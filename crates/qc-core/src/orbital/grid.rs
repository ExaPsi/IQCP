//! MO grid evaluator
//!
//! Evaluates molecular orbital wavefunctions on a uniform 3D grid for
//! isosurface extraction and visualization.
//!
//! # Algorithm
//!
//! For each grid point r = (x, y, z), the MO value is:
//!
//! ```text
//! psi_i(r) = sum_mu { C_{mu,i} * chi_mu(r) }
//! ```
//!
//! where chi_mu is the contracted Gaussian basis function:
//!
//! ```text
//! chi_mu(r) = sum_k { c_k * N(alpha_k, l,m,n) * (x-Ax)^l (y-Ay)^m (z-Az)^n * exp(-alpha_k |r-A|^2) }
//! ```
//!
//! # Performance
//!
//! Two key optimizations are used:
//!
//! 1. **Gaussian screening**: Skip shell contributions where the most diffuse
//!    primitive gives exp(-alpha_min * r^2) < 1e-15, i.e., alpha_min * r^2 > 34.54.
//!
//! 2. **Shell-batched evaluation**: The outer loop iterates over shells (keeping
//!    shell data in cache), and the inner loop iterates over grid points.
//!
//! # Normalization
//!
//! The normalization follows the same convention as the SCF integral engine
//! (see `integrals::overlap::cartesian_gaussian_normalization`):
//!
//! ```text
//! N(alpha, i, j, k) = (2*alpha/pi)^(3/4) * (4*alpha)^((i+j+k)/2)
//!                    / sqrt((2i-1)!! * (2j-1)!! * (2k-1)!!)
//! ```
//!
//! The contraction coefficients stored in `GaussianPrimitive::coefficient` are
//! the raw coefficients from the basis set library (e.g., PySCF), and the
//! primitive normalization is applied explicitly during evaluation.
//!
//! # References
//!
//! - Szabo & Ostlund, Eq. 3.203 (Gaussian normalization)
//! - PySCF `mol.eval_gto('GTOval', coords)` for validation

use crate::basis::{AngularMomentum, BasisSet, ContractedShell};
use crate::integrals::{cartesian_components, CartesianPower};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use thiserror::Error;

// =============================================================================
// Error Types
// =============================================================================

/// Errors from orbital grid evaluation
#[derive(Debug, Error, Clone, PartialEq)]
pub enum OrbitalError {
    /// MO coefficient vector has wrong length
    #[error("MO coefficients length {actual} does not match basis size {expected}")]
    CoefficientLengthMismatch { expected: usize, actual: usize },

    /// Grid dimensions are invalid
    #[error("Grid dimension must be >= 2, got {0}")]
    InvalidGridDimension(usize),

    /// Grid spacing is invalid
    #[error("Grid spacing must be positive, got {0}")]
    InvalidGridSpacing(f64),

    /// Basis set is empty
    #[error("Basis set has no shells")]
    EmptyBasis,

    /// Density matrix has wrong size
    #[error("Density matrix size {actual} does not match expected {expected} (nbf^2)")]
    DensityMatrixSizeMismatch { expected: usize, actual: usize },

    /// Element not supported for promolecule density (Z > 18 or Z < 1)
    #[error("No atomic density profile available for element Z={0} (only H-Ar supported)")]
    UnsupportedElement(u32),

    /// Two grids have different sizes
    #[error(
        "Grid size mismatch: molecular density has {expected} points, promolecule has {actual}"
    )]
    GridSizeMismatch { expected: usize, actual: usize },
}

// =============================================================================
// Configuration and Result Types
// =============================================================================

/// Configuration for grid auto-determination
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GridConfig {
    /// Padding factor for grid extent beyond outermost nucleus.
    /// Default: 4.0 (extends 4 * alpha_min^{-1/2} beyond atoms).
    pub padding_factor: f64,

    /// Grid spacing in Bohr.
    /// If None, auto-determined from grid_dims.
    pub spacing: Option<f64>,

    /// Grid dimensions [nx, ny, nz].
    /// If None with spacing, auto-determined from bounds.
    pub dims: Option<[usize; 3]>,
}

impl Default for GridConfig {
    fn default() -> Self {
        Self {
            padding_factor: 4.0,
            spacing: None,
            dims: None,
        }
    }
}

/// Input for MO grid evaluation (WASM-friendly)
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MoGridInput {
    /// MO coefficient vector (C_{mu,i} for basis function mu)
    pub mo_coefficients: Vec<f64>,

    /// Grid origin [x, y, z] in Bohr
    pub grid_origin: [f64; 3],

    /// Grid spacing in Bohr (uniform in all directions)
    pub grid_spacing: f64,

    /// Grid dimensions [nx, ny, nz]
    pub grid_dims: [usize; 3],
}

/// Result from MO grid evaluation
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MoGridResult {
    /// Grid values (flattened, C-order: x-slowest, z-fastest)
    /// index = ix * ny * nz + iy * nz + iz
    pub values: Vec<f64>,

    /// Grid origin [x, y, z] in Bohr
    pub grid_origin: [f64; 3],

    /// Grid spacing in Bohr
    pub grid_spacing: f64,

    /// Grid dimensions [nx, ny, nz]
    pub grid_dims: [usize; 3],

    /// Maximum absolute value in the grid
    pub max_abs_value: f64,

    /// Approximate norm-squared integral (sum(psi^2) * dV)
    pub norm_sq_integral: f64,

    /// Computation time in milliseconds
    pub compute_time_ms: f64,
}

// =============================================================================
// Gaussian Screening Threshold
// =============================================================================

/// Screening threshold: skip when exp(-alpha * r^2) < 1e-15
/// Equivalent to alpha * r^2 > -ln(1e-15) = 34.5388...
pub(crate) const GAUSSIAN_SCREENING_THRESHOLD: f64 = 34.5388;

// =============================================================================
// Grid Bounds Auto-Determination
// =============================================================================

/// Compute automatic grid bounds extending beyond the molecule.
///
/// The grid extends `padding_factor * alpha_min^{-1/2}` beyond the outermost
/// nucleus in each Cartesian direction, where `alpha_min` is the global minimum
/// Gaussian exponent across all shells.
///
/// # Arguments
///
/// * `basis` - The basis set containing atoms and shells
/// * `padding_factor` - Multiplicative factor (default 4.0)
///
/// # Returns
///
/// `(min_corner, max_corner)` - The grid bounding box in Bohr
///
/// # Reference
///
/// Phase 2 TDD: "Extend 4 * alpha_min^{-1/2} beyond the outermost nucleus"
pub fn auto_grid_bounds(basis: &BasisSet, padding_factor: f64) -> ([f64; 3], [f64; 3]) {
    // Find bounding box of all atom positions
    let mut min_pos = [f64::INFINITY; 3];
    let mut max_pos = [f64::NEG_INFINITY; 3];

    for atom in &basis.atoms {
        for d in 0..3 {
            min_pos[d] = min_pos[d].min(atom.position[d]);
            max_pos[d] = max_pos[d].max(atom.position[d]);
        }
    }

    // Find global minimum exponent across all shells
    let alpha_min = basis
        .shells
        .iter()
        .flat_map(|s| s.primitives.iter().map(|p| p.exponent))
        .fold(f64::INFINITY, f64::min);

    // Compute padding: padding_factor / sqrt(alpha_min)
    let padding = if alpha_min > 0.0 && alpha_min.is_finite() {
        padding_factor / alpha_min.sqrt()
    } else {
        // Fallback: 10 Bohr padding
        10.0
    };

    let min_corner = [
        min_pos[0] - padding,
        min_pos[1] - padding,
        min_pos[2] - padding,
    ];
    let max_corner = [
        max_pos[0] + padding,
        max_pos[1] + padding,
        max_pos[2] + padding,
    ];

    (min_corner, max_corner)
}

// =============================================================================
// Core Grid Evaluation
// =============================================================================

/// Evaluate a molecular orbital on a uniform 3D grid.
///
/// Computes psi_i(r) = sum_mu { C_{mu,i} * chi_mu(r) } at every grid point,
/// where chi_mu(r) are the contracted Gaussian basis functions.
///
/// # Grid Indexing Convention
///
/// C-order (x-slowest, z-fastest):
/// ```text
/// index = ix * ny * nz + iy * nz + iz
/// ```
///
/// # Arguments
///
/// * `mo_coefficients` - MO coefficient vector C_{mu,i} for all basis functions mu
/// * `basis` - The basis set (provides shell centers, exponents, coefficients)
/// * `grid_origin` - Origin [x, y, z] of the grid in Bohr
/// * `grid_spacing` - Uniform grid spacing in Bohr
/// * `grid_dims` - Grid dimensions [nx, ny, nz]
///
/// # Returns
///
/// A flat vector of grid values in C-order (length nx * ny * nz).
///
/// # Errors
///
/// Returns error if MO coefficient length does not match basis size,
/// or if grid dimensions are invalid.
///
/// # Performance
///
/// Uses shell-batched evaluation with Gaussian screening for efficiency.
/// A 50^3 grid for H2O STO-3G completes in under 2 seconds.
pub fn evaluate_mo_on_grid(
    mo_coefficients: &[f64],
    basis: &BasisSet,
    grid_origin: [f64; 3],
    grid_spacing: f64,
    grid_dims: [usize; 3],
) -> Result<Vec<f64>, OrbitalError> {
    let [nx, ny, nz] = grid_dims;
    let n_basis = basis.n_basis;

    // Validate inputs
    if mo_coefficients.len() != n_basis {
        return Err(OrbitalError::CoefficientLengthMismatch {
            expected: n_basis,
            actual: mo_coefficients.len(),
        });
    }
    if nx < 2 || ny < 2 || nz < 2 {
        return Err(OrbitalError::InvalidGridDimension(nx.min(ny).min(nz)));
    }
    if grid_spacing <= 0.0 {
        return Err(OrbitalError::InvalidGridSpacing(grid_spacing));
    }
    if basis.shells.is_empty() {
        return Err(OrbitalError::EmptyBasis);
    }

    let total_points = nx * ny * nz;
    let mut grid_values = vec![0.0f64; total_points];

    // Shell-batched evaluation:
    // Outer loop: iterate over shells (keeping shell data in cache)
    // Inner loop: iterate over grid points
    let mut basis_offset = 0;

    for shell in &basis.shells {
        let n_funcs = shell.n_basis_functions();
        let l = shell.l_value();

        // Get Cartesian components for this angular momentum
        // This should never fail for supported angular momenta (S, P, D)
        let components = match cartesian_components(l) {
            Ok(c) => c,
            Err(_) => {
                basis_offset += n_funcs;
                continue;
            }
        };

        // Extract MO coefficients for this shell's basis functions
        let shell_mo_coeffs: Vec<f64> = (0..n_funcs)
            .map(|i| mo_coefficients[basis_offset + i])
            .collect();

        // Quick check: if all MO coefficients for this shell are zero, skip
        let max_coeff = shell_mo_coeffs
            .iter()
            .map(|c| c.abs())
            .fold(0.0f64, f64::max);
        if max_coeff < 1e-20 {
            basis_offset += n_funcs;
            continue;
        }

        // Get minimum exponent for screening
        let alpha_min_shell = shell.min_exponent();

        // Shell center
        let [ax, ay, az] = shell.center;

        // Evaluate this shell's contribution at each grid point
        evaluate_shell_on_grid(
            shell,
            &components,
            &shell_mo_coeffs,
            alpha_min_shell,
            ax,
            ay,
            az,
            grid_origin,
            grid_spacing,
            grid_dims,
            &mut grid_values,
        );

        basis_offset += n_funcs;
    }

    Ok(grid_values)
}

/// Evaluate a single shell's contribution to the MO grid.
///
/// This is the inner kernel of the grid evaluator. For each grid point,
/// it computes the contracted Gaussian basis functions for this shell
/// and adds the weighted contribution to the grid values.
///
/// Uses a general approach that handles all angular momenta (S, P, D)
/// with proper per-component normalization.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn evaluate_shell_on_grid(
    shell: &ContractedShell,
    components: &[CartesianPower],
    shell_mo_coeffs: &[f64],
    alpha_min_shell: f64,
    ax: f64,
    ay: f64,
    az: f64,
    grid_origin: [f64; 3],
    grid_spacing: f64,
    grid_dims: [usize; 3],
    grid_values: &mut [f64],
) {
    let [nx, ny, nz] = grid_dims;
    let [x0, y0, z0] = grid_origin;

    // For s and p shells, all components share the same normalization
    // (because the normalization only depends on total L for these cases).
    // For d shells (and higher), each component has a different normalization.
    // We use the general approach for all, which is clean and correct.
    // The inner loop is optimized for the common s and p cases.

    let is_s = shell.angular_momentum == AngularMomentum::S;
    let is_p = shell.angular_momentum == AngularMomentum::P;

    // Precompute normalization factors for each primitive and component
    // to avoid recomputation inside the innermost loop.
    // norm_coef[prim_idx][comp_idx] = coefficient * normalization
    let n_prims = shell.primitives.len();
    let n_comps = components.len();
    let mut norm_coef: Vec<Vec<f64>> = Vec::with_capacity(n_prims);
    let mut exponents: Vec<f64> = Vec::with_capacity(n_prims);

    for prim in &shell.primitives {
        exponents.push(prim.exponent);
        let mut nc = Vec::with_capacity(n_comps);
        for comp in components {
            nc.push(prim.coefficient * cartesian_norm(prim.exponent, comp));
        }
        norm_coef.push(nc);
    }

    // Iterate over grid points
    for ix in 0..nx {
        let gx = x0 + ix as f64 * grid_spacing;
        let dx = gx - ax;
        let dx2 = dx * dx;

        // X-direction screening
        if alpha_min_shell * dx2 > GAUSSIAN_SCREENING_THRESHOLD {
            continue;
        }

        for iy in 0..ny {
            let gy = y0 + iy as f64 * grid_spacing;
            let dy = gy - ay;
            let dy2 = dy * dy;
            let dxy2 = dx2 + dy2;

            // XY screening
            if alpha_min_shell * dxy2 > GAUSSIAN_SCREENING_THRESHOLD {
                continue;
            }

            for iz in 0..nz {
                let gz = z0 + iz as f64 * grid_spacing;
                let dz = gz - az;
                let dist_sq = dxy2 + dz * dz;

                // Full Gaussian screening
                if alpha_min_shell * dist_sq > GAUSSIAN_SCREENING_THRESHOLD {
                    continue;
                }

                let flat_idx = ix * ny * nz + iy * nz + iz;

                if is_s {
                    // S-shell fast path: single component, no angular part
                    let mut val = 0.0;
                    for (p, &alpha) in exponents.iter().enumerate() {
                        val += norm_coef[p][0] * (-alpha * dist_sq).exp();
                    }
                    grid_values[flat_idx] += shell_mo_coeffs[0] * val;
                } else if is_p {
                    // P-shell fast path: all 3 components share the same
                    // normalization (since (1,0,0), (0,1,0), (0,0,1) all have L=1
                    // with the same double factorial structure).
                    // Compute shared radial sum once.
                    let mut radial = 0.0;
                    for (p, &alpha) in exponents.iter().enumerate() {
                        radial += norm_coef[p][0] * (-alpha * dist_sq).exp();
                    }
                    let disps = [dx, dy, dz];
                    for comp_idx in 0..3 {
                        grid_values[flat_idx] +=
                            shell_mo_coeffs[comp_idx] * disps[comp_idx] * radial;
                    }
                } else {
                    // General case (D-shells and above): each component has
                    // different normalization.
                    let disps = [dx, dy, dz];
                    for (comp_idx, powers) in components.iter().enumerate() {
                        let angular = angular_factor(&disps, powers);
                        if angular.abs() < 1e-30 {
                            continue;
                        }
                        let mut radial_sum = 0.0;
                        for (p, &alpha) in exponents.iter().enumerate() {
                            radial_sum += norm_coef[p][comp_idx] * (-alpha * dist_sq).exp();
                        }
                        grid_values[flat_idx] += shell_mo_coeffs[comp_idx] * angular * radial_sum;
                    }
                }
            }
        }
    }
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Compute angular factor (x-Ax)^i * (y-Ay)^j * (z-Az)^k
#[inline]
pub(crate) fn angular_factor(disps: &[f64; 3], powers: &CartesianPower) -> f64 {
    let mut result = 1.0;
    if powers.i > 0 {
        result *= disps[0].powi(powers.i as i32);
    }
    if powers.j > 0 {
        result *= disps[1].powi(powers.j as i32);
    }
    if powers.k > 0 {
        result *= disps[2].powi(powers.k as i32);
    }
    result
}

/// Double factorial (2n-1)!! = 1 * 3 * 5 * ... * (2n-1)
/// Convention: (-1)!! = 0!! = 1
#[inline]
pub(crate) fn double_factorial(n: u32) -> u64 {
    if n <= 1 {
        1
    } else {
        let mut result = 1u64;
        let mut k = n;
        while k > 1 {
            result *= k as u64;
            k -= 2;
        }
        result
    }
}

/// Compute Cartesian Gaussian normalization factor.
///
/// ```text
/// N(alpha, i, j, k) = (2*alpha/pi)^(3/4) * (4*alpha)^((i+j+k)/2)
///                    / sqrt((2i-1)!! * (2j-1)!! * (2k-1)!!)
/// ```
///
/// Matches `integrals::overlap::cartesian_gaussian_normalization` exactly.
///
/// Reference: Szabo & Ostlund, Eq. 3.203
#[inline]
pub(crate) fn cartesian_norm(alpha: f64, powers: &CartesianPower) -> f64 {
    let base = (2.0 * alpha / PI).powf(0.75);
    let l = powers.i + powers.j + powers.k;

    if l == 0 {
        return base;
    }

    let four_alpha = 4.0 * alpha;
    let ang_num = four_alpha.powf(l as f64 / 2.0);

    let denom_i = double_factorial(if powers.i > 0 { 2 * powers.i - 1 } else { 0 }) as f64;
    let denom_j = double_factorial(if powers.j > 0 { 2 * powers.j - 1 } else { 0 }) as f64;
    let denom_k = double_factorial(if powers.k > 0 { 2 * powers.k - 1 } else { 0 }) as f64;

    let ang_denom = (denom_i * denom_j * denom_k).sqrt();

    base * ang_num / ang_denom
}

// =============================================================================
// Evaluate MO on grid (full result struct version)
// =============================================================================

/// Evaluate a molecular orbital on an auto-determined or specified grid.
///
/// This is the high-level entry point that returns a full result struct
/// with metadata including norm-squared integral and timing.
///
/// # Arguments
///
/// * `mo_coefficients` - MO coefficient vector
/// * `basis` - The molecular basis set
/// * `grid_origin` - Grid origin in Bohr
/// * `grid_spacing` - Grid spacing in Bohr
/// * `grid_dims` - Grid dimensions [nx, ny, nz]
///
/// # Returns
///
/// A `MoGridResult` containing the grid values and metadata.
pub fn evaluate_mo_on_grid_full(
    mo_coefficients: &[f64],
    basis: &BasisSet,
    grid_origin: [f64; 3],
    grid_spacing: f64,
    grid_dims: [usize; 3],
) -> Result<MoGridResult, OrbitalError> {
    let values = evaluate_mo_on_grid(mo_coefficients, basis, grid_origin, grid_spacing, grid_dims)?;

    // Compute metadata
    let dv = grid_spacing.powi(3);
    let norm_sq_integral: f64 = values.iter().map(|v| v * v).sum::<f64>() * dv;
    let max_abs_value = values.iter().map(|v| v.abs()).fold(0.0f64, f64::max);

    Ok(MoGridResult {
        values,
        grid_origin,
        grid_spacing,
        grid_dims,
        max_abs_value,
        norm_sq_integral,
        // Timing is handled by the caller (WASM layer uses js_sys::Date::now()
        // since std::time::Instant is not supported on wasm32-unknown-unknown)
        compute_time_ms: 0.0,
    })
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::{Atom, BasisSet};
    use approx::assert_abs_diff_eq;

    // -------------------------------------------------------------------------
    // Helper: Build H2 STO-3G basis
    // -------------------------------------------------------------------------

    fn h2_basis() -> BasisSet {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        BasisSet::build(vec![h1, h2], "sto-3g").unwrap()
    }

    fn h2o_basis() -> BasisSet {
        let o = Atom::new(8, [0.0, 0.0, 0.221665]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430901, -0.886659]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430901, -0.886659]).unwrap();
        BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap()
    }

    // -------------------------------------------------------------------------
    // Normalization tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_cartesian_norm_s_type() {
        // s-type: N = (2*alpha/pi)^(3/4)
        let powers = CartesianPower::new(0, 0, 0);
        let alpha = 3.42525091;
        let norm = cartesian_norm(alpha, &powers);
        let expected = (2.0 * alpha / PI).powf(0.75);
        assert_abs_diff_eq!(norm, expected, epsilon = 1e-12);
    }

    #[test]
    fn test_cartesian_norm_p_type() {
        // p-type: N = (2*alpha/pi)^(3/4) * (4*alpha)^(1/2)
        let powers_x = CartesianPower::new(1, 0, 0);
        let alpha = 5.0331513;
        let norm = cartesian_norm(alpha, &powers_x);
        let expected = (2.0 * alpha / PI).powf(0.75) * (4.0 * alpha).sqrt();
        assert_abs_diff_eq!(norm, expected, epsilon = 1e-12);

        // All p components should have the same normalization
        let norm_y = cartesian_norm(alpha, &CartesianPower::new(0, 1, 0));
        let norm_z = cartesian_norm(alpha, &CartesianPower::new(0, 0, 1));
        assert_abs_diff_eq!(norm, norm_y, epsilon = 1e-15);
        assert_abs_diff_eq!(norm, norm_z, epsilon = 1e-15);
    }

    #[test]
    fn test_cartesian_norm_d_type() {
        let alpha = 1.0;

        // dxx (2,0,0): N = (2a/pi)^(3/4) * (4a)^1 / sqrt(3!!)
        //            = (2a/pi)^(3/4) * 4a / sqrt(3)
        let norm_xx = cartesian_norm(alpha, &CartesianPower::new(2, 0, 0));
        let expected_xx = (2.0 * alpha / PI).powf(0.75) * 4.0 * alpha / 3.0f64.sqrt();
        assert_abs_diff_eq!(norm_xx, expected_xx, epsilon = 1e-12);

        // dxy (1,1,0): N = (2a/pi)^(3/4) * (4a)^1 / sqrt(1*1*1) = (2a/pi)^(3/4) * 4a
        let norm_xy = cartesian_norm(alpha, &CartesianPower::new(1, 1, 0));
        let expected_xy = (2.0 * alpha / PI).powf(0.75) * 4.0 * alpha;
        assert_abs_diff_eq!(norm_xy, expected_xy, epsilon = 1e-12);

        // dxy normalization should be sqrt(3) times dxx normalization
        assert_abs_diff_eq!(norm_xy / norm_xx, 3.0f64.sqrt(), epsilon = 1e-12);
    }

    // -------------------------------------------------------------------------
    // Grid bounds tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_auto_grid_bounds_h2() {
        let basis = h2_basis();
        let (min_corner, max_corner) = auto_grid_bounds(&basis, 4.0);

        // H2 atoms at (0,0,0) and (0,0,1.4)
        // Min exponent in STO-3G H: 0.1688554
        // padding = 4 / sqrt(0.1688554) = 9.734
        let alpha_min: f64 = 0.1688554;
        let padding = 4.0 / alpha_min.sqrt();

        assert_abs_diff_eq!(min_corner[0], -padding, epsilon = 0.01);
        assert_abs_diff_eq!(min_corner[1], -padding, epsilon = 0.01);
        assert_abs_diff_eq!(min_corner[2], -padding, epsilon = 0.01);
        assert_abs_diff_eq!(max_corner[0], padding, epsilon = 0.01);
        assert_abs_diff_eq!(max_corner[1], padding, epsilon = 0.01);
        assert_abs_diff_eq!(max_corner[2], 1.4 + padding, epsilon = 0.01);
    }

    #[test]
    fn test_auto_grid_bounds_h2o() {
        let basis = h2o_basis();
        let (min_corner, max_corner) = auto_grid_bounds(&basis, 4.0);

        // H2O: atoms span y=[-1.431, 1.431], z=[-0.887, 0.222]
        // padding should be at least 7.9 (from O 2sp exponent 0.255580)
        // but global min is H's 0.1688554, so padding = 9.734

        // Verify bounds are reasonable
        assert!(min_corner[0] < -9.0);
        assert!(max_corner[0] > 9.0);
        assert!(min_corner[1] < -10.0); // -1.431 - 9.734
        assert!(max_corner[1] > 10.0); // 1.431 + 9.734
    }

    // -------------------------------------------------------------------------
    // Input validation tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_coefficient_length_mismatch() {
        let basis = h2_basis();
        let wrong_coeffs = vec![1.0]; // should be 2
        let result = evaluate_mo_on_grid(&wrong_coeffs, &basis, [0.0; 3], 1.0, [10, 10, 10]);
        assert!(matches!(
            result,
            Err(OrbitalError::CoefficientLengthMismatch { .. })
        ));
    }

    #[test]
    fn test_invalid_grid_dimension() {
        let basis = h2_basis();
        let coeffs = vec![0.5, 0.5];
        let result = evaluate_mo_on_grid(&coeffs, &basis, [0.0; 3], 1.0, [1, 10, 10]);
        assert!(matches!(result, Err(OrbitalError::InvalidGridDimension(_))));
    }

    #[test]
    fn test_invalid_grid_spacing() {
        let basis = h2_basis();
        let coeffs = vec![0.5, 0.5];
        let result = evaluate_mo_on_grid(&coeffs, &basis, [0.0; 3], 0.0, [10, 10, 10]);
        assert!(matches!(result, Err(OrbitalError::InvalidGridSpacing(_))));
    }

    // -------------------------------------------------------------------------
    // H2 sigma_g golden tests (PySCF reference)
    // -------------------------------------------------------------------------

    #[test]
    fn test_h2_sigma_g_point_values() {
        // PySCF reference values from tests/golden/orbital/h2_point_values.json
        // MO coefficients: [0.5489340403996354, 0.5489340403996352]
        let basis = h2_basis();
        let mo_coeffs = vec![0.5489340403996354, 0.5489340403996352];

        // Test at bond midpoint r=(0,0,0.7)
        // PySCF: MO0 = 3.577196327424507e-01
        let values = evaluate_mo_on_grid(
            &mo_coeffs,
            &basis,
            [-0.5, -0.5, 0.2], // small grid around midpoint
            0.5,
            [3, 3, 3],
        )
        .unwrap();

        // The grid point at (ix=1, iy=1, iz=1) = (-0.5+0.5, -0.5+0.5, 0.2+0.5) = (0.0, 0.0, 0.7)
        let midpoint_value = values[1 * 3 * 3 + 1 * 3 + 1];
        assert_abs_diff_eq!(midpoint_value, 3.577196327424507e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2_sigma_g_on_atom() {
        // PySCF: r=(0,0,0), MO0 = 4.212434409465830e-01
        let basis = h2_basis();
        let mo_coeffs = vec![0.5489340403996354, 0.5489340403996352];

        let values =
            evaluate_mo_on_grid(&mo_coeffs, &basis, [0.0, 0.0, 0.0], 0.7, [3, 3, 3]).unwrap();

        // Grid point (0,0,0) = origin = atom 1
        let atom_value = values[0];
        assert_abs_diff_eq!(atom_value, 4.212434409465830e-01, epsilon = 1e-8);

        // Grid point (0,0,2) = (0, 0, 1.4) = atom 2
        let atom2_value = values[0 * 3 * 3 + 0 * 3 + 2];
        assert_abs_diff_eq!(atom2_value, 4.212434409465828e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2_sigma_g_off_axis() {
        // PySCF: r=(1,0,0.7), MO0 = 1.891668849482989e-01
        let basis = h2_basis();
        let mo_coeffs = vec![0.5489340403996354, 0.5489340403996352];

        let values =
            evaluate_mo_on_grid(&mo_coeffs, &basis, [1.0, 0.0, 0.7], 0.5, [2, 2, 2]).unwrap();

        let val = values[0]; // (0,0,0) = (1.0, 0.0, 0.7)
        assert_abs_diff_eq!(val, 1.891668849482989e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2_sigma_g_general_point() {
        // PySCF: r=(0.5,0.5,0.5), MO0 = 2.474670861940863e-01
        let basis = h2_basis();
        let mo_coeffs = vec![0.5489340403996354, 0.5489340403996352];

        let values =
            evaluate_mo_on_grid(&mo_coeffs, &basis, [0.5, 0.5, 0.5], 1.0, [2, 2, 2]).unwrap();

        let val = values[0]; // (0,0,0) = (0.5, 0.5, 0.5)
        assert_abs_diff_eq!(val, 2.474670861940863e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2_sigma_g_far_point() {
        // PySCF: r=(2,0,0), MO0 = 5.575974312159111e-02
        let basis = h2_basis();
        let mo_coeffs = vec![0.5489340403996354, 0.5489340403996352];

        let values =
            evaluate_mo_on_grid(&mo_coeffs, &basis, [2.0, 0.0, 0.0], 1.0, [2, 2, 2]).unwrap();

        let val = values[0];
        assert_abs_diff_eq!(val, 5.575974312159111e-02, epsilon = 1e-8);
    }

    #[test]
    fn test_h2_sigma_g_negative_z() {
        // PySCF: r=(0,0,-2), MO0 = 4.205664508227352e-02
        let basis = h2_basis();
        let mo_coeffs = vec![0.5489340403996354, 0.5489340403996352];

        let values =
            evaluate_mo_on_grid(&mo_coeffs, &basis, [0.0, 0.0, -2.0], 1.0, [2, 2, 2]).unwrap();

        let val = values[0];
        assert_abs_diff_eq!(val, 4.205664508227352e-02, epsilon = 1e-8);
    }

    // -------------------------------------------------------------------------
    // H2O tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_h2o_homo_at_x_axis() {
        // H2O HOMO is pure px orbital.
        // PySCF: r=(1,0,0), HOMO = 4.293649798818524e-01
        // MO coefficients: [~0, ~0, 1.0, ~0, ~0, ~0, ~0]
        let basis = h2o_basis();
        let mo_coeffs = vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0];

        let values =
            evaluate_mo_on_grid(&mo_coeffs, &basis, [1.0, 0.0, 0.0], 1.0, [2, 2, 2]).unwrap();

        let val = values[0]; // (1.0, 0.0, 0.0)
        assert_abs_diff_eq!(val, 4.293649798818524e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2o_homo_on_yz_plane() {
        // H2O HOMO (px): should be zero on the yz-plane (x=0)
        // PySCF confirms: r=(0,0,0.222), HOMO ~ 2e-16 (zero)
        let basis = h2o_basis();
        let mo_coeffs = vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0];

        let values =
            evaluate_mo_on_grid(&mo_coeffs, &basis, [0.0, 0.0, 0.221665], 1.0, [2, 2, 2]).unwrap();

        // At x=0, px orbital is zero
        let val = values[0];
        assert!(val.abs() < 1e-14);
    }

    #[test]
    fn test_h2o_mo3_mixed_character() {
        // H2O MO3 (lone pair) has mixed s, pz, H character
        // PySCF MO coefficients (full precision from PySCF 2.11.0):
        // PySCF: r=(0,0,0.222), MO3 = -8.6806237006611986e-01
        let basis = h2o_basis();
        let mo_coeffs = vec![
            -1.0317195713388164e-01,
            5.3679807124188206e-01,
            0.0,
            0.0,
            7.7649134735796044e-01,
            -2.7810301986702651e-01,
            -2.7810301986702557e-01,
        ];

        let values =
            evaluate_mo_on_grid(&mo_coeffs, &basis, [0.0, 0.0, 0.221665], 1.0, [2, 2, 2]).unwrap();

        let val = values[0];
        assert_abs_diff_eq!(val, -8.6806237006611986e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2o_mo3_general_point() {
        // PySCF: r=(0.5,0.5,0.5), MO3 = 3.1052942251890669e-01
        let basis = h2o_basis();
        let mo_coeffs = vec![
            -1.0317195713388164e-01,
            5.3679807124188206e-01,
            0.0,
            0.0,
            7.7649134735796044e-01,
            -2.7810301986702651e-01,
            -2.7810301986702557e-01,
        ];

        let values =
            evaluate_mo_on_grid(&mo_coeffs, &basis, [0.5, 0.5, 0.5], 1.0, [2, 2, 2]).unwrap();

        let val = values[0];
        assert_abs_diff_eq!(val, 3.1052942251890669e-01, epsilon = 1e-8);
    }

    // -------------------------------------------------------------------------
    // Norm-squared integral tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_h2_norm_squared_integral() {
        // For the H2 sigma_g orbital, the norm-squared integral should be ~1.0
        // on a sufficiently large grid.
        let basis = h2_basis();
        let mo_coeffs = vec![0.5489340403996354, 0.5489340403996352];

        let (min_corner, max_corner) = auto_grid_bounds(&basis, 4.0);
        let nx = 40;
        let spacing = (max_corner[0] - min_corner[0]) / (nx as f64 - 1.0);

        // Use uniform spacing in all directions for simplicity
        let ny = ((max_corner[1] - min_corner[1]) / spacing).ceil() as usize + 1;
        let nz = ((max_corner[2] - min_corner[2]) / spacing).ceil() as usize + 1;

        let result =
            evaluate_mo_on_grid_full(&mo_coeffs, &basis, min_corner, spacing, [nx, ny, nz])
                .unwrap();

        // Norm should be approximately 1.0 within 5% grid quadrature error
        assert!(
            result.norm_sq_integral > 0.95 && result.norm_sq_integral < 1.05,
            "Norm^2 = {:.6}, expected ~1.0",
            result.norm_sq_integral
        );
    }

    // -------------------------------------------------------------------------
    // Gaussian screening tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_screening_very_far_point() {
        // A point very far from the molecule should have ~zero orbital value
        let basis = h2_basis();
        let mo_coeffs = vec![0.5489340403996354, 0.5489340403996352];

        let values =
            evaluate_mo_on_grid(&mo_coeffs, &basis, [50.0, 50.0, 50.0], 1.0, [2, 2, 2]).unwrap();

        // All values should be essentially zero due to screening
        for val in &values {
            assert!(val.abs() < 1e-14, "Expected ~0 at far point, got {}", val);
        }
    }

    // -------------------------------------------------------------------------
    // Grid result metadata tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_mo_grid_result_metadata() {
        let basis = h2_basis();
        let mo_coeffs = vec![0.5489340403996354, 0.5489340403996352];

        let result =
            evaluate_mo_on_grid_full(&mo_coeffs, &basis, [-5.0, -5.0, -5.0], 1.0, [11, 11, 13])
                .unwrap();

        assert_eq!(result.values.len(), 11 * 11 * 13);
        assert_eq!(result.grid_dims, [11, 11, 13]);
        assert!(result.max_abs_value > 0.0);
        assert!(result.compute_time_ms >= 0.0);
    }
}
