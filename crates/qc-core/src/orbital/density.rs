//! Total electron density grid evaluator
//!
//! Evaluates the total electron density rho(r) on a uniform 3D grid from
//! the density matrix and basis set:
//!
//! ```text
//! rho(r) = sum_{mu,nu} D_{mu,nu} * chi_mu(r) * chi_nu(r)
//! ```
//!
//! where D is the RHF density matrix (already includes the factor of 2,
//! Szabo & Ostlund Eq. 3.145) and chi_mu(r) are contracted Gaussian
//! basis functions.
//!
//! # Algorithm
//!
//! The evaluator uses a grid-point-centric approach with shell-batched
//! basis function evaluation:
//!
//! 1. For each grid point, evaluate all basis functions chi_mu(r)
//! 2. Contract with the density matrix using the symmetry-optimized formula:
//!    ```text
//!    rho = sum_mu D[mu,mu] * chi[mu]^2 + 2 * sum_{mu < nu} D[mu,nu] * chi[mu] * chi[nu]
//!    ```
//!
//! # Performance
//!
//! - Gaussian screening: skip contributions where exp(-alpha * r^2) < 1e-15
//! - The chi vector is allocated once and zeroed each iteration
//! - D is symmetric, so only upper-triangle + diagonal is needed
//!
//! # References
//!
//! - Szabo & Ostlund (1996). "Modern Quantum Chemistry". Dover, Eq. 3.145.
//! - PySCF density evaluation: `pyscf.dft.numint.eval_rho`

use crate::basis::{AngularMomentum, BasisSet};
use crate::integrals::{cartesian_components, CartesianPower};
use serde::{Deserialize, Serialize};

use super::grid::{angular_factor, cartesian_norm, OrbitalError, GAUSSIAN_SCREENING_THRESHOLD};

// =============================================================================
// Result Type
// =============================================================================

/// Result of evaluating total electron density on a 3D grid.
///
/// rho(r) = sum_{mu,nu} D_{mu,nu} chi_mu(r) chi_nu(r)
/// where D already contains the factor of 2 (RHF convention, Szabo & Ostlund Eq. 3.145).
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct DensityGridResult {
    /// Flat density values rho(r) at each grid point (C-order: x-slowest, z-fastest)
    pub values: Vec<f64>,
    /// Grid dimensions [nx, ny, nz]
    pub grid_dims: [usize; 3],
    /// Grid origin [x, y, z] in bohr
    pub grid_origin: [f64; 3],
    /// Grid spacing in bohr
    pub grid_spacing: f64,
    /// Integrated density: N_integrated = sum rho(r_i) * dV
    pub integrated_density: f64,
    /// Expected number of electrons (from basis set / SCF input)
    pub n_electrons_expected: usize,
    /// Maximum density value in the grid
    pub max_density: f64,
    /// Computation time in milliseconds
    pub compute_time_ms: f64,
}

// =============================================================================
// Core Density Evaluation
// =============================================================================

/// Evaluate total electron density on a uniform 3D grid.
///
/// Computes rho(r) = sum_{mu,nu} D_{mu,nu} * chi_mu(r) * chi_nu(r) at every
/// grid point, where D is the RHF density matrix and chi_mu(r) are the
/// contracted Gaussian basis functions.
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
/// * `density_matrix` - Row-major density matrix D_{mu,nu} of size nbf*nbf.
///   This is the full RHF density matrix that includes the factor of 2.
/// * `nbf` - Number of basis functions
/// * `n_electrons` - Expected number of electrons (for validation metadata)
/// * `basis` - The basis set (provides shell centers, exponents, coefficients)
/// * `grid_origin` - Origin [x, y, z] of the grid in Bohr
/// * `grid_spacing` - Uniform grid spacing in Bohr
/// * `grid_dims` - Grid dimensions [nx, ny, nz]
///
/// # Returns
///
/// A `DensityGridResult` containing grid values and integration metadata.
///
/// # Errors
///
/// Returns error if density matrix size does not match nbf^2,
/// if grid dimensions are invalid, or if the basis set is empty.
///
/// # References
///
/// Szabo & Ostlund, Eq. 3.145: D_mu,nu = 2 * sum_a^{N/2} C_{mu,a} C_{nu,a}
pub fn evaluate_density_on_grid(
    density_matrix: &[f64],
    nbf: usize,
    n_electrons: usize,
    basis: &BasisSet,
    grid_origin: [f64; 3],
    grid_spacing: f64,
    grid_dims: [usize; 3],
) -> Result<DensityGridResult, OrbitalError> {
    let [nx, ny, nz] = grid_dims;

    // ---- Input validation ----
    if density_matrix.len() != nbf * nbf {
        return Err(OrbitalError::DensityMatrixSizeMismatch {
            expected: nbf * nbf,
            actual: density_matrix.len(),
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

    // ---- Precompute per-shell data ----
    // For each shell: normalization factors, exponents, center, min exponent,
    // cartesian components, and basis offset.
    let shell_data: Vec<ShellPrecomputed> = {
        let mut data = Vec::with_capacity(basis.shells.len());
        let mut offset = 0usize;
        for shell in &basis.shells {
            let n_funcs = shell.n_basis_functions();
            let l = shell.l_value();
            let components = match cartesian_components(l) {
                Ok(c) => c,
                Err(_) => {
                    offset += n_funcs;
                    continue;
                }
            };

            let n_prims = shell.primitives.len();
            let n_comps = components.len();
            let mut norm_coef = Vec::with_capacity(n_prims);
            let mut exponents = Vec::with_capacity(n_prims);

            for prim in &shell.primitives {
                exponents.push(prim.exponent);
                let mut nc = Vec::with_capacity(n_comps);
                for comp in &components {
                    nc.push(prim.coefficient * cartesian_norm(prim.exponent, comp));
                }
                norm_coef.push(nc);
            }

            let alpha_min = shell.min_exponent();

            data.push(ShellPrecomputed {
                center: shell.center,
                angular_momentum: shell.angular_momentum,
                components,
                norm_coef,
                exponents,
                alpha_min,
                basis_offset: offset,
            });

            offset += n_funcs;
        }
        data
    };

    // ---- Evaluate density on grid ----
    let total_points = nx * ny * nz;
    let mut grid_values = vec![0.0f64; total_points];

    // Allocate chi vector once, reused for each grid point
    let mut chi = vec![0.0f64; nbf];

    let [x0, y0, z0] = grid_origin;

    for ix in 0..nx {
        let gx = x0 + ix as f64 * grid_spacing;

        for iy in 0..ny {
            let gy = y0 + iy as f64 * grid_spacing;

            for iz in 0..nz {
                let gz = z0 + iz as f64 * grid_spacing;

                // Zero the chi vector for this grid point
                for c in chi.iter_mut() {
                    *c = 0.0;
                }

                // Evaluate all basis functions at this grid point
                for sd in &shell_data {
                    let dx = gx - sd.center[0];
                    let dy = gy - sd.center[1];
                    let dz = gz - sd.center[2];
                    let dist_sq = dx * dx + dy * dy + dz * dz;

                    // Gaussian screening
                    if sd.alpha_min * dist_sq > GAUSSIAN_SCREENING_THRESHOLD {
                        continue;
                    }

                    let is_s = sd.angular_momentum == AngularMomentum::S;
                    let is_p = sd.angular_momentum == AngularMomentum::P;

                    if is_s {
                        // S-shell: single component, no angular part
                        let mut val = 0.0;
                        for (p, &alpha) in sd.exponents.iter().enumerate() {
                            val += sd.norm_coef[p][0] * (-alpha * dist_sq).exp();
                        }
                        chi[sd.basis_offset] = val;
                    } else if is_p {
                        // P-shell: all 3 components share the same radial sum
                        let mut radial = 0.0;
                        for (p, &alpha) in sd.exponents.iter().enumerate() {
                            radial += sd.norm_coef[p][0] * (-alpha * dist_sq).exp();
                        }
                        let disps = [dx, dy, dz];
                        for comp_idx in 0..3 {
                            chi[sd.basis_offset + comp_idx] = disps[comp_idx] * radial;
                        }
                    } else {
                        // General case (D-shells and above): per-component normalization
                        let disps = [dx, dy, dz];
                        for (comp_idx, powers) in sd.components.iter().enumerate() {
                            let angular = angular_factor(&disps, powers);
                            if angular.abs() < 1e-30 {
                                continue;
                            }
                            let mut radial_sum = 0.0;
                            for (p, &alpha) in sd.exponents.iter().enumerate() {
                                radial_sum += sd.norm_coef[p][comp_idx] * (-alpha * dist_sq).exp();
                            }
                            chi[sd.basis_offset + comp_idx] = angular * radial_sum;
                        }
                    }
                }

                // Contract chi with D using symmetry:
                // rho = sum_mu D[mu,mu] * chi[mu]^2
                //     + 2 * sum_{mu < nu} D[mu,nu] * chi[mu] * chi[nu]
                let mut rho = 0.0;
                for mu in 0..nbf {
                    let chi_mu = chi[mu];
                    if chi_mu.abs() < 1e-30 {
                        continue;
                    }
                    // Diagonal term
                    rho += density_matrix[mu * nbf + mu] * chi_mu * chi_mu;
                    // Off-diagonal terms (mu < nu)
                    for nu in (mu + 1)..nbf {
                        let chi_nu = chi[nu];
                        if chi_nu.abs() < 1e-30 {
                            continue;
                        }
                        rho += 2.0 * density_matrix[mu * nbf + nu] * chi_mu * chi_nu;
                    }
                }

                let flat_idx = ix * ny * nz + iy * nz + iz;
                grid_values[flat_idx] = rho;
            }
        }
    }

    // ---- Compute metadata ----
    let dv = grid_spacing.powi(3);
    let integrated_density: f64 = grid_values.iter().sum::<f64>() * dv;
    let max_density = grid_values.iter().cloned().fold(0.0f64, f64::max);

    Ok(DensityGridResult {
        values: grid_values,
        grid_dims,
        grid_origin,
        grid_spacing,
        integrated_density,
        n_electrons_expected: n_electrons,
        max_density,
        // Timing is handled by the caller (WASM layer uses js_sys::Date::now()
        // since std::time::Instant is not supported on wasm32-unknown-unknown)
        compute_time_ms: 0.0,
    })
}

// =============================================================================
// Internal Data Structures
// =============================================================================

/// Precomputed data for a single shell, avoiding recomputation in the inner loop.
struct ShellPrecomputed {
    center: [f64; 3],
    angular_momentum: AngularMomentum,
    components: Vec<CartesianPower>,
    /// norm_coef[prim_idx][comp_idx] = contraction_coeff * normalization
    norm_coef: Vec<Vec<f64>>,
    exponents: Vec<f64>,
    alpha_min: f64,
    basis_offset: usize,
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::{Atom, BasisSet};
    use approx::assert_abs_diff_eq;

    // =========================================================================
    // Helper: Build test basis sets
    // =========================================================================

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

    // PySCF 2.11.0 reference values
    //
    // H2 STO-3G (R=1.4 bohr), SCF energy = -1.11671432506255
    // Density matrix (includes factor of 2):
    // D = [[0.602657161419, 0.602657161419],
    //      [0.602657161419, 0.602657161419]]
    fn h2_density_matrix() -> Vec<f64> {
        vec![
            6.026571614189372e-01,
            6.026571614189370e-01,
            6.026571614189370e-01,
            6.026571614189368e-01,
        ]
    }

    // H2O STO-3G, SCF energy = -74.9630231335037
    // Density matrix (7x7, row-major):
    fn h2o_density_matrix() -> Vec<f64> {
        vec![
            2.106277830064741e+00,
            -4.461535227711436e-01,
            1.322835007549169e-17,
            1.629217466013592e-16,
            -1.085881938229459e-01,
            -2.835657164921020e-02,
            -2.835657164921002e-02,
            -4.461535227711436e-01,
            1.967805166046474e+00,
            -1.154397687301236e-16,
            -8.489433926532169e-16,
            6.175785252220366e-01,
            -3.431351566262539e-02,
            -3.431351566262618e-02,
            1.322835007549169e-17,
            -1.154397687301236e-16,
            2.000000000000000e+00,
            1.324013041190176e-16,
            1.209941783305419e-16,
            1.138339817945259e-16,
            1.606778345197832e-16,
            1.629217466013592e-16,
            -8.489433926532169e-16,
            1.324013041190176e-16,
            7.356315442784612e-01,
            -7.644540111146455e-16,
            5.398419093925200e-01,
            -5.398419093925191e-01,
            -1.085881938229459e-01,
            6.175785252220366e-01,
            1.209941783305419e-16,
            -7.644540111146455e-16,
            1.239425489277904e+00,
            -4.729157305763869e-01,
            -4.729157305763884e-01,
            -2.835657164921020e-02,
            -3.431351566262539e-02,
            1.138339817945259e-16,
            5.398419093925200e-01,
            -4.729157305763869e-01,
            6.012715980391495e-01,
            -1.910524652863499e-01,
            -2.835657164921002e-02,
            -3.431351566262618e-02,
            1.606778345197832e-16,
            -5.398419093925191e-01,
            -4.729157305763884e-01,
            -1.910524652863499e-01,
            6.012715980391500e-01,
        ]
    }

    /// Evaluate density at a single point by placing a tiny grid with the
    /// first grid point at the desired coordinate.
    fn density_at_point(
        density_matrix: &[f64],
        nbf: usize,
        n_electrons: usize,
        basis: &BasisSet,
        point: [f64; 3],
    ) -> f64 {
        let result = evaluate_density_on_grid(
            density_matrix,
            nbf,
            n_electrons,
            basis,
            point,
            0.001, // very small spacing so other grid points are nearby
            [2, 2, 2],
        )
        .unwrap();
        result.values[0]
    }

    // =========================================================================
    // T1: Input Validation Tests
    // =========================================================================

    #[test]
    fn test_wrong_density_matrix_size() {
        let basis = h2_basis();
        let wrong_d = vec![1.0, 2.0, 3.0]; // should be 4 (2x2)
        let result = evaluate_density_on_grid(&wrong_d, 2, 2, &basis, [0.0; 3], 1.0, [3, 3, 3]);
        assert!(matches!(
            result,
            Err(OrbitalError::DensityMatrixSizeMismatch { .. })
        ));
    }

    #[test]
    fn test_invalid_grid_dims() {
        let basis = h2_basis();
        let d = h2_density_matrix();
        let result = evaluate_density_on_grid(&d, 2, 2, &basis, [0.0; 3], 1.0, [1, 3, 3]);
        assert!(matches!(result, Err(OrbitalError::InvalidGridDimension(_))));
    }

    #[test]
    fn test_zero_spacing() {
        let basis = h2_basis();
        let d = h2_density_matrix();
        let result = evaluate_density_on_grid(&d, 2, 2, &basis, [0.0; 3], 0.0, [3, 3, 3]);
        assert!(matches!(result, Err(OrbitalError::InvalidGridSpacing(_))));
    }

    #[test]
    fn test_empty_basis() {
        // Build an empty basis set manually
        let empty_basis = BasisSet {
            atoms: vec![],
            shells: vec![],
            name: String::new(),
            n_basis: 0,
            n_electrons: 0,
            nuclear_repulsion: 0.0,
        };
        let d = vec![0.0; 0];
        let result = evaluate_density_on_grid(&d, 0, 0, &empty_basis, [0.0; 3], 1.0, [3, 3, 3]);
        assert!(matches!(result, Err(OrbitalError::EmptyBasis)));
    }

    // =========================================================================
    // T2: H2 Golden Value Tests (PySCF 2.11.0 reference)
    // =========================================================================

    #[test]
    fn test_h2_density_at_atom1() {
        // PySCF: rho at (0,0,0) = 3.548920730810347e-01
        let basis = h2_basis();
        let d = h2_density_matrix();
        let rho = density_at_point(&d, 2, 2, &basis, [0.0, 0.0, 0.0]);
        assert_abs_diff_eq!(rho, 3.548920730810347e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2_density_at_midpoint() {
        // PySCF: rho at (0,0,0.7) = 2.559266712987876e-01
        let basis = h2_basis();
        let d = h2_density_matrix();
        let rho = density_at_point(&d, 2, 2, &basis, [0.0, 0.0, 0.7]);
        assert_abs_diff_eq!(rho, 2.559266712987876e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2_density_at_atom2() {
        // PySCF: rho at (0,0,1.4) = 3.548920730810345e-01
        let basis = h2_basis();
        let d = h2_density_matrix();
        let rho = density_at_point(&d, 2, 2, &basis, [0.0, 0.0, 1.4]);
        assert_abs_diff_eq!(rho, 3.548920730810345e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2_density_off_axis() {
        // PySCF: rho at (1,0,0.7) = 7.156822072208588e-02
        let basis = h2_basis();
        let d = h2_density_matrix();
        let rho = density_at_point(&d, 2, 2, &basis, [1.0, 0.0, 0.7]);
        assert_abs_diff_eq!(rho, 7.156822072208588e-02, epsilon = 1e-8);
    }

    // =========================================================================
    // T3: H2O Golden Value Tests (PySCF 2.11.0 reference)
    // =========================================================================

    #[test]
    fn test_h2o_density_near_oxygen() {
        // PySCF: rho at (0,0,0.221665) = 1.933138767520558e+02 (near O nucleus)
        let basis = h2o_basis();
        let d = h2o_density_matrix();
        let rho = density_at_point(&d, 7, 10, &basis, [0.0, 0.0, 0.221665]);
        assert_abs_diff_eq!(rho, 1.933138767520558e+02, epsilon = 1e-6);
    }

    #[test]
    fn test_h2o_density_near_h1() {
        // PySCF: rho at (0,1.430901,-0.886659) = 3.627094851968106e-01
        let basis = h2o_basis();
        let d = h2o_density_matrix();
        let rho = density_at_point(&d, 7, 10, &basis, [0.0, 1.430901, -0.886659]);
        assert_abs_diff_eq!(rho, 3.627094851968106e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2o_density_off_axis() {
        // PySCF: rho at (1,0,0) = 4.911442048926800e-01
        let basis = h2o_basis();
        let d = h2o_density_matrix();
        let rho = density_at_point(&d, 7, 10, &basis, [1.0, 0.0, 0.0]);
        assert_abs_diff_eq!(rho, 4.911442048926800e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_h2o_density_general_point() {
        // PySCF: rho at (0.5,0.5,0.5) = 7.637351101341471e-01
        let basis = h2o_basis();
        let d = h2o_density_matrix();
        let rho = density_at_point(&d, 7, 10, &basis, [0.5, 0.5, 0.5]);
        assert_abs_diff_eq!(rho, 7.637351101341471e-01, epsilon = 1e-8);
    }

    // =========================================================================
    // T4: Integrated Density Tests
    // =========================================================================

    #[test]
    fn test_h2_integrated_density() {
        // H2 has 2 electrons. The integrated density on a sufficiently large
        // grid should approximate 2.0 within 5% grid quadrature error.
        let basis = h2_basis();
        let d = h2_density_matrix();

        let (min_corner, max_corner) = super::super::grid::auto_grid_bounds(&basis, 4.0);
        let spacing = 0.3;
        let nx = ((max_corner[0] - min_corner[0]) / spacing).ceil() as usize + 1;
        let ny = ((max_corner[1] - min_corner[1]) / spacing).ceil() as usize + 1;
        let nz = ((max_corner[2] - min_corner[2]) / spacing).ceil() as usize + 1;

        let result =
            evaluate_density_on_grid(&d, 2, 2, &basis, min_corner, spacing, [nx, ny, nz]).unwrap();

        let n_expected = 2.0;
        let rel_error = (result.integrated_density - n_expected).abs() / n_expected;
        assert!(
            rel_error < 0.05,
            "H2 integrated density = {:.4}, expected ~2.0, relative error = {:.4}",
            result.integrated_density,
            rel_error
        );
        assert_eq!(result.n_electrons_expected, 2);
    }

    #[test]
    fn test_h2o_integrated_density() {
        // H2O has 10 electrons. The integrated density on a sufficiently large
        // grid should approximate 10.0 within 5% grid quadrature error.
        let basis = h2o_basis();
        let d = h2o_density_matrix();

        let (min_corner, max_corner) = super::super::grid::auto_grid_bounds(&basis, 4.0);
        let spacing = 0.3;
        let nx = ((max_corner[0] - min_corner[0]) / spacing).ceil() as usize + 1;
        let ny = ((max_corner[1] - min_corner[1]) / spacing).ceil() as usize + 1;
        let nz = ((max_corner[2] - min_corner[2]) / spacing).ceil() as usize + 1;

        let result =
            evaluate_density_on_grid(&d, 7, 10, &basis, min_corner, spacing, [nx, ny, nz]).unwrap();

        let n_expected = 10.0;
        let rel_error = (result.integrated_density - n_expected).abs() / n_expected;
        assert!(
            rel_error < 0.05,
            "H2O integrated density = {:.4}, expected ~10.0, relative error = {:.4}",
            result.integrated_density,
            rel_error
        );
        assert_eq!(result.n_electrons_expected, 10);
    }

    // =========================================================================
    // Properties Tests
    // =========================================================================

    #[test]
    fn test_density_non_negativity() {
        // Electron density must be non-negative everywhere.
        let basis = h2_basis();
        let d = h2_density_matrix();

        let result =
            evaluate_density_on_grid(&d, 2, 2, &basis, [-3.0, -3.0, -3.0], 0.5, [13, 13, 15])
                .unwrap();

        for (i, &rho) in result.values.iter().enumerate() {
            assert!(
                rho >= -1e-14,
                "Density should be non-negative, got {rho} at index {i}"
            );
        }
    }

    #[test]
    fn test_density_decay_at_distance() {
        // Density should decrease as we move far from the molecule.
        let basis = h2_basis();
        let d = h2_density_matrix();

        let rho_near = density_at_point(&d, 2, 2, &basis, [0.0, 0.0, 0.7]);
        let rho_far = density_at_point(&d, 2, 2, &basis, [5.0, 0.0, 0.0]);

        assert!(
            rho_near > rho_far,
            "Density near molecule ({rho_near}) should exceed density far away ({rho_far})"
        );
    }

    #[test]
    fn test_h2_symmetry_about_midpoint() {
        // H2 density should be symmetric about the bond midpoint (0,0,0.7).
        // rho(0,0,0) should equal rho(0,0,1.4) (by symmetry of the two H atoms).
        let basis = h2_basis();
        let d = h2_density_matrix();

        let rho_h1 = density_at_point(&d, 2, 2, &basis, [0.0, 0.0, 0.0]);
        let rho_h2 = density_at_point(&d, 2, 2, &basis, [0.0, 0.0, 1.4]);

        assert_abs_diff_eq!(rho_h1, rho_h2, epsilon = 1e-10);
    }

    // =========================================================================
    // Grid Metadata Tests
    // =========================================================================

    #[test]
    fn test_grid_metadata_correct() {
        let basis = h2_basis();
        let d = h2_density_matrix();

        let result =
            evaluate_density_on_grid(&d, 2, 2, &basis, [-2.0, -2.0, -2.0], 1.0, [5, 5, 7]).unwrap();

        assert_eq!(result.values.len(), 5 * 5 * 7);
        assert_eq!(result.grid_dims, [5, 5, 7]);
        assert_eq!(result.grid_origin, [-2.0, -2.0, -2.0]);
        assert_abs_diff_eq!(result.grid_spacing, 1.0, epsilon = 1e-15);
        assert_eq!(result.n_electrons_expected, 2);
        assert!(result.max_density > 0.0);
        assert!(result.integrated_density > 0.0);
        assert!(result.compute_time_ms >= 0.0);
    }
}
