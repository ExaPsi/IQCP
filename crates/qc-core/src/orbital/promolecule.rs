//! Promolecule density and difference density computation
//!
//! This module provides the infrastructure for computing the promolecule
//! density (superposition of spherically-averaged isolated atomic densities)
//! and the difference density (deformation density):
//!
//! ```text
//! Delta-rho(r) = rho_molecule(r) - rho_promolecule(r)
//! ```
//!
//! where:
//! - `rho_molecule(r)` is the converged SCF electron density
//! - `rho_promolecule(r) = sum_A rho_atom_A(|r - R_A|)` is the superposition
//!   of spherically-averaged atomic densities placed at the nuclear positions
//!
//! # Physical Interpretation
//!
//! - **Delta-rho > 0** (accumulation): electrons move INTO this region during
//!   bond formation. For H2, this is the bonding region between the nuclei.
//! - **Delta-rho < 0** (depletion): electrons move OUT of this region. For H2,
//!   this is the region behind each nucleus.
//!
//! # Data Source
//!
//! Atomic density radial profiles are pre-tabulated from UHF/STO-3G calculations
//! using PySCF 2.11.0. The profiles are embedded at compile time from JSON files
//! in `content/presets/atomic_densities/`. Each profile contains 1001 points on a
//! uniform radial grid from r=0 to r=10 bohr with spacing 0.01 bohr.
//!
//! # References
//!
//! - Coppens (1997). "X-Ray Charge Densities and Chemical Bonding." Oxford.
//! - Phase 3 TDD Section 8.3, FR-DEN-04.

use serde::{Deserialize, Serialize};
use std::sync::OnceLock;

use super::grid::OrbitalError;

// =============================================================================
// Embedded Atomic Density Data (compile-time)
// =============================================================================

// Include JSON data for all 18 elements (H through Ar) at compile time.
// Total size: ~18 * 1001 * ~12 bytes/value = ~216 KB of JSON strings.
// Parsed on first access via OnceLock.
static H_DATA: &str = include_str!("../../../../content/presets/atomic_densities/h.json");
static HE_DATA: &str = include_str!("../../../../content/presets/atomic_densities/he.json");
static LI_DATA: &str = include_str!("../../../../content/presets/atomic_densities/li.json");
static BE_DATA: &str = include_str!("../../../../content/presets/atomic_densities/be.json");
static B_DATA: &str = include_str!("../../../../content/presets/atomic_densities/b.json");
static C_DATA: &str = include_str!("../../../../content/presets/atomic_densities/c.json");
static N_DATA: &str = include_str!("../../../../content/presets/atomic_densities/n.json");
static O_DATA: &str = include_str!("../../../../content/presets/atomic_densities/o.json");
static F_DATA: &str = include_str!("../../../../content/presets/atomic_densities/f.json");
static NE_DATA: &str = include_str!("../../../../content/presets/atomic_densities/ne.json");
static NA_DATA: &str = include_str!("../../../../content/presets/atomic_densities/na.json");
static MG_DATA: &str = include_str!("../../../../content/presets/atomic_densities/mg.json");
static AL_DATA: &str = include_str!("../../../../content/presets/atomic_densities/al.json");
static SI_DATA: &str = include_str!("../../../../content/presets/atomic_densities/si.json");
static P_DATA: &str = include_str!("../../../../content/presets/atomic_densities/p.json");
static S_DATA: &str = include_str!("../../../../content/presets/atomic_densities/s.json");
static CL_DATA: &str = include_str!("../../../../content/presets/atomic_densities/cl.json");
static AR_DATA: &str = include_str!("../../../../content/presets/atomic_densities/ar.json");

/// Returns the JSON string for element with atomic number Z.
fn json_for_z(z: u32) -> Option<&'static str> {
    match z {
        1 => Some(H_DATA),
        2 => Some(HE_DATA),
        3 => Some(LI_DATA),
        4 => Some(BE_DATA),
        5 => Some(B_DATA),
        6 => Some(C_DATA),
        7 => Some(N_DATA),
        8 => Some(O_DATA),
        9 => Some(F_DATA),
        10 => Some(NE_DATA),
        11 => Some(NA_DATA),
        12 => Some(MG_DATA),
        13 => Some(AL_DATA),
        14 => Some(SI_DATA),
        15 => Some(P_DATA),
        16 => Some(S_DATA),
        17 => Some(CL_DATA),
        18 => Some(AR_DATA),
        _ => None,
    }
}

// =============================================================================
// JSON Deserialization Format
// =============================================================================

/// Raw JSON format matching the files in `content/presets/atomic_densities/`.
///
/// JSON uses snake_case keys: atomic_number, n_electrons, density_values, etc.
#[derive(Deserialize)]
struct AtomicDensityJson {
    #[allow(dead_code)]
    element: String,
    atomic_number: u32,
    n_electrons: u32,
    #[allow(dead_code)]
    method: String,
    #[allow(dead_code)]
    basis: String,
    radial_grid: RadialGridJson,
    density_values: Vec<f64>,
}

#[derive(Deserialize)]
struct RadialGridJson {
    #[allow(dead_code)]
    r_min: f64,
    r_max: f64,
    spacing: f64,
    #[allow(dead_code)]
    n_points: usize,
}

// =============================================================================
// Atomic Density Profile
// =============================================================================

/// Pre-tabulated spherically-averaged atomic density radial profile.
///
/// The density is stored as rho(r) on a uniform radial grid from r=0
/// to r=r_max. At distances beyond r_max, the density is zero.
///
/// These profiles are generated from UHF/STO-3G calculations using PySCF 2.11.0
/// with 6-point octahedral angular averaging for spherical symmetry.
///
/// Reference: Phase 3 TDD Section 8.3, FR-DEN-04.
#[derive(Debug, Clone)]
pub struct AtomicDensityProfile {
    /// Atomic number
    pub atomic_number: u32,
    /// Number of electrons in the atom
    pub n_electrons: u32,
    /// Maximum radius of the radial grid (bohr)
    pub r_max: f64,
    /// Uniform spacing of the radial grid (bohr)
    pub r_spacing: f64,
    /// Radial density values rho(r) at r = 0, r_spacing, 2*r_spacing, ...
    pub density_values: Vec<f64>,
}

impl AtomicDensityProfile {
    /// Evaluate the atomic density at distance r from the nucleus
    /// using linear interpolation of the radial profile.
    ///
    /// Returns 0.0 for r >= r_max or r < 0 (density negligible beyond
    /// the tabulated range of 10 bohr).
    ///
    /// # Linear Interpolation
    ///
    /// Given uniform grid spacing h and distance r:
    /// ```text
    /// k = floor(r / h)           // integer index
    /// t = r / h - k              // fractional part, t in [0, 1)
    /// rho(r) = (1-t) * rho[k] + t * rho[k+1]
    /// ```
    ///
    /// With 0.01 bohr spacing, the interpolation error is negligible for
    /// the smooth STO-3G radial density functions (~1e-6 relative error).
    #[inline]
    pub fn evaluate(&self, r: f64) -> f64 {
        if r < 0.0 || r >= self.r_max || self.density_values.is_empty() {
            return 0.0;
        }

        let t = r / self.r_spacing;
        let k = t as usize;
        let frac = t - k as f64;

        if k + 1 >= self.density_values.len() {
            // At or beyond last tabulated point
            return *self.density_values.last().unwrap_or(&0.0);
        }

        // Linear interpolation
        (1.0 - frac) * self.density_values[k] + frac * self.density_values[k + 1]
    }
}

// =============================================================================
// Global Profile Cache (parsed once from embedded JSON)
// =============================================================================

/// Global storage for all parsed atomic density profiles.
/// Parsed on first access from compile-time embedded JSON strings.
static PROFILES: OnceLock<Vec<Option<AtomicDensityProfile>>> = OnceLock::new();

/// Initialize and return the global profile cache.
fn get_profiles() -> &'static Vec<Option<AtomicDensityProfile>> {
    PROFILES.get_or_init(|| {
        let mut profiles: Vec<Option<AtomicDensityProfile>> = Vec::with_capacity(19);
        // Index 0 is unused (no element with Z=0)
        profiles.push(None);

        for z in 1..=18u32 {
            let profile = json_for_z(z).and_then(|json_str| {
                let raw: AtomicDensityJson = serde_json::from_str(json_str).ok()?;
                Some(AtomicDensityProfile {
                    atomic_number: raw.atomic_number,
                    n_electrons: raw.n_electrons,
                    r_max: raw.radial_grid.r_max,
                    r_spacing: raw.radial_grid.spacing,
                    density_values: raw.density_values,
                })
            });
            profiles.push(profile);
        }

        profiles
    })
}

/// Get the atomic density profile for element with atomic number Z.
///
/// Returns `None` if Z < 1 or Z > 18 (only H through Ar are available).
/// The profile is lazily parsed from embedded JSON on first access.
pub fn get_atomic_density(z: u32) -> Option<&'static AtomicDensityProfile> {
    let profiles = get_profiles();
    profiles.get(z as usize)?.as_ref()
}

/// Check if all atoms in a system have available atomic density profiles.
///
/// Returns a vector of atomic numbers that are NOT supported (empty = all supported).
pub fn unsupported_elements(atomic_numbers: &[u32]) -> Vec<u32> {
    atomic_numbers
        .iter()
        .copied()
        .filter(|&z| get_atomic_density(z).is_none())
        .collect()
}

// =============================================================================
// Promolecule Grid Evaluation
// =============================================================================

/// Evaluate the promolecule density on a uniform 3D grid.
///
/// The promolecule is the superposition of spherically-averaged atomic densities:
///
/// ```text
/// rho_promolecule(r) = sum_A rho_atom_A(|r - R_A|)
/// ```
///
/// Each atomic density is evaluated by linear interpolation of the pre-tabulated
/// radial profile.
///
/// # Arguments
///
/// * `atoms` - Atom specifications: `(Z, [x, y, z])` with positions in Bohr
/// * `grid_origin` - Grid origin `[x, y, z]` in Bohr
/// * `grid_spacing` - Uniform grid spacing in Bohr
/// * `grid_dims` - Grid dimensions `[nx, ny, nz]`
///
/// # Grid Indexing Convention
///
/// C-order (x-slowest, z-fastest):
/// ```text
/// index = ix * ny * nz + iy * nz + iz
/// ```
///
/// # Returns
///
/// Flat vector of promolecule density values (same indexing as molecular density grid).
///
/// # Errors
///
/// Returns error if any atom has an unsupported atomic number (Z > 18 or Z < 1),
/// if grid dimensions are invalid, or if grid spacing is non-positive.
pub fn evaluate_promolecule_on_grid(
    atoms: &[(u32, [f64; 3])],
    grid_origin: [f64; 3],
    grid_spacing: f64,
    grid_dims: [usize; 3],
) -> Result<Vec<f64>, OrbitalError> {
    let [nx, ny, nz] = grid_dims;

    // Input validation
    if nx < 2 || ny < 2 || nz < 2 {
        return Err(OrbitalError::InvalidGridDimension(nx.min(ny).min(nz)));
    }
    if grid_spacing <= 0.0 {
        return Err(OrbitalError::InvalidGridSpacing(grid_spacing));
    }

    // Look up atomic density profiles for all atoms
    let profiles: Vec<(&AtomicDensityProfile, [f64; 3])> = atoms
        .iter()
        .map(|&(z, pos)| {
            get_atomic_density(z)
                .ok_or(OrbitalError::UnsupportedElement(z))
                .map(|p| (p, pos))
        })
        .collect::<Result<Vec<_>, _>>()?;

    // Evaluate promolecule density on the grid
    let total_points = nx * ny * nz;
    let mut grid_values = vec![0.0f64; total_points];

    let [x0, y0, z0] = grid_origin;

    for ix in 0..nx {
        let gx = x0 + ix as f64 * grid_spacing;

        for iy in 0..ny {
            let gy = y0 + iy as f64 * grid_spacing;

            for iz in 0..nz {
                let gz = z0 + iz as f64 * grid_spacing;

                let flat_idx = ix * ny * nz + iy * nz + iz;

                // Sum contributions from all atoms
                let mut rho_pro = 0.0;
                for &(profile, center) in &profiles {
                    let dx = gx - center[0];
                    let dy = gy - center[1];
                    let dz = gz - center[2];
                    let r = (dx * dx + dy * dy + dz * dz).sqrt();
                    rho_pro += profile.evaluate(r);
                }

                grid_values[flat_idx] = rho_pro;
            }
        }
    }

    Ok(grid_values)
}

// =============================================================================
// Difference Density
// =============================================================================

/// Result of difference density computation.
///
/// Contains the difference density grid values and summary statistics
/// useful for validation and visualization.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct DifferenceDensityResult {
    /// Flat Delta-rho values at each grid point (C-order: x-slowest, z-fastest)
    pub values: Vec<f64>,
    /// Integrated Delta-rho: sum(Delta-rho) * dV.
    /// Should be approximately 0 for density conservation.
    pub integrated_delta_rho: f64,
    /// Maximum positive Delta-rho value (electron accumulation peak)
    pub max_accumulation: f64,
    /// Maximum negative Delta-rho value (electron depletion trough, stored as negative)
    pub max_depletion: f64,
}

/// Compute difference density: Delta-rho = rho_molecule - rho_promolecule.
///
/// The difference density (deformation density) reveals where electrons
/// redistribute when isolated atoms form a molecule:
///
/// ```text
/// Delta-rho(r) = rho_molecule(r) - rho_promolecule(r)
/// ```
///
/// # Density Conservation
///
/// The integral of Delta-rho over all space must be zero (same total
/// number of electrons in molecule and promolecule). On a finite grid:
///
/// ```text
/// |sum_i Delta-rho(r_i) * dV| / N_electrons < 0.05
/// ```
///
/// # Arguments
///
/// * `molecular_density` - Total molecular density grid values from
///   `evaluate_density_on_grid()` (C-order flat array)
/// * `promolecule_density` - Promolecule density grid values from
///   `evaluate_promolecule_on_grid()` (same grid specification)
/// * `grid_spacing` - Uniform grid spacing in Bohr (for integration)
///
/// # Errors
///
/// Returns error if the two density grids have different lengths.
pub fn compute_difference_density(
    molecular_density: &[f64],
    promolecule_density: &[f64],
    grid_spacing: f64,
) -> Result<DifferenceDensityResult, OrbitalError> {
    if molecular_density.len() != promolecule_density.len() {
        return Err(OrbitalError::GridSizeMismatch {
            expected: molecular_density.len(),
            actual: promolecule_density.len(),
        });
    }

    let dv = grid_spacing.powi(3);

    let mut values = Vec::with_capacity(molecular_density.len());
    let mut sum = 0.0f64;
    let mut max_accumulation = 0.0f64;
    let mut max_depletion = 0.0f64;

    for (&rho_mol, &rho_pro) in molecular_density.iter().zip(promolecule_density.iter()) {
        let delta = rho_mol - rho_pro;
        values.push(delta);
        sum += delta;

        if delta > max_accumulation {
            max_accumulation = delta;
        }
        if delta < max_depletion {
            max_depletion = delta;
        }
    }

    let integrated_delta_rho = sum * dv;

    Ok(DifferenceDensityResult {
        values,
        integrated_delta_rho,
        max_accumulation,
        max_depletion,
    })
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    // =========================================================================
    // AtomicDensityProfile unit tests
    // =========================================================================

    fn make_test_profile() -> AtomicDensityProfile {
        AtomicDensityProfile {
            atomic_number: 99,
            n_electrons: 1,
            r_max: 3.0,
            r_spacing: 1.0,
            density_values: vec![1.0, 0.5, 0.2, 0.05],
        }
    }

    #[test]
    fn test_evaluate_at_grid_points() {
        let p = make_test_profile();
        assert_abs_diff_eq!(p.evaluate(0.0), 1.0, epsilon = 1e-15);
        assert_abs_diff_eq!(p.evaluate(1.0), 0.5, epsilon = 1e-15);
        assert_abs_diff_eq!(p.evaluate(2.0), 0.2, epsilon = 1e-15);
    }

    #[test]
    fn test_linear_interpolation_midpoints() {
        let p = make_test_profile();
        // Between 0 and 1: (1.0 + 0.5) / 2 = 0.75
        assert_abs_diff_eq!(p.evaluate(0.5), 0.75, epsilon = 1e-15);
        // Between 1 and 2: (0.5 + 0.2) / 2 = 0.35
        assert_abs_diff_eq!(p.evaluate(1.5), 0.35, epsilon = 1e-15);
        // Between 2 and 3: (0.2 + 0.05) / 2 = 0.125
        assert_abs_diff_eq!(p.evaluate(2.5), 0.125, epsilon = 1e-15);
    }

    #[test]
    fn test_interpolation_quarter_points() {
        let p = make_test_profile();
        // At r=0.25: (1-0.25)*1.0 + 0.25*0.5 = 0.875
        assert_abs_diff_eq!(p.evaluate(0.25), 0.875, epsilon = 1e-15);
        // At r=0.75: (1-0.75)*1.0 + 0.75*0.5 = 0.625
        assert_abs_diff_eq!(p.evaluate(0.75), 0.625, epsilon = 1e-15);
    }

    #[test]
    fn test_beyond_rmax_returns_zero() {
        let p = make_test_profile();
        assert_eq!(p.evaluate(3.0), 0.0);
        assert_eq!(p.evaluate(5.0), 0.0);
        assert_eq!(p.evaluate(100.0), 0.0);
    }

    #[test]
    fn test_negative_r_returns_zero() {
        let p = make_test_profile();
        assert_eq!(p.evaluate(-0.1), 0.0);
        assert_eq!(p.evaluate(-5.0), 0.0);
    }

    #[test]
    fn test_empty_profile_returns_zero() {
        let p = AtomicDensityProfile {
            atomic_number: 1,
            n_electrons: 1,
            r_max: 10.0,
            r_spacing: 0.01,
            density_values: vec![],
        };
        assert_eq!(p.evaluate(0.0), 0.0);
        assert_eq!(p.evaluate(1.0), 0.0);
    }

    // =========================================================================
    // Embedded atomic density data tests
    // =========================================================================

    #[test]
    fn test_all_18_elements_load() {
        for z in 1..=18u32 {
            let profile = get_atomic_density(z);
            assert!(
                profile.is_some(),
                "Atomic density profile for Z={z} should be available"
            );
            let p = profile.unwrap();
            assert_eq!(p.atomic_number, z);
            assert_eq!(p.density_values.len(), 1001);
            assert_abs_diff_eq!(p.r_max, 10.0, epsilon = 1e-10);
            assert_abs_diff_eq!(p.r_spacing, 0.01, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_z_out_of_range_returns_none() {
        assert!(get_atomic_density(0).is_none());
        assert!(get_atomic_density(19).is_none());
        assert!(get_atomic_density(100).is_none());
    }

    #[test]
    fn test_hydrogen_density_at_origin() {
        // PySCF 2.11.0 UHF/STO-3G: rho_H(0) = 3.946941451474386e-01
        let h = get_atomic_density(1).unwrap();
        assert_abs_diff_eq!(h.evaluate(0.0), 3.946941451474e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_hydrogen_density_at_0p7() {
        // PySCF 2.11.0 UHF/STO-3G: rho_H(0.7) = 1.061659462803928e-01
        let h = get_atomic_density(1).unwrap();
        assert_abs_diff_eq!(h.evaluate(0.7), 1.061659462803928e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_hydrogen_density_at_1p0() {
        // PySCF 2.11.0 UHF/STO-3G: rho_H(1.0) = 4.974498017996373e-02
        let h = get_atomic_density(1).unwrap();
        assert_abs_diff_eq!(h.evaluate(1.0), 4.974498017996e-02, epsilon = 1e-8);
    }

    #[test]
    fn test_oxygen_density_at_origin() {
        // PySCF 2.11.0 UHF/STO-3G: rho_O(0) = 194.230308933455
        let o = get_atomic_density(8).unwrap();
        assert_abs_diff_eq!(o.evaluate(0.0), 194.230308933, epsilon = 1e-3);
    }

    #[test]
    fn test_helium_density_at_origin() {
        // He has 2 electrons, density at origin should be ~2.0
        let he = get_atomic_density(2).unwrap();
        assert!(
            he.evaluate(0.0) > 1.0,
            "He density at origin should be > 1.0"
        );
        assert_eq!(he.n_electrons, 2);
    }

    #[test]
    fn test_electron_counts_correct() {
        let expected: [(u32, u32); 18] = [
            (1, 1),
            (2, 2),
            (3, 3),
            (4, 4),
            (5, 5),
            (6, 6),
            (7, 7),
            (8, 8),
            (9, 9),
            (10, 10),
            (11, 11),
            (12, 12),
            (13, 13),
            (14, 14),
            (15, 15),
            (16, 16),
            (17, 17),
            (18, 18),
        ];
        for (z, n_e) in expected {
            let p = get_atomic_density(z).unwrap();
            assert_eq!(
                p.n_electrons, n_e,
                "Element Z={z} should have {n_e} electrons"
            );
        }
    }

    #[test]
    fn test_density_non_negative() {
        // All atomic densities must be non-negative
        for z in 1..=18u32 {
            let p = get_atomic_density(z).unwrap();
            for (i, &rho) in p.density_values.iter().enumerate() {
                assert!(
                    rho >= -1e-14,
                    "Z={z}: density at index {i} should be non-negative, got {rho}"
                );
            }
        }
    }

    #[test]
    fn test_density_monotonically_decreasing_far_from_nucleus() {
        // Beyond ~1 bohr, atomic densities should generally decrease monotonically.
        // We test from index 200 (r=2.0) onward where all atoms are smooth.
        for z in 1..=18u32 {
            let p = get_atomic_density(z).unwrap();
            for i in 200..p.density_values.len() - 1 {
                assert!(
                    p.density_values[i] >= p.density_values[i + 1] - 1e-14,
                    "Z={z}: density should decrease at large r, but index {i} ({}) < index {} ({})",
                    p.density_values[i],
                    i + 1,
                    p.density_values[i + 1]
                );
            }
        }
    }

    // =========================================================================
    // unsupported_elements tests
    // =========================================================================

    #[test]
    fn test_unsupported_elements_empty_for_valid() {
        assert!(unsupported_elements(&[1, 6, 8]).is_empty());
    }

    #[test]
    fn test_unsupported_elements_detects_invalid() {
        let unsup = unsupported_elements(&[1, 26, 8, 92]);
        assert_eq!(unsup, vec![26, 92]);
    }

    // =========================================================================
    // Promolecule grid evaluation tests
    // =========================================================================

    #[test]
    fn test_promolecule_single_h_at_origin() {
        // Single H at origin: promolecule = H atomic density
        // At (0,0,0): rho_pro = rho_H(0) = 0.3947
        let atoms = vec![(1u32, [0.0, 0.0, 0.0])];
        let grid = evaluate_promolecule_on_grid(&atoms, [0.0, 0.0, 0.0], 0.001, [2, 2, 2]).unwrap();
        // grid[0] is at (0,0,0)
        assert_abs_diff_eq!(grid[0], 3.946941451474e-01, epsilon = 1e-8);
    }

    #[test]
    fn test_promolecule_h2_at_midpoint() {
        // H2 at (0,0,0) and (0,0,1.4): promolecule at midpoint (0,0,0.7)
        // rho_pro(midpoint) = rho_H(0.7) + rho_H(0.7) = 2 * 0.10617
        // PySCF: rho_pro at midpoint = 2.123318925607856e-01
        let atoms = vec![(1u32, [0.0, 0.0, 0.0]), (1u32, [0.0, 0.0, 1.4])];
        let grid = evaluate_promolecule_on_grid(&atoms, [0.0, 0.0, 0.7], 0.001, [2, 2, 2]).unwrap();
        assert_abs_diff_eq!(grid[0], 2.123318925607856e-01, epsilon = 1e-6);
    }

    #[test]
    fn test_promolecule_h2_symmetry() {
        // rho_pro at (0,0,0) should equal rho_pro at (0,0,1.4) by symmetry
        let atoms = vec![(1u32, [0.0, 0.0, 0.0]), (1u32, [0.0, 0.0, 1.4])];

        let grid_at_0 =
            evaluate_promolecule_on_grid(&atoms, [0.0, 0.0, 0.0], 0.001, [2, 2, 2]).unwrap();
        let grid_at_14 =
            evaluate_promolecule_on_grid(&atoms, [0.0, 0.0, 1.4], 0.001, [2, 2, 2]).unwrap();

        assert_abs_diff_eq!(grid_at_0[0], grid_at_14[0], epsilon = 1e-10);
    }

    #[test]
    fn test_promolecule_h2_at_nuclei() {
        // PySCF: rho_pro at H1 = 4.140534049811135e-01
        // = rho_H(0.0) + rho_H(1.4) = 0.3947 + 0.0194 = 0.4141
        let atoms = vec![(1u32, [0.0, 0.0, 0.0]), (1u32, [0.0, 0.0, 1.4])];
        let grid = evaluate_promolecule_on_grid(&atoms, [0.0, 0.0, 0.0], 0.001, [2, 2, 2]).unwrap();
        assert_abs_diff_eq!(grid[0], 4.140534049811135e-01, epsilon = 1e-6);
    }

    #[test]
    fn test_promolecule_invalid_element() {
        let atoms = vec![(1u32, [0.0, 0.0, 0.0]), (26u32, [0.0, 0.0, 1.0])];
        let result = evaluate_promolecule_on_grid(&atoms, [0.0, 0.0, 0.0], 1.0, [3, 3, 3]);
        assert!(result.is_err());
    }

    #[test]
    fn test_promolecule_invalid_grid_dims() {
        let atoms = vec![(1u32, [0.0, 0.0, 0.0])];
        let result = evaluate_promolecule_on_grid(&atoms, [0.0, 0.0, 0.0], 1.0, [1, 3, 3]);
        assert!(matches!(result, Err(OrbitalError::InvalidGridDimension(_))));
    }

    #[test]
    fn test_promolecule_invalid_spacing() {
        let atoms = vec![(1u32, [0.0, 0.0, 0.0])];
        let result = evaluate_promolecule_on_grid(&atoms, [0.0, 0.0, 0.0], 0.0, [3, 3, 3]);
        assert!(matches!(result, Err(OrbitalError::InvalidGridSpacing(_))));
    }

    // =========================================================================
    // Difference density tests
    // =========================================================================

    #[test]
    fn test_difference_density_identical_grids() {
        // If mol == pro, delta should be all zeros
        let mol = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let result = compute_difference_density(&mol, &mol, 1.0).unwrap();
        for &v in &result.values {
            assert_abs_diff_eq!(v, 0.0, epsilon = 1e-15);
        }
        assert_abs_diff_eq!(result.integrated_delta_rho, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(result.max_accumulation, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(result.max_depletion, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_difference_density_size_mismatch() {
        let mol = vec![1.0, 2.0, 3.0];
        let pro = vec![1.0, 2.0];
        let result = compute_difference_density(&mol, &pro, 1.0);
        assert!(matches!(result, Err(OrbitalError::GridSizeMismatch { .. })));
    }

    #[test]
    fn test_difference_density_accumulation_depletion() {
        // mol > pro at midpoint (accumulation), mol < pro at edges (depletion)
        let mol = vec![0.5, 1.0, 0.5];
        let pro = vec![0.8, 0.6, 0.8];
        let result = compute_difference_density(&mol, &pro, 1.0).unwrap();

        assert_abs_diff_eq!(result.values[0], -0.3, epsilon = 1e-15); // depletion
        assert_abs_diff_eq!(result.values[1], 0.4, epsilon = 1e-15); // accumulation
        assert_abs_diff_eq!(result.values[2], -0.3, epsilon = 1e-15); // depletion

        assert_abs_diff_eq!(result.max_accumulation, 0.4, epsilon = 1e-15);
        assert_abs_diff_eq!(result.max_depletion, -0.3, epsilon = 1e-15);

        // Integrated: (-0.3 + 0.4 - 0.3) * 1.0 = -0.2
        assert_abs_diff_eq!(result.integrated_delta_rho, -0.2, epsilon = 1e-15);
    }

    // =========================================================================
    // H2 golden value tests (PySCF 2.11.0 reference)
    // =========================================================================
    //
    // System: H2 at (0,0,0) and (0,0,1.4) in bohr, RHF/STO-3G
    // Atomic reference: H atom, UHF/STO-3G, spherically averaged
    //
    // Delta-rho = rho_mol - rho_pro where rho_pro = rho_H(r1) + rho_H(r2)
    //
    // PySCF 2.11.0 golden values (mol.unit='bohr'):
    //   H1 (0,0,0):     rho_mol=3.548920731e-01  rho_pro=4.140534050e-01  Delta=-5.916133190e-02
    //   midpoint:        rho_mol=2.559266713e-01  rho_pro=2.123318926e-01  Delta=+4.359477874e-02
    //   H2 (0,0,1.4):   rho_mol=3.548920731e-01  rho_pro=4.140534050e-01  Delta=-5.916133190e-02
    //   behind H1:       rho_mol=4.135653945e-02  rho_pro=5.126016313e-02  Delta=-9.903623680e-03
    //   behind H2:       rho_mol=4.135653945e-02  rho_pro=5.126016313e-02  Delta=-9.903623680e-03
    //   off-axis:        rho_mol=7.156822072e-02  rho_pro=5.937722581e-02  Delta=+1.219099491e-02

    #[test]
    fn test_h2_promolecule_golden_values() {
        // Verify the promolecule density at specific points matches PySCF
        let h_profile = get_atomic_density(1).unwrap();

        // At H1 (0,0,0): distance to H1=0, distance to H2=1.4
        let rho_pro_h1 = h_profile.evaluate(0.0) + h_profile.evaluate(1.4);
        assert_abs_diff_eq!(rho_pro_h1, 4.140534049811135e-01, epsilon = 1e-6);

        // At midpoint (0,0,0.7): distance to both = 0.7
        let rho_pro_mid = h_profile.evaluate(0.7) + h_profile.evaluate(0.7);
        assert_abs_diff_eq!(rho_pro_mid, 2.123318925607856e-01, epsilon = 1e-6);

        // Behind H1 (0,0,-1): distance to H1=1.0, distance to H2=2.4
        let rho_pro_behind1 = h_profile.evaluate(1.0) + h_profile.evaluate(2.4);
        assert_abs_diff_eq!(rho_pro_behind1, 5.126016313240526e-02, epsilon = 1e-6);

        // Off-axis (1,0,0.7): distances = sqrt(1+0.49) = sqrt(1.49)
        let r_off = (1.0f64 + 0.49).sqrt();
        let rho_pro_off = h_profile.evaluate(r_off) + h_profile.evaluate(r_off);
        assert_abs_diff_eq!(rho_pro_off, 5.937722581241780e-02, epsilon = 1e-6);
    }

    #[test]
    fn test_h2_difference_density_accumulation_at_midpoint() {
        // At the H2 bond midpoint, electrons accumulate: Delta-rho > 0
        // PySCF golden: Delta = +4.359477873800197e-02

        // We test the promolecule contribution only (molecular density
        // would require running SCF, which is outside this unit test scope)
        let h_profile = get_atomic_density(1).unwrap();
        let rho_pro_mid = h_profile.evaluate(0.7) + h_profile.evaluate(0.7);

        // The molecular density at midpoint is 2.559266712987876e-01
        // (from density.rs golden tests with the same H2 system)
        let rho_mol_mid = 2.559266712987876e-01;
        let delta_mid = rho_mol_mid - rho_pro_mid;

        assert!(
            delta_mid > 0.0,
            "Delta-rho at H2 midpoint should be positive (accumulation), got {delta_mid}"
        );
        assert_abs_diff_eq!(delta_mid, 4.359477873800197e-02, epsilon = 1e-5);
    }

    #[test]
    fn test_h2_difference_density_depletion_at_nuclei() {
        // At H atom positions, electrons deplete: Delta-rho < 0
        // PySCF golden: Delta = -5.916133190007877e-02

        let h_profile = get_atomic_density(1).unwrap();
        let rho_pro_h1 = h_profile.evaluate(0.0) + h_profile.evaluate(1.4);

        // Molecular density at H1 (0,0,0) from density.rs golden tests
        let rho_mol_h1 = 3.548920730810347e-01;
        let delta_h1 = rho_mol_h1 - rho_pro_h1;

        assert!(
            delta_h1 < 0.0,
            "Delta-rho at H1 should be negative (depletion), got {delta_h1}"
        );
        assert_abs_diff_eq!(delta_h1, -5.916133190007877e-02, epsilon = 1e-5);
    }

    #[test]
    fn test_h2_difference_density_symmetry() {
        // Delta-rho at (0,0,0) should equal Delta-rho at (0,0,1.4)
        let h_profile = get_atomic_density(1).unwrap();

        let rho_pro_h1 = h_profile.evaluate(0.0) + h_profile.evaluate(1.4);
        let rho_pro_h2 = h_profile.evaluate(1.4) + h_profile.evaluate(0.0);

        // Promolecule densities must be equal
        assert_abs_diff_eq!(rho_pro_h1, rho_pro_h2, epsilon = 1e-15);

        // With equal molecular densities (from density.rs golden tests):
        let rho_mol_h1 = 3.548920730810347e-01;
        let rho_mol_h2 = 3.548920730810345e-01;
        let delta_h1 = rho_mol_h1 - rho_pro_h1;
        let delta_h2 = rho_mol_h2 - rho_pro_h2;

        assert_abs_diff_eq!(delta_h1, delta_h2, epsilon = 1e-10);
    }

    #[test]
    fn test_h2_behind_nuclei_golden() {
        // Behind H1 (0,0,-1) and behind H2 (0,0,2.4)
        // PySCF golden: Delta = -9.903623679900213e-03
        let h_profile = get_atomic_density(1).unwrap();

        let rho_pro_behind1 = h_profile.evaluate(1.0) + h_profile.evaluate(2.4);
        let rho_mol_behind1 = 4.135653945250505e-02;
        let delta_behind1 = rho_mol_behind1 - rho_pro_behind1;

        assert!(delta_behind1 < 0.0, "Delta behind H1 should be negative");
        assert_abs_diff_eq!(delta_behind1, -9.903623679900e-03, epsilon = 1e-5);

        // Behind H2 should be symmetric
        let rho_pro_behind2 = h_profile.evaluate(2.4) + h_profile.evaluate(1.0);
        let rho_mol_behind2 = 4.135653945250502e-02;
        let delta_behind2 = rho_mol_behind2 - rho_pro_behind2;

        assert_abs_diff_eq!(delta_behind1, delta_behind2, epsilon = 1e-10);
    }

    #[test]
    fn test_h2_off_axis_golden() {
        // Off-axis at (1,0,0.7)
        // PySCF golden: Delta = +1.219099490966809e-02
        let h_profile = get_atomic_density(1).unwrap();

        let r_off = (1.0f64 + 0.49).sqrt(); // distance from each H to (1,0,0.7)
        let rho_pro_off = h_profile.evaluate(r_off) + h_profile.evaluate(r_off);
        let rho_mol_off = 7.156822072208589e-02;
        let delta_off = rho_mol_off - rho_pro_off;

        assert!(
            delta_off > 0.0,
            "Delta off-axis near midpoint should be positive"
        );
        assert_abs_diff_eq!(delta_off, 1.219099490966809e-02, epsilon = 1e-5);
    }

    // =========================================================================
    // compute_difference_density integration test
    // =========================================================================

    #[test]
    fn test_compute_difference_density_conservation() {
        // Construct mock molecular and promolecular density grids
        // that satisfy density conservation (same total integral).
        let n = 27; // 3x3x3
        let spacing: f64 = 1.0;

        // Both grids integrate to 2.0 electrons
        let total_integral = 2.0;
        let per_point = total_integral / (n as f64 * spacing.powi(3));

        let mol: Vec<f64> = (0..n)
            .map(|i| per_point + 0.01 * (i as f64 - 13.0))
            .collect();
        let pro: Vec<f64> = (0..n)
            .map(|i| per_point - 0.01 * (i as f64 - 13.0))
            .collect();

        let result = compute_difference_density(&mol, &pro, spacing).unwrap();

        // The sums differ, so integrated_delta_rho won't be exactly 0
        // But it should be small relative to total electrons
        let mol_sum: f64 = mol.iter().sum::<f64>() * spacing.powi(3);
        let pro_sum: f64 = pro.iter().sum::<f64>() * spacing.powi(3);
        let expected_delta = mol_sum - pro_sum;
        assert_abs_diff_eq!(result.integrated_delta_rho, expected_delta, epsilon = 1e-12);
    }
}
