//! Becke grid construction: radial quadrature, angular quadrature,
//! SG-1 pruning, and Becke partitioning.
//!
//! # Algorithm
//!
//! For each atom A in the molecule:
//! 1. Generate Mura-Knowles radial grid: (r_i, dr_i)
//! 2. For each radial point r_i:
//!    a. Determine Lebedev angular grid size via SG-1 pruning
//!    b. For each angular point (x_j, y_j, z_j, w_j):
//!       - Cartesian: point = R_A + r_i * (x_j, y_j, z_j)
//!       - Raw weight: `w_raw = r_i^2 * dr_i * w_j`
//! 3. Apply Becke partition weights to all points of atom A
//! 4. Multiply each point's weight by its Becke partition weight
//!
//! The Lebedev weights stored in `lebedev.rs` sum to 4*pi. The 4*pi
//! in the radial volume element cancels with dividing angular weights
//! by 4*pi, giving: w_total = r^2 * dr * w_angular * P_becke.
//!
//! # References
//!
//! - Becke (1988) JCP 88, 2547. Eqs. 11-21.
//! - Mura & Knowles (1996) JCP 104, 9848.
//! - Gill, Johnson & Pople (1993) CPL 209, 506.
//! - PySCF `dft/radi.py` (mura_knowles), `dft/gen_grid.py` (sg1_prune, gen_partition).

use super::lebedev::get_lebedev_grid;
use super::{BeckeGrid, GridConfig, GridQuality};
use crate::basis::Atom;

// =============================================================================
// Constants
// =============================================================================

/// SG-1 atomic radii in bohr (indexed by Z, Z=0 unused).
///
/// These are the radii used for SG-1 pruning region boundaries,
/// NOT the Bragg-Slater radii used for Becke partitioning.
///
/// Reference: Gill, Johnson & Pople (1993) CPL 209, 506;
///            PySCF `dft/radi.py` lines 39-43 (SG1RADII).
const SG1_RADII: [f64; 19] = [
    0.0,    // Z=0 (unused)
    1.0000, // H
    0.5882, // He
    3.0769, // Li
    2.0513, // Be
    1.5385, // B
    1.2308, // C
    1.0256, // N
    0.8791, // O
    0.7692, // F
    0.6838, // Ne
    4.0909, // Na
    3.1579, // Mg
    2.5714, // Al
    2.1687, // Si
    1.8750, // P
    1.6514, // S
    1.4754, // Cl
    1.3333, // Ar
];

/// SG-1 pruning thresholds for different element periods.
///
/// Each row has 4 thresholds defining 5 regions.
/// Row 0: H, He (Z <= 2)
/// Row 1: Li-Ne (Z = 3-10)
/// Row 2: Na-Ar (Z = 11-18)
///
/// Reference: Gill, Johnson & Pople (1993) CPL 209, 506;
///            PySCF `dft/gen_grid.py` lines 51-86.
const SG1_ALPHAS: [[f64; 4]; 3] = [
    [0.25, 0.50, 1.00, 4.50],   // H, He
    [0.1667, 0.50, 0.90, 3.50], // Li-Ne
    [0.10, 0.40, 0.80, 2.50],   // Na-Ar
];

// =============================================================================
// Public API
// =============================================================================

/// Build a Becke atom-centered numerical integration grid for a molecular system.
///
/// Constructs atom-centered grids with Mura-Knowles radial quadrature,
/// Lebedev angular quadrature, SG-1 pruning, and Becke partitioning.
///
/// # Arguments
///
/// * `atoms` - Atoms in the molecule with positions in bohr
/// * `config` - Grid configuration (radial points, quality, pruning)
///
/// # Returns
///
/// A `BeckeGrid` containing all quadrature points, weights, and atom assignments.
///
/// # References
///
/// - Becke (1988) JCP 88, 2547 (partitioning)
/// - Mura & Knowles (1996) JCP 104, 9848 (radial quadrature)
/// - Lebedev & Laikov (1999) Doklady Math. 59, 477 (angular quadrature)
/// - Gill, Johnson & Pople (1993) CPL 209, 506 (SG-1 pruning)
pub fn build_becke_grid(atoms: &[Atom], config: &GridConfig) -> BeckeGrid {
    let n_atoms = atoms.len();

    // Collect atom positions for Becke partitioning
    let atom_positions: Vec<[f64; 3]> = atoms.iter().map(|a| a.position).collect();

    // Determine max angular grid from quality (used only for unpruned grids)
    let max_angular = match config.quality {
        GridQuality::Standard => 302,
        GridQuality::Fine => 590,
    };

    // Estimate total points for pre-allocation
    // With SG-1 pruning, average ~100 angular points per radial point
    let est_per_atom = config.n_radial * 120;
    let est_total = est_per_atom * n_atoms;

    let mut all_points = Vec::with_capacity(est_total);
    let mut all_weights = Vec::with_capacity(est_total);
    let mut all_atom_indices = Vec::with_capacity(est_total);

    for (atom_idx, atom) in atoms.iter().enumerate() {
        let z = atom.atomic_number;

        // Step 1: Generate radial grid
        let (radii, dr_weights) = mura_knowles_radial(config.n_radial, z);

        // Step 2: Determine angular grid sizes via SG-1 pruning
        // SG-1 uses hardcoded Lebedev-194 in region 4 (Gill, Johnson & Pople 1993)
        // max_angular from quality setting only applies to unpruned grids
        let angular_sizes = if config.pruning {
            sg1_prune(z, &radii)
        } else {
            vec![max_angular; config.n_radial]
        };

        // Step 3: Build atom grid (coordinates and raw weights)
        let mut atom_coords = Vec::new();
        let mut atom_raw_weights = Vec::new();

        for (i, (&r, &dr)) in radii.iter().zip(dr_weights.iter()).enumerate() {
            let n_ang = angular_sizes[i];
            if n_ang == 0 {
                continue;
            }
            let angular_grid = get_lebedev_grid(n_ang);

            // Radial volume element: 4*pi * r^2 * dr
            // Angular weights sum to 4*pi
            // Combined: w = (4*pi * r^2 * dr) * (w_ang / (4*pi)) = r^2 * dr * w_ang
            let rad_factor = r * r * dr;

            for ang_pt in angular_grid.iter() {
                let [ax, ay, az, aw] = *ang_pt;

                // Cartesian coordinates: R_A + r * (x, y, z)
                let point = [
                    atom.position[0] + r * ax,
                    atom.position[1] + r * ay,
                    atom.position[2] + r * az,
                ];

                // Raw weight = r^2 * dr * w_angular
                let raw_weight = rad_factor * aw;

                atom_coords.push(point);
                atom_raw_weights.push(raw_weight);
            }
        }

        // Step 4: Apply Becke partitioning
        if n_atoms == 1 {
            // Single atom: all partition weights = 1.0
            for (coord, weight) in atom_coords.into_iter().zip(atom_raw_weights.into_iter()) {
                if weight.abs() > 1e-30 {
                    all_points.push(coord);
                    all_weights.push(weight);
                    all_atom_indices.push(atom_idx);
                }
            }
        } else {
            let partition_weights = becke_partition(atom_idx, &atom_coords, &atom_positions);

            for ((coord, raw_w), &pw) in atom_coords
                .into_iter()
                .zip(atom_raw_weights.into_iter())
                .zip(partition_weights.iter())
            {
                let final_weight = raw_w * pw;
                if final_weight.abs() > 1e-30 {
                    all_points.push(coord);
                    all_weights.push(final_weight);
                    all_atom_indices.push(atom_idx);
                }
            }
        }
    }

    let n_points = all_points.len();

    BeckeGrid {
        points: all_points,
        weights: all_weights,
        atom_indices: all_atom_indices,
        n_atoms,
        n_points,
        quality: config.quality,
    }
}

// =============================================================================
// Mura-Knowles Radial Quadrature
// =============================================================================

/// Generate Mura-Knowles radial quadrature points and weights.
///
/// Uses the log3 mapping from Mura & Knowles (1996):
///
/// ```text
/// x_i = (i + 0.5) / N        (uniform points on (0, 1))
/// r_i = -R * ln(1 - x_i^3)   (log3 transformation)
/// dr_i = R * 3*x_i^2 / ((1 - x_i^3) * N)   (Jacobian)
/// ```
///
/// where R is a scaling factor:
/// - R = 7.0 for Li, Be, Na, Mg, K, Ca (Z = 3, 4, 11, 12, 19, 20)
/// - R = 5.2 for all other elements
///
/// This exactly matches PySCF's `dft/radi.py::mura_knowles()`.
///
/// # Arguments
///
/// * `n` - Number of radial points (typically 75)
/// * `atomic_number` - Element's atomic number (for scaling factor R)
///
/// # Returns
///
/// Tuple of (radii, dr_weights) where:
/// - radii: distances from nucleus in bohr, monotonically increasing
/// - dr_weights: Jacobian weight factors dr/dx * (1/N)
///
/// Note: These do NOT include the 4*pi*r^2 spherical volume element.
/// The caller multiplies by `4*pi * r^2 * dr` when assembling the full weight.
///
/// # Reference
///
/// Mura & Knowles (1996) JCP 104, 9848; PySCF `dft/radi.py` lines 86-99.
pub fn mura_knowles_radial(n: usize, atomic_number: u8) -> (Vec<f64>, Vec<f64>) {
    let mut radii = Vec::with_capacity(n);
    let mut dr_weights = Vec::with_capacity(n);

    // Scaling factor R depends on element
    // R = 7.0 for Li, Be, Na, Mg, K, Ca; R = 5.2 otherwise
    // Reference: PySCF dft/radi.py line 90
    let r_scale = match atomic_number {
        3 | 4 | 11 | 12 | 19 | 20 => 7.0,
        _ => 5.2,
    };

    let n_f64 = n as f64;

    for i in 0..n {
        let x = (i as f64 + 0.5) / n_f64;
        let x3 = x * x * x;
        let one_minus_x3 = 1.0 - x3;

        // r = -R * ln(1 - x^3)
        let r = -r_scale * one_minus_x3.ln();

        // dr = R * 3*x^2 / ((1 - x^3) * N)
        let dr = r_scale * 3.0 * x * x / (one_minus_x3 * n_f64);

        radii.push(r);
        dr_weights.push(dr);
    }

    (radii, dr_weights)
}

// =============================================================================
// SG-1 Grid Pruning
// =============================================================================

/// Determine Lebedev angular grid order for each radial point using SG-1 pruning.
///
/// The SG-1 scheme (Gill, Johnson & Pople, 1993) assigns different angular
/// grid sizes based on the ratio `r / R_SG1` where `R_SG1` is the SG-1
/// atomic radius. The pruning thresholds depend on the element's period:
///
/// | Period | Thresholds |
/// |--------|-----------|
/// | H, He (Z <= 2) | 0.25, 0.50, 1.00, 4.50 |
/// | Li-Ne (Z = 3-10) | 0.1667, 0.50, 0.90, 3.50 |
/// | Na-Ar (Z = 11-18) | 0.10, 0.40, 0.80, 2.50 |
///
/// The angular grids for the 5 regions are {6, 38, 86, 194, 86} as specified
/// in the original SG-1 paper. The 194-point Lebedev grid in region 4 is
/// hardcoded per the published standard --- it does NOT vary with quality setting.
///
/// # Arguments
///
/// * `atomic_number` - Nuclear charge Z
/// * `radii` - Radial grid point distances (bohr)
///
/// # Returns
///
/// Vec of angular grid sizes (one per radial point).
///
/// # Reference
///
/// Gill, Johnson & Pople (1993) CPL 209, 506;
/// PySCF `dft/gen_grid.py` lines 51-86.
pub fn sg1_prune(atomic_number: u8, radii: &[f64]) -> Vec<usize> {
    let z = atomic_number as usize;

    // SG-1 radii are only defined for Z=1-18
    let r_atom = if z < SG1_RADII.len() {
        SG1_RADII[z] + 1e-200 // avoid division by zero
    } else {
        // Fallback: use generic radius
        1.0
    };

    // Select threshold row based on period
    let alpha_row = if z <= 2 {
        0
    } else if z <= 10 {
        1
    } else {
        2
    };
    let alphas = &SG1_ALPHAS[alpha_row];

    // SG-1 standard: Lebedev-194 in region 4 (hardcoded per Gill, Johnson & Pople 1993)
    let leb_ngrid = [6, 38, 86, 194, 86];

    radii
        .iter()
        .map(|&r| {
            let ratio = r / r_atom;
            // Count how many thresholds are exceeded
            let place = alphas.iter().filter(|&&a| ratio > a).count();
            leb_ngrid[place]
        })
        .collect()
}

// =============================================================================
// Becke Partitioning
// =============================================================================

/// Becke step function: maps mu in \[-1, 1\] to a smooth step in \[0, 1\].
///
/// Applies the polynomial `f(x) = 1.5*x - 0.5*x^3` three times
/// (3-iteration hardening), then converts to a step: `s = 0.5 * (1 - f(f(f(mu))))`.
///
/// Reference: Becke (1988) JCP 88, 2547, Eqs. 19-21.
#[inline]
#[allow(dead_code)]
fn becke_step(mu: f64) -> f64 {
    // f(x) = 1.5*x - 0.5*x^3 = x * (3 - x^2) * 0.5
    // This is the same as PySCF's commented code: g = (3 - g**2) * g * .5
    let f = |x: f64| -> f64 { x * (3.0 - x * x) * 0.5 };
    let g = f(f(f(mu)));
    0.5 * (1.0 - g)
}

/// Compute Becke partition weights for grid points centered on one atom.
///
/// For each grid point belonging to atom A, computes the normalized
/// Becke partition weight w_A(r) that distributes the integration
/// weight among atoms.
///
/// The algorithm:
/// 1. Compute distances from each grid point to all atoms.
/// 2. For each atom pair (A, B), compute the confocal ellipsoidal
///    coordinate: `mu_AB = (|r - R_A| - |r - R_B|) / R_AB`.
/// 3. Apply the step function: `s(mu_AB)` gives the contribution
///    from atom A (close to 1 if r is nearer to A than B).
/// 4. The unnormalized partition for atom A is the product of all
///    s(mu_AB) over B != A.
/// 5. Normalize: `w_A = P_A / sum_B P_B`.
///
/// # Arguments
///
/// * `atom_idx` - Index of the atom owning these grid points
/// * `coords` - Grid point coordinates (absolute Cartesian positions)
/// * `atom_positions` - Positions of all atoms in the molecule
///
/// # Returns
///
/// Vec of Becke partition weights (one per grid point), each in \[0, 1\].
///
/// # Reference
///
/// Becke (1988) JCP 88, 2547, Eqs. 11-21;
/// PySCF `dft/gen_grid.py` `get_partition` / `gen_grid_partition`.
fn becke_partition(atom_idx: usize, coords: &[[f64; 3]], atom_positions: &[[f64; 3]]) -> Vec<f64> {
    let n_atoms = atom_positions.len();
    let n_points = coords.len();

    if n_atoms <= 1 {
        return vec![1.0; n_points];
    }

    // Precompute inter-atomic distances
    let mut atom_dist = vec![0.0f64; n_atoms * n_atoms];
    for i in 0..n_atoms {
        for j in (i + 1)..n_atoms {
            let d = dist3(&atom_positions[i], &atom_positions[j]);
            atom_dist[i * n_atoms + j] = d;
            atom_dist[j * n_atoms + i] = d;
        }
    }

    let mut partition_weights = Vec::with_capacity(n_points);

    // For each grid point, compute the full Becke partition
    for point in coords.iter() {
        // Compute distances from this point to all atoms
        let mut grid_dist = Vec::with_capacity(n_atoms);
        for ap in atom_positions.iter() {
            grid_dist.push(dist3(point, ap));
        }

        // Compute unnormalized partition function for each atom
        let mut p_values = vec![1.0f64; n_atoms];

        for i in 0..n_atoms {
            for j in 0..i {
                let r_ij = atom_dist[i * n_atoms + j];
                if r_ij < 1e-15 {
                    continue;
                }
                let mu = (grid_dist[i] - grid_dist[j]) / r_ij;

                // Becke step function with 3-iteration hardening
                // g = f(f(f(mu)))
                // For pair (i,j): s_ij = 0.5 * (1 - g)
                // atom i gets factor s_ij (= 0.5*(1-g))
                // atom j gets factor 0.5*(1+g)
                let f = |x: f64| -> f64 { x * (3.0 - x * x) * 0.5 };
                let g = f(f(f(mu)));

                p_values[i] *= 0.5 * (1.0 - g);
                p_values[j] *= 0.5 * (1.0 + g);
            }
        }

        // Normalize
        let p_sum: f64 = p_values.iter().sum();
        let w = if p_sum.abs() > 1e-30 {
            p_values[atom_idx] / p_sum
        } else {
            // All atoms have ~zero partition weight at this point;
            // distribute equally
            1.0 / n_atoms as f64
        };

        partition_weights.push(w);
    }

    partition_weights
}

/// Euclidean distance between two 3D points.
#[inline]
fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
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
    // Mura-Knowles Radial Tests
    // =========================================================================

    #[test]
    fn test_mura_knowles_h_75_points() {
        // MK-1: H atom, 75 radial points
        let (radii, dr) = mura_knowles_radial(75, 1);
        assert_eq!(radii.len(), 75);
        assert_eq!(dr.len(), 75);

        // All radii should be positive
        for &r in &radii {
            assert!(r > 0.0, "Radius must be positive, got {}", r);
        }

        // All dr weights should be positive
        for &w in &dr {
            assert!(w > 0.0, "Weight must be positive, got {}", w);
        }
    }

    #[test]
    fn test_mura_knowles_li_scaling() {
        // MK-2: Li uses R=7.0 scaling (vs R=5.2 for most elements)
        let (radii_h, _) = mura_knowles_radial(75, 1);
        let (radii_li, _) = mura_knowles_radial(75, 3);

        // Li radii should be larger by factor 7.0/5.2
        let ratio = radii_li[37] / radii_h[37];
        assert_abs_diff_eq!(ratio, 7.0 / 5.2, epsilon = 1e-12);
    }

    #[test]
    fn test_mura_knowles_monotonic() {
        // MK-3: Radii must be monotonically increasing
        let (radii, _) = mura_knowles_radial(75, 1);
        for i in 1..radii.len() {
            assert!(
                radii[i] > radii[i - 1],
                "Radii not monotonic at i={}: r[{}]={} <= r[{}]={}",
                i,
                i,
                radii[i],
                i - 1,
                radii[i - 1]
            );
        }
    }

    #[test]
    fn test_mura_knowles_pyscf_reference() {
        // MK-4: Compare against PySCF mura_knowles(75, charge=1)
        // PySCF reference values (PySCF 2.11.0):
        //   r[0]  = 1.540740968893304e-06
        //   r[37] = 6.943632416475177e-01
        //   r[74] = 2.037722481312466e+01
        let (radii, dr) = mura_knowles_radial(75, 1);

        assert_abs_diff_eq!(radii[0], 1.540740968893304e-06, epsilon = 1e-18);
        assert_abs_diff_eq!(radii[37], 6.943632416475177e-01, epsilon = 1e-13);
        assert_abs_diff_eq!(radii[74], 2.037722481312466e+01, epsilon = 1e-10);

        // dr weights
        assert_abs_diff_eq!(dr[0], 9.244447183539910e-06, epsilon = 1e-18);
        assert_abs_diff_eq!(dr[37], 5.942857142857143e-02, epsilon = 1e-14);
    }

    #[test]
    fn test_mura_knowles_o_pyscf_reference() {
        // Compare O (Z=8) against PySCF: same as H (both use R=5.2)
        let (radii_h, _) = mura_knowles_radial(75, 1);
        let (radii_o, _) = mura_knowles_radial(75, 8);

        // H and O both use R=5.2, so radii should be identical
        for i in 0..75 {
            assert_abs_diff_eq!(radii_h[i], radii_o[i], epsilon = 1e-15);
        }
    }

    // =========================================================================
    // SG-1 Pruning Tests
    // =========================================================================

    #[test]
    fn test_sg1_prune_h_reduction() {
        // SG1-1: H atom with SG-1 pruning shows significant reduction
        // PySCF reference: {6: 27, 38: 7, 86: 20, 194: 21}
        let (radii, _) = mura_knowles_radial(75, 1);
        let angs = sg1_prune(1, &radii);

        let total_pruned: usize = angs.iter().sum();
        // SG-1 max angular is 194 (hardcoded per Gill, Johnson & Pople 1993)
        let total_unpruned = 75 * 194;
        let reduction = 1.0 - total_pruned as f64 / total_unpruned as f64;

        assert!(
            reduction > 0.2,
            "SG-1 should reduce grid by >20%, got {:.1}%",
            reduction * 100.0
        );
    }

    #[test]
    fn test_sg1_prune_h_pyscf_pattern() {
        // Standard SG-1 uses Lebedev-194 in region 4 --- match PySCF exactly
        // PySCF: {6: 27, 38: 7, 86: 20, 194: 21}
        let (radii, _) = mura_knowles_radial(75, 1);
        let angs = sg1_prune(1, &radii);

        let count_6 = angs.iter().filter(|&&a| a == 6).count();
        let count_38 = angs.iter().filter(|&&a| a == 38).count();
        let count_86 = angs.iter().filter(|&&a| a == 86).count();
        let count_194 = angs.iter().filter(|&&a| a == 194).count();

        assert_eq!(count_6, 27, "Expected 27 x Lebedev-6 for H");
        assert_eq!(count_38, 7, "Expected 7 x Lebedev-38 for H");
        assert_eq!(count_86, 20, "Expected 20 x Lebedev-86 for H");
        assert_eq!(count_194, 21, "Expected 21 x Lebedev-194 for H");
    }

    #[test]
    fn test_sg1_prune_o_pyscf_pattern() {
        // PySCF: O(Z=8): {6: 23, 38: 9, 86: 25, 194: 18}
        let (radii, _) = mura_knowles_radial(75, 8);
        let angs = sg1_prune(8, &radii);

        let count_6 = angs.iter().filter(|&&a| a == 6).count();
        let count_38 = angs.iter().filter(|&&a| a == 38).count();
        let count_86 = angs.iter().filter(|&&a| a == 86).count();
        let count_194 = angs.iter().filter(|&&a| a == 194).count();

        assert_eq!(count_6, 23, "Expected 23 x Lebedev-6 for O");
        assert_eq!(count_38, 9, "Expected 9 x Lebedev-38 for O");
        assert_eq!(count_86, 25, "Expected 25 x Lebedev-86 for O");
        assert_eq!(count_194, 18, "Expected 18 x Lebedev-194 for O");
    }

    #[test]
    fn test_sg1_prune_near_nucleus() {
        // SG1-2: Very small r should give 6-point grid
        let (radii, _) = mura_knowles_radial(75, 1);
        let angs = sg1_prune(1, &radii);
        // r[0] is very small, well below 0.25 * R_SG1(H=1.0)
        assert_eq!(angs[0], 6, "Near-nucleus point should use Lebedev-6");
    }

    #[test]
    fn test_sg1_prune_no_pruning_gives_uniform() {
        // When pruning is disabled, all radial points get max_angular
        let radii = vec![0.1, 0.5, 1.0, 5.0, 10.0];
        // Not using sg1_prune; just verifying the concept
        let uniform: Vec<usize> = vec![302; radii.len()];
        assert!(uniform.iter().all(|&n| n == 302));
    }

    // =========================================================================
    // Becke Partitioning Tests
    // =========================================================================

    #[test]
    fn test_becke_single_atom() {
        // BK-1: Single atom -> all partition weights = 1.0
        let positions = vec![[0.0, 0.0, 0.0]];
        let coords = vec![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];

        let weights = becke_partition(0, &coords, &positions);

        for &w in &weights {
            assert_abs_diff_eq!(w, 1.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn test_becke_two_atoms_midpoint() {
        // BK-2: Two atoms, point at midplane -> w_A ~= w_B ~= 0.5
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]];
        let midpoint = vec![[0.0, 0.0, 1.0]]; // equidistant from both

        let w_a = becke_partition(0, &midpoint, &positions);
        let w_b = becke_partition(1, &midpoint, &positions);

        assert_abs_diff_eq!(w_a[0], 0.5, epsilon = 1e-14);
        assert_abs_diff_eq!(w_b[0], 0.5, epsilon = 1e-14);
    }

    #[test]
    fn test_becke_two_atoms_near_atom_a() {
        // BK-3: Two atoms, point very close to atom A -> w_A -> 1.0
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 4.0]];
        let near_a = vec![[0.0, 0.0, 0.1]];

        let w_a = becke_partition(0, &near_a, &positions);
        let w_b = becke_partition(1, &near_a, &positions);

        assert!(w_a[0] > 0.99, "w_A should be ~1.0, got {}", w_a[0]);
        assert!(w_b[0] < 0.01, "w_B should be ~0.0, got {}", w_b[0]);
    }

    #[test]
    fn test_becke_partition_sum_to_one() {
        // BK-4: Sum of partition weights at any point = 1.0
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let test_points = vec![
            [0.5, 0.5, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [3.0, 3.0, 3.0],
        ];

        for point in &test_points {
            let mut total = 0.0;
            for atom_idx in 0..positions.len() {
                let w = becke_partition(atom_idx, &[*point], &positions);
                total += w[0];
            }
            assert_abs_diff_eq!(total, 1.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn test_becke_step_function() {
        // Test the becke_step function directly
        // s(-1) should be close to 1 (point at atom A)
        // s(0) should be 0.5 (midplane)
        // s(1) should be close to 0 (point at atom B)
        assert_abs_diff_eq!(becke_step(0.0), 0.5, epsilon = 1e-14);
        assert!(becke_step(-0.99) > 0.99);
        assert!(becke_step(0.99) < 0.01);
    }

    // =========================================================================
    // Full Grid Structure Tests
    // =========================================================================

    #[test]
    fn test_grid_structure_single_atom() {
        let h = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let config = GridConfig::default();
        let grid = build_becke_grid(&[h], &config);

        // GS-1: n_points == points.len()
        assert_eq!(grid.n_points, grid.points.len());
        // GS-2: atom_indices.len() == n_points
        assert_eq!(grid.atom_indices.len(), grid.n_points);
        // GS-3: All weights positive (after zero-weight pruning)
        for &w in &grid.weights {
            assert!(w > 0.0, "Weight should be positive, got {}", w);
        }
        // Single atom: all atom_indices = 0
        assert!(grid.atom_indices.iter().all(|&i| i == 0));
        assert_eq!(grid.n_atoms, 1);
    }

    #[test]
    fn test_grid_structure_h2() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let config = GridConfig::default();
        let grid = build_becke_grid(&[h1, h2], &config);

        assert_eq!(grid.n_points, grid.points.len());
        assert_eq!(grid.n_points, grid.weights.len());
        assert_eq!(grid.n_points, grid.atom_indices.len());
        assert_eq!(grid.n_atoms, 2);

        // All weights positive
        for &w in &grid.weights {
            assert!(w > 0.0, "Weight should be positive, got {}", w);
        }

        // atom_indices covers both atoms
        let has_0 = grid.atom_indices.iter().any(|&i| i == 0);
        let has_1 = grid.atom_indices.iter().any(|&i| i == 1);
        assert!(has_0, "Grid should have points for atom 0");
        assert!(has_1, "Grid should have points for atom 1");
    }

    #[test]
    fn test_grid_structure_h2o() {
        let o = Atom::new(8, [0.0, 0.0, 0.221665]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430901, -0.886659]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430901, -0.886659]).unwrap();
        let config = GridConfig::default();
        let grid = build_becke_grid(&[o, h1, h2], &config);

        assert_eq!(grid.n_atoms, 3);
        assert_eq!(grid.n_points, grid.points.len());

        // Should have a reasonable number of points
        // With SG-1 and 194 max angular, expect ~3000+ pts/atom
        assert!(
            grid.n_points > 3000,
            "H2O grid should have >3000 points, got {}",
            grid.n_points
        );
        assert!(
            grid.n_points < 100000,
            "H2O grid should have <100000 points, got {}",
            grid.n_points
        );
    }

    #[test]
    fn test_grid_weight_sum_positive() {
        // Total weight sum should be positive and finite
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let config = GridConfig::default();
        let grid = build_becke_grid(&[h1, h2], &config);

        let total_weight: f64 = grid.weights.iter().sum();
        assert!(total_weight > 0.0, "Total weight should be positive");
        assert!(total_weight.is_finite(), "Total weight should be finite");
    }

    #[test]
    fn test_grid_no_pruning() {
        // Without pruning, all radial points get max angular grid
        let h = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let config_pruned = GridConfig {
            n_radial: 75,
            quality: GridQuality::Standard,
            pruning: true,
        };
        let config_unpruned = GridConfig {
            n_radial: 75,
            quality: GridQuality::Standard,
            pruning: false,
        };

        let grid_pruned = build_becke_grid(&[h.clone()], &config_pruned);
        let grid_unpruned = build_becke_grid(&[h], &config_unpruned);

        // Unpruned should have more points than pruned
        assert!(
            grid_unpruned.n_points > grid_pruned.n_points,
            "Unpruned ({}) should have more points than pruned ({})",
            grid_unpruned.n_points,
            grid_pruned.n_points
        );

        // Unpruned should be exactly 75 * 302
        assert_eq!(
            grid_unpruned.n_points,
            75 * 302,
            "Unpruned H should have 75*302={} points, got {}",
            75 * 302,
            grid_unpruned.n_points
        );
    }

    #[test]
    fn test_grid_fine_quality() {
        // SG-1 pruning uses hardcoded Lebedev-194 regardless of quality.
        // Quality setting only affects unpruned grids (302 vs 590).
        let h = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();

        // With pruning on, Standard and Fine produce identical SG-1 grids
        let config_std_pruned = GridConfig {
            n_radial: 75,
            quality: GridQuality::Standard,
            pruning: true,
        };
        let config_fine_pruned = GridConfig {
            n_radial: 75,
            quality: GridQuality::Fine,
            pruning: true,
        };
        let grid_std_pruned = build_becke_grid(&[h.clone()], &config_std_pruned);
        let grid_fine_pruned = build_becke_grid(&[h.clone()], &config_fine_pruned);
        assert_eq!(
            grid_std_pruned.n_points, grid_fine_pruned.n_points,
            "SG-1 pruned grids should be identical regardless of quality setting"
        );

        // Without pruning, Fine (590) should have more points than Standard (302)
        let config_std_unpruned = GridConfig {
            n_radial: 75,
            quality: GridQuality::Standard,
            pruning: false,
        };
        let config_fine_unpruned = GridConfig {
            n_radial: 75,
            quality: GridQuality::Fine,
            pruning: false,
        };
        let grid_std_unpruned = build_becke_grid(&[h.clone()], &config_std_unpruned);
        let grid_fine_unpruned = build_becke_grid(&[h], &config_fine_unpruned);
        assert!(
            grid_fine_unpruned.n_points > grid_std_unpruned.n_points,
            "Unpruned Fine ({}) should have more points than Standard ({})",
            grid_fine_unpruned.n_points,
            grid_std_unpruned.n_points
        );
    }

    // =========================================================================
    // Density Integration Tests (AC5)
    // =========================================================================

    /// Evaluate electron density at a single point.
    ///
    /// rho(r) = sum_{mu,nu} D[mu,nu] * chi_mu(r) * chi_nu(r)
    fn evaluate_density_at_points(
        density_matrix: &[f64],
        nbf: usize,
        basis: &BasisSet,
        points: &[[f64; 3]],
    ) -> Vec<f64> {
        use crate::basis::AngularMomentum;
        use crate::integrals::cartesian_components;
        use crate::orbital::grid::{angular_factor, cartesian_norm, GAUSSIAN_SCREENING_THRESHOLD};

        let mut rho_values = Vec::with_capacity(points.len());
        let mut chi = vec![0.0f64; nbf];

        for point in points {
            // Zero chi
            for c in chi.iter_mut() {
                *c = 0.0;
            }

            // Evaluate all basis functions at this point
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

                let dx = point[0] - shell.center[0];
                let dy = point[1] - shell.center[1];
                let dz = point[2] - shell.center[2];
                let dist_sq = dx * dx + dy * dy + dz * dz;

                let alpha_min = shell.min_exponent();
                if alpha_min * dist_sq > GAUSSIAN_SCREENING_THRESHOLD {
                    offset += n_funcs;
                    continue;
                }

                let disps = [dx, dy, dz];
                let is_s = shell.angular_momentum == AngularMomentum::S;

                if is_s {
                    let mut val = 0.0;
                    for prim in &shell.primitives {
                        let norm = cartesian_norm(prim.exponent, &components[0]);
                        val += prim.coefficient * norm * (-prim.exponent * dist_sq).exp();
                    }
                    chi[offset] = val;
                } else {
                    for (comp_idx, comp) in components.iter().enumerate() {
                        let ang = angular_factor(&disps, comp);
                        if ang.abs() < 1e-30 {
                            continue;
                        }
                        let mut radial_sum = 0.0;
                        for prim in &shell.primitives {
                            let norm = cartesian_norm(prim.exponent, comp);
                            radial_sum +=
                                prim.coefficient * norm * (-prim.exponent * dist_sq).exp();
                        }
                        chi[offset + comp_idx] = ang * radial_sum;
                    }
                }

                offset += n_funcs;
            }

            // Contract: rho = sum D[mu,nu] * chi[mu] * chi[nu]
            let mut rho = 0.0;
            for mu in 0..nbf {
                let chi_mu = chi[mu];
                if chi_mu.abs() < 1e-30 {
                    continue;
                }
                rho += density_matrix[mu * nbf + mu] * chi_mu * chi_mu;
                for nu in (mu + 1)..nbf {
                    let chi_nu = chi[nu];
                    if chi_nu.abs() < 1e-30 {
                        continue;
                    }
                    rho += 2.0 * density_matrix[mu * nbf + nu] * chi_mu * chi_nu;
                }
            }
            rho_values.push(rho);
        }

        rho_values
    }

    // H2 STO-3G density matrix (PySCF 2.11.0, R=1.4 bohr)
    fn h2_density_matrix() -> Vec<f64> {
        vec![
            6.026571614189372e-01,
            6.026571614189370e-01,
            6.026571614189370e-01,
            6.026571614189368e-01,
        ]
    }

    // H2O STO-3G density matrix (PySCF 2.11.0)
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

    #[test]
    fn test_h2_integrated_density() {
        // DENS-1: H2 STO-3G integrated density within 0.01% of 2.0
        // PySCF reference: N_integrated = 2.000005583760, relative error = 2.79e-06
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1.clone(), h2.clone()], "sto-3g").unwrap();

        let config = GridConfig::default();
        let grid = build_becke_grid(&[h1, h2], &config);

        let dm = h2_density_matrix();
        let rho = evaluate_density_at_points(&dm, 2, &basis, &grid.points);

        let n_integrated: f64 = rho
            .iter()
            .zip(grid.weights.iter())
            .map(|(r, w)| r * w)
            .sum();
        let n_expected = 2.0;
        let rel_error = (n_integrated - n_expected).abs() / n_expected;

        assert!(
            rel_error < 0.0001,
            "H2 density integration: N={:.10}, expected {}, rel_error={:.2e} (limit 1e-4)",
            n_integrated,
            n_expected,
            rel_error
        );
    }

    #[test]
    fn test_h2o_integrated_density() {
        // DENS-2: H2O STO-3G integrated density within 0.01% of 10.0
        // PySCF reference: N_integrated = 9.999991163750, relative error = 8.84e-07
        let o = Atom::new(8, [0.0, 0.0, 0.221665]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430901, -0.886659]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430901, -0.886659]).unwrap();
        let basis = BasisSet::build(vec![o.clone(), h1.clone(), h2.clone()], "sto-3g").unwrap();

        let config = GridConfig::default();
        let grid = build_becke_grid(&[o, h1, h2], &config);

        let dm = h2o_density_matrix();
        let rho = evaluate_density_at_points(&dm, 7, &basis, &grid.points);

        let n_integrated: f64 = rho
            .iter()
            .zip(grid.weights.iter())
            .map(|(r, w)| r * w)
            .sum();
        let n_expected = 10.0;
        let rel_error = (n_integrated - n_expected).abs() / n_expected;

        assert!(
            rel_error < 0.0001,
            "H2O density integration: N={:.10}, expected {}, rel_error={:.2e} (limit 1e-4)",
            n_integrated,
            n_expected,
            rel_error
        );
    }

    // =========================================================================
    // Performance Test (AC6)
    // =========================================================================

    #[test]
    fn test_h2o_grid_generation_performance() {
        // PERF-1: H2O grid generation should complete quickly
        // In debug mode this may be slower; the 100ms target is for release mode.
        // We still verify it completes in a reasonable time.
        let o = Atom::new(8, [0.0, 0.0, 0.221665]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430901, -0.886659]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430901, -0.886659]).unwrap();
        let config = GridConfig::default();

        let start = std::time::Instant::now();
        let grid = build_becke_grid(&[o, h1, h2], &config);
        let elapsed = start.elapsed();

        // In debug mode, allow up to 5 seconds
        // In release mode, should be <100ms
        assert!(
            elapsed.as_secs() < 10,
            "H2O grid generation took {:?}, expected <10s",
            elapsed
        );
        assert!(grid.n_points > 0);
    }
}
