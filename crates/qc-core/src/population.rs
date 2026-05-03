//! Mulliken and Lowdin population analysis
//!
//! Partitions the total electron density among atoms in a molecule,
//! yielding atomic charges that reveal charge transfer and polarity.
//!
//! # Methods
//!
//! ## Mulliken (1955)
//!
//! The simplest partitioning: assign the diagonal elements of the
//! density-overlap product DS to each atom's basis functions.
//!
//! ```text
//! N_A = sum_{mu in A} (DS)_{mu,mu}
//! q_A = Z_A - N_A
//! ```
//!
//! ## Lowdin (1950)
//!
//! Uses the symmetric orthogonalization S^{1/2} to distribute shared
//! overlap equally between atoms:
//!
//! ```text
//! D_L = S^{1/2} D S^{1/2}
//! N_A = sum_{mu in A} (D_L)_{mu,mu}
//! q_A = Z_A - N_A
//! ```
//!
//! # Conservation
//!
//! Both methods conserve total electron count:
//! ```text
//! sum_A N_A = Tr(DS) = n_electrons
//! sum_A q_A = sum_A Z_A - n_electrons = Q_molecular
//! ```
//!
//! # References
//!
//! - Mulliken, R.S. (1955). J. Chem. Phys. 23, 1833.
//! - Lowdin, P.-O. (1950). J. Chem. Phys. 18, 365.
//! - PySCF implementation: `pyscf/scf/hf.py` (mulliken_pop)

use crate::basis::atomic_number_to_symbol;
use nalgebra::{DMatrix, SymmetricEigen};
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Module version (matches crate version)
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// ============================================================================
// Error Types
// ============================================================================

/// Errors from population analysis
#[derive(Debug, Error, Clone, PartialEq)]
pub enum PopulationError {
    /// Matrix dimension mismatch
    #[error("Matrix dimension mismatch: expected {expected} elements, got {actual}")]
    DimensionMismatch {
        /// Expected number of elements (nbf * nbf)
        expected: usize,
        /// Actual number of elements
        actual: usize,
    },

    /// Basis function count mismatch
    #[error("Basis function count mismatch: atoms sum to {atom_sum}, expected {nbf}")]
    BasisCountMismatch {
        /// Sum of basis functions across atoms
        atom_sum: usize,
        /// Expected total
        nbf: usize,
    },

    /// Overlap matrix has non-positive eigenvalue
    #[error("Overlap matrix has non-positive eigenvalue: {0:.6e}")]
    NonPositiveOverlap(f64),
}

// ============================================================================
// Result Types
// ============================================================================

/// Population and charge for a single atom.
///
/// Contains both Mulliken and Lowdin results for side-by-side comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AtomCharge {
    /// Index of the atom in the molecule (0-based)
    pub atom_index: usize,

    /// Atomic number (Z)
    pub atomic_number: u8,

    /// Element symbol (e.g., "O", "H")
    pub element: String,

    /// Mulliken charge: q_A = Z_A - sum_{mu in A} (DS)_{mu,mu}
    pub mulliken_charge: f64,

    /// Lowdin charge: q_A = Z_A - sum_{mu in A} (S^{1/2} D S^{1/2})_{mu,mu}
    pub lowdin_charge: f64,

    /// Mulliken gross atomic population
    pub mulliken_population: f64,

    /// Lowdin atomic population
    pub lowdin_population: f64,
}

/// Combined Mulliken and Lowdin population analysis result.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PopulationResult {
    /// Per-atom charges and populations
    pub atoms: Vec<AtomCharge>,

    /// Total Mulliken charge (sum of all atomic charges)
    pub total_mulliken_charge: f64,

    /// Total Lowdin charge (sum of all atomic charges)
    pub total_lowdin_charge: f64,

    /// Computation time in microseconds
    pub compute_time_us: f64,
}

// ============================================================================
// Core Algorithm
// ============================================================================

/// Compute Mulliken and Lowdin population analysis.
///
/// Given the converged density matrix D, the overlap matrix S, and the
/// atom-to-basis-function mapping, this computes both Mulliken and Lowdin
/// atomic charges and populations.
///
/// # Arguments
///
/// * `density_matrix` - Row-major D (nbf x nbf). Must include the factor
///   of 2 for RHF (i.e., D = 2 * C_occ @ C_occ^T). Tr(DS) = n_electrons.
/// * `overlap_matrix` - Row-major S (nbf x nbf). Symmetric, positive definite.
/// * `nbf` - Number of basis functions.
/// * `atoms` - Slice of (atomic_number, n_basis_on_atom) tuples. The sum of
///   all n_basis_on_atom must equal nbf.
///
/// # Returns
///
/// `PopulationResult` with per-atom Mulliken and Lowdin charges and populations.
///
/// # Errors
///
/// - `PopulationError::DimensionMismatch` if matrix sizes are wrong
/// - `PopulationError::BasisCountMismatch` if atom basis counts don't sum to nbf
/// - `PopulationError::NonPositiveOverlap` if S has non-positive eigenvalues
///
/// # Algorithm
///
/// ## Mulliken
/// Computes only the diagonal of DS in O(N^2):
/// ```text
/// (DS)_{mu,mu} = sum_nu D_{mu,nu} * S_{nu,mu}
/// N_A = sum_{mu in A} (DS)_{mu,mu}
/// q_A = Z_A - N_A
/// ```
///
/// ## Lowdin
/// Eigendecomposes S, builds S^{1/2}, computes diagonal of S^{1/2} D S^{1/2}:
/// ```text
/// S = V * diag(sigma) * V^T
/// S^{1/2} = V * diag(sqrt(sigma)) * V^T
/// (D_L)_{mu,mu} = sum_{a,b} S^{1/2}_{mu,a} * D_{a,b} * S^{1/2}_{b,mu}
/// ```
///
/// # References
///
/// - Mulliken, R.S. (1955). J. Chem. Phys. 23, 1833.
/// - Lowdin, P.-O. (1950). J. Chem. Phys. 18, 365.
pub fn population_analysis(
    density_matrix: &[f64],
    overlap_matrix: &[f64],
    nbf: usize,
    atoms: &[(u8, usize)],
) -> Result<PopulationResult, PopulationError> {
    // std::time::Instant is not available on wasm32-unknown-unknown
    #[cfg(not(target_arch = "wasm32"))]
    let start = std::time::Instant::now();

    // --- Input validation ---
    let expected = nbf * nbf;
    if density_matrix.len() != expected {
        return Err(PopulationError::DimensionMismatch {
            expected,
            actual: density_matrix.len(),
        });
    }
    if overlap_matrix.len() != expected {
        return Err(PopulationError::DimensionMismatch {
            expected,
            actual: overlap_matrix.len(),
        });
    }
    let atom_sum: usize = atoms.iter().map(|(_, n)| n).sum();
    if atom_sum != nbf {
        return Err(PopulationError::BasisCountMismatch { atom_sum, nbf });
    }

    // --- Mulliken: compute diagonal of DS in O(N^2) ---
    // (DS)_{mu,mu} = sum_nu D_{mu,nu} * S_{nu,mu}
    // S is symmetric so S_{nu,mu} = S_{mu,nu}, but we use S_{nu,mu} directly
    let mut mulliken_per_ao = vec![0.0f64; nbf];
    for mu in 0..nbf {
        let mut ds_mu_mu = 0.0;
        for nu in 0..nbf {
            // D[mu,nu] is at index mu*nbf + nu (row-major)
            // S[nu,mu] is at index nu*nbf + mu (row-major)
            ds_mu_mu += density_matrix[mu * nbf + nu] * overlap_matrix[nu * nbf + mu];
        }
        mulliken_per_ao[mu] = ds_mu_mu;
    }

    // --- Lowdin: eigendecompose S, build S^{1/2}, compute diagonal of S^{1/2} D S^{1/2} ---
    let s_mat = DMatrix::from_row_slice(nbf, nbf, overlap_matrix);
    let eigen = SymmetricEigen::new(s_mat);

    // Verify all eigenvalues are positive
    let min_eigenvalue = eigen.eigenvalues.min();
    if min_eigenvalue <= 0.0 {
        return Err(PopulationError::NonPositiveOverlap(min_eigenvalue));
    }

    // Build S^{1/2} = V * diag(sqrt(lambda)) * V^T
    let mut s_sqrt = DMatrix::zeros(nbf, nbf);
    for i in 0..nbf {
        let sqrt_lambda = eigen.eigenvalues[i].sqrt();
        for mu in 0..nbf {
            let v_mu_i = eigen.eigenvectors[(mu, i)];
            if v_mu_i.abs() < 1e-16 {
                continue;
            }
            let scaled = v_mu_i * sqrt_lambda;
            for nu in 0..nbf {
                s_sqrt[(mu, nu)] += scaled * eigen.eigenvectors[(nu, i)];
            }
        }
    }

    // Compute diagonal of S^{1/2} D S^{1/2}
    // (SDS)_{mu,mu} = sum_{a,b} S_sqrt[mu,a] * D[a,b] * S_sqrt[b,mu]
    // Since S_sqrt is symmetric: S_sqrt[b,mu] = S_sqrt[mu,b]
    // So: (SDS)_{mu,mu} = sum_{a,b} S_sqrt[mu,a] * D[a,b] * S_sqrt[mu,b]
    //                    = sum_a S_sqrt[mu,a] * (sum_b D[a,b] * S_sqrt[mu,b])
    let mut lowdin_per_ao = vec![0.0f64; nbf];
    for mu in 0..nbf {
        let mut sds_mu_mu = 0.0;
        for a in 0..nbf {
            let s_mu_a: f64 = s_sqrt[(mu, a)];
            if s_mu_a.abs() < 1e-16 {
                continue;
            }
            // Compute sum_b D[a,b] * S_sqrt[mu,b]
            let mut ds_sum = 0.0;
            for b in 0..nbf {
                ds_sum += density_matrix[a * nbf + b] * s_sqrt[(mu, b)];
            }
            sds_mu_mu += s_mu_a * ds_sum;
        }
        lowdin_per_ao[mu] = sds_mu_mu;
    }

    // --- Aggregate per-atom populations and charges ---
    let mut atom_charges = Vec::with_capacity(atoms.len());
    let mut mu_offset = 0;

    for (atom_idx, &(z, n_basis)) in atoms.iter().enumerate() {
        let mut mulliken_pop = 0.0;
        let mut lowdin_pop = 0.0;

        for mu in mu_offset..mu_offset + n_basis {
            mulliken_pop += mulliken_per_ao[mu];
            lowdin_pop += lowdin_per_ao[mu];
        }

        let z_f64 = z as f64;
        let element = atomic_number_to_symbol(z).unwrap_or("?").to_string();

        atom_charges.push(AtomCharge {
            atom_index: atom_idx,
            atomic_number: z,
            element,
            mulliken_charge: z_f64 - mulliken_pop,
            lowdin_charge: z_f64 - lowdin_pop,
            mulliken_population: mulliken_pop,
            lowdin_population: lowdin_pop,
        });

        mu_offset += n_basis;
    }

    let total_mulliken: f64 = atom_charges.iter().map(|a| a.mulliken_charge).sum();
    let total_lowdin: f64 = atom_charges.iter().map(|a| a.lowdin_charge).sum();

    #[cfg(not(target_arch = "wasm32"))]
    let compute_time_us = start.elapsed().as_secs_f64() * 1_000_000.0;
    #[cfg(target_arch = "wasm32")]
    let compute_time_us = 0.0;

    Ok(PopulationResult {
        atoms: atom_charges,
        total_mulliken_charge: total_mulliken,
        total_lowdin_charge: total_lowdin,
        compute_time_us,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    /// Helper: build H2 STO-3G density and overlap matrices from the PresetSystem test data.
    fn h2_sto3g_data() -> (Vec<f64>, Vec<f64>, usize, Vec<(u8, usize)>) {
        // From crates/qc-core/src/scf/mod.rs PresetSystem::h2_sto3g_test()
        // Run SCF to get density matrix
        use crate::scf::{rhf_scf, ConvergenceProfile, PresetSystem, ScfConfig};

        let system = PresetSystem::h2_sto3g_test();
        let config = ScfConfig::new(ConvergenceProfile::Tight);
        let result = rhf_scf(&system, &config).expect("H2 SCF failed");

        let density = result.density_matrix.clone();
        let overlap = system.s_matrix.clone();
        let nbf = system.nbf;
        let atoms = vec![(1u8, 1usize), (1u8, 1usize)]; // Two H atoms, 1 basis each

        (density, overlap, nbf, atoms)
    }

    #[test]
    fn test_h2_symmetric_charges() {
        let (density, overlap, nbf, atoms) = h2_sto3g_data();
        let result = population_analysis(&density, &overlap, nbf, &atoms)
            .expect("Population analysis failed");

        assert_eq!(result.atoms.len(), 2);

        // By symmetry, both H atoms should have zero charge
        assert_abs_diff_eq!(result.atoms[0].mulliken_charge, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.atoms[1].mulliken_charge, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.atoms[0].lowdin_charge, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.atoms[1].lowdin_charge, 0.0, epsilon = 1e-10);

        // Total charge should be zero (neutral molecule)
        assert_abs_diff_eq!(result.total_mulliken_charge, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.total_lowdin_charge, 0.0, epsilon = 1e-10);

        // Populations should each be 1.0 (2 electrons, 2 atoms)
        assert_abs_diff_eq!(result.atoms[0].mulliken_population, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.atoms[1].mulliken_population, 1.0, epsilon = 1e-10);

        // Element symbols
        assert_eq!(result.atoms[0].element, "H");
        assert_eq!(result.atoms[1].element, "H");
    }

    #[test]
    fn test_h2o_mulliken_charges() {
        // H2O STO-3G with exact IQCP geometry
        // PySCF 2.11.0 reference:
        //   O: Mulliken charge = -0.3657447338
        //   H: Mulliken charge = +0.1828723669 (each)
        let (density, overlap, nbf, atoms) = h2o_sto3g_data();
        let result = population_analysis(&density, &overlap, nbf, &atoms)
            .expect("Population analysis failed");

        assert_eq!(result.atoms.len(), 3);

        // Mulliken charges
        assert_abs_diff_eq!(
            result.atoms[0].mulliken_charge,
            -0.3657447338,
            epsilon = 1e-4
        );
        assert_abs_diff_eq!(
            result.atoms[1].mulliken_charge,
            0.1828723669,
            epsilon = 1e-4
        );
        assert_abs_diff_eq!(
            result.atoms[2].mulliken_charge,
            0.1828723669,
            epsilon = 1e-4
        );

        // Charge conservation: sum = 0 for neutral molecule
        assert_abs_diff_eq!(result.total_mulliken_charge, 0.0, epsilon = 1e-10);

        // O should be negative, H should be positive
        assert!(
            result.atoms[0].mulliken_charge < 0.0,
            "O should have negative charge"
        );
        assert!(
            result.atoms[1].mulliken_charge > 0.0,
            "H should have positive charge"
        );
    }

    #[test]
    fn test_h2o_lowdin_charges() {
        // PySCF 2.11.0 reference:
        //   O: Lowdin charge = -0.2530126139
        //   H: Lowdin charge = +0.1265063069 (each)
        let (density, overlap, nbf, atoms) = h2o_sto3g_data();
        let result = population_analysis(&density, &overlap, nbf, &atoms)
            .expect("Population analysis failed");

        // Lowdin charges
        assert_abs_diff_eq!(result.atoms[0].lowdin_charge, -0.2530126139, epsilon = 1e-4);
        assert_abs_diff_eq!(result.atoms[1].lowdin_charge, 0.1265063069, epsilon = 1e-4);
        assert_abs_diff_eq!(result.atoms[2].lowdin_charge, 0.1265063069, epsilon = 1e-4);

        // Charge conservation
        assert_abs_diff_eq!(result.total_lowdin_charge, 0.0, epsilon = 1e-10);

        // Lowdin charges should be smaller in magnitude than Mulliken
        assert!(
            result.atoms[0].lowdin_charge.abs() < result.atoms[0].mulliken_charge.abs(),
            "Lowdin O charge magnitude should be smaller than Mulliken"
        );
    }

    #[test]
    fn test_h2o_populations() {
        let (density, overlap, nbf, atoms) = h2o_sto3g_data();
        let result = population_analysis(&density, &overlap, nbf, &atoms)
            .expect("Population analysis failed");

        // Mulliken populations should sum to n_electrons = 10
        let mull_total: f64 = result.atoms.iter().map(|a| a.mulliken_population).sum();
        assert_abs_diff_eq!(mull_total, 10.0, epsilon = 1e-10);

        // Lowdin populations should sum to n_electrons = 10
        let low_total: f64 = result.atoms.iter().map(|a| a.lowdin_population).sum();
        assert_abs_diff_eq!(low_total, 10.0, epsilon = 1e-10);
    }

    #[test]
    fn test_error_dimension_mismatch() {
        let density = vec![1.0; 4]; // 2x2
        let overlap = vec![1.0; 9]; // 3x3
        let atoms = vec![(1u8, 2usize)];

        let result = population_analysis(&density, &overlap, 2, &atoms);
        assert!(result.is_err());
        match result.unwrap_err() {
            PopulationError::DimensionMismatch { expected, actual } => {
                assert_eq!(expected, 4);
                assert_eq!(actual, 9);
            }
            other => panic!("Expected DimensionMismatch, got: {:?}", other),
        }
    }

    #[test]
    fn test_error_basis_count_mismatch() {
        let density = vec![1.0; 4]; // 2x2
        let overlap = vec![1.0, 0.0, 0.0, 1.0]; // identity
        let atoms = vec![(1u8, 1usize)]; // Only 1 basis function, but nbf=2

        let result = population_analysis(&density, &overlap, 2, &atoms);
        assert!(result.is_err());
        match result.unwrap_err() {
            PopulationError::BasisCountMismatch { atom_sum, nbf } => {
                assert_eq!(atom_sum, 1);
                assert_eq!(nbf, 2);
            }
            other => panic!("Expected BasisCountMismatch, got: {:?}", other),
        }
    }

    #[test]
    fn test_single_atom_identity_overlap() {
        // Trivial case: 1 atom with identity overlap
        // D = [[2.0]] (1 electron pair), S = [[1.0]]
        let density = vec![2.0];
        let overlap = vec![1.0];
        let atoms = vec![(1u8, 1usize)]; // H with 1 basis

        let result =
            population_analysis(&density, &overlap, 1, &atoms).expect("Population analysis failed");

        assert_eq!(result.atoms.len(), 1);
        // Population = DS[0,0] = 2.0, charge = Z - pop = 1 - 2 = -1
        assert_abs_diff_eq!(result.atoms[0].mulliken_population, 2.0, epsilon = 1e-15);
        assert_abs_diff_eq!(result.atoms[0].mulliken_charge, -1.0, epsilon = 1e-15);
        assert_abs_diff_eq!(result.atoms[0].lowdin_population, 2.0, epsilon = 1e-15);
        assert_abs_diff_eq!(result.atoms[0].lowdin_charge, -1.0, epsilon = 1e-15);
    }

    #[test]
    fn test_performance() {
        let (density, overlap, nbf, atoms) = h2o_sto3g_data();

        let start = std::time::Instant::now();
        for _ in 0..100 {
            let _ = population_analysis(&density, &overlap, nbf, &atoms);
        }
        let elapsed = start.elapsed();

        // 100 iterations should take well under 100ms total
        // Single call should be < 1ms for nbf=7
        assert!(
            elapsed.as_millis() < 1000,
            "100 population analyses took {:?}, expected < 1000ms",
            elapsed
        );
    }

    // ========================================================================
    // H2O test data helper
    // ========================================================================

    /// Generate H2O STO-3G test data using on-the-fly integrals + SCF.
    fn h2o_sto3g_data() -> (Vec<f64>, Vec<f64>, usize, Vec<(u8, usize)>) {
        use crate::basis::{Atom, BasisSet};
        use crate::integrals;
        use crate::scf::{rhf_scf, ConvergenceProfile, PresetSystem, ScfConfig};

        // Build basis set from IQCP's exact H2O geometry (bohr)
        let o = Atom::new(8, [0.0, 0.0, 0.2217282]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.4305447, -0.8869128]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.4305447, -0.8869128]).unwrap();

        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();
        let nbf = basis.n_basis;

        // Compute integrals
        let s_matrix = integrals::overlap_matrix(&basis);
        let t_matrix = integrals::kinetic_matrix(&basis);
        let v_matrix = integrals::nuclear_matrix(&basis);
        let h_core: Vec<f64> = t_matrix
            .iter()
            .zip(v_matrix.iter())
            .map(|(t, v)| t + v)
            .collect();
        let eri = integrals::eri_compressed(&basis);

        // Build PresetSystem for SCF
        let system = PresetSystem {
            system_id: "h2o_test".to_string(),
            label: "H2O test".to_string(),
            nbf,
            nelec: basis.n_electrons,
            e_nuc: basis.nuclear_repulsion,
            s_matrix: s_matrix.clone(),
            h_core,
            eri_compressed: eri,
        };

        let config = ScfConfig::new(ConvergenceProfile::Tight);
        let result = rhf_scf(&system, &config).expect("H2O SCF failed");

        // Build atom-to-basis mapping
        let atoms: Vec<(u8, usize)> = basis
            .atoms
            .iter()
            .enumerate()
            .map(|(idx, atom)| {
                let n_basis: usize = basis
                    .shells
                    .iter()
                    .filter(|s| s.atom_idx == idx)
                    .map(|s| s.n_basis_functions())
                    .sum();
                (atom.atomic_number, n_basis)
            })
            .collect();

        (result.density_matrix, s_matrix, nbf, atoms)
    }
}
