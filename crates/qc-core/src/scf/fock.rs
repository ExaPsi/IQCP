//! Fock matrix decomposition for educational inspection
//!
//! Decomposes the Fock matrix into its physical components:
//!   F = H^core + G, where G = J - 0.5*K
//!
//! The Coulomb (J) and exchange (K) matrices are computed separately so that
//! students can inspect each contribution independently. This function is NOT
//! used in the SCF hot path -- it is called once after convergence for
//! pedagogical inspection.
//!
//! # Theory
//!
//! The Fock matrix is the effective one-electron operator in Hartree-Fock theory.
//! It encodes both the one-electron core Hamiltonian (kinetic + nuclear attraction)
//! and the mean-field electron-electron interaction:
//!
//! ```text
//! F_{mn} = H^core_{mn} + G_{mn}
//! G_{mn} = J_{mn} - 0.5 * K_{mn}
//!
//! J_{mn} = sum_{ls} P_{ls} (mn|ls)     [Coulomb: classical repulsion]
//! K_{mn} = sum_{ls} P_{ls} (ml|ns)     [Exchange: quantum effect, no classical analogue]
//! ```
//!
//! where P is the density matrix with factor of 2 for RHF closed-shell:
//!   P = 2 * C_occ * C_occ^T
//!
//! # References
//!
//! - Szabo & Ostlund (1996), *Modern Quantum Chemistry*, Eq. 3.154
//! - PySCF: `hf.py` get_jk (separate J, K computation)
//! - IQCP: `scf/mod.rs` build_fock() (combined G computation)

use serde::{Deserialize, Serialize};

use super::eri_get;

/// Result of decomposing the Fock matrix into its physical components.
///
/// All matrices are stored as flat `Vec<f64>` in row-major order, dimension `nbf x nbf`.
///
/// Convention: Szabo & Ostlund (1996), Eq. 3.154
///   F = H^core + G, where G = J - 0.5*K
///   J_{mn} = sum_{ls} P_{ls} (mn|ls)     [Coulomb]
///   K_{mn} = sum_{ls} P_{ls} (ml|ns)     [Exchange]
///   P = D = 2*C_occ*C_occ^T              [density matrix with factor of 2]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FockDecomposition {
    /// Core Hamiltonian H^core = T + V (one-electron part)
    pub h_core: Vec<f64>,
    /// Coulomb matrix J: J_{mn} = sum_{ls} P_{ls} (mn|ls)
    pub j_matrix: Vec<f64>,
    /// Exchange matrix K: K_{mn} = sum_{ls} P_{ls} (ml|ns)
    pub k_matrix: Vec<f64>,
    /// Two-electron contribution G = J - 0.5*K
    pub g_matrix: Vec<f64>,
    /// Full Fock matrix F = H^core + G
    pub f_matrix: Vec<f64>,
    /// Density matrix P (input, echoed for display; includes factor of 2)
    pub density: Vec<f64>,
    /// Number of basis functions
    pub nbf: usize,
    /// Basis function labels (e.g., ["H1 1s", "H2 1s"])
    pub labels: Vec<String>,
}

/// Compute Fock matrix decomposition with separate J and K matrices.
///
/// This function is designed for educational inspection, not for the
/// SCF hot path. It performs the same computation as `build_fock()`
/// in `scf/mod.rs` but stores J and K in separate matrices for transparency.
///
/// # Arguments
///
/// * `h_core` - Core Hamiltonian matrix (nbf x nbf, flat row-major)
/// * `eri_compressed` - 8-fold symmetry compressed ERI tensor
/// * `density` - Density matrix P = 2*C_occ*C_occ^T (nbf x nbf, flat row-major)
/// * `nbf` - Number of basis functions
/// * `labels` - Basis function labels
///
/// # Returns
///
/// `FockDecomposition` with H^core, J, K, G, F matrices
///
/// # Algorithm
///
/// For each pair (mu, nu), accumulate over all (lambda, sigma):
///
/// ```text
/// J_{mu,nu} += P_{lambda,sigma} * (mu nu | lambda sigma)
/// K_{mu,nu} += P_{lambda,sigma} * (mu lambda | nu sigma)
/// ```
///
/// Then compute G = J - 0.5*K and F = H^core + G.
///
/// This is identical to `build_fock()` but separates J and K storage.
///
/// # References
///
/// - Szabo & Ostlund (1996), Eq. 3.154
/// - PySCF: hf.py get_jk (separate J,K computation)
/// - IQCP: scf/mod.rs build_fock() (combined G computation)
pub fn compute_fock_decomposition(
    h_core: &[f64],
    eri_compressed: &[f64],
    density: &[f64],
    nbf: usize,
    labels: Vec<String>,
) -> FockDecomposition {
    let n2 = nbf * nbf;

    // Initialize J and K matrices to zero
    let mut j_matrix = vec![0.0_f64; n2];
    let mut k_matrix = vec![0.0_f64; n2];

    // Accumulate J and K contributions.
    //
    // Exploits μ >= ν and λ >= σ symmetry (~4x speedup).
    // Same algorithm as build_fock() in scf/mod.rs, but stores J and K
    // separately for educational inspection.
    //
    // Szabo & Ostlund Eq. 3.154:
    //   J_{mu,nu} = sum_{lambda,sigma} P_{lambda,sigma} * (mu nu | lambda sigma)
    //   K_{mu,nu} = sum_{lambda,sigma} P_{lambda,sigma} * (mu lambda | nu sigma)
    for mu in 0..nbf {
        for nu in 0..=mu {
            let mut j_mn = 0.0;
            let mut k_mn = 0.0;

            for lambda in 0..nbf {
                // Diagonal: σ = λ
                {
                    let d_ll = density[lambda * nbf + lambda];
                    let j_integral = eri_get(eri_compressed, mu, nu, lambda, lambda);
                    let k_integral = eri_get(eri_compressed, mu, lambda, nu, lambda);
                    j_mn += d_ll * j_integral;
                    k_mn += d_ll * k_integral;
                }
                // Off-diagonal: σ < λ
                for sigma in 0..lambda {
                    let d_ls = density[lambda * nbf + sigma];
                    let j_integral = eri_get(eri_compressed, mu, nu, lambda, sigma);
                    j_mn += 2.0 * d_ls * j_integral;
                    let k_integral_1 = eri_get(eri_compressed, mu, lambda, nu, sigma);
                    let k_integral_2 = eri_get(eri_compressed, mu, sigma, nu, lambda);
                    k_mn += d_ls * (k_integral_1 + k_integral_2);
                }
            }

            j_matrix[mu * nbf + nu] = j_mn;
            k_matrix[mu * nbf + nu] = k_mn;
            if mu != nu {
                j_matrix[nu * nbf + mu] = j_mn;
                k_matrix[nu * nbf + mu] = k_mn;
            }
        }
    }

    // Compute G = J - 0.5*K and F = H^core + G
    let mut g_matrix = vec![0.0_f64; n2];
    let mut f_matrix = vec![0.0_f64; n2];

    for i in 0..n2 {
        g_matrix[i] = j_matrix[i] - 0.5 * k_matrix[i];
        f_matrix[i] = h_core[i] + g_matrix[i];
    }

    FockDecomposition {
        h_core: h_core.to_vec(),
        j_matrix,
        k_matrix,
        g_matrix,
        f_matrix,
        density: density.to_vec(),
        nbf,
        labels,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scf::PresetSystem;

    // ========================================================================
    // PySCF Golden Reference Values (H2 STO-3G, R=1.4 bohr)
    // Generated: PySCF 2.11.0, conv_tol=1e-12
    // ========================================================================

    /// PySCF reference J matrix for H2 STO-3G at R=1.4 bohr
    const PYSCF_J: [f64; 4] = [
        1.345430416741971e+00,
        8.933020749953189e-01,
        8.933020749953189e-01,
        1.345430416741970e+00,
    ];

    /// PySCF reference K matrix for H2 STO-3G at R=1.4 bohr
    const PYSCF_K: [f64; 4] = [
        1.181117517432615e+00,
        1.057614974304675e+00,
        1.057614974304675e+00,
        1.181117517432615e+00,
    ];

    /// PySCF reference G = J - 0.5*K for H2 STO-3G at R=1.4 bohr
    const PYSCF_G: [f64; 4] = [
        7.548716580256631e-01,
        3.644945878429816e-01,
        3.644945878429816e-01,
        7.548716580256630e-01,
    ];

    /// PySCF reference F = H^core + G for H2 STO-3G at R=1.4 bohr
    const PYSCF_F: [f64; 4] = [
        -3.655373508811572e-01,
        -5.938853765466354e-01,
        -5.938853765466354e-01,
        -3.655373508811574e-01,
    ];

    /// PySCF reference density matrix P (includes factor of 2)
    const PYSCF_D: [f64; 4] = [
        6.026571614189372e-01,
        6.026571614189370e-01,
        6.026571614189370e-01,
        6.026571614189368e-01,
    ];

    /// Helper: Get converged density matrix for H2 STO-3G test system.
    ///
    /// Runs SCF to convergence and returns the density matrix as flat row-major.
    fn h2_converged_density() -> (PresetSystem, Vec<f64>) {
        use crate::scf::{rhf_scf, ConvergenceProfile, ScfConfig};

        let system = PresetSystem::h2_sto3g_test();
        let config = ScfConfig::new(ConvergenceProfile::Tight);
        let output = rhf_scf(&system, &config).expect("H2 SCF should converge");
        assert!(
            output.converged,
            "H2 SCF must converge for Fock decomposition test"
        );

        let density = output.density_matrix;
        (system, density)
    }

    // ========================================================================
    // FB-R1: J matrix symmetry
    // ========================================================================

    #[test]
    fn test_j_matrix_symmetry() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        let nbf = decomp.nbf;
        for i in 0..nbf {
            for j in 0..nbf {
                let j_ij = decomp.j_matrix[i * nbf + j];
                let j_ji = decomp.j_matrix[j * nbf + i];
                assert!(
                    (j_ij - j_ji).abs() < 1e-15,
                    "J[{},{}] = {:.15e} != J[{},{}] = {:.15e}",
                    i,
                    j,
                    j_ij,
                    j,
                    i,
                    j_ji
                );
            }
        }
    }

    // ========================================================================
    // FB-R2: K matrix symmetry
    // ========================================================================

    #[test]
    fn test_k_matrix_symmetry() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        let nbf = decomp.nbf;
        for i in 0..nbf {
            for j in 0..nbf {
                let k_ij = decomp.k_matrix[i * nbf + j];
                let k_ji = decomp.k_matrix[j * nbf + i];
                assert!(
                    (k_ij - k_ji).abs() < 1e-15,
                    "K[{},{}] = {:.15e} != K[{},{}] = {:.15e}",
                    i,
                    j,
                    k_ij,
                    j,
                    i,
                    k_ji
                );
            }
        }
    }

    // ========================================================================
    // FB-R3: G = J - 0.5*K identity
    // ========================================================================

    #[test]
    fn test_g_equals_j_minus_half_k() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        for i in 0..decomp.nbf * decomp.nbf {
            let expected = decomp.j_matrix[i] - 0.5 * decomp.k_matrix[i];
            assert!(
                (decomp.g_matrix[i] - expected).abs() < 1e-15,
                "G[{}] = {:.15e} != J[{}] - 0.5*K[{}] = {:.15e}",
                i,
                decomp.g_matrix[i],
                i,
                i,
                expected
            );
        }
    }

    // ========================================================================
    // FB-R4: F = H^core + G identity
    // ========================================================================

    #[test]
    fn test_f_equals_hcore_plus_g() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        for i in 0..decomp.nbf * decomp.nbf {
            let expected = decomp.h_core[i] + decomp.g_matrix[i];
            assert!(
                (decomp.f_matrix[i] - expected).abs() < 1e-15,
                "F[{}] = {:.15e} != H_core[{}] + G[{}] = {:.15e}",
                i,
                decomp.f_matrix[i],
                i,
                i,
                expected
            );
        }
    }

    // ========================================================================
    // FB-R5: F matches build_fock() from SCF engine
    // ========================================================================

    #[test]
    fn test_f_matches_build_fock() {
        use crate::scf::{build_fock, rhf_scf, ConvergenceProfile, ScfConfig};
        use nalgebra::DMatrix;

        let system = PresetSystem::h2_sto3g_test();
        let config = ScfConfig::new(ConvergenceProfile::Tight);
        let output = rhf_scf(&system, &config).expect("H2 SCF should converge");

        let nbf = system.nbf;
        let density = &output.density_matrix;

        // Build Fock using the SCF engine's build_fock
        let h_core_mat = DMatrix::from_row_slice(nbf, nbf, &system.h_core);
        let density_mat = DMatrix::from_row_slice(nbf, nbf, density);
        let fock_scf = build_fock(&h_core_mat, &density_mat, &system.eri_compressed, nbf);

        // Build Fock using decomposition
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];
        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            density,
            nbf,
            labels,
        );

        // Compare element by element
        for i in 0..nbf {
            for j in 0..nbf {
                let f_decomp = decomp.f_matrix[i * nbf + j];
                // build_fock uses nalgebra column-major: (row, col) indexing
                let f_scf = fock_scf[(i, j)];
                assert!(
                    (f_decomp - f_scf).abs() < 1e-12,
                    "F_decomp[{},{}] = {:.15e} != F_scf[{},{}] = {:.15e}, diff = {:.2e}",
                    i,
                    j,
                    f_decomp,
                    i,
                    j,
                    f_scf,
                    (f_decomp - f_scf).abs()
                );
            }
        }
    }

    // ========================================================================
    // FB-R6: Correct dimensions
    // ========================================================================

    #[test]
    fn test_correct_dimensions() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        let n2 = system.nbf * system.nbf;
        assert_eq!(decomp.nbf, 2);
        assert_eq!(decomp.h_core.len(), n2);
        assert_eq!(decomp.j_matrix.len(), n2);
        assert_eq!(decomp.k_matrix.len(), n2);
        assert_eq!(decomp.g_matrix.len(), n2);
        assert_eq!(decomp.f_matrix.len(), n2);
        assert_eq!(decomp.density.len(), n2);
    }

    // ========================================================================
    // FB-R7: Labels correct count
    // ========================================================================

    #[test]
    fn test_labels_correct_count() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        assert_eq!(decomp.labels.len(), decomp.nbf);
    }

    // ========================================================================
    // FB-G1: J matrix matches PySCF (golden test)
    // ========================================================================

    #[test]
    fn test_j_matrix_matches_pyscf() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        for i in 0..4 {
            assert!(
                (decomp.j_matrix[i] - PYSCF_J[i]).abs() < 1e-10,
                "J[{}]: IQCP = {:.15e}, PySCF = {:.15e}, diff = {:.2e}",
                i,
                decomp.j_matrix[i],
                PYSCF_J[i],
                (decomp.j_matrix[i] - PYSCF_J[i]).abs()
            );
        }
    }

    // ========================================================================
    // FB-G2: K matrix matches PySCF (golden test)
    // ========================================================================

    #[test]
    fn test_k_matrix_matches_pyscf() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        for i in 0..4 {
            assert!(
                (decomp.k_matrix[i] - PYSCF_K[i]).abs() < 1e-10,
                "K[{}]: IQCP = {:.15e}, PySCF = {:.15e}, diff = {:.2e}",
                i,
                decomp.k_matrix[i],
                PYSCF_K[i],
                (decomp.k_matrix[i] - PYSCF_K[i]).abs()
            );
        }
    }

    // ========================================================================
    // FB-G3: F matrix matches PySCF (golden test)
    // ========================================================================

    #[test]
    fn test_f_matrix_matches_pyscf() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        for i in 0..4 {
            assert!(
                (decomp.f_matrix[i] - PYSCF_F[i]).abs() < 1e-10,
                "F[{}]: IQCP = {:.15e}, PySCF = {:.15e}, diff = {:.2e}",
                i,
                decomp.f_matrix[i],
                PYSCF_F[i],
                (decomp.f_matrix[i] - PYSCF_F[i]).abs()
            );
        }
    }

    // ========================================================================
    // FB-P1/P2: All elements finite
    // ========================================================================

    #[test]
    fn test_all_elements_finite() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        for (name, matrix) in &[
            ("J", &decomp.j_matrix),
            ("K", &decomp.k_matrix),
            ("G", &decomp.g_matrix),
            ("F", &decomp.f_matrix),
        ] {
            for (i, val) in matrix.iter().enumerate() {
                assert!(val.is_finite(), "{}[{}] = {} is not finite", name, i, val);
            }
        }
    }

    // ========================================================================
    // FB-P3: J diagonal >= 0 (Coulomb self-interaction is non-negative)
    // ========================================================================

    #[test]
    fn test_j_diagonal_non_negative() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        for i in 0..decomp.nbf {
            let j_ii = decomp.j_matrix[i * decomp.nbf + i];
            assert!(
                j_ii >= 0.0,
                "J[{},{}] = {:.15e} should be non-negative",
                i,
                i,
                j_ii
            );
        }
    }

    // ========================================================================
    // FB-P4: Tr(P*J) > 0 (total Coulomb energy is positive/repulsive)
    // ========================================================================

    #[test]
    fn test_coulomb_energy_positive() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        let n2 = decomp.nbf * decomp.nbf;
        let tr_pj: f64 = (0..n2)
            .map(|i| decomp.density[i] * decomp.j_matrix[i])
            .sum();
        assert!(
            tr_pj > 0.0,
            "Tr(P*J) = {:.15e} should be positive (Coulomb repulsion)",
            tr_pj
        );
    }

    // ========================================================================
    // FB-P5: Tr(P*K) > 0 (exchange energy is positive before -0.5 factor)
    // ========================================================================

    #[test]
    fn test_exchange_energy_positive() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        let n2 = decomp.nbf * decomp.nbf;
        let tr_pk: f64 = (0..n2)
            .map(|i| decomp.density[i] * decomp.k_matrix[i])
            .sum();
        assert!(tr_pk > 0.0, "Tr(P*K) = {:.15e} should be positive", tr_pk);
    }

    // ========================================================================
    // FB-G*: G and density matrix match PySCF (additional golden checks)
    // ========================================================================

    #[test]
    fn test_g_matrix_matches_pyscf() {
        let (system, density) = h2_converged_density();
        let labels = vec!["H1 1s".to_string(), "H2 1s".to_string()];

        let decomp = compute_fock_decomposition(
            &system.h_core,
            &system.eri_compressed,
            &density,
            system.nbf,
            labels,
        );

        for i in 0..4 {
            assert!(
                (decomp.g_matrix[i] - PYSCF_G[i]).abs() < 1e-10,
                "G[{}]: IQCP = {:.15e}, PySCF = {:.15e}, diff = {:.2e}",
                i,
                decomp.g_matrix[i],
                PYSCF_G[i],
                (decomp.g_matrix[i] - PYSCF_G[i]).abs()
            );
        }
    }

    #[test]
    fn test_density_matches_pyscf() {
        let (_system, density) = h2_converged_density();

        // Our SCF density should match PySCF's within tolerance
        for i in 0..4 {
            assert!(
                (density[i] - PYSCF_D[i]).abs() < 1e-10,
                "D[{}]: IQCP = {:.15e}, PySCF = {:.15e}, diff = {:.2e}",
                i,
                density[i],
                PYSCF_D[i],
                (density[i] - PYSCF_D[i]).abs()
            );
        }
    }
}
