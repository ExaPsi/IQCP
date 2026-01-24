//! Core Hamiltonian assembly: H^core = T + V
//!
//! This module assembles the one-electron core Hamiltonian matrix by combining
//! the kinetic energy (T) and nuclear attraction (V) matrices:
//!
//! ```text
//! H^core_μν = T_μν + V_μν
//!           = ⟨μ|-½∇²|ν⟩ + Σ_A (-Z_A) ⟨μ|1/|r-R_A||ν⟩
//! ```
//!
//! # Physical Meaning
//!
//! The core Hamiltonian represents the interaction of electrons with:
//! - Their own kinetic energy (T): momentum of the electron
//! - The nuclear electrostatic potential (V): attraction to all nuclei
//!
//! It does NOT include electron-electron repulsion, which is handled separately
//! in the Fock matrix construction via the two-electron integral contributions.
//!
//! # Properties
//!
//! - **Symmetric:** H^core_μν = H^core_νμ (inherits from T and V symmetry)
//! - **Hermitian:** Real-valued for real basis functions
//! - **Negative off-diagonal:** Typically, off-diagonal elements are negative due to V
//! - **Diagonal dominance:** Diagonal elements have large negative contributions from V
//!
//! # References
//!
//! - Szabo & Ostlund, Modern Quantum Chemistry, Eq. 3.148
//! - Helgaker, Jorgensen & Olsen (2000), Molecular Electronic-Structure Theory, Ch. 10
//! - PySCF: `references/pyscf/pyscf/scf/hf.py` (get_hcore function)
//!
//! # Example
//!
//! ```rust
//! use qc_core::basis::{Atom, BasisSet};
//! use qc_core::integrals::{hcore_matrix, kinetic_matrix, nuclear_matrix};
//!
//! // Build H2 molecule
//! let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
//! let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
//! let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
//!
//! // Compute H^core matrix
//! let hcore = hcore_matrix(&basis);
//!
//! // H^core is a 2x2 symmetric matrix (flattened row-major)
//! assert_eq!(hcore.len(), 4);
//!
//! // Verify H^core = T + V
//! let t = kinetic_matrix(&basis);
//! let v = nuclear_matrix(&basis);
//! for i in 0..4 {
//!     assert!((hcore[i] - (t[i] + v[i])).abs() < 1e-15);
//! }
//! ```

use super::{kinetic_matrix, kinetic_matrix_spherical, nuclear_matrix, nuclear_matrix_spherical};
use crate::basis::BasisSet;

// =============================================================================
// Core Hamiltonian Matrix
// =============================================================================

/// Compute the full core Hamiltonian matrix for a basis set
///
/// Returns the N x N core Hamiltonian matrix H^core where N is the total number
/// of basis functions. The matrix is symmetric: H^core_μν = H^core_νμ.
///
/// The core Hamiltonian is the sum of kinetic and nuclear attraction matrices:
/// ```text
/// H^core = T + V
/// ```
///
/// # Arguments
///
/// * `basis` - The molecular basis set (contains atoms with positions and charges)
///
/// # Returns
///
/// Vector of length N*N containing the flattened row-major core Hamiltonian matrix
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{Atom, BasisSet};
/// use qc_core::integrals::hcore_matrix;
///
/// // H2 molecule
/// let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
/// let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
/// let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
///
/// let hcore = hcore_matrix(&basis);
/// assert_eq!(hcore.len(), 4);  // 2x2 matrix
///
/// // Check symmetry
/// assert!((hcore[1] - hcore[2]).abs() < 1e-15);
///
/// // Diagonal elements are typically negative (V dominates T)
/// // For small molecules like H2, this may vary
/// ```
pub fn hcore_matrix(basis: &BasisSet) -> Vec<f64> {
    // Compute kinetic energy matrix T
    let t_matrix = kinetic_matrix(basis);

    // Compute nuclear attraction matrix V
    let v_matrix = nuclear_matrix(basis);

    // Compute H^core = T + V (element-wise addition)
    t_matrix
        .iter()
        .zip(v_matrix.iter())
        .map(|(t, v)| t + v)
        .collect()
}

/// Compute the core Hamiltonian matrix in spherical harmonic basis
///
/// This computes H^core = T + V in the spherical harmonic basis.
/// For basis sets without d-orbitals or higher, this is identical to
/// `hcore_matrix()`.
///
/// # Arguments
///
/// * `basis` - The molecular basis set
///
/// # Returns
///
/// Vector of length N_sph * N_sph containing the flattened row-major
/// spherical core Hamiltonian matrix, where N_sph is the number of
/// spherical basis functions.
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{Atom, BasisSet};
/// use qc_core::integrals::hcore_matrix_spherical;
///
/// // H2O with 6-31G* (has d-polarization on O)
/// let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
/// let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();
/// let h2 = Atom::new(1, [0.0, -1.43, -0.47]).unwrap();
///
/// let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();
/// let hcore_sph = hcore_matrix_spherical(&basis);
///
/// assert_eq!(hcore_sph.len(), basis.n_basis_spherical() * basis.n_basis_spherical());
/// ```
pub fn hcore_matrix_spherical(basis: &BasisSet) -> Vec<f64> {
    // Compute kinetic energy matrix T in spherical basis
    let t_matrix = kinetic_matrix_spherical(basis);

    // Compute nuclear attraction matrix V in spherical basis
    let v_matrix = nuclear_matrix_spherical(basis);

    // Compute H^core = T + V (element-wise addition)
    t_matrix
        .iter()
        .zip(v_matrix.iter())
        .map(|(t, v)| t + v)
        .collect()
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::Atom;
    use crate::integrals::{kinetic_matrix as t_matrix_fn, nuclear_matrix as v_matrix_fn};
    use approx::assert_abs_diff_eq;

    // Tolerance for golden test comparisons
    // H^core inherits tolerances from T and V matrices.
    // The nuclear attraction integral (V) uses recursive VRR/HTR with Boys functions,
    // which accumulates errors resulting in ~1e-5 to 1e-6 relative error.
    // For H2 with s-orbitals, the error is around 1e-6.
    // For H2O with p-orbitals, the accumulated error is ~1e-5.
    const TOL_H2: f64 = 1e-5;
    const TOL_H2O: f64 = 1e-5;

    // -------------------------------------------------------------------------
    // Golden tests: H2 STO-3G
    // -------------------------------------------------------------------------

    #[test]
    fn test_golden_h2_sto3g_hcore() {
        // Reference: PySCF 2.11.0
        //
        // H2 at 1.3984 Bohr separation
        // Expected H^core matrix (2x2):
        //
        // Generated with:
        // ```python
        // from pyscf import gto, scf
        // mol = gto.Mole()
        // mol.atom = 'H 0 0 0; H 0 0 1.3984'
        // mol.basis = 'sto-3g'
        // mol.unit = 'B'
        // mol.build()
        // mf = scf.RHF(mol)
        // print(mf.get_hcore())
        // ```
        //
        // Reference values:
        // H[0,0] = -1.120958539482569e+00
        // H[0,1] = -9.593741137393389e-01
        // H[1,0] = -9.593741137393389e-01
        // H[1,1] = -1.120958539482569e+00

        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let hcore = hcore_matrix(&basis);

        // Reference values from PySCF
        let h_00 = -1.120958539482569e+00;
        let h_01 = -9.593741137393389e-01;

        assert_abs_diff_eq!(hcore[0], h_00, epsilon = TOL_H2); // H[0,0]
        assert_abs_diff_eq!(hcore[1], h_01, epsilon = TOL_H2); // H[0,1]
        assert_abs_diff_eq!(hcore[2], h_01, epsilon = TOL_H2); // H[1,0] (symmetry)
        assert_abs_diff_eq!(hcore[3], h_00, epsilon = TOL_H2); // H[1,1]
    }

    // -------------------------------------------------------------------------
    // Golden tests: H2O STO-3G
    // -------------------------------------------------------------------------

    #[test]
    fn test_golden_h2o_sto3g_hcore() {
        // Reference: PySCF 2.11.0
        //
        // H2O geometry (Bohr):
        // O: [0, 0, 0.2216656]
        // H1: [0, 1.4309295, -0.8866625]
        // H2: [0, -1.4309295, -0.8866625]
        //
        // 7x7 H^core matrix
        //
        // Generated with:
        // ```python
        // from pyscf import gto, scf
        // mol = gto.Mole()
        // mol.atom = '''
        // O   0.0000000   0.0000000   0.2216656303019316
        // H   0.0000000   1.4309295343303710  -0.8866625212077263
        // H   0.0000000  -1.4309295343303710  -0.8866625212077263
        // '''
        // mol.basis = 'sto-3g'
        // mol.unit = 'B'
        // mol.build()
        // mf = scf.RHF(mol)
        // print(repr(mf.get_hcore().flatten()))
        // ```

        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        let hcore = hcore_matrix(&basis);
        assert_eq!(hcore.len(), 49); // 7x7

        // Reference values from PySCF 2.11.0 (flattened row-major)
        // H^core = T + V
        #[rustfmt::skip]
        let reference = [
            -32.72025252747679, -7.612647658092316, 0.0, 0.0, 0.01898834225829542,
            -1.7470967486879878, -1.7470967486879878,
            -7.612647658092317, -9.334210191717451, 0.0, 0.0, 0.22339810110550098,
            -3.737626039294491, -3.737626039294491,
            0.0, 0.0, -7.457158706515188, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -7.613191451727299, 0.0, -2.02824154714686, 2.02824154714686,
            0.018988342258295425, 0.223398101105501, 0.0, 0.0, -7.550767356765422,
            1.6438933044362245, 1.6438933044362245,
            -1.7470967486879876, -3.7376260392944913, 0.0, -2.0282415471468607,
            1.643893304436224, -5.074340305149788, -1.6064807356381123,
            -1.7470967486879876, -3.7376260392944913, 0.0, 2.0282415471468607,
            1.643893304436224, -1.6064807356381123, -5.074340305149788,
        ];

        // Check all elements
        for (i, &ref_val) in reference.iter().enumerate() {
            assert!(
                (hcore[i] - ref_val).abs() < TOL_H2O,
                "Mismatch at index {}: computed {} vs reference {} (error: {:.2e})",
                i,
                hcore[i],
                ref_val,
                (hcore[i] - ref_val).abs()
            );
        }
    }

    // -------------------------------------------------------------------------
    // H^core = T + V verification
    // -------------------------------------------------------------------------

    #[test]
    fn test_hcore_equals_t_plus_v_h2() {
        // Explicitly verify H^core = T + V for H2
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let hcore = hcore_matrix(&basis);
        let t = t_matrix_fn(&basis);
        let v = v_matrix_fn(&basis);

        for i in 0..hcore.len() {
            let expected = t[i] + v[i];
            assert_abs_diff_eq!(hcore[i], expected, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_hcore_equals_t_plus_v_h2o() {
        // Explicitly verify H^core = T + V for H2O
        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        let hcore = hcore_matrix(&basis);
        let t = t_matrix_fn(&basis);
        let v = v_matrix_fn(&basis);

        for i in 0..hcore.len() {
            let expected = t[i] + v[i];
            assert_abs_diff_eq!(hcore[i], expected, epsilon = 1e-15);
        }
    }

    // -------------------------------------------------------------------------
    // Matrix properties tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_hcore_matrix_symmetry() {
        // H^core matrix should be symmetric: H^core_μν = H^core_νμ
        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.43, -0.47]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        let hcore = hcore_matrix(&basis);
        let n = basis.n_basis;

        for i in 0..n {
            for j in 0..n {
                assert_abs_diff_eq!(hcore[i * n + j], hcore[j * n + i], epsilon = 1e-14);
            }
        }
    }

    #[test]
    fn test_hcore_matrix_size() {
        // Verify matrix size matches basis set size
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let hcore = hcore_matrix(&basis);
        let n = basis.n_basis;

        assert_eq!(hcore.len(), n * n);
    }

    #[test]
    fn test_hcore_h2o_matrix_size() {
        // Verify H2O matrix size: 7 basis functions for STO-3G
        // O: 1s, 2s, 2px, 2py, 2pz = 5 functions
        // H: 1s = 1 function each, 2 H atoms = 2 functions
        // Total = 7
        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.43, -0.47]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        assert_eq!(basis.n_basis, 7);

        let hcore = hcore_matrix(&basis);
        assert_eq!(hcore.len(), 49); // 7x7
    }

    #[test]
    fn test_hcore_diagonal_properties() {
        // The diagonal elements of H^core are typically negative because
        // the nuclear attraction (V, negative) usually dominates kinetic energy (T, positive)
        // For most basis functions, especially core orbitals, this holds.
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let hcore = hcore_matrix(&basis);
        let n = basis.n_basis;

        // For H2 with STO-3G, diagonal elements should be negative
        for i in 0..n {
            let diag = hcore[i * n + i];
            assert!(
                diag < 0.0,
                "Diagonal element H[{},{}] should be negative, got {}",
                i,
                i,
                diag
            );
        }
    }

    // -------------------------------------------------------------------------
    // Spherical harmonic transformation tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_hcore_spherical_dimension_reduction_6_31g_star() {
        // Test dimension reduction: 6-31G* has d-orbitals which reduce from 6 to 5 functions
        //
        // For H2O with 6-31G*:
        // O: 1s(1), 2s(1), 3s(1), 2px(1), 2py(1), 2pz(1), 3px(1), 3py(1), 3pz(1), d(6 cart/5 sph)
        //    = 9 + 6 = 15 Cartesian, 9 + 5 = 14 spherical
        // H: 1s(1), 2s(1) = 2 each, total 4
        // Total: Cartesian = 15 + 4 = 19, Spherical = 14 + 4 = 18

        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();

        // Verify Cartesian vs spherical dimension
        let n_cart = basis.n_basis;
        let n_sph = basis.n_basis_spherical();

        assert_eq!(
            n_cart, 19,
            "H2O 6-31G* should have 19 Cartesian basis functions"
        );
        assert_eq!(
            n_sph, 18,
            "H2O 6-31G* should have 18 spherical basis functions"
        );
        assert!(
            basis.has_spherical_difference(),
            "Basis should have d-orbitals"
        );

        // Compute matrices
        let hcore_cart = hcore_matrix(&basis);
        let hcore_sph = hcore_matrix_spherical(&basis);

        // Verify dimensions
        assert_eq!(hcore_cart.len(), n_cart * n_cart);
        assert_eq!(hcore_sph.len(), n_sph * n_sph);
    }

    #[test]
    fn test_hcore_spherical_no_change_for_sp_only() {
        // For STO-3G (s and p orbitals only), spherical and Cartesian are identical
        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        // STO-3G has no d-orbitals
        assert!(!basis.has_spherical_difference());
        assert_eq!(basis.n_basis, basis.n_basis_spherical());

        let hcore_cart = hcore_matrix(&basis);
        let hcore_sph = hcore_matrix_spherical(&basis);

        // Matrices should be identical
        assert_eq!(hcore_cart.len(), hcore_sph.len());
        for (cart, sph) in hcore_cart.iter().zip(hcore_sph.iter()) {
            assert_abs_diff_eq!(cart, sph, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_hcore_spherical_symmetry() {
        // Spherical H^core matrix should remain symmetric
        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();

        let hcore_sph = hcore_matrix_spherical(&basis);
        let n = basis.n_basis_spherical();

        for i in 0..n {
            for j in 0..n {
                let diff = (hcore_sph[i * n + j] - hcore_sph[j * n + i]).abs();
                assert!(
                    diff < 1e-14,
                    "Spherical H^core not symmetric at ({}, {}): {} vs {}",
                    i,
                    j,
                    hcore_sph[i * n + j],
                    hcore_sph[j * n + i]
                );
            }
        }
    }

    #[test]
    fn test_hcore_spherical_equals_t_plus_v_spherical() {
        // Verify H^core_sph = T_sph + V_sph
        use crate::integrals::{kinetic_matrix_spherical, nuclear_matrix_spherical};

        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();

        let hcore_sph = hcore_matrix_spherical(&basis);
        let t_sph = kinetic_matrix_spherical(&basis);
        let v_sph = nuclear_matrix_spherical(&basis);

        for i in 0..hcore_sph.len() {
            let expected = t_sph[i] + v_sph[i];
            assert_abs_diff_eq!(hcore_sph[i], expected, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_hcore_spherical_diagonal_negative() {
        // Diagonal elements should remain negative after spherical transformation
        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();

        let hcore_sph = hcore_matrix_spherical(&basis);
        let n = basis.n_basis_spherical();

        for i in 0..n {
            let diag = hcore_sph[i * n + i];
            assert!(
                diag < 0.0,
                "Spherical H^core diagonal [{},{}] should be negative, got {}",
                i,
                i,
                diag
            );
        }
    }
}
