//! Molecular integral computation
//!
//! This module implements the computation of molecular integrals over
//! contracted Gaussian basis functions.
//!
//! # One-Electron Integrals (Obara-Saika)
//!
//! - **Overlap** `S_ij = ⟨i|j⟩` (US-027a) - Implemented
//! - **Kinetic** `T_ij = ⟨i|-½∇²|j⟩` (US-027b) - Implemented
//! - **Nuclear** `V_ij = ⟨i|-Z_C/|r-R_C||j⟩` (US-027b) - Implemented
//! - **H_core** `H^core_ij = T_ij + V_ij` (US-027b) - Implemented
//!
//! # Two-Electron Integrals (Rys Quadrature)
//!
//! - **ERI** `(ij|kl)` (US-028a) - Implemented for s,p angular momentum
//!
//! # Algorithm Overview
//!
//! The Obara-Saika method computes integrals through recurrence relations:
//!
//! 1. **Gaussian Product Theorem**: Two Gaussians at centers A and B combine to
//!    form a Gaussian at center P with exponent p = α + β
//!
//! 2. **Vertical Recurrence (VRR)**: Build angular momentum on center A
//!    ```text
//!    [a+1_i|b] = (P_i - A_i)[a|b] + (a_i/2p)[a-1_i|b] + (b_i/2p)[a|b-1_i]
//!    ```
//!
//! 3. **Horizontal Transfer (HTR)**: Transfer angular momentum from A to B
//!    ```text
//!    [a|b+1_i] = [a+1_i|b] + (A_i - B_i)[a|b]
//!    ```
//!
//! # Kinetic Energy Integrals
//!
//! The kinetic energy integral is computed using the relation:
//! ```text
//! T_1D(a, b) = -2β² S_1D(a, b+2) + β(2b+1) S_1D(a, b) - ½b(b-1) S_1D(a, b-2)
//! ```
//!
//! The 3D integral separates as: `T = T_x S_y S_z + S_x T_y S_z + S_x S_y T_z`
//!
//! # References
//!
//! - Obara & Saika (1986), J. Chem. Phys. 84, 3963
//! - Head-Gordon & Pople (1988), J. Chem. Phys. 89, 5777
//! - Helgaker, Jorgensen & Olsen (2000), Molecular Electronic-Structure Theory, Ch. 9
//! - libcint implementation: `references/libcint/src/g1e.c`
//!
//! # Example
//!
//! ```rust
//! use qc_core::integrals::{overlap_matrix, kinetic_matrix, nuclear_matrix, hcore_matrix};
//! use qc_core::basis::{Atom, BasisSet};
//!
//! // Build basis set for H2
//! let atoms = vec![
//!     Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
//!     Atom::new(1, [0.0, 0.0, 1.4]).unwrap(),
//! ];
//! let basis = BasisSet::build(atoms, "sto-3g").unwrap();
//!
//! // Compute overlap matrix
//! let s_matrix = overlap_matrix(&basis);
//! assert_eq!(s_matrix.len(), 4); // 2x2 flattened
//!
//! // Compute kinetic matrix
//! let t_matrix = kinetic_matrix(&basis);
//! assert_eq!(t_matrix.len(), 4); // 2x2 flattened
//!
//! // Compute nuclear attraction matrix
//! let v_matrix = nuclear_matrix(&basis);
//! assert_eq!(v_matrix.len(), 4); // 2x2 flattened
//!
//! // Compute core Hamiltonian H^core = T + V
//! let hcore = hcore_matrix(&basis);
//! assert_eq!(hcore.len(), 4); // 2x2 flattened
//! ```

pub mod breakdown;
mod cartesian;
pub mod deriv1;
pub mod dipole;
pub mod eri;
mod gaussian_product;
mod hcore;
mod kinetic;
mod nuclear;
pub(crate) mod overlap;
pub mod spherical;

pub use breakdown::{
    eri_with_breakdown, integral_with_breakdown, EriBreakdown, EriMethod, EriPrimitiveContribution,
    IntegralBreakdown, PrimitiveContribution,
};
pub use cartesian::{cartesian_components, CartesianPower};
pub use dipole::{dipole_matrix, shell_dipole_first_deriv, shell_dipole_integrals};
pub use eri::{
    eri_compressed, eri_compressed_spherical, eri_compressed_spherical_with_progress,
    eri_compressed_with_progress, eri_get, eri_index, pair_index, primitive_eri, shell_eri,
    shell_eri_gradient_direct, shell_eri_spherical, shell_eri_with_derivatives,
    shell_eri_with_second_derivatives, EriDerivResult, EriError, EriResult, EriSecondDerivResult,
    EriSphericalResult, GaussianProduct2e, RysCoefficients, MAX_ERI_ANGULAR_MOMENTUM,
    MAX_TOTAL_ANGULAR_MOMENTUM,
};

#[cfg(feature = "parallel")]
pub use eri::{eri_compressed_parallel, eri_compressed_spherical_parallel};
pub use gaussian_product::GaussianProduct;
pub use hcore::{hcore_matrix, hcore_matrix_spherical};
pub(crate) use kinetic::kinetic_1d;
pub use kinetic::{kinetic_matrix, kinetic_matrix_spherical, primitive_kinetic, shell_kinetic};
pub use nuclear::{nuclear_matrix, nuclear_matrix_spherical, primitive_nuclear, shell_nuclear};
pub use overlap::{
    evaluate_overlap_vs_distance, overlap_matrix, overlap_matrix_spherical, primitive_overlap,
    pyscf_eri_normalization, radial_gto_norm, shell_overlap, shell_renorm_factor,
};
pub use spherical::{
    transform_one_electron_matrix, EriSphericalTransform, OneElectronSphericalTransform,
    SphericalTransform, CART2SPH_D,
};

use thiserror::Error;

/// Errors that can occur during integral computation
#[derive(Error, Debug, Clone, PartialEq)]
pub enum IntegralError {
    /// Angular momentum exceeds supported maximum
    #[error("Angular momentum {0} exceeds maximum supported value {1}")]
    AngularMomentumTooHigh(u32, u32),

    /// Invalid basis set configuration
    #[error("Invalid basis set: {0}")]
    InvalidBasis(String),

    /// Numerical instability detected
    #[error("Numerical instability: {0}")]
    NumericalInstability(String),
}

/// Maximum supported angular momentum (d orbitals = 2)
pub const MAX_ANGULAR_MOMENTUM: u32 = 2;

/// Module version for compatibility tracking
pub const VERSION: &str = "0.1.0";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_valid() {
        // Verify VERSION matches expected semver pattern
        assert!(VERSION.starts_with("0."), "VERSION should start with '0.'");
    }

    #[test]
    fn max_angular_momentum_supports_d_orbitals() {
        assert_eq!(MAX_ANGULAR_MOMENTUM, 2);
    }

    /// Comprehensive accuracy test comparing IQCP on-the-fly integrals against PySCF 2.11.0.
    ///
    /// This test was instrumental in diagnosing and fixing two critical numerical issues:
    /// 1. The Cody (1969) erf approximation had only ~7-digit accuracy (replaced with libm::erf)
    /// 2. The Boys-function auxiliary-index nuclear attraction VRR suffered from catastrophic
    ///    cancellation when the nucleus was near a Gaussian center (replaced with Rys quadrature)
    ///
    /// All integrals now match PySCF within ~1e-9, limited by the 10-digit precision of the
    /// stored basis set coefficients.
    #[test]
    fn diagnostic_h2o_sto3g_integral_accuracy() {
        // Compare IQCP on-the-fly integrals against PySCF 2.11.0 reference values
        // at maximum precision to identify which integrals contribute to SCF energy error.
        //
        // PySCF reference generated with:
        //   mol.atom = 'O 0 0 0.2216656303019316; H 0 1.430929534330371 -0.8866625212077263; ...'
        //   mol.basis = 'sto-3g'; mol.unit = 'B'; mol.build()
        use crate::basis::{Atom, BasisSet};

        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();
        let n = basis.n_basis;
        assert_eq!(n, 7);

        let s = overlap_matrix(&basis);
        let t = kinetic_matrix(&basis);
        let v = nuclear_matrix(&basis);
        let h = hcore_matrix(&basis);
        let eri = eri_compressed(&basis);

        // Track max errors
        let mut max_s_err = 0.0_f64;
        let mut max_t_err = 0.0_f64;
        let mut max_v_err = 0.0_f64;
        let mut max_h_err = 0.0_f64;
        let mut max_eri_err = 0.0_f64;

        // ===== OVERLAP =====
        // PySCF reference (only non-zero upper triangle)
        let s_refs: &[(usize, usize, f64)] = &[
            (0, 0, 1.000000000000000000e+00),
            (0, 1, 2.367039365108475901e-01),
            (0, 5, 5.390224430387011156e-02),
            (0, 6, 5.390224430387011156e-02),
            (1, 1, 1.000000000000000000e+00),
            (1, 5, 4.744202424058437950e-01),
            (1, 6, 4.744202424058437950e-01),
            (2, 2, 1.000000000000000000e+00),
            (3, 3, 1.000000000000000000e+00),
            (3, 5, 3.108959396001613196e-01),
            (3, 6, -3.108959396001613196e-01),
            (4, 4, 1.000000000000000000e+00),
            (4, 5, -2.408048151792149538e-01),
            (4, 6, -2.408048151792149538e-01),
            (5, 5, 1.000000000000000222e+00),
            (5, 6, 2.515251126945288740e-01),
            (6, 6, 1.000000000000000222e+00),
        ];
        for &(i, j, ref_val) in s_refs {
            let err = (s[i * n + j] - ref_val).abs();
            if err > max_s_err {
                max_s_err = err;
            }
        }

        // ===== HCORE =====
        let h_refs: &[(usize, usize, f64)] = &[
            (0, 0, -3.272025252747678792e+01),
            (0, 1, -7.612647658092315694e+00),
            (0, 4, 1.898834225829542158e-02),
            (0, 5, -1.747096748687987811e+00),
            (0, 6, -1.747096748687987811e+00),
            (1, 1, -9.334210191717451366e+00),
            (1, 4, 2.233981011055009847e-01),
            (1, 5, -3.737626039294490887e+00),
            (1, 6, -3.737626039294490887e+00),
            (2, 2, -7.457158706515188307e+00),
            (3, 3, -7.613191451727298720e+00),
            (3, 5, -2.028241547146859780e+00),
            (3, 6, 2.028241547146859780e+00),
            (4, 4, -7.550767356765422100e+00),
            (4, 5, 1.643893304436224501e+00),
            (4, 6, 1.643893304436224501e+00),
            (5, 5, -5.074340305149788399e+00),
            (5, 6, -1.606480735638112334e+00),
            (6, 6, -5.074340305149788399e+00),
        ];
        for &(i, j, ref_val) in h_refs {
            let err = (h[i * n + j] - ref_val).abs();
            if err > max_h_err {
                max_h_err = err;
            }
        }

        // Also check T and V separately
        let t_ref_00 = 2.900319994553956700e+01;
        let v_ref_00 = -6.172345247301635400e+01;
        max_t_err = (t[0] - t_ref_00).abs();
        max_v_err = (v[0] - v_ref_00).abs();

        // ===== ERIs =====
        let eri_refs: &[(usize, usize, usize, usize, f64)] = &[
            (0, 0, 0, 0, 4.785065404705503234e+00),
            (1, 1, 1, 1, 8.172063215260582103e-01),
            (2, 2, 2, 2, 8.801590933750453871e-01),
            (3, 3, 3, 3, 8.801590933750453871e-01),
            (4, 4, 4, 4, 8.801590933750451651e-01),
            (5, 5, 5, 5, 7.746059439198977881e-01),
            (6, 6, 6, 6, 7.746059439198977881e-01),
        ];
        for &(i, j, k, l, ref_val) in eri_refs {
            let val = eri_get(&eri, n, i, j, k, l);
            let err = (val - ref_val).abs();
            if err > max_eri_err {
                max_eri_err = err;
            }
        }

        // Nuclear repulsion
        let e_nuc_ref = 9.189403756920761168e+00;
        let e_nuc_err = (basis.nuclear_repulsion - e_nuc_ref).abs();

        // Assert accuracy: all integrals should match PySCF within 5e-9
        // (limited by 10-digit precision of stored basis set coefficients)
        let integral_tol = 5e-9;
        assert!(
            max_s_err < integral_tol,
            "Overlap max error {:.2e} exceeds tolerance {:.2e}",
            max_s_err,
            integral_tol
        );
        assert!(
            max_h_err < integral_tol,
            "Hcore max error {:.2e} exceeds tolerance {:.2e}",
            max_h_err,
            integral_tol
        );
        assert!(
            max_eri_err < integral_tol,
            "ERI max error {:.2e} exceeds tolerance {:.2e}",
            max_eri_err,
            integral_tol
        );
        assert!(
            e_nuc_err < 1e-12,
            "E_nuc error {:.2e} exceeds tolerance 1e-12",
            e_nuc_err
        );
    }
}
