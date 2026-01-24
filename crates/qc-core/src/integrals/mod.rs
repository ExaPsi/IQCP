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

mod cartesian;
pub mod eri;
mod gaussian_product;
mod hcore;
mod kinetic;
mod nuclear;
mod overlap;
pub mod spherical;

pub use cartesian::{cartesian_components, CartesianPower};
pub use eri::{
    eri_compressed, eri_compressed_spherical, eri_get, eri_index, pair_index, primitive_eri,
    shell_eri, shell_eri_spherical, EriError, EriResult, EriSphericalResult, GaussianProduct2e,
    RysCoefficients, MAX_ERI_ANGULAR_MOMENTUM, MAX_TOTAL_ANGULAR_MOMENTUM,
};

#[cfg(feature = "parallel")]
pub use eri::{eri_compressed_parallel, eri_compressed_spherical_parallel};
pub use gaussian_product::GaussianProduct;
pub use hcore::{hcore_matrix, hcore_matrix_spherical};
pub use kinetic::{kinetic_matrix, kinetic_matrix_spherical, primitive_kinetic, shell_kinetic};
pub use nuclear::{nuclear_matrix, nuclear_matrix_spherical, primitive_nuclear, shell_nuclear};
pub use overlap::{
    overlap_matrix, overlap_matrix_spherical, primitive_overlap, pyscf_eri_normalization,
    radial_gto_norm, shell_overlap, shell_renorm_factor,
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
}
