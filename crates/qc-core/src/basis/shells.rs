//! Contracted shell data structures
//!
//! This module defines the ContractedShell struct which represents a set of
//! basis functions sharing the same angular momentum and center but with
//! different contraction coefficients.
//!
//! # Background
//!
//! A contracted shell is a linear combination of primitive Gaussian functions:
//!
//! ```text
//! phi_mu(r) = sum_i c_i * g_i(alpha_i, r - R)
//! ```
//!
//! where `c_i` are contraction coefficients, `alpha_i` are exponents, and `R`
//! is the center (typically an atomic position).
//!
//! # References
//!
//! - Szabo & Ostlund (1996). "Modern Quantum Chemistry". Dover, Chapter 3.
//! - PySCF shell implementation: `references/pyscf/pyscf/gto/mole.py`

use super::primitives::{AngularMomentum, GaussianPrimitive};
use serde::{Deserialize, Serialize};
use std::fmt;

// =============================================================================
// ContractedShell
// =============================================================================

/// A contracted shell of Gaussian functions
///
/// A contracted shell consists of:
/// - A set of primitive Gaussians with different exponents
/// - Contraction coefficients for each primitive
/// - A center point in 3D space (typically an atomic nucleus)
/// - An angular momentum quantum number
///
/// All primitives in a shell share the same center and angular momentum.
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{ContractedShell, AngularMomentum, GaussianPrimitive};
///
/// // Create hydrogen 1s shell (STO-3G)
/// let primitives = vec![
///     GaussianPrimitive::new(3.4252509100, 0.1543289707),
///     GaussianPrimitive::new(0.6239137300, 0.5353281424),
///     GaussianPrimitive::new(0.1688554000, 0.4446345420),
/// ];
///
/// let shell = ContractedShell::new(
///     AngularMomentum::S,
///     primitives,
///     [0.0, 0.0, 0.0],  // Center at origin
///     0,                 // Atom index
/// );
///
/// assert_eq!(shell.n_basis_functions(), 1);  // S shell has 1 function
/// assert_eq!(shell.n_primitives(), 3);        // 3 primitives in contraction
/// ```
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ContractedShell {
    /// Angular momentum (S, P, D)
    pub angular_momentum: AngularMomentum,

    /// Primitive Gaussians in this contraction
    pub primitives: Vec<GaussianPrimitive>,

    /// Center coordinates in Bohr [x, y, z]
    pub center: [f64; 3],

    /// Index of the atom this shell belongs to
    pub atom_idx: usize,
}

impl ContractedShell {
    /// Create a new contracted shell
    ///
    /// # Arguments
    ///
    /// * `angular_momentum` - Angular momentum type (S, P, D)
    /// * `primitives` - Vector of primitive Gaussians
    /// * `center` - Center coordinates in Bohr [x, y, z]
    /// * `atom_idx` - Index of the owning atom
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::{ContractedShell, AngularMomentum, GaussianPrimitive};
    ///
    /// let prims = vec![
    ///     GaussianPrimitive::new(1.0, 0.5),
    ///     GaussianPrimitive::new(0.3, 0.5),
    /// ];
    /// let shell = ContractedShell::new(AngularMomentum::P, prims, [1.0, 0.0, 0.0], 0);
    ///
    /// assert_eq!(shell.n_basis_functions(), 3);
    /// ```
    pub fn new(
        angular_momentum: AngularMomentum,
        primitives: Vec<GaussianPrimitive>,
        center: [f64; 3],
        atom_idx: usize,
    ) -> Self {
        Self {
            angular_momentum,
            primitives,
            center,
            atom_idx,
        }
    }

    /// Number of basis functions in this shell (Cartesian)
    ///
    /// This is determined by the angular momentum:
    /// - S: 1 function
    /// - P: 3 functions (px, py, pz)
    /// - D: 6 functions (dxx, dxy, dxz, dyy, dyz, dzz)
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::{ContractedShell, AngularMomentum, GaussianPrimitive};
    ///
    /// let prims = vec![GaussianPrimitive::new(1.0, 1.0)];
    ///
    /// let s_shell = ContractedShell::new(AngularMomentum::S, prims.clone(), [0.0; 3], 0);
    /// let p_shell = ContractedShell::new(AngularMomentum::P, prims.clone(), [0.0; 3], 0);
    /// let d_shell = ContractedShell::new(AngularMomentum::D, prims, [0.0; 3], 0);
    ///
    /// assert_eq!(s_shell.n_basis_functions(), 1);
    /// assert_eq!(p_shell.n_basis_functions(), 3);
    /// assert_eq!(d_shell.n_basis_functions(), 6);
    /// ```
    #[inline]
    pub fn n_basis_functions(&self) -> usize {
        self.angular_momentum.n_cartesian()
    }

    /// Number of basis functions in this shell (spherical harmonics)
    ///
    /// This is determined by the angular momentum using 2l+1:
    /// - S: 1 function
    /// - P: 3 functions
    /// - D: 5 functions (d-2, d-1, d0, d+1, d+2)
    ///
    /// For S and P, this equals the Cartesian count. The difference
    /// appears at D (5 vs 6) and higher angular momenta.
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::{ContractedShell, AngularMomentum, GaussianPrimitive};
    ///
    /// let prims = vec![GaussianPrimitive::new(1.0, 1.0)];
    ///
    /// let s_shell = ContractedShell::new(AngularMomentum::S, prims.clone(), [0.0; 3], 0);
    /// let p_shell = ContractedShell::new(AngularMomentum::P, prims.clone(), [0.0; 3], 0);
    /// let d_shell = ContractedShell::new(AngularMomentum::D, prims, [0.0; 3], 0);
    ///
    /// assert_eq!(s_shell.n_basis_functions_spherical(), 1);
    /// assert_eq!(p_shell.n_basis_functions_spherical(), 3);
    /// assert_eq!(d_shell.n_basis_functions_spherical(), 5);
    /// ```
    #[inline]
    pub fn n_basis_functions_spherical(&self) -> usize {
        self.angular_momentum.n_spherical()
    }

    /// Number of primitive Gaussians in the contraction
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::{ContractedShell, AngularMomentum, GaussianPrimitive};
    ///
    /// let prims = vec![
    ///     GaussianPrimitive::new(1.0, 0.5),
    ///     GaussianPrimitive::new(0.3, 0.5),
    /// ];
    /// let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0; 3], 0);
    ///
    /// assert_eq!(shell.n_primitives(), 2);
    /// ```
    #[inline]
    pub fn n_primitives(&self) -> usize {
        self.primitives.len()
    }

    /// Get the l value (angular momentum quantum number)
    #[inline]
    pub fn l_value(&self) -> u32 {
        self.angular_momentum.l_value()
    }

    /// Check if this is an S shell
    #[inline]
    pub fn is_s_shell(&self) -> bool {
        self.angular_momentum == AngularMomentum::S
    }

    /// Check if this is a P shell
    #[inline]
    pub fn is_p_shell(&self) -> bool {
        self.angular_momentum == AngularMomentum::P
    }

    /// Check if this is a D shell
    #[inline]
    pub fn is_d_shell(&self) -> bool {
        self.angular_momentum == AngularMomentum::D
    }

    /// Check if the shell is valid (has primitives with positive exponents)
    pub fn is_valid(&self) -> bool {
        !self.primitives.is_empty() && self.primitives.iter().all(|p| p.is_valid())
    }

    /// Get the center x coordinate
    #[inline]
    pub fn x(&self) -> f64 {
        self.center[0]
    }

    /// Get the center y coordinate
    #[inline]
    pub fn y(&self) -> f64 {
        self.center[1]
    }

    /// Get the center z coordinate
    #[inline]
    pub fn z(&self) -> f64 {
        self.center[2]
    }

    /// Compute the maximum exponent in this shell
    ///
    /// Useful for screening in integral evaluation.
    pub fn max_exponent(&self) -> f64 {
        self.primitives
            .iter()
            .map(|p| p.exponent)
            .fold(0.0, f64::max)
    }

    /// Compute the minimum exponent in this shell
    ///
    /// Useful for estimating extent of the shell.
    pub fn min_exponent(&self) -> f64 {
        self.primitives
            .iter()
            .map(|p| p.exponent)
            .fold(f64::INFINITY, f64::min)
    }
}

impl fmt::Display for ContractedShell {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}({}) @ [{:.4}, {:.4}, {:.4}] (atom {})",
            self.angular_momentum,
            self.n_primitives(),
            self.center[0],
            self.center[1],
            self.center[2],
            self.atom_idx
        )
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn make_test_primitives() -> Vec<GaussianPrimitive> {
        vec![
            GaussianPrimitive::new(3.4252509100, 0.1543289707),
            GaussianPrimitive::new(0.6239137300, 0.5353281424),
            GaussianPrimitive::new(0.1688554000, 0.4446345420),
        ]
    }

    #[test]
    fn test_contracted_shell_creation() {
        let prims = make_test_primitives();
        let shell = ContractedShell::new(AngularMomentum::S, prims.clone(), [0.0, 0.0, 0.74], 0);

        assert_eq!(shell.angular_momentum, AngularMomentum::S);
        assert_eq!(shell.primitives.len(), 3);
        assert_abs_diff_eq!(shell.center[0], 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(shell.center[1], 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(shell.center[2], 0.74, epsilon = 1e-10);
        assert_eq!(shell.atom_idx, 0);
    }

    #[test]
    fn test_shell_n_basis_functions() {
        let prims = vec![GaussianPrimitive::new(1.0, 1.0)];

        let s_shell = ContractedShell::new(AngularMomentum::S, prims.clone(), [0.0; 3], 0);
        let p_shell = ContractedShell::new(AngularMomentum::P, prims.clone(), [0.0; 3], 0);
        let d_shell = ContractedShell::new(AngularMomentum::D, prims, [0.0; 3], 0);

        assert_eq!(s_shell.n_basis_functions(), 1);
        assert_eq!(p_shell.n_basis_functions(), 3);
        assert_eq!(d_shell.n_basis_functions(), 6);
    }

    #[test]
    fn test_shell_n_primitives() {
        let prims = make_test_primitives();
        let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0; 3], 0);

        assert_eq!(shell.n_primitives(), 3);
    }

    #[test]
    fn test_shell_l_value() {
        let prims = vec![GaussianPrimitive::new(1.0, 1.0)];

        let s_shell = ContractedShell::new(AngularMomentum::S, prims.clone(), [0.0; 3], 0);
        let p_shell = ContractedShell::new(AngularMomentum::P, prims.clone(), [0.0; 3], 0);
        let d_shell = ContractedShell::new(AngularMomentum::D, prims, [0.0; 3], 0);

        assert_eq!(s_shell.l_value(), 0);
        assert_eq!(p_shell.l_value(), 1);
        assert_eq!(d_shell.l_value(), 2);
    }

    #[test]
    fn test_shell_type_predicates() {
        let prims = vec![GaussianPrimitive::new(1.0, 1.0)];

        let s_shell = ContractedShell::new(AngularMomentum::S, prims.clone(), [0.0; 3], 0);
        let p_shell = ContractedShell::new(AngularMomentum::P, prims.clone(), [0.0; 3], 0);
        let d_shell = ContractedShell::new(AngularMomentum::D, prims, [0.0; 3], 0);

        assert!(s_shell.is_s_shell());
        assert!(!s_shell.is_p_shell());
        assert!(!s_shell.is_d_shell());

        assert!(!p_shell.is_s_shell());
        assert!(p_shell.is_p_shell());
        assert!(!p_shell.is_d_shell());

        assert!(!d_shell.is_s_shell());
        assert!(!d_shell.is_p_shell());
        assert!(d_shell.is_d_shell());
    }

    #[test]
    fn test_shell_is_valid() {
        // Valid shell
        let valid_prims = vec![
            GaussianPrimitive::new(1.0, 0.5),
            GaussianPrimitive::new(0.3, 0.5),
        ];
        let valid_shell = ContractedShell::new(AngularMomentum::S, valid_prims, [0.0; 3], 0);
        assert!(valid_shell.is_valid());

        // Invalid: negative exponent
        let invalid_prims = vec![
            GaussianPrimitive::new(1.0, 0.5),
            GaussianPrimitive::new(-0.3, 0.5), // Invalid
        ];
        let invalid_shell = ContractedShell::new(AngularMomentum::S, invalid_prims, [0.0; 3], 0);
        assert!(!invalid_shell.is_valid());

        // Invalid: empty primitives
        let empty_shell = ContractedShell::new(AngularMomentum::S, vec![], [0.0; 3], 0);
        assert!(!empty_shell.is_valid());
    }

    #[test]
    fn test_shell_coordinates() {
        let prims = vec![GaussianPrimitive::new(1.0, 1.0)];
        let shell = ContractedShell::new(AngularMomentum::S, prims, [1.0, 2.0, 3.0], 0);

        assert_abs_diff_eq!(shell.x(), 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(shell.y(), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(shell.z(), 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_shell_max_exponent() {
        let prims = make_test_primitives();
        let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0; 3], 0);

        assert_abs_diff_eq!(shell.max_exponent(), 3.4252509100, epsilon = 1e-10);
    }

    #[test]
    fn test_shell_min_exponent() {
        let prims = make_test_primitives();
        let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0; 3], 0);

        assert_abs_diff_eq!(shell.min_exponent(), 0.1688554000, epsilon = 1e-10);
    }

    #[test]
    fn test_shell_display() {
        let prims = vec![GaussianPrimitive::new(1.0, 1.0)];
        let shell = ContractedShell::new(AngularMomentum::P, prims, [1.0, 2.0, 3.0], 5);

        let display = format!("{}", shell);
        assert!(display.contains("P"));
        assert!(display.contains("(1)")); // 1 primitive
        assert!(display.contains("atom 5"));
    }

    #[test]
    fn test_shell_clone() {
        let prims = make_test_primitives();
        let shell1 = ContractedShell::new(AngularMomentum::S, prims, [0.0; 3], 0);
        let shell2 = shell1.clone();

        assert_eq!(shell1, shell2);
        assert_eq!(shell1.n_primitives(), shell2.n_primitives());
    }

    #[test]
    fn test_shell_serde() {
        let prims = make_test_primitives();
        let shell = ContractedShell::new(AngularMomentum::P, prims, [1.0, 2.0, 3.0], 1);

        let serialized = serde_json::to_string(&shell).unwrap();
        let deserialized: ContractedShell = serde_json::from_str(&serialized).unwrap();

        assert_eq!(shell, deserialized);
    }

    // -------------------------------------------------------------------------
    // Reference test: H2O STO-3G shell structure
    // -------------------------------------------------------------------------

    #[test]
    fn test_h2o_sto3g_shell_structure() {
        // H2O STO-3G has 5 shells:
        // O: 1s (S), 2s (S), 2p (P)
        // H1: 1s (S)
        // H2: 1s (S)

        // This tests that the shell types give correct basis function counts
        let o_1s_prims = vec![GaussianPrimitive::new(130.709, 0.15433)];
        let o_2s_prims = vec![GaussianPrimitive::new(5.033, -0.09997)];
        let o_2p_prims = vec![GaussianPrimitive::new(5.033, 0.15592)];
        let h_1s_prims = vec![GaussianPrimitive::new(3.425, 0.15433)];

        let o_1s = ContractedShell::new(AngularMomentum::S, o_1s_prims, [0.0, 0.0, 0.0], 0);
        let o_2s = ContractedShell::new(AngularMomentum::S, o_2s_prims, [0.0, 0.0, 0.0], 0);
        let o_2p = ContractedShell::new(AngularMomentum::P, o_2p_prims, [0.0, 0.0, 0.0], 0);
        let h1_1s =
            ContractedShell::new(AngularMomentum::S, h_1s_prims.clone(), [1.0, 0.0, 0.0], 1);
        let h2_1s = ContractedShell::new(AngularMomentum::S, h_1s_prims, [-1.0, 0.0, 0.0], 2);

        // Total basis functions: 1 + 1 + 3 + 1 + 1 = 7
        let total_basis: usize = [&o_1s, &o_2s, &o_2p, &h1_1s, &h2_1s]
            .iter()
            .map(|s| s.n_basis_functions())
            .sum();

        assert_eq!(total_basis, 7);
    }
}
