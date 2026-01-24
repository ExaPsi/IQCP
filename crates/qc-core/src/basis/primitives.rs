//! Primitive Gaussian functions and angular momentum types
//!
//! This module defines the fundamental building blocks for Gaussian-type
//! basis functions used in quantum chemistry calculations.
//!
//! # Background
//!
//! Gaussian-type orbitals (GTOs) are used as basis functions in most modern
//! quantum chemistry programs due to their computational efficiency in integral
//! evaluation. A primitive Gaussian is:
//!
//! ```text
//! g(r) = N * x^l * y^m * z^n * exp(-alpha * r^2)
//! ```
//!
//! where `N` is a normalization constant, `alpha` is the exponent,
//! and `(l, m, n)` define the angular momentum.
//!
//! # References
//!
//! - Szabo & Ostlund (1996). "Modern Quantum Chemistry". Dover, Chapter 3.
//! - PySCF basis set implementation: `references/pyscf/pyscf/gto/basis/`

use serde::{Deserialize, Serialize};
use std::fmt;

// =============================================================================
// Angular Momentum
// =============================================================================

/// Angular momentum quantum number for basis functions
///
/// In Cartesian Gaussian functions, angular momentum determines the
/// polynomial prefactor and the number of components:
///
/// - S (l=0): 1 component (x^0 y^0 z^0)
/// - P (l=1): 3 components (x, y, z)
/// - D (l=2): 6 components (xx, xy, xz, yy, yz, zz)
///
/// # Formula
///
/// Number of Cartesian components: `(l+1)(l+2)/2`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AngularMomentum {
    /// S shell (l=0): 1 Cartesian component
    S = 0,
    /// P shell (l=1): 3 Cartesian components (px, py, pz)
    P = 1,
    /// D shell (l=2): 6 Cartesian components (dxx, dxy, dxz, dyy, dyz, dzz)
    D = 2,
}

impl AngularMomentum {
    /// Number of Cartesian basis functions for this angular momentum
    ///
    /// Formula: `(l+1)(l+2)/2`
    /// - S (l=0): 1
    /// - P (l=1): 3
    /// - D (l=2): 6
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::AngularMomentum;
    ///
    /// assert_eq!(AngularMomentum::S.n_cartesian(), 1);
    /// assert_eq!(AngularMomentum::P.n_cartesian(), 3);
    /// assert_eq!(AngularMomentum::D.n_cartesian(), 6);
    /// ```
    #[inline]
    pub fn n_cartesian(&self) -> usize {
        let l = *self as usize;
        (l + 1) * (l + 2) / 2
    }

    /// Get the angular momentum quantum number (l value)
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::AngularMomentum;
    ///
    /// assert_eq!(AngularMomentum::S.l_value(), 0);
    /// assert_eq!(AngularMomentum::P.l_value(), 1);
    /// assert_eq!(AngularMomentum::D.l_value(), 2);
    /// ```
    #[inline]
    pub fn l_value(&self) -> u32 {
        *self as u32
    }

    /// Convert from l value to AngularMomentum
    ///
    /// # Returns
    ///
    /// - `Some(AngularMomentum)` if l is 0, 1, or 2
    /// - `None` for unsupported l values
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::AngularMomentum;
    ///
    /// assert_eq!(AngularMomentum::from_l(0), Some(AngularMomentum::S));
    /// assert_eq!(AngularMomentum::from_l(1), Some(AngularMomentum::P));
    /// assert_eq!(AngularMomentum::from_l(2), Some(AngularMomentum::D));
    /// assert_eq!(AngularMomentum::from_l(3), None);
    /// ```
    #[inline]
    pub fn from_l(l: u32) -> Option<Self> {
        match l {
            0 => Some(Self::S),
            1 => Some(Self::P),
            2 => Some(Self::D),
            _ => None,
        }
    }

    /// Convert from shell symbol to AngularMomentum
    ///
    /// Accepts both uppercase and lowercase symbols.
    ///
    /// # Returns
    ///
    /// - `Some(AngularMomentum)` for "S", "s", "P", "p", "D", "d"
    /// - `None` for other symbols
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::AngularMomentum;
    ///
    /// assert_eq!(AngularMomentum::from_symbol("S"), Some(AngularMomentum::S));
    /// assert_eq!(AngularMomentum::from_symbol("p"), Some(AngularMomentum::P));
    /// assert_eq!(AngularMomentum::from_symbol("D"), Some(AngularMomentum::D));
    /// assert_eq!(AngularMomentum::from_symbol("F"), None);
    /// ```
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        match symbol.to_uppercase().as_str() {
            "S" => Some(Self::S),
            "P" => Some(Self::P),
            "D" => Some(Self::D),
            _ => None,
        }
    }

    /// Get the shell symbol
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::AngularMomentum;
    ///
    /// assert_eq!(AngularMomentum::S.symbol(), "S");
    /// assert_eq!(AngularMomentum::P.symbol(), "P");
    /// assert_eq!(AngularMomentum::D.symbol(), "D");
    /// ```
    #[inline]
    pub fn symbol(&self) -> &'static str {
        match self {
            Self::S => "S",
            Self::P => "P",
            Self::D => "D",
        }
    }

    /// Number of spherical harmonic basis functions for this angular momentum
    ///
    /// Formula: `2l + 1`
    /// - S (l=0): 1
    /// - P (l=1): 3
    /// - D (l=2): 5
    ///
    /// Note: For S and P orbitals, spherical equals Cartesian.
    /// The difference appears at D (5 vs 6) and higher.
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::AngularMomentum;
    ///
    /// assert_eq!(AngularMomentum::S.n_spherical(), 1);
    /// assert_eq!(AngularMomentum::P.n_spherical(), 3);
    /// assert_eq!(AngularMomentum::D.n_spherical(), 5);
    /// ```
    #[inline]
    pub fn n_spherical(&self) -> usize {
        let l = *self as usize;
        2 * l + 1
    }
}

impl fmt::Display for AngularMomentum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

// =============================================================================
// Gaussian Primitive
// =============================================================================

/// A primitive Gaussian function
///
/// Represents a single Gaussian function:
/// ```text
/// g(r) = N * exp(-alpha * |r - R|^2)
/// ```
///
/// where `alpha` is the exponent and `N` is a normalization constant
/// (combined with the coefficient in contracted basis sets).
///
/// # Fields
///
/// - `exponent`: The Gaussian exponent (alpha) in atomic units (Bohr^-2)
/// - `coefficient`: The contraction coefficient (pre-normalized from basis set library)
///
/// # Note on Normalization
///
/// The coefficients stored here are "pre-normalized" as they appear in
/// standard basis set libraries like PySCF. This means the normalization
/// factor for each primitive is already incorporated into the coefficient.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GaussianPrimitive {
    /// Exponent (alpha) in atomic units (Bohr^-2)
    ///
    /// Determines how "tight" or "diffuse" the Gaussian is:
    /// - Large exponent: tight, localized function (core electrons)
    /// - Small exponent: diffuse, spread-out function (valence electrons)
    pub exponent: f64,

    /// Contraction coefficient (pre-normalized)
    ///
    /// This coefficient is multiplied by the primitive Gaussian when
    /// forming a contracted shell. The value includes the normalization
    /// factor from the original basis set definition.
    pub coefficient: f64,
}

impl GaussianPrimitive {
    /// Create a new Gaussian primitive
    ///
    /// # Arguments
    ///
    /// * `exponent` - Gaussian exponent (alpha) in Bohr^-2
    /// * `coefficient` - Contraction coefficient (pre-normalized)
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::GaussianPrimitive;
    ///
    /// // First primitive of hydrogen STO-3G basis
    /// let prim = GaussianPrimitive::new(3.4252509100, 0.1543289707);
    /// assert!((prim.exponent - 3.4252509100).abs() < 1e-10);
    /// assert!((prim.coefficient - 0.1543289707).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn new(exponent: f64, coefficient: f64) -> Self {
        Self {
            exponent,
            coefficient,
        }
    }

    /// Check if the primitive is valid (positive exponent)
    ///
    /// # Returns
    ///
    /// `true` if the exponent is positive, `false` otherwise
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::GaussianPrimitive;
    ///
    /// let valid = GaussianPrimitive::new(1.0, 0.5);
    /// assert!(valid.is_valid());
    ///
    /// let invalid = GaussianPrimitive::new(-1.0, 0.5);
    /// assert!(!invalid.is_valid());
    /// ```
    #[inline]
    pub fn is_valid(&self) -> bool {
        self.exponent > 0.0
    }
}

impl fmt::Display for GaussianPrimitive {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "GaussianPrimitive(exp={:.6e}, coef={:.10})",
            self.exponent, self.coefficient
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

    // -------------------------------------------------------------------------
    // AngularMomentum tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_angular_momentum_n_cartesian() {
        // S shell: (0+1)(0+2)/2 = 1
        assert_eq!(AngularMomentum::S.n_cartesian(), 1);
        // P shell: (1+1)(1+2)/2 = 3
        assert_eq!(AngularMomentum::P.n_cartesian(), 3);
        // D shell: (2+1)(2+2)/2 = 6
        assert_eq!(AngularMomentum::D.n_cartesian(), 6);
    }

    #[test]
    fn test_angular_momentum_n_spherical() {
        // S shell: 2*0+1 = 1
        assert_eq!(AngularMomentum::S.n_spherical(), 1);
        // P shell: 2*1+1 = 3
        assert_eq!(AngularMomentum::P.n_spherical(), 3);
        // D shell: 2*2+1 = 5
        assert_eq!(AngularMomentum::D.n_spherical(), 5);
    }

    #[test]
    fn test_angular_momentum_spherical_vs_cartesian() {
        // S and P are the same
        assert_eq!(
            AngularMomentum::S.n_cartesian(),
            AngularMomentum::S.n_spherical()
        );
        assert_eq!(
            AngularMomentum::P.n_cartesian(),
            AngularMomentum::P.n_spherical()
        );
        // D differs: 6 Cartesian vs 5 spherical
        assert_ne!(
            AngularMomentum::D.n_cartesian(),
            AngularMomentum::D.n_spherical()
        );
        assert_eq!(AngularMomentum::D.n_cartesian(), 6);
        assert_eq!(AngularMomentum::D.n_spherical(), 5);
    }

    #[test]
    fn test_angular_momentum_l_value() {
        assert_eq!(AngularMomentum::S.l_value(), 0);
        assert_eq!(AngularMomentum::P.l_value(), 1);
        assert_eq!(AngularMomentum::D.l_value(), 2);
    }

    #[test]
    fn test_angular_momentum_from_l() {
        assert_eq!(AngularMomentum::from_l(0), Some(AngularMomentum::S));
        assert_eq!(AngularMomentum::from_l(1), Some(AngularMomentum::P));
        assert_eq!(AngularMomentum::from_l(2), Some(AngularMomentum::D));
        assert_eq!(AngularMomentum::from_l(3), None);
        assert_eq!(AngularMomentum::from_l(100), None);
    }

    #[test]
    fn test_angular_momentum_from_symbol() {
        // Uppercase
        assert_eq!(AngularMomentum::from_symbol("S"), Some(AngularMomentum::S));
        assert_eq!(AngularMomentum::from_symbol("P"), Some(AngularMomentum::P));
        assert_eq!(AngularMomentum::from_symbol("D"), Some(AngularMomentum::D));

        // Lowercase
        assert_eq!(AngularMomentum::from_symbol("s"), Some(AngularMomentum::S));
        assert_eq!(AngularMomentum::from_symbol("p"), Some(AngularMomentum::P));
        assert_eq!(AngularMomentum::from_symbol("d"), Some(AngularMomentum::D));

        // Invalid
        assert_eq!(AngularMomentum::from_symbol("F"), None);
        assert_eq!(AngularMomentum::from_symbol("G"), None);
        assert_eq!(AngularMomentum::from_symbol(""), None);
        assert_eq!(AngularMomentum::from_symbol("SP"), None);
    }

    #[test]
    fn test_angular_momentum_symbol() {
        assert_eq!(AngularMomentum::S.symbol(), "S");
        assert_eq!(AngularMomentum::P.symbol(), "P");
        assert_eq!(AngularMomentum::D.symbol(), "D");
    }

    #[test]
    fn test_angular_momentum_display() {
        assert_eq!(format!("{}", AngularMomentum::S), "S");
        assert_eq!(format!("{}", AngularMomentum::P), "P");
        assert_eq!(format!("{}", AngularMomentum::D), "D");
    }

    #[test]
    fn test_angular_momentum_equality_and_hash() {
        use std::collections::HashSet;

        assert_eq!(AngularMomentum::S, AngularMomentum::S);
        assert_ne!(AngularMomentum::S, AngularMomentum::P);

        // Test hashability
        let mut set = HashSet::new();
        set.insert(AngularMomentum::S);
        set.insert(AngularMomentum::P);
        set.insert(AngularMomentum::S); // duplicate

        assert_eq!(set.len(), 2);
        assert!(set.contains(&AngularMomentum::S));
        assert!(set.contains(&AngularMomentum::P));
    }

    #[test]
    fn test_angular_momentum_serde() {
        let am = AngularMomentum::P;
        let serialized = serde_json::to_string(&am).unwrap();
        let deserialized: AngularMomentum = serde_json::from_str(&serialized).unwrap();
        assert_eq!(am, deserialized);
    }

    // -------------------------------------------------------------------------
    // GaussianPrimitive tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_gaussian_primitive_creation() {
        let prim = GaussianPrimitive::new(3.4252509100, 0.1543289707);

        assert_abs_diff_eq!(prim.exponent, 3.4252509100, epsilon = 1e-10);
        assert_abs_diff_eq!(prim.coefficient, 0.1543289707, epsilon = 1e-10);
    }

    #[test]
    fn test_gaussian_primitive_is_valid() {
        // Valid primitive
        let valid = GaussianPrimitive::new(1.0, 0.5);
        assert!(valid.is_valid());

        // Valid with large exponent
        let large = GaussianPrimitive::new(1000.0, 0.001);
        assert!(large.is_valid());

        // Valid with small exponent
        let small = GaussianPrimitive::new(0.01, 0.9);
        assert!(small.is_valid());

        // Invalid: negative exponent
        let negative = GaussianPrimitive::new(-1.0, 0.5);
        assert!(!negative.is_valid());

        // Invalid: zero exponent
        let zero = GaussianPrimitive::new(0.0, 0.5);
        assert!(!zero.is_valid());

        // Valid: negative coefficient is allowed (can occur in contracted shells)
        let neg_coef = GaussianPrimitive::new(1.0, -0.1);
        assert!(neg_coef.is_valid());
    }

    #[test]
    fn test_gaussian_primitive_display() {
        let prim = GaussianPrimitive::new(3.425251, 0.1543289707);
        let display = format!("{}", prim);
        assert!(display.contains("GaussianPrimitive"));
        assert!(display.contains("3.425251"));
    }

    #[test]
    fn test_gaussian_primitive_equality() {
        let p1 = GaussianPrimitive::new(1.0, 0.5);
        let p2 = GaussianPrimitive::new(1.0, 0.5);
        let p3 = GaussianPrimitive::new(1.0, 0.6);
        let p4 = GaussianPrimitive::new(2.0, 0.5);

        assert_eq!(p1, p2);
        assert_ne!(p1, p3);
        assert_ne!(p1, p4);
    }

    #[test]
    fn test_gaussian_primitive_copy() {
        let p1 = GaussianPrimitive::new(1.0, 0.5);
        let p2 = p1; // Copy
        assert_eq!(p1.exponent, p2.exponent);
        assert_eq!(p1.coefficient, p2.coefficient);
    }

    #[test]
    fn test_gaussian_primitive_serde() {
        let prim = GaussianPrimitive::new(3.4252509100, 0.1543289707);
        let serialized = serde_json::to_string(&prim).unwrap();
        let deserialized: GaussianPrimitive = serde_json::from_str(&serialized).unwrap();

        assert_abs_diff_eq!(prim.exponent, deserialized.exponent, epsilon = 1e-15);
        assert_abs_diff_eq!(prim.coefficient, deserialized.coefficient, epsilon = 1e-15);
    }

    // -------------------------------------------------------------------------
    // Reference values from STO-3G (PySCF 2.11.0)
    // -------------------------------------------------------------------------

    #[test]
    fn test_hydrogen_sto3g_primitives() {
        // Hydrogen 1s shell in STO-3G from PySCF
        let h_primitives = [
            GaussianPrimitive::new(3.4252509100, 0.1543289707),
            GaussianPrimitive::new(0.6239137300, 0.5353281424),
            GaussianPrimitive::new(0.1688554000, 0.4446345420),
        ];

        // Check exponents decrease (tighter to diffuse)
        assert!(h_primitives[0].exponent > h_primitives[1].exponent);
        assert!(h_primitives[1].exponent > h_primitives[2].exponent);

        // All coefficients should be positive for H 1s
        for prim in &h_primitives {
            assert!(prim.coefficient > 0.0);
            assert!(prim.is_valid());
        }
    }

    #[test]
    fn test_carbon_sto3g_valence_primitives() {
        // Carbon 2s shell has a negative coefficient (inner/outer mixing)
        let c_2s_primitives = [
            GaussianPrimitive::new(2.9412494000, -0.0999672301),
            GaussianPrimitive::new(0.6834831000, 0.3995128303),
            GaussianPrimitive::new(0.2222899000, 0.7001154705),
        ];

        // First primitive has negative coefficient (valid)
        assert!(c_2s_primitives[0].coefficient < 0.0);
        assert!(c_2s_primitives[0].is_valid());

        // All exponents should be positive
        for prim in &c_2s_primitives {
            assert!(prim.exponent > 0.0);
        }
    }
}
