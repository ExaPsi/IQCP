//! Atom data structures and element utilities
//!
//! This module provides atomic data types and element symbol/atomic number
//! conversion utilities for the IQCP basis set infrastructure.
//!
//! # Supported Elements
//!
//! For the IQCP paper scope (educational quantum chemistry), we support
//! elements H through Ne (Z = 1-10). This covers the first two periods
//! and includes all common small molecule examples:
//!
//! | Z | Symbol | Name |
//! |---|--------|------|
//! | 1 | H | Hydrogen |
//! | 2 | He | Helium |
//! | 3 | Li | Lithium |
//! | 4 | Be | Beryllium |
//! | 5 | B | Boron |
//! | 6 | C | Carbon |
//! | 7 | N | Nitrogen |
//! | 8 | O | Oxygen |
//! | 9 | F | Fluorine |
//! | 10 | Ne | Neon |
//!
//! # Coordinate System
//!
//! All atomic coordinates are in Bohr (atomic units). The conversion from
//! Angstrom is: 1 Bohr = 0.529177210903 Angstrom.
//!
//! # References
//!
//! - NIST atomic data: https://physics.nist.gov/PhysRefData/Elements/

use serde::{Deserialize, Serialize};
use std::fmt;
use thiserror::Error;

// =============================================================================
// Constants
// =============================================================================

/// Conversion factor from Bohr to Angstrom (CODATA 2014, matches PySCF)
/// Reference: PySCF lib/param.py BOHR = 0.52917721092
pub const BOHR_TO_ANGSTROM: f64 = 0.52917721092;

/// Conversion factor from Angstrom to Bohr
/// This is 1 / BOHR_TO_ANGSTROM
pub const ANGSTROM_TO_BOHR: f64 = 1.8897261245650618;

/// Minimum supported atomic number (Hydrogen)
pub const MIN_ATOMIC_NUMBER: u8 = 1;

/// Maximum supported atomic number (Neon)
pub const MAX_ATOMIC_NUMBER: u8 = 10;

/// Element symbols for Z = 1 to 10
const ELEMENT_SYMBOLS: [&str; 10] = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"];

/// Element names for Z = 1 to 10
const ELEMENT_NAMES: [&str; 10] = [
    "Hydrogen",
    "Helium",
    "Lithium",
    "Beryllium",
    "Boron",
    "Carbon",
    "Nitrogen",
    "Oxygen",
    "Fluorine",
    "Neon",
];

// =============================================================================
// Error Types
// =============================================================================

/// Errors that can occur when working with elements
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum ElementError {
    /// Unknown element symbol
    #[error("Unknown element symbol: '{0}'")]
    UnknownSymbol(String),

    /// Unsupported atomic number (outside Z=1-10 range)
    #[error("Unsupported atomic number: {0} (must be 1-10 for IQCP)")]
    UnsupportedElement(u8),
}

// =============================================================================
// Element Conversion Functions
// =============================================================================

/// Convert element symbol to atomic number
///
/// Accepts both uppercase and lowercase symbols (e.g., "H", "h", "HE", "He").
///
/// # Arguments
///
/// * `symbol` - Element symbol (case insensitive)
///
/// # Returns
///
/// - `Ok(u8)` - Atomic number (1-10)
/// - `Err(ElementError::UnknownSymbol)` - Symbol not recognized
///
/// # Example
///
/// ```rust
/// use qc_core::basis::symbol_to_atomic_number;
///
/// assert_eq!(symbol_to_atomic_number("H").unwrap(), 1);
/// assert_eq!(symbol_to_atomic_number("He").unwrap(), 2);
/// assert_eq!(symbol_to_atomic_number("O").unwrap(), 8);
/// assert_eq!(symbol_to_atomic_number("ne").unwrap(), 10);
///
/// assert!(symbol_to_atomic_number("X").is_err());
/// ```
pub fn symbol_to_atomic_number(symbol: &str) -> Result<u8, ElementError> {
    let upper = symbol.trim().to_uppercase();

    match upper.as_str() {
        "H" => Ok(1),
        "HE" => Ok(2),
        "LI" => Ok(3),
        "BE" => Ok(4),
        "B" => Ok(5),
        "C" => Ok(6),
        "N" => Ok(7),
        "O" => Ok(8),
        "F" => Ok(9),
        "NE" => Ok(10),
        _ => Err(ElementError::UnknownSymbol(symbol.to_string())),
    }
}

/// Convert atomic number to element symbol
///
/// # Arguments
///
/// * `z` - Atomic number (1-10)
///
/// # Returns
///
/// - `Ok(&str)` - Element symbol
/// - `Err(ElementError::UnsupportedElement)` - Atomic number outside supported range
///
/// # Example
///
/// ```rust
/// use qc_core::basis::atomic_number_to_symbol;
///
/// assert_eq!(atomic_number_to_symbol(1).unwrap(), "H");
/// assert_eq!(atomic_number_to_symbol(6).unwrap(), "C");
/// assert_eq!(atomic_number_to_symbol(10).unwrap(), "Ne");
///
/// assert!(atomic_number_to_symbol(0).is_err());
/// assert!(atomic_number_to_symbol(11).is_err());
/// ```
pub fn atomic_number_to_symbol(z: u8) -> Result<&'static str, ElementError> {
    if !(MIN_ATOMIC_NUMBER..=MAX_ATOMIC_NUMBER).contains(&z) {
        return Err(ElementError::UnsupportedElement(z));
    }
    Ok(ELEMENT_SYMBOLS[(z - 1) as usize])
}

/// Convert atomic number to element name
///
/// # Arguments
///
/// * `z` - Atomic number (1-10)
///
/// # Returns
///
/// - `Ok(&str)` - Element name (e.g., "Hydrogen")
/// - `Err(ElementError::UnsupportedElement)` - Atomic number outside supported range
///
/// # Example
///
/// ```rust
/// use qc_core::basis::atomic_number_to_name;
///
/// assert_eq!(atomic_number_to_name(1).unwrap(), "Hydrogen");
/// assert_eq!(atomic_number_to_name(6).unwrap(), "Carbon");
/// assert_eq!(atomic_number_to_name(8).unwrap(), "Oxygen");
/// ```
pub fn atomic_number_to_name(z: u8) -> Result<&'static str, ElementError> {
    if !(MIN_ATOMIC_NUMBER..=MAX_ATOMIC_NUMBER).contains(&z) {
        return Err(ElementError::UnsupportedElement(z));
    }
    Ok(ELEMENT_NAMES[(z - 1) as usize])
}

/// Check if an atomic number is supported
///
/// # Example
///
/// ```rust
/// use qc_core::basis::is_supported_element;
///
/// assert!(is_supported_element(1));   // H
/// assert!(is_supported_element(10));  // Ne
/// assert!(!is_supported_element(0));  // Invalid
/// assert!(!is_supported_element(11)); // Na (not supported)
/// ```
#[inline]
pub fn is_supported_element(z: u8) -> bool {
    (MIN_ATOMIC_NUMBER..=MAX_ATOMIC_NUMBER).contains(&z)
}

// =============================================================================
// Atom
// =============================================================================

/// An atom in a molecular system
///
/// Represents an atomic nucleus with its element type and position in space.
/// Coordinates are always stored in atomic units (Bohr).
///
/// # Example
///
/// ```rust
/// use qc_core::basis::Atom;
///
/// // Create a hydrogen atom at the origin
/// let h = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
/// assert_eq!(h.symbol(), "H");
/// assert_eq!(h.atomic_number, 1);
///
/// // Create an oxygen atom
/// let o = Atom::from_symbol("O", [0.0, 0.0, 0.0]).unwrap();
/// assert_eq!(o.atomic_number, 8);
/// ```
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    /// Atomic number (Z = 1-10)
    pub atomic_number: u8,

    /// Element symbol (e.g., "H", "C", "O")
    pub symbol: String,

    /// Cartesian coordinates in Bohr [x, y, z]
    pub position: [f64; 3],
}

impl Atom {
    /// Create a new atom from atomic number
    ///
    /// # Arguments
    ///
    /// * `atomic_number` - Atomic number (1-10)
    /// * `position` - Cartesian coordinates in Bohr [x, y, z]
    ///
    /// # Returns
    ///
    /// - `Ok(Atom)` - Successfully created atom
    /// - `Err(ElementError::UnsupportedElement)` - Invalid atomic number
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::Atom;
    ///
    /// let carbon = Atom::new(6, [0.0, 0.0, 0.0]).unwrap();
    /// assert_eq!(carbon.symbol(), "C");
    /// assert_eq!(carbon.atomic_number, 6);
    /// ```
    pub fn new(atomic_number: u8, position: [f64; 3]) -> Result<Self, ElementError> {
        let symbol = atomic_number_to_symbol(atomic_number)?.to_string();
        Ok(Self {
            atomic_number,
            symbol,
            position,
        })
    }

    /// Create a new atom from element symbol
    ///
    /// # Arguments
    ///
    /// * `symbol` - Element symbol (case insensitive)
    /// * `position` - Cartesian coordinates in Bohr [x, y, z]
    ///
    /// # Returns
    ///
    /// - `Ok(Atom)` - Successfully created atom
    /// - `Err(ElementError::UnknownSymbol)` - Unknown element symbol
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::Atom;
    ///
    /// let oxygen = Atom::from_symbol("O", [0.0, 0.0, 0.117]).unwrap();
    /// assert_eq!(oxygen.atomic_number, 8);
    ///
    /// // Case insensitive
    /// let nitrogen = Atom::from_symbol("n", [0.0, 0.0, 0.0]).unwrap();
    /// assert_eq!(nitrogen.atomic_number, 7);
    /// ```
    pub fn from_symbol(symbol: &str, position: [f64; 3]) -> Result<Self, ElementError> {
        let atomic_number = symbol_to_atomic_number(symbol)?;
        let symbol_str = atomic_number_to_symbol(atomic_number)?.to_string();
        Ok(Self {
            atomic_number,
            symbol: symbol_str,
            position,
        })
    }

    /// Create an atom with position in Angstrom (will be converted to Bohr)
    ///
    /// # Arguments
    ///
    /// * `atomic_number` - Atomic number (1-10)
    /// * `position_angstrom` - Cartesian coordinates in Angstrom [x, y, z]
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::Atom;
    ///
    /// // H2 with bond length 0.74 Angstrom
    /// let h1 = Atom::new_angstrom(1, [0.0, 0.0, 0.0]).unwrap();
    /// let h2 = Atom::new_angstrom(1, [0.0, 0.0, 0.74]).unwrap();
    ///
    /// // Positions are stored in Bohr
    /// assert!(h2.position[2] > 1.0); // 0.74 Angstrom > 1 Bohr
    /// ```
    pub fn new_angstrom(
        atomic_number: u8,
        position_angstrom: [f64; 3],
    ) -> Result<Self, ElementError> {
        let position_bohr = [
            position_angstrom[0] * ANGSTROM_TO_BOHR,
            position_angstrom[1] * ANGSTROM_TO_BOHR,
            position_angstrom[2] * ANGSTROM_TO_BOHR,
        ];
        Self::new(atomic_number, position_bohr)
    }

    /// Get the element symbol
    #[inline]
    pub fn symbol(&self) -> &str {
        &self.symbol
    }

    /// Get the element name (e.g., "Hydrogen")
    pub fn name(&self) -> &'static str {
        // Safe because we validate atomic_number on construction
        atomic_number_to_name(self.atomic_number).unwrap_or("Unknown")
    }

    /// Get the nuclear charge (same as atomic number for neutral atoms)
    #[inline]
    pub fn nuclear_charge(&self) -> f64 {
        self.atomic_number as f64
    }

    /// Get the x coordinate in Bohr
    #[inline]
    pub fn x(&self) -> f64 {
        self.position[0]
    }

    /// Get the y coordinate in Bohr
    #[inline]
    pub fn y(&self) -> f64 {
        self.position[1]
    }

    /// Get the z coordinate in Bohr
    #[inline]
    pub fn z(&self) -> f64 {
        self.position[2]
    }

    /// Get position as a slice
    #[inline]
    pub fn position_slice(&self) -> &[f64; 3] {
        &self.position
    }

    /// Calculate distance to another atom in Bohr
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::Atom;
    ///
    /// let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
    /// let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
    ///
    /// let dist = h1.distance_to(&h2);
    /// assert!((dist - 1.4).abs() < 1e-10);
    /// ```
    pub fn distance_to(&self, other: &Atom) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} [{:.6}, {:.6}, {:.6}]",
            self.symbol, self.position[0], self.position[1], self.position[2]
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
    // Element conversion tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_symbol_to_atomic_number() {
        // Standard cases
        assert_eq!(symbol_to_atomic_number("H").unwrap(), 1);
        assert_eq!(symbol_to_atomic_number("He").unwrap(), 2);
        assert_eq!(symbol_to_atomic_number("Li").unwrap(), 3);
        assert_eq!(symbol_to_atomic_number("Be").unwrap(), 4);
        assert_eq!(symbol_to_atomic_number("B").unwrap(), 5);
        assert_eq!(symbol_to_atomic_number("C").unwrap(), 6);
        assert_eq!(symbol_to_atomic_number("N").unwrap(), 7);
        assert_eq!(symbol_to_atomic_number("O").unwrap(), 8);
        assert_eq!(symbol_to_atomic_number("F").unwrap(), 9);
        assert_eq!(symbol_to_atomic_number("Ne").unwrap(), 10);
    }

    #[test]
    fn test_symbol_to_atomic_number_case_insensitive() {
        assert_eq!(symbol_to_atomic_number("h").unwrap(), 1);
        assert_eq!(symbol_to_atomic_number("HE").unwrap(), 2);
        assert_eq!(symbol_to_atomic_number("he").unwrap(), 2);
        assert_eq!(symbol_to_atomic_number("C").unwrap(), 6);
        assert_eq!(symbol_to_atomic_number("c").unwrap(), 6);
        assert_eq!(symbol_to_atomic_number("NE").unwrap(), 10);
        assert_eq!(symbol_to_atomic_number("ne").unwrap(), 10);
    }

    #[test]
    fn test_symbol_to_atomic_number_with_whitespace() {
        assert_eq!(symbol_to_atomic_number(" H ").unwrap(), 1);
        assert_eq!(symbol_to_atomic_number("  O  ").unwrap(), 8);
    }

    #[test]
    fn test_symbol_to_atomic_number_unknown() {
        let result = symbol_to_atomic_number("X");
        assert!(matches!(result, Err(ElementError::UnknownSymbol(_))));

        let result = symbol_to_atomic_number("Na"); // Sodium not supported
        assert!(matches!(result, Err(ElementError::UnknownSymbol(_))));

        let result = symbol_to_atomic_number("");
        assert!(matches!(result, Err(ElementError::UnknownSymbol(_))));
    }

    #[test]
    fn test_atomic_number_to_symbol() {
        assert_eq!(atomic_number_to_symbol(1).unwrap(), "H");
        assert_eq!(atomic_number_to_symbol(2).unwrap(), "He");
        assert_eq!(atomic_number_to_symbol(3).unwrap(), "Li");
        assert_eq!(atomic_number_to_symbol(4).unwrap(), "Be");
        assert_eq!(atomic_number_to_symbol(5).unwrap(), "B");
        assert_eq!(atomic_number_to_symbol(6).unwrap(), "C");
        assert_eq!(atomic_number_to_symbol(7).unwrap(), "N");
        assert_eq!(atomic_number_to_symbol(8).unwrap(), "O");
        assert_eq!(atomic_number_to_symbol(9).unwrap(), "F");
        assert_eq!(atomic_number_to_symbol(10).unwrap(), "Ne");
    }

    #[test]
    fn test_atomic_number_to_symbol_invalid() {
        assert!(matches!(
            atomic_number_to_symbol(0),
            Err(ElementError::UnsupportedElement(0))
        ));
        assert!(matches!(
            atomic_number_to_symbol(11),
            Err(ElementError::UnsupportedElement(11))
        ));
        assert!(matches!(
            atomic_number_to_symbol(100),
            Err(ElementError::UnsupportedElement(100))
        ));
    }

    #[test]
    fn test_atomic_number_to_name() {
        assert_eq!(atomic_number_to_name(1).unwrap(), "Hydrogen");
        assert_eq!(atomic_number_to_name(6).unwrap(), "Carbon");
        assert_eq!(atomic_number_to_name(8).unwrap(), "Oxygen");
        assert_eq!(atomic_number_to_name(10).unwrap(), "Neon");
    }

    #[test]
    fn test_is_supported_element() {
        assert!(!is_supported_element(0));
        assert!(is_supported_element(1));
        assert!(is_supported_element(5));
        assert!(is_supported_element(10));
        assert!(!is_supported_element(11));
        assert!(!is_supported_element(100));
    }

    // -------------------------------------------------------------------------
    // Atom tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_atom_creation() {
        let atom = Atom::new(6, [1.0, 2.0, 3.0]).unwrap();

        assert_eq!(atom.atomic_number, 6);
        assert_eq!(atom.symbol, "C");
        assert_abs_diff_eq!(atom.position[0], 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(atom.position[1], 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(atom.position[2], 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_atom_from_symbol() {
        let atom = Atom::from_symbol("O", [0.0, 0.0, 0.117]).unwrap();

        assert_eq!(atom.atomic_number, 8);
        assert_eq!(atom.symbol, "O");
        assert_abs_diff_eq!(atom.position[2], 0.117, epsilon = 1e-10);
    }

    #[test]
    fn test_atom_from_symbol_case_insensitive() {
        let atom = Atom::from_symbol("c", [0.0, 0.0, 0.0]).unwrap();
        assert_eq!(atom.atomic_number, 6);
        assert_eq!(atom.symbol, "C"); // Normalized to uppercase
    }

    #[test]
    fn test_atom_invalid() {
        assert!(Atom::new(0, [0.0, 0.0, 0.0]).is_err());
        assert!(Atom::new(11, [0.0, 0.0, 0.0]).is_err());
        assert!(Atom::from_symbol("Na", [0.0, 0.0, 0.0]).is_err());
        assert!(Atom::from_symbol("", [0.0, 0.0, 0.0]).is_err());
    }

    #[test]
    fn test_atom_new_angstrom() {
        let atom = Atom::new_angstrom(1, [0.0, 0.0, 0.74]).unwrap();

        // 0.74 Angstrom = 0.74 * 1.8897259886 = 1.3984 Bohr
        assert_abs_diff_eq!(atom.position[2], 0.74 * ANGSTROM_TO_BOHR, epsilon = 1e-10);
    }

    #[test]
    fn test_atom_getters() {
        let atom = Atom::new(8, [1.0, 2.0, 3.0]).unwrap();

        assert_eq!(atom.symbol(), "O");
        assert_eq!(atom.name(), "Oxygen");
        assert_abs_diff_eq!(atom.nuclear_charge(), 8.0, epsilon = 1e-10);
        assert_abs_diff_eq!(atom.x(), 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(atom.y(), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(atom.z(), 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_atom_distance() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();

        assert_abs_diff_eq!(h1.distance_to(&h2), 1.4, epsilon = 1e-10);
        assert_abs_diff_eq!(h2.distance_to(&h1), 1.4, epsilon = 1e-10);
    }

    #[test]
    fn test_atom_distance_3d() {
        let a1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let a2 = Atom::new(1, [1.0, 2.0, 2.0]).unwrap();

        // Distance = sqrt(1 + 4 + 4) = 3
        assert_abs_diff_eq!(a1.distance_to(&a2), 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_atom_display() {
        let atom = Atom::new(6, [1.0, 2.0, 3.0]).unwrap();
        let display = format!("{}", atom);

        assert!(display.contains("C"));
        assert!(display.contains("1.0"));
        assert!(display.contains("2.0"));
        assert!(display.contains("3.0"));
    }

    #[test]
    fn test_atom_clone_and_equality() {
        let a1 = Atom::new(6, [1.0, 2.0, 3.0]).unwrap();
        let a2 = a1.clone();

        assert_eq!(a1, a2);

        let a3 = Atom::new(6, [1.0, 2.0, 3.1]).unwrap();
        assert_ne!(a1, a3);
    }

    #[test]
    fn test_atom_serde() {
        let atom = Atom::new(8, [1.0, 2.0, 3.0]).unwrap();
        let serialized = serde_json::to_string(&atom).unwrap();
        let deserialized: Atom = serde_json::from_str(&serialized).unwrap();

        assert_eq!(atom, deserialized);
    }

    // -------------------------------------------------------------------------
    // Reference system tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_h2_system() {
        // H2 with bond length 0.74 Angstrom
        let h1 = Atom::new_angstrom(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new_angstrom(1, [0.0, 0.0, 0.74]).unwrap();

        let bond_length_bohr = h1.distance_to(&h2);
        let bond_length_angstrom = bond_length_bohr * BOHR_TO_ANGSTROM;

        // Use reasonable tolerance for round-trip conversion
        assert_abs_diff_eq!(bond_length_angstrom, 0.74, epsilon = 1e-8);
    }

    #[test]
    fn test_h2o_atoms() {
        // Build H2O geometry
        let o = Atom::from_symbol("O", [0.0, 0.0, 0.0]).unwrap();
        let h1 = Atom::from_symbol("H", [1.43, 1.11, 0.0]).unwrap();
        let h2 = Atom::from_symbol("H", [-1.43, 1.11, 0.0]).unwrap();

        assert_eq!(o.atomic_number, 8);
        assert_eq!(h1.atomic_number, 1);
        assert_eq!(h2.atomic_number, 1);

        // Total nuclear charge: 8 + 1 + 1 = 10
        let total_charge: f64 = [&o, &h1, &h2].iter().map(|a| a.nuclear_charge()).sum();
        assert_abs_diff_eq!(total_charge, 10.0, epsilon = 1e-10);
    }
}
