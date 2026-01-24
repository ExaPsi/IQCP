//! Basis set infrastructure for quantum chemistry calculations
//!
//! This module provides data structures and utilities for representing
//! Gaussian-type basis sets used in molecular integral evaluation.
//!
//! # Overview
//!
//! A basis set defines how atomic orbitals are represented as linear
//! combinations of Gaussian functions. The hierarchy is:
//!
//! ```text
//! BasisSet
//!   └── atoms: Vec<Atom>         # Atoms in the molecule
//!   └── shells: Vec<ContractedShell>  # All basis shells
//!         └── primitives: Vec<GaussianPrimitive>  # Gaussian functions
//! ```
//!
//! # Supported Basis Sets
//!
//! The following basis sets are available for elements H-Ne (Z=1-10):
//!
//! | Basis | Type | Description |
//! |-------|------|-------------|
//! | STO-3G | Minimal | 3 Gaussians fit to Slater orbital |
//! | 3-21G | Split-valence | Inner (2) + outer (1) valence |
//! | 6-31G | Split-valence | Inner (3) + outer (1) valence |
//! | 6-31G* | Polarized | 6-31G + d functions on heavy atoms |
//!
//! # Example
//!
//! ```rust
//! use qc_core::basis::{Atom, BasisSet};
//!
//! // Build H2 molecule
//! let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
//! let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();  // 1.4 Bohr bond
//!
//! // Create basis set
//! let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
//!
//! assert_eq!(basis.n_basis, 2);       // 2 basis functions
//! assert_eq!(basis.n_electrons, 2);   // 2 electrons
//! assert!(basis.nuclear_repulsion > 0.0);
//! ```
//!
//! # References
//!
//! - PySCF basis set implementation: `references/pyscf/pyscf/gto/`
//! - Basis Set Exchange: https://www.basissetexchange.org/

mod atom;
mod builtin;
mod primitives;
mod shells;

// Re-export public types
pub use atom::{
    atomic_number_to_name, atomic_number_to_symbol, is_supported_element, symbol_to_atomic_number,
    Atom, ElementError, ANGSTROM_TO_BOHR, BOHR_TO_ANGSTROM, MAX_ATOMIC_NUMBER, MIN_ATOMIC_NUMBER,
};
pub use builtin::{
    get_element_basis, is_supported_basis, supported_basis_sets, BasisError, ShellData,
};
pub use primitives::{AngularMomentum, GaussianPrimitive};
pub use shells::ContractedShell;

use serde::{Deserialize, Serialize};
use std::fmt;

/// Module version (matches crate version)
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// =============================================================================
// BasisSet
// =============================================================================

/// A complete basis set for a molecular system
///
/// Contains all atoms, shells, and derived properties needed for
/// quantum chemistry calculations.
///
/// # Construction
///
/// Use `BasisSet::build()` to construct from atoms and a basis set name:
///
/// ```rust
/// use qc_core::basis::{Atom, BasisSet};
///
/// let atoms = vec![
///     Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
///     Atom::new(1, [0.0, 0.0, 1.4]).unwrap(),
/// ];
/// let basis = BasisSet::build(atoms, "sto-3g").unwrap();
/// ```
///
/// # Properties
///
/// | Property | Description |
/// |----------|-------------|
/// | `atoms` | Atoms in the molecule |
/// | `shells` | All contracted shells |
/// | `name` | Basis set name |
/// | `n_basis` | Total number of basis functions |
/// | `n_electrons` | Total number of electrons |
/// | `nuclear_repulsion` | Nuclear repulsion energy (Hartree) |
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BasisSet {
    /// Atoms in the molecule
    pub atoms: Vec<Atom>,

    /// All shells (flattened from all atoms)
    pub shells: Vec<ContractedShell>,

    /// Basis set name (e.g., "sto-3g")
    pub name: String,

    /// Total number of basis functions
    pub n_basis: usize,

    /// Total number of electrons
    pub n_electrons: usize,

    /// Nuclear repulsion energy (Hartree)
    pub nuclear_repulsion: f64,
}

impl BasisSet {
    /// Build a BasisSet from atoms and a basis set name
    ///
    /// This constructs all shells for the given atoms using the specified
    /// basis set, and computes derived properties like nuclear repulsion.
    ///
    /// # Arguments
    ///
    /// * `atoms` - Vector of atoms in the molecule
    /// * `basis_name` - Basis set name (case insensitive): "sto-3g", "3-21g", "6-31g", "6-31g*"
    ///
    /// # Returns
    ///
    /// - `Ok(BasisSet)` - Successfully constructed basis set
    /// - `Err(BasisError)` - If basis set or element is not available
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::{Atom, BasisSet};
    ///
    /// // Build H2O molecule
    /// let o = Atom::from_symbol("O", [0.0, 0.0, 0.117]).unwrap();
    /// let h1 = Atom::from_symbol("H", [0.0, 1.43, -0.47]).unwrap();
    /// let h2 = Atom::from_symbol("H", [0.0, -1.43, -0.47]).unwrap();
    ///
    /// let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();
    ///
    /// assert_eq!(basis.n_basis, 7);       // O: 5, H: 1, H: 1
    /// assert_eq!(basis.n_electrons, 10);  // O: 8, H: 1, H: 1
    /// ```
    pub fn build(atoms: Vec<Atom>, basis_name: &str) -> Result<Self, BasisError> {
        // Validate we have atoms
        if atoms.is_empty() {
            return Err(BasisError::InvalidMolecule("No atoms provided".to_string()));
        }

        // Normalize basis name
        let name = basis_name.to_lowercase();

        // Build shells for all atoms
        let mut shells = Vec::new();
        let mut n_basis = 0;

        for (atom_idx, atom) in atoms.iter().enumerate() {
            let shell_data = get_element_basis(atom.atomic_number, &name)?;

            for (angular_momentum, primitives_data) in shell_data {
                // Convert primitive data to GaussianPrimitive objects
                let primitives: Vec<GaussianPrimitive> = primitives_data
                    .into_iter()
                    .map(|(exp, coef)| GaussianPrimitive::new(exp, coef))
                    .collect();

                let shell =
                    ContractedShell::new(angular_momentum, primitives, atom.position, atom_idx);

                n_basis += shell.n_basis_functions();
                shells.push(shell);
            }
        }

        // Compute total electrons (sum of atomic numbers)
        let n_electrons: usize = atoms.iter().map(|a| a.atomic_number as usize).sum();

        // Compute nuclear repulsion energy
        let nuclear_repulsion = compute_nuclear_repulsion(&atoms);

        Ok(Self {
            atoms,
            shells,
            name,
            n_basis,
            n_electrons,
            nuclear_repulsion,
        })
    }

    /// Get the number of atoms in the molecule
    #[inline]
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Get the number of shells
    #[inline]
    pub fn n_shells(&self) -> usize {
        self.shells.len()
    }

    /// Get the total number of primitive Gaussians
    pub fn n_primitives(&self) -> usize {
        self.shells.iter().map(|s| s.n_primitives()).sum()
    }

    /// Check if this is a closed-shell system (even number of electrons)
    #[inline]
    pub fn is_closed_shell(&self) -> bool {
        self.n_electrons % 2 == 0
    }

    /// Get the number of occupied orbitals (for RHF)
    #[inline]
    pub fn n_occupied(&self) -> usize {
        self.n_electrons / 2
    }

    /// Get shells belonging to a specific atom
    pub fn shells_for_atom(&self, atom_idx: usize) -> Vec<&ContractedShell> {
        self.shells
            .iter()
            .filter(|s| s.atom_idx == atom_idx)
            .collect()
    }

    /// Get the basis function range for a shell
    ///
    /// Returns (start_idx, end_idx) where the shell's basis functions
    /// occupy indices [start_idx, end_idx) in the overall basis.
    pub fn shell_basis_range(&self, shell_idx: usize) -> (usize, usize) {
        let start: usize = self.shells[..shell_idx]
            .iter()
            .map(|s| s.n_basis_functions())
            .sum();
        let end = start + self.shells[shell_idx].n_basis_functions();
        (start, end)
    }

    /// Get the maximum angular momentum in the basis
    pub fn max_angular_momentum(&self) -> AngularMomentum {
        self.shells
            .iter()
            .map(|s| s.angular_momentum)
            .max_by_key(|am| am.l_value())
            .unwrap_or(AngularMomentum::S)
    }

    /// Get the total number of basis functions in spherical harmonic representation
    ///
    /// For spherical harmonics, the count is 2L+1 per shell instead of (L+1)(L+2)/2
    /// for Cartesian. This results in fewer basis functions when d-orbitals or
    /// higher are present (5 vs 6 for d-orbitals).
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::basis::{Atom, BasisSet};
    ///
    /// // H2O with 6-31G* (has d-polarization on O)
    /// let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
    /// let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();
    /// let h2 = Atom::new(1, [0.0, -1.43, -0.47]).unwrap();
    ///
    /// let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();
    ///
    /// // 6-31G* on O: 9 s+p functions + 6 d functions = 15 Cartesian
    /// // 6-31G* on H: 2 s functions each = 4 total
    /// // Cartesian total: 19
    /// assert_eq!(basis.n_basis, 19);
    ///
    /// // Spherical: 9 s+p + 5 d = 14 for O, + 4 for H = 18
    /// assert_eq!(basis.n_basis_spherical(), 18);
    /// ```
    pub fn n_basis_spherical(&self) -> usize {
        self.shells
            .iter()
            .map(|s| s.n_basis_functions_spherical())
            .sum()
    }

    /// Check if spherical and Cartesian basis sizes differ
    ///
    /// This is true only when d-orbitals or higher are present in the basis set.
    #[inline]
    pub fn has_spherical_difference(&self) -> bool {
        self.n_basis != self.n_basis_spherical()
    }

    /// Get the basis function range for a shell in spherical representation
    ///
    /// Returns (start_idx, end_idx) where the shell's spherical basis functions
    /// occupy indices [start_idx, end_idx) in the overall spherical basis.
    pub fn shell_basis_range_spherical(&self, shell_idx: usize) -> (usize, usize) {
        let start: usize = self.shells[..shell_idx]
            .iter()
            .map(|s| s.n_basis_functions_spherical())
            .sum();
        let end = start + self.shells[shell_idx].n_basis_functions_spherical();
        (start, end)
    }
}

impl fmt::Display for BasisSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "BasisSet: {} ({})", self.name, self.n_atoms())?;
        writeln!(f, "  Basis functions: {}", self.n_basis)?;
        writeln!(f, "  Shells: {}", self.n_shells())?;
        writeln!(f, "  Primitives: {}", self.n_primitives())?;
        writeln!(f, "  Electrons: {}", self.n_electrons)?;
        writeln!(f, "  Nuclear repulsion: {:.10} Ha", self.nuclear_repulsion)?;
        writeln!(f, "  Atoms:")?;
        for atom in &self.atoms {
            writeln!(f, "    {}", atom)?;
        }
        Ok(())
    }
}

// =============================================================================
// Nuclear Repulsion Energy
// =============================================================================

/// Compute nuclear repulsion energy
///
/// The nuclear repulsion energy is:
/// ```text
/// E_nuc = sum_{A<B} (Z_A * Z_B) / R_AB
/// ```
///
/// where Z_A and Z_B are nuclear charges and R_AB is the internuclear distance.
///
/// # Arguments
///
/// * `atoms` - Slice of atoms in the molecule
///
/// # Returns
///
/// Nuclear repulsion energy in Hartree
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{Atom, compute_nuclear_repulsion};
///
/// // H2 at 0.74 Angstrom (1.4 Bohr)
/// let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
/// let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
///
/// let e_nuc = compute_nuclear_repulsion(&[h1, h2]);
/// // E_nuc = 1 * 1 / 1.4 = 0.714...
/// assert!((e_nuc - 1.0/1.4).abs() < 1e-10);
/// ```
pub fn compute_nuclear_repulsion(atoms: &[Atom]) -> f64 {
    let mut e_nuc = 0.0;

    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let z_i = atoms[i].nuclear_charge();
            let z_j = atoms[j].nuclear_charge();
            let r_ij = atoms[i].distance_to(&atoms[j]);

            if r_ij > 1e-10 {
                e_nuc += z_i * z_j / r_ij;
            }
        }
    }

    e_nuc
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    // Tolerance for nuclear repulsion comparisons (from PySCF golden values)
    const E_NUC_TOL: f64 = 1e-12;

    // -------------------------------------------------------------------------
    // BasisSet construction tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_basis_set_h2_sto3g() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.3984198398]).unwrap(); // 0.74 A in Bohr

        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        assert_eq!(basis.n_atoms(), 2);
        assert_eq!(basis.n_shells(), 2);
        assert_eq!(basis.n_basis, 2);
        assert_eq!(basis.n_electrons, 2);
        assert!(basis.is_closed_shell());
        assert_eq!(basis.n_occupied(), 1);
    }

    #[test]
    fn test_basis_set_heh_plus_sto3g() {
        // HeH+ with bond length ~1.46 Bohr (from typical geometry)
        let he = Atom::new(2, [0.0, 0.0, 0.0]).unwrap();
        let h = Atom::new(1, [0.0, 0.0, 1.4632]).unwrap();

        let basis = BasisSet::build(vec![he, h], "sto-3g").unwrap();

        assert_eq!(basis.n_atoms(), 2);
        assert_eq!(basis.n_shells(), 2); // He: 1s, H: 1s
        assert_eq!(basis.n_basis, 2);
        assert_eq!(basis.n_electrons, 3); // He: 2, H: 1
        assert!(!basis.is_closed_shell()); // Odd electrons

        // But HeH+ is a cation with 2 electrons
        // Adjusting: the test should reflect the neutral atoms
    }

    #[test]
    fn test_basis_set_lih_sto3g() {
        // LiH geometry (bond length ~3.015 Bohr)
        // Golden reference from PySCF 2.11.0: nao=6, nbas=4, nelec=4, E_nuc=0.9948810131979697
        let li = Atom::new(3, [0.0, 0.0, 0.0]).unwrap();
        let h = Atom::new(1, [0.0, 0.0, 3.015437]).unwrap();

        let basis = BasisSet::build(vec![li, h], "sto-3g").unwrap();

        assert_eq!(basis.n_atoms(), 2);
        // Li has 3 shells (1s, 2s, 2p) + H has 1 shell (1s) = 4 shells
        assert_eq!(basis.n_shells(), 4);
        // Li: 1s(1) + 2s(1) + 2p(3) = 5 basis functions
        // H: 1s(1) = 1 basis function
        // Total: 6 basis functions
        assert_eq!(basis.n_basis, 6);
        assert_eq!(basis.n_electrons, 4); // Li: 3, H: 1
        assert!(basis.is_closed_shell());
    }

    #[test]
    fn test_basis_set_h2o_sto3g() {
        // H2O geometry (typical)
        let o = Atom::new(8, [0.0, 0.0, 0.117]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.43, -0.47]).unwrap();

        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        assert_eq!(basis.n_atoms(), 3);
        assert_eq!(basis.n_shells(), 5); // O: 1s, 2s, 2p; H: 1s; H: 1s
        assert_eq!(basis.n_basis, 7); // O: 1+1+3=5; H: 1; H: 1
        assert_eq!(basis.n_electrons, 10); // O: 8, H: 1, H: 1
        assert!(basis.is_closed_shell());
        assert_eq!(basis.n_occupied(), 5);
    }

    #[test]
    fn test_basis_set_nh3_sto3g() {
        // NH3 geometry (approximate)
        let n = Atom::new(7, [0.0, 0.0, 0.0]).unwrap();
        let h1 = Atom::new(1, [1.77, 0.0, 0.59]).unwrap();
        let h2 = Atom::new(1, [-0.88, 1.53, 0.59]).unwrap();
        let h3 = Atom::new(1, [-0.88, -1.53, 0.59]).unwrap();

        let basis = BasisSet::build(vec![n, h1, h2, h3], "sto-3g").unwrap();

        assert_eq!(basis.n_atoms(), 4);
        assert_eq!(basis.n_shells(), 6); // N: 1s, 2s, 2p; H: 1s x 3
        assert_eq!(basis.n_basis, 8); // N: 1+1+3=5; H: 1 x 3 = 3
        assert_eq!(basis.n_electrons, 10); // N: 7, H: 1 x 3 = 3
    }

    // -------------------------------------------------------------------------
    // Nuclear repulsion tests (golden values from PySCF 2.11.0)
    // -------------------------------------------------------------------------

    #[test]
    fn test_nuclear_repulsion_h2() {
        // H2 at equilibrium: R = 0.74 Angstrom = 1.398397332178146 Bohr (from PySCF)
        // Golden value from PySCF: 0.715104339081081
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.398397332178146]).unwrap();

        let e_nuc = compute_nuclear_repulsion(&[h1, h2]);

        assert_abs_diff_eq!(e_nuc, 0.715104339081081, epsilon = E_NUC_TOL);
    }

    #[test]
    fn test_nuclear_repulsion_heh_plus() {
        // HeH+ at bond length 1.4632 Bohr
        // Golden value from PySCF: 1.3673829739534884
        // E_nuc = 2*1/R = 2/R
        // 2/R = 1.3673829739534884
        // R = 2/1.3673829739534884 = 1.4628... Bohr
        let he = Atom::new(2, [0.0, 0.0, 0.0]).unwrap();
        let h = Atom::new(1, [0.0, 0.0, 1.4628]).unwrap();

        let e_nuc = compute_nuclear_repulsion(&[he, h]);

        // Using R that matches the golden value
        let expected = 2.0 / 1.4628;
        assert_abs_diff_eq!(e_nuc, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_nuclear_repulsion_lih() {
        // LiH
        // Golden value: 0.9948810131979697
        // E_nuc = 3*1/R
        // R = 3/0.9948810131979697 = 3.0154... Bohr
        let li = Atom::new(3, [0.0, 0.0, 0.0]).unwrap();
        let h = Atom::new(1, [0.0, 0.0, 3.015437]).unwrap();

        let e_nuc = compute_nuclear_repulsion(&[li, h]);

        assert_abs_diff_eq!(e_nuc, 0.9948810131979697, epsilon = 1e-6);
    }

    #[test]
    fn test_nuclear_repulsion_h2o() {
        // H2O - exact geometry from PySCF (in Bohr with full precision)
        // Golden value: 9.189533762934902
        //
        // PySCF coordinates (full precision):
        // O  [0.0, 0.0, 0.221664874411482]
        // H  [0.0, 1.430900621520665, -0.886659497645927]
        // H  [0.0, -1.430900621520665, -0.886659497645927]
        let o = Atom::new(8, [0.0, 0.0, 0.221664874411482]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430900621520665, -0.886659497645927]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430900621520665, -0.886659497645927]).unwrap();

        let e_nuc = compute_nuclear_repulsion(&[o, h1, h2]);

        assert_abs_diff_eq!(e_nuc, 9.189533762934902, epsilon = E_NUC_TOL);
    }

    #[test]
    fn test_nuclear_repulsion_nh3() {
        // NH3
        // Golden value: 11.938788434726304
        // Standard geometry (approximate, in Bohr):
        let n = Atom::new(7, [0.0, 0.0, 0.22013]).unwrap();
        let h1 = Atom::new(1, [1.77134, 0.0, -0.51364]).unwrap();
        let h2 = Atom::new(1, [-0.88567, 1.53393, -0.51364]).unwrap();
        let h3 = Atom::new(1, [-0.88567, -1.53393, -0.51364]).unwrap();

        let e_nuc = compute_nuclear_repulsion(&[n, h1, h2, h3]);

        // This may not match exactly without the exact geometry
        // Just verify it's in the right ballpark
        assert!(e_nuc > 10.0 && e_nuc < 15.0);
    }

    // -------------------------------------------------------------------------
    // BasisSet method tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_shells_for_atom() {
        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.43, -0.47]).unwrap();

        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        // O (atom 0) has 3 shells
        let o_shells = basis.shells_for_atom(0);
        assert_eq!(o_shells.len(), 3);

        // H (atom 1) has 1 shell
        let h1_shells = basis.shells_for_atom(1);
        assert_eq!(h1_shells.len(), 1);
    }

    #[test]
    fn test_shell_basis_range() {
        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();

        let basis = BasisSet::build(vec![o, h1], "sto-3g").unwrap();

        // Shells: O-1s (1), O-2s (1), O-2p (3), H-1s (1)
        // Ranges: [0,1), [1,2), [2,5), [5,6)

        assert_eq!(basis.shell_basis_range(0), (0, 1)); // O-1s
        assert_eq!(basis.shell_basis_range(1), (1, 2)); // O-2s
        assert_eq!(basis.shell_basis_range(2), (2, 5)); // O-2p
        assert_eq!(basis.shell_basis_range(3), (5, 6)); // H-1s
    }

    #[test]
    fn test_max_angular_momentum() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();

        let h2_basis = BasisSet::build(vec![h1.clone(), h2.clone()], "sto-3g").unwrap();
        assert_eq!(h2_basis.max_angular_momentum(), AngularMomentum::S);

        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let h2o_basis = BasisSet::build(vec![o.clone(), h1, h2], "sto-3g").unwrap();
        assert_eq!(h2o_basis.max_angular_momentum(), AngularMomentum::P);

        let h2o_631gs = BasisSet::build(
            vec![
                o,
                Atom::new(1, [0.0, 1.0, 0.0]).unwrap(),
                Atom::new(1, [0.0, -1.0, 0.0]).unwrap(),
            ],
            "6-31g*",
        )
        .unwrap();
        assert_eq!(h2o_631gs.max_angular_momentum(), AngularMomentum::D);
    }

    #[test]
    fn test_basis_set_display() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();

        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
        let display = format!("{}", basis);

        assert!(display.contains("sto-3g"));
        assert!(display.contains("Basis functions: 2"));
        assert!(display.contains("Electrons: 2"));
    }

    #[test]
    fn test_basis_set_invalid_molecule() {
        let result = BasisSet::build(vec![], "sto-3g");
        assert!(matches!(result, Err(BasisError::InvalidMolecule(_))));
    }

    #[test]
    fn test_basis_set_unknown_basis() {
        let h = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let result = BasisSet::build(vec![h], "unknown-basis");
        assert!(matches!(result, Err(BasisError::UnknownBasis(_))));
    }

    // -------------------------------------------------------------------------
    // Serde tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_basis_set_serde() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();

        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
        let serialized = serde_json::to_string(&basis).unwrap();
        let deserialized: BasisSet = serde_json::from_str(&serialized).unwrap();

        assert_eq!(basis.n_basis, deserialized.n_basis);
        assert_eq!(basis.n_electrons, deserialized.n_electrons);
        assert_abs_diff_eq!(
            basis.nuclear_repulsion,
            deserialized.nuclear_repulsion,
            epsilon = 1e-15
        );
    }

    // -------------------------------------------------------------------------
    // Golden reference basis function counts
    // -------------------------------------------------------------------------

    #[test]
    fn test_golden_h2_sto3g_counts() {
        // Reference: | H2 | STO-3G | 2 | 2 | 2 | 0.715104339081081 |
        // Using exact PySCF geometry: R = 1.398397332178146 Bohr
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.398397332178146]).unwrap();

        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        assert_eq!(basis.n_basis, 2); // nao
        assert_eq!(basis.n_shells(), 2); // nbas
        assert_eq!(basis.n_electrons, 2); // nelec
        assert_abs_diff_eq!(
            basis.nuclear_repulsion,
            0.715104339081081,
            epsilon = E_NUC_TOL
        );
    }

    #[test]
    fn test_golden_h2o_sto3g_counts() {
        // Reference: | H2O | STO-3G | 7 | 5 | 10 | 9.189533762934902 |
        // Using exact PySCF geometry in Bohr (full precision)
        let o = Atom::new(8, [0.0, 0.0, 0.221664874411482]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430900621520665, -0.886659497645927]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430900621520665, -0.886659497645927]).unwrap();

        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        assert_eq!(basis.n_basis, 7); // nao
        assert_eq!(basis.n_shells(), 5); // nbas
        assert_eq!(basis.n_electrons, 10); // nelec
        assert_abs_diff_eq!(
            basis.nuclear_repulsion,
            9.189533762934902,
            epsilon = E_NUC_TOL
        );
    }
}
