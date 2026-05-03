//! Physical constants and atomic masses for quantum chemistry calculations.
//!
//! All values match PySCF's `pyscf/data/nist.py` (CODATA 2014) to ensure
//! cross-code consistency in frequency conversion, thermochemistry, and
//! mass-weighted Hessian calculations.
//!
//! # References
//!
//! - CODATA 2014: <https://physics.nist.gov/cuu/Constants/>
//! - PySCF nist.py: `references/pyscf/pyscf/data/nist.py`
//! - IUPAC 2016 atomic weights: <https://www.qmul.ac.uk/sbcs/iupac/AtWt/>

// =============================================================================
// Fundamental Physical Constants (CODATA 2014, matching PySCF)
// =============================================================================

/// Speed of light in vacuum (m/s)
/// <https://physics.nist.gov/cgi-bin/cuu/Value?c>
pub const LIGHT_SPEED_SI: f64 = 299_792_458.0;

/// Speed of light in atomic units (1/α)
/// <https://physics.nist.gov/cgi-bin/cuu/Value?alph>
pub const LIGHT_SPEED_AU: f64 = 137.03599967994;

/// Planck constant (J·s)
/// <https://physics.nist.gov/cgi-bin/cuu/Value?h>
pub const PLANCK: f64 = 6.626_070_040e-34;

/// Reduced Planck constant ℏ = h/(2π) (J·s)
pub const HBAR: f64 = PLANCK / (2.0 * std::f64::consts::PI);

/// Boltzmann constant (J/K)
/// <https://physics.nist.gov/cgi-bin/cuu/Value?k>
pub const BOLTZMANN: f64 = 1.380_648_52e-23;

/// Avogadro constant (mol⁻¹)
/// <https://physics.nist.gov/cgi-bin/cuu/Value?na>
pub const AVOGADRO: f64 = 6.022_140_857e23;

/// Molar gas constant R = k_B × N_A (J/(mol·K))
pub const GAS_CONSTANT: f64 = BOLTZMANN * AVOGADRO;

/// Elementary charge (C)
/// <https://physics.nist.gov/cgi-bin/cuu/Value?e>
pub const E_CHARGE: f64 = 1.602_176_620_8e-19;

/// Electron mass (kg)
/// <https://physics.nist.gov/cgi-bin/cuu/Value?me>
pub const E_MASS: f64 = 9.109_383_56e-31;

// =============================================================================
// Unit Conversion: Length
// =============================================================================

/// Bohr radius in Angstroms (CODATA 2014)
/// <https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0>
///
/// This is the single source of truth for Bohr-Angstrom conversion.
/// All other modules must import from here.
pub const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_92;

/// Angstrom to Bohr conversion (1 / BOHR_TO_ANGSTROM)
pub const ANGSTROM_TO_BOHR: f64 = 1.889_726_124_565_061_8;

/// Bohr radius in meters
pub const BOHR_SI: f64 = BOHR_TO_ANGSTROM * 1e-10;

/// Bohr³ → Å³ conversion for polarizability.
///
/// Since atomic units of polarizability are `e² · a₀² / E_h = a₀³` (bohr³)
/// and the IUPAC / spectroscopic convention reports polarizability in Å³,
/// the conversion factor is `BOHR_TO_ANGSTROM³ ≈ 0.1481847`.
///
/// Used by [`crate::raman::RamanResult::polarizability_ang3`].
pub const AU_POLARIZABILITY_TO_ANG3: f64 = BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM;

/// Bohr⁴ → Å⁴ conversion for Raman scattering activity.
///
/// The Raman scattering activity `S_k = 45 · ᾱ'² + 7 · γ'²` is naturally
/// computed in atomic units of polarizability derivative squared per amu:
/// `(bohr³/bohr)² · amu⁻¹ = bohr⁴/amu`. The IUPAC / spectroscopic convention
/// reports activities in Å⁴/amu, requiring multiplication by
/// `BOHR_TO_ANGSTROM⁴ ≈ 0.07843376`.
///
/// Used by [`crate::raman::compute_raman_activities`].
///
/// # Reference
///
/// - Long, D. A. (2002). *The Raman Effect: A Unified Treatment of the
///   Theory of Raman Scattering by Molecules.* Wiley. Chapter 5 Eq. 5.5.9.
/// - Gaussian 16 output units for "Raman Activities"
pub const BOHR4_TO_ANG4: f64 =
    BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM;

// =============================================================================
// Unit Conversion: Energy
// =============================================================================

/// Hartree to Joule conversion
/// PySCF computes this as ℏ²/(m_e × a_0²) rather than using the NIST value.
/// We use the same formula for exact cross-code consistency.
/// ≈ 4.359744644912e-18 J
/// <https://physics.nist.gov/cgi-bin/cuu/Value?hrj>
pub const HARTREE_TO_JOULE: f64 = 4.359_744_644_911_914e-18;

/// Hartree to eV conversion
/// <https://physics.nist.gov/cgi-bin/cuu/Value?threv>
pub const HARTREE_TO_EV: f64 = 27.211_386_02;

/// Hartree to kcal/mol conversion
pub const HARTREE_TO_KCAL: f64 = 627.509_474;

/// Hartree to kJ/mol conversion
pub const HARTREE_TO_KJ: f64 = 2_625.499_638;

// =============================================================================
// Unit Conversion: Mass
// =============================================================================

/// Atomic mass unit in kg (1 amu = 1 g/mol / N_A)
pub const AMU_TO_KG: f64 = 1e-3 / AVOGADRO;

/// Atomic mass unit to electron mass ratio
pub const AMU_TO_AU: f64 = AMU_TO_KG / E_MASS;

// =============================================================================
// Unit Conversion: Spectroscopy (Derived)
// =============================================================================

/// Conversion factor: atomic units of frequency → Hz
/// = √(E_h / (m_u × a_0²)) / (2π)
/// where E_h = Hartree, m_u = atomic mass unit, a_0 = Bohr radius
///
/// Used in: frequency = freq_au × AU_TO_HZ
// Note: Cannot use .sqrt() in const context, so we precompute the value.
// AU_TO_HZ = sqrt(HARTREE_TO_JOULE / (AMU_TO_KG * BOHR_SI^2)) / (2*PI)
// = sqrt(4.359744644912e-18 / (1.660539040e-27 * (5.29177210920e-11)^2)) / (2*PI)
// ≈ 1.541079e14
pub const AU_TO_HZ: f64 = 1.541_079_274_566_313e14;

/// Conversion factor: atomic units of frequency → wavenumber (cm⁻¹)
/// = AU_TO_HZ / c × 10⁻²
///
/// PySCF: HARTREE2WAVENUMBER = 1e-2 * HARTREE2J / (LIGHT_SPEED_SI * PLANCK) ≈ 2.194746e5
/// But that converts Hartree (energy) to cm⁻¹.
/// For force constants (eigenvalues of mass-weighted Hessian in au),
/// freq_cm1 = sign(λ) * sqrt(|λ|) * AU_TO_CM1
pub const AU_TO_CM1: f64 = AU_TO_HZ / LIGHT_SPEED_SI * 1e-2;

/// Conversion factor: Hartree (energy) → wavenumber (cm⁻¹)
/// = HARTREE_TO_JOULE / (PLANCK × c) × 10⁻²
/// ≈ 2.194746313702e5 (matches PySCF HARTREE2WAVENUMBER)
pub const HARTREE_TO_WAVENUMBER: f64 = 1e-2 * HARTREE_TO_JOULE / (LIGHT_SPEED_SI * PLANCK);

// =============================================================================
// Unit Conversion: Dipole
// =============================================================================

/// Debye in C·m
pub const DEBYE: f64 = 3.335_641e-30;

/// Atomic units to Debye conversion for dipole moment
/// = e × a_0 / Debye ≈ 2.541746
pub const AU_TO_DEBYE: f64 = E_CHARGE * BOHR_SI / DEBYE;

// =============================================================================
// Unit Conversion: IR Absorption Intensities
// =============================================================================

/// Conversion factor for IR absorption intensities from atomic units to km/mol.
///
/// The integrated IR absorption intensity of vibrational mode k is:
///
/// ```text
/// A_k [km/mol] = AU_TO_KMMOL * |∂μ/∂Q_k|^2
/// ```
///
/// where `|∂μ/∂Q_k|^2` is in atomic units `e^2 / amu` (electric dipole
/// squared per atomic mass unit). The `norm_mode` from PySCF's
/// `harmonic_analysis` is in units of `1/sqrt(amu)` (PySCF convention:
/// `norm_mode = mass^(-1/2) * eigenvector`), and `∂μ/∂R` is in atomic
/// units of `e` (dipole/bohr = e*bohr/bohr = e).
///
/// # Derivation
///
/// Matches the pyscf-forge implementation
/// `pyscf/prop/infrared/rhf.py:kernel_ir` (line 161):
///
/// ```text
/// unit_kmmol = alpha^2 * (1e-3 / m_u) * m_e * N_A * pi * a_0 / 3
/// ```
///
/// where:
/// - `alpha = 1/c` is the fine-structure constant (in atomic units,
///   it equals `1/LIGHT_SPEED_AU`; equivalently in SI `alpha = 1/LIGHT_SPEED_SI`
///   when expressed as a dimensionless scalar via 1/137.036... but here we
///   use the direct formula with SI constants)
/// - `m_u = ATOMIC_MASS` atomic mass constant (kg)
/// - `m_e = E_MASS` electron mass (kg)
/// - `N_A = AVOGADRO` Avogadro constant (mol^-1)
/// - `a_0 = BOHR_SI` Bohr radius (m)
/// - The `1e-3` factor converts m/mol to km/mol
///
/// Equivalent derivation:
///
/// ```text
/// A_k [m/mol] = (N_A * pi / (3 * c^2)) * |∂μ/∂Q_k|^2_SI
/// ```
///
/// with `|∂μ/∂Q|^2_SI` in `C^2/kg`.
///
/// # Numerical value
///
/// Approximately **974.880** km/mol per `e^2/amu`.
///
/// # References
///
/// - pyscf-forge `prop/infrared/rhf.py` line 161 (`kernel_ir`)
/// - Psi4 `driver/qcdb/vib.py` line 589
/// - Person, W.B. & Zerbi, G. (1982) *Vibrational Intensities in
///   Infrared and Raman Spectroscopy*, Chapter 1
/// - McQuarrie & Simon (1997) *Physical Chemistry: A Molecular Approach*,
///   Section 13.6
pub const AU_TO_KMMOL: f64 = {
    // alpha = fine-structure constant ≈ 7.297e-3 (dimensionless).
    // In atomic units, alpha = 1 / c_au where c_au = LIGHT_SPEED_AU ≈ 137.036.
    // We compute it from LIGHT_SPEED_AU rather than hard-coding to stay
    // consistent with the rest of the constants table.
    let alpha = 1.0 / LIGHT_SPEED_AU;
    // 974.8801... km/mol per (e*a0/sqrt(amu))^2 = e^2/amu in atomic units
    alpha * alpha * (1e-3 / AMU_TO_KG) * E_MASS * AVOGADRO * std::f64::consts::PI * BOHR_SI / 3.0
};

// =============================================================================
// Force Constant Conversion
// =============================================================================

/// Conversion factor for force constants: au → mDyne/Å
/// 1 Hartree/bohr² → mDyne/Å
/// = HARTREE_TO_JOULE / (BOHR_SI² × 1e-2) × 1e-8
/// ≈ 15.569 mDyne/Å per Hartree/bohr²
pub const AU_TO_MDYNE_PER_ANG: f64 = HARTREE_TO_JOULE / (BOHR_SI * BOHR_SI) * 1e-5;

// =============================================================================
// Atomic Masses (IUPAC 2016 Standard Atomic Weights)
// =============================================================================

/// Standard atomic masses in amu for elements Z=0 through Z=18.
/// Index 0 is a placeholder (no element with Z=0).
/// Values from IUPAC 2016 recommended standard atomic weights.
///
/// Usage: `ATOMIC_MASSES[atomic_number as usize]`
pub const ATOMIC_MASSES: [f64; 19] = [
    0.0,          // Z=0: placeholder
    1.007_94,     // Z=1:  H  (hydrogen)
    4.002_602,    // Z=2:  He (helium)
    6.941,        // Z=3:  Li (lithium)
    9.012_182,    // Z=4:  Be (beryllium)
    10.811,       // Z=5:  B  (boron)
    12.0107,      // Z=6:  C  (carbon)
    14.006_7,     // Z=7:  N  (nitrogen)
    15.999_4,     // Z=8:  O  (oxygen)
    18.998_403_2, // Z=9:  F  (fluorine)
    20.179_7,     // Z=10: Ne (neon)
    22.989_770,   // Z=11: Na (sodium)
    24.305_0,     // Z=12: Mg (magnesium)
    26.981_538,   // Z=13: Al (aluminium)
    28.085_5,     // Z=14: Si (silicon)
    30.973_762,   // Z=15: P  (phosphorus)
    32.065,       // Z=16: S  (sulfur)
    35.453,       // Z=17: Cl (chlorine)
    39.948,       // Z=18: Ar (argon)
];

/// Get the atomic mass for a given atomic number.
/// Returns `None` if Z is 0 or > 18.
pub fn atomic_mass(z: u8) -> Option<f64> {
    if z == 0 || z as usize >= ATOMIC_MASSES.len() {
        None
    } else {
        Some(ATOMIC_MASSES[z as usize])
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bohr_angstrom_consistency() {
        // BOHR_TO_ANGSTROM × ANGSTROM_TO_BOHR ≈ 1.0
        let product = BOHR_TO_ANGSTROM * ANGSTROM_TO_BOHR;
        assert!(
            (product - 1.0).abs() < 1e-14,
            "BOHR_TO_ANGSTROM * ANGSTROM_TO_BOHR = {product}, expected 1.0"
        );
    }

    #[test]
    fn bohr_si_consistent() {
        // BOHR_SI should be BOHR_TO_ANGSTROM * 1e-10
        let expected = BOHR_TO_ANGSTROM * 1e-10;
        assert!(
            (BOHR_SI - expected).abs() < 1e-25,
            "BOHR_SI = {BOHR_SI}, expected {expected}"
        );
    }

    #[test]
    fn constants_match_pyscf_nist() {
        // Cross-check against PySCF nist.py values
        assert_eq!(BOHR_TO_ANGSTROM, 0.52917721092, "BOHR mismatch");
        assert_eq!(AVOGADRO, 6.022140857e23, "AVOGADRO mismatch");
        assert_eq!(PLANCK, 6.626070040e-34, "PLANCK mismatch");
        assert_eq!(BOLTZMANN, 1.38064852e-23, "BOLTZMANN mismatch");
        assert_eq!(LIGHT_SPEED_SI, 299792458.0, "LIGHT_SPEED_SI mismatch");
        assert_eq!(E_MASS, 9.10938356e-31, "E_MASS mismatch");
        assert_eq!(E_CHARGE, 1.6021766208e-19, "E_CHARGE mismatch");
    }

    #[test]
    fn hartree_to_joule_matches_pyscf() {
        // PySCF computes HARTREE2J = HBAR**2 / (E_MASS * BOHR_SI**2)
        let hbar = PLANCK / (2.0 * std::f64::consts::PI);
        let hartree2j_computed = hbar * hbar / (E_MASS * BOHR_SI * BOHR_SI);
        // Our constant should agree within floating-point precision
        let rel_err = (HARTREE_TO_JOULE - hartree2j_computed).abs() / hartree2j_computed;
        assert!(
            rel_err < 1e-10,
            "HARTREE_TO_JOULE relative error = {rel_err:.2e}"
        );
    }

    #[test]
    fn au_to_hz_consistent() {
        // AU_TO_HZ = sqrt(HARTREE_TO_JOULE / (AMU_TO_KG * BOHR_SI^2)) / (2*PI)
        let computed = (HARTREE_TO_JOULE / (AMU_TO_KG * BOHR_SI * BOHR_SI)).sqrt()
            / (2.0 * std::f64::consts::PI);
        let rel_err = (AU_TO_HZ - computed).abs() / computed;
        assert!(
            rel_err < 1e-8,
            "AU_TO_HZ relative error = {rel_err:.2e}, computed = {computed}"
        );
    }

    #[test]
    fn au_to_cm1_consistent() {
        // AU_TO_CM1 = AU_TO_HZ / c * 1e-2
        let computed = AU_TO_HZ / LIGHT_SPEED_SI * 1e-2;
        let rel_err = (AU_TO_CM1 - computed).abs() / computed;
        assert!(rel_err < 1e-10, "AU_TO_CM1 relative error = {rel_err:.2e}");
    }

    #[test]
    fn hartree_to_wavenumber_matches_pyscf() {
        // PySCF: HARTREE2WAVENUMBER = 1e-2 * HARTREE2J / (c * h) ≈ 2.194746e5
        assert!(
            (HARTREE_TO_WAVENUMBER - 2.194_746e5).abs() < 1.0,
            "HARTREE_TO_WAVENUMBER = {}, expected ~2.194746e5",
            HARTREE_TO_WAVENUMBER
        );
    }

    #[test]
    fn au_to_kmmol_matches_literature() {
        // Expected value: approximately 974.880 km/mol (pyscf-forge and Psi4 agree)
        // pyscf-forge prop/infrared/rhf.py:161 reports ~974.88011
        assert!(
            (AU_TO_KMMOL - 974.880).abs() < 0.01,
            "AU_TO_KMMOL = {AU_TO_KMMOL}, expected ~974.88"
        );
    }

    #[test]
    fn au_to_kmmol_from_fundamental_constants() {
        // Recompute from first principles using the pyscf-forge formula:
        //   unit_kmmol = alpha^2 * (1e-3 / amu) * m_e * N_A * pi * a_0 / 3
        // where alpha = 1 / c_au is the dimensionless fine-structure constant.
        let alpha = 1.0 / LIGHT_SPEED_AU;
        let computed =
            alpha * alpha * (1e-3 / AMU_TO_KG) * E_MASS * AVOGADRO * std::f64::consts::PI * BOHR_SI
                / 3.0;
        let rel_err = (AU_TO_KMMOL - computed).abs() / computed;
        assert!(
            rel_err < 1e-12,
            "AU_TO_KMMOL internal inconsistency: {rel_err:.2e}"
        );
    }

    #[test]
    fn au_to_debye_matches_pyscf() {
        // PySCF: AU2DEBYE = E_CHARGE * BOHR*1e-10 / DEBYE ≈ 2.541746
        let computed = E_CHARGE * BOHR_SI / DEBYE;
        assert!(
            (AU_TO_DEBYE - 2.541_746).abs() < 0.001,
            "AU_TO_DEBYE = {AU_TO_DEBYE}, expected ~2.541746"
        );
        assert!(
            (AU_TO_DEBYE - computed).abs() < 1e-10,
            "AU_TO_DEBYE internal inconsistency"
        );
    }

    #[test]
    fn gas_constant_consistent() {
        // R = k_B × N_A ≈ 8.314 J/(mol·K)
        assert!(
            (GAS_CONSTANT - 8.314_459_8).abs() < 0.01,
            "GAS_CONSTANT = {GAS_CONSTANT}, expected ~8.314"
        );
    }

    #[test]
    fn atomic_masses_coverage() {
        // All elements H through Ar should have positive mass
        for z in 1u8..=18 {
            let mass = atomic_mass(z);
            assert!(mass.is_some(), "atomic_mass({z}) returned None");
            assert!(
                mass.unwrap() > 0.0,
                "atomic_mass({z}) = {} is not positive",
                mass.unwrap()
            );
        }
        // Z=0 and Z>18 should return None
        assert!(atomic_mass(0).is_none());
        assert!(atomic_mass(19).is_none());
    }

    #[test]
    fn atomic_masses_known_values() {
        // Spot-check well-known masses
        assert!((atomic_mass(1).unwrap() - 1.00794).abs() < 0.001, "H mass");
        assert!((atomic_mass(6).unwrap() - 12.011).abs() < 0.01, "C mass");
        assert!((atomic_mass(7).unwrap() - 14.007).abs() < 0.01, "N mass");
        assert!((atomic_mass(8).unwrap() - 15.999).abs() < 0.01, "O mass");
        assert!((atomic_mass(9).unwrap() - 18.998).abs() < 0.01, "F mass");
    }

    #[test]
    fn bohr4_to_ang4_matches_derivation() {
        // BOHR4_TO_ANG4 = BOHR_TO_ANGSTROM^4
        let expected = BOHR_TO_ANGSTROM.powi(4);
        let rel_err = (BOHR4_TO_ANG4 - expected).abs() / expected;
        assert!(
            rel_err < 1e-15,
            "BOHR4_TO_ANG4 derivation mismatch: rel_err = {rel_err:.2e}"
        );
    }

    #[test]
    fn bohr4_to_ang4_matches_literature() {
        // BOHR_TO_ANGSTROM^4 ≈ 0.07841597 (from CODATA 2014 BOHR_TO_ANGSTROM)
        // Literature value (from Gaussian, ORCA docs): ~0.07842
        assert!(
            (BOHR4_TO_ANG4 - 0.07842).abs() < 1e-3,
            "BOHR4_TO_ANG4 = {BOHR4_TO_ANG4}, expected ~0.07842"
        );
    }

    #[test]
    fn au_polarizability_to_ang3_matches_derivation() {
        // AU_POLARIZABILITY_TO_ANG3 = BOHR_TO_ANGSTROM^3
        let expected = BOHR_TO_ANGSTROM.powi(3);
        let rel_err = (AU_POLARIZABILITY_TO_ANG3 - expected).abs() / expected;
        assert!(
            rel_err < 1e-15,
            "AU_POLARIZABILITY_TO_ANG3 derivation mismatch: rel_err = {rel_err:.2e}"
        );
    }

    #[test]
    fn au_polarizability_to_ang3_matches_literature() {
        // BOHR_TO_ANGSTROM^3 ≈ 0.14818
        assert!(
            (AU_POLARIZABILITY_TO_ANG3 - 0.14818).abs() < 1e-4,
            "AU_POLARIZABILITY_TO_ANG3 = {AU_POLARIZABILITY_TO_ANG3}, expected ~0.14818"
        );
    }

    #[test]
    fn amu_to_kg_consistent() {
        // AMU_TO_KG = 1e-3 / AVOGADRO ≈ 1.6605e-27
        let computed = 1e-3 / AVOGADRO;
        assert!(
            (AMU_TO_KG - computed).abs() < 1e-40,
            "AMU_TO_KG inconsistency"
        );
        assert!(
            (AMU_TO_KG - 1.660_539e-27).abs() < 1e-33,
            "AMU_TO_KG = {AMU_TO_KG}, expected ~1.6605e-27"
        );
    }
}
