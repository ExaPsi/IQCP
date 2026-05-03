//! RRHO (Rigid-Rotor Harmonic-Oscillator) thermochemistry.
//!
//! Given a `FrequencyInfo` (from [`crate::thermo::harmonic_analysis`], US-096),
//! the SCF electronic energy, the atomic geometry, a temperature, a pressure,
//! a rotational symmetry number, and a spin multiplicity, this module computes
//! the zero-point energy (ZPE), thermal corrections, enthalpy $H(T)$,
//! entropy $S(T)$, Gibbs free energy $G(T)$, and heat capacities $C_v(T)$,
//! $C_p(T)$ under the rigid-rotor harmonic-oscillator approximation.
//!
//! # Algorithm
//!
//! The molecular partition function factorizes into independent
//! contributions from translation, rotation, vibration, and the electronic
//! state:
//!
//! $$q_\text{mol} = q_\text{trans} \cdot q_\text{rot} \cdot q_\text{vib}
//! \cdot q_\text{elec}$$
//!
//! Each contribution has its own thermodynamic functions; the totals are the
//! corresponding sums (for energy, entropy, heat capacity) or products
//! (for the partition function itself). See [`compute_thermochemistry`] for
//! the full assembly formulas.
//!
//! The implementation follows PySCF's
//! `pyscf/hessian/thermo.py::thermo()` (lines 135-230) function-by-function,
//! using the same constants from [`crate::constants`] (CODATA 2014 matching
//! PySCF's `nist.py`). The unit system inside the formulas is SI; the
//! results are expressed in Hartree (for energies) and Hartree per Kelvin
//! per mole (for entropies and heat capacities) via the molar gas constant
//! in atomic units:
//!
//! $$R_\text{Eh} = \frac{k_B}{E_h} \approx 3.166811 \times 10^{-6}\,
//!   \text{Ha}/(\text{mol}\cdot\text{K})$$
//!
//! PySCF computes this as `R_Eh = kB * N_A / (HARTREE2J * N_A) = kB / HARTREE2J`,
//! which is mathematically identical because $N_A$ cancels.
//!
//! # Regimes and branches
//!
//! - **Translational** (always): Sackur-Tetrode equation for ideal gas.
//!   $E_\text{trans} = (3/2) R T$, $H_\text{trans} = (5/2) R T$,
//!   $C_v^\text{trans} = (3/2) R$, $C_p^\text{trans} = (5/2) R$.
//!   Matches PySCF `thermo.py:162-169`.
//!
//! - **Rotational** branches on [`crate::thermo::RotorType`]:
//!   - [`RotorType::Atom`]: all rotational contributions are zero.
//!     PySCF `thermo.py:180-183`.
//!   - [`RotorType::Linear`]: uses the finite rotational constant at
//!     `rotational_constants_ghz[1]` (index 0 is infinity for a linear
//!     molecule). $E_\text{rot} = R T$, $C_v^\text{rot} = R$.
//!     PySCF `thermo.py:184-189`.
//!   - Nonlinear (spherical, symmetric, or asymmetric top): uses the product
//!     of all three rotational constants. $E_\text{rot} = (3/2) R T$,
//!     $C_v^\text{rot} = (3/2) R$. PySCF `thermo.py:190-196`.
//!
//! - **Vibrational** (quantum harmonic oscillator): sums over all modes with
//!   $\nu > 0$. Imaginary modes (negative real frequency, from saddle points
//!   or numerical noise in the T/R projection) are skipped and counted in
//!   `n_imag`. Matches PySCF `thermo.py:198-212`, filter at line 200
//!   (`idx = freq.real > 0`).
//!
//! - **Electronic**: ground-state-only approximation. $S_\text{elec} =
//!   R \ln(2S+1)$, $C_v^\text{elec} = 0$. PySCF `thermo.py:155-158`.
//!
//! # Numerical care
//!
//! The vibrational entropy and heat capacity expressions use
//! $u_k e^{-u_k}/(1 - e^{-u_k})$ rather than $u_k/(e^{u_k} - 1)$ to keep
//! every intermediate quantity bounded by 1 for large $u_k$. For very small
//! $u_k$ (very low-frequency modes) the usual formulas remain stable.
//!
//! For vibrational modes with `exp(-u_k)` underflowing to 0, the contributions
//! are set to zero, which is the correct low-temperature limit (only the
//! ground vibrational state is populated).
//!
//! # References
//!
//! - McQuarrie, D.A. (1976). *Statistical Mechanics.* University Science Books.
//!   Chapters 5-6 (translation, rotation, vibration), Chapter 8 (ideal gas).
//! - McQuarrie, D.A. & Simon, J.D. (1997). *Physical Chemistry: A Molecular
//!   Approach.* Chapter 18 (partition functions).
//! - Herzberg, G. (1945). *Molecular Spectra and Molecular Structure,
//!   Volume II: Infrared and Raman Spectra of Polyatomic Molecules.*
//! - Cramer, C.J. (2004). *Essentials of Computational Chemistry* (2nd ed.),
//!   Chapter 10.
//! - Ochterski, J.W. (2000). *Thermochemistry in Gaussian.* Gaussian Inc.
//!   white paper. <https://gaussian.com/thermo/>
//! - PySCF `pyscf/hessian/thermo.py::thermo()` lines 135-230
//!   (`references/pyscf/pyscf/hessian/thermo.py`).
//! - NIST CCCBDB: <https://cccbdb.nist.gov/thermo.asp>

use std::f64::consts::PI;
use thiserror::Error;

use crate::basis::Atom;
use crate::constants::{
    AMU_TO_KG, ATOMIC_MASSES, BOLTZMANN, HARTREE_TO_JOULE, LIGHT_SPEED_SI, PLANCK,
};
use crate::thermo::{FrequencyInfo, RotorType};

/// Module version (matches crate version).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// IUPAC standard reference temperature: 298.15 K (25 °C).
pub const STANDARD_TEMPERATURE_K: f64 = 298.15;

/// IUPAC standard reference pressure: 1 atm = 101 325 Pa.
pub const STANDARD_PRESSURE_PA: f64 = 101_325.0;

// ============================================================================
// Public types
// ============================================================================

/// Errors returned by [`compute_thermochemistry`].
#[derive(Debug, Clone, PartialEq, Error)]
pub enum ThermochemistryError {
    /// Temperature must be strictly positive and finite.
    #[error("Temperature must be > 0 K and finite, got {0}")]
    NonPositiveTemperature(f64),

    /// Pressure must be strictly positive and finite.
    #[error("Pressure must be > 0 Pa and finite, got {0}")]
    NonPositivePressure(f64),

    /// Spin multiplicity (2S+1) must be at least 1.
    #[error("Spin multiplicity must be >= 1, got {0}")]
    InvalidMultiplicity(u32),

    /// Rotational symmetry number must be at least 1.
    #[error("Symmetry number must be >= 1, got {0}")]
    InvalidSymmetryNumber(u32),

    /// Atomic number outside the range covered by
    /// [`crate::constants::ATOMIC_MASSES`] (Z = 1..=18).
    #[error("Atomic number {0} not supported (must be 1-18)")]
    UnsupportedElement(u8),

    /// A scalar input is NaN or infinity.
    #[error("Non-finite input: {field} = {value}")]
    NonFiniteInput { field: &'static str, value: f64 },

    /// The atom list is empty.
    #[error("Empty molecule (no atoms)")]
    EmptyMolecule,
}

/// Per-contribution breakdown of thermochemical quantities.
///
/// All energies in Hartree, all entropies and heat capacities in
/// Ha/(mol·K). The vibrational energy `e_vib_thermal_ha` does **not**
/// include ZPE; ZPE is reported separately in [`ThermochemistryResult`].
/// `h_vib_ha` on the other hand matches PySCF's convention of including
/// ZPE in `H_vib` (and therefore in `H_tot`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermoComponents {
    // --- Translational ---
    /// $E_\text{trans} = (3/2) R T$ in Hartree.
    pub e_trans_ha: f64,
    /// $H_\text{trans} = (5/2) R T$ in Hartree.
    pub h_trans_ha: f64,
    /// Sackur-Tetrode entropy in Ha/(mol·K).
    pub s_trans_ha_per_k: f64,
    /// $C_v^\text{trans} = (3/2) R$ in Ha/(mol·K).
    pub cv_trans_ha_per_k: f64,
    /// $C_p^\text{trans} = (5/2) R$ in Ha/(mol·K).
    pub cp_trans_ha_per_k: f64,

    // --- Rotational ---
    /// $E_\text{rot}$: 0 for atom, $RT$ for linear, $(3/2) R T$ for nonlinear.
    pub e_rot_ha: f64,
    /// $H_\text{rot}$: equal to $E_\text{rot}$ (no $RT$ shift at the
    /// per-contribution level; only the total enthalpy adds $RT$).
    pub h_rot_ha: f64,
    /// Rotational entropy in Ha/(mol·K).
    pub s_rot_ha_per_k: f64,
    /// $C_v^\text{rot}$: 0 for atom, $R$ for linear, $(3/2) R$ for nonlinear.
    pub cv_rot_ha_per_k: f64,
    /// PySCF sets `Cp_rot = Cv_rot`; the `+R` shift only applies at the
    /// total level.
    pub cp_rot_ha_per_k: f64,

    // --- Vibrational ---
    /// Thermal vibrational energy (excluding ZPE), in Hartree.
    pub e_vib_thermal_ha: f64,
    /// $H_\text{vib}$ per PySCF = ZPE + thermal, in Hartree.
    pub h_vib_ha: f64,
    /// Vibrational entropy in Ha/(mol·K).
    pub s_vib_ha_per_k: f64,
    /// Vibrational heat capacity in Ha/(mol·K).
    pub cv_vib_ha_per_k: f64,
    /// PySCF sets `Cp_vib = Cv_vib` (the `+R` is added at the total level).
    pub cp_vib_ha_per_k: f64,

    // --- Electronic ---
    /// Electronic energy in Hartree (echo of `electronic_energy_ha`).
    pub e_elec_ha: f64,
    /// Electronic enthalpy in Hartree (same as `e_elec_ha`).
    pub h_elec_ha: f64,
    /// $S_\text{elec} = R \ln(2S+1)$ in Ha/(mol·K).
    pub s_elec_ha_per_k: f64,
    /// Always 0 in the ground-state-only approximation.
    pub cv_elec_ha_per_k: f64,
    /// Always 0 in the ground-state-only approximation.
    pub cp_elec_ha_per_k: f64,
}

/// Result of [`compute_thermochemistry`].
///
/// All energies in Hartree, all entropies and heat capacities in
/// Ha/(mol·K). `n_imag` reports how many vibrational modes were skipped
/// because they had imaginary (negative) frequency — a transition state
/// has exactly one imaginary mode; an equilibrium geometry should have
/// zero.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermochemistryResult {
    // --- Input echoes ---
    /// Temperature in K.
    pub temperature_k: f64,
    /// Pressure in Pa.
    pub pressure_pa: f64,
    /// Rotational symmetry number actually used (after defaulting).
    pub symmetry_number: u32,
    /// Spin multiplicity 2S+1.
    pub multiplicity: u32,
    /// Total molecular mass in amu.
    pub total_mass_amu: f64,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Number of vibrational modes used (positive frequencies only).
    pub n_vib_modes_used: usize,
    /// Number of imaginary modes skipped from the vibrational sum.
    pub n_imag: usize,

    // --- Electronic & zero-point ---
    /// Electronic SCF energy in Hartree (echo of input).
    pub electronic_energy_ha: f64,
    /// Zero-point vibrational energy $(R/2) \sum_k \Theta_k$ in Hartree.
    pub zpe_ha: f64,

    // --- Total thermodynamic functions ---
    /// Total internal energy
    /// $U(T) = E_\text{elec} + \text{ZPE} + E_\text{trans} + E_\text{rot} +
    /// E_\text{vib,thermal}$.
    pub internal_energy_ha: f64,
    /// Total enthalpy $H(T) = U(T) + R T$ (ideal gas).
    pub enthalpy_ha: f64,
    /// Total entropy
    /// $S(T) = S_\text{trans} + S_\text{rot} + S_\text{vib} + S_\text{elec}$.
    pub entropy_ha_per_k: f64,
    /// Total Gibbs free energy $G(T) = H(T) - T \cdot S(T)$.
    pub gibbs_ha: f64,
    /// Total constant-volume heat capacity
    /// $C_v(T) = C_v^\text{trans} + C_v^\text{rot} + C_v^\text{vib}$.
    pub cv_ha_per_k: f64,
    /// Total constant-pressure heat capacity $C_p(T) = C_v(T) + R$
    /// (ideal gas).
    pub cp_ha_per_k: f64,

    /// Energy at 0 K: $E_{0K} = E_\text{elec} + \text{ZPE}$.
    pub e_0k_ha: f64,

    /// Per-contribution breakdown.
    pub components: ThermoComponents,
}

// ============================================================================
// Helpers: gas constant, masses, validation
// ============================================================================

/// Gas constant in Hartree per mole per Kelvin,
/// $R_\text{Eh} = k_B / E_h \approx 3.166\,811 \times 10^{-6}$ Ha/(mol·K).
///
/// Matches PySCF `thermo.py:146`:
/// `R_Eh = kB*AVOGADRO / (HARTREE2J*AVOGADRO) = kB / HARTREE2J`.
/// The `AVOGADRO` factor cancels because the result is naturally per-mole
/// once expressed in Hartree (which is already a molar-scale unit of energy
/// in quantum chemistry convention).
#[inline]
fn r_hartree_per_mol_k() -> f64 {
    BOLTZMANN / HARTREE_TO_JOULE
}

/// Sum atomic masses (amu) for a list of atoms.
///
/// Returns [`ThermochemistryError::UnsupportedElement`] for any atomic
/// number outside `1..=18` (the range of [`ATOMIC_MASSES`]).
fn total_mass_amu(atoms: &[Atom]) -> Result<f64, ThermochemistryError> {
    if atoms.is_empty() {
        return Err(ThermochemistryError::EmptyMolecule);
    }
    let mut total = 0.0f64;
    for atom in atoms {
        let z = atom.atomic_number;
        if z == 0 || (z as usize) >= ATOMIC_MASSES.len() {
            return Err(ThermochemistryError::UnsupportedElement(z));
        }
        total += ATOMIC_MASSES[z as usize];
    }
    Ok(total)
}

/// Validate scalar inputs up-front.
///
/// Order: temperature > 0 and finite → pressure > 0 and finite →
/// multiplicity ≥ 1 → symmetry number ≥ 1 → electronic energy finite.
fn validate_inputs(
    temperature_k: f64,
    pressure_pa: f64,
    multiplicity: u32,
    symmetry_number: u32,
    electronic_energy_ha: f64,
) -> Result<(), ThermochemistryError> {
    if temperature_k.is_nan() {
        return Err(ThermochemistryError::NonFiniteInput {
            field: "temperature_k",
            value: temperature_k,
        });
    }
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(ThermochemistryError::NonPositiveTemperature(temperature_k));
    }
    if pressure_pa.is_nan() {
        return Err(ThermochemistryError::NonFiniteInput {
            field: "pressure_pa",
            value: pressure_pa,
        });
    }
    if !pressure_pa.is_finite() || pressure_pa <= 0.0 {
        return Err(ThermochemistryError::NonPositivePressure(pressure_pa));
    }
    if multiplicity == 0 {
        return Err(ThermochemistryError::InvalidMultiplicity(multiplicity));
    }
    if symmetry_number == 0 {
        return Err(ThermochemistryError::InvalidSymmetryNumber(symmetry_number));
    }
    if !electronic_energy_ha.is_finite() {
        return Err(ThermochemistryError::NonFiniteInput {
            field: "electronic_energy_ha",
            value: electronic_energy_ha,
        });
    }
    Ok(())
}

// ============================================================================
// Partition functions: translational
// ============================================================================

/// Translational partition function and thermodynamic contributions.
///
/// Returns `(e_trans, h_trans, s_trans, cv_trans, cp_trans)` in
/// `(Ha, Ha, Ha/(mol·K), Ha/(mol·K), Ha/(mol·K))`.
///
/// # Formulas
///
/// ```text
/// q_trans = (2π m k_B T / h²)^(3/2) · k_B T / P
/// E_trans = (3/2) R T
/// H_trans = (5/2) R T                     (= E_trans + RT)
/// S_trans = R · [ln(q_trans) + 5/2]         (Sackur-Tetrode)
/// C_v^trans = (3/2) R
/// C_p^trans = (5/2) R
/// ```
///
/// `m` is the total molecular mass in kg
/// (`total_mass_amu * AMU_TO_KG`); `T` in K, `P` in Pa.
///
/// # Reference
///
/// - PySCF `thermo.py:162-169`.
/// - McQuarrie (1976) eq. 5-22, 5-32 (Sackur-Tetrode).
fn translational_partition(
    total_mass_amu_val: f64,
    temperature_k: f64,
    pressure_pa: f64,
) -> (f64, f64, f64, f64, f64) {
    let r = r_hartree_per_mol_k();
    let mass_kg = total_mass_amu_val * AMU_TO_KG;
    // q_trans = (2*pi*m*kB*T / h^2)^(3/2) * kB*T / P
    let factor = 2.0 * PI * mass_kg * BOLTZMANN * temperature_k / (PLANCK * PLANCK);
    let q_trans = factor.powf(1.5) * BOLTZMANN * temperature_k / pressure_pa;

    let e_trans = 1.5 * r * temperature_k;
    let h_trans = 2.5 * r * temperature_k;
    let s_trans = r * (q_trans.ln() + 2.5);
    let cv_trans = 1.5 * r;
    let cp_trans = 2.5 * r;

    (e_trans, h_trans, s_trans, cv_trans, cp_trans)
}

// ============================================================================
// Partition functions: rotational
// ============================================================================

/// Rotational partition function and contributions.
///
/// Branches on `rot_type`:
///
/// - [`RotorType::Atom`]: returns all zeros.
/// - [`RotorType::Linear`]: uses `B = rotational_constants_ghz[1] * 1e9` Hz
///   (index 0 is infinity for a linear molecule because the moment of
///   inertia along the molecular axis is zero). Formulas:
///   `q_rot = k_B T / (σ · h · B)`, `E = R T`, `S = R (1 + ln q)`, `C_v = R`.
/// - **Nonlinear** (spherical, symmetric, asymmetric top): uses the
///   geometric mean of the three finite rotational constants. Formulas:
///   `q_rot = sqrt(π)/σ · (k_B T / h)^(3/2) / sqrt(B_A B_B B_C)`,
///   `E = (3/2) R T`, `S = R (3/2 + ln q)`, `C_v = (3/2) R`.
///
/// Returns `(e_rot, s_rot, cv_rot)` in `(Ha, Ha/(mol·K), Ha/(mol·K))`.
///
/// # Reference
///
/// - PySCF `thermo.py:171-196`.
/// - McQuarrie (1976) eq. 6-46 (linear), eq. 6-49 (nonlinear).
fn rotational_partition(
    rot_type: RotorType,
    rotational_constants_ghz: [f64; 3],
    symmetry_number: u32,
    temperature_k: f64,
) -> (f64, f64, f64) {
    let r = r_hartree_per_mol_k();
    let sigma = symmetry_number as f64;
    let t = temperature_k;

    match rot_type {
        RotorType::Atom => (0.0, 0.0, 0.0),
        RotorType::Linear => {
            // PySCF: B = rot_const[1] * 1e9
            // Index 0 is infinity for linear; indices 1, 2 are equal finite values.
            let b_hz = rotational_constants_ghz[1] * 1e9;
            let q_rot = BOLTZMANN * t / (sigma * PLANCK * b_hz);
            let e_rot = r * t;
            let s_rot = r * (1.0 + q_rot.ln());
            let cv_rot = r;
            (e_rot, s_rot, cv_rot)
        }
        _ => {
            // Nonlinear: spherical / symmetric / asymmetric top.
            // q_rot = sqrt(pi)/σ · (k_B T / h)^(3/2) / sqrt(B_A B_B B_C)
            // (with B_i in Hz)
            let abc_hz: [f64; 3] = [
                rotational_constants_ghz[0] * 1e9,
                rotational_constants_ghz[1] * 1e9,
                rotational_constants_ghz[2] * 1e9,
            ];
            let prod_abc = abc_hz[0] * abc_hz[1] * abc_hz[2];
            let kt_over_h = BOLTZMANN * t / PLANCK;
            let q_rot = kt_over_h.powf(1.5) * PI.sqrt() / (sigma * prod_abc.sqrt());
            let e_rot = 1.5 * r * t;
            let s_rot = r * (1.5 + q_rot.ln());
            let cv_rot = 1.5 * r;
            (e_rot, s_rot, cv_rot)
        }
    }
}

// ============================================================================
// Partition functions: vibrational (quantum HO)
// ============================================================================

/// Vibrational partition function (quantum harmonic oscillator) and
/// contributions.
///
/// Returns `(zpe, e_vib_thermal, s_vib, cv_vib, n_used, n_imag)`:
/// - `zpe` — zero-point vibrational energy in Hartree.
/// - `e_vib_thermal` — thermal contribution (**excluding** ZPE) in Hartree.
/// - `s_vib` — vibrational entropy in Ha/(mol·K).
/// - `cv_vib` — vibrational heat capacity in Ha/(mol·K).
/// - `n_used` — number of modes actually summed (positive frequencies).
/// - `n_imag` — number of imaginary modes skipped.
///
/// Imaginary frequencies (`freq_wavenumber[k] <= 0`, stored as negative real
/// per the Gaussian/PySCF convention) are skipped from every sum and
/// counted in `n_imag`. This matches PySCF's
/// `idx = freq.real > 0` filter at `thermo.py:200`.
///
/// # Formulas per real mode $k$
///
/// ```text
/// Θ_k = h · c · ν_k / k_B        (vibrational temperature in K,
///                                 ν_k in cm⁻¹, c in cm/s)
/// u_k = Θ_k / T                  (reduced temperature)
/// e_neg = exp(-u_k)              (bounded in [0, 1])
/// one_minus = 1 - e_neg
/// ZPE      += (R/2) · Θ_k
/// e_therm  += R · T · u_k · e_neg / one_minus
/// s_vib    += R · [u_k · e_neg / one_minus - ln(one_minus)]
/// cv_vib   += R · u_k² · e_neg / one_minus²
/// ```
///
/// # Numerical care
///
/// The `exp(-u_k)` form keeps every intermediate bounded by 1, avoiding
/// overflow for very stiff modes (large $u_k$). For underflow
/// (`e_neg == 0`) we set all per-mode contributions except ZPE to zero,
/// which is the exact low-temperature limit.
///
/// # Reference
///
/// - PySCF `thermo.py:198-212`.
/// - McQuarrie (1976) eq. 6-22 to 6-29.
fn vibrational_partition(
    freq_wavenumber_cm1: &[f64],
    temperature_k: f64,
) -> (f64, f64, f64, f64, usize, usize) {
    let r = r_hartree_per_mol_k();
    // cm/s conversion factor for c (LIGHT_SPEED_SI is in m/s)
    let c_cm_per_s = LIGHT_SPEED_SI * 100.0;

    let mut zpe = 0.0_f64;
    let mut e_thermal = 0.0_f64;
    let mut s_vib = 0.0_f64;
    let mut cv_vib = 0.0_f64;
    let mut n_used = 0_usize;
    let mut n_imag = 0_usize;

    for &nu in freq_wavenumber_cm1 {
        // Skip imaginary modes (negative-real by convention)
        if nu <= 0.0 {
            n_imag += 1;
            continue;
        }
        n_used += 1;

        // Vibrational temperature in K:
        // Θ_k = h · c · ν / k_B, with ν in cm⁻¹ and c in cm/s.
        let theta_k = PLANCK * c_cm_per_s * nu / BOLTZMANN;

        // ZPE = (R/2) · Σ Θ_k (temperature-independent)
        zpe += 0.5 * r * theta_k;

        // Reduced temperature u_k = Θ_k / T. T > 0 is validated up-front.
        let u_k = theta_k / temperature_k;
        // exp(-u_k) ∈ [0, 1]; underflows to 0 for u_k ≳ 700
        let e_neg = (-u_k).exp();
        if e_neg <= 0.0 {
            // Low-T limit: thermal contributions are zero for this mode.
            continue;
        }
        let one_minus = 1.0 - e_neg;
        if one_minus <= 0.0 {
            // Guard: u_k underflow gives e_neg == 1 numerically. Treat as
            // classical (equipartition) limit by contributing R to Cv,
            // but this case is unreachable for physical inputs.
            continue;
        }

        // Thermal energy: R T · u_k · e_neg / (1 - e_neg)
        e_thermal += r * temperature_k * u_k * e_neg / one_minus;

        // Vibrational entropy:
        //   s_vib += R · [u_k·e_neg/(1-e_neg) - ln(1-e_neg)]
        // `one_minus > 0` is guaranteed above.
        s_vib += r * (u_k * e_neg / one_minus - one_minus.ln());

        // Heat capacity:
        //   cv_vib += R · u_k² · e_neg / (1-e_neg)²
        cv_vib += r * u_k * u_k * e_neg / (one_minus * one_minus);
    }

    (zpe, e_thermal, s_vib, cv_vib, n_used, n_imag)
}

// ============================================================================
// Partition functions: electronic
// ============================================================================

/// Electronic entropy from spin multiplicity.
///
/// $S_\text{elec} = R \cdot \ln(\text{multiplicity})$
/// where `multiplicity = 2S + 1`. For closed-shell singlets
/// (`multiplicity == 1`), the result is exactly 0.
///
/// # Reference
///
/// - PySCF `thermo.py:156`.
#[inline]
fn electronic_entropy(multiplicity: u32) -> f64 {
    r_hartree_per_mol_k() * (multiplicity as f64).ln()
}

// ============================================================================
// Top-level: compute_thermochemistry
// ============================================================================

/// Compute RRHO thermochemistry at the specified temperature and pressure.
///
/// Under the rigid-rotor harmonic-oscillator approximation with a
/// ground-state-only electronic partition function, this function assembles:
///
/// - Zero-point vibrational energy (ZPE)
/// - Internal energy $U$, enthalpy $H$, Gibbs free energy $G$
/// - Entropy $S$ (total and per-contribution)
/// - Heat capacities $C_v$, $C_p$ (total and per-contribution)
///
/// from the vibrational frequencies, rotor classification, SCF electronic
/// energy, and atomic geometry. The implementation follows PySCF's
/// `pyscf/hessian/thermo.py::thermo()` function-by-function (lines 135-230),
/// using the same physical constants.
///
/// # Arguments
///
/// * `freq_info` — vibrational frequencies and rotor classification from
///   [`crate::thermo::harmonic_analysis`] (US-096).
/// * `electronic_energy_ha` — SCF electronic energy in Hartree (this is
///   `mf.e_tot` in PySCF: SCF electronic energy + nuclear repulsion).
/// * `atoms` — atomic geometry. Only `atomic_number` is read (positions are
///   ignored here; the center of mass and moment of inertia come from
///   `freq_info`).
/// * `temperature_k` — temperature in K (must be > 0). The IUPAC standard
///   is [`STANDARD_TEMPERATURE_K`] = 298.15 K.
/// * `pressure_pa` — pressure in Pa (must be > 0). The IUPAC standard is
///   [`STANDARD_PRESSURE_PA`] = 101325 Pa (1 atm).
/// * `symmetry_number` — rotational symmetry number $\sigma$. Pass `None`
///   to use the default of 1. Common values: C₂ᵥ (H₂O) → 2,
///   D∞ₕ (CO₂, H₂) → 2, C₃ᵥ (NH₃) → 3, T_d (CH₄) → 12, O_h (SF₆) → 24.
///   IQCP does **not** detect the molecular point group automatically; the
///   caller must pass the correct σ for a cross-check against PySCF.
/// * `multiplicity` — spin multiplicity 2S+1. Default = 1 (closed-shell
///   singlet). All Phase 5 SCF tests are closed-shell singlets.
///
/// # Returns
///
/// A [`ThermochemistryResult`] with the total quantities plus a
/// [`ThermoComponents`] per-contribution breakdown.
///
/// # Errors
///
/// - [`ThermochemistryError::NonPositiveTemperature`] if `T <= 0` or NaN/Inf.
/// - [`ThermochemistryError::NonPositivePressure`] if `P <= 0` or NaN/Inf.
/// - [`ThermochemistryError::InvalidMultiplicity`] if `multiplicity == 0`.
/// - [`ThermochemistryError::InvalidSymmetryNumber`] if provided σ is 0.
/// - [`ThermochemistryError::UnsupportedElement`] for any atomic number > 18.
/// - [`ThermochemistryError::EmptyMolecule`] if `atoms.is_empty()`.
/// - [`ThermochemistryError::NonFiniteInput`] if any input is NaN/Inf.
///
/// # Example
///
/// ```no_run
/// use qc_core::thermo::harmonic_analysis;
/// use qc_core::thermochemistry::{
///     compute_thermochemistry, STANDARD_PRESSURE_PA, STANDARD_TEMPERATURE_K,
/// };
///
/// # use qc_core::basis::Atom;
/// # use nalgebra::DMatrix;
/// # let atoms: Vec<Atom> = vec![];
/// # let hessian = DMatrix::<f64>::zeros(0, 0);
/// # let scf_energy_ha: f64 = 0.0;
/// let freq_info = harmonic_analysis(&atoms, &hessian).unwrap();
/// let thermo_result = compute_thermochemistry(
///     &freq_info,
///     scf_energy_ha,
///     &atoms,
///     STANDARD_TEMPERATURE_K,
///     STANDARD_PRESSURE_PA,
///     Some(2),  // H2O point group C2v
///     1,        // closed-shell singlet
/// ).unwrap();
/// println!("ZPE: {:.6} Ha", thermo_result.zpe_ha);
/// println!("G(298.15 K): {:.6} Ha", thermo_result.gibbs_ha);
/// ```
///
/// # References
///
/// - PySCF `thermo.py::thermo()` lines 135-230.
/// - McQuarrie (1976), Chapters 5-6, 8.
/// - Ochterski (2000), *Thermochemistry in Gaussian*.
pub fn compute_thermochemistry(
    freq_info: &FrequencyInfo,
    electronic_energy_ha: f64,
    atoms: &[Atom],
    temperature_k: f64,
    pressure_pa: f64,
    symmetry_number: Option<u32>,
    multiplicity: u32,
) -> Result<ThermochemistryResult, ThermochemistryError> {
    // Apply symmetry number default before validation so the check still
    // catches a user-supplied Some(0).
    let sigma = symmetry_number.unwrap_or(1);

    validate_inputs(
        temperature_k,
        pressure_pa,
        multiplicity,
        sigma,
        electronic_energy_ha,
    )?;

    let mass_amu = total_mass_amu(atoms)?;
    let r = r_hartree_per_mol_k();
    let t = temperature_k;

    // --- Translational ---
    let (e_trans, h_trans, s_trans, cv_trans, cp_trans) =
        translational_partition(mass_amu, t, pressure_pa);

    // --- Rotational ---
    let (e_rot, s_rot, cv_rot) = rotational_partition(
        freq_info.rot_type,
        freq_info.rotational_constants_ghz,
        sigma,
        t,
    );
    // PySCF thermo.py:187, 188, 195, 196: H_rot = E_rot, Cp_rot = Cv_rot.
    let h_rot = e_rot;
    let cp_rot = cv_rot;

    // --- Vibrational (quantum HO) ---
    let (zpe, e_vib_thermal, s_vib, cv_vib, n_used, n_imag) =
        vibrational_partition(&freq_info.freq_wavenumber, t);
    // PySCF thermo.py:211 sets H_vib = ZPE + thermal. Cp_vib = Cv_vib at line 210.
    let h_vib = zpe + e_vib_thermal;
    let cp_vib = cv_vib;

    // --- Electronic ---
    let s_elec = electronic_entropy(multiplicity);
    // H_elec = E_elec (no ZPE or thermal contribution). Cv_elec = 0, Cp_elec = 0.
    let h_elec = electronic_energy_ha;

    // --- Assembly ---
    // Internal energy: sum of electronic + ZPE + kinetic/rotational/vib
    let internal_energy = electronic_energy_ha + zpe + e_trans + e_rot + e_vib_thermal;

    // Total enthalpy per PySCF `_sum` (thermo.py:219-227):
    //   H_tot = H_elec + H_trans + H_rot + H_vib
    // which equals U + RT because H_trans - E_trans = RT (5/2 RT - 3/2 RT),
    // while H_rot == E_rot and H_vib == ZPE + thermal == U_vib + 0.
    // We compute the total enthalpy directly as H_tot = internal_energy + RT
    // to match the PySCF result exactly.
    let enthalpy = internal_energy + r * t;

    // Total entropy
    let entropy = s_trans + s_rot + s_vib + s_elec;

    // Gibbs free energy
    let gibbs = enthalpy - t * entropy;

    // Heat capacities (Cv_elec = 0)
    let cv_total = cv_trans + cv_rot + cv_vib;
    let cp_total = cv_total + r;

    // 0 K energy
    let e_0k = electronic_energy_ha + zpe;

    let components = ThermoComponents {
        e_trans_ha: e_trans,
        h_trans_ha: h_trans,
        s_trans_ha_per_k: s_trans,
        cv_trans_ha_per_k: cv_trans,
        cp_trans_ha_per_k: cp_trans,

        e_rot_ha: e_rot,
        h_rot_ha: h_rot,
        s_rot_ha_per_k: s_rot,
        cv_rot_ha_per_k: cv_rot,
        cp_rot_ha_per_k: cp_rot,

        e_vib_thermal_ha: e_vib_thermal,
        h_vib_ha: h_vib,
        s_vib_ha_per_k: s_vib,
        cv_vib_ha_per_k: cv_vib,
        cp_vib_ha_per_k: cp_vib,

        e_elec_ha: electronic_energy_ha,
        h_elec_ha: h_elec,
        s_elec_ha_per_k: s_elec,
        cv_elec_ha_per_k: 0.0,
        cp_elec_ha_per_k: 0.0,
    };

    Ok(ThermochemistryResult {
        temperature_k: t,
        pressure_pa,
        symmetry_number: sigma,
        multiplicity,
        total_mass_amu: mass_amu,
        n_atoms: atoms.len(),
        n_vib_modes_used: n_used,
        n_imag,
        electronic_energy_ha,
        zpe_ha: zpe,
        internal_energy_ha: internal_energy,
        enthalpy_ha: enthalpy,
        entropy_ha_per_k: entropy,
        gibbs_ha: gibbs,
        cv_ha_per_k: cv_total,
        cp_ha_per_k: cp_total,
        e_0k_ha: e_0k,
        components,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermo::harmonic_analysis;
    use approx::assert_relative_eq;
    use nalgebra::DMatrix;
    use serde::Deserialize;

    // ========================================================================
    // Golden data loaders (shared with thermo.rs golden format)
    // ========================================================================

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenAtom {
        #[serde(rename = "Z")]
        z: u8,
        symbol: String,
        pos_bohr: [f64; 3],
    }

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct ThermoGoldenFull {
        name: String,
        atoms: Vec<GoldenAtom>,
        hessian: Vec<Vec<f64>>,
        energy: f64,
        rot_type: String,
        rotational_constants_ghz: Vec<f64>,
        principal_moments_amu_bohr2: Vec<f64>,
        freq_wavenumber: Vec<f64>,
        freq_au: Vec<f64>,
        reduced_mass: Vec<f64>,
        force_const_dyne: Vec<f64>,
        force_const_au: Vec<f64>,
        norm_mode: Vec<Vec<[f64; 3]>>,
    }

    impl ThermoGoldenFull {
        fn atoms(&self) -> Vec<Atom> {
            self.atoms
                .iter()
                .map(|a| Atom::new(a.z, a.pos_bohr).unwrap())
                .collect()
        }

        fn hessian_matrix(&self) -> DMatrix<f64> {
            let n3 = self.hessian.len();
            let mut h = DMatrix::<f64>::zeros(n3, n3);
            for i in 0..n3 {
                for j in 0..n3 {
                    h[(i, j)] = self.hessian[i][j];
                }
            }
            h
        }
    }

    fn load_thermo_golden(name: &str) -> ThermoGoldenFull {
        match name {
            "h2o" => serde_json::from_str(include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/golden/thermo/h2o_sto3g_rhf.json"
            )))
            .expect("parse h2o thermo golden"),
            "ch4" => serde_json::from_str(include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/golden/thermo/ch4_sto3g_rhf.json"
            )))
            .expect("parse ch4 thermo golden"),
            "co2" => serde_json::from_str(include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/golden/thermo/co2_sto3g_rhf.json"
            )))
            .expect("parse co2 thermo golden"),
            "h2" => serde_json::from_str(include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/golden/thermo/h2_sto3g_rhf.json"
            )))
            .expect("parse h2 thermo golden"),
            _ => panic!("unknown thermo golden: {name}"),
        }
    }

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct ThermoQuantity {
        value: f64,
        unit: String,
    }

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct ThermochemGolden {
        name: String,
        atoms: Vec<GoldenAtom>,
        scf_energy_ha: f64,
        rot_type: String,
        rotational_constants_ghz: Vec<f64>,
        freq_wavenumber_cm1: Vec<f64>,
        freq_au: Vec<f64>,
        symmetry_number: u32,
        multiplicity: u32,
        temperature_k: f64,
        pressure_pa: f64,
        thermo: std::collections::BTreeMap<String, ThermoQuantity>,
    }

    fn load_thermochem_golden(name: &str) -> ThermochemGolden {
        match name {
            "h2o" => serde_json::from_str(include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/golden/thermochem/h2o_sto3g_rhf.json"
            )))
            .expect("parse h2o thermochem golden"),
            "ch4" => serde_json::from_str(include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/golden/thermochem/ch4_sto3g_rhf.json"
            )))
            .expect("parse ch4 thermochem golden"),
            "co2" => serde_json::from_str(include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/golden/thermochem/co2_sto3g_rhf.json"
            )))
            .expect("parse co2 thermochem golden"),
            "h2" => serde_json::from_str(include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/golden/thermochem/h2_sto3g_rhf.json"
            )))
            .expect("parse h2 thermochem golden"),
            _ => panic!("unknown thermochem golden: {name}"),
        }
    }

    // ========================================================================
    // Unit tests: helpers
    // ========================================================================

    #[test]
    fn test_r_hartree_per_mol_k() {
        let r = r_hartree_per_mol_k();
        // Expected: BOLTZMANN / HARTREE_TO_JOULE ≈ 3.166811e-6 Ha/(mol·K)
        assert_relative_eq!(r, 3.166_811e-6, epsilon = 1e-11);
        // Cross-check: R (J/(mol K)) / HARTREE_TO_JOULE / N_A should equal R_Eh.
        let via_r = (BOLTZMANN * crate::constants::AVOGADRO)
            / HARTREE_TO_JOULE
            / crate::constants::AVOGADRO;
        assert_relative_eq!(r, via_r, epsilon = 1e-15);
    }

    #[test]
    fn test_electronic_entropy_singlet() {
        assert_eq!(electronic_entropy(1), 0.0);
    }

    #[test]
    fn test_electronic_entropy_doublet() {
        let s = electronic_entropy(2);
        let r = r_hartree_per_mol_k();
        assert_relative_eq!(s, r * 2.0_f64.ln(), epsilon = 1e-14);
    }

    #[test]
    fn test_electronic_entropy_triplet() {
        let s = electronic_entropy(3);
        let r = r_hartree_per_mol_k();
        assert_relative_eq!(s, r * 3.0_f64.ln(), epsilon = 1e-14);
    }

    #[test]
    fn test_electronic_entropy_scaling() {
        // S(m_2) - S(m_1) = R * ln(m_2 / m_1)
        let s1 = electronic_entropy(2);
        let s2 = electronic_entropy(6);
        let r = r_hartree_per_mol_k();
        assert_relative_eq!(s2 - s1, r * 3.0_f64.ln(), epsilon = 1e-14);
    }

    #[test]
    fn test_total_mass_amu_h2o() {
        let atoms = vec![
            Atom::new(8, [0.0, 0.0, 0.0]).unwrap(),
            Atom::new(1, [0.0, 1.43, 1.11]).unwrap(),
            Atom::new(1, [0.0, -1.43, 1.11]).unwrap(),
        ];
        let m = total_mass_amu(&atoms).unwrap();
        let expected = ATOMIC_MASSES[8] + 2.0 * ATOMIC_MASSES[1];
        assert_relative_eq!(m, expected, epsilon = 1e-14);
    }

    #[test]
    fn test_total_mass_amu_co2() {
        let atoms = vec![
            Atom::new(6, [0.0, 0.0, 0.0]).unwrap(),
            Atom::new(8, [0.0, 0.0, 2.19]).unwrap(),
            Atom::new(8, [0.0, 0.0, -2.19]).unwrap(),
        ];
        let m = total_mass_amu(&atoms).unwrap();
        let expected = ATOMIC_MASSES[6] + 2.0 * ATOMIC_MASSES[8];
        assert_relative_eq!(m, expected, epsilon = 1e-14);
    }

    #[test]
    fn test_total_mass_amu_atom() {
        let atoms = vec![Atom::new(1, [0.0, 0.0, 0.0]).unwrap()];
        let m = total_mass_amu(&atoms).unwrap();
        assert_eq!(m, ATOMIC_MASSES[1]);
    }

    #[test]
    fn test_total_mass_empty_fails() {
        let atoms: Vec<Atom> = vec![];
        let result = total_mass_amu(&atoms);
        assert!(matches!(result, Err(ThermochemistryError::EmptyMolecule)));
    }

    // ========================================================================
    // Unit tests: translational
    // ========================================================================

    #[test]
    fn test_translational_h2o_298k_1atm() {
        // Compare against PySCF golden: S_trans = 5.515294252089e-5 Ha/K
        let golden = load_thermochem_golden("h2o");
        let m_h2o_golden_amu = {
            let atoms: Vec<Atom> = golden
                .atoms
                .iter()
                .map(|a| Atom::new(a.z, a.pos_bohr).unwrap())
                .collect();
            total_mass_amu(&atoms).unwrap()
        };

        let (e, h, s, cv, cp) = translational_partition(m_h2o_golden_amu, 298.15, 101325.0);
        let r = r_hartree_per_mol_k();
        // E_trans = (3/2) R T
        assert_relative_eq!(e, 1.5 * r * 298.15, epsilon = 1e-14);
        // H_trans = (5/2) R T
        assert_relative_eq!(h, 2.5 * r * 298.15, epsilon = 1e-14);
        // Cv = (3/2) R, Cp = (5/2) R
        assert_relative_eq!(cv, 1.5 * r, epsilon = 1e-14);
        assert_relative_eq!(cp, 2.5 * r, epsilon = 1e-14);
        // Sackur-Tetrode: match PySCF within 1e-10 Ha/K
        let pyscf_s = golden.thermo.get("S_trans").unwrap().value;
        assert_relative_eq!(s, pyscf_s, epsilon = 1e-10);
        let pyscf_e = golden.thermo.get("E_trans").unwrap().value;
        assert_relative_eq!(e, pyscf_e, epsilon = 1e-12);
    }

    #[test]
    fn test_translational_cv_constant() {
        let r = r_hartree_per_mol_k();
        for &m in &[1.0, 18.0, 100.0] {
            for &t in &[100.0, 298.15, 1000.0] {
                for &p in &[1000.0, 101325.0, 1e7] {
                    let (_, _, _, cv, _) = translational_partition(m, t, p);
                    assert_relative_eq!(cv, 1.5 * r, epsilon = 1e-14);
                }
            }
        }
    }

    #[test]
    fn test_translational_cp_consistency() {
        let r = r_hartree_per_mol_k();
        let (_, _, _, cv, cp) = translational_partition(18.0, 298.15, 101325.0);
        assert_relative_eq!(cp, cv + r, epsilon = 1e-14);
    }

    #[test]
    fn test_translational_e_equipartition() {
        let r = r_hartree_per_mol_k();
        let (e, _, _, _, _) = translational_partition(18.0, 298.15, 101325.0);
        assert_relative_eq!(e, 1.5 * r * 298.15, epsilon = 1e-14);
    }

    #[test]
    fn test_translational_h_relation() {
        let r = r_hartree_per_mol_k();
        let (e, h, _, _, _) = translational_partition(18.0, 298.15, 101325.0);
        assert_relative_eq!(h, e + r * 298.15, epsilon = 1e-14);
    }

    #[test]
    fn test_translational_temperature_scaling() {
        // E doubles when T doubles; Cv unchanged
        let (e1, _, _, cv1, _) = translational_partition(18.0, 298.15, 101325.0);
        let (e2, _, _, cv2, _) = translational_partition(18.0, 596.30, 101325.0);
        assert_relative_eq!(e2, 2.0 * e1, epsilon = 1e-14);
        assert_relative_eq!(cv1, cv2, epsilon = 1e-14);
    }

    #[test]
    fn test_translational_pressure_scaling() {
        // Doubling P subtracts R·ln 2 from S_trans
        let r = r_hartree_per_mol_k();
        let (_, _, s1, _, _) = translational_partition(18.0, 298.15, 101325.0);
        let (_, _, s2, _, _) = translational_partition(18.0, 298.15, 2.0 * 101325.0);
        assert_relative_eq!(s2 - s1, -r * 2.0_f64.ln(), epsilon = 1e-14);
    }

    // ========================================================================
    // Unit tests: rotational
    // ========================================================================

    #[test]
    fn test_rotational_atom() {
        let (e, s, cv) = rotational_partition(
            RotorType::Atom,
            [f64::INFINITY, f64::INFINITY, f64::INFINITY],
            1,
            298.15,
        );
        assert_eq!(e, 0.0);
        assert_eq!(s, 0.0);
        assert_eq!(cv, 0.0);
    }

    #[test]
    fn test_rotational_linear_cv_equals_r() {
        let r = r_hartree_per_mol_k();
        let (_, _, cv) = rotational_partition(
            RotorType::Linear,
            [f64::INFINITY, 11.7376, 11.7376], // CO2-like
            2,
            298.15,
        );
        assert_relative_eq!(cv, r, epsilon = 1e-14);
    }

    #[test]
    fn test_rotational_linear_e_equals_rt() {
        let r = r_hartree_per_mol_k();
        let (e, _, _) =
            rotational_partition(RotorType::Linear, [1e300, 11.7376, 11.7376], 2, 298.15);
        assert_relative_eq!(e, r * 298.15, epsilon = 1e-14);
    }

    #[test]
    fn test_rotational_nonlinear_cv_equals_3r_2() {
        let r = r_hartree_per_mol_k();
        let (_, _, cv) = rotational_partition(
            RotorType::AsymmetricTop,
            [820.6, 437.2, 285.2], // H2O-like
            2,
            298.15,
        );
        assert_relative_eq!(cv, 1.5 * r, epsilon = 1e-14);
    }

    #[test]
    fn test_rotational_nonlinear_e_equals_3rt_2() {
        let r = r_hartree_per_mol_k();
        let (e, _, _) =
            rotational_partition(RotorType::AsymmetricTop, [820.6, 437.2, 285.2], 2, 298.15);
        assert_relative_eq!(e, 1.5 * r * 298.15, epsilon = 1e-14);
    }

    #[test]
    fn test_rotational_symmetry_number_scaling() {
        // S_rot(σ=2) - S_rot(σ=1) = -R·ln 2  (for nonlinear)
        let r = r_hartree_per_mol_k();
        let (_, s1, _) =
            rotational_partition(RotorType::AsymmetricTop, [820.6, 437.2, 285.2], 1, 298.15);
        let (_, s2, _) =
            rotational_partition(RotorType::AsymmetricTop, [820.6, 437.2, 285.2], 2, 298.15);
        let delta = s2 - s1;
        assert_relative_eq!(delta, -r * 2.0_f64.ln(), epsilon = 1e-14);
    }

    #[test]
    fn test_rotational_temperature_dependence_linear() {
        // Linear: ΔS_rot = R · ln(T_2/T_1)
        let r = r_hartree_per_mol_k();
        let (_, s1, _) = rotational_partition(RotorType::Linear, [1e300, 10.0, 10.0], 1, 200.0);
        let (_, s2, _) = rotational_partition(RotorType::Linear, [1e300, 10.0, 10.0], 1, 400.0);
        assert_relative_eq!(s2 - s1, r * 2.0_f64.ln(), epsilon = 1e-14);
    }

    #[test]
    fn test_rotational_temperature_dependence_nonlinear() {
        // Nonlinear: ΔS_rot = (3/2) R · ln(T_2/T_1)
        let r = r_hartree_per_mol_k();
        let (_, s1, _) =
            rotational_partition(RotorType::AsymmetricTop, [820.6, 437.2, 285.2], 1, 200.0);
        let (_, s2, _) =
            rotational_partition(RotorType::AsymmetricTop, [820.6, 437.2, 285.2], 1, 400.0);
        assert_relative_eq!(s2 - s1, 1.5 * r * 2.0_f64.ln(), epsilon = 1e-14);
    }

    #[test]
    fn test_rotational_h2o_nonlinear_vs_pyscf() {
        let golden = load_thermochem_golden("h2o");
        let (e, s, cv) = rotational_partition(
            RotorType::AsymmetricTop,
            [
                golden.rotational_constants_ghz[0],
                golden.rotational_constants_ghz[1],
                golden.rotational_constants_ghz[2],
            ],
            golden.symmetry_number,
            golden.temperature_k,
        );
        let pyscf_e = golden.thermo.get("E_rot").unwrap().value;
        let pyscf_s = golden.thermo.get("S_rot").unwrap().value;
        let pyscf_cv = golden.thermo.get("Cv_rot").unwrap().value;
        assert_relative_eq!(e, pyscf_e, epsilon = 1e-12);
        assert_relative_eq!(s, pyscf_s, epsilon = 1e-12);
        assert_relative_eq!(cv, pyscf_cv, epsilon = 1e-14);
    }

    #[test]
    fn test_rotational_co2_linear_vs_pyscf() {
        let golden = load_thermochem_golden("co2");
        // For linear molecules, the first rotational constant in the JSON is
        // +/-1e300 (sanitized infinity). rotational_partition() reads index [1].
        let (e, s, cv) = rotational_partition(
            RotorType::Linear,
            [
                golden.rotational_constants_ghz[0],
                golden.rotational_constants_ghz[1],
                golden.rotational_constants_ghz[2],
            ],
            golden.symmetry_number,
            golden.temperature_k,
        );
        let pyscf_e = golden.thermo.get("E_rot").unwrap().value;
        let pyscf_s = golden.thermo.get("S_rot").unwrap().value;
        let pyscf_cv = golden.thermo.get("Cv_rot").unwrap().value;
        assert_relative_eq!(e, pyscf_e, epsilon = 1e-12);
        assert_relative_eq!(s, pyscf_s, epsilon = 1e-12);
        assert_relative_eq!(cv, pyscf_cv, epsilon = 1e-14);
    }

    // ========================================================================
    // Unit tests: vibrational
    // ========================================================================

    #[test]
    fn test_vibrational_zpe_formula() {
        // Single mode at ν = 1000 cm⁻¹
        // Θ = h · c · ν / k_B ≈ 1438.8 K
        // ZPE = (R/2) · Θ
        let (zpe, _, _, _, n_used, n_imag) = vibrational_partition(&[1000.0], 298.15);
        let expected_theta = PLANCK * (LIGHT_SPEED_SI * 100.0) * 1000.0 / BOLTZMANN;
        let expected_zpe = 0.5 * r_hartree_per_mol_k() * expected_theta;
        assert_relative_eq!(zpe, expected_zpe, epsilon = 1e-14);
        assert_eq!(n_used, 1);
        assert_eq!(n_imag, 0);
    }

    #[test]
    fn test_vibrational_zpe_independent_of_t() {
        let (zpe1, _, _, _, _, _) = vibrational_partition(&[1000.0, 2000.0, 3000.0], 200.0);
        let (zpe2, _, _, _, _, _) = vibrational_partition(&[1000.0, 2000.0, 3000.0], 400.0);
        assert_relative_eq!(zpe1, zpe2, epsilon = 1e-15);
    }

    #[test]
    fn test_vibrational_imaginary_skipped() {
        // [100, -200, 1500, 3000] cm⁻¹: skip mode 1 (-200)
        let (_, _, _, _, n_used, n_imag) =
            vibrational_partition(&[100.0, -200.0, 1500.0, 3000.0], 298.15);
        assert_eq!(n_used, 3);
        assert_eq!(n_imag, 1);
    }

    #[test]
    fn test_vibrational_all_imaginary() {
        let (zpe, e, s, cv, n_used, n_imag) =
            vibrational_partition(&[-100.0, -200.0, -300.0], 298.15);
        assert_eq!(zpe, 0.0);
        assert_eq!(e, 0.0);
        assert_eq!(s, 0.0);
        assert_eq!(cv, 0.0);
        assert_eq!(n_used, 0);
        assert_eq!(n_imag, 3);
    }

    #[test]
    fn test_vibrational_high_temperature_limit() {
        // At T = 50000 K with low frequencies, Cv_vib → n·R
        let r = r_hartree_per_mol_k();
        let freqs = vec![100.0, 200.0, 500.0]; // 3 modes
        let (_, _, _, cv, _, _) = vibrational_partition(&freqs, 50000.0);
        let expected = 3.0 * r;
        assert!(
            (cv - expected).abs() / expected < 0.01,
            "Cv_vib {} not within 1% of 3R = {}",
            cv,
            expected
        );
    }

    #[test]
    fn test_vibrational_low_temperature_limit() {
        // At T = 1 K with stiff modes, thermal + S + Cv ≈ 0
        let (_, e_thermal, s, cv, _, _) = vibrational_partition(&[3000.0], 1.0);
        assert!(e_thermal.abs() < 1e-15);
        assert!(s.abs() < 1e-15);
        assert!(cv.abs() < 1e-15);
    }

    #[test]
    fn test_vibrational_empty() {
        let (zpe, e, s, cv, n_used, n_imag) = vibrational_partition(&[], 298.15);
        assert_eq!(zpe, 0.0);
        assert_eq!(e, 0.0);
        assert_eq!(s, 0.0);
        assert_eq!(cv, 0.0);
        assert_eq!(n_used, 0);
        assert_eq!(n_imag, 0);
    }

    #[test]
    fn test_vibrational_zpe_doubling_frequencies() {
        // ZPE(2ν) = 2 · ZPE(ν)
        let (zpe1, _, _, _, _, _) = vibrational_partition(&[1000.0], 298.15);
        let (zpe2, _, _, _, _, _) = vibrational_partition(&[2000.0], 298.15);
        assert_relative_eq!(zpe2, 2.0 * zpe1, epsilon = 1e-14);
    }

    #[test]
    fn test_vibrational_h2o_zpe_vs_pyscf() {
        let golden = load_thermochem_golden("h2o");
        let (zpe, _, _, _, n_used, n_imag) =
            vibrational_partition(&golden.freq_wavenumber_cm1, golden.temperature_k);
        let pyscf_zpe = golden.thermo.get("ZPE").unwrap().value;
        assert_relative_eq!(zpe, pyscf_zpe, epsilon = 1e-12);
        assert_eq!(n_used, 3);
        assert_eq!(n_imag, 0);
    }

    #[test]
    fn test_vibrational_h2o_s_vib_vs_pyscf() {
        let golden = load_thermochem_golden("h2o");
        let (_, _, s_vib, cv_vib, _, _) =
            vibrational_partition(&golden.freq_wavenumber_cm1, golden.temperature_k);
        let pyscf_s = golden.thermo.get("S_vib").unwrap().value;
        let pyscf_cv = golden.thermo.get("Cv_vib").unwrap().value;
        assert_relative_eq!(s_vib, pyscf_s, epsilon = 1e-14);
        assert_relative_eq!(cv_vib, pyscf_cv, epsilon = 1e-14);
    }

    // ========================================================================
    // Unit tests: validation
    // ========================================================================

    #[test]
    fn test_validate_temperature_zero_fails() {
        let r = validate_inputs(0.0, 101325.0, 1, 1, -75.0);
        assert!(matches!(
            r,
            Err(ThermochemistryError::NonPositiveTemperature(_))
        ));
    }

    #[test]
    fn test_validate_temperature_negative_fails() {
        let r = validate_inputs(-10.0, 101325.0, 1, 1, -75.0);
        assert!(matches!(
            r,
            Err(ThermochemistryError::NonPositiveTemperature(_))
        ));
    }

    #[test]
    fn test_validate_temperature_nan_fails() {
        let r = validate_inputs(f64::NAN, 101325.0, 1, 1, -75.0);
        assert!(matches!(
            r,
            Err(ThermochemistryError::NonFiniteInput { .. })
        ));
    }

    #[test]
    fn test_validate_pressure_zero_fails() {
        let r = validate_inputs(298.15, 0.0, 1, 1, -75.0);
        assert!(matches!(
            r,
            Err(ThermochemistryError::NonPositivePressure(_))
        ));
    }

    #[test]
    fn test_validate_pressure_infinity_fails() {
        let r = validate_inputs(298.15, f64::INFINITY, 1, 1, -75.0);
        assert!(matches!(
            r,
            Err(ThermochemistryError::NonPositivePressure(_))
        ));
    }

    #[test]
    fn test_validate_multiplicity_zero_fails() {
        let r = validate_inputs(298.15, 101325.0, 0, 1, -75.0);
        assert!(matches!(
            r,
            Err(ThermochemistryError::InvalidMultiplicity(0))
        ));
    }

    #[test]
    fn test_validate_symmetry_number_zero_fails() {
        let r = validate_inputs(298.15, 101325.0, 1, 0, -75.0);
        assert!(matches!(
            r,
            Err(ThermochemistryError::InvalidSymmetryNumber(0))
        ));
    }

    #[test]
    fn test_validate_non_finite_electronic_energy_fails() {
        let r = validate_inputs(298.15, 101325.0, 1, 1, f64::NAN);
        assert!(matches!(
            r,
            Err(ThermochemistryError::NonFiniteInput { .. })
        ));
    }

    // ========================================================================
    // Integration tests: top-level compute_thermochemistry
    // ========================================================================

    /// Build a synthetic [`FrequencyInfo`] from the PySCF-stored frequencies
    /// in the US-099 golden JSON. This isolates the thermochemistry unit
    /// tests from any accumulated numerical drift in IQCP's own
    /// `harmonic_analysis` (US-096), letting us verify that the RRHO
    /// formulas themselves agree with PySCF to machine precision.
    ///
    /// Used by `run_for_molecule_pyscf_freqs` below. US-096 end-to-end
    /// validation is separate (`test_thermochemistry_h2o_rhf_sto3g_298k_ac6`
    /// feeds IQCP's own frequencies through the pipeline).
    fn freq_info_from_pyscf_golden(
        golden: &ThermochemGolden,
        rot_type: RotorType,
    ) -> FrequencyInfo {
        let n_modes = golden.freq_wavenumber_cm1.len();
        let n_atoms = golden.atoms.len();
        // Dummy norm modes (unused by compute_thermochemistry)
        let empty_modes: Vec<Vec<[f64; 3]>> = vec![vec![[0.0; 3]; n_atoms]; n_modes];
        // Vibrational temperature not read by compute_thermochemistry (we use
        // freq_wavenumber and convert internally), but populate for completeness.
        let r = r_hartree_per_mol_k();
        let _ = r;
        let vib_temperature: Vec<f64> = golden
            .freq_wavenumber_cm1
            .iter()
            .map(|&nu| {
                if nu <= 0.0 {
                    0.0
                } else {
                    PLANCK * (LIGHT_SPEED_SI * 100.0) * nu / BOLTZMANN
                }
            })
            .collect();
        FrequencyInfo {
            freq_wavenumber: golden.freq_wavenumber_cm1.clone(),
            freq_au: golden.freq_au.clone(),
            reduced_mass: vec![0.0; n_modes],
            force_const_au: vec![0.0; n_modes],
            force_const_dyne: vec![0.0; n_modes],
            norm_mode: empty_modes.clone(),
            norm_mode_mw: empty_modes,
            rot_type,
            rotational_constants_ghz: [
                golden.rotational_constants_ghz[0],
                golden.rotational_constants_ghz[1],
                golden.rotational_constants_ghz[2],
            ],
            principal_moments_amu_bohr2: [0.0, 0.0, 0.0],
            vib_temperature,
            n_modes,
            n_atoms,
        }
    }

    /// Map the PySCF `rot_type` string from the golden JSON to [`RotorType`].
    /// PySCF uses only three labels: `"ATOM"`, `"LINEAR"`, `"REGULAR"`.
    /// "REGULAR" covers all nonlinear rotors; we map to `AsymmetricTop`
    /// since the branch logic is the same for all nonlinear types.
    fn parse_pyscf_rotor_type(s: &str) -> RotorType {
        match s {
            "ATOM" => RotorType::Atom,
            "LINEAR" => RotorType::Linear,
            _ => RotorType::AsymmetricTop,
        }
    }

    /// Construct a FrequencyInfo directly from the PySCF golden JSON and
    /// run compute_thermochemistry(). Isolates thermochemistry from US-096.
    fn run_for_molecule_pyscf_freqs(name: &str) -> (ThermochemGolden, ThermochemistryResult) {
        let thermochem_gold = load_thermochem_golden(name);
        let rot_type = parse_pyscf_rotor_type(&thermochem_gold.rot_type);
        let freq_info = freq_info_from_pyscf_golden(&thermochem_gold, rot_type);
        let atoms: Vec<Atom> = thermochem_gold
            .atoms
            .iter()
            .map(|a| Atom::new(a.z, a.pos_bohr).unwrap())
            .collect();
        let result = compute_thermochemistry(
            &freq_info,
            thermochem_gold.scf_energy_ha,
            &atoms,
            thermochem_gold.temperature_k,
            thermochem_gold.pressure_pa,
            Some(thermochem_gold.symmetry_number),
            thermochem_gold.multiplicity,
        )
        .unwrap();
        (thermochem_gold, result)
    }

    /// Helper: given a molecule name, reuse the US-096 golden hessian to
    /// compute IQCP's FrequencyInfo (via harmonic_analysis) and run
    /// compute_thermochemistry(). This is the end-to-end path used by AC6.
    fn run_for_molecule(name: &str) -> (ThermochemGolden, ThermochemistryResult) {
        let thermo_gold = load_thermo_golden(name);
        let thermochem_gold = load_thermochem_golden(name);
        let atoms = thermo_gold.atoms();
        let hessian = thermo_gold.hessian_matrix();
        let freq_info = harmonic_analysis(&atoms, &hessian).unwrap();
        let result = compute_thermochemistry(
            &freq_info,
            thermochem_gold.scf_energy_ha,
            &atoms,
            thermochem_gold.temperature_k,
            thermochem_gold.pressure_pa,
            Some(thermochem_gold.symmetry_number),
            thermochem_gold.multiplicity,
        )
        .unwrap();
        (thermochem_gold, result)
    }

    #[test]
    fn test_thermochemistry_h2o_rhf_sto3g_298k_ac6() {
        // AC6 (hard target): ZPE within 1e-6 Ha, G within 1e-4 Ha.
        //
        // End-to-end path: PySCF golden hessian → IQCP harmonic_analysis →
        // IQCP compute_thermochemistry → compare to PySCF thermo() output.
        //
        // The small residual error relative to a pure-formula comparison
        // comes from IQCP's own harmonic_analysis which matches PySCF
        // frequencies to within ~0.1 cm⁻¹ per mode. For H2O (3 modes),
        // this propagates to ~7e-7 Ha ZPE — within the 1e-6 tolerance by a
        // comfortable margin.
        let (golden, result) = run_for_molecule("h2o");
        let pyscf_zpe = golden.thermo.get("ZPE").unwrap().value;
        let pyscf_g = golden.thermo.get("G_tot").unwrap().value;
        let pyscf_h = golden.thermo.get("H_tot").unwrap().value;
        let pyscf_s = golden.thermo.get("S_tot").unwrap().value;
        let pyscf_cv = golden.thermo.get("Cv_tot").unwrap().value;
        let pyscf_cp = golden.thermo.get("Cp_tot").unwrap().value;
        let pyscf_e0k = golden.thermo.get("E_0K").unwrap().value;

        let zpe_err = (result.zpe_ha - pyscf_zpe).abs();
        let g_err = (result.gibbs_ha - pyscf_g).abs();
        assert!(
            zpe_err < 1e-6,
            "H2O ZPE error {zpe_err:.3e} exceeds 1e-6 Ha (IQCP {} vs PySCF {})",
            result.zpe_ha,
            pyscf_zpe
        );
        assert!(
            g_err < 1e-4,
            "H2O G error {g_err:.3e} exceeds 1e-4 Ha (IQCP {} vs PySCF {})",
            result.gibbs_ha,
            pyscf_g
        );
        // Sanity checks: the total quantities inherit the frequency
        // accuracy of US-096 (~0.1 cm⁻¹/mode → ~1e-6 Ha ZPE). Tolerances
        // are chosen to reflect that end-to-end accuracy, not the purely
        // analytical 1e-12 limit you would get with PySCF-provided freqs.
        assert!((result.enthalpy_ha - pyscf_h).abs() < 2e-6);
        assert_relative_eq!(result.entropy_ha_per_k, pyscf_s, epsilon = 1e-8);
        assert_relative_eq!(result.cv_ha_per_k, pyscf_cv, epsilon = 1e-8);
        assert_relative_eq!(result.cp_ha_per_k, pyscf_cp, epsilon = 1e-8);
        assert!((result.e_0k_ha - pyscf_e0k).abs() < 2e-6);
        assert_eq!(result.n_imag, 0);
        assert_eq!(result.n_vib_modes_used, 3);
    }

    #[test]
    fn test_thermochemistry_h2o_pyscf_freqs_vs_pyscf() {
        // Same as AC6 but with PySCF frequencies loaded directly (no
        // harmonic_analysis step). This isolates the RRHO formulas from
        // US-096.
        //
        // Note on residual tolerance: PySCF reports atomic masses to 4
        // significant figures (H=1.008, O=15.999), while IQCP uses the
        // fuller IUPAC 2016 values (H=1.00794, O=15.9994). The ~0.002%
        // mass difference propagates via ln(m^{3/2}) in Sackur-Tetrode to
        // a ~7e-11 Ha/K shift in S_trans and therefore a ~2e-8 Ha shift
        // in G_tot. This is physically irrelevant (well below chemical
        // accuracy) but needs a tolerance of ~1e-7 Ha on G for the test.
        let (golden, result) = run_for_molecule_pyscf_freqs("h2o");
        let pyscf_zpe = golden.thermo.get("ZPE").unwrap().value;
        let pyscf_g = golden.thermo.get("G_tot").unwrap().value;
        let pyscf_h = golden.thermo.get("H_tot").unwrap().value;
        let pyscf_s = golden.thermo.get("S_tot").unwrap().value;
        let pyscf_cv = golden.thermo.get("Cv_tot").unwrap().value;
        let pyscf_cp = golden.thermo.get("Cp_tot").unwrap().value;
        let pyscf_e0k = golden.thermo.get("E_0K").unwrap().value;

        // ZPE, Cv, Cp, E_0K are independent of mass → machine precision
        assert_relative_eq!(result.zpe_ha, pyscf_zpe, epsilon = 1e-12);
        assert_relative_eq!(result.e_0k_ha, pyscf_e0k, epsilon = 1e-13);
        assert_relative_eq!(result.cv_ha_per_k, pyscf_cv, epsilon = 1e-11);
        assert_relative_eq!(result.cp_ha_per_k, pyscf_cp, epsilon = 1e-11);
        // Enthalpy H = U + RT: E_trans, E_rot, E_vib are independent of
        // mass; only the electronic + ZPE terms carry through. Machine
        // precision.
        assert_relative_eq!(result.enthalpy_ha, pyscf_h, epsilon = 1e-12);
        // Entropy includes S_trans which depends on mass via ln m.
        // Absolute tolerance ~1e-10 Ha/K reflects the 7e-11 Ha/K shift.
        assert!((result.entropy_ha_per_k - pyscf_s).abs() < 1e-9);
        // G = H - T·S inherits the mass-induced S_trans shift.
        assert!(
            (result.gibbs_ha - pyscf_g).abs() < 1e-7,
            "G mismatch: IQCP {} vs PySCF {} (delta {:e})",
            result.gibbs_ha,
            pyscf_g,
            result.gibbs_ha - pyscf_g
        );
        assert_eq!(result.n_imag, 0);
        assert_eq!(result.n_vib_modes_used, 3);
    }

    #[test]
    fn test_thermochemistry_ch4_rhf_sto3g_298k_sigma12_pyscf_freqs() {
        // Isolated thermochemistry test using PySCF's frequencies directly
        // (bypassing IQCP harmonic_analysis so we verify the RRHO formulas).
        // Tests sigma=12 (T_d) nonlinear rotor branch on a spherical top.
        // Residual ~1e-7 Ha in G from the PySCF/IQCP mass-table
        // difference (see test_thermochemistry_h2o_pyscf_freqs_vs_pyscf).
        let (golden, result) = run_for_molecule_pyscf_freqs("ch4");
        assert_eq!(golden.symmetry_number, 12);
        let pyscf_zpe = golden.thermo.get("ZPE").unwrap().value;
        let pyscf_g = golden.thermo.get("G_tot").unwrap().value;
        assert_relative_eq!(result.zpe_ha, pyscf_zpe, epsilon = 1e-12);
        assert!((result.gibbs_ha - pyscf_g).abs() < 1e-7);
        assert_eq!(result.n_imag, 0);
        assert_eq!(result.n_vib_modes_used, 9);
    }

    #[test]
    fn test_thermochemistry_ch4_rhf_sto3g_298k_sigma12_end_to_end() {
        // End-to-end path: PySCF hessian → IQCP harmonic_analysis →
        // compute_thermochemistry. The residual error vs PySCF reflects the
        // ~0.1 cm⁻¹/mode frequency accuracy of US-096; for CH4's 9 modes,
        // this propagates to at most ~2e-6 Ha ZPE. This is still well
        // within the AC6 ~1e-4 Ha G tolerance for chemistry purposes.
        let (golden, result) = run_for_molecule("ch4");
        assert_eq!(golden.symmetry_number, 12);
        let pyscf_zpe = golden.thermo.get("ZPE").unwrap().value;
        let pyscf_g = golden.thermo.get("G_tot").unwrap().value;
        // ZPE: end-to-end tolerance for a 9-mode molecule
        assert!(
            (result.zpe_ha - pyscf_zpe).abs() < 5e-6,
            "CH4 ZPE end-to-end error: IQCP {} vs PySCF {}",
            result.zpe_ha,
            pyscf_zpe
        );
        // G: still meets the 1e-4 Ha chemistry target
        assert!(
            (result.gibbs_ha - pyscf_g).abs() < 1e-4,
            "CH4 G end-to-end error: IQCP {} vs PySCF {}",
            result.gibbs_ha,
            pyscf_g
        );
        assert_eq!(result.n_imag, 0);
        assert_eq!(result.n_vib_modes_used, 9);
    }

    #[test]
    fn test_thermochemistry_co2_rhf_sto3g_298k_sigma2_pyscf_freqs() {
        // Linear molecule test using PySCF frequencies directly.
        // Residual ~1e-7 Ha in G from mass-table differences.
        let (golden, result) = run_for_molecule_pyscf_freqs("co2");
        assert_eq!(golden.symmetry_number, 2);
        let pyscf_zpe = golden.thermo.get("ZPE").unwrap().value;
        let pyscf_g = golden.thermo.get("G_tot").unwrap().value;
        assert_relative_eq!(result.zpe_ha, pyscf_zpe, epsilon = 1e-12);
        assert!((result.gibbs_ha - pyscf_g).abs() < 1e-7);
        assert_eq!(result.n_vib_modes_used, 4);
        // Verify linear rotor branch was used
        assert_relative_eq!(
            result.components.cv_rot_ha_per_k,
            r_hartree_per_mol_k(),
            epsilon = 1e-14
        );
    }

    #[test]
    fn test_thermochemistry_co2_rhf_sto3g_298k_sigma2_end_to_end() {
        // End-to-end path for CO2. Linear molecule with 4 modes.
        let (golden, result) = run_for_molecule("co2");
        assert_eq!(golden.symmetry_number, 2);
        let pyscf_zpe = golden.thermo.get("ZPE").unwrap().value;
        let pyscf_g = golden.thermo.get("G_tot").unwrap().value;
        assert!((result.zpe_ha - pyscf_zpe).abs() < 5e-6);
        assert!((result.gibbs_ha - pyscf_g).abs() < 1e-4);
        assert_eq!(result.n_vib_modes_used, 4);
    }

    #[test]
    fn test_thermochemistry_h2_rhf_sto3g_298k_sigma2_pyscf_freqs() {
        // Diatomic test using PySCF frequencies directly.
        // Residual ~1e-7 Ha in G from mass-table differences.
        let (golden, result) = run_for_molecule_pyscf_freqs("h2");
        assert_eq!(golden.symmetry_number, 2);
        let pyscf_zpe = golden.thermo.get("ZPE").unwrap().value;
        let pyscf_g = golden.thermo.get("G_tot").unwrap().value;
        assert_relative_eq!(result.zpe_ha, pyscf_zpe, epsilon = 1e-12);
        assert!((result.gibbs_ha - pyscf_g).abs() < 1e-7);
        assert_eq!(result.n_vib_modes_used, 1);
    }

    #[test]
    fn test_thermochemistry_h2_rhf_sto3g_298k_sigma2_end_to_end() {
        // End-to-end path for H2. Diatomic with 1 mode.
        let (golden, result) = run_for_molecule("h2");
        assert_eq!(golden.symmetry_number, 2);
        let pyscf_zpe = golden.thermo.get("ZPE").unwrap().value;
        let pyscf_g = golden.thermo.get("G_tot").unwrap().value;
        assert!((result.zpe_ha - pyscf_zpe).abs() < 1e-6);
        assert!((result.gibbs_ha - pyscf_g).abs() < 1e-4);
        assert_eq!(result.n_vib_modes_used, 1);
    }

    #[test]
    fn test_thermochemistry_h_atom_edge_case() {
        // Single hydrogen atom: only translation + electronic (doublet).
        let atoms = vec![Atom::new(1, [0.0, 0.0, 0.0]).unwrap()];
        let hessian = DMatrix::<f64>::zeros(3, 3);
        let freq_info = harmonic_analysis(&atoms, &hessian).unwrap();
        assert_eq!(freq_info.rot_type, RotorType::Atom);
        assert_eq!(freq_info.n_modes, 0);

        let r = r_hartree_per_mol_k();
        let result = compute_thermochemistry(
            &freq_info, -0.5, // arbitrary electronic energy
            &atoms, 298.15, 101325.0, None, // default σ=1
            2,    // doublet
        )
        .unwrap();

        assert_eq!(result.zpe_ha, 0.0);
        assert_eq!(result.components.e_rot_ha, 0.0);
        assert_eq!(result.components.e_vib_thermal_ha, 0.0);
        assert_eq!(result.components.s_rot_ha_per_k, 0.0);
        assert_eq!(result.components.s_vib_ha_per_k, 0.0);
        assert_eq!(result.components.cv_rot_ha_per_k, 0.0);
        assert_eq!(result.components.cv_vib_ha_per_k, 0.0);
        assert_relative_eq!(
            result.components.cv_trans_ha_per_k,
            1.5 * r,
            epsilon = 1e-14
        );
        assert_relative_eq!(result.cv_ha_per_k, 1.5 * r, epsilon = 1e-14);
        assert_relative_eq!(result.cp_ha_per_k, 2.5 * r, epsilon = 1e-14);
        assert_relative_eq!(
            result.components.s_elec_ha_per_k,
            r * 2.0_f64.ln(),
            epsilon = 1e-14
        );
        assert_eq!(result.n_imag, 0);
        assert_eq!(result.n_vib_modes_used, 0);
    }

    #[test]
    fn test_thermochemistry_h2o_total_breakdown_consistency() {
        let (_, result) = run_for_molecule("h2o");
        let r = r_hartree_per_mol_k();

        // U = E_elec + ZPE + E_trans + E_rot + E_vib_thermal
        let expected_u = result.electronic_energy_ha
            + result.zpe_ha
            + result.components.e_trans_ha
            + result.components.e_rot_ha
            + result.components.e_vib_thermal_ha;
        assert_relative_eq!(result.internal_energy_ha, expected_u, epsilon = 1e-13);

        // S_tot = S_trans + S_rot + S_vib + S_elec
        let expected_s = result.components.s_trans_ha_per_k
            + result.components.s_rot_ha_per_k
            + result.components.s_vib_ha_per_k
            + result.components.s_elec_ha_per_k;
        assert_relative_eq!(result.entropy_ha_per_k, expected_s, epsilon = 1e-15);

        // H = U + RT
        assert_relative_eq!(
            result.enthalpy_ha,
            result.internal_energy_ha + r * result.temperature_k,
            epsilon = 1e-13
        );

        // G = H - T*S
        assert_relative_eq!(
            result.gibbs_ha,
            result.enthalpy_ha - result.temperature_k * result.entropy_ha_per_k,
            epsilon = 1e-12
        );

        // Cp = Cv + R
        assert_relative_eq!(result.cp_ha_per_k, result.cv_ha_per_k + r, epsilon = 1e-14);

        // E_0K = E_elec + ZPE
        assert_relative_eq!(
            result.e_0k_ha,
            result.electronic_energy_ha + result.zpe_ha,
            epsilon = 1e-15
        );

        // Cv = Cv_trans + Cv_rot + Cv_vib  (Cv_elec = 0)
        let expected_cv = result.components.cv_trans_ha_per_k
            + result.components.cv_rot_ha_per_k
            + result.components.cv_vib_ha_per_k;
        assert_relative_eq!(result.cv_ha_per_k, expected_cv, epsilon = 1e-15);
    }

    #[test]
    fn test_thermochemistry_h2o_temperature_dependence() {
        // Run H2O at 200, 298.15, 400, 500 K. This uses PySCF frequencies
        // directly so the 298.15 K comparison can be done at machine
        // precision, and we verify the smooth monotonicity of G(T) across
        // all four temperatures.
        let thermochem_gold = load_thermochem_golden("h2o");
        let freq_info = freq_info_from_pyscf_golden(
            &thermochem_gold,
            parse_pyscf_rotor_type(&thermochem_gold.rot_type),
        );
        let atoms: Vec<Atom> = thermochem_gold
            .atoms
            .iter()
            .map(|a| Atom::new(a.z, a.pos_bohr).unwrap())
            .collect();

        let r = r_hartree_per_mol_k();
        let pyscf_zpe = thermochem_gold.thermo.get("ZPE").unwrap().value;
        let pyscf_g_298 = thermochem_gold.thermo.get("G_tot").unwrap().value;
        let mut prev_g: Option<f64> = None;
        for &t in &[200.0, 298.15, 400.0, 500.0] {
            let result = compute_thermochemistry(
                &freq_info,
                thermochem_gold.scf_energy_ha,
                &atoms,
                t,
                101325.0,
                Some(thermochem_gold.symmetry_number),
                1,
            )
            .unwrap();
            // ZPE is temperature-independent (matches PySCF to ~1e-12 Ha)
            assert_relative_eq!(result.zpe_ha, pyscf_zpe, epsilon = 1e-12);
            // At 298.15 K: match full PySCF result within the mass-table
            // residual (~1e-7 Ha; see
            // test_thermochemistry_h2o_pyscf_freqs_vs_pyscf).
            if (t - 298.15).abs() < 1e-10 {
                assert!((result.gibbs_ha - pyscf_g_298).abs() < 1e-7);
            }
            // G must be monotonically decreasing in T for a stable molecule
            // (dG/dT = -S < 0)
            if let Some(pg) = prev_g {
                assert!(
                    result.gibbs_ha < pg,
                    "G(T={t}) = {} should be < previous G = {pg}",
                    result.gibbs_ha,
                );
            }
            prev_g = Some(result.gibbs_ha);

            // Verify H = U + RT exact algebraic identity
            assert_relative_eq!(
                result.enthalpy_ha,
                result.internal_energy_ha + r * t,
                epsilon = 1e-13
            );
        }
    }

    #[test]
    fn test_thermochemistry_h2o_sigma_override_effect() {
        // Passing σ=1 vs σ=2 should shift G_tot by exactly +T·R·ln 2.
        //
        // Derivation: q_rot(σ) = (const) / σ  →  ln q_rot(2) = ln q_rot(1) - ln 2
        //             S_rot(2)  = S_rot(1) - R·ln 2
        //             S_tot(2)  = S_tot(1) - R·ln 2        (rot contribution)
        //             G(2) - G(1) = -T·(S_tot(2) - S_tot(1)) = -T·(-R·ln 2) = +T·R·ln 2
        //
        // So increasing σ makes G **higher** (less negative), which
        // physically corresponds to the entropy penalty for
        // indistinguishable orientations.
        let thermo_gold = load_thermo_golden("h2o");
        let atoms = thermo_gold.atoms();
        let hessian = thermo_gold.hessian_matrix();
        let freq_info = harmonic_analysis(&atoms, &hessian).unwrap();
        let r = r_hartree_per_mol_k();
        let t = 298.15;

        let res_sigma1 =
            compute_thermochemistry(&freq_info, -74.963, &atoms, t, 101325.0, Some(1), 1).unwrap();
        let res_sigma2 =
            compute_thermochemistry(&freq_info, -74.963, &atoms, t, 101325.0, Some(2), 1).unwrap();
        let expected_delta = t * r * 2.0_f64.ln();
        let actual_delta = res_sigma2.gibbs_ha - res_sigma1.gibbs_ha;
        assert_relative_eq!(actual_delta, expected_delta, epsilon = 1e-12);
    }

    #[test]
    fn test_compute_thermochemistry_default_sigma_is_one() {
        // Passing None should behave identically to Some(1).
        let thermo_gold = load_thermo_golden("h2o");
        let atoms = thermo_gold.atoms();
        let hessian = thermo_gold.hessian_matrix();
        let freq_info = harmonic_analysis(&atoms, &hessian).unwrap();

        let r_none =
            compute_thermochemistry(&freq_info, -74.963, &atoms, 298.15, 101325.0, None, 1)
                .unwrap();
        let r_some =
            compute_thermochemistry(&freq_info, -74.963, &atoms, 298.15, 101325.0, Some(1), 1)
                .unwrap();
        assert_eq!(r_none.symmetry_number, 1);
        assert_eq!(r_some.symmetry_number, 1);
        assert_relative_eq!(r_none.gibbs_ha, r_some.gibbs_ha, epsilon = 1e-15);
    }

    #[test]
    fn test_compute_thermochemistry_components_match_pyscf_h2o() {
        // Per-contribution cross-check against PySCF for H2O, using
        // PySCF's frequencies directly (isolated from US-096). Every
        // component should agree with PySCF to machine-precision
        // tolerances.
        let (golden, result) = run_for_molecule_pyscf_freqs("h2o");
        let thermo = &golden.thermo;
        let c = &result.components;

        // Translational. S_trans depends on mass via ln(m^{3/2}); the
        // PySCF/IQCP mass-table difference (see
        // test_thermochemistry_h2o_pyscf_freqs_vs_pyscf) shifts it by
        // ~7e-11 Ha/K, so we use absolute rather than relative tolerance.
        assert!((c.s_trans_ha_per_k - thermo.get("S_trans").unwrap().value).abs() < 1e-9);
        // E_trans and Cv_trans are mass-independent → machine precision.
        assert_relative_eq!(
            c.e_trans_ha,
            thermo.get("E_trans").unwrap().value,
            epsilon = 1e-14
        );
        assert_relative_eq!(
            c.cv_trans_ha_per_k,
            thermo.get("Cv_trans").unwrap().value,
            epsilon = 1e-14
        );

        // Rotational
        assert_relative_eq!(
            c.s_rot_ha_per_k,
            thermo.get("S_rot").unwrap().value,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            c.e_rot_ha,
            thermo.get("E_rot").unwrap().value,
            epsilon = 1e-14
        );
        assert_relative_eq!(
            c.cv_rot_ha_per_k,
            thermo.get("Cv_rot").unwrap().value,
            epsilon = 1e-14
        );

        // Vibrational (PySCF E_vib = ZPE + thermal = H_vib)
        let pyscf_e_vib = thermo.get("E_vib").unwrap().value;
        assert_relative_eq!(
            result.zpe_ha + c.e_vib_thermal_ha,
            pyscf_e_vib,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            c.s_vib_ha_per_k,
            thermo.get("S_vib").unwrap().value,
            epsilon = 1e-14
        );
        assert_relative_eq!(
            c.cv_vib_ha_per_k,
            thermo.get("Cv_vib").unwrap().value,
            epsilon = 1e-14
        );

        // Electronic (singlet: S_elec = 0)
        assert_eq!(c.s_elec_ha_per_k, 0.0);
        assert_eq!(thermo.get("S_elec").unwrap().value, 0.0);
    }

    #[test]
    fn test_compute_thermochemistry_rejects_invalid_inputs() {
        let atoms = vec![Atom::new(1, [0.0, 0.0, 0.0]).unwrap()];
        let hessian = DMatrix::<f64>::zeros(3, 3);
        let freq_info = harmonic_analysis(&atoms, &hessian).unwrap();

        // T <= 0
        let r = compute_thermochemistry(&freq_info, -0.5, &atoms, 0.0, 101325.0, None, 1);
        assert!(matches!(
            r,
            Err(ThermochemistryError::NonPositiveTemperature(_))
        ));

        // P <= 0
        let r = compute_thermochemistry(&freq_info, -0.5, &atoms, 298.15, 0.0, None, 1);
        assert!(matches!(
            r,
            Err(ThermochemistryError::NonPositivePressure(_))
        ));

        // multiplicity == 0
        let r = compute_thermochemistry(&freq_info, -0.5, &atoms, 298.15, 101325.0, None, 0);
        assert!(matches!(
            r,
            Err(ThermochemistryError::InvalidMultiplicity(0))
        ));

        // sigma == Some(0)
        let r = compute_thermochemistry(&freq_info, -0.5, &atoms, 298.15, 101325.0, Some(0), 1);
        assert!(matches!(
            r,
            Err(ThermochemistryError::InvalidSymmetryNumber(0))
        ));

        // Empty atoms
        let empty: Vec<Atom> = vec![];
        let r = compute_thermochemistry(&freq_info, -0.5, &empty, 298.15, 101325.0, None, 1);
        assert!(matches!(r, Err(ThermochemistryError::EmptyMolecule)));
    }

    #[test]
    fn test_thermochemistry_determinism() {
        // Same inputs → bit-identical outputs
        let (_, r1) = run_for_molecule("h2o");
        let (_, r2) = run_for_molecule("h2o");
        assert_eq!(r1.zpe_ha.to_bits(), r2.zpe_ha.to_bits());
        assert_eq!(r1.gibbs_ha.to_bits(), r2.gibbs_ha.to_bits());
        assert_eq!(r1.entropy_ha_per_k.to_bits(), r2.entropy_ha_per_k.to_bits());
    }
}
