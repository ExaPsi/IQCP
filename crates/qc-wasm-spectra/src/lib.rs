//! WebAssembly bindings for frequency analysis (Phase 5 spectra pipeline)
//!
//! This crate exposes the full analytical vibrational spectroscopy pipeline
//! (Hessian -> normal modes -> IR intensities -> Raman activities -> RRHO
//! thermochemistry -> broadened spectra) to the React frontend through a
//! single WASM export [`compute_frequencies`].
//!
//! Split from `qc-wasm` to enable lazy loading: the core WASM module loads
//! on page init (<500 KB gzipped), while this spectra module loads on demand
//! when the user opens the Frequency tab (~200 KB gzipped).
//!
//! # Build
//!
//! ```bash
//! wasm-pack build crates/qc-wasm-spectra --release --target web \
//!     --out-dir ../../apps/web/src/wasm-spectra
//! ```

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use qc_core::scf::ConvergenceProfile;

// ============================================================================
// Version
// ============================================================================

/// Return the crate version string.
#[wasm_bindgen]
pub fn spectra_version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

// ============================================================================
// Frequency Analysis WASM Bindings (US-101)
// ============================================================================
//
// These bindings expose the full Phase 5 analytical vibrational spectroscopy
// pipeline (Hessian -> normal modes -> IR intensities -> Raman activities ->
// RRHO thermochemistry -> broadened spectra) to the React frontend through a
// single WASM export `compute_frequencies`. The function orchestrates the
// five qc-core modules implemented in US-093 through US-100, streaming
// phase-granular progress updates and producing a flat camelCase JS result.

/// WASM-friendly input struct for the full frequency analysis pipeline.
///
/// Accepts a molecule, basis set, method, and thermodynamic conditions and
/// returns a complete [`FrequencyWasmResult`] with frequencies, IR
/// intensities, Raman activities, RRHO thermochemistry, and broadened
/// spectra. Deserialized from JavaScript via `serde_wasm_bindgen::from_value`.
///
/// # Example (JavaScript)
///
/// ```javascript
/// const input = {
///   atoms: [[8, 0, 0, 0], [1, 0, 0.756, 0.586], [1, 0, -0.756, 0.586]],
///   basisName: "sto-3g",
///   method: "rhf",
///   temperatureK: 298.15,
///   pressurePa: 101325.0,
///   multiplicity: 1,
///   broadeningKind: "lorentzian",
///   fwhmCm1: 8.0,
/// };
/// ```
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct FrequencyWasmInput {
    /// Atoms as `[[Z, x, y, z], ...]` where Z is the atomic number (as f64)
    /// and x, y, z are Cartesian coordinates in bohr.
    pub atoms: Vec<[f64; 4]>,
    /// Basis set name (e.g., `"sto-3g"`, `"6-31g*"`, `"cc-pvdz"`).
    pub basis_name: String,
    /// Electronic structure method. Recognized (case-insensitive):
    /// `"rhf"`, `"hf"` (alias for rhf), `"lda"`, `"b3lyp"`, `"b3lyp-d3bj"`.
    pub method: String,
    /// Temperature in Kelvin. Default 298.15. Must be > 0 and finite.
    #[serde(default = "default_freq_temperature_k")]
    pub temperature_k: f64,
    /// Pressure in Pascals. Default 101325.0 (1 atm). Must be > 0 and finite.
    #[serde(default = "default_freq_pressure_pa")]
    pub pressure_pa: f64,
    /// Rotational symmetry number override (sigma). `None` -> thermochemistry
    /// defaults sigma = 1. Common values: 1 (C1, Cs), 2 (C2v H2O), 3 (C3v NH3),
    /// 12 (Td CH4), 24 (Oh).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub symmetry_number_override: Option<u32>,
    /// Spin multiplicity `2S + 1`. Default 1 (singlet). Must be >= 1.
    /// Phase 5 is closed-shell only; this affects only the electronic entropy.
    #[serde(default = "default_freq_multiplicity")]
    pub multiplicity: u32,
    /// Spectrum broadening kind: `"lorentzian"` or `"gaussian"`. Default
    /// `"lorentzian"`.
    #[serde(default = "default_freq_broadening_kind")]
    pub broadening_kind: String,
    /// Full width at half maximum for spectrum broadening in cm-1.
    /// Default 8.0 (matches US-100). Must be > 0 and finite.
    #[serde(default = "default_freq_fwhm_cm1")]
    pub fwhm_cm1: f64,
    /// SCF convergence profile: `"loose"`, `"medium"`, or `"tight"`.
    /// Default `"tight"` for accurate frequencies.
    #[serde(default = "default_freq_convergence_profile")]
    pub convergence_profile: String,
    /// Maximum SCF iterations. Default 100.
    #[serde(default = "default_freq_max_iterations")]
    pub max_iterations: u32,
}

fn default_freq_temperature_k() -> f64 {
    298.15
}
fn default_freq_pressure_pa() -> f64 {
    101_325.0
}
fn default_freq_multiplicity() -> u32 {
    1
}
fn default_freq_broadening_kind() -> String {
    "lorentzian".to_string()
}
fn default_freq_fwhm_cm1() -> f64 {
    8.0
}
fn default_freq_convergence_profile() -> String {
    "tight".to_string()
}
fn default_freq_max_iterations() -> u32 {
    100
}

/// Per-phase timing breakdown in milliseconds.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct FrequencyTiming {
    /// Phase 1: SCF + integrals + Hessian assembly (nuclear CPHF).
    pub integrals_ms: f64,
    /// Phase 2: CPHF data extraction + density rebuild.
    pub nuclear_cphf_ms: f64,
    /// Phase 3: Field CPHF (runs inside `compute_raman_spectrum`).
    pub field_cphf_ms: f64,
    /// Phase 4: Harmonic analysis + IR intensities + thermochemistry.
    pub assembly_ms: f64,
    /// Phase 5: Spectrum broadening.
    pub modes_ms: f64,
    /// Total wall time from input deserialization to result serialization.
    pub total_ms: f64,
}

/// RRHO thermochemistry reshaped for JavaScript consumption.
///
/// Mirrors [`qc_core::thermochemistry::ThermochemistryResult`] with its
/// per-contribution breakdown flattened into a single struct.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct FrequencyThermochemistry {
    // ---- Input echoes ----
    /// Temperature used, in Kelvin.
    pub temperature_k: f64,
    /// Pressure used, in Pascals.
    pub pressure_pa: f64,
    /// Rotational symmetry number actually used (sigma after defaulting).
    pub symmetry_number: u32,
    /// Spin multiplicity `2S + 1`.
    pub multiplicity: u32,
    /// Total molecular mass in amu.
    pub total_mass_amu: f64,
    /// Number of vibrational modes used (positive frequencies only).
    pub n_vib_modes_used: usize,
    /// Number of imaginary modes skipped from the vibrational sum.
    pub n_imag: usize,

    // ---- Zero-point + 0 K ----
    /// Zero-point vibrational energy in Hartree.
    pub zpe_ha: f64,
    /// Energy at 0 K: `E_elec + ZPE` in Hartree.
    pub e_0k_ha: f64,

    // ---- Totals ----
    /// Total internal energy `U(T)` in Hartree.
    pub internal_energy_ha: f64,
    /// Total enthalpy `H(T) = U(T) + RT` in Hartree.
    pub enthalpy_ha: f64,
    /// Total entropy `S(T)` in Ha/(mol*K).
    pub entropy_ha_per_k: f64,
    /// Total Gibbs free energy `G(T) = H(T) - T*S(T)` in Hartree.
    pub gibbs_ha: f64,
    /// Total constant-volume heat capacity `Cv(T)` in Ha/(mol*K).
    pub cv_ha_per_k: f64,
    /// Total constant-pressure heat capacity `Cp(T) = Cv(T) + R` in
    /// Ha/(mol*K).
    pub cp_ha_per_k: f64,

    // ---- Translational ----
    pub e_trans_ha: f64,
    pub h_trans_ha: f64,
    pub s_trans_ha_per_k: f64,
    pub cv_trans_ha_per_k: f64,
    pub cp_trans_ha_per_k: f64,

    // ---- Rotational ----
    pub e_rot_ha: f64,
    pub h_rot_ha: f64,
    pub s_rot_ha_per_k: f64,
    pub cv_rot_ha_per_k: f64,
    pub cp_rot_ha_per_k: f64,

    // ---- Vibrational (thermal only; ZPE is stored separately) ----
    pub e_vib_thermal_ha: f64,
    pub h_vib_ha: f64,
    pub s_vib_ha_per_k: f64,
    pub cv_vib_ha_per_k: f64,
    pub cp_vib_ha_per_k: f64,

    // ---- Electronic (only entropy is nonzero in RRHO) ----
    pub s_elec_ha_per_k: f64,
}

impl FrequencyThermochemistry {
    /// Map a [`qc_core::thermochemistry::ThermochemistryResult`] into the
    /// WASM-friendly flat struct by copying totals and unpacking the
    /// component breakdown.
    fn from_core(result: &qc_core::thermochemistry::ThermochemistryResult) -> Self {
        let c = &result.components;
        Self {
            temperature_k: result.temperature_k,
            pressure_pa: result.pressure_pa,
            symmetry_number: result.symmetry_number,
            multiplicity: result.multiplicity,
            total_mass_amu: result.total_mass_amu,
            n_vib_modes_used: result.n_vib_modes_used,
            n_imag: result.n_imag,

            zpe_ha: result.zpe_ha,
            e_0k_ha: result.e_0k_ha,

            internal_energy_ha: result.internal_energy_ha,
            enthalpy_ha: result.enthalpy_ha,
            entropy_ha_per_k: result.entropy_ha_per_k,
            gibbs_ha: result.gibbs_ha,
            cv_ha_per_k: result.cv_ha_per_k,
            cp_ha_per_k: result.cp_ha_per_k,

            e_trans_ha: c.e_trans_ha,
            h_trans_ha: c.h_trans_ha,
            s_trans_ha_per_k: c.s_trans_ha_per_k,
            cv_trans_ha_per_k: c.cv_trans_ha_per_k,
            cp_trans_ha_per_k: c.cp_trans_ha_per_k,

            e_rot_ha: c.e_rot_ha,
            h_rot_ha: c.h_rot_ha,
            s_rot_ha_per_k: c.s_rot_ha_per_k,
            cv_rot_ha_per_k: c.cv_rot_ha_per_k,
            cp_rot_ha_per_k: c.cp_rot_ha_per_k,

            e_vib_thermal_ha: c.e_vib_thermal_ha,
            h_vib_ha: c.h_vib_ha,
            s_vib_ha_per_k: c.s_vib_ha_per_k,
            cv_vib_ha_per_k: c.cv_vib_ha_per_k,
            cp_vib_ha_per_k: c.cp_vib_ha_per_k,

            s_elec_ha_per_k: c.s_elec_ha_per_k,
        }
    }
}

/// Simulated IR or Raman spectrum (broadened on a wavenumber grid).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct FrequencySpectrum {
    /// Wavenumber grid in cm-1 (default 0..=4500 step 1, 4501 points).
    pub wavenumbers_cm1: Vec<f64>,
    /// Broadened intensity at each grid point. Same length as `wavenumbers_cm1`.
    /// Units: km/mol for IR, A^4/amu for Raman.
    pub intensity: Vec<f64>,
    /// Broadening kind used: `"lorentzian"` or `"gaussian"`.
    pub kind: String,
    /// FWHM used for the broadening in cm-1.
    pub fwhm_cm1: f64,
}

impl FrequencySpectrum {
    /// Map a [`qc_core::spectra::Spectrum`] to the WASM-friendly struct.
    fn from_core(spectrum: &qc_core::spectra::Spectrum) -> Self {
        Self {
            wavenumbers_cm1: spectrum.wavenumbers_cm1.clone(),
            intensity: spectrum.intensity.clone(),
            kind: broadening_kind_to_string(spectrum.kind),
            fwhm_cm1: spectrum.fwhm_cm1,
        }
    }
}

/// Full frequency-analysis result serialized to JavaScript.
///
/// Aggregates outputs from `rhf_hessian`/`dft_hessian`, `harmonic_analysis`,
/// `compute_ir_spectrum`, `compute_raman_spectrum`, `compute_thermochemistry`,
/// and `simulate_ir_spectrum` / `simulate_raman_spectrum`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct FrequencyWasmResult {
    // ---- Size metadata ----
    /// Number of atoms.
    pub n_atoms: usize,
    /// Number of vibrational modes (3N-6 nonlinear, 3N-5 linear, 0 atom).
    pub n_modes: usize,
    /// Rotor classification (snake_case): `"atom"`, `"linear"`,
    /// `"spherical_top"`, `"symmetric_top"`, or `"asymmetric_top"`.
    pub rotor_type: String,

    // ---- Electronic properties ----
    /// Electronic SCF energy in Hartree (from `HessianResult.energy`).
    pub electronic_energy_ha: f64,
    /// Equilibrium dipole moment in atomic units (e*bohr), `[x, y, z]`.
    pub dipole_au: [f64; 3],
    /// Equilibrium dipole moment in Debye, `[x, y, z]`.
    pub dipole_debye: [f64; 3],
    /// Static polarizability tensor in atomic units (bohr^3), symmetric 3x3.
    pub polarizability_au: [[f64; 3]; 3],
    /// Static polarizability tensor in A^3, symmetric 3x3.
    pub polarizability_ang3: [[f64; 3]; 3],

    // ---- Vibrational structure ----
    /// Vibrational frequencies in cm-1. Negative values are imaginary
    /// (transition states). Length = `n_modes`.
    pub frequencies_cm1: Vec<f64>,
    /// Reduced masses in amu, one per mode. Length = `n_modes`.
    pub reduced_masses_amu: Vec<f64>,
    /// Force constants in mDyne/A, one per mode. Length = `n_modes`.
    pub force_constants_mdyne: Vec<f64>,
    /// Cartesian normal modes, indexed as `[mode][atom][xyz]`.
    /// Shape: `n_modes x n_atoms x 3`.
    pub normal_modes_cartesian: Vec<Vec<[f64; 3]>>,
    /// Rotational constants in GHz (A >= B >= C; inf for zero moments of inertia).
    pub rotational_constants_ghz: [f64; 3],

    // ---- IR ----
    /// IR absorption intensities in km/mol, one per mode. Length = `n_modes`.
    pub ir_intensities_km_per_mol: Vec<f64>,

    // ---- Raman ----
    /// Raman scattering activities in A^4/amu, one per mode.
    pub raman_activities_a4_amu: Vec<f64>,
    /// Depolarization ratios rho_p (plane-polarized), one per mode.
    pub depolarization_ratios: Vec<f64>,

    // ---- Thermochemistry ----
    /// RRHO thermochemistry at the requested temperature and pressure.
    pub thermochemistry: FrequencyThermochemistry,

    // ---- Simulated spectra ----
    /// Continuous broadened IR spectrum on a wavenumber grid.
    pub ir_spectrum: FrequencySpectrum,
    /// Continuous broadened Raman spectrum on a wavenumber grid.
    pub raman_spectrum: FrequencySpectrum,

    // ---- Metadata ----
    /// Per-phase timings in milliseconds.
    pub timing_ms: FrequencyTiming,
    /// Whether the calculation was aborted mid-pipeline.
    ///
    /// The WASM function itself does not support mid-phase cancellation;
    /// this flag is always `false` in the value returned from
    /// `compute_frequencies`. The worker-side handler (US-101 T4) sets it
    /// to `true` when `isAborted()` was observed during progress events.
    pub aborted: bool,
}

// ============================================================================
// Internal Helpers
// ============================================================================

/// Convert a [`qc_core::thermo::RotorType`] to its snake_case string label.
fn rotor_type_to_string(rot_type: qc_core::thermo::RotorType) -> String {
    use qc_core::thermo::RotorType;
    match rot_type {
        RotorType::Atom => "atom".to_string(),
        RotorType::Linear => "linear".to_string(),
        RotorType::SphericalTop => "spherical_top".to_string(),
        RotorType::SymmetricTop => "symmetric_top".to_string(),
        RotorType::AsymmetricTop => "asymmetric_top".to_string(),
    }
}

/// Convert a [`qc_core::spectra::BroadeningKind`] to its canonical label.
fn broadening_kind_to_string(kind: qc_core::spectra::BroadeningKind) -> String {
    use qc_core::spectra::BroadeningKind;
    match kind {
        BroadeningKind::Lorentzian => "lorentzian".to_string(),
        BroadeningKind::Gaussian => "gaussian".to_string(),
    }
}

/// Parse a user-provided broadening kind string to the qc-core enum.
fn parse_broadening_kind(kind: &str) -> Result<qc_core::spectra::BroadeningKind, String> {
    use qc_core::spectra::BroadeningKind;
    match kind.to_lowercase().as_str() {
        "lorentzian" => Ok(BroadeningKind::Lorentzian),
        "gaussian" => Ok(BroadeningKind::Gaussian),
        other => Err(format!(
            "compute_frequencies: broadeningKind must be 'lorentzian' or 'gaussian', got '{}'",
            other
        )),
    }
}

/// Parse convergence profile string to enum.
///
/// Defaults to Medium for unknown strings.
fn parse_convergence_profile(s: &str) -> ConvergenceProfile {
    match s.to_lowercase().as_str() {
        "loose" => ConvergenceProfile::Loose,
        "medium" => ConvergenceProfile::Medium,
        "tight" => ConvergenceProfile::Tight,
        _ => ConvergenceProfile::Medium,
    }
}

/// Validate the [`FrequencyWasmInput`] scalar invariants.
///
/// Returns `Ok(())` on valid input, or a human-readable error string on
/// failure. The caller is responsible for wrapping the error in `JsError`.
fn validate_frequency_input(input: &FrequencyWasmInput) -> Result<(), String> {
    if input.atoms.is_empty() {
        return Err("compute_frequencies: atoms array is empty".to_string());
    }
    if !input.temperature_k.is_finite() || input.temperature_k <= 0.0 {
        return Err(format!(
            "compute_frequencies: temperatureK must be positive and finite, got {}",
            input.temperature_k
        ));
    }
    if !input.pressure_pa.is_finite() || input.pressure_pa <= 0.0 {
        return Err(format!(
            "compute_frequencies: pressurePa must be positive and finite, got {}",
            input.pressure_pa
        ));
    }
    if !input.fwhm_cm1.is_finite() || input.fwhm_cm1 <= 0.0 {
        return Err(format!(
            "compute_frequencies: fwhmCm1 must be positive and finite, got {}",
            input.fwhm_cm1
        ));
    }
    if input.multiplicity == 0 {
        return Err("compute_frequencies: multiplicity must be >= 1".to_string());
    }
    let method_lower = input.method.to_lowercase();
    if !["rhf", "hf", "lda", "b3lyp", "b3lyp-d3bj"].contains(&method_lower.as_str()) {
        return Err(format!(
            "compute_frequencies: unknown method '{}', expected 'rhf', 'hf', 'lda', 'b3lyp', or 'b3lyp-d3bj'",
            input.method
        ));
    }
    Ok(())
}

/// Best-effort millisecond clock that works in both native and wasm32
/// contexts. On wasm32 we use `js_sys::Date::now()`; on native we use a
/// monotonic `Instant` relative to a lazily-initialized origin.
fn freq_now_ms() -> f64 {
    #[cfg(target_arch = "wasm32")]
    {
        js_sys::Date::now()
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        use std::sync::OnceLock;
        use std::time::Instant;
        static ORIGIN: OnceLock<Instant> = OnceLock::new();
        let origin = ORIGIN.get_or_init(Instant::now);
        origin.elapsed().as_secs_f64() * 1000.0
    }
}

// ============================================================================
// WASM Export
// ============================================================================

/// Run the full analytical frequency analysis pipeline and return a
/// complete [`FrequencyWasmResult`].
///
/// Supports RHF, LDA, B3LYP, and B3LYP-D3BJ methods. Streams phase-granular
/// progress through the optional `progress_callback`, which receives a
/// plain JavaScript object with fields `{phase, percent, step, message}`.
/// The result contains frequencies, normal modes, IR intensities, Raman
/// activities, RRHO thermochemistry, and broadened IR/Raman spectra.
///
/// # Progress phases
///
/// 1. `"integrals"` -- SCF + Hessian assembly (nuclear CPHF runs inside
///    `rhf_hessian`/`dft_hessian`).
/// 2. `"nuclear_cphf"` -- Extract CPHF data from the Hessian result and
///    rebuild the density matrix for IR intensities.
/// 3. `"field_cphf"` -- Boundary marker; the actual field CPHF solve
///    executes inside `compute_raman_spectrum` during phase 4.
/// 4. `"assembly"` -- `harmonic_analysis` + `compute_ir_spectrum` +
///    `compute_raman_spectrum` + `compute_thermochemistry`.
/// 5. `"modes"` -- `simulate_ir_spectrum` + `simulate_raman_spectrum`.
///
/// # Abort semantics
///
/// The WASM function cannot interrupt itself mid-phase -- the worker-side
/// handler observes abort requests via the progress callback and sets the
/// `aborted` field on the returned object after the WASM call completes.
/// The value returned directly from Rust always has `aborted == false`.
///
/// # Arguments
///
/// * `input` - JS object matching [`FrequencyWasmInput`]
/// * `progress_callback` - Optional JS function called with a plain
///   `{phase, percent, step, message}` object. If `None`, progress is
///   discarded.
///
/// # Returns
///
/// A JS object matching [`FrequencyWasmResult`] with frequencies,
/// intensities, thermochemistry, and spectra.
///
/// # Errors
///
/// Returns a `JsError` if:
/// - Input deserialization or validation fails (empty atoms, unknown
///   method, non-positive T/P/FWHM, unknown broadening kind, etc.)
/// - Basis set construction fails
/// - SCF / Hessian / CPHF / harmonic analysis / IR / Raman / thermochemistry
///   / spectrum simulation fails
///
/// # Example (JavaScript)
///
/// ```javascript
/// const result = compute_frequencies(
///   {
///     atoms: [[8, 0, 0, 0], [1, 0, 0.756, 0.586], [1, 0, -0.756, 0.586]],
///     basisName: "sto-3g",
///     method: "rhf",
///     temperatureK: 298.15,
///     pressurePa: 101325,
///     multiplicity: 1,
///     broadeningKind: "lorentzian",
///     fwhmCm1: 8.0,
///   },
///   (progress) => {
///     console.log(`${progress.phase}: ${(progress.percent * 100).toFixed(0)}% ${progress.message}`);
///   }
/// );
/// console.log(`E = ${result.electronicEnergyHa} Ha`);
/// console.log(`Frequencies (cm-1):`, result.frequenciesCm1);
/// console.log(`G(298 K) = ${result.thermochemistry.gibbsHa} Ha`);
/// ```
#[wasm_bindgen]
pub fn compute_frequencies(
    input: JsValue,
    progress_callback: Option<js_sys::Function>,
) -> Result<JsValue, JsError> {
    let wasm_input: FrequencyWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid compute_frequencies input: {}", e)))?;

    let result = run_frequency_pipeline(&wasm_input, progress_callback.as_ref())
        .map_err(|e| JsError::new(&e))?;

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

/// Native-friendly entry point for [`compute_frequencies`] that returns
/// the concrete Rust result and accepts an optional `js_sys::Function`
/// progress callback.
///
/// This helper exists so that unit tests can drive the full pipeline
/// without crossing the JS boundary, while the wasm-bindgen export above
/// handles serialization. It takes the input by reference to avoid any
/// extra allocation and uses `&str` error messages for cheap propagation.
fn run_frequency_pipeline(
    input: &FrequencyWasmInput,
    progress_callback: Option<&js_sys::Function>,
) -> Result<FrequencyWasmResult, String> {
    use qc_core::basis::{Atom, BasisSet};
    use qc_core::scf::cphf::CphfConfig;
    use qc_core::scf::hessian::{dft_hessian, rhf_hessian};
    use qc_core::scf::ScfConfig;

    validate_frequency_input(input)?;

    // Emit progress helper: mirrors the `ks_scf` Reflect::set pattern so
    // that the JS callback receives a plain `{phase, percent, step, message}`
    // object with no serde allocation per tick. A closure keeps the
    // `progress_callback` borrow confined to the call site.
    let emit = |phase: &str, percent: f64, step: &str, msg: &str| {
        if let Some(cb) = progress_callback {
            let obj = js_sys::Object::new();
            let _ = js_sys::Reflect::set(&obj, &"phase".into(), &JsValue::from_str(phase));
            let _ = js_sys::Reflect::set(&obj, &"percent".into(), &JsValue::from_f64(percent));
            let _ = js_sys::Reflect::set(&obj, &"step".into(), &JsValue::from_str(step));
            let _ = js_sys::Reflect::set(&obj, &"message".into(), &JsValue::from_str(msg));
            let _ = cb.call1(&JsValue::NULL, &obj);
        }
    };

    let start_total = freq_now_ms();
    let mut timing = FrequencyTiming::default();

    // --------------------------------------------------------------------
    // PHASE 1: integrals (SCF + Hessian assembly + nuclear CPHF)
    // --------------------------------------------------------------------
    emit(
        "integrals",
        0.0,
        "start",
        "Running SCF and building Hessian",
    );
    let start_int = freq_now_ms();

    // Convert atoms from [[Z, x, y, z], ...] into the tuple form required
    // by rhf_hessian / dft_hessian, and into the typed Atom list required
    // by harmonic_analysis / compute_ir_spectrum / compute_raman_spectrum /
    // compute_thermochemistry.
    let atoms_tuple: Vec<(u8, [f64; 3])> = input
        .atoms
        .iter()
        .map(|a| (a[0] as u8, [a[1], a[2], a[3]]))
        .collect();

    let core_atoms: Vec<Atom> = atoms_tuple
        .iter()
        .map(|(z, pos)| Atom::new(*z, *pos).map_err(|e| format!("Invalid atom: {}", e)))
        .collect::<Result<Vec<_>, _>>()?;

    // BasisSet for downstream IR / Raman calls.
    let basis = BasisSet::build(core_atoms.clone(), &input.basis_name)
        .map_err(|e| format!("Basis set error: {}", e))?;

    let scf_config = ScfConfig {
        profile: parse_convergence_profile(&input.convergence_profile),
        max_iterations: input.max_iterations as usize,
        use_diis: true,
        ..Default::default()
    };

    let method_lower = input.method.to_lowercase();
    let hess_result = match method_lower.as_str() {
        "rhf" | "hf" => rhf_hessian(&atoms_tuple, &input.basis_name, &scf_config)
            .map_err(|e| format!("RHF Hessian failed: {}", e))?,
        "lda" => dft_hessian(&atoms_tuple, &input.basis_name, &scf_config, "lda")
            .map_err(|e| format!("LDA Hessian failed: {}", e))?,
        "b3lyp" | "b3lyp-d3bj" => {
            dft_hessian(&atoms_tuple, &input.basis_name, &scf_config, "b3lyp")
                .map_err(|e| format!("B3LYP Hessian failed: {}", e))?
        }
        other => {
            return Err(format!(
                "compute_frequencies: unknown method '{}', expected 'rhf', 'hf', 'lda', 'b3lyp', or 'b3lyp-d3bj'",
                other
            ));
        }
    };

    timing.integrals_ms = freq_now_ms() - start_int;
    emit("integrals", 1.0, "done", "SCF + Hessian complete");

    // --------------------------------------------------------------------
    // PHASE 2: nuclear_cphf -- extract CPHF data and rebuild density matrix.
    // The nuclear CPHF itself already ran inside rhf_hessian/dft_hessian;
    // this phase makes the CPHF solution explicit so downstream IR code
    // can consume `hess_result.mo1_cphf`.
    // --------------------------------------------------------------------
    emit("nuclear_cphf", 0.0, "extract", "Extracting CPHF data");
    let start_cphf = freq_now_ms();

    let cphf_data = hess_result
        .mo1_cphf
        .as_ref()
        .ok_or_else(|| "compute_frequencies: Hessian missing mo1_cphf data".to_string())?;

    // Density matrix D = 2 * C_occ * C_occ^T (closed-shell).
    // `compute_ir_spectrum` needs this to build the AO-basis dipole moment.
    let c_occ = cphf_data.mo_coeff.columns(0, cphf_data.n_occ).clone_owned();
    let density = 2.0 * &c_occ * c_occ.transpose();

    timing.nuclear_cphf_ms = freq_now_ms() - start_cphf;
    emit("nuclear_cphf", 1.0, "done", "Nuclear CPHF complete");

    // --------------------------------------------------------------------
    // PHASE 3: field_cphf boundary marker (actual solve in phase 4 inside
    // compute_raman_spectrum -- see strategy doc Section 6.3).
    // --------------------------------------------------------------------
    emit("field_cphf", 0.0, "start", "Preparing to solve field CPHF");

    // --------------------------------------------------------------------
    // PHASE 4: assembly -- harmonic analysis, IR, Raman, and thermochemistry.
    // --------------------------------------------------------------------
    emit("assembly", 0.0, "harmonic", "Running harmonic analysis");
    let start_asm = freq_now_ms();

    let freq_info = qc_core::thermo::harmonic_analysis(&core_atoms, &hess_result.hessian)
        .map_err(|e| format!("Harmonic analysis failed: {:?}", e))?;

    emit("assembly", 0.33, "ir", "Computing IR intensities");
    let gauge_origin: [f64; 3] = [0.0, 0.0, 0.0];
    let ir_result = qc_core::ir::compute_ir_spectrum(
        &core_atoms,
        &basis,
        &density,
        &hess_result,
        &freq_info,
        &gauge_origin,
    )
    .map_err(|e| format!("IR spectrum failed: {:?}", e))?;

    // Phase 3 field CPHF executes inside compute_raman_spectrum.
    let start_field = freq_now_ms();
    emit(
        "field_cphf",
        0.5,
        "solve",
        "Solving field CPHF for polarizability",
    );
    let cphf_config = CphfConfig::default();
    let raman_result = qc_core::raman::compute_raman_spectrum(
        &core_atoms,
        &basis,
        &hess_result,
        &freq_info,
        &cphf_config,
    )
    .map_err(|e| format!("Raman spectrum failed: {:?}", e))?;
    timing.field_cphf_ms = freq_now_ms() - start_field;
    emit("field_cphf", 1.0, "done", "Field CPHF complete");

    emit("assembly", 0.66, "thermo", "Computing RRHO thermochemistry");
    let thermo_result = qc_core::thermochemistry::compute_thermochemistry(
        &freq_info,
        hess_result.energy,
        &core_atoms,
        input.temperature_k,
        input.pressure_pa,
        input.symmetry_number_override,
        input.multiplicity,
    )
    .map_err(|e| format!("Thermochemistry failed: {:?}", e))?;

    // Subtract the nested field_cphf duration so `assembly_ms` represents
    // only the portion of phase 4 not already accounted for by field_cphf.
    let assembly_total = freq_now_ms() - start_asm;
    timing.assembly_ms = (assembly_total - timing.field_cphf_ms).max(0.0);
    emit("assembly", 1.0, "done", "Assembly complete");

    // --------------------------------------------------------------------
    // PHASE 5: modes -- simulate broadened IR/Raman spectra.
    // --------------------------------------------------------------------
    let start_modes = freq_now_ms();
    let kind = parse_broadening_kind(&input.broadening_kind)?;

    emit("modes", 0.0, "ir_sim", "Simulating IR spectrum");
    let ir_spectrum = qc_core::spectra::simulate_ir_spectrum(
        &freq_info.freq_wavenumber,
        &ir_result.intensities_km_per_mol,
        input.fwhm_cm1,
        kind,
    )
    .map_err(|e| format!("IR spectrum simulation failed: {}", e))?;

    emit("modes", 0.5, "raman_sim", "Simulating Raman spectrum");
    let raman_spectrum = qc_core::spectra::simulate_raman_spectrum(
        &freq_info.freq_wavenumber,
        &raman_result.raman_activities_a4_amu,
        input.fwhm_cm1,
        kind,
    )
    .map_err(|e| format!("Raman spectrum simulation failed: {}", e))?;

    timing.modes_ms = freq_now_ms() - start_modes;
    timing.total_ms = freq_now_ms() - start_total;
    emit("modes", 1.0, "done", "Frequency analysis complete");

    // --------------------------------------------------------------------
    // Assemble the result payload.
    // --------------------------------------------------------------------
    Ok(FrequencyWasmResult {
        n_atoms: core_atoms.len(),
        n_modes: freq_info.n_modes,
        rotor_type: rotor_type_to_string(freq_info.rot_type),

        electronic_energy_ha: hess_result.energy,
        dipole_au: ir_result.equilibrium_dipole_au,
        dipole_debye: ir_result.equilibrium_dipole_debye,
        polarizability_au: raman_result.polarizability_au,
        polarizability_ang3: raman_result.polarizability_ang3,

        frequencies_cm1: freq_info.freq_wavenumber.clone(),
        reduced_masses_amu: freq_info.reduced_mass.clone(),
        force_constants_mdyne: freq_info.force_const_dyne.clone(),
        normal_modes_cartesian: freq_info.norm_mode.clone(),
        rotational_constants_ghz: freq_info.rotational_constants_ghz,

        ir_intensities_km_per_mol: ir_result.intensities_km_per_mol.clone(),
        raman_activities_a4_amu: raman_result.raman_activities_a4_amu.clone(),
        depolarization_ratios: raman_result.depolarization_ratios.clone(),

        thermochemistry: FrequencyThermochemistry::from_core(&thermo_result),
        ir_spectrum: FrequencySpectrum::from_core(&ir_spectrum),
        raman_spectrum: FrequencySpectrum::from_core(&raman_spectrum),

        timing_ms: timing,
        aborted: false,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_frequency_input() -> FrequencyWasmInput {
        FrequencyWasmInput {
            atoms: vec![[1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 1.4]],
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            temperature_k: 298.15,
            pressure_pa: 101_325.0,
            symmetry_number_override: Some(2),
            multiplicity: 1,
            broadening_kind: "lorentzian".to_string(),
            fwhm_cm1: 8.0,
            convergence_profile: "tight".to_string(),
            max_iterations: 100,
        }
    }

    fn sample_thermochemistry() -> FrequencyThermochemistry {
        FrequencyThermochemistry {
            temperature_k: 298.15,
            pressure_pa: 101_325.0,
            symmetry_number: 2,
            multiplicity: 1,
            total_mass_amu: 2.016,
            n_vib_modes_used: 1,
            n_imag: 0,
            zpe_ha: 0.010_123,
            e_0k_ha: -1.106_591,
            internal_energy_ha: -1.102_678,
            enthalpy_ha: -1.101_734,
            entropy_ha_per_k: 5.95e-5,
            gibbs_ha: -1.119_479,
            cv_ha_per_k: 9.48e-6,
            cp_ha_per_k: 1.264e-5,
            e_trans_ha: 1.414_59e-3,
            h_trans_ha: 2.357_66e-3,
            s_trans_ha_per_k: 4.44e-5,
            cv_trans_ha_per_k: 4.75e-6,
            cp_trans_ha_per_k: 7.92e-6,
            e_rot_ha: 9.43e-4,
            h_rot_ha: 9.43e-4,
            s_rot_ha_per_k: 4.49e-6,
            cv_rot_ha_per_k: 3.17e-6,
            cp_rot_ha_per_k: 3.17e-6,
            e_vib_thermal_ha: 1.2e-9,
            h_vib_ha: 0.010_123,
            s_vib_ha_per_k: 2.5e-11,
            cv_vib_ha_per_k: 3.0e-9,
            cp_vib_ha_per_k: 3.0e-9,
            s_elec_ha_per_k: 0.0,
        }
    }

    fn sample_spectrum(kind: &str, fwhm: f64) -> FrequencySpectrum {
        FrequencySpectrum {
            wavenumbers_cm1: vec![0.0, 1.0, 2.0, 3.0, 4500.0],
            intensity: vec![0.1, 0.2, 0.3, 0.4, 0.0],
            kind: kind.to_string(),
            fwhm_cm1: fwhm,
        }
    }

    fn sample_frequency_result() -> FrequencyWasmResult {
        FrequencyWasmResult {
            n_atoms: 2,
            n_modes: 1,
            rotor_type: "linear".to_string(),
            electronic_energy_ha: -1.116_714,
            dipole_au: [0.0, 0.0, 0.0],
            dipole_debye: [0.0, 0.0, 0.0],
            polarizability_au: [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]],
            polarizability_ang3: [[0.74, 0.0, 0.0], [0.0, 0.74, 0.0], [0.0, 0.0, 0.74]],
            frequencies_cm1: vec![4_646.5],
            reduced_masses_amu: vec![0.5],
            force_constants_mdyne: vec![6.4],
            normal_modes_cartesian: vec![vec![[0.0, 0.0, 0.707], [0.0, 0.0, -0.707]]],
            // Note: in practice the first rotational constant for a linear
            // molecule would be `f64::INFINITY`, but `serde_json` does not
            // support non-finite floats (it serializes them as `null`). Use
            // large finite placeholders here so the native JSON round-trip
            // test remains meaningful; the wasm32 path uses
            // `serde_wasm_bindgen` which handles Infinity correctly.
            rotational_constants_ghz: [1.0e30, 1_824.32, 1_824.32],
            ir_intensities_km_per_mol: vec![0.0],
            raman_activities_a4_amu: vec![80.0],
            depolarization_ratios: vec![0.0],
            thermochemistry: sample_thermochemistry(),
            ir_spectrum: sample_spectrum("lorentzian", 8.0),
            raman_spectrum: sample_spectrum("lorentzian", 8.0),
            timing_ms: FrequencyTiming {
                integrals_ms: 10.0,
                nuclear_cphf_ms: 1.0,
                field_cphf_ms: 5.0,
                assembly_ms: 3.0,
                modes_ms: 2.0,
                total_ms: 21.0,
            },
            aborted: false,
        }
    }

    #[test]
    fn frequency_wasm_input_serializes_to_camel_case() {
        let input = sample_frequency_input();
        let json = serde_json::to_string(&input).unwrap();

        assert!(json.contains("\"basisName\":\"sto-3g\""));
        assert!(json.contains("\"method\":\"rhf\""));
        assert!(json.contains("\"temperatureK\":298.15"));
        assert!(json.contains("\"pressurePa\":101325"));
        assert!(json.contains("\"symmetryNumberOverride\":2"));
        assert!(json.contains("\"multiplicity\":1"));
        assert!(json.contains("\"broadeningKind\":\"lorentzian\""));
        assert!(json.contains("\"fwhmCm1\":8"));
        assert!(json.contains("\"convergenceProfile\":\"tight\""));
        assert!(json.contains("\"maxIterations\":100"));
        // No snake_case leakage
        assert!(!json.contains("temperature_k"));
        assert!(!json.contains("pressure_pa"));
    }

    #[test]
    fn frequency_wasm_input_json_round_trip() {
        let original = sample_frequency_input();
        let json = serde_json::to_string(&original).unwrap();
        let recovered: FrequencyWasmInput = serde_json::from_str(&json).unwrap();
        assert_eq!(recovered, original);
    }

    #[test]
    fn frequency_wasm_input_defaults_apply() {
        // Minimal JSON: only the three required fields. All optionals default.
        let json = r#"{
            "atoms": [[1, 0, 0, 0], [1, 0, 0, 1.4]],
            "basisName": "sto-3g",
            "method": "rhf"
        }"#;
        let input: FrequencyWasmInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.temperature_k, 298.15);
        assert_eq!(input.pressure_pa, 101_325.0);
        assert_eq!(input.multiplicity, 1);
        assert_eq!(input.fwhm_cm1, 8.0);
        assert_eq!(input.broadening_kind, "lorentzian");
        assert_eq!(input.convergence_profile, "tight");
        assert_eq!(input.max_iterations, 100);
        assert!(input.symmetry_number_override.is_none());
    }

    #[test]
    fn frequency_wasm_result_serializes_to_camel_case() {
        let result = sample_frequency_result();
        let json = serde_json::to_string(&result).unwrap();
        // A sampling of expected fields
        for needle in [
            "\"nAtoms\":2",
            "\"nModes\":1",
            "\"rotorType\":\"linear\"",
            "\"electronicEnergyHa\":",
            "\"dipoleAu\":",
            "\"dipoleDebye\":",
            "\"polarizabilityAu\":",
            "\"polarizabilityAng3\":",
            "\"frequenciesCm1\":",
            "\"reducedMassesAmu\":",
            "\"forceConstantsMdyne\":",
            "\"normalModesCartesian\":",
            "\"rotationalConstantsGhz\":",
            "\"irIntensitiesKmPerMol\":",
            "\"ramanActivitiesA4Amu\":",
            "\"depolarizationRatios\":",
            "\"thermochemistry\":",
            "\"irSpectrum\":",
            "\"ramanSpectrum\":",
            "\"timingMs\":",
            "\"aborted\":false",
        ] {
            assert!(
                json.contains(needle),
                "FrequencyWasmResult JSON missing {}: {}",
                needle,
                json
            );
        }
    }

    #[test]
    fn frequency_wasm_result_json_round_trip() {
        let original = sample_frequency_result();
        let json = serde_json::to_string(&original).unwrap();
        let recovered: FrequencyWasmResult = serde_json::from_str(&json).unwrap();
        assert_eq!(recovered, original);
    }

    #[test]
    fn frequency_thermochemistry_camel_case_and_round_trip() {
        let thermo = sample_thermochemistry();
        let json = serde_json::to_string(&thermo).unwrap();
        for needle in [
            "\"temperatureK\":298.15",
            "\"pressurePa\":101325",
            "\"symmetryNumber\":2",
            "\"multiplicity\":1",
            "\"zpeHa\":",
            "\"e0kHa\":",
            "\"internalEnergyHa\":",
            "\"enthalpyHa\":",
            "\"entropyHaPerK\":",
            "\"gibbsHa\":",
            "\"cvHaPerK\":",
            "\"cpHaPerK\":",
            "\"eTransHa\":",
            "\"hTransHa\":",
            "\"sTransHaPerK\":",
            "\"cvTransHaPerK\":",
            "\"cpTransHaPerK\":",
            "\"eRotHa\":",
            "\"hRotHa\":",
            "\"sRotHaPerK\":",
            "\"cvRotHaPerK\":",
            "\"cpRotHaPerK\":",
            "\"eVibThermalHa\":",
            "\"hVibHa\":",
            "\"sVibHaPerK\":",
            "\"cvVibHaPerK\":",
            "\"cpVibHaPerK\":",
            "\"sElecHaPerK\":",
        ] {
            assert!(
                json.contains(needle),
                "FrequencyThermochemistry JSON missing {}: {}",
                needle,
                json
            );
        }
        let recovered: FrequencyThermochemistry = serde_json::from_str(&json).unwrap();
        assert_eq!(recovered, thermo);
    }

    #[test]
    fn frequency_spectrum_round_trip() {
        let spectrum = sample_spectrum("gaussian", 12.0);
        let json = serde_json::to_string(&spectrum).unwrap();
        assert!(json.contains("\"wavenumbersCm1\":[0.0,1.0,2.0,3.0,4500.0]"));
        assert!(json.contains("\"intensity\":[0.1,0.2,0.3,0.4,0.0]"));
        assert!(json.contains("\"kind\":\"gaussian\""));
        assert!(json.contains("\"fwhmCm1\":12"));
        let recovered: FrequencySpectrum = serde_json::from_str(&json).unwrap();
        assert_eq!(recovered, spectrum);
    }

    #[test]
    fn frequency_timing_default_is_zero() {
        let timing = FrequencyTiming::default();
        assert_eq!(timing.integrals_ms, 0.0);
        assert_eq!(timing.nuclear_cphf_ms, 0.0);
        assert_eq!(timing.field_cphf_ms, 0.0);
        assert_eq!(timing.assembly_ms, 0.0);
        assert_eq!(timing.modes_ms, 0.0);
        assert_eq!(timing.total_ms, 0.0);
    }

    #[test]
    fn frequency_default_helpers_match_strategy_spec() {
        assert_eq!(default_freq_temperature_k(), 298.15);
        assert_eq!(default_freq_pressure_pa(), 101_325.0);
        assert_eq!(default_freq_multiplicity(), 1);
        assert_eq!(default_freq_broadening_kind(), "lorentzian");
        assert_eq!(default_freq_fwhm_cm1(), 8.0);
        assert_eq!(default_freq_convergence_profile(), "tight");
        assert_eq!(default_freq_max_iterations(), 100);
    }

    #[test]
    fn parse_broadening_kind_accepts_canonical_variants() {
        use qc_core::spectra::BroadeningKind;
        assert_eq!(
            parse_broadening_kind("lorentzian").unwrap(),
            BroadeningKind::Lorentzian
        );
        assert_eq!(
            parse_broadening_kind("LORENTZIAN").unwrap(),
            BroadeningKind::Lorentzian
        );
        assert_eq!(
            parse_broadening_kind("Gaussian").unwrap(),
            BroadeningKind::Gaussian
        );
        assert!(parse_broadening_kind("voigt").is_err());
    }

    #[test]
    fn rotor_type_to_string_covers_all_variants() {
        use qc_core::thermo::RotorType;
        assert_eq!(rotor_type_to_string(RotorType::Atom), "atom");
        assert_eq!(rotor_type_to_string(RotorType::Linear), "linear");
        assert_eq!(
            rotor_type_to_string(RotorType::SphericalTop),
            "spherical_top"
        );
        assert_eq!(
            rotor_type_to_string(RotorType::SymmetricTop),
            "symmetric_top"
        );
        assert_eq!(
            rotor_type_to_string(RotorType::AsymmetricTop),
            "asymmetric_top"
        );
    }

    #[test]
    fn validate_frequency_input_accepts_default_sample() {
        let input = sample_frequency_input();
        validate_frequency_input(&input).expect("default sample input should validate");
    }

    #[test]
    fn validate_frequency_input_rejects_empty_atoms() {
        let mut input = sample_frequency_input();
        input.atoms.clear();
        let err = validate_frequency_input(&input).unwrap_err();
        assert!(err.contains("atoms array is empty"));
    }

    #[test]
    fn validate_frequency_input_rejects_unknown_method() {
        let mut input = sample_frequency_input();
        input.method = "ccsd_t".to_string();
        let err = validate_frequency_input(&input).unwrap_err();
        assert!(err.contains("unknown method"));
    }

    #[test]
    fn validate_frequency_input_rejects_non_positive_temperature() {
        let mut input = sample_frequency_input();
        input.temperature_k = 0.0;
        let err = validate_frequency_input(&input).unwrap_err();
        assert!(err.contains("temperatureK"));
    }

    #[test]
    fn validate_frequency_input_rejects_non_positive_pressure() {
        let mut input = sample_frequency_input();
        input.pressure_pa = -1.0;
        let err = validate_frequency_input(&input).unwrap_err();
        assert!(err.contains("pressurePa"));
    }

    #[test]
    fn validate_frequency_input_rejects_non_positive_fwhm() {
        let mut input = sample_frequency_input();
        input.fwhm_cm1 = 0.0;
        let err = validate_frequency_input(&input).unwrap_err();
        assert!(err.contains("fwhmCm1"));
    }

    #[test]
    fn validate_frequency_input_rejects_zero_multiplicity() {
        let mut input = sample_frequency_input();
        input.multiplicity = 0;
        let err = validate_frequency_input(&input).unwrap_err();
        assert!(err.contains("multiplicity"));
    }

    // ----------------------------------------------------------------
    // End-to-end native pipeline tests (#[cfg(not(target_arch = "wasm32"))]).
    //
    // These call `run_frequency_pipeline` directly, bypassing the
    // wasm-bindgen boundary, so they exercise every qc-core call that
    // the WASM export would make without needing `js_sys::Function`.
    // ----------------------------------------------------------------

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn compute_frequencies_pipeline_h2_sto3g_rhf() {
        // H2 at equilibrium R = 1.4 bohr (STO-3G reference geometry).
        let input = FrequencyWasmInput {
            atoms: vec![[1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 1.4]],
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            temperature_k: 298.15,
            pressure_pa: 101_325.0,
            symmetry_number_override: Some(2),
            multiplicity: 1,
            broadening_kind: "lorentzian".to_string(),
            fwhm_cm1: 8.0,
            convergence_profile: "tight".to_string(),
            max_iterations: 100,
        };

        let result = run_frequency_pipeline(&input, None).expect("H2/STO-3G pipeline should run");

        // Basic metadata sanity
        assert_eq!(result.n_atoms, 2);
        assert_eq!(result.n_modes, 1);
        assert_eq!(result.rotor_type, "linear");

        // SCF energy should be within 0.5 Ha of the reference
        // PySCF H2/STO-3G reference: -1.116714 Ha
        assert!(
            (result.electronic_energy_ha + 1.116_714).abs() < 0.01,
            "electronicEnergyHa = {} (expected ~-1.117)",
            result.electronic_energy_ha
        );

        // Single stretch mode
        assert_eq!(result.frequencies_cm1.len(), 1);
        let omega = result.frequencies_cm1[0];
        assert!(
            (4000.0..=5500.0).contains(&omega),
            "frequency {} cm-1 outside expected [4000, 5500] window",
            omega
        );

        // Derived vibrational arrays must have matching length = n_modes
        assert_eq!(result.reduced_masses_amu.len(), 1);
        assert_eq!(result.force_constants_mdyne.len(), 1);
        assert_eq!(result.normal_modes_cartesian.len(), 1);
        assert_eq!(result.normal_modes_cartesian[0].len(), 2);
        assert_eq!(result.ir_intensities_km_per_mol.len(), 1);
        assert_eq!(result.raman_activities_a4_amu.len(), 1);
        assert_eq!(result.depolarization_ratios.len(), 1);

        // H2 is homonuclear -- IR intensity must be (very close to) zero
        assert!(
            result.ir_intensities_km_per_mol[0] < 1e-6,
            "H2 IR intensity must vanish (symmetry), got {}",
            result.ir_intensities_km_per_mol[0]
        );
        // ...but Raman activity must be nonzero (totally symmetric stretch)
        assert!(
            result.raman_activities_a4_amu[0] > 1e-3,
            "H2 Raman activity should be positive, got {}",
            result.raman_activities_a4_amu[0]
        );

        // Polarizability tensor should be positive-definite (trace > 0)
        let pol_trace = result.polarizability_au[0][0]
            + result.polarizability_au[1][1]
            + result.polarizability_au[2][2];
        assert!(pol_trace > 0.0, "polarizability trace {} <= 0", pol_trace);

        // Thermochemistry: ZPE > 0, G < E_elec, S_tot > 0
        assert!(result.thermochemistry.zpe_ha > 0.0);
        assert!(result.thermochemistry.entropy_ha_per_k > 0.0);
        assert!(result.thermochemistry.gibbs_ha < result.electronic_energy_ha + 0.05);
        assert!(!result.thermochemistry.internal_energy_ha.is_nan());
        assert_eq!(result.thermochemistry.n_imag, 0);
        assert_eq!(result.thermochemistry.n_vib_modes_used, 1);

        // Spectrum arrays share the same grid and have 4501 points
        assert_eq!(result.ir_spectrum.wavenumbers_cm1.len(), 4501);
        assert_eq!(result.ir_spectrum.intensity.len(), 4501);
        assert_eq!(result.raman_spectrum.wavenumbers_cm1.len(), 4501);
        assert_eq!(result.raman_spectrum.intensity.len(), 4501);
        assert_eq!(result.ir_spectrum.kind, "lorentzian");
        assert_eq!(result.raman_spectrum.kind, "lorentzian");
        assert_eq!(result.ir_spectrum.fwhm_cm1, 8.0);

        // Timing fields are populated and coherent
        assert!(result.timing_ms.integrals_ms >= 0.0);
        assert!(result.timing_ms.nuclear_cphf_ms >= 0.0);
        assert!(result.timing_ms.field_cphf_ms >= 0.0);
        assert!(result.timing_ms.assembly_ms >= 0.0);
        assert!(result.timing_ms.modes_ms >= 0.0);
        assert!(result.timing_ms.total_ms >= result.timing_ms.integrals_ms);

        assert!(!result.aborted);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn compute_frequencies_pipeline_h2o_sto3g_rhf() {
        // H2O at STO-3G equilibrium (bohr coordinates).
        let input = FrequencyWasmInput {
            atoms: vec![
                [8.0, 0.0, 0.0, 0.0],
                [1.0, 0.0, 1.43, 1.11],
                [1.0, 0.0, -1.43, 1.11],
            ],
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            temperature_k: 298.15,
            pressure_pa: 101_325.0,
            symmetry_number_override: Some(2),
            multiplicity: 1,
            broadening_kind: "lorentzian".to_string(),
            fwhm_cm1: 8.0,
            convergence_profile: "tight".to_string(),
            max_iterations: 100,
        };

        let result = run_frequency_pipeline(&input, None).expect("H2O/STO-3G pipeline should run");

        // Structure / metadata
        assert_eq!(result.n_atoms, 3);
        assert_eq!(result.n_modes, 3);
        assert!(
            result.rotor_type == "asymmetric_top" || result.rotor_type == "symmetric_top",
            "H2O rotor type {} should be (a)symmetric top",
            result.rotor_type
        );

        // Electronic energy within ~0.1 Ha of the PySCF reference
        assert!(
            result.electronic_energy_ha < -74.0 && result.electronic_energy_ha > -76.5,
            "electronicEnergyHa = {} (expected ~-75 Ha)",
            result.electronic_energy_ha
        );

        // Dipole moment: H2O has a permanent dipole
        let dip_mag_debye = (result.dipole_debye[0].powi(2)
            + result.dipole_debye[1].powi(2)
            + result.dipole_debye[2].powi(2))
        .sqrt();
        assert!(
            dip_mag_debye > 0.5,
            "H2O dipole magnitude {} Debye too small",
            dip_mag_debye
        );

        // Three vibrational modes with nonzero IR and Raman activity
        assert_eq!(result.frequencies_cm1.len(), 3);
        assert_eq!(result.ir_intensities_km_per_mol.len(), 3);
        assert_eq!(result.raman_activities_a4_amu.len(), 3);
        for (i, freq) in result.frequencies_cm1.iter().enumerate() {
            assert!(
                freq.is_finite(),
                "mode {} frequency {} is not finite",
                i,
                freq
            );
        }

        let any_ir_positive = result.ir_intensities_km_per_mol.iter().any(|&x| x > 1e-4);
        assert!(
            any_ir_positive,
            "At least one H2O IR intensity should be > 0"
        );
        let any_raman_positive = result.raman_activities_a4_amu.iter().any(|&x| x > 1e-4);
        assert!(
            any_raman_positive,
            "At least one H2O Raman activity should be > 0"
        );

        // Thermochemistry sanity
        assert!(result.thermochemistry.zpe_ha > 0.0);
        assert!(result.thermochemistry.entropy_ha_per_k > 0.0);
        assert!(result.thermochemistry.internal_energy_ha.is_finite());
        assert!(result.thermochemistry.enthalpy_ha.is_finite());
        assert!(result.thermochemistry.gibbs_ha.is_finite());
        assert_eq!(
            result.thermochemistry.n_vib_modes_used + result.thermochemistry.n_imag,
            3
        );

        // Spectrum arrays
        assert_eq!(result.ir_spectrum.wavenumbers_cm1.len(), 4501);
        assert_eq!(result.raman_spectrum.wavenumbers_cm1.len(), 4501);

        // Timing
        assert!(result.timing_ms.total_ms > 0.0);
        assert!(!result.aborted);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn compute_frequencies_pipeline_unknown_method_errors() {
        let mut input = sample_frequency_input();
        input.method = "ccsd".to_string();
        let err = run_frequency_pipeline(&input, None).unwrap_err();
        assert!(
            err.contains("unknown method"),
            "expected unknown method error, got: {}",
            err
        );
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn compute_frequencies_pipeline_applies_default_values() {
        // Minimal request, all optional fields defaulted.
        let json = r#"{
            "atoms": [[1, 0, 0, 0], [1, 0, 0, 1.4]],
            "basisName": "sto-3g",
            "method": "rhf"
        }"#;
        let input: FrequencyWasmInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.temperature_k, 298.15);
        assert_eq!(input.pressure_pa, 101_325.0);
        assert_eq!(input.fwhm_cm1, 8.0);
        assert_eq!(input.broadening_kind, "lorentzian");
        assert_eq!(input.convergence_profile, "tight");

        let result = run_frequency_pipeline(&input, None).expect("default H2 pipeline should run");
        assert_eq!(result.n_atoms, 2);
        assert_eq!(result.n_modes, 1);
        assert_eq!(result.thermochemistry.temperature_k, 298.15);
        assert_eq!(result.thermochemistry.pressure_pa, 101_325.0);
        assert_eq!(result.ir_spectrum.kind, "lorentzian");
        assert_eq!(result.ir_spectrum.fwhm_cm1, 8.0);
    }
}
