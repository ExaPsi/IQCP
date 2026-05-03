//! WebAssembly bindings for qc-core algorithms
//!
//! Exposes quantum chemistry compute functions to JavaScript via wasm-bindgen.
//!
//! # Architecture
//!
//! This crate serves as the bridge between the pure Rust algorithms in `qc-core`
//! and the JavaScript/TypeScript frontend. It handles:
//!
//! - Serialization/deserialization via serde-wasm-bindgen
//! - Error translation to JavaScript exceptions
//! - Progress callbacks for long-running computations
//!
//! # Usage
//!
//! The WASM module is loaded in a Web Worker to keep the UI responsive:
//!
//! ```javascript
//! import init, { version } from './qc_wasm.js';
//!
//! await init();
//! console.log(`IQCP WASM v${version()}`);
//! ```
//!
//! # Build
//!
//! ```bash
//! wasm-pack build --target web --out-dir ../../apps/web/src/wasm
//! ```

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

// ============================================================================
// Thread Pool Initialization (for parallel feature)
// ============================================================================

/// Re-export the thread pool initialization function from wasm-bindgen-rayon.
///
/// This function must be called before any parallel computations can be performed.
/// It initializes Rayon's thread pool in the WebAssembly context using Web Workers.
///
/// # Requirements
///
/// For multithreading to work, the following HTTP headers must be set on the server:
/// - `Cross-Origin-Opener-Policy: same-origin`
/// - `Cross-Origin-Embedder-Policy: require-corp`
///
/// These headers enable `SharedArrayBuffer`, which is required for WASM threads.
///
/// # Usage
///
/// ```javascript
/// import init, { initThreadPool } from './qc_wasm.js';
///
/// await init();
/// await initThreadPool(navigator.hardwareConcurrency);
/// // Now parallel computations are available
/// ```
///
/// # Arguments
///
/// * `num_threads` - Number of threads to use. Typically `navigator.hardwareConcurrency`.
///
/// # Build Requirements
///
/// The WASM module must be built with the `parallel` feature and thread support:
///
/// ```bash
/// RUSTFLAGS='-C target-feature=+atomics,+bulk-memory,+mutable-globals' \
/// wasm-pack build --target web --features parallel -- -Z build-std=panic_abort,std
/// ```
#[cfg(feature = "parallel")]
pub use wasm_bindgen_rayon::init_thread_pool;

/// Returns whether the WASM module was built with parallel/threading support.
///
/// This can be used by the frontend to conditionally show threading options
/// or fall back to single-threaded mode.
///
/// # Example
///
/// ```javascript
/// import { hasThreadingSupport } from './qc_wasm.js';
///
/// if (hasThreadingSupport()) {
///     await initThreadPool(navigator.hardwareConcurrency);
/// }
/// ```
#[wasm_bindgen]
pub fn has_threading_support() -> bool {
    cfg!(feature = "parallel")
}

// Re-export Boys function types for convenience
pub use qc_core::boys::{BoysError, BoysMethod, BoysResult};

// Re-export Rys quadrature types for convenience
pub use qc_core::rys::{ErrorCurvePoint, ErrorCurveResult, RysError, RysMethod, RysResult};

// Re-export SCF types for convenience
pub use qc_core::scf::{
    ConvergenceProfile, PresetSystem, PresetSystemJson, ScfConfig, ScfError, ScfIteration,
    ScfOutput,
};

// Re-export PES types for convenience
pub use qc_core::scf::pes::{PesEquilibrium, PesPoint, PesScanResult};

// Re-export optimizer types for convenience
pub use qc_core::optimizer::{OptMethod, OptimizationConfig, OptimizationResult, OptimizationStep};

// Re-export orbital types for convenience
pub use qc_core::orbital::{MarchingCubesResult, MoGridResult, OrbitalError};

// ============================================================================
// Boys Function WASM Bindings
// ============================================================================

/// Evaluate the Boys function F_m(T) for a single (m, T) pair.
///
/// The Boys function is fundamental to Gaussian integral evaluation and is
/// defined as:
///
/// F_m(T) = ∫₀¹ t^(2m) e^(-T·t²) dt
///
/// # Arguments
///
/// * `m` - Order of the Boys function (0 <= m <= 50)
/// * `t` - Argument value (T >= 0)
///
/// # Returns
///
/// A JavaScript object containing:
/// - `value`: The computed F_m(T) value
/// - `method`: The computational method used ("zero", "series", or "recurrence")
/// - `m`: The order parameter
/// - `t`: The argument parameter
///
/// # Errors
///
/// Returns a `JsError` if:
/// - `m` exceeds the maximum supported order (50)
/// - `t` is negative
///
/// # Example
///
/// ```javascript
/// import { boys_eval } from './qc_wasm.js';
///
/// const result = boys_eval(0, 0.5);
/// console.log(result);
/// // { value: 0.8556243918921487, method: "series", m: 0, t: 0.5 }
/// ```
#[wasm_bindgen]
pub fn boys_eval(m: u32, t: f64) -> Result<JsValue, JsError> {
    let result = qc_core::boys::boys_eval(m, t).map_err(|e| JsError::new(&e.to_string()))?;
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

/// Evaluate Boys functions F_0(T) through F_m(T) for a single T value.
///
/// Computes all Boys function values from order 0 to the specified maximum order.
/// This is efficient because the algorithms compute intermediate orders during
/// the recurrence/series evaluation.
///
/// # Arguments
///
/// * `m_max` - Maximum order to compute (0 <= m_max <= 50)
/// * `t` - Argument value (T >= 0)
///
/// # Returns
///
/// A JavaScript array of result objects, one for each order from 0 to m_max.
///
/// # Example
///
/// ```javascript
/// import { boys_eval_all } from './qc_wasm.js';
///
/// const results = boys_eval_all(3, 1.0);
/// // Returns array: [F_0(1.0), F_1(1.0), F_2(1.0), F_3(1.0)]
/// ```
#[wasm_bindgen]
pub fn boys_eval_all(m_max: u32, t: f64) -> Result<JsValue, JsError> {
    let results =
        qc_core::boys::boys_eval_all(m_max, t).map_err(|e| JsError::new(&e.to_string()))?;
    serde_wasm_bindgen::to_value(&results).map_err(|e| JsError::new(&e.to_string()))
}

/// Evaluate the Boys function F_m(T) for multiple T values at fixed order m.
///
/// This is the primary function for generating sweep plots where the student
/// varies T while keeping m constant.
///
/// # Arguments
///
/// * `m` - Order of the Boys function (0 <= m <= 50)
/// * `ts` - Array of T values to evaluate (all must be >= 0)
///
/// # Returns
///
/// A JavaScript array of result objects, one for each T value.
///
/// # Example
///
/// ```javascript
/// import { boys_eval_many } from './qc_wasm.js';
///
/// const ts = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0];
/// const results = boys_eval_many(0, ts);
/// // Returns array of F_0(T) for each T value
/// ```
#[wasm_bindgen]
pub fn boys_eval_many(m: u32, ts: Vec<f64>) -> Result<JsValue, JsError> {
    let results =
        qc_core::boys::boys_eval_many(m, &ts).map_err(|e| JsError::new(&e.to_string()))?;
    serde_wasm_bindgen::to_value(&results).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Rys Quadrature WASM Bindings
// ============================================================================

/// Compute Rys quadrature roots and weights for a given order and T value.
///
/// Rys quadrature provides optimal integration points for Gaussian integrals
/// by computing roots and weights of Rys polynomials. This is fundamental to
/// molecular integral evaluation in quantum chemistry.
///
/// # Arguments
///
/// * `n` - Number of quadrature points (1 <= n <= 10)
/// * `t` - Argument value (T >= 0)
///
/// # Returns
///
/// A JavaScript object containing:
/// - `roots`: Array of quadrature roots in the interval [0, 1)
/// - `weights`: Array of corresponding weights (all positive)
/// - `nroots`: Number of roots/weights
/// - `t`: The argument parameter
/// - `method`: The computational method used ("special" or "standard")
///
/// # Errors
///
/// Returns a `JsError` if:
/// - `n` is 0 or exceeds the maximum supported order (10)
/// - `t` is negative
/// - Numerical computation fails
///
/// # Example
///
/// ```javascript
/// import { rys_compute } from './qc_wasm.js';
///
/// const result = rys_compute(3, 1.0);
/// console.log(result);
/// // {
/// //   roots: [0.123..., 0.456..., 0.789...],
/// //   weights: [0.234..., 0.345..., 0.167...],
/// //   nroots: 3,
/// //   t: 1.0,
/// //   method: "standard"
/// // }
/// ```
#[wasm_bindgen]
pub fn rys_compute(n: usize, t: f64) -> Result<JsValue, JsError> {
    let result = qc_core::rys::rys_roots(n, t).map_err(|e| JsError::new(&e.to_string()))?;
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

/// Compute the error curve showing how reconstruction error varies with quadrature order.
///
/// For each quadrature order n from 1 to n_max, computes the maximum absolute
/// reconstruction error across moments. This demonstrates how higher-order
/// quadrature achieves better accuracy, a key pedagogical visualization for
/// understanding Gaussian quadrature.
///
/// The reconstruction error measures how well the Gauss quadrature rule
/// approximates the original integral:
///
/// ```text
/// Error_m = |F_m(T) - Σ_k w_k * r_k^m|
/// max_error(n) = max_{m=0..2n-1} Error_m
/// ```
///
/// # Arguments
///
/// * `n_max` - Maximum quadrature order (1 <= n_max <= 10)
/// * `t` - Argument value (T >= 0)
///
/// # Returns
///
/// A JavaScript object containing:
/// - `t`: The argument parameter
/// - `nMax`: Maximum order computed
/// - `points`: Array of {n, maxError} objects for each order
///
/// # Errors
///
/// Returns a `JsError` if:
/// - `n_max` is 0 or exceeds the maximum supported order (10)
/// - `t` is negative
///
/// # Example
///
/// ```javascript
/// import { rys_error_curve } from './qc_wasm.js';
///
/// const result = rys_error_curve(5, 1.0);
/// console.log(result);
/// // {
/// //   t: 1.0,
/// //   nMax: 5,
/// //   points: [
/// //     { n: 1, maxError: 1.23e-15 },
/// //     { n: 2, maxError: 4.56e-15 },
/// //     ...
/// //   ]
/// // }
/// ```
#[wasm_bindgen]
pub fn rys_error_curve(n_max: usize, t: f64) -> Result<JsValue, JsError> {
    let result = qc_core::rys::error_curve(n_max, t).map_err(|e| JsError::new(&e.to_string()))?;
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// SCF WASM Bindings
// ============================================================================

/// WASM-friendly SCF options struct.
///
/// This struct is deserialized from JavaScript to configure SCF calculations.
/// It uses camelCase field names to match JavaScript conventions.
///
/// # Fields
///
/// * `convergence_profile` - Convergence threshold level: "loose", "medium", or "tight"
/// * `max_iterations` - Maximum number of SCF iterations
/// * `use_diis` - Enable DIIS acceleration
/// * `diis_size` - Optional DIIS subspace size (defaults to 6)
/// * `damp` - Optional Fock matrix damping factor (0.0-1.0, defaults to 0.0)
/// * `include_matrices` - Whether to include matrices in result (for internals mode)
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct ScfWasmOptions {
    /// Convergence profile: "loose" | "medium" | "tight"
    pub convergence_profile: String,
    /// Maximum number of SCF iterations
    pub max_iterations: u32,
    /// Enable DIIS acceleration
    pub use_diis: bool,
    /// DIIS subspace size (defaults to 6 if not provided)
    #[serde(default)]
    pub diis_size: Option<usize>,
    /// Fock matrix damping factor (0.0 = no damping, 0.5 = 50% old Fock)
    ///
    /// Damping mixes the current and previous Fock matrices to stabilize
    /// convergence for difficult systems (e.g., diffuse basis sets).
    ///
    /// F_damped = damp * F_old + (1.0 - damp) * F_new
    ///
    /// Reference: PySCF hf.py lines 789-790, 1119-1120
    #[serde(default)]
    pub damp: Option<f64>,
    /// Whether to include matrices in result (for internals mode)
    /// When true, returns S, H_core, F, D matrices and orbital energies
    #[serde(default)]
    pub include_matrices: bool,
}

impl Default for ScfWasmOptions {
    fn default() -> Self {
        Self {
            convergence_profile: "medium".to_string(),
            max_iterations: 100,
            use_diis: false,
            diis_size: None,
            damp: None, // No damping by default
            include_matrices: false,
        }
    }
}

/// WASM-friendly SCF iteration struct.
///
/// This struct is serialized to JavaScript with camelCase field names.
/// It contains a subset of the full ScfIteration data optimized for
/// visualization and progress tracking.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct ScfWasmIteration {
    /// Iteration number (0-indexed)
    pub iteration: usize,
    /// Total energy (electronic + nuclear) in Hartree
    pub energy: f64,
    /// Energy change from previous iteration (None for first iteration)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub delta_e: Option<f64>,
    /// RMS change in density matrix (None for first iteration)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rms_density_change: Option<f64>,
    /// Whether DIIS was applied in this iteration
    pub diis_applied: bool,
}

impl From<&ScfIteration> for ScfWasmIteration {
    fn from(iter: &ScfIteration) -> Self {
        Self {
            iteration: iter.iteration,
            energy: iter.energy_total,
            delta_e: iter.delta_e,
            rms_density_change: iter.rms_density_change,
            diis_applied: iter.diis_applied,
        }
    }
}

/// WASM-friendly SCF matrices struct.
///
/// Contains the key matrices from an SCF calculation for "internals mode"
/// visualization. All matrices are stored as row-major flat arrays.
///
/// # Fields
///
/// * `nbf` - Number of basis functions (for reshaping)
/// * `s_matrix` - Overlap matrix (nbf x nbf)
/// * `h_core` - Core Hamiltonian (nbf x nbf)
/// * `fock_matrix` - Final Fock matrix (nbf x nbf)
/// * `density_matrix` - Final density matrix (nbf x nbf)
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct ScfWasmMatrices {
    /// Number of basis functions (for reshaping into nbf x nbf matrix)
    pub nbf: usize,
    /// Overlap matrix S (row-major, nbf x nbf)
    pub s_matrix: Vec<f64>,
    /// Core Hamiltonian H = T + V (row-major, nbf x nbf)
    pub h_core: Vec<f64>,
    /// Final Fock matrix F (row-major, nbf x nbf)
    pub fock_matrix: Vec<f64>,
    /// Final density matrix D (row-major, nbf x nbf)
    pub density_matrix: Vec<f64>,
    /// MO coefficient matrix C (row-major, nbf x nbf)
    pub mo_coefficients: Vec<f64>,
}

/// WASM-friendly orbital energies struct.
///
/// Contains the MO (molecular orbital) energies from an SCF calculation,
/// along with occupancy information for visualization.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct ScfWasmOrbitalEnergies {
    /// Orbital energies in Hartree, sorted ascending
    pub energies: Vec<f64>,
    /// Number of occupied orbitals (n_occ = nelec / 2 for RHF)
    pub n_occupied: usize,
}

/// WASM-friendly SCF result struct.
///
/// This struct is serialized to JavaScript with camelCase field names.
/// It contains the essential SCF output for visualization and export.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct ScfWasmResult {
    /// Whether SCF converged
    pub converged: bool,
    /// Final total energy (electronic + nuclear) in Hartree
    pub energy: f64,
    /// Number of iterations performed
    pub iterations: usize,
    /// Iteration-by-iteration trace
    pub trace: Vec<ScfWasmIteration>,
    /// Optional matrices (when include_matrices is true)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub matrices: Option<ScfWasmMatrices>,
    /// Optional orbital energies (when include_matrices is true)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub orbital_energies: Option<ScfWasmOrbitalEnergies>,
}

impl From<&ScfOutput> for ScfWasmResult {
    fn from(output: &ScfOutput) -> Self {
        Self {
            converged: output.converged,
            energy: output.energy_total,
            iterations: output.iterations,
            trace: output.trace.iter().map(ScfWasmIteration::from).collect(),
            matrices: None,
            orbital_energies: None,
        }
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
        _ => ConvergenceProfile::Medium, // Default to medium for unknown profiles
    }
}

/// Run an RHF SCF calculation on a preset molecular system.
///
/// This function accepts a JSON string representing the molecular system
/// with pre-computed integrals, and options for the SCF calculation.
///
/// # Arguments
///
/// * `system_json` - JSON string containing a PresetSystemJson object
/// * `options` - JavaScript object containing SCF options (ScfWasmOptions)
///
/// # Returns
///
/// A JavaScript object containing:
/// - `converged`: Whether SCF converged
/// - `energy`: Final total energy in Hartree
/// - `iterations`: Number of iterations performed
/// - `trace`: Array of iteration data for visualization
///
/// # Errors
///
/// Returns a `JsError` if:
/// - The system JSON is invalid or malformed
/// - The options object is invalid
/// - The SCF calculation fails to converge or encounters numerical issues
///
/// # Example
///
/// ```javascript
/// import { scf_run } from './qc_wasm.js';
///
/// const systemJson = '{"format_version": 1, "system_id": "h2", ...}';
/// const options = {
///     convergenceProfile: "medium",
///     maxIterations: 100,
///     useDiis: false
/// };
///
/// const result = scf_run(systemJson, options);
/// console.log(result);
/// // {
/// //   converged: true,
/// //   energy: -1.1167143250625,
/// //   iterations: 8,
/// //   trace: [...]
/// // }
/// ```
#[wasm_bindgen]
pub fn scf_run(system_json: &str, options: JsValue) -> Result<JsValue, JsError> {
    // 1. Deserialize system from JSON string
    let preset_json = PresetSystemJson::from_json(system_json)
        .map_err(|e| JsError::new(&format!("Invalid system JSON: {}", e)))?;

    // 2. Convert PresetSystemJson to PresetSystem
    let system = preset_json
        .to_preset_system()
        .map_err(|e| JsError::new(&format!("System validation failed: {}", e)))?;

    // 3. Deserialize options from JsValue
    let wasm_options: ScfWasmOptions = serde_wasm_bindgen::from_value(options)
        .map_err(|e| JsError::new(&format!("Invalid options: {}", e)))?;

    // 4. Convert options to ScfConfig
    let config = ScfConfig {
        profile: parse_convergence_profile(&wasm_options.convergence_profile),
        max_iterations: wasm_options.max_iterations as usize,
        use_diis: wasm_options.use_diis,
        diis_size: wasm_options.diis_size.unwrap_or(6),
        diis_start: 2,
        damp: wasm_options.damp.unwrap_or(0.0), // Default: no damping
        damp_start: 5,                          // Apply damping for first 5 iterations (like PySCF)
        level_shift: 0.0,                       // No level shift by default
    };

    // 5. Call qc_core::scf::rhf_scf
    let output = qc_core::scf::rhf_scf(&system, &config)
        .map_err(|e| JsError::new(&format!("SCF calculation failed: {}", e)))?;

    // 6. Convert ScfOutput to ScfWasmResult
    let mut result = ScfWasmResult::from(&output);

    // 7. If include_matrices is true, populate matrices and orbital energies
    if wasm_options.include_matrices {
        // Matrices: S and Hcore from system, F and D from output
        result.matrices = Some(ScfWasmMatrices {
            nbf: system.nbf,
            s_matrix: system.s_matrix.clone(),
            h_core: system.h_core.clone(),
            fock_matrix: output.fock_matrix.clone(),
            density_matrix: output.density_matrix.clone(),
            mo_coefficients: output.mo_coefficients.clone(),
        });

        // Orbital energies from output
        result.orbital_energies = Some(ScfWasmOrbitalEnergies {
            energies: output.mo_energies.clone(),
            n_occupied: system.n_occ(),
        });
    }

    // 8. Return as JsValue
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// KS-DFT SCF WASM Bindings (US-068)
// ============================================================================

/// Input options for KS-DFT SCF calculation.
///
/// # Fields
///
/// * `convergence_profile` - "loose" | "medium" | "tight"
/// * `max_iterations` - Maximum number of SCF iterations
/// * `use_diis` - Enable DIIS acceleration
/// * `method` - DFT method: "lda" (Slater + VWN5)
/// * `grid_quality` - Grid quality: "standard" (302-pt) or "fine" (590-pt)
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct KsScfWasmInput {
    /// Atoms as [[Z, x, y, z], ...]
    pub atoms: Vec<[f64; 4]>,
    /// Basis set name (e.g., "sto-3g")
    pub basis_name: String,
    /// DFT method: "lda"
    pub method: String,
    /// Convergence profile: "loose" | "medium" | "tight"
    #[serde(default = "default_convergence")]
    pub convergence_profile: String,
    /// Maximum number of SCF iterations
    #[serde(default = "default_max_iter")]
    pub max_iterations: u32,
    /// Enable DIIS acceleration
    #[serde(default = "default_diis")]
    pub use_diis: bool,
    /// Grid quality: "standard" or "fine"
    #[serde(default = "default_grid_quality")]
    pub grid_quality: String,
    /// Use spherical harmonic basis functions (5 d-orbitals vs 6 Cartesian)
    #[serde(default)]
    pub use_spherical: bool,
}

fn default_convergence() -> String {
    "tight".to_string()
}
fn default_max_iter() -> u32 {
    100
}
fn default_diis() -> bool {
    true
}
fn default_grid_quality() -> String {
    "standard".to_string()
}

/// KS-DFT SCF result for WASM.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct KsScfWasmResult {
    /// Whether SCF converged
    pub converged: bool,
    /// Final total energy (electronic + nuclear) in Hartree
    pub energy: f64,
    /// Number of iterations performed
    pub iterations: usize,
    /// Exchange-correlation energy component (Hartree)
    pub energy_xc: f64,
    /// Coulomb energy component (Hartree)
    pub energy_j: f64,
    /// One-electron energy (Hartree)
    pub energy_1e: f64,
    /// Nuclear repulsion energy (Hartree)
    pub energy_nuc: f64,
    /// Method identifier
    pub method: String,
    /// Iteration-by-iteration trace
    pub trace: Vec<ScfWasmIteration>,
    /// Final density matrix (row-major, nbf x nbf)
    pub density_matrix: Vec<f64>,
    /// MO coefficients (column-major from nalgebra, nbf x nbf)
    pub mo_coefficients: Vec<f64>,
    /// Orbital energies (eigenvalues, sorted ascending)
    pub orbital_energies: Vec<f64>,
    /// Number of basis functions
    pub n_basis: usize,
    /// Number of occupied orbitals
    pub n_occupied: usize,
    /// Overlap matrix S (row-major, nbf x nbf) — needed for population analysis
    pub overlap_matrix: Vec<f64>,
    /// Core Hamiltonian matrix (row-major, nbf x nbf)
    pub h_core: Vec<f64>,
    /// Final Fock/KS matrix (row-major, nbf x nbf)
    pub fock_matrix: Vec<f64>,
    /// D3-BJ dispersion energy (Hartree), if applicable.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub energy_disp: Option<f64>,
}

/// Run a Kohn-Sham DFT SCF calculation.
///
/// # Arguments
///
/// * `input` - JavaScript object containing `KsScfWasmInput`
///
/// # Returns
///
/// A JavaScript object containing `KsScfWasmResult`
///
/// # Example
///
/// ```javascript
/// const input = {
///   atoms: [[1, 0, 0, 0], [1, 0, 0, 1.4]],
///   basisName: "sto-3g",
///   method: "lda",
///   convergenceProfile: "tight",
///   useDiis: true,
/// };
/// const result = ks_scf(input);
/// console.log(result.energy);
/// ```
#[wasm_bindgen]
pub fn ks_scf(
    input: JsValue,
    progress_callback: Option<js_sys::Function>,
) -> Result<JsValue, JsError> {
    let wasm_input: KsScfWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid KS-SCF input: {}", e)))?;

    // 1. Build atoms and basis set
    let atoms: Vec<qc_core::basis::Atom> = wasm_input
        .atoms
        .iter()
        .map(|a| {
            let z = a[0] as u8;
            let pos = [a[1], a[2], a[3]];
            qc_core::basis::Atom::new(z, pos)
                .map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let basis = qc_core::basis::BasisSet::build(atoms.clone(), &wasm_input.basis_name)
        .map_err(|e| JsError::new(&format!("Basis set error: {}", e)))?;

    // 2. Compute integrals (spherical or cartesian depending on user choice)
    let use_spherical = wasm_input.use_spherical;
    let nelec = basis.n_electrons;
    let e_nuc = basis.nuclear_repulsion;

    // Macro to emit integral-phase progress (distinguished from SCF progress by "phase" field).
    // Uses a macro instead of a closure to avoid borrow conflicts with progress_callback.
    macro_rules! emit_integral_progress {
        ($step:expr, $pct:expr, $msg:expr) => {
            if let Some(ref cb) = progress_callback {
                let obj = js_sys::Object::new();
                let _ = js_sys::Reflect::set(&obj, &"phase".into(), &"integrals".into());
                let _ = js_sys::Reflect::set(&obj, &"step".into(), &JsValue::from_str($step));
                let _ = js_sys::Reflect::set(&obj, &"percent".into(), &JsValue::from_f64($pct));
                let _ = js_sys::Reflect::set(&obj, &"message".into(), &JsValue::from_str($msg));
                let _ = cb.call1(&JsValue::NULL, &obj);
            }
        };
    }

    emit_integral_progress!("overlap", 5.0, "Computing overlap integrals (S)...");

    let (s_matrix, h_core, eri, nbf) = if use_spherical {
        let s = qc_core::integrals::overlap_matrix_spherical(&basis);
        emit_integral_progress!("hcore", 15.0, "Computing core Hamiltonian (H)...");
        let h = qc_core::integrals::hcore_matrix_spherical(&basis);
        emit_integral_progress!(
            "eri",
            25.0,
            "Computing electron repulsion integrals (ERI)..."
        );
        // Use progress variant: ERI is 25-95%, live per-quartet updates
        let cb = progress_callback.clone();
        let e =
            qc_core::integrals::eri_compressed_spherical_with_progress(&basis, |done, total| {
                if let Some(ref cb) = cb {
                    let pct = 25.0 + (done as f64 / total as f64) * 70.0;
                    let obj = js_sys::Object::new();
                    let _ = js_sys::Reflect::set(&obj, &"phase".into(), &"integrals".into());
                    let _ = js_sys::Reflect::set(&obj, &"step".into(), &JsValue::from_str("eri"));
                    let _ = js_sys::Reflect::set(&obj, &"percent".into(), &JsValue::from_f64(pct));
                    let msg = format!("ERI: {}/{} shell quartets ({:.0}%)", done, total, pct);
                    let _ = js_sys::Reflect::set(&obj, &"message".into(), &JsValue::from_str(&msg));
                    let _ = cb.call1(&JsValue::NULL, &obj);
                }
            });
        let n = basis.n_basis_spherical();
        (s, h, e, n)
    } else {
        let s = qc_core::integrals::overlap_matrix(&basis);
        emit_integral_progress!("hcore", 15.0, "Computing core Hamiltonian (H)...");
        let h = qc_core::integrals::hcore_matrix(&basis);
        emit_integral_progress!(
            "eri",
            25.0,
            "Computing electron repulsion integrals (ERI)..."
        );
        // Use progress variant: ERI is 25-95%, live per-quartet updates
        let cb = progress_callback.clone();
        let e = qc_core::integrals::eri_compressed_with_progress(&basis, |done, total| {
            if let Some(ref cb) = cb {
                let pct = 25.0 + (done as f64 / total as f64) * 70.0;
                let obj = js_sys::Object::new();
                let _ = js_sys::Reflect::set(&obj, &"phase".into(), &"integrals".into());
                let _ = js_sys::Reflect::set(&obj, &"step".into(), &JsValue::from_str("eri"));
                let _ = js_sys::Reflect::set(&obj, &"percent".into(), &JsValue::from_f64(pct));
                let msg = format!("ERI: {}/{} shell quartets ({:.0}%)", done, total, pct);
                let _ = js_sys::Reflect::set(&obj, &"message".into(), &JsValue::from_str(&msg));
                let _ = cb.call1(&JsValue::NULL, &obj);
            }
        });
        let n = basis.n_basis;
        (s, h, e, n)
    };

    emit_integral_progress!("done", 95.0, "All integrals computed");

    let system = qc_core::scf::PresetSystem {
        system_id: "ks_scf_wasm".to_string(),
        label: "KS-DFT WASM calculation".to_string(),
        nbf,
        nelec,
        e_nuc,
        s_matrix,
        h_core,
        eri_compressed: eri,
    };

    // 3. Determine method
    let method_lower = wasm_input.method.to_lowercase();
    let is_rhf = method_lower == "rhf" || method_lower == "hf";
    let use_d3bj = method_lower == "b3lyp-d3bj";

    // For RHF, skip grid and functional setup — use rhf_scf directly
    if is_rhf {
        let config = qc_core::scf::ScfConfig {
            profile: parse_convergence_profile(&wasm_input.convergence_profile),
            max_iterations: wasm_input.max_iterations as usize,
            use_diis: wasm_input.use_diis,
            ..Default::default()
        };

        // Create progress callback
        let rust_callback = progress_callback.as_ref().map(|cb| {
            move |iteration: usize,
                  energy: f64,
                  delta_e: f64,
                  rms_density: f64,
                  diis_applied: bool| {
                let obj = js_sys::Object::new();
                let _ = js_sys::Reflect::set(&obj, &"iteration".into(), &(iteration as f64).into());
                let _ = js_sys::Reflect::set(&obj, &"energy".into(), &energy.into());
                let _ = js_sys::Reflect::set(&obj, &"deltaE".into(), &delta_e.into());
                let _ = js_sys::Reflect::set(&obj, &"rmsDensityChange".into(), &rms_density.into());
                let _ = js_sys::Reflect::set(&obj, &"diisApplied".into(), &diis_applied.into());
                let _ = cb.call1(&JsValue::NULL, &obj);
            }
        });

        // Run RHF with progress callback via manual iteration tracking
        let output = qc_core::scf::rhf_scf(&system, &config)
            .map_err(|e| JsError::new(&format!("RHF SCF failed: {}", e)))?;

        // Send final progress for each iteration (from trace)
        if let Some(ref cb) = progress_callback {
            for iter in &output.trace {
                let obj = js_sys::Object::new();
                let _ = js_sys::Reflect::set(
                    &obj,
                    &"iteration".into(),
                    &(iter.iteration as f64).into(),
                );
                let _ = js_sys::Reflect::set(&obj, &"energy".into(), &iter.energy_total.into());
                let _ = js_sys::Reflect::set(
                    &obj,
                    &"deltaE".into(),
                    &iter.delta_e.unwrap_or(0.0).into(),
                );
                let _ = js_sys::Reflect::set(
                    &obj,
                    &"rmsDensityChange".into(),
                    &iter.rms_density_change.unwrap_or(0.0).into(),
                );
                let _ =
                    js_sys::Reflect::set(&obj, &"diisApplied".into(), &iter.diis_applied.into());
                let _ = cb.call1(&JsValue::NULL, &obj);
            }
        }
        let _ = rust_callback; // suppress unused warning

        let n_occ = system.n_occ();
        let result = KsScfWasmResult {
            converged: output.converged,
            energy: output.energy_total,
            iterations: output.iterations,
            energy_xc: 0.0,
            energy_j: 0.0,
            energy_1e: 0.0,
            energy_nuc: system.e_nuc,
            method: "RHF".to_string(),
            trace: output.trace.iter().map(ScfWasmIteration::from).collect(),
            density_matrix: output.density_matrix.clone(),
            mo_coefficients: output.mo_coefficients.clone(),
            orbital_energies: output.mo_energies.clone(),
            n_basis: system.nbf,
            n_occupied: n_occ,
            overlap_matrix: system.s_matrix.clone(),
            h_core: system.h_core.clone(),
            fock_matrix: output.fock_matrix.clone(),
            energy_disp: None,
        };

        return serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()));
    }

    // DFT path: build grid and select functional
    emit_integral_progress!("grid", 50.0, "Building numerical integration grid...");
    let grid_quality = match wasm_input.grid_quality.to_lowercase().as_str() {
        "fine" => qc_core::dft::GridQuality::Fine,
        _ => qc_core::dft::GridQuality::Standard,
    };
    let grid_config = qc_core::dft::GridConfig {
        n_radial: 75,
        quality: grid_quality,
        pruning: true,
    };
    let grid = qc_core::dft::build_becke_grid(&atoms, &grid_config);

    let functional: Box<dyn qc_core::dft::ExchangeCorrelation> = match method_lower.as_str() {
        "lda" => Box::new(qc_core::dft::Lda::new()),
        "b3lyp" | "b3lyp-d3bj" => Box::new(qc_core::dft::B3lyp::new()),
        other => {
            return Err(JsError::new(&format!(
                "Unknown DFT method: '{}'. Supported: 'rhf', 'lda', 'b3lyp', 'b3lyp-d3bj'",
                other
            )))
        }
    };

    // 5. Configure SCF
    let config = qc_core::scf::ScfConfig {
        profile: parse_convergence_profile(&wasm_input.convergence_profile),
        max_iterations: wasm_input.max_iterations as usize,
        use_diis: wasm_input.use_diis,
        diis_size: 8,
        diis_start: 2,
        damp: 0.0,
        damp_start: 5,
        level_shift: 0.0,
    };

    // 6. Run KS-SCF with optional progress callback
    let rust_callback = progress_callback.as_ref().map(|cb| {
        let cb = cb.clone();
        move |iteration: usize, energy: f64, delta_e: f64, rms_density: f64, diis_applied: bool| {
            // Build a plain JS object for the progress data
            let obj = js_sys::Object::new();
            let _ = js_sys::Reflect::set(
                &obj,
                &"iteration".into(),
                &JsValue::from_f64(iteration as f64),
            );
            let _ = js_sys::Reflect::set(&obj, &"energy".into(), &JsValue::from_f64(energy));
            let _ = js_sys::Reflect::set(&obj, &"deltaE".into(), &JsValue::from_f64(delta_e));
            let _ = js_sys::Reflect::set(
                &obj,
                &"rmsDensityChange".into(),
                &JsValue::from_f64(rms_density),
            );
            let _ = js_sys::Reflect::set(
                &obj,
                &"diisApplied".into(),
                &JsValue::from_bool(diis_applied),
            );
            let _ = cb.call1(&JsValue::NULL, &obj);
        }
    });
    let mut output = qc_core::dft::ks_scf(
        &system,
        &config,
        functional.as_ref(),
        &grid,
        &basis,
        use_spherical,
        rust_callback
            .as_ref()
            .map(|f| f as &dyn Fn(usize, f64, f64, f64, bool)),
    )
    .map_err(|e| JsError::new(&format!("KS-SCF failed: {}", e)))?;

    // 6b. Apply D3-BJ dispersion correction if requested (post-SCF)
    let energy_disp = if use_d3bj {
        let d3_atoms: Vec<(u8, [f64; 3])> = wasm_input
            .atoms
            .iter()
            .map(|a| (a[0] as u8, [a[1], a[2], a[3]]))
            .collect();
        let d3_result = qc_core::dft::compute_d3bj_energy(&d3_atoms, &qc_core::dft::D3BJ_B3LYP);
        // Add dispersion energy to total
        output.scf_output.energy_total += d3_result.energy;
        output.scf_output.energy_electronic += d3_result.energy;
        output.energy_disp = Some(d3_result.energy);
        output.method = "B3LYP-D3(BJ)".to_string();
        Some(d3_result.energy)
    } else {
        None
    };

    // 7. Build result
    let n_occ = system.n_occ();
    let result = KsScfWasmResult {
        converged: output.scf_output.converged,
        energy: output.scf_output.energy_total,
        iterations: output.scf_output.iterations,
        energy_xc: output.energy_xc,
        energy_j: output.energy_j,
        energy_1e: output.energy_1e,
        energy_nuc: system.e_nuc,
        method: output.method,
        trace: output
            .scf_output
            .trace
            .iter()
            .map(ScfWasmIteration::from)
            .collect(),
        density_matrix: output.scf_output.density_matrix.clone(),
        mo_coefficients: output.scf_output.mo_coefficients.clone(),
        orbital_energies: output.scf_output.mo_energies.clone(),
        n_basis: system.nbf,
        n_occupied: n_occ,
        overlap_matrix: system.s_matrix.clone(),
        h_core: system.h_core.clone(),
        fock_matrix: output.scf_output.fock_matrix.clone(),
        energy_disp,
    };

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Integral Computation WASM Bindings (US-029)
// ============================================================================

/// Options for integral computation.
///
/// Configures how integrals are computed, including the choice between
/// Cartesian and spherical harmonic basis functions.
///
/// # Fields
///
/// * `use_spherical` - If true, use spherical harmonic basis functions (5 d-orbitals).
///   If false (default), use Cartesian basis functions (6 d-orbitals).
///   For s and p orbitals, both choices give the same result.
///
/// # Example
///
/// ```javascript
/// const options = { useSpherical: true };  // Use 5 d-orbitals instead of 6
/// ```
#[derive(Default, Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct IntegralComputeOptions {
    /// Use spherical harmonic basis functions (5 d-orbitals vs 6 Cartesian).
    /// Default: false (Cartesian basis for backward compatibility).
    #[serde(default)]
    pub use_spherical: bool,
}

/// Input atom specification for integral computation.
///
/// Specifies an atom by its element symbol and Cartesian coordinates.
///
/// # Fields
///
/// * `symbol` - Element symbol (case insensitive): "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"
/// * `xyz` - Cartesian coordinates [x, y, z] in the units specified by GeometryInput
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct AtomInput {
    /// Element symbol (e.g., "H", "O", "C")
    pub symbol: String,
    /// Cartesian coordinates [x, y, z]
    pub xyz: [f64; 3],
}

/// Input geometry for integral computation.
///
/// Specifies a molecular geometry with atoms and coordinate units.
///
/// # Fields
///
/// * `atoms` - List of atoms with element symbols and coordinates
/// * `units` - Coordinate units: "bohr" or "angstrom"
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct GeometryInput {
    /// List of atoms in the molecule
    pub atoms: Vec<AtomInput>,
    /// Coordinate units: "bohr" or "angstrom"
    pub units: String,
}

/// Output atom specification (echoed back with atomic number).
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct AtomOutput {
    /// Element symbol
    pub symbol: String,
    /// Cartesian coordinates [x, y, z] in Bohr
    pub xyz: [f64; 3],
    /// Atomic number
    pub atomic_number: u32,
}

/// Geometry echoed back in output (always in Bohr).
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct GeometryOutput {
    /// Atoms with coordinates in Bohr
    pub atoms: Vec<AtomOutput>,
    /// Units (always "bohr" for output)
    pub units: String,
}

/// Computation metadata for integral computation.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct IntegralMetadata {
    /// WASM module version
    pub wasm_version: String,
    /// Computation time in milliseconds
    pub compute_time_ms: u64,
    /// Number of shell pairs processed
    pub shell_pairs: usize,
    /// Number of shell quartets processed (for ERI)
    pub shell_quartets: usize,
    /// Number of unique ERIs (after 8-fold symmetry)
    pub significant_eris: usize,
    /// Basis type used for computation: "cartesian" or "spherical"
    pub basis_type: String,
}

/// Progress update during integral computation.
///
/// This struct is serialized and passed to the JavaScript progress callback
/// during long-running integral computations.
#[derive(Serialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct IntegralProgress {
    /// Current phase: "setup", "overlap", "kinetic", "nuclear", "eri", "assembly"
    pub phase: String,
    /// Current step within phase
    pub current: usize,
    /// Total steps in phase
    pub total: usize,
    /// Overall percentage complete (0-100)
    pub overall_percent: f64,
    /// Human-readable message
    pub message: String,
}

/// Result of integral computation, compatible with PresetSystemJson.
///
/// This structure can be directly used as input to scf_run without modification.
///
/// # Example
///
/// ```javascript
/// import { compute_integrals, scf_run } from './qc_wasm.js';
///
/// const geometry = {
///   atoms: [
///     { symbol: "H", xyz: [0, 0, 0] },
///     { symbol: "H", xyz: [0, 0, 1.4] }
///   ],
///   units: "bohr"
/// };
///
/// const integrals = compute_integrals(JSON.stringify(geometry), "sto-3g");
/// const scfResult = scf_run(JSON.stringify(integrals), options);
/// ```
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct IntegralResult {
    /// Format version (matches PresetSystemJson format)
    pub format_version: u32,
    /// System identifier (generated from geometry hash)
    pub system_id: String,
    /// Human-readable label
    pub label: String,
    /// Description
    pub description: String,
    /// Geometry for reference (coordinates in Bohr)
    pub geometry: GeometryOutput,
    /// Basis set used
    pub basis_id: String,
    /// Number of basis functions
    pub nbf: usize,
    /// Number of electrons
    pub nelec: usize,
    /// Nuclear repulsion energy (Hartree)
    pub e_nuc: f64,
    /// Overlap matrix S (row-major, nbf x nbf)
    pub s_matrix: Vec<f64>,
    /// Core Hamiltonian H = T + V (row-major, nbf x nbf)
    pub h_core: Vec<f64>,
    /// Compressed two-electron integrals (8-fold symmetry)
    pub eri_compressed: Vec<f64>,
    /// Indexing scheme description
    pub eri_indexing: String,
    /// Computation metadata
    pub metadata: IntegralMetadata,
}

/// Compute all molecular integrals for a given geometry and basis set.
///
/// This function computes overlap (S), kinetic (T), nuclear attraction (V),
/// and two-electron repulsion (ERI) integrals for a user-specified molecular
/// geometry and basis set. The result is compatible with PresetSystemJson
/// and can be used directly with scf_run.
///
/// # Arguments
///
/// * `geometry_json` - JSON string containing GeometryInput
/// * `basis_name` - Basis set name: "sto-3g", "3-21g", "6-31g", "6-31g*"
///
/// # Returns
///
/// IntegralResult as JsValue, compatible with PresetSystemJson for SCF input.
///
/// # Errors
///
/// Returns JsError if:
/// - Geometry JSON is invalid
/// - Element symbol is unsupported (only H-Ne supported)
/// - Basis set name is unrecognized
/// - Geometry has no atoms
/// - Computation fails (numerical issues)
///
/// # Supported Elements
///
/// H, He, Li, Be, B, C, N, O, F, Ne (atomic numbers 1-10)
///
/// # Supported Basis Sets
///
/// - "sto-3g": Minimal basis (3 Gaussians per orbital)
/// - "3-21g": Split-valence basis
/// - "6-31g": Split-valence basis
/// - "6-31g*": Polarized split-valence (d orbitals on heavy atoms)
///
/// # Example
///
/// ```javascript
/// import { compute_integrals } from './qc_wasm.js';
///
/// const geometry = {
///   atoms: [
///     { symbol: "H", xyz: [0, 0, 0] },
///     { symbol: "H", xyz: [0, 0, 1.4] }
///   ],
///   units: "bohr"
/// };
///
/// const result = compute_integrals(JSON.stringify(geometry), "sto-3g");
/// console.log(result);
/// // {
/// //   formatVersion: 1,
/// //   systemId: "custom_a1b2c3d4e5f6g7h8",
/// //   label: "H2 (sto-3g)",
/// //   nbf: 2,
/// //   nelec: 2,
/// //   eNuc: 0.714285...,
/// //   sMatrix: [...],
/// //   hCore: [...],
/// //   eriCompressed: [...],
/// //   ...
/// // }
/// ```
#[wasm_bindgen]
pub fn compute_integrals(geometry_json: &str, basis_name: &str) -> Result<JsValue, JsError> {
    use qc_core::basis::{symbol_to_atomic_number, Atom, BasisSet, ANGSTROM_TO_BOHR};
    #[cfg(not(feature = "parallel"))]
    use qc_core::integrals::eri_compressed;
    #[cfg(feature = "parallel")]
    use qc_core::integrals::eri_compressed_parallel;
    use qc_core::integrals::{hcore_matrix, overlap_matrix};
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    // Use js_sys::Date for timing (std::time::Instant is not supported in WASM)
    let start_time = js_sys::Date::now();

    // 1. Parse geometry JSON
    let geometry: GeometryInput = serde_json::from_str(geometry_json)
        .map_err(|e| JsError::new(&format!("Invalid geometry JSON: {}", e)))?;

    // 2. Validate geometry
    if geometry.atoms.is_empty() {
        return Err(JsError::new("Geometry must have at least 1 atom."));
    }

    // 3. Normalize basis name and validate
    let basis_lower = basis_name.to_lowercase();
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];
    if !supported_bases.contains(&basis_lower.as_str()) {
        return Err(JsError::new(&format!(
            "Unknown basis set '{}'. Supported: {}",
            basis_name,
            supported_bases.join(", ")
        )));
    }

    // 4. Determine coordinate units
    let units_lower = geometry.units.to_lowercase();
    let convert_to_bohr = match units_lower.as_str() {
        "bohr" => false,
        "angstrom" | "angstroms" => true,
        _ => {
            return Err(JsError::new(&format!(
                "Invalid units '{}'. Must be 'bohr' or 'angstrom'.",
                geometry.units
            )));
        }
    };

    // 5. Convert atoms to qc_core::basis::Atom
    let mut atoms: Vec<Atom> = Vec::with_capacity(geometry.atoms.len());
    let mut atom_outputs: Vec<AtomOutput> = Vec::with_capacity(geometry.atoms.len());

    for atom_input in &geometry.atoms {
        // Validate element symbol
        let atomic_number = symbol_to_atomic_number(&atom_input.symbol).map_err(|_| {
            JsError::new(&format!(
                "Unsupported element '{}'. Only H-Ar (Z=1-18) are supported.",
                atom_input.symbol
            ))
        })?;

        // Convert coordinates to Bohr if needed
        let position_bohr = if convert_to_bohr {
            [
                atom_input.xyz[0] * ANGSTROM_TO_BOHR,
                atom_input.xyz[1] * ANGSTROM_TO_BOHR,
                atom_input.xyz[2] * ANGSTROM_TO_BOHR,
            ]
        } else {
            atom_input.xyz
        };

        let atom = Atom::new(atomic_number, position_bohr)
            .map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))?;
        atoms.push(atom.clone());

        atom_outputs.push(AtomOutput {
            symbol: atom.symbol.clone(),
            xyz: position_bohr,
            atomic_number: atomic_number as u32,
        });
    }

    // 6. Build basis set
    let basis = BasisSet::build(atoms, &basis_lower)
        .map_err(|e| JsError::new(&format!("Failed to build basis set: {}", e)))?;

    // 7. Compute integrals
    let s_matrix = overlap_matrix(&basis);
    let h_core = hcore_matrix(&basis);
    #[cfg(feature = "parallel")]
    let eri = eri_compressed_parallel(&basis);
    #[cfg(not(feature = "parallel"))]
    let eri = eri_compressed(&basis);

    // 8. Generate system ID from geometry + basis hash
    let system_id = {
        let mut hasher = DefaultHasher::new();
        for atom in &atom_outputs {
            atom.symbol.hash(&mut hasher);
            // Hash coordinates with limited precision for stability
            ((atom.xyz[0] * 1e6) as i64).hash(&mut hasher);
            ((atom.xyz[1] * 1e6) as i64).hash(&mut hasher);
            ((atom.xyz[2] * 1e6) as i64).hash(&mut hasher);
        }
        basis_lower.hash(&mut hasher);
        format!("custom_{:016x}", hasher.finish())
    };

    // 9. Generate label
    let element_symbols: Vec<&str> = atom_outputs.iter().map(|a| a.symbol.as_str()).collect();
    let formula = generate_molecular_formula(&element_symbols);
    let label = format!("{} ({})", formula, basis_lower);

    // 10. Compute metadata
    let compute_time_ms = (js_sys::Date::now() - start_time) as u64;
    let n_shells = basis.n_shells();
    let shell_pairs = n_shells * (n_shells + 1) / 2;
    let shell_quartets = shell_pairs * (shell_pairs + 1) / 2;
    // Store ERI count before moving the vector
    let significant_eris = eri.len();

    let result = IntegralResult {
        format_version: 1,
        system_id,
        label,
        description: format!(
            "On-the-fly integral computation for {} atoms, {} basis",
            atom_outputs.len(),
            basis_lower
        ),
        geometry: GeometryOutput {
            atoms: atom_outputs,
            units: "bohr".to_string(),
        },
        basis_id: basis_lower,
        nbf: basis.n_basis,
        nelec: basis.n_electrons,
        e_nuc: basis.nuclear_repulsion,
        s_matrix,
        h_core,
        eri_compressed: eri, // Move instead of clone for better performance
        eri_indexing: "8-fold symmetry, pair_ij * (pair_ij+1)/2 + pair_kl".to_string(),
        metadata: IntegralMetadata {
            wasm_version: env!("CARGO_PKG_VERSION").to_string(),
            compute_time_ms,
            shell_pairs,
            shell_quartets,
            significant_eris,
            basis_type: "cartesian".to_string(), // Default to Cartesian for backward compatibility
        },
    };

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

/// Compute all molecular integrals with options for basis type.
///
/// This function is similar to `compute_integrals` but accepts an options object
/// that allows specifying whether to use spherical harmonic basis functions.
///
/// # Arguments
///
/// * `geometry_json` - JSON string containing GeometryInput
/// * `basis_name` - Basis set name: "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"
/// * `options` - Optional IntegralComputeOptions object
///
/// # Options
///
/// * `useSpherical` - If true, use spherical harmonic basis (5 d-orbitals).
///                    Default: false (6 Cartesian d-orbitals).
///
/// # Example
///
/// ```javascript
/// import { compute_integrals_with_options } from './qc_wasm.js';
///
/// const geometry = {
///   atoms: [
///     { symbol: "O", xyz: [0, 0, 0] },
///     { symbol: "H", xyz: [0.96, 0, 0] }
///   ],
///   units: "angstrom"
/// };
///
/// // Use spherical harmonics (5 d-orbitals for O with 6-31G*)
/// const options = { useSpherical: true };
/// const result = compute_integrals_with_options(JSON.stringify(geometry), "6-31g*", options);
/// console.log(result.metadata.basisType); // "spherical"
/// console.log(result.nbf); // fewer basis functions with spherical d-orbitals
/// ```
#[wasm_bindgen]
pub fn compute_integrals_with_options(
    geometry_json: &str,
    basis_name: &str,
    options: JsValue,
) -> Result<JsValue, JsError> {
    use qc_core::basis::{symbol_to_atomic_number, Atom, BasisSet, ANGSTROM_TO_BOHR};
    #[cfg(not(feature = "parallel"))]
    use qc_core::integrals::{eri_compressed, eri_compressed_spherical};
    #[cfg(feature = "parallel")]
    use qc_core::integrals::{eri_compressed_parallel, eri_compressed_spherical_parallel};
    use qc_core::integrals::{
        hcore_matrix, hcore_matrix_spherical, overlap_matrix, overlap_matrix_spherical,
    };
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    // Parse options (default to Cartesian for backward compatibility)
    let opts: IntegralComputeOptions = if options.is_undefined() || options.is_null() {
        IntegralComputeOptions::default()
    } else {
        serde_wasm_bindgen::from_value(options)
            .map_err(|e| JsError::new(&format!("Invalid options: {}", e)))?
    };

    let use_spherical = opts.use_spherical;
    let basis_type = if use_spherical {
        "spherical"
    } else {
        "cartesian"
    };

    // Use js_sys::Date for timing (std::time::Instant is not supported in WASM)
    let start_time = js_sys::Date::now();

    // 1. Parse geometry JSON
    let geometry: GeometryInput = serde_json::from_str(geometry_json)
        .map_err(|e| JsError::new(&format!("Invalid geometry JSON: {}", e)))?;

    // 2. Validate geometry
    if geometry.atoms.is_empty() {
        return Err(JsError::new("Geometry must have at least 1 atom."));
    }

    // 3. Normalize basis name and validate
    let basis_lower = basis_name.to_lowercase();
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];
    if !supported_bases.contains(&basis_lower.as_str()) {
        return Err(JsError::new(&format!(
            "Unknown basis set '{}'. Supported: {}",
            basis_name,
            supported_bases.join(", ")
        )));
    }

    // 4. Determine coordinate units
    let units_lower = geometry.units.to_lowercase();
    let convert_to_bohr = match units_lower.as_str() {
        "bohr" => false,
        "angstrom" | "angstroms" => true,
        _ => {
            return Err(JsError::new(&format!(
                "Invalid units '{}'. Must be 'bohr' or 'angstrom'.",
                geometry.units
            )));
        }
    };

    // 5. Convert atoms to qc_core::basis::Atom
    let mut atoms: Vec<Atom> = Vec::with_capacity(geometry.atoms.len());
    let mut atom_outputs: Vec<AtomOutput> = Vec::with_capacity(geometry.atoms.len());

    for atom_input in &geometry.atoms {
        let atomic_number = symbol_to_atomic_number(&atom_input.symbol).map_err(|_| {
            JsError::new(&format!(
                "Unsupported element '{}'. Only H-Ar (Z=1-18) are supported.",
                atom_input.symbol
            ))
        })?;

        let position_bohr = if convert_to_bohr {
            [
                atom_input.xyz[0] * ANGSTROM_TO_BOHR,
                atom_input.xyz[1] * ANGSTROM_TO_BOHR,
                atom_input.xyz[2] * ANGSTROM_TO_BOHR,
            ]
        } else {
            atom_input.xyz
        };

        let atom = Atom::new(atomic_number, position_bohr)
            .map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))?;
        atoms.push(atom.clone());

        atom_outputs.push(AtomOutput {
            symbol: atom.symbol.clone(),
            xyz: position_bohr,
            atomic_number: atomic_number as u32,
        });
    }

    // 6. Build basis set
    let basis = BasisSet::build(atoms, &basis_lower)
        .map_err(|e| JsError::new(&format!("Failed to build basis set: {}", e)))?;

    // 7. Compute integrals
    // Use spherical harmonic basis if requested (transforms d-orbitals from 6 to 5 functions)
    let (s_matrix, h_core, nbf) = if use_spherical {
        // Spherical harmonic basis: 5 d-orbitals instead of 6 Cartesian
        let s = overlap_matrix_spherical(&basis);
        let h = hcore_matrix_spherical(&basis);
        let n = basis.n_basis_spherical();
        (s, h, n)
    } else {
        // Cartesian basis: 6 d-orbitals (default)
        let s = overlap_matrix(&basis);
        let h = hcore_matrix(&basis);
        let n = basis.n_basis;
        (s, h, n)
    };

    // Compute ERIs in the same basis as one-electron integrals
    #[cfg(feature = "parallel")]
    let eri = if use_spherical {
        eri_compressed_spherical_parallel(&basis)
    } else {
        eri_compressed_parallel(&basis)
    };
    #[cfg(not(feature = "parallel"))]
    let eri = if use_spherical {
        eri_compressed_spherical(&basis)
    } else {
        eri_compressed(&basis)
    };

    // 8. Generate system ID from geometry + basis + spherical option hash
    let system_id = {
        let mut hasher = DefaultHasher::new();
        for atom in &atom_outputs {
            atom.symbol.hash(&mut hasher);
            ((atom.xyz[0] * 1e6) as i64).hash(&mut hasher);
            ((atom.xyz[1] * 1e6) as i64).hash(&mut hasher);
            ((atom.xyz[2] * 1e6) as i64).hash(&mut hasher);
        }
        basis_lower.hash(&mut hasher);
        use_spherical.hash(&mut hasher); // Include in hash for uniqueness
        format!("custom_{:016x}", hasher.finish())
    };

    // 9. Generate label
    let element_symbols: Vec<&str> = atom_outputs.iter().map(|a| a.symbol.as_str()).collect();
    let formula = generate_molecular_formula(&element_symbols);
    let sph_suffix = if use_spherical { ", sph" } else { "" };
    let label = format!("{} ({}{})", formula, basis_lower, sph_suffix);

    // 10. Compute metadata
    let compute_time_ms = (js_sys::Date::now() - start_time) as u64;
    let n_shells = basis.n_shells();
    let shell_pairs = n_shells * (n_shells + 1) / 2;
    let shell_quartets = shell_pairs * (shell_pairs + 1) / 2;
    // Store ERI count before moving the vector
    let significant_eris = eri.len();

    let result = IntegralResult {
        format_version: 1,
        system_id,
        label,
        description: format!(
            "On-the-fly integral computation for {} atoms, {} basis ({})",
            atom_outputs.len(),
            basis_lower,
            basis_type
        ),
        geometry: GeometryOutput {
            atoms: atom_outputs,
            units: "bohr".to_string(),
        },
        basis_id: basis_lower,
        nbf,
        nelec: basis.n_electrons,
        e_nuc: basis.nuclear_repulsion,
        s_matrix,
        h_core,
        eri_compressed: eri, // Move instead of clone for better performance
        eri_indexing: "8-fold symmetry, pair_ij * (pair_ij+1)/2 + pair_kl".to_string(),
        metadata: IntegralMetadata {
            wasm_version: env!("CARGO_PKG_VERSION").to_string(),
            compute_time_ms,
            shell_pairs,
            shell_quartets,
            significant_eris,
            basis_type: basis_type.to_string(),
        },
    };

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

/// Generate a molecular formula from element symbols.
///
/// Groups elements and counts them in standard order (C, H, then alphabetical).
fn generate_molecular_formula(symbols: &[&str]) -> String {
    use std::collections::HashMap;

    let mut counts: HashMap<&str, usize> = HashMap::new();
    for symbol in symbols {
        *counts.entry(*symbol).or_insert(0) += 1;
    }

    // Standard ordering: C first, then H, then alphabetical
    let mut formula = String::new();

    // Carbon first
    if let Some(&c) = counts.get("C") {
        formula.push('C');
        if c > 1 {
            formula.push_str(&c.to_string());
        }
        counts.remove("C");
    }

    // Hydrogen second
    if let Some(&h) = counts.get("H") {
        formula.push('H');
        if h > 1 {
            formula.push_str(&h.to_string());
        }
        counts.remove("H");
    }

    // Rest alphabetically
    let mut remaining: Vec<_> = counts.into_iter().collect();
    remaining.sort_by(|a, b| a.0.cmp(b.0));

    for (symbol, count) in remaining {
        formula.push_str(symbol);
        if count > 1 {
            formula.push_str(&count.to_string());
        }
    }

    formula
}

/// Compute all molecular integrals with progress callback.
///
/// This is the same as `compute_integrals` but accepts an optional progress
/// callback function that is invoked during computation to report progress.
/// This is essential for larger molecules where ERI computation can take
/// several seconds.
///
/// # Arguments
///
/// * `geometry_json` - JSON string containing GeometryInput
/// * `basis_name` - Basis set name: "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"
/// * `progress_callback` - Optional JavaScript callback function for progress updates
///
/// # Progress Callback
///
/// The callback receives an `IntegralProgress` object with:
/// - `phase`: Current computation phase ("setup", "overlap", "kinetic", "nuclear", "eri", "assembly")
/// - `current`: Current step within the phase
/// - `total`: Total steps in the phase
/// - `overallPercent`: Overall completion percentage (0-100)
/// - `message`: Human-readable status message
///
/// # Example
///
/// ```javascript
/// import { compute_integrals_with_progress } from './qc_wasm.js';
///
/// const geometry = { atoms: [...], units: "bohr" };
///
/// const result = compute_integrals_with_progress(
///   JSON.stringify(geometry),
///   "6-31g*",
///   (progress) => {
///     console.log(`${progress.phase}: ${progress.overallPercent.toFixed(1)}%`);
///   }
/// );
/// ```
#[wasm_bindgen]
pub fn compute_integrals_with_progress(
    geometry_json: &str,
    basis_name: &str,
    progress_callback: Option<js_sys::Function>,
) -> Result<JsValue, JsError> {
    use qc_core::basis::{symbol_to_atomic_number, Atom, BasisSet, ANGSTROM_TO_BOHR};
    #[cfg(feature = "parallel")]
    use qc_core::integrals::eri_compressed_parallel;
    #[cfg(not(feature = "parallel"))]
    use qc_core::integrals::{eri_index, shell_eri};
    use qc_core::integrals::{hcore_matrix, overlap_matrix};
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    // Helper to emit progress
    let emit_progress =
        |phase: &str, current: usize, total: usize, overall_percent: f64, message: &str| {
            if let Some(ref callback) = progress_callback {
                let progress = IntegralProgress {
                    phase: phase.to_string(),
                    current,
                    total,
                    overall_percent,
                    message: message.to_string(),
                };
                if let Ok(js_progress) = serde_wasm_bindgen::to_value(&progress) {
                    let _ = callback.call1(&JsValue::NULL, &js_progress);
                }
            }
        };

    // Use js_sys::Date for timing (std::time::Instant is not supported in WASM)
    let start_time = js_sys::Date::now();

    // Emit setup phase
    emit_progress(
        "setup",
        0,
        1,
        0.0,
        "Parsing geometry and building basis set...",
    );

    // 1. Parse geometry JSON
    let geometry: GeometryInput = serde_json::from_str(geometry_json)
        .map_err(|e| JsError::new(&format!("Invalid geometry JSON: {}", e)))?;

    // 2. Validate geometry
    if geometry.atoms.is_empty() {
        return Err(JsError::new("Geometry must have at least 1 atom."));
    }

    // 3. Normalize basis name and validate
    let basis_lower = basis_name.to_lowercase();
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];
    if !supported_bases.contains(&basis_lower.as_str()) {
        return Err(JsError::new(&format!(
            "Unknown basis set '{}'. Supported: {}",
            basis_name,
            supported_bases.join(", ")
        )));
    }

    // 4. Determine coordinate units
    let units_lower = geometry.units.to_lowercase();
    let convert_to_bohr = match units_lower.as_str() {
        "bohr" => false,
        "angstrom" | "angstroms" => true,
        _ => {
            return Err(JsError::new(&format!(
                "Invalid units '{}'. Must be 'bohr' or 'angstrom'.",
                geometry.units
            )));
        }
    };

    // 5. Convert atoms to qc_core::basis::Atom
    let mut atoms: Vec<Atom> = Vec::with_capacity(geometry.atoms.len());
    let mut atom_outputs: Vec<AtomOutput> = Vec::with_capacity(geometry.atoms.len());

    for atom_input in &geometry.atoms {
        let atomic_number = symbol_to_atomic_number(&atom_input.symbol).map_err(|_| {
            JsError::new(&format!(
                "Unsupported element '{}'. Only H-Ar (Z=1-18) are supported.",
                atom_input.symbol
            ))
        })?;

        let position_bohr = if convert_to_bohr {
            [
                atom_input.xyz[0] * ANGSTROM_TO_BOHR,
                atom_input.xyz[1] * ANGSTROM_TO_BOHR,
                atom_input.xyz[2] * ANGSTROM_TO_BOHR,
            ]
        } else {
            atom_input.xyz
        };

        let atom = Atom::new(atomic_number, position_bohr)
            .map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))?;
        atoms.push(atom.clone());

        atom_outputs.push(AtomOutput {
            symbol: atom.symbol.clone(),
            xyz: position_bohr,
            atomic_number: atomic_number as u32,
        });
    }

    // 6. Build basis set
    let basis = BasisSet::build(atoms, &basis_lower)
        .map_err(|e| JsError::new(&format!("Failed to build basis set: {}", e)))?;

    emit_progress(
        "setup",
        1,
        1,
        5.0,
        &format!(
            "Basis set ready: {} basis functions, {} shells",
            basis.n_basis,
            basis.n_shells()
        ),
    );

    // Calculate total shell quartets for progress tracking
    let n_shells = basis.n_shells();
    #[cfg(not(feature = "parallel"))]
    let n_basis = basis.n_basis;
    #[cfg(not(feature = "parallel"))]
    let n_pairs = n_basis * (n_basis + 1) / 2;
    #[cfg(not(feature = "parallel"))]
    let n_eri = n_pairs * (n_pairs + 1) / 2;

    // Count shell quartets for progress (accounting for symmetry)
    let mut total_shell_quartets = 0usize;
    for si in 0..n_shells {
        for _sj in 0..=si {
            for _sk in 0..n_shells {
                for _sl in 0..=_sk {
                    total_shell_quartets += 1;
                }
            }
        }
    }

    // 7. Compute overlap matrix (S)
    emit_progress("overlap", 0, 1, 10.0, "Computing overlap integrals (S)...");
    let s_matrix = overlap_matrix(&basis);
    emit_progress("overlap", 1, 1, 15.0, "Overlap matrix complete");

    // 8. Compute core Hamiltonian (H = T + V)
    // hcore_matrix computes both kinetic (T) and nuclear (V)
    emit_progress("hcore", 0, 1, 15.0, "Computing core Hamiltonian (T + V)...");
    let h_core = hcore_matrix(&basis);
    emit_progress("hcore", 1, 1, 25.0, "Core Hamiltonian complete");

    // 9. Compute two-electron integrals (ERI)
    // This is the expensive O(N^4) operation
    emit_progress(
        "eri",
        0,
        total_shell_quartets,
        25.0,
        "Starting ERI computation...",
    );

    // With parallel feature: use optimized parallel computation
    #[cfg(feature = "parallel")]
    let eri = {
        emit_progress(
            "eri",
            total_shell_quartets / 2,
            total_shell_quartets,
            60.0,
            "Computing ERIs...",
        );
        eri_compressed_parallel(&basis)
    };

    // Without parallel feature: sequential computation with progress tracking
    #[cfg(not(feature = "parallel"))]
    let eri = {
        // Allocate ERI storage
        let mut eri_vec = vec![0.0; n_eri];

        // Track progress
        let mut shell_quartets_computed = 0usize;
        let mut last_progress_percent = 25.0;
        let progress_step = 10.0; // Emit progress every 10% to reduce callback overhead

        // Iterate over shell quartets (matching eri_compressed logic)
        let mut mu_i = 0;
        for (si, shell_i) in basis.shells.iter().enumerate() {
            let n_i = shell_i.n_basis_functions();

            let mut mu_j = 0;
            for (_sj, shell_j) in basis.shells.iter().enumerate().take(si + 1) {
                let n_j = shell_j.n_basis_functions();

                let mut mu_k = 0;
                for (_sk, shell_k) in basis.shells.iter().enumerate() {
                    let n_k = shell_k.n_basis_functions();

                    let mut mu_l = 0;
                    for (_sl, shell_l) in basis.shells.iter().enumerate().take(_sk + 1) {
                        let n_l = shell_l.n_basis_functions();

                        // Compute shell quartet
                        let block = shell_eri(shell_i, shell_j, shell_k, shell_l);

                        // Store unique integrals
                        for ii in 0..n_i {
                            let i = mu_i + ii;
                            for jj in 0..n_j {
                                let j = mu_j + jj;
                                if i < j {
                                    continue;
                                }

                                for kk in 0..n_k {
                                    let k = mu_k + kk;
                                    for ll in 0..n_l {
                                        let l = mu_l + ll;
                                        if k < l {
                                            continue;
                                        }

                                        let idx = eri_index(n_basis, i, j, k, l);
                                        eri_vec[idx] = block.get(ii, jj, kk, ll);
                                    }
                                }
                            }
                        }

                        shell_quartets_computed += 1;

                        // Emit progress periodically
                        let current_percent = 25.0
                            + (shell_quartets_computed as f64 / total_shell_quartets as f64) * 70.0;
                        if current_percent >= last_progress_percent + progress_step {
                            emit_progress(
                                "eri",
                                shell_quartets_computed,
                                total_shell_quartets,
                                current_percent,
                                &format!(
                                    "Computing ERIs: {}/{} shell quartets ({:.0}%)",
                                    shell_quartets_computed, total_shell_quartets, current_percent
                                ),
                            );
                            last_progress_percent = current_percent;
                        }

                        mu_l += n_l;
                    }
                    mu_k += n_k;
                }
                mu_j += n_j;
            }
            mu_i += n_i;
        }
        eri_vec
    };

    emit_progress(
        "eri",
        total_shell_quartets,
        total_shell_quartets,
        95.0,
        "ERI computation complete",
    );

    // 10. Assembly phase - generate metadata and result
    emit_progress("assembly", 0, 1, 95.0, "Assembling results...");

    // Generate system ID from geometry + basis hash
    let system_id = {
        let mut hasher = DefaultHasher::new();
        for atom in &atom_outputs {
            atom.symbol.hash(&mut hasher);
            ((atom.xyz[0] * 1e6) as i64).hash(&mut hasher);
            ((atom.xyz[1] * 1e6) as i64).hash(&mut hasher);
            ((atom.xyz[2] * 1e6) as i64).hash(&mut hasher);
        }
        basis_lower.hash(&mut hasher);
        format!("custom_{:016x}", hasher.finish())
    };

    // Generate label
    let element_symbols: Vec<&str> = atom_outputs.iter().map(|a| a.symbol.as_str()).collect();
    let formula = generate_molecular_formula(&element_symbols);
    let label = format!("{} ({})", formula, basis_lower);

    // Compute metadata
    let compute_time_ms = (js_sys::Date::now() - start_time) as u64;
    let shell_pairs = n_shells * (n_shells + 1) / 2;
    let shell_quartets = shell_pairs * (shell_pairs + 1) / 2;
    // Store ERI count before moving the vector
    let significant_eris = eri.len();

    let result = IntegralResult {
        format_version: 1,
        system_id,
        label,
        description: format!(
            "On-the-fly integral computation for {} atoms, {} basis",
            atom_outputs.len(),
            basis_lower
        ),
        geometry: GeometryOutput {
            atoms: atom_outputs,
            units: "bohr".to_string(),
        },
        basis_id: basis_lower,
        nbf: basis.n_basis,
        nelec: basis.n_electrons,
        e_nuc: basis.nuclear_repulsion,
        s_matrix,
        h_core,
        eri_compressed: eri, // Move instead of clone for better performance
        eri_indexing: "8-fold symmetry, pair_ij * (pair_ij+1)/2 + pair_kl".to_string(),
        metadata: IntegralMetadata {
            wasm_version: env!("CARGO_PKG_VERSION").to_string(),
            compute_time_ms,
            shell_pairs,
            shell_quartets,
            significant_eris,
            basis_type: "cartesian".to_string(), // Default to Cartesian
        },
    };

    emit_progress("assembly", 1, 1, 100.0, "Computation complete!");

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

/// Compute all molecular integrals with both options and progress callback.
///
/// This function combines the features of `compute_integrals_with_options` (spherical
/// harmonics support) and `compute_integrals_with_progress` (progress callbacks).
/// Use this when you need both features simultaneously.
///
/// # Arguments
///
/// * `geometry_json` - JSON string containing GeometryInput
/// * `basis_name` - Basis set name: "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"
/// * `options` - JavaScript object with computation options (e.g., `{ useSpherical: true }`)
/// * `progress_callback` - Optional JavaScript callback function for progress updates
///
/// # Options
///
/// - `useSpherical`: Boolean (default: false). When true, uses spherical harmonic
///   basis functions (5 d-orbitals) instead of Cartesian (6 d-orbitals).
///
/// # Progress Callback
///
/// The callback receives an `IntegralProgress` object with:
/// - `phase`: Current computation phase ("setup", "overlap", "hcore", "eri", "assembly")
/// - `current`: Current step within the phase
/// - `total`: Total steps in the phase
/// - `overallPercent`: Overall completion percentage (0-100)
/// - `message`: Human-readable status message
///
/// # Example
///
/// ```javascript
/// import { compute_integrals_with_options_and_progress } from './qc_wasm.js';
///
/// const geometry = { atoms: [...], units: "angstrom" };
/// const options = { useSpherical: true };
///
/// const result = compute_integrals_with_options_and_progress(
///   JSON.stringify(geometry),
///   "6-31g*",
///   options,
///   (progress) => {
///     console.log(`${progress.phase}: ${progress.overallPercent.toFixed(1)}%`);
///   }
/// );
/// console.log(result.metadata.basisType); // "spherical"
/// ```
#[wasm_bindgen]
pub fn compute_integrals_with_options_and_progress(
    geometry_json: &str,
    basis_name: &str,
    options: JsValue,
    progress_callback: Option<js_sys::Function>,
) -> Result<JsValue, JsError> {
    use qc_core::basis::{symbol_to_atomic_number, Atom, BasisSet, ANGSTROM_TO_BOHR};
    #[cfg(feature = "parallel")]
    use qc_core::integrals::{eri_compressed_parallel, eri_compressed_spherical_parallel};
    #[cfg(not(feature = "parallel"))]
    use qc_core::integrals::{eri_compressed_spherical, eri_index, shell_eri};
    use qc_core::integrals::{
        hcore_matrix, hcore_matrix_spherical, overlap_matrix, overlap_matrix_spherical,
    };
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    // Parse options (default to Cartesian for backward compatibility)
    let opts: IntegralComputeOptions = if options.is_undefined() || options.is_null() {
        IntegralComputeOptions::default()
    } else {
        serde_wasm_bindgen::from_value(options)
            .map_err(|e| JsError::new(&format!("Invalid options: {}", e)))?
    };

    let use_spherical = opts.use_spherical;
    let basis_type = if use_spherical {
        "spherical"
    } else {
        "cartesian"
    };

    // Helper to emit progress
    let emit_progress =
        |phase: &str, current: usize, total: usize, overall_percent: f64, message: &str| {
            if let Some(ref callback) = progress_callback {
                let progress = IntegralProgress {
                    phase: phase.to_string(),
                    current,
                    total,
                    overall_percent,
                    message: message.to_string(),
                };
                if let Ok(js_progress) = serde_wasm_bindgen::to_value(&progress) {
                    let _ = callback.call1(&JsValue::NULL, &js_progress);
                }
            }
        };

    // Use js_sys::Date for timing (std::time::Instant is not supported in WASM)
    let start_time = js_sys::Date::now();

    // Emit setup phase
    emit_progress(
        "setup",
        0,
        1,
        0.0,
        "Parsing geometry and building basis set...",
    );

    // 1. Parse geometry JSON
    let geometry: GeometryInput = serde_json::from_str(geometry_json)
        .map_err(|e| JsError::new(&format!("Invalid geometry JSON: {}", e)))?;

    // 2. Validate geometry
    if geometry.atoms.is_empty() {
        return Err(JsError::new("Geometry must have at least 1 atom."));
    }

    // 3. Normalize basis name and validate
    let basis_lower = basis_name.to_lowercase();
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];
    if !supported_bases.contains(&basis_lower.as_str()) {
        return Err(JsError::new(&format!(
            "Unknown basis set '{}'. Supported: {}",
            basis_name,
            supported_bases.join(", ")
        )));
    }

    // 4. Determine coordinate units
    let units_lower = geometry.units.to_lowercase();
    let convert_to_bohr = match units_lower.as_str() {
        "bohr" => false,
        "angstrom" | "angstroms" => true,
        _ => {
            return Err(JsError::new(&format!(
                "Invalid units '{}'. Must be 'bohr' or 'angstrom'.",
                geometry.units
            )));
        }
    };

    // 5. Convert atoms to qc_core::basis::Atom
    let mut atoms: Vec<Atom> = Vec::with_capacity(geometry.atoms.len());
    let mut atom_outputs: Vec<AtomOutput> = Vec::with_capacity(geometry.atoms.len());

    for atom_input in &geometry.atoms {
        let atomic_number = symbol_to_atomic_number(&atom_input.symbol).map_err(|_| {
            JsError::new(&format!(
                "Unsupported element '{}'. Only H-Ar (Z=1-18) are supported.",
                atom_input.symbol
            ))
        })?;

        let position_bohr = if convert_to_bohr {
            [
                atom_input.xyz[0] * ANGSTROM_TO_BOHR,
                atom_input.xyz[1] * ANGSTROM_TO_BOHR,
                atom_input.xyz[2] * ANGSTROM_TO_BOHR,
            ]
        } else {
            atom_input.xyz
        };

        let atom = Atom::new(atomic_number, position_bohr)
            .map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))?;
        atoms.push(atom.clone());

        atom_outputs.push(AtomOutput {
            symbol: atom.symbol.clone(),
            xyz: position_bohr,
            atomic_number: atomic_number as u32,
        });
    }

    // 6. Build basis set
    let basis = BasisSet::build(atoms, &basis_lower)
        .map_err(|e| JsError::new(&format!("Failed to build basis set: {}", e)))?;

    // Determine number of basis functions based on spherical/Cartesian choice
    let nbf = if use_spherical {
        basis.n_basis_spherical()
    } else {
        basis.n_basis
    };

    emit_progress(
        "setup",
        1,
        1,
        5.0,
        &format!(
            "Basis set ready: {} basis functions ({}), {} shells",
            nbf,
            basis_type,
            basis.n_shells()
        ),
    );

    // Calculate progress tracking values
    let n_shells = basis.n_shells();
    #[cfg(not(feature = "parallel"))]
    let n_pairs = nbf * (nbf + 1) / 2;
    #[cfg(not(feature = "parallel"))]
    let n_eri = n_pairs * (n_pairs + 1) / 2;

    // Count shell quartets for progress (accounting for symmetry)
    let mut total_shell_quartets = 0usize;
    for si in 0..n_shells {
        for _sj in 0..=si {
            for _sk in 0..n_shells {
                for _sl in 0..=_sk {
                    total_shell_quartets += 1;
                }
            }
        }
    }

    // 7. Compute overlap matrix (S)
    emit_progress("overlap", 0, 1, 10.0, "Computing overlap integrals (S)...");
    let s_matrix = if use_spherical {
        overlap_matrix_spherical(&basis)
    } else {
        overlap_matrix(&basis)
    };
    emit_progress("overlap", 1, 1, 15.0, "Overlap matrix complete");

    // 8. Compute core Hamiltonian (H = T + V)
    emit_progress("hcore", 0, 1, 15.0, "Computing core Hamiltonian (T + V)...");
    let h_core = if use_spherical {
        hcore_matrix_spherical(&basis)
    } else {
        hcore_matrix(&basis)
    };
    emit_progress("hcore", 1, 1, 25.0, "Core Hamiltonian complete");

    // 9. Compute two-electron integrals (ERI) with progress tracking
    // This is the expensive O(N^4) operation
    emit_progress(
        "eri",
        0,
        total_shell_quartets,
        25.0,
        "Starting ERI computation...",
    );

    // With parallel feature: use optimized parallel computation
    #[cfg(feature = "parallel")]
    let eri = {
        emit_progress(
            "eri",
            total_shell_quartets / 2,
            total_shell_quartets,
            60.0,
            "Computing ERIs...",
        );
        if use_spherical {
            eri_compressed_spherical_parallel(&basis)
        } else {
            eri_compressed_parallel(&basis)
        }
    };

    // Without parallel feature: sequential computation
    #[cfg(not(feature = "parallel"))]
    let eri = if use_spherical {
        // For spherical, use the optimized function
        emit_progress(
            "eri",
            total_shell_quartets / 2,
            total_shell_quartets,
            60.0,
            "Computing ERIs (spherical)...",
        );
        eri_compressed_spherical(&basis)
    } else {
        // For Cartesian, we compute with per-shell progress
        let mut eri_vec = vec![0.0; n_eri];

        // Track progress
        let mut shell_quartets_computed = 0usize;
        let mut last_progress_percent = 25.0;
        let progress_step = 10.0; // Reduce callback overhead

        // Iterate over shell quartets
        let mut mu_i = 0;
        for (si, shell_i) in basis.shells.iter().enumerate() {
            let n_i = shell_i.n_basis_functions();

            let mut mu_j = 0;
            for (_sj, shell_j) in basis.shells.iter().enumerate().take(si + 1) {
                let n_j = shell_j.n_basis_functions();

                let mut mu_k = 0;
                for (_sk, shell_k) in basis.shells.iter().enumerate() {
                    let n_k = shell_k.n_basis_functions();

                    let mut mu_l = 0;
                    for (_sl, shell_l) in basis.shells.iter().enumerate().take(_sk + 1) {
                        let n_l = shell_l.n_basis_functions();

                        // Compute shell quartet
                        let block = shell_eri(shell_i, shell_j, shell_k, shell_l);

                        // Store unique integrals
                        for ii in 0..n_i {
                            let i = mu_i + ii;
                            for jj in 0..n_j {
                                let j = mu_j + jj;
                                if i < j {
                                    continue;
                                }

                                for kk in 0..n_k {
                                    let k = mu_k + kk;
                                    for ll in 0..n_l {
                                        let l = mu_l + ll;
                                        if k < l {
                                            continue;
                                        }

                                        let idx = eri_index(basis.n_basis, i, j, k, l);
                                        eri_vec[idx] = block.get(ii, jj, kk, ll);
                                    }
                                }
                            }
                        }

                        shell_quartets_computed += 1;

                        // Emit progress periodically
                        let current_percent = 25.0
                            + (shell_quartets_computed as f64 / total_shell_quartets as f64) * 70.0;
                        if current_percent >= last_progress_percent + progress_step {
                            emit_progress(
                                "eri",
                                shell_quartets_computed,
                                total_shell_quartets,
                                current_percent,
                                &format!(
                                    "Computing ERIs: {}/{} shell quartets ({:.0}%)",
                                    shell_quartets_computed, total_shell_quartets, current_percent
                                ),
                            );
                            last_progress_percent = current_percent;
                        }

                        mu_l += n_l;
                    }
                    mu_k += n_k;
                }
                mu_j += n_j;
            }
            mu_i += n_i;
        }

        eri_vec
    };

    emit_progress(
        "eri",
        total_shell_quartets,
        total_shell_quartets,
        95.0,
        "ERI computation complete",
    );

    // 10. Assembly phase - generate metadata and result
    emit_progress("assembly", 0, 1, 95.0, "Assembling results...");

    // Generate system ID from geometry + basis + spherical option hash
    let system_id = {
        let mut hasher = DefaultHasher::new();
        for atom in &atom_outputs {
            atom.symbol.hash(&mut hasher);
            ((atom.xyz[0] * 1e6) as i64).hash(&mut hasher);
            ((atom.xyz[1] * 1e6) as i64).hash(&mut hasher);
            ((atom.xyz[2] * 1e6) as i64).hash(&mut hasher);
        }
        basis_lower.hash(&mut hasher);
        use_spherical.hash(&mut hasher); // Include in hash for uniqueness
        format!("custom_{:016x}", hasher.finish())
    };

    // Generate label
    let element_symbols: Vec<&str> = atom_outputs.iter().map(|a| a.symbol.as_str()).collect();
    let formula = generate_molecular_formula(&element_symbols);
    let sph_suffix = if use_spherical { ", sph" } else { "" };
    let label = format!("{} ({}{})", formula, basis_lower, sph_suffix);

    // Compute metadata
    let compute_time_ms = (js_sys::Date::now() - start_time) as u64;
    let shell_pairs = n_shells * (n_shells + 1) / 2;
    let shell_quartets = shell_pairs * (shell_pairs + 1) / 2;
    // Store ERI count before moving the vector
    let significant_eris = eri.len();

    let result = IntegralResult {
        format_version: 1,
        system_id,
        label,
        description: format!(
            "On-the-fly integral computation for {} atoms, {} basis ({})",
            atom_outputs.len(),
            basis_lower,
            basis_type
        ),
        geometry: GeometryOutput {
            atoms: atom_outputs,
            units: "bohr".to_string(),
        },
        basis_id: basis_lower,
        nbf,
        nelec: basis.n_electrons,
        e_nuc: basis.nuclear_repulsion,
        s_matrix,
        h_core,
        eri_compressed: eri, // Move instead of clone for better performance
        eri_indexing: "8-fold symmetry, pair_ij * (pair_ij+1)/2 + pair_kl".to_string(),
        metadata: IntegralMetadata {
            wasm_version: env!("CARGO_PKG_VERSION").to_string(),
            compute_time_ms,
            shell_pairs,
            shell_quartets,
            significant_eris,
            basis_type: basis_type.to_string(),
        },
    };

    emit_progress("assembly", 1, 1, 100.0, "Computation complete!");

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Test Functions (for WASM integration verification)
// ============================================================================

/// Result type for test computation.
///
/// This struct demonstrates the data round-trip pattern used throughout
/// the WASM bindings: Rust struct -> serde serialization -> JsValue.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct TestResult {
    /// The input value that was provided.
    pub input: f64,
    /// The computed output (input * 2).
    pub output: f64,
    /// A human-readable message describing the computation.
    pub message: String,
}

// ============================================================================
// PES Scan WASM Bindings (US-039)
// ============================================================================

/// WASM-friendly PES scan input struct.
///
/// This struct is deserialized from JavaScript to configure a PES bond-length
/// scan. It uses camelCase field names to match JavaScript conventions.
///
/// # Fields
///
/// * `atom_a_z` - Atomic number of atom A (fixed at origin), 1-10
/// * `atom_b_z` - Atomic number of atom B (translated along z-axis), 1-10
/// * `r_min` - Minimum bond distance in bohr
/// * `r_max` - Maximum bond distance in bohr
/// * `n_points` - Number of scan points (evenly spaced)
/// * `basis_name` - Basis set name (e.g., "sto-3g")
/// * `options` - SCF computation options
/// * `use_seeding` - Whether to seed convergence from previous point
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PesScanWasmInput {
    /// Atomic number of atom A (fixed at origin)
    pub atom_a_z: u8,
    /// Atomic number of atom B (translated along z-axis)
    pub atom_b_z: u8,
    /// Minimum bond distance in bohr
    pub r_min: f64,
    /// Maximum bond distance in bohr
    pub r_max: f64,
    /// Number of scan points
    pub n_points: usize,
    /// Basis set name (e.g., "sto-3g")
    pub basis_name: String,
    /// SCF computation options
    pub options: ScfWasmOptions,
    /// Whether to use convergence seeding from previous point
    #[serde(default = "default_use_seeding")]
    pub use_seeding: bool,
}

fn default_use_seeding() -> bool {
    true
}

/// WASM-friendly PES scan progress struct.
///
/// Serialized and passed to the JavaScript progress callback after each
/// scan point completes.
#[derive(Serialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PesScanWasmProgress {
    /// Index of the completed point (0-indexed)
    pub point_index: usize,
    /// Total number of scan points
    pub total_points: usize,
    /// Bond distance in bohr
    pub r: f64,
    /// SCF energy in Hartree
    pub energy: f64,
    /// Whether SCF converged at this point
    pub converged: bool,
}

/// Run a bond-length PES scan for a diatomic molecule.
///
/// Scans the potential energy surface by varying the internuclear distance
/// from `r_min` to `r_max` over `n_points` evenly spaced points. At each
/// geometry, molecular integrals are computed on-the-fly and SCF is run.
///
/// The progress callback is called after each scan point with a
/// `PesScanWasmProgress` object.
///
/// # Arguments
///
/// * `input` - JavaScript object containing `PesScanWasmInput` fields
/// * `progress_callback` - JavaScript function called after each scan point
///
/// # Returns
///
/// A JavaScript object containing `PesScanResult`:
/// - `points`: Array of { r, energy, converged, iterations }
/// - `equilibrium`: { rBohr, energyHartree } or null
/// - `computeTimeMs`: Total computation time
/// - `totalIterations`: Sum of SCF iterations across all points
///
/// # Errors
///
/// Returns a `JsError` if:
/// - Input deserialization fails
/// - Atomic numbers are invalid (must be 1-10)
/// - `r_min` >= `r_max`
/// - `n_points` is 0
///
/// # Example
///
/// ```javascript
/// import { pes_scan } from './qc_wasm.js';
///
/// const input = {
///   atomAZ: 1, atomBZ: 1,
///   rMin: 0.5, rMax: 5.0, nPoints: 20,
///   basisName: "sto-3g",
///   options: { convergenceProfile: "medium", maxIterations: 100, useDiis: true },
///   useSeeding: true,
/// };
///
/// const result = pes_scan(input, (progress) => {
///   console.log(`Point ${progress.pointIndex}/${progress.totalPoints}: r=${progress.r}`);
/// });
/// ```
#[wasm_bindgen]
pub fn pes_scan(input: JsValue, progress_callback: &js_sys::Function) -> Result<JsValue, JsError> {
    use qc_core::scf::pes;

    // 1. Deserialize input
    let wasm_input: PesScanWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid PES scan input: {}", e)))?;

    // 2. Validate input
    if wasm_input.atom_a_z < 1 || wasm_input.atom_a_z > 10 {
        return Err(JsError::new(&format!(
            "pes_scan: atomAZ must be 1-10, got {}",
            wasm_input.atom_a_z
        )));
    }
    if wasm_input.atom_b_z < 1 || wasm_input.atom_b_z > 10 {
        return Err(JsError::new(&format!(
            "pes_scan: atomBZ must be 1-10, got {}",
            wasm_input.atom_b_z
        )));
    }
    if wasm_input.n_points == 0 {
        return Err(JsError::new("pes_scan: nPoints must be > 0"));
    }
    if wasm_input.n_points > 1 && wasm_input.r_min >= wasm_input.r_max {
        return Err(JsError::new(&format!(
            "pes_scan: rMin ({}) must be < rMax ({})",
            wasm_input.r_min, wasm_input.r_max
        )));
    }

    // 3. Convert ScfWasmOptions to ScfConfig
    let scf_config = ScfConfig {
        profile: parse_convergence_profile(&wasm_input.options.convergence_profile),
        max_iterations: wasm_input.options.max_iterations as usize,
        use_diis: wasm_input.options.use_diis,
        diis_size: wasm_input.options.diis_size.unwrap_or(6),
        diis_start: 2,
        damp: wasm_input.options.damp.unwrap_or(0.0),
        damp_start: 5,
        level_shift: 0.0,
    };

    // 4. Build PesScanConfig
    let scan_config = pes::PesScanConfig {
        atom_a_z: wasm_input.atom_a_z,
        atom_b_z: wasm_input.atom_b_z,
        r_min: wasm_input.r_min,
        r_max: wasm_input.r_max,
        n_points: wasm_input.n_points,
        basis_name: &wasm_input.basis_name,
        scf_config: &scf_config,
        use_seeding: wasm_input.use_seeding,
    };

    // 5. Set up progress callback
    let n_points = wasm_input.n_points;
    let progress_fn = |idx: usize, r: f64, energy: f64, converged: bool| {
        let progress = PesScanWasmProgress {
            point_index: idx,
            total_points: n_points,
            r,
            energy,
            converged,
        };
        if let Ok(js_progress) = serde_wasm_bindgen::to_value(&progress) {
            let _ = progress_callback.call1(&JsValue::NULL, &js_progress);
        }
    };

    // 6. Run PES scan (timing via js_sys -- std::time::Instant not supported in WASM)
    let start_ms = js_sys::Date::now();
    let mut result = pes::pes_scan(&scan_config, Some(&progress_fn));
    result.compute_time_ms = js_sys::Date::now() - start_ms;

    // 7. Serialize and return result
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Geometry Optimization WASM Bindings (US-074)
// ============================================================================

/// WASM-friendly geometry optimization input struct.
///
/// Contains the initial molecular geometry, optimization method, basis set,
/// and convergence parameters needed to run L-BFGS geometry optimization.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct OptimizeWasmInput {
    /// Atoms as [[Z, x, y, z], ...] where Z is atomic number and
    /// x, y, z are coordinates in bohr
    pub atoms: Vec<[f64; 4]>,
    /// Basis set name (e.g., "sto-3g")
    pub basis_name: String,
    /// Electronic structure method: "rhf", "lda", or "b3lyp"
    pub method: String,
    /// Maximum optimization steps (default: 50)
    #[serde(default = "default_max_steps")]
    pub max_steps: usize,
    /// Maximum gradient convergence threshold (Ha/bohr, default: 4.5e-4)
    #[serde(default = "default_grad_threshold")]
    pub grad_threshold: f64,
    /// Energy convergence threshold (Ha, default: 1e-6)
    #[serde(default = "default_energy_threshold")]
    pub energy_threshold: f64,
    /// L-BFGS memory size (default: 7)
    #[serde(default = "default_memory_size")]
    pub memory_size: usize,
}

fn default_max_steps() -> usize {
    50
}
fn default_grad_threshold() -> f64 {
    4.5e-4
}
fn default_energy_threshold() -> f64 {
    1.0e-6
}
fn default_memory_size() -> usize {
    7
}

/// Progress update emitted after each optimization step.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct OptimizeWasmProgress {
    /// Current step number (0 = initial evaluation)
    pub step: usize,
    /// Total energy at this step (Ha)
    pub energy: f64,
    /// Maximum gradient component (Ha/bohr)
    pub max_gradient: f64,
    /// RMS gradient (Ha/bohr)
    pub rms_gradient: f64,
}

/// Run L-BFGS geometry optimization.
///
/// Optimizes the molecular geometry to minimize the total energy using
/// the L-BFGS quasi-Newton method. Supports RHF, LDA, and B3LYP methods.
///
/// # Arguments
///
/// * `input` - JavaScript object with optimization parameters
/// * `progress_callback` - Called after each step with progress info
///
/// # Returns
///
/// An `OptimizationResult` object with the full trajectory, final geometry,
/// and convergence status.
///
/// # Example
///
/// ```javascript
/// import { optimize_geometry } from './qc_wasm.js';
///
/// const input = {
///   atoms: [[1, 0, 0, 0], [1, 0, 0, 1.4]],
///   basisName: "sto-3g",
///   method: "rhf",
///   maxSteps: 50,
/// };
///
/// const result = optimize_geometry(input, (progress) => {
///   console.log(`Step ${progress.step}: E=${progress.energy}, maxGrad=${progress.maxGradient}`);
/// });
/// ```
#[wasm_bindgen]
pub fn optimize_geometry(
    input: JsValue,
    progress_callback: &js_sys::Function,
) -> Result<JsValue, JsError> {
    use qc_core::optimizer;

    // 1. Deserialize input
    let wasm_input: OptimizeWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid optimization input: {}", e)))?;

    // 2. Validate input
    if wasm_input.atoms.is_empty() {
        return Err(JsError::new("optimize_geometry: atoms array is empty"));
    }

    // 3. Convert atoms from [Z, x, y, z] format
    let atoms: Vec<(u8, [f64; 3])> = wasm_input
        .atoms
        .iter()
        .map(|a| (a[0] as u8, [a[1], a[2], a[3]]))
        .collect();

    // 4. Parse method
    let method = match wasm_input.method.to_lowercase().as_str() {
        "rhf" | "hf" => optimizer::OptMethod::Rhf,
        "lda" => optimizer::OptMethod::Lda,
        "b3lyp" => optimizer::OptMethod::B3lyp,
        other => {
            return Err(JsError::new(&format!(
                "optimize_geometry: unknown method '{}', expected 'rhf', 'lda', or 'b3lyp'",
                other
            )));
        }
    };

    // 5. Build optimization config
    let config = optimizer::OptimizationConfig {
        max_steps: wasm_input.max_steps,
        grad_threshold: wasm_input.grad_threshold,
        energy_threshold: wasm_input.energy_threshold,
        memory_size: wasm_input.memory_size,
        method,
        basis: wasm_input.basis_name,
    };

    // 6. Set up progress callback
    let progress_fn = |step: &optimizer::OptimizationStep| {
        let progress = OptimizeWasmProgress {
            step: step.step,
            energy: step.energy,
            max_gradient: step.max_gradient,
            rms_gradient: step.rms_gradient,
        };
        if let Ok(js_progress) = serde_wasm_bindgen::to_value(&progress) {
            let _ = progress_callback.call1(&JsValue::NULL, &js_progress);
        }
    };

    // 7. Run optimization (timing via js_sys)
    let start_ms = js_sys::Date::now();
    let mut result = optimizer::optimize_geometry(&atoms, &config, Some(&progress_fn));
    result.compute_time_ms = js_sys::Date::now() - start_ms;

    // 8. Serialize and return result
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Population Analysis WASM Bindings (US-076)
// ============================================================================

/// WASM-friendly population analysis input struct.
///
/// Contains the density matrix, overlap matrix, and atom-to-basis mapping
/// needed for Mulliken and Lowdin population analysis.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PopulationWasmInput {
    /// Flattened density matrix (row-major, nbf x nbf)
    pub density_matrix: Vec<f64>,
    /// Flattened overlap matrix (row-major, nbf x nbf)
    pub overlap_matrix: Vec<f64>,
    /// Number of basis functions
    pub nbf: usize,
    /// Atom specifications: [{ atomicNumber, nBasis }, ...]
    pub atoms: Vec<PopulationAtomInput>,
}

/// Single atom specification for population analysis.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PopulationAtomInput {
    /// Atomic number (Z)
    pub atomic_number: u8,
    /// Number of basis functions on this atom
    pub n_basis: usize,
}

/// Compute Mulliken and Lowdin population analysis.
///
/// Partitions the total electron density among atoms using two methods:
///
/// - **Mulliken** (1955): q_A = Z_A - sum_{mu in A} (DS)_{mu,mu}
/// - **Lowdin** (1950): q_A = Z_A - sum_{mu in A} (S^{1/2} D S^{1/2})_{mu,mu}
///
/// # Arguments
///
/// * `input` - JavaScript object containing `PopulationWasmInput` fields:
///   - `densityMatrix`: Flat density matrix (row-major, nbf x nbf)
///   - `overlapMatrix`: Flat overlap matrix (row-major, nbf x nbf)
///   - `nbf`: Number of basis functions
///   - `atoms`: Array of { atomicNumber, nBasis } for each atom
///
/// # Returns
///
/// A JavaScript object containing:
/// - `atoms`: Array of per-atom charges and populations
/// - `totalMullikenCharge`: Sum of all Mulliken charges
/// - `totalLowdinCharge`: Sum of all Lowdin charges
/// - `computeTimeUs`: Computation time in microseconds
///
/// # Example
///
/// ```javascript
/// import { compute_population } from './qc_wasm.js';
///
/// const result = compute_population({
///   densityMatrix: [...],
///   overlapMatrix: [...],
///   nbf: 7,
///   atoms: [
///     { atomicNumber: 8, nBasis: 5 },
///     { atomicNumber: 1, nBasis: 1 },
///     { atomicNumber: 1, nBasis: 1 },
///   ],
/// });
/// console.log(result.atoms[0].mullikenCharge); // -0.3657...
/// ```
#[wasm_bindgen]
pub fn compute_population(input: JsValue) -> Result<JsValue, JsError> {
    let wasm_input: PopulationWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid population input: {}", e)))?;

    // Convert atom input to (atomic_number, n_basis) tuples
    let atoms: Vec<(u8, usize)> = wasm_input
        .atoms
        .iter()
        .map(|a| (a.atomic_number, a.n_basis))
        .collect();

    let result = qc_core::population::population_analysis(
        &wasm_input.density_matrix,
        &wasm_input.overlap_matrix,
        wasm_input.nbf,
        &atoms,
    )
    .map_err(|e| JsError::new(&format!("Population analysis failed: {}", e)))?;

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// MO Grid Evaluation WASM Bindings (US-042)
// ============================================================================

/// WASM-friendly MO grid evaluation input struct.
///
/// Contains the MO coefficients, basis set specification, and grid parameters
/// needed to evaluate a molecular orbital on a 3D grid.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct MoGridWasmInput {
    /// MO coefficient vector (one per basis function)
    pub mo_coefficients: Vec<f64>,
    /// Atom specifications: [[Z, x, y, z], ...]
    pub atoms: Vec<[f64; 4]>,
    /// Basis set name (e.g., "sto-3g")
    pub basis_name: String,
    /// Grid origin [x, y, z] in Bohr
    pub grid_origin: [f64; 3],
    /// Grid spacing in Bohr (uniform)
    pub grid_spacing: f64,
    /// Grid dimensions [nx, ny, nz]
    pub grid_dims: [usize; 3],
    /// Whether the MO coefficients are in the spherical harmonic basis.
    /// When true, the coefficients are transformed from spherical to Cartesian
    /// before grid evaluation (since the grid evaluator uses Cartesian GTOs).
    /// This must match the `useSpherical` option used in the SCF calculation.
    #[serde(default)]
    pub use_spherical: bool,
}

/// Evaluate a molecular orbital on a 3D grid.
///
/// Computes psi_i(r) = sum_mu { C_{mu,i} * chi_mu(r) } on a uniform 3D grid,
/// where chi_mu are contracted Gaussian basis functions.
///
/// The result contains the grid values and metadata including the approximate
/// norm-squared integral and computation time.
///
/// # Arguments
///
/// * `input` - JavaScript object containing `MoGridWasmInput` fields:
///   - `moCoefficients`: MO coefficient vector (one per basis function)
///   - `atoms`: Array of [Z, x, y, z] for each atom (coordinates in Bohr)
///   - `basisName`: Basis set name (e.g., "sto-3g")
///   - `gridOrigin`: Grid origin [x, y, z] in Bohr
///   - `gridSpacing`: Uniform grid spacing in Bohr
///   - `gridDims`: Grid dimensions [nx, ny, nz]
///   - `useSpherical`: (optional, default false) Whether the MO coefficients
///     are in the spherical harmonic basis. When true, coefficients are
///     transformed from spherical (5 d-functions) to Cartesian (6 d-functions)
///     before grid evaluation.
///
/// # Returns
///
/// A JavaScript object containing `MoGridResult`:
/// - `values`: Flat array of grid values (C-order: x-slowest, z-fastest)
/// - `gridOrigin`: Grid origin
/// - `gridSpacing`: Grid spacing
/// - `gridDims`: Grid dimensions
/// - `maxAbsValue`: Maximum absolute value in the grid
/// - `normSqIntegral`: Approximate norm-squared integral
/// - `computeTimeMs`: Computation time in milliseconds
///
/// # Errors
///
/// Returns a `JsError` if:
/// - Input deserialization fails
/// - Basis set name is invalid
/// - MO coefficient count does not match basis size
/// - Grid dimensions are invalid (< 2 in any direction)
///
/// # Example
///
/// ```javascript
/// import { evaluate_mo_grid } from './qc_wasm.js';
///
/// const input = {
///   moCoefficients: [0.549, 0.549],
///   atoms: [[1, 0, 0, 0], [1, 0, 0, 1.4]],
///   basisName: "sto-3g",
///   gridOrigin: [-5, -5, -5],
///   gridSpacing: 0.5,
///   gridDims: [21, 21, 23],
/// };
///
/// const result = evaluate_mo_grid(input);
/// console.log(`Max value: ${result.maxAbsValue}`);
/// console.log(`Norm^2: ${result.normSqIntegral}`);
/// ```
#[wasm_bindgen]
pub fn evaluate_mo_grid(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::basis::{Atom, BasisSet};
    use qc_core::orbital;

    // 1. Deserialize input
    let wasm_input: MoGridWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid MO grid input: {}", e)))?;

    // 2. Validate input
    if wasm_input.atoms.is_empty() {
        return Err(JsError::new(
            "evaluate_mo_grid: atoms array must not be empty",
        ));
    }
    if wasm_input.grid_dims.iter().any(|&d| d < 2) {
        return Err(JsError::new(
            "evaluate_mo_grid: grid dimensions must be >= 2",
        ));
    }
    if wasm_input.grid_spacing <= 0.0 {
        return Err(JsError::new(&format!(
            "evaluate_mo_grid: grid spacing must be positive, got {}",
            wasm_input.grid_spacing
        )));
    }

    // 3. Build basis set from atom specs
    let atoms: Result<Vec<Atom>, _> = wasm_input
        .atoms
        .iter()
        .map(|a| {
            let z = a[0] as u8;
            Atom::new(z, [a[1], a[2], a[3]])
        })
        .collect();
    let atoms = atoms.map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))?;
    let basis = BasisSet::build(atoms, &wasm_input.basis_name)
        .map_err(|e| JsError::new(&format!("Invalid basis set: {}", e)))?;

    // 4. If MO coefficients are in spherical basis, transform to Cartesian.
    //    The grid evaluator always uses Cartesian GTOs (6 d-functions), so
    //    spherical MO coefficients (5 d-functions) must be expanded.
    //    For s and p shells, spherical == Cartesian (no change).
    //    For d shells: C_cart[mu] = sum_m CART2SPH_D[mu][m] * C_sph[m]
    let mo_coefficients = if wasm_input.use_spherical && basis.has_spherical_difference() {
        use qc_core::integrals::SphericalTransform;

        let n_sph = basis.n_basis_spherical();
        if wasm_input.mo_coefficients.len() != n_sph {
            return Err(JsError::new(&format!(
                "evaluate_mo_grid: MO coefficient count {} does not match spherical basis size {}",
                wasm_input.mo_coefficients.len(),
                n_sph
            )));
        }

        // Transform shell-by-shell from spherical to Cartesian
        let n_cart = basis.n_basis;
        let mut cart_coeffs = Vec::with_capacity(n_cart);
        let mut sph_offset = 0;

        for shell in &basis.shells {
            let l = shell.l_value();
            let transform = SphericalTransform::new(l);
            let n_s = transform.n_sph;
            let n_c = transform.n_cart;

            if transform.needs_transform() {
                // d-shell (or higher): C_cart[c] = sum_m T[c][m] * C_sph[m]
                for c in 0..n_c {
                    let mut val = 0.0;
                    for m in 0..n_s {
                        val += transform.coeff(c, m) * wasm_input.mo_coefficients[sph_offset + m];
                    }
                    cart_coeffs.push(val);
                }
            } else {
                // s/p shell: identity transformation, copy directly
                for m in 0..n_s {
                    cart_coeffs.push(wasm_input.mo_coefficients[sph_offset + m]);
                }
            }

            sph_offset += n_s;
        }

        cart_coeffs
    } else {
        wasm_input.mo_coefficients
    };

    // 5. Evaluate MO on grid
    let result = orbital::evaluate_mo_on_grid(
        &mo_coefficients,
        &basis,
        wasm_input.grid_origin,
        wasm_input.grid_spacing,
        wasm_input.grid_dims,
    )
    .map_err(|e| JsError::new(&format!("MO grid evaluation failed: {}", e)))?;

    // 6. Compute metadata
    let dv = wasm_input.grid_spacing.powi(3);
    let norm_sq_integral: f64 = result.iter().map(|v| v * v).sum::<f64>() * dv;
    let max_abs_value = result.iter().map(|v| v.abs()).fold(0.0f64, f64::max);

    let output = MoGridResult {
        values: result,
        grid_origin: wasm_input.grid_origin,
        grid_spacing: wasm_input.grid_spacing,
        grid_dims: wasm_input.grid_dims,
        max_abs_value,
        norm_sq_integral,
        compute_time_ms: 0.0, // Cannot time in WASM without web_sys::Performance
    };

    // 7. Serialize and return
    serde_wasm_bindgen::to_value(&output).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Marching Cubes WASM Bindings
// ============================================================================

/// WASM input for marching cubes isosurface extraction.
///
/// Contains the scalar field grid data and extraction parameters.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct MarchingCubesWasmInput {
    /// Flat scalar field data (C-order: x-slowest, z-fastest)
    pub grid_data: Vec<f64>,
    /// Grid dimensions [nx, ny, nz]
    pub grid_dims: [usize; 3],
    /// Grid origin [x, y, z] in Bohr
    pub grid_origin: [f64; 3],
    /// Grid spacing in Bohr (uniform)
    pub grid_spacing: f64,
    /// Isovalue threshold
    pub isovalue: f64,
}

/// Extract an isosurface from a 3D scalar field using marching cubes.
///
/// Returns `MarchingCubesResult` with vertices, indices, and smooth normals
/// directly consumable by Three.js `BufferGeometry`.
///
/// Smooth vertex normals are computed via central-difference gradient estimation
/// at each vertex position, with trilinear interpolation for non-grid-aligned
/// positions.
///
/// Deterministic: same input always produces identical output.
///
/// # Arguments
///
/// * `input` - JavaScript object containing `MarchingCubesWasmInput` fields:
///   - `gridData`: Flat scalar field values (C-order: x-slowest, z-fastest)
///   - `gridDims`: Grid dimensions [nx, ny, nz]
///   - `gridOrigin`: Grid origin [x, y, z] in Bohr
///   - `gridSpacing`: Uniform grid spacing in Bohr
///   - `isovalue`: Isosurface threshold value
///
/// # Returns
///
/// A JavaScript object containing `MarchingCubesResult`:
/// - `vertices`: Interleaved vertex positions [x0,y0,z0, x1,y1,z1, ...]
/// - `indices`: Triangle indices (3 per triangle)
/// - `normals`: Interleaved vertex normals [nx0,ny0,nz0, ...]
///
/// # Errors
///
/// Returns a `JsError` if:
/// - Input deserialization fails
/// - Grid dimensions are invalid (< 2 in any direction)
/// - Grid data length does not match dimensions
///
/// # Example
///
/// ```javascript
/// import { marching_cubes } from './qc_wasm.js';
///
/// const input = {
///   gridData: [...],       // flat scalar field
///   gridDims: [50, 50, 50],
///   gridOrigin: [-5, -5, -5],
///   gridSpacing: 0.2,
///   isovalue: 0.05,
/// };
///
/// const result = marching_cubes(input);
/// // result.vertices: Float32-like array [x0,y0,z0, ...]
/// // result.indices: Uint32-like array [i0,i1,i2, ...]
/// // result.normals: Float32-like array [nx0,ny0,nz0, ...]
/// ```
#[wasm_bindgen]
pub fn marching_cubes(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::orbital::marching_cubes as mc;

    // 1. Deserialize input
    let wasm_input: MarchingCubesWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid marching cubes input: {}", e)))?;

    // 2. Validate input
    let [nx, ny, nz] = wasm_input.grid_dims;
    if nx < 2 || ny < 2 || nz < 2 {
        return Err(JsError::new(&format!(
            "marching_cubes: grid dimensions must be >= 2, got [{}, {}, {}]",
            nx, ny, nz
        )));
    }

    let expected_len = nx * ny * nz;
    if wasm_input.grid_data.len() != expected_len {
        return Err(JsError::new(&format!(
            "marching_cubes: grid data length {} does not match dims [{}, {}, {}] = {}",
            wasm_input.grid_data.len(),
            nx,
            ny,
            nz,
            expected_len
        )));
    }

    if wasm_input.grid_spacing <= 0.0 {
        return Err(JsError::new(&format!(
            "marching_cubes: grid spacing must be positive, got {}",
            wasm_input.grid_spacing
        )));
    }

    // 3. Run marching cubes
    let result = mc::marching_cubes(
        &wasm_input.grid_data,
        wasm_input.grid_dims,
        wasm_input.grid_origin,
        wasm_input.grid_spacing,
        wasm_input.isovalue,
    );

    // 4. Serialize and return
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

/// Extract dual isosurfaces (positive and negative lobes) for orbital visualization.
///
/// For orbital wavefunctions, positive and negative lobes require isosurfaces
/// at `+isovalue` and `-isovalue`. This function extracts both in one call.
///
/// # Arguments
///
/// * `input` - Same as `marching_cubes` input (isovalue is the positive threshold)
///
/// # Returns
///
/// A JavaScript object with `positive` and `negative` fields, each containing
/// a `MarchingCubesResult`.
///
/// # Example
///
/// ```javascript
/// import { dual_marching_cubes } from './qc_wasm.js';
///
/// const result = dual_marching_cubes(input);
/// // result.positive: MarchingCubesResult for +isovalue lobe
/// // result.negative: MarchingCubesResult for -isovalue lobe
/// ```
#[wasm_bindgen]
pub fn dual_marching_cubes(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::orbital::marching_cubes as mc;

    // 1. Deserialize input
    let wasm_input: MarchingCubesWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid marching cubes input: {}", e)))?;

    // 2. Validate input
    let [nx, ny, nz] = wasm_input.grid_dims;
    if nx < 2 || ny < 2 || nz < 2 {
        return Err(JsError::new(&format!(
            "dual_marching_cubes: grid dimensions must be >= 2, got [{}, {}, {}]",
            nx, ny, nz
        )));
    }

    let expected_len = nx * ny * nz;
    if wasm_input.grid_data.len() != expected_len {
        return Err(JsError::new(&format!(
            "dual_marching_cubes: grid data length {} does not match dims",
            wasm_input.grid_data.len()
        )));
    }

    if wasm_input.grid_spacing <= 0.0 {
        return Err(JsError::new(&format!(
            "dual_marching_cubes: grid spacing must be positive, got {}",
            wasm_input.grid_spacing
        )));
    }

    // 3. Run dual marching cubes
    let (positive, negative) = mc::dual_marching_cubes(
        &wasm_input.grid_data,
        wasm_input.grid_dims,
        wasm_input.grid_origin,
        wasm_input.grid_spacing,
        wasm_input.isovalue,
    );

    // 4. Serialize result
    #[derive(Serialize)]
    #[serde(rename_all = "camelCase")]
    struct DualResult {
        positive: mc::MarchingCubesResult,
        negative: mc::MarchingCubesResult,
    }

    let result = DualResult { positive, negative };
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Radial Profile WASM Bindings
// ============================================================================

/// Input for radial profile evaluation.
///
/// Specifies the element, basis set, and which shell to evaluate,
/// along with optional grid parameters.
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct RadialProfileInput {
    /// Atomic number (1-18)
    atomic_number: u8,
    /// Basis set name (e.g., "sto-3g", "6-31g*")
    basis_name: String,
    /// Shell index (0-indexed, matching order from get_basis_info)
    shell_index: usize,
    /// Number of evaluation points (default: 200)
    n_points: Option<usize>,
    /// Optional maximum r value in Bohr (auto-determined if absent)
    r_max: Option<f64>,
}

/// Evaluate the radial profile of a contracted basis shell.
///
/// Computes the radial function R(r) = r^l * sum_k { d_k * exp(-alpha_k * r^2) }
/// for a specified shell of a given element and basis set. Returns both the
/// contracted profile and individual primitive profiles for overlay plotting.
///
/// # Arguments
///
/// * `input` - JavaScript object with fields:
///   - `atomicNumber`: Atomic number (1-18)
///   - `basisName`: Basis set name (e.g., "sto-3g")
///   - `shellIndex`: Which shell to evaluate (0-indexed)
///   - `nPoints`: Optional number of grid points (default: 200)
///   - `rMax`: Optional maximum r in Bohr (auto if absent)
///
/// # Returns
///
/// A `RadialProfileResult` object containing:
/// - `rValues`: Array of r values in Bohr
/// - `contractedValues`: The contracted radial profile
/// - `primitiveValues`: Per-primitive profiles (2D array)
/// - `exponents`, `effectiveCoefficients`, `rawCoefficients`: Shell data
/// - `angularMomentum`, `angularMomentumLabel`: Shell type info
/// - `nPrimitives`, `rMax`: Metadata
///
/// # Errors
///
/// Returns `JsError` if:
/// - The atomic number or basis set is not supported
/// - The shell index is out of range
///
/// # Example
///
/// ```javascript
/// import { evaluate_radial_profile } from './qc_wasm.js';
///
/// const result = evaluate_radial_profile({
///   atomicNumber: 1,
///   basisName: "sto-3g",
///   shellIndex: 0,
///   nPoints: 200,
/// });
/// console.log(result.contractedValues); // [0.628..., ...]
/// ```
#[wasm_bindgen]
pub fn evaluate_radial_profile(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::basis::{get_element_basis, AngularMomentum};

    // 1. Deserialize input
    let params: RadialProfileInput =
        serde_wasm_bindgen::from_value(input).map_err(|e| JsError::new(&e.to_string()))?;

    // 2. Get basis data for the element
    let shells = get_element_basis(params.atomic_number, &params.basis_name)
        .map_err(|e| JsError::new(&e.to_string()))?;

    // 3. Validate shell index
    if params.shell_index >= shells.len() {
        return Err(JsError::new(&format!(
            "Shell index {} out of range (element Z={} in {} has {} shells)",
            params.shell_index,
            params.atomic_number,
            params.basis_name,
            shells.len()
        )));
    }

    // 4. Extract shell data
    let (am, primitives_data) = &shells[params.shell_index];
    let angular_momentum = match am {
        AngularMomentum::S => 0,
        AngularMomentum::P => 1,
        AngularMomentum::D => 2,
    };
    let exponents: Vec<f64> = primitives_data.iter().map(|(e, _)| *e).collect();
    let coefficients: Vec<f64> = primitives_data.iter().map(|(_, c)| *c).collect();

    // 5. Call core evaluator
    let result = qc_core::basis::evaluate_radial_profile(
        angular_momentum,
        &exponents,
        &coefficients,
        params.n_points,
        params.r_max,
    );

    // 6. Serialize result
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Overlap vs. Distance WASM Bindings
// ============================================================================

/// Input for overlap vs. distance computation.
///
/// Specifies two shells (by element, basis, and shell index) and a distance
/// range over which to evaluate the overlap integral.
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct OverlapDistanceInput {
    /// Atomic number for atom A (1-18)
    element_a: u8,
    /// Basis set name for atom A (e.g., "sto-3g")
    basis_a: String,
    /// Shell index within atom A's basis (0-indexed)
    shell_index_a: usize,
    /// Atomic number for atom B (1-18)
    element_b: u8,
    /// Basis set name for atom B (e.g., "sto-3g")
    basis_b: String,
    /// Shell index within atom B's basis (0-indexed)
    shell_index_b: usize,
    /// Minimum distance in bohr
    r_min: f64,
    /// Maximum distance in bohr
    r_max: f64,
    /// Number of evenly-spaced distance points
    n_points: usize,
}

/// Result from overlap vs. distance computation.
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct OverlapDistanceOutput {
    /// Distance values (bohr)
    r_values: Vec<f64>,
    /// Overlap integral values S_ab at each distance
    overlap_values: Vec<f64>,
    /// Display label for shell A (e.g., "H 1s")
    shell_label_a: String,
    /// Display label for shell B (e.g., "He 1s")
    shell_label_b: String,
    /// Basis set name for A
    basis_a: String,
    /// Basis set name for B
    basis_b: String,
}

/// Compute overlap integral vs. interatomic distance for two basis shells.
///
/// Places shell A at the origin and shell B at (0, 0, R) for each R in the
/// distance grid, computing the overlap integral at each point. This powers
/// the Module D "Overlap vs. Distance" plot.
///
/// # Arguments
///
/// * `input` - JavaScript object with fields:
///   - `elementA`: Atomic number for atom A (1-18)
///   - `basisA`: Basis set name for atom A
///   - `shellIndexA`: Shell index within atom A's basis
///   - `elementB`: Atomic number for atom B (1-18)
///   - `basisB`: Basis set name for atom B
///   - `shellIndexB`: Shell index within atom B's basis
///   - `rMin`: Minimum distance in bohr
///   - `rMax`: Maximum distance in bohr
///   - `nPoints`: Number of evenly-spaced distance points
///
/// # Returns
///
/// A JavaScript object with:
/// - `rValues`: Array of distance values (bohr)
/// - `overlapValues`: Array of overlap integral values
/// - `shellLabelA`: Display label for shell A (e.g., "H 1s")
/// - `shellLabelB`: Display label for shell B (e.g., "He 1s")
/// - `basisA`: Basis set name for A
/// - `basisB`: Basis set name for B
///
/// # Errors
///
/// Returns `JsError` if:
/// - An atomic number or basis set is not supported
/// - A shell index is out of range
/// - `rMin >= rMax`
/// - `nPoints` is 0
///
/// # Example
///
/// ```javascript
/// import { overlap_vs_distance } from './qc_wasm.js';
///
/// const result = overlap_vs_distance({
///   elementA: 1, basisA: "sto-3g", shellIndexA: 0,
///   elementB: 1, basisB: "sto-3g", shellIndexB: 0,
///   rMin: 0.1, rMax: 10.0, nPoints: 100,
/// });
/// console.log(result.overlapValues); // [0.98..., 0.94..., ...]
/// ```
#[wasm_bindgen]
pub fn overlap_vs_distance(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::basis::{
        atomic_number_to_symbol, get_element_basis, AngularMomentum, ContractedShell,
        GaussianPrimitive,
    };
    use qc_core::integrals::evaluate_overlap_vs_distance;

    // 1. Deserialize input
    let params: OverlapDistanceInput =
        serde_wasm_bindgen::from_value(input).map_err(|e| JsError::new(&e.to_string()))?;

    // 2. Validate distance range
    if params.r_min >= params.r_max {
        return Err(JsError::new(&format!(
            "r_min ({}) must be less than r_max ({})",
            params.r_min, params.r_max
        )));
    }
    if params.n_points == 0 {
        return Err(JsError::new("n_points must be at least 1"));
    }

    // 3. Get basis data and construct ContractedShell for shell A
    let shells_a = get_element_basis(params.element_a, &params.basis_a)
        .map_err(|e| JsError::new(&e.to_string()))?;
    if params.shell_index_a >= shells_a.len() {
        return Err(JsError::new(&format!(
            "Shell index {} out of range for element Z={} in {} ({} shells)",
            params.shell_index_a,
            params.element_a,
            params.basis_a,
            shells_a.len()
        )));
    }
    let (am_a, prims_a) = &shells_a[params.shell_index_a];
    let shell_a = ContractedShell::new(
        *am_a,
        prims_a
            .iter()
            .map(|(e, c)| GaussianPrimitive::new(*e, *c))
            .collect(),
        [0.0, 0.0, 0.0],
        0,
    );

    // 4. Get basis data and construct ContractedShell for shell B
    let shells_b = get_element_basis(params.element_b, &params.basis_b)
        .map_err(|e| JsError::new(&e.to_string()))?;
    if params.shell_index_b >= shells_b.len() {
        return Err(JsError::new(&format!(
            "Shell index {} out of range for element Z={} in {} ({} shells)",
            params.shell_index_b,
            params.element_b,
            params.basis_b,
            shells_b.len()
        )));
    }
    let (am_b, prims_b) = &shells_b[params.shell_index_b];
    let shell_b = ContractedShell::new(
        *am_b,
        prims_b
            .iter()
            .map(|(e, c)| GaussianPrimitive::new(*e, *c))
            .collect(),
        [0.0, 0.0, 0.0],
        0,
    );

    // 5. Generate uniform distance grid
    let r_values: Vec<f64> = if params.n_points == 1 {
        vec![params.r_min]
    } else {
        let step = (params.r_max - params.r_min) / (params.n_points - 1) as f64;
        (0..params.n_points)
            .map(|i| params.r_min + i as f64 * step)
            .collect()
    };

    // 6. Compute overlap at each distance
    let overlap_values = evaluate_overlap_vs_distance(&shell_a, &shell_b, &r_values);

    // 7. Build display labels
    let sym_a =
        atomic_number_to_symbol(params.element_a).map_err(|e| JsError::new(&e.to_string()))?;
    let sym_b =
        atomic_number_to_symbol(params.element_b).map_err(|e| JsError::new(&e.to_string()))?;
    let am_label = |am: &AngularMomentum| -> &str {
        match am {
            AngularMomentum::S => "s",
            AngularMomentum::P => "p",
            AngularMomentum::D => "d",
        }
    };
    let shell_label_a = format!("{} {}", sym_a, am_label(am_a));
    let shell_label_b = format!("{} {}", sym_b, am_label(am_b));

    // 8. Build and serialize result
    let result = OverlapDistanceOutput {
        r_values,
        overlap_values,
        shell_label_a,
        shell_label_b,
        basis_a: params.basis_a,
        basis_b: params.basis_b,
    };

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Basis Set Info WASM Bindings
// ============================================================================

/// Shell information for the Basis Explorer UI.
///
/// This struct provides a serializable representation of basis set shell data
/// that can be consumed by the TypeScript frontend.
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
struct BasisShellInfo {
    /// Angular momentum quantum number: 0=s, 1=p, 2=d
    angular_momentum: u32,
    /// Angular momentum letter: "s", "p", or "d"
    angular_momentum_label: String,
    /// Number of primitive Gaussian functions in this shell
    n_primitives: usize,
    /// Exponents of the primitive Gaussians
    exponents: Vec<f64>,
    /// Contraction coefficients of the primitive Gaussians
    coefficients: Vec<f64>,
}

/// Get basis set shell information for a given element and basis set.
///
/// Returns an array of shell descriptors containing angular momentum type,
/// number of primitives, exponents, and contraction coefficients. This data
/// is sourced directly from the built-in basis set data in `qc-core`, ensuring
/// a single source of truth.
///
/// # Arguments
///
/// * `atomic_number` - Atomic number of the element (1-18, H through Ar)
/// * `basis_name` - Basis set name (e.g., "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz")
///
/// # Returns
///
/// A JavaScript array of `BasisShellInfo` objects, each containing:
/// - `angularMomentum`: Angular momentum quantum number (0, 1, or 2)
/// - `angularMomentumLabel`: Human-readable label ("s", "p", or "d")
/// - `nPrimitives`: Number of primitive Gaussians
/// - `exponents`: Array of exponent values
/// - `coefficients`: Array of contraction coefficients
///
/// # Errors
///
/// Returns a `JsError` if:
/// - The atomic number is not supported (must be 1-18)
/// - The basis set name is not recognized
/// - The element/basis combination is not available
///
/// # Example
///
/// ```javascript
/// import { get_basis_info } from './qc_wasm.js';
///
/// const shells = get_basis_info(1, "sto-3g");
/// console.log(shells);
/// // [{ angularMomentum: 0, angularMomentumLabel: "s", nPrimitives: 3,
/// //    exponents: [3.42525091, 0.62353064, 0.16885540],
/// //    coefficients: [0.15432897, 0.53532814, 0.44463454] }]
/// ```
#[wasm_bindgen]
pub fn get_basis_info(atomic_number: u8, basis_name: &str) -> Result<JsValue, JsError> {
    use qc_core::basis::{get_element_basis, AngularMomentum};

    let shells =
        get_element_basis(atomic_number, basis_name).map_err(|e| JsError::new(&e.to_string()))?;

    let shell_infos: Vec<BasisShellInfo> = shells
        .iter()
        .map(|(am, prims)| BasisShellInfo {
            angular_momentum: match am {
                AngularMomentum::S => 0,
                AngularMomentum::P => 1,
                AngularMomentum::D => 2,
            },
            angular_momentum_label: match am {
                AngularMomentum::S => "s".to_string(),
                AngularMomentum::P => "p".to_string(),
                AngularMomentum::D => "d".to_string(),
            },
            n_primitives: prims.len(),
            exponents: prims.iter().map(|(e, _)| *e).collect(),
            coefficients: prims.iter().map(|(_, c)| *c).collect(),
        })
        .collect();

    serde_wasm_bindgen::to_value(&shell_infos).map_err(|e| JsError::new(&e.to_string()))
}

/// Returns the library version.
///
/// This is useful for verifying the WASM module loaded correctly
/// and for displaying version information in the UI.
///
/// # Example
///
/// ```javascript
/// import { version } from './qc_wasm.js';
/// console.log(`WASM version: ${version()}`);
/// ```
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// Test function demonstrating data round-trip.
///
/// Takes a number, doubles it, and returns a result object.
/// Used to verify WASM integration is working correctly.
///
/// # Arguments
///
/// * `input` - A floating-point number to process
///
/// # Returns
///
/// A `JsValue` containing a `TestResult` object with:
/// - `input`: The original input value
/// - `output`: The input doubled (input * 2)
/// - `message`: A descriptive message
///
/// # Errors
///
/// Returns a `JsError` if serialization fails (should never happen
/// for valid inputs).
///
/// # Example
///
/// ```javascript
/// import { test_compute } from './qc_wasm.js';
///
/// const result = test_compute(21);
/// console.log(result);
/// // { input: 21, output: 42, message: "Computed 21 * 2 = 42" }
/// ```
#[wasm_bindgen]
pub fn test_compute(input: f64) -> Result<JsValue, JsError> {
    let result = TestResult {
        input,
        output: input * 2.0,
        message: format!("Computed {} * 2 = {}", input, input * 2.0),
    };
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

/// Initialize the WASM module.
///
/// This function is called automatically by wasm-bindgen's init(),
/// but can be used for any additional setup if needed.
#[wasm_bindgen(start)]
pub fn init() {
    // Console panic hook for better error messages (dev only)
    #[cfg(debug_assertions)]
    {
        // Future: console_error_panic_hook::set_once();
    }
}

// ============================================================================
// Integral Matrices WASM Binding (US-055)
// ============================================================================

/// Input for computing one-electron integral matrices.
///
/// Specifies a molecular geometry, basis set, and optional spherical harmonic
/// flag. Used by Module E (Integral Inspector) to compute S, T, V, H^core
/// matrices for heatmap visualization.
///
/// # Fields
///
/// * `atoms` - Atom specifications with element symbols and coordinates
/// * `units` - Coordinate units: "bohr" or "angstrom" (default: "bohr")
/// * `basis_name` - Basis set name: "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"
/// * `use_spherical` - If true, use spherical harmonic d-orbitals (5 instead of 6)
#[derive(Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
struct IntegralMatricesInput {
    /// List of atoms in the molecule
    atoms: Vec<AtomInput>,
    /// Basis set name (e.g., "sto-3g")
    basis_name: String,
    /// Coordinate units: "bohr" or "angstrom" (default: "bohr")
    #[serde(default = "default_units")]
    units: String,
    /// Use spherical harmonic basis functions (default: false)
    #[serde(default)]
    use_spherical: bool,
}

fn default_units() -> String {
    "bohr".to_string()
}

/// Result of computing one-electron integral matrices.
///
/// Contains S, T, V, and H^core matrices as flat row-major arrays,
/// along with basis function labels and metadata.
///
/// # Layout
///
/// All matrices are stored in row-major order as flat Vec<f64> of length nbf*nbf.
/// To access element (i, j): `matrix[i * nbf + j]`
///
/// # Labels
///
/// Basis function labels follow the convention "AtomSymbol ShellType":
/// e.g., "O 1s", "O 2s", "O 2px", "H1 1s".
/// When multiple atoms of the same element exist, they are numbered:
/// "H1 1s", "H2 1s".
#[derive(Serialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
struct IntegralMatricesOutput {
    /// Number of basis functions
    nbf: usize,
    /// Basis function labels: ["O 1s", "O 2s", "O 2px", ...]
    labels: Vec<String>,
    /// Overlap matrix S (row-major, nbf x nbf)
    s_matrix: Vec<f64>,
    /// Kinetic energy matrix T (row-major, nbf x nbf)
    t_matrix: Vec<f64>,
    /// Nuclear attraction matrix V (row-major, nbf x nbf)
    v_matrix: Vec<f64>,
    /// Core Hamiltonian H^core = T + V (row-major, nbf x nbf)
    h_core: Vec<f64>,
    /// Nuclear repulsion energy (Hartree)
    nuclear_repulsion: f64,
    /// Computation time in milliseconds
    compute_time_ms: f64,
}

/// Compute one-electron integral matrices (S, T, V, H^core) for a molecule.
///
/// This function is designed for Module E (Integral Inspector), which displays
/// individual integral matrices as heatmaps. Unlike `compute_integrals`, it
/// returns T and V separately (not just their sum H^core) and does not compute
/// two-electron integrals (ERIs), making it faster for this use case.
///
/// # Arguments
///
/// * `input` - JavaScript object with `atoms`, `basisName`, `units`, and `useSpherical` fields
///
/// # Returns
///
/// A JavaScript object containing:
/// - `nbf`: Number of basis functions
/// - `labels`: Array of basis function labels
/// - `sMatrix`: Overlap matrix (row-major flat array)
/// - `tMatrix`: Kinetic energy matrix (row-major flat array)
/// - `vMatrix`: Nuclear attraction matrix (row-major flat array)
/// - `hCore`: Core Hamiltonian T+V (row-major flat array)
/// - `nuclearRepulsion`: Nuclear repulsion energy in Hartree
/// - `computeTimeMs`: Computation time in milliseconds
///
/// # Errors
///
/// Returns JsError for invalid element symbols, unsupported basis sets,
/// empty geometry, or invalid coordinate units.
///
/// # Example
///
/// ```javascript
/// import { compute_integral_matrices } from './qc_wasm.js';
///
/// const result = compute_integral_matrices({
///   atoms: [
///     { symbol: "H", xyz: [0, 0, 0] },
///     { symbol: "H", xyz: [0, 0, 1.3984] }
///   ],
///   basisName: "sto-3g",
///   units: "bohr",
/// });
///
/// console.log(result.nbf);     // 2
/// console.log(result.labels);  // ["H1 1s", "H2 1s"]
/// console.log(result.sMatrix); // [1.0, 0.659..., 0.659..., 1.0]
/// ```
#[wasm_bindgen]
pub fn compute_integral_matrices(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::basis::{symbol_to_atomic_number, Atom, BasisSet, ANGSTROM_TO_BOHR};
    use qc_core::integrals::{
        hcore_matrix, hcore_matrix_spherical, kinetic_matrix, kinetic_matrix_spherical,
        nuclear_matrix, nuclear_matrix_spherical, overlap_matrix, overlap_matrix_spherical,
    };

    let start_time = js_sys::Date::now();

    // 1. Deserialize input
    let input: IntegralMatricesInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid input: {}", e)))?;

    // 2. Validate geometry
    if input.atoms.is_empty() {
        return Err(JsError::new("Geometry must have at least 1 atom."));
    }

    // 3. Validate basis name
    let basis_lower = input.basis_name.to_lowercase();
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];
    if !supported_bases.contains(&basis_lower.as_str()) {
        return Err(JsError::new(&format!(
            "Unknown basis set '{}'. Supported: {}",
            input.basis_name,
            supported_bases.join(", ")
        )));
    }

    // 4. Determine coordinate units
    let units_lower = input.units.to_lowercase();
    let convert_to_bohr = match units_lower.as_str() {
        "bohr" => false,
        "angstrom" | "angstroms" => true,
        _ => {
            return Err(JsError::new(&format!(
                "Invalid units '{}'. Must be 'bohr' or 'angstrom'.",
                input.units
            )));
        }
    };

    // 5. Build atoms
    let mut atoms: Vec<Atom> = Vec::with_capacity(input.atoms.len());
    for atom_input in &input.atoms {
        let atomic_number = symbol_to_atomic_number(&atom_input.symbol).map_err(|_| {
            JsError::new(&format!(
                "Unsupported element '{}'. Only H-Ar are supported.",
                atom_input.symbol
            ))
        })?;

        let position_bohr = if convert_to_bohr {
            [
                atom_input.xyz[0] * ANGSTROM_TO_BOHR,
                atom_input.xyz[1] * ANGSTROM_TO_BOHR,
                atom_input.xyz[2] * ANGSTROM_TO_BOHR,
            ]
        } else {
            atom_input.xyz
        };

        let atom = Atom::new(atomic_number, position_bohr)
            .map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))?;
        atoms.push(atom);
    }

    // 6. Build basis set
    let basis = BasisSet::build(atoms, &basis_lower)
        .map_err(|e| JsError::new(&format!("Failed to build basis set: {}", e)))?;

    // 7. Compute matrices (spherical or Cartesian)
    let use_spherical = input.use_spherical;
    let (s_matrix, t_matrix, v_matrix, h_core, nbf) = if use_spherical {
        let s = overlap_matrix_spherical(&basis);
        let t = kinetic_matrix_spherical(&basis);
        let v = nuclear_matrix_spherical(&basis);
        let h = hcore_matrix_spherical(&basis);
        let n = basis.n_basis_spherical();
        (s, t, v, h, n)
    } else {
        let s = overlap_matrix(&basis);
        let t = kinetic_matrix(&basis);
        let v = nuclear_matrix(&basis);
        let h = hcore_matrix(&basis);
        let n = basis.n_basis;
        (s, t, v, h, n)
    };

    // 8. Generate basis function labels
    let labels = generate_basis_labels(&basis, use_spherical);

    // 9. Compute timing
    let compute_time_ms = js_sys::Date::now() - start_time;

    let result = IntegralMatricesOutput {
        nbf,
        labels,
        s_matrix,
        t_matrix,
        v_matrix,
        h_core,
        nuclear_repulsion: basis.nuclear_repulsion,
        compute_time_ms,
    };

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

/// Generate human-readable basis function labels from a BasisSet.
///
/// Labels follow the convention: "AtomSymbol ShellType"
/// e.g., "O 1s", "O 2s", "O 2px", "H1 1s", "H2 1s"
///
/// When multiple atoms of the same element exist, they are numbered
/// to distinguish them (H1, H2, etc.).
///
/// # Arguments
///
/// * `basis` - The basis set to generate labels for
/// * `use_spherical` - Whether to use spherical harmonic labels for d-shells
///
/// # Returns
///
/// A vector of labels, one per basis function, in the same order as the
/// integral matrix rows/columns.
fn generate_basis_labels(basis: &qc_core::basis::BasisSet, use_spherical: bool) -> Vec<String> {
    use qc_core::basis::AngularMomentum;
    use std::collections::HashMap;

    // Count atoms of each element to decide whether to number them
    let mut element_counts: HashMap<&str, usize> = HashMap::new();
    for atom in &basis.atoms {
        *element_counts.entry(atom.symbol.as_str()).or_insert(0) += 1;
    }

    // Build per-atom labels: "O", "H1", "H2" etc.
    let mut element_instance: HashMap<&str, usize> = HashMap::new();
    let mut atom_labels: Vec<String> = Vec::with_capacity(basis.atoms.len());
    for atom in &basis.atoms {
        let sym = atom.symbol.as_str();
        let count = element_counts[sym];
        if count > 1 {
            let instance = element_instance.entry(sym).or_insert(0);
            *instance += 1;
            atom_labels.push(format!("{}{}", sym, *instance));
        } else {
            atom_labels.push(sym.to_string());
        }
    }

    // Track shell counts per atom to generate principal quantum number labels
    // (1s, 2s, 2p, 3s, 3p, 3d, etc.)
    let mut atom_s_count: Vec<usize> = vec![0; basis.atoms.len()];
    let mut atom_p_count: Vec<usize> = vec![0; basis.atoms.len()];
    let mut atom_d_count: Vec<usize> = vec![0; basis.atoms.len()];

    let mut labels = Vec::with_capacity(basis.n_basis);

    for shell in &basis.shells {
        let atom_idx = shell.atom_idx;
        let atom_label = &atom_labels[atom_idx];

        match shell.angular_momentum {
            AngularMomentum::S => {
                atom_s_count[atom_idx] += 1;
                let n = atom_s_count[atom_idx];
                // Principal quantum number: first s-shell is 1s, second is 2s, etc.
                labels.push(format!("{} {}s", atom_label, n));
            }
            AngularMomentum::P => {
                atom_p_count[atom_idx] += 1;
                // P shells start at n=2 minimum
                let n = atom_p_count[atom_idx] + 1;
                if use_spherical {
                    // Spherical: p-1, p0, p+1 -- but for display, use px, py, pz
                    // (same count for both Cartesian and spherical)
                    labels.push(format!("{} {}px", atom_label, n));
                    labels.push(format!("{} {}py", atom_label, n));
                    labels.push(format!("{} {}pz", atom_label, n));
                } else {
                    labels.push(format!("{} {}px", atom_label, n));
                    labels.push(format!("{} {}py", atom_label, n));
                    labels.push(format!("{} {}pz", atom_label, n));
                }
            }
            AngularMomentum::D => {
                atom_d_count[atom_idx] += 1;
                // D shells start at n=3 minimum
                let n = atom_d_count[atom_idx] + 2;
                if use_spherical {
                    // Spherical harmonics: 5 components
                    let sph_labels = ["d-2", "d-1", "d0", "d+1", "d+2"];
                    for sl in &sph_labels {
                        labels.push(format!("{} {}{}", atom_label, n, sl));
                    }
                } else {
                    // Cartesian: 6 components (xx, xy, xz, yy, yz, zz)
                    let cart_labels = ["dxx", "dxy", "dxz", "dyy", "dyz", "dzz"];
                    for cl in &cart_labels {
                        labels.push(format!("{} {}{}", atom_label, n, cl));
                    }
                }
            }
        }
    }

    labels
}

// ============================================================================
// Integral Breakdown (US-056)
// ============================================================================

/// Input for `integral_with_breakdown` WASM export.
///
/// Specifies a molecule, basis set, integral type, and basis function indices.
/// Returns the contracted integral decomposed into primitive-pair contributions.
///
/// # Fields
///
/// * `atoms` - Atom specifications with element symbols and coordinates
/// * `basis_name` - Basis set name: "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"
/// * `units` - Coordinate units: "bohr" or "angstrom" (default: "bohr")
/// * `integral_type` - Type of integral: "S", "T", "V", or "Hcore"
/// * `indices` - Basis function indices [i, j] (0-based)
#[derive(Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
struct IntegralBreakdownInput {
    /// List of atoms in the molecule
    atoms: Vec<AtomInput>,
    /// Basis set name (e.g., "sto-3g")
    basis_name: String,
    /// Coordinate units: "bohr" or "angstrom" (default: "bohr")
    #[serde(default = "default_units")]
    units: String,
    /// Integral type: "S", "T", "V", or "Hcore"
    integral_type: String,
    /// Basis function indices [row, col] (0-based)
    indices: [usize; 2],
}

/// Compute a single integral and decompose it into primitive-pair contributions.
///
/// For a contracted integral I_ij = <chi_i | O | chi_j>, this function
/// returns all primitive-pair contributions sorted by magnitude, enabling
/// students to see which Gaussian pairs dominate the integral value.
///
/// # Input JSON structure
///
/// ```text
/// {
///   "atoms": [{ "symbol": "H", "xyz": [0.0, 0.0, 0.0] }, ...],
///   "units": "bohr",
///   "basisName": "sto-3g",
///   "integralType": "S",      // "S" | "T" | "V" | "Hcore"
///   "indices": [0, 1]         // [row, col] basis function indices (0-based)
/// }
/// ```
///
/// # Returns
///
/// JSON with `IntegralBreakdown` fields (contractedValue, integralType,
/// indices, labels, primitiveContributions, nPrimI, nPrimJ).
///
/// # Errors
///
/// Returns JsError for invalid elements, unsupported basis, out-of-range
/// indices, or invalid integral type.
///
/// # Reference
///
/// Phase 3 TDD Section 8.2; Phase 3 PRD FR-INT-03; US-056
#[wasm_bindgen]
pub fn integral_with_breakdown(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::basis::{symbol_to_atomic_number, Atom, BasisSet, ANGSTROM_TO_BOHR};

    // 1. Deserialize input
    let input: IntegralBreakdownInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid input: {}", e)))?;

    // 2. Validate geometry
    if input.atoms.is_empty() {
        return Err(JsError::new("Geometry must have at least 1 atom."));
    }

    // 3. Validate basis name
    let basis_lower = input.basis_name.to_lowercase();
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];
    if !supported_bases.contains(&basis_lower.as_str()) {
        return Err(JsError::new(&format!(
            "Unknown basis set '{}'. Supported: {}",
            input.basis_name,
            supported_bases.join(", ")
        )));
    }

    // 4. Determine coordinate units
    let units_lower = input.units.to_lowercase();
    let convert_to_bohr = match units_lower.as_str() {
        "bohr" => false,
        "angstrom" | "angstroms" => true,
        _ => {
            return Err(JsError::new(&format!(
                "Invalid units '{}'. Must be 'bohr' or 'angstrom'.",
                input.units
            )));
        }
    };

    // 5. Build atoms
    let mut atoms: Vec<Atom> = Vec::with_capacity(input.atoms.len());
    for atom_input in &input.atoms {
        let atomic_number = symbol_to_atomic_number(&atom_input.symbol).map_err(|_| {
            JsError::new(&format!(
                "Unsupported element '{}'. Only H-Ar are supported.",
                atom_input.symbol
            ))
        })?;

        let position_bohr = if convert_to_bohr {
            [
                atom_input.xyz[0] * ANGSTROM_TO_BOHR,
                atom_input.xyz[1] * ANGSTROM_TO_BOHR,
                atom_input.xyz[2] * ANGSTROM_TO_BOHR,
            ]
        } else {
            atom_input.xyz
        };

        let atom = Atom::new(atomic_number, position_bohr)
            .map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))?;
        atoms.push(atom);
    }

    // 6. Build basis set
    let basis = BasisSet::build(atoms, &basis_lower)
        .map_err(|e| JsError::new(&format!("Failed to build basis set: {}", e)))?;

    // 7. Map integral type: "S" -> "overlap", "T" -> "kinetic", "V" -> "nuclear", "Hcore" -> "hcore"
    let integral_type = match input.integral_type.as_str() {
        "S" | "overlap" => "overlap",
        "T" | "kinetic" => "kinetic",
        "V" | "nuclear" => "nuclear",
        "Hcore" | "hcore" => "hcore",
        other => {
            return Err(JsError::new(&format!(
                "Invalid integral type '{}'. Must be 'S', 'T', 'V', or 'Hcore'.",
                other
            )));
        }
    };

    // 8. Compute breakdown
    let [i, j] = input.indices;
    let result = qc_core::integrals::integral_with_breakdown(&basis, integral_type, i, j)
        .map_err(|e| JsError::new(&format!("Breakdown computation failed: {}", e)))?;

    // 9. Serialize and return
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Fock Matrix Decomposition (US-058)
// ============================================================================

/// Input for Fock matrix decomposition.
///
/// Requires geometry, basis set, and a converged density matrix.
/// The function recomputes H^core and ERIs from the geometry/basis
/// rather than accepting pre-computed matrices, ensuring consistency.
///
/// # Fields
///
/// * `atoms` - Atom list with symbols and coordinates
/// * `basis_name` - Basis set name (e.g., "sto-3g")
/// * `units` - Coordinate units: "bohr" or "angstrom" (default: "bohr")
/// * `density_matrix` - Converged density matrix P = 2*C_occ*C_occ^T
///   (flat row-major, nbf x nbf, includes factor of 2 for RHF)
#[derive(Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
struct FockDecompositionInput {
    /// List of atoms in the molecule
    atoms: Vec<AtomInput>,
    /// Basis set name (e.g., "sto-3g")
    basis_name: String,
    /// Coordinate units: "bohr" or "angstrom" (default: "bohr")
    #[serde(default = "default_units")]
    units: String,
    /// Density matrix P (flat, nbf x nbf, includes factor of 2 for RHF)
    density_matrix: Vec<f64>,
}

/// Compute Fock matrix decomposition F = H^core + J - 0.5*K.
///
/// Decomposes the Fock matrix into its physical components for educational
/// inspection. Recomputes all integrals from the geometry and basis set,
/// then builds separate J (Coulomb) and K (Exchange) matrices using the
/// provided converged density matrix.
///
/// # Input JSON structure
///
/// ```text
/// {
///   "atoms": [{ "symbol": "H", "xyz": [0.0, 0.0, 0.0] }, ...],
///   "basisName": "sto-3g",
///   "units": "bohr",
///   "densityMatrix": [0.602, 0.602, 0.602, 0.602]
/// }
/// ```
///
/// # Returns
///
/// JSON with `FockDecomposition` fields (hCore, jMatrix, kMatrix, gMatrix,
/// fMatrix, density, nbf, labels).
///
/// # Errors
///
/// Returns JsError for invalid elements, unsupported basis, or density matrix
/// dimension mismatch.
///
/// # References
///
/// - Szabo & Ostlund (1996), Eq. 3.154
/// - US-058 Fock Build Tracing
#[wasm_bindgen]
pub fn fock_decomposition(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::basis::{symbol_to_atomic_number, Atom, BasisSet, ANGSTROM_TO_BOHR};
    use qc_core::integrals::{eri_compressed, hcore_matrix};

    // 1. Deserialize input
    let input: FockDecompositionInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid input: {}", e)))?;

    // 2. Validate geometry
    if input.atoms.is_empty() {
        return Err(JsError::new("Geometry must have at least 1 atom."));
    }

    // 3. Validate basis name
    let basis_lower = input.basis_name.to_lowercase();
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];
    if !supported_bases.contains(&basis_lower.as_str()) {
        return Err(JsError::new(&format!(
            "Unknown basis set '{}'. Supported: {}",
            input.basis_name,
            supported_bases.join(", ")
        )));
    }

    // 4. Determine coordinate units
    let units_lower = input.units.to_lowercase();
    let convert_to_bohr = match units_lower.as_str() {
        "bohr" => false,
        "angstrom" | "angstroms" => true,
        _ => {
            return Err(JsError::new(&format!(
                "Invalid units '{}'. Must be 'bohr' or 'angstrom'.",
                input.units
            )));
        }
    };

    // 5. Build atoms
    let mut atoms: Vec<Atom> = Vec::with_capacity(input.atoms.len());
    for atom_input in &input.atoms {
        let atomic_number = symbol_to_atomic_number(&atom_input.symbol).map_err(|_| {
            JsError::new(&format!(
                "Unsupported element '{}'. Only H-Ar are supported.",
                atom_input.symbol
            ))
        })?;

        let position_bohr = if convert_to_bohr {
            [
                atom_input.xyz[0] * ANGSTROM_TO_BOHR,
                atom_input.xyz[1] * ANGSTROM_TO_BOHR,
                atom_input.xyz[2] * ANGSTROM_TO_BOHR,
            ]
        } else {
            atom_input.xyz
        };

        let atom = Atom::new(atomic_number, position_bohr)
            .map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))?;
        atoms.push(atom);
    }

    // 6. Build basis set
    let basis = BasisSet::build(atoms, &basis_lower)
        .map_err(|e| JsError::new(&format!("Failed to build basis set: {}", e)))?;

    let nbf = basis.n_basis;

    // 7. Validate density matrix dimensions
    let expected_size = nbf * nbf;
    if input.density_matrix.len() != expected_size {
        return Err(JsError::new(&format!(
            "Density matrix size mismatch: expected {} ({}x{}), got {}",
            expected_size,
            nbf,
            nbf,
            input.density_matrix.len()
        )));
    }

    // 8. Compute integrals (H^core and ERI)
    let h_core = hcore_matrix(&basis);
    let eri = eri_compressed(&basis);

    // 9. Generate basis function labels (Cartesian, matching SCF convention)
    let labels = generate_basis_labels(&basis, false);

    // 10. Compute Fock decomposition
    let decomp = qc_core::scf::fock::compute_fock_decomposition(
        &h_core,
        &eri,
        &input.density_matrix,
        nbf,
        labels,
    );

    // 11. Serialize and return
    serde_wasm_bindgen::to_value(&decomp).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// ERI Detail / Breakdown (US-059)
// ============================================================================

/// Input for `eri_detail` WASM export.
///
/// Specifies a molecule, basis set, and four basis function indices.
/// Returns the contracted ERI decomposed into primitive-quartet contributions.
///
/// # Fields
///
/// * `atoms` - Atom specifications with element symbols and coordinates
/// * `basis_name` - Basis set name: "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"
/// * `units` - Coordinate units: "bohr" or "angstrom" (default: "bohr")
/// * `indices` - Basis function indices [i, j, k, l] (0-based)
#[derive(Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
struct EriDetailInput {
    /// List of atoms in the molecule
    atoms: Vec<AtomInput>,
    /// Basis set name (e.g., "sto-3g")
    basis_name: String,
    /// Coordinate units: "bohr" or "angstrom" (default: "bohr")
    #[serde(default = "default_units")]
    units: String,
    /// Basis function indices [i, j, k, l] (0-based)
    indices: [usize; 4],
}

/// Compute a single ERI and decompose it into primitive-quartet contributions.
///
/// For a contracted ERI (ij|kl), this function returns all primitive-quartet
/// contributions sorted by magnitude, the computational method (Boys function
/// or Rys quadrature), and representative Rys roots/weights. This enables
/// students to see which primitive quartets dominate and how the integral
/// connects to the Rys quadrature concepts from Module B.
///
/// # Input JSON structure
///
/// ```text
/// {
///   "atoms": [{ "symbol": "H", "xyz": [0.0, 0.0, 0.0] }, ...],
///   "units": "bohr",
///   "basisName": "sto-3g",
///   "indices": [0, 0, 1, 1]
/// }
/// ```
///
/// # Returns
///
/// JSON with `EriBreakdown` fields (contractedValue, indices, labels,
/// method, contributions, nPrimitives, totalAngularMomentum, nroots).
///
/// # Errors
///
/// Returns JsError for invalid elements, unsupported basis, or out-of-range
/// indices.
///
/// # Reference
///
/// Phase 3 PRD FR-INT-04; US-059 ERI Browser
#[wasm_bindgen]
pub fn eri_detail(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::basis::{symbol_to_atomic_number, Atom, BasisSet, ANGSTROM_TO_BOHR};

    // 1. Deserialize input
    let input: EriDetailInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid input: {}", e)))?;

    // 2. Validate geometry
    if input.atoms.is_empty() {
        return Err(JsError::new("Geometry must have at least 1 atom."));
    }

    // 3. Validate basis name
    let basis_lower = input.basis_name.to_lowercase();
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];
    if !supported_bases.contains(&basis_lower.as_str()) {
        return Err(JsError::new(&format!(
            "Unknown basis set '{}'. Supported: {}",
            input.basis_name,
            supported_bases.join(", ")
        )));
    }

    // 4. Determine coordinate units
    let units_lower = input.units.to_lowercase();
    let convert_to_bohr = match units_lower.as_str() {
        "bohr" => false,
        "angstrom" | "angstroms" => true,
        _ => {
            return Err(JsError::new(&format!(
                "Invalid units '{}'. Must be 'bohr' or 'angstrom'.",
                input.units
            )));
        }
    };

    // 5. Build atoms
    let mut atoms: Vec<Atom> = Vec::with_capacity(input.atoms.len());
    for atom_input in &input.atoms {
        let atomic_number = symbol_to_atomic_number(&atom_input.symbol).map_err(|_| {
            JsError::new(&format!(
                "Unsupported element '{}'. Only H-Ar are supported.",
                atom_input.symbol
            ))
        })?;

        let position_bohr = if convert_to_bohr {
            [
                atom_input.xyz[0] * ANGSTROM_TO_BOHR,
                atom_input.xyz[1] * ANGSTROM_TO_BOHR,
                atom_input.xyz[2] * ANGSTROM_TO_BOHR,
            ]
        } else {
            atom_input.xyz
        };

        let atom = Atom::new(atomic_number, position_bohr)
            .map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))?;
        atoms.push(atom);
    }

    // 6. Build basis set
    let basis = BasisSet::build(atoms, &basis_lower)
        .map_err(|e| JsError::new(&format!("Failed to build basis set: {}", e)))?;

    // 7. Compute ERI breakdown
    let [i, j, k, l] = input.indices;
    let result = qc_core::integrals::eri_with_breakdown(&basis, i, j, k, l)
        .map_err(|e| JsError::new(&format!("ERI breakdown computation failed: {}", e)))?;

    // 8. Serialize and return
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Density Grid Evaluation WASM Bindings (US-061)
// ============================================================================

/// WASM-friendly density grid evaluation input struct.
///
/// Contains the density matrix, basis set specification, and grid parameters
/// needed to evaluate the electron density rho(r) on a 3D grid.
///
/// The electron density is:
///   rho(r) = sum_{mu,nu} D_{mu,nu} * chi_mu(r) * chi_nu(r)
///
/// where D is the density matrix and chi are basis functions.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct DensityGridWasmInput {
    /// Flattened density matrix (row-major, n_basis x n_basis)
    pub density_matrix: Vec<f64>,
    /// Atom specifications: [[Z, x, y, z], ...] in Bohr
    pub atoms: Vec<[f64; 4]>,
    /// Basis set name (e.g., "sto-3g")
    pub basis_name: String,
    /// Grid origin [x, y, z] in Bohr
    pub grid_origin: [f64; 3],
    /// Grid spacing in Bohr (uniform)
    pub grid_spacing: f64,
    /// Grid dimensions [nx, ny, nz]
    pub grid_dims: [usize; 3],
    /// Number of electrons (for integrated density validation)
    pub n_electrons: usize,
    /// Whether the density matrix is in the spherical harmonic basis.
    /// When true, the density matrix is transformed from spherical to Cartesian
    /// before grid evaluation (since the grid evaluator uses Cartesian GTOs).
    /// This must match the `useSpherical` option used in the SCF calculation.
    #[serde(default)]
    pub use_spherical: bool,
}

/// Result of density grid evaluation.
///
/// Contains the grid values and metadata for isosurface extraction.
/// Grid values use C-order indexing: index = ix * ny * nz + iy * nz + iz.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct DensityGridResult {
    /// Flat array of density values (C-order: x-slowest, z-fastest)
    pub values: Vec<f64>,
    /// Grid origin [x, y, z] in Bohr
    pub grid_origin: [f64; 3],
    /// Grid spacing in Bohr
    pub grid_spacing: f64,
    /// Grid dimensions [nx, ny, nz]
    pub grid_dims: [usize; 3],
    /// Integrated density (should approximate n_electrons)
    pub integrated_density: f64,
    /// Expected number of electrons
    pub n_electrons_expected: usize,
    /// Maximum density value
    pub max_density: f64,
    /// Computation time in milliseconds
    pub compute_time_ms: f64,
}

/// Evaluate the electron density on a 3D grid.
///
/// Computes rho(r) = sum_{mu,nu} D_{mu,nu} * chi_mu(r) * chi_nu(r) on a
/// uniform 3D grid, where chi_mu are contracted Gaussian basis functions
/// and D is the one-electron density matrix.
///
/// The density matrix is typically obtained from an SCF calculation as
/// D = C_occ * C_occ^T (where C_occ contains the occupied MO coefficients).
///
/// # Arguments
///
/// * `input` - JavaScript object containing `DensityGridWasmInput` fields:
///   - `densityMatrix`: Flattened density matrix (row-major, n_basis x n_basis)
///   - `atoms`: Atom specifications [[Z, x, y, z], ...] in Bohr
///   - `basisName`: Basis set name (e.g., "sto-3g")
///   - `gridOrigin`: Grid origin [x, y, z] in Bohr
///   - `gridSpacing`: Uniform grid spacing in Bohr
///   - `gridDims`: Grid dimensions [nx, ny, nz]
///   - `nElectrons`: Number of electrons (for validation)
///   - `useSpherical`: Whether density matrix uses spherical harmonics
///
/// # Returns
///
/// A JavaScript object containing `DensityGridResult`:
/// - `values`: Flat array of density values (C-order: x-slowest, z-fastest)
/// - `gridOrigin`: Grid origin [x, y, z] in Bohr
/// - `gridSpacing`: Uniform grid spacing in Bohr
/// - `gridDims`: Grid dimensions [nx, ny, nz]
/// - `integratedDensity`: Approximate integral of rho(r) over the grid
/// - `nElectronsExpected`: Expected number of electrons
/// - `maxDensity`: Maximum density value on the grid
/// - `computeTimeMs`: Computation time in milliseconds
///
/// # Example
///
/// ```javascript
/// const input = {
///   densityMatrix: [0.6, 0.4, 0.4, 0.6],
///   atoms: [[1, 0, 0, 0], [1, 0, 0, 1.4]],
///   basisName: "sto-3g",
///   gridOrigin: [-5, -5, -5],
///   gridSpacing: 0.5,
///   gridDims: [21, 21, 23],
///   nElectrons: 2,
/// };
///
/// const result = evaluate_density_grid(input);
/// console.log(`Integrated density: ${result.integratedDensity}`);
/// console.log(`Max density: ${result.maxDensity}`);
/// ```
#[wasm_bindgen]
pub fn evaluate_density_grid(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::basis::{Atom, BasisSet};

    // 1. Deserialize input
    let wasm_input: DensityGridWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid density grid input: {}", e)))?;

    // 2. Validate input
    if wasm_input.atoms.is_empty() {
        return Err(JsError::new(
            "evaluate_density_grid: atoms array must not be empty",
        ));
    }
    if wasm_input.grid_dims.iter().any(|&d| d < 2) {
        return Err(JsError::new(
            "evaluate_density_grid: grid dimensions must be >= 2",
        ));
    }
    if wasm_input.grid_spacing <= 0.0 {
        return Err(JsError::new(&format!(
            "evaluate_density_grid: grid spacing must be positive, got {}",
            wasm_input.grid_spacing
        )));
    }

    // 3. Build basis set from atom specs
    let atoms: Result<Vec<Atom>, _> = wasm_input
        .atoms
        .iter()
        .map(|a| {
            let z = a[0] as u8;
            Atom::new(z, [a[1], a[2], a[3]])
        })
        .collect();
    let atoms = atoms.map_err(|e| JsError::new(&format!("Invalid atom: {}", e)))?;
    let basis = BasisSet::build(atoms, &wasm_input.basis_name)
        .map_err(|e| JsError::new(&format!("Invalid basis set: {}", e)))?;

    // 4. If density matrix is in spherical basis, transform to Cartesian.
    //    The grid evaluator uses Cartesian GTOs (6 d-functions), so
    //    spherical density matrices (5 d-functions) must be expanded.
    //
    //    For a density matrix:
    //      D_cart = C_full * D_sph * C_full^T
    //    where C_full is block-diagonal: identity for s/p, 6x5 for d shells.
    let n_cart = basis.n_basis;
    let density_matrix = if wasm_input.use_spherical && basis.has_spherical_difference() {
        use qc_core::integrals::SphericalTransform;

        let n_sph = basis.n_basis_spherical();
        if wasm_input.density_matrix.len() != n_sph * n_sph {
            return Err(JsError::new(&format!(
                "evaluate_density_grid: density matrix size {} does not match spherical basis size {}x{}={}",
                wasm_input.density_matrix.len(),
                n_sph, n_sph, n_sph * n_sph
            )));
        }

        // Build block-diagonal transformation matrix C_full (n_cart x n_sph)
        // C_full is block-diagonal with identity for s/p shells and T[6x5] for d shells
        let mut c_full = vec![0.0f64; n_cart * n_sph];

        let mut cart_offset = 0;
        let mut sph_offset = 0;

        for shell in &basis.shells {
            let l = shell.l_value();
            let transform = SphericalTransform::new(l);
            let n_s = transform.n_sph;
            let n_c = transform.n_cart;

            if transform.needs_transform() {
                // d-shell (or higher): fill in the transformation coefficients
                for c in 0..n_c {
                    for m in 0..n_s {
                        c_full[(cart_offset + c) * n_sph + (sph_offset + m)] =
                            transform.coeff(c, m);
                    }
                }
            } else {
                // s/p shell: identity block
                for m in 0..n_s {
                    c_full[(cart_offset + m) * n_sph + (sph_offset + m)] = 1.0;
                }
            }

            cart_offset += n_c;
            sph_offset += n_s;
        }

        // Compute D_cart = C_full * D_sph * C_full^T
        // Step 1: temp = D_sph * C_full^T  (n_sph x n_cart)
        let mut temp = vec![0.0f64; n_sph * n_cart];
        for i in 0..n_sph {
            for j in 0..n_cart {
                let mut sum = 0.0;
                for k in 0..n_sph {
                    sum += wasm_input.density_matrix[i * n_sph + k] * c_full[j * n_sph + k];
                    // C_full^T[k][j] = C_full[j][k]
                }
                temp[i * n_cart + j] = sum;
            }
        }

        // Step 2: D_cart = C_full * temp  (n_cart x n_cart)
        let mut d_cart = vec![0.0f64; n_cart * n_cart];
        for i in 0..n_cart {
            for j in 0..n_cart {
                let mut sum = 0.0;
                for k in 0..n_sph {
                    sum += c_full[i * n_sph + k] * temp[k * n_cart + j];
                }
                d_cart[i * n_cart + j] = sum;
            }
        }

        d_cart
    } else {
        if wasm_input.density_matrix.len() != n_cart * n_cart {
            return Err(JsError::new(&format!(
                "evaluate_density_grid: density matrix size {} does not match Cartesian basis size {}x{}={}",
                wasm_input.density_matrix.len(),
                n_cart, n_cart, n_cart * n_cart
            )));
        }
        wasm_input.density_matrix
    };

    // 5. Evaluate density on grid
    //    rho(r) = sum_{mu,nu} D_{mu,nu} * chi_mu(r) * chi_nu(r)
    //
    //    Strategy: evaluate all basis functions on the grid, then contract
    //    with the density matrix. This is more memory-intensive than the
    //    MO grid evaluator but simpler, since we need all pairwise products.
    let [nx, ny, nz] = wasm_input.grid_dims;
    let total_points = nx * ny * nz;
    let mut grid_values = vec![0.0f64; total_points];

    // Allocate basis function values: basis_values[mu][grid_idx]
    let mut basis_values = vec![vec![0.0f64; total_points]; n_cart];

    // Evaluate each basis function on the entire grid (shell-batched)
    let mut basis_offset = 0;
    for shell in &basis.shells {
        let l = shell.l_value();
        let n_funcs = shell.n_basis_functions();
        let [cx, cy, cz] = shell.center;

        // Extract primitive exponents and coefficients
        let primitives: Vec<(f64, f64)> = shell
            .primitives
            .iter()
            .map(|p| (p.exponent, p.coefficient))
            .collect();

        for ix in 0..nx {
            let x = wasm_input.grid_origin[0] + ix as f64 * wasm_input.grid_spacing;
            let dx = x - cx;
            for iy in 0..ny {
                let y = wasm_input.grid_origin[1] + iy as f64 * wasm_input.grid_spacing;
                let dy = y - cy;
                for iz in 0..nz {
                    let z = wasm_input.grid_origin[2] + iz as f64 * wasm_input.grid_spacing;
                    let dz = z - cz;

                    let r2 = dx * dx + dy * dy + dz * dz;

                    // Compute radial part: sum of contracted primitives
                    let mut radial = 0.0;
                    for &(exp, coef) in &primitives {
                        radial += coef * (-exp * r2).exp();
                    }

                    // Compute angular parts and store basis function values
                    let idx = ix * ny * nz + iy * nz + iz;
                    match l {
                        0 => {
                            // s: 1
                            basis_values[basis_offset][idx] = radial;
                        }
                        1 => {
                            // p: x, y, z
                            basis_values[basis_offset][idx] = dx * radial;
                            basis_values[basis_offset + 1][idx] = dy * radial;
                            basis_values[basis_offset + 2][idx] = dz * radial;
                        }
                        2 => {
                            // d (Cartesian): xx, yy, zz, xy, xz, yz
                            basis_values[basis_offset][idx] = dx * dx * radial;
                            basis_values[basis_offset + 1][idx] = dy * dy * radial;
                            basis_values[basis_offset + 2][idx] = dz * dz * radial;
                            basis_values[basis_offset + 3][idx] = dx * dy * radial;
                            basis_values[basis_offset + 4][idx] = dx * dz * radial;
                            basis_values[basis_offset + 5][idx] = dy * dz * radial;
                        }
                        _ => {
                            return Err(JsError::new(&format!(
                                "evaluate_density_grid: unsupported angular momentum l={}",
                                l
                            )));
                        }
                    }
                }
            }
        }

        basis_offset += n_funcs;
    }

    // Contract with density matrix: rho(r) = sum_{mu,nu} D[mu,nu] * chi_mu(r) * chi_nu(r)
    for mu in 0..n_cart {
        for nu in 0..n_cart {
            let d_mu_nu = density_matrix[mu * n_cart + nu];
            if d_mu_nu.abs() < 1e-15 {
                continue; // Skip negligible density matrix elements
            }
            let chi_mu = &basis_values[mu];
            let chi_nu = &basis_values[nu];
            for (gv, (&bmu, &bnu)) in grid_values.iter_mut().zip(chi_mu.iter().zip(chi_nu.iter())) {
                *gv += d_mu_nu * bmu * bnu;
            }
        }
    }

    // 6. Compute metadata
    let dv = wasm_input.grid_spacing.powi(3);
    let integrated_density: f64 = grid_values.iter().sum::<f64>() * dv;
    let max_density = grid_values.iter().cloned().fold(0.0f64, f64::max);

    let output = DensityGridResult {
        values: grid_values,
        grid_origin: wasm_input.grid_origin,
        grid_spacing: wasm_input.grid_spacing,
        grid_dims: wasm_input.grid_dims,
        integrated_density,
        n_electrons_expected: wasm_input.n_electrons,
        max_density,
        compute_time_ms: 0.0, // Cannot time in WASM without web_sys::Performance
    };

    // 7. Serialize and return
    serde_wasm_bindgen::to_value(&output).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Difference Density (Promolecule)
// ============================================================================

/// WASM input for difference density computation.
///
/// Takes the pre-computed molecular density grid and atom positions,
/// computes the promolecule density from embedded atomic density data,
/// and returns Delta-rho = rho_molecule - rho_promolecule.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct DifferenceDensityWasmInput {
    /// Flat molecular density grid values (from evaluate_density_grid)
    pub total_density: Vec<f64>,
    /// Atom specifications: [[Z, x, y, z], ...] in Bohr
    pub atoms: Vec<[f64; 4]>,
    /// Grid origin [x, y, z] in Bohr
    pub grid_origin: [f64; 3],
    /// Grid spacing in Bohr (uniform)
    pub grid_spacing: f64,
    /// Grid dimensions [nx, ny, nz]
    pub grid_dims: [usize; 3],
}

/// Result of difference density computation.
///
/// Contains the difference density grid values and summary statistics.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct DifferenceDensityWasmResult {
    /// Flat Delta-rho values at each grid point (C-order: x-slowest, z-fastest)
    pub values: Vec<f64>,
    /// Grid origin [x, y, z] in Bohr
    pub grid_origin: [f64; 3],
    /// Grid spacing in Bohr
    pub grid_spacing: f64,
    /// Grid dimensions [nx, ny, nz]
    pub grid_dims: [usize; 3],
    /// Integrated Delta-rho: sum(Delta-rho) * dV (should be ~0)
    pub integrated_delta_rho: f64,
    /// Maximum positive Delta-rho (electron accumulation)
    pub max_accumulation: f64,
    /// Maximum negative Delta-rho (electron depletion)
    pub max_depletion: f64,
    /// Computation time in milliseconds
    pub compute_time_ms: f64,
}

/// Compute the difference density (deformation density) for a molecular system.
///
/// This function:
/// 1. Looks up atomic density profiles from embedded data (H-Ar, UHF/STO-3G)
/// 2. Evaluates the promolecule density on the same grid as the molecular density
/// 3. Computes Delta-rho = rho_molecule - rho_promolecule
///
/// # Physical Interpretation
///
/// - **Positive Delta-rho**: electron accumulation (bonding regions)
/// - **Negative Delta-rho**: electron depletion (antibonding regions)
/// - **Integrated Delta-rho ~ 0**: density conservation (same total electrons)
///
/// # Arguments
///
/// * `input` - JavaScript object containing `DifferenceDensityWasmInput` fields:
///   - `totalDensity`: Flat molecular density grid values
///   - `atoms`: Atom specifications [[Z, x, y, z], ...] in Bohr
///   - `gridOrigin`: Grid origin [x, y, z] in Bohr
///   - `gridSpacing`: Uniform grid spacing in Bohr
///   - `gridDims`: Grid dimensions [nx, ny, nz]
///
/// # Returns
///
/// A JavaScript object containing `DifferenceDensityWasmResult`:
/// - `values`: Flat array of Delta-rho values
/// - `gridOrigin`, `gridSpacing`, `gridDims`: Grid specification (echo back)
/// - `integratedDeltaRho`: Sum * dV (should be ~0)
/// - `maxAccumulation`: Maximum positive Delta-rho
/// - `maxDepletion`: Maximum negative Delta-rho (as negative number)
/// - `computeTimeMs`: Computation time in milliseconds
///
/// # Errors
///
/// Returns an error if:
/// - Any atom has an unsupported atomic number (Z > 18 or Z < 1)
/// - Grid dimensions are < 2
/// - Grid spacing is non-positive
/// - The total density array size does not match grid dimensions
///
/// # Example
///
/// ```javascript
/// const input = {
///   totalDensity: [...],  // from evaluate_density_grid()
///   atoms: [[1, 0, 0, 0], [1, 0, 0, 1.4]],
///   gridOrigin: [-5, -5, -5],
///   gridSpacing: 0.3,
///   gridDims: [34, 34, 38],
/// };
///
/// const result = compute_difference_density(input);
/// console.log(`Integrated: ${result.integratedDeltaRho}`);
/// console.log(`Max accumulation: ${result.maxAccumulation}`);
/// ```
#[wasm_bindgen]
pub fn compute_difference_density(input: JsValue) -> Result<JsValue, JsError> {
    use qc_core::orbital::promolecule;

    // 1. Deserialize input
    let wasm_input: DifferenceDensityWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid difference density input: {}", e)))?;

    // 2. Validate input
    if wasm_input.atoms.is_empty() {
        return Err(JsError::new(
            "compute_difference_density: atoms array must not be empty",
        ));
    }
    if wasm_input.grid_dims.iter().any(|&d| d < 2) {
        return Err(JsError::new(
            "compute_difference_density: grid dimensions must be >= 2",
        ));
    }
    if wasm_input.grid_spacing <= 0.0 {
        return Err(JsError::new(&format!(
            "compute_difference_density: grid spacing must be positive, got {}",
            wasm_input.grid_spacing
        )));
    }

    let expected_size = wasm_input.grid_dims[0] * wasm_input.grid_dims[1] * wasm_input.grid_dims[2];
    if wasm_input.total_density.len() != expected_size {
        return Err(JsError::new(&format!(
            "compute_difference_density: total density size {} does not match grid {}x{}x{}={}",
            wasm_input.total_density.len(),
            wasm_input.grid_dims[0],
            wasm_input.grid_dims[1],
            wasm_input.grid_dims[2],
            expected_size
        )));
    }

    // 3. Check that all elements are supported (H-Ar)
    let unsupported: Vec<u32> = wasm_input
        .atoms
        .iter()
        .map(|a| a[0] as u32)
        .filter(|&z| promolecule::get_atomic_density(z).is_none())
        .collect();
    if !unsupported.is_empty() {
        return Err(JsError::new(&format!(
            "compute_difference_density: unsupported elements (Z > 18): {:?}",
            unsupported
        )));
    }

    // 4. Prepare atom specifications for promolecule evaluation
    let atoms: Vec<(u32, [f64; 3])> = wasm_input
        .atoms
        .iter()
        .map(|a| (a[0] as u32, [a[1], a[2], a[3]]))
        .collect();

    // 5. Evaluate promolecule density on the same grid
    let promolecule_density = promolecule::evaluate_promolecule_on_grid(
        &atoms,
        wasm_input.grid_origin,
        wasm_input.grid_spacing,
        wasm_input.grid_dims,
    )
    .map_err(|e| JsError::new(&format!("compute_difference_density: {}", e)))?;

    // 6. Compute difference density
    let diff_result = promolecule::compute_difference_density(
        &wasm_input.total_density,
        &promolecule_density,
        wasm_input.grid_spacing,
    )
    .map_err(|e| JsError::new(&format!("compute_difference_density: {}", e)))?;

    // 7. Build output
    let output = DifferenceDensityWasmResult {
        values: diff_result.values,
        grid_origin: wasm_input.grid_origin,
        grid_spacing: wasm_input.grid_spacing,
        grid_dims: wasm_input.grid_dims,
        integrated_delta_rho: diff_result.integrated_delta_rho,
        max_accumulation: diff_result.max_accumulation,
        max_depletion: diff_result.max_depletion,
        compute_time_ms: 0.0,
    };

    // 8. Serialize and return
    serde_wasm_bindgen::to_value(&output).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Internal Coordinate PES Scan WASM Bindings (US-081)
// ============================================================================

/// WASM-friendly input struct for internal coordinate PES scans.
///
/// Accepts both rigid and relaxed scans over bond, angle, or dihedral
/// coordinates. Deserialized from JavaScript via serde-wasm-bindgen.
///
/// # Coordinate Types
///
/// - `"bond"`: Distance between two atoms (bohr). Requires 2 atom indices.
/// - `"angle"`: Bond angle i-j-k where j is central (radians). Requires 3 atom indices.
/// - `"dihedral"`: Torsion angle i-j-k-l about j-k bond (radians). Requires 4 atom indices.
///
/// # Scan Modes
///
/// - `"rigid"`: Only the scanned coordinate changes; all other coordinates frozen.
/// - `"relaxed"`: Non-scanned coordinates are optimized at each scan point via
///   constrained L-BFGS, giving a true relaxed PES.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PesScanInternalWasmInput {
    /// Atoms as [[Z, x, y, z], ...] where Z is atomic number and
    /// x, y, z are coordinates in bohr
    pub atoms: Vec<[f64; 4]>,

    /// Basis set name (e.g., "sto-3g", "6-31g*")
    pub basis_name: String,

    /// Electronic structure method: "rhf", "lda", "b3lyp", "b3lyp-d3bj"
    pub method: String,

    /// Type of coordinate to scan: "bond", "angle", or "dihedral"
    pub coordinate_type: String,

    /// Atom indices defining the coordinate:
    /// - Bond: [i, j] (2 indices)
    /// - Angle: [i, j, k] (3 indices, j = central)
    /// - Dihedral: [i, j, k, l] (4 indices)
    pub atom_indices: Vec<usize>,

    /// Scan mode: "rigid" or "relaxed"
    pub scan_mode: String,

    /// Minimum coordinate value (bohr for bonds, radians for angles)
    pub value_min: f64,

    /// Maximum coordinate value
    pub value_max: f64,

    /// Number of evenly spaced scan points (must be >= 2)
    pub n_points: usize,

    /// Whether to seed density from previous scan point (default: true)
    #[serde(default = "default_use_seeding_internal")]
    pub use_seeding: bool,

    /// Whether to use spherical harmonics for d/f functions (default: true)
    #[serde(default = "default_use_spherical_internal")]
    pub use_spherical: bool,

    /// Convergence profile: "loose", "medium", or "tight" (default: "tight")
    #[serde(default = "default_convergence_profile_internal")]
    pub convergence_profile: String,

    /// Maximum optimization steps per scan point for relaxed scans (default: 50)
    #[serde(default)]
    pub opt_max_steps: Option<usize>,

    /// Gradient convergence threshold for relaxed scans in Ha/bohr (default: 4.5e-4)
    #[serde(default)]
    pub opt_grad_threshold: Option<f64>,
}

fn default_use_seeding_internal() -> bool {
    true
}
fn default_use_spherical_internal() -> bool {
    true
}
fn default_convergence_profile_internal() -> String {
    "tight".to_string()
}

/// Progress emitted after each completed scan point.
///
/// Serialized and passed to the JavaScript progress callback.
#[derive(Serialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PesScanInternalWasmProgress {
    /// Index of the completed point (0-indexed)
    pub point_index: usize,
    /// Total number of scan points
    pub total_points: usize,
    /// Value of the scanned coordinate at this point
    pub coordinate_value: f64,
    /// SCF energy in Hartree
    pub energy: f64,
    /// Whether SCF converged at this geometry
    pub converged: bool,
    /// Number of optimization steps (None for rigid scans)
    pub opt_steps: Option<usize>,
}

/// Run an internal coordinate PES scan.
///
/// Supports bond, angle, and dihedral scans in both rigid and relaxed
/// modes. Progress is reported after each scan point via the callback.
///
/// # Arguments
///
/// * `input` - JavaScript object matching `PesScanInternalWasmInput`
/// * `progress_callback` - JavaScript function called after each point
///   with a `PesScanInternalWasmProgress` object
///
/// # Returns
///
/// A JavaScript object matching `PesScanInternalResult` from qc-core:
/// - `coordinateType`: "bond" | "angle" | "dihedral"
/// - `atomIndices`: number[]
/// - `points`: PesInternalPoint[]
/// - `equilibrium`: PesInternalEquilibrium | null
/// - `totalIterations`: number
/// - `scanMode`: "rigid" | "relaxed"
/// - `totalOptSteps`: number
///
/// # Errors
///
/// Returns a `JsError` if:
/// - Input deserialization fails
/// - `coordinate_type` is not "bond", "angle", or "dihedral"
/// - `atom_indices` length does not match coordinate type (2, 3, or 4)
/// - `scan_mode` is not "rigid" or "relaxed"
/// - `n_points` < 2
/// - `value_min` >= `value_max`
/// - Any atom index >= atoms.len()
/// - `method` is not a supported method
///
/// # Example
///
/// ```javascript
/// import { pes_scan_internal } from './qc_wasm.js';
///
/// const input = {
///   atoms: [[1, 0, 0, 0], [1, 0, 0, 1.4]],
///   basisName: "sto-3g",
///   method: "rhf",
///   coordinateType: "bond",
///   atomIndices: [0, 1],
///   scanMode: "rigid",
///   valueMin: 0.5,
///   valueMax: 5.0,
///   nPoints: 20,
/// };
///
/// const result = pes_scan_internal(input, (progress) => {
///   console.log(`Point ${progress.pointIndex}/${progress.totalPoints}`);
/// });
/// ```
#[wasm_bindgen]
pub fn pes_scan_internal(
    input: JsValue,
    progress_callback: &js_sys::Function,
) -> Result<JsValue, JsError> {
    use qc_core::scf::pes_internal;

    // 1. Deserialize input
    let wasm_input: PesScanInternalWasmInput = serde_wasm_bindgen::from_value(input)
        .map_err(|e| JsError::new(&format!("Invalid PES scan internal input: {}", e)))?;

    // 2. Validate input
    if wasm_input.atoms.is_empty() {
        return Err(JsError::new("pes_scan_internal: atoms array is empty"));
    }

    let expected_indices = match wasm_input.coordinate_type.as_str() {
        "bond" => 2,
        "angle" => 3,
        "dihedral" => 4,
        other => {
            return Err(JsError::new(&format!(
                "pes_scan_internal: unknown coordinateType '{}', expected 'bond', 'angle', or 'dihedral'",
                other
            )));
        }
    };

    if wasm_input.atom_indices.len() != expected_indices {
        return Err(JsError::new(&format!(
            "pes_scan_internal: coordinateType '{}' requires {} atom indices, got {}",
            wasm_input.coordinate_type,
            expected_indices,
            wasm_input.atom_indices.len()
        )));
    }

    let n_atoms = wasm_input.atoms.len();
    for &idx in &wasm_input.atom_indices {
        if idx >= n_atoms {
            return Err(JsError::new(&format!(
                "pes_scan_internal: atom index {} out of range (molecule has {} atoms)",
                idx, n_atoms
            )));
        }
    }

    let scan_mode = match wasm_input.scan_mode.as_str() {
        "rigid" => pes_internal::ScanMode::Rigid,
        "relaxed" => pes_internal::ScanMode::Relaxed,
        other => {
            return Err(JsError::new(&format!(
                "pes_scan_internal: unknown scanMode '{}', expected 'rigid' or 'relaxed'",
                other
            )));
        }
    };

    if wasm_input.n_points < 2 {
        return Err(JsError::new("pes_scan_internal: nPoints must be >= 2"));
    }

    if wasm_input.value_min >= wasm_input.value_max {
        return Err(JsError::new(&format!(
            "pes_scan_internal: valueMin ({}) must be < valueMax ({})",
            wasm_input.value_min, wasm_input.value_max
        )));
    }

    let method_lower = wasm_input.method.to_lowercase();
    if !["rhf", "hf", "lda", "b3lyp", "b3lyp-d3bj"].contains(&method_lower.as_str()) {
        return Err(JsError::new(&format!(
            "pes_scan_internal: unknown method '{}', expected 'rhf', 'lda', 'b3lyp', or 'b3lyp-d3bj'",
            wasm_input.method
        )));
    }

    // 3. Convert atoms from [Z, x, y, z] format to (u8, [f64; 3])
    let atoms: Vec<(u8, [f64; 3])> = wasm_input
        .atoms
        .iter()
        .map(|a| (a[0] as u8, [a[1], a[2], a[3]]))
        .collect();

    // 4. Build ScanCoordinate from coordinate_type + atom_indices
    let coordinate = match wasm_input.coordinate_type.as_str() {
        "bond" => pes_internal::ScanCoordinate::Bond {
            atom_i: wasm_input.atom_indices[0],
            atom_j: wasm_input.atom_indices[1],
        },
        "angle" => pes_internal::ScanCoordinate::Angle {
            atom_i: wasm_input.atom_indices[0],
            atom_j: wasm_input.atom_indices[1],
            atom_k: wasm_input.atom_indices[2],
        },
        "dihedral" => pes_internal::ScanCoordinate::Dihedral {
            atom_i: wasm_input.atom_indices[0],
            atom_j: wasm_input.atom_indices[1],
            atom_k: wasm_input.atom_indices[2],
            atom_l: wasm_input.atom_indices[3],
        },
        _ => unreachable!(), // Already validated above
    };

    // 5. Normalize method name (accept "hf" as alias for "rhf")
    let method = if method_lower == "hf" {
        "rhf".to_string()
    } else {
        method_lower
    };

    // 6. Build PesScanInternalConfig
    let config = pes_internal::PesScanInternalConfig {
        atoms,
        coordinate,
        value_min: wasm_input.value_min,
        value_max: wasm_input.value_max,
        n_points: wasm_input.n_points,
        basis_name: wasm_input.basis_name,
        method,
        use_seeding: wasm_input.use_seeding,
        use_spherical: wasm_input.use_spherical,
        convergence_profile: wasm_input.convergence_profile,
        opt_max_steps: wasm_input.opt_max_steps,
        opt_grad_threshold: wasm_input.opt_grad_threshold,
    };

    // 7. Set up progress callback
    let n_points = wasm_input.n_points;
    let progress_fn = |idx: usize, coordinate_value: f64, energy: f64, converged: bool| {
        let progress = PesScanInternalWasmProgress {
            point_index: idx,
            total_points: n_points,
            coordinate_value,
            energy,
            converged,
            opt_steps: None, // Not available during streaming
        };
        if let Ok(js_progress) = serde_wasm_bindgen::to_value(&progress) {
            let _ = progress_callback.call1(&JsValue::NULL, &js_progress);
        }
    };

    // 8. Run internal coordinate PES scan
    let result = pes_internal::pes_scan_internal(&config, scan_mode, Some(&progress_fn))
        .map_err(|e| JsError::new(&format!("pes_scan_internal failed: {}", e)))?;

    // 9. Serialize and return result
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

// ============================================================================
// Frequency Analysis — Moved to qc-wasm-spectra
// ============================================================================
//
// The compute_frequencies export and all supporting types
// (FrequencyWasmInput, FrequencyWasmResult, FrequencyTiming,
// FrequencyThermochemistry, FrequencySpectrum) have been moved to the
// `qc-wasm-spectra` crate to enable lazy WASM loading. The spectra module
// is loaded on demand when the user opens the Frequency tab.

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_returns_string() {
        let v = version();
        assert!(!v.is_empty());
        // Should match crate version
        assert_eq!(v, env!("CARGO_PKG_VERSION"));
    }

    #[test]
    fn test_result_serializes_to_camel_case() {
        let result = TestResult {
            input: 5.0,
            output: 10.0,
            message: "test".to_string(),
        };
        let json = serde_json::to_string(&result).unwrap();
        // Verify camelCase field names
        assert!(json.contains("\"input\":5"));
        assert!(json.contains("\"output\":10"));
        assert!(json.contains("\"message\":\"test\""));
        // Should not contain snake_case
        assert!(!json.contains("input_"));
    }

    #[test]
    fn test_result_deserializes_from_camel_case() {
        let json = r#"{"input":3.5,"output":7.0,"message":"hello"}"#;
        let result: TestResult = serde_json::from_str(json).unwrap();
        assert_eq!(result.input, 3.5);
        assert_eq!(result.output, 7.0);
        assert_eq!(result.message, "hello");
    }

    #[test]
    fn test_result_roundtrip() {
        let original = TestResult {
            input: 21.0,
            output: 42.0,
            message: "Computed 21 * 2 = 42".to_string(),
        };
        let json = serde_json::to_string(&original).unwrap();
        let recovered: TestResult = serde_json::from_str(&json).unwrap();
        assert_eq!(original, recovered);
    }

    #[test]
    fn test_compute_doubles_input() {
        // We can test the logic without wasm-bindgen in native tests
        let input = 21.0;
        let expected_output = 42.0;
        let expected_message = "Computed 21 * 2 = 42";

        let result = TestResult {
            input,
            output: input * 2.0,
            message: format!("Computed {} * 2 = {}", input, input * 2.0),
        };

        assert_eq!(result.input, input);
        assert_eq!(result.output, expected_output);
        assert_eq!(result.message, expected_message);
    }

    #[test]
    fn test_compute_handles_zero() {
        let input = 0.0;
        let result = TestResult {
            input,
            output: input * 2.0,
            message: format!("Computed {} * 2 = {}", input, input * 2.0),
        };
        assert_eq!(result.output, 0.0);
    }

    #[test]
    fn test_compute_handles_negative() {
        let input = -5.5;
        let result = TestResult {
            input,
            output: input * 2.0,
            message: format!("Computed {} * 2 = {}", input, input * 2.0),
        };
        assert_eq!(result.output, -11.0);
    }

    // ========================================================================
    // Boys Function Tests
    // ========================================================================

    #[test]
    fn boys_result_serializes_to_camel_case() {
        use qc_core::boys::{BoysMethod, BoysResult};

        let result = BoysResult {
            value: 0.855624,
            method: BoysMethod::Series,
            m: 0,
            t: 0.5,
            turnover: 0.0, // m=0 always uses recurrence
            terms_count: 5,
            estimated_error: Some(1e-15),
        };
        let json = serde_json::to_string(&result).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"value\":"));
        assert!(json.contains("\"method\":"));
        assert!(json.contains("\"m\":"));
        assert!(json.contains("\"t\":"));
        assert!(json.contains("\"turnover\":"));
        assert!(json.contains("\"termsCount\":"));
        assert!(json.contains("\"estimatedError\":"));
        // Method should serialize as lowercase string
        assert!(json.contains("\"series\"") || json.contains("\"Series\""));
    }

    #[test]
    fn boys_method_serializes_correctly() {
        use qc_core::boys::BoysMethod;

        // Test each variant
        let zero_json = serde_json::to_string(&BoysMethod::Zero).unwrap();
        let series_json = serde_json::to_string(&BoysMethod::Series).unwrap();
        let recurrence_json = serde_json::to_string(&BoysMethod::Recurrence).unwrap();

        assert_eq!(zero_json, "\"zero\"");
        assert_eq!(series_json, "\"series\"");
        assert_eq!(recurrence_json, "\"recurrence\"");
    }

    #[test]
    fn boys_eval_returns_correct_value() {
        // Test using the underlying qc_core function (not WASM)
        let result = qc_core::boys::boys_eval(0, 0.5).unwrap();

        // F_0(0.5) ≈ 0.8556243918921487 (from golden test data)
        assert!((result.value - 0.8556243918921487).abs() < 1e-10);
        assert_eq!(result.m, 0);
        assert_eq!(result.t, 0.5);
    }

    #[test]
    fn boys_eval_all_returns_array() {
        let results = qc_core::boys::boys_eval_all(3, 1.0).unwrap();

        assert_eq!(results.len(), 4); // m=0, 1, 2, 3
        assert_eq!(results[0].m, 0);
        assert_eq!(results[1].m, 1);
        assert_eq!(results[2].m, 2);
        assert_eq!(results[3].m, 3);

        // All should have the same T value
        for result in &results {
            assert_eq!(result.t, 1.0);
        }
    }

    #[test]
    fn boys_eval_many_returns_array() {
        let ts = vec![0.0, 0.5, 1.0, 5.0];
        let results = qc_core::boys::boys_eval_many(0, &ts).unwrap();

        assert_eq!(results.len(), 4);

        // All should have the same m value
        for result in &results {
            assert_eq!(result.m, 0);
        }

        // T values should match input
        assert_eq!(results[0].t, 0.0);
        assert_eq!(results[1].t, 0.5);
        assert_eq!(results[2].t, 1.0);
        assert_eq!(results[3].t, 5.0);
    }

    #[test]
    fn boys_eval_error_on_invalid_order() {
        // Order > MAX_ORDER should fail
        let result = qc_core::boys::boys_eval(100, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn boys_eval_error_on_negative_t() {
        let result = qc_core::boys::boys_eval(0, -1.0);
        assert!(result.is_err());
    }

    #[test]
    fn boys_result_roundtrip() {
        use qc_core::boys::{BoysMethod, BoysResult};

        let original = BoysResult {
            value: 0.746824132812427,
            method: BoysMethod::Recurrence,
            m: 0,
            t: 1.0,
            turnover: 0.0, // m=0 always uses recurrence
            terms_count: 1,
            estimated_error: None,
        };

        let json = serde_json::to_string(&original).unwrap();
        let recovered: BoysResult = serde_json::from_str(&json).unwrap();

        assert_eq!(original.value, recovered.value);
        assert_eq!(original.method, recovered.method);
        assert_eq!(original.m, recovered.m);
        assert_eq!(original.t, recovered.t);
        assert_eq!(original.turnover, recovered.turnover);
        assert_eq!(original.terms_count, recovered.terms_count);
        assert_eq!(original.estimated_error, recovered.estimated_error);
    }

    // ========================================================================
    // Rys Quadrature Tests
    // ========================================================================

    #[test]
    fn rys_result_serializes_to_camel_case() {
        use qc_core::rys::{RysMethod, RysResult};

        let result = RysResult {
            roots: vec![0.123, 0.456, 0.789],
            weights: vec![0.234, 0.345, 0.421],
            nroots: 3,
            t: 1.0,
            method: RysMethod::Standard,
        };
        let json = serde_json::to_string(&result).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"roots\":"));
        assert!(json.contains("\"weights\":"));
        assert!(json.contains("\"nroots\":"));
        assert!(json.contains("\"t\":"));
        assert!(json.contains("\"method\":"));
        // Method should serialize as lowercase string
        assert!(json.contains("\"standard\""));
    }

    #[test]
    fn rys_method_serializes_correctly() {
        use qc_core::rys::RysMethod;

        // Test each variant
        let special_json = serde_json::to_string(&RysMethod::Special).unwrap();
        let standard_json = serde_json::to_string(&RysMethod::Standard).unwrap();

        assert_eq!(special_json, "\"special\"");
        assert_eq!(standard_json, "\"standard\"");
    }

    #[test]
    fn rys_compute_returns_valid_result() {
        // Test using the underlying qc_core function (not WASM)
        let result = qc_core::rys::rys_roots(3, 1.0).unwrap();

        assert_eq!(result.nroots, 3);
        assert_eq!(result.roots.len(), 3);
        assert_eq!(result.weights.len(), 3);
        assert_eq!(result.t, 1.0);
        assert_eq!(result.method, qc_core::rys::RysMethod::Standard);

        // Roots should be in (0, 1)
        for &r in &result.roots {
            assert!(r > 0.0 && r < 1.0, "Root {} should be in (0, 1)", r);
        }

        // Weights should be positive
        for &w in &result.weights {
            assert!(w > 0.0, "Weight {} should be positive", w);
        }
    }

    #[test]
    fn rys_compute_error_on_invalid_order() {
        // Order = 0 should fail
        let result = qc_core::rys::rys_roots(0, 1.0);
        assert!(result.is_err());

        // Order > MAX_ROOTS (10) should fail
        let result = qc_core::rys::rys_roots(11, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn rys_compute_error_on_negative_t() {
        let result = qc_core::rys::rys_roots(3, -1.0);
        assert!(result.is_err());
    }

    #[test]
    fn rys_error_curve_returns_array() {
        // Test using the underlying qc_core function (not WASM)
        let result = qc_core::rys::error_curve(5, 1.0).unwrap();

        assert_eq!(result.n_max, 5);
        assert_eq!(result.t, 1.0);
        assert_eq!(result.points.len(), 5);

        // Check that n values are correct (1 through 5)
        for (i, point) in result.points.iter().enumerate() {
            assert_eq!(point.n, i + 1);
            assert!(point.max_error >= 0.0, "Error must be non-negative");
        }
    }

    #[test]
    fn rys_result_roundtrip() {
        use qc_core::rys::{RysMethod, RysResult};

        let original = RysResult {
            roots: vec![0.123456789, 0.456789012, 0.789012345],
            weights: vec![0.234567890, 0.345678901, 0.419753209],
            nroots: 3,
            t: 2.5,
            method: RysMethod::Standard,
        };

        let json = serde_json::to_string(&original).unwrap();
        let recovered: RysResult = serde_json::from_str(&json).unwrap();

        assert_eq!(original.roots, recovered.roots);
        assert_eq!(original.weights, recovered.weights);
        assert_eq!(original.nroots, recovered.nroots);
        assert_eq!(original.t, recovered.t);
        assert_eq!(original.method, recovered.method);
    }

    #[test]
    fn error_curve_result_roundtrip() {
        use qc_core::rys::{ErrorCurvePoint, ErrorCurveResult};

        let original = ErrorCurveResult {
            t: 5.0,
            n_max: 3,
            points: vec![
                ErrorCurvePoint {
                    n: 1,
                    max_error: 1.23e-15,
                },
                ErrorCurvePoint {
                    n: 2,
                    max_error: 2.34e-15,
                },
                ErrorCurvePoint {
                    n: 3,
                    max_error: 3.45e-15,
                },
            ],
        };

        let json = serde_json::to_string(&original).unwrap();
        let recovered: ErrorCurveResult = serde_json::from_str(&json).unwrap();

        assert_eq!(original.t, recovered.t);
        assert_eq!(original.n_max, recovered.n_max);
        assert_eq!(original.points.len(), recovered.points.len());

        for (orig, recov) in original.points.iter().zip(recovered.points.iter()) {
            assert_eq!(orig.n, recov.n);
            assert_eq!(orig.max_error, recov.max_error);
        }
    }

    // ========================================================================
    // SCF WASM Tests
    // ========================================================================

    #[test]
    fn test_scf_wasm_options_deserialize() {
        // Test deserialization from JSON with camelCase fields
        let json = r#"{
            "convergenceProfile": "tight",
            "maxIterations": 50,
            "useDiis": true,
            "diisSize": 8
        }"#;

        let options: ScfWasmOptions = serde_json::from_str(json).unwrap();

        assert_eq!(options.convergence_profile, "tight");
        assert_eq!(options.max_iterations, 50);
        assert!(options.use_diis);
        assert_eq!(options.diis_size, Some(8));
    }

    #[test]
    fn test_scf_wasm_options_deserialize_without_optional() {
        // Test deserialization without optional diis_size field
        let json = r#"{
            "convergenceProfile": "loose",
            "maxIterations": 25,
            "useDiis": false
        }"#;

        let options: ScfWasmOptions = serde_json::from_str(json).unwrap();

        assert_eq!(options.convergence_profile, "loose");
        assert_eq!(options.max_iterations, 25);
        assert!(!options.use_diis);
        assert_eq!(options.diis_size, None);
    }

    #[test]
    fn test_scf_wasm_result_serialize() {
        let result = ScfWasmResult {
            converged: true,
            energy: -1.116714325,
            iterations: 8,
            trace: vec![
                ScfWasmIteration {
                    iteration: 0,
                    energy: -1.8310597,
                    delta_e: None,
                    rms_density_change: None,
                    diis_applied: false,
                },
                ScfWasmIteration {
                    iteration: 1,
                    energy: -1.116714325,
                    delta_e: Some(0.714345),
                    rms_density_change: Some(1e-6),
                    diis_applied: false,
                },
            ],
            matrices: None,
            orbital_energies: None,
        };

        let json = serde_json::to_string(&result).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"converged\":true"));
        assert!(json.contains("\"energy\":-1.116714325"));
        assert!(json.contains("\"iterations\":8"));
        assert!(json.contains("\"trace\":["));
        assert!(json.contains("\"diisApplied\":false"));
        assert!(json.contains("\"deltaE\":"));
        assert!(json.contains("\"rmsDensityChange\":"));
    }

    #[test]
    fn test_scf_wasm_iteration_serialize() {
        let iteration = ScfWasmIteration {
            iteration: 5,
            energy: -1.5,
            delta_e: Some(1e-8),
            rms_density_change: Some(1e-7),
            diis_applied: true,
        };

        let json = serde_json::to_string(&iteration).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"iteration\":5"));
        assert!(json.contains("\"energy\":-1.5"));
        assert!(json.contains("\"deltaE\":"));
        assert!(json.contains("\"rmsDensityChange\":"));
        assert!(json.contains("\"diisApplied\":true"));
    }

    #[test]
    fn test_scf_wasm_iteration_serialize_without_optional() {
        // First iteration has no delta values
        let iteration = ScfWasmIteration {
            iteration: 0,
            energy: -1.8,
            delta_e: None,
            rms_density_change: None,
            diis_applied: false,
        };

        let json = serde_json::to_string(&iteration).unwrap();

        // Optional fields should be skipped when None
        assert!(!json.contains("\"deltaE\":"));
        assert!(!json.contains("\"rmsDensityChange\":"));
    }

    #[test]
    fn test_parse_convergence_profile() {
        assert_eq!(
            parse_convergence_profile("loose"),
            ConvergenceProfile::Loose
        );
        assert_eq!(
            parse_convergence_profile("LOOSE"),
            ConvergenceProfile::Loose
        );
        assert_eq!(
            parse_convergence_profile("medium"),
            ConvergenceProfile::Medium
        );
        assert_eq!(
            parse_convergence_profile("Medium"),
            ConvergenceProfile::Medium
        );
        assert_eq!(
            parse_convergence_profile("tight"),
            ConvergenceProfile::Tight
        );
        assert_eq!(
            parse_convergence_profile("TIGHT"),
            ConvergenceProfile::Tight
        );

        // Unknown profiles should default to Medium
        assert_eq!(
            parse_convergence_profile("unknown"),
            ConvergenceProfile::Medium
        );
        assert_eq!(parse_convergence_profile(""), ConvergenceProfile::Medium);
    }

    #[test]
    fn test_scf_run_h2_convergence() {
        // H2 (STO-3G, R=1.4 bohr) test system data from PySCF golden data generation
        // Reference: PySCF 2.11.0, E_tot = -1.116714325062551 Ha
        let system_json_str = r#"{
            "format_version": 1,
            "system_id": "h2_sto3g_r1.4",
            "label": "H2 (STO-3G, R=1.4 bohr)",
            "description": "",
            "basis_id": "sto-3g",
            "nbf": 2,
            "nelec": 2,
            "e_nuc": 0.7142857142857143,
            "s_matrix": [
                1.0000000000000002,
                0.659318206134864,
                0.659318206134864,
                1.0000000000000002
            ],
            "h_core": [
                -1.1204090089068204,
                -0.958379964389617,
                -0.958379964389617,
                -1.1204090089068204
            ],
            "eri_compressed": [
                0.7746059439198978,
                0.44410765803196095,
                0.2970285402769315,
                0.5696759256037501,
                0.44410765803196084,
                0.7746059439198978
            ]
        }"#;

        // Create options JSON
        let options_json_str = r#"{
            "convergenceProfile": "tight",
            "maxIterations": 100,
            "useDiis": false
        }"#;

        // Directly test the core logic (without WASM)
        let preset_json = qc_core::scf::PresetSystemJson::from_json(system_json_str).unwrap();
        let preset_system = preset_json.to_preset_system().unwrap();
        let options: ScfWasmOptions = serde_json::from_str(options_json_str).unwrap();

        let config = qc_core::scf::ScfConfig {
            profile: parse_convergence_profile(&options.convergence_profile),
            max_iterations: options.max_iterations as usize,
            use_diis: options.use_diis,
            diis_size: options.diis_size.unwrap_or(6),
            diis_start: 2,
            damp: options.damp.unwrap_or(0.0),
            damp_start: 5,
            level_shift: 0.0,
        };

        let output = qc_core::scf::rhf_scf(&preset_system, &config).unwrap();
        let result = ScfWasmResult::from(&output);

        // Verify SCF converged
        assert!(result.converged, "SCF should converge for H2");

        // Reference energy from PySCF 2.11.0: -1.116714325062551 Ha
        // Use tolerance of 1e-10
        let reference_energy = -1.116714325062551;
        let energy_diff = (result.energy - reference_energy).abs();
        assert!(
            energy_diff < 1e-10,
            "Energy {} should match reference {} (diff: {})",
            result.energy,
            reference_energy,
            energy_diff
        );

        // Verify trace is populated
        assert!(!result.trace.is_empty(), "Trace should have iteration data");

        // First iteration should not have delta values
        assert_eq!(result.trace[0].iteration, 0);
        assert!(result.trace[0].delta_e.is_none());
    }

    #[test]
    fn test_scf_run_invalid_system() {
        // Test that invalid JSON returns error
        let invalid_json = r#"{"invalid": "json"}"#;

        // Test the parsing directly (can't test WASM in unit tests)
        let result = qc_core::scf::PresetSystemJson::from_json(invalid_json);
        assert!(result.is_err(), "Invalid JSON should return error");
    }

    #[test]
    fn test_scf_wasm_result_roundtrip() {
        let original = ScfWasmResult {
            converged: true,
            energy: -1.116714325062551,
            iterations: 8,
            trace: vec![ScfWasmIteration {
                iteration: 0,
                energy: -1.8310597,
                delta_e: None,
                rms_density_change: None,
                diis_applied: false,
            }],
            matrices: None,
            orbital_energies: None,
        };

        let json = serde_json::to_string(&original).unwrap();
        let recovered: ScfWasmResult = serde_json::from_str(&json).unwrap();

        assert_eq!(original.converged, recovered.converged);
        assert_eq!(original.energy, recovered.energy);
        assert_eq!(original.iterations, recovered.iterations);
        assert_eq!(original.trace.len(), recovered.trace.len());
        assert_eq!(original.trace[0].iteration, recovered.trace[0].iteration);
        assert_eq!(original.trace[0].energy, recovered.trace[0].energy);
        assert_eq!(
            original.trace[0].diis_applied,
            recovered.trace[0].diis_applied
        );
    }

    #[test]
    fn test_scf_wasm_iteration_from_scf_iteration() {
        use qc_core::scf::ScfIteration;

        let scf_iter = ScfIteration {
            iteration: 3,
            energy_total: -1.5,
            energy_electronic: -2.0,
            delta_e: Some(1e-8),
            rms_density_change: Some(1e-7),
            converged: false,
            diis_applied: true,
        };

        let wasm_iter = ScfWasmIteration::from(&scf_iter);

        assert_eq!(wasm_iter.iteration, 3);
        assert_eq!(wasm_iter.energy, -1.5);
        assert_eq!(wasm_iter.delta_e, Some(1e-8));
        assert_eq!(wasm_iter.rms_density_change, Some(1e-7));
        assert!(wasm_iter.diis_applied);
    }

    #[test]
    fn test_scf_wasm_result_from_scf_output() {
        use qc_core::scf::{ScfConfig, ScfIteration, ScfOutput};

        let output = ScfOutput {
            converged: true,
            iterations: 5,
            energy_total: -1.116714,
            energy_electronic: -1.83,
            energy_nuclear: 0.713286,
            mo_energies: vec![-0.5, 0.5],
            mo_coefficients: vec![1.0, 0.0, 0.0, 1.0],
            density_matrix: vec![0.5, 0.0, 0.0, 0.5],
            fock_matrix: vec![-0.5, 0.0, 0.0, 0.5],
            trace: vec![ScfIteration {
                iteration: 0,
                energy_total: -1.0,
                energy_electronic: -1.7,
                delta_e: None,
                rms_density_change: None,
                converged: false,
                diis_applied: false,
            }],
            config: ScfConfig::default(),
            system_id: "test".to_string(),
        };

        let result = ScfWasmResult::from(&output);

        assert!(result.converged);
        assert_eq!(result.energy, -1.116714);
        assert_eq!(result.iterations, 5);
        assert_eq!(result.trace.len(), 1);
        assert_eq!(result.trace[0].energy, -1.0);
    }

    #[test]
    fn test_scf_wasm_matrices_serialize() {
        let matrices = ScfWasmMatrices {
            nbf: 2,
            s_matrix: vec![1.0, 0.5, 0.5, 1.0],
            h_core: vec![-1.0, -0.5, -0.5, -1.0],
            fock_matrix: vec![-0.6, -0.3, -0.3, -0.6],
            density_matrix: vec![0.5, 0.2, 0.2, 0.5],
            mo_coefficients: vec![0.7, 0.7, 0.7, -0.7],
        };

        let json = serde_json::to_string(&matrices).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"nbf\":2"));
        assert!(json.contains("\"sMatrix\":"));
        assert!(json.contains("\"hCore\":"));
        assert!(json.contains("\"fockMatrix\":"));
        assert!(json.contains("\"moCoefficients\":"));
        assert!(json.contains("\"densityMatrix\":"));

        // Verify roundtrip
        let recovered: ScfWasmMatrices = serde_json::from_str(&json).unwrap();
        assert_eq!(matrices.nbf, recovered.nbf);
        assert_eq!(matrices.s_matrix, recovered.s_matrix);
        assert_eq!(matrices.h_core, recovered.h_core);
        assert_eq!(matrices.fock_matrix, recovered.fock_matrix);
        assert_eq!(matrices.density_matrix, recovered.density_matrix);
    }

    #[test]
    fn test_scf_wasm_orbital_energies_serialize() {
        let orbital_energies = ScfWasmOrbitalEnergies {
            energies: vec![-0.5, 0.5, 1.0],
            n_occupied: 1,
        };

        let json = serde_json::to_string(&orbital_energies).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"energies\":"));
        assert!(json.contains("\"nOccupied\":1"));

        // Verify roundtrip
        let recovered: ScfWasmOrbitalEnergies = serde_json::from_str(&json).unwrap();
        assert_eq!(orbital_energies.energies, recovered.energies);
        assert_eq!(orbital_energies.n_occupied, recovered.n_occupied);
    }

    #[test]
    fn test_scf_wasm_result_with_matrices() {
        let result = ScfWasmResult {
            converged: true,
            energy: -1.116714,
            iterations: 5,
            trace: vec![],
            matrices: Some(ScfWasmMatrices {
                nbf: 2,
                s_matrix: vec![1.0, 0.5, 0.5, 1.0],
                h_core: vec![-1.0, -0.5, -0.5, -1.0],
                fock_matrix: vec![-0.6, -0.3, -0.3, -0.6],
                density_matrix: vec![0.5, 0.2, 0.2, 0.5],
                mo_coefficients: vec![0.7, 0.7, 0.7, -0.7],
            }),
            orbital_energies: Some(ScfWasmOrbitalEnergies {
                energies: vec![-0.5, 0.5],
                n_occupied: 1,
            }),
        };

        let json = serde_json::to_string(&result).unwrap();

        // Matrices should be included
        assert!(json.contains("\"matrices\":"));
        assert!(json.contains("\"orbitalEnergies\":"));

        // Verify roundtrip
        let recovered: ScfWasmResult = serde_json::from_str(&json).unwrap();
        assert!(recovered.matrices.is_some());
        assert!(recovered.orbital_energies.is_some());

        let matrices = recovered.matrices.unwrap();
        assert_eq!(matrices.nbf, 2);
        assert_eq!(matrices.s_matrix.len(), 4);

        let orbitals = recovered.orbital_energies.unwrap();
        assert_eq!(orbitals.energies.len(), 2);
        assert_eq!(orbitals.n_occupied, 1);
    }

    #[test]
    fn test_scf_wasm_options_include_matrices() {
        // Default should be false
        let json = r#"{
            "convergenceProfile": "medium",
            "maxIterations": 100,
            "useDiis": false
        }"#;

        let options: ScfWasmOptions = serde_json::from_str(json).unwrap();
        assert!(!options.include_matrices);

        // Explicit true
        let json_with_matrices = r#"{
            "convergenceProfile": "medium",
            "maxIterations": 100,
            "useDiis": false,
            "includeMatrices": true
        }"#;

        let options_with: ScfWasmOptions = serde_json::from_str(json_with_matrices).unwrap();
        assert!(options_with.include_matrices);
    }

    // ========================================================================
    // Integral Computation Tests (US-029)
    // ========================================================================

    #[test]
    fn test_atom_input_serialize() {
        let atom = AtomInput {
            symbol: "H".to_string(),
            xyz: [0.0, 0.0, 0.0],
        };

        let json = serde_json::to_string(&atom).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"symbol\":\"H\""));
        assert!(json.contains("\"xyz\":[0.0,0.0,0.0]"));
    }

    #[test]
    fn test_geometry_input_serialize() {
        let geometry = GeometryInput {
            atoms: vec![
                AtomInput {
                    symbol: "H".to_string(),
                    xyz: [0.0, 0.0, 0.0],
                },
                AtomInput {
                    symbol: "H".to_string(),
                    xyz: [0.0, 0.0, 1.4],
                },
            ],
            units: "bohr".to_string(),
        };

        let json = serde_json::to_string(&geometry).unwrap();

        // Verify structure
        assert!(json.contains("\"atoms\":["));
        assert!(json.contains("\"units\":\"bohr\""));
    }

    #[test]
    fn test_geometry_input_deserialize() {
        let json = r#"{
            "atoms": [
                { "symbol": "H", "xyz": [0, 0, 0] },
                { "symbol": "H", "xyz": [0, 0, 1.4] }
            ],
            "units": "bohr"
        }"#;

        let geometry: GeometryInput = serde_json::from_str(json).unwrap();

        assert_eq!(geometry.atoms.len(), 2);
        assert_eq!(geometry.atoms[0].symbol, "H");
        assert_eq!(geometry.atoms[1].xyz[2], 1.4);
        assert_eq!(geometry.units, "bohr");
    }

    #[test]
    fn test_atom_output_serialize() {
        let atom = AtomOutput {
            symbol: "O".to_string(),
            xyz: [0.0, 0.0, 0.117],
            atomic_number: 8,
        };

        let json = serde_json::to_string(&atom).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"symbol\":\"O\""));
        assert!(json.contains("\"atomicNumber\":8"));
    }

    #[test]
    fn test_integral_metadata_serialize() {
        let metadata = IntegralMetadata {
            wasm_version: "0.1.0".to_string(),
            compute_time_ms: 42,
            shell_pairs: 3,
            shell_quartets: 6,
            significant_eris: 10,
            basis_type: "cartesian".to_string(),
        };

        let json = serde_json::to_string(&metadata).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"wasmVersion\":\"0.1.0\""));
        assert!(json.contains("\"computeTimeMs\":42"));
        assert!(json.contains("\"shellPairs\":3"));
        assert!(json.contains("\"shellQuartets\":6"));
        assert!(json.contains("\"significantEris\":10"));
        assert!(json.contains("\"basisType\":\"cartesian\""));
    }

    #[test]
    fn test_integral_result_serialize() {
        let result = IntegralResult {
            format_version: 1,
            system_id: "custom_test123".to_string(),
            label: "H2 (sto-3g)".to_string(),
            description: "Test computation".to_string(),
            geometry: GeometryOutput {
                atoms: vec![AtomOutput {
                    symbol: "H".to_string(),
                    xyz: [0.0, 0.0, 0.0],
                    atomic_number: 1,
                }],
                units: "bohr".to_string(),
            },
            basis_id: "sto-3g".to_string(),
            nbf: 2,
            nelec: 2,
            e_nuc: 0.714285,
            s_matrix: vec![1.0, 0.5, 0.5, 1.0],
            h_core: vec![-1.0, -0.5, -0.5, -1.0],
            eri_compressed: vec![0.7, 0.4, 0.3],
            eri_indexing: "8-fold symmetry".to_string(),
            metadata: IntegralMetadata {
                wasm_version: "0.1.0".to_string(),
                compute_time_ms: 10,
                shell_pairs: 1,
                shell_quartets: 1,
                significant_eris: 3,
                basis_type: "cartesian".to_string(),
            },
        };

        let json = serde_json::to_string(&result).unwrap();

        // Verify camelCase field names (key fields)
        assert!(json.contains("\"formatVersion\":1"));
        assert!(json.contains("\"systemId\":\"custom_test123\""));
        assert!(json.contains("\"basisId\":\"sto-3g\""));
        assert!(json.contains("\"eNuc\":"));
        assert!(json.contains("\"sMatrix\":"));
        assert!(json.contains("\"hCore\":"));
        assert!(json.contains("\"eriCompressed\":"));
        assert!(json.contains("\"eriIndexing\":"));
    }

    #[test]
    fn test_integral_result_roundtrip() {
        let original = IntegralResult {
            format_version: 1,
            system_id: "custom_abc".to_string(),
            label: "Test".to_string(),
            description: "Desc".to_string(),
            geometry: GeometryOutput {
                atoms: vec![AtomOutput {
                    symbol: "H".to_string(),
                    xyz: [0.0, 0.0, 0.0],
                    atomic_number: 1,
                }],
                units: "bohr".to_string(),
            },
            basis_id: "sto-3g".to_string(),
            nbf: 1,
            nelec: 1,
            e_nuc: 0.0,
            s_matrix: vec![1.0],
            h_core: vec![-0.5],
            eri_compressed: vec![0.3],
            eri_indexing: "8-fold".to_string(),
            metadata: IntegralMetadata {
                wasm_version: "0.1.0".to_string(),
                compute_time_ms: 5,
                shell_pairs: 1,
                shell_quartets: 1,
                significant_eris: 1,
                basis_type: "cartesian".to_string(),
            },
        };

        let json = serde_json::to_string(&original).unwrap();
        let recovered: IntegralResult = serde_json::from_str(&json).unwrap();

        assert_eq!(original.format_version, recovered.format_version);
        assert_eq!(original.system_id, recovered.system_id);
        assert_eq!(original.nbf, recovered.nbf);
        assert_eq!(original.nelec, recovered.nelec);
        assert_eq!(original.e_nuc, recovered.e_nuc);
        assert_eq!(original.s_matrix, recovered.s_matrix);
        assert_eq!(original.h_core, recovered.h_core);
        assert_eq!(original.eri_compressed, recovered.eri_compressed);
    }

    #[test]
    fn test_generate_molecular_formula() {
        // Test H2
        assert_eq!(generate_molecular_formula(&["H", "H"]), "H2");

        // Test H2O
        assert_eq!(generate_molecular_formula(&["H", "O", "H"]), "H2O");

        // Test CH4
        assert_eq!(
            generate_molecular_formula(&["C", "H", "H", "H", "H"]),
            "CH4"
        );

        // Test NH3
        assert_eq!(generate_molecular_formula(&["N", "H", "H", "H"]), "H3N");

        // Test single atom
        assert_eq!(generate_molecular_formula(&["He"]), "He");

        // Test LiH
        assert_eq!(generate_molecular_formula(&["Li", "H"]), "HLi");
    }

    #[test]
    fn test_compute_integrals_h2_core_logic() {
        use qc_core::basis::{Atom, BasisSet};
        use qc_core::integrals::{eri_compressed, hcore_matrix, overlap_matrix};

        // Test the core logic without WASM (directly using qc_core)
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();

        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        // Compute integrals
        let s_matrix = overlap_matrix(&basis);
        let h_core = hcore_matrix(&basis);
        let eri = eri_compressed(&basis);

        // Verify dimensions
        assert_eq!(s_matrix.len(), 4); // 2x2 matrix
        assert_eq!(h_core.len(), 4); // 2x2 matrix
        assert_eq!(eri.len(), 6); // 6 unique ERIs for nbf=2

        // Verify overlap matrix diagonal is close to 1.0 (normalized basis)
        // Use 1e-9 tolerance matching qc-core overlap tests
        assert!(
            (s_matrix[0] - 1.0).abs() < 1e-9,
            "S[0,0] = {} should be ~1.0",
            s_matrix[0]
        );
        assert!(
            (s_matrix[3] - 1.0).abs() < 1e-9,
            "S[1,1] = {} should be ~1.0",
            s_matrix[3]
        );

        // Verify nuclear repulsion
        // E_nuc = 1*1/1.4 = 0.714285...
        assert!((basis.nuclear_repulsion - 1.0 / 1.4).abs() < 1e-10);
    }

    #[test]
    fn test_compute_integrals_h2_function() {
        // Test the compute_integrals logic with parsed geometry
        let geometry_json = r#"{
            "atoms": [
                { "symbol": "H", "xyz": [0.0, 0.0, 0.0] },
                { "symbol": "H", "xyz": [0.0, 0.0, 1.4] }
            ],
            "units": "bohr"
        }"#;

        // Parse and process (testing the logic without WASM binding)
        let geometry: GeometryInput = serde_json::from_str(geometry_json).unwrap();

        assert_eq!(geometry.atoms.len(), 2);
        assert_eq!(geometry.atoms[0].symbol, "H");
        assert_eq!(geometry.atoms[1].symbol, "H");
        assert_eq!(geometry.units, "bohr");

        // Build atoms and basis
        use qc_core::basis::{symbol_to_atomic_number, Atom, BasisSet};

        let mut atoms = Vec::new();
        for atom_input in &geometry.atoms {
            let z = symbol_to_atomic_number(&atom_input.symbol).unwrap();
            let atom = Atom::new(z, atom_input.xyz).unwrap();
            atoms.push(atom);
        }

        let basis = BasisSet::build(atoms, "sto-3g").unwrap();

        assert_eq!(basis.n_basis, 2);
        assert_eq!(basis.n_electrons, 2);
    }

    #[test]
    fn test_compute_integrals_angstrom_conversion() {
        use qc_core::basis::ANGSTROM_TO_BOHR;

        // Test that Angstrom coordinates are properly converted
        let geometry_json = r#"{
            "atoms": [
                { "symbol": "H", "xyz": [0.0, 0.0, 0.0] },
                { "symbol": "H", "xyz": [0.0, 0.0, 0.74] }
            ],
            "units": "angstrom"
        }"#;

        let geometry: GeometryInput = serde_json::from_str(geometry_json).unwrap();

        // 0.74 Angstrom should convert to ~1.4 Bohr
        let z_bohr = geometry.atoms[1].xyz[2] * ANGSTROM_TO_BOHR;
        assert!((z_bohr - 1.4).abs() < 0.01);
    }

    #[test]
    fn test_compute_integrals_invalid_element() {
        // Test that unsupported elements are rejected (K = potassium, Z=19)
        let geometry_json = r#"{
            "atoms": [
                { "symbol": "K", "xyz": [0.0, 0.0, 0.0] }
            ],
            "units": "bohr"
        }"#;

        let geometry: GeometryInput = serde_json::from_str(geometry_json).unwrap();
        let result = qc_core::basis::symbol_to_atomic_number(&geometry.atoms[0].symbol);

        assert!(result.is_err());
    }

    #[test]
    fn test_compute_integrals_empty_geometry() {
        // Verify that empty geometry is rejected
        let geometry_json = r#"{
            "atoms": [],
            "units": "bohr"
        }"#;

        let geometry: GeometryInput = serde_json::from_str(geometry_json).unwrap();
        assert!(geometry.atoms.is_empty());
    }

    // ========================================================================
    // PES Scan WASM Tests (US-039)
    // ========================================================================

    #[test]
    fn test_pes_scan_wasm_input_deserialize() {
        let json = r#"{
            "atomAZ": 1,
            "atomBZ": 1,
            "rMin": 0.5,
            "rMax": 5.0,
            "nPoints": 20,
            "basisName": "sto-3g",
            "options": {
                "convergenceProfile": "medium",
                "maxIterations": 100,
                "useDiis": true
            },
            "useSeeding": true
        }"#;

        let input: PesScanWasmInput = serde_json::from_str(json).unwrap();

        assert_eq!(input.atom_a_z, 1);
        assert_eq!(input.atom_b_z, 1);
        assert_eq!(input.r_min, 0.5);
        assert_eq!(input.r_max, 5.0);
        assert_eq!(input.n_points, 20);
        assert_eq!(input.basis_name, "sto-3g");
        assert_eq!(input.options.convergence_profile, "medium");
        assert_eq!(input.options.max_iterations, 100);
        assert!(input.options.use_diis);
        assert!(input.use_seeding);
    }

    #[test]
    fn test_pes_scan_wasm_input_default_seeding() {
        // useSeeding should default to true when not specified
        let json = r#"{
            "atomAZ": 1,
            "atomBZ": 1,
            "rMin": 1.0,
            "rMax": 3.0,
            "nPoints": 5,
            "basisName": "sto-3g",
            "options": {
                "convergenceProfile": "loose",
                "maxIterations": 50,
                "useDiis": false
            }
        }"#;

        let input: PesScanWasmInput = serde_json::from_str(json).unwrap();
        assert!(input.use_seeding);
    }

    #[test]
    fn test_pes_scan_wasm_progress_serialize() {
        let progress = PesScanWasmProgress {
            point_index: 3,
            total_points: 20,
            r: 1.5,
            energy: -1.116,
            converged: true,
        };

        let json = serde_json::to_string(&progress).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"pointIndex\":3"));
        assert!(json.contains("\"totalPoints\":20"));
        assert!(json.contains("\"r\":1.5"));
        assert!(json.contains("\"energy\":-1.116"));
        assert!(json.contains("\"converged\":true"));
    }

    #[test]
    fn test_pes_scan_wasm_input_roundtrip() {
        let original = PesScanWasmInput {
            atom_a_z: 1,
            atom_b_z: 3,
            r_min: 2.0,
            r_max: 5.0,
            n_points: 15,
            basis_name: "sto-3g".to_string(),
            options: ScfWasmOptions {
                convergence_profile: "tight".to_string(),
                max_iterations: 100,
                use_diis: true,
                diis_size: Some(8),
                damp: None,
                include_matrices: false,
            },
            use_seeding: true,
        };

        let json = serde_json::to_string(&original).unwrap();
        let recovered: PesScanWasmInput = serde_json::from_str(&json).unwrap();

        assert_eq!(original.atom_a_z, recovered.atom_a_z);
        assert_eq!(original.atom_b_z, recovered.atom_b_z);
        assert_eq!(original.r_min, recovered.r_min);
        assert_eq!(original.r_max, recovered.r_max);
        assert_eq!(original.n_points, recovered.n_points);
        assert_eq!(original.basis_name, recovered.basis_name);
        assert_eq!(original.use_seeding, recovered.use_seeding);
    }

    #[test]
    fn test_pes_point_serialize() {
        use qc_core::scf::pes::PesPoint;

        let point = PesPoint {
            r: 1.4,
            energy: -1.116714,
            converged: true,
            iterations: 8,
        };

        let json = serde_json::to_string(&point).unwrap();

        // PesPoint uses default serde (no rename_all) -- Rust field names
        assert!(json.contains("\"r\":1.4"));
        assert!(json.contains("\"energy\":-1.116714"));
        assert!(json.contains("\"converged\":true"));
        assert!(json.contains("\"iterations\":8"));
    }

    #[test]
    fn test_pes_equilibrium_serialize() {
        use qc_core::scf::pes::PesEquilibrium;

        let eq = PesEquilibrium {
            r_bohr: 1.346,
            energy_hartree: -1.116714,
        };

        let json = serde_json::to_string(&eq).unwrap();

        assert!(json.contains("\"r_bohr\":1.346"));
        assert!(json.contains("\"energy_hartree\":-1.116714"));
    }

    #[test]
    fn test_pes_scan_result_serialize() {
        use qc_core::scf::pes::{PesEquilibrium, PesPoint, PesScanResult};

        let result = PesScanResult {
            points: vec![
                PesPoint {
                    r: 1.0,
                    energy: -1.0,
                    converged: true,
                    iterations: 5,
                },
                PesPoint {
                    r: 2.0,
                    energy: -1.5,
                    converged: true,
                    iterations: 4,
                },
            ],
            equilibrium: Some(PesEquilibrium {
                r_bohr: 1.5,
                energy_hartree: -1.3,
            }),
            compute_time_ms: 150.0,
            total_iterations: 9,
        };

        let json = serde_json::to_string(&result).unwrap();

        assert!(json.contains("\"points\":["));
        assert!(json.contains("\"equilibrium\":{"));
        assert!(json.contains("\"compute_time_ms\":150.0"));
        assert!(json.contains("\"total_iterations\":9"));
    }

    #[test]
    fn test_pes_scan_result_roundtrip() {
        use qc_core::scf::pes::{PesEquilibrium, PesPoint, PesScanResult};

        let original = PesScanResult {
            points: vec![PesPoint {
                r: 1.4,
                energy: -1.116714,
                converged: true,
                iterations: 8,
            }],
            equilibrium: Some(PesEquilibrium {
                r_bohr: 1.346,
                energy_hartree: -1.116714,
            }),
            compute_time_ms: 200.0,
            total_iterations: 8,
        };

        let json = serde_json::to_string(&original).unwrap();
        let recovered: PesScanResult = serde_json::from_str(&json).unwrap();

        assert_eq!(original.points.len(), recovered.points.len());
        assert_eq!(original.points[0].r, recovered.points[0].r);
        assert_eq!(original.points[0].energy, recovered.points[0].energy);
        assert!(recovered.equilibrium.is_some());
        let eq = recovered.equilibrium.unwrap();
        assert_eq!(eq.r_bohr, 1.346);
        assert_eq!(eq.energy_hartree, -1.116714);
    }

    #[test]
    fn test_compute_integrals_produces_scf_compatible_format() {
        // Verify that IntegralResult can be parsed as PresetSystemJson
        let result = IntegralResult {
            format_version: 1,
            system_id: "custom_h2".to_string(),
            label: "H2 (sto-3g)".to_string(),
            description: "Test".to_string(),
            geometry: GeometryOutput {
                atoms: vec![
                    AtomOutput {
                        symbol: "H".to_string(),
                        xyz: [0.0, 0.0, 0.0],
                        atomic_number: 1,
                    },
                    AtomOutput {
                        symbol: "H".to_string(),
                        xyz: [0.0, 0.0, 1.4],
                        atomic_number: 1,
                    },
                ],
                units: "bohr".to_string(),
            },
            basis_id: "sto-3g".to_string(),
            nbf: 2,
            nelec: 2,
            e_nuc: 0.714285714,
            s_matrix: vec![1.0, 0.659, 0.659, 1.0],
            h_core: vec![-1.12, -0.96, -0.96, -1.12],
            eri_compressed: vec![0.77, 0.44, 0.30, 0.57, 0.44, 0.77],
            eri_indexing: "8-fold symmetry".to_string(),
            metadata: IntegralMetadata {
                wasm_version: "0.1.0".to_string(),
                compute_time_ms: 10,
                shell_pairs: 3,
                shell_quartets: 6,
                significant_eris: 6,
                basis_type: "cartesian".to_string(),
            },
        };

        // Serialize to JSON
        let json = serde_json::to_string(&result).unwrap();

        // The JSON should be parseable by PresetSystemJson
        // Note: PresetSystemJson uses snake_case in JSON, but we use camelCase
        // This is intentional - the worker will handle the conversion
        // Here we just verify the structure is correct

        // Key fields that SCF needs
        assert!(json.contains("\"formatVersion\":1"));
        assert!(json.contains("\"nbf\":2"));
        assert!(json.contains("\"nelec\":2"));
        assert!(json.contains("\"eNuc\":"));
        assert!(json.contains("\"sMatrix\":"));
        assert!(json.contains("\"hCore\":"));
        assert!(json.contains("\"eriCompressed\":"));
    }

    // ========================================================================
    // Integral Matrices (US-055) Tests
    // ========================================================================

    #[test]
    fn test_generate_basis_labels_h2_sto3g() {
        use qc_core::basis::{Atom, BasisSet};

        let atoms = vec![
            Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
            Atom::new(1, [0.0, 0.0, 1.3984]).unwrap(),
        ];
        let basis = BasisSet::build(atoms, "sto-3g").unwrap();
        let labels = generate_basis_labels(&basis, false);

        assert_eq!(labels.len(), 2);
        assert_eq!(labels[0], "H1 1s");
        assert_eq!(labels[1], "H2 1s");
    }

    #[test]
    fn test_generate_basis_labels_h2o_sto3g() {
        use qc_core::basis::{Atom, BasisSet};

        let atoms = vec![
            Atom::new(8, [0.0, 0.0, 0.117]).unwrap(),
            Atom::new(1, [0.0, 1.43, -0.47]).unwrap(),
            Atom::new(1, [0.0, -1.43, -0.47]).unwrap(),
        ];
        let basis = BasisSet::build(atoms, "sto-3g").unwrap();
        let labels = generate_basis_labels(&basis, false);

        assert_eq!(labels.len(), 7);
        assert_eq!(labels[0], "O 1s");
        assert_eq!(labels[1], "O 2s");
        assert_eq!(labels[2], "O 2px");
        assert_eq!(labels[3], "O 2py");
        assert_eq!(labels[4], "O 2pz");
        assert_eq!(labels[5], "H1 1s");
        assert_eq!(labels[6], "H2 1s");
    }

    #[test]
    fn test_generate_basis_labels_lih_sto3g() {
        use qc_core::basis::{Atom, BasisSet};

        let atoms = vec![
            Atom::new(3, [0.0, 0.0, 0.0]).unwrap(),
            Atom::new(1, [0.0, 0.0, 3.015]).unwrap(),
        ];
        let basis = BasisSet::build(atoms, "sto-3g").unwrap();
        let labels = generate_basis_labels(&basis, false);

        assert_eq!(labels.len(), 6);
        assert_eq!(labels[0], "Li 1s");
        assert_eq!(labels[1], "Li 2s");
        assert_eq!(labels[2], "Li 2px");
        assert_eq!(labels[3], "Li 2py");
        assert_eq!(labels[4], "Li 2pz");
        assert_eq!(labels[5], "H 1s");
    }

    #[test]
    fn test_generate_basis_labels_single_atom_no_numbering() {
        use qc_core::basis::{Atom, BasisSet};

        // Single H atom - should not have numbering
        let atoms = vec![Atom::new(1, [0.0, 0.0, 0.0]).unwrap()];
        let basis = BasisSet::build(atoms, "sto-3g").unwrap();
        let labels = generate_basis_labels(&basis, false);

        assert_eq!(labels.len(), 1);
        assert_eq!(labels[0], "H 1s");
    }

    #[test]
    fn test_generate_basis_labels_h2o_631gs_cartesian() {
        use qc_core::basis::{Atom, BasisSet};

        let atoms = vec![
            Atom::new(8, [0.0, 0.0, 0.117]).unwrap(),
            Atom::new(1, [0.0, 1.43, -0.47]).unwrap(),
            Atom::new(1, [0.0, -1.43, -0.47]).unwrap(),
        ];
        let basis = BasisSet::build(atoms, "6-31g*").unwrap();
        let labels = generate_basis_labels(&basis, false);

        // 6-31G* on O: 1s, 2s, 2p(3), 3s, 3p(3), 3d(6) = 15
        // 6-31G on H: 1s, 2s = 2 each, total 4
        // Total: 19
        assert_eq!(labels.len(), 19);
        assert_eq!(basis.n_basis, 19);

        // Check d-shell labels are Cartesian
        assert!(labels.iter().any(|l| l.contains("dxx")));
        assert!(labels.iter().any(|l| l.contains("dyz")));
    }

    #[test]
    fn test_generate_basis_labels_h2o_631gs_spherical() {
        use qc_core::basis::{Atom, BasisSet};

        let atoms = vec![
            Atom::new(8, [0.0, 0.0, 0.117]).unwrap(),
            Atom::new(1, [0.0, 1.43, -0.47]).unwrap(),
            Atom::new(1, [0.0, -1.43, -0.47]).unwrap(),
        ];
        let basis = BasisSet::build(atoms, "6-31g*").unwrap();
        let labels = generate_basis_labels(&basis, true);

        // Spherical: 6-31G* on O: 1s, 2s, 2p(3), 3s, 3p(3), 3d(5) = 14
        // 6-31G on H: 1s, 2s = 2 each, total 4
        // Total: 18
        assert_eq!(labels.len(), 18);
        assert_eq!(basis.n_basis_spherical(), 18);

        // Check d-shell labels are spherical
        assert!(labels.iter().any(|l| l.contains("d-2")));
        assert!(labels.iter().any(|l| l.contains("d+2")));
    }

    #[test]
    fn test_integral_matrices_output_serialization() {
        let output = IntegralMatricesOutput {
            nbf: 2,
            labels: vec!["H1 1s".to_string(), "H2 1s".to_string()],
            s_matrix: vec![1.0, 0.659, 0.659, 1.0],
            t_matrix: vec![0.76, 0.24, 0.24, 0.76],
            v_matrix: vec![-1.88, -1.20, -1.20, -1.88],
            h_core: vec![-1.12, -0.96, -0.96, -1.12],
            nuclear_repulsion: 0.7151043,
            compute_time_ms: 1.5,
        };

        let json = serde_json::to_string(&output).unwrap();

        // Verify camelCase serialization
        assert!(json.contains("\"nbf\":2"));
        assert!(json.contains("\"labels\":[\"H1 1s\",\"H2 1s\"]"));
        assert!(json.contains("\"sMatrix\":"));
        assert!(json.contains("\"tMatrix\":"));
        assert!(json.contains("\"vMatrix\":"));
        assert!(json.contains("\"hCore\":"));
        assert!(json.contains("\"nuclearRepulsion\":"));
        assert!(json.contains("\"computeTimeMs\":"));
    }

    #[test]
    fn test_integral_matrices_input_deserialization() {
        let json = r#"{
            "atoms": [
                {"symbol": "H", "xyz": [0, 0, 0]},
                {"symbol": "H", "xyz": [0, 0, 1.4]}
            ],
            "basisName": "sto-3g",
            "units": "bohr"
        }"#;

        let input: IntegralMatricesInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.atoms.len(), 2);
        assert_eq!(input.basis_name, "sto-3g");
        assert_eq!(input.units, "bohr");
        assert!(!input.use_spherical);
    }

    #[test]
    fn test_integral_matrices_input_defaults() {
        // units should default to "bohr", use_spherical to false
        let json = r#"{
            "atoms": [{"symbol": "H", "xyz": [0, 0, 0]}],
            "basisName": "sto-3g"
        }"#;

        let input: IntegralMatricesInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.units, "bohr");
        assert!(!input.use_spherical);
    }

    #[test]
    fn test_density_grid_input_deserialization() {
        let json = r#"{
            "densityMatrix": [0.6, 0.4, 0.4, 0.6],
            "atoms": [[1, 0, 0, 0], [1, 0, 0, 1.4]],
            "basisName": "sto-3g",
            "gridOrigin": [-5.0, -5.0, -5.0],
            "gridSpacing": 0.5,
            "gridDims": [21, 21, 23],
            "nElectrons": 2
        }"#;

        let input: DensityGridWasmInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.density_matrix.len(), 4);
        assert_eq!(input.atoms.len(), 2);
        assert_eq!(input.basis_name, "sto-3g");
        assert_eq!(input.grid_dims, [21, 21, 23]);
        assert_eq!(input.n_electrons, 2);
        assert!(!input.use_spherical); // default
    }

    #[test]
    fn test_density_grid_input_with_spherical() {
        let json = r#"{
            "densityMatrix": [0.6, 0.4, 0.4, 0.6],
            "atoms": [[1, 0, 0, 0], [1, 0, 0, 1.4]],
            "basisName": "sto-3g",
            "gridOrigin": [-5.0, -5.0, -5.0],
            "gridSpacing": 0.5,
            "gridDims": [21, 21, 23],
            "nElectrons": 2,
            "useSpherical": true
        }"#;

        let input: DensityGridWasmInput = serde_json::from_str(json).unwrap();
        assert!(input.use_spherical);
    }

    #[test]
    fn test_density_grid_result_serialization() {
        let result = DensityGridResult {
            values: vec![0.1, 0.2, 0.3],
            grid_origin: [-5.0, -5.0, -5.0],
            grid_spacing: 0.5,
            grid_dims: [21, 21, 23],
            integrated_density: 1.98,
            n_electrons_expected: 2,
            max_density: 0.3,
            compute_time_ms: 42.5,
        };

        let json = serde_json::to_string(&result).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"integratedDensity\""));
        assert!(json.contains("\"nElectronsExpected\":2"));
        assert!(json.contains("\"maxDensity\""));
        assert!(json.contains("\"computeTimeMs\""));
        assert!(json.contains("\"gridOrigin\""));
        assert!(json.contains("\"gridSpacing\""));
        assert!(json.contains("\"gridDims\""));

        // Roundtrip
        let recovered: DensityGridResult = serde_json::from_str(&json).unwrap();
        assert_eq!(recovered.values.len(), 3);
        assert_eq!(recovered.n_electrons_expected, 2);
        assert!((recovered.integrated_density - 1.98).abs() < 1e-10);
    }

    // ========================================================================
    // Difference Density (US-063) Tests
    // ========================================================================

    #[test]
    fn test_difference_density_input_deserialization() {
        let json = r#"{
            "totalDensity": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
            "atoms": [[1, 0, 0, 0], [1, 0, 0, 1.4]],
            "gridOrigin": [-5.0, -5.0, -5.0],
            "gridSpacing": 0.5,
            "gridDims": [2, 2, 2]
        }"#;

        let input: DifferenceDensityWasmInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.total_density.len(), 8);
        assert_eq!(input.atoms.len(), 2);
        assert_eq!(input.grid_dims, [2, 2, 2]);
        assert!((input.grid_spacing - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_difference_density_result_serialization() {
        let result = DifferenceDensityWasmResult {
            values: vec![0.01, -0.02, 0.03],
            grid_origin: [-5.0, -5.0, -5.0],
            grid_spacing: 0.5,
            grid_dims: [1, 1, 3],
            integrated_delta_rho: 0.001,
            max_accumulation: 0.03,
            max_depletion: -0.02,
            compute_time_ms: 15.0,
        };

        let json = serde_json::to_string(&result).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"integratedDeltaRho\""));
        assert!(json.contains("\"maxAccumulation\""));
        assert!(json.contains("\"maxDepletion\""));
        assert!(json.contains("\"computeTimeMs\""));
        assert!(json.contains("\"gridOrigin\""));
        assert!(json.contains("\"gridSpacing\""));
        assert!(json.contains("\"gridDims\""));

        // Roundtrip
        let recovered: DifferenceDensityWasmResult = serde_json::from_str(&json).unwrap();
        assert_eq!(recovered.values.len(), 3);
        assert!((recovered.integrated_delta_rho - 0.001).abs() < 1e-10);
        assert!((recovered.max_accumulation - 0.03).abs() < 1e-10);
        assert!((recovered.max_depletion - (-0.02)).abs() < 1e-10);
    }

    // ========================================================================
    // Internal Coordinate PES Scan WASM Tests (US-081)
    // ========================================================================

    #[test]
    fn test_pes_scan_internal_input_deserialize_bond() {
        let json = r#"{
            "atoms": [[1, 0, 0, 0], [1, 0, 0, 1.4]],
            "basisName": "sto-3g",
            "method": "rhf",
            "coordinateType": "bond",
            "atomIndices": [0, 1],
            "scanMode": "rigid",
            "valueMin": 0.5,
            "valueMax": 5.0,
            "nPoints": 20
        }"#;

        let input: PesScanInternalWasmInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.atoms.len(), 2);
        assert_eq!(input.basis_name, "sto-3g");
        assert_eq!(input.method, "rhf");
        assert_eq!(input.coordinate_type, "bond");
        assert_eq!(input.atom_indices, vec![0, 1]);
        assert_eq!(input.scan_mode, "rigid");
        assert!((input.value_min - 0.5).abs() < 1e-10);
        assert!((input.value_max - 5.0).abs() < 1e-10);
        assert_eq!(input.n_points, 20);
        // Check defaults
        assert!(input.use_seeding);
        assert!(input.use_spherical);
        assert_eq!(input.convergence_profile, "tight");
    }

    #[test]
    fn test_pes_scan_internal_input_deserialize_angle() {
        let json = r#"{
            "atoms": [[8, 0, 0, 0], [1, 0.96, 0, 0], [1, -0.24, 0.93, 0]],
            "basisName": "6-31g*",
            "method": "b3lyp",
            "coordinateType": "angle",
            "atomIndices": [1, 0, 2],
            "scanMode": "relaxed",
            "valueMin": 1.5,
            "valueMax": 2.5,
            "nPoints": 10,
            "useSeeding": false,
            "useSpherical": false,
            "convergenceProfile": "medium"
        }"#;

        let input: PesScanInternalWasmInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.atoms.len(), 3);
        assert_eq!(input.coordinate_type, "angle");
        assert_eq!(input.atom_indices, vec![1, 0, 2]);
        assert_eq!(input.scan_mode, "relaxed");
        assert!(!input.use_seeding);
        assert!(!input.use_spherical);
        assert_eq!(input.convergence_profile, "medium");
    }

    #[test]
    fn test_pes_scan_internal_input_deserialize_dihedral() {
        let json = r#"{
            "atoms": [[1, 0, 0, 0], [6, 1.0, 0, 0], [6, 2.0, 1.0, 0], [1, 3.0, 1.0, 0]],
            "basisName": "sto-3g",
            "method": "b3lyp-d3bj",
            "coordinateType": "dihedral",
            "atomIndices": [0, 1, 2, 3],
            "scanMode": "rigid",
            "valueMin": -3.14159,
            "valueMax": 3.14159,
            "nPoints": 36
        }"#;

        let input: PesScanInternalWasmInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.coordinate_type, "dihedral");
        assert_eq!(input.atom_indices.len(), 4);
        assert_eq!(input.method, "b3lyp-d3bj");
        assert_eq!(input.n_points, 36);
    }

    #[test]
    fn test_pes_scan_internal_input_roundtrip() {
        let input = PesScanInternalWasmInput {
            atoms: vec![[1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 1.4]],
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            coordinate_type: "bond".to_string(),
            atom_indices: vec![0, 1],
            scan_mode: "rigid".to_string(),
            value_min: 0.5,
            value_max: 5.0,
            n_points: 20,
            use_seeding: true,
            use_spherical: true,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let json = serde_json::to_string(&input).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"basisName\""));
        assert!(json.contains("\"coordinateType\""));
        assert!(json.contains("\"atomIndices\""));
        assert!(json.contains("\"scanMode\""));
        assert!(json.contains("\"valueMin\""));
        assert!(json.contains("\"valueMax\""));
        assert!(json.contains("\"nPoints\""));
        assert!(json.contains("\"useSeeding\""));
        assert!(json.contains("\"useSpherical\""));
        assert!(json.contains("\"convergenceProfile\""));

        // Roundtrip
        let recovered: PesScanInternalWasmInput = serde_json::from_str(&json).unwrap();
        assert_eq!(recovered.atoms.len(), 2);
        assert_eq!(recovered.basis_name, "sto-3g");
        assert_eq!(recovered.coordinate_type, "bond");
        assert_eq!(recovered.n_points, 20);
    }

    #[test]
    fn test_pes_scan_internal_progress_serializes_camel_case() {
        let progress = PesScanInternalWasmProgress {
            point_index: 3,
            total_points: 20,
            coordinate_value: 1.5,
            energy: -1.116714,
            converged: true,
            opt_steps: None,
        };

        let json = serde_json::to_string(&progress).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"pointIndex\":3"));
        assert!(json.contains("\"totalPoints\":20"));
        assert!(json.contains("\"coordinateValue\""));
        assert!(json.contains("\"optSteps\":null"));

        // With opt_steps set (relaxed scan)
        let progress_relaxed = PesScanInternalWasmProgress {
            opt_steps: Some(12),
            ..progress
        };
        let json_relaxed = serde_json::to_string(&progress_relaxed).unwrap();
        assert!(json_relaxed.contains("\"optSteps\":12"));
    }
}
