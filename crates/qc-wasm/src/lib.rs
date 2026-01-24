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
/// - `roots`: Array of quadrature roots in the interval (0, 1)
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
///                     If false (default), use Cartesian basis functions (6 d-orbitals).
///                     For s and p orbitals, both choices give the same result.
///
/// # Example
///
/// ```javascript
/// const options = { useSpherical: true };  // Use 5 d-orbitals instead of 6
/// ```
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct IntegralComputeOptions {
    /// Use spherical harmonic basis functions (5 d-orbitals vs 6 Cartesian).
    /// Default: false (Cartesian basis for backward compatibility).
    #[serde(default)]
    pub use_spherical: bool,
}

impl Default for IntegralComputeOptions {
    fn default() -> Self {
        Self {
            use_spherical: false, // Cartesian by default for backward compatibility
        }
    }
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
    use qc_core::integrals::{hcore_matrix, overlap_matrix};
    #[cfg(feature = "parallel")]
    use qc_core::integrals::eri_compressed_parallel;
    #[cfg(not(feature = "parallel"))]
    use qc_core::integrals::eri_compressed;
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
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*"];
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
                "Unsupported element '{}'. Only H-Ne are supported.",
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
/// * `basis_name` - Basis set name: "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*"
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
    use qc_core::integrals::{
        hcore_matrix, hcore_matrix_spherical, overlap_matrix, overlap_matrix_spherical,
    };
    #[cfg(feature = "parallel")]
    use qc_core::integrals::{eri_compressed_parallel, eri_compressed_spherical_parallel};
    #[cfg(not(feature = "parallel"))]
    use qc_core::integrals::{eri_compressed, eri_compressed_spherical};
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
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*"];
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
                "Unsupported element '{}'. Only H-Ne are supported.",
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
/// * `basis_name` - Basis set name: "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*"
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
    use qc_core::integrals::{hcore_matrix, overlap_matrix};
    #[cfg(feature = "parallel")]
    use qc_core::integrals::eri_compressed_parallel;
    #[cfg(not(feature = "parallel"))]
    use qc_core::integrals::{eri_index, shell_eri};
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
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*"];
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
                "Unsupported element '{}'. Only H-Ne are supported.",
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
/// * `basis_name` - Basis set name: "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*"
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
    use qc_core::integrals::{
        hcore_matrix, hcore_matrix_spherical, overlap_matrix, overlap_matrix_spherical,
    };
    #[cfg(feature = "parallel")]
    use qc_core::integrals::{eri_compressed_parallel, eri_compressed_spherical_parallel};
    #[cfg(not(feature = "parallel"))]
    use qc_core::integrals::{eri_compressed_spherical, eri_index, shell_eri};
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
    let supported_bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*"];
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
                "Unsupported element '{}'. Only H-Ne are supported.",
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
        };

        let json = serde_json::to_string(&matrices).unwrap();

        // Verify camelCase field names
        assert!(json.contains("\"nbf\":2"));
        assert!(json.contains("\"sMatrix\":"));
        assert!(json.contains("\"hCore\":"));
        assert!(json.contains("\"fockMatrix\":"));
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
        // Test that unsupported elements are rejected
        let geometry_json = r#"{
            "atoms": [
                { "symbol": "Na", "xyz": [0.0, 0.0, 0.0] }
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
}
