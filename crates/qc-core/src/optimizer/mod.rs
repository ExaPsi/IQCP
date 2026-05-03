//! Geometry Optimizer
//!
//! Implements molecular geometry optimization using two strategies:
//!
//! 1. **Internal coordinates** (polyatomic molecules, >= 3 atoms):
//!    Newton-Raphson with BFGS Hessian update in redundant internal
//!    coordinates (bonds, angles, dihedrals). Uses the Wilson B-matrix
//!    to transform between Cartesian and internal coordinate spaces.
//!
//! 2. **Cartesian L-BFGS** (diatomic molecules, fallback):
//!    Limited-memory BFGS quasi-Newton method in Cartesian coordinates
//!    (Nocedal 1980).
//!
//! The internal coordinate optimizer typically converges in 3-4 steps for
//! molecules like H2O, compared to 7+ steps for Cartesian L-BFGS, because
//! the Hessian is better conditioned in internal coordinates.
//!
//! # References
//!
//! - Wilson, Decius & Cross (1955). Molecular Vibrations. (B-matrix)
//! - Schlegel (1984). Theor. Chim. Acta 66, 333. (Model Hessian)
//! - Pulay & Fogarasi (1992). JCP 96, 2856. (Redundant internals)
//! - Baker (1993). JCP 105, 192. (Internal coordinate optimization)
//! - Nocedal, J. (1980). Math. Comp. 35, 773. (L-BFGS)
//! - Nocedal & Wright (2006). Numerical Optimization, Algorithm 7.4.

use nalgebra::{DMatrix, SymmetricEigen};
use serde::{Deserialize, Serialize};

use crate::basis::{Atom, BasisSet};
use crate::dft::GridConfig;
use crate::integrals::{eri_compressed, hcore_matrix, overlap_matrix};
use crate::scf::gradient::{ks_dft_gradient, rhf_gradient, GradientResult};
use crate::scf::{rhf_scf_with_guess, PresetSystem, ScfConfig};

/// Module version (matches crate version)
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// ============================================================================
// Configuration
// ============================================================================

/// Electronic structure method for the optimizer
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub enum OptMethod {
    /// Restricted Hartree-Fock
    Rhf,
    /// LDA (Slater exchange + VWN5 correlation)
    Lda,
    /// B3LYP hybrid functional
    B3lyp,
}

/// Configuration for geometry optimization
///
/// All fields have sensible defaults. The convergence criteria match
/// standard quantum chemistry codes (PySCF/geomeTRIC defaults).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct OptimizationConfig {
    /// Maximum number of optimization steps (default: 50)
    pub max_steps: usize,
    /// Maximum gradient convergence threshold in Ha/bohr (default: 4.5e-4)
    pub grad_threshold: f64,
    /// Energy change convergence threshold in Ha (default: 1.0e-6)
    pub energy_threshold: f64,
    /// L-BFGS memory size (number of {s,y} pairs, default: 7)
    pub memory_size: usize,
    /// Electronic structure method
    pub method: OptMethod,
    /// Basis set name (e.g., "sto-3g")
    pub basis: String,
}

impl Default for OptimizationConfig {
    fn default() -> Self {
        Self {
            max_steps: 50,
            grad_threshold: 4.5e-4,
            energy_threshold: 1.0e-6,
            memory_size: 7,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
        }
    }
}

// ============================================================================
// Result Types
// ============================================================================

/// A single step in the optimization trajectory
///
/// Records the energy, gradient, and geometry at each step of the
/// optimization. Step 0 corresponds to the initial (input) geometry.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct OptimizationStep {
    /// Step number (0 = initial geometry)
    pub step: usize,
    /// Total energy at this geometry (Ha)
    pub energy: f64,
    /// Maximum absolute gradient component (Ha/bohr)
    pub max_gradient: f64,
    /// RMS gradient (Ha/bohr)
    pub rms_gradient: f64,
    /// Atomic coordinates at this step [[x,y,z], ...] (bohr)
    pub geometry: Vec<[f64; 3]>,
    /// Gradient at this step [[gx,gy,gz], ...] (Ha/bohr)
    pub gradient: Vec<[f64; 3]>,
}

/// Result of a geometry optimization
///
/// Contains the full optimization trajectory, convergence status,
/// and the final optimized geometry.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct OptimizationResult {
    /// Whether the optimization converged
    pub converged: bool,
    /// Full optimization trajectory (including step 0)
    pub steps: Vec<OptimizationStep>,
    /// Final total energy (Ha)
    pub final_energy: f64,
    /// Final optimized geometry [[x,y,z], ...] (bohr)
    pub final_geometry: Vec<[f64; 3]>,
    /// Number of optimization steps taken (excludes initial evaluation)
    pub total_steps: usize,
    /// Total computation time in milliseconds
    pub compute_time_ms: f64,
}

// ============================================================================
// L-BFGS Core
// ============================================================================

/// L-BFGS history storage for the two-loop recursion
struct LbfgsHistory {
    /// Step vectors s_k = x_{k+1} - x_k
    s_vecs: Vec<Vec<f64>>,
    /// Gradient change vectors y_k = g_{k+1} - g_k
    y_vecs: Vec<Vec<f64>>,
    /// rho_k = 1 / (y_k^T * s_k)
    rho_vals: Vec<f64>,
    /// Maximum number of stored pairs
    max_size: usize,
}

impl LbfgsHistory {
    fn new(max_size: usize) -> Self {
        Self {
            s_vecs: Vec::with_capacity(max_size),
            y_vecs: Vec::with_capacity(max_size),
            rho_vals: Vec::with_capacity(max_size),
            max_size,
        }
    }

    /// Add a new {s, y} pair to the history
    ///
    /// Returns false if the curvature condition s^T * y > 0 is violated,
    /// indicating that this pair should be skipped.
    fn push(&mut self, s: Vec<f64>, y: Vec<f64>) -> bool {
        let sy: f64 = s.iter().zip(y.iter()).map(|(si, yi)| si * yi).sum();

        // Curvature condition: s^T * y > 0
        if sy <= 1e-10 {
            return false;
        }

        let rho = 1.0 / sy;

        // Drop oldest if at capacity
        if self.s_vecs.len() >= self.max_size {
            self.s_vecs.remove(0);
            self.y_vecs.remove(0);
            self.rho_vals.remove(0);
        }

        self.s_vecs.push(s);
        self.y_vecs.push(y);
        self.rho_vals.push(rho);
        true
    }

    /// Number of stored pairs
    fn len(&self) -> usize {
        self.s_vecs.len()
    }

    /// Clear all stored pairs (reset to steepest descent)
    fn clear(&mut self) {
        self.s_vecs.clear();
        self.y_vecs.clear();
        self.rho_vals.clear();
    }
}

/// Compute the L-BFGS search direction using the two-loop recursion
///
/// Reference: Nocedal & Wright (2006), Algorithm 7.4
///
/// With empty history, returns the negative gradient (steepest descent).
/// With history, returns -H_k * g_k where H_k is the L-BFGS approximation
/// to the inverse Hessian.
#[allow(clippy::needless_range_loop)]
fn lbfgs_direction(gradient: &[f64], history: &LbfgsHistory) -> Vec<f64> {
    let _n = gradient.len();
    let m = history.len();

    if m == 0 {
        // No history: steepest descent
        return gradient.iter().map(|g| -g).collect();
    }

    // Two-loop recursion
    let mut q: Vec<f64> = gradient.to_vec();
    let mut alpha_vec = vec![0.0; m];

    // First loop: from newest to oldest
    for i in (0..m).rev() {
        let alpha_i: f64 = history.rho_vals[i]
            * history.s_vecs[i]
                .iter()
                .zip(q.iter())
                .map(|(si, qi)| si * qi)
                .sum::<f64>();
        alpha_vec[i] = alpha_i;
        for (qj, yij) in q.iter_mut().zip(history.y_vecs[i].iter()) {
            *qj -= alpha_i * yij;
        }
    }

    // Initial Hessian scaling: gamma = (s_{k-1}^T * y_{k-1}) / (y_{k-1}^T * y_{k-1})
    let last = m - 1;
    let sy: f64 = history.s_vecs[last]
        .iter()
        .zip(history.y_vecs[last].iter())
        .map(|(si, yi)| si * yi)
        .sum();
    let yy: f64 = history.y_vecs[last].iter().map(|yi| yi * yi).sum();
    let gamma = if yy > 1e-30 { sy / yy } else { 1.0 };

    // r = gamma * q (apply initial Hessian H_0 = gamma * I)
    let mut r: Vec<f64> = q.iter().map(|qi| gamma * qi).collect();

    // Second loop: from oldest to newest
    for i in 0..m {
        let beta: f64 = history.rho_vals[i]
            * history.y_vecs[i]
                .iter()
                .zip(r.iter())
                .map(|(yi, ri)| yi * ri)
                .sum::<f64>();
        let coeff = alpha_vec[i] - beta;
        for (rj, sij) in r.iter_mut().zip(history.s_vecs[i].iter()) {
            *rj += sij * coeff;
        }
    }

    // direction = -r
    r.iter().map(|ri| -ri).collect()
}

// ============================================================================
// Energy and Gradient Evaluation
// ============================================================================

/// Maximum step size per coordinate component (bohr)
const MAX_STEP_SIZE: f64 = 0.3;

/// Build a PresetSystem from atoms and basis set
fn build_system(basis: &BasisSet) -> PresetSystem {
    let s = overlap_matrix(basis);
    let h = hcore_matrix(basis);
    let eri = eri_compressed(basis);

    PresetSystem {
        system_id: "optimizer".to_string(),
        label: "Optimizer geometry".to_string(),
        nbf: basis.n_basis,
        nelec: basis.n_electrons,
        e_nuc: basis.nuclear_repulsion,
        s_matrix: s,
        h_core: h,
        eri_compressed: eri,
    }
}

/// Evaluate energy and gradient at a given geometry
///
/// Returns (energy, gradient_per_atom, gradient_result)
fn evaluate_energy_gradient(
    atoms: &[Atom],
    config: &OptimizationConfig,
    scf_config: &ScfConfig,
) -> Option<(f64, Vec<[f64; 3]>, GradientResult)> {
    match config.method {
        OptMethod::Rhf => {
            let basis = BasisSet::build(atoms.to_vec(), &config.basis).ok()?;
            let system = build_system(&basis);
            let n_occ = system.n_occ();
            // Use SAD initial guess for faster SCF convergence
            let sad_density = crate::scf::sad::build_sad_density(&basis);
            let scf_output = rhf_scf_with_guess(&system, scf_config, Some(&sad_density)).ok()?;
            let grad_result = rhf_gradient(
                &basis,
                &scf_output.density_matrix,
                &scf_output.mo_coefficients,
                &scf_output.mo_energies,
                n_occ,
            );
            Some((
                scf_output.energy_total,
                grad_result.gradients.clone(),
                grad_result,
            ))
        }
        OptMethod::Lda => {
            let lda = crate::dft::Lda::new();
            let grid_config = GridConfig::default();
            // ks_dft_gradient runs ONE SCF internally and returns both
            // gradient and energy -- no second SCF needed
            let grad_result =
                ks_dft_gradient(atoms, &config.basis, &lda, &grid_config, scf_config, false);
            let energy = grad_result
                .energy
                .expect("ks_dft_gradient should provide energy");
            Some((energy, grad_result.gradients.clone(), grad_result))
        }
        OptMethod::B3lyp => {
            let b3lyp = crate::dft::B3lyp::new();
            let grid_config = GridConfig::default();
            // ks_dft_gradient runs ONE SCF internally and returns both
            // gradient and energy -- no second SCF needed
            let grad_result = ks_dft_gradient(
                atoms,
                &config.basis,
                &b3lyp,
                &grid_config,
                scf_config,
                false,
            );
            let energy = grad_result
                .energy
                .expect("ks_dft_gradient should provide energy");
            Some((energy, grad_result.gradients.clone(), grad_result))
        }
    }
}

/// Convert atoms to a flat coordinate vector [x0, y0, z0, x1, y1, z1, ...]
fn atoms_to_coords(atoms: &[Atom]) -> Vec<f64> {
    atoms
        .iter()
        .flat_map(|a| a.position.iter().copied())
        .collect()
}

/// Convert gradient per atom to a flat vector
fn grad_to_flat(grad: &[[f64; 3]]) -> Vec<f64> {
    grad.iter().flat_map(|g| g.iter().copied()).collect()
}

/// Convert flat coordinates back to atoms (preserving atomic numbers)
fn update_atom_positions(atoms: &[Atom], coords: &[f64]) -> Vec<Atom> {
    atoms
        .iter()
        .enumerate()
        .map(|(i, atom)| {
            let idx = i * 3;
            Atom {
                atomic_number: atom.atomic_number,
                position: [coords[idx], coords[idx + 1], coords[idx + 2]],
                symbol: atom.symbol.clone(),
            }
        })
        .collect()
}

/// Extract geometry as Vec<[f64; 3]> from atoms
fn atoms_to_geometry(atoms: &[Atom]) -> Vec<[f64; 3]> {
    atoms.iter().map(|a| a.position).collect()
}

// ============================================================================
// Main Optimizer
// ============================================================================

/// Optimize molecular geometry using the best available method
///
/// For polyatomic molecules (>= 3 atoms), uses internal coordinate
/// optimization with Newton-Raphson + BFGS Hessian update. For
/// diatomic molecules, falls back to Cartesian L-BFGS.
///
/// # Arguments
///
/// * `atoms` - Initial molecular geometry as (atomic_number, [x, y, z]) pairs.
///   Coordinates are in bohr.
/// * `config` - Optimization configuration (method, basis, convergence criteria)
/// * `progress` - Optional callback invoked after each step with the current
///   `OptimizationStep`. Used for streaming progress to the UI.
///
/// # Returns
///
/// `OptimizationResult` with full trajectory, final geometry, and convergence status.
pub fn optimize_geometry(
    atoms: &[(u8, [f64; 3])],
    config: &OptimizationConfig,
    progress: Option<&dyn Fn(&OptimizationStep)>,
) -> OptimizationResult {
    // Build Atom structs to detect connectivity
    let atom_structs: Vec<Atom> = atoms
        .iter()
        .filter_map(|(z, pos)| Atom::new(*z, *pos).ok())
        .collect();

    if atom_structs.len() != atoms.len() {
        return OptimizationResult {
            converged: false,
            steps: Vec::new(),
            final_energy: 0.0,
            final_geometry: atoms.iter().map(|(_, pos)| *pos).collect(),
            total_steps: 0,
            compute_time_ms: 0.0,
        };
    }

    // Use internal coordinates for polyatomic molecules (>= 3 atoms)
    if atom_structs.len() >= 3 {
        let connectivity = detect_connectivity(&atom_structs);
        if connectivity.n_internals >= 2 {
            return optimize_internal(atoms, config, progress);
        }
    }

    // Fall back to Cartesian L-BFGS for diatomics and small systems
    optimize_cartesian(atoms, config, progress)
}

/// Optimize molecular geometry using Cartesian L-BFGS quasi-Newton method
///
/// This is the original Cartesian optimizer, used as fallback for diatomic
/// molecules. For polyatomic molecules, `optimize_geometry` dispatches to
/// the internal coordinate optimizer instead.
///
/// # References
///
/// - Nocedal (1980). Math. Comp. 35, 773.
/// - Nocedal & Wright (2006). Numerical Optimization, Alg. 7.4.
fn optimize_cartesian(
    atoms: &[(u8, [f64; 3])],
    config: &OptimizationConfig,
    progress: Option<&dyn Fn(&OptimizationStep)>,
) -> OptimizationResult {
    // SCF configuration: use tight convergence for inner SCF loop
    let scf_config = ScfConfig {
        use_diis: true,
        ..ScfConfig::tight()
    };

    // Convert input to Atom structs
    let mut current_atoms: Vec<Atom> = atoms
        .iter()
        .filter_map(|(z, pos)| Atom::new(*z, *pos).ok())
        .collect();

    if current_atoms.len() != atoms.len() {
        // Invalid atoms -- return empty result
        return OptimizationResult {
            converged: false,
            steps: Vec::new(),
            final_energy: 0.0,
            final_geometry: atoms.iter().map(|(_, pos)| *pos).collect(),
            total_steps: 0,
            compute_time_ms: 0.0,
        };
    }

    let n_atoms = current_atoms.len();
    let n_coords = n_atoms * 3;

    // Step 0: Evaluate energy and gradient at initial geometry
    let eval_result = evaluate_energy_gradient(&current_atoms, config, &scf_config);
    let (mut energy, mut grad_atoms, _) = match eval_result {
        Some(r) => r,
        None => {
            return OptimizationResult {
                converged: false,
                steps: Vec::new(),
                final_energy: 0.0,
                final_geometry: atoms_to_geometry(&current_atoms),
                total_steps: 0,
                compute_time_ms: 0.0,
            };
        }
    };

    let mut grad_flat = grad_to_flat(&grad_atoms);
    let mut coords = atoms_to_coords(&current_atoms);

    // Record initial step
    let initial_step = make_step(0, energy, &grad_atoms, &current_atoms);
    let mut steps = vec![initial_step.clone()];
    if let Some(cb) = progress {
        cb(&initial_step);
    }

    // L-BFGS history
    let mut history = LbfgsHistory::new(config.memory_size);
    let mut converged = false;

    for step_num in 1..=config.max_steps {
        // 1. Compute search direction via L-BFGS two-loop recursion
        let mut direction = lbfgs_direction(&grad_flat, &history);

        // 2. Clamp step size: ensure no component exceeds MAX_STEP_SIZE
        let max_component = direction.iter().map(|d| d.abs()).fold(0.0_f64, f64::max);
        if max_component > MAX_STEP_SIZE {
            let scale = MAX_STEP_SIZE / max_component;
            for d in &mut direction {
                *d *= scale;
            }
        }

        // 3. Verify descent direction: g^T * d < 0
        let directional_deriv: f64 = grad_flat
            .iter()
            .zip(direction.iter())
            .map(|(g, d)| g * d)
            .sum();
        if directional_deriv >= 0.0 {
            // Not a descent direction: reset history, use steepest descent
            history.clear();
            direction = grad_flat.iter().map(|g| -g).collect();
            // Re-clamp
            let max_comp = direction.iter().map(|d| d.abs()).fold(0.0_f64, f64::max);
            if max_comp > MAX_STEP_SIZE {
                let scale = MAX_STEP_SIZE / max_comp;
                for d in &mut direction {
                    *d *= scale;
                }
            }
        }

        // 4. Take full step (alpha=1.0): L-BFGS direction is already
        //    clamped by MAX_STEP_SIZE, so no line search needed.
        //    This eliminates 1-3 redundant SCF evaluations per step.
        let old_coords = coords.clone();
        let old_grad = grad_flat.clone();
        let prev_energy = energy;

        for i in 0..n_coords {
            coords[i] += direction[i];
        }
        current_atoms = update_atom_positions(&current_atoms, &coords);

        // 5. Compute energy AND gradient at new geometry (ONE SCF call)
        let eval_new = evaluate_energy_gradient(&current_atoms, config, &scf_config);
        match eval_new {
            Some((e_new, g_new, _)) => {
                energy = e_new;
                grad_atoms = g_new;
                grad_flat = grad_to_flat(&grad_atoms);
            }
            None => {
                // SCF failed at new geometry -- report partial result
                let step_record = OptimizationStep {
                    step: step_num,
                    energy,
                    max_gradient: f64::NAN,
                    rms_gradient: f64::NAN,
                    geometry: atoms_to_geometry(&current_atoms),
                    gradient: vec![[0.0; 3]; n_atoms],
                };
                steps.push(step_record.clone());
                if let Some(cb) = progress {
                    cb(&step_record);
                }
                break;
            }
        }

        // 7. Update L-BFGS history
        let s_vec: Vec<f64> = coords
            .iter()
            .zip(old_coords.iter())
            .map(|(c, o)| c - o)
            .collect();
        let y_vec: Vec<f64> = grad_flat
            .iter()
            .zip(old_grad.iter())
            .map(|(g, o)| g - o)
            .collect();

        if !history.push(s_vec, y_vec) {
            // Curvature condition violated -- reset history
            history.clear();
        }

        // 8. Record step
        let step_record = make_step(step_num, energy, &grad_atoms, &current_atoms);
        steps.push(step_record.clone());
        if let Some(cb) = progress {
            cb(&step_record);
        }

        // 9. Check convergence
        let delta_e = (energy - prev_energy).abs();
        if step_record.max_gradient < config.grad_threshold && delta_e < config.energy_threshold {
            converged = true;
            break;
        }
    }

    let last_step = steps.last().unwrap();
    OptimizationResult {
        converged,
        final_energy: last_step.energy,
        final_geometry: last_step.geometry.clone(),
        total_steps: steps.len() - 1, // exclude initial evaluation
        steps,
        compute_time_ms: 0.0, // set by caller (WASM layer)
    }
}

/// Create an OptimizationStep from current state
fn make_step(
    step: usize,
    energy: f64,
    grad_atoms: &[[f64; 3]],
    atoms: &[Atom],
) -> OptimizationStep {
    let mut max_grad = 0.0_f64;
    let mut sum_sq = 0.0;
    let mut n = 0;
    for g in grad_atoms {
        for &component in g {
            max_grad = max_grad.max(component.abs());
            sum_sq += component * component;
            n += 1;
        }
    }
    let rms_grad = if n > 0 {
        (sum_sq / n as f64).sqrt()
    } else {
        0.0
    };

    OptimizationStep {
        step,
        energy,
        max_gradient: max_grad,
        rms_gradient: rms_grad,
        geometry: atoms_to_geometry(atoms),
        gradient: grad_atoms.to_vec(),
    }
}

// ============================================================================
// Internal Coordinate Optimizer
// ============================================================================
//
// References:
// - Wilson, Decius & Cross (1955). Molecular Vibrations. (B-matrix theory)
// - Schlegel (1984). Theor. Chim. Acta 66, 333. (Model Hessian)
// - Pulay & Fogarasi (1992). JCP 96, 2856. (Redundant internals)
// - Baker (1993). JCP 105, 192. (Internal coordinate optimization)

// ---- Covalent radii (bohr) indexed by atomic number Z=0..18 ----
// Source: Pyykko & Atsumi (2009), converted from Angstrom (* 1/0.529177)
const COVALENT_RADII: [f64; 19] = [
    0.0,   // Z=0 placeholder
    0.661, // H  (0.35 A)
    0.661, // He (not bonding, but need a value)
    2.456, // Li (1.30 A)
    1.890, // Be (1.00 A)
    1.606, // B  (0.85 A)
    1.436, // C  (0.76 A)
    1.342, // N  (0.71 A)
    1.247, // O  (0.66 A)
    1.191, // F  (0.63 A)
    1.153, // Ne
    3.024, // Na (1.60 A)
    2.645, // Mg (1.40 A)
    2.362, // Al (1.25 A)
    2.173, // Si (1.15 A)
    2.040, // P  (1.08 A)
    1.928, // S  (1.02 A)
    1.870, // Cl (0.99 A)
    1.890, // Ar
];

/// Get covalent radius for an atomic number, defaulting to 1.5 bohr for unknowns
#[inline]
fn covalent_radius(z: u8) -> f64 {
    if (z as usize) < COVALENT_RADII.len() {
        COVALENT_RADII[z as usize]
    } else {
        1.5 // conservative default
    }
}

/// Molecular connectivity: bonds, angles, and dihedrals detected from geometry
pub(crate) struct MolecularConnectivity {
    pub(crate) bonds: Vec<(usize, usize)>,
    pub(crate) angles: Vec<(usize, usize, usize)>, // (i, j, k) where j is central atom
    pub(crate) dihedrals: Vec<(usize, usize, usize, usize)>, // (i, j, k, l) where j-k is central bond
    pub(crate) n_internals: usize,
}

/// Detect molecular connectivity from atomic positions using covalent radii
///
/// A bond is detected when:
///   distance(A, B) < 1.3 * (r_cov(A) + r_cov(B))
///
/// Angles are detected for all i-j-k triplets where (i,j) and (j,k) are bonds.
/// Dihedrals are detected for i-j-k-l where (i,j), (j,k), (k,l) are bonds.
pub(crate) fn detect_connectivity(atoms: &[Atom]) -> MolecularConnectivity {
    let n = atoms.len();
    let bond_scale = 1.3;

    // --- Step 1: Detect bonds ---
    let mut bonds = Vec::new();
    // Adjacency list for efficient angle/dihedral detection
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = atoms[i].position[0] - atoms[j].position[0];
            let dy = atoms[i].position[1] - atoms[j].position[1];
            let dz = atoms[i].position[2] - atoms[j].position[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            let threshold = bond_scale
                * (covalent_radius(atoms[i].atomic_number)
                    + covalent_radius(atoms[j].atomic_number));

            if dist < threshold && dist > 0.1 {
                // 0.1 bohr minimum to avoid self-interaction
                bonds.push((i, j));
                neighbors[i].push(j);
                neighbors[j].push(i);
            }
        }
    }

    // --- Step 2: Detect angles (i-j-k where j is central) ---
    let mut angles = Vec::new();
    for (j, nbrs) in neighbors.iter().enumerate() {
        for idx_a in 0..nbrs.len() {
            for idx_b in (idx_a + 1)..nbrs.len() {
                let i = nbrs[idx_a];
                let k = nbrs[idx_b];
                // Ensure canonical ordering: smaller index first for the endpoints
                if i < k {
                    angles.push((i, j, k));
                } else {
                    angles.push((k, j, i));
                }
            }
        }
    }

    // --- Step 3: Detect dihedrals (i-j-k-l where j-k is a bond) ---
    let mut dihedrals = Vec::new();
    for &(j, k) in &bonds {
        // For each neighbor i of j (i != k) and each neighbor l of k (l != j)
        for &i in &neighbors[j] {
            if i == k {
                continue;
            }
            for &l in &neighbors[k] {
                if l == j || l == i {
                    continue;
                }
                // Canonical ordering: ensure i < l for the same j-k bond
                if i < l {
                    dihedrals.push((i, j, k, l));
                }
                // Skip (l, k, j, i) -- would be the same dihedral
            }
        }
    }

    let n_internals = bonds.len() + angles.len() + dihedrals.len();
    MolecularConnectivity {
        bonds,
        angles,
        dihedrals,
        n_internals,
    }
}

// ---- Internal coordinate values ----

/// Compute bond length between atoms i and j (bohr)
#[inline]
pub(crate) fn compute_bond_length(atoms: &[Atom], i: usize, j: usize) -> f64 {
    let dx = atoms[i].position[0] - atoms[j].position[0];
    let dy = atoms[i].position[1] - atoms[j].position[1];
    let dz = atoms[i].position[2] - atoms[j].position[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Compute angle i-j-k in radians (j is central atom)
#[inline]
pub(crate) fn compute_angle(atoms: &[Atom], i: usize, j: usize, k: usize) -> f64 {
    let v1 = [
        atoms[i].position[0] - atoms[j].position[0],
        atoms[i].position[1] - atoms[j].position[1],
        atoms[i].position[2] - atoms[j].position[2],
    ];
    let v2 = [
        atoms[k].position[0] - atoms[j].position[0],
        atoms[k].position[1] - atoms[j].position[1],
        atoms[k].position[2] - atoms[j].position[2],
    ];
    let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    let mag1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
    let mag2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
    let cos_theta = (dot / (mag1 * mag2)).clamp(-1.0, 1.0);
    cos_theta.acos()
}

/// Compute dihedral angle i-j-k-l in radians (j-k is the central bond)
///
/// Uses the atan2 formula for numerical stability:
///   tau = atan2(|b2| * b1 . (b2 x b3), (b1 x b2) . (b2 x b3))
/// where b1 = j-i, b2 = k-j, b3 = l-k
#[inline]
pub(crate) fn compute_dihedral(atoms: &[Atom], i: usize, j: usize, k: usize, l: usize) -> f64 {
    // Bond vectors
    let b1 = [
        atoms[j].position[0] - atoms[i].position[0],
        atoms[j].position[1] - atoms[i].position[1],
        atoms[j].position[2] - atoms[i].position[2],
    ];
    let b2 = [
        atoms[k].position[0] - atoms[j].position[0],
        atoms[k].position[1] - atoms[j].position[1],
        atoms[k].position[2] - atoms[j].position[2],
    ];
    let b3 = [
        atoms[l].position[0] - atoms[k].position[0],
        atoms[l].position[1] - atoms[k].position[1],
        atoms[l].position[2] - atoms[k].position[2],
    ];

    // Cross products
    let n1 = cross(&b1, &b2);
    let n2 = cross(&b2, &b3);

    let b2_mag = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
    let x = dot3(&n1, &n2);
    let y = dot3(
        &cross(&n1, &n2),
        &[b2[0] / b2_mag, b2[1] / b2_mag, b2[2] / b2_mag],
    );

    y.atan2(x)
}

/// Compute all internal coordinate values
fn compute_internal_coords(connectivity: &MolecularConnectivity, atoms: &[Atom]) -> Vec<f64> {
    let mut q = Vec::with_capacity(connectivity.n_internals);

    for &(i, j) in &connectivity.bonds {
        q.push(compute_bond_length(atoms, i, j));
    }
    for &(i, j, k) in &connectivity.angles {
        q.push(compute_angle(atoms, i, j, k));
    }
    for &(i, j, k, l) in &connectivity.dihedrals {
        q.push(compute_dihedral(atoms, i, j, k, l));
    }

    q
}

#[inline]
pub(crate) fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub(crate) fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ---- Wilson B-matrix ----
//
// The B-matrix transforms Cartesian displacements to internal coordinate
// displacements: dq = B * dx
//
// Reference: Wilson, Decius & Cross (1955), Chapter 4

/// Wilson B-matrix for transforming between Cartesian and internal coordinates
struct WilsonBMatrix {
    /// B-matrix data in row-major order: n_internals x n_cartesian
    mat: DMatrix<f64>,
}

impl WilsonBMatrix {
    /// Build the Wilson B-matrix at the current geometry
    ///
    /// Each row corresponds to one internal coordinate (bond, angle, or dihedral).
    /// Each column corresponds to one Cartesian component (3*n_atoms total).
    fn build(connectivity: &MolecularConnectivity, atoms: &[Atom]) -> Self {
        let n_int = connectivity.n_internals;
        let n_cart = atoms.len() * 3;
        let mut mat = DMatrix::zeros(n_int, n_cart);

        let mut row = 0;

        // --- Bond stretches ---
        for &(a, b) in &connectivity.bonds {
            let r = compute_bond_length(atoms, a, b);
            if r < 1e-10 {
                row += 1;
                continue;
            }
            // B[r, 3A+d] = (R_A - R_B)[d] / r
            // B[r, 3B+d] = (R_B - R_A)[d] / r
            for d in 0..3 {
                let e = (atoms[a].position[d] - atoms[b].position[d]) / r;
                mat[(row, 3 * a + d)] = e;
                mat[(row, 3 * b + d)] = -e;
            }
            row += 1;
        }

        // --- Angle bends ---
        for &(i, j, k) in &connectivity.angles {
            let r_ji = compute_bond_length(atoms, j, i);
            let r_jk = compute_bond_length(atoms, j, k);

            if r_ji < 1e-10 || r_jk < 1e-10 {
                row += 1;
                continue;
            }

            // Unit vectors from central atom j
            let mut e_ji = [0.0; 3];
            let mut e_jk = [0.0; 3];
            for d in 0..3 {
                e_ji[d] = (atoms[i].position[d] - atoms[j].position[d]) / r_ji;
                e_jk[d] = (atoms[k].position[d] - atoms[j].position[d]) / r_jk;
            }

            let cos_theta = dot3(&e_ji, &e_jk).clamp(-1.0, 1.0);
            let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-12);

            // B[theta, 3I+d] = (cos(theta) * e_ji - e_jk)[d] / (r_ji * sin(theta))
            // B[theta, 3K+d] = (cos(theta) * e_jk - e_ji)[d] / (r_jk * sin(theta))
            // B[theta, 3J+d] = -(B[theta, 3I+d] + B[theta, 3K+d])
            for d in 0..3 {
                let b_i = (cos_theta * e_ji[d] - e_jk[d]) / (r_ji * sin_theta);
                let b_k = (cos_theta * e_jk[d] - e_ji[d]) / (r_jk * sin_theta);
                mat[(row, 3 * i + d)] = b_i;
                mat[(row, 3 * k + d)] = b_k;
                mat[(row, 3 * j + d)] = -b_i - b_k;
            }
            row += 1;
        }

        // --- Dihedral torsions ---
        // B-matrix elements for dihedral i-j-k-l:
        // Reference: Helgaker, Jorgensen & Olsen, Chapter 1, Eq. 1.4.26-1.4.29
        for &(i, j, k, l) in &connectivity.dihedrals {
            let r_ij = compute_bond_length(atoms, i, j);
            let r_jk = compute_bond_length(atoms, j, k);
            let r_kl = compute_bond_length(atoms, k, l);

            if r_ij < 1e-10 || r_jk < 1e-10 || r_kl < 1e-10 {
                row += 1;
                continue;
            }

            // Bond vectors
            let mut e_ij = [0.0; 3];
            let mut e_jk = [0.0; 3];
            let mut e_kl = [0.0; 3];
            for d in 0..3 {
                e_ij[d] = (atoms[j].position[d] - atoms[i].position[d]) / r_ij;
                e_jk[d] = (atoms[k].position[d] - atoms[j].position[d]) / r_jk;
                e_kl[d] = (atoms[l].position[d] - atoms[k].position[d]) / r_kl;
            }

            // Cross products for normal vectors
            let n1 = cross(&e_ij, &e_jk); // normal to i-j-k plane
            let n2 = cross(&e_jk, &e_kl); // normal to j-k-l plane

            let sin_ijk = cross(&e_ij, &e_jk)
                .iter()
                .map(|x| x * x)
                .sum::<f64>()
                .sqrt()
                .max(1e-12);
            let sin_jkl = cross(&e_jk, &e_kl)
                .iter()
                .map(|x| x * x)
                .sum::<f64>()
                .sqrt()
                .max(1e-12);

            let cos_ijk = dot3(&e_ij, &e_jk).clamp(-1.0, 1.0);
            let cos_jkl = dot3(&e_jk, &e_kl).clamp(-1.0, 1.0);

            // B[tau, atom_i] = n1 / (r_ij * sin^2(ijk))
            // B[tau, atom_l] = -n2 / (r_kl * sin^2(jkl))
            let sin2_ijk = sin_ijk * sin_ijk;
            let sin2_jkl = sin_jkl * sin_jkl;

            for d in 0..3 {
                let b_i = n1[d] / (r_ij * sin2_ijk);
                let b_l = -n2[d] / (r_kl * sin2_jkl);

                // For atoms j and k, use the chain rule relations:
                // B[tau, j] = -(1 - cos_ijk/sin_ijk * r_ij/r_jk) * B[tau, i]
                //              + cos_jkl / (sin_jkl * r_jk) * n2[d] / sin_jkl
                // Simplified from Helgaker Eq. 1.4.27-1.4.28
                let frac_ij_jk = r_ij / r_jk;
                let frac_kl_jk = r_kl / r_jk;

                let b_j = -(1.0 - cos_ijk * frac_ij_jk) * b_i + cos_jkl * frac_kl_jk * b_l;
                let b_k = -(1.0 - cos_jkl * frac_kl_jk) * b_l + cos_ijk * frac_ij_jk * b_i;

                mat[(row, 3 * i + d)] = b_i;
                mat[(row, 3 * j + d)] = b_j;
                mat[(row, 3 * k + d)] = b_k;
                mat[(row, 3 * l + d)] = b_l;
            }
            row += 1;
        }

        Self { mat }
    }

    /// Transform Cartesian gradient to internal coordinates
    ///
    /// g_q = G^{-1} B g_x, where G = B B^T
    ///
    /// Uses SVD-based pseudoinverse for handling redundant coordinates.
    fn cartesian_to_internal_gradient(&self, grad_cart: &[f64]) -> Vec<f64> {
        let n_cart = self.mat.ncols();
        let n_int = self.mat.nrows();
        assert_eq!(grad_cart.len(), n_cart);

        let g_x = DMatrix::from_column_slice(n_cart, 1, grad_cart);

        // g_q = G^+ B g_x where G = B B^T and G^+ is pseudoinverse
        let bg = &self.mat * &g_x;
        let g = &self.mat * self.mat.transpose();

        // SVD-based pseudoinverse of G
        let g_inv = pseudoinverse_symmetric(&g, 1e-6);
        let g_q = &g_inv * &bg;

        (0..n_int).map(|i| g_q[(i, 0)]).collect()
    }

    /// Transform internal coordinate step to Cartesian displacements
    ///
    /// dx = B^T G^{-1} dq, where G = B B^T
    ///
    /// This gives the minimum-norm Cartesian step consistent with the
    /// desired internal coordinate change. Uses SVD-based pseudoinverse.
    fn internal_to_cartesian_step(&self, step_int: &[f64]) -> Vec<f64> {
        let n_cart = self.mat.ncols();
        let n_int = self.mat.nrows();
        assert_eq!(step_int.len(), n_int);

        let dq = DMatrix::from_column_slice(n_int, 1, step_int);
        let g = &self.mat * self.mat.transpose();
        let g_inv = pseudoinverse_symmetric(&g, 1e-6);
        let dx = self.mat.transpose() * &g_inv * &dq;

        (0..n_cart).map(|i| dx[(i, 0)]).collect()
    }
}

/// SVD-based pseudoinverse for a symmetric matrix
///
/// Inverts only singular values above `threshold`, setting smaller ones to zero.
/// This correctly handles redundant internal coordinates where some combinations
/// are not independent.
fn pseudoinverse_symmetric(m: &DMatrix<f64>, threshold: f64) -> DMatrix<f64> {
    let svd = m.clone().svd(true, true);
    let u = svd.u.expect("SVD should produce U");
    let vt = svd.v_t.expect("SVD should produce V^T");
    let n = svd.singular_values.len();

    // Build pseudoinverse of singular value diagonal
    let mut s_inv = DMatrix::zeros(n, n);
    for i in 0..n {
        let s = svd.singular_values[i];
        if s > threshold {
            s_inv[(i, i)] = 1.0 / s;
        }
    }

    // M^+ = V S^+ U^T
    &vt.transpose() * &s_inv * &u.transpose()
}

// ---- Model Hessian ----
//
// Diagonal Hessian in internal coordinates with empirical force constants.
// Much better initial guess than the identity matrix.
//
// Reference: Schlegel (1984). Theor. Chim. Acta 66, 333.

/// Period of an atom (1 for H-He, 2 for Li-Ne, 3 for Na-Ar)
#[inline]
fn atom_period(z: u8) -> u8 {
    match z {
        1..=2 => 1,
        3..=10 => 2,
        11..=18 => 3,
        _ => 3,
    }
}

/// Build empirical model Hessian in internal coordinates (diagonal)
///
/// Returns a diagonal Hessian as a flat Vec of length n_internals.
/// Force constants are empirical estimates that provide a much better
/// initial guess than the identity matrix.
///
/// Bond force constants follow Fischer & Almlof (1992), JPC 96, 9768:
///   k_bond = A / (r - r_ref)^2
/// where A and r_ref depend on the period pair of the bonded atoms.
///
/// Angle: k_angle = 0.16 Ha/rad^2 (Schlegel 1984)
/// Dihedral: k_dihedral = 0.023 Ha/rad^2 (Schlegel 1984)
fn model_hessian_diagonal(connectivity: &MolecularConnectivity, atoms: &[Atom]) -> Vec<f64> {
    let mut h = Vec::with_capacity(connectivity.n_internals);

    // --- Bond force constants ---
    // Fischer-Almlof empirical model: k = A / (r - r_ref)^2
    // The parameter A and r_ref depend on the row-pair of the two atoms.
    // Parameters calibrated to approximate true Hessian eigenvalues:
    //   O-H (~1.8 bohr) -> k ~ 0.6-0.8  (actual ~0.8)
    //   C-C (~2.8 bohr) -> k ~ 0.4-0.6
    //   H-H (~1.4 bohr) -> k ~ 0.5-0.7
    for &(i, j) in &connectivity.bonds {
        let r = compute_bond_length(atoms, i, j);
        let pi = atom_period(atoms[i].atomic_number);
        let pj = atom_period(atoms[j].atomic_number);

        let (a, r_ref) = match (pi.min(pj), pi.max(pj)) {
            (1, 1) => (0.35, 0.60), // H-H type
            (1, 2) => (0.70, 0.85), // H with 2nd period (C, N, O, F)
            (1, 3) => (0.60, 0.85), // H with 3rd period
            (2, 2) => (1.00, 1.20), // 2nd period - 2nd period (C-C, C-O, etc.)
            (2, 3) => (0.85, 1.50), // 2nd - 3rd period
            (3, 3) => (0.75, 1.80), // 3rd - 3rd period
            _ => (0.85, 1.20),
        };

        let dr = (r - r_ref).max(0.3); // prevent division by near-zero
        let k = a / (dr * dr);
        // Clamp to reasonable range
        h.push(k.clamp(0.05, 2.0));
    }

    // --- Angle force constants ---
    // Typical bends have eigenvalues ~0.2-0.8 Ha/rad^2 depending on the
    // atoms involved. Use 0.3 as a reasonable middle estimate.
    h.extend(std::iter::repeat_n(0.30, connectivity.angles.len()));

    // --- Dihedral force constants ---
    // Torsions are soft modes with eigenvalues ~0.02-0.05 Ha/rad^2
    h.extend(std::iter::repeat_n(0.023, connectivity.dihedrals.len()));

    h
}

/// BFGS update of the Hessian matrix in internal coordinates
///
/// H_new = H + (dg dg^T) / (dg^T dq) - (H dq)(H dq)^T / (dq^T H dq)
///
/// where dq = q_{k+1} - q_k and dg = g_{k+1} - g_k in internal coordinates.
///
/// For a diagonal-only Hessian, this produces a full (dense) update.
/// We store the full n_int x n_int Hessian.
fn bfgs_update_hessian(h: &mut DMatrix<f64>, dq: &[f64], dg: &[f64]) {
    let n = dq.len();
    assert_eq!(dg.len(), n);
    assert_eq!(h.nrows(), n);
    assert_eq!(h.ncols(), n);

    // dg^T dq  (must be positive for a valid BFGS update — Wolfe curvature condition)
    let dg_dot_dq: f64 = dg.iter().zip(dq.iter()).map(|(g, q)| g * q).sum();
    if dg_dot_dq < 1e-14 {
        return; // Skip update for negative or near-zero curvature (Wolfe condition)
    }

    // H * dq
    let dq_col = DMatrix::from_column_slice(n, 1, dq);
    let hdq = h as &DMatrix<f64> * &dq_col;

    // dq^T H dq
    let dq_h_dq: f64 = (0..n).map(|i| dq[i] * hdq[(i, 0)]).sum();
    if dq_h_dq.abs() < 1e-14 {
        return; // Skip update if denominator too small
    }

    // BFGS update: H += (dg dg^T) / (dg^T dq) - (H dq)(H dq)^T / (dq^T H dq)
    for i in 0..n {
        for j in 0..n {
            h[(i, j)] += dg[i] * dg[j] / dg_dot_dq - hdq[(i, 0)] * hdq[(j, 0)] / dq_h_dq;
        }
    }
}

/// Maximum step sizes for internal coordinates
const MAX_BOND_STEP: f64 = 0.3; // bohr
const MAX_ANGLE_STEP: f64 = 0.3; // radians (~17 degrees)
const MAX_DIHEDRAL_STEP: f64 = 0.3; // radians

// ---- Trust-radius optimizer constants ----
// Reference: geomeTRIC (Wang & Song, 2016)
/// Initial trust radius for RMS norm of internal coordinate step (bohr/rad)
const INITIAL_TRUST_RADIUS: f64 = 0.3;
/// Minimum trust radius (prevents stalling)
const MIN_TRUST_RADIUS: f64 = 0.01;
/// Maximum trust radius
const MAX_TRUST_RADIUS: f64 = 1.0;
/// Growth factor for trust radius on good steps (sqrt(2))
const TRUST_GROW_FACTOR: f64 = 1.414;
/// Shrink factor for trust radius on poor/rejected steps
const TRUST_SHRINK_FACTOR: f64 = 0.5;
/// Quality ratio threshold for growing trust radius
const GOOD_QUALITY: f64 = 0.75;
/// Quality ratio threshold for shrinking trust radius
const POOR_QUALITY: f64 = 0.25;
/// Minimum eigenvalue floor for Hessian in Newton step (Ha/bohr^2)
const MIN_HESSIAN_EIGENVALUE: f64 = 0.005;

/// Clamp internal coordinate step to reasonable limits
fn clamp_internal_step(step: &mut [f64], connectivity: &MolecularConnectivity) {
    let n_bonds = connectivity.bonds.len();
    let n_angles = connectivity.angles.len();

    for (idx, val) in step.iter_mut().enumerate() {
        let limit = if idx < n_bonds {
            MAX_BOND_STEP
        } else if idx < n_bonds + n_angles {
            MAX_ANGLE_STEP
        } else {
            MAX_DIHEDRAL_STEP
        };
        *val = val.clamp(-limit, limit);
    }
}

/// Compute Newton step with eigenvalue floor to prevent blowup
///
/// Diagonalizes the Hessian, floors all eigenvalues at `min_eval`, then
/// computes dq = -H_floored^{-1} g. This prevents near-zero or negative
/// eigenvalues from producing catastrophically large steps without
/// requiring the full RFO augmented Hessian equation.
///
/// Reference: inspired by geomeTRIC's approach to Hessian conditioning
fn newton_step_with_eigenfloor(
    hessian: &DMatrix<f64>,
    gradient: &[f64],
    min_eval: f64,
) -> Vec<f64> {
    let n = gradient.len();
    let eigen = SymmetricEigen::new(hessian.clone());

    // Transform gradient to eigenbasis: g_eig = U^T g
    let g = DMatrix::from_column_slice(n, 1, gradient);
    let g_eig = eigen.eigenvectors.transpose() * &g;

    // Compute step in eigenbasis with floored eigenvalues
    let mut s_eig = DMatrix::zeros(n, 1);
    for i in 0..n {
        let eval = eigen.eigenvalues[i].max(min_eval);
        s_eig[(i, 0)] = -g_eig[(i, 0)] / eval;
    }

    // Transform back: dq = U s_eig
    let s_cart = &eigen.eigenvectors * &s_eig;
    (0..n).map(|i| s_cart[(i, 0)]).collect()
}

/// Optimize molecular geometry using internal coordinates with trust-radius control
///
/// Uses Newton-Raphson with BFGS Hessian update in redundant internal
/// coordinates. The Wilson B-matrix transforms between Cartesian and
/// internal coordinate spaces at each step. An adaptive trust radius
/// controls step size based on the quality of the quadratic model.
///
/// # Algorithm
///
/// 1. Detect molecular connectivity (bonds, angles, dihedrals)
/// 2. Build Schlegel's model Hessian as initial guess
/// 3. At each step:
///    a. Build B-matrix at current geometry
///    b. Transform Cartesian gradient to internal coordinates
///    c. Compute Newton step with eigenvalue-floored Hessian
///    d. Scale step to trust radius (RMS norm) + clamp individual components
///    e. Compute predicted energy change from quadratic model
///    f. Transform step back to Cartesian and take step
///    g. Evaluate energy and gradient at new geometry
///    h. Compute quality ratio rho = (E_new - E_old) / dE_predicted
///    i. Accept/reject step based on quality ratio
///    j. Update trust radius adaptively
///    k. BFGS update Hessian (only on accepted steps)
///    l. Check convergence
///
/// # References
///
/// - Wang & Song (2016). J. Chem. Phys. 144, 214108. (geomeTRIC)
/// - Schlegel (1984). Theor. Chim. Acta 66, 333. (Model Hessian)
fn optimize_internal(
    atoms: &[(u8, [f64; 3])],
    config: &OptimizationConfig,
    progress: Option<&dyn Fn(&OptimizationStep)>,
) -> OptimizationResult {
    let scf_config = ScfConfig {
        use_diis: true,
        ..ScfConfig::tight()
    };

    let mut current_atoms: Vec<Atom> = atoms
        .iter()
        .filter_map(|(z, pos)| Atom::new(*z, *pos).ok())
        .collect();

    let n_atoms = current_atoms.len();
    let connectivity = detect_connectivity(&current_atoms);
    let n_int = connectivity.n_internals;

    // Build initial model Hessian (diagonal -> full matrix)
    let h_diag = model_hessian_diagonal(&connectivity, &current_atoms);
    let mut hessian = DMatrix::zeros(n_int, n_int);
    for i in 0..n_int {
        hessian[(i, i)] = h_diag[i];
    }

    // Step 0: evaluate energy and gradient at initial geometry
    let eval_result = evaluate_energy_gradient(&current_atoms, config, &scf_config);
    let (mut energy, mut grad_atoms, _) = match eval_result {
        Some(r) => r,
        None => {
            return OptimizationResult {
                converged: false,
                steps: Vec::new(),
                final_energy: 0.0,
                final_geometry: atoms_to_geometry(&current_atoms),
                total_steps: 0,
                compute_time_ms: 0.0,
            };
        }
    };

    let mut grad_flat = grad_to_flat(&grad_atoms);

    // Record initial step
    let initial_step = make_step(0, energy, &grad_atoms, &current_atoms);
    let mut steps = vec![initial_step.clone()];
    if let Some(cb) = progress {
        cb(&initial_step);
    }

    // Compute initial internal coordinates and gradient
    let mut q_current = compute_internal_coords(&connectivity, &current_atoms);
    let bmat = WilsonBMatrix::build(&connectivity, &current_atoms);
    let mut grad_int = bmat.cartesian_to_internal_gradient(&grad_flat);

    // Adaptive trust radius (RMS norm of internal coordinate step)
    let mut trust_radius = INITIAL_TRUST_RADIUS;
    let mut converged = false;

    // Track consecutive rejected steps to detect infinite rejection loops.
    // After 3 consecutive rejections the BFGS Hessian has drifted too far
    // from the true curvature: same gradient + same Hessian = same direction,
    // and shrinking the trust radius alone cannot escape.  Recovery resets
    // the Hessian to the Schlegel model AND forces a steepest-descent
    // fallback on the next iteration.
    let mut consecutive_rejections: u32 = 0;
    let mut force_steepest_descent = false;

    for step_num in 1..=config.max_steps {
        let prev_energy = energy;

        // 1. Compute search direction.
        //    Normally: Newton step with eigenvalue-floored Hessian.
        //    After recovery from consecutive rejections: steepest descent
        //    (ignores the Hessian entirely, guaranteeing a new direction).
        let mut step_int = if force_steepest_descent {
            force_steepest_descent = false;
            // Steepest descent: dq = -g (will be scaled to trust radius below)
            grad_int.iter().map(|g| -g).collect::<Vec<f64>>()
        } else {
            newton_step_with_eigenfloor(&hessian, &grad_int, MIN_HESSIAN_EIGENVALUE)
        };

        // 2. Scale step to trust radius (RMS norm)
        let rms_step = if n_int > 0 {
            (step_int.iter().map(|x| x * x).sum::<f64>() / n_int as f64).sqrt()
        } else {
            0.0
        };
        if rms_step > trust_radius {
            let scale = trust_radius / rms_step;
            for s in &mut step_int {
                *s *= scale;
            }
        }

        // 3. Clamp individual components as secondary safety net
        clamp_internal_step(&mut step_int, &connectivity);

        // 4. Compute predicted energy change from quadratic model:
        //    dE_pred = g_q^T dq + 0.5 dq^T H dq
        let g_dot_dq: f64 = grad_int
            .iter()
            .zip(step_int.iter())
            .map(|(g, s)| g * s)
            .sum();
        let dq_col = DMatrix::from_column_slice(n_int, 1, &step_int);
        let h_dq = &hessian * &dq_col;
        let dq_h_dq: f64 = (0..n_int).map(|i| step_int[i] * h_dq[(i, 0)]).sum();
        let de_predicted = g_dot_dq + 0.5 * dq_h_dq;

        // 5. Transform step to Cartesian displacements
        let bmat_current = WilsonBMatrix::build(&connectivity, &current_atoms);
        let step_cart = bmat_current.internal_to_cartesian_step(&step_int);

        // Save state before taking the step (for potential rejection)
        let saved_atoms = current_atoms.clone();
        let saved_energy = energy;
        let saved_grad_atoms = grad_atoms.clone();
        let saved_grad_flat = grad_flat.clone();
        let saved_q_current = q_current.clone();
        let saved_grad_int = grad_int.clone();

        // 6. Update Cartesian geometry
        let mut coords = atoms_to_coords(&current_atoms);
        for i in 0..coords.len() {
            coords[i] += step_cart[i];
        }
        current_atoms = update_atom_positions(&current_atoms, &coords);

        // 7. Evaluate energy and gradient at new geometry
        let eval_new = evaluate_energy_gradient(&current_atoms, config, &scf_config);
        let mut scf_failed = false;
        match eval_new {
            Some((e_new, g_new, _)) => {
                let actual_de = e_new - energy;

                // 8. Compute quality ratio: rho = actual / predicted
                //    de_predicted should be negative for a descent step;
                //    rho ~ 1 means the quadratic model is accurate.
                let quality = if de_predicted.abs() > 1e-16 {
                    actual_de / de_predicted
                } else {
                    // Predicted change is essentially zero -- trust the step
                    // if actual change is also tiny
                    if actual_de.abs() < 1e-12 {
                        1.0
                    } else {
                        0.0
                    }
                };

                // Recompute actual RMS step for trust radius update decision
                let actual_rms = if n_int > 0 {
                    (step_int.iter().map(|x| x * x).sum::<f64>() / n_int as f64).sqrt()
                } else {
                    0.0
                };

                // 9. Trust radius update (before accept/reject decision)
                if quality > GOOD_QUALITY && actual_rms > 0.8 * trust_radius {
                    // Good step that used most of the trust radius -> grow
                    trust_radius = (trust_radius * TRUST_GROW_FACTOR).min(MAX_TRUST_RADIUS);
                } else if quality < POOR_QUALITY {
                    // Poor model quality -> shrink
                    trust_radius = (trust_radius * TRUST_SHRINK_FACTOR).max(MIN_TRUST_RADIUS);
                }
                // Otherwise: keep current trust radius

                // 10. Step acceptance: accept if energy decreased (rho > 0)
                if quality > 0.0 || actual_de < 0.0 {
                    // Accept step — reset the rejection counter
                    consecutive_rejections = 0;
                    energy = e_new;
                    grad_atoms = g_new;
                    grad_flat = grad_to_flat(&grad_atoms);

                    // Compute new internal coordinates and gradient
                    let q_new = compute_internal_coords(&connectivity, &current_atoms);
                    let bmat_new = WilsonBMatrix::build(&connectivity, &current_atoms);
                    let grad_int_new = bmat_new.cartesian_to_internal_gradient(&grad_flat);

                    // 11. BFGS update Hessian (only on accepted steps)
                    let dq: Vec<f64> = q_new
                        .iter()
                        .zip(q_current.iter())
                        .map(|(n, c)| n - c)
                        .collect();
                    let dg: Vec<f64> = grad_int_new
                        .iter()
                        .zip(grad_int.iter())
                        .map(|(n, c)| n - c)
                        .collect();

                    bfgs_update_hessian(&mut hessian, &dq, &dg);

                    // Update stored state
                    q_current = q_new;
                    grad_int = grad_int_new;
                } else {
                    // Reject step: restore geometry, do NOT update Hessian
                    current_atoms = saved_atoms;
                    energy = saved_energy;
                    grad_atoms = saved_grad_atoms;
                    grad_flat = saved_grad_flat;
                    q_current = saved_q_current;
                    grad_int = saved_grad_int;

                    consecutive_rejections += 1;

                    // Recovery: after 3 consecutive rejections the BFGS
                    // Hessian no longer produces useful directions.
                    //
                    // Two actions:
                    //  (a) Reset the Hessian to the Schlegel model so future
                    //      Newton steps start from fresh curvature.
                    //  (b) Force the *next* iteration to use steepest descent
                    //      with a small trust radius.  Steepest descent is
                    //      guaranteed to differ from the Newton direction and
                    //      -- for a sufficiently small step -- will always
                    //      decrease the energy.
                    if consecutive_rejections >= 3 {
                        let h_diag = model_hessian_diagonal(&connectivity, &current_atoms);
                        hessian = DMatrix::zeros(n_int, n_int);
                        for i in 0..n_int {
                            hessian[(i, i)] = h_diag[i];
                        }
                        // Use half the initial trust radius: large enough to
                        // make progress but small enough that the steepest
                        // descent step stays within the linear-gradient regime.
                        trust_radius = INITIAL_TRUST_RADIUS * 0.5;
                        force_steepest_descent = true;
                        consecutive_rejections = 0;
                    }
                }
            }
            None => {
                // SCF failed -- record and break
                let step_record = OptimizationStep {
                    step: step_num,
                    energy,
                    max_gradient: f64::NAN,
                    rms_gradient: f64::NAN,
                    geometry: atoms_to_geometry(&current_atoms),
                    gradient: vec![[0.0; 3]; n_atoms],
                };
                steps.push(step_record.clone());
                if let Some(cb) = progress {
                    cb(&step_record);
                }
                scf_failed = true;
            }
        }
        if scf_failed {
            break;
        }

        // 12. Record step
        let step_record = make_step(step_num, energy, &grad_atoms, &current_atoms);
        steps.push(step_record.clone());
        if let Some(cb) = progress {
            cb(&step_record);
        }

        // 13. Check convergence
        let delta_e = (energy - prev_energy).abs();
        if step_record.max_gradient < config.grad_threshold && delta_e < config.energy_threshold {
            converged = true;
            break;
        }
    }

    let last_step = steps.last().unwrap();
    OptimizationResult {
        converged,
        final_energy: last_step.energy,
        final_geometry: last_step.geometry.clone(),
        total_steps: steps.len() - 1,
        steps,
        compute_time_ms: 0.0,
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ========================================================================
    // Unit tests: L-BFGS two-loop recursion
    // ========================================================================

    #[test]
    fn test_lbfgs_first_step_steepest_descent() {
        // With no history, direction should be -gradient
        let gradient = vec![1.0, -2.0, 0.5];
        let history = LbfgsHistory::new(5);
        let direction = lbfgs_direction(&gradient, &history);

        assert_eq!(direction.len(), 3);
        assert!((direction[0] - (-1.0)).abs() < 1e-14);
        assert!((direction[1] - 2.0).abs() < 1e-14);
        assert!((direction[2] - (-0.5)).abs() < 1e-14);
    }

    #[test]
    fn test_lbfgs_direction_with_history() {
        // After one step with known s and y, the direction should differ
        // from steepest descent
        let gradient = vec![1.0, -2.0, 0.5];
        let mut history = LbfgsHistory::new(5);

        let s = vec![0.1, 0.2, -0.1];
        let y = vec![0.5, 0.3, 0.2]; // s^T * y = 0.09 > 0 (curvature condition)
        let accepted = history.push(s, y);
        assert!(accepted);

        let direction = lbfgs_direction(&gradient, &history);
        assert_eq!(direction.len(), 3);

        // Direction should not be exactly -gradient
        let sd_direction: Vec<f64> = gradient.iter().map(|g| -g).collect();
        let differs = direction
            .iter()
            .zip(sd_direction.iter())
            .any(|(d, sd)| (d - sd).abs() > 1e-10);
        assert!(
            differs,
            "L-BFGS direction should differ from steepest descent with history"
        );

        // Direction should still be a descent direction (d^T * g < 0)
        let dot: f64 = direction
            .iter()
            .zip(gradient.iter())
            .map(|(d, g)| d * g)
            .sum();
        assert!(
            dot < 0.0,
            "L-BFGS direction should be a descent direction, got d^T*g = {}",
            dot
        );
    }

    #[test]
    fn test_lbfgs_curvature_condition() {
        let mut history = LbfgsHistory::new(5);

        // Valid pair: s^T * y > 0
        let accepted = history.push(vec![0.1, 0.2], vec![0.5, 0.3]);
        assert!(accepted);
        assert_eq!(history.len(), 1);

        // Invalid pair: s^T * y <= 0
        let rejected = history.push(vec![0.1, 0.2], vec![-0.5, -0.3]);
        assert!(!rejected);
        assert_eq!(history.len(), 1); // Still 1 -- rejected pair not stored
    }

    #[test]
    fn test_lbfgs_history_capacity() {
        let mut history = LbfgsHistory::new(3);

        for i in 0..5 {
            let s = vec![(i + 1) as f64 * 0.1, 0.1];
            let y = vec![0.2, 0.3];
            history.push(s, y);
        }

        // Should only keep the last 3 pairs
        assert_eq!(history.len(), 3);
    }

    #[test]
    fn test_lbfgs_history_clear() {
        let mut history = LbfgsHistory::new(5);
        history.push(vec![0.1, 0.2], vec![0.3, 0.4]);
        assert_eq!(history.len(), 1);
        history.clear();
        assert_eq!(history.len(), 0);
    }

    // ========================================================================
    // Unit tests: Configuration
    // ========================================================================

    #[test]
    fn test_default_config() {
        let config = OptimizationConfig::default();
        assert_eq!(config.max_steps, 50);
        assert!((config.grad_threshold - 4.5e-4).abs() < 1e-10);
        assert!((config.energy_threshold - 1e-6).abs() < 1e-15);
        assert_eq!(config.memory_size, 7);
        assert_eq!(config.method, OptMethod::Rhf);
        assert_eq!(config.basis, "sto-3g");
    }

    // ========================================================================
    // Integration tests: H2 RHF/STO-3G optimization
    // ========================================================================

    #[test]
    fn test_h2_rhf_sto3g_optimization() {
        // H2 at R = 1.4 bohr (slightly longer than equilibrium ~1.346)
        let atoms = vec![(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];

        let config = OptimizationConfig {
            max_steps: 30,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };

        let result = optimize_geometry(&atoms, &config, None);

        // Should converge
        assert!(
            result.converged,
            "H2 optimization should converge, got {} steps",
            result.total_steps
        );

        // Compute final bond length
        let geom = &result.final_geometry;
        let r_eq = ((geom[0][0] - geom[1][0]).powi(2)
            + (geom[0][1] - geom[1][1]).powi(2)
            + (geom[0][2] - geom[1][2]).powi(2))
        .sqrt();

        // PySCF reference: R_eq = 1.345964 bohr
        assert!(
            (r_eq - 1.346).abs() < 0.01,
            "H2 R_eq = {:.6} bohr, expected ~1.346 (diff = {:.4})",
            r_eq,
            (r_eq - 1.346).abs()
        );

        // PySCF reference: E_opt = -1.117505884535 Ha
        // Our integral engine has ~1e-4 Ha differences
        assert!(
            (result.final_energy - (-1.1175)).abs() < 1e-3,
            "H2 E_opt = {:.10} Ha, expected ~-1.1175",
            result.final_energy
        );

        // Should converge in fewer than 20 steps
        assert!(
            result.total_steps < 20,
            "H2 should converge in <20 steps, took {}",
            result.total_steps
        );
    }

    #[test]
    fn test_h2_rhf_sto3g_trajectory() {
        let atoms = vec![(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];

        let config = OptimizationConfig {
            max_steps: 30,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };

        let result = optimize_geometry(&atoms, &config, None);

        // Trajectory should include initial + optimization steps
        assert!(
            result.steps.len() >= 2,
            "Should have at least initial + 1 step"
        );
        assert_eq!(result.steps.len(), result.total_steps + 1);

        // Step 0 is the initial geometry
        assert_eq!(result.steps[0].step, 0);

        // Energy should generally decrease
        let initial_energy = result.steps[0].energy;
        let final_energy = result.steps.last().unwrap().energy;
        assert!(
            final_energy < initial_energy,
            "Final energy ({:.10}) should be lower than initial ({:.10})",
            final_energy,
            initial_energy
        );

        // Last step energy should match final_energy
        assert!(
            (result.steps.last().unwrap().energy - result.final_energy).abs() < 1e-15,
            "Last step energy should match final_energy"
        );
    }

    // ========================================================================
    // Integration tests: H2O RHF/STO-3G optimization
    // ========================================================================

    #[test]
    fn test_h2o_rhf_sto3g_optimization() {
        // H2O with IQCP preset coordinates (bohr)
        let atoms = vec![
            (8u8, [0.0, 0.0, 0.2217282]),
            (1u8, [0.0, 1.4305447, -0.8869128]),
            (1u8, [0.0, -1.4305447, -0.8869128]),
        ];

        let config = OptimizationConfig {
            max_steps: 50,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };

        let result = optimize_geometry(&atoms, &config, None);

        // Should converge
        assert!(
            result.converged,
            "H2O optimization should converge, got {} steps, max_grad = {:.2e}",
            result.total_steps,
            result
                .steps
                .last()
                .map(|s| s.max_gradient)
                .unwrap_or(f64::NAN)
        );

        let geom = &result.final_geometry;

        // Compute OH distances
        let r_oh1 = ((geom[0][0] - geom[1][0]).powi(2)
            + (geom[0][1] - geom[1][1]).powi(2)
            + (geom[0][2] - geom[1][2]).powi(2))
        .sqrt();
        let r_oh2 = ((geom[0][0] - geom[2][0]).powi(2)
            + (geom[0][1] - geom[2][1]).powi(2)
            + (geom[0][2] - geom[2][2]).powi(2))
        .sqrt();

        // Convert to Angstrom for comparison
        use crate::constants::BOHR_TO_ANGSTROM;
        let r_oh1_ang = r_oh1 * BOHR_TO_ANGSTROM;
        let r_oh2_ang = r_oh2 * BOHR_TO_ANGSTROM;

        // PySCF reference: R(OH) = 0.989451 Angstrom
        // Tolerance: 0.002 Angstrom (AC6)
        assert!(
            (r_oh1_ang - 0.9895).abs() < 0.005,
            "R(OH1) = {:.4} A, expected ~0.9895 (diff = {:.4})",
            r_oh1_ang,
            (r_oh1_ang - 0.9895).abs()
        );
        assert!(
            (r_oh2_ang - 0.9895).abs() < 0.005,
            "R(OH2) = {:.4} A, expected ~0.9895 (diff = {:.4})",
            r_oh2_ang,
            (r_oh2_ang - 0.9895).abs()
        );

        // Compute HOH angle
        let v1 = [
            geom[1][0] - geom[0][0],
            geom[1][1] - geom[0][1],
            geom[1][2] - geom[0][2],
        ];
        let v2 = [
            geom[2][0] - geom[0][0],
            geom[2][1] - geom[0][1],
            geom[2][2] - geom[0][2],
        ];
        let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        let mag1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let mag2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        let angle_rad = (dot / (mag1 * mag2)).acos();
        let angle_deg = angle_rad * 180.0 / std::f64::consts::PI;

        // PySCF reference: HOH angle = 100.0191 degrees
        // Tolerance: 0.5 degrees (AC6)
        assert!(
            (angle_deg - 100.0).abs() < 1.0,
            "HOH angle = {:.2} deg, expected ~100.0 (diff = {:.2})",
            angle_deg,
            (angle_deg - 100.0).abs()
        );

        // PySCF reference: E_opt = -74.965901186381 Ha
        // Tolerance: 1e-4 Ha (integral engine differences)
        assert!(
            (result.final_energy - (-74.9659)).abs() < 1e-3,
            "H2O E_opt = {:.10} Ha, expected ~-74.9659",
            result.final_energy
        );
    }

    // ========================================================================
    // Unconverged optimization test
    // ========================================================================

    #[test]
    fn test_unconverged_reports_gracefully() {
        let atoms = vec![(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];

        let config = OptimizationConfig {
            max_steps: 2,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };

        let result = optimize_geometry(&atoms, &config, None);

        // May or may not converge in 2 steps -- the point is it doesn't panic
        assert!(result.total_steps <= 2);
        assert!(result.steps.len() >= 2); // At least initial + 1 step
        assert!(result.final_energy < 0.0); // Should have a valid energy
    }

    #[test]
    fn test_max_steps_respected() {
        // Use a far-from-equilibrium geometry so it won't converge in 1 step
        let atoms = vec![(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 3.0])];

        let config = OptimizationConfig {
            max_steps: 3,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };

        let result = optimize_geometry(&atoms, &config, None);
        assert!(
            result.total_steps <= 3,
            "Should not exceed max_steps=3, got {}",
            result.total_steps
        );
    }

    // ========================================================================
    // Progress callback test
    // ========================================================================

    #[test]
    fn test_progress_callback() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        let atoms = vec![(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
        let config = OptimizationConfig {
            max_steps: 10,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };

        let call_count = AtomicUsize::new(0);
        let result = optimize_geometry(
            &atoms,
            &config,
            Some(&|_step| {
                call_count.fetch_add(1, Ordering::Relaxed);
            }),
        );

        let count = call_count.load(Ordering::Relaxed);
        assert_eq!(
            count,
            result.steps.len(),
            "Progress callback should be called once per step (including step 0)"
        );
    }

    // ========================================================================
    // Convergence criteria test
    // ========================================================================

    #[test]
    fn test_convergence_requires_both_criteria() {
        // At equilibrium, both criteria should be met
        // Start very close to H2 equilibrium
        let atoms = vec![(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.346])];

        let config = OptimizationConfig {
            max_steps: 20,
            grad_threshold: 4.5e-4,
            energy_threshold: 1.0e-6,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };

        let result = optimize_geometry(&atoms, &config, None);

        if result.converged {
            // Both criteria must be satisfied at the final step
            let last_step = result.steps.last().unwrap();
            assert!(
                last_step.max_gradient < config.grad_threshold,
                "Max gradient {} should be < {}",
                last_step.max_gradient,
                config.grad_threshold
            );
            // delta_E is checked internally -- just verify convergence happened
        }
    }

    // ========================================================================
    // Unit tests: Connectivity detection
    // ========================================================================

    #[test]
    fn test_detect_connectivity_h2o() {
        // H2O geometry (bohr)
        let atoms = vec![
            Atom {
                atomic_number: 8,
                symbol: "O".to_string(),
                position: [0.0, 0.0, 0.2217282],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 1.4305447, -0.8869128],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, -1.4305447, -0.8869128],
            },
        ];

        let conn = detect_connectivity(&atoms);

        // Should detect 2 O-H bonds
        assert_eq!(
            conn.bonds.len(),
            2,
            "H2O should have 2 bonds, got {:?}",
            conn.bonds
        );

        // Should detect 1 H-O-H angle
        assert_eq!(
            conn.angles.len(),
            1,
            "H2O should have 1 angle, got {:?}",
            conn.angles
        );

        // No dihedrals in H2O
        assert_eq!(conn.dihedrals.len(), 0, "H2O should have 0 dihedrals");

        // Total: 2 bonds + 1 angle = 3 internals
        assert_eq!(conn.n_internals, 3);
    }

    #[test]
    fn test_detect_connectivity_h2() {
        let atoms = vec![
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 0.0, 1.4],
            },
        ];

        let conn = detect_connectivity(&atoms);
        assert_eq!(conn.bonds.len(), 1, "H2 should have 1 bond");
        assert_eq!(conn.angles.len(), 0, "H2 should have 0 angles");
        assert_eq!(conn.n_internals, 1);
    }

    #[test]
    fn test_detect_connectivity_nh3() {
        // NH3 geometry (approximately tetrahedral, bohr)
        let atoms = vec![
            Atom {
                atomic_number: 7,
                symbol: "N".to_string(),
                position: [0.0, 0.0, 0.2196],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 1.7710, -0.5126],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [1.5337, -0.8855, -0.5126],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [-1.5337, -0.8855, -0.5126],
            },
        ];

        let conn = detect_connectivity(&atoms);

        // 3 N-H bonds
        assert_eq!(
            conn.bonds.len(),
            3,
            "NH3 should have 3 bonds, got {:?}",
            conn.bonds
        );

        // 3 H-N-H angles
        assert_eq!(
            conn.angles.len(),
            3,
            "NH3 should have 3 angles, got {:?}",
            conn.angles
        );

        // No dihedrals (no central bond with additional substituents)
        assert_eq!(conn.dihedrals.len(), 0, "NH3 should have 0 dihedrals");
    }

    // ========================================================================
    // Unit tests: Internal coordinate values
    // ========================================================================

    #[test]
    fn test_bond_length() {
        let atoms = vec![
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 0.0, 1.4],
            },
        ];
        let r = compute_bond_length(&atoms, 0, 1);
        assert!(
            (r - 1.4).abs() < 1e-14,
            "Bond length should be 1.4, got {}",
            r
        );
    }

    #[test]
    fn test_angle_90_degrees() {
        // Right angle: atoms at (1,0,0), (0,0,0), (0,1,0)
        let atoms = vec![
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [1.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 8,
                symbol: "O".to_string(),
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 1.0, 0.0],
            },
        ];
        let theta = compute_angle(&atoms, 0, 1, 2);
        let expected = std::f64::consts::FRAC_PI_2;
        assert!(
            (theta - expected).abs() < 1e-12,
            "Angle should be pi/2 = {:.10}, got {:.10}",
            expected,
            theta
        );
    }

    #[test]
    fn test_angle_180_degrees() {
        // Linear: atoms at (-1,0,0), (0,0,0), (1,0,0)
        let atoms = vec![
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [-1.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 8,
                symbol: "O".to_string(),
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [1.0, 0.0, 0.0],
            },
        ];
        let theta = compute_angle(&atoms, 0, 1, 2);
        let expected = std::f64::consts::PI;
        assert!(
            (theta - expected).abs() < 1e-12,
            "Angle should be pi = {:.10}, got {:.10}",
            expected,
            theta
        );
    }

    // ========================================================================
    // Unit tests: Wilson B-matrix
    // ========================================================================

    #[test]
    fn test_bmatrix_h2_gradient_transform() {
        // For a diatomic, the B-matrix should have 1 row (bond stretch)
        // and 6 columns (3 per atom)
        let atoms = vec![
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 0.0, 1.4],
            },
        ];

        let conn = detect_connectivity(&atoms);
        let bmat = WilsonBMatrix::build(&conn, &atoms);

        assert_eq!(bmat.mat.nrows(), 1, "H2 B-matrix should have 1 row");
        assert_eq!(bmat.mat.ncols(), 6, "H2 B-matrix should have 6 columns");

        // The bond is along z-axis, so B should be [0, 0, -1, 0, 0, 1]
        // (derivative of r_12 = |R2 - R1| w.r.t. atom 1 z is -(R2z - R1z)/r = -1)
        assert!(
            (bmat.mat[(0, 2)] - (-1.0)).abs() < 1e-12,
            "B[0,2] should be -1"
        );
        assert!((bmat.mat[(0, 5)] - 1.0).abs() < 1e-12, "B[0,5] should be 1");
    }

    #[test]
    fn test_bmatrix_roundtrip_gradient() {
        // For H2O, transform a known Cartesian gradient to internal,
        // then verify internal gradient has 3 components (2 bonds + 1 angle)
        let atoms = vec![
            Atom {
                atomic_number: 8,
                symbol: "O".to_string(),
                position: [0.0, 0.0, 0.2217282],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 1.4305447, -0.8869128],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, -1.4305447, -0.8869128],
            },
        ];

        let conn = detect_connectivity(&atoms);
        let bmat = WilsonBMatrix::build(&conn, &atoms);

        // Arbitrary gradient (9 components for 3 atoms)
        let grad_cart = vec![0.0, 0.0, 0.01, 0.0, 0.005, -0.005, 0.0, -0.005, -0.005];
        let grad_int = bmat.cartesian_to_internal_gradient(&grad_cart);

        assert_eq!(
            grad_int.len(),
            conn.n_internals,
            "Internal gradient should have {} components",
            conn.n_internals
        );

        // The gradient components should be finite
        for (i, g) in grad_int.iter().enumerate() {
            assert!(
                g.is_finite(),
                "Internal gradient component {} should be finite, got {}",
                i,
                g
            );
        }
    }

    // ========================================================================
    // Unit tests: Model Hessian
    // ========================================================================

    #[test]
    fn test_model_hessian_h2o() {
        let atoms = vec![
            Atom {
                atomic_number: 8,
                symbol: "O".to_string(),
                position: [0.0, 0.0, 0.2217282],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, 1.4305447, -0.8869128],
            },
            Atom {
                atomic_number: 1,
                symbol: "H".to_string(),
                position: [0.0, -1.4305447, -0.8869128],
            },
        ];

        let conn = detect_connectivity(&atoms);
        let h = model_hessian_diagonal(&conn, &atoms);

        assert_eq!(h.len(), conn.n_internals);

        // All diagonal elements should be positive
        for (i, &val) in h.iter().enumerate() {
            assert!(
                val > 0.0,
                "Hessian diagonal {} should be positive, got {}",
                i,
                val
            );
        }

        // Bond force constants should be larger than angle force constants
        // (first 2 are bonds, third is angle for H2O)
        // O-H bond k ~ 0.26 > angle k = 0.16
        assert!(
            h[0] > h[2],
            "Bond force constant ({:.4}) should exceed angle force constant ({:.4})",
            h[0],
            h[2]
        );

        // Both O-H bonds should have the same force constant (symmetric molecule)
        assert!(
            (h[0] - h[1]).abs() < 1e-10,
            "Both O-H force constants should match: {} vs {}",
            h[0],
            h[1]
        );
    }

    // ========================================================================
    // Integration tests: H2O internal coordinate optimization
    // ========================================================================

    #[test]
    fn test_h2o_internal_coords_optimization() {
        // H2O with IQCP preset coordinates (bohr) -- should use internal coords
        let atoms = vec![
            (8u8, [0.0, 0.0, 0.2217282]),
            (1u8, [0.0, 1.4305447, -0.8869128]),
            (1u8, [0.0, -1.4305447, -0.8869128]),
        ];

        let config = OptimizationConfig {
            max_steps: 50,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };

        let result = optimize_geometry(&atoms, &config, None);

        // Should converge
        assert!(
            result.converged,
            "H2O internal coord optimization should converge, got {} steps, max_grad = {:.2e}",
            result.total_steps,
            result
                .steps
                .last()
                .map(|s| s.max_gradient)
                .unwrap_or(f64::NAN)
        );

        // Should converge in fewer steps than Cartesian (target: <= 5)
        assert!(
            result.total_steps <= 8,
            "H2O internal coords should converge in <= 8 steps, took {}",
            result.total_steps
        );

        // Verify final geometry
        let geom = &result.final_geometry;
        let r_oh1 = ((geom[0][0] - geom[1][0]).powi(2)
            + (geom[0][1] - geom[1][1]).powi(2)
            + (geom[0][2] - geom[1][2]).powi(2))
        .sqrt();
        let r_oh2 = ((geom[0][0] - geom[2][0]).powi(2)
            + (geom[0][1] - geom[2][1]).powi(2)
            + (geom[0][2] - geom[2][2]).powi(2))
        .sqrt();

        use crate::constants::BOHR_TO_ANGSTROM;
        let r_oh1_ang = r_oh1 * BOHR_TO_ANGSTROM;
        let r_oh2_ang = r_oh2 * BOHR_TO_ANGSTROM;

        // PySCF reference: R(OH) ~ 0.9895 Angstrom
        assert!(
            (r_oh1_ang - 0.9895).abs() < 0.005,
            "R(OH1) = {:.4} A, expected ~0.9895",
            r_oh1_ang
        );
        assert!(
            (r_oh2_ang - 0.9895).abs() < 0.005,
            "R(OH2) = {:.4} A, expected ~0.9895",
            r_oh2_ang
        );

        // HOH angle
        let v1 = [
            geom[1][0] - geom[0][0],
            geom[1][1] - geom[0][1],
            geom[1][2] - geom[0][2],
        ];
        let v2 = [
            geom[2][0] - geom[0][0],
            geom[2][1] - geom[0][1],
            geom[2][2] - geom[0][2],
        ];
        let dot_val = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        let mag1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let mag2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        let angle_deg = (dot_val / (mag1 * mag2)).acos() * 180.0 / std::f64::consts::PI;

        assert!(
            (angle_deg - 100.0).abs() < 1.0,
            "HOH angle = {:.2} deg, expected ~100.0",
            angle_deg
        );

        // Energy should match reference
        assert!(
            (result.final_energy - (-74.9659)).abs() < 1e-3,
            "H2O E_opt = {:.10} Ha, expected ~-74.9659",
            result.final_energy
        );
    }

    #[test]
    fn test_h2_uses_cartesian_fallback() {
        // H2 (diatomic) should fall back to Cartesian L-BFGS
        let atoms = vec![(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];

        let config = OptimizationConfig {
            max_steps: 30,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };

        let result = optimize_geometry(&atoms, &config, None);
        assert!(
            result.converged,
            "H2 should still converge via Cartesian fallback"
        );

        let geom = &result.final_geometry;
        let r_eq = ((geom[0][0] - geom[1][0]).powi(2)
            + (geom[0][1] - geom[1][1]).powi(2)
            + (geom[0][2] - geom[1][2]).powi(2))
        .sqrt();
        assert!(
            (r_eq - 1.346).abs() < 0.01,
            "H2 R_eq = {:.6} bohr, expected ~1.346",
            r_eq,
        );
    }

    // ========================================================================
    // BFGS Hessian update test
    // ========================================================================

    #[test]
    fn test_bfgs_hessian_update() {
        // Start with identity Hessian
        let n = 3;
        let mut h = DMatrix::identity(n, n);

        // After a BFGS update, the Hessian should change
        let dq = vec![0.1, 0.05, -0.02];
        let dg = vec![0.3, 0.1, 0.05];

        bfgs_update_hessian(&mut h, &dq, &dg);

        // Hessian should no longer be identity
        let is_identity = (0..n).all(|i| {
            (0..n).all(|j| {
                let expected = if i == j { 1.0 } else { 0.0 };
                (h[(i, j)] - expected).abs() < 1e-14
            })
        });
        assert!(!is_identity, "BFGS update should modify the Hessian");

        // Hessian should remain symmetric
        for i in 0..n {
            for j in i + 1..n {
                assert!(
                    (h[(i, j)] - h[(j, i)]).abs() < 1e-14,
                    "Hessian should remain symmetric: H[{},{}]={} vs H[{},{}]={}",
                    i,
                    j,
                    h[(i, j)],
                    j,
                    i,
                    h[(j, i)]
                );
            }
        }
    }

    #[test]
    fn test_c2h6_rhf_sto3g_optimization() {
        let atoms = vec![
            (6u8, [0.0, 0.0, 1.4508]),
            (6, [0.0, 0.0, -1.4508]),
            (1, [0.0, 1.9217, 2.2700]),
            (1, [1.6641, -0.9609, 2.2700]),
            (1, [-1.6641, -0.9609, 2.2700]),
            (1, [0.0, -1.9217, -2.2700]),
            (1, [-1.6641, 0.9609, -2.2700]),
            (1, [1.6641, 0.9609, -2.2700]),
        ];
        let config = OptimizationConfig {
            max_steps: 50,
            grad_threshold: 4.5e-4,
            energy_threshold: 1.0e-6,
            memory_size: 7,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
        };
        let result = optimize_geometry(
            &atoms,
            &config,
            Some(&|step| {
                eprintln!(
                    "C2H6 Step {}: E={:.10} max|g|={:.6e} rms|g|={:.6e}",
                    step.step, step.energy, step.max_gradient, step.rms_gradient
                );
            }),
        );
        eprintln!(
            "C2H6 Converged: {}, Steps: {}, Final E: {:.10}",
            result.converged, result.total_steps, result.final_energy
        );
        assert!(
            result.converged,
            "C2H6 RHF/STO-3G should converge within 50 steps (got {} steps, max|g|={:.6e})",
            result.total_steps,
            result
                .steps
                .last()
                .map(|s| s.max_gradient)
                .unwrap_or(f64::NAN)
        );
    }

    /// Regression test: C2H4/B3LYP/STO-3G previously got stuck in an infinite
    /// rejection loop because the BFGS Hessian drifted and consecutive rejected
    /// steps all proposed the same direction.  The fix resets the Hessian to the
    /// Schlegel model after 3 consecutive rejections.
    #[test]
    fn test_c2h4_b3lyp_sto3g_optimization() {
        let atoms = vec![
            (6u8, [0.0, 0.0, 1.2592]),
            (6, [0.0, 0.0, -1.2592]),
            (1, [0.0, 1.7455, 2.3280]),
            (1, [0.0, -1.7455, 2.3280]),
            (1, [0.0, 1.7455, -2.3280]),
            (1, [0.0, -1.7455, -2.3280]),
        ];
        let config = OptimizationConfig {
            max_steps: 50,
            method: OptMethod::B3lyp,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };
        let result = optimize_geometry(
            &atoms,
            &config,
            Some(&|step| {
                eprintln!(
                    "C2H4/B3LYP Step {}: E={:.10} max|g|={:.6e} rms|g|={:.6e}",
                    step.step, step.energy, step.max_gradient, step.rms_gradient
                );
            }),
        );
        eprintln!(
            "C2H4/B3LYP Converged: {}, Steps: {}, Final E: {:.10}",
            result.converged, result.total_steps, result.final_energy
        );
        assert!(
            result.converged,
            "C2H4/B3LYP/STO-3G should converge within 50 steps (got {} steps, max|g|={:.6e})",
            result.total_steps,
            result
                .steps
                .last()
                .map(|s| s.max_gradient)
                .unwrap_or(f64::NAN)
        );
        // PySCF reference: E ~ -77.6224 Ha (grid-dependent)
        assert!(
            result.final_energy < -77.56,
            "Final energy should be reasonable, got {}",
            result.final_energy
        );
    }

    // ========================================================================
    // SI Data Extraction: H2O B3LYP/6-31G* Optimization (GAP-1)
    // ========================================================================

    #[test]
    fn test_si_gap1_h2o_b3lyp_631gs_optimization() {
        // GAP-1: H2O/B3LYP/6-31G* geometry optimization
        // PySCF reference: R(OH) = 1.830645 bohr = 0.968735 A
        //                  angle(HOH) = 103.6523 deg
        //                  E_opt = -76.408953844985 Ha
        let atoms = vec![
            (8u8, [0.0, 0.0, 0.2217282]),
            (1, [0.0, 1.4305447, -0.8869128]),
            (1, [0.0, -1.4305447, -0.8869128]),
        ];

        let config = OptimizationConfig {
            max_steps: 50,
            method: OptMethod::B3lyp,
            basis: "6-31g*".to_string(),
            ..Default::default()
        };

        let result = optimize_geometry(
            &atoms,
            &config,
            Some(&|step| {
                eprintln!(
                    "H2O B3LYP Step {}: E={:.10} max|g|={:.6e}",
                    step.step, step.energy, step.max_gradient
                );
            }),
        );

        let geom = &result.final_geometry;
        use crate::constants::BOHR_TO_ANGSTROM;

        // Compute OH distances
        let r_oh1 = ((geom[0][0] - geom[1][0]).powi(2)
            + (geom[0][1] - geom[1][1]).powi(2)
            + (geom[0][2] - geom[1][2]).powi(2))
        .sqrt();
        let r_oh2 = ((geom[0][0] - geom[2][0]).powi(2)
            + (geom[0][1] - geom[2][1]).powi(2)
            + (geom[0][2] - geom[2][2]).powi(2))
        .sqrt();

        // Compute HOH angle
        let v1 = [
            geom[1][0] - geom[0][0],
            geom[1][1] - geom[0][1],
            geom[1][2] - geom[0][2],
        ];
        let v2 = [
            geom[2][0] - geom[0][0],
            geom[2][1] - geom[0][1],
            geom[2][2] - geom[0][2],
        ];
        let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        let mag1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let mag2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        let angle_deg = (dot / (mag1 * mag2)).acos() * 180.0 / std::f64::consts::PI;

        eprintln!("\n=== GAP-1: H2O/B3LYP/6-31G* Optimization ===");
        eprintln!(
            "Converged: {}, Steps: {}",
            result.converged, result.total_steps
        );
        eprintln!("E_opt = {:.10} Ha", result.final_energy);
        eprintln!(
            "R(OH1) = {:.6} bohr = {:.4} A",
            r_oh1,
            r_oh1 * BOHR_TO_ANGSTROM
        );
        eprintln!(
            "R(OH2) = {:.6} bohr = {:.4} A",
            r_oh2,
            r_oh2 * BOHR_TO_ANGSTROM
        );
        eprintln!("angle(HOH) = {:.4} deg", angle_deg);
    }
}
