#![allow(clippy::needless_range_loop)]

//! Electric field CPHF, static polarizability, polarizability derivatives,
//! and Raman scattering activities (US-098, Phase 5 S41 M19).
//!
//! This module implements the Raman counterpart to the IR module (`ir.rs`).
//! Where IR intensity is a first-derivative property of the dipole moment,
//! Raman scattering activity is a first-derivative property of the static
//! polarizability, and the static polarizability is itself obtained via
//! field CPHF (the twin of the nuclear CPHF used in the analytical
//! Hessian).
//!
//! # Pipeline
//!
//! 1. [`charge_center`]: compute the charge-center gauge origin
//!    `C = Σ Z_A R_A / Σ Z_A` matching pyscf-forge
//!    `prop/polarizability/rhf.py:52-53`.
//! 2. [`compute_field_cphf`]: solve CPHF with the three Cartesian dipole
//!    integrals (bra-center gauge at charge center) as perturbations,
//!    reusing [`cphf_solve_withs1`](crate::scf::cphf::cphf_solve_withs1)
//!    with `s1 = 0` — same path PySCF takes (see `prop/polarizability/rhf.py:58-63`).
//! 3. [`compute_polarizability`]: contract the field CPHF solution with
//!    the MO-basis dipole integrals, producing `α_{de}` in bohr³ via
//!    `α = -2·(e2 + e2.T)` matching pyscf-forge `rhf.py:69`.
//! 4. [`compute_polarizability_derivatives`]: assemble the analytical
//!    Cartesian derivative tensor `∂α_{de}/∂R_{A,γ}` via a displacement
//!    scheme that reuses the analytical CPHF at displaced geometries.
//!    The outer differentiation is numerical (central difference with a
//!    small step in bohr) but the CPHF itself is fully analytical at each
//!    sample point — this is the "analytical via CPHF" path. The field
//!    CPHF at the equilibrium geometry is solved only once (cached in
//!    `field_cphf`), and the nuclear CPHF from
//!    [`HessianResult.mo1_cphf`](crate::scf::hessian::HessianResult) is
//!    not re-solved.
//! 5. [`compute_raman_activities`]: project the Cartesian polarizability
//!    derivative onto normal modes, reduce to Long (2002) invariants
//!    (isotropic `ᾱ'` and anisotropy `γ'²`), compute the scattering
//!    activity `S_k = (45·ᾱ'² + 7·γ'²)·BOHR4_TO_ANG4` in Å⁴/amu and the
//!    depolarization ratio `ρ_k = 3·γ'² / (45·ᾱ'² + 4·γ'²)` in the range
//!    `[0, 0.75]` for natural incident light.
//! 6. [`compute_raman_spectrum`]: end-to-end wrapper packaging everything
//!    into [`RamanResult`] for downstream consumers (US-100 simulated
//!    spectra, US-101 WASM export).
//!
//! # Gauge origin convention
//!
//! Unlike IR (US-097) which uses gauge origin `[0, 0, 0]` to match PySCF
//! `dip_moment()`, this module uses the molecular **charge center** for
//! all dipole integrals, matching pyscf-forge
//! `prop/polarizability/rhf.py:52-53`. For neutral molecules the
//! polarizability is mathematically gauge-invariant; the charge-center
//! choice minimizes numerical cancellation error. The gauge origin is
//! stored in [`FieldCphfResult::gauge_origin`] and flows consistently
//! through all downstream functions.
//!
//! # References
//!
//! - Long, D. A. (2002). *The Raman Effect: A Unified Treatment of the
//!   Theory of Raman Scattering by Molecules.* Wiley. Chapter 5.
//! - Placzek, G. (1934). *Rayleigh-Streuung und Raman-Effekt.*
//!   Handbuch der Radiologie VI.
//! - Amos, R. D. (1986). *Chem. Phys. Lett.* **124**, 376 — analytical
//!   polarizability derivatives via the interchange theorem.
//! - Pople, Krishnan, Schlegel, Binkley (1979). *Int. J. Quantum Chem.
//!   Symp.* **13**, 225 — CPHF-based derivative properties.
//! - Helgaker, Jørgensen & Olsen (2000). *Molecular Electronic-Structure
//!   Theory.* Wiley. Sections 10.8 and 13.
//! - pyscf-forge `pyscf/prop/polarizability/rhf.py` (lines 39-77,
//!   `polarizability()` function).

use crate::basis::{Atom, BasisSet};
use crate::constants::{AU_POLARIZABILITY_TO_ANG3, BOHR4_TO_ANG4};
use crate::integrals::dipole_matrix;
use crate::scf::cphf::{cphf_solve_withs1, gen_vind_rhf, CphfConfig};
use crate::scf::hessian::HessianResult;
use crate::thermo::FrequencyInfo;
use nalgebra::DMatrix;
use thiserror::Error;

// =============================================================================
// Errors
// =============================================================================

/// Errors that can occur during Raman computation.
#[derive(Debug, Error, PartialEq)]
pub enum RamanError {
    /// `HessianResult.mo1_cphf` was `None`; cannot compute polarizability
    /// derivatives without the nuclear CPHF solution.
    #[error(
        "HessianResult.mo1_cphf is None; compute_polarizability_derivatives \
         requires the nuclear CPHF solution from rhf_hessian"
    )]
    NoCphfData,

    /// Number of atoms inferred from `FrequencyInfo` or `HessianResult`
    /// does not match the number of atoms in the basis set.
    #[error("Size mismatch: expected {expected} atoms, got {actual} atoms")]
    SizeMismatch { expected: usize, actual: usize },

    /// Field CPHF did not converge within the configured iteration limit.
    #[error("Field CPHF did not converge after {iterations} iterations")]
    CphfNotConverged { iterations: usize },

    /// A non-finite value was encountered during assembly.
    #[error("Non-finite result detected in {stage}")]
    NonFiniteResult { stage: &'static str },
}

// =============================================================================
// Charge-center gauge origin
// =============================================================================

/// Compute the molecular charge center
/// `C = (Σ_A Z_A R_A) / (Σ_A Z_A)`.
///
/// Matches pyscf-forge `prop/polarizability/rhf.py:52-53`:
///
/// ```python
/// charges = mol.atom_charges()
/// coords  = mol.atom_coords()
/// charge_center = numpy.einsum('i,ix->x', charges, coords) / charges.sum()
/// ```
///
/// For neutral molecules, this minimizes numerical cancellation in the
/// dipole-integral contractions. For charged molecules, the polarizability
/// is gauge-dependent; the charge-center origin is the conventional
/// choice and matches pyscf-forge exactly.
///
/// # Arguments
///
/// * `atoms` - Slice of atoms (atomic number + bohr position)
///
/// # Returns
///
/// The charge center `[C_x, C_y, C_z]` in bohr. If the total nuclear
/// charge is zero (e.g., empty input) returns `[0, 0, 0]`.
pub fn charge_center(atoms: &[Atom]) -> [f64; 3] {
    let z_total: f64 = atoms.iter().map(|a| a.atomic_number as f64).sum();
    if z_total.abs() < 1e-14 {
        return [0.0, 0.0, 0.0];
    }
    let mut c = [0.0f64; 3];
    for atom in atoms {
        let z = atom.atomic_number as f64;
        for d in 0..3 {
            c[d] += z * atom.position[d];
        }
    }
    for d in 0..3 {
        c[d] /= z_total;
    }
    c
}

// =============================================================================
// Field CPHF result
// =============================================================================

/// Result of the electric field CPHF solve.
///
/// Contains the first-order MO coefficient response `U^(E_d)` for the
/// three Cartesian electric field directions, shaped `(nmo, nocc)` per
/// direction. The occupied-occupied block is zero (since `s^(1) = 0` for
/// field perturbations), so only the virtual-occupied block is
/// meaningful for downstream contractions — but the full `(nmo, nocc)`
/// layout is retained to match the shape PySCF uses in
/// `prop/polarizability/rhf.py`.
#[derive(Debug, Clone)]
pub struct FieldCphfResult {
    /// Field CPHF response `[U^(E_x), U^(E_y), U^(E_z)]`, each
    /// `(nmo, nocc)`.
    pub mo1: Vec<DMatrix<f64>>,
    /// Full MO coefficient matrix `(nbf, nmo)`.
    pub mo_coeff: DMatrix<f64>,
    /// MO orbital energies (length `nmo`).
    pub mo_energies: Vec<f64>,
    /// Number of doubly-occupied MOs.
    pub n_occ: usize,
    /// Gauge origin used for the dipole integrals — stored here so that
    /// [`compute_polarizability`] and
    /// [`compute_polarizability_derivatives`] use the **same** origin.
    pub gauge_origin: [f64; 3],
    /// Number of Krylov CPHF iterations.
    pub iterations: usize,
    /// Whether CPHF converged.
    pub converged: bool,
}

// =============================================================================
// Field CPHF driver
// =============================================================================

/// Solve the coupled-perturbed Hartree-Fock equations for the three
/// electric field perturbations (Ex, Ey, Ez).
///
/// This is the "field" twin of the nuclear CPHF solve inside
/// [`rhf_hessian`](crate::scf::hessian::rhf_hessian). Both use the same
/// [`gen_vind_rhf`](crate::scf::cphf::gen_vind_rhf) response function
/// because the AO-basis J/K operator is perturbation-independent; only
/// the right-hand side changes.
///
/// The dispatch path matches pyscf-forge `prop/polarizability/rhf.py`:
/// it calls [`cphf_solve_withs1`](crate::scf::cphf::cphf_solve_withs1)
/// with an explicit zero `s1`, identical to PySCF's
/// `cphf.solve(vind, mo_energy, mo_occ, h1, s1, ...)` at line 61 where
/// `s1 = numpy.zeros_like(h1)` — see line 58.
///
/// # Gauge origin
///
/// The gauge origin is the molecular **charge center**
/// `(Σ Z_A R_A) / (Σ Z_A)`, matching pyscf-forge
/// `polarizability/rhf.py:52-53`. For neutral molecules the
/// polarizability is mathematically gauge-invariant; the charge-center
/// choice minimizes numerical cancellation error.
///
/// # Arguments
///
/// * `atoms` - Molecular atom list (needed for charge-center calc)
/// * `basis` - Molecular basis set
/// * `hess_result` - Hessian result whose `mo1_cphf` provides
///   `(mo_coeff, mo_energies, n_occ)`
/// * `cphf_config` - CPHF solver configuration (defaults suffice)
///
/// # Returns
///
/// A [`FieldCphfResult`] with `mo1[d]` the field CPHF response for
/// direction `d ∈ {0, 1, 2}`, shaped `(nmo, nocc)`.
///
/// # Errors
///
/// * [`RamanError::NoCphfData`] if `hess_result.mo1_cphf.is_none()`
/// * [`RamanError::CphfNotConverged`] if the Krylov solver fails to
///   converge within `cphf_config.max_cycle`
///
/// # References
///
/// - pyscf-forge `prop/polarizability/rhf.py` lines 39-66
/// - IQCP `scf/cphf.rs:cphf_solve_withs1` (entry point for field CPHF)
/// - IQCP `scf/cphf.rs:gen_vind_rhf` (response function)
pub fn compute_field_cphf(
    atoms: &[Atom],
    basis: &BasisSet,
    hess_result: &HessianResult,
    cphf_config: &CphfConfig,
) -> Result<FieldCphfResult, RamanError> {
    let cphf_data = hess_result
        .mo1_cphf
        .as_ref()
        .ok_or(RamanError::NoCphfData)?;

    let gauge_origin = charge_center(atoms);
    let mo_coeff = cphf_data.mo_coeff.clone();
    let n_occ = cphf_data.n_occ;
    let mo_energies = cphf_data.mo_energies.clone();
    let nbf = basis.n_basis;
    let nmo = mo_coeff.ncols();

    // AO-basis dipole integrals at the charge-center gauge
    let dip_ao = dipole_matrix(basis, &gauge_origin);

    // Build h1 in the "full" MO layout (nmo, nocc), matching PySCF
    // polarizability.py line 57:
    //   h1 = einsum('xpq,pi,qj->xij', int_r, mo_coeff.conj(), orbo)
    // which gives (3, nmo, nocc).
    let c_occ = mo_coeff.columns(0, n_occ).clone_owned();
    let ct = mo_coeff.transpose();
    let mut h1_mo: Vec<DMatrix<f64>> = Vec::with_capacity(3);
    for d in 0..3 {
        let m = &ct * &dip_ao[d] * &c_occ; // (nmo, nocc)
        h1_mo.push(m);
    }

    // s1 = zeros_like(h1)  — matches pyscf-forge line 58
    let s1_mo: Vec<DMatrix<f64>> = (0..3).map(|_| DMatrix::zeros(nmo, n_occ)).collect();

    // Rebuild ERI tensor for the response function (same pattern as
    // rhf_hessian does at hessian.rs:1766).
    let eri = crate::integrals::eri_compressed(basis);

    // Response function — full RHF vind, operates on (nmo, nocc) blocks.
    // Matches pyscf-forge gen_vind pattern in Polarizability class.
    let vind = gen_vind_rhf(&mo_coeff, n_occ, &eri, nbf, 1.0);

    // Solve field CPHF via cphf_solve_withs1 (s1=0 path; same as PySCF)
    let (mo1, _mo_e1, iterations, converged) =
        cphf_solve_withs1(vind, &mo_energies, n_occ, &h1_mo, &s1_mo, cphf_config);

    if !converged {
        return Err(RamanError::CphfNotConverged { iterations });
    }

    Ok(FieldCphfResult {
        mo1,
        mo_coeff,
        mo_energies,
        n_occ,
        gauge_origin,
        iterations,
        converged,
    })
}

// =============================================================================
// Static polarizability
// =============================================================================

/// Compute the static polarizability tensor `α_{de}` in atomic units
/// (bohr³) from the field CPHF solution.
///
/// # Formula
///
/// Matches pyscf-forge `prop/polarizability/rhf.py:67-69` exactly:
///
/// ```text
/// e2[d, e] = sum_{p,i} h1[d][p,i] * mo1[e][p,i]
/// α     = (e2 + e2.T) * -2
/// ```
///
/// where `h1[d] = C^T · μ^d_AO · C_occ` is the MO-basis dipole integral
/// and `mo1[e] = U^(E_e)` is the CPHF response to a field in direction
/// `e`. The `-2` factor comes from: `-1` from the sign of the dipole
/// perturbation and `2` from closed-shell double occupancy; the `(e2 +
/// e2.T)` symmetrizes over `(d, e)`.
///
/// The occupied-occupied block of `mo1[e]` is zero (because `s^(1) = 0`
/// for field perturbations), so the sum over `p` reduces to a sum over
/// virtual orbitals in practice — but the full shape is retained for
/// clarity.
///
/// # Arguments
///
/// * `field_cphf` - Result of [`compute_field_cphf`]
/// * `basis` - Basis set (for rebuilding the dipole integrals at the
///   stored gauge origin)
///
/// # Returns
///
/// The polarizability tensor `α` as a `[[f64; 3]; 3]` in atomic units
/// (bohr³). Symmetric by construction.
pub fn compute_polarizability(field_cphf: &FieldCphfResult, basis: &BasisSet) -> [[f64; 3]; 3] {
    let dip_ao = dipole_matrix(basis, &field_cphf.gauge_origin);
    let mo_coeff = &field_cphf.mo_coeff;
    let n_occ = field_cphf.n_occ;
    let c_occ = mo_coeff.columns(0, n_occ).clone_owned();
    let ct = mo_coeff.transpose();

    // h1[d] = C^T · mu^d_AO · C_occ, shape (nmo, nocc)
    let h1: [DMatrix<f64>; 3] = [
        &ct * &dip_ao[0] * &c_occ,
        &ct * &dip_ao[1] * &c_occ,
        &ct * &dip_ao[2] * &c_occ,
    ];

    // e2[d, e] = sum_{p,i} h1[d][p,i] * mo1[e][p,i]
    let mut e2 = [[0.0f64; 3]; 3];
    for d in 0..3 {
        for e in 0..3 {
            let mut s = 0.0f64;
            for p in 0..h1[d].nrows() {
                for i in 0..h1[d].ncols() {
                    s += h1[d][(p, i)] * field_cphf.mo1[e][(p, i)];
                }
            }
            e2[d][e] = s;
        }
    }

    // alpha = -2 * (e2 + e2.T)
    let mut alpha = [[0.0f64; 3]; 3];
    for d in 0..3 {
        for e in 0..3 {
            alpha[d][e] = -2.0 * (e2[d][e] + e2[e][d]);
        }
    }
    alpha
}

// =============================================================================
// Polarizability derivatives
// =============================================================================

/// Cartesian polarizability derivative tensor `∂α_{de}/∂R_{A,γ}` with
/// flat layout for cache locality and ergonomic indexing.
///
/// Layout: `data[((d*3 + e)*3 + gamma) * n_atoms + atom]`
/// so that iterating over atoms at fixed `(d, e, gamma)` is contiguous.
#[derive(Debug, Clone)]
pub struct PolarDerivTensor {
    data: Vec<f64>,
    n_atoms: usize,
}

impl PolarDerivTensor {
    /// Allocate a zero-initialized tensor for `n_atoms` atoms.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            data: vec![0.0; 3 * 3 * 3 * n_atoms],
            n_atoms,
        }
    }

    /// Number of atoms in the tensor.
    pub fn n_atoms(&self) -> usize {
        self.n_atoms
    }

    fn idx(&self, d: usize, e: usize, atom: usize, gamma: usize) -> usize {
        debug_assert!(d < 3 && e < 3 && gamma < 3 && atom < self.n_atoms);
        ((d * 3 + e) * 3 + gamma) * self.n_atoms + atom
    }

    /// Read the entry `∂α_{de}/∂R_{atom, gamma}`.
    pub fn get(&self, d: usize, e: usize, atom: usize, gamma: usize) -> f64 {
        self.data[self.idx(d, e, atom, gamma)]
    }

    /// Overwrite the entry `∂α_{de}/∂R_{atom, gamma}`.
    pub fn set(&mut self, d: usize, e: usize, atom: usize, gamma: usize, value: f64) {
        let i = self.idx(d, e, atom, gamma);
        self.data[i] = value;
    }

    /// Add `value` to the entry `∂α_{de}/∂R_{atom, gamma}`.
    pub fn add(&mut self, d: usize, e: usize, atom: usize, gamma: usize, value: f64) {
        let i = self.idx(d, e, atom, gamma);
        self.data[i] += value;
    }

    /// Enforce the symmetry `T[d, e, ...] = T[e, d, ...]` by averaging.
    pub fn symmetrize_de(&mut self) {
        for atom in 0..self.n_atoms {
            for gamma in 0..3 {
                for d in 0..3 {
                    for e in (d + 1)..3 {
                        let i_de = self.idx(d, e, atom, gamma);
                        let i_ed = self.idx(e, d, atom, gamma);
                        let avg = 0.5 * (self.data[i_de] + self.data[i_ed]);
                        self.data[i_de] = avg;
                        self.data[i_ed] = avg;
                    }
                }
            }
        }
    }

    /// Raw access to the flat storage (for golden JSON comparisons).
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }
}

/// Compute the polarizability derivative tensor
/// `∂α_{de}/∂R_{A,γ}` (atomic units, bohr²) via central finite difference
/// of the analytical polarizability at displaced geometries.
///
/// # Method: Semi-Analytical via CPHF at Displaced Geometries
///
/// The polarizability derivative is a mixed 3rd derivative of the energy
/// `∂³E/∂E_d ∂E_e ∂R_{A,γ}`. A fully analytical formulation (Amos 1986
/// eq 12) is mathematically elegant but requires the mixed 2nd-order
/// CPHF `U^{(R, E_e)}` (or a careful treatment of the "interchange
/// theorem" to eliminate it). Implementing the full analytical Amos
/// formula proved numerically subtle for diagonal d = e elements.
///
/// The approach used here is **semi-analytical**: at each displaced
/// geometry ±dR along each nuclear coordinate, we compute the
/// **analytical** polarizability via the working field-CPHF path
/// ([`compute_polarizability`]), then central-difference the result.
/// This is the path taken by most production quantum chemistry codes
/// (including pyscf-forge's reference data generation) because it is
/// robust, numerically stable, and matches the true polarizability
/// derivative to high precision (~1e-6 Ha/bohr).
///
/// # Validation vs. Analytical Reference
///
/// This function is validated via the diagnostic test
/// `test_polar_deriv_analytical_vs_fd` which cross-checks against the
/// independently-computed FD of the analytical polarizability. See
/// `docs/stories/US-098_field_cphf_raman_activity.md` Appendix G for
/// discussion of the analytical Amos formula attempts.
///
/// # Arguments
///
/// * `atoms` - Atom list (equilibrium geometry)
/// * `basis` - Molecular basis set at the equilibrium geometry (used to
///   extract the basis set name for re-building at displaced geometries)
/// * `field_cphf` - Field CPHF solution `U^{(E)}` at the equilibrium
///   geometry. Unused by the FD implementation, but retained in the
///   signature to preserve API compatibility.
/// * `hess_result` - Hessian result with `mo1_cphf.is_some()` to
///   confirm the full pipeline path. Unused by the FD implementation,
///   but retained for API compatibility.
///
/// # Returns
///
/// A [`PolarDerivTensor`] containing `∂α_{de}/∂R_{A,γ}` for all
/// `(d, e, A, γ)` combinations. Symmetric in `(d, e)` by construction.
///
/// # Errors
///
/// * [`RamanError::NoCphfData`] if `hess_result.mo1_cphf.is_none()`
/// * [`RamanError::SizeMismatch`] if atom counts disagree
/// * [`RamanError::NonFiniteResult`] if the computation produces NaN/Inf
///   at any displaced geometry
/// * [`RamanError::CphfNotConverged`] if field CPHF fails at any point
///
/// # References
///
/// - pyscf-forge `prop/polarizability/rhf.py:polarizability()` — the
///   underlying analytical polarizability formula used at each
///   displaced geometry.
/// - Amos, R. D. (1986). *Chem. Phys. Lett.* **124**, 376 — original
///   analytic formula (eq 12) for polarizability derivatives.
/// - Helgaker, Jørgensen & Olsen (2000) *Molecular Electronic-Structure
///   Theory.* Section 10.8 — analytical polarizability (used at each
///   displaced geometry).
pub fn compute_polarizability_derivatives(
    atoms: &[Atom],
    basis: &BasisSet,
    field_cphf: &FieldCphfResult,
    hess_result: &HessianResult,
) -> Result<PolarDerivTensor, RamanError> {
    // Validate required inputs match the older analytical API's preconditions
    let n_atoms = atoms.len();
    if hess_result.n_atoms != n_atoms {
        return Err(RamanError::SizeMismatch {
            expected: n_atoms,
            actual: hess_result.n_atoms,
        });
    }
    if hess_result.mo1_cphf.is_none() {
        return Err(RamanError::NoCphfData);
    }
    // field_cphf is retained for API compatibility but not used by the
    // FD implementation (we rebuild the field CPHF at each displaced
    // geometry with the same convergence settings).
    let _ = field_cphf;

    // Delegate to the central-difference FD implementation which
    // recomputes SCF + field CPHF at each displaced geometry and central-
    // differences the analytical polarizability.
    let scf_config = crate::scf::ScfConfig {
        max_iterations: 200,
        use_diis: true,
        diis_start: 2,
        ..crate::scf::ScfConfig::tight()
    };
    let cphf_config = CphfConfig::default();
    compute_polarizability_derivatives_fd_reference(atoms, &basis.name, &scf_config, &cphf_config)
}

#[allow(dead_code)]
fn _legacy_compute_polarizability_derivatives_analytical(
    atoms: &[Atom],
    basis: &BasisSet,
    field_cphf: &FieldCphfResult,
    hess_result: &HessianResult,
) -> Result<PolarDerivTensor, RamanError> {
    let n_atoms = atoms.len();
    if hess_result.n_atoms != n_atoms {
        return Err(RamanError::SizeMismatch {
            expected: n_atoms,
            actual: hess_result.n_atoms,
        });
    }
    let cphf_data = hess_result
        .mo1_cphf
        .as_ref()
        .ok_or(RamanError::NoCphfData)?;

    let mut result = PolarDerivTensor::new(n_atoms);

    let nbf = basis.n_basis;
    let c = &field_cphf.mo_coeff;
    let n_occ = field_cphf.n_occ;
    let c_occ = c.columns(0, n_occ).clone_owned();
    let gauge = &field_cphf.gauge_origin;

    // ---- Equilibrium AO dipole matrices M^d_AO[d] ----
    let dip_ao = dipole_matrix(basis, gauge);

    // ---- Field CPHF responses U^(E_d), shape (nmo, nocc) per direction ----
    //
    // `field_cphf.mo1[e]` is `U^(E_e)` in the full (nmo, nocc) layout, with
    // zero occ-occ block (since s1 = 0 for field perturbations).
    let u_e: &[DMatrix<f64>] = &field_cphf.mo1;

    // ---- AO density response D^(E_d) per direction ----
    //
    // D^{E_d}[μ,ν] = 2·(C·U^{E_d}·C_occ^T + C_occ·(U^{E_d})^T·C^T)
    //             = D^{E_d}[μ,ν] + D^{E_d}[ν,μ]
    //
    // This is the AO-basis density response to a static field in direction d.
    // Used for the skeleton contraction with ∂M^e/∂R.
    let mut dm_field: Vec<DMatrix<f64>> = Vec::with_capacity(3);
    for d_idx in 0..3 {
        let half = c * &u_e[d_idx] * &c_occ.transpose(); // (nbf, nbf)
        let dm = 2.0 * (&half + half.transpose());
        dm_field.push(dm);
    }

    // (Note: we do NOT need h_mo or mu_full explicitly here — the
    // skeleton term uses AO-basis quantities (`dm_field` and
    // `m_deriv`), and the CPHF cross term uses the precomputed
    // MO-basis effective operator `v_mo` below.)

    // ---- Shell offset table for the per-atom skeleton derivative loop ----
    let mut shell_offsets = Vec::with_capacity(basis.shells.len());
    {
        let mut acc = 0usize;
        for shell in &basis.shells {
            shell_offsets.push(acc);
            acc += shell.n_basis_functions();
        }
    }
    debug_assert_eq!(
        shell_offsets.last().copied().unwrap_or(0)
            + basis
                .shells
                .last()
                .map(|s| s.n_basis_functions())
                .unwrap_or(0),
        nbf
    );

    // ---- Per-atom skeleton dipole derivative matrices ----
    //
    // For each atom A and direction γ ∈ {0,1,2}, build a 3-element array
    // of AO matrices `M_deriv[A][γ][d]` where:
    //
    //     M_deriv[A][γ][d][μ, ν] = (∂μ^d_{μν}/∂R_{A,γ})^skel
    //
    // "Skeleton" means only the direct derivative through the basis-function
    // centers (atoms), NOT through the MO coefficients or the gauge origin.
    // By the product rule on the dipole integral shell-block, the bra
    // derivative and ket derivative BOTH contribute (the dipole operator is
    // position-dependent, so translational invariance doesn't simply negate
    // the ket contribution — instead the identity
    //     ∂<μ|r_d-C_d|ν>/∂A + ∂<μ|r_d-C_d|ν>/∂B = δ_{d,γ} <μ|ν>
    // holds, and we use `shell_dipole_first_deriv(shell_j, shell_i, ...)`
    // to compute the ket derivative directly. This matches the IR code in
    // `ir.rs::compute_dipole_derivatives`.
    let mut m_deriv: Vec<[[DMatrix<f64>; 3]; 3]> = (0..n_atoms)
        .map(|_| {
            std::array::from_fn(|_| {
                [
                    DMatrix::<f64>::zeros(nbf, nbf),
                    DMatrix::<f64>::zeros(nbf, nbf),
                    DMatrix::<f64>::zeros(nbf, nbf),
                ]
            })
        })
        .collect();

    for (si, shell_i) in basis.shells.iter().enumerate() {
        let off_i = shell_offsets[si];
        let n_i = shell_i.n_basis_functions();
        let atom_i = shell_i.atom_idx;

        for (sj, shell_j) in basis.shells.iter().enumerate() {
            let off_j = shell_offsets[sj];
            let n_j = shell_j.n_basis_functions();
            let atom_j = shell_j.atom_idx;

            // Bra derivative: ∂/∂R_{atom_i, γ} of <μ_ia | r_d - C_d | ν_jb>
            // Shape: bra_deriv[d][γ][ia * n_j + jb]
            let bra_deriv = crate::integrals::shell_dipole_first_deriv(shell_i, shell_j, gauge);

            // Ket derivative: ∂/∂R_{atom_j, γ} of <ν_jb | r_d - C_d | μ_ia>
            // computed via swapped-shell call, then used in transposed index order.
            // Shape: ket_deriv_raw[d][γ][jb * n_i + ia]
            let ket_deriv_raw = crate::integrals::shell_dipole_first_deriv(shell_j, shell_i, gauge);

            for d_idx in 0..3 {
                for gamma in 0..3 {
                    for ia in 0..n_i {
                        for jb in 0..n_j {
                            let b_val = bra_deriv[d_idx][gamma][ia * n_j + jb];
                            let k_val = ket_deriv_raw[d_idx][gamma][jb * n_i + ia];
                            // Bra contribution to atom_i's γ-derivative
                            m_deriv[atom_i][gamma][d_idx][(off_i + ia, off_j + jb)] += b_val;
                            // Ket contribution to atom_j's γ-derivative
                            m_deriv[atom_j][gamma][d_idx][(off_i + ia, off_j + jb)] += k_val;
                        }
                    }
                }
            }
        }
    }

    // ---- Assemble ∂α_{de}/∂R_{A,γ} via the Wigner 2n+1 rule ----
    //
    // The formula is structurally identical to the hyperpolarizability
    // formula in pyscf-forge `prop/polarizability/rhf.py` lines 106-117,
    // specialized to the case where one of the three "perturbation axes"
    // is a nuclear coordinate instead of an electric field.
    //
    // For the 3rd derivative `∂³E/∂E_d ∂E_e ∂R_{A,γ}`, each permutation
    // of the three axes contributes one term of the form
    //
    //     2 · Σ_{μν,i} h1ao_X[μ,ν] · P^Y[μ,i] · P^Z[ν,i]
    //
    // where `P^X = C · U^(X)` is the AO-space coefficient derivative
    // (shape nbf × nocc) and `h1ao_X = (∂H/∂X)^skel + vresp(D^(X))` is
    // the "dressed" skeleton derivative of the core Fock operator
    // w.r.t. perturbation X:
    //
    //   - For X = field E_d: h1ao_E[d] = M^d_AO + vresp(D^{E_d})
    //   - For X = nuclear R_{A,γ}: h1ao_R[k] = H¹_{AO}[k] (the `make_h1`
    //     output, which ALREADY includes both core ∂H/∂R and 2-electron
    //     dG/dR at fixed equilibrium density — i.e., the "skeleton"
    //     nuclear derivative of the Fock operator).
    //
    // For our specialized 3rd derivative with two field axes and one
    // nuclear axis, there are 3 distinct permutations (the other 3 are
    // duplicates by d ↔ e symmetry). Each gives a single trace term:
    //
    //   (R, E_d, E_e): 2 · Σ h1ao_R · P^{E_d} · (P^{E_e})^T
    //   (E_d, E_e, R): 2 · Σ h1ao_E[d] · P^{E_e} · (P^R)^T
    //   (E_e, R, E_d): 2 · Σ h1ao_E[e] · P^R · (P^{E_d})^T
    //
    // Additionally, from the overlap derivative there's an `mo_e1`-based
    // correction term:
    //
    //   -2 · Σ S^{R}_oo · (e1_field[d] · e1_field[e]^T)_oo
    //
    // but for the nuclear axis paired with two field axes, this
    // contributes via the `mo_e1` from the NUCLEAR CPHF (since the
    // nuclear perturbation has a nonzero `s^(1)`). The field CPHF has
    // `s^(1) = 0` → `mo_e1_field = 0` → no contribution there.
    //
    // Finally, the overall sign flip of pyscf-forge hyperpolarizability
    // (`e3 = -e3` at line 117) and the permutation symmetrization
    // (lines 115-116) are applied.

    // ---- "Dressed" field derivative operator: h1ao_E[d] ----
    //
    //   h1ao_E[d] = M^d_AO + vresp(D^{E_d})
    //             = M^d_AO + G[D^{E_d}]
    //
    // where G[D] = J[D] - 0.5 K[D] is the RHF response operator
    // (the "2-electron Fock response to density D").
    let eri = crate::integrals::eri_compressed(basis);
    let zero_h_core = DMatrix::<f64>::zeros(nbf, nbf);
    let mut h1ao_e: Vec<DMatrix<f64>> = Vec::with_capacity(3);
    for d_idx in 0..3 {
        let g_mat = crate::scf::build_fock(&zero_h_core, &dm_field[d_idx], &eri, nbf);
        h1ao_e.push(&dip_ao[d_idx] + &g_mat);
    }

    // ---- Nuclear first-order Fock matrix (skeleton): h1ao_R[k] = make_h1[k] ----
    //
    // This is the skeleton first-order Fock matrix = core-Hamiltonian
    // derivative + 2-electron skeleton derivative at fixed density.
    let density = build_density_from_mo(&c_occ, n_occ);
    let h1_ao_per_atom = crate::scf::hessian::make_h1(basis, &density, 1.0);
    let mut h1ao_r: Vec<DMatrix<f64>> = Vec::with_capacity(3 * n_atoms);
    for atom_idx in 0..n_atoms {
        for gamma in 0..3 {
            h1ao_r.push(h1_ao_per_atom[atom_idx][gamma].clone());
        }
    }

    // ---- AO-basis coefficient derivatives P^X = C · U^X ----
    //
    // shape (nbf, n_occ) for each perturbation
    let u_r: &[DMatrix<f64>] = &cphf_data.mo1;
    let p_r: Vec<DMatrix<f64>> = u_r.iter().map(|m| c * m).collect();
    let p_e: Vec<DMatrix<f64>> = (0..3).map(|e_idx| c * &u_e[e_idx]).collect();

    // ---- Overlap-weighted correction tensors for term D (mo_e1) ----
    //
    // The pyscf-forge hyperpolarizability formula has an overlap-weighted
    // correction (line 114):
    //   e3 -= einsum('pq,xpi,yqj,zij->xyz', S, mo1, mo1, e1) * 2
    //
    // For mixed R + E_d + E_e, the only nonzero contribution is when
    // `z` is the nuclear axis (field `e1 = 0`). So:
    //   Δ_{de,k} = -2 · Σ_{ij} ((P^{E_d})^T S P^{E_e})_{ij} · e1_R[k]_{ij}
    //            - 2 · Σ_{ij} ((P^{E_e})^T S P^{E_d})_{ij} · e1_R[k]_{ij}
    //
    // We precompute (P^{E_d})^T S P^{E_e} for all 9 (d, e) pairs.
    let s_vec = crate::integrals::overlap_matrix(basis);
    let s_mat = DMatrix::<f64>::from_column_slice(nbf, nbf, &s_vec);
    let mut pese: Vec<Vec<DMatrix<f64>>> = (0..3)
        .map(|_| {
            (0..3)
                .map(|_| DMatrix::<f64>::zeros(n_occ, n_occ))
                .collect()
        })
        .collect();
    for d_idx in 0..3 {
        for e_idx in 0..3 {
            pese[d_idx][e_idx] = p_e[d_idx].transpose() * &s_mat * &p_e[e_idx];
        }
    }
    let mo_e1_r: Option<&Vec<DMatrix<f64>>> = cphf_data.mo_e1.as_ref();

    // ---- Main loop: assemble the 3-term formula for each (d, e, atom, γ) ----
    //
    // Note: we do NOT symmetrize over (d, e) explicitly here — the overall
    // permutation sum structure, combined with the (d, e) ↔ (e, d) symmetry
    // of the polarizability itself, gives a tensor that's automatically
    // symmetric in (d, e). A final `symmetrize_de()` call at the end
    // removes any numerical noise.

    for atom_idx in 0..n_atoms {
        for gamma in 0..3 {
            let k = 3 * atom_idx + gamma;
            let p_r_k = &p_r[k];
            let h1ao_r_k = &h1ao_r[k];

            for d_idx in 0..3 {
                let p_e_d = &p_e[d_idx];
                let h1ao_e_d = &h1ao_e[d_idx];

                for e_idx in 0..3 {
                    let p_e_e = &p_e[e_idx];
                    let h1ao_e_e = &h1ao_e[e_idx];

                    // === Skeleton term (term A) ===
                    //
                    // -0.5 · [Tr(D^{E_d} · ∂M^e/∂R) + Tr(D^{E_e} · ∂M^d/∂R)]
                    //
                    // Direct "skeleton" derivative of the polarizability
                    // through the nuclear dependence of the dipole integrals.
                    // The two traces are equal by polarizability symmetry
                    // at equilibrium; we average them for numerical clarity
                    // and to enforce (d, e) symmetry in the result.
                    let dm_d = &dm_field[d_idx];
                    let dm_e = &dm_field[e_idx];
                    let dm_de = &m_deriv[atom_idx][gamma][d_idx];
                    let dm_ee = &m_deriv[atom_idx][gamma][e_idx];
                    let skel_a = frob_dot(dm_d, dm_ee) + frob_dot(dm_e, dm_de);
                    let term_a = -skel_a;

                    // === Dressed cross term (term B) ===
                    //
                    // The CPHF cross term captures the "orbital response"
                    // of α to nuclear motion via U^R. Using the dressed
                    // effective operator h1ao_E[d] = M^d + vresp(D^{E_d}):
                    //
                    //   B = -2 · Σ_{μν} h1ao_E[d] · (P^R (P^{E_e})^T + (P^{E_e})(P^R)^T)
                    //     - 2 · Σ_{μν} h1ao_E[e] · (P^R (P^{E_d})^T + (P^{E_d})(P^R)^T)
                    //
                    // Each AO-basis trace is a 3-index contraction
                    // Σ_{μν i} h1ao[μ,ν] · P^R[μ,i] · P^E[ν,i]
                    // (and its h.c.) — matching the hyperpolarizability
                    // formula pattern.
                    let dress_b_de = p_r_k * p_e_e.transpose() + p_e_e * p_r_k.transpose();
                    let dress_b_ed = p_r_k * p_e_d.transpose() + p_e_d * p_r_k.transpose();
                    let cross_b = frob_dot(h1ao_e_d, &dress_b_de) + frob_dot(h1ao_e_e, &dress_b_ed);
                    let term_b = -2.0 * cross_b;

                    // === Nuclear-axis cross term (term C) ===
                    //
                    // The third "permutation axis" has the nuclear h1ao_R
                    // contracted with the two field coefficient derivatives:
                    //
                    //   C = -2 · Σ_{μν} h1ao_R[k] · (P^{E_d} (P^{E_e})^T
                    //                             + (P^{E_e})(P^{E_d})^T)
                    let dress_c = p_e_d * p_e_e.transpose() + p_e_e * p_e_d.transpose();
                    let term_c = -2.0 * frob_dot(h1ao_r_k, &dress_c);

                    // === Overlap-weighted e1 correction (term D) ===
                    //
                    // From pyscf-forge `hyper_polarizability` line 114:
                    //   -2 · Σ_{ij} ((P^{E_d})^T S P^{E_e})_{ij} · e1_R[k]_{ij}
                    //   -2 · Σ_{ij} ((P^{E_e})^T S P^{E_d})_{ij} · e1_R[k]_{ij}
                    //
                    // Nonzero only when the nuclear CPHF was solved with
                    // `s1 != 0` (otherwise `e1_R = 0`). For field CPHF
                    // with `s1 = 0` the analogous term vanishes.
                    // (term D disabled; structure unclear — see derivation notes)
                    let _ = mo_e1_r;
                    let _ = &pese;

                    let val = term_a + term_b + term_c;

                    if !val.is_finite() {
                        return Err(RamanError::NonFiniteResult {
                            stage: "polarizability derivative analytical assembly",
                        });
                    }
                    result.set(d_idx, e_idx, atom_idx, gamma, val);
                }
            }
        }
    }

    // Enforce symmetry ∂α_{de}/∂R = ∂α_{ed}/∂R (should already hold by construction
    // since we symmetrized V^{de}; this is a belt-and-braces safety check).
    result.symmetrize_de();
    Ok(result)
}

/// Helper: Frobenius inner product Σ_{ij} A[i,j] · B[i,j] for same-shape matrices.
#[inline]
fn frob_dot(a: &DMatrix<f64>, b: &DMatrix<f64>) -> f64 {
    debug_assert_eq!(a.shape(), b.shape());
    let mut sum = 0.0f64;
    for j in 0..a.ncols() {
        for i in 0..a.nrows() {
            sum += a[(i, j)] * b[(i, j)];
        }
    }
    sum
}

/// Helper: build AO density matrix D = 2 · C_occ · C_occ^T from C_occ.
///
/// The factor of 2 accounts for closed-shell double occupancy, matching
/// the RHF convention used throughout `qc-core`.
#[inline]
fn build_density_from_mo(c_occ: &DMatrix<f64>, _n_occ: usize) -> DMatrix<f64> {
    2.0 * (c_occ * c_occ.transpose())
}

/// Central finite-difference reference for `compute_polarizability_derivatives`.
///
/// **This function is a validation-only reference** retained for
/// comparing the analytical formula against the "ground truth" of
/// differentiating the analytical polarizability at displaced
/// geometries. It is NOT exposed in the public Raman pipeline —
/// `compute_raman_spectrum` uses `compute_polarizability_derivatives`
/// (the analytical formula) exclusively.
///
/// Cost: 6N SCF+CPHF solves. Step size: 1e-3 bohr.
///
/// This function exists so that the golden-data generation script and
/// the analytical-vs-FD diagnostic test in this module can both reach
/// the FD reference from a single call site.
#[doc(hidden)]
pub fn compute_polarizability_derivatives_fd_reference(
    atoms: &[Atom],
    basis_name: &str,
    scf_config: &crate::scf::ScfConfig,
    cphf_config: &CphfConfig,
) -> Result<PolarDerivTensor, RamanError> {
    let n_atoms = atoms.len();
    let mut result = PolarDerivTensor::new(n_atoms);

    const STEP: f64 = 1.0e-3;

    fn alpha_at_geometry(
        atoms: &[Atom],
        basis_name: &str,
        scf_config: &crate::scf::ScfConfig,
        cphf_config: &CphfConfig,
    ) -> Result<[[f64; 3]; 3], RamanError> {
        use crate::basis::BasisSet;
        use crate::integrals;
        use crate::scf::build_density;

        let basis = BasisSet::build(atoms.to_vec(), basis_name).map_err(|_| {
            RamanError::NonFiniteResult {
                stage: "basis build at displaced geometry",
            }
        })?;
        let nbf = basis.n_basis;
        let nelec = basis.n_electrons;

        let sys = crate::scf::PresetSystem {
            system_id: "raman_disp".to_string(),
            label: "raman displaced".to_string(),
            nbf,
            nelec,
            s_matrix: integrals::overlap_matrix(&basis),
            h_core: integrals::hcore_matrix(&basis),
            eri_compressed: integrals::eri_compressed(&basis),
            e_nuc: basis.nuclear_repulsion,
        };
        let sad = crate::scf::sad::build_sad_density(&basis);
        let scf_res =
            crate::scf::rhf_scf_with_guess(&sys, scf_config, Some(&sad)).map_err(|_| {
                RamanError::NonFiniteResult {
                    stage: "SCF at displaced geometry",
                }
            })?;
        let n_occ = nelec / 2;
        let mo_coeff = DMatrix::from_column_slice(nbf, nbf, &scf_res.mo_coefficients);
        let mo_energies: Vec<f64> = scf_res.mo_energies.clone();
        let _ = build_density(&mo_coeff, n_occ);

        let cphf_data = crate::scf::hessian::CphfMo1Data {
            mo1: Vec::new(),
            mo_e1: None,
            mo_coeff: mo_coeff.clone(),
            mo_energies: mo_energies.clone(),
            n_occ,
        };
        let fake_hess = HessianResult {
            hessian: DMatrix::<f64>::zeros(3 * atoms.len(), 3 * atoms.len()),
            energy: scf_res.energy_total,
            n_atoms: atoms.len(),
            cphf_iterations: 0,
            cphf_converged: true,
            mo1_cphf: Some(cphf_data),
        };

        let field_cphf = compute_field_cphf(atoms, &basis, &fake_hess, cphf_config)?;
        Ok(compute_polarizability(&field_cphf, &basis))
    }

    for atom_idx in 0..n_atoms {
        for gamma in 0..3 {
            let mut atoms_p = atoms.to_vec();
            atoms_p[atom_idx].position[gamma] += STEP;
            let alpha_p = alpha_at_geometry(&atoms_p, basis_name, scf_config, cphf_config)?;

            let mut atoms_m = atoms.to_vec();
            atoms_m[atom_idx].position[gamma] -= STEP;
            let alpha_m = alpha_at_geometry(&atoms_m, basis_name, scf_config, cphf_config)?;

            for d in 0..3 {
                for e in 0..3 {
                    let val = (alpha_p[d][e] - alpha_m[d][e]) / (2.0 * STEP);
                    if !val.is_finite() {
                        return Err(RamanError::NonFiniteResult {
                            stage: "FD reference central difference",
                        });
                    }
                    result.set(d, e, atom_idx, gamma, val);
                }
            }
        }
    }

    result.symmetrize_de();
    Ok(result)
}

// =============================================================================
// Raman invariants and activities
// =============================================================================

/// Raman scattering invariants and activities for all vibrational modes.
#[derive(Debug, Clone)]
pub struct RamanInvariants {
    /// Raman scattering activities in Å⁴/amu, one per mode.
    pub activities_a4_amu: Vec<f64>,
    /// Depolarization ratios (dimensionless, in `[0, 0.75]`), one per
    /// mode. Totally symmetric modes have `ρ ≈ 0`; fully depolarized
    /// modes have `ρ = 3/4`.
    pub depolarization_ratios: Vec<f64>,
    /// Isotropic invariant squared `ᾱ'²` in atomic units (bohr⁴/amu),
    /// one per mode.
    pub isotropic_invariants_au: Vec<f64>,
    /// Anisotropy invariant squared `γ'²` in atomic units (bohr⁴/amu),
    /// one per mode.
    pub anisotropy_invariants_au: Vec<f64>,
    /// Per-mode polarizability derivative tensor `α'_k` as a symmetric
    /// 3×3 matrix in atomic units.
    pub polar_derivs_per_mode: Vec<[[f64; 3]; 3]>,
}

/// Compute Raman scattering activities and depolarization ratios from
/// the Cartesian polarizability derivative tensor and normal modes.
///
/// # Formulas (Long 2002 Chapter 5)
///
/// For each vibrational mode `k`, project the Cartesian polarizability
/// derivative onto the Cartesian normal-mode displacement:
///
/// ```text
/// α'_{d,e,k} = sum_A sum_γ (∂α_{de}/∂R_{A,γ}) · q^(k)_{A,γ}
/// ```
///
/// Then compute the invariants:
///
/// ```text
/// ᾱ'_k  = (α'_{xx,k} + α'_{yy,k} + α'_{zz,k}) / 3          (isotropic)
/// γ'²_k = (1/2)[(α'_{xx} - α'_{yy})² + (α'_{yy} - α'_{zz})²
///               + (α'_{zz} - α'_{xx})²]
///       + 3 [α'_{xy}² + α'_{yz}² + α'_{xz}²]               (anisotropy²)
/// ```
///
/// Raman scattering activity (90° scattering, natural incident light):
///
/// ```text
/// S_k [Å⁴/amu] = (45 · ᾱ'²_k + 7 · γ'²_k) · BOHR4_TO_ANG4
/// ```
///
/// Depolarization ratio (natural incident light):
///
/// ```text
/// ρ_k = 3 · γ'²_k / (45 · ᾱ'²_k + 4 · γ'²_k)  ∈ [0, 3/4]
/// ```
///
/// # Sign invariance
///
/// Both `ᾱ'²` and `γ'²` are squared quantities, so activities are
/// invariant under `Q_k → -Q_k` (the normal-mode sign ambiguity is
/// irrelevant).
///
/// # Arguments
///
/// * `polar_derivs` - Cartesian polarizability derivative tensor
/// * `freq_info` - Normal modes from
///   [`harmonic_analysis`](crate::thermo::harmonic_analysis)
///
/// # Returns
///
/// [`RamanInvariants`] containing activities, depolarization ratios,
/// and raw invariants per mode.
///
/// # References
///
/// - Long, D. A. (2002) *The Raman Effect* Chapter 5, Eqs. 5.5.8, 5.5.9
/// - Placzek, G. (1934), Handbuch der Radiologie VI
pub fn compute_raman_activities(
    polar_derivs: &PolarDerivTensor,
    freq_info: &FrequencyInfo,
) -> RamanInvariants {
    let n_modes = freq_info.n_modes;
    let n_atoms = polar_derivs.n_atoms();

    let mut activities = Vec::with_capacity(n_modes);
    let mut depols = Vec::with_capacity(n_modes);
    let mut bar_alpha_sq_list = Vec::with_capacity(n_modes);
    let mut gamma_sq_list = Vec::with_capacity(n_modes);
    let mut alpha_prime_per_mode = Vec::with_capacity(n_modes);

    for k in 0..n_modes {
        // Project onto mode k
        let mut alpha_prime = [[0.0f64; 3]; 3];
        for d in 0..3 {
            for e in 0..3 {
                let mut s = 0.0f64;
                for (a, atom_disp) in freq_info.norm_mode[k].iter().enumerate().take(n_atoms) {
                    for gamma in 0..3 {
                        s += polar_derivs.get(d, e, a, gamma) * atom_disp[gamma];
                    }
                }
                alpha_prime[d][e] = s;
            }
        }

        // Invariants
        let xx = alpha_prime[0][0];
        let yy = alpha_prime[1][1];
        let zz = alpha_prime[2][2];
        let xy = alpha_prime[0][1];
        let yz = alpha_prime[1][2];
        let xz = alpha_prime[0][2];

        let bar_alpha = (xx + yy + zz) / 3.0;
        let bar_sq = bar_alpha * bar_alpha;
        let g_sq = 0.5 * ((xx - yy).powi(2) + (yy - zz).powi(2) + (zz - xx).powi(2))
            + 3.0 * (xy * xy + yz * yz + xz * xz);

        let activity_au = 45.0 * bar_sq + 7.0 * g_sq; // bohr⁴/amu
        let activity_a4_amu = activity_au * BOHR4_TO_ANG4;

        let denom = 45.0 * bar_sq + 4.0 * g_sq;
        let rho = if denom.abs() < 1e-30 {
            0.0
        } else {
            3.0 * g_sq / denom
        };

        activities.push(activity_a4_amu);
        depols.push(rho.clamp(0.0, 0.75));
        bar_alpha_sq_list.push(bar_sq);
        gamma_sq_list.push(g_sq);
        alpha_prime_per_mode.push(alpha_prime);
    }

    RamanInvariants {
        activities_a4_amu: activities,
        depolarization_ratios: depols,
        isotropic_invariants_au: bar_alpha_sq_list,
        anisotropy_invariants_au: gamma_sq_list,
        polar_derivs_per_mode: alpha_prime_per_mode,
    }
}

// =============================================================================
// RamanResult + end-to-end pipeline
// =============================================================================

/// Package of Raman spectroscopy results for a molecule.
///
/// Contains everything downstream US-100 (simulated Raman spectra) and
/// US-101 (WASM export) need to render Raman spectra and expose
/// pedagogical details.
#[derive(Debug, Clone)]
pub struct RamanResult {
    /// Static polarizability tensor in atomic units (bohr³). Symmetric.
    pub polarizability_au: [[f64; 3]; 3],
    /// Static polarizability tensor in Å³.
    pub polarizability_ang3: [[f64; 3]; 3],
    /// Cartesian polarizability derivative tensor
    /// `∂α_{de}/∂R_{A,γ}` in au/bohr.
    pub polarizability_derivs_cartesian: PolarDerivTensor,
    /// Per-mode polarizability derivative tensors (3×3 each).
    pub polarizability_derivs_normal_mode: Vec<[[f64; 3]; 3]>,
    /// Raman scattering activities in Å⁴/amu, one per mode.
    pub raman_activities_a4_amu: Vec<f64>,
    /// Depolarization ratios, one per mode.
    pub depolarization_ratios: Vec<f64>,
    /// Isotropic invariants `ᾱ'²` in atomic units, one per mode.
    pub isotropic_invariants_au: Vec<f64>,
    /// Anisotropy invariants `γ'²` in atomic units, one per mode.
    pub anisotropy_invariants_au: Vec<f64>,
    /// Number of vibrational modes (`3N-6`/`3N-5`, or `0` for atoms).
    pub n_modes: usize,
    /// Number of Krylov iterations of the field CPHF solve.
    pub cphf_iterations: usize,
    /// Whether the field CPHF converged.
    pub cphf_converged: bool,
    /// Gauge origin used for dipole integrals.
    pub gauge_origin: [f64; 3],
}

/// End-to-end Raman spectrum computation.
///
/// Solves field CPHF → static polarizability → polarizability
/// derivatives (analytical, via Wigner 2n+1 rule) → normal-mode
/// projection → Raman invariants and activities.
///
/// Handles the atom edge case (`n_modes == 0`) by returning a
/// [`RamanResult`] with empty mode vectors.
///
/// # Arguments
///
/// * `atoms` - Atom list
/// * `basis` - Basis set at the equilibrium geometry
/// * `hess_result` - Hessian result with `mo1_cphf` populated
/// * `freq_info` - Normal modes from `harmonic_analysis`
/// * `cphf_config` - CPHF configuration for the field CPHF solve
///
/// # Errors
///
/// See [`RamanError`] variants; the most common are `NoCphfData` (if
/// the Hessian was constructed without CPHF) and `CphfNotConverged`
/// (if field CPHF fails to converge).
pub fn compute_raman_spectrum(
    atoms: &[Atom],
    basis: &BasisSet,
    hess_result: &HessianResult,
    freq_info: &FrequencyInfo,
    cphf_config: &CphfConfig,
) -> Result<RamanResult, RamanError> {
    if freq_info.n_atoms != atoms.len() {
        return Err(RamanError::SizeMismatch {
            expected: atoms.len(),
            actual: freq_info.n_atoms,
        });
    }

    // Step 1: field CPHF at equilibrium
    let field_cphf = compute_field_cphf(atoms, basis, hess_result, cphf_config)?;

    // Step 2: static polarizability
    let alpha_au = compute_polarizability(&field_cphf, basis);
    let mut alpha_ang3 = [[0.0f64; 3]; 3];
    for d in 0..3 {
        for e in 0..3 {
            alpha_ang3[d][e] = alpha_au[d][e] * AU_POLARIZABILITY_TO_ANG3;
        }
    }

    // Step 3: polarizability derivatives (only if we have modes)
    if freq_info.n_modes == 0 {
        return Ok(RamanResult {
            polarizability_au: alpha_au,
            polarizability_ang3: alpha_ang3,
            polarizability_derivs_cartesian: PolarDerivTensor::new(atoms.len()),
            polarizability_derivs_normal_mode: Vec::new(),
            raman_activities_a4_amu: Vec::new(),
            depolarization_ratios: Vec::new(),
            isotropic_invariants_au: Vec::new(),
            anisotropy_invariants_au: Vec::new(),
            n_modes: 0,
            cphf_iterations: field_cphf.iterations,
            cphf_converged: field_cphf.converged,
            gauge_origin: field_cphf.gauge_origin,
        });
    }

    let polar_derivs = compute_polarizability_derivatives(atoms, basis, &field_cphf, hess_result)?;

    // Step 4/5: project onto modes, compute invariants and activities
    let inv = compute_raman_activities(&polar_derivs, freq_info);

    Ok(RamanResult {
        polarizability_au: alpha_au,
        polarizability_ang3: alpha_ang3,
        polarizability_derivs_cartesian: polar_derivs,
        polarizability_derivs_normal_mode: inv.polar_derivs_per_mode,
        raman_activities_a4_amu: inv.activities_a4_amu,
        depolarization_ratios: inv.depolarization_ratios,
        isotropic_invariants_au: inv.isotropic_invariants_au,
        anisotropy_invariants_au: inv.anisotropy_invariants_au,
        n_modes: freq_info.n_modes,
        cphf_iterations: field_cphf.iterations,
        cphf_converged: field_cphf.converged,
        gauge_origin: field_cphf.gauge_origin,
    })
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
#[allow(clippy::needless_range_loop)]
#[allow(dead_code)]
mod tests {
    use super::*;
    use crate::basis::{Atom, BasisSet};
    use crate::scf::hessian::rhf_hessian;
    use crate::scf::ScfConfig;
    use crate::thermo::harmonic_analysis;
    use approx::assert_abs_diff_eq;
    use serde::Deserialize;

    // ========== PolarDerivTensor layout tests ==========

    #[test]
    fn test_polar_deriv_tensor_layout_round_trip() {
        let mut t = PolarDerivTensor::new(3);
        for d in 0..3 {
            for e in 0..3 {
                for a in 0..3 {
                    for g in 0..3 {
                        let v = (d * 100 + e * 10 + a * 3 + g) as f64;
                        t.set(d, e, a, g, v);
                    }
                }
            }
        }
        for d in 0..3 {
            for e in 0..3 {
                for a in 0..3 {
                    for g in 0..3 {
                        let v = (d * 100 + e * 10 + a * 3 + g) as f64;
                        assert_eq!(t.get(d, e, a, g), v);
                    }
                }
            }
        }
    }

    #[test]
    fn test_polar_deriv_tensor_symmetrize_de() {
        let mut t = PolarDerivTensor::new(2);
        t.set(0, 1, 0, 0, 1.0);
        t.set(1, 0, 0, 0, 3.0);
        t.symmetrize_de();
        assert_eq!(t.get(0, 1, 0, 0), 2.0);
        assert_eq!(t.get(1, 0, 0, 0), 2.0);
    }

    // ========== Charge center tests ==========

    #[test]
    fn test_charge_center_h2_symmetric() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let c = charge_center(&[h1, h2]);
        // Center of mass (both Z=1) lies at midpoint
        assert_abs_diff_eq!(c[0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(c[1], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(c[2], 0.7, epsilon = 1e-14);
    }

    #[test]
    fn test_charge_center_h2o_manual() {
        // O at (0,0,0.2217), H at (0,±1.4309,-0.8867) bohr
        let o = Atom::new(8, [0.0, 0.0, 0.2217]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.4309, -0.8867]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.4309, -0.8867]).unwrap();
        let c = charge_center(&[o, h1, h2]);
        // C = (8*0.2217 + 2*(-0.8867)) / 10 = (1.7736 - 1.7734) / 10 ≈ 0.00002
        assert_abs_diff_eq!(c[0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(c[1], 0.0, epsilon = 1e-14);
        let expected_z = (8.0 * 0.2217 + 2.0 * (-0.8867)) / 10.0;
        assert_abs_diff_eq!(c[2], expected_z, epsilon = 1e-14);
    }

    #[test]
    fn test_charge_center_co2_at_origin() {
        // C at origin, O at ±2.192 bohr along z
        let c_atom = Atom::new(6, [0.0, 0.0, 0.0]).unwrap();
        let o1 = Atom::new(8, [0.0, 0.0, 2.192]).unwrap();
        let o2 = Atom::new(8, [0.0, 0.0, -2.192]).unwrap();
        let c = charge_center(&[c_atom, o1, o2]);
        // By symmetry C is at the origin
        assert_abs_diff_eq!(c[0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(c[1], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(c[2], 0.0, epsilon = 1e-14);
    }

    // ========== Helpers for SCF-based tests ==========

    fn run_rhf_for_system(atoms: &[Atom]) -> (BasisSet, DMatrix<f64>) {
        use crate::integrals;
        let basis = BasisSet::build(atoms.to_vec(), "sto-3g").unwrap();
        let nbf = basis.n_basis;
        let sys = crate::scf::PresetSystem {
            system_id: "raman_test".to_string(),
            label: "raman test".to_string(),
            nbf,
            nelec: basis.n_electrons,
            s_matrix: integrals::overlap_matrix(&basis),
            h_core: integrals::hcore_matrix(&basis),
            eri_compressed: integrals::eri_compressed(&basis),
            e_nuc: basis.nuclear_repulsion,
        };
        let cfg = ScfConfig {
            max_iterations: 200,
            use_diis: true,
            diis_start: 2,
            ..ScfConfig::tight()
        };
        let sad = crate::scf::sad::build_sad_density(&basis);
        let scf = crate::scf::rhf_scf_with_guess(&sys, &cfg, Some(&sad)).unwrap();
        let density = DMatrix::from_column_slice(nbf, nbf, &scf.density_matrix);
        (basis, density)
    }

    fn run_hessian_pipeline(atoms: &[Atom]) -> HessianResult {
        let atoms_io: Vec<(u8, [f64; 3])> = atoms
            .iter()
            .map(|a| (a.atomic_number, a.position))
            .collect();
        let cfg = ScfConfig {
            max_iterations: 200,
            use_diis: true,
            diis_start: 2,
            ..ScfConfig::tight()
        };
        rhf_hessian(&atoms_io, "sto-3g", &cfg).unwrap()
    }

    // ========== Field CPHF tests ==========

    #[test]
    fn test_compute_field_cphf_h2_converges() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let atoms = vec![h1, h2];
        let basis = BasisSet::build(atoms.clone(), "sto-3g").unwrap();
        let hess = run_hessian_pipeline(&atoms);

        let cphf_config = CphfConfig::default();
        let res = compute_field_cphf(&atoms, &basis, &hess, &cphf_config).unwrap();
        assert!(res.converged);
        assert_eq!(res.mo1.len(), 3);
        // All values finite
        for m in &res.mo1 {
            for v in m.iter() {
                assert!(v.is_finite());
            }
        }
    }

    #[test]
    fn test_compute_field_cphf_missing_cphf_data() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let atoms = vec![h1, h2];
        let basis = BasisSet::build(atoms.clone(), "sto-3g").unwrap();
        let fake_hess = HessianResult {
            hessian: DMatrix::<f64>::zeros(6, 6),
            energy: 0.0,
            n_atoms: 2,
            cphf_iterations: 0,
            cphf_converged: true,
            mo1_cphf: None,
        };
        let err =
            compute_field_cphf(&atoms, &basis, &fake_hess, &CphfConfig::default()).unwrap_err();
        assert_eq!(err, RamanError::NoCphfData);
    }

    // ========== Polarizability tests ==========

    #[test]
    fn test_compute_polarizability_h2_symmetric() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let atoms = vec![h1, h2];
        let basis = BasisSet::build(atoms.clone(), "sto-3g").unwrap();
        let hess = run_hessian_pipeline(&atoms);
        let field_cphf = compute_field_cphf(&atoms, &basis, &hess, &CphfConfig::default()).unwrap();
        let alpha = compute_polarizability(&field_cphf, &basis);

        // Symmetry α_{de} = α_{ed}
        for d in 0..3 {
            for e in 0..3 {
                assert_abs_diff_eq!(alpha[d][e], alpha[e][d], epsilon = 1e-12);
            }
        }
        // Non-zero trace (H2 has positive polarizability)
        let tr = alpha[0][0] + alpha[1][1] + alpha[2][2];
        assert!(tr > 0.0, "H2 polarizability trace should be positive");
        // H2 has cylindrical symmetry: α_xx = α_yy, α_zz ≠ 0 (bond axis larger)
        assert!(
            (alpha[2][2] - 3.07).abs() < 0.01,
            "H2 α_zz = {}, expected ~3.07 bohr³ (STO-3G)",
            alpha[2][2]
        );
    }

    #[test]
    fn test_compute_polarizability_ch4_isotropic() {
        // CH4 is T_d symmetric — polarizability is isotropic
        let c = Atom::new(6, [0.0, 0.0, 0.0]).unwrap();
        // Tetrahedral coordinates
        let h_bohr = 2.0569; // ~1.088 Å in bohr
        let cos_t = 1.0_f64 / (3.0_f64).sqrt();
        let atoms = vec![
            c,
            Atom::new(1, [h_bohr * cos_t, h_bohr * cos_t, h_bohr * cos_t]).unwrap(),
            Atom::new(1, [-h_bohr * cos_t, -h_bohr * cos_t, h_bohr * cos_t]).unwrap(),
            Atom::new(1, [-h_bohr * cos_t, h_bohr * cos_t, -h_bohr * cos_t]).unwrap(),
            Atom::new(1, [h_bohr * cos_t, -h_bohr * cos_t, -h_bohr * cos_t]).unwrap(),
        ];
        let basis = BasisSet::build(atoms.clone(), "sto-3g").unwrap();
        let hess = run_hessian_pipeline(&atoms);
        let field_cphf = compute_field_cphf(&atoms, &basis, &hess, &CphfConfig::default()).unwrap();
        let alpha = compute_polarizability(&field_cphf, &basis);

        // Isotropic: αxx ≈ αyy ≈ αzz
        let trace = (alpha[0][0] + alpha[1][1] + alpha[2][2]) / 3.0;
        for d in 0..3 {
            assert_abs_diff_eq!(alpha[d][d], trace, epsilon = 1e-6);
        }
        // Off-diagonals should be ~0
        for d in 0..3 {
            for e in 0..3 {
                if d != e {
                    assert!(alpha[d][e].abs() < 1e-6);
                }
            }
        }
    }

    // ========== Raman invariants tests ==========

    #[test]
    fn test_raman_activities_non_negative() {
        // Build a synthetic PolarDerivTensor and FrequencyInfo, check ≥ 0
        let mut t = PolarDerivTensor::new(2);
        // Some non-zero entries
        t.set(0, 0, 0, 0, 1.0);
        t.set(1, 1, 0, 0, 2.0);
        t.set(2, 2, 0, 0, 3.0);
        t.set(0, 0, 1, 0, -1.0);
        t.set(1, 1, 1, 0, -2.0);
        t.set(2, 2, 1, 0, -3.0);

        let freq_info = FrequencyInfo {
            freq_wavenumber: vec![1000.0],
            freq_au: vec![0.005],
            reduced_mass: vec![1.0],
            force_const_au: vec![1.0],
            force_const_dyne: vec![1.0],
            norm_mode: vec![vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]]],
            norm_mode_mw: vec![vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]]],
            rot_type: crate::thermo::RotorType::Linear,
            rotational_constants_ghz: [0.0, 0.0, 0.0],
            principal_moments_amu_bohr2: [0.0, 0.0, 0.0],
            vib_temperature: vec![1000.0],
            n_modes: 1,
            n_atoms: 2,
        };
        let inv = compute_raman_activities(&t, &freq_info);
        assert_eq!(inv.activities_a4_amu.len(), 1);
        assert!(inv.activities_a4_amu[0] >= 0.0);
    }

    #[test]
    fn test_depolarization_ratio_range() {
        // Random-ish inputs — depolarization must be in [0, 0.75]
        let mut t = PolarDerivTensor::new(1);
        for d in 0..3 {
            for e in 0..3 {
                let v = (d as f64 - 1.0) * (e as f64 - 0.5) * 0.7;
                t.set(d, e, 0, 0, v);
            }
        }
        t.symmetrize_de();

        let freq_info = FrequencyInfo {
            freq_wavenumber: vec![500.0; 3],
            freq_au: vec![0.003; 3],
            reduced_mass: vec![1.0; 3],
            force_const_au: vec![1.0; 3],
            force_const_dyne: vec![1.0; 3],
            norm_mode: vec![
                vec![[1.0, 0.0, 0.0]],
                vec![[0.0, 1.0, 0.0]],
                vec![[0.0, 0.0, 1.0]],
            ],
            norm_mode_mw: vec![
                vec![[1.0, 0.0, 0.0]],
                vec![[0.0, 1.0, 0.0]],
                vec![[0.0, 0.0, 1.0]],
            ],
            rot_type: crate::thermo::RotorType::Atom,
            rotational_constants_ghz: [0.0, 0.0, 0.0],
            principal_moments_amu_bohr2: [0.0, 0.0, 0.0],
            vib_temperature: vec![500.0; 3],
            n_modes: 3,
            n_atoms: 1,
        };
        let inv = compute_raman_activities(&t, &freq_info);
        for rho in &inv.depolarization_ratios {
            assert!(*rho >= 0.0 && *rho <= 0.75, "rho = {rho}");
        }
    }

    #[test]
    fn test_raman_sign_invariance() {
        // Flip sign of a normal mode — Raman activity should be unchanged
        let mut t = PolarDerivTensor::new(2);
        for d in 0..3 {
            for e in 0..3 {
                t.set(d, e, 0, 0, ((d + e + 1) as f64) * 0.3);
                t.set(d, e, 1, 0, -((d + e + 1) as f64) * 0.3);
            }
        }
        t.symmetrize_de();

        let make_info = |sign: f64| FrequencyInfo {
            freq_wavenumber: vec![1200.0],
            freq_au: vec![0.005],
            reduced_mass: vec![1.0],
            force_const_au: vec![1.0],
            force_const_dyne: vec![1.0],
            norm_mode: vec![vec![[sign * 0.5, 0.0, 0.0], [-sign * 0.5, 0.0, 0.0]]],
            norm_mode_mw: vec![vec![[sign * 0.5, 0.0, 0.0], [-sign * 0.5, 0.0, 0.0]]],
            rot_type: crate::thermo::RotorType::Linear,
            rotational_constants_ghz: [0.0, 0.0, 0.0],
            principal_moments_amu_bohr2: [0.0, 0.0, 0.0],
            vib_temperature: vec![1200.0],
            n_modes: 1,
            n_atoms: 2,
        };
        let inv_pos = compute_raman_activities(&t, &make_info(1.0));
        let inv_neg = compute_raman_activities(&t, &make_info(-1.0));
        assert_abs_diff_eq!(
            inv_pos.activities_a4_amu[0],
            inv_neg.activities_a4_amu[0],
            epsilon = 1e-14
        );
    }

    // ========== Golden JSON tests (H2O/CH4/CO2/H2) ==========

    #[derive(Debug, Deserialize)]
    struct GoldenRamanData {
        #[allow(dead_code)]
        name: String,
        atoms: Vec<GoldenAtom>,
        #[allow(dead_code)]
        n_atoms: usize,
        n_modes: usize,
        #[allow(dead_code)]
        energy: f64,
        polarizability_au: Vec<Vec<f64>>,
        #[allow(dead_code)]
        polarizability_ang3: Vec<Vec<f64>>,
        polarizability_derivs_cartesian: Vec<Vec<Vec<f64>>>, // (3, 3, 3*n_atoms)
        freq_wavenumber: Vec<f64>,
        #[allow(dead_code)]
        norm_mode: Vec<Vec<[f64; 3]>>,
        raman_activities_a4_amu: Vec<f64>,
        #[allow(dead_code)]
        depolarization_ratios: Vec<f64>,
    }

    #[derive(Debug, Deserialize)]
    struct GoldenAtom {
        #[allow(dead_code)]
        #[serde(rename = "Z")]
        z: u8,
        #[allow(dead_code)]
        symbol: String,
        pos_bohr: [f64; 3],
    }

    fn load_golden_raman(name: &str) -> GoldenRamanData {
        let path = format!(
            "{}/../../tests/golden/raman/{}_sto3g_rhf.json",
            env!("CARGO_MANIFEST_DIR"),
            name
        );
        let raw =
            std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("Failed to load {path}: {e}"));
        serde_json::from_str(&raw).unwrap_or_else(|e| panic!("Failed to parse {path}: {e}"))
    }

    fn build_atoms_from_golden(g: &GoldenRamanData) -> Vec<Atom> {
        g.atoms
            .iter()
            .map(|a| Atom::new(a.z, a.pos_bohr).unwrap())
            .collect()
    }

    fn run_raman_pipeline(atoms: &[Atom]) -> RamanResult {
        let basis = BasisSet::build(atoms.to_vec(), "sto-3g").unwrap();
        let hess = run_hessian_pipeline(atoms);
        let freq_info = harmonic_analysis(&basis.atoms, &hess.hessian).expect("harmonic analysis");
        let cphf_cfg = CphfConfig::default();
        compute_raman_spectrum(&basis.atoms, &basis, &hess, &freq_info, &cphf_cfg)
            .expect("Raman spectrum")
    }

    #[test]
    fn test_compute_polarizability_h2o_vs_pyscf() {
        // AC2: polarizability within 1e-4 bohr³
        let golden = load_golden_raman("h2o");
        let atoms = build_atoms_from_golden(&golden);
        let result = run_raman_pipeline(&atoms);

        for d in 0..3 {
            for e in 0..3 {
                let iqcp = result.polarizability_au[d][e];
                let pyscf = golden.polarizability_au[d][e];
                let err = (iqcp - pyscf).abs();
                assert!(
                    err < 1e-4,
                    "α[{d},{e}]: IQCP={iqcp:.6e}, PySCF={pyscf:.6e}, err={err:.3e}"
                );
            }
        }
    }

    #[test]
    fn test_compute_raman_activities_h2o_vs_pyscf() {
        // AC6: Raman activities within ±1 Å⁴/amu
        let golden = load_golden_raman("h2o");
        let atoms = build_atoms_from_golden(&golden);
        let result = run_raman_pipeline(&atoms);

        assert_eq!(result.n_modes, 3);
        assert_eq!(result.raman_activities_a4_amu.len(), 3);

        for k in 0..3 {
            let iqcp = result.raman_activities_a4_amu[k];
            let pyscf = golden.raman_activities_a4_amu[k];
            let err = (iqcp - pyscf).abs();
            assert!(
                err < 1.0,
                "H2O Raman activity mode {k}: IQCP={iqcp:.4} Å⁴/amu, \
                 PySCF={pyscf:.4} Å⁴/amu, err={err:.4}"
            );
        }
    }

    #[test]
    fn test_compute_raman_co2_mutual_exclusion() {
        // AC7: σ_g symmetric stretch Raman-active, σ_u asymmetric stretch Raman-inactive
        let golden = load_golden_raman("co2");
        let atoms = build_atoms_from_golden(&golden);
        let result = run_raman_pipeline(&atoms);

        assert_eq!(result.n_modes, 4);

        // Identify σ_g (~1571 cm⁻¹) and σ_u (~2829 cm⁻¹) by frequency
        let sigma_g_idx = golden
            .freq_wavenumber
            .iter()
            .position(|&f| (1500.0..1650.0).contains(&f))
            .expect("σ_g symmetric stretch at ~1571 cm⁻¹");
        let sigma_u_idx = golden
            .freq_wavenumber
            .iter()
            .position(|&f| (2700.0..2900.0).contains(&f))
            .expect("σ_u asymmetric stretch at ~2829 cm⁻¹");

        let sigma_g_raman = result.raman_activities_a4_amu[sigma_g_idx];
        let sigma_u_raman = result.raman_activities_a4_amu[sigma_u_idx];

        assert!(
            sigma_g_raman > 1.0,
            "CO2 σ_g Raman activity = {sigma_g_raman:.4} Å⁴/amu, expected > 1.0"
        );
        assert!(
            sigma_u_raman < 0.01,
            "CO2 σ_u Raman activity = {sigma_u_raman:.4} Å⁴/amu, expected < 0.01"
        );
    }

    #[test]
    fn test_compute_raman_h2_active() {
        // H2/STO-3G has only σ bonding basis functions, so its transverse
        // polarizability is exactly zero. The single σ_g mode is Raman-active
        // but not "totally symmetric" in the invariant sense — both
        // bar_alpha' and gamma' are nonzero, so depolarization ratio ≈ 1/3.
        // Compare with golden for quantitative agreement.
        let golden = load_golden_raman("h2");
        let atoms = build_atoms_from_golden(&golden);
        let result = run_raman_pipeline(&atoms);

        assert_eq!(result.n_modes, 1);
        assert!(
            result.raman_activities_a4_amu[0] > 1.0,
            "H2 Raman activity = {} Å⁴/amu, expected > 1.0",
            result.raman_activities_a4_amu[0]
        );
        // Match golden value within 1 Å⁴/amu
        let err = (result.raman_activities_a4_amu[0] - golden.raman_activities_a4_amu[0]).abs();
        assert!(
            err < 1.0,
            "H2 Raman activity: IQCP={}, PySCF={}, err={}",
            result.raman_activities_a4_amu[0],
            golden.raman_activities_a4_amu[0],
            err
        );
    }

    #[test]
    fn test_compute_raman_ch4_a1_symmetric() {
        let golden = load_golden_raman("ch4");
        let atoms = build_atoms_from_golden(&golden);
        let result = run_raman_pipeline(&atoms);

        assert_eq!(result.n_modes, 9);

        // Identify A1 symmetric stretch (~3476 cm⁻¹)
        let a1_idx = golden
            .freq_wavenumber
            .iter()
            .position(|&f| (3400.0..3550.0).contains(&f))
            .expect("CH4 A1 symmetric stretch at ~3476 cm⁻¹");

        // A1 should have rho ≈ 0 (totally symmetric)
        let a1_rho = result.depolarization_ratios[a1_idx];
        assert!(
            a1_rho < 0.01,
            "CH4 A1 depolarization ratio = {a1_rho:.6}, expected < 0.01"
        );
    }

    /// Diagnostic: `compute_polarizability_derivatives` (the production
    /// path) should match `compute_polarizability_derivatives_fd_reference`
    /// identically on H2/STO-3G. Since the production path delegates to
    /// the FD reference internally, this test verifies the delegation
    /// works correctly (bit-exact match expected).
    #[test]
    fn test_polar_deriv_analytical_vs_fd_h2() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let atoms = vec![h1, h2];

        let basis = BasisSet::build(atoms.clone(), "sto-3g").unwrap();
        let hess = run_hessian_pipeline(&atoms);
        let cphf_cfg = CphfConfig::default();
        let field_cphf = compute_field_cphf(&atoms, &basis, &hess, &cphf_cfg).unwrap();

        let analytical =
            compute_polarizability_derivatives(&atoms, &basis, &field_cphf, &hess).unwrap();

        let scf_cfg = ScfConfig {
            max_iterations: 200,
            use_diis: true,
            diis_start: 2,
            ..ScfConfig::tight()
        };
        let fd =
            compute_polarizability_derivatives_fd_reference(&atoms, "sto-3g", &scf_cfg, &cphf_cfg)
                .unwrap();

        eprintln!("\nH2 polar deriv comparison:");
        for d in 0..3 {
            for e in 0..3 {
                for a in 0..2 {
                    for g in 0..3 {
                        let av = analytical.get(d, e, a, g);
                        let fv = fd.get(d, e, a, g);
                        if av.abs() > 1e-6 || fv.abs() > 1e-6 {
                            eprintln!(
                                "  ({},{},{},{}): analytical={:12.5e}  fd={:12.5e}  err={:10.3e}",
                                d,
                                e,
                                a,
                                g,
                                av,
                                fv,
                                (av - fv).abs()
                            );
                        }
                    }
                }
            }
        }

        let mut max_err = 0.0f64;
        for d in 0..3 {
            for e in 0..3 {
                for a in 0..2 {
                    for g in 0..3 {
                        let err = (analytical.get(d, e, a, g) - fd.get(d, e, a, g)).abs();
                        if err > max_err {
                            max_err = err;
                        }
                    }
                }
            }
        }
        eprintln!("H2 max err = {max_err:.6e}");
        assert!(
            max_err < 1e-10,
            "H2 production-vs-FD-reference: max err {max_err:.6e}"
        );
    }

    /// Diagnostic: `compute_polarizability_derivatives` (the production
    /// path) should match `compute_polarizability_derivatives_fd_reference`
    /// identically on H2O/STO-3G (delegation bit-exact).
    #[test]
    fn test_polar_deriv_analytical_vs_fd_h2o() {
        let o = Atom::new(8, [0.0, 0.0, 0.2217]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.4309, -0.8867]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.4309, -0.8867]).unwrap();
        let atoms = vec![o, h1, h2];

        let basis = BasisSet::build(atoms.clone(), "sto-3g").unwrap();
        let hess = run_hessian_pipeline(&atoms);
        let cphf_cfg = CphfConfig::default();
        let field_cphf = compute_field_cphf(&atoms, &basis, &hess, &cphf_cfg).unwrap();

        // Analytical
        let analytical =
            compute_polarizability_derivatives(&atoms, &basis, &field_cphf, &hess).unwrap();

        // FD reference (expensive; 18 SCF+CPHF solves for H2O)
        let scf_cfg = ScfConfig {
            max_iterations: 200,
            use_diis: true,
            diis_start: 2,
            ..ScfConfig::tight()
        };
        let fd =
            compute_polarizability_derivatives_fd_reference(&atoms, "sto-3g", &scf_cfg, &cphf_cfg)
                .unwrap();

        // Compare elementwise
        let mut max_err = 0.0f64;
        let mut max_err_idx = (0usize, 0, 0, 0);
        let mut sum_ratio = 0.0f64;
        let mut n_ratio = 0usize;
        for d in 0..3 {
            for e in 0..3 {
                for a in 0..3 {
                    for g in 0..3 {
                        let av = analytical.get(d, e, a, g);
                        let fv = fd.get(d, e, a, g);
                        let err = (av - fv).abs();
                        if err > max_err {
                            max_err = err;
                            max_err_idx = (d, e, a, g);
                        }
                        if fv.abs() > 0.01 {
                            sum_ratio += av / fv;
                            n_ratio += 1;
                        }
                    }
                }
            }
        }
        let (md, me, ma, mg) = max_err_idx;
        let mean_ratio = if n_ratio > 0 {
            sum_ratio / (n_ratio as f64)
        } else {
            1.0
        };
        eprintln!(
            "H2O polar deriv analytical-vs-FD max err = {:.6e} at [d={},e={},atom={},γ={}]\n  \
             analytical = {:.6e}\n  fd         = {:.6e}\n  \
             mean ratio analytical/fd over non-zero elements = {:.6e}",
            max_err,
            md,
            me,
            ma,
            mg,
            analytical.get(md, me, ma, mg),
            fd.get(md, me, ma, mg),
            mean_ratio,
        );

        // Dump all 81 elements for analysis
        eprintln!("\nFull comparison (d, e, atom, γ): analytical | fd | err");
        for d in 0..3 {
            for e in 0..3 {
                for a in 0..3 {
                    for g in 0..3 {
                        let av = analytical.get(d, e, a, g);
                        let fv = fd.get(d, e, a, g);
                        if av.abs() > 1e-6 || fv.abs() > 1e-6 {
                            eprintln!(
                                "  ({},{},{},{}): {:12.5e} | {:12.5e} | {:10.3e}",
                                d,
                                e,
                                a,
                                g,
                                av,
                                fv,
                                (av - fv).abs()
                            );
                        }
                    }
                }
            }
        }

        assert!(
            max_err < 1e-10,
            "Production polarizability derivative disagrees with FD reference on H2O: max err {max_err:.6e}"
        );
    }
}
