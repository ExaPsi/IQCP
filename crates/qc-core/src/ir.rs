//! Infrared absorption intensities via analytical dipole derivatives
//!
//! This module implements the top-level IR pipeline for US-097:
//!
//! 1. **Equilibrium dipole moment** — [`compute_dipole_moment`] evaluates
//!    `μ_d = sum_A Z_A * (R_A,d - O_d) - tr(D * dipole^d)` matching
//!    PySCF `scf/hf.py:dip_moment()` (lines 1371-1428).
//!
//! 2. **Analytical dipole derivatives** — [`compute_dipole_derivatives`]
//!    builds the 3 × (3*n_atoms) matrix `∂μ_d/∂R_{A,e}` by assembling
//!    three contributions:
//!    - **Nuclear**: `Z_A * δ_{d,e}`
//!    - **Skeleton**: `-tr(D * ∂dipole^d/∂R_{A,e})` using
//!      [`shell_dipole_first_deriv`](crate::integrals::shell_dipole_first_deriv)
//!    - **CPHF response**: `-4 * sum_{pi} U^(R_{A,e})_{pi} * dipole^d_MO[p,i]`
//!      using `HessianResult.mo1_cphf` (reused, no re-solve of CPHF)
//!
//! 3. **Normal mode projection + km/mol** — [`compute_ir_intensities`]
//!    projects the Cartesian dipole derivatives onto the `FrequencyInfo`
//!    normal modes and applies the [`AU_TO_KMMOL`](crate::constants::AU_TO_KMMOL)
//!    conversion factor.
//!
//! The output is packaged in [`IrSpectrumResult`] for downstream US-100
//! (simulated IR spectra) and WASM export.
//!
//! # References
//!
//! - Helgaker, Jørgensen & Olsen (2000) Chapters 9 & 13
//! - Obara & Saika (1986) J. Chem. Phys. 84, 3963
//! - Pulay (1969) Mol. Phys. 17, 197
//! - PySCF `scf/hf.py:dip_moment` (line 1371), `hessian/rhf.py:solve_mo1` (line 286)
//! - pyscf-forge `prop/infrared/rhf.py:kernel_dipderiv` and `kernel_ir`
//! - Person, W.B. & Zerbi, G. (1982) *Vibrational Intensities in
//!   Infrared and Raman Spectroscopy*

use crate::basis::{Atom, BasisSet};
use crate::constants::{AU_TO_DEBYE, AU_TO_KMMOL};
use crate::integrals::{dipole_matrix, shell_dipole_first_deriv};
use crate::scf::hessian::{CphfMo1Data, HessianResult};
use crate::thermo::FrequencyInfo;
use nalgebra::DMatrix;
use thiserror::Error;

// =============================================================================
// Errors
// =============================================================================

/// Errors that can occur during IR intensity computation.
#[derive(Debug, Error, PartialEq)]
pub enum IrError {
    /// `HessianResult.mo1_cphf` was `None`; cannot compute dipole derivatives
    /// without the CPHF solution.
    #[error(
        "HessianResult.mo1_cphf is None; compute_dipole_derivatives requires \
         the CPHF solution from rhf_hessian or dft_hessian"
    )]
    NoCphfData,

    /// Number of atoms inferred from `FrequencyInfo` or `HessianResult` does
    /// not match the number of atoms in the basis set.
    #[error("Size mismatch: expected {expected} atoms, got {actual} atoms")]
    SizeMismatch { expected: usize, actual: usize },

    /// A non-finite value was encountered in the dipole derivative or
    /// intensity computation.
    #[error("Non-finite result detected in {stage}")]
    NonFiniteResult { stage: &'static str },
}

// =============================================================================
// Equilibrium dipole moment
// =============================================================================

/// Compute the equilibrium electric dipole moment of a closed-shell
/// molecule.
///
/// # Formula
///
/// Matches PySCF `scf/hf.py:dip_moment()` lines 1419-1421 exactly:
///
/// ```text
/// mu_d = sum_A Z_A * (R_{A,d} - O_d) - tr(D * dipole^d)
/// ```
///
/// where `O` is the gauge origin (default `[0, 0, 0]`), `Z_A` is the
/// nuclear charge of atom `A`, `R_{A,d}` is the `d`-th coordinate of
/// atom `A` (bohr), `D` is the closed-shell density matrix (with the
/// factor of 2 for double occupancy baked in), and `dipole^d` is the
/// AO-basis dipole integral matrix for direction `d`.
///
/// # Arguments
///
/// * `atoms` - Vector of atoms (atomic number + bohr position)
/// * `basis` - Molecular basis set
/// * `density` - Closed-shell density matrix `(nbf, nbf)`
/// * `gauge_origin` - Gauge origin `C` for the dipole operator (default
///   `[0, 0, 0]` matches PySCF)
///
/// # Returns
///
/// `[mu_x, mu_y, mu_z]` in atomic units (electron·bohr).
///
/// Multiply by [`AU_TO_DEBYE`](crate::constants::AU_TO_DEBYE) to convert to Debye.
///
/// # Reference
///
/// - PySCF `scf/hf.py:1371-1428` (`dip_moment`)
/// - PySCF `scf/hf.py:1418` (`intor_symmetric('int1e_r', comp=3)`)
/// - Helgaker, Jørgensen & Olsen (2000), Eq. 13.3.2
pub fn compute_dipole_moment(
    atoms: &[Atom],
    basis: &BasisSet,
    density: &DMatrix<f64>,
    gauge_origin: &[f64; 3],
) -> [f64; 3] {
    let dip_ao = dipole_matrix(basis, gauge_origin);
    let mut mu = [0.0f64; 3];
    for d in 0..3 {
        // Nuclear: sum_A Z_A * (R_{A,d} - O_d)
        let nuclear: f64 = atoms
            .iter()
            .map(|a| (a.atomic_number as f64) * (a.position[d] - gauge_origin[d]))
            .sum();
        // Electronic: -tr(D * dipole^d) = -sum_{i,j} D[i,j] * dipole^d[i,j]
        let mut electronic = 0.0f64;
        for i in 0..density.nrows() {
            for j in 0..density.ncols() {
                electronic += density[(i, j)] * dip_ao[d][(i, j)];
            }
        }
        mu[d] = nuclear - electronic;
    }
    mu
}

// =============================================================================
// Analytical dipole derivatives
// =============================================================================

/// Compute the analytical dipole derivative matrix `∂μ_d/∂R_{A,e}` in
/// atomic units, reusing the CPHF solution from a prior Hessian run.
///
/// # Formula
///
/// ```text
/// ∂μ_d/∂R_{A,e} =
///     Z_A * δ_{d,e}                              [nuclear]
///   - sum_{μν} D_{μν} * ∂dipole^d_{μν}/∂R_{A,e}  [skeleton, from shell_dipole_first_deriv]
///   - 4 * sum_{pi} U^(R_{A,e})_{pi} * dipole^d_MO[p,i]   [CPHF response]
/// ```
///
/// The CPHF term uses the identity (for a closed-shell density and a
/// real basis):
///
/// ```text
/// -tr(∂D/∂R · μ^d) = -4 * sum_{pi} mo1[p,i] * (C^T μ^d_AO C)[p,i]
/// ```
///
/// where `p` runs over all MOs and `i` over occupied MOs. The factor
/// of 4 comes from:
/// - Factor 2 for closed-shell double occupancy
/// - Factor 2 for the bra/ket symmetry of the density-response trace
///   over a symmetric dipole integral matrix
///
/// The occupied-occupied block of `mo1` is set to `-s1_MO/2` by
/// `cphf_solve_withs1` (enforced by the orthonormality constraint),
/// so this formula automatically includes the `s^(1)/2` correction.
///
/// # Arguments
///
/// * `atoms` - Vector of atoms (atomic number + bohr position)
/// * `basis` - Molecular basis set
/// * `density` - Closed-shell density matrix `(nbf, nbf)`
/// * `hess_result` - Result of `rhf_hessian` or `dft_hessian`, which
///   must have `mo1_cphf.is_some()`
/// * `gauge_origin` - Gauge origin `C` for the dipole operator
///
/// # Returns
///
/// A `(3, 3*n_atoms)` matrix with layout
/// `result[(d, 3*A + e)] = ∂μ_d/∂R_{A,e}` in atomic units (`e` = elementary charge).
///
/// # Errors
///
/// * [`IrError::NoCphfData`] if `hess_result.mo1_cphf.is_none()`
/// * [`IrError::SizeMismatch`] if the number of atoms in
///   `hess_result` doesn't match `atoms.len()`
/// * [`IrError::NonFiniteResult`] if the computation produces NaN or Inf
///
/// # Reference
///
/// - pyscf-forge `prop/infrared/rhf.py:kernel_dipderiv` (function body)
/// - Helgaker, Jørgensen & Olsen (2000), Section 13.3, Eq. 11.1.11
pub fn compute_dipole_derivatives(
    atoms: &[Atom],
    basis: &BasisSet,
    density: &DMatrix<f64>,
    hess_result: &HessianResult,
    gauge_origin: &[f64; 3],
) -> Result<DMatrix<f64>, IrError> {
    let n_atoms = atoms.len();
    if hess_result.n_atoms != n_atoms {
        return Err(IrError::SizeMismatch {
            expected: n_atoms,
            actual: hess_result.n_atoms,
        });
    }
    let cphf_data = hess_result.mo1_cphf.as_ref().ok_or(IrError::NoCphfData)?;

    let nbf = basis.n_basis;
    let mut result = DMatrix::<f64>::zeros(3, 3 * n_atoms);

    // ------------------------------------------------------------------
    // (1) Nuclear contribution: Z_A * delta_{d,e}
    // ------------------------------------------------------------------
    for (a, atom) in atoms.iter().enumerate() {
        let z_a = atom.atomic_number as f64;
        for d in 0..3 {
            result[(d, 3 * a + d)] += z_a;
        }
    }

    // ------------------------------------------------------------------
    // (2) Skeleton electronic contribution:
    //     -sum_{μν} D_{μν} * d(dipole^d_{μν})/dR_{A,e}
    //
    // For dipole integrals, translational invariance is NOT the simple
    // `d/dR_j = -d/dR_i` that holds for overlap/kinetic/nuclear
    // integrals. The dipole operator `(r_d - C_d)` explicitly depends
    // on position, so:
    //
    //     d/dA_e <mu|r_d - C_d|nu> + d/dB_e <mu|r_d - C_d|nu>
    //         = δ_{d,e} <mu|nu>
    //
    // (derived from the overlap-shift identity by differentiating both
    // sides and using d/dA_e (A_d - C_d) = δ_{d,e}).
    //
    // To correctly compute the ket-center derivative, we call
    // `shell_dipole_first_deriv(shell_j, shell_i, gauge_origin)` which
    // returns `d/dR_{atom_j, e} <nu|r_d - C_d|mu>`, then transpose the
    // (ia, jb) indices to match our (shell_i, shell_j) block layout.
    // ------------------------------------------------------------------
    let mut shell_offsets = Vec::with_capacity(basis.shells.len());
    let mut acc = 0usize;
    for shell in &basis.shells {
        shell_offsets.push(acc);
        acc += shell.n_basis_functions();
    }

    for (si, shell_i) in basis.shells.iter().enumerate() {
        let off_i = shell_offsets[si];
        let n_i = shell_i.n_basis_functions();
        let atom_i = shell_i.atom_idx;

        for (sj, shell_j) in basis.shells.iter().enumerate() {
            let off_j = shell_offsets[sj];
            let n_j = shell_j.n_basis_functions();
            let atom_j = shell_j.atom_idx;

            // Bra-center derivative of the dipole integral:
            // bra_deriv[d][e][ia*n_j + jb] = d/dR_{atom_i, e} <mu_ia | r_d - C_d | nu_jb>
            let bra_deriv = shell_dipole_first_deriv(shell_i, shell_j, gauge_origin);

            // Ket-center derivative: compute via the reversed shell order.
            // shell_dipole_first_deriv(shell_j, shell_i) returns
            //     d/dR_{atom_j, e} <nu | r_d - C_d | mu>
            // with shape (n_j, n_i). To match our (ia, jb) block we need
            // (transpose) the (ia, jb) indices: ket_deriv[d][e][ia * n_j + jb] =
            //     d/dR_{atom_j, e} <nu_jb | r_d - C_d | mu_ia>.
            // Since the dipole operator is Hermitian on a real basis,
            //     <nu_jb | r_d - C_d | mu_ia> = <mu_ia | r_d - C_d | nu_jb>
            // so the TRANSPOSED ket_deriv equals d/dR_{atom_j, e}
            // applied to the original (ia, jb) bra-ket labeling.
            let ket_deriv_raw = shell_dipole_first_deriv(shell_j, shell_i, gauge_origin);
            // ket_deriv_raw[d][e][jb * n_i + ia] = d/dR_{atom_j, e} <nu_jb | r_d - C_d | mu_ia>

            for d in 0..3 {
                for e in 0..3 {
                    let mut bra_acc = 0.0f64;
                    let mut ket_acc = 0.0f64;
                    for ia in 0..n_i {
                        for jb in 0..n_j {
                            let d_ij = density[(off_i + ia, off_j + jb)];
                            bra_acc += d_ij * bra_deriv[d][e][ia * n_j + jb];
                            // ket_deriv_raw has (jb, ia) layout for the reversed call
                            ket_acc += d_ij * ket_deriv_raw[d][e][jb * n_i + ia];
                        }
                    }
                    // Electronic skeleton contribution: -D * d(dipole)/dR
                    result[(d, 3 * atom_i + e)] -= bra_acc;
                    result[(d, 3 * atom_j + e)] -= ket_acc;
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // (3) CPHF response contribution:
    //     -4 * sum_{p,i} mo1[k][p,i] * dipole^d_MO[p,i]
    //
    // where k = 3*A + e and dipole^d_MO[p,i] = (C^T * dipole^d_AO * C)[p,i]
    // (sliced to occupied columns 0..n_occ for i).
    // ------------------------------------------------------------------
    let dip_ao = dipole_matrix(basis, gauge_origin);
    let c = &cphf_data.mo_coeff;
    let n_occ = cphf_data.n_occ;

    // Pre-compute the MO-basis dipole matrices sliced to (nmo × n_occ)
    let dip_mo: Vec<DMatrix<f64>> = (0..3)
        .map(|d| {
            let ct_dip_c = c.transpose() * &dip_ao[d] * c;
            ct_dip_c.columns(0, n_occ).clone_owned()
        })
        .collect();

    for a in 0..n_atoms {
        for e in 0..3 {
            let k = 3 * a + e;
            let mo1_k = &cphf_data.mo1[k];
            // mo1_k has shape (nmo, n_occ). Contract with dip_mo[d] element-wise.
            for d in 0..3 {
                let mut sum = 0.0f64;
                // dip_mo[d] is also (nmo, n_occ).
                for p in 0..mo1_k.nrows() {
                    for i in 0..n_occ {
                        sum += mo1_k[(p, i)] * dip_mo[d][(p, i)];
                    }
                }
                result[(d, 3 * a + e)] -= 4.0 * sum;
            }
        }
    }

    // ------------------------------------------------------------------
    // Sanity check: non-finite values
    // ------------------------------------------------------------------
    for val in result.iter() {
        if !val.is_finite() {
            return Err(IrError::NonFiniteResult {
                stage: "compute_dipole_derivatives",
            });
        }
    }
    // Also sanity check nbf matches
    debug_assert_eq!(density.nrows(), nbf);
    debug_assert_eq!(density.ncols(), nbf);

    Ok(result)
}

// =============================================================================
// IR intensities
// =============================================================================

/// Package of IR spectrum results for a molecule.
///
/// Contains all quantities needed for downstream US-100 spectrum
/// visualization (broadened Lorentzian/Gaussian plots) as well as raw
/// per-mode dipole derivative information for pedagogical inspection.
#[derive(Debug, Clone)]
pub struct IrSpectrumResult {
    /// Equilibrium dipole moment in atomic units (`e*bohr`).
    pub equilibrium_dipole_au: [f64; 3],
    /// Equilibrium dipole moment in Debye (for display).
    pub equilibrium_dipole_debye: [f64; 3],
    /// Cartesian dipole derivatives, shape `(3, 3*n_atoms)`,
    /// layout `[d, 3*A + e] = ∂μ_d/∂R_{A,e}` in au/bohr.
    pub dipole_derivs_cartesian: DMatrix<f64>,
    /// Per-mode normal-coordinate dipole derivatives: one `[f64; 3]`
    /// vector per mode, in `e/sqrt(amu)` (atomic units).
    pub dipole_derivs_normal_mode: Vec<[f64; 3]>,
    /// Integrated IR intensities in km/mol, one value per mode.
    /// Ordering matches `freq_info.freq_wavenumber`.
    pub intensities_km_per_mol: Vec<f64>,
    /// Number of vibrational modes (`3N-6` or `3N-5` depending on rotor
    /// type; `0` for atoms).
    pub n_modes: usize,
}

/// Compute per-mode IR absorption intensities from Cartesian dipole
/// derivatives and vibrational normal modes.
///
/// # Formula
///
/// For each vibrational mode `k`:
///
/// ```text
/// ∂μ/∂Q_k[d] = sum_{A, e} (∂μ_d/∂R_{A,e}) * norm_mode[k][A][e]
/// A_k [km/mol] = AU_TO_KMMOL * (|∂μ/∂Q_k|^2)
/// ```
///
/// where `|∂μ/∂Q_k|^2 = ∂μ_x^2 + ∂μ_y^2 + ∂μ_z^2` is the squared
/// magnitude and [`AU_TO_KMMOL`](crate::constants::AU_TO_KMMOL) is the
/// conversion factor (~974.88 km/mol per `e^2/amu`).
///
/// The sign ambiguity of normal modes (PySCF vs IQCP may differ in the
/// sign of `norm_mode`) is irrelevant because the intensity depends on
/// the squared magnitude.
///
/// # Arguments
///
/// * `dipole_derivs` - `(3, 3*n_atoms)` matrix in au/bohr (`e`)
/// * `freq_info` - `FrequencyInfo` from `harmonic_analysis` with
///   `norm_mode` field
///
/// # Returns
///
/// Ordered list of IR intensities (one per mode) in km/mol, plus
/// per-mode `∂μ/∂Q` vectors (in `e/sqrt(amu)`, atomic units).
///
/// # Reference
///
/// - pyscf-forge `prop/infrared/rhf.py:kernel_ir` (line 161)
/// - Psi4 `driver/qcdb/vib.py` line 589
/// - Person & Zerbi (1982), *Vibrational Intensities*
pub fn compute_ir_intensities(
    dipole_derivs: &DMatrix<f64>,
    freq_info: &FrequencyInfo,
) -> (Vec<f64>, Vec<[f64; 3]>) {
    let n_modes = freq_info.n_modes;
    let mut intensities = Vec::with_capacity(n_modes);
    let mut dmu_dq = Vec::with_capacity(n_modes);

    for k in 0..n_modes {
        let mut dqk = [0.0f64; 3];
        for (a, atom_displacement) in freq_info.norm_mode[k].iter().enumerate() {
            for e in 0..3 {
                for (d, dqk_d) in dqk.iter_mut().enumerate() {
                    *dqk_d += dipole_derivs[(d, 3 * a + e)] * atom_displacement[e];
                }
            }
        }
        let mag_sq = dqk[0] * dqk[0] + dqk[1] * dqk[1] + dqk[2] * dqk[2];
        intensities.push(mag_sq * AU_TO_KMMOL);
        dmu_dq.push(dqk);
    }

    (intensities, dmu_dq)
}

/// Top-level convenience wrapper: compute the full IR spectrum
/// (equilibrium dipole + derivatives + intensities) in a single call.
///
/// This is the preferred entry point for US-097 downstream consumers
/// (US-100 simulated spectra, WASM export).
///
/// # Arguments
///
/// * `atoms` - Atom list (`n_atoms`)
/// * `basis` - Basis set
/// * `density` - Closed-shell density matrix `(nbf, nbf)`
/// * `hess_result` - Hessian result with `mo1_cphf.is_some()`
/// * `freq_info` - Normal mode analysis output
/// * `gauge_origin` - Gauge origin for the dipole operator (default
///   `[0, 0, 0]`)
///
/// # Errors
///
/// Propagates [`IrError`] from [`compute_dipole_derivatives`] and
/// returns [`IrError::SizeMismatch`] if `freq_info.n_atoms != atoms.len()`.
pub fn compute_ir_spectrum(
    atoms: &[Atom],
    basis: &BasisSet,
    density: &DMatrix<f64>,
    hess_result: &HessianResult,
    freq_info: &FrequencyInfo,
    gauge_origin: &[f64; 3],
) -> Result<IrSpectrumResult, IrError> {
    if freq_info.n_atoms != atoms.len() {
        return Err(IrError::SizeMismatch {
            expected: atoms.len(),
            actual: freq_info.n_atoms,
        });
    }

    let dipole_au = compute_dipole_moment(atoms, basis, density, gauge_origin);
    let dipole_debye = [
        dipole_au[0] * AU_TO_DEBYE,
        dipole_au[1] * AU_TO_DEBYE,
        dipole_au[2] * AU_TO_DEBYE,
    ];

    // Handle atom edge case: no modes, no derivatives to compute
    if freq_info.n_modes == 0 {
        return Ok(IrSpectrumResult {
            equilibrium_dipole_au: dipole_au,
            equilibrium_dipole_debye: dipole_debye,
            dipole_derivs_cartesian: DMatrix::<f64>::zeros(3, 3 * atoms.len()),
            dipole_derivs_normal_mode: Vec::new(),
            intensities_km_per_mol: Vec::new(),
            n_modes: 0,
        });
    }

    let dipole_derivs =
        compute_dipole_derivatives(atoms, basis, density, hess_result, gauge_origin)?;
    let (intensities, dmu_dq) = compute_ir_intensities(&dipole_derivs, freq_info);

    Ok(IrSpectrumResult {
        equilibrium_dipole_au: dipole_au,
        equilibrium_dipole_debye: dipole_debye,
        dipole_derivs_cartesian: dipole_derivs,
        dipole_derivs_normal_mode: dmu_dq,
        intensities_km_per_mol: intensities,
        n_modes: freq_info.n_modes,
    })
}

// Suppress unused-import warnings in non-test compilation
#[allow(dead_code)]
const _: fn() = || {
    let _ = CphfMo1Data {
        mo1: Vec::new(),
        mo_e1: None,
        mo_coeff: DMatrix::<f64>::zeros(0, 0),
        mo_energies: Vec::new(),
        n_occ: 0,
    };
};

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
    use crate::scf::{rhf_scf_with_guess, ScfConfig};
    use crate::thermo::harmonic_analysis;
    use approx::assert_abs_diff_eq;
    use serde::Deserialize;

    // ========== Helpers ==========

    #[derive(Debug, Deserialize)]
    struct GoldenIrData {
        #[allow(dead_code)]
        name: String,
        atoms: Vec<GoldenAtom>,
        #[allow(dead_code)]
        n_atoms: usize,
        n_modes: usize,
        #[allow(dead_code)]
        energy: f64,
        dipole_au: Vec<f64>,
        dipole_debye: Vec<f64>,
        freq_wavenumber: Vec<f64>,
        #[allow(dead_code)]
        norm_mode: Vec<Vec<[f64; 3]>>,
        dipole_derivs_cartesian: Vec<Vec<f64>>,
        #[allow(dead_code)]
        #[serde(rename = "dmu_dQ")]
        dmu_d_q: Vec<[f64; 3]>,
        ir_intensities_km_mol: Vec<f64>,
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

    fn load_golden_ir(name: &str) -> GoldenIrData {
        let path = format!(
            "{}/../../tests/golden/ir/{}_sto3g_rhf.json",
            env!("CARGO_MANIFEST_DIR"),
            name
        );
        let raw =
            std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("Failed to load {path}: {e}"));
        serde_json::from_str(&raw).unwrap_or_else(|e| panic!("Failed to parse {path}: {e}"))
    }

    fn build_atoms(golden: &GoldenIrData) -> Vec<Atom> {
        golden
            .atoms
            .iter()
            .map(|a| Atom::new(a.z, a.pos_bohr).unwrap())
            .collect()
    }

    fn run_rhf_for_system(atoms: &[Atom]) -> (BasisSet, DMatrix<f64>) {
        use crate::integrals;
        let basis = BasisSet::build(atoms.to_vec(), "sto-3g").unwrap();
        let nbf = basis.n_basis;
        let sys = crate::scf::PresetSystem {
            system_id: "ir_test".to_string(),
            label: "IR test system".to_string(),
            nbf,
            nelec: basis.n_electrons,
            s_matrix: integrals::overlap_matrix(&basis),
            h_core: integrals::hcore_matrix(&basis),
            eri_compressed: integrals::eri_compressed(&basis),
            e_nuc: basis.nuclear_repulsion,
        };
        let scf_config = ScfConfig {
            max_iterations: 200,
            use_diis: true,
            diis_start: 2,
            ..ScfConfig::tight()
        };
        let sad = crate::scf::sad::build_sad_density(&basis);
        let scf = rhf_scf_with_guess(&sys, &scf_config, Some(&sad)).unwrap();
        let density = DMatrix::from_column_slice(nbf, nbf, &scf.density_matrix);
        (basis, density)
    }

    // ========== Unit tests ==========

    #[test]
    fn test_compute_dipole_moment_h2_zero() {
        // H2 homonuclear: dipole moment must be exactly zero by inversion symmetry
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let atoms = vec![h1, h2];
        let (basis, density) = run_rhf_for_system(&atoms);
        let mu = compute_dipole_moment(&atoms, &basis, &density, &[0.0, 0.0, 0.0]);
        // By symmetry all components = 0
        for d in 0..3 {
            assert!(
                mu[d].abs() < 1e-10,
                "H2 dipole component {d} = {}, expected 0",
                mu[d]
            );
        }
    }

    #[test]
    fn test_compute_dipole_moment_h2o_vs_pyscf() {
        // AC2: H2O dipole moment must match PySCF within 1e-4 Debye
        let golden = load_golden_ir("h2o");
        let atoms = build_atoms(&golden);
        let (basis, density) = run_rhf_for_system(&atoms);
        let mu_au = compute_dipole_moment(&atoms, &basis, &density, &[0.0, 0.0, 0.0]);
        let mu_debye = [
            mu_au[0] * AU_TO_DEBYE,
            mu_au[1] * AU_TO_DEBYE,
            mu_au[2] * AU_TO_DEBYE,
        ];
        for d in 0..3 {
            let err = (mu_debye[d] - golden.dipole_debye[d]).abs();
            assert!(
                err < 1e-4,
                "H2O dipole[{d}] = {} Debye, PySCF = {} Debye, err = {}",
                mu_debye[d],
                golden.dipole_debye[d],
                err
            );
        }
    }

    #[test]
    fn test_compute_dipole_moment_ch4_zero() {
        // CH4 T_d symmetry: dipole moment = 0
        let golden = load_golden_ir("ch4");
        let atoms = build_atoms(&golden);
        let (basis, density) = run_rhf_for_system(&atoms);
        let mu = compute_dipole_moment(&atoms, &basis, &density, &[0.0, 0.0, 0.0]);
        for d in 0..3 {
            assert!(
                mu[d].abs() < 1e-9,
                "CH4 dipole component {d} = {}, expected 0",
                mu[d]
            );
        }
    }

    #[test]
    fn test_compute_dipole_moment_co2_zero() {
        // CO2 D∞h symmetry: dipole moment = 0
        let golden = load_golden_ir("co2");
        let atoms = build_atoms(&golden);
        let (basis, density) = run_rhf_for_system(&atoms);
        let mu = compute_dipole_moment(&atoms, &basis, &density, &[0.0, 0.0, 0.0]);
        for d in 0..3 {
            assert!(
                mu[d].abs() < 1e-9,
                "CO2 dipole component {d} = {}, expected 0",
                mu[d]
            );
        }
    }

    // ========== Error path tests ==========

    #[test]
    fn test_compute_dipole_derivatives_missing_cphf() {
        // Build a mock HessianResult with mo1_cphf = None and confirm
        // IrError::NoCphfData is returned.
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let atoms = vec![h1, h2];
        let (basis, density) = run_rhf_for_system(&atoms);
        let fake_hess = HessianResult {
            hessian: DMatrix::<f64>::zeros(6, 6),
            energy: 0.0,
            n_atoms: 2,
            cphf_iterations: 0,
            cphf_converged: true,
            mo1_cphf: None,
        };
        let err =
            compute_dipole_derivatives(&atoms, &basis, &density, &fake_hess, &[0.0, 0.0, 0.0])
                .unwrap_err();
        assert_eq!(err, IrError::NoCphfData);
    }

    #[test]
    fn test_compute_dipole_derivatives_size_mismatch() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let atoms = vec![h1, h2];
        let (basis, density) = run_rhf_for_system(&atoms);
        // Mock: hessian reports wrong n_atoms
        let fake_hess = HessianResult {
            hessian: DMatrix::<f64>::zeros(9, 9),
            energy: 0.0,
            n_atoms: 3, // wrong!
            cphf_iterations: 0,
            cphf_converged: true,
            mo1_cphf: None,
        };
        let err =
            compute_dipole_derivatives(&atoms, &basis, &density, &fake_hess, &[0.0, 0.0, 0.0])
                .unwrap_err();
        assert!(matches!(err, IrError::SizeMismatch { .. }));
    }

    // ========== Full pipeline golden tests ==========

    /// Helper: run the full RHF + Hessian + IR pipeline for a molecule
    /// and return the IrSpectrumResult.
    fn run_ir_pipeline(atoms: &[Atom]) -> IrSpectrumResult {
        use crate::scf::build_density;
        let basis = BasisSet::build(atoms.to_vec(), "sto-3g").unwrap();
        let atoms_io: Vec<(u8, [f64; 3])> = atoms
            .iter()
            .map(|a| (a.atomic_number, a.position))
            .collect();

        let scf_config = ScfConfig {
            max_iterations: 200,
            use_diis: true,
            diis_start: 2,
            ..ScfConfig::tight()
        };

        let hess = rhf_hessian(&atoms_io, "sto-3g", &scf_config).unwrap();

        // Run the same re-diagonalization + density rebuild as rhf_hessian does
        // internally, so we get the fully self-consistent density.
        let cphf_data = hess.mo1_cphf.as_ref().expect("mo1_cphf");
        let density = build_density(&cphf_data.mo_coeff, cphf_data.n_occ);

        let freq_info = harmonic_analysis(&basis.atoms, &hess.hessian).expect("harmonic analysis");

        compute_ir_spectrum(
            &basis.atoms,
            &basis,
            &density,
            &hess,
            &freq_info,
            &[0.0, 0.0, 0.0],
        )
        .expect("IR spectrum")
    }

    #[test]
    fn test_ir_intensities_h2o_rhf_vs_pyscf() {
        // AC5: H2O IR intensities must match PySCF within 1 km/mol
        let golden = load_golden_ir("h2o");
        let atoms = build_atoms(&golden);
        let result = run_ir_pipeline(&atoms);

        assert_eq!(result.n_modes, 3);
        assert_eq!(result.intensities_km_per_mol.len(), 3);

        // Sort both (IQCP and golden) by frequency for robust comparison
        // (both should already be sorted but double-check)
        let iqcp_intensities = &result.intensities_km_per_mol;
        let pyscf_intensities = &golden.ir_intensities_km_mol;

        for k in 0..3 {
            let err = (iqcp_intensities[k] - pyscf_intensities[k]).abs();
            assert!(
                err < 1.0,
                "H2O IR intensity mode {k}: IQCP = {} km/mol, PySCF = {} km/mol, err = {}",
                iqcp_intensities[k],
                pyscf_intensities[k],
                err
            );
        }
    }

    #[test]
    fn test_ir_intensities_h2o_dipole_moment_also_matches() {
        // Cross-check: the equilibrium dipole in IrSpectrumResult should match
        // the PySCF reference to 1e-4 Debye (AC2 via integrated test)
        let golden = load_golden_ir("h2o");
        let atoms = build_atoms(&golden);
        let result = run_ir_pipeline(&atoms);
        for d in 0..3 {
            let err = (result.equilibrium_dipole_debye[d] - golden.dipole_debye[d]).abs();
            assert!(
                err < 1e-4,
                "H2O dipole[{d}] = {} Debye (IQCP) vs {} Debye (PySCF), err = {}",
                result.equilibrium_dipole_debye[d],
                golden.dipole_debye[d],
                err
            );
        }
    }

    #[test]
    fn test_ir_intensities_ch4_a1_inactive() {
        // AC6: CH4 A1 symmetric stretch (~3476 cm-1) has intensity ~0
        let golden = load_golden_ir("ch4");
        let atoms = build_atoms(&golden);
        let result = run_ir_pipeline(&atoms);

        assert_eq!(result.n_modes, 9);

        // CH4 has 9 modes:
        //   3 T2 bend (~1688 cm-1): IR active (~5.94 km/mol each)
        //   2 E bend (~1908 cm-1): IR inactive (0 by symmetry)
        //   1 A1 stretch (~3476 cm-1): IR inactive
        //   3 T2 stretch (~3733 cm-1): IR active (~0.67 km/mol each)
        //
        // Identify the A1 stretch by frequency: the singleton mode in the
        // range 3400-3550 cm-1. Use a PySCF golden-sorted index search.
        let pyscf_a1_freq = golden
            .freq_wavenumber
            .iter()
            .find(|&&f| (3400.0..3550.0).contains(&f))
            .copied()
            .expect("A1 symmetric stretch at ~3476 cm-1");
        let a1_idx = golden
            .freq_wavenumber
            .iter()
            .position(|&f| (f - pyscf_a1_freq).abs() < 1.0)
            .unwrap();

        let a1_intensity = result.intensities_km_per_mol[a1_idx];
        assert!(
            a1_intensity.abs() < 0.01,
            "CH4 A1 symmetric stretch at idx {} (freq {} cm-1) has intensity {} km/mol, expected < 0.01",
            a1_idx,
            pyscf_a1_freq,
            a1_intensity
        );

        // Also verify: max deviation from PySCF for all modes < 1 km/mol
        for k in 0..9 {
            let err = (result.intensities_km_per_mol[k] - golden.ir_intensities_km_mol[k]).abs();
            assert!(
                err < 1.0,
                "CH4 IR intensity mode {k}: IQCP = {}, PySCF = {}, err = {}",
                result.intensities_km_per_mol[k],
                golden.ir_intensities_km_mol[k],
                err
            );
        }
    }

    #[test]
    fn test_ir_intensities_co2_mutual_exclusion() {
        // CO2 σ_g symmetric stretch (~1571 cm-1) has intensity ~0
        //     σ_u asymmetric stretch (~2829 cm-1) has intensity > 100 km/mol
        // This tests the D∞h mutual exclusion preview for US-098
        let golden = load_golden_ir("co2");
        let atoms = build_atoms(&golden);
        let result = run_ir_pipeline(&atoms);

        assert_eq!(result.n_modes, 4);

        let sigma_g_idx = golden
            .freq_wavenumber
            .iter()
            .position(|&f| (1500.0..1650.0).contains(&f))
            .expect("σ_g symmetric stretch at ~1571 cm-1");
        let sigma_u_idx = golden
            .freq_wavenumber
            .iter()
            .position(|&f| (2700.0..2900.0).contains(&f))
            .expect("σ_u asymmetric stretch at ~2829 cm-1");

        let sigma_g_intensity = result.intensities_km_per_mol[sigma_g_idx];
        let sigma_u_intensity = result.intensities_km_per_mol[sigma_u_idx];

        assert!(
            sigma_g_intensity.abs() < 0.01,
            "CO2 σ_g symmetric stretch intensity = {} km/mol, expected < 0.01",
            sigma_g_intensity
        );
        assert!(
            sigma_u_intensity > 100.0,
            "CO2 σ_u asymmetric stretch intensity = {} km/mol, expected > 100",
            sigma_u_intensity
        );
    }

    #[test]
    fn test_ir_intensities_h2_zero() {
        // H2 homonuclear: single mode, zero intensity
        let golden = load_golden_ir("h2");
        let atoms = build_atoms(&golden);
        let result = run_ir_pipeline(&atoms);
        assert_eq!(result.n_modes, 1);
        assert!(
            result.intensities_km_per_mol[0] < 0.01,
            "H2 single-mode intensity = {} km/mol, expected < 0.01",
            result.intensities_km_per_mol[0]
        );
        // Also: equilibrium dipole is zero
        for d in 0..3 {
            assert_abs_diff_eq!(result.equilibrium_dipole_au[d], 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_dipole_derivatives_h2o_vs_pyscf() {
        // Component-wise comparison of dipole derivatives against PySCF
        let golden = load_golden_ir("h2o");
        let atoms = build_atoms(&golden);
        let result = run_ir_pipeline(&atoms);

        let dd = &result.dipole_derivs_cartesian;
        assert_eq!(dd.nrows(), 3);
        assert_eq!(dd.ncols(), 3 * 3);

        for d in 0..3 {
            for k in 0..9 {
                let iqcp_val = dd[(d, k)];
                let pyscf_val = golden.dipole_derivs_cartesian[d][k];
                let err = (iqcp_val - pyscf_val).abs();
                assert!(
                    err < 1e-5,
                    "dipole_derivs[{d}][{k}] = {} (IQCP) vs {} (PySCF), err = {}",
                    iqcp_val,
                    pyscf_val,
                    err
                );
            }
        }
    }

    #[test]
    fn test_dipole_derivatives_h2o_translational_invariance() {
        // Sum over atoms of ∂μ/∂R_A should be zero for a neutral molecule.
        let golden = load_golden_ir("h2o");
        let atoms = build_atoms(&golden);
        let result = run_ir_pipeline(&atoms);
        let dd = &result.dipole_derivs_cartesian;
        let n_atoms = atoms.len();
        for d in 0..3 {
            for e in 0..3 {
                let mut sum = 0.0;
                for a in 0..n_atoms {
                    sum += dd[(d, 3 * a + e)];
                }
                assert!(
                    sum.abs() < 1e-8,
                    "Translational sum d={d} e={e} = {}, expected 0",
                    sum
                );
            }
        }
    }
}
