//! Coupled-Perturbed Hartree-Fock (CPHF) solver
//!
//! Solves the CPHF equations to determine the first-order response of MO
//! coefficients to external perturbations (nuclear displacements or electric
//! fields). The CPHF solver produces the first-order density matrix response
//! needed for analytical Hessian, IR intensity, and Raman activity calculations.
//!
//! # Two solver paths
//!
//! - `cphf_solve_withs1()`: Nuclear perturbation path (field-dependent basis,
//!   nonzero S^(1)). The occupied-occupied block is fixed at -S^(1)/2.
//!   Returns both U^(1) and first-order orbital energies.
//!
//! - `cphf_solve_nos1()`: Electric field perturbation path (field-independent
//!   basis, S^(1) = 0). Only the virtual-occupied block is solved.
//!
//! Both paths use a custom Lanczos-like Krylov subspace solver matching
//! PySCF's `lib.krylov()` (linalg_helper.py lines 1222-1376).
//!
//! # Algorithm
//!
//! The CPHF equations (Helgaker et al. 2000, Eq. 10.6.15):
//! ```text
//! (e_a - e_i) U^(1)_{ai} + sum_{bj} A_{ai,bj} U^(1)_{bj}
//!     = -h^(1)_{ai} + e_i * s^(1)_{ai}
//! ```
//!
//! are reformulated as `(1 + A_tilde) x = b` where the orbital energy
//! denominator is absorbed into the operator, and solved by Krylov iteration.
//!
//! # References
//!
//! - PySCF `cphf.py` lines 29-148: solve_withs1, solve_nos1
//! - PySCF `linalg_helper.py` lines 1222-1376: krylov()
//! - PySCF `hessian/rhf.py` lines 286-360: solve_mo1, gen_vind
//! - Pople et al. (1979). Int. J. Quantum Chem. Symp. 13, 225.
//! - Helgaker, Jorgensen & Olsen (2000). Ch. 10: CPHF equations.

use nalgebra::DMatrix;

use super::eri_get;

// =============================================================================
// Data Types
// =============================================================================

/// Configuration for the CPHF solver.
///
/// # Fields
///
/// * `max_cycle` - Maximum Krylov iterations (default: 50)
/// * `tol` - Convergence tolerance; solver converges when `||r||^2 < tol^2`
///   (default: 1e-9, matching PySCF's `conv_tol_cpscf`)
/// * `level_shift` - Added to orbital energy denominator for near-degenerate
///   systems (default: 0.0)
///
/// # Reference
///
/// PySCF `cphf.py` line 31: max_cycle=50, tol=1e-9, level_shift=0
#[derive(Debug, Clone)]
pub struct CphfConfig {
    pub max_cycle: usize,
    pub tol: f64,
    pub level_shift: f64,
}

impl Default for CphfConfig {
    fn default() -> Self {
        Self {
            max_cycle: 50,
            tol: 1e-9,
            level_shift: 0.0,
        }
    }
}

/// Result of the CPHF solver.
///
/// Contains the first-order MO coefficients and orbital energies for
/// all perturbations.
///
/// # Layout
///
/// * `mo1[k]` - nmo x nocc DMatrix (for solve_withs1) or nvir x nocc
///   (for solve_nos1), where k indexes perturbation (e.g., atom*3+dir)
/// * `mo_e1[k]` - nocc x nocc DMatrix of first-order orbital energies
///   (only for solve_withs1; None for solve_nos1)
pub struct CphfResult {
    /// First-order MO coefficients U^(1)_{pq} in MO basis.
    /// Indexed as [perturbation_index] -> DMatrix.
    pub mo1: Vec<DMatrix<f64>>,
    /// First-order orbital energies (occupied block).
    /// Indexed as [perturbation_index] -> DMatrix (nocc x nocc).
    /// None for solve_nos1 (field-independent basis).
    pub mo_e1: Option<Vec<DMatrix<f64>>>,
    /// Number of Krylov iterations to convergence.
    pub iterations: usize,
    /// Whether the CPHF equations converged within max_cycle iterations.
    pub converged: bool,
}

// =============================================================================
// Krylov Solver
// =============================================================================

/// Linear dependence threshold for Krylov subspace.
///
/// When the squared norm of an orthogonalized Krylov vector falls below
/// this threshold, the vector is discarded to avoid numerical instability.
///
/// Reference: PySCF linalg_helper.py line 1223 (DSOLVE_LINDEP)
const KRYLOV_LINDEP: f64 = 1e-14;

/// Custom Lanczos-like Krylov subspace solver for `(1+A)x = b`.
///
/// Matches PySCF's `lib.krylov()` (linalg_helper.py lines 1222-1376).
///
/// # Algorithm
///
/// 1. Initialize: x1 = b; QR-orthogonalize (not normalize); track innerprod
/// 2. Iterate: axt = A(x1); orthogonalize against all previous xs; check convergence
/// 3. Project: build Gram matrix h; solve h*c = g; reconstruct x = sum c_i * xs_i
///
/// # Arguments
///
/// * `aop` - The operator A: given flat x (nroots * ndim_per_root), returns A*x
/// * `b` - Right-hand side: flat f64 array of shape (nroots * ndim_per_root)
/// * `nroots` - Number of simultaneous RHS vectors
/// * `tol` - Convergence tolerance (convergence when max_innerprod < tol^2)
/// * `max_cycle` - Maximum iterations
///
/// # Returns
///
/// `(solution, n_iterations, converged)` where solution is a flat Vec<f64>
/// of the same shape as b.
///
/// # Reference
///
/// PySCF linalg_helper.py lines 1266-1376
#[allow(clippy::needless_range_loop, clippy::collapsible_if)]
fn krylov_solve<F>(
    aop: &mut F,
    b: &[f64],
    nroots: usize,
    tol: f64,
    max_cycle: usize,
) -> (Vec<f64>, usize, bool)
where
    F: FnMut(&[f64]) -> Vec<f64>,
{
    let ndim = b.len() / nroots;
    assert_eq!(b.len(), nroots * ndim, "b length must be nroots * ndim");

    // Initialize x1 = b (each root is one row of shape ndim)
    // PySCF linalg_helper.py line 1276: x1 = b
    let mut x1: Vec<Vec<f64>> = (0..nroots)
        .map(|r| b[r * ndim..(r + 1) * ndim].to_vec())
        .collect();

    // For multiple roots: QR orthogonalize but DON'T normalize.
    // PySCF line 1281-1285: x1, rmat = _qr(x1, dot); x1 *= rmat.diagonal()[:,None]
    let mut innerprod: Vec<f64>;
    let max_innerprod: f64;

    if nroots > 1 {
        let (q_vecs, r_diag) = krylov_qr(&x1);
        // Scale back by r diagonal: vectors are orthogonal but not normalized
        x1 = q_vecs
            .iter()
            .zip(r_diag.iter())
            .map(|(q, &d)| q.iter().map(|&v| v * d).collect())
            .collect();
        innerprod = r_diag.iter().map(|&d| d * d).collect();
        max_innerprod = innerprod.iter().cloned().fold(0.0f64, f64::max);
    } else {
        let ip: f64 = x1[0].iter().map(|&v| v * v).sum();
        innerprod = vec![ip];
        max_innerprod = ip;
    }

    // PySCF line 1294: if max_innerprod < lindep or < tol**2, already converged
    if max_innerprod < KRYLOV_LINDEP || max_innerprod < tol * tol {
        return (vec![0.0; b.len()], 0, true);
    }

    // Limit max_cycle to ndim (PySCF line 1309)
    let max_cycle = max_cycle.min(ndim);

    // Storage for trial vectors and their A-products
    let mut xs: Vec<Vec<f64>> = Vec::new();
    let mut ax: Vec<Vec<f64>> = Vec::new();

    let mut converged = false;
    let mut n_iters = 0;

    for cycle in 0..max_cycle {
        n_iters = cycle + 1;

        // Apply operator: axt = A(x1) -- PySCF line 1311
        let x1_flat: Vec<f64> = x1.iter().flat_map(|v| v.iter().copied()).collect();
        let axt_flat = aop(&x1_flat);
        let mut axt: Vec<Vec<f64>> = (0..x1.len())
            .map(|r| axt_flat[r * ndim..(r + 1) * ndim].to_vec())
            .collect();

        // Store current x1 and axt -- PySCF lines 1314-1315
        xs.extend(x1.iter().cloned());
        ax.extend(axt.iter().cloned());

        // Orthogonalize: x1 = axt - sum_i (axt . xs[i]) / innerprod[i] * xs[i]
        // PySCF lines 1319-1323 (modified Gram-Schmidt against ALL previous vectors)
        let mut x1_new = axt.clone();
        for i in 0..xs.len() {
            // w[r] = dot(axt[r], xs[i]) / innerprod[i]
            for r in 0..x1_new.len() {
                let dot_val: f64 = axt[r].iter().zip(xs[i].iter()).map(|(&a, &x)| a * x).sum();
                let w = dot_val / innerprod[i];
                for j in 0..ndim {
                    x1_new[r][j] -= xs[i][j] * w;
                }
            }
        }
        // Clear axt reference
        axt.clear();

        // QR and innerprod for new vectors -- PySCF lines 1326-1331
        let innerprod1: Vec<f64>;
        if nroots > 1 || x1_new.len() > 1 {
            let (q_vecs, r_diag) = krylov_qr(&x1_new);
            x1_new = q_vecs
                .iter()
                .zip(r_diag.iter())
                .map(|(q, &d)| q.iter().map(|&v| v * d).collect())
                .collect();
            innerprod1 = r_diag.iter().map(|&d| d * d).collect();
        } else {
            let ip: f64 = x1_new[0].iter().map(|&v| v * v).sum();
            innerprod1 = vec![ip];
        }

        let max_innerprod = innerprod1.iter().cloned().fold(0.0f64, f64::max);

        // Convergence check -- PySCF line 1335
        if max_innerprod < KRYLOV_LINDEP || max_innerprod < tol * tol {
            converged = true;
            break;
        }

        // Filter linearly independent vectors -- PySCF lines 1338-1343
        let tol_sq = tol * tol;
        let mut filtered_x1 = Vec::new();
        for (idx, ip) in innerprod1.iter().enumerate() {
            if *ip > KRYLOV_LINDEP && *ip > tol_sq {
                if idx < x1_new.len() {
                    filtered_x1.push(x1_new[idx].clone());
                    innerprod.push(*ip);
                }
            }
        }
        x1 = filtered_x1;

        if x1.is_empty() {
            converged = true;
            break;
        }
    }

    // Projection step: build Gram matrix and solve -- PySCF lines 1348-1369
    let nd = xs.len();
    if nd == 0 {
        return (vec![0.0; b.len()], n_iters, converged);
    }

    // h[i,j] = dot(xs[i], ax[j]) -- PySCF line 1357-1358
    let mut h = DMatrix::zeros(nd, nd);
    for i in 0..nd {
        for j in 0..nd {
            let dot_val: f64 = xs[i].iter().zip(ax[j].iter()).map(|(&x, &a)| x * a).sum();
            h[(i, j)] = dot_val;
        }
    }

    // Add identity contribution: h[i,i] += innerprod[i] -- PySCF lines 1361-1362
    for i in 0..nd {
        h[(i, i)] += innerprod[i];
    }

    // g[i, r] = dot(b_r, xs[i]) -- PySCF lines 1364-1366
    let mut g = DMatrix::zeros(nd, nroots);
    for i in 0..nd {
        for r in 0..nroots {
            let b_r = &b[r * ndim..(r + 1) * ndim];
            let dot_val: f64 = b_r.iter().zip(xs[i].iter()).map(|(&b, &x)| b * x).sum();
            g[(i, r)] = dot_val;
        }
    }

    // Solve h * c = g -- PySCF line 1368
    let c = match h.clone().lu().solve(&g) {
        Some(c) => c,
        None => {
            // Fallback: try pseudoinverse
            match h.clone().pseudo_inverse(1e-14) {
                Ok(h_inv) => &h_inv * &g,
                Err(_) => {
                    return (vec![0.0; b.len()], n_iters, false);
                }
            }
        }
    };

    // Reconstruct x = sum_i c[i,r] * xs[i] -- PySCF line 1369
    let mut result = vec![0.0; b.len()];
    for r in 0..nroots {
        for i in 0..nd {
            let coeff = c[(i, r)];
            for j in 0..ndim {
                result[r * ndim + j] += coeff * xs[i][j];
            }
        }
    }

    (result, n_iters, converged)
}

/// QR decomposition for a set of vectors (orthogonal but not normalized).
///
/// Returns `(Q_normalized_vectors, R_diagonal_elements)` where Q has rows
/// that are unit-normalized, and R diagonal gives the original norms.
/// Vectors with norm below KRYLOV_LINDEP are discarded.
///
/// This matches PySCF's `_qr()` in linalg_helper.py lines 1412-1433,
/// but returns the diagonal of the R matrix (the norms) directly.
///
/// # Reference
///
/// PySCF linalg_helper.py lines 1412-1433
#[allow(clippy::needless_range_loop)]
fn krylov_qr(xs: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<f64>) {
    let nvec = xs.len();
    if nvec == 0 {
        return (Vec::new(), Vec::new());
    }

    let ndim = xs[0].len();
    let mut qs: Vec<Vec<f64>> = Vec::with_capacity(nvec);
    let mut norms: Vec<f64> = Vec::with_capacity(nvec);

    for i in 0..nvec {
        let mut xi = xs[i].clone();

        // Orthogonalize against existing basis
        for j in 0..qs.len() {
            let prod: f64 = qs[j].iter().zip(xi.iter()).map(|(&q, &x)| q * x).sum();
            for k in 0..ndim {
                xi[k] -= qs[j][k] * prod;
            }
        }

        // Check linear independence
        let innerprod: f64 = xi.iter().map(|&v| v * v).sum();
        let norm = innerprod.sqrt();

        if innerprod > KRYLOV_LINDEP {
            // Normalize
            for k in 0..ndim {
                xi[k] /= norm;
            }
            qs.push(xi);
            norms.push(norm);
        }
    }

    (qs, norms)
}

// =============================================================================
// CPHF Solver: solve_withs1 (nuclear perturbation path)
// =============================================================================

/// Nuclear perturbation path: solve CPHF with nonzero S^(1).
///
/// The occupied-occupied block is fixed at -S^(1)/2 and not iterated.
/// Returns both U^(1) and first-order orbital energies.
///
/// # Algorithm (PySCF cphf.py lines 86-148)
///
/// 1. Compute `hs = h1 - s1 * eps_i` (RHS before denominator)
/// 2. Initial guess: virtual block `= -(h1 - s1*eps_i) * e_ai`
///    occupied block `= -s1_oo / 2`
/// 3. Define `vind_vo` wrapper: apply e_ai to virtual block, zero occupied
/// 4. Call krylov_solve on the batched system
/// 5. Post-process: restore occupied block, recompute virtual exactly,
///    compute mo_e1
///
/// # Arguments
///
/// * `vind` - Response function: given flat mo1 `(npert * nmo * nocc)`,
///   returns induced potential in MO basis (same shape)
/// * `mo_energy` - Orbital energies (nmo-length slice)
/// * `n_occ` - Number of occupied orbitals
/// * `h1_mo` - First-order Hamiltonian in MO basis (one nmo x nocc per perturbation)
/// * `s1_mo` - First-order overlap in MO basis (one nmo x nocc per perturbation)
/// * `config` - CPHF configuration
///
/// # Returns
///
/// `(mo1, mo_e1, iterations, converged)` where:
/// - `mo1[k]` is nmo x nocc DMatrix for perturbation k
/// - `mo_e1[k]` is nocc x nocc DMatrix of first-order orbital energies
///
/// # Reference
///
/// PySCF cphf.py lines 86-148
#[allow(clippy::needless_range_loop)]
pub fn cphf_solve_withs1<F>(
    vind: F,
    mo_energy: &[f64],
    n_occ: usize,
    h1_mo: &[DMatrix<f64>],
    s1_mo: &[DMatrix<f64>],
    config: &CphfConfig,
) -> (Vec<DMatrix<f64>>, Vec<DMatrix<f64>>, usize, bool)
where
    F: FnMut(&[f64]) -> Vec<f64>,
{
    let npert = h1_mo.len();
    assert_eq!(
        s1_mo.len(),
        npert,
        "h1_mo and s1_mo must have same number of perturbations"
    );
    assert!(!h1_mo.is_empty(), "Must have at least one perturbation");

    let nmo = h1_mo[0].nrows();
    let nocc = h1_mo[0].ncols();
    assert_eq!(nocc, n_occ);
    let nvir = nmo - nocc;

    // Orbital energies
    let e_i: Vec<f64> = mo_energy[..nocc].to_vec();
    let e_a: Vec<f64> = mo_energy[nocc..nmo].to_vec();

    // e_ai = 1 / (e_a + level_shift - e_i)
    // PySCF cphf.py line 111
    let mut e_ai = DMatrix::zeros(nvir, nocc);
    for a in 0..nvir {
        for i in 0..nocc {
            e_ai[(a, i)] = 1.0 / (e_a[a] + config.level_shift - e_i[i]);
        }
    }

    // Step 1: hs = h1 - s1 * e_i  (PySCF line 117)
    let mut hs: Vec<DMatrix<f64>> = Vec::with_capacity(npert);
    for k in 0..npert {
        let mut hs_k = h1_mo[k].clone();
        for p in 0..nmo {
            for i in 0..nocc {
                hs_k[(p, i)] -= s1_mo[k][(p, i)] * e_i[i];
            }
        }
        hs.push(hs_k);
    }

    // Step 2: Initial guess (PySCF lines 118-120)
    let mut mo1base: Vec<DMatrix<f64>> = Vec::with_capacity(npert);
    for k in 0..npert {
        let mut base = hs[k].clone();
        // Virtual block: *= -e_ai
        for a in 0..nvir {
            for i in 0..nocc {
                base[(nocc + a, i)] *= -e_ai[(a, i)];
            }
        }
        // Occupied block: = -s1_oo / 2
        for oi in 0..nocc {
            for oj in 0..nocc {
                base[(oi, oj)] = -s1_mo[k][(oi, oj)] * 0.5;
            }
        }
        mo1base.push(base);
    }

    // Flatten mo1base for krylov solver
    let ndim_per_root = nmo * nocc;
    let nroots = npert;
    let b_flat: Vec<f64> = mo1base
        .iter()
        .flat_map(|m| (0..nmo).flat_map(move |p| (0..nocc).map(move |i| m[(p, i)])))
        .collect();

    // Step 3: vind_vo wrapper (PySCF lines 122-129)
    // The wrapper applies e_ai to virtual block and zeros occupied block.
    // We use a RefCell to share vind between the closure and post-processing.
    let level_shift = config.level_shift;
    let e_ai_clone = e_ai.clone();
    let vind_cell = std::cell::RefCell::new(vind);
    let mut vind_vo = |x: &[f64]| -> Vec<f64> {
        // Call the raw vind
        let v_flat = vind_cell.borrow_mut()(x);
        let n = x.len() / ndim_per_root;
        let mut result = vec![0.0; v_flat.len()];

        for r in 0..n {
            let offset = r * ndim_per_root;
            // If level_shift != 0, subtract mo1 * level_shift
            // PySCF cphf.py lines 126-127
            for p in 0..nmo {
                for i in 0..nocc {
                    let idx = offset + p * nocc + i;
                    let mut v = v_flat[idx];
                    if level_shift.abs() > 1e-15 {
                        v -= x[idx] * level_shift;
                    }
                    // Virtual block: multiply by e_ai (PySCF line 128)
                    if p >= nocc {
                        result[idx] = v * e_ai_clone[(p - nocc, i)];
                    }
                    // Occupied block: zero (PySCF line 129)
                    // result[idx] = 0.0 is already initialized
                }
            }
        }
        result
    };

    // Step 4: Krylov solve (PySCF lines 130-131)
    let (mo1_flat, iterations, converged) =
        krylov_solve(&mut vind_vo, &b_flat, nroots, config.tol, config.max_cycle);
    // Release the vind_vo closure to drop the borrow on vind_cell
    let _ = vind_vo;

    // Step 5: Post-processing (PySCF lines 132-148)
    // Reshape to matrices
    let mut mo1: Vec<DMatrix<f64>> = Vec::with_capacity(npert);
    for k in 0..npert {
        let mut m = DMatrix::zeros(nmo, nocc);
        for p in 0..nmo {
            for i in 0..nocc {
                m[(p, i)] = mo1_flat[k * ndim_per_root + p * nocc + i];
            }
        }
        mo1.push(m);
    }

    // Restore occupied block = mo1base occupied block = -s1/2
    // PySCF line 133: mo1[:,occidx] = mo1base[:,occidx]
    for k in 0..npert {
        for oi in 0..nocc {
            for oj in 0..nocc {
                mo1[k][(oi, oj)] = mo1base[k][(oi, oj)];
            }
        }
    }

    // Compute final response with converged mo1
    // PySCF line 136: hs += fvind(mo1)
    let mo1_for_vind: Vec<f64> = mo1
        .iter()
        .flat_map(|m| (0..nmo).flat_map(move |p| (0..nocc).map(move |i| m[(p, i)])))
        .collect();
    let v_final = vind_cell.borrow_mut()(&mo1_for_vind);

    for k in 0..npert {
        let offset = k * ndim_per_root;
        for p in 0..nmo {
            for i in 0..nocc {
                hs[k][(p, i)] += v_final[offset + p * nocc + i];
            }
        }
    }

    // Recompute virtual block exactly from final response
    // PySCF line 137: mo1[:,viridx] = hs[:,viridx] / (e_i - e_a[:,None])
    for k in 0..npert {
        for a in 0..nvir {
            for i in 0..nocc {
                mo1[k][(nocc + a, i)] = hs[k][(nocc + a, i)] / (e_i[i] - e_a[a]);
            }
        }
    }

    // Compute first-order orbital energies
    // PySCF lines 141-142:
    // mo_e1 = hs[:,occidx,:]
    // mo_e1 += mo1[:,occidx] * (e_i[:,None] - e_i)
    let mut mo_e1: Vec<DMatrix<f64>> = Vec::with_capacity(npert);
    for k in 0..npert {
        let mut e1 = DMatrix::zeros(nocc, nocc);
        for oi in 0..nocc {
            for oj in 0..nocc {
                e1[(oi, oj)] = hs[k][(oi, oj)] + mo1[k][(oi, oj)] * (e_i[oi] - e_i[oj]);
            }
        }
        mo_e1.push(e1);
    }

    (mo1, mo_e1, iterations, converged)
}

// =============================================================================
// CPHF Solver: solve_nos1 (electric field perturbation path)
// =============================================================================

/// Electric field perturbation path: solve CPHF without S^(1).
///
/// Only the virtual-occupied block is solved. No occupied-occupied block
/// or first-order orbital energies are computed.
///
/// # Algorithm (PySCF cphf.py lines 53-83)
///
/// 1. Compute `e_ai = 1 / (e_a + level_shift - e_i)`
/// 2. Initial guess: `mo1base = -h1 * e_ai`
/// 3. Define `vind_vo` wrapper: apply e_ai, subtract level_shift
/// 4. Call krylov_solve
/// 5. Return mo1 in (nvir x nocc) shape
///
/// # Arguments
///
/// * `vind` - Response function (operates on virtual-occupied block only)
/// * `mo_energy` - Orbital energies (nmo-length slice)
/// * `n_occ` - Number of occupied orbitals
/// * `h1_mo` - First-order Hamiltonian in MO basis (nvir x nocc per perturbation)
/// * `config` - CPHF configuration
///
/// # Returns
///
/// `(mo1, iterations, converged)` where mo1[k] is nvir x nocc DMatrix
///
/// # Reference
///
/// PySCF cphf.py lines 53-83
#[allow(clippy::needless_range_loop)]
pub fn cphf_solve_nos1<F>(
    mut vind: F,
    mo_energy: &[f64],
    n_occ: usize,
    h1_mo: &[DMatrix<f64>],
    config: &CphfConfig,
) -> (Vec<DMatrix<f64>>, usize, bool)
where
    F: FnMut(&[f64]) -> Vec<f64>,
{
    let npert = h1_mo.len();
    assert!(!h1_mo.is_empty(), "Must have at least one perturbation");

    let nvir = h1_mo[0].nrows();
    let nocc = h1_mo[0].ncols();
    assert_eq!(nocc, n_occ);
    let nmo = nvir + nocc;

    // Orbital energies
    let e_i: Vec<f64> = mo_energy[..nocc].to_vec();
    let e_a: Vec<f64> = mo_energy[nocc..nmo].to_vec();

    // e_ai = 1 / (e_a + level_shift - e_i)
    // PySCF cphf.py line 69
    let mut e_ai = DMatrix::zeros(nvir, nocc);
    for a in 0..nvir {
        for i in 0..nocc {
            e_ai[(a, i)] = 1.0 / (e_a[a] + config.level_shift - e_i[i]);
        }
    }

    // Initial guess: mo1base = h1 * -e_ai (PySCF line 70)
    let mut mo1base: Vec<DMatrix<f64>> = Vec::with_capacity(npert);
    for k in 0..npert {
        let mut base = DMatrix::zeros(nvir, nocc);
        for a in 0..nvir {
            for i in 0..nocc {
                base[(a, i)] = h1_mo[k][(a, i)] * (-e_ai[(a, i)]);
            }
        }
        mo1base.push(base);
    }

    // Flatten for krylov solver
    let ndim_per_root = nvir * nocc;
    let nroots = npert;
    let b_flat: Vec<f64> = mo1base
        .iter()
        .flat_map(|m| (0..nvir).flat_map(move |a| (0..nocc).map(move |i| m[(a, i)])))
        .collect();

    // vind_vo wrapper (PySCF lines 73-79)
    let level_shift = config.level_shift;
    let e_ai_clone = e_ai.clone();
    let mut vind_vo = move |x: &[f64]| -> Vec<f64> {
        let v_flat = vind(x);
        let n = x.len() / ndim_per_root;
        let mut result = vec![0.0; v_flat.len()];

        for r in 0..n {
            let offset = r * ndim_per_root;
            for a in 0..nvir {
                for i in 0..nocc {
                    let idx = offset + a * nocc + i;
                    let mut v = v_flat[idx];
                    if level_shift.abs() > 1e-15 {
                        v -= x[idx] * level_shift;
                    }
                    result[idx] = v * e_ai_clone[(a, i)];
                }
            }
        }
        result
    };

    // Krylov solve
    let (mo1_flat, iterations, converged) =
        krylov_solve(&mut vind_vo, &b_flat, nroots, config.tol, config.max_cycle);

    // Reshape to matrices
    let mut mo1: Vec<DMatrix<f64>> = Vec::with_capacity(npert);
    for k in 0..npert {
        let mut m = DMatrix::zeros(nvir, nocc);
        for a in 0..nvir {
            for i in 0..nocc {
                m[(a, i)] = mo1_flat[k * ndim_per_root + a * nocc + i];
            }
        }
        mo1.push(m);
    }

    (mo1, iterations, converged)
}

// =============================================================================
// CPHF Solver: dispatcher
// =============================================================================

/// Solve the CPHF equations (dispatcher).
///
/// Routes to `cphf_solve_withs1` if `s1_mo` is provided (nuclear perturbation),
/// or `cphf_solve_nos1` if `s1_mo` is None (electric field perturbation).
///
/// The perturbation matrices are provided per-atom as `[dX, dY, dZ]` arrays.
/// Internally, these are flattened to a batch for the Krylov solver, and
/// the tolerance is scaled by the number of atoms (PySCF convention).
///
/// # Arguments
///
/// * `vind` - Response function
/// * `mo_energy` - Orbital energies
/// * `n_occ` - Number of occupied orbitals
/// * `h1_mo` - First-order Hamiltonian `[atom][dir]` arrays of nmo x nocc
/// * `s1_mo` - First-order overlap `[atom][dir]` arrays, or None
/// * `config` - CPHF configuration
///
/// # Reference
///
/// PySCF cphf.py lines 29-50, hessian/rhf.py lines 316-332 (batching)
pub fn cphf_solve<F>(
    vind: F,
    mo_energy: &[f64],
    n_occ: usize,
    h1_mo: &[[DMatrix<f64>; 3]],
    s1_mo: Option<&[[DMatrix<f64>; 3]]>,
    config: &CphfConfig,
) -> CphfResult
where
    F: FnMut(&[f64]) -> Vec<f64>,
{
    let n_atoms = h1_mo.len();

    // Flatten [atom][3] -> batch of 3*n_atoms perturbations
    let h1_batch: Vec<DMatrix<f64>> = h1_mo.iter().flat_map(|dirs| dirs.iter().cloned()).collect();

    // Scale tolerance by number of atoms (PySCF rhf.py line 330)
    let mut scaled_config = config.clone();
    scaled_config.tol = config.tol * n_atoms as f64;

    if let Some(s1) = s1_mo {
        assert_eq!(s1.len(), n_atoms);
        let s1_batch: Vec<DMatrix<f64>> = s1.iter().flat_map(|dirs| dirs.iter().cloned()).collect();

        let (mo1_flat, mo_e1_flat, iterations, converged) =
            cphf_solve_withs1(vind, mo_energy, n_occ, &h1_batch, &s1_batch, &scaled_config);

        // Reshape back to [atom][3]
        let mo1: Vec<DMatrix<f64>> = mo1_flat;
        let mo_e1: Vec<DMatrix<f64>> = mo_e1_flat;

        CphfResult {
            mo1,
            mo_e1: Some(mo_e1),
            iterations,
            converged,
        }
    } else {
        let (mo1_flat, iterations, converged) =
            cphf_solve_nos1(vind, mo_energy, n_occ, &h1_batch, &scaled_config);

        CphfResult {
            mo1: mo1_flat,
            mo_e1: None,
            iterations,
            converged,
        }
    }
}

// =============================================================================
// Response Function: RHF vind
// =============================================================================

/// Generate the RHF response function (vind).
///
/// Returns a closure that, given trial U^(1) as a flat vector, computes
/// the induced potential `J[D^(1)] - c_hf*K[D^(1)]` in MO basis.
///
/// For each trial U^(1) (nmo x nocc or nvir x nocc):
/// 1. Build trial density in AO basis: `D^(1) = C * mo1 * 2 * C_occ^T + transpose`
/// 2. Compute `J[D^(1)] - c_hf * K[D^(1)]`
/// 3. Transform back to MO basis: `V_mo = C^T * V_AO * C_occ`
///
/// # Arguments
///
/// * `mo_coeff` - Full MO coefficient matrix C (nbf x nmo)
/// * `n_occ` - Number of occupied orbitals
/// * `eri` - Two-electron integrals in compressed storage
/// * `nbf` - Number of basis functions
/// * `hf_exchange_fraction` - Fraction of HF exchange (1.0 for RHF, 0.2 for B3LYP)
///
/// # Reference
///
/// PySCF hessian/rhf.py lines 343-360 (gen_vind)
pub fn gen_vind_rhf(
    mo_coeff: &DMatrix<f64>,
    n_occ: usize,
    eri: &[f64],
    nbf: usize,
    hf_exchange_fraction: f64,
) -> impl FnMut(&[f64]) -> Vec<f64> {
    let c = mo_coeff.clone();
    let nmo = c.ncols();
    let c_occ = c.columns(0, n_occ).clone_owned();
    let eri = eri.to_vec();
    let ct = c.transpose();

    move |mo1_flat: &[f64]| -> Vec<f64> {
        let ndim_per_root = nmo * n_occ;
        let nset = mo1_flat.len() / ndim_per_root;
        let mut result = vec![0.0; mo1_flat.len()];

        for s in 0..nset {
            let offset = s * ndim_per_root;

            // Reshape mo1 to nmo x nocc DMatrix
            let mut mo1 = DMatrix::zeros(nmo, n_occ);
            for p in 0..nmo {
                for i in 0..n_occ {
                    mo1[(p, i)] = mo1_flat[offset + p * n_occ + i];
                }
            }

            // Build trial density: dm = C * mo1 * 2 * C_occ^T
            // PySCF line 353: dm = reduce(dot, (mo_coeff, x*2, mocc.T))
            let dm_half = &c * &mo1 * 2.0 * c_occ.transpose();
            // Symmetrize: dm1 = dm + dm.T (PySCF line 354)
            let dm1 = &dm_half + &dm_half.transpose();

            // Compute J[D^(1)] and K[D^(1)]
            // Use the same algorithm as build_fock but with D^(1) instead of D
            let mut j_matrix = DMatrix::zeros(nbf, nbf);
            let mut k_matrix = DMatrix::zeros(nbf, nbf);

            for mu in 0..nbf {
                for nu in 0..=mu {
                    let mut j_mn = 0.0;
                    let mut k_mn = 0.0;

                    for lambda in 0..nbf {
                        // Diagonal: sigma = lambda
                        {
                            let d_ll = dm1[(lambda, lambda)];
                            let j_int = eri_get(&eri, mu, nu, lambda, lambda);
                            let k_int = eri_get(&eri, mu, lambda, nu, lambda);
                            j_mn += d_ll * j_int;
                            k_mn += d_ll * k_int;
                        }
                        // Off-diagonal: sigma < lambda
                        for sigma in 0..lambda {
                            let d_ls = dm1[(lambda, sigma)];
                            let j_int = eri_get(&eri, mu, nu, lambda, sigma);
                            j_mn += 2.0 * d_ls * j_int;
                            let k_int_1 = eri_get(&eri, mu, lambda, nu, sigma);
                            let k_int_2 = eri_get(&eri, mu, sigma, nu, lambda);
                            k_mn += d_ls * (k_int_1 + k_int_2);
                        }
                    }

                    j_matrix[(mu, nu)] = j_mn;
                    k_matrix[(mu, nu)] = k_mn;
                    if mu != nu {
                        j_matrix[(nu, mu)] = j_mn;
                        k_matrix[(nu, mu)] = k_mn;
                    }
                }
            }

            // V = J - c_hf * 0.5 * K
            let v_ao = &j_matrix - &k_matrix * (hf_exchange_fraction * 0.5);

            // Transform to MO basis: v_mo = C^T * V * C_occ
            // PySCF line 358: v1vo[i] = reduce(dot, (mo_coeff.T, x, mocc))
            let v_mo = &ct * &v_ao * &c_occ;

            // Store result
            for p in 0..nmo {
                for i in 0..n_occ {
                    result[offset + p * n_occ + i] = v_mo[(p, i)];
                }
            }
        }

        result
    }
}

/// Generate the RHF response function that operates on virtual-occupied
/// block only (for solve_nos1).
///
/// Same as `gen_vind_rhf` but the mo1 input/output is nvir x nocc.
pub fn gen_vind_rhf_vo(
    mo_coeff: &DMatrix<f64>,
    n_occ: usize,
    eri: &[f64],
    nbf: usize,
    hf_exchange_fraction: f64,
) -> impl FnMut(&[f64]) -> Vec<f64> {
    let c = mo_coeff.clone();
    let nmo = c.ncols();
    let nvir = nmo - n_occ;
    let c_occ = c.columns(0, n_occ).clone_owned();
    let eri = eri.to_vec();
    let ct = c.transpose();

    move |mo1_flat: &[f64]| -> Vec<f64> {
        let ndim_per_root = nvir * n_occ;
        let nset = mo1_flat.len() / ndim_per_root;
        let mut result = vec![0.0; mo1_flat.len()];

        for s in 0..nset {
            let offset = s * ndim_per_root;

            // Reshape mo1 to nmo x nocc (embed virtual-occupied block)
            let mut mo1_full = DMatrix::zeros(nmo, n_occ);
            for a in 0..nvir {
                for i in 0..n_occ {
                    mo1_full[(n_occ + a, i)] = mo1_flat[offset + a * n_occ + i];
                }
            }

            // Build trial density
            let dm_half = &c * &mo1_full * 2.0 * c_occ.transpose();
            let dm1 = &dm_half + &dm_half.transpose();

            // Compute J[D^(1)] and K[D^(1)]
            let mut j_matrix = DMatrix::zeros(nbf, nbf);
            let mut k_matrix = DMatrix::zeros(nbf, nbf);

            for mu in 0..nbf {
                for nu in 0..=mu {
                    let mut j_mn = 0.0;
                    let mut k_mn = 0.0;

                    for lambda in 0..nbf {
                        {
                            let d_ll = dm1[(lambda, lambda)];
                            let j_int = eri_get(&eri, mu, nu, lambda, lambda);
                            let k_int = eri_get(&eri, mu, lambda, nu, lambda);
                            j_mn += d_ll * j_int;
                            k_mn += d_ll * k_int;
                        }
                        for sigma in 0..lambda {
                            let d_ls = dm1[(lambda, sigma)];
                            let j_int = eri_get(&eri, mu, nu, lambda, sigma);
                            j_mn += 2.0 * d_ls * j_int;
                            let k_int_1 = eri_get(&eri, mu, lambda, nu, sigma);
                            let k_int_2 = eri_get(&eri, mu, sigma, nu, lambda);
                            k_mn += d_ls * (k_int_1 + k_int_2);
                        }
                    }

                    j_matrix[(mu, nu)] = j_mn;
                    k_matrix[(mu, nu)] = k_mn;
                    if mu != nu {
                        j_matrix[(nu, mu)] = j_mn;
                        k_matrix[(nu, mu)] = k_mn;
                    }
                }
            }

            // V = J - c_hf * 0.5 * K
            let v_ao = &j_matrix - &k_matrix * (hf_exchange_fraction * 0.5);

            // Transform to MO basis
            let v_mo = &ct * &v_ao * &c_occ;

            // Extract virtual-occupied block
            for a in 0..nvir {
                for i in 0..n_occ {
                    result[offset + a * n_occ + i] = v_mo[(n_occ + a, i)];
                }
            }
        }

        result
    }
}

// =============================================================================
// Response Function: DFT vind (AC6: DFT XC response framework)
// =============================================================================

/// Type alias for XC response callback.
///
/// Given a first-order density matrix D^(1) in AO basis (nbf x nbf),
/// returns the XC response contribution V_xc^(1) in AO basis (nbf x nbf):
///
/// ```text
/// V_xc^(1)_{μν} = ∫ f_xc(r) * ρ^(1)(r) * χ_μ(r) χ_ν(r) dr   (LDA)
///               + GGA gradient terms involving f''_{ρσ}, f''_{σσ}  (GGA)
/// ```
///
/// The 0.5 scale factor for GGA sigma terms is applied inside the callback
/// (PySCF convention: vsigma terms are multiplied by 0.5).
///
/// This callback is provided by the XC functional module (US-094).
/// For pure HF, pass `None`.
///
/// # Reference
///
/// PySCF `_response_functions.py` lines 66-98: DFT vind adds `nr_rks_fxc` to J-K
pub type XcResponseFn = dyn Fn(&DMatrix<f64>) -> DMatrix<f64>;

/// Generate a DFT response function (vind) with optional XC kernel.
///
/// Extends `gen_vind_rhf` to support DFT methods where:
/// - The HF exchange is scaled by `hf_exchange_fraction` (e.g., 0.2 for B3LYP)
/// - An XC response callback adds the exchange-correlation kernel contribution
///
/// The total response for each trial density D^(1) is:
/// ```text
/// V^(1) = J[D^(1)] - (c_hf / 2) * K[D^(1)] + V_xc[D^(1)]
/// ```
///
/// For pure HF: `hf_exchange_fraction = 1.0`, `xc_response = None`
/// For pure DFT (e.g., SVWN): `hf_exchange_fraction = 0.0`, `xc_response = Some(fxc_callback)`
/// For hybrid DFT (e.g., B3LYP): `hf_exchange_fraction = 0.2`, `xc_response = Some(fxc_callback)`
///
/// # Arguments
///
/// * `mo_coeff` - Full MO coefficient matrix C (nbf x nmo)
/// * `n_occ` - Number of occupied orbitals
/// * `eri` - Two-electron integrals in compressed storage
/// * `nbf` - Number of basis functions
/// * `hf_exchange_fraction` - Fraction of HF exchange (1.0 for RHF, 0.2 for B3LYP, 0.0 for pure DFT)
/// * `xc_response` - Optional XC response callback (None for pure HF)
///
/// # Reference
///
/// PySCF `_response_functions.py` lines 46-98: DFT gen_response
/// PySCF `hessian/rhf.py` lines 343-360: gen_vind wraps gen_response
pub fn gen_vind_dft<'a>(
    mo_coeff: &DMatrix<f64>,
    n_occ: usize,
    eri: &[f64],
    nbf: usize,
    hf_exchange_fraction: f64,
    xc_response: Option<&'a XcResponseFn>,
) -> impl FnMut(&[f64]) -> Vec<f64> + 'a {
    let c = mo_coeff.clone();
    let nmo = c.ncols();
    let c_occ = c.columns(0, n_occ).clone_owned();
    let eri = eri.to_vec();
    let ct = c.transpose();

    move |mo1_flat: &[f64]| -> Vec<f64> {
        let ndim_per_root = nmo * n_occ;
        let nset = mo1_flat.len() / ndim_per_root;
        let mut result = vec![0.0; mo1_flat.len()];

        for s in 0..nset {
            let offset = s * ndim_per_root;

            // Reshape mo1 to nmo x nocc DMatrix
            let mut mo1 = DMatrix::zeros(nmo, n_occ);
            for p in 0..nmo {
                for i in 0..n_occ {
                    mo1[(p, i)] = mo1_flat[offset + p * n_occ + i];
                }
            }

            // Build trial density: dm = C * mo1 * 2 * C_occ^T
            // PySCF line 353: dm = reduce(dot, (mo_coeff, x*2, mocc.T))
            let dm_half = &c * &mo1 * 2.0 * c_occ.transpose();
            // Symmetrize: dm1 = dm + dm.T (PySCF line 354)
            let dm1 = &dm_half + &dm_half.transpose();

            // Start with XC response if provided
            // PySCF _response_functions.py line 71-72: v1 = ni.nr_rks_fxc(...)
            let mut v_ao = if let Some(xc_fn) = &xc_response {
                xc_fn(&dm1)
            } else {
                DMatrix::zeros(nbf, nbf)
            };

            // Add J and K contributions (only if hf_exchange_fraction > 0 or pure HF)
            // PySCF _response_functions.py lines 76-94
            let mut j_matrix = DMatrix::zeros(nbf, nbf);
            let mut k_matrix = DMatrix::zeros(nbf, nbf);

            for mu in 0..nbf {
                for nu in 0..=mu {
                    let mut j_mn = 0.0;
                    let mut k_mn = 0.0;

                    for lambda in 0..nbf {
                        {
                            let d_ll = dm1[(lambda, lambda)];
                            let j_int = eri_get(&eri, mu, nu, lambda, lambda);
                            let k_int = eri_get(&eri, mu, lambda, nu, lambda);
                            j_mn += d_ll * j_int;
                            k_mn += d_ll * k_int;
                        }
                        for sigma in 0..lambda {
                            let d_ls = dm1[(lambda, sigma)];
                            let j_int = eri_get(&eri, mu, nu, lambda, sigma);
                            j_mn += 2.0 * d_ls * j_int;
                            let k_int_1 = eri_get(&eri, mu, lambda, nu, sigma);
                            let k_int_2 = eri_get(&eri, mu, sigma, nu, lambda);
                            k_mn += d_ls * (k_int_1 + k_int_2);
                        }
                    }

                    j_matrix[(mu, nu)] = j_mn;
                    k_matrix[(mu, nu)] = k_mn;
                    if mu != nu {
                        j_matrix[(nu, mu)] = j_mn;
                        k_matrix[(nu, mu)] = k_mn;
                    }
                }
            }

            // V += J - c_hf * 0.5 * K
            // PySCF _response_functions.py line 93: v1 += vj - .5 * vk
            v_ao += &j_matrix;
            v_ao -= &k_matrix * (hf_exchange_fraction * 0.5);

            // Transform to MO basis: v_mo = C^T * V * C_occ
            let v_mo = &ct * &v_ao * &c_occ;

            // Store result
            for p in 0..nmo {
                for i in 0..n_occ {
                    result[offset + p * n_occ + i] = v_mo[(p, i)];
                }
            }
        }

        result
    }
}

/// Generate a DFT response function that operates on virtual-occupied
/// block only (for solve_nos1).
///
/// Same as `gen_vind_dft` but the mo1 input/output is nvir x nocc.
pub fn gen_vind_dft_vo<'a>(
    mo_coeff: &DMatrix<f64>,
    n_occ: usize,
    eri: &[f64],
    nbf: usize,
    hf_exchange_fraction: f64,
    xc_response: Option<&'a XcResponseFn>,
) -> impl FnMut(&[f64]) -> Vec<f64> + 'a {
    let c = mo_coeff.clone();
    let nmo = c.ncols();
    let nvir = nmo - n_occ;
    let c_occ = c.columns(0, n_occ).clone_owned();
    let eri = eri.to_vec();
    let ct = c.transpose();

    move |mo1_flat: &[f64]| -> Vec<f64> {
        let ndim_per_root = nvir * n_occ;
        let nset = mo1_flat.len() / ndim_per_root;
        let mut result = vec![0.0; mo1_flat.len()];

        for s in 0..nset {
            let offset = s * ndim_per_root;

            // Reshape mo1 to nmo x nocc (embed virtual-occupied block)
            let mut mo1_full = DMatrix::zeros(nmo, n_occ);
            for a in 0..nvir {
                for i in 0..n_occ {
                    mo1_full[(n_occ + a, i)] = mo1_flat[offset + a * n_occ + i];
                }
            }

            // Build trial density
            let dm_half = &c * &mo1_full * 2.0 * c_occ.transpose();
            let dm1 = &dm_half + &dm_half.transpose();

            // XC response
            let mut v_ao = if let Some(xc_fn) = &xc_response {
                xc_fn(&dm1)
            } else {
                DMatrix::zeros(nbf, nbf)
            };

            // Compute J[D^(1)] and K[D^(1)]
            let mut j_matrix = DMatrix::zeros(nbf, nbf);
            let mut k_matrix = DMatrix::zeros(nbf, nbf);

            for mu in 0..nbf {
                for nu in 0..=mu {
                    let mut j_mn = 0.0;
                    let mut k_mn = 0.0;

                    for lambda in 0..nbf {
                        {
                            let d_ll = dm1[(lambda, lambda)];
                            let j_int = eri_get(&eri, mu, nu, lambda, lambda);
                            let k_int = eri_get(&eri, mu, lambda, nu, lambda);
                            j_mn += d_ll * j_int;
                            k_mn += d_ll * k_int;
                        }
                        for sigma in 0..lambda {
                            let d_ls = dm1[(lambda, sigma)];
                            let j_int = eri_get(&eri, mu, nu, lambda, sigma);
                            j_mn += 2.0 * d_ls * j_int;
                            let k_int_1 = eri_get(&eri, mu, lambda, nu, sigma);
                            let k_int_2 = eri_get(&eri, mu, sigma, nu, lambda);
                            k_mn += d_ls * (k_int_1 + k_int_2);
                        }
                    }

                    j_matrix[(mu, nu)] = j_mn;
                    k_matrix[(mu, nu)] = k_mn;
                    if mu != nu {
                        j_matrix[(nu, mu)] = j_mn;
                        k_matrix[(nu, mu)] = k_mn;
                    }
                }
            }

            // V += J - c_hf * 0.5 * K
            v_ao += &j_matrix;
            v_ao -= &k_matrix * (hf_exchange_fraction * 0.5);

            // Transform to MO basis
            let v_mo = &ct * &v_ao * &c_occ;

            // Extract virtual-occupied block
            for a in 0..nvir {
                for i in 0..n_occ {
                    result[offset + a * n_occ + i] = v_mo[(n_occ + a, i)];
                }
            }
        }

        result
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::{Atom, BasisSet};
    use crate::integrals;
    use crate::scf::hessian::{ao_to_mo, make_h1, make_s1};
    use crate::scf::{rhf_scf, ConvergenceProfile, PresetSystem, ScfConfig, ScfOutput};

    // ========================================================================
    // Test helpers: basis set construction (mirrors hessian.rs patterns)
    // ========================================================================

    /// Create H₂ STO-3G basis set at a given bond length (in bohr).
    fn h2_sto3g_basis(bond_length: f64) -> BasisSet {
        let atoms = vec![
            Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
            Atom::new(1, [0.0, 0.0, bond_length]).unwrap(),
        ];
        BasisSet::build(atoms, "sto-3g").expect("Failed to build H2 STO-3G basis")
    }

    /// Run SCF on H₂ and return output + basis.
    fn h2_scf(bond_length: f64) -> (ScfOutput, BasisSet) {
        let basis = h2_sto3g_basis(bond_length);
        let s_matrix = integrals::overlap_matrix(&basis);
        let h_core = integrals::hcore_matrix(&basis);
        let eri = integrals::eri_compressed(&basis);
        let e_nuc = 1.0 / bond_length; // Z_A*Z_B/R for H₂

        let system = PresetSystem {
            system_id: "h2_cphf_test".to_string(),
            label: "H2".to_string(),
            nbf: basis.n_basis,
            nelec: 2,
            s_matrix,
            h_core,
            eri_compressed: eri,
            e_nuc,
        };

        let config = ScfConfig {
            use_diis: true,
            ..ScfConfig::new(ConvergenceProfile::Tight)
        };
        let output = rhf_scf(&system, &config).expect("SCF should converge for H2");
        (output, basis)
    }

    /// H₂O reference geometry (bohr).
    const H2O_O: [f64; 3] = [0.0, 0.0, 0.0];
    const H2O_H1: [f64; 3] = [0.0, 1.43, 1.11];
    const H2O_H2: [f64; 3] = [0.0, -1.43, 1.11];

    /// Create H₂O STO-3G basis set.
    fn h2o_sto3g_basis() -> BasisSet {
        let atoms = vec![
            Atom::new(8, H2O_O).unwrap(),
            Atom::new(1, H2O_H1).unwrap(),
            Atom::new(1, H2O_H2).unwrap(),
        ];
        BasisSet::build(atoms, "sto-3g").expect("Failed to build H2O STO-3G basis")
    }

    /// Run SCF on H₂O and return output + basis.
    fn h2o_scf() -> (ScfOutput, BasisSet) {
        let basis = h2o_sto3g_basis();
        let s_matrix = integrals::overlap_matrix(&basis);
        let h_core = integrals::hcore_matrix(&basis);
        let eri = integrals::eri_compressed(&basis);

        let system = PresetSystem {
            system_id: "h2o_cphf_test".to_string(),
            label: "H2O".to_string(),
            nbf: basis.n_basis,
            nelec: basis.n_electrons,
            s_matrix,
            h_core,
            eri_compressed: eri,
            e_nuc: basis.nuclear_repulsion,
        };

        let config = ScfConfig {
            use_diis: true,
            ..ScfConfig::new(ConvergenceProfile::Tight)
        };
        let output = rhf_scf(&system, &config).expect("SCF should converge for H2O");
        (output, basis)
    }

    /// Extract CPHF inputs from SCF output.
    ///
    /// Returns (mo_coeff, mo_energies, density, eri_compressed, n_occ)
    /// following the same convention as hessian.rs tests.
    fn cphf_inputs_from_scf(
        output: &ScfOutput,
        basis: &BasisSet,
    ) -> (DMatrix<f64>, Vec<f64>, DMatrix<f64>, Vec<f64>, usize) {
        let nbf = basis.n_basis;
        let n_occ = basis.n_electrons / 2;

        // Reconstruct DMatrix from ScfOutput (column-major from nalgebra as_slice)
        // Using from_row_slice to match convention in hessian.rs tests
        let mo_coeff = DMatrix::from_row_slice(nbf, nbf, &output.mo_coefficients);
        let density = DMatrix::from_row_slice(nbf, nbf, &output.density_matrix);

        let eri = integrals::eri_compressed(basis);

        (mo_coeff, output.mo_energies.clone(), density, eri, n_occ)
    }

    /// Truncate full nmo x nmo MO perturbation matrices to nmo x nocc.
    ///
    /// The CPHF solver needs only the occupied columns.
    fn truncate_to_nocc(
        matrices: &[[DMatrix<f64>; 3]],
        _nmo: usize,
        n_occ: usize,
    ) -> Vec<[DMatrix<f64>; 3]> {
        matrices
            .iter()
            .map(|dirs| {
                [
                    dirs[0].columns(0, n_occ).clone_owned(),
                    dirs[1].columns(0, n_occ).clone_owned(),
                    dirs[2].columns(0, n_occ).clone_owned(),
                ]
            })
            .collect()
    }

    // ========================================================================
    // Test 1: Krylov solver on simple known system
    // ========================================================================

    #[test]
    fn test_krylov_simple_system() {
        // Solve (1+A)x = b where A is a small 5x5 symmetric positive matrix
        let n = 5;
        let a_data = vec![
            0.01, 0.002, 0.003, 0.001, 0.002, 0.002, 0.015, 0.001, 0.003, 0.001, 0.003, 0.001,
            0.02, 0.002, 0.003, 0.001, 0.003, 0.002, 0.012, 0.002, 0.002, 0.001, 0.003, 0.002,
            0.018,
        ];
        let a_mat = DMatrix::from_row_slice(n, n, &a_data);

        let b = vec![1.0, 2.0, 3.0, 4.0, 5.0];

        // Direct solve for reference: (I + A)x = b
        let ia = DMatrix::identity(n, n) + &a_mat;
        let x_ref = ia
            .lu()
            .solve(&DMatrix::from_column_slice(n, 1, &b))
            .unwrap();

        // Krylov solve
        let a_for_krylov = a_mat.clone();
        let mut aop = move |x: &[f64]| -> Vec<f64> {
            let xv = DMatrix::from_column_slice(n, 1, x);
            let ax = &a_for_krylov * &xv;
            ax.as_slice().to_vec()
        };

        let (x_krylov, iters, conv) = krylov_solve(&mut aop, &b, 1, 1e-10, 30);
        assert!(conv, "Krylov should converge");
        assert!(iters <= 10, "Should converge quickly: iters = {}", iters);

        for i in 0..n {
            assert!(
                (x_krylov[i] - x_ref[(i, 0)]).abs() < 1e-6,
                "Krylov x[{}] = {:.10e}, ref = {:.10e}",
                i,
                x_krylov[i],
                x_ref[(i, 0)]
            );
        }
    }

    #[test]
    fn test_krylov_multiple_roots() {
        // Solve (1+A)x = b with 3 RHS vectors simultaneously
        let n = 5;
        let a_data = vec![
            0.01, 0.002, 0.003, 0.001, 0.002, 0.002, 0.015, 0.001, 0.003, 0.001, 0.003, 0.001,
            0.02, 0.002, 0.003, 0.001, 0.003, 0.002, 0.012, 0.002, 0.002, 0.001, 0.003, 0.002,
            0.018,
        ];
        let a_mat = DMatrix::from_row_slice(n, n, &a_data);
        let ia = DMatrix::identity(n, n) + &a_mat;

        let nroots = 3;
        let b_all = vec![
            1.0, 2.0, 3.0, 4.0, 5.0, // root 0
            5.0, 4.0, 3.0, 2.0, 1.0, // root 1
            1.0, 0.0, 1.0, 0.0, 1.0, // root 2
        ];

        let a_for_krylov = a_mat.clone();
        let mut aop = move |x: &[f64]| -> Vec<f64> {
            let mut result = vec![0.0; x.len()];
            for r in 0..x.len() / n {
                let xv = DMatrix::from_column_slice(n, 1, &x[r * n..(r + 1) * n]);
                let ax = &a_for_krylov * &xv;
                result[r * n..(r + 1) * n].copy_from_slice(ax.as_slice());
            }
            result
        };

        let (x_krylov, _iters, conv) = krylov_solve(&mut aop, &b_all, nroots, 1e-10, 30);
        assert!(conv, "Krylov should converge");

        for r in 0..nroots {
            let bv = DMatrix::from_column_slice(n, 1, &b_all[r * n..(r + 1) * n]);
            let x_ref = ia.clone().lu().solve(&bv).unwrap();
            for i in 0..n {
                assert!(
                    (x_krylov[r * n + i] - x_ref[(i, 0)]).abs() < 1e-8,
                    "Root {} x[{}] = {:.10e}, ref = {:.10e}",
                    r,
                    i,
                    x_krylov[r * n + i],
                    x_ref[(i, 0)]
                );
            }
        }
    }

    // ========================================================================
    // Test 2: CPHF convergence for H₂/STO-3G (AC7: <=2 iterations)
    // ========================================================================

    #[test]
    fn test_cphf_h2_sto3g_convergence() {
        let (output, basis) = h2_scf(1.4);
        let (mo_coeff, mo_energies, density, eri, n_occ) = cphf_inputs_from_scf(&output, &basis);
        let nbf = basis.n_basis;

        // Build H¹ and S¹
        let h1_ao = make_h1(&basis, &density, 1.0);
        let s1_ao = make_s1(&basis);

        // Transform to MO basis
        let h1_mo = ao_to_mo(&h1_ao, &mo_coeff);
        let s1_mo = ao_to_mo(&s1_ao, &mo_coeff);

        // Truncate to nmo x nocc
        let h1_cphf = truncate_to_nocc(&h1_mo, nbf, n_occ);
        let s1_cphf = truncate_to_nocc(&s1_mo, nbf, n_occ);

        // Build vind
        let mut vind = gen_vind_rhf(&mo_coeff, n_occ, &eri, nbf, 1.0);

        // Flatten perturbations
        let h1_flat: Vec<DMatrix<f64>> = h1_cphf.iter().flat_map(|d| d.iter().cloned()).collect();
        let s1_flat: Vec<DMatrix<f64>> = s1_cphf.iter().flat_map(|d| d.iter().cloned()).collect();

        let config = CphfConfig {
            tol: 1e-9 * basis.atoms.len() as f64,
            ..CphfConfig::default()
        };

        let (mo1, _mo_e1, iterations, converged) =
            cphf_solve_withs1(&mut vind, &mo_energies, n_occ, &h1_flat, &s1_flat, &config);

        assert!(converged, "CPHF should converge for H2");
        assert!(
            iterations <= 2,
            "H2/STO-3G CPHF should converge in <=2 iterations, got {}",
            iterations
        );
        eprintln!("H2/STO-3G CPHF converged in {} iterations", iterations);

        // Verify results are finite
        for (k, m) in mo1.iter().enumerate() {
            assert!(
                m.iter().all(|&v| v.is_finite()),
                "mo1[{}] has non-finite values",
                k
            );
        }
    }

    // ========================================================================
    // Test 3: CPHF convergence for H₂O/STO-3G (AC8: <=20 iterations)
    // ========================================================================

    #[test]
    fn test_cphf_h2o_sto3g_convergence() {
        let (output, basis) = h2o_scf();
        let (mo_coeff, mo_energies, density, eri, n_occ) = cphf_inputs_from_scf(&output, &basis);
        let nbf = basis.n_basis;

        let h1_ao = make_h1(&basis, &density, 1.0);
        let s1_ao = make_s1(&basis);
        let h1_mo = ao_to_mo(&h1_ao, &mo_coeff);
        let s1_mo = ao_to_mo(&s1_ao, &mo_coeff);

        let h1_cphf = truncate_to_nocc(&h1_mo, nbf, n_occ);
        let s1_cphf = truncate_to_nocc(&s1_mo, nbf, n_occ);

        let mut vind = gen_vind_rhf(&mo_coeff, n_occ, &eri, nbf, 1.0);

        let h1_flat: Vec<DMatrix<f64>> = h1_cphf.iter().flat_map(|d| d.iter().cloned()).collect();
        let s1_flat: Vec<DMatrix<f64>> = s1_cphf.iter().flat_map(|d| d.iter().cloned()).collect();

        let config = CphfConfig {
            tol: 1e-9 * basis.atoms.len() as f64,
            ..CphfConfig::default()
        };

        let (_mo1, _mo_e1, iterations, converged) =
            cphf_solve_withs1(&mut vind, &mo_energies, n_occ, &h1_flat, &s1_flat, &config);

        assert!(converged, "CPHF should converge for H2O/STO-3G");
        assert!(
            iterations <= 20,
            "H2O/STO-3G CPHF should converge within 20 iterations, got {}",
            iterations
        );
        eprintln!("H2O/STO-3G CPHF converged in {} iterations", iterations);
    }

    // ========================================================================
    // Test 4: Occupied block of mo1 = -s1/2 exactly (AC1 sub-check)
    // ========================================================================

    #[test]
    fn test_cphf_h2_occupied_block() {
        let (output, basis) = h2_scf(1.4);
        let (mo_coeff, mo_energies, density, eri, n_occ) = cphf_inputs_from_scf(&output, &basis);
        let nbf = basis.n_basis;

        let h1_ao = make_h1(&basis, &density, 1.0);
        let s1_ao = make_s1(&basis);
        let h1_mo = ao_to_mo(&h1_ao, &mo_coeff);
        let s1_mo = ao_to_mo(&s1_ao, &mo_coeff);

        let h1_cphf = truncate_to_nocc(&h1_mo, nbf, n_occ);
        let s1_cphf = truncate_to_nocc(&s1_mo, nbf, n_occ);

        let mut vind = gen_vind_rhf(&mo_coeff, n_occ, &eri, nbf, 1.0);

        let h1_flat: Vec<DMatrix<f64>> = h1_cphf.iter().flat_map(|d| d.iter().cloned()).collect();
        let s1_flat: Vec<DMatrix<f64>> = s1_cphf.iter().flat_map(|d| d.iter().cloned()).collect();

        let config = CphfConfig {
            tol: 1e-9 * basis.atoms.len() as f64,
            ..CphfConfig::default()
        };

        let (mo1, _mo_e1, _iters, converged) =
            cphf_solve_withs1(&mut vind, &mo_energies, n_occ, &h1_flat, &s1_flat, &config);
        assert!(converged);

        // Verify occupied block: U^(1)_ij = -S^(1)_ij / 2
        for k in 0..mo1.len() {
            for oi in 0..n_occ {
                for oj in 0..n_occ {
                    let expected = -s1_flat[k][(oi, oj)] * 0.5;
                    let actual = mo1[k][(oi, oj)];
                    assert!(
                        (actual - expected).abs() < 1e-14,
                        "mo1[{}] occ block ({},{}) = {:.6e}, expected = {:.6e}",
                        k,
                        oi,
                        oj,
                        actual,
                        expected
                    );
                }
            }
        }
    }

    // ========================================================================
    // Test 5: solve_nos1 produces virtual-occupied block (AC2)
    // ========================================================================

    #[test]
    fn test_cphf_nos1_h2() {
        let (output, basis) = h2_scf(1.4);
        let (mo_coeff, mo_energies, _density, eri, n_occ) = cphf_inputs_from_scf(&output, &basis);
        let nbf = basis.n_basis;
        let nmo = nbf;
        let nvir = nmo - n_occ;

        // Create simple perturbation in MO basis (virtual-occupied block)
        let mut h1_vo = vec![DMatrix::zeros(nvir, n_occ); 3];
        h1_vo[0][(0, 0)] = 0.1;
        h1_vo[1][(0, 0)] = 0.2;
        h1_vo[2][(0, 0)] = 0.3;

        let mut vind = gen_vind_rhf_vo(&mo_coeff, n_occ, &eri, nbf, 1.0);

        let config = CphfConfig::default();
        let (mo1, iterations, converged) =
            cphf_solve_nos1(&mut vind, &mo_energies, n_occ, &h1_vo, &config);

        assert!(converged, "CPHF nos1 should converge for H2");
        assert!(
            iterations <= 5,
            "H2/STO-3G nos1 should converge quickly, got {}",
            iterations
        );
        eprintln!("H2/STO-3G CPHF nos1 converged in {} iterations", iterations);

        // Verify shape
        assert_eq!(mo1.len(), 3);
        for m in &mo1 {
            assert_eq!(m.nrows(), nvir);
            assert_eq!(m.ncols(), n_occ);
        }
    }

    // ========================================================================
    // Test 6: Batched vs individual solve (AC10: max diff < 1e-8)
    // ========================================================================

    #[test]
    fn test_cphf_batched_vs_individual() {
        let (output, basis) = h2o_scf();
        let (mo_coeff, mo_energies, density, eri, n_occ) = cphf_inputs_from_scf(&output, &basis);
        let nbf = basis.n_basis;
        let n_atoms = basis.atoms.len();

        let h1_ao = make_h1(&basis, &density, 1.0);
        let s1_ao = make_s1(&basis);
        let h1_mo = ao_to_mo(&h1_ao, &mo_coeff);
        let s1_mo = ao_to_mo(&s1_ao, &mo_coeff);
        let h1_cphf = truncate_to_nocc(&h1_mo, nbf, n_occ);
        let s1_cphf = truncate_to_nocc(&s1_mo, nbf, n_occ);

        // Batched solve: all atoms at once
        let config_batched = CphfConfig::default();
        let mut vind_b = gen_vind_rhf(&mo_coeff, n_occ, &eri, nbf, 1.0);
        let batched = cphf_solve(
            &mut vind_b,
            &mo_energies,
            n_occ,
            &h1_cphf,
            Some(&s1_cphf),
            &config_batched,
        );
        assert!(batched.converged, "Batched CPHF should converge");
        eprintln!(
            "Batched H2O/STO-3G CPHF converged in {} iterations",
            batched.iterations
        );

        // Individual solves: one atom at a time
        let mut individual_mo1: Vec<DMatrix<f64>> = Vec::new();
        for ia in 0..n_atoms {
            let h1_single = vec![h1_cphf[ia].clone()];
            let s1_single = vec![s1_cphf[ia].clone()];

            let h1_flat: Vec<DMatrix<f64>> =
                h1_single.iter().flat_map(|d| d.iter().cloned()).collect();
            let s1_flat: Vec<DMatrix<f64>> =
                s1_single.iter().flat_map(|d| d.iter().cloned()).collect();

            let mut vind_i = gen_vind_rhf(&mo_coeff, n_occ, &eri, nbf, 1.0);
            let config_single = CphfConfig {
                tol: 1e-9,
                ..CphfConfig::default()
            };

            let (mo1_i, _e1_i, _iters, conv) = cphf_solve_withs1(
                &mut vind_i,
                &mo_energies,
                n_occ,
                &h1_flat,
                &s1_flat,
                &config_single,
            );
            assert!(conv, "Individual CPHF should converge for atom {}", ia);
            individual_mo1.extend(mo1_i);
        }

        // Compare
        assert_eq!(batched.mo1.len(), individual_mo1.len());
        let mut max_diff = 0.0f64;
        for (k, (b_m, i_m)) in batched.mo1.iter().zip(individual_mo1.iter()).enumerate() {
            let diff = (b_m - i_m).abs().max();
            if diff > max_diff {
                max_diff = diff;
            }
            assert!(
                diff < 1e-8,
                "Batched vs individual mismatch at perturbation {}: max_diff = {:.2e}",
                k,
                diff
            );
        }
        eprintln!("Batched vs individual max diff: {:.2e}", max_diff);
    }

    // ========================================================================
    // Test 7: DFT vind framework (AC6)
    // ========================================================================

    #[test]
    fn test_cphf_dft_vind_framework() {
        // Verify that gen_vind_dft with xc_response=None gives the same result
        // as gen_vind_rhf, and that providing a mock XC callback modifies the
        // response as expected.

        let (output, basis) = h2_scf(1.4);
        let (mo_coeff, mo_energies, density, eri, n_occ) = cphf_inputs_from_scf(&output, &basis);
        let nbf = basis.n_basis;

        // Part 1: gen_vind_dft with no XC response should equal gen_vind_rhf
        {
            let h1_ao = make_h1(&basis, &density, 1.0);
            let s1_ao = make_s1(&basis);
            let h1_mo = ao_to_mo(&h1_ao, &mo_coeff);
            let s1_mo = ao_to_mo(&s1_ao, &mo_coeff);
            let h1_cphf = truncate_to_nocc(&h1_mo, nbf, n_occ);
            let s1_cphf = truncate_to_nocc(&s1_mo, nbf, n_occ);

            let h1_flat: Vec<DMatrix<f64>> =
                h1_cphf.iter().flat_map(|d| d.iter().cloned()).collect();
            let s1_flat: Vec<DMatrix<f64>> =
                s1_cphf.iter().flat_map(|d| d.iter().cloned()).collect();

            // Solve with gen_vind_rhf (reference)
            let mut vind_rhf = gen_vind_rhf(&mo_coeff, n_occ, &eri, nbf, 1.0);
            let config = CphfConfig {
                tol: 1e-9 * basis.atoms.len() as f64,
                ..CphfConfig::default()
            };
            let (mo1_rhf, _, _, conv_rhf) = cphf_solve_withs1(
                &mut vind_rhf,
                &mo_energies,
                n_occ,
                &h1_flat,
                &s1_flat,
                &config,
            );
            assert!(conv_rhf, "RHF vind should converge");

            // Solve with gen_vind_dft, xc_response=None, hf_exchange=1.0
            let mut vind_dft = gen_vind_dft(&mo_coeff, n_occ, &eri, nbf, 1.0, None);
            let (mo1_dft, _, _, conv_dft) = cphf_solve_withs1(
                &mut vind_dft,
                &mo_energies,
                n_occ,
                &h1_flat,
                &s1_flat,
                &config,
            );
            assert!(conv_dft, "DFT vind (no XC) should converge");

            // They should give identical results
            for (k, (rhf_m, dft_m)) in mo1_rhf.iter().zip(mo1_dft.iter()).enumerate() {
                let diff = (rhf_m - dft_m).abs().max();
                assert!(
                    diff < 1e-12,
                    "RHF vs DFT(no XC) mismatch at perturbation {}: max_diff = {:.2e}",
                    k,
                    diff
                );
            }
            eprintln!("Part 1 PASSED: gen_vind_dft(xc=None) matches gen_vind_rhf");
        }

        // Part 2: gen_vind_dft with a mock XC callback modifies the response
        {
            // Create a mock XC response: a simple scaled-identity kernel
            // V_xc^(1) = alpha * trace(D^(1)) * I
            // This is a crude mock that adds a constant diagonal shift proportional
            // to the density response trace — just to verify the framework calls it.
            let mock_alpha = 0.01;
            let nbf_copy = nbf;
            let mock_xc: Box<XcResponseFn> = Box::new(move |dm1: &DMatrix<f64>| {
                let trace = dm1.trace();
                DMatrix::from_diagonal_element(nbf_copy, nbf_copy, trace * mock_alpha)
            });

            // Build a simple test input: single perturbation
            let nmo = nbf;
            let ndim = nmo * n_occ;
            let mut test_mo1 = vec![0.0; ndim];
            // Put a nonzero value in the virtual-occupied block
            test_mo1[n_occ * n_occ] = 0.5; // first virtual, first occ

            // Evaluate both vind functions with same input
            let mut vind_no_xc = gen_vind_dft(&mo_coeff, n_occ, &eri, nbf, 1.0, None);
            let mut vind_with_xc = gen_vind_dft(&mo_coeff, n_occ, &eri, nbf, 1.0, Some(&*mock_xc));

            let result_no_xc = vind_no_xc(&test_mo1);
            let result_with_xc = vind_with_xc(&test_mo1);

            // Results should differ due to XC contribution
            let diff_norm: f64 = result_no_xc
                .iter()
                .zip(result_with_xc.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!(
                diff_norm > 1e-10,
                "XC response should modify vind output; diff_norm = {:.2e}",
                diff_norm
            );
            eprintln!(
                "Part 2 PASSED: XC callback modifies response (diff_norm = {:.6e})",
                diff_norm
            );
        }

        // Part 3: gen_vind_dft with scaled HF exchange (hybrid DFT scenario)
        {
            let nmo = nbf;
            let ndim = nmo * n_occ;
            let mut test_mo1 = vec![0.0; ndim];
            test_mo1[n_occ * n_occ] = 0.5;

            // Full HF exchange
            let mut vind_full_hf = gen_vind_dft(&mo_coeff, n_occ, &eri, nbf, 1.0, None);
            let result_full = vind_full_hf(&test_mo1);

            // Scaled HF exchange (B3LYP-like: 0.2)
            let mut vind_hybrid = gen_vind_dft(&mo_coeff, n_occ, &eri, nbf, 0.2, None);
            let result_hybrid = vind_hybrid(&test_mo1);

            // No HF exchange (pure DFT)
            let mut vind_pure = gen_vind_dft(&mo_coeff, n_occ, &eri, nbf, 0.0, None);
            let result_pure = vind_pure(&test_mo1);

            // All should differ (K term scales differently)
            let diff_full_hybrid: f64 = result_full
                .iter()
                .zip(result_hybrid.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            let diff_full_pure: f64 = result_full
                .iter()
                .zip(result_pure.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();

            assert!(
                diff_full_hybrid > 1e-10,
                "Full vs hybrid should differ: {:.2e}",
                diff_full_hybrid
            );
            assert!(
                diff_full_pure > 1e-10,
                "Full vs pure DFT should differ: {:.2e}",
                diff_full_pure
            );
            // pure DFT has more difference from full HF than hybrid
            assert!(
                diff_full_pure > diff_full_hybrid,
                "Pure DFT should differ more from full HF than hybrid"
            );
            eprintln!(
                "Part 3 PASSED: HF exchange scaling works \
                 (full-hybrid = {:.6e}, full-pure = {:.6e})",
                diff_full_hybrid, diff_full_pure
            );
        }

        // Part 4: gen_vind_dft_vo (virtual-occupied variant)
        {
            let nmo = nbf;
            let nvir = nmo - n_occ;
            let ndim = nvir * n_occ;
            let mut test_mo1 = vec![0.0; ndim];
            test_mo1[0] = 0.5;

            // Compare RHF vo variant with DFT vo variant (no XC)
            let mut vind_rhf_vo = gen_vind_rhf_vo(&mo_coeff, n_occ, &eri, nbf, 1.0);
            let mut vind_dft_vo = gen_vind_dft_vo(&mo_coeff, n_occ, &eri, nbf, 1.0, None);

            let result_rhf = vind_rhf_vo(&test_mo1);
            let result_dft = vind_dft_vo(&test_mo1);

            let diff: f64 = result_rhf
                .iter()
                .zip(result_dft.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!(
                diff < 1e-12,
                "RHF_vo vs DFT_vo(no XC) should be identical: diff = {:.2e}",
                diff
            );
            eprintln!("Part 4 PASSED: gen_vind_dft_vo matches gen_vind_rhf_vo");
        }
    }

    // ========================================================================
    // Test 8: H₂O/6-31G* CPHF convergence (AC9: <=30 iterations)
    // ========================================================================

    #[test]
    fn test_cphf_h2o_631gs_convergence() {
        // Build H₂O at the specified geometry
        let atoms = vec![
            Atom::new(8, [0.0, 0.0, 0.0]).unwrap(),    // O
            Atom::new(1, [0.0, 1.43, 1.11]).unwrap(),  // H1
            Atom::new(1, [0.0, -1.43, 1.11]).unwrap(), // H2
        ];
        let basis = BasisSet::build(atoms, "6-31g*").expect("Failed to build H2O/6-31G* basis");

        eprintln!(
            "H2O/6-31G*: nbf={}, n_electrons={}, n_occ={}",
            basis.n_basis,
            basis.n_electrons,
            basis.n_occupied()
        );

        // Compute integrals
        let s_matrix = integrals::overlap_matrix(&basis);
        let h_core = integrals::hcore_matrix(&basis);
        let eri = integrals::eri_compressed(&basis);

        let system = PresetSystem {
            system_id: "h2o_631gs_cphf_test".to_string(),
            label: "H2O/6-31G*".to_string(),
            nbf: basis.n_basis,
            nelec: basis.n_electrons,
            s_matrix,
            h_core,
            eri_compressed: eri.clone(),
            e_nuc: basis.nuclear_repulsion,
        };

        // Run SCF
        let config = ScfConfig {
            use_diis: true,
            ..ScfConfig::new(ConvergenceProfile::Tight)
        };
        let output = rhf_scf(&system, &config).expect("SCF should converge for H2O/6-31G*");
        assert!(output.converged, "H2O/6-31G* SCF should converge");
        eprintln!(
            "H2O/6-31G* SCF converged in {} iterations, E = {:.10}",
            output.iterations, output.energy_total
        );

        // Extract CPHF inputs
        let nbf = basis.n_basis;
        let n_occ = basis.n_electrons / 2;
        let mo_coeff = DMatrix::from_row_slice(nbf, nbf, &output.mo_coefficients);
        let density = DMatrix::from_row_slice(nbf, nbf, &output.density_matrix);

        // Build H¹ and S¹
        let h1_ao = make_h1(&basis, &density, 1.0);
        let s1_ao = make_s1(&basis);

        // Transform to MO basis
        let h1_mo = ao_to_mo(&h1_ao, &mo_coeff);
        let s1_mo = ao_to_mo(&s1_ao, &mo_coeff);

        // Truncate to nmo x nocc
        let h1_cphf = truncate_to_nocc(&h1_mo, nbf, n_occ);
        let s1_cphf = truncate_to_nocc(&s1_mo, nbf, n_occ);

        // Build vind
        let mut vind = gen_vind_rhf(&mo_coeff, n_occ, &eri, nbf, 1.0);

        // Flatten perturbations
        let h1_flat: Vec<DMatrix<f64>> = h1_cphf.iter().flat_map(|d| d.iter().cloned()).collect();
        let s1_flat: Vec<DMatrix<f64>> = s1_cphf.iter().flat_map(|d| d.iter().cloned()).collect();

        // Scale tolerance by number of atoms (as in batched solve)
        let cphf_config = CphfConfig {
            tol: 1e-9 * basis.atoms.len() as f64,
            ..CphfConfig::default()
        };

        let (mo1, _mo_e1, iterations, converged) = cphf_solve_withs1(
            &mut vind,
            &output.mo_energies,
            n_occ,
            &h1_flat,
            &s1_flat,
            &cphf_config,
        );

        eprintln!(
            "H2O/6-31G* CPHF: converged={}, iterations={}",
            converged, iterations
        );

        assert!(converged, "CPHF should converge for H2O/6-31G*");
        assert!(
            iterations <= 30,
            "H2O/6-31G* CPHF should converge within 30 iterations, got {}",
            iterations
        );

        // Verify results are finite and non-trivial
        let mut max_abs_mo1 = 0.0f64;
        for (k, m) in mo1.iter().enumerate() {
            assert!(
                m.iter().all(|&v| v.is_finite()),
                "mo1[{}] has non-finite values",
                k
            );
            let this_max = m.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
            if this_max > max_abs_mo1 {
                max_abs_mo1 = this_max;
            }
        }
        assert!(
            max_abs_mo1 > 1e-10,
            "mo1 values should be non-trivial, max = {:.2e}",
            max_abs_mo1
        );
        eprintln!("max |mo1| = {:.6e}", max_abs_mo1);
    }
}
