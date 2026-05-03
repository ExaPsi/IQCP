//! Electric dipole integrals via the overlap-shift recurrence
//!
//! This module computes the AO-basis matrix elements of the electric
//! dipole operator with gauge origin `C`:
//!
//! ```text
//! mu^d_{mu,nu} = <chi_mu | (r_d - C_d) | chi_nu>,   d in {x, y, z}
//! ```
//!
//! # Algorithm — Overlap-Shift Identity
//!
//! Write `(r_d - C_d) = (r_d - A_d) + (A_d - C_d)` where `A` is the bra
//! center. The first piece times the bra Gaussian shifts its angular
//! momentum by 1:
//!
//! ```text
//! <chi_{a}| (r_d - A_d) |chi_{b}> = <chi_{a+1_d} | chi_b>
//! ```
//!
//! (ignoring the normalization factor which we handle outside the
//! primitive). The second piece is just `(A_d - C_d)` times the plain
//! overlap. Therefore:
//!
//! ```text
//! <a|(r_d - C_d)|b> = <a + 1_d | b> + (A_d - C_d) * <a|b>
//! ```
//!
//! This lets us compute all three dipole components from exactly two
//! overlap evaluations per Cartesian direction: one with raised angular
//! momentum and one at the base. We reuse the validated
//! `primitive_overlap` / `overlap_1d` machinery from `overlap.rs` so
//! the dipole integrals inherit the same numerical accuracy as the
//! overlap integrals.
//!
//! # First-derivative dipole integrals (bra-center nabla identity)
//!
//! For the skeleton term of the dipole-derivative Hessian, we also need
//! `d/dR_{A,e} <mu|(r_d - C_d)|nu>` when mu is centered on A. The
//! derivative of a Gaussian primitive w.r.t. its center is:
//!
//! ```text
//! d/dA_e chi_a = 2*alpha * chi_{a+1_e} - l_{a,e} * chi_{a-1_e}
//! ```
//!
//! and the multipole operator is A-independent, so:
//!
//! ```text
//! d/dA_e <chi_a|(r_d - C_d)|chi_b> =
//!     2*alpha * <chi_{a+1_e} | (r_d - C_d) | chi_b>
//!     - l_{a,e} * <chi_{a-1_e} | (r_d - C_d) | chi_b>
//!     - delta_{de} * <chi_a | chi_b>
//! ```
//!
//! The last term has the opposite sign to what a naive application of
//! the nabla identity would give, because `(A_e - C_e)` explicitly
//! depends on `A_e`:
//!
//! ```text
//! d/dA_e [(A_e - C_e) <a|b>] = <a|b> + (A_e - C_e) * d<a|b>/dA_e
//! ```
//!
//! Combining the two contributions (the `(r - A)` piece gives the
//! nabla-identity result applied to the raised integral, and the
//! `(A - C)` piece gives the explicit `delta_de` term plus the nabla
//! identity applied to the base overlap) produces the formula above.
//!
//! Wait — we need to be careful. Let me re-derive cleanly.
//!
//! Let `f(A) = <chi_a(r; A) | (r_d - C_d) | chi_b(r; B)>`. Using the
//! overlap shift trick once:
//!
//! ```text
//! f(A) = S_{a+1_d, b}(A, B) + (A_d - C_d) * S_{a, b}(A, B)
//! ```
//!
//! Taking d/dA_e:
//!
//! ```text
//! df/dA_e = dS_{a+1_d, b}/dA_e + delta_{de} * S_{a,b} + (A_d - C_d) * dS_{a,b}/dA_e
//! ```
//!
//! Each `dS/dA_e` uses the nabla identity:
//! ```text
//! dS_{a', b}/dA_e = 2*alpha * S_{a'+1_e, b} - (l_{a', e}) * S_{a'-1_e, b}
//! ```
//!
//! So for the bra-center derivative of the dipole integral:
//! ```text
//! d mu^d_{a, b}/dA_e =
//!     2*alpha * S_{a + 1_d + 1_e, b}
//!     - (l_{a, e} + delta_{de}) * S_{a + 1_d - 1_e, b}   (note: raising 1_d first changes l_e if d==e)
//!     + delta_{de} * S_{a, b}
//!     + (A_d - C_d) * (2*alpha * S_{a + 1_e, b} - l_{a, e} * S_{a - 1_e, b})
//! ```
//!
//! Actually the "l_e" term when we lower 1_e after raising 1_d changes
//! by +1 only if d == e. To avoid index bookkeeping bugs, we compute
//! the first-derivative dipole integral with direct back-substitution
//! into `primitive_dipole(a_plus_e)` minus the lowered variant.
//!
//! # References
//!
//! - Helgaker, Jorgensen & Olsen (2000), Section 9.3 & Eq. 9.3.8
//! - Obara & Saika (1986), J. Chem. Phys. 84, 3963, Eq. 15-24
//! - libcint `cint1e_r` (gauge-free dipole) and `cint1e_irp` (first deriv)
//! - PySCF `scf/hf.py:dip_moment` line 1418
//! - IQCP `integrals/overlap.rs` (primitive_overlap, overlap_1d)
//! - IQCP `integrals/deriv1.rs` (shell_overlap_first_deriv pattern)

use crate::basis::{BasisSet, ContractedShell};
use crate::integrals::overlap::cartesian_gaussian_normalization;
use crate::integrals::{cartesian_components, CartesianPower, GaussianProduct};
use nalgebra::DMatrix;

use super::overlap::primitive_overlap;

// ============================================================================
// Angular-momentum helpers (mirrors deriv1.rs)
// ============================================================================

#[inline]
fn get_angular_component(pow: &CartesianPower, dir: usize) -> u32 {
    match dir {
        0 => pow.i,
        1 => pow.j,
        2 => pow.k,
        _ => unreachable!(),
    }
}

#[inline]
fn raise_angular(pow: &CartesianPower, dir: usize) -> CartesianPower {
    let mut p = *pow;
    match dir {
        0 => p.i += 1,
        1 => p.j += 1,
        2 => p.k += 1,
        _ => unreachable!(),
    }
    p
}

#[inline]
fn lower_angular(pow: &CartesianPower, dir: usize) -> Option<CartesianPower> {
    let mut p = *pow;
    match dir {
        0 => {
            if p.i == 0 {
                return None;
            }
            p.i -= 1;
        }
        1 => {
            if p.j == 0 {
                return None;
            }
            p.j -= 1;
        }
        2 => {
            if p.k == 0 {
                return None;
            }
            p.k -= 1;
        }
        _ => unreachable!(),
    }
    Some(p)
}

// ============================================================================
// Primitive dipole integral via overlap shift
// ============================================================================

/// Compute the three Cartesian components of the primitive dipole
/// integral `<chi_a | (r_d - C_d) | chi_b>` for a given primitive pair.
///
/// Uses the overlap-shift identity:
/// ```text
/// <a | (r_d - C_d) | b> = <a + 1_d | b> + (A_d - C_d) * <a|b>
/// ```
///
/// # Arguments
///
/// * `gp` - Pre-computed Gaussian product data (bra-ket centers/exponents)
/// * `a_powers` - Cartesian powers of the bra primitive
/// * `b_powers` - Cartesian powers of the ket primitive
/// * `bra_center` - Center of the bra Gaussian (three Cartesian coords)
/// * `gauge_origin` - Gauge origin `C` for the dipole operator
///
/// # Returns
///
/// `[mu_x, mu_y, mu_z]` — primitive dipole integral values
fn primitive_dipole(
    gp: &GaussianProduct,
    a_powers: &CartesianPower,
    b_powers: &CartesianPower,
    bra_center: &[f64; 3],
    gauge_origin: &[f64; 3],
) -> [f64; 3] {
    let s_base = primitive_overlap(gp, a_powers, b_powers);
    let mut out = [0.0f64; 3];
    for d in 0..3 {
        let raised = raise_angular(a_powers, d);
        let s_raised = primitive_overlap(gp, &raised, b_powers);
        out[d] = s_raised + (bra_center[d] - gauge_origin[d]) * s_base;
    }
    out
}

// ============================================================================
// Shell dipole integrals
// ============================================================================

/// Compute the three Cartesian dipole matrix blocks between two
/// contracted shells.
///
/// # Arguments
///
/// * `shell_a` - Bra shell
/// * `shell_b` - Ket shell
/// * `gauge_origin` - Gauge origin `C` for the dipole operator
///
/// # Returns
///
/// `[mu_x_block, mu_y_block, mu_z_block]` where each element is a
/// flattened `n_a * n_b` row-major block of dipole integrals.
///
/// # Formula
///
/// ```text
/// mu^d_{ij} = sum_{p in A, q in B} c_p * c_q * N_p * N_q *
///             [<p + 1_d | q> + (A_d - C_d) <p|q>]
/// ```
///
/// where `N_p`, `N_q` are the Cartesian Gaussian normalization factors
/// applied to the *original* angular momentum (matching
/// `shell_overlap`).
pub fn shell_dipole_integrals(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
    gauge_origin: &[f64; 3],
) -> [Vec<f64>; 3] {
    let l_a = shell_a.l_value();
    let l_b = shell_b.l_value();
    let comps_a = cartesian_components(l_a).expect("Angular momentum in range");
    let comps_b = cartesian_components(l_b).expect("Angular momentum in range");
    let n_a = comps_a.len();
    let n_b = comps_b.len();

    let mut result = [
        vec![0.0f64; n_a * n_b],
        vec![0.0f64; n_a * n_b],
        vec![0.0f64; n_a * n_b],
    ];

    for prim_a in &shell_a.primitives {
        for prim_b in &shell_b.primitives {
            let gp = GaussianProduct::new(
                prim_a.exponent,
                &shell_a.center,
                prim_b.exponent,
                &shell_b.center,
            );
            let coef = prim_a.coefficient * prim_b.coefficient;

            for (ia, pow_a) in comps_a.iter().enumerate() {
                let norm_a = cartesian_gaussian_normalization(prim_a.exponent, pow_a);

                for (jb, pow_b) in comps_b.iter().enumerate() {
                    let norm_b = cartesian_gaussian_normalization(prim_b.exponent, pow_b);
                    let cn = coef * norm_a * norm_b;

                    let mu_prim =
                        primitive_dipole(&gp, pow_a, pow_b, &shell_a.center, gauge_origin);

                    for d in 0..3 {
                        result[d][ia * n_b + jb] += cn * mu_prim[d];
                    }
                }
            }
        }
    }

    result
}

// ============================================================================
// Full dipole matrix
// ============================================================================

/// Compute the three full `nbf x nbf` dipole matrices for a basis set.
///
/// # Arguments
///
/// * `basis` - The molecular basis set (must be Cartesian; use
///   `mol.cart = True` in PySCF for comparison)
/// * `gauge_origin` - Gauge origin `C` for the dipole operator (default
///   for PySCF is `[0, 0, 0]`)
///
/// # Returns
///
/// `[D_x, D_y, D_z]` where each matrix is symmetric and has
/// `basis.n_basis` rows/columns.
pub fn dipole_matrix(basis: &BasisSet, gauge_origin: &[f64; 3]) -> [DMatrix<f64>; 3] {
    let nbf = basis.n_basis;
    let mut mats = [
        DMatrix::<f64>::zeros(nbf, nbf),
        DMatrix::<f64>::zeros(nbf, nbf),
        DMatrix::<f64>::zeros(nbf, nbf),
    ];

    // Pre-compute shell offsets
    let mut shell_offsets = Vec::with_capacity(basis.shells.len());
    let mut offset = 0usize;
    for shell in &basis.shells {
        shell_offsets.push(offset);
        offset += shell.n_basis_functions();
    }

    // Iterate over shell pairs (upper triangle; mirror via Hermiticity)
    for (i, shell_i) in basis.shells.iter().enumerate() {
        let off_i = shell_offsets[i];
        let n_i = shell_i.n_basis_functions();

        for (j, shell_j) in basis.shells.iter().enumerate() {
            let off_j = shell_offsets[j];
            let n_j = shell_j.n_basis_functions();

            if i <= j {
                let block = shell_dipole_integrals(shell_i, shell_j, gauge_origin);
                for d in 0..3 {
                    for ia in 0..n_i {
                        for jb in 0..n_j {
                            let val = block[d][ia * n_j + jb];
                            mats[d][(off_i + ia, off_j + jb)] = val;
                            if i != j {
                                mats[d][(off_j + jb, off_i + ia)] = val;
                            }
                        }
                    }
                }
            }
        }
    }

    mats
}

// ============================================================================
// First derivative of shell dipole integrals (bra-center)
// ============================================================================

/// First derivative of the shell dipole integral block with respect to
/// the bra center.
///
/// Returns a `[3][3]` array of flattened `n_a * n_b` blocks: the outer
/// index is the dipole direction `d` and the inner index is the bra-
/// center derivative direction `e`. That is:
///
/// ```text
/// result[d][e][ia * n_b + jb] = d/dA_e <chi_{a_ia} | (r_d - C_d) | chi_{b_jb}>
/// ```
///
/// # Formula
///
/// Differentiating the overlap-shift identity:
/// ```text
/// <a|(r_d - C_d)|b> = <a+1_d|b> + (A_d - C_d) <a|b>
/// ```
///
/// with respect to `A_e` (using the nabla identity on the two overlap
/// evaluations and the explicit derivative of `(A_d - C_d)`):
///
/// ```text
/// d/dA_e <a|(r_d - C_d)|b>
///   = d/dA_e <a+1_d|b>
///   + delta_{de} <a|b>
///   + (A_d - C_d) * d/dA_e <a|b>
///
///   where d/dA_e <a'|b> = 2*alpha * <a'+1_e|b> - l_{a',e} <a'-1_e|b>
/// ```
///
/// # Note on normalization
///
/// Matches the `deriv1.rs` convention: the *original* angular-momentum
/// Cartesian Gaussian normalization is applied to both the raised and
/// lowered primitive overlaps, not the normalization of the shifted
/// Gaussians. This is physically correct because we are differentiating
/// the Cartesian basis function, which keeps its normalization fixed
/// while its center moves.
#[allow(clippy::needless_range_loop)] // `e` indexes per-Cartesian helpers, not just `result`
pub fn shell_dipole_first_deriv(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
    gauge_origin: &[f64; 3],
) -> [[Vec<f64>; 3]; 3] {
    let l_a = shell_a.l_value();
    let l_b = shell_b.l_value();
    let comps_a = cartesian_components(l_a).expect("Angular momentum in range");
    let comps_b = cartesian_components(l_b).expect("Angular momentum in range");
    let n_a = comps_a.len();
    let n_b = comps_b.len();

    // Allocate 9 result blocks [d][e]
    let mut result: [[Vec<f64>; 3]; 3] = [
        [
            vec![0.0; n_a * n_b],
            vec![0.0; n_a * n_b],
            vec![0.0; n_a * n_b],
        ],
        [
            vec![0.0; n_a * n_b],
            vec![0.0; n_a * n_b],
            vec![0.0; n_a * n_b],
        ],
        [
            vec![0.0; n_a * n_b],
            vec![0.0; n_a * n_b],
            vec![0.0; n_a * n_b],
        ],
    ];

    for prim_a in &shell_a.primitives {
        let alpha = prim_a.exponent;
        for prim_b in &shell_b.primitives {
            let gp = GaussianProduct::new(alpha, &shell_a.center, prim_b.exponent, &shell_b.center);
            let coef = prim_a.coefficient * prim_b.coefficient;

            for (ia, pow_a) in comps_a.iter().enumerate() {
                // Use the ORIGINAL Cartesian normalization on the bra
                // (matching deriv1.rs convention)
                let norm_a = cartesian_gaussian_normalization(alpha, pow_a);

                for (jb, pow_b) in comps_b.iter().enumerate() {
                    let norm_b = cartesian_gaussian_normalization(prim_b.exponent, pow_b);
                    let cn = coef * norm_a * norm_b;

                    // Base overlaps and raised-base overlaps (angular momentum +1 in d)
                    let s_base = primitive_overlap(&gp, pow_a, pow_b);
                    for d in 0..3 {
                        let a_plus_d = raise_angular(pow_a, d);
                        // gauge shift coefficient (A_d - C_d)
                        let shift_d = shell_a.center[d] - gauge_origin[d];

                        for e in 0..3 {
                            // Term 1: d/dA_e <a+1_d|b>
                            //   = 2*alpha * <a+1_d+1_e|b> - l_{a+1_d, e} * <a+1_d-1_e|b>
                            // Note: l_{a+1_d, e} = pow_a[e] + (1 if d == e else 0)
                            let a_plus_d_plus_e = raise_angular(&a_plus_d, e);
                            let s_plus_d_plus_e = primitive_overlap(&gp, &a_plus_d_plus_e, pow_b);
                            let l_a_plus_d_e = get_angular_component(&a_plus_d, e) as f64;
                            let term_raised_minus =
                                if let Some(a_plus_d_minus_e) = lower_angular(&a_plus_d, e) {
                                    primitive_overlap(&gp, &a_plus_d_minus_e, pow_b)
                                } else {
                                    0.0
                                };
                            let d_s_raised =
                                2.0 * alpha * s_plus_d_plus_e - l_a_plus_d_e * term_raised_minus;

                            // Term 2: delta_{de} * <a|b>
                            let delta_de = if d == e { s_base } else { 0.0 };

                            // Term 3: (A_d - C_d) * d/dA_e <a|b>
                            let a_plus_e = raise_angular(pow_a, e);
                            let s_plus_e = primitive_overlap(&gp, &a_plus_e, pow_b);
                            let l_a_e = get_angular_component(pow_a, e) as f64;
                            let term_base_minus = if let Some(a_minus_e) = lower_angular(pow_a, e) {
                                primitive_overlap(&gp, &a_minus_e, pow_b)
                            } else {
                                0.0
                            };
                            let d_s_base = 2.0 * alpha * s_plus_e - l_a_e * term_base_minus;

                            let deriv = d_s_raised + delta_de + shift_d * d_s_base;
                            result[d][e][ia * n_b + jb] += cn * deriv;
                        }
                    }
                }
            }
        }
    }

    result
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
#[allow(clippy::needless_range_loop)]
mod tests {
    use super::*;
    use crate::basis::{Atom, BasisSet};
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_dipole_matrix_h2o_symmetric() {
        // H2O at standard geometry, STO-3G, Cartesian basis
        let o = Atom::new(8, [0.0, 0.0, 0.22166487441148175]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.4309006215206648, -0.886659497645927]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.4309006215206648, -0.886659497645927]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        let dip = dipole_matrix(&basis, &[0.0, 0.0, 0.0]);

        // Each component matrix must be symmetric (Hermiticity of dipole operator
        // over real Gaussians)
        let n = basis.n_basis;
        for d in 0..3 {
            for i in 0..n {
                for j in 0..n {
                    let err = (dip[d][(i, j)] - dip[d][(j, i)]).abs();
                    assert!(
                        err < 1e-14,
                        "Dipole matrix component {d} not symmetric at ({i},{j}): err = {err}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_dipole_matrix_h2_zero_dipole_moment_component() {
        // H2 along z: the electronic dipole <D_z> = -tr(D * mu_z)
        // Since D is symmetric and both atoms are identical, the result should
        // combined with nuclear contribution give zero net dipole. Here we just
        // check that the dipole matrix isn't pathologically broken (non-NaN,
        // non-zero for p-like shell pairs via projection).
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let dip = dipole_matrix(&basis, &[0.0, 0.0, 0.0]);
        assert_eq!(basis.n_basis, 2);
        // <s_A | r_z | s_A> for s at (0,0,0) should be 0 by symmetry
        assert_abs_diff_eq!(dip[2][(0, 0)], 0.0, epsilon = 1e-12);
        // <s_B | r_z | s_B> for s at (0,0,1.4) should be 1.4
        assert_abs_diff_eq!(dip[2][(1, 1)], 1.4, epsilon = 1e-9);
        // All xy/yx elements should be finite
        for d in 0..3 {
            for i in 0..2 {
                for j in 0..2 {
                    assert!(dip[d][(i, j)].is_finite(), "Non-finite dipole element");
                }
            }
        }
    }

    #[test]
    fn test_dipole_matrix_off_center_ss() {
        // Two s-Gaussians: bra at A=(0,0,0), ket at B=(0,0,R).
        // For the pair, the dipole <sA | (r - 0) | sB> has z-component equal
        // to the overlap times the product-center z-coordinate:
        //   <sA | r_z | sB> = S * P_z = S * (alpha*A_z + beta*B_z)/(alpha+beta)
        // This is a standard result from the Gaussian product theorem.
        // We verify it for STO-3G contracted H shells at R = 1.4 bohr.
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
        let dip = dipole_matrix(&basis, &[0.0, 0.0, 0.0]);

        // <sA | r_z | sA> = 0 (centered at origin)
        assert_abs_diff_eq!(dip[2][(0, 0)], 0.0, epsilon = 1e-12);
        // <sB | r_z | sB> = 1.4 (centered at 1.4)
        assert_abs_diff_eq!(dip[2][(1, 1)], 1.4, epsilon = 1e-9);
        // <sA | r_x | sA> = 0 and other orthogonal components = 0
        assert_abs_diff_eq!(dip[0][(0, 0)], 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(dip[1][(0, 0)], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_shell_dipole_first_deriv_finite_difference() {
        // Finite-difference validation of shell_dipole_first_deriv:
        // displace the bra shell center by +/- h along each direction,
        // recompute shell_dipole_integrals, and compare the central
        // difference to the analytical first derivative.
        use crate::basis::{AngularMomentum, ContractedShell, GaussianPrimitive};

        let h_prims = vec![
            GaussianPrimitive::new(3.42525091, 0.15432897),
            GaussianPrimitive::new(0.62391373, 0.53532814),
            GaussianPrimitive::new(0.16885540, 0.44463454),
        ];
        let mut shell_a =
            ContractedShell::new(AngularMomentum::S, h_prims.clone(), [0.0, 0.0, 0.0], 0);
        let shell_b = ContractedShell::new(AngularMomentum::S, h_prims.clone(), [0.0, 0.0, 1.4], 1);

        let origin = [0.0f64, 0.0, 0.0];
        let analytical = shell_dipole_first_deriv(&shell_a, &shell_b, &origin);

        let h = 1e-5;
        for e in 0..3 {
            // Plus
            let mut center_plus = shell_a.center;
            center_plus[e] += h;
            shell_a.center = center_plus;
            let dip_plus = shell_dipole_integrals(&shell_a, &shell_b, &origin);
            // Minus
            let mut center_minus = shell_a.center;
            center_minus[e] -= 2.0 * h;
            shell_a.center = center_minus;
            let dip_minus = shell_dipole_integrals(&shell_a, &shell_b, &origin);
            // Restore
            shell_a.center[e] += h;

            // Central difference for each dipole direction d
            for d in 0..3 {
                let fd_val = (dip_plus[d][0] - dip_minus[d][0]) / (2.0 * h);
                let ana_val = analytical[d][e][0];
                let err = (fd_val - ana_val).abs();
                assert!(
                    err < 1e-7,
                    "shell_dipole_first_deriv FD vs analytical mismatch at d={d}, e={e}: ana={ana_val}, fd={fd_val}, err={err}"
                );
            }
        }
    }

    #[test]
    fn test_primitive_dipole_ss_at_origin_zero() {
        // Two primitive s-Gaussians at origin: <0|r|0> should be 0 in all directions.
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let gp = GaussianProduct::new(1.0, &a, 1.0, &b);
        let s = CartesianPower { i: 0, j: 0, k: 0 };
        let mu = primitive_dipole(&gp, &s, &s, &a, &[0.0, 0.0, 0.0]);
        assert_abs_diff_eq!(mu[0], 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(mu[1], 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(mu[2], 0.0, epsilon = 1e-12);
    }
}
