//! First-derivative one-electron integrals
//!
//! Standalone shell-pair-level functions for computing first derivatives of
//! overlap, kinetic, and nuclear attraction integrals with respect to nuclear
//! (basis function center) positions.
//!
//! These functions use the angular momentum raising/lowering identity
//! (nabla identity) from Helgaker, Jorgensen & Olsen (2000), Eq. 9.3.32:
//!
//! ```text
//! d/dA_x [N * g(alpha, A, l)] = N * [2*alpha * g(alpha, A, l+1_x) - l_x * g(alpha, A, l-1_x)]
//! ```
//!
//! The derivative with respect to center A (bra side) produces two terms:
//! - Term 1: +2*alpha * INTEGRAL(raised_angular, ket_angular) [angular momentum raised by 1]
//! - Term 2: -l_dir * INTEGRAL(lowered_angular, ket_angular) [angular momentum lowered by 1]
//!
//! The ORIGINAL normalization factors are used (not the shifted normalization).
//!
//! # Design
//!
//! These functions are factored out from `scf::gradient` to enable reuse by:
//! - The Hessian code (H^1 construction for second derivatives)
//! - Finite-difference validation of second derivatives
//! - Standalone integral derivative inspection
//!
//! # References
//!
//! - Helgaker, Jorgensen & Olsen (2000), Ch. 9, Eq. 9.3.32
//! - Pulay (1969), Mol. Phys. 17, 197
//! - PySCF `grad/rhf.py` lines 33-76

use crate::basis::ContractedShell;
use crate::integrals::overlap::{cartesian_gaussian_normalization, overlap_1d};
use crate::integrals::{cartesian_components, CartesianPower, GaussianProduct};

// ============================================================================
// Angular Momentum Helpers (duplicated from gradient.rs for standalone use)
// ============================================================================

/// Get the angular momentum component for a specific direction (0=x, 1=y, 2=z)
#[inline]
fn get_angular_component(pow: &CartesianPower, dir: usize) -> u32 {
    match dir {
        0 => pow.i,
        1 => pow.j,
        2 => pow.k,
        _ => unreachable!(),
    }
}

/// Raise angular momentum by 1 in the specified direction
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

/// Lower angular momentum by 1 in the specified direction (caller must ensure l_dir > 0)
#[inline]
fn lower_angular(pow: &CartesianPower, dir: usize) -> CartesianPower {
    let mut p = *pow;
    match dir {
        0 => {
            debug_assert!(p.i > 0);
            p.i -= 1;
        }
        1 => {
            debug_assert!(p.j > 0);
            p.j -= 1;
        }
        2 => {
            debug_assert!(p.k > 0);
            p.k -= 1;
        }
        _ => unreachable!(),
    }
    p
}

// ============================================================================
// Primitive integral wrappers
// ============================================================================

/// Compute primitive overlap integral from pre-computed GaussianProduct.
/// Wrapper that calls the 1D overlap functions and combines.
///
/// Reference: gradient.rs `primitive_overlap_from_gp`
#[inline]
fn primitive_overlap_from_gp(
    gp: &GaussianProduct,
    a_powers: &CartesianPower,
    b_powers: &CartesianPower,
) -> f64 {
    let s_x = overlap_1d(
        gp.pa[0],
        gp.ab[0],
        gp.one_over_2p,
        a_powers.i as i32,
        b_powers.i as i32,
    );
    let s_y = overlap_1d(
        gp.pa[1],
        gp.ab[1],
        gp.one_over_2p,
        a_powers.j as i32,
        b_powers.j as i32,
    );
    let s_z = overlap_1d(
        gp.pa[2],
        gp.ab[2],
        gp.one_over_2p,
        a_powers.k as i32,
        b_powers.k as i32,
    );
    gp.ss_integral * s_x * s_y * s_z
}

/// Compute primitive kinetic integral from pre-computed GaussianProduct.
/// Wrapper that calls the 1D overlap/kinetic functions and combines.
///
/// Reference: gradient.rs `primitive_kinetic_from_gp`
#[inline]
fn primitive_kinetic_from_gp(
    gp: &GaussianProduct,
    a_powers: &CartesianPower,
    b_powers: &CartesianPower,
    beta: f64,
) -> f64 {
    use crate::integrals::kinetic_1d;

    let a_x = a_powers.i as i32;
    let a_y = a_powers.j as i32;
    let a_z = a_powers.k as i32;
    let b_x = b_powers.i as i32;
    let b_y = b_powers.j as i32;
    let b_z = b_powers.k as i32;

    let s_x = overlap_1d(gp.pa[0], gp.ab[0], gp.one_over_2p, a_x, b_x);
    let s_y = overlap_1d(gp.pa[1], gp.ab[1], gp.one_over_2p, a_y, b_y);
    let s_z = overlap_1d(gp.pa[2], gp.ab[2], gp.one_over_2p, a_z, b_z);

    let t_x = kinetic_1d(gp.pa[0], gp.ab[0], gp.one_over_2p, a_x, b_x, beta);
    let t_y = kinetic_1d(gp.pa[1], gp.ab[1], gp.one_over_2p, a_y, b_y, beta);
    let t_z = kinetic_1d(gp.pa[2], gp.ab[2], gp.one_over_2p, a_z, b_z, beta);

    gp.ss_integral * (t_x * s_y * s_z + s_x * t_y * s_z + s_x * s_y * t_z)
}

// ============================================================================
// Public API: First-derivative one-electron integrals
// ============================================================================

/// Compute the first derivative of the overlap integral with respect to the
/// bra center (shell_a).
///
/// Returns `[dS/dx, dS/dy, dS/dz]` where each element is a flattened
/// `n_a x n_b` matrix (row-major).
///
/// Uses the nabla identity (Helgaker et al., Eq. 9.3.32):
/// ```text
/// dS_{mu,nu}/dA_x = sum_{prim} c_a * c_b * N_a * N_b *
///     [2*alpha * <g_{a+1_x} | g_b> - l_{a_x} * <g_{a-1_x} | g_b>]
/// ```
///
/// The ORIGINAL normalization N_a is preserved (not the shifted normalization
/// of the raised/lowered Gaussian).
///
/// # Arguments
///
/// * `shell_a` - Bra shell (the center being differentiated)
/// * `shell_b` - Ket shell
///
/// # Returns
///
/// Array of three vectors `[dS/dx, dS/dy, dS/dz]`, each of length `n_a * n_b`.
pub fn shell_overlap_first_deriv(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
) -> [Vec<f64>; 3] {
    let l_a = shell_a.l_value();
    let l_b = shell_b.l_value();
    let comps_a = cartesian_components(l_a).expect("Angular momentum in range");
    let comps_b = cartesian_components(l_b).expect("Angular momentum in range");
    let n_a = comps_a.len();
    let n_b = comps_b.len();

    let mut result = [
        vec![0.0; n_a * n_b],
        vec![0.0; n_a * n_b],
        vec![0.0; n_a * n_b],
    ];

    for prim_a in &shell_a.primitives {
        for prim_b in &shell_b.primitives {
            let alpha = prim_a.exponent;
            let gp = GaussianProduct::new(alpha, &shell_a.center, prim_b.exponent, &shell_b.center);
            let coef = prim_a.coefficient * prim_b.coefficient;

            for (ia, pow_a) in comps_a.iter().enumerate() {
                let norm_a = cartesian_gaussian_normalization(alpha, pow_a);

                for (jb, pow_b) in comps_b.iter().enumerate() {
                    let norm_b = cartesian_gaussian_normalization(prim_b.exponent, pow_b);
                    let cn = coef * norm_a * norm_b;

                    for (dir, res_dir) in result.iter_mut().enumerate() {
                        let l_dir = get_angular_component(pow_a, dir);

                        // Term 1: +2*alpha * <g_{a+1_dir} | g_b>
                        let pow_a_plus = raise_angular(pow_a, dir);
                        let s_plus = primitive_overlap_from_gp(&gp, &pow_a_plus, pow_b);
                        let term1 = 2.0 * alpha * cn * s_plus;

                        // Term 2: -l_dir * <g_{a-1_dir} | g_b>
                        let term2 = if l_dir > 0 {
                            let pow_a_minus = lower_angular(pow_a, dir);
                            let s_minus = primitive_overlap_from_gp(&gp, &pow_a_minus, pow_b);
                            -(l_dir as f64) * cn * s_minus
                        } else {
                            0.0
                        };

                        res_dir[ia * n_b + jb] += term1 + term2;
                    }
                }
            }
        }
    }

    result
}

/// Compute the first derivative of the kinetic energy integral with respect
/// to the bra center (shell_a).
///
/// Returns `[dT/dx, dT/dy, dT/dz]` where each element is a flattened
/// `n_a x n_b` matrix (row-major).
///
/// Uses the same nabla identity as the overlap derivative, but applied to
/// the kinetic energy integral T = Tx*Sy*Sz + Sx*Ty*Sz + Sx*Sy*Tz.
///
/// # Arguments
///
/// * `shell_a` - Bra shell (the center being differentiated)
/// * `shell_b` - Ket shell
///
/// # Returns
///
/// Array of three vectors `[dT/dx, dT/dy, dT/dz]`, each of length `n_a * n_b`.
pub fn shell_kinetic_first_deriv(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
) -> [Vec<f64>; 3] {
    let comps_a = cartesian_components(shell_a.l_value()).expect("Angular momentum in range");
    let comps_b = cartesian_components(shell_b.l_value()).expect("Angular momentum in range");
    let n_a = comps_a.len();
    let n_b = comps_b.len();

    let mut result = [
        vec![0.0; n_a * n_b],
        vec![0.0; n_a * n_b],
        vec![0.0; n_a * n_b],
    ];

    for prim_a in &shell_a.primitives {
        for prim_b in &shell_b.primitives {
            let alpha = prim_a.exponent;
            let beta = prim_b.exponent;
            let gp = GaussianProduct::new(alpha, &shell_a.center, beta, &shell_b.center);
            let coef = prim_a.coefficient * prim_b.coefficient;

            for (ia, pow_a) in comps_a.iter().enumerate() {
                let norm_a = cartesian_gaussian_normalization(alpha, pow_a);

                for (jb, pow_b) in comps_b.iter().enumerate() {
                    let norm_b = cartesian_gaussian_normalization(beta, pow_b);
                    let cn = coef * norm_a * norm_b;

                    for (dir, res_dir) in result.iter_mut().enumerate() {
                        let l_dir = get_angular_component(pow_a, dir);

                        // Term 1: +2*alpha * <g_{a+1_dir} | T | g_b>
                        let pow_a_plus = raise_angular(pow_a, dir);
                        let t_plus = primitive_kinetic_from_gp(&gp, &pow_a_plus, pow_b, beta);
                        let term1 = 2.0 * alpha * cn * t_plus;

                        // Term 2: -l_dir * <g_{a-1_dir} | T | g_b>
                        let term2 = if l_dir > 0 {
                            let pow_a_minus = lower_angular(pow_a, dir);
                            let t_minus = primitive_kinetic_from_gp(&gp, &pow_a_minus, pow_b, beta);
                            -(l_dir as f64) * cn * t_minus
                        } else {
                            0.0
                        };

                        res_dir[ia * n_b + jb] += term1 + term2;
                    }
                }
            }
        }
    }

    result
}

/// Compute the first derivative of the nuclear attraction integral with
/// respect to the bra center (basis function derivative).
///
/// Returns `[dV/dx, dV/dy, dV/dz]` where each element is a flattened
/// `n_a x n_b` matrix (row-major).
///
/// This computes the derivative arising from differentiating the basis
/// function on center A (bra side). The nuclear charge Z is included
/// in the result.
///
/// # Arguments
///
/// * `shell_a` - Bra shell (the center being differentiated)
/// * `shell_b` - Ket shell
/// * `nuclear_center` - Position of the nucleus [x, y, z] in bohr
/// * `nuclear_charge` - Nuclear charge Z (as f64)
///
/// # Returns
///
/// Array of three vectors `[dV/dx, dV/dy, dV/dz]`, each of length `n_a * n_b`.
pub fn shell_nuclear_first_deriv_basis(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
    nuclear_center: &[f64; 3],
    nuclear_charge: f64,
) -> [Vec<f64>; 3] {
    let comps_a = cartesian_components(shell_a.l_value()).expect("Angular momentum in range");
    let comps_b = cartesian_components(shell_b.l_value()).expect("Angular momentum in range");
    let n_a = comps_a.len();
    let n_b = comps_b.len();

    let mut result = [
        vec![0.0; n_a * n_b],
        vec![0.0; n_a * n_b],
        vec![0.0; n_a * n_b],
    ];

    for prim_a in &shell_a.primitives {
        for prim_b in &shell_b.primitives {
            let alpha = prim_a.exponent;
            let gp = GaussianProduct::new(alpha, &shell_a.center, prim_b.exponent, &shell_b.center);
            let coef = prim_a.coefficient * prim_b.coefficient;

            for (ia, pow_a) in comps_a.iter().enumerate() {
                let norm_a = cartesian_gaussian_normalization(alpha, pow_a);

                for (jb, pow_b) in comps_b.iter().enumerate() {
                    let norm_b = cartesian_gaussian_normalization(prim_b.exponent, pow_b);
                    let cn = coef * norm_a * norm_b * nuclear_charge;

                    for (dir, res_dir) in result.iter_mut().enumerate() {
                        let l_dir = get_angular_component(pow_a, dir);

                        // Term 1: +2*alpha * <g_{a+1_dir} | V | g_b>
                        let pow_a_plus = raise_angular(pow_a, dir);
                        let v_plus = crate::integrals::primitive_nuclear(
                            &gp,
                            &pow_a_plus,
                            pow_b,
                            nuclear_center,
                        );
                        let term1 = 2.0 * alpha * cn * v_plus;

                        // Term 2: -l_dir * <g_{a-1_dir} | V | g_b>
                        let term2 = if l_dir > 0 {
                            let pow_a_minus = lower_angular(pow_a, dir);
                            let v_minus = crate::integrals::primitive_nuclear(
                                &gp,
                                &pow_a_minus,
                                pow_b,
                                nuclear_center,
                            );
                            -(l_dir as f64) * cn * v_minus
                        } else {
                            0.0
                        };

                        res_dir[ia * n_b + jb] += term1 + term2;
                    }
                }
            }
        }
    }

    result
}

/// Compute the first derivative of the nuclear attraction integral with
/// respect to the nuclear center position.
///
/// Returns `[dV/dx, dV/dy, dV/dz]` for the nuclear center C.
///
/// Uses translational invariance of the 3-center integral:
/// ```text
/// dV^C/dA + dV^C/dB + dV^C/dC = 0
/// ```
/// Therefore:
/// ```text
/// dV^C/dC = -(dV^C/dA + dV^C/dB)
/// ```
/// where `dV^C/dA` is the bra-center derivative (computed by
/// `shell_nuclear_first_deriv_basis(a, b)`) and `dV^C/dB` is the ket-center
/// derivative (obtained by computing `shell_nuclear_first_deriv_basis(b, a)`
/// and transposing the AO indices).
///
/// This is fully analytical (no finite differences) and matches the precision
/// of the bra-center derivative used in the gradient and Hessian codes.
///
/// # Arguments
///
/// * `shell_a` - Bra shell
/// * `shell_b` - Ket shell
/// * `nuclear_center` - Position of the nucleus [x, y, z] in bohr
/// * `nuclear_charge` - Nuclear charge Z (as f64)
///
/// # Returns
///
/// Array of three vectors `[dV/dx, dV/dy, dV/dz]`, each of length `n_a * n_b`.
///
/// # References
///
/// - Translational invariance: Helgaker, Jorgensen & Olsen (2000), Ch. 9
/// - PySCF `grad/rhf.py` — nuclear center derivative via translational invariance
pub fn shell_nuclear_first_deriv_center(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
    nuclear_center: &[f64; 3],
    nuclear_charge: f64,
) -> [Vec<f64>; 3] {
    let n_a = shell_a.n_basis_functions();
    let n_b = shell_b.n_basis_functions();

    // dV^C/dA: derivative of <A|V^C|B> wrt center A (bra)
    let dv_da = shell_nuclear_first_deriv_basis(shell_a, shell_b, nuclear_center, nuclear_charge);

    // dV^C/dB: derivative of <A|V^C|B> wrt center B (ket)
    // Compute via: <B|V^C|A> is the same integral with swapped bra/ket,
    // so shell_nuclear_first_deriv_basis(b, a, C, Z) gives d<B|V^C|A>/dB
    // which is (n_b x n_a). Transpose to get (n_a x n_b) for d<A|V^C|B>/dB.
    let dv_db_transposed =
        shell_nuclear_first_deriv_basis(shell_b, shell_a, nuclear_center, nuclear_charge);

    // dV^C/dC = -(dV^C/dA + dV^C/dB)
    let mut result = [
        vec![0.0; n_a * n_b],
        vec![0.0; n_a * n_b],
        vec![0.0; n_a * n_b],
    ];

    for dir in 0..3 {
        for ia in 0..n_a {
            for jb in 0..n_b {
                let idx_ab = ia * n_b + jb;
                let idx_ba = jb * n_a + ia; // transposed index
                result[dir][idx_ab] = -(dv_da[dir][idx_ab] + dv_db_transposed[dir][idx_ba]);
            }
        }
    }

    result
}

// ============================================================================
// Public API: Second-derivative one-electron integrals
// ============================================================================

/// Compute the second derivative of the overlap integral with respect to
/// both the bra center A (direction d) and ket center B (direction e).
///
/// Returns 9 components `[d²S/dA_x·dB_x, d²S/dA_x·dB_y, ..., d²S/dA_z·dB_z]`
/// where `result[d*3 + e]` = ∂²S/∂A_d ∂B_e.
///
/// Each component is a flattened `n_a x n_b` matrix (row-major).
///
/// # Algorithm
///
/// Uses the double nabla identity (Helgaker et al., Ch. 9):
/// ```text
/// ∂²[μν]/∂A_d ∂B_e = 4·α·β · [μ_d+1, ν_e+1]      (Term 1)
///                   - 2·α·ν_e · [μ_d+1, ν_e-1]      (Term 2, only if ν_e > 0)
///                   - 2·μ_d·β · [μ_d-1, ν_e+1]      (Term 3, only if μ_d > 0)
///                   + μ_d·ν_e · [μ_d-1, ν_e-1]      (Term 4, only if both > 0)
/// ```
///
/// Apply nabla on ket first (∂/∂B_e), then on bra (∂/∂A_d).
/// The ORIGINAL normalization N(α, pow_a) and N(β, pow_b) is preserved.
///
/// # References
///
/// - Helgaker, Jorgensen & Olsen (2000), Ch. 9, Eq. 9.3.32 (extended to 2nd order)
/// - Phase 5 plan: Section 4a, second-derivative nabla identity
/// - PySCF: `int1e_ipovlpip` (∂²/∂A_d ∂B_e pattern)
pub fn shell_overlap_second_deriv(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
) -> [Vec<f64>; 9] {
    let comps_a = cartesian_components(shell_a.l_value()).expect("Angular momentum in range");
    let comps_b = cartesian_components(shell_b.l_value()).expect("Angular momentum in range");
    let n_a = comps_a.len();
    let n_b = comps_b.len();

    let mut result: [Vec<f64>; 9] = std::array::from_fn(|_| vec![0.0; n_a * n_b]);

    for prim_a in &shell_a.primitives {
        for prim_b in &shell_b.primitives {
            let alpha = prim_a.exponent;
            let beta = prim_b.exponent;
            let gp = GaussianProduct::new(alpha, &shell_a.center, beta, &shell_b.center);
            let coef = prim_a.coefficient * prim_b.coefficient;

            for (ia, pow_a) in comps_a.iter().enumerate() {
                let norm_a = cartesian_gaussian_normalization(alpha, pow_a);

                for (jb, pow_b) in comps_b.iter().enumerate() {
                    let norm_b = cartesian_gaussian_normalization(beta, pow_b);
                    let cn = coef * norm_a * norm_b;

                    for d in 0..3 {
                        let mu_d = get_angular_component(pow_a, d);

                        for e in 0..3 {
                            let nu_e = get_angular_component(pow_b, e);

                            // Term 1: +4·α·β · <g_{a+1_d} | g_{b+1_e}>
                            let pow_a_up = raise_angular(pow_a, d);
                            let pow_b_up = raise_angular(pow_b, e);
                            let s_pp = primitive_overlap_from_gp(&gp, &pow_a_up, &pow_b_up);
                            let mut val = 4.0 * alpha * beta * cn * s_pp;

                            // Term 2: -2·α·ν_e · <g_{a+1_d} | g_{b-1_e}>
                            if nu_e > 0 {
                                let pow_b_down = lower_angular(pow_b, e);
                                let s_pm = primitive_overlap_from_gp(&gp, &pow_a_up, &pow_b_down);
                                val -= 2.0 * alpha * (nu_e as f64) * cn * s_pm;
                            }

                            // Term 3: -2·μ_d·β · <g_{a-1_d} | g_{b+1_e}>
                            if mu_d > 0 {
                                let pow_a_down = lower_angular(pow_a, d);
                                let s_mp = primitive_overlap_from_gp(&gp, &pow_a_down, &pow_b_up);
                                val -= 2.0 * (mu_d as f64) * beta * cn * s_mp;
                            }

                            // Term 4: +μ_d·ν_e · <g_{a-1_d} | g_{b-1_e}>
                            if mu_d > 0 && nu_e > 0 {
                                let pow_a_down = lower_angular(pow_a, d);
                                let pow_b_down = lower_angular(pow_b, e);
                                let s_mm = primitive_overlap_from_gp(&gp, &pow_a_down, &pow_b_down);
                                val += (mu_d as f64) * (nu_e as f64) * cn * s_mm;
                            }

                            result[d * 3 + e][ia * n_b + jb] += val;
                        }
                    }
                }
            }
        }
    }

    result
}

/// Compute the second derivative of the kinetic energy integral with respect to
/// both the bra center A (direction d) and ket center B (direction e).
///
/// Returns 9 components `[d²T/dA_x·dB_x, d²T/dA_x·dB_y, ..., d²T/dA_z·dB_z]`
/// where `result[d*3 + e]` = ∂²T/∂A_d ∂B_e.
///
/// Each component is a flattened `n_a x n_b` matrix (row-major).
///
/// # Algorithm
///
/// Same double nabla identity as the overlap second derivative, but applied to
/// the kinetic energy integral T = ⟨μ|-½∇²|ν⟩. The nabla raises/lowers the
/// angular momentum of the *basis functions*, not the kinetic operator:
///
/// ```text
/// ∂²T_{μν}/∂A_d ∂B_e = 4·α·β · T(μ_d+1, ν_e+1)
///                     - 2·α·ν_e · T(μ_d+1, ν_e-1)
///                     - 2·μ_d·β · T(μ_d-1, ν_e+1)
///                     + μ_d·ν_e · T(μ_d-1, ν_e-1)
/// ```
///
/// # References
///
/// - Helgaker, Jorgensen & Olsen (2000), Ch. 9
/// - PySCF: `int1e_ipkinip` (∂²/∂A_d ∂B_e pattern)
pub fn shell_kinetic_second_deriv(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
) -> [Vec<f64>; 9] {
    let comps_a = cartesian_components(shell_a.l_value()).expect("Angular momentum in range");
    let comps_b = cartesian_components(shell_b.l_value()).expect("Angular momentum in range");
    let n_a = comps_a.len();
    let n_b = comps_b.len();

    let mut result: [Vec<f64>; 9] = std::array::from_fn(|_| vec![0.0; n_a * n_b]);

    for prim_a in &shell_a.primitives {
        for prim_b in &shell_b.primitives {
            let alpha = prim_a.exponent;
            let beta = prim_b.exponent;
            let gp = GaussianProduct::new(alpha, &shell_a.center, beta, &shell_b.center);
            let coef = prim_a.coefficient * prim_b.coefficient;

            for (ia, pow_a) in comps_a.iter().enumerate() {
                let norm_a = cartesian_gaussian_normalization(alpha, pow_a);

                for (jb, pow_b) in comps_b.iter().enumerate() {
                    let norm_b = cartesian_gaussian_normalization(beta, pow_b);
                    let cn = coef * norm_a * norm_b;

                    for d in 0..3 {
                        let mu_d = get_angular_component(pow_a, d);

                        for e in 0..3 {
                            let nu_e = get_angular_component(pow_b, e);

                            // Term 1: +4·α·β · T(a+1_d, b+1_e)
                            let pow_a_up = raise_angular(pow_a, d);
                            let pow_b_up = raise_angular(pow_b, e);
                            let t_pp = primitive_kinetic_from_gp(&gp, &pow_a_up, &pow_b_up, beta);
                            let mut val = 4.0 * alpha * beta * cn * t_pp;

                            // Term 2: -2·α·ν_e · T(a+1_d, b-1_e)
                            if nu_e > 0 {
                                let pow_b_down = lower_angular(pow_b, e);
                                let t_pm =
                                    primitive_kinetic_from_gp(&gp, &pow_a_up, &pow_b_down, beta);
                                val -= 2.0 * alpha * (nu_e as f64) * cn * t_pm;
                            }

                            // Term 3: -2·μ_d·β · T(a-1_d, b+1_e)
                            if mu_d > 0 {
                                let pow_a_down = lower_angular(pow_a, d);
                                let t_mp =
                                    primitive_kinetic_from_gp(&gp, &pow_a_down, &pow_b_up, beta);
                                val -= 2.0 * (mu_d as f64) * beta * cn * t_mp;
                            }

                            // Term 4: +μ_d·ν_e · T(a-1_d, b-1_e)
                            if mu_d > 0 && nu_e > 0 {
                                let pow_a_down = lower_angular(pow_a, d);
                                let pow_b_down = lower_angular(pow_b, e);
                                let t_mm =
                                    primitive_kinetic_from_gp(&gp, &pow_a_down, &pow_b_down, beta);
                                val += (mu_d as f64) * (nu_e as f64) * cn * t_mm;
                            }

                            result[d * 3 + e][ia * n_b + jb] += val;
                        }
                    }
                }
            }
        }
    }

    result
}

/// Compute the second derivative of the nuclear attraction integral with
/// respect to both the bra center A (direction d) and ket center B (direction e).
///
/// Returns 9 components `[d²V/dA_x·dB_x, d²V/dA_x·dB_y, ..., d²V/dA_z·dB_z]`
/// where `result[d*3 + e]` = ∂²V/∂A_d ∂B_e.
///
/// Each component is a flattened `n_a x n_b` matrix (row-major).
///
/// The nuclear charge Z is included in the result.
///
/// # Algorithm
///
/// Same double nabla identity as the overlap second derivative, but applied to
/// the nuclear attraction integral V = ⟨μ|-Z/|r-C||ν⟩:
///
/// ```text
/// ∂²V_{μν}/∂A_d ∂B_e = 4·α·β · V(μ_d+1, ν_e+1)
///                     - 2·α·ν_e · V(μ_d+1, ν_e-1)
///                     - 2·μ_d·β · V(μ_d-1, ν_e+1)
///                     + μ_d·ν_e · V(μ_d-1, ν_e-1)
/// ```
///
/// # Arguments
///
/// * `shell_a` - Bra shell (center A being differentiated)
/// * `shell_b` - Ket shell (center B being differentiated)
/// * `nuclear_center` - Position of the nucleus [x, y, z] in bohr
/// * `nuclear_charge` - Nuclear charge Z (as f64)
///
/// # References
///
/// - Helgaker, Jorgensen & Olsen (2000), Ch. 9
/// - PySCF: `int1e_ipnucip` (∂²/∂A_d ∂B_e pattern)
pub fn shell_nuclear_second_deriv(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
    nuclear_center: &[f64; 3],
    nuclear_charge: f64,
) -> [Vec<f64>; 9] {
    let comps_a = cartesian_components(shell_a.l_value()).expect("Angular momentum in range");
    let comps_b = cartesian_components(shell_b.l_value()).expect("Angular momentum in range");
    let n_a = comps_a.len();
    let n_b = comps_b.len();

    let mut result: [Vec<f64>; 9] = std::array::from_fn(|_| vec![0.0; n_a * n_b]);

    for prim_a in &shell_a.primitives {
        for prim_b in &shell_b.primitives {
            let alpha = prim_a.exponent;
            let beta = prim_b.exponent;
            let gp = GaussianProduct::new(alpha, &shell_a.center, beta, &shell_b.center);
            let coef = prim_a.coefficient * prim_b.coefficient;

            for (ia, pow_a) in comps_a.iter().enumerate() {
                let norm_a = cartesian_gaussian_normalization(alpha, pow_a);

                for (jb, pow_b) in comps_b.iter().enumerate() {
                    let norm_b = cartesian_gaussian_normalization(beta, pow_b);
                    let cn = coef * norm_a * norm_b * nuclear_charge;

                    for d in 0..3 {
                        let mu_d = get_angular_component(pow_a, d);

                        for e in 0..3 {
                            let nu_e = get_angular_component(pow_b, e);

                            // Term 1: +4·α·β · V(a+1_d, b+1_e)
                            let pow_a_up = raise_angular(pow_a, d);
                            let pow_b_up = raise_angular(pow_b, e);
                            let v_pp = crate::integrals::primitive_nuclear(
                                &gp,
                                &pow_a_up,
                                &pow_b_up,
                                nuclear_center,
                            );
                            let mut val = 4.0 * alpha * beta * cn * v_pp;

                            // Term 2: -2·α·ν_e · V(a+1_d, b-1_e)
                            if nu_e > 0 {
                                let pow_b_down = lower_angular(pow_b, e);
                                let v_pm = crate::integrals::primitive_nuclear(
                                    &gp,
                                    &pow_a_up,
                                    &pow_b_down,
                                    nuclear_center,
                                );
                                val -= 2.0 * alpha * (nu_e as f64) * cn * v_pm;
                            }

                            // Term 3: -2·μ_d·β · V(a-1_d, b+1_e)
                            if mu_d > 0 {
                                let pow_a_down = lower_angular(pow_a, d);
                                let v_mp = crate::integrals::primitive_nuclear(
                                    &gp,
                                    &pow_a_down,
                                    &pow_b_up,
                                    nuclear_center,
                                );
                                val -= 2.0 * (mu_d as f64) * beta * cn * v_mp;
                            }

                            // Term 4: +μ_d·ν_e · V(a-1_d, b-1_e)
                            if mu_d > 0 && nu_e > 0 {
                                let pow_a_down = lower_angular(pow_a, d);
                                let pow_b_down = lower_angular(pow_b, e);
                                let v_mm = crate::integrals::primitive_nuclear(
                                    &gp,
                                    &pow_a_down,
                                    &pow_b_down,
                                    nuclear_center,
                                );
                                val += (mu_d as f64) * (nu_e as f64) * cn * v_mm;
                            }

                            result[d * 3 + e][ia * n_b + jb] += val;
                        }
                    }
                }
            }
        }
    }

    result
}

/// Compute the same-center (bra-bra) second derivative of the nuclear attraction
/// integral: ∂²V^C/∂A_d ∂A_e.
///
/// Both derivatives act on the BRA center A. This is needed for nuclear-center
/// contributions to the analytical Hessian, where translational invariance
/// relates the nuclear-center derivative to basis-center derivatives:
///
/// ```text
/// ∂²V^C/∂R_C ∂R_B = -(∂²V^C/∂R_A ∂R_B + ∂²V^C/∂R_B ∂R_B)
/// ```
///
/// # Algorithm
///
/// Applies the nabla identity twice on the bra center:
///
/// ```text
/// ∂²V^C/∂A_d ∂A_e = 4·α² · V(a+1_d+1_e, b)
///                   - 2·α·(μ_e + δ_{de}) · V(a+1_d-1_e, b)
///                   - 2·α·μ_d · V(a-1_d+1_e, b)
///                   + μ_d·(μ_e - δ_{de}) · V(a-1_d-1_e, b)
/// ```
///
/// Note the δ_{de} correction: when d=e, the angular momentum raised/lowered in
/// direction d also modifies the angular momentum that the second derivative in
/// direction e=d sees.
///
/// # Arguments
///
/// * `shell_a` - Bra shell (both derivatives act on this center)
/// * `shell_b` - Ket shell (unchanged)
/// * `nuclear_center` - Position of nucleus C
/// * `nuclear_charge` - Charge Z of nucleus C
///
/// # Returns
///
/// 9-component array `[∂²V/∂A_x·∂A_x, ∂²V/∂A_x·∂A_y, ..., ∂²V/∂A_z·∂A_z]`
/// indexed as `[d*3 + e]`.
///
/// # References
///
/// - Helgaker, Jorgensen & Olsen (2000), Ch. 9 (double nabla on same center)
/// - PySCF: `int1e_ipipnuc` (∂²/∂A_d ∂A_e pattern)
pub fn shell_nuclear_second_deriv_bra_bra(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
    nuclear_center: &[f64; 3],
    nuclear_charge: f64,
) -> [Vec<f64>; 9] {
    let comps_a = cartesian_components(shell_a.l_value()).expect("Angular momentum in range");
    let comps_b = cartesian_components(shell_b.l_value()).expect("Angular momentum in range");
    let n_a = comps_a.len();
    let n_b = comps_b.len();

    let mut result: [Vec<f64>; 9] = std::array::from_fn(|_| vec![0.0; n_a * n_b]);

    for prim_a in &shell_a.primitives {
        for prim_b in &shell_b.primitives {
            let alpha = prim_a.exponent;
            let beta = prim_b.exponent;
            let gp = GaussianProduct::new(alpha, &shell_a.center, beta, &shell_b.center);
            let coef = prim_a.coefficient * prim_b.coefficient;

            for (ia, pow_a) in comps_a.iter().enumerate() {
                let norm_a = cartesian_gaussian_normalization(alpha, pow_a);

                for (jb, pow_b) in comps_b.iter().enumerate() {
                    let norm_b = cartesian_gaussian_normalization(beta, pow_b);
                    let cn = coef * norm_a * norm_b * nuclear_charge;

                    for d in 0..3 {
                        let mu_d = get_angular_component(pow_a, d);

                        for e in 0..3 {
                            // After raising in direction d, the angular momentum in direction e is:
                            // mu_e_after_raise_d = mu_e + delta(d,e)
                            let mu_e_orig = get_angular_component(pow_a, e);
                            let delta_de: u32 = if d == e { 1 } else { 0 };

                            // Term 1: +4·α² · V(a+1_d+1_e, b)
                            let pow_a_up_d = raise_angular(pow_a, d);
                            let pow_a_up_d_up_e = raise_angular(&pow_a_up_d, e);
                            let v_pp = crate::integrals::primitive_nuclear(
                                &gp,
                                &pow_a_up_d_up_e,
                                pow_b,
                                nuclear_center,
                            );
                            let mut val = 4.0 * alpha * alpha * cn * v_pp;

                            // Term 2: -2·α·(μ_e + δ_{de}) · V(a+1_d-1_e, b)
                            // Only if μ_e + δ_{de} > 0 (which is always true since δ_{de}=1 when d=e)
                            let mu_e_raised = mu_e_orig + delta_de;
                            if mu_e_raised > 0 {
                                let pow_a_up_d_down_e = lower_angular(&pow_a_up_d, e);
                                let v_pm = crate::integrals::primitive_nuclear(
                                    &gp,
                                    &pow_a_up_d_down_e,
                                    pow_b,
                                    nuclear_center,
                                );
                                val -= 2.0 * alpha * (mu_e_raised as f64) * cn * v_pm;
                            }

                            // Term 3: -2·α·μ_d · V(a-1_d+1_e, b)
                            if mu_d > 0 {
                                let pow_a_down_d = lower_angular(pow_a, d);
                                let pow_a_down_d_up_e = raise_angular(&pow_a_down_d, e);
                                let v_mp = crate::integrals::primitive_nuclear(
                                    &gp,
                                    &pow_a_down_d_up_e,
                                    pow_b,
                                    nuclear_center,
                                );
                                val -= 2.0 * alpha * (mu_d as f64) * cn * v_mp;
                            }

                            // Term 4: +μ_d · (μ_e_after_lower_d) · V(a-1_d-1_e, b)
                            // After lowering d, the angular momentum in direction e is:
                            //   if d == e: μ_d - 1  (since lowering d also lowers e=d)
                            //   if d != e: μ_e      (unchanged)
                            // Need μ_d > 0 to lower d, then need the angular momentum
                            // in direction e of the lowered function to be > 0.
                            if mu_d > 0 {
                                let pow_a_down_d = lower_angular(pow_a, d);
                                let ang_e_after_lower_d = get_angular_component(&pow_a_down_d, e);
                                if ang_e_after_lower_d > 0 {
                                    let pow_a_down_d_down_e = lower_angular(&pow_a_down_d, e);
                                    let v_mm = crate::integrals::primitive_nuclear(
                                        &gp,
                                        &pow_a_down_d_down_e,
                                        pow_b,
                                        nuclear_center,
                                    );
                                    val += (mu_d as f64) * (ang_e_after_lower_d as f64) * cn * v_mm;
                                }
                            }

                            result[d * 3 + e][ia * n_b + jb] += val;
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
mod tests {
    use super::*;
    use crate::basis::{Atom, BasisSet};
    use crate::integrals;
    use approx::assert_abs_diff_eq;

    /// Build H2/STO-3G basis set
    fn h2_sto3g() -> BasisSet {
        let atoms = vec![
            Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
            Atom::new(1, [0.0, 0.0, 1.4]).unwrap(),
        ];
        BasisSet::build(atoms, "sto-3g").unwrap()
    }

    /// Build HF/STO-3G basis set
    fn hf_sto3g() -> BasisSet {
        let atoms = vec![
            Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
            Atom::new(9, [0.0, 0.0, 1.7328]).unwrap(), // HF bond length ~0.917 A
        ];
        BasisSet::build(atoms, "sto-3g").unwrap()
    }

    /// Build H2O/STO-3G basis set
    fn h2o_sto3g() -> BasisSet {
        let atoms = vec![
            Atom::new(8, [0.0, 0.0, 0.2217282]).unwrap(),
            Atom::new(1, [0.0, 1.4305447, -0.8869128]).unwrap(),
            Atom::new(1, [0.0, -1.4305447, -0.8869128]).unwrap(),
        ];
        BasisSet::build(atoms, "sto-3g").unwrap()
    }

    // =========================================================================
    // Test 1: Overlap first derivative (s|s) pair - H2/STO-3G
    // =========================================================================

    #[test]
    fn test_overlap_first_deriv_ss() {
        let basis = h2_sto3g();

        // Shell 0: 1s on H(0,0,0), Shell 1: 1s on H(0,0,1.4)
        let shell_a = &basis.shells[0];
        let shell_b = &basis.shells[1];

        let ds = shell_overlap_first_deriv(shell_a, shell_b);

        // For (s|s) pair: each direction gives a 1x1 matrix
        assert_eq!(ds[0].len(), 1);
        assert_eq!(ds[1].len(), 1);
        assert_eq!(ds[2].len(), 1);

        // x and y derivatives should be zero (bond along z)
        assert_abs_diff_eq!(ds[0][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(ds[1][0], 0.0, epsilon = 1e-14);

        // z derivative should be non-zero (atoms separated along z)
        assert!(
            ds[2][0].abs() > 1e-5,
            "dS/dz for (s|s) should be non-zero, got {}",
            ds[2][0]
        );

        // Translational invariance: dS/dA + dS/dB = 0
        // dS/dB is obtained by swapping shell_a and shell_b
        let ds_swap = shell_overlap_first_deriv(shell_b, shell_a);
        for dir in 0..3 {
            let sum = ds[dir][0] + ds_swap[dir][0];
            assert_abs_diff_eq!(sum, 0.0, epsilon = 1e-12,);
        }
    }

    // =========================================================================
    // Test 2: Overlap first derivative (s|p) pair - HF/STO-3G
    // =========================================================================

    #[test]
    fn test_overlap_first_deriv_sp() {
        let basis = hf_sto3g();

        // Find s shell on H and p shell on F
        let shell_s = &basis.shells[0]; // H 1s
                                        // F has shells: 1s(idx 1), 2s(idx 2), 2p(idx 3)
        let shell_p = &basis.shells[3]; // F 2p

        assert_eq!(shell_s.n_basis_functions(), 1);
        assert_eq!(shell_p.n_basis_functions(), 3);

        let ds = shell_overlap_first_deriv(shell_s, shell_p);

        // (s|p) pair: each direction gives a 1x3 matrix
        for dir in 0..3 {
            assert_eq!(ds[dir].len(), 3, "dir {} should give 1x3 = 3 elements", dir);
        }
    }

    // =========================================================================
    // Test 3: Kinetic first derivative (s|s) pair - H2/STO-3G
    // =========================================================================

    #[test]
    fn test_kinetic_first_deriv_ss() {
        let basis = h2_sto3g();

        let shell_a = &basis.shells[0];
        let shell_b = &basis.shells[1];

        let dt = shell_kinetic_first_deriv(shell_a, shell_b);

        // x and y derivatives should be zero
        assert_abs_diff_eq!(dt[0][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(dt[1][0], 0.0, epsilon = 1e-14);

        // z derivative should be non-zero
        assert!(
            dt[2][0].abs() > 1e-5,
            "dT/dz for (s|s) should be non-zero, got {}",
            dt[2][0]
        );

        // Translational invariance: dT/dA + dT/dB = 0
        let dt_swap = shell_kinetic_first_deriv(shell_b, shell_a);
        for dir in 0..3 {
            let sum = dt[dir][0] + dt_swap[dir][0];
            assert_abs_diff_eq!(sum, 0.0, epsilon = 1e-12,);
        }
    }

    // =========================================================================
    // Test 4: Nuclear first derivative (s|s) pair - H2/STO-3G
    // =========================================================================

    #[test]
    fn test_nuclear_first_deriv_ss() {
        let basis = h2_sto3g();

        let shell_a = &basis.shells[0]; // 1s on H(0,0,0)
        let shell_b = &basis.shells[1]; // 1s on H(0,0,1.4)
        let c = &basis.atoms[0].position; // Nucleus at H1
        let z = basis.atoms[0].atomic_number as f64;

        // Basis function derivative
        let dv_basis = shell_nuclear_first_deriv_basis(shell_a, shell_b, c, z);
        assert_eq!(dv_basis[0].len(), 1);
        assert_eq!(dv_basis[1].len(), 1);
        assert_eq!(dv_basis[2].len(), 1);

        // Nuclear center derivative
        let dv_center = shell_nuclear_first_deriv_center(shell_a, shell_b, c, z);
        assert_eq!(dv_center[0].len(), 1);
        assert_eq!(dv_center[1].len(), 1);
        assert_eq!(dv_center[2].len(), 1);

        // x and y basis derivatives should be zero (axial symmetry)
        assert_abs_diff_eq!(dv_basis[0][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(dv_basis[1][0], 0.0, epsilon = 1e-14);

        // z basis derivative should be non-zero
        assert!(
            dv_basis[2][0].abs() > 1e-6,
            "dV_basis/dz should be non-zero, got {}",
            dv_basis[2][0]
        );

        // Nuclear center derivative: x and y should be zero
        assert_abs_diff_eq!(dv_center[0][0], 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(dv_center[1][0], 0.0, epsilon = 1e-10);

        // z center derivative should be non-zero
        assert!(
            dv_center[2][0].abs() > 1e-6,
            "dV_center/dz should be non-zero, got {}",
            dv_center[2][0]
        );
    }

    // =========================================================================
    // Test 5: Overlap finite-difference validation
    // =========================================================================

    #[test]
    fn test_overlap_finite_diff() {
        // Verify dS/dA matches [S(A+h) - S(A-h)]/(2h) with h=1e-5
        let basis = h2_sto3g();
        let step = 1e-5;

        let shell_a = &basis.shells[0]; // 1s on H(0,0,0)
        let shell_b = &basis.shells[1]; // 1s on H(0,0,1.4)

        let ds_analytic = shell_overlap_first_deriv(shell_a, shell_b);

        for dir in 0..3 {
            // Displace shell_a center
            let mut center_plus = shell_a.center;
            center_plus[dir] += step;
            let mut center_minus = shell_a.center;
            center_minus[dir] -= step;

            let shell_a_plus = ContractedShell::new(
                shell_a.angular_momentum,
                shell_a.primitives.clone(),
                center_plus,
                shell_a.atom_idx,
            );
            let shell_a_minus = ContractedShell::new(
                shell_a.angular_momentum,
                shell_a.primitives.clone(),
                center_minus,
                shell_a.atom_idx,
            );

            // Compute overlap at displaced geometries
            let s_plus = compute_shell_overlap(&shell_a_plus, shell_b);
            let s_minus = compute_shell_overlap(&shell_a_minus, shell_b);

            for idx in 0..ds_analytic[dir].len() {
                let fd = (s_plus[idx] - s_minus[idx]) / (2.0 * step);
                assert_abs_diff_eq!(ds_analytic[dir][idx], fd, epsilon = 1e-8,);
            }
        }
    }

    /// Helper: compute overlap integral for a shell pair (contracted)
    fn compute_shell_overlap(shell_a: &ContractedShell, shell_b: &ContractedShell) -> Vec<f64> {
        let comps_a = cartesian_components(shell_a.l_value()).unwrap();
        let comps_b = cartesian_components(shell_b.l_value()).unwrap();
        let n_a = comps_a.len();
        let n_b = comps_b.len();
        let mut result = vec![0.0; n_a * n_b];

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
                        let s = primitive_overlap_from_gp(&gp, pow_a, pow_b);
                        result[ia * n_b + jb] += coef * norm_a * norm_b * s;
                    }
                }
            }
        }

        result
    }

    // =========================================================================
    // Test 6: Gradient equivalence - H2O/STO-3G
    //
    // Build full dS and dHcore matrices for H2O/STO-3G using deriv1 functions,
    // compare against finite-difference of the full matrices. Then contract
    // with a converged SCF density to verify the one-electron gradient
    // contribution matches.
    // =========================================================================

    /// Assemble the full dS/dR_A derivative matrix for a given atom and direction
    /// using the shell-pair deriv1 functions.
    fn assemble_ds_matrix(basis: &BasisSet, atom_idx: usize, dir: usize) -> Vec<f64> {
        let nbf = basis.n_basis;
        let mut ds = vec![0.0; nbf * nbf];

        let mut mu_offset = 0;
        for shell_a in &basis.shells {
            let n_a = shell_a.n_basis_functions();
            let mut nu_offset = 0;
            for shell_b in &basis.shells {
                let n_b = shell_b.n_basis_functions();

                // Bra side: differentiate shell_a
                if shell_a.atom_idx == atom_idx {
                    let ds_block = shell_overlap_first_deriv(shell_a, shell_b);
                    for ia in 0..n_a {
                        for ib in 0..n_b {
                            ds[(mu_offset + ia) * nbf + (nu_offset + ib)] +=
                                ds_block[dir][ia * n_b + ib];
                        }
                    }
                }

                // Ket side: differentiate shell_b (when it sits on atom_idx)
                // d<mu|nu>/dB_x where B is the center of shell_b
                // = d<nu|mu>/dB_x (overlap is symmetric under bra/ket swap)
                // = shell_overlap_first_deriv(shell_b, shell_a)[dir][ib * n_a + ia]
                if shell_b.atom_idx == atom_idx {
                    let ds_ket = shell_overlap_first_deriv(shell_b, shell_a);
                    for ia in 0..n_a {
                        for ib in 0..n_b {
                            ds[(mu_offset + ia) * nbf + (nu_offset + ib)] +=
                                ds_ket[dir][ib * n_a + ia];
                        }
                    }
                }

                nu_offset += n_b;
            }
            mu_offset += n_a;
        }

        ds
    }

    /// Assemble the full dH/dR_A derivative matrix (kinetic + nuclear)
    /// for a given atom and direction using deriv1 functions.
    fn assemble_dh_matrix(basis: &BasisSet, atom_idx: usize, dir: usize) -> Vec<f64> {
        let nbf = basis.n_basis;
        let mut dh = vec![0.0; nbf * nbf];

        let mut mu_offset = 0;
        for shell_a in &basis.shells {
            let n_a = shell_a.n_basis_functions();
            let mut nu_offset = 0;
            for shell_b in &basis.shells {
                let n_b = shell_b.n_basis_functions();

                // --- Kinetic derivative ---
                if shell_a.atom_idx == atom_idx {
                    let dt = shell_kinetic_first_deriv(shell_a, shell_b);
                    for ia in 0..n_a {
                        for ib in 0..n_b {
                            dh[(mu_offset + ia) * nbf + (nu_offset + ib)] += dt[dir][ia * n_b + ib];
                        }
                    }
                }
                if shell_b.atom_idx == atom_idx {
                    let dt_ket = shell_kinetic_first_deriv(shell_b, shell_a);
                    for ia in 0..n_a {
                        for ib in 0..n_b {
                            dh[(mu_offset + ia) * nbf + (nu_offset + ib)] +=
                                dt_ket[dir][ib * n_a + ia];
                        }
                    }
                }

                // --- Nuclear attraction derivative: basis function side ---
                if shell_a.atom_idx == atom_idx {
                    for atom_c in &basis.atoms {
                        let dv = shell_nuclear_first_deriv_basis(
                            shell_a,
                            shell_b,
                            &atom_c.position,
                            atom_c.atomic_number as f64,
                        );
                        for ia in 0..n_a {
                            for ib in 0..n_b {
                                dh[(mu_offset + ia) * nbf + (nu_offset + ib)] +=
                                    dv[dir][ia * n_b + ib];
                            }
                        }
                    }
                }
                if shell_b.atom_idx == atom_idx {
                    for atom_c in &basis.atoms {
                        let dv_ket = shell_nuclear_first_deriv_basis(
                            shell_b,
                            shell_a,
                            &atom_c.position,
                            atom_c.atomic_number as f64,
                        );
                        for ia in 0..n_a {
                            for ib in 0..n_b {
                                dh[(mu_offset + ia) * nbf + (nu_offset + ib)] +=
                                    dv_ket[dir][ib * n_a + ia];
                            }
                        }
                    }
                }

                nu_offset += n_b;
            }
            mu_offset += n_a;
        }

        // --- Nuclear center derivative ---
        let z_c = basis.atoms[atom_idx].atomic_number;
        if z_c > 0 {
            let mut mu_offset = 0;
            for shell_a in &basis.shells {
                let n_a = shell_a.n_basis_functions();
                let mut nu_offset = 0;
                for shell_b in &basis.shells {
                    let n_b = shell_b.n_basis_functions();

                    let dv_center = shell_nuclear_first_deriv_center(
                        shell_a,
                        shell_b,
                        &basis.atoms[atom_idx].position,
                        z_c as f64,
                    );
                    for ia in 0..n_a {
                        for ib in 0..n_b {
                            dh[(mu_offset + ia) * nbf + (nu_offset + ib)] +=
                                dv_center[dir][ia * n_b + ib];
                        }
                    }

                    nu_offset += n_b;
                }
                mu_offset += n_a;
            }
        }

        dh
    }

    #[test]
    fn test_gradient_equivalence_h2o() {
        // Validate the full dS and dH matrices assembled from deriv1 functions
        // against central finite-difference of the full S and H matrices.
        //
        // This is the definitive test: if the assembled derivative matrices match
        // finite-difference, then the deriv1 functions are correct AND the
        // assembly logic (bra/ket contributions, nuclear center term) is correct.

        let basis = h2o_sto3g();
        let nbf = basis.n_basis;
        let n_atoms = basis.atoms.len();
        let step = 1e-5;

        for atom_idx in 0..n_atoms {
            for dir in 0..3 {
                // Analytical derivative matrices from deriv1
                let ds_an = assemble_ds_matrix(&basis, atom_idx, dir);
                let dh_an = assemble_dh_matrix(&basis, atom_idx, dir);

                // Finite-difference derivative matrices
                let mut atoms_plus = basis.atoms.clone();
                atoms_plus[atom_idx].position[dir] += step;
                let bp = BasisSet::build(atoms_plus, "sto-3g").unwrap();
                let s_plus = integrals::overlap_matrix(&bp);
                let h_plus = integrals::hcore_matrix(&bp);

                let mut atoms_minus = basis.atoms.clone();
                atoms_minus[atom_idx].position[dir] -= step;
                let bm = BasisSet::build(atoms_minus, "sto-3g").unwrap();
                let s_minus = integrals::overlap_matrix(&bm);
                let h_minus = integrals::hcore_matrix(&bm);

                for mu in 0..nbf {
                    for nu in 0..nbf {
                        let ds_fd = (s_plus[mu * nbf + nu] - s_minus[mu * nbf + nu]) / (2.0 * step);
                        let dh_fd = (h_plus[mu * nbf + nu] - h_minus[mu * nbf + nu]) / (2.0 * step);

                        assert!(
                            (ds_an[mu * nbf + nu] - ds_fd).abs() < 1e-6,
                            "dS[{},{}] atom {} dir {}: analytic {:.6e} vs FD {:.6e} (err {:.2e})",
                            mu,
                            nu,
                            atom_idx,
                            dir,
                            ds_an[mu * nbf + nu],
                            ds_fd,
                            (ds_an[mu * nbf + nu] - ds_fd).abs(),
                        );
                        assert!(
                            (dh_an[mu * nbf + nu] - dh_fd).abs() < 1e-6,
                            "dH[{},{}] atom {} dir {}: analytic {:.6e} vs FD {:.6e} (err {:.2e})",
                            mu,
                            nu,
                            atom_idx,
                            dir,
                            dh_an[mu * nbf + nu],
                            dh_fd,
                            (dh_an[mu * nbf + nu] - dh_fd).abs(),
                        );
                    }
                }
            }
        }
    }

    // =========================================================================
    // Test 7: Overlap second derivative (s|s) pair - H2/STO-3G
    // =========================================================================

    #[test]
    fn test_overlap_second_deriv_ss() {
        let basis = h2_sto3g();

        let shell_a = &basis.shells[0]; // 1s on H(0,0,0)
        let shell_b = &basis.shells[1]; // 1s on H(0,0,1.4)

        let d2s = shell_overlap_second_deriv(shell_a, shell_b);

        // For (s|s) pair: each component gives a 1x1 matrix
        for comp in 0..9 {
            assert_eq!(d2s[comp].len(), 1, "Component {} should be 1x1", comp);
        }

        // d²S/dA_z·dB_z should be nonzero (bond along z axis)
        assert!(
            d2s[8][0].abs() > 1e-5,
            "d²S/dA_z·dB_z for (s|s) should be nonzero, got {}",
            d2s[8][0]
        );

        // Off-diagonal symmetry: d²S/dA_x·dB_y = d²S/dA_y·dB_x
        assert_abs_diff_eq!(d2s[0 * 3 + 1][0], d2s[1 * 3 + 0][0], epsilon = 1e-14);

        // Components involving x or y only should be zero (axial symmetry along z)
        // d²S/dA_x·dB_x: x-x should be nonzero (Gaussian width effect)
        // d²S/dA_x·dB_y should be zero
        assert_abs_diff_eq!(d2s[0 * 3 + 1][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(d2s[0 * 3 + 2][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(d2s[1 * 3 + 0][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(d2s[1 * 3 + 2][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(d2s[2 * 3 + 0][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(d2s[2 * 3 + 1][0], 0.0, epsilon = 1e-14);

        // d²S/dA_x·dB_x = d²S/dA_y·dB_y (cylindrical symmetry)
        assert_abs_diff_eq!(d2s[0][0], d2s[4][0], epsilon = 1e-14);
    }

    // =========================================================================
    // Test 8: Overlap second derivative (s|p) pair - HF/STO-3G
    // =========================================================================

    #[test]
    fn test_overlap_second_deriv_sp() {
        let basis = hf_sto3g();

        let shell_s = &basis.shells[0]; // H 1s
        let shell_p = &basis.shells[3]; // F 2p

        assert_eq!(shell_s.n_basis_functions(), 1);
        assert_eq!(shell_p.n_basis_functions(), 3);

        let d2s = shell_overlap_second_deriv(shell_s, shell_p);

        // (s|p) pair: each component gives a 1x3 matrix
        for comp in 0..9 {
            assert_eq!(
                d2s[comp].len(),
                3,
                "Component {} should be 1x3 = 3 elements",
                comp
            );
        }
    }

    // =========================================================================
    // Test 9: Overlap second derivative (p|p) pair - H2O/STO-3G
    // =========================================================================

    #[test]
    fn test_overlap_second_deriv_pp() {
        let basis = h2o_sto3g();

        // O has shells: 1s(idx 0), 2s(idx 1), 2p(idx 2)
        let shell_p = &basis.shells[2]; // O 2p

        assert_eq!(shell_p.n_basis_functions(), 3);

        // Self-pair: (p|p)
        let d2s = shell_overlap_second_deriv(shell_p, shell_p);

        // (p|p) pair: each component gives a 3x3 matrix
        for comp in 0..9 {
            assert_eq!(
                d2s[comp].len(),
                9,
                "Component {} should be 3x3 = 9 elements",
                comp
            );
        }

        // For same-center pair (p|p), d²S/dA_d·dB_e = d²S/dA_d·dA_e with sign
        // since A = B means displacement in A is same as displacement in B
        // For same-center: ∂²/∂A_d ∂A_e and ∂²/∂A_d ∂B_e differ by a sign pattern
    }

    // =========================================================================
    // Test 10: Kinetic second derivative (s|s) pair - H2/STO-3G
    // =========================================================================

    #[test]
    fn test_kinetic_second_deriv_ss() {
        let basis = h2_sto3g();

        let shell_a = &basis.shells[0];
        let shell_b = &basis.shells[1];

        let d2t = shell_kinetic_second_deriv(shell_a, shell_b);

        // For (s|s): 1x1 per component
        for comp in 0..9 {
            assert_eq!(d2t[comp].len(), 1);
        }

        // d²T/dA_z·dB_z should be nonzero
        assert!(
            d2t[8][0].abs() > 1e-5,
            "d²T/dA_z·dB_z for (s|s) should be nonzero, got {}",
            d2t[8][0]
        );

        // Axial symmetry: off-axis cross terms zero
        assert_abs_diff_eq!(d2t[0 * 3 + 1][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(d2t[0 * 3 + 2][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(d2t[1 * 3 + 0][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(d2t[1 * 3 + 2][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(d2t[2 * 3 + 0][0], 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(d2t[2 * 3 + 1][0], 0.0, epsilon = 1e-14);

        // Cylindrical symmetry: d²T/dA_x·dB_x = d²T/dA_y·dB_y
        assert_abs_diff_eq!(d2t[0][0], d2t[4][0], epsilon = 1e-14);
    }

    // =========================================================================
    // Test 11: Nuclear second derivative (s|s) pair - H2/STO-3G
    // =========================================================================

    #[test]
    fn test_nuclear_second_deriv_ss() {
        let basis = h2_sto3g();

        let shell_a = &basis.shells[0];
        let shell_b = &basis.shells[1];
        let c = &basis.atoms[0].position;
        let z = basis.atoms[0].atomic_number as f64;

        let d2v = shell_nuclear_second_deriv(shell_a, shell_b, c, z);

        // For (s|s): 1x1 per component
        for comp in 0..9 {
            assert_eq!(d2v[comp].len(), 1);
        }

        // d²V/dA_z·dB_z should be nonzero
        assert!(
            d2v[8][0].abs() > 1e-5,
            "d²V/dA_z·dB_z for (s|s) should be nonzero, got {}",
            d2v[8][0]
        );
    }

    // =========================================================================
    // Test 12: Overlap second derivative finite-difference validation
    //
    // CRITICAL: This is the definitive correctness test. For each (d,e),
    // compare analytical d²S/dA_d·dB_e against finite-difference of
    // dS/dA_d evaluated at displaced B positions:
    //
    //   d²S/dA_d·dB_e ≈ [dS/dA_d(B_e+h) - dS/dA_d(B_e-h)] / (2h)
    // =========================================================================

    #[test]
    fn test_overlap_second_deriv_finite_diff() {
        let basis = h2_sto3g();
        let step = 1e-5;

        let shell_a = &basis.shells[0]; // 1s on H(0,0,0)
        let shell_b = &basis.shells[1]; // 1s on H(0,0,1.4)

        let d2s = shell_overlap_second_deriv(shell_a, shell_b);

        for d in 0..3 {
            for e in 0..3 {
                // Displace shell_b center in direction e
                let mut center_plus = shell_b.center;
                center_plus[e] += step;
                let mut center_minus = shell_b.center;
                center_minus[e] -= step;

                let shell_b_plus = ContractedShell::new(
                    shell_b.angular_momentum,
                    shell_b.primitives.clone(),
                    center_plus,
                    shell_b.atom_idx,
                );
                let shell_b_minus = ContractedShell::new(
                    shell_b.angular_momentum,
                    shell_b.primitives.clone(),
                    center_minus,
                    shell_b.atom_idx,
                );

                // Compute dS/dA_d at displaced B positions
                let ds_plus = shell_overlap_first_deriv(shell_a, &shell_b_plus);
                let ds_minus = shell_overlap_first_deriv(shell_a, &shell_b_minus);

                for idx in 0..d2s[d * 3 + e].len() {
                    let fd = (ds_plus[d][idx] - ds_minus[d][idx]) / (2.0 * step);
                    let analytical = d2s[d * 3 + e][idx];
                    assert!(
                        (analytical - fd).abs() < 1e-6,
                        "d2S/dA_{}_dB_{} idx {}: analytical {:.8e} vs FD {:.8e} (err {:.2e})",
                        d,
                        e,
                        idx,
                        analytical,
                        fd,
                        (analytical - fd).abs(),
                    );
                }
            }
        }
    }

    // =========================================================================
    // Test 13: Overlap second derivative finite-diff for (s|p) pair
    // =========================================================================

    #[test]
    fn test_overlap_second_deriv_sp_finite_diff() {
        let basis = hf_sto3g();
        let step = 1e-5;

        let shell_s = &basis.shells[0]; // H 1s
        let shell_p = &basis.shells[3]; // F 2p

        let d2s = shell_overlap_second_deriv(shell_s, shell_p);

        for d in 0..3 {
            for e in 0..3 {
                let mut center_plus = shell_p.center;
                center_plus[e] += step;
                let mut center_minus = shell_p.center;
                center_minus[e] -= step;

                let shell_p_plus = ContractedShell::new(
                    shell_p.angular_momentum,
                    shell_p.primitives.clone(),
                    center_plus,
                    shell_p.atom_idx,
                );
                let shell_p_minus = ContractedShell::new(
                    shell_p.angular_momentum,
                    shell_p.primitives.clone(),
                    center_minus,
                    shell_p.atom_idx,
                );

                let ds_plus = shell_overlap_first_deriv(shell_s, &shell_p_plus);
                let ds_minus = shell_overlap_first_deriv(shell_s, &shell_p_minus);

                for idx in 0..d2s[d * 3 + e].len() {
                    let fd = (ds_plus[d][idx] - ds_minus[d][idx]) / (2.0 * step);
                    let analytical = d2s[d * 3 + e][idx];
                    assert!(
                        (analytical - fd).abs() < 1e-6,
                        "(s|p) d2S/dA_{}_dB_{} idx {}: analytical {:.8e} vs FD {:.8e} (err {:.2e})",
                        d,
                        e,
                        idx,
                        analytical,
                        fd,
                        (analytical - fd).abs(),
                    );
                }
            }
        }
    }

    // =========================================================================
    // Test 14: Kinetic second derivative finite-diff validation
    // =========================================================================

    #[test]
    fn test_kinetic_second_deriv_finite_diff() {
        let basis = h2_sto3g();
        let step = 1e-5;

        let shell_a = &basis.shells[0];
        let shell_b = &basis.shells[1];

        let d2t = shell_kinetic_second_deriv(shell_a, shell_b);

        for d in 0..3 {
            for e in 0..3 {
                let mut center_plus = shell_b.center;
                center_plus[e] += step;
                let mut center_minus = shell_b.center;
                center_minus[e] -= step;

                let shell_b_plus = ContractedShell::new(
                    shell_b.angular_momentum,
                    shell_b.primitives.clone(),
                    center_plus,
                    shell_b.atom_idx,
                );
                let shell_b_minus = ContractedShell::new(
                    shell_b.angular_momentum,
                    shell_b.primitives.clone(),
                    center_minus,
                    shell_b.atom_idx,
                );

                let dt_plus = shell_kinetic_first_deriv(shell_a, &shell_b_plus);
                let dt_minus = shell_kinetic_first_deriv(shell_a, &shell_b_minus);

                for idx in 0..d2t[d * 3 + e].len() {
                    let fd = (dt_plus[d][idx] - dt_minus[d][idx]) / (2.0 * step);
                    let analytical = d2t[d * 3 + e][idx];
                    assert!(
                        (analytical - fd).abs() < 1e-6,
                        "d2T/dA_{}_dB_{} idx {}: analytical {:.8e} vs FD {:.8e} (err {:.2e})",
                        d,
                        e,
                        idx,
                        analytical,
                        fd,
                        (analytical - fd).abs(),
                    );
                }
            }
        }
    }

    // =========================================================================
    // Test 15: Nuclear second derivative finite-diff validation
    // =========================================================================

    #[test]
    fn test_nuclear_second_deriv_finite_diff() {
        let basis = h2_sto3g();
        let step = 1e-5;

        let shell_a = &basis.shells[0];
        let shell_b = &basis.shells[1];
        let c = &basis.atoms[0].position;
        let z = basis.atoms[0].atomic_number as f64;

        let d2v = shell_nuclear_second_deriv(shell_a, shell_b, c, z);

        for d in 0..3 {
            for e in 0..3 {
                let mut center_plus = shell_b.center;
                center_plus[e] += step;
                let mut center_minus = shell_b.center;
                center_minus[e] -= step;

                let shell_b_plus = ContractedShell::new(
                    shell_b.angular_momentum,
                    shell_b.primitives.clone(),
                    center_plus,
                    shell_b.atom_idx,
                );
                let shell_b_minus = ContractedShell::new(
                    shell_b.angular_momentum,
                    shell_b.primitives.clone(),
                    center_minus,
                    shell_b.atom_idx,
                );

                let dv_plus = shell_nuclear_first_deriv_basis(shell_a, &shell_b_plus, c, z);
                let dv_minus = shell_nuclear_first_deriv_basis(shell_a, &shell_b_minus, c, z);

                for idx in 0..d2v[d * 3 + e].len() {
                    let fd = (dv_plus[d][idx] - dv_minus[d][idx]) / (2.0 * step);
                    let analytical = d2v[d * 3 + e][idx];
                    assert!(
                        (analytical - fd).abs() < 1e-6,
                        "d2V/dA_{}_dB_{} idx {}: analytical {:.8e} vs FD {:.8e} (err {:.2e})",
                        d,
                        e,
                        idx,
                        analytical,
                        fd,
                        (analytical - fd).abs(),
                    );
                }
            }
        }
    }

    // =========================================================================
    // Test 16: Overlap second derivative finite-diff for (p|p) pair
    //
    // Tests on H2O/STO-3G with O 2p shell to validate higher angular
    // momentum second derivatives.
    // =========================================================================

    #[test]
    fn test_overlap_second_deriv_pp_finite_diff() {
        let basis = h2o_sto3g();
        let step = 1e-5;

        // O 2p shell paired with H 1s shell (different centers)
        let shell_p = &basis.shells[2]; // O 2p
        let shell_h = &basis.shells[3]; // H1 1s

        let d2s = shell_overlap_second_deriv(shell_p, shell_h);

        for d in 0..3 {
            for e in 0..3 {
                let mut center_plus = shell_h.center;
                center_plus[e] += step;
                let mut center_minus = shell_h.center;
                center_minus[e] -= step;

                let shell_h_plus = ContractedShell::new(
                    shell_h.angular_momentum,
                    shell_h.primitives.clone(),
                    center_plus,
                    shell_h.atom_idx,
                );
                let shell_h_minus = ContractedShell::new(
                    shell_h.angular_momentum,
                    shell_h.primitives.clone(),
                    center_minus,
                    shell_h.atom_idx,
                );

                let ds_plus = shell_overlap_first_deriv(shell_p, &shell_h_plus);
                let ds_minus = shell_overlap_first_deriv(shell_p, &shell_h_minus);

                for idx in 0..d2s[d * 3 + e].len() {
                    let fd = (ds_plus[d][idx] - ds_minus[d][idx]) / (2.0 * step);
                    let analytical = d2s[d * 3 + e][idx];
                    assert!(
                        (analytical - fd).abs() < 1e-6,
                        "(p|s) d2S/dA_{}_dB_{} idx {}: analytical {:.8e} vs FD {:.8e} (err {:.2e})",
                        d,
                        e,
                        idx,
                        analytical,
                        fd,
                        (analytical - fd).abs(),
                    );
                }
            }
        }
    }

    // =========================================================================
    // Test 17: Nuclear second derivative (cross-center) for H2O with p-shells
    //
    // Tests shell_nuclear_second_deriv for all shell pairs and all nuclei
    // in H2O/STO-3G against central finite differences of
    // shell_nuclear_first_deriv_basis.
    //
    // This is the KEY diagnostic for the Hessian accuracy issue.
    // =========================================================================

    #[test]
    fn test_nuclear_second_deriv_h2o_cross_center_fd() {
        let basis = h2o_sto3g();
        let step = 1e-5;
        let n_shells = basis.shells.len();
        let mut max_err_overall = 0.0f64;

        for si in 0..n_shells {
            let shell_a = &basis.shells[si];
            let n_a = shell_a.n_basis_functions();

            for sj in 0..n_shells {
                let shell_b = &basis.shells[sj];
                let n_b = shell_b.n_basis_functions();

                for atom_c in &basis.atoms {
                    let z_c = atom_c.atomic_number as f64;
                    if z_c == 0.0 {
                        continue;
                    }
                    let r_c = &atom_c.position;

                    // Analytical cross-center second derivative: d²V/dA_d dB_e
                    let d2v = shell_nuclear_second_deriv(shell_a, shell_b, r_c, z_c);

                    for d in 0..3 {
                        for e in 0..3 {
                            // Finite difference: displace ket center (B) in direction e
                            let mut center_plus = shell_b.center;
                            center_plus[e] += step;
                            let mut center_minus = shell_b.center;
                            center_minus[e] -= step;

                            let shell_b_plus = ContractedShell::new(
                                shell_b.angular_momentum,
                                shell_b.primitives.clone(),
                                center_plus,
                                shell_b.atom_idx,
                            );
                            let shell_b_minus = ContractedShell::new(
                                shell_b.angular_momentum,
                                shell_b.primitives.clone(),
                                center_minus,
                                shell_b.atom_idx,
                            );

                            // Compute dV/dA_d at displaced B positions
                            let dv_plus =
                                shell_nuclear_first_deriv_basis(shell_a, &shell_b_plus, r_c, z_c);
                            let dv_minus =
                                shell_nuclear_first_deriv_basis(shell_a, &shell_b_minus, r_c, z_c);

                            for idx in 0..n_a * n_b {
                                let fd = (dv_plus[d][idx] - dv_minus[d][idx]) / (2.0 * step);
                                let analytical = d2v[d * 3 + e][idx];
                                let err = (analytical - fd).abs();
                                max_err_overall = max_err_overall.max(err);

                                if err > 1e-5 {
                                    eprintln!(
                                        "CROSS d2V/dA{}_dB{} shells({},{}) atom_c={} idx={}: \
                                         anal={:.8e} fd={:.8e} err={:.2e}",
                                        d,
                                        e,
                                        si,
                                        sj,
                                        atom_c.atomic_number,
                                        idx,
                                        analytical,
                                        fd,
                                        err,
                                    );
                                }

                                assert!(
                                    err < 1e-5,
                                    "d2V/dA{}_dB{} shells({},{}) atom_c={} idx={}: \
                                     anal={:.8e} fd={:.8e} err={:.2e}",
                                    d,
                                    e,
                                    si,
                                    sj,
                                    atom_c.atomic_number,
                                    idx,
                                    analytical,
                                    fd,
                                    err,
                                );
                            }
                        }
                    }
                }
            }
        }
        eprintln!(
            "H2O cross-center nuclear second deriv max error: {:.2e}",
            max_err_overall
        );
    }

    // =========================================================================
    // Test 18: Nuclear second derivative (bra-bra) for H2O with p-shells
    //
    // Tests shell_nuclear_second_deriv_bra_bra for all shell pairs and all
    // nuclei in H2O/STO-3G against central finite differences.
    //
    // FD: d²V/dA_d dA_e ≈ [dV/dA_d(A_e+h) - dV/dA_d(A_e-h)] / (2h)
    // where we displace the bra center in direction e and compute the
    // bra first derivative in direction d.
    // =========================================================================

    #[test]
    fn test_nuclear_second_deriv_h2o_bra_bra_fd() {
        let basis = h2o_sto3g();
        let step = 1e-5;
        let n_shells = basis.shells.len();
        let mut max_err_overall = 0.0f64;

        for si in 0..n_shells {
            let shell_a = &basis.shells[si];
            let n_a = shell_a.n_basis_functions();

            for sj in 0..n_shells {
                let shell_b = &basis.shells[sj];
                let n_b = shell_b.n_basis_functions();

                for atom_c in &basis.atoms {
                    let z_c = atom_c.atomic_number as f64;
                    if z_c == 0.0 {
                        continue;
                    }
                    let r_c = &atom_c.position;

                    // Analytical bra-bra second derivative: d²V/dA_d dA_e
                    let d2v = shell_nuclear_second_deriv_bra_bra(shell_a, shell_b, r_c, z_c);

                    for d in 0..3 {
                        for e in 0..3 {
                            // Finite difference: displace bra center (A) in direction e
                            let mut center_plus = shell_a.center;
                            center_plus[e] += step;
                            let mut center_minus = shell_a.center;
                            center_minus[e] -= step;

                            let shell_a_plus = ContractedShell::new(
                                shell_a.angular_momentum,
                                shell_a.primitives.clone(),
                                center_plus,
                                shell_a.atom_idx,
                            );
                            let shell_a_minus = ContractedShell::new(
                                shell_a.angular_momentum,
                                shell_a.primitives.clone(),
                                center_minus,
                                shell_a.atom_idx,
                            );

                            // Compute dV/dA_d at displaced A positions
                            let dv_plus =
                                shell_nuclear_first_deriv_basis(&shell_a_plus, shell_b, r_c, z_c);
                            let dv_minus =
                                shell_nuclear_first_deriv_basis(&shell_a_minus, shell_b, r_c, z_c);

                            for idx in 0..n_a * n_b {
                                let fd = (dv_plus[d][idx] - dv_minus[d][idx]) / (2.0 * step);
                                let analytical = d2v[d * 3 + e][idx];
                                let err = (analytical - fd).abs();
                                max_err_overall = max_err_overall.max(err);

                                // Use relative tolerance for large values
                                let denom = analytical.abs().max(1.0);
                                let rel_err = err / denom;

                                if err > 1e-4 || rel_err > 1e-6 {
                                    eprintln!(
                                        "BRABRA d2V/dA{}_dA{} shells({},{}) atom_c={} idx={}: \
                                         anal={:.8e} fd={:.8e} err={:.2e} rel={:.2e}",
                                        d,
                                        e,
                                        si,
                                        sj,
                                        atom_c.atomic_number,
                                        idx,
                                        analytical,
                                        fd,
                                        err,
                                        rel_err,
                                    );
                                }

                                assert!(
                                    err < 1e-4 && rel_err < 1e-5,
                                    "d2V/dA{}_dA{} shells({},{}) atom_c={} idx={}: \
                                     anal={:.8e} fd={:.8e} err={:.2e} rel={:.2e}",
                                    d,
                                    e,
                                    si,
                                    sj,
                                    atom_c.atomic_number,
                                    idx,
                                    analytical,
                                    fd,
                                    err,
                                    rel_err,
                                );
                            }
                        }
                    }
                }
            }
        }
        eprintln!(
            "H2O bra-bra nuclear second deriv max error: {:.2e}",
            max_err_overall
        );
    }
}
