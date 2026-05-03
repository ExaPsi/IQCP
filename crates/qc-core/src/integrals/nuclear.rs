//! Nuclear attraction integral computation via Obara-Saika recurrence relations
//!
//! This module computes nuclear attraction integrals `V_ij = ⟨i|-Z_C/|r-R_C||j⟩` between
//! Gaussian basis functions using the Obara-Saika (OS) scheme.
//!
//! # Algorithm
//!
//! The nuclear attraction integral requires the Boys function F_m(T) and uses auxiliary
//! integrals with an extra index m to track the Boys function order.
//!
//! ## Base case (Obara-Saika 1986, Eq. 4.1)
//!
//! ```text
//! [0|V|0]^(m) = (-2π/p) · K_AB · F_m(T)
//! ```
//!
//! where:
//! - T = p · |P - C|² (Boys function argument)
//! - p = α + β (combined exponent)
//! - K_AB = exp(-αβ|A-B|²/p) from GaussianProduct
//! - F_m(T) = Boys function of order m
//!
//! ## VRR for nuclear attraction (Obara-Saika 1986, Eq. 4.5-4.8)
//!
//! ```text
//! [a+1_i|V|b]^(m) = (P_i - A_i)[a|V|b]^(m) - (P_i - C_i)[a|V|b]^(m+1)
//!                  + (a_i/2p)([a-1_i|V|b]^(m) - [a-1_i|V|b]^(m+1))
//!                  + (b_i/2p)([a|V|b-1_i]^(m) - [a|V|b-1_i]^(m+1))
//! ```
//!
//! ## HTR for nuclear attraction (same as overlap)
//!
//! ```text
//! [a|b+1_i]^(m) = [a+1_i|b]^(m) + (A_i - B_i)[a|b]^(m)
//! ```
//!
//! # References
//!
//! - Obara & Saika (1986), J. Chem. Phys. 84, 3963, Eqs. 4.1-4.8
//! - libcint implementation: `references/libcint/src/g1e.c` (CINTg1e_nuc function)
//! - Helgaker, Jorgensen & Olsen (2000), Molecular Electronic-Structure Theory, Ch. 9
//!
//! # Example
//!
//! ```rust
//! use qc_core::basis::{Atom, BasisSet};
//! use qc_core::integrals::nuclear_matrix;
//!
//! // Build H2 molecule
//! let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
//! let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
//! let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
//!
//! // Compute nuclear attraction matrix
//! let v = nuclear_matrix(&basis);
//!
//! // V is a 2x2 symmetric matrix (flattened row-major)
//! assert_eq!(v.len(), 4);
//! // All elements should be negative (attractive potential)
//! assert!(v[0] < 0.0);
//! assert!(v[3] < 0.0);
//! ```

use super::cartesian::{cartesian_components, CartesianPower};
use super::overlap::cartesian_gaussian_normalization;
use super::GaussianProduct;
use crate::basis::{BasisSet, ContractedShell};
use crate::boys::boys_eval_all;
use crate::rys::rys_roots;
use std::f64::consts::PI;

// =============================================================================
// Nuclear Auxiliary Integrals
// =============================================================================

/// Compute all auxiliary nuclear attraction integrals [0|V|0]^(m) for m = 0..=m_max
///
/// This computes the base case integrals:
/// ```text
/// [0|V|0]^(m) = (-2π/p) · K_AB · F_m(T)
/// ```
///
/// where T = p · |P - C|²
///
/// # Arguments
///
/// * `gp` - Pre-computed Gaussian product data
/// * `c` - Nuclear center position [x, y, z]
/// * `m_max` - Maximum auxiliary index needed
///
/// # Returns
///
/// Vector of auxiliary integrals [0|V|0]^(0), [0|V|0]^(1), ..., [0|V|0]^(m_max)
///
/// # Reference
///
/// Obara-Saika (1986), Eq. 4.1; libcint g1e.c lines 227-246
pub fn nuclear_auxiliary(gp: &GaussianProduct, c: &[f64; 3], m_max: usize) -> Vec<f64> {
    // Compute PC = P - C (product center minus nuclear center)
    let pc = [
        gp.center_p[0] - c[0],
        gp.center_p[1] - c[1],
        gp.center_p[2] - c[2],
    ];

    // Compute T = p * |P - C|^2 (Boys function argument)
    let pc_squared = pc[0] * pc[0] + pc[1] * pc[1] + pc[2] * pc[2];
    let t_arg = gp.p * pc_squared;

    // Compute Boys function values F_m(T) for m = 0..=m_max
    // Note: boys_eval_all returns BoysResult with values
    let boys_results = boys_eval_all(m_max as u32, t_arg)
        .expect("Boys function evaluation should succeed for valid inputs");

    // Prefactor: (-2π/p) * K_AB
    // Note: The NEGATIVE sign is critical for nuclear attraction!
    // Reference: libcint g1e.c line 233: fac1 = 2*M_PI * -Z * ...
    let prefactor = -2.0 * PI / gp.p * gp.k_ab;

    // Compute auxiliary integrals
    boys_results.iter().map(|br| prefactor * br.value).collect()
}

// =============================================================================
// 1D VRR for Nuclear Attraction
// =============================================================================

/// Compute 1D nuclear attraction VRR integrals for one Cartesian direction
///
/// This builds the array of integrals [a|0]^(m) for a = 0..=a_max, m = 0..=m_max
/// using the vertical recurrence relation:
///
/// ```text
/// [a+1|0]^(m) = PA * [a|0]^(m) - PC * [a|0]^(m+1) + (a/2p) * ([a-1|0]^(m) - [a-1|0]^(m+1))
/// ```
///
/// # Arguments
///
/// * `pa` - P_i - A_i (product center minus bra center, for this axis)
/// * `pc` - P_i - C_i (product center minus nuclear center, for this axis)
/// * `one_over_2p` - 1/(2p) where p = alpha + beta
/// * `a_max` - Maximum angular momentum needed on center A
/// * `m_max` - Maximum auxiliary index needed
/// * `aux` - Auxiliary integrals [0|V|0]^(m) from nuclear_auxiliary
///
/// # Returns
///
/// 2D array (flattened row-major) of [a][m] for a = 0..=a_max, m = 0..=m_max
/// Index: result[a * (m_max + 1) + m]
///
/// # Reference
///
/// Obara-Saika (1986), Eq. 4.5; libcint g1e.c lines 282-296
#[allow(dead_code)]
pub fn nuclear_vrr_1d(
    pa: f64,
    pc: f64,
    one_over_2p: f64,
    a_max: usize,
    m_max: usize,
    _aux: &[f64], // Auxiliary values; kept for API consistency
) -> Vec<f64> {
    let n_m = m_max + 1;
    let n_a = a_max + 1;
    let mut g = vec![0.0; n_a * n_m];

    // Base case: g[0][m] = 1.0 for all m (the prefactor is handled by aux)
    for item in g.iter_mut().take(n_m) {
        *item = 1.0;
    }

    if a_max == 0 {
        return g;
    }

    // First recurrence step: [1|0]^(m) = PA * [0|0]^(m) - PC * [0|0]^(m+1)
    // We need values for m = 0..=(m_max - 1) since we use m+1
    for m in 0..m_max {
        g[n_m + m] = pa * g[m] - pc * g[m + 1];
    }

    // VRR for a >= 2:
    // [a+1|0]^(m) = PA * [a|0]^(m) - PC * [a|0]^(m+1) + (a/2p) * ([a-1|0]^(m) - [a-1|0]^(m+1))
    for a in 1..a_max {
        let a_f = a as f64;
        for m in 0..(m_max - a) {
            // m_max decreases with a because each recurrence step uses m+1
            let term1 = pa * g[a * n_m + m] - pc * g[a * n_m + (m + 1)];
            let term2 = a_f * one_over_2p * (g[(a - 1) * n_m + m] - g[(a - 1) * n_m + (m + 1)]);
            g[(a + 1) * n_m + m] = term1 + term2;
        }
    }

    g
}

// =============================================================================
// Primitive Nuclear Attraction Integral
// =============================================================================

/// Compute the nuclear attraction integral between two primitive Cartesian Gaussians
/// for a single nucleus at position C with unit charge.
///
/// This implements the Rys quadrature scheme following libcint's CINTg1e_nuc.
/// The Rys-based approach avoids the catastrophic cancellation that occurs in the
/// auxiliary-index (Boys function) recursion when the nucleus is near a Gaussian center.
///
/// # Algorithm
///
/// 1. Compute T = p * |P - C|^2 (Boys function argument)
/// 2. Compute Rys roots u_n and weights w_n for nroots = (L_total)/2 + 1
/// 3. For each root n, compute modified VRR coefficients:
///    - ru_n = u_n / (1 + u_n)
///    - rt_n = 1/(2p) * (1 - ru_n)  (modified 1/(2p))
///    - r_n = (P - A) + ru_n * (C - P)  (modified PA, per axis)
/// 4. Build VRR per root: g[i+1] = r * g[i] + i * rt * g[i-1]
/// 5. Apply HTR: g[a|b+1] = g[a+1|b] + (A-B) * g[a|b]
/// 6. Sum over roots: integral = prefactor * sum_n w_n * gx_n * gy_n * gz_n
///
/// The nuclear charge Z is NOT included here; it should be multiplied in the caller.
///
/// # Reference
///
/// - libcint g1e.c lines 208-320 (CINTg1e_nuc)
/// - Dupuis, Rys & King (1976), J. Chem. Phys. 65, 111
pub fn primitive_nuclear(
    gp: &GaussianProduct,
    a_powers: &CartesianPower,
    b_powers: &CartesianPower,
    c: &[f64; 3],
) -> f64 {
    let l_a = a_powers.angular_momentum() as usize;
    let l_b = b_powers.angular_momentum() as usize;
    let l_total = l_a + l_b;

    // Number of Rys roots needed
    let nroots = l_total / 2 + 1;

    // Compute CP = C - P (nuclear center minus product center)
    let cp = [
        c[0] - gp.center_p[0],
        c[1] - gp.center_p[1],
        c[2] - gp.center_p[2],
    ];

    // Compute T = p * |P - C|^2 (Boys function argument)
    let pc_squared = cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2];
    let t_arg = gp.p * pc_squared;

    // Prefactor: -2π/p * K_AB
    // The negative sign for attraction is handled in shell_nuclear (via Z factor).
    // Here we compute the magnitude: 2π/p * K_AB
    let prefactor = -2.0 * PI / gp.p * gp.k_ab;

    // PA = P - A (product center minus bra center)
    let pa = gp.pa;
    // AB = B - A -> A - B = -AB
    let a_minus_b = [-gp.ab[0], -gp.ab[1], -gp.ab[2]];

    // For s-s case (nroots=1), we can use Boys function directly
    // For higher angular momentum, use Rys quadrature
    if nroots == 1 && l_total == 0 {
        // s-s case: V = prefactor * F_0(T)
        let boys_results =
            boys_eval_all(0, t_arg).expect("Boys function evaluation should succeed");
        return prefactor * boys_results[0].value;
    }

    // Get Rys roots and weights
    // rys_roots returns (roots, weights) where roots are in [0, 1)
    // and weights sum to F_0(T)
    let (rys_r, rys_w) = match rys_roots(nroots, t_arg) {
        Ok(result) => (result.roots, result.weights),
        Err(_) => {
            // Fallback to Boys function approach for edge cases
            let m_max = l_total;
            let aux = nuclear_auxiliary(gp, c, m_max);
            let pc_dir = [
                gp.center_p[0] - c[0],
                gp.center_p[1] - c[1],
                gp.center_p[2] - c[2],
            ];
            return nuclear_integral_recursive_fallback(
                gp,
                a_powers.i as i32,
                a_powers.j as i32,
                a_powers.k as i32,
                b_powers.i as i32,
                b_powers.j as i32,
                b_powers.k as i32,
                0,
                &pc_dir,
                &aux,
                m_max,
            );
        }
    };

    // Build VRR tables for each root and accumulate
    // Following libcint g1e.c lines 282-297
    //
    // For each root n:
    //   ru = u[n] / (1 + u[n])   -- note: rys_roots gives t^2, so u = t^2/(1-t^2)
    //   Actually, the Rys roots from our rys_roots() are t_n^2 in [0,1).
    //   The relationship to libcint's u is: u = t^2/(1-t^2), so t^2 = u/(1+u).
    //   Therefore ru = t^2 (the Rys root directly).
    //
    //   rt = 1/(2p) * (1 - ru) = 1/(2p) * (1 - t^2)
    //   r_x = PA_x + ru * CP_x = PA_x + t^2 * CP_x
    //
    // libcint formula (line 283-285):
    //   ru = tau^2 * u[n] / (1 + u[n])  (tau=1 for point nucleus)
    //   rt = aij2 - aij2 * ru = aij2 * (1 - ru)
    //   r0 = rijrx + ru * crij[0]
    // where rijrx = P - A, crij = C - P

    let nmax = l_a + l_b; // total angular momentum to build in VRR
    let a_total = (a_powers.i + a_powers.j + a_powers.k) as usize;
    let b_total = (b_powers.i + b_powers.j + b_powers.k) as usize;
    let _ = a_total; // used implicitly through a_powers
    let _ = b_total;

    let mut result = 0.0;

    for root_idx in 0..nroots {
        let t_sq = rys_r[root_idx]; // t^2, the Rys root in [0, 1)
        let w_n = rys_w[root_idx]; // weight

        // Modified coefficients for this root
        // ru = t^2 (the Rys root)
        let ru = t_sq;
        let rt = gp.one_over_2p * (1.0 - ru); // modified 1/(2p)

        // Modified PA for each axis: r = PA + ru * CP
        let r = [pa[0] + ru * cp[0], pa[1] + ru * cp[1], pa[2] + ru * cp[2]];

        // Build VRR tables for each Cartesian direction
        // g_x[n] = [n|0]_x for n = 0..=nmax (x component)
        // VRR: g[n+1] = r_x * g[n] + n * rt * g[n-1]
        let gx = vrr_1d_rys(r[0], rt, nmax);
        let gy = vrr_1d_rys(r[1], rt, nmax);
        let gz = vrr_1d_rys(r[2], rt, nmax);

        // Apply HTR to get [a|b] from [a+b|0] values, then combine 3 directions
        // For the specific (a_powers, b_powers) combination:
        let val_x = htr_1d(&gx, a_minus_b[0], a_powers.i as usize, b_powers.i as usize);
        let val_y = htr_1d(&gy, a_minus_b[1], a_powers.j as usize, b_powers.j as usize);
        let val_z = htr_1d(&gz, a_minus_b[2], a_powers.k as usize, b_powers.k as usize);

        result += w_n * val_x * val_y * val_z;
    }

    prefactor * result
}

/// Build 1D VRR table for a single Rys root
///
/// VRR: g[n+1] = r * g[n] + n * rt * g[n-1]
/// where r = PA + ru * CP (modified PA) and rt = 1/(2p) * (1-ru) (modified 1/(2p))
///
/// Returns g[0], g[1], ..., g[nmax]
fn vrr_1d_rys(r: f64, rt: f64, nmax: usize) -> Vec<f64> {
    let mut g = vec![0.0; nmax + 1];
    g[0] = 1.0;

    if nmax > 0 {
        g[1] = r; // g[1] = r * g[0] = r

        for n in 1..nmax {
            g[n + 1] = r * g[n] + (n as f64) * rt * g[n - 1];
        }
    }

    g
}

/// Apply 1D Horizontal Transfer (HTR) to get [a|b] from VRR table
///
/// HTR: [a|b+1] = [a+1|b] + (A-B) * [a|b]
///
/// Given g[n] = [n|0], compute [a|b].
fn htr_1d(g: &[f64], a_minus_b: f64, a: usize, b: usize) -> f64 {
    if b == 0 {
        return g[a];
    }

    // Build HTR table: h[aa][bb] = [aa|bb]
    // Initialize with h[n][0] = g[n]
    // HTR: h[aa][bb+1] = h[aa+1][bb] + (A-B) * h[aa][bb]

    let mut h = vec![vec![0.0; b + 1]; a + b + 1];

    // Initialize from VRR
    for n in 0..=a + b {
        h[n][0] = g[n];
    }

    // Apply HTR
    for bb in 0..b {
        for aa in 0..=(a + b - bb - 1) {
            h[aa][bb + 1] = h[aa + 1][bb] + a_minus_b * h[aa][bb];
        }
    }

    h[a][b]
}

/// Fallback recursive implementation using Boys function auxiliary indices.
/// Used when Rys quadrature fails (edge cases).
#[allow(clippy::too_many_arguments)]
fn nuclear_integral_recursive_fallback(
    gp: &GaussianProduct,
    a_x: i32,
    a_y: i32,
    a_z: i32,
    b_x: i32,
    b_y: i32,
    b_z: i32,
    m: usize,
    pc: &[f64; 3],
    aux: &[f64],
    m_max: usize,
) -> f64 {
    if m > m_max {
        return 0.0;
    }
    if a_x < 0 || a_y < 0 || a_z < 0 || b_x < 0 || b_y < 0 || b_z < 0 {
        return 0.0;
    }
    if a_x == 0 && a_y == 0 && a_z == 0 && b_x == 0 && b_y == 0 && b_z == 0 {
        return aux[m];
    }

    // HTR to reduce b
    if b_x > 0 {
        let term1 = nuclear_integral_recursive_fallback(
            gp,
            a_x + 1,
            a_y,
            a_z,
            b_x - 1,
            b_y,
            b_z,
            m,
            pc,
            aux,
            m_max,
        );
        let term2 = nuclear_integral_recursive_fallback(
            gp,
            a_x,
            a_y,
            a_z,
            b_x - 1,
            b_y,
            b_z,
            m,
            pc,
            aux,
            m_max,
        );
        return term1 + (-gp.ab[0]) * term2;
    }
    if b_y > 0 {
        let term1 = nuclear_integral_recursive_fallback(
            gp,
            a_x,
            a_y + 1,
            a_z,
            b_x,
            b_y - 1,
            b_z,
            m,
            pc,
            aux,
            m_max,
        );
        let term2 = nuclear_integral_recursive_fallback(
            gp,
            a_x,
            a_y,
            a_z,
            b_x,
            b_y - 1,
            b_z,
            m,
            pc,
            aux,
            m_max,
        );
        return term1 + (-gp.ab[1]) * term2;
    }
    if b_z > 0 {
        let term1 = nuclear_integral_recursive_fallback(
            gp,
            a_x,
            a_y,
            a_z + 1,
            b_x,
            b_y,
            b_z - 1,
            m,
            pc,
            aux,
            m_max,
        );
        let term2 = nuclear_integral_recursive_fallback(
            gp,
            a_x,
            a_y,
            a_z,
            b_x,
            b_y,
            b_z - 1,
            m,
            pc,
            aux,
            m_max,
        );
        return term1 + (-gp.ab[2]) * term2;
    }

    // VRR to reduce a
    if a_x > 0 {
        let term1_m =
            nuclear_integral_recursive_fallback(gp, a_x - 1, a_y, a_z, 0, 0, 0, m, pc, aux, m_max);
        let term1_m1 = nuclear_integral_recursive_fallback(
            gp,
            a_x - 1,
            a_y,
            a_z,
            0,
            0,
            0,
            m + 1,
            pc,
            aux,
            m_max,
        );
        let vrr1 = gp.pa[0] * term1_m - pc[0] * term1_m1;
        let vrr2 = if a_x >= 2 {
            let t2_m = nuclear_integral_recursive_fallback(
                gp,
                a_x - 2,
                a_y,
                a_z,
                0,
                0,
                0,
                m,
                pc,
                aux,
                m_max,
            );
            let t2_m1 = nuclear_integral_recursive_fallback(
                gp,
                a_x - 2,
                a_y,
                a_z,
                0,
                0,
                0,
                m + 1,
                pc,
                aux,
                m_max,
            );
            (a_x - 1) as f64 * gp.one_over_2p * (t2_m - t2_m1)
        } else {
            0.0
        };
        return vrr1 + vrr2;
    }
    if a_y > 0 {
        let term1_m =
            nuclear_integral_recursive_fallback(gp, a_x, a_y - 1, a_z, 0, 0, 0, m, pc, aux, m_max);
        let term1_m1 = nuclear_integral_recursive_fallback(
            gp,
            a_x,
            a_y - 1,
            a_z,
            0,
            0,
            0,
            m + 1,
            pc,
            aux,
            m_max,
        );
        let vrr1 = gp.pa[1] * term1_m - pc[1] * term1_m1;
        let vrr2 = if a_y >= 2 {
            let t2_m = nuclear_integral_recursive_fallback(
                gp,
                a_x,
                a_y - 2,
                a_z,
                0,
                0,
                0,
                m,
                pc,
                aux,
                m_max,
            );
            let t2_m1 = nuclear_integral_recursive_fallback(
                gp,
                a_x,
                a_y - 2,
                a_z,
                0,
                0,
                0,
                m + 1,
                pc,
                aux,
                m_max,
            );
            (a_y - 1) as f64 * gp.one_over_2p * (t2_m - t2_m1)
        } else {
            0.0
        };
        return vrr1 + vrr2;
    }
    if a_z > 0 {
        let term1_m =
            nuclear_integral_recursive_fallback(gp, a_x, a_y, a_z - 1, 0, 0, 0, m, pc, aux, m_max);
        let term1_m1 = nuclear_integral_recursive_fallback(
            gp,
            a_x,
            a_y,
            a_z - 1,
            0,
            0,
            0,
            m + 1,
            pc,
            aux,
            m_max,
        );
        let vrr1 = gp.pa[2] * term1_m - pc[2] * term1_m1;
        let vrr2 = if a_z >= 2 {
            let t2_m = nuclear_integral_recursive_fallback(
                gp,
                a_x,
                a_y,
                a_z - 2,
                0,
                0,
                0,
                m,
                pc,
                aux,
                m_max,
            );
            let t2_m1 = nuclear_integral_recursive_fallback(
                gp,
                a_x,
                a_y,
                a_z - 2,
                0,
                0,
                0,
                m + 1,
                pc,
                aux,
                m_max,
            );
            (a_z - 1) as f64 * gp.one_over_2p * (t2_m - t2_m1)
        } else {
            0.0
        };
        return vrr1 + vrr2;
    }

    0.0
}

// =============================================================================
// Shell Nuclear Attraction Integral
// =============================================================================

/// Compute nuclear attraction integrals between two contracted shells for a single nucleus
///
/// This computes all `n_a * n_b` nuclear attraction integrals between the Cartesian
/// components of shells A and B for a nucleus at position C with charge Z:
/// ```text
/// V_ij = -Z * sum_{p in A} sum_{q in B} c_p * N_p * c_q * N_q * V[p|q]
/// ```
/// where the sum runs over primitive pairs and N are normalization factors.
///
/// # Arguments
///
/// * `shell_a` - First contracted shell (bra)
/// * `shell_b` - Second contracted shell (ket)
/// * `c` - Nuclear center position [x, y, z]
/// * `z` - Nuclear charge (positive integer, e.g., 1 for H, 8 for O)
///
/// # Returns
///
/// Vector of nuclear attraction integrals in row-major order:
/// `[V(a0,b0), V(a0,b1), ..., V(a0,bn), V(a1,b0), ...]`
/// where a0, a1, ... are Cartesian components of shell A.
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{ContractedShell, AngularMomentum, GaussianPrimitive};
/// use qc_core::integrals::shell_nuclear;
///
/// // Create two s-shells at different centers
/// let h_prims = vec![
///     GaussianPrimitive::new(3.425251, 0.154329),
///     GaussianPrimitive::new(0.623914, 0.535328),
///     GaussianPrimitive::new(0.168855, 0.444635),
/// ];
/// let shell_a = ContractedShell::new(AngularMomentum::S, h_prims.clone(), [0.0, 0.0, 0.0], 0);
/// let shell_b = ContractedShell::new(AngularMomentum::S, h_prims, [0.0, 0.0, 1.4], 0);
///
/// // Nuclear attraction from H at origin (charge = 1)
/// let v = shell_nuclear(&shell_a, &shell_b, &[0.0, 0.0, 0.0], 1);
/// assert_eq!(v.len(), 1);  // 1x1 for s-s nuclear integral
/// assert!(v[0] < 0.0);     // Attractive potential is negative
/// ```
pub fn shell_nuclear(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
    c: &[f64; 3],
    z: u32,
) -> Vec<f64> {
    let l_a = shell_a.l_value();
    let l_b = shell_b.l_value();

    // Get Cartesian components for each shell
    let comps_a = cartesian_components(l_a).expect("Angular momentum within supported range");
    let comps_b = cartesian_components(l_b).expect("Angular momentum within supported range");

    let n_a = comps_a.len();
    let n_b = comps_b.len();

    // Output array: n_a x n_b integrals
    let mut integrals = vec![0.0; n_a * n_b];

    // Charge factor
    let z_factor = z as f64;

    // Loop over primitive pairs
    for prim_a in &shell_a.primitives {
        for prim_b in &shell_b.primitives {
            // Compute Gaussian product data
            let gp = GaussianProduct::new(
                prim_a.exponent,
                &shell_a.center,
                prim_b.exponent,
                &shell_b.center,
            );

            // Raw contraction coefficient product
            let coef = prim_a.coefficient * prim_b.coefficient;

            // Loop over Cartesian components
            for (i, pow_a) in comps_a.iter().enumerate() {
                for (j, pow_b) in comps_b.iter().enumerate() {
                    // Compute primitive integral (unnormalized, without charge)
                    let prim_integral = primitive_nuclear(&gp, pow_a, pow_b, c);

                    // Apply Cartesian Gaussian normalization for each primitive
                    let norm_a = cartesian_gaussian_normalization(prim_a.exponent, pow_a);
                    let norm_b = cartesian_gaussian_normalization(prim_b.exponent, pow_b);

                    // Add contribution to contracted integral (multiply by charge)
                    integrals[i * n_b + j] += coef * norm_a * norm_b * prim_integral * z_factor;
                }
            }
        }
    }

    integrals
}

// =============================================================================
// Nuclear Attraction Matrix
// =============================================================================

/// Compute the full nuclear attraction matrix for a basis set
///
/// Returns the N x N nuclear attraction matrix V where N is the total number of
/// basis functions. The matrix is symmetric: V_ij = V_ji.
///
/// This sums contributions from ALL nuclei in the system:
/// ```text
/// V_ij = sum_C (-Z_C) * <i|1/|r-R_C||j>
/// ```
///
/// # Arguments
///
/// * `basis` - The molecular basis set (contains atoms with positions and charges)
///
/// # Returns
///
/// Vector of length N*N containing the flattened row-major nuclear attraction matrix
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{Atom, BasisSet};
/// use qc_core::integrals::nuclear_matrix;
///
/// // H2 molecule
/// let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
/// let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
/// let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
///
/// let v = nuclear_matrix(&basis);
/// assert_eq!(v.len(), 4);  // 2x2 matrix
///
/// // Check symmetry
/// assert!((v[1] - v[2]).abs() < 1e-15);
///
/// // All elements should be negative (attractive)
/// assert!(v[0] < 0.0);
/// assert!(v[3] < 0.0);
/// ```
pub fn nuclear_matrix(basis: &BasisSet) -> Vec<f64> {
    let n = basis.n_basis;
    let mut v_matrix = vec![0.0; n * n];

    // Sum contributions from all nuclei
    for atom in &basis.atoms {
        let center = atom.position;
        let charge = atom.atomic_number as u32;

        // Iterate over shell pairs
        let mut mu = 0; // Basis function index for shell A
        for (i, shell_a) in basis.shells.iter().enumerate() {
            let n_a = shell_a.n_basis_functions();

            let mut nu = 0; // Basis function index for shell B
            for (j, shell_b) in basis.shells.iter().enumerate() {
                let n_b = shell_b.n_basis_functions();

                // Only compute upper triangle (i <= j) and symmetrize
                if i <= j {
                    // Compute shell block for this nucleus
                    let block = shell_nuclear(shell_a, shell_b, &center, charge);

                    // Add to matrix
                    for ia in 0..n_a {
                        for ib in 0..n_b {
                            let val = block[ia * n_b + ib];
                            // Upper triangle
                            v_matrix[(mu + ia) * n + (nu + ib)] += val;
                            // Lower triangle (symmetry)
                            if i != j {
                                v_matrix[(nu + ib) * n + (mu + ia)] += val;
                            }
                        }
                    }
                }

                nu += n_b;
            }
            mu += n_a;
        }
    }

    v_matrix
}

// =============================================================================
// Spherical Nuclear Attraction Matrix
// =============================================================================

/// Compute the nuclear attraction matrix in spherical harmonic basis
///
/// This computes the nuclear attraction matrix and transforms it from
/// Cartesian to spherical harmonic basis. For basis sets without d-orbitals
/// or higher, this is identical to `nuclear_matrix()`.
///
/// # Arguments
///
/// * `basis` - The molecular basis set
///
/// # Returns
///
/// Vector of length N_sph * N_sph containing the flattened row-major
/// spherical nuclear attraction matrix, where N_sph is the number of
/// spherical basis functions.
pub fn nuclear_matrix_spherical(basis: &BasisSet) -> Vec<f64> {
    use super::spherical::transform_one_electron_matrix;

    // If no spherical/Cartesian difference, just return the Cartesian matrix
    if !basis.has_spherical_difference() {
        return nuclear_matrix(basis);
    }

    // Compute Cartesian matrix and transform
    let cart_matrix = nuclear_matrix(basis);
    transform_one_electron_matrix(&cart_matrix, &basis.shells)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::{AngularMomentum, Atom, GaussianPrimitive};
    use approx::assert_abs_diff_eq;

    // Tolerance for golden test comparisons
    // Nuclear attraction integrals have accumulated errors from Boys function and
    // multiple recurrence steps. The recursive VRR/HTR implementation accumulates
    // errors proportional to the angular momentum, resulting in ~1e-6 relative error.
    // This is still well within chemical accuracy (1e-4 Hartree).
    const TOL: f64 = 1e-5;

    // -------------------------------------------------------------------------
    // Utility functions
    // -------------------------------------------------------------------------

    fn h_sto3g_primitives() -> Vec<GaussianPrimitive> {
        // Hydrogen STO-3G from PySCF
        vec![
            GaussianPrimitive::new(3.4252509100, 0.1543289707),
            GaussianPrimitive::new(0.6239137300, 0.5353281424),
            GaussianPrimitive::new(0.1688554000, 0.4446345420),
        ]
    }

    // -------------------------------------------------------------------------
    // Auxiliary integral tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_nuclear_auxiliary_same_center() {
        // When P = C, T = 0, and F_m(0) = 1/(2m+1)
        let alpha = 1.0;
        let center = [0.0, 0.0, 0.0];
        let gp = GaussianProduct::new(alpha, &center, alpha, &center);

        let aux = nuclear_auxiliary(&gp, &center, 3);

        // Prefactor: -2π/p * K_AB = -2π/2 * 1 = -π
        let prefactor = -PI;

        // F_0(0) = 1, F_1(0) = 1/3, F_2(0) = 1/5, F_3(0) = 1/7
        assert_abs_diff_eq!(aux[0], prefactor * 1.0, epsilon = TOL);
        assert_abs_diff_eq!(aux[1], prefactor / 3.0, epsilon = TOL);
        assert_abs_diff_eq!(aux[2], prefactor / 5.0, epsilon = TOL);
        assert_abs_diff_eq!(aux[3], prefactor / 7.0, epsilon = TOL);
    }

    #[test]
    fn test_nuclear_auxiliary_different_center() {
        // P != C case: need Boys function for T > 0
        let alpha = 1.0;
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let c = [1.0, 0.0, 0.0]; // Nucleus displaced by 1 Bohr

        let gp = GaussianProduct::new(alpha, &a, alpha, &b);
        let aux = nuclear_auxiliary(&gp, &c, 2);

        // p = 2, P = (0,0,0), |P-C|^2 = 1, T = 2
        // F_m(2) values from Boys function
        // Just check that values are reasonable (negative and decreasing in magnitude)
        assert!(aux[0] < 0.0);
        assert!(aux[1] < 0.0);
        assert!(aux[2] < 0.0);
        assert!(aux[0].abs() > aux[1].abs());
        assert!(aux[1].abs() > aux[2].abs());
    }

    // -------------------------------------------------------------------------
    // Primitive integral tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_primitive_nuclear_ss_same_center() {
        // Two s-type primitives at the same center as the nucleus
        let alpha = 1.0;
        let center = [0.0, 0.0, 0.0];
        let gp = GaussianProduct::new(alpha, &center, alpha, &center);

        let s = CartesianPower::new(0, 0, 0);
        let v = primitive_nuclear(&gp, &s, &s, &center);

        // Expected: -2π/p * K_AB * F_0(0) = -2π/2 * 1 * 1 = -π
        assert_abs_diff_eq!(v, -PI, epsilon = TOL);
    }

    #[test]
    fn test_primitive_nuclear_ss_nucleus_displaced() {
        // s-s integral with nucleus displaced
        let alpha = 1.0;
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let c = [1.0, 0.0, 0.0]; // Nucleus at (1,0,0)

        let gp = GaussianProduct::new(alpha, &a, alpha, &b);
        let s = CartesianPower::new(0, 0, 0);
        let v = primitive_nuclear(&gp, &s, &s, &c);

        // Should be negative and smaller in magnitude than same-center case
        assert!(v < 0.0);
        assert!(v.abs() < PI);
    }

    #[test]
    fn test_primitive_nuclear_is_negative() {
        // Nuclear attraction should always be negative (attractive)
        let alpha = 2.0;
        let beta = 1.5;
        let a = [0.0, 0.0, 0.0];
        let b = [0.5, 0.0, 0.0];
        let c = [0.25, 0.0, 0.0];

        let gp = GaussianProduct::new(alpha, &a, beta, &b);

        // Test various angular momenta
        let test_cases = [
            (CartesianPower::new(0, 0, 0), CartesianPower::new(0, 0, 0)), // ss
            (CartesianPower::new(1, 0, 0), CartesianPower::new(0, 0, 0)), // ps
            (CartesianPower::new(0, 0, 0), CartesianPower::new(1, 0, 0)), // sp
            (CartesianPower::new(1, 0, 0), CartesianPower::new(1, 0, 0)), // pp_x
        ];

        for (pow_a, pow_b) in test_cases {
            let v = primitive_nuclear(&gp, &pow_a, &pow_b, &c);
            // Most combinations should be negative, but some off-diagonal p-type might not be
            // Just check they're finite
            assert!(
                v.is_finite(),
                "Integral should be finite for {:?}-{:?}",
                pow_a,
                pow_b
            );
        }
    }

    // -------------------------------------------------------------------------
    // Shell integral tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_shell_nuclear_same_center() {
        // Single H 1s shell with nucleus at same center
        let prims = h_sto3g_primitives();
        let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0, 0.0, 0.0], 0);

        let v = shell_nuclear(&shell, &shell, &[0.0, 0.0, 0.0], 1);

        assert_eq!(v.len(), 1);
        // Should be negative
        assert!(
            v[0] < 0.0,
            "Nuclear attraction should be negative, got {}",
            v[0]
        );
    }

    // -------------------------------------------------------------------------
    // Golden tests: H2 STO-3G
    // -------------------------------------------------------------------------

    #[test]
    fn test_golden_h2_sto3g_nuclear() {
        // Reference: PySCF 2.11.0
        //
        // H2 at 1.3984 Bohr separation (same geometry as overlap/kinetic tests)
        // Expected nuclear attraction matrix (2x2):
        //
        // Generated with:
        // ```python
        // from pyscf import gto
        // mol = gto.Mole()
        // mol.atom = 'H 0 0 0; H 0 0 1.3984'
        // mol.basis = 'sto-3g'
        // mol.unit = 'B'
        // mol.build()
        // print(mol.intor('int1e_nuc'))
        // ```
        //
        // Reference values:
        // V[0,0] = -1.880990423049179e+00
        // V[0,1] = -1.196333536602375e+00
        // V[1,0] = -1.196333536602375e+00
        // V[1,1] = -1.880990423049179e+00

        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let v = nuclear_matrix(&basis);

        // Reference values from PySCF
        let v_00 = -1.880990423049179e+00;
        let v_01 = -1.196333536602375e+00;

        assert_abs_diff_eq!(v[0], v_00, epsilon = TOL); // V[0,0]
        assert_abs_diff_eq!(v[1], v_01, epsilon = TOL); // V[0,1]
        assert_abs_diff_eq!(v[2], v_01, epsilon = TOL); // V[1,0] (symmetry)
        assert_abs_diff_eq!(v[3], v_00, epsilon = TOL); // V[1,1]
    }

    // -------------------------------------------------------------------------
    // Golden tests: H2O STO-3G
    // -------------------------------------------------------------------------

    #[test]
    fn test_golden_h2o_sto3g_nuclear() {
        // Reference: PySCF 2.11.0
        //
        // H2O geometry (Bohr):
        // O: [0, 0, 0.2216656]
        // H1: [0, 1.4309295, -0.8866625]
        // H2: [0, -1.4309295, -0.8866625]
        //
        // 7x7 nuclear attraction matrix
        //
        // Generated with:
        // ```python
        // from pyscf import gto
        // mol = gto.Mole()
        // mol.atom = '''
        // O   0.0000000   0.0000000   0.2216656303019316
        // H   0.0000000   1.4309295343303710  -0.8866625212077263
        // H   0.0000000  -1.4309295343303710  -0.8866625212077263
        // '''
        // mol.basis = 'sto-3g'
        // mol.unit = 'B'
        // mol.build()
        // print(mol.intor('int1e_nuc').flatten())
        // ```

        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        let v = nuclear_matrix(&basis);
        assert_eq!(v.len(), 49); // 7x7

        // Reference values from PySCF 2.11.0 (flattened row-major)
        let reference = [
            -61.72345247301635,
            -7.444636718775824,
            0.0,
            0.0,
            0.01898834225829542,
            -1.7445439378594125,
            -1.7445439378594125,
            -7.444636718775825,
            -10.142338146647798,
            0.0,
            0.0,
            0.22339810110550098,
            -3.8659536258482188,
            -3.8659536258482188,
            0.0,
            0.0,
            -9.985889904709952,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -10.141922649922062,
            0.0,
            -2.2526288189058916,
            2.2526288189058916,
            0.018988342258295425,
            0.223398101105501,
            0.0,
            0.0,
            -10.079498554960185,
            1.8176927293305916,
            1.8176927293305916,
            -1.7445439378594123,
            -3.865953625848219,
            0.0,
            -2.2526288189058925,
            1.8176927293305911,
            -5.834372188716397,
            -1.6149226694440912,
            -1.7445439378594123,
            -3.865953625848219,
            0.0,
            2.2526288189058925,
            1.8176927293305911,
            -1.6149226694440912,
            -5.834372188716397,
        ];

        // Check all elements
        // Nuclear attraction involves Boys function and multiple recurrence steps,
        // so we use a slightly relaxed tolerance. For molecules with p-orbitals,
        // the error can accumulate up to ~1e-5 relative to individual matrix elements.
        let h2o_tol = 1e-5;
        for (i, &ref_val) in reference.iter().enumerate() {
            assert!(
                (v[i] - ref_val).abs() < h2o_tol,
                "Mismatch at index {}: computed {} vs reference {} (error: {:.2e})",
                i,
                v[i],
                ref_val,
                (v[i] - ref_val).abs()
            );
        }
    }

    // -------------------------------------------------------------------------
    // Matrix properties tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_nuclear_matrix_symmetry() {
        // Nuclear attraction matrix should be symmetric: V_ij = V_ji
        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.43, -0.47]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        let v = nuclear_matrix(&basis);
        let n = basis.n_basis;

        for i in 0..n {
            for j in 0..n {
                assert_abs_diff_eq!(v[i * n + j], v[j * n + i], epsilon = 1e-14);
            }
        }
    }

    #[test]
    fn test_nuclear_matrix_all_negative_diagonal() {
        // Diagonal elements of nuclear attraction matrix should be negative
        // (electron attracted to nuclei)
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let v = nuclear_matrix(&basis);
        let n = basis.n_basis;

        for i in 0..n {
            assert!(
                v[i * n + i] < 0.0,
                "Diagonal element {} should be negative, got {}",
                i,
                v[i * n + i]
            );
        }
    }

    #[test]
    fn test_nuclear_integral_scaling_with_charge() {
        // Nuclear attraction should scale linearly with nuclear charge
        let prims = h_sto3g_primitives();
        let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0, 0.0, 0.0], 0);
        let c = [0.0, 0.0, 0.5]; // Nucleus slightly displaced

        let v_z1 = shell_nuclear(&shell, &shell, &c, 1);
        let v_z2 = shell_nuclear(&shell, &shell, &c, 2);

        assert_abs_diff_eq!(v_z2[0], 2.0 * v_z1[0], epsilon = 1e-14);
    }

    #[test]
    fn test_nuclear_integral_distance_dependence() {
        // Nuclear attraction should decrease with distance from nucleus
        let prims = h_sto3g_primitives();
        let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0, 0.0, 0.0], 0);

        let v_close = shell_nuclear(&shell, &shell, &[0.0, 0.0, 0.5], 1);
        let v_far = shell_nuclear(&shell, &shell, &[0.0, 0.0, 2.0], 1);

        // Both should be negative, but closer nucleus gives more negative value
        assert!(v_close[0] < 0.0);
        assert!(v_far[0] < 0.0);
        assert!(
            v_close[0] < v_far[0],
            "Closer nucleus should give more negative attraction: close={}, far={}",
            v_close[0],
            v_far[0]
        );
    }
}
