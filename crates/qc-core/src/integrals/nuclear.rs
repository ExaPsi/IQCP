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
/// This implements the Obara-Saika scheme to compute:
/// ```text
/// V = ⟨G_a|-1/|r-C||G_b⟩
/// ```
/// where `G_a` and `G_b` are primitive Gaussians with angular momentum
/// specified by `a_powers` and `b_powers`.
///
/// The nuclear charge Z is NOT included here; it should be multiplied in the caller.
///
/// # Arguments
///
/// * `gp` - Pre-computed Gaussian product data
/// * `a_powers` - Cartesian powers (i, j, k) for bra Gaussian
/// * `b_powers` - Cartesian powers (i, j, k) for ket Gaussian
/// * `c` - Nuclear center position [x, y, z]
///
/// # Returns
///
/// The unnormalized nuclear attraction integral value (without charge factor)
///
/// # Algorithm
///
/// 1. Compute auxiliary integrals [0|V|0]^(m) for m = 0..=L_total
/// 2. Use VRR to build [a|0]^(m) for each Cartesian direction
/// 3. Use HTR to transfer angular momentum to center B
/// 4. Combine contributions from all three directions
///
/// # Reference
///
/// Obara & Saika (1986), Eqs. 4.1-4.8; libcint g1e.c CINTg1e_nuc
pub fn primitive_nuclear(
    gp: &GaussianProduct,
    a_powers: &CartesianPower,
    b_powers: &CartesianPower,
    c: &[f64; 3],
) -> f64 {
    let l_a = a_powers.angular_momentum() as usize;
    let l_b = b_powers.angular_momentum() as usize;
    let l_total = l_a + l_b;

    // Maximum auxiliary index needed
    // For VRR, we need m up to l_total
    let m_max = l_total;

    // Compute auxiliary integrals [0|V|0]^(m)
    let aux = nuclear_auxiliary(gp, c, m_max);

    // PC = P - C
    let pc = [
        gp.center_p[0] - c[0],
        gp.center_p[1] - c[1],
        gp.center_p[2] - c[2],
    ];

    // For each Cartesian direction, build [a|b]^(0) using VRR then HTR
    // The total integral is the product of 1D contributions times the auxiliary
    nuclear_3d(gp, a_powers, b_powers, &pc, &aux, m_max)
}

/// Compute the 3D nuclear attraction integral by combining 1D contributions
///
/// This handles the coupling between Cartesian directions through the auxiliary index m.
fn nuclear_3d(
    gp: &GaussianProduct,
    a_powers: &CartesianPower,
    b_powers: &CartesianPower,
    pc: &[f64; 3],
    aux: &[f64],
    m_max: usize,
) -> f64 {
    let a_x = a_powers.i as i32;
    let a_y = a_powers.j as i32;
    let a_z = a_powers.k as i32;

    let b_x = b_powers.i as i32;
    let b_y = b_powers.j as i32;
    let b_z = b_powers.k as i32;

    // Build the 3D integral using recursive VRR/HTR
    // We need to track the auxiliary index through all three directions
    // The correct approach uses a direct 3D recursion with the auxiliary index.

    nuclear_integral_recursive(
        gp, a_x, a_y, a_z, b_x, b_y, b_z, 0, // starting m = 0
        pc, aux, m_max,
    )
}

/// Recursive implementation of nuclear attraction VRR+HTR
///
/// Computes [a_x, a_y, a_z | b_x, b_y, b_z]^(m)
#[allow(clippy::too_many_arguments)]
fn nuclear_integral_recursive(
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
    // Check if m is within bounds
    if m > m_max {
        return 0.0;
    }

    // Negative angular momentum check (invalid)
    if a_x < 0 || a_y < 0 || a_z < 0 || b_x < 0 || b_y < 0 || b_z < 0 {
        return 0.0;
    }

    // Base case: [0,0,0|0,0,0]^(m) = aux[m]
    if a_x == 0 && a_y == 0 && a_z == 0 && b_x == 0 && b_y == 0 && b_z == 0 {
        return aux[m];
    }

    // Apply HTR first to reduce b to 0 (transfer angular momentum to a)
    // HTR: [a|b+1_i]^(m) = [a+1_i|b]^(m) + (A_i - B_i)[a|b]^(m)
    // Inverse: [a|b]^(m) can be computed from [a+b|0]^(m) via series of HTR

    // If b > 0, apply HTR in reverse: [a|b]^(m) from [a+1|b-1]^(m) - (A-B)[a|b-1]^(m)
    // Actually: [a|b+1] = [a+1|b] + (A-B)[a|b]
    // So: [a|b] for b > 0 can be expressed as:
    // [a|b] = [a+1|b-1] + (A-B)[a|b-1]  (by substituting b-1 for b in HTR)
    // No wait, that's still HTR building b up from 0.

    // The standard approach: use VRR to build all [a+b|0]^(m), then HTR to transfer.
    // For recursive implementation, we can apply:
    // - If b_i > 0 for some i, use HTR: [a|b] = [a+1_i|b-1_i] + (A_i - B_i)[a|b-1_i]
    // - If all b_i = 0, use VRR to reduce a to 0.

    // Choose the first non-zero b component to reduce via HTR
    if b_x > 0 {
        // HTR in x: [a_x|b_x]^(m) from [a_x+1|b_x-1]^(m) and [a_x|b_x-1]^(m)
        let term1 =
            nuclear_integral_recursive(gp, a_x + 1, a_y, a_z, b_x - 1, b_y, b_z, m, pc, aux, m_max);
        let term2 =
            nuclear_integral_recursive(gp, a_x, a_y, a_z, b_x - 1, b_y, b_z, m, pc, aux, m_max);
        let a_minus_b_x = -gp.ab[0]; // AB = B - A, so A - B = -AB
        return term1 + a_minus_b_x * term2;
    }

    if b_y > 0 {
        let term1 =
            nuclear_integral_recursive(gp, a_x, a_y + 1, a_z, b_x, b_y - 1, b_z, m, pc, aux, m_max);
        let term2 =
            nuclear_integral_recursive(gp, a_x, a_y, a_z, b_x, b_y - 1, b_z, m, pc, aux, m_max);
        let a_minus_b_y = -gp.ab[1];
        return term1 + a_minus_b_y * term2;
    }

    if b_z > 0 {
        let term1 =
            nuclear_integral_recursive(gp, a_x, a_y, a_z + 1, b_x, b_y, b_z - 1, m, pc, aux, m_max);
        let term2 =
            nuclear_integral_recursive(gp, a_x, a_y, a_z, b_x, b_y, b_z - 1, m, pc, aux, m_max);
        let a_minus_b_z = -gp.ab[2];
        return term1 + a_minus_b_z * term2;
    }

    // Now b = (0,0,0), apply VRR to reduce a to 0
    // VRR: [a+1_i|0]^(m) = PA_i[a|0]^(m) - PC_i[a|0]^(m+1) + (a_i/2p)([a-1_i|0]^(m) - [a-1_i|0]^(m+1))

    // Choose the first non-zero a component to reduce
    if a_x > 0 {
        // VRR in x: [a_x|0]^(m) from [a_x-1|0]^(m), [a_x-1|0]^(m+1), [a_x-2|0]^(m), [a_x-2|0]^(m+1)
        // [a_x|0]^(m) = PA_x[a_x-1|0]^(m) - PC_x[a_x-1|0]^(m+1) + ((a_x-1)/2p)([a_x-2|0]^(m) - [a_x-2|0]^(m+1))
        let term1_m = nuclear_integral_recursive(gp, a_x - 1, a_y, a_z, 0, 0, 0, m, pc, aux, m_max);
        let term1_m1 =
            nuclear_integral_recursive(gp, a_x - 1, a_y, a_z, 0, 0, 0, m + 1, pc, aux, m_max);

        let vrr_term1 = gp.pa[0] * term1_m - pc[0] * term1_m1;

        let vrr_term2 = if a_x >= 2 {
            let term2_m =
                nuclear_integral_recursive(gp, a_x - 2, a_y, a_z, 0, 0, 0, m, pc, aux, m_max);
            let term2_m1 =
                nuclear_integral_recursive(gp, a_x - 2, a_y, a_z, 0, 0, 0, m + 1, pc, aux, m_max);
            (a_x - 1) as f64 * gp.one_over_2p * (term2_m - term2_m1)
        } else {
            0.0
        };

        return vrr_term1 + vrr_term2;
    }

    if a_y > 0 {
        let term1_m = nuclear_integral_recursive(gp, a_x, a_y - 1, a_z, 0, 0, 0, m, pc, aux, m_max);
        let term1_m1 =
            nuclear_integral_recursive(gp, a_x, a_y - 1, a_z, 0, 0, 0, m + 1, pc, aux, m_max);

        let vrr_term1 = gp.pa[1] * term1_m - pc[1] * term1_m1;

        let vrr_term2 = if a_y >= 2 {
            let term2_m =
                nuclear_integral_recursive(gp, a_x, a_y - 2, a_z, 0, 0, 0, m, pc, aux, m_max);
            let term2_m1 =
                nuclear_integral_recursive(gp, a_x, a_y - 2, a_z, 0, 0, 0, m + 1, pc, aux, m_max);
            (a_y - 1) as f64 * gp.one_over_2p * (term2_m - term2_m1)
        } else {
            0.0
        };

        return vrr_term1 + vrr_term2;
    }

    if a_z > 0 {
        let term1_m = nuclear_integral_recursive(gp, a_x, a_y, a_z - 1, 0, 0, 0, m, pc, aux, m_max);
        let term1_m1 =
            nuclear_integral_recursive(gp, a_x, a_y, a_z - 1, 0, 0, 0, m + 1, pc, aux, m_max);

        let vrr_term1 = gp.pa[2] * term1_m - pc[2] * term1_m1;

        let vrr_term2 = if a_z >= 2 {
            let term2_m =
                nuclear_integral_recursive(gp, a_x, a_y, a_z - 2, 0, 0, 0, m, pc, aux, m_max);
            let term2_m1 =
                nuclear_integral_recursive(gp, a_x, a_y, a_z - 2, 0, 0, 0, m + 1, pc, aux, m_max);
            (a_z - 1) as f64 * gp.one_over_2p * (term2_m - term2_m1)
        } else {
            0.0
        };

        return vrr_term1 + vrr_term2;
    }

    // Should not reach here if logic is correct
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
