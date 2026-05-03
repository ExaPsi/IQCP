//! Overlap integral computation via Obara-Saika recurrence relations
//!
//! This module computes overlap integrals `S_ij = <i|j>` between Gaussian basis
//! functions using the efficient Obara-Saika (OS) recurrence relations.
//!
//! # Algorithm
//!
//! The computation proceeds in three steps:
//!
//! 1. **Base integral**: `[0|0] = (pi/p)^{3/2} * K_AB` (from Gaussian product theorem)
//!
//! 2. **Vertical Recurrence (VRR)**: Build angular momentum on center A
//!    ```text
//!    [a+1_i|b] = (P_i - A_i)[a|b] + (a_i/2p)[a-1_i|b] + (b_i/2p)[a|b-1_i]
//!    ```
//!
//! 3. **Horizontal Transfer (HTR)**: Transfer angular momentum from A to B
//!    ```text
//!    [a|b+1_i] = [a+1_i|b] + (A_i - B_i)[a|b]
//!    ```
//!
//! # References
//!
//! - Obara & Saika (1986), J. Chem. Phys. 84, 3963
//! - Head-Gordon & Pople (1988), J. Chem. Phys. 89, 5777
//! - libcint implementation: `references/libcint/src/g1e.c` lines 125-184
//!
//! # Example
//!
//! ```rust
//! use qc_core::basis::{Atom, BasisSet};
//! use qc_core::integrals::overlap_matrix;
//!
//! // Build H2 molecule
//! let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
//! let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
//! let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
//!
//! // Compute overlap matrix
//! let s = overlap_matrix(&basis);
//!
//! // S is a 2x2 symmetric matrix (flattened row-major)
//! assert_eq!(s.len(), 4);
//! // Diagonal elements are ~1.0 (normalized basis functions)
//! assert!((s[0] - 1.0).abs() < 1e-9);  // S[0,0]
//! assert!((s[3] - 1.0).abs() < 1e-9);  // S[1,1]
//! ```

use super::cartesian::{cartesian_components, CartesianPower};
use super::GaussianProduct;
use crate::basis::{BasisSet, ContractedShell};
use std::f64::consts::PI;

// =============================================================================
// Primitive Overlap via Obara-Saika Recurrence
// =============================================================================

/// Compute the overlap integral between two primitive Cartesian Gaussians
///
/// This implements the Obara-Saika recurrence relations to compute:
/// ```text
/// [a|b] = <G_a|G_b>
/// ```
/// where `G_a` and `G_b` are primitive Gaussians with angular momentum
/// specified by `a_powers` and `b_powers`.
///
/// # Arguments
///
/// * `gp` - Pre-computed Gaussian product data
/// * `a_powers` - Cartesian powers (i, j, k) for bra Gaussian
/// * `b_powers` - Cartesian powers (i, j, k) for ket Gaussian
///
/// # Returns
///
/// The overlap integral value
///
/// # Algorithm
///
/// For each Cartesian direction x, y, z:
/// 1. Use VRR to build `[a|0]` from `[0|0]`
/// 2. Use HTR to build `[a|b]` from `[a+b|0]`
/// 3. Multiply contributions from all three directions
///
/// # Reference
///
/// libcint `g1e.c` lines 164-182: VRR builds all angular momentum on one center,
/// then HTR transfers to the other center.
pub fn primitive_overlap(
    gp: &GaussianProduct,
    a_powers: &CartesianPower,
    b_powers: &CartesianPower,
) -> f64 {
    // Base case: [0|0] = (pi/p)^{3/2} * K_AB (already stored in gp.ss_integral)
    // For the OS scheme, we compute each Cartesian direction separately

    // Compute 1D overlap integrals for each direction
    let s_x = overlap_1d(
        gp.pa[0],          // P_x - A_x
        gp.ab[0],          // A_x - B_x (note: AB = B - A, so we negate)
        gp.one_over_2p,    // 1/(2p)
        a_powers.i as i32, // angular momentum on A
        b_powers.i as i32, // angular momentum on B
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

    // Total primitive overlap = [s|s] * S_x * S_y * S_z
    gp.ss_integral * s_x * s_y * s_z
}

/// Compute 1D overlap integral component using Obara-Saika recurrence
///
/// This computes the 1D component of the overlap integral:
/// ```text
/// S_1D(a, b) = <x^a * exp(-alpha*(x-A)^2) | x^b * exp(-beta*(x-B)^2)>
/// ```
///
/// # Arguments
///
/// * `pa` - P - A (product center minus bra center)
/// * `ab` - B - A (ket center minus bra center)
/// * `one_over_2p` - 1/(2p) where p = alpha + beta
/// * `a` - Angular momentum on center A
/// * `b` - Angular momentum on center B
///
/// # Algorithm
///
/// 1. VRR to build [n|0] for n = 0, 1, ..., a+b
/// 2. HTR to build [a|b] from [a+b|0], [a+b-1|0], etc.
///
/// Reference: libcint g1e.c lines 164-182
pub(crate) fn overlap_1d(pa: f64, ab: f64, one_over_2p: f64, a: i32, b: i32) -> f64 {
    // Handle base cases
    if a < 0 || b < 0 {
        return 0.0;
    }

    let total = (a + b) as usize;

    // Build VRR table: g[n] = [n|0] for n = 0, 1, ..., a+b
    // VRR: [n+1|0] = PA * [n|0] + n/(2p) * [n-1|0]
    // Reference: libcint g1e.c lines 169-173

    // Pre-allocate array for VRR values
    let mut g = vec![0.0; total + 1];

    // Base case: [0|0] = 1 (the prefactor is handled separately)
    g[0] = 1.0;

    if total > 0 {
        // [1|0] = PA * [0|0] = PA
        g[1] = pa;

        // VRR: [n+1|0] = PA * [n|0] + n/(2p) * [n-1|0]
        // Reference: libcint line 170-173:
        //   gx[(i+1)*di] = i * aij2 * gx[(i-1)*di] + rijrx[0] * gx[i*di];
        for n in 1..total {
            g[n + 1] = pa * g[n] + (n as f64) * one_over_2p * g[n - 1];
        }
    }

    // Now apply HTR to get [a|b] from the g[n] = [n|0] values
    // HTR: [a|b+1] = [a+1|b] + AB * [a|b]
    // where AB = B - A = ab
    //
    // Reference: libcint g1e.c lines 175-182:
    //   gx[n] = gx[n+di-dj] + rirj[0] * gx[n-dj];
    // Note: rirj = ri - rj = A - B = -ab, but libcint uses ri - rj
    //
    // We need to be careful: the HTR relation is
    //   [a|b+1] = [a+1|b] + (A-B)[a|b]
    //
    // Since ab = B - A, we have A - B = -ab

    if b == 0 {
        // No HTR needed, just return [a|0]
        return g[a as usize];
    }

    // Build HTR table
    // We need [a|b] where a and b are the input angular momenta
    // Start with [a+b|0], [a+b-1|0], etc. and work down to [a|b]
    //
    // Create 2D table: h[aa][bb] = [aa|bb]
    // Initialize with h[n][0] = g[n] for n = 0..=a+b
    // Then fill using: h[aa][bb+1] = h[aa+1][bb] + (A-B) * h[aa][bb]

    let a_u = a as usize;
    let b_u = b as usize;

    // For HTR, we iterate from b=0 up to b
    // At each step, we need h[aa][bb] for aa = a-bb+1 down to a-bb
    // Simplify: use 1D array, iterating over b and updating

    // h[i] represents [a - b + i | b_current]
    // Actually, let's be more explicit with 2D for clarity

    let mut h = vec![vec![0.0; b_u + 1]; a_u + b_u + 1];

    // Initialize [n|0] values from VRR
    for n in 0..=a_u + b_u {
        h[n][0] = g[n];
    }

    // HTR: [aa|bb+1] = [aa+1|bb] + (A-B)[aa|bb]
    // A - B = -ab
    let a_minus_b = -ab;

    for bb in 0..b_u {
        for aa in 0..=(a_u + b_u - bb - 1) {
            h[aa][bb + 1] = h[aa + 1][bb] + a_minus_b * h[aa][bb];
        }
    }

    h[a_u][b_u]
}

// =============================================================================
// Shell Overlap
// =============================================================================

/// Compute overlap integrals between two contracted shells
///
/// This computes all `n_a * n_b` overlap integrals between the Cartesian
/// components of shells A and B:
/// ```text
/// S_ij = sum_{p in A} sum_{q in B} c_p * N_p * c_q * N_q * [p|q]
/// ```
/// where the sum runs over primitive pairs and N are normalization factors.
///
/// # Arguments
///
/// * `shell_a` - First contracted shell (bra)
/// * `shell_b` - Second contracted shell (ket)
///
/// # Returns
///
/// Vector of overlap integrals in row-major order:
/// `[S(a0,b0), S(a0,b1), ..., S(a0,bn), S(a1,b0), ...]`
/// where a0, a1, ... are Cartesian components of shell A.
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{ContractedShell, AngularMomentum, GaussianPrimitive};
/// use qc_core::integrals::shell_overlap;
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
/// let s = shell_overlap(&shell_a, &shell_b);
/// assert_eq!(s.len(), 1);  // 1x1 for s-s overlap
/// ```
pub fn shell_overlap(shell_a: &ContractedShell, shell_b: &ContractedShell) -> Vec<f64> {
    let l_a = shell_a.l_value();
    let l_b = shell_b.l_value();

    // Get Cartesian components for each shell
    let comps_a = cartesian_components(l_a).expect("Angular momentum within supported range");
    let comps_b = cartesian_components(l_b).expect("Angular momentum within supported range");

    let n_a = comps_a.len();
    let n_b = comps_b.len();

    // Output array: n_a x n_b integrals
    let mut integrals = vec![0.0; n_a * n_b];

    // The basis set coefficients from standard libraries (like PySCF's STO-3G)
    // are raw contraction coefficients, NOT including the primitive Gaussian
    // normalization factor. We must include the normalization:
    //
    // N(alpha, i, j, k) = (2*alpha/pi)^(3/4) * (4*alpha)^((i+j+k)/2)
    //                    / sqrt((2i-1)!! * (2j-1)!! * (2k-1)!!)
    //
    // For s-type (i=j=k=0): N = (2*alpha/pi)^(3/4)

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
                    // Compute primitive integral (unnormalized)
                    let prim_integral = primitive_overlap(&gp, pow_a, pow_b);

                    // Apply Cartesian Gaussian normalization for each primitive
                    let norm_a = cartesian_gaussian_normalization(prim_a.exponent, pow_a);
                    let norm_b = cartesian_gaussian_normalization(prim_b.exponent, pow_b);

                    // Add contribution to contracted integral
                    integrals[i * n_b + j] += coef * norm_a * norm_b * prim_integral;
                }
            }
        }
    }

    integrals
}

/// Compute normalization factor for a Cartesian Gaussian primitive
///
/// For a primitive Gaussian with exponent alpha and Cartesian powers (i, j, k):
/// ```text
/// N(alpha, i, j, k) = (2*alpha/pi)^(3/4) * (4*alpha)^((i+j+k)/2)
///                    / sqrt((2i-1)!! * (2j-1)!! * (2k-1)!!)
/// ```
///
/// For s-type (i=j=k=0): N = (2*alpha/pi)^(3/4)
///
/// Reference: Szabo & Ostlund, Eq. 3.203; PySCF normalization conventions
#[inline]
pub(crate) fn cartesian_gaussian_normalization(alpha: f64, powers: &CartesianPower) -> f64 {
    // Base normalization: (2*alpha/pi)^(3/4)
    // Rewrite as ((2α/π)^3)^(1/4) = sqrt(sqrt((2α/π)^3)) to avoid powf()
    let two_a_pi = 2.0 * alpha / PI;
    let base = (two_a_pi * two_a_pi * two_a_pi).sqrt().sqrt();

    let i = powers.i;
    let j = powers.j;
    let k = powers.k;
    let l = i + j + k;

    if l == 0 {
        // s-type: just the base normalization
        return base;
    }

    // Angular part: (4*alpha)^(L/2) / sqrt((2i-1)!! * (2j-1)!! * (2k-1)!!)
    // Compute (4α)^(L/2) without powf() using explicit cases
    let four_alpha = 4.0 * alpha;
    let ang_num = match l {
        1 => four_alpha.sqrt(),
        2 => four_alpha,
        3 => four_alpha * four_alpha.sqrt(),
        4 => four_alpha * four_alpha,
        _ => four_alpha.powf(l as f64 / 2.0), // fallback for L > 4
    };

    // Double factorial denominator for each axis
    let denom_i = double_factorial(if i > 0 { 2 * i - 1 } else { 0 }) as f64;
    let denom_j = double_factorial(if j > 0 { 2 * j - 1 } else { 0 }) as f64;
    let denom_k = double_factorial(if k > 0 { 2 * k - 1 } else { 0 }) as f64;

    let ang_denom = (denom_i * denom_j * denom_k).sqrt();

    base * ang_num / ang_denom
}

/// Compute double factorial (2n-1)!! = 1 * 3 * 5 * ... * (2n-1)
///
/// Special cases:
/// - (-1)!! = 1 (by convention)
/// - 0!! = 1
/// - 1!! = 1
#[inline]
pub(crate) fn double_factorial(n: u32) -> u64 {
    if n <= 1 {
        1
    } else {
        let mut result = 1u64;
        let mut k = n;
        while k > 1 {
            result *= k as u64;
            k -= 2;
        }
        result
    }
}

/// Compute PySCF-compatible primitive normalization for ERI computation.
///
/// PySCF uses different normalization conventions for different angular momenta:
/// - For s and p orbitals (l <= 1): Cartesian Gaussian normalization
/// - For d and higher orbitals (l >= 2): gto_norm(l, alpha) / 2
///
/// This function returns the normalization factor that makes IQCP's ERIs
/// match PySCF's output exactly.
///
/// # Arguments
///
/// * `l` - Angular momentum quantum number
/// * `alpha` - Gaussian exponent
/// * `powers` - Cartesian powers (i, j, k) for the orbital
///
/// # Returns
///
/// Normalization factor for ERI computation that matches PySCF conventions.
///
/// # Reference
///
/// - PySCF gto_norm: `pyscf/gto/mole.py`
/// - Verified against PySCF 2.11.0 with mol.cart=True
#[inline]
pub fn pyscf_eri_normalization(l: u32, alpha: f64, powers: &CartesianPower) -> f64 {
    if l <= 1 {
        // For s and p orbitals, use standard Cartesian normalization
        // This gives overlap = 1.0 for all s and p components
        cartesian_gaussian_normalization(alpha, powers)
    } else {
        // For d and higher orbitals, PySCF uses gto_norm which is UNIFORM
        // for all Cartesian components of the shell.
        // This does NOT give overlap = 1.0, but matches PySCF's convention.
        // Reference: PySCF gto_norm in pyscf/gto/mole.py
        radial_gto_norm(l, alpha)
    }
}

/// Compute radial GTO normalization factor (PySCF gto_norm convention)
///
/// This is PySCF's gto_norm function, which normalizes the radial part
/// `g = r^l * exp(-alpha * r^2)` such that the integral of g^2 * r^2
/// over all space equals 1.
///
/// Formula:
/// ```text
/// gto_norm(l, alpha) = 2^((l+2)/2) * (2*alpha)^((l+1.5)/2) / sqrt((2l+1)!!) / pi^0.25
/// ```
///
/// # Arguments
///
/// * `l` - Angular momentum quantum number
/// * `alpha` - Gaussian exponent
///
/// # Returns
///
/// Normalization factor that is the same for all Cartesian components of shell L.
///
/// # Reference
///
/// - PySCF gto_norm: `pyscf/gto/mole.py`
/// - H. B. Schlegel and M. J. Frisch, Int. J. Quant. Chem., 54(1995), 83-87
#[inline]
pub fn radial_gto_norm(l: u32, alpha: f64) -> f64 {
    // Match PySCF's gto_norm exactly:
    // gto_norm(l, alpha) = 1 / sqrt(gaussian_int(2*l + 2, 2*alpha))
    //
    // where gaussian_int(n, alpha) = gamma((n+1)/2) / (2 * alpha^((n+1)/2))
    //
    // Using gamma(l + 1.5) = (2l+1)!! * sqrt(pi) / 2^(l+1):
    // gto_norm = 2^((l+2)/2) * (2*alpha)^((l+1.5)/2) / sqrt((2l+1)!!) / pi^0.25

    let l_f64 = l as f64;
    let two_alpha = 2.0 * alpha;
    let double_fact = double_factorial(2 * l + 1) as f64;

    // 2^((l+2)/2)
    let factor1 = 2.0_f64.powf((l_f64 + 2.0) / 2.0);

    // (2*alpha)^((l+1.5)/2)
    let factor2 = two_alpha.powf((l_f64 + 1.5) / 2.0);

    // 1 / (sqrt((2l+1)!!) * pi^0.25)
    let denominator = double_fact.sqrt() * PI.powf(0.25);

    factor1 * factor2 / denominator
}

/// Compute the Gaussian integral formula used by PySCF
///
/// This computes:
/// ```text
/// gaussian_int(n, alpha) = gamma((n+1)/2) / (2 * alpha^((n+1)/2))
/// ```
///
/// For even n: gamma((n+1)/2) = sqrt(pi) * (n-1)!! / 2^(n/2)
/// For odd n:  gamma((n+1)/2) = ((n-1)/2)!
///
/// # Arguments
///
/// * `n` - The power parameter (typically 2*l + 2 for contracted overlap)
/// * `alpha` - The sum of two exponents
///
/// # Reference
///
/// PySCF gaussian_int in pyscf/gto/mole.py
#[inline]
pub fn gaussian_int(n: u32, alpha: f64) -> f64 {
    // gamma((n+1)/2) / (2 * alpha^((n+1)/2))
    let half_n_plus_1 = (n as f64 + 1.0) / 2.0;

    // Use std::f64::consts::PI and gamma via exp(lgamma)
    // For numerical stability, we compute gamma using the relation:
    // gamma(x) = (x-1)! for positive integers, or use log-gamma for half-integers

    // For half-integer x = k + 0.5:
    // gamma(k + 0.5) = (2k-1)!! * sqrt(pi) / 2^k
    // For integer x = k:
    // gamma(k) = (k-1)!

    let n_mod_2 = n % 2;
    let gamma_val = if n_mod_2 == 0 {
        // n is even: (n+1)/2 = k + 0.5 for k = n/2
        // gamma(k + 0.5) = (2k-1)!! * sqrt(pi) / 2^k
        let k = n / 2;
        let double_fact = if k == 0 {
            1.0
        } else {
            double_factorial(2 * k - 1) as f64
        };
        double_fact * PI.sqrt() / (1u64 << k) as f64
    } else {
        // n is odd: (n+1)/2 = k (integer) for k = (n+1)/2
        // gamma(k) = (k-1)!
        let k = n.div_ceil(2);
        factorial(k - 1) as f64
    };

    gamma_val / (2.0 * alpha.powf(half_n_plus_1))
}

/// Compute factorial n!
#[inline]
fn factorial(n: u32) -> u64 {
    if n <= 1 {
        1
    } else {
        (2..=n as u64).product()
    }
}

/// Compute the contracted shell renormalization factor for PySCF compatibility
///
/// PySCF's `_nomalize_contracted_ao` function computes a renormalization factor
/// that makes the contracted self-overlap exactly 1.0 in the gto_norm convention.
///
/// The contracted self-overlap using gto_norm'd coefficients is:
/// ```text
/// S_contracted = sum_{p,q} c[p] * gto_norm(l, alpha[p]) * gaussian_int(2l+2, alpha[p]+alpha[q])
///                        * c[q] * gto_norm(l, alpha[q])
/// ```
///
/// The renormalization factor is: `1 / sqrt(S_contracted)`
///
/// For perfectly normalized primitives, S_contracted would be exactly 1.0, but due to
/// the gto_norm convention, it's typically very close to but not exactly 1.0 (e.g., 1.0000002).
/// This small difference creates a ~1e-7 relative error per coefficient, which accumulates
/// to ~4e-7 in ERIs (four coefficients).
///
/// # Arguments
///
/// * `shell` - The contracted shell
///
/// # Returns
///
/// The shell renormalization factor to multiply with integrals.
///
/// # Reference
///
/// - PySCF `_nomalize_contracted_ao` in pyscf/gto/mole.py (lines ~995-1030)
/// - This factor ensures IQCP ERIs match PySCF at 1e-10 or better tolerance.
pub fn shell_renorm_factor(shell: &crate::basis::ContractedShell) -> f64 {
    let l = shell.l_value();
    let n = 2 * l + 2;

    // Compute the contracted self-overlap in gto_norm convention
    let mut s_contracted = 0.0;

    for p in &shell.primitives {
        let gto_norm_p = radial_gto_norm(l, p.exponent);
        let c_p = p.coefficient * gto_norm_p;

        for q in &shell.primitives {
            let gto_norm_q = radial_gto_norm(l, q.exponent);
            let c_q = q.coefficient * gto_norm_q;

            // gaussian_int(2l+2, alpha_p + alpha_q)
            let ee = p.exponent + q.exponent;
            let integral = gaussian_int(n, ee);

            s_contracted += c_p * integral * c_q;
        }
    }

    // Renormalization factor: 1 / sqrt(S_contracted)
    1.0 / s_contracted.sqrt()
}

// =============================================================================
// Overlap Matrix
// =============================================================================

/// Compute the full overlap matrix for a basis set
///
/// Returns the N x N overlap matrix S where N is the total number of
/// basis functions. The matrix is symmetric: S_ij = S_ji.
///
/// # Arguments
///
/// * `basis` - The molecular basis set
///
/// # Returns
///
/// Vector of length N*N containing the flattened row-major overlap matrix
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{Atom, BasisSet};
/// use qc_core::integrals::overlap_matrix;
///
/// // H2 molecule
/// let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
/// let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
/// let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
///
/// let s = overlap_matrix(&basis);
/// assert_eq!(s.len(), 4);  // 2x2 matrix
///
/// // Check symmetry
/// assert!((s[1] - s[2]).abs() < 1e-15);
///
/// // Check normalization (diagonal should be ~1)
/// assert!((s[0] - 1.0).abs() < 1e-9);
/// assert!((s[3] - 1.0).abs() < 1e-9);
/// ```
pub fn overlap_matrix(basis: &BasisSet) -> Vec<f64> {
    let n = basis.n_basis;
    let mut s_matrix = vec![0.0; n * n];

    // Iterate over shell pairs
    let mut mu = 0; // Basis function index for shell A
    for (i, shell_a) in basis.shells.iter().enumerate() {
        let n_a = shell_a.n_basis_functions();

        let mut nu = 0; // Basis function index for shell B
        for (j, shell_b) in basis.shells.iter().enumerate() {
            let n_b = shell_b.n_basis_functions();

            // Only compute upper triangle (i <= j) and symmetrize
            if i <= j {
                // Compute shell block
                let block = shell_overlap(shell_a, shell_b);

                // Copy to matrix
                for ia in 0..n_a {
                    for ib in 0..n_b {
                        let val = block[ia * n_b + ib];
                        // Upper triangle
                        s_matrix[(mu + ia) * n + (nu + ib)] = val;
                        // Lower triangle (symmetry)
                        if i != j {
                            s_matrix[(nu + ib) * n + (mu + ia)] = val;
                        }
                    }
                }
            }

            nu += n_b;
        }
        mu += n_a;
    }

    s_matrix
}

// =============================================================================
// Spherical Overlap Matrix
// =============================================================================

/// Compute the overlap matrix in spherical harmonic basis
///
/// This computes the overlap matrix and transforms it from Cartesian to
/// spherical harmonic basis. For basis sets without d-orbitals or higher,
/// this is identical to `overlap_matrix()`.
///
/// # Arguments
///
/// * `basis` - The molecular basis set
///
/// # Returns
///
/// Vector of length N_sph * N_sph containing the flattened row-major
/// spherical overlap matrix, where N_sph is the number of spherical
/// basis functions.
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{Atom, BasisSet};
/// use qc_core::integrals::overlap_matrix_spherical;
///
/// // H2O with 6-31G* (has d-polarization on O)
/// let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
/// let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();
/// let h2 = Atom::new(1, [0.0, -1.43, -0.47]).unwrap();
///
/// let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();
/// let s_sph = overlap_matrix_spherical(&basis);
///
/// // Spherical basis has fewer functions than Cartesian for d-orbitals
/// assert_eq!(s_sph.len(), basis.n_basis_spherical() * basis.n_basis_spherical());
/// ```
pub fn overlap_matrix_spherical(basis: &BasisSet) -> Vec<f64> {
    use super::spherical::transform_one_electron_matrix;

    // If no spherical/Cartesian difference, just return the Cartesian matrix
    if !basis.has_spherical_difference() {
        return overlap_matrix(basis);
    }

    // Compute Cartesian matrix and transform
    let cart_matrix = overlap_matrix(basis);
    transform_one_electron_matrix(&cart_matrix, &basis.shells)
}

// =============================================================================
// Overlap vs. Distance Scan
// =============================================================================

/// Evaluate the overlap integral between two contracted shells at multiple
/// interatomic distances.
///
/// For each distance R in `r_values`, this function places `shell_a` at the
/// origin `[0, 0, 0]` and `shell_b` at `[0, 0, R]`, then computes their
/// overlap integral using `shell_overlap`.
///
/// # Arguments
///
/// * `shell_a` - First contracted shell (will be placed at origin)
/// * `shell_b` - Second contracted shell (will be placed at `[0, 0, R]`)
/// * `r_values` - Interatomic distances (in Bohr) at which to evaluate
///
/// # Returns
///
/// Vector of overlap integral values `S_ab(R)`, one per distance in `r_values`.
/// For s-s overlaps, this is the single overlap element. For higher angular
/// momentum, this returns the overlap of the first Cartesian component of each
/// shell (e.g., `px-px` for p-p overlap, `s-px` for s-p overlap).
///
/// # Postconditions
///
/// * Return vector has length `r_values.len()`
/// * At R=0 with identical shells: `S_ab ~ 1.0` (within normalization precision)
/// * At R -> inf: `S_ab -> 0`
///
/// # References
///
/// * Obara & Saika (1986), J. Chem. Phys. 84, 3963
/// * Existing `shell_overlap` implementation in this module
pub fn evaluate_overlap_vs_distance(
    shell_a: &ContractedShell,
    shell_b: &ContractedShell,
    r_values: &[f64],
) -> Vec<f64> {
    r_values
        .iter()
        .map(|&r| {
            // Place shell_a at origin, shell_b at [0, 0, R]
            let mut a = shell_a.clone();
            a.center = [0.0, 0.0, 0.0];

            let mut b = shell_b.clone();
            b.center = [0.0, 0.0, r];

            // Compute shell overlap block and extract [0, 0] element
            // (first Cartesian component of each shell)
            let block = shell_overlap(&a, &b);
            block[0]
        })
        .collect()
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
    // The 1e-9 tolerance accounts for accumulated floating-point errors
    // in normalization and contraction from coefficients stored at different
    // precision levels. This still represents high accuracy for quantum chemistry.
    // Implementation-level tests use stricter tolerances where appropriate.
    const TOL: f64 = 1e-9;

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
    // Double factorial tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_double_factorial() {
        assert_eq!(double_factorial(0), 1);
        assert_eq!(double_factorial(1), 1);
        assert_eq!(double_factorial(3), 3); // 3!! = 3
        assert_eq!(double_factorial(5), 15); // 5!! = 5 * 3 = 15
        assert_eq!(double_factorial(7), 105); // 7!! = 7 * 5 * 3 = 105
    }

    // -------------------------------------------------------------------------
    // 1D overlap tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_overlap_1d_ss() {
        // [0|0] should be 1.0 (base case, prefactor handled separately)
        let val = overlap_1d(0.5, 1.0, 0.25, 0, 0);
        assert_abs_diff_eq!(val, 1.0, epsilon = TOL);
    }

    #[test]
    fn test_overlap_1d_ps() {
        // [1|0] = PA
        let pa = 0.3;
        let val = overlap_1d(pa, 1.0, 0.25, 1, 0);
        assert_abs_diff_eq!(val, pa, epsilon = TOL);
    }

    #[test]
    fn test_overlap_1d_sp() {
        // [0|1] uses HTR: [0|1] = [1|0] + (A-B)[0|0]
        // With PA = 0.3, AB = B-A = 1.0, so A-B = -1.0
        // [0|1] = 0.3 + (-1.0) * 1.0 = -0.7
        let pa = 0.3;
        let ab = 1.0;
        let val = overlap_1d(pa, ab, 0.25, 0, 1);
        assert_abs_diff_eq!(val, pa - ab, epsilon = TOL);
    }

    #[test]
    fn test_overlap_1d_pp() {
        // Test [1|1] case
        // [2|0] = PA * [1|0] + 1/(2p) * [0|0]
        // [1|1] = [2|0] + (A-B)[1|0]
        let pa = 0.3;
        let ab = 1.0;
        let one_over_2p = 0.25;

        let val = overlap_1d(pa, ab, one_over_2p, 1, 1);

        // Manual calculation
        let g0 = 1.0;
        let g1 = pa; // 0.3
        let g2 = pa * g1 + one_over_2p * g0; // 0.09 + 0.25 = 0.34
        let a_minus_b = -ab; // -1.0
        let h_1_1 = g2 + a_minus_b * g1; // 0.34 - 0.3 = 0.04

        assert_abs_diff_eq!(val, h_1_1, epsilon = TOL);
    }

    // -------------------------------------------------------------------------
    // Primitive overlap tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_primitive_overlap_ss_same_center() {
        // Two s-type primitives at the same center
        let alpha = 1.0;
        let center = [0.0, 0.0, 0.0];
        let gp = GaussianProduct::new(alpha, &center, alpha, &center);

        let s = CartesianPower::new(0, 0, 0);
        let val = primitive_overlap(&gp, &s, &s);

        // [s|s] = (pi/p)^(3/2) for same center
        let p = 2.0 * alpha;
        let expected = (PI / p).powf(1.5);
        assert_abs_diff_eq!(val, expected, epsilon = TOL);
    }

    #[test]
    fn test_primitive_overlap_ss_different_centers() {
        // Two s-type primitives at different centers
        let alpha = 1.0;
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 1.0];
        let gp = GaussianProduct::new(alpha, &a, alpha, &b);

        let s = CartesianPower::new(0, 0, 0);
        let val = primitive_overlap(&gp, &s, &s);

        // [s|s] = (pi/p)^(3/2) * K_AB
        let p = 2.0 * alpha;
        let mu = alpha * alpha / p;
        let ab_sq = 1.0;
        let k_ab = (-mu * ab_sq).exp();
        let expected = (PI / p).powf(1.5) * k_ab;

        assert_abs_diff_eq!(val, expected, epsilon = TOL);
    }

    // -------------------------------------------------------------------------
    // Shell overlap tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_shell_overlap_same_center() {
        // Single H 1s shell overlapping with itself
        let prims = h_sto3g_primitives();
        let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0, 0.0, 0.0], 0);

        let s = shell_overlap(&shell, &shell);

        assert_eq!(s.len(), 1);
        // Should be normalized to 1.0 (within numerical precision of contraction)
        assert_abs_diff_eq!(s[0], 1.0, epsilon = TOL);
    }

    // -------------------------------------------------------------------------
    // Golden tests: H2 STO-3G
    // -------------------------------------------------------------------------

    #[test]
    fn test_golden_h2_sto3g_overlap() {
        // Reference: tests/golden/integrals/h2_sto3g_overlap.json
        // Generated from PySCF 2.11.0
        //
        // H2 at 1.3984 Bohr separation
        // Expected overlap matrix (2x2):
        // [1.0, 0.6598721980070731]
        // [0.6598721980070731, 1.0]

        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let s = overlap_matrix(&basis);

        // Reference values from PySCF
        let s_00 = 1.0000000000000002;
        let s_01 = 0.6598721980070731;

        assert_abs_diff_eq!(s[0], s_00, epsilon = TOL); // S[0,0]
        assert_abs_diff_eq!(s[1], s_01, epsilon = TOL); // S[0,1]
        assert_abs_diff_eq!(s[2], s_01, epsilon = TOL); // S[1,0] (symmetry)
        assert_abs_diff_eq!(s[3], s_00, epsilon = TOL); // S[1,1]
    }

    // -------------------------------------------------------------------------
    // Golden tests: H2O STO-3G
    // -------------------------------------------------------------------------

    #[test]
    fn test_golden_h2o_sto3g_overlap() {
        // Reference: tests/golden/integrals/h2o_sto3g_overlap.json
        // Generated from PySCF 2.11.0
        //
        // H2O geometry (Bohr):
        // O: [0, 0, 0.2216656]
        // H1: [0, 1.4309295, -0.8866625]
        // H2: [0, -1.4309295, -0.8866625]
        //
        // 7x7 overlap matrix

        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        let s = overlap_matrix(&basis);
        assert_eq!(s.len(), 49); // 7x7

        // Reference values from PySCF (flattened row-major)
        let reference = [
            1.0,
            0.2367039365108476,
            0.0,
            0.0,
            0.0,
            0.053902244303870125,
            0.053902244303870125,
            0.23670393651084762,
            1.0,
            0.0,
            0.0,
            0.0,
            0.47442024240584385,
            0.47442024240584385,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.3108959396001613,
            -0.3108959396001613,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            -0.24080481517921498,
            -0.24080481517921498,
            0.05390224430387014,
            0.47442024240584385,
            0.0,
            0.3108959396001614,
            -0.24080481517921493,
            1.0000000000000002,
            0.251525112694529,
            0.05390224430387014,
            0.47442024240584385,
            0.0,
            -0.3108959396001614,
            -0.24080481517921493,
            0.251525112694529,
            1.0000000000000002,
        ];

        // Check all elements
        for (i, &ref_val) in reference.iter().enumerate() {
            assert!(
                (s[i] - ref_val).abs() < TOL,
                "Mismatch at index {}: computed {} vs reference {}",
                i,
                s[i],
                ref_val
            );
        }
    }

    // -------------------------------------------------------------------------
    // Matrix properties tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_overlap_matrix_symmetry() {
        // Overlap matrix should be symmetric: S_ij = S_ji
        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.43, -0.47]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        let s = overlap_matrix(&basis);
        let n = basis.n_basis;

        for i in 0..n {
            for j in 0..n {
                assert_abs_diff_eq!(s[i * n + j], s[j * n + i], epsilon = 1e-14);
            }
        }
    }

    #[test]
    fn test_overlap_matrix_positive_definite() {
        // Overlap matrix should be positive definite (all eigenvalues > 0)
        // A quick check: diagonal elements should be positive
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let s = overlap_matrix(&basis);
        let n = basis.n_basis;

        for i in 0..n {
            assert!(
                s[i * n + i] > 0.0,
                "Diagonal element {} should be positive",
                i
            );
        }
    }

    #[test]
    fn test_overlap_matrix_normalized_diagonal() {
        // For normalized basis functions, diagonal should be close to 1.0
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let s = overlap_matrix(&basis);
        let n = basis.n_basis;

        for i in 0..n {
            assert_abs_diff_eq!(s[i * n + i], 1.0, epsilon = TOL);
        }
    }

    // -------------------------------------------------------------------------
    // gaussian_int and shell_renorm_factor tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_gaussian_int_values() {
        // Reference values from PySCF/scipy:
        // gaussian_int(n, alpha) = gamma((n+1)/2) / (2 * alpha^((n+1)/2))

        // n=2, alpha=1.0: gamma(1.5) / (2 * 1^1.5) = sqrt(pi)/2 / 2 = sqrt(pi)/4
        let val = super::gaussian_int(2, 1.0);
        let expected = PI.sqrt() / 4.0; // 0.4431134627263791
        assert_abs_diff_eq!(val, expected, epsilon = 1e-12);

        // n=2, alpha=2.0: gamma(1.5) / (2 * 2^1.5) = sqrt(pi)/2 / (2*2.828...)
        let val = super::gaussian_int(2, 2.0);
        let expected = PI.sqrt() / 2.0 / (2.0 * 2.0_f64.powf(1.5));
        assert_abs_diff_eq!(val, expected, epsilon = 1e-12);

        // n=4, alpha=1.0: gamma(2.5) / (2 * 1^2.5) = (3/2)*sqrt(pi)/2 / 2
        // gamma(2.5) = 1.5 * gamma(1.5) = 1.5 * sqrt(pi)/2 = 0.75*sqrt(pi)
        let val = super::gaussian_int(4, 1.0);
        let expected = 0.75 * PI.sqrt() / 2.0; // = 3*sqrt(pi)/8
        assert_abs_diff_eq!(val, expected, epsilon = 1e-12);
    }

    #[test]
    fn test_shell_renorm_h_631g_inner() {
        // Test shell_renorm_factor for H 6-31G* inner s shell (3 primitives)
        // Reference from Python: shell_renorm = 0.999999985249434

        use crate::basis::ContractedShell;

        let prims = vec![
            GaussianPrimitive::new(18.7311370, 0.03349460434),
            GaussianPrimitive::new(2.8253937, 0.23472695355),
            GaussianPrimitive::new(0.6401217, 0.81375732610),
        ];
        let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0, 0.0, 0.0], 0);

        let renorm = super::shell_renorm_factor(&shell);

        // The shell renorm factor should be very close to 1.0
        // Python gives: 0.999999985249434
        assert_abs_diff_eq!(renorm, 0.999999985249434, epsilon = 1e-12);
    }

    #[test]
    fn test_shell_renorm_h_631g_outer() {
        // Test shell_renorm_factor for H 6-31G* outer s shell (1 primitive)
        // For a single-primitive shell, the renorm factor should be exactly 1.0

        use crate::basis::ContractedShell;

        let prims = vec![GaussianPrimitive::new(0.1612778, 1.0)];
        let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0, 0.0, 0.0], 0);

        let renorm = super::shell_renorm_factor(&shell);

        // For a single primitive with coefficient 1.0, the renorm should be exactly 1.0
        assert_abs_diff_eq!(renorm, 1.0, epsilon = 1e-12);
    }

    // -------------------------------------------------------------------------
    // Spherical harmonic transformation tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_overlap_spherical_dimension_reduction() {
        // Test dimension reduction: 6-31G* has d-orbitals which reduce from 6 to 5 functions
        //
        // For H2O with 6-31G*:
        // O: 1s(1), 2s(1), 3s(1), 2px(1), 2py(1), 2pz(1), 3px(1), 3py(1), 3pz(1), d(6 cart/5 sph)
        //    = 9 + 6 = 15 Cartesian, 9 + 5 = 14 spherical
        // H: 1s(1), 2s(1) = 2 each, total 4
        // Total: Cartesian = 15 + 4 = 19, Spherical = 14 + 4 = 18

        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();

        // Verify Cartesian vs spherical dimension
        let n_cart = basis.n_basis;
        let n_sph = basis.n_basis_spherical();

        assert_eq!(
            n_cart, 19,
            "H2O 6-31G* should have 19 Cartesian basis functions"
        );
        assert_eq!(
            n_sph, 18,
            "H2O 6-31G* should have 18 spherical basis functions"
        );
        assert!(
            basis.has_spherical_difference(),
            "Basis should have d-orbitals"
        );

        // Compute matrices
        let s_cart = overlap_matrix(&basis);
        let s_sph = overlap_matrix_spherical(&basis);

        // Verify dimensions
        assert_eq!(s_cart.len(), n_cart * n_cart);
        assert_eq!(s_sph.len(), n_sph * n_sph);
    }

    #[test]
    fn test_overlap_spherical_no_change_for_sp_only() {
        // For STO-3G (s and p orbitals only), spherical and Cartesian are identical
        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        // STO-3G has no d-orbitals
        assert!(!basis.has_spherical_difference());
        assert_eq!(basis.n_basis, basis.n_basis_spherical());

        let s_cart = overlap_matrix(&basis);
        let s_sph = overlap_matrix_spherical(&basis);

        // Matrices should be identical
        assert_eq!(s_cart.len(), s_sph.len());
        for (_i, (cart, sph)) in s_cart.iter().zip(s_sph.iter()).enumerate() {
            assert_abs_diff_eq!(cart, sph, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_overlap_spherical_symmetry() {
        // Spherical overlap matrix should remain symmetric
        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();

        let s_sph = overlap_matrix_spherical(&basis);
        let n = basis.n_basis_spherical();

        for i in 0..n {
            for j in 0..n {
                let diff = (s_sph[i * n + j] - s_sph[j * n + i]).abs();
                assert!(
                    diff < 1e-14,
                    "Spherical overlap not symmetric at ({}, {}): {} vs {}",
                    i,
                    j,
                    s_sph[i * n + j],
                    s_sph[j * n + i]
                );
            }
        }
    }

    #[test]
    fn test_overlap_spherical_positive_diagonal() {
        // Spherical overlap matrix diagonal elements should be positive
        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();

        let s_sph = overlap_matrix_spherical(&basis);
        let n = basis.n_basis_spherical();

        for i in 0..n {
            let diag = s_sph[i * n + i];
            assert!(
                diag > 0.0,
                "Spherical overlap diagonal [{},{}] should be positive, got {}",
                i,
                i,
                diag
            );
        }
    }

    // -------------------------------------------------------------------------
    // evaluate_overlap_vs_distance tests (US-054)
    // -------------------------------------------------------------------------
    //
    // Golden values generated from PySCF 2.11.0:
    //   mol.atom = f'H 0 0 0; H 0 0 {R}'
    //   mol.unit = 'bohr'; mol.basis = 'sto-3g'; mol.build(verbose=0)
    //   S = mol.intor('int1e_ovlp'); S[0,1]

    /// Helper: construct H(1s) STO-3G contracted shell at a given center
    fn h_1s_sto3g_shell(center: [f64; 3]) -> ContractedShell {
        ContractedShell::new(AngularMomentum::S, h_sto3g_primitives(), center, 0)
    }

    /// Helper: construct He(1s) STO-3G contracted shell at a given center
    fn he_1s_sto3g_shell(center: [f64; 3]) -> ContractedShell {
        ContractedShell::new(
            AngularMomentum::S,
            vec![
                GaussianPrimitive::new(6.3624213900, 0.1543289704),
                GaussianPrimitive::new(1.1589230000, 0.5353281416),
                GaussianPrimitive::new(0.3136497900, 0.4446345413),
            ],
            center,
            0,
        )
    }

    // OD-R1 through OD-R8: H(1s)-H(1s) STO-3G overlap at specific distances
    // Reference: PySCF 2.11.0 mol.intor('int1e_ovlp')[0,1]
    //
    // Note: IQCP's Obara-Saika implementation matches PySCF within ~1e-9.
    // The small discrepancy arises from accumulated floating-point rounding
    // in the recurrence relations and normalization conventions. This is
    // consistent with the existing overlap module tolerance (TOL = 1e-9).
    #[test]
    fn test_overlap_vs_distance_h_h_sto3g_golden() {
        let shell_a = h_1s_sto3g_shell([0.0, 0.0, 0.0]);
        let shell_b = h_1s_sto3g_shell([0.0, 0.0, 0.0]);

        let r_values = vec![0.5, 1.0, 1.3984, 2.0, 3.0, 5.0, 8.0, 10.0];
        let expected = vec![
            9.405268653093932e-01, // R=0.5
            7.965883006970911e-01, // R=1.0
            6.598721980070731e-01, // R=1.3984
            4.627776954301663e-01, // R=2.0
            2.261896476964506e-01, // R=3.0
            3.747005011027010e-02, // R=5.0
            9.626096800740930e-04, // R=8.0
            4.319600325118806e-05, // R=10.0
        ];

        let result = evaluate_overlap_vs_distance(&shell_a, &shell_b, &r_values);

        assert_eq!(result.len(), expected.len());
        for (computed, reference) in result.iter().zip(expected.iter()) {
            // Use relative tolerance for values spanning many orders of magnitude.
            // For small absolute values (< 1e-10), use proportional absolute tolerance.
            let r: f64 = *reference;
            let c: f64 = *computed;
            let tol = if r.abs() < 1e-10 {
                1e-14
            } else {
                r.abs() * 1e-9
            };
            assert!(
                (c - r).abs() < tol,
                "Overlap mismatch: computed={:.15e}, reference={:.15e}, diff={:.2e}, tol={:.2e}",
                c,
                r,
                (c - r).abs(),
                tol,
            );
        }
    }

    // OD-R9: Self-overlap at R=0
    // The diagonal overlap S(R=0) for identical normalized shells should be
    // very close to 1.0. We use the same tolerance as the overlap matrix
    // diagonal tests (TOL = 1e-9).
    #[test]
    fn test_overlap_vs_distance_self_overlap_r0() {
        let shell = h_1s_sto3g_shell([0.0, 0.0, 0.0]);
        let result = evaluate_overlap_vs_distance(&shell, &shell, &[0.0]);

        assert_eq!(result.len(), 1);
        assert_abs_diff_eq!(result[0], 1.0, epsilon = TOL);
    }

    // OD-R10: Symmetry: S(H,He) = S(He,H)
    #[test]
    fn test_overlap_vs_distance_h_he_symmetry() {
        let h_shell = h_1s_sto3g_shell([0.0, 0.0, 0.0]);
        let he_shell = he_1s_sto3g_shell([0.0, 0.0, 0.0]);

        let r = [2.0];
        let s_h_he = evaluate_overlap_vs_distance(&h_shell, &he_shell, &r);
        let s_he_h = evaluate_overlap_vs_distance(&he_shell, &h_shell, &r);

        assert_abs_diff_eq!(s_h_he[0], s_he_h[0], epsilon = 1e-14);
    }

    // OD-R11: Pairwise symmetry for shell_a, shell_b swap
    #[test]
    fn test_overlap_vs_distance_pairwise_symmetry() {
        let h_shell = h_1s_sto3g_shell([0.0, 0.0, 0.0]);
        let he_shell = he_1s_sto3g_shell([0.0, 0.0, 0.0]);

        let r_values: Vec<f64> = (1..=10).map(|i| i as f64 * 0.5).collect();
        let s_ab = evaluate_overlap_vs_distance(&h_shell, &he_shell, &r_values);
        let s_ba = evaluate_overlap_vs_distance(&he_shell, &h_shell, &r_values);

        for (a, b) in s_ab.iter().zip(s_ba.iter()) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-14);
        }
    }

    // OD-R12 through OD-R14: H(1s)-He(1s) STO-3G at specific distances
    // Reference: PySCF 2.11.0 (mol.spin=1 for odd electron count)
    #[test]
    fn test_overlap_vs_distance_h_he_sto3g_golden() {
        let h_shell = h_1s_sto3g_shell([0.0, 0.0, 0.0]);
        let he_shell = he_1s_sto3g_shell([0.0, 0.0, 0.0]);

        let r_values = vec![1.0, 2.0, 5.0];
        let expected = vec![
            7.142217424910036e-01, // R=1.0
            3.590320602452557e-01, // R=2.0
            1.651056913935484e-02, // R=5.0
        ];

        let result = evaluate_overlap_vs_distance(&h_shell, &he_shell, &r_values);

        assert_eq!(result.len(), expected.len());
        for (computed, reference) in result.iter().zip(expected.iter()) {
            assert_abs_diff_eq!(computed, reference, epsilon = 1e-10);
        }
    }

    // OD-R15: Monotonic decrease for s-s overlap as R increases
    #[test]
    fn test_overlap_vs_distance_monotonic_decrease() {
        let shell = h_1s_sto3g_shell([0.0, 0.0, 0.0]);

        let r_values: Vec<f64> = (1..=100).map(|i| i as f64 * 0.1).collect();
        let overlaps = evaluate_overlap_vs_distance(&shell, &shell, &r_values);

        for i in 1..overlaps.len() {
            assert!(
                overlaps[i] < overlaps[i - 1],
                "Overlap not monotonically decreasing at R={:.1}: S[{}]={} >= S[{}]={}",
                r_values[i],
                i,
                overlaps[i],
                i - 1,
                overlaps[i - 1],
            );
        }
    }

    // OD-R16: Non-negativity for s-s overlap
    #[test]
    fn test_overlap_vs_distance_non_negative() {
        let shell = h_1s_sto3g_shell([0.0, 0.0, 0.0]);

        let r_values: Vec<f64> = (0..=200).map(|i| i as f64 * 0.1).collect();
        let overlaps = evaluate_overlap_vs_distance(&shell, &shell, &r_values);

        for (i, &s) in overlaps.iter().enumerate() {
            assert!(
                s >= 0.0,
                "Overlap should be non-negative at R={:.1}, got {}",
                r_values[i],
                s,
            );
        }
    }

    // OD-R17: Output length equals input length
    #[test]
    fn test_overlap_vs_distance_output_length() {
        let shell = h_1s_sto3g_shell([0.0, 0.0, 0.0]);

        let r_values: Vec<f64> = (0..50).map(|i| i as f64 * 0.2).collect();
        let result = evaluate_overlap_vs_distance(&shell, &shell, &r_values);

        assert_eq!(result.len(), r_values.len());
    }

    // OD-R18: Empty r_values returns empty vec
    #[test]
    fn test_overlap_vs_distance_empty_input() {
        let shell = h_1s_sto3g_shell([0.0, 0.0, 0.0]);
        let result = evaluate_overlap_vs_distance(&shell, &shell, &[]);

        assert!(result.is_empty());
    }
}
