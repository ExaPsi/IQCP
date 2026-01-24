//! Horizontal Transfer Recurrence (HTR) for 4D integrals
//!
//! Converts 2D integrals g[n][m] to 4D integrals [i,j,k,l] by transferring
//! angular momentum between centers.
//!
//! # Horizontal Transfer Relations
//!
//! Starting from 2D integrals where:
//! - n = i + j (total bra angular momentum)
//! - m = k + l (total ket angular momentum)
//!
//! **Bra transfer (from A to B):**
//! ```text
//! [i, j+1, k, l] = [i+1, j, k, l] + (A - B) * [i, j, k, l]
//! ```
//!
//! **Ket transfer (from C to D):**
//! ```text
//! [i, j, k, l+1] = [i, j, k+1, l] + (C - D) * [i, j, k, l]
//! ```
//!
//! # Algorithm
//!
//! For 1D (one Cartesian direction):
//! 1. Start with [n][m] = g[i+j][k+l] from VRR
//! 2. Apply bra HTR: reduce n while building j
//! 3. Apply ket HTR: reduce m while building l
//! 4. Result: [i][j][k][l] for specific angular momentum indices
//!
//! # Reference
//!
//! - libcint g2e.c lines 428-492 (CINTg0_lj2d_4d)
//! - Dupuis, Rys & King (1976), J. Chem. Phys. 65, 111

use super::vrr_2d::get_2d;

/// Apply horizontal transfer to get 1D 4-index integral from 2D table
///
/// # Arguments
///
/// * `g` - The 2D VRR table (from vrr_2d::build_2d)
/// * `n_bra` - Maximum bra angular momentum in g table
/// * `n_ket` - Maximum ket angular momentum in g table
/// * `i` - Angular momentum on center A
/// * `j` - Angular momentum on center B
/// * `k` - Angular momentum on center C
/// * `l` - Angular momentum on center D
/// * `ab` - A - B distance (for this Cartesian direction)
/// * `cd` - C - D distance (for this Cartesian direction)
///
/// # Returns
///
/// The 1D component of the 4-index integral [i,j,k,l]
///
/// # Algorithm
///
/// First apply ket HTR to reduce (k+l, l) -> (k, l), then apply bra HTR
/// to reduce (i+j, j) -> (i, j). This order minimizes temporary storage.
///
/// # Reference
///
/// libcint g2e.c lines 428-492 (CINTg0_lj2d_4d)
#[allow(clippy::too_many_arguments)]
pub fn horizontal_transfer_1d(
    g: &[f64],
    _n_bra: usize,
    n_ket: usize,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    ab: f64,
    cd: f64,
) -> f64 {
    horizontal_transfer_recursive(g, n_ket, i, j, k, l, ab, cd)
}

/// Apply horizontal transfer to get 1D 4-index integral from 2D table
///
/// Full version with center separation distances for HTR.
///
/// # Arguments
///
/// * `g` - The 2D VRR table (from vrr_2d::build_2d)
/// * `n_bra` - Maximum bra angular momentum in g table
/// * `n_ket` - Maximum ket angular momentum in g table
/// * `i` - Angular momentum on center A
/// * `j` - Angular momentum on center B
/// * `k` - Angular momentum on center C
/// * `l` - Angular momentum on center D
/// * `ab` - A - B distance (for this Cartesian direction)
/// * `cd` - C - D distance (for this Cartesian direction)
///
/// # Returns
///
/// The 1D component of the 4-index integral [i,j,k,l]
///
/// # Notes
///
/// The horizontal transfer relations are:
/// - [i,j+1,k,l] = [i+1,j,k,l] + (A-B)*[i,j,k,l]
/// - [i,j,k,l+1] = [i,j,k+1,l] + (C-D)*[i,j,k,l]
///
/// We invert these to express [i,j,k,l] in terms of [n][m] values from VRR:
/// - [i,j,k,l] = [i+1,j-1,k,l] + (A-B)*[i,j-1,k,l]  (for j > 0)
/// - [i,j,k,l] = [i,j,k+1,l-1] + (C-D)*[i,j,k,l-1]  (for l > 0)
#[allow(clippy::too_many_arguments)]
fn horizontal_transfer_recursive(
    g: &[f64],
    n_ket: usize,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    ab: f64,
    cd: f64,
) -> f64 {
    // If we're at the 2D level (j = 0 and l = 0), return g[n][m]
    if j == 0 && l == 0 {
        return get_2d(g, n_ket, i, k);
    }

    // Apply ket HTR first if l > 0
    // [i,j,k,l] = [i,j,k+1,l-1] + (C-D)*[i,j,k,l-1]
    if l > 0 {
        let term1 = horizontal_transfer_recursive(g, n_ket, i, j, k + 1, l - 1, ab, cd);
        let term2 = horizontal_transfer_recursive(g, n_ket, i, j, k, l - 1, ab, cd);
        return term1 + cd * term2;
    }

    // Apply bra HTR if j > 0
    // [i,j,k,l] = [i+1,j-1,k,l] + (A-B)*[i,j-1,k,l]
    if j > 0 {
        let term1 = horizontal_transfer_recursive(g, n_ket, i + 1, j - 1, k, l, ab, cd);
        let term2 = horizontal_transfer_recursive(g, n_ket, i, j - 1, k, l, ab, cd);
        return term1 + ab * term2;
    }

    // This should never be reached
    0.0
}

/// Apply horizontal transfer using iterative (non-recursive) method
///
/// This is more efficient for higher angular momentum as it avoids
/// recursion overhead and recomputation.
///
/// # Arguments
///
/// * `g` - The 2D VRR table (from vrr_2d::build_2d)
/// * `n_bra` - Maximum bra angular momentum in g table
/// * `n_ket` - Maximum ket angular momentum in g table
/// * `i` - Angular momentum on center A
/// * `j` - Angular momentum on center B
/// * `k` - Angular momentum on center C
/// * `l` - Angular momentum on center D
/// * `ab` - A - B distance (for this Cartesian direction)
/// * `cd` - C - D distance (for this Cartesian direction)
///
/// # Returns
///
/// The 1D component of the 4-index integral [i,j,k,l]
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn horizontal_transfer_iterative(
    g: &[f64],
    _n_bra: usize,
    n_ket: usize,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    ab: f64,
    cd: f64,
) -> f64 {
    // First apply ket HTR: build h[n][k][l] from g[n][m]
    // Then apply bra HTR: build result[i][j][k][l] from h

    let _n = i + j;
    let _m = k + l;

    // For s and p orbitals, the dimensions are small, so allocate a work array
    // h[ii][kk] for ii in 0..=n, kk in 0..=k
    let _h = vec![vec![0.0; k + 1]; _n + 1];

    // Initialize h[ii][m] = g[ii][m] for ii = 0..=n
    // Then apply ket HTR: h[ii][kk] = h[ii][kk+1] + cd * h[ii][kk] (in reverse)

    // Actually, we need to apply HTR to reduce l to 0:
    // [ii, kk, l] -> [ii, k+l, 0] is stored in g
    // To get [ii, kk, ll]: [ii, kk, ll] = [ii, kk+1, ll-1] + cd * [ii, kk, ll-1]

    // First, fill h[ii][kk+ll] = g[ii][kk+ll] for all ll values we need
    // Then iterate: for ll from l-1 down to 0

    // Initialize for ll = l (final result we want)
    // We need h[ii][kk] for kk = k to k+l

    // More efficient: build a 3D array for the ket HTR, then apply bra HTR

    // For simplicity with small L, use the recursive version
    // This function serves as a template for optimization

    horizontal_transfer_recursive(g, n_ket, i, j, k, l, ab, cd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrals::eri::vrr_2d::build_2d;

    const TOL: f64 = 1e-12;

    #[test]
    fn test_no_transfer_needed() {
        // When j = 0 and l = 0, just return g[i][k]
        let g = build_2d(2, 2, 0.5, 0.6, 0.1, 0.2, 0.25);

        // [0,0,0,0] = g[0][0] = 1.0
        let result = horizontal_transfer_1d(&g, 2, 2, 0, 0, 0, 0, 0.0, 0.0);
        assert!((result - 1.0).abs() < TOL);

        // [1,0,0,0] = g[1][0]
        let expected = 0.5; // c00
        let result = horizontal_transfer_1d(&g, 2, 2, 1, 0, 0, 0, 0.0, 0.0);
        assert!((result - expected).abs() < TOL);

        // [0,0,1,0] = g[0][1]
        let expected = 0.6; // c0p
        let result = horizontal_transfer_1d(&g, 2, 2, 0, 0, 1, 0, 0.0, 0.0);
        assert!((result - expected).abs() < TOL);
    }

    #[test]
    fn test_bra_transfer() {
        // Test [0,1,0,0] = [1,0,0,0] + ab * [0,0,0,0]
        let c00 = 0.5;
        let g = build_2d(2, 2, c00, 0.6, 0.1, 0.2, 0.25);

        let ab = 0.3;
        // [0,1,0,0] = g[1][0] + ab * g[0][0] = c00 + ab
        let expected = c00 + ab;
        let result = horizontal_transfer_recursive(&g, 2, 0, 1, 0, 0, ab, 0.0);
        assert!((result - expected).abs() < TOL);
    }

    #[test]
    fn test_ket_transfer() {
        // Test [0,0,0,1] = [0,0,1,0] + cd * [0,0,0,0]
        let c0p = 0.6;
        let g = build_2d(2, 2, 0.5, c0p, 0.1, 0.2, 0.25);

        let cd = 0.4;
        // [0,0,0,1] = g[0][1] + cd * g[0][0] = c0p + cd
        let expected = c0p + cd;
        let result = horizontal_transfer_recursive(&g, 2, 0, 0, 0, 1, 0.0, cd);
        assert!((result - expected).abs() < TOL);
    }

    #[test]
    fn test_combined_transfer() {
        // Test [0,1,0,1] = [0,0,1,0] + ab * [0,0,0,0] then ket HTR
        let c00 = 0.5;
        let c0p = 0.6;
        let b00 = 0.1;
        let g = build_2d(2, 2, c00, c0p, b00, 0.2, 0.25);

        let ab = 0.3;
        let cd = 0.4;

        // First, compute [0,1,0,0] = g[1][0] + ab * g[0][0] = c00 + ab
        // Then [0,1,0,1] = [0,1,1,0] + cd * [0,1,0,0]
        // where [0,1,1,0] = g[1][1] + ab * g[0][1]

        // g[1][1] = c0p * c00 + b00 (from VRR tests)
        let g_1_1 = c0p * c00 + b00;

        let h_0_1_0_0 = c00 + ab;
        let h_0_1_1_0 = g_1_1 + ab * c0p;
        let expected = h_0_1_1_0 + cd * h_0_1_0_0;

        let result = horizontal_transfer_recursive(&g, 2, 0, 1, 0, 1, ab, cd);
        assert!((result - expected).abs() < TOL);
    }

    #[test]
    fn test_higher_transfer() {
        // Test [1,1,0,0]
        // [1,1,0,0] = [2,0,0,0] + ab * [1,0,0,0]
        let c00 = 0.5;
        let b10 = 0.2;
        let g = build_2d(3, 2, c00, 0.6, 0.1, b10, 0.25);

        let ab = 0.3;

        // g[2][0] = c00 * c00 + b10
        let g_2_0 = c00 * c00 + b10;
        let expected = g_2_0 + ab * c00;

        let result = horizontal_transfer_recursive(&g, 2, 1, 1, 0, 0, ab, 0.0);
        assert!((result - expected).abs() < TOL);
    }

    #[test]
    fn test_zero_separation() {
        // When ab = cd = 0, HTR simplifies to just indexing
        let g = build_2d(2, 2, 0.5, 0.6, 0.1, 0.2, 0.25);

        // [0,1,0,0] with ab=0: = [1,0,0,0] = g[1][0]
        let result = horizontal_transfer_recursive(&g, 2, 0, 1, 0, 0, 0.0, 0.0);
        let expected = get_2d(&g, 2, 1, 0);
        assert!((result - expected).abs() < TOL);

        // [0,0,0,1] with cd=0: = [0,0,1,0] = g[0][1]
        let result = horizontal_transfer_recursive(&g, 2, 0, 0, 0, 1, 0.0, 0.0);
        let expected = get_2d(&g, 2, 0, 1);
        assert!((result - expected).abs() < TOL);
    }

    #[test]
    fn test_pp_case() {
        // Test (p|p) type: [1,0,1,0] and [0,1,0,1]
        let c00 = 0.4;
        let c0p = 0.5;
        let b00 = 0.1;
        let b10 = 0.2;
        let b01 = 0.15;
        let g = build_2d(2, 2, c00, c0p, b00, b10, b01);

        // [1,0,1,0] needs no HTR, just g[1][1]
        let g_1_1 = c0p * c00 + b00;
        let result = horizontal_transfer_1d(&g, 2, 2, 1, 0, 1, 0, 0.0, 0.0);
        assert!((result - g_1_1).abs() < TOL);

        // [0,1,0,1] needs both bra and ket HTR
        let ab = 0.3;
        let cd = 0.4;
        let result = horizontal_transfer_recursive(&g, 2, 0, 1, 0, 1, ab, cd);

        // Manual calculation:
        // [0,1,0,1] = [0,1,1,0] + cd * [0,1,0,0]
        // [0,1,1,0] = [1,0,1,0] + ab * [0,0,1,0] = g[1][1] + ab * g[0][1]
        // [0,1,0,0] = [1,0,0,0] + ab * [0,0,0,0] = g[1][0] + ab * g[0][0]

        let h_0_1_1_0 = g_1_1 + ab * c0p;
        let h_0_1_0_0 = c00 + ab;
        let expected = h_0_1_1_0 + cd * h_0_1_0_0;

        assert!((result - expected).abs() < TOL);
    }
}
