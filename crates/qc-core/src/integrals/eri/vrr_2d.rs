//! 2D Vertical Recurrence Relations (VRR) for ERI computation
//!
//! Builds the 2D integral table g[n][m] for a single Cartesian direction
//! using the vertical recurrence relations.
//!
//! # 2D Recurrence Relations
//!
//! The recurrence builds integrals in a 2D grid indexed by (n, m) where:
//! - n = total angular momentum on bra side (A + B)
//! - m = total angular momentum on ket side (C + D)
//!
//! **Bra recurrence (build n from n-1):**
//! ```text
//! [n+1, m] = c00 * [n, m] + n * b10 * [n-1, m] + m * b00 * [n, m-1]
//! ```
//!
//! **Ket recurrence (build m from m-1):**
//! ```text
//! [n, m+1] = c0p * [n, m] + m * b01 * [n, m-1] + n * b00 * [n-1, m]
//! ```
//!
//! # Base Case
//!
//! The base case [0, 0] = 1.0 (prefactors are handled separately).
//!
//! # Reference
//!
//! - libcint g2e.c lines 272-421 (CINTg0_2e_2d)
//! - Dupuis, Rys & King (1976), J. Chem. Phys. 65, 111

/// Build the 2D VRR integral table for one Cartesian direction
///
/// # Arguments
///
/// * `n_bra` - Maximum angular momentum needed on bra side (i + j)
/// * `n_ket` - Maximum angular momentum needed on ket side (k + l)
/// * `c00` - Center offset coefficient for bra recurrence
/// * `c0p` - Center offset coefficient for ket recurrence
/// * `b00` - Cross-term coupling coefficient
/// * `b10` - Bra recurrence coefficient
/// * `b01` - Ket recurrence coefficient
///
/// # Returns
///
/// A 2D array (flattened row-major) of size (n_bra + 1) x (n_ket + 1).
/// Access: g[n * (n_ket + 1) + m] for g[n][m].
///
/// # Algorithm
///
/// We use a mixed recurrence strategy:
/// 1. Build column [n][0] using pure bra recurrence
/// 2. Build each row [n][m] using ket recurrence with cross-term
///
/// This approach minimizes numerical error accumulation.
///
/// # Reference
///
/// libcint g2e.c lines 272-421 (CINTg0_2e_2d)
pub fn build_2d(
    n_bra: usize,
    n_ket: usize,
    c00: f64,
    c0p: f64,
    b00: f64,
    b10: f64,
    b01: f64,
) -> Vec<f64> {
    let rows = n_bra + 1;
    let cols = n_ket + 1;

    // Allocate 2D array (row-major: g[n][m] at index n * cols + m)
    let mut g = vec![0.0; rows * cols];

    // Base case: [0, 0] = 1.0
    g[0] = 1.0;

    // Handle the case where no angular momentum is needed
    if n_bra == 0 && n_ket == 0 {
        return g;
    }

    // First, build the first column (m = 0) using bra recurrence
    // [n+1, 0] = c00 * [n, 0] + n * b10 * [n-1, 0]
    if n_bra > 0 {
        // [1, 0] = c00 * [0, 0] = c00
        g[cols] = c00;

        for n in 1..n_bra {
            // [n+1, 0] = c00 * [n, 0] + n * b10 * [n-1, 0]
            g[(n + 1) * cols] = c00 * g[n * cols] + (n as f64) * b10 * g[(n - 1) * cols];
        }
    }

    // Now build rows (m > 0) using ket recurrence
    // [n, m+1] = c0p * [n, m] + m * b01 * [n, m-1] + n * b00 * [n-1, m]
    for m in 0..n_ket {
        // First entry of the new column: n = 0
        // [0, m+1] = c0p * [0, m] + m * b01 * [0, m-1]
        // (no cross-term since n = 0)
        if m == 0 {
            g[m + 1] = c0p * g[0];
        } else {
            g[m + 1] = c0p * g[m] + (m as f64) * b01 * g[m - 1];
        }

        // Build remaining entries: n > 0
        // [n, m+1] = c0p * [n, m] + m * b01 * [n, m-1] + n * b00 * [n-1, m]
        for n in 1..=n_bra {
            let g_n_m = g[n * cols + m];
            let g_n_m_minus_1 = if m > 0 { g[n * cols + m - 1] } else { 0.0 };
            let g_n_minus_1_m = g[(n - 1) * cols + m];

            g[n * cols + m + 1] =
                c0p * g_n_m + (m as f64) * b01 * g_n_m_minus_1 + (n as f64) * b00 * g_n_minus_1_m;
        }
    }

    g
}

/// Build the 2D VRR integral table into a pre-allocated buffer.
///
/// Same algorithm as `build_2d` but avoids heap allocation by writing
/// into the caller's buffer. The buffer must have at least
/// `(n_bra + 1) * (n_ket + 1)` elements; only those elements are written.
///
/// # Arguments
///
/// * `g` - Pre-allocated output buffer (must be large enough)
/// * `n_bra` - Maximum angular momentum on bra side
/// * `n_ket` - Maximum angular momentum on ket side
/// * `c00` - Center offset coefficient for bra recurrence
/// * `c0p` - Center offset coefficient for ket recurrence
/// * `b00` - Cross-term coupling coefficient
/// * `b10` - Bra recurrence coefficient
/// * `b01` - Ket recurrence coefficient
///
/// # Returns
///
/// The number of elements written: `(n_bra + 1) * (n_ket + 1)`
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn build_2d_into(
    g: &mut [f64],
    n_bra: usize,
    n_ket: usize,
    c00: f64,
    c0p: f64,
    b00: f64,
    b10: f64,
    b01: f64,
) -> usize {
    let rows = n_bra + 1;
    let cols = n_ket + 1;
    let size = rows * cols;

    // Zero out the used portion
    for v in &mut g[..size] {
        *v = 0.0;
    }

    // Base case: [0, 0] = 1.0
    g[0] = 1.0;

    if n_bra == 0 && n_ket == 0 {
        return size;
    }

    // Build the first column (m = 0) using bra recurrence
    if n_bra > 0 {
        g[cols] = c00;
        for n in 1..n_bra {
            g[(n + 1) * cols] = c00 * g[n * cols] + (n as f64) * b10 * g[(n - 1) * cols];
        }
    }

    // Build rows (m > 0) using ket recurrence
    for m in 0..n_ket {
        if m == 0 {
            g[m + 1] = c0p * g[0];
        } else {
            g[m + 1] = c0p * g[m] + (m as f64) * b01 * g[m - 1];
        }

        for n in 1..=n_bra {
            let g_n_m = g[n * cols + m];
            let g_n_m_minus_1 = if m > 0 { g[n * cols + m - 1] } else { 0.0 };
            let g_n_minus_1_m = g[(n - 1) * cols + m];

            g[n * cols + m + 1] =
                c0p * g_n_m + (m as f64) * b01 * g_n_m_minus_1 + (n as f64) * b00 * g_n_minus_1_m;
        }
    }

    size
}

/// Get value from 2D VRR table
///
/// # Arguments
///
/// * `g` - The 2D VRR table (from build_2d)
/// * `n_ket` - Number of columns - 1 (to compute stride)
/// * `n` - Bra angular momentum index
/// * `m` - Ket angular momentum index
///
/// # Returns
///
/// The value g[n][m]
#[inline]
pub fn get_2d(g: &[f64], n_ket: usize, n: usize, m: usize) -> f64 {
    g[n * (n_ket + 1) + m]
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    #[test]
    fn test_base_case() {
        // When n_bra = n_ket = 0, just get [0,0] = 1
        let g = build_2d(0, 0, 0.3, 0.4, 0.1, 0.2, 0.25);
        assert_eq!(g.len(), 1);
        assert!((g[0] - 1.0).abs() < TOL);
    }

    #[test]
    fn test_first_bra_recurrence() {
        // [1, 0] = c00 * [0, 0] = c00
        let c00 = 0.5;
        let g = build_2d(1, 0, c00, 0.4, 0.1, 0.2, 0.25);

        assert_eq!(g.len(), 2); // 2 rows, 1 column
        assert!((g[0] - 1.0).abs() < TOL); // [0,0]
        assert!((g[1] - c00).abs() < TOL); // [1,0]
    }

    #[test]
    fn test_first_ket_recurrence() {
        // [0, 1] = c0p * [0, 0] = c0p
        let c0p = 0.6;
        let g = build_2d(0, 1, 0.5, c0p, 0.1, 0.2, 0.25);

        assert_eq!(g.len(), 2); // 1 row, 2 columns
        assert!((g[0] - 1.0).abs() < TOL); // [0,0]
        assert!((g[1] - c0p).abs() < TOL); // [0,1]
    }

    #[test]
    fn test_second_bra_recurrence() {
        // [2, 0] = c00 * [1, 0] + 1 * b10 * [0, 0]
        //        = c00 * c00 + b10
        let c00 = 0.5;
        let b10 = 0.3;
        let g = build_2d(2, 0, c00, 0.4, 0.1, b10, 0.25);

        let expected = c00 * c00 + b10;
        assert!((g[2] - expected).abs() < TOL); // [2,0]
    }

    #[test]
    fn test_second_ket_recurrence() {
        // [0, 2] = c0p * [0, 1] + 1 * b01 * [0, 0]
        //        = c0p * c0p + b01
        let c0p = 0.6;
        let b01 = 0.25;
        let g = build_2d(0, 2, 0.5, c0p, 0.1, 0.3, b01);

        let expected = c0p * c0p + b01;
        assert!((g[2] - expected).abs() < TOL); // [0,2]
    }

    #[test]
    fn test_cross_term() {
        // [1, 1] = c0p * [1, 0] + 0 * b01 * [1, -1] + 1 * b00 * [0, 0]
        //        = c0p * c00 + b00
        let c00 = 0.5;
        let c0p = 0.6;
        let b00 = 0.1;

        let g = build_2d(1, 1, c00, c0p, b00, 0.3, 0.25);

        // [1,0] = c00
        // [1,1] = c0p * c00 + b00
        let expected = c0p * c00 + b00;
        let cols = 2; // n_ket + 1
        assert!((g[cols + 1] - expected).abs() < TOL);
    }

    #[test]
    fn test_full_grid() {
        // Build a 3x3 grid and verify consistency
        let c00 = 0.3;
        let c0p = 0.4;
        let b00 = 0.1;
        let b10 = 0.2;
        let b01 = 0.15;

        let g = build_2d(2, 2, c00, c0p, b00, b10, b01);
        let cols = 3;

        // Check [0,0] = 1
        assert!((g[0] - 1.0).abs() < TOL);

        // Check [1,0] = c00
        assert!((g[cols] - c00).abs() < TOL);

        // Check [0,1] = c0p
        assert!((g[1] - c0p).abs() < TOL);

        // Check [2,0] = c00 * [1,0] + 1 * b10 * [0,0]
        let g_2_0 = c00 * c00 + b10;
        assert!((g[2 * cols] - g_2_0).abs() < TOL);

        // Check [0,2] = c0p * [0,1] + 1 * b01 * [0,0]
        let g_0_2 = c0p * c0p + b01;
        assert!((g[2] - g_0_2).abs() < TOL);

        // Check [1,1] = c0p * [1,0] + 0 + 1 * b00 * [0,0]
        let g_1_1 = c0p * c00 + b00;
        assert!((g[cols + 1] - g_1_1).abs() < TOL);

        // Check [2,1] = c0p * [2,0] + 0 + 2 * b00 * [1,0]
        let g_2_1 = c0p * g_2_0 + 2.0 * b00 * c00;
        assert!((g[2 * cols + 1] - g_2_1).abs() < TOL);

        // Check [1,2] = c0p * [1,1] + 1 * b01 * [1,0] + 1 * b00 * [0,1]
        let g_1_2 = c0p * g_1_1 + b01 * c00 + b00 * c0p;
        assert!((g[cols + 2] - g_1_2).abs() < TOL);
    }

    #[test]
    fn test_get_2d() {
        let c00 = 0.5;
        let c0p = 0.6;
        let b00 = 0.1;
        let b10 = 0.2;
        let b01 = 0.25;

        let g = build_2d(2, 2, c00, c0p, b00, b10, b01);

        // Verify get_2d gives same values as direct indexing
        let n_ket = 2;
        for n in 0..=2 {
            for m in 0..=2 {
                let expected = g[n * 3 + m];
                let actual = get_2d(&g, n_ket, n, m);
                assert!(
                    (actual - expected).abs() < TOL,
                    "Mismatch at ({}, {}): {} vs {}",
                    n,
                    m,
                    actual,
                    expected
                );
            }
        }
    }

    #[test]
    fn test_zero_coefficients() {
        // With all coefficients = 0, only [0,0] = 1
        let g = build_2d(2, 2, 0.0, 0.0, 0.0, 0.0, 0.0);

        assert!((g[0] - 1.0).abs() < TOL); // [0,0] = 1

        // All others should be 0
        for (i, &val) in g.iter().enumerate().skip(1) {
            assert!(val.abs() < TOL, "g[{}] = {} should be 0", i, val);
        }
    }
}
