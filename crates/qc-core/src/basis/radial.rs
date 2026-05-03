//! Radial profile evaluation for basis function visualization.
//!
//! Evaluates the radial part of contracted Gaussian basis functions at
//! arbitrary r values for Module D's radial profile plot.
//!
//! # Radial Function
//!
//! The radial profile of a contracted shell is:
//!
//! ```text
//! R(r) = r^l * sum_k { d_k * exp(-alpha_k * r^2) }
//! ```
//!
//! where:
//! - `l` is the angular momentum quantum number (0=s, 1=p, 2=d)
//! - `alpha_k` are the primitive Gaussian exponents
//! - `d_k = c_k * N_k` is the effective coefficient (raw contraction coefficient
//!   times primitive normalization)
//!
//! # Normalization Convention
//!
//! Uses the same normalization as `cartesian_gaussian_normalization` in
//! `integrals/overlap.rs` with a representative Cartesian component for
//! each angular momentum:
//! - s (l=0): powers (0,0,0)
//! - p (l=1): powers (0,0,1)
//! - d (l=2): powers (1,1,0)
//!
//! This ensures visual profiles are consistent with the SCF integral engine.
//!
//! # References
//!
//! - Szabo & Ostlund (1996), "Modern Quantum Chemistry", Section 3.6, Eq. 3.203-3.210

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Result of evaluating a radial profile for a contracted shell.
///
/// Contains the contracted profile, per-primitive profiles, and metadata
/// needed for plotting and educational tooltips.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RadialProfileResult {
    /// r values in Bohr at which the profile was evaluated
    pub r_values: Vec<f64>,

    /// Contracted radial profile R(r) = r^l * sum_k d_k * exp(-alpha_k * r^2)
    pub contracted_values: Vec<f64>,

    /// Individual primitive profiles R_k(r) = r^l * d_k * exp(-alpha_k * r^2)
    ///
    /// Outer Vec: one per primitive. Inner Vec: values at each r point.
    pub primitive_values: Vec<Vec<f64>>,

    /// Exponents of each primitive (for tooltips)
    pub exponents: Vec<f64>,

    /// Effective coefficients d_k = c_k * N_k (for tooltips)
    pub effective_coefficients: Vec<f64>,

    /// Raw contraction coefficients c_k (for display)
    pub raw_coefficients: Vec<f64>,

    /// Angular momentum quantum number l (0=s, 1=p, 2=d)
    pub angular_momentum: u32,

    /// Angular momentum label ("s", "p", "d")
    pub angular_momentum_label: String,

    /// Number of primitives
    pub n_primitives: usize,

    /// Auto-determined r_max in Bohr
    pub r_max: f64,
}

/// Compute the normalization factor for the radial part of a Gaussian primitive.
///
/// Uses a representative Cartesian component for each angular momentum to
/// match the integral engine normalization:
/// - s (l=0): (0,0,0) -> N = (2*alpha/pi)^(3/4)
/// - p (l=1): (0,0,1) -> N = (2*alpha/pi)^(3/4) * (4*alpha)^(1/2)
/// - d (l=2): (1,1,0) -> N = (2*alpha/pi)^(3/4) * (4*alpha)
///
/// # Arguments
///
/// * `alpha` - Primitive Gaussian exponent
/// * `l` - Angular momentum quantum number (0, 1, or 2)
///
/// # References
///
/// Szabo & Ostlund Eq. 3.203. Consistent with
/// `crate::integrals::overlap::cartesian_gaussian_normalization`.
#[inline]
pub fn primitive_normalization(alpha: f64, l: u32) -> f64 {
    // Base normalization: (2*alpha/pi)^(3/4)
    let base = (2.0 * alpha / PI).powf(0.75);

    match l {
        0 => {
            // s-type: N = (2*alpha/pi)^(3/4)
            // Representative Cartesian: (0,0,0), no angular factor
            base
        }
        1 => {
            // p-type: N = (2*alpha/pi)^(3/4) * (4*alpha)^(1/2)
            // Representative Cartesian: (0,0,1)
            // Double factorials: (-1)!!=1, (-1)!!=1, (1)!!=1, product=1
            base * (4.0 * alpha).sqrt()
        }
        2 => {
            // d-type: N = (2*alpha/pi)^(3/4) * (4*alpha)
            // Representative Cartesian: (1,1,0)
            // Powers (1,1,0): (2*1-1)!!=1, (2*1-1)!!=1, (-1)!!=1
            // Angular factor: (4*alpha)^(2/2) / sqrt(1*1*1) = 4*alpha
            base * 4.0 * alpha
        }
        _ => {
            // Unsupported angular momentum; fall back to base
            base
        }
    }
}

/// Compute the auto-range r_max for a set of exponents.
///
/// Returns the r value where the most diffuse Gaussian (smallest exponent)
/// falls below a visual threshold, ensuring the full tail is visible.
///
/// Formula: `r_max = min(5.0 / sqrt(alpha_min), 15.0)`
///
/// This places the cutoff where `exp(-alpha_min * r_max^2) = exp(-25) ~ 1.4e-11`,
/// well below visual resolution.
///
/// # Arguments
///
/// * `exponents` - Slice of primitive Gaussian exponents
///
/// # Returns
///
/// The computed r_max in Bohr, clamped to [0.5, 15.0].
///
/// # Panics
///
/// Panics if `exponents` is empty.
pub fn auto_r_max(exponents: &[f64]) -> f64 {
    assert!(!exponents.is_empty(), "exponents must not be empty");

    let alpha_min = exponents.iter().copied().fold(f64::INFINITY, f64::min);

    // r_max = 5.0 / sqrt(alpha_min), clamped to [0.5, 15.0]
    let r = 5.0 / alpha_min.sqrt();
    r.clamp(0.5, 15.0)
}

/// Evaluate the radial profile of a contracted shell.
///
/// Computes R(r) = r^l * sum_k { d_k * exp(-alpha_k * r^2) }
/// where d_k = c_k * N_k (contraction coefficient * primitive normalization).
///
/// The evaluation produces both the contracted profile and individual
/// primitive profiles for overlay plotting.
///
/// # Arguments
///
/// * `angular_momentum` - Shell angular momentum (0=s, 1=p, 2=d)
/// * `exponents` - Primitive Gaussian exponents alpha_k
/// * `coefficients` - Raw contraction coefficients c_k
/// * `n_points` - Number of r values to evaluate (default: 200)
/// * `r_max_opt` - Optional maximum r value in Bohr (auto-determined if None)
///
/// # Returns
///
/// `RadialProfileResult` containing r values, contracted profile,
/// individual primitive profiles, and metadata for plotting.
///
/// # Panics
///
/// Panics if `exponents` and `coefficients` have different lengths,
/// or if either is empty.
pub fn evaluate_radial_profile(
    angular_momentum: u32,
    exponents: &[f64],
    coefficients: &[f64],
    n_points: Option<usize>,
    r_max_opt: Option<f64>,
) -> RadialProfileResult {
    assert_eq!(
        exponents.len(),
        coefficients.len(),
        "exponents and coefficients must have the same length"
    );
    assert!(!exponents.is_empty(), "must have at least one primitive");

    let n_prims = exponents.len();
    let n_pts = n_points.unwrap_or(200);
    let r_max = r_max_opt.unwrap_or_else(|| auto_r_max(exponents));

    // Compute effective coefficients: d_k = c_k * N_k
    let effective_coeffs: Vec<f64> = exponents
        .iter()
        .zip(coefficients.iter())
        .map(|(&alpha, &c)| c * primitive_normalization(alpha, angular_momentum))
        .collect();

    // Generate linearly-spaced r values from 0 to r_max
    let r_values: Vec<f64> = if n_pts <= 1 {
        vec![0.0]
    } else {
        (0..n_pts)
            .map(|i| r_max * (i as f64) / ((n_pts - 1) as f64))
            .collect()
    };

    // Evaluate per-primitive and contracted profiles
    let mut primitive_values: Vec<Vec<f64>> = vec![vec![0.0; n_pts]; n_prims];
    let mut contracted_values = vec![0.0; n_pts];

    for (k, (&alpha, &d_k)) in exponents.iter().zip(effective_coeffs.iter()).enumerate() {
        for (i, &r) in r_values.iter().enumerate() {
            // R_k(r) = r^l * d_k * exp(-alpha_k * r^2)
            let r_l = r.powi(angular_momentum as i32);
            let gauss = (-alpha * r * r).exp();
            let val = r_l * d_k * gauss;

            primitive_values[k][i] = val;
            contracted_values[i] += val;
        }
    }

    let angular_momentum_label = match angular_momentum {
        0 => "s".to_string(),
        1 => "p".to_string(),
        2 => "d".to_string(),
        _ => format!("l={}", angular_momentum),
    };

    RadialProfileResult {
        r_values,
        contracted_values,
        primitive_values,
        exponents: exponents.to_vec(),
        effective_coefficients: effective_coeffs,
        raw_coefficients: coefficients.to_vec(),
        angular_momentum,
        angular_momentum_label,
        n_primitives: n_prims,
        r_max,
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    // =========================================================================
    // H 1s STO-3G reference data (from PySCF 2.11.0)
    // =========================================================================

    /// Exponents for hydrogen 1s STO-3G (PySCF 2.11.0)
    const H_1S_EXPONENTS: [f64; 3] = [3.4252509100, 0.6239137300, 0.1688554000];
    /// Raw contraction coefficients for hydrogen 1s STO-3G (PySCF 2.11.0)
    const H_1S_COEFFICIENTS: [f64; 3] = [0.1543289707, 0.5353281424, 0.4446345420];

    // O 2p STO-3G reference data (from PySCF 2.11.0)
    const O_2P_EXPONENTS: [f64; 3] = [5.0331513000, 1.1695961000, 0.3803890000];
    const O_2P_COEFFICIENTS: [f64; 3] = [0.1559162685, 0.6076837141, 0.3919573862];

    // =========================================================================
    // PySCF cross-validation (golden values from eval_gto)
    // =========================================================================

    #[test]
    fn test_h_1s_sto3g_at_origin() {
        // PySCF: r=0.0  chi(r)=6.282468823221000e-01
        // For s-type at r=0: R(0) = sum_k d_k * exp(0) = sum_k d_k
        let result =
            evaluate_radial_profile(0, &H_1S_EXPONENTS, &H_1S_COEFFICIENTS, Some(10), None);
        let r0_value = result.contracted_values[0];
        let expected_sum = result.effective_coefficients.iter().sum::<f64>();

        // R(0) must equal sum of effective coefficients
        assert_abs_diff_eq!(r0_value, expected_sum, epsilon = 1e-14);

        // Cross-validate against PySCF eval_gto value
        assert_abs_diff_eq!(r0_value, 6.282468823221000e-01, epsilon = 1e-10);
    }

    #[test]
    fn test_h_1s_sto3g_at_r1() {
        // PySCF: r=1.0  chi(r)=2.230358271219306e-01
        let result = evaluate_radial_profile(
            0,
            &H_1S_EXPONENTS,
            &H_1S_COEFFICIENTS,
            Some(201),
            Some(10.0),
        );

        // Find the value closest to r=1.0
        // With n=201, spacing=10/200=0.05, index for r=1.0 is 20
        let idx = 20;
        assert_abs_diff_eq!(result.r_values[idx], 1.0, epsilon = 1e-14);
        assert_abs_diff_eq!(
            result.contracted_values[idx],
            2.230358271219306e-01,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_h_1s_sto3g_at_multiple_r() {
        // PySCF golden values for H 1s STO-3G along z-axis
        let reference: &[(f64, f64)] = &[
            (0.0, 6.282468823221000e-01),
            (0.5, 4.268012209132965e-01),
            (1.0, 2.230358271219306e-01),
            (1.5, 1.230106274332335e-01),
            (2.0, 6.456483925938930e-02),
            (3.0, 1.923765504354119e-02),
            (4.0, 5.612811274401157e-03),
            (5.0, 1.225296708600061e-03),
        ];

        // Also generate a full profile to verify the first point
        let result = evaluate_radial_profile(0, &H_1S_EXPONENTS, &H_1S_COEFFICIENTS, None, None);

        // For each reference point, find the closest r in the grid or evaluate directly
        for &(r_ref, expected) in reference {
            // Compute value directly for exact r
            let mut val = 0.0;
            for (&alpha, &coef) in H_1S_EXPONENTS.iter().zip(H_1S_COEFFICIENTS.iter()) {
                let norm = primitive_normalization(alpha, 0);
                val += coef * norm * (-alpha * r_ref * r_ref).exp();
            }
            assert_abs_diff_eq!(val, expected, epsilon = 1e-10,);
        }

        // Also verify first point of the generated result
        assert_abs_diff_eq!(result.contracted_values[0], reference[0].1, epsilon = 1e-10);
    }

    // =========================================================================
    // s-type profile shape tests
    // =========================================================================

    #[test]
    fn test_h_1s_sto3g_monotonic_decay() {
        // s-type with all positive coefficients should decay monotonically from r=0
        let result =
            evaluate_radial_profile(0, &H_1S_EXPONENTS, &H_1S_COEFFICIENTS, Some(200), None);
        for i in 1..result.contracted_values.len() {
            assert!(
                result.contracted_values[i] <= result.contracted_values[i - 1] + 1e-15,
                "s-type H 1s profile must be monotonically non-increasing, violated at i={}: {} > {}",
                i,
                result.contracted_values[i],
                result.contracted_values[i - 1],
            );
        }
    }

    #[test]
    fn test_h_1s_sto3g_positive() {
        // s-type with all positive coefficients should be non-negative everywhere
        let result =
            evaluate_radial_profile(0, &H_1S_EXPONENTS, &H_1S_COEFFICIENTS, Some(200), None);
        for (i, &v) in result.contracted_values.iter().enumerate() {
            assert!(
                v >= -1e-15,
                "s-type profile should be non-negative at index {}, got {}",
                i,
                v
            );
        }
    }

    // =========================================================================
    // p-type profile shape tests
    // =========================================================================

    #[test]
    fn test_p_type_zero_at_origin() {
        // p-type: R(0) = 0^1 * sum_k d_k * exp(0) = 0
        let result =
            evaluate_radial_profile(1, &O_2P_EXPONENTS, &O_2P_COEFFICIENTS, Some(200), None);
        assert_abs_diff_eq!(result.contracted_values[0], 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_p_type_has_maximum() {
        // p-type should have a maximum at r > 0, then decay
        let result =
            evaluate_radial_profile(1, &O_2P_EXPONENTS, &O_2P_COEFFICIENTS, Some(200), None);
        let max_val = result
            .contracted_values
            .iter()
            .copied()
            .fold(0.0_f64, f64::max);
        assert!(max_val > 0.0, "p-type profile must have a positive maximum");
        let max_idx = result
            .contracted_values
            .iter()
            .position(|&v| v == max_val)
            .unwrap();
        assert!(
            max_idx > 0,
            "p-type maximum must not be at r=0 (found at index {})",
            max_idx
        );
    }

    #[test]
    fn test_p_type_pyscf_crossvalidation() {
        // PySCF: O 2pz at r=0.5 -> 7.071411779390142e-01
        // PySCF: O 2pz at r=1.0 -> 4.521398203032004e-01
        // The pz AO along z-axis: chi_pz(0,0,z) = z * R_radial(z) = R(z) for our convention
        // since R(r) = r^l * sum_k d_k exp(-alpha_k r^2) and for pz along z, z=r
        let r_vals = [0.5, 1.0, 1.5, 2.0];
        let expected = [
            7.071411779390142e-01,
            4.521398203032004e-01,
            2.201369187979187e-01,
            9.247754891306965e-02,
        ];

        for (&r, &exp_val) in r_vals.iter().zip(expected.iter()) {
            let mut val = 0.0;
            for (&alpha, &coef) in O_2P_EXPONENTS.iter().zip(O_2P_COEFFICIENTS.iter()) {
                let norm = primitive_normalization(alpha, 1);
                val += r * coef * norm * (-alpha * r * r).exp();
            }
            assert_abs_diff_eq!(val, exp_val, epsilon = 1e-10,);
        }
    }

    // =========================================================================
    // d-type profile shape tests
    // =========================================================================

    #[test]
    fn test_d_type_zero_at_origin() {
        // d-type: R(0) = 0^2 * ... = 0
        // Use some arbitrary d-type exponents and coefficients
        let d_exponents = [5.0, 1.0, 0.3];
        let d_coefficients = [0.2, 0.5, 0.3];
        let result = evaluate_radial_profile(2, &d_exponents, &d_coefficients, Some(200), None);
        assert_abs_diff_eq!(result.contracted_values[0], 0.0, epsilon = 1e-15);
    }

    // =========================================================================
    // Contracted equals sum of primitives
    // =========================================================================

    #[test]
    fn test_contracted_equals_sum_of_primitives() {
        let result =
            evaluate_radial_profile(0, &H_1S_EXPONENTS, &H_1S_COEFFICIENTS, Some(50), None);
        for (r_idx, &contracted) in result.contracted_values.iter().enumerate() {
            let sum: f64 = result.primitive_values.iter().map(|pv| pv[r_idx]).sum();
            assert_abs_diff_eq!(contracted, sum, epsilon = 1e-14,);
        }
    }

    #[test]
    fn test_contracted_equals_sum_p_type() {
        let result =
            evaluate_radial_profile(1, &O_2P_EXPONENTS, &O_2P_COEFFICIENTS, Some(50), None);
        for (r_idx, &contracted) in result.contracted_values.iter().enumerate() {
            let sum: f64 = result.primitive_values.iter().map(|pv| pv[r_idx]).sum();
            assert_abs_diff_eq!(contracted, sum, epsilon = 1e-14,);
        }
    }

    // =========================================================================
    // Auto-range tests
    // =========================================================================

    #[test]
    fn test_auto_range_h_1s_sto3g() {
        // H 1s STO-3G most diffuse exponent: 0.168855
        // r_max = min(5.0 / sqrt(0.168855), 15.0) = min(12.17, 15.0) = 12.17
        let r_max = auto_r_max(&H_1S_EXPONENTS);
        assert!(
            r_max > 10.0 && r_max < 15.0,
            "r_max = {} should be ~12.2 for H STO-3G",
            r_max
        );
        // More precise check
        let expected = 5.0 / H_1S_EXPONENTS[2].sqrt();
        assert_abs_diff_eq!(r_max, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_auto_range_tight_exponent() {
        // Very tight Gaussian: alpha = 100.0
        // r_max = min(5.0 / sqrt(100.0), 15.0) = min(0.5, 15.0) = 0.5
        let r_max = auto_r_max(&[100.0]);
        assert_abs_diff_eq!(r_max, 0.5, epsilon = 0.01);
    }

    #[test]
    fn test_auto_range_very_diffuse() {
        // Very diffuse Gaussian: alpha = 0.01
        // r_max = min(5.0 / sqrt(0.01), 15.0) = min(50.0, 15.0) = 15.0
        let r_max = auto_r_max(&[0.01]);
        assert_abs_diff_eq!(r_max, 15.0, epsilon = 1e-10);
    }

    #[test]
    fn test_auto_range_uses_most_diffuse() {
        // With multiple exponents, should use the smallest (most diffuse)
        let r_max_all = auto_r_max(&[100.0, 10.0, 0.5]);
        let r_max_diffuse = auto_r_max(&[0.5]);
        assert_abs_diff_eq!(r_max_all, r_max_diffuse, epsilon = 1e-14);
    }

    // =========================================================================
    // Normalization tests
    // =========================================================================

    #[test]
    fn test_primitive_normalization_s_type() {
        // N(alpha, l=0) = (2*alpha/pi)^(3/4)
        let alpha = 1.0;
        let n = primitive_normalization(alpha, 0);
        let expected = (2.0 * alpha / PI).powf(0.75);
        assert_abs_diff_eq!(n, expected, epsilon = 1e-14);
    }

    #[test]
    fn test_primitive_normalization_p_type() {
        // N(alpha, l=1) = (2*alpha/pi)^(3/4) * (4*alpha)^(1/2)
        let alpha = 1.0;
        let n = primitive_normalization(alpha, 1);
        let expected = (2.0 * alpha / PI).powf(0.75) * (4.0 * alpha).sqrt();
        assert_abs_diff_eq!(n, expected, epsilon = 1e-14);
    }

    #[test]
    fn test_primitive_normalization_d_type() {
        // N(alpha, l=2) = (2*alpha/pi)^(3/4) * 4*alpha
        // For d_{xy} (powers 1,1,0): no double-factorial correction
        let alpha = 1.0;
        let n = primitive_normalization(alpha, 2);
        let expected = (2.0 * alpha / PI).powf(0.75) * 4.0 * alpha;
        assert_abs_diff_eq!(n, expected, epsilon = 1e-14);
    }

    #[test]
    fn test_normalization_matches_integral_engine() {
        // Verify that primitive_normalization for each l
        // matches cartesian_gaussian_normalization with representative powers
        use crate::integrals::overlap::cartesian_gaussian_normalization;
        use crate::integrals::CartesianPower;

        let alpha = 2.5;

        // s-type: powers (0,0,0)
        let n_radial = primitive_normalization(alpha, 0);
        let n_integral = cartesian_gaussian_normalization(alpha, &CartesianPower::new(0, 0, 0));
        assert_abs_diff_eq!(n_radial, n_integral, epsilon = 1e-14);

        // p-type: powers (0,0,1)
        let n_radial = primitive_normalization(alpha, 1);
        let n_integral = cartesian_gaussian_normalization(alpha, &CartesianPower::new(0, 0, 1));
        assert_abs_diff_eq!(n_radial, n_integral, epsilon = 1e-14);

        // d-type: powers (1,1,0)
        let n_radial = primitive_normalization(alpha, 2);
        let n_integral = cartesian_gaussian_normalization(alpha, &CartesianPower::new(1, 1, 0));
        assert_abs_diff_eq!(n_radial, n_integral, epsilon = 1e-14);
    }

    #[test]
    fn test_normalization_various_alphas() {
        // Test normalization is positive for multiple alpha values
        for &alpha in &[0.1, 0.5, 1.0, 5.0, 50.0, 130.0] {
            let n_s = primitive_normalization(alpha, 0);
            assert!(
                n_s > 0.0,
                "s-type norm must be positive for alpha={}",
                alpha
            );

            let n_p = primitive_normalization(alpha, 1);
            assert!(
                n_p > 0.0,
                "p-type norm must be positive for alpha={}",
                alpha
            );

            let n_d = primitive_normalization(alpha, 2);
            assert!(
                n_d > 0.0,
                "d-type norm must be positive for alpha={}",
                alpha
            );

            // Verify ratios match analytical formulas:
            // N_p / N_s = (4*alpha)^(1/2) = 2*sqrt(alpha)
            let ratio_ps = n_p / n_s;
            assert_abs_diff_eq!(ratio_ps, 2.0 * alpha.sqrt(), epsilon = 1e-12);

            // N_d / N_s = 4*alpha
            let ratio_ds = n_d / n_s;
            assert_abs_diff_eq!(ratio_ds, 4.0 * alpha, epsilon = 1e-12);
        }
    }

    // =========================================================================
    // Output structure tests
    // =========================================================================

    #[test]
    fn test_output_dimensions() {
        let n_points = 100;
        let result =
            evaluate_radial_profile(0, &H_1S_EXPONENTS, &H_1S_COEFFICIENTS, Some(n_points), None);
        assert_eq!(result.r_values.len(), n_points);
        assert_eq!(result.contracted_values.len(), n_points);
        assert_eq!(result.primitive_values.len(), 3); // 3 primitives in STO-3G
        for pv in &result.primitive_values {
            assert_eq!(pv.len(), n_points);
        }
        assert_eq!(result.exponents.len(), 3);
        assert_eq!(result.effective_coefficients.len(), 3);
        assert_eq!(result.raw_coefficients.len(), 3);
        assert_eq!(result.angular_momentum, 0);
        assert_eq!(result.angular_momentum_label, "s");
        assert_eq!(result.n_primitives, 3);
    }

    #[test]
    fn test_r_values_start_at_zero() {
        let result =
            evaluate_radial_profile(0, &H_1S_EXPONENTS, &H_1S_COEFFICIENTS, Some(50), None);
        assert_abs_diff_eq!(result.r_values[0], 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_r_values_end_at_r_max() {
        let r_max = 5.0;
        let result = evaluate_radial_profile(
            0,
            &H_1S_EXPONENTS,
            &H_1S_COEFFICIENTS,
            Some(50),
            Some(r_max),
        );
        assert_abs_diff_eq!(*result.r_values.last().unwrap(), r_max, epsilon = 1e-10);
        assert_abs_diff_eq!(result.r_max, r_max, epsilon = 1e-10);
    }

    #[test]
    fn test_r_values_uniformly_spaced() {
        let n_points = 101;
        let r_max = 10.0;
        let result = evaluate_radial_profile(
            0,
            &H_1S_EXPONENTS,
            &H_1S_COEFFICIENTS,
            Some(n_points),
            Some(r_max),
        );
        let spacing = r_max / (n_points as f64 - 1.0);
        for i in 1..result.r_values.len() {
            let delta = result.r_values[i] - result.r_values[i - 1];
            assert_abs_diff_eq!(delta, spacing, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_metadata_p_type() {
        let result =
            evaluate_radial_profile(1, &O_2P_EXPONENTS, &O_2P_COEFFICIENTS, Some(50), None);
        assert_eq!(result.angular_momentum, 1);
        assert_eq!(result.angular_momentum_label, "p");
        assert_eq!(result.n_primitives, 3);
    }

    #[test]
    fn test_metadata_d_type() {
        let result = evaluate_radial_profile(2, &[5.0, 1.0], &[0.5, 0.5], Some(50), None);
        assert_eq!(result.angular_momentum, 2);
        assert_eq!(result.angular_momentum_label, "d");
        assert_eq!(result.n_primitives, 2);
    }

    #[test]
    fn test_raw_coefficients_preserved() {
        let result =
            evaluate_radial_profile(0, &H_1S_EXPONENTS, &H_1S_COEFFICIENTS, Some(10), None);
        for (i, &c) in result.raw_coefficients.iter().enumerate() {
            assert_abs_diff_eq!(c, H_1S_COEFFICIENTS[i], epsilon = 1e-15);
        }
    }

    #[test]
    fn test_exponents_preserved() {
        let result =
            evaluate_radial_profile(0, &H_1S_EXPONENTS, &H_1S_COEFFICIENTS, Some(10), None);
        for (i, &e) in result.exponents.iter().enumerate() {
            assert_abs_diff_eq!(e, H_1S_EXPONENTS[i], epsilon = 1e-15);
        }
    }

    // =========================================================================
    // Serde round-trip
    // =========================================================================

    #[test]
    fn test_serde_roundtrip() {
        let result =
            evaluate_radial_profile(0, &H_1S_EXPONENTS, &H_1S_COEFFICIENTS, Some(20), None);
        let json = serde_json::to_string(&result).unwrap();
        let deserialized: RadialProfileResult = serde_json::from_str(&json).unwrap();

        assert_eq!(result.r_values.len(), deserialized.r_values.len());
        assert_eq!(result.angular_momentum, deserialized.angular_momentum);
        assert_eq!(result.n_primitives, deserialized.n_primitives);
        for (a, b) in result
            .contracted_values
            .iter()
            .zip(deserialized.contracted_values.iter())
        {
            assert_abs_diff_eq!(a, b, epsilon = 1e-15);
        }
    }
}
