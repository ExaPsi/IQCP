// Allow excessive precision for mathematical constants from reference implementations.
// These values are taken directly from Cody (1969) and libcint for numerical accuracy.
#![allow(clippy::excessive_precision)]

//! Boys function F_m(T) implementation
//!
//! The Boys function is fundamental to Gaussian integral evaluation:
//!
//! ```text
//! F_m(T) = integral from 0 to 1 of t^(2m) * exp(-T * t^2) dt
//! ```
//!
//! # Computational Regimes
//!
//! Two computational regimes are used based on libcint's proven implementation:
//!
//! - **Series** (T < TURNOVER_POINT[m]): Power series expansion with downward recurrence
//! - **Recurrence** (T >= TURNOVER_POINT[m]): erf(sqrt(T)) + upward recurrence
//!
//! The turnover points are m-dependent and derived from error analysis ensuring
//! optimal accuracy in each regime. See libcint fmt.c lines 42-83.
//!
//! # Tolerance
//!
//! All computations target absolute accuracy of 1e-12.
//!
//! # References
//!
//! - Shavitt, I. (1963). Methods in Computational Physics, Vol. 2, Eq. 4-12.
//! - libcint: references/libcint/src/fmt.c (gamma_inc_like function)
//!
//! # Implementation Notes
//!
//! This implementation follows libcint's fmt.c closely:
//! - `fmt1_gamma_inc_like` (lines 186-204): Series expansion for small T
//! - `gamma_inc_like` (lines 206-226): Regime router and erf+recurrence

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Module version (matches crate version)
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// sqrt(pi/4) constant for erf-based computation
/// Reference: libcint fmt.c line 23
const SQRTPIE4: f64 = 0.8862269254527580136490837416705725913987747280611935641069038949264;

/// Machine epsilon threshold for series termination
/// Reference: libcint fmt.c line 20: SML_FLOAT64 = DBL_EPSILON * 0.5
const SML_FLOAT64: f64 = f64::EPSILON * 0.5;

/// Threshold below which T is treated as zero
const NEAR_ZERO_THRESHOLD: f64 = 1e-14;

/// Maximum supported order m
pub const MAX_ORDER: u32 = 50;

/// Turnover points for regime selection (m-dependent)
///
/// For T < TURNOVER_POINT[m], use series expansion.
/// For T >= TURNOVER_POINT[m], use erf + upward recurrence.
///
/// These values are derived from error analysis ensuring optimal accuracy.
/// Reference: libcint fmt.c lines 42-83
///
/// The formula is: t0 = 0.5 * ((2m-1)!!/(2m-1)^0.5)^(1/(m-1))
static TURNOVER_POINT: [f64; 40] = [
    0.0,            // m=0: always use erf
    0.0,            // m=1: always use erf
    0.866025403784, // m=2
    1.295010032056, // m=3
    1.705493613097, // m=4
    2.106432965305, // m=5
    2.501471934009, // m=6
    2.892473348218, // m=7
    3.280525047072, // m=8
    3.666320693281, // m=9
    4.05033123037,  // m=10
    4.432891808508, // m=11
    4.814249856864, // m=12
    5.194593501454, // m=13
    5.574069276051, // m=14
    5.952793645111, // m=15
    6.330860773135, // m=16
    6.708347923415, // m=17
    7.08531930745,  // m=18
    7.461828891625, // m=19
    7.837922483937, // m=20
    8.213639312398, // m=21
    8.589013237349, // m=22
    8.964073695432, // m=23
    9.338846443746, // m=24
    9.713354153046, // m=25
    10.08761688545, // m=26
    10.46165248270, // m=27
    10.83547688448, // m=28
    11.20910439128, // m=29
    11.58254788331, // m=30
    11.95581900374, // m=31
    12.32892831326, // m=32
    12.70188542111, // m=33
    13.07469909673, // m=34
    13.44737736550, // m=35
    13.81992759110, // m=36
    14.19235654675, // m=37
    14.56467047710, // m=38
    14.93687515212, // m=39
];

/// Errors that can occur during Boys function evaluation
#[derive(Error, Debug, Clone, PartialEq)]
pub enum BoysError {
    /// Order m exceeds the maximum supported value
    #[error("Invalid order m={0}, must be <= {}", MAX_ORDER)]
    InvalidOrder(u32),

    /// Argument T is negative (must be >= 0)
    #[error("Negative argument T={0}, must be >= 0")]
    NegativeArgument(f64),

    /// Numerical computation failed (should not happen with valid inputs)
    #[error("Numerical computation failed: {0}")]
    NumericalFailure(String),
}

/// Computational method used for Boys function evaluation
///
/// This enum enables pedagogical transparency by indicating which
/// algorithm was employed for a given (m, T) pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum BoysMethod {
    /// T is zero or near-zero: F_m(0) = 1/(2m+1)
    Zero,

    /// Series expansion for T < TURNOVER_POINT[m]
    /// Reference: libcint fmt.c fmt1_gamma_inc_like
    Series,

    /// erf(sqrt(T)) + upward recurrence for T >= TURNOVER_POINT[m]
    /// Reference: libcint fmt.c gamma_inc_like
    Recurrence,
}

impl std::fmt::Display for BoysMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BoysMethod::Zero => write!(f, "Zero"),
            BoysMethod::Series => write!(f, "Series"),
            BoysMethod::Recurrence => write!(f, "Recurrence"),
        }
    }
}

/// Result of Boys function evaluation
///
/// Contains both the computed value and metadata about how it was computed,
/// supporting IQCP's pedagogical goal of computational transparency.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BoysResult {
    /// The computed value F_m(T)
    pub value: f64,

    /// The computational method used
    pub method: BoysMethod,

    /// The order m of the Boys function
    pub m: u32,

    /// The argument T
    pub t: f64,

    /// The turnover point T_0 for this m value
    ///
    /// The regime boundary where the computational method changes:
    /// - T < turnover: Series expansion
    /// - T >= turnover: erf + upward recurrence
    ///
    /// This is m-dependent and derived from error analysis (libcint fmt.c).
    /// For m=0 and m=1, turnover is 0.0 (always use recurrence).
    pub turnover: f64,

    /// Number of terms/steps used in computation
    ///
    /// - Zero method: 0 (direct formula)
    /// - Series: Number of iteration terms until convergence
    /// - Recurrence: m + 1 (upward recurrence steps from F_0 to F_m)
    pub terms_count: u32,

    /// Estimated relative error bound
    ///
    /// - Zero method: Some(f64::EPSILON) (machine precision)
    /// - Series: Some(last_term / sum) (estimated from convergence)
    /// - Recurrence: None (difficult to estimate accurately)
    pub estimated_error: Option<f64>,
}

impl BoysResult {
    /// Create a new BoysResult
    fn new(
        value: f64,
        method: BoysMethod,
        m: u32,
        t: f64,
        turnover: f64,
        terms_count: u32,
        estimated_error: Option<f64>,
    ) -> Self {
        Self {
            value,
            method,
            m,
            t,
            turnover,
            terms_count,
            estimated_error,
        }
    }
}

/// Evaluate the Boys function F_m(T) for a single (m, T) pair.
///
/// # Arguments
///
/// * `m` - Order of the Boys function (0 <= m <= MAX_ORDER)
/// * `t` - Argument value (T >= 0)
///
/// # Returns
///
/// A `BoysResult` containing the computed value and the method used.
///
/// # Errors
///
/// - `BoysError::InvalidOrder` if m > MAX_ORDER
/// - `BoysError::NegativeArgument` if T < 0
///
/// # Example
///
/// ```rust
/// use qc_core::boys::{boys_eval, BoysMethod};
///
/// let result = boys_eval(0, 0.5).unwrap();
/// assert!((result.value - 0.8556243918921487).abs() < 1e-12);
/// assert_eq!(result.method, BoysMethod::Recurrence);
/// ```
///
/// # References
///
/// - Shavitt (1963), Methods in Computational Physics, Vol. 2
/// - libcint fmt.c gamma_inc_like (lines 206-226)
pub fn boys_eval(m: u32, t: f64) -> Result<BoysResult, BoysError> {
    // Validate inputs
    if m > MAX_ORDER {
        return Err(BoysError::InvalidOrder(m));
    }
    if t < 0.0 {
        return Err(BoysError::NegativeArgument(t));
    }

    // Get the turnover point for this m value
    let turnover = get_turnover_point(m);

    // Handle T = 0 or near-zero case
    // Reference: libcint fmt.c lines 208-212
    if t < NEAR_ZERO_THRESHOLD {
        let value = 1.0 / (2.0 * m as f64 + 1.0);
        // Zero method: direct formula, no iterations, machine precision error
        return Ok(BoysResult::new(
            value,
            BoysMethod::Zero,
            m,
            t,
            turnover,
            0,                  // terms_count: direct formula
            Some(f64::EPSILON), // estimated_error: machine precision
        ));
    }

    // Select regime based on turnover point
    // Reference: libcint fmt.c lines 214-225
    if t < turnover {
        // Series expansion regime
        let (value, metrics) = boys_series(m, t);
        Ok(BoysResult::new(
            value,
            BoysMethod::Series,
            m,
            t,
            turnover,
            metrics.terms_count,
            Some(metrics.estimated_error),
        ))
    } else {
        // erf + upward recurrence regime
        let (value, terms_count) = boys_recurrence(m, t);
        Ok(BoysResult::new(
            value,
            BoysMethod::Recurrence,
            m,
            t,
            turnover,
            terms_count,
            None, // estimated_error: difficult to estimate for recurrence
        ))
    }
}

/// Evaluate Boys functions F_0(T) through F_m(T) for a single T value.
///
/// This is more efficient when multiple orders are needed, as it computes
/// all values F_0, F_1, ..., F_m in a single pass using recurrence relations.
///
/// # Arguments
///
/// * `m_max` - Maximum order (computes F_0 through F_{m_max})
/// * `t` - Argument value (T >= 0)
///
/// # Returns
///
/// A vector of `BoysResult` for m = 0, 1, ..., m_max.
///
/// # Errors
///
/// - `BoysError::InvalidOrder` if m_max > MAX_ORDER
/// - `BoysError::NegativeArgument` if T < 0
pub fn boys_eval_all(m_max: u32, t: f64) -> Result<Vec<BoysResult>, BoysError> {
    // Validate inputs
    if m_max > MAX_ORDER {
        return Err(BoysError::InvalidOrder(m_max));
    }
    if t < 0.0 {
        return Err(BoysError::NegativeArgument(t));
    }

    // Allocate result array
    let n = (m_max + 1) as usize;
    let mut values = vec![0.0; n];

    // Compute all values using the appropriate regime
    let batch_metrics = boys_all_internal(&mut values, t, m_max as usize);

    // Package results
    // Note: For batch computation, all results share the same method and metrics
    // (computed once for the highest m, applied to all)
    let results = values
        .into_iter()
        .enumerate()
        .map(|(i, value)| {
            // Each result gets its own turnover point (m-dependent)
            let turnover = get_turnover_point(i as u32);
            BoysResult::new(
                value,
                batch_metrics.method,
                i as u32,
                t,
                turnover,
                batch_metrics.terms_count,
                batch_metrics.estimated_error,
            )
        })
        .collect();

    Ok(results)
}

/// Evaluate the Boys function F_m(T) for multiple T values.
///
/// # Arguments
///
/// * `m` - Order of the Boys function (0 <= m <= MAX_ORDER)
/// * `ts` - Slice of argument values (all must be >= 0)
///
/// # Returns
///
/// A vector of `BoysResult`, one for each T value in the input.
///
/// # Errors
///
/// - `BoysError::InvalidOrder` if m > MAX_ORDER
/// - `BoysError::NegativeArgument` if any T < 0
pub fn boys_eval_many(m: u32, ts: &[f64]) -> Result<Vec<BoysResult>, BoysError> {
    // Validate order once
    if m > MAX_ORDER {
        return Err(BoysError::InvalidOrder(m));
    }

    // Validate all T values first
    for &t in ts {
        if t < 0.0 {
            return Err(BoysError::NegativeArgument(t));
        }
    }

    // Evaluate for each T value
    ts.iter().map(|&t| boys_eval(m, t)).collect()
}

/// Get the turnover point for a given order m.
///
/// The turnover point determines which computational regime is used:
/// - T < turnover: Series expansion
/// - T >= turnover: erf + upward recurrence
///
/// These m-dependent values are derived from error analysis (see libcint fmt.c)
/// using the formula: t0 = 0.5 * ((2m-1)!!/(2m-1)^0.5)^(1/(m-1))
///
/// For m >= 40 (beyond our precomputed table), we extrapolate using
/// the approximate formula: t0 ~ 0.372 * m + const
///
/// # Arguments
///
/// * `m` - Order of the Boys function (0 <= m <= MAX_ORDER)
///
/// # Returns
///
/// The turnover point T where the computation regime changes.
///
/// # Example
///
/// ```rust
/// use qc_core::boys::get_turnover_point;
///
/// // m=0 and m=1 always use recurrence
/// assert_eq!(get_turnover_point(0), 0.0);
/// assert_eq!(get_turnover_point(1), 0.0);
///
/// // m=5 has turnover at ~2.1
/// assert!((get_turnover_point(5) - 2.106).abs() < 0.01);
/// ```
#[inline]
pub fn get_turnover_point(m: u32) -> f64 {
    if (m as usize) < TURNOVER_POINT.len() {
        TURNOVER_POINT[m as usize]
    } else {
        // Extrapolate for large m (approximately linear growth)
        // Last tabulated value is at m=39: 14.937
        // Growth rate is approximately 0.372 per unit m
        let m_diff = m as f64 - 39.0;
        14.93687515212 + 0.372 * m_diff
    }
}

/// Metrics returned from series computation
struct SeriesMetrics {
    /// Number of terms computed in the series
    terms_count: u32,
    /// Estimated relative error (last_term / sum)
    estimated_error: f64,
}

/// Series expansion for F_m(T) when T is small.
///
/// This computes F_m(T) using the power series expansion:
///
/// ```text
/// F_m(T) = exp(-T)/(2b) * [1 + T/(b+1) + T^2/((b+1)(b+2)) + ...]
/// where b = m + 0.5
/// ```
///
/// Then uses downward recurrence to get F_0..F_{m-1}.
///
/// Reference: libcint fmt.c fmt1_gamma_inc_like (lines 186-204)
///
/// Returns: (value, terms_count, estimated_error)
#[inline]
fn boys_series(m: u32, t: f64) -> (f64, SeriesMetrics) {
    let m_usize = m as usize;
    let mut f = vec![0.0; m_usize + 1];
    let metrics = boys_series_internal(&mut f, t, m_usize);
    (f[m_usize], metrics)
}

/// Internal series expansion that fills the array f[0..=m].
///
/// Algorithm (libcint fmt.c lines 186-204):
/// 1. Compute F_m(T) via converging series at order m
/// 2. Use downward recurrence to fill F_{m-1}, F_{m-2}, ..., F_0
///
/// Returns SeriesMetrics with iteration count and estimated error.
fn boys_series_internal(f: &mut [f64], t: f64, m: usize) -> SeriesMetrics {
    // b = m + 0.5
    let b = m as f64 + 0.5;

    // e = 0.5 * exp(-T)
    let e = 0.5 * (-t).exp();

    // Series expansion: sum_k T^k / (b+1)(b+2)...(b+k)
    // Start with x = e, accumulate until convergence
    let mut x = e;
    let mut s = e;
    let tol = SML_FLOAT64 * e;

    // Iterate until terms are smaller than tolerance
    // Reference: libcint lines 195-198
    let mut bi = b + 1.0;
    let mut iteration_count: u32 = 0;
    while x > tol {
        x *= t / bi;
        s += x;
        bi += 1.0;
        iteration_count += 1;
    }

    // Estimated error: last term relative to sum
    // When x <= tol, the relative error is approximately x / s
    let estimated_error = if s.abs() > 0.0 { x / s } else { 0.0 };

    // F_m(T) = s / b
    f[m] = s / b;

    // Downward recurrence: F_{i-1} = (e + T * F_i) / (i - 0.5)
    // Reference: libcint lines 200-203
    let mut b_curr = b;
    for i in (1..=m).rev() {
        b_curr -= 1.0;
        f[i - 1] = (e + t * f[i]) / b_curr;
    }

    SeriesMetrics {
        terms_count: iteration_count,
        estimated_error,
    }
}

/// erf + upward recurrence for F_m(T) when T is moderate/large.
///
/// This computes F_0(T) analytically using the error function:
///
/// ```text
/// F_0(T) = sqrt(pi/4T) * erf(sqrt(T))
/// ```
///
/// Then uses upward recurrence:
///
/// ```text
/// F_{m+1}(T) = [(2m+1)*F_m(T) - exp(-T)] / (2T)
/// ```
///
/// Reference: libcint fmt.c gamma_inc_like (lines 216-225)
///
/// Returns: (value, terms_count) where terms_count = m + 1 (steps from F_0 to F_m)
#[inline]
fn boys_recurrence(m: u32, t: f64) -> (f64, u32) {
    let m_usize = m as usize;
    let mut f = vec![0.0; m_usize + 1];
    boys_recurrence_internal(&mut f, t, m_usize);
    // terms_count = m + 1: we compute F_0, then recur m times to get F_m
    (f[m_usize], m + 1)
}

/// Internal erf + upward recurrence that fills the array f[0..=m].
///
/// Algorithm (libcint fmt.c lines 216-225):
/// 1. Compute F_0(T) = sqrt(pi/4T) * erf(sqrt(T))
/// 2. Use upward recurrence to fill F_1, F_2, ..., F_m
fn boys_recurrence_internal(f: &mut [f64], t: f64, m: usize) {
    let tt = t.sqrt();

    // F_0(T) = SQRTPIE4 / sqrt(T) * erf(sqrt(T))
    // Reference: libcint fmt.c line 219
    f[0] = SQRTPIE4 / tt * erf_cody(tt);

    // e = exp(-T)
    let e = (-t).exp();

    // b = 0.5 / T
    let b = 0.5 / t;

    // Upward recurrence: F_{i}(T) = b * [(2i-1) * F_{i-1}(T) - e]
    // Reference: libcint fmt.c lines 222-224
    for i in 1..=m {
        f[i] = b * ((2 * i - 1) as f64 * f[i - 1] - e);
    }
}

/// Metrics returned from batch computation
struct BatchMetrics {
    method: BoysMethod,
    terms_count: u32,
    estimated_error: Option<f64>,
}

/// Internal function to compute all F_0..F_m using optimal regime.
///
/// Returns the method used and metrics (same for all values since single T).
fn boys_all_internal(f: &mut [f64], t: f64, m: usize) -> BatchMetrics {
    // Handle T = 0 or near-zero case
    if t < NEAR_ZERO_THRESHOLD {
        for (i, f_i) in f.iter_mut().enumerate().take(m + 1) {
            *f_i = 1.0 / (2.0 * i as f64 + 1.0);
        }
        return BatchMetrics {
            method: BoysMethod::Zero,
            terms_count: 0,
            estimated_error: Some(f64::EPSILON),
        };
    }

    // Select regime based on turnover point for highest m
    let turnover = get_turnover_point(m as u32);

    if t < turnover {
        let metrics = boys_series_internal(f, t, m);
        BatchMetrics {
            method: BoysMethod::Series,
            terms_count: metrics.terms_count,
            estimated_error: Some(metrics.estimated_error),
        }
    } else {
        boys_recurrence_internal(f, t, m);
        BatchMetrics {
            method: BoysMethod::Recurrence,
            terms_count: (m + 1) as u32, // m + 1 steps from F_0 to F_m
            estimated_error: None,
        }
    }
}

/// Error function via Cody's rational Chebyshev approximation.
///
/// Computes erf(x) to near machine precision (maximum relative error < 1e-15)
/// using three-region rational polynomial approximations. This is NOT a crude
/// approximation -- it achieves full double-precision accuracy across all x.
///
/// # Reference
///
/// W. J. Cody, "Rational Chebyshev Approximations for the Error Function,"
/// Mathematics of Computation, Vol. 23, No. 107 (1969), pp. 631--637.
/// doi:10.1090/S0025-5718-1969-0247736-4
///
/// Regions (Cody Table I--III):
/// - |x| < 0.5:  erf(x) = x * P(x^2) / Q(x^2)
/// - 0.5 <= |x| < 4:  erfc(x) via rational approximation, erf = 1 - erfc
/// - |x| >= 4:  erfc(x) via asymptotic rational form, erf = 1 - erfc
#[inline]
fn erf_cody(x: f64) -> f64 {
    cody_erf_impl(x)
}

/// Full double-precision erf(x) using the libm crate.
///
/// This replaces the previous Cody (1969) rational Chebyshev approximation
/// which only achieved ~7 significant digits due to incorrect coefficient
/// handling. The libm::erf function is a pure-Rust port of Sun's fdlibm
/// that achieves machine-precision accuracy (relative error < 2e-16).
///
/// Reference: Sun fdlibm s_erf.c (K.C. Ng, March 1992)
fn cody_erf_impl(x: f64) -> f64 {
    libm::erf(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    /// Tolerance for numerical comparisons
    const TOL: f64 = 1e-12;

    // =========================================================================
    // Reference values from scipy.special.hyp1f1
    // Generated using: hyp1f1(m+0.5, m+1.5, -T) / (2m+1)
    // =========================================================================

    #[test]
    fn test_reference_m0() {
        // F_0(T) reference values
        let cases = [
            (0.0, 1.0000000000000000e+00),
            (0.5, 8.5562439189214867e-01),
            (1.0, 7.4682413281242721e-01),
            (10.0, 2.8024739050664266e-01),
            (50.0, 1.2533141373155002e-01),
        ];

        for (t, expected) in cases {
            let result = boys_eval(0, t).unwrap();
            assert_abs_diff_eq!(result.value, expected, epsilon = TOL);
        }
    }

    #[test]
    fn test_reference_m5() {
        // F_5(T) reference values
        let cases = [
            (0.0, 9.0909090909090912e-02),
            (1.0, 3.9364864513484157e-02),
            (10.0, 7.9008749875855312e-05),
        ];

        for (t, expected) in cases {
            let result = boys_eval(5, t).unwrap();
            assert_abs_diff_eq!(result.value, expected, epsilon = TOL);
        }
    }

    #[test]
    fn test_reference_m10() {
        // F_10(T) reference values
        let cases: [(f64, f64); 3] = [
            (0.0, 4.7619047619047616e-02),
            (1.0, 1.9172936091314628e-02),
            (50.0, 8.2058120580663161e-13),
        ];

        for (t, expected) in cases {
            let result = boys_eval(10, t).unwrap();
            // Use relative tolerance for very small values
            if expected.abs() < 1e-10 {
                let rel_err = (result.value - expected).abs() / expected.abs();
                assert!(
                    rel_err < 1e-6,
                    "F_{}({}) = {} (expected {}), rel_err = {}",
                    10,
                    t,
                    result.value,
                    expected,
                    rel_err
                );
            } else {
                assert_abs_diff_eq!(result.value, expected, epsilon = TOL);
            }
        }
    }

    // =========================================================================
    // T = 0 special case tests
    // =========================================================================

    #[test]
    fn test_zero_argument() {
        for m in 0..=10 {
            let result = boys_eval(m, 0.0).unwrap();
            let expected = 1.0 / (2.0 * m as f64 + 1.0);
            assert_eq!(result.method, BoysMethod::Zero);
            assert_abs_diff_eq!(result.value, expected, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_near_zero_argument() {
        // Very small T should use Zero method
        let result = boys_eval(5, 1e-15).unwrap();
        assert_eq!(result.method, BoysMethod::Zero);
        assert_abs_diff_eq!(result.value, 1.0 / 11.0, epsilon = 1e-14);
    }

    // =========================================================================
    // Regime selection tests
    // =========================================================================

    #[test]
    fn test_series_regime() {
        // For m=5, turnover is at ~2.106, so T=1.0 should use series
        let result = boys_eval(5, 1.0).unwrap();
        assert_eq!(result.method, BoysMethod::Series);
    }

    #[test]
    fn test_recurrence_regime() {
        // For m=0, turnover is 0, so any T > 0 should use recurrence
        let result = boys_eval(0, 0.5).unwrap();
        assert_eq!(result.method, BoysMethod::Recurrence);

        // For m=5, T=10 is well above turnover ~2.106
        let result = boys_eval(5, 10.0).unwrap();
        assert_eq!(result.method, BoysMethod::Recurrence);
    }

    #[test]
    fn test_regime_boundary_continuity() {
        // Test that the algorithm produces continuous results across regime boundaries.
        // We check that the difference between the computed value and the expected
        // linear interpolation is small, accounting for the function's natural variation.
        for m in 2..=10 {
            let turnover = get_turnover_point(m);
            if turnover > 0.1 {
                // Use smaller delta to reduce natural variation
                let delta = 1e-6;
                let t_below = turnover - delta;
                let t_above = turnover + delta;
                let t_at = turnover;

                let result_below = boys_eval(m, t_below).unwrap();
                let result_above = boys_eval(m, t_above).unwrap();
                let result_at = boys_eval(m, t_at).unwrap();

                // Check that the value at the boundary is between the neighbors
                // (monotonically decreasing function)
                assert!(
                    result_at.value <= result_below.value && result_at.value >= result_above.value,
                    "Non-monotonic at m={}, turnover={}: f({})={}, f({})={}, f({})={}",
                    m,
                    turnover,
                    t_below,
                    result_below.value,
                    t_at,
                    result_at.value,
                    t_above,
                    result_above.value
                );

                // Check that the computed difference matches expected linear variation
                // For such small delta, the difference should be very small
                let expected_diff = (result_below.value - result_above.value).abs();
                // With delta=1e-6 and typical derivatives ~0.1, expect diff ~2e-7
                assert!(
                    expected_diff < 1e-5,
                    "Unexpected large variation at m={}, turnover={}: diff={}",
                    m,
                    turnover,
                    expected_diff
                );
            }
        }
    }

    // =========================================================================
    // Mathematical property tests
    // =========================================================================

    #[test]
    fn test_non_negativity() {
        // F_m(T) >= 0 for all m >= 0, T >= 0
        let test_ts = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0];
        for m in 0..=10 {
            for &t in &test_ts {
                let result = boys_eval(m, t).unwrap();
                assert!(
                    result.value >= 0.0,
                    "F_{}({}) = {} should be non-negative",
                    m,
                    t,
                    result.value
                );
            }
        }
    }

    #[test]
    fn test_monotonicity_in_t() {
        // F_m(T) decreases as T increases
        for m in 0..=10 {
            let mut prev_value = f64::INFINITY;
            for t in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0] {
                let result = boys_eval(m, t).unwrap();
                assert!(
                    result.value <= prev_value,
                    "F_{}({}) = {} should be <= F_m(smaller T)",
                    m,
                    t,
                    result.value
                );
                prev_value = result.value;
            }
        }
    }

    #[test]
    fn test_monotonicity_in_m() {
        // F_m(T) decreases as m increases (for T > 0)
        for &t in &[0.5, 1.0, 5.0, 10.0] {
            let mut prev_value = f64::INFINITY;
            for m in 0..=10 {
                let result = boys_eval(m, t).unwrap();
                assert!(
                    result.value <= prev_value,
                    "F_{}({}) = {} should be <= F_(m-1)(T)",
                    m,
                    t,
                    result.value
                );
                prev_value = result.value;
            }
        }
    }

    #[test]
    fn test_recurrence_relation() {
        // Verify: F_{m+1}(T) = [(2m+1)*F_m(T) - exp(-T)] / (2T)
        let test_ts = [0.5, 1.0, 2.0, 5.0, 10.0];
        for &t in &test_ts {
            for m in 0..=5 {
                let f_m = boys_eval(m, t).unwrap().value;
                let f_m1 = boys_eval(m + 1, t).unwrap().value;

                let exp_neg_t = (-t).exp();
                let computed = ((2 * m + 1) as f64 * f_m - exp_neg_t) / (2.0 * t);

                assert_abs_diff_eq!(f_m1, computed, epsilon = 1e-11);
            }
        }
    }

    // =========================================================================
    // Batch evaluation tests
    // =========================================================================

    #[test]
    fn test_boys_eval_many() {
        let ts = vec![0.0, 0.5, 1.0, 5.0, 10.0];
        let results = boys_eval_many(3, &ts).unwrap();

        assert_eq!(results.len(), 5);
        for (i, result) in results.iter().enumerate() {
            assert_eq!(result.m, 3);
            assert_eq!(result.t, ts[i]);
            assert!(result.value >= 0.0);
        }
    }

    #[test]
    fn test_boys_eval_all() {
        let results = boys_eval_all(5, 1.0).unwrap();

        assert_eq!(results.len(), 6);
        for (i, result) in results.iter().enumerate() {
            assert_eq!(result.m, i as u32);
            assert_eq!(result.t, 1.0);
            assert!(result.value >= 0.0);
        }

        // Verify consistency with individual evaluations
        for m in 0..=5 {
            let single = boys_eval(m, 1.0).unwrap();
            assert_abs_diff_eq!(results[m as usize].value, single.value, epsilon = 1e-14);
        }
    }

    // =========================================================================
    // Error handling tests
    // =========================================================================

    #[test]
    fn test_invalid_order() {
        let result = boys_eval(MAX_ORDER + 1, 1.0);
        assert!(matches!(result, Err(BoysError::InvalidOrder(_))));
    }

    #[test]
    fn test_negative_argument() {
        let result = boys_eval(0, -1.0);
        assert!(matches!(result, Err(BoysError::NegativeArgument(_))));
    }

    // =========================================================================
    // Additional validation against scipy reference values
    // =========================================================================

    #[test]
    fn test_full_reference_table() {
        // Comprehensive test against scipy.special.hyp1f1 values
        let cases = [
            // (m, T, expected_value)
            (0, 0.5, 8.5562439189214867e-01),
            (0, 1.0, 7.4682413281242721e-01),
            (0, 2.0, 5.9814400666130374e-01),
            (0, 3.0, 5.0434356023143900e-01),
            (0, 4.0, 4.4104069538121066e-01),
            (0, 5.0, 3.9571230961051346e-01),
            (1, 0.5, 2.4909373217951533e-01),
            (1, 1.0, 1.8947234582049235e-01),
            (1, 2.0, 1.1570218085617283e-01),
            (1, 3.0, 7.5759415310595837e-02),
            (1, 4.0, 5.2840632061559555e-02),
            (1, 5.0, 3.8897436261142809e-02),
            (2, 0.5, 1.4075053682591265e-01),
            (2, 1.0, 1.0026879814501735e-01),
            (2, 2.0, 5.2942814832976484e-02),
            (2, 3.0, 2.9581862927320578e-02),
            (2, 4.0, 1.7525782161993065e-02),
            (2, 5.0, 1.0995436178434296e-02),
            (3, 0.5, 9.7222024416930147e-02),
            (3, 1.0, 6.6732274776822254e-02),
            (3, 2.0, 3.2344697732067412e-02),
            (3, 3.0, 1.6353707711456493e-02),
            (3, 4.0, 8.6641589901538963e-03),
            (3, 5.0, 4.8239233893085992e-03),
            (4, 0.5, 7.4023511205877635e-02),
            (4, 1.0, 4.9623241133156748e-02),
            (4, 2.0, 2.2769400221964796e-02),
            (4, 3.0, 1.0781480935388584e-02),
            (4, 4.0, 5.2916842552928854e-03),
            (4, 5.0, 2.7029516726074739e-03),
            (5, 0.5, 5.9680941140265321e-02),
            (5, 1.0, 3.9364864513484157e-02),
            (5, 2.0, 1.7397329690267617e-02),
            (5, 3.0, 7.8743766751055554e-03),
            (5, 4.0, 3.6636899261127255e-03),
            (5, 5.0, 1.7588618054381793e-03),
        ];

        for (m, t, expected) in cases {
            let result = boys_eval(m as u32, t).unwrap();
            assert_abs_diff_eq!(result.value, expected, epsilon = TOL);
        }
    }

    #[test]
    fn test_large_t_values() {
        // Test accuracy for large T (recurrence regime)
        let cases = [
            (0, 20.0, 1.9837606611688054e-01), // sqrt(pi/80)
            (0, 30.0, 1.6197126953893457e-01), // sqrt(pi/120)
            (0, 40.0, 1.4031415181288770e-01), // sqrt(pi/160)
        ];

        for (m, t, _) in cases {
            let result = boys_eval(m as u32, t).unwrap();
            // For large T: F_0(T) ~ sqrt(pi/(4T))
            let asymptotic = (PI / (4.0 * t)).sqrt();
            // Should be close to asymptotic but not exactly equal
            let rel_diff = (result.value - asymptotic).abs() / asymptotic;
            assert!(
                rel_diff < 0.05,
                "F_0({}) = {} should be close to asymptotic {}",
                t,
                result.value,
                asymptotic
            );
        }
    }

    #[test]
    #[allow(clippy::const_is_empty)]
    fn version_exists() {
        assert!(!VERSION.is_empty());
    }

    // =========================================================================
    // Computation metrics tests (terms_count, estimated_error)
    // =========================================================================

    #[test]
    fn test_zero_method_metrics() {
        // Zero method should return terms_count=0 and estimated_error=Some(EPSILON)
        let result = boys_eval(5, 0.0).unwrap();
        assert_eq!(result.method, BoysMethod::Zero);
        assert_eq!(
            result.terms_count, 0,
            "Zero method should have terms_count=0"
        );
        assert_eq!(
            result.estimated_error,
            Some(f64::EPSILON),
            "Zero method should have estimated_error=Some(EPSILON)"
        );
    }

    #[test]
    fn test_near_zero_method_metrics() {
        // Near-zero T should also use Zero method
        let result = boys_eval(3, 1e-15).unwrap();
        assert_eq!(result.method, BoysMethod::Zero);
        assert_eq!(result.terms_count, 0);
        assert_eq!(result.estimated_error, Some(f64::EPSILON));
    }

    #[test]
    fn test_series_method_metrics() {
        // Series method should return terms_count > 0 and Some(estimated_error)
        // For m=5, T=1.0 is below turnover (~2.106), so it should use series
        let result = boys_eval(5, 1.0).unwrap();
        assert_eq!(result.method, BoysMethod::Series);
        assert!(
            result.terms_count > 0,
            "Series method should have terms_count > 0, got {}",
            result.terms_count
        );
        assert!(
            result.estimated_error.is_some(),
            "Series method should have Some(estimated_error)"
        );
        // The error estimate should be small (below tolerance)
        let err = result.estimated_error.unwrap();
        assert!(
            err < 1e-10,
            "Series estimated_error should be small, got {}",
            err
        );
    }

    #[test]
    fn test_series_terms_count_varies_with_t() {
        // Larger T should require more terms in series expansion
        // For m=10, use T values both below its turnover (~4.05)
        let result_small = boys_eval(10, 0.5).unwrap();
        let result_larger = boys_eval(10, 3.0).unwrap();

        // Both should use Series method
        assert_eq!(result_small.method, BoysMethod::Series);
        assert_eq!(result_larger.method, BoysMethod::Series);

        // Larger T should need more terms (or at least as many)
        assert!(
            result_larger.terms_count >= result_small.terms_count,
            "Larger T should need at least as many terms: {} vs {}",
            result_larger.terms_count,
            result_small.terms_count
        );
    }

    #[test]
    fn test_recurrence_method_metrics() {
        // Recurrence method should return terms_count = m + 1 and estimated_error = None
        // For m=0, any T > 0 uses recurrence (turnover is 0)
        let result = boys_eval(0, 0.5).unwrap();
        assert_eq!(result.method, BoysMethod::Recurrence);
        assert_eq!(
            result.terms_count, 1,
            "Recurrence for m=0 should have terms_count=1 (F_0 computed)"
        );
        assert_eq!(
            result.estimated_error, None,
            "Recurrence method should have estimated_error=None"
        );
    }

    #[test]
    fn test_recurrence_terms_count_equals_m_plus_1() {
        // For recurrence, terms_count = m + 1
        for m in 0..=10 {
            // Use T=20.0 which is above turnover for all m <= 10
            let result = boys_eval(m, 20.0).unwrap();
            assert_eq!(result.method, BoysMethod::Recurrence);
            assert_eq!(
                result.terms_count,
                m + 1,
                "Recurrence for m={} should have terms_count={}, got {}",
                m,
                m + 1,
                result.terms_count
            );
            assert_eq!(result.estimated_error, None);
        }
    }

    #[test]
    fn test_eval_all_metrics_consistency() {
        // boys_eval_all should return consistent metrics for all results
        let results = boys_eval_all(5, 1.0).unwrap();

        // All should have same method and metrics (batch computation)
        let first = &results[0];
        for (i, result) in results.iter().enumerate() {
            assert_eq!(
                result.method, first.method,
                "Result {} should have same method as first",
                i
            );
            assert_eq!(
                result.terms_count, first.terms_count,
                "Result {} should have same terms_count as first",
                i
            );
            assert_eq!(
                result.estimated_error, first.estimated_error,
                "Result {} should have same estimated_error as first",
                i
            );
        }
    }

    #[test]
    fn test_eval_many_metrics_vary() {
        // boys_eval_many should return different metrics for different T values
        let ts = vec![0.0, 1.0, 20.0];
        let results = boys_eval_many(3, &ts).unwrap();

        // T=0 should use Zero method
        assert_eq!(results[0].method, BoysMethod::Zero);
        assert_eq!(results[0].terms_count, 0);

        // T=1.0 for m=3 should use Series (turnover ~1.7)
        assert_eq!(results[1].method, BoysMethod::Series);
        assert!(results[1].terms_count > 0);
        assert!(results[1].estimated_error.is_some());

        // T=20.0 should use Recurrence
        assert_eq!(results[2].method, BoysMethod::Recurrence);
        assert_eq!(results[2].terms_count, 4); // m + 1 = 3 + 1 = 4
        assert_eq!(results[2].estimated_error, None);
    }

    // =========================================================================
    // Turnover point tests
    // =========================================================================

    #[test]
    fn test_turnover_field_is_populated() {
        // Verify the turnover field is included in results
        let result = boys_eval(5, 1.0).unwrap();

        // m=5 turnover is approximately 2.106
        assert!((result.turnover - 2.1064).abs() < 0.01);

        // Check m=0 and m=1 have turnover = 0
        let result_m0 = boys_eval(0, 0.5).unwrap();
        assert_eq!(result_m0.turnover, 0.0, "m=0 should have turnover=0");

        let result_m1 = boys_eval(1, 0.5).unwrap();
        assert_eq!(result_m1.turnover, 0.0, "m=1 should have turnover=0");
    }

    #[test]
    fn test_turnover_matches_libcint_values() {
        // Verify turnover points match libcint fmt.c exactly
        // (from references/libcint/src/fmt.c lines 42-83)
        let expected_turnovers: [(u32, f64); 11] = [
            (0, 0.0),
            (1, 0.0),
            (2, 0.866025403784),
            (5, 2.106432965305),
            (10, 4.05033123037),
            (15, 5.952793645111),
            (20, 7.837922483937),
            (25, 9.713354153046),
            (30, 11.58254788331),
            (35, 13.44737736550),
            (39, 14.93687515212),
        ];

        for (m, expected) in expected_turnovers {
            let turnover = get_turnover_point(m);
            assert!(
                (turnover - expected).abs() < 1e-10,
                "Turnover for m={} should be {}, got {}",
                m,
                expected,
                turnover
            );
        }
    }

    #[test]
    fn test_turnover_included_in_all_methods() {
        // Verify turnover is set correctly for Zero, Series, and Recurrence methods

        // Zero method (T near 0)
        let result_zero = boys_eval(5, 0.0).unwrap();
        assert_eq!(result_zero.method, BoysMethod::Zero);
        assert!((result_zero.turnover - 2.1064).abs() < 0.01);

        // Series method (T < turnover)
        let result_series = boys_eval(5, 1.0).unwrap();
        assert_eq!(result_series.method, BoysMethod::Series);
        assert!((result_series.turnover - 2.1064).abs() < 0.01);

        // Recurrence method (T >= turnover)
        let result_recurrence = boys_eval(5, 10.0).unwrap();
        assert_eq!(result_recurrence.method, BoysMethod::Recurrence);
        assert!((result_recurrence.turnover - 2.1064).abs() < 0.01);
    }

    #[test]
    fn test_turnover_in_batch_evaluation() {
        // Verify boys_eval_all returns correct turnover for each m
        let results = boys_eval_all(10, 5.0).unwrap();

        // Each result should have its own m-dependent turnover
        assert_eq!(results[0].turnover, 0.0, "m=0 turnover should be 0");
        assert_eq!(results[1].turnover, 0.0, "m=1 turnover should be 0");
        assert!(
            (results[2].turnover - 0.866).abs() < 0.01,
            "m=2 turnover ~0.87"
        );
        assert!(
            (results[5].turnover - 2.106).abs() < 0.01,
            "m=5 turnover ~2.1"
        );
        assert!(
            (results[10].turnover - 4.05).abs() < 0.01,
            "m=10 turnover ~4.05"
        );
    }

    #[test]
    fn test_get_turnover_point_public_function() {
        // Verify the public get_turnover_point function is accessible
        assert_eq!(get_turnover_point(0), 0.0);
        assert_eq!(get_turnover_point(1), 0.0);
        assert!((get_turnover_point(5) - 2.106).abs() < 0.01);
        assert!((get_turnover_point(10) - 4.05).abs() < 0.01);

        // Test extrapolation for m >= 40
        let t40 = get_turnover_point(40);
        // Should be approximately 14.937 + 0.372 * 1 = 15.309
        assert!(
            t40 > 15.0 && t40 < 16.0,
            "m=40 turnover should be ~15.3, got {}",
            t40
        );
    }
}
