/**
 * Boys Function Constants and Utilities
 *
 * Contains m-dependent turnover points and helper functions that match
 * the qc-core implementation. These are used for accurate regime
 * visualization in the frontend.
 *
 * IMPORTANT: There are TWO sets of thresholds:
 *
 * 1. THEORETICAL REGIMES (from lecture notes, Chapter 4):
 *    - Small T (T < 25): Series expansion
 *    - Moderate T (25 <= T < 30+5m): erf formula + upward recurrence
 *    - Large T (T >= 30+5m): Asymptotic expansion
 *
 * 2. IMPLEMENTATION REGIMES (from libcint, what qc-core actually uses):
 *    - Zero (T < 1e-14): Direct formula F_m(0) = 1/(2m+1)
 *    - Series (T < TURNOVER_POINT[m]): Power series + DOWNWARD recurrence
 *    - Recurrence (T >= TURNOVER_POINT[m]): erf + UPWARD recurrence
 *    (Note: Implementation does NOT have a separate asymptotic method)
 *
 * The key insight is that m=0 and m=1 have turnover=0, meaning they ALWAYS
 * use the recurrence method (erf-based) for any T > 0. Only m >= 2 can use
 * the series method, and only when T is below their specific turnover point.
 *
 * NUMERICAL STABILITY:
 * - Upward recurrence: F_{m+1} = [(2m+1)F_m - e^{-T}]/(2T)
 *   Division by 2T is problematic for small T (catastrophic cancellation)
 * - Downward recurrence: F_{m-1} = [2(m-1/2)F_m + e^{-T}]/(2T-1/(m-1/2))
 *   Stable for small T because terms reinforce rather than cancel
 *
 * Reference: crates/qc-core/src/boys/mod.rs (TURNOVER_POINT array)
 *
 * @module lib/boysConstants
 */

/**
 * Threshold below which T is treated as zero.
 * Matches qc-core NEAR_ZERO_THRESHOLD constant.
 */
export const NEAR_ZERO_THRESHOLD = 1e-14;

/**
 * M-dependent turnover points for regime selection (IMPLEMENTATION).
 *
 * For T < TURNOVER_POINTS[m], the series expansion + downward recurrence is used.
 * For T >= TURNOVER_POINTS[m], the erf + upward recurrence is used.
 *
 * CRITICAL: m=0 and m=1 have turnover=0, meaning:
 * - For m=0: F_0(T) is computed directly from erf(sqrt(T)), no recurrence needed
 * - For m=1: Only 1 recurrence step from F_0, very stable
 * - For m>=2 with small T: Series + downward recurrence is more stable
 *
 * These values are derived from error analysis ensuring optimal accuracy.
 * Reference: libcint fmt.c lines 42-83
 *
 * Formula: t0 = 0.5 * ((2m-1)!!/(2m-1)^0.5)^(1/(m-1))
 */
export const TURNOVER_POINTS: readonly number[] = [
  0.0,            // m=0: always use erf (no recurrence needed, direct formula)
  0.0,            // m=1: always use erf + 1 recurrence step (very stable)
  0.866025403784, // m=2: use series for T < 0.87, recurrence for T >= 0.87
  1.295010032056, // m=3
  1.705493613097, // m=4
  2.106432965305, // m=5
  2.501471934009, // m=6
  2.892473348218, // m=7
  3.280525047072, // m=8
  3.666320693281, // m=9
  4.05033123037,  // m=10
] as const;

/**
 * Get the turnover point for a given order m.
 *
 * For m values beyond the table (m > 10), extrapolates using
 * the approximate linear growth rate.
 *
 * @param m - Order of the Boys function
 * @returns The T value at which regime switches from series to recurrence
 */
export function getTurnoverPoint(m: number): number {
  if (m < 0) {
    return 0.0;
  }
  if (m < TURNOVER_POINTS.length) {
    return TURNOVER_POINTS[m];
  }
  // Extrapolate for large m (approximately linear growth ~0.38 per unit m)
  // Based on pattern in TURNOVER_POINTS: growth rate is ~0.38
  const lastIdx = TURNOVER_POINTS.length - 1;
  const lastValue = TURNOVER_POINTS[lastIdx];
  const mDiff = m - lastIdx;
  return lastValue + 0.38 * mDiff;
}

/**
 * Determine which computational method will be used for given (m, T).
 *
 * This function mirrors the exact logic in qc-core/src/boys/mod.rs.
 *
 * The three methods are:
 * 1. **zero**: T < 1e-14, uses direct formula F_m(0) = 1/(2m+1)
 * 2. **series**: T < TURNOVER_POINT[m], uses series expansion + downward recurrence
 * 3. **recurrence**: T >= TURNOVER_POINT[m], uses erf(sqrt(T)) + upward recurrence
 *
 * Key insight for m=0 and m=1:
 * - Their turnover points are 0.0, so for any T > 0, recurrence is always used
 * - This is numerically safe because:
 *   - m=0: F_0 is computed directly from erf, no recurrence at all
 *   - m=1: Only 1 recurrence step, very stable even for small T
 *
 * @param m - Order of the Boys function
 * @param T - Argument value
 * @returns The method: 'zero', 'series', or 'recurrence'
 */
export function getImplementationMethod(m: number, T: number): 'zero' | 'series' | 'recurrence' {
  if (T < NEAR_ZERO_THRESHOLD) {
    return 'zero';
  }

  const turnover = getTurnoverPoint(m);
  if (T < turnover) {
    return 'series';
  }
  return 'recurrence';
}

/**
 * Alias for getImplementationMethod for backward compatibility.
 * @deprecated Use getImplementationMethod instead for clarity.
 */
export function getBoysRegime(m: number, T: number): 'zero' | 'series' | 'recurrence' {
  return getImplementationMethod(m, T);
}

/**
 * Get the asymptotic threshold for a given order m.
 *
 * For T > 30 + 5*m, the asymptotic expansion achieves machine precision.
 * This is the theoretical "Large T" regime boundary from the lecture notes.
 *
 * Formula: T_asymptotic = 30 + 5 * m
 *
 * Reference: Lecture notes ch04_s06.tex, Table of regimes
 *
 * @param m - Order of the Boys function
 * @returns The T value above which asymptotic expansion is accurate to machine precision
 */
export function getAsymptoticThreshold(m: number): number {
  return 30 + 5 * m;
}

/**
 * Theoretical regime classification for educational purposes.
 *
 * The lecture notes describe 3 theoretical regimes:
 * - Small T (T < 25): Series expansion
 * - Moderate T (25 <= T < 30+5m): erf formula + upward recurrence
 * - Large T (T > 30+5m): Asymptotic expansion
 *
 * The implementation (following libcint) uses only 2 practical methods:
 * - Series: T < TURNOVER_POINT[m]
 * - Recurrence: T >= TURNOVER_POINT[m] (combines moderate + large T)
 *
 * @param m - Order of the Boys function
 * @param T - Argument value
 * @returns The theoretical regime: 'small', 'moderate', or 'large'
 */
export function getTheoreticalRegime(m: number, T: number): 'small' | 'moderate' | 'large' {
  const SMALL_T_THRESHOLD = 25;
  const asymptoticThreshold = getAsymptoticThreshold(m);

  if (T < SMALL_T_THRESHOLD) {
    return 'small';
  }
  if (T < asymptoticThreshold) {
    return 'moderate';
  }
  return 'large';
}

/**
 * Maximum T value displayed in charts.
 */
export const MAX_T_VALUE = 50;

/**
 * Maximum m value supported in the UI.
 */
export const MAX_M_VALUE = 10;

/**
 * Theoretical threshold where series expansion transitions to recurrence.
 * This is the "Small T" boundary from the lecture notes (Chapter 4).
 *
 * Theoretical Regimes:
 * - Small T: T < 25 (use series expansion)
 * - Moderate T: 25 <= T < 30+5m (use erf + upward recurrence)
 * - Large T: T >= 30+5m (use asymptotic expansion)
 *
 * Note: The actual implementation uses m-dependent turnover points from libcint,
 * which are much smaller than 25 (e.g., m=5 has turnover ~2.1).
 */
export const THEORETICAL_SMALL_T_THRESHOLD = 25;

/**
 * Theoretical regime formulas from lecture notes (Chapter 4).
 *
 * These are the PEDAGOGICAL formulas shown to students.
 * The implementation may use equivalent but different formulations.
 */
export const THEORETICAL_FORMULAS = {
  /**
   * Series expansion (Small T regime, T < 25)
   * F_m(T) = e^{-T} * sum_{k=0}^{inf} [(2m-1)!! * T^k / (2m+2k+1)!!]
   */
  series: 'F_m(T) = e^{-T} \\sum_{k=0}^{\\infty} \\frac{(2m-1)!! \\cdot T^k}{(2m+2k+1)!!}',

  /**
   * erf + Upward Recurrence (Moderate T regime, 25 <= T < 30+5m)
   * F_0(T) = (sqrt(pi) / (2*sqrt(T))) * erf(sqrt(T))
   * F_{m+1}(T) = [(2m+1)*F_m(T) - e^{-T}] / (2T)
   */
  recurrence: 'F_0(T) = \\frac{\\sqrt{\\pi}}{2\\sqrt{T}} \\cdot \\text{erf}(\\sqrt{T}), \\quad F_{m+1}(T) = \\frac{(2m+1) F_m(T) - e^{-T}}{2T}',

  /**
   * Asymptotic expansion (Large T regime, T >= 30+5m)
   * F_m(T) ~ [(2m-1)!! / 2^{m+1}] * sqrt(pi / T^{2m+1})
   */
  asymptotic: 'F_m(T) \\approx \\frac{(2m-1)!!}{2^{m+1}} \\sqrt{\\frac{\\pi}{T^{2m+1}}}',
} as const;

/**
 * Implementation formulas from libcint (what qc-core actually uses).
 *
 * These may differ from the theoretical formulas but are mathematically equivalent.
 */
export const IMPLEMENTATION_FORMULAS = {
  /**
   * Zero method (T near 0)
   * F_m(0) = 1 / (2m+1)
   */
  zero: 'F_m(0) = \\frac{1}{2m+1}',

  /**
   * Series expansion + Downward Recurrence (T < TURNOVER_POINT[m], m >= 2)
   *
   * First compute F_{m_high} using series:
   * F_m(T) = (e^{-T}/(2b)) * [1 + T/(b+1) + T^2/((b+1)(b+2)) + ...]
   * where b = m + 0.5
   *
   * Then use DOWNWARD recurrence to get lower orders:
   * F_{m-1}(T) = (2T * F_m(T) + e^{-T}) / (2m - 1)
   *
   * Downward recurrence is stable for small T because terms reinforce rather than cancel.
   */
  series: 'F_m(T) = \\frac{e^{-T}}{2b} \\left[ 1 + \\frac{T}{b+1} + \\cdots \\right], \\text{ then } F_{m-1} = \\frac{2T \\cdot F_m + e^{-T}}{2m-1}',

  /**
   * erf + Upward Recurrence (T >= TURNOVER_POINT[m])
   *
   * First compute F_0 from the error function:
   * F_0(T) = (sqrt(pi) / (2*sqrt(T))) * erf(sqrt(T))
   *
   * Then use UPWARD recurrence to get higher orders:
   * F_{m+1}(T) = [(2m+1)*F_m(T) - e^{-T}] / (2T)
   *
   * Upward recurrence is stable when T is large enough that subtraction doesn't
   * cause catastrophic cancellation (the turnover point ensures this).
   */
  recurrence: 'F_0(T) = \\frac{\\sqrt{\\pi}}{2\\sqrt{T}} \\cdot \\text{erf}(\\sqrt{T}), \\quad F_{m+1} = \\frac{(2m+1) F_m - e^{-T}}{2T}',

  /**
   * Special case explanations for m=0 and m=1
   */
  m0_explanation: 'F_0(T) computed directly from erf formula (no recurrence needed)',
  m1_explanation: 'F_1(T) = (1 - e^{-T}) / (2T) \\text{ (only 1 recurrence step from } F_0\\text{)}',
} as const;

/**
 * Detailed explanation of why each method is used.
 */
export const METHOD_EXPLANATIONS = {
  zero: {
    short: 'T is effectively zero',
    long: 'When T < 1e-14, the exponential e^{-T} is indistinguishable from 1, and the exact formula F_m(0) = 1/(2m+1) is used directly.',
  },
  series: {
    short: 'Series + downward recurrence for small T',
    long: 'For small T (below the m-dependent turnover), upward recurrence would suffer from catastrophic cancellation. Instead, we compute a high-order F_m using series expansion, then use numerically stable downward recurrence to get lower orders.',
  },
  recurrence: {
    short: 'erf + upward recurrence for moderate/large T',
    long: 'For T above the turnover point, we compute F_0 directly from erf(sqrt(T)), then use upward recurrence. This is stable because (2m+1)*F_m is much larger than e^{-T} when T is large enough.',
  },
  m0_special: {
    short: 'Direct erf formula (no recurrence)',
    long: 'For m=0, F_0(T) = sqrt(pi)/(2*sqrt(T)) * erf(sqrt(T)) is computed directly. No recurrence is needed, so there is no stability concern.',
  },
  m1_special: {
    short: 'Only 1 recurrence step (very stable)',
    long: 'For m=1, we compute F_0 from erf, then just one recurrence step to get F_1. This single step is stable even for relatively small T.',
  },
} as const;
