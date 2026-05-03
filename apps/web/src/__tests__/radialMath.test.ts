/**
 * Unit tests for radialMath.ts -- client-side radial profile evaluation.
 *
 * Cross-validates the TypeScript evaluation against golden values from
 * the Rust `crates/qc-core/src/basis/radial.rs` tests, which were themselves
 * validated against PySCF 2.11.0 eval_gto output.
 *
 * @module __tests__/radialMath.test
 */

import { describe, it, expect } from 'vitest';
import {
  evaluateRadialProfileTS,
  computeOverlapMetric,
} from '../components/basisExplorer/radialMath';

// =============================================================================
// Reference data (from PySCF 2.11.0 via Rust golden tests)
// =============================================================================

/** H 1s STO-3G exponents (PySCF 2.11.0) */
const H_1S_EXPONENTS = [3.4252509100, 0.6239137300, 0.1688554000];
/** H 1s STO-3G raw contraction coefficients (PySCF 2.11.0) */
const H_1S_COEFFICIENTS = [0.1543289707, 0.5353281424, 0.4446345420];

/** O 2p STO-3G exponents (PySCF 2.11.0) */
const O_2P_EXPONENTS = [5.0331513000, 1.1695961000, 0.3803890000];
/** O 2p STO-3G raw contraction coefficients (PySCF 2.11.0) */
const O_2P_COEFFICIENTS = [0.1559162685, 0.6076837141, 0.3919573862];

// =============================================================================
// RM: evaluateRadialProfileTS tests
// =============================================================================

describe('evaluateRadialProfileTS', () => {
  // RM-1: TS evaluation matches WASM for H 1s STO-3G
  it('RM-1: matches Rust golden values for H 1s STO-3G', () => {
    const rValues = [0.0, 0.5, 1.0, 2.0];
    const result = evaluateRadialProfileTS(
      H_1S_EXPONENTS,
      H_1S_COEFFICIENTS,
      0,
      rValues,
    );

    // Golden values from PySCF eval_gto (validated in Rust tests)
    const expected = [
      6.282468823221000e-01,  // r=0.0
      4.268012209132965e-01,  // r=0.5
      2.230358271219306e-01,  // r=1.0
      6.456483925938930e-02,  // r=2.0
    ];

    for (let i = 0; i < rValues.length; i++) {
      expect(Math.abs(result[i] - expected[i])).toBeLessThan(1e-10);
    }
  });

  // RM-2: TS evaluation matches WASM for O 2p STO-3G
  it('RM-2: matches Rust golden values for O 2p STO-3G', () => {
    const rValues = [0.0, 0.5, 1.0, 1.5];
    const result = evaluateRadialProfileTS(
      O_2P_EXPONENTS,
      O_2P_COEFFICIENTS,
      1,
      rValues,
    );

    // Golden values from PySCF eval_gto (validated in Rust tests)
    const expected = [
      0.0,                    // r=0.0, p-type is zero at origin
      7.071411779390142e-01,  // r=0.5
      4.521398203032004e-01,  // r=1.0
      2.201369187979187e-01,  // r=1.5
    ];

    for (let i = 0; i < rValues.length; i++) {
      expect(Math.abs(result[i] - expected[i])).toBeLessThan(1e-10);
    }
  });

  // RM-3: s-type R(0) equals sum of effective coefficients
  it('RM-3: s-type R(0) equals sum of effective coefficients', () => {
    const rValues = [0.0];
    const result = evaluateRadialProfileTS(
      H_1S_EXPONENTS,
      H_1S_COEFFICIENTS,
      0,
      rValues,
    );

    // At r=0 for s-type: R(0) = sum_k c_k * N(alpha_k) * exp(0) = sum_k d_k
    // This is 6.282468823221000e-01 from PySCF
    expect(result[0]).toBeCloseTo(6.282468823221000e-01, 10);
  });

  // RM-4: p-type R(0) equals zero
  it('RM-4: p-type R(0) is exactly zero', () => {
    const rValues = [0.0];
    const result = evaluateRadialProfileTS(
      O_2P_EXPONENTS,
      O_2P_COEFFICIENTS,
      1,
      rValues,
    );
    expect(Math.abs(result[0])).toBeLessThan(1e-15);
  });

  // RM-5: d-type R(0) equals zero
  it('RM-5: d-type R(0) is exactly zero', () => {
    const rValues = [0.0];
    const result = evaluateRadialProfileTS(
      [5.0, 1.0, 0.3],
      [0.2, 0.5, 0.3],
      2,
      rValues,
    );
    expect(Math.abs(result[0])).toBeLessThan(1e-15);
  });

  // RM-6: Modified exponents produce different profile
  it('RM-6: modified exponents change the profile', () => {
    const rValues = [0.0, 0.5, 1.0, 2.0, 3.0];

    const original = evaluateRadialProfileTS(
      H_1S_EXPONENTS,
      H_1S_COEFFICIENTS,
      0,
      rValues,
    );

    // Double the first exponent
    const modifiedExponents = [H_1S_EXPONENTS[0] * 2, H_1S_EXPONENTS[1], H_1S_EXPONENTS[2]];
    const modified = evaluateRadialProfileTS(
      modifiedExponents,
      H_1S_COEFFICIENTS,
      0,
      rValues,
    );

    // At least one point must differ by > 1e-6
    let maxDiff = 0;
    for (let i = 0; i < rValues.length; i++) {
      maxDiff = Math.max(maxDiff, Math.abs(original[i] - modified[i]));
    }
    expect(maxDiff).toBeGreaterThan(1e-6);
  });

  // RM-7: 200-point evaluation completes in < 5 ms
  it('RM-7: 200-point evaluation completes in < 5 ms', () => {
    // Generate 200 r-points
    const nPts = 200;
    const rMax = 12.0;
    const rValues = Array.from({ length: nPts }, (_, i) => (rMax * i) / (nPts - 1));

    const start = performance.now();
    evaluateRadialProfileTS(H_1S_EXPONENTS, H_1S_COEFFICIENTS, 0, rValues);
    const elapsed = performance.now() - start;

    expect(elapsed).toBeLessThan(5.0);
  });

  // Additional: contracted equals sum of individual primitives at each r
  it('contracted profile equals sum of individual primitives', () => {
    const rValues = [0.0, 0.5, 1.0, 2.0, 5.0];
    const contracted = evaluateRadialProfileTS(
      H_1S_EXPONENTS,
      H_1S_COEFFICIENTS,
      0,
      rValues,
    );

    // Evaluate each primitive separately and sum
    for (let i = 0; i < rValues.length; i++) {
      let sum = 0;
      for (let k = 0; k < H_1S_EXPONENTS.length; k++) {
        const singlePrim = evaluateRadialProfileTS(
          [H_1S_EXPONENTS[k]],
          [H_1S_COEFFICIENTS[k]],
          0,
          [rValues[i]],
        );
        sum += singlePrim[0];
      }
      expect(Math.abs(contracted[i] - sum)).toBeLessThan(1e-14);
    }
  });

  // Additional: output length matches input
  it('output length matches rValues length', () => {
    const rValues = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    const result = evaluateRadialProfileTS(
      H_1S_EXPONENTS,
      H_1S_COEFFICIENTS,
      0,
      rValues,
    );
    expect(result.length).toBe(rValues.length);
  });
});

// =============================================================================
// OV: computeOverlapMetric tests
// =============================================================================

describe('computeOverlapMetric', () => {
  /** Generate a uniform r grid */
  function rGrid(rMax: number, n: number): number[] {
    return Array.from({ length: n }, (_, i) => (rMax * i) / (n - 1));
  }

  // OV-1: Identical profiles yield S=1.0
  it('OV-1: identical profiles yield S = 1.0', () => {
    const r = rGrid(10.0, 200);
    const profile = evaluateRadialProfileTS(
      H_1S_EXPONENTS,
      H_1S_COEFFICIENTS,
      0,
      r,
    );
    const S = computeOverlapMetric(profile, profile, r);
    expect(Math.abs(S - 1.0)).toBeLessThan(1e-10);
  });

  // OV-2: Orthogonal profiles yield S ~ 0.0
  it('OV-2: non-overlapping narrow Gaussians yield S near 0', () => {
    const r = rGrid(20.0, 1000);
    // Two narrow Gaussians centered far apart (in terms of peak location)
    // One very tight (large alpha, peaked near 0), one very diffuse (small alpha, peaked far)
    const profile1 = evaluateRadialProfileTS([1000.0], [1.0], 0, r);
    const profile2 = evaluateRadialProfileTS([0.001], [1.0], 0, r);
    const S = computeOverlapMetric(profile1, profile2, r);
    // These won't be perfectly orthogonal on a finite grid, but overlap should be very low
    expect(Math.abs(S)).toBeLessThan(0.1);
  });

  // OV-3: Slightly modified profile yields S close to 1.0
  it('OV-3: 1% exponent change gives S > 0.99', () => {
    const r = rGrid(12.0, 200);
    const original = evaluateRadialProfileTS(
      H_1S_EXPONENTS,
      H_1S_COEFFICIENTS,
      0,
      r,
    );

    // Change first exponent by 1%
    const modExponents = [
      H_1S_EXPONENTS[0] * 1.01,
      H_1S_EXPONENTS[1],
      H_1S_EXPONENTS[2],
    ];
    const modified = evaluateRadialProfileTS(
      modExponents,
      H_1S_COEFFICIENTS,
      0,
      r,
    );

    const S = computeOverlapMetric(original, modified, r);
    expect(S).toBeGreaterThan(0.99);
    expect(S).toBeLessThanOrEqual(1.0);
  });

  // OV-4: S is symmetric
  it('OV-4: overlap is symmetric (S(a,b) = S(b,a))', () => {
    const r = rGrid(12.0, 200);
    const profile1 = evaluateRadialProfileTS(
      H_1S_EXPONENTS,
      H_1S_COEFFICIENTS,
      0,
      r,
    );

    const modExponents = [
      H_1S_EXPONENTS[0] * 2.0,
      H_1S_EXPONENTS[1],
      H_1S_EXPONENTS[2],
    ];
    const profile2 = evaluateRadialProfileTS(
      modExponents,
      H_1S_COEFFICIENTS,
      0,
      r,
    );

    const S12 = computeOverlapMetric(profile1, profile2, r);
    const S21 = computeOverlapMetric(profile2, profile1, r);
    expect(Math.abs(S12 - S21)).toBeLessThan(1e-14);
  });

  // OV-5: S is bounded in [-1, 1]
  it('OV-5: overlap is bounded in [-1, 1]', () => {
    const r = rGrid(12.0, 200);
    const profile1 = evaluateRadialProfileTS(
      H_1S_EXPONENTS,
      H_1S_COEFFICIENTS,
      0,
      r,
    );

    // Large exponent modification
    const modExponents = [
      H_1S_EXPONENTS[0] * 5.0,
      H_1S_EXPONENTS[1] * 0.3,
      H_1S_EXPONENTS[2] * 8.0,
    ];
    const profile2 = evaluateRadialProfileTS(
      modExponents,
      H_1S_COEFFICIENTS,
      0,
      r,
    );

    const S = computeOverlapMetric(profile1, profile2, r);
    expect(S).toBeGreaterThanOrEqual(-1.0);
    expect(S).toBeLessThanOrEqual(1.0 + 1e-10);
  });

  // OV-6: Trapezoidal rule accuracy for analytical overlap
  it('OV-6: trapezoidal overlap matches analytical within 1%', () => {
    // For two s-type single-primitive Gaussians:
    //   integral_0^inf exp(-a*r^2) * exp(-b*r^2) dr = sqrt(pi/(a+b))/2
    // The normalized overlap:
    //   S = integral / sqrt(integral_a * integral_b)
    //   integral_a = sqrt(pi/(2a))/2, integral_b = sqrt(pi/(2b))/2
    //   S = [sqrt(pi/(a+b))/2] / sqrt([sqrt(pi/(2a))/2] * [sqrt(pi/(2b))/2])
    //   S = sqrt(pi/(a+b)) / (pi/(2*sqrt(a*b)))^(1/2) ... simplifying:
    //   S = [2*sqrt(a*b)/(a+b)]^(1/2)

    const a = 1.0;
    const b = 2.0;

    // Analytical normalized overlap of two bare Gaussians (without normalization prefactors)
    const S_analytical = Math.sqrt(2.0 * Math.sqrt(a * b) / (a + b));

    // Numerical: evaluate profiles on r grid and compute overlap
    const rMax = 10.0;
    const n = 500;
    const r = rGrid(rMax, n);

    // Use the raw Gaussians exp(-a*r^2), not normalized, to match the formula
    const profileA = r.map(ri => Math.exp(-a * ri * ri));
    const profileB = r.map(ri => Math.exp(-b * ri * ri));

    const S_numerical = computeOverlapMetric(profileA, profileB, r);

    const relError = Math.abs(S_numerical - S_analytical) / S_analytical;
    expect(relError).toBeLessThan(0.01); // Within 1%
  });

  // Additional: two-point grid edge case
  it('handles minimum grid size (2 points)', () => {
    const r = [0.0, 1.0];
    const p1 = [1.0, 0.5];
    const p2 = [1.0, 0.5];
    const S = computeOverlapMetric(p1, p2, r);
    expect(Math.abs(S - 1.0)).toBeLessThan(1e-10);
  });

  // Additional: single-point returns 1.0
  it('single-point grid returns 1.0', () => {
    const S = computeOverlapMetric([1.0], [2.0], [0.0]);
    expect(S).toBe(1.0);
  });
});
