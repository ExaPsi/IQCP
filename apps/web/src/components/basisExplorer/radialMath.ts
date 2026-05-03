/**
 * Client-side radial profile evaluation for real-time exponent adjustment.
 *
 * This is a TypeScript port of the formula used in
 * `crates/qc-core/src/basis/radial.rs::evaluate_radial_profile`,
 * but evaluated on the client for instant (<100ms) slider response.
 *
 * Formula: R(r) = r^l * sum_k { c_k * N(alpha_k) * exp(-alpha_k * r^2) }
 * where N(alpha_k) is the primitive normalization from normalization.ts.
 *
 * Both functions are pure, with no dependencies on React, Zustand, or the worker.
 *
 * @module components/basisExplorer/radialMath
 */

import { primitiveNormalization } from './normalization';

/**
 * Evaluate a contracted radial profile at given r values.
 *
 * @param exponents - Primitive Gaussian exponents alpha_k
 * @param coefficients - Raw contraction coefficients c_k
 * @param angularMomentum - Shell angular momentum (0, 1, or 2)
 * @param rValues - Array of r values (in Bohr) at which to evaluate
 * @returns Array of R(r) values at each r point
 *
 * Preconditions:
 *   - exponents.length === coefficients.length
 *   - exponents.length >= 1
 *   - All exponents > 0
 *   - rValues.length >= 1
 *
 * Postconditions:
 *   - Return array length === rValues.length
 *   - For s-type with all positive coefficients: R(r) >= 0 for all r >= 0
 *   - R(0) = 0 for p-type (l=1) and d-type (l=2)
 */
export function evaluateRadialProfileTS(
  exponents: number[],
  coefficients: number[],
  angularMomentum: number,
  rValues: number[],
): number[] {
  const nPrims = exponents.length;
  const nPts = rValues.length;

  // Pre-compute effective coefficients: d_k = c_k * N(alpha_k)
  const effectiveCoeffs = new Array<number>(nPrims);
  for (let k = 0; k < nPrims; k++) {
    effectiveCoeffs[k] = coefficients[k] * primitiveNormalization(exponents[k], angularMomentum);
  }

  // Evaluate contracted profile at each r value
  const result = new Array<number>(nPts);
  for (let i = 0; i < nPts; i++) {
    const r = rValues[i];
    const r2 = r * r;

    // r^l prefactor
    let rL: number;
    switch (angularMomentum) {
      case 0:
        rL = 1.0;
        break;
      case 1:
        rL = r;
        break;
      case 2:
        rL = r2;
        break;
      default:
        rL = Math.pow(r, angularMomentum);
    }

    // Sum over primitives: sum_k { d_k * exp(-alpha_k * r^2) }
    let sum = 0.0;
    for (let k = 0; k < nPrims; k++) {
      sum += effectiveCoeffs[k] * Math.exp(-exponents[k] * r2);
    }

    result[i] = rL * sum;
  }

  return result;
}

/**
 * Compute the normalized overlap between two radial profiles.
 *
 * Uses trapezoidal rule over the r grid to compute:
 *   S = integral(R1 * R2 dr) / sqrt(integral(R1^2 dr) * integral(R2^2 dr))
 *
 * @param profile1 - First profile R(r) values
 * @param profile2 - Second profile R(r) values
 * @param rValues - r grid points (must be uniformly spaced)
 * @returns Normalized overlap S in [-1, 1]; typically [0, 1]
 *
 * Preconditions:
 *   - profile1.length === profile2.length === rValues.length
 *   - rValues.length >= 2
 *
 * Postconditions:
 *   - If profile1 === profile2 (element-wise), S = 1.0 (within floating-point tolerance)
 *   - |S| <= 1.0
 */
export function computeOverlapMetric(
  profile1: number[],
  profile2: number[],
  rValues: number[],
): number {
  const n = rValues.length;
  if (n < 2) return 1.0;

  // Uniform grid spacing
  const dr = rValues[1] - rValues[0];

  // Trapezoidal rule: integral ~ dr * (f[0]/2 + f[1] + f[2] + ... + f[n-2] + f[n-1]/2)
  let cross = 0.0;
  let norm1 = 0.0;
  let norm2 = 0.0;

  // First point (weight 0.5)
  cross += 0.5 * profile1[0] * profile2[0];
  norm1 += 0.5 * profile1[0] * profile1[0];
  norm2 += 0.5 * profile2[0] * profile2[0];

  // Interior points (weight 1.0)
  for (let i = 1; i < n - 1; i++) {
    cross += profile1[i] * profile2[i];
    norm1 += profile1[i] * profile1[i];
    norm2 += profile2[i] * profile2[i];
  }

  // Last point (weight 0.5)
  cross += 0.5 * profile1[n - 1] * profile2[n - 1];
  norm1 += 0.5 * profile1[n - 1] * profile1[n - 1];
  norm2 += 0.5 * profile2[n - 1] * profile2[n - 1];

  // Multiply by dr
  cross *= dr;
  norm1 *= dr;
  norm2 *= dr;

  // Normalized overlap
  const denominator = Math.sqrt(norm1 * norm2);
  if (denominator < 1e-30) return 1.0; // Both profiles are essentially zero

  return cross / denominator;
}
