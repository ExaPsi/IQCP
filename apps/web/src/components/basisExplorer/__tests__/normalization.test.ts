/**
 * Unit tests for primitiveNormalization utility.
 *
 * Validates that the TypeScript implementation matches the Rust
 * `primitive_normalization` in `crates/qc-core/src/basis/radial.rs`.
 *
 * @module basisExplorer/__tests__/normalization
 */

import { describe, it, expect } from 'vitest';
import { primitiveNormalization } from '../normalization';

// =============================================================================
// s-type normalization
// =============================================================================

describe('primitiveNormalization - s-type (l=0)', () => {
  it('computes N = (2*alpha/pi)^(3/4) for alpha=1.0', () => {
    const alpha = 1.0;
    const expected = Math.pow((2.0 * alpha) / Math.PI, 0.75);
    expect(primitiveNormalization(alpha, 0)).toBeCloseTo(expected, 14);
  });

  it('computes correctly for H 1s STO-3G first primitive (alpha=3.4252509100)', () => {
    const alpha = 3.42525091;
    const expected = Math.pow((2.0 * alpha) / Math.PI, 0.75);
    expect(primitiveNormalization(alpha, 0)).toBeCloseTo(expected, 12);
  });
});

// =============================================================================
// p-type normalization
// =============================================================================

describe('primitiveNormalization - p-type (l=1)', () => {
  it('computes N = (2*alpha/pi)^(3/4) * sqrt(4*alpha) for alpha=1.0', () => {
    const alpha = 1.0;
    const expected = Math.pow((2.0 * alpha) / Math.PI, 0.75) * Math.sqrt(4.0 * alpha);
    expect(primitiveNormalization(alpha, 1)).toBeCloseTo(expected, 14);
  });

  it('computes correctly for a range of exponents', () => {
    for (const alpha of [0.5, 2.0, 10.0]) {
      const expected = Math.pow((2.0 * alpha) / Math.PI, 0.75) * Math.sqrt(4.0 * alpha);
      expect(primitiveNormalization(alpha, 1)).toBeCloseTo(expected, 12);
    }
  });
});

// =============================================================================
// d-type normalization (d_xy representative)
// =============================================================================

describe('primitiveNormalization - d-type (l=2)', () => {
  it('computes N = (2*alpha/pi)^(3/4) * 4*alpha for alpha=1.0', () => {
    const alpha = 1.0;
    const expected = Math.pow((2.0 * alpha) / Math.PI, 0.75) * 4.0 * alpha;
    expect(primitiveNormalization(alpha, 2)).toBeCloseTo(expected, 14);
  });

  it('computes correctly for typical d-type exponent (alpha=0.8)', () => {
    const alpha = 0.8;
    const expected = Math.pow((2.0 * alpha) / Math.PI, 0.75) * 4.0 * alpha;
    expect(primitiveNormalization(alpha, 2)).toBeCloseTo(expected, 12);
  });
});

// =============================================================================
// Positivity
// =============================================================================

describe('primitiveNormalization - positivity', () => {
  it('returns positive values for all angular momenta and alpha > 0', () => {
    for (const alpha of [0.1, 0.5, 1.0, 5.0, 50.0, 130.0]) {
      expect(primitiveNormalization(alpha, 0)).toBeGreaterThan(0);
      expect(primitiveNormalization(alpha, 1)).toBeGreaterThan(0);
      expect(primitiveNormalization(alpha, 2)).toBeGreaterThan(0);
    }
  });
});

// =============================================================================
// Ratio tests (analytical relationships between angular momenta)
// =============================================================================

describe('primitiveNormalization - ratio tests', () => {
  it('N_p/N_s = 2*sqrt(alpha)', () => {
    for (const alpha of [0.5, 1.0, 5.0]) {
      const ns = primitiveNormalization(alpha, 0);
      const np = primitiveNormalization(alpha, 1);
      expect(np / ns).toBeCloseTo(2 * Math.sqrt(alpha), 12);
    }
  });

  it('N_d/N_s = 4*alpha', () => {
    for (const alpha of [0.5, 1.0, 5.0]) {
      const ns = primitiveNormalization(alpha, 0);
      const nd = primitiveNormalization(alpha, 2);
      expect(nd / ns).toBeCloseTo(4 * alpha, 12);
    }
  });

  it('N_d/N_p = 2*sqrt(alpha)', () => {
    for (const alpha of [0.5, 1.0, 5.0]) {
      const np = primitiveNormalization(alpha, 1);
      const nd = primitiveNormalization(alpha, 2);
      expect(nd / np).toBeCloseTo(2 * Math.sqrt(alpha), 12);
    }
  });
});

// =============================================================================
// Effective coefficient identity
// =============================================================================

describe('primitiveNormalization - effective coefficient', () => {
  it('d_k = c_k * N_k identity holds for H 1s STO-3G', () => {
    const exponents = [3.42525091, 0.62391373, 0.1688554];
    const coefficients = [0.1543289707, 0.5353281424, 0.4446345420];

    for (let i = 0; i < exponents.length; i++) {
      const nk = primitiveNormalization(exponents[i], 0);
      const dk = coefficients[i] * nk;
      expect(dk).toBeCloseTo(coefficients[i] * nk, 12);
      expect(nk).toBeGreaterThan(0);
      expect(dk).toBeGreaterThan(0);
    }
  });
});

// =============================================================================
// Default/unsupported angular momentum
// =============================================================================

describe('primitiveNormalization - unsupported angular momentum', () => {
  it('falls back to base normalization for l > 2', () => {
    const alpha = 1.0;
    const base = Math.pow((2.0 * alpha) / Math.PI, 0.75);
    expect(primitiveNormalization(alpha, 3)).toBeCloseTo(base, 14);
    expect(primitiveNormalization(alpha, 4)).toBeCloseTo(base, 14);
  });
});
