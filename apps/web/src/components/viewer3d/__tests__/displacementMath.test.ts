/**
 * Unit tests for DisplacementArrows helper functions (US-102).
 *
 * Pure-logic tests for the exported helpers in `displacementMath.ts`:
 *   - displacementFieldNorm
 *   - normalizeDisplacement
 *   - animationPhase
 *   - displacementEndpoint
 *   - wrapModeIndex
 *
 * The `<DisplacementArrows>` component itself (R3F + Three.js) cannot be
 * unit-mounted in a node-test environment without WebGL. It is exercised
 * via the manual smoke test documented in VP-US-102 Section 17.
 *
 * @module components/viewer3d/__tests__/displacementMath.test
 */

import { describe, it, expect } from 'vitest';
import {
  displacementFieldNorm,
  normalizeDisplacement,
  animationPhase,
  displacementEndpoint,
  wrapModeIndex,
  type DisplacementVec,
} from '../displacementMath';

// ============================================================================
// displacementFieldNorm
// ============================================================================

describe('displacementFieldNorm', () => {
  it('returns 0 for an empty field', () => {
    expect(displacementFieldNorm([])).toBe(0);
  });

  it('returns 0 for a fully-zero field', () => {
    expect(
      displacementFieldNorm([
        [0, 0, 0],
        [0, 0, 0],
      ])
    ).toBe(0);
  });

  it('computes the L2 norm of a single-atom vector', () => {
    expect(displacementFieldNorm([[3, 4, 0]])).toBeCloseTo(5, 15);
  });

  it('concatenates across atoms (H₂ antisymmetric stretch)', () => {
    // Two atoms with opposite unit displacements along x → total norm = √2
    expect(
      displacementFieldNorm([
        [1, 0, 0],
        [-1, 0, 0],
      ])
    ).toBeCloseTo(Math.sqrt(2), 15);
  });
});

// ============================================================================
// normalizeDisplacement
// ============================================================================

describe('normalizeDisplacement', () => {
  it('returns an all-zero field unchanged when the input is zero', () => {
    const input: DisplacementVec[] = [
      [0, 0, 0],
      [0, 0, 0],
    ];
    const result = normalizeDisplacement(input);
    expect(result).toEqual(input);
  });

  it('rescales to unit L2 norm', () => {
    const input: DisplacementVec[] = [
      [3, 4, 0],
      [0, 0, 0],
    ];
    const out = normalizeDisplacement(input);
    expect(displacementFieldNorm(out)).toBeCloseTo(1, 14);
  });

  it('preserves the relative direction of each atom', () => {
    const input: DisplacementVec[] = [
      [1, 0, 0],
      [-1, 0, 0],
    ];
    const out = normalizeDisplacement(input);
    // Original ratio: atom 1 is opposite to atom 0. Must still hold after scaling.
    expect(out[1][0]).toBeCloseTo(-out[0][0], 15);
  });

  it('does not mutate the input array', () => {
    const input: DisplacementVec[] = [
      [3, 4, 0],
      [1, 2, 2],
    ];
    const snapshot = JSON.parse(JSON.stringify(input));
    normalizeDisplacement(input);
    expect(input).toEqual(snapshot);
  });

  it('already-normalized input is stable (idempotent up to rounding)', () => {
    const input: DisplacementVec[] = [
      [1 / Math.sqrt(2), 1 / Math.sqrt(2), 0],
    ];
    const out = normalizeDisplacement(input);
    expect(out[0][0]).toBeCloseTo(input[0][0], 15);
  });
});

// ============================================================================
// animationPhase
// ============================================================================

describe('animationPhase', () => {
  it('returns 0 at t = 0', () => {
    expect(animationPhase(0, 1)).toBe(0);
  });

  it('returns 1 at t = 0.25 with speed = 1 (quarter period)', () => {
    expect(animationPhase(0.25, 1)).toBeCloseTo(1, 14);
  });

  it('returns 0 again at t = 0.5 with speed = 1 (half period)', () => {
    expect(animationPhase(0.5, 1)).toBeCloseTo(0, 14);
  });

  it('returns -1 at t = 0.75 with speed = 1 (three-quarter period)', () => {
    expect(animationPhase(0.75, 1)).toBeCloseTo(-1, 14);
  });

  it('is bounded in [-1, 1]', () => {
    for (const t of [0, 0.1, 0.5, 1.3, 42, 137.5]) {
      const p = animationPhase(t, 1.5);
      expect(p).toBeGreaterThanOrEqual(-1);
      expect(p).toBeLessThanOrEqual(1);
    }
  });

  it('doubling speed halves the period', () => {
    // At speed = 2, the quarter-period (peak) is at t = 0.125
    expect(animationPhase(0.125, 2)).toBeCloseTo(1, 14);
  });
});

// ============================================================================
// displacementEndpoint
// ============================================================================

describe('displacementEndpoint', () => {
  it('returns start when phase = 0', () => {
    expect(displacementEndpoint([1, 2, 3], [1, 0, 0], 0.5, 0)).toEqual([
      1, 2, 3,
    ]);
  });

  it('returns start + A·d when phase = 1', () => {
    expect(displacementEndpoint([0, 0, 0], [1, 0, 0], 0.5, 1)).toEqual([
      0.5, 0, 0,
    ]);
  });

  it('returns start − A·d when phase = -1', () => {
    expect(displacementEndpoint([0, 0, 0], [1, 0, 0], 0.5, -1)).toEqual([
      -0.5, 0, 0,
    ]);
  });

  it('scales linearly with amplitude', () => {
    const e1 = displacementEndpoint([0, 0, 0], [0, 1, 0], 1.0, 1);
    const e2 = displacementEndpoint([0, 0, 0], [0, 1, 0], 2.0, 1);
    expect(e2[1]).toBeCloseTo(2 * e1[1], 15);
  });

  it('preserves the start position for a zero displacement', () => {
    expect(displacementEndpoint([5, 6, 7], [0, 0, 0], 1.0, 1)).toEqual([
      5, 6, 7,
    ]);
  });
});

// ============================================================================
// wrapModeIndex
// ============================================================================

describe('wrapModeIndex', () => {
  it('returns 0 for an empty mode list', () => {
    expect(wrapModeIndex(0, 0)).toBe(0);
    expect(wrapModeIndex(5, 0)).toBe(0);
    expect(wrapModeIndex(-1, 0)).toBe(0);
  });

  it('returns the index unchanged when in range', () => {
    expect(wrapModeIndex(0, 6)).toBe(0);
    expect(wrapModeIndex(3, 6)).toBe(3);
    expect(wrapModeIndex(5, 6)).toBe(5);
  });

  it('wraps positive overflow', () => {
    expect(wrapModeIndex(6, 6)).toBe(0);
    expect(wrapModeIndex(7, 6)).toBe(1);
    expect(wrapModeIndex(13, 6)).toBe(1);
  });

  it('wraps negative underflow', () => {
    expect(wrapModeIndex(-1, 6)).toBe(5);
    expect(wrapModeIndex(-7, 6)).toBe(5);
  });
});
