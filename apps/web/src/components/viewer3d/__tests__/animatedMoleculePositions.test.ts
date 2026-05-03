/**
 * Unit tests for AnimatedMolecule helper functions (US-102 UX refinement).
 *
 * Pure-logic tests for the exported helpers in `animatedMoleculePositions.ts`:
 *   - `computeDisplacedPositions`
 *   - `computeBondGeometry`
 *
 * The `<AnimatedMolecule>` component itself (R3F + Three.js) cannot be
 * unit-mounted in a node-test environment without WebGL. It is exercised
 * via the manual smoke test documented in the US-102 verification report.
 *
 * @module components/viewer3d/__tests__/animatedMoleculePositions.test
 */

import { describe, it, expect } from 'vitest';
import {
  computeDisplacedPositions,
  computeBondGeometry,
  type AnimAtom,
  type AnimDisplacement,
} from '../animatedMoleculePositions';

// Float32Array stores values in IEEE 754 single precision (~7 decimal digits).
// Use ~6 decimal digits of tolerance for any value read back from a Float32Array,
// while values that round-trip through Number (Float64) can use ~14.
const F32_DIGITS = 6;
const F64_DIGITS = 14;

// ============================================================================
// computeDisplacedPositions
// ============================================================================

describe('computeDisplacedPositions', () => {
  // ---- H₂ along z, antisymmetric stretch -----------------------------------
  // Two H atoms at z = ±0.7 bohr, displacement is +z on atom 1, -z on atom 0.
  const h2Atoms: AnimAtom[] = [
    [1, 0, 0, -0.7],
    [1, 0, 0, +0.7],
  ];
  const h2Stretch: AnimDisplacement[] = [
    [0, 0, -1],
    [0, 0, +1],
  ];

  // ---- H₂O bend mode ------------------------------------------------------
  // O at origin, two H atoms in the xz plane. Bend mode pushes both H atoms
  // perpendicular to their bonds (toward each other along ±x in this idealized
  // version).
  const h2oAtoms: AnimAtom[] = [
    [8, 0, 0, 0],
    [1, 0.78, 0, 0.59],
    [1, -0.78, 0, 0.59],
  ];
  const h2oBend: AnimDisplacement[] = [
    [0, 0, -0.05], // O moves slightly down to conserve momentum
    [-0.4, 0, 0.3], // H1 moves toward H2 and slightly up
    [+0.4, 0, 0.3], // H2 moves toward H1 and slightly up
  ];

  it('returns equilibrium positions when phase = 0 (within Float32 precision)', () => {
    // Note: -0.7 and 0.7 are not exactly representable in IEEE 754 single
    // precision, so the round-trip through Float32Array introduces a small
    // (~1e-8) error. The pattern is still "no displacement applied" — the
    // contract is that the equilibrium x/y/z are preserved up to F32 precision.
    const out = computeDisplacedPositions(h2Atoms, h2Stretch, 0.5, 0);
    expect(out[0]).toBe(0); // exact zeros remain exact
    expect(out[1]).toBe(0);
    expect(out[2]).toBeCloseTo(-0.7, F32_DIGITS); // z atom 0
    expect(out[3]).toBe(0);
    expect(out[4]).toBe(0);
    expect(out[5]).toBeCloseTo(0.7, F32_DIGITS); // z atom 1
  });

  it('returns equilibrium positions when amplitude = 0 (any phase)', () => {
    const out1 = computeDisplacedPositions(h2Atoms, h2Stretch, 0, 1);
    const out2 = computeDisplacedPositions(h2Atoms, h2Stretch, 0, -1);
    const out3 = computeDisplacedPositions(h2Atoms, h2Stretch, 0, 0.42);
    for (const out of [out1, out2, out3]) {
      expect(out[2]).toBeCloseTo(-0.7, F32_DIGITS);
      expect(out[5]).toBeCloseTo(0.7, F32_DIGITS);
    }
  });

  it('returns equilibrium + amplitude·displacement at phase = 1', () => {
    const out = computeDisplacedPositions(h2Atoms, h2Stretch, 0.3, 1);
    // Atom 0: z = -0.7 + 1·0.3·(-1) = -1.0
    expect(out[2]).toBeCloseTo(-1.0, F32_DIGITS);
    // Atom 1: z = +0.7 + 1·0.3·(+1) = +1.0
    expect(out[5]).toBeCloseTo(1.0, F32_DIGITS);
    // x and y components should be unchanged (exact zeros)
    expect(out[0]).toBe(0);
    expect(out[1]).toBe(0);
    expect(out[3]).toBe(0);
    expect(out[4]).toBe(0);
  });

  it('returns equilibrium − amplitude·displacement at phase = -1', () => {
    const out = computeDisplacedPositions(h2Atoms, h2Stretch, 0.3, -1);
    // Atom 0: z = -0.7 + (-1)·0.3·(-1) = -0.4
    expect(out[2]).toBeCloseTo(-0.4, F32_DIGITS);
    // Atom 1: z = +0.7 + (-1)·0.3·(+1) = +0.4
    expect(out[5]).toBeCloseTo(0.4, F32_DIGITS);
  });

  it('produces a symmetric H₂ stretch trajectory between phase = +1 and -1', () => {
    const peakPlus = computeDisplacedPositions(h2Atoms, h2Stretch, 0.3, +1);
    const peakMinus = computeDisplacedPositions(h2Atoms, h2Stretch, 0.3, -1);
    // The two atoms swap their roles when phase reverses sign:
    // - At phase = +1, atom 0 is at -1.0, atom 1 at +1.0 (stretched)
    // - At phase = -1, atom 0 is at -0.4, atom 1 at +0.4 (compressed)
    // The midpoint between atoms is invariant (mode is non-translational):
    const mp1 = (peakPlus[2] + peakPlus[5]) * 0.5;
    const mp2 = (peakMinus[2] + peakMinus[5]) * 0.5;
    expect(mp1).toBeCloseTo(0, F32_DIGITS);
    expect(mp2).toBeCloseTo(0, F32_DIGITS);
  });

  it('moves H₂O bend-mode atoms by their expected per-atom displacements', () => {
    const out = computeDisplacedPositions(h2oAtoms, h2oBend, 0.5, 1);
    // O: z = 0 + 1·0.5·(-0.05) = -0.025
    expect(out[0]).toBeCloseTo(0, F32_DIGITS);
    expect(out[1]).toBeCloseTo(0, F32_DIGITS);
    expect(out[2]).toBeCloseTo(-0.025, F32_DIGITS);
    // H1: x = 0.78 + 1·0.5·(-0.4) = 0.58, z = 0.59 + 1·0.5·0.3 = 0.74
    expect(out[3]).toBeCloseTo(0.58, F32_DIGITS);
    expect(out[4]).toBeCloseTo(0, F32_DIGITS);
    expect(out[5]).toBeCloseTo(0.74, F32_DIGITS);
    // H2: x = -0.78 + 1·0.5·0.4 = -0.58, z = 0.59 + 1·0.5·0.3 = 0.74
    expect(out[6]).toBeCloseTo(-0.58, F32_DIGITS);
    expect(out[7]).toBeCloseTo(0, F32_DIGITS);
    expect(out[8]).toBeCloseTo(0.74, F32_DIGITS);
  });

  it('returns a Float32Array of length 3·atoms.length', () => {
    const out = computeDisplacedPositions(h2oAtoms, h2oBend, 0.5, 0.5);
    expect(out).toBeInstanceOf(Float32Array);
    expect(out.length).toBe(9);
  });

  it('treats missing displacement entries as zero (no exception)', () => {
    // Only one displacement vector for two atoms — atom 1 should stay put.
    const out = computeDisplacedPositions(
      h2Atoms,
      [[0, 0, -1]] as AnimDisplacement[],
      0.3,
      1
    );
    // Atom 0: z = -0.7 + 1·0.3·(-1) = -1.0
    expect(out[2]).toBeCloseTo(-1.0, F32_DIGITS);
    // Atom 1: should stay at +0.7
    expect(out[5]).toBeCloseTo(0.7, F32_DIGITS);
  });

  it('does not mutate the input atoms or displacement arrays', () => {
    const atomsSnap = JSON.parse(JSON.stringify(h2Atoms));
    const dispSnap = JSON.parse(JSON.stringify(h2Stretch));
    computeDisplacedPositions(h2Atoms, h2Stretch, 0.42, 0.7);
    expect(h2Atoms).toEqual(atomsSnap);
    expect(h2Stretch).toEqual(dispSnap);
  });

  it('handles an empty atoms list', () => {
    const out = computeDisplacedPositions([], [], 0.5, 1);
    expect(out).toBeInstanceOf(Float32Array);
    expect(out.length).toBe(0);
  });

  it('writes into a caller-supplied scratch buffer (zero allocation path)', () => {
    const scratch = new Float32Array(6);
    const returned = computeDisplacedPositions(
      h2Atoms,
      h2Stretch,
      0.3,
      1,
      scratch
    );
    // Same buffer object — no copy
    expect(returned).toBe(scratch);
    expect(scratch[2]).toBeCloseTo(-1.0, F32_DIGITS);
    expect(scratch[5]).toBeCloseTo(1.0, F32_DIGITS);
  });

  it('falls back to a fresh allocation if the scratch buffer is too small', () => {
    // Caller passes a buffer too small to hold all atoms — must NOT write
    // out of bounds. The contract is to allocate a new buffer in that case.
    const scratch = new Float32Array(3); // only one atom's worth
    const returned = computeDisplacedPositions(
      h2Atoms,
      h2Stretch,
      0.3,
      1,
      scratch
    );
    expect(returned).not.toBe(scratch);
    expect(returned.length).toBe(6);
    // The original (too-small) buffer must be left untouched.
    expect(scratch[0]).toBe(0);
    expect(scratch[1]).toBe(0);
    expect(scratch[2]).toBe(0);
  });

  it('scales linearly with amplitude (twice amplitude → twice displacement)', () => {
    const out1 = computeDisplacedPositions(h2Atoms, h2Stretch, 0.3, 1);
    const out2 = computeDisplacedPositions(h2Atoms, h2Stretch, 0.6, 1);
    // Δz of atom 0: 0.3 vs 0.6 below equilibrium
    const dz1 = out1[2] - h2Atoms[0][3];
    const dz2 = out2[2] - h2Atoms[0][3];
    expect(dz2).toBeCloseTo(2 * dz1, F32_DIGITS);
  });
});

// ============================================================================
// computeBondGeometry
// ============================================================================

describe('computeBondGeometry', () => {
  it('computes midpoint, length, and direction for a unit-z bond', () => {
    // Atom 0 at origin, atom 1 at (0, 0, 1).
    const positions = new Float32Array([0, 0, 0, 0, 0, 1]);
    const bonds: Array<readonly [number, number]> = [[0, 1]];
    const out = computeBondGeometry(positions, bonds);
    expect(out).toHaveLength(1);
    expect(out[0].midpoint[0]).toBeCloseTo(0, F32_DIGITS);
    expect(out[0].midpoint[1]).toBeCloseTo(0, F32_DIGITS);
    expect(out[0].midpoint[2]).toBeCloseTo(0.5, F32_DIGITS);
    expect(out[0].length).toBeCloseTo(1.0, F32_DIGITS);
    expect(out[0].direction[0]).toBeCloseTo(0, F32_DIGITS);
    expect(out[0].direction[1]).toBeCloseTo(0, F32_DIGITS);
    expect(out[0].direction[2]).toBeCloseTo(1, F32_DIGITS);
  });

  it('computes midpoint, length, and direction for a unit-y bond', () => {
    // Atom 0 at origin, atom 1 at (0, 1, 0).
    const positions = new Float32Array([0, 0, 0, 0, 1, 0]);
    const bonds: Array<readonly [number, number]> = [[0, 1]];
    const out = computeBondGeometry(positions, bonds);
    expect(out[0].direction[1]).toBeCloseTo(1, F32_DIGITS);
    expect(out[0].length).toBeCloseTo(1, F32_DIGITS);
  });

  it('preserves bond length under rigid translation of both atoms', () => {
    // Translate both atoms by (5, -3, 7); length must be invariant.
    const positions = new Float32Array([5, -3, 7, 5, -3, 8]);
    const bonds: Array<readonly [number, number]> = [[0, 1]];
    const out = computeBondGeometry(positions, bonds);
    expect(out[0].length).toBeCloseTo(1.0, F32_DIGITS);
    expect(out[0].midpoint[0]).toBeCloseTo(5, F32_DIGITS);
    expect(out[0].midpoint[1]).toBeCloseTo(-3, F32_DIGITS);
    expect(out[0].midpoint[2]).toBeCloseTo(7.5, F32_DIGITS);
  });

  it('stretches bond length when atoms separate along the bond axis', () => {
    // Original: 1 bohr; stretched: 2 bohr.
    const stretched = new Float32Array([0, 0, -1, 0, 0, 1]);
    const out = computeBondGeometry(stretched, [[0, 1]]);
    expect(out[0].length).toBeCloseTo(2.0, F32_DIGITS);
    // Direction is still +z
    expect(out[0].direction[2]).toBeCloseTo(1, F32_DIGITS);
    // Midpoint is at origin
    expect(out[0].midpoint[2]).toBeCloseTo(0, F32_DIGITS);
  });

  it('contracts bond length when atoms approach each other', () => {
    const compressed = new Float32Array([0, 0, -0.25, 0, 0, 0.25]);
    const out = computeBondGeometry(compressed, [[0, 1]]);
    expect(out[0].length).toBeCloseTo(0.5, F32_DIGITS);
  });

  it('handles a degenerate bond (zero length) without producing NaN', () => {
    const positions = new Float32Array([1, 2, 3, 1, 2, 3]);
    const out = computeBondGeometry(positions, [[0, 1]]);
    expect(out[0].length).toBeLessThan(1e-9);
    // Direction defaults to the cylinder default axis (+Y)
    expect(Number.isNaN(out[0].direction[0])).toBe(false);
    expect(Number.isNaN(out[0].direction[1])).toBe(false);
    expect(Number.isNaN(out[0].direction[2])).toBe(false);
    expect(out[0].direction).toEqual([0, 1, 0]);
  });

  it('returns one BondGeometry per input bond', () => {
    // 3 atoms in a line: A — B — C
    const positions = new Float32Array([0, 0, 0, 0, 0, 1, 0, 0, 2]);
    const bonds: Array<readonly [number, number]> = [
      [0, 1],
      [1, 2],
    ];
    const out = computeBondGeometry(positions, bonds);
    expect(out).toHaveLength(2);
    expect(out[0].midpoint[2]).toBeCloseTo(0.5, F32_DIGITS);
    expect(out[1].midpoint[2]).toBeCloseTo(1.5, F32_DIGITS);
    expect(out[0].length).toBeCloseTo(1, F32_DIGITS);
    expect(out[1].length).toBeCloseTo(1, F32_DIGITS);
  });

  it('computes a unit direction vector (length = 1) for any non-degenerate bond', () => {
    // An arbitrary diagonal bond.
    const positions = new Float32Array([0, 0, 0, 1, 2, 2]);
    const out = computeBondGeometry(positions, [[0, 1]]);
    const [dx, dy, dz] = out[0].direction;
    const norm = Math.sqrt(dx * dx + dy * dy + dz * dz);
    expect(norm).toBeCloseTo(1.0, F32_DIGITS);
    // And the length is the original Euclidean distance √(1+4+4) = 3
    expect(out[0].length).toBeCloseTo(3.0, F32_DIGITS);
  });

  it('returns an empty array when no bonds are passed', () => {
    const positions = new Float32Array([0, 0, 0, 0, 0, 1]);
    const out = computeBondGeometry(positions, []);
    expect(out).toEqual([]);
  });
});

// ============================================================================
// Integration: end-to-end simulation of a single frame
// ============================================================================

describe('animatedMoleculePositions: full-frame integration', () => {
  it('replicates the H₂ stretch frame trajectory bond-by-bond', () => {
    const atoms: AnimAtom[] = [
      [1, 0, 0, -0.7],
      [1, 0, 0, +0.7],
    ];
    const stretch: AnimDisplacement[] = [
      [0, 0, -1],
      [0, 0, +1],
    ];
    // At phase = 0, the bond should have its equilibrium length 1.4.
    const equil = computeDisplacedPositions(atoms, stretch, 0.3, 0);
    const equilBond = computeBondGeometry(equil, [[0, 1]]);
    expect(equilBond[0].length).toBeCloseTo(1.4, F32_DIGITS);

    // At phase = +1, atoms have moved away by 0.3 each → bond length = 2.0
    const stretched = computeDisplacedPositions(atoms, stretch, 0.3, 1);
    const stretchedBond = computeBondGeometry(stretched, [[0, 1]]);
    expect(stretchedBond[0].length).toBeCloseTo(2.0, F32_DIGITS);

    // At phase = -1, atoms have moved toward each other by 0.3 each → bond length = 0.8
    const compressed = computeDisplacedPositions(atoms, stretch, 0.3, -1);
    const compressedBond = computeBondGeometry(compressed, [[0, 1]]);
    expect(compressedBond[0].length).toBeCloseTo(0.8, F32_DIGITS);

    // Bond direction is still +z throughout (sign of dz never flips at these
    // amplitudes — atoms do not pass through each other)
    expect(equilBond[0].direction[2]).toBeCloseTo(1, F32_DIGITS);
    expect(stretchedBond[0].direction[2]).toBeCloseTo(1, F32_DIGITS);
    expect(compressedBond[0].direction[2]).toBeCloseTo(1, F32_DIGITS);
  });

  it('preserves the H₂ midpoint at the origin throughout the trajectory', () => {
    const atoms: AnimAtom[] = [
      [1, 0, 0, -0.7],
      [1, 0, 0, +0.7],
    ];
    const stretch: AnimDisplacement[] = [
      [0, 0, -1],
      [0, 0, +1],
    ];
    for (const phase of [-1, -0.5, 0, 0.5, 1]) {
      const pos = computeDisplacedPositions(atoms, stretch, 0.3, phase);
      const bond = computeBondGeometry(pos, [[0, 1]]);
      expect(bond[0].midpoint[2]).toBeCloseTo(0, F32_DIGITS);
    }
  });
});

// Reference: F64_DIGITS is reserved for any future tests that operate purely on
// non-Float32 values. Currently all assertions use F32_DIGITS because the data
// flows through a Float32Array at some point.
void F64_DIGITS;
