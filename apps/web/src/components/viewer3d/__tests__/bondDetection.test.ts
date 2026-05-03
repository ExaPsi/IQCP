/**
 * Unit tests for bond detection algorithm.
 *
 * Tests use actual coordinates from the preset system JSON files
 * to validate correct bond counts for all 5 preset molecules.
 *
 * @module viewer3d/__tests__/bondDetection
 */

import { describe, it, expect } from 'vitest';
import { detectBonds } from '../bondDetection';
import type { Atom3D } from '../AtomSpheres';

// =============================================================================
// Group 1: Edge Cases
// =============================================================================

describe('detectBonds - edge cases', () => {
  it('returns empty array for empty input', () => {
    expect(detectBonds([])).toHaveLength(0);
  });

  it('returns empty array for single atom', () => {
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
    ];
    expect(detectBonds(atoms)).toHaveLength(0);
  });

  it('returns empty array for two atoms very far apart', () => {
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
      { symbol: 'H', position: [0, 0, 100] },
    ];
    expect(detectBonds(atoms)).toHaveLength(0);
  });

  it('returns empty array for two atoms at same position (degenerate)', () => {
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
      { symbol: 'H', position: [0, 0, 0] },
    ];
    expect(detectBonds(atoms)).toHaveLength(0);
  });

  it('handles unknown element symbols with default radius', () => {
    // DEFAULT_COVALENT_RADIUS_BOHR = 1.0
    // threshold = 1.2 * (1.0 + 1.0) = 2.4
    // distance = 1.5 < 2.4, so bonded
    const atoms: Atom3D[] = [
      { symbol: 'Zz', position: [0, 0, 0] },
      { symbol: 'Zz', position: [0, 0, 1.5] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(1);
  });
});

// =============================================================================
// Group 2: Threshold Boundary
// =============================================================================

describe('detectBonds - threshold boundary', () => {
  // H-H threshold = 1.2 * (0.586 + 0.586) = 1.4064 bohr

  it('detects bond just below threshold (H-H at 1.40 bohr)', () => {
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
      { symbol: 'H', position: [0, 0, 1.40] },
    ];
    expect(detectBonds(atoms)).toHaveLength(1);
  });

  it('detects no bond just above threshold (H-H at 1.41 bohr)', () => {
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
      { symbol: 'H', position: [0, 0, 1.41] },
    ];
    expect(detectBonds(atoms)).toHaveLength(0);
  });

  it('detects no bond exactly at threshold (strict less-than)', () => {
    // threshold = 1.2 * (0.586 + 0.586) = 1.4064
    const threshold = 1.2 * (0.586 + 0.586);
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
      { symbol: 'H', position: [0, 0, threshold] },
    ];
    expect(detectBonds(atoms)).toHaveLength(0);
  });
});

// =============================================================================
// Group 3: All 5 Presets (AC4, AC6)
// Coordinates taken directly from content/presets/systems/ JSON files
// =============================================================================

describe('detectBonds - preset molecules', () => {
  it('detects 1 bond in H2 (h2_sto3g_r1.4)', () => {
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0.0, 0.0, 0.0] },
      { symbol: 'H', position: [0.0, 0.0, 1.4] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(1);
  });

  it('detects 1 bond in HeH+ (heh_plus_sto3g)', () => {
    const atoms: Atom3D[] = [
      { symbol: 'He', position: [0.0, 0.0, 0.0] },
      { symbol: 'H', position: [0.0, 0.0, 1.4632] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(1);
  });

  it('detects 1 bond in LiH (lih_sto3g)', () => {
    const atoms: Atom3D[] = [
      { symbol: 'Li', position: [0.0, 0.0, 0.0] },
      { symbol: 'H', position: [0.0, 0.0, 3.0139] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(1);
  });

  it('detects 2 bonds in H2O (h2o_sto3g)', () => {
    const atoms: Atom3D[] = [
      { symbol: 'O', position: [0.0, 0.0, 0.2217282] },
      { symbol: 'H', position: [0.0, 1.4305447, -0.8869128] },
      { symbol: 'H', position: [0.0, -1.4305447, -0.8869128] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(2);
  });

  it('detects 3 bonds in NH3 (nh3_sto3g)', () => {
    const atoms: Atom3D[] = [
      { symbol: 'N', position: [0.0, 0.0, 0.219705] },
      { symbol: 'H', position: [0.0, 1.7714918, -0.512645] },
      { symbol: 'H', position: [1.5342036, -0.8857459, -0.512645] },
      { symbol: 'H', position: [-1.5342036, -0.8857459, -0.512645] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(3);
  });
});

// =============================================================================
// Group 4: Bond Properties Validation
// =============================================================================

describe('detectBonds - bond properties', () => {
  it('computes correct properties for H2 bond', () => {
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
      { symbol: 'H', position: [0, 0, 1.4] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(1);

    const bond = bonds[0];
    expect(bond.atomIndexA).toBe(0);
    expect(bond.atomIndexB).toBe(1);
    expect(bond.length).toBeCloseTo(1.4, 10);
    expect(bond.midpoint[0]).toBeCloseTo(0, 10);
    expect(bond.midpoint[1]).toBeCloseTo(0, 10);
    expect(bond.midpoint[2]).toBeCloseTo(0.7, 10);
    expect(bond.direction[0]).toBeCloseTo(0, 10);
    expect(bond.direction[1]).toBeCloseTo(0, 10);
    expect(bond.direction[2]).toBeCloseTo(1.0, 10);
  });

  it('computes correct bond lengths for H2O', () => {
    const atoms: Atom3D[] = [
      { symbol: 'O', position: [0.0, 0.0, 0.2217282] },
      { symbol: 'H', position: [0.0, 1.4305447, -0.8869128] },
      { symbol: 'H', position: [0.0, -1.4305447, -0.8869128] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(2);

    // Both O-H bonds should have approximately the same length (~1.8098 bohr)
    expect(bonds[0].length).toBeCloseTo(1.8098, 2);
    expect(bonds[1].length).toBeCloseTo(1.8098, 2);
  });

  it('NH3 bonds all connect N (index 0) to H atoms', () => {
    const atoms: Atom3D[] = [
      { symbol: 'N', position: [0.0, 0.0, 0.219705] },
      { symbol: 'H', position: [0.0, 1.7714918, -0.512645] },
      { symbol: 'H', position: [1.5342036, -0.8857459, -0.512645] },
      { symbol: 'H', position: [-1.5342036, -0.8857459, -0.512645] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(3);

    // All bonds should have atomIndexA = 0 (nitrogen)
    for (const bond of bonds) {
      expect(bond.atomIndexA).toBe(0);
    }
    // Atom indices B should be 1, 2, 3
    expect(bonds[0].atomIndexB).toBe(1);
    expect(bonds[1].atomIndexB).toBe(2);
    expect(bonds[2].atomIndexB).toBe(3);
  });
});

// =============================================================================
// Group 5: No False Positives
// =============================================================================

describe('detectBonds - no false positives', () => {
  it('H2O has no H-H bond', () => {
    const atoms: Atom3D[] = [
      { symbol: 'O', position: [0.0, 0.0, 0.2217282] },
      { symbol: 'H', position: [0.0, 1.4305447, -0.8869128] },
      { symbol: 'H', position: [0.0, -1.4305447, -0.8869128] },
    ];
    const bonds = detectBonds(atoms);

    // No bond should have both indices among 1 and 2 (H atoms)
    const hhBonds = bonds.filter(
      (b) => b.atomIndexA >= 1 && b.atomIndexB >= 1
    );
    expect(hhBonds).toHaveLength(0);
  });

  it('NH3 has no H-H bonds', () => {
    const atoms: Atom3D[] = [
      { symbol: 'N', position: [0.0, 0.0, 0.219705] },
      { symbol: 'H', position: [0.0, 1.7714918, -0.512645] },
      { symbol: 'H', position: [1.5342036, -0.8857459, -0.512645] },
      { symbol: 'H', position: [-1.5342036, -0.8857459, -0.512645] },
    ];
    const bonds = detectBonds(atoms);

    // No bond should have both indices among 1, 2, 3 (H atoms)
    const hhBonds = bonds.filter(
      (b) => b.atomIndexA >= 1 && b.atomIndexB >= 1
    );
    expect(hhBonds).toHaveLength(0);
  });
});

// =============================================================================
// Group 6: Custom Geometry (AC7)
// =============================================================================

describe('detectBonds - custom geometry', () => {
  it('detects bond between two C atoms at typical bond length', () => {
    // C-C threshold = 1.2 * (1.455 + 1.455) = 3.492 bohr
    // 2.83 < 3.492 --> bonded
    const atoms: Atom3D[] = [
      { symbol: 'C', position: [0, 0, 0] },
      { symbol: 'C', position: [0, 0, 2.83] },
    ];
    expect(detectBonds(atoms)).toHaveLength(1);
  });

  it('detects O-H bond at typical water distance', () => {
    // O-H threshold = 1.2 * (1.304 + 0.586) = 2.268 bohr
    // 1.81 < 2.268 --> bonded
    const atoms: Atom3D[] = [
      { symbol: 'O', position: [0, 0, 0] },
      { symbol: 'H', position: [0, 0, 1.81] },
    ];
    expect(detectBonds(atoms)).toHaveLength(1);
  });

  it('detects bond with mixed known and unknown elements', () => {
    // H threshold: 0.586, "X" uses DEFAULT_COVALENT_RADIUS_BOHR = 1.0
    // threshold = 1.2 * (0.586 + 1.0) = 1.9032 bohr
    // 1.5 < 1.9032 --> bonded
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
      { symbol: 'X', position: [0, 0, 1.5] },
    ];
    expect(detectBonds(atoms)).toHaveLength(1);
  });
});

// =============================================================================
// Group 7: Bond Disappearance (PES scanning scenario, T8)
// =============================================================================

describe('detectBonds - bond disappearance at large distances', () => {
  // H-H threshold = 1.2 * (0.586 + 0.586) = 1.4064 bohr
  const testDistances: [number, boolean][] = [
    [1.0, true],
    [1.3, true],
    [1.4, true],   // 1.4 < 1.4064 --> bonded
    [1.41, false],  // 1.41 > 1.4064 --> not bonded
    [1.5, false],
    [5.0, false],
    [100.0, false],
  ];

  for (const [distance, expectedBonded] of testDistances) {
    it(`H-H at ${distance} bohr: ${expectedBonded ? 'bonded' : 'not bonded'}`, () => {
      const atoms: Atom3D[] = [
        { symbol: 'H', position: [0, 0, 0] },
        { symbol: 'H', position: [0, 0, distance] },
      ];
      const bonds = detectBonds(atoms);
      expect(bonds).toHaveLength(expectedBonded ? 1 : 0);
    });
  }
});

// =============================================================================
// Group 8: Direction vector validation
// =============================================================================

describe('detectBonds - direction vectors', () => {
  it('computes unit direction vector along x-axis', () => {
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
      { symbol: 'H', position: [1.3, 0, 0] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(1);
    expect(bonds[0].direction[0]).toBeCloseTo(1.0, 10);
    expect(bonds[0].direction[1]).toBeCloseTo(0.0, 10);
    expect(bonds[0].direction[2]).toBeCloseTo(0.0, 10);
  });

  it('computes unit direction vector along y-axis', () => {
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
      { symbol: 'H', position: [0, 1.3, 0] },
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(1);
    expect(bonds[0].direction[0]).toBeCloseTo(0.0, 10);
    expect(bonds[0].direction[1]).toBeCloseTo(1.0, 10);
    expect(bonds[0].direction[2]).toBeCloseTo(0.0, 10);
  });

  it('computes unit direction vector for diagonal bond', () => {
    const atoms: Atom3D[] = [
      { symbol: 'H', position: [0, 0, 0] },
      { symbol: 'H', position: [0.7, 0.7, 0] }, // distance = sqrt(0.98) ~ 0.99
    ];
    const bonds = detectBonds(atoms);
    expect(bonds).toHaveLength(1);

    const len = Math.sqrt(0.7 * 0.7 + 0.7 * 0.7);
    expect(bonds[0].direction[0]).toBeCloseTo(0.7 / len, 10);
    expect(bonds[0].direction[1]).toBeCloseTo(0.7 / len, 10);
    expect(bonds[0].direction[2]).toBeCloseTo(0.0, 10);

    // Direction should be a unit vector
    const mag = Math.sqrt(
      bonds[0].direction[0] ** 2 +
        bonds[0].direction[1] ** 2 +
        bonds[0].direction[2] ** 2
    );
    expect(mag).toBeCloseTo(1.0, 10);
  });
});
