/**
 * Bond detection utility for 3D molecular visualization.
 *
 * Pure function with no React or Three.js dependencies. Detects covalent
 * bonds between atoms using the standard distance-based criterion:
 *
 *   bonded if dist(A, B) < BOND_TOLERANCE * (r_cov(A) + r_cov(B))
 *
 * where BOND_TOLERANCE = 1.2 accounts for bonds slightly longer than
 * ideal covalent distances while avoiding false positives.
 *
 * @module viewer3d/bondDetection
 */

import type { Atom3D } from './AtomSpheres';
import {
  COVALENT_RADII_BOHR,
  DEFAULT_COVALENT_RADIUS_BOHR,
  BOND_TOLERANCE,
} from './constants';

/**
 * A detected bond between two atoms.
 *
 * Contains the atom indices into the source Atom3D[] array
 * and precomputed geometric data needed for cylinder rendering.
 */
export interface Bond {
  /** Index of the first atom in the source Atom3D[] array */
  atomIndexA: number;
  /** Index of the second atom in the source Atom3D[] array */
  atomIndexB: number;
  /** Midpoint [x, y, z] between the two atom centers (bohr) */
  midpoint: [number, number, number];
  /** Unit direction vector from atom A to atom B */
  direction: [number, number, number];
  /** Euclidean distance between atom centers (bohr) */
  length: number;
}

/**
 * Detect bonds between atoms using covalent radius criterion.
 *
 * Two atoms A, B are bonded if:
 *   dist(A, B) < BOND_TOLERANCE * (r_cov(A) + r_cov(B))
 *
 * The algorithm is O(n^2/2) pairwise comparison, which is more than
 * sufficient for our molecule sizes (2-10 atoms). Each unique pair
 * (i, j) with i < j is checked exactly once, ensuring deterministic
 * output order and no self-bonds.
 *
 * @param atoms - Array of atoms with positions in bohr
 * @returns Array of detected bonds, ordered by atom index pairs
 */
export function detectBonds(atoms: Atom3D[]): Bond[] {
  const bonds: Bond[] = [];
  const n = atoms.length;

  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const a = atoms[i];
      const b = atoms[j];

      // Compute Euclidean distance in bohr
      const dx = b.position[0] - a.position[0];
      const dy = b.position[1] - a.position[1];
      const dz = b.position[2] - a.position[2];
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);

      // Skip zero-distance (degenerate geometry) to prevent division by zero
      if (dist < 1e-10) continue;

      // Look up covalent radii with fallback for unknown elements
      const rA = COVALENT_RADII_BOHR[a.symbol] ?? DEFAULT_COVALENT_RADIUS_BOHR;
      const rB = COVALENT_RADII_BOHR[b.symbol] ?? DEFAULT_COVALENT_RADIUS_BOHR;

      // Bond criterion: strict less-than (not <=)
      const threshold = BOND_TOLERANCE * (rA + rB);

      if (dist < threshold) {
        const invDist = 1.0 / dist;
        bonds.push({
          atomIndexA: i,
          atomIndexB: j,
          midpoint: [
            (a.position[0] + b.position[0]) / 2,
            (a.position[1] + b.position[1]) / 2,
            (a.position[2] + b.position[2]) / 2,
          ],
          direction: [dx * invDist, dy * invDist, dz * invDist],
          length: dist,
        });
      }
    }
  }

  return bonds;
}
