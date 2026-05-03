/**
 * Pure math helpers for AnimatedMolecule.tsx (US-102 UX refinement).
 *
 * Isolated in a no-R3F module so the geometric logic can be unit-tested in
 * Node without mounting a Three.js Canvas (which requires WebGL).
 *
 * Two helpers:
 *   - `computeDisplacedPositions` — given equilibrium atom positions, a
 *     displacement field, an amplitude, and a phase factor, returns the
 *     atom positions at that animation phase.
 *   - `computeBondGeometry` — given a flat array of atom positions and an
 *     atom-pair list, returns the midpoint, length, and unit direction of
 *     each bond. Used by AnimatedMolecule to update bond cylinder transforms
 *     each frame as the atoms move.
 *
 * Both helpers operate in the same length unit as the atom positions
 * (the existing viewer uses bohr; the displacement vectors are also in
 * bohr — see `viewer3d/DisplacementArrows.tsx` and `crates/qc-core` for
 * provenance).
 *
 * @module components/viewer3d/animatedMoleculePositions
 * @see US-102 Frequency Tab UI
 */

/**
 * One atom in a normal-mode animation, in the same flat shape used by the
 * Frequency tab and worker protocol: `[Z, x, y, z]`. Z is the atomic number
 * (unused for geometric computation but kept for API parity with the
 * existing worker payload). x, y, z are equilibrium coordinates.
 */
export type AnimAtom = readonly [number, number, number, number];

/**
 * One displacement vector for an atom, `[dx, dy, dz]`, in the same length
 * unit as the equilibrium positions.
 */
export type AnimDisplacement = readonly [number, number, number];

/**
 * Geometry of a single bond between two atoms.
 *
 * `length` is the Euclidean distance between the atom centers; `direction`
 * is the unit vector from atom A to atom B (or `[0, 1, 0]` — the cylinder
 * default axis — if the two atoms are degenerate). `midpoint` is the
 * centroid.
 */
export interface BondGeometry {
  midpoint: [number, number, number];
  length: number;
  direction: [number, number, number];
}

/**
 * Compute displaced atom positions for a normal mode at a given animation
 * phase.
 *
 * For each atom A:
 *
 *   r_A(t) = R_A + amplitude · phase · q^(k)_A
 *
 * where R_A is the equilibrium position and q^(k)_A is the displacement
 * vector for atom A in mode k (already normalized by the caller; this
 * function performs no normalization).
 *
 * The output is a flat `Float32Array` of length `3 · atoms.length` for
 * direct use with Three.js InstancedMesh / BufferAttribute APIs.
 *
 * Special cases:
 *   - `phase === 0` → returns equilibrium positions (round-trip through
 *     `Float32Array` introduces ~1e-7 error, but the contract is "no
 *     displacement applied")
 *   - `amplitude === 0` → returns equilibrium positions for any phase
 *   - displacement length shorter than atoms → missing entries treated as
 *     zero (atom stays at equilibrium)
 *
 * **Allocation control**: Pass `out` to reuse a pre-allocated scratch
 * buffer (the per-frame `useFrame` callback in `AnimatedMolecule` does
 * this). When `out` is omitted, the function allocates a new array — used
 * by tests and by code paths where allocation is not on the hot path.
 *
 * Pure function: does not mutate inputs.
 *
 * @param atoms        - Equilibrium atoms `[Z, x, y, z][]`
 * @param displacement - Per-atom displacement vectors `[dx, dy, dz][]`
 * @param amplitude    - Peak amplitude (slider value, in atom-position units)
 * @param phase        - Animation phase factor in `[-1, 1]` (typically
 *                       `Math.sin(2π · speed · t)`)
 * @param out          - Optional pre-allocated output buffer. Must have
 *                       length `≥ 3 · atoms.length`. If omitted, a new
 *                       array is allocated.
 * @returns A `Float32Array` (the same `out` buffer if provided) containing
 *          interleaved (x0, y0, z0, x1, y1, z1, …) displaced positions.
 */
export function computeDisplacedPositions(
  atoms: readonly AnimAtom[],
  displacement: readonly AnimDisplacement[],
  amplitude: number,
  phase: number,
  out?: Float32Array
): Float32Array {
  const target =
    out && out.length >= atoms.length * 3
      ? out
      : new Float32Array(atoms.length * 3);
  const scale = amplitude * phase;
  for (let i = 0; i < atoms.length; i++) {
    const atom = atoms[i];
    const x = atom[1];
    const y = atom[2];
    const z = atom[3];
    const d = displacement[i];
    if (d) {
      target[i * 3 + 0] = x + scale * d[0];
      target[i * 3 + 1] = y + scale * d[1];
      target[i * 3 + 2] = z + scale * d[2];
    } else {
      target[i * 3 + 0] = x;
      target[i * 3 + 1] = y;
      target[i * 3 + 2] = z;
    }
  }
  return target;
}

/**
 * Compute the bond geometry (midpoint, length, unit direction) for every
 * bond in `bonds`, given the current atom positions.
 *
 * `positions` is the flat `[x0, y0, z0, x1, y1, z1, …]` form returned by
 * `computeDisplacedPositions`. `bonds` is a list of `[atomIndexA,
 * atomIndexB]` pairs (typically the output of `detectBonds()` projected
 * down to indices, since detected bonds at equilibrium can be reused
 * across frames as long as no bond breaks during the animation — which is
 * physically impossible at the small amplitudes in scope here).
 *
 * For a degenerate bond (two atoms at the same position), `length` is set
 * to a tiny epsilon and `direction` defaults to the cylinder's natural
 * `+Y` axis so Three.js does not produce NaN quaternions.
 *
 * @param positions - Flat (x0, y0, z0, x1, y1, z1, …) array
 * @param bonds     - Atom-index pairs `[A, B]`
 * @returns A new array of `BondGeometry` records, one per input bond.
 */
export function computeBondGeometry(
  positions: Float32Array,
  bonds: ReadonlyArray<readonly [number, number]>
): BondGeometry[] {
  const result: BondGeometry[] = new Array(bonds.length);
  for (let i = 0; i < bonds.length; i++) {
    const [a, b] = bonds[i];
    const ax = positions[a * 3 + 0];
    const ay = positions[a * 3 + 1];
    const az = positions[a * 3 + 2];
    const bx = positions[b * 3 + 0];
    const by = positions[b * 3 + 1];
    const bz = positions[b * 3 + 2];

    const dx = bx - ax;
    const dy = by - ay;
    const dz = bz - az;
    const length = Math.sqrt(dx * dx + dy * dy + dz * dz);

    const midpoint: [number, number, number] = [
      (ax + bx) * 0.5,
      (ay + by) * 0.5,
      (az + bz) * 0.5,
    ];

    let direction: [number, number, number];
    if (length < 1e-10) {
      // Degenerate: pick the cylinder's default axis to avoid NaN quaternion.
      direction = [0, 1, 0];
    } else {
      const inv = 1.0 / length;
      direction = [dx * inv, dy * inv, dz * inv];
    }

    result[i] = { midpoint, length, direction };
  }
  return result;
}
