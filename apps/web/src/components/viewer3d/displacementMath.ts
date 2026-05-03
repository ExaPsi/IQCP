/**
 * Pure math helpers for DisplacementArrows.tsx.
 *
 * Isolated in a no-R3F module so that the logic can be unit-tested in Node
 * without mounting a Three.js Canvas (which requires WebGL).
 *
 * @module components/viewer3d/displacementMath
 * @see US-102 Frequency Tab UI, AC3
 */

/**
 * A single Cartesian atomic displacement (Δx, Δy, Δz) in whatever units the
 * caller supplies. `normalizeDisplacement` below rescales to a unit-norm
 * vector field which the viewer then multiplies by the user-chosen amplitude.
 */
export type DisplacementVec = [number, number, number];

/**
 * Compute the L2 norm of the full displacement field (concatenating all atom
 * vectors). Used for `normalizeDisplacement`.
 */
export function displacementFieldNorm(
  displacement: readonly DisplacementVec[]
): number {
  let sumSq = 0;
  for (const [dx, dy, dz] of displacement) {
    sumSq += dx * dx + dy * dy + dz * dz;
  }
  return Math.sqrt(sumSq);
}

/**
 * Rescale a cartesian normal mode vector so that the concatenated field has
 * unit L2 norm. Returns a new array; the input is not mutated. For a
 * zero-norm input (all atoms fixed), returns the input unchanged.
 *
 * The resulting unit-norm field is ideal for display: the viewer multiplies
 * by the user amplitude slider (in Å) so that every mode has a consistent
 * visual magnitude regardless of how qc-core normalized the raw
 * `normalModesCartesian` output.
 */
export function normalizeDisplacement(
  displacement: readonly DisplacementVec[]
): DisplacementVec[] {
  const norm = displacementFieldNorm(displacement);
  if (norm === 0) return displacement.map((v) => [v[0], v[1], v[2]]);
  return displacement.map((v) => [v[0] / norm, v[1] / norm, v[2] / norm]);
}

/**
 * Compute the animation phase factor `sin(2π · speed · t)` for a given
 * elapsed time and speed multiplier. Returned value ∈ [-1, 1].
 *
 * Exported so the animation loop's math can be unit-tested separately
 * from `useFrame`.
 */
export function animationPhase(elapsedSeconds: number, speed: number): number {
  return Math.sin(2 * Math.PI * speed * elapsedSeconds);
}

/**
 * Compute the endpoint of an atom's displacement arrow at a given animation
 * phase.
 *
 *   end = start + phase · amplitude · displacement
 *
 * @param start              - Atom position (x, y, z) in Angstroms
 * @param displacement       - Unit-norm displacement vector for the atom
 * @param amplitudeBohr  - User-chosen amplitude
 * @param phase              - `animationPhase(...)` ∈ [-1, 1]
 * @returns                    Endpoint [x, y, z] in Angstroms
 */
export function displacementEndpoint(
  start: readonly [number, number, number],
  displacement: readonly [number, number, number],
  amplitudeBohr: number,
  phase: number
): [number, number, number] {
  return [
    start[0] + phase * amplitudeBohr * displacement[0],
    start[1] + phase * amplitudeBohr * displacement[1],
    start[2] + phase * amplitudeBohr * displacement[2],
  ];
}

/**
 * Wrap a mode index modulo the total mode count, handling negative indices
 * (e.g., from a "previous" button). Returns 0 for an empty list.
 */
export function wrapModeIndex(index: number, total: number): number {
  if (total <= 0) return 0;
  return ((index % total) + total) % total;
}
