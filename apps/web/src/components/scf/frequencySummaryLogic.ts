/**
 * Pure helpers for `FrequencySummaryPanel.tsx`.
 *
 * Extracted into its own module so the component file exports only JSX
 * components (satisfies the `react-refresh/only-export-components` rule)
 * and so tests can import the helpers without pulling in React.
 *
 * @module components/scf/frequencySummaryLogic
 * @see US-102 Frequency Tab UI, AC7
 */

import type { RotorType } from '../../worker/protocol';

/**
 * Compute the magnitude of a 3-vector (e.g. dipole in e·bohr or Debye).
 */
export function dipoleMagnitude(
  vec: readonly [number, number, number]
): number {
  return Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

/**
 * Compute the isotropic polarizability (trace/3) from a symmetric 3×3 tensor.
 */
export function isotropicPolarizability(
  tensor: readonly [
    readonly [number, number, number],
    readonly [number, number, number],
    readonly [number, number, number],
  ]
): number {
  return (tensor[0][0] + tensor[1][1] + tensor[2][2]) / 3;
}

/**
 * Human-readable label for a RotorType variant.
 */
export function formatRotorType(t: RotorType): string {
  switch (t) {
    case 'atom':
      return 'Atom';
    case 'linear':
      return 'Linear';
    case 'spherical_top':
      return 'Spherical top';
    case 'symmetric_top':
      return 'Symmetric top';
    case 'asymmetric_top':
      return 'Asymmetric top';
  }
}

/**
 * Format a rotational constant for display.
 *
 *   Infinity  → "—"  (atom all three; linear first)
 *   Finite    → value.toFixed(3)
 *   NaN       → "—"
 */
export function formatRotConst(x: number): string {
  if (!Number.isFinite(x)) return '—';
  return x.toFixed(3);
}
