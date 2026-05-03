/**
 * Pure helpers for `FrequencyPanel.tsx`.
 *
 * Extracted into its own module so the component file exports only JSX
 * components (satisfies the `react-refresh/only-export-components` rule)
 * and so tests can import the helpers without pulling in React.
 *
 * @module components/scf/frequencyPanelLogic
 * @see US-102 Frequency Tab UI, AC1
 */

import type { FrequencyResult } from '../../worker/protocol';

/**
 * Row data for a single vibrational mode. Assembled from a FrequencyResult
 * by `buildFrequencyRows`.
 */
export interface FrequencyRow {
  modeIndex: number;         // 0-based — row order in result.frequenciesCm1
  frequency: number;         // cm⁻¹ (negative for imaginary)
  irIntensity: number;       // km/mol
  ramanActivity: number;     // Å⁴/amu
  depolarization: number;    // dimensionless
  reducedMass: number;       // amu
  forceConstant: number;     // mDyne/Å
  isImaginary: boolean;      // frequency < 0
}

/** Sort column key — union of FrequencyRow numeric fields. */
export type FrequencyColumnKey =
  | 'modeIndex'
  | 'frequency'
  | 'irIntensity'
  | 'ramanActivity'
  | 'depolarization'
  | 'reducedMass'
  | 'forceConstant';

/** Ascending or descending sort direction. */
export type SortDirection = 'asc' | 'desc';

/**
 * Assemble FrequencyRow[] from a FrequencyResult. Preserves input order
 * (mode 0 = result.frequenciesCm1[0]).
 */
export function buildFrequencyRows(
  result: FrequencyResult | null
): FrequencyRow[] {
  if (!result || result.frequenciesCm1.length === 0) return [];
  const n = result.frequenciesCm1.length;
  const rows: FrequencyRow[] = [];
  for (let i = 0; i < n; i++) {
    const freq = result.frequenciesCm1[i];
    rows.push({
      modeIndex: i,
      frequency: freq,
      irIntensity: result.irIntensitiesKmPerMol[i] ?? 0,
      ramanActivity: result.ramanActivitiesA4Amu[i] ?? 0,
      depolarization: result.depolarizationRatios[i] ?? 0,
      reducedMass: result.reducedMassesAmu[i] ?? 0,
      forceConstant: result.forceConstantsMdyne[i] ?? 0,
      isImaginary: freq < 0,
    });
  }
  return rows;
}

/**
 * Return a new sorted copy of `rows` by the given column and direction.
 * Stable: rows with equal keys preserve original order.
 */
export function sortFrequencyRows(
  rows: readonly FrequencyRow[],
  key: FrequencyColumnKey,
  direction: SortDirection
): FrequencyRow[] {
  const indexed = rows.map((r, i) => ({ row: r, origIdx: i }));
  indexed.sort((a, b) => {
    const va = a.row[key];
    const vb = b.row[key];
    if (va < vb) return direction === 'asc' ? -1 : 1;
    if (va > vb) return direction === 'asc' ? 1 : -1;
    return a.origIdx - b.origIdx;
  });
  return indexed.map((x) => x.row);
}

/**
 * Count imaginary frequencies in a FrequencyResult. Returns 0 for null.
 */
export function countImaginary(result: FrequencyResult | null): number {
  if (!result) return 0;
  let count = 0;
  for (const f of result.frequenciesCm1) {
    if (f < 0) count += 1;
  }
  return count;
}

/**
 * Find the index of the first imaginary mode (< 0), or -1 if none.
 */
export function firstImaginaryIndex(result: FrequencyResult | null): number {
  if (!result) return -1;
  for (let i = 0; i < result.frequenciesCm1.length; i++) {
    if (result.frequenciesCm1[i] < 0) return i;
  }
  return -1;
}

/**
 * Format a frequency value for display. Returns the absolute value; the
 * sign channel is carried by a separate "imaginary" label alongside the
 * value to satisfy WCAG 1.4.1 (color is not the only indicator).
 */
export function formatFrequency(value: number): string {
  const abs = Math.abs(value);
  return abs.toFixed(1);
}
