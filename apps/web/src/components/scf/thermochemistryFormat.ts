/**
 * Pure formatting helpers for `ThermochemistryPanel.tsx`.
 *
 * Extracted into its own module so the component file exports only JSX
 * components (satisfies the `react-refresh/only-export-components` rule)
 * and so tests can import the helpers without pulling in React.
 *
 * @module components/scf/thermochemistryFormat
 * @see US-102 Frequency Tab UI, AC4
 */

import {
  convertEnergy,
  convertEntropyPerK,
  type EnergyUnitsMode,
} from '../../lib/units';

/**
 * Format an energy value (Hartree) in the requested display units.
 *
 * Hartree uses 6 decimal places; kcal/mol and kJ/mol use 3 decimal places.
 */
export function formatEnergyValue(
  valueHa: number,
  mode: EnergyUnitsMode
): string {
  const converted = convertEnergy(valueHa, mode);
  if (!Number.isFinite(converted)) return '—';
  const decimals = mode === 'hartree' ? 6 : 3;
  return converted.toFixed(decimals);
}

/**
 * Format an entropy or heat-capacity value in Ha/(mol·K) in the requested
 * display units.
 *
 * Hartree uses scientific notation (4 sig figs) because the raw values are
 * O(1e-5); kcal/(mol·K) and kJ/(mol·K) use fixed notation (5 decimal places).
 */
export function formatEntropyValue(
  valueHaPerK: number,
  mode: EnergyUnitsMode
): string {
  const converted = convertEntropyPerK(valueHaPerK, mode);
  if (!Number.isFinite(converted)) return '—';
  if (mode === 'hartree') {
    return converted.toExponential(4);
  }
  return converted.toFixed(5);
}
