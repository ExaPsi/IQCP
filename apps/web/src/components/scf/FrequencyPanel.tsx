/**
 * FrequencyPanel — sortable 7-column table of vibrational modes.
 *
 * Part of the Frequency tab workflow in Module C (SCF Sandbox). Renders one
 * row per normal mode with: mode index, frequency (cm⁻¹), IR intensity
 * (km/mol), Raman activity (Å⁴/amu), depolarization ratio, reduced mass
 * (amu), and force constant (mDyne/Å).
 *
 * Sortable by clicking any column header. Imaginary frequencies are
 * highlighted in red and labeled "imaginary" (dual-signal to satisfy
 * WCAG 1.4.1 — color is not the only indicator).
 *
 * Clicking a row calls `onSelectMode(modeIndex)` so the spectrum plot and
 * normal-mode viewer can sync to the selected mode.
 *
 * Pure helpers (`buildFrequencyRows`, `sortFrequencyRows`, `countImaginary`,
 * `firstImaginaryIndex`, `formatFrequency`) are defined in
 * `./frequencyPanelLogic.ts` so this file exports only the React component.
 *
 * @module components/scf/FrequencyPanel
 * @see US-102 Frequency Tab UI, AC1
 */

import { useMemo, useState, useCallback } from 'react';
import type { FrequencyResult } from '../../worker/protocol';
import {
  buildFrequencyRows,
  sortFrequencyRows,
  countImaginary,
  formatFrequency,
  type FrequencyRow,
  type FrequencyColumnKey,
  type SortDirection,
} from './frequencyPanelLogic';

// ============================================================================
// Props
// ============================================================================

export interface FrequencyPanelProps {
  /** Frequency result from WASM, or null for empty state. */
  result: FrequencyResult | null;
  /** Currently-selected mode index (0-based), or null. */
  selectedMode: number | null;
  /** Callback invoked when the user clicks a row. */
  onSelectMode: (modeIndex: number) => void;
}

// ============================================================================
// Column definitions
// ============================================================================

interface FrequencyColumn {
  key: FrequencyColumnKey;
  label: string;
  ariaLabel: string;
  format: (row: FrequencyRow) => string;
  align: 'left' | 'right';
}

const COLUMNS: readonly FrequencyColumn[] = [
  {
    key: 'modeIndex',
    label: '#',
    ariaLabel: 'Mode number',
    format: (r) => `${r.modeIndex + 1}`,
    align: 'left',
  },
  {
    key: 'frequency',
    label: 'ν̃ (cm⁻¹)',
    ariaLabel: 'Vibrational frequency in wavenumbers',
    format: (r) =>
      r.isImaginary
        ? `−${formatFrequency(r.frequency)}i`
        : formatFrequency(r.frequency),
    align: 'right',
  },
  {
    key: 'irIntensity',
    label: 'I_IR (km/mol)',
    ariaLabel: 'IR absorption intensity in km per mole',
    format: (r) => r.irIntensity.toFixed(3),
    align: 'right',
  },
  {
    key: 'ramanActivity',
    label: 'I_Raman (Å⁴/amu)',
    ariaLabel: 'Raman scattering activity in Å⁴ per amu',
    format: (r) => r.ramanActivity.toFixed(3),
    align: 'right',
  },
  {
    key: 'depolarization',
    label: 'ρ',
    ariaLabel: 'Depolarization ratio',
    format: (r) => r.depolarization.toFixed(3),
    align: 'right',
  },
  {
    key: 'reducedMass',
    label: 'μ (amu)',
    ariaLabel: 'Reduced mass in amu',
    format: (r) => r.reducedMass.toFixed(4),
    align: 'right',
  },
  {
    key: 'forceConstant',
    label: 'k (mDyne/Å)',
    ariaLabel: 'Force constant in mDyne per Angstrom',
    format: (r) => r.forceConstant.toFixed(4),
    align: 'right',
  },
];

// ============================================================================
// Component
// ============================================================================

/**
 * Sortable seven-column table of vibrational modes with click-to-select.
 *
 * @see US-102 Frequency Tab UI, AC1
 */
export function FrequencyPanel({
  result,
  selectedMode,
  onSelectMode,
}: FrequencyPanelProps): JSX.Element {
  const [sortKey, setSortKey] = useState<FrequencyColumnKey>('modeIndex');
  const [sortDirection, setSortDirection] = useState<SortDirection>('asc');

  const rawRows = useMemo(() => buildFrequencyRows(result), [result]);
  const sortedRows = useMemo(
    () => sortFrequencyRows(rawRows, sortKey, sortDirection),
    [rawRows, sortKey, sortDirection]
  );

  const handleSort = useCallback(
    (key: FrequencyColumnKey) => {
      if (key === sortKey) {
        setSortDirection((prev) => (prev === 'asc' ? 'desc' : 'asc'));
      } else {
        setSortKey(key);
        setSortDirection('asc');
      }
    },
    [sortKey]
  );

  // Empty-state card
  if (!result) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 text-center">
        <div className="text-sm font-semibold text-slate-700 mb-1">
          Vibrational Modes
        </div>
        <div className="text-xs text-slate-500">
          Run frequency analysis to see vibrational modes.
        </div>
      </div>
    );
  }

  if (rawRows.length === 0) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 text-center">
        <div className="text-sm font-semibold text-slate-700 mb-1">
          Vibrational Modes
        </div>
        <div className="text-xs text-slate-500">
          No vibrational modes (atomic system).
        </div>
      </div>
    );
  }

  const imaginaryCount = countImaginary(result);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      <div className="px-4 py-3 bg-slate-50 border-b border-slate-200 flex items-center justify-between">
        <h3 className="text-sm font-semibold text-slate-700">
          Vibrational Modes
        </h3>
        <span className="text-xs text-slate-500">
          {rawRows.length} mode{rawRows.length !== 1 ? 's' : ''}
          {imaginaryCount > 0 && (
            <span className="ml-2 text-red-600 font-semibold">
              ({imaginaryCount} imaginary)
            </span>
          )}
        </span>
      </div>

      <div className="overflow-x-auto">
        <table
          className="w-full text-xs"
          role="grid"
          aria-label="Vibrational modes table"
        >
          <thead className="bg-slate-50 border-b border-slate-200 sticky top-0">
            <tr>
              {COLUMNS.map((col) => {
                const isActive = sortKey === col.key;
                const arrow = isActive
                  ? sortDirection === 'asc'
                    ? '↑'
                    : '↓'
                  : '';
                return (
                  <th
                    key={col.key}
                    scope="col"
                    className={`px-3 py-2 font-semibold text-slate-700 ${
                      col.align === 'right' ? 'text-right' : 'text-left'
                    }`}
                  >
                    <button
                      type="button"
                      onClick={() => handleSort(col.key)}
                      aria-label={`Sort by ${col.ariaLabel}`}
                      aria-sort={
                        isActive
                          ? sortDirection === 'asc'
                            ? 'ascending'
                            : 'descending'
                          : 'none'
                      }
                      className="hover:text-blue-600 focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 rounded"
                    >
                      {col.label} {arrow}
                    </button>
                  </th>
                );
              })}
            </tr>
          </thead>
          <tbody>
            {sortedRows.map((row) => {
              const isSelected = selectedMode === row.modeIndex;
              const classes = [
                'border-b border-slate-100 cursor-pointer transition-colors',
                isSelected
                  ? 'bg-blue-100 border-l-4 border-l-blue-500'
                  : row.isImaginary
                    ? 'bg-red-50 hover:bg-red-100'
                    : 'hover:bg-slate-50',
                row.isImaginary ? 'text-red-600 font-semibold' : 'text-slate-700',
              ].join(' ');
              return (
                <tr
                  key={row.modeIndex}
                  className={classes}
                  onClick={() => onSelectMode(row.modeIndex)}
                  onKeyDown={(e) => {
                    if (e.key === 'Enter' || e.key === ' ') {
                      e.preventDefault();
                      onSelectMode(row.modeIndex);
                    }
                  }}
                  tabIndex={0}
                  role="row"
                  aria-selected={isSelected}
                  aria-label={
                    row.isImaginary
                      ? `Mode ${row.modeIndex + 1}, imaginary frequency ${formatFrequency(row.frequency)} cm⁻¹`
                      : `Mode ${row.modeIndex + 1}, frequency ${formatFrequency(row.frequency)} cm⁻¹`
                  }
                >
                  {COLUMNS.map((col) => (
                    <td
                      key={col.key}
                      className={`px-3 py-2 tabular-nums ${
                        col.align === 'right' ? 'text-right' : 'text-left'
                      }`}
                    >
                      {col.key === 'frequency' && row.isImaginary ? (
                        <span>
                          {col.format(row)}{' '}
                          <em className="not-italic text-[10px] text-red-500 ml-1">
                            imaginary
                          </em>
                        </span>
                      ) : (
                        col.format(row)
                      )}
                    </td>
                  ))}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}

export default FrequencyPanel;
