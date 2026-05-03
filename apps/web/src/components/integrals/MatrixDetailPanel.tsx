/**
 * MatrixDetailPanel - Detail panel for a selected matrix element.
 *
 * Shows the value of the selected (i, j) element along with:
 * - Basis function labels for row i and column j
 * - Numerical value with full precision
 * - Physical interpretation text for the matrix type
 * - Primitive breakdown table (US-056)
 *
 * @module components/integrals/MatrixDetailPanel
 * @see US-055 Integral Matrix Heatmap
 * @see US-056 Primitive Breakdown Panel
 */

import type { MatrixTab } from '../../stores/integralStore';
import { PrimitiveBreakdownTable } from './PrimitiveBreakdownTable';
import type { IntegralBreakdownResult } from './PrimitiveBreakdownTable';

// ============================================================================
// Types
// ============================================================================

export interface MatrixDetailPanelProps {
  /** Row index */
  row: number;
  /** Column index */
  col: number;
  /** Basis function label for row */
  rowLabel: string;
  /** Basis function label for column */
  colLabel: string;
  /** Matrix element value */
  value: number;
  /** Matrix type for interpretation text */
  matrixType: MatrixTab;
  /** Primitive breakdown data (null if not yet computed) */
  breakdownData?: IntegralBreakdownResult | null;
  /** Whether breakdown is loading */
  breakdownLoading?: boolean;
  /** Breakdown computation error */
  breakdownError?: string | null;
  /** Whether all primitive pairs are shown (vs. top-N) */
  showAllPrimitives?: boolean;
  /** Callback to toggle show all / top-N */
  onToggleShowAll?: () => void;
}

// ============================================================================
// Interpretation Text
// ============================================================================

/**
 * Get the physical interpretation of a matrix element based on type and indices.
 */
function getInterpretation(
  matrixType: MatrixTab,
  row: number,
  col: number,
  rowLabel: string,
  colLabel: string,
  value: number,
): string {
  const isDiagonal = row === col;
  const subscript = `${row + 1},${col + 1}`;

  switch (matrixType) {
    case 'S':
      if (isDiagonal) {
        return `S_${subscript} = 1.0: The overlap of ${rowLabel} with itself is always 1 for normalized basis functions.`;
      }
      if (Math.abs(value) < 1e-10) {
        return `S_${subscript} = 0: ${rowLabel} and ${colLabel} are orthogonal -- they have no spatial overlap. This often occurs between functions of different symmetry (e.g., s and p_y in the molecular plane).`;
      }
      return `S_${subscript}: Overlap between ${rowLabel} and ${colLabel}. A value of ${value.toFixed(4)} indicates ${Math.abs(value) > 0.5 ? 'significant' : 'moderate'} spatial overlap between these basis functions.`;

    case 'T':
      if (isDiagonal) {
        return `T_${subscript}: Kinetic energy of ${rowLabel}. Positive because kinetic energy is always positive for any wavefunction.`;
      }
      return `T_${subscript}: Kinetic energy coupling between ${rowLabel} and ${colLabel}. Arises from the second derivative in the kinetic energy operator.`;

    case 'V':
      if (isDiagonal) {
        return `V_${subscript}: Nuclear attraction of ${rowLabel}. Negative because the Coulomb potential between electrons and nuclei is attractive.`;
      }
      return `V_${subscript}: Nuclear attraction between ${rowLabel} and ${colLabel}. Negative values reflect the attractive electron-nuclear interaction. Larger magnitudes indicate stronger attraction.`;

    case 'Hcore':
      if (isDiagonal) {
        return `H^core_${subscript} = T_${subscript} + V_${subscript}: One-electron energy of ${rowLabel}. Typically negative because nuclear attraction (V) dominates kinetic energy (T).`;
      }
      return `H^core_${subscript} = T_${subscript} + V_${subscript}: One-electron coupling between ${rowLabel} and ${colLabel}. This is the one-electron part of the Fock matrix, combining kinetic energy and nuclear attraction.`;
  }
}

// ============================================================================
// Component
// ============================================================================

/**
 * Detail panel showing information about a selected matrix element.
 *
 * Provides the numerical value with full precision, basis function labels,
 * and a physical interpretation to help students understand what each
 * matrix element represents.
 */
export function MatrixDetailPanel({
  row,
  col,
  rowLabel,
  colLabel,
  value,
  matrixType,
  breakdownData,
  breakdownLoading,
  breakdownError,
  showAllPrimitives = false,
  onToggleShowAll,
}: MatrixDetailPanelProps) {
  const interpretation = getInterpretation(matrixType, row, col, rowLabel, colLabel, value);
  const isDiagonal = row === col;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
      <h3 className="text-sm font-semibold text-slate-700 mb-3">
        Selected Element
      </h3>

      {/* Element identification */}
      <div className="grid grid-cols-2 gap-3 mb-4">
        <div className="bg-slate-50 rounded-lg px-3 py-2">
          <p className="text-xs text-slate-500">Row ({row + 1})</p>
          <p className="text-sm font-medium text-slate-800">{rowLabel}</p>
        </div>
        <div className="bg-slate-50 rounded-lg px-3 py-2">
          <p className="text-xs text-slate-500">Column ({col + 1})</p>
          <p className="text-sm font-medium text-slate-800">{colLabel}</p>
        </div>
      </div>

      {/* Value display */}
      <div className="bg-slate-50 rounded-lg px-3 py-2 mb-4">
        <p className="text-xs text-slate-500 mb-1">
          {matrixType === 'Hcore' ? 'H^core' : matrixType}
          <sub>{row + 1},{col + 1}</sub>
          {isDiagonal && ' (diagonal)'}
        </p>
        <p className="text-lg font-mono font-semibold text-slate-900">
          {value.toExponential(10)}
        </p>
        <p className="text-xs text-slate-400 font-mono mt-0.5">
          = {value.toFixed(12)}
        </p>
      </div>

      {/* Physical interpretation */}
      <div className="bg-blue-50 rounded-lg px-3 py-2 border border-blue-100 mb-4">
        <p className="text-xs font-medium text-blue-800 mb-1">
          Physical Interpretation
        </p>
        <p className="text-xs text-blue-700 leading-relaxed">
          {interpretation}
        </p>
      </div>

      {/* Symmetry note */}
      {!isDiagonal && (
        <div className="text-xs text-slate-500">
          <p>
            This matrix is symmetric: {matrixType === 'Hcore' ? 'H^core' : matrixType}
            <sub>{row + 1},{col + 1}</sub> = {matrixType === 'Hcore' ? 'H^core' : matrixType}
            <sub>{col + 1},{row + 1}</sub>
          </p>
        </div>
      )}

      {/* Primitive breakdown (US-056) */}
      {breakdownLoading && (
        <div className="mt-4 border-t border-slate-100 pt-3">
          <div className="flex items-center gap-2">
            <div className="animate-spin rounded-full h-3.5 w-3.5 border-b-2 border-primary-600" />
            <p className="text-xs text-slate-500">Loading primitive breakdown...</p>
          </div>
        </div>
      )}
      {breakdownError && !breakdownLoading && (
        <div className="mt-4 border-t border-slate-100 pt-3">
          <p className="text-xs text-red-600">
            Breakdown error: {breakdownError}
          </p>
        </div>
      )}
      {breakdownData && !breakdownLoading && (
        <div className="mt-4 border-t border-slate-100 pt-3">
          <PrimitiveBreakdownTable
            breakdown={breakdownData}
            showAll={showAllPrimitives}
            onToggleShowAll={onToggleShowAll ?? (() => {})}
          />
        </div>
      )}
    </div>
  );
}
