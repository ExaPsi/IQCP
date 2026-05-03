/**
 * DftComparisonPanel - HF vs DFT side-by-side energy comparison.
 *
 * Provides a "Compare RHF vs DFT" button that triggers sequential
 * computations: first RHF, then the selected DFT method. Results
 * are displayed in a compact table showing energy, iterations,
 * convergence status, and energy difference.
 *
 * Only available when a DFT method is selected (not RHF).
 *
 * @module components/scf/DftComparisonPanel
 * @see US-070 DFT UI + Method Selector + Deep Links (AC4)
 */

import { useScfStore } from '../../stores/scfStore';
import type { DftCompareState } from '../../stores/scfStore';
import { getDftMethodInfo } from '../../types/dft';

/**
 * Props for the DftComparisonPanel component.
 */
export interface DftComparisonPanelProps {
  /** Start a comparison run */
  onCompare: () => void;
  /** Whether the worker is ready */
  isReady: boolean;
  /** Whether a computation is currently running */
  isRunning: boolean;
}

/**
 * Format an energy value compactly.
 */
function fmtE(e: number): string {
  return e.toFixed(8);
}

/**
 * Format energy difference with sign and color indicator.
 */
function DeltaE({ value }: { value: number }) {
  const sign = value < 0 ? '' : '+';
  const color = value < 0 ? 'text-green-700' : value > 0 ? 'text-red-700' : 'text-slate-700';
  return (
    <span className={`font-mono ${color}`}>
      {sign}{value.toFixed(8)}
    </span>
  );
}

/**
 * Phase progress indicator.
 */
function PhaseIndicator({ phase }: { phase: DftCompareState['phase'] }) {
  switch (phase) {
    case 'running_rhf':
      return (
        <div className="flex items-center gap-2 text-sm text-blue-600">
          <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-600" />
          Running RHF...
        </div>
      );
    case 'running_dft':
      return (
        <div className="flex items-center gap-2 text-sm text-teal-600">
          <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-teal-600" />
          Running DFT...
        </div>
      );
    default:
      return null;
  }
}

/**
 * DFT vs HF comparison panel with sequential computation and results table.
 *
 * @example
 * ```tsx
 * <DftComparisonPanel
 *   onCompare={handleCompare}
 *   isReady={isReady}
 *   isRunning={isRunning}
 * />
 * ```
 */
export function DftComparisonPanel({ onCompare, isReady, isRunning }: DftComparisonPanelProps) {
  const method = useScfStore((s) => s.method);
  const compareState = useScfStore((s) => s.dftCompareState);
  const resetDftCompare = useScfStore((s) => s.resetDftCompare);

  const methodInfo = getDftMethodInfo(method);
  const isComparing = compareState.comparing;
  const hasResults = compareState.phase === 'complete';
  const hasError = compareState.phase === 'error';

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      <h3 className="text-sm font-semibold text-slate-700 mb-3">
        RHF vs {methodInfo.label} Comparison
      </h3>

      {/* Compare button */}
      <div className="flex gap-2 mb-3">
        <button
          type="button"
          onClick={onCompare}
          disabled={!isReady || isRunning || isComparing}
          className="flex-1 px-3 py-2 rounded-lg text-sm font-medium bg-teal-600 text-white hover:bg-teal-700 transition-colors disabled:opacity-50 disabled:cursor-not-allowed focus:outline-none focus-visible:ring-2 focus-visible:ring-teal-500"
          aria-label={`Compare RHF vs ${methodInfo.label}`}
        >
          {isComparing ? 'Comparing...' : `Compare RHF vs ${methodInfo.label}`}
        </button>
        {hasResults && (
          <button
            type="button"
            onClick={resetDftCompare}
            className="px-3 py-2 rounded-lg text-sm font-medium bg-slate-100 text-slate-700 hover:bg-slate-200 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-slate-500"
            aria-label="Reset comparison"
          >
            Reset
          </button>
        )}
      </div>

      {/* Progress indicator */}
      {isComparing && <PhaseIndicator phase={compareState.phase} />}

      {/* Error display */}
      {hasError && compareState.error && (
        <div className="flex items-center gap-2 text-sm text-red-600 mb-3">
          <svg className="w-4 h-4 flex-shrink-0" fill="currentColor" viewBox="0 0 20 20">
            <path
              fillRule="evenodd"
              d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z"
              clipRule="evenodd"
            />
          </svg>
          {compareState.error}
        </div>
      )}

      {/* Results table */}
      {hasResults && compareState.rhfResult && compareState.dftResult && (
        <>
          <table className="w-full text-xs border-collapse">
            <thead>
              <tr className="border-b border-slate-200">
                <th className="py-1.5 text-left font-medium text-slate-500"></th>
                <th className="py-1.5 text-right font-medium text-slate-500">RHF</th>
                <th className="py-1.5 text-right font-medium text-slate-500">{methodInfo.label}</th>
                <th className="py-1.5 text-right font-medium text-slate-500">Difference</th>
              </tr>
            </thead>
            <tbody className="text-slate-700">
              <tr className="border-b border-slate-100">
                <td className="py-1.5 font-medium">Energy (Ha)</td>
                <td className="py-1.5 text-right font-mono">{fmtE(compareState.rhfResult.energy)}</td>
                <td className="py-1.5 text-right font-mono">{fmtE(compareState.dftResult.energy)}</td>
                <td className="py-1.5 text-right">
                  <DeltaE value={compareState.dftResult.energy - compareState.rhfResult.energy} />
                </td>
              </tr>
              <tr className="border-b border-slate-100">
                <td className="py-1.5 font-medium">Iterations</td>
                <td className="py-1.5 text-right font-mono">{compareState.rhfResult.iterations}</td>
                <td className="py-1.5 text-right font-mono">{compareState.dftResult.iterations}</td>
                <td className="py-1.5 text-right font-mono text-slate-500">
                  {compareState.dftResult.iterations - compareState.rhfResult.iterations >= 0 ? '+' : ''}
                  {compareState.dftResult.iterations - compareState.rhfResult.iterations}
                </td>
              </tr>
              <tr>
                <td className="py-1.5 font-medium">Converged</td>
                <td className="py-1.5 text-right">
                  {compareState.rhfResult.converged ? (
                    <span className="text-green-600">Yes</span>
                  ) : (
                    <span className="text-red-600">No</span>
                  )}
                </td>
                <td className="py-1.5 text-right">
                  {compareState.dftResult.converged ? (
                    <span className="text-green-600">Yes</span>
                  ) : (
                    <span className="text-red-600">No</span>
                  )}
                </td>
                <td className="py-1.5 text-right text-slate-400">--</td>
              </tr>
            </tbody>
          </table>

          {/* Pedagogical note */}
          <p className="mt-3 text-xs text-slate-500 italic">
            Correlation energy captured by DFT but missing from Hartree-Fock typically
            lowers the total energy. A negative difference indicates DFT recovers
            additional electron correlation.
          </p>
        </>
      )}
    </div>
  );
}

export default DftComparisonPanel;
