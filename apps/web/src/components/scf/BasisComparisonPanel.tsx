/**
 * BasisComparisonPanel - Basis set comparison controls and progress.
 *
 * Provides a collapsible panel for selecting basis sets and running
 * multi-basis SCF comparisons on the same molecule.
 *
 * Features:
 * - Checkbox list for available basis sets (STO-3G, 3-21G, 6-31G, 6-31G*)
 * - Compare button (disabled until >= 2 basis sets are selected)
 * - Cancel button during computation
 * - Progress indicator showing which basis is currently being computed
 *
 * Only shown when the user has entered a custom geometry (custom mode).
 *
 * @module components/scf/BasisComparisonPanel
 * @see US-045 Basis Set Comparison Mode
 */

import { useState, useCallback, useMemo } from 'react';
import { useScfStore } from '../../stores/scfStore';
import { useBasisComparison } from '../../hooks/useBasisComparison';
import type { GeometryInput, BasisSetName } from '../../worker/protocol';

// ============================================================================
// Constants
// ============================================================================

/**
 * Available basis sets for comparison with display labels.
 *
 * Ordered from smallest (minimal) to largest (polarized split-valence).
 */
const AVAILABLE_BASES: { id: BasisSetName; label: string; description: string }[] = [
  { id: 'sto-3g', label: 'STO-3G', description: 'Minimal basis' },
  { id: '3-21g', label: '3-21G', description: 'Split-valence' },
  { id: '6-31g', label: '6-31G', description: 'Split-valence (larger)' },
  { id: '6-31g*', label: '6-31G*', description: 'Polarized split-valence' },
  { id: 'cc-pvdz', label: 'cc-pVDZ', description: 'Correlation-consistent DZ' },
];

// ============================================================================
// Props
// ============================================================================

/**
 * Props for BasisComparisonPanel.
 */
interface BasisComparisonPanelProps {
  /** Geometry input for the current molecule */
  geometry: GeometryInput;
  /** Disable all controls (e.g., worker not ready) */
  disabled?: boolean;
}

// ============================================================================
// Component
// ============================================================================

/**
 * Basis set comparison configuration and progress panel.
 *
 * Rendered in the Module E sidebar when custom geometry is active.
 * Styled consistently with PesScanPanel (white card, slate-50 header).
 *
 * @example
 * ```tsx
 * {inputMode.mode === 'custom' && geometry.atoms.length > 0 && (
 *   <BasisComparisonPanel geometry={geometry} disabled={!isReady} />
 * )}
 * ```
 */
export function BasisComparisonPanel({ geometry, disabled = false }: BasisComparisonPanelProps) {
  // Store state
  const basisCompareState = useScfStore((state) => state.basisCompareState);
  const setBasisCompareSelected = useScfStore((state) => state.setBasisCompareSelected);

  // Hook for comparison execution
  const { isReady, isComparing, compare, cancel, clear } = useBasisComparison();

  // Collapsible panel state
  const [isExpanded, setIsExpanded] = useState(true);

  // Derived state
  const selectedBases = basisCompareState.selectedBases;
  const hasEnoughBases = selectedBases.length >= 2;
  const controlsDisabled = disabled || !isReady || isComparing;

  // Toggle a basis set selection
  const handleToggleBasis = useCallback(
    (basisId: string) => {
      if (controlsDisabled) return;

      const newBases = selectedBases.includes(basisId)
        ? selectedBases.filter((b) => b !== basisId)
        : [...selectedBases, basisId];

      // Limit to 4 basis sets maximum
      if (newBases.length <= 4) {
        setBasisCompareSelected(newBases);
      }
    },
    [selectedBases, controlsDisabled, setBasisCompareSelected]
  );

  // Start comparison
  const handleCompare = useCallback(() => {
    if (!hasEnoughBases || controlsDisabled) return;
    compare(selectedBases, geometry);
  }, [hasEnoughBases, controlsDisabled, compare, selectedBases, geometry]);

  // Reset results
  const handleReset = useCallback(() => {
    clear();
  }, [clear]);

  // Progress info
  const progress = basisCompareState.progress;
  const progressPercent = useMemo(() => {
    if (!progress) return 0;
    return (progress.current / progress.total) * 100;
  }, [progress]);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      {/* Header - clickable to expand/collapse */}
      <button
        type="button"
        onClick={() => setIsExpanded(!isExpanded)}
        className="w-full px-4 py-3 bg-slate-50 border-b border-slate-200 flex items-center justify-between text-left hover:bg-slate-100 transition-colors"
        aria-expanded={isExpanded}
        aria-controls="basis-compare-content"
      >
        <div className="flex items-center gap-2">
          <h3 className="text-sm font-semibold text-slate-700">
            Basis Set Comparison
          </h3>
          {basisCompareState.results.length > 0 && !isComparing && (
            <span className="text-xs text-slate-500">
              ({basisCompareState.results.length} results)
            </span>
          )}
        </div>
        <svg
          className={`w-4 h-4 text-slate-500 transition-transform ${isExpanded ? 'rotate-180' : ''}`}
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {/* Content */}
      {isExpanded && (
        <div id="basis-compare-content" className="p-4 space-y-4">
          {/* Description */}
          <p className="text-xs text-slate-500">
            Compare SCF convergence and energy across different basis sets for the same molecule.
          </p>

          {/* Basis Set Checkboxes */}
          <fieldset disabled={controlsDisabled}>
            <legend className="block text-xs font-medium text-slate-700 mb-2">
              Select basis sets (2-4):
            </legend>
            <div className="space-y-2">
              {AVAILABLE_BASES.map((basis) => {
                const isChecked = selectedBases.includes(basis.id);
                const isDisabledByLimit = !isChecked && selectedBases.length >= 4;

                return (
                  <label
                    key={basis.id}
                    className={`flex items-center gap-2.5 px-2.5 py-1.5 rounded-lg cursor-pointer transition-colors ${
                      isChecked
                        ? 'bg-blue-50 border border-blue-200'
                        : 'border border-transparent hover:bg-slate-50'
                    } ${(controlsDisabled || isDisabledByLimit) ? 'opacity-50 cursor-not-allowed' : ''}`}
                  >
                    <input
                      type="checkbox"
                      checked={isChecked}
                      onChange={() => handleToggleBasis(basis.id)}
                      disabled={controlsDisabled || isDisabledByLimit}
                      className="rounded border-slate-300 text-blue-600 focus:ring-blue-500 h-3.5 w-3.5"
                    />
                    <div className="flex-1 min-w-0">
                      <span className="text-sm font-medium text-slate-800">
                        {basis.label}
                      </span>
                      <span className="text-xs text-slate-500 ml-1.5">
                        {basis.description}
                      </span>
                    </div>
                  </label>
                );
              })}
            </div>
          </fieldset>

          {/* Validation message */}
          {!hasEnoughBases && (
            <p className="text-xs text-amber-600">
              Select at least 2 basis sets to compare.
            </p>
          )}

          {/* Progress Bar (during comparison) */}
          {isComparing && progress && (
            <div className="space-y-1">
              <div className="flex items-center justify-between text-xs text-slate-600">
                <span>
                  Computing{' '}
                  <span className="font-medium">{progress.currentBasis.toUpperCase()}</span>...
                </span>
                <span>{progress.current} / {progress.total}</span>
              </div>
              <div
                className="w-full bg-slate-200 rounded-full h-2 overflow-hidden"
                role="progressbar"
                aria-valuenow={progress.current}
                aria-valuemin={0}
                aria-valuemax={progress.total}
                aria-label={`Basis comparison progress: ${progress.current} of ${progress.total} basis sets`}
              >
                <div
                  className="bg-blue-500 h-full rounded-full transition-all duration-300"
                  style={{ width: `${progressPercent}%` }}
                />
              </div>
            </div>
          )}

          {/* Error Display */}
          {basisCompareState.error && !isComparing && (
            <div
              className="text-xs text-red-600 bg-red-50 border border-red-200 rounded-lg px-3 py-2"
              role="alert"
            >
              {basisCompareState.error}
            </div>
          )}

          {/* Completion Summary */}
          {!isComparing && basisCompareState.results.length > 0 && !basisCompareState.error && (
            <div className="text-xs text-green-700 bg-green-50 border border-green-200 rounded-lg px-3 py-2">
              <div className="flex items-center gap-1.5">
                <svg className="w-3.5 h-3.5 flex-shrink-0" fill="currentColor" viewBox="0 0 20 20">
                  <path
                    fillRule="evenodd"
                    d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z"
                    clipRule="evenodd"
                  />
                </svg>
                <span>
                  Comparison complete: {basisCompareState.results.length} basis sets
                </span>
              </div>
            </div>
          )}

          {/* Action Buttons */}
          <div className="flex gap-2">
            {!isComparing ? (
              <>
                <button
                  type="button"
                  onClick={handleCompare}
                  disabled={!hasEnoughBases || controlsDisabled}
                  className="flex-1 px-3 py-2 rounded-lg text-sm font-medium bg-blue-600 text-white hover:bg-blue-700 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
                  aria-label="Start basis set comparison"
                >
                  Compare
                </button>
                {basisCompareState.results.length > 0 && (
                  <button
                    type="button"
                    onClick={handleReset}
                    disabled={controlsDisabled}
                    className="px-3 py-2 rounded-lg text-sm font-medium bg-slate-100 text-slate-700 hover:bg-slate-200 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-slate-500 disabled:opacity-50 disabled:cursor-not-allowed"
                    aria-label="Reset basis comparison results"
                  >
                    Reset
                  </button>
                )}
              </>
            ) : (
              <button
                type="button"
                onClick={cancel}
                className="flex-1 px-3 py-2 rounded-lg text-sm font-medium bg-red-100 text-red-700 hover:bg-red-200 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-red-500"
                aria-label="Cancel basis comparison"
              >
                <span className="flex items-center justify-center gap-2">
                  <span className="animate-spin rounded-full h-3.5 w-3.5 border-b-2 border-red-600" />
                  Cancel
                </span>
              </button>
            )}
          </div>
        </div>
      )}
    </div>
  );
}

export type { BasisComparisonPanelProps };
