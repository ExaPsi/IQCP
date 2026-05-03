/**
 * OverlapDistanceControls - Controls for the overlap vs. distance feature.
 *
 * Provides:
 * - "Overlap vs Distance" toggle button to enter/exit overlap mode
 * - Atom B element selector (defaults to same element as atom A)
 * - Shell B selector (defaults to same shell type as selected shell)
 * - "Add Pair" button to overlay additional traces (up to 4)
 * - "Clear Pairs" button to remove all overlay traces
 *
 * The primary pair always uses atom A's current element/basis/shell.
 * Overlap mode and comparison mode are mutually exclusive.
 *
 * @module components/basisExplorer/OverlapDistanceControls
 * @see US-054 Overlap vs. Distance Plot
 */

import { useState, useMemo, memo, useCallback } from 'react';
import { ELEMENTS, BASIS_SETS } from './basisData';
import { useBasisInfo } from '../../hooks/useBasisInfo';
import { deriveShellDisplayInfo } from './shellDisplayInfo';

// ============================================================================
// Types
// ============================================================================

/**
 * Configuration for one overlap trace (the "B" side of a pair).
 */
export interface OverlapTraceConfig {
  /** Atomic number for the second atom */
  elementB: number;
  /** Basis set name for the second atom */
  basisB: string;
  /** Shell index for the second atom */
  shellIndexB: number;
}

export interface OverlapDistanceControlsProps {
  /** Whether overlap distance mode is active */
  isActive: boolean;
  /** Toggle overlap distance mode on/off */
  onToggle: () => void;
  /** Current overlap traces (B-side configurations) */
  overlapTraces: OverlapTraceConfig[];
  /** Add a new overlap trace */
  onAddTrace: (trace: OverlapTraceConfig) => void;
  /** Remove a trace by index */
  onRemoveTrace: (index: number) => void;
  /** Clear all traces */
  onClearTraces: () => void;
  /** Primary element (atom A) from current selection */
  primaryElement: number;
  /** Primary basis set (atom A) from current selection */
  primaryBasis: string;
  /** Primary shell index (atom A) from current selection */
  primaryShell: number | null;
  /** Whether a shell is selected (required for overlap mode) */
  hasShellSelected: boolean;
  /** Whether comparison mode is active (mutually exclusive) */
  comparisonModeActive: boolean;
}

// ============================================================================
// Constants
// ============================================================================

/** Maximum number of overlay traces (including the default primary pair) */
const MAX_TRACES = 4;

// ============================================================================
// Add Pair Form Sub-component
// ============================================================================

/**
 * Form for selecting the B-side of a new overlap pair.
 */
function AddPairForm({
  defaultElement,
  defaultBasis,
  onAdd,
  disabled,
}: {
  defaultElement: number;
  defaultBasis: string;
  onAdd: (config: OverlapTraceConfig) => void;
  disabled: boolean;
}) {
  const [elementB, setElementB] = useState(defaultElement);
  const [basisB, setBasisB] = useState(defaultBasis);
  const [shellIndexB, setShellIndexB] = useState(0);

  // Fetch shell info for the selected element B and basis B
  const { shells: shellsB } = useBasisInfo(elementB, basisB);

  // Derive display labels for B-side shells
  const shellLabelsB = useMemo(() => {
    if (!shellsB) return [];
    return deriveShellDisplayInfo(shellsB, elementB);
  }, [shellsB, elementB]);

  // Reset shell index when element or basis changes
  const handleElementChange = useCallback((z: number) => {
    setElementB(z);
    setShellIndexB(0);
  }, []);

  const handleBasisChange = useCallback((b: string) => {
    setBasisB(b);
    setShellIndexB(0);
  }, []);

  const handleAdd = useCallback(() => {
    onAdd({ elementB, basisB, shellIndexB });
  }, [elementB, basisB, shellIndexB, onAdd]);

  // Get supported elements
  const supportedElements = useMemo(
    () => ELEMENTS.filter((el) => el.supported),
    [],
  );

  return (
    <div className="space-y-2 bg-slate-50 rounded-lg p-3 border border-slate-200">
      <div className="text-xs font-medium text-slate-600 mb-1">
        Add comparison pair (atom B)
      </div>

      {/* Element B selector */}
      <div className="grid grid-cols-2 gap-2">
        <div>
          <label
            htmlFor="overlap-element-b"
            className="block text-xs text-slate-500 mb-0.5"
          >
            Element B
          </label>
          <select
            id="overlap-element-b"
            value={elementB}
            onChange={(e) => handleElementChange(parseInt(e.target.value, 10))}
            className="w-full rounded-md border border-slate-300 px-2 py-1 text-sm
              focus:outline-none focus:ring-2 focus:ring-primary-500 focus:border-primary-500"
          >
            {supportedElements.map((el) => (
              <option key={el.z} value={el.z}>
                {el.symbol} ({el.name})
              </option>
            ))}
          </select>
        </div>

        {/* Basis B selector */}
        <div>
          <label
            htmlFor="overlap-basis-b"
            className="block text-xs text-slate-500 mb-0.5"
          >
            Basis B
          </label>
          <select
            id="overlap-basis-b"
            value={basisB}
            onChange={(e) => handleBasisChange(e.target.value)}
            className="w-full rounded-md border border-slate-300 px-2 py-1 text-sm
              focus:outline-none focus:ring-2 focus:ring-primary-500 focus:border-primary-500"
          >
            {BASIS_SETS.map((bs) => (
              <option key={bs.value} value={bs.value}>
                {bs.label}
              </option>
            ))}
          </select>
        </div>
      </div>

      {/* Shell B selector */}
      <div>
        <label
          htmlFor="overlap-shell-b"
          className="block text-xs text-slate-500 mb-0.5"
        >
          Shell B
        </label>
        <select
          id="overlap-shell-b"
          value={shellIndexB}
          onChange={(e) => setShellIndexB(parseInt(e.target.value, 10))}
          disabled={shellLabelsB.length === 0}
          className="w-full rounded-md border border-slate-300 px-2 py-1 text-sm
            focus:outline-none focus:ring-2 focus:ring-primary-500 focus:border-primary-500
            disabled:bg-slate-100 disabled:text-slate-400"
        >
          {shellLabelsB.length === 0 ? (
            <option value={0}>Loading...</option>
          ) : (
            shellLabelsB.map((si, idx) => (
              <option key={idx} value={idx}>
                {si.label} ({si.angularMomentumLabel}-type, {si.nPrimitives} primitives)
              </option>
            ))
          )}
        </select>
      </div>

      {/* Add button */}
      <button
        type="button"
        onClick={handleAdd}
        disabled={disabled || shellLabelsB.length === 0}
        className="w-full px-3 py-1.5 text-xs font-medium rounded-md transition-colors
          focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
          bg-primary-50 text-primary-700 border border-primary-200 hover:bg-primary-100
          disabled:bg-slate-100 disabled:text-slate-400 disabled:cursor-not-allowed disabled:border-slate-200"
      >
        Add Pair
      </button>
    </div>
  );
}

// ============================================================================
// Main Component
// ============================================================================

/**
 * Controls panel for overlap vs. distance mode.
 *
 * Renders a toggle button and, when active, controls for selecting
 * the B-side of overlap pairs and managing overlay traces.
 */
export const OverlapDistanceControls = memo(function OverlapDistanceControls({
  isActive,
  onToggle,
  overlapTraces,
  onAddTrace,
  onRemoveTrace,
  onClearTraces,
  primaryElement,
  primaryBasis,
  primaryShell,
  hasShellSelected,
  comparisonModeActive,
}: OverlapDistanceControlsProps) {
  const atMaxTraces = overlapTraces.length >= MAX_TRACES;
  const canActivate = hasShellSelected && !comparisonModeActive;

  // Element symbol lookup
  const getElementSymbol = useCallback((z: number) => {
    const el = ELEMENTS.find((e) => e.z === z);
    return el?.symbol ?? `Z${z}`;
  }, []);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      {/* Toggle button */}
      <div className="flex items-center justify-between mb-3">
        <h3 className="text-sm font-semibold text-slate-700">
          Overlap vs. Distance
        </h3>
        <button
          type="button"
          onClick={onToggle}
          disabled={!canActivate}
          aria-pressed={isActive}
          title={
            comparisonModeActive
              ? 'Exit comparison mode first to use overlap distance'
              : !hasShellSelected
                ? 'Select a shell to enable overlap distance'
                : isActive
                  ? 'Exit overlap distance mode'
                  : 'View overlap integral vs interatomic distance'
          }
          className={`px-3 py-1.5 text-xs font-medium rounded-lg transition-colors
            focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
            ${
              !canActivate
                ? 'bg-slate-100 text-slate-400 cursor-not-allowed'
                : isActive
                  ? 'bg-primary-600 text-white hover:bg-primary-700'
                  : 'bg-primary-50 text-primary-700 border border-primary-200 hover:bg-primary-100'
            }`}
        >
          {isActive ? 'Exit Overlap' : 'Overlap Plot'}
        </button>
      </div>

      {/* Hint when disabled */}
      {!canActivate && (
        <p className="text-xs text-slate-400 italic">
          {comparisonModeActive
            ? 'Exit comparison mode to use overlap distance mode.'
            : 'Select a shell from the list above to enable overlap distance.'}
        </p>
      )}

      {/* Controls when active */}
      {isActive && hasShellSelected && primaryShell !== null && (
        <div className="space-y-3">
          {/* Current primary pair info */}
          <div className="bg-blue-50 rounded-lg px-3 py-2 border border-blue-100">
            <p className="text-xs text-blue-800">
              <strong>Primary pair (atom A):</strong>{' '}
              {getElementSymbol(primaryElement)} shell #{primaryShell}{' '}
              ({primaryBasis.toUpperCase()})
            </p>
          </div>

          {/* Existing overlay traces */}
          {overlapTraces.length > 0 && (
            <div className="space-y-1">
              <div className="flex items-center justify-between">
                <span className="text-xs font-medium text-slate-600">
                  Overlay traces ({overlapTraces.length}/{MAX_TRACES})
                </span>
                <button
                  type="button"
                  onClick={onClearTraces}
                  className="text-xs text-red-600 hover:text-red-700 transition-colors
                    focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-red-500"
                >
                  Clear all
                </button>
              </div>
              {overlapTraces.map((trace, idx) => (
                <div
                  key={idx}
                  className="flex items-center justify-between bg-slate-50 rounded-md px-3 py-1.5"
                >
                  <span className="text-xs text-slate-700">
                    {getElementSymbol(trace.elementB)} shell #{trace.shellIndexB}{' '}
                    ({trace.basisB.toUpperCase()})
                  </span>
                  <button
                    type="button"
                    onClick={() => onRemoveTrace(idx)}
                    className="text-xs text-slate-400 hover:text-red-500 transition-colors
                      focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-red-500"
                    aria-label={`Remove trace ${idx + 1}`}
                  >
                    Remove
                  </button>
                </div>
              ))}
            </div>
          )}

          {/* Add pair form */}
          {!atMaxTraces && (
            <AddPairForm
              defaultElement={primaryElement}
              defaultBasis={primaryBasis}
              onAdd={onAddTrace}
              disabled={atMaxTraces}
            />
          )}

          {atMaxTraces && (
            <p className="text-xs text-slate-400 italic text-center">
              Maximum {MAX_TRACES} traces reached. Remove a trace to add another.
            </p>
          )}
        </div>
      )}
    </div>
  );
});

export default OverlapDistanceControls;
