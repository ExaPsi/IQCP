/**
 * EriBrowser - Interactive browser for two-electron repulsion integrals.
 *
 * Provides a searchable quartet selector with four numeric inputs for
 * (i, j, k, l) indices, a "Top ERIs" quick-select dropdown showing
 * the 20 most significant ERIs by magnitude, and an "Inspect" button
 * to fetch the full primitive-quartet breakdown.
 *
 * The top-N list is computed entirely in JavaScript from the cached
 * eriCompressed array (no WASM call needed). Only the detailed
 * breakdown triggers a WASM worker call.
 *
 * @module components/integrals/EriBrowser
 * @see US-059 ERI Browser
 * @see FR-INT-04
 */

import { useMemo, useCallback, memo, useState, useEffect } from 'react';
import { useIntegralStore } from '../../stores/integralStore';
import type { EriCacheState } from '../../stores/integralStore';
import { useWorker } from '../../hooks/useWorker';
import { useEriDetail } from '../../hooks/useEriDetail';
import type { EriDetailResult } from '../../hooks/useEriDetail';
import { EriQuartetDetail } from './EriQuartetDetail';
import type { GeometryInput } from '../../worker/protocol';

// ============================================================================
// Types
// ============================================================================

export interface EriBrowserProps {
  /** Current geometry for ERI detail requests */
  geometry: GeometryInput;
  /** Current basis set name */
  basisName: string;
  /** Number of basis functions */
  nbf: number;
  /** ERI cache data (compressed flat array with 8-fold symmetry) */
  eriCache: Extract<EriCacheState, { status: 'success' }>;
}

/**
 * An ERI value with its indices, used for the top-N dropdown.
 */
interface EriEntry {
  indices: [number, number, number, number];
  value: number;
  absMagnitude: number;
}

// ============================================================================
// ERI Index Utilities
// ============================================================================

/**
 * Compute the index into the 8-fold symmetry compressed ERI array.
 *
 * Uses the standard pair index formula:
 *   pair(i,j) = max(i,j) * (max(i,j) + 1) / 2 + min(i,j)
 *   eri_index = max(pair_ij, pair_kl) * (max(pair_ij, pair_kl) + 1) / 2 + min(pair_ij, pair_kl)
 */
function eriIndex(i: number, j: number, k: number, l: number): number {
  const maxIJ = Math.max(i, j);
  const minIJ = Math.min(i, j);
  const pairIJ = (maxIJ * (maxIJ + 1)) / 2 + minIJ;

  const maxKL = Math.max(k, l);
  const minKL = Math.min(k, l);
  const pairKL = (maxKL * (maxKL + 1)) / 2 + minKL;

  const maxPair = Math.max(pairIJ, pairKL);
  const minPair = Math.min(pairIJ, pairKL);
  return (maxPair * (maxPair + 1)) / 2 + minPair;
}

/**
 * Extract the value of ERI (i,j,k,l) from compressed storage.
 */
function eriGet(
  eriCompressed: number[],
  i: number,
  j: number,
  k: number,
  l: number,
): number {
  const idx = eriIndex(i, j, k, l);
  return eriCompressed[idx] ?? 0;
}

/**
 * Build the top-N ERIs by magnitude from the compressed array.
 *
 * Iterates over all unique (i,j,k,l) quartets with canonical ordering
 * (ij >= kl, i >= j, k >= l) and sorts by absolute magnitude.
 */
function computeTopEris(eriCompressed: number[], nbf: number, topN: number): EriEntry[] {
  const entries: EriEntry[] = [];

  for (let i = 0; i < nbf; i++) {
    for (let j = 0; j <= i; j++) {
      const pairIJ = (i * (i + 1)) / 2 + j;
      for (let k = 0; k < nbf; k++) {
        for (let l = 0; l <= k; l++) {
          const pairKL = (k * (k + 1)) / 2 + l;
          if (pairIJ < pairKL) continue; // only canonical ordering

          const idx = eriIndex(i, j, k, l);
          const value = eriCompressed[idx] ?? 0;
          entries.push({
            indices: [i, j, k, l],
            value,
            absMagnitude: Math.abs(value),
          });
        }
      }
    }
  }

  // Sort by magnitude descending
  entries.sort((a, b) => b.absMagnitude - a.absMagnitude);

  return entries.slice(0, topN);
}

// ============================================================================
// Constants
// ============================================================================

/** Number of top ERIs to show in the dropdown */
const TOP_N = 20;

// ============================================================================
// Sub-components
// ============================================================================

/**
 * Four numeric inputs for quartet index selection.
 */
const QuartetInput = memo(function QuartetInput({
  values,
  maxIndex,
  onChange,
}: {
  values: [number, number, number, number];
  maxIndex: number;
  onChange: (quartet: [number, number, number, number]) => void;
}) {
  const labels = ['i', 'j', 'k', 'l'] as const;

  const handleChange = useCallback(
    (pos: number, val: string) => {
      const num = parseInt(val, 10);
      if (isNaN(num) || num < 0 || num >= maxIndex) return;
      const next = [...values] as [number, number, number, number];
      next[pos] = num;
      onChange(next);
    },
    [values, maxIndex, onChange],
  );

  return (
    <div className="flex items-center gap-1">
      <span className="text-xs text-slate-500 mr-1">(</span>
      {labels.map((label, idx) => (
        <div key={label} className="flex items-center">
          {idx === 2 && <span className="text-xs text-slate-400 mx-0.5">|</span>}
          <div className="flex flex-col items-center">
            <label
              htmlFor={`eri-idx-${label}`}
              className="text-[9px] text-slate-400 font-mono mb-0.5"
            >
              {label}
            </label>
            <input
              id={`eri-idx-${label}`}
              type="number"
              min={0}
              max={maxIndex - 1}
              value={values[idx]}
              onChange={(e) => handleChange(idx, e.target.value)}
              className="w-10 px-1 py-1 text-xs text-center font-mono border border-slate-300 rounded bg-white
                focus:outline-none focus:ring-1 focus:ring-primary-500 focus:border-primary-500"
              aria-label={`Basis function index ${label}`}
            />
          </div>
        </div>
      ))}
      <span className="text-xs text-slate-500 ml-1">)</span>
    </div>
  );
});

// ============================================================================
// Main Component
// ============================================================================

/**
 * Interactive ERI browser with quartet search, top-N dropdown, and
 * detailed breakdown panel.
 *
 * Features:
 * - Four numeric inputs for (i,j,k,l) indices
 * - Quick-select dropdown of top-20 ERIs by magnitude
 * - "Inspect" button triggers WASM breakdown computation
 * - Displays EriQuartetDetail when data is available
 * - Inline ERI value preview before inspection
 */
export const EriBrowser = memo(function EriBrowser({
  geometry,
  basisName,
  nbf,
  eriCache,
}: EriBrowserProps) {
  const { isReady } = useWorker();

  // Store state
  const eriSelectedQuartet = useIntegralStore((s) => s.eriSelectedQuartet);
  const eriShowAllPrimitives = useIntegralStore((s) => s.eriShowAllPrimitives);
  const setEriSelectedQuartet = useIntegralStore((s) => s.setEriSelectedQuartet);
  const toggleEriShowAllPrimitives = useIntegralStore((s) => s.toggleEriShowAllPrimitives);

  // Local state for the input fields (allows editing before committing)
  const [inputValues, setInputValues] = useState<[number, number, number, number]>([0, 0, 0, 0]);

  // Sync input values from store selection
  useEffect(() => {
    if (eriSelectedQuartet) {
      setInputValues(eriSelectedQuartet);
    }
  }, [eriSelectedQuartet]);

  // Compute top-N ERIs from the compressed cache
  const topEris = useMemo(() => {
    return computeTopEris(eriCache.eriCompressed, nbf, TOP_N);
  }, [eriCache.eriCompressed, nbf]);

  // Preview value for the currently entered quartet
  const previewValue = useMemo(() => {
    const [i, j, k, l] = inputValues;
    if (i < 0 || i >= nbf || j < 0 || j >= nbf || k < 0 || k >= nbf || l < 0 || l >= nbf) {
      return null;
    }
    return eriGet(eriCache.eriCompressed, i, j, k, l);
  }, [inputValues, eriCache.eriCompressed, nbf]);

  // ERI detail hook
  const eriDetail = useEriDetail({
    geometry,
    basisName,
    quartet: eriSelectedQuartet,
    isReady,
  });

  // Handle quartet selection from input fields
  const handleInspect = useCallback(() => {
    const [i, j, k, l] = inputValues;
    if (i >= 0 && i < nbf && j >= 0 && j < nbf && k >= 0 && k < nbf && l >= 0 && l < nbf) {
      setEriSelectedQuartet([i, j, k, l]);
    }
  }, [inputValues, nbf, setEriSelectedQuartet]);

  // Trigger compute when quartet changes (via store) and is valid
  useEffect(() => {
    if (eriSelectedQuartet && isReady) {
      // Small delay so state is consistent before compute
      const timer = setTimeout(() => {
        eriDetail.compute();
      }, 50);
      return () => clearTimeout(timer);
    }
    // Note: eriDetail.compute is intentionally excluded from deps
    // to avoid infinite loops. The compute function itself is stable
    // but references eriSelectedQuartet, so we rely on the quartet
    // change to trigger this effect.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [eriSelectedQuartet, isReady]);

  // Handle dropdown selection
  const handleDropdownSelect = useCallback(
    (e: React.ChangeEvent<HTMLSelectElement>) => {
      const val = e.target.value;
      if (!val) return;
      const parts = val.split(',').map(Number) as [number, number, number, number];
      if (parts.length === 4 && parts.every((p) => !isNaN(p))) {
        setInputValues(parts);
        setEriSelectedQuartet(parts);
      }
    },
    [setEriSelectedQuartet],
  );

  // Format a dropdown label for an ERI entry
  const formatDropdownLabel = useCallback(
    (entry: EriEntry, rank: number): string => {
      const [i, j, k, l] = entry.indices;
      const sign = entry.value >= 0 ? '+' : '';
      return `#${rank}: (${i}${j}|${k}${l}) = ${sign}${entry.value.toFixed(6)}`;
    },
    [],
  );

  return (
    <div className="space-y-4">
      {/* Header */}
      <div>
        <h2 className="text-lg font-bold text-slate-800">
          ERI Browser
        </h2>
        <p className="text-xs text-slate-600 mt-1">
          Browse two-electron repulsion integrals (ij|kl). Select a quartet
          by entering indices or choosing from the significant ERIs list,
          then inspect the primitive-quartet breakdown.
        </p>
      </div>

      {/* Controls row */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
        <div className="flex flex-col sm:flex-row gap-4">
          {/* Quartet input */}
          <div className="flex-1">
            <h3 className="text-xs font-semibold text-slate-400 uppercase tracking-wider mb-2">
              Select Quartet
            </h3>
            <div className="flex items-end gap-3">
              <QuartetInput
                values={inputValues}
                maxIndex={nbf}
                onChange={setInputValues}
              />
              <button
                onClick={handleInspect}
                disabled={!isReady || previewValue === null}
                className="px-4 py-1.5 text-xs font-medium bg-primary-600 text-white rounded-lg
                  hover:bg-primary-700 disabled:bg-slate-300 disabled:cursor-not-allowed
                  transition-colors focus-visible:outline-none focus-visible:ring-2
                  focus-visible:ring-primary-500 focus-visible:ring-offset-1"
              >
                Inspect
              </button>
            </div>
            {previewValue !== null && (
              <p className="text-[10px] text-slate-500 mt-1.5 font-mono">
                ({inputValues[0]}{inputValues[1]}|{inputValues[2]}{inputValues[3]}) = {previewValue.toExponential(6)}
              </p>
            )}
          </div>

          {/* Top ERIs dropdown */}
          <div className="flex-1">
            <h3 className="text-xs font-semibold text-slate-400 uppercase tracking-wider mb-2">
              Significant ERIs (Top {topEris.length})
            </h3>
            <select
              onChange={handleDropdownSelect}
              value={eriSelectedQuartet ? eriSelectedQuartet.join(',') : ''}
              className="w-full px-3 py-1.5 text-xs font-mono border border-slate-300 rounded-lg
                bg-white focus:outline-none focus:ring-2 focus:ring-primary-500 focus:border-transparent"
              aria-label="Select a significant ERI by magnitude"
            >
              <option value="">-- Select an ERI --</option>
              {topEris.map((entry, idx) => {
                const key = entry.indices.join(',');
                return (
                  <option key={key} value={key}>
                    {formatDropdownLabel(entry, idx + 1)}
                  </option>
                );
              })}
            </select>
            <p className="text-[10px] text-slate-400 mt-1">
              {eriCache.eriCount} unique ERIs (8-fold symmetry), computed in{' '}
              {eriCache.computeTimeMs.toFixed(1)} ms
            </p>
          </div>
        </div>
      </div>

      {/* Loading state */}
      {eriDetail.loading && (
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <div className="flex items-center justify-center gap-3">
            <div className="animate-spin rounded-full h-5 w-5 border-b-2 border-primary-600" />
            <span className="text-sm text-slate-600">
              Computing ERI breakdown...
            </span>
          </div>
        </div>
      )}

      {/* Error state */}
      {eriDetail.error && !eriDetail.loading && (
        <div className="bg-white rounded-xl shadow-sm border border-red-200 p-5">
          <div className="flex items-start gap-3">
            <div className="flex-shrink-0 w-7 h-7 bg-red-100 rounded-full flex items-center justify-center">
              <span className="text-red-600 text-xs font-bold">!</span>
            </div>
            <div className="flex-1">
              <h3 className="text-sm font-semibold text-red-800">ERI Detail Error</h3>
              <p className="text-sm text-red-700 mt-1">{eriDetail.error}</p>
            </div>
          </div>
        </div>
      )}

      {/* ERI detail panel */}
      {eriDetail.data && !eriDetail.loading && (
        <EriQuartetDetail
          data={eriDetail.data as EriDetailResult}
          showAll={eriShowAllPrimitives}
          onToggleShowAll={toggleEriShowAllPrimitives}
        />
      )}
    </div>
  );
});
