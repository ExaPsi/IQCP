/**
 * BoysResultDisplay - Display area for Boys function computation results.
 *
 * Shows the computed value, method indicator, and explanatory text for
 * single point evaluation, or a chart for sweep view.
 * Handles all computation states (idle, computing, success, error).
 *
 * @module components/boys/BoysResultDisplay
 */

import { useBoysStore } from '../../stores/boysStore';
import { useBoysSweep } from '../../hooks/useBoysSweep';
import { useBoysEvalAll } from '../../hooks/useBoysEvalAll';
import { BoysChart } from './BoysChart';
import { BoysInternalsPanel } from './BoysInternalsPanel';
import { BoysMultiOrderPanel } from './BoysMultiOrderPanel';
import { Math as MathInline } from '../common';
import type { BoysMethod } from '../../worker/protocol';
import {
  getTurnoverPoint,
  getAsymptoticThreshold,
  THEORETICAL_SMALL_T_THRESHOLD,
  getMaxTValue,
  THEORETICAL_FORMULAS,
} from '../../lib/boysConstants';

/**
 * Format a number in scientific notation for display.
 *
 * Shows 12 significant digits for precision while remaining readable.
 */
function formatScientific(value: number): string {
  // For very small values, use exponential
  if (Math.abs(value) < 1e-4 && value !== 0) {
    return value.toExponential(12);
  }
  // For values close to 1, show decimal
  if (Math.abs(value) >= 0.0001 && Math.abs(value) < 1e6) {
    return value.toPrecision(13);
  }
  // For large values, use exponential
  return value.toExponential(12);
}

/**
 * Method badge colors and labels.
 * Note: The label for 'recurrence' may be overridden for m=0 and m=1.
 */
const METHOD_STYLES: Record<BoysMethod, { bg: string; text: string; label: string }> = {
  zero: {
    bg: 'bg-amber-100',
    text: 'text-amber-800',
    label: 'Zero',
  },
  series: {
    bg: 'bg-blue-100',
    text: 'text-blue-800',
    label: 'Series',
  },
  recurrence: {
    bg: 'bg-green-100',
    text: 'text-green-800',
    label: 'Recurrence',
  },
};

/**
 * Get m-aware method style. For m=0, recurrence is shown as "erf (Direct)".
 * For m=1, it's shown as "erf + 1 Step".
 */
function getMethodStyleWithM(method: BoysMethod | undefined, m: number): { bg: string; text: string; label: string } {
  if (!method || !METHOD_STYLES[method]) {
    return { bg: 'bg-gray-100', text: 'text-gray-800', label: 'Unknown' };
  }

  // Special labels for recurrence when m=0 or m=1
  if (method === 'recurrence') {
    if (m === 0) {
      return { bg: 'bg-green-100', text: 'text-green-800', label: 'erf (Direct)' };
    }
    if (m === 1) {
      return { bg: 'bg-green-100', text: 'text-green-800', label: 'erf + 1 Step' };
    }
  }

  return METHOD_STYLES[method];
}

/**
 * Method explanation text with inline math components.
 * Uses m-aware explanations for special cases.
 */
function getMethodExplanationContent(method: BoysMethod | undefined, m: number, T: number): React.ReactNode {
  const turnover = getTurnoverPoint(m);

  switch (method) {
    case 'zero':
      return (
        <>
          <MathInline>T</MathInline> is effectively zero (<MathInline>{`< 10^{-14}`}</MathInline>),
          so <MathInline>{`F_m(0) = \\frac{1}{2m+1}`}</MathInline> is used directly.
        </>
      );
    case 'series':
      return (
        <>
          <MathInline>{`T = ${T.toFixed(2)}`}</MathInline> is below the turnover point
          (<MathInline>{`${turnover.toFixed(2)}`}</MathInline>), so series expansion with
          <strong> downward recurrence</strong> is used. This is stable for small <MathInline>T</MathInline>.
        </>
      );
    case 'recurrence':
      // Special explanations for m=0 and m=1
      if (m === 0) {
        return (
          <>
            For <MathInline>m=0</MathInline>, <MathInline>{`F_0`}</MathInline> is computed
            <strong> directly</strong> from the error function. No recurrence needed!
          </>
        );
      }
      if (m === 1) {
        return (
          <>
            For <MathInline>m=1</MathInline>, compute <MathInline>{`F_0`}</MathInline> from erf,
            then just <strong>one recurrence step</strong> to get <MathInline>{`F_1`}</MathInline>.
            This is always stable.
          </>
        );
      }
      // General case: m >= 2
      return (
        <>
          <MathInline>{`T = ${T.toFixed(2)}`}</MathInline> is at or above the turnover point
          (<MathInline>{`${turnover.toFixed(2)}`}</MathInline>), so the erf-based formula with
          <strong> upward recurrence</strong> is used. This is stable when <MathInline>T</MathInline> is large enough.
        </>
      );
    default:
      return 'Unknown computational method.';
  }
}



/**
 * Loading spinner component.
 */
function LoadingSpinner() {
  return (
    <div className="flex items-center justify-center py-8">
      <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600" />
      <span className="ml-3 text-slate-600">Computing...</span>
    </div>
  );
}

/**
 * Single point result display.
 */
function SingleResultDisplay() {
  const m = useBoysStore((state) => state.m);
  const T = useBoysStore((state) => state.T);
  const compute = useBoysStore((state) => state.compute);
  const mode = useBoysStore((state) => state.mode);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 min-h-[200px]">
      <h2 className="text-lg font-semibold text-slate-800 mb-4">Result</h2>

      {/* Idle state */}
      {compute.status === 'idle' && (
        <div className="text-center py-8 text-slate-500">
          <p className="text-lg">Waiting for computation...</p>
          <p className="text-sm mt-1">Adjust parameters to see results</p>
        </div>
      )}

      {/* Computing state */}
      {compute.status === 'computing' && <LoadingSpinner />}

      {/* Error state */}
      {compute.status === 'error' && (
        <div className="bg-red-50 border border-red-200 rounded-lg p-4">
          <p className="text-red-800 font-medium">Computation Error</p>
          <p className="text-red-600 text-sm mt-1">{compute.error}</p>
        </div>
      )}

      {/* Success state */}
      {compute.status === 'success' && (
        <div className="space-y-4">
          {/* Main result */}
          <div className="bg-slate-50 rounded-lg p-4">
            <div className="text-sm text-slate-600 mb-1">
              <MathInline>{`F_{${m}}(${Number.isInteger(T) ? T.toFixed(1) : T}) =`}</MathInline>
            </div>
            <div className="text-2xl font-mono text-slate-900 break-all">
              {formatScientific(compute.result.value)}
            </div>
          </div>

          {/* Method badge */}
          <div className="flex items-center gap-2">
            <span className="text-sm text-slate-600">Method:</span>
            {(() => {
              const style = getMethodStyleWithM(compute.result?.method, compute.result?.m ?? m);
              return (
                <span
                  className={`px-3 py-1 rounded-full text-sm font-medium ${style.bg} ${style.text}`}
                >
                  {style.label}
                </span>
              );
            })()}
          </div>

          {/* Conditional content based on display mode */}
          {mode === 'explain' ? (
            /* Explain mode - educational content */
            <div className="text-sm text-slate-600 bg-slate-50 rounded-lg p-3">
              <span className="font-medium">Why this method? </span>
              {getMethodExplanationContent(compute.result?.method, compute.result?.m ?? m, compute.result?.t ?? T)}
            </div>
          ) : (
            /* Internals mode - computational details */
            <BoysInternalsPanel result={compute.result} />
          )}

          {/* Additional details */}
          <div className="text-xs text-slate-500 border-t border-slate-200 pt-3 mt-3">
            <div className="grid grid-cols-2 gap-2">
              <div>
                <span className="font-medium">Order (<MathInline>m</MathInline>):</span> {compute.result.m}
              </div>
              <div>
                <span className="font-medium">Argument (<MathInline>T</MathInline>):</span>{' '}
                {compute.result.t.toFixed(6)}
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

/**
 * Sweep chart display with regime legend or internals.
 */
function SweepResultDisplay() {
  const m = useBoysStore((state) => state.m);
  const T = useBoysStore((state) => state.T);
  const logScale = useBoysStore((state) => state.logScale);
  const sweepData = useBoysStore((state) => state.sweepData);
  const sweepStatus = useBoysStore((state) => state.sweepStatus);
  const sweepError = useBoysStore((state) => state.sweepError);
  const mode = useBoysStore((state) => state.mode);
  const compute = useBoysStore((state) => state.compute);

  // Use the sweep hook to trigger computations
  const { isComputing } = useBoysSweep();

  // Get the turnover point for the current m
  const turnover = getTurnoverPoint(m);

  return (
    <div className="space-y-4">
      {/* Chart */}
      <BoysChart
        sweepData={sweepData}
        logScale={logScale}
        currentT={T}
        m={m}
        loading={isComputing || sweepStatus === 'pending'}
        error={sweepError}
      />

      {/* Conditional content based on display mode */}
      {mode === 'explain' ? (
        /* Explain mode - 3-regime legend with theoretical explanation */
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
          <h3 className="text-sm font-semibold text-slate-700 mb-3">
            Theoretical Computational Regimes (Lecture Notes, Chapter 4)
          </h3>

          {/* 3-regime legend */}
          <div className="grid grid-cols-2 md:grid-cols-4 gap-3 text-xs mb-4">
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-blue-100 border border-blue-300" />
              <span className="text-slate-600">
                <span className="font-medium">Small <MathInline>T</MathInline></span>
                <span className="block text-slate-500">(<MathInline>{`T < ${THEORETICAL_SMALL_T_THRESHOLD}`}</MathInline>)</span>
              </span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-green-100 border border-green-300" />
              <span className="text-slate-600">
                <span className="font-medium">Moderate <MathInline>T</MathInline></span>
                <span className="block text-slate-500">(<MathInline>{`${THEORETICAL_SMALL_T_THRESHOLD} \\leq T < ${getAsymptoticThreshold(m)}`}</MathInline>)</span>
              </span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-violet-100 border border-violet-300" />
              <span className="text-slate-600">
                <span className="font-medium">Large <MathInline>T</MathInline></span>
                <span className="block text-slate-500">(<MathInline>{`T \\geq ${getAsymptoticThreshold(m)}`}</MathInline>)</span>
              </span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-4 h-2 border-t-2 border-red-500" />
              <span className="text-slate-600">
                <span className="font-medium">Current <MathInline>T</MathInline></span>
              </span>
            </div>
          </div>

          {/* Theoretical formulas from lecture notes */}
          <div className="bg-slate-50 rounded-lg p-3 mb-3 space-y-3">
            <h4 className="text-xs font-semibold text-slate-600 mb-2">Theoretical Formulas (Lecture Notes)</h4>

            <div className="space-y-2">
              <div>
                <span className="text-xs text-blue-600 font-medium">Small <MathInline>T</MathInline> (Series Expansion):</span>
                <div className="bg-white px-2 py-1 rounded mt-1 overflow-x-auto">
                  <MathInline className="text-xs">{THEORETICAL_FORMULAS.series}</MathInline>
                </div>
              </div>

              <div>
                <span className="text-xs text-green-600 font-medium">Moderate <MathInline>T</MathInline> (erf + Upward Recurrence):</span>
                <div className="bg-white px-2 py-1 rounded mt-1 overflow-x-auto">
                  <MathInline className="text-xs">{THEORETICAL_FORMULAS.recurrence}</MathInline>
                </div>
              </div>

              <div>
                <span className="text-xs text-violet-600 font-medium">Large <MathInline>T</MathInline> (Asymptotic Expansion):</span>
                <div className="bg-white px-2 py-1 rounded mt-1 overflow-x-auto">
                  <MathInline className="text-xs">{THEORETICAL_FORMULAS.asymptotic}</MathInline>
                </div>
              </div>
            </div>
          </div>

          {/* Implementation note */}
          <div className="bg-amber-50 border border-amber-200 rounded-lg p-3">
            <h4 className="text-xs font-semibold text-amber-700 mb-1">Implementation (qc-core/libcint)</h4>
            <p className="text-xs text-amber-700">
              The implementation uses 2 methods with an <MathInline>m</MathInline>-dependent turnover point:
            </p>

            {/* Method details based on current m */}
            {m === 0 && (
              <div className="text-xs text-amber-600 mt-1 bg-amber-100 rounded p-1.5">
                <span className="font-medium">For m=0:</span>{' '}
                <MathInline>{`F_0`}</MathInline> is computed <strong>directly</strong> from erf (no recurrence).
                Turnover is 0 because there is no stability concern.
              </div>
            )}
            {m === 1 && (
              <div className="text-xs text-amber-600 mt-1 bg-amber-100 rounded p-1.5">
                <span className="font-medium">For m=1:</span>{' '}
                Only <strong>1 recurrence step</strong> from <MathInline>{`F_0`}</MathInline>. Turnover is 0
                because one step is always stable.
              </div>
            )}
            {m >= 2 && (
              <ul className="text-xs text-amber-600 mt-1 space-y-0.5">
                <li>
                  <span className="font-medium">Series + Downward:</span>{' '}
                  <MathInline>{`T < ${turnover.toFixed(2)}`}</MathInline> (stable for small T)
                </li>
                <li>
                  <span className="font-medium">erf + Upward:</span>{' '}
                  <MathInline>{`T \\geq ${turnover.toFixed(2)}`}</MathInline> (stable when T large enough)
                </li>
              </ul>
            )}

            {turnover > 0 && turnover < getMaxTValue(m) && (
              <p className="text-xs text-amber-600 mt-1">
                The solid gray line at <MathInline>{`T=${turnover.toFixed(2)}`}</MathInline> shows this boundary.
              </p>
            )}

            <p className="text-xs text-amber-600 mt-1 italic">
              No separate asymptotic method - the recurrence formula is accurate for all large <MathInline>T</MathInline>.
            </p>
          </div>

          <p className="text-xs text-slate-500 mt-3">
            The asymptotic threshold (<MathInline>{`30+5m`}</MathInline>) for <MathInline>{`m=${m}`}</MathInline> is <MathInline>{`T=${getAsymptoticThreshold(m)}`}</MathInline>.
            {getAsymptoticThreshold(m) > getMaxTValue(m) && (
              <> This is beyond the chart range (<MathInline>{`T < ${getMaxTValue(m)}`}</MathInline>), so only Small and Moderate regimes are visible.</>
            )}
            {' '}The red line marks the current <MathInline>T</MathInline> value ({T.toFixed(1)}).
          </p>
        </div>
      ) : (
        /* Internals mode - computational details for current T */
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
          <h3 className="text-sm font-semibold text-slate-700 mb-3">
            Computational Details at <MathInline>{`T = ${T.toFixed(1)}`}</MathInline>
          </h3>
          {compute.status === 'success' ? (
            <BoysInternalsPanel result={compute.result} />
          ) : compute.status === 'computing' ? (
            <div className="text-sm text-slate-500">Computing...</div>
          ) : compute.status === 'error' ? (
            <div className="text-sm text-red-600">{compute.error}</div>
          ) : (
            <div className="text-sm text-slate-500">
              Adjust parameters to see computational details.
            </div>
          )}
        </div>
      )}
    </div>
  );
}

/**
 * Multi-order result display (F_0..F_m at fixed T).
 *
 * Shows a bar chart and table comparing Boys function values
 * across all orders from 0 to m at the current T value.
 */
function MultiOrderResultDisplay() {
  const T = useBoysStore((state) => state.T);
  const multiOrderData = useBoysStore((state) => state.multiOrderData);
  const multiOrderStatus = useBoysStore((state) => state.multiOrderStatus);
  const multiOrderError = useBoysStore((state) => state.multiOrderError);

  // Use the multi-order hook to trigger computations
  const { isComputing } = useBoysEvalAll();

  return (
    <BoysMultiOrderPanel
      data={multiOrderData}
      loading={isComputing || multiOrderStatus === 'pending'}
      error={multiOrderError}
      currentT={T}
    />
  );
}

/**
 * Display area for Boys function computation results.
 *
 * Shows:
 * - Single view: Loading state, error state, success with value and method
 * - Sweep view: Interactive chart with regime annotations
 * - Multi-order mode: Bar chart + table of F_0..F_m at fixed T
 *
 * @example
 * ```tsx
 * <BoysResultDisplay />
 * ```
 */
export function BoysResultDisplay() {
  const view = useBoysStore((state) => state.view);
  const mode = useBoysStore((state) => state.mode);

  // Multi-order mode overrides view selection
  if (mode === 'multi-order') {
    return <MultiOrderResultDisplay />;
  }

  if (view === 'sweep') {
    return <SweepResultDisplay />;
  }

  return <SingleResultDisplay />;
}

export default BoysResultDisplay;
