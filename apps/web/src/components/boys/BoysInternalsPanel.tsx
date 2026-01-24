/**
 * BoysInternalsPanel - Display computational internals for Boys function.
 *
 * Shows detailed information about the computational method, formula,
 * iteration steps, error bounds, and regime selection thresholds.
 * Used in "Internals" display mode for advanced users.
 *
 * @module components/boys/BoysInternalsPanel
 */

import type { BoysEvalResult, BoysMethod } from '../../worker/protocol';
import { Math as MathInline, MathDisplay } from '../common';
import {
  getTurnoverPoint,
  getAsymptoticThreshold,
  getTheoreticalRegime,
  TURNOVER_POINTS,
  THEORETICAL_SMALL_T_THRESHOLD,
  THEORETICAL_FORMULAS,
  IMPLEMENTATION_FORMULAS,
} from '../../lib/boysConstants';

/**
 * Method details for each computational regime (IMPLEMENTATION details).
 */
interface MethodDetail {
  label: string;
  formula: string; // LaTeX formula
  description: React.ReactNode;
}

/**
 * Get method details with m-aware descriptions.
 * For m=0 and m=1, the "recurrence" method has special behavior.
 */
function getMethodDetails(method: BoysMethod, m: number): MethodDetail {
  switch (method) {
    case 'zero':
      return {
        label: 'Zero (Direct)',
        formula: IMPLEMENTATION_FORMULAS.zero,
        description: (
          <>
            <MathInline>T</MathInline> is effectively zero (<MathInline>{`T < 10^{-14}`}</MathInline>), so the exact formula <MathInline>{`F_m(0) = 1/(2m+1)`}</MathInline> is used directly.
          </>
        ),
      };
    case 'series':
      return {
        label: 'Series + Downward Recurrence',
        formula: IMPLEMENTATION_FORMULAS.series,
        description: (
          <>
            <MathInline>T</MathInline> is below the turnover point (<MathInline>{`T < ${getTurnoverPoint(m).toFixed(2)}`}</MathInline>).
            Series expansion computes <MathInline>{`F_m`}</MathInline>, then <strong>downward</strong> recurrence
            calculates lower orders. This is stable for small <MathInline>T</MathInline> because terms reinforce rather than cancel.
          </>
        ),
      };
    case 'recurrence':
      // Special cases for m=0 and m=1
      if (m === 0) {
        return {
          label: 'erf (Direct)',
          formula: 'F_0(T) = \\frac{\\sqrt{\\pi}}{2\\sqrt{T}} \\cdot \\text{erf}(\\sqrt{T})',
          description: (
            <>
              For <MathInline>m=0</MathInline>, <MathInline>{`F_0`}</MathInline> is computed <strong>directly</strong> from
              the error function. No recurrence is needed, so there is no stability concern.
              This is why the turnover is 0 for <MathInline>m=0</MathInline>.
            </>
          ),
        };
      }
      if (m === 1) {
        return {
          label: 'erf + 1 Recurrence Step',
          formula: 'F_0 = \\frac{\\sqrt{\\pi}}{2\\sqrt{T}} \\text{erf}(\\sqrt{T}), \\quad F_1 = \\frac{F_0 - e^{-T}}{2T}',
          description: (
            <>
              For <MathInline>m=1</MathInline>, compute <MathInline>{`F_0`}</MathInline> from erf,
              then just <strong>one</strong> recurrence step to get <MathInline>{`F_1`}</MathInline>.
              A single step is stable even for small <MathInline>T</MathInline>, which is why
              the turnover is 0 for <MathInline>m=1</MathInline>.
            </>
          ),
        };
      }
      // General case: m >= 2 using recurrence
      return {
        label: 'erf + Upward Recurrence',
        formula: IMPLEMENTATION_FORMULAS.recurrence,
        description: (
          <>
            <MathInline>T</MathInline> is at or above the turnover point (<MathInline>{`T \\geq ${getTurnoverPoint(m).toFixed(2)}`}</MathInline>).
            Computed <MathInline>{`F_0`}</MathInline> via error function, then used <strong>upward</strong> recurrence
            for <MathInline>{`${m}`}</MathInline> steps. This is stable because when <MathInline>T</MathInline> is large enough,
            <MathInline>{`(2m+1)F_m \\gg e^{-T}`}</MathInline>, avoiding catastrophic cancellation.
          </>
        ),
      };
    default:
      return {
        label: 'Unknown',
        formula: '',
        description: 'Unknown computational method.',
      };
  }
}


/**
 * Theoretical regime details from lecture notes.
 */
interface TheoreticalRegimeDetail {
  label: string;
  formula: string;
  color: string;
  description: string;
}

const THEORETICAL_REGIME_DETAILS: Record<'small' | 'moderate' | 'large', TheoreticalRegimeDetail> = {
  small: {
    label: 'Small T',
    formula: THEORETICAL_FORMULAS.series,
    color: 'text-blue-600',
    description: 'Series expansion is efficient and convergent for small arguments.',
  },
  moderate: {
    label: 'Moderate T',
    formula: THEORETICAL_FORMULAS.recurrence,
    color: 'text-green-600',
    description: 'Error function formula with upward recurrence provides numerical stability.',
  },
  large: {
    label: 'Large T (Asymptotic)',
    formula: THEORETICAL_FORMULAS.asymptotic,
    color: 'text-violet-600',
    description: 'Asymptotic expansion achieves machine precision for large arguments.',
  },
};

/**
 * Props for BoysInternalsPanel.
 */
interface BoysInternalsPanelProps {
  /** The computation result to display internals for */
  result: BoysEvalResult;
}

/**
 * Formats an error bound value for display.
 *
 * @param error - The error value, null, or undefined
 * @returns Formatted string representation
 */
function formatErrorBound(error: number | null | undefined): string {
  if (error === null || error === undefined) {
    return 'Not estimated (recurrence method)';
  }
  if (error < 1e-15) {
    return '< 1e-15 (machine precision)';
  }
  return `~ ${error.toExponential(2)}`;
}

/**
 * Generates description of computation steps based on method.
 *
 * @param method - The computational method used
 * @param termsCount - Number of terms/steps used
 * @param m - Order of the Boys function
 * @returns Human-readable description of steps
 */
function getStepsDescription(
  method: BoysMethod,
  termsCount: number,
  m: number
): React.ReactNode {
  switch (method) {
    case 'zero':
      return 'Direct calculation (no iteration needed)';
    case 'series':
      return (
        <>
          {termsCount} series terms until convergence, then {m > 0 ? `${m} downward recurrence steps` : 'no recurrence needed'}
        </>
      );
    case 'recurrence':
      if (m === 0) {
        return (
          <>
            Direct erf formula (no recurrence needed for <MathInline>m=0</MathInline>)
          </>
        );
      }
      if (m === 1) {
        return (
          <>
            1 upward recurrence step (<MathInline>{`F_0 \\to F_1`}</MathInline>)
          </>
        );
      }
      return (
        <>
          {m} upward recurrence steps (<MathInline>{`F_0 \\to F_{${m}}`}</MathInline>)
        </>
      );
  }
}

/**
 * Panel displaying computational internals for Boys function evaluation.
 *
 * Shows:
 * - Computational regime and method details
 * - Mathematical formula used
 * - Number of computation steps
 * - Estimated error bound
 * - Regime selection thresholds
 *
 * @example
 * ```tsx
 * <BoysInternalsPanel result={boysResult} />
 * ```
 */
export function BoysInternalsPanel({ result }: BoysInternalsPanelProps) {
  // Defensive check: ensure result is valid before accessing properties
  if (!result || !result.method) {
    return (
      <div className="bg-slate-50 rounded-lg p-4 border border-slate-200">
        <p className="text-sm text-slate-500">No computation result available.</p>
      </div>
    );
  }

  // Get m-aware method details (special handling for m=0 and m=1)
  const methodDetails = getMethodDetails(result.method, result.m);

  // Get the theoretical regime for the current (m, T)
  const theoreticalRegime = getTheoreticalRegime(result.m, result.t);
  const theoreticalDetails = THEORETICAL_REGIME_DETAILS[theoreticalRegime];
  const asymptoticThreshold = getAsymptoticThreshold(result.m);
  const turnover = getTurnoverPoint(result.m);

  return (
    <div className="bg-slate-50 rounded-lg p-4 border border-slate-200 space-y-4">
      {/* Theoretical Regime (from lecture notes) */}
      <div className="bg-white rounded-lg p-3 border border-slate-200">
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-2">
          Theoretical Regime (Lecture Notes)
        </h4>
        <div className="flex items-center gap-2 mb-2">
          <span className={`px-2 py-1 rounded text-sm font-medium ${
            theoreticalRegime === 'small' ? 'bg-blue-100 text-blue-800' :
            theoreticalRegime === 'moderate' ? 'bg-green-100 text-green-800' :
            'bg-violet-100 text-violet-800'
          }`}>
            {theoreticalDetails.label}
          </span>
          <span className="text-xs text-slate-500">
            {theoreticalRegime === 'small' && <MathInline>{`T < ${THEORETICAL_SMALL_T_THRESHOLD}`}</MathInline>}
            {theoreticalRegime === 'moderate' && <MathInline>{`${THEORETICAL_SMALL_T_THRESHOLD} \\leq T < ${asymptoticThreshold}`}</MathInline>}
            {theoreticalRegime === 'large' && <MathInline>{`T \\geq ${asymptoticThreshold}`}</MathInline>}
          </span>
        </div>
        <div className="bg-slate-50 px-3 py-2 rounded overflow-x-auto mb-2">
          <MathDisplay className="!my-0 text-sm">{theoreticalDetails.formula}</MathDisplay>
        </div>
        <p className="text-xs text-slate-600">{theoreticalDetails.description}</p>
      </div>

      {/* Implementation Method (what qc-core actually uses) */}
      <div>
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-2">
          Implementation Method (qc-core)
        </h4>
        <div className="flex items-center gap-2 mb-2">
          <span className={`px-2 py-1 rounded text-sm font-medium ${
            result.method === 'zero' ? 'bg-amber-100 text-amber-800' :
            result.method === 'series' ? 'bg-blue-100 text-blue-800' :
            'bg-green-100 text-green-800'
          }`}>
            {methodDetails.label}
          </span>
        </div>
        <div className="bg-slate-100 px-3 py-2 rounded overflow-x-auto mb-2">
          <MathDisplay className="!my-0 text-sm">{methodDetails.formula}</MathDisplay>
        </div>
        <p className="text-xs text-slate-600">{methodDetails.description}</p>
      </div>

      {/* Computation Steps */}
      <div>
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-1">
          Computation Steps
        </h4>
        <p className="text-sm text-slate-700">
          {getStepsDescription(result.method, result.termsCount, result.m)}
        </p>
      </div>

      {/* Estimated Error Bound */}
      <div>
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-1">
          Estimated Error Bound
        </h4>
        <p className="text-sm font-mono text-slate-700">
          {formatErrorBound(result.estimatedError)}
        </p>
        {result.method === 'series' && result.estimatedError !== null && (
          <p className="text-xs text-slate-500 mt-1">
            Relative error based on series convergence tolerance
          </p>
        )}
      </div>

      {/* Regime Thresholds Comparison */}
      <div className="pt-3 border-t border-slate-200">
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-3">
          Regime Thresholds Comparison
        </h4>

        {/* Theoretical Regimes Table */}
        <div className="mb-4">
          <p className="text-xs font-medium text-slate-700 mb-2">Theoretical (Lecture Notes, Chapter 4):</p>
          <div className="text-xs text-slate-600 space-y-1 pl-2">
            <div className="flex justify-between">
              <span className="text-blue-600">Small <MathInline>T</MathInline> (Series):</span>
              <MathInline>{`T < ${THEORETICAL_SMALL_T_THRESHOLD}`}</MathInline>
            </div>
            <div className="flex justify-between">
              <span className="text-green-600">Moderate <MathInline>T</MathInline> (erf + Recurrence):</span>
              <MathInline>{`${THEORETICAL_SMALL_T_THRESHOLD} \\leq T < ${asymptoticThreshold}`}</MathInline>
            </div>
            <div className="flex justify-between">
              <span className="text-violet-600">Large <MathInline>T</MathInline> (Asymptotic):</span>
              <MathInline>{`T \\geq ${asymptoticThreshold}`}</MathInline>
            </div>
          </div>
        </div>

        {/* Implementation Regimes Table */}
        <div className="mb-3">
          <p className="text-xs font-medium text-slate-700 mb-2">Implementation (libcint/qc-core):</p>
          <div className="text-xs text-slate-600 space-y-1 pl-2">
            <div className="flex justify-between">
              <span className="text-amber-600">Zero:</span>
              <MathInline>{`T < 10^{-14}`}</MathInline>
            </div>
            <div className="flex justify-between">
              <span className="text-blue-600">Series:</span>
              {getTurnoverPoint(result.m) > 0
                ? <MathInline>{`T < ${getTurnoverPoint(result.m).toFixed(2)}`}</MathInline>
                : <span className="text-slate-500">(not used for <MathInline>{`m=0,1`}</MathInline>)</span>}
            </div>
            <div className="flex justify-between">
              <span className="text-green-600">Recurrence:</span>
              {getTurnoverPoint(result.m) > 0
                ? <MathInline>{`T \\geq ${getTurnoverPoint(result.m).toFixed(2)}`}</MathInline>
                : <span>all <MathInline>{`T > 0`}</MathInline></span>}
            </div>
          </div>
        </div>

        {/* Explanation Box */}
        <div className="bg-amber-50 border border-amber-200 rounded p-2">
          <p className="text-xs text-amber-700">
            <span className="font-medium">Key Difference:</span> The implementation uses <MathInline>m</MathInline>-dependent
            turnover points (much smaller than 25) optimized for double precision. The "Recurrence" method handles
            both moderate and large <MathInline>T</MathInline> theoretical regimes. There is no separate asymptotic method
            in the implementation because the recurrence formula remains stable and accurate for all large <MathInline>T</MathInline>.
          </p>
        </div>
      </div>

      {/* Current Values Summary */}
      <div className="pt-2 border-t border-slate-200">
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-2">
          Current Parameters
        </h4>
        <div className="text-xs text-slate-600 space-y-1">
          <p>
            <MathInline>{`m = ${result.m}`}</MathInline>,{' '}
            <MathInline>{`T = ${result.t.toFixed(4)}`}</MathInline>
          </p>
          <p>
            Implementation turnover for <MathInline>{`m=${result.m}`}</MathInline>:{' '}
            {turnover > 0 ? (
              <MathInline>{`T_0 = ${turnover.toFixed(4)}`}</MathInline>
            ) : (
              <span className="text-amber-600 font-medium">0 (always recurrence)</span>
            )}
          </p>
          <p>
            Theoretical asymptotic threshold (<MathInline>{`30+5m`}</MathInline>):{' '}
            <MathInline>{`T = ${asymptoticThreshold}`}</MathInline>
          </p>
        </div>

        {/* Special explanation for m=0 and m=1 */}
        {result.m === 0 && (
          <div className="bg-blue-50 border border-blue-200 rounded p-2 mt-2">
            <p className="text-xs text-blue-700">
              <span className="font-medium">Why m=0 always uses recurrence:</span>{' '}
              <MathInline>{`F_0`}</MathInline> is computed directly from the error function formula.
              No recurrence steps are needed, so there is no numerical stability concern.
            </p>
          </div>
        )}
        {result.m === 1 && (
          <div className="bg-blue-50 border border-blue-200 rounded p-2 mt-2">
            <p className="text-xs text-blue-700">
              <span className="font-medium">Why m=1 always uses recurrence:</span>{' '}
              Only one recurrence step is needed (<MathInline>{`F_0 \\to F_1`}</MathInline>).
              A single step is stable even for relatively small <MathInline>T</MathInline>.
            </p>
          </div>
        )}
        {result.m >= 2 && (
          <div className="bg-slate-100 rounded p-2 mt-2">
            <p className="text-xs text-slate-600">
              For <MathInline>{`m=${result.m}`}</MathInline>:{' '}
              {result.t < turnover ? (
                <>
                  <MathInline>{`T=${result.t.toFixed(2)} < ${turnover.toFixed(2)}`}</MathInline> so <strong>Series + downward recurrence</strong> is used.
                </>
              ) : (
                <>
                  <MathInline>{`T=${result.t.toFixed(2)} \\geq ${turnover.toFixed(2)}`}</MathInline> so <strong>erf + upward recurrence</strong> is used.
                </>
              )}
            </p>
          </div>
        )}

        <p className="text-xs text-slate-400 mt-2 italic">
          Turnover increases with <MathInline>m</MathInline>:{' '}
          <MathInline>{`m=2`}</MathInline> at {TURNOVER_POINTS[2].toFixed(2)},{' '}
          <MathInline>{`m=5`}</MathInline> at {TURNOVER_POINTS[5].toFixed(2)},{' '}
          <MathInline>{`m=10`}</MathInline> at {TURNOVER_POINTS[10].toFixed(2)}.
        </p>
      </div>

      {/* Numerical Stability Note */}
      <div className="pt-2 border-t border-slate-200">
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-2">
          Numerical Stability
        </h4>
        <div className="text-xs text-slate-600 space-y-2">
          <div>
            <span className="font-medium text-green-600">Upward recurrence:</span>{' '}
            <MathInline>{`F_{m+1} = [(2m+1)F_m - e^{-T}]/(2T)`}</MathInline>
            <p className="text-slate-500 mt-0.5">
              Division by <MathInline>2T</MathInline> is problematic for small <MathInline>T</MathInline> (subtraction causes cancellation).
            </p>
          </div>
          <div>
            <span className="font-medium text-blue-600">Downward recurrence:</span>{' '}
            <MathInline>{`F_{m-1} = [2T \\cdot F_m + e^{-T}]/(2m-1)`}</MathInline>
            <p className="text-slate-500 mt-0.5">
              Stable for small <MathInline>T</MathInline> because terms add rather than subtract.
            </p>
          </div>
        </div>
      </div>
    </div>
  );
}

export default BoysInternalsPanel;
