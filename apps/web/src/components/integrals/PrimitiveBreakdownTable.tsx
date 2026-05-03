/**
 * PrimitiveBreakdownTable - Table showing primitive-pair contributions
 * to a contracted integral.
 *
 * Displays all primitive Gaussian pairs sorted by magnitude, with
 * expandable "show all" when more than topN pairs exist. Includes
 * a sum verification row confirming the invariant:
 *   sum(weighted_contribution) == contracted_value within 1e-12.
 *
 * @module components/integrals/PrimitiveBreakdownTable
 * @see US-056 Primitive Breakdown Panel
 * @see FR-INT-03
 */

import { useMemo, memo } from 'react';

// ============================================================================
// Types
// ============================================================================

/**
 * A single primitive-pair contribution to a contracted integral.
 *
 * This type is defined locally to decouple from protocol.ts, which
 * may be updated concurrently by another agent. It matches the
 * IntegralBreakdownResult/PrimitiveContributionResult from the
 * worker protocol.
 */
export interface PrimitiveContributionResult {
  /** Primitive indices within shells [p, q] */
  primIndices: [number, number];
  /** Primitive exponents [alpha_p, alpha_q] */
  exponents: [number, number];
  /** Raw contraction coefficients [c_p, c_q] */
  coefficients: [number, number];
  /** Coefficients with normalization [c_p * N_p, c_q * N_q] */
  normCoefficients: [number, number];
  /** Bare primitive integral value */
  primitiveValue: number;
  /** Weighted contribution to contracted integral */
  weightedContribution: number;
}

export interface IntegralBreakdownResult {
  /** Contracted integral value */
  contractedValue: number;
  /** Integral type */
  integralType: string;
  /** Basis function indices [i, j] */
  indices: [number, number];
  /** Basis function labels */
  labels: [string, string];
  /** Number of primitives in shell i */
  nPrimI: number;
  /** Number of primitives in shell j */
  nPrimJ: number;
  /** Primitive contributions sorted by |weightedContribution| descending */
  primitiveContributions: PrimitiveContributionResult[];
}

export interface PrimitiveBreakdownTableProps {
  /** Breakdown data from integral_with_breakdown */
  breakdown: IntegralBreakdownResult;
  /** Whether to show all contributions (vs. top-N) */
  showAll: boolean;
  /** Callback to toggle show all */
  onToggleShowAll: () => void;
  /** Number of contributions to show by default */
  topN?: number;
}

// ============================================================================
// Constants
// ============================================================================

/** Default number of top contributions to show */
const DEFAULT_TOP_N = 5;

/** Threshold for sum verification (should match Rust invariant) */
const SUM_VERIFICATION_TOLERANCE = 1e-12;

// ============================================================================
// Helpers
// ============================================================================

/**
 * Format a floating-point number with a fixed number of decimal places.
 */
function formatFixed(value: number, decimals: number): string {
  return value.toFixed(decimals);
}

/**
 * Format a number in scientific notation with a specified number of
 * significant digits after the decimal point.
 */
function formatSci(value: number, decimals: number): string {
  return value.toExponential(decimals);
}

/**
 * Map the WASM integral type string to a display symbol.
 *
 * The WASM side uses "overlap", "kinetic", "nuclear", "hcore"
 * while the UI display uses S, T, V, H^core.
 */
function integralTypeToSymbol(integralType: string): string {
  switch (integralType.toLowerCase()) {
    case 'overlap':
    case 's':
      return 'S';
    case 'kinetic':
    case 't':
      return 'T';
    case 'nuclear':
    case 'v':
      return 'V';
    case 'hcore':
      return 'H^core';
    default:
      return integralType;
  }
}

// ============================================================================
// Sub-components
// ============================================================================

/**
 * Horizontal percentage bar rendered inline in the table cell.
 * Width is proportional to the percentage, color indicates sign.
 */
const PercentageBar = memo(function PercentageBar({
  percentage,
  isPositive,
}: {
  percentage: number;
  isPositive: boolean;
}) {
  const clampedPct = Math.min(Math.abs(percentage), 100);
  return (
    <div className="flex items-center gap-1.5">
      <div className="flex-1 h-3 bg-slate-100 rounded-sm overflow-hidden">
        <div
          className={`h-full rounded-sm transition-all ${
            isPositive ? 'bg-blue-400' : 'bg-amber-400'
          }`}
          style={{ width: `${clampedPct}%` }}
        />
      </div>
      <span className="text-[10px] font-mono text-slate-600 w-12 text-right">
        {formatFixed(percentage, 1)}%
      </span>
    </div>
  );
});

// ============================================================================
// Main Component
// ============================================================================

/**
 * Table showing primitive-pair contributions to a contracted integral.
 *
 * Features:
 * - Sorted by |weighted_contribution| descending (pre-sorted from WASM)
 * - Top-N with "Show all N pairs" expand button
 * - Sum verification row with green check / red warning
 * - Percentage bars for visual dominance comparison
 * - Zebra-striped rows with monospace numbers
 */
export const PrimitiveBreakdownTable = memo(function PrimitiveBreakdownTable({
  breakdown,
  showAll,
  onToggleShowAll,
  topN = DEFAULT_TOP_N,
}: PrimitiveBreakdownTableProps) {
  const {
    contractedValue,
    integralType,
    indices,
    labels,
    nPrimI,
    nPrimJ,
    primitiveContributions,
  } = breakdown;

  const totalPairs = primitiveContributions.length;
  const hasMore = totalPairs > topN;

  // Visible rows based on showAll toggle
  const visibleContributions = useMemo(() => {
    if (showAll || totalPairs <= topN) {
      return primitiveContributions;
    }
    return primitiveContributions.slice(0, topN);
  }, [primitiveContributions, showAll, topN, totalPairs]);

  // Compute sum and percentages
  const { sumOfContributions, sumDiff, isVerified, percentages } = useMemo(() => {
    const sum = primitiveContributions.reduce(
      (acc, pc) => acc + pc.weightedContribution,
      0,
    );
    const diff = Math.abs(sum - contractedValue);
    const verified = diff < SUM_VERIFICATION_TOLERANCE;

    // Compute percentage for each contribution
    const absContracted = Math.abs(contractedValue);
    const pcts = primitiveContributions.map((pc) =>
      absContracted > 1e-30
        ? (100 * Math.abs(pc.weightedContribution)) / absContracted
        : 0,
    );

    return {
      sumOfContributions: sum,
      sumDiff: diff,
      isVerified: verified,
      percentages: pcts,
    };
  }, [primitiveContributions, contractedValue]);

  // Determine the visible percentages (sliced like contributions)
  const visiblePercentages = useMemo(() => {
    if (showAll || totalPairs <= topN) return percentages;
    return percentages.slice(0, topN);
  }, [percentages, showAll, topN, totalPairs]);

  const symbol = integralTypeToSymbol(integralType);

  return (
    <div>
      {/* Header */}
      <h4 className="text-sm font-semibold text-slate-700 mb-1">
        Primitive Breakdown
      </h4>
      <p className="text-xs text-slate-500 mb-3">
        {symbol}
        <sub>
          {indices[0] + 1},{indices[1] + 1}
        </sub>{' '}
        = {formatSci(contractedValue, 10)}{' '}
        <span className="text-slate-400">
          ({nPrimI} x {nPrimJ} = {totalPairs} primitive pairs)
        </span>
      </p>
      <p className="text-[10px] text-slate-400 mb-2">
        {labels[0]} ⟨{symbol === 'H^core' ? 'T+V' : symbol}⟩{' '}
        {labels[1]}
      </p>

      {/* Table */}
      <div className="overflow-x-auto">
        <table
          className="w-full text-xs border-collapse"
          aria-label={`Primitive pair contributions to ${symbol} element (${indices[0] + 1}, ${indices[1] + 1})`}
        >
          <thead>
            <tr className="border-b border-slate-200">
              <th className="px-1.5 py-1.5 text-left text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                #
              </th>
              <th className="px-1.5 py-1.5 text-left text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                Pair
              </th>
              <th className="px-1.5 py-1.5 text-right text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                &alpha;<sub>p</sub>
              </th>
              <th className="px-1.5 py-1.5 text-right text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                &alpha;<sub>q</sub>
              </th>
              <th className="px-1.5 py-1.5 text-right text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                c<sub>p</sub>
              </th>
              <th className="px-1.5 py-1.5 text-right text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                c<sub>q</sub>
              </th>
              <th className="px-1.5 py-1.5 text-right text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                Prim. Integral
              </th>
              <th className="px-1.5 py-1.5 text-right text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                Contribution
              </th>
              <th className="px-1.5 py-1.5 text-left text-[10px] font-semibold text-slate-500 uppercase tracking-wider w-28">
                % of Total
              </th>
            </tr>
          </thead>
          <tbody>
            {visibleContributions.map((pc, idx) => {
              const rank = idx + 1;
              const isEven = idx % 2 === 0;
              const pct = visiblePercentages[idx];
              const isPositive = pc.weightedContribution >= 0;

              return (
                <tr
                  key={`${pc.primIndices[0]}-${pc.primIndices[1]}`}
                  className={`border-b border-slate-50 ${
                    isEven ? 'bg-white' : 'bg-slate-50/50'
                  } hover:bg-blue-50/40 transition-colors`}
                >
                  <td className="px-1.5 py-1 font-mono text-slate-400 text-right">
                    {rank}
                  </td>
                  <td className="px-1.5 py-1 font-mono text-slate-700">
                    ({pc.primIndices[0]},{pc.primIndices[1]})
                  </td>
                  <td className="px-1.5 py-1 font-mono text-right text-slate-700">
                    {formatFixed(pc.exponents[0], 6)}
                  </td>
                  <td className="px-1.5 py-1 font-mono text-right text-slate-700">
                    {formatFixed(pc.exponents[1], 6)}
                  </td>
                  <td className="px-1.5 py-1 font-mono text-right text-slate-700">
                    {formatFixed(pc.coefficients[0], 6)}
                  </td>
                  <td className="px-1.5 py-1 font-mono text-right text-slate-700">
                    {formatFixed(pc.coefficients[1], 6)}
                  </td>
                  <td className="px-1.5 py-1 font-mono text-right text-slate-700">
                    {formatSci(pc.primitiveValue, 6)}
                  </td>
                  <td className="px-1.5 py-1 font-mono text-right text-slate-700">
                    {formatSci(pc.weightedContribution, 6)}
                  </td>
                  <td className="px-1.5 py-1">
                    <PercentageBar percentage={pct} isPositive={isPositive} />
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      {/* Show all / collapse toggle */}
      {hasMore && (
        <button
          onClick={onToggleShowAll}
          className="mt-2 text-xs text-primary-600 hover:text-primary-800 font-medium transition-colors"
          aria-expanded={showAll}
        >
          {showAll
            ? `Show top ${topN} only`
            : `Show all ${totalPairs} pairs`}
        </button>
      )}

      {/* Sum verification row */}
      <div
        className={`mt-3 px-3 py-2 rounded-lg border text-xs font-mono ${
          isVerified
            ? 'bg-green-50 border-green-200 text-green-800'
            : 'bg-red-50 border-red-200 text-red-800'
        }`}
        role="status"
        aria-label={`Sum verification: ${isVerified ? 'passed' : 'failed'}`}
      >
        <div className="flex items-center gap-2">
          <span
            className={`flex-shrink-0 w-4 h-4 rounded-full flex items-center justify-center text-[10px] font-bold ${
              isVerified
                ? 'bg-green-200 text-green-700'
                : 'bg-red-200 text-red-700'
            }`}
            aria-hidden="true"
          >
            {isVerified ? '\u2713' : '!'}
          </span>
          <div className="flex-1">
            <p>
              Sum = {formatSci(sumOfContributions, 12)}
            </p>
            <p className="text-[10px] opacity-75 mt-0.5">
              Contracted = {formatSci(contractedValue, 12)} | diff ={' '}
              {formatSci(sumDiff, 2)}
            </p>
          </div>
        </div>
      </div>

      {/* Educational note */}
      <div className="mt-3 bg-slate-50 rounded-lg px-3 py-2 border border-slate-100">
        <p className="text-[10px] text-slate-500 leading-relaxed">
          Each contracted integral is a sum over primitive Gaussian pairs.
          The pairs are sorted by magnitude: the top rows contribute most
          to the final value. Diffuse primitives (small exponents) tend to
          dominate overlap integrals at moderate bond distances, while tight
          primitives (large exponents) contribute more to kinetic energy.
        </p>
      </div>
    </div>
  );
});
