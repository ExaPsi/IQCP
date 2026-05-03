/**
 * EriQuartetDetail - Full breakdown of a contracted two-electron
 * repulsion integral (ERI) into its primitive-quartet contributions.
 *
 * Displays:
 * 1. Header: "(ij|kl) = value" with basis function labels
 * 2. Method badge: Boys Function (blue) or Rys Quadrature (green)
 * 3. T parameter and cross-reference to Module D (for Rys method)
 * 4. Primitive-quartet table sorted by |weightedContribution|
 * 5. Sum verification badge
 *
 * The primitive table mirrors PrimitiveBreakdownTable but with 4 indices
 * instead of 2. It defaults to top-10 with a "Show all" toggle.
 *
 * @module components/integrals/EriQuartetDetail
 * @see US-059 ERI Browser
 * @see Szabo & Ostlund, Eq. 3.155
 * @see Dupuis, Rys & King (1976)
 */

import { useMemo, memo } from 'react';
import type { EriDetailResult } from '../../hooks/useEriDetail';
import { encodeRunState } from '../../lib/deeplink';
import { ROUTES } from '../../lib/routes';
import { APP_VERSION } from '../../types/run-state';

// ============================================================================
// Types
// ============================================================================

export interface EriQuartetDetailProps {
  /** ERI breakdown result data */
  data: EriDetailResult;
  /** Whether to show all primitives or top-N */
  showAll: boolean;
  /** Toggle show all */
  onToggleShowAll: () => void;
}

// ============================================================================
// Constants
// ============================================================================

/** Default number of top contributions to show */
const DEFAULT_TOP_N = 10;

/** Threshold for sum verification (matches Rust invariant) */
const SUM_VERIFICATION_TOLERANCE = 1e-12;

// ============================================================================
// Helpers
// ============================================================================

function formatFixed(value: number, decimals: number): string {
  return value.toFixed(decimals);
}

function formatSci(value: number, decimals: number): string {
  return value.toExponential(decimals);
}

/**
 * Build a deep link URL to Module D (Rys Quadrature Lab)
 * with the given n and T parameters pre-set.
 */
function buildRysDeepLink(n: number, T: number): string {
  const state = {
    schema_version: 'run_state_v1' as const,
    app_version: APP_VERSION,
    module: 'rys' as const,
    rys: {
      n,
      T,
      target: '1e-8' as const,
    },
    ui: { mode: 'explain' as const },
  };
  const encoded = encodeRunState(state);
  return `${ROUTES.RYS}?run=${encoded}`;
}

// ============================================================================
// Sub-components
// ============================================================================

/**
 * Percentage bar for the primitive contribution table.
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

/**
 * Method badge indicating Boys Function or Rys Quadrature.
 */
const MethodBadge = memo(function MethodBadge({
  method,
}: {
  method: EriDetailResult['method'];
}) {
  if (method.type === 'boysFunction') {
    return (
      <span className="inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs font-medium bg-blue-100 text-blue-800 border border-blue-200">
        <span className="w-2 h-2 rounded-full bg-blue-500" aria-hidden="true" />
        Boys Function F_0(T)
      </span>
    );
  }

  return (
    <span className="inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs font-medium bg-emerald-100 text-emerald-800 border border-emerald-200">
      <span className="w-2 h-2 rounded-full bg-emerald-500" aria-hidden="true" />
      Rys Quadrature (n={method.nroots})
    </span>
  );
});

/**
 * Rys quadrature info section with roots/weights table
 * and "Open in Module B" cross-reference link.
 */
const RysInfoSection = memo(function RysInfoSection({
  method,
}: {
  method: Extract<EriDetailResult['method'], { type: 'rysQuadrature' }>;
}) {
  const deepLink = buildRysDeepLink(method.nroots, method.tParameter);

  return (
    <div className="bg-emerald-50 rounded-lg border border-emerald-200 p-3 space-y-3">
      <div className="flex items-center justify-between">
        <h5 className="text-xs font-semibold text-emerald-800 uppercase tracking-wider">
          Rys Quadrature Parameters
        </h5>
        <span className="text-xs font-mono text-emerald-700">
          T = {formatFixed(method.tParameter, 6)}
        </span>
      </div>

      {/* Roots and weights table */}
      <table
        className="w-full text-xs border-collapse"
        aria-label="Rys quadrature roots and weights"
      >
        <thead>
          <tr className="border-b border-emerald-200">
            <th className="px-2 py-1 text-left text-[10px] font-semibold text-emerald-600 uppercase">
              Root #
            </th>
            <th className="px-2 py-1 text-right text-[10px] font-semibold text-emerald-600 uppercase">
              t_i
            </th>
            <th className="px-2 py-1 text-right text-[10px] font-semibold text-emerald-600 uppercase">
              w_i
            </th>
          </tr>
        </thead>
        <tbody>
          {method.roots.map((root, idx) => (
            <tr key={idx} className="border-b border-emerald-100">
              <td className="px-2 py-1 font-mono text-emerald-700">{idx + 1}</td>
              <td className="px-2 py-1 font-mono text-right text-emerald-800">
                {formatFixed(root, 10)}
              </td>
              <td className="px-2 py-1 font-mono text-right text-emerald-800">
                {formatFixed(method.weights[idx], 10)}
              </td>
            </tr>
          ))}
        </tbody>
      </table>

      {/* Module B cross-reference */}
      <div className="flex items-center gap-2 pt-1">
        <a
          href={deepLink}
          className="text-xs text-emerald-700 hover:text-emerald-900 font-medium underline decoration-emerald-400 underline-offset-2 transition-colors"
          title={`Open Rys Quadrature Lab with n=${method.nroots}, T=${formatFixed(method.tParameter, 4)}`}
        >
          Explore in Module D (Rys Quadrature Lab) &rarr;
        </a>
      </div>

      <p className="text-[10px] text-emerald-600 leading-relaxed">
        These roots and weights come from the representative primitive quartet
        (most diffuse, smallest T). Each primitive quartet has its own T value;
        see the T column in the breakdown table below.
      </p>
    </div>
  );
});

/**
 * Boys function info section with T parameter display.
 */
const BoysInfoSection = memo(function BoysInfoSection({
  method,
}: {
  method: Extract<EriDetailResult['method'], { type: 'boysFunction' }>;
}) {
  return (
    <div className="bg-blue-50 rounded-lg border border-blue-200 p-3">
      <div className="flex items-center justify-between mb-2">
        <h5 className="text-xs font-semibold text-blue-800 uppercase tracking-wider">
          Boys Function Parameters
        </h5>
        <span className="text-xs font-mono text-blue-700">
          T = {formatFixed(method.tParameter, 6)}
        </span>
      </div>

      <p className="text-[10px] text-blue-700 leading-relaxed">
        This integral involves only s-type shells (total angular momentum L = 0),
        so it is computed via a single quadrature point equivalent to the Boys
        function F_0(T). For multi-root Rys quadrature, select a quartet involving
        p-type or higher angular momentum functions.
      </p>
    </div>
  );
});

// ============================================================================
// Main Component
// ============================================================================

/**
 * Full breakdown panel for a selected ERI quartet.
 *
 * Features:
 * - Header with (ij|kl) = value and basis function labels
 * - Method badge (Boys / Rys)
 * - Rys info section with roots/weights table and Module D link
 * - Primitive-quartet table sorted by contribution magnitude
 * - Top-N with "Show all" toggle
 * - Sum verification badge
 */
export const EriQuartetDetail = memo(function EriQuartetDetail({
  data,
  showAll,
  onToggleShowAll,
}: EriQuartetDetailProps) {
  const {
    contractedValue,
    indices,
    labels,
    method,
    contributions,
    nPrimitives,
    totalAngularMomentum,
    nroots,
  } = data;

  const totalQuartets = contributions.length;
  const hasMore = totalQuartets > DEFAULT_TOP_N;

  // Visible rows based on showAll toggle
  const visibleContributions = useMemo(() => {
    if (showAll || totalQuartets <= DEFAULT_TOP_N) {
      return contributions;
    }
    return contributions.slice(0, DEFAULT_TOP_N);
  }, [contributions, showAll, totalQuartets]);

  // Sum verification and percentages
  const { sumOfContributions, sumDiff, isVerified, percentages } = useMemo(() => {
    const sum = contributions.reduce(
      (acc, pc) => acc + pc.weightedContribution,
      0,
    );
    const diff = Math.abs(sum - contractedValue);
    const verified = diff < SUM_VERIFICATION_TOLERANCE;

    const absContracted = Math.abs(contractedValue);
    const pcts = contributions.map((pc) =>
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
  }, [contributions, contractedValue]);

  const visiblePercentages = useMemo(() => {
    if (showAll || totalQuartets <= DEFAULT_TOP_N) return percentages;
    return percentages.slice(0, DEFAULT_TOP_N);
  }, [percentages, showAll, totalQuartets]);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5 space-y-4">
      {/* Header */}
      <div>
        <div className="flex items-start justify-between gap-3">
          <div>
            <h3 className="text-sm font-semibold text-slate-800">
              ERI Detail: ({indices[0]}{indices[1]}|{indices[2]}{indices[3]})
            </h3>
            <p className="text-xs text-slate-500 mt-0.5">
              ⟨{labels[0]} {labels[1]} | {labels[2]} {labels[3]}⟩
            </p>
          </div>
          <MethodBadge method={method} />
        </div>

        {/* Value display */}
        <div className="mt-2 px-3 py-2 bg-slate-50 rounded-lg border border-slate-100">
          <span className="text-xs text-slate-500">Value: </span>
          <span className="text-sm font-mono font-semibold text-slate-800">
            {formatSci(contractedValue, 10)}
          </span>
          <span className="text-[10px] text-slate-400 ml-2">Ha</span>
        </div>

        {/* Metadata */}
        <div className="flex flex-wrap gap-x-4 gap-y-1 mt-2 text-[10px] text-slate-500">
          <span>L_total = {totalAngularMomentum}</span>
          <span>nroots = {nroots}</span>
          <span>
            Primitives: {nPrimitives[0]}x{nPrimitives[1]}x{nPrimitives[2]}x{nPrimitives[3]} = {totalQuartets}
          </span>
        </div>
      </div>

      {/* Method info section */}
      {method.type === 'rysQuadrature' ? (
        <RysInfoSection method={method} />
      ) : (
        <BoysInfoSection method={method} />
      )}

      {/* Primitive breakdown table */}
      <div>
        <h4 className="text-sm font-semibold text-slate-700 mb-2">
          Primitive-Quartet Breakdown
        </h4>

        <div className="overflow-x-auto">
          <table
            className="w-full text-xs border-collapse"
            aria-label={`Primitive quartet contributions to ERI (${indices[0]}${indices[1]}|${indices[2]}${indices[3]})`}
          >
            <thead>
              <tr className="border-b border-slate-200">
                <th className="px-1.5 py-1.5 text-left text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                  #
                </th>
                <th className="px-1.5 py-1.5 text-left text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                  (p,q,r,s)
                </th>
                <th className="px-1.5 py-1.5 text-right text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                  Prim. ERI
                </th>
                <th className="px-1.5 py-1.5 text-right text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                  Contribution
                </th>
                <th className="px-1.5 py-1.5 text-right text-[10px] font-semibold text-slate-500 uppercase tracking-wider">
                  T
                </th>
                <th className="px-1.5 py-1.5 text-left text-[10px] font-semibold text-slate-500 uppercase tracking-wider w-24">
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
                    key={`${pc.primIndices[0]}-${pc.primIndices[1]}-${pc.primIndices[2]}-${pc.primIndices[3]}`}
                    className={`border-b border-slate-50 ${
                      isEven ? 'bg-white' : 'bg-slate-50/50'
                    } hover:bg-blue-50/40 transition-colors`}
                  >
                    <td className="px-1.5 py-1 font-mono text-slate-400 text-right">
                      {rank}
                    </td>
                    <td className="px-1.5 py-1 font-mono text-slate-700">
                      ({pc.primIndices[0]},{pc.primIndices[1]},{pc.primIndices[2]},{pc.primIndices[3]})
                    </td>
                    <td className="px-1.5 py-1 font-mono text-right text-slate-700">
                      {formatSci(pc.primitiveValue, 4)}
                    </td>
                    <td className="px-1.5 py-1 font-mono text-right text-slate-700">
                      {formatSci(pc.weightedContribution, 4)}
                    </td>
                    <td className="px-1.5 py-1 font-mono text-right text-slate-500">
                      {formatFixed(pc.tParameter, 4)}
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
              ? `Show top ${DEFAULT_TOP_N} only`
              : `Show all ${totalQuartets} quartets`}
          </button>
        )}
      </div>

      {/* Sum verification */}
      <div
        className={`px-3 py-2 rounded-lg border text-xs font-mono ${
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
      <div className="bg-slate-50 rounded-lg px-3 py-2 border border-slate-100">
        <p className="text-[10px] text-slate-500 leading-relaxed">
          Each contracted ERI is a sum over primitive Gaussian quartets.
          For STO-3G (3 primitives per shell), an (ss|ss) integral has
          3^4 = 81 primitive quartets. The quartets are sorted by
          contribution magnitude: the top rows dominate the sum while
          quartets involving tight primitives far apart contribute
          negligibly -- foreshadowing Schwarz screening in production codes.
        </p>
      </div>
    </div>
  );
});
