/**
 * BasisComparisonTable - Summary table and discussion prompts.
 *
 * Displays a tabular summary of SCF results across different basis sets
 * with columns for basis name, number of basis functions, energy,
 * iterations, and convergence status. Sorted by energy ascending
 * (lowest energy = most accurate per the variational principle).
 *
 * Below the table, guided discussion prompts target common student
 * misconceptions about basis set quality and the variational principle.
 *
 * @module components/scf/BasisComparisonTable
 * @see US-045 Basis Set Comparison Mode
 */

import { useMemo } from 'react';
import { useScfStore } from '../../stores/scfStore';
import type { BasisComparisonResult } from '../../stores/scfStore';

// ============================================================================
// Constants
// ============================================================================

/**
 * Human-readable basis set label mapping.
 */
const BASIS_LABELS: Record<string, string> = {
  'sto-3g': 'STO-3G',
  '3-21g': '3-21G',
  '6-31g': '6-31G',
  '6-31g*': '6-31G*',
  '6-31+g*': '6-31+G*',
  'cc-pvdz': 'cc-pVDZ',
};

/**
 * Color dots for basis sets matching the chart colors.
 */
const BASIS_COLORS: string[] = [
  '#3b82f6', // blue
  '#10b981', // green
  '#f59e0b', // amber
  '#8b5cf6', // purple
];

// ============================================================================
// Discussion Prompts
// ============================================================================

/**
 * Guided discussion prompts targeting common misconceptions about
 * basis set quality and the variational principle.
 */
const DISCUSSION_PROMPTS = [
  {
    id: 'variational',
    question:
      'Why is the energy with a larger basis set (e.g., 6-31G*) lower than with STO-3G? Does this mean STO-3G gives a "wrong" answer?',
    hint:
      'Hint: The variational principle states that any approximate wavefunction gives an energy that is an upper bound to the true energy. A larger basis set provides more flexibility to lower the energy, but even the STO-3G result is valid -- it is simply a less tight upper bound.',
  },
  {
    id: 'cost',
    question:
      'How does the number of basis functions scale with basis set size? How might this affect computation time for a larger molecule?',
    hint:
      'Hint: The number of two-electron integrals scales as O(N^4) where N is the number of basis functions. Doubling N increases the integral count by roughly 16x. Consider the trade-off between accuracy and computational cost.',
  },
  {
    id: 'practical',
    question:
      'If you needed to compute energies for 1000 molecules quickly, which basis set would you choose and why?',
    hint:
      'Hint: In computational chemistry practice, the choice depends on the property of interest and the required accuracy. For screening or relative comparisons, a smaller basis may suffice. For quantitative thermochemistry, larger basis sets are essential.',
  },
];

// ============================================================================
// Sub-Components
// ============================================================================

/**
 * Single row in the comparison table.
 */
function ComparisonRow({
  result,
  rank,
  colorIndex,
  isLowest,
}: {
  result: BasisComparisonResult;
  rank: number;
  colorIndex: number;
  isLowest: boolean;
}) {
  const label = BASIS_LABELS[result.basisName] ?? result.basisName.toUpperCase();
  const color = BASIS_COLORS[colorIndex % BASIS_COLORS.length];

  return (
    <tr
      className={`border-b border-slate-100 ${
        result.error ? 'bg-red-50' : isLowest ? 'bg-blue-50' : ''
      }`}
    >
      {/* Rank */}
      <td className="py-2 px-3 text-sm text-slate-500">{rank}</td>

      {/* Basis name with color indicator */}
      <td className="py-2 px-3">
        <div className="flex items-center gap-2">
          <span
            className="w-2.5 h-2.5 rounded-full flex-shrink-0"
            style={{ backgroundColor: color }}
            aria-hidden="true"
          />
          <span className="text-sm font-medium text-slate-800">{label}</span>
        </div>
      </td>

      {/* Number of basis functions */}
      <td className="py-2 px-3 text-right font-mono text-sm text-slate-800">
        {result.error ? '--' : result.nbf}
      </td>

      {/* Energy */}
      <td className="py-2 px-3 text-right font-mono text-sm text-slate-800">
        {result.error ? (
          <span className="text-red-600 text-xs">{result.error}</span>
        ) : (
          result.energy.toFixed(10)
        )}
      </td>

      {/* Iterations */}
      <td className="py-2 px-3 text-right font-mono text-sm text-slate-800">
        {result.error ? '--' : result.iterations}
      </td>

      {/* Converged */}
      <td className="py-2 px-3 text-center">
        {result.error ? (
          <span className="text-red-600" aria-label="error">
            &#10007;
          </span>
        ) : result.converged ? (
          <span className="text-green-600" aria-label="converged">
            &#10003;
          </span>
        ) : (
          <span className="text-amber-600" aria-label="unconverged">
            &#8226;
          </span>
        )}
      </td>

      {/* Compute time */}
      <td className="py-2 px-3 text-right text-xs text-slate-500">
        {result.error
          ? '--'
          : result.computeTimeMs < 1000
            ? `${result.computeTimeMs} ms`
            : `${(result.computeTimeMs / 1000).toFixed(1)} s`}
      </td>
    </tr>
  );
}

// ============================================================================
// Main Component
// ============================================================================

/**
 * Summary table of basis set comparison results with discussion prompts.
 *
 * Features:
 * - Sorted by energy ascending (lowest = most accurate)
 * - Color-coded rows matching the convergence chart
 * - Highlighted lowest-energy row (blue background)
 * - Discussion prompts below the table for educational engagement
 *
 * Only renders when comparison results exist.
 *
 * @example
 * ```tsx
 * <BasisComparisonTable />
 * ```
 */
export function BasisComparisonTable() {
  const results = useScfStore((state) => state.basisCompareState.results);
  const comparing = useScfStore((state) => state.basisCompareState.comparing);

  // Sort results by energy ascending (errors at the end)
  const sortedResults = useMemo(() => {
    return [...results].sort((a, b) => {
      if (a.error && b.error) return 0;
      if (a.error) return 1;
      if (b.error) return -1;
      return a.energy - b.energy;
    });
  }, [results]);

  // Build a map from basisName to original index (for color assignment)
  const basisColorMap = useMemo(() => {
    const map = new Map<string, number>();
    results.forEach((r, i) => map.set(r.basisName, i));
    return map;
  }, [results]);

  // Find lowest energy for highlighting
  const lowestEnergy = useMemo(() => {
    const converged = sortedResults.filter((r) => r.converged && !r.error);
    if (converged.length === 0) return null;
    return converged[0].energy;
  }, [sortedResults]);

  // Don't render if no results
  if (results.length === 0 && !comparing) {
    return null;
  }

  return (
    <div className="space-y-4">
      {/* Summary Table */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
        <div className="px-4 py-3 bg-slate-50 border-b border-slate-200">
          <h3 className="text-sm font-semibold text-slate-700">
            Basis Set Comparison Summary
          </h3>
        </div>

        <div className="overflow-x-auto">
          <table className="w-full">
            <thead>
              <tr className="border-b border-slate-200 bg-slate-50/50">
                <th className="text-left py-2.5 px-3 text-xs font-medium text-slate-600 w-10">
                  #
                </th>
                <th className="text-left py-2.5 px-3 text-xs font-medium text-slate-600">
                  Basis Set
                </th>
                <th className="text-right py-2.5 px-3 text-xs font-medium text-slate-600">
                  # Functions
                </th>
                <th className="text-right py-2.5 px-3 text-xs font-medium text-slate-600">
                  Energy (Ha)
                </th>
                <th className="text-right py-2.5 px-3 text-xs font-medium text-slate-600">
                  Iterations
                </th>
                <th className="text-center py-2.5 px-3 text-xs font-medium text-slate-600">
                  Converged
                </th>
                <th className="text-right py-2.5 px-3 text-xs font-medium text-slate-600">
                  Time
                </th>
              </tr>
            </thead>
            <tbody>
              {sortedResults.map((result, idx) => (
                <ComparisonRow
                  key={result.basisName}
                  result={result}
                  rank={idx + 1}
                  colorIndex={basisColorMap.get(result.basisName) ?? idx}
                  isLowest={lowestEnergy !== null && result.energy === lowestEnergy && !result.error}
                />
              ))}
            </tbody>
          </table>
        </div>

        {/* Variational principle note */}
        {sortedResults.length >= 2 && !comparing && (
          <div className="px-4 py-2.5 bg-blue-50 border-t border-blue-100 text-xs text-blue-800">
            <span className="font-medium">Variational principle:</span>{' '}
            The lowest energy is the best upper bound to the true ground-state energy.
            Larger basis sets provide more variational flexibility.
          </div>
        )}
      </div>

      {/* Discussion Prompts */}
      {sortedResults.length >= 2 && !comparing && (
        <div className="bg-amber-50 border border-amber-200 rounded-xl p-4">
          <h4 className="text-sm font-semibold text-amber-900 mb-3">
            Discussion Questions
          </h4>
          <ol className="space-y-3">
            {DISCUSSION_PROMPTS.map((prompt, idx) => (
              <li key={prompt.id} className="text-sm text-amber-900">
                <p className="font-medium mb-0.5">
                  {idx + 1}. {prompt.question}
                </p>
                <p className="text-xs text-amber-700 italic">
                  {prompt.hint}
                </p>
              </li>
            ))}
          </ol>
        </div>
      )}
    </div>
  );
}

export default BasisComparisonTable;
