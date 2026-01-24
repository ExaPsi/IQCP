/**
 * ScfIterationTable - Display table for SCF iteration history.
 *
 * Shows iteration number, energy, energy change (delta), and DIIS error
 * in a scrollable table format. Updates live during computation and
 * highlights the converged iteration.
 *
 * @module components/scf/ScfIterationTable
 */

import { useRef, useEffect } from 'react';
import { useScfStore } from '../../stores/scfStore';
import { Math as MathRender } from '../common/Math';
import type { ScfIterationHistory } from '../../worker/protocol';

/**
 * Format energy value with appropriate precision.
 *
 * Shows 10 decimal places for Hartree values.
 */
function formatEnergy(value: number): string {
  return value.toFixed(10);
}

/**
 * Format delta or error value in scientific notation.
 *
 * Shows dash for first iteration or if value is not available.
 */
function formatDelta(value: number | undefined, isFirst: boolean): string {
  if (isFirst || value === undefined) return '-';
  if (Math.abs(value) < 1e-15) return '0.00e+00';
  return value.toExponential(2);
}

/**
 * Props for ScfIterationTable.
 */
interface ScfIterationTableProps {
  /** Maximum height before scrolling (default: 320px = ~8 rows) */
  maxHeight?: number;
}

/**
 * Individual iteration row component.
 */
function IterationRow({
  iteration,
  isFirst,
  isConverged,
}: {
  iteration: ScfIterationHistory;
  isFirst: boolean;
  isConverged: boolean;
}) {
  return (
    <tr
      className={`border-b border-slate-100 ${
        isConverged
          ? 'bg-green-50 text-green-900'
          : 'hover:bg-slate-50'
      }`}
    >
      <td className="px-3 py-2 text-center font-mono">
        {iteration.iteration}
        {isConverged && <span className="text-green-600 ml-1">*</span>}
      </td>
      <td className="px-3 py-2 font-mono text-right">
        {formatEnergy(iteration.energy)}
      </td>
      <td className="px-3 py-2 font-mono text-right">
        {formatDelta(iteration.delta, isFirst)}
      </td>
      <td className="px-3 py-2 font-mono text-right">
        {formatDelta(iteration.diisError, isFirst)}
      </td>
      <td className="px-3 py-2 text-center text-xs">
        {isConverged && (
          <span className="text-green-700 font-medium">Converged</span>
        )}
      </td>
    </tr>
  );
}

/**
 * SCF iteration history table.
 *
 * Features:
 * - Live updates during computation
 * - Auto-scrolls to show latest iteration
 * - Sticky header when scrollable
 * - Highlights converged row with green background
 * - Shows formatted energy, delta E, and DIIS error
 *
 * @example
 * ```tsx
 * <ScfIterationTable maxHeight={400} />
 * ```
 */
export function ScfIterationTable({ maxHeight = 320 }: ScfIterationTableProps) {
  const history = useScfStore((state) => state.history);
  const compute = useScfStore((state) => state.compute);

  // Ref for auto-scrolling to bottom
  const tableContainerRef = useRef<HTMLDivElement>(null);
  const isRunning = compute.status === 'running';

  // Determine if we should show converged indicator
  const isConverged = compute.status === 'success' && compute.result.converged;
  const convergedIteration = isConverged ? history.length : -1;

  // Auto-scroll to bottom when new iterations arrive (during running)
  useEffect(() => {
    if (isRunning && tableContainerRef.current) {
      tableContainerRef.current.scrollTop = tableContainerRef.current.scrollHeight;
    }
  }, [history.length, isRunning]);

  // Empty state
  if (history.length === 0) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
        <h3 className="text-sm font-semibold text-slate-700 mb-4">
          Iteration History
        </h3>
        <div className="text-center py-8 text-slate-500">
          <p className="text-sm">No iterations yet</p>
          <p className="text-xs mt-1">Click "Run SCF" to start computation</p>
        </div>
      </div>
    );
  }

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      <div className="px-4 py-3 border-b border-slate-200">
        <h3 className="text-sm font-semibold text-slate-700">
          Iteration History
          <span className="text-slate-500 font-normal ml-2">
            ({history.length} iteration{history.length !== 1 ? 's' : ''})
          </span>
        </h3>
      </div>

      <div
        ref={tableContainerRef}
        className="overflow-auto"
        style={{ maxHeight: `${maxHeight}px` }}
      >
        <table className="w-full text-sm">
          <thead className="bg-slate-50 sticky top-0 z-10">
            <tr className="text-slate-600 text-xs tracking-wider">
              <th className="px-3 py-2 text-center font-medium w-16">#</th>
              <th className="px-3 py-2 text-right font-medium">
                <span className="flex items-center justify-end gap-1">
                  <MathRender>{String.raw`E`}</MathRender>
                  <span className="text-slate-400">(Ha)</span>
                </span>
              </th>
              <th className="px-3 py-2 text-right font-medium w-24">
                <MathRender>{String.raw`\Delta E`}</MathRender>
              </th>
              <th className="px-3 py-2 text-right font-medium w-28">
                <span className="flex items-center justify-end gap-1">
                  <MathRender>{String.raw`\|\mathbf{e}\|`}</MathRender>
                  <span className="text-slate-400 text-[10px]">(DIIS)</span>
                </span>
              </th>
              <th className="px-3 py-2 text-center font-medium w-20">Status</th>
            </tr>
          </thead>
          <tbody className="text-slate-800">
            {history.map((iter, index) => (
              <IterationRow
                key={iter.iteration}
                iteration={iter}
                isFirst={index === 0}
                isConverged={iter.iteration === convergedIteration}
              />
            ))}
          </tbody>
        </table>
      </div>

      {/* Summary footer */}
      {compute.status === 'success' && (
        <div className="px-4 py-3 bg-slate-50 border-t border-slate-200 text-xs text-slate-600">
          <div className="flex justify-between">
            <span>
              Final energy:{' '}
              <span className="font-mono font-medium text-slate-800">
                {formatEnergy(compute.result.energy)} Ha
              </span>
            </span>
            <span>
              Status:{' '}
              <span className={compute.result.converged ? 'text-green-700 font-medium' : 'text-amber-700 font-medium'}>
                {compute.result.converged ? 'Converged' : 'Did not converge'}
              </span>
            </span>
          </div>
        </div>
      )}

      {/* Running indicator */}
      {isRunning && (
        <div className="px-4 py-2 bg-blue-50 border-t border-blue-100 text-xs text-blue-700 flex items-center gap-2">
          <span className="animate-spin rounded-full h-3 w-3 border-b-2 border-blue-600" />
          Computing iteration {history.length + 1}...
        </div>
      )}
    </div>
  );
}

export default ScfIterationTable;
