/**
 * RysDataTable - Table displaying Rys quadrature roots and weights.
 *
 * Shows roots (x_i) and weights (w_i) in a formatted table with
 * high precision display (15 significant digits).
 *
 * @module components/rys/RysDataTable
 */

import type { RysComputeResult } from '../../worker/protocol';

/**
 * Props for RysDataTable component.
 */
export interface RysDataTableProps {
  /** Computed roots and weights from rys_compute */
  result: RysComputeResult;
}

/**
 * Format a number with 15 significant digits.
 *
 * Uses toPrecision for consistent display across different magnitudes.
 */
function formatPrecision(value: number): string {
  // Handle very small values
  if (Math.abs(value) < 1e-10 && value !== 0) {
    return value.toExponential(14);
  }
  // For regular values, show high precision
  return value.toPrecision(15);
}

/**
 * RysDataTable - Displays Rys quadrature roots and weights.
 *
 * Features:
 * - Index column (1-based)
 * - Root column with subscript notation (x_i)
 * - Weight column with subscript notation (w_i)
 * - Footer showing weight sum for verification
 * - High precision (15 digits) for numerical accuracy
 *
 * @example
 * ```tsx
 * {compute.status === 'success' && (
 *   <RysDataTable result={compute.result} />
 * )}
 * ```
 */
export function RysDataTable({ result }: RysDataTableProps) {
  const { roots, weights, nroots } = result;

  // Calculate weight sum for footer
  const weightSum = weights.reduce((sum, w) => sum + w, 0);

  return (
    <div className="overflow-x-auto">
      <table className="w-full text-sm" aria-label="Rys quadrature roots and weights">
        <thead>
          <tr className="border-b border-slate-200">
            <th className="px-3 py-2 text-left font-medium text-slate-600 w-16">
              i
            </th>
            <th className="px-3 py-2 text-left font-medium text-slate-600">
              Root (x<sub>i</sub>)
            </th>
            <th className="px-3 py-2 text-left font-medium text-slate-600">
              Weight (w<sub>i</sub>)
            </th>
          </tr>
        </thead>
        <tbody>
          {roots.map((root, index) => (
            <tr
              key={index}
              className="border-b border-slate-100 hover:bg-slate-50 transition-colors"
            >
              <td className="px-3 py-2 font-mono text-slate-500">
                {index + 1}
              </td>
              <td className="px-3 py-2 font-mono text-slate-800 break-all">
                {formatPrecision(root)}
              </td>
              <td className="px-3 py-2 font-mono text-slate-800 break-all">
                {formatPrecision(weights[index])}
              </td>
            </tr>
          ))}
        </tbody>
        <tfoot>
          <tr className="border-t border-slate-200 bg-slate-50">
            <td className="px-3 py-2 font-medium text-slate-600" colSpan={2}>
              Sum of weights
            </td>
            <td className="px-3 py-2 font-mono text-slate-800">
              {formatPrecision(weightSum)}
            </td>
          </tr>
        </tfoot>
      </table>

      {/* Summary text */}
      <div className="mt-3 text-xs text-slate-500">
        <p>
          {nroots} quadrature {nroots === 1 ? 'point' : 'points'} at T = {result.t.toFixed(4)}
        </p>
        <p className="mt-1">
          Roots are in the interval (0, 1) and weights are positive.
          The weight sum equals F<sub>0</sub>(T), the zeroth-order Boys function.
        </p>
      </div>
    </div>
  );
}

export default RysDataTable;
