/**
 * PopulationTable - Mulliken and Lowdin atomic charge display.
 *
 * Shows a table of per-atom charges and populations from population analysis
 * (Mulliken and Lowdin methods). Includes color-coded charges and a sum row.
 *
 * @module components/scf/PopulationTable
 * @see US-076 Mulliken/Lowdin Population Analysis
 */

import type { PopulationAnalysisResult } from '../../worker/protocol';

/**
 * Props for PopulationTable.
 */
export interface PopulationTableProps {
  /** Population analysis result from WASM */
  result: PopulationAnalysisResult;
  /** Additional CSS classes */
  className?: string;
}

/**
 * Format a charge value with sign and 4 decimal places.
 */
function formatCharge(charge: number): string {
  const sign = charge >= 0 ? '+' : '';
  return `${sign}${charge.toFixed(4)}`;
}

/**
 * Get CSS class for charge color coding.
 * Negative charges (electron-rich) get blue, positive (electron-poor) get red.
 */
function chargeColorClass(charge: number): string {
  if (Math.abs(charge) < 1e-6) return 'text-slate-500';
  if (charge < 0) return 'text-blue-600 font-medium';
  return 'text-red-600 font-medium';
}

/**
 * PopulationTable component.
 *
 * Displays Mulliken and Lowdin atomic charges in a compact table format.
 * Charges are color-coded: blue for negative (electron-rich), red for
 * positive (electron-poor).
 */
export function PopulationTable({ result, className = '' }: PopulationTableProps) {
  return (
    <div className={`bg-white rounded-lg border border-slate-200 ${className}`}>
      <div className="px-4 py-3 border-b border-slate-200">
        <h3 className="text-sm font-semibold text-slate-800">
          Population Analysis
        </h3>
        <p className="text-xs text-slate-500 mt-0.5">
          Mulliken &amp; Lowdin atomic charges
        </p>
      </div>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-slate-200 bg-slate-50">
              <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase tracking-wider">
                Atom
              </th>
              <th className="px-3 py-2 text-center text-xs font-medium text-slate-500 uppercase tracking-wider">
                Element
              </th>
              <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase tracking-wider">
                Mulliken
              </th>
              <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase tracking-wider">
                Lowdin
              </th>
            </tr>
          </thead>
          <tbody>
            {result.atoms.map((atom) => (
              <tr
                key={atom.atomIndex}
                className="border-b border-slate-100 hover:bg-slate-50"
              >
                <td className="px-3 py-1.5 text-slate-600 font-mono text-xs">
                  {atom.atomIndex}
                </td>
                <td className="px-3 py-1.5 text-center font-medium text-slate-800">
                  {atom.element}
                </td>
                <td className={`px-3 py-1.5 text-right font-mono text-xs ${chargeColorClass(atom.mullikenCharge)}`}>
                  {formatCharge(atom.mullikenCharge)}
                </td>
                <td className={`px-3 py-1.5 text-right font-mono text-xs ${chargeColorClass(atom.lowdinCharge)}`}>
                  {formatCharge(atom.lowdinCharge)}
                </td>
              </tr>
            ))}

            {/* Sum row */}
            <tr className="bg-slate-50 font-medium">
              <td className="px-3 py-1.5 text-slate-600 text-xs" colSpan={2}>
                Sum
              </td>
              <td className="px-3 py-1.5 text-right font-mono text-xs text-slate-600">
                {formatCharge(result.totalMullikenCharge)}
              </td>
              <td className="px-3 py-1.5 text-right font-mono text-xs text-slate-600">
                {formatCharge(result.totalLowdinCharge)}
              </td>
            </tr>
          </tbody>
        </table>
      </div>

      <div className="px-4 py-2 text-xs text-slate-400 border-t border-slate-100">
        Computed in {result.computeTimeUs < 1000
          ? `${result.computeTimeUs.toFixed(0)} us`
          : `${(result.computeTimeUs / 1000).toFixed(1)} ms`}
      </div>
    </div>
  );
}
