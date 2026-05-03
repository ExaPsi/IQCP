/**
 * DftInfoPanel - Information display for DFT methods.
 *
 * Shows DFT-specific computation details when a DFT method is selected:
 * - Functional name and composition
 * - Grid type and approximate size
 * - Energy decomposition after convergence (E_1e, E_J, E_xc, E_nuc)
 *
 * Only renders when the selected method is not RHF.
 *
 * @module components/scf/DftInfoPanel
 * @see US-070 DFT UI + Method Selector + Deep Links
 */

import { getDftMethodInfo, isDftMethod } from '../../types/dft';
import type { DftMethod } from '../../types/dft';
import type { KsScfResult } from '../../worker/protocol';

/**
 * Props for the DftInfoPanel component.
 */
export interface DftInfoPanelProps {
  /** Selected DFT method */
  method: DftMethod;
  /** KS-DFT result (when available, for energy decomposition) */
  ksResult?: KsScfResult | null;
}

/**
 * Format an energy value in Hartree with 8 decimal places.
 */
function formatEnergy(e: number): string {
  return e.toFixed(8);
}

/**
 * DFT information panel with functional details and energy decomposition.
 *
 * @example
 * ```tsx
 * <DftInfoPanel method="b3lyp" ksResult={result} />
 * ```
 */
export function DftInfoPanel({ method, ksResult }: DftInfoPanelProps) {
  // Only render for DFT methods
  if (!isDftMethod(method)) return null;

  const info = getDftMethodInfo(method);

  return (
    <div className="bg-teal-50 rounded-xl border border-teal-200 p-4">
      <h3 className="text-sm font-semibold text-teal-800 mb-2">
        DFT Method Info
      </h3>

      {/* Functional name */}
      <div className="space-y-2 text-xs text-teal-700">
        <div>
          <span className="font-medium">Functional: </span>
          <span>{info.functionalName ?? info.description}</span>
        </div>

        {/* Grid info */}
        {info.usesGrid && (
          <div>
            <span className="font-medium">Grid: </span>
            <span>Becke grid, ~3,000 points/atom (Lebedev-302)</span>
          </div>
        )}

        {/* Method description for custom geometry */}
        {method !== 'rhf' && (
          <div className="bg-teal-100/50 rounded-lg p-2 text-xs text-teal-600">
            DFT methods compute integrals on-the-fly. No separate integral computation step needed.
          </div>
        )}
      </div>

      {/* Energy decomposition (after convergence) */}
      {ksResult && (
        <div className="mt-3 pt-3 border-t border-teal-200">
          <h4 className="text-xs font-semibold text-teal-800 mb-2">
            Energy Decomposition
          </h4>
          <table className="w-full text-xs">
            <tbody className="text-teal-700">
              <tr>
                <td className="py-0.5 font-medium">E_1e</td>
                <td className="py-0.5 text-right font-mono">{formatEnergy(ksResult.energy1e)} Ha</td>
              </tr>
              <tr>
                <td className="py-0.5 font-medium">E_J</td>
                <td className="py-0.5 text-right font-mono">{formatEnergy(ksResult.energyJ)} Ha</td>
              </tr>
              <tr>
                <td className="py-0.5 font-medium">E_xc</td>
                <td className="py-0.5 text-right font-mono">{formatEnergy(ksResult.energyXc)} Ha</td>
              </tr>
              <tr>
                <td className="py-0.5 font-medium">E_nuc</td>
                <td className="py-0.5 text-right font-mono">{formatEnergy(ksResult.energyNuc)} Ha</td>
              </tr>
              <tr className="border-t border-teal-200 font-semibold">
                <td className="py-1">E_total</td>
                <td className="py-1 text-right font-mono">{formatEnergy(ksResult.energy)} Ha</td>
              </tr>
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}

export default DftInfoPanel;
