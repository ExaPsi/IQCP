/**
 * FrequencySummaryPanel — 4-card summary of equilibrium electronic properties.
 *
 * Rendered above the FrequencyPanel table in the Frequency tab. Shows:
 *   1. Dipole moment: |μ| in Debye with components tooltip
 *   2. Isotropic polarizability α_iso (a.u. + Å³)
 *   3. Rotor type (human-readable label)
 *   4. Rotational constants A / B / C (GHz — Infinity rendered as "—")
 *
 * Pure helpers (`dipoleMagnitude`, `isotropicPolarizability`, `formatRotorType`,
 * `formatRotConst`) are defined in `./frequencySummaryLogic.ts` so this file
 * exports only the React component.
 *
 * @module components/scf/FrequencySummaryPanel
 * @see US-102 Frequency Tab UI, AC7
 */

import type { FrequencyResult } from '../../worker/protocol';
import {
  dipoleMagnitude,
  isotropicPolarizability,
  formatRotorType,
  formatRotConst,
} from './frequencySummaryLogic';

// ============================================================================
// Props
// ============================================================================

export interface FrequencySummaryPanelProps {
  /** Frequency result from WASM. Null → panel renders nothing. */
  result: FrequencyResult | null;
}

// ============================================================================
// Component
// ============================================================================

/**
 * 4-card summary panel for the Frequency tab.
 *
 * @see US-102 Frequency Tab UI, AC7
 */
export function FrequencySummaryPanel({
  result,
}: FrequencySummaryPanelProps): JSX.Element | null {
  if (!result) return null;

  const dipoleDebye = dipoleMagnitude(result.dipoleDebye);
  const alphaIsoAu = isotropicPolarizability(result.polarizabilityAu);
  const alphaIsoAng3 = isotropicPolarizability(result.polarizabilityAng3);
  const rotorLabel = formatRotorType(result.rotorType);
  const [a, b, c] = result.rotationalConstantsGhz;
  const dipoleComponents = result.dipoleDebye.map((d) => d.toFixed(3)).join(', ');

  return (
    <div className="grid grid-cols-2 md:grid-cols-4 gap-3 mb-4">
      <div
        className="bg-white border border-slate-200 rounded-lg px-4 py-3 shadow-sm"
        title={`Components: [${dipoleComponents}] D`}
      >
        <div className="text-[11px] uppercase tracking-wide font-semibold text-slate-500 mb-1">
          Dipole moment
        </div>
        <div className="text-base font-semibold text-slate-800 tabular-nums">
          {dipoleDebye.toFixed(3)} D
        </div>
        <div className="text-[11px] text-slate-500 mt-0.5">|μ| (Debye)</div>
      </div>

      <div
        className="bg-white border border-slate-200 rounded-lg px-4 py-3 shadow-sm"
        title="Isotropic polarizability = (α_xx + α_yy + α_zz) / 3"
      >
        <div className="text-[11px] uppercase tracking-wide font-semibold text-slate-500 mb-1">
          α isotropic
        </div>
        <div className="text-base font-semibold text-slate-800 tabular-nums">
          {alphaIsoAng3.toFixed(3)} Å³
        </div>
        <div className="text-[11px] text-slate-500 mt-0.5 tabular-nums">
          {alphaIsoAu.toFixed(3)} a.u.
        </div>
      </div>

      <div className="bg-white border border-slate-200 rounded-lg px-4 py-3 shadow-sm">
        <div className="text-[11px] uppercase tracking-wide font-semibold text-slate-500 mb-1">
          Rotor type
        </div>
        <div className="text-base font-semibold text-slate-800">
          {rotorLabel}
        </div>
        <div className="text-[11px] text-slate-500 mt-0.5">
          {result.nAtoms} atom{result.nAtoms !== 1 ? 's' : ''}
        </div>
      </div>

      <div className="bg-white border border-slate-200 rounded-lg px-4 py-3 shadow-sm">
        <div className="text-[11px] uppercase tracking-wide font-semibold text-slate-500 mb-1">
          Rotational constants
        </div>
        <div className="text-base font-semibold text-slate-800 tabular-nums">
          A {formatRotConst(a)} GHz
        </div>
        <div className="text-[11px] text-slate-500 mt-0.5 tabular-nums">
          B {formatRotConst(b)} · C {formatRotConst(c)} GHz
        </div>
      </div>
    </div>
  );
}

export default FrequencySummaryPanel;
