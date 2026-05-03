/**
 * ThermochemistryPanel — T/P slider + RRHO thermochemistry display.
 *
 * Shows all RRHO totals (ZPE, U, H, S, G, Cv, Cp) plus a 4-row
 * component-breakdown table (translational / rotational / vibrational /
 * electronic). The temperature slider (100–1000 K) and pressure input are
 * debounced at 100 ms; the parent component performs the client-side
 * recompute via `recomputeThermochemistry` (see `lib/recomputeThermochemistry.ts`).
 *
 * Units toggle (Hartree / kcal·mol⁻¹ / kJ·mol⁻¹) rescales every displayed
 * value via `lib/units.ts` helpers.
 *
 * @module components/scf/ThermochemistryPanel
 * @see US-102 Frequency Tab UI, AC4
 */

import { useCallback, useEffect, useRef, useState } from 'react';
import type {
  FrequencyResult,
  FrequencyThermochemistry,
} from '../../worker/protocol';
import {
  energyUnitsLabel,
  entropyUnitsLabel,
  type EnergyUnitsMode,
} from '../../lib/units';
import {
  formatEnergyValue,
  formatEntropyValue,
} from './thermochemistryFormat';

// ============================================================================
// Props
// ============================================================================

export interface ThermochemistryPanelProps {
  /** Frequency result (needed for mode count / empty-state detection). */
  result: FrequencyResult | null;
  /** Thermochemistry to display (typically recomputed client-side). */
  displayThermo: FrequencyThermochemistry | null;
  /** Current user-controlled temperature in K. */
  temperatureK: number;
  /** Current user-controlled pressure in Pa. */
  pressurePa: number;
  /** Current energy units for the display. */
  unitsMode: EnergyUnitsMode;
  /** T change callback (parent is responsible for debounce coordination). */
  onChangeTemperature: (temperatureK: number) => void;
  /** P change callback. */
  onChangePressure: (pressurePa: number) => void;
  /** Units toggle callback. */
  onChangeUnits: (mode: EnergyUnitsMode) => void;
}

// ============================================================================
// Constants
// ============================================================================

const T_MIN = 100;
const T_MAX = 1000;
const DEBOUNCE_MS = 100;

// ============================================================================
// Debounce helper
// ============================================================================

function debounce<T extends (...args: Parameters<T>) => void>(
  fn: T,
  delay: number
): T & { cancel: () => void } {
  let timeoutId: ReturnType<typeof setTimeout> | null = null;
  const debounced = ((...args: Parameters<T>) => {
    if (timeoutId) clearTimeout(timeoutId);
    timeoutId = setTimeout(() => {
      fn(...args);
      timeoutId = null;
    }, delay);
  }) as T & { cancel: () => void };
  debounced.cancel = () => {
    if (timeoutId) {
      clearTimeout(timeoutId);
      timeoutId = null;
    }
  };
  return debounced;
}

// ============================================================================
// Component
// ============================================================================

/**
 * Thermochemistry display panel.
 *
 * @see US-102 Frequency Tab UI, AC4
 */
export function ThermochemistryPanel({
  result,
  displayThermo,
  temperatureK,
  pressurePa,
  unitsMode,
  onChangeTemperature,
  onChangePressure,
  onChangeUnits,
}: ThermochemistryPanelProps): JSX.Element {
  // Local T for responsive slider; debounced commit to parent.
  const [localT, setLocalT] = useState(temperatureK);
  const [localP, setLocalP] = useState(pressurePa);

  useEffect(() => setLocalT(temperatureK), [temperatureK]);
  useEffect(() => setLocalP(pressurePa), [pressurePa]);

  // Debounced forwards.
  const debouncedTRef = useRef<((v: number) => void) | null>(null);
  const debouncedPRef = useRef<((v: number) => void) | null>(null);
  useEffect(() => {
    debouncedTRef.current = debounce(
      (v: number) => onChangeTemperature(v),
      DEBOUNCE_MS
    );
    debouncedPRef.current = debounce(
      (v: number) => onChangePressure(v),
      DEBOUNCE_MS
    );
  }, [onChangeTemperature, onChangePressure]);

  const handleTChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const v = Number(e.target.value);
      setLocalT(v);
      debouncedTRef.current?.(v);
    },
    []
  );

  const handlePChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const v = Number(e.target.value);
      setLocalP(v);
      debouncedPRef.current?.(v);
    },
    []
  );

  // Empty state
  if (!result || !displayThermo) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
        <h3 className="text-sm font-semibold text-slate-700 mb-1">
          Thermochemistry
        </h3>
        <div className="text-xs text-slate-500">
          Run frequency analysis to compute thermochemistry.
        </div>
      </div>
    );
  }

  const hasVibModes = displayThermo.nVibModesUsed > 0;
  const ePost = energyUnitsLabel(unitsMode);
  const sPost = entropyUnitsLabel(unitsMode);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4 space-y-3">
      <h3 className="text-sm font-semibold text-slate-700">Thermochemistry</h3>

      {/* Controls row */}
      <div className="space-y-2 text-xs">
        <label className="flex items-center gap-2 text-slate-600">
          <span className="font-medium w-20">Temperature</span>
          <input
            type="range"
            min={T_MIN}
            max={T_MAX}
            step={1}
            value={localT}
            onChange={handleTChange}
            aria-label="Temperature in Kelvin"
            className="flex-1 accent-blue-600"
          />
          <span className="tabular-nums font-semibold w-16 text-right">
            {localT.toFixed(2)} K
          </span>
        </label>

        <label className="flex items-center gap-2 text-slate-600">
          <span className="font-medium w-20">Pressure</span>
          <input
            type="number"
            min={1}
            step={1000}
            value={localP}
            onChange={handlePChange}
            aria-label="Pressure in Pascals"
            className="flex-1 px-2 py-1 border border-slate-300 rounded tabular-nums"
          />
          <span className="font-medium w-10">Pa</span>
        </label>

        {/* Units segmented control */}
        <div className="flex items-center gap-2 text-slate-600">
          <span className="font-medium w-20">Units</span>
          <div
            className="flex-1 bg-slate-200 rounded p-0.5 flex gap-0.5"
            role="group"
            aria-label="Energy units"
          >
            {(['hartree', 'kcal_mol', 'kj_mol'] as const).map((mode) => (
              <button
                key={mode}
                type="button"
                onClick={() => onChangeUnits(mode)}
                aria-pressed={unitsMode === mode}
                className={`flex-1 px-2 py-1 rounded text-[11px] font-medium transition ${
                  unitsMode === mode
                    ? 'bg-white text-slate-800 shadow-sm'
                    : 'text-slate-500 hover:text-slate-700'
                }`}
              >
                {energyUnitsLabel(mode)}
              </button>
            ))}
          </div>
        </div>
      </div>

      {/* Totals */}
      <div className="space-y-1 text-xs border-t border-slate-100 pt-2">
        {hasVibModes ? (
          <Row
            label="ZPE"
            value={formatEnergyValue(displayThermo.zpeHa, unitsMode)}
            unit={ePost}
          />
        ) : (
          <Row label="ZPE" value="—" unit="no vib modes" />
        )}
        <Row
          label="U(T)"
          value={formatEnergyValue(displayThermo.internalEnergyHa, unitsMode)}
          unit={ePost}
        />
        <Row
          label="H(T)"
          value={formatEnergyValue(displayThermo.enthalpyHa, unitsMode)}
          unit={ePost}
        />
        <Row
          label="S(T)"
          value={formatEntropyValue(displayThermo.entropyHaPerK, unitsMode)}
          unit={sPost}
        />
        <Row
          label="G(T)"
          value={formatEnergyValue(displayThermo.gibbsHa, unitsMode)}
          unit={ePost}
        />
        <Row
          label="Cv(T)"
          value={formatEntropyValue(displayThermo.cvHaPerK, unitsMode)}
          unit={sPost}
        />
        <Row
          label="Cp(T)"
          value={formatEntropyValue(displayThermo.cpHaPerK, unitsMode)}
          unit={sPost}
        />
      </div>

      {/* Component breakdown */}
      <div className="border-t border-slate-100 pt-2">
        <table className="w-full text-[11px]">
          <thead>
            <tr className="text-slate-500">
              <th className="text-left font-semibold py-1">Component</th>
              <th className="text-right font-semibold py-1">E ({ePost})</th>
              <th className="text-right font-semibold py-1">S ({sPost})</th>
              <th className="text-right font-semibold py-1">Cv ({sPost})</th>
            </tr>
          </thead>
          <tbody className="tabular-nums">
            <BreakdownRow
              label="Translational"
              e={displayThermo.eTransHa}
              s={displayThermo.sTransHaPerK}
              cv={displayThermo.cvTransHaPerK}
              mode={unitsMode}
            />
            <BreakdownRow
              label="Rotational"
              e={displayThermo.eRotHa}
              s={displayThermo.sRotHaPerK}
              cv={displayThermo.cvRotHaPerK}
              mode={unitsMode}
            />
            <BreakdownRow
              label={
                displayThermo.nImag > 0
                  ? `Vibrational (⚠ ${displayThermo.nImag} imag)`
                  : 'Vibrational'
              }
              e={displayThermo.eVibThermalHa}
              s={displayThermo.sVibHaPerK}
              cv={displayThermo.cvVibHaPerK}
              mode={unitsMode}
              warn={displayThermo.nImag > 0}
            />
            <BreakdownRow
              label="Electronic"
              e={0}
              s={displayThermo.sElecHaPerK}
              cv={0}
              mode={unitsMode}
            />
          </tbody>
        </table>
      </div>
    </div>
  );
}

// ============================================================================
// Internal row primitives
// ============================================================================

function Row({
  label,
  value,
  unit,
}: {
  label: string;
  value: string;
  unit: string;
}): JSX.Element {
  return (
    <div className="flex items-center justify-between text-slate-700">
      <span className="font-medium">{label}</span>
      <span className="tabular-nums">
        {value} <span className="text-[10px] text-slate-500">{unit}</span>
      </span>
    </div>
  );
}

function BreakdownRow({
  label,
  e,
  s,
  cv,
  mode,
  warn,
}: {
  label: string;
  e: number;
  s: number;
  cv: number;
  mode: EnergyUnitsMode;
  warn?: boolean;
}): JSX.Element {
  return (
    <tr
      className={
        warn
          ? 'text-red-600 font-semibold'
          : 'text-slate-700 border-t border-slate-100'
      }
    >
      <td className="py-1">{label}</td>
      <td className="text-right py-1">{formatEnergyValue(e, mode)}</td>
      <td className="text-right py-1">{formatEntropyValue(s, mode)}</td>
      <td className="text-right py-1">{formatEntropyValue(cv, mode)}</td>
    </tr>
  );
}

export default ThermochemistryPanel;
