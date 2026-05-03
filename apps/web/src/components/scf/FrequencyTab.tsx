/**
 * FrequencyTab — top-level panel for the Frequency workflow in Module C.
 *
 * Reads frequency state from the Zustand `scfStore` (US-103 lift from
 * component-local `useState`). Sub-components retain their existing
 * prop-based APIs (container/presenter pattern).
 *
 * Responsibilities:
 *   - Render the two-column Module C layout (sidebar + main area)
 *   - Dispatch the frequency analysis request via `useFrequency`
 *   - Coordinate a single `selectedMode` across the 4 data views
 *     (FrequencyPanel, SpectrumPlot, NormalModeViewer, ThermochemistryPanel)
 *   - Recompute thermochemistry client-side whenever T or P changes
 *     (via `lib/recomputeThermochemistry.ts`)
 *   - Show the ImaginaryWarningBanner when any frequency is negative
 *
 * @module components/scf/FrequencyTab
 * @see US-102 Frequency Tab UI
 * @see US-103 Frequency State + Deep Links (Zustand lift)
 */

import { useCallback, useEffect } from 'react';
import type {
  FrequencyProgress,
  BroadeningKind,
} from '../../worker/protocol';
import type { DftMethod } from '../../types/dft';
import type { Atom3D } from '../viewer3d/AtomSpheres';
import { useFrequency } from '../../hooks/useFrequency';
import { recomputeThermochemistry } from '../../lib/recomputeThermochemistry';
import type { EnergyUnitsMode } from '../../lib/units';
import { useScfStore } from '../../stores/scfStore';
import { FrequencyPanel } from './FrequencyPanel';
import { FrequencySummaryPanel } from './FrequencySummaryPanel';
import { SpectrumPlot } from './SpectrumPlot';
import { NormalModeViewer } from './NormalModeViewer';
import { ThermochemistryPanel } from './ThermochemistryPanel';
import { ImaginaryWarningBanner } from './ImaginaryWarningBanner';

// ============================================================================
// Props
// ============================================================================

export interface FrequencyTabProps {
  /** Whether the worker is ready. */
  isReady: boolean;
  /** Current atoms in worker format, or null. */
  workerAtoms: [number, number, number, number][] | null;
  /** Atoms3D for the 3D viewer. */
  atoms3D: Atom3D[] | null;
  /** Current basis set. */
  basisName: string;
  /** Current method. */
  method: DftMethod;
  /** Callback when the user clicks "Re-optimize" — switches to Optimize tab. */
  onNavigateToOptimize: () => void;
  /** Reused viewer3d label-toggle state. */
  showLabels: boolean;
  /** Reused viewer3d label-toggle setter. */
  onToggleLabels: (show: boolean) => void;
}

// ============================================================================
// Internal: progress bar
// ============================================================================

/** Human-readable label for each frequency pipeline phase. */
const PHASE_LABELS: Record<FrequencyProgress['phase'], string> = {
  integrals: 'Integrals + Hessian',
  nuclear_cphf: 'Nuclear CPHF',
  field_cphf: 'Field CPHF (polarizability)',
  assembly: 'Harmonic analysis',
  modes: 'Spectrum broadening',
};

function RunFrequencyPanel({
  isReady,
  hasAtoms,
  isComputing,
  progress,
  error,
  onRun,
  onCancel,
}: {
  isReady: boolean;
  hasAtoms: boolean;
  isComputing: boolean;
  progress: FrequencyProgress | null;
  error: string | null;
  onRun: () => void;
  onCancel: () => void;
}): JSX.Element {
  const percent = progress ? Math.round(progress.percent * 100) : 0;
  const phaseLabel = progress ? PHASE_LABELS[progress.phase] : null;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4 space-y-3">
      <h3 className="text-sm font-semibold text-slate-700">
        Frequency Analysis
      </h3>
      <div className="flex gap-2">
        <button
          type="button"
          onClick={isComputing ? undefined : onRun}
          disabled={!isReady || !hasAtoms || isComputing}
          aria-label={
            isComputing
              ? 'Frequency analysis running'
              : 'Run frequency analysis'
          }
          className={`flex-1 px-4 py-2.5 rounded-lg text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-emerald-500 ${
            isComputing
              ? 'bg-emerald-400 text-white cursor-not-allowed'
              : 'bg-emerald-600 text-white hover:bg-emerald-700'
          } disabled:opacity-50 disabled:cursor-not-allowed`}
        >
          {isComputing ? (
            <span className="flex items-center justify-center gap-2">
              <span className="animate-spin rounded-full h-3.5 w-3.5 border-b-2 border-white" />
              Running...
            </span>
          ) : (
            'Run Frequency Analysis'
          )}
        </button>
        {isComputing && (
          <button
            type="button"
            onClick={onCancel}
            aria-label="Cancel frequency analysis"
            className="px-3 py-2 bg-red-600 text-white rounded-lg text-xs font-semibold hover:bg-red-700 focus:outline-none focus-visible:ring-2 focus-visible:ring-red-500"
          >
            Cancel
          </button>
        )}
      </div>

      {/* Progress bar */}
      {isComputing && progress && (
        <div className="space-y-1">
          <div className="flex justify-between text-xs text-slate-600">
            <span className="font-medium">{phaseLabel}</span>
            <span className="tabular-nums">{percent}%</span>
          </div>
          <div className="h-1.5 bg-slate-200 rounded-full overflow-hidden">
            <div
              className="h-full bg-blue-500 transition-[width] duration-150"
              style={{ width: `${percent}%` }}
            />
          </div>
        </div>
      )}

      {/* Error card */}
      {error && !isComputing && (
        <div
          role="alert"
          className="bg-red-50 border border-red-300 rounded px-3 py-2 text-xs text-red-700"
        >
          <span className="font-semibold">Error:</span> {error}
        </div>
      )}

      {!hasAtoms && !isComputing && (
        <div className="text-[11px] text-slate-500">
          Select a molecule to enable frequency analysis.
        </div>
      )}
    </div>
  );
}

// ============================================================================
// Main component
// ============================================================================

/**
 * Frequency workflow tab — the top-level panel rendered by ModuleC when
 * `workflowTab === 'frequency'`.
 *
 * US-103: Reads from Zustand `scfStore.frequencyState` instead of local
 * `useState<LocalFrequencyState>`. State persists across tab switches.
 *
 * @see US-102 Frequency Tab UI
 * @see US-103 Frequency State + Deep Links
 */
export function FrequencyTab({
  isReady,
  workerAtoms,
  atoms3D,
  basisName,
  method,
  onNavigateToOptimize,
  showLabels,
  onToggleLabels,
}: FrequencyTabProps): JSX.Element {
  // ---- Zustand store selectors (US-103) -----------------------------------

  const frequencyState = useScfStore((state) => state.frequencyState);
  const setFrequencySelectedMode = useScfStore((state) => state.setFrequencySelectedMode);
  const setFrequencyTemperature = useScfStore((state) => state.setFrequencyTemperature);
  const setFrequencyPressure = useScfStore((state) => state.setFrequencyPressure);
  const setFrequencySpectrumTab = useScfStore((state) => state.setFrequencySpectrumTab);
  const setFrequencyBroadeningKind = useScfStore((state) => state.setFrequencyBroadeningKind);
  const setFrequencyFwhmCm1 = useScfStore((state) => state.setFrequencyFwhmCm1);
  const setFrequencyAmplitudeBohr = useScfStore((state) => state.setFrequencyAmplitudeBohr);
  const setFrequencyAnimationSpeed = useScfStore((state) => state.setFrequencyAnimationSpeed);
  const setFrequencyIsAnimating = useScfStore((state) => state.setFrequencyIsAnimating);
  const setFrequencyShowDisplacementArrows = useScfStore((state) => state.setFrequencyShowDisplacementArrows);
  const setFrequencyUnitsMode = useScfStore((state) => state.setFrequencyUnitsMode);
  const setFrequencyDisplayThermo = useScfStore((state) => state.setFrequencyDisplayThermo);

  // ---- Worker dispatch (US-103: refactored hook, no setter callbacks) -----

  const { run, cancel } = useFrequency({
    atoms: workerAtoms,
    basisName,
    method,
    isReady,
  });

  // ---- Client-side thermochemistry recompute ----------------------------

  // When the user drags T or P, recompute the display thermo from the
  // frozen WASM result. This is fast (<5 ms for 30 modes) so we can run
  // synchronously on every (T, P, result) change.
  useEffect(() => {
    if (!frequencyState.result) return;
    try {
      const newThermo = recomputeThermochemistry({
        result: frequencyState.result,
        temperatureK: frequencyState.temperatureK,
        pressurePa: frequencyState.pressurePa,
      });
      setFrequencyDisplayThermo(newThermo);
    } catch (err) {
      // Keep the previous displayThermo on error (e.g., invalid T/P
      // input). The error is displayed alongside the sliders.
      console.warn('[FrequencyTab] recomputeThermochemistry failed:', err);
    }
    // frequencyState.result is a reference — recompute only when it or T/P changes.
  }, [frequencyState.result, frequencyState.temperatureK, frequencyState.pressurePa, setFrequencyDisplayThermo]);

  // ---- Coordinated callbacks --------------------------------------------

  const handleSelectMode = useCallback((modeIndex: number) => {
    setFrequencySelectedMode(modeIndex);
  }, [setFrequencySelectedMode]);

  const handleChangeTab = useCallback((tab: 'ir' | 'raman') => {
    setFrequencySpectrumTab(tab);
  }, [setFrequencySpectrumTab]);

  const handleChangeBroadening = useCallback(
    (kind: BroadeningKind, fwhm: number) => {
      setFrequencyBroadeningKind(kind);
      setFrequencyFwhmCm1(fwhm);
    },
    [setFrequencyBroadeningKind, setFrequencyFwhmCm1]
  );

  const handleChangeAmplitude = useCallback((value: number) => {
    setFrequencyAmplitudeBohr(value);
  }, [setFrequencyAmplitudeBohr]);

  const handleChangeSpeed = useCallback((value: number) => {
    setFrequencyAnimationSpeed(value);
  }, [setFrequencyAnimationSpeed]);

  const handleToggleAnimation = useCallback(() => {
    setFrequencyIsAnimating(!frequencyState.isAnimating);
  }, [setFrequencyIsAnimating, frequencyState.isAnimating]);

  const handleToggleDisplacementArrows = useCallback((show: boolean) => {
    setFrequencyShowDisplacementArrows(show);
  }, [setFrequencyShowDisplacementArrows]);

  const handleChangeTemperature = useCallback((t: number) => {
    setFrequencyTemperature(t);
  }, [setFrequencyTemperature]);

  const handleChangePressure = useCallback((p: number) => {
    setFrequencyPressure(p);
  }, [setFrequencyPressure]);

  const handleChangeUnits = useCallback((mode: EnergyUnitsMode) => {
    setFrequencyUnitsMode(mode);
  }, [setFrequencyUnitsMode]);

  // Imaginary warning: jump to first imaginary mode
  const handleShowImaginary = useCallback(
    (modeIndex: number) => {
      setFrequencySelectedMode(modeIndex);
    },
    [setFrequencySelectedMode]
  );

  // ---- Derived values ---------------------------------------------------

  const hasAtoms = (workerAtoms?.length ?? 0) > 0;

  return (
    <>
      {/* Left sidebar */}
      <div className="lg:col-span-1">
        <div className="lg:sticky lg:top-4 lg:max-h-[calc(100vh-7rem)] lg:overflow-y-auto lg:overflow-x-hidden scrollbar-thin scrollbar-on-hover lg:pr-1 space-y-4">
          <RunFrequencyPanel
            isReady={isReady}
            hasAtoms={hasAtoms}
            isComputing={frequencyState.isComputing}
            progress={frequencyState.progress}
            error={frequencyState.error}
            onRun={run}
            onCancel={cancel}
          />
          <ThermochemistryPanel
            result={frequencyState.result}
            displayThermo={frequencyState.displayThermo}
            temperatureK={frequencyState.temperatureK}
            pressurePa={frequencyState.pressurePa}
            unitsMode={frequencyState.unitsMode}
            onChangeTemperature={handleChangeTemperature}
            onChangePressure={handleChangePressure}
            onChangeUnits={handleChangeUnits}
          />
        </div>
      </div>

      {/* Right main area */}
      <div className="lg:col-span-3 space-y-4">
        <ImaginaryWarningBanner
          result={frequencyState.result}
          onShowMode={handleShowImaginary}
          onNavigateToOptimize={onNavigateToOptimize}
        />

        <FrequencySummaryPanel result={frequencyState.result} />

        <NormalModeViewer
          result={frequencyState.result}
          selectedMode={frequencyState.selectedMode}
          onSelectMode={handleSelectMode}
          atoms={atoms3D}
          showLabels={showLabels}
          onToggleLabels={onToggleLabels}
          amplitudeAng={frequencyState.amplitudeBohr}
          animationSpeed={frequencyState.animationSpeed}
          isAnimating={frequencyState.isAnimating}
          showDisplacementArrows={frequencyState.showDisplacementArrows}
          onChangeAmplitude={handleChangeAmplitude}
          onChangeSpeed={handleChangeSpeed}
          onToggleAnimation={handleToggleAnimation}
          onToggleDisplacementArrows={handleToggleDisplacementArrows}
        />

        <SpectrumPlot
          result={frequencyState.result}
          selectedMode={frequencyState.selectedMode}
          onSelectMode={handleSelectMode}
          spectrumTab={frequencyState.spectrumTab}
          onChangeTab={handleChangeTab}
          broadeningKind={frequencyState.broadeningKind}
          fwhmCm1={frequencyState.fwhmCm1}
          onChangeBroadening={handleChangeBroadening}
        />

        <FrequencyPanel
          result={frequencyState.result}
          selectedMode={frequencyState.selectedMode}
          onSelectMode={handleSelectMode}
        />
      </div>
    </>
  );
}

export default FrequencyTab;
