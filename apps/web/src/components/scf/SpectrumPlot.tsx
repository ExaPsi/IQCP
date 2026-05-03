/**
 * SpectrumPlot — tabbed IR/Raman spectrum plot with broadening controls.
 *
 * Renders a Plotly scatter plot with:
 *   - Reversed x-axis (4500 → 0 cm⁻¹, conventional orientation)
 *   - Stick spectrum trace (one bar per mode)
 *   - Broadened continuous trace (Lorentzian or Gaussian, FWHM 2–20 cm⁻¹)
 *   - Tab bar: IR / Raman
 *   - Lorentzian / Gaussian toggle
 *   - FWHM slider
 *   - Click-to-select: clicking near a stick selects that mode
 *
 * The broadened curve is recomputed **client-side** from the stick amplitudes
 * whenever the FWHM or kernel changes — avoids round-tripping to WASM for
 * slider drags.
 *
 * @module components/scf/SpectrumPlot
 * @see US-102 Frequency Tab UI, AC2
 */

import { useCallback, useMemo, useRef, useState, useEffect } from 'react';
import type { Data, Layout, PlotMouseEvent } from 'plotly.js';
import { PlotPanel } from '../common/PlotPanel';
import type {
  FrequencyResult,
  BroadeningKind,
} from '../../worker/protocol';
import {
  computeBroadenedSpectrum,
  findNearestMode,
} from './spectrumMath';

// ============================================================================
// Public types
// ============================================================================

export interface SpectrumPlotProps {
  /** Frequency result from WASM, or null for empty state. */
  result: FrequencyResult | null;
  /** Currently-selected mode index, for annotation. */
  selectedMode: number | null;
  /** Click-to-select callback (nearest mode within threshold). */
  onSelectMode: (modeIndex: number) => void;
  /** Visible spectrum tab (IR or Raman). */
  spectrumTab: 'ir' | 'raman';
  /** Tab change callback. */
  onChangeTab: (tab: 'ir' | 'raman') => void;
  /** Broadening kernel choice. */
  broadeningKind: BroadeningKind;
  /** Current FWHM in cm⁻¹. */
  fwhmCm1: number;
  /** Broadening or FWHM change callback (called with both fields). */
  onChangeBroadening: (kind: BroadeningKind, fwhmCm1: number) => void;
}

// ============================================================================
// Internal debounce helper
// ============================================================================

/**
 * Simple debounce — same pattern as `useBoys.ts` / `ModuleC.tsx`.
 */
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
// Constants
// ============================================================================

const IR_COLOR = '#3b82f6'; // blue-500
const RAMAN_COLOR = '#10b981'; // emerald-500
const SELECTED_COLOR = '#ef4444'; // red-500

const MIN_FWHM = 2;
const MAX_FWHM = 20;
const DEBOUNCE_MS = 100;

// ============================================================================
// Component
// ============================================================================

/**
 * Tabbed IR/Raman spectrum plot.
 *
 * @see US-102 Frequency Tab UI, AC2
 */
export function SpectrumPlot({
  result,
  selectedMode,
  onSelectMode,
  spectrumTab,
  onChangeTab,
  broadeningKind,
  fwhmCm1,
  onChangeBroadening,
}: SpectrumPlotProps): JSX.Element {
  // Local FWHM value for responsive slider drag; debounced-forwarded to parent.
  const [localFwhm, setLocalFwhm] = useState(fwhmCm1);

  // Keep local FWHM in sync with prop if the parent changes it externally.
  useEffect(() => {
    setLocalFwhm(fwhmCm1);
  }, [fwhmCm1]);

  // Debounced callback to commit FWHM changes to the parent.
  const debouncedRef = useRef<((kind: BroadeningKind, fwhm: number) => void) | null>(
    null
  );
  useEffect(() => {
    debouncedRef.current = debounce(
      (kind: BroadeningKind, fwhm: number) => onChangeBroadening(kind, fwhm),
      DEBOUNCE_MS
    );
    return () => {
      // No cancel() accessor here since we rebuild on deps change.
    };
  }, [onChangeBroadening]);

  // Broadened trace recompute (pure client-side)
  const broadenedTrace = useMemo(() => {
    if (!result) return null;
    const freqs = result.frequenciesCm1;
    const amps =
      spectrumTab === 'ir'
        ? result.irIntensitiesKmPerMol
        : result.ramanActivitiesA4Amu;
    return computeBroadenedSpectrum(freqs, amps, broadeningKind, localFwhm);
  }, [result, spectrumTab, broadeningKind, localFwhm]);

  // Stick trace — one vertical bar per real mode
  const stickTrace = useMemo(() => {
    if (!result) return null;
    const freqs = result.frequenciesCm1;
    const amps =
      spectrumTab === 'ir'
        ? result.irIntensitiesKmPerMol
        : result.ramanActivitiesA4Amu;
    const xs: number[] = [];
    const ys: number[] = [];
    const colors: string[] = [];
    const base = spectrumTab === 'ir' ? IR_COLOR : RAMAN_COLOR;
    for (let i = 0; i < freqs.length; i++) {
      if (freqs[i] <= 0) continue; // skip imaginary
      xs.push(freqs[i]);
      ys.push(amps[i] ?? 0);
      colors.push(i === selectedMode ? SELECTED_COLOR : base);
    }
    return { xs, ys, colors };
  }, [result, selectedMode, spectrumTab]);

  // Plotly data traces
  const plotData: Data[] = useMemo(() => {
    if (!stickTrace || !broadenedTrace) return [];
    const base = spectrumTab === 'ir' ? IR_COLOR : RAMAN_COLOR;
    return [
      {
        type: 'scatter',
        mode: 'lines',
        name: 'Broadened',
        x: broadenedTrace.grid,
        y: broadenedTrace.intensity,
        line: { color: base, width: 1.5 },
        hoverinfo: 'x+y',
      },
      {
        type: 'bar',
        name: 'Sticks',
        x: stickTrace.xs,
        y: stickTrace.ys,
        marker: {
          color: stickTrace.colors,
          line: { width: 0 },
        },
        width: 4,
        hoverinfo: 'x+y',
      },
    ];
  }, [stickTrace, broadenedTrace, spectrumTab]);

  // Plotly layout — reversed x-axis is AC2
  const layout: Partial<Layout> = useMemo(() => {
    const yLabel =
      spectrumTab === 'ir'
        ? 'IR intensity (km/mol)'
        : 'Raman activity (Å⁴/amu)';
    return {
      xaxis: {
        title: { text: 'Wavenumber (cm⁻¹)' },
        // Reversed: high frequencies on the left, low on the right.
        range: [4500, 0],
        autorange: false,
        zeroline: false,
      },
      yaxis: {
        title: { text: yLabel },
        rangemode: 'tozero',
      },
      showlegend: false,
      margin: { t: 20, r: 20, b: 50, l: 70 },
    };
  }, [spectrumTab]);

  // Click handler — find the nearest stick and select it
  const handlePlotClick = useCallback(
    (event: PlotMouseEvent) => {
      if (!result) return;
      const point = event.points[0];
      if (!point) return;
      const clickedX = point.x;
      if (typeof clickedX !== 'number') return;
      const idx = findNearestMode(result.frequenciesCm1, clickedX);
      if (idx >= 0) onSelectMode(idx);
    },
    [result, onSelectMode]
  );

  // FWHM slider handler (debounced)
  const handleFwhmChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const value = Number(e.target.value);
      setLocalFwhm(value);
      debouncedRef.current?.(broadeningKind, value);
    },
    [broadeningKind]
  );

  // Kind toggle — immediate (not debounced)
  const handleKindChange = useCallback(
    (kind: BroadeningKind) => {
      setLocalFwhm((cur) => {
        onChangeBroadening(kind, cur);
        return cur;
      });
    },
    [onChangeBroadening]
  );

  // Empty state
  if (!result) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 text-center">
        <div className="text-sm font-semibold text-slate-700 mb-1">
          IR / Raman Spectrum
        </div>
        <div className="text-xs text-slate-500">
          Run frequency analysis to see simulated spectra.
        </div>
      </div>
    );
  }

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      {/* Header: tab bar + broadening controls */}
      <div className="px-4 py-3 bg-slate-50 border-b border-slate-200 flex flex-col lg:flex-row lg:items-center gap-3 lg:justify-between">
        <h3 className="text-sm font-semibold text-slate-700">
          Simulated Spectrum
        </h3>
        <div className="flex flex-wrap items-center gap-3">
          {/* IR / Raman tab bar */}
          <div
            className="bg-slate-200 rounded-lg p-0.5 flex gap-0.5"
            role="tablist"
            aria-label="Spectrum type"
          >
            <button
              type="button"
              role="tab"
              aria-selected={spectrumTab === 'ir'}
              onClick={() => onChangeTab('ir')}
              className={`px-3 py-1 rounded text-xs font-semibold transition ${
                spectrumTab === 'ir'
                  ? 'bg-white text-slate-800 shadow-sm'
                  : 'text-slate-500 hover:text-slate-700'
              }`}
            >
              IR
            </button>
            <button
              type="button"
              role="tab"
              aria-selected={spectrumTab === 'raman'}
              onClick={() => onChangeTab('raman')}
              className={`px-3 py-1 rounded text-xs font-semibold transition ${
                spectrumTab === 'raman'
                  ? 'bg-white text-slate-800 shadow-sm'
                  : 'text-slate-500 hover:text-slate-700'
              }`}
            >
              Raman
            </button>
          </div>

          {/* Lorentzian / Gaussian segmented control */}
          <div
            className="bg-slate-200 rounded-lg p-0.5 flex gap-0.5"
            role="group"
            aria-label="Broadening kernel"
          >
            <button
              type="button"
              onClick={() => handleKindChange('lorentzian')}
              aria-pressed={broadeningKind === 'lorentzian'}
              className={`px-2.5 py-1 rounded text-xs font-medium transition ${
                broadeningKind === 'lorentzian'
                  ? 'bg-white text-slate-800 shadow-sm'
                  : 'text-slate-500 hover:text-slate-700'
              }`}
            >
              Lorentzian
            </button>
            <button
              type="button"
              onClick={() => handleKindChange('gaussian')}
              aria-pressed={broadeningKind === 'gaussian'}
              className={`px-2.5 py-1 rounded text-xs font-medium transition ${
                broadeningKind === 'gaussian'
                  ? 'bg-white text-slate-800 shadow-sm'
                  : 'text-slate-500 hover:text-slate-700'
              }`}
            >
              Gaussian
            </button>
          </div>

          {/* FWHM slider */}
          <label className="flex items-center gap-2 text-xs text-slate-600">
            <span className="font-medium">FWHM</span>
            <input
              type="range"
              min={MIN_FWHM}
              max={MAX_FWHM}
              step={0.5}
              value={localFwhm}
              onChange={handleFwhmChange}
              aria-label="Full width at half maximum, in cm⁻¹"
              className="w-32 accent-blue-600"
            />
            <span className="tabular-nums font-semibold w-12 text-right">
              {localFwhm.toFixed(1)} cm⁻¹
            </span>
          </label>
        </div>
      </div>

      {/* Plot */}
      <div className="p-2">
        <PlotPanel
          data={plotData}
          layout={layout}
          onClick={handlePlotClick}
          ariaLabel={`Simulated ${spectrumTab === 'ir' ? 'IR' : 'Raman'} spectrum`}
          minHeight={320}
        />
      </div>
    </div>
  );
}

export default SpectrumPlot;
