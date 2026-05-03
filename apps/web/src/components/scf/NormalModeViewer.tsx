/**
 * NormalModeViewer — 3D normal-mode animation viewer for the Frequency tab.
 *
 * Renders the lazy-loaded `AnimatedMolecule` (US-102 UX refinement) — atoms
 * physically oscillate along the selected normal mode, with bonds following
 * the atoms each frame. This replaced the original DisplacementArrows
 * pulse-overlay design after user testing showed the animated-atoms
 * approach (used by every mature QC viewer: GaussView, Avogadro, Jmol,
 * IQmol, ChemCraft) is much more intuitive for novices.
 *
 * Below the canvas, a control panel exposes:
 *   - Amplitude slider (0.1–2.0 Å)
 *   - Animation speed slider (0.5–3.0 Hz multiplier)
 *   - Play / Pause button
 *   - Prev / Next mode buttons (wrap around)
 *   - "Show direction arrows" toggle (off by default — opt-in for advanced
 *     users who want to see the static displacement vector field overlaid on
 *     the animation)
 *   - Current mode label "Mode i / N   ν̃ = xxx cm⁻¹"
 *
 * `AnimatedMolecule` is lazy-loaded so the Three.js chunk is only pulled in
 * when the Frequency tab is first activated.
 *
 * @module components/scf/NormalModeViewer
 * @see US-102 Frequency Tab UI, AC3
 */

import React, { Suspense, useCallback, useEffect, useMemo } from 'react';
import type { FrequencyResult } from '../../worker/protocol';
import type { Atom3D } from '../viewer3d/AtomSpheres';
import { wrapModeIndex } from '../viewer3d/displacementMath';

// ============================================================================
// Lazy viewer chunk — shares chunk with the existing Module C viewer3d entry.
// ============================================================================

const LazyAnimatedMolecule = React.lazy(() =>
  import('../viewer3d').then((m) => ({ default: m.AnimatedMolecule }))
);

// ============================================================================
// Props
// ============================================================================

export interface NormalModeViewerProps {
  /** Frequency result (needed for normal modes array). null → empty state. */
  result: FrequencyResult | null;
  /** Currently-selected mode index. null → atoms render at equilibrium. */
  selectedMode: number | null;
  /** Callback when the user uses prev/next buttons (or pressing a mode). */
  onSelectMode: (modeIndex: number) => void;
  /** Atoms to render in the 3D scene. */
  atoms: Atom3D[] | null;
  /** Display label toggle (reused from viewer3d state). */
  showLabels: boolean;
  /** Label toggle callback. */
  onToggleLabels: (show: boolean) => void;
  /** Oscillation amplitude in bohr. */
  amplitudeAng: number;
  /** Animation speed multiplier (cycles per second). */
  animationSpeed: number;
  /** Is animation playing? */
  isAnimating: boolean;
  /**
   * Whether to overlay static direction-vector arrows on the animation
   * (US-102 UX refinement). Default `false` — animated atoms are the
   * primary visualization; arrows are an opt-in extra for advanced users.
   */
  showDisplacementArrows: boolean;
  /** Amplitude change callback. */
  onChangeAmplitude: (value: number) => void;
  /** Speed change callback. */
  onChangeSpeed: (value: number) => void;
  /** Play/pause toggle callback. */
  onToggleAnimation: () => void;
  /** Show-arrows toggle callback. */
  onToggleDisplacementArrows: (show: boolean) => void;
}

// ============================================================================
// Component
// ============================================================================

/**
 * 3D normal-mode viewer with displacement arrows + control panel.
 *
 * @see US-102 Frequency Tab UI, AC3
 */
export function NormalModeViewer({
  result,
  selectedMode,
  onSelectMode,
  atoms,
  showLabels,
  onToggleLabels,
  amplitudeAng,
  animationSpeed,
  isAnimating,
  showDisplacementArrows,
  onChangeAmplitude,
  onChangeSpeed,
  onToggleAnimation,
  onToggleDisplacementArrows,
}: NormalModeViewerProps): JSX.Element {
  const nModes = result?.frequenciesCm1.length ?? 0;

  // Respect prefers-reduced-motion: if the user has the system setting on,
  // pause animation when the component first mounts. We only trigger this
  // once (via the parent callback) since the parent owns isAnimating.
  useEffect(() => {
    if (typeof window === 'undefined' || typeof window.matchMedia !== 'function') {
      return;
    }
    const mq = window.matchMedia('(prefers-reduced-motion: reduce)');
    if (mq.matches && isAnimating) {
      onToggleAnimation();
    }
    // Deliberately empty deps: only run on mount.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Selected mode label
  const modeLabel = useMemo(() => {
    if (!result || selectedMode === null || nModes === 0) {
      return 'No mode selected';
    }
    const freq = result.frequenciesCm1[selectedMode];
    if (typeof freq !== 'number') return `Mode ${selectedMode + 1} / ${nModes}`;
    const sign = freq < 0 ? '−' : '';
    return `Mode ${selectedMode + 1} / ${nModes}   ν̃ = ${sign}${Math.abs(freq).toFixed(1)} cm⁻¹${
      freq < 0 ? ' (imaginary)' : ''
    }`;
  }, [result, selectedMode, nModes]);

  // Prev / Next handlers
  const handlePrev = useCallback(() => {
    if (nModes === 0 || selectedMode === null) return;
    onSelectMode(wrapModeIndex(selectedMode - 1, nModes));
  }, [selectedMode, nModes, onSelectMode]);

  const handleNext = useCallback(() => {
    if (nModes === 0 || selectedMode === null) return;
    onSelectMode(wrapModeIndex(selectedMode + 1, nModes));
  }, [selectedMode, nModes, onSelectMode]);

  // Current mode's cartesian displacement vectors
  const modeDisplacement = useMemo<[number, number, number][]>(() => {
    if (!result || selectedMode === null) return [];
    const modes = result.normalModesCartesian;
    const mode = modes[selectedMode];
    if (!mode) return [];
    return mode.map((v) => [v[0], v[1], v[2]] as [number, number, number]);
  }, [result, selectedMode]);

  // Empty state — no atoms
  if (!atoms || atoms.length === 0) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 text-center">
        <div className="text-sm font-semibold text-slate-700 mb-1">
          Normal Mode Viewer
        </div>
        <div className="text-xs text-slate-500">
          Select a molecule and run frequency analysis to visualize normal modes.
        </div>
      </div>
    );
  }

  // Accessible aria-label for the canvas
  const canvasAriaLabel =
    result && selectedMode !== null
      ? `3D view of molecule showing displacement arrows for ${modeLabel}`
      : '3D view of molecule';

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      {/* Header */}
      <div className="px-4 py-3 bg-slate-50 border-b border-slate-200 flex items-center justify-between">
        <h3 className="text-sm font-semibold text-slate-700">
          Normal Mode Viewer
        </h3>
        <label className="flex items-center gap-1.5 text-xs text-slate-600 cursor-pointer select-none">
          <input
            type="checkbox"
            checked={showLabels}
            onChange={(e) => onToggleLabels(e.target.checked)}
            className="rounded border-slate-300 text-blue-600 focus:ring-blue-500 h-3.5 w-3.5"
          />
          Labels
        </label>
      </div>

      {/* 3D canvas area */}
      <div
        className="h-[24rem]"
        role="img"
        aria-label={canvasAriaLabel}
      >
        <Suspense
          fallback={
            <div className="flex items-center justify-center h-full bg-slate-50">
              <div className="text-center">
                <div className="animate-spin rounded-full h-6 w-6 border-b-2 border-blue-600 mx-auto mb-2" />
                <span className="text-xs text-slate-500">Loading 3D viewer...</span>
              </div>
            </div>
          }
        >
          <LazyAnimatedMolecule
            atoms={atoms}
            displacement={modeDisplacement}
            amplitude={amplitudeAng}
            speed={animationSpeed}
            isAnimating={
              isAnimating && result !== null && selectedMode !== null
            }
            showArrows={showDisplacementArrows}
            showLabels={showLabels}
            ariaLabel={canvasAriaLabel}
          />
        </Suspense>
      </div>

      {/* Control panel */}
      <div className="px-4 py-3 bg-slate-50 border-t border-slate-200 space-y-3">
        {/* Mode label + prev/next */}
        <div className="flex items-center justify-between gap-3">
          <div className="text-xs font-semibold text-slate-700 tabular-nums">
            {modeLabel}
          </div>
          <div className="flex gap-1">
            <button
              type="button"
              onClick={handlePrev}
              disabled={nModes === 0 || selectedMode === null}
              aria-label="Previous mode"
              className="px-2 py-1 bg-white border border-slate-300 rounded text-xs font-medium hover:bg-slate-50 disabled:opacity-40 disabled:cursor-not-allowed focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500"
            >
              ← Prev
            </button>
            <button
              type="button"
              onClick={handleNext}
              disabled={nModes === 0 || selectedMode === null}
              aria-label="Next mode"
              className="px-2 py-1 bg-white border border-slate-300 rounded text-xs font-medium hover:bg-slate-50 disabled:opacity-40 disabled:cursor-not-allowed focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500"
            >
              Next →
            </button>
          </div>
        </div>

        {/* Amplitude + speed sliders */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
          <label className="flex items-center gap-2 text-xs text-slate-600">
            <span className="font-medium w-20">Amplitude</span>
            <input
              type="range"
              min={0.1}
              max={2.0}
              step={0.1}
              value={amplitudeAng}
              onChange={(e) => onChangeAmplitude(Number(e.target.value))}
              aria-label="Displacement amplitude in Angstroms"
              className="flex-1 accent-blue-600"
            />
            <span className="tabular-nums font-semibold w-12 text-right">
              {amplitudeAng.toFixed(1)} Å
            </span>
          </label>
          <label className="flex items-center gap-2 text-xs text-slate-600">
            <span className="font-medium w-20">Speed</span>
            <input
              type="range"
              min={0.5}
              max={3.0}
              step={0.1}
              value={animationSpeed}
              onChange={(e) => onChangeSpeed(Number(e.target.value))}
              aria-label="Animation speed multiplier"
              className="flex-1 accent-blue-600"
            />
            <span className="tabular-nums font-semibold w-12 text-right">
              {animationSpeed.toFixed(1)}×
            </span>
          </label>
        </div>

        {/* Play / pause + show-arrows toggle */}
        <div className="flex flex-wrap items-center justify-center gap-3">
          <button
            type="button"
            onClick={onToggleAnimation}
            aria-pressed={isAnimating}
            aria-label={isAnimating ? 'Pause animation' : 'Play animation'}
            className={`px-4 py-1.5 rounded text-xs font-semibold transition ${
              isAnimating
                ? 'bg-amber-100 text-amber-800 border border-amber-300 hover:bg-amber-200'
                : 'bg-emerald-600 text-white hover:bg-emerald-700'
            } focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500`}
          >
            {isAnimating ? '⏸ Pause' : '▶ Play'}
          </button>
          <label
            className="flex items-center gap-1.5 text-xs text-slate-600 cursor-pointer select-none"
            title="Overlay static direction-vector arrows on the animation (advanced)"
          >
            <input
              type="checkbox"
              checked={showDisplacementArrows}
              onChange={(e) => onToggleDisplacementArrows(e.target.checked)}
              aria-label="Show static displacement direction arrows"
              className="rounded border-slate-300 text-blue-600 focus:ring-blue-500 h-3.5 w-3.5"
            />
            Show arrows
          </label>
        </div>
      </div>
    </div>
  );
}

export default NormalModeViewer;
