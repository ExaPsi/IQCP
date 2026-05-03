/**
 * Local state shape and default values for `FrequencyTab.tsx`.
 *
 * Isolated in its own module (no React / JSX / Plotly / Three.js imports)
 * so the shape can be unit-tested and imported from the future US-103
 * Zustand lift without pulling in the full component tree.
 *
 * @module components/scf/frequencyTabState
 * @see US-102 Frequency Tab UI
 * @see US-103 Frequency State + Deep Links (will lift this shape to Zustand)
 */

import type {
  FrequencyResult,
  FrequencyThermochemistry,
  FrequencyProgress,
  BroadeningKind,
} from '../../worker/protocol';
import type { EnergyUnitsMode } from '../../lib/units';

/**
 * Local state for the Frequency tab.
 *
 * Fields are grouped by responsibility and documented inline. US-102 stores
 * this inside a React `useState` hook on the `FrequencyTab` component.
 * US-103 will lift this shape into `scfStore.frequencyState` verbatim (or
 * nearly so — the `isComputing`, `progress`, and `error` fields will move
 * into a discriminated-union loading state then).
 *
 * @see US-102 Section 7.1 — authoritative shape for US-103 lift
 */
export interface LocalFrequencyState {
  // ---- Result + loading ----
  /** Frozen result from the last successful WASM call. null if not yet run. */
  result: FrequencyResult | null;
  /** Whether the WASM call is currently running. */
  isComputing: boolean;
  /** Progress from the most recent phase event. null if not running. */
  progress: FrequencyProgress | null;
  /** Error message from the last WASM call, or null on success / no run. */
  error: string | null;

  // ---- Shared mode selection across the 3 data views ----
  /**
   * Currently-highlighted mode index (0-based), or null if no selection.
   * Reset to 0 on each new successful result.
   */
  selectedMode: number | null;

  // ---- Client-side thermochemistry recompute ----
  /**
   * User-controlled temperature in K (default 298.15). Separate from
   * `result.thermochemistry.temperatureK` which is frozen to the original
   * WASM call's T.
   */
  temperatureK: number;
  /** User-controlled pressure in Pa (default 101325). */
  pressurePa: number;
  /**
   * Thermochemistry recomputed at the current (temperatureK, pressurePa)
   * via the TypeScript port of `crates/qc-core/src/thermochemistry.rs`.
   * Equal to `result.thermochemistry` immediately after a fresh WASM result,
   * then diverges as the user drags the T slider.
   */
  displayThermo: FrequencyThermochemistry | null;

  // ---- Spectrum plot UI state ----
  /** Currently-visible spectrum tab. */
  spectrumTab: 'ir' | 'raman';
  /** Broadening kernel choice (default 'lorentzian'). */
  broadeningKind: BroadeningKind;
  /** FWHM in cm⁻¹ for the broadening (default 8.0, range 2–20). */
  fwhmCm1: number;

  // ---- Normal-mode viewer UI state ----
  /** Displacement arrow amplitude in Ångströms (default 0.5, range 0.1–2.0). */
  amplitudeBohr: number;
  /** Animation speed multiplier (default 1.0, range 0.5–3.0). */
  animationSpeed: number;
  /**
   * Whether the displacement animation is playing. Overridden to `false`
   * at mount time if `prefers-reduced-motion: reduce` is set.
   */
  isAnimating: boolean;
  /**
   * Whether to overlay static displacement-vector arrows on the animated
   * molecule (US-102 UX refinement). Default `false` — novices see the
   * physically intuitive animated atoms first; advanced users can toggle
   * the arrows on for a directional-field view.
   */
  showDisplacementArrows: boolean;

  // ---- Units toggle ----
  /** Energy units for the thermochemistry display. */
  unitsMode: EnergyUnitsMode;
}

/**
 * Default state on first mount.
 *
 * `selectedMode` becomes 0 once a result arrives; `displayThermo` is
 * populated from the result.
 */
export const DEFAULT_LOCAL_FREQUENCY_STATE: LocalFrequencyState = {
  result: null,
  isComputing: false,
  progress: null,
  error: null,
  selectedMode: null,
  temperatureK: 298.15,
  pressurePa: 101325,
  displayThermo: null,
  spectrumTab: 'ir',
  broadeningKind: 'lorentzian',
  fwhmCm1: 8.0,
  amplitudeBohr: 0.5,
  animationSpeed: 1.0,
  isAnimating: true,
  showDisplacementArrows: false,
  unitsMode: 'kcal_mol',
};
