/**
 * Basis Explorer State Store
 *
 * Zustand store managing the state for Module D (Basis Explorer).
 * Handles element selection, basis set selection, and shell selection.
 *
 * @module stores/basisStore
 */

import { create } from 'zustand';
import { devtools } from 'zustand/middleware';

/**
 * Y-axis scale mode for plots.
 */
export type AxisScale = 'linear' | 'log';

/**
 * Configuration for one overlap trace (the B-side of a pair).
 *
 * @see US-054 Overlap vs. Distance Plot
 */
export interface OverlapTraceConfig {
  /** Atomic number for the second atom */
  elementB: number;
  /** Basis set name for the second atom */
  basisB: string;
  /** Shell index for the second atom */
  shellIndexB: number;
}

/**
 * Basis store state interface.
 */
export interface BasisState {
  /** Atomic number of the selected element (1-18) */
  selectedElement: number;
  /** Name of the selected basis set */
  selectedBasis: string;
  /** Index of the selected shell (null = none selected) */
  selectedShell: number | null;
  /** Y-axis scale mode for the radial profile plot */
  radialProfileScale: AxisScale;
  /** Whether basis set comparison mode is active */
  comparisonMode: boolean;
  /** Basis sets selected for comparison (max 5, includes primary) */
  comparisonBases: string[];
  /**
   * Exponent overrides keyed by shell index.
   * Each entry maps a shell index to an array of modified exponents.
   * If a shell index is not present, the default exponents are used.
   *
   * @see US-053 Exponent Sliders
   */
  exponentOverrides: Record<number, number[]>;
  /**
   * Whether overlap distance mode is active.
   * Mutually exclusive with comparisonMode.
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  overlapDistanceMode: boolean;
  /**
   * Overlap distance comparison traces (up to 4).
   * Each trace specifies the B-side of a shell pair.
   * The A-side always uses the current element/basis/shell selection.
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  overlapTraces: OverlapTraceConfig[];
  /**
   * Y-axis scale mode for the overlap distance plot.
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  overlapPlotScale: AxisScale;
}

/**
 * Basis store actions interface.
 */
export interface BasisActions {
  /** Set the selected element by atomic number */
  setElement: (z: number) => void;
  /** Set the selected basis set by name */
  setBasis: (name: string) => void;
  /** Set the selected shell index (or null to deselect) */
  setShell: (index: number | null) => void;
  /** Set the y-axis scale mode for the radial profile plot */
  setRadialProfileScale: (scale: AxisScale) => void;
  /** Toggle comparison mode on/off */
  toggleComparisonMode: () => void;
  /** Set comparison mode explicitly */
  setComparisonMode: (active: boolean) => void;
  /** Set the full list of comparison basis sets */
  setComparisonBases: (bases: string[]) => void;
  /** Add a basis set to comparison (max 5 total) */
  addComparisonBasis: (basis: string) => void;
  /** Remove a basis set from comparison */
  removeComparisonBasis: (basis: string) => void;
  /**
   * Set a single exponent override for a specific shell and primitive.
   * Creates the override entry if it does not yet exist (copies default exponents).
   *
   * @param shellIndex - Shell index (matching shell list order)
   * @param primitiveIndex - Index of the primitive within the shell
   * @param newAlpha - New exponent value (must be > 0)
   * @param defaultExponents - Default exponents for the shell (used to initialize the override)
   *
   * @see US-053 Exponent Sliders
   */
  setExponentOverride: (
    shellIndex: number,
    primitiveIndex: number,
    newAlpha: number,
    defaultExponents: number[],
  ) => void;
  /**
   * Reset exponent overrides for a specific shell (restores defaults).
   *
   * @see US-053 Exponent Sliders
   */
  resetExponentOverrides: (shellIndex: number) => void;
  /**
   * Reset all exponent overrides across all shells.
   *
   * @see US-053 Exponent Sliders
   */
  resetAllExponentOverrides: () => void;
  /**
   * Toggle overlap distance mode on/off.
   * Deactivates comparison mode when activating overlap mode.
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  toggleOverlapDistanceMode: () => void;
  /**
   * Set overlap distance mode explicitly.
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  setOverlapDistanceMode: (active: boolean) => void;
  /**
   * Add an overlap trace for comparison (max 4 total).
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  addOverlapTrace: (trace: OverlapTraceConfig) => void;
  /**
   * Remove an overlap trace by index.
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  removeOverlapTrace: (index: number) => void;
  /**
   * Clear all overlap traces.
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  clearOverlapTraces: () => void;
  /**
   * Set the y-axis scale for the overlap distance plot.
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  setOverlapPlotScale: (scale: AxisScale) => void;
  /** Reset to default state */
  reset: () => void;
}

/** Maximum number of basis sets in a comparison */
const MAX_COMPARISON_BASES = 5;

/** Maximum number of overlap traces */
const MAX_OVERLAP_TRACES = 4;

/**
 * Default state values.
 */
const DEFAULT_STATE: BasisState = {
  selectedElement: 1,
  selectedBasis: 'sto-3g',
  selectedShell: null,
  radialProfileScale: 'linear',
  comparisonMode: false,
  comparisonBases: [],
  exponentOverrides: {},
  overlapDistanceMode: false,
  overlapTraces: [],
  overlapPlotScale: 'linear',
};

/**
 * Basis Explorer state store.
 *
 * Manages element, basis set, and shell selection for Module D.
 *
 * @example
 * ```typescript
 * import { useBasisStore } from '@/stores/basisStore';
 *
 * function BasisControls() {
 *   const { selectedElement, selectedBasis, setElement, setBasis } = useBasisStore();
 *   // ...
 * }
 * ```
 */
export const useBasisStore = create<BasisState & BasisActions>()(
  devtools(
    (set) => ({
      // Initial state
      ...DEFAULT_STATE,

      // Actions
      setElement: (z: number) => {
        const clamped = Math.max(1, Math.min(18, Math.round(z)));
        set(
          {
            selectedElement: clamped,
            selectedShell: null,
            comparisonBases: [],
            comparisonMode: false,
            exponentOverrides: {},
            overlapDistanceMode: false,
            overlapTraces: [],
          },
          false,
          'setElement'
        );
      },

      setBasis: (name: string) => {
        set(
          { selectedBasis: name, selectedShell: null, exponentOverrides: {} },
          false,
          'setBasis'
        );
      },

      setShell: (index: number | null) => {
        set({ selectedShell: index }, false, 'setShell');
      },

      setRadialProfileScale: (scale: AxisScale) => {
        set({ radialProfileScale: scale }, false, 'setRadialProfileScale');
      },

      toggleComparisonMode: () => {
        set(
          (state) => ({ comparisonMode: !state.comparisonMode }),
          false,
          'toggleComparisonMode'
        );
      },

      setComparisonMode: (active: boolean) => {
        set({ comparisonMode: active }, false, 'setComparisonMode');
      },

      setComparisonBases: (bases: string[]) => {
        const unique = [...new Set(bases)].slice(0, MAX_COMPARISON_BASES);
        set({ comparisonBases: unique }, false, 'setComparisonBases');
      },

      addComparisonBasis: (basis: string) => {
        set(
          (state) => {
            if (state.comparisonBases.includes(basis)) return state;
            if (state.comparisonBases.length >= MAX_COMPARISON_BASES) return state;
            return { comparisonBases: [...state.comparisonBases, basis] };
          },
          false,
          'addComparisonBasis'
        );
      },

      removeComparisonBasis: (basis: string) => {
        set(
          (state) => ({
            comparisonBases: state.comparisonBases.filter((b) => b !== basis),
          }),
          false,
          'removeComparisonBasis'
        );
      },

      setExponentOverride: (
        shellIndex: number,
        primitiveIndex: number,
        newAlpha: number,
        defaultExponents: number[],
      ) => {
        set(
          (state) => {
            const existing = state.exponentOverrides[shellIndex];
            // Copy default exponents if this shell has no overrides yet
            const exponents = existing ? [...existing] : [...defaultExponents];
            exponents[primitiveIndex] = newAlpha;
            return {
              exponentOverrides: {
                ...state.exponentOverrides,
                [shellIndex]: exponents,
              },
            };
          },
          false,
          'setExponentOverride'
        );
      },

      resetExponentOverrides: (shellIndex: number) => {
        set(
          (state) => {
            const newOverrides = { ...state.exponentOverrides };
            delete newOverrides[shellIndex];
            return { exponentOverrides: newOverrides };
          },
          false,
          'resetExponentOverrides'
        );
      },

      resetAllExponentOverrides: () => {
        set({ exponentOverrides: {} }, false, 'resetAllExponentOverrides');
      },

      toggleOverlapDistanceMode: () => {
        set(
          (state) => ({
            overlapDistanceMode: !state.overlapDistanceMode,
            // Deactivate comparison mode when activating overlap mode
            comparisonMode: state.overlapDistanceMode ? state.comparisonMode : false,
          }),
          false,
          'toggleOverlapDistanceMode'
        );
      },

      setOverlapDistanceMode: (active: boolean) => {
        set(
          {
            overlapDistanceMode: active,
            // Deactivate comparison mode when activating overlap mode
            ...(active ? { comparisonMode: false } : {}),
          },
          false,
          'setOverlapDistanceMode'
        );
      },

      addOverlapTrace: (trace: OverlapTraceConfig) => {
        set(
          (state) => {
            if (state.overlapTraces.length >= MAX_OVERLAP_TRACES) return state;
            return { overlapTraces: [...state.overlapTraces, trace] };
          },
          false,
          'addOverlapTrace'
        );
      },

      removeOverlapTrace: (index: number) => {
        set(
          (state) => ({
            overlapTraces: state.overlapTraces.filter((_, i) => i !== index),
          }),
          false,
          'removeOverlapTrace'
        );
      },

      clearOverlapTraces: () => {
        set({ overlapTraces: [] }, false, 'clearOverlapTraces');
      },

      setOverlapPlotScale: (scale: AxisScale) => {
        set({ overlapPlotScale: scale }, false, 'setOverlapPlotScale');
      },

      reset: () => {
        set({ ...DEFAULT_STATE }, false, 'reset');
      },
    }),
    { name: 'basisStore' }
  )
);
