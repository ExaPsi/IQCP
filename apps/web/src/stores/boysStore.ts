/**
 * Boys Module State Store
 *
 * Zustand store managing the state for Module A (Boys Function Lab).
 * Handles parameter inputs, computation status, and URL synchronization.
 *
 * @module stores/boysStore
 */

import { create } from 'zustand';
import { devtools } from 'zustand/middleware';
import type { BoysEvalResult, BoysSweepResult } from '../worker/protocol';
import type { BoysParams } from '../types/run-state';

/**
 * Computation status discriminated union.
 *
 * Provides type-safe handling of async computation states.
 */
export type ComputeStatus =
  | { status: 'idle' }
  | { status: 'computing' }
  | { status: 'success'; result: BoysEvalResult }
  | { status: 'error'; error: string };

/**
 * Sweep computation status discriminated union.
 *
 * Handles the async state for sweep computations (plotting F_m(T) curves).
 */
export type SweepStatus = 'idle' | 'pending' | 'success' | 'error';

/**
 * Display mode for Boys function results.
 * - 'explain': Educational explanations of methods
 * - 'internals': Raw computational details and metrics
 */
export type DisplayMode = 'explain' | 'internals';

/**
 * Boys store state interface.
 */
export interface BoysState {
  /** Order m of the Boys function F_m(T), must be 0-10 */
  m: number;
  /** Argument T (must be >= 0), typically 0-50 */
  T: number;
  /** View mode: single point or sweep over T range */
  view: 'single' | 'sweep';
  /** Display mode: explain (educational) or internals (raw details) */
  mode: DisplayMode;
  /** Computation status (idle, computing, success, error) */
  compute: ComputeStatus;
  /** Flag indicating if state has been initialized from URL */
  urlInitialized: boolean;
  /** Log scale toggle for sweep chart (UI preference, not persisted in URL) */
  logScale: boolean;
  /** Sweep computation data */
  sweepData: BoysSweepResult | null;
  /** Sweep computation status */
  sweepStatus: SweepStatus;
  /** Sweep computation error message */
  sweepError: string | null;
}

/**
 * Boys store actions interface.
 */
export interface BoysActions {
  /** Set the order m (clamped to 0-10) */
  setM: (m: number) => void;
  /** Set the argument T (clamped to >= 0) */
  setT: (T: number) => void;
  /** Set the view mode */
  setView: (view: 'single' | 'sweep') => void;
  /** Set the display mode (explain or internals) */
  setMode: (mode: DisplayMode) => void;
  /** Initialize state from URL parameters */
  initializeFromURL: (params: BoysParams) => void;
  /** Mark computation as started */
  startCompute: () => void;
  /** Store successful computation result */
  setComputeResult: (result: BoysEvalResult) => void;
  /** Store computation error */
  setComputeError: (error: string) => void;
  /** Reset to default state */
  reset: () => void;
  /** Toggle log scale for sweep chart */
  setLogScale: (logScale: boolean) => void;
  /** Set sweep computation data */
  setSweepData: (data: BoysSweepResult | null) => void;
  /** Set sweep computation status */
  setSweepStatus: (status: SweepStatus) => void;
  /** Set sweep computation error */
  setSweepError: (error: string | null) => void;
}

/**
 * Default state values matching DEFAULT_BOYS_STATE from run-state.ts
 */
const DEFAULT_STATE: BoysState = {
  m: 0,
  T: 1.0,
  view: 'single',
  mode: 'explain',
  compute: { status: 'idle' },
  urlInitialized: false,
  logScale: false,
  sweepData: null,
  sweepStatus: 'idle',
  sweepError: null,
};

/**
 * Boys function state store.
 *
 * Manages parameter inputs (m, T, view), computation status,
 * and URL synchronization for Module A.
 *
 * @example
 * ```typescript
 * import { useBoysStore } from '@/stores/boysStore';
 *
 * function BoysControls() {
 *   const { m, T, setM, setT, compute } = useBoysStore();
 *
 *   return (
 *     <>
 *       <input value={m} onChange={e => setM(Number(e.target.value))} />
 *       <input value={T} onChange={e => setT(Number(e.target.value))} />
 *       {compute.status === 'success' && <span>{compute.result.value}</span>}
 *     </>
 *   );
 * }
 * ```
 */
export const useBoysStore = create<BoysState & BoysActions>()(
  devtools(
    (set) => ({
      // Initial state
      ...DEFAULT_STATE,

      // Actions
      setM: (m: number) => {
        // Clamp m to valid range [0, 10]
        const clampedM = Math.max(0, Math.min(10, Math.round(m)));
        set({ m: clampedM }, false, 'setM');
      },

      setT: (T: number) => {
        // Clamp T to valid range [0, 50]
        const clampedT = Math.max(0, Math.min(50, T));
        set({ T: clampedT }, false, 'setT');
      },

      setView: (view: 'single' | 'sweep') => {
        set({ view }, false, 'setView');
      },

      setMode: (mode: DisplayMode) => {
        set({ mode }, false, 'setMode');
      },

      initializeFromURL: (params: BoysParams) => {
        set(
          {
            m: Math.max(0, Math.min(10, Math.round(params.m))),
            T: Math.max(0, Math.min(50, params.T)),
            view: params.view,
            urlInitialized: true,
          },
          false,
          'initializeFromURL'
        );
      },

      startCompute: () => {
        set({ compute: { status: 'computing' } }, false, 'startCompute');
      },

      setComputeResult: (result: BoysEvalResult) => {
        set(
          { compute: { status: 'success', result } },
          false,
          'setComputeResult'
        );
      },

      setComputeError: (error: string) => {
        set({ compute: { status: 'error', error } }, false, 'setComputeError');
      },

      reset: () => {
        set({ ...DEFAULT_STATE, urlInitialized: true }, false, 'reset');
      },

      setLogScale: (logScale: boolean) => {
        set({ logScale }, false, 'setLogScale');
      },

      setSweepData: (data: BoysSweepResult | null) => {
        set({ sweepData: data }, false, 'setSweepData');
      },

      setSweepStatus: (status: SweepStatus) => {
        set({ sweepStatus: status }, false, 'setSweepStatus');
      },

      setSweepError: (error: string | null) => {
        set({ sweepError: error }, false, 'setSweepError');
      },
    }),
    { name: 'boysStore' }
  )
);

/**
 * Get current Boys parameters in BoysParams format.
 *
 * Useful for URL synchronization and deep link generation.
 *
 * @returns Current m, T, view as BoysParams
 */
export function getBoysParams(): BoysParams {
  const { m, T, view } = useBoysStore.getState();
  return { m, T, view };
}
