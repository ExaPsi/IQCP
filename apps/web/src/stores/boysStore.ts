/**
 * Boys Module State Store
 *
 * Zustand store managing the state for Module C (Boys Function Lab).
 * Handles parameter inputs, computation status, and URL synchronization.
 *
 * @module stores/boysStore
 */

import { create } from 'zustand';
import { devtools } from 'zustand/middleware';
import type { BoysEvalResult, BoysSweepResult, BoysEvalAllResult } from '../worker/protocol';
import type { BoysParams } from '../types/run-state';
import { getMaxTValue } from '../lib/boysConstants';

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
 * - 'multi-order': Multi-order comparison showing F_0..F_m at fixed T
 */
export type DisplayMode = 'explain' | 'internals' | 'multi-order';

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
  /** Multi-order computation data (F_0..F_m at fixed T) */
  multiOrderData: BoysEvalAllResult | null;
  /** Multi-order computation status */
  multiOrderStatus: SweepStatus;
  /** Multi-order computation error message */
  multiOrderError: string | null;
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
  /** Set multi-order computation data */
  setMultiOrderData: (data: BoysEvalAllResult | null) => void;
  /** Set multi-order computation status */
  setMultiOrderStatus: (status: SweepStatus) => void;
  /** Set multi-order computation error */
  setMultiOrderError: (error: string | null) => void;
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
  multiOrderData: null,
  multiOrderStatus: 'idle',
  multiOrderError: null,
};

/**
 * Boys function state store.
 *
 * Manages parameter inputs (m, T, view), computation status,
 * and URL synchronization for Module C.
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
    (set, get) => ({
      // Initial state
      ...DEFAULT_STATE,

      // Actions
      setM: (m: number) => {
        // Clamp m to valid range [0, 10]
        const clampedM = Math.max(0, Math.min(10, Math.round(m)));
        set({ m: clampedM }, false, 'setM');
      },

      setT: (T: number) => {
        // Clamp T to valid range [0, maxT(m)]
        const currentM = get().m;
        const maxT = getMaxTValue(currentM);
        const clampedT = Math.max(0, Math.min(maxT, T));
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
            T: Math.max(0, Math.min(getMaxTValue(Math.round(params.m)), params.T)),
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

      setMultiOrderData: (data: BoysEvalAllResult | null) => {
        set({ multiOrderData: data }, false, 'setMultiOrderData');
      },

      setMultiOrderStatus: (status: SweepStatus) => {
        set({ multiOrderStatus: status }, false, 'setMultiOrderStatus');
      },

      setMultiOrderError: (error: string | null) => {
        set({ multiOrderError: error }, false, 'setMultiOrderError');
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
