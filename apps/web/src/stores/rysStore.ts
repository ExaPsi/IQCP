/**
 * Rys Module State Store
 *
 * Zustand store managing the state for Module C (Rys Quadrature Lab).
 * Handles parameter inputs, computation status, error curve status,
 * and URL synchronization.
 *
 * @module stores/rysStore
 */

import { create } from 'zustand';
import { devtools } from 'zustand/middleware';
import type { RysComputeResult, RysErrorCurveResult } from '../worker/protocol';
import type { RysParams } from '../types/run-state';

/**
 * Target accuracy options for Rys quadrature.
 * Represented as string literals for URL encoding.
 */
export type TargetAccuracy = '1e-4' | '1e-6' | '1e-8';

/**
 * Display mode for Rys function results.
 * - 'explain': Educational explanations of quadrature
 * - 'internals': Raw computational details and sanity checks
 */
export type DisplayMode = 'explain' | 'internals';

/**
 * Computation status discriminated union for rys_compute.
 *
 * Provides type-safe handling of async computation states.
 */
export type ComputeStatus =
  | { status: 'idle' }
  | { status: 'computing' }
  | { status: 'success'; result: RysComputeResult }
  | { status: 'error'; error: string };

/**
 * Error curve computation status discriminated union for rys_error_curve.
 */
export type ErrorCurveStatus =
  | { status: 'idle' }
  | { status: 'computing' }
  | { status: 'success'; result: RysErrorCurveResult }
  | { status: 'error'; error: string };

/**
 * Rys store state interface.
 */
export interface RysState {
  /** Number of quadrature points (order), must be 1-10 */
  n: number;
  /** Parameter T (must be >= 0), typically 0-50 */
  T: number;
  /** Target accuracy for order recommendation */
  target: TargetAccuracy;
  /** Display mode: explain (educational) or internals (raw details) */
  mode: DisplayMode;
  /** Computation status for rys_compute (idle, computing, success, error) */
  compute: ComputeStatus;
  /** Error curve computation status */
  errorCurve: ErrorCurveStatus;
  /** Flag indicating if state has been initialized from URL */
  urlInitialized: boolean;
}

/**
 * Rys store actions interface.
 */
export interface RysActions {
  /** Set the order n (clamped to 1-10) */
  setN: (n: number) => void;
  /** Set the parameter T (clamped to 0-50) */
  setT: (T: number) => void;
  /** Set the target accuracy */
  setTarget: (target: TargetAccuracy) => void;
  /** Set the display mode (explain or internals) */
  setMode: (mode: DisplayMode) => void;
  /** Initialize state from URL parameters */
  initializeFromURL: (params: RysParams) => void;
  /** Mark computation as started */
  startCompute: () => void;
  /** Store successful computation result */
  setComputeResult: (result: RysComputeResult) => void;
  /** Store computation error */
  setComputeError: (error: string) => void;
  /** Mark error curve computation as started */
  startErrorCurve: () => void;
  /** Store successful error curve result */
  setErrorCurveResult: (result: RysErrorCurveResult) => void;
  /** Store error curve computation error */
  setErrorCurveError: (error: string) => void;
  /** Reset to default state */
  reset: () => void;
}

/**
 * Default state values matching DEFAULT_RYS_STATE from run-state.ts
 */
const DEFAULT_STATE: RysState = {
  n: 3,
  T: 1.0,
  target: '1e-6',
  mode: 'explain',
  compute: { status: 'idle' },
  errorCurve: { status: 'idle' },
  urlInitialized: false,
};

/**
 * Rys quadrature state store.
 *
 * Manages parameter inputs (n, T, target), computation status,
 * error curve status, and URL synchronization for Module C.
 *
 * @example
 * ```typescript
 * import { useRysStore } from '@/stores/rysStore';
 *
 * function RysControls() {
 *   const { n, T, setN, setT, compute } = useRysStore();
 *
 *   return (
 *     <>
 *       <input value={n} onChange={e => setN(Number(e.target.value))} />
 *       <input value={T} onChange={e => setT(Number(e.target.value))} />
 *       {compute.status === 'success' && (
 *         <span>Roots: {compute.result.roots.length}</span>
 *       )}
 *     </>
 *   );
 * }
 * ```
 */
export const useRysStore = create<RysState & RysActions>()(
  devtools(
    (set) => ({
      // Initial state
      ...DEFAULT_STATE,

      // Actions
      setN: (n: number) => {
        // Clamp n to valid range [1, 10]
        const clampedN = Math.max(1, Math.min(10, Math.round(n)));
        set({ n: clampedN }, false, 'setN');
      },

      setT: (T: number) => {
        // Clamp T to valid range [0, 50]
        const clampedT = Math.max(0, Math.min(50, T));
        set({ T: clampedT }, false, 'setT');
      },

      setTarget: (target: TargetAccuracy) => {
        set({ target }, false, 'setTarget');
      },

      setMode: (mode: DisplayMode) => {
        set({ mode }, false, 'setMode');
      },

      initializeFromURL: (params: RysParams) => {
        set(
          {
            n: Math.max(1, Math.min(10, Math.round(params.n))),
            T: Math.max(0, Math.min(50, params.T)),
            target: params.target,
            urlInitialized: true,
          },
          false,
          'initializeFromURL'
        );
      },

      startCompute: () => {
        set({ compute: { status: 'computing' } }, false, 'startCompute');
      },

      setComputeResult: (result: RysComputeResult) => {
        set(
          { compute: { status: 'success', result } },
          false,
          'setComputeResult'
        );
      },

      setComputeError: (error: string) => {
        set({ compute: { status: 'error', error } }, false, 'setComputeError');
      },

      startErrorCurve: () => {
        set({ errorCurve: { status: 'computing' } }, false, 'startErrorCurve');
      },

      setErrorCurveResult: (result: RysErrorCurveResult) => {
        set(
          { errorCurve: { status: 'success', result } },
          false,
          'setErrorCurveResult'
        );
      },

      setErrorCurveError: (error: string) => {
        set(
          { errorCurve: { status: 'error', error } },
          false,
          'setErrorCurveError'
        );
      },

      reset: () => {
        set({ ...DEFAULT_STATE, urlInitialized: true }, false, 'reset');
      },
    }),
    { name: 'rysStore' }
  )
);

/**
 * Get current Rys parameters in RysParams format.
 *
 * Useful for URL synchronization and deep link generation.
 *
 * @returns Current n, T, target as RysParams
 */
export function getRysParams(): RysParams {
  const { n, T, target } = useRysStore.getState();
  return { n, T, target };
}

/**
 * Convert target accuracy string to numeric value.
 *
 * @param target - Target accuracy string ('1e-4', '1e-6', '1e-8')
 * @returns Numeric value (0.0001, 0.000001, or 0.00000001)
 */
export function targetToNumber(target: TargetAccuracy): number {
  switch (target) {
    case '1e-4':
      return 1e-4;
    case '1e-6':
      return 1e-6;
    case '1e-8':
      return 1e-8;
  }
}

/**
 * Find recommended order that meets target accuracy from error curve.
 *
 * @param errorCurve - Error curve result
 * @param target - Target accuracy string
 * @returns Recommended order n, or null if none meets target
 */
export function findRecommendedOrder(
  errorCurve: RysErrorCurveResult,
  target: TargetAccuracy
): number | null {
  const targetValue = targetToNumber(target);

  for (const point of errorCurve.points) {
    if (point.maxError <= targetValue) {
      return point.n;
    }
  }

  return null; // No order meets target accuracy
}
