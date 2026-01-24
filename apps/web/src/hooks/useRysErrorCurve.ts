/**
 * useRysErrorCurve - React hook for Rys error curve computation.
 *
 * Provides on-demand computation of error vs. quadrature order curve.
 * Unlike useRys, this is not auto-triggered on parameter changes.
 *
 * @module hooks/useRysErrorCurve
 */

import { useCallback } from 'react';
import { useWorker } from './useWorker';
import { useRysStore } from '../stores/rysStore';
import type { RysErrorCurveResult, RysErrorCurveRequest } from '../worker/protocol';

/**
 * Maximum quadrature order to test in error curve.
 */
const MAX_ORDER = 10;

/**
 * Return type of useRysErrorCurve hook.
 */
export interface UseRysErrorCurveResult {
  /** Whether the worker is ready to receive requests */
  isReady: boolean;
  /** Whether error curve computation is in progress */
  isComputing: boolean;
  /** Trigger error curve computation for current T */
  computeErrorCurve: () => Promise<void>;
}

/**
 * React hook for Rys error curve computation.
 *
 * Computes error vs. quadrature order for the current T value.
 * This is an on-demand computation triggered by user action,
 * not automatically on parameter changes.
 *
 * @returns Object with isReady state, isComputing flag, and compute function
 *
 * @example
 * ```typescript
 * function ErrorCurveSection() {
 *   const { isReady, isComputing, computeErrorCurve } = useRysErrorCurve();
 *   const { errorCurve } = useRysStore();
 *
 *   return (
 *     <div>
 *       <button
 *         onClick={computeErrorCurve}
 *         disabled={!isReady || isComputing}
 *       >
 *         {isComputing ? 'Computing...' : 'Compute Error Curve'}
 *       </button>
 *       {errorCurve.status === 'success' && (
 *         <RysErrorChart data={errorCurve.result} />
 *       )}
 *     </div>
 *   );
 * }
 * ```
 */
export function useRysErrorCurve(): UseRysErrorCurveResult {
  const { send, isReady } = useWorker();

  // Get store state and actions via selectors
  const T = useRysStore((state) => state.T);
  const errorCurveStatus = useRysStore((state) => state.errorCurve.status);
  const startErrorCurve = useRysStore((state) => state.startErrorCurve);
  const setErrorCurveResult = useRysStore((state) => state.setErrorCurveResult);
  const setErrorCurveError = useRysStore((state) => state.setErrorCurveError);

  const isComputing = errorCurveStatus === 'computing';

  // Compute error curve function
  const computeErrorCurve = useCallback(async () => {
    if (!isReady || isComputing) return;

    startErrorCurve();

    try {
      const request: Omit<RysErrorCurveRequest, 'requestId'> = {
        type: 'rys_error_curve',
        T: T,
        max_order: MAX_ORDER,
      };
      const result = await send<RysErrorCurveResult>(request);
      setErrorCurveResult(result);
    } catch (err) {
      const message = err instanceof Error ? err.message : 'Error curve computation failed';
      setErrorCurveError(message);
    }
  }, [isReady, isComputing, T, send, startErrorCurve, setErrorCurveResult, setErrorCurveError]);

  return {
    isReady,
    isComputing,
    computeErrorCurve,
  };
}
