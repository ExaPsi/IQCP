/**
 * useBoys - React hook for Boys function computation.
 *
 * Combines the Zustand store with the Web Worker to provide
 * automatic debounced computation on parameter changes.
 *
 * @module hooks/useBoys
 */

import { useEffect, useMemo, useRef, useCallback } from 'react';
import { useWorker } from './useWorker';
import { useBoysStore } from '../stores/boysStore';
import type { BoysEvalResult, BoysEvalRequest } from '../worker/protocol';

/**
 * Simple debounce implementation.
 *
 * Returns a debounced version of the function that delays
 * invocation until after `delay` milliseconds have passed
 * since the last call.
 */
function debounce<T extends (...args: Parameters<T>) => void>(
  fn: T,
  delay: number
): T & { cancel: () => void } {
  let timeoutId: ReturnType<typeof setTimeout> | null = null;

  const debounced = ((...args: Parameters<T>) => {
    if (timeoutId) {
      clearTimeout(timeoutId);
    }
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

/**
 * Return type of useBoys hook.
 */
export interface UseBoysResult {
  /** Whether the worker is ready to receive requests */
  isReady: boolean;
  /** Worker error (if any) */
  workerError: Error | null;
  /** Manually trigger a computation (normally auto-triggered on param change) */
  compute: () => void;
}

/**
 * Debounce delay for computation triggers (ms).
 *
 * Per IQCP spec, slider debouncing is 100ms to balance
 * responsiveness with avoiding excessive computation.
 */
const COMPUTE_DEBOUNCE_MS = 100;

/**
 * React hook for Boys function computation.
 *
 * Automatically triggers debounced computation when m or T changes.
 * Manages the computation lifecycle (computing → success/error).
 *
 * @returns Object with isReady state, workerError, and compute function
 *
 * @example
 * ```typescript
 * function BoysComponent() {
 *   const { isReady, workerError } = useBoys();
 *   const { m, T, compute } = useBoysStore();
 *
 *   if (workerError) return <div>Worker error: {workerError.message}</div>;
 *   if (!isReady) return <div>Initializing...</div>;
 *
 *   return (
 *     <div>
 *       <p>m={m}, T={T}</p>
 *       {compute.status === 'success' && <p>Result: {compute.result.value}</p>}
 *     </div>
 *   );
 * }
 * ```
 */
export function useBoys(): UseBoysResult {
  const { send, isReady, error: workerError } = useWorker();

  // Get store state and actions via selectors to minimize re-renders
  const m = useBoysStore((state) => state.m);
  const T = useBoysStore((state) => state.T);
  const startCompute = useBoysStore((state) => state.startCompute);
  const setComputeResult = useBoysStore((state) => state.setComputeResult);
  const setComputeError = useBoysStore((state) => state.setComputeError);

  // Track if we've done initial computation
  const hasComputedRef = useRef(false);

  // Create stable compute function
  const doCompute = useCallback(
    async (computeM: number, computeT: number) => {
      if (!isReady) return;

      startCompute();

      try {
        // Type assertion needed because TypeScript can't narrow Omit<UnionType, 'field'>
        const request: Omit<BoysEvalRequest, 'requestId'> = {
          type: 'boys_eval',
          m: computeM,
          T: computeT,
        };
        const result = await send<BoysEvalResult>(request);
        setComputeResult(result);
      } catch (err) {
        const message = err instanceof Error ? err.message : 'Computation failed';
        setComputeError(message);
      }
    },
    [isReady, send, startCompute, setComputeResult, setComputeError]
  );

  // Create debounced compute function
  const debouncedCompute = useMemo(
    () => debounce(doCompute, COMPUTE_DEBOUNCE_MS),
    [doCompute]
  );

  // Trigger computation on parameter change
  useEffect(() => {
    if (!isReady) return;

    // Always compute when parameters change
    debouncedCompute(m, T);

    return () => {
      debouncedCompute.cancel();
    };
  }, [isReady, m, T, debouncedCompute]);

  // Track that we've computed at least once
  useEffect(() => {
    if (isReady && !hasComputedRef.current) {
      hasComputedRef.current = true;
    }
  }, [isReady]);

  // Manual compute function (bypasses debounce)
  const compute = useCallback(() => {
    doCompute(m, T);
  }, [doCompute, m, T]);

  return {
    isReady,
    workerError,
    compute,
  };
}
