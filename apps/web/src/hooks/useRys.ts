/**
 * useRys - React hook for Rys quadrature computation.
 *
 * Combines the Zustand store with the Web Worker to provide
 * automatic debounced computation on parameter changes.
 *
 * @module hooks/useRys
 */

import { useEffect, useMemo, useRef, useCallback } from 'react';
import { useWorker } from './useWorker';
import { useRysStore } from '../stores/rysStore';
import type { RysComputeResult, RysComputeRequest } from '../worker/protocol';

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
 * Return type of useRys hook.
 */
export interface UseRysResult {
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
 * React hook for Rys quadrature computation.
 *
 * Automatically triggers debounced computation when n or T changes.
 * Manages the computation lifecycle (computing -> success/error).
 *
 * @returns Object with isReady state, workerError, and compute function
 *
 * @example
 * ```typescript
 * function RysComponent() {
 *   const { isReady, workerError } = useRys();
 *   const { n, T, compute } = useRysStore();
 *
 *   if (workerError) return <div>Worker error: {workerError.message}</div>;
 *   if (!isReady) return <div>Initializing...</div>;
 *
 *   return (
 *     <div>
 *       <p>n={n}, T={T}</p>
 *       {compute.status === 'success' && (
 *         <p>Roots: {compute.result.roots.length}</p>
 *       )}
 *     </div>
 *   );
 * }
 * ```
 */
export function useRys(): UseRysResult {
  const { send, isReady, error: workerError } = useWorker();

  // Get store state and actions via selectors to minimize re-renders
  const n = useRysStore((state) => state.n);
  const T = useRysStore((state) => state.T);
  const startCompute = useRysStore((state) => state.startCompute);
  const setComputeResult = useRysStore((state) => state.setComputeResult);
  const setComputeError = useRysStore((state) => state.setComputeError);

  // Track if we've done initial computation
  const hasComputedRef = useRef(false);

  // Create stable compute function
  const doCompute = useCallback(
    async (computeN: number, computeT: number) => {
      if (!isReady) return;

      startCompute();

      try {
        const request: Omit<RysComputeRequest, 'requestId'> = {
          type: 'rys_compute',
          n: computeN,
          T: computeT,
        };
        const result = await send<RysComputeResult>(request);
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
    debouncedCompute(n, T);

    return () => {
      debouncedCompute.cancel();
    };
  }, [isReady, n, T, debouncedCompute]);

  // Track that we've computed at least once
  useEffect(() => {
    if (isReady && !hasComputedRef.current) {
      hasComputedRef.current = true;
    }
  }, [isReady]);

  // Manual compute function (bypasses debounce)
  const compute = useCallback(() => {
    doCompute(n, T);
  }, [doCompute, n, T]);

  return {
    isReady,
    workerError,
    compute,
  };
}
