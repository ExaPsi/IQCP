/**
 * useBoysEvalAll - React hook for multi-order Boys function computation.
 *
 * Computes F_0(T) through F_m(T) at a single T value for the
 * multi-order comparison view. Uses the same debouncing pattern
 * as useBoys and useBoysSweep for consistency.
 *
 * @module hooks/useBoysEvalAll
 */

import { useEffect, useMemo, useCallback } from 'react';
import { useWorker } from './useWorker';
import { useBoysStore } from '../stores/boysStore';
import type { BoysEvalAllResult, BoysEvalAllRequest } from '../worker/protocol';

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
 * Return type of useBoysEvalAll hook.
 */
export interface UseBoysEvalAllResult {
  /** Whether multi-order evaluation is currently computing */
  isComputing: boolean;
  /** Error message if computation failed */
  error: string | null;
  /** Manually trigger a computation */
  compute: () => void;
}

/**
 * Debounce delay for computation triggers (ms).
 *
 * Per IQCP spec, computation debouncing is 100ms.
 */
const COMPUTE_DEBOUNCE_MS = 100;

/**
 * React hook for multi-order Boys function computation.
 *
 * Automatically triggers debounced computation when m or T changes,
 * but only when in 'multi-order' display mode. Computes F_0(T)
 * through F_m(T) at the current T value.
 *
 * @returns Object with isComputing, error, and compute function
 *
 * @example
 * ```typescript
 * function BoysMultiOrderView() {
 *   const { isComputing, error } = useBoysEvalAll();
 *   const multiOrderData = useBoysStore((state) => state.multiOrderData);
 *
 *   if (error) return <div>Error: {error}</div>;
 *   if (isComputing) return <div>Computing...</div>;
 *   if (!multiOrderData) return <div>No data</div>;
 *
 *   return <BoysMultiOrderPanel data={multiOrderData} />;
 * }
 * ```
 */
export function useBoysEvalAll(): UseBoysEvalAllResult {
  const { send, isReady } = useWorker();

  // Get store state and actions via selectors
  const m = useBoysStore((state) => state.m);
  const T = useBoysStore((state) => state.T);
  const mode = useBoysStore((state) => state.mode);
  const multiOrderStatus = useBoysStore((state) => state.multiOrderStatus);
  const multiOrderError = useBoysStore((state) => state.multiOrderError);
  const setMultiOrderData = useBoysStore((state) => state.setMultiOrderData);
  const setMultiOrderStatus = useBoysStore((state) => state.setMultiOrderStatus);
  const setMultiOrderError = useBoysStore((state) => state.setMultiOrderError);

  // Create stable compute function
  const doCompute = useCallback(
    async (computeM: number, computeT: number) => {
      if (!isReady) return;

      setMultiOrderStatus('pending');
      setMultiOrderError(null);

      try {
        const request: Omit<BoysEvalAllRequest, 'requestId'> = {
          type: 'boys_eval_all',
          mMax: computeM,
          T: computeT,
        };
        const result = await send<BoysEvalAllResult>(request);
        setMultiOrderData(result);
        setMultiOrderStatus('success');
      } catch (err) {
        const message = err instanceof Error ? err.message : 'Multi-order computation failed';
        setMultiOrderError(message);
        setMultiOrderStatus('error');
      }
    },
    [isReady, send, setMultiOrderData, setMultiOrderStatus, setMultiOrderError]
  );

  // Create debounced compute function
  const debouncedCompute = useMemo(
    () => debounce(doCompute, COMPUTE_DEBOUNCE_MS),
    [doCompute]
  );

  // Trigger computation when m or T changes (only in multi-order mode)
  useEffect(() => {
    if (!isReady || mode !== 'multi-order') return;

    debouncedCompute(m, T);

    return () => {
      debouncedCompute.cancel();
    };
  }, [isReady, m, T, mode, debouncedCompute]);

  // Manual compute function (bypasses debounce)
  const compute = useCallback(() => {
    if (mode === 'multi-order') {
      doCompute(m, T);
    }
  }, [doCompute, m, T, mode]);

  return {
    isComputing: multiOrderStatus === 'pending',
    error: multiOrderError,
    compute,
  };
}
