/**
 * useBoysSweep - React hook for Boys function sweep computation.
 *
 * Computes F_m(T) over a range of T values for plotting curves.
 * Uses the same debouncing pattern as useBoys for consistency.
 *
 * @module hooks/useBoysSweep
 */

import { useEffect, useMemo, useCallback } from 'react';
import { useWorker } from './useWorker';
import { useBoysStore } from '../stores/boysStore';
import { getMaxTValue } from '../lib/boysConstants';
import type { BoysSweepResult, BoysSweepRequest } from '../worker/protocol';

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
 * Return type of useBoysSweep hook.
 */
export interface UseBoysSweepResult {
  /** Whether sweep is currently computing */
  isComputing: boolean;
  /** Error message if sweep failed */
  error: string | null;
  /** Manually trigger a sweep computation */
  compute: () => void;
}

/**
 * Debounce delay for sweep computation (ms).
 *
 * Per IQCP spec, computation debouncing is 100ms.
 */
const SWEEP_DEBOUNCE_MS = 100;

/**
 * Default sweep parameters.
 *
 * Number of points for smooth curves (fixed).
 * T range is dynamic based on m via getMaxTValue(m).
 */
const SWEEP_POINTS = 200;

/**
 * React hook for Boys function sweep computation.
 *
 * Automatically triggers debounced sweep computation when m changes.
 * Only computes when in 'sweep' view mode.
 *
 * @returns Object with isComputing, error, and compute function
 *
 * @example
 * ```typescript
 * function BoysSweepView() {
 *   const { isComputing, error } = useBoysSweep();
 *   const sweepData = useBoysStore((state) => state.sweepData);
 *
 *   if (error) return <div>Error: {error}</div>;
 *   if (isComputing) return <div>Computing...</div>;
 *   if (!sweepData) return <div>No data</div>;
 *
 *   return <BoysChart data={sweepData} />;
 * }
 * ```
 */
export function useBoysSweep(): UseBoysSweepResult {
  const { send, isReady } = useWorker();

  // Get store state and actions via selectors
  const m = useBoysStore((state) => state.m);
  const view = useBoysStore((state) => state.view);
  const sweepStatus = useBoysStore((state) => state.sweepStatus);
  const sweepError = useBoysStore((state) => state.sweepError);
  const setSweepData = useBoysStore((state) => state.setSweepData);
  const setSweepStatus = useBoysStore((state) => state.setSweepStatus);
  const setSweepError = useBoysStore((state) => state.setSweepError);

  // Create stable compute function
  const doCompute = useCallback(
    async (computeM: number) => {
      if (!isReady) return;

      setSweepStatus('pending');
      setSweepError(null);

      try {
        const request: Omit<BoysSweepRequest, 'requestId'> = {
          type: 'boys_sweep',
          m: computeM,
          T_range: [0, getMaxTValue(computeM)],
          points: SWEEP_POINTS,
        };
        const result = await send<BoysSweepResult>(request);
        setSweepData(result);
        setSweepStatus('success');
      } catch (err) {
        const message = err instanceof Error ? err.message : 'Sweep computation failed';
        setSweepError(message);
        setSweepStatus('error');
      }
    },
    [isReady, send, setSweepData, setSweepStatus, setSweepError]
  );

  // Create debounced compute function
  const debouncedCompute = useMemo(
    () => debounce(doCompute, SWEEP_DEBOUNCE_MS),
    [doCompute]
  );

  // Trigger computation when m changes (only in sweep view)
  useEffect(() => {
    if (!isReady || view !== 'sweep') return;

    debouncedCompute(m);

    return () => {
      debouncedCompute.cancel();
    };
  }, [isReady, m, view, debouncedCompute]);

  // Manual compute function (bypasses debounce)
  const compute = useCallback(() => {
    if (view === 'sweep') {
      doCompute(m);
    }
  }, [doCompute, m, view]);

  return {
    isComputing: sweepStatus === 'pending',
    error: sweepError,
    compute,
  };
}
