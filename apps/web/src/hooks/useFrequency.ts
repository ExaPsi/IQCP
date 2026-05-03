/**
 * useFrequency — React hook for running the frequency-analysis WASM pipeline.
 *
 * Thin lifecycle hook mirroring the pattern in `useOptimizer.ts`. Reads
 * frequency state from the Zustand `scfStore` and calls store actions
 * directly for state updates (US-103 lift from component-local useState).
 *
 * @module hooks/useFrequency
 * @see US-102 Frequency Tab UI (original local state design)
 * @see US-103 Frequency State + Deep Links (Zustand lift)
 */

import { useCallback, useRef } from 'react';
import { useWorker, createRequestId } from './useWorker';
import { useScfStore } from '../stores/scfStore';
import type {
  FrequencyRequest,
  FrequencyResult,
  FrequencyProgress,
  WorkerProgress,
  RequestId,
} from '../worker/protocol';
import type { DftMethod } from '../types/dft';

// ============================================================================
// Types
// ============================================================================

/**
 * Parameters for the useFrequency hook.
 *
 * US-103: Setter callbacks removed — the hook reads broadening/temperature
 * from the Zustand store and calls store actions directly.
 */
export interface UseFrequencyParams {
  /** Atom coordinates as [Z, x, y, z][] (null if not available). */
  atoms: [number, number, number, number][] | null;
  /** Basis set name (e.g., "sto-3g", "6-31g*"). */
  basisName: string;
  /** Electronic structure method. */
  method: DftMethod;
  /** Whether the worker is ready to accept requests. */
  isReady: boolean;
}

/**
 * Return type of useFrequency.
 */
export interface UseFrequencyResult {
  /** Dispatch a frequency analysis request to the worker. */
  run: () => void;
  /** Cancel a running frequency analysis (if any). */
  cancel: () => void;
  /** Whether frequency analysis is running (from store). */
  loading: boolean;
  /** Error message (from store, null if no error). */
  error: string | null;
}

// ============================================================================
// Helpers
// ============================================================================

/**
 * Type guard for frequency-phase progress events.
 */
function isFrequencyProgress(p: WorkerProgress): p is FrequencyProgress {
  return p.module === 'frequency';
}

/**
 * Map DftMethod to the method string expected by compute_frequencies.
 *
 * Matches the mapping in `useOptimizer.toOptMethod`; D3(BJ) dispersion is
 * not yet supported in the frequency pipeline, so we fall back to plain
 * B3LYP for "b3lyp-d3bj".
 */
function toFreqMethod(method: DftMethod): string {
  switch (method) {
    case 'rhf':
      return 'rhf';
    case 'lda':
      return 'lda';
    case 'b3lyp':
      return 'b3lyp';
    case 'b3lyp-d3bj':
      // D3(BJ) dispersion correction is not wired into the frequency
      // pipeline as of US-101; fall back to plain B3LYP.
      return 'b3lyp';
  }
}

// ============================================================================
// Hook implementation
// ============================================================================

/**
 * Zustand-integrated React hook for frequency analysis.
 *
 * Reads `temperatureK`, `pressurePa`, `broadeningKind`, `fwhmCm1` from
 * the store and calls `setFrequencyIsComputing`, `setFrequencyProgress`,
 * `setFrequencyResult`, `setFrequencyError`, `setFrequencySelectedMode`
 * store actions directly. Follows the `useOptimizer.ts` pattern.
 *
 * @see US-103 Frequency State + Deep Links
 */
export function useFrequency(params: UseFrequencyParams): UseFrequencyResult {
  const { atoms, basisName, method, isReady } = params;

  const { send, cancel: cancelRequest, isReady: workerReady } = useWorker();
  const requestIdRef = useRef<RequestId | null>(null);

  // Read frequency state from Zustand store
  const frequencyState = useScfStore((state) => state.frequencyState);
  const setFrequencyIsComputing = useScfStore((state) => state.setFrequencyIsComputing);
  const setFrequencyProgress = useScfStore((state) => state.setFrequencyProgress);
  const setFrequencyResult = useScfStore((state) => state.setFrequencyResult);
  const setFrequencyError = useScfStore((state) => state.setFrequencyError);
  const setFrequencySelectedMode = useScfStore((state) => state.setFrequencySelectedMode);
  const setFrequencyDisplayThermo = useScfStore((state) => state.setFrequencyDisplayThermo);

  // Wrap progress callback in a type guard that only forwards frequency events.
  const handleProgress = useCallback(
    (progress: WorkerProgress) => {
      if (!isFrequencyProgress(progress)) return;
      setFrequencyProgress(progress);
    },
    [setFrequencyProgress]
  );

  const run = useCallback(async () => {
    if (!isReady || !workerReady) return;
    if (!atoms || atoms.length === 0) return;

    const requestId = createRequestId();
    requestIdRef.current = requestId;

    // Start: clear previous state
    setFrequencyIsComputing(true);
    setFrequencyError(null);
    setFrequencyProgress(null);

    try {
      // Read current store values at dispatch time
      const { temperatureK, pressurePa, broadeningKind, fwhmCm1 } =
        useScfStore.getState().frequencyState;

      const request: Omit<FrequencyRequest, 'requestId'> = {
        type: 'frequency',
        atoms,
        basisName,
        method: toFreqMethod(method),
        temperatureK,
        pressurePa,
        broadeningKind,
        fwhmCm1,
      };

      const result = await send<FrequencyResult>(request, {
        onProgress: handleProgress,
      });

      // Success: update store with result
      setFrequencyIsComputing(false);
      setFrequencyProgress(null);
      setFrequencyError(null);
      setFrequencyResult(result);
      setFrequencySelectedMode(result.frequenciesCm1.length > 0 ? 0 : null);
      setFrequencyDisplayThermo(result.thermochemistry);
    } catch (err) {
      const message =
        err instanceof Error
          ? err.message
          : 'Frequency analysis failed (unknown error)';
      setFrequencyIsComputing(false);
      setFrequencyProgress(null);
      setFrequencyError(message);
    } finally {
      requestIdRef.current = null;
    }
  }, [
    isReady,
    workerReady,
    atoms,
    basisName,
    method,
    send,
    handleProgress,
    setFrequencyIsComputing,
    setFrequencyProgress,
    setFrequencyResult,
    setFrequencyError,
    setFrequencySelectedMode,
    setFrequencyDisplayThermo,
  ]);

  const cancel = useCallback(() => {
    const id = requestIdRef.current;
    if (id) cancelRequest(id);
  }, [cancelRequest]);

  return {
    run,
    cancel,
    loading: frequencyState.isComputing,
    error: frequencyState.error,
  };
}
