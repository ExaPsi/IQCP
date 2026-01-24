/**
 * useIntegralCompute - React hook for on-the-fly integral computation.
 *
 * This hook handles integral computation for custom geometry mode in the SCF
 * module. It sends geometry and basis set data to the Web Worker, receives
 * progress updates during computation, and returns an IntegralComputeResult
 * that can be used as a "virtual preset" for SCF calculations.
 *
 * The hook integrates with the scfStore for state management, allowing the
 * UI to reflect computation progress and access results.
 *
 * @module hooks/useIntegralCompute
 */

import { useCallback, useRef } from 'react';
import { useWorker, createRequestId } from './useWorker';
import { useScfStore } from '../stores/scfStore';
import type {
  IntegralComputeResult,
  IntegralProgress,
  IntegralComputeRequest,
  WorkerProgress,
  RequestId,
} from '../worker/protocol';

/**
 * Return type of useIntegralCompute hook.
 */
export interface UseIntegralComputeResult {
  /** Whether the worker is ready to receive requests */
  isReady: boolean;
  /** Worker error (if any) */
  workerError: Error | null;
  /** Whether integral computation is currently running */
  isComputing: boolean;
  /** Current progress (0-100 percent) */
  progress: number;
  /** Current computation phase */
  phase: string | null;
  /** Computation error (if any) */
  computeError: string | null;
  /** Start integral computation for current custom geometry */
  compute: () => void;
  /** Cancel the running computation */
  cancel: () => void;
}

/**
 * Type guard to check if a progress update is for integral computation.
 */
function isIntegralProgress(progress: WorkerProgress): progress is IntegralProgress {
  return progress.module === 'integral';
}

/**
 * React hook for integral computation with progress streaming.
 *
 * This hook manages the lifecycle of integral computation for custom geometries:
 * 1. Reads geometry and basis set from scfStore
 * 2. Sends 'integral_compute' request to worker
 * 3. Streams progress updates to store
 * 4. On success, calls `completeIntegralCompute` with result
 *
 * The result becomes a "virtual preset" that can be used for SCF calculations
 * exactly like pre-computed system presets.
 *
 * @returns Object with isReady, workerError, isComputing, progress, phase, error, compute, and cancel
 *
 * @example
 * ```typescript
 * function CustomGeometryPanel() {
 *   const { isReady, isComputing, progress, phase, compute, cancel } = useIntegralCompute();
 *   const integralCompute = useScfStore((s) => s.integralCompute);
 *
 *   return (
 *     <div>
 *       <button onClick={compute} disabled={!isReady || isComputing}>
 *         Compute Integrals
 *       </button>
 *       <button onClick={cancel} disabled={!isComputing}>
 *         Cancel
 *       </button>
 *       {isComputing && (
 *         <div>
 *           <p>Computing {phase}... {progress.toFixed(0)}%</p>
 *           <progress value={progress} max={100} />
 *         </div>
 *       )}
 *       {integralCompute.status === 'success' && (
 *         <p>Ready for SCF: {integralCompute.systemData.nbf} basis functions</p>
 *       )}
 *     </div>
 *   );
 * }
 * ```
 */
export function useIntegralCompute(): UseIntegralComputeResult {
  const { send, cancel: workerCancel, isReady, error: workerError } = useWorker();

  // Store selectors
  const inputMode = useScfStore((state) => state.inputMode);
  const integralCompute = useScfStore((state) => state.integralCompute);
  const integralComputeRequestId = useScfStore((state) => state.integralComputeRequestId);
  const useSpherical = useScfStore((state) => state.useSpherical);

  // Store actions
  const startIntegralCompute = useScfStore((state) => state.startIntegralCompute);
  const updateIntegralProgress = useScfStore((state) => state.updateIntegralProgress);
  const completeIntegralCompute = useScfStore((state) => state.completeIntegralCompute);
  const failIntegralCompute = useScfStore((state) => state.failIntegralCompute);
  const cancelIntegralCompute = useScfStore((state) => state.cancelIntegralCompute);

  // Ref to track current request ID for cleanup
  const requestIdRef = useRef<RequestId | null>(null);

  /**
   * Handle progress updates from the worker.
   *
   * Called during each phase of integral computation with overall progress.
   */
  const handleProgress = useCallback(
    (progress: WorkerProgress) => {
      if (!isIntegralProgress(progress)) return;

      updateIntegralProgress(progress.phase, progress.overallPercent);
    },
    [updateIntegralProgress]
  );

  /**
   * Start integral computation.
   *
   * Does nothing if:
   * - Worker is not ready
   * - Computation is already running
   * - Not in custom mode
   * - Geometry is empty
   */
  const compute = useCallback(async () => {
    // Validate preconditions
    if (!isReady) {
      console.warn('[useIntegralCompute] Worker not ready');
      return;
    }

    if (integralCompute.status === 'computing') {
      console.warn('[useIntegralCompute] Computation already running');
      return;
    }

    if (inputMode.mode !== 'custom') {
      console.warn('[useIntegralCompute] Not in custom mode');
      return;
    }

    const { geometry, basisSet } = inputMode;

    if (geometry.atoms.length === 0) {
      failIntegralCompute('No atoms in geometry');
      return;
    }

    const requestId = createRequestId();
    requestIdRef.current = requestId;
    startIntegralCompute(requestId);

    try {
      const request: Omit<IntegralComputeRequest, 'requestId'> = {
        type: 'integral_compute',
        geometry,
        basisSet,
        useSpherical,
      };

      const result = await send<IntegralComputeResult>(request, { onProgress: handleProgress });

      completeIntegralCompute(result);
    } catch (err) {
      const message = err instanceof Error ? err.message : 'Integral computation failed';
      failIntegralCompute(message);
    } finally {
      requestIdRef.current = null;
    }
  }, [
    isReady,
    integralCompute.status,
    inputMode,
    useSpherical,
    send,
    handleProgress,
    startIntegralCompute,
    completeIntegralCompute,
    failIntegralCompute,
  ]);

  /**
   * Cancel the running computation.
   *
   * Sends a cancel request to the worker for the current request ID.
   */
  const cancel = useCallback(() => {
    if (integralComputeRequestId) {
      workerCancel(integralComputeRequestId as RequestId);
      cancelIntegralCompute();
    }
  }, [integralComputeRequestId, workerCancel, cancelIntegralCompute]);

  // Derive progress and phase from store state
  const isComputing = integralCompute.status === 'computing';
  const progress = integralCompute.status === 'computing' ? integralCompute.progress : 0;
  const phase = integralCompute.status === 'computing' ? integralCompute.phase : null;
  const computeError = integralCompute.status === 'error' ? integralCompute.error : null;

  return {
    isReady,
    workerError,
    isComputing,
    progress,
    phase,
    computeError,
    compute,
    cancel,
  };
}
