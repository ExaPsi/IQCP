/**
 * useScf - React hook for SCF computation.
 *
 * Unlike useBoys/useRys which auto-compute on parameter changes,
 * this hook provides explicit run() and cancelRun() methods for
 * controlling SCF computation. It handles progress streaming via
 * the worker's onProgress callback.
 *
 * Always requests matrix data (S, Hcore, F, D) and orbital energies
 * to support the internals mode (US-018). The payload is small (~2KB max)
 * so including matrices has negligible impact on performance.
 *
 * @module hooks/useScf
 */

import { useCallback, useRef } from 'react';
import { useWorker, createRequestId } from './useWorker';
import { useScfStore } from '../stores/scfStore';
import type {
  ScfRunResult,
  ScfIterationProgress,
  ScfRunRequest,
  WorkerProgress,
  RequestId,
} from '../worker/protocol';

/**
 * Return type of useScf hook.
 */
export interface UseScfResult {
  /** Whether the worker is ready to receive requests */
  isReady: boolean;
  /** Worker error (if any) */
  workerError: Error | null;
  /** Whether SCF is currently running */
  isRunning: boolean;
  /** Start an SCF computation */
  run: () => void;
  /** Cancel the running computation */
  cancelRun: () => void;
}

/**
 * Type guard to check if a progress update is for SCF.
 */
function isScfProgress(progress: WorkerProgress): progress is ScfIterationProgress {
  return progress.module === 'scf';
}

/**
 * React hook for SCF computation with progress streaming.
 *
 * Unlike useBoys/useRys, this does NOT auto-compute on parameter change.
 * Instead, computation is triggered by the explicit run() method.
 * Progress updates are streamed via the worker's onProgress callback
 * and stored in the scfStore for live UI updates.
 *
 * @returns Object with isReady, workerError, isRunning, run, and cancelRun
 *
 * @example
 * ```typescript
 * function ScfComponent() {
 *   const { isReady, isRunning, run, cancelRun } = useScf();
 *   const history = useScfStore((s) => s.history);
 *   const compute = useScfStore((s) => s.compute);
 *
 *   return (
 *     <div>
 *       <button onClick={run} disabled={!isReady || isRunning}>
 *         Run SCF
 *       </button>
 *       <button onClick={cancelRun} disabled={!isRunning}>
 *         Cancel
 *       </button>
 *       {isRunning && <p>Iteration {history.length}...</p>}
 *       {compute.status === 'success' && (
 *         <p>Converged: {compute.result.energy} Ha</p>
 *       )}
 *     </div>
 *   );
 * }
 * ```
 */
export function useScf(): UseScfResult {
  const { send, cancel, isReady, error: workerError } = useWorker();

  // Store selectors - get individual values to minimize re-renders
  const systemId = useScfStore((state) => state.systemId);
  const convergenceProfile = useScfStore((state) => state.convergenceProfile);
  const maxIterations = useScfStore((state) => state.maxIterations);
  const useDiis = useScfStore((state) => state.useDiis);
  const damp = useScfStore((state) => state.damp);
  const computeStatus = useScfStore((state) => state.compute.status);
  const runningRequestId = useScfStore((state) => state.runningRequestId);

  // Store actions
  const startRun = useScfStore((state) => state.startRun);
  const addIteration = useScfStore((state) => state.addIteration);
  const setRunResult = useScfStore((state) => state.setRunResult);
  const setRunError = useScfStore((state) => state.setRunError);
  const setRunCancelled = useScfStore((state) => state.setRunCancelled);

  // Ref to track current request ID for cleanup
  const requestIdRef = useRef<RequestId | null>(null);

  /**
   * Handle progress updates from the worker.
   *
   * Called for each SCF iteration with energy and convergence info.
   */
  const handleProgress = useCallback(
    (progress: WorkerProgress) => {
      if (!isScfProgress(progress)) return;

      addIteration({
        iteration: progress.iteration,
        energy: progress.energy,
        delta: progress.delta,
        diisError: progress.diisError,
      });
    },
    [addIteration]
  );

  /**
   * Start an SCF computation.
   *
   * Does nothing if worker is not ready or computation is already running.
   */
  const run = useCallback(async () => {
    if (!isReady || computeStatus === 'running') return;

    const requestId = createRequestId();
    requestIdRef.current = requestId;
    startRun(requestId);

    try {
      // Type assertion needed because TypeScript can't narrow Omit<UnionType, 'field'>
      // Always include matrices to support internals mode (US-018).
      // Payload is small (~2KB max for our preset systems) so this has
      // negligible impact on performance.
      const request: Omit<ScfRunRequest, 'requestId'> = {
        type: 'scf_run',
        systemId,
        options: {
          convergenceProfile,
          maxIterations,
          useDiis,
          damp: damp > 0 ? damp : undefined,
          includeMatrices: true,
        },
      };
      const result = await send<ScfRunResult>(request, { onProgress: handleProgress });

      // Check if computation was aborted
      if (result.aborted) {
        setRunCancelled();
      } else {
        setRunResult(result);
      }
    } catch (err) {
      const message = err instanceof Error ? err.message : 'SCF computation failed';
      setRunError(message);
    } finally {
      requestIdRef.current = null;
    }
  }, [
    isReady,
    computeStatus,
    systemId,
    convergenceProfile,
    maxIterations,
    useDiis,
    damp,
    send,
    handleProgress,
    startRun,
    setRunResult,
    setRunCancelled,
    setRunError,
  ]);

  /**
   * Cancel the running computation.
   *
   * Sends a cancel request to the worker for the current request ID.
   */
  const cancelRun = useCallback(() => {
    if (runningRequestId) {
      cancel(runningRequestId as RequestId);
    }
  }, [runningRequestId, cancel]);

  return {
    isReady,
    workerError,
    isRunning: computeStatus === 'running',
    run,
    cancelRun,
  };
}
