/**
 * useOptimizer - React hook for geometry optimization.
 *
 * Orchestrates L-BFGS geometry optimization via the WASM optimize_geometry
 * function. Progress is streamed step-by-step and stored in scfStore's
 * optimizationState.
 *
 * Pattern follows usePesScan: explicit trigger, progress streaming, cancellation.
 *
 * @module hooks/useOptimizer
 * @see US-075 Optimization UI + Trajectory Animation
 */

import { useCallback, useRef } from 'react';
import { useWorker, createRequestId } from './useWorker';
import { useScfStore } from '../stores/scfStore';
import type {
  OptimizeGeometryRequest,
  OptimizationResult,
  OptimizationProgress,
  WorkerProgress,
  RequestId,
} from '../worker/protocol';
import type { DftMethod } from '../types/dft';

// ============================================================================
// Types
// ============================================================================

/**
 * Parameters for the useOptimizer hook.
 */
export interface UseOptimizerParams {
  /** Atom coordinates as [Z, x, y, z][] (null if not available) */
  atoms: [number, number, number, number][] | null;
  /** Basis set name */
  basisName: string;
  /** Electronic structure method */
  method: DftMethod;
  /** Whether the worker is ready */
  isReady: boolean;
}

/**
 * Return type of useOptimizer hook.
 */
export interface UseOptimizerResult {
  /** Start geometry optimization */
  optimize: () => void;
  /** Cancel running optimization */
  cancel: () => void;
  /** Whether optimization is running */
  loading: boolean;
  /** Error message (null if no error) */
  error: string | null;
}

/**
 * Type guard to check if a progress update is for optimization.
 */
function isOptimizationProgress(
  progress: WorkerProgress
): progress is OptimizationProgress {
  return progress.module === 'optimization';
}

/**
 * Map DftMethod to the method string for optimize_geometry.
 */
function toOptMethod(method: DftMethod): string {
  switch (method) {
    case 'rhf':
      return 'rhf';
    case 'lda':
      return 'lda';
    case 'b3lyp':
      return 'b3lyp';
    case 'b3lyp-d3bj':
      return 'b3lyp'; // D3(BJ) dispersion not yet supported in optimizer
  }
}

// ============================================================================
// Hook Implementation
// ============================================================================

/**
 * React hook for geometry optimization with progress streaming.
 *
 * Unlike useScf which runs a single-point calculation, this hook runs
 * the full L-BFGS optimization loop. Progress is streamed after each
 * optimization step.
 *
 * @param params - Atoms, basis, method, and worker readiness
 * @returns Object with optimize, cancel, loading, and error
 */
export function useOptimizer({
  atoms,
  basisName,
  method,
  isReady,
}: UseOptimizerParams): UseOptimizerResult {
  const { send, cancel: cancelRequest, isReady: workerReady } = useWorker();

  // Store selectors
  const running = useScfStore((state) => state.optimizationState.running);
  const optimizationRequestId = useScfStore(
    (state) => state.optimizationState.requestId
  );

  // Store actions
  const startOptimization = useScfStore((state) => state.startOptimization);
  const addOptimizationStep = useScfStore((state) => state.addOptimizationStep);
  const completeOptimization = useScfStore((state) => state.completeOptimization);
  const failOptimization = useScfStore((state) => state.failOptimization);
  const cancelOptimization = useScfStore((state) => state.cancelOptimization);

  // Ref to track current request ID for cleanup
  const requestIdRef = useRef<RequestId | null>(null);

  /**
   * Handle progress updates from the worker.
   */
  const handleProgress = useCallback(
    (progress: WorkerProgress) => {
      if (!isOptimizationProgress(progress)) return;

      // We receive progress with step/energy/gradient, but the full
      // OptimizationStepResult (with geometry/gradient arrays) comes
      // in the final result. For the progress plot, we create a
      // partial step result.
      addOptimizationStep({
        step: progress.step,
        energy: progress.energy,
        maxGradient: progress.maxGradient,
        rmsGradient: progress.rmsGradient,
        geometry: [],  // Not available from progress callback
        gradient: [],  // Not available from progress callback
      });
    },
    [addOptimizationStep]
  );

  /**
   * Start geometry optimization.
   */
  const optimize = useCallback(async () => {
    if (!isReady || !workerReady || running) return;
    if (!atoms || atoms.length === 0) return;

    const requestId = createRequestId();
    requestIdRef.current = requestId;
    startOptimization(requestId);

    try {
      const request = {
        type: 'optimize_geometry' as const,
        atoms,
        basisName,
        method: toOptMethod(method),
      } as Omit<OptimizeGeometryRequest, 'requestId'>;

      const result = await send<OptimizationResult>(request, {
        onProgress: handleProgress,
      });

      completeOptimization(result);
    } catch (err) {
      const message =
        err instanceof Error ? err.message : 'Geometry optimization failed';
      if (
        err instanceof Error &&
        (err.message.includes('cancelled') ||
          err.message.includes('COMPUTATION_CANCELLED'))
      ) {
        cancelOptimization();
      } else {
        failOptimization(message);
      }
    } finally {
      requestIdRef.current = null;
    }
  }, [
    isReady,
    workerReady,
    running,
    atoms,
    basisName,
    method,
    send,
    handleProgress,
    startOptimization,
    completeOptimization,
    failOptimization,
    cancelOptimization,
  ]);

  /**
   * Cancel running optimization.
   */
  const cancel = useCallback(() => {
    if (optimizationRequestId) {
      cancelRequest(optimizationRequestId as RequestId);
      cancelOptimization();
    }
  }, [optimizationRequestId, cancelRequest, cancelOptimization]);

  return {
    optimize,
    cancel,
    loading: running,
    error: useScfStore.getState().optimizationState.error,
  };
}
