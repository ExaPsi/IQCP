/**
 * useBasisComparison - React hook for multi-basis SCF comparison.
 *
 * Orchestrates sequential integral_compute + scf_run for each selected
 * basis set. Results accumulate in the scfStore's basisCompareState.
 *
 * The hook runs independently from the main SCF compute pipeline --
 * it does not modify the main compute/history state.
 *
 * @module hooks/useBasisComparison
 * @see US-045 Basis Set Comparison Mode
 */

import { useCallback, useRef } from 'react';
import { useWorker } from './useWorker';
import { useScfStore } from '../stores/scfStore';
import type { BasisComparisonResult } from '../stores/scfStore';
import type {
  GeometryInput,
  IntegralComputeRequest,
  IntegralComputeResult,
  ScfRunRequest,
  ScfRunResult,
  ConvergenceProfile,
  BasisSetName,
} from '../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * SCF options for comparison runs.
 */
export interface BasisCompareOptions {
  /** Convergence profile (loose/medium/tight) */
  convergenceProfile: ConvergenceProfile;
  /** Maximum number of SCF iterations */
  maxIterations: number;
  /** Enable DIIS acceleration */
  useDiis: boolean;
  /** Fock matrix damping factor */
  damp: number;
  /** Use spherical harmonic basis functions */
  useSpherical: boolean;
}

/**
 * Return type of useBasisComparison hook.
 */
export interface UseBasisComparisonReturn {
  /** Whether the worker is ready */
  isReady: boolean;
  /** Whether a comparison run is in progress */
  isComparing: boolean;
  /** Start a comparison run for the given bases and geometry */
  compare: (bases: string[], geometry: GeometryInput) => Promise<void>;
  /** Cancel the running comparison */
  cancel: () => void;
  /** Clear all comparison results */
  clear: () => void;
}

// ============================================================================
// Hook Implementation
// ============================================================================

/**
 * React hook for multi-basis SCF comparison.
 *
 * For each basis set in the list, the hook:
 * 1. Sends an `integral_compute` request to compute integrals
 * 2. Sends an `scf_run` request using the computed system
 * 3. Stores the result in the scfStore's basisCompareState
 *
 * Progress updates are posted after each basis completes.
 * Supports cancellation via an abort flag.
 *
 * @returns Object with isReady, isComparing, compare, cancel, and clear
 *
 * @example
 * ```typescript
 * const { isReady, isComparing, compare, cancel, clear } = useBasisComparison();
 *
 * const handleCompare = async () => {
 *   await compare(['sto-3g', '6-31g', '6-31g*'], geometry);
 * };
 * ```
 */
export function useBasisComparison(): UseBasisComparisonReturn {
  const { send, isReady } = useWorker();

  // Store selectors
  const comparing = useScfStore((state) => state.basisCompareState.comparing);

  // Store actions
  const startBasisCompare = useScfStore((state) => state.startBasisCompare);
  const updateBasisCompareProgress = useScfStore((state) => state.updateBasisCompareProgress);
  const addBasisCompareResult = useScfStore((state) => state.addBasisCompareResult);
  const completeBasisCompare = useScfStore((state) => state.completeBasisCompare);
  const failBasisCompare = useScfStore((state) => state.failBasisCompare);
  const resetBasisCompare = useScfStore((state) => state.resetBasisCompare);

  // Abort flag for cancellation
  const abortRef = useRef(false);

  /**
   * Run the multi-basis comparison.
   *
   * @param bases - Array of basis set names to compare
   * @param geometry - Molecular geometry (same for all bases)
   */
  const compare = useCallback(
    async (bases: string[], geometry: GeometryInput) => {
      if (!isReady || comparing || bases.length < 2) return;

      abortRef.current = false;
      startBasisCompare();

      try {
        for (let i = 0; i < bases.length; i++) {
          if (abortRef.current) break;

          const basisName = bases[i];

          // Update progress
          updateBasisCompareProgress({
            current: i + 1,
            total: bases.length,
            currentBasis: basisName,
          });

          const startTime = performance.now();

          try {
            // Read SCF options from the store at the time of the run
            const storeState = useScfStore.getState();
            const options: BasisCompareOptions = {
              convergenceProfile: storeState.convergenceProfile,
              maxIterations: storeState.maxIterations,
              useDiis: storeState.useDiis,
              damp: storeState.damp,
              useSpherical: storeState.useSpherical,
            };

            // Step 1: Compute integrals for this basis
            // Type assertion needed because TypeScript can't narrow Omit<UnionType, 'field'>
            const integralRequest: Omit<IntegralComputeRequest, 'requestId'> = {
              type: 'integral_compute',
              geometry,
              basisSet: basisName as BasisSetName,
              useSpherical: options.useSpherical,
            };
            const integralResult = await send<IntegralComputeResult>(integralRequest);

            if (abortRef.current) break;

            // Step 2: Run SCF with the computed integrals
            const scfRequest: Omit<ScfRunRequest, 'requestId'> = {
              type: 'scf_run',
              systemId: integralResult.systemId,
              options: {
                convergenceProfile: options.convergenceProfile,
                maxIterations: options.maxIterations,
                useDiis: options.useDiis,
                damp: options.damp > 0 ? options.damp : undefined,
                includeMatrices: false,
              },
            };
            const scfResult = await send<ScfRunResult>(scfRequest);

            if (abortRef.current) break;

            const elapsedMs = performance.now() - startTime;

            // Store the result
            const result: BasisComparisonResult = {
              basisName,
              nbf: integralResult.nbf,
              energy: scfResult.energy,
              converged: scfResult.converged,
              iterations: scfResult.iterations,
              computeTimeMs: Math.round(elapsedMs),
              history: scfResult.history,
            };

            addBasisCompareResult(result);
          } catch (err) {
            if (abortRef.current) break;

            const elapsedMs = performance.now() - startTime;
            const errorMessage = err instanceof Error ? err.message : 'Unknown error';

            // Store failed result
            const result: BasisComparisonResult = {
              basisName,
              nbf: 0,
              energy: 0,
              converged: false,
              iterations: 0,
              computeTimeMs: Math.round(elapsedMs),
              history: [],
              error: errorMessage,
            };

            addBasisCompareResult(result);
          }
        }

        if (abortRef.current) {
          failBasisCompare('Comparison cancelled');
        } else {
          completeBasisCompare();
        }
      } catch (err) {
        const message = err instanceof Error ? err.message : 'Comparison failed';
        failBasisCompare(message);
      }
    },
    [
      isReady,
      comparing,
      send,
      startBasisCompare,
      updateBasisCompareProgress,
      addBasisCompareResult,
      completeBasisCompare,
      failBasisCompare,
    ]
  );

  /**
   * Cancel the running comparison.
   */
  const cancel = useCallback(() => {
    abortRef.current = true;
  }, []);

  /**
   * Clear all comparison results.
   */
  const clear = useCallback(() => {
    resetBasisCompare();
  }, [resetBasisCompare]);

  return {
    isReady,
    isComparing: comparing,
    compare,
    cancel,
    clear,
  };
}
