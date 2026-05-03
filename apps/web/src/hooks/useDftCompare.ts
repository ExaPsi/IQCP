/**
 * useDftCompare - Hook for DFT vs HF comparison mode.
 *
 * Orchestrates two sequential computations:
 * 1. RHF via scf_run (using preset system data or custom integrals)
 * 2. Selected DFT method via ks_scf (using atom coordinates)
 *
 * Results are stored in scfStore.dftCompareState for display
 * in the DftComparisonPanel component.
 *
 * @module hooks/useDftCompare
 * @see US-070 DFT UI + Method Selector + Deep Links (AC4)
 */

import { useCallback } from 'react';
import { useWorker } from './useWorker';
import { useScfStore } from '../stores/scfStore';
import { isDftMethod, toKsMethod } from '../types/dft';
import type { KsScfResult, KsScfRequest } from '../worker/protocol';

/**
 * Parameters for the useDftCompare hook.
 */
export interface UseDftCompareParams {
  /** Worker atoms for DFT computation: [Z, x, y, z][] */
  atoms: [number, number, number, number][] | null;
  /** Basis set name */
  basisName: string;
  /** Whether the worker is ready */
  isReady: boolean;
}

/**
 * Return type of useDftCompare hook.
 */
export interface UseDftCompareResult {
  /** Start a comparison run */
  compare: () => Promise<void>;
  /** Whether comparison is in progress */
  isComparing: boolean;
}

/**
 * React hook for DFT vs HF comparison mode.
 *
 * Runs two sequential computations (RHF then DFT) for the same molecule
 * and stores the results in the scfStore for side-by-side display.
 *
 * @param params - Atoms, basis name, and worker readiness
 * @returns Object with compare function and isComparing flag
 *
 * @example
 * ```typescript
 * const { compare, isComparing } = useDftCompare({
 *   atoms: workerAtoms,
 *   basisName: 'sto-3g',
 *   isReady: true,
 * });
 *
 * <button onClick={compare} disabled={isComparing}>Compare</button>
 * ```
 */
export function useDftCompare({ atoms, basisName, isReady }: UseDftCompareParams): UseDftCompareResult {
  const { send } = useWorker();

  // Store selectors
  const method = useScfStore((s) => s.method);
  const convergenceProfile = useScfStore((s) => s.convergenceProfile);
  const maxIterations = useScfStore((s) => s.maxIterations);
  const useDiis = useScfStore((s) => s.useDiis);
  const isComparing = useScfStore((s) => s.dftCompareState.comparing);

  // Store actions
  const startDftCompare = useScfStore((s) => s.startDftCompare);
  const setDftCompareRhfResult = useScfStore((s) => s.setDftCompareRhfResult);
  const setDftCompareDftResult = useScfStore((s) => s.setDftCompareDftResult);
  const completeDftCompare = useScfStore((s) => s.completeDftCompare);
  const failDftCompare = useScfStore((s) => s.failDftCompare);

  const compare = useCallback(async () => {
    if (!isReady || isComparing || !isDftMethod(method)) return;
    if (!atoms || atoms.length === 0) {
      failDftCompare('No atom geometry available for DFT comparison');
      return;
    }

    startDftCompare();

    try {
      // Phase 1: Run RHF via ks_scf (method='rhf') using the actual molecule geometry.
      // Previously used scf_run with systemId, which broke for non-preset molecules
      // (e.g., C6H6 would show H2 energy because systemId was stale).
      const rhfRequest = {
        type: 'ks_scf' as const,
        atoms,
        basisName,
        method: 'rhf',
        convergenceProfile,
        maxIterations,
        useDiis,
      } as Omit<KsScfRequest, 'requestId'>;
      const rhfResult = await send<KsScfResult>(rhfRequest);

      setDftCompareRhfResult({
        energy: rhfResult.energy,
        iterations: rhfResult.iterations,
        converged: rhfResult.converged,
      });

      // Phase 2: Run DFT via ks_scf
      const ksRequest = {
        type: 'ks_scf' as const,
        atoms,
        basisName,
        method: toKsMethod(method),
        convergenceProfile,
        maxIterations,
        useDiis,
      } as Omit<KsScfRequest, 'requestId'>;
      const ksResult = await send<KsScfResult>(ksRequest);

      setDftCompareDftResult({
        energy: ksResult.energy,
        iterations: ksResult.iterations,
        converged: ksResult.converged,
        method: ksResult.method,
      });

      completeDftCompare();
    } catch (err) {
      const message = err instanceof Error ? err.message : 'Comparison failed';
      failDftCompare(message);
    }
  }, [
    isReady, isComparing, method, atoms, basisName,
    convergenceProfile, maxIterations, useDiis,
    send, startDftCompare, setDftCompareRhfResult, setDftCompareDftResult,
    completeDftCompare, failDftCompare,
  ]);

  return { compare, isComparing };
}
