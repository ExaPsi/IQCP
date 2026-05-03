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
 * US-070: Routes computation to scf_run (RHF) or ks_scf (DFT)
 * based on the selected method. The DFT path requires atom coordinates
 * and basis name, passed via RunOptions.
 *
 * @module hooks/useScf
 */

import { useCallback, useRef } from 'react';
import { useWorker, createRequestId } from './useWorker';
import { useScfStore } from '../stores/scfStore';
import { isDftMethod, toKsMethod } from '../types/dft';
import type {
  ScfRunResult,
  KsScfResult,
  KsScfRequest,
  ScfIterationHistory,
  ScfIterationProgress,
  ScfIntegralProgress,
  ScfRunRequest,
  PopulationAnalysisResult,
  WorkerProgress,
  RequestId,
} from '../worker/protocol';

/**
 * Options for the run() method.
 *
 * When using DFT methods, atoms and basisName must be provided
 * since the ks_scf WASM path requires atom coordinates.
 *
 * @see US-070 DFT UI + Method Selector + Deep Links
 */
export interface RunOptions {
  /** Atom coordinates as [Z, x, y, z] arrays (required for DFT methods) */
  atoms?: [number, number, number, number][];
  /** Basis set name (required for DFT methods) */
  basisName?: string;
}

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
  run: (options?: RunOptions) => void;
  /** Cancel the running computation */
  cancelRun: () => void;
}

/**
 * Type guard to check if a progress update is for SCF iterations.
 */
function isScfProgress(progress: WorkerProgress): progress is ScfIterationProgress {
  return progress.module === 'scf';
}

/**
 * Type guard to check if a progress update is for pre-SCF integral computation.
 */
function isScfIntegralProgress(progress: WorkerProgress): progress is ScfIntegralProgress {
  return progress.module === 'scf_integrals';
}

/**
 * Atom-to-basis-function mapping for preset STO-3G systems.
 *
 * Each entry is { atomicNumber, nBasis } where nBasis is the number
 * of basis functions centered on that atom in STO-3G:
 * - H: 1 (1s)
 * - He: 1 (1s)
 * - Li: 2 (1s, 2s)
 * - N: 5 (1s, 2s, 2px, 2py, 2pz)
 * - O: 5 (1s, 2s, 2px, 2py, 2pz)
 *
 * @see US-076 Population Analysis
 */
/**
 * STO-3G basis function count per atom by atomic number.
 * Used to construct atom-to-basis mapping for population analysis.
 */
const STO3G_NBASIS: Record<number, number> = {
  1: 1,  // H: 1s
  2: 1,  // He: 1s
  3: 2,  // Li: 1s, 2s
  4: 2,  // Be: 1s, 2s
  5: 5,  // B: 1s, 2s, 2p
  6: 5,  // C: 1s, 2s, 2p
  7: 5,  // N: 1s, 2s, 2p
  8: 5,  // O: 1s, 2s, 2p
  9: 5,  // F: 1s, 2s, 2p
  10: 5, // Ne: 1s, 2s, 2p
};

/**
 * Extract atom-to-basis mapping from geometry atoms (for custom geometry systems).
 * Assumes STO-3G basis for now (can be extended for other basis sets).
 */
function getAtomBasisFromGeometry(
  atoms: { atomicNumber: number }[]
): { atomicNumber: number; nBasis: number }[] | undefined {
  const result: { atomicNumber: number; nBasis: number }[] = [];
  for (const atom of atoms) {
    const nBasis = STO3G_NBASIS[atom.atomicNumber];
    if (nBasis === undefined) return undefined; // Unknown element
    result.push({ atomicNumber: atom.atomicNumber, nBasis });
  }
  return result;
}

/**
 * Basis function counts per atom for supported basis sets.
 * Key: "Z:basis" for cartesian, "Z:basis:sph" for spherical.
 * Values cross-validated against PySCF 2.11.0.
 * Spherical entries only needed for basis sets with d-shells (6-31G*, 6-31+G*).
 */
const BASIS_NBASIS: Record<string, number> = {
  // STO-3G: minimal basis (no d-shells, cart == sph)
  '1:sto-3g': 1, '2:sto-3g': 1, '3:sto-3g': 5, '4:sto-3g': 5, '5:sto-3g': 5,
  '6:sto-3g': 5, '7:sto-3g': 5, '8:sto-3g': 5, '9:sto-3g': 5, '10:sto-3g': 5,
  '11:sto-3g': 9, '12:sto-3g': 9, '13:sto-3g': 9, '14:sto-3g': 9, '15:sto-3g': 9,
  '16:sto-3g': 9, '17:sto-3g': 9, '18:sto-3g': 9,
  // 3-21G (no d-shells, cart == sph)
  '1:3-21g': 2, '2:3-21g': 2, '3:3-21g': 9, '4:3-21g': 9, '5:3-21g': 9,
  '6:3-21g': 9, '7:3-21g': 9, '8:3-21g': 9, '9:3-21g': 9, '10:3-21g': 9,
  '11:3-21g': 13, '12:3-21g': 13, '13:3-21g': 13, '14:3-21g': 13, '15:3-21g': 13,
  '16:3-21g': 13, '17:3-21g': 13, '18:3-21g': 13,
  // 6-31G (no d-shells, cart == sph)
  '1:6-31g': 2, '2:6-31g': 2, '3:6-31g': 9, '4:6-31g': 9, '5:6-31g': 9,
  '6:6-31g': 9, '7:6-31g': 9, '8:6-31g': 9, '9:6-31g': 9, '10:6-31g': 9,
  '11:6-31g': 13, '12:6-31g': 13, '13:6-31g': 13, '14:6-31g': 13, '15:6-31g': 13,
  '16:6-31g': 13, '17:6-31g': 13, '18:6-31g': 13,
  // 6-31G* Cartesian (d = 6 functions): all Li-Ar get d-shell
  '1:6-31g*': 2, '2:6-31g*': 2, '3:6-31g*': 15, '4:6-31g*': 15, '5:6-31g*': 15,
  '6:6-31g*': 15, '7:6-31g*': 15, '8:6-31g*': 15, '9:6-31g*': 15, '10:6-31g*': 15,
  '11:6-31g*': 19, '12:6-31g*': 19, '13:6-31g*': 19, '14:6-31g*': 19, '15:6-31g*': 19,
  '16:6-31g*': 19, '17:6-31g*': 19, '18:6-31g*': 19,
  // 6-31G* Spherical (d = 5 functions): each d-shell → 1 fewer BF
  '1:6-31g*:sph': 2, '2:6-31g*:sph': 2, '3:6-31g*:sph': 14, '4:6-31g*:sph': 14, '5:6-31g*:sph': 14,
  '6:6-31g*:sph': 14, '7:6-31g*:sph': 14, '8:6-31g*:sph': 14, '9:6-31g*:sph': 14, '10:6-31g*:sph': 14,
  '11:6-31g*:sph': 18, '12:6-31g*:sph': 18, '13:6-31g*:sph': 18, '14:6-31g*:sph': 18, '15:6-31g*:sph': 18,
  '16:6-31g*:sph': 18, '17:6-31g*:sph': 18, '18:6-31g*:sph': 18,
  // 6-31+G* Cartesian (diffuse sp + d polarization on Li-Ar)
  '1:6-31+g*': 2, '2:6-31+g*': 2, '3:6-31+g*': 19, '4:6-31+g*': 19, '5:6-31+g*': 19,
  '6:6-31+g*': 19, '7:6-31+g*': 19, '8:6-31+g*': 19, '9:6-31+g*': 19, '10:6-31+g*': 19,
  // 6-31+G* Spherical (each d-shell → 1 fewer)
  '1:6-31+g*:sph': 2, '2:6-31+g*:sph': 2, '3:6-31+g*:sph': 18, '4:6-31+g*:sph': 18, '5:6-31+g*:sph': 18,
  '6:6-31+g*:sph': 18, '7:6-31+g*:sph': 18, '8:6-31+g*:sph': 18, '9:6-31+g*:sph': 18, '10:6-31+g*:sph': 18,
};

/**
 * Get atom-to-basis mapping from atom specs and basis set name.
 * Works for any basis set (not just STO-3G).
 * Tries spherical key first (e.g., "8:6-31g*:sph"), falls back to cartesian.
 */
function getAtomBasisFromAtomSpecs(
  atoms: [number, number, number, number][],
  basisName: string,
  totalNbf: number,
): { atomicNumber: number; nBasis: number }[] | undefined {
  // Try spherical keys first, then cartesian
  for (const suffix of [':sph', '']) {
    const result: { atomicNumber: number; nBasis: number }[] = [];
    let sum = 0;
    let allFound = true;
    for (const [z] of atoms) {
      const key = `${z}:${basisName}${suffix}`;
      const nb = BASIS_NBASIS[key];
      if (nb === undefined) { allFound = false; break; }
      result.push({ atomicNumber: z, nBasis: nb });
      sum += nb;
    }
    if (allFound && sum === totalNbf) return result;
  }
  return undefined;
}

const PRESET_ATOM_BASIS: Record<string, { atomicNumber: number; nBasis: number }[]> = {
  'h2_sto3g_r1.4': [
    { atomicNumber: 1, nBasis: 1 },
    { atomicNumber: 1, nBasis: 1 },
  ],
  'heh_plus_sto3g': [
    { atomicNumber: 2, nBasis: 1 },
    { atomicNumber: 1, nBasis: 1 },
  ],
  'lih_sto3g': [
    { atomicNumber: 3, nBasis: 2 },
    { atomicNumber: 1, nBasis: 1 },
  ],
  'h2o_sto3g': [
    { atomicNumber: 8, nBasis: 5 },
    { atomicNumber: 1, nBasis: 1 },
    { atomicNumber: 1, nBasis: 1 },
  ],
  'nh3_sto3g': [
    { atomicNumber: 7, nBasis: 5 },
    { atomicNumber: 1, nBasis: 1 },
    { atomicNumber: 1, nBasis: 1 },
    { atomicNumber: 1, nBasis: 1 },
  ],
};

/**
 * Adapt a KsScfResult to the ScfRunResult shape for store compatibility.
 *
 * The existing result display components (ScfResultDisplay, ScfEnergyPlot,
 * ScfResidualPlot, ScfIterationTable) all consume ScfRunResult. By adapting
 * the KS result, we avoid duplicating these components.
 *
 * Limitation: DFT results do not include matrices or orbital energies
 * in the current WASM export. These fields are left undefined.
 */
function adaptKsResult(ksResult: KsScfResult): ScfRunResult {
  // The WASM trace entries have camelCase field names from Rust serde:
  //   { iteration, energy, deltaE, rmsDensityChange, diisApplied }
  // But ScfIterationHistory expects:
  //   { iteration, energy, delta, diisError? }
  // The RHF handler (scf.ts) does this mapping; we must do the same here.
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const history: ScfIterationHistory[] = (ksResult.trace as any[]).map((iter) => {
    const entry: ScfIterationHistory = {
      iteration: iter.iteration ?? 0,
      energy: iter.energy ?? 0,
      // WASM serde camelCase: deltaE; ScfIterationHistory expects: delta
      delta: iter.deltaE ?? iter.delta ?? 0,
    };
    // WASM serde camelCase: rmsDensityChange; ScfIterationHistory expects: diisError
    const rms = iter.rmsDensityChange ?? iter.diisError;
    if (rms !== null && rms !== undefined) {
      entry.diisError = rms;
    }
    return entry;
  });

  return {
    energy: ksResult.energy,
    converged: ksResult.converged,
    iterations: ksResult.iterations,
    aborted: false,
    history,
    matrices: ksResult.densityMatrix ? {
      nbf: ksResult.nBasis,
      sMatrix: ksResult.overlapMatrix ?? [],
      hCore: ksResult.hCore ?? [],
      fockMatrix: ksResult.fockMatrix ?? [],
      densityMatrix: ksResult.densityMatrix,
      moCoefficients: ksResult.moCoefficients,
    } : undefined,
    orbitalEnergies: ksResult.orbitalEnergies ? {
      energies: ksResult.orbitalEnergies,
      nOccupied: ksResult.nOccupied,
    } : undefined,
  };
}

/**
 * React hook for SCF computation with progress streaming.
 *
 * Unlike useBoys/useRys, this does NOT auto-compute on parameter change.
 * Instead, computation is triggered by the explicit run() method.
 * Progress updates are streamed via the worker's onProgress callback
 * and stored in the scfStore for live UI updates.
 *
 * US-070: Routes to scf_run for RHF or ks_scf for DFT methods.
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
 *       <button onClick={() => run({ atoms: workerAtoms, basisName })} disabled={!isReady || isRunning}>
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
  const useSpherical = useScfStore((state) => state.useSpherical);
  const gridQuality = useScfStore((state) => state.gridQuality);
  const levelShift = useScfStore((state) => state.levelShift);
  const method = useScfStore((state) => state.method);
  const computeStatus = useScfStore((state) => state.compute.status);
  const runningRequestId = useScfStore((state) => state.runningRequestId);

  // Store actions
  const startRun = useScfStore((state) => state.startRun);
  const addIteration = useScfStore((state) => state.addIteration);
  const updateScfIntegralProgress = useScfStore((state) => state.updateScfIntegralProgress);
  const setRunResult = useScfStore((state) => state.setRunResult);
  const setRunError = useScfStore((state) => state.setRunError);
  const setRunCancelled = useScfStore((state) => state.setRunCancelled);
  const setKsResult = useScfStore((state) => state.setKsResult);
  const setPopulationResult = useScfStore((state) => state.setPopulationResult);

  // Ref to track current request ID for cleanup
  const requestIdRef = useRef<RequestId | null>(null);

  /**
   * Handle progress updates from the worker.
   *
   * Handles both integral-phase progress (S, Hcore, ERI, grid)
   * and SCF iteration progress.
   */
  const handleProgress = useCallback(
    (progress: WorkerProgress) => {
      if (isScfIntegralProgress(progress)) {
        updateScfIntegralProgress(progress.step, progress.percent);
        return;
      }
      if (!isScfProgress(progress)) return;

      addIteration({
        iteration: progress.iteration,
        energy: progress.energy,
        delta: progress.delta,
        diisError: progress.diisError,
      });
    },
    [addIteration, updateScfIntegralProgress]
  );

  /**
   * Start an SCF computation.
   *
   * Routes to scf_run (RHF) or ks_scf (DFT) based on the selected method.
   * For DFT methods, atoms and basisName must be provided in options.
   *
   * Does nothing if worker is not ready or computation is already running.
   */
  const run = useCallback(async (options?: RunOptions) => {
    if (!isReady || computeStatus === 'running') return;

    const requestId = createRequestId();
    requestIdRef.current = requestId;
    startRun(requestId);

    // Clear previous results
    setKsResult(null);
    setPopulationResult(null);

    try {
      // Use ks_scf path when: (1) DFT method selected, OR (2) RHF with on-the-fly atoms
      const hasOnTheFlyAtoms = options?.atoms && options.atoms.length > 0 && options?.basisName;
      if (isDftMethod(method) || hasOnTheFlyAtoms) {
        // === On-the-fly path: use ks_scf (supports RHF, LDA, B3LYP, B3LYP-D3BJ) ===
        const atoms = options?.atoms;
        const basisName = options?.basisName;

        if (!atoms || atoms.length === 0) {
          throw new Error('Atom geometry is required for on-the-fly calculations');
        }
        if (!basisName) {
          throw new Error('Basis set name is required for on-the-fly calculations');
        }

        // Type assertion needed because TypeScript can't narrow Omit<UnionType, 'field'>
        const ksRequest = {
          type: 'ks_scf' as const,
          atoms,
          basisName,
          method: toKsMethod(method),
          convergenceProfile,
          maxIterations,
          useDiis,
          useSpherical,
          gridQuality,
        } as Omit<KsScfRequest, 'requestId'>;
        const ksResult = await send<KsScfResult>(
          ksRequest,
          { onProgress: handleProgress }
        );

        // Store raw KS result for DftInfoPanel energy decomposition
        setKsResult(ksResult);

        // Adapt to ScfRunResult shape for existing display components
        const adapted = adaptKsResult(ksResult);
        setRunResult(adapted);

        // Population analysis for DFT (US-076)
        // KS result now includes overlap matrix S directly
        if (ksResult.densityMatrix && ksResult.densityMatrix.length > 0 &&
            ksResult.overlapMatrix && ksResult.overlapMatrix.length > 0) {
          // Get atom-basis mapping from the molecule preset or custom geometry
          let atomBasis: { atomicNumber: number; nBasis: number }[] | undefined;
          const mol = options?.atoms;
          if (mol && mol.length > 0) {
            // Build from the atoms we sent to KS-SCF
            const basisName = options?.basisName ?? 'sto-3g';
            atomBasis = getAtomBasisFromAtomSpecs(mol, basisName, ksResult.nBasis);
          }
          if (!atomBasis) {
            // Fallback: try preset lookup
            atomBasis = PRESET_ATOM_BASIS[systemId];
          }
          if (atomBasis) {
            try {
              const popReq: Omit<import('../worker/protocol').PopulationAnalysisRequest, 'requestId'> = {
                type: 'population_analysis',
                densityMatrix: ksResult.densityMatrix,
                overlapMatrix: ksResult.overlapMatrix,
                nbf: ksResult.nBasis,
                atoms: atomBasis,
              };
              const popResult = await send<PopulationAnalysisResult>(popReq);
              setPopulationResult(popResult);
            } catch (popErr) {
              console.warn('[useScf] DFT population analysis failed:', popErr);
            }
          }
        }
      } else {
        // === RHF path: use scf_run (existing behavior) ===
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
            levelShift: levelShift > 0 ? levelShift : undefined,
            includeMatrices: true,
          },
        };
        const result = await send<ScfRunResult>(request, { onProgress: handleProgress });

        // Check if computation was aborted
        if (result.aborted) {
          setRunCancelled();
        } else {
          setRunResult(result);

          // Trigger population analysis if matrices are available (US-076)
          if (result.matrices && result.matrices.sMatrix.length > 0) {
            // Try preset lookup first, then custom geometry
            let atomBasis: { atomicNumber: number; nBasis: number }[] | undefined = PRESET_ATOM_BASIS[systemId];
            if (!atomBasis) {
              // Custom geometry: extract atom info from integral compute result
              const integralState = useScfStore.getState().integralCompute;
              if (integralState.status === 'success' && integralState.systemData.geometry) {
                atomBasis = getAtomBasisFromGeometry(integralState.systemData.geometry.atoms);
              }
            }
            if (atomBasis) {
              try {
                const popReq: Omit<import('../worker/protocol').PopulationAnalysisRequest, 'requestId'> = {
                  type: 'population_analysis',
                  densityMatrix: result.matrices.densityMatrix,
                  overlapMatrix: result.matrices.sMatrix,
                  nbf: result.matrices.nbf,
                  atoms: atomBasis,
                };
                const popResult = await send<PopulationAnalysisResult>(popReq);
                setPopulationResult(popResult);
              } catch (popErr) {
                // Population analysis is non-critical; don't fail the SCF
                console.warn('[useScf] Population analysis failed:', popErr);
              }
            }
          }
        }
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
    method,
    systemId,
    convergenceProfile,
    maxIterations,
    useDiis,
    useSpherical,
    gridQuality,
    levelShift,
    damp,
    send,
    handleProgress,
    startRun,
    setRunResult,
    setRunCancelled,
    setRunError,
    setKsResult,
    setPopulationResult,
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
