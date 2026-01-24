/**
 * SCF (Self-Consistent Field) handlers.
 *
 * Handles SCF computation requests by loading preset systems and
 * calling the WASM SCF function. Supports progress streaming and
 * cancellation.
 *
 * @module worker/handlers/scf
 */

import type {
  ScfRunRequest,
  ScfRunResult,
  WorkerProgress,
  ErrorResponse,
  ScfIterationHistory,
} from '../protocol';
import { getPreset, listAllSystemIds, getCustomSystemCount } from '../presets';

/**
 * Progress callback type for SCF iteration updates.
 */
export type ProgressCallback = (progress: WorkerProgress) => void;

/**
 * Abort flag checker type.
 */
export type AbortChecker = () => boolean;

/**
 * WASM SCF iteration structure (matches ScfWasmIteration in qc-wasm).
 *
 * Note: Field names use camelCase as they come from serde serialization.
 */
interface ScfWasmIteration {
  /** Iteration number (0-indexed) */
  iteration: number;
  /** Total energy (electronic + nuclear) in Hartree */
  energy: number;
  /** Energy change from previous iteration */
  deltaE: number | null;
  /** RMS change in density matrix */
  rmsDensityChange: number | null;
  /** Whether DIIS was applied in this iteration */
  diisApplied: boolean;
}

/**
 * WASM SCF matrices structure (matches ScfWasmMatrices in qc-wasm).
 *
 * Note: Field names use camelCase as they come from serde serialization.
 */
interface ScfWasmMatrices {
  /** Number of basis functions (for reshaping) */
  nbf: number;
  /** Overlap matrix S (row-major, nbf x nbf) */
  sMatrix: number[];
  /** Core Hamiltonian H = T + V (row-major, nbf x nbf) */
  hCore: number[];
  /** Final Fock matrix F (row-major, nbf x nbf) */
  fockMatrix: number[];
  /** Final density matrix D (row-major, nbf x nbf) */
  densityMatrix: number[];
}

/**
 * WASM orbital energies structure (matches ScfWasmOrbitalEnergies in qc-wasm).
 *
 * Note: Field names use camelCase as they come from serde serialization.
 */
interface ScfWasmOrbitalEnergies {
  /** Orbital energies in Hartree, sorted ascending */
  energies: number[];
  /** Number of occupied orbitals (n_occ = nelec / 2 for RHF) */
  nOccupied: number;
}

/**
 * WASM SCF result structure (matches ScfWasmResult in qc-wasm).
 *
 * Note: Field names use camelCase as they come from serde serialization.
 */
interface ScfWasmResult {
  /** Whether SCF converged */
  converged: boolean;
  /** Final total energy in Hartree */
  energy: number;
  /** Number of iterations performed */
  iterations: number;
  /** Iteration-by-iteration trace */
  trace: ScfWasmIteration[];
  /** Optional matrices (when includeMatrices is true) */
  matrices?: ScfWasmMatrices;
  /** Optional orbital energies (when includeMatrices is true) */
  orbitalEnergies?: ScfWasmOrbitalEnergies;
}

/**
 * Type for the WASM scf_run function.
 *
 * This is the function signature exported from qc-wasm.
 */
export type ScfRunFn = (systemJson: string, options: unknown) => ScfWasmResult;

// Store reference to the WASM scf_run function (set during worker initialization)
let wasmScfRun: ScfRunFn | null = null;

/**
 * Set the WASM scf_run function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function is stored and used by handleScfRun.
 *
 * @param fn - The WASM scf_run function
 */
export function setScfRunFn(fn: ScfRunFn): void {
  wasmScfRun = fn;
}

/**
 * Handle scf_run request.
 *
 * Loads a preset molecular system and runs an RHF SCF calculation.
 * Progress is emitted for each iteration, and the computation can be
 * cancelled via the abort checker.
 *
 * @param request - The SCF run request
 * @param onProgress - Callback for progress updates
 * @param isAborted - Function to check if computation should abort
 * @returns SCF result or error response
 *
 * @example
 * ```typescript
 * const result = handleScfRun(
 *   { type: 'scf_run', requestId: id, systemId: 'h2_sto3g_r1.4', options },
 *   (progress) => postMessage({ type: 'progress', requestId: id, progress }),
 *   () => abortFlags.get(id) === true
 * );
 * ```
 */
export function handleScfRun(
  request: ScfRunRequest,
  onProgress: ProgressCallback,
  isAborted: AbortChecker
): ScfRunResult | ErrorResponse {
  const { requestId, systemId, options } = request;

  // Check if WASM function is available
  if (!wasmScfRun) {
    return {
      type: 'error',
      requestId,
      code: 'WORKER_NOT_READY',
      message: 'SCF WASM function not initialized',
    };
  }

  // 1. Load preset system (or custom computed system)
  const preset = getPreset(systemId);

  if (!preset) {
    const allSystems = listAllSystemIds();
    const customCount = getCustomSystemCount();
    const systemList =
      customCount > 0
        ? `${allSystems.slice(0, 10).join(', ')}${allSystems.length > 10 ? `, ... (${allSystems.length} total, ${customCount} custom)` : ''}`
        : 'h2_sto3g_r1.4, heh_plus_sto3g, lih_sto3g, h2o_sto3g, nh3_sto3g';
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `Unknown system: ${systemId}. Available systems: ${systemList}`,
    };
  }

  // Check abort before starting
  if (isAborted()) {
    return {
      energy: 0,
      converged: false,
      iterations: 0,
      aborted: true,
      history: [],
    };
  }

  // 2. Prepare WASM options (camelCase for JavaScript/WASM)
  const wasmOptions = {
    convergenceProfile: options.convergenceProfile,
    maxIterations: options.maxIterations,
    useDiis: options.useDiis,
    diisSize: options.diisSize,
    includeMatrices: options.includeMatrices ?? false,
    damp: options.damp,
  };

  // 3. Call WASM SCF function
  let result: ScfWasmResult;
  try {
    // Pass the preset JSON as a string
    result = wasmScfRun(JSON.stringify(preset), wasmOptions);
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'SCF computation failed',
    };
  }

  // 4. Process trace and emit progress for each iteration
  const history: ScfIterationHistory[] = [];

  for (let i = 0; i < result.trace.length; i++) {
    const iter = result.trace[i];

    // Check abort between iterations
    if (isAborted()) {
      // Return partial results up to this point
      return {
        energy: iter.energy,
        converged: false,
        iterations: i + 1,
        aborted: true,
        history,
      };
    }

    // Add to history
    const historyEntry: ScfIterationHistory = {
      iteration: iter.iteration,
      energy: iter.energy,
      delta: iter.deltaE ?? 0,
    };

    // Include DIIS error if available
    if (iter.rmsDensityChange !== null && iter.rmsDensityChange !== undefined) {
      historyEntry.diisError = iter.rmsDensityChange;
    }

    history.push(historyEntry);

    // Determine if this is the final converged iteration
    const isFinal = i === result.trace.length - 1;
    const iterConverged = isFinal && result.converged;

    // Emit progress
    onProgress({
      module: 'scf',
      iteration: iter.iteration,
      energy: iter.energy,
      delta: iter.deltaE ?? 0,
      diisError: iter.rmsDensityChange ?? undefined,
      converged: iterConverged,
      current: i + 1,
      total: result.iterations,
      message: iterConverged
        ? `Converged after ${i + 1} iterations`
        : `Iteration ${i + 1}/${result.iterations}: E = ${iter.energy.toFixed(10)} Ha`,
    });
  }

  // 5. Build final result
  const finalResult: ScfRunResult = {
    energy: result.energy,
    converged: result.converged,
    iterations: result.iterations,
    aborted: false,
    history,
  };

  // 6. Include matrices and orbital energies if available
  if (result.matrices) {
    finalResult.matrices = {
      nbf: result.matrices.nbf,
      sMatrix: result.matrices.sMatrix,
      hCore: result.matrices.hCore,
      fockMatrix: result.matrices.fockMatrix,
      densityMatrix: result.matrices.densityMatrix,
    };
  }

  if (result.orbitalEnergies) {
    finalResult.orbitalEnergies = {
      energies: result.orbitalEnergies.energies,
      nOccupied: result.orbitalEnergies.nOccupied,
    };
  }

  return finalResult;
}
