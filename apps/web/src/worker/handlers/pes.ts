/**
 * PES (Potential Energy Surface) scan handlers.
 *
 * Handles PES bond-length scan requests by calling the WASM pes_scan
 * function. Supports progress streaming and cancellation.
 *
 * @module worker/handlers/pes
 * @see US-039 PES WASM Export & Worker
 */

import type {
  PesScanRequest,
  PesScanResult,
  WorkerProgress,
  ErrorResponse,
} from '../protocol';

/**
 * Progress callback type for PES scan updates.
 */
export type ProgressCallback = (progress: WorkerProgress) => void;

/**
 * Abort flag checker type.
 */
export type AbortChecker = () => boolean;

/**
 * WASM PES scan progress structure (matches PesScanWasmProgress in qc-wasm).
 *
 * Note: Field names use camelCase as they come from serde serialization.
 */
interface PesScanWasmProgress {
  /** Index of the completed point (0-indexed) */
  pointIndex: number;
  /** Total number of scan points */
  totalPoints: number;
  /** Bond distance in bohr */
  r: number;
  /** SCF energy in Hartree */
  energy: number;
  /** Whether SCF converged at this geometry */
  converged: boolean;
}

/**
 * WASM PES point structure (matches PesPoint in qc-core).
 *
 * Note: Field names use snake_case as PesPoint in qc-core does NOT
 * use serde rename_all = "camelCase".
 */
interface PesWasmPoint {
  r: number;
  energy: number;
  converged: boolean;
  iterations: number;
}

/**
 * WASM PES equilibrium structure (matches PesEquilibrium in qc-core).
 *
 * Note: Field names use snake_case.
 */
interface PesWasmEquilibrium {
  r_bohr: number;
  energy_hartree: number;
}

/**
 * WASM PES scan result structure (matches PesScanResult in qc-core).
 *
 * Note: Field names use snake_case as PesScanResult in qc-core does NOT
 * use serde rename_all = "camelCase".
 */
interface PesWasmResult {
  points: PesWasmPoint[];
  equilibrium: PesWasmEquilibrium | null;
  compute_time_ms: number;
  total_iterations: number;
}

/**
 * Type for the WASM pes_scan function.
 *
 * This is the function signature exported from qc-wasm.
 */
export type PesScanFn = (
  input: unknown,
  progressCallback: (progress: PesScanWasmProgress) => void
) => PesWasmResult;

// Store reference to the WASM pes_scan function (set during worker initialization)
let wasmPesScan: PesScanFn | null = null;

/**
 * Set the WASM pes_scan function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function is stored and used by handlePesScan.
 *
 * @param fn - The WASM pes_scan function
 */
export function setPesScanFn(fn: PesScanFn): void {
  wasmPesScan = fn;
}

/**
 * Handle pes_scan request.
 *
 * Runs a PES bond-length scan by calling the WASM pes_scan function.
 * Progress is emitted for each scan point, and the computation can be
 * cancelled via the abort checker.
 *
 * @param request - The PES scan request
 * @param onProgress - Callback for progress updates
 * @param isAborted - Function to check if computation should abort
 * @returns PES scan result or error response
 *
 * @example
 * ```typescript
 * const result = handlePesScan(
 *   { type: 'pes_scan', requestId: id, atomAZ: 1, atomBZ: 1, ... },
 *   (progress) => postMessage({ type: 'progress', requestId: id, progress }),
 *   () => abortFlags.get(id) === true
 * );
 * ```
 */
export function handlePesScan(
  request: PesScanRequest,
  onProgress: ProgressCallback,
  isAborted: AbortChecker
): PesScanResult | ErrorResponse {
  const { requestId, atomAZ, atomBZ, rMin, rMax, nPoints, basisName, options, useSeeding } =
    request;

  // Check if WASM function is available
  if (!wasmPesScan) {
    return {
      type: 'error',
      requestId,
      code: 'WORKER_NOT_READY',
      message: 'PES scan WASM function not initialized',
    };
  }

  // Validate input parameters
  if (atomAZ < 1 || atomAZ > 10) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `atomAZ must be 1-10, got ${atomAZ}`,
    };
  }
  if (atomBZ < 1 || atomBZ > 10) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `atomBZ must be 1-10, got ${atomBZ}`,
    };
  }
  if (nPoints < 1) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `nPoints must be >= 1, got ${nPoints}`,
    };
  }
  if (nPoints > 1 && rMin >= rMax) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `rMin (${rMin}) must be < rMax (${rMax})`,
    };
  }

  // Check abort before starting
  if (isAborted()) {
    return {
      points: [],
      equilibrium: null,
      compute_time_ms: 0,
      total_iterations: 0,
    };
  }

  // Build WASM input (camelCase for serde deserialization)
  const wasmInput = {
    atomAZ,
    atomBZ,
    rMin,
    rMax,
    nPoints,
    basisName,
    options: {
      convergenceProfile: options.convergenceProfile,
      maxIterations: options.maxIterations,
      useDiis: options.useDiis,
      diisSize: options.diisSize,
      damp: options.damp,
      includeMatrices: false,
    },
    useSeeding: useSeeding ?? true,
  };

  // Track whether abort was requested during the scan
  let abortedDuringProgress = false;
  const collectedPoints: PesWasmPoint[] = [];

  // Set up progress callback that wraps WASM progress into worker progress
  const wasmProgressCallback = (wasmProgress: PesScanWasmProgress): void => {
    // Check abort -- if aborted, record partial results
    if (isAborted()) {
      abortedDuringProgress = true;
      // We cannot stop the WASM execution, but we stop emitting progress
      return;
    }

    // Collect points for partial results on abort
    collectedPoints.push({
      r: wasmProgress.r,
      energy: wasmProgress.energy,
      converged: wasmProgress.converged,
      iterations: 0, // Not available from progress callback
    });

    onProgress({
      module: 'pes',
      pointIndex: wasmProgress.pointIndex,
      totalPoints: wasmProgress.totalPoints,
      r: wasmProgress.r,
      energy: wasmProgress.energy,
      converged: wasmProgress.converged,
      current: wasmProgress.pointIndex + 1,
      total: wasmProgress.totalPoints,
      message: wasmProgress.converged
        ? `Point ${wasmProgress.pointIndex + 1}/${wasmProgress.totalPoints}: r=${wasmProgress.r.toFixed(3)} E=${wasmProgress.energy.toFixed(8)} Ha`
        : `Point ${wasmProgress.pointIndex + 1}/${wasmProgress.totalPoints}: r=${wasmProgress.r.toFixed(3)} (not converged)`,
    });
  };

  // Call WASM pes_scan function
  let result: PesWasmResult;
  try {
    result = wasmPesScan(wasmInput, wasmProgressCallback);
  } catch (e) {
    // If aborted during execution, return partial results
    if (abortedDuringProgress && collectedPoints.length > 0) {
      return {
        points: collectedPoints,
        equilibrium: null,
        compute_time_ms: 0,
        total_iterations: 0,
      };
    }

    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'PES scan computation failed',
    };
  }

  // If aborted during computation, return partial results from the WASM result
  // (WASM still runs to completion, but we only report what we have)
  if (abortedDuringProgress) {
    return {
      points: result.points,
      equilibrium: null,
      compute_time_ms: result.compute_time_ms,
      total_iterations: result.total_iterations,
    };
  }

  // Return the full result (structure matches PesScanResult directly)
  return {
    points: result.points,
    equilibrium: result.equilibrium,
    compute_time_ms: result.compute_time_ms,
    total_iterations: result.total_iterations,
  };
}
