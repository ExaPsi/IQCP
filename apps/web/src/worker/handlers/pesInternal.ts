/**
 * Internal coordinate PES scan handler.
 *
 * Handles pes_scan_internal requests by calling the WASM pes_scan_internal
 * function. Supports progress streaming and cancellation. Supports both
 * rigid and relaxed scan modes over bonds, angles, and dihedrals.
 *
 * @module worker/handlers/pesInternal
 * @see US-081 PES Scan WASM Export + Worker Handler
 */

import type {
  PesScanInternalRequest,
  PesScanInternalResult,
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
 * WASM PES scan internal progress structure (matches PesScanInternalWasmProgress in qc-wasm).
 *
 * Note: Field names use camelCase as they come from serde serialization
 * with rename_all = "camelCase".
 */
interface PesScanInternalWasmProgress {
  /** Index of the completed point (0-indexed) */
  pointIndex: number;
  /** Total number of scan points */
  totalPoints: number;
  /** Value of the scanned coordinate */
  coordinateValue: number;
  /** SCF energy in Hartree */
  energy: number;
  /** Whether SCF converged at this geometry */
  converged: boolean;
  /** Optimization steps (null for rigid scans) */
  optSteps: number | null;
}

/**
 * Type for the WASM pes_scan_internal function.
 *
 * This is the function signature exported from qc-wasm.
 */
export type PesScanInternalFn = (
  input: unknown,
  progressCallback: (progress: PesScanInternalWasmProgress) => void
) => PesScanInternalResult;

// Store reference to the WASM pes_scan_internal function (set during worker initialization)
let wasmPesScanInternal: PesScanInternalFn | null = null;

/**
 * Set the WASM pes_scan_internal function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function is stored and used by handlePesScanInternal.
 *
 * @param fn - The WASM pes_scan_internal function
 */
export function setPesScanInternalFn(fn: PesScanInternalFn): void {
  wasmPesScanInternal = fn;
}

/**
 * Expected number of atom indices for each coordinate type.
 */
const COORDINATE_ATOM_COUNTS: Record<string, number> = {
  bond: 2,
  angle: 3,
  dihedral: 4,
};

/**
 * Handle pes_scan_internal request.
 *
 * Runs an internal coordinate PES scan by calling the WASM pes_scan_internal
 * function. Progress is emitted for each scan point, and the computation can
 * be cancelled via the abort checker.
 *
 * @param request - The PES scan internal request
 * @param onProgress - Callback for progress updates
 * @param isAborted - Function to check if computation should abort
 * @returns PES scan internal result or error response
 *
 * @example
 * ```typescript
 * const result = handlePesScanInternal(
 *   { type: 'pes_scan_internal', requestId: id, atoms: [...], ... },
 *   (progress) => postMessage({ type: 'progress', requestId: id, progress }),
 *   () => abortFlags.get(id) === true
 * );
 * ```
 */
export function handlePesScanInternal(
  request: PesScanInternalRequest,
  onProgress: ProgressCallback,
  isAborted: AbortChecker
): PesScanInternalResult | ErrorResponse {
  const {
    requestId,
    atoms,
    basisName,
    method,
    coordinateType,
    atomIndices,
    scanMode,
    valueMin,
    valueMax,
    nPoints,
    useSeeding,
    useSpherical,
    convergenceProfile,
    optMaxSteps,
    optGradThreshold,
  } = request;

  // Check if WASM function is available
  if (!wasmPesScanInternal) {
    return {
      type: 'error',
      requestId,
      code: 'NOT_IMPLEMENTED',
      message: 'pes_scan_internal WASM function not available (rebuild WASM module)',
    };
  }

  // Validate input parameters
  if (!atoms || atoms.length === 0) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: 'atoms array is empty',
    };
  }

  const validMethods = ['rhf', 'hf', 'lda', 'b3lyp', 'b3lyp-d3bj'];
  if (!validMethods.includes(method.toLowerCase())) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `Invalid method '${method}', expected one of: ${validMethods.join(', ')}`,
    };
  }

  const validCoordinateTypes = ['bond', 'angle', 'dihedral'];
  if (!validCoordinateTypes.includes(coordinateType)) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `Invalid coordinateType '${coordinateType}', expected one of: ${validCoordinateTypes.join(', ')}`,
    };
  }

  const expectedAtomCount = COORDINATE_ATOM_COUNTS[coordinateType];
  if (atomIndices.length !== expectedAtomCount) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `coordinateType '${coordinateType}' requires ${expectedAtomCount} atom indices, got ${atomIndices.length}`,
    };
  }

  const validScanModes = ['rigid', 'relaxed'];
  if (!validScanModes.includes(scanMode)) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `Invalid scanMode '${scanMode}', expected one of: ${validScanModes.join(', ')}`,
    };
  }

  if (nPoints < 2) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `nPoints must be >= 2, got ${nPoints}`,
    };
  }

  if (valueMin >= valueMax) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `valueMin (${valueMin}) must be < valueMax (${valueMax})`,
    };
  }

  // Check abort before starting
  if (isAborted()) {
    return {
      coordinate_type: coordinateType,
      atom_indices: atomIndices,
      points: [],
      equilibrium: null,
      total_iterations: 0,
      scan_mode: scanMode,
      total_opt_steps: 0,
    };
  }

  // Build WASM input (camelCase for serde deserialization)
  const wasmInput: Record<string, unknown> = {
    atoms,
    basisName,
    method: method.toLowerCase(),
    coordinateType,
    atomIndices,
    scanMode,
    valueMin,
    valueMax,
    nPoints,
    useSeeding: useSeeding ?? true,
    useSpherical: useSpherical ?? true,
    convergenceProfile: convergenceProfile ?? 'tight',
  };

  // Only include relaxed scan optimization parameters when provided
  if (optMaxSteps !== undefined) {
    wasmInput.optMaxSteps = optMaxSteps;
  }
  if (optGradThreshold !== undefined) {
    wasmInput.optGradThreshold = optGradThreshold;
  }

  // Track whether abort was requested during the scan
  let abortedDuringProgress = false;

  // Set up progress callback that wraps WASM progress into worker progress
  const wasmProgressCallback = (wasmProgress: PesScanInternalWasmProgress): void => {
    // Check abort -- if aborted, record and stop emitting
    if (isAborted()) {
      abortedDuringProgress = true;
      return;
    }

    onProgress({
      module: 'pes_internal',
      pointIndex: wasmProgress.pointIndex,
      totalPoints: wasmProgress.totalPoints,
      coordinateValue: wasmProgress.coordinateValue,
      energy: wasmProgress.energy,
      converged: wasmProgress.converged,
      optSteps: wasmProgress.optSteps,
      current: wasmProgress.pointIndex + 1,
      total: wasmProgress.totalPoints,
      message: wasmProgress.converged
        ? `Point ${wasmProgress.pointIndex + 1}/${wasmProgress.totalPoints}: E = ${wasmProgress.energy.toFixed(8)} Ha`
        : `Point ${wasmProgress.pointIndex + 1}/${wasmProgress.totalPoints}: (not converged)`,
    });
  };

  // Call WASM pes_scan_internal function with timing
  let result: PesScanInternalResult;
  const startTime = performance.now();
  try {
    result = wasmPesScanInternal(wasmInput, wasmProgressCallback);
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'Internal coordinate PES scan failed',
    };
  }
  const elapsedMs = performance.now() - startTime;

  // Inject compute time into result
  (result as PesScanInternalResult & { compute_time_ms: number }).compute_time_ms = elapsedMs;

  // If aborted during computation, return result but nullify equilibrium
  if (abortedDuringProgress) {
    return {
      ...result,
      equilibrium: null,
    };
  }

  return result;
}
