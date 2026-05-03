/**
 * Geometry optimization handler.
 *
 * Handles optimize_geometry requests by calling the WASM optimize_geometry
 * function. Supports progress streaming per optimization step.
 *
 * @module worker/handlers/optimizeGeometry
 * @see US-075 Optimization UI + Trajectory Animation
 */

import type {
  OptimizeGeometryRequest,
  OptimizationResult,
  WorkerProgress,
  ErrorResponse,
} from '../protocol';

/**
 * Progress callback type for optimization updates.
 */
export type ProgressCallback = (progress: WorkerProgress) => void;

/**
 * Abort flag checker type.
 */
export type AbortChecker = () => boolean;

/**
 * WASM optimization progress structure (matches OptimizeWasmProgress in qc-wasm).
 *
 * Note: Field names use camelCase as they come from serde serialization.
 */
interface OptimizeWasmProgress {
  /** Current step number (0 = initial evaluation) */
  step: number;
  /** Total energy at this step (Ha) */
  energy: number;
  /** Maximum gradient component (Ha/bohr) */
  maxGradient: number;
  /** RMS gradient (Ha/bohr) */
  rmsGradient: number;
}

/**
 * Type for the WASM optimize_geometry function.
 *
 * This is the function signature exported from qc-wasm.
 */
export type OptimizeGeometryFn = (
  input: unknown,
  progressCallback: (progress: OptimizeWasmProgress) => void
) => OptimizationResult;

// Store reference to the WASM optimize_geometry function (set during worker initialization)
let wasmOptimizeGeometry: OptimizeGeometryFn | null = null;

/**
 * Set the WASM optimize_geometry function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function is stored and used by handleOptimizeGeometry.
 *
 * @param fn - The WASM optimize_geometry function
 */
export function setOptimizeGeometryFn(fn: OptimizeGeometryFn): void {
  wasmOptimizeGeometry = fn;
}

/**
 * Handle optimize_geometry request.
 *
 * Runs L-BFGS geometry optimization by calling the WASM optimize_geometry
 * function. Progress is emitted for each optimization step.
 *
 * @param request - The optimization request
 * @param onProgress - Callback for progress updates
 * @param isAborted - Function to check if computation should abort
 * @returns Optimization result or error response
 */
export function handleOptimizeGeometry(
  request: OptimizeGeometryRequest,
  onProgress: ProgressCallback,
  isAborted: AbortChecker
): OptimizationResult | ErrorResponse {
  const { requestId, atoms, basisName, method, maxSteps, gradThreshold, energyThreshold } = request;

  // Check if WASM function is available
  if (!wasmOptimizeGeometry) {
    return {
      type: 'error',
      requestId,
      code: 'NOT_IMPLEMENTED',
      message: 'optimize_geometry WASM function not available (rebuild WASM module)',
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
      message: `Invalid method '${method}', expected one of: rhf, lda, b3lyp, b3lyp-d3bj`,
    };
  }

  // Check abort before starting
  if (isAborted()) {
    return {
      converged: false,
      steps: [],
      finalEnergy: 0,
      finalGeometry: [],
      totalSteps: 0,
      computeTimeMs: 0,
    };
  }

  // Build WASM input (camelCase for serde deserialization)
  const wasmInput = {
    atoms,
    basisName,
    method: method.toLowerCase(),
    maxSteps: maxSteps ?? 50,
    gradThreshold: gradThreshold ?? 4.5e-4,
    energyThreshold: energyThreshold ?? 1.0e-6,
    memorySize: 7,
  };

  // Track abort state during progress
  let abortedDuringProgress = false;

  // Set up progress callback that wraps WASM progress into worker progress
  const wasmProgressCallback = (wasmProgress: OptimizeWasmProgress): void => {
    // Check abort
    if (isAborted()) {
      abortedDuringProgress = true;
      return;
    }

    onProgress({
      module: 'optimization',
      step: wasmProgress.step,
      energy: wasmProgress.energy,
      maxGradient: wasmProgress.maxGradient,
      rmsGradient: wasmProgress.rmsGradient,
      current: wasmProgress.step,
      total: 0, // Unknown total for optimization
      message: `Step ${wasmProgress.step}: E=${wasmProgress.energy.toFixed(8)} Ha, max|g|=${wasmProgress.maxGradient.toExponential(2)}`,
    });
  };

  // Call WASM optimize_geometry function
  try {
    const result = wasmOptimizeGeometry(wasmInput, wasmProgressCallback);

    // If aborted during execution, mark as not converged but return trajectory
    if (abortedDuringProgress) {
      return {
        ...result,
        converged: false,
      };
    }

    return result;
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'Geometry optimization failed',
    };
  }
}
