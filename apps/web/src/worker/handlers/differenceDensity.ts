/**
 * Difference density grid evaluation handler.
 *
 * Handles difference_density requests by calling the WASM
 * compute_difference_density function. The handler receives the cached
 * total molecular density grid and atom positions, evaluates the
 * promolecule density internally (using embedded atomic density data),
 * and returns Delta-rho = rho_molecule - rho_promolecule.
 *
 * @module worker/handlers/differenceDensity
 * @see US-063 Difference Density
 */

import type { DifferenceDensityRequest, DifferenceDensityResult, ErrorResponse } from '../protocol';

/**
 * Type for the WASM compute_difference_density function.
 *
 * This is the function signature exported from qc-wasm.
 * It takes a JsValue input and returns a JsValue result.
 */
export type DifferenceDensityFn = (input: unknown) => DifferenceDensityResult;

// Store reference to the WASM function (set during worker initialization)
let wasmDifferenceDensity: DifferenceDensityFn | null = null;

/**
 * Set the WASM compute_difference_density function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function is stored and used by handleDifferenceDensity.
 *
 * @param fn - The WASM compute_difference_density function
 */
export function setDifferenceDensityFn(fn: DifferenceDensityFn): void {
  wasmDifferenceDensity = fn;
}

// ============================================================================
// Handler
// ============================================================================

/**
 * Handle difference_density request.
 *
 * Computes Delta-rho = rho_molecule - rho_promolecule on the provided grid.
 * The molecular density values are passed in from a cached density_grid result.
 * The promolecule is evaluated inside WASM using embedded atomic density data.
 *
 * @param request - The difference density request
 * @returns Difference density result or error response
 */
export function handleDifferenceDensity(
  request: DifferenceDensityRequest
): DifferenceDensityResult | ErrorResponse {
  const {
    requestId,
    totalDensity,
    atoms,
    gridOrigin,
    gridSpacing,
    gridDims,
  } = request;

  // Check if WASM function is available
  if (!wasmDifferenceDensity) {
    return {
      type: 'error',
      requestId,
      code: 'WORKER_NOT_READY',
      message: 'Difference density WASM function not initialized. Rebuild WASM?',
    };
  }

  // Validate input parameters
  if (!totalDensity || totalDensity.length === 0) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: 'totalDensity must be a non-empty array',
    };
  }

  if (!atoms || atoms.length === 0) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: 'atoms must be a non-empty array',
    };
  }

  if (gridSpacing <= 0) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `gridSpacing must be positive, got ${gridSpacing}`,
    };
  }

  if (gridDims.some((d) => d < 2)) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `gridDims must all be >= 2, got [${gridDims.join(', ')}]`,
    };
  }

  // Verify grid size matches total density length
  const expectedSize = gridDims[0] * gridDims[1] * gridDims[2];
  if (totalDensity.length !== expectedSize) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `totalDensity length (${totalDensity.length}) does not match gridDims product (${expectedSize})`,
    };
  }

  // Build WASM input (camelCase for serde deserialization)
  const wasmInput = {
    totalDensity,
    atoms,
    gridOrigin,
    gridSpacing,
    gridDims,
  };

  // Call WASM compute_difference_density function
  try {
    const result = wasmDifferenceDensity(wasmInput);
    return result;
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'Difference density evaluation failed',
    };
  }
}
