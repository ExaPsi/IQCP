/**
 * MO grid evaluation handler.
 *
 * Handles MO grid evaluation requests by calling the WASM evaluate_mo_grid
 * function. The result contains a flat array of orbital values on a 3D grid.
 *
 * @module worker/handlers/moGrid
 * @see US-042 MO Grid Evaluation
 */

import type { MoGridRequest, MoGridResult, ErrorResponse } from '../protocol';

/**
 * Type for the WASM evaluate_mo_grid function.
 *
 * This is the function signature exported from qc-wasm.
 */
export type MoGridFn = (input: unknown) => MoGridResult;

// Store reference to the WASM evaluate_mo_grid function (set during worker initialization)
let wasmMoGrid: MoGridFn | null = null;

/**
 * Set the WASM evaluate_mo_grid function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function is stored and used by handleMoGrid.
 *
 * @param fn - The WASM evaluate_mo_grid function
 */
export function setMoGridFn(fn: MoGridFn): void {
  wasmMoGrid = fn;
}

/**
 * Handle mo_grid request.
 *
 * Evaluates a molecular orbital on a 3D grid by calling the WASM
 * evaluate_mo_grid function. Returns the grid values and metadata.
 *
 * @param request - The MO grid request
 * @returns MO grid result or error response
 *
 * @example
 * ```typescript
 * const result = handleMoGrid({
 *   type: 'mo_grid',
 *   requestId: id,
 *   moCoefficients: [0.549, 0.549],
 *   atoms: [[1, 0, 0, 0], [1, 0, 0, 1.4]],
 *   basisName: 'sto-3g',
 *   gridOrigin: [-5, -5, -5],
 *   gridSpacing: 0.5,
 *   gridDims: [21, 21, 23],
 * });
 * ```
 */
export function handleMoGrid(
  request: MoGridRequest
): MoGridResult | ErrorResponse {
  const { requestId, moCoefficients, atoms, basisName, gridOrigin, gridSpacing, gridDims, useSpherical } =
    request;

  // Check if WASM function is available
  if (!wasmMoGrid) {
    return {
      type: 'error',
      requestId,
      code: 'WORKER_NOT_READY',
      message: 'MO grid evaluation WASM function not initialized',
    };
  }

  // Validate input parameters
  if (!moCoefficients || moCoefficients.length === 0) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: 'moCoefficients must be a non-empty array',
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

  // Build WASM input (camelCase for serde deserialization)
  const wasmInput: Record<string, unknown> = {
    moCoefficients,
    atoms,
    basisName,
    gridOrigin,
    gridSpacing,
    gridDims,
  };
  // Only include useSpherical when true (serde default is false)
  if (useSpherical) {
    wasmInput.useSpherical = true;
  }

  // Call WASM evaluate_mo_grid function
  try {
    const result = wasmMoGrid(wasmInput);
    return result;
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'MO grid evaluation failed',
    };
  }
}
