/**
 * Marching cubes isosurface extraction handler.
 *
 * Handles marching cubes and dual marching cubes requests by calling
 * the WASM marching_cubes / dual_marching_cubes functions.
 *
 * @module worker/handlers/marchingCubes
 * @see US-043 Marching Cubes Isosurface
 */

import type {
  MarchingCubesRequest,
  DualMarchingCubesRequest,
  MarchingCubesResult,
  DualMarchingCubesResult,
  ErrorResponse,
} from '../protocol';

/**
 * Type for the WASM marching_cubes function.
 */
export type MarchingCubesFn = (input: unknown) => MarchingCubesResult;

/**
 * Type for the WASM dual_marching_cubes function.
 */
export type DualMarchingCubesFn = (input: unknown) => DualMarchingCubesResult;

// Store references to the WASM functions (set during worker initialization)
let wasmMarchingCubes: MarchingCubesFn | null = null;
let wasmDualMarchingCubes: DualMarchingCubesFn | null = null;

/**
 * Set the WASM marching_cubes function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function is stored and used by handleMarchingCubes.
 *
 * @param fn - The WASM marching_cubes function
 */
export function setMarchingCubesFn(fn: MarchingCubesFn): void {
  wasmMarchingCubes = fn;
}

/**
 * Set the WASM dual_marching_cubes function reference.
 *
 * @param fn - The WASM dual_marching_cubes function
 */
export function setDualMarchingCubesFn(fn: DualMarchingCubesFn): void {
  wasmDualMarchingCubes = fn;
}

/**
 * Handle marching_cubes request.
 *
 * Extracts an isosurface from a 3D scalar field using the standard
 * marching cubes algorithm with smooth vertex normals.
 *
 * @param request - The marching cubes request
 * @returns Marching cubes result or error response
 */
export function handleMarchingCubes(
  request: MarchingCubesRequest
): MarchingCubesResult | ErrorResponse {
  const { requestId, gridData, gridDims, gridOrigin, gridSpacing, isovalue } = request;

  // Check if WASM function is available
  if (!wasmMarchingCubes) {
    return {
      type: 'error',
      requestId,
      code: 'WORKER_NOT_READY',
      message: 'Marching cubes WASM function not initialized',
    };
  }

  // Validate input parameters
  if (!gridData || gridData.length === 0) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: 'gridData must be a non-empty array',
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

  const expectedLen = gridDims[0] * gridDims[1] * gridDims[2];
  if (gridData.length !== expectedLen) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `gridData length ${gridData.length} does not match dims [${gridDims.join(', ')}] = ${expectedLen}`,
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

  // Build WASM input (camelCase for serde deserialization)
  const wasmInput = {
    gridData,
    gridDims,
    gridOrigin,
    gridSpacing,
    isovalue,
  };

  // Call WASM marching_cubes function
  try {
    const result = wasmMarchingCubes(wasmInput);
    return result;
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'Marching cubes extraction failed',
    };
  }
}

/**
 * Handle dual_marching_cubes request.
 *
 * Extracts dual isosurfaces (positive and negative orbital lobes)
 * from a 3D scalar field.
 *
 * @param request - The dual marching cubes request
 * @returns Dual marching cubes result or error response
 */
export function handleDualMarchingCubes(
  request: DualMarchingCubesRequest
): DualMarchingCubesResult | ErrorResponse {
  const { requestId, gridData, gridDims, gridOrigin, gridSpacing, isovalue } = request;

  // Check if WASM function is available
  if (!wasmDualMarchingCubes) {
    return {
      type: 'error',
      requestId,
      code: 'WORKER_NOT_READY',
      message: 'Dual marching cubes WASM function not initialized',
    };
  }

  // Validate input parameters
  if (!gridData || gridData.length === 0) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: 'gridData must be a non-empty array',
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

  if (gridSpacing <= 0) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `gridSpacing must be positive, got ${gridSpacing}`,
    };
  }

  // Build WASM input
  const wasmInput = {
    gridData,
    gridDims,
    gridOrigin,
    gridSpacing,
    isovalue,
  };

  // Call WASM dual_marching_cubes function
  try {
    const result = wasmDualMarchingCubes(wasmInput);
    return result;
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'Dual marching cubes extraction failed',
    };
  }
}
