/**
 * Density grid evaluation handler.
 *
 * Handles density grid evaluation requests by calling the WASM evaluate_density_grid
 * function. The result contains a flat array of electron density values on a 3D grid.
 *
 * Includes a simple single-entry cache keyed by density matrix fingerprint + grid
 * parameters, since the same density is typically evaluated multiple times with
 * different isovalues (only the marching cubes step changes, not the grid).
 *
 * @module worker/handlers/densityGrid
 * @see US-061 Density Isosurface Visualization
 */

import type { DensityGridRequest, DensityGridResult, ErrorResponse } from '../protocol';

/**
 * Type for the WASM evaluate_density_grid function.
 *
 * This is the function signature exported from qc-wasm.
 */
export type DensityGridFn = (input: unknown) => DensityGridResult;

// Store reference to the WASM evaluate_density_grid function (set during worker initialization)
let wasmDensityGrid: DensityGridFn | null = null;

/**
 * Set the WASM evaluate_density_grid function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function is stored and used by handleDensityGrid.
 *
 * @param fn - The WASM evaluate_density_grid function
 */
export function setDensityGridFn(fn: DensityGridFn): void {
  wasmDensityGrid = fn;
}

// ============================================================================
// Density Grid Caching (US-061-T7)
// ============================================================================

/**
 * Simple single-entry cache for density grid results.
 *
 * Caches the last computed result keyed by a fingerprint of the density matrix
 * and grid parameters. This avoids recomputation when only the isovalue changes
 * (which only affects the marching cubes step, not the grid evaluation).
 */
let cachedResult: DensityGridResult | null = null;
let cachedKey: string | null = null;

/**
 * Compute a cache key from the density grid request.
 *
 * Uses a lightweight fingerprint: density matrix length, first element,
 * last element, and sum (not a full comparison, but sufficient for the
 * common case where the density matrix changes between SCF runs).
 *
 * @param request - The density grid request
 * @returns A string cache key
 */
function densityCacheKey(request: DensityGridRequest): string {
  const dLen = request.densityMatrix.length;
  const dFirst = dLen > 0 ? request.densityMatrix[0] : 0;
  const dLast = dLen > 0 ? request.densityMatrix[dLen - 1] : 0;
  // Use a running sum for fingerprinting; cheaper than hashing all elements
  let dSum = 0;
  for (let i = 0; i < dLen; i++) {
    dSum += request.densityMatrix[i];
  }
  return `${dLen}:${dFirst}:${dLast}:${dSum}:${request.gridDims.join(',')}:${request.gridSpacing}:${request.gridOrigin.join(',')}:${request.basisName}:${request.useSpherical ?? false}`;
}

/**
 * Invalidate the density grid cache.
 *
 * Call this when external state changes might invalidate cached results
 * (e.g., after a new SCF run).
 */
export function invalidateDensityGridCache(): void {
  cachedResult = null;
  cachedKey = null;
}

// ============================================================================
// Handler
// ============================================================================

/**
 * Handle density_grid request.
 *
 * Evaluates the electron density rho(r) on a 3D grid by calling the WASM
 * evaluate_density_grid function. Returns the grid values and metadata.
 *
 * Uses a single-entry cache to avoid recomputation when the density matrix
 * and grid parameters haven't changed.
 *
 * @param request - The density grid request
 * @returns Density grid result or error response
 *
 * @example
 * ```typescript
 * const result = handleDensityGrid({
 *   type: 'density_grid',
 *   requestId: id,
 *   densityMatrix: [0.6, 0.4, 0.4, 0.6],
 *   atoms: [[1, 0, 0, 0], [1, 0, 0, 1.4]],
 *   basisName: 'sto-3g',
 *   gridOrigin: [-5, -5, -5],
 *   gridSpacing: 0.5,
 *   gridDims: [21, 21, 23],
 *   nElectrons: 2,
 * });
 * ```
 */
export function handleDensityGrid(
  request: DensityGridRequest
): DensityGridResult | ErrorResponse {
  const {
    requestId,
    densityMatrix,
    atoms,
    basisName,
    gridOrigin,
    gridSpacing,
    gridDims,
    nElectrons,
    useSpherical,
  } = request;

  // Check if WASM function is available
  if (!wasmDensityGrid) {
    return {
      type: 'error',
      requestId,
      code: 'WORKER_NOT_READY',
      message: 'Density grid evaluation WASM function not initialized',
    };
  }

  // Validate input parameters
  if (!densityMatrix || densityMatrix.length === 0) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: 'densityMatrix must be a non-empty array',
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

  if (nElectrons < 1) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `nElectrons must be >= 1, got ${nElectrons}`,
    };
  }

  // Check cache
  const key = densityCacheKey(request);
  if (cachedKey === key && cachedResult) {
    // Return cached result (shallow copy is fine since values are immutable once computed)
    return cachedResult;
  }

  // Build WASM input (camelCase for serde deserialization)
  const wasmInput: Record<string, unknown> = {
    densityMatrix,
    atoms,
    basisName,
    gridOrigin,
    gridSpacing,
    gridDims,
    nElectrons,
  };
  // Only include useSpherical when true (serde default is false)
  if (useSpherical) {
    wasmInput.useSpherical = true;
  }

  // Call WASM evaluate_density_grid function
  try {
    const result = wasmDensityGrid(wasmInput);

    // Cache the result
    cachedKey = key;
    cachedResult = result;

    return result;
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'Density grid evaluation failed',
    };
  }
}
