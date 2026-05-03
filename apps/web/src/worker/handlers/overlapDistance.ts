/**
 * Overlap vs. distance handler for computing overlap integrals at varying
 * interatomic distances from WASM.
 *
 * Wraps the `overlap_vs_distance` WASM function to compute the overlap
 * integral between two basis shells as a function of distance for Module D's
 * overlap vs. distance plot.
 *
 * @module worker/handlers/overlapDistance
 * @see US-054 Overlap vs. Distance Plot
 */

import type {
  OverlapDistanceRequest,
  OverlapDistanceResult,
  ErrorResponse,
} from '../protocol';

// ============================================================================
// WASM Function Types
// ============================================================================

/**
 * Type signature for the WASM overlap_vs_distance function.
 *
 * @param input - Serialized OverlapDistanceInput object
 * @returns OverlapDistanceResult with distance values, overlap values, and labels
 */
export type OverlapVsDistanceFn = (input: unknown) => OverlapDistanceResult;

// ============================================================================
// Module State
// ============================================================================

/** Injected WASM function reference */
let wasmOverlapVsDistance: OverlapVsDistanceFn | null = null;

// ============================================================================
// Setup
// ============================================================================

/**
 * Inject the WASM overlap_vs_distance function.
 *
 * Called during worker initialization to wire up the WASM dependency.
 *
 * @param fn - The overlap_vs_distance function from the WASM module
 */
export function setOverlapVsDistanceFn(fn: OverlapVsDistanceFn): void {
  wasmOverlapVsDistance = fn;
}

// ============================================================================
// Handler
// ============================================================================

/**
 * Handle an overlap_distance request.
 *
 * Evaluates the overlap integral between two basis shells at multiple
 * interatomic distances for the specified element/basis/shell pairs.
 *
 * @param request - The overlap distance request
 * @returns OverlapDistanceResult or error response
 */
export function handleOverlapDistance(
  request: OverlapDistanceRequest
): OverlapDistanceResult | ErrorResponse {
  if (!wasmOverlapVsDistance) {
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'WORKER_NOT_READY',
      message: 'WASM overlap_vs_distance function not available. Rebuild WASM?',
    };
  }

  try {
    // Build the input matching the Rust OverlapDistanceInput struct
    const input = {
      elementA: request.elementA,
      basisA: request.basisA,
      shellIndexA: request.shellIndexA,
      elementB: request.elementB,
      basisB: request.basisB,
      shellIndexB: request.shellIndexB,
      rMin: request.rMin,
      rMax: request.rMax,
      nPoints: request.nPoints,
    };

    return wasmOverlapVsDistance(input);
  } catch (error) {
    const message = error instanceof Error ? error.message : 'Unknown error in overlap_vs_distance';
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'HANDLER_ERROR',
      message,
    };
  }
}
