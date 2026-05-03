/**
 * Integral matrices handler for computing S, T, V, H^core via WASM.
 *
 * Wraps the `compute_integral_matrices` WASM function to provide
 * one-electron integral matrices for Module E's heatmap visualization.
 * Unlike the integral_compute handler, this returns T and V separately
 * and does not compute ERIs, making it faster for display purposes.
 *
 * @module worker/handlers/integralMatrices
 * @see US-055 Integral Matrix Heatmap
 */

import type {
  IntegralMatricesRequest,
  IntegralMatricesResult,
  ErrorResponse,
} from '../protocol';

// ============================================================================
// WASM Function Types
// ============================================================================

/**
 * Type signature for the WASM compute_integral_matrices function.
 *
 * Takes a single input object (deserialized by serde-wasm-bindgen)
 * and returns the IntegralMatricesResult.
 *
 * @param input - Object with atoms, basisName, units, useSpherical
 * @returns IntegralMatricesResult with matrices, labels, nbf
 */
export type ComputeIntegralMatricesFn = (input: unknown) => IntegralMatricesResult;

// ============================================================================
// Module State
// ============================================================================

/** Injected WASM function reference */
let wasmComputeIntegralMatrices: ComputeIntegralMatricesFn | null = null;

// ============================================================================
// Setup
// ============================================================================

/**
 * Inject the WASM compute_integral_matrices function.
 *
 * Called during worker initialization to wire up the WASM dependency.
 *
 * @param fn - The compute_integral_matrices function from the WASM module
 */
export function setComputeIntegralMatricesFn(fn: ComputeIntegralMatricesFn): void {
  wasmComputeIntegralMatrices = fn;
}

// ============================================================================
// Handler
// ============================================================================

/**
 * Handle an integral_matrices request.
 *
 * Builds the WASM input from the request fields and calls
 * compute_integral_matrices. Returns the result containing S, T, V,
 * H^core matrices, basis function labels, and nuclear repulsion energy.
 *
 * @param request - The integral matrices request
 * @returns IntegralMatricesResult or error response
 */
export function handleIntegralMatrices(
  request: IntegralMatricesRequest
): IntegralMatricesResult | ErrorResponse {
  if (!wasmComputeIntegralMatrices) {
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'WORKER_NOT_READY',
      message: 'WASM compute_integral_matrices function not available. Rebuild WASM?',
    };
  }

  try {
    // Build the input matching the Rust IntegralMatricesInput struct (camelCase).
    // The request uses GeometryInput (with atoms + units), but WASM expects
    // atoms and units as top-level fields alongside basisName.
    const input: Record<string, unknown> = {
      atoms: request.geometry.atoms,
      basisName: request.basisName,
      units: request.geometry.units,
    };

    // Only include useSpherical when true (serde default is false)
    if (request.useSpherical) {
      input.useSpherical = true;
    }

    return wasmComputeIntegralMatrices(input);
  } catch (error) {
    const message =
      error instanceof Error
        ? error.message
        : 'Unknown error in compute_integral_matrices';
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'HANDLER_ERROR',
      message,
    };
  }
}
