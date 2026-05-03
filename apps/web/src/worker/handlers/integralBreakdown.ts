/**
 * Integral breakdown handler for computing primitive-pair decompositions via WASM.
 *
 * Wraps the `integral_with_breakdown` WASM function to decompose a contracted
 * integral into its constituent primitive-pair contributions. Used by Module E's
 * PrimitiveBreakdownTable to show students which Gaussian pairs dominate.
 *
 * @module worker/handlers/integralBreakdown
 * @see US-056 Primitive Breakdown Panel
 */

import type {
  IntegralBreakdownRequest,
  IntegralBreakdownResult,
  ErrorResponse,
} from '../protocol';

// ============================================================================
// WASM Function Types
// ============================================================================

/**
 * Type signature for the WASM integral_with_breakdown function.
 *
 * Takes a single input object (deserialized by serde-wasm-bindgen)
 * and returns the IntegralBreakdownResult.
 *
 * @param input - Object with atoms, basisName, units, integralType, indices
 * @returns IntegralBreakdownResult with contracted value and primitive contributions
 */
export type IntegralBreakdownFn = (input: unknown) => IntegralBreakdownResult;

// ============================================================================
// Module State
// ============================================================================

/** Injected WASM function reference */
let wasmIntegralBreakdown: IntegralBreakdownFn | null = null;

// ============================================================================
// Setup
// ============================================================================

/**
 * Inject the WASM integral_with_breakdown function.
 *
 * Called during worker initialization to wire up the WASM dependency.
 *
 * @param fn - The integral_with_breakdown function from the WASM module
 */
export function setIntegralBreakdownFn(fn: IntegralBreakdownFn): void {
  wasmIntegralBreakdown = fn;
}

// ============================================================================
// Handler
// ============================================================================

/**
 * Handle an integral_breakdown request.
 *
 * Builds the WASM input from the request fields and calls
 * integral_with_breakdown. Returns the contracted integral value
 * and all primitive-pair contributions sorted by magnitude.
 *
 * @param request - The integral breakdown request
 * @returns IntegralBreakdownResult or error response
 */
export function handleIntegralBreakdown(
  request: IntegralBreakdownRequest
): IntegralBreakdownResult | ErrorResponse {
  if (!wasmIntegralBreakdown) {
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'WORKER_NOT_READY',
      message: 'WASM integral_with_breakdown function not available. Rebuild WASM?',
    };
  }

  try {
    // Build the input matching the Rust IntegralBreakdownInput struct (camelCase).
    const input = {
      atoms: request.geometry.atoms,
      basisName: request.basisName,
      units: request.geometry.units,
      integralType: request.integralType,
      indices: request.indices,
    };

    return wasmIntegralBreakdown(input);
  } catch (error) {
    const message =
      error instanceof Error
        ? error.message
        : 'Unknown error in integral_with_breakdown';
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'HANDLER_ERROR',
      message,
    };
  }
}
