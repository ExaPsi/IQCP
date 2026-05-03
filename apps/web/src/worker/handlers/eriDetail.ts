/**
 * ERI detail handler for computing primitive-quartet decompositions via WASM.
 *
 * Wraps the `eri_detail` WASM function to decompose a contracted ERI (ij|kl)
 * into its constituent primitive-quartet contributions. Used by Module E's
 * ERI Browser to show students the structure of two-electron integrals.
 *
 * @module worker/handlers/eriDetail
 * @see US-059 ERI Browser
 */

import type {
  EriDetailRequest,
  EriDetailResult,
  ErrorResponse,
} from '../protocol';

// ============================================================================
// WASM Function Types
// ============================================================================

/**
 * Type signature for the WASM eri_detail function.
 *
 * Takes a single input object (deserialized by serde-wasm-bindgen)
 * and returns the EriDetailResult.
 *
 * @param input - Object with atoms, basisName, units, indices
 * @returns EriDetailResult with contracted value and primitive contributions
 */
export type EriDetailFn = (input: unknown) => EriDetailResult;

// ============================================================================
// Module State
// ============================================================================

/** Injected WASM function reference */
let wasmEriDetail: EriDetailFn | null = null;

// ============================================================================
// Setup
// ============================================================================

/**
 * Inject the WASM eri_detail function.
 *
 * Called during worker initialization to wire up the WASM dependency.
 *
 * @param fn - The eri_detail function from the WASM module
 */
export function setEriDetailFn(fn: EriDetailFn): void {
  wasmEriDetail = fn;
}

// ============================================================================
// Handler
// ============================================================================

/**
 * Handle an eri_detail request.
 *
 * Builds the WASM input from the request fields and calls eri_detail.
 * Returns the contracted ERI value and all primitive-quartet contributions
 * sorted by magnitude, plus method metadata (Boys/Rys).
 *
 * @param request - The ERI detail request
 * @returns EriDetailResult or error response
 */
export function handleEriDetail(
  request: EriDetailRequest
): EriDetailResult | ErrorResponse {
  if (!wasmEriDetail) {
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'WORKER_NOT_READY',
      message: 'WASM eri_detail function not available. Rebuild WASM?',
    };
  }

  try {
    // Build the input matching the Rust EriDetailInput struct (camelCase).
    const input = {
      atoms: request.geometry.atoms,
      basisName: request.basisName,
      units: request.geometry.units,
      indices: request.indices,
    };

    return wasmEriDetail(input);
  } catch (error) {
    const message =
      error instanceof Error
        ? error.message
        : 'Unknown error in eri_detail';
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'HANDLER_ERROR',
      message,
    };
  }
}
