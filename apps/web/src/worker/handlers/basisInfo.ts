/**
 * Basis info handler for querying basis set shell data from WASM.
 *
 * Wraps the `get_basis_info` WASM function to provide basis set shell
 * information sourced from qc-core's builtin.rs, ensuring a single
 * source of truth for basis set data.
 *
 * @module worker/handlers/basisInfo
 * @see US-049 AC6
 */

import type {
  BasisInfoRequest,
  BasisInfoResult,
  ErrorResponse,
} from '../protocol';

// ============================================================================
// WASM Function Types
// ============================================================================

/**
 * Type signature for the WASM get_basis_info function.
 *
 * @param atomicNumber - Atomic number (1-18)
 * @param basisName - Basis set name (e.g., "sto-3g")
 * @returns Array of shell info objects
 */
export type GetBasisInfoFn = (atomicNumber: number, basisName: string) => BasisInfoResult;

// ============================================================================
// Module State
// ============================================================================

/** Injected WASM function reference */
let wasmGetBasisInfo: GetBasisInfoFn | null = null;

// ============================================================================
// Setup
// ============================================================================

/**
 * Inject the WASM get_basis_info function.
 *
 * Called during worker initialization to wire up the WASM dependency.
 *
 * @param fn - The get_basis_info function from the WASM module
 */
export function setGetBasisInfoFn(fn: GetBasisInfoFn): void {
  wasmGetBasisInfo = fn;
}

// ============================================================================
// Handler
// ============================================================================

/**
 * Handle a basis_info request.
 *
 * Queries the WASM module for shell structure data for the given
 * element and basis set combination.
 *
 * @param request - The basis info request
 * @returns Shell info array or error response
 */
export function handleBasisInfo(
  request: BasisInfoRequest
): BasisInfoResult | ErrorResponse {
  if (!wasmGetBasisInfo) {
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'WORKER_NOT_READY',
      message: 'WASM get_basis_info function not available. Rebuild WASM?',
    };
  }

  try {
    return wasmGetBasisInfo(request.atomicNumber, request.basisName);
  } catch (error) {
    const message = error instanceof Error ? error.message : 'Unknown error in get_basis_info';
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'HANDLER_ERROR',
      message,
    };
  }
}
