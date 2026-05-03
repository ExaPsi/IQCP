/**
 * Fock decomposition handler for computing J, K, G, F matrix decomposition via WASM.
 *
 * Wraps the `fock_decomposition` WASM function to decompose the Fock matrix
 * into H^core, J (Coulomb), K (Exchange), G = J - 0.5*K, and F = H^core + G.
 * Used by Module E's FockBuildTracer to show students the physical components
 * of the Fock matrix.
 *
 * @module worker/handlers/fockDecomposition
 * @see US-058 Fock Build Tracing
 */

import type {
  FockDecompositionRequest,
  FockDecompositionResult,
  ErrorResponse,
} from '../protocol';

// ============================================================================
// WASM Function Types
// ============================================================================

/**
 * Type signature for the WASM fock_decomposition function.
 *
 * Takes a single input object (deserialized by serde-wasm-bindgen)
 * and returns the FockDecompositionResult.
 *
 * @param input - Object with atoms, basisName, units, densityMatrix
 * @returns FockDecompositionResult with hCore, jMatrix, kMatrix, gMatrix, fMatrix, density, nbf, labels
 */
export type FockDecompositionFn = (input: unknown) => FockDecompositionResult;

// ============================================================================
// Module State
// ============================================================================

/** Injected WASM function reference */
let wasmFockDecomposition: FockDecompositionFn | null = null;

// ============================================================================
// Setup
// ============================================================================

/**
 * Inject the WASM fock_decomposition function.
 *
 * Called during worker initialization to wire up the WASM dependency.
 *
 * @param fn - The fock_decomposition function from the WASM module
 */
export function setFockDecompositionFn(fn: FockDecompositionFn): void {
  wasmFockDecomposition = fn;
}

// ============================================================================
// Handler
// ============================================================================

/**
 * Handle a fock_decomposition request.
 *
 * Builds the WASM input from the request fields and calls
 * fock_decomposition. Returns the decomposed Fock matrix with
 * separate J (Coulomb) and K (Exchange) contributions.
 *
 * @param request - The Fock decomposition request
 * @returns FockDecompositionResult or error response
 */
export function handleFockDecomposition(
  request: FockDecompositionRequest
): FockDecompositionResult | ErrorResponse {
  if (!wasmFockDecomposition) {
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'WORKER_NOT_READY',
      message: 'WASM fock_decomposition function not available. Rebuild WASM?',
    };
  }

  try {
    // Build the input matching the Rust FockDecompositionInput struct (camelCase).
    const input = {
      atoms: request.geometry.atoms,
      basisName: request.basisSet,
      units: request.geometry.units,
      densityMatrix: request.densityMatrix,
    };

    return wasmFockDecomposition(input);
  } catch (error) {
    const message =
      error instanceof Error
        ? error.message
        : 'Unknown error in fock_decomposition';
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'HANDLER_ERROR',
      message,
    };
  }
}
