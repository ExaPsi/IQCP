/**
 * Radial profile handler for evaluating basis function radial profiles from WASM.
 *
 * Wraps the `evaluate_radial_profile` WASM function to compute the radial
 * part of contracted Gaussian basis functions for Module D's radial profile plot.
 *
 * @module worker/handlers/radialProfile
 * @see US-050 Radial Profile Visualization
 */

import type {
  RadialProfileRequest,
  RadialProfileResult,
  ErrorResponse,
} from '../protocol';

// ============================================================================
// WASM Function Types
// ============================================================================

/**
 * Type signature for the WASM evaluate_radial_profile function.
 *
 * @param input - Serialized RadialProfileInput object
 * @returns RadialProfileResult with r values, profiles, and metadata
 */
export type EvaluateRadialProfileFn = (input: unknown) => RadialProfileResult;

// ============================================================================
// Module State
// ============================================================================

/** Injected WASM function reference */
let wasmEvaluateRadialProfile: EvaluateRadialProfileFn | null = null;

// ============================================================================
// Setup
// ============================================================================

/**
 * Inject the WASM evaluate_radial_profile function.
 *
 * Called during worker initialization to wire up the WASM dependency.
 *
 * @param fn - The evaluate_radial_profile function from the WASM module
 */
export function setEvaluateRadialProfileFn(fn: EvaluateRadialProfileFn): void {
  wasmEvaluateRadialProfile = fn;
}

// ============================================================================
// Handler
// ============================================================================

/**
 * Handle a radial_profile request.
 *
 * Evaluates the radial profile of a contracted basis shell for the
 * specified element, basis set, and shell index.
 *
 * @param request - The radial profile request
 * @returns RadialProfileResult or error response
 */
export function handleRadialProfile(
  request: RadialProfileRequest
): RadialProfileResult | ErrorResponse {
  if (!wasmEvaluateRadialProfile) {
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'WORKER_NOT_READY',
      message: 'WASM evaluate_radial_profile function not available. Rebuild WASM?',
    };
  }

  try {
    // Build the input matching the Rust RadialProfileInput struct
    const input: Record<string, unknown> = {
      atomicNumber: request.atomicNumber,
      basisName: request.basisName,
      shellIndex: request.shellIndex,
    };

    // Only include optional fields if provided
    if (request.nPoints !== undefined) {
      input.nPoints = request.nPoints;
    }
    if (request.rMax !== undefined) {
      input.rMax = request.rMax;
    }

    return wasmEvaluateRadialProfile(input);
  } catch (error) {
    const message = error instanceof Error ? error.message : 'Unknown error in evaluate_radial_profile';
    return {
      type: 'error',
      requestId: request.requestId,
      code: 'HANDLER_ERROR',
      message,
    };
  }
}
