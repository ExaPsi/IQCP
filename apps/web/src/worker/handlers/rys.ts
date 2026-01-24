/**
 * Rys quadrature handlers.
 *
 * Handle rys_compute and rys_error_curve requests by calling the WASM module.
 *
 * @module worker/handlers/rys
 */

import { rys_compute, rys_error_curve } from '../../wasm/qc_wasm';
import type {
  RysComputeRequest,
  RysErrorCurveRequest,
  RysComputeResult,
  RysErrorCurveResult,
  ErrorResponse,
} from '../protocol';

/**
 * Create an error response.
 */
function errorResponse(
  requestId: string,
  code: ErrorResponse['code'],
  message: string
): ErrorResponse {
  return {
    type: 'error',
    requestId: requestId as ErrorResponse['requestId'],
    code,
    message,
  };
}

/**
 * Handle rys_compute request (roots and weights).
 *
 * Calls the WASM rys_compute function to compute Rys quadrature roots and weights
 * for a given order n and parameter T.
 *
 * @param request - The rys_compute request containing n and T parameters
 * @returns The computed roots and weights, or an error response
 */
export function handleRysCompute(
  request: RysComputeRequest
): RysComputeResult | ErrorResponse {
  const { requestId, n, T } = request;

  // Validate parameters
  if (n < 1 || !Number.isInteger(n)) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid order n=${n}: must be a positive integer`);
  }
  if (T < 0) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid argument T=${T}: must be >= 0`);
  }

  try {
    // Call WASM function - note: WASM uses n (usize) and t (f64)
    const result = rys_compute(n, T) as RysComputeResult;
    return result;
  } catch (error) {
    const message = error instanceof Error ? error.message : 'Unknown WASM error';
    return errorResponse(requestId, 'HANDLER_ERROR', `Rys compute failed: ${message}`);
  }
}

/**
 * Handle rys_error_curve request (order vs error).
 *
 * Calls the WASM rys_error_curve function to compute reconstruction error
 * as a function of quadrature order for a given T parameter.
 *
 * @param request - The rys_error_curve request containing T and max_order parameters
 * @returns Arrays of orders and corresponding errors, or an error response
 */
export function handleRysErrorCurve(
  request: RysErrorCurveRequest
): RysErrorCurveResult | ErrorResponse {
  const { requestId, T, max_order } = request;

  // Validate parameters
  if (max_order < 1 || !Number.isInteger(max_order)) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid max_order=${max_order}: must be a positive integer`);
  }
  if (T < 0) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid argument T=${T}: must be >= 0`);
  }

  try {
    // Call WASM function
    const result = rys_error_curve(max_order, T) as RysErrorCurveResult;
    return result;
  } catch (error) {
    const message = error instanceof Error ? error.message : 'Unknown WASM error';
    return errorResponse(requestId, 'HANDLER_ERROR', `Rys error curve failed: ${message}`);
  }
}
