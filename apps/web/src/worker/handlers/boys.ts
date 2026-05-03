/**
 * Boys function handlers.
 *
 * Handle boys_eval and boys_sweep requests by calling the WASM module.
 *
 * @module worker/handlers/boys
 */

import { boys_eval, boys_eval_many, boys_eval_all } from '../../wasm/qc_wasm';
import type {
  BoysEvalRequest,
  BoysSweepRequest,
  BoysEvalAllRequest,
  BoysEvalResult,
  BoysSweepResult,
  BoysEvalAllResult,
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
 * Handle boys_eval request (single point evaluation).
 *
 * Calls the WASM boys_eval function and returns the result.
 */
export function handleBoysEval(
  request: BoysEvalRequest
): BoysEvalResult | ErrorResponse {
  const { requestId, m, T } = request;

  // Validate parameters
  if (m < 0 || !Number.isInteger(m)) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid order m=${m}: must be a non-negative integer`);
  }
  if (T < 0) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid argument T=${T}: must be >= 0`);
  }

  try {
    // Call WASM function
    const result = boys_eval(m, T) as BoysEvalResult;
    return result;
  } catch (error) {
    const message = error instanceof Error ? error.message : 'Unknown WASM error';
    return errorResponse(requestId, 'HANDLER_ERROR', `Boys evaluation failed: ${message}`);
  }
}

/**
 * Handle boys_sweep request (range evaluation).
 *
 * Generates T values across the specified range and evaluates
 * the Boys function at each point.
 */
export function handleBoysSweep(
  request: BoysSweepRequest
): BoysSweepResult | ErrorResponse {
  const { requestId, m, T_range, points } = request;

  // Validate parameters
  if (m < 0 || !Number.isInteger(m)) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid order m=${m}: must be a non-negative integer`);
  }
  if (T_range[0] < 0) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid T_min=${T_range[0]}: must be >= 0`);
  }
  if (T_range[1] < T_range[0]) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid range: T_max=${T_range[1]} must be >= T_min=${T_range[0]}`);
  }
  if (points < 1 || !Number.isInteger(points)) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid points=${points}: must be a positive integer`);
  }

  try {
    // Generate T values
    const [T_min, T_max] = T_range;
    const T_values: number[] = [];

    if (points === 1) {
      T_values.push(T_min);
    } else {
      const step = (T_max - T_min) / (points - 1);
      for (let i = 0; i < points; i++) {
        T_values.push(T_min + i * step);
      }
    }

    // Call WASM function with array
    const ts = new Float64Array(T_values);
    const results = boys_eval_many(m, ts) as BoysEvalResult[];

    return {
      results,
      m,
      T_range,
      points,
    };
  } catch (error) {
    const message = error instanceof Error ? error.message : 'Unknown WASM error';
    return errorResponse(requestId, 'HANDLER_ERROR', `Boys sweep failed: ${message}`);
  }
}

/**
 * Handle boys_eval_all request (multi-order evaluation).
 *
 * Computes F_0(T) through F_mMax(T) at a single T value.
 * This is efficient because the WASM function computes
 * intermediate orders during the recurrence/series evaluation.
 */
export function handleBoysEvalAll(
  request: BoysEvalAllRequest
): BoysEvalAllResult | ErrorResponse {
  const { requestId, mMax, T } = request;

  // Validate parameters
  if (mMax < 0 || !Number.isInteger(mMax)) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid mMax=${mMax}: must be a non-negative integer`);
  }
  if (T < 0) {
    return errorResponse(requestId, 'INVALID_PARAMS', `Invalid argument T=${T}: must be >= 0`);
  }

  try {
    const results = boys_eval_all(mMax, T) as BoysEvalResult[];
    return {
      results,
      mMax,
      T,
    };
  } catch (error) {
    const message = error instanceof Error ? error.message : 'Unknown WASM error';
    return errorResponse(requestId, 'HANDLER_ERROR', `Boys eval_all failed: ${message}`);
  }
}
