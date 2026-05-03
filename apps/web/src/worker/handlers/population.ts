/**
 * Population analysis handler.
 *
 * Calls the WASM compute_population function to perform Mulliken and Lowdin
 * population analysis from density and overlap matrices.
 *
 * @module worker/handlers/population
 * @see US-076 Mulliken/Lowdin Population Analysis
 */

import type {
  PopulationAnalysisRequest,
  PopulationAnalysisResult,
  ErrorResponse,
} from '../protocol';

/**
 * Type for the WASM compute_population function.
 */
export type ComputePopulationFn = (input: unknown) => PopulationAnalysisResult;

// Store reference to the WASM function
let wasmComputePopulation: ComputePopulationFn | null = null;

/**
 * Set the WASM compute_population function reference.
 *
 * Must be called during worker initialization after the WASM module loads.
 */
export function setComputePopulationFn(fn: ComputePopulationFn): void {
  wasmComputePopulation = fn;
}

/**
 * Handle population_analysis request.
 *
 * Computes Mulliken and Lowdin atomic charges from the density and overlap
 * matrices by calling the WASM population analysis function.
 *
 * @param request - The population analysis request
 * @returns Population result or error response
 */
export function handlePopulation(
  request: PopulationAnalysisRequest
): PopulationAnalysisResult | ErrorResponse {
  const { requestId } = request;

  if (!wasmComputePopulation) {
    return {
      type: 'error',
      requestId,
      code: 'WORKER_NOT_READY',
      message: 'Population analysis WASM function not initialized',
    };
  }

  try {
    const input = {
      densityMatrix: request.densityMatrix,
      overlapMatrix: request.overlapMatrix,
      nbf: request.nbf,
      atoms: request.atoms,
    };

    return wasmComputePopulation(input);
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'Population analysis failed',
    };
  }
}
