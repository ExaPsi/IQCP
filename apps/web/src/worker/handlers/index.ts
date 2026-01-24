/**
 * Handler registry for worker message processing.
 *
 * This module exports all handlers and provides a unified interface
 * for the worker to route messages to the appropriate handler.
 *
 * @module worker/handlers
 */

// Export individual handlers
export { handlePing } from './ping';
export type { ThreadingStatus } from './ping';
export { handleBoysEval, handleBoysSweep } from './boys';
export { handleRysCompute, handleRysErrorCurve } from './rys';
export { handleScfRun, setScfRunFn } from './scf';
export {
  handleIntegralCompute,
  setComputeIntegralsFn,
  setComputeIntegralsWithOptionsFn,
  setComputeIntegralsWithOptionsAndProgressFn,
} from './integral';

// Re-export types for convenience
export type { ProgressCallback, AbortChecker, ScfRunFn } from './scf';
export type {
  ComputeIntegralsFn,
  ComputeIntegralsWithOptionsFn,
  ComputeIntegralsWithOptionsAndProgressFn,
} from './integral';
