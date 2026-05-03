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
export { handleBoysEval, handleBoysSweep, handleBoysEvalAll } from './boys';
export { handleRysCompute, handleRysErrorCurve } from './rys';
export { handleScfRun, setScfRunFn } from './scf';
export {
  handleIntegralCompute,
  setComputeIntegralsFn,
  setComputeIntegralsWithOptionsFn,
  setComputeIntegralsWithOptionsAndProgressFn,
} from './integral';
export { handlePesScan, setPesScanFn } from './pes';
export { handlePesScanInternal, setPesScanInternalFn } from './pesInternal';
export { handleMoGrid, setMoGridFn } from './moGrid';
export {
  handleMarchingCubes,
  handleDualMarchingCubes,
  setMarchingCubesFn,
  setDualMarchingCubesFn,
} from './marchingCubes';
export { handleBasisInfo, setGetBasisInfoFn } from './basisInfo';
export { handleRadialProfile, setEvaluateRadialProfileFn } from './radialProfile';
export { handleOverlapDistance, setOverlapVsDistanceFn } from './overlapDistance';
export { handleIntegralMatrices, setComputeIntegralMatricesFn } from './integralMatrices';
export { handleIntegralBreakdown, setIntegralBreakdownFn } from './integralBreakdown';
export { handleFockDecomposition, setFockDecompositionFn } from './fockDecomposition';
export { handleEriDetail, setEriDetailFn } from './eriDetail';
export { handleDensityGrid, setDensityGridFn } from './densityGrid';
export { handleDifferenceDensity, setDifferenceDensityFn } from './differenceDensity';
export { handleOptimizeGeometry, setOptimizeGeometryFn } from './optimizeGeometry';
export { handlePopulation, setComputePopulationFn } from './population';
export { handleFrequency, setComputeFrequenciesFn } from './hessian';

// Re-export types for convenience
export type { ProgressCallback, AbortChecker, ScfRunFn } from './scf';
export type {
  ComputeIntegralsFn,
  ComputeIntegralsWithOptionsFn,
  ComputeIntegralsWithOptionsAndProgressFn,
} from './integral';
export type { PesScanFn } from './pes';
export type { PesScanInternalFn } from './pesInternal';
export type { MoGridFn } from './moGrid';
export type { MarchingCubesFn, DualMarchingCubesFn } from './marchingCubes';
export type { GetBasisInfoFn } from './basisInfo';
export type { EvaluateRadialProfileFn } from './radialProfile';
export type { OverlapVsDistanceFn } from './overlapDistance';
export type { ComputeIntegralMatricesFn } from './integralMatrices';
export type { IntegralBreakdownFn } from './integralBreakdown';
export type { FockDecompositionFn } from './fockDecomposition';
export type { EriDetailFn } from './eriDetail';
export type { DensityGridFn } from './densityGrid';
export type { DifferenceDensityFn } from './differenceDensity';
export type { OptimizeGeometryFn } from './optimizeGeometry';
export type { ComputePopulationFn } from './population';
export type { ComputeFrequenciesFn } from './hessian';
