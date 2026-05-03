/**
 * Compute Web Worker
 *
 * This worker handles all heavy computation (Boys functions, Rys quadrature,
 * SCF iterations) off the main thread to keep the UI responsive.
 *
 * The worker initializes the WASM module on startup and routes incoming
 * messages to the appropriate handler.
 *
 * @module worker/compute.worker
 */

import init, { version, scf_run, has_threading_support } from '../wasm/qc_wasm';
// Note: compute_integrals will be available after US-029-T3 is implemented
// For now, we import it dynamically to avoid TypeScript errors
import * as wasmModule from '../wasm/qc_wasm';

// Type assertion for init_thread_pool which is only available with parallel feature
type InitThreadPoolFn = (numThreads: number) => Promise<void>;
import type {
  WorkerRequest,
  WorkerResponse,
  RequestId,
  ErrorResponse,
  ProgressResponse,
  WorkerProgress,
  ResultResponse,
  CancelResult,
  KsScfResult,
} from './protocol';
import { assertNever } from './protocol';
import {
  handlePing,
  handleBoysEval,
  handleBoysSweep,
  handleBoysEvalAll,
  handleRysCompute,
  handleRysErrorCurve,
  handleScfRun,
  setScfRunFn,
  handleIntegralCompute,
  setComputeIntegralsFn,
  setComputeIntegralsWithOptionsFn,
  setComputeIntegralsWithOptionsAndProgressFn,
  handlePesScan,
  setPesScanFn,
  handlePesScanInternal,
  setPesScanInternalFn,
  handleMoGrid,
  setMoGridFn,
  handleMarchingCubes,
  handleDualMarchingCubes,
  setMarchingCubesFn,
  setDualMarchingCubesFn,
  handleBasisInfo,
  setGetBasisInfoFn,
  handleRadialProfile,
  setEvaluateRadialProfileFn,
  handleOverlapDistance,
  setOverlapVsDistanceFn,
  handleIntegralMatrices,
  setComputeIntegralMatricesFn,
  handleIntegralBreakdown,
  setIntegralBreakdownFn,
  handleFockDecomposition,
  setFockDecompositionFn,
  handleEriDetail,
  setEriDetailFn,
  handleDensityGrid,
  setDensityGridFn,
  handleDifferenceDensity,
  setDifferenceDensityFn,
  handleOptimizeGeometry,
  setOptimizeGeometryFn,
  handlePopulation,
  setComputePopulationFn,
  handleFrequency,
  setComputeFrequenciesFn,
} from './handlers';
import type { ThreadingStatus } from './handlers';
import { registerCustomSystem } from './presets';

// ============================================================================
// Worker State
// ============================================================================

/**
 * Worker initialization state.
 */
type WorkerState =
  | { status: 'pending' }
  | { status: 'ready'; wasmVersion: string; threadsAvailable: boolean; numThreads: number }
  | { status: 'error'; message: string };

/** Current worker state */
let state: WorkerState = { status: 'pending' };

/** Abort flags for cancellable operations, keyed by request ID */
const abortFlags = new Map<RequestId, boolean>();

/** Whether WASM initialization has completed (success or failure) */
let initializationComplete = false;

/** Whether the spectra WASM module (qc-wasm-spectra) has been lazy-loaded */
let spectraModuleLoaded = false;

/** Queue of messages received before initialization completed */
const pendingMessages: WorkerRequest[] = [];

// ============================================================================
// Response Helpers
// ============================================================================

/**
 * Post a response to the main thread.
 */
function respond(response: WorkerResponse): void {
  self.postMessage(response);
}

/**
 * Create an error response.
 */
function errorResponse(
  requestId: RequestId,
  code: ErrorResponse['code'],
  message: string
): ErrorResponse {
  return { type: 'error', requestId, code, message };
}

/**
 * Create a result response.
 */
function resultResponse(requestId: RequestId, data: unknown): ResultResponse {
  return { type: 'result', requestId, data };
}

/**
 * Create a progress response.
 */
function progressResponse(requestId: RequestId, progress: WorkerProgress): ProgressResponse {
  return { type: 'progress', requestId, progress };
}

// ============================================================================
// Initialization
// ============================================================================

/**
 * Initialize the WASM module.
 *
 * This is called automatically when the worker starts.
 * Sets up the WASM module, initializes thread pool if available,
 * and wires up handler dependencies.
 */
async function initialize(): Promise<void> {
  try {
    await init();
    const wasmVersion = version();

    // Wire up the SCF handler with the WASM function
    setScfRunFn(scf_run);

    // Wire up the integral compute handler with the WASM functions
    // Priority: compute_integrals_with_options_and_progress (best) > others
    //
    // Available functions in qc-wasm:
    //   - compute_integrals: basic (no options, no progress)
    //   - compute_integrals_with_progress: has progress callback, no options
    //   - compute_integrals_with_options: has options (spherical), no progress
    //   - compute_integrals_with_options_and_progress: has both (preferred)
    const computeIntegralsWithOptionsAndProgress = (wasmModule as Record<string, unknown>)[
      'compute_integrals_with_options_and_progress'
    ];
    const computeIntegralsWithProgress = (wasmModule as Record<string, unknown>)[
      'compute_integrals_with_progress'
    ];
    const computeIntegralsWithOptions = (wasmModule as Record<string, unknown>)[
      'compute_integrals_with_options'
    ];
    const computeIntegrals = (wasmModule as Record<string, unknown>)['compute_integrals'];

    // Set up compute_integrals_with_options_and_progress (best: has both options and progress)
    if (typeof computeIntegralsWithOptionsAndProgress === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setComputeIntegralsWithOptionsAndProgressFn(computeIntegralsWithOptionsAndProgress as any);
    }

    // Set up compute_integrals_with_progress (for Cartesian with progress fallback)
    if (typeof computeIntegralsWithProgress === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setComputeIntegralsFn(computeIntegralsWithProgress as any);
    } else if (typeof computeIntegrals === 'function') {
      // Fallback to non-progress version (wrap to match expected signature)
      setComputeIntegralsFn(
        ((geometryJson: string, basisName: string) =>
          // eslint-disable-next-line @typescript-eslint/no-explicit-any
          (computeIntegrals as (a: string, b: string) => any)(geometryJson, basisName)) as any
      );
    }

    // Set up compute_integrals_with_options for spherical harmonic support (no progress fallback)
    if (typeof computeIntegralsWithOptions === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setComputeIntegralsWithOptionsFn(computeIntegralsWithOptions as any);
    }

    // Wire up the PES scan handler with the WASM function
    const pesScanFn = (wasmModule as Record<string, unknown>)['pes_scan'];
    if (typeof pesScanFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setPesScanFn(pesScanFn as any);
    } else {
      console.warn('[Worker] WASM function "pes_scan" not found — PES scan will be unavailable. Rebuild WASM?');
    }

    // Wire up the internal coordinate PES scan handler with the WASM function (US-081)
    const pesScanInternalFn = (wasmModule as Record<string, unknown>)['pes_scan_internal'];
    if (typeof pesScanInternalFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setPesScanInternalFn(pesScanInternalFn as any);
    } else {
      console.warn('[Worker] WASM function "pes_scan_internal" not found — internal coordinate PES scan will be unavailable. Rebuild WASM?');
    }

    // Wire up the MO grid evaluation handler with the WASM function
    const moGridFn = (wasmModule as Record<string, unknown>)['evaluate_mo_grid'];
    if (typeof moGridFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setMoGridFn(moGridFn as any);
    } else {
      console.warn('[Worker] WASM function "evaluate_mo_grid" not found — orbital visualization will be unavailable. Rebuild WASM?');
    }

    // Wire up the marching cubes handler with the WASM function
    const marchingCubesFn = (wasmModule as Record<string, unknown>)['marching_cubes'];
    if (typeof marchingCubesFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setMarchingCubesFn(marchingCubesFn as any);
    } else {
      console.warn('[Worker] WASM function "marching_cubes" not found — isosurface extraction will be unavailable. Rebuild WASM?');
    }

    // Wire up the dual marching cubes handler with the WASM function
    const dualMarchingCubesFn = (wasmModule as Record<string, unknown>)['dual_marching_cubes'];
    if (typeof dualMarchingCubesFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setDualMarchingCubesFn(dualMarchingCubesFn as any);
    } else {
      console.warn('[Worker] WASM function "dual_marching_cubes" not found — dual isosurface extraction will be unavailable. Rebuild WASM?');
    }

    // Wire up the basis info handler with the WASM function
    const getBasisInfoFn = (wasmModule as Record<string, unknown>)['get_basis_info'];
    if (typeof getBasisInfoFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setGetBasisInfoFn(getBasisInfoFn as any);
    } else {
      console.warn('[Worker] WASM function "get_basis_info" not found — basis info queries will be unavailable. Rebuild WASM?');
    }

    // Wire up the radial profile handler with the WASM function
    const evaluateRadialProfileFn = (wasmModule as Record<string, unknown>)['evaluate_radial_profile'];
    if (typeof evaluateRadialProfileFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setEvaluateRadialProfileFn(evaluateRadialProfileFn as any);
    } else {
      console.warn('[Worker] WASM function "evaluate_radial_profile" not found — radial profile evaluation will be unavailable. Rebuild WASM?');
    }

    // Wire up the overlap vs distance handler with the WASM function
    const overlapVsDistanceFn = (wasmModule as Record<string, unknown>)['overlap_vs_distance'];
    if (typeof overlapVsDistanceFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setOverlapVsDistanceFn(overlapVsDistanceFn as any);
    } else {
      console.warn('[Worker] WASM function "overlap_vs_distance" not found — overlap distance plot will be unavailable. Rebuild WASM?');
    }

    // Wire up the integral matrices handler with the WASM function
    const computeIntegralMatricesFn = (wasmModule as Record<string, unknown>)['compute_integral_matrices'];
    if (typeof computeIntegralMatricesFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setComputeIntegralMatricesFn(computeIntegralMatricesFn as any);
    } else {
      console.warn('[Worker] WASM function "compute_integral_matrices" not found — integral matrices will be unavailable. Rebuild WASM?');
    }

    // Wire up integral_with_breakdown (US-056)
    const integralBreakdownFn = (wasmModule as Record<string, unknown>)['integral_with_breakdown'];
    if (typeof integralBreakdownFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setIntegralBreakdownFn(integralBreakdownFn as any);
    } else {
      console.warn('[Worker] WASM function "integral_with_breakdown" not found — primitive breakdown will be unavailable. Rebuild WASM?');
    }

    // Wire up fock_decomposition (US-058)
    const fockDecompositionFn = (wasmModule as Record<string, unknown>)['fock_decomposition'];
    if (typeof fockDecompositionFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setFockDecompositionFn(fockDecompositionFn as any);
    } else {
      console.warn('[Worker] WASM function "fock_decomposition" not found — Fock build tracing will be unavailable. Rebuild WASM?');
    }

    // Wire up eri_detail (US-059)
    const eriDetailFn = (wasmModule as Record<string, unknown>)['eri_detail'];
    if (typeof eriDetailFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setEriDetailFn(eriDetailFn as any);
    } else {
      console.warn('[Worker] WASM function "eri_detail" not found — ERI browser will be unavailable. Rebuild WASM?');
    }

    // Wire up evaluate_density_grid (US-061)
    const evaluateDensityGridFn = (wasmModule as Record<string, unknown>)['evaluate_density_grid'];
    if (typeof evaluateDensityGridFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setDensityGridFn(evaluateDensityGridFn as any);
    } else {
      console.warn('[Worker] WASM function "evaluate_density_grid" not found — density grid evaluation will be unavailable. Rebuild WASM?');
    }

    // Wire up compute_difference_density (US-063)
    const computeDifferenceDensityFn = (wasmModule as Record<string, unknown>)['compute_difference_density'];
    if (typeof computeDifferenceDensityFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setDifferenceDensityFn(computeDifferenceDensityFn as any);
    } else {
      console.warn('[Worker] WASM function "compute_difference_density" not found — difference density will be unavailable. Rebuild WASM?');
    }

    // Wire up optimize_geometry (US-075)
    const optimizeGeometryFn = (wasmModule as Record<string, unknown>)['optimize_geometry'];
    if (typeof optimizeGeometryFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setOptimizeGeometryFn(optimizeGeometryFn as any);
    } else {
      console.warn('[Worker] WASM function "optimize_geometry" not found — geometry optimization will be unavailable. Rebuild WASM?');
    }

    // Wire compute_population (US-076)
    const computePopulationFn = (wasmModule as Record<string, unknown>)['compute_population'];
    if (typeof computePopulationFn === 'function') {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      setComputePopulationFn(computePopulationFn as any);
    } else {
      console.warn('[Worker] WASM function "compute_population" not found — population analysis will be unavailable. Rebuild WASM?');
    }

    // NOTE: compute_frequencies has been moved to qc-wasm-spectra and is
    // lazy-loaded on first frequency request (see case 'frequency' below).

    // Check if threading is available and initialize thread pool
    let threadsAvailable = false;
    let numThreads = 0;

    // First check if SharedArrayBuffer is available (requires COOP/COEP headers)
    const hasSharedArrayBuffer = typeof SharedArrayBuffer !== 'undefined';

    try {
      // has_threading_support() returns true if WASM was built with parallel feature
      threadsAvailable = has_threading_support();

      if (threadsAvailable && hasSharedArrayBuffer) {
        // Get the init_thread_pool function from the module
        const initThreadPool = (wasmModule as Record<string, unknown>)[
          'initThreadPool'
        ] as InitThreadPoolFn | undefined;

        if (initThreadPool) {
          // Use hardware concurrency or default to 4 threads
          numThreads = navigator.hardwareConcurrency || 4;

          try {
            await initThreadPool(numThreads);
          } catch {
            // Thread pool init can fail if SharedArrayBuffer is not available
            // (missing COOP/COEP headers) or other browser restrictions
            threadsAvailable = false;
            numThreads = 0;
          }
        } else {
          threadsAvailable = false;
        }
      } else {
        threadsAvailable = false;
      }
    } catch {
      // has_threading_support() itself failed - treat as no threading
      threadsAvailable = false;
      numThreads = 0;
    }

    state = { status: 'ready', wasmVersion, threadsAvailable, numThreads };
  } catch (error) {
    const message = error instanceof Error ? error.message : 'Unknown initialization error';
    state = { status: 'error', message };
    console.error('[Worker] WASM initialization failed:', message);
  }

  // Mark initialization as complete and drain any queued messages
  initializationComplete = true;
  while (pendingMessages.length > 0) {
    const queued = pendingMessages.shift()!;
    void handleMessage(queued);
  }
}

// ============================================================================
// Message Handler
// ============================================================================

/**
 * Handle incoming messages from the main thread.
 */
async function handleMessage(request: WorkerRequest): Promise<void> {
  const { requestId } = request;

  // Ping is special - works even if WASM not ready
  if (request.type === 'ping') {
    const wasmVersion = state.status === 'ready' ? state.wasmVersion : 'not_initialized';
    const threadingStatus: ThreadingStatus =
      state.status === 'ready'
        ? { threadsAvailable: state.threadsAvailable, numThreads: state.numThreads }
        : { threadsAvailable: false, numThreads: 0 };
    respond(handlePing(request, wasmVersion, threadingStatus));
    return;
  }

  // All other requests require WASM to be ready
  if (state.status !== 'ready') {
    respond(
      errorResponse(
        requestId,
        'WORKER_NOT_READY',
        state.status === 'pending'
          ? 'Worker is still initializing WASM module'
          : `WASM initialization failed: ${state.message}`
      )
    );
    return;
  }

  // Clear any previous abort flag for this request
  abortFlags.delete(requestId);

  try {
    switch (request.type) {
      case 'boys_eval': {
        const result = handleBoysEval(request);
        if ('code' in result) {
          respond(result); // It's an error response
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'boys_sweep': {
        const result = handleBoysSweep(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'boys_eval_all': {
        const result = handleBoysEvalAll(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'rys_compute': {
        const result = handleRysCompute(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'rys_error_curve': {
        const result = handleRysErrorCurve(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'scf_run': {
        // Set up progress callback
        const onProgress = (progress: WorkerProgress): void => {
          respond(progressResponse(requestId, progress));
        };

        // Set up abort checker
        const isAborted = (): boolean => {
          return abortFlags.get(requestId) === true;
        };

        const result = handleScfRun(request, onProgress, isAborted);

        // Clean up abort flag
        abortFlags.delete(requestId);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'integral_compute': {
        // Set up progress callback
        const onProgress = (progress: WorkerProgress): void => {
          respond(progressResponse(requestId, progress));
        };

        // Set up abort checker
        const isAborted = (): boolean => {
          return abortFlags.get(requestId) === true;
        };

        const result = handleIntegralCompute(request, onProgress, isAborted);

        // Clean up abort flag
        abortFlags.delete(requestId);

        if ('code' in result) {
          respond(result);
        } else {
          // Register the computed system in the cache for SCF to find
          registerCustomSystem(result.systemId, result);
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'pes_scan': {
        // Set up progress callback
        const onProgress = (progress: WorkerProgress): void => {
          respond(progressResponse(requestId, progress));
        };

        // Set up abort checker
        const isAborted = (): boolean => {
          return abortFlags.get(requestId) === true;
        };

        const result = handlePesScan(request, onProgress, isAborted);

        // Clean up abort flag
        abortFlags.delete(requestId);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'pes_scan_internal': {
        // Internal coordinate PES scan (US-081)
        const onProgress = (progress: WorkerProgress): void => {
          respond(progressResponse(requestId, progress));
        };

        const isAborted = (): boolean => {
          return abortFlags.get(requestId) === true;
        };

        const result = handlePesScanInternal(request, onProgress, isAborted);

        abortFlags.delete(requestId);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'mo_grid': {
        const result = handleMoGrid(request);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'marching_cubes': {
        const result = handleMarchingCubes(request);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'dual_marching_cubes': {
        const result = handleDualMarchingCubes(request);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'basis_info': {
        const result = handleBasisInfo(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'radial_profile': {
        const result = handleRadialProfile(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'overlap_distance': {
        const result = handleOverlapDistance(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'integral_matrices': {
        const result = handleIntegralMatrices(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'integral_breakdown': {
        const result = handleIntegralBreakdown(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'fock_decomposition': {
        const result = handleFockDecomposition(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'eri_detail': {
        const result = handleEriDetail(request);
        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'density_grid': {
        const result = handleDensityGrid(request);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'difference_density': {
        const result = handleDifferenceDensity(request);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'optimize_geometry': {
        // Geometry optimization (US-075)
        const onProgress = (progress: WorkerProgress): void => {
          respond(progressResponse(requestId, progress));
        };

        const isAborted = (): boolean => {
          return abortFlags.get(requestId) === true;
        };

        const result = handleOptimizeGeometry(request, onProgress, isAborted);

        abortFlags.delete(requestId);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'population_analysis': {
        // Mulliken/Lowdin population analysis (US-076)
        const result = handlePopulation(request);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'frequency': {
        // Frequency analysis (US-101): Hessian -> normal modes -> IR/Raman ->
        // RRHO thermochemistry -> broadened spectra. Streams phase-granular
        // progress and supports cooperative cancellation.
        //
        // The spectra WASM module (qc-wasm-spectra) is lazy-loaded on the
        // first frequency request to keep the initial WASM bundle < 500 KB.
        const onProgress = (progress: WorkerProgress): void => {
          respond(progressResponse(requestId, progress));
        };

        const isAborted = (): boolean => {
          return abortFlags.get(requestId) === true;
        };

        // Lazy-load spectra WASM on first frequency request
        try {
          if (!spectraModuleLoaded) {
            const spectraModule = await import('../wasm-spectra/qc_wasm_spectra');
            await spectraModule.default(); // init WASM
            // eslint-disable-next-line @typescript-eslint/no-explicit-any
            setComputeFrequenciesFn(spectraModule.compute_frequencies as any);
            spectraModuleLoaded = true;
          }
        } catch (err) {
          const message = err instanceof Error ? err.message : String(err);
          respond(
            errorResponse(
              requestId,
              'HANDLER_ERROR',
              `Failed to load spectra WASM module: ${message}`,
            ),
          );
          break;
        }

        const result = handleFrequency(request, onProgress, isAborted);

        abortFlags.delete(requestId);

        if ('code' in result) {
          respond(result);
        } else {
          respond(resultResponse(requestId, result));
        }
        break;
      }

      case 'ks_scf': {
        // KS-DFT SCF calculation (US-068) with iteration-by-iteration progress streaming
        try {
          const ksScfFn = (wasmModule as Record<string, unknown>)['ks_scf'] as
            | ((input: unknown, progressCallback?: (progress: unknown) => void) => KsScfResult)
            | undefined;
          if (!ksScfFn) {
            respond(
              errorResponse(
                requestId,
                'NOT_IMPLEMENTED',
                'ks_scf WASM function not available (rebuild WASM module)',
              ),
            );
            break;
          }
          const maxIter = request.maxIterations ?? 100;
          const ksInput: Record<string, unknown> = {
              atoms: request.atoms,
              basisName: request.basisName,
              method: request.method,
              convergenceProfile: request.convergenceProfile ?? 'tight',
              maxIterations: maxIter,
              useDiis: request.useDiis ?? true,
              gridQuality: request.gridQuality ?? 'standard',
          };
          if (request.useSpherical) {
            ksInput.useSpherical = true;
          }
          const result = ksScfFn(
            ksInput,
            (progressData: unknown) => {
              const raw = progressData as Record<string, unknown>;
              if (raw.phase === 'integrals') {
                // Integral-phase progress (S, Hcore, ERI, grid)
                respond(
                  progressResponse(requestId, {
                    module: 'scf_integrals' as const,
                    step: String(raw.step),
                    percent: Number(raw.percent),
                    message: String(raw.message),
                    current: 0,
                    total: 100,
                  }),
                );
              } else {
                // SCF iteration progress
                const p = raw as unknown as {
                  iteration: number;
                  energy: number;
                  deltaE: number;
                  rmsDensityChange: number;
                  diisApplied: boolean;
                };
                respond(
                  progressResponse(requestId, {
                    module: 'scf',
                    iteration: p.iteration,
                    energy: p.energy,
                    delta: p.deltaE ?? 0,
                    diisError: p.rmsDensityChange,
                    converged: false,
                    current: p.iteration,
                    total: maxIter,
                    message: `Iteration ${p.iteration}: E = ${p.energy.toFixed(10)} Ha`,
                  }),
                );
              }
            },
          );
          respond(resultResponse(requestId, result));
        } catch (err) {
          const message = err instanceof Error ? err.message : String(err);
          respond(errorResponse(requestId, 'HANDLER_ERROR', `KS-SCF error: ${message}`));
        }
        break;
      }

      case 'cancel': {
        const { targetRequestId } = request;
        // Set abort flag for the target request
        abortFlags.set(targetRequestId, true);

        const cancelResult: CancelResult = {
          cancelled: true,
          targetRequestId,
        };
        respond(resultResponse(requestId, cancelResult));
        break;
      }

      default:
        // This ensures all cases are handled at compile time
        assertNever(request);
    }
  } catch (error) {
    // Catch any unexpected errors from handlers
    const message = error instanceof Error ? error.message : 'Unknown handler error';
    respond(errorResponse(requestId, 'HANDLER_ERROR', message));
  }
}

// ============================================================================
// Worker Entry Point
// ============================================================================

/**
 * Message event handler.
 *
 * Messages received before WASM initialization completes are queued
 * and processed once initialization finishes. This prevents race
 * conditions where the singleton's init ping arrives before WASM
 * functions are wired up, which would cause subsequent requests
 * to fail with "Worker is still initializing WASM module".
 *
 * Note: We don't validate the message structure here because TypeScript
 * guarantees the shape at the call site. If needed, use `isWorkerRequest`
 * from protocol.ts for runtime validation.
 */
self.onmessage = (event: MessageEvent<WorkerRequest>) => {
  if (!initializationComplete) {
    pendingMessages.push(event.data);
    return;
  }
  void handleMessage(event.data);
};

/**
 * Error event handler.
 *
 * Note: In worker context, onerror receives different parameters than in window context.
 */
self.onerror = (message: string | Event) => {
  const errorMessage = typeof message === 'string' ? message : 'Unknown worker error';
  console.error('[Worker] Unhandled error:', errorMessage);
};

/**
 * Unhandled rejection handler.
 */
self.onunhandledrejection = (event: PromiseRejectionEvent) => {
  console.error('[Worker] Unhandled promise rejection:', event.reason);
};

// Initialize WASM on worker start
void initialize();

// Export for testing (if needed)
export { initialize, handleMessage, state };
