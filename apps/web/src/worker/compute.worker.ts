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
} from './protocol';
import { assertNever } from './protocol';
import {
  handlePing,
  handleBoysEval,
  handleBoysSweep,
  handleRysCompute,
  handleRysErrorCurve,
  handleScfRun,
  setScfRunFn,
  handleIntegralCompute,
  setComputeIntegralsFn,
  setComputeIntegralsWithOptionsFn,
  setComputeIntegralsWithOptionsAndProgressFn,
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
 * Note: We don't validate the message structure here because TypeScript
 * guarantees the shape at the call site. If needed, use `isWorkerRequest`
 * from protocol.ts for runtime validation.
 */
self.onmessage = (event: MessageEvent<WorkerRequest>) => {
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
