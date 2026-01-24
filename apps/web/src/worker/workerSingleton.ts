/**
 * Worker Singleton
 *
 * Provides a shared Web Worker instance for all compute operations.
 * This ensures that custom systems registered during integral computation
 * are available for subsequent SCF calculations.
 *
 * The worker is created lazily on first access and shared across all
 * hooks (useWorker, useScf, useIntegralCompute, etc.).
 *
 * @module worker/workerSingleton
 */

import type {
  WorkerResponse,
  WorkerRequest,
  RequestId,
  WorkerProgress,
  WorkerErrorCode,
} from './protocol';
import { createRequestId } from './protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * Pending request handlers.
 */
interface PendingRequest<T = unknown> {
  /** Resolve the promise with result data */
  resolve: (value: T) => void;
  /** Reject the promise with an error */
  reject: (error: WorkerError) => void;
  /** Optional progress callback */
  onProgress?: (progress: WorkerProgress) => void;
}

/**
 * Error thrown by worker operations.
 */
export class WorkerError extends Error {
  /** Machine-readable error code */
  readonly code: WorkerErrorCode;

  constructor(code: WorkerErrorCode, message: string) {
    super(message);
    this.name = 'WorkerError';
    this.code = code;
  }
}

/**
 * Listener callback for state changes.
 */
export type StateChangeListener = (isReady: boolean, error: Error | null) => void;

// ============================================================================
// Singleton State
// ============================================================================

/** The shared worker instance */
let worker: Worker | null = null;

/** Pending requests map */
const pendingRequests = new Map<RequestId, PendingRequest>();

/** Whether the worker is ready (WASM initialized) */
let isReady = false;

/** Worker initialization error, if any */
let initError: Error | null = null;

/** Whether initialization is in progress */
let isInitializing = false;

/** Promise for initialization (for async awaiting) */
let initPromise: Promise<void> | null = null;

/** Listeners for state changes */
const stateListeners = new Set<StateChangeListener>();

// ============================================================================
// State Management
// ============================================================================

/**
 * Notify all listeners of state changes.
 */
function notifyStateChange(): void {
  stateListeners.forEach((listener) => {
    try {
      listener(isReady, initError);
    } catch (e) {
      console.error('[WorkerSingleton] Listener error:', e);
    }
  });
}

/**
 * Subscribe to worker state changes.
 *
 * The listener is called immediately with current state and
 * again whenever the state changes.
 *
 * @param listener - Callback receiving (isReady, error)
 * @returns Unsubscribe function
 */
export function subscribeToWorkerState(listener: StateChangeListener): () => void {
  stateListeners.add(listener);
  // Immediately notify with current state
  listener(isReady, initError);
  return () => {
    stateListeners.delete(listener);
  };
}

/**
 * Get current worker state synchronously.
 *
 * @returns Object with isReady and error
 */
export function getWorkerState(): { isReady: boolean; error: Error | null } {
  return { isReady, error: initError };
}

// ============================================================================
// Worker Management
// ============================================================================

/**
 * Initialize the shared worker.
 *
 * This is called automatically on first request, but can be called
 * earlier to pre-warm the worker.
 *
 * @returns Promise that resolves when worker is ready
 */
export function initializeWorker(): Promise<void> {
  // If already ready, resolve immediately
  if (isReady && worker) {
    return Promise.resolve();
  }

  // If initialization failed, reject with cached error
  if (initError) {
    return Promise.reject(initError);
  }

  // If already initializing, return the existing promise
  if (isInitializing && initPromise) {
    return initPromise;
  }

  // Start initialization
  isInitializing = true;

  initPromise = new Promise<void>((resolve, reject) => {
    try {
      // Create worker using Vite's worker import syntax
      worker = new Worker(new URL('./compute.worker.ts', import.meta.url), {
        type: 'module',
      });

      // Handle messages from worker
      worker.onmessage = (event: MessageEvent<WorkerResponse>) => {
        handleWorkerMessage(event.data);
      };

      // Handle worker errors
      worker.onerror = (event: ErrorEvent) => {
        console.error('[WorkerSingleton] Worker error:', event.message);
        const error = new Error(`Worker error: ${event.message}`);

        if (!isReady) {
          // Initialization error
          initError = error;
          isInitializing = false;
          notifyStateChange();
          reject(error);
        } else {
          // Runtime error - reject all pending requests
          pendingRequests.forEach((pending, requestId) => {
            pending.reject(new WorkerError('HANDLER_ERROR', `Worker error: ${event.message}`));
            pendingRequests.delete(requestId);
          });
        }
      };

      // Handle worker message errors
      worker.onmessageerror = (event: MessageEvent) => {
        console.error('[WorkerSingleton] Message error:', event);
        const error = new Error('Worker message serialization error');

        if (!isReady) {
          initError = error;
          isInitializing = false;
          notifyStateChange();
          reject(error);
        }
      };

      // Send initial ping to verify worker is ready
      const initRequestId = createRequestId();
      pendingRequests.set(initRequestId, {
        resolve: () => {
          isReady = true;
          isInitializing = false;
          notifyStateChange();
          resolve();
        },
        reject: (err) => {
          initError = err;
          isInitializing = false;
          notifyStateChange();
          reject(err);
        },
      });

      worker.postMessage({ type: 'ping', requestId: initRequestId });
    } catch (e) {
      const error = e instanceof Error ? e : new Error('Failed to create worker');
      initError = error;
      isInitializing = false;
      notifyStateChange();
      reject(error);
    }
  });

  return initPromise;
}

/**
 * Handle incoming messages from the worker.
 */
function handleWorkerMessage(response: WorkerResponse): void {
  const pending = pendingRequests.get(response.requestId);

  if (!pending) {
    // This can happen if a request was cancelled or timed out
    console.warn('[WorkerSingleton] Received response for unknown request:', response.requestId);
    return;
  }

  switch (response.type) {
    case 'pong':
      pending.resolve(response.wasmVersion);
      pendingRequests.delete(response.requestId);
      break;

    case 'result':
      pending.resolve(response.data);
      pendingRequests.delete(response.requestId);
      break;

    case 'error':
      pending.reject(new WorkerError(response.code, response.message));
      pendingRequests.delete(response.requestId);
      break;

    case 'progress':
      // Progress update - call callback but keep pending
      pending.onProgress?.(response.progress);
      // Don't delete from pending - more messages may come
      break;
  }
}

// ============================================================================
// Public API
// ============================================================================

/**
 * Request options for sendToWorker().
 */
export interface SendOptions {
  /** Callback for progress updates (e.g., SCF iterations) */
  onProgress?: (progress: WorkerProgress) => void;
}

/**
 * Send a request to the shared worker.
 *
 * Automatically initializes the worker if not already done.
 *
 * @param request - The request to send (without requestId)
 * @param options - Optional configuration including progress callback
 * @returns Promise that resolves with the result data
 * @throws WorkerError if the worker returns an error response
 */
export async function sendToWorker<T>(
  request: Omit<WorkerRequest, 'requestId'>,
  options?: SendOptions
): Promise<T> {
  // Ensure worker is initialized
  await initializeWorker();

  if (!worker) {
    throw new WorkerError('WORKER_NOT_READY', 'Worker not initialized');
  }

  return new Promise<T>((resolve, reject) => {
    const requestId = createRequestId();

    // Store pending request handlers
    pendingRequests.set(requestId, {
      resolve: resolve as (value: unknown) => void,
      reject,
      onProgress: options?.onProgress,
    });

    // Send message to worker
    const fullRequest = { ...request, requestId } as WorkerRequest;
    worker!.postMessage(fullRequest);
  });
}

/**
 * Cancel a running request on the shared worker.
 *
 * @param targetRequestId - The request ID to cancel
 */
export function cancelWorkerRequest(targetRequestId: RequestId): void {
  if (!worker) {
    console.warn('[WorkerSingleton] Cannot cancel: worker not initialized');
    return;
  }

  const requestId = createRequestId();
  worker.postMessage({
    type: 'cancel',
    requestId,
    targetRequestId,
  });
}

/**
 * Ping the worker and get WASM version.
 *
 * @returns Promise that resolves with WASM version string
 */
export async function pingWorker(): Promise<string> {
  return sendToWorker<string>({ type: 'ping' });
}

/**
 * Terminate the shared worker.
 *
 * This should only be called during cleanup (e.g., hot module reload).
 * After calling this, the worker will be re-created on next request.
 */
export function terminateWorker(): void {
  if (worker) {
    // Reject all pending requests
    pendingRequests.forEach((pending) => {
      pending.reject(new WorkerError('WORKER_NOT_READY', 'Worker was terminated'));
    });
    pendingRequests.clear();

    worker.terminate();
    worker = null;
    isReady = false;
    isInitializing = false;
    initPromise = null;
    initError = null;
    notifyStateChange();
  }
}

// Re-export types
export type { WorkerProgress, RequestId };
export { createRequestId };
