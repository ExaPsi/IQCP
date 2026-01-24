/**
 * useWorker - React hook for Web Worker communication.
 *
 * Provides a type-safe interface for sending requests to the shared
 * compute worker and receiving responses.
 *
 * Uses a singleton worker instance to ensure that custom systems
 * registered during integral computation are available for subsequent
 * SCF calculations.
 *
 * @module hooks/useWorker
 *
 * @example
 * ```typescript
 * function MyComponent() {
 *   const { isReady, error, send, ping } = useWorker();
 *
 *   useEffect(() => {
 *     if (isReady) {
 *       ping().then(version => console.log('WASM:', version));
 *     }
 *   }, [isReady, ping]);
 *
 *   const handleCompute = async () => {
 *     const result = await send<BoysEvalResult>({ type: 'boys_eval', m: 0, T: 1.0 });
 *     console.log('Result:', result);
 *   };
 *
 *   return <button onClick={handleCompute} disabled={!isReady}>Compute</button>;
 * }
 * ```
 */

import { useState, useEffect, useCallback } from 'react';
import type { WorkerRequest, WorkerProgress, RequestId } from '../worker/protocol';
import { createRequestId } from '../worker/protocol';
import {
  sendToWorker,
  cancelWorkerRequest,
  pingWorker,
  subscribeToWorkerState,
  initializeWorker,
  WorkerError,
  type SendOptions,
} from '../worker/workerSingleton';

// ============================================================================
// Types
// ============================================================================

/**
 * Return type of useWorker hook.
 */
export interface UseWorkerResult {
  /**
   * Whether the worker is ready to receive requests.
   * The worker is ready once it has successfully initialized WASM.
   */
  isReady: boolean;

  /**
   * Error from worker initialization or fatal errors.
   * null if no error has occurred.
   */
  error: Error | null;

  /**
   * Send a request to the worker.
   *
   * @param request - The request to send (without requestId)
   * @param options - Optional configuration including progress callback
   * @returns Promise that resolves with the result data
   * @throws WorkerError if the worker returns an error response
   *
   * @example
   * ```typescript
   * const result = await send<BoysEvalResult>(
   *   { type: 'boys_eval', m: 0, T: 1.0 }
   * );
   * ```
   */
  send: <T>(request: Omit<WorkerRequest, 'requestId'>, options?: SendOptions) => Promise<T>;

  /**
   * Cancel a running request.
   *
   * @param targetRequestId - The request ID to cancel
   *
   * @example
   * ```typescript
   * // Store request ID when sending
   * const requestId = createRequestId();
   * send({ type: 'scf_run', ... });
   *
   * // Later, cancel if needed
   * cancel(requestId);
   * ```
   */
  cancel: (targetRequestId: RequestId) => void;

  /**
   * Ping the worker and get WASM version.
   *
   * This is a convenience method for verifying worker health.
   *
   * @returns Promise that resolves with WASM version string
   *
   * @example
   * ```typescript
   * const version = await ping();
   * console.log(`WASM version: ${version}`);
   * ```
   */
  ping: () => Promise<string>;
}

// ============================================================================
// Hook Implementation
// ============================================================================

/**
 * React hook for communicating with the shared compute Web Worker.
 *
 * Uses a singleton worker instance that is shared across all components.
 * This ensures that custom systems registered during integral computation
 * are available for subsequent SCF calculations.
 *
 * @returns Object with send, cancel, ping methods and isReady/error state
 */
export function useWorker(): UseWorkerResult {
  // Track worker state
  const [isReady, setIsReady] = useState(false);
  const [error, setError] = useState<Error | null>(null);

  // ============================================================================
  // Subscribe to Worker State
  // ============================================================================

  useEffect(() => {
    // Subscribe to state changes from the singleton
    const unsubscribe = subscribeToWorkerState((ready, err) => {
      setIsReady(ready);
      setError(err);
    });

    // Eagerly initialize the worker
    initializeWorker().catch((err) => {
      console.error('[useWorker] Failed to initialize worker:', err);
    });

    return unsubscribe;
  }, []);

  // ============================================================================
  // Send Method
  // ============================================================================

  const send = useCallback(
    <T>(request: Omit<WorkerRequest, 'requestId'>, options?: SendOptions): Promise<T> => {
      return sendToWorker<T>(request, options);
    },
    []
  );

  // ============================================================================
  // Cancel Method
  // ============================================================================

  const cancel = useCallback((targetRequestId: RequestId): void => {
    cancelWorkerRequest(targetRequestId);
  }, []);

  // ============================================================================
  // Ping Method
  // ============================================================================

  const ping = useCallback((): Promise<string> => {
    return pingWorker();
  }, []);

  // ============================================================================
  // Return Value
  // ============================================================================

  return {
    isReady,
    error,
    send,
    cancel,
    ping,
  };
}

// Re-export types for convenience
export type { WorkerProgress, RequestId, SendOptions };
export { createRequestId, WorkerError };
