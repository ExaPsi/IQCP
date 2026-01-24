/**
 * Ping handler for worker health checks.
 *
 * @module worker/handlers/ping
 */

import type { PingRequest, PongResponse } from '../protocol';

/**
 * Threading status information.
 */
export interface ThreadingStatus {
  /** Whether threading is available in this WASM build */
  threadsAvailable: boolean;
  /** Number of threads initialized (0 if threading not available) */
  numThreads: number;
}

/**
 * Handle ping request.
 *
 * Returns pong response with WASM version string and threading status.
 *
 * @param request - The ping request
 * @param wasmVersion - Version string from WASM module
 * @param threadingStatus - Threading availability and thread count
 * @returns Pong response with WASM version and threading info
 */
export function handlePing(
  request: PingRequest,
  wasmVersion: string,
  threadingStatus: ThreadingStatus
): PongResponse {
  return {
    type: 'pong',
    requestId: request.requestId,
    wasmVersion,
    threadsAvailable: threadingStatus.threadsAvailable,
    numThreads: threadingStatus.numThreads,
  };
}
