/**
 * WASM module initialization and re-exports.
 *
 * This module provides a safe way to initialize the WASM module
 * and access its exported functions.
 *
 * @example
 * ```typescript
 * import { initWasm, version, test_compute } from '@/lib/wasm';
 *
 * await initWasm();
 * console.log('Version:', version());
 * console.log('Test:', test_compute(21));
 * ```
 */

import init, { version, test_compute } from '../wasm/qc_wasm';

/** Track initialization state */
let initialized = false;

/**
 * Initialize the WASM module.
 *
 * Safe to call multiple times; only initializes once.
 * Must be called before using any WASM functions.
 *
 * @example
 * ```typescript
 * await initWasm();
 * const v = version();
 * ```
 */
export async function initWasm(): Promise<void> {
  if (!initialized) {
    await init();
    initialized = true;
    console.log(`WASM module loaded: qc-wasm v${version()}`);
  }
}

/**
 * Check if WASM module is initialized.
 */
export function isWasmInitialized(): boolean {
  return initialized;
}

/**
 * Get WASM module version.
 *
 * @throws Error if WASM not initialized
 */
export function getWasmVersion(): string {
  if (!initialized) {
    throw new Error('WASM not initialized. Call initWasm() first.');
  }
  return version();
}

// Re-export WASM functions for convenience
export { version, test_compute };

/**
 * TestResult type matching the Rust struct.
 * This is the return type of test_compute().
 */
export interface TestResult {
  /** The input value that was provided */
  input: number;
  /** The computed output (input * 2) */
  output: number;
  /** A human-readable message describing the computation */
  message: string;
}
