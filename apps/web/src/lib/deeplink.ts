/**
 * Deep Link Encoding and Decoding
 *
 * Provides URL-safe serialization of application state for deep linking.
 * Uses lz-string compression to minimize URL length while maintaining
 * deterministic encoding.
 *
 * Encoding pipeline:
 * 1. JSON.stringify (deterministic for same input)
 * 2. lz-string compressToEncodedURIComponent (URL-safe)
 *
 * Decoding pipeline:
 * 1. lz-string decompressFromEncodedURIComponent
 * 2. JSON.parse
 * 3. Zod validation
 *
 * @module deeplink
 */

import LZString from 'lz-string';
import { RunStateV1Schema } from '../types/run-state.schema';
import type { RunStateV1 } from '../types/run-state';


/**
 * Encode a RunStateV1 object to a URL-safe string.
 *
 * The encoding process:
 * 1. Serializes the state to JSON using JSON.stringify
 * 2. Compresses the JSON using lz-string's URL-safe encoding
 *
 * The output is deterministic for the same input, making it
 * suitable for comparison and caching.
 *
 * @param state - The state object to encode
 * @returns URL-safe encoded string suitable for use in query parameters
 *
 * @example
 * ```typescript
 * import { encodeRunState } from '@/lib/deeplink';
 * import { DEFAULT_BOYS_STATE } from '@/types/run-state';
 *
 * const encoded = encodeRunState(DEFAULT_BOYS_STATE);
 * // Use in URL: `?run=${encoded}`
 * console.log('Encoded length:', encoded.length); // Typically ~150 chars
 * ```
 */
export function encodeRunState(state: RunStateV1): string {
  const json = JSON.stringify(state);
  const compressed = LZString.compressToEncodedURIComponent(json);
  return compressed;
}

/**
 * Decode a URL-safe string to a RunStateV1 object.
 *
 * The decoding process:
 * 1. Decompresses the string using lz-string
 * 2. Parses the resulting JSON
 * 3. Validates the object against the RunStateV1 Zod schema
 *
 * Returns null if any step fails (does not throw). This allows
 * callers to gracefully handle invalid URLs without try-catch blocks.
 *
 * @param encoded - The URL-safe encoded string, or null
 * @returns Validated RunStateV1 object, or null if invalid/missing
 *
 * @example
 * ```typescript
 * import { decodeRunState } from '@/lib/deeplink';
 *
 * const params = new URLSearchParams(window.location.search);
 * const encoded = params.get('run');
 * const state = decodeRunState(encoded);
 *
 * if (state) {
 *   // Apply the decoded state
 *   applyState(state);
 * } else {
 *   // Show error or use default state
 *   showInvalidLinkError();
 * }
 * ```
 */
export function decodeRunState(encoded: string | null): RunStateV1 | null {
  if (!encoded) {
    return null;
  }

  try {
    // Step 1: Restore + characters that browsers convert to spaces
    // In URL query strings, + is interpreted as space by URLSearchParams
    // lz-string's encodeURIComponent output uses + as a valid character
    // So we need to convert any spaces back to + before decompressing
    const restored = encoded.replace(/ /g, '+');

    // Step 2: Decompress
    const json = LZString.decompressFromEncodedURIComponent(restored);
    if (!json) {
      console.warn('[deeplink] Decompression failed: encoded string may be corrupted');
      return null;
    }

    // Step 3: Parse JSON
    let data: unknown;
    try {
      data = JSON.parse(json);
    } catch (parseError) {
      console.warn('[deeplink] JSON parse failed:', parseError);
      return null;
    }

    // Step 4: Validate with Zod
    const result = RunStateV1Schema.safeParse(data);
    if (!result.success) {
      console.warn('[deeplink] Validation failed:', result.error.format());
      return null;
    }

    return result.data;
  } catch (error) {
    // Catch any unexpected errors
    console.warn('[deeplink] Unexpected decode error:', error);
    return null;
  }
}

/**
 * Calculate the encoded length of a state object.
 *
 * Useful for checking if a state will fit within URL limits before
 * actually updating the URL.
 *
 * @param state - The state object to measure
 * @returns Encoded length in characters
 *
 * @example
 * ```typescript
 * import { estimateEncodedLength } from '@/lib/deeplink';
 *
 * const length = estimateEncodedLength(state);
 * if (length > 500) {
 *   console.warn('State may be too large for some browsers');
 * }
 * ```
 */
export function estimateEncodedLength(state: RunStateV1): number {
  return encodeRunState(state).length;
}

/**
 * Check if a state can be safely encoded to URL.
 *
 * Browser URL limits vary:
 * - Most modern browsers: ~2000 characters
 * - IE11 (legacy): ~2083 characters
 * - Server limits vary widely
 *
 * This function uses a conservative limit of 500 characters for the
 * encoded state parameter, leaving room for the base URL, path,
 * and other query parameters.
 *
 * @param state - The state object to check
 * @returns true if state encodes to less than 500 characters
 *
 * @example
 * ```typescript
 * import { isStateUrlSafe } from '@/lib/deeplink';
 *
 * if (!isStateUrlSafe(state)) {
 *   showWarning('Link may be too long for some browsers');
 * }
 * ```
 */
export function isStateUrlSafe(state: RunStateV1): boolean {
  return estimateEncodedLength(state) < 500;
}

/**
 * Verify that encode/decode round-trip preserves state.
 *
 * This is primarily useful for testing, but can also be used
 * to verify state integrity before sharing.
 *
 * @param state - The state to verify
 * @returns true if encode(decode(state)) equals state
 *
 * @example
 * ```typescript
 * import { verifyRoundTrip } from '@/lib/deeplink';
 *
 * if (!verifyRoundTrip(state)) {
 *   console.error('State fails round-trip verification');
 * }
 * ```
 */
export function verifyRoundTrip(state: RunStateV1): boolean {
  const encoded = encodeRunState(state);
  const decoded = decodeRunState(encoded);

  if (!decoded) {
    return false;
  }

  // Compare JSON representations for deep equality
  return JSON.stringify(state) === JSON.stringify(decoded);
}
