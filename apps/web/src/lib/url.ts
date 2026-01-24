/**
 * URL State Management
 *
 * Provides utilities for reading and writing application state to/from
 * the browser URL. Uses the History API to update URLs without page
 * reloads, preserving a smooth user experience.
 *
 * @module url
 */

import { encodeRunState, decodeRunState } from './deeplink';
import type { RunStateV1 } from '../types/run-state';

/**
 * URL parameter name for encoded state.
 *
 * The parameter name is intentionally short to minimize URL length.
 * Example URL: /boys?run=N4IgdghgtgpiBc...
 */
export const STATE_PARAM = 'run';

/**
 * Update the URL with encoded state without page reload.
 *
 * Uses history.replaceState to avoid polluting browser history
 * with every state change. This means the back button won't
 * step through every parameter adjustment.
 *
 * @param state - The state to encode into URL
 *
 * @example
 * ```typescript
 * import { updateURL } from '@/lib/url';
 *
 * // User changes T value
 * const newState = { ...state, boys: { ...state.boys, T: 5.0 } };
 * updateURL(newState);
 * // URL becomes: /boys?run=<encoded>
 * ```
 */
export function updateURL(state: RunStateV1): void {
  const encoded = encodeRunState(state);
  const url = new URL(window.location.href);
  url.searchParams.set(STATE_PARAM, encoded);
  window.history.replaceState(null, '', url.toString());
}

/**
 * Read and decode state from the current URL.
 *
 * Reads the 'run' query parameter and decodes it to a RunStateV1 object.
 * Returns null if no state parameter exists or if the state is invalid.
 *
 * @returns Decoded state, or null if no valid state in URL
 *
 * @example
 * ```typescript
 * import { getStateFromURL } from '@/lib/url';
 * import { DEFAULT_BOYS_STATE } from '@/types/run-state';
 *
 * // On page mount
 * const urlState = getStateFromURL();
 * if (urlState) {
 *   store.setState(urlState);
 * } else {
 *   store.setState(DEFAULT_BOYS_STATE);
 * }
 * ```
 */
export function getStateFromURL(): RunStateV1 | null {
  const params = new URLSearchParams(window.location.search);
  const encoded = params.get(STATE_PARAM);
  return decodeRunState(encoded);
}

/**
 * Remove the state parameter from URL.
 *
 * Useful after showing an error or when user clicks "reset".
 * Uses replaceState to avoid adding history entries.
 *
 * @example
 * ```typescript
 * import { clearURL } from '@/lib/url';
 * import { DEFAULT_BOYS_STATE } from '@/types/run-state';
 *
 * function handleReset() {
 *   clearURL();
 *   store.setState(DEFAULT_BOYS_STATE);
 * }
 * ```
 */
export function clearURL(): void {
  const url = new URL(window.location.href);
  url.searchParams.delete(STATE_PARAM);
  window.history.replaceState(null, '', url.toString());
}

/**
 * Generate a shareable URL for a state.
 *
 * Creates a complete URL that can be copied and shared with others.
 * The URL includes the encoded state as a query parameter.
 *
 * @param state - The state to encode
 * @param baseUrl - Optional base URL (defaults to current origin + module path)
 * @returns Full URL string with encoded state
 *
 * @example
 * ```typescript
 * import { generateShareURL } from '@/lib/url';
 *
 * const shareUrl = generateShareURL(state);
 * await navigator.clipboard.writeText(shareUrl);
 * showToast('Link copied to clipboard!');
 * ```
 *
 * @example
 * ```typescript
 * // With custom base URL for production
 * const shareUrl = generateShareURL(state, 'https://iqcp.example.edu/boys');
 * ```
 */
export function generateShareURL(state: RunStateV1, baseUrl?: string): string {
  const base = baseUrl ?? `${window.location.origin}/${state.module}`;
  const encoded = encodeRunState(state);
  const url = new URL(base);
  url.searchParams.set(STATE_PARAM, encoded);
  return url.toString();
}

/**
 * Check if the current URL has a state parameter.
 *
 * Useful for detecting whether the user arrived via a deep link
 * vs navigating directly to the page.
 *
 * @returns true if URL contains a state parameter (whether valid or not)
 *
 * @example
 * ```typescript
 * import { hasStateInURL, getStateFromURL } from '@/lib/url';
 *
 * if (hasStateInURL()) {
 *   const state = getStateFromURL();
 *   if (state) {
 *     applyState(state);
 *   } else {
 *     showInvalidLinkError();
 *   }
 * } else {
 *   useDefaultState();
 * }
 * ```
 */
export function hasStateInURL(): boolean {
  const params = new URLSearchParams(window.location.search);
  return params.has(STATE_PARAM);
}

/**
 * Push a new history entry with state.
 *
 * Unlike updateURL which uses replaceState, this creates a new
 * history entry, allowing the user to navigate back to previous
 * states using the browser back button.
 *
 * Use sparingly - typically only for significant state changes
 * like switching modules or loading a preset.
 *
 * @param state - The state to encode into URL
 *
 * @example
 * ```typescript
 * import { pushStateToURL } from '@/lib/url';
 *
 * // User loads a preset - create history entry
 * function loadPreset(presetState: RunStateV1) {
 *   pushStateToURL(presetState);
 *   applyState(presetState);
 * }
 * ```
 */
export function pushStateToURL(state: RunStateV1): void {
  const encoded = encodeRunState(state);
  const url = new URL(window.location.href);
  url.searchParams.set(STATE_PARAM, encoded);
  window.history.pushState(null, '', url.toString());
}

/**
 * Get the raw encoded state string from URL without decoding.
 *
 * Useful for debugging or when you need to pass the encoded
 * string to another function.
 *
 * @returns The encoded state string, or null if not present
 *
 * @example
 * ```typescript
 * import { getRawStateFromURL } from '@/lib/url';
 *
 * const encoded = getRawStateFromURL();
 * console.log('Encoded state:', encoded);
 * console.log('Length:', encoded?.length ?? 0);
 * ```
 */
export function getRawStateFromURL(): string | null {
  const params = new URLSearchParams(window.location.search);
  return params.get(STATE_PARAM);
}
