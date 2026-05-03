/**
 * useRadialProfile - React hook for computing radial profile data via WASM.
 *
 * Sends a `radial_profile` request to the Web Worker when the selected
 * shell changes, and returns the profile data for plotting. Results are
 * cached to avoid redundant WASM calls when revisiting the same shell.
 *
 * Cache is invalidated when the element or basis set changes.
 *
 * @module hooks/useRadialProfile
 */

import { useState, useEffect, useRef, useCallback } from 'react';
import { useWorker } from './useWorker';
import type {
  RadialProfileRequest,
  RadialProfileResult,
} from '../worker/protocol';

// Re-export the result type for convenience
export type { RadialProfileResult } from '../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * Parameters for the useRadialProfile hook.
 */
export interface UseRadialProfileParams {
  /** Atomic number of the element (1-18) */
  atomicNumber: number;
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** Shell index (0-indexed), or null if no shell is selected */
  shellIndex: number | null;
  /** Whether the worker is ready */
  isReady: boolean;
}

/**
 * Return type of useRadialProfile hook.
 */
export interface UseRadialProfileReturn {
  /** Whether a computation is in progress */
  loading: boolean;
  /** Error message if the computation failed */
  error: string | null;
  /** Radial profile data, or null if not yet computed */
  data: RadialProfileResult | null;
}

// ============================================================================
// Cache
// ============================================================================

/**
 * In-memory cache for radial profile results.
 * Key format: `${atomicNumber}-${basisName}-${shellIndex}`
 */
const cache = new Map<string, RadialProfileResult>();

function cacheKey(z: number, basis: string, shellIdx: number): string {
  return `${z}-${basis}-${shellIdx}`;
}

/**
 * Clear all cached entries for a given element and basis combination.
 * Called when the element or basis changes to avoid stale data.
 */
function invalidateCache(z: number, basis: string): void {
  const prefix = `${z}-${basis}-`;
  for (const key of cache.keys()) {
    if (key.startsWith(prefix)) {
      cache.delete(key);
    }
  }
}

// ============================================================================
// Hook
// ============================================================================

/**
 * React hook for computing radial profile data from WASM.
 *
 * When the shellIndex changes (and is not null), sends a `radial_profile`
 * worker request. Results are cached to avoid re-computing the same shell.
 * Cache is invalidated when the element or basis set changes.
 *
 * @param params - Hook parameters
 * @returns Object with loading state, error, and profile data
 *
 * @example
 * ```typescript
 * function RadialProfilePanel() {
 *   const { isReady } = useWorker();
 *   const profile = useRadialProfile({
 *     atomicNumber: 8,
 *     basisName: 'sto-3g',
 *     shellIndex: 0,
 *     isReady,
 *   });
 *
 *   if (profile.loading) return <p>Computing...</p>;
 *   if (profile.error) return <p>Error: {profile.error}</p>;
 *   if (profile.data) return <RadialProfilePlot profileData={profile.data} />;
 * }
 * ```
 */
export function useRadialProfile({
  atomicNumber,
  basisName,
  shellIndex,
  isReady,
}: UseRadialProfileParams): UseRadialProfileReturn {
  const { send } = useWorker();
  const [data, setData] = useState<RadialProfileResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Track the current request to avoid stale updates
  const requestRef = useRef(0);

  // Track previous element/basis to detect changes for cache invalidation
  const prevParamsRef = useRef({ atomicNumber, basisName });

  // Invalidate cache when element or basis changes
  useEffect(() => {
    const prev = prevParamsRef.current;
    if (prev.atomicNumber !== atomicNumber || prev.basisName !== basisName) {
      invalidateCache(prev.atomicNumber, prev.basisName);
      prevParamsRef.current = { atomicNumber, basisName };
      // Clear current data when switching element/basis
      setData(null);
      setError(null);
    }
  }, [atomicNumber, basisName]);

  const fetchProfile = useCallback(
    async (z: number, basis: string, shellIdx: number) => {
      // Check cache first
      const key = cacheKey(z, basis, shellIdx);
      const cached = cache.get(key);
      if (cached) {
        setData(cached);
        setLoading(false);
        setError(null);
        return;
      }

      const requestId = ++requestRef.current;
      setLoading(true);
      setError(null);

      try {
        const request: Omit<RadialProfileRequest, 'requestId'> = {
          type: 'radial_profile',
          atomicNumber: z,
          basisName: basis,
          shellIndex: shellIdx,
        };
        const result = await send<RadialProfileResult>(request);

        // Only update if this is still the latest request
        if (requestRef.current === requestId) {
          cache.set(key, result);
          setData(result);
          setLoading(false);
        }
      } catch (err) {
        if (requestRef.current === requestId) {
          const message =
            err instanceof Error ? err.message : 'Failed to compute radial profile';
          setError(message);
          setData(null);
          setLoading(false);
        }
      }
    },
    [send]
  );

  useEffect(() => {
    if (!isReady || shellIndex === null) {
      // Clear data when no shell is selected
      if (shellIndex === null) {
        setData(null);
        setError(null);
        setLoading(false);
      }
      return;
    }

    fetchProfile(atomicNumber, basisName, shellIndex);
  }, [isReady, atomicNumber, basisName, shellIndex, fetchProfile]);

  return { loading, error, data };
}
