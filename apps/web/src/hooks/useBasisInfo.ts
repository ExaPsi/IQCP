/**
 * useBasisInfo - React hook for querying basis set shell data from WASM.
 *
 * Fetches shell information (angular momentum, primitives, exponents,
 * coefficients) for a given element and basis set combination. Data is
 * sourced from qc-core's builtin.rs via the WASM bridge, ensuring
 * a single source of truth (US-049 AC6).
 *
 * Results are cached to avoid redundant WASM calls when revisiting
 * the same element/basis combination.
 *
 * @module hooks/useBasisInfo
 */

import { useState, useEffect, useRef, useCallback } from 'react';
import { useWorker } from './useWorker';
import type { BasisInfoRequest, BasisInfoResult, BasisShellInfo } from '../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * Return type of useBasisInfo hook.
 */
export interface UseBasisInfoResult {
  /** Shell data from WASM, or null if not yet loaded */
  shells: BasisShellInfo[] | null;
  /** Whether data is currently being fetched */
  loading: boolean;
  /** Error message if the query failed */
  error: string | null;
}

// ============================================================================
// Cache
// ============================================================================

/**
 * In-memory cache for basis info results.
 * Key format: `${atomicNumber}-${basisName}`
 */
const cache = new Map<string, BasisShellInfo[]>();

function cacheKey(z: number, basis: string): string {
  return `${z}-${basis}`;
}

// ============================================================================
// Hook
// ============================================================================

/**
 * React hook for querying basis set shell information from WASM.
 *
 * Sends a `basis_info` request to the Web Worker when the element
 * or basis set changes, and returns the shell data.
 *
 * @param atomicNumber - Atomic number of the element (1-18)
 * @param basisName - Basis set name (e.g., "sto-3g", "6-31g*")
 * @returns Object with shells data, loading state, and error
 *
 * @example
 * ```typescript
 * function ShellDisplay() {
 *   const { shells, loading, error } = useBasisInfo(8, 'sto-3g');
 *
 *   if (loading) return <p>Loading...</p>;
 *   if (error) return <p>Error: {error}</p>;
 *   if (!shells) return null;
 *
 *   return (
 *     <ul>
 *       {shells.map((s, i) => (
 *         <li key={i}>{s.angularMomentumLabel} shell: {s.nPrimitives} primitives</li>
 *       ))}
 *     </ul>
 *   );
 * }
 * ```
 */
export function useBasisInfo(atomicNumber: number, basisName: string): UseBasisInfoResult {
  const { send, isReady } = useWorker();
  const [shells, setShells] = useState<BasisShellInfo[] | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Track the current request to avoid stale updates
  const requestRef = useRef(0);

  const fetchBasisInfo = useCallback(
    async (z: number, basis: string) => {
      // Check cache first
      const key = cacheKey(z, basis);
      const cached = cache.get(key);
      if (cached) {
        setShells(cached);
        setLoading(false);
        setError(null);
        return;
      }

      const requestId = ++requestRef.current;
      setLoading(true);
      setError(null);

      try {
        const request: Omit<BasisInfoRequest, 'requestId'> = {
          type: 'basis_info',
          atomicNumber: z,
          basisName: basis,
        };
        const result = await send<BasisInfoResult>(request);

        // Only update if this is still the latest request
        if (requestRef.current === requestId) {
          cache.set(key, result);
          setShells(result);
          setLoading(false);
        }
      } catch (err) {
        if (requestRef.current === requestId) {
          const message = err instanceof Error ? err.message : 'Failed to fetch basis info';
          setError(message);
          setShells(null);
          setLoading(false);
        }
      }
    },
    [send]
  );

  useEffect(() => {
    if (!isReady) return;
    fetchBasisInfo(atomicNumber, basisName);
  }, [isReady, atomicNumber, basisName, fetchBasisInfo]);

  return { shells, loading, error };
}
