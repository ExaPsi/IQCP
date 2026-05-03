/**
 * useIntegralBreakdown - React hook for computing primitive-pair
 * breakdown of a single integral via WASM.
 *
 * Sends an `integral_breakdown` request to the Web Worker when a
 * cell is selected in the integral heatmap. Results are cached by
 * (geometry hash, basisName, integralType, i, j) to avoid redundant
 * WASM calls when revisiting the same cell.
 *
 * Cache is invalidated when the molecule or basis set changes.
 *
 * @module hooks/useIntegralBreakdown
 * @see US-056 Primitive Breakdown Panel
 * @see FR-INT-03
 */

import { useState, useEffect, useRef, useCallback } from 'react';
import { useWorker } from './useWorker';
import type { IntegralBreakdownResult } from '../components/integrals/PrimitiveBreakdownTable';
import type { GeometryInput, BasisSetName } from '../worker/protocol';

// Re-export the result type for convenience
export type { IntegralBreakdownResult } from '../components/integrals/PrimitiveBreakdownTable';

// ============================================================================
// Types
// ============================================================================

/**
 * Parameters for the useIntegralBreakdown hook.
 */
export interface UseIntegralBreakdownParams {
  /** Atom specifications for the molecule */
  geometry: GeometryInput | null;
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** Integral type matching the active matrix tab */
  integralType: 'S' | 'T' | 'V' | 'Hcore';
  /** Basis function row index (0-based) */
  indexI: number;
  /** Basis function column index (0-based) */
  indexJ: number;
  /** Whether the worker is ready */
  isReady: boolean;
  /** Whether the hook should compute (false = skip, e.g., no cell selected) */
  enabled: boolean;
}

/**
 * Return type of useIntegralBreakdown hook.
 */
export interface UseIntegralBreakdownReturn {
  /** Whether a computation is in progress */
  loading: boolean;
  /** Error message if the computation failed */
  error: string | null;
  /** Breakdown data, or null if not yet computed */
  data: IntegralBreakdownResult | null;
}

// ============================================================================
// Cache
// ============================================================================

/**
 * In-memory cache for integral breakdown results.
 *
 * Key format: `${geometryHash}-${basisName}-${integralType}-${i}-${j}`
 *
 * Using a module-level Map ensures cache persists across re-renders
 * and component remounts, but is cleared on page reload.
 */
const cache = new Map<string, IntegralBreakdownResult>();

/**
 * Maximum number of cached breakdown results.
 * Each result is small (a few KB for 9-36 primitives), so 100 entries
 * is well within memory budget.
 */
const MAX_CACHE_SIZE = 100;

/**
 * Generate a simple hash string from a geometry for cache keying.
 *
 * Uses atom symbols and coordinates rounded to 6 decimal places.
 * This is not cryptographic -- it only needs to distinguish different
 * molecular geometries for cache invalidation.
 */
function geometryHash(geometry: GeometryInput): string {
  const atomStr = geometry.atoms
    .map(
      (a) =>
        `${a.symbol}:${a.xyz[0].toFixed(6)},${a.xyz[1].toFixed(6)},${a.xyz[2].toFixed(6)}`,
    )
    .join('|');
  return `${atomStr}:${geometry.units}`;
}

/**
 * Build the cache key from all parameters that affect the result.
 */
function buildCacheKey(
  geometry: GeometryInput,
  basisName: string,
  integralType: string,
  i: number,
  j: number,
): string {
  return `${geometryHash(geometry)}-${basisName}-${integralType}-${i}-${j}`;
}

/**
 * Evict oldest entries if cache exceeds maximum size.
 */
function evictIfNeeded(): void {
  if (cache.size <= MAX_CACHE_SIZE) return;

  // Remove the first (oldest) entry
  const firstKey = cache.keys().next().value;
  if (firstKey !== undefined) {
    cache.delete(firstKey);
  }
}

// ============================================================================
// Hook
// ============================================================================

/**
 * React hook for computing the primitive-pair breakdown of a single integral.
 *
 * When enabled and the cell parameters change, sends an `integral_breakdown`
 * worker request. Results are cached to avoid re-computing the same cell.
 *
 * @param params - Hook parameters
 * @returns Object with loading state, error, and breakdown data
 *
 * @example
 * ```typescript
 * function DetailPanel({ row, col }: { row: number; col: number }) {
 *   const { isReady } = useWorker();
 *   const breakdown = useIntegralBreakdown({
 *     geometry: currentGeometry,
 *     basisName: 'sto-3g',
 *     integralType: 'S',
 *     indexI: row,
 *     indexJ: col,
 *     isReady,
 *     enabled: true,
 *   });
 *
 *   if (breakdown.loading) return <p>Loading breakdown...</p>;
 *   if (breakdown.error) return <p>Error: {breakdown.error}</p>;
 *   if (breakdown.data) return <PrimitiveBreakdownTable breakdown={breakdown.data} />;
 * }
 * ```
 */
export function useIntegralBreakdown({
  geometry,
  basisName,
  integralType,
  indexI,
  indexJ,
  isReady,
  enabled,
}: UseIntegralBreakdownParams): UseIntegralBreakdownReturn {
  const { send } = useWorker();
  const [data, setData] = useState<IntegralBreakdownResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Track the current request to avoid stale updates
  const requestRef = useRef(0);

  const fetchBreakdown = useCallback(
    async (
      geom: GeometryInput,
      basis: string,
      type: 'S' | 'T' | 'V' | 'Hcore',
      i: number,
      j: number,
    ) => {
      // Check cache first
      const key = buildCacheKey(geom, basis, type, i, j);
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
        // Send the integral_breakdown request to the worker.
        // The worker handler and protocol types are being added by a
        // concurrent agent. Use a type assertion until the protocol
        // union is updated.
        const request = {
          type: 'integral_breakdown' as const,
          geometry: geom,
          basisName: basis as BasisSetName,
          integralType: type,
          indices: [i, j] as [number, number],
        };
        const result = await send<IntegralBreakdownResult>(request as never);

        // Only update if this is still the latest request
        if (requestRef.current === requestId) {
          evictIfNeeded();
          cache.set(key, result);
          setData(result);
          setLoading(false);
        }
      } catch (err) {
        if (requestRef.current === requestId) {
          const message =
            err instanceof Error
              ? err.message
              : 'Failed to compute primitive breakdown';
          setError(message);
          setData(null);
          setLoading(false);
        }
      }
    },
    [send],
  );

  // Clear data when parameters change significantly (new molecule or basis)
  // but keep showing stale data during loading of a new cell in the same system
  const prevGeomRef = useRef<string | null>(null);
  const prevBasisRef = useRef<string | null>(null);

  useEffect(() => {
    if (!enabled || !isReady || !geometry) {
      // If disabled, clear everything
      if (!enabled) {
        setData(null);
        setLoading(false);
        setError(null);
      }
      return;
    }

    const currentGeomHash = geometryHash(geometry);

    // If the molecule or basis changed, clear stale data immediately
    if (
      prevGeomRef.current !== null &&
      (prevGeomRef.current !== currentGeomHash ||
        prevBasisRef.current !== basisName)
    ) {
      setData(null);
    }

    prevGeomRef.current = currentGeomHash;
    prevBasisRef.current = basisName;

    fetchBreakdown(geometry, basisName, integralType, indexI, indexJ);
  }, [
    enabled,
    isReady,
    geometry,
    basisName,
    integralType,
    indexI,
    indexJ,
    fetchBreakdown,
  ]);

  return { loading, error, data };
}
