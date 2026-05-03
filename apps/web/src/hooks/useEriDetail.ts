/**
 * useEriDetail - React hook for computing ERI primitive-quartet breakdown
 * via the WASM Web Worker.
 *
 * Sends an `eri_detail` request to the worker when the user selects a
 * quartet and clicks "Inspect". Unlike useIntegralBreakdown, this hook
 * does NOT auto-compute -- it exposes a `compute()` function that the
 * caller triggers explicitly.
 *
 * Results are cached by (geometry hash, basisName, i, j, k, l) to avoid
 * redundant WASM calls when revisiting the same quartet.
 *
 * @module hooks/useEriDetail
 * @see US-059 ERI Browser
 * @see FR-INT-04
 */

import { useState, useRef, useCallback, useEffect } from 'react';
import { useWorker } from './useWorker';
import { useIntegralStore } from '../stores/integralStore';
import type { GeometryInput, BasisSetName } from '../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * A single primitive-quartet contribution to a contracted ERI.
 *
 * Defined locally to decouple from protocol.ts, which may be
 * updated concurrently by another agent.
 */
export interface EriPrimitiveContribution {
  /** Primitive indices [p, q, r, s] within their shells (0-based) */
  primIndices: [number, number, number, number];
  /** Primitive exponents [alpha_p, alpha_q, alpha_r, alpha_s] */
  exponents: [number, number, number, number];
  /** Raw contraction coefficients (without normalization) */
  coefficients: [number, number, number, number];
  /** Contraction coefficients with normalization */
  normCoefficients: [number, number, number, number];
  /** Bare primitive ERI value */
  primitiveValue: number;
  /** Weighted contribution = product(normCoefficients) * primitiveValue */
  weightedContribution: number;
  /** T parameter for this primitive quartet */
  tParameter: number;
}

/**
 * Method used to compute the ERI.
 *
 * Discriminated union matching Rust's EriMethod enum.
 */
export type EriMethod =
  | { type: 'boysFunction'; tParameter: number }
  | { type: 'rysQuadrature'; nroots: number; tParameter: number; roots: number[]; weights: number[] };

/**
 * Result from eri_detail worker request.
 *
 * Contains the contracted ERI value and all primitive-quartet
 * contributions sorted by |weightedContribution| descending.
 */
export interface EriDetailResult {
  /** Contracted ERI value */
  contractedValue: number;
  /** Basis function indices [i, j, k, l] */
  indices: [number, number, number, number];
  /** Basis function labels */
  labels: [string, string, string, string];
  /** Computation method (Boys function or Rys quadrature) */
  method: EriMethod;
  /** Primitive-quartet contributions sorted by |weightedContribution| descending */
  contributions: EriPrimitiveContribution[];
  /** Number of primitives per shell [n_i, n_j, n_k, n_l] */
  nPrimitives: [number, number, number, number];
  /** Total angular momentum L_total */
  totalAngularMomentum: number;
  /** Number of Rys quadrature roots */
  nroots: number;
}

/**
 * Parameters for the useEriDetail hook.
 */
export interface UseEriDetailParams {
  /** Atom specifications for the molecule */
  geometry: GeometryInput | null;
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** Selected quartet [i, j, k, l] (null if none) */
  quartet: [number, number, number, number] | null;
  /** Whether the worker is ready */
  isReady: boolean;
}

/**
 * Return type of useEriDetail hook.
 */
export interface UseEriDetailReturn {
  /** ERI breakdown data, or null if not yet computed */
  data: EriDetailResult | null;
  /** Whether a computation is in progress */
  loading: boolean;
  /** Error message if the computation failed */
  error: string | null;
  /** Trigger an ERI detail computation for the current quartet */
  compute: () => void;
}

// ============================================================================
// Cache
// ============================================================================

/**
 * In-memory cache for ERI detail results.
 *
 * Key format: `${geometryHash}-${basisName}-${i}-${j}-${k}-${l}`
 *
 * Module-level Map persists across re-renders and component remounts,
 * cleared on page reload.
 */
const cache = new Map<string, EriDetailResult>();

/** Maximum cached entries (each result is a few KB for STO-3G). */
const MAX_CACHE_SIZE = 50;

/**
 * Generate a simple hash string from a geometry for cache keying.
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
  quartet: [number, number, number, number],
): string {
  return `${geometryHash(geometry)}-${basisName}-${quartet[0]}-${quartet[1]}-${quartet[2]}-${quartet[3]}`;
}

/**
 * Evict oldest entries if cache exceeds maximum size.
 */
function evictIfNeeded(): void {
  if (cache.size <= MAX_CACHE_SIZE) return;
  const firstKey = cache.keys().next().value;
  if (firstKey !== undefined) {
    cache.delete(firstKey);
  }
}

// ============================================================================
// Hook
// ============================================================================

/**
 * React hook for computing the primitive-quartet breakdown of a single ERI.
 *
 * Unlike useIntegralBreakdown (which auto-computes on cell selection),
 * this hook exposes a `compute()` function that the caller triggers
 * explicitly, typically when the user clicks "Inspect".
 *
 * Results are cached to avoid re-computing the same quartet.
 *
 * @param params - Hook parameters
 * @returns Object with loading state, error, data, and compute trigger
 *
 * @example
 * ```typescript
 * function EriBrowserPanel() {
 *   const { isReady } = useWorker();
 *   const eri = useEriDetail({
 *     geometry: currentGeometry,
 *     basisName: 'sto-3g',
 *     quartet: [0, 0, 1, 1],
 *     isReady,
 *   });
 *
 *   return (
 *     <button onClick={eri.compute} disabled={eri.loading}>
 *       Inspect
 *     </button>
 *   );
 * }
 * ```
 */
export function useEriDetail({
  geometry,
  basisName,
  quartet,
  isReady,
}: UseEriDetailParams): UseEriDetailReturn {
  const { send } = useWorker();
  const [data, setData] = useState<EriDetailResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Store actions for syncing ERI detail state
  const setEriDetail = useIntegralStore((s) => s.setEriDetail);

  // Track the current request to avoid stale updates
  const requestRef = useRef(0);

  // Reset when quartet or geometry changes
  const prevKeyRef = useRef<string | null>(null);

  useEffect(() => {
    if (!geometry || !quartet) {
      setData(null);
      setLoading(false);
      setError(null);
      setEriDetail({ status: 'idle' });
      prevKeyRef.current = null;
      return;
    }

    const key = buildCacheKey(geometry, basisName, quartet);
    if (prevKeyRef.current !== null && prevKeyRef.current !== key) {
      // Quartet or system changed -- reset results but don't auto-compute
      setData(null);
      setError(null);
      setEriDetail({ status: 'idle' });
    }
    prevKeyRef.current = key;
  }, [geometry, basisName, quartet, setEriDetail]);

  const compute = useCallback(() => {
    if (!isReady || !geometry || !quartet) return;

    // Check cache first
    const key = buildCacheKey(geometry, basisName, quartet);
    const cached = cache.get(key);
    if (cached) {
      setData(cached);
      setLoading(false);
      setError(null);
      setEriDetail({ status: 'success', data: cached });
      return;
    }

    const requestId = ++requestRef.current;
    setLoading(true);
    setError(null);
    setEriDetail({ status: 'computing' });

    const request = {
      type: 'eri_detail' as const,
      geometry,
      basisName: basisName as BasisSetName,
      indices: quartet,
    };

    send<EriDetailResult>(request as never)
      .then((result) => {
        if (requestRef.current === requestId) {
          evictIfNeeded();
          cache.set(key, result);
          setData(result);
          setLoading(false);
          setEriDetail({ status: 'success', data: result });
        }
      })
      .catch((err) => {
        if (requestRef.current === requestId) {
          const message =
            err instanceof Error
              ? err.message
              : 'Failed to compute ERI breakdown';
          setError(message);
          setData(null);
          setLoading(false);
          setEriDetail({ status: 'error', error: message });
        }
      });
  }, [isReady, geometry, basisName, quartet, send, setEriDetail]);

  return { data, loading, error, compute };
}
