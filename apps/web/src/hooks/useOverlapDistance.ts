/**
 * useOverlapDistance - React hook for computing overlap vs distance data via WASM.
 *
 * Sends an `overlap_distance` request to the Web Worker when the selected
 * shell pair changes, and returns the overlap curve data for plotting.
 * Results are cached to avoid redundant WASM calls when revisiting the same
 * shell pair configuration.
 *
 * Cache is invalidated when the element or basis set changes.
 *
 * @module hooks/useOverlapDistance
 * @see US-054 Overlap vs. Distance Plot
 */

import { useState, useEffect, useRef, useCallback } from 'react';
import { useWorker } from './useWorker';
import type { BasisSetName, OverlapDistanceResult } from '../worker/protocol';

// Re-export the result type for convenience
export type { OverlapDistanceResult } from '../worker/protocol';

/**
 * Parameters for the useOverlapDistance hook.
 */
export interface UseOverlapDistanceParams {
  /** Atomic number for atom A */
  elementA: number;
  /** Basis set name for atom A */
  basisNameA: string;
  /** Shell index within atom A's basis */
  shellIndexA: number;
  /** Atomic number for atom B */
  elementB: number;
  /** Basis set name for atom B */
  basisNameB: string;
  /** Shell index within atom B's basis */
  shellIndexB: number;
  /** Number of distance points (default: 100) */
  nPoints?: number;
  /** Maximum distance in bohr (default: auto-scaled from covalent radii) */
  rMax?: number;
  /** Whether the hook should compute (false = skip) */
  isReady: boolean;
}

/**
 * Return type of useOverlapDistance hook.
 */
export interface UseOverlapDistanceReturn {
  /** Whether a computation is in progress */
  loading: boolean;
  /** Error message if the computation failed */
  error: string | null;
  /** Overlap distance data, or null if not yet computed */
  data: OverlapDistanceResult | null;
}

// ============================================================================
// Covalent radii (bohr) for auto-scaling rMax
// ============================================================================

/**
 * Approximate covalent radii in bohr for H-Ar.
 * Used to auto-scale the maximum interatomic distance for overlap plots.
 */
const COVALENT_RADII_BOHR: Record<number, number> = {
  1: 0.59,   // H
  2: 0.56,   // He
  3: 2.42,   // Li
  4: 1.70,   // Be
  5: 1.55,   // B
  6: 1.44,   // C
  7: 1.32,   // N
  8: 1.21,   // O
  9: 1.12,   // F
  10: 1.12,  // Ne
  11: 3.02,  // Na
  12: 2.57,  // Mg
  13: 2.23,  // Al
  14: 2.10,  // Si
  15: 1.98,  // P
  16: 1.89,  // S
  17: 1.87,  // Cl
  18: 1.87,  // Ar
};

/**
 * Compute an appropriate rMax based on the covalent radii of two elements.
 *
 * Uses the heuristic: rMax = 3 * (r_cov_A + r_cov_B), clamped to [6.0, 15.0] bohr.
 * This ensures that the overlap curve is shown over a physically meaningful range:
 * smaller atoms (H-H) get a shorter range, while larger atoms (Na-Na) get a longer range.
 */
export function autoScaleRMax(elementA: number, elementB: number): number {
  const rA = COVALENT_RADII_BOHR[elementA] ?? 1.5;
  const rB = COVALENT_RADII_BOHR[elementB] ?? 1.5;
  const rMax = 3 * (rA + rB);
  return Math.max(6.0, Math.min(15.0, rMax));
}

// ============================================================================
// Cache
// ============================================================================

/**
 * In-memory cache for overlap distance results.
 * Key format: `${elA}-${basisA}-${shellA}-${elB}-${basisB}-${shellB}-${nPts}-${rMax}`
 */
const cache = new Map<string, OverlapDistanceResult>();

function cacheKey(
  elA: number,
  basisA: string,
  shellA: number,
  elB: number,
  basisB: string,
  shellB: number,
  nPoints: number,
  rMax: number,
): string {
  return `${elA}-${basisA}-${shellA}-${elB}-${basisB}-${shellB}-${nPoints}-${rMax}`;
}

// ============================================================================
// Hook
// ============================================================================

/**
 * React hook for computing overlap vs distance data from WASM.
 *
 * When the shell pair parameters change (and isReady is true), sends an
 * `overlap_distance` worker request. Results are cached to avoid re-computing
 * the same pair.
 *
 * @param params - Hook parameters
 * @returns Object with loading state, error, and overlap distance data
 *
 * @example
 * ```typescript
 * function OverlapPanel() {
 *   const { isReady } = useWorker();
 *   const overlap = useOverlapDistance({
 *     elementA: 1, basisNameA: 'sto-3g', shellIndexA: 0,
 *     elementB: 1, basisNameB: 'sto-3g', shellIndexB: 0,
 *     isReady,
 *   });
 *
 *   if (overlap.loading) return <p>Computing...</p>;
 *   if (overlap.error) return <p>Error: {overlap.error}</p>;
 *   if (overlap.data) return <OverlapDistancePlot traces={[overlap.data]} />;
 * }
 * ```
 */
export function useOverlapDistance({
  elementA,
  basisNameA,
  shellIndexA,
  elementB,
  basisNameB,
  shellIndexB,
  nPoints = 100,
  rMax: rMaxProp,
  isReady,
}: UseOverlapDistanceParams): UseOverlapDistanceReturn {
  // Auto-scale rMax from covalent radii if not explicitly provided
  const rMax = rMaxProp ?? autoScaleRMax(elementA, elementB);
  const { send } = useWorker();
  const [data, setData] = useState<OverlapDistanceResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Track the current request to avoid stale updates
  const requestRef = useRef(0);

  const fetchOverlap = useCallback(
    async (
      elA: number,
      bA: string,
      sA: number,
      elB: number,
      bB: string,
      sB: number,
      pts: number,
      rm: number,
    ) => {
      // Check cache first
      const key = cacheKey(elA, bA, sA, elB, bB, sB, pts, rm);
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
        const request = {
          type: 'overlap_distance' as const,
          elementA: elA,
          basisA: bA as BasisSetName,
          shellIndexA: sA,
          elementB: elB,
          basisB: bB as BasisSetName,
          shellIndexB: sB,
          rMin: 0.1,
          rMax: rm,
          nPoints: pts,
        };
        const result = await send<OverlapDistanceResult>(request);

        // Only update if this is still the latest request
        if (requestRef.current === requestId) {
          cache.set(key, result);
          setData(result);
          setLoading(false);
        }
      } catch (err) {
        if (requestRef.current === requestId) {
          const message =
            err instanceof Error ? err.message : 'Failed to compute overlap vs distance';
          setError(message);
          setData(null);
          setLoading(false);
        }
      }
    },
    [send],
  );

  useEffect(() => {
    if (!isReady) {
      return;
    }

    fetchOverlap(
      elementA,
      basisNameA,
      shellIndexA,
      elementB,
      basisNameB,
      shellIndexB,
      nPoints,
      rMax,
    );
  }, [
    isReady,
    elementA,
    basisNameA,
    shellIndexA,
    elementB,
    basisNameB,
    shellIndexB,
    nPoints,
    rMax,
    fetchOverlap,
  ]);

  return { loading, error, data };
}
