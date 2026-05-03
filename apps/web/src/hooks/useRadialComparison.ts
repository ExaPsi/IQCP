/**
 * useRadialComparison - React hook for comparing radial profiles across basis sets.
 *
 * Enables side-by-side comparison of the same shell type (e.g., "1s") across
 * different basis sets (e.g., STO-3G, 3-21G, 6-31G) for Module D's Basis
 * Explorer. For each selected basis set, the hook:
 *
 * 1. Fetches shell structure via the WASM worker
 * 2. Matches shells by the target shell type label
 * 3. Fetches radial profiles for each matching shell
 * 4. Returns grouped results with role annotations
 *
 * Split-valence bases return multiple profiles per basis set (inner + outer
 * valence), which is the key pedagogical insight this hook enables.
 *
 * @module hooks/useRadialComparison
 * @see US-052 Basis Set Comparison Mode (Module D)
 */

import { useState, useEffect, useRef, useCallback } from 'react';
import { useWorker } from './useWorker';
import { deriveShellDisplayInfo } from '../components/basisExplorer/shellDisplayInfo';
import { findMatchingShells } from '../components/basisExplorer/shellMatching';
import type {
  BasisInfoRequest,
  BasisShellInfo,
  RadialProfileRequest,
  RadialProfileResult,
} from '../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * A single radial profile associated with a basis set and shell role.
 */
export interface ComparisonProfile {
  /** Basis set internal name (e.g., "sto-3g") */
  basisName: string;
  /** Basis set display label (e.g., "STO-3G") */
  basisLabel: string;
  /** Shell role within the basis (e.g., "valence", "inner valence", "outer valence") */
  shellRole: string;
  /** The radial profile data from WASM */
  profile: RadialProfileResult;
}

/**
 * Parameters for the useRadialComparison hook.
 */
export interface UseRadialComparisonParams {
  /** Atomic number (1-18) */
  atomicNumber: number;
  /** Shell type label to compare (e.g., "1s", "2p") */
  shellTypeLabel: string;
  /** List of basis set names to compare */
  selectedBases: string[];
  /** Whether comparison mode is active */
  isActive: boolean;
}

/**
 * Return value of the useRadialComparison hook.
 */
export interface UseRadialComparisonReturn {
  /** Map of basis name -> ComparisonProfile[] */
  profiles: Map<string, ComparisonProfile[]>;
  /** Whether any fetch is in progress */
  loading: boolean;
  /** Error message if any fetch failed */
  error: string | null;
}

// ============================================================================
// Display Label Mapping
// ============================================================================

const BASIS_DISPLAY_LABELS: Record<string, string> = {
  'sto-3g': 'STO-3G',
  '3-21g': '3-21G',
  '6-31g': '6-31G',
  '6-31g*': '6-31G*',
  '6-31+g*': '6-31+G*',
};

// ============================================================================
// Cache
// ============================================================================

/**
 * Cache for basis info results.
 * Key: `${atomicNumber}-${basisName}`
 */
const basisInfoCache = new Map<string, BasisShellInfo[]>();

/**
 * Cache for radial profile results.
 * Key: `${atomicNumber}-${basisName}-${shellIndex}`
 */
const profileCache = new Map<string, RadialProfileResult>();

// ============================================================================
// Hook
// ============================================================================

/**
 * React hook for comparing radial profiles across multiple basis sets.
 *
 * @param params - Configuration for the comparison
 * @returns Object with profiles map, loading state, and error
 *
 * @example
 * ```typescript
 * const { profiles, loading, error } = useRadialComparison({
 *   atomicNumber: 1,
 *   shellTypeLabel: '1s',
 *   selectedBases: ['sto-3g', '3-21g', '6-31g'],
 *   isActive: true,
 * });
 *
 * // profiles.get('3-21g') returns:
 * // [
 * //   { basisName: '3-21g', basisLabel: '3-21G', shellRole: 'inner valence', profile: {...} },
 * //   { basisName: '3-21g', basisLabel: '3-21G', shellRole: 'outer valence', profile: {...} },
 * // ]
 * ```
 */
export function useRadialComparison({
  atomicNumber,
  shellTypeLabel,
  selectedBases,
  isActive,
}: UseRadialComparisonParams): UseRadialComparisonReturn {
  const { send, isReady } = useWorker();
  const [profiles, setProfiles] = useState<Map<string, ComparisonProfile[]>>(new Map());
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Track request generation to avoid stale updates
  const requestGenRef = useRef(0);

  /**
   * Fetch basis info for a single element+basis, using cache.
   */
  const fetchBasisInfo = useCallback(
    async (z: number, basis: string): Promise<BasisShellInfo[]> => {
      const key = `${z}-${basis}`;
      const cached = basisInfoCache.get(key);
      if (cached) return cached;

      const request: Omit<BasisInfoRequest, 'requestId'> = {
        type: 'basis_info',
        atomicNumber: z,
        basisName: basis,
      };
      const result = await send<BasisShellInfo[]>(request);
      basisInfoCache.set(key, result);
      return result;
    },
    [send],
  );

  /**
   * Fetch a single radial profile, using cache.
   */
  const fetchProfile = useCallback(
    async (z: number, basis: string, shellIndex: number): Promise<RadialProfileResult> => {
      const key = `${z}-${basis}-${shellIndex}`;
      const cached = profileCache.get(key);
      if (cached) return cached;

      const request: Omit<RadialProfileRequest, 'requestId'> = {
        type: 'radial_profile',
        atomicNumber: z,
        basisName: basis,
        shellIndex,
      };
      const result = await send<RadialProfileResult>(request);
      profileCache.set(key, result);
      return result;
    },
    [send],
  );

  // Stable serialization of selectedBases for dependency tracking
  const basesKey = selectedBases.join(',');

  useEffect(() => {
    // Guard: don't fetch if inactive, not ready, or no bases selected
    if (!isActive || !isReady || selectedBases.length === 0 || !shellTypeLabel) {
      if (!isActive) {
        setProfiles(new Map());
        setError(null);
        setLoading(false);
      }
      return;
    }

    const gen = ++requestGenRef.current;

    async function fetchAll() {
      setLoading(true);
      setError(null);

      try {
        const resultMap = new Map<string, ComparisonProfile[]>();

        // Fetch all basis sets in parallel
        const basisPromises = selectedBases.map(async (basis) => {
          try {
            // Step 1: Get shell structure for this element+basis
            const wasmShells = await fetchBasisInfo(atomicNumber, basis);
            const displayInfo = deriveShellDisplayInfo(wasmShells, atomicNumber);

            // Step 2: Find matching shells for the target shell type
            const matches = findMatchingShells(wasmShells, displayInfo, shellTypeLabel);

            if (matches.length === 0) {
              // No matching shells (e.g., asking for "3d" from STO-3G)
              return;
            }

            // Step 3: Fetch radial profiles for each matching shell
            const profilePromises = matches.map(async (match) => {
              const profile = await fetchProfile(atomicNumber, basis, match.shellIndex);
              return {
                basisName: basis,
                basisLabel: BASIS_DISPLAY_LABELS[basis] ?? basis.toUpperCase(),
                shellRole: match.role,
                profile,
              } satisfies ComparisonProfile;
            });

            const basisProfiles = await Promise.all(profilePromises);
            resultMap.set(basis, basisProfiles);
          } catch (err) {
            // Individual basis failure: log and continue with others
            console.warn(`[useRadialComparison] Failed for ${basis}:`, err);
          }
        });

        await Promise.all(basisPromises);

        // Only update if this is still the latest request
        if (requestGenRef.current === gen) {
          setProfiles(resultMap);
          setLoading(false);

          // Partial failure: some bases failed
          if (resultMap.size < selectedBases.length && resultMap.size > 0) {
            setError(
              `Failed to load ${selectedBases.length - resultMap.size} of ${selectedBases.length} basis sets`,
            );
          }
        }
      } catch (err) {
        if (requestGenRef.current === gen) {
          const message = err instanceof Error ? err.message : 'Failed to fetch comparison data';
          setError(message);
          setProfiles(new Map());
          setLoading(false);
        }
      }
    }

    fetchAll();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [isActive, isReady, atomicNumber, shellTypeLabel, basesKey, fetchBasisInfo, fetchProfile]);

  return { profiles, loading, error };
}
