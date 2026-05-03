/**
 * useExponentSliders - Hook for managing exponent slider state and computation.
 *
 * Encapsulates the logic for:
 * 1. Reading default exponents from the WASM-computed RadialProfileResult
 * 2. Reading current overrides from basisStore.exponentOverrides[shellIndex]
 * 3. Computing the modified profile client-side via evaluateRadialProfileTS
 * 4. Computing the overlap metric via computeOverlapMetric
 * 5. Debouncing: slider onChange fires frequently, but profile recomputation
 *    is debounced at 50ms (well under the 100ms AC2 budget)
 *
 * @see US-053 Exponent Sliders
 * @module hooks/useExponentSliders
 */

import { useMemo, useCallback, useState, useEffect, useRef } from 'react';
import { useBasisStore } from '../stores/basisStore';
import {
  evaluateRadialProfileTS,
  computeOverlapMetric,
} from '../components/basisExplorer/radialMath';
import type { RadialProfileResult } from '../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * Parameters for the useExponentSliders hook.
 */
export interface UseExponentSlidersParams {
  /** Shell index (0-indexed), or null if no shell is selected */
  shellIndex: number | null;
  /** WASM-computed radial profile for the selected shell (the "original") */
  profileData: RadialProfileResult | null;
}

/**
 * Return type of useExponentSliders hook.
 */
export interface UseExponentSlidersReturn {
  /** Default (original) exponents from the basis set */
  defaultExponents: number[];
  /** Current (possibly modified) exponents */
  currentExponents: number[];
  /** Modified profile R(r) values (null if not modified) */
  modifiedProfile: number[] | null;
  /** Normalized overlap metric S (1.0 when unmodified) */
  overlapMetric: number;
  /** Whether any exponent has been modified */
  isModified: boolean;
  /** Callback to change a single exponent */
  handleExponentChange: (primitiveIndex: number, newAlpha: number) => void;
  /** Callback to reset all exponents for the current shell */
  handleReset: () => void;
}

/**
 * Debounce delay for profile recomputation (ms).
 * 50ms is well under the 100ms AC2 budget since the actual computation
 * takes < 0.1ms for typical shells (3 primitives x 200 r-points).
 */
const RECOMPUTE_DEBOUNCE_MS = 50;

// ============================================================================
// Hook
// ============================================================================

/**
 * React hook for managing exponent slider state and real-time profile recomputation.
 *
 * Connects the basisStore exponent overrides, the WASM-computed original profile,
 * and the client-side TypeScript recomputation for instant slider response.
 *
 * @param params - Hook parameters
 * @returns Object with exponents, modified profile, overlap metric, and handlers
 *
 * @example
 * ```typescript
 * const {
 *   defaultExponents,
 *   currentExponents,
 *   modifiedProfile,
 *   overlapMetric,
 *   isModified,
 *   handleExponentChange,
 *   handleReset,
 * } = useExponentSliders({
 *   shellIndex: selectedShell,
 *   profileData,
 * });
 * ```
 */
export function useExponentSliders({
  shellIndex,
  profileData,
}: UseExponentSlidersParams): UseExponentSlidersReturn {
  // Store selectors
  const exponentOverrides = useBasisStore((state) => state.exponentOverrides);
  const setExponentOverride = useBasisStore((state) => state.setExponentOverride);
  const resetExponentOverrides = useBasisStore((state) => state.resetExponentOverrides);

  // Default exponents from the original profile
  const defaultExponents = useMemo(() => {
    if (!profileData) return [];
    return profileData.exponents;
  }, [profileData]);

  // Current exponents: overrides if present, else defaults
  const currentExponents = useMemo(() => {
    if (shellIndex === null || !profileData) return defaultExponents;
    const overrides = exponentOverrides[shellIndex];
    if (overrides && overrides.length === defaultExponents.length) {
      return overrides;
    }
    return defaultExponents;
  }, [shellIndex, profileData, exponentOverrides, defaultExponents]);

  // Whether any exponent is modified
  const isModified = useMemo(() => {
    if (defaultExponents.length === 0 || currentExponents.length === 0) return false;
    if (defaultExponents.length !== currentExponents.length) return false;
    return defaultExponents.some(
      (d, i) => Math.abs(d - currentExponents[i]) > 1e-12
    );
  }, [defaultExponents, currentExponents]);

  // Debounced profile recomputation state
  const [computedState, setComputedState] = useState<{
    modifiedProfile: number[] | null;
    overlapMetric: number;
  }>({ modifiedProfile: null, overlapMetric: 1.0 });

  const debounceTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  // Recompute modified profile when exponents change (debounced)
  useEffect(() => {
    if (debounceTimerRef.current) {
      clearTimeout(debounceTimerRef.current);
    }

    if (!isModified || !profileData || shellIndex === null) {
      setComputedState({ modifiedProfile: null, overlapMetric: 1.0 });
      return;
    }

    debounceTimerRef.current = setTimeout(() => {
      // Evaluate modified profile using the same r grid
      const modifiedValues = evaluateRadialProfileTS(
        currentExponents,
        profileData.rawCoefficients,
        profileData.angularMomentum,
        profileData.rValues,
      );

      // Compute overlap with original
      const S = computeOverlapMetric(
        profileData.contractedValues,
        modifiedValues,
        profileData.rValues,
      );

      setComputedState({
        modifiedProfile: modifiedValues,
        overlapMetric: S,
      });
    }, RECOMPUTE_DEBOUNCE_MS);

    return () => {
      if (debounceTimerRef.current) {
        clearTimeout(debounceTimerRef.current);
      }
    };
  }, [isModified, currentExponents, profileData, shellIndex]);

  // Handler for changing a single exponent
  const handleExponentChange = useCallback(
    (primitiveIndex: number, newAlpha: number) => {
      if (shellIndex === null) return;
      setExponentOverride(shellIndex, primitiveIndex, newAlpha, defaultExponents);
    },
    [shellIndex, setExponentOverride, defaultExponents],
  );

  // Handler for resetting all exponents
  const handleReset = useCallback(() => {
    if (shellIndex === null) return;
    resetExponentOverrides(shellIndex);
  }, [shellIndex, resetExponentOverrides]);

  return {
    defaultExponents,
    currentExponents,
    modifiedProfile: computedState.modifiedProfile,
    overlapMetric: computedState.overlapMetric,
    isModified,
    handleExponentChange,
    handleReset,
  };
}
