/**
 * useIntegralMatrices - React hook for computing integral matrices via WASM.
 *
 * Sends an `integral_matrices` request to the Web Worker when the system
 * selection changes and stores the result in the integral store.
 *
 * After one-electron integrals (S, T, V, Hcore) succeed, fires a background
 * `integral_compute` request to compute and cache ERI data for use by
 * US-058 (Fock Build Tracing) and US-059 (ERI Browser).
 *
 * Auto-computes when inputs change -- one-electron integral computation
 * is fast enough (sub-second for all preset systems) that no manual
 * "Compute" button is needed.
 *
 * @module hooks/useIntegralMatrices
 * @see US-055 Integral Matrix Heatmap
 * @see US-057 System Selection (AC3 - ERI caching)
 */

import { useEffect, useRef, useCallback } from 'react';
import { useWorker } from './useWorker';
import { useIntegralStore } from '../stores/integralStore';
import { MOLECULE_PRESETS, getPresetIdForMoleculeBasis } from '../stores/scfStore';
import type { IntegralMatricesResult } from '../stores/integralStore';
import type { GeometryInput, BasisSetName, IntegralComputeResult } from '../worker/protocol';

// ============================================================================
// Preset System Geometries
// ============================================================================

/**
 * Geometry definitions for preset systems.
 *
 * These match the systems in content/presets/systems.json.
 * Coordinates are in Bohr.
 */
/**
 * Preset geometries generated from the shared MOLECULE_PRESETS.
 *
 * Includes ALL molecules available in the SCF Sandbox (Module E),
 * each paired with STO-3G as the default basis set. Keyed by the
 * pre-computed preset ID when available, otherwise by molecule ID.
 */
export const PRESET_GEOMETRIES: Record<string, { geometry: GeometryInput; basisName: BasisSetName }> =
  Object.fromEntries(
    MOLECULE_PRESETS.map((mol) => {
      const presetId = getPresetIdForMoleculeBasis(mol.id, 'sto-3g');
      const key = presetId ?? mol.id;
      return [
        key,
        {
          geometry: {
            atoms: mol.atoms.map((a) => ({ symbol: a.symbol, xyz: a.xyz as [number, number, number] })),
            units: 'bohr' as const,
          },
          basisName: 'sto-3g' as BasisSetName,
        },
      ];
    })
  );

// ============================================================================
// Hook
// ============================================================================

/**
 * React hook for computing integral matrices from WASM.
 *
 * Watches the integral store for system/basis changes and automatically
 * sends computation requests to the worker. Updates the store with results.
 *
 * After one-electron integrals succeed, automatically fires a background
 * `integral_compute` request to compute and cache ERI data. The ERI
 * request uses a separate staleness counter (`eriRequestRef`) to handle
 * the case where the system changes while ERI is still computing.
 *
 * @example
 * ```typescript
 * function IntegralPanel() {
 *   useIntegralMatrices();
 *   const compute = useIntegralStore(s => s.compute);
 *   const eriCache = useIntegralStore(s => s.eriCache);
 *
 *   if (compute.status === 'computing') return <p>Computing...</p>;
 *   if (compute.status === 'success') return <Heatmap data={compute.data} />;
 * }
 * ```
 */
export function useIntegralMatrices(): void {
  const { send, isReady } = useWorker();

  // Store selectors
  const systemId = useIntegralStore((s) => s.systemId);
  const inputMode = useIntegralStore((s) => s.inputMode);
  const customGeometry = useIntegralStore((s) => s.customGeometry);
  const customBasis = useIntegralStore((s) => s.customBasis);
  const useSpherical = useIntegralStore((s) => s.useSpherical);
  const setCompute = useIntegralStore((s) => s.setCompute);
  const setEriCache = useIntegralStore((s) => s.setEriCache);

  // Track the current request to avoid stale updates (one-electron integrals)
  const requestRef = useRef(0);
  // Separate staleness counter for ERI background requests
  const eriRequestRef = useRef(0);

  /**
   * Compute ERI in the background using the existing `integral_compute`
   * worker request. This request returns S, Hcore, AND eriCompressed.
   * We only need the ERI data from the result.
   */
  const computeEriBackground = useCallback(
    async (geometry: GeometryInput, basisName: BasisSetName, spherical: boolean) => {
      const eriId = ++eriRequestRef.current;
      setEriCache({ status: 'computing' });

      try {
        const request = {
          type: 'integral_compute' as const,
          geometry,
          basisSet: basisName,
          useSpherical: spherical,
        };
        const result = await send<IntegralComputeResult>(request as never);

        // Only update if this is still the latest ERI request
        if (eriRequestRef.current === eriId) {
          setEriCache({
            status: 'success',
            eriCompressed: result.eriCompressed,
            eriIndexing: result.eriIndexing,
            eriCount: result.eriCompressed.length,
            computeTimeMs: result.metadata.computeTimeMs,
          });
        }
      } catch (err) {
        if (eriRequestRef.current === eriId) {
          const message =
            err instanceof Error ? err.message : 'Failed to compute ERI';
          setEriCache({ status: 'error', error: message });
        }
      }
    },
    [send, setEriCache],
  );

  const computeMatrices = useCallback(
    async (geometry: GeometryInput, basisName: BasisSetName, spherical: boolean) => {
      const requestId = ++requestRef.current;
      // Also invalidate any in-flight ERI request
      ++eriRequestRef.current;
      setCompute({ status: 'computing' });
      setEriCache({ status: 'idle' });

      try {
        const request = {
          type: 'integral_matrices' as const,
          geometry,
          basisName,
          useSpherical: spherical,
        };
        const result = await send<IntegralMatricesResult>(request as never);

        // Only update if this is still the latest request
        if (requestRef.current === requestId) {
          setCompute({ status: 'success', data: result });

          // Fire background ERI computation now that one-electron integrals succeeded.
          // This runs asynchronously and does not block the heatmap display.
          void computeEriBackground(geometry, basisName, spherical);
        }
      } catch (err) {
        if (requestRef.current === requestId) {
          const message =
            err instanceof Error ? err.message : 'Failed to compute integral matrices';
          setCompute({ status: 'error', error: message });
        }
      }
    },
    [send, setCompute, setEriCache, computeEriBackground],
  );

  useEffect(() => {
    if (!isReady) return;

    if (inputMode === 'preset') {
      const preset = PRESET_GEOMETRIES[systemId];
      if (!preset) {
        setCompute({ status: 'error', error: `Unknown preset system: ${systemId}` });
        return;
      }
      computeMatrices(preset.geometry, preset.basisName, useSpherical);
    } else if (inputMode === 'custom' && customGeometry) {
      computeMatrices(customGeometry, customBasis, useSpherical);
    }
  }, [isReady, inputMode, systemId, customGeometry, customBasis, useSpherical, computeMatrices, setCompute]);
}
