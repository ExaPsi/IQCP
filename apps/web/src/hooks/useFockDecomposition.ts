/**
 * useFockDecomposition - React hook for computing the Fock matrix
 * decomposition via WASM.
 *
 * Sends a `fock_decomposition` request to the Web Worker, which
 * recomputes integrals (H^core, ERI) from the geometry and basis
 * set, then decomposes the Fock matrix into:
 *   F = H^core + G, where G = J - 0.5*K
 *
 * Unlike auto-computing hooks (useIntegralMatrices), this hook
 * exposes a manual `compute()` trigger because the Fock decomposition
 * requires a density matrix from a prior SCF run. The user must
 * click "Compute Fock Decomposition" to initiate.
 *
 * @module hooks/useFockDecomposition
 * @see US-058 Fock Build Tracing
 * @see Szabo & Ostlund (1996), Eq. 3.154
 */

import { useCallback, useRef } from 'react';
import { useWorker } from './useWorker';
import { useIntegralStore } from '../stores/integralStore';
import type { FockDecompositionResult } from '../stores/integralStore';
import type { GeometryInput, BasisSetName } from '../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * Parameters for the useFockDecomposition hook.
 */
export interface UseFockDecompositionParams {
  /** Molecule geometry (null if not available) */
  geometry: GeometryInput | null;
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** Density matrix from converged SCF (flat, nbf x nbf, includes factor of 2) */
  densityMatrix: number[] | null;
  /** Whether the worker is ready */
  isReady: boolean;
}

/**
 * Return type of useFockDecomposition hook.
 */
export interface UseFockDecompositionReturn {
  /** The decomposition result data, or null if not yet computed */
  data: FockDecompositionResult | null;
  /** Whether computation is in progress */
  loading: boolean;
  /** Error message if computation failed */
  error: string | null;
  /** Trigger the Fock decomposition computation */
  compute: () => void;
}

// ============================================================================
// Hook
// ============================================================================

/**
 * React hook for computing the Fock matrix decomposition.
 *
 * Does NOT auto-compute on mount or parameter change. Exposes a
 * `compute()` function triggered by user interaction (button click).
 *
 * Updates the integral store's `fockDecomposition` state on result/error.
 * Handles staleness: if parameters change while a request is in-flight,
 * the stale result is discarded.
 *
 * @param params - Hook parameters
 * @returns Object with loading state, error, data, and compute trigger
 *
 * @example
 * ```typescript
 * function FockPanel() {
 *   const { isReady } = useWorker();
 *   const { data, loading, error, compute } = useFockDecomposition({
 *     geometry: currentGeometry,
 *     basisName: 'sto-3g',
 *     densityMatrix: scfDensity,
 *     isReady,
 *   });
 *
 *   return (
 *     <div>
 *       <button onClick={compute} disabled={loading}>
 *         Compute Fock Decomposition
 *       </button>
 *       {data && <FockBuildTracer decomposition={data} />}
 *     </div>
 *   );
 * }
 * ```
 */
export function useFockDecomposition({
  geometry,
  basisName,
  densityMatrix,
  isReady,
}: UseFockDecompositionParams): UseFockDecompositionReturn {
  const { send } = useWorker();

  // Store selectors
  const fockDecomposition = useIntegralStore((s) => s.fockDecomposition);
  const setFockDecomposition = useIntegralStore((s) => s.setFockDecomposition);

  // Track the current request to avoid stale updates
  const requestRef = useRef(0);

  const compute = useCallback(() => {
    if (!isReady || !geometry || !densityMatrix || densityMatrix.length === 0) {
      return;
    }

    const requestId = ++requestRef.current;
    setFockDecomposition({ status: 'computing' });

    // Send the fock_decomposition request to the worker.
    // The worker handler and protocol types are being added by a
    // concurrent agent. Use a type assertion until the protocol
    // union is updated.
    const request = {
      type: 'fock_decomposition' as const,
      geometry,
      basisSet: basisName as BasisSetName,
      densityMatrix,
    };

    send<FockDecompositionResult>(request as never)
      .then((result) => {
        // Only update if this is still the latest request
        if (requestRef.current === requestId) {
          setFockDecomposition({ status: 'success', data: result });
        }
      })
      .catch((err) => {
        if (requestRef.current === requestId) {
          const message =
            err instanceof Error
              ? err.message
              : 'Failed to compute Fock decomposition';
          setFockDecomposition({ status: 'error', error: message });
        }
      });
  }, [isReady, geometry, basisName, densityMatrix, send, setFockDecomposition]);

  // Derive loading/error/data from the store state
  const loading = fockDecomposition.status === 'computing';
  const error = fockDecomposition.status === 'error' ? fockDecomposition.error : null;
  const data = fockDecomposition.status === 'success' ? fockDecomposition.data : null;

  return { data, loading, error, compute };
}
