/**
 * usePesScan - React hook for PES scan computation.
 *
 * Provides explicit startScan() (legacy diatomic), startInternalScan()
 * (generalized internal coordinate), and cancelScan() methods for
 * controlling PES scans. Progress updates are streamed via the worker's
 * onProgress callback and stored in the scfStore's pesState.
 *
 * @module hooks/usePesScan
 * @see US-040 PES Scan UI (legacy diatomic)
 * @see US-082 Coordinate Selector + Scan Mode UI (internal coordinate)
 */

import { useCallback, useRef } from 'react';
import { useWorker, createRequestId } from './useWorker';
import { useScfStore } from '../stores/scfStore';
import type {
  PesScanRequest,
  PesScanResult,
  PesScanProgress,
  PesScanInternalRequest,
  PesScanInternalResult,
  PesScanInternalProgress,
  WorkerProgress,
  RequestId,
} from '../worker/protocol';
import type { CoordinateType } from '../components/scf/CoordinateSelector';
import type { ScanMode } from '../components/scf/ScanModeToggle';

// Re-export from shared utility for backward compatibility
export { ELEMENT_TO_Z } from '../lib/elements';

// ============================================================================
// Types
// ============================================================================

/**
 * Parameters for starting an internal coordinate PES scan.
 */
export interface InternalScanParams {
  atoms: [number, number, number, number][];
  basisName: string;
  method: string;
  coordinateType: CoordinateType;
  atomIndices: number[];
  scanMode: ScanMode;
  /** Min value (bohr for bonds, radians for angles/dihedrals) */
  valueMin: number;
  /** Max value (bohr for bonds, radians for angles/dihedrals) */
  valueMax: number;
  nPoints: number;
  useSpherical?: boolean;
  /** Max optimization steps per scan point for relaxed scans (default: 50) */
  optMaxSteps?: number;
  /** Gradient convergence threshold for relaxed scans in Ha/bohr (default: 4.5e-4) */
  optGradThreshold?: number;
}

/**
 * Return type of usePesScan hook.
 */
export interface UsePesScanResult {
  /** Whether the worker is ready to receive requests */
  isReady: boolean;
  /** Whether a PES scan is currently running */
  isScanning: boolean;
  /** Start a legacy diatomic PES scan with the current configuration */
  startScan: (atomAZ: number, atomBZ: number, basisName: string) => void;
  /** Start an internal coordinate PES scan */
  startInternalScan: (params: InternalScanParams) => void;
  /** Cancel the running PES scan */
  cancelScan: () => void;
}

/**
 * Type guard to check if a progress update is for legacy PES scan.
 */
function isPesProgress(progress: WorkerProgress): progress is PesScanProgress {
  return progress.module === 'pes';
}

/**
 * Type guard to check if a progress update is for internal coordinate PES scan.
 */
function isPesInternalProgress(progress: WorkerProgress): progress is PesScanInternalProgress {
  return progress.module === 'pes_internal';
}

// ============================================================================
// Hook Implementation
// ============================================================================

/**
 * React hook for PES scan computation.
 *
 * Supports two scan modes:
 * 1. Legacy diatomic scan via `startScan()` (PesScanRequest)
 * 2. Internal coordinate scan via `startInternalScan()` (PesScanInternalRequest)
 *
 * Both stream incremental progress and accumulate results in scfStore.pesState.
 *
 * @returns Object with isReady, isScanning, startScan, startInternalScan, cancelScan
 */
export function usePesScan(): UsePesScanResult {
  const { send, cancel, isReady } = useWorker();

  // Store selectors
  const scanning = useScfStore((state) => state.pesState.scanning);
  const rMin = useScfStore((state) => state.pesState.rMin);
  const rMax = useScfStore((state) => state.pesState.rMax);
  const nPoints = useScfStore((state) => state.pesState.nPoints);
  const convergenceProfile = useScfStore((state) => state.convergenceProfile);
  const maxIterations = useScfStore((state) => state.maxIterations);
  const useDiis = useScfStore((state) => state.useDiis);
  const damp = useScfStore((state) => state.damp);
  const scanRequestId = useScfStore((state) => state.pesState.scanRequestId);

  // Store actions
  const startPesScan = useScfStore((state) => state.startPesScan);
  const updatePesProgress = useScfStore((state) => state.updatePesProgress);
  const completePesScan = useScfStore((state) => state.completePesScan);
  const completePesInternalScan = useScfStore((state) => state.completePesInternalScan);
  const failPesScan = useScfStore((state) => state.failPesScan);
  const cancelPesScanStore = useScfStore((state) => state.cancelPesScan);

  // Ref to track current request ID for cleanup
  const requestIdRef = useRef<RequestId | null>(null);

  /**
   * Handle progress updates from the worker (legacy diatomic scan).
   */
  const handleProgress = useCallback(
    (progress: WorkerProgress) => {
      if (!isPesProgress(progress)) return;

      updatePesProgress({
        r: progress.r,
        energy: progress.energy,
        converged: progress.converged,
        iterations: 0,
      });
    },
    [updatePesProgress],
  );

  /**
   * Handle progress updates from the worker (internal coordinate scan).
   *
   * Maps internal coordinate progress to the PesPoint format used by
   * pesState.results so the PesCurvePlot can display points progressively.
   */
  const handleInternalProgress = useCallback(
    (progress: WorkerProgress) => {
      if (!isPesInternalProgress(progress)) return;

      updatePesProgress({
        r: progress.coordinateValue,
        energy: progress.energy,
        converged: progress.converged,
        iterations: 0,
      });
    },
    [updatePesProgress],
  );

  /**
   * Start a legacy diatomic PES scan.
   *
   * @param atomAZ - Atomic number of atom A (fixed at origin)
   * @param atomBZ - Atomic number of atom B (translated along z-axis)
   * @param basisName - Basis set name (e.g., "sto-3g")
   */
  const startScan = useCallback(
    async (atomAZ: number, atomBZ: number, basisName: string) => {
      if (!isReady || scanning) return;

      const requestId = createRequestId();
      requestIdRef.current = requestId;
      startPesScan(requestId);

      try {
        const request: Omit<PesScanRequest, 'requestId'> = {
          type: 'pes_scan',
          atomAZ,
          atomBZ,
          rMin,
          rMax,
          nPoints,
          basisName,
          options: {
            convergenceProfile,
            maxIterations,
            useDiis,
            damp: damp > 0 ? damp : undefined,
            includeMatrices: false,
          },
          useSeeding: true,
        };
        const result = await send<PesScanResult>(request, { onProgress: handleProgress });

        completePesScan(result);
      } catch (err) {
        const message = err instanceof Error ? err.message : 'PES scan failed';
        if (
          err instanceof Error &&
          (err.message.includes('cancelled') || err.message.includes('COMPUTATION_CANCELLED'))
        ) {
          cancelPesScanStore();
        } else {
          failPesScan(message);
        }
      } finally {
        requestIdRef.current = null;
      }
    },
    [
      isReady, scanning, rMin, rMax, nPoints, convergenceProfile, maxIterations,
      useDiis, damp, send, handleProgress, startPesScan, completePesScan,
      failPesScan, cancelPesScanStore,
    ],
  );

  /**
   * Start an internal coordinate PES scan.
   *
   * Sends a PesScanInternalRequest to the worker. Progress updates are
   * mapped to PesPoint format for compatibility with PesCurvePlot.
   *
   * When complete, the result is adapted to the PesScanResult format
   * expected by completePesScan().
   *
   * @param params - Internal scan parameters
   */
  const startInternalScan = useCallback(
    async (params: InternalScanParams) => {
      if (!isReady || scanning) return;

      const requestId = createRequestId();
      requestIdRef.current = requestId;
      startPesScan(requestId);

      try {
        const request: Omit<PesScanInternalRequest, 'requestId'> = {
          type: 'pes_scan_internal',
          atoms: params.atoms,
          basisName: params.basisName,
          method: params.method,
          coordinateType: params.coordinateType,
          atomIndices: params.atomIndices,
          scanMode: params.scanMode,
          valueMin: params.valueMin,
          valueMax: params.valueMax,
          nPoints: params.nPoints,
          useSeeding: true,
          useSpherical: params.useSpherical ?? true,
          convergenceProfile: 'tight',
          optMaxSteps: params.optMaxSteps,
          optGradThreshold: params.optGradThreshold,
        };

        const result = await send<PesScanInternalResult>(request, {
          onProgress: handleInternalProgress,
        });

        // Preserve full internal result for coordinate tracking (US-083)
        completePesInternalScan(result);

        // Adapt PesScanInternalResult to PesScanResult format for the store
        const adaptedResult: PesScanResult = {
          points: result.points.map((p) => ({
            r: p.coordinate_value,
            energy: p.energy,
            converged: p.converged,
            iterations: p.scf_iterations,
          })),
          equilibrium: result.equilibrium
            ? {
                r_bohr: result.equilibrium.value,
                energy_hartree: result.equilibrium.energy,
              }
            : null,
          compute_time_ms: (result as PesScanInternalResult & { compute_time_ms?: number }).compute_time_ms ?? 0,
          total_iterations: result.total_iterations,
        };

        completePesScan(adaptedResult);
      } catch (err) {
        const message = err instanceof Error ? err.message : 'PES scan failed';
        if (
          err instanceof Error &&
          (err.message.includes('cancelled') || err.message.includes('COMPUTATION_CANCELLED'))
        ) {
          cancelPesScanStore();
        } else {
          failPesScan(message);
        }
      } finally {
        requestIdRef.current = null;
      }
    },
    [
      isReady, scanning, send, handleInternalProgress, startPesScan,
      completePesScan, completePesInternalScan, failPesScan, cancelPesScanStore,
    ],
  );

  /**
   * Cancel the running PES scan.
   */
  const cancelScan = useCallback(() => {
    if (scanRequestId) {
      cancel(scanRequestId as RequestId);
      cancelPesScanStore();
    }
  }, [scanRequestId, cancel, cancelPesScanStore]);

  return {
    isReady,
    isScanning: scanning,
    startScan,
    startInternalScan,
    cancelScan,
  };
}
