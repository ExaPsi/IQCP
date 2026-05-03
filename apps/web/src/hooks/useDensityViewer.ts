/**
 * useDensityViewer - React hook for density isosurface computation.
 *
 * Orchestrates the two-step pipeline for density visualization:
 * 1. Density grid evaluation: compute electron density on a 3D grid
 * 2. Marching cubes: extract single isosurface mesh (density >= 0)
 *
 * For difference density mode (US-063):
 * 1. Density grid evaluation (same as total mode -- may use cached grid)
 * 2. Difference density request: compute Delta-rho = rho_mol - rho_promolecule
 * 3. Dual marching cubes: extract accumulation (+) and depletion (-) isosurfaces
 *
 * Grid data is cached so that isovalue changes (step 2/3 only)
 * and cross-section slicing do not require recomputation of step 1.
 * This makes isovalue slider interaction fast (<=200ms target).
 *
 * Unlike useOrbitalViewer:
 * - Uses `density_grid` request (not `mo_grid`)
 * - Total mode: single `marching_cubes` (density is always non-negative)
 * - Difference mode: `dual_marching_cubes` for +/- isosurfaces
 * - Exposes raw gridData for client-side cross-section slicing
 *
 * Isovalue changes are debounced at 100ms per project convention.
 *
 * @module hooks/useDensityViewer
 * @see US-062 Density Isosurface & Cross-Sections
 * @see US-063 Difference Density
 */

import { useState, useEffect, useRef, useCallback, useMemo } from 'react';
import { useWorker } from './useWorker';
import type {
  DensityGridResult,
  MarchingCubesResult,
  MarchingCubesRequest,
  DensityGridRequest,
  DifferenceDensityRequest,
  DifferenceDensityResult,
  DualMarchingCubesRequest,
  DualMarchingCubesResult,
} from '../worker/protocol';
import type { DensityMode } from '../types/density';

// ============================================================================
// Types
// ============================================================================

/**
 * Parameters for the useDensityViewer hook.
 */
export interface UseDensityViewerParams {
  /** Whether density visualization is active (user has opened density tab) */
  active: boolean;
  /** Density matrix from converged SCF (flat, row-major) */
  densityMatrix: number[] | null;
  /** Atom specifications for grid bounds: [Z, x, y, z] */
  atoms: [number, number, number, number][] | null;
  /** Basis set name */
  basisName: string;
  /** Number of electrons */
  nElectrons: number;
  /** Whether spherical harmonics were used */
  useSpherical: boolean;
  /** Isovalue for total density isosurface extraction */
  isovalue: number;
  /** Density visualization mode: 'total' or 'difference' */
  densityMode: DensityMode;
  /** Isovalue for difference density isosurface extraction (US-063) */
  diffIsovalue: number;
  /** Whether the worker is ready */
  isReady: boolean;
}

/**
 * Return value of the useDensityViewer hook.
 */
export interface DensityViewerState {
  /** Whether a computation is in progress */
  loading: boolean;
  /** Error message, or null */
  error: string | null;
  /** Cached total density grid data (reused across isovalue/slice changes) */
  gridData: DensityGridResult | null;
  /** Current total density isosurface mesh data (single surface) */
  meshData: MarchingCubesResult | null;
  /** Cached difference density grid data (US-063) */
  diffGridData: DifferenceDensityResult | null;
  /** Current difference density dual isosurface mesh (accumulation/depletion) (US-063) */
  diffMeshData: DualMarchingCubesResult | null;
}

// ============================================================================
// Constants
// ============================================================================

/**
 * Padding in Bohr beyond outermost atom position for grid bounds.
 *
 * 6 bohr provides sufficient spatial extent for STO-3G and small
 * split-valence basis sets used in the educational context.
 */
const GRID_PADDING_BOHR = 6.0;

/**
 * Grid spacing in Bohr for density evaluation.
 *
 * 0.2 bohr gives reasonable visual quality for isosurfaces.
 */
const GRID_SPACING_BOHR = 0.2;

/**
 * Debounce delay for isovalue changes (ms).
 */
const ISOVALUE_DEBOUNCE_MS = 100;

// ============================================================================
// Helpers
// ============================================================================

/**
 * Compute auto-grid bounds from atom positions.
 *
 * Extends GRID_PADDING_BOHR beyond the outermost atom position in each
 * direction, with uniform spacing of GRID_SPACING_BOHR.
 *
 * @param atoms - Atom positions as [Z, x, y, z] arrays
 * @returns Grid parameters: origin, spacing, and dimensions
 */
function computeGridBounds(atoms: [number, number, number, number][]): {
  origin: [number, number, number];
  spacing: number;
  dims: [number, number, number];
} {
  if (atoms.length === 0) {
    return {
      origin: [-GRID_PADDING_BOHR, -GRID_PADDING_BOHR, -GRID_PADDING_BOHR],
      spacing: GRID_SPACING_BOHR,
      dims: [2, 2, 2],
    };
  }

  let xMin = Infinity, xMax = -Infinity;
  let yMin = Infinity, yMax = -Infinity;
  let zMin = Infinity, zMax = -Infinity;

  for (const [, x, y, z] of atoms) {
    if (x < xMin) xMin = x;
    if (x > xMax) xMax = x;
    if (y < yMin) yMin = y;
    if (y > yMax) yMax = y;
    if (z < zMin) zMin = z;
    if (z > zMax) zMax = z;
  }

  const ox = xMin - GRID_PADDING_BOHR;
  const oy = yMin - GRID_PADDING_BOHR;
  const oz = zMin - GRID_PADDING_BOHR;

  const lx = (xMax + GRID_PADDING_BOHR) - ox;
  const ly = (yMax + GRID_PADDING_BOHR) - oy;
  const lz = (zMax + GRID_PADDING_BOHR) - oz;

  const nx = Math.max(2, Math.ceil(lx / GRID_SPACING_BOHR) + 1);
  const ny = Math.max(2, Math.ceil(ly / GRID_SPACING_BOHR) + 1);
  const nz = Math.max(2, Math.ceil(lz / GRID_SPACING_BOHR) + 1);

  return {
    origin: [ox, oy, oz],
    spacing: GRID_SPACING_BOHR,
    dims: [nx, ny, nz],
  };
}

// ============================================================================
// Hook
// ============================================================================

/**
 * useDensityViewer - orchestrates density grid evaluation and marching cubes.
 *
 * Pipeline (total mode):
 * 1. When active becomes true (or densityMatrix changes): compute density
 *    grid via worker, cache result, then run marching cubes.
 * 2. When isovalue changes (with cached grid): skip grid evaluation,
 *    send MarchingCubesRequest with cached grid data.
 * 3. Cross-section slicing uses gridData directly -- no worker call.
 *
 * Pipeline (difference mode, US-063):
 * 1. Compute total density grid (same as above, may use cache).
 * 2. Send difference_density request with cached total density values.
 * 3. Run dual_marching_cubes on the Delta-rho grid for +/- surfaces.
 * 4. When diffIsovalue changes, re-run dual marching cubes on cached diff grid.
 *
 * @param params - Hook parameters
 * @returns Current density viewer state
 */
export function useDensityViewer(params: UseDensityViewerParams): DensityViewerState {
  const {
    active,
    densityMatrix,
    atoms,
    basisName,
    nElectrons,
    useSpherical,
    isovalue,
    densityMode,
    diffIsovalue,
    isReady,
  } = params;

  const { send } = useWorker();

  // State -- total density
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [gridData, setGridData] = useState<DensityGridResult | null>(null);
  const [meshData, setMeshData] = useState<MarchingCubesResult | null>(null);

  // State -- difference density (US-063)
  const [diffGridData, setDiffGridData] = useState<DifferenceDensityResult | null>(null);
  const [diffMeshData, setDiffMeshData] = useState<DualMarchingCubesResult | null>(null);

  // Track which density matrix the current grid cache is for
  // (use a fingerprint to detect changes without deep comparison)
  const cachedDensityRef = useRef<number[] | null>(null);

  // Abort controller for cancelling in-flight requests on parameter changes
  const abortRef = useRef(0);

  // Debounce timer for isovalue changes
  const isovalueTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  // ============================================================================
  // Total Density: Marching Cubes Step (runs with cached grid data)
  // ============================================================================

  /**
   * Run single marching cubes on cached grid data.
   *
   * This is the fast path used when only the isovalue changes.
   * Unlike orbital viewer, we use single marching cubes (not dual)
   * because density is always non-negative.
   */
  const runMarchingCubes = useCallback(
    async (grid: DensityGridResult, iso: number, generation: number) => {
      try {
        const mcRequest: Omit<MarchingCubesRequest, 'requestId'> = {
          type: 'marching_cubes',
          gridData: grid.values,
          gridDims: grid.gridDims,
          gridOrigin: grid.gridOrigin,
          gridSpacing: grid.gridSpacing,
          isovalue: iso,
        };
        const result = await send<MarchingCubesResult>(mcRequest);

        if (abortRef.current === generation) {
          setMeshData(result);
          setLoading(false);
          setError(null);
        }
      } catch (err) {
        if (abortRef.current === generation) {
          setError(err instanceof Error ? err.message : 'Marching cubes failed');
          setLoading(false);
        }
      }
    },
    [send]
  );

  // ============================================================================
  // Difference Density: Dual Marching Cubes Step (US-063)
  // ============================================================================

  /**
   * Run dual marching cubes on cached difference density grid data.
   *
   * This is the fast path used when only the diffIsovalue changes.
   * Uses dual_marching_cubes because Delta-rho has both positive
   * (accumulation) and negative (depletion) values.
   */
  const runDiffMarchingCubes = useCallback(
    async (diffGrid: DifferenceDensityResult, iso: number, generation: number) => {
      try {
        const dmcRequest: Omit<DualMarchingCubesRequest, 'requestId'> = {
          type: 'dual_marching_cubes',
          gridData: diffGrid.values,
          gridDims: diffGrid.gridDims,
          gridOrigin: diffGrid.gridOrigin,
          gridSpacing: diffGrid.gridSpacing,
          isovalue: iso,
        };
        const result = await send<DualMarchingCubesResult>(dmcRequest);

        if (abortRef.current === generation) {
          setDiffMeshData(result);
          setLoading(false);
          setError(null);
        }
      } catch (err) {
        if (abortRef.current === generation) {
          setError(err instanceof Error ? err.message : 'Dual marching cubes failed');
          setLoading(false);
        }
      }
    },
    [send]
  );

  // ============================================================================
  // Difference Density Pipeline (US-063)
  // ============================================================================

  /**
   * Run the difference density pipeline:
   * 1. Use cached total density grid (must already exist)
   * 2. Compute difference density via WASM
   * 3. Run dual marching cubes for +/- isosurfaces
   */
  const runDiffPipeline = useCallback(
    async (totalGrid: DensityGridResult, iso: number, generation: number) => {
      if (!atoms || atoms.length === 0) return;

      try {
        // Send difference_density request with cached total density
        const diffRequest: Omit<DifferenceDensityRequest, 'requestId'> = {
          type: 'difference_density',
          totalDensity: totalGrid.values,
          atoms,
          gridOrigin: totalGrid.gridOrigin,
          gridSpacing: totalGrid.gridSpacing,
          gridDims: totalGrid.gridDims,
        };
        const diffResult = await send<DifferenceDensityResult>(diffRequest);

        // Check if we've been superseded
        if (abortRef.current !== generation) return;

        // Cache the difference density grid
        setDiffGridData(diffResult);

        // Run dual marching cubes for +/- isosurfaces
        await runDiffMarchingCubes(diffResult, iso, generation);
      } catch (err) {
        if (abortRef.current === generation) {
          setError(err instanceof Error ? err.message : 'Difference density evaluation failed');
          setLoading(false);
        }
      }
    },
    [atoms, send, runDiffMarchingCubes]
  );

  // ============================================================================
  // Full Pipeline (density grid eval + marching cubes or diff density)
  // ============================================================================

  /**
   * Run the full pipeline: density grid evaluation followed by marching cubes
   * (total mode) or difference density + dual marching cubes (difference mode).
   *
   * Used when the density matrix changes and the grid cache is invalidated.
   */
  const runFullPipeline = useCallback(
    async (dMatrix: number[], iso: number, mode: DensityMode, diffIso: number, generation: number) => {
      if (!atoms || atoms.length === 0) {
        return;
      }

      setLoading(true);
      setError(null);

      try {
        // 1. Compute grid bounds from atom positions
        const { origin, spacing, dims } = computeGridBounds(atoms);

        // 2. Send DensityGridRequest to worker (always needed)
        const gridRequest: Omit<DensityGridRequest, 'requestId'> = {
          type: 'density_grid',
          densityMatrix: dMatrix,
          atoms,
          basisName,
          gridOrigin: origin,
          gridSpacing: spacing,
          gridDims: dims,
          nElectrons,
          useSpherical,
        };
        const gridResult = await send<DensityGridResult>(gridRequest);

        // Check if we've been superseded
        if (abortRef.current !== generation) return;

        // 3. Cache the grid result
        setGridData(gridResult);
        cachedDensityRef.current = dMatrix;

        // 4. Depending on mode, run appropriate pipeline
        if (mode === 'difference') {
          await runDiffPipeline(gridResult, diffIso, generation);
        } else {
          await runMarchingCubes(gridResult, iso, generation);
        }
      } catch (err) {
        if (abortRef.current === generation) {
          setError(err instanceof Error ? err.message : 'Density grid evaluation failed');
          setLoading(false);
        }
      }
    },
    [atoms, basisName, nElectrons, useSpherical, send, runMarchingCubes, runDiffPipeline]
  );

  // ============================================================================
  // Effect: Respond to activation / density matrix changes (full pipeline)
  // ============================================================================

  useEffect(() => {
    // Clear any pending isovalue debounce
    if (isovalueTimerRef.current) {
      clearTimeout(isovalueTimerRef.current);
      isovalueTimerRef.current = null;
    }

    if (!active || !isReady || !densityMatrix || !atoms) {
      // Not active or not ready -- clear state
      setGridData(null);
      setMeshData(null);
      setDiffGridData(null);
      setDiffMeshData(null);
      setError(null);
      setLoading(false);
      cachedDensityRef.current = null;
      return;
    }

    // Check if density matrix has changed (reference equality)
    if (cachedDensityRef.current === densityMatrix && gridData !== null) {
      // Grid is still valid -- no need for full pipeline
      return;
    }

    // Increment generation to invalidate in-flight requests
    const generation = ++abortRef.current;

    // Run the full pipeline
    runFullPipeline(densityMatrix, isovalue, densityMode, diffIsovalue, generation);

    // Note: isovalue/diffIsovalue are NOT in the dependency array here.
    // Isovalue-only changes are handled by the effects below.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [active, isReady, densityMatrix, atoms, basisName, runFullPipeline]);

  // ============================================================================
  // Effect: Respond to density mode switch (with cached total grid)
  // ============================================================================

  useEffect(() => {
    // Only run if we have a cached total grid and the mode just changed
    if (!active || gridData === null) return;

    // Clear stale mesh data from the other mode
    if (densityMode === 'difference') {
      // Switching to difference mode: need to compute diff grid if not cached
      if (diffGridData === null && atoms) {
        const generation = ++abortRef.current;
        setLoading(true);
        runDiffPipeline(gridData, diffIsovalue, generation);
      }
    } else {
      // Switching to total mode: need to re-run marching cubes if mesh is stale
      if (meshData === null) {
        const generation = ++abortRef.current;
        setLoading(true);
        runMarchingCubes(gridData, isovalue, generation);
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [densityMode]);

  // ============================================================================
  // Effect: Respond to total density isovalue changes (marching cubes only)
  // ============================================================================

  useEffect(() => {
    // Only run if we have cached grid data and in total mode
    if (gridData === null || !active || densityMode !== 'total') {
      return;
    }

    // Clear any previous debounce timer
    if (isovalueTimerRef.current) {
      clearTimeout(isovalueTimerRef.current);
    }

    // Debounce isovalue changes
    isovalueTimerRef.current = setTimeout(() => {
      isovalueTimerRef.current = null;
      const generation = ++abortRef.current;
      setLoading(true);
      runMarchingCubes(gridData, isovalue, generation);
    }, ISOVALUE_DEBOUNCE_MS);

    return () => {
      if (isovalueTimerRef.current) {
        clearTimeout(isovalueTimerRef.current);
        isovalueTimerRef.current = null;
      }
    };
  }, [isovalue, gridData, active, densityMode, runMarchingCubes]);

  // ============================================================================
  // Effect: Respond to diff isovalue changes (dual marching cubes only, US-063)
  // ============================================================================

  useEffect(() => {
    // Only run if we have cached diff grid data and in difference mode
    if (diffGridData === null || !active || densityMode !== 'difference') {
      return;
    }

    // Clear any previous debounce timer
    if (isovalueTimerRef.current) {
      clearTimeout(isovalueTimerRef.current);
    }

    // Debounce diff isovalue changes
    isovalueTimerRef.current = setTimeout(() => {
      isovalueTimerRef.current = null;
      const generation = ++abortRef.current;
      setLoading(true);
      runDiffMarchingCubes(diffGridData, diffIsovalue, generation);
    }, ISOVALUE_DEBOUNCE_MS);

    return () => {
      if (isovalueTimerRef.current) {
        clearTimeout(isovalueTimerRef.current);
        isovalueTimerRef.current = null;
      }
    };
  }, [diffIsovalue, diffGridData, active, densityMode, runDiffMarchingCubes]);

  // ============================================================================
  // Cleanup on unmount
  // ============================================================================

  useEffect(() => {
    const abortControl = abortRef;
    const timerControl = isovalueTimerRef;

    return () => {
      // Invalidate any in-flight requests
      abortControl.current++;

      // Clear debounce timer
      if (timerControl.current) {
        clearTimeout(timerControl.current);
        timerControl.current = null;
      }
    };
  }, []);

  // ============================================================================
  // Return
  // ============================================================================

  return useMemo(
    () => ({ loading, error, gridData, meshData, diffGridData, diffMeshData }),
    [loading, error, gridData, meshData, diffGridData, diffMeshData]
  );
}
