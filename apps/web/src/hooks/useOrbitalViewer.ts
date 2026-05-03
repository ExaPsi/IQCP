/**
 * useOrbitalViewer - React hook for orbital isosurface computation.
 *
 * Orchestrates the two-step pipeline for orbital visualization:
 * 1. MO grid evaluation: evaluate the selected MO on a 3D grid
 * 2. Dual marching cubes: extract positive and negative lobe meshes
 *
 * The grid data is cached so that isovalue changes (step 2 only)
 * do not require recomputation of the grid (step 1). This makes
 * isovalue slider interaction fast (<=200ms target).
 *
 * Isovalue changes are debounced at 100ms per project convention.
 *
 * @module hooks/useOrbitalViewer
 * @see US-044 Orbital Viewer UI
 */

import { useState, useEffect, useRef, useCallback, useMemo } from 'react';
import { useWorker } from './useWorker';
import type {
  MoGridResult,
  DualMarchingCubesResult,
  DualMarchingCubesRequest,
  MoGridRequest,
} from '../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * Parameters for the useOrbitalViewer hook.
 */
export interface UseOrbitalViewerParams {
  /** Index of the selected MO (0-based), or null for none */
  selectedOrbital: number | null;
  /** Isovalue for isosurface extraction */
  isovalue: number;
  /** Full C matrix (nbf x nbf), column-major flat array from WASM */
  moCoefficients: number[] | null;
  /** Number of basis functions */
  nbf: number;
  /** Atom specifications: [Z, x, y, z] in Bohr */
  atoms: [number, number, number, number][] | null;
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** Whether the worker is ready to receive requests */
  isReady: boolean;
  /**
   * Whether the MO coefficients are in the spherical harmonic basis.
   * Must match the `useSpherical` option used in the SCF/integral computation.
   * @default false
   */
  useSpherical?: boolean;
}

/**
 * Return value of the useOrbitalViewer hook.
 */
export interface OrbitalViewerState {
  /** Whether a computation is in progress */
  loading: boolean;
  /** Error message from computation, or null */
  error: string | null;
  /** Cached grid data (reused across isovalue changes) */
  gridData: MoGridResult | null;
  /** Current isosurface mesh data (positive and negative lobes) */
  meshData: DualMarchingCubesResult | null;
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
 * Grid spacing in Bohr for MO evaluation.
 *
 * 0.2 bohr gives reasonable visual quality for isosurfaces.
 * Finer grids are more accurate but slower (grid points scale as 1/h^3).
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
 * Extract a single MO coefficient vector from the full C matrix.
 *
 * The C matrix from Rust/WASM is stored in column-major order
 * (nalgebra convention). Column i (MO i) consists of elements
 * at indices i*nbf through (i+1)*nbf - 1 in the flat array.
 *
 * @param cMatrix - Full C matrix as a flat column-major array
 * @param orbitalIndex - 0-based index of the MO to extract
 * @param nbf - Number of basis functions
 * @returns MO coefficient vector of length nbf
 */
function extractMoCoefficients(
  cMatrix: number[],
  orbitalIndex: number,
  nbf: number
): number[] {
  const start = orbitalIndex * nbf;
  const end = start + nbf;
  return cMatrix.slice(start, end);
}

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
    // Fallback for empty molecule (shouldn't happen in practice)
    return {
      origin: [-GRID_PADDING_BOHR, -GRID_PADDING_BOHR, -GRID_PADDING_BOHR],
      spacing: GRID_SPACING_BOHR,
      dims: [2, 2, 2],
    };
  }

  // Find bounding box of atom positions
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

  // Extend by padding
  const ox = xMin - GRID_PADDING_BOHR;
  const oy = yMin - GRID_PADDING_BOHR;
  const oz = zMin - GRID_PADDING_BOHR;

  const lx = (xMax + GRID_PADDING_BOHR) - ox;
  const ly = (yMax + GRID_PADDING_BOHR) - oy;
  const lz = (zMax + GRID_PADDING_BOHR) - oz;

  // Compute grid dimensions (at least 2 points per dimension)
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
 * useOrbitalViewer - orchestrates MO grid evaluation and marching cubes.
 *
 * Pipeline:
 * 1. When selectedOrbital changes: extract MO coefficients, compute grid
 *    bounds, send MoGridRequest, cache result, then run marching cubes.
 * 2. When isovalue changes (with cached grid): skip grid evaluation,
 *    send DualMarchingCubesRequest with cached grid data.
 *
 * @param params - Hook parameters
 * @returns Current orbital viewer state
 *
 * @example
 * ```typescript
 * const { loading, error, meshData } = useOrbitalViewer({
 *   selectedOrbital: 4,
 *   isovalue: 0.03,
 *   moCoefficients: result.matrices.moCoefficients,
 *   nbf: 7,
 *   atoms: [[8, 0, 0, 0.22], [1, 0, 1.43, -0.88], [1, 0, -1.43, -0.88]],
 *   basisName: 'sto-3g',
 *   isReady: true,
 * });
 * ```
 */
export function useOrbitalViewer(params: UseOrbitalViewerParams): OrbitalViewerState {
  const {
    selectedOrbital,
    isovalue,
    moCoefficients,
    nbf,
    atoms,
    basisName,
    isReady,
    useSpherical = false,
  } = params;

  const { send } = useWorker();

  // State
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [gridData, setGridData] = useState<MoGridResult | null>(null);
  const [meshData, setMeshData] = useState<DualMarchingCubesResult | null>(null);

  // Track which orbital the current grid cache is for
  const cachedOrbitalRef = useRef<number | null>(null);

  // Abort controller for cancelling in-flight requests on parameter changes
  const abortRef = useRef(0);

  // Debounce timer for isovalue changes
  const isovalueTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  // ============================================================================
  // Marching Cubes Step (runs with cached grid data)
  // ============================================================================

  /**
   * Run dual marching cubes on cached grid data.
   *
   * This is the fast path used when only the isovalue changes.
   */
  const runMarchingCubes = useCallback(
    async (grid: MoGridResult, iso: number, generation: number) => {
      try {
        // Type assertion needed because TypeScript can't narrow Omit<UnionType, 'field'>
        const mcRequest: Omit<DualMarchingCubesRequest, 'requestId'> = {
          type: 'dual_marching_cubes',
          gridData: grid.values,
          gridDims: grid.gridDims,
          gridOrigin: grid.gridOrigin,
          gridSpacing: grid.gridSpacing,
          isovalue: iso,
        };
        const result = await send<DualMarchingCubesResult>(mcRequest);

        // Only apply if this is still the current generation
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
  // Full Pipeline (grid eval + marching cubes)
  // ============================================================================

  /**
   * Run the full pipeline: MO grid evaluation followed by marching cubes.
   *
   * Used when the selected orbital changes and the grid cache is invalidated.
   */
  const runFullPipeline = useCallback(
    async (orbitalIdx: number, iso: number, generation: number) => {
      if (!moCoefficients || !atoms || atoms.length === 0 || nbf <= 0) {
        return;
      }

      setLoading(true);
      setError(null);

      try {
        // 1. Extract the MO coefficient vector for the selected orbital
        const moCoeffs = extractMoCoefficients(moCoefficients, orbitalIdx, nbf);

        // 2. Compute grid bounds from atom positions
        const { origin, spacing, dims } = computeGridBounds(atoms);

        // 3. Send MoGridRequest to worker
        // Type assertion needed because TypeScript can't narrow Omit<UnionType, 'field'>
        const gridRequest: Omit<MoGridRequest, 'requestId'> = {
          type: 'mo_grid',
          moCoefficients: moCoeffs,
          atoms,
          basisName,
          gridOrigin: origin,
          gridSpacing: spacing,
          gridDims: dims,
          useSpherical,
        };
        const gridResult = await send<MoGridResult>(gridRequest);

        // Check if we've been superseded
        if (abortRef.current !== generation) return;

        // 4. Cache the grid result
        setGridData(gridResult);
        cachedOrbitalRef.current = orbitalIdx;

        // 5. Run marching cubes on the new grid
        await runMarchingCubes(gridResult, iso, generation);
      } catch (err) {
        if (abortRef.current === generation) {
          setError(err instanceof Error ? err.message : 'Orbital grid evaluation failed');
          setLoading(false);
        }
      }
    },
    [moCoefficients, atoms, nbf, basisName, useSpherical, send, runMarchingCubes]
  );

  // ============================================================================
  // Effect: Respond to selectedOrbital changes (full pipeline)
  // ============================================================================

  useEffect(() => {
    // Clear any pending isovalue debounce
    if (isovalueTimerRef.current) {
      clearTimeout(isovalueTimerRef.current);
      isovalueTimerRef.current = null;
    }

    if (selectedOrbital === null || !isReady || !moCoefficients || !atoms) {
      // No orbital selected or not ready -- clear state
      setGridData(null);
      setMeshData(null);
      setError(null);
      setLoading(false);
      cachedOrbitalRef.current = null;
      return;
    }

    // Validate orbital index
    if (selectedOrbital < 0 || selectedOrbital >= nbf) {
      setError(`Invalid orbital index: ${selectedOrbital} (nbf=${nbf})`);
      return;
    }

    // Increment generation to invalidate in-flight requests
    const generation = ++abortRef.current;

    // Run the full pipeline
    runFullPipeline(selectedOrbital, isovalue, generation);

    // Note: we do NOT include isovalue in the dependency array here.
    // Isovalue-only changes are handled by the effect below.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedOrbital, isReady, moCoefficients, atoms, nbf, basisName, runFullPipeline]);

  // ============================================================================
  // Effect: Respond to isovalue changes (marching cubes only, debounced)
  // ============================================================================

  useEffect(() => {
    // Only run if we have cached grid data for the current orbital
    if (
      gridData === null ||
      selectedOrbital === null ||
      cachedOrbitalRef.current !== selectedOrbital
    ) {
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
  }, [isovalue, gridData, selectedOrbital, runMarchingCubes]);

  // ============================================================================
  // Cleanup on unmount
  // ============================================================================

  useEffect(() => {
    // Capture refs for cleanup (satisfies react-hooks/exhaustive-deps)
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
    () => ({ loading, error, gridData, meshData }),
    [loading, error, gridData, meshData]
  );
}
