/**
 * Density Slice Extraction Utility
 *
 * Pure utility function that extracts a 2D slice from a cached 3D
 * density grid. Runs synchronously on the main thread -- no WASM
 * or worker call needed. For a 50x50x50 grid, a slice is just
 * 2500 array reads, completing in <1 ms.
 *
 * The grid uses C-order indexing: x-slowest, z-fastest.
 * Index into flat array: values[ix * ny * nz + iy * nz + iz]
 *
 * @module lib/densitySlice
 * @see US-062 Density Isosurface & Cross-Sections
 */

import type { PlaneAxis, DensitySliceResult } from '../types/density';

/**
 * Extract a 2D slice from a 3D density grid.
 *
 * Given the flat grid values (C-order: x-slowest, z-fastest) and grid
 * metadata, extracts a 2D array of density values at a specified plane.
 *
 * The plane position is snapped to the nearest grid point.
 *
 * @param values - Flat density grid values (C-order)
 * @param gridDims - Grid dimensions [nx, ny, nz]
 * @param gridOrigin - Grid origin [ox, oy, oz] in Bohr
 * @param gridSpacing - Uniform grid spacing in Bohr
 * @param planeAxis - Which plane to slice ('xy', 'xz', 'yz')
 * @param planePosition - Position along the perpendicular axis (Bohr)
 * @returns Slice result with 2D data and axis metadata
 */
export function extractDensitySlice(
  values: number[],
  gridDims: [number, number, number],
  gridOrigin: [number, number, number],
  gridSpacing: number,
  planeAxis: PlaneAxis,
  planePosition: number
): DensitySliceResult {
  const [nx, ny, nz] = gridDims;
  const [ox, oy, oz] = gridOrigin;

  switch (planeAxis) {
    case 'xy': {
      // Fixed z: plane contains x and y axes
      const iz = clampIndex(Math.round((planePosition - oz) / gridSpacing), nz);
      const actualPosition = oz + iz * gridSpacing;

      const xCoords = makeCoords(ox, nx, gridSpacing);
      const yCoords = makeCoords(oy, ny, gridSpacing);

      // data[iy][ix] = values[ix * ny * nz + iy * nz + iz]
      const data: number[][] = [];
      for (let iy = 0; iy < ny; iy++) {
        const row: number[] = [];
        for (let ix = 0; ix < nx; ix++) {
          row.push(values[ix * ny * nz + iy * nz + iz]);
        }
        data.push(row);
      }

      return {
        data,
        xCoords,
        yCoords,
        xLabel: 'x (bohr)',
        yLabel: 'y (bohr)',
        perpendicularIndex: iz,
        actualPosition,
      };
    }

    case 'xz': {
      // Fixed y: plane contains x and z axes
      const iy = clampIndex(Math.round((planePosition - oy) / gridSpacing), ny);
      const actualPosition = oy + iy * gridSpacing;

      const xCoords = makeCoords(ox, nx, gridSpacing);
      const zCoords = makeCoords(oz, nz, gridSpacing);

      // data[iz][ix] = values[ix * ny * nz + iy * nz + iz]
      const data: number[][] = [];
      for (let iz = 0; iz < nz; iz++) {
        const row: number[] = [];
        for (let ix = 0; ix < nx; ix++) {
          row.push(values[ix * ny * nz + iy * nz + iz]);
        }
        data.push(row);
      }

      return {
        data,
        xCoords,
        yCoords: zCoords,
        xLabel: 'x (bohr)',
        yLabel: 'z (bohr)',
        perpendicularIndex: iy,
        actualPosition,
      };
    }

    case 'yz': {
      // Fixed x: plane contains y and z axes
      const ix = clampIndex(Math.round((planePosition - ox) / gridSpacing), nx);
      const actualPosition = ox + ix * gridSpacing;

      const yCoords = makeCoords(oy, ny, gridSpacing);
      const zCoords = makeCoords(oz, nz, gridSpacing);

      // data[iz][iy] = values[ix * ny * nz + iy * nz + iz]
      const data: number[][] = [];
      for (let iz = 0; iz < nz; iz++) {
        const row: number[] = [];
        for (let iy = 0; iy < ny; iy++) {
          row.push(values[ix * ny * nz + iy * nz + iz]);
        }
        data.push(row);
      }

      return {
        data,
        xCoords: yCoords,
        yCoords: zCoords,
        xLabel: 'y (bohr)',
        yLabel: 'z (bohr)',
        perpendicularIndex: ix,
        actualPosition,
      };
    }
  }
}

/**
 * Clamp an index to [0, maxExclusive - 1].
 */
function clampIndex(idx: number, maxExclusive: number): number {
  return Math.max(0, Math.min(maxExclusive - 1, idx));
}

/**
 * Generate coordinate array for a grid axis.
 *
 * @param origin - Starting position in Bohr
 * @param n - Number of grid points along this axis
 * @param spacing - Grid spacing in Bohr
 * @returns Array of coordinate values
 */
function makeCoords(origin: number, n: number, spacing: number): number[] {
  const coords: number[] = [];
  for (let i = 0; i < n; i++) {
    coords.push(origin + i * spacing);
  }
  return coords;
}

/**
 * Get the range of the perpendicular axis for a given plane.
 *
 * Useful for setting slider bounds on the plane position control.
 *
 * @param planeAxis - Which plane ('xy', 'xz', 'yz')
 * @param gridOrigin - Grid origin [ox, oy, oz] in Bohr
 * @param gridDims - Grid dimensions [nx, ny, nz]
 * @param gridSpacing - Grid spacing in Bohr
 * @returns [min, max] range of the perpendicular axis in Bohr
 */
export function getPerpendicularRange(
  planeAxis: PlaneAxis,
  gridOrigin: [number, number, number],
  gridDims: [number, number, number],
  gridSpacing: number
): [number, number] {
  const [ox, oy, oz] = gridOrigin;
  const [nx, ny, nz] = gridDims;

  switch (planeAxis) {
    case 'xy':
      return [oz, oz + (nz - 1) * gridSpacing];
    case 'xz':
      return [oy, oy + (ny - 1) * gridSpacing];
    case 'yz':
      return [ox, ox + (nx - 1) * gridSpacing];
  }
}
