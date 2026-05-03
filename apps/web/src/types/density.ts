/**
 * Density Visualization Type Definitions
 *
 * Shared types for the density isosurface and cross-section
 * visualization components (US-062).
 *
 * @module types/density
 */

/**
 * Axis for the cross-section plane.
 *
 * The plane is perpendicular to the named axis:
 * - 'xy': z = constant (plane contains x and y axes)
 * - 'xz': y = constant (plane contains x and z axes)
 * - 'yz': x = constant (plane contains y and z axes)
 */
export type PlaneAxis = 'xy' | 'xz' | 'yz';

/**
 * Color scale mode for density cross-section.
 */
export type ColorScale = 'linear' | 'log';

/**
 * Density visualization mode.
 *
 * - 'total': Total electron density rho(r)
 * - 'difference': Difference density (D - D_atomic), planned for US-063
 */
export type DensityMode = 'total' | 'difference';

/**
 * Result of extracting a 2D slice from a 3D density grid.
 */
export interface DensitySliceResult {
  /** 2D array of density values [row][col] */
  data: number[][];
  /** Coordinate values along the horizontal axis of the slice */
  xCoords: number[];
  /** Coordinate values along the vertical axis of the slice */
  yCoords: number[];
  /** Axis label for the horizontal axis (e.g., "x (bohr)") */
  xLabel: string;
  /** Axis label for the vertical axis (e.g., "y (bohr)") */
  yLabel: string;
  /** Index along the perpendicular axis (nearest grid point) */
  perpendicularIndex: number;
  /** Actual position along the perpendicular axis (snapped to grid point, in Bohr) */
  actualPosition: number;
}

/**
 * Density visualization state managed by the SCF store.
 *
 * Controls the density isosurface and cross-section panels.
 * Only meaningful after SCF convergence.
 */
export interface DensityViewerState {
  /** Density mode: 'total' electron density or 'difference' (D - D_atomic) */
  densityMode: DensityMode;
  /** Isovalue for density isosurface in e/bohr^3 */
  densityIsovalue: number;
  /**
   * Isovalue for difference density isosurface in e/bohr^3.
   * Separate from densityIsovalue because difference density magnitudes
   * are much smaller than total density (typical: 0.005 vs 0.05).
   * @see US-063 Difference Density (AC3)
   */
  diffIsovalue: number;
  /** Plane axis for cross-section */
  planeAxis: PlaneAxis;
  /** Position along the perpendicular axis (Bohr) */
  planePosition: number;
  /** Color scale for cross-section */
  colorScale: ColorScale;
  /** Whether density visualization is active */
  showDensity: boolean;
}

/**
 * Default density viewer state values.
 */
export const DEFAULT_DENSITY_STATE: DensityViewerState = {
  densityMode: 'total',
  densityIsovalue: 0.05,
  diffIsovalue: 0.005,
  planeAxis: 'xz',
  planePosition: 0.0,
  colorScale: 'linear',
  showDensity: false,
};
