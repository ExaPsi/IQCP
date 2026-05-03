/**
 * DensityCrossSection - 2D heatmap of electron density cross-section.
 *
 * Renders a Plotly heatmap showing a planar slice through the 3D density
 * grid. The slice is extracted client-side from the cached grid data using
 * extractDensitySlice() -- no worker call needed, as slice extraction is
 * <1 ms for a 50x50 grid.
 *
 * Supports linear and log color scales. The log scale uses log10(rho + eps)
 * to handle zero values, revealing the exponential decay of electron density
 * away from nuclei -- a useful pedagogical view connecting to Gaussian
 * basis function tails.
 *
 * For difference density mode (US-063): uses a diverging RdBu colorscale
 * centered at zero, with symmetric range so that accumulation (positive)
 * and depletion (negative) are visually balanced.
 *
 * @module components/scf/DensityCrossSection
 * @see US-062 Density Isosurface & Cross-Sections
 * @see US-063 Difference Density
 */

import { useMemo } from 'react';
import { PlotPanel } from '../common/PlotPanel';
import { extractDensitySlice } from '../../lib/densitySlice';
import type { PlaneAxis, ColorScale } from '../../types/density';
import type { Data, Layout } from 'plotly.js';

// ============================================================================
// Types
// ============================================================================

/**
 * Props for the DensityCrossSection component.
 */
export interface DensityCrossSectionProps {
  /** Cached density grid values (flat array, C-order) */
  gridValues: number[] | null;
  /** Grid dimensions [nx, ny, nz] */
  gridDims: [number, number, number] | null;
  /** Grid origin [x, y, z] in Bohr */
  gridOrigin: [number, number, number] | null;
  /** Grid spacing in Bohr */
  gridSpacing: number | null;
  /** Which plane to slice */
  planeAxis: PlaneAxis;
  /** Position along the perpendicular axis (in Bohr) */
  planePosition: number;
  /** Color scale mode */
  colorScale: ColorScale;
  /** Whether data is loading */
  loading?: boolean;
  /**
   * Whether to use a diverging colorscale (for difference density).
   * When true, uses RdBu reversed colorscale centered at 0 with
   * symmetric range (red=depletion, blue=accumulation).
   * @see US-063 Difference Density
   */
  diverging?: boolean;
}

// ============================================================================
// Constants
// ============================================================================

/**
 * Small epsilon added to density before log10 to handle zeros.
 * log10(1e-15) = -15, which is below any physical density value.
 */
const LOG_EPSILON = 1e-15;

/**
 * Perpendicular axis label for each plane type.
 */
const PERP_AXIS_LABEL: Record<PlaneAxis, string> = {
  xy: 'z',
  xz: 'y',
  yz: 'x',
};

// ============================================================================
// Component
// ============================================================================

/**
 * DensityCrossSection - 2D planar cross-section of electron density.
 *
 * Displays a heatmap of density values on a selected plane (XY, XZ, or YZ)
 * at a specified position. Supports linear and logarithmic color scales.
 *
 * For difference density (diverging=true), uses a blue-white-red colorscale
 * centered at zero, showing accumulation (blue) and depletion (red) regions.
 *
 * The cursor hover shows coordinate position and density value for
 * quantitative inspection.
 *
 * @example
 * ```tsx
 * <DensityCrossSection
 *   gridValues={gridData.values}
 *   gridDims={gridData.gridDims}
 *   gridOrigin={gridData.gridOrigin}
 *   gridSpacing={gridData.gridSpacing}
 *   planeAxis="xz"
 *   planePosition={0.0}
 *   colorScale="linear"
 *   diverging={false}
 * />
 * ```
 */
export function DensityCrossSection({
  gridValues,
  gridDims,
  gridOrigin,
  gridSpacing,
  planeAxis,
  planePosition,
  colorScale,
  loading = false,
  diverging = false,
}: DensityCrossSectionProps) {
  // Extract the 2D slice from the cached 3D grid
  const sliceResult = useMemo(() => {
    if (!gridValues || !gridDims || !gridOrigin || gridSpacing === null) {
      return null;
    }

    return extractDensitySlice(
      gridValues,
      gridDims,
      gridOrigin,
      gridSpacing,
      planeAxis,
      planePosition
    );
  }, [gridValues, gridDims, gridOrigin, gridSpacing, planeAxis, planePosition]);

  // Build Plotly data and layout
  const { plotData, plotLayout } = useMemo(() => {
    if (!sliceResult) {
      return { plotData: [] as Data[], plotLayout: {} as Partial<Layout> };
    }

    const perpLabel = PERP_AXIS_LABEL[planeAxis];
    const isDiverging = diverging;

    // Choose colorscale and data transform based on mode
    let zData: number[][];
    let plotColorscale: string;
    let zmin: number | undefined;
    let zmax: number | undefined;
    let colorbarTitle: string;
    let hoverTemplate: string;

    if (isDiverging) {
      // Diverging mode for difference density (US-063):
      // Symmetric range centered at 0, RdBu reversed (blue=positive, red=negative)
      if (colorScale === 'log') {
        // For log scale with diverging data, use signed log transform:
        // sign(v) * log10(|v| + eps)
        zData = sliceResult.data.map((row) =>
          row.map((v) => {
            const sign = v >= 0 ? 1 : -1;
            return sign * Math.log10(Math.abs(v) + LOG_EPSILON);
          })
        );
        colorbarTitle = 'sign * log\u2081\u2080(|\u0394\u03C1|)';
      } else {
        zData = sliceResult.data;
        colorbarTitle = '\u0394\u03C1 (e/bohr\u00B3)';
      }

      // Compute symmetric range for diverging colorscale
      let absMax = 0;
      for (const row of zData) {
        for (const v of row) {
          const abs = Math.abs(v);
          if (abs > absMax) absMax = abs;
        }
      }
      zmin = -absMax;
      zmax = absMax;

      // RdBu reversed: blue for positive (accumulation), red for negative (depletion)
      plotColorscale = 'RdBu_r';

      hoverTemplate = `${sliceResult.xLabel.split(' ')[0]}=%{x:.2f}, ${sliceResult.yLabel.split(' ')[0]}=%{y:.2f}<br>\u0394\u03C1=%{z:.4e} e/bohr\u00B3<extra></extra>`;
    } else {
      // Standard mode for total density
      if (colorScale === 'log') {
        zData = sliceResult.data.map((row) =>
          row.map((v) => Math.log10(Math.max(v, LOG_EPSILON)))
        );
        colorbarTitle = 'log\u2081\u2080(\u03C1)';
        hoverTemplate = `${sliceResult.xLabel.split(' ')[0]}=%{x:.2f}, ${sliceResult.yLabel.split(' ')[0]}=%{y:.2f}<br>log\u2081\u2080(\u03C1)=%{z:.2f}<extra></extra>`;
      } else {
        zData = sliceResult.data;
        colorbarTitle = '\u03C1 (e/bohr\u00B3)';
        hoverTemplate = `${sliceResult.xLabel.split(' ')[0]}=%{x:.2f}, ${sliceResult.yLabel.split(' ')[0]}=%{y:.2f}<br>\u03C1=%{z:.4e} e/bohr\u00B3<extra></extra>`;
      }

      plotColorscale = 'Viridis';
    }

    const data: Data[] = [
      {
        type: 'heatmap',
        x: sliceResult.xCoords,
        y: sliceResult.yCoords,
        z: zData,
        colorscale: plotColorscale,
        ...(zmin !== undefined ? { zmin } : {}),
        ...(zmax !== undefined ? { zmax } : {}),
        colorbar: {
          title: {
            text: colorbarTitle,
            side: 'right',
          },
          tickformat: colorScale === 'log' ? '.1f' : '.2e',
        },
        hovertemplate: hoverTemplate,
        zsmooth: 'best',
      },
    ];

    const titlePrefix = isDiverging ? 'Difference density' : 'Density';
    const layout: Partial<Layout> = {
      title: {
        text: `${titlePrefix} cross-section (${planeAxis.toUpperCase()} plane, ${perpLabel} = ${sliceResult.actualPosition.toFixed(2)} bohr)`,
        font: { size: 14 },
      },
      xaxis: {
        title: { text: sliceResult.xLabel },
        scaleanchor: 'y',
        scaleratio: 1,
      },
      yaxis: {
        title: { text: sliceResult.yLabel },
      },
      margin: { t: 50, r: 80, b: 50, l: 60 },
    };

    return { plotData: data, plotLayout: layout };
  }, [sliceResult, colorScale, planeAxis, diverging]);

  if (!gridValues && !loading) {
    return (
      <div className="bg-slate-50 rounded-lg border border-slate-200 p-8 text-center">
        <p className="text-sm text-slate-500">
          Compute density grid to view cross-sections
        </p>
      </div>
    );
  }

  return (
    <PlotPanel
      data={plotData}
      layout={plotLayout}
      loading={loading}
      ariaLabel={`${diverging ? 'Difference density' : 'Electron density'} cross-section on ${planeAxis.toUpperCase()} plane`}
      className="w-full"
      minHeight={350}
    />
  );
}
