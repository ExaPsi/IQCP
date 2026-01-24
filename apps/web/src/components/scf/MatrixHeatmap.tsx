/**
 * MatrixHeatmap - Plotly heatmap for displaying SCF matrices.
 *
 * Visualizes symmetric matrices (overlap, Hcore, Fock, density) as
 * interactive heatmaps with appropriate colorscales for each type.
 *
 * Features:
 * - Reshapes flat row-major arrays into 2D matrices
 * - Uses diverging colorscales (RdBu) for signed values
 * - Uses sequential colorscales (Blues/Greens) for positive values
 * - Interactive hover with cell values
 * - 1-indexed axis labels for student clarity
 *
 * @module components/scf/MatrixHeatmap
 */

import { useMemo } from 'react';
import { PlotPanel } from '../common/PlotPanel';
import type { Data, Layout } from 'plotly.js';

/**
 * Supported colorscale types for matrix visualization.
 *
 * - 'RdBu': Diverging red-white-blue, centered at 0 (for signed values like Hcore, Fock)
 * - 'Blues': Sequential blue gradient (for overlap matrix, 0-1 range)
 * - 'Greens': Sequential green gradient (for density matrix, 0-2 range)
 */
export type MatrixColorscale = 'RdBu' | 'Blues' | 'Greens';

/**
 * Props for MatrixHeatmap component.
 */
export interface MatrixHeatmapProps {
  /** Flat array of matrix values (row-major order) */
  data: number[];
  /** Matrix dimension (creates nbf x nbf matrix) */
  nbf: number;
  /** Title shown above the heatmap */
  title: string;
  /** Plotly colorscale (default: 'RdBu' for signed values) */
  colorscale?: MatrixColorscale;
  /** Additional CSS classes */
  className?: string;
  /** Minimum height in pixels (default: 320) */
  minHeight?: number;
}

/**
 * Reshape a flat row-major array into a 2D matrix.
 *
 * @param data - Flat array of length nbf * nbf
 * @param nbf - Number of basis functions (matrix dimension)
 * @returns 2D array where z[i][j] = data[i * nbf + j]
 */
function reshapeToMatrix(data: number[], nbf: number): number[][] {
  const matrix: number[][] = [];
  for (let i = 0; i < nbf; i++) {
    const row: number[] = [];
    for (let j = 0; j < nbf; j++) {
      row.push(data[i * nbf + j]);
    }
    matrix.push(row);
  }
  return matrix;
}

/**
 * Generate axis tick labels (1-indexed for student clarity).
 *
 * @param nbf - Number of basis functions
 * @returns Array of tick labels ["1", "2", ..., "nbf"]
 */
function generateTickLabels(nbf: number): string[] {
  return Array.from({ length: nbf }, (_, i) => String(i + 1));
}

/**
 * Get colorscale configuration based on type.
 *
 * For diverging scales (RdBu), we center at 0 by setting zmid.
 * For sequential scales, we don't set zmid.
 */
function getColorscaleConfig(
  colorscale: MatrixColorscale,
  data: number[]
): { zmid?: number; zmin?: number; zmax?: number } {
  if (colorscale === 'RdBu') {
    // Diverging scale: center at 0
    const absMax = Math.max(...data.map(Math.abs));
    return {
      zmid: 0,
      zmin: -absMax,
      zmax: absMax,
    };
  }

  // Sequential scales: use natural data range
  return {};
}

/**
 * MatrixHeatmap - Interactive heatmap visualization for SCF matrices.
 *
 * Displays symmetric matrices (S, Hcore, F, D) as color-coded heatmaps
 * with hover information showing cell values.
 *
 * @example
 * ```tsx
 * <MatrixHeatmap
 *   data={scfMatrices.sMatrix}
 *   nbf={scfMatrices.nbf}
 *   title="Overlap Matrix (S)"
 *   colorscale="Blues"
 * />
 * ```
 */
export function MatrixHeatmap({
  data,
  nbf,
  title,
  colorscale = 'RdBu',
  className = '',
  minHeight = 320,
}: MatrixHeatmapProps) {
  // Reshape flat array to 2D matrix
  const matrix = useMemo(() => reshapeToMatrix(data, nbf), [data, nbf]);

  // Generate tick labels (1-indexed)
  const tickLabels = useMemo(() => generateTickLabels(nbf), [nbf]);

  // Get colorscale configuration
  const colorConfig = useMemo(
    () => getColorscaleConfig(colorscale, data),
    [colorscale, data]
  );

  // Prepare Plotly data
  const plotData = useMemo((): Data[] => {
    return [
      {
        z: matrix,
        type: 'heatmap',
        colorscale: colorscale,
        x: tickLabels,
        y: tickLabels,
        // Reverse y-axis so (1,1) is in top-left corner
        // This matches standard matrix display convention
        ...colorConfig,
        hovertemplate:
          'Row: %{y}<br>Col: %{x}<br>Value: %{z:.6e}<extra></extra>',
        showscale: true,
        colorbar: {
          title: {
            text: 'Value',
            font: { size: 11 },
          },
          thickness: 15,
          len: 0.9,
          tickformat: '.2e',
          tickfont: { size: 10 },
        },
      },
    ];
  }, [matrix, tickLabels, colorscale, colorConfig]);

  // Prepare layout
  const layout = useMemo((): Partial<Layout> => {
    return {
      title: {
        text: title,
        font: { size: 14 },
      },
      xaxis: {
        title: { text: 'Column (basis function)', font: { size: 11 } },
        tickmode: 'array',
        tickvals: tickLabels,
        ticktext: tickLabels,
        tickfont: { size: 10 },
        side: 'bottom',
        constrain: 'domain',
      },
      yaxis: {
        title: { text: 'Row (basis function)', font: { size: 11 } },
        tickmode: 'array',
        tickvals: tickLabels,
        ticktext: tickLabels,
        tickfont: { size: 10 },
        autorange: 'reversed', // Top-to-bottom for matrix convention
        scaleanchor: 'x', // Keep aspect ratio square
        constrain: 'domain',
      },
      margin: { t: 50, r: 80, b: 60, l: 70 },
    };
  }, [title, tickLabels]);

  // Handle empty or invalid data
  if (data.length === 0 || nbf === 0) {
    return (
      <div
        className={`bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center ${className}`}
        style={{ minHeight: `${minHeight}px` }}
      >
        <div className="text-center text-slate-500">
          <p className="text-sm">No matrix data available</p>
        </div>
      </div>
    );
  }

  // Validate data length
  if (data.length !== nbf * nbf) {
    return (
      <div
        className={`bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center ${className}`}
        style={{ minHeight: `${minHeight}px` }}
      >
        <div className="text-center text-red-500">
          <p className="text-sm">Invalid matrix data</p>
          <p className="text-xs mt-1">
            Expected {nbf * nbf} elements, got {data.length}
          </p>
        </div>
      </div>
    );
  }

  return (
    <PlotPanel
      data={plotData}
      layout={layout}
      minHeight={minHeight}
      className={className}
      ariaLabel={`${title} heatmap visualization`}
    />
  );
}

export default MatrixHeatmap;
