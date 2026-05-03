/**
 * IntegralMatrixHeatmap - Interactive Plotly heatmap for one-electron integral matrices.
 *
 * Displays a one-electron integral matrix (S, T, V, or H^core) as a color-coded
 * heatmap with basis function labels on both axes, cell click selection, and
 * keyboard navigation for accessibility.
 *
 * Colormap selection follows the physical properties of each matrix:
 * - S (Overlap): Blues (all non-negative, diagonal = 1.0)
 * - T (Kinetic): Blues (all positive, kinetic energy is positive-definite)
 * - V (Nuclear): Reds reversed (all negative, attraction)
 * - H^core: RdBu diverging centered at 0 (mixed positive T and negative V)
 *
 * @module components/integrals/IntegralMatrixHeatmap
 * @see US-055 Integral Matrix Heatmap
 */

import { useMemo, useCallback, useRef, useEffect } from 'react';
import Plotly from 'plotly.js-cartesian-dist';
import createPlotlyComponent from 'react-plotly.js/factory';
import type { Data, Layout, Config, PlotMouseEvent } from 'plotly.js';
import type { MatrixTab, SelectedCell } from '../../stores/integralStore';

const Plot = createPlotlyComponent(Plotly);

// ============================================================================
// Types
// ============================================================================

export interface IntegralMatrixHeatmapProps {
  /** Flat row-major matrix data */
  matrix: number[];
  /** Number of basis functions (matrix is nbf x nbf) */
  nbf: number;
  /** Basis function labels for axes */
  labels: string[];
  /** Matrix type for colormap selection and title */
  matrixType: MatrixTab;
  /** Currently selected cell (null if none) */
  selectedCell: SelectedCell | null;
  /** Callback when cell is clicked or deselected (null on Escape) */
  onCellSelect: (cell: SelectedCell | null) => void;
  /** Minimum height in pixels */
  minHeight?: number;
}

// ============================================================================
// Matrix Type Configuration
// ============================================================================

interface MatrixConfig {
  title: string;
  colorscale: string;
  diverging: boolean;
}

const MATRIX_CONFIGS: Record<MatrixTab, MatrixConfig> = {
  S: {
    title: 'Overlap Matrix (S)',
    colorscale: 'Blues',
    diverging: false,
  },
  T: {
    title: 'Kinetic Energy Matrix (T)',
    colorscale: 'Blues',
    diverging: false,
  },
  V: {
    title: 'Nuclear Attraction Matrix (V)',
    colorscale: 'Reds',
    diverging: false,
  },
  Hcore: {
    title: 'Core Hamiltonian (H<sup>core</sup>)',
    colorscale: 'RdBu',
    diverging: true,
  },
};

// ============================================================================
// Helpers
// ============================================================================

/**
 * Reshape a flat row-major array into a 2D matrix.
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

// ============================================================================
// Component
// ============================================================================

/**
 * Interactive heatmap for one-electron integral matrices.
 *
 * Features:
 * - Color-coded cells with matrix-type-appropriate colormap
 * - Basis function labels on both axes
 * - Click-to-select cell interaction
 * - Hover template with full-precision values
 * - Keyboard arrow navigation for accessibility
 * - Selected cell highlighted with annotation overlay
 */
export function IntegralMatrixHeatmap({
  matrix,
  nbf,
  labels,
  matrixType,
  selectedCell,
  onCellSelect,
  minHeight = 400,
}: IntegralMatrixHeatmapProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const config = MATRIX_CONFIGS[matrixType];

  // Reshape flat array to 2D
  const matrix2D = useMemo(() => reshapeToMatrix(matrix, nbf), [matrix, nbf]);

  // Build colorscale-specific configuration
  const colorConfig = useMemo(() => {
    if (config.diverging) {
      const absMax = Math.max(...matrix.map(Math.abs));
      return { zmid: 0, zmin: -absMax, zmax: absMax };
    }
    // For V matrix (all negative), reverse the scale
    if (matrixType === 'V') {
      return { reversescale: true };
    }
    return {};
  }, [matrix, config.diverging, matrixType]);

  // Plotly data traces
  const plotData = useMemo((): Data[] => {
    const traces: Data[] = [
      {
        z: matrix2D,
        x: labels,
        y: labels,
        type: 'heatmap',
        colorscale: config.colorscale,
        ...colorConfig,
        hovertemplate:
          '%{y} | %{x}<br>Value: %{z:.8e}<extra></extra>',
        showscale: true,
        colorbar: {
          title: { text: 'Value', font: { size: 11 } },
          thickness: 15,
          len: 0.9,
          tickformat: '.3e',
          tickfont: { size: 10 },
        },
      } as Data,
    ];

    // Add selection highlight via a scatter trace
    if (selectedCell) {
      traces.push({
        x: [labels[selectedCell.col]],
        y: [labels[selectedCell.row]],
        type: 'scatter',
        mode: 'markers',
        marker: {
          size: 20,
          color: 'rgba(0,0,0,0)',
          line: {
            color: '#000000',
            width: 3,
          },
          symbol: 'square',
        },
        hoverinfo: 'skip',
        showlegend: false,
      } as Data);
    }

    return traces;
  }, [matrix2D, labels, config.colorscale, colorConfig, selectedCell]);

  // Plotly layout
  const layout = useMemo((): Partial<Layout> => {
    return {
      title: {
        text: config.title,
        font: { size: 14 },
      },
      xaxis: {
        title: { text: 'Basis Function', font: { size: 11 } },
        tickmode: 'array' as const,
        tickvals: labels,
        ticktext: labels,
        tickfont: { size: nbf > 10 ? 8 : 10 },
        tickangle: nbf > 7 ? -45 : 0,
        side: 'bottom' as const,
        constrain: 'domain' as const,
      },
      yaxis: {
        title: { text: 'Basis Function', font: { size: 11 } },
        tickmode: 'array' as const,
        tickvals: labels,
        ticktext: labels,
        tickfont: { size: nbf > 10 ? 8 : 10 },
        autorange: 'reversed' as const,
        scaleanchor: 'x',
        constrain: 'domain' as const,
      },
      margin: { t: 50, r: 90, b: nbf > 7 ? 80 : 60, l: nbf > 7 ? 80 : 60 },
      autosize: true,
      paper_bgcolor: 'transparent',
      plot_bgcolor: 'white',
      hovermode: 'closest' as const,
    };
  }, [config.title, labels, nbf]);

  // Plotly config
  const plotConfig = useMemo((): Partial<Config> => ({
    responsive: true,
    displayModeBar: 'hover' as const,
    modeBarButtonsToRemove: [
      'lasso2d',
      'select2d',
      'autoScale2d',
      'hoverClosestCartesian',
      'hoverCompareCartesian',
    ],
    displaylogo: false,
    toImageButtonOptions: {
      format: 'png' as const,
      filename: `iqcp-${matrixType}-matrix`,
      scale: 2,
    },
  }), [matrixType]);

  // Click handler: extract cell coordinates from Plotly event
  const handleClick = useCallback(
    (event: PlotMouseEvent) => {
      if (!event.points || event.points.length === 0) return;
      const point = event.points[0];

      // Map the tick label back to an index
      const colLabel = point.x as string;
      const rowLabel = point.y as string;
      const col = labels.indexOf(colLabel);
      const row = labels.indexOf(rowLabel);

      if (row >= 0 && col >= 0) {
        onCellSelect({ row, col });
      }
    },
    [labels, onCellSelect],
  );

  // Keyboard navigation for accessibility
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const handleKeyDown = (e: KeyboardEvent) => {
      if (!selectedCell) {
        // If no cell is selected, select (0, 0) on any arrow key
        if (['ArrowUp', 'ArrowDown', 'ArrowLeft', 'ArrowRight'].includes(e.key)) {
          e.preventDefault();
          onCellSelect({ row: 0, col: 0 });
        }
        return;
      }

      let { row, col } = selectedCell;

      switch (e.key) {
        case 'ArrowUp':
          e.preventDefault();
          row = Math.max(0, row - 1);
          break;
        case 'ArrowDown':
          e.preventDefault();
          row = Math.min(nbf - 1, row + 1);
          break;
        case 'ArrowLeft':
          e.preventDefault();
          col = Math.max(0, col - 1);
          break;
        case 'ArrowRight':
          e.preventDefault();
          col = Math.min(nbf - 1, col + 1);
          break;
        case 'Escape':
          e.preventDefault();
          onCellSelect(null);
          return;
        default:
          return;
      }

      onCellSelect({ row, col });
    };

    container.addEventListener('keydown', handleKeyDown);
    return () => container.removeEventListener('keydown', handleKeyDown);
  }, [selectedCell, nbf, onCellSelect]);

  // Validation
  if (matrix.length === 0 || nbf === 0) {
    return (
      <div
        className="bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center"
        style={{ minHeight: `${minHeight}px` }}
      >
        <p className="text-sm text-slate-500">No matrix data available</p>
      </div>
    );
  }

  if (matrix.length !== nbf * nbf) {
    return (
      <div
        className="bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center"
        style={{ minHeight: `${minHeight}px` }}
      >
        <div className="text-center text-red-500">
          <p className="text-sm">Invalid matrix data</p>
          <p className="text-xs mt-1">
            Expected {nbf * nbf} elements, got {matrix.length}
          </p>
        </div>
      </div>
    );
  }

  const plotDivId = `integral-${matrixType}-heatmap`;

  return (
    <div className="space-y-0">
      <div
        ref={containerRef}
        className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden"
        style={{ minHeight: `${minHeight}px` }}
        role="img"
        aria-label={`${config.title} heatmap with ${nbf} basis functions`}
        tabIndex={0}
      >
        <Plot
          divId={plotDivId}
          data={plotData}
          layout={layout}
          config={plotConfig}
          onClick={handleClick}
          useResizeHandler={true}
          style={{ width: '100%', height: '100%' }}
        />
      </div>
      <div className="flex justify-end mt-1">
        <button
          onClick={() => {
            const el = document.getElementById(plotDivId);
            if (el && (window as unknown as Record<string, unknown>).Plotly) {
              (
                (window as unknown as Record<string, unknown>).Plotly as {
                  downloadImage: (
                    el: HTMLElement,
                    opts: Record<string, unknown>,
                  ) => void;
                }
              ).downloadImage(el, {
                format: 'svg',
                filename: `iqcp-${matrixType}-matrix`,
                width: 600,
                height: 500,
              });
            }
          }}
          className="text-xs text-slate-500 hover:text-slate-700 flex items-center gap-1"
          title="Export heatmap as SVG"
        >
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
          </svg>
          Export SVG
        </button>
      </div>
    </div>
  );
}
