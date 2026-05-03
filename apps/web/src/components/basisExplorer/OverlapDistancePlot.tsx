/**
 * OverlapDistancePlot - Interactive overlap integral vs. distance visualization.
 *
 * Displays one or more overlap curves showing how the overlap integral S_ab
 * between two basis functions decays with interatomic distance R. Supports
 * overlaying up to 4 traces for comparison (e.g., STO-3G vs 3-21G,
 * or H-H vs H-He pairs).
 *
 * Features:
 * - Plotly scatter/line plot with R (bohr) x-axis, S_ab y-axis
 * - Linear/log y-axis toggle (log shows exponential decay as a line)
 * - Hover info with R and S_ab values
 * - Up to 4 overlaid traces with distinct colors and legend
 * - Educational info panel connecting to Module E overlap matrix
 *
 * @module components/basisExplorer/OverlapDistancePlot
 * @see US-054 Overlap vs. Distance Plot
 */

import { useMemo, useCallback, memo } from 'react';
import Plotly from 'plotly.js-cartesian-dist';
import { PlotPanel } from '../common/PlotPanel';
import type { Data, Layout } from 'plotly.js';
import type { OverlapDistanceResult } from '../../hooks/useOverlapDistance';

// ============================================================================
// Types
// ============================================================================

/**
 * Props for OverlapDistancePlot.
 */
export interface OverlapDistancePlotProps {
  /** Array of overlap distance results (one per trace, max 4) */
  traces: OverlapDistanceResult[];
  /** Whether data is currently being loaded */
  isLoading: boolean;
  /** Y-axis scale mode */
  yAxisScale: 'linear' | 'log';
  /** Callback when scale mode changes */
  onScaleChange: (scale: 'linear' | 'log') => void;
  /** Minimum plot height in pixels (default: 350) */
  minHeight?: number;
}

// ============================================================================
// Constants
// ============================================================================

/**
 * Color palette for overlap traces.
 * Four visually distinct, colorblind-accessible colors.
 */
const TRACE_COLORS = [
  '#2563eb', // blue-600
  '#dc2626', // red-600
  '#16a34a', // green-600
  '#9333ea', // purple-600
];

/**
 * Build a display label for a trace from its result metadata.
 * Format: "H 1s / H 1s (STO-3G)" or "H 1s / He 1s (STO-3G / STO-3G)"
 */
function traceLabel(result: OverlapDistanceResult): string {
  const pairLabel = `${result.shellLabelA} / ${result.shellLabelB}`;
  if (result.basisA === result.basisB) {
    return `${pairLabel} (${result.basisA.toUpperCase()})`;
  }
  return `${pairLabel} (${result.basisA.toUpperCase()} / ${result.basisB.toUpperCase()})`;
}

// ============================================================================
// Scale Toggle Component
// ============================================================================

/**
 * Toggle button for linear/log y-axis scale.
 * Replicates the pattern from RadialProfilePlot.
 */
function ScaleToggle({
  scale,
  onChange,
}: {
  scale: 'linear' | 'log';
  onChange: (scale: 'linear' | 'log') => void;
}) {
  return (
    <div className="flex items-center gap-2">
      <span className="text-xs text-slate-500">Y-axis:</span>
      <div className="inline-flex rounded-md border border-slate-300" role="group">
        <button
          type="button"
          onClick={() => onChange('linear')}
          aria-pressed={scale === 'linear'}
          className={`px-3 py-1 text-xs font-medium rounded-l-md transition-colors
            focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
            ${
              scale === 'linear'
                ? 'bg-primary-100 text-primary-800 border-primary-300'
                : 'bg-white text-slate-600 hover:bg-slate-50'
            }`}
        >
          Linear
        </button>
        <button
          type="button"
          onClick={() => onChange('log')}
          aria-pressed={scale === 'log'}
          className={`px-3 py-1 text-xs font-medium rounded-r-md border-l border-slate-300 transition-colors
            focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
            ${
              scale === 'log'
                ? 'bg-primary-100 text-primary-800 border-primary-300'
                : 'bg-white text-slate-600 hover:bg-slate-50'
            }`}
        >
          Log
        </button>
      </div>
    </div>
  );
}

// ============================================================================
// Main Component
// ============================================================================

/**
 * Interactive overlap vs. distance plot.
 *
 * Each trace is a line plot of S_ab(R) for a pair of basis functions.
 * Up to 4 traces can be overlaid for comparison with distinct colors.
 *
 * @example
 * ```tsx
 * <OverlapDistancePlot
 *   traces={[overlapResult1, overlapResult2]}
 *   isLoading={false}
 *   yAxisScale="linear"
 *   onScaleChange={setScale}
 * />
 * ```
 */
export const OverlapDistancePlot = memo(function OverlapDistancePlot({
  traces,
  isLoading,
  yAxisScale,
  onScaleChange,
  minHeight = 350,
}: OverlapDistancePlotProps) {
  // Stable divId for Plotly SVG export
  const plotDivId = 'overlap-distance-plot';

  // SVG export handler
  const handleExportSvg = useCallback(() => {
    const plotEl = document.getElementById(plotDivId);
    if (plotEl) {
      Plotly.downloadImage(plotEl as unknown as Plotly.PlotlyHTMLElement, {
        format: 'svg',
        filename: `overlap-distance-${traces.length}traces`,
        width: 800,
        height: 500,
      });
    }
  }, [traces.length]);

  // Build Plotly data traces
  const plotData = useMemo((): Data[] => {
    if (traces.length === 0) return [];

    return traces.map((result, idx) => {
      const color = TRACE_COLORS[idx % TRACE_COLORS.length];
      const label = traceLabel(result);

      // For log scale, filter out zero/negative values
      let rValues = result.rValues;
      let overlapValues = result.overlapValues;

      if (yAxisScale === 'log') {
        const filtered = filterPositive(rValues, overlapValues);
        rValues = filtered.x;
        overlapValues = filtered.y;
      }

      return {
        x: rValues,
        y: overlapValues,
        type: 'scatter' as const,
        mode: 'lines' as const,
        name: label,
        line: {
          color,
          width: 2.5,
        },
        hovertemplate:
          `<b>${label}</b><br>` +
          'R = %{x:.3f} bohr<br>' +
          'S<sub>ab</sub> = %{y:.6e}' +
          '<extra></extra>',
      };
    });
  }, [traces, yAxisScale]);

  // Build Plotly layout
  const layout = useMemo((): Partial<Layout> => {
    return {
      title: {
        text: 'Overlap S<sub>ab</sub> vs. Distance',
        font: { size: 14 },
      },
      xaxis: {
        title: { text: 'Interatomic Distance R (bohr)', font: { size: 12 } },
        rangemode: 'tozero',
        zeroline: true,
      },
      yaxis: {
        title: { text: 'Overlap Integral S<sub>ab</sub>', font: { size: 12 } },
        type: yAxisScale === 'log' ? 'log' : 'linear',
        exponentformat: 'e',
        zeroline: yAxisScale === 'linear',
      },
      margin: { t: 50, r: 20, b: 50, l: 70 },
      showlegend: traces.length > 1,
      legend: {
        x: 1,
        y: 1,
        xanchor: 'right',
        yanchor: 'top',
        font: { size: 10 },
        bgcolor: 'rgba(255,255,255,0.8)',
        borderwidth: 1,
        bordercolor: '#e2e8f0',
      },
      hovermode: 'closest',
    };
  }, [yAxisScale, traces.length]);

  // Empty state: no traces and not loading
  if (traces.length === 0 && !isLoading) {
    return (
      <div
        className="bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center"
        style={{ minHeight: `${minHeight}px` }}
        role="img"
        aria-label="Overlap vs distance plot placeholder"
      >
        <div className="text-center text-slate-500 px-6">
          <div className="inline-flex items-center justify-center w-12 h-12 rounded-full bg-slate-100 mb-3">
            <svg
              className="w-6 h-6 text-slate-400"
              fill="none"
              viewBox="0 0 24 24"
              strokeWidth={1.5}
              stroke="currentColor"
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                d="M7.5 14.25v2.25m3-4.5v4.5m3-6.75v6.75m3-9v9M6 20.25h12A2.25 2.25 0 0020.25 18V6A2.25 2.25 0 0018 3.75H6A2.25 2.25 0 003.75 6v12A2.25 2.25 0 006 20.25z"
              />
            </svg>
          </div>
          <p className="text-sm font-medium">Overlap vs. Distance</p>
          <p className="text-xs mt-1">
            Select shells to view how the overlap integral decays with distance
          </p>
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-2">
      {/* Header with scale toggle */}
      <div className="flex items-center justify-between">
        <h3 className="text-sm font-semibold text-slate-700">
          Overlap vs. Distance
        </h3>
        <ScaleToggle scale={yAxisScale} onChange={onScaleChange} />
      </div>

      {/* Plot */}
      <PlotPanel
        data={plotData}
        layout={layout}
        loading={isLoading}
        minHeight={minHeight}
        divId={plotDivId}
        ariaLabel={
          traces.length > 0
            ? `Overlap integral vs interatomic distance plot with ${traces.length} trace${traces.length > 1 ? 's' : ''}`
            : 'Overlap vs distance plot'
        }
      />

      {/* Export SVG button */}
      {traces.length > 0 && (
        <div className="flex justify-end">
          <button
            type="button"
            onClick={handleExportSvg}
            className="inline-flex items-center gap-1.5 px-2.5 py-1 text-xs font-medium text-slate-500 hover:text-slate-700 transition-colors
              focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500 rounded"
            aria-label="Export plot as SVG"
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" />
            </svg>
            Export SVG
          </button>
        </div>
      )}

      {/* Educational cross-reference to Module E */}
      {traces.length > 0 && (
        <div className="bg-amber-50 rounded-lg px-4 py-3 border border-amber-100">
          <div className="flex items-start gap-2">
            <span
              className="inline-flex items-center justify-center w-5 h-5 rounded-full bg-amber-200 text-amber-800 text-xs font-bold flex-shrink-0 mt-0.5"
              aria-hidden="true"
            >
              E
            </span>
            <p className="text-xs text-amber-800">
              <strong>Connection to Module E (SCF):</strong> This overlap integral S<sub>ab</sub> is
              the same quantity that appears as the off-diagonal element of the overlap
              matrix <strong>S</strong> in the SCF procedure. When S<sub>ab</sub> is large,
              the two basis functions have significant spatial overlap, indicating potential
              for covalent bonding. The overlap matrix measures the non-orthogonality of
              the basis set.
            </p>
          </div>
        </div>
      )}

      {/* Log scale insight */}
      {traces.length > 0 && yAxisScale === 'log' && (
        <div className="bg-blue-50 rounded-lg px-4 py-3 border border-blue-100">
          <p className="text-xs text-blue-800">
            <strong>Log scale insight:</strong> On a log scale, the exponential decay of the
            overlap integral appears as an approximately straight line. The slope is governed
            by the smallest exponent in the shell pair: more diffuse functions (smaller exponents)
            decay more slowly, extending the overlap to larger distances.
          </p>
        </div>
      )}
    </div>
  );
});

// ============================================================================
// Utilities
// ============================================================================

/**
 * Filter out non-positive values for log-scale display.
 * Plotly cannot display zero or negative values on a log axis.
 */
function filterPositive(
  x: number[],
  y: number[],
): { x: number[]; y: number[] } {
  const filteredX: number[] = [];
  const filteredY: number[] = [];
  for (let i = 0; i < x.length; i++) {
    if (y[i] > 0) {
      filteredX.push(x[i]);
      filteredY.push(y[i]);
    }
  }
  return { x: filteredX, y: filteredY };
}

export default OverlapDistancePlot;
