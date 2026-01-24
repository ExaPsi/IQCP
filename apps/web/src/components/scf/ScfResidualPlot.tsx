/**
 * ScfResidualPlot - Residual vs iteration plot for SCF computation.
 *
 * Shows the convergence residuals (delta E and DIIS error) on a
 * logarithmic scale, visualizing the convergence rate. Includes
 * a horizontal line at the convergence threshold.
 *
 * @module components/scf/ScfResidualPlot
 */

import { useMemo } from 'react';
import { PlotPanel } from '../common/PlotPanel';
import { useScfStore, CONVERGENCE_THRESHOLDS } from '../../stores/scfStore';
import type { Data, Layout } from 'plotly.js';

/**
 * Props for ScfResidualPlot.
 */
interface ScfResidualPlotProps {
  /** Minimum height in pixels (default: 280) */
  minHeight?: number;
}

/**
 * Get numeric value of convergence threshold.
 */
function getThresholdValue(thresholdStr: string): number {
  // Parse strings like "1e-6 Ha" or "1e-6"
  const match = thresholdStr.match(/(\d+e[+-]?\d+)/i);
  if (match) {
    return parseFloat(match[1]);
  }
  return 1e-6; // default fallback
}

/**
 * SCF residual convergence plot.
 *
 * Features:
 * - Log-scale Y-axis for residuals
 * - Dual traces: Delta E (red) and DIIS error (orange)
 * - Horizontal line at convergence threshold
 * - Live updates during computation
 *
 * @example
 * ```tsx
 * <ScfResidualPlot minHeight={300} />
 * ```
 */
export function ScfResidualPlot({ minHeight = 280 }: ScfResidualPlotProps) {
  const history = useScfStore((state) => state.history);
  const compute = useScfStore((state) => state.compute);
  const convergenceProfile = useScfStore((state) => state.convergenceProfile);

  const isRunning = compute.status === 'running';
  // Need at least 2 iterations to show residuals (first has no delta)
  const hasData = history.length > 1;

  // Get convergence threshold for the horizontal line
  const threshold = getThresholdValue(CONVERGENCE_THRESHOLDS[convergenceProfile].energy);

  // Prepare plot data
  const plotData = useMemo((): Data[] => {
    if (!hasData) return [];

    // Skip first iteration (no delta)
    const iterations = history.slice(1).map((h) => h.iteration);
    const deltas = history.slice(1).map((h) => Math.abs(h.delta));
    const diisErrors = history
      .slice(1)
      .map((h) => (h.diisError !== undefined ? h.diisError : null));

    const traces: Data[] = [
      {
        x: iterations,
        y: deltas,
        type: 'scatter',
        mode: 'lines+markers',
        name: '|Delta E|',
        line: {
          color: '#ef4444', // red-500
          width: 2,
        },
        marker: {
          color: '#ef4444',
          size: 6,
        },
        hovertemplate: 'Iteration %{x}<br>|Delta E|: %{y:.2e}<extra></extra>',
      },
    ];

    // Add DIIS error trace if available
    const hasDiisData = diisErrors.some((e) => e !== null);
    if (hasDiisData) {
      traces.push({
        x: iterations,
        y: diisErrors,
        type: 'scatter',
        mode: 'lines+markers',
        name: 'DIIS Error',
        line: {
          color: '#f97316', // orange-500
          width: 2,
        },
        marker: {
          color: '#f97316',
          size: 6,
        },
        connectgaps: true,
        hovertemplate: 'Iteration %{x}<br>DIIS Error: %{y:.2e}<extra></extra>',
      });
    }

    return traces;
  }, [history, hasData]);

  // Prepare layout
  const layout = useMemo((): Partial<Layout> => {
    // Calculate Y-axis range based on data
    let yMin = threshold / 100; // Show at least 2 orders below threshold
    let yMax = 1; // Start at 1 (typical initial residual scale)

    if (hasData) {
      const allValues = history
        .slice(1)
        .flatMap((h) => [Math.abs(h.delta), h.diisError])
        .filter((v): v is number => v !== undefined && v !== null && v > 0);

      if (allValues.length > 0) {
        const minVal = Math.min(...allValues);
        const maxVal = Math.max(...allValues);
        yMin = Math.min(yMin, minVal / 10);
        yMax = Math.max(yMax, maxVal * 10);
      }
    }

    const baseLayout: Partial<Layout> = {
      title: {
        text: 'Residual vs Iteration',
        font: { size: 14 },
      },
      xaxis: {
        title: { text: 'Iteration', font: { size: 12 } },
        dtick: 1,
        tickmode: history.length > 20 ? 'auto' : 'linear',
      },
      yaxis: {
        title: { text: 'Residual', font: { size: 12 } },
        type: 'log',
        range: [Math.log10(yMin), Math.log10(yMax)],
        exponentformat: 'e',
      },
      margin: { t: 50, r: 20, b: 50, l: 80 },
      showlegend: true,
      legend: {
        x: 1,
        xanchor: 'right',
        y: 1,
        bgcolor: 'rgba(255,255,255,0.8)',
        font: { size: 10 },
      },
      // Horizontal line at convergence threshold
      shapes: [
        {
          type: 'line',
          x0: 0,
          x1: Math.max(history.length + 1, 10),
          y0: threshold,
          y1: threshold,
          line: {
            color: '#22c55e', // green-500
            width: 2,
            dash: 'dash',
          },
        },
      ],
      annotations: [
        {
          x: 1,
          y: Math.log10(threshold),
          xref: 'paper',
          xanchor: 'right',
          yanchor: 'bottom',
          text: `Threshold (${threshold.toExponential(0)})`,
          showarrow: false,
          font: { size: 9, color: '#22c55e' },
        },
      ],
    };

    return baseLayout;
  }, [history, hasData, threshold]);

  // Empty state
  if (!hasData && !isRunning) {
    return (
      <div
        className="bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center"
        style={{ minHeight: `${minHeight}px` }}
      >
        <div className="text-center text-slate-500">
          <p className="text-sm">Residual plot will appear here</p>
          <p className="text-xs mt-1">Shows convergence rate on log scale</p>
        </div>
      </div>
    );
  }

  // Show minimal loading if running but no data yet
  if (isRunning && !hasData) {
    return (
      <div
        className="bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center"
        style={{ minHeight: `${minHeight}px` }}
      >
        <div className="flex items-center gap-3 text-slate-500">
          <span className="animate-spin rounded-full h-5 w-5 border-b-2 border-blue-600" />
          <span className="text-sm">Waiting for iterations...</span>
        </div>
      </div>
    );
  }

  return (
    <PlotPanel
      data={plotData}
      layout={layout}
      loading={false}
      minHeight={minHeight}
      ariaLabel="SCF residual convergence plot"
    />
  );
}

export default ScfResidualPlot;
