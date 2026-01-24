/**
 * ScfEnergyPlot - Energy vs iteration plot for SCF computation.
 *
 * Shows the total electronic energy at each SCF iteration,
 * visualizing the convergence trajectory. Updates live during
 * computation.
 *
 * @module components/scf/ScfEnergyPlot
 */

import { useMemo } from 'react';
import { PlotPanel } from '../common/PlotPanel';
import { useScfStore } from '../../stores/scfStore';
import type { Data, Layout } from 'plotly.js';

/**
 * Props for ScfEnergyPlot.
 */
interface ScfEnergyPlotProps {
  /** Minimum height in pixels (default: 280) */
  minHeight?: number;
}

/**
 * SCF energy convergence plot.
 *
 * Features:
 * - Energy (Hartree) vs iteration number
 * - Live updates during computation
 * - Horizontal annotation at final converged energy
 * - Clean scientific styling
 *
 * @example
 * ```tsx
 * <ScfEnergyPlot minHeight={300} />
 * ```
 */
export function ScfEnergyPlot({ minHeight = 280 }: ScfEnergyPlotProps) {
  const history = useScfStore((state) => state.history);
  const compute = useScfStore((state) => state.compute);

  const isRunning = compute.status === 'running';
  const hasData = history.length > 0;

  // Prepare plot data
  const plotData = useMemo((): Data[] => {
    if (!hasData) return [];

    const iterations = history.map((h) => h.iteration);
    const energies = history.map((h) => h.energy);

    return [
      {
        x: iterations,
        y: energies,
        type: 'scatter',
        mode: 'lines+markers',
        name: 'Energy',
        line: {
          color: '#3b82f6', // blue-500
          width: 2,
        },
        marker: {
          color: '#3b82f6',
          size: 6,
        },
        hovertemplate: 'Iteration %{x}<br>Energy: %{y:.10f} Ha<extra></extra>',
      },
    ];
  }, [history, hasData]);

  // Prepare layout with annotations
  const layout = useMemo((): Partial<Layout> => {
    const baseLayout: Partial<Layout> = {
      title: {
        text: 'Energy vs Iteration',
        font: { size: 14 },
      },
      xaxis: {
        title: { text: 'Iteration', font: { size: 12 } },
        dtick: 1,
        tickmode: history.length > 20 ? 'auto' : 'linear',
      },
      yaxis: {
        title: { text: 'Energy (Hartree)', font: { size: 12 } },
        exponentformat: 'e',
      },
      margin: { t: 50, r: 20, b: 50, l: 80 },
      showlegend: false,
    };

    // Add horizontal line at final energy if converged
    if (compute.status === 'success' && compute.result.converged && history.length > 0) {
      const finalEnergy = compute.result.energy;
      baseLayout.shapes = [
        {
          type: 'line',
          x0: 0,
          x1: history.length + 0.5,
          y0: finalEnergy,
          y1: finalEnergy,
          line: {
            color: '#22c55e', // green-500
            width: 1,
            dash: 'dot',
          },
        },
      ];
      baseLayout.annotations = [
        {
          x: history.length,
          y: finalEnergy,
          xanchor: 'left',
          yanchor: 'bottom',
          text: ` ${finalEnergy.toFixed(8)} Ha`,
          showarrow: false,
          font: { size: 10, color: '#22c55e' },
        },
      ];
    }

    return baseLayout;
  }, [history, compute]);

  // Empty state
  if (!hasData && !isRunning) {
    return (
      <div
        className="bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center"
        style={{ minHeight: `${minHeight}px` }}
      >
        <div className="text-center text-slate-500">
          <p className="text-sm">Energy plot will appear here</p>
          <p className="text-xs mt-1">Run an SCF calculation to see results</p>
        </div>
      </div>
    );
  }

  return (
    <PlotPanel
      data={plotData}
      layout={layout}
      loading={isRunning && history.length === 0}
      minHeight={minHeight}
      ariaLabel="SCF energy convergence plot"
    />
  );
}

export default ScfEnergyPlot;
