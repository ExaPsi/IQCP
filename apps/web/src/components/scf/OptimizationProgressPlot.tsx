/**
 * OptimizationProgressPlot - Live-updating dual-axis plot for geometry optimization.
 *
 * Shows energy (left y-axis) and max gradient (right y-axis, log scale) vs
 * optimization step number. Includes a horizontal dashed line at the gradient
 * convergence threshold (4.5e-4 Ha/bohr).
 *
 * @module components/scf/OptimizationProgressPlot
 * @see US-075 Optimization UI + Trajectory Animation (AC2)
 */

import React, { useMemo } from 'react';
import Plot from 'react-plotly.js';
import type { OptimizationStepResult } from '../../worker/protocol';

/**
 * Gradient convergence threshold in Ha/bohr.
 * Matches OptimizationConfig default in qc-core.
 */
const GRAD_THRESHOLD = 4.5e-4;

/**
 * Props for OptimizationProgressPlot.
 */
export interface OptimizationProgressPlotProps {
  /** Optimization steps (grows during optimization) */
  steps: OptimizationStepResult[];
  /** Whether optimization has converged (null = still running) */
  converged: boolean | null;
}

/**
 * Dual-axis plot showing energy and gradient convergence during optimization.
 *
 * Live-updates as new steps arrive. Uses Plotly with two y-axes:
 * - Left: total energy in Hartree
 * - Right: max gradient in Ha/bohr (log scale)
 */
export const OptimizationProgressPlot = React.memo(function OptimizationProgressPlot({
  steps,
  converged,
}: OptimizationProgressPlotProps) {
  const { stepNumbers, energies, gradients } = useMemo(() => {
    const nums: number[] = [];
    const ens: number[] = [];
    const grads: number[] = [];
    for (const s of steps) {
      nums.push(s.step);
      ens.push(s.energy);
      grads.push(s.maxGradient);
    }
    return { stepNumbers: nums, energies: ens, gradients: grads };
  }, [steps]);

  if (steps.length === 0) return null;

  const titleText = converged === true
    ? 'Optimization Converged'
    : converged === false && steps.length > 0
      ? 'Optimization (not converged)'
      : 'Optimization Progress';

  const data: Plotly.Data[] = [
    {
      x: stepNumbers,
      y: energies,
      type: 'scatter',
      mode: 'lines+markers',
      name: 'Energy (Ha)',
      line: { color: '#3b82f6', width: 2 },
      marker: { size: 5 },
      yaxis: 'y',
    },
    {
      x: stepNumbers,
      y: gradients,
      type: 'scatter',
      mode: 'lines+markers',
      name: 'Max |grad| (Ha/bohr)',
      line: { color: '#ef4444', width: 2, dash: 'dot' },
      marker: { size: 5 },
      yaxis: 'y2',
    },
    // Gradient threshold line
    {
      x: [stepNumbers[0], stepNumbers[stepNumbers.length - 1]],
      y: [GRAD_THRESHOLD, GRAD_THRESHOLD],
      type: 'scatter',
      mode: 'lines',
      name: `Threshold (${GRAD_THRESHOLD})`,
      line: { color: '#9ca3af', width: 1, dash: 'dash' },
      yaxis: 'y2',
      showlegend: true,
    },
  ];

  const layout: Partial<Plotly.Layout> = {
    title: { text: titleText, font: { size: 14 } },
    xaxis: {
      title: { text: 'Step' },
      dtick: 1,
    },
    yaxis: {
      title: { text: 'Energy (Ha)', font: { color: '#3b82f6' } },
      tickfont: { color: '#3b82f6' },
      side: 'left',
    },
    yaxis2: {
      title: { text: 'Max |gradient| (Ha/bohr)', font: { color: '#ef4444' } },
      tickfont: { color: '#ef4444' },
      type: 'log',
      overlaying: 'y',
      side: 'right',
    },
    legend: {
      x: 0.01,
      y: 0.99,
      bgcolor: 'rgba(255,255,255,0.8)',
      bordercolor: '#e2e8f0',
      borderwidth: 1,
      font: { size: 11 },
    },
    margin: { t: 40, r: 60, b: 50, l: 60 },
    autosize: true,
    height: 280,
  };

  const config: Partial<Plotly.Config> = {
    responsive: true,
    displayModeBar: 'hover' as const,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'] as Plotly.ModeBarDefaultButtons[],
  };

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      <Plot
        data={data}
        layout={layout}
        config={config}
        className="w-full"
        useResizeHandler
        style={{ width: '100%' }}
      />
    </div>
  );
});
