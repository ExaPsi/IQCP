/**
 * BasisComparisonChart - Overlaid convergence plot for basis set comparison.
 *
 * Displays SCF energy vs. iteration with one color-coded trace per basis set.
 * Only renders when there are comparison results to show.
 *
 * @module components/scf/BasisComparisonChart
 * @see US-045 Basis Set Comparison Mode
 */

import { useMemo } from 'react';
import { PlotPanel } from '../common/PlotPanel';
import { useScfStore } from '../../stores/scfStore';
import type { Data, Layout } from 'plotly.js';

// ============================================================================
// Constants
// ============================================================================

/**
 * Color palette for up to 4 basis sets.
 *
 * Each color is visually distinct and accessible:
 * - Blue (STO-3G default): #3b82f6 (Tailwind blue-500)
 * - Green (3-21G default): #10b981 (Tailwind emerald-500)
 * - Amber (6-31G default): #f59e0b (Tailwind amber-500)
 * - Purple (6-31G* default): #8b5cf6 (Tailwind violet-500)
 */
const BASIS_COLORS: string[] = [
  '#3b82f6', // blue
  '#10b981', // green
  '#f59e0b', // amber
  '#8b5cf6', // purple
];

/**
 * Human-readable basis set label mapping.
 */
const BASIS_LABELS: Record<string, string> = {
  'sto-3g': 'STO-3G',
  '3-21g': '3-21G',
  '6-31g': '6-31G',
  '6-31g*': '6-31G*',
  '6-31+g*': '6-31+G*',
  'cc-pvdz': 'cc-pVDZ',
};

// ============================================================================
// Component
// ============================================================================

/**
 * Props for BasisComparisonChart.
 */
interface BasisComparisonChartProps {
  /** Minimum height of the plot in pixels (default: 320) */
  minHeight?: number;
}

/**
 * Overlaid convergence plot comparing SCF energy trajectories
 * across multiple basis sets.
 *
 * Features:
 * - One trace per basis set, color-coded
 * - Energy (Ha) on y-axis, iteration number on x-axis
 * - Hover info showing basis set name, iteration, and energy
 * - Only rendered when comparison results exist
 *
 * @example
 * ```tsx
 * <BasisComparisonChart minHeight={350} />
 * ```
 */
export function BasisComparisonChart({ minHeight = 320 }: BasisComparisonChartProps) {
  const results = useScfStore((state) => state.basisCompareState.results);
  const comparing = useScfStore((state) => state.basisCompareState.comparing);

  const hasData = results.length > 0;

  // Build Plotly traces - one per basis set
  const plotData = useMemo((): Data[] => {
    if (!hasData) return [];

    return results
      .filter((r) => r.history.length > 0 && !r.error)
      .map((result, index) => {
        const color = BASIS_COLORS[index % BASIS_COLORS.length];
        const label = BASIS_LABELS[result.basisName] ?? result.basisName.toUpperCase();

        return {
          x: result.history.map((h) => h.iteration),
          y: result.history.map((h) => h.energy),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: label,
          line: {
            color,
            width: 2,
          },
          marker: {
            color,
            size: 5,
            symbol: 'circle',
          },
          hovertemplate:
            `<b>${label}</b><br>` +
            'Iteration %{x}<br>' +
            'Energy: %{y:.10f} Ha' +
            '<extra></extra>',
        };
      });
  }, [results, hasData]);

  // Build layout
  const layout = useMemo((): Partial<Layout> => ({
    title: {
      text: 'SCF Convergence by Basis Set',
      font: { size: 14 },
    },
    xaxis: {
      title: { text: 'Iteration', font: { size: 12 } },
      dtick: 1,
    },
    yaxis: {
      title: { text: 'Energy (Hartree)', font: { size: 12 } },
      exponentformat: 'e',
    },
    margin: { t: 50, r: 30, b: 50, l: 80 },
    showlegend: true,
    legend: {
      x: 1,
      y: 1,
      xanchor: 'right',
      yanchor: 'top',
      font: { size: 11 },
      bgcolor: 'rgba(255,255,255,0.85)',
      bordercolor: '#e2e8f0',
      borderwidth: 1,
    },
  }), []);

  // Don't render if no data and not comparing
  if (!hasData && !comparing) {
    return null;
  }

  return (
    <PlotPanel
      data={plotData}
      layout={layout}
      loading={comparing && results.length === 0}
      minHeight={minHeight}
      ariaLabel="Overlaid SCF convergence plot comparing energy trajectories across different basis sets"
    />
  );
}

export default BasisComparisonChart;
