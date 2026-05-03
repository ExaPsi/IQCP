/**
 * BoysChart - Plotly chart for Boys function F_m(T) curves.
 *
 * Displays the computed sweep data as a line chart with:
 * - F_m(T) vs T curve
 * - Vertical marker at current T value
 * - Dynamic regime boundary based on m-dependent turnover point
 * - Optional log scale for y-axis
 *
 * The regime boundary varies with order m:
 * - m=0,1: turnover=0 (always recurrence)
 * - m=5: turnover~2.1
 * - m=10: turnover~4.05
 *
 * @module components/boys/BoysChart
 */

import { useMemo, useCallback } from 'react';
import Plotly from 'plotly.js-cartesian-dist';
import { PlotPanel } from '../common';
import type { Data, Layout, Shape, Annotations } from 'plotly.js';
import type { BoysSweepResult } from '../../worker/protocol';
import {
  getTurnoverPoint,
  getAsymptoticThreshold,
  getMaxTValue,
  THEORETICAL_SMALL_T_THRESHOLD,
} from '../../lib/boysConstants';

// Use Plotly's Annotations type for individual annotation objects
type PlotlyAnnotation = Annotations;

/**
 * Plotly div ID for SVG export.
 */
const SWEEP_PLOT_DIV_ID = 'boys-sweep-plot';

/**
 * Props for BoysChart component.
 */
export interface BoysChartProps {
  /** Sweep computation results */
  sweepData: BoysSweepResult | null;
  /** Use logarithmic y-axis */
  logScale: boolean;
  /** Current T value for vertical marker */
  currentT: number;
  /** Current order m (for dynamic T range) */
  m?: number;
  /** Whether sweep is loading */
  loading?: boolean;
  /** Error message to display */
  error?: string | null;
}

/**
 * Colors for the 3-regime theoretical view.
 * Maps to the lecture notes description:
 * - Small T (< 25): Series expansion (blue)
 * - Moderate T (25 <= T < 30+5m): erf + recurrence (green)
 * - Large T (> 30+5m): Asymptotic expansion (purple)
 */
const THEORETICAL_REGIME_COLORS = {
  small: 'rgba(59, 130, 246, 0.15)', // Blue - Series
  moderate: 'rgba(34, 197, 94, 0.15)', // Green - erf + recurrence
  large: 'rgba(139, 92, 246, 0.15)', // Purple/Violet - Asymptotic
};

/**
 * Create regime boundary shapes for the plot background.
 *
 * Shows 3 theoretical regimes from the lecture notes:
 * - Small T (< 25): Series expansion (blue)
 * - Moderate T (25 <= T < 30+5m): erf + upward recurrence (green)
 * - Large T (> 30+5m): Asymptotic expansion (purple)
 *
 * Also shows the implementation turnover point (m-dependent).
 *
 * @param turnover - The m-dependent turnover point (implementation boundary)
 * @param m - The current order m
 */
function createRegimeShapes(turnover: number, m: number, maxT: number): Partial<Shape>[] {
  const shapes: Partial<Shape>[] = [];

  // Calculate thresholds
  const smallTThreshold = THEORETICAL_SMALL_T_THRESHOLD; // 25
  const asymptoticThreshold = getAsymptoticThreshold(m); // 30 + 5*m

  // Small T regime (T < 25) - Series (blue)
  if (smallTThreshold > 0) {
    shapes.push({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: 0,
      x1: Math.min(smallTThreshold, maxT),
      y0: 0,
      y1: 1,
      fillcolor: THEORETICAL_REGIME_COLORS.small,
      line: { width: 0 },
      layer: 'below',
    });
  }

  // Moderate T regime (25 <= T < 30+5m) - erf + recurrence (green)
  if (smallTThreshold < maxT) {
    const moderateEnd = Math.min(asymptoticThreshold, maxT);
    if (smallTThreshold < moderateEnd) {
      shapes.push({
        type: 'rect',
        xref: 'x',
        yref: 'paper',
        x0: smallTThreshold,
        x1: moderateEnd,
        y0: 0,
        y1: 1,
        fillcolor: THEORETICAL_REGIME_COLORS.moderate,
        line: { width: 0 },
        layer: 'below',
      });
    }
  }

  // Large T regime (T > 30+5m) - Asymptotic (purple)
  // Only show if asymptotic threshold is within visible range
  if (asymptoticThreshold < maxT) {
    shapes.push({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: asymptoticThreshold,
      x1: maxT,
      y0: 0,
      y1: 1,
      fillcolor: THEORETICAL_REGIME_COLORS.large,
      line: { width: 0 },
      layer: 'below',
    });
  }

  // Vertical line at theoretical Small T threshold (T=25)
  if (smallTThreshold > 0 && smallTThreshold < maxT) {
    shapes.push({
      type: 'line',
      xref: 'x',
      yref: 'paper',
      x0: smallTThreshold,
      x1: smallTThreshold,
      y0: 0,
      y1: 1,
      line: {
        color: 'rgba(59, 130, 246, 0.5)', // Blue
        width: 1,
        dash: 'dot',
      },
    });
  }

  // Vertical line at asymptotic threshold (30+5m) if within visible range
  if (asymptoticThreshold < maxT) {
    shapes.push({
      type: 'line',
      xref: 'x',
      yref: 'paper',
      x0: asymptoticThreshold,
      x1: asymptoticThreshold,
      y0: 0,
      y1: 1,
      line: {
        color: 'rgba(139, 92, 246, 0.5)', // Purple
        width: 1,
        dash: 'dot',
      },
    });
  }

  // Also show implementation turnover point if > 0 and visible
  if (turnover > 0 && turnover < maxT) {
    shapes.push({
      type: 'line',
      xref: 'x',
      yref: 'paper',
      x0: turnover,
      x1: turnover,
      y0: 0,
      y1: 1,
      line: {
        color: 'rgba(100, 116, 139, 0.7)', // Gray
        width: 2,
        dash: 'solid',
      },
    });
  }

  return shapes;
}

/**
 * Create regime boundary annotations.
 *
 * Shows annotations for:
 * - Text labels inside each regime region (accessibility: not color-only)
 * - Theoretical Small T threshold (T=25)
 * - Asymptotic threshold (30+5m) if visible
 * - Implementation turnover point if > 0
 *
 * @param turnover - The m-dependent turnover point (implementation)
 * @param m - The current order m
 */
function createRegimeAnnotations(turnover: number, m: number, maxT: number): Partial<PlotlyAnnotation>[] {
  const annotations: Partial<PlotlyAnnotation>[] = [];
  const asymptoticThreshold = getAsymptoticThreshold(m);

  // --- Regime text labels (inside each colored region, for accessibility) ---

  // "Series" label in the blue region (centered between 0 and min(25, maxT))
  const seriesEnd = Math.min(THEORETICAL_SMALL_T_THRESHOLD, maxT);
  if (seriesEnd > 0) {
    annotations.push({
      x: seriesEnd / 2,
      y: 0.97,
      xref: 'x',
      yref: 'paper',
      text: 'Series',
      showarrow: false,
      font: { size: 10, color: 'rgba(59, 130, 246, 0.7)' },
      xanchor: 'center',
      yanchor: 'top',
    });
  }

  // "Recurrence" label in the green region (centered between 25 and min(30+5m, maxT))
  if (THEORETICAL_SMALL_T_THRESHOLD < maxT) {
    const moderateEnd = Math.min(asymptoticThreshold, maxT);
    if (THEORETICAL_SMALL_T_THRESHOLD < moderateEnd) {
      annotations.push({
        x: (THEORETICAL_SMALL_T_THRESHOLD + moderateEnd) / 2,
        y: 0.97,
        xref: 'x',
        yref: 'paper',
        text: 'Recurrence',
        showarrow: false,
        font: { size: 10, color: 'rgba(34, 197, 94, 0.7)' },
        xanchor: 'center',
        yanchor: 'top',
      });
    }
  }

  // "Asymptotic" label in the violet region (centered between 30+5m and maxT)
  if (asymptoticThreshold < maxT) {
    annotations.push({
      x: (asymptoticThreshold + maxT) / 2,
      y: 0.97,
      xref: 'x',
      yref: 'paper',
      text: 'Asymptotic',
      showarrow: false,
      font: { size: 10, color: 'rgba(139, 92, 246, 0.7)' },
      xanchor: 'center',
      yanchor: 'top',
    });
  }

  // --- Boundary threshold annotations ---

  // Annotation at T=25 (theoretical series/moderate boundary)
  if (THEORETICAL_SMALL_T_THRESHOLD < maxT) {
    annotations.push({
      x: THEORETICAL_SMALL_T_THRESHOLD,
      y: 1.02,
      xref: 'x',
      yref: 'paper',
      text: 'T=25',
      showarrow: false,
      font: { size: 9, color: '#3b82f6' },
      xanchor: 'center',
    });
  }

  // Annotation at asymptotic threshold (30+5m) if visible
  if (asymptoticThreshold < maxT) {
    annotations.push({
      x: asymptoticThreshold,
      y: 1.02,
      xref: 'x',
      yref: 'paper',
      text: `30+5(${m})=${asymptoticThreshold}`,
      showarrow: false,
      font: { size: 9, color: '#8b5cf6' },
      xanchor: 'center',
    });
  }

  // Annotation for implementation turnover point if > 0 and visible
  if (turnover > 0 && turnover < maxT) {
    annotations.push({
      x: turnover,
      y: 0.95,
      xref: 'x',
      yref: 'paper',
      text: `impl(${m})=${turnover.toFixed(1)}`,
      showarrow: false,
      font: { size: 9, color: '#64748b' },
      xanchor: 'left',
      bgcolor: 'rgba(255, 255, 255, 0.8)',
      borderpad: 2,
    });
  }

  return annotations;
}

/**
 * Create vertical marker shape at current T value.
 */
function createTMarkerShape(currentT: number): Partial<Shape> {
  return {
    type: 'line',
    xref: 'x',
    yref: 'paper',
    x0: currentT,
    x1: currentT,
    y0: 0,
    y1: 1,
    line: {
      color: '#dc2626', // Red
      width: 2,
      dash: 'solid',
    },
  };
}

/**
 * Create annotation for current T marker.
 */
function createTMarkerAnnotation(currentT: number, currentValue: number | undefined): Partial<PlotlyAnnotation> {
  return {
    x: currentT,
    y: 1.06,
    xref: 'x',
    yref: 'paper',
    text: currentValue !== undefined
      ? `T=${currentT}<br>F=${currentValue.toExponential(4)}`
      : `T=${currentT}`,
    showarrow: false,
    font: { size: 10, color: '#dc2626' },
    xanchor: 'center',
    bgcolor: 'white',
    bordercolor: '#dc2626',
    borderwidth: 1,
    borderpad: 2,
  };
}

/**
 * BoysChart - Plotly chart for Boys function curves.
 *
 * Displays F_m(T) vs T with regime boundaries and current T marker.
 *
 * @example
 * ```tsx
 * <BoysChart
 *   sweepData={sweepData}
 *   logScale={false}
 *   currentT={5.0}
 *   loading={isComputing}
 * />
 * ```
 */
export function BoysChart({
  sweepData,
  logScale,
  currentT,
  m = 0,
  loading = false,
  error = null,
}: BoysChartProps) {
  const maxT = getMaxTValue(m);
  // Extract x and y values from sweep data
  const { xValues, yValues, currentValue } = useMemo(() => {
    if (!sweepData || !sweepData.results.length) {
      return { xValues: [], yValues: [], currentValue: undefined };
    }

    const x: number[] = [];
    const y: number[] = [];
    let closestValue: number | undefined;
    let closestDist = Infinity;

    for (const result of sweepData.results) {
      x.push(result.t);
      y.push(result.value);

      // Find the value closest to currentT
      const dist = Math.abs(result.t - currentT);
      if (dist < closestDist) {
        closestDist = dist;
        closestValue = result.value;
      }
    }

    return { xValues: x, yValues: y, currentValue: closestValue };
  }, [sweepData, currentT]);

  // Calculate y-axis range for consistent display
  const yRange = useMemo(() => {
    if (!yValues.length) return undefined;

    const minY = Math.min(...yValues);
    const maxY = Math.max(...yValues);

    if (logScale) {
      // For log scale, ensure we don't go below a minimum value
      const logMin = Math.max(minY, 1e-15);
      return [Math.log10(logMin) - 0.5, Math.log10(maxY) + 0.5];
    }

    // Add 10% padding for linear scale
    const padding = (maxY - minY) * 0.1;
    return [Math.max(0, minY - padding), maxY + padding];
  }, [yValues, logScale]);

  // Create Plotly data traces
  const data: Data[] = useMemo(() => {
    if (!xValues.length) return [];

    return [
      {
        x: xValues,
        y: yValues,
        type: 'scatter',
        mode: 'lines',
        name: sweepData ? `F_${sweepData.m}(T)` : 'F_m(T)',
        line: {
          color: '#3b82f6',
          width: 2,
        },
        hovertemplate: 'T = %{x:.2f}<br>F = %{y:.6e}<extra></extra>',
      },
    ];
  }, [xValues, yValues, sweepData]);

  // Get the current order m and its turnover point
  const currentM = sweepData?.m ?? 0;
  const turnover = getTurnoverPoint(currentM);

  // Create Plotly layout
  const layout: Partial<Layout> = useMemo(() => {
    const shapes: Partial<Shape>[] = [
      ...createRegimeShapes(turnover, currentM, maxT),
      createTMarkerShape(currentT),
    ];

    const annotations: Partial<PlotlyAnnotation>[] = [
      ...createRegimeAnnotations(turnover, currentM, maxT),
      createTMarkerAnnotation(currentT, currentValue),
    ];

    return {
      title: {
        text: sweepData ? `Boys Function F<sub>${sweepData.m}</sub>(T)` : 'Boys Function',
        font: { size: 16 },
      },
      xaxis: {
        title: { text: 'T' },
        range: [0, maxT],
        showgrid: true,
        gridcolor: 'rgba(0, 0, 0, 0.1)',
      },
      yaxis: {
        title: { text: sweepData ? `F<sub>${sweepData.m}</sub>(T)` : 'F<sub>m</sub>(T)' },
        type: logScale ? 'log' : 'linear',
        range: yRange,
        showgrid: true,
        gridcolor: 'rgba(0, 0, 0, 0.1)',
        exponentformat: 'e',
      },
      shapes,
      annotations,
      margin: { t: 60, r: 30, b: 50, l: 70 },
      showlegend: false,
      hovermode: 'closest',
    };
  }, [sweepData, logScale, currentT, currentValue, yRange, turnover, currentM, maxT]);

  // SVG export handler
  const handleExportSvg = useCallback(() => {
    const plotEl = document.getElementById(SWEEP_PLOT_DIV_ID);
    if (plotEl) {
      Plotly.downloadImage(plotEl as unknown as Plotly.PlotlyHTMLElement, {
        format: 'svg',
        filename: `boys-F${sweepData?.m ?? m}-sweep`,
        width: 800,
        height: 500,
      });
    }
  }, [sweepData, m]);

  return (
    <div>
      <PlotPanel
        data={data}
        layout={layout}
        loading={loading}
        error={error}
        ariaLabel={`Boys function chart for order m=${sweepData?.m ?? 0}`}
        minHeight={400}
        className="w-full"
        divId={SWEEP_PLOT_DIV_ID}
      />

      {/* Export SVG button */}
      {sweepData && sweepData.results.length > 0 && (
        <div className="flex justify-end mt-1">
          <button
            type="button"
            onClick={handleExportSvg}
            className="text-xs text-slate-500 hover:text-slate-700 flex items-center gap-1"
            aria-label="Export sweep plot as SVG"
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" />
            </svg>
            Export SVG
          </button>
        </div>
      )}
    </div>
  );
}

export default BoysChart;
