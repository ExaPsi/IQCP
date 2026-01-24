/**
 * RysErrorChart - Plotly chart for Rys error vs. quadrature order.
 *
 * Displays the maximum reconstruction error as a function of quadrature order,
 * with a horizontal line at the target accuracy threshold.
 *
 * @module components/rys/RysErrorChart
 */

import { useMemo } from 'react';
import { PlotPanel } from '../common';
import type { Data, Layout, Shape, Annotations } from 'plotly.js';
import type { RysErrorCurveResult } from '../../worker/protocol';
import { type TargetAccuracy, targetToNumber, findRecommendedOrder } from '../../stores/rysStore';

// Use Plotly's Annotations type for individual annotation objects
type PlotlyAnnotation = Annotations;

/**
 * Props for RysErrorChart component.
 */
export interface RysErrorChartProps {
  /** Error curve computation results */
  errorCurve: RysErrorCurveResult | null;
  /** Target accuracy for threshold line */
  target: TargetAccuracy;
  /** Whether computation is loading */
  loading?: boolean;
  /** Error message to display */
  error?: string | null;
}

/**
 * Create horizontal threshold line at target accuracy.
 */
function createThresholdLine(target: TargetAccuracy): Partial<Shape> {
  const targetValue = targetToNumber(target);
  return {
    type: 'line',
    xref: 'paper',
    yref: 'y',
    x0: 0,
    x1: 1,
    y0: targetValue,
    y1: targetValue,
    line: {
      color: '#dc2626', // Red
      width: 2,
      dash: 'dash',
    },
  };
}

/**
 * Create annotation for threshold line label.
 */
function createThresholdAnnotation(target: TargetAccuracy): Partial<PlotlyAnnotation> {
  const targetValue = targetToNumber(target);
  return {
    x: 1,
    y: Math.log10(targetValue),
    xref: 'paper',
    yref: 'y',
    text: `Target: ${target}`,
    showarrow: false,
    font: { size: 10, color: '#dc2626' },
    xanchor: 'right',
    yanchor: 'bottom',
    bgcolor: 'white',
    borderpad: 2,
  };
}

/**
 * Create marker annotation for recommended order.
 */
function createRecommendedAnnotation(
  n: number,
  errorValue: number
): Partial<PlotlyAnnotation> {
  return {
    x: n,
    y: Math.log10(errorValue),
    xref: 'x',
    yref: 'y',
    text: `n=${n}`,
    showarrow: true,
    arrowhead: 2,
    arrowsize: 1,
    arrowcolor: '#22c55e',
    ax: 30,
    ay: -30,
    font: { size: 11, color: '#22c55e' },
    bgcolor: 'white',
    bordercolor: '#22c55e',
    borderwidth: 1,
    borderpad: 3,
  };
}

/**
 * RysErrorChart - Plotly chart for error vs. quadrature order.
 *
 * Displays:
 * - Max reconstruction error vs n (log scale y-axis)
 * - Horizontal threshold line at target accuracy
 * - Marker highlighting the recommended order
 *
 * @example
 * ```tsx
 * <RysErrorChart
 *   errorCurve={errorCurve}
 *   target="1e-6"
 *   loading={isComputing}
 * />
 * ```
 */
export function RysErrorChart({
  errorCurve,
  target,
  loading = false,
  error = null,
}: RysErrorChartProps) {
  // Find recommended order
  const recommendedOrder = useMemo(() => {
    if (!errorCurve) return null;
    return findRecommendedOrder(errorCurve, target);
  }, [errorCurve, target]);

  // Extract data points
  const { xValues, yValues, recommendedError } = useMemo(() => {
    if (!errorCurve || !errorCurve.points.length) {
      return { xValues: [], yValues: [], recommendedError: undefined };
    }

    const x: number[] = [];
    const y: number[] = [];
    let recError: number | undefined;

    for (const point of errorCurve.points) {
      x.push(point.n);
      y.push(point.maxError);

      if (recommendedOrder && point.n === recommendedOrder) {
        recError = point.maxError;
      }
    }

    return { xValues: x, yValues: y, recommendedError: recError };
  }, [errorCurve, recommendedOrder]);

  // Create Plotly data traces
  const data: Data[] = useMemo(() => {
    if (!xValues.length) return [];

    const traces: Data[] = [
      // Main error line
      {
        x: xValues,
        y: yValues,
        type: 'scatter',
        mode: 'lines+markers',
        name: 'Max Error',
        line: {
          color: '#3b82f6',
          width: 2,
        },
        marker: {
          color: '#3b82f6',
          size: 8,
        },
        hovertemplate: 'n = %{x}<br>Error = %{y:.2e}<extra></extra>',
      },
    ];

    // Add marker for recommended order
    if (recommendedOrder && recommendedError !== undefined) {
      traces.push({
        x: [recommendedOrder],
        y: [recommendedError],
        type: 'scatter',
        mode: 'markers',
        name: 'Recommended',
        marker: {
          color: '#22c55e',
          size: 14,
          symbol: 'circle',
          line: {
            color: 'white',
            width: 2,
          },
        },
        hovertemplate: `Recommended: n = ${recommendedOrder}<br>Error = ${recommendedError.toExponential(2)}<extra></extra>`,
      });
    }

    return traces;
  }, [xValues, yValues, recommendedOrder, recommendedError]);

  // Create Plotly layout
  const layout: Partial<Layout> = useMemo(() => {
    const shapes: Partial<Shape>[] = [createThresholdLine(target)];

    const annotations: Partial<PlotlyAnnotation>[] = [
      createThresholdAnnotation(target),
    ];

    // Add annotation for recommended order
    if (recommendedOrder && recommendedError !== undefined) {
      annotations.push(createRecommendedAnnotation(recommendedOrder, recommendedError));
    }

    // Calculate y-axis range
    let yMin = targetToNumber(target) / 100; // Show some space below target
    let yMax = 1; // Max error shouldn't exceed 1

    if (yValues.length) {
      const dataMin = Math.min(...yValues);
      const dataMax = Math.max(...yValues);
      yMin = Math.min(yMin, dataMin / 10);
      yMax = Math.max(yMax, dataMax * 10);
    }

    return {
      title: {
        text: 'Reconstruction Error vs Quadrature Order',
        font: { size: 14 },
      },
      xaxis: {
        title: { text: 'Quadrature Order (n)' },
        range: [0.5, 10.5],
        tickmode: 'linear',
        tick0: 1,
        dtick: 1,
        showgrid: true,
        gridcolor: 'rgba(0, 0, 0, 0.1)',
      },
      yaxis: {
        title: { text: 'Max Reconstruction Error' },
        type: 'log',
        range: [Math.log10(yMin), Math.log10(yMax)],
        showgrid: true,
        gridcolor: 'rgba(0, 0, 0, 0.1)',
        exponentformat: 'e',
      },
      shapes,
      annotations,
      margin: { t: 50, r: 30, b: 50, l: 70 },
      showlegend: false,
      hovermode: 'closest',
    };
  }, [target, yValues, recommendedOrder, recommendedError]);

  return (
    <PlotPanel
      data={data}
      layout={layout}
      loading={loading}
      error={error}
      ariaLabel={`Rys error curve chart showing error vs quadrature order`}
      minHeight={350}
      className="w-full"
    />
  );
}

export default RysErrorChart;
