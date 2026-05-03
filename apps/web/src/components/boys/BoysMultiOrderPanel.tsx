/**
 * BoysMultiOrderPanel - Multi-order comparison view for Boys function.
 *
 * Displays F_0(T) through F_m(T) at the current T value in both
 * table and bar chart form. This view helps students observe the
 * key property that F_m(T) decreases monotonically with increasing m.
 *
 * @module components/boys/BoysMultiOrderPanel
 */

import { useMemo, useCallback } from 'react';
import Plotly from 'plotly.js-cartesian-dist';
import { PlotPanel } from '../common';
import type { Data, Layout } from 'plotly.js';
import type { BoysEvalAllResult, BoysMethod } from '../../worker/protocol';

/**
 * Props for BoysMultiOrderPanel component.
 */
export interface BoysMultiOrderPanelProps {
  /** Multi-order computation results */
  data: BoysEvalAllResult | null;
  /** Whether computation is in progress */
  loading?: boolean;
  /** Error message to display */
  error?: string | null;
  /** Current T value (for display) */
  currentT: number;
}

/**
 * Method badge colors for the table.
 */
const METHOD_COLORS: Record<BoysMethod, { bg: string; text: string }> = {
  zero: { bg: 'bg-amber-100', text: 'text-amber-800' },
  series: { bg: 'bg-blue-100', text: 'text-blue-800' },
  recurrence: { bg: 'bg-green-100', text: 'text-green-800' },
};

/**
 * Bar colors corresponding to each computation method.
 */
const METHOD_BAR_COLORS: Record<BoysMethod, string> = {
  zero: '#f59e0b',     // amber
  series: '#3b82f6',   // blue
  recurrence: '#22c55e', // green
};

/**
 * Plotly div ID for SVG export.
 */
const PLOT_DIV_ID = 'boys-multi-order-plot';

/**
 * Format a number for display in the table.
 */
function formatValue(value: number): string {
  if (Math.abs(value) < 1e-4 && value !== 0) {
    return value.toExponential(10);
  }
  if (Math.abs(value) >= 0.0001 && Math.abs(value) < 1e6) {
    return value.toPrecision(11);
  }
  return value.toExponential(10);
}

/**
 * BoysMultiOrderPanel - Table and bar chart showing F_m(T) for all orders.
 *
 * Pedagogically valuable because students can observe:
 * 1. F_m(T) decreases monotonically with increasing m
 * 2. Different computation methods are used for different m at the same T
 * 3. The rate of decrease depends on T
 *
 * @example
 * ```tsx
 * <BoysMultiOrderPanel
 *   data={multiOrderData}
 *   loading={isComputing}
 *   currentT={5.0}
 * />
 * ```
 */
export function BoysMultiOrderPanel({
  data,
  loading = false,
  error = null,
  currentT,
}: BoysMultiOrderPanelProps) {
  // Build bar chart data
  const plotData: Data[] = useMemo(() => {
    if (!data || !data.results.length) return [];

    const mValues = data.results.map((r) => `m=${r.m}`);
    const values = data.results.map((r) => r.value);
    const colors = data.results.map((r) => METHOD_BAR_COLORS[r.method]);

    return [
      {
        x: mValues,
        y: values,
        type: 'bar',
        marker: {
          color: colors,
        },
        hovertemplate: '%{x}<br>F = %{y:.6e}<extra></extra>',
      },
    ];
  }, [data]);

  const plotLayout: Partial<Layout> = useMemo(() => {
    return {
      title: {
        text: `Multi-Order Comparison at T = ${currentT.toFixed(2)}`,
        font: { size: 14 },
      },
      xaxis: {
        title: { text: 'Order m' },
      },
      yaxis: {
        title: { text: 'F_m(T)' },
        exponentformat: 'e',
      },
      margin: { t: 50, r: 30, b: 50, l: 70 },
      showlegend: false,
      hovermode: 'closest',
    };
  }, [currentT]);

  // SVG export handler
  const handleExportSvg = useCallback(() => {
    const plotEl = document.getElementById(PLOT_DIV_ID);
    if (plotEl) {
      Plotly.downloadImage(plotEl as unknown as Plotly.PlotlyHTMLElement, {
        format: 'svg',
        filename: `boys-multi-order-T${currentT.toFixed(2)}`,
        width: 800,
        height: 500,
      });
    }
  }, [currentT]);

  return (
    <div className="space-y-4">
      {/* Bar chart */}
      <PlotPanel
        data={plotData}
        layout={plotLayout}
        loading={loading}
        error={error}
        ariaLabel={`Boys function multi-order comparison bar chart at T=${currentT.toFixed(2)}`}
        minHeight={300}
        className="w-full"
        divId={PLOT_DIV_ID}
      />

      {/* Export SVG button */}
      {data && data.results.length > 0 && (
        <div className="flex justify-end">
          <button
            type="button"
            onClick={handleExportSvg}
            className="text-xs text-slate-500 hover:text-slate-700 flex items-center gap-1"
            aria-label="Export multi-order plot as SVG"
          >
            <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" />
            </svg>
            Export SVG
          </button>
        </div>
      )}

      {/* Data table */}
      {data && data.results.length > 0 && (
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
          <div className="p-4 border-b border-slate-100">
            <h3 className="text-sm font-semibold text-slate-700">
              F<sub>m</sub>(T) values at T = {currentT.toFixed(4)}
            </h3>
            <p className="text-xs text-slate-500 mt-1">
              The Boys function decreases monotonically with increasing order m.
            </p>
          </div>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="bg-slate-50 text-slate-600">
                  <th className="px-4 py-2 text-left font-medium">m</th>
                  <th className="px-4 py-2 text-right font-medium">F<sub>m</sub>(T)</th>
                  <th className="px-4 py-2 text-center font-medium">Method</th>
                  <th className="px-4 py-2 text-right font-medium">Terms</th>
                </tr>
              </thead>
              <tbody>
                {data.results.map((result, index) => (
                  <tr
                    key={result.m}
                    className={`border-t border-slate-100 ${
                      index % 2 === 0 ? 'bg-white' : 'bg-slate-50/50'
                    }`}
                  >
                    <td className="px-4 py-2 font-mono font-medium text-slate-800">
                      {result.m}
                    </td>
                    <td className="px-4 py-2 font-mono text-right text-slate-700">
                      {formatValue(result.value)}
                    </td>
                    <td className="px-4 py-2 text-center">
                      <span
                        className={`inline-block px-2 py-0.5 rounded text-xs font-medium ${
                          METHOD_COLORS[result.method].bg
                        } ${METHOD_COLORS[result.method].text}`}
                      >
                        {result.method}
                      </span>
                    </td>
                    <td className="px-4 py-2 font-mono text-right text-slate-600">
                      {result.termsCount}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          {/* Legend */}
          <div className="px-4 py-3 bg-slate-50 border-t border-slate-100">
            <div className="flex flex-wrap gap-4 text-xs text-slate-500">
              <div className="flex items-center gap-1.5">
                <div className="w-3 h-3 rounded bg-amber-100 border border-amber-300" />
                <span>Zero (T = 0)</span>
              </div>
              <div className="flex items-center gap-1.5">
                <div className="w-3 h-3 rounded bg-blue-100 border border-blue-300" />
                <span>Series expansion</span>
              </div>
              <div className="flex items-center gap-1.5">
                <div className="w-3 h-3 rounded bg-green-100 border border-green-300" />
                <span>erf + recurrence</span>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Empty state */}
      {!loading && !error && (!data || data.results.length === 0) && (
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 text-center">
          <p className="text-slate-500">
            Select multi-order mode to compare F<sub>0</sub>(T) through F<sub>m</sub>(T)
            at the current T value.
          </p>
        </div>
      )}
    </div>
  );
}

export default BoysMultiOrderPanel;
