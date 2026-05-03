/**
 * RadialProfilePlot - Interactive radial profile visualization for contracted shells.
 *
 * Displays the radial function R(r) = r^l * sum_k { d_k * exp(-alpha_k * r^2) }
 * for a contracted Gaussian shell. Shows the contracted sum as a solid line
 * and individual primitive contributions as dashed lines.
 *
 * Students can toggle between linear and log y-axis scales to examine
 * both the peak structure and the tail behavior of diffuse Gaussians.
 *
 * @module components/basisExplorer/RadialProfilePlot
 */

import { useMemo, useCallback, memo } from 'react';
import Plotly from 'plotly.js-cartesian-dist';
import { PlotPanel } from '../common/PlotPanel';
import type { Data, Layout } from 'plotly.js';
import type { RadialProfileResult } from '../../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * Props for RadialProfilePlot.
 */
export interface RadialProfilePlotProps {
  /** Radial profile data from worker, or null if not yet computed */
  profileData: RadialProfileResult | null;
  /** Whether data is currently being computed */
  isLoading: boolean;
  /** Y-axis scale mode */
  yAxisScale: 'linear' | 'log';
  /** Callback when scale mode changes */
  onScaleChange: (scale: 'linear' | 'log') => void;
  /** Minimum height in pixels (default: 350) */
  minHeight?: number;
  /**
   * Modified profile values from client-side evaluation (when exponents are adjusted).
   * If provided and exponentsModified is true, shown as an overlaid trace alongside the original.
   *
   * @see US-053 Exponent Sliders
   */
  modifiedProfileValues?: number[] | null;
  /**
   * Whether exponents have been modified (controls whether overlay is shown).
   *
   * @see US-053 Exponent Sliders
   */
  exponentsModified?: boolean;
}

// ============================================================================
// Constants
// ============================================================================

/**
 * Color palette for primitive Gaussian traces.
 * Six visually distinct colors for differentiating primitives.
 */
const PRIMITIVE_COLORS = [
  '#94a3b8', // slate-400
  '#f97316', // orange-500
  '#22c55e', // green-500
  '#a855f7', // purple-500
  '#ef4444', // red-500
  '#06b6d4', // cyan-500
];

/** Color for the contracted sum trace */
const CONTRACTED_COLOR = '#2563eb'; // blue-600

/** Color for the modified profile overlay trace (US-053) */
const MODIFIED_COLOR = '#ea580c'; // orange-600

// ============================================================================
// Scale Toggle Component
// ============================================================================

/**
 * Toggle button for linear/log y-axis scale.
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
 * Radial profile plot for contracted Gaussian basis functions.
 *
 * Features:
 * - Contracted sum shown as solid blue line (width 3)
 * - Individual primitives shown as dashed lines with distinct colors
 * - Tooltips display exponent (alpha) and effective coefficient (d) for each primitive
 * - Linear/log toggle for y-axis (log useful for examining tail behavior)
 * - d-type shells include a note that only the radial part is shown
 *
 * @example
 * ```tsx
 * <RadialProfilePlot
 *   profileData={profileResult}
 *   isLoading={loading}
 *   yAxisScale="linear"
 *   onScaleChange={setScale}
 * />
 * ```
 */
export const RadialProfilePlot = memo(function RadialProfilePlot({
  profileData,
  isLoading,
  yAxisScale,
  onScaleChange,
  minHeight = 350,
  modifiedProfileValues,
  exponentsModified = false,
}: RadialProfilePlotProps) {
  // Stable divId for Plotly SVG export
  const plotDivId = 'radial-profile-plot';

  // Whether the modified profile overlay should be shown
  const showModified = exponentsModified && modifiedProfileValues != null && profileData != null;

  // SVG export handler
  const handleExportSvg = useCallback(() => {
    const plotEl = document.getElementById(plotDivId);
    if (plotEl) {
      const filename = profileData
        ? `basis-radial-${profileData.angularMomentumLabel}-${profileData.nPrimitives}p`
        : 'basis-radial';
      Plotly.downloadImage(plotEl as unknown as Plotly.PlotlyHTMLElement, {
        format: 'svg',
        filename,
        width: 800,
        height: 500,
      });
    }
  }, [profileData]);

  // Prepare Plotly data traces
  const plotData = useMemo((): Data[] => {
    if (!profileData) return [];

    const traces: Data[] = [];

    // Add individual primitive traces (dashed, behind contracted)
    for (let k = 0; k < profileData.nPrimitives; k++) {
      const alpha = profileData.exponents[k];
      const dEff = profileData.effectiveCoefficients[k];
      const dRaw = profileData.rawCoefficients[k];
      const color = PRIMITIVE_COLORS[k % PRIMITIVE_COLORS.length];

      // For log scale, filter out zero/negative values to avoid -Infinity
      let rValues = profileData.rValues;
      let yValues = profileData.primitiveValues[k];

      if (yAxisScale === 'log') {
        const filtered = filterPositive(rValues, yValues);
        rValues = filtered.x;
        yValues = filtered.y;
      }

      traces.push({
        x: rValues,
        y: yValues,
        type: 'scatter',
        mode: 'lines',
        name: `Primitive ${k + 1}`,
        line: {
          color,
          width: 1.5,
          dash: 'dash',
        },
        hovertemplate:
          `<b>Primitive ${k + 1}</b><br>` +
          `r = %{x:.3f} bohr<br>` +
          `R(r) = %{y:.6e}<br>` +
          `<i>\u03B1</i> = ${alpha.toExponential(4)}<br>` +
          `d<sub>eff</sub> = ${dEff.toFixed(6)}<br>` +
          `c<sub>raw</sub> = ${dRaw.toFixed(6)}` +
          '<extra></extra>',
      });
    }

    // Add contracted sum trace (solid, on top)
    // When modified profile is shown, relabel this as "Original"
    let rValues = profileData.rValues;
    let yValues = profileData.contractedValues;

    if (yAxisScale === 'log') {
      const filtered = filterPositive(rValues, yValues);
      rValues = filtered.x;
      yValues = filtered.y;
    }

    traces.push({
      x: rValues,
      y: yValues,
      type: 'scatter',
      mode: 'lines',
      name: showModified ? 'Original' : 'Contracted',
      line: {
        color: CONTRACTED_COLOR,
        width: showModified ? 2 : 3,
        dash: showModified ? 'dash' : 'solid',
      },
      hovertemplate:
        `<b>${showModified ? 'Original' : 'Contracted Sum'}</b><br>` +
        'r = %{x:.3f} bohr<br>' +
        'R(r) = %{y:.6e}' +
        '<extra></extra>',
    });

    // Add modified profile overlay (US-053)
    if (showModified && modifiedProfileValues) {
      let modR = profileData.rValues;
      let modY = modifiedProfileValues;

      if (yAxisScale === 'log') {
        const filtered = filterPositive(modR, modY);
        modR = filtered.x;
        modY = filtered.y;
      }

      traces.push({
        x: modR,
        y: modY,
        type: 'scatter',
        mode: 'lines',
        name: 'Modified',
        line: {
          color: MODIFIED_COLOR,
          width: 3,
        },
        hovertemplate:
          '<b>Modified</b><br>' +
          'r = %{x:.3f} bohr<br>' +
          'R(r) = %{y:.6e}' +
          '<extra></extra>',
      });
    }

    return traces;
  }, [profileData, yAxisScale, showModified, modifiedProfileValues]);

  // Prepare Plotly layout
  const layout = useMemo((): Partial<Layout> => {
    let titleText = profileData
      ? `Radial Profile: ${profileData.angularMomentumLabel}-type shell (${profileData.nPrimitives} primitives)`
      : 'Radial Profile';

    // Add note for d-type shells
    if (profileData && profileData.angularMomentum === 2) {
      titleText += '<br><span style="font-size:10px;color:#64748b">Radial part only (angular structure not shown)</span>';
    }

    return {
      title: {
        text: titleText,
        font: { size: 14 },
      },
      xaxis: {
        title: { text: 'r (bohr)', font: { size: 12 } },
        rangemode: 'tozero',
        zeroline: true,
      },
      yaxis: {
        title: { text: 'R(r) (a.u.)', font: { size: 12 } },
        type: yAxisScale === 'log' ? 'log' : 'linear',
        exponentformat: 'e',
        zeroline: yAxisScale === 'linear',
      },
      margin: { t: profileData?.angularMomentum === 2 ? 65 : 50, r: 20, b: 50, l: 70 },
      showlegend: true,
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
  }, [profileData, yAxisScale]);

  // Empty state: no data and not loading
  if (!profileData && !isLoading) {
    return (
      <div
        className="bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center"
        style={{ minHeight: `${minHeight}px` }}
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
                d="M3 13.125C3 12.504 3.504 12 4.125 12h2.25c.621 0 1.125.504 1.125 1.125v6.75C7.5 20.496 6.996 21 6.375 21h-2.25A1.125 1.125 0 013 19.875v-6.75zM9.75 8.625c0-.621.504-1.125 1.125-1.125h2.25c.621 0 1.125.504 1.125 1.125v11.25c0 .621-.504 1.125-1.125 1.125h-2.25a1.125 1.125 0 01-1.125-1.125V8.625zM16.5 4.125c0-.621.504-1.125 1.125-1.125h2.25C20.496 3 21 3.504 21 4.125v15.75c0 .621-.504 1.125-1.125 1.125h-2.25a1.125 1.125 0 01-1.125-1.125V4.125z"
              />
            </svg>
          </div>
          <p className="text-sm font-medium">Radial Profile</p>
          <p className="text-xs mt-1">
            Select a shell from the list to view its radial profile
          </p>
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-2">
      {/* Scale toggle */}
      <div className="flex items-center justify-between">
        <h3 className="text-sm font-semibold text-slate-700">Radial Profile</h3>
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
          profileData
            ? `Radial profile of ${profileData.angularMomentumLabel}-type shell with ${profileData.nPrimitives} primitives`
            : 'Radial profile plot'
        }
      />

      {/* Export SVG button */}
      {profileData && (
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

      {/* Educational note */}
      {profileData && (
        <div className="bg-blue-50 rounded-lg px-4 py-3 border border-blue-100">
          <p className="text-xs text-blue-800">
            {profileData.angularMomentum === 0 && (
              <>
                <strong>s-type (l=0):</strong> The radial function is maximal at the nucleus (r=0)
                and decays monotonically. Each primitive contributes a Gaussian with different
                width (controlled by exponent <i>&alpha;</i>). Tight primitives (large <i>&alpha;</i>)
                capture the cusp near the nucleus, while diffuse primitives (small <i>&alpha;</i>)
                describe the tail.
              </>
            )}
            {profileData.angularMomentum === 1 && (
              <>
                <strong>p-type (l=1):</strong> The r<sup>1</sup> prefactor forces R(0)&nbsp;=&nbsp;0.
                The profile rises from zero, reaches a maximum, then decays. This reflects the
                nodal structure of p orbitals, which have zero electron density at the nucleus.
              </>
            )}
            {profileData.angularMomentum === 2 && (
              <>
                <strong>d-type (l=2):</strong> The r<sup>2</sup> prefactor pushes the maximum farther
                from the nucleus. Only the radial part R(r) = r<sup>2</sup> &middot; &sum; d<sub>k</sub>
                exp(-&alpha;<sub>k</sub>r<sup>2</sup>) is shown; the full angular structure (e.g., d<sub>xy</sub>,
                d<sub>z<sup>2</sup></sub>) is not visualized.
              </>
            )}
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
 *
 * Plotly cannot display zero or negative values on a log axis.
 * This filters out such points to prevent rendering artifacts.
 */
function filterPositive(
  x: number[],
  y: number[]
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

export default RadialProfilePlot;
