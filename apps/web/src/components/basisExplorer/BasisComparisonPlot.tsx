/**
 * BasisComparisonPlot - Overlaid radial profiles from multiple basis sets.
 *
 * Displays the contracted radial function R(r) for the same shell type
 * across 2-5 selected basis sets. Each basis set gets a distinct color.
 * For split-valence bases, multiple shells (inner + outer valence) are
 * shown as separate traces with the same color but different dash patterns.
 *
 * This plot is the core visualization for AC1 (side-by-side profiles)
 * and AC2 (split-valence inner/outer distinction).
 *
 * @module components/basisExplorer/BasisComparisonPlot
 * @see US-052 Basis Set Comparison Mode
 */

import { useMemo, useCallback, memo } from 'react';
import Plotly from 'plotly.js-cartesian-dist';
import { PlotPanel } from '../common/PlotPanel';
import { COMPARISON_COLORS, SHELL_DASH_PATTERNS } from './comparisonColors';
import type { Data, Layout } from 'plotly.js';
import type { ComparisonProfile } from '../../hooks/useRadialComparison';
import { ELEMENTS } from './basisData';

// ============================================================================
// Types
// ============================================================================

export interface BasisComparisonPlotProps {
  /** Map of basis name -> radial profile result(s) for matching shell type */
  profiles: Map<string, ComparisonProfile[]>;
  /** Shell type being compared (e.g., "1s", "2p") */
  shellTypeLabel: string;
  /** Atomic number of the element being examined */
  atomicNumber: number;
  /** Whether any profile is still loading */
  isLoading: boolean;
  /** Y-axis scale mode */
  yAxisScale: 'linear' | 'log';
  /** Callback when scale mode changes */
  onScaleChange: (scale: 'linear' | 'log') => void;
  /** Minimum height in pixels */
  minHeight?: number;
}

// ============================================================================
// Component
// ============================================================================

/**
 * Overlaid comparison plot for radial profiles across basis sets.
 *
 * Features:
 * - One trace per (basis, shell) combination
 * - Basis sets identified by color from COMPARISON_COLORS
 * - Inner valence: solid line; outer valence: dashed line
 * - Hover shows basis name, role, r, and R(r)
 * - Legend grouped by basis set
 * - Linear/log toggle for y-axis
 */
export const BasisComparisonPlot = memo(function BasisComparisonPlot({
  profiles,
  shellTypeLabel,
  atomicNumber,
  isLoading,
  yAxisScale,
  onScaleChange,
  minHeight = 380,
}: BasisComparisonPlotProps) {
  // Stable divId for Plotly SVG export
  const plotDivId = 'basis-comparison-plot';

  // Get element symbol for title
  const elementSymbol = ELEMENTS.find((e) => e.z === atomicNumber)?.symbol ?? `Z=${atomicNumber}`;

  // SVG export handler
  const handleExportSvg = useCallback(() => {
    const plotEl = document.getElementById(plotDivId);
    if (plotEl) {
      Plotly.downloadImage(plotEl as unknown as Plotly.PlotlyHTMLElement, {
        format: 'svg',
        filename: `basis-comparison-${elementSymbol}-${shellTypeLabel}`,
        width: 800,
        height: 500,
      });
    }
  }, [elementSymbol, shellTypeLabel]);

  // Build Plotly traces
  const plotData = useMemo((): Data[] => {
    if (profiles.size === 0) return [];

    const traces: Data[] = [];

    for (const [, basisProfiles] of profiles) {
      for (let k = 0; k < basisProfiles.length; k++) {
        const cp = basisProfiles[k];
        const color = COMPARISON_COLORS[cp.basisName] ?? '#6b7280';
        const dash = SHELL_DASH_PATTERNS[k] ?? 'solid';

        // Build trace name: "STO-3G" or "3-21G inner valence"
        let traceName = cp.basisLabel;
        if (basisProfiles.length > 1) {
          traceName += ` (${cp.shellRole})`;
        }

        // Filter for log scale
        let rValues = cp.profile.rValues;
        let yValues = cp.profile.contractedValues;

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
          name: traceName,
          line: {
            color,
            width: 2.5,
            dash: dash as 'solid' | 'dash' | 'dot' | 'dashdot' | 'longdash',
          },
          hovertemplate:
            `<b>${traceName}</b><br>` +
            `r = %{x:.3f} bohr<br>` +
            `R(r) = %{y:.6e}` +
            '<extra></extra>',
        });
      }
    }

    return traces;
  }, [profiles, yAxisScale]);

  // Plotly layout
  const layout = useMemo((): Partial<Layout> => {
    return {
      title: {
        text: `${elementSymbol} ${shellTypeLabel} Radial Profile Comparison`,
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
      margin: { t: 50, r: 20, b: 50, l: 70 },
      showlegend: true,
      legend: {
        x: 1,
        y: 1,
        xanchor: 'right',
        yanchor: 'top',
        font: { size: 10 },
        bgcolor: 'rgba(255,255,255,0.9)',
        borderwidth: 1,
        bordercolor: '#e2e8f0',
      },
      hovermode: 'closest',
    };
  }, [elementSymbol, shellTypeLabel, yAxisScale]);

  // Empty state: no profiles and not loading
  if (profiles.size === 0 && !isLoading) {
    return (
      <div
        className="bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center"
        style={{ minHeight: `${minHeight}px` }}
      >
        <div className="text-center text-slate-500 px-6">
          <p className="text-sm font-medium">No comparison data</p>
          <p className="text-xs mt-1">
            Select at least two basis sets and a shell type to compare.
          </p>
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-2">
      {/* Scale toggle */}
      <div className="flex items-center justify-between">
        <h3 className="text-sm font-semibold text-slate-700">
          Basis Set Comparison
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
        ariaLabel={`Radial profile comparison of ${shellTypeLabel} shell for ${elementSymbol} across ${profiles.size} basis sets`}
      />

      {/* Export SVG button */}
      {profiles.size > 0 && (
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
    </div>
  );
});

// ============================================================================
// Scale Toggle (reused pattern from RadialProfilePlot)
// ============================================================================

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
// Utilities
// ============================================================================

/**
 * Filter out non-positive values for log-scale display.
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
