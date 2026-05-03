/**
 * PesCurvePlot - PES energy vs. coordinate visualization.
 *
 * Displays the potential energy surface (PES) as an interactive Plotly
 * scatter+line plot with:
 * - Progressive update as scan points arrive
 * - Dynamic x-axis label based on coordinate type (bond/angle/dihedral)
 * - Dual-trace overlay for rigid (dashed) + relaxed (solid) comparison
 * - Optimization step tooltip for relaxed scan points
 * - Equilibrium marker (vertical dashed line + annotation)
 * - Separate visual treatment for unconverged points
 * - Summary panel with equilibrium data
 * - Collapsible computational notes with RHF dissociation warning
 *
 * Backward compatible: when only store data is used (no rigidPoints/relaxedPoints
 * props), renders the original single-trace behavior from US-041.
 *
 * @module components/scf/PesCurvePlot
 * @see US-041 PES Curve Visualization
 * @see US-083 Extended PES Plot + Coordinate Tracking
 */

import { useMemo, useState, useCallback } from 'react';
import { PlotPanel } from '../common/PlotPanel';
import { useScfStore } from '../../stores/scfStore';
import type { Data, Layout, Shape, Annotations, PlotHoverEvent, PlotMouseEvent } from 'plotly.js';
import type { PesPoint, PesInternalPoint } from '../../worker/protocol';
import type { CoordinateType } from './CoordinateSelector';

// ============================================================================
// Constants
// ============================================================================

/** Blue color for converged points (Tailwind blue-500) */
const COLOR_CONVERGED = '#3b82f6';

/** Red color for unconverged points (Tailwind red-500) */
const COLOR_UNCONVERGED = '#ef4444';

/** Green color for equilibrium marker (Tailwind green-500) */
const COLOR_EQUILIBRIUM = '#22c55e';

/** Amber color for selected point marker (Tailwind amber-500) */
const COLOR_SELECTED = '#f59e0b';

/** Darker amber for selected point marker border (Tailwind amber-600) */
const COLOR_SELECTED_BORDER = '#d97706';

/** Bond distance threshold above which the RHF dissociation note appears (bohr) */
const RHF_DISSOCIATION_NOTE_THRESHOLD = 3.0;

/** Lighter blue for rigid scan trace (Tailwind blue-300) */
const COLOR_RIGID = '#93c5fd';

/** Darker blue for relaxed scan trace (Tailwind blue-600) */
const COLOR_RELAXED = '#2563eb';

/**
 * Dynamic x-axis labels based on coordinate type (US-083 AC1).
 */
const X_AXIS_LABELS: Record<CoordinateType, string> = {
  bond: 'Bond Length (bohr)',
  angle: 'Angle (degrees)',
  dihedral: 'Dihedral Angle (degrees)',
};

/**
 * Convert coordinate values from internal units to display units.
 *
 * Bond lengths are kept in bohr (identity). Angles and dihedrals are
 * converted from radians to degrees for display.
 *
 * @param values - Raw coordinate values (bohr or radians)
 * @param coordinateType - Type of coordinate for unit selection
 * @returns Converted values suitable for the x-axis
 */
function convertXValues(
  values: number[],
  coordinateType: CoordinateType,
): number[] {
  if (coordinateType === 'bond') return values;
  return values.map((v) => v * (180 / Math.PI));
}

// ============================================================================
// Sub-Components
// ============================================================================

/**
 * Summary panel showing equilibrium data and point statistics.
 */
function PesSummaryPanel() {
  const results = useScfStore((state) => state.pesState.results);
  const equilibrium = useScfStore((state) => state.pesState.equilibrium);
  const computeTimeMs = useScfStore((state) => state.pesState.computeTimeMs);
  const pesCoordinateConfig = useScfStore((state) => state.pesState.pesCoordinateConfig);
  const coordType = pesCoordinateConfig?.coordinateType ?? 'bond';

  const convergedCount = useMemo(
    () => results.filter((p) => p.converged).length,
    [results]
  );
  const unconvergedCount = results.length - convergedCount;

  if (results.length === 0) return null;

  const eqLabel = coordType === 'angle'
    ? 'Equilibrium Angle'
    : coordType === 'dihedral'
      ? 'Equilibrium Dihedral'
      : 'Equilibrium R';

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      <h3 className="text-sm font-semibold text-slate-700 mb-3">
        PES Scan Summary
      </h3>
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
        {equilibrium && (
          <>
            <div>
              <span className="text-slate-500 block text-xs">
                {eqLabel}
              </span>
              <span className="font-mono text-slate-800">
                {coordType === 'bond'
                  ? `${equilibrium.r_bohr.toFixed(4)} bohr`
                  : `${(equilibrium.r_bohr * (180 / Math.PI)).toFixed(1)} deg`}
              </span>
            </div>
            <div>
              <span className="text-slate-500 block text-xs">
                Equilibrium E
              </span>
              <span className="font-mono text-slate-800">
                {equilibrium.energy_hartree.toFixed(8)} Ha
              </span>
            </div>
          </>
        )}
        <div>
          <span className="text-slate-500 block text-xs">Points</span>
          <span className="font-mono text-slate-800">
            {results.length} total
            {unconvergedCount > 0 && (
              <span className="text-red-600 ml-1">
                ({unconvergedCount} unconverged)
              </span>
            )}
          </span>
        </div>
        {computeTimeMs !== null && (
          <div>
            <span className="text-slate-500 block text-xs">
              Compute Time
            </span>
            <span className="font-mono text-slate-800">
              {computeTimeMs < 1000
                ? `${computeTimeMs.toFixed(1)} ms`
                : `${(computeTimeMs / 1000).toFixed(2)} s`}
            </span>
          </div>
        )}
      </div>
    </div>
  );
}

/**
 * Collapsible computational notes panel with point detail table
 * and RHF dissociation warning.
 */
function PesNotesPanel() {
  const results = useScfStore((state) => state.pesState.results);
  const pesCoordinateConfig = useScfStore((state) => state.pesState.pesCoordinateConfig);
  const coordType = pesCoordinateConfig?.coordinateType ?? 'bond';
  const [isExpanded, setIsExpanded] = useState(false);

  // RHF dissociation warning only applies to bond length scans
  const hasLargeR = useMemo(
    () => coordType === 'bond' && results.some((p) => p.r > RHF_DISSOCIATION_NOTE_THRESHOLD),
    [results, coordType]
  );

  const toggleExpanded = useCallback(() => {
    setIsExpanded((prev) => !prev);
  }, []);

  if (results.length === 0) return null;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      {/* Header - clickable to expand/collapse */}
      <button
        type="button"
        onClick={toggleExpanded}
        className="w-full px-4 py-3 bg-slate-50 border-b border-slate-200 flex items-center justify-between text-left hover:bg-slate-100 transition-colors"
        aria-expanded={isExpanded}
        aria-controls="pes-notes-content"
      >
        <h3 className="text-sm font-semibold text-slate-700">
          Computational Notes
        </h3>
        <svg
          className={`w-4 h-4 text-slate-500 transition-transform ${isExpanded ? 'rotate-180' : ''}`}
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {/* Content */}
      {isExpanded && (
        <div id="pes-notes-content" className="p-4 space-y-4">
          {/* Point detail table */}
          <div className="overflow-x-auto">
            <table className="w-full text-xs">
              <thead>
                <tr className="border-b border-slate-200">
                  <th className="text-left py-2 px-2 font-medium text-slate-600">#</th>
                  <th className="text-right py-2 px-2 font-medium text-slate-600">
                    {coordType === 'bond' ? 'R (bohr)' : coordType === 'angle' ? 'Angle (deg)' : 'Dihedral (deg)'}
                  </th>
                  <th className="text-right py-2 px-2 font-medium text-slate-600">Energy (Ha)</th>
                  <th className="text-center py-2 px-2 font-medium text-slate-600">Converged</th>
                  <th className="text-right py-2 px-2 font-medium text-slate-600">Iterations</th>
                </tr>
              </thead>
              <tbody>
                {results.map((point, idx) => (
                  <tr
                    key={idx}
                    className={`border-b border-slate-100 ${
                      !point.converged ? 'bg-red-50' : ''
                    }`}
                  >
                    <td className="py-1.5 px-2 text-slate-500">{idx + 1}</td>
                    <td className="py-1.5 px-2 text-right font-mono text-slate-800">
                      {coordType === 'bond'
                        ? point.r.toFixed(4)
                        : (point.r * (180 / Math.PI)).toFixed(1)}
                    </td>
                    <td className="py-1.5 px-2 text-right font-mono text-slate-800">
                      {point.energy.toFixed(10)}
                    </td>
                    <td className="py-1.5 px-2 text-center">
                      {point.converged ? (
                        <span className="text-green-600" aria-label="converged">
                          &#10003;
                        </span>
                      ) : (
                        <span className="text-red-600" aria-label="unconverged">
                          &#10007;
                        </span>
                      )}
                    </td>
                    <td className="py-1.5 px-2 text-right font-mono text-slate-800">
                      {point.iterations}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          {/* RHF Dissociation Note */}
          {hasLargeR && (
            <div
              className="bg-amber-50 border border-amber-200 rounded-lg px-3 py-2 text-xs text-amber-800"
              role="note"
            >
              <p className="font-medium mb-1">
                RHF Dissociation Limit
              </p>
              <p>
                At large bond distances (R &gt; ~3 bohr), RHF produces an incorrect
                dissociation limit. This is a known limitation of single-determinant
                methods -- the energy should approach the sum of isolated atom
                energies but does not. Correlated methods (e.g., CISD, CCSD, FCI)
                are needed for correct dissociation behavior.
              </p>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

// ============================================================================
// Helper: Split points into converged/unconverged
// ============================================================================

interface SplitPoints {
  converged: PesPoint[];
  unconverged: PesPoint[];
}

function splitByConvergence(results: PesPoint[]): SplitPoints {
  const converged: PesPoint[] = [];
  const unconverged: PesPoint[] = [];

  for (const p of results) {
    if (p.converged) {
      converged.push(p);
    } else {
      unconverged.push(p);
    }
  }

  return { converged, unconverged };
}

// ============================================================================
// Main Component
// ============================================================================

/**
 * Props for PesCurvePlot.
 *
 * All new props are optional for backward compatibility (US-083).
 * When `rigidPoints` and `relaxedPoints` are both provided,
 * the plot renders a dual-trace overlay instead of the single-trace
 * store-based mode.
 */
interface PesCurvePlotProps {
  /** Minimum height of the plot in pixels (default: 320) */
  minHeight?: number;
  /** Coordinate type for dynamic x-axis labeling (AC1) */
  coordinateType?: CoordinateType;
  /** Rigid scan points for dual-trace overlay (AC2) */
  rigidPoints?: PesInternalPoint[];
  /** Relaxed scan points for dual-trace overlay (AC2) */
  relaxedPoints?: PesInternalPoint[];
  /**
   * Called when the user hovers over a PES data point.
   * Emits the index into pesInternalResult.points or store results,
   * or null when the user moves away from all points.
   * Used to drive 3D geometry animation from the PES curve.
   */
  onHoverPoint?: (pointIndex: number | null) => void;
  /**
   * Called when the user clicks a PES data point to lock selection.
   * Emits the point index for the clicked point.
   */
  onClickPoint?: (pointIndex: number) => void;
  /**
   * Currently selected (locked) point index.
   * Shown with an amber diamond highlight marker on the plot.
   */
  selectedPointIndex?: number | null;
}

/**
 * Build traces for dual-trace overlay mode (rigid + relaxed).
 *
 * Rigid scan is shown with dashed line and lighter color.
 * Relaxed scan is shown with solid line, darker color, and opt_steps tooltip.
 * Unconverged points in either trace get red open markers.
 *
 * @param rigidPoints - Points from a rigid PES scan
 * @param relaxedPoints - Points from a relaxed PES scan
 * @param coordinateType - Coordinate type for x-axis conversion
 * @param xLabel - X-axis unit label for hover template
 * @returns Plotly data traces
 */
function buildDualTraces(
  rigidPoints: PesInternalPoint[],
  relaxedPoints: PesInternalPoint[],
  coordinateType: CoordinateType,
  xLabel: string,
): Data[] {
  const traces: Data[] = [];

  // --- Rigid trace (dashed, lighter blue) ---
  const rigidConverged = rigidPoints.filter((p) => p.converged);
  const rigidUnconverged = rigidPoints.filter((p) => !p.converged);

  if (rigidConverged.length > 0) {
    traces.push({
      x: convertXValues(rigidConverged.map((p) => p.coordinate_value), coordinateType),
      y: rigidConverged.map((p) => p.energy),
      type: 'scatter',
      mode: 'lines+markers',
      name: 'Rigid',
      line: {
        color: COLOR_RIGID,
        width: 2,
        dash: 'dash',
      },
      marker: {
        color: COLOR_RIGID,
        size: 5,
        symbol: 'circle',
      },
      hovertemplate:
        `x = %{x:.4f} ${xLabel}<br>E = %{y:.10f} Ha<extra>Rigid</extra>`,
    });
  }

  if (rigidUnconverged.length > 0) {
    traces.push({
      x: convertXValues(rigidUnconverged.map((p) => p.coordinate_value), coordinateType),
      y: rigidUnconverged.map((p) => p.energy),
      type: 'scatter',
      mode: 'markers',
      name: 'Rigid (unconverged)',
      marker: {
        color: COLOR_UNCONVERGED,
        size: 7,
        symbol: 'circle-open',
        line: { color: COLOR_UNCONVERGED, width: 2 },
      },
      hovertemplate:
        `x = %{x:.4f} ${xLabel}<br>E = %{y:.10f} Ha<extra>Rigid (unconverged)</extra>`,
      showlegend: false,
    });
  }

  // --- Relaxed trace (solid, darker blue, with opt_steps tooltip) ---
  const relaxedConverged = relaxedPoints.filter((p) => p.converged);
  const relaxedUnconverged = relaxedPoints.filter((p) => !p.converged);

  if (relaxedConverged.length > 0) {
    traces.push({
      x: convertXValues(relaxedConverged.map((p) => p.coordinate_value), coordinateType),
      y: relaxedConverged.map((p) => p.energy),
      type: 'scatter',
      mode: 'lines+markers',
      name: 'Relaxed',
      line: {
        color: COLOR_RELAXED,
        width: 2.5,
      },
      marker: {
        color: COLOR_RELAXED,
        size: 7,
        symbol: 'circle',
      },
      customdata: relaxedConverged.map((p) => [p.opt_steps ?? 'N/A']),
      hovertemplate:
        `x = %{x:.4f} ${xLabel}<br>` +
        'E = %{y:.10f} Ha<br>' +
        'Opt steps: %{customdata[0]}<extra>Relaxed</extra>',
    } as Data);
  }

  if (relaxedUnconverged.length > 0) {
    traces.push({
      x: convertXValues(relaxedUnconverged.map((p) => p.coordinate_value), coordinateType),
      y: relaxedUnconverged.map((p) => p.energy),
      type: 'scatter',
      mode: 'markers',
      name: 'Relaxed (unconverged)',
      marker: {
        color: COLOR_UNCONVERGED,
        size: 9,
        symbol: 'circle-open',
        line: { color: COLOR_UNCONVERGED, width: 2 },
      },
      hovertemplate:
        `x = %{x:.4f} ${xLabel}<br>E = %{y:.10f} Ha<extra>Relaxed (unconverged)</extra>`,
      showlegend: false,
    });
  }

  return traces;
}

/**
 * PES energy vs. coordinate plot with equilibrium marker and notes.
 *
 * Features:
 * - Progressive update as scan points arrive from the worker
 * - Dynamic x-axis label based on coordinate type (bohr/degrees) (AC1)
 * - Dual-trace overlay for rigid vs relaxed comparison (AC2)
 * - Optimization step count in tooltip for relaxed points (AC3)
 * - Solid blue markers + line for converged points
 * - Red hollow markers for unconverged points
 * - Vertical dashed line at equilibrium
 * - Summary panel with equilibrium data and point statistics
 * - Collapsible computational notes with detail table
 * - RHF dissociation limit note for large R
 *
 * Backward compatible: when rigidPoints/relaxedPoints are not provided,
 * falls back to the original single-trace store-based mode.
 *
 * @example
 * ```tsx
 * // Single-trace mode (backward compatible, reads from store)
 * <PesCurvePlot minHeight={350} />
 *
 * // Dynamic axis label
 * <PesCurvePlot coordinateType="angle" />
 *
 * // Dual-trace mode
 * <PesCurvePlot
 *   coordinateType="bond"
 *   rigidPoints={rigidData}
 *   relaxedPoints={relaxedData}
 * />
 * ```
 */
export function PesCurvePlot({
  minHeight = 320,
  coordinateType,
  rigidPoints,
  relaxedPoints,
  onHoverPoint,
  onClickPoint,
  selectedPointIndex,
}: PesCurvePlotProps) {
  const results = useScfStore((state) => state.pesState.results);
  const scanning = useScfStore((state) => state.pesState.scanning);
  const equilibrium = useScfStore((state) => state.pesState.equilibrium);
  const pesInternalResult = useScfStore((state) => state.pesState.pesInternalResult);

  // Determine if we are in dual-trace mode
  const isDualTrace = rigidPoints !== undefined && relaxedPoints !== undefined;

  // Effective x-axis label
  const xAxisLabel = X_AXIS_LABELS[coordinateType ?? 'bond'];
  const xUnitLabel = coordinateType === 'angle' || coordinateType === 'dihedral'
    ? 'deg'
    : 'bohr';

  const hasData = isDualTrace
    ? (rigidPoints.length > 0 || relaxedPoints.length > 0)
    : results.length > 0;

  // Split results into converged and unconverged for single-trace mode
  const { converged, unconverged } = useMemo(
    () => isDualTrace ? { converged: [], unconverged: [] } : splitByConvergence(results),
    [results, isDualTrace]
  );

  // Build Plotly traces
  const plotData = useMemo((): Data[] => {
    if (!hasData) return [];

    // Dynamic coordinate label and unit for hover templates (Fix 2)
    const coordLabel = coordinateType === 'angle'
      ? 'Angle'
      : coordinateType === 'dihedral'
        ? 'Dihedral'
        : 'R';
    const coordUnit = coordinateType === 'angle' || coordinateType === 'dihedral'
      ? 'deg'
      : 'bohr';

    // --- Dual-trace mode (US-083 AC2) ---
    if (isDualTrace) {
      const traces = buildDualTraces(
        rigidPoints,
        relaxedPoints,
        coordinateType ?? 'bond',
        xUnitLabel,
      );

      // Equilibrium marker
      if (equilibrium) {
        const eqX = coordinateType && coordinateType !== 'bond'
          ? equilibrium.r_bohr * (180 / Math.PI)
          : equilibrium.r_bohr;
        traces.push({
          x: [eqX],
          y: [equilibrium.energy_hartree],
          type: 'scatter',
          mode: 'markers',
          name: 'Equilibrium',
          marker: {
            color: COLOR_EQUILIBRIUM,
            size: 12,
            symbol: 'star',
            line: { color: '#16a34a', width: 1 },
          },
          hovertemplate:
            `Eq = %{x:.4f} ${xUnitLabel}<br>E_eq = %{y:.10f} Ha<extra>Equilibrium</extra>`,
          showlegend: true,
        });
      }

      // Selected point marker (amber diamond)
      if (selectedPointIndex != null) {
        const internalPoints = pesInternalResult?.points;
        const selectedPt = internalPoints?.[selectedPointIndex];
        if (selectedPt) {
          const selX = convertXValues([selectedPt.coordinate_value], coordinateType ?? 'bond')[0];
          traces.push({
            x: [selX],
            y: [selectedPt.energy],
            type: 'scatter',
            mode: 'markers',
            name: 'Selected',
            marker: {
              color: COLOR_SELECTED,
              size: 14,
              symbol: 'diamond',
              line: { color: COLOR_SELECTED_BORDER, width: 2 },
            },
            hovertemplate:
              `<b>Selected</b><br>${coordinateType === 'angle' ? 'Angle' : coordinateType === 'dihedral' ? 'Dihedral' : 'R'} = %{x:.4f} ${xUnitLabel}<br>E = %{y:.10f} Ha<extra></extra>`,
            showlegend: false,
          });
        }
      }

      return traces;
    }

    // --- Single-trace mode (backward compatible) ---
    const traces: Data[] = [];

    // Trace 1: Converged points (solid blue line + markers)
    if (converged.length > 0) {
      traces.push({
        x: convertXValues(converged.map((p) => p.r), coordinateType ?? 'bond'),
        y: converged.map((p) => p.energy),
        type: 'scatter',
        mode: 'lines+markers',
        name: 'Converged',
        line: {
          color: COLOR_CONVERGED,
          width: 2,
        },
        marker: {
          color: COLOR_CONVERGED,
          size: 7,
          symbol: 'circle',
        },
        hovertemplate:
          `${coordLabel} = %{x:.4f} ${coordUnit}<br>E = %{y:.10f} Ha<extra>Converged</extra>`,
      });
    }

    // Trace 2: Unconverged points (red open markers, no line)
    if (unconverged.length > 0) {
      traces.push({
        x: convertXValues(unconverged.map((p) => p.r), coordinateType ?? 'bond'),
        y: unconverged.map((p) => p.energy),
        type: 'scatter',
        mode: 'markers',
        name: 'Unconverged',
        marker: {
          color: COLOR_UNCONVERGED,
          size: 9,
          symbol: 'circle-open',
          line: {
            color: COLOR_UNCONVERGED,
            width: 2,
          },
        },
        hovertemplate:
          `${coordLabel} = %{x:.4f} ${coordUnit}<br>E = %{y:.10f} Ha<extra>Unconverged</extra>`,
      });
    }

    // Trace 3: Equilibrium point marker (star)
    if (equilibrium) {
      const eqX = coordinateType && coordinateType !== 'bond'
        ? equilibrium.r_bohr * (180 / Math.PI)
        : equilibrium.r_bohr;
      traces.push({
        x: [eqX],
        y: [equilibrium.energy_hartree],
        type: 'scatter',
        mode: 'markers',
        name: 'Equilibrium',
        marker: {
          color: COLOR_EQUILIBRIUM,
          size: 12,
          symbol: 'star',
          line: {
            color: '#16a34a', // green-600
            width: 1,
          },
        },
        hovertemplate:
          `${coordLabel}_eq = %{x:.4f} ${coordUnit}<br>E_eq = %{y:.10f} Ha<extra>Equilibrium</extra>`,
        showlegend: true,
      });
    }

    // Selected point marker (amber diamond)
    if (selectedPointIndex != null) {
      // Try internal points first (has coordinate_value), fall back to store results
      const internalPoints = pesInternalResult?.points;
      const selectedPt = internalPoints?.[selectedPointIndex];
      if (selectedPt) {
        const selX = convertXValues([selectedPt.coordinate_value], coordinateType ?? 'bond')[0];
        traces.push({
          x: [selX],
          y: [selectedPt.energy],
          type: 'scatter',
          mode: 'markers',
          name: 'Selected',
          marker: {
            color: COLOR_SELECTED,
            size: 14,
            symbol: 'diamond',
            line: { color: COLOR_SELECTED_BORDER, width: 2 },
          },
          hovertemplate:
            `<b>Selected</b><br>${coordLabel} = %{x:.4f} ${coordUnit}<br>E = %{y:.10f} Ha<extra></extra>`,
          showlegend: false,
        });
      } else if (results[selectedPointIndex]) {
        const fallbackPt = results[selectedPointIndex];
        const selX = convertXValues([fallbackPt.r], coordinateType ?? 'bond')[0];
        traces.push({
          x: [selX],
          y: [fallbackPt.energy],
          type: 'scatter',
          mode: 'markers',
          name: 'Selected',
          marker: {
            color: COLOR_SELECTED,
            size: 14,
            symbol: 'diamond',
            line: { color: COLOR_SELECTED_BORDER, width: 2 },
          },
          hovertemplate:
            `<b>Selected</b><br>${coordLabel} = %{x:.4f} ${coordUnit}<br>E = %{y:.10f} Ha<extra></extra>`,
          showlegend: false,
        });
      }
    }

    return traces;
  }, [converged, unconverged, equilibrium, hasData, isDualTrace, rigidPoints, relaxedPoints, coordinateType, xUnitLabel, selectedPointIndex, pesInternalResult, results]);

  // Build layout with equilibrium marker shapes and annotations
  const layout = useMemo((): Partial<Layout> => {
    const baseLayout: Partial<Layout> = {
      title: {
        text: 'Potential Energy Surface',
        font: { size: 14 },
      },
      xaxis: {
        title: { text: xAxisLabel, font: { size: 12 } },
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
        font: { size: 10 },
        bgcolor: 'rgba(255,255,255,0.8)',
        bordercolor: '#e2e8f0',
        borderwidth: 1,
      },
    };

    // Add equilibrium vertical line and annotation
    if (equilibrium) {
      const eqX = coordinateType && coordinateType !== 'bond'
        ? equilibrium.r_bohr * (180 / Math.PI)
        : equilibrium.r_bohr;

      const eqLabel = coordinateType === 'angle' || coordinateType === 'dihedral'
        ? `  ${eqX.toFixed(1)} deg`
        : `  R<sub>eq</sub> = ${eqX.toFixed(3)} bohr`;

      const shapes: Partial<Shape>[] = [
        {
          type: 'line',
          x0: eqX,
          x1: eqX,
          y0: 0,
          y1: 1,
          yref: 'paper',
          line: {
            color: COLOR_EQUILIBRIUM,
            width: 1.5,
            dash: 'dash',
          },
        },
      ];

      const annotations: Partial<Annotations>[] = [
        {
          x: eqX,
          y: 1,
          yref: 'paper',
          xanchor: 'left',
          yanchor: 'top',
          text: eqLabel,
          showarrow: false,
          font: { size: 11, color: '#16a34a' },
          bgcolor: 'rgba(255,255,255,0.85)',
          borderpad: 3,
        },
      ];

      baseLayout.shapes = shapes;
      baseLayout.annotations = annotations;
    }

    return baseLayout;
  }, [equilibrium, xAxisLabel, coordinateType]);

  // Hover handler: map hovered x-value to the closest point index in pesInternalResult.
  // We match by x-value (coordinate_value) rather than Plotly pointIndex because the
  // plot has multiple traces (converged/unconverged/equilibrium) with different indices.
  const handlePlotHover = useCallback((event: PlotHoverEvent) => {
    if (!onHoverPoint || !event.points[0]) return;
    const hoveredX = event.points[0].x;
    if (typeof hoveredX !== 'number') return;

    // Skip hover on the equilibrium marker trace (it's synthetic, not a scan point)
    const traceName = event.points[0].data?.name;
    if (traceName === 'Equilibrium') return;

    // Try to match against pesInternalResult.points (preferred - has geometry data)
    const internalPoints = pesInternalResult?.points;
    if (internalPoints && internalPoints.length > 0) {
      // Convert hovered x back to internal units (display may be in degrees)
      const internalX = (coordinateType === 'angle' || coordinateType === 'dihedral')
        ? hoveredX * (Math.PI / 180)
        : hoveredX;

      let bestIdx = 0;
      let bestDist = Math.abs(internalPoints[0].coordinate_value - internalX);
      for (let i = 1; i < internalPoints.length; i++) {
        const dist = Math.abs(internalPoints[i].coordinate_value - internalX);
        if (dist < bestDist) {
          bestDist = dist;
          bestIdx = i;
        }
      }
      onHoverPoint(bestIdx);
      return;
    }

    // Fallback: match against store results (legacy single-trace mode)
    if (results.length > 0) {
      const internalX = (coordinateType === 'angle' || coordinateType === 'dihedral')
        ? hoveredX * (Math.PI / 180)
        : hoveredX;

      let bestIdx = 0;
      let bestDist = Math.abs(results[0].r - internalX);
      for (let i = 1; i < results.length; i++) {
        const dist = Math.abs(results[i].r - internalX);
        if (dist < bestDist) {
          bestDist = dist;
          bestIdx = i;
        }
      }
      onHoverPoint(bestIdx);
    }
  }, [onHoverPoint, pesInternalResult, results, coordinateType]);

  const handlePlotUnhover = useCallback(() => {
    onHoverPoint?.(null);
  }, [onHoverPoint]);

  // Click handler: map clicked x-value to closest point index (same logic as hover)
  // but calls onClickPoint to lock the selection.
  const handlePlotClick = useCallback((event: PlotMouseEvent) => {
    if (!onClickPoint || !event.points[0]) return;
    const clickedX = event.points[0].x;
    if (typeof clickedX !== 'number') return;

    // Skip click on equilibrium/selected marker traces
    const traceName = (event.points[0] as { data?: { name?: string } }).data?.name;
    if (traceName === 'Equilibrium' || traceName === 'Selected') return;

    // Try internal points first
    const internalPoints = pesInternalResult?.points;
    if (internalPoints && internalPoints.length > 0) {
      const internalX = (coordinateType === 'angle' || coordinateType === 'dihedral')
        ? clickedX * (Math.PI / 180)
        : clickedX;

      let bestIdx = 0;
      let bestDist = Math.abs(internalPoints[0].coordinate_value - internalX);
      for (let i = 1; i < internalPoints.length; i++) {
        const dist = Math.abs(internalPoints[i].coordinate_value - internalX);
        if (dist < bestDist) {
          bestDist = dist;
          bestIdx = i;
        }
      }
      onClickPoint(bestIdx);
      return;
    }

    // Fallback: store results
    if (results.length > 0) {
      const internalX = (coordinateType === 'angle' || coordinateType === 'dihedral')
        ? clickedX * (Math.PI / 180)
        : clickedX;

      let bestIdx = 0;
      let bestDist = Math.abs(results[0].r - internalX);
      for (let i = 1; i < results.length; i++) {
        const dist = Math.abs(results[i].r - internalX);
        if (dist < bestDist) {
          bestDist = dist;
          bestIdx = i;
        }
      }
      onClickPoint(bestIdx);
    }
  }, [onClickPoint, pesInternalResult, results, coordinateType]);

  // Empty state (no data and not scanning)
  if (!hasData && !scanning) {
    return null;
  }

  return (
    <div className="space-y-4">
      {/* Main PES plot */}
      <PlotPanel
        data={plotData}
        layout={layout}
        loading={scanning && results.length === 0}
        minHeight={minHeight}
        ariaLabel={`Potential energy surface plot showing energy versus ${coordinateType ?? 'bond length'}`}
        onHover={onHoverPoint ? handlePlotHover : undefined}
        onUnhover={onHoverPoint ? handlePlotUnhover : undefined}
        onClick={onClickPoint ? handlePlotClick : undefined}
      />

      {/* Summary panel */}
      {!scanning && hasData && <PesSummaryPanel />}

      {/* Computational notes (collapsible) */}
      {!scanning && hasData && <PesNotesPanel />}
    </div>
  );
}

export default PesCurvePlot;
