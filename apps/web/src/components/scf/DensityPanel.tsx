/**
 * DensityPanel - Controls for density isosurface and cross-section visualization.
 *
 * Progressive disclosure: only appears after SCF convergence and when the
 * density tab is active. Contains controls for:
 * - Density mode selector (total / difference) (US-062, US-063)
 * - Isovalue slider (AC2)
 * - Separate difference density isovalue slider (US-063 AC3)
 * - Plane axis selector (AC5)
 * - Plane position slider (AC6)
 * - Color scale toggle (AC5)
 * - Integrated density display (AC7)
 * - Contextual tooltips (AC4, US-063 AC6)
 *
 * @module components/scf/DensityPanel
 * @see US-062 Density Isosurface & Cross-Sections
 * @see US-063 Difference Density
 */

import { useMemo } from 'react';
import { useScfStore } from '../../stores/scfStore';
import type { PlaneAxis, ColorScale, DensityMode } from '../../types/density';
import type { DifferenceDensityResult } from '../../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/**
 * Props for the DensityPanel component.
 */
export interface DensityPanelProps {
  /** Whether density grid data is available */
  hasGrid: boolean;
  /** Whether density grid is being computed */
  loading: boolean;
  /** Error message from density computation */
  error: string | null;
  /** Integrated density from the grid */
  integratedDensity: number | null;
  /** Expected number of electrons */
  nElectronsExpected: number | null;
  /** Maximum density value in the grid */
  maxDensity: number | null;
  /** Grid bounds for plane position slider */
  gridBounds: {
    min: [number, number, number];
    max: [number, number, number];
  } | null;
  /** Difference density grid result (US-063) */
  diffGridData: DifferenceDensityResult | null;
}

// ============================================================================
// Sub-components
// ============================================================================

/**
 * Info icon SVG for contextual tooltips.
 */
function InfoIcon() {
  return (
    <svg
      className="w-4 h-4 text-slate-400"
      viewBox="0 0 20 20"
      fill="currentColor"
      aria-hidden="true"
    >
      <path
        fillRule="evenodd"
        d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a.75.75 0 000 1.5h.253a.25.25 0 01.244.304l-.459 2.066A1.75 1.75 0 0010.747 15H11a.75.75 0 000-1.5h-.253a.25.25 0 01-.244-.304l.459-2.066A1.75 1.75 0 009.253 9H9z"
        clipRule="evenodd"
      />
    </svg>
  );
}

/**
 * Integrated density display with quality indicator (AC7).
 */
function IntegratedDensityDisplay({
  integratedDensity,
  nElectronsExpected,
}: {
  integratedDensity: number | null;
  nElectronsExpected: number | null;
}) {
  const { deviation, quality, qualityColor } = useMemo(() => {
    if (integratedDensity === null || nElectronsExpected === null || nElectronsExpected === 0) {
      return { deviation: null, quality: null, qualityColor: '' };
    }

    const dev = ((integratedDensity - nElectronsExpected) / nElectronsExpected) * 100;
    const absDev = Math.abs(dev);

    let q: string;
    let color: string;
    if (absDev < 3) {
      q = 'Good';
      color = 'text-green-600';
    } else if (absDev < 5) {
      q = 'Acceptable';
      color = 'text-yellow-600';
    } else {
      q = 'Low';
      color = 'text-red-600';
    }

    return { deviation: dev, quality: q, qualityColor: color };
  }, [integratedDensity, nElectronsExpected]);

  if (integratedDensity === null || nElectronsExpected === null) {
    return null;
  }

  return (
    <div className="bg-slate-50 rounded-lg p-3 space-y-1.5">
      <div className="flex items-center gap-1">
        <h4 className="text-xs font-semibold text-slate-600">Integrated Density</h4>
        <span
          className="cursor-help"
          title="This is a numerical integration over a finite grid. Deviations from N_electrons are due to grid discretization, not SCF errors."
          aria-label="Integrated density explanation"
        >
          <InfoIcon />
        </span>
      </div>
      <div className="grid grid-cols-2 gap-x-3 gap-y-0.5 text-xs">
        <span className="text-slate-500">N_integrated:</span>
        <span className="text-slate-800 font-mono">{integratedDensity.toFixed(2)}</span>
        <span className="text-slate-500">N_electrons:</span>
        <span className="text-slate-800 font-mono">{nElectronsExpected}</span>
        {deviation !== null && (
          <>
            <span className="text-slate-500">Deviation:</span>
            <span className={`font-mono ${qualityColor}`}>
              {deviation >= 0 ? '+' : ''}{deviation.toFixed(1)}%
            </span>
          </>
        )}
        {quality && (
          <>
            <span className="text-slate-500">Grid quality:</span>
            <span className={`font-medium ${qualityColor}`}>{quality}</span>
          </>
        )}
      </div>
    </div>
  );
}

/**
 * Integrated Delta-rho display for difference density mode (US-063 AC5).
 */
function IntegratedDeltaDisplay({
  diffGridData,
}: {
  diffGridData: DifferenceDensityResult | null;
}) {
  if (!diffGridData) return null;

  return (
    <div className="bg-slate-50 rounded-lg p-3 space-y-1.5">
      <div className="flex items-center gap-1">
        <h4 className="text-xs font-semibold text-slate-600">Density Conservation</h4>
        <span
          className="cursor-help"
          title="The integrated difference density should be approximately zero, since the molecule has the same total number of electrons as the sum of its constituent atoms. Small deviations are due to grid discretization."
          aria-label="Density conservation explanation"
        >
          <InfoIcon />
        </span>
      </div>
      <div className="grid grid-cols-2 gap-x-3 gap-y-0.5 text-xs">
        <span className="text-slate-500">Integrated Delta-rho:</span>
        <span className="text-slate-800 font-mono">{diffGridData.integratedDeltaRho.toFixed(4)}</span>
        <span className="text-slate-500">Max accumulation:</span>
        <span className="text-green-700 font-mono">+{diffGridData.maxAccumulation.toFixed(4)}</span>
        <span className="text-slate-500">Max depletion:</span>
        <span className="text-red-700 font-mono">{diffGridData.maxDepletion.toFixed(4)}</span>
      </div>
    </div>
  );
}

// ============================================================================
// Main Component
// ============================================================================

/**
 * DensityPanel - density visualization controls.
 *
 * Renders controls for the density isosurface and cross-section visualization.
 * Designed for the sidebar area of Module E, appearing only after SCF
 * convergence (progressive disclosure).
 *
 * @example
 * ```tsx
 * <DensityPanel
 *   hasGrid={densityViewer.gridData !== null}
 *   loading={densityViewer.loading}
 *   error={densityViewer.error}
 *   integratedDensity={densityViewer.gridData?.integratedDensity ?? null}
 *   nElectronsExpected={densityViewer.gridData?.nElectronsExpected ?? null}
 *   maxDensity={densityViewer.gridData?.maxDensity ?? null}
 *   gridBounds={gridBounds}
 *   diffGridData={densityViewer.diffGridData}
 * />
 * ```
 */
export function DensityPanel({
  hasGrid,
  loading,
  error,
  integratedDensity,
  nElectronsExpected,
  // maxDensity reserved for future adaptive isovalue range
  gridBounds,
  diffGridData,
}: DensityPanelProps) {
  // Read density state from store
  const densityIsovalue = useScfStore((s) => s.viewerState.densityIsovalue);
  const diffIsovalue = useScfStore((s) => s.viewerState.diffIsovalue);
  const densityMode = useScfStore((s) => s.viewerState.densityMode);
  const planeAxis = useScfStore((s) => s.viewerState.planeAxis);
  const planePosition = useScfStore((s) => s.viewerState.planePosition);
  const colorScale = useScfStore((s) => s.viewerState.colorScale);
  const showDensity = useScfStore((s) => s.viewerState.showDensity);

  // Store actions
  const setDensityIsovalue = useScfStore((s) => s.setDensityIsovalue);
  const setDiffIsovalue = useScfStore((s) => s.setDiffIsovalue);
  const setDensityMode = useScfStore((s) => s.setDensityMode);
  const setPlaneAxis = useScfStore((s) => s.setPlaneAxis);
  const setPlanePosition = useScfStore((s) => s.setPlanePosition);
  const setColorScale = useScfStore((s) => s.setColorScale);
  const setShowDensity = useScfStore((s) => s.setShowDensity);

  const isDifference = densityMode === 'difference';

  // Compute perpendicular axis range for plane position slider
  const perpRange = useMemo(() => {
    if (!gridBounds) return { min: -5, max: 5, step: 0.2 };
    const { min: gMin, max: gMax } = gridBounds;

    switch (planeAxis) {
      case 'xy':
        return { min: gMin[2], max: gMax[2], step: 0.2 };
      case 'xz':
        return { min: gMin[1], max: gMax[1], step: 0.2 };
      case 'yz':
        return { min: gMin[0], max: gMax[0], step: 0.2 };
    }
  }, [gridBounds, planeAxis]);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      {/* Header */}
      <div className="px-4 py-3 bg-slate-50 border-b border-slate-200 flex items-center justify-between">
        <h3 className="text-sm font-semibold text-slate-700">
          Density Visualization
        </h3>
        <label className="flex items-center gap-1.5 text-xs text-slate-600 cursor-pointer select-none">
          <input
            type="checkbox"
            checked={showDensity}
            onChange={(e) => setShowDensity(e.target.checked)}
            className="rounded border-slate-300 text-teal-600 focus:ring-teal-500 h-3.5 w-3.5"
          />
          Show
        </label>
      </div>

      <div className="p-4 space-y-4">
        {/* Loading indicator */}
        {loading && (
          <div className="flex items-center gap-2 text-xs text-slate-500">
            <div className="animate-spin rounded-full h-3 w-3 border-b-2 border-teal-600" />
            {isDifference ? 'Computing difference density...' : 'Computing density grid...'}
          </div>
        )}

        {/* Error display */}
        {error && (
          <div className="bg-red-50 border border-red-200 rounded-lg p-2">
            <p className="text-xs text-red-600">{error}</p>
          </div>
        )}

        {/* Density mode selector (US-063) */}
        <div className="space-y-1.5">
          <div className="flex items-center gap-1">
            <label className="text-xs font-medium text-slate-700">
              Density Mode
            </label>
            {isDifference && (
              <span
                className="cursor-help"
                title="The difference density shows how electrons rearrange when atoms form a molecule. The 'reference' is a promolecule -- the sum of spherically-averaged isolated atom densities placed at the molecular geometry. Regions where Delta-rho > 0 (solid green) show where electrons have accumulated relative to separated atoms; regions where Delta-rho < 0 (translucent red) show depletion. This is not an absolute physical quantity -- it depends on the choice of atomic reference."
                aria-label="Difference density explanation"
              >
                <InfoIcon />
              </span>
            )}
          </div>
          <div className="flex gap-1">
            {(['total', 'difference'] as DensityMode[]).map((mode) => (
              <button
                key={mode}
                onClick={() => setDensityMode(mode)}
                disabled={!showDensity || loading}
                className={`flex-1 px-2 py-1 text-xs font-medium rounded transition-colors ${
                  densityMode === mode
                    ? 'bg-teal-600 text-white'
                    : 'bg-slate-100 text-slate-600 hover:bg-slate-200'
                } disabled:opacity-50 disabled:cursor-not-allowed`}
                aria-pressed={densityMode === mode}
              >
                {mode === 'total' ? 'Total' : 'Difference'}
              </button>
            ))}
          </div>
        </div>

        {/* Element availability note for difference mode (US-063 AC7) */}
        {isDifference && (
          <p className="text-[10px] text-amber-600 bg-amber-50 rounded px-2 py-1">
            Difference density is available for molecules containing H-Ar only.
          </p>
        )}

        {/* Total density isovalue slider (AC2) -- shown in total mode */}
        {!isDifference && (
          <div className="space-y-1.5">
            <div className="flex items-center gap-1">
              <label
                htmlFor="density-isovalue"
                className="text-xs font-medium text-slate-700"
              >
                Density Isovalue
              </label>
              <span
                className="cursor-help"
                title="This surface encloses all points where the electron density exceeds the chosen threshold."
                aria-label="Density isosurface explanation"
              >
                <InfoIcon />
              </span>
            </div>
            <div className="flex items-center gap-2">
              <input
                id="density-isovalue"
                type="range"
                min={0.01}
                max={0.50}
                step={0.005}
                value={densityIsovalue}
                onChange={(e) => setDensityIsovalue(parseFloat(e.target.value))}
                disabled={!showDensity || loading}
                className="flex-1 h-1.5 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-teal-600 disabled:opacity-50"
                aria-label="Density isovalue slider"
              />
              <span className="text-xs font-mono text-slate-600 w-12 text-right">
                {densityIsovalue.toFixed(3)}
              </span>
            </div>
            <p className="text-[10px] text-slate-400">e/bohr^3</p>
          </div>
        )}

        {/* Difference density isovalue slider (US-063 AC3) -- shown in difference mode */}
        {isDifference && (
          <div className="space-y-1.5">
            <div className="flex items-center gap-1">
              <label
                htmlFor="diff-isovalue"
                className="text-xs font-medium text-slate-700"
              >
                Difference Isovalue
              </label>
              <span
                className="cursor-help"
                title="Controls the threshold for accumulation (green, positive) and depletion (red, negative) isosurfaces. Lower values reveal more subtle charge rearrangement; higher values show only the strongest effects."
                aria-label="Difference density isovalue explanation"
              >
                <InfoIcon />
              </span>
            </div>
            <div className="flex items-center gap-2">
              <input
                id="diff-isovalue"
                type="range"
                min={0.001}
                max={0.05}
                step={0.001}
                value={diffIsovalue}
                onChange={(e) => setDiffIsovalue(parseFloat(e.target.value))}
                disabled={!showDensity || loading}
                className="flex-1 h-1.5 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-teal-600 disabled:opacity-50"
                aria-label="Difference density isovalue slider"
              />
              <span className="text-xs font-mono text-slate-600 w-14 text-right">
                {diffIsovalue.toFixed(3)}
              </span>
            </div>
            <p className="text-[10px] text-slate-400">e/bohr^3</p>
          </div>
        )}

        {/* Plane axis selector (AC5) */}
        <div className="space-y-1.5">
          <label className="text-xs font-medium text-slate-700">
            Cross-Section Plane
          </label>
          <div className="flex gap-1">
            {(['xy', 'xz', 'yz'] as PlaneAxis[]).map((axis) => (
              <button
                key={axis}
                onClick={() => setPlaneAxis(axis)}
                disabled={!hasGrid}
                className={`flex-1 px-2 py-1 text-xs font-medium rounded transition-colors ${
                  planeAxis === axis
                    ? 'bg-teal-600 text-white'
                    : 'bg-slate-100 text-slate-600 hover:bg-slate-200'
                } disabled:opacity-50 disabled:cursor-not-allowed`}
                aria-pressed={planeAxis === axis}
              >
                {axis.toUpperCase()}
              </button>
            ))}
          </div>
        </div>

        {/* Plane position slider (AC6) */}
        <div className="space-y-1.5">
          <label
            htmlFor="plane-position"
            className="text-xs font-medium text-slate-700"
          >
            Plane Position
          </label>
          <div className="flex items-center gap-2">
            <input
              id="plane-position"
              type="range"
              min={perpRange.min}
              max={perpRange.max}
              step={perpRange.step}
              value={planePosition}
              onChange={(e) => setPlanePosition(parseFloat(e.target.value))}
              disabled={!hasGrid}
              className="flex-1 h-1.5 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-teal-600 disabled:opacity-50"
              aria-label="Plane position slider"
            />
            <span className="text-xs font-mono text-slate-600 w-16 text-right">
              {planePosition.toFixed(2)} bohr
            </span>
          </div>
        </div>

        {/* Color scale toggle (AC5) */}
        <div className="space-y-1.5">
          <label className="text-xs font-medium text-slate-700">
            Color Scale
          </label>
          <div className="flex gap-1">
            {(['linear', 'log'] as ColorScale[]).map((scale) => (
              <button
                key={scale}
                onClick={() => setColorScale(scale)}
                disabled={!hasGrid}
                className={`flex-1 px-2 py-1 text-xs font-medium rounded transition-colors ${
                  colorScale === scale
                    ? 'bg-teal-600 text-white'
                    : 'bg-slate-100 text-slate-600 hover:bg-slate-200'
                } disabled:opacity-50 disabled:cursor-not-allowed`}
                aria-pressed={colorScale === scale}
              >
                {scale === 'linear' ? 'Linear' : 'Log'}
              </button>
            ))}
          </div>
        </div>

        {/* Integrated density / Delta-rho display */}
        {isDifference ? (
          <IntegratedDeltaDisplay diffGridData={diffGridData} />
        ) : (
          <IntegratedDensityDisplay
            integratedDensity={integratedDensity}
            nElectronsExpected={nElectronsExpected}
          />
        )}
      </div>
    </div>
  );
}
