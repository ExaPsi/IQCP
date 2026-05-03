/**
 * FockStepView - Matrix display for each step of the Fock build process.
 *
 * Renders the appropriate matrix heatmap for the current step:
 *   Step 1: H^core matrix with title and educational text
 *   Step 2: G, J, or K matrix (selected by step2View toggle)
 *   Step 3: F = H^core + G matrix with clickable elements
 *
 * Uses a general-purpose Plotly heatmap (not the IntegralMatrixHeatmap,
 * which is coupled to MatrixTab = S|T|V|Hcore). The Fock tracer needs
 * different colorscales and titles for J, K, G, F, density matrices.
 *
 * @module components/integrals/FockStepView
 * @see US-058 Fock Build Tracing
 */

import { useMemo, useCallback, useRef, useEffect, memo } from 'react';
import Plotly from 'plotly.js-cartesian-dist';
import createPlotlyComponent from 'react-plotly.js/factory';
import type { Data, Layout, Config, PlotMouseEvent } from 'plotly.js';
import type {
  FockDecompositionResult,
  FockBuildStep,
  FockStep2View,
} from '../../stores/integralStore';
import { Math as MathDisplay } from '../common/Math';

const Plot = createPlotlyComponent(Plotly);

// ============================================================================
// Types
// ============================================================================

export interface FockStepViewProps {
  /** Fock decomposition result data */
  decomposition: FockDecompositionResult;
  /** Current step (1, 2, or 3) */
  step: FockBuildStep;
  /** Sub-view within step 2 (G, J, or K) */
  step2View: FockStep2View;
  /** Selected element [row, col] for highlighting */
  selectedElement: [number, number] | null;
  /** Element selection handler (only active in step 3) */
  onElementSelect: (element: [number, number] | null) => void;
}

// ============================================================================
// Configuration
// ============================================================================

interface MatrixViewConfig {
  title: string;
  colorscale: string;
  diverging: boolean;
  educationalText: string;
  formula?: React.ReactNode;
}

function getStepConfig(
  step: FockBuildStep,
  step2View: FockStep2View,
): MatrixViewConfig {
  switch (step) {
    case 1:
      return {
        title: 'Step 1: Core Hamiltonian H<sup>core</sup> = T + V',
        colorscale: 'RdBu',
        diverging: true,
        educationalText:
          'The core Hamiltonian represents the one-electron kinetic energy (T) and nuclear attraction (V), computed before electron-electron interactions are included. This is the same H\u1D9C\u1D52\u02B3\u1D49 matrix from the one-electron integral inspector.',
      };
    case 2:
      switch (step2View) {
        case 'G':
          return {
            title: 'Step 2: Two-Electron Contribution G(P)',
            colorscale: 'RdBu',
            diverging: true,
            formula: <MathDisplay>{'G = J - 0.5 \\cdot K'}</MathDisplay>,
            educationalText:
              'G combines Coulomb repulsion (J) and exchange stabilization (K). Coulomb repulsion is a classical effect: electrons repel each other. Exchange is a purely quantum effect arising from the Pauli exclusion principle. The factor of 0.5 in G = J - 0.5*K accounts for double-counting in the RHF density matrix convention (P = 2*C*C^T).',
          };
        case 'J':
          return {
            title: 'Coulomb Matrix J',
            colorscale: 'Blues',
            diverging: false,
            formula: <MathDisplay>{'J_{\\mu\\nu} = \\sum_{\\lambda\\sigma} P_{\\lambda\\sigma} (\\mu\\nu|\\lambda\\sigma)'}</MathDisplay>,
            educationalText:
              'The Coulomb matrix represents classical electrostatic repulsion between the charge distribution of basis functions mu,nu and the total electron density. J is always positive on the diagonal: an electron always repels the electron density described by the same basis function.',
          };
        case 'K':
          return {
            title: 'Exchange Matrix K',
            colorscale: 'Oranges',
            diverging: false,
            formula: <MathDisplay>{'K_{\\mu\\nu} = \\sum_{\\lambda\\sigma} P_{\\lambda\\sigma} (\\mu\\lambda|\\nu\\sigma)'}</MathDisplay>,
            educationalText:
              'The exchange matrix has no classical analogue. It arises from the antisymmetry requirement of the many-electron wavefunction (Pauli exclusion principle). Exchange reduces electron repulsion for same-spin electrons, stabilizing the system. Notice the different ERI index ordering: (ml|ns) vs. Coulomb (mn|ls).',
          };
      }
      break; // unreachable but TypeScript needs it
    case 3:
      return {
        title: 'Step 3: Fock Matrix F = H<sup>core</sup> + G',
        colorscale: 'RdBu',
        diverging: true,
        formula: <MathDisplay>{'F = H^{\\text{core}} + G = H^{\\text{core}} + J - 0.5K'}</MathDisplay>,
        educationalText:
          'The Fock matrix is the effective one-electron Hamiltonian in Hartree-Fock theory. It combines the fixed one-electron terms (H\u1D9C\u1D52\u02B3\u1D49) with the mean-field electron-electron interaction (G). Click any element to see the numerical breakdown.',
      };
  }
  // Fallback (should not be reached)
  return {
    title: 'Matrix',
    colorscale: 'RdBu',
    diverging: true,
    educationalText: '',
  };
}

// ============================================================================
// Helpers
// ============================================================================

/**
 * Get the matrix data for the current step and sub-view.
 */
function getMatrixForStep(
  decomposition: FockDecompositionResult,
  step: FockBuildStep,
  step2View: FockStep2View,
): number[] {
  switch (step) {
    case 1:
      return decomposition.hCore;
    case 2:
      switch (step2View) {
        case 'G':
          return decomposition.gMatrix;
        case 'J':
          return decomposition.jMatrix;
        case 'K':
          return decomposition.kMatrix;
      }
      break; // unreachable
    case 3:
      return decomposition.fMatrix;
  }
  return decomposition.fMatrix; // fallback
}

/**
 * Reshape a flat row-major array into a 2D matrix.
 */
function reshapeToMatrix(data: number[], nbf: number): number[][] {
  const matrix: number[][] = [];
  for (let i = 0; i < nbf; i++) {
    const row: number[] = [];
    for (let j = 0; j < nbf; j++) {
      row.push(data[i * nbf + j]);
    }
    matrix.push(row);
  }
  return matrix;
}

// ============================================================================
// Component
// ============================================================================

/**
 * Renders the matrix heatmap for the current Fock build step.
 *
 * Features:
 * - Step-appropriate colorscale and title
 * - Educational text and formula display
 * - Selected element highlighting (step 3 only)
 * - Click interaction for element drill-down (step 3 only)
 * - Hover template with full-precision values
 * - Keyboard Escape to deselect
 */
export const FockStepView = memo(function FockStepView({
  decomposition,
  step,
  step2View,
  selectedElement,
  onElementSelect,
}: FockStepViewProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const { nbf, labels } = decomposition;
  const config = getStepConfig(step, step2View);
  const matrixData = getMatrixForStep(decomposition, step, step2View);

  // Reshape flat array to 2D
  const matrix2D = useMemo(() => reshapeToMatrix(matrixData, nbf), [matrixData, nbf]);

  // Color configuration
  const colorConfig = useMemo(() => {
    if (config.diverging) {
      const absMax = Math.max(...matrixData.map(Math.abs));
      return { zmid: 0, zmin: -absMax, zmax: absMax };
    }
    return {};
  }, [matrixData, config.diverging]);

  // Plotly data traces
  const plotData = useMemo((): Data[] => {
    const traces: Data[] = [
      {
        z: matrix2D,
        x: labels,
        y: labels,
        type: 'heatmap',
        colorscale: config.colorscale,
        ...colorConfig,
        hovertemplate:
          '%{y} | %{x}<br>Value: %{z:.8e}<extra></extra>',
        showscale: true,
        colorbar: {
          title: { text: 'Value', font: { size: 11 } },
          thickness: 15,
          len: 0.9,
          tickformat: '.3e',
          tickfont: { size: 10 },
        },
      } as Data,
    ];

    // Selection highlight (step 3 only)
    if (step === 3 && selectedElement) {
      const [row, col] = selectedElement;
      traces.push({
        x: [labels[col]],
        y: [labels[row]],
        type: 'scatter',
        mode: 'markers',
        marker: {
          size: 20,
          color: 'rgba(0,0,0,0)',
          line: {
            color: '#000000',
            width: 3,
          },
          symbol: 'square',
        },
        hoverinfo: 'skip',
        showlegend: false,
      } as Data);
    }

    return traces;
  }, [matrix2D, labels, config.colorscale, colorConfig, step, selectedElement]);

  // Plotly layout
  const layout = useMemo((): Partial<Layout> => {
    return {
      title: {
        text: config.title,
        font: { size: 14 },
      },
      xaxis: {
        title: { text: 'Basis Function', font: { size: 11 } },
        tickmode: 'array' as const,
        tickvals: labels,
        ticktext: labels,
        tickfont: { size: nbf > 10 ? 8 : 10 },
        tickangle: nbf > 7 ? -45 : 0,
        side: 'bottom' as const,
        constrain: 'domain' as const,
      },
      yaxis: {
        title: { text: 'Basis Function', font: { size: 11 } },
        tickmode: 'array' as const,
        tickvals: labels,
        ticktext: labels,
        tickfont: { size: nbf > 10 ? 8 : 10 },
        autorange: 'reversed' as const,
        scaleanchor: 'x',
        constrain: 'domain' as const,
      },
      margin: { t: 50, r: 90, b: nbf > 7 ? 80 : 60, l: nbf > 7 ? 80 : 60 },
      autosize: true,
      paper_bgcolor: 'transparent',
      plot_bgcolor: 'white',
      hovermode: 'closest' as const,
    };
  }, [config.title, labels, nbf]);

  // Plotly config
  const plotConfig = useMemo((): Partial<Config> => ({
    responsive: true,
    displayModeBar: 'hover' as const,
    modeBarButtonsToRemove: [
      'lasso2d',
      'select2d',
      'autoScale2d',
      'hoverClosestCartesian',
      'hoverCompareCartesian',
    ],
    displaylogo: false,
    toImageButtonOptions: {
      format: 'png' as const,
      filename: `iqcp-fock-step${step}`,
      scale: 2,
    },
  }), [step]);

  // Click handler (only active in step 3)
  const handleClick = useCallback(
    (event: PlotMouseEvent) => {
      if (step !== 3) return;
      if (!event.points || event.points.length === 0) return;
      const point = event.points[0];

      const colLabel = point.x as string;
      const rowLabel = point.y as string;
      const col = labels.indexOf(colLabel);
      const row = labels.indexOf(rowLabel);

      if (row >= 0 && col >= 0) {
        onElementSelect([row, col]);
      }
    },
    [step, labels, onElementSelect],
  );

  // Keyboard: Escape to deselect
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Escape' && selectedElement) {
        e.preventDefault();
        onElementSelect(null);
      }
    };

    container.addEventListener('keydown', handleKeyDown);
    return () => container.removeEventListener('keydown', handleKeyDown);
  }, [selectedElement, onElementSelect]);

  // Build a label for the current step/sub-view (used for SVG filename)
  const stepLabel = step === 2 ? `step2-${step2View}` : `step${step}`;
  const plotDivId = `fock-step-heatmap-${stepLabel}`;

  return (
    <div className="space-y-3">
      {/* Formula display */}
      {config.formula && (
        <div className="bg-slate-50 rounded-lg px-4 py-2 border border-slate-200">
          <p className="text-sm font-mono text-slate-800 text-center">
            {config.formula}
          </p>
        </div>
      )}

      {/* Heatmap */}
      <div
        ref={containerRef}
        className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden"
        style={{ minHeight: '400px' }}
        role="img"
        aria-label={`${config.title} heatmap with ${nbf} basis functions`}
        tabIndex={0}
      >
        <Plot
          divId={plotDivId}
          data={plotData}
          layout={layout}
          config={plotConfig}
          onClick={handleClick}
          useResizeHandler={true}
          style={{ width: '100%', height: '100%' }}
        />
      </div>

      {/* Export SVG + click instruction row */}
      <div className="flex items-center justify-between">
        <div>
          {step === 3 && !selectedElement && (
            <p className="text-xs text-slate-500">
              Click any element in the Fock matrix to see its H<sup>core</sup> + J − 0.5K breakdown.
            </p>
          )}
        </div>
        <button
          onClick={() => {
            const el = document.getElementById(plotDivId);
            if (el && (window as unknown as Record<string, unknown>).Plotly) {
              (
                (window as unknown as Record<string, unknown>).Plotly as {
                  downloadImage: (
                    el: HTMLElement,
                    opts: Record<string, unknown>,
                  ) => void;
                }
              ).downloadImage(el, {
                format: 'svg',
                filename: `iqcp-fock-${stepLabel}-matrix`,
                width: 600,
                height: 500,
              });
            }
          }}
          className="text-xs text-slate-500 hover:text-slate-700 flex items-center gap-1"
          title="Export heatmap as SVG"
        >
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
          </svg>
          Export SVG
        </button>
      </div>

      {/* Educational text */}
      <div className="bg-blue-50 rounded-lg px-4 py-3 border border-blue-100">
        <p className="text-xs text-blue-800 leading-relaxed">
          {config.educationalText}
        </p>
      </div>
    </div>
  );
});
