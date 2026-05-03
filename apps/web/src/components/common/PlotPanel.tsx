/**
 * PlotPanel - Responsive Plotly.js wrapper component.
 *
 * Wraps react-plotly.js with:
 * - Responsive sizing using ResizeObserver
 * - Loading overlay while data is being computed
 * - Error state display
 * - Consistent styling with the IQCP design system
 *
 * @module components/common/PlotPanel
 */

import { useRef, useEffect, useState, useMemo } from 'react';
// Use cartesian distribution to include heatmap support
// (basic-dist only has scatter, bar, pie)
import Plotly from 'plotly.js-cartesian-dist';
import createPlotlyComponent from 'react-plotly.js/factory';
import type { Data, Layout, Config, PlotHoverEvent, PlotMouseEvent } from 'plotly.js';

// Create the Plot component using the cartesian Plotly distribution
const Plot = createPlotlyComponent(Plotly);

/**
 * Props for PlotPanel component.
 */
export interface PlotPanelProps {
  /** Plotly data traces */
  data: Data[];
  /** Plotly layout configuration */
  layout?: Partial<Layout>;
  /** Plotly config options */
  config?: Partial<Config>;
  /** Show loading overlay when true */
  loading?: boolean;
  /** Error message to display (overrides chart) */
  error?: string | null;
  /** Additional CSS classes for the container */
  className?: string;
  /** Accessible label for the plot */
  ariaLabel?: string;
  /** Minimum height in pixels (default: 300) */
  minHeight?: number;
  /** Called when a data point is hovered */
  onHover?: (event: PlotHoverEvent) => void;
  /** Called when a data point is unhovered */
  onUnhover?: (event: PlotMouseEvent) => void;
  /** Called when user clicks on plot data */
  onClick?: (event: PlotMouseEvent) => void;
  /** Optional HTML id for the plot div (used by Plotly.downloadImage) */
  divId?: string;
}

/**
 * Default Plotly config for IQCP plots.
 *
 * Optimized for educational use:
 * - Mode bar appears on hover
 * - Removed less useful tools
 * - Responsive sizing enabled
 */
const DEFAULT_CONFIG: Partial<Config> = {
  responsive: true,
  displayModeBar: 'hover',
  modeBarButtonsToRemove: [
    'lasso2d',
    'select2d',
    'autoScale2d',
    'hoverClosestCartesian',
    'hoverCompareCartesian',
  ],
  displaylogo: false,
  toImageButtonOptions: {
    format: 'png',
    filename: 'iqcp-plot',
    scale: 2,
  },
};

/**
 * Default Plotly layout for IQCP plots.
 *
 * Clean, scientific appearance with:
 * - Minimal margins
 * - Clean font choices
 * - Automatic sizing
 */
const DEFAULT_LAYOUT: Partial<Layout> = {
  autosize: true,
  margin: { t: 40, r: 30, b: 50, l: 60 },
  font: { family: 'system-ui, -apple-system, sans-serif', size: 12 },
  paper_bgcolor: 'transparent',
  plot_bgcolor: 'white',
  hovermode: 'closest',
};

/**
 * Loading spinner component for overlay.
 */
function LoadingOverlay() {
  return (
    <div className="absolute inset-0 bg-white/80 flex items-center justify-center z-10">
      <div className="flex items-center gap-3">
        <div className="animate-spin rounded-full h-6 w-6 border-b-2 border-blue-600" />
        <span className="text-slate-600 text-sm">Computing...</span>
      </div>
    </div>
  );
}

/**
 * Error display component.
 */
function ErrorDisplay({ message }: { message: string }) {
  return (
    <div className="absolute inset-0 flex items-center justify-center">
      <div className="bg-red-50 border border-red-200 rounded-lg p-4 max-w-md text-center">
        <p className="text-red-800 font-medium">Plot Error</p>
        <p className="text-red-600 text-sm mt-1">{message}</p>
      </div>
    </div>
  );
}

/**
 * PlotPanel - Responsive Plotly wrapper with loading and error states.
 *
 * Features:
 * - Uses ResizeObserver for responsive sizing
 * - Shows loading overlay during computation
 * - Displays error state when data fails to load
 * - Consistent styling with IQCP design system
 *
 * @example
 * ```tsx
 * <PlotPanel
 *   data={[{ x: [1, 2, 3], y: [1, 4, 9], type: 'scatter' }]}
 *   layout={{ title: 'My Plot' }}
 *   loading={isComputing}
 *   error={computeError}
 *   ariaLabel="Boys function curve"
 * />
 * ```
 */
export function PlotPanel({
  data,
  layout,
  config,
  loading = false,
  error = null,
  className = '',
  ariaLabel = 'Interactive plot',
  minHeight = 300,
  onHover,
  onUnhover,
  onClick,
  divId,
}: PlotPanelProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [dimensions, setDimensions] = useState({ width: 0, height: 0 });

  // Set up ResizeObserver for responsive sizing
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const resizeObserver = new ResizeObserver((entries) => {
      const entry = entries[0];
      if (entry) {
        const { width, height } = entry.contentRect;
        setDimensions({
          width: Math.floor(width),
          height: Math.floor(Math.max(height, minHeight)),
        });
      }
    });

    resizeObserver.observe(container);

    // Initial measurement
    const rect = container.getBoundingClientRect();
    setDimensions({
      width: Math.floor(rect.width),
      height: Math.floor(Math.max(rect.height, minHeight)),
    });

    return () => resizeObserver.disconnect();
  }, [minHeight]);

  // Merge layouts with responsive dimensions
  const mergedLayout = useMemo((): Partial<Layout> => {
    return {
      ...DEFAULT_LAYOUT,
      ...layout,
      width: dimensions.width > 0 ? dimensions.width : undefined,
      height: dimensions.height > 0 ? dimensions.height : undefined,
    };
  }, [layout, dimensions]);

  // Merge configs
  const mergedConfig = useMemo((): Partial<Config> => {
    return {
      ...DEFAULT_CONFIG,
      ...config,
    };
  }, [config]);

  return (
    <div
      ref={containerRef}
      className={`relative bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden ${className}`}
      style={{ minHeight: `${minHeight}px` }}
      role="img"
      aria-label={ariaLabel}
    >
      {/* Loading overlay */}
      {loading && <LoadingOverlay />}

      {/* Error display */}
      {error && <ErrorDisplay message={error} />}

      {/* Plot (hidden when error) */}
      {!error && dimensions.width > 0 && (
        <Plot
          data={data}
          layout={mergedLayout}
          config={mergedConfig}
          useResizeHandler={false} // We handle resizing ourselves
          style={{ width: '100%', height: '100%' }}
          onHover={onHover}
          onUnhover={onUnhover}
          onClick={onClick}
          divId={divId}
        />
      )}

      {/* Placeholder when dimensions not yet measured */}
      {!error && dimensions.width === 0 && (
        <div className="absolute inset-0 flex items-center justify-center">
          <div className="animate-pulse text-slate-400 text-sm">
            Initializing plot...
          </div>
        </div>
      )}
    </div>
  );
}

export default PlotPanel;
