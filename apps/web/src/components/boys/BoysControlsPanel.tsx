/**
 * BoysControlsPanel - Parameter input controls for Boys function.
 *
 * Provides inputs for m (order), T (argument), and view mode.
 * Connected to the boysStore for state management.
 *
 * @module components/boys/BoysControlsPanel
 */

import { useBoysStore } from '../../stores/boysStore';

/**
 * Props for BoysControlsPanel.
 */
interface BoysControlsPanelProps {
  /** Disable all controls (e.g., during computation) */
  disabled?: boolean;
}

/**
 * Controls panel for Boys function parameters.
 *
 * Includes:
 * - Order (m): Integer selector 0-10 with +/- buttons
 * - Argument (T): Range slider 0-50
 * - View mode: Toggle between 'single' and 'sweep'
 *
 * @example
 * ```tsx
 * <BoysControlsPanel />
 * ```
 */
export function BoysControlsPanel({ disabled = false }: BoysControlsPanelProps) {
  const m = useBoysStore((state) => state.m);
  const T = useBoysStore((state) => state.T);
  const view = useBoysStore((state) => state.view);
  const mode = useBoysStore((state) => state.mode);
  const logScale = useBoysStore((state) => state.logScale);
  const setM = useBoysStore((state) => state.setM);
  const setT = useBoysStore((state) => state.setT);
  const setView = useBoysStore((state) => state.setView);
  const setMode = useBoysStore((state) => state.setMode);
  const setLogScale = useBoysStore((state) => state.setLogScale);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
      <h2 className="text-lg font-semibold text-slate-800 mb-4">Parameters</h2>

      {/* Order (m) input */}
      <div className="mb-6">
        <label
          htmlFor="boys-m-input"
          className="block text-sm font-medium text-slate-700 mb-2"
        >
          Order (m)
        </label>
        <div className="flex items-center gap-3">
          <button
            type="button"
            onClick={() => setM(m - 1)}
            disabled={disabled || m <= 0}
            className="w-10 h-10 flex items-center justify-center rounded-lg border border-slate-300 bg-white text-slate-700 hover:bg-slate-50 disabled:opacity-50 disabled:cursor-not-allowed focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500"
            aria-label="Decrease m"
          >
            <span className="text-xl">−</span>
          </button>
          <input
            id="boys-m-input"
            type="number"
            min={0}
            max={10}
            step={1}
            value={m}
            onChange={(e) => setM(Number(e.target.value))}
            disabled={disabled}
            className="w-20 h-10 text-center text-lg font-mono border border-slate-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50"
            aria-label="Order m value"
          />
          <button
            type="button"
            onClick={() => setM(m + 1)}
            disabled={disabled || m >= 10}
            className="w-10 h-10 flex items-center justify-center rounded-lg border border-slate-300 bg-white text-slate-700 hover:bg-slate-50 disabled:opacity-50 disabled:cursor-not-allowed focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500"
            aria-label="Increase m"
          >
            <span className="text-xl">+</span>
          </button>
        </div>
        <p className="mt-1 text-xs text-slate-500">
          Integer from 0 to 10
        </p>
      </div>

      {/* Argument (T) slider */}
      <div className="mb-6">
        <label
          htmlFor="boys-t-slider"
          className="block text-sm font-medium text-slate-700 mb-2"
        >
          Argument (T)
        </label>
        <input
          id="boys-t-slider"
          type="range"
          min={0}
          max={50}
          step={0.1}
          value={T}
          onChange={(e) => setT(Number(e.target.value))}
          disabled={disabled}
          className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600 disabled:opacity-50 disabled:cursor-not-allowed"
          aria-label="Argument T slider"
          aria-valuemin={0}
          aria-valuemax={50}
          aria-valuenow={T}
        />
        <div className="flex justify-between mt-2">
          <span className="text-xs text-slate-500">0</span>
          <span className="text-sm font-mono text-slate-700">
            T = {T.toFixed(1)}
          </span>
          <span className="text-xs text-slate-500">50</span>
        </div>
      </div>

      {/* View mode toggle */}
      <div className="mb-6">
        <label className="block text-sm font-medium text-slate-700 mb-2">
          View Mode
        </label>
        <div className="flex gap-2">
          <button
            type="button"
            onClick={() => setView('single')}
            disabled={disabled}
            className={`flex-1 px-4 py-2 rounded-lg text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
              view === 'single'
                ? 'bg-blue-600 text-white'
                : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
            } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
            aria-pressed={view === 'single'}
          >
            Single
          </button>
          <button
            type="button"
            onClick={() => setView('sweep')}
            disabled={disabled}
            className={`flex-1 px-4 py-2 rounded-lg text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
              view === 'sweep'
                ? 'bg-blue-600 text-white'
                : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
            } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
            aria-pressed={view === 'sweep'}
          >
            Sweep
          </button>
        </div>
        <p className="mt-1 text-xs text-slate-500">
          {view === 'single'
            ? 'Evaluate at single T value'
            : 'Plot F_m(T) over T range'}
        </p>
      </div>

      {/* Display mode toggle */}
      <div className="mb-6">
        <label className="block text-sm font-medium text-slate-700 mb-2">
          Display Mode
        </label>
        <div className="flex gap-2">
          <button
            type="button"
            onClick={() => setMode('explain')}
            disabled={disabled}
            className={`flex-1 px-4 py-2 rounded-lg text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
              mode === 'explain'
                ? 'bg-blue-600 text-white'
                : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
            } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
            aria-pressed={mode === 'explain'}
          >
            Explain
          </button>
          <button
            type="button"
            onClick={() => setMode('internals')}
            disabled={disabled}
            className={`flex-1 px-4 py-2 rounded-lg text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-purple-500 ${
              mode === 'internals'
                ? 'bg-purple-600 text-white'
                : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
            } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
            aria-pressed={mode === 'internals'}
          >
            Internals
          </button>
        </div>
        <p className="mt-1 text-xs text-slate-500">
          {mode === 'explain'
            ? 'Educational explanations of computation methods'
            : 'Raw computational details and metrics'}
        </p>
      </div>

      {/* Log scale toggle (only shown in sweep view) */}
      {view === 'sweep' && (
        <div>
          <label className="flex items-center gap-3 cursor-pointer">
            <input
              type="checkbox"
              checked={logScale}
              onChange={(e) => setLogScale(e.target.checked)}
              disabled={disabled}
              className="w-4 h-4 rounded border-slate-300 text-blue-600 focus:ring-2 focus:ring-blue-500 disabled:opacity-50"
              aria-label="Toggle logarithmic y-axis scale"
            />
            <span className="text-sm font-medium text-slate-700">
              Log scale (y-axis)
            </span>
          </label>
          <p className="mt-1 text-xs text-slate-500 ml-7">
            Use logarithmic scale for better visualization of small values
          </p>
        </div>
      )}
    </div>
  );
}

export default BoysControlsPanel;
