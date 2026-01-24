/**
 * RysControlsPanel - Parameter input controls for Rys quadrature.
 *
 * Provides inputs for n (order), T (parameter), target accuracy,
 * and display mode. Connected to the rysStore for state management.
 *
 * @module components/rys/RysControlsPanel
 */

import { useRysStore, type TargetAccuracy } from '../../stores/rysStore';
import { Math } from '../common/Math';

/**
 * Props for RysControlsPanel.
 */
interface RysControlsPanelProps {
  /** Disable all controls (e.g., during computation) */
  disabled?: boolean;
}

/**
 * Target accuracy options with display labels.
 */
const TARGET_OPTIONS: { value: TargetAccuracy; latex: string }[] = [
  { value: '1e-4', latex: '10^{-4}' },
  { value: '1e-6', latex: '10^{-6}' },
  { value: '1e-8', latex: '10^{-8}' },
];

/**
 * Shell quartet options for quick selection.
 * Each defines: label, total angular momentum L, and computed n_r.
 * n_r = floor(L/2) + 1
 */
interface ShellQuartet {
  label: string;
  L: number;
  nRoots: number;
  description: string;
}

const SHELL_QUARTETS: ShellQuartet[] = [
  { label: '(ss|ss)', L: 0, nRoots: 1, description: 's-type only' },
  { label: '(ps|ss)', L: 1, nRoots: 1, description: 'p and s orbitals' },
  { label: '(pp|ss)', L: 2, nRoots: 2, description: 'p-p overlap' },
  { label: '(pp|pp)', L: 4, nRoots: 3, description: 'p-p on both' },
  { label: '(dd|pp)', L: 6, nRoots: 4, description: 'd orbitals' },
  { label: '(dd|dd)', L: 8, nRoots: 5, description: 'd-d overlap' },
  { label: '(ff|ff)', L: 12, nRoots: 7, description: 'f orbitals' },
];

/**
 * Controls panel for Rys quadrature parameters.
 *
 * Includes:
 * - Order (n): Integer selector 1-10 with +/- buttons
 * - Parameter (T): Range slider 0-50
 * - Target accuracy: Segmented buttons (10^-4, 10^-6, 10^-8)
 * - Display mode: Toggle between 'explain' and 'internals'
 *
 * @example
 * ```tsx
 * <RysControlsPanel />
 * ```
 */
export function RysControlsPanel({ disabled = false }: RysControlsPanelProps) {
  const n = useRysStore((state) => state.n);
  const T = useRysStore((state) => state.T);
  const target = useRysStore((state) => state.target);
  const mode = useRysStore((state) => state.mode);
  const setN = useRysStore((state) => state.setN);
  const setT = useRysStore((state) => state.setT);
  const setTarget = useRysStore((state) => state.setTarget);
  const setMode = useRysStore((state) => state.setMode);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
      <h2 className="text-lg font-semibold text-slate-800 mb-4">Parameters</h2>

      {/* Order (n) input */}
      <div className="mb-6">
        <label
          htmlFor="rys-n-input"
          className="block text-sm font-medium text-slate-700 mb-2"
        >
          Order (n)
        </label>
        <div className="flex items-center gap-3">
          <button
            type="button"
            onClick={() => setN(n - 1)}
            disabled={disabled || n <= 1}
            className="w-10 h-10 flex items-center justify-center rounded-lg border border-slate-300 bg-white text-slate-700 hover:bg-slate-50 disabled:opacity-50 disabled:cursor-not-allowed focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500"
            aria-label="Decrease n"
          >
            <span className="text-xl">-</span>
          </button>
          <input
            id="rys-n-input"
            type="number"
            min={1}
            max={10}
            step={1}
            value={n}
            onChange={(e) => setN(Number(e.target.value))}
            disabled={disabled}
            className="w-20 h-10 text-center text-lg font-mono border border-slate-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50"
            aria-label="Order n value"
          />
          <button
            type="button"
            onClick={() => setN(n + 1)}
            disabled={disabled || n >= 10}
            className="w-10 h-10 flex items-center justify-center rounded-lg border border-slate-300 bg-white text-slate-700 hover:bg-slate-50 disabled:opacity-50 disabled:cursor-not-allowed focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500"
            aria-label="Increase n"
          >
            <span className="text-xl">+</span>
          </button>
        </div>
        <p className="mt-1 text-xs text-slate-500">
          Quadrature points (1 to 10)
        </p>
      </div>

      {/* Shell Quartet Quick Select */}
      <div className="mb-6">
        <label
          htmlFor="rys-shell-quartet"
          className="block text-sm font-medium text-slate-700 mb-2"
        >
          Quick Select: Shell Quartet
        </label>
        <select
          id="rys-shell-quartet"
          value=""
          onChange={(e) => {
            const quartet = SHELL_QUARTETS.find(q => q.label === e.target.value);
            if (quartet && quartet.nRoots <= 10) {
              setN(quartet.nRoots);
            }
          }}
          disabled={disabled}
          className="w-full h-10 px-3 border border-slate-300 rounded-lg bg-white text-slate-700 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
        >
          <option value="">Select shell quartet...</option>
          {SHELL_QUARTETS.map((quartet) => (
            <option
              key={quartet.label}
              value={quartet.label}
              disabled={quartet.nRoots > 10}
            >
              {quartet.label} (L={quartet.L}, n={quartet.nRoots})
            </option>
          ))}
        </select>
        <p className="mt-1 text-xs text-slate-500">
          Sets <Math>{'n'}</Math> based on root count rule: <Math>{'n_r = \\lfloor L/2 \\rfloor + 1'}</Math>
        </p>
      </div>

      {/* Parameter (T) slider */}
      <div className="mb-6">
        <label
          htmlFor="rys-t-slider"
          className="block text-sm font-medium text-slate-700 mb-2"
        >
          Parameter <Math>{'T'}</Math>
        </label>
        <input
          id="rys-t-slider"
          type="range"
          min={0}
          max={50}
          step={0.1}
          value={T}
          onChange={(e) => setT(Number(e.target.value))}
          disabled={disabled}
          className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600 disabled:opacity-50 disabled:cursor-not-allowed"
          aria-label="Parameter T slider"
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

      {/* Target accuracy selector */}
      <div className="mb-6">
        <label className="block text-sm font-medium text-slate-700 mb-2">
          Target Accuracy
        </label>
        <div className="flex gap-2">
          {TARGET_OPTIONS.map((option) => (
            <button
              key={option.value}
              type="button"
              onClick={() => setTarget(option.value)}
              disabled={disabled}
              className={`flex-1 px-3 py-2 rounded-lg text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
                target === option.value
                  ? 'bg-blue-600 text-white'
                  : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
              } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
              aria-pressed={target === option.value}
            >
              <Math>{option.latex}</Math>
            </button>
          ))}
        </div>
        <p className="mt-1 text-xs text-slate-500">
          Maximum acceptable reconstruction error
        </p>
      </div>

      {/* Display mode toggle */}
      <div>
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
            ? 'Educational explanations of quadrature'
            : 'Raw computational details and sanity checks'}
        </p>
      </div>
    </div>
  );
}

export default RysControlsPanel;
