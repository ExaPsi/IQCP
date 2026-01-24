/**
 * ScfResultDisplay - Container for SCF result visualization.
 *
 * Combines the iteration table, convergence plots, and internals panel
 * into a cohesive display area. Supports two display modes:
 * - "explain": Standard view with plots and summary (default)
 * - "internals": Detailed matrix inspection for US-018
 *
 * @module components/scf/ScfResultDisplay
 */

import { ScfIterationTable } from './ScfIterationTable';
import { ScfEnergyPlot } from './ScfEnergyPlot';
import { ScfResidualPlot } from './ScfResidualPlot';
import { ScfInternalsPanel } from './ScfInternalsPanel';
import { Math } from '../common/Math';
import { useScfStore, type DisplayMode } from '../../stores/scfStore';

/**
 * Props for ScfResultDisplay.
 */
interface ScfResultDisplayProps {
  /** Additional CSS classes */
  className?: string;
}

/**
 * Mode toggle button component.
 */
function ModeToggleButton({
  mode,
  currentMode,
  label,
  icon,
  onClick,
}: {
  mode: DisplayMode;
  currentMode: DisplayMode;
  label: string;
  icon: React.ReactNode;
  onClick: () => void;
}) {
  const isSelected = mode === currentMode;

  return (
    <button
      type="button"
      onClick={onClick}
      className={`flex items-center gap-2 px-4 py-2 text-sm font-medium transition-colors rounded-lg ${
        isSelected
          ? 'bg-blue-600 text-white'
          : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
      }`}
      aria-pressed={isSelected}
      aria-label={`${label} view`}
    >
      {icon}
      {label}
    </button>
  );
}

/**
 * Mode toggle bar component.
 */
function ModeToggle({
  currentMode,
  onModeChange,
  disabled,
}: {
  currentMode: DisplayMode;
  onModeChange: (mode: DisplayMode) => void;
  disabled: boolean;
}) {
  return (
    <div
      className={`flex gap-2 p-1 bg-slate-100 rounded-xl ${
        disabled ? 'opacity-50 pointer-events-none' : ''
      }`}
      role="group"
      aria-label="Result display mode"
    >
      <ModeToggleButton
        mode="explain"
        currentMode={currentMode}
        label="Explanation"
        icon={
          <svg
            className="w-4 h-4"
            fill="none"
            viewBox="0 0 24 24"
            stroke="currentColor"
            aria-hidden="true"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z"
            />
          </svg>
        }
        onClick={() => onModeChange('explain')}
      />
      <ModeToggleButton
        mode="internals"
        currentMode={currentMode}
        label="Inspect Internals"
        icon={
          <svg
            className="w-4 h-4"
            fill="none"
            viewBox="0 0 24 24"
            stroke="currentColor"
            aria-hidden="true"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M4 6a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2H6a2 2 0 01-2-2V6zM14 6a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2h-2a2 2 0 01-2-2V6zM4 16a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2H6a2 2 0 01-2-2v-2zM14 16a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2h-2a2 2 0 01-2-2v-2z"
            />
          </svg>
        }
        onClick={() => onModeChange('internals')}
      />
    </div>
  );
}

/**
 * Standard explanation view with plots and summary.
 */
function ExplanationView({ className = '' }: { className?: string }) {
  const compute = useScfStore((state) => state.compute);

  return (
    <div className={`space-y-4 ${className}`}>
      {/* Iteration Table */}
      <ScfIterationTable maxHeight={280} />

      {/* Plots Grid */}
      <div className="grid grid-cols-1 xl:grid-cols-2 gap-4">
        <ScfEnergyPlot minHeight={280} />
        <ScfResidualPlot minHeight={280} />
      </div>

      {/* Summary Section (when converged) */}
      {compute.status === 'success' && (
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
          <h3 className="text-sm font-semibold text-slate-700 mb-3">
            Computation Summary
          </h3>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
            <div>
              <span className="text-slate-500 block text-xs flex items-center gap-1">
                Final Energy <Math>{String.raw`E`}</Math>
              </span>
              <span className="font-mono text-slate-800">
                {compute.result.energy.toFixed(10)} Ha
              </span>
            </div>
            <div>
              <span className="text-slate-500 block text-xs">Iterations</span>
              <span className="font-mono text-slate-800">
                {compute.result.iterations}
              </span>
            </div>
            <div>
              <span className="text-slate-500 block text-xs">Status</span>
              <span
                className={`font-medium ${
                  compute.result.converged ? 'text-green-700' : 'text-amber-700'
                }`}
              >
                {compute.result.converged ? 'Converged' : 'Did not converge'}
              </span>
            </div>
            <div>
              <span className="text-slate-500 block text-xs">Aborted</span>
              <span className="text-slate-800">
                {compute.result.aborted ? 'Yes' : 'No'}
              </span>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

/**
 * Internals view with matrices and detailed diagnostics.
 */
function InternalsView({ className = '' }: { className?: string }) {
  const compute = useScfStore((state) => state.compute);
  const convergenceProfile = useScfStore((state) => state.convergenceProfile);
  const useDiis = useScfStore((state) => state.useDiis);

  // Internals view requires successful computation with matrix data
  if (compute.status !== 'success') {
    return (
      <div className={`bg-white rounded-xl shadow-sm border border-slate-200 p-8 text-center ${className}`}>
        <div className="text-slate-500">
          <svg
            className="w-12 h-12 mx-auto mb-4 text-slate-300"
            fill="none"
            viewBox="0 0 24 24"
            stroke="currentColor"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={1.5}
              d="M4 6a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2H6a2 2 0 01-2-2V6zM14 6a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2h-2a2 2 0 01-2-2V6zM4 16a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2H6a2 2 0 01-2-2v-2zM14 16a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2h-2a2 2 0 01-2-2v-2z"
            />
          </svg>
          <p className="text-sm font-medium">No internals data available</p>
          <p className="text-xs mt-1">
            Run an SCF calculation to inspect matrices and orbitals
          </p>
        </div>
      </div>
    );
  }

  const { result } = compute;

  // Check if matrix data is available
  if (!result.matrices || !result.orbitalEnergies) {
    return (
      <div className={`bg-amber-50 border border-amber-200 rounded-xl p-6 ${className}`}>
        <h3 className="text-sm font-semibold text-amber-800 mb-2">
          Matrix Data Not Available
        </h3>
        <p className="text-amber-700 text-sm">
          Matrix data was not included in this computation result.
          Re-run the SCF calculation to enable internals inspection.
        </p>
      </div>
    );
  }

  return (
    <ScfInternalsPanel
      matrices={result.matrices}
      orbitalEnergies={result.orbitalEnergies}
      converged={result.converged}
      iterations={result.iterations}
      energy={result.energy}
      convergenceProfile={convergenceProfile}
      diisEnabled={useDiis}
      aborted={result.aborted}
      className={className}
    />
  );
}

/**
 * Container for SCF computation results.
 *
 * Supports two display modes controlled by the store:
 * - "explain": Standard view with iteration table, energy/residual plots
 * - "internals": Matrix heatmaps, orbital energies, and diagnostics
 *
 * @example
 * ```tsx
 * <ScfResultDisplay />
 * ```
 */
export function ScfResultDisplay({ className = '' }: ScfResultDisplayProps) {
  const compute = useScfStore((state) => state.compute);
  const mode = useScfStore((state) => state.mode);
  const setMode = useScfStore((state) => state.setMode);

  // Determine if mode toggle should be enabled
  // Only enable when we have results (not during idle/running/error/cancelled)
  const canToggleMode = compute.status === 'success';

  return (
    <div className={`space-y-4 ${className}`}>
      {/* Mode Toggle Header */}
      <div className="flex items-center justify-between">
        <h2 className="text-lg font-semibold text-slate-800">Results</h2>
        <ModeToggle
          currentMode={mode}
          onModeChange={setMode}
          disabled={!canToggleMode}
        />
      </div>

      {/* Mode-specific content */}
      {mode === 'explain' ? (
        <ExplanationView />
      ) : (
        <InternalsView />
      )}

      {/* Error Display (shown in both modes) */}
      {compute.status === 'error' && (
        <div className="bg-red-50 border border-red-200 rounded-xl p-4">
          <h3 className="text-sm font-semibold text-red-800 mb-1">
            Computation Error
          </h3>
          <p className="text-red-600 text-sm">{compute.error}</p>
        </div>
      )}

      {/* Cancelled Display (shown in both modes) */}
      {compute.status === 'cancelled' && (
        <div className="bg-amber-50 border border-amber-200 rounded-xl p-4">
          <h3 className="text-sm font-semibold text-amber-800 mb-1">
            Computation Cancelled
          </h3>
          <p className="text-amber-600 text-sm">
            The SCF calculation was stopped before completion.
            Partial results (if any) are shown above.
          </p>
        </div>
      )}
    </div>
  );
}

export default ScfResultDisplay;
