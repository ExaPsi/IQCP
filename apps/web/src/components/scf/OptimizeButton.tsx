/**
 * OptimizeButton - Trigger geometry optimization with progress bar.
 *
 * Shows "Optimize Geometry" when idle, a progress bar with step count,
 * energy, and gradient when running, and a cancel button during optimization.
 *
 * @module components/scf/OptimizeButton
 * @see US-075 Optimization UI + Trajectory Animation (AC1)
 */

import { useScfStore } from '../../stores/scfStore';

/**
 * Props for OptimizeButton.
 */
export interface OptimizeButtonProps {
  /** Trigger optimization */
  onOptimize: () => void;
  /** Cancel running optimization */
  onCancel: () => void;
  /** Whether optimization is running */
  loading: boolean;
  /** Whether the button should be disabled */
  disabled: boolean;
  /** Maximum optimization steps (default: 50) */
  maxSteps?: number;
}

/**
 * Button that triggers geometry optimization with live progress bar.
 *
 * Displays a multi-phase progress bar during optimization:
 * - Step 0 (initial evaluation): indeterminate bar
 * - Step 1+: determinate bar with energy and gradient readout
 */
export function OptimizeButton({
  onOptimize,
  onCancel,
  loading,
  disabled,
  maxSteps = 50,
}: OptimizeButtonProps) {
  const progress = useScfStore(
    (state) => state.optimizationState.progress
  );
  const stepCount = progress.length;
  const lastStep = stepCount > 0 ? progress[stepCount - 1] : null;
  const stepProgress = stepCount > 0 ? (stepCount / maxSteps) * 100 : 0;

  return (
    <div className="space-y-2">
      {/* Button row */}
      <div className="flex gap-2">
        <button
          type="button"
          onClick={loading ? undefined : onOptimize}
          disabled={disabled || loading}
          className={`flex-1 px-4 py-2.5 rounded-lg text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-amber-500 ${
            loading
              ? 'bg-amber-400 text-white cursor-not-allowed'
              : 'bg-amber-600 text-white hover:bg-amber-700'
          } disabled:opacity-50 disabled:cursor-not-allowed`}
          aria-label={
            loading
              ? `Optimizing geometry, step ${stepCount}`
              : 'Optimize geometry'
          }
        >
          {loading ? (
            <span className="flex items-center justify-center gap-2">
              <span className="animate-spin rounded-full h-3.5 w-3.5 border-b-2 border-white" />
              Optimizing...
            </span>
          ) : (
            <span className="flex items-center justify-center gap-2">
              <svg
                className="w-4 h-4"
                fill="none"
                stroke="currentColor"
                viewBox="0 0 24 24"
                aria-hidden="true"
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M13 7h8m0 0v8m0-8l-8 8-4-4-6 6"
                />
              </svg>
              Optimize Geometry
            </span>
          )}
        </button>
        {loading && (
          <button
            type="button"
            onClick={onCancel}
            className="px-3 py-2.5 rounded-lg text-sm font-medium transition-colors bg-slate-100 text-slate-700 hover:bg-slate-200 focus:outline-none focus-visible:ring-2 focus-visible:ring-slate-500"
            aria-label="Cancel optimization"
          >
            Cancel
          </button>
        )}
      </div>

      {/* Progress bar (during optimization) */}
      {loading && (
        <div className="space-y-1.5">
          {/* Phase label */}
          <div className="flex items-center justify-between text-xs">
            <span className="text-amber-700 font-medium flex items-center gap-1.5">
              <span className="w-1.5 h-1.5 rounded-full bg-amber-600 animate-pulse" />
              {stepCount === 0
                ? 'Evaluating initial geometry...'
                : `Step ${stepCount}/${maxSteps}`}
            </span>
            {stepCount > 0 && (
              <span className="text-slate-500 font-mono">{Math.min(stepProgress, 100).toFixed(0)}%</span>
            )}
          </div>

          {/* Progress bar */}
          <div className="h-1.5 bg-slate-200 rounded-full overflow-hidden">
            {stepCount === 0 ? (
              <div className="h-full bg-amber-500 rounded-full animate-progress-indeterminate" />
            ) : (
              <div
                className="h-full bg-amber-500 rounded-full transition-all duration-300"
                style={{ width: `${Math.min(stepProgress, 100)}%` }}
              />
            )}
          </div>

          {/* Live energy & gradient readout */}
          {lastStep && (
            <p className="text-xs text-slate-500 font-mono">
              E = {lastStep.energy.toFixed(10)} Ha
              <span className="ml-2">
                |g|<sub>max</sub> = {lastStep.maxGradient.toExponential(2)}
              </span>
            </p>
          )}
        </div>
      )}
    </div>
  );
}
