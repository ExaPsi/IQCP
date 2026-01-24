/**
 * IntegralComputeButton - Button to trigger integral computation.
 *
 * Provides:
 * - "Compute Integrals" button (disabled when geometry invalid)
 * - Progress display during computation (phase + percentage)
 * - Cancel button during computation
 *
 * @module components/scf/custom/IntegralComputeButton
 */

import type { IntegralProgress } from '../../../worker/protocol';

/**
 * Props for IntegralComputeButton component.
 */
export interface IntegralComputeButtonProps {
  /** Callback to start computation */
  onCompute: () => void;
  /** Callback to cancel computation */
  onCancel: () => void;
  /** Disable the button (e.g., geometry invalid) */
  disabled?: boolean;
  /** Whether computation is currently running */
  isComputing: boolean;
  /** Current progress (when computing) */
  progress?: IntegralProgress;
}

/**
 * Human-readable labels for integral computation phases.
 */
const PHASE_LABELS: Record<string, string> = {
  overlap: 'Computing overlap integrals (S)',
  kinetic: 'Computing kinetic integrals (T)',
  nuclear: 'Computing nuclear attraction (V)',
  eri: 'Computing two-electron integrals (ERI)',
  assembly: 'Assembling integral matrices',
};

/**
 * Get human-readable label for a phase.
 */
function getPhaseLabel(phase: string): string {
  return PHASE_LABELS[phase] || `Computing ${phase}`;
}

/**
 * Progress bar component.
 */
function ProgressBar({ percent }: { percent: number }) {
  const clampedPercent = Math.max(0, Math.min(100, percent));

  return (
    <div
      className="w-full h-2 bg-slate-200 rounded-full overflow-hidden"
      role="progressbar"
      aria-valuenow={clampedPercent}
      aria-valuemin={0}
      aria-valuemax={100}
    >
      <div
        className="h-full bg-blue-600 transition-all duration-200"
        style={{ width: `${clampedPercent}%` }}
      />
    </div>
  );
}

/**
 * Button to compute integrals with progress display.
 *
 * Shows a progress bar and phase description during computation,
 * with an option to cancel the operation.
 *
 * @example
 * ```tsx
 * const [isComputing, setIsComputing] = useState(false);
 * const [progress, setProgress] = useState<IntegralProgress>();
 *
 * <IntegralComputeButton
 *   onCompute={handleCompute}
 *   onCancel={handleCancel}
 *   disabled={!isGeometryValid}
 *   isComputing={isComputing}
 *   progress={progress}
 * />
 * ```
 */
export function IntegralComputeButton({
  onCompute,
  onCancel,
  disabled = false,
  isComputing,
  progress,
}: IntegralComputeButtonProps) {
  if (isComputing) {
    return (
      <div>
        {/* Progress Display */}
        <div className="mb-3">
          <div className="flex items-center justify-between mb-1">
            <span className="text-sm font-medium text-slate-700">
              {progress ? getPhaseLabel(progress.phase) : 'Starting...'}
            </span>
            <span className="text-sm text-slate-500">
              {progress ? `${Math.round(progress.overallPercent)}%` : '0%'}
            </span>
          </div>
          <ProgressBar percent={progress?.overallPercent || 0} />
          {progress?.message && (
            <p className="mt-1 text-xs text-slate-500">{progress.message}</p>
          )}
        </div>

        {/* Cancel Button */}
        <button
          type="button"
          onClick={onCancel}
          className="w-full px-4 py-3 rounded-lg font-medium transition-colors bg-red-100 text-red-700 hover:bg-red-200 focus:outline-none focus-visible:ring-2 focus-visible:ring-red-500"
          aria-label="Cancel integral computation"
        >
          <span className="flex items-center justify-center gap-2">
            <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
              <path
                fillRule="evenodd"
                d="M10 18a8 8 0 100-16 8 8 0 000 16zM8 7a1 1 0 00-1 1v4a1 1 0 001 1h4a1 1 0 001-1V8a1 1 0 00-1-1H8z"
                clipRule="evenodd"
              />
            </svg>
            Cancel Computation
          </span>
        </button>
      </div>
    );
  }

  return (
    <div>
      <button
        type="button"
        onClick={onCompute}
        disabled={disabled}
        className={`w-full px-4 py-3 rounded-lg font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
          disabled
            ? 'bg-slate-200 text-slate-400 cursor-not-allowed'
            : 'bg-green-600 text-white hover:bg-green-700'
        }`}
        aria-label="Compute integrals for custom geometry"
      >
        <span className="flex items-center justify-center gap-2">
          <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M9 7h6m0 10v-3m-3 3h.01M9 17h.01M9 14h.01M12 14h.01M15 11h.01M12 11h.01M9 11h.01M7 21h10a2 2 0 002-2V5a2 2 0 00-2-2H7a2 2 0 00-2 2v14a2 2 0 002 2z"
            />
          </svg>
          Compute Integrals
        </span>
      </button>
      {disabled && (
        <p className="mt-1 text-xs text-slate-500 text-center">
          Enter a valid geometry to enable integral computation
        </p>
      )}
    </div>
  );
}

export default IntegralComputeButton;
