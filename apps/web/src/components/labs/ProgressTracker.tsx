/**
 * ProgressTracker - Visual indicator of worksheet completion progress.
 *
 * Displays a progress bar and percentage, with options to view
 * detailed progress or clear all progress.
 *
 * @module components/labs/ProgressTracker
 */

import { useState, useCallback } from 'react';

/**
 * Props for ProgressTracker component.
 */
export interface ProgressTrackerProps {
  /** Completion percentage (0-100) */
  percent: number;
  /** Number of completed sections */
  completedCount: number;
  /** Total number of trackable sections */
  totalCount: number;
  /** Callback to clear all progress */
  onClearProgress?: () => void;
  /** Additional CSS classes */
  className?: string;
}

/**
 * ProgressTracker - Shows worksheet completion progress.
 *
 * Features:
 * - Visual progress bar
 * - Percentage and section count display
 * - Optional clear progress button with confirmation
 * - Accessible progress indicator
 *
 * @example
 * ```tsx
 * <ProgressTracker
 *   percent={60}
 *   completedCount={3}
 *   totalCount={5}
 *   onClearProgress={() => clearProgress()}
 * />
 * ```
 */
export function ProgressTracker({
  percent,
  completedCount,
  totalCount,
  onClearProgress,
  className = '',
}: ProgressTrackerProps) {
  const [showConfirm, setShowConfirm] = useState(false);

  const handleClearClick = useCallback(() => {
    setShowConfirm(true);
  }, []);

  const handleConfirmClear = useCallback(() => {
    onClearProgress?.();
    setShowConfirm(false);
  }, [onClearProgress]);

  const handleCancelClear = useCallback(() => {
    setShowConfirm(false);
  }, []);

  // Determine progress bar color based on completion
  const progressColor =
    percent === 100
      ? 'bg-green-500'
      : percent >= 50
        ? 'bg-blue-500'
        : 'bg-blue-400';

  return (
    <div className={`${className}`}>
      {/* Progress bar */}
      <div className="mb-2">
        <div className="flex items-center justify-between text-sm mb-1">
          <span className="font-medium text-slate-700">Progress</span>
          <span className="text-slate-600">
            {completedCount} / {totalCount} sections
          </span>
        </div>
        <div
          className="h-2 bg-slate-200 rounded-full overflow-hidden"
          role="progressbar"
          aria-valuenow={percent}
          aria-valuemin={0}
          aria-valuemax={100}
          aria-label={`Worksheet progress: ${percent}%`}
        >
          <div
            className={`h-full ${progressColor} transition-all duration-300`}
            style={{ width: `${percent}%` }}
          />
        </div>
      </div>

      {/* Percentage and clear button */}
      <div className="flex items-center justify-between">
        <span className="text-sm text-slate-600">
          {percent}% complete
          {percent === 100 && (
            <span className="ml-2 text-green-600 font-medium">
              Lab complete!
            </span>
          )}
        </span>

        {onClearProgress && completedCount > 0 && !showConfirm && (
          <button
            type="button"
            onClick={handleClearClick}
            className="
              text-xs text-slate-500 hover:text-slate-700
              underline transition-colors
            "
          >
            Clear progress
          </button>
        )}

        {showConfirm && (
          <div className="flex items-center gap-2">
            <span className="text-xs text-slate-500">Clear all progress?</span>
            <button
              type="button"
              onClick={handleConfirmClear}
              className="
                px-2 py-0.5 text-xs rounded
                bg-red-100 text-red-700 hover:bg-red-200
                transition-colors
              "
            >
              Yes
            </button>
            <button
              type="button"
              onClick={handleCancelClear}
              className="
                px-2 py-0.5 text-xs rounded
                bg-slate-100 text-slate-700 hover:bg-slate-200
                transition-colors
              "
            >
              No
            </button>
          </div>
        )}
      </div>
    </div>
  );
}

export default ProgressTracker;
