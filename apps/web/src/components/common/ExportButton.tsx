/**
 * ExportButton - Button to export computation results as JSON artifact.
 *
 * Used for downloading run artifacts containing state, results, and metadata
 * for sharing, grading, and reproducibility.
 * Shows "Exported!" feedback after successful export.
 *
 * @module components/common/ExportButton
 */

import { useState, useCallback } from 'react';

/**
 * Props for ExportButton component.
 */
export interface ExportButtonProps {
  /** Function to create the artifact when clicked */
  onExport: () => void;
  /** Whether the button should be disabled (no results available) */
  disabled?: boolean;
  /** Additional CSS classes */
  className?: string;
  /** Button label (default: "Export Results") */
  label?: string;
}

/**
 * Duration to show "Exported!" feedback (ms).
 */
const FEEDBACK_DURATION_MS = 2000;

/**
 * ExportButton - Downloads artifact as JSON with visual feedback.
 *
 * Features:
 * - Download icon (SVG)
 * - Shows "Exported!" feedback for 2 seconds after click
 * - Disabled state when no results available
 * - Uses same button style as CopyLinkButton for consistency
 * - Accessible with proper ARIA attributes
 *
 * @example
 * ```tsx
 * // Basic usage
 * <ExportButton
 *   onExport={handleExport}
 *   disabled={!hasResults}
 * />
 *
 * // With custom label
 * <ExportButton
 *   onExport={handleExport}
 *   label="Download JSON"
 *   disabled={isLoading}
 * />
 * ```
 */
export function ExportButton({
  onExport,
  disabled = false,
  className = '',
  label = 'Export Results',
}: ExportButtonProps) {
  const [exported, setExported] = useState(false);
  const [error, setError] = useState(false);

  const handleExport = useCallback(() => {
    try {
      onExport();
      setExported(true);
      setError(false);

      // Reset after feedback duration
      setTimeout(() => {
        setExported(false);
      }, FEEDBACK_DURATION_MS);
    } catch (err) {
      console.error('Failed to export artifact:', err);
      setError(true);

      // Reset error state after feedback duration
      setTimeout(() => {
        setError(false);
      }, FEEDBACK_DURATION_MS);
    }
  }, [onExport]);

  // Determine button text and style based on state
  const buttonText = exported ? 'Exported!' : error ? 'Failed' : label;
  const buttonStyle = exported
    ? 'bg-green-600 hover:bg-green-700 text-white'
    : error
      ? 'bg-red-600 hover:bg-red-700 text-white'
      : 'bg-slate-100 hover:bg-slate-200 text-slate-700';

  return (
    <button
      type="button"
      onClick={handleExport}
      disabled={disabled || exported}
      className={`
        inline-flex items-center gap-2 px-4 py-2 rounded-lg
        text-sm font-medium transition-colors
        focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500
        disabled:opacity-50 disabled:cursor-not-allowed
        ${buttonStyle}
        ${className}
      `}
      aria-label={exported ? 'Results exported successfully' : 'Export results as JSON file'}
    >
      {/* Icon */}
      {exported ? (
        // Checkmark icon
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
            d="M5 13l4 4L19 7"
          />
        </svg>
      ) : error ? (
        // X icon
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
            d="M6 18L18 6M6 6l12 12"
          />
        </svg>
      ) : (
        // Download icon
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
            d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4"
          />
        </svg>
      )}
      {buttonText}
    </button>
  );
}

export default ExportButton;
