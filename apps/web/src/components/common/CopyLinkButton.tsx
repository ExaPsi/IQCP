/**
 * CopyLinkButton - Button to copy the current URL to clipboard.
 *
 * Used for sharing deep links to specific application states.
 * Shows "Copied!" feedback after successful copy.
 *
 * @module components/common/CopyLinkButton
 */

import { useState, useCallback } from 'react';

/**
 * Props for CopyLinkButton component.
 */
export interface CopyLinkButtonProps {
  /** Optional custom URL to copy (defaults to current URL) */
  url?: string;
  /** Additional CSS classes */
  className?: string;
  /** Button label (default: "Copy Link") */
  label?: string;
  /** Disabled state */
  disabled?: boolean;
}

/**
 * Duration to show "Copied!" feedback (ms).
 */
const FEEDBACK_DURATION_MS = 2000;

/**
 * CopyLinkButton - Copies URL to clipboard with visual feedback.
 *
 * Features:
 * - Uses navigator.clipboard API for secure copying
 * - Shows "Copied!" feedback for 2 seconds after successful copy
 * - Accessible with proper aria labels
 * - Falls back gracefully if clipboard API unavailable
 *
 * @example
 * ```tsx
 * // Copy current page URL
 * <CopyLinkButton />
 *
 * // Copy specific URL with custom label
 * <CopyLinkButton
 *   url="https://example.com/state=xyz"
 *   label="Share State"
 * />
 * ```
 */
export function CopyLinkButton({
  url,
  className = '',
  label = 'Copy Link',
  disabled = false,
}: CopyLinkButtonProps) {
  const [copied, setCopied] = useState(false);
  const [error, setError] = useState(false);

  const handleCopy = useCallback(async () => {
    const textToCopy = url ?? window.location.href;

    try {
      await navigator.clipboard.writeText(textToCopy);
      setCopied(true);
      setError(false);

      // Reset after feedback duration
      setTimeout(() => {
        setCopied(false);
      }, FEEDBACK_DURATION_MS);
    } catch (err) {
      console.error('Failed to copy to clipboard:', err);
      setError(true);

      // Reset error state after feedback duration
      setTimeout(() => {
        setError(false);
      }, FEEDBACK_DURATION_MS);
    }
  }, [url]);

  // Determine button text and style based on state
  const buttonText = copied ? 'Copied!' : error ? 'Failed' : label;
  const buttonStyle = copied
    ? 'bg-green-600 hover:bg-green-700 text-white'
    : error
      ? 'bg-red-600 hover:bg-red-700 text-white'
      : 'bg-slate-100 hover:bg-slate-200 text-slate-700';

  return (
    <button
      type="button"
      onClick={handleCopy}
      disabled={disabled || copied}
      className={`
        inline-flex items-center gap-2 px-4 py-2 rounded-lg
        text-sm font-medium transition-colors
        focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500
        disabled:opacity-50 disabled:cursor-not-allowed
        ${buttonStyle}
        ${className}
      `}
      aria-label={copied ? 'Link copied to clipboard' : `Copy link to clipboard`}
    >
      {/* Icon */}
      {copied ? (
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
        // Link icon
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
            d="M13.828 10.172a4 4 0 00-5.656 0l-4 4a4 4 0 105.656 5.656l1.102-1.101m-.758-4.899a4 4 0 005.656 0l4-4a4 4 0 00-5.656-5.656l-1.1 1.1"
          />
        </svg>
      )}
      {buttonText}
    </button>
  );
}

export default CopyLinkButton;
