/**
 * Deep Link Error Component
 *
 * Displays a user-friendly error message when a deep link cannot be decoded,
 * along with a button to reset to the default state.
 *
 * @module DeepLinkError
 */

import { clearURL } from '../lib/url';

/**
 * Props for the DeepLinkError component
 */
export interface DeepLinkErrorProps {
  /**
   * Callback function invoked when user clicks the reset button.
   * The component will clear the URL automatically; this callback
   * should reset the application state to defaults.
   */
  onReset: () => void;

  /**
   * Optional custom error message.
   * @default 'The link you followed contains invalid data.'
   */
  message?: string;

  /**
   * Optional additional CSS classes for the container.
   */
  className?: string;
}

/**
 * Error UI shown when a deep link cannot be decoded.
 *
 * Provides:
 * - A warning icon and heading
 * - A descriptive error message
 * - Explanation of possible causes
 * - A reset button to clear the invalid URL
 *
 * The component follows accessibility best practices:
 * - Uses role="alert" for screen reader announcements
 * - Uses aria-live="polite" for non-intrusive updates
 * - Button is keyboard-accessible with focus styles
 *
 * @example
 * ```tsx
 * import { DeepLinkError } from '@/components/DeepLinkError';
 * import { DEFAULT_BOYS_STATE } from '@/types/run-state';
 *
 * function BoysPage() {
 *   const [state, setState] = useState<RunStateV1 | null>(null);
 *   const [showError, setShowError] = useState(false);
 *
 *   useEffect(() => {
 *     const urlState = getStateFromURL();
 *     if (hasStateInURL() && !urlState) {
 *       setShowError(true);
 *     } else {
 *       setState(urlState ?? DEFAULT_BOYS_STATE);
 *     }
 *   }, []);
 *
 *   if (showError) {
 *     return (
 *       <DeepLinkError
 *         onReset={() => {
 *           setState(DEFAULT_BOYS_STATE);
 *           setShowError(false);
 *         }}
 *       />
 *     );
 *   }
 *
 *   // ... render normal content
 * }
 * ```
 */
export function DeepLinkError({
  onReset,
  message = 'The link you followed contains invalid data.',
  className = '',
}: DeepLinkErrorProps): React.ReactElement {
  /**
   * Handle reset button click.
   * Clears the URL parameter first, then invokes the callback.
   */
  const handleReset = () => {
    clearURL();
    onReset();
  };

  return (
    <div
      className={`rounded-lg border border-amber-200 bg-amber-50 p-4 ${className}`}
      role="alert"
      aria-live="polite"
    >
      <div className="flex items-start gap-3">
        {/* Warning Icon */}
        <div className="flex-shrink-0">
          <svg
            className="h-5 w-5 text-amber-500"
            viewBox="0 0 20 20"
            fill="currentColor"
            aria-hidden="true"
          >
            <path
              fillRule="evenodd"
              d="M8.485 2.495c.673-1.167 2.357-1.167 3.03 0l6.28 10.875c.673 1.167-.17 2.625-1.516 2.625H3.72c-1.347 0-2.189-1.458-1.515-2.625L8.485 2.495zM10 5a.75.75 0 01.75.75v3.5a.75.75 0 01-1.5 0v-3.5A.75.75 0 0110 5zm0 9a1 1 0 100-2 1 1 0 000 2z"
              clipRule="evenodd"
            />
          </svg>
        </div>

        {/* Content */}
        <div className="flex-1">
          <h3 className="text-sm font-medium text-amber-800">Invalid Deep Link</h3>

          <p className="mt-1 text-sm text-amber-700">{message}</p>

          <p className="mt-1 text-sm text-amber-600">
            This may happen if the link was truncated, modified, or created with an incompatible
            version of the application.
          </p>

          {/* Reset Button */}
          <div className="mt-3">
            <button
              type="button"
              onClick={handleReset}
              className="rounded-md bg-amber-100 px-3 py-2 text-sm font-medium text-amber-800 hover:bg-amber-200 focus:outline-none focus:ring-2 focus:ring-amber-500 focus:ring-offset-2"
            >
              Reset to Default State
            </button>
          </div>
        </div>
      </div>
    </div>
  );
}

export default DeepLinkError;
