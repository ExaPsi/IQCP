/**
 * ImaginaryWarningBanner — red warning banner for imaginary vibrational modes.
 *
 * Rendered above the Frequency tab's summary panel whenever any frequency is
 * negative (by convention, imaginary modes are reported as negative reals).
 * Imaginary modes typically indicate:
 *   - A transition state (exactly 1 imaginary)
 *   - An unconverged or non-equilibrium geometry (1+ imaginary)
 *
 * Provides two actions:
 *   - "Show mode" jumps the selected mode to the first imaginary frequency
 *     (so the normal-mode viewer displays its displacement arrows).
 *   - "Re-optimize" switches the Module C workflow tab to the Optimize tab.
 *
 * @module components/scf/ImaginaryWarningBanner
 * @see US-102 Frequency Tab UI, AC6
 */

import type { FrequencyResult } from '../../worker/protocol';
import {
  countImaginary,
  firstImaginaryIndex,
} from './frequencyPanelLogic';

// ============================================================================
// Props
// ============================================================================

export interface ImaginaryWarningBannerProps {
  /** Frequency result — read `frequenciesCm1` to count imaginary modes. */
  result: FrequencyResult | null;
  /** Callback when user clicks "Show mode" — selects the first imaginary mode. */
  onShowMode: (modeIndex: number) => void;
  /** Callback when user clicks "Re-optimize" — switches to the Optimize tab. */
  onNavigateToOptimize: () => void;
}

// ============================================================================
// Component
// ============================================================================

/**
 * Red warning banner shown when any frequency is imaginary.
 *
 * Returns `null` if there is no result or if all frequencies are real.
 *
 * @see US-102 Frequency Tab UI, AC6
 */
export function ImaginaryWarningBanner({
  result,
  onShowMode,
  onNavigateToOptimize,
}: ImaginaryWarningBannerProps): JSX.Element | null {
  const count = countImaginary(result);
  if (count === 0) return null;

  const firstIdx = firstImaginaryIndex(result);
  const plural = count !== 1;

  return (
    <div
      role="alert"
      aria-live="polite"
      className="bg-red-100 border border-red-400 rounded-lg px-4 py-3 mb-4 flex flex-col sm:flex-row sm:items-center sm:justify-between gap-3"
    >
      <div className="flex items-start gap-2">
        <svg
          className="w-5 h-5 flex-shrink-0 text-red-600 mt-0.5"
          fill="none"
          stroke="currentColor"
          strokeWidth={2}
          viewBox="0 0 24 24"
          aria-hidden="true"
        >
          <path
            strokeLinecap="round"
            strokeLinejoin="round"
            d="M12 9v2m0 4h.01M10.29 3.86L1.82 18a2 2 0 001.71 3h16.94a2 2 0 001.71-3L13.71 3.86a2 2 0 00-3.42 0z"
          />
        </svg>
        <div className="text-sm text-red-800">
          <span className="font-semibold">
            {count} imaginary frequenc{plural ? 'ies' : 'y'} detected
          </span>{' '}
          — likely transition state or unconverged geometry.
        </div>
      </div>
      <div className="flex gap-2 flex-shrink-0">
        {firstIdx >= 0 && (
          <button
            type="button"
            onClick={() => onShowMode(firstIdx)}
            className="px-3 py-1.5 bg-white border border-red-300 text-red-700 rounded text-xs font-semibold hover:bg-red-50 focus:outline-none focus-visible:ring-2 focus-visible:ring-red-500"
            aria-label="Show first imaginary mode in the 3D viewer"
          >
            Show mode
          </button>
        )}
        <button
          type="button"
          onClick={onNavigateToOptimize}
          className="px-3 py-1.5 bg-red-600 text-white rounded text-xs font-semibold hover:bg-red-700 focus:outline-none focus-visible:ring-2 focus-visible:ring-red-500"
          aria-label="Switch to the Optimize tab to re-run geometry optimization"
        >
          Re-optimize
        </button>
      </div>
    </div>
  );
}

export default ImaginaryWarningBanner;
