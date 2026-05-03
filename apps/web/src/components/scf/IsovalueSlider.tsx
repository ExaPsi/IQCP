/**
 * IsovalueSlider - Isovalue control for orbital isosurface rendering.
 *
 * Provides a range slider for adjusting the isovalue threshold used
 * in marching cubes extraction. The isovalue determines which points
 * in space are considered "inside" the orbital -- higher values show
 * a smaller, tighter surface; lower values show a larger, more diffuse
 * surface.
 *
 * Includes a contextual tooltip (info icon) explaining the concept
 * for educational purposes.
 *
 * Range: 0.005 to 0.10 a.u., step 0.005
 * Default: 0.03 a.u.
 *
 * @module components/scf/IsovalueSlider
 * @see US-044 Orbital Viewer UI
 */

import { useState, useCallback } from 'react';

// ============================================================================
// Types
// ============================================================================

/**
 * Props for the IsovalueSlider component.
 */
export interface IsovalueSliderProps {
  /** Current isovalue in atomic units */
  isovalue: number;
  /** Callback when isovalue changes */
  onChange: (val: number) => void;
  /** Disable the slider (e.g., when no orbital is selected) */
  disabled?: boolean;
}

// ============================================================================
// Constants
// ============================================================================

/** Minimum isovalue (a.u.) */
const ISO_MIN = 0.005;

/** Maximum isovalue (a.u.) */
const ISO_MAX = 0.1;

/** Step size for isovalue slider */
const ISO_STEP = 0.005;

/** Educational tooltip text */
const TOOLTIP_TEXT =
  'The isovalue is a threshold \u2014 you are seeing all points where ' +
  '|\u03c8| exceeds this value. The orbital extends beyond this surface. ' +
  'Lower values show a larger surface; higher values show a tighter core.';

// ============================================================================
// Component
// ============================================================================

/**
 * IsovalueSlider - range input for orbital isosurface threshold.
 *
 * Styled to match existing Module E slider controls (Tailwind, slate
 * color scheme). Shows the current value formatted to 3 decimal places
 * and includes an info tooltip explaining the concept.
 *
 * @example
 * ```tsx
 * <IsovalueSlider
 *   isovalue={0.03}
 *   onChange={(val) => store.setIsovalue(val)}
 * />
 * ```
 */
export function IsovalueSlider({
  isovalue,
  onChange,
  disabled = false,
}: IsovalueSliderProps) {
  const [showTooltip, setShowTooltip] = useState(false);

  const handleChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      onChange(parseFloat(e.target.value));
    },
    [onChange]
  );

  return (
    <div>
      {/* Label row with info icon */}
      <div className="flex items-center justify-between mb-1">
        <label
          htmlFor="isovalue-slider"
          className="text-xs font-medium text-slate-600"
        >
          Isovalue (a.u.)
        </label>
        <div className="flex items-center gap-2">
          {/* Current value display */}
          <span className="text-xs font-mono text-slate-700">
            {isovalue.toFixed(3)}
          </span>

          {/* Info icon with tooltip */}
          <div className="relative">
            <button
              type="button"
              className="text-slate-400 hover:text-slate-600 focus:outline-none focus:text-blue-500 transition-colors"
              onMouseEnter={() => setShowTooltip(true)}
              onMouseLeave={() => setShowTooltip(false)}
              onFocus={() => setShowTooltip(true)}
              onBlur={() => setShowTooltip(false)}
              aria-label="What is the isovalue?"
              tabIndex={0}
            >
              <svg
                className="w-3.5 h-3.5"
                fill="none"
                stroke="currentColor"
                viewBox="0 0 24 24"
                aria-hidden="true"
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"
                />
              </svg>
            </button>

            {/* Tooltip popover */}
            {showTooltip && (
              <div
                className="absolute right-0 bottom-full mb-2 w-64 p-3 bg-slate-800 text-white text-xs rounded-lg shadow-lg z-50 leading-relaxed"
                role="tooltip"
              >
                {TOOLTIP_TEXT}
                <div className="absolute bottom-0 right-3 translate-y-1/2 rotate-45 w-2 h-2 bg-slate-800" />
              </div>
            )}
          </div>
        </div>
      </div>

      {/* Range slider */}
      <input
        id="isovalue-slider"
        type="range"
        min={ISO_MIN}
        max={ISO_MAX}
        step={ISO_STEP}
        value={isovalue}
        onChange={handleChange}
        disabled={disabled}
        className="w-full h-1.5 rounded-lg appearance-none cursor-pointer bg-slate-200 accent-blue-600 disabled:opacity-50 disabled:cursor-not-allowed"
        aria-label={`Isovalue: ${isovalue.toFixed(3)} atomic units`}
        aria-valuemin={ISO_MIN}
        aria-valuemax={ISO_MAX}
        aria-valuenow={isovalue}
      />

      {/* Min/Max labels */}
      <div className="flex justify-between text-[10px] text-slate-400 mt-0.5">
        <span>{ISO_MIN.toFixed(3)}</span>
        <span>{ISO_MAX.toFixed(3)}</span>
      </div>
    </div>
  );
}

export default IsovalueSlider;
