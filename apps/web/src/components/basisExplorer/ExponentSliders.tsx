/**
 * ExponentSliders - Interactive sliders for adjusting primitive Gaussian exponents.
 *
 * Provides one logarithmic-scale slider per primitive in the selected shell.
 * Each slider adjusts the exponent within a factor of 10 from its default value.
 * The overlap metric quantifies how much the modification changed the function.
 *
 * Pedagogical purpose: Students discover that tight primitives control the
 * near-nucleus region, diffuse primitives control the tail, and the contraction
 * coefficients are optimized for specific exponent values.
 *
 * @see US-053 Exponent Sliders
 * @module components/basisExplorer/ExponentSliders
 */

import { memo, useCallback, useId } from 'react';

// ============================================================================
// Types
// ============================================================================

/**
 * Props for ExponentSliders component.
 */
export interface ExponentSlidersProps {
  /** Default (original) exponents from the basis set */
  defaultExponents: number[];
  /** Current (possibly modified) exponents */
  currentExponents: number[];
  /** Angular momentum label ("s", "p", "d") */
  angularMomentumLabel: string;
  /** Shell label (e.g., "1s", "2p", "2s'") */
  shellLabel: string;
  /** Callback when an exponent is changed */
  onExponentChange: (primitiveIndex: number, newAlpha: number) => void;
  /** Callback to reset all exponents to defaults */
  onReset: () => void;
  /** Whether any exponent has been modified */
  isModified: boolean;
  /** Normalized overlap with the original function (0-1) */
  overlapMetric: number;
}

// ============================================================================
// Helpers
// ============================================================================

/**
 * Format an exponent value for display.
 * Uses scientific notation for values >= 10 or <= 0.1, fixed otherwise.
 */
function formatAlpha(alpha: number): string {
  if (alpha >= 10 || alpha < 0.01) {
    return alpha.toExponential(3);
  }
  return alpha.toFixed(4);
}

/**
 * Get the CSS color class for the overlap metric value.
 */
function overlapColor(S: number): string {
  if (S >= 0.90) return 'text-green-700';
  if (S >= 0.50) return 'text-amber-600';
  return 'text-red-600';
}

/**
 * Get the background gradient color for the overlap bar.
 */
function overlapBarColor(S: number): string {
  if (S >= 0.90) return 'bg-green-500';
  if (S >= 0.50) return 'bg-amber-500';
  return 'bg-red-500';
}

// ============================================================================
// Slider Row Component
// ============================================================================

interface PrimitiveSliderProps {
  index: number;
  defaultAlpha: number;
  currentAlpha: number;
  onChange: (index: number, value: number) => void;
  sliderId: string;
}

/**
 * Single primitive exponent slider with logarithmic scale.
 *
 * The slider operates on log10(alpha) internally, where:
 * - min = log10(alpha_default) - 1  (i.e., alpha_default / 10)
 * - max = log10(alpha_default) + 1  (i.e., alpha_default * 10)
 * - midpoint = log10(alpha_default) (i.e., the default value)
 *
 * This gives uniform visual weight per decade of change.
 */
const PrimitiveSlider = memo(function PrimitiveSlider({
  index,
  defaultAlpha,
  currentAlpha,
  onChange,
  sliderId,
}: PrimitiveSliderProps) {
  const logDefault = Math.log10(defaultAlpha);
  const logMin = logDefault - 1; // alpha_default / 10
  const logMax = logDefault + 1; // alpha_default * 10
  const logCurrent = Math.log10(currentAlpha);

  // Ratio to default for display
  const ratio = currentAlpha / defaultAlpha;
  const isModified = Math.abs(ratio - 1.0) > 1e-6;

  const handleChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const logValue = parseFloat(e.target.value);
      const alpha = Math.pow(10, logValue);
      onChange(index, alpha);
    },
    [index, onChange],
  );

  return (
    <div className="flex items-center gap-3">
      {/* Label */}
      <label
        htmlFor={sliderId}
        className="text-xs font-mono text-slate-600 w-6 text-right shrink-0"
      >
        {'\u03B1'}<sub>{index + 1}</sub>
      </label>

      {/* Slider with range labels */}
      <div className="flex-1 flex flex-col">
        <input
          id={sliderId}
          type="range"
          min={logMin}
          max={logMax}
          step={0.005}
          value={logCurrent}
          onChange={handleChange}
          className="w-full h-2 accent-primary-600 cursor-pointer"
          aria-label={`Exponent for primitive ${index + 1}`}
          aria-valuemin={Math.pow(10, logMin)}
          aria-valuemax={Math.pow(10, logMax)}
          aria-valuenow={currentAlpha}
        />
        <div className="flex justify-between mt-0.5">
          <span className="text-[10px] text-slate-400 font-mono">{formatAlpha(Math.pow(10, logMin))}</span>
          <span className="text-[10px] text-slate-400 font-mono">{formatAlpha(Math.pow(10, logMax))}</span>
        </div>
      </div>

      {/* Value display */}
      <span
        className={`text-xs font-mono w-24 text-right shrink-0 ${
          isModified ? 'text-amber-700 font-semibold' : 'text-slate-500'
        }`}
      >
        {formatAlpha(currentAlpha)}
      </span>

      {/* Ratio to default */}
      <span
        className={`text-xs w-14 text-right shrink-0 ${
          isModified ? 'text-amber-600' : 'text-slate-400'
        }`}
      >
        {isModified ? `\u00D7${ratio.toFixed(2)}` : '\u00D71.00'}
      </span>
    </div>
  );
});

// ============================================================================
// Main Component
// ============================================================================

/**
 * Panel of logarithmic sliders for adjusting individual primitive Gaussian exponents.
 *
 * Features:
 * - One slider per primitive with log-scale mapping (equal slider distance = equal multiplicative change)
 * - Current exponent value and ratio-to-default displayed alongside each slider
 * - "Reset to Defaults" button to restore original values
 * - Overlap metric (similarity S) with numerical display and visual bar indicator
 * - Keyboard accessible with aria labels
 *
 * @example
 * ```tsx
 * <ExponentSliders
 *   defaultExponents={[3.425, 0.624, 0.169]}
 *   currentExponents={[3.425, 0.624, 0.169]}
 *   angularMomentumLabel="s"
 *   shellLabel="1s"
 *   onExponentChange={(idx, val) => store.setExponentOverride(shell, idx, val)}
 *   onReset={() => store.resetExponentOverrides(shell)}
 *   isModified={false}
 *   overlapMetric={1.0}
 * />
 * ```
 */
export const ExponentSliders = memo(function ExponentSliders({
  defaultExponents,
  currentExponents,
  angularMomentumLabel,
  shellLabel,
  onExponentChange,
  onReset,
  isModified,
  overlapMetric,
}: ExponentSlidersProps) {
  const baseId = useId();

  // Clamp S to [0, 1] for display (theoretical range is [-1, 1] but practically always positive)
  const displayS = Math.max(0, Math.min(1, overlapMetric));
  const barWidthPercent = displayS * 100;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      {/* Header row */}
      <div className="flex items-center justify-between mb-3">
        <div>
          <h3 className="text-sm font-semibold text-slate-700">
            Modify Gaussian Exponents
          </h3>
          <p className="text-xs text-slate-500 mt-0.5">
            {shellLabel} shell ({angularMomentumLabel}-type,{' '}
            {defaultExponents.length} primitive{defaultExponents.length !== 1 ? 's' : ''})
          </p>
        </div>
        <button
          type="button"
          onClick={onReset}
          disabled={!isModified}
          className={`px-3 py-1.5 text-xs font-medium rounded-md transition-colors
            focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
            ${
              isModified
                ? 'bg-slate-100 text-slate-700 hover:bg-slate-200 border border-slate-300'
                : 'bg-slate-50 text-slate-400 border border-slate-200 cursor-not-allowed'
            }`}
          aria-label="Reset exponents to defaults"
        >
          Reset to Defaults
        </button>
      </div>

      {/* Tooltip */}
      <div className="bg-slate-50 rounded-lg px-3 py-2 mb-3 border border-slate-100">
        <p className="text-xs text-slate-600">
          Each slider adjusts a primitive Gaussian exponent on a logarithmic scale.
          Increasing {'\u03B1'} makes the Gaussian tighter (concentrated near the nucleus);
          decreasing {'\u03B1'} makes it more diffuse (spread farther from the nucleus).
          Range: {'\u03B1'}<sub>default</sub>/10 to {'\u03B1'}<sub>default</sub>{'\u00D7'}10.
        </p>
      </div>

      {/* Slider column headers */}
      <div className="flex items-center gap-3 mb-1 px-0.5">
        <span className="text-xs text-slate-400 w-6 text-right shrink-0"></span>
        <span className="flex-1 text-xs text-slate-400 text-center">tighter &larr; slider &rarr; diffuse</span>
        <span className="text-xs text-slate-400 w-24 text-right shrink-0">value</span>
        <span className="text-xs text-slate-400 w-14 text-right shrink-0">ratio</span>
      </div>

      {/* Primitive sliders */}
      <div className="space-y-2 mb-4">
        {defaultExponents.map((defaultAlpha, idx) => (
          <PrimitiveSlider
            key={idx}
            index={idx}
            defaultAlpha={defaultAlpha}
            currentAlpha={currentExponents[idx]}
            onChange={onExponentChange}
            sliderId={`${baseId}-prim-${idx}`}
          />
        ))}
      </div>

      {/* Overlap metric */}
      <div className="border-t border-slate-100 pt-3">
        <div className="flex items-center justify-between mb-1.5">
          <span className="text-xs text-slate-600 font-medium">
            Similarity to original
          </span>
          <span className={`text-sm font-mono font-semibold ${overlapColor(displayS)}`}>
            S = {displayS.toFixed(4)}
          </span>
        </div>

        {/* Visual bar */}
        <div className="w-full h-2 bg-slate-100 rounded-full overflow-hidden">
          <div
            className={`h-full rounded-full transition-all duration-150 ${overlapBarColor(displayS)}`}
            style={{ width: `${barWidthPercent}%` }}
            role="progressbar"
            aria-valuenow={displayS}
            aria-valuemin={0}
            aria-valuemax={1}
            aria-label={`Overlap metric: ${displayS.toFixed(4)}`}
          />
        </div>

        {/* Educational note when modified */}
        {isModified && displayS < 0.95 && (
          <p className="text-xs text-slate-500 mt-2">
            The contraction coefficients were optimized for the default exponents.
            Changing exponents disrupts this balance, reducing the similarity to the
            original basis function.
          </p>
        )}
      </div>
    </div>
  );
});

export default ExponentSliders;
