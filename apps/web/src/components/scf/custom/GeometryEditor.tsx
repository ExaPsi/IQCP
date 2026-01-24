/**
 * GeometryEditor - Textarea editor for XYZ molecular geometry input.
 *
 * Provides:
 * - Monospace textarea for XYZ coordinate entry
 * - Real-time validation with inline error/warning display
 * - Unit toggle (Bohr/Angstrom)
 * - Statistics display (atom count, electron count, molecular formula)
 *
 * @module components/scf/custom/GeometryEditor
 */

import { useMemo } from 'react';
import type { CoordinateUnits } from '../../../worker/protocol';
import {
  parseAndValidate,
  type ValidationResult,
} from '../../../lib/geometryValidation';

/**
 * Props for GeometryEditor component.
 */
export interface GeometryEditorProps {
  /** Current XYZ text value */
  value: string;
  /** Callback when text changes */
  onChange: (text: string) => void;
  /** Current coordinate units */
  units: CoordinateUnits;
  /** Callback when units change */
  onUnitsChange: (units: CoordinateUnits) => void;
  /** Disable editing */
  disabled?: boolean;
  /** Placeholder text */
  placeholder?: string;
}

/**
 * Unit options for coordinate display.
 */
const UNIT_OPTIONS: { value: CoordinateUnits; label: string; symbol: string }[] = [
  { value: 'angstrom', label: 'Angstrom', symbol: '\u00C5' },
  { value: 'bohr', label: 'Bohr', symbol: 'a\u2080' },
];

/**
 * Editor for XYZ molecular geometry with live validation.
 *
 * Accepts XYZ format text input:
 * ```
 * H  0.0  0.0  0.0
 * H  0.0  0.0  1.4
 * ```
 *
 * Shows validation errors in red, warnings in amber, and
 * geometry statistics when valid.
 *
 * @example
 * ```tsx
 * const [xyz, setXyz] = useState('H  0.0  0.0  0.0\nH  0.0  0.0  1.4');
 * const [units, setUnits] = useState<CoordinateUnits>('bohr');
 *
 * <GeometryEditor
 *   value={xyz}
 *   onChange={setXyz}
 *   units={units}
 *   onUnitsChange={setUnits}
 * />
 * ```
 */
export function GeometryEditor({
  value,
  onChange,
  units,
  onUnitsChange,
  disabled = false,
  placeholder = 'H  0.0  0.0  0.0\nH  0.0  0.0  1.4',
}: GeometryEditorProps) {
  // Validate geometry on every change
  const validation: ValidationResult = useMemo(() => {
    return parseAndValidate(value);
  }, [value]);

  const { valid, errors, warnings, stats } = validation;

  // Determine border color based on validation state
  const borderColor = !value.trim()
    ? 'border-slate-300'
    : errors.length > 0
      ? 'border-red-400'
      : warnings.length > 0
        ? 'border-amber-400'
        : 'border-green-400';

  return (
    <div>
      {/* Label and Unit Toggle */}
      <div className="flex items-center justify-between mb-2">
        <label
          htmlFor="geometry-editor"
          className="block text-sm font-medium text-slate-700"
        >
          Molecular Geometry (XYZ Format)
        </label>
        <div className="flex items-center gap-2">
          <span className="text-xs text-slate-500">Units:</span>
          <div
            role="radiogroup"
            aria-label="Coordinate units"
            className="flex gap-1"
          >
            {UNIT_OPTIONS.map(({ value: unitValue, label, symbol }) => (
              <button
                key={unitValue}
                type="button"
                role="radio"
                aria-checked={units === unitValue}
                onClick={() => onUnitsChange(unitValue)}
                disabled={disabled}
                className={`px-2 py-1 text-xs rounded transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
                  units === unitValue
                    ? 'bg-blue-600 text-white'
                    : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
                } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
                title={label}
              >
                {symbol}
              </button>
            ))}
          </div>
        </div>
      </div>

      {/* Textarea */}
      <textarea
        id="geometry-editor"
        value={value}
        onChange={(e) => onChange(e.target.value)}
        disabled={disabled}
        placeholder={placeholder}
        rows={6}
        spellCheck={false}
        className={`w-full px-3 py-2 font-mono text-sm border rounded-lg bg-white text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 resize-y ${borderColor} ${
          disabled ? 'opacity-50 cursor-not-allowed bg-slate-50' : ''
        }`}
        aria-label="XYZ geometry input"
        aria-describedby="geometry-validation geometry-stats"
      />

      {/* Validation Messages */}
      <div id="geometry-validation" className="mt-2">
        {/* Errors */}
        {errors.map((error, index) => (
          <div
            key={`error-${index}`}
            className="flex items-start gap-2 text-red-600 text-sm mb-1"
            role="alert"
          >
            <svg className="w-4 h-4 mt-0.5 flex-shrink-0" fill="currentColor" viewBox="0 0 20 20">
              <path
                fillRule="evenodd"
                d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z"
                clipRule="evenodd"
              />
            </svg>
            <span>
              {error.line !== undefined && (
                <span className="font-medium">Line {error.line}: </span>
              )}
              {error.message}
            </span>
          </div>
        ))}

        {/* Warnings */}
        {warnings.map((warning, index) => (
          <div
            key={`warning-${index}`}
            className="flex items-start gap-2 text-amber-600 text-sm mb-1"
            role="alert"
          >
            <svg className="w-4 h-4 mt-0.5 flex-shrink-0" fill="currentColor" viewBox="0 0 20 20">
              <path
                fillRule="evenodd"
                d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z"
                clipRule="evenodd"
              />
            </svg>
            <span>{warning.message}</span>
          </div>
        ))}
      </div>

      {/* Statistics Panel - shown when geometry has atoms */}
      {stats.atomCount > 0 && (
        <div
          id="geometry-stats"
          className={`mt-2 px-3 py-2 rounded-lg text-sm ${
            valid
              ? 'bg-green-50 border border-green-200'
              : 'bg-slate-50 border border-slate-200'
          }`}
        >
          <div className="flex flex-wrap gap-x-4 gap-y-1">
            <div>
              <span className="text-slate-500">Formula: </span>
              <span className="font-medium text-slate-800">{stats.molecularFormula}</span>
            </div>
            <div>
              <span className="text-slate-500">Atoms: </span>
              <span className="font-medium text-slate-800">{stats.atomCount}</span>
            </div>
            <div>
              <span className="text-slate-500">Electrons: </span>
              <span className="font-medium text-slate-800">{stats.electronCount}</span>
            </div>
            <div>
              <span className="text-slate-500">Elements: </span>
              <span className="font-medium text-slate-800">{stats.elements.join(', ')}</span>
            </div>
          </div>
        </div>
      )}

      {/* Help Text */}
      <p className="mt-2 text-xs text-slate-500">
        Enter one atom per line: <code className="bg-slate-100 px-1 rounded">Element X Y Z</code>
        {' '}(e.g., <code className="bg-slate-100 px-1 rounded">H 0.0 0.0 0.0</code>).
        Supported elements: H-Ne. Maximum 15 atoms.
      </p>
    </div>
  );
}

export default GeometryEditor;
