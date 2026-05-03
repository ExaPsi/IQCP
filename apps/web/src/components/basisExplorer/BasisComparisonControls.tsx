/**
 * BasisComparisonControls - Controls for basis set comparison mode.
 *
 * Provides:
 * - "Compare Basis Sets" toggle button to enter/exit comparison mode
 * - Checkbox list for selecting 2-5 basis sets to compare
 * - Shell type selector dropdown (e.g., "1s", "2s", "2p")
 *
 * The primary basis (from BasisSelector) is always included and cannot
 * be unchecked.
 *
 * @module components/basisExplorer/BasisComparisonControls
 * @see US-052 Basis Set Comparison Mode
 */

import { useMemo, memo } from 'react';
import { BASIS_SETS } from './basisData';
import { COMPARISON_COLORS } from './comparisonColors';
import type { ShellTypeOption } from './shellMatching';

// ============================================================================
// Types
// ============================================================================

export interface BasisComparisonControlsProps {
  /** Whether comparison mode is active */
  isActive: boolean;
  /** Toggle comparison mode on/off */
  onToggle: () => void;
  /** Currently selected basis sets for comparison */
  selectedBases: string[];
  /** Callback to add a basis to comparison */
  onAddBasis: (basis: string) => void;
  /** Callback to remove a basis from comparison */
  onRemoveBasis: (basis: string) => void;
  /** The primary basis set (always included, cannot be unchecked) */
  primaryBasis: string;
  /** Available shell types for the current element */
  shellTypes: ShellTypeOption[];
  /** Currently selected shell type for comparison */
  selectedShellType: string;
  /** Callback when shell type changes */
  onShellTypeChange: (label: string) => void;
  /** Whether a shell is selected (comparison requires shell selection) */
  hasShellSelected: boolean;
}

// ============================================================================
// Component
// ============================================================================

/**
 * Controls panel for basis set comparison mode.
 *
 * Renders a toggle button and, when active, a checkbox list of basis sets
 * and a shell type selector. The primary basis is always checked and
 * disabled. Maximum 5 basis sets can be compared simultaneously.
 */
export const BasisComparisonControls = memo(function BasisComparisonControls({
  isActive,
  onToggle,
  selectedBases,
  onAddBasis,
  onRemoveBasis,
  primaryBasis,
  shellTypes,
  selectedShellType,
  onShellTypeChange,
  hasShellSelected,
}: BasisComparisonControlsProps) {
  // Ensure primary basis is always in the selection
  const effectiveBases = useMemo(() => {
    if (!selectedBases.includes(primaryBasis)) {
      return [primaryBasis, ...selectedBases];
    }
    return selectedBases;
  }, [selectedBases, primaryBasis]);

  const atMaxBases = effectiveBases.length >= 5;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      {/* Toggle button */}
      <div className="flex items-center justify-between mb-3">
        <h3 className="text-sm font-semibold text-slate-700">
          Basis Set Comparison
        </h3>
        <button
          type="button"
          onClick={onToggle}
          disabled={!hasShellSelected}
          aria-pressed={isActive}
          title={
            !hasShellSelected
              ? 'Select a shell to enable comparison'
              : isActive
                ? 'Exit comparison mode'
                : 'Compare this shell across basis sets'
          }
          className={`px-3 py-1.5 text-xs font-medium rounded-lg transition-colors
            focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
            ${
              !hasShellSelected
                ? 'bg-slate-100 text-slate-400 cursor-not-allowed'
                : isActive
                  ? 'bg-primary-600 text-white hover:bg-primary-700'
                  : 'bg-primary-50 text-primary-700 border border-primary-200 hover:bg-primary-100'
            }`}
        >
          {isActive ? 'Exit Comparison' : 'Compare Bases'}
        </button>
      </div>

      {/* Hint when no shell selected */}
      {!hasShellSelected && (
        <p className="text-xs text-slate-400 italic">
          Select a shell from the list above to enable basis set comparison.
        </p>
      )}

      {/* Comparison controls (visible when active) */}
      {isActive && hasShellSelected && (
        <div className="space-y-4">
          {/* Shell type selector */}
          {shellTypes.length > 0 && (
            <div>
              <label
                htmlFor="shell-type-select"
                className="block text-xs font-medium text-slate-600 mb-1"
              >
                Shell type to compare
              </label>
              <select
                id="shell-type-select"
                value={selectedShellType}
                onChange={(e) => onShellTypeChange(e.target.value)}
                className="w-full rounded-md border border-slate-300 px-3 py-1.5 text-sm
                  focus:outline-none focus:ring-2 focus:ring-primary-500 focus:border-primary-500"
              >
                {shellTypes.map((st) => (
                  <option key={st.label} value={st.label}>
                    {st.label} ({['s', 'p', 'd'][st.angularMomentum]}-type)
                  </option>
                ))}
              </select>
            </div>
          )}

          {/* Basis set checkboxes */}
          <div>
            <div className="flex items-center justify-between mb-2">
              <span className="text-xs font-medium text-slate-600">
                Select basis sets
              </span>
              <span className="text-xs text-slate-400">
                {effectiveBases.length} of 5 selected
              </span>
            </div>
            <div className="space-y-1.5">
              {BASIS_SETS.map((bs) => {
                const isPrimary = bs.value === primaryBasis;
                const isChecked = effectiveBases.includes(bs.value);
                const isDisabled = isPrimary || (!isChecked && atMaxBases);
                const color = COMPARISON_COLORS[bs.value] ?? '#6b7280';

                return (
                  <label
                    key={bs.value}
                    className={`flex items-center gap-2 p-1.5 rounded-md transition-colors cursor-pointer
                      ${isDisabled && !isPrimary ? 'opacity-50 cursor-not-allowed' : ''}
                      ${isChecked ? 'bg-slate-50' : 'hover:bg-slate-50'}`}
                  >
                    <input
                      type="checkbox"
                      checked={isChecked}
                      disabled={isDisabled}
                      onChange={() => {
                        if (isPrimary) return;
                        if (isChecked) {
                          onRemoveBasis(bs.value);
                        } else {
                          onAddBasis(bs.value);
                        }
                      }}
                      className="rounded border-slate-300 text-primary-600
                        focus:ring-primary-500 focus:ring-offset-0"
                    />
                    {/* Color indicator */}
                    <span
                      className="w-3 h-3 rounded-sm flex-shrink-0"
                      style={{ backgroundColor: color }}
                      aria-hidden="true"
                    />
                    <span className="text-sm text-slate-700 flex-1">
                      {bs.label}
                    </span>
                    {isPrimary && (
                      <span className="text-xs text-slate-400 italic">
                        primary
                      </span>
                    )}
                  </label>
                );
              })}
            </div>
          </div>
        </div>
      )}
    </div>
  );
});
