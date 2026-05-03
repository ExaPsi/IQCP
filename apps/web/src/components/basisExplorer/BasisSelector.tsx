/**
 * BasisSelector - Dropdown for selecting a Gaussian basis set.
 *
 * Provides a styled select element with all 5 supported basis sets.
 * Connected to the basisStore for state management.
 *
 * @module components/basisExplorer/BasisSelector
 */

import { useBasisStore } from '../../stores/basisStore';
import { BASIS_SETS } from './basisData';

/**
 * Basis set selector dropdown.
 *
 * Displays all available basis sets with their descriptions.
 * Styled to match existing Module C controls.
 */
export function BasisSelector() {
  const selectedBasis = useBasisStore((state) => state.selectedBasis);
  const setBasis = useBasisStore((state) => state.setBasis);

  const currentBasis = BASIS_SETS.find((b) => b.value === selectedBasis);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
      <h2 className="text-lg font-semibold text-slate-800 mb-4">Basis Set</h2>

      <div className="space-y-3">
        <label htmlFor="basis-select" className="block text-sm font-medium text-slate-700">
          Select basis set
        </label>
        <select
          id="basis-select"
          value={selectedBasis}
          onChange={(e) => setBasis(e.target.value)}
          className="
            w-full px-3 py-2 rounded-lg border border-slate-300
            bg-white text-slate-800 text-sm
            focus:outline-none focus:ring-2 focus:ring-primary-500 focus:border-primary-500
            transition-colors
          "
          aria-label="Basis set"
        >
          {BASIS_SETS.map((basis) => (
            <option key={basis.value} value={basis.value}>
              {basis.label}
            </option>
          ))}
        </select>

        {/* Description of selected basis */}
        {currentBasis && (
          <p className="text-xs text-slate-500 mt-2">
            {currentBasis.description}
          </p>
        )}
      </div>
    </div>
  );
}
