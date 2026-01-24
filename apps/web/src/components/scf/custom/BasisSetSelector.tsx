/**
 * BasisSetSelector - Dropdown selector for basis sets.
 *
 * Provides a dropdown to select from supported basis sets for
 * custom geometry integral computation.
 *
 * @module components/scf/custom/BasisSetSelector
 */

import type { BasisSetName } from '../../../worker/protocol';

/**
 * Props for BasisSetSelector component.
 */
export interface BasisSetSelectorProps {
  /** Currently selected basis set */
  value: BasisSetName;
  /** Callback when basis set changes */
  onChange: (basis: BasisSetName) => void;
  /** Disable the selector */
  disabled?: boolean;
}

/**
 * Information about each supported basis set.
 */
interface BasisSetInfo {
  value: BasisSetName;
  label: string;
  description: string;
  /** Relative computational cost (1-5 scale) */
  cost: number;
}

/**
 * Available basis sets with descriptions.
 */
const BASIS_SETS: BasisSetInfo[] = [
  {
    value: 'sto-3g',
    label: 'STO-3G',
    description: 'Minimal basis set, fast but low accuracy',
    cost: 1,
  },
  {
    value: '3-21g',
    label: '3-21G',
    description: 'Split-valence basis, improved accuracy',
    cost: 2,
  },
  {
    value: '6-31g',
    label: '6-31G',
    description: 'Double-zeta split-valence, good balance',
    cost: 3,
  },
  {
    value: '6-31g*',
    label: '6-31G*',
    description: 'With d polarization on heavy atoms',
    cost: 4,
  },
  {
    value: '6-31+g*',
    label: '6-31+G*',
    description: 'With diffuse sp + d polarization (good for anions)',
    cost: 5,
  },
];

/**
 * Render cost indicator dots.
 */
function CostIndicator({ cost }: { cost: number }) {
  return (
    <span className="inline-flex gap-0.5 ml-2" aria-label={`Cost: ${cost} of 5`}>
      {[1, 2, 3, 4, 5].map((level) => (
        <span
          key={level}
          className={`w-1.5 h-1.5 rounded-full ${
            level <= cost ? 'bg-blue-500' : 'bg-slate-200'
          }`}
        />
      ))}
    </span>
  );
}

/**
 * Dropdown selector for basis sets.
 *
 * Shows basis set name, description, and relative computational cost.
 *
 * @example
 * ```tsx
 * const [basis, setBasis] = useState<BasisSetName>('sto-3g');
 * <BasisSetSelector value={basis} onChange={setBasis} />
 * ```
 */
export function BasisSetSelector({
  value,
  onChange,
  disabled = false,
}: BasisSetSelectorProps) {
  const currentBasis = BASIS_SETS.find((b) => b.value === value);

  return (
    <div>
      <label
        htmlFor="basis-set-select"
        className="block text-sm font-medium text-slate-700 mb-2"
      >
        Basis Set
      </label>
      <select
        id="basis-set-select"
        value={value}
        onChange={(e) => onChange(e.target.value as BasisSetName)}
        disabled={disabled}
        className="w-full px-3 py-2 border border-slate-300 rounded-lg bg-white text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
        aria-label="Select basis set"
        aria-describedby="basis-set-description"
      >
        {BASIS_SETS.map(({ value: basisValue, label }) => (
          <option key={basisValue} value={basisValue}>
            {label}
          </option>
        ))}
      </select>
      {currentBasis && (
        <div id="basis-set-description" className="mt-2 text-xs text-slate-500">
          <div className="flex items-center">
            <span>{currentBasis.description}</span>
            <CostIndicator cost={currentBasis.cost} />
          </div>
        </div>
      )}
    </div>
  );
}

export default BasisSetSelector;
