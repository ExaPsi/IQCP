/**
 * MethodSelector - Dropdown for selecting computation method.
 *
 * Displays a dropdown with RHF, LDA, B3LYP, and B3LYP-D3(BJ) options.
 * B3LYP-D3(BJ) is shown but disabled with a "coming soon" note since
 * the D3 dispersion correction is not yet implemented.
 *
 * Below the dropdown, shows a brief description of the selected method.
 *
 * @module components/scf/MethodSelector
 * @see US-070 DFT UI + Method Selector + Deep Links
 */

import { DFT_METHODS, getDftMethodInfo } from '../../types/dft';
import type { DftMethod } from '../../types/dft';

/**
 * Props for the MethodSelector component.
 */
export interface MethodSelectorProps {
  /** Currently selected method */
  method: DftMethod;
  /** Callback when method changes */
  onChange: (method: DftMethod) => void;
  /** Whether the selector is disabled */
  disabled?: boolean;
}

/**
 * Method selector dropdown for choosing between RHF and DFT methods.
 *
 * @example
 * ```tsx
 * <MethodSelector
 *   method={method}
 *   onChange={setMethod}
 *   disabled={isRunning}
 * />
 * ```
 */
export function MethodSelector({ method, onChange, disabled = false }: MethodSelectorProps) {
  const info = getDftMethodInfo(method);

  return (
    <div>
      <label
        htmlFor="scf-method-select"
        className="block text-sm font-medium text-slate-700 mb-2"
      >
        Method
      </label>
      <select
        id="scf-method-select"
        value={method}
        onChange={(e) => {
          const selected = e.target.value as DftMethod;
          // Only allow available methods
          const methodInfo = getDftMethodInfo(selected);
          if (methodInfo.available) {
            onChange(selected);
          }
        }}
        disabled={disabled}
        className="w-full px-3 py-2 border border-slate-300 rounded-lg bg-white text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
        aria-label="Select computation method"
      >
        {DFT_METHODS.map((m) => (
          <option
            key={m.id}
            value={m.id}
            disabled={!m.available}
          >
            {m.label}{!m.available ? ' (coming soon)' : ''}
          </option>
        ))}
      </select>
      <p className="mt-1 text-xs text-slate-500">
        {info.description}
      </p>
    </div>
  );
}

export default MethodSelector;
