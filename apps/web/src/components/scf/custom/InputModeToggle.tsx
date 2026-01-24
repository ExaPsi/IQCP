/**
 * InputModeToggle - Toggle between preset and custom geometry modes.
 *
 * Provides a button group for switching between pre-computed preset
 * systems and custom geometry input mode.
 *
 * Uses ARIA radiogroup pattern for accessibility.
 *
 * @module components/scf/custom/InputModeToggle
 */

/**
 * Input mode for SCF module.
 */
export type InputMode = 'preset' | 'custom';

/**
 * Props for InputModeToggle component.
 */
export interface InputModeToggleProps {
  /** Current input mode */
  mode: InputMode;
  /** Callback when mode changes */
  onChange: (mode: InputMode) => void;
  /** Disable the toggle */
  disabled?: boolean;
}

/**
 * Mode options for the toggle.
 */
const MODE_OPTIONS: { value: InputMode; label: string; description: string }[] = [
  {
    value: 'preset',
    label: 'Preset',
    description: 'Use pre-computed molecular systems with optimized integrals',
  },
  {
    value: 'custom',
    label: 'Custom',
    description: 'Enter custom molecular geometry and compute integrals on-the-fly',
  },
];

/**
 * Toggle between preset and custom input modes.
 *
 * Renders as a button group with ARIA radiogroup semantics for
 * keyboard navigation and screen reader support.
 *
 * Note: The parent component should provide the label and wrapper
 * for consistent spacing with other controls.
 *
 * @example
 * ```tsx
 * const [mode, setMode] = useState<InputMode>('preset');
 * <InputModeToggle mode={mode} onChange={setMode} />
 * ```
 */
export function InputModeToggle({
  mode,
  onChange,
  disabled = false,
}: InputModeToggleProps) {
  const handleKeyDown = (event: React.KeyboardEvent, optionValue: InputMode) => {
    // Handle arrow key navigation for radiogroup
    if (event.key === 'ArrowLeft' || event.key === 'ArrowRight') {
      event.preventDefault();
      const newValue = optionValue === 'preset' ? 'custom' : 'preset';
      onChange(newValue);
    }
  };

  return (
    <div className="space-y-1">
      <div
        role="radiogroup"
        aria-label="Input Mode"
        className="flex gap-2"
      >
        {MODE_OPTIONS.map(({ value, label, description }) => (
          <button
            key={value}
            type="button"
            role="radio"
            aria-checked={mode === value}
            aria-label={`${label}: ${description}`}
            tabIndex={mode === value ? 0 : -1}
            onClick={() => onChange(value)}
            onKeyDown={(e) => handleKeyDown(e, value)}
            disabled={disabled}
            className={`flex-1 px-4 py-2 rounded-lg text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
              mode === value
                ? 'bg-blue-600 text-white'
                : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
            } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
          >
            {label}
          </button>
        ))}
      </div>
      <p className="text-xs text-slate-500">
        {MODE_OPTIONS.find((opt) => opt.value === mode)?.description}
      </p>
    </div>
  );
}

export default InputModeToggle;
