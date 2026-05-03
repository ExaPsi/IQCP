/**
 * OrbitalSelector - Dropdown for selecting which molecular orbital to visualize.
 *
 * Displays a styled <select> element listing all molecular orbitals from
 * the converged SCF calculation. Each option shows the orbital index,
 * energy, and type (occupied/virtual), with HOMO and LUMO specially marked.
 *
 * @module components/scf/OrbitalSelector
 * @see US-044 Orbital Viewer UI
 */

// ============================================================================
// Types
// ============================================================================

/**
 * Props for the OrbitalSelector component.
 */
export interface OrbitalSelectorProps {
  /** Orbital energies in ascending order (Hartree) */
  energies: number[];
  /** Number of occupied orbitals (HOMO index = nOccupied - 1) */
  nOccupied: number;
  /** Currently selected orbital index, or null for none */
  selectedOrbital: number | null;
  /** Callback when a new orbital is selected (null to deselect) */
  onSelect: (idx: number | null) => void;
  /** Disable the selector (e.g., during computation) */
  disabled?: boolean;
}

// ============================================================================
// Helpers
// ============================================================================

/**
 * Format the label for a molecular orbital option.
 *
 * Shows orbital number (1-indexed), energy, and designation:
 * - Core occupied orbitals: no special suffix
 * - HOMO: marked with [HOMO]
 * - LUMO: marked with [LUMO]
 * - Other virtual orbitals: marked as (virtual)
 *
 * @param index - 0-based orbital index
 * @param energy - Orbital energy in Hartree
 * @param nOccupied - Number of occupied orbitals
 * @returns Formatted label string
 */
function formatOrbitalLabel(
  index: number,
  energy: number,
  nOccupied: number
): string {
  const orbitalNumber = index + 1;
  const energyStr = energy.toFixed(4);

  if (index === nOccupied - 1) {
    return `MO ${orbitalNumber}: ${energyStr} Ha [HOMO]`;
  }
  if (index === nOccupied) {
    return `MO ${orbitalNumber}: ${energyStr} Ha [LUMO]`;
  }
  if (index < nOccupied) {
    // Core occupied orbital (not HOMO)
    return `MO ${orbitalNumber}: ${energyStr} Ha`;
  }
  // Virtual orbital (not LUMO)
  return `MO ${orbitalNumber}: ${energyStr} Ha (virtual)`;
}

// ============================================================================
// Component
// ============================================================================

/**
 * OrbitalSelector - dropdown to choose which MO to visualize as an isosurface.
 *
 * Styled to match existing Module E controls (Tailwind, rounded borders,
 * slate color scheme). The "None" option at top clears the orbital selection.
 *
 * @example
 * ```tsx
 * <OrbitalSelector
 *   energies={[-20.55, -1.34, -0.71, -0.57, -0.39, 0.58, 0.69]}
 *   nOccupied={5}
 *   selectedOrbital={4}
 *   onSelect={(idx) => store.setSelectedOrbital(idx)}
 * />
 * ```
 */
export function OrbitalSelector({
  energies,
  nOccupied,
  selectedOrbital,
  onSelect,
  disabled = false,
}: OrbitalSelectorProps) {
  const handleChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    const value = e.target.value;
    if (value === '') {
      onSelect(null);
    } else {
      onSelect(parseInt(value, 10));
    }
  };

  return (
    <div>
      <label
        htmlFor="orbital-selector"
        className="block text-xs font-medium text-slate-600 mb-1"
      >
        Orbital
      </label>
      <select
        id="orbital-selector"
        value={selectedOrbital !== null ? String(selectedOrbital) : ''}
        onChange={handleChange}
        disabled={disabled}
        className="w-full text-sm border border-slate-300 rounded-lg px-3 py-2 bg-white text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-400 focus:border-blue-400 disabled:bg-slate-100 disabled:text-slate-400 disabled:cursor-not-allowed"
        aria-label="Select molecular orbital for 3D visualization"
      >
        <option value="">None</option>
        {energies.map((energy, index) => (
          <option key={index} value={String(index)}>
            {formatOrbitalLabel(index, energy, nOccupied)}
          </option>
        ))}
      </select>
    </div>
  );
}

export default OrbitalSelector;
