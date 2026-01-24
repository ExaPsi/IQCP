/**
 * OrbitalEnergiesTable - Display table for molecular orbital energies.
 *
 * Shows orbital energies with occupation numbers and HOMO/LUMO highlighting.
 * Designed for educational use, clearly showing which orbitals are occupied
 * and which are virtual (unoccupied).
 *
 * Features:
 * - 1-indexed orbital numbers for student clarity
 * - HOMO row highlighted in green
 * - LUMO row highlighted in orange
 * - Occupation column (2 for doubly occupied, 0 for virtual)
 * - Type column (Occupied/Virtual)
 *
 * @module components/scf/OrbitalEnergiesTable
 */

import { Math } from '../common/Math';

/**
 * Props for OrbitalEnergiesTable component.
 */
export interface OrbitalEnergiesTableProps {
  /** Orbital energies in ascending order (Hartree) */
  energies: number[];
  /** Number of occupied orbitals (HOMO index = nOccupied - 1) */
  nOccupied: number;
  /** Additional CSS classes */
  className?: string;
}

/**
 * Format orbital energy with 6 decimal places.
 *
 * @param energy - Orbital energy in Hartree
 * @returns Formatted string
 */
function formatEnergy(energy: number): string {
  return energy.toFixed(6);
}

/**
 * Determine orbital type based on index and occupation.
 *
 * @param index - 0-based orbital index
 * @param nOccupied - Number of occupied orbitals
 * @returns Type label
 */
function getOrbitalType(index: number, nOccupied: number): string {
  if (index < nOccupied) {
    return 'Occupied';
  }
  return 'Virtual';
}

/**
 * Get special designation for HOMO/LUMO orbitals.
 *
 * @param index - 0-based orbital index
 * @param nOccupied - Number of occupied orbitals
 * @returns Special designation or null
 */
function getSpecialDesignation(
  index: number,
  nOccupied: number
): string | null {
  if (index === nOccupied - 1) {
    return 'HOMO';
  }
  if (index === nOccupied) {
    return 'LUMO';
  }
  return null;
}

/**
 * Get row style classes based on orbital type.
 *
 * @param index - 0-based orbital index
 * @param nOccupied - Number of occupied orbitals
 * @returns Tailwind CSS classes
 */
function getRowClasses(index: number, nOccupied: number): string {
  if (index === nOccupied - 1) {
    // HOMO: green background
    return 'bg-green-50 text-green-900';
  }
  if (index === nOccupied) {
    // LUMO: orange background
    return 'bg-orange-50 text-orange-900';
  }
  return 'hover:bg-slate-50';
}

/**
 * Individual orbital row component.
 */
function OrbitalRow({
  index,
  energy,
  nOccupied,
}: {
  index: number;
  energy: number;
  nOccupied: number;
}) {
  const specialDesignation = getSpecialDesignation(index, nOccupied);
  const isOccupied = index < nOccupied;

  return (
    <tr className={`border-b border-slate-100 ${getRowClasses(index, nOccupied)}`}>
      {/* Orbital number (1-indexed) */}
      <td className="px-3 py-2 text-center font-mono">
        {index + 1}
        {specialDesignation && (
          <span
            className={`ml-2 text-xs font-medium ${
              specialDesignation === 'HOMO' ? 'text-green-600' : 'text-orange-600'
            }`}
          >
            ({specialDesignation})
          </span>
        )}
      </td>

      {/* Energy in Hartree */}
      <td className="px-3 py-2 font-mono text-right">
        {formatEnergy(energy)}
      </td>

      {/* Occupation (2 for doubly occupied in RHF, 0 for virtual) */}
      <td className="px-3 py-2 text-center font-mono">
        {isOccupied ? '2' : '0'}
      </td>

      {/* Type */}
      <td className="px-3 py-2 text-center text-xs">
        <span
          className={`font-medium ${
            isOccupied ? 'text-green-700' : 'text-slate-500'
          }`}
        >
          {getOrbitalType(index, nOccupied)}
        </span>
      </td>
    </tr>
  );
}

/**
 * OrbitalEnergiesTable - Molecular orbital energies display.
 *
 * Shows orbital energies with clear HOMO/LUMO identification for
 * educational purposes in the SCF sandbox.
 *
 * @example
 * ```tsx
 * <OrbitalEnergiesTable
 *   energies={[-0.5, -0.2, 0.3, 0.8]}
 *   nOccupied={2}
 * />
 * ```
 */
export function OrbitalEnergiesTable({
  energies,
  nOccupied,
  className = '',
}: OrbitalEnergiesTableProps) {
  // Empty state
  if (energies.length === 0) {
    return (
      <div className={`bg-white rounded-xl shadow-sm border border-slate-200 p-6 ${className}`}>
        <h3 className="text-sm font-semibold text-slate-700 mb-4">
          Orbital Energies
        </h3>
        <div className="text-center py-4 text-slate-500">
          <p className="text-sm">No orbital energy data available</p>
        </div>
      </div>
    );
  }

  // Calculate HOMO/LUMO gap if both exist
  const homoEnergy = nOccupied > 0 ? energies[nOccupied - 1] : null;
  const lumoEnergy = nOccupied < energies.length ? energies[nOccupied] : null;
  const homoLumoGap =
    homoEnergy !== null && lumoEnergy !== null ? lumoEnergy - homoEnergy : null;

  return (
    <div className={`bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden ${className}`}>
      {/* Header */}
      <div className="px-4 py-3 border-b border-slate-200">
        <h3 className="text-sm font-semibold text-slate-700 flex flex-wrap items-center gap-2">
          <span>Orbital Energies</span>
          <Math>{String.raw`\varepsilon_i`}</Math>
          <span className="text-slate-500 font-normal">
            ({energies.length} orbital{energies.length !== 1 ? 's' : ''},{' '}
            {nOccupied} occupied)
          </span>
        </h3>
      </div>

      {/* Table with scrolling */}
      <div className="overflow-auto max-h-64">
        <table className="w-full text-sm">
          <thead className="bg-slate-50 sticky top-0 z-10">
            <tr className="text-slate-600 text-xs uppercase tracking-wider">
              <th className="px-3 py-2 text-center font-medium w-28">Orbital #</th>
              <th className="px-3 py-2 text-right font-medium">Energy (Ha)</th>
              <th className="px-3 py-2 text-center font-medium w-20">Occ.</th>
              <th className="px-3 py-2 text-center font-medium w-24">Type</th>
            </tr>
          </thead>
          <tbody className="text-slate-800">
            {energies.map((energy, index) => (
              <OrbitalRow
                key={index}
                index={index}
                energy={energy}
                nOccupied={nOccupied}
              />
            ))}
          </tbody>
        </table>
      </div>

      {/* HOMO-LUMO Gap Footer */}
      {homoLumoGap !== null && (
        <div className="px-4 py-3 bg-slate-50 border-t border-slate-200">
          <div className="flex flex-wrap justify-between items-center gap-2 text-xs">
            <span className="text-slate-600 flex items-center gap-1">
              <Math>{String.raw`\varepsilon_{\text{LUMO}} - \varepsilon_{\text{HOMO}}`}</Math>
              <span>=</span>
            </span>
            <span className="font-mono font-medium text-slate-800">
              {homoLumoGap.toFixed(6)} Ha ({(homoLumoGap * 27.2114).toFixed(3)} eV)
            </span>
          </div>
        </div>
      )}

      {/* Legend */}
      <div className="px-4 py-2 bg-slate-50 border-t border-slate-100 text-xs text-slate-500 flex gap-4">
        <span className="flex items-center gap-1">
          <span className="w-3 h-3 rounded bg-green-100 border border-green-300" />
          HOMO
        </span>
        <span className="flex items-center gap-1">
          <span className="w-3 h-3 rounded bg-orange-100 border border-orange-300" />
          LUMO
        </span>
      </div>
    </div>
  );
}

export default OrbitalEnergiesTable;
