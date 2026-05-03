/**
 * Energy unit conversion helpers.
 *
 * Used by the thermochemistry display in the Frequency tab (US-102) to
 * convert between Hartree, kcal/mol, and kJ/mol. Constants match the
 * CODATA 2018 values used in `crates/qc-core/src/constants.rs`.
 *
 * @module lib/units
 * @see US-102 Frequency Tab UI
 */

/**
 * Hartree → kcal/mol conversion factor.
 *
 * 1 Ha = 627.5094740631 kcal/mol (CODATA 2018).
 */
export const HA_TO_KCAL_MOL = 627.5094740631;

/**
 * Hartree → kJ/mol conversion factor.
 *
 * 1 Ha = 2625.4996394799 kJ/mol (CODATA 2018).
 */
export const HA_TO_KJ_MOL = 2625.4996394799;

/**
 * Bohr → Angstrom conversion factor.
 *
 * 1 bohr = 0.529177210903 Angstrom (CODATA 2018).
 */
export const BOHR_TO_ANGSTROM_PRECISE = 0.529177210903;

/**
 * Energy display mode used by the Frequency tab UI.
 */
export type EnergyUnitsMode = 'hartree' | 'kcal_mol' | 'kj_mol';

/**
 * Convert an energy in Hartree to the requested display units.
 *
 * @param valueHa - Energy value in Hartree
 * @param mode    - Target units
 * @returns Converted value in the target units
 */
export function convertEnergy(valueHa: number, mode: EnergyUnitsMode): number {
  switch (mode) {
    case 'hartree':
      return valueHa;
    case 'kcal_mol':
      return valueHa * HA_TO_KCAL_MOL;
    case 'kj_mol':
      return valueHa * HA_TO_KJ_MOL;
  }
}

/**
 * Convert an entropy (or heat capacity) in Ha/(mol·K) to the requested display units.
 *
 * Uses the same per-mole factor as bulk energy (Ha → kcal/mol → kJ/mol), since
 * Ha/(mol·K) scales by the same conversion constants as Hartree.
 *
 * @param valueHaPerK - Entropy / heat capacity in Ha/(mol·K)
 * @param mode        - Target units
 * @returns Converted value
 */
export function convertEntropyPerK(
  valueHaPerK: number,
  mode: EnergyUnitsMode
): number {
  // Entropies and heat capacities share the same scaling as energies when
  // written in Ha/(mol·K) → kcal/(mol·K) or kJ/(mol·K).
  return convertEnergy(valueHaPerK, mode);
}

/**
 * Human-readable label for an energy unit mode.
 */
export function energyUnitsLabel(mode: EnergyUnitsMode): string {
  switch (mode) {
    case 'hartree':
      return 'Ha';
    case 'kcal_mol':
      return 'kcal/mol';
    case 'kj_mol':
      return 'kJ/mol';
  }
}

/**
 * Human-readable label for an entropy / heat-capacity unit.
 */
export function entropyUnitsLabel(mode: EnergyUnitsMode): string {
  switch (mode) {
    case 'hartree':
      return 'Ha/(mol·K)';
    case 'kcal_mol':
      return 'kcal/(mol·K)';
    case 'kj_mol':
      return 'kJ/(mol·K)';
  }
}
