/**
 * PES scan input validation.
 *
 * Provides validation for PES bond-length scan configuration parameters.
 *
 * @module components/scf/pesValidation
 * @see US-040 PES Scan UI
 */

/**
 * Validate PES scan configuration parameters.
 *
 * Checks that the scan range and point count form a valid configuration:
 * - rMin must be at least 0.3 bohr (avoids nuclear overlap)
 * - rMin must be strictly less than rMax
 * - nPoints must be in [2, 100]
 * - All inputs must be valid numbers
 *
 * @param rMin - Minimum bond distance in bohr
 * @param rMax - Maximum bond distance in bohr
 * @param nPoints - Number of scan points
 * @returns Error message string, or null if valid
 *
 * @example
 * ```typescript
 * const error = validatePesConfig(0.5, 5.0, 20);
 * // null (valid)
 *
 * const error2 = validatePesConfig(3.0, 2.0, 20);
 * // "Minimum distance must be less than maximum"
 * ```
 */
export function validatePesConfig(
  rMin: number,
  rMax: number,
  nPoints: number
): string | null {
  if (Number.isNaN(rMin) || Number.isNaN(rMax) || Number.isNaN(nPoints)) {
    return 'All inputs must be valid numbers';
  }
  if (rMin < 0.3) return 'Minimum distance must be at least 0.3 bohr';
  if (rMin >= rMax) return 'Minimum distance must be less than maximum';
  if (nPoints < 2) return 'Need at least 2 scan points';
  if (nPoints > 100) return 'Maximum 100 scan points';
  return null; // valid
}
