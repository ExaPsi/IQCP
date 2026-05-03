/**
 * Shared element-to-atomic-number mapping.
 *
 * Single source of truth for element symbol -> Z conversion
 * across the entire IQCP frontend. Supports elements H through Ar (Z=1-18),
 * matching the basis set library limitations.
 *
 * @module lib/elements
 */

// ============================================================================
// Element Data
// ============================================================================

/** Map element symbols to atomic numbers (Z) for H-Ar. */
export const ELEMENT_TO_Z: Record<string, number> = {
  H: 1, He: 2, Li: 3, Be: 4, B: 5,
  C: 6, N: 7, O: 8, F: 9, Ne: 10,
  Na: 11, Mg: 12, Al: 13, Si: 14, P: 15,
  S: 16, Cl: 17, Ar: 18,
};

/** Reverse mapping: atomic number (Z) to element symbol. */
export const Z_TO_ELEMENT: Record<number, string> = {
  1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B',
  6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
  11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
  16: 'S', 17: 'Cl', 18: 'Ar',
};

/** Supported element symbols as a readonly array. */
export const SUPPORTED_ELEMENTS = [
  'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
  'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
] as const;

/** Type for supported element symbols. */
export type SupportedElement = (typeof SUPPORTED_ELEMENTS)[number];

/**
 * Check if an element symbol is supported (Z=1-18).
 */
export function isElementSupported(symbol: string): symbol is SupportedElement {
  return symbol in ELEMENT_TO_Z;
}
