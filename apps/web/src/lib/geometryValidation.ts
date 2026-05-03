/**
 * Geometry Validation Utilities
 *
 * Provides parsing and validation for XYZ molecular geometry input.
 * Supports elements H-Ne (matching basis set library limitations).
 *
 * @module lib/geometryValidation
 */

import type { GeometryInput, AtomInput } from '../worker/protocol';
export type { CoordinateUnits } from '../worker/protocol';

// Import for local use and re-export for backward compatibility
import {
  ELEMENT_TO_Z as _ELEMENT_TO_Z,
  isElementSupported as _isElementSupported,
} from './elements';
const isElementSupported = _isElementSupported;
export {
  ELEMENT_TO_Z,
  SUPPORTED_ELEMENTS,
} from './elements';
export { isElementSupported } from './elements';
export type { SupportedElement } from './elements';

/**
 * Element symbol to electron count mapping (neutral atom).
 * For neutral atoms, electron count equals atomic number.
 */
export const ELEMENT_ELECTRONS: Record<string, number> = _ELEMENT_TO_Z;

// ============================================================================
// Validation Types
// ============================================================================

/**
 * Types of validation errors.
 */
export type ValidationErrorType =
  | 'empty_geometry'
  | 'invalid_element'
  | 'parse_error'
  | 'invalid_coordinate'
  | 'too_many_atoms';

/**
 * Types of validation warnings.
 */
export type ValidationWarningType = 'odd_electrons';

/**
 * A validation error with details.
 */
export interface ValidationError {
  type: ValidationErrorType;
  message: string;
  line?: number;
}

/**
 * A validation warning with details.
 */
export interface ValidationWarning {
  type: ValidationWarningType;
  message: string;
}

/**
 * Statistics about the parsed geometry.
 */
export interface GeometryStats {
  /** Number of atoms in the geometry */
  atomCount: number;
  /** Total electron count (sum of atomic numbers) */
  electronCount: number;
  /** Hill notation molecular formula (e.g., "H2O") */
  molecularFormula: string;
  /** Unique elements present */
  elements: string[];
}

/**
 * Result of geometry validation.
 */
export interface ValidationResult {
  /** Whether the geometry is valid for computation */
  valid: boolean;
  /** List of errors (empty if valid) */
  errors: ValidationError[];
  /** List of warnings (may exist even if valid) */
  warnings: ValidationWarning[];
  /** Statistics about the geometry */
  stats: GeometryStats;
  /** Parsed atoms (empty if invalid) */
  atoms: AtomInput[];
}

// ============================================================================
// Constants
// ============================================================================

/** Maximum number of atoms allowed */
export const MAX_ATOMS = 15;

/** Bohr to Angstrom conversion factor */
export const BOHR_TO_ANGSTROM = 0.529177210903;

/** Angstrom to Bohr conversion factor */
export const ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM;

// ============================================================================
// Parsing Functions
// ============================================================================

/**
 * Parse a single line of XYZ format.
 *
 * Expected format: "Symbol X Y Z" with whitespace separation.
 *
 * @param line - The line to parse
 * @param lineNumber - Line number for error reporting (1-indexed)
 * @returns Parsed atom or error
 */
export function parseXyzLine(
  line: string,
  lineNumber: number
): { atom: AtomInput } | { error: ValidationError } {
  const trimmed = line.trim();

  // Skip empty lines and comments
  if (trimmed === '' || trimmed.startsWith('#')) {
    return { atom: null as unknown as AtomInput }; // Skip marker
  }

  // Split by whitespace
  const parts = trimmed.split(/\s+/);

  if (parts.length < 4) {
    return {
      error: {
        type: 'parse_error',
        message: `Could not parse line ${lineNumber}: expected 'Symbol X Y Z'`,
        line: lineNumber,
      },
    };
  }

  const [symbol, xStr, yStr, zStr] = parts;

  // Normalize element symbol (first letter uppercase, rest lowercase)
  const normalizedSymbol = symbol.charAt(0).toUpperCase() + symbol.slice(1).toLowerCase();

  // Check if element is supported
  if (!isElementSupported(normalizedSymbol)) {
    return {
      error: {
        type: 'invalid_element',
        message: `Unsupported element '${symbol}'. Only H-Ar are supported.`,
        line: lineNumber,
      },
    };
  }

  // Parse coordinates
  const x = parseFloat(xStr);
  const y = parseFloat(yStr);
  const z = parseFloat(zStr);

  if (isNaN(x) || isNaN(y) || isNaN(z)) {
    return {
      error: {
        type: 'invalid_coordinate',
        message: `Invalid coordinate on line ${lineNumber}: expected number`,
        line: lineNumber,
      },
    };
  }

  // Check for infinity (can occur with very large numbers)
  if (!isFinite(x) || !isFinite(y) || !isFinite(z)) {
    return {
      error: {
        type: 'invalid_coordinate',
        message: `Coordinate out of range on line ${lineNumber}`,
        line: lineNumber,
      },
    };
  }

  return {
    atom: {
      symbol: normalizedSymbol,
      xyz: [x, y, z],
    },
  };
}

/**
 * Parse XYZ text input into atoms.
 *
 * @param text - Multi-line XYZ format text
 * @returns Array of atoms and any errors encountered
 */
export function parseXyzText(text: string): {
  atoms: AtomInput[];
  errors: ValidationError[];
} {
  const lines = text.split('\n');
  const atoms: AtomInput[] = [];
  const errors: ValidationError[] = [];

  for (let i = 0; i < lines.length; i++) {
    const lineNumber = i + 1;
    const line = lines[i];

    // Skip empty lines
    if (line.trim() === '' || line.trim().startsWith('#')) {
      continue;
    }

    const result = parseXyzLine(line, lineNumber);

    if ('error' in result) {
      errors.push(result.error);
    } else if (result.atom) {
      atoms.push(result.atom);
    }
  }

  return { atoms, errors };
}

// ============================================================================
// Molecular Formula
// ============================================================================

/**
 * Generate Hill notation molecular formula.
 *
 * Hill notation: Carbon first, then Hydrogen, then alphabetical.
 * For molecules without Carbon: alphabetical order.
 *
 * @param atoms - List of atoms
 * @returns Molecular formula string (e.g., "H2O", "CH4")
 */
export function generateMolecularFormula(atoms: AtomInput[]): string {
  if (atoms.length === 0) return '';

  // Count elements
  const counts: Record<string, number> = {};
  for (const atom of atoms) {
    counts[atom.symbol] = (counts[atom.symbol] || 0) + 1;
  }

  // Sort according to Hill notation
  const elements = Object.keys(counts);
  const hasCarbon = elements.includes('C');

  let sortedElements: string[];
  if (hasCarbon) {
    // Carbon first, then Hydrogen, then alphabetical
    sortedElements = [];
    if (counts['C']) {
      sortedElements.push('C');
    }
    if (counts['H']) {
      sortedElements.push('H');
    }
    elements
      .filter((e) => e !== 'C' && e !== 'H')
      .sort()
      .forEach((e) => sortedElements.push(e));
  } else {
    // Alphabetical order
    sortedElements = elements.sort();
  }

  // Build formula string
  let formula = '';
  for (const element of sortedElements) {
    const count = counts[element];
    formula += element;
    if (count > 1) {
      formula += count.toString();
    }
  }

  return formula;
}

// ============================================================================
// Validation
// ============================================================================

/**
 * Validate a parsed geometry.
 *
 * Checks:
 * - Non-empty
 * - Maximum atom count
 * - Even electron count (RHF requirement)
 *
 * @param atoms - Parsed atoms
 * @returns Validation result with errors, warnings, and stats
 */
export function validateGeometry(atoms: AtomInput[]): ValidationResult {
  const errors: ValidationError[] = [];
  const warnings: ValidationWarning[] = [];

  // Check for empty geometry
  if (atoms.length === 0) {
    errors.push({
      type: 'empty_geometry',
      message: 'At least one atom is required',
    });

    return {
      valid: false,
      errors,
      warnings,
      stats: {
        atomCount: 0,
        electronCount: 0,
        molecularFormula: '',
        elements: [],
      },
      atoms: [],
    };
  }

  // Check for too many atoms
  if (atoms.length > MAX_ATOMS) {
    errors.push({
      type: 'too_many_atoms',
      message: `Maximum ${MAX_ATOMS} atoms supported`,
    });
  }

  // Calculate electron count
  let electronCount = 0;
  const elements = new Set<string>();

  for (const atom of atoms) {
    if (isElementSupported(atom.symbol)) {
      electronCount += ELEMENT_ELECTRONS[atom.symbol];
      elements.add(atom.symbol);
    }
  }

  // Check for odd electrons (RHF requirement)
  if (electronCount % 2 !== 0) {
    warnings.push({
      type: 'odd_electrons',
      message: `Warning: Odd electron count (${electronCount}) requires UHF (not supported)`,
    });
  }

  // Generate stats
  const stats: GeometryStats = {
    atomCount: atoms.length,
    electronCount,
    molecularFormula: generateMolecularFormula(atoms),
    elements: Array.from(elements).sort(),
  };

  return {
    valid: errors.length === 0,
    errors,
    warnings,
    stats,
    atoms,
  };
}

/**
 * Full validation pipeline: parse + validate.
 *
 * @param text - XYZ format text input
 * @returns Complete validation result
 */
export function parseAndValidate(text: string): ValidationResult {
  const { atoms, errors: parseErrors } = parseXyzText(text);
  const validationResult = validateGeometry(atoms);

  return {
    ...validationResult,
    errors: [...parseErrors, ...validationResult.errors],
    valid: parseErrors.length === 0 && validationResult.valid,
  };
}

// ============================================================================
// Unit Conversion
// ============================================================================

/**
 * Convert coordinates from Angstrom to Bohr.
 *
 * @param atoms - Atoms with coordinates in Angstrom
 * @returns Atoms with coordinates in Bohr
 */
export function angstromToBohr(atoms: AtomInput[]): AtomInput[] {
  return atoms.map((atom) => ({
    ...atom,
    xyz: [
      atom.xyz[0] * ANGSTROM_TO_BOHR,
      atom.xyz[1] * ANGSTROM_TO_BOHR,
      atom.xyz[2] * ANGSTROM_TO_BOHR,
    ] as [number, number, number],
  }));
}

/**
 * Convert coordinates from Bohr to Angstrom.
 *
 * @param atoms - Atoms with coordinates in Bohr
 * @returns Atoms with coordinates in Angstrom
 */
export function bohrToAngstrom(atoms: AtomInput[]): AtomInput[] {
  return atoms.map((atom) => ({
    ...atom,
    xyz: [
      atom.xyz[0] * BOHR_TO_ANGSTROM,
      atom.xyz[1] * BOHR_TO_ANGSTROM,
      atom.xyz[2] * BOHR_TO_ANGSTROM,
    ] as [number, number, number],
  }));
}

/**
 * Convert geometry to Bohr if needed.
 *
 * @param geometry - Geometry with units
 * @returns Geometry with all coordinates in Bohr
 */
export function ensureBohr(geometry: GeometryInput): GeometryInput {
  if (geometry.units === 'bohr') {
    return geometry;
  }

  return {
    atoms: angstromToBohr(geometry.atoms),
    units: 'bohr',
  };
}

// ============================================================================
// XYZ Text Generation
// ============================================================================

/**
 * Convert atoms to XYZ format text.
 *
 * @param atoms - Atoms to convert
 * @param precision - Decimal precision for coordinates (default: 6)
 * @returns XYZ format text
 */
export function atomsToXyzText(atoms: AtomInput[], precision = 6): string {
  return atoms
    .map(
      (atom) =>
        `${atom.symbol}  ${atom.xyz[0].toFixed(precision)}  ${atom.xyz[1].toFixed(precision)}  ${atom.xyz[2].toFixed(precision)}`
    )
    .join('\n');
}

/**
 * Convert a GeometryInput to XYZ format text.
 *
 * @param geometry - Geometry to convert
 * @param precision - Decimal precision for coordinates
 * @returns XYZ format text
 */
export function geometryToXyzText(geometry: GeometryInput, precision = 6): string {
  return atomsToXyzText(geometry.atoms, precision);
}
