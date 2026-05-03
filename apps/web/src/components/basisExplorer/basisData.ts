/**
 * Static basis set shell data for the Basis Explorer UI.
 *
 * This data mirrors the shell structure defined in `crates/qc-core/src/basis/builtin.rs`.
 * It is used for display purposes in Module D. The actual numerical coefficients
 * live in the Rust/WASM layer -- this module only stores structural metadata
 * (angular momentum type, primitive count, shell label).
 *
 * Elements H-Ar (Z=1-18) are fully supported.
 *
 * @module components/basisExplorer/basisData
 */

/**
 * Angular momentum type identifier.
 */
export type AngularMomentumType = 's' | 'p' | 'd';

/**
 * Shell information for display in the UI.
 */
export interface ShellInfo {
  /** Human-readable label, e.g. "1s", "2sp", "3d" */
  label: string;
  /** Angular momentum quantum number: 0=s, 1=p, 2=d */
  angularMomentum: number;
  /** Angular momentum letter */
  angularMomentumLabel: AngularMomentumType;
  /** Number of primitive Gaussian functions in this shell */
  nPrimitives: number;
  /** Brief description of the shell role */
  description: string;
}

/**
 * Element information for the periodic table selector.
 */
export interface ElementInfo {
  /** Atomic number */
  z: number;
  /** Element symbol */
  symbol: string;
  /** Element name */
  name: string;
  /** Row in periodic table (1-3) */
  row: number;
  /** Column in periodic table (1-18) */
  col: number;
  /** Whether basis data is available for this element */
  supported: boolean;
}

/**
 * Basis set option for the selector.
 */
export interface BasisSetOption {
  /** Internal name (lowercase, used as key) */
  value: string;
  /** Display label */
  label: string;
  /** Short description */
  description: string;
}

// =============================================================================
// Element Data
// =============================================================================

export const ELEMENTS: ElementInfo[] = [
  // Row 1
  { z: 1, symbol: 'H', name: 'Hydrogen', row: 1, col: 1, supported: true },
  { z: 2, symbol: 'He', name: 'Helium', row: 1, col: 18, supported: true },
  // Row 2
  { z: 3, symbol: 'Li', name: 'Lithium', row: 2, col: 1, supported: true },
  { z: 4, symbol: 'Be', name: 'Beryllium', row: 2, col: 2, supported: true },
  { z: 5, symbol: 'B', name: 'Boron', row: 2, col: 13, supported: true },
  { z: 6, symbol: 'C', name: 'Carbon', row: 2, col: 14, supported: true },
  { z: 7, symbol: 'N', name: 'Nitrogen', row: 2, col: 15, supported: true },
  { z: 8, symbol: 'O', name: 'Oxygen', row: 2, col: 16, supported: true },
  { z: 9, symbol: 'F', name: 'Fluorine', row: 2, col: 17, supported: true },
  { z: 10, symbol: 'Ne', name: 'Neon', row: 2, col: 18, supported: true },
  // Row 3
  { z: 11, symbol: 'Na', name: 'Sodium', row: 3, col: 1, supported: true },
  { z: 12, symbol: 'Mg', name: 'Magnesium', row: 3, col: 2, supported: true },
  { z: 13, symbol: 'Al', name: 'Aluminium', row: 3, col: 13, supported: true },
  { z: 14, symbol: 'Si', name: 'Silicon', row: 3, col: 14, supported: true },
  { z: 15, symbol: 'P', name: 'Phosphorus', row: 3, col: 15, supported: true },
  { z: 16, symbol: 'S', name: 'Sulfur', row: 3, col: 16, supported: true },
  { z: 17, symbol: 'Cl', name: 'Chlorine', row: 3, col: 17, supported: true },
  { z: 18, symbol: 'Ar', name: 'Argon', row: 3, col: 18, supported: true },
];

// =============================================================================
// Basis Set Options
// =============================================================================

export const BASIS_SETS: BasisSetOption[] = [
  {
    value: 'sto-3g',
    label: 'STO-3G',
    description: 'Minimal basis (3 Gaussians per Slater orbital)',
  },
  {
    value: '3-21g',
    label: '3-21G',
    description: 'Split-valence (inner/outer valence splitting)',
  },
  {
    value: '6-31g',
    label: '6-31G',
    description: 'Split-valence (6 core + 3/1 valence primitives)',
  },
  {
    value: '6-31g*',
    label: '6-31G*',
    description: 'Split-valence + d polarization on heavy atoms',
  },
  {
    value: '6-31+g*',
    label: '6-31+G*',
    description: 'Split-valence + diffuse sp + d polarization',
  },
  {
    value: 'cc-pvdz',
    label: 'cc-pVDZ',
    description: 'Correlation-consistent double-zeta with polarization (Dunning)',
  },
];

// =============================================================================
// Shell Data
// =============================================================================

/**
 * Shell data for each (element Z, basis set) combination.
 *
 * Key format: `${z}-${basisName}`, e.g. "1-sto-3g", "8-6-31g*"
 *
 * Shell labels follow standard spectroscopic notation:
 * - Number = principal quantum number (shell layer)
 * - Letter = angular momentum type (s, p, d)
 *
 * For split-valence bases, the inner/outer designation indicates
 * which part of the split-valence description the shell belongs to.
 */
const SHELL_DATA: Record<string, ShellInfo[]> = {
  // ===========================================================================
  // STO-3G
  // H, He: 1 shell (1s)
  // Li-Ne: 3 shells (1s, 2s, 2p)
  // ===========================================================================
  '1-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core/valence s orbital' },
  ],
  '2-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core/valence s orbital' },
  ],
  '3-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '4-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '5-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '6-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '7-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '8-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '9-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '10-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  // Na-Ar STO-3G: 5 shells (3S + 2P, all 3 primitives)
  '11-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner s orbital' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner p orbital' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '12-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner s orbital' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner p orbital' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '13-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner s orbital' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner p orbital' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '14-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner s orbital' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner p orbital' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '15-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner s orbital' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner p orbital' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '16-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner s orbital' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner p orbital' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '17-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner s orbital' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner p orbital' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],
  '18-sto-3g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner s orbital' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Valence s orbital' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner p orbital' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p orbital' },
  ],

  // ===========================================================================
  // 3-21G
  // H, He: 2 shells (inner s, outer s)
  // Li, Be: 3 shells (core s, inner s, outer s)
  // B-Ne: 5 shells (core s, inner s, outer s, inner p, outer p)
  // ===========================================================================
  '1-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '1s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
  ],
  '2-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '1s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
  ],
  '3-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
  ],
  '4-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
  ],
  '5-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '6-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '7-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '8-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '9-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '10-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '11-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '12-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '13-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '14-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '15-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '16-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '17-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '18-3-21g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 2, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 2, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],

  // ===========================================================================
  // 6-31G
  // H, He: 2 shells (inner s [3 prim], outer s [1 prim])
  // Li-Ne: 5 shells (core s [6 prim], inner valence s [3 prim], outer s [1],
  //                   inner valence p [3 prim], outer p [1 prim])
  // ===========================================================================
  '1-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '1s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
  ],
  '2-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '1s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
  ],
  '3-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '4-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '5-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '6-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '7-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '8-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '9-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '10-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '11-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '12-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '13-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '14-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '15-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '16-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '17-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],
  '18-6-31g': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
  ],

  // ===========================================================================
  // 6-31G*
  // Same as 6-31G, but with an added d polarization shell on Li-Ne (Z >= 3)
  // H, He: 2 shells (same as 6-31G)
  // Li-Ne: 6 shells (6-31G shells + 1 d polarization [1 prim])
  // ===========================================================================
  '1-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '1s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
  ],
  '2-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '1s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
  ],
  '3-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '4-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '5-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '6-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '7-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '8-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '9-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '10-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '11-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '12-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '13-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '14-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '15-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '16-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '17-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '18-6-31g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],

  // ===========================================================================
  // 6-31+G*
  // H, He: 2 shells (same as 6-31G -- no diffuse functions)
  // Li-Ne: 8 shells (core s [6], inner valence s [3], inner valence p [3],
  //                   outer s [1], outer p [1], d polarization [1],
  //                   diffuse s [1], diffuse p [1])
  // ===========================================================================
  '1-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '1s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
  ],
  '2-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '1s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
  ],
  '3-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
    { label: '2s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
  ],
  '4-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
    { label: '2s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
  ],
  '5-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
    { label: '2s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
  ],
  '6-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
    { label: '2s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
  ],
  '7-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
    { label: '2s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
  ],
  '8-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
    { label: '2s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
  ],
  '9-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
    { label: '2s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
  ],
  '10-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '2s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '2p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
    { label: '2s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
  ],
  '11-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '3s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '12-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '3s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '13-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '3s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '14-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '3s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '15-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '3s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '16-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '3s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '17-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '3s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],
  '18-6-31+g*': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Core s orbital' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 6, description: 'Inner core s' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Inner valence s' },
    { label: '3s\'', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Outer valence s' },
    { label: '3s+', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Diffuse s function' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 6, description: 'Inner core p' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Inner valence p' },
    { label: '3p\'', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Outer valence p' },
    { label: '3p+', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Diffuse p function' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d function' },
  ],

  // ===========================================================================
  // cc-pVDZ (Dunning correlation-consistent)
  // H, He: (4s1p) -> [2s1p]: 3 shells, 5 Cartesian basis functions
  // Li-Ne: (9s4p1d) -> [3s2p1d]: 6 shells, 15 Cartesian basis functions
  // Na-Ar: (12s8p1d) -> [4s3p1d]: 8 shells, 19 Cartesian basis functions
  // ===========================================================================
  '1-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s shell' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Polarization p shell' },
  ],
  '2-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 3, description: 'Core s shell' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Polarization p shell' },
  ],
  '3-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell (outer)' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '4-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell (outer)' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '5-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell (outer)' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '6-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell (outer)' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '7-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell (outer)' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '8-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell (outer)' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '9-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell (outer)' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '10-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 8, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 3, description: 'Valence p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell (outer)' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  // Na-Ar cc-pVDZ: (12s8p1d) -> [4s3p1d]: 8 shells, 19 Cartesian basis functions
  '11-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 3)' },
    { label: '4s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (outer)' },
    { label: '4p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '12-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 3)' },
    { label: '4s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (outer)' },
    { label: '4p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '13-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 3)' },
    { label: '4s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (outer)' },
    { label: '4p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '14-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 3)' },
    { label: '4s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (outer)' },
    { label: '4p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '15-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 3)' },
    { label: '4s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (outer)' },
    { label: '4p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '16-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 3)' },
    { label: '4s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (outer)' },
    { label: '4p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '17-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 3)' },
    { label: '4s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (outer)' },
    { label: '4p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
  '18-cc-pvdz': [
    { label: '1s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 1)' },
    { label: '2s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 2)' },
    { label: '3s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 11, description: 'Core s shell (contraction 3)' },
    { label: '4s', angularMomentum: 0, angularMomentumLabel: 's', nPrimitives: 1, description: 'Valence s shell' },
    { label: '2p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (inner)' },
    { label: '3p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 7, description: 'Core p shell (outer)' },
    { label: '4p', angularMomentum: 1, angularMomentumLabel: 'p', nPrimitives: 1, description: 'Valence p shell' },
    { label: '3d', angularMomentum: 2, angularMomentumLabel: 'd', nPrimitives: 1, description: 'Polarization d shell' },
  ],
};

// =============================================================================
// Lookup Functions
// =============================================================================

/**
 * Get the shell list for a given element and basis set.
 *
 * @param z - Atomic number (1-18)
 * @param basisName - Basis set name (e.g., 'sto-3g', '6-31g*')
 * @returns Array of shell info, or null if the combination is not supported
 */
export function getShells(z: number, basisName: string): ShellInfo[] | null {
  const key = `${z}-${basisName}`;
  return SHELL_DATA[key] ?? null;
}

/**
 * Get element info by atomic number.
 *
 * @param z - Atomic number (1-18)
 * @returns Element info or undefined if not in the table
 */
export function getElement(z: number): ElementInfo | undefined {
  return ELEMENTS.find((el) => el.z === z);
}

/**
 * Check if an element is supported (has basis data available).
 *
 * @param z - Atomic number
 * @returns true if basis data exists for this element
 */
export function isElementSupported(z: number): boolean {
  return z >= 1 && z <= 18;
}
