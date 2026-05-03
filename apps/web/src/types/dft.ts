/**
 * DFT Method Types and Registry
 *
 * Type definitions and metadata for DFT method selection in Module E.
 * Provides a registry of available methods (RHF, LDA, B3LYP, B3LYP-D3BJ)
 * with display metadata and utility functions for routing computations
 * to the correct WASM path.
 *
 * @module types/dft
 * @see US-070 DFT UI + Method Selector + Deep Links
 */

/**
 * DFT method selection for Module E.
 *
 * - 'rhf': Restricted Hartree-Fock (no DFT, uses existing scf_run)
 * - 'lda': Local Density Approximation (Slater exchange + VWN5 correlation)
 * - 'b3lyp': B3LYP hybrid GGA (20% HF exchange + Becke88 + LYP)
 * - 'b3lyp-d3bj': B3LYP with D3(BJ) dispersion correction
 */
export type DftMethod = 'rhf' | 'lda' | 'b3lyp' | 'b3lyp-d3bj';

/**
 * Display metadata for each DFT method.
 */
export interface DftMethodInfo {
  /** Method identifier */
  id: DftMethod;
  /** Display label for the dropdown */
  label: string;
  /** Short description */
  description: string;
  /** Functional name for the info panel (e.g., "Slater + VWN5") */
  functionalName: string | null;
  /** Whether this method uses a numerical integration grid */
  usesGrid: boolean;
  /** Whether this method is currently available */
  available: boolean;
  /** Reason if not available */
  unavailableReason?: string;
}

/**
 * Method info registry.
 *
 * B3LYP-D3(BJ) is listed but disabled until US-071b implements the
 * D3-BJ dispersion correction.
 */
export const DFT_METHODS: DftMethodInfo[] = [
  {
    id: 'rhf',
    label: 'RHF',
    description: 'Restricted Hartree-Fock',
    functionalName: null,
    usesGrid: false,
    available: true,
  },
  {
    id: 'lda',
    label: 'LDA',
    description: 'Local Density Approximation',
    functionalName: 'Slater + VWN5',
    usesGrid: true,
    available: true,
  },
  {
    id: 'b3lyp',
    label: 'B3LYP',
    description: 'Becke 3-parameter Lee-Yang-Parr',
    functionalName: 'B3LYP (0.20 HF + 0.72 B88 + 0.81 LYP + 0.19 VWN5)',
    usesGrid: true,
    available: true,
  },
  {
    id: 'b3lyp-d3bj',
    label: 'B3LYP-D3(BJ)',
    description: 'B3LYP with Grimme D3 dispersion (Becke-Johnson damping)',
    functionalName: 'B3LYP + D3(BJ)',
    usesGrid: true,
    available: true,
  },
];

/**
 * Look up method info by ID.
 */
export function getDftMethodInfo(method: DftMethod): DftMethodInfo {
  return DFT_METHODS.find((m) => m.id === method) ?? DFT_METHODS[0];
}

/**
 * Check if a method requires the ks_scf WASM path.
 */
export function isDftMethod(method: DftMethod): boolean {
  return method !== 'rhf';
}

/**
 * Map DftMethod to the `method` string expected by KsScfRequest.
 *
 * The KsScfRequest.method field accepts the string identifiers used
 * by the Rust ks_scf function (e.g., "lda", "b3lyp").
 * RHF does not use this path at all.
 */
export function toKsMethod(method: DftMethod): string {
  switch (method) {
    case 'rhf': return 'rhf';
    case 'lda': return 'lda';
    case 'b3lyp': return 'b3lyp';
    case 'b3lyp-d3bj': return 'b3lyp-d3bj';
    default: return 'rhf';
  }
}

/**
 * Get the run button label for a given method.
 */
export function getRunButtonLabel(method: DftMethod): string {
  switch (method) {
    case 'rhf': return 'Run SCF';
    case 'lda': return 'Run KS-DFT (LDA)';
    case 'b3lyp': return 'Run KS-DFT (B3LYP)';
    case 'b3lyp-d3bj': return 'Run KS-DFT (B3LYP-D3)';
  }
}

/**
 * Default method for new sessions.
 */
export const DEFAULT_METHOD: DftMethod = 'rhf';
