/**
 * RunStateV1 - URL-encodable application state
 *
 * This module defines the complete state schema needed to reproduce a user's view
 * of any module in the IQCP application. It includes schema versioning for
 * future migration support.
 *
 * @module run-state
 */

/** Module identifier for routing and type discrimination */
export type ModuleId = 'boys' | 'rys' | 'scf';

/**
 * Boys module parameters
 *
 * Controls the visualization of the Boys function F_m(T) with
 * support for single-point evaluation or parameter sweeps.
 */
export interface BoysParams {
  /** Order m of the Boys function F_m(T), must be >= 0 */
  m: number;
  /** Argument T (must be >= 0) */
  T: number;
  /** View mode: single point evaluation or sweep over T range */
  view: 'single' | 'sweep';
}

/**
 * Rys module parameters
 *
 * Controls the Rys quadrature visualization showing roots, weights,
 * and accuracy analysis for different orders.
 */
export interface RysParams {
  /** Number of quadrature points (order), typically 1-10 */
  n: number;
  /** Parameter T for the Boys function moments */
  T: number;
  /** Target accuracy for order recommendation */
  target: '1e-4' | '1e-6' | '1e-8';
}

/**
 * Custom geometry atom for deep link serialization.
 */
export interface ScfCustomAtom {
  /** Element symbol (e.g., "H", "O", "C") */
  symbol: string;
  /** Cartesian coordinates [x, y, z] */
  xyz: [number, number, number];
}

/**
 * Custom geometry data for deep link serialization.
 */
export interface ScfCustomGeometry {
  /** List of atoms in the molecule */
  atoms: ScfCustomAtom[];
  /** Coordinate units: "bohr" or "angstrom" */
  units: 'bohr' | 'angstrom';
}

/**
 * SCF module parameters
 *
 * Controls the Self-Consistent Field computation including
 * convergence settings and DIIS acceleration.
 *
 * Supports both preset mode (using pre-computed integrals) and
 * custom geometry mode (on-the-fly integral computation).
 */
export interface ScfParams {
  /** Identifier for the curated molecular system (or custom_* for custom) */
  system_id: string;
  /** Maximum number of SCF iterations */
  max_iter: number;
  /** Convergence profile determining energy/density thresholds */
  conv: 'loose' | 'medium' | 'tight';
  /** Enable DIIS (Direct Inversion in the Iterative Subspace) acceleration */
  diis: boolean;
  /**
   * Fock matrix damping factor (0.0 = no damping, 0.7 = heavy damping).
   * Optional - defaults to 0.0 if not specified.
   */
  damp?: number;

  // === Custom Geometry Mode (optional) ===

  /** Input mode: "preset" (default) or "custom" */
  input_mode?: 'preset' | 'custom';
  /** Custom geometry data (when input_mode === 'custom') */
  custom_geometry?: ScfCustomGeometry;
  /** Basis set for custom geometry (e.g., "sto-3g", "6-31g") */
  custom_basis?: string;
}

/**
 * Preset reference for curated systems
 *
 * Optional references to pre-configured systems used for
 * demonstration and educational purposes.
 */
export interface PresetRef {
  /** SCF curated system ID (e.g., "h2_sto3g_r1.4") */
  system_id?: string;
  /** Molecule ID for display in Boys/Rys modules */
  molecule_id?: string;
  /** Basis set ID for display in Boys/Rys modules */
  basis_id?: string;
}

/**
 * UI settings
 *
 * Controls the display mode and visual theme of the application.
 */
export interface UiSettings {
  /** Display mode: educational explanations or raw computational internals */
  mode: 'explain' | 'internals';
  /** Color theme (optional, defaults to system preference) */
  theme?: 'light' | 'dark';
}

/**
 * RunStateV1 - URL-encodable application state
 *
 * This is the complete state needed to reproduce a user's view
 * of any module in the application. It includes schema versioning
 * for future migration support.
 *
 * The state is designed to be:
 * - JSON-serializable with deterministic output
 * - Compressible to under 500 characters for URL embedding
 * - Validatable at runtime via Zod schema
 * - Version-aware for future schema migrations
 *
 * @example
 * ```typescript
 * const state: RunStateV1 = {
 *   schema_version: 'run_state_v1',
 *   app_version: '1.0.0-study',
 *   module: 'boys',
 *   boys: { m: 3, T: 5.0, view: 'single' },
 *   ui: { mode: 'explain' }
 * };
 * ```
 */
export interface RunStateV1 {
  /** Schema version for migration support (always 'run_state_v1' for this version) */
  schema_version: 'run_state_v1';
  /** Application version that created this state */
  app_version: string;
  /** Active module determining which view is displayed */
  module: ModuleId;

  /** Preset references (optional) for curated system configurations */
  preset?: PresetRef;

  /** Boys module parameters (present when module === 'boys') */
  boys?: BoysParams;
  /** Rys module parameters (present when module === 'rys') */
  rys?: RysParams;
  /** SCF module parameters (present when module === 'scf') */
  scf?: ScfParams;

  /** UI settings controlling display mode and theme */
  ui: UiSettings;
}

/** Current application version for state creation */
export const APP_VERSION = '0.1.0';

/**
 * Default Boys state for new users
 *
 * Provides sensible starting values for exploring the Boys function:
 * - Order m=0 (simplest case)
 * - T=1.0 (mid-range value in series regime)
 * - Single point view (focused evaluation)
 */
export const DEFAULT_BOYS_STATE: RunStateV1 = {
  schema_version: 'run_state_v1',
  app_version: APP_VERSION,
  module: 'boys',
  boys: { m: 0, T: 1.0, view: 'single' },
  ui: { mode: 'explain' },
};

/**
 * Default Rys state for new users
 *
 * Provides sensible starting values for exploring Rys quadrature:
 * - Order n=3 (common practical value)
 * - T=1.0 (mid-range parameter)
 * - Medium accuracy target (balanced precision)
 */
export const DEFAULT_RYS_STATE: RunStateV1 = {
  schema_version: 'run_state_v1',
  app_version: APP_VERSION,
  module: 'rys',
  rys: { n: 3, T: 1.0, target: '1e-6' },
  ui: { mode: 'explain' },
};

/**
 * Default SCF state for new users
 *
 * Provides sensible starting values for exploring SCF convergence:
 * - H2 molecule with STO-3G basis at R=1.4 bohr (equilibrium)
 * - 50 maximum iterations (sufficient for most cases)
 * - Medium convergence profile
 * - DIIS enabled for faster convergence
 */
export const DEFAULT_SCF_STATE: RunStateV1 = {
  schema_version: 'run_state_v1',
  app_version: APP_VERSION,
  module: 'scf',
  scf: { system_id: 'h2_sto3g_r1.4', max_iter: 50, conv: 'medium', diis: true },
  ui: { mode: 'explain' },
};

/**
 * Get the default state for a given module
 *
 * @param module - The module identifier
 * @returns The default RunStateV1 for that module
 *
 * @example
 * ```typescript
 * const state = getDefaultState('boys');
 * // Returns DEFAULT_BOYS_STATE
 * ```
 */
export function getDefaultState(module: ModuleId): RunStateV1 {
  switch (module) {
    case 'boys':
      return DEFAULT_BOYS_STATE;
    case 'rys':
      return DEFAULT_RYS_STATE;
    case 'scf':
      return DEFAULT_SCF_STATE;
  }
}
