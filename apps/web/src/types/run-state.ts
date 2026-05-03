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
export type ModuleId = 'boys' | 'rys' | 'scf' | 'basis' | 'integrals';

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
  /**
   * Selected molecule preset ID (e.g., 'h2', 'h2o').
   * Used with `basis` to determine whether a pre-computed preset exists.
   * Optional for backward compatibility.
   */
  molecule?: string;
  /**
   * Selected basis set ID (e.g., 'sto-3g', '6-31g*').
   * Used with `molecule` to determine whether a pre-computed preset exists.
   * Optional for backward compatibility.
   */
  basis?: string;
  /** Convergence profile determining energy/density thresholds */
  conv: 'loose' | 'medium' | 'tight';
  /** Enable DIIS (Direct Inversion in the Iterative Subspace) acceleration */
  diis: boolean;
  /**
   * Fock matrix damping factor (0.0 = no damping, 0.7 = heavy damping).
   * Optional - defaults to 0.0 if not specified.
   */
  damp?: number;

  /**
   * Computation method: 'rhf', 'lda', 'b3lyp', or 'b3lyp-d3bj'.
   * Optional - defaults to 'rhf' if not specified.
   * @see US-070 DFT UI + Method Selector + Deep Links
   */
  method?: 'rhf' | 'lda' | 'b3lyp' | 'b3lyp-d3bj';

  // === Custom Geometry Mode (optional) ===

  /** Input mode: "preset" (default) or "custom" */
  input_mode?: 'preset' | 'custom';
  /** Custom geometry data (when input_mode === 'custom') */
  custom_geometry?: ScfCustomGeometry;
  /** Basis set for custom geometry (e.g., "sto-3g", "6-31g") */
  custom_basis?: string;

  // === PES Scan Configuration (optional, Phase 2) ===

  /** PES scan configuration (only config, not results -- those are too large for URLs) */
  pes?: PesScanConfig;

  // === Orbital Viewer Configuration (optional, Phase 2 US-044) ===

  /** Selected orbital index for isosurface visualization (0-based) */
  orbital?: number;
  /** Isovalue for isosurface extraction (a.u., default 0.03) */
  isovalue?: number;

  // === Basis Set Comparison (optional, Phase 2 US-045) ===

  /** Selected basis sets for comparison mode */
  comparison_bases?: string[];

  // === Density Visualization (optional, Phase 3 US-062, US-063) ===

  /** Density visualization mode: 'total' or 'difference' */
  density_mode?: 'total' | 'difference';
  /** Density isosurface isovalue in e/bohr^3 */
  density_isovalue?: number;
  /** Difference density isosurface isovalue in e/bohr^3 (US-063) */
  diff_isovalue?: number;
  /** Cross-section plane axis */
  plane_axis?: 'xy' | 'xz' | 'yz';
  /** Cross-section plane position in Bohr */
  plane_position?: number;
  /** Cross-section color scale */
  color_scale?: 'linear' | 'log';

  // === Frequency Analysis (optional, Phase 5 US-103) ===

  /** Selected normal mode index (0-based). Omitted when null. */
  freq_mode?: number;
  /** Temperature in K for thermochemistry display. Omitted when 298.15 (default). */
  freq_temp?: number;
  /** Pressure in Pa for thermochemistry display. Omitted when 101325 (default). */
  freq_pres?: number;
  /** Broadening kind for spectrum: 'l' = Lorentzian, 'g' = Gaussian. Omitted when 'l' (default). */
  freq_broad?: 'l' | 'g';
  /** FWHM in cm^-1 for spectrum broadening. Omitted when 8.0 (default). */
  freq_fwhm?: number;
  /** Spectrum tab: 'ir' or 'raman'. Omitted when 'ir' (default). */
  freq_tab?: 'ir' | 'raman';
  /** Show displacement arrows on normal mode viewer. Omitted when false (default). */
  freq_arrows?: boolean;
  /** Energy units mode for thermo: 'ha' = hartree, 'kc' = kcal/mol, 'kj' = kJ/mol. Omitted when 'kc' (default). */
  freq_units?: 'ha' | 'kc' | 'kj';
}

/**
 * PES scan configuration for deep link serialization.
 *
 * Phase 2 (legacy): Stores only r_min, r_max, n_points for diatomic
 * bond-length scans.
 *
 * Phase 4 (extended): Adds coordinate_type, atom_indices, scan_mode
 * for generalized internal coordinate scans. These fields are optional
 * for backward compatibility -- their absence indicates a legacy
 * diatomic scan.
 *
 * @see US-040 PES Scan UI (original diatomic)
 * @see US-085 Deep Links + Backward Compatibility
 */
export interface PesScanConfig {
  /** Minimum scan value: bohr for bonds, degrees for angles/dihedrals */
  r_min: number;
  /** Maximum scan value */
  r_max: number;
  /** Number of scan points */
  n_points: number;

  // === Phase 4: Internal coordinate fields (optional, backward compat) ===

  /**
   * Coordinate type for internal coordinate scans.
   * Absent in legacy diatomic URLs.
   * @see US-085 Deep Links
   */
  coordinate_type?: 'bond' | 'angle' | 'dihedral';

  /**
   * Atom indices defining the scanned coordinate (0-based).
   * Length matches coordinate type: 2 (bond), 3 (angle), 4 (dihedral).
   * Absent in legacy diatomic URLs.
   */
  atom_indices?: number[];

  /**
   * Scan mode: rigid (frozen geometry) or relaxed (constrained optimization).
   * Absent in legacy diatomic URLs (implies rigid diatomic scan).
   */
  scan_mode?: 'rigid' | 'relaxed';
}

/**
 * Basis Explorer module parameters (Module D)
 *
 * Controls the Basis Explorer view including element selection,
 * basis set selection, shell selection, and comparison mode state.
 *
 * @see US-052 Basis Set Comparison Mode
 */
export interface BasisParams {
  /** Selected element atomic number (1-18) */
  element: number;
  /** Selected basis set name */
  basis: string;
  /** Selected shell index (null = none) */
  shell: number | null;
  /** Comparison mode active */
  compare?: boolean;
  /** Basis sets selected for comparison */
  compare_bases?: string[];
  /** Shell type label for comparison (e.g., "1s", "2p") */
  compare_shell_type?: string;
  /**
   * Exponent overrides: map of shell index -> modified exponents array.
   * Only shells with non-default exponents are included.
   * Omitted if no exponents have been modified.
   *
   * @see US-053 Exponent Sliders
   */
  exponent_overrides?: Record<string, number[]>;
  /**
   * Whether overlap distance mode is active.
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  overlap_mode?: boolean;
  /**
   * Overlap trace configurations for comparison.
   * Each entry defines the B-side of a shell pair.
   *
   * @see US-054 Overlap vs. Distance Plot
   */
  overlap_traces?: Array<{
    element_b: number;
    basis_b: string;
    shell_index_b: number;
  }>;
}

/**
 * Integral Inspector module parameters (Module E)
 *
 * Controls the Integral Inspector view including molecule selection,
 * matrix tab selection, cell selection state, Fock build step, and
 * ERI quartet selection.
 *
 * @see US-055 Integral Matrix Heatmap
 * @see US-058 Fock Build Tracing
 * @see US-059 ERI Browser
 * @see US-060 Module E Deep Links
 */
export interface IntegralParams {
  /** System ID (preset name) or 'custom' for user-entered geometry */
  system_id?: string;
  /** Active matrix tab */
  tab?: 'S' | 'T' | 'V' | 'Hcore';
  /** Selected cell [row, col] */
  cell?: [number, number];
  /** Custom geometry (if not preset) */
  custom?: {
    atoms: Array<{ symbol: string; xyz: [number, number, number] }>;
    units: 'bohr' | 'angstrom';
    basis: string;
  };

  // === Fock Build Tracing (US-058, deep-linked via US-060) ===

  /** Fock build step: 1 = H^core, 2 = G(P) with J/K, 3 = F = H^core + G */
  fock_step?: 1 | 2 | 3;
  /** Selected Fock matrix element [row, col] for decomposition drill-down */
  fock_element?: [number, number];

  // === ERI Browser (US-059, deep-linked via US-060) ===

  /** Selected ERI quartet [i, j, k, l] (0-based basis function indices) */
  eri_quartet?: [number, number, number, number];
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
  /** Display mode: educational explanations, raw internals, or multi-order comparison */
  mode: 'explain' | 'internals' | 'multi-order';
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

  /** Basis Explorer parameters (present when module === 'basis') */
  basis?: BasisParams;

  /** Integral Inspector parameters (present when module === 'integrals') */
  integrals?: IntegralParams;

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
 * Default Basis Explorer state for new users
 *
 * Provides sensible starting values for exploring basis sets:
 * - Hydrogen (simplest element, 1 shell in minimal basis)
 * - STO-3G (minimal basis, good starting point)
 * - No shell selected initially
 */
export const DEFAULT_BASIS_STATE: RunStateV1 = {
  schema_version: 'run_state_v1',
  app_version: APP_VERSION,
  module: 'basis',
  basis: { element: 1, basis: 'sto-3g', shell: null },
  ui: { mode: 'explain' },
};

/**
 * Default Integral Inspector state for new users
 *
 * Provides sensible starting values:
 * - H2 (STO-3G) system (simplest molecule)
 * - Overlap matrix tab (most intuitive starting point)
 * - No cell selected
 */
export const DEFAULT_INTEGRALS_STATE: RunStateV1 = {
  schema_version: 'run_state_v1',
  app_version: APP_VERSION,
  module: 'integrals',
  integrals: { system_id: 'h2_sto3g_r1.4', tab: 'S' },
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
    case 'basis':
      return DEFAULT_BASIS_STATE;
    case 'integrals':
      return DEFAULT_INTEGRALS_STATE;
  }
}
