/**
 * Zod Validation Schema for RunStateV1
 *
 * Provides runtime validation for URL-decoded state objects with
 * structured error messages for debugging invalid deep links.
 *
 * @module run-state.schema
 */

import { z } from 'zod';
import type { RunStateV1 } from './run-state';

/**
 * Module ID validation schema
 *
 * Validates the module identifier used for routing.
 */
export const ModuleIdSchema = z.enum(['boys', 'rys', 'scf', 'basis', 'integrals']);

/**
 * Boys parameters validation schema
 *
 * Validates Boys module parameters with sensible bounds:
 * - m: 0-50 (practical range for Boys functions)
 * - T: 0-1000 (covers all computational regimes)
 * - view: single or sweep mode
 */
export const BoysParamsSchema = z.object({
  m: z.number().int().min(0).max(50),
  T: z.number().min(0).max(1000),
  view: z.enum(['single', 'sweep']),
});

/**
 * Rys parameters validation schema
 *
 * Validates Rys quadrature parameters:
 * - n: 1-10 (practical quadrature orders)
 * - T: 0-100 (typical moment parameter range)
 * - target: predefined accuracy levels
 */
export const RysParamsSchema = z.object({
  n: z.number().int().min(1).max(10),
  T: z.number().min(0).max(100),
  target: z.enum(['1e-4', '1e-6', '1e-8']),
});

/**
 * Coordinate tuple for atom positions [x, y, z].
 */
export const CoordinateTupleSchema = z.tuple([z.number(), z.number(), z.number()]);

/**
 * Custom geometry atom validation schema
 *
 * Validates individual atoms in custom geometry:
 * - symbol: element symbol (H-Ne supported)
 * - xyz: 3-element coordinate tuple
 */
export const ScfCustomAtomSchema = z.object({
  symbol: z.string().min(1).max(2),
  xyz: CoordinateTupleSchema,
});

/**
 * Custom geometry validation schema
 *
 * Validates the complete geometry specification:
 * - atoms: 1-50 atoms (practical limit for educational use)
 * - units: bohr or angstrom
 */
export const ScfCustomGeometrySchema = z.object({
  atoms: z.array(ScfCustomAtomSchema).min(1).max(50),
  units: z.enum(['bohr', 'angstrom']),
});

/**
 * PES scan configuration validation schema
 *
 * Validates PES scan parameters for both legacy diatomic scans and
 * Phase 4 internal coordinate scans (bond/angle/dihedral).
 *
 * Legacy format: { r_min, r_max, n_points }
 * Extended format: { r_min, r_max, n_points, coordinate_type, atom_indices, scan_mode }
 *
 * All Phase 4 fields are optional for backward compatibility.
 * r_min/r_max max raised to 360 to accommodate angle values in degrees.
 *
 * @see US-040 PES Scan UI (original diatomic)
 * @see US-085 Deep Links + Backward Compatibility
 */
export const PesScanConfigSchema = z.object({
  r_min: z.number().min(0.3).max(360),
  r_max: z.number().min(0.3).max(360),
  n_points: z.number().int().min(2).max(100),
  // Phase 4 extensions (optional for backward compat, US-085)
  coordinate_type: z.enum(['bond', 'angle', 'dihedral']).optional(),
  atom_indices: z.array(z.number().int().min(0).max(49)).min(2).max(4).optional(),
  scan_mode: z.enum(['rigid', 'relaxed']).optional(),
});

/**
 * SCF parameters validation schema
 *
 * Validates SCF module parameters:
 * - system_id: non-empty identifier for curated systems
 * - max_iter: 1-500 (practical iteration limits)
 * - conv: predefined convergence profiles
 * - diis: boolean for DIIS acceleration toggle
 * - input_mode: optional preset/custom mode selector
 * - custom_geometry: optional geometry for custom mode
 * - custom_basis: optional basis set for custom mode
 * - pes: optional PES scan configuration (Phase 2)
 */
export const ScfParamsSchema = z
  .object({
    system_id: z.string().min(1).max(100),
    max_iter: z.number().int().min(1).max(500),
    conv: z.enum(['loose', 'medium', 'tight']),
    diis: z.boolean(),
    // Method selection (optional, US-070)
    method: z.enum(['rhf', 'lda', 'b3lyp', 'b3lyp-d3bj']).optional(),
    // Custom mode fields (optional)
    input_mode: z.enum(['preset', 'custom']).optional(),
    custom_geometry: ScfCustomGeometrySchema.optional(),
    custom_basis: z.string().min(1).max(20).optional(),
    // PES scan configuration (optional, Phase 2)
    pes: PesScanConfigSchema.optional(),
    // Orbital viewer configuration (optional, Phase 2 US-044)
    orbital: z.number().optional(),
    isovalue: z.number().optional(),
    // Fock matrix damping factor (optional, Phase 2)
    damp: z.number().optional(),
    // Basis set comparison selection (optional, Phase 2 US-045)
    comparison_bases: z.array(z.string()).optional(),
    // Density visualization parameters (optional, Phase 3 US-062, US-063)
    density_mode: z.enum(['total', 'difference']).optional(),
    density_isovalue: z.number().min(0.01).max(0.50).optional(),
    diff_isovalue: z.number().min(0.001).max(0.05).optional(),
    plane_axis: z.enum(['xy', 'xz', 'yz']).optional(),
    plane_position: z.number().optional(),
    color_scale: z.enum(['linear', 'log']).optional(),
    // Frequency analysis parameters (optional, Phase 5 US-103)
    freq_mode: z.number().int().min(0).max(999).optional(),
    freq_temp: z.number().min(10).max(2000).optional(),
    freq_pres: z.number().min(100).max(1e8).optional(),
    freq_broad: z.enum(['l', 'g']).optional(),
    freq_fwhm: z.number().min(0.5).max(100).optional(),
    freq_tab: z.enum(['ir', 'raman']).optional(),
    freq_arrows: z.boolean().optional(),
    freq_units: z.enum(['ha', 'kc', 'kj']).optional(),
  })
  .refine(
    (data) => {
      // If custom mode, require custom_geometry and custom_basis
      if (data.input_mode === 'custom') {
        return data.custom_geometry !== undefined && data.custom_basis !== undefined;
      }
      return true;
    },
    { message: 'Custom mode requires custom_geometry and custom_basis' }
  );

/**
 * Basis Explorer parameters validation schema
 *
 * Validates Module D parameters:
 * - element: atomic number 1-18
 * - basis: basis set name
 * - shell: shell index (nullable)
 * - compare: comparison mode flag (optional)
 * - compare_bases: basis sets for comparison (optional)
 * - compare_shell_type: shell type label for comparison (optional)
 *
 * @see US-052 Basis Set Comparison Mode
 */
export const BasisParamsSchema = z.object({
  element: z.number().int().min(1).max(18),
  basis: z.string().min(1).max(20),
  shell: z.number().int().min(0).nullable(),
  compare: z.boolean().optional(),
  compare_bases: z.array(z.string().min(1).max(20)).max(5).optional(),
  compare_shell_type: z.string().min(1).max(10).optional(),
  exponent_overrides: z.record(
    z.string(),
    z.array(z.number().positive()),
  ).optional(),
  // US-054: Overlap vs. Distance Plot
  overlap_mode: z.boolean().optional(),
  overlap_traces: z.array(z.object({
    element_b: z.number().int().min(1).max(18),
    basis_b: z.string().min(1).max(20),
    shell_index_b: z.number().int().min(0),
  })).max(4).optional(),
});

/**
 * Integral Inspector parameters validation schema
 *
 * Validates Module E parameters:
 * - system_id: preset system identifier (optional)
 * - tab: active matrix tab S/T/V/Hcore (optional)
 * - cell: selected cell [row, col] (optional)
 * - custom: custom geometry specification (optional)
 * - fock_step: Fock build step 1/2/3 (optional, US-060)
 * - fock_element: selected Fock element [row, col] (optional, US-060)
 * - eri_quartet: selected ERI quartet [i, j, k, l] (optional, US-060)
 *
 * @see US-055 Integral Matrix Heatmap
 * @see US-060 Module E Deep Links
 */
export const IntegralParamsSchema = z.object({
  system_id: z.string().min(1).max(100).optional(),
  tab: z.enum(['S', 'T', 'V', 'Hcore']).optional(),
  cell: z.tuple([z.number().int().min(0), z.number().int().min(0)]).optional(),
  custom: z.object({
    atoms: z.array(z.object({
      symbol: z.string().min(1).max(2),
      xyz: z.tuple([z.number(), z.number(), z.number()]),
    })).min(1).max(50),
    units: z.enum(['bohr', 'angstrom']),
    basis: z.string().min(1).max(20),
  }).optional(),
  // Fock Build Tracing (US-058, deep-linked via US-060)
  fock_step: z.union([z.literal(1), z.literal(2), z.literal(3)]).optional(),
  fock_element: z.tuple([z.number().int().min(0), z.number().int().min(0)]).optional(),
  // ERI Browser (US-059, deep-linked via US-060)
  eri_quartet: z.tuple([
    z.number().int().min(0),
    z.number().int().min(0),
    z.number().int().min(0),
    z.number().int().min(0),
  ]).optional(),
});

/**
 * Preset reference validation schema
 *
 * Validates optional preset references for curated configurations.
 * All fields are optional strings.
 */
export const PresetRefSchema = z
  .object({
    system_id: z.string().optional(),
    molecule_id: z.string().optional(),
    basis_id: z.string().optional(),
  })
  .optional();

/**
 * UI settings validation schema
 *
 * Validates UI display settings:
 * - mode: explain (educational) or internals (raw data)
 * - theme: optional light/dark preference
 */
export const UiSettingsSchema = z.object({
  mode: z.enum(['explain', 'internals']),
  theme: z.enum(['light', 'dark']).optional(),
});

/**
 * RunStateV1 validation schema
 *
 * Complete schema for validating URL-decoded state objects.
 * Includes a refinement to ensure the active module has its
 * corresponding parameters defined.
 *
 * Use `.safeParse()` for validation without throwing:
 *
 * @example
 * ```typescript
 * const result = RunStateV1Schema.safeParse(data);
 * if (result.success) {
 *   const state: RunStateV1 = result.data;
 *   // Use validated state
 * } else {
 *   console.error('Validation failed:', result.error.format());
 *   // Handle invalid state
 * }
 * ```
 */
export const RunStateV1Schema = z
  .object({
    schema_version: z.literal('run_state_v1'),
    app_version: z.string().min(1).max(50),
    module: ModuleIdSchema,
    preset: PresetRefSchema,
    boys: BoysParamsSchema.optional(),
    rys: RysParamsSchema.optional(),
    scf: ScfParamsSchema.optional(),
    basis: BasisParamsSchema.optional(),
    integrals: IntegralParamsSchema.optional(),
    ui: UiSettingsSchema,
  })
  .refine(
    (data) => {
      // Ensure the active module has its params defined
      switch (data.module) {
        case 'boys':
          return data.boys !== undefined;
        case 'rys':
          return data.rys !== undefined;
        case 'scf':
          return data.scf !== undefined;
        case 'basis':
          return data.basis !== undefined;
        case 'integrals':
          return data.integrals !== undefined;
        default:
          return false;
      }
    },
    { message: 'Active module must have its parameters defined' }
  );

/**
 * Type guard for RunStateV1
 *
 * Checks if an unknown value is a valid RunStateV1 object.
 * Uses Zod schema validation internally.
 *
 * @param value - The value to check
 * @returns true if value is a valid RunStateV1
 *
 * @example
 * ```typescript
 * const data: unknown = JSON.parse(someString);
 * if (isRunStateV1(data)) {
 *   // data is now typed as RunStateV1
 *   console.log(data.module);
 * }
 * ```
 */
export function isRunStateV1(value: unknown): value is RunStateV1 {
  return RunStateV1Schema.safeParse(value).success;
}

/** Result type from validateRunState when validation succeeds */
export interface ValidationSuccess {
  success: true;
  data: RunStateV1;
}

/** Result type from validateRunState when validation fails */
export interface ValidationFailure {
  success: false;
  error: z.ZodError;
}

/** Union type for validateRunState result */
export type ValidationResult = ValidationSuccess | ValidationFailure;

/**
 * Validate a value against the RunStateV1 schema
 *
 * Returns a result object with either the validated data or error information.
 * This is a thin wrapper around Zod's safeParse for convenience.
 *
 * @param value - The value to validate
 * @returns Object with success boolean and either data or error
 *
 * @example
 * ```typescript
 * const result = validateRunState(data);
 * if (result.success) {
 *   applyState(result.data);
 * } else {
 *   showError(result.error.format());
 * }
 * ```
 */
export function validateRunState(value: unknown): ValidationResult {
  const result = RunStateV1Schema.safeParse(value);
  if (result.success) {
    // Type assertion is safe because the Zod schema validates the structure
    // The assertion is needed because Zod's type inference doesn't perfectly match
    // our TypeScript interfaces (e.g., tuple types from z.tuple vs [number, number, number])
    return { success: true, data: result.data as RunStateV1 };
  }
  return { success: false, error: result.error };
}
