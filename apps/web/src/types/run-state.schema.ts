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
export const ModuleIdSchema = z.enum(['boys', 'rys', 'scf']);

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
 */
export const ScfParamsSchema = z
  .object({
    system_id: z.string().min(1).max(100),
    max_iter: z.number().int().min(1).max(500),
    conv: z.enum(['loose', 'medium', 'tight']),
    diis: z.boolean(),
    // Custom mode fields (optional)
    input_mode: z.enum(['preset', 'custom']).optional(),
    custom_geometry: ScfCustomGeometrySchema.optional(),
    custom_basis: z.string().min(1).max(20).optional(),
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
