/**
 * RunArtifactV1 - Exportable computation artifacts
 *
 * This module defines the complete artifact schema for exporting computation
 * results. An artifact contains the input state (RunStateV1), computed results,
 * and metadata for reproducibility and grading.
 *
 * The artifact is designed to be:
 * - JSON-serializable for file export/import
 * - Validatable at runtime via Zod schema
 * - Version-aware for future schema migrations
 * - Self-documenting with all provenance information
 *
 * @module run-artifact
 * @see TDD Section 5.4 for specification
 */

import { z } from 'zod';
import type { RunStateV1 } from './run-state';
import { RunStateV1Schema } from './run-state.schema';

// ============================================================================
// Module-Specific Result Data Types
// ============================================================================

/**
 * Boys function computation result data
 *
 * Contains the computed value of F_m(T) along with diagnostic information
 * about the computational method used.
 */
export interface BoysArtifactData {
  /** Computed value F_m(T) */
  value: number;
  /** Computational method used: "zero", "series", or "recurrence" */
  method: string;
  /** Number of terms/steps used in computation */
  terms_count: number;
  /** Estimated relative error bound (null if not available) */
  estimated_error: number | null;
}

/**
 * Rys quadrature computation result data
 *
 * Contains the computed roots and weights for Gaussian quadrature,
 * along with optional error curve data for accuracy analysis.
 */
export interface RysArtifactData {
  /** Quadrature roots in the interval (0, 1) */
  roots: number[];
  /** Corresponding weights (all strictly positive) */
  weights: number[];
  /** Error curve data showing max error vs quadrature order */
  error_curve?: ErrorCurvePoint[];
  /** Maximum reconstruction error across all moments */
  reconstruction_error: number;
}

/**
 * A single point on the error curve
 */
export interface ErrorCurvePoint {
  /** Quadrature order (number of roots/weights) */
  n: number;
  /** Maximum absolute reconstruction error at this order */
  max_error: number;
}

/**
 * SCF computation result data
 *
 * Contains the final energy, convergence status, and full iteration
 * history for the Self-Consistent Field calculation.
 */
export interface ScfArtifactData {
  /** Final total electronic energy in Hartree */
  energy: number;
  /** Whether SCF achieved convergence */
  converged: boolean;
  /** Number of iterations performed */
  iterations: number;
  /** Whether computation was aborted by user */
  aborted: boolean;
  /** Full iteration history for plotting */
  trace: ScfIterationHistory[];
  /** Orbital energies in Hartree (if internals mode) */
  orbital_energies?: number[];
  /** HOMO-LUMO gap in Hartree (if internals mode) */
  homo_lumo_gap?: number;
}

/**
 * Single iteration in SCF history
 *
 * Records the energy and convergence metrics at each SCF iteration.
 */
export interface ScfIterationHistory {
  /** Iteration number (1-indexed) */
  iteration: number;
  /** Total electronic energy at this iteration (Hartree) */
  energy: number;
  /** Energy change from previous iteration */
  delta: number;
  /** DIIS error metric (if DIIS was used) */
  diis_error?: number;
}

// ============================================================================
// Discriminated Union for Results
// ============================================================================

/**
 * Boys module result wrapper
 */
export interface BoysArtifactResult {
  type: 'boys';
  data: BoysArtifactData;
}

/**
 * Rys module result wrapper
 */
export interface RysArtifactResult {
  type: 'rys';
  data: RysArtifactData;
}

/**
 * SCF module result wrapper
 */
export interface ScfArtifactResult {
  type: 'scf';
  data: ScfArtifactData;
}

/**
 * Discriminated union of all module-specific results
 *
 * Use the `type` field to discriminate between modules:
 *
 * @example
 * ```typescript
 * function processResults(results: ArtifactResults) {
 *   switch (results.type) {
 *     case 'boys':
 *       console.log('Boys value:', results.data.value);
 *       break;
 *     case 'rys':
 *       console.log('Rys roots:', results.data.roots);
 *       break;
 *     case 'scf':
 *       console.log('SCF energy:', results.data.energy);
 *       break;
 *   }
 * }
 * ```
 */
export type ArtifactResults =
  | BoysArtifactResult
  | RysArtifactResult
  | ScfArtifactResult;

// ============================================================================
// Artifact Metadata
// ============================================================================

/**
 * Metadata about the artifact and computation environment
 *
 * Provides provenance information for reproducibility verification.
 */
export interface ArtifactMetadata {
  /** ISO 8601 timestamp when artifact was created */
  timestamp: string;
  /** Application version that created this artifact */
  app_version: string;
  /** Version string from the WASM module */
  wasm_version: string;
  /** Total computation time in milliseconds */
  compute_time_ms: number;
  /** Unique identifier for this run (hash of state + results) */
  run_id: string;
}

/**
 * @deprecated Use ArtifactMetadata instead
 * Alias for backward compatibility
 */
export type ComputeMetadata = ArtifactMetadata;

// ============================================================================
// Main RunArtifactV1 Interface
// ============================================================================

/**
 * RunArtifactV1 - Complete exportable computation artifact
 *
 * This is the full artifact exported when users download their results.
 * It contains all information needed for:
 * - Reproducing the computation (via state)
 * - Grading and verification (via results)
 * - Provenance tracking (via metadata)
 *
 * @example
 * ```typescript
 * const artifact: RunArtifactV1 = {
 *   schema_version: 'run_artifact_v1',
 *   state: { ... },
 *   results: { type: 'scf', data: { ... } },
 *   metadata: {
 *     timestamp: '2026-01-17T14:30:45.123Z',
 *     app_version: '0.1.0',
 *     wasm_version: '0.1.0',
 *     compute_time_ms: 45,
 *     run_id: 'a1b2c3d4'
 *   }
 * };
 * ```
 */
export interface RunArtifactV1 {
  /** Schema version for migration support (always 'run_artifact_v1' for this version) */
  schema_version: 'run_artifact_v1';
  /** Input state used for the computation */
  state: RunStateV1;
  /** Computed results (discriminated by module type) */
  results: ArtifactResults;
  /** Metadata about the artifact and computation environment */
  metadata: ArtifactMetadata;
}

// ============================================================================
// Zod Validation Schemas
// ============================================================================

/**
 * Zod schema for ScfIterationHistory
 */
export const ScfIterationHistorySchema = z.object({
  iteration: z.number().int().min(1),
  energy: z.number(),
  delta: z.number(),
  diis_error: z.number().optional(),
});

/**
 * Zod schema for BoysArtifactData
 */
export const BoysArtifactDataSchema = z.object({
  value: z.number(),
  method: z.string().min(1),
  terms_count: z.number().int().min(0),
  estimated_error: z.number().nullable(),
});

/**
 * Zod schema for ErrorCurvePoint
 */
export const ErrorCurvePointSchema = z.object({
  n: z.number().int().min(1),
  max_error: z.number().min(0),
});

/**
 * Zod schema for RysArtifactData
 */
export const RysArtifactDataSchema = z.object({
  roots: z.array(z.number()),
  weights: z.array(z.number()),
  error_curve: z.array(ErrorCurvePointSchema).optional(),
  reconstruction_error: z.number().min(0),
});

/**
 * Zod schema for ScfArtifactData
 */
export const ScfArtifactDataSchema = z.object({
  energy: z.number(),
  converged: z.boolean(),
  iterations: z.number().int().min(0),
  aborted: z.boolean(),
  trace: z.array(ScfIterationHistorySchema),
  orbital_energies: z.array(z.number()).optional(),
  homo_lumo_gap: z.number().optional(),
});

/**
 * Zod schema for BoysArtifactResult
 */
export const BoysArtifactResultSchema = z.object({
  type: z.literal('boys'),
  data: BoysArtifactDataSchema,
});

/**
 * Zod schema for RysArtifactResult
 */
export const RysArtifactResultSchema = z.object({
  type: z.literal('rys'),
  data: RysArtifactDataSchema,
});

/**
 * Zod schema for ScfArtifactResult
 */
export const ScfArtifactResultSchema = z.object({
  type: z.literal('scf'),
  data: ScfArtifactDataSchema,
});

/**
 * Zod schema for ArtifactResults discriminated union
 */
export const ArtifactResultsSchema = z.discriminatedUnion('type', [
  BoysArtifactResultSchema,
  RysArtifactResultSchema,
  ScfArtifactResultSchema,
]);

/**
 * Zod schema for ArtifactMetadata
 */
export const ArtifactMetadataSchema = z.object({
  timestamp: z.string().datetime(),
  app_version: z.string().min(1).max(50),
  wasm_version: z.string().min(1),
  compute_time_ms: z.number().min(0),
  run_id: z.string().min(1).max(100),
});

/**
 * @deprecated Use ArtifactMetadataSchema instead
 * Alias for backward compatibility
 */
export const ComputeMetadataSchema = ArtifactMetadataSchema;

/**
 * Zod schema for RunArtifactV1
 *
 * Complete validation schema for imported artifacts.
 * Use `.safeParse()` for validation without throwing:
 *
 * @example
 * ```typescript
 * const result = RunArtifactV1Schema.safeParse(data);
 * if (result.success) {
 *   const artifact: RunArtifactV1 = result.data;
 *   // Use validated artifact
 * } else {
 *   console.error('Validation failed:', result.error.format());
 *   // Handle invalid artifact
 * }
 * ```
 */
export const RunArtifactV1Schema = z.object({
  schema_version: z.literal('run_artifact_v1'),
  state: RunStateV1Schema,
  results: ArtifactResultsSchema,
  metadata: ArtifactMetadataSchema,
});

// ============================================================================
// Validation Utilities
// ============================================================================

/**
 * Type guard for RunArtifactV1
 *
 * Checks if an unknown value is a valid RunArtifactV1 object.
 * Uses Zod schema validation internally.
 *
 * @param value - The value to check
 * @returns true if value is a valid RunArtifactV1
 *
 * @example
 * ```typescript
 * const data: unknown = JSON.parse(fileContents);
 * if (isRunArtifactV1(data)) {
 *   // data is now typed as RunArtifactV1
 *   console.log(data.results.type);
 * }
 * ```
 */
export function isRunArtifactV1(value: unknown): value is RunArtifactV1 {
  return RunArtifactV1Schema.safeParse(value).success;
}

/** Result type from validateRunArtifact when validation succeeds */
export interface ArtifactValidationSuccess {
  success: true;
  data: RunArtifactV1;
}

/** Result type from validateRunArtifact when validation fails */
export interface ArtifactValidationFailure {
  success: false;
  error: z.ZodError;
}

/** Union type for validateRunArtifact result */
export type ArtifactValidationResult =
  | ArtifactValidationSuccess
  | ArtifactValidationFailure;

/**
 * Validate a value against the RunArtifactV1 schema
 *
 * Returns a result object with either the validated data or error information.
 * This is a thin wrapper around Zod's safeParse for convenience.
 *
 * @param value - The value to validate
 * @returns Object with success boolean and either data or error
 *
 * @example
 * ```typescript
 * const result = validateRunArtifact(data);
 * if (result.success) {
 *   restoreArtifact(result.data);
 * } else {
 *   showError(result.error.format());
 * }
 * ```
 */
export function validateRunArtifact(value: unknown): ArtifactValidationResult {
  const result = RunArtifactV1Schema.safeParse(value);
  if (result.success) {
    return { success: true, data: result.data };
  }
  return { success: false, error: result.error };
}
