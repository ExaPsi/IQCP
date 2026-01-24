/**
 * Type Definitions Index
 *
 * Re-exports all public types from the types module for convenient imports.
 *
 * @module types
 */

// ============================================================================
// Run State Types
// ============================================================================

export {
  // Core types
  type ModuleId,
  type BoysParams,
  type RysParams,
  type ScfParams,
  type PresetRef,
  type UiSettings,
  type RunStateV1,
  // Constants
  APP_VERSION,
  // Default states
  DEFAULT_BOYS_STATE,
  DEFAULT_RYS_STATE,
  DEFAULT_SCF_STATE,
  // Utilities
  getDefaultState,
} from './run-state';

export {
  // Schemas
  ModuleIdSchema,
  BoysParamsSchema,
  RysParamsSchema,
  ScfParamsSchema,
  PresetRefSchema,
  UiSettingsSchema,
  RunStateV1Schema,
  // Type guards and validation
  isRunStateV1,
  validateRunState,
  // Validation result types
  type ValidationSuccess,
  type ValidationFailure,
  type ValidationResult,
} from './run-state.schema';

// ============================================================================
// Run Artifact Types
// ============================================================================

export {
  // Module-specific data types
  type BoysArtifactData,
  type RysArtifactData,
  type ErrorCurvePoint,
  type ScfArtifactData,
  type ScfIterationHistory,
  // Result wrapper types
  type BoysArtifactResult,
  type RysArtifactResult,
  type ScfArtifactResult,
  type ArtifactResults,
  // Metadata types
  type ArtifactMetadata,
  type ComputeMetadata,
  // Main artifact type
  type RunArtifactV1,
  // Schemas
  ScfIterationHistorySchema,
  BoysArtifactDataSchema,
  ErrorCurvePointSchema,
  RysArtifactDataSchema,
  ScfArtifactDataSchema,
  BoysArtifactResultSchema,
  RysArtifactResultSchema,
  ScfArtifactResultSchema,
  ArtifactResultsSchema,
  ArtifactMetadataSchema,
  ComputeMetadataSchema,
  RunArtifactV1Schema,
  // Type guards and validation
  isRunArtifactV1,
  validateRunArtifact,
  // Validation result types
  type ArtifactValidationSuccess,
  type ArtifactValidationFailure,
  type ArtifactValidationResult,
} from './run-artifact';
