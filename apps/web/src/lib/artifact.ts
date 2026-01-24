/**
 * Run Artifact Export Utilities
 *
 * Provides functions for creating, identifying, and downloading run artifacts.
 * Artifacts capture a complete snapshot of a computation for sharing and
 * reproducibility.
 *
 * @module artifact
 */

import type {
  RunArtifactV1,
  ArtifactResults,
  ArtifactMetadata,
  BoysArtifactData,
  RysArtifactData,
  ScfArtifactData,
} from '../types/run-artifact';
import type { RunStateV1 } from '../types/run-state';
import { APP_VERSION } from '../types/run-state';
import type { BoysEvalResult, RysComputeResult, ScfRunResult, ScfIterationHistory } from '../worker/protocol';
import { useBoysStore } from '../stores/boysStore';
import { useRysStore, type TargetAccuracy } from '../stores/rysStore';
import { useScfStore } from '../stores/scfStore';

/**
 * DJB2 hash algorithm implementation.
 *
 * A simple, fast hash function suitable for generating short identifiers.
 * Not cryptographically secure, but produces good distribution for
 * string inputs.
 *
 * @param str - Input string to hash
 * @returns 32-bit unsigned integer hash
 */
function djb2Hash(str: string): number {
  let hash = 5381;
  for (let i = 0; i < str.length; i++) {
    // hash * 33 + char
    hash = (hash * 33) ^ str.charCodeAt(i);
  }
  // Convert to unsigned 32-bit integer
  return hash >>> 0;
}

/**
 * Generate a unique run ID from state and results.
 *
 * Uses a simple hash of the JSON representation of the state and results
 * to create a deterministic identifier. The same inputs will always
 * produce the same run ID.
 *
 * @param state - The run state containing parameters
 * @param results - The computation results
 * @returns Hexadecimal string (8 characters)
 *
 * @example
 * ```typescript
 * import { generateRunId } from '@/lib/artifact';
 *
 * const runId = generateRunId(state, results);
 * console.log(runId); // e.g., "a1b2c3d4"
 * ```
 */
export function generateRunId(state: RunStateV1, results: ArtifactResults): string {
  const combined = JSON.stringify({ state, results });
  const hash = djb2Hash(combined);
  return hash.toString(16).padStart(8, '0');
}

/**
 * Generate ISO 8601 timestamp for artifact metadata.
 *
 * @returns ISO 8601 formatted timestamp string
 */
function generateTimestamp(): string {
  return new Date().toISOString();
}

/**
 * Create a RunArtifactV1 from current state and results.
 *
 * Assembles a complete artifact with metadata including timestamps,
 * version information, and a unique run identifier.
 *
 * @param state - The complete run state
 * @param results - The computation results
 * @param options - Optional metadata overrides
 * @param options.computeTimeMs - Computation time in milliseconds
 * @param options.wasmVersion - WASM module version string
 * @returns Complete artifact ready for export
 *
 * @example
 * ```typescript
 * import { createArtifact } from '@/lib/artifact';
 *
 * const artifact = createArtifact(state, results, {
 *   computeTimeMs: 12.5,
 *   wasmVersion: '0.1.0',
 * });
 * ```
 */
export function createArtifact(
  state: RunStateV1,
  results: ArtifactResults,
  options?: {
    computeTimeMs?: number;
    wasmVersion?: string;
  }
): RunArtifactV1 {
  const timestamp = generateTimestamp();
  const runId = generateRunId(state, results);

  const metadata: ArtifactMetadata = {
    timestamp,
    app_version: APP_VERSION,
    wasm_version: options?.wasmVersion ?? 'unknown',
    compute_time_ms: options?.computeTimeMs ?? 0,
    run_id: runId,
  };

  return {
    schema_version: 'run_artifact_v1',
    state,
    results,
    metadata,
  };
}

/**
 * Format a date to compact timestamp string.
 *
 * Produces format: 20260117T143045Z (compact ISO 8601 basic format)
 *
 * @param date - Date to format
 * @returns Compact timestamp string
 */
function formatCompactTimestamp(date: Date): string {
  const year = date.getUTCFullYear();
  const month = String(date.getUTCMonth() + 1).padStart(2, '0');
  const day = String(date.getUTCDate()).padStart(2, '0');
  const hours = String(date.getUTCHours()).padStart(2, '0');
  const minutes = String(date.getUTCMinutes()).padStart(2, '0');
  const seconds = String(date.getUTCSeconds()).padStart(2, '0');

  return `${year}${month}${day}T${hours}${minutes}${seconds}Z`;
}

/**
 * Generate a descriptive filename for artifact export.
 *
 * Format varies by module:
 * - Boys: iqcp_boys_m{m}_T{T}_20260117T143045Z.json
 * - Rys: iqcp_rys_n{n}_T{T}_20260117T143045Z.json
 * - SCF: iqcp_scf_{systemId}_{conv}_{diis}_20260117T143045Z.json
 *
 * @param artifact - The artifact to generate filename for
 * @returns Descriptive filename with .json extension
 *
 * @example
 * ```typescript
 * import { generateFilename } from '@/lib/artifact';
 *
 * const filename = generateFilename(artifact);
 * // "iqcp_boys_m3_T5.00_20260117T143045Z.json"
 * ```
 */
export function generateFilename(artifact: RunArtifactV1): string {
  const { state, metadata } = artifact;
  const module = state.module;

  // Parse timestamp from ISO 8601 format
  const timestamp = formatCompactTimestamp(new Date(metadata.timestamp));

  let params: string;

  switch (module) {
    case 'boys': {
      const boysParams = state.boys;
      if (boysParams) {
        const m = boysParams.m;
        const T = boysParams.T.toFixed(2);
        params = `m${m}_T${T}`;
      } else {
        params = 'unknown';
      }
      break;
    }

    case 'rys': {
      const rysParams = state.rys;
      if (rysParams) {
        const n = rysParams.n;
        const T = rysParams.T.toFixed(2);
        params = `n${n}_T${T}`;
      } else {
        params = 'unknown';
      }
      break;
    }

    case 'scf': {
      const scfParams = state.scf;
      if (scfParams) {
        const systemId = scfParams.system_id;
        const conv = scfParams.conv;
        const diis = scfParams.diis ? 'diis' : 'nodiis';
        params = `${systemId}_${conv}_${diis}`;
      } else {
        params = 'unknown';
      }
      break;
    }

    default:
      params = 'unknown';
  }

  return `iqcp_${module}_${params}_${timestamp}.json`;
}

/**
 * Download an artifact as a JSON file using the browser Blob API.
 *
 * Creates a temporary anchor element to trigger a download without
 * any external dependencies. The file is formatted with 2-space
 * indentation for human readability.
 *
 * @param artifact - The artifact to download
 *
 * @example
 * ```typescript
 * import { downloadArtifact, createArtifact } from '@/lib/artifact';
 *
 * const artifact = createArtifact(state, results);
 * downloadArtifact(artifact);
 * // Browser will prompt to save the file
 * ```
 */
export function downloadArtifact(artifact: RunArtifactV1): void {
  // Serialize to formatted JSON
  const json = JSON.stringify(artifact, null, 2);

  // Create Blob with JSON MIME type
  const blob = new Blob([json], { type: 'application/json' });

  // Create object URL for the blob
  const url = URL.createObjectURL(blob);

  // Generate filename
  const filename = generateFilename(artifact);

  // Create temporary anchor element
  const anchor = document.createElement('a');
  anchor.href = url;
  anchor.download = filename;

  // Append to DOM, click, and remove
  document.body.appendChild(anchor);
  anchor.click();
  document.body.removeChild(anchor);

  // Clean up object URL
  URL.revokeObjectURL(url);
}

// ============================================================================
// Artifact Import (Restoration)
// ============================================================================

/**
 * Result from restoreArtifact function.
 */
export interface RestoreResult {
  /** Whether restoration was successful */
  success: boolean;
  /** Any warnings encountered during restoration */
  warnings: string[];
}

/**
 * Convert Boys artifact data to BoysEvalResult format.
 *
 * Maps snake_case artifact fields to camelCase store format.
 *
 * @param data - Artifact data in snake_case format
 * @param params - Boys parameters from artifact state
 * @returns BoysEvalResult compatible with store
 */
function convertBoysArtifactToResult(
  data: BoysArtifactData,
  params: { m: number; T: number }
): BoysEvalResult {
  return {
    value: data.value,
    method: data.method as BoysEvalResult['method'],
    m: params.m,
    t: params.T,
    termsCount: data.terms_count,
    estimatedError: data.estimated_error,
  };
}

/**
 * Convert Rys artifact data to RysComputeResult format.
 *
 * Maps snake_case artifact fields to camelCase store format.
 *
 * @param data - Artifact data in snake_case format
 * @param params - Rys parameters from artifact state
 * @returns RysComputeResult compatible with store
 */
function convertRysArtifactToResult(
  data: RysArtifactData,
  params: { n: number; T: number }
): RysComputeResult {
  return {
    roots: data.roots,
    weights: data.weights,
    nroots: data.roots.length,
    t: params.T,
    // Artifact doesn't store method, default to 'standard'
    method: params.T === 0 || params.n === 1 ? 'special' : 'standard',
  };
}

/**
 * Convert SCF artifact data to ScfRunResult format.
 *
 * Maps snake_case artifact fields to camelCase store format.
 *
 * @param data - Artifact data in snake_case format
 * @returns ScfRunResult compatible with store
 */
function convertScfArtifactToResult(data: ScfArtifactData): ScfRunResult {
  // Convert trace (snake_case) to history (camelCase)
  const history: ScfIterationHistory[] = data.trace.map((item) => ({
    iteration: item.iteration,
    energy: item.energy,
    delta: item.delta,
    diisError: item.diis_error,
  }));

  return {
    energy: data.energy,
    converged: data.converged,
    iterations: data.iterations,
    aborted: data.aborted,
    history,
    // Note: orbital energies and matrices are optional and may not be in artifact
    orbitalEnergies: data.orbital_energies
      ? {
          energies: data.orbital_energies,
          // We don't store nOccupied in artifact, would need to recalculate
          // For now, estimate as half the electrons (RHF)
          nOccupied: Math.floor(data.orbital_energies.length / 2),
        }
      : undefined,
  };
}

/**
 * Restore application state from an imported artifact.
 *
 * This function restores both the input parameters and computed results
 * from an artifact, enabling session restoration without re-computation.
 *
 * The function:
 * 1. Identifies the module type from the artifact
 * 2. Restores the module-specific parameters via store actions
 * 3. Restores the computed results to the store
 *
 * @param artifact - Validated RunArtifactV1 object
 * @returns Object indicating success and any warnings
 *
 * @example
 * ```typescript
 * import { restoreArtifact } from '@/lib/artifact';
 * import { validateRunArtifact } from '@/types/run-artifact';
 *
 * const result = validateRunArtifact(jsonData);
 * if (result.success) {
 *   const { success, warnings } = restoreArtifact(result.data);
 *   if (success) {
 *     console.log('Artifact restored!');
 *     warnings.forEach(w => console.warn(w));
 *   }
 * }
 * ```
 */
export function restoreArtifact(artifact: RunArtifactV1): RestoreResult {
  const warnings: string[] = [];

  // Check schema version
  if (artifact.schema_version !== 'run_artifact_v1') {
    warnings.push(
      `Schema version mismatch: expected 'run_artifact_v1', got '${artifact.schema_version}'`
    );
  }

  const { state, results } = artifact;
  const module = state.module;

  // Verify module type matches results type
  if (module !== results.type) {
    warnings.push(
      `Module mismatch: state.module='${module}' but results.type='${results.type}'`
    );
  }

  try {
    switch (results.type) {
      case 'boys': {
        const boysParams = state.boys;
        if (!boysParams) {
          return {
            success: false,
            warnings: [...warnings, 'Missing boys parameters in artifact state'],
          };
        }

        // Restore parameters
        const boysStore = useBoysStore.getState();
        boysStore.setM(boysParams.m);
        boysStore.setT(boysParams.T);
        boysStore.setView(boysParams.view);

        // Restore UI mode if present
        if (state.ui?.mode) {
          boysStore.setMode(state.ui.mode);
        }

        // Restore results
        const boysResult = convertBoysArtifactToResult(results.data, boysParams);
        boysStore.setComputeResult(boysResult);

        break;
      }

      case 'rys': {
        const rysParams = state.rys;
        if (!rysParams) {
          return {
            success: false,
            warnings: [...warnings, 'Missing rys parameters in artifact state'],
          };
        }

        // Restore parameters
        const rysStore = useRysStore.getState();
        rysStore.setN(rysParams.n);
        rysStore.setT(rysParams.T);
        rysStore.setTarget(rysParams.target as TargetAccuracy);

        // Restore UI mode if present
        if (state.ui?.mode) {
          rysStore.setMode(state.ui.mode);
        }

        // Restore results
        const rysResult = convertRysArtifactToResult(results.data, rysParams);
        rysStore.setComputeResult(rysResult);

        // Restore error curve if present in artifact
        if (results.data.error_curve && results.data.error_curve.length > 0) {
          rysStore.setErrorCurveResult({
            t: rysParams.T,
            nMax: results.data.error_curve.length,
            points: results.data.error_curve.map((p) => ({
              n: p.n,
              maxError: p.max_error,
            })),
          });
        }

        break;
      }

      case 'scf': {
        const scfParams = state.scf;
        if (!scfParams) {
          return {
            success: false,
            warnings: [...warnings, 'Missing scf parameters in artifact state'],
          };
        }

        // Restore parameters
        const scfStore = useScfStore.getState();
        scfStore.setSystemId(scfParams.system_id);
        scfStore.setConvergenceProfile(scfParams.conv);
        scfStore.setMaxIterations(scfParams.max_iter);
        scfStore.setUseDiis(scfParams.diis);

        // Restore UI mode if present
        if (state.ui?.mode) {
          scfStore.setMode(state.ui.mode);
        }

        // Restore results
        const scfResult = convertScfArtifactToResult(results.data);
        scfStore.setRunResult(scfResult);

        break;
      }

      default: {
        // Should not happen with validated artifact, but handle gracefully
        const exhaustiveCheck: never = results;
        return {
          success: false,
          warnings: [
            ...warnings,
            `Unknown results type: ${(exhaustiveCheck as ArtifactResults).type}`,
          ],
        };
      }
    }

    return { success: true, warnings };
  } catch (error) {
    const errorMessage =
      error instanceof Error ? error.message : 'Unknown error during restore';
    return {
      success: false,
      warnings: [...warnings, errorMessage],
    };
  }
}

/**
 * Get the expected module for an artifact.
 *
 * Useful for checking if user is on the correct module page before import.
 *
 * @param artifact - Validated artifact
 * @returns Module name ('boys', 'rys', or 'scf')
 */
export function getArtifactModule(artifact: RunArtifactV1): 'boys' | 'rys' | 'scf' {
  return artifact.results.type;
}
