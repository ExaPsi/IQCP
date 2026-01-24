/**
 * Integral computation handlers.
 *
 * Handles on-the-fly computation of molecular integrals (overlap, kinetic,
 * nuclear attraction, and two-electron) for user-specified geometries.
 * Results are compatible with PresetSystemJson for seamless SCF integration.
 *
 * @module worker/handlers/integral
 * @see US-029 Integral WASM Integration
 */

import type {
  IntegralComputeRequest,
  IntegralComputeResult,
  WorkerProgress,
  ErrorResponse,
  IntegralPhase,
  GeometryInput,
  BasisSetName,
} from '../protocol';

/**
 * Progress callback type for integral computation updates.
 */
export type ProgressCallback = (progress: WorkerProgress) => void;

/**
 * Abort flag checker type.
 */
export type AbortChecker = () => boolean;

/**
 * WASM integral progress structure (matches IntegralWasmProgress in qc-wasm).
 *
 * Note: Field names use camelCase as they come from serde serialization.
 */
interface IntegralWasmProgress {
  /** Current phase: "overlap", "kinetic", "nuclear", "eri", "assembly" */
  phase: string;
  /** Current step within phase */
  current: number;
  /** Total steps in phase */
  total: number;
  /** Overall percentage complete (0-100) */
  overallPercent: number;
  /** Human-readable message */
  message: string;
}

/**
 * WASM atom output structure (matches AtomOutput in qc-wasm).
 */
interface WasmAtomOutput {
  symbol: string;
  xyz: [number, number, number];
  atomicNumber: number;
}

/**
 * WASM geometry output structure (matches GeometryOutput in qc-wasm).
 */
interface WasmGeometryOutput {
  atoms: WasmAtomOutput[];
  units: string;
}

/**
 * WASM integral metadata structure (matches IntegralMetadata in qc-wasm).
 */
interface WasmIntegralMetadata {
  wasmVersion: string;
  computeTimeMs: number;
  shellPairs: number;
  shellQuartets: number;
  significantEris: number;
  /** Basis type used: "cartesian" or "spherical" */
  basisType: string;
}

/**
 * WASM integral result structure (matches IntegralResult in qc-wasm).
 *
 * Note: Field names use camelCase as they come from serde serialization.
 */
interface IntegralWasmResult {
  formatVersion: number;
  systemId: string;
  label: string;
  description: string;
  geometry: WasmGeometryOutput;
  basisId: string;
  nbf: number;
  nelec: number;
  eNuc: number;
  sMatrix: number[];
  hCore: number[];
  eriCompressed: number[];
  eriIndexing: string;
  metadata: WasmIntegralMetadata;
}

/**
 * Options for WASM integral computation.
 * Matches IntegralComputeOptions in qc-wasm.
 */
interface IntegralComputeOptions {
  /** Use spherical harmonic basis functions (5 d-orbitals vs 6 Cartesian) */
  useSpherical?: boolean;
}

/**
 * Type for the WASM compute_integrals function.
 *
 * This is the function signature exported from qc-wasm.
 */
export type ComputeIntegralsFn = (
  geometryJson: string,
  basisName: string,
  progressCallback?: (progress: IntegralWasmProgress) => void
) => IntegralWasmResult;

/**
 * Type for the WASM compute_integrals_with_options function.
 *
 * This is the function signature for the version that accepts options.
 */
export type ComputeIntegralsWithOptionsFn = (
  geometryJson: string,
  basisName: string,
  options: IntegralComputeOptions
) => IntegralWasmResult;

/**
 * Type for the WASM compute_integrals_with_options_and_progress function.
 *
 * This is the function signature for the version that accepts both options
 * and a progress callback. This is the preferred function when available
 * as it supports both spherical harmonics and progress reporting.
 */
export type ComputeIntegralsWithOptionsAndProgressFn = (
  geometryJson: string,
  basisName: string,
  options: IntegralComputeOptions,
  progressCallback: (progress: IntegralWasmProgress) => void
) => IntegralWasmResult;

// Store reference to the WASM compute_integrals function (set during worker initialization)
let wasmComputeIntegrals: ComputeIntegralsFn | null = null;

// Store reference to the WASM compute_integrals_with_options function
let wasmComputeIntegralsWithOptions: ComputeIntegralsWithOptionsFn | null = null;

// Store reference to the WASM compute_integrals_with_options_and_progress function
let wasmComputeIntegralsWithOptionsAndProgress: ComputeIntegralsWithOptionsAndProgressFn | null =
  null;

/**
 * Set the WASM compute_integrals function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function is stored and used by handleIntegralCompute.
 *
 * @param fn - The WASM compute_integrals function
 */
export function setComputeIntegralsFn(fn: ComputeIntegralsFn): void {
  wasmComputeIntegrals = fn;
}

/**
 * Set the WASM compute_integrals_with_options function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function enables spherical harmonic basis support.
 *
 * @param fn - The WASM compute_integrals_with_options function
 */
export function setComputeIntegralsWithOptionsFn(fn: ComputeIntegralsWithOptionsFn): void {
  wasmComputeIntegralsWithOptions = fn;
}

/**
 * Set the WASM compute_integrals_with_options_and_progress function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. This function supports both spherical harmonics AND progress
 * callbacks, making it the preferred option when available.
 *
 * @param fn - The WASM compute_integrals_with_options_and_progress function
 */
export function setComputeIntegralsWithOptionsAndProgressFn(
  fn: ComputeIntegralsWithOptionsAndProgressFn
): void {
  wasmComputeIntegralsWithOptionsAndProgress = fn;
}

/**
 * Supported element symbols (H through Ne).
 */
const SUPPORTED_ELEMENTS = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'];

/**
 * Supported basis set names.
 */
const SUPPORTED_BASIS_SETS: BasisSetName[] = ['sto-3g', '3-21g', '6-31g', '6-31g*', '6-31+g*'];

/**
 * Validate geometry input.
 *
 * @param geometry - The geometry to validate
 * @returns Error message if invalid, null if valid
 */
function validateGeometry(geometry: GeometryInput): string | null {
  if (!geometry.atoms || geometry.atoms.length === 0) {
    return 'Geometry must have at least 1 atom.';
  }

  for (let i = 0; i < geometry.atoms.length; i++) {
    const atom = geometry.atoms[i];
    if (!SUPPORTED_ELEMENTS.includes(atom.symbol)) {
      return `Unsupported element '${atom.symbol}' at position ${i}. Only H-Ne are supported.`;
    }
    if (!Array.isArray(atom.xyz) || atom.xyz.length !== 3) {
      return `Invalid coordinates for atom ${i} (${atom.symbol}). Expected [x, y, z] array.`;
    }
    for (let j = 0; j < 3; j++) {
      if (typeof atom.xyz[j] !== 'number' || !Number.isFinite(atom.xyz[j])) {
        return `Invalid coordinate value for atom ${i} (${atom.symbol}) at index ${j}.`;
      }
    }
  }

  if (geometry.units !== 'bohr' && geometry.units !== 'angstrom') {
    return `Invalid units '${geometry.units}'. Must be 'bohr' or 'angstrom'.`;
  }

  return null;
}

/**
 * Validate basis set name.
 *
 * @param basisSet - The basis set name to validate
 * @returns Error message if invalid, null if valid
 */
function validateBasisSet(basisSet: string): string | null {
  if (!SUPPORTED_BASIS_SETS.includes(basisSet as BasisSetName)) {
    return `Unknown basis set '${basisSet}'. Supported: ${SUPPORTED_BASIS_SETS.join(', ')}`;
  }
  return null;
}

/**
 * Handle integral_compute request.
 *
 * Computes all molecular integrals (S, T, V, ERI) for a user-specified
 * geometry and basis set. Progress updates are emitted for each computation
 * phase, and the operation can be cancelled via the abort checker.
 *
 * The result is compatible with PresetSystemJson, enabling direct use
 * with the existing SCF workflow.
 *
 * @param request - The integral compute request
 * @param onProgress - Callback for progress updates
 * @param isAborted - Function to check if computation should abort
 * @returns Integral result or error response
 *
 * @example
 * ```typescript
 * const result = handleIntegralCompute(
 *   {
 *     type: 'integral_compute',
 *     requestId: id,
 *     geometry: { atoms: [{ symbol: 'H', xyz: [0, 0, 0] }, { symbol: 'H', xyz: [0, 0, 1.4] }], units: 'bohr' },
 *     basisSet: 'sto-3g'
 *   },
 *   (progress) => postMessage({ type: 'progress', requestId: id, progress }),
 *   () => abortFlags.get(id) === true
 * );
 * ```
 */
export function handleIntegralCompute(
  request: IntegralComputeRequest,
  onProgress: ProgressCallback,
  isAborted: AbortChecker
): IntegralComputeResult | ErrorResponse {
  const { requestId, geometry, basisSet, useSpherical } = request;

  // Check if WASM function is available
  // Priority: options+progress > progress-only (if not spherical) > options-only > basic
  // If useSpherical is requested, we need either options+progress or options version
  const hasOptionsAndProgress = wasmComputeIntegralsWithOptionsAndProgress !== null;
  const hasOptions = wasmComputeIntegralsWithOptions !== null;
  const hasProgress = wasmComputeIntegrals !== null;

  if (useSpherical && !hasOptionsAndProgress && !hasOptions) {
    return {
      type: 'error',
      requestId,
      code: 'WORKER_NOT_READY',
      message: 'Spherical harmonics require compute_integrals_with_options (not available)',
    };
  }
  if (!hasOptionsAndProgress && !hasOptions && !hasProgress) {
    return {
      type: 'error',
      requestId,
      code: 'WORKER_NOT_READY',
      message: 'Integral computation WASM function not initialized',
    };
  }

  // Validate geometry
  const geometryError = validateGeometry(geometry);
  if (geometryError) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: geometryError,
    };
  }

  // Validate basis set
  const basisError = validateBasisSet(basisSet);
  if (basisError) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: basisError,
    };
  }

  // Check abort before starting
  if (isAborted()) {
    return {
      type: 'error',
      requestId,
      code: 'COMPUTATION_CANCELLED',
      message: 'Computation cancelled before starting',
    };
  }

  // Set up progress callback that wraps WASM progress into worker progress
  let lastPhase: IntegralPhase = 'overlap';
  const wasmProgressCallback = (wasmProgress: IntegralWasmProgress): void => {
    // Check abort during computation
    if (isAborted()) {
      // Note: We can't actually stop WASM execution mid-computation,
      // but we can stop emitting progress updates
      return;
    }

    // Map WASM phase string to IntegralPhase type
    const phase = wasmProgress.phase as IntegralPhase;
    lastPhase = phase;

    onProgress({
      module: 'integral',
      phase,
      current: wasmProgress.current,
      total: wasmProgress.total,
      overallPercent: wasmProgress.overallPercent,
      message: wasmProgress.message,
    });
  };

  // Call WASM compute_integrals function
  // Priority order:
  //   1. compute_integrals_with_options_and_progress - best, has both options and progress
  //   2. compute_integrals_with_progress - if not using spherical (has progress, no options)
  //   3. compute_integrals_with_options - if spherical needed but no progress version
  //   4. compute_integrals - legacy fallback (no options, no progress)
  let result: IntegralWasmResult;
  try {
    // Serialize geometry to JSON for WASM
    const geometryJson = JSON.stringify(geometry);
    const options: IntegralComputeOptions = { useSpherical: useSpherical ?? false };

    if (wasmComputeIntegralsWithOptionsAndProgress) {
      // Best option: supports both options (spherical) AND progress callbacks
      result = wasmComputeIntegralsWithOptionsAndProgress(
        geometryJson,
        basisSet,
        options,
        wasmProgressCallback
      );
    } else if (!useSpherical && wasmComputeIntegrals) {
      // Second best for Cartesian: has progress callback support
      result = wasmComputeIntegrals(geometryJson, basisSet, wasmProgressCallback);
    } else if (wasmComputeIntegralsWithOptions) {
      // Fall back to options version (supports spherical, no progress)
      result = wasmComputeIntegralsWithOptions(geometryJson, basisSet, options);
    } else if (wasmComputeIntegrals) {
      // Legacy fallback (Cartesian only with progress)
      result = wasmComputeIntegrals(geometryJson, basisSet, wasmProgressCallback);
    } else {
      throw new Error('No WASM integral function available');
    }
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'Integral computation failed',
    };
  }

  // Check abort after computation (user may have cancelled during computation)
  if (isAborted()) {
    return {
      type: 'error',
      requestId,
      code: 'COMPUTATION_CANCELLED',
      message: `Computation cancelled during ${lastPhase} phase`,
    };
  }

  // Return WASM result directly - the structure already matches IntegralComputeResult
  // since serde serializes with camelCase. This avoids unnecessary object recreation
  // and copying of large arrays (sMatrix, hCore, eriCompressed).
  //
  // The only type coercion needed is for the string literal types.
  return {
    formatVersion: result.formatVersion,
    systemId: result.systemId,
    label: result.label,
    description: result.description,
    // Geometry: pass through directly, just type-narrow the units
    geometry: {
      atoms: result.geometry.atoms,
      units: result.geometry.units as 'bohr' | 'angstrom',
    },
    basisId: result.basisId,
    nbf: result.nbf,
    nelec: result.nelec,
    eNuc: result.eNuc,
    // Large arrays: pass through directly without copying
    sMatrix: result.sMatrix,
    hCore: result.hCore,
    eriCompressed: result.eriCompressed,
    eriIndexing: result.eriIndexing,
    // Metadata: pass through, just type-narrow basisType
    metadata: {
      wasmVersion: result.metadata.wasmVersion,
      computeTimeMs: result.metadata.computeTimeMs,
      shellPairs: result.metadata.shellPairs,
      shellQuartets: result.metadata.shellQuartets,
      significantEris: result.metadata.significantEris,
      basisType: (result.metadata.basisType as 'cartesian' | 'spherical') || 'cartesian',
    },
  };
}
