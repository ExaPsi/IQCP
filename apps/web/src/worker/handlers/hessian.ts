/**
 * Frequency analysis handler.
 *
 * Handles `frequency` requests by calling the WASM `compute_frequencies`
 * function. Streams phase-granular progress events and supports
 * cancellation via an abort flag (cooperative — the WASM function cannot
 * interrupt itself mid-phase; see the `aborted` semantics in
 * `FrequencyResult`).
 *
 * Follows the structural pattern established by
 * `handlers/pesInternal.ts` and `handlers/scf.ts`.
 *
 * @module worker/handlers/hessian
 * @see US-101 Frequency WASM Export + Worker Handler
 */

import type {
  FrequencyRequest,
  FrequencyResult,
  FrequencyPhase,
  WorkerProgress,
  ErrorResponse,
} from '../protocol';

/**
 * Progress callback type for frequency analysis updates.
 */
export type ProgressCallback = (progress: WorkerProgress) => void;

/**
 * Abort flag checker type.
 */
export type AbortChecker = () => boolean;

/**
 * WASM progress structure emitted by `compute_frequencies`.
 *
 * Built on the Rust side via `js_sys::Reflect::set` as a plain JS object
 * with these four camelCase fields. Mirrors the `emit(...)` helper in
 * `run_frequency_pipeline`.
 */
interface FrequencyWasmProgress {
  /** Current pipeline phase */
  phase: FrequencyPhase;
  /** Fractional completion within the current phase (0.0 – 1.0) */
  percent: number;
  /** Sub-step identifier (e.g., "start", "done", "harmonic", "ir_sim") */
  step: string;
  /** Human-readable message */
  message: string;
}

/**
 * Type for the WASM `compute_frequencies` function.
 *
 * The second argument is an optional progress callback. The return value
 * is typed as `unknown` because `serde_wasm_bindgen` produces a plain
 * JavaScript object that we narrow to `FrequencyResult` after the call.
 */
export type ComputeFrequenciesFn = (
  input: unknown,
  progressCallback?: (progress: FrequencyWasmProgress) => void,
) => unknown;

// Store reference to the WASM compute_frequencies function
// (set during worker initialization).
let wasmComputeFrequencies: ComputeFrequenciesFn | null = null;

/**
 * Set the WASM `compute_frequencies` function reference.
 *
 * This must be called during worker initialization after the WASM module
 * is loaded. The function is stored and used by `handleFrequency`.
 *
 * @param fn - The WASM `compute_frequencies` function
 */
export function setComputeFrequenciesFn(fn: ComputeFrequenciesFn): void {
  wasmComputeFrequencies = fn;
}

/** Valid electronic-structure methods accepted by `compute_frequencies`. */
const VALID_METHODS = ['rhf', 'hf', 'lda', 'b3lyp', 'b3lyp-d3bj'] as const;

/** Valid broadening kinds accepted by `compute_frequencies`. */
const VALID_BROADENING_KINDS = ['lorentzian', 'gaussian'] as const;

/** Valid SCF convergence profiles. */
const VALID_CONVERGENCE_PROFILES = ['loose', 'medium', 'tight'] as const;

/**
 * Handle a frequency analysis request.
 *
 * Runs the full analytical vibrational spectroscopy pipeline by calling
 * the WASM `compute_frequencies` function. Streams progress for each
 * pipeline phase and supports cancellation via the abort checker.
 *
 * Cancellation is cooperative: the WASM function cannot interrupt itself
 * mid-phase. If the abort flag is set during the run, we stop forwarding
 * progress events to the main thread and mark the final result with
 * `aborted: true` once the WASM call returns. Results computed up to the
 * abort point are preserved.
 *
 * @param request - The frequency analysis request
 * @param onProgress - Callback for progress updates
 * @param isAborted - Function to check if computation should abort
 * @returns `FrequencyResult` on success, `ErrorResponse` on failure
 *
 * @example
 * ```typescript
 * const result = handleFrequency(
 *   { type: 'frequency', requestId: id, atoms: [...], basisName: 'sto-3g', method: 'rhf' },
 *   (progress) => postMessage({ type: 'progress', requestId: id, progress }),
 *   () => abortFlags.get(id) === true,
 * );
 * ```
 */
export function handleFrequency(
  request: FrequencyRequest,
  onProgress: ProgressCallback,
  isAborted: AbortChecker,
): FrequencyResult | ErrorResponse {
  const {
    requestId,
    atoms,
    basisName,
    method,
    temperatureK,
    pressurePa,
    symmetryNumberOverride,
    multiplicity,
    broadeningKind,
    fwhmCm1,
    convergenceProfile,
    maxIterations,
  } = request;

  // ---- Availability check -------------------------------------------------
  if (!wasmComputeFrequencies) {
    return {
      type: 'error',
      requestId,
      code: 'NOT_IMPLEMENTED',
      message:
        'compute_frequencies WASM function not available (rebuild WASM module)',
    };
  }

  // ---- Input validation ---------------------------------------------------
  if (!atoms || atoms.length === 0) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: 'FrequencyRequest: atoms array must not be empty',
    };
  }

  if (!basisName || basisName.trim() === '') {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: 'FrequencyRequest: basisName is required',
    };
  }

  const methodLower = method.toLowerCase();
  if (!(VALID_METHODS as readonly string[]).includes(methodLower)) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `FrequencyRequest: method must be one of ${VALID_METHODS.join(', ')}, got '${method}'`,
    };
  }

  if (temperatureK !== undefined && !(Number.isFinite(temperatureK) && temperatureK > 0)) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `FrequencyRequest: temperatureK must be positive and finite, got ${temperatureK}`,
    };
  }

  if (pressurePa !== undefined && !(Number.isFinite(pressurePa) && pressurePa > 0)) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `FrequencyRequest: pressurePa must be positive and finite, got ${pressurePa}`,
    };
  }

  if (fwhmCm1 !== undefined && !(Number.isFinite(fwhmCm1) && fwhmCm1 > 0)) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `FrequencyRequest: fwhmCm1 must be positive and finite, got ${fwhmCm1}`,
    };
  }

  if (multiplicity !== undefined && !(Number.isInteger(multiplicity) && multiplicity >= 1)) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `FrequencyRequest: multiplicity must be an integer >= 1, got ${multiplicity}`,
    };
  }

  if (
    symmetryNumberOverride !== undefined &&
    !(Number.isInteger(symmetryNumberOverride) && symmetryNumberOverride >= 1)
  ) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `FrequencyRequest: symmetryNumberOverride must be an integer >= 1, got ${symmetryNumberOverride}`,
    };
  }

  if (
    broadeningKind !== undefined &&
    !(VALID_BROADENING_KINDS as readonly string[]).includes(broadeningKind)
  ) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `FrequencyRequest: broadeningKind must be one of ${VALID_BROADENING_KINDS.join(', ')}, got '${broadeningKind}'`,
    };
  }

  if (
    convergenceProfile !== undefined &&
    !(VALID_CONVERGENCE_PROFILES as readonly string[]).includes(convergenceProfile)
  ) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `FrequencyRequest: convergenceProfile must be one of ${VALID_CONVERGENCE_PROFILES.join(', ')}, got '${convergenceProfile}'`,
    };
  }

  if (
    maxIterations !== undefined &&
    !(Number.isInteger(maxIterations) && maxIterations >= 1)
  ) {
    return {
      type: 'error',
      requestId,
      code: 'INVALID_PARAMS',
      message: `FrequencyRequest: maxIterations must be a positive integer, got ${maxIterations}`,
    };
  }

  // ---- Early abort check --------------------------------------------------
  if (isAborted()) {
    return buildAbortedPlaceholder(atoms.length);
  }

  // ---- Build WASM input (camelCase for serde deserialization) -------------
  // Note: omit optional fields when undefined so that Rust serde defaults
  // apply (default_freq_temperature_k, default_freq_fwhm_cm1, etc.).
  const wasmInput: Record<string, unknown> = {
    atoms,
    basisName,
    method: methodLower,
  };
  if (temperatureK !== undefined) wasmInput.temperatureK = temperatureK;
  if (pressurePa !== undefined) wasmInput.pressurePa = pressurePa;
  if (symmetryNumberOverride !== undefined)
    wasmInput.symmetryNumberOverride = symmetryNumberOverride;
  if (multiplicity !== undefined) wasmInput.multiplicity = multiplicity;
  if (broadeningKind !== undefined) wasmInput.broadeningKind = broadeningKind;
  if (fwhmCm1 !== undefined) wasmInput.fwhmCm1 = fwhmCm1;
  if (convergenceProfile !== undefined)
    wasmInput.convergenceProfile = convergenceProfile;
  if (maxIterations !== undefined) wasmInput.maxIterations = maxIterations;

  // ---- Progress translator + cooperative abort tracking ------------------
  // The WASM function cannot interrupt itself. We track whether an abort
  // was observed during any progress event and, if so, mark the returned
  // result with `aborted: true` once the WASM call completes.
  let abortedDuringProgress = false;

  const wasmProgressCallback = (wasmProgress: FrequencyWasmProgress): void => {
    if (isAborted()) {
      abortedDuringProgress = true;
      return; // Stop emitting further events; WASM keeps running.
    }

    const percent = Number.isFinite(wasmProgress.percent) ? wasmProgress.percent : 0;
    onProgress({
      module: 'frequency',
      phase: wasmProgress.phase,
      percent,
      step: wasmProgress.step ?? '',
      message: wasmProgress.message ?? '',
      current: Math.round(percent * 100),
      total: 100,
    });
  };

  // ---- Invoke WASM --------------------------------------------------------
  const startTime = performance.now();
  let rawResult: unknown;
  try {
    rawResult = wasmComputeFrequencies(wasmInput, wasmProgressCallback);
  } catch (e) {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: e instanceof Error ? e.message : 'Frequency analysis failed',
    };
  }
  const elapsedMs = performance.now() - startTime;

  if (!rawResult || typeof rawResult !== 'object') {
    return {
      type: 'error',
      requestId,
      code: 'HANDLER_ERROR',
      message: 'compute_frequencies returned an unexpected value',
    };
  }

  // ---- Tag with handler-side metadata ------------------------------------
  // The WASM function always returns aborted=false; the handler overrides
  // it when the abort flag was observed during a progress event. We also
  // patch timingMs.totalMs with the handler-measured elapsed time to
  // include any JS-side overhead.
  const typedResult = rawResult as FrequencyResult;
  const patchedTiming: FrequencyResult['timingMs'] = {
    ...typedResult.timingMs,
    totalMs: elapsedMs,
  };

  return {
    ...typedResult,
    timingMs: patchedTiming,
    aborted: abortedDuringProgress,
  };
}

/**
 * Build a minimal aborted-placeholder result for the early-abort case
 * (abort requested before we ever called the WASM function).
 *
 * Mirrors the zeroed-out `FrequencyResult` shape with empty arrays and
 * zero totals so that downstream consumers can still destructure safely.
 */
function buildAbortedPlaceholder(nAtoms: number): FrequencyResult {
  const zeroTiming: FrequencyResult['timingMs'] = {
    integralsMs: 0,
    nuclearCphfMs: 0,
    fieldCphfMs: 0,
    assemblyMs: 0,
    modesMs: 0,
    totalMs: 0,
  };

  const zeroThermo: FrequencyResult['thermochemistry'] = {
    temperatureK: 0,
    pressurePa: 0,
    symmetryNumber: 0,
    multiplicity: 0,
    totalMassAmu: 0,
    nVibModesUsed: 0,
    nImag: 0,
    zpeHa: 0,
    e0kHa: 0,
    internalEnergyHa: 0,
    enthalpyHa: 0,
    entropyHaPerK: 0,
    gibbsHa: 0,
    cvHaPerK: 0,
    cpHaPerK: 0,
    eTransHa: 0,
    hTransHa: 0,
    sTransHaPerK: 0,
    cvTransHaPerK: 0,
    cpTransHaPerK: 0,
    eRotHa: 0,
    hRotHa: 0,
    sRotHaPerK: 0,
    cvRotHaPerK: 0,
    cpRotHaPerK: 0,
    eVibThermalHa: 0,
    hVibHa: 0,
    sVibHaPerK: 0,
    cvVibHaPerK: 0,
    cpVibHaPerK: 0,
    sElecHaPerK: 0,
  };

  const zeroSpectrum: FrequencyResult['irSpectrum'] = {
    wavenumbersCm1: [],
    intensity: [],
    kind: 'lorentzian',
    fwhmCm1: 0,
  };

  return {
    nAtoms,
    nModes: 0,
    rotorType: 'atom',
    electronicEnergyHa: 0,
    dipoleAu: [0, 0, 0],
    dipoleDebye: [0, 0, 0],
    polarizabilityAu: [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ],
    polarizabilityAng3: [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ],
    frequenciesCm1: [],
    reducedMassesAmu: [],
    forceConstantsMdyne: [],
    normalModesCartesian: [],
    rotationalConstantsGhz: [0, 0, 0],
    irIntensitiesKmPerMol: [],
    ramanActivitiesA4Amu: [],
    depolarizationRatios: [],
    thermochemistry: zeroThermo,
    irSpectrum: zeroSpectrum,
    ramanSpectrum: { ...zeroSpectrum },
    timingMs: zeroTiming,
    aborted: true,
  };
}
