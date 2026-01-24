/**
 * Web Worker Message Protocol
 *
 * Defines the typed message protocol for communication between
 * the main thread (React UI) and the compute Web Worker.
 *
 * All messages use discriminated unions with a `type` field for
 * type-safe switching and runtime validation.
 *
 * @module worker/protocol
 */

// ============================================================================
// Request ID
// ============================================================================

/**
 * Branded type for request IDs to prevent mixing with other strings.
 *
 * Using a branded type ensures type safety at compile time while
 * still being a string at runtime for serialization.
 */
export type RequestId = string & { readonly __brand: 'RequestId' };

/**
 * Generate a unique request ID.
 *
 * Uses `crypto.randomUUID()` for cryptographically random IDs
 * that are globally unique without coordination.
 *
 * @returns A branded RequestId
 *
 * @example
 * ```typescript
 * const id = createRequestId();
 * // id is typed as RequestId, not just string
 * ```
 */
export function createRequestId(): RequestId {
  return crypto.randomUUID() as RequestId;
}

// ============================================================================
// Geometry Input Types (for integral_compute request)
// ============================================================================

/**
 * Supported basis sets for on-the-fly integral computation.
 *
 * Currently supports minimal and small split-valence basis sets
 * suitable for educational use.
 *
 * - sto-3g: Minimal basis (3 Gaussians fit to Slater)
 * - 3-21g: Split-valence
 * - 6-31g: Split-valence (larger)
 * - 6-31g*: Split-valence + d polarization on heavy atoms
 * - 6-31+g*: Split-valence + diffuse sp + d polarization on heavy atoms
 */
export type BasisSetName = 'sto-3g' | '3-21g' | '6-31g' | '6-31g*' | '6-31+g*';

/**
 * Coordinate unit system.
 *
 * - "bohr": Atomic units (1 bohr = 0.529177 Angstrom)
 * - "angstrom": SI-derived (1 Angstrom = 1e-10 m)
 */
export type CoordinateUnits = 'bohr' | 'angstrom';

/**
 * Single atom specification for geometry input.
 *
 * @see IntegralComputeRequest
 */
export interface AtomInput {
  /** Element symbol (H-Ne supported) */
  symbol: string;
  /** Cartesian coordinates [x, y, z] */
  xyz: [number, number, number];
}

/**
 * Input geometry for integral computation.
 *
 * @see IntegralComputeRequest
 */
export interface GeometryInput {
  /** List of atoms in the molecule */
  atoms: AtomInput[];
  /** Coordinate units: "bohr" or "angstrom" */
  units: CoordinateUnits;
}

// ============================================================================
// SCF Options (for scf_run request)
// ============================================================================

/**
 * SCF convergence profile options.
 */
export type ConvergenceProfile = 'loose' | 'medium' | 'tight';

/**
 * SCF computation options.
 */
export interface ScfOptions {
  /** Convergence profile (loose/medium/tight) */
  convergenceProfile: ConvergenceProfile;
  /** Maximum number of iterations */
  maxIterations: number;
  /** Enable DIIS acceleration */
  useDiis: boolean;
  /** DIIS subspace size (if enabled) */
  diisSize?: number;
  /**
   * Fock matrix damping factor (0.0 = no damping, 0.5 = 50% old Fock).
   *
   * Damping mixes the current and previous Fock matrices to stabilize
   * convergence for difficult systems (e.g., diffuse basis sets like 6-31+G*).
   *
   * F_damped = damp * F_old + (1.0 - damp) * F_new
   *
   * Reference: PySCF hf.py damping implementation
   * @default 0.0
   */
  damp?: number;
  /**
   * Whether to include matrices in result (for internals mode).
   * When true, returns S, H_core, F, D matrices and orbital energies.
   * @default false
   */
  includeMatrices?: boolean;
}

// ============================================================================
// Worker Requests
// ============================================================================

/**
 * Base fields present in all requests.
 */
interface BaseRequest {
  /** Unique identifier for correlating responses */
  requestId: RequestId;
}

/**
 * Ping request to verify worker is alive and get WASM version.
 */
export interface PingRequest extends BaseRequest {
  type: 'ping';
}

/**
 * Boys function evaluation at a single point.
 *
 * @see Module B - Boys Function Lab
 */
export interface BoysEvalRequest extends BaseRequest {
  type: 'boys_eval';
  /** Order m of the Boys function F_m(T) */
  m: number;
  /** Argument T (must be >= 0) */
  T: number;
}

/**
 * Boys function sweep over a range of T values.
 *
 * Used for plotting F_m(T) curves.
 */
export interface BoysSweepRequest extends BaseRequest {
  type: 'boys_sweep';
  /** Order m of the Boys function */
  m: number;
  /** Range of T values [min, max] */
  T_range: [number, number];
  /** Number of points to compute */
  points: number;
}

/**
 * Rys quadrature roots and weights computation.
 *
 * @see Module C - Rys Quadrature Lab
 */
export interface RysComputeRequest extends BaseRequest {
  type: 'rys_compute';
  /** Number of quadrature points (order) */
  n: number;
  /** Parameter T for the Boys function moments */
  T: number;
}

/**
 * Rys error curve computation for order selection.
 *
 * Computes reconstruction error vs quadrature order.
 */
export interface RysErrorCurveRequest extends BaseRequest {
  type: 'rys_error_curve';
  /** Parameter T */
  T: number;
  /** Maximum quadrature order to test */
  max_order: number;
}

/**
 * SCF (Self-Consistent Field) computation request.
 *
 * @see Module E - SCF Sandbox
 */
export interface ScfRunRequest extends BaseRequest {
  type: 'scf_run';
  /** Identifier for the molecular system preset */
  systemId: string;
  /** SCF computation options */
  options: ScfOptions;
}

/**
 * Integral computation request.
 *
 * Computes overlap (S), kinetic (T), nuclear (V), and two-electron (ERI)
 * integrals for a user-specified geometry and basis set.
 *
 * The result is compatible with PresetSystemJson, allowing seamless
 * integration with the existing SCF workflow.
 *
 * @see US-029 Integral WASM Integration
 */
export interface IntegralComputeRequest extends BaseRequest {
  type: 'integral_compute';
  /** Geometry specification with atoms and coordinate units */
  geometry: GeometryInput;
  /** Basis set name for the computation */
  basisSet: BasisSetName;
  /**
   * Use spherical harmonic basis functions (5 d-orbitals) instead of Cartesian (6 d-orbitals).
   * For s and p orbitals, both choices give the same result.
   * @default false (Cartesian basis for backward compatibility)
   */
  useSpherical?: boolean;
}

/**
 * Cancel a running computation.
 *
 * The worker will set an abort flag that handlers can check.
 */
export interface CancelRequest extends BaseRequest {
  type: 'cancel';
  /** The request ID of the computation to cancel */
  targetRequestId: RequestId;
}

/**
 * Union of all worker request types.
 *
 * Use the `type` discriminator to narrow the type:
 *
 * @example
 * ```typescript
 * function handleRequest(request: WorkerRequest) {
 *   switch (request.type) {
 *     case 'ping':
 *       // request is narrowed to PingRequest
 *       break;
 *     case 'boys_eval':
 *       // request is narrowed to BoysEvalRequest
 *       console.log(request.m, request.T);
 *       break;
 *     // ... handle other cases
 *   }
 * }
 * ```
 */
export type WorkerRequest =
  | PingRequest
  | BoysEvalRequest
  | BoysSweepRequest
  | RysComputeRequest
  | RysErrorCurveRequest
  | ScfRunRequest
  | IntegralComputeRequest
  | CancelRequest;

/**
 * All valid request type strings.
 */
export type WorkerRequestType = WorkerRequest['type'];

// ============================================================================
// Worker Responses
// ============================================================================

/**
 * Base fields present in all responses.
 */
interface BaseResponse {
  /** Request ID this response correlates to */
  requestId: RequestId;
}

/**
 * Successful pong response to ping request.
 */
export interface PongResponse extends BaseResponse {
  type: 'pong';
  /** Version string from the WASM module */
  wasmVersion: string;
  /** Whether threading is available in this WASM build */
  threadsAvailable: boolean;
  /** Number of threads initialized (0 if threading not available) */
  numThreads: number;
}

/**
 * Successful result response.
 *
 * The `data` field type depends on the request type.
 * Use type inference from the request to narrow the data type.
 */
export interface ResultResponse extends BaseResponse {
  type: 'result';
  /** Result data (type varies by request) */
  data: unknown;
}

/**
 * Error response.
 *
 * All errors include a machine-readable code and human-readable message.
 */
export interface ErrorResponse extends BaseResponse {
  type: 'error';
  /** Machine-readable error code */
  code: WorkerErrorCode;
  /** Human-readable error message */
  message: string;
}

/**
 * Progress update response for streaming operations.
 *
 * Used primarily for SCF iteration updates.
 */
export interface ProgressResponse extends BaseResponse {
  type: 'progress';
  /** Progress details */
  progress: WorkerProgress;
}

/**
 * Union of all worker response types.
 *
 * Use the `type` discriminator to narrow the type:
 *
 * @example
 * ```typescript
 * function handleResponse(response: WorkerResponse) {
 *   switch (response.type) {
 *     case 'pong':
 *       console.log('WASM version:', response.wasmVersion);
 *       break;
 *     case 'result':
 *       console.log('Data:', response.data);
 *       break;
 *     case 'error':
 *       console.error(`[${response.code}] ${response.message}`);
 *       break;
 *     case 'progress':
 *       updateProgress(response.progress);
 *       break;
 *   }
 * }
 * ```
 */
export type WorkerResponse =
  | PongResponse
  | ResultResponse
  | ErrorResponse
  | ProgressResponse;

/**
 * All valid response type strings.
 */
export type WorkerResponseType = WorkerResponse['type'];

// ============================================================================
// Error Codes
// ============================================================================

/**
 * Worker error codes.
 *
 * These codes help the UI determine appropriate error handling and
 * user-facing messages.
 */
export type WorkerErrorCode =
  /** Worker or WASM not yet initialized */
  | 'WORKER_NOT_READY'
  /** WASM module failed to initialize */
  | 'WASM_INIT_FAILED'
  /** Unrecognized request type */
  | 'UNKNOWN_REQUEST_TYPE'
  /** Handler exists but not yet implemented */
  | 'NOT_IMPLEMENTED'
  /** Handler threw an exception during execution */
  | 'HANDLER_ERROR'
  /** Invalid parameters provided */
  | 'INVALID_PARAMS'
  /** Computation was cancelled by user */
  | 'COMPUTATION_CANCELLED';

// ============================================================================
// Progress Types
// ============================================================================

/**
 * Module identifier for progress updates.
 */
export type ProgressModule = 'boys' | 'rys' | 'scf' | 'integral';

/**
 * Base progress fields common to all modules.
 */
interface BaseProgress {
  /** Current step number */
  current: number;
  /** Total expected steps (0 if unknown) */
  total: number;
  /** Human-readable progress message */
  message: string;
}

/**
 * SCF-specific iteration progress.
 *
 * Posted after each SCF iteration with convergence information.
 */
export interface ScfIterationProgress extends BaseProgress {
  module: 'scf';
  /** Current iteration number (1-indexed) */
  iteration: number;
  /** Current total electronic energy (Hartree) */
  energy: number;
  /** Energy change from previous iteration */
  delta: number;
  /** DIIS error metric (if DIIS enabled) */
  diisError?: number;
  /** Whether convergence criteria met */
  converged: boolean;
}

/**
 * Generic progress for Boys function operations.
 *
 * Used for sweep operations with many points.
 */
export interface BoysProgress extends BaseProgress {
  module: 'boys';
}

/**
 * Generic progress for Rys quadrature operations.
 */
export interface RysProgress extends BaseProgress {
  module: 'rys';
}

/**
 * Integral computation phase identifier.
 *
 * Phases are processed in order, with ERI being the slowest.
 * Weight distribution: overlap (5%), kinetic (5%), nuclear (10%),
 * eri (75%), assembly (5%).
 */
export type IntegralPhase = 'overlap' | 'kinetic' | 'nuclear' | 'eri' | 'assembly';

/**
 * Integral computation progress update.
 *
 * Emitted during each phase of integral computation to provide
 * feedback on long-running computations.
 *
 * @see US-029 Integral WASM Integration
 */
export interface IntegralProgress extends BaseProgress {
  module: 'integral';
  /** Current computation phase */
  phase: IntegralPhase;
  /** Overall completion percentage (0-100) */
  overallPercent: number;
}

/**
 * Union of all progress types.
 *
 * Use the `module` discriminator to narrow the type.
 */
export type WorkerProgress =
  | ScfIterationProgress
  | BoysProgress
  | RysProgress
  | IntegralProgress;

// ============================================================================
// Result Types (for type inference)
// ============================================================================

/**
 * Result from ping request.
 */
export interface PingResult {
  wasmVersion: string;
  /** Whether threading is available in this WASM build */
  threadsAvailable: boolean;
  /** Number of threads initialized (0 if threading not available) */
  numThreads: number;
}

/**
 * Computational method used for Boys function evaluation.
 * Matches the enum from qc-core.
 */
export type BoysMethod = 'zero' | 'series' | 'recurrence';

/**
 * Result from boys_eval request.
 * Matches the BoysResult struct from qc-core/qc-wasm.
 */
export interface BoysEvalResult {
  /** Computed value F_m(T) */
  value: number;
  /** Which computational method was used */
  method: BoysMethod;
  /** The order m of the Boys function */
  m: number;
  /** The argument T */
  t: number;
  /**
   * Number of terms/steps used in computation.
   * - Zero method: 0 (direct formula)
   * - Series: number of iteration terms until convergence
   * - Recurrence: m + 1 (upward recurrence steps)
   */
  termsCount: number;
  /**
   * Estimated relative error bound (if available).
   * - Zero method: machine epsilon
   * - Series: estimated from convergence tolerance
   * - Recurrence: null (difficult to estimate)
   */
  estimatedError: number | null;
}

/**
 * Result from boys_sweep request.
 */
export interface BoysSweepResult {
  /** Array of individual point results */
  results: BoysEvalResult[];
  /** The order m used for all evaluations */
  m: number;
  /** Range of T values [min, max] */
  T_range: [number, number];
  /** Number of points computed */
  points: number;
}

/**
 * Computational method used for Rys quadrature.
 * Matches the RysMethod enum from qc-core.
 */
export type RysMethod = 'special' | 'standard';

/**
 * Result from rys_compute request.
 * Matches the RysResult struct from qc-core.
 */
export interface RysComputeResult {
  /** Quadrature roots in the interval (0, 1) */
  roots: number[];
  /** Corresponding weights (all strictly positive) */
  weights: number[];
  /** Number of roots/weights (== roots.length == weights.length) */
  nroots: number;
  /** The argument T used for computation */
  t: number;
  /** Computational method used ("special" for T=0 or n=1, "standard" otherwise) */
  method: RysMethod;
}

/**
 * A single point on the error curve showing max reconstruction error for a given order.
 * Matches the ErrorCurvePoint struct from qc-core.
 */
export interface ErrorCurvePoint {
  /** Quadrature order (number of roots/weights) */
  n: number;
  /** Maximum absolute reconstruction error across all moments 0..2n-1 */
  maxError: number;
}

/**
 * Result from rys_error_curve request.
 * Matches the ErrorCurveResult struct from qc-core.
 */
export interface RysErrorCurveResult {
  /** The argument T used for computation */
  t: number;
  /** Maximum order computed (n_max) */
  nMax: number;
  /** Error curve points for n = 1, 2, ..., n_max */
  points: ErrorCurvePoint[];
}

/**
 * SCF matrices for internals mode visualization.
 *
 * Contains the key matrices from an SCF calculation stored as
 * row-major flat arrays. Use `nbf` to reshape into square matrices.
 *
 * @see US-018 SCF Internals Mode
 */
export interface ScfMatrices {
  /** Number of basis functions (for reshaping into nbf x nbf matrix) */
  nbf: number;
  /** Overlap matrix S (row-major, nbf x nbf) */
  sMatrix: number[];
  /** Core Hamiltonian H = T + V (row-major, nbf x nbf) */
  hCore: number[];
  /** Final Fock matrix F (row-major, nbf x nbf) */
  fockMatrix: number[];
  /** Final density matrix D (row-major, nbf x nbf) */
  densityMatrix: number[];
}

/**
 * Orbital energies from SCF calculation.
 *
 * Contains MO (molecular orbital) energies with occupancy information
 * for energy level diagram visualization.
 *
 * @see US-018 SCF Internals Mode
 */
export interface OrbitalEnergies {
  /** Orbital energies in Hartree, sorted ascending */
  energies: number[];
  /** Number of occupied orbitals (n_occ = nelec / 2 for RHF) */
  nOccupied: number;
}

/**
 * Result from scf_run request.
 */
export interface ScfRunResult {
  /** Final total electronic energy (Hartree) */
  energy: number;
  /** Whether SCF converged */
  converged: boolean;
  /** Number of iterations performed */
  iterations: number;
  /** Whether computation was aborted */
  aborted: boolean;
  /** Iteration history for plotting */
  history: ScfIterationHistory[];
  /** Optional matrices (when includeMatrices option is true) */
  matrices?: ScfMatrices;
  /** Optional orbital energies (when includeMatrices option is true) */
  orbitalEnergies?: OrbitalEnergies;
}

/**
 * Single iteration in SCF history.
 */
export interface ScfIterationHistory {
  iteration: number;
  energy: number;
  delta: number;
  diisError?: number;
}

/**
 * Result from cancel request.
 */
export interface CancelResult {
  /** Whether cancellation was acknowledged */
  cancelled: boolean;
  /** The request ID that was cancelled */
  targetRequestId: RequestId;
}

/**
 * Atom output in integral result geometry.
 *
 * Echoes back the input atoms with added atomic number information.
 */
export interface AtomOutput {
  /** Element symbol */
  symbol: string;
  /** Cartesian coordinates [x, y, z] in Bohr */
  xyz: [number, number, number];
  /** Atomic number (Z) */
  atomicNumber: number;
}

/**
 * Geometry output echoed back in integral result.
 */
export interface GeometryOutput {
  /** List of atoms with atomic numbers */
  atoms: AtomOutput[];
  /** Coordinate units (always "bohr" in output) */
  units: CoordinateUnits;
}

/**
 * Basis type used for computation.
 * - "cartesian": 6 Cartesian d-orbitals (xx, xy, xz, yy, yz, zz)
 * - "spherical": 5 spherical harmonic d-orbitals (d-2, d-1, d0, d+1, d+2)
 */
export type BasisType = 'cartesian' | 'spherical';

/**
 * Computation metadata from integral calculation.
 */
export interface IntegralMetadata {
  /** WASM module version */
  wasmVersion: string;
  /** Computation time in milliseconds */
  computeTimeMs: number;
  /** Number of shell pairs processed */
  shellPairs: number;
  /** Number of shell quartets processed (for ERI) */
  shellQuartets: number;
  /** Number of significant ERIs (after screening) */
  significantEris: number;
  /** Basis type used: "cartesian" (6 d-orbitals) or "spherical" (5 d-orbitals) */
  basisType: BasisType;
}

/**
 * Result from integral_compute request.
 *
 * This structure is compatible with PresetSystemJson, allowing direct
 * use with the existing SCF workflow without modification.
 *
 * The result contains all pre-computed integrals needed for an SCF
 * calculation: overlap (S), core Hamiltonian (H = T + V), and
 * two-electron integrals (ERI).
 *
 * @see US-029 Integral WASM Integration
 */
export interface IntegralComputeResult {
  /** Format version (1) */
  formatVersion: number;
  /** Generated system ID (e.g., "custom_abc123") */
  systemId: string;
  /** Human-readable label */
  label: string;
  /** Description of the system */
  description: string;
  /** Input geometry (echoed back, always in Bohr) */
  geometry: GeometryOutput;
  /** Basis set used */
  basisId: string;
  /** Number of basis functions */
  nbf: number;
  /** Number of electrons */
  nelec: number;
  /** Nuclear repulsion energy (Hartree) */
  eNuc: number;
  /** Overlap matrix S (row-major, nbf x nbf) */
  sMatrix: number[];
  /** Core Hamiltonian H = T + V (row-major, nbf x nbf) */
  hCore: number[];
  /** Compressed two-electron integrals (8-fold symmetry) */
  eriCompressed: number[];
  /** Indexing scheme description for ERI */
  eriIndexing: string;
  /** Computation metadata */
  metadata: IntegralMetadata;
}

// ============================================================================
// Type Helpers
// ============================================================================

/**
 * Extract the result type for a given request type.
 *
 * This enables type-safe responses from the worker.
 *
 * @example
 * ```typescript
 * // ResultFor<'ping'> = PingResult
 * // ResultFor<'boys_eval'> = BoysEvalResult
 * ```
 */
export type ResultFor<T extends WorkerRequestType> = T extends 'ping'
  ? PingResult
  : T extends 'boys_eval'
    ? BoysEvalResult
    : T extends 'boys_sweep'
      ? BoysSweepResult
      : T extends 'rys_compute'
        ? RysComputeResult
        : T extends 'rys_error_curve'
          ? RysErrorCurveResult
          : T extends 'scf_run'
            ? ScfRunResult
            : T extends 'integral_compute'
              ? IntegralComputeResult
              : T extends 'cancel'
                ? CancelResult
                : never;

/**
 * Type guard to check if a value is a valid WorkerRequest.
 */
export function isWorkerRequest(value: unknown): value is WorkerRequest {
  if (typeof value !== 'object' || value === null) {
    return false;
  }
  const obj = value as Record<string, unknown>;
  if (typeof obj.type !== 'string' || typeof obj.requestId !== 'string') {
    return false;
  }
  const validTypes: WorkerRequestType[] = [
    'ping',
    'boys_eval',
    'boys_sweep',
    'rys_compute',
    'rys_error_curve',
    'scf_run',
    'integral_compute',
    'cancel',
  ];
  return validTypes.includes(obj.type as WorkerRequestType);
}

/**
 * Type guard to check if a value is a valid WorkerResponse.
 */
export function isWorkerResponse(value: unknown): value is WorkerResponse {
  if (typeof value !== 'object' || value === null) {
    return false;
  }
  const obj = value as Record<string, unknown>;
  if (typeof obj.type !== 'string' || typeof obj.requestId !== 'string') {
    return false;
  }
  const validTypes: WorkerResponseType[] = ['pong', 'result', 'error', 'progress'];
  return validTypes.includes(obj.type as WorkerResponseType);
}

/**
 * Assert that a switch statement is exhaustive.
 *
 * Use this in the default case of a switch on a discriminated union
 * to get a compile-time error if a case is missed.
 *
 * @example
 * ```typescript
 * switch (request.type) {
 *   case 'ping': // ...
 *   case 'boys_eval': // ...
 *   // If you forget a case, TypeScript will error
 *   default:
 *     assertNever(request);
 * }
 * ```
 */
export function assertNever(value: never): never {
  throw new Error(`Unexpected value: ${JSON.stringify(value)}`);
}
