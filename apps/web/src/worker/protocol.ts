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
export type BasisSetName = 'sto-3g' | '3-21g' | '6-31g' | '6-31g*' | '6-31+g*' | 'cc-pvdz';

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
  /**
   * Level shift for virtual orbitals (Hartree).
   *
   * Shifts virtual orbital energies up by this amount to widen the
   * HOMO-LUMO gap and stabilize SCF convergence.
   *
   * @default 0.0 (no level shift)
   */
  levelShift?: number;
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
 * @see Module C - Boys Function Lab
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
 * Boys function evaluation for all orders 0..mMax at a single T.
 *
 * Used for multi-order comparison: shows how F_m(T) varies with m
 * at a fixed T value. This is pedagogically valuable because students
 * can observe the monotonic decrease of F_m(T) with increasing m.
 *
 * @see Module C - Boys Function Lab
 */
export interface BoysEvalAllRequest extends BaseRequest {
  type: 'boys_eval_all';
  /** Maximum order m to compute (evaluates F_0(T) through F_mMax(T)) */
  mMax: number;
  /** Argument T (must be >= 0) */
  T: number;
}

/**
 * Result from boys_eval_all request.
 */
export interface BoysEvalAllResult {
  /** Array of results for m=0..mMax */
  results: BoysEvalResult[];
  /** Maximum order computed */
  mMax: number;
  /** The T value used */
  T: number;
}

/**
 * Rys quadrature roots and weights computation.
 *
 * @see Module D - Rys Quadrature Lab
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
 * KS-DFT SCF computation request.
 *
 * Runs a Kohn-Sham DFT calculation using the specified functional
 * and numerical integration grid.
 *
 * @see US-068 KS-SCF Loop + V_xc Matrix
 */
export interface KsScfRequest extends BaseRequest {
  type: 'ks_scf';
  /** Atoms as [Z, x, y, z] arrays (coordinates in bohr) */
  atoms: [number, number, number, number][];
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** DFT method: "lda" */
  method: string;
  /** Convergence profile: "loose" | "medium" | "tight" */
  convergenceProfile?: string;
  /** Maximum SCF iterations */
  maxIterations?: number;
  /** Enable DIIS acceleration */
  useDiis?: boolean;
  /** Grid quality: "standard" | "fine" */
  gridQuality?: string;
  /** Use spherical harmonic basis functions (5 d-orbitals vs 6 Cartesian) */
  useSpherical?: boolean;
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
 * PES (Potential Energy Surface) bond-length scan request.
 *
 * Scans the potential energy surface of a diatomic molecule by varying
 * the internuclear distance. Emits progress after each scan point.
 *
 * @see US-039 PES WASM Export & Worker
 */
export interface PesScanRequest extends BaseRequest {
  type: 'pes_scan';
  /** Atomic number of atom A (fixed at origin), 1-10 */
  atomAZ: number;
  /** Atomic number of atom B (translated along z-axis), 1-10 */
  atomBZ: number;
  /** Minimum bond distance in bohr */
  rMin: number;
  /** Maximum bond distance in bohr */
  rMax: number;
  /** Number of scan points (evenly spaced) */
  nPoints: number;
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** SCF computation options */
  options: ScfOptions;
  /** Whether to use convergence seeding from previous point (default: true) */
  useSeeding?: boolean;
}

/**
 * MO grid evaluation request.
 *
 * Evaluates a molecular orbital on a 3D grid for isosurface extraction.
 *
 * @see US-042 MO Grid Evaluation
 */
export interface MoGridRequest extends BaseRequest {
  type: 'mo_grid';
  /** MO coefficient vector (one per basis function) */
  moCoefficients: number[];
  /** Atom specifications: [[Z, x, y, z], ...] in Bohr */
  atoms: [number, number, number, number][];
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** Grid origin [x, y, z] in Bohr */
  gridOrigin: [number, number, number];
  /** Grid spacing in Bohr (uniform) */
  gridSpacing: number;
  /** Grid dimensions [nx, ny, nz] */
  gridDims: [number, number, number];
  /**
   * Whether the MO coefficients are in the spherical harmonic basis.
   * Must match the `useSpherical` option used in the SCF/integral computation.
   * When true, WASM transforms coefficients from spherical (5 d-functions) to
   * Cartesian (6 d-functions) before grid evaluation.
   * @default false
   */
  useSpherical?: boolean;
}

/**
 * Marching cubes isosurface extraction request.
 *
 * Extracts an isosurface from a 3D scalar field (e.g., from MO grid evaluation).
 * Returns triangle mesh data directly consumable by Three.js BufferGeometry.
 *
 * @see US-043 Marching Cubes Isosurface
 */
export interface MarchingCubesRequest extends BaseRequest {
  type: 'marching_cubes';
  /** Flat scalar field data (C-order: x-slowest, z-fastest) */
  gridData: number[];
  /** Grid dimensions [nx, ny, nz] */
  gridDims: [number, number, number];
  /** Grid origin [x, y, z] in Bohr */
  gridOrigin: [number, number, number];
  /** Grid spacing in Bohr (uniform) */
  gridSpacing: number;
  /** Isovalue threshold */
  isovalue: number;
}

/**
 * Dual marching cubes request for orbital visualization.
 *
 * Extracts isosurfaces at both +isovalue and -isovalue for positive
 * and negative orbital lobes.
 *
 * @see US-043 Marching Cubes Isosurface
 */
export interface DualMarchingCubesRequest extends BaseRequest {
  type: 'dual_marching_cubes';
  /** Flat scalar field data (C-order: x-slowest, z-fastest) */
  gridData: number[];
  /** Grid dimensions [nx, ny, nz] */
  gridDims: [number, number, number];
  /** Grid origin [x, y, z] in Bohr */
  gridOrigin: [number, number, number];
  /** Grid spacing in Bohr (uniform) */
  gridSpacing: number;
  /** Isovalue threshold (positive; negative lobe uses -isovalue) */
  isovalue: number;
}

/**
 * Basis set shell information request.
 *
 * Queries the WASM module for the shell structure of a given element
 * in a specific basis set. Returns exponents and contraction coefficients
 * sourced directly from the built-in basis data in qc-core.
 *
 * @see Module D - Basis Explorer (US-049 AC6)
 */
export interface BasisInfoRequest extends BaseRequest {
  type: 'basis_info';
  /** Atomic number of the element (1-18) */
  atomicNumber: number;
  /** Basis set name (e.g., "sto-3g", "6-31g*") */
  basisName: string;
}

/**
 * Radial profile evaluation request.
 *
 * Evaluates the radial part of a contracted basis shell for
 * visualization in Module D's radial profile plot.
 *
 * @see US-050 Radial Profile Visualization
 */
export interface RadialProfileRequest extends BaseRequest {
  type: 'radial_profile';
  /** Atomic number (1-18) */
  atomicNumber: number;
  /** Basis set name (e.g., "sto-3g", "6-31g*") */
  basisName: string;
  /** Shell index (0-indexed, matching order from basis_info) */
  shellIndex: number;
  /** Number of evaluation points (default: 200) */
  nPoints?: number;
  /** Optional maximum r in Bohr (auto-determined if absent) */
  rMax?: number;
}

/**
 * Overlap vs. distance computation request.
 *
 * Computes the overlap integral S_ab(R) for two selected basis shells
 * at n_points evenly spaced distances from r_min to r_max.
 *
 * @see US-054 Overlap vs. Distance Plot (Module D)
 */
export interface OverlapDistanceRequest extends BaseRequest {
  type: 'overlap_distance';
  /** Atomic number for atom A (1-18) */
  elementA: number;
  /** Basis set name for atom A */
  basisA: BasisSetName;
  /** Shell index within atom A's basis (0-indexed) */
  shellIndexA: number;
  /** Atomic number for atom B (1-18) */
  elementB: number;
  /** Basis set name for atom B */
  basisB: BasisSetName;
  /** Shell index within atom B's basis (0-indexed) */
  shellIndexB: number;
  /** Minimum distance in bohr */
  rMin: number;
  /** Maximum distance in bohr */
  rMax: number;
  /** Number of distance points */
  nPoints: number;
}

/**
 * Integral matrices computation request.
 *
 * Computes one-electron integral matrices (S, T, V, H^core) for a molecule.
 * Unlike integral_compute, this returns T and V separately and does NOT
 * compute two-electron integrals (ERIs), making it faster for heatmap display.
 *
 * @see US-055 Integral Matrix Heatmap
 */
export interface IntegralMatricesRequest extends BaseRequest {
  type: 'integral_matrices';
  /** Geometry specification with atoms and coordinate units */
  geometry: GeometryInput;
  /** Basis set name (e.g., "sto-3g") */
  basisName: BasisSetName;
  /** Whether to use spherical harmonics (default: false) */
  useSpherical?: boolean;
}

/**
 * Request to compute a single integral with primitive-pair decomposition.
 *
 * Returns the contracted integral value and all primitive-pair contributions
 * sorted by magnitude, enabling inspection of which Gaussian pairs dominate.
 *
 * @see US-056 Primitive Breakdown Panel
 * @see FR-INT-03
 */
export interface IntegralBreakdownRequest extends BaseRequest {
  type: 'integral_breakdown';
  /** Geometry specification with atoms and coordinate units */
  geometry: GeometryInput;
  /** Basis set name */
  basisName: BasisSetName;
  /** Integral type */
  integralType: 'S' | 'T' | 'V' | 'Hcore';
  /** Basis function indices [row, col] (0-based) */
  indices: [number, number];
}

/**
 * Fock matrix decomposition request.
 *
 * Decomposes F = H^core + G(P) into separate J (Coulomb) and K (Exchange)
 * contributions for educational inspection. Requires a converged density
 * matrix from a prior SCF run.
 *
 * @see US-058 Fock Build Tracing
 * @see FR-INT-05
 */
export interface FockDecompositionRequest extends BaseRequest {
  type: 'fock_decomposition';
  /** Molecule geometry */
  geometry: GeometryInput;
  /** Basis set name */
  basisSet: BasisSetName;
  /** Density matrix P (flat, nbf x nbf, includes factor of 2 for RHF) */
  densityMatrix: number[];
}

/**
 * ERI detail request.
 *
 * Decomposes a contracted ERI (ij|kl) into its primitive-quartet
 * contributions with method metadata (Boys vs Rys).
 *
 * @see US-059 ERI Browser
 * @see FR-INT-04
 */
export interface EriDetailRequest extends BaseRequest {
  type: 'eri_detail';
  /** Molecule geometry */
  geometry: GeometryInput;
  /** Basis set name */
  basisName: BasisSetName;
  /** Basis function indices [i, j, k, l] (0-based) */
  indices: [number, number, number, number];
}

/**
 * Density grid evaluation request.
 *
 * Evaluates the electron density rho(r) on a 3D grid for isosurface extraction.
 * The density is computed as rho(r) = sum_{mu,nu} D_{mu,nu} * chi_mu(r) * chi_nu(r).
 *
 * @see US-061 Density Isosurface Visualization
 */
export interface DensityGridRequest extends BaseRequest {
  type: 'density_grid';
  /** Flattened density matrix (row-major, n_basis x n_basis) */
  densityMatrix: number[];
  /** Atom specifications: [[Z, x, y, z], ...] in Bohr */
  atoms: [number, number, number, number][];
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** Grid origin [x, y, z] in Bohr */
  gridOrigin: [number, number, number];
  /** Grid spacing in Bohr (uniform) */
  gridSpacing: number;
  /** Grid dimensions [nx, ny, nz] */
  gridDims: [number, number, number];
  /** Number of electrons (for integrated density validation) */
  nElectrons: number;
  /**
   * Whether the density matrix is in the spherical harmonic basis.
   * Must match the `useSpherical` option used in the SCF/integral computation.
   * When true, WASM transforms the density matrix from spherical to
   * Cartesian before grid evaluation.
   * @default false
   */
  useSpherical?: boolean;
}

/**
 * Difference density grid evaluation request.
 *
 * Computes Delta-rho = rho_molecule - rho_promolecule on a 3D grid.
 * The molecular density grid is passed in (cached from a prior density_grid
 * request), and the promolecule is evaluated internally using embedded
 * atomic density profiles.
 *
 * The result contains the difference density grid values, integrated
 * Delta-rho (should be approximately zero for density conservation),
 * and max accumulation/depletion values.
 *
 * @see US-063 Difference Density
 */
export interface DifferenceDensityRequest extends BaseRequest {
  type: 'difference_density';
  /** Cached total molecular density grid from a prior density_grid request */
  totalDensity: number[];
  /** Atom specifications: [[Z, x, y, z], ...] in Bohr */
  atoms: [number, number, number, number][];
  /** Grid origin [x, y, z] in Bohr */
  gridOrigin: [number, number, number];
  /** Grid spacing in Bohr (uniform) */
  gridSpacing: number;
  /** Grid dimensions [nx, ny, nz] */
  gridDims: [number, number, number];
}

/**
 * Result from difference density grid evaluation.
 *
 * @see US-063 Difference Density
 */
export interface DifferenceDensityResult {
  /** Flat array of Delta-rho values (C-order: x-slowest, z-fastest) */
  values: number[];
  /** Grid dimensions [nx, ny, nz] */
  gridDims: [number, number, number];
  /** Grid origin [x, y, z] in Bohr */
  gridOrigin: [number, number, number];
  /** Grid spacing in Bohr */
  gridSpacing: number;
  /** Integrated Delta-rho (sum * dV, should be approximately 0 for density conservation) */
  integratedDeltaRho: number;
  /** Maximum positive value (accumulation peak) */
  maxAccumulation: number;
  /** Maximum negative value (depletion peak, stored as negative number) */
  maxDepletion: number;
  /** Computation time in milliseconds */
  computeTimeMs: number;
}

/**
 * Geometry optimization request.
 *
 * Runs L-BFGS geometry optimization using the specified method and basis set.
 * Emits progress after each optimization step with energy and gradient.
 *
 * @see US-075 Optimization UI + Trajectory Animation
 */
export interface OptimizeGeometryRequest extends BaseRequest {
  type: 'optimize_geometry';
  /** Atoms as [Z, x, y, z] arrays (coordinates in bohr) */
  atoms: [number, number, number, number][];
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** Electronic structure method: "rhf", "lda", or "b3lyp" */
  method: string;
  /** Maximum optimization steps (default: 50) */
  maxSteps?: number;
  /** Maximum gradient convergence threshold in Ha/bohr (default: 4.5e-4) */
  gradThreshold?: number;
  /** Energy convergence threshold in Ha (default: 1e-6) */
  energyThreshold?: number;
}

/**
 * Internal coordinate PES scan request.
 *
 * Scans the potential energy surface along a bond, angle, or dihedral
 * coordinate. Supports both rigid and relaxed scan modes.
 *
 * @see US-081 PES Scan WASM Export + Worker Handler
 */
export interface PesScanInternalRequest extends BaseRequest {
  type: 'pes_scan_internal';
  /** Atoms as [Z, x, y, z] arrays (coordinates in bohr) */
  atoms: [number, number, number, number][];
  /** Basis set name (e.g., "sto-3g") */
  basisName: string;
  /** Electronic structure method: "rhf", "lda", "b3lyp", "b3lyp-d3bj" */
  method: string;
  /** Type of coordinate to scan */
  coordinateType: 'bond' | 'angle' | 'dihedral';
  /** Atom indices defining the coordinate (2, 3, or 4 indices) */
  atomIndices: number[];
  /** Scan mode: rigid (frozen geometry) or relaxed (constrained optimization) */
  scanMode: 'rigid' | 'relaxed';
  /** Minimum coordinate value (bohr for bonds, radians for angles) */
  valueMin: number;
  /** Maximum coordinate value */
  valueMax: number;
  /** Number of evenly spaced scan points (>= 2) */
  nPoints: number;
  /** Whether to seed density from previous point (default: true) */
  useSeeding?: boolean;
  /** Whether to use spherical harmonics (default: true) */
  useSpherical?: boolean;
  /** Convergence profile (default: "tight") */
  convergenceProfile?: string;
  /** Max optimization steps per scan point for relaxed scans (default: 50) */
  optMaxSteps?: number;
  /** Gradient convergence threshold for relaxed scans in Ha/bohr (default: 4.5e-4) */
  optGradThreshold?: number;
}

/**
 * Population analysis request.
 *
 * Computes Mulliken and Lowdin atomic charges from the density and
 * overlap matrices. Requires SCF to have converged with includeMatrices.
 *
 * @see US-076 Mulliken/Lowdin Population Analysis
 */
export interface PopulationAnalysisRequest extends BaseRequest {
  type: 'population_analysis';
  /** Flattened density matrix (row-major, nbf x nbf) */
  densityMatrix: number[];
  /** Flattened overlap matrix (row-major, nbf x nbf) */
  overlapMatrix: number[];
  /** Number of basis functions */
  nbf: number;
  /** Atom specifications with basis function counts */
  atoms: { atomicNumber: number; nBasis: number }[];
}

// ============================================================================
// Frequency Analysis (US-101)
// ============================================================================

/**
 * Broadening kernel for simulated IR / Raman spectra.
 *
 * Matches `qc_core::spectra::BroadeningKind` (lowercased at the WASM boundary).
 *
 * @see US-101 Frequency WASM Export + Worker Handler
 */
export type BroadeningKind = 'lorentzian' | 'gaussian';

/**
 * Rotor classification from the harmonic analysis pipeline.
 *
 * Matches `qc_core::thermo::RotorType` serialized as snake_case. The enum
 * has exactly 5 variants — there is no oblate/prolate split; both are
 * reported as `'symmetric_top'`.
 *
 * @see crates/qc-core/src/thermo.rs RotorType
 * @see US-101 Frequency WASM Export + Worker Handler
 */
export type RotorType =
  | 'atom'
  | 'linear'
  | 'spherical_top'
  | 'symmetric_top'
  | 'asymmetric_top';

/**
 * Frequency analysis request.
 *
 * Runs the full analytical vibrational spectroscopy pipeline:
 * Hessian → normal modes → IR intensities → Raman activities →
 * RRHO thermochemistry → broadened IR/Raman spectra.
 *
 * // MUST match FrequencyWasmInput in crates/qc-wasm/src/lib.rs
 *
 * @see US-101 Frequency WASM Export + Worker Handler
 */
export interface FrequencyRequest extends BaseRequest {
  type: 'frequency';
  /** Atoms as [Z, x, y, z] arrays (atomic number as number, coordinates in bohr) */
  atoms: [number, number, number, number][];
  /** Basis set name (e.g., "sto-3g", "6-31g*", "cc-pvdz") */
  basisName: string;
  /**
   * Electronic structure method (case-insensitive):
   * `"rhf"`, `"hf"` (alias for rhf), `"lda"`, `"b3lyp"`, or `"b3lyp-d3bj"`.
   */
  method: string;
  /** Temperature in Kelvin (default 298.15). Must be > 0 if provided. */
  temperatureK?: number;
  /** Pressure in Pascals (default 101325 = 1 atm). Must be > 0 if provided. */
  pressurePa?: number;
  /**
   * Rotational symmetry number override (σ).
   * If omitted, thermochemistry defaults σ = 1.
   * Common values: 1 (C₁, Cs), 2 (C₂ᵥ H₂O), 3 (C₃ᵥ NH₃), 12 (Tᵈ CH₄).
   */
  symmetryNumberOverride?: number;
  /** Spin multiplicity 2S+1 (default 1 for singlet). Must be >= 1. */
  multiplicity?: number;
  /** Spectrum broadening kind (default "lorentzian") */
  broadeningKind?: BroadeningKind;
  /** FWHM in cm⁻¹ for spectrum broadening (default 8.0). Must be > 0. */
  fwhmCm1?: number;
  /** SCF convergence profile (default "tight" for accurate frequencies) */
  convergenceProfile?: ConvergenceProfile;
  /** Maximum SCF iterations (default 100) */
  maxIterations?: number;
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
  | BoysEvalAllRequest
  | RysComputeRequest
  | RysErrorCurveRequest
  | ScfRunRequest
  | KsScfRequest
  | IntegralComputeRequest
  | PesScanRequest
  | MoGridRequest
  | MarchingCubesRequest
  | DualMarchingCubesRequest
  | BasisInfoRequest
  | RadialProfileRequest
  | OverlapDistanceRequest
  | IntegralMatricesRequest
  | IntegralBreakdownRequest
  | FockDecompositionRequest
  | EriDetailRequest
  | DensityGridRequest
  | DifferenceDensityRequest
  | OptimizeGeometryRequest
  | PesScanInternalRequest
  | PopulationAnalysisRequest
  | FrequencyRequest
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
export type ProgressModule =
  | 'boys'
  | 'rys'
  | 'scf'
  | 'integral'
  | 'pes'
  | 'pes_internal'
  | 'optimization'
  | 'frequency';

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
 * Progress during integral computation phase of ks_scf.
 *
 * Emitted before SCF iterations begin, during S/Hcore/ERI/grid computation.
 */
export interface ScfIntegralProgress extends BaseProgress {
  module: 'scf_integrals';
  /** Current integral step: "overlap", "hcore", "eri", "grid", "done" */
  step: string;
  /** Overall completion percentage (0-100) */
  percent: number;
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
 * PES scan progress update.
 *
 * Emitted after each scan point completes with the SCF result
 * for that geometry.
 *
 * @see US-039 PES WASM Export & Worker
 */
export interface PesScanProgress extends BaseProgress {
  module: 'pes';
  /** Index of the completed point (0-indexed) */
  pointIndex: number;
  /** Total number of scan points */
  totalPoints: number;
  /** Bond distance in bohr */
  r: number;
  /** SCF energy in Hartree */
  energy: number;
  /** Whether SCF converged at this geometry */
  converged: boolean;
}

/**
 * Geometry optimization progress update.
 *
 * Emitted after each L-BFGS optimization step with energy and gradient info.
 *
 * @see US-075 Optimization UI + Trajectory Animation
 */
export interface OptimizationProgress extends BaseProgress {
  module: 'optimization';
  /** Current step number (0 = initial evaluation) */
  step: number;
  /** Total energy at this step (Ha) */
  energy: number;
  /** Maximum absolute gradient component (Ha/bohr) */
  maxGradient: number;
  /** RMS gradient (Ha/bohr) */
  rmsGradient: number;
}

/**
 * Internal coordinate PES scan progress update.
 *
 * Emitted after each scan point completes with energy and convergence info.
 *
 * @see US-081 PES Scan WASM Export + Worker Handler
 */
export interface PesScanInternalProgress extends BaseProgress {
  module: 'pes_internal';
  /** Index of the completed point (0-indexed) */
  pointIndex: number;
  /** Total number of scan points */
  totalPoints: number;
  /** Value of the scanned coordinate */
  coordinateValue: number;
  /** SCF energy in Hartree */
  energy: number;
  /** Whether SCF converged */
  converged: boolean;
  /** Optimization steps (null for rigid scans) */
  optSteps: number | null;
}

/**
 * Frequency analysis pipeline phase identifier.
 *
 * Phases are emitted in order by the WASM `compute_frequencies` pipeline:
 * 1. `integrals` — SCF + Hessian assembly (includes nuclear CPHF).
 * 2. `nuclear_cphf` — Extract CPHF data + rebuild density matrix for IR.
 * 3. `field_cphf` — Solve field CPHF for polarizability (inside Raman).
 * 4. `assembly` — Harmonic analysis + IR + Raman + thermochemistry.
 * 5. `modes` — Broadened IR/Raman spectrum simulation.
 *
 * @see US-101 Frequency WASM Export + Worker Handler
 */
export type FrequencyPhase =
  | 'integrals'
  | 'nuclear_cphf'
  | 'field_cphf'
  | 'assembly'
  | 'modes';

/**
 * Frequency analysis progress update.
 *
 * Emitted at phase boundaries (and optionally intermediate steps within a
 * phase) by the WASM `compute_frequencies` function. The UI uses `phase` to
 * drive a 5-segment progress bar and `message` for the human-readable label.
 *
 * @see US-101 Frequency WASM Export + Worker Handler
 */
export interface FrequencyProgress extends BaseProgress {
  module: 'frequency';
  /** Current pipeline phase */
  phase: FrequencyPhase;
  /** Fractional completion within the current phase (0.0 – 1.0) */
  percent: number;
  /** Sub-step identifier (e.g., "start", "done", "harmonic", "ir", "raman_sim") */
  step: string;
}

/**
 * Union of all progress types.
 *
 * Use the `module` discriminator to narrow the type.
 */
export type WorkerProgress =
  | ScfIterationProgress
  | ScfIntegralProgress
  | BoysProgress
  | RysProgress
  | IntegralProgress
  | PesScanProgress
  | PesScanInternalProgress
  | OptimizationProgress
  | FrequencyProgress;

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
  /** MO coefficient matrix C (row-major, nbf x nbf) */
  moCoefficients: number[];
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
 * Result from ks_scf request (KS-DFT calculation).
 *
 * @see US-068 KS-SCF Loop + V_xc Matrix
 */
export interface KsScfResult {
  /** Final total energy in Hartree */
  energy: number;
  /** Whether SCF converged */
  converged: boolean;
  /** Number of iterations performed */
  iterations: number;
  /** Exchange-correlation energy (Hartree) */
  energyXc: number;
  /** Coulomb energy (Hartree) */
  energyJ: number;
  /** One-electron energy (Hartree) */
  energy1e: number;
  /** Nuclear repulsion energy (Hartree) */
  energyNuc: number;
  /** Method identifier (e.g., "LDA (Slater + VWN5)") */
  method: string;
  /** Iteration history for plotting */
  trace: ScfIterationHistory[];
  /** Final density matrix (row-major, nbf x nbf) */
  densityMatrix: number[];
  /** MO coefficients (column-major, nbf x nbf) */
  moCoefficients: number[];
  /** Orbital energies (eigenvalues, sorted ascending) */
  orbitalEnergies: number[];
  /** Number of basis functions */
  nBasis: number;
  /** Number of occupied orbitals */
  nOccupied: number;
  /** Overlap matrix S (row-major, nbf x nbf) — for population analysis */
  overlapMatrix: number[];
  /** Core Hamiltonian matrix (row-major, nbf x nbf) */
  hCore: number[];
  /** Final Fock/KS matrix (row-major, nbf x nbf) */
  fockMatrix: number[];
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
 * Shell information returned from WASM basis set query.
 *
 * Contains the angular momentum type, primitive count, exponents,
 * and contraction coefficients for a single contracted shell.
 *
 * @see US-049 AC6 - Data sourced from builtin.rs via WASM
 */
export interface BasisShellInfo {
  /** Angular momentum quantum number: 0=s, 1=p, 2=d */
  angularMomentum: number;
  /** Angular momentum letter: "s", "p", or "d" */
  angularMomentumLabel: string;
  /** Number of primitive Gaussian functions */
  nPrimitives: number;
  /** Exponents of the primitive Gaussians */
  exponents: number[];
  /** Contraction coefficients */
  coefficients: number[];
}

/**
 * Result from basis_info request.
 *
 * Array of shell descriptors for a given element and basis set combination.
 */
export type BasisInfoResult = BasisShellInfo[];

/**
 * Result from radial_profile request.
 *
 * Contains r values, contracted profile, individual primitive profiles,
 * and metadata for plotting and educational tooltips.
 *
 * @see US-050 Radial Profile Visualization
 */
export interface RadialProfileResult {
  /** r values in Bohr */
  rValues: number[];
  /** Contracted profile values R(r) = r^l * sum_k d_k * exp(-alpha_k * r^2) */
  contractedValues: number[];
  /** Individual primitive profiles [primitive_k][r_point] */
  primitiveValues: number[][];
  /** Exponents for each primitive */
  exponents: number[];
  /** Effective coefficients d_k = c_k * N_k */
  effectiveCoefficients: number[];
  /** Raw contraction coefficients c_k */
  rawCoefficients: number[];
  /** Angular momentum l (0=s, 1=p, 2=d) */
  angularMomentum: number;
  /** Angular momentum label ("s", "p", "d") */
  angularMomentumLabel: string;
  /** Number of primitives */
  nPrimitives: number;
  /** Auto-determined r_max in Bohr */
  rMax: number;
}

/**
 * Result from overlap_distance request.
 *
 * Contains the overlap curve data for plotting in Module D.
 *
 * @see US-054 Overlap vs. Distance Plot
 */
export interface OverlapDistanceResult {
  /** Distance values (bohr) */
  rValues: number[];
  /** Overlap integral values S_ab at each distance */
  overlapValues: number[];
  /** Display label for shell A (e.g., "H 1s") */
  shellLabelA: string;
  /** Display label for shell B (e.g., "H 1s") */
  shellLabelB: string;
  /** Basis set name for A */
  basisA: string;
  /** Basis set name for B */
  basisB: string;
}

/**
 * Result from integral_matrices request.
 *
 * Contains all four one-electron integral matrices and basis function labels.
 * Matrices are stored as flat row-major arrays of length nbf*nbf.
 * Access element (i, j) as: matrix[i * nbf + j].
 *
 * @see US-055 Integral Matrix Heatmap
 */
export interface IntegralMatricesResult {
  /** Number of basis functions */
  nbf: number;
  /** Basis function labels: ["O 1s", "O 2s", "O 2px", ...] */
  labels: string[];
  /** Overlap matrix S (row-major, nbf x nbf) */
  sMatrix: number[];
  /** Kinetic energy matrix T (row-major, nbf x nbf) */
  tMatrix: number[];
  /** Nuclear attraction matrix V (row-major, nbf x nbf) */
  vMatrix: number[];
  /** Core Hamiltonian H^core = T + V (row-major, nbf x nbf) */
  hCore: number[];
  /** Nuclear repulsion energy in Hartree */
  nuclearRepulsion: number;
  /** Computation time in milliseconds */
  computeTimeMs: number;
}

/**
 * A single primitive-pair contribution to a contracted integral.
 *
 * Represents one (p, q) term in the double sum over primitive pairs.
 *
 * @see US-056 Primitive Breakdown Panel
 */
export interface PrimitiveContributionResult {
  /** Primitive indices within shells [p, q] */
  primIndices: [number, number];
  /** Primitive exponents [alpha_p, alpha_q] */
  exponents: [number, number];
  /** Raw contraction coefficients [c_p, c_q] */
  coefficients: [number, number];
  /** Coefficients with normalization [c_p * N_p, c_q * N_q] */
  normCoefficients: [number, number];
  /** Bare primitive integral value */
  primitiveValue: number;
  /** Weighted contribution to contracted integral */
  weightedContribution: number;
}

/**
 * Result from integral_breakdown request.
 *
 * Contains the contracted integral value and all primitive-pair contributions
 * sorted by |weightedContribution| descending.
 *
 * @see US-056 Primitive Breakdown Panel
 */
export interface IntegralBreakdownResult {
  /** Contracted integral value */
  contractedValue: number;
  /** Integral type */
  integralType: string;
  /** Basis function indices [i, j] */
  indices: [number, number];
  /** Basis function labels */
  labels: [string, string];
  /** Number of primitives in shell i */
  nPrimI: number;
  /** Number of primitives in shell j */
  nPrimJ: number;
  /** Primitive contributions sorted by |weightedContribution| descending */
  primitiveContributions: PrimitiveContributionResult[];
}

/**
 * Result from fock_decomposition request.
 *
 * Contains the decomposed Fock matrix: F = H^core + J - 0.5*K.
 * All matrices are flat arrays, nbf x nbf, row-major.
 *
 * @see US-058 Fock Build Tracing
 */
export interface FockDecompositionResult {
  /** Core Hamiltonian H^core = T + V (row-major, nbf x nbf) */
  hCore: number[];
  /** Coulomb matrix J_{mn} = sum_{ls} P_{ls} (mn|ls) (row-major, nbf x nbf) */
  jMatrix: number[];
  /** Exchange matrix K_{mn} = sum_{ls} P_{ls} (ml|ns) (row-major, nbf x nbf) */
  kMatrix: number[];
  /** Two-electron contribution G = J - 0.5*K (row-major, nbf x nbf) */
  gMatrix: number[];
  /** Full Fock matrix F = H^core + G (row-major, nbf x nbf) */
  fMatrix: number[];
  /** Density matrix P (echoed, row-major, nbf x nbf; includes factor of 2) */
  density: number[];
  /** Number of basis functions */
  nbf: number;
  /** Basis function labels (e.g., ["H1 1s", "H2 1s"]) */
  labels: string[];
}

/**
 * A single primitive-quartet contribution to a contracted ERI.
 *
 * @see US-059 ERI Browser
 */
export interface EriPrimitiveContribution {
  /** Primitive indices [p, q, r, s] within their shells (0-based) */
  primIndices: [number, number, number, number];
  /** Primitive exponents [alpha_p, alpha_q, alpha_r, alpha_s] */
  exponents: [number, number, number, number];
  /** Raw contraction coefficients (without normalization) */
  coefficients: [number, number, number, number];
  /** Contraction coefficients with normalization */
  normCoefficients: [number, number, number, number];
  /** Bare primitive ERI value */
  primitiveValue: number;
  /** Weighted contribution = product(normCoefficients) * primitiveValue */
  weightedContribution: number;
  /** T parameter for this primitive quartet */
  tParameter: number;
}

/**
 * Method used to compute the ERI.
 *
 * Discriminated union matching Rust's EriMethod enum.
 *
 * @see US-059 ERI Browser
 */
export type EriMethod =
  | { type: 'boysFunction'; tParameter: number }
  | {
      type: 'rysQuadrature';
      nroots: number;
      tParameter: number;
      roots: number[];
      weights: number[];
    };

/**
 * Result from eri_detail request.
 *
 * Contains the contracted ERI value and all primitive-quartet
 * contributions sorted by |weightedContribution| descending.
 *
 * @see US-059 ERI Browser
 */
export interface EriDetailResult {
  /** Contracted ERI value */
  contractedValue: number;
  /** Basis function indices [i, j, k, l] */
  indices: [number, number, number, number];
  /** Basis function labels */
  labels: [string, string, string, string];
  /** Computation method (Boys function or Rys quadrature) */
  method: EriMethod;
  /** Primitive-quartet contributions sorted by |weightedContribution| descending */
  contributions: EriPrimitiveContribution[];
  /** Number of primitives per shell [n_i, n_j, n_k, n_l] */
  nPrimitives: [number, number, number, number];
  /** Total angular momentum L_total */
  totalAngularMomentum: number;
  /** Number of Rys quadrature roots */
  nroots: number;
}

// ============================================================================
// Density Grid Result Types (US-061)
// ============================================================================

/**
 * Result of density grid evaluation.
 *
 * Contains the electron density values on a 3D grid and metadata for
 * validation and visualization.
 *
 * Note: Field names use camelCase (matching Rust serde rename_all = "camelCase").
 *
 * @see US-061 Density Isosurface Visualization
 */
export interface DensityGridResult {
  /** Flat array of density values (C-order: x-slowest, z-fastest) */
  values: number[];
  /** Grid origin [x, y, z] in Bohr */
  gridOrigin: [number, number, number];
  /** Grid spacing in Bohr */
  gridSpacing: number;
  /** Grid dimensions [nx, ny, nz] */
  gridDims: [number, number, number];
  /** Integrated density (numerical integral, should approximate nElectronsExpected) */
  integratedDensity: number;
  /** Expected number of electrons */
  nElectronsExpected: number;
  /** Maximum density value on the grid */
  maxDensity: number;
  /** Computation time in milliseconds */
  computeTimeMs: number;
}

// ============================================================================
// Population Analysis Result Types (US-076)
// ============================================================================

/**
 * Per-atom charge and population from Mulliken/Lowdin analysis.
 *
 * @see US-076 Mulliken/Lowdin Population Analysis
 */
export interface AtomChargeResult {
  /** Atom index (0-based) */
  atomIndex: number;
  /** Atomic number (Z) */
  atomicNumber: number;
  /** Element symbol (e.g., "O", "H") */
  element: string;
  /** Mulliken charge: q_A = Z_A - N_A */
  mullikenCharge: number;
  /** Lowdin charge */
  lowdinCharge: number;
  /** Mulliken gross atomic population */
  mullikenPopulation: number;
  /** Lowdin atomic population */
  lowdinPopulation: number;
}

/**
 * Result from population_analysis request.
 *
 * Contains per-atom Mulliken and Lowdin charges and populations.
 *
 * @see US-076 Mulliken/Lowdin Population Analysis
 */
export interface PopulationAnalysisResult {
  /** Per-atom charges and populations */
  atoms: AtomChargeResult[];
  /** Total Mulliken charge (sum of all atomic charges) */
  totalMullikenCharge: number;
  /** Total Lowdin charge (sum of all atomic charges) */
  totalLowdinCharge: number;
  /** Computation time in microseconds */
  computeTimeUs: number;
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

// ============================================================================
// Optimization Result Types (US-075)
// ============================================================================

/**
 * A single step in the optimization trajectory.
 *
 * Matches OptimizationStep from qc-core (serde rename_all = "camelCase").
 *
 * @see US-075 Optimization UI + Trajectory Animation
 */
export interface OptimizationStepResult {
  /** Step number (0 = initial evaluation) */
  step: number;
  /** Total energy at this geometry (Ha) */
  energy: number;
  /** Maximum absolute gradient component (Ha/bohr) */
  maxGradient: number;
  /** RMS gradient (Ha/bohr) */
  rmsGradient: number;
  /** Atomic coordinates at this step [[x,y,z], ...] in bohr */
  geometry: [number, number, number][];
  /** Gradient at this step [[gx,gy,gz], ...] in Ha/bohr */
  gradient: [number, number, number][];
}

/**
 * Result of a geometry optimization.
 *
 * Matches OptimizationResult from qc-core (serde rename_all = "camelCase").
 *
 * @see US-075 Optimization UI + Trajectory Animation
 */
export interface OptimizationResult {
  /** Whether the optimization converged */
  converged: boolean;
  /** Full optimization trajectory (including step 0) */
  steps: OptimizationStepResult[];
  /** Final total energy (Ha) */
  finalEnergy: number;
  /** Final optimized geometry [[x,y,z], ...] in bohr */
  finalGeometry: [number, number, number][];
  /** Number of optimization steps taken (excludes initial evaluation) */
  totalSteps: number;
  /** Total computation time in milliseconds */
  computeTimeMs: number;
}

// ============================================================================
// Internal Coordinate PES Scan Result Types (US-081)
// ============================================================================

/**
 * Snapshot of all internal coordinates at a geometry.
 *
 * Contains bonds, angles, and dihedrals measured at the current scan point.
 * Note: Field names use snake_case (matching Rust serde defaults).
 *
 * @see US-081 PES Scan WASM Export + Worker Handler
 */
export interface InternalCoordinateSnapshot {
  /** Bonds as [atom_i, atom_j, distance_bohr] */
  bonds: [number, number, number][];
  /** Angles as [atom_i, atom_j_central, atom_k, angle_radians] */
  angles: [number, number, number, number][];
  /** Dihedrals as [atom_i, atom_j, atom_k, atom_l, angle_radians] */
  dihedrals: [number, number, number, number, number][];
}

/**
 * Single point on the scanned PES (internal coordinates).
 *
 * Note: Field names use snake_case (matching Rust serde defaults).
 *
 * @see US-081 PES Scan WASM Export + Worker Handler
 */
export interface PesInternalPoint {
  /** Value of the scanned coordinate */
  coordinate_value: number;
  /** SCF energy in Hartree */
  energy: number;
  /** Whether SCF converged */
  converged: boolean;
  /** Number of SCF iterations */
  scf_iterations: number;
  /** Cartesian geometry at this point ([x,y,z] per atom) */
  geometry: [number, number, number][];
  /** Optimization steps (null for rigid scans) */
  opt_steps: number | null;
  /** Internal coordinate snapshot (bonds, angles, dihedrals) */
  internal_coordinates: InternalCoordinateSnapshot | null;
}

/**
 * Equilibrium from parabolic interpolation of lowest-energy points.
 *
 * Note: Field names use snake_case (matching Rust serde defaults).
 *
 * @see US-081 PES Scan WASM Export + Worker Handler
 */
export interface PesInternalEquilibrium {
  /** Interpolated equilibrium coordinate value */
  value: number;
  /** Interpolated equilibrium energy in Hartree */
  energy: number;
}

/**
 * Result of an internal coordinate PES scan.
 *
 * Matches the Rust PesScanInternalResult struct serialized via serde.
 * Note: Field names use snake_case as PesScanInternalResult in qc-core
 * does NOT use serde rename_all = "camelCase".
 *
 * @see US-081 PES Scan WASM Export + Worker Handler
 */
export interface PesScanInternalResult {
  /** Type of coordinate scanned */
  coordinate_type: string;
  /** Atom indices defining the scanned coordinate */
  atom_indices: number[];
  /** Scan points ordered by coordinate value */
  points: PesInternalPoint[];
  /** Equilibrium from parabolic interpolation (null if not found) */
  equilibrium: PesInternalEquilibrium | null;
  /** Total SCF iterations across all points */
  total_iterations: number;
  /** Scan mode: "rigid" or "relaxed" */
  scan_mode: string;
  /** Total optimization steps (relaxed only; 0 for rigid) */
  total_opt_steps: number;
}

// ============================================================================
// PES Scan Result Types (US-039)
// ============================================================================

/**
 * A single point on the potential energy surface.
 *
 * Note: Field names match Rust's default serde serialization (snake_case)
 * since PesPoint in qc-core does not use rename_all = "camelCase".
 *
 * @see US-039 PES WASM Export & Worker
 */
export interface PesPoint {
  /** Bond distance in bohr */
  r: number;
  /** Total SCF energy in Hartree (electronic + nuclear) */
  energy: number;
  /** Whether SCF converged at this geometry */
  converged: boolean;
  /** Number of SCF iterations at this point */
  iterations: number;
}

/**
 * Equilibrium geometry from parabolic interpolation.
 *
 * Note: Field names use snake_case (matching Rust serde defaults).
 */
export interface PesEquilibrium {
  /** Interpolated equilibrium bond distance in bohr */
  r_bohr: number;
  /** Interpolated equilibrium energy in Hartree */
  energy_hartree: number;
}

/**
 * Result of a PES bond-length scan.
 *
 * Note: Field names use snake_case (matching Rust serde defaults).
 */
export interface PesScanResult {
  /** Energy at each scanned geometry */
  points: PesPoint[];
  /** Equilibrium from parabolic interpolation (null if not found) */
  equilibrium: PesEquilibrium | null;
  /** Total computation time in milliseconds */
  compute_time_ms: number;
  /** Total SCF iterations across all scan points */
  total_iterations: number;
}

// ============================================================================
// Frequency Analysis Result Types (US-101)
// ============================================================================

/**
 * RRHO thermochemistry reshaped for JavaScript consumption.
 *
 * Flattens `ThermochemistryResult` totals plus the `ThermoComponents`
 * breakdown into a single struct. All energies are in Hartree; entropies
 * and heat capacities are in Ha/(mol·K).
 *
 * // MUST match FrequencyThermochemistry in crates/qc-wasm/src/lib.rs
 *
 * @see US-101 Frequency WASM Export + Worker Handler
 * @see qc_core::thermochemistry::ThermochemistryResult
 */
export interface FrequencyThermochemistry {
  // ---- Input echoes ----
  /** Temperature used, in Kelvin. */
  temperatureK: number;
  /** Pressure used, in Pascals. */
  pressurePa: number;
  /** Rotational symmetry number actually used (σ after defaulting). */
  symmetryNumber: number;
  /** Spin multiplicity 2S+1. */
  multiplicity: number;
  /** Total molecular mass in amu. */
  totalMassAmu: number;
  /** Number of vibrational modes used (positive frequencies only). */
  nVibModesUsed: number;
  /** Number of imaginary modes skipped from the vibrational sum. */
  nImag: number;

  // ---- Zero-point + 0 K ----
  /** Zero-point vibrational energy in Hartree. */
  zpeHa: number;
  /** Energy at 0 K: `E_elec + ZPE` in Hartree. */
  e0kHa: number;

  // ---- Totals ----
  /** Total internal energy `U(T)` in Hartree. */
  internalEnergyHa: number;
  /** Total enthalpy `H(T) = U(T) + RT` in Hartree. */
  enthalpyHa: number;
  /** Total entropy `S(T)` in Ha/(mol·K). */
  entropyHaPerK: number;
  /** Total Gibbs free energy `G(T) = H(T) - T·S(T)` in Hartree. */
  gibbsHa: number;
  /** Total constant-volume heat capacity `Cv(T)` in Ha/(mol·K). */
  cvHaPerK: number;
  /** Total constant-pressure heat capacity `Cp(T) = Cv(T) + R` in Ha/(mol·K). */
  cpHaPerK: number;

  // ---- Translational contributions ----
  eTransHa: number;
  hTransHa: number;
  sTransHaPerK: number;
  cvTransHaPerK: number;
  cpTransHaPerK: number;

  // ---- Rotational contributions ----
  eRotHa: number;
  hRotHa: number;
  sRotHaPerK: number;
  cvRotHaPerK: number;
  cpRotHaPerK: number;

  // ---- Vibrational contributions (thermal only; ZPE stored separately) ----
  eVibThermalHa: number;
  hVibHa: number;
  sVibHaPerK: number;
  cvVibHaPerK: number;
  cpVibHaPerK: number;

  // ---- Electronic contribution (only entropy is nonzero in RRHO) ----
  sElecHaPerK: number;
}

/**
 * Simulated IR or Raman spectrum (broadened on a wavenumber grid).
 *
 * // MUST match FrequencySpectrum in crates/qc-wasm/src/lib.rs
 *
 * @see US-101 Frequency WASM Export + Worker Handler
 */
export interface FrequencySpectrum {
  /** Wavenumber grid in cm⁻¹ (default 0..=4500 step 1, 4501 points). */
  wavenumbersCm1: number[];
  /**
   * Broadened intensity at each grid point. Same length as `wavenumbersCm1`.
   * Units: km/mol for IR, Å⁴/amu for Raman.
   */
  intensity: number[];
  /** Broadening kind used. */
  kind: BroadeningKind;
  /** FWHM used for the broadening in cm⁻¹. */
  fwhmCm1: number;
}

/**
 * Per-phase timing breakdown in milliseconds.
 *
 * // MUST match FrequencyTiming in crates/qc-wasm/src/lib.rs
 *
 * @see US-101 Frequency WASM Export + Worker Handler
 */
export interface FrequencyTiming {
  /** Phase 1: SCF + integrals + Hessian assembly (nuclear CPHF). */
  integralsMs: number;
  /** Phase 2: CPHF data extraction + density rebuild. */
  nuclearCphfMs: number;
  /** Phase 3: Field CPHF (runs inside `compute_raman_spectrum`). */
  fieldCphfMs: number;
  /** Phase 4: Harmonic analysis + IR intensities + thermochemistry. */
  assemblyMs: number;
  /** Phase 5: Spectrum broadening. */
  modesMs: number;
  /** Total wall time from input deserialization to result serialization. */
  totalMs: number;
}

/**
 * Full frequency-analysis result returned by `compute_frequencies`.
 *
 * Aggregates outputs from `rhf_hessian`/`dft_hessian`, `harmonic_analysis`,
 * `compute_ir_spectrum`, `compute_raman_spectrum`, `compute_thermochemistry`,
 * and `simulate_ir_spectrum` / `simulate_raman_spectrum`.
 *
 * // MUST match FrequencyWasmResult in crates/qc-wasm/src/lib.rs
 *
 * @see US-101 Frequency WASM Export + Worker Handler
 */
export interface FrequencyResult {
  // ---- Size metadata ----
  /** Number of atoms. */
  nAtoms: number;
  /** Number of vibrational modes (3N-6 nonlinear, 3N-5 linear, 0 atom). */
  nModes: number;
  /** Rotor classification (one of the 5 `RotorType` variants). */
  rotorType: RotorType;

  // ---- Electronic properties ----
  /** Electronic SCF energy in Hartree (from `HessianResult.energy`). */
  electronicEnergyHa: number;
  /** Equilibrium dipole moment in atomic units (e·bohr), [x, y, z]. */
  dipoleAu: [number, number, number];
  /** Equilibrium dipole moment in Debye, [x, y, z]. */
  dipoleDebye: [number, number, number];
  /** Static polarizability tensor in atomic units (bohr³), symmetric 3×3. */
  polarizabilityAu: [
    [number, number, number],
    [number, number, number],
    [number, number, number],
  ];
  /** Static polarizability tensor in Å³, symmetric 3×3. */
  polarizabilityAng3: [
    [number, number, number],
    [number, number, number],
    [number, number, number],
  ];

  // ---- Vibrational structure ----
  /**
   * Vibrational frequencies in cm⁻¹. Negative values are imaginary
   * (transition states). Length = `nModes`.
   */
  frequenciesCm1: number[];
  /** Reduced masses in amu, one per mode. Length = `nModes`. */
  reducedMassesAmu: number[];
  /** Force constants in mDyne/Å, one per mode. Length = `nModes`. */
  forceConstantsMdyne: number[];
  /**
   * Cartesian normal modes, indexed as `[mode][atom][xyz]`.
   * Shape: `nModes × nAtoms × 3`.
   */
  normalModesCartesian: [number, number, number][][];
  /** Rotational constants in GHz (A ≥ B ≥ C; Infinity for zero moments). */
  rotationalConstantsGhz: [number, number, number];

  // ---- IR ----
  /** IR absorption intensities in km/mol, one per mode. Length = `nModes`. */
  irIntensitiesKmPerMol: number[];

  // ---- Raman ----
  /** Raman scattering activities in Å⁴/amu, one per mode. */
  ramanActivitiesA4Amu: number[];
  /** Depolarization ratios ρ_p (plane-polarized), one per mode. */
  depolarizationRatios: number[];

  // ---- Thermochemistry ----
  /** RRHO thermochemistry at the requested temperature and pressure. */
  thermochemistry: FrequencyThermochemistry;

  // ---- Simulated spectra ----
  /** Continuous broadened IR spectrum on a wavenumber grid. */
  irSpectrum: FrequencySpectrum;
  /** Continuous broadened Raman spectrum on a wavenumber grid. */
  ramanSpectrum: FrequencySpectrum;

  // ---- Metadata ----
  /** Per-phase timings in milliseconds. */
  timingMs: FrequencyTiming;
  /**
   * Whether the calculation was aborted mid-pipeline.
   *
   * The WASM function always sets this to `false`; the worker-side handler
   * sets it to `true` when `isAborted()` was observed during progress events.
   */
  aborted: boolean;
}

// ============================================================================
// MO Grid Result Types (US-042)
// ============================================================================

/**
 * Result of MO grid evaluation.
 *
 * Contains the grid values and metadata for isosurface extraction.
 * Grid values use C-order indexing: index = ix * ny * nz + iy * nz + iz.
 *
 * Note: Field names use camelCase (matching Rust serde rename_all = "camelCase").
 *
 * @see US-042 MO Grid Evaluation
 */
export interface MoGridResult {
  /** Flat array of grid values (C-order: x-slowest, z-fastest) */
  values: number[];
  /** Grid origin [x, y, z] in Bohr */
  gridOrigin: [number, number, number];
  /** Grid spacing in Bohr */
  gridSpacing: number;
  /** Grid dimensions [nx, ny, nz] */
  gridDims: [number, number, number];
  /** Maximum absolute value in the grid */
  maxAbsValue: number;
  /** Approximate norm-squared integral (sum(psi^2) * dV) */
  normSqIntegral: number;
  /** Computation time in milliseconds */
  computeTimeMs: number;
}

// ============================================================================
// Marching Cubes Result Types (US-043)
// ============================================================================

/**
 * Result from marching cubes isosurface extraction.
 *
 * Contains triangle mesh data directly consumable by Three.js BufferGeometry.
 * Vertices and normals use f32 precision (GPU-native); indices use u32.
 *
 * Note: Field names use camelCase (matching Rust serde rename_all = "camelCase").
 *
 * @see US-043 Marching Cubes Isosurface
 */
export interface MarchingCubesResult {
  /** Interleaved vertex positions [x0,y0,z0, x1,y1,z1, ...] */
  vertices: number[];
  /** Triangle indices (3 per triangle) */
  indices: number[];
  /** Interleaved vertex normals [nx0,ny0,nz0, nx1,ny1,nz1, ...] */
  normals: number[];
}

/**
 * Result from dual marching cubes (positive and negative orbital lobes).
 *
 * @see US-043 Marching Cubes Isosurface
 */
export interface DualMarchingCubesResult {
  /** Positive lobe isosurface (field > +isovalue) */
  positive: MarchingCubesResult;
  /** Negative lobe isosurface (field < -isovalue) */
  negative: MarchingCubesResult;
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
      : T extends 'boys_eval_all'
        ? BoysEvalAllResult
        : T extends 'rys_compute'
        ? RysComputeResult
        : T extends 'rys_error_curve'
          ? RysErrorCurveResult
          : T extends 'scf_run'
            ? ScfRunResult
            : T extends 'ks_scf'
              ? KsScfResult
              : T extends 'integral_compute'
                ? IntegralComputeResult
                : T extends 'pes_scan'
                  ? PesScanResult
                  : T extends 'mo_grid'
                    ? MoGridResult
                    : T extends 'marching_cubes'
                      ? MarchingCubesResult
                      : T extends 'dual_marching_cubes'
                        ? DualMarchingCubesResult
                        : T extends 'basis_info'
                          ? BasisInfoResult
                          : T extends 'radial_profile'
                            ? RadialProfileResult
                            : T extends 'overlap_distance'
                              ? OverlapDistanceResult
                              : T extends 'integral_matrices'
                                ? IntegralMatricesResult
                                : T extends 'integral_breakdown'
                                  ? IntegralBreakdownResult
                                  : T extends 'fock_decomposition'
                                    ? FockDecompositionResult
                                    : T extends 'eri_detail'
                                      ? EriDetailResult
                                      : T extends 'density_grid'
                                        ? DensityGridResult
                                        : T extends 'difference_density'
                                          ? DifferenceDensityResult
                                          : T extends 'optimize_geometry'
                                            ? OptimizationResult
                                            : T extends 'pes_scan_internal'
                                              ? PesScanInternalResult
                                              : T extends 'population_analysis'
                                                ? PopulationAnalysisResult
                                                : T extends 'frequency'
                                                  ? FrequencyResult
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
    'boys_eval_all',
    'rys_compute',
    'rys_error_curve',
    'scf_run',
    'ks_scf',
    'integral_compute',
    'pes_scan',
    'mo_grid',
    'marching_cubes',
    'dual_marching_cubes',
    'basis_info',
    'radial_profile',
    'overlap_distance',
    'integral_matrices',
    'integral_breakdown',
    'fock_decomposition',
    'eri_detail',
    'density_grid',
    'difference_density',
    'optimize_geometry',
    'pes_scan_internal',
    'population_analysis',
    'frequency',
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
