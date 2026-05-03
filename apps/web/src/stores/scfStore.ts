/**
 * SCF Module State Store
 *
 * Zustand store managing the state for Module E (SCF Sandbox).
 * Handles parameter inputs, computation status, iteration history,
 * and URL synchronization.
 *
 * Unlike Modules B/C which auto-compute on parameter changes,
 * Module E uses explicit Run/Cancel buttons due to longer computation times.
 *
 * @module stores/scfStore
 */

import { create } from 'zustand';
import { devtools } from 'zustand/middleware';
import type {
  ScfRunResult,
  ScfIterationHistory,
  ConvergenceProfile,
  GeometryInput,
  BasisSetName,
  CoordinateUnits,
  IntegralPhase,
  IntegralComputeResult,
  PesPoint,
  PesEquilibrium,
  PesScanResult,
  PesScanInternalResult,
  OptimizationStepResult,
  OptimizationResult,
  FrequencyResult,
  FrequencyThermochemistry,
  FrequencyProgress,
  BroadeningKind,
} from '../worker/protocol';
import type { EnergyUnitsMode } from '../lib/units';
import type { ScfParams } from '../types/run-state';
import type { PlaneAxis, ColorScale, DensityMode } from '../types/density';
import { DEFAULT_DENSITY_STATE } from '../types/density';
import type { DftMethod } from '../types/dft';
import { DEFAULT_METHOD } from '../types/dft';

// ============================================================================
// Custom Geometry Types
// ============================================================================

/**
 * Input mode for SCF parameters.
 *
 * Determines whether the SCF calculation uses a pre-computed preset system
 * or custom geometry with on-the-fly integral computation.
 */
export type ScfInputMode =
  | { mode: 'preset' }
  | { mode: 'custom'; geometry: GeometryInput; basisSet: BasisSetName };

/**
 * Integral computation status discriminated union.
 *
 * Tracks the progress of on-the-fly integral computation for custom geometries.
 */
export type IntegralComputeStatus =
  | { status: 'idle' }
  | { status: 'computing'; phase: IntegralPhase; progress: number }
  | { status: 'success'; systemData: IntegralComputeResult }
  | { status: 'error'; error: string };

/**
 * SCF computation status discriminated union.
 *
 * Provides type-safe handling of async computation states.
 */
export type ScfComputeStatus =
  | { status: 'idle' }
  | { status: 'running' }
  | { status: 'success'; result: ScfRunResult }
  | { status: 'error'; error: string }
  | { status: 'cancelled' };

/**
 * Display mode for SCF results.
 * - 'explain': Educational explanations
 * - 'internals': Raw matrices and details (US-018)
 */
export type DisplayMode = 'explain' | 'internals';

/**
 * Available system preset info for dropdown.
 */
export interface SystemPreset {
  /** System ID matching worker presets */
  id: string;
  /** Display label */
  label: string;
  /** Brief description */
  description: string;
  /** Number of basis functions */
  nbf: number;
  /** Number of electrons */
  nelec: number;
}

/**
 * Phase 2: 3D viewer state for Module E integration (US-035).
 *
 * Controls the 3D molecular viewer panel: label visibility,
 * orbital selection for isosurface rendering, and isovalue threshold.
 */
export interface ViewerState {
  /** Whether to show atom labels in the 3D viewer */
  showLabels: boolean;
  /** Index of the selected molecular orbital for isosurface display (null = none) */
  selectedOrbital: number | null;
  /** Isovalue for orbital isosurface rendering (atomic units) */
  isovalue: number;

  // === Phase 3: Density visualization state (US-062) ===

  /** Density mode: 'total' electron density or 'difference' (D - D_atomic) */
  densityMode: DensityMode;
  /** Isovalue for density isosurface in e/bohr^3 */
  densityIsovalue: number;
  /**
   * Isovalue for difference density isosurface in e/bohr^3.
   * Separate from densityIsovalue because difference density magnitudes
   * are much smaller than total density (typical: 0.005 vs 0.05).
   * @see US-063 Difference Density (AC3)
   */
  diffIsovalue: number;
  /** Plane axis for cross-section */
  planeAxis: PlaneAxis;
  /** Position along the perpendicular axis (Bohr) */
  planePosition: number;
  /** Color scale for cross-section */
  colorScale: ColorScale;
  /** Whether density visualization is active */
  showDensity: boolean;
}

/**
 * Default viewer state values.
 */
export const DEFAULT_VIEWER_STATE: ViewerState = {
  showLabels: false,
  selectedOrbital: null,
  isovalue: 0.03,
  // Phase 3: Density defaults (US-062)
  ...DEFAULT_DENSITY_STATE,
};

/**
 * Phase 2: PES scan state for Module E integration (US-040).
 *
 * Tracks the configuration, progress, and results of a PES bond-length scan.
 * Results accumulate incrementally as each scan point completes.
 */
export interface PesState {
  /** Whether a PES scan is currently running */
  scanning: boolean;
  /** Minimum bond distance in bohr */
  rMin: number;
  /** Maximum bond distance in bohr */
  rMax: number;
  /** Number of scan points */
  nPoints: number;
  /** Accumulated scan results (grows incrementally during scan) */
  results: PesPoint[];
  /** Equilibrium geometry from parabolic interpolation (set on completion) */
  equilibrium: PesEquilibrium | null;
  /** Error message from the scan (null if no error) */
  error: string | null;
  /** Request ID for the running scan (for cancellation) */
  scanRequestId: string | null;
  /** Total computation time in milliseconds (set on completion) */
  computeTimeMs: number | null;
  /** Total SCF iterations across all scan points (set on completion) */
  totalIterations: number | null;
  /** Full internal coordinate scan result (preserved for coordinate tracking, US-083) */
  pesInternalResult: PesScanInternalResult | null;
  /** Keys of coordinates selected for tracking (US-083) */
  trackedCoordinates: string[];

  // === Phase 4: URL-initialized internal coordinate config (US-085) ===

  /** Coordinate type restored from URL (null if not from URL or legacy format) */
  urlCoordinateType: 'bond' | 'angle' | 'dihedral' | null;
  /** Atom indices restored from URL (null if not from URL) */
  urlAtomIndices: number[] | null;
  /** Scan mode restored from URL (null if not from URL) */
  urlScanMode: 'rigid' | 'relaxed' | null;
  /** Min value restored from URL in display units (degrees for angles) */
  urlValueMin: number | null;
  /** Max value restored from URL in display units */
  urlValueMax: number | null;

  /**
   * Current coordinate configuration synced from PesScanPanel for URL encoding.
   * Updated whenever the user changes coordinate settings in the panel.
   * Null when no internal coordinate config has been set (legacy diatomic mode).
   * @see US-085 Deep Links
   */
  pesCoordinateConfig: {
    coordinateType: 'bond' | 'angle' | 'dihedral';
    atomIndices: number[];
    scanMode: 'rigid' | 'relaxed';
    valueMin: number;
    valueMax: number;
    nPoints: number;
  } | null;
}

/**
 * Default PES scan state values.
 */
export const DEFAULT_PES_STATE: PesState = {
  scanning: false,
  rMin: 0.5,
  rMax: 5.0,
  nPoints: 20,
  results: [],
  equilibrium: null,
  error: null,
  scanRequestId: null,
  computeTimeMs: null,
  totalIterations: null,
  pesInternalResult: null,
  trackedCoordinates: [],
  // Phase 4: URL-initialized internal coordinate config (US-085)
  urlCoordinateType: null,
  urlAtomIndices: null,
  urlScanMode: null,
  urlValueMin: null,
  urlValueMax: null,
  pesCoordinateConfig: null,
};

/**
 * Phase 2: Basis set comparison result for a single basis (US-045).
 *
 * Contains SCF results and metadata for one basis set within
 * a multi-basis comparison run.
 */
export interface BasisComparisonResult {
  /** Basis set name (e.g., "sto-3g", "6-31g*") */
  basisName: string;
  /** Number of basis functions */
  nbf: number;
  /** Final SCF energy in Hartree */
  energy: number;
  /** Whether SCF converged */
  converged: boolean;
  /** Number of SCF iterations */
  iterations: number;
  /** Computation time in milliseconds (integrals + SCF) */
  computeTimeMs: number;
  /** Full iteration history for convergence plotting */
  history: ScfIterationHistory[];
  /** Error message if this basis failed */
  error?: string;
}

/**
 * Phase 2: Basis set comparison state for Module E integration (US-045).
 *
 * Tracks the configuration, progress, and results of a multi-basis
 * comparison run. Results accumulate as each basis completes.
 */
export interface BasisCompareState {
  /** Selected basis sets to compare (2-4) */
  selectedBases: string[];
  /** Whether a comparison run is in progress */
  comparing: boolean;
  /** Accumulated results (grows as each basis completes) */
  results: BasisComparisonResult[];
  /** Progress indicator showing current basis */
  progress: { current: number; total: number; currentBasis: string } | null;
  /** Error message from the comparison run */
  error: string | null;
}

/**
 * Default basis comparison state values.
 */
export const DEFAULT_BASIS_COMPARE_STATE: BasisCompareState = {
  selectedBases: ['sto-3g', '6-31g'],
  comparing: false,
  results: [],
  progress: null,
  error: null,
};

// ============================================================================
// DFT Comparison Types (US-070)
// ============================================================================

/**
 * DFT vs HF comparison state.
 *
 * Tracks the progress and results of a sequential RHF + DFT comparison run.
 *
 * @see US-070 DFT UI + Method Selector + Deep Links (AC4)
 */
export interface DftCompareState {
  /** Whether a comparison run is in progress */
  comparing: boolean;
  /** Which phase of comparison is active */
  phase: 'idle' | 'running_rhf' | 'running_dft' | 'complete' | 'error';
  /** RHF result */
  rhfResult: { energy: number; iterations: number; converged: boolean } | null;
  /** DFT result */
  dftResult: { energy: number; iterations: number; converged: boolean; method: string } | null;
  /** Error message */
  error: string | null;
}

/**
 * Geometry optimization state (US-075).
 *
 * Tracks the progress, results, and trajectory animation state
 * for L-BFGS geometry optimization.
 */
export interface OptimizationState {
  /** Whether optimization is currently running */
  running: boolean;
  /** Accumulated optimization step results (grows during optimization) */
  progress: OptimizationStepResult[];
  /** Final optimization result (set on completion) */
  result: OptimizationResult | null;
  /** Error message from optimization (null if no error) */
  error: string | null;
  /** Request ID for the running optimization (for cancellation) */
  requestId: string | null;
  /** Current trajectory animation step (for slider) */
  trajectoryStep: number;
  /** Whether to show ghost overlay of initial geometry */
  showGhostOverlay: boolean;
}

/**
 * Default optimization state values.
 */
export const DEFAULT_OPTIMIZATION_STATE: OptimizationState = {
  running: false,
  progress: [],
  result: null,
  error: null,
  requestId: null,
  trajectoryStep: 0,
  showGhostOverlay: false,
};

// ============================================================================
// Frequency Analysis State (US-103)
// ============================================================================

/**
 * Phase 5: Frequency analysis state for Module C.
 *
 * Lifted from `LocalFrequencyState` (US-102) into Zustand by US-103.
 * Replaces the React `useState<LocalFrequencyState>` in FrequencyTab.
 *
 * @see US-102 frequencyTabState.ts (original shape)
 * @see US-103 Frequency State + Deep Links
 */
export interface FrequencyState {
  /** Frozen result from the last successful WASM call. null if not yet run. */
  result: FrequencyResult | null;
  /** Whether the WASM call is currently running. */
  isComputing: boolean;
  /** Progress from the most recent phase event. null if not running. */
  progress: FrequencyProgress | null;
  /** Error message from the last WASM call, or null on success / no run. */
  error: string | null;
  /** Currently-highlighted mode index (0-based), or null if no selection. */
  selectedMode: number | null;
  /** User-controlled temperature in K (default 298.15). */
  temperatureK: number;
  /** User-controlled pressure in Pa (default 101325). */
  pressurePa: number;
  /** Thermochemistry recomputed at current (T, P). */
  displayThermo: FrequencyThermochemistry | null;
  /** Currently-visible spectrum tab. */
  spectrumTab: 'ir' | 'raman';
  /** Broadening kernel choice (default 'lorentzian'). */
  broadeningKind: BroadeningKind;
  /** FWHM in cm^-1 for broadening (default 8.0, range 2-20). */
  fwhmCm1: number;
  /** Displacement arrow amplitude (default 0.5, range 0.1-2.0). Renamed from amplitudeAngstrom (US-103). */
  amplitudeBohr: number;
  /** Animation speed multiplier (default 1.0, range 0.5-3.0). */
  animationSpeed: number;
  /** Whether the displacement animation is playing. */
  isAnimating: boolean;
  /** Whether to overlay static displacement-vector arrows. */
  showDisplacementArrows: boolean;
  /** Energy units for the thermochemistry display. */
  unitsMode: EnergyUnitsMode;
}

/**
 * Default frequency state values.
 *
 * Matches US-102 `DEFAULT_LOCAL_FREQUENCY_STATE` with the
 * `amplitudeAngstrom` -> `amplitudeBohr` rename (US-103).
 */
export const DEFAULT_FREQUENCY_STATE: FrequencyState = {
  result: null,
  isComputing: false,
  progress: null,
  error: null,
  selectedMode: null,
  temperatureK: 298.15,
  pressurePa: 101325,
  displayThermo: null,
  spectrumTab: 'ir',
  broadeningKind: 'lorentzian',
  fwhmCm1: 8.0,
  amplitudeBohr: 0.5,
  animationSpeed: 1.0,
  isAnimating: true,
  showDisplacementArrows: false,
  unitsMode: 'kcal_mol',
};

/**
 * Default DFT comparison state values.
 */
export const DEFAULT_DFT_COMPARE_STATE: DftCompareState = {
  comparing: false,
  phase: 'idle',
  rhfResult: null,
  dftResult: null,
  error: null,
};

/**
 * SCF store state interface.
 */
export interface ScfState {
  // === Parameters (user inputs) ===
  /** Selected molecular system ID */
  systemId: string;
  /** Selected molecule preset ID (independent of basis) */
  selectedMolecule: string;
  /** Selected basis set ID (independent of molecule) */
  selectedBasis: BasisSetName;
  /** Convergence profile */
  convergenceProfile: ConvergenceProfile;
  /** Maximum iterations */
  maxIterations: number;
  /** DIIS enabled/disabled */
  useDiis: boolean;
  /**
   * Fock matrix damping factor (0.0 = no damping, 0.7 = heavy damping).
   *
   * Damping mixes the current and previous Fock matrices to stabilize
   * convergence for difficult systems (e.g., diffuse basis sets like 6-31+G*).
   *
   * F_damped = damp * F_old + (1.0 - damp) * F_new
   */
  damp: number;
  /**
   * Use spherical harmonic basis functions (5D, 7F) vs Cartesian (6D, 10F).
   *
   * Spherical harmonics are the quantum chemistry standard because they:
   * - Avoid linear dependencies in larger basis sets
   * - Result in smaller matrices (faster computation)
   * - Produce cleaner orbital analysis (pure angular momentum states)
   *
   * For s and p orbitals, there is no difference. The distinction only
   * matters for d, f, or higher angular momentum functions.
   *
   * @default false (Cartesian for backward compatibility)
   */
  useSpherical: boolean;
  /**
   * DFT numerical integration grid quality.
   *
   * Controls the number of grid points per atom for DFT calculations.
   * Higher quality gives more accurate DFT energies but is slower.
   *
   * - 'standard': ~4,000 points per atom
   * - 'fine': ~8,000 points per atom
   *
   * Only affects DFT methods (LDA, B3LYP, B3LYP-D3BJ); ignored for RHF.
   *
   * @default 'standard'
   */
  gridQuality: 'standard' | 'fine';
  /**
   * Level shift for virtual orbitals (Hartree).
   *
   * Shifts virtual orbital energies up by this amount to stabilize
   * SCF convergence for systems with small HOMO-LUMO gaps.
   *
   * - 0.0: Off (no level shift)
   * - 0.25: Mild shift
   * - 0.5: Standard shift
   * - 1.0: Strong shift
   *
   * @default 0.0
   */
  levelShift: number;
  /** Display mode */
  mode: DisplayMode;

  // === Custom geometry inputs ===
  /** Current input mode (preset or custom) */
  inputMode: ScfInputMode;
  /** Integral computation status (for custom geometry) */
  integralCompute: IntegralComputeStatus;
  /** Raw text input for custom geometry (XYZ format) */
  customGeometryText: string;
  /** Coordinate units for custom geometry */
  customUnits: CoordinateUnits;

  // === Computation state ===
  /** Current compute status */
  compute: ScfComputeStatus;
  /** Live iteration history (updates during computation) */
  history: ScfIterationHistory[];
  /** Integral computation progress (0-100) during pre-SCF phase */
  integralPercent: number;
  /** Current integral step name (e.g., "overlap", "hcore", "eri", "grid") */
  integralStep: string;
  /** Current running request ID (for cancellation) */
  runningRequestId: string | null;
  /** Running integral compute request ID (for cancellation) */
  integralComputeRequestId: string | null;

  // === Phase 2: 3D viewer state (US-035) ===
  /** 3D molecular viewer state */
  viewerState: ViewerState;

  // === Phase 2: PES scan state (US-040) ===
  /** PES bond-length scan state */
  pesState: PesState;

  // === Phase 2: Basis set comparison state (US-045) ===
  /** Multi-basis comparison state */
  basisCompareState: BasisCompareState;

  // === Phase 3: DFT method selection (US-070) ===
  /** Currently selected computation method (RHF or DFT) */
  method: DftMethod;
  /** DFT vs HF comparison state */
  dftCompareState: DftCompareState;
  /** Raw KS-DFT result for energy decomposition display (null if no DFT result) */
  ksResult: import('../worker/protocol').KsScfResult | null;

  // === Phase 3: Geometry optimization state (US-075) ===
  /** Geometry optimization state */
  optimizationState: OptimizationState;

  // === Phase 5: Frequency analysis state (US-103) ===
  /** Frequency analysis state (lifted from component-local useState) */
  frequencyState: FrequencyState;

  // === Phase 3: Population analysis state (US-076) ===
  /** Population analysis result (null if not computed) */
  populationResult: import('../worker/protocol').PopulationAnalysisResult | null;
  /** Whether to show charge labels on 3D viewer atoms */
  showChargeLabels: boolean;
  /**
   * Whether to color atoms by Mulliken charge in the 3D viewer.
   *
   * When enabled, atom sphere colors use a diverging red-white-blue
   * scale based on Mulliken partial charges instead of CPK element colors.
   * Requires populationResult to be available (post-SCF convergence).
   */
  showChargeOverlay: boolean;

  // === URL sync ===
  /** Flag indicating if state has been initialized from URL */
  urlInitialized: boolean;
}

/**
 * SCF store actions interface.
 */
export interface ScfActions {
  // Parameter setters
  setSystemId: (systemId: string) => void;
  /** Set molecule (clears results, auto-selects preset or switches to custom mode) */
  setSelectedMolecule: (id: string) => void;
  /** Set basis set (clears results, auto-selects preset or switches to custom mode) */
  setSelectedBasis: (basis: BasisSetName) => void;
  setConvergenceProfile: (profile: ConvergenceProfile) => void;
  setMaxIterations: (max: number) => void;
  setUseDiis: (useDiis: boolean) => void;
  setDamp: (damp: number) => void;
  setUseSpherical: (useSpherical: boolean) => void;
  /** Set DFT grid quality (only affects DFT methods) */
  setGridQuality: (quality: 'standard' | 'fine') => void;
  /** Set level shift for virtual orbitals (Hartree) */
  setLevelShift: (shift: number) => void;
  setMode: (mode: DisplayMode) => void;

  // URL initialization
  initializeFromURL: (params: ScfParams) => void;

  // Input mode switching
  /** Switch to preset mode (using pre-computed integrals) */
  switchToPresetMode: () => void;
  /** Switch to custom geometry mode */
  switchToCustomMode: () => void;

  // Custom geometry setters
  /** Set the parsed geometry input */
  setCustomGeometry: (geometry: GeometryInput) => void;
  /** Set the raw geometry text (XYZ format) */
  setCustomGeometryText: (text: string) => void;
  /** Set coordinate units for custom geometry */
  setCustomUnits: (units: CoordinateUnits) => void;
  /** Set basis set for custom geometry */
  setCustomBasisSet: (basis: BasisSetName) => void;

  // Integral computation lifecycle
  /** Mark integral computation as started */
  startIntegralCompute: (requestId: string) => void;
  /** Update integral computation progress */
  updateIntegralProgress: (phase: IntegralPhase, progress: number) => void;
  /** Mark integral computation as complete */
  completeIntegralCompute: (result: IntegralComputeResult) => void;
  /** Mark integral computation as failed */
  failIntegralCompute: (error: string) => void;
  /** Cancel running integral computation */
  cancelIntegralCompute: () => void;

  // Computation lifecycle
  startRun: (requestId: string) => void;
  addIteration: (iteration: ScfIterationHistory) => void;
  /** Update pre-SCF integral computation progress */
  updateScfIntegralProgress: (step: string, percent: number) => void;
  setRunResult: (result: ScfRunResult) => void;
  setRunError: (error: string) => void;
  setRunCancelled: () => void;

  // Phase 2: 3D viewer actions (US-035)
  /** Toggle atom label visibility in 3D viewer */
  setShowLabels: (show: boolean) => void;
  /** Select a molecular orbital for isosurface display */
  setSelectedOrbital: (idx: number | null) => void;
  /** Set isovalue for orbital isosurface rendering */
  setIsovalue: (val: number) => void;

  // Phase 3: Density visualization actions (US-062, US-063)
  /** Set density visualization mode (total/difference) */
  setDensityMode: (mode: DensityMode) => void;
  /** Set density isosurface isovalue (clamped to [0.01, 0.50]) */
  setDensityIsovalue: (val: number) => void;
  /** Set difference density isosurface isovalue (clamped to [0.001, 0.05]) */
  setDiffIsovalue: (val: number) => void;
  /** Set cross-section plane axis */
  setPlaneAxis: (axis: PlaneAxis) => void;
  /** Set cross-section plane position along perpendicular axis */
  setPlanePosition: (pos: number) => void;
  /** Set cross-section color scale */
  setColorScale: (scale: ColorScale) => void;
  /** Toggle density visualization active state */
  setShowDensity: (show: boolean) => void;

  // Phase 2: PES scan actions (US-040)
  /** Set PES scan configuration parameters */
  setPesConfig: (config: Partial<{ rMin: number; rMax: number; nPoints: number }>) => void;
  /** Mark PES scan as started */
  startPesScan: (requestId: string) => void;
  /** Cancel a running PES scan */
  cancelPesScan: () => void;
  /** Add a single PES progress point to the results array */
  updatePesProgress: (point: PesPoint) => void;
  /** Mark PES scan as complete with full result */
  completePesScan: (result: PesScanResult) => void;
  /** Mark PES scan as failed with error message */
  failPesScan: (error: string) => void;
  /** Reset PES scan state (clear results, keep config) */
  resetPes: () => void;
  /** Store full internal coordinate scan result for coordinate tracking (US-083) */
  completePesInternalScan: (result: PesScanInternalResult) => void;
  /** Set tracked coordinate keys for coordinate tracking panel (US-083) */
  setTrackedCoordinates: (keys: string[]) => void;
  /** Sync internal coordinate config from PesScanPanel for URL encoding (US-085) */
  setPesCoordinateConfig: (config: PesState['pesCoordinateConfig']) => void;

  // Phase 2: Basis set comparison actions (US-045)
  /** Set the selected basis sets for comparison */
  setBasisCompareSelected: (bases: string[]) => void;
  /** Mark basis comparison as started */
  startBasisCompare: () => void;
  /** Update basis comparison progress indicator */
  updateBasisCompareProgress: (progress: BasisCompareState['progress']) => void;
  /** Add a single basis comparison result */
  addBasisCompareResult: (result: BasisComparisonResult) => void;
  /** Mark basis comparison as complete */
  completeBasisCompare: () => void;
  /** Mark basis comparison as failed */
  failBasisCompare: (error: string) => void;
  /** Reset basis comparison state (clear results, keep selection) */
  resetBasisCompare: () => void;

  // Phase 3: DFT method selection actions (US-070)
  /** Set the computation method (clears previous results) */
  setMethod: (method: DftMethod) => void;
  /** Mark DFT comparison as started */
  startDftCompare: () => void;
  /** Store RHF result from comparison */
  setDftCompareRhfResult: (result: { energy: number; iterations: number; converged: boolean }) => void;
  /** Store DFT result from comparison */
  setDftCompareDftResult: (result: { energy: number; iterations: number; converged: boolean; method: string }) => void;
  /** Mark comparison as complete */
  completeDftCompare: () => void;
  /** Mark comparison as failed */
  failDftCompare: (error: string) => void;
  /** Reset comparison state */
  resetDftCompare: () => void;
  /** Store raw KS-DFT result for energy decomposition display */
  setKsResult: (result: import('../worker/protocol').KsScfResult | null) => void;

  // Phase 3: Geometry optimization actions (US-075)
  /** Mark optimization as started */
  startOptimization: (requestId: string) => void;
  /** Add a single optimization step to the progress array */
  addOptimizationStep: (step: OptimizationStepResult) => void;
  /** Mark optimization as complete with full result */
  completeOptimization: (result: OptimizationResult) => void;
  /** Mark optimization as failed with error message */
  failOptimization: (error: string) => void;
  /** Cancel a running optimization */
  cancelOptimization: () => void;
  /** Set the current trajectory animation step */
  setTrajectoryStep: (step: number) => void;
  /** Toggle ghost overlay of initial geometry */
  setShowGhostOverlay: (show: boolean) => void;
  /** Apply optimized geometry as the active custom geometry */
  applyOptimizedGeometry: (atoms: [number, number, number, number][], basisName: string) => void;
  /** Reset optimization state */
  resetOptimization: () => void;

  // Phase 3: Population analysis actions (US-076)
  /** Set population analysis result */
  setPopulationResult: (result: import('../worker/protocol').PopulationAnalysisResult | null) => void;
  /** Toggle charge labels on 3D viewer */
  setShowChargeLabels: (show: boolean) => void;
  /** Toggle charge-based coloring of atom spheres in 3D viewer */
  setShowChargeOverlay: (show: boolean) => void;

  // Phase 5: Frequency analysis actions (US-103)
  /** Set frequency result */
  setFrequencyResult: (result: FrequencyResult | null) => void;
  /** Set frequency computing state */
  setFrequencyIsComputing: (isComputing: boolean) => void;
  /** Set frequency progress */
  setFrequencyProgress: (progress: FrequencyProgress | null) => void;
  /** Set frequency error */
  setFrequencyError: (error: string | null) => void;
  /** Set selected normal mode index */
  setFrequencySelectedMode: (mode: number | null) => void;
  /** Set thermochemistry temperature */
  setFrequencyTemperature: (temperatureK: number) => void;
  /** Set thermochemistry pressure */
  setFrequencyPressure: (pressurePa: number) => void;
  /** Set display thermochemistry (recomputed at current T, P) */
  setFrequencyDisplayThermo: (thermo: FrequencyThermochemistry | null) => void;
  /** Set spectrum tab (ir or raman) */
  setFrequencySpectrumTab: (tab: 'ir' | 'raman') => void;
  /** Set broadening kernel */
  setFrequencyBroadeningKind: (kind: BroadeningKind) => void;
  /** Set FWHM in cm^-1 */
  setFrequencyFwhmCm1: (fwhm: number) => void;
  /** Set displacement amplitude */
  setFrequencyAmplitudeBohr: (amplitude: number) => void;
  /** Set animation speed multiplier */
  setFrequencyAnimationSpeed: (speed: number) => void;
  /** Set whether animation is playing */
  setFrequencyIsAnimating: (isAnimating: boolean) => void;
  /** Set whether displacement arrows are shown */
  setFrequencyShowDisplacementArrows: (show: boolean) => void;
  /** Set energy units mode */
  setFrequencyUnitsMode: (mode: EnergyUnitsMode) => void;
  /** Reset frequency state to defaults */
  resetFrequencyState: () => void;

  // Reset
  reset: () => void;
  clearHistory: () => void;
}

/**
 * Available system presets for the dropdown.
 *
 * These must match the system IDs in the worker/backend presets.
 */
export const SYSTEM_PRESETS: SystemPreset[] = [
  {
    id: 'h2_sto3g_r1.4',
    label: 'H2 / STO-3G',
    description: 'Hydrogen molecule at equilibrium (R=1.4 bohr)',
    nbf: 2,
    nelec: 2,
  },
  {
    id: 'heh_plus_sto3g',
    label: 'HeH+ / STO-3G',
    description: 'Helium hydride cation',
    nbf: 2,
    nelec: 2,
  },
  {
    id: 'lih_sto3g',
    label: 'LiH / STO-3G',
    description: 'Lithium hydride',
    nbf: 6,
    nelec: 4,
  },
  {
    id: 'h2o_sto3g',
    label: 'H2O / STO-3G',
    description: 'Water molecule',
    nbf: 7,
    nelec: 10,
  },
  {
    id: 'nh3_sto3g',
    label: 'NH3 / STO-3G',
    description: 'Ammonia molecule',
    nbf: 8,
    nelec: 10,
  },
];

// ============================================================================
// Molecule + Basis Set Presets (Independent Selectors)
// ============================================================================

/**
 * A molecule preset for the molecule selector dropdown.
 *
 * Contains the geometry (in bohr) and metadata for each molecule.
 * When combined with a basis set selection, determines how the
 * SCF calculation is configured.
 */
export interface MoleculePreset {
  /** Unique molecule identifier */
  id: string;
  /** Display label with subscript Unicode (e.g., "H\u2082") */
  label: string;
  /** Chemical formula for display */
  formula: string;
  /** Brief description */
  description: string;
  /** Atom positions in bohr */
  atoms: { symbol: string; xyz: [number, number, number] }[];
  /** Molecular charge (0 = neutral, 1 = cation) */
  charge: number;
  /** Total number of electrons */
  nelec: number;
}

/**
 * Molecule presets with geometries in bohr.
 *
 * All geometries match the pre-computed STO-3G preset systems.
 */
export const MOLECULE_PRESETS: MoleculePreset[] = [
  {
    id: 'h2',
    label: 'H\u2082',
    formula: 'H\u2082',
    description: 'Hydrogen molecule',
    charge: 0,
    nelec: 2,
    atoms: [
      { symbol: 'H', xyz: [0, 0, 0] },
      { symbol: 'H', xyz: [0, 0, 1.4] },
    ],
  },
  {
    id: 'heh_plus',
    label: 'HeH\u207A',
    formula: 'HeH\u207A',
    description: 'Helium hydride cation',
    charge: 1,
    nelec: 2,
    atoms: [
      { symbol: 'He', xyz: [0, 0, 0] },
      { symbol: 'H', xyz: [0, 0, 1.4632] },
    ],
  },
  {
    id: 'lih',
    label: 'LiH',
    formula: 'LiH',
    description: 'Lithium hydride',
    charge: 0,
    nelec: 4,
    atoms: [
      { symbol: 'Li', xyz: [0, 0, 0] },
      { symbol: 'H', xyz: [0, 0, 3.015] },
    ],
  },
  {
    id: 'h2o',
    label: 'H\u2082O',
    formula: 'H\u2082O',
    description: 'Water molecule',
    charge: 0,
    nelec: 10,
    atoms: [
      { symbol: 'O', xyz: [0, 0, 0.2217282] },
      { symbol: 'H', xyz: [0, 1.4305447, -0.8869128] },
      { symbol: 'H', xyz: [0, -1.4305447, -0.8869128] },
    ],
  },
  {
    id: 'nh3',
    label: 'NH\u2083',
    formula: 'NH\u2083',
    description: 'Ammonia molecule',
    charge: 0,
    nelec: 10,
    atoms: [
      { symbol: 'N', xyz: [0, 0, 0.219705] },
      { symbol: 'H', xyz: [0, 1.7714918, -0.512645] },
      { symbol: 'H', xyz: [1.5342036, -0.8857459, -0.512645] },
      { symbol: 'H', xyz: [-1.5342036, -0.8857459, -0.512645] },
    ],
  },
  {
    id: 'hf',
    label: 'HF',
    formula: 'HF',
    description: 'Hydrogen fluoride',
    charge: 0,
    nelec: 10,
    atoms: [
      { symbol: 'H', xyz: [0, 0, 0] },
      { symbol: 'F', xyz: [0, 0, 1.7328] },
    ],
  },
  {
    id: 'ch4',
    label: 'CH\u2084',
    formula: 'CH\u2084',
    description: 'Methane (tetrahedral)',
    charge: 0,
    nelec: 10,
    atoms: [
      { symbol: 'C', xyz: [0, 0, 0] },
      { symbol: 'H', xyz: [1.186, 1.186, 1.186] },
      { symbol: 'H', xyz: [-1.186, -1.186, 1.186] },
      { symbol: 'H', xyz: [-1.186, 1.186, -1.186] },
      { symbol: 'H', xyz: [1.186, -1.186, -1.186] },
    ],
  },
  {
    id: 'c2h4',
    label: 'C\u2082H\u2084',
    formula: 'C\u2082H\u2084',
    description: 'Ethylene (\u03C0 bond)',
    charge: 0,
    nelec: 16,
    atoms: [
      { symbol: 'C', xyz: [0, 0, 1.2592] },
      { symbol: 'C', xyz: [0, 0, -1.2592] },
      { symbol: 'H', xyz: [0, 1.7455, 2.3280] },
      { symbol: 'H', xyz: [0, -1.7455, 2.3280] },
      { symbol: 'H', xyz: [0, 1.7455, -2.3280] },
      { symbol: 'H', xyz: [0, -1.7455, -2.3280] },
    ],
  },
  {
    id: 'c2h6',
    label: 'C\u2082H\u2086',
    formula: 'C\u2082H\u2086',
    description: 'Ethane',
    charge: 0,
    nelec: 18,
    atoms: [
      { symbol: 'C', xyz: [0, 0, 1.4508] },
      { symbol: 'C', xyz: [0, 0, -1.4508] },
      { symbol: 'H', xyz: [0, 1.9217, 2.2700] },
      { symbol: 'H', xyz: [1.6641, -0.9609, 2.2700] },
      { symbol: 'H', xyz: [-1.6641, -0.9609, 2.2700] },
      { symbol: 'H', xyz: [0, -1.9217, -2.2700] },
      { symbol: 'H', xyz: [-1.6641, 0.9609, -2.2700] },
      { symbol: 'H', xyz: [1.6641, 0.9609, -2.2700] },
    ],
  },
  {
    id: 'ch3oh',
    label: 'CH\u2083OH',
    formula: 'CH\u2083OH',
    description: 'Methanol',
    charge: 0,
    nelec: 18,
    atoms: [
      { symbol: 'C', xyz: [-1.2712, 0, 0.3316] },
      { symbol: 'O', xyz: [1.1616, 0, -0.2536] },
      { symbol: 'H', xyz: [-1.8256, 1.6856, 1.0332] },
      { symbol: 'H', xyz: [-1.8256, -1.6856, 1.0332] },
      { symbol: 'H', xyz: [-1.8256, 0, -1.6248] },
      { symbol: 'H', xyz: [2.4344, 0, 0.8344] },
    ],
  },
  {
    id: 'co',
    label: 'CO',
    formula: 'CO',
    description: 'Carbon monoxide (triple bond, C\u207B O\u207A dipole)',
    charge: 0,
    nelec: 14,
    atoms: [
      { symbol: 'C', xyz: [0, 0, 0] },
      { symbol: 'O', xyz: [0, 0, 2.1322] },
    ],
  },
  {
    id: 'hcn',
    label: 'HCN',
    formula: 'HCN',
    description: 'Hydrogen cyanide (triple bond)',
    charge: 0,
    nelec: 14,
    atoms: [
      { symbol: 'H', xyz: [0, 0, -1.9106] },
      { symbol: 'C', xyz: [0, 0, 0] },
      { symbol: 'N', xyz: [0, 0, 2.1793] },
    ],
  },
  {
    id: 'co2',
    label: 'CO\u2082',
    formula: 'CO\u2082',
    description: 'Carbon dioxide (linear)',
    charge: 0,
    nelec: 22,
    atoms: [
      { symbol: 'C', xyz: [0, 0, 0] },
      { symbol: 'O', xyz: [0, 0, 2.1944] },
      { symbol: 'O', xyz: [0, 0, -2.1944] },
    ],
  },
  {
    id: 'hcl',
    label: 'HCl',
    formula: 'HCl',
    description: 'Hydrogen chloride',
    charge: 0,
    nelec: 18,
    atoms: [
      { symbol: 'H', xyz: [0, 0, 0] },
      { symbol: 'Cl', xyz: [0, 0, 2.4086] },
    ],
  },
  {
    id: 'h2s',
    label: 'H\u2082S',
    formula: 'H\u2082S',
    description: 'Hydrogen sulfide',
    charge: 0,
    nelec: 18,
    atoms: [
      { symbol: 'S', xyz: [0, 0, 0.2037] },
      { symbol: 'H', xyz: [0, 1.8273, -0.8149] },
      { symbol: 'H', xyz: [0, -1.8273, -0.8149] },
    ],
  },
  {
    id: 'nacl',
    label: 'NaCl',
    formula: 'NaCl',
    description: 'Sodium chloride (ionic bond)',
    charge: 0,
    nelec: 28,
    atoms: [
      { symbol: 'Na', xyz: [0, 0, 0] },
      { symbol: 'Cl', xyz: [0, 0, 4.4610] },
    ],
  },
  {
    id: 'c6h6',
    label: 'C\u2086H\u2086',
    formula: 'C\u2086H\u2086',
    description: 'Benzene (aromatic \u03C0 system)',
    charge: 0,
    nelec: 42,
    atoms: [
      { symbol: 'C', xyz: [2.640000, 0.000000, 0] },
      { symbol: 'H', xyz: [4.688300, 0.000000, 0] },
      { symbol: 'C', xyz: [1.320000, 2.286307, 0] },
      { symbol: 'H', xyz: [2.344150, 4.060187, 0] },
      { symbol: 'C', xyz: [-1.320000, 2.286307, 0] },
      { symbol: 'H', xyz: [-2.344150, 4.060187, 0] },
      { symbol: 'C', xyz: [-2.640000, 0.000000, 0] },
      { symbol: 'H', xyz: [-4.688300, 0.000000, 0] },
      { symbol: 'C', xyz: [-1.320000, -2.286307, 0] },
      { symbol: 'H', xyz: [-2.344150, -4.060187, 0] },
      { symbol: 'C', xyz: [1.320000, -2.286307, 0] },
      { symbol: 'H', xyz: [2.344150, -4.060187, 0] },
    ],
  },
];

/**
 * Basis set option for the basis set selector dropdown.
 */
export interface BasisSetOption {
  /** Basis set identifier matching BasisSetName */
  id: BasisSetName;
  /** Display label */
  label: string;
}

/**
 * Available basis sets for the basis set selector dropdown.
 */
export const BASIS_SET_OPTIONS: BasisSetOption[] = [
  { id: 'sto-3g', label: 'STO-3G' },
  { id: '3-21g', label: '3-21G' },
  { id: '6-31g', label: '6-31G' },
  { id: '6-31g*', label: '6-31G*' },
  { id: '6-31+g*', label: '6-31+G*' },
  { id: 'cc-pvdz', label: 'cc-pVDZ' },
];

/**
 * Mapping from (moleculeId, basisId) to pre-computed preset system ID.
 *
 * Only combinations that have pre-computed integral JSON files are included.
 * All other combinations require on-the-fly integral computation.
 */
const MOLECULE_BASIS_TO_PRESET: Record<string, string> = {
  'h2:sto-3g': 'h2_sto3g_r1.4',
  'heh_plus:sto-3g': 'heh_plus_sto3g',
  'lih:sto-3g': 'lih_sto3g',
  'h2o:sto-3g': 'h2o_sto3g',
  'nh3:sto-3g': 'nh3_sto3g',
};

/**
 * Reverse mapping from preset system ID to (moleculeId, basisId).
 */
const PRESET_TO_MOLECULE_BASIS: Record<string, { molecule: string; basis: BasisSetName }> = {
  'h2_sto3g_r1.4': { molecule: 'h2', basis: 'sto-3g' },
  'heh_plus_sto3g': { molecule: 'heh_plus', basis: 'sto-3g' },
  'lih_sto3g': { molecule: 'lih', basis: 'sto-3g' },
  'h2o_sto3g': { molecule: 'h2o', basis: 'sto-3g' },
  'nh3_sto3g': { molecule: 'nh3', basis: 'sto-3g' },
};

/**
 * Look up the pre-computed preset ID for a molecule + basis combination.
 *
 * @returns The preset system ID if a pre-computed preset exists, undefined otherwise.
 */
export function getPresetIdForMoleculeBasis(moleculeId: string, basisId: string): string | undefined {
  return MOLECULE_BASIS_TO_PRESET[`${moleculeId}:${basisId}`];
}

/**
 * Convert molecule atoms to XYZ text format for the custom geometry editor.
 */
function moleculeToXyzText(atoms: { symbol: string; xyz: [number, number, number] }[]): string {
  return atoms
    .map((a) => `${a.symbol}  ${a.xyz[0].toFixed(6)}  ${a.xyz[1].toFixed(6)}  ${a.xyz[2].toFixed(6)}`)
    .join('\n');
}

/**
 * Get molecule preset by ID.
 */
export function getMoleculePreset(id: string): MoleculePreset | undefined {
  return MOLECULE_PRESETS.find((m) => m.id === id);
}

/**
 * Convergence profile threshold descriptions for tooltips.
 */
export const CONVERGENCE_THRESHOLDS: Record<ConvergenceProfile, { energy: string; density: string }> = {
  loose: { energy: '1e-4 Ha', density: '1e-4' },
  medium: { energy: '1e-6 Ha', density: '1e-6' },
  tight: { energy: '1e-8 Ha', density: '1e-8' },
};

/**
 * Default basis set for custom geometry mode.
 */
export const DEFAULT_CUSTOM_BASIS: BasisSetName = 'sto-3g';

/**
 * Default state values matching DEFAULT_SCF_STATE from run-state.ts
 */
const DEFAULT_STATE: ScfState = {
  systemId: 'h2_sto3g_r1.4',
  selectedMolecule: 'h2',
  selectedBasis: 'sto-3g',
  convergenceProfile: 'medium',
  maxIterations: 50,
  useDiis: true,
  damp: 0.0,
  useSpherical: false, // Cartesian for backward compatibility
  gridQuality: 'standard',
  levelShift: 0.0,
  mode: 'explain',
  // Custom geometry defaults
  inputMode: { mode: 'preset' },
  integralCompute: { status: 'idle' },
  customGeometryText: '',
  customUnits: 'angstrom',
  // Computation state
  compute: { status: 'idle' },
  history: [],
  integralPercent: 0,
  integralStep: '',
  runningRequestId: null,
  integralComputeRequestId: null,
  // Phase 2: 3D viewer state (US-035)
  viewerState: { ...DEFAULT_VIEWER_STATE },
  // Phase 2: PES scan state (US-040)
  pesState: { ...DEFAULT_PES_STATE },
  // Phase 2: Basis set comparison state (US-045)
  basisCompareState: { ...DEFAULT_BASIS_COMPARE_STATE },
  // Phase 3: DFT method selection (US-070)
  method: DEFAULT_METHOD,
  dftCompareState: { ...DEFAULT_DFT_COMPARE_STATE },
  ksResult: null,
  // Phase 3: Geometry optimization (US-075)
  optimizationState: { ...DEFAULT_OPTIMIZATION_STATE },
  // Phase 5: Frequency analysis (US-103)
  frequencyState: { ...DEFAULT_FREQUENCY_STATE },
  populationResult: null,
  showChargeLabels: false,
  showChargeOverlay: false,
  urlInitialized: false,
};

/**
 * SCF state store.
 *
 * Manages parameter inputs, computation status, iteration history,
 * and URL synchronization for Module E.
 *
 * @example
 * ```typescript
 * import { useScfStore, SYSTEM_PRESETS } from '@/stores/scfStore';
 *
 * function ScfControls() {
 *   const { systemId, setSystemId, compute, history } = useScfStore();
 *
 *   return (
 *     <>
 *       <select value={systemId} onChange={e => setSystemId(e.target.value)}>
 *         {SYSTEM_PRESETS.map(p => (
 *           <option key={p.id} value={p.id}>{p.label}</option>
 *         ))}
 *       </select>
 *       {compute.status === 'running' && <p>Running... {history.length} iterations</p>}
 *     </>
 *   );
 * }
 * ```
 */
export const useScfStore = create<ScfState & ScfActions>()(
  devtools(
    (set, get) => ({
      // Initial state
      ...DEFAULT_STATE,

      // Actions
      setSystemId: (systemId: string) => {
        set(
          (state) => ({
            systemId,
            // Reset orbital selection when switching systems -- old MO data is invalidated
            viewerState: {
              ...state.viewerState,
              selectedOrbital: null,
              isovalue: DEFAULT_VIEWER_STATE.isovalue,
            },
          }),
          false,
          'setSystemId'
        );
      },

      setSelectedMolecule: (moleculeId: string) => {
        const { selectedBasis } = get();
        const molecule = getMoleculePreset(moleculeId);
        if (!molecule) return;

        const presetId = getPresetIdForMoleculeBasis(moleculeId, selectedBasis);

        if (presetId) {
          // Pre-computed preset exists: use fast preset path
          set(
            (state) => ({
              selectedMolecule: moleculeId,
              systemId: presetId,
              inputMode: { mode: 'preset' as const },
              // Clear previous results
              compute: { status: 'idle' as const },
              history: [],
              integralCompute: { status: 'idle' as const },
              populationResult: null,
              ksResult: null,
              viewerState: {
                ...state.viewerState,
                selectedOrbital: null,
                isovalue: DEFAULT_VIEWER_STATE.isovalue,
              },
            }),
            false,
            'setSelectedMolecule'
          );
        } else {
          // No preset: switch to custom mode with geometry from molecule definition
          const geometry: GeometryInput = {
            atoms: molecule.atoms.map((a) => ({ symbol: a.symbol, xyz: a.xyz })),
            units: 'bohr',
          };
          set(
            (state) => ({
              selectedMolecule: moleculeId,
              inputMode: {
                mode: 'custom' as const,
                geometry,
                basisSet: selectedBasis,
              },
              customGeometryText: moleculeToXyzText(molecule.atoms),
              customUnits: 'bohr',
              // Clear previous results
              compute: { status: 'idle' as const },
              history: [],
              integralCompute: { status: 'idle' as const },
              populationResult: null,
              ksResult: null,
              viewerState: {
                ...state.viewerState,
                selectedOrbital: null,
                isovalue: DEFAULT_VIEWER_STATE.isovalue,
              },
            }),
            false,
            'setSelectedMolecule'
          );
        }
      },

      setSelectedBasis: (basisId: BasisSetName) => {
        const { selectedMolecule } = get();
        const molecule = getMoleculePreset(selectedMolecule);

        const presetId = getPresetIdForMoleculeBasis(selectedMolecule, basisId);

        if (presetId) {
          // Pre-computed preset exists: use fast preset path
          set(
            (state) => ({
              selectedBasis: basisId,
              systemId: presetId,
              inputMode: { mode: 'preset' as const },
              // Clear previous results
              compute: { status: 'idle' as const },
              history: [],
              integralCompute: { status: 'idle' as const },
              populationResult: null,
              ksResult: null,
              viewerState: {
                ...state.viewerState,
                selectedOrbital: null,
                isovalue: DEFAULT_VIEWER_STATE.isovalue,
              },
            }),
            false,
            'setSelectedBasis'
          );
        } else if (molecule) {
          // No preset: switch to custom mode with geometry from molecule definition
          const geometry: GeometryInput = {
            atoms: molecule.atoms.map((a) => ({ symbol: a.symbol, xyz: a.xyz })),
            units: 'bohr',
          };
          set(
            (state) => ({
              selectedBasis: basisId,
              inputMode: {
                mode: 'custom' as const,
                geometry,
                basisSet: basisId,
              },
              customGeometryText: moleculeToXyzText(molecule.atoms),
              customUnits: 'bohr',
              // Clear previous results
              compute: { status: 'idle' as const },
              history: [],
              integralCompute: { status: 'idle' as const },
              populationResult: null,
              ksResult: null,
              viewerState: {
                ...state.viewerState,
                selectedOrbital: null,
                isovalue: DEFAULT_VIEWER_STATE.isovalue,
              },
            }),
            false,
            'setSelectedBasis'
          );
        }
      },

      setConvergenceProfile: (convergenceProfile: ConvergenceProfile) => {
        set({ convergenceProfile }, false, 'setConvergenceProfile');
      },

      setMaxIterations: (maxIterations: number) => {
        // Clamp to valid range [1, 200]
        const clamped = Math.max(1, Math.min(200, Math.round(maxIterations)));
        set({ maxIterations: clamped }, false, 'setMaxIterations');
      },

      setUseDiis: (useDiis: boolean) => {
        set({ useDiis }, false, 'setUseDiis');
      },

      setDamp: (damp: number) => {
        // Clamp to valid range [0.0, 0.9]
        const clamped = Math.max(0.0, Math.min(0.9, damp));
        set({ damp: clamped }, false, 'setDamp');
      },

      setUseSpherical: (useSpherical: boolean) => {
        set(
          (state) => ({
            useSpherical,
            // Reset integral compute status when basis type changes (for custom mode)
            integralCompute: state.inputMode.mode === 'custom' ? { status: 'idle' } : state.integralCompute,
          }),
          false,
          'setUseSpherical'
        );
      },

      setGridQuality: (quality: 'standard' | 'fine') => {
        set({ gridQuality: quality }, false, 'setGridQuality');
      },

      setLevelShift: (shift: number) => {
        set({ levelShift: shift }, false, 'setLevelShift');
      },

      setMode: (mode: DisplayMode) => {
        set({ mode }, false, 'setMode');
      },

      initializeFromURL: (params: ScfParams) => {
        // Handle custom mode initialization from URL
        const dampValue = params.damp !== undefined ? Math.max(0.0, Math.min(0.9, params.damp)) : 0.0;

        // Restore method from URL (US-070); fall back to 'rhf' if unavailable
        const methodFromUrl = params.method ?? DEFAULT_METHOD;
        // If method is b3lyp-d3bj (unavailable), fall back to 'rhf'
        const method: DftMethod = methodFromUrl === 'b3lyp-d3bj' ? DEFAULT_METHOD : methodFromUrl;

        // Restore molecule + basis from URL, or infer from system_id
        const moleculeFromUrl = params.molecule;
        const basisFromUrl = params.basis;
        let selectedMolecule = DEFAULT_STATE.selectedMolecule;
        let selectedBasis = DEFAULT_STATE.selectedBasis;
        if (moleculeFromUrl && getMoleculePreset(moleculeFromUrl)) {
          selectedMolecule = moleculeFromUrl;
        }
        if (basisFromUrl && BASIS_SET_OPTIONS.some((b) => b.id === basisFromUrl)) {
          selectedBasis = basisFromUrl as BasisSetName;
        }
        // If molecule/basis not in URL, try to infer from system_id via reverse mapping
        if (!moleculeFromUrl || !basisFromUrl) {
          const inferred = PRESET_TO_MOLECULE_BASIS[params.system_id];
          if (inferred) {
            if (!moleculeFromUrl) selectedMolecule = inferred.molecule;
            if (!basisFromUrl) selectedBasis = inferred.basis;
          }
        }

        // Restore PES config from URL if present (US-085: handles both legacy and internal coordinate formats)
        const pesUpdate = params.pes
          ? {
              pesState: {
                ...DEFAULT_PES_STATE,
                rMin: params.pes.r_min,
                rMax: params.pes.r_max,
                nPoints: params.pes.n_points,
                // Phase 4: restore internal coordinate config if present (US-085)
                ...(params.pes.coordinate_type
                  ? {
                      urlCoordinateType: params.pes.coordinate_type,
                      urlAtomIndices: params.pes.atom_indices ?? null,
                      urlScanMode: params.pes.scan_mode ?? null,
                      urlValueMin: params.pes.r_min,
                      urlValueMax: params.pes.r_max,
                    }
                  : {}),
              },
            }
          : {};

        // Restore orbital + density viewer state from URL if present (US-044, US-062)
        const hasDensityParams =
          params.density_mode !== undefined ||
          params.density_isovalue !== undefined ||
          params.diff_isovalue !== undefined ||
          params.plane_axis !== undefined ||
          params.plane_position !== undefined ||
          params.color_scale !== undefined;

        const hasViewerParams =
          params.orbital !== undefined ||
          params.isovalue !== undefined ||
          hasDensityParams;

        const viewerUpdate = hasViewerParams
          ? {
              viewerState: {
                ...DEFAULT_VIEWER_STATE,
                selectedOrbital: params.orbital ?? null,
                isovalue: params.isovalue ?? DEFAULT_VIEWER_STATE.isovalue,
                // Density state from URL (US-062, US-063)
                densityMode: params.density_mode ?? DEFAULT_VIEWER_STATE.densityMode,
                densityIsovalue: params.density_isovalue ?? DEFAULT_VIEWER_STATE.densityIsovalue,
                diffIsovalue: params.diff_isovalue ?? DEFAULT_VIEWER_STATE.diffIsovalue,
                planeAxis: params.plane_axis ?? DEFAULT_VIEWER_STATE.planeAxis,
                planePosition: params.plane_position ?? DEFAULT_VIEWER_STATE.planePosition,
                colorScale: params.color_scale ?? DEFAULT_VIEWER_STATE.colorScale,
                showDensity: hasDensityParams,
              },
            }
          : {};

        // Restore basis comparison selection from URL if present (US-045)
        const basisCompareUpdate = params.comparison_bases
          ? {
              basisCompareState: {
                ...DEFAULT_BASIS_COMPARE_STATE,
                selectedBases: params.comparison_bases,
              },
            }
          : {};

        // Restore frequency state from URL if present (US-103)
        const hasFreqParams =
          params.freq_mode !== undefined ||
          params.freq_temp !== undefined ||
          params.freq_pres !== undefined ||
          params.freq_broad !== undefined ||
          params.freq_fwhm !== undefined ||
          params.freq_tab !== undefined ||
          params.freq_arrows !== undefined ||
          params.freq_units !== undefined;

        const freqUpdate = hasFreqParams
          ? {
              frequencyState: {
                ...DEFAULT_FREQUENCY_STATE,
                selectedMode: params.freq_mode ?? null,
                temperatureK: params.freq_temp ?? DEFAULT_FREQUENCY_STATE.temperatureK,
                pressurePa: params.freq_pres ?? DEFAULT_FREQUENCY_STATE.pressurePa,
                broadeningKind: params.freq_broad === 'g' ? 'gaussian' as const : 'lorentzian' as const,
                fwhmCm1: params.freq_fwhm ?? DEFAULT_FREQUENCY_STATE.fwhmCm1,
                spectrumTab: params.freq_tab ?? DEFAULT_FREQUENCY_STATE.spectrumTab,
                showDisplacementArrows: params.freq_arrows ?? DEFAULT_FREQUENCY_STATE.showDisplacementArrows,
                unitsMode:
                  params.freq_units === 'ha' ? 'hartree' as const :
                  params.freq_units === 'kj' ? 'kj_mol' as const :
                  'kcal_mol' as const,
              },
            }
          : {};

        if (params.input_mode === 'custom' && params.custom_geometry && params.custom_basis) {
          const geometry: GeometryInput = {
            atoms: params.custom_geometry.atoms.map((a) => ({
              symbol: a.symbol,
              xyz: a.xyz,
            })),
            units: params.custom_geometry.units,
          };
          set(
            {
              systemId: params.system_id,
              selectedMolecule,
              selectedBasis,
              convergenceProfile: params.conv,
              maxIterations: Math.max(1, Math.min(200, params.max_iter)),
              useDiis: params.diis,
              damp: dampValue,
              method,
              inputMode: { mode: 'custom', geometry, basisSet: params.custom_basis as BasisSetName },
              customUnits: params.custom_geometry.units,
              urlInitialized: true,
              ...pesUpdate,
              ...viewerUpdate,
              ...basisCompareUpdate,
              ...freqUpdate,
            },
            false,
            'initializeFromURL'
          );
        } else {
          // Preset mode (default)
          set(
            {
              systemId: params.system_id,
              selectedMolecule,
              selectedBasis,
              convergenceProfile: params.conv,
              maxIterations: Math.max(1, Math.min(200, params.max_iter)),
              useDiis: params.diis,
              damp: dampValue,
              method,
              inputMode: { mode: 'preset' },
              urlInitialized: true,
              ...pesUpdate,
              ...viewerUpdate,
              ...basisCompareUpdate,
              ...freqUpdate,
            },
            false,
            'initializeFromURL'
          );
        }
      },

      // === Input Mode Switching ===

      switchToPresetMode: () => {
        const { selectedMolecule, selectedBasis } = get();
        const presetId = getPresetIdForMoleculeBasis(selectedMolecule, selectedBasis);

        if (presetId) {
          // Current molecule+basis has a preset
          set(
            {
              inputMode: { mode: 'preset' },
              integralCompute: { status: 'idle' },
              systemId: presetId,
            },
            false,
            'switchToPresetMode'
          );
        } else {
          // No preset for current combo: reset basis to STO-3G
          const fallbackPresetId = getPresetIdForMoleculeBasis(selectedMolecule, 'sto-3g');
          set(
            {
              inputMode: { mode: 'preset' },
              integralCompute: { status: 'idle' },
              selectedBasis: 'sto-3g',
              systemId: fallbackPresetId ?? 'h2_sto3g_r1.4',
            },
            false,
            'switchToPresetMode'
          );
        }
      },

      switchToCustomMode: () => {
        set(
          (state) => ({
            inputMode: {
              mode: 'custom',
              geometry: { atoms: [], units: state.customUnits },
              basisSet: DEFAULT_CUSTOM_BASIS,
            },
          }),
          false,
          'switchToCustomMode'
        );
      },

      // === Custom Geometry Setters ===

      setCustomGeometry: (geometry: GeometryInput) => {
        set(
          (state) => {
            if (state.inputMode.mode === 'custom') {
              return {
                inputMode: { ...state.inputMode, geometry },
                // Reset integral compute status when geometry changes
                integralCompute: { status: 'idle' },
              };
            }
            return {};
          },
          false,
          'setCustomGeometry'
        );
      },

      setCustomGeometryText: (text: string) => {
        set({ customGeometryText: text }, false, 'setCustomGeometryText');
      },

      setCustomUnits: (units: CoordinateUnits) => {
        set(
          (state) => {
            if (state.inputMode.mode === 'custom') {
              return {
                customUnits: units,
                inputMode: {
                  ...state.inputMode,
                  geometry: { ...state.inputMode.geometry, units },
                },
                // Reset integral compute status when units change
                integralCompute: { status: 'idle' },
              };
            }
            return { customUnits: units };
          },
          false,
          'setCustomUnits'
        );
      },

      setCustomBasisSet: (basis: BasisSetName) => {
        set(
          (state) => {
            if (state.inputMode.mode === 'custom') {
              return {
                inputMode: { ...state.inputMode, basisSet: basis },
                // Reset integral compute status when basis changes
                integralCompute: { status: 'idle' },
              };
            }
            return {};
          },
          false,
          'setCustomBasisSet'
        );
      },

      // === Integral Computation Lifecycle ===

      startIntegralCompute: (requestId: string) => {
        set(
          {
            integralCompute: { status: 'computing', phase: 'overlap', progress: 0 },
            integralComputeRequestId: requestId,
          },
          false,
          'startIntegralCompute'
        );
      },

      updateIntegralProgress: (phase: IntegralPhase, progress: number) => {
        set(
          {
            integralCompute: { status: 'computing', phase, progress },
          },
          false,
          'updateIntegralProgress'
        );
      },

      completeIntegralCompute: (result: IntegralComputeResult) => {
        set(
          {
            integralCompute: { status: 'success', systemData: result },
            integralComputeRequestId: null,
            // Update systemId to the generated custom system ID
            systemId: result.systemId,
          },
          false,
          'completeIntegralCompute'
        );
      },

      failIntegralCompute: (error: string) => {
        set(
          {
            integralCompute: { status: 'error', error },
            integralComputeRequestId: null,
          },
          false,
          'failIntegralCompute'
        );
      },

      cancelIntegralCompute: () => {
        set(
          {
            integralCompute: { status: 'idle' },
            integralComputeRequestId: null,
          },
          false,
          'cancelIntegralCompute'
        );
      },

      // === SCF Computation Lifecycle ===

      startRun: (requestId: string) => {
        set(
          (state) => ({
            compute: { status: 'running' },
            history: [],
            integralPercent: 0,
            integralStep: '',
            runningRequestId: requestId,
            // Clear orbital selection -- old MO data is invalidated by new SCF run
            viewerState: {
              ...state.viewerState,
              selectedOrbital: null,
            },
          }),
          false,
          'startRun'
        );
      },

      addIteration: (iteration: ScfIterationHistory) => {
        set(
          (state) => ({
            history: [...state.history, iteration],
          }),
          false,
          'addIteration'
        );
      },

      updateScfIntegralProgress: (step: string, percent: number) => {
        set(
          { integralStep: step, integralPercent: percent },
          false,
          'updateScfIntegralProgress'
        );
      },

      setRunResult: (result: ScfRunResult) => {
        // If history is empty (no streaming progress, e.g. DFT path),
        // populate from result's history so convergence plots render.
        const currentHistory = get().history;
        const history =
          currentHistory.length === 0 && result.history && result.history.length > 0
            ? result.history
            : currentHistory;
        set(
          {
            compute: { status: 'success', result },
            runningRequestId: null,
            integralPercent: 0,
            integralStep: '',
            history,
          },
          false,
          'setRunResult'
        );
      },

      setRunError: (error: string) => {
        set(
          {
            compute: { status: 'error', error },
            runningRequestId: null,
            integralPercent: 0,
            integralStep: '',
          },
          false,
          'setRunError'
        );
      },

      setRunCancelled: () => {
        set(
          {
            compute: { status: 'cancelled' },
            runningRequestId: null,
            integralPercent: 0,
            integralStep: '',
          },
          false,
          'setRunCancelled'
        );
      },

      // === Phase 2: 3D Viewer Actions (US-035) ===

      setShowLabels: (show: boolean) => {
        set(
          (state) => ({
            viewerState: { ...state.viewerState, showLabels: show },
          }),
          false,
          'setShowLabels'
        );
      },

      setSelectedOrbital: (idx: number | null) => {
        set(
          (state) => ({
            viewerState: { ...state.viewerState, selectedOrbital: idx },
          }),
          false,
          'setSelectedOrbital'
        );
      },

      setIsovalue: (val: number) => {
        // Clamp to valid range [0.005, 0.10]
        const clamped = Math.max(0.005, Math.min(0.10, val));
        set(
          (state) => ({
            viewerState: { ...state.viewerState, isovalue: clamped },
          }),
          false,
          'setIsovalue'
        );
      },

      // === Phase 3: Density Visualization Actions (US-062) ===

      setDensityMode: (mode: DensityMode) => {
        set(
          (state) => ({
            viewerState: { ...state.viewerState, densityMode: mode },
          }),
          false,
          'setDensityMode'
        );
      },

      setDensityIsovalue: (val: number) => {
        // Clamp to valid range [0.01, 0.50]
        const clamped = Math.max(0.01, Math.min(0.50, val));
        set(
          (state) => ({
            viewerState: { ...state.viewerState, densityIsovalue: clamped },
          }),
          false,
          'setDensityIsovalue'
        );
      },

      setDiffIsovalue: (val: number) => {
        // Clamp to valid range [0.001, 0.05] (AC3)
        const clamped = Math.max(0.001, Math.min(0.05, val));
        set(
          (state) => ({
            viewerState: { ...state.viewerState, diffIsovalue: clamped },
          }),
          false,
          'setDiffIsovalue'
        );
      },

      setPlaneAxis: (axis: PlaneAxis) => {
        set(
          (state) => ({
            viewerState: { ...state.viewerState, planeAxis: axis },
          }),
          false,
          'setPlaneAxis'
        );
      },

      setPlanePosition: (pos: number) => {
        set(
          (state) => ({
            viewerState: { ...state.viewerState, planePosition: pos },
          }),
          false,
          'setPlanePosition'
        );
      },

      setColorScale: (scale: ColorScale) => {
        set(
          (state) => ({
            viewerState: { ...state.viewerState, colorScale: scale },
          }),
          false,
          'setColorScale'
        );
      },

      setShowDensity: (show: boolean) => {
        set(
          (state) => ({
            viewerState: {
              ...state.viewerState,
              showDensity: show,
              // When enabling density, deselect orbital to avoid visual clutter
              selectedOrbital: show ? null : state.viewerState.selectedOrbital,
            },
          }),
          false,
          'setShowDensity'
        );
      },

      // === Phase 2: PES Scan Actions (US-040) ===

      setPesConfig: (config) => {
        set(
          (state) => ({
            pesState: { ...state.pesState, ...config },
          }),
          false,
          'setPesConfig'
        );
      },

      startPesScan: (requestId: string) => {
        set(
          (state) => ({
            pesState: {
              ...state.pesState,
              scanning: true,
              results: [],
              equilibrium: null,
              error: null,
              scanRequestId: requestId,
              pesInternalResult: null,
              trackedCoordinates: [],
            },
          }),
          false,
          'startPesScan'
        );
      },

      cancelPesScan: () => {
        set(
          (state) => ({
            pesState: {
              ...state.pesState,
              scanning: false,
              scanRequestId: null,
            },
          }),
          false,
          'cancelPesScan'
        );
      },

      updatePesProgress: (point: PesPoint) => {
        set(
          (state) => ({
            pesState: {
              ...state.pesState,
              results: [...state.pesState.results, point],
            },
          }),
          false,
          'updatePesProgress'
        );
      },

      completePesScan: (result: PesScanResult) => {
        set(
          (state) => ({
            pesState: {
              ...state.pesState,
              scanning: false,
              results: result.points,
              equilibrium: result.equilibrium,
              error: null,
              scanRequestId: null,
              computeTimeMs: result.compute_time_ms,
              totalIterations: result.total_iterations,
            },
          }),
          false,
          'completePesScan'
        );
      },

      failPesScan: (error: string) => {
        set(
          (state) => ({
            pesState: {
              ...state.pesState,
              scanning: false,
              error,
              scanRequestId: null,
            },
          }),
          false,
          'failPesScan'
        );
      },

      resetPes: () => {
        set(
          (state) => ({
            pesState: {
              ...DEFAULT_PES_STATE,
              rMin: state.pesState.rMin,
              rMax: state.pesState.rMax,
              nPoints: state.pesState.nPoints,
            },
          }),
          false,
          'resetPes'
        );
      },

      // === Phase 4: PES Coordinate Tracking Actions (US-083) ===

      completePesInternalScan: (result: PesScanInternalResult) => {
        set(
          (state) => ({
            pesState: {
              ...state.pesState,
              pesInternalResult: result,
            },
          }),
          false,
          'completePesInternalScan'
        );
      },

      setTrackedCoordinates: (keys: string[]) => {
        set(
          (state) => ({
            pesState: {
              ...state.pesState,
              trackedCoordinates: keys,
            },
          }),
          false,
          'setTrackedCoordinates'
        );
      },

      setPesCoordinateConfig: (config) => {
        set(
          (state) => ({
            pesState: { ...state.pesState, pesCoordinateConfig: config },
          }),
          false,
          'setPesCoordinateConfig'
        );
      },

      // === Phase 2: Basis Set Comparison Actions (US-045) ===

      setBasisCompareSelected: (bases: string[]) => {
        set(
          (state) => ({
            basisCompareState: { ...state.basisCompareState, selectedBases: bases },
          }),
          false,
          'setBasisCompareSelected'
        );
      },

      startBasisCompare: () => {
        set(
          (state) => ({
            basisCompareState: {
              ...state.basisCompareState,
              comparing: true,
              results: [],
              progress: null,
              error: null,
            },
          }),
          false,
          'startBasisCompare'
        );
      },

      updateBasisCompareProgress: (progress) => {
        set(
          (state) => ({
            basisCompareState: { ...state.basisCompareState, progress },
          }),
          false,
          'updateBasisCompareProgress'
        );
      },

      addBasisCompareResult: (result) => {
        set(
          (state) => ({
            basisCompareState: {
              ...state.basisCompareState,
              results: [...state.basisCompareState.results, result],
            },
          }),
          false,
          'addBasisCompareResult'
        );
      },

      completeBasisCompare: () => {
        set(
          (state) => ({
            basisCompareState: {
              ...state.basisCompareState,
              comparing: false,
              progress: null,
            },
          }),
          false,
          'completeBasisCompare'
        );
      },

      failBasisCompare: (error: string) => {
        set(
          (state) => ({
            basisCompareState: {
              ...state.basisCompareState,
              comparing: false,
              progress: null,
              error,
            },
          }),
          false,
          'failBasisCompare'
        );
      },

      resetBasisCompare: () => {
        set(
          (state) => ({
            basisCompareState: {
              ...DEFAULT_BASIS_COMPARE_STATE,
              selectedBases: state.basisCompareState.selectedBases,
            },
          }),
          false,
          'resetBasisCompare'
        );
      },

      // === Phase 3: DFT Method Selection Actions (US-070) ===

      setMethod: (method: DftMethod) => {
        set(
          (state) => ({
            method,
            // Clear results -- old data is from a different method
            compute: { status: 'idle' },
            history: [],
            runningRequestId: null,
            ksResult: null,
            // Clear DFT comparison state
            dftCompareState: { ...DEFAULT_DFT_COMPARE_STATE },
            // Clear orbital/density selection since results are invalidated
            viewerState: {
              ...state.viewerState,
              selectedOrbital: null,
              showDensity: false,
            },
          }),
          false,
          'setMethod'
        );
      },

      startDftCompare: () => {
        set(
          {
            dftCompareState: {
              comparing: true,
              phase: 'running_rhf',
              rhfResult: null,
              dftResult: null,
              error: null,
            },
          },
          false,
          'startDftCompare'
        );
      },

      setDftCompareRhfResult: (result) => {
        set(
          (state) => ({
            dftCompareState: {
              ...state.dftCompareState,
              phase: 'running_dft',
              rhfResult: result,
            },
          }),
          false,
          'setDftCompareRhfResult'
        );
      },

      setDftCompareDftResult: (result) => {
        set(
          (state) => ({
            dftCompareState: {
              ...state.dftCompareState,
              dftResult: result,
            },
          }),
          false,
          'setDftCompareDftResult'
        );
      },

      completeDftCompare: () => {
        set(
          (state) => ({
            dftCompareState: {
              ...state.dftCompareState,
              comparing: false,
              phase: 'complete',
            },
          }),
          false,
          'completeDftCompare'
        );
      },

      failDftCompare: (error: string) => {
        set(
          (state) => ({
            dftCompareState: {
              ...state.dftCompareState,
              comparing: false,
              phase: 'error',
              error,
            },
          }),
          false,
          'failDftCompare'
        );
      },

      resetDftCompare: () => {
        set(
          { dftCompareState: { ...DEFAULT_DFT_COMPARE_STATE } },
          false,
          'resetDftCompare'
        );
      },

      setKsResult: (result) => {
        set({ ksResult: result }, false, 'setKsResult');
      },

      // === Phase 3: Geometry Optimization Actions (US-075) ===

      startOptimization: (requestId: string) => {
        set(
          {
            optimizationState: {
              ...DEFAULT_OPTIMIZATION_STATE,
              running: true,
              requestId,
            },
          },
          false,
          'startOptimization'
        );
      },

      addOptimizationStep: (step: OptimizationStepResult) => {
        set(
          (state) => ({
            optimizationState: {
              ...state.optimizationState,
              progress: [...state.optimizationState.progress, step],
              // Update trajectory step to latest during optimization
              trajectoryStep: step.step,
            },
          }),
          false,
          'addOptimizationStep'
        );
      },

      completeOptimization: (result: OptimizationResult) => {
        set(
          (state) => ({
            optimizationState: {
              ...state.optimizationState,
              running: false,
              result,
              requestId: null,
              trajectoryStep: result.totalSteps,
            },
          }),
          false,
          'completeOptimization'
        );
      },

      failOptimization: (error: string) => {
        set(
          (state) => ({
            optimizationState: {
              ...state.optimizationState,
              running: false,
              error,
              requestId: null,
            },
          }),
          false,
          'failOptimization'
        );
      },

      cancelOptimization: () => {
        set(
          (state) => ({
            optimizationState: {
              ...state.optimizationState,
              running: false,
              requestId: null,
            },
          }),
          false,
          'cancelOptimization'
        );
      },

      setTrajectoryStep: (step: number) => {
        set(
          (state) => ({
            optimizationState: { ...state.optimizationState, trajectoryStep: step },
          }),
          false,
          'setTrajectoryStep'
        );
      },

      setShowGhostOverlay: (show: boolean) => {
        set(
          (state) => ({
            optimizationState: { ...state.optimizationState, showGhostOverlay: show },
          }),
          false,
          'setShowGhostOverlay'
        );
      },

      applyOptimizedGeometry: (
        atoms: [number, number, number, number][],
        basisName: string
      ) => {
        // Convert [Z, x, y, z] to GeometryInput format
        const Z_TO_SYMBOL: Record<number, string> = {
          1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B',
          6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
        };
        const geometry = {
          atoms: atoms.map((a) => ({
            symbol: Z_TO_SYMBOL[a[0]] ?? 'X',
            xyz: [a[1], a[2], a[3]] as [number, number, number],
          })),
          units: 'bohr' as const,
        };

        set(
          (state) => ({
            inputMode: {
              mode: 'custom' as const,
              geometry,
              basisSet: basisName as BasisSetName,
            },
            customGeometryText: geometry.atoms
              .map((a) => `${a.symbol}  ${a.xyz[0].toFixed(6)}  ${a.xyz[1].toFixed(6)}  ${a.xyz[2].toFixed(6)}`)
              .join('\n'),
            customUnits: 'bohr' as const,
            // Reset SCF results since geometry changed
            compute: { status: 'idle' },
            history: [],
            runningRequestId: null,
            ksResult: null,
            // Reset integral compute (new geometry needs new integrals for RHF)
            integralCompute: { status: 'idle' },
            // Clear orbital/density selection
            viewerState: {
              ...state.viewerState,
              selectedOrbital: null,
              showDensity: false,
            },
          }),
          false,
          'applyOptimizedGeometry'
        );
      },

      resetOptimization: () => {
        set(
          { optimizationState: { ...DEFAULT_OPTIMIZATION_STATE } },
          false,
          'resetOptimization'
        );
      },

      // Phase 3: Population analysis (US-076)
      setPopulationResult: (result) => {
        set({ populationResult: result }, false, 'setPopulationResult');
      },

      setShowChargeLabels: (show) => {
        set({ showChargeLabels: show }, false, 'setShowChargeLabels');
      },

      setShowChargeOverlay: (show) => {
        set({ showChargeOverlay: show }, false, 'setShowChargeOverlay');
      },

      // Phase 5: Frequency analysis actions (US-103)

      setFrequencyResult: (result) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, result } }),
          false,
          'setFrequencyResult'
        );
      },

      setFrequencyIsComputing: (isComputing) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, isComputing } }),
          false,
          'setFrequencyIsComputing'
        );
      },

      setFrequencyProgress: (progress) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, progress } }),
          false,
          'setFrequencyProgress'
        );
      },

      setFrequencyError: (error) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, error } }),
          false,
          'setFrequencyError'
        );
      },

      setFrequencySelectedMode: (selectedMode) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, selectedMode } }),
          false,
          'setFrequencySelectedMode'
        );
      },

      setFrequencyTemperature: (temperatureK) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, temperatureK } }),
          false,
          'setFrequencyTemperature'
        );
      },

      setFrequencyPressure: (pressurePa) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, pressurePa } }),
          false,
          'setFrequencyPressure'
        );
      },

      setFrequencyDisplayThermo: (displayThermo) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, displayThermo } }),
          false,
          'setFrequencyDisplayThermo'
        );
      },

      setFrequencySpectrumTab: (spectrumTab) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, spectrumTab } }),
          false,
          'setFrequencySpectrumTab'
        );
      },

      setFrequencyBroadeningKind: (broadeningKind) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, broadeningKind } }),
          false,
          'setFrequencyBroadeningKind'
        );
      },

      setFrequencyFwhmCm1: (fwhmCm1) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, fwhmCm1 } }),
          false,
          'setFrequencyFwhmCm1'
        );
      },

      setFrequencyAmplitudeBohr: (amplitudeBohr) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, amplitudeBohr } }),
          false,
          'setFrequencyAmplitudeBohr'
        );
      },

      setFrequencyAnimationSpeed: (animationSpeed) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, animationSpeed } }),
          false,
          'setFrequencyAnimationSpeed'
        );
      },

      setFrequencyIsAnimating: (isAnimating) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, isAnimating } }),
          false,
          'setFrequencyIsAnimating'
        );
      },

      setFrequencyShowDisplacementArrows: (showDisplacementArrows) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, showDisplacementArrows } }),
          false,
          'setFrequencyShowDisplacementArrows'
        );
      },

      setFrequencyUnitsMode: (unitsMode) => {
        set(
          (state) => ({ frequencyState: { ...state.frequencyState, unitsMode } }),
          false,
          'setFrequencyUnitsMode'
        );
      },

      resetFrequencyState: () => {
        set(
          { frequencyState: { ...DEFAULT_FREQUENCY_STATE } },
          false,
          'resetFrequencyState'
        );
      },

      reset: () => {
        set({ ...DEFAULT_STATE, urlInitialized: true }, false, 'reset');
      },

      clearHistory: () => {
        set({ history: [], integralPercent: 0, integralStep: '', compute: { status: 'idle' } }, false, 'clearHistory');
      },
    }),
    { name: 'scfStore' }
  )
);

/**
 * Get current SCF parameters in ScfParams format.
 *
 * Useful for URL synchronization and deep link generation.
 * Includes custom geometry information when in custom mode.
 *
 * @returns Current parameters as ScfParams
 */
export function getScfParams(): ScfParams {
  const { systemId, selectedMolecule, selectedBasis, convergenceProfile, maxIterations, useDiis, damp, method, inputMode, pesState, basisCompareState, viewerState, frequencyState } = useScfStore.getState();

  const baseParams: ScfParams = {
    system_id: systemId,
    conv: convergenceProfile,
    max_iter: maxIterations,
    diis: useDiis,
    damp: damp > 0 ? damp : undefined, // Only include if non-zero
    // Only include method when non-default (US-070)
    method: method !== 'rhf' ? method : undefined,
    // Include molecule + basis when non-default for deep link
    molecule: selectedMolecule !== DEFAULT_STATE.selectedMolecule ? selectedMolecule : undefined,
    basis: selectedBasis !== DEFAULT_STATE.selectedBasis ? selectedBasis as string : undefined,
  };

  // Include PES scan config if non-default (US-085: includes internal coordinate fields)
  const coordConfig = pesState.pesCoordinateConfig;
  if (coordConfig) {
    // Internal coordinate config available: check if it differs from defaults
    const isDefault =
      coordConfig.coordinateType === 'bond' &&
      coordConfig.atomIndices.length === 2 &&
      coordConfig.atomIndices[0] === 0 &&
      coordConfig.atomIndices[1] === 1 &&
      coordConfig.scanMode === 'rigid' &&
      coordConfig.valueMin === DEFAULT_PES_STATE.rMin &&
      coordConfig.valueMax === DEFAULT_PES_STATE.rMax &&
      coordConfig.nPoints === DEFAULT_PES_STATE.nPoints;

    if (!isDefault) {
      baseParams.pes = {
        r_min: coordConfig.valueMin,
        r_max: coordConfig.valueMax,
        n_points: coordConfig.nPoints,
        coordinate_type: coordConfig.coordinateType,
        atom_indices: coordConfig.atomIndices,
        scan_mode: coordConfig.scanMode,
      };
    }
  } else {
    // Legacy path: use rMin/rMax/nPoints from pesState
    const hasNonDefault =
      pesState.rMin !== DEFAULT_PES_STATE.rMin ||
      pesState.rMax !== DEFAULT_PES_STATE.rMax ||
      pesState.nPoints !== DEFAULT_PES_STATE.nPoints;

    if (hasNonDefault) {
      baseParams.pes = {
        r_min: pesState.rMin,
        r_max: pesState.rMax,
        n_points: pesState.nPoints,
      };
    }
  }

  // Include basis comparison selection if non-default (US-045)
  const defaultBases = DEFAULT_BASIS_COMPARE_STATE.selectedBases;
  const currentBases = basisCompareState.selectedBases;
  const basesChanged =
    currentBases.length !== defaultBases.length ||
    currentBases.some((b, i) => b !== defaultBases[i]);
  if (basesChanged) {
    baseParams.comparison_bases = currentBases;
  }

  // Include density visualization params if non-default (US-062, US-063)
  if (viewerState.densityMode !== DEFAULT_DENSITY_STATE.densityMode) {
    baseParams.density_mode = viewerState.densityMode;
  }
  if (viewerState.densityIsovalue !== DEFAULT_DENSITY_STATE.densityIsovalue) {
    baseParams.density_isovalue = viewerState.densityIsovalue;
  }
  if (viewerState.diffIsovalue !== DEFAULT_DENSITY_STATE.diffIsovalue) {
    baseParams.diff_isovalue = viewerState.diffIsovalue;
  }
  if (viewerState.planeAxis !== DEFAULT_DENSITY_STATE.planeAxis) {
    baseParams.plane_axis = viewerState.planeAxis;
  }
  if (viewerState.planePosition !== DEFAULT_DENSITY_STATE.planePosition) {
    baseParams.plane_position = viewerState.planePosition;
  }
  if (viewerState.colorScale !== DEFAULT_DENSITY_STATE.colorScale) {
    baseParams.color_scale = viewerState.colorScale;
  }

  // Include frequency state params if non-default (US-103)
  if (frequencyState.selectedMode !== null) {
    baseParams.freq_mode = frequencyState.selectedMode;
  }
  if (frequencyState.temperatureK !== DEFAULT_FREQUENCY_STATE.temperatureK) {
    baseParams.freq_temp = frequencyState.temperatureK;
  }
  if (frequencyState.pressurePa !== DEFAULT_FREQUENCY_STATE.pressurePa) {
    baseParams.freq_pres = frequencyState.pressurePa;
  }
  if (frequencyState.broadeningKind !== DEFAULT_FREQUENCY_STATE.broadeningKind) {
    baseParams.freq_broad = frequencyState.broadeningKind === 'gaussian' ? 'g' : 'l';
  }
  if (frequencyState.fwhmCm1 !== DEFAULT_FREQUENCY_STATE.fwhmCm1) {
    baseParams.freq_fwhm = frequencyState.fwhmCm1;
  }
  if (frequencyState.spectrumTab !== DEFAULT_FREQUENCY_STATE.spectrumTab) {
    baseParams.freq_tab = frequencyState.spectrumTab;
  }
  if (frequencyState.showDisplacementArrows !== DEFAULT_FREQUENCY_STATE.showDisplacementArrows) {
    baseParams.freq_arrows = frequencyState.showDisplacementArrows;
  }
  if (frequencyState.unitsMode !== DEFAULT_FREQUENCY_STATE.unitsMode) {
    baseParams.freq_units = frequencyState.unitsMode === 'hartree' ? 'ha' : frequencyState.unitsMode === 'kj_mol' ? 'kj' : 'kc';
  }

  // Include custom geometry data when in custom mode
  if (inputMode.mode === 'custom') {
    return {
      ...baseParams,
      input_mode: 'custom',
      custom_geometry: {
        atoms: inputMode.geometry.atoms.map((a) => ({
          symbol: a.symbol,
          xyz: a.xyz,
        })),
        units: inputMode.geometry.units,
      },
      custom_basis: inputMode.basisSet,
    };
  }

  return baseParams;
}

/**
 * Get system preset by ID.
 *
 * @param id - System preset ID
 * @returns SystemPreset or undefined if not found
 */
export function getSystemPreset(id: string): SystemPreset | undefined {
  return SYSTEM_PRESETS.find((p) => p.id === id);
}
