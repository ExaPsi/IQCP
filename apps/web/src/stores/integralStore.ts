/**
 * Integral Inspector State Store
 *
 * Zustand store managing the state for Module E (Integral Inspector).
 * Handles system selection, matrix tab navigation, cell selection,
 * and computation state for one-electron integral matrices.
 *
 * The store does not use `persist` middleware because integral state
 * is transient -- deep links capture the configuration, not the results.
 *
 * @module stores/integralStore
 * @see US-055 Integral Matrix Heatmap
 */

import { create } from 'zustand';
import { devtools } from 'zustand/middleware';
import type { GeometryInput, BasisSetName } from '../worker/protocol';

// ============================================================================
// Types
// ============================================================================

/** Which matrix is currently being viewed */
export type MatrixTab = 'S' | 'T' | 'V' | 'Hcore';

/** Selected cell in the heatmap */
export interface SelectedCell {
  row: number;
  col: number;
}

/**
 * Result from integral_matrices worker request.
 *
 * This type mirrors the IntegralMatricesResult that will be returned
 * by the WASM module. The field names use camelCase matching Rust's
 * serde(rename_all = "camelCase") convention.
 */
export interface IntegralMatricesResult {
  /** Number of basis functions */
  nbf: number;
  /** Basis function labels: ["O 1s", "O 2s", ...] */
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
 * Async state for integral computation.
 * Uses a discriminated union for type-safe status handling.
 */
export type IntegralComputeState =
  | { status: 'idle' }
  | { status: 'computing' }
  | { status: 'success'; data: IntegralMatricesResult }
  | { status: 'error'; error: string };

/**
 * Async state for ERI (electron repulsion integral) background caching.
 *
 * ERI computation runs in the background after one-electron integrals
 * succeed. The cached ERI data is needed by US-058 (Fock Build Tracing)
 * and US-059 (ERI Browser).
 *
 * Uses 8-fold symmetry compressed storage: the ERI tensor (mu nu | la si)
 * is stored as a flat array with pair_ij*(pair_ij+1)/2 + pair_kl indexing.
 *
 * @see US-057 System Selection (AC3)
 */
export type EriCacheState =
  | { status: 'idle' }
  | { status: 'computing' }
  | { status: 'success'; eriCompressed: number[]; eriIndexing: string; eriCount: number; computeTimeMs: number }
  | { status: 'error'; error: string };

/** Input mode for molecule specification */
export type IntegralInputMode = 'preset' | 'custom';

// ============================================================================
// Fock Build Tracing Types (US-058)
// ============================================================================

/**
 * Result from fock_decomposition worker request.
 *
 * Contains the Fock matrix decomposed into its physical components:
 * F = H^core + G, where G = J - 0.5*K.
 *
 * All matrices are flat arrays, nbf x nbf, row-major order.
 * Field names use camelCase matching Rust's serde(rename_all = "camelCase").
 *
 * @see US-058 Fock Build Tracing
 * @see Szabo & Ostlund (1996), Eq. 3.154
 */
export interface FockDecompositionResult {
  /** Core Hamiltonian H^core = T + V (row-major, nbf x nbf) */
  hCore: number[];
  /** Coulomb matrix J: J_{mn} = sum_{ls} P_{ls} (mn|ls) (row-major, nbf x nbf) */
  jMatrix: number[];
  /** Exchange matrix K: K_{mn} = sum_{ls} P_{ls} (ml|ns) (row-major, nbf x nbf) */
  kMatrix: number[];
  /** Two-electron contribution G = J - 0.5*K (row-major, nbf x nbf) */
  gMatrix: number[];
  /** Full Fock matrix F = H^core + G (row-major, nbf x nbf) */
  fMatrix: number[];
  /** Density matrix P (echoed, row-major, nbf x nbf) */
  density: number[];
  /** Number of basis functions */
  nbf: number;
  /** Basis function labels */
  labels: string[];
}

/**
 * Async state for Fock decomposition computation.
 * Uses a discriminated union for type-safe status handling.
 *
 * @see US-058 Fock Build Tracing
 */
export type FockDecompositionState =
  | { status: 'idle' }
  | { status: 'computing' }
  | { status: 'success'; data: FockDecompositionResult }
  | { status: 'error'; error: string };

/** Fock build step: 1 = H^core, 2 = G(P) with J/K, 3 = F = H^core + G */
export type FockBuildStep = 1 | 2 | 3;

/** Sub-view within step 2 (G decomposition) */
export type FockStep2View = 'G' | 'J' | 'K';

// ============================================================================
// ERI Browser Types (US-059)
// ============================================================================

/**
 * Import EriDetailResult from the hook to avoid circular dependency.
 * The type is re-declared inline here for store independence.
 */
export interface EriDetailResultStore {
  contractedValue: number;
  indices: [number, number, number, number];
  labels: [string, string, string, string];
  method:
    | { type: 'boysFunction'; tParameter: number }
    | { type: 'rysQuadrature'; nroots: number; tParameter: number; roots: number[]; weights: number[] };
  contributions: Array<{
    primIndices: [number, number, number, number];
    exponents: [number, number, number, number];
    coefficients: [number, number, number, number];
    normCoefficients: [number, number, number, number];
    primitiveValue: number;
    weightedContribution: number;
    tParameter: number;
  }>;
  nPrimitives: [number, number, number, number];
  totalAngularMomentum: number;
  nroots: number;
}

/**
 * Async state for ERI detail computation.
 * Uses a discriminated union for type-safe status handling.
 *
 * @see US-059 ERI Browser
 */
export type EriDetailState =
  | { status: 'idle' }
  | { status: 'computing' }
  | { status: 'success'; data: EriDetailResultStore }
  | { status: 'error'; error: string };

// ============================================================================
// State Interface
// ============================================================================

/**
 * Integral Inspector state.
 */
export interface IntegralState {
  /** Selected preset system ID (e.g., "h2_sto3g_r1.4") */
  systemId: string;
  /** Input mode: preset or custom geometry */
  inputMode: IntegralInputMode;
  /** Custom geometry text (XYZ format) */
  customGeometryText: string;
  /** Parsed custom geometry (null if invalid or not custom mode) */
  customGeometry: GeometryInput | null;
  /** Basis set for custom geometry */
  customBasis: BasisSetName;
  /** Whether to use spherical harmonics (default true) */
  useSpherical: boolean;

  /** Active matrix tab */
  activeTab: MatrixTab;
  /** Selected cell in the heatmap (null if none) */
  selectedCell: SelectedCell | null;

  /** Computation state for one-electron integrals (S, T, V, Hcore) */
  compute: IntegralComputeState;

  /**
   * Background ERI cache state.
   * Computed after one-electron integrals succeed, cached for US-058/US-059.
   * @see US-057 AC3
   */
  eriCache: EriCacheState;

  /** Whether URL state has been initialized */
  urlInitialized: boolean;

  /** Whether all primitive pairs are shown in the breakdown table (US-056) */
  showAllPrimitives: boolean;

  // Fock Build Tracing (US-058)

  /** Current Fock build step (1=H^core, 2=G with J/K, 3=F=H^core+G) */
  fockBuildStep: FockBuildStep;
  /** Sub-view within step 2 (G, J, or K matrix) */
  fockStep2View: FockStep2View;
  /** Selected matrix element for Fock drill-down [row, col] */
  fockSelectedElement: [number, number] | null;
  /** Fock decomposition computation state */
  fockDecomposition: FockDecompositionState;

  // ERI Browser (US-059)

  /** Selected ERI quartet [i, j, k, l] (null if none) */
  eriSelectedQuartet: [number, number, number, number] | null;
  /** ERI detail computation state */
  eriDetail: EriDetailState;
  /** Whether to show all primitive quartets in the ERI breakdown table */
  eriShowAllPrimitives: boolean;
}

/**
 * Integral Inspector actions.
 */
export interface IntegralActions {
  /** Set the preset system ID (clears compute state) */
  setSystemId: (id: string) => void;
  /** Set the input mode */
  setInputMode: (mode: IntegralInputMode) => void;
  /** Set custom geometry text */
  setCustomGeometryText: (text: string) => void;
  /** Set parsed custom geometry */
  setCustomGeometry: (geom: GeometryInput | null) => void;
  /** Set basis set for custom geometry */
  setCustomBasis: (basis: BasisSetName) => void;
  /** Set spherical harmonics toggle */
  setUseSpherical: (use: boolean) => void;
  /** Switch matrix tab (preserves cell selection) */
  setActiveTab: (tab: MatrixTab) => void;
  /** Select a cell (or null to clear) */
  selectCell: (cell: SelectedCell | null) => void;
  /** Update computation state */
  setCompute: (state: IntegralComputeState) => void;
  /** Update ERI cache state */
  setEriCache: (state: EriCacheState) => void;
  /** Mark URL state as initialized */
  setUrlInitialized: () => void;
  /** Toggle show all primitives in breakdown table (US-056) */
  toggleShowAllPrimitives: () => void;

  // Fock Build Tracing (US-058)

  /** Set the current Fock build step */
  setFockBuildStep: (step: FockBuildStep) => void;
  /** Set the sub-view within step 2 */
  setFockStep2View: (view: FockStep2View) => void;
  /** Select a Fock matrix element for drill-down (or null to clear) */
  setFockSelectedElement: (element: [number, number] | null) => void;
  /** Update the Fock decomposition computation state */
  setFockDecomposition: (state: FockDecompositionState) => void;

  // ERI Browser (US-059)

  /** Set the selected ERI quartet (or null to clear) */
  setEriSelectedQuartet: (quartet: [number, number, number, number] | null) => void;
  /** Update the ERI detail computation state */
  setEriDetail: (state: EriDetailState) => void;
  /** Toggle show all primitive quartets in the ERI breakdown table */
  toggleEriShowAllPrimitives: () => void;

  /** Reset to default state */
  reset: () => void;
}

// ============================================================================
// Default State
// ============================================================================

/**
 * Default H2 geometry in XYZ format (Bohr units).
 */
export const DEFAULT_GEOMETRY_TEXT = `H  0.000000  0.000000  0.000000
H  0.000000  0.000000  1.398400`;

const DEFAULT_STATE: IntegralState = {
  systemId: 'h2_sto3g_r1.4',
  inputMode: 'preset',
  customGeometryText: DEFAULT_GEOMETRY_TEXT,
  customGeometry: null,
  customBasis: 'sto-3g',
  useSpherical: true,

  activeTab: 'S',
  selectedCell: null,

  compute: { status: 'idle' },
  eriCache: { status: 'idle' },
  urlInitialized: false,
  showAllPrimitives: false,

  // Fock Build Tracing (US-058)
  fockBuildStep: 1,
  fockStep2View: 'G',
  fockSelectedElement: null,
  fockDecomposition: { status: 'idle' },

  // ERI Browser (US-059)
  eriSelectedQuartet: null,
  eriDetail: { status: 'idle' },
  eriShowAllPrimitives: false,
};

// ============================================================================
// Store
// ============================================================================

/**
 * Integral Inspector state store.
 *
 * @example
 * ```typescript
 * import { useIntegralStore } from '@/stores/integralStore';
 *
 * function IntegralControls() {
 *   const systemId = useIntegralStore(s => s.systemId);
 *   const setSystemId = useIntegralStore(s => s.setSystemId);
 *   // ...
 * }
 * ```
 */
export const useIntegralStore = create<IntegralState & IntegralActions>()(
  devtools(
    (set) => ({
      ...DEFAULT_STATE,

      setSystemId: (id: string) => {
        set(
          {
            systemId: id,
            compute: { status: 'idle' },
            eriCache: { status: 'idle' },
            selectedCell: null,
            fockDecomposition: { status: 'idle' },
            fockBuildStep: 1,
            fockStep2View: 'G',
            fockSelectedElement: null,
            eriSelectedQuartet: null,
            eriDetail: { status: 'idle' },
            eriShowAllPrimitives: false,
          },
          false,
          'setSystemId',
        );
      },

      setInputMode: (mode: IntegralInputMode) => {
        set(
          {
            inputMode: mode,
            compute: { status: 'idle' },
            eriCache: { status: 'idle' },
            selectedCell: null,
            fockDecomposition: { status: 'idle' },
            fockSelectedElement: null,
            eriSelectedQuartet: null,
            eriDetail: { status: 'idle' },
            eriShowAllPrimitives: false,
          },
          false,
          'setInputMode',
        );
      },

      setCustomGeometryText: (text: string) => {
        set({ customGeometryText: text }, false, 'setCustomGeometryText');
      },

      setCustomGeometry: (geom: GeometryInput | null) => {
        set(
          {
            customGeometry: geom,
            compute: { status: 'idle' },
            eriCache: { status: 'idle' },
            selectedCell: null,
            fockDecomposition: { status: 'idle' },
            fockSelectedElement: null,
            eriSelectedQuartet: null,
            eriDetail: { status: 'idle' },
            eriShowAllPrimitives: false,
          },
          false,
          'setCustomGeometry',
        );
      },

      setCustomBasis: (basis: BasisSetName) => {
        set(
          {
            customBasis: basis,
            compute: { status: 'idle' },
            eriCache: { status: 'idle' },
            selectedCell: null,
            fockDecomposition: { status: 'idle' },
            fockSelectedElement: null,
            eriSelectedQuartet: null,
            eriDetail: { status: 'idle' },
            eriShowAllPrimitives: false,
          },
          false,
          'setCustomBasis',
        );
      },

      setUseSpherical: (use: boolean) => {
        set(
          {
            useSpherical: use,
            compute: { status: 'idle' },
            eriCache: { status: 'idle' },
            selectedCell: null,
          },
          false,
          'setUseSpherical',
        );
      },

      setActiveTab: (tab: MatrixTab) => {
        set({ activeTab: tab, showAllPrimitives: false }, false, 'setActiveTab');
      },

      selectCell: (cell: SelectedCell | null) => {
        set({ selectedCell: cell, showAllPrimitives: false }, false, 'selectCell');
      },

      setCompute: (state: IntegralComputeState) => {
        set({ compute: state }, false, 'setCompute');
      },

      setEriCache: (state: EriCacheState) => {
        set({ eriCache: state }, false, 'setEriCache');
      },

      setUrlInitialized: () => {
        set({ urlInitialized: true }, false, 'setUrlInitialized');
      },

      toggleShowAllPrimitives: () => {
        set(
          (state) => ({ showAllPrimitives: !state.showAllPrimitives }),
          false,
          'toggleShowAllPrimitives',
        );
      },

      // Fock Build Tracing (US-058)

      setFockBuildStep: (step: FockBuildStep) => {
        set({ fockBuildStep: step, fockSelectedElement: null }, false, 'setFockBuildStep');
      },

      setFockStep2View: (view: FockStep2View) => {
        set({ fockStep2View: view }, false, 'setFockStep2View');
      },

      setFockSelectedElement: (element: [number, number] | null) => {
        set({ fockSelectedElement: element }, false, 'setFockSelectedElement');
      },

      setFockDecomposition: (state: FockDecompositionState) => {
        set({ fockDecomposition: state }, false, 'setFockDecomposition');
      },

      // ERI Browser (US-059)

      setEriSelectedQuartet: (quartet: [number, number, number, number] | null) => {
        set(
          {
            eriSelectedQuartet: quartet,
            eriDetail: { status: 'idle' },
            eriShowAllPrimitives: false,
          },
          false,
          'setEriSelectedQuartet',
        );
      },

      setEriDetail: (state: EriDetailState) => {
        set({ eriDetail: state }, false, 'setEriDetail');
      },

      toggleEriShowAllPrimitives: () => {
        set(
          (state) => ({ eriShowAllPrimitives: !state.eriShowAllPrimitives }),
          false,
          'toggleEriShowAllPrimitives',
        );
      },

      reset: () => {
        set({ ...DEFAULT_STATE }, false, 'reset');
      },
    }),
    { name: 'integralStore' },
  ),
);
