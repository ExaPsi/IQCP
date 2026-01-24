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
} from '../worker/protocol';
import type { ScfParams } from '../types/run-state';

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
 * SCF store state interface.
 */
export interface ScfState {
  // === Parameters (user inputs) ===
  /** Selected molecular system ID */
  systemId: string;
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
  /** Current running request ID (for cancellation) */
  runningRequestId: string | null;
  /** Running integral compute request ID (for cancellation) */
  integralComputeRequestId: string | null;

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
  setConvergenceProfile: (profile: ConvergenceProfile) => void;
  setMaxIterations: (max: number) => void;
  setUseDiis: (useDiis: boolean) => void;
  setDamp: (damp: number) => void;
  setUseSpherical: (useSpherical: boolean) => void;
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
  setRunResult: (result: ScfRunResult) => void;
  setRunError: (error: string) => void;
  setRunCancelled: () => void;

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
  convergenceProfile: 'medium',
  maxIterations: 50,
  useDiis: true,
  damp: 0.0,
  useSpherical: false, // Cartesian for backward compatibility
  mode: 'explain',
  // Custom geometry defaults
  inputMode: { mode: 'preset' },
  integralCompute: { status: 'idle' },
  customGeometryText: '',
  customUnits: 'angstrom',
  // Computation state
  compute: { status: 'idle' },
  history: [],
  runningRequestId: null,
  integralComputeRequestId: null,
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
    (set) => ({
      // Initial state
      ...DEFAULT_STATE,

      // Actions
      setSystemId: (systemId: string) => {
        set({ systemId }, false, 'setSystemId');
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

      setMode: (mode: DisplayMode) => {
        set({ mode }, false, 'setMode');
      },

      initializeFromURL: (params: ScfParams) => {
        // Handle custom mode initialization from URL
        const dampValue = params.damp !== undefined ? Math.max(0.0, Math.min(0.9, params.damp)) : 0.0;
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
              convergenceProfile: params.conv,
              maxIterations: Math.max(1, Math.min(200, params.max_iter)),
              useDiis: params.diis,
              damp: dampValue,
              inputMode: { mode: 'custom', geometry, basisSet: params.custom_basis as BasisSetName },
              customUnits: params.custom_geometry.units,
              urlInitialized: true,
            },
            false,
            'initializeFromURL'
          );
        } else {
          // Preset mode (default)
          set(
            {
              systemId: params.system_id,
              convergenceProfile: params.conv,
              maxIterations: Math.max(1, Math.min(200, params.max_iter)),
              useDiis: params.diis,
              damp: dampValue,
              inputMode: { mode: 'preset' },
              urlInitialized: true,
            },
            false,
            'initializeFromURL'
          );
        }
      },

      // === Input Mode Switching ===

      switchToPresetMode: () => {
        set(
          {
            inputMode: { mode: 'preset' },
            integralCompute: { status: 'idle' },
          },
          false,
          'switchToPresetMode'
        );
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
          {
            compute: { status: 'running' },
            history: [],
            runningRequestId: requestId,
          },
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

      setRunResult: (result: ScfRunResult) => {
        set(
          {
            compute: { status: 'success', result },
            runningRequestId: null,
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
          },
          false,
          'setRunCancelled'
        );
      },

      reset: () => {
        set({ ...DEFAULT_STATE, urlInitialized: true }, false, 'reset');
      },

      clearHistory: () => {
        set({ history: [], compute: { status: 'idle' } }, false, 'clearHistory');
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
  const { systemId, convergenceProfile, maxIterations, useDiis, damp, inputMode } = useScfStore.getState();

  const baseParams: ScfParams = {
    system_id: systemId,
    conv: convergenceProfile,
    max_iter: maxIterations,
    diis: useDiis,
    damp: damp > 0 ? damp : undefined, // Only include if non-zero
  };

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
