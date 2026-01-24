/**
 * ScfControlsPanel - Parameter input controls for SCF computation.
 *
 * Provides inputs for:
 * - Input mode toggle (Preset/Custom)
 * - System preset selector (dropdown) - in preset mode
 * - Custom geometry editor with validation - in custom mode
 * - Basis set selector - in custom mode
 * - Convergence profile (button group)
 * - DIIS toggle (checkbox)
 * - Max iterations (number input)
 * - Compute Integrals button - in custom mode
 * - Run/Cancel buttons
 *
 * Connected to the scfStore for state management.
 *
 * @module components/scf/ScfControlsPanel
 */

import { useCallback, useMemo } from 'react';
import {
  useScfStore,
  SYSTEM_PRESETS,
  CONVERGENCE_THRESHOLDS,
  getSystemPreset,
} from '../../stores/scfStore';
import { useScf } from '../../hooks/useScf';
import { useIntegralCompute } from '../../hooks/useIntegralCompute';
import { parseAndValidate } from '../../lib/geometryValidation';
import { Math } from '../common/Math';
import type { ConvergenceProfile, IntegralProgress, BasisSetName } from '../../worker/protocol';
import type { ExampleGeometry } from '../../lib/exampleGeometries';
import { geometryToXyzText } from '../../lib/geometryValidation';

// Import custom geometry components
import {
  InputModeToggle,
  GeometryEditor,
  BasisSetSelector,
  BasisTypeToggle,
  IntegralComputeButton,
  ExampleGeometrySelector,
} from './custom';

/**
 * Props for ScfControlsPanel.
 */
interface ScfControlsPanelProps {
  /** Disable all controls */
  disabled?: boolean;
}

/**
 * Convergence profile options for the button group.
 */
const CONVERGENCE_OPTIONS: { value: ConvergenceProfile; label: string }[] = [
  { value: 'loose', label: 'Loose' },
  { value: 'medium', label: 'Medium' },
  { value: 'tight', label: 'Tight' },
];

/**
 * Controls panel for SCF parameters.
 *
 * Includes:
 * - Input mode toggle (Preset/Custom)
 * - System selector dropdown (preset mode) or geometry editor (custom mode)
 * - Convergence profile button group (loose/medium/tight)
 * - DIIS toggle checkbox
 * - Max iterations number input
 * - Compute Integrals button (custom mode)
 * - Run/Cancel action buttons
 *
 * @example
 * ```tsx
 * <ScfControlsPanel />
 * ```
 */
export function ScfControlsPanel({ disabled = false }: ScfControlsPanelProps) {
  // Store state
  const systemId = useScfStore((state) => state.systemId);
  const convergenceProfile = useScfStore((state) => state.convergenceProfile);
  const maxIterations = useScfStore((state) => state.maxIterations);
  const useDiis = useScfStore((state) => state.useDiis);
  const damp = useScfStore((state) => state.damp);
  const useSpherical = useScfStore((state) => state.useSpherical);
  const compute = useScfStore((state) => state.compute);
  const inputMode = useScfStore((state) => state.inputMode);
  const integralCompute = useScfStore((state) => state.integralCompute);
  const customGeometryText = useScfStore((state) => state.customGeometryText);
  const customUnits = useScfStore((state) => state.customUnits);

  // Store actions
  const setSystemId = useScfStore((state) => state.setSystemId);
  const setConvergenceProfile = useScfStore((state) => state.setConvergenceProfile);
  const setMaxIterations = useScfStore((state) => state.setMaxIterations);
  const setUseDiis = useScfStore((state) => state.setUseDiis);
  const setDamp = useScfStore((state) => state.setDamp);
  const setUseSpherical = useScfStore((state) => state.setUseSpherical);
  const switchToPresetMode = useScfStore((state) => state.switchToPresetMode);
  const switchToCustomMode = useScfStore((state) => state.switchToCustomMode);
  const setCustomGeometryText = useScfStore((state) => state.setCustomGeometryText);
  const setCustomGeometry = useScfStore((state) => state.setCustomGeometry);
  const setCustomUnits = useScfStore((state) => state.setCustomUnits);
  const setCustomBasisSet = useScfStore((state) => state.setCustomBasisSet);

  // Hooks for computation
  const { isReady, isRunning, run, cancelRun } = useScf();
  const {
    isReady: integralReady,
    isComputing: isComputingIntegrals,
    compute: computeIntegrals,
    cancel: cancelIntegrals,
  } = useIntegralCompute();

  // Derived state
  const isPresetMode = inputMode.mode === 'preset';
  const isCustomMode = inputMode.mode === 'custom';

  // Get current system info for display (preset mode)
  const currentSystem = isPresetMode ? getSystemPreset(systemId) : null;

  // Validate custom geometry text
  const validationResult = useMemo(
    () => (isCustomMode ? parseAndValidate(customGeometryText) : null),
    [isCustomMode, customGeometryText]
  );

  const isGeometryValid = validationResult?.valid ?? false;
  const hasIntegrals = integralCompute.status === 'success';

  // In custom mode, SCF can only run after integrals are computed
  const canRunScf = isPresetMode || hasIntegrals;

  // Determine if controls should be disabled
  const controlsDisabled = disabled || isRunning || isComputingIntegrals;

  // Handle mode change
  const handleModeChange = useCallback(
    (mode: 'preset' | 'custom') => {
      if (mode === 'preset') {
        switchToPresetMode();
      } else {
        switchToCustomMode();
      }
    },
    [switchToPresetMode, switchToCustomMode]
  );

  // Handle geometry text change
  const handleGeometryTextChange = useCallback(
    (text: string) => {
      setCustomGeometryText(text);
      // Also parse and set the structured geometry
      const result = parseAndValidate(text);
      if (result.valid) {
        setCustomGeometry({
          atoms: result.atoms,
          units: customUnits,
        });
      }
    },
    [setCustomGeometryText, setCustomGeometry, customUnits]
  );

  // Handle units change
  const handleUnitsChange = useCallback(
    (units: 'bohr' | 'angstrom') => {
      setCustomUnits(units);
      // Re-parse with new units
      const result = parseAndValidate(customGeometryText);
      if (result.valid) {
        setCustomGeometry({
          atoms: result.atoms,
          units,
        });
      }
    },
    [setCustomUnits, setCustomGeometry, customGeometryText]
  );

  // Handle example geometry selection
  const handleExampleSelect = useCallback(
    (example: ExampleGeometry) => {
      // Convert geometry to XYZ text
      const text = geometryToXyzText(example.geometry);
      setCustomGeometryText(text);
      setCustomGeometry(example.geometry);
      setCustomUnits(example.geometry.units);
      setCustomBasisSet(example.recommendedBasis);
    },
    [setCustomGeometryText, setCustomGeometry, setCustomUnits, setCustomBasisSet]
  );

  // Build progress object for IntegralComputeButton
  const progress: IntegralProgress | undefined =
    integralCompute.status === 'computing'
      ? {
          module: 'integral' as const,
          phase: integralCompute.phase,
          overallPercent: integralCompute.progress,
          current: 0,
          total: 0,
          message: `Computing ${integralCompute.phase}...`,
        }
      : undefined;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
      <h2 className="text-lg font-semibold text-slate-800 mb-4">Parameters</h2>

      {/* Input Mode Toggle */}
      <div className="mb-6">
        <label className="block text-sm font-medium text-slate-700 mb-2">Input Mode</label>
        <InputModeToggle
          mode={inputMode.mode}
          onChange={handleModeChange}
          disabled={controlsDisabled}
        />
      </div>

      {/* Preset Mode: System Selector */}
      {isPresetMode && (
        <div className="mb-6">
          <label
            htmlFor="scf-system-select"
            className="block text-sm font-medium text-slate-700 mb-2"
          >
            Molecular System
          </label>
          <select
            id="scf-system-select"
            value={systemId}
            onChange={(e) => setSystemId(e.target.value)}
            disabled={controlsDisabled}
            className="w-full px-3 py-2 border border-slate-300 rounded-lg bg-white text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
            aria-label="Select molecular system"
          >
            {SYSTEM_PRESETS.map((preset) => (
              <option key={preset.id} value={preset.id}>
                {preset.label}
              </option>
            ))}
          </select>
          {currentSystem && (
            <div className="mt-2 text-xs text-slate-500">
              <p>{currentSystem.description}</p>
              <p className="mt-1">
                <span className="font-medium">{currentSystem.nbf}</span> basis functions,{' '}
                <span className="font-medium">{currentSystem.nelec}</span> electrons
              </p>
            </div>
          )}
        </div>
      )}

      {/* Custom Mode: Geometry Editor & Controls */}
      {isCustomMode && (
        <>
          {/* Example Geometry Selector */}
          <div className="mb-4">
            <ExampleGeometrySelector onSelect={handleExampleSelect} disabled={controlsDisabled} />
          </div>

          {/* Geometry Editor */}
          <div className="mb-4">
            <GeometryEditor
              value={customGeometryText}
              onChange={handleGeometryTextChange}
              units={customUnits}
              onUnitsChange={handleUnitsChange}
              disabled={controlsDisabled}
              placeholder="H  0.0  0.0  0.0&#10;H  0.0  0.0  1.4"
            />
          </div>

          {/* Basis Set Selector */}
          <div className="mb-4">
            <BasisSetSelector
              value={inputMode.mode === 'custom' ? inputMode.basisSet : 'sto-3g'}
              onChange={(basis: BasisSetName) => setCustomBasisSet(basis)}
              disabled={controlsDisabled}
            />
          </div>

          {/* Basis Type Toggle (Cartesian vs Spherical) */}
          <div className="mb-4">
            <BasisTypeToggle
              useSpherical={useSpherical}
              onChange={setUseSpherical}
              disabled={controlsDisabled}
            />
          </div>

          {/* Compute Integrals Button */}
          <div className="mb-6">
            <IntegralComputeButton
              onCompute={computeIntegrals}
              onCancel={cancelIntegrals}
              disabled={!integralReady || !isGeometryValid || controlsDisabled}
              isComputing={isComputingIntegrals}
              progress={progress}
            />

            {/* Integral computation status */}
            {integralCompute.status === 'success' && (
              <div className="mt-2 flex items-center gap-2 text-green-700 text-sm">
                <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
                  <path
                    fillRule="evenodd"
                    d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z"
                    clipRule="evenodd"
                  />
                </svg>
                <span>
                  Integrals ready: {integralCompute.systemData.nbf} basis functions,{' '}
                  {integralCompute.systemData.nelec} electrons
                </span>
              </div>
            )}
            {integralCompute.status === 'error' && (
              <div className="mt-2 flex items-center gap-2 text-red-700 text-sm">
                <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
                  <path
                    fillRule="evenodd"
                    d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z"
                    clipRule="evenodd"
                  />
                </svg>
                <span>{integralCompute.error}</span>
              </div>
            )}
          </div>
        </>
      )}

      {/* Convergence Profile */}
      <div className="mb-6">
        <label className="block text-sm font-medium text-slate-700 mb-2">Convergence Profile</label>
        <div className="flex gap-2">
          {CONVERGENCE_OPTIONS.map(({ value, label }) => (
            <button
              key={value}
              type="button"
              onClick={() => setConvergenceProfile(value)}
              disabled={controlsDisabled}
              className={`flex-1 px-3 py-2 rounded-lg text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
                convergenceProfile === value
                  ? 'bg-blue-600 text-white'
                  : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
              } ${controlsDisabled ? 'opacity-50 cursor-not-allowed' : ''}`}
              aria-pressed={convergenceProfile === value}
              title={`Energy: ${CONVERGENCE_THRESHOLDS[value].energy}, Density: ${CONVERGENCE_THRESHOLDS[value].density}`}
            >
              {label}
            </button>
          ))}
        </div>
        <p className="mt-1 text-xs text-slate-500">
          <Math>{String.raw`\Delta E`}</Math>: {CONVERGENCE_THRESHOLDS[convergenceProfile].energy},{' '}
          <Math>{String.raw`\|\Delta\mathbf{P}\|`}</Math>: {CONVERGENCE_THRESHOLDS[convergenceProfile].density}
        </p>
      </div>

      {/* DIIS Toggle */}
      <div className="mb-6">
        <label className="flex items-center gap-3 cursor-pointer">
          <input
            type="checkbox"
            checked={useDiis}
            onChange={(e) => setUseDiis(e.target.checked)}
            disabled={controlsDisabled}
            className="w-4 h-4 rounded border-slate-300 text-blue-600 focus:ring-2 focus:ring-blue-500 disabled:opacity-50"
            aria-label="Enable DIIS acceleration"
          />
          <span className="text-sm font-medium text-slate-700">Enable DIIS Acceleration</span>
        </label>
        <p className="mt-1 text-xs text-slate-500 ml-7">
          Direct Inversion in the Iterative Subspace: extrapolates{' '}
          <Math>{String.raw`\mathbf{F}' = \sum_i c_i \mathbf{F}_i`}</Math>
        </p>
      </div>

      {/* Damping Control */}
      <div className="mb-6">
        <div className="flex items-center justify-between mb-2">
          <label
            htmlFor="scf-damp-select"
            className="block text-sm font-medium text-slate-700"
          >
            Damping
          </label>
          <select
            id="scf-damp-select"
            value={damp}
            onChange={(e) => setDamp(parseFloat(e.target.value))}
            disabled={controlsDisabled}
            className="px-3 py-1.5 border border-slate-300 rounded-lg bg-white text-slate-800 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
            aria-label="Select damping level"
          >
            <option value={0}>Off</option>
            <option value={0.3}>Light (0.3)</option>
            <option value={0.5}>Medium (0.5)</option>
            <option value={0.7}>Heavy (0.7)</option>
          </select>
        </div>
        <p className="text-xs text-slate-500">
          Mixes old and new Fock matrices to stabilize convergence for difficult systems (e.g., diffuse basis sets)
        </p>
      </div>

      {/* Max Iterations */}
      <div className="mb-6">
        <label
          htmlFor="scf-max-iter"
          className="block text-sm font-medium text-slate-700 mb-2"
        >
          Maximum Iterations
        </label>
        <input
          id="scf-max-iter"
          type="number"
          min={1}
          max={200}
          step={1}
          value={maxIterations}
          onChange={(e) => setMaxIterations(Number(e.target.value))}
          disabled={controlsDisabled}
          className="w-full px-3 py-2 border border-slate-300 rounded-lg text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
          aria-label="Maximum number of SCF iterations"
        />
        <p className="mt-1 text-xs text-slate-500">Range: 1-200 iterations</p>
      </div>

      {/* Run/Cancel Buttons */}
      <div className="flex gap-3">
        <button
          type="button"
          onClick={run}
          disabled={disabled || !isReady || isRunning || !canRunScf}
          className={`flex-1 px-4 py-3 rounded-lg font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
            isRunning
              ? 'bg-blue-400 text-white cursor-not-allowed'
              : 'bg-blue-600 text-white hover:bg-blue-700'
          } disabled:opacity-50 disabled:cursor-not-allowed`}
          aria-label="Run SCF calculation"
          title={
            isCustomMode && !hasIntegrals
              ? 'Compute integrals first before running SCF'
              : undefined
          }
        >
          {isRunning ? (
            <span className="flex items-center justify-center gap-2">
              <span className="animate-spin rounded-full h-4 w-4 border-b-2 border-white" />
              Running...
            </span>
          ) : (
            'Run SCF'
          )}
        </button>
        <button
          type="button"
          onClick={cancelRun}
          disabled={!isRunning}
          className="px-4 py-3 rounded-lg font-medium transition-colors bg-slate-100 text-slate-700 hover:bg-slate-200 focus:outline-none focus-visible:ring-2 focus-visible:ring-slate-500 disabled:opacity-50 disabled:cursor-not-allowed"
          aria-label="Cancel SCF calculation"
        >
          Cancel
        </button>
      </div>

      {/* Custom mode hint when integrals not computed */}
      {isCustomMode && !hasIntegrals && !isComputingIntegrals && isGeometryValid && (
        <p className="mt-2 text-xs text-amber-600">
          Click &quot;Compute Integrals&quot; before running SCF
        </p>
      )}

      {/* Status Indicator */}
      {compute.status !== 'idle' && compute.status !== 'running' && (
        <div className="mt-4">
          {compute.status === 'success' && (
            <div className="flex items-center gap-2 text-green-700 text-sm">
              <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
                <path
                  fillRule="evenodd"
                  d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z"
                  clipRule="evenodd"
                />
              </svg>
              {compute.result.converged ? 'Converged' : 'Did not converge'} in{' '}
              {compute.result.iterations} iterations
            </div>
          )}
          {compute.status === 'error' && (
            <div className="flex items-center gap-2 text-red-700 text-sm">
              <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
                <path
                  fillRule="evenodd"
                  d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z"
                  clipRule="evenodd"
                />
              </svg>
              Error: {compute.error}
            </div>
          )}
          {compute.status === 'cancelled' && (
            <div className="flex items-center gap-2 text-amber-700 text-sm">
              <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
                <path
                  fillRule="evenodd"
                  d="M10 18a8 8 0 100-16 8 8 0 000 16zM8 7a1 1 0 00-1 1v4a1 1 0 001 1h4a1 1 0 001-1V8a1 1 0 00-1-1H8z"
                  clipRule="evenodd"
                />
              </svg>
              Calculation cancelled
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export default ScfControlsPanel;
