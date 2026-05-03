/**
 * ScfControlsPanel - Parameter input controls for SCF computation.
 *
 * Tab-based wizard layout that guides novice users step-by-step
 * while remaining compact for experienced users.
 *
 * **Tab 1 "System":** Molecule, basis set, and method selectors.
 *   In custom mode, also shows geometry editor and basis set selector.
 *
 * **Tab 2 "Options":** Convergence profile, DIIS, damping, max iterations,
 *   and basis type toggle. All have sensible defaults so novices can skip.
 *
 * **Run/Cancel** buttons are always visible at the bottom regardless of tab.
 *
 * Connected to the scfStore for state management.
 *
 * @module components/scf/ScfControlsPanel
 */

import { useCallback, useMemo, useState } from 'react';
import {
  useScfStore,
  CONVERGENCE_THRESHOLDS,
  MOLECULE_PRESETS,
  BASIS_SET_OPTIONS,
  getMoleculePreset,
  getPresetIdForMoleculeBasis,
} from '../../stores/scfStore';
import { useScf } from '../../hooks/useScf';
import { parseAndValidate } from '../../lib/geometryValidation';
import { Math } from '../common/Math';
import { MethodSelector } from './MethodSelector';
import { isDftMethod, getRunButtonLabel } from '../../types/dft';
import type { ConvergenceProfile, BasisSetName, ScfIterationHistory } from '../../worker/protocol';
import type { ExampleGeometry } from '../../lib/exampleGeometries';
import { geometryToXyzText } from '../../lib/geometryValidation';
import type { RunOptions } from '../../hooks/useScf';

// ============================================================================
// SCF Progress Bar
// ============================================================================

/** Human-readable labels for integral computation steps. */
const INTEGRAL_STEP_LABELS: Record<string, string> = {
  overlap: 'Overlap (S)',
  hcore: 'Core Hamiltonian (H)',
  eri: 'Electron repulsion integrals (ERI)',
  grid: 'Numerical integration grid',
  done: 'Integrals complete',
};

/**
 * SCF progress bar showing computation phases with live progress.
 *
 * Two-phase display:
 * 1. Integrals: determinate bar showing S → Hcore → ERI → Grid progress
 * 2. SCF iterations: determinate bar with iteration count and live energy
 */
function ScfProgressBar({
  history,
  maxIterations,
  integralPercent,
  integralStep,
}: {
  history: ScfIterationHistory[];
  maxIterations: number;
  integralPercent: number;
  integralStep: string;
}) {
  const inScfPhase = history.length > 0;
  const iterProgress = inScfPhase ? (history.length / maxIterations) * 100 : 0;
  const stepLabel = INTEGRAL_STEP_LABELS[integralStep] || integralStep;

  return (
    <div className="mt-3 space-y-1.5">
      {/* Phase label row */}
      <div className="flex items-center justify-between text-xs">
        {inScfPhase ? (
          <>
            <span className="text-blue-700 font-medium flex items-center gap-1.5">
              <span className="w-1.5 h-1.5 rounded-full bg-blue-600 animate-pulse" />
              SCF iteration
              <span className="font-mono">{history.length}/{maxIterations}</span>
            </span>
            <span className="text-green-600 flex items-center gap-1">
              <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20"><path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" /></svg>
              Integrals
            </span>
          </>
        ) : (
          <>
            <span className="text-blue-700 font-medium flex items-center gap-1.5">
              <span className="w-1.5 h-1.5 rounded-full bg-blue-600 animate-pulse" />
              {integralPercent > 0 ? stepLabel : 'Initializing...'}
            </span>
            {integralPercent > 0 && (
              <span className="text-slate-500 font-mono">{integralPercent.toFixed(0)}%</span>
            )}
          </>
        )}
      </div>

      {/* Progress bar */}
      <div className="h-1.5 bg-slate-200 rounded-full overflow-hidden">
        <div
          className={`h-full rounded-full transition-all duration-300 ${
            inScfPhase ? 'bg-blue-600' : integralPercent > 0 ? 'bg-amber-500' : 'bg-blue-500 animate-progress-indeterminate'
          }`}
          style={
            integralPercent > 0 || inScfPhase
              ? { width: `${globalThis.Math.min(inScfPhase ? iterProgress : integralPercent, 100)}%` }
              : undefined
          }
        />
      </div>

      {/* Live energy readout (during SCF phase) */}
      {inScfPhase && history.length > 0 && (
        <p className="text-xs text-slate-500 font-mono">
          E = {history[history.length - 1].energy.toFixed(10)} Ha
          {history.length >= 2 && (
            <span className="ml-2">
              ΔE = {history[history.length - 1].delta.toExponential(2)}
            </span>
          )}
        </p>
      )}
    </div>
  );
}

// Import custom geometry components
import {
  InputModeToggle,
  GeometryEditor,
  BasisSetSelector,
  BasisTypeToggle,
  ExampleGeometrySelector,
} from './custom';

/**
 * Props for ScfControlsPanel.
 */
interface ScfControlsPanelProps {
  /** Disable all controls */
  disabled?: boolean;
  /** Atom coordinates for DFT path [Z, x, y, z][] */
  workerAtoms?: [number, number, number, number][] | null;
  /** Basis set name for DFT path */
  basisName?: string;
  /** Hide the Run/Cancel buttons section (used in Optimize tab where OptimizeButton replaces it) */
  hideRunButton?: boolean;
}

/**
 * Active tab in the Parameters panel.
 */
type ControlsTab = 'system' | 'options';

/**
 * Convergence profile options for the button group.
 */
const CONVERGENCE_OPTIONS: { value: ConvergenceProfile; label: string }[] = [
  { value: 'loose', label: 'Loose' },
  { value: 'medium', label: 'Medium' },
  { value: 'tight', label: 'Tight' },
];

/**
 * Props for MoleculeBasisSelector.
 */
interface MoleculeBasisSelectorProps {
  selectedMolecule: string;
  selectedBasis: string;
  onMoleculeChange: (id: string) => void;
  onBasisChange: (basis: BasisSetName) => void;
  disabled?: boolean;
}

/**
 * Independent Molecule + Basis Set selector for preset mode.
 *
 * When a pre-computed preset exists for the combination (e.g., H2 + STO-3G),
 * it uses the fast preset path. Otherwise, it switches to on-the-fly
 * integral computation automatically.
 *
 * HeH+ (charged system) is disabled for non-STO-3G basis sets because
 * charged systems require the pre-computed STO-3G preset.
 */
function MoleculeBasisSelector({
  selectedMolecule,
  selectedBasis,
  onMoleculeChange,
  onBasisChange,
  disabled = false,
}: MoleculeBasisSelectorProps) {
  const currentMolecule = getMoleculePreset(selectedMolecule);
  const hasPreset = !!getPresetIdForMoleculeBasis(selectedMolecule, selectedBasis);

  // HeH+ is charged -- only works with pre-computed STO-3G preset
  const isChargedMolecule = currentMolecule ? currentMolecule.charge !== 0 : false;

  return (
    <div className="space-y-3">
      {/* Molecule selector */}
      <div>
        <label
          htmlFor="scf-molecule-select"
          className="block text-sm font-medium text-slate-700 mb-1"
        >
          Molecule
        </label>
        <select
          id="scf-molecule-select"
          value={selectedMolecule}
          onChange={(e) => onMoleculeChange(e.target.value)}
          disabled={disabled}
          className="w-full px-3 py-2 border border-slate-300 rounded-lg bg-white text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
          aria-label="Select molecule"
        >
          {MOLECULE_PRESETS.map((mol) => (
            <option key={mol.id} value={mol.id}>
              {mol.label}
            </option>
          ))}
        </select>
      </div>

      {/* Basis set selector */}
      <div>
        <label
          htmlFor="scf-basis-select"
          className="block text-sm font-medium text-slate-700 mb-1"
        >
          Basis Set
        </label>
        <select
          id="scf-basis-select"
          value={selectedBasis}
          onChange={(e) => onBasisChange(e.target.value as BasisSetName)}
          disabled={disabled}
          className="w-full px-3 py-2 border border-slate-300 rounded-lg bg-white text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
          aria-label="Select basis set"
        >
          {BASIS_SET_OPTIONS.map((basis) => {
            // HeH+ only works with STO-3G preset
            const basisDisabled = isChargedMolecule && basis.id !== 'sto-3g';
            return (
              <option
                key={basis.id}
                value={basis.id}
                disabled={basisDisabled}
              >
                {basis.label}{basisDisabled ? ' (preset only for charged systems)' : ''}
              </option>
            );
          })}
        </select>
      </div>

      {/* System info */}
      {currentMolecule && (
        <div className="text-xs text-slate-500 space-y-0.5">
          <p>
            {currentMolecule.description} ({currentMolecule.atoms.length} atom{currentMolecule.atoms.length !== 1 ? 's' : ''}, {currentMolecule.nelec} electrons)
          </p>
          {!hasPreset && (
            <p className="text-amber-600 flex items-center gap-1">
              <svg className="w-3.5 h-3.5 flex-shrink-0" fill="currentColor" viewBox="0 0 20 20">
                <path fillRule="evenodd" d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z" clipRule="evenodd" />
              </svg>
              On-the-fly integrals (may be slower)
            </p>
          )}
        </div>
      )}
    </div>
  );
}

/**
 * Checkmark icon shown on tabs with valid configuration.
 */
function TabCheckmark() {
  return (
    <svg
      className="w-3.5 h-3.5 text-green-500"
      fill="currentColor"
      viewBox="0 0 20 20"
      aria-hidden="true"
    >
      <path
        fillRule="evenodd"
        d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
        clipRule="evenodd"
      />
    </svg>
  );
}

/**
 * Controls panel for SCF parameters.
 *
 * Uses a tab-based layout with two tabs:
 * - **System**: molecule, basis, method, and custom geometry controls
 * - **Options**: convergence, DIIS, damping, max iterations
 *
 * Run/Cancel buttons are always visible at the bottom.
 *
 * @example
 * ```tsx
 * <ScfControlsPanel />
 * ```
 */
export function ScfControlsPanel({ disabled = false, workerAtoms, basisName: basisNameProp, hideRunButton = false }: ScfControlsPanelProps) {
  // Local UI state for active tab
  const [activeTab, setActiveTab] = useState<ControlsTab>('system');

  // Store state
  const selectedMolecule = useScfStore((state) => state.selectedMolecule);
  const selectedBasis = useScfStore((state) => state.selectedBasis);
  const convergenceProfile = useScfStore((state) => state.convergenceProfile);
  const maxIterations = useScfStore((state) => state.maxIterations);
  const useDiis = useScfStore((state) => state.useDiis);
  const damp = useScfStore((state) => state.damp);
  const useSpherical = useScfStore((state) => state.useSpherical);
  const compute = useScfStore((state) => state.compute);
  const history = useScfStore((state) => state.history);
  const integralPercent = useScfStore((state) => state.integralPercent);
  const integralStep = useScfStore((state) => state.integralStep);
  const inputMode = useScfStore((state) => state.inputMode);
  const customGeometryText = useScfStore((state) => state.customGeometryText);
  const customUnits = useScfStore((state) => state.customUnits);
  const method = useScfStore((state) => state.method);
  const setMethod = useScfStore((state) => state.setMethod);
  const gridQuality = useScfStore((state) => state.gridQuality);
  const setGridQuality = useScfStore((state) => state.setGridQuality);
  const levelShift = useScfStore((state) => state.levelShift);
  const setLevelShift = useScfStore((state) => state.setLevelShift);

  // Store actions
  const setSelectedMolecule = useScfStore((state) => state.setSelectedMolecule);
  const setSelectedBasis = useScfStore((state) => state.setSelectedBasis);
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

  // Derived state
  const isPresetMode = inputMode.mode === 'preset';
  const isCustomMode = inputMode.mode === 'custom';

  // Validate custom geometry text
  const validationResult = useMemo(
    () => (isCustomMode ? parseAndValidate(customGeometryText) : null),
    [isCustomMode, customGeometryText]
  );

  const isGeometryValid = validationResult?.valid ?? false;

  // All SCF paths (RHF, LDA, B3LYP) now compute integrals on-the-fly via
  // the ks_scf WASM path when workerAtoms are available, so custom mode only
  // requires valid geometry.
  const canRunScf = isPresetMode || (isCustomMode && isGeometryValid);

  // Determine if controls should be disabled
  const controlsDisabled = disabled || isRunning;

  // System tab validity: a valid molecule is selected (always true for preset mode,
  // requires valid geometry for custom mode)
  const systemConfigured = isPresetMode || (isCustomMode && isGeometryValid);

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

  // Tab definitions
  const tabs: { id: ControlsTab; label: string; showCheck: boolean }[] = [
    { id: 'system', label: 'System', showCheck: systemConfigured },
    { id: 'options', label: 'Options', showCheck: false },
  ];

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
      <h2 className="text-lg font-semibold text-slate-800 mb-3">Parameters</h2>

      {/* Tab Bar */}
      <div
        className="flex gap-1 mb-4 bg-slate-100 rounded-lg p-1"
        role="tablist"
        aria-label="Parameter sections"
      >
        {tabs.map(({ id, label, showCheck }) => (
          <button
            key={id}
            type="button"
            role="tab"
            id={`tab-${id}`}
            aria-selected={activeTab === id}
            aria-controls={`panel-${id}`}
            onClick={() => setActiveTab(id)}
            className={`flex-1 flex items-center justify-center gap-1.5 px-3 py-1.5 rounded-md text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
              activeTab === id
                ? 'bg-white text-blue-700 shadow-sm'
                : 'text-slate-600 hover:text-slate-800'
            }`}
          >
            {label}
            {showCheck && <TabCheckmark />}
          </button>
        ))}
      </div>

      {/* Tab Panels */}

      {/* System Tab */}
      <div
        id="panel-system"
        role="tabpanel"
        aria-labelledby="tab-system"
        hidden={activeTab !== 'system'}
      >
        {activeTab === 'system' && (
          <div className="space-y-4">
            {/* Input Mode Toggle */}
            <InputModeToggle
              mode={inputMode.mode}
              onChange={handleModeChange}
              disabled={controlsDisabled}
            />

            {/* Method Selector */}
            <MethodSelector
              method={method}
              onChange={setMethod}
              disabled={controlsDisabled}
            />

            {/* Molecule + Basis Selectors (preset mode only) */}
            {isPresetMode && (
              <MoleculeBasisSelector
                selectedMolecule={selectedMolecule}
                selectedBasis={selectedBasis}
                onMoleculeChange={setSelectedMolecule}
                onBasisChange={setSelectedBasis}
                disabled={controlsDisabled}
              />
            )}

            {/* Custom Mode: Geometry Editor & Controls */}
            {isCustomMode && (
              <div className="space-y-4 pt-2 border-t border-slate-200">
                {/* Example Geometry Selector */}
                <ExampleGeometrySelector onSelect={handleExampleSelect} disabled={controlsDisabled} />

                {/* Geometry Editor */}
                <GeometryEditor
                  value={customGeometryText}
                  onChange={handleGeometryTextChange}
                  units={customUnits}
                  onUnitsChange={handleUnitsChange}
                  disabled={controlsDisabled}
                  placeholder="H  0.0  0.0  0.0&#10;H  0.0  0.0  1.4"
                />

                {/* Basis Set Selector */}
                <BasisSetSelector
                  value={inputMode.mode === 'custom' ? inputMode.basisSet : 'sto-3g'}
                  onChange={(basis: BasisSetName) => setCustomBasisSet(basis)}
                  disabled={controlsDisabled}
                />
              </div>
            )}
          </div>
        )}
      </div>

      {/* Options Tab */}
      <div
        id="panel-options"
        role="tabpanel"
        aria-labelledby="tab-options"
        hidden={activeTab !== 'options'}
      >
        {activeTab === 'options' && (
          <div className="space-y-5">
            {/* Convergence Profile */}
            <div>
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
            <div>
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
            <div>
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

            {/* Level Shift */}
            <div>
              <div className="flex items-center justify-between mb-2">
                <label
                  htmlFor="scf-level-shift"
                  className="block text-sm font-medium text-slate-700"
                >
                  Level Shift
                </label>
                <select
                  id="scf-level-shift"
                  value={levelShift}
                  onChange={(e) => setLevelShift(parseFloat(e.target.value))}
                  disabled={controlsDisabled}
                  className="px-3 py-1.5 border border-slate-300 rounded-lg bg-white text-slate-800 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
                  aria-label="Select level shift"
                >
                  <option value={0}>Off</option>
                  <option value={0.25}>Mild (0.25)</option>
                  <option value={0.5}>Standard (0.5)</option>
                  <option value={1.0}>Strong (1.0)</option>
                </select>
              </div>
              <p className="text-xs text-slate-500">
                Shifts virtual orbitals up to stabilize convergence for systems with small HOMO-LUMO gaps
              </p>
            </div>

            {/* Max Iterations */}
            <div>
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

            {/* Basis Type Toggle (Cartesian vs Spherical) */}
            <BasisTypeToggle
              useSpherical={useSpherical}
              onChange={setUseSpherical}
              disabled={controlsDisabled}
            />

            {/* DFT Grid Quality (only shown for DFT methods) */}
            {isDftMethod(method) && (
              <div>
                <label
                  htmlFor="scf-grid-quality"
                  className="block text-sm font-medium text-slate-700 mb-2"
                >
                  DFT Grid Quality
                </label>
                <select
                  id="scf-grid-quality"
                  value={gridQuality}
                  onChange={(e) => setGridQuality(e.target.value as 'standard' | 'fine')}
                  disabled={controlsDisabled}
                  className="w-full px-3 py-2 border border-slate-300 rounded-lg bg-white text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
                  aria-label="Select DFT grid quality"
                >
                  <option value="standard">Standard (~4,000 pts/atom)</option>
                  <option value="fine">Fine (~8,000 pts/atom)</option>
                </select>
                <p className="mt-1 text-xs text-slate-500">
                  Higher quality grids give more accurate DFT energies but are slower
                </p>
              </div>
            )}
          </div>
        )}
      </div>

      {/* Divider before Run buttons (hidden when used in Optimize tab) */}
      {!hideRunButton && (
      <div className="border-t border-slate-200 mt-5 pt-4">
        {/* Run/Cancel Buttons - always visible */}
        <div className="flex gap-3">
          <button
            type="button"
            onClick={() => {
              const options: RunOptions = {};
              // Pass atoms/basisName for on-the-fly computation:
              // - Always for DFT methods (require ks_scf path)
              // - For RHF when no pre-computed preset exists (on-the-fly integrals)
              const presetExists = !!getPresetIdForMoleculeBasis(selectedMolecule, selectedBasis);
              if (workerAtoms && (isDftMethod(method) || !presetExists)) {
                options.atoms = workerAtoms;
                options.basisName = basisNameProp;
              }
              run(options);
            }}
            disabled={disabled || !isReady || isRunning || !canRunScf}
            className={`flex-1 px-4 py-3 rounded-lg font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
              isRunning
                ? 'bg-blue-400 text-white cursor-not-allowed'
                : isDftMethod(method)
                  ? 'bg-teal-600 text-white hover:bg-teal-700'
                  : 'bg-blue-600 text-white hover:bg-blue-700'
            } disabled:opacity-50 disabled:cursor-not-allowed`}
            aria-label={`Run ${isDftMethod(method) ? 'KS-DFT' : 'SCF'} calculation`}
          >
            {isRunning ? (
              <span className="flex items-center justify-center gap-2">
                <span className="animate-spin rounded-full h-4 w-4 border-b-2 border-white" />
                Running...
              </span>
            ) : (
              getRunButtonLabel(method)
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

        {/* Progress Bar (during computation) */}
        {compute.status === 'running' && (
          <ScfProgressBar
            history={history}
            maxIterations={maxIterations}
            integralPercent={integralPercent}
            integralStep={integralStep}
          />
        )}

        {/* Status Indicator */}
        {compute.status !== 'idle' && compute.status !== 'running' && (
          <div className="mt-3">
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
      )}
    </div>
  );
}

export default ScfControlsPanel;
