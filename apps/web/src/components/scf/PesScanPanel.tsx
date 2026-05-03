/**
 * PesScanPanel - PES scan controls and progress.
 *
 * Provides a collapsible panel for configuring and running PES
 * (Potential Energy Surface) scans. Supports both the legacy diatomic
 * bond-length scan path and the generalized internal coordinate scan
 * (bond/angle/dihedral) on arbitrary molecules.
 *
 * Features:
 * - CoordinateSelector for choosing coordinate type and atoms
 * - ScanModeToggle for rigid vs relaxed scan mode
 * - Scan range inputs (min/max/nPoints)
 * - Input validation with inline error messages
 * - Start/Cancel buttons with loading state
 * - Progress bar showing completed/total points
 *
 * For diatomic molecules (2 atoms), automatically selects Bond with
 * atoms [0, 1] and hides the coordinate selector for simplicity.
 *
 * @module components/scf/PesScanPanel
 * @see US-040 PES Scan UI (original diatomic)
 * @see US-082 Coordinate Selector + Scan Mode UI
 */

import { useState, useCallback, useMemo, useEffect } from 'react';
import { useScfStore } from '../../stores/scfStore';
import { usePesScan } from '../../hooks/usePesScan';
import { validatePesConfig } from './pesValidation';
import {
  CoordinateSelector,
  computeCoordinateValue,
  COORDINATE_LABELS,
  ELEMENT_SYMBOLS,
} from './CoordinateSelector';
import { ScanModeToggle } from './ScanModeToggle';
import type { CoordinateType } from './CoordinateSelector';
import type { ScanMode } from './ScanModeToggle';

// ============================================================================
// Props
// ============================================================================

/**
 * Atom info extracted from the current system for PES scan.
 */
export interface DiatomicInfo {
  /** Element symbol of atom A (fixed at origin) */
  symbolA: string;
  /** Element symbol of atom B (translated along z-axis) */
  symbolB: string;
  /** Atomic number of atom A */
  atomAZ: number;
  /** Atomic number of atom B */
  atomBZ: number;
  /** Basis set name for the current system */
  basisName: string;
}

/**
 * Props for PesScanPanel.
 */
interface PesScanPanelProps {
  /** Diatomic molecule information (null if not diatomic or data unavailable) */
  diatomic: DiatomicInfo | null;
  /** Atoms from the current geometry as [Z, x, y, z] arrays */
  workerAtoms: [number, number, number, number][] | null;
  /** Basis set name for the current system */
  basisName: string;
  /** Electronic structure method */
  method: string;
  /** Whether to use spherical harmonics */
  useSpherical: boolean;
  /** Disable all controls (e.g., worker not ready) */
  disabled?: boolean;
}

// ============================================================================
// Helpers
// ============================================================================

/**
 * Compute default scan range for a given coordinate type and current value.
 *
 * Bond: current +/- 1 bohr (min clamped to 0.5)
 * Angle: current +/- 30 degrees (clamped to [1, 179])
 * Dihedral: full rotation -180 to 180
 */
function getDefaultRange(
  coordinateType: CoordinateType,
  currentValue: number | null,
): { min: number; max: number; nPoints: number } {
  if (coordinateType === 'bond') {
    const cv = currentValue ?? 2.0;
    return {
      min: Math.max(0.5, cv - 1.0),
      max: cv + 1.0,
      nPoints: 20,
    };
  }
  if (coordinateType === 'angle') {
    const cv = currentValue ?? 109.5;
    return {
      min: Math.max(1, cv - 30),
      max: Math.min(179, cv + 30),
      nPoints: 20,
    };
  }
  // dihedral
  return {
    min: -180,
    max: 180,
    nPoints: 36,
  };
}

/**
 * Validate scan range for internal coordinate scans.
 *
 * Extends the basic PES validation with coordinate-type-specific checks.
 */
function validateInternalScanConfig(
  coordinateType: CoordinateType,
  valueMin: number,
  valueMax: number,
  nPoints: number,
  currentValue: number | null,
): string | null {
  if (Number.isNaN(valueMin) || Number.isNaN(valueMax) || Number.isNaN(nPoints)) {
    return 'All inputs must be valid numbers';
  }
  if (nPoints < 2) return 'Need at least 2 scan points';
  if (nPoints > 100) return 'Maximum 100 scan points';

  if (coordinateType === 'bond') {
    if (valueMin < 0.3) return 'Minimum distance must be at least 0.3 bohr';
    if (valueMin >= valueMax) return 'Minimum must be less than maximum';
  } else if (coordinateType === 'angle') {
    if (valueMin < 0 || valueMax > 180) return 'Angle must be between 0 and 180 degrees';
    if (valueMin >= valueMax) return 'Minimum must be less than maximum';
  } else {
    // Dihedral: min can be > max for wrapping, but they cannot be equal
    if (valueMin === valueMax) return 'Min and max cannot be equal';
  }

  if (currentValue === null) {
    return 'Select distinct atoms to define the coordinate';
  }

  return null;
}

// ============================================================================
// Component
// ============================================================================

/**
 * PES scan configuration and progress panel.
 *
 * Supports two paths:
 * 1. **Legacy diatomic** (2 atoms): uses the old PesScanRequest (bond only)
 * 2. **Internal coordinate** (2+ atoms): uses PesScanInternalRequest
 *
 * For diatomics, the coordinate selector is simplified (auto-selects bond
 * with atoms [0, 1]). For polyatomics, the full coordinate selector is shown.
 *
 * @example
 * ```tsx
 * <PesScanPanel
 *   diatomic={diatomicInfo}
 *   workerAtoms={workerAtoms}
 *   basisName={basisName}
 *   method={method}
 *   useSpherical={useSpherical}
 *   disabled={!isReady}
 * />
 * ```
 */
export function PesScanPanel({
  diatomic,
  workerAtoms,
  basisName,
  method,
  useSpherical,
  disabled = false,
}: PesScanPanelProps) {
  const atoms = workerAtoms ?? [];
  const nAtoms = atoms.length;
  const isDiatomic = nAtoms === 2;

  // Store state (legacy PES)
  const pesState = useScfStore((state) => state.pesState);
  const setPesConfig = useScfStore((state) => state.setPesConfig);
  const setPesCoordinateConfig = useScfStore((state) => state.setPesCoordinateConfig);
  const resetPes = useScfStore((state) => state.resetPes);

  // Hook for worker communication
  const { isReady, isScanning, startScan, startInternalScan, cancelScan } = usePesScan();

  // Collapsible panel state
  const [isExpanded, setIsExpanded] = useState(true);

  // Internal coordinate state (local to this panel, seeded from URL if available, US-085)
  const [coordinateType, setCoordinateType] = useState<CoordinateType>(
    pesState.urlCoordinateType ?? 'bond'
  );
  const [atomIndices, setAtomIndices] = useState<number[]>(
    pesState.urlAtomIndices ?? [0, 1]
  );
  const [scanMode, setScanMode] = useState<ScanMode>(
    pesState.urlScanMode ?? 'rigid'
  );

  // Relaxed scan optimization settings
  const [optMaxSteps, setOptMaxSteps] = useState(50);
  const [optGradThreshold, setOptGradThreshold] = useState(4.5e-4);
  const [localMin, setLocalMin] = useState<number>(
    pesState.urlValueMin ?? 0.5
  );
  const [localMax, setLocalMax] = useState<number>(
    pesState.urlValueMax ?? 5.0
  );
  const [localNPoints, setLocalNPoints] = useState<number>(pesState.nPoints);

  // For legacy diatomic path
  const [legacyRMin, setLegacyRMin] = useState(pesState.rMin);
  const [legacyRMax, setLegacyRMax] = useState(pesState.rMax);
  const [legacyNPoints, setLegacyNPoints] = useState(pesState.nPoints);

  // Clear URL-seeded values after consuming them (one-time, US-085)
  useEffect(() => {
    if (pesState.urlCoordinateType !== null) {
      useScfStore.setState((state) => ({
        pesState: {
          ...state.pesState,
          urlCoordinateType: null,
          urlAtomIndices: null,
          urlScanMode: null,
          urlValueMin: null,
          urlValueMax: null,
        },
      }));
    }
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  // Sync coordinate config to store for URL encoding (US-085)
  useEffect(() => {
    if (!isDiatomic) {
      setPesCoordinateConfig({
        coordinateType,
        atomIndices,
        scanMode,
        valueMin: localMin,
        valueMax: localMax,
        nPoints: localNPoints,
      });
    } else {
      // Diatomic mode: clear coordinate config so legacy path is used
      setPesCoordinateConfig(null);
    }
  }, [coordinateType, atomIndices, scanMode, localMin, localMax, localNPoints, isDiatomic, setPesCoordinateConfig]);

  // Current coordinate value
  const currentValue = useMemo(
    () => computeCoordinateValue(atoms, coordinateType, atomIndices),
    [atoms, coordinateType, atomIndices],
  );

  // Initialize defaults when atoms or coordinate type changes
  // Skip if URL-seeded values were applied (avoid overwriting them)
  const urlSeeded = pesState.urlCoordinateType !== null;
  useEffect(() => {
    if (atoms.length < 2) return;
    if (urlSeeded) return; // Don't override URL-seeded values on first render
    const defaults = getDefaultRange(coordinateType, currentValue);
    setLocalMin(defaults.min);
    setLocalMax(defaults.max);
    setLocalNPoints(defaults.nPoints);
  }, [coordinateType, atoms.length]); // eslint-disable-line react-hooks/exhaustive-deps

  // Auto-select bond [0,1] for diatomics
  useEffect(() => {
    if (isDiatomic) {
      setCoordinateType('bond');
      setAtomIndices([0, 1]);
    }
  }, [isDiatomic]);

  // Validation
  const validationError = useMemo(() => {
    if (isDiatomic) {
      return validatePesConfig(legacyRMin, legacyRMax, legacyNPoints);
    }
    return validateInternalScanConfig(coordinateType, localMin, localMax, localNPoints, currentValue);
  }, [isDiatomic, legacyRMin, legacyRMax, legacyNPoints, coordinateType, localMin, localMax, localNPoints, currentValue]);

  const isValid = validationError === null;

  // Controls disabled during scan or when worker not ready
  const controlsDisabled = disabled || !isReady || isScanning;

  // Legacy diatomic handlers
  const handleLegacyRMinChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const value = parseFloat(e.target.value);
      setLegacyRMin(value);
      if (!Number.isNaN(value)) setPesConfig({ rMin: value });
    },
    [setPesConfig],
  );

  const handleLegacyRMaxChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const value = parseFloat(e.target.value);
      setLegacyRMax(value);
      if (!Number.isNaN(value)) setPesConfig({ rMax: value });
    },
    [setPesConfig],
  );

  const handleLegacyNPointsChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const value = parseInt(e.target.value, 10);
      setLegacyNPoints(value);
      if (!Number.isNaN(value)) setPesConfig({ nPoints: value });
    },
    [setPesConfig],
  );

  // Start scan handler
  const handleStartScan = useCallback(() => {
    if (!isValid || controlsDisabled) return;

    if (isDiatomic && diatomic) {
      // Use legacy diatomic path
      startScan(diatomic.atomAZ, diatomic.atomBZ, diatomic.basisName);
    } else if (atoms.length >= 2) {
      // Use internal coordinate path
      startInternalScan({
        atoms,
        basisName,
        method,
        coordinateType,
        atomIndices,
        scanMode,
        valueMin: coordinateType === 'bond' ? localMin : localMin * (Math.PI / 180),
        valueMax: coordinateType === 'bond' ? localMax : localMax * (Math.PI / 180),
        nPoints: localNPoints,
        useSpherical,
        ...(scanMode === 'relaxed' ? {
          optMaxSteps,
          optGradThreshold,
        } : {}),
      });
    }
  }, [
    isValid, controlsDisabled, isDiatomic, diatomic, atoms, basisName, method,
    coordinateType, atomIndices, scanMode, localMin, localMax, localNPoints,
    useSpherical, optMaxSteps, optGradThreshold, startScan, startInternalScan,
  ]);

  // Reset handler
  const handleReset = useCallback(() => {
    resetPes();
    if (isDiatomic) {
      setLegacyRMin(pesState.rMin);
      setLegacyRMax(pesState.rMax);
      setLegacyNPoints(pesState.nPoints);
    }
  }, [resetPes, pesState.rMin, pesState.rMax, pesState.nPoints, isDiatomic]);

  // Progress calculation
  const completedPoints = pesState.results.length;
  const totalPoints = isDiatomic ? pesState.nPoints : localNPoints;
  const progressPercent = totalPoints > 0 ? (completedPoints / totalPoints) * 100 : 0;

  // Header subtitle
  const headerSubtitle = useMemo(() => {
    if (isDiatomic && diatomic) {
      return `(${diatomic.symbolA}-${diatomic.symbolB})`;
    }
    if (atoms.length >= 2) {
      const selectedSymbols = atomIndices
        .slice(0, coordinateType === 'bond' ? 2 : coordinateType === 'angle' ? 3 : 4)
        .map((i) => ELEMENT_SYMBOLS[atoms[i]?.[0] ?? 0] ?? 'X')
        .join('-');
      return `(${COORDINATE_LABELS[coordinateType]}: ${selectedSymbols})`;
    }
    return '';
  }, [isDiatomic, diatomic, atoms, atomIndices, coordinateType]);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      {/* Header - clickable to expand/collapse */}
      <button
        type="button"
        onClick={() => setIsExpanded(!isExpanded)}
        className="w-full px-4 py-3 bg-slate-50 border-b border-slate-200 flex items-center justify-between text-left hover:bg-slate-100 transition-colors"
        aria-expanded={isExpanded}
        aria-controls="pes-scan-content"
      >
        <div className="flex items-center gap-2">
          <h3 className="text-sm font-semibold text-slate-700">
            PES Scan
          </h3>
          <span className="text-xs text-slate-500">
            {headerSubtitle}
          </span>
        </div>
        <svg
          className={`w-4 h-4 text-slate-500 transition-transform ${isExpanded ? 'rotate-180' : ''}`}
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {/* Content */}
      {isExpanded && (
        <div id="pes-scan-content" className="p-4 space-y-4">
          {/* Scan mode toggle (always shown for polyatomics, hidden for diatomics) */}
          {!isDiatomic && (
            <>
              <ScanModeToggle
                scanMode={scanMode}
                onScanModeChange={setScanMode}
                disabled={controlsDisabled}
              />

              {scanMode === 'relaxed' && (
                <details className="mt-2 group">
                  <summary className="text-xs text-slate-500 cursor-pointer hover:text-slate-700 select-none flex items-center gap-1">
                    <svg className="w-3 h-3 transition-transform group-open:rotate-90" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" /></svg>
                    Advanced options
                  </summary>
                  <div className="mt-2 pl-4 grid grid-cols-[1fr_auto] gap-x-3 gap-y-1.5 items-center">
                    <label htmlFor="pes-opt-max-steps" className="text-xs text-slate-500">
                      Max opt steps per point
                    </label>
                    <input
                      id="pes-opt-max-steps"
                      type="number"
                      min={5}
                      max={200}
                      step={1}
                      value={optMaxSteps}
                      onChange={(e) => setOptMaxSteps(Number(e.target.value))}
                      disabled={controlsDisabled}
                      className="w-16 px-1.5 py-0.5 text-xs border border-slate-300 rounded focus:outline-none focus:ring-1 focus:ring-blue-500 disabled:opacity-50"
                      aria-label="Maximum optimization steps per scan point"
                    />
                    <label htmlFor="pes-opt-grad-threshold" className="text-xs text-slate-500">
                      Convergence (Ha/bohr)
                    </label>
                    <input
                      id="pes-opt-grad-threshold"
                      type="number"
                      step={0.0001}
                      min={0.00001}
                      max={0.01}
                      value={optGradThreshold}
                      onChange={(e) => setOptGradThreshold(Number(e.target.value))}
                      disabled={controlsDisabled}
                      className="w-20 px-1.5 py-0.5 text-xs border border-slate-300 rounded focus:outline-none focus:ring-1 focus:ring-blue-500 disabled:opacity-50"
                      aria-label="Gradient convergence threshold"
                    />
                  </div>
                </details>
              )}

              <CoordinateSelector
                coordinateType={coordinateType}
                onCoordinateTypeChange={setCoordinateType}
                atomIndices={atomIndices}
                onAtomIndicesChange={setAtomIndices}
                atoms={atoms}
                valueMin={localMin}
                onValueMinChange={setLocalMin}
                valueMax={localMax}
                onValueMaxChange={setLocalMax}
                nPoints={localNPoints}
                onNPointsChange={setLocalNPoints}
                disabled={controlsDisabled}
              />
            </>
          )}

          {/* Legacy diatomic controls */}
          {isDiatomic && (
            <>
              {/* Description */}
              <p className="text-xs text-slate-500">
                Scan the potential energy surface by varying the {diatomic?.symbolA ?? '?'}-{diatomic?.symbolB ?? '?'} bond length.
              </p>

              {/* Range Inputs */}
              <div className="grid grid-cols-2 gap-3">
                <div>
                  <label
                    htmlFor="pes-rmin"
                    className="block text-xs font-medium text-slate-700 mb-1"
                  >
                    Min Distance (bohr)
                  </label>
                  <input
                    id="pes-rmin"
                    type="number"
                    min={0.3}
                    max={20.0}
                    step={0.1}
                    value={legacyRMin}
                    onChange={handleLegacyRMinChange}
                    disabled={controlsDisabled}
                    className="w-full px-2 py-1.5 border border-slate-300 rounded-lg text-sm text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
                    aria-label="Minimum bond distance in bohr"
                  />
                </div>

                <div>
                  <label
                    htmlFor="pes-rmax"
                    className="block text-xs font-medium text-slate-700 mb-1"
                  >
                    Max Distance (bohr)
                  </label>
                  <input
                    id="pes-rmax"
                    type="number"
                    min={0.3}
                    max={20.0}
                    step={0.1}
                    value={legacyRMax}
                    onChange={handleLegacyRMaxChange}
                    disabled={controlsDisabled}
                    className="w-full px-2 py-1.5 border border-slate-300 rounded-lg text-sm text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
                    aria-label="Maximum bond distance in bohr"
                  />
                </div>
              </div>

              <div>
                <label
                  htmlFor="pes-npoints"
                  className="block text-xs font-medium text-slate-700 mb-1"
                >
                  Number of Points
                </label>
                <input
                  id="pes-npoints"
                  type="number"
                  min={2}
                  max={100}
                  step={1}
                  value={legacyNPoints}
                  onChange={handleLegacyNPointsChange}
                  disabled={controlsDisabled}
                  className="w-full px-2 py-1.5 border border-slate-300 rounded-lg text-sm text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
                  aria-label="Number of PES scan points"
                />
                <p className="mt-0.5 text-xs text-slate-400">Range: 2-100</p>
              </div>
            </>
          )}

          {/* Validation Error */}
          {validationError && (
            <div
              className="text-xs text-red-600 bg-red-50 border border-red-200 rounded-lg px-3 py-2"
              role="alert"
            >
              {validationError}
            </div>
          )}

          {/* Progress Bar (during scan) */}
          {isScanning && (
            <div className="space-y-1">
              <div className="flex items-center justify-between text-xs text-slate-600">
                <span>Scanning...</span>
                <span>{completedPoints} / {totalPoints} points</span>
              </div>
              <div
                className="w-full bg-slate-200 rounded-full h-2 overflow-hidden"
                role="progressbar"
                aria-valuenow={completedPoints}
                aria-valuemin={0}
                aria-valuemax={totalPoints}
                aria-label={`PES scan progress: ${completedPoints} of ${totalPoints} points`}
              >
                <div
                  className="bg-blue-500 h-full rounded-full transition-all duration-300"
                  style={{ width: `${progressPercent}%` }}
                />
              </div>
            </div>
          )}

          {/* Error Display */}
          {pesState.error && !isScanning && (
            <div
              className="text-xs text-red-600 bg-red-50 border border-red-200 rounded-lg px-3 py-2"
              role="alert"
            >
              {pesState.error}
            </div>
          )}

          {/* Completion Summary */}
          {!isScanning && pesState.results.length > 0 && !pesState.error && (
            <div className="text-xs text-green-700 bg-green-50 border border-green-200 rounded-lg px-3 py-2">
              <div className="flex items-center gap-1.5">
                <svg className="w-3.5 h-3.5 flex-shrink-0" fill="currentColor" viewBox="0 0 20 20">
                  <path
                    fillRule="evenodd"
                    d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z"
                    clipRule="evenodd"
                  />
                </svg>
                <span>
                  Scan complete: {pesState.results.length} points
                  {pesState.equilibrium && (
                    <>, R<sub>eq</sub> = {pesState.equilibrium.r_bohr.toFixed(3)} bohr</>
                  )}
                </span>
              </div>
            </div>
          )}

          {/* Action Buttons */}
          <div className="flex gap-2">
            {!isScanning ? (
              <>
                <button
                  type="button"
                  onClick={handleStartScan}
                  disabled={!isValid || controlsDisabled}
                  className="flex-1 px-3 py-2 rounded-lg text-sm font-medium bg-blue-600 text-white hover:bg-blue-700 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
                  aria-label="Start PES scan"
                >
                  Start Scan
                </button>
                {pesState.results.length > 0 && (
                  <button
                    type="button"
                    onClick={handleReset}
                    disabled={controlsDisabled}
                    className="px-3 py-2 rounded-lg text-sm font-medium bg-slate-100 text-slate-700 hover:bg-slate-200 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-slate-500 disabled:opacity-50 disabled:cursor-not-allowed"
                    aria-label="Reset PES scan results"
                  >
                    Reset
                  </button>
                )}
              </>
            ) : (
              <button
                type="button"
                onClick={cancelScan}
                className="flex-1 px-3 py-2 rounded-lg text-sm font-medium bg-red-100 text-red-700 hover:bg-red-200 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-red-500"
                aria-label="Cancel PES scan"
              >
                <span className="flex items-center justify-center gap-2">
                  <span className="animate-spin rounded-full h-3.5 w-3.5 border-b-2 border-red-600" />
                  Cancel Scan
                </span>
              </button>
            )}
          </div>
        </div>
      )}
    </div>
  );
}

export type { PesScanPanelProps };
