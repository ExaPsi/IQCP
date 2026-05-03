/**
 * ModuleE - Module B: Integral Inspector page.
 *
 * Interactive exploration of one-electron integral matrices (S, T, V, H^core)
 * and step-by-step Fock matrix construction visualization for 17 preset
 * molecules with 6 basis sets, plus custom geometry input.
 *
 * Select a molecule and basis set, then inspect the integral matrices as
 * heatmaps with drill-down to individual element values. Clicking a cell
 * shows the primitive-pair breakdown (US-056).
 *
 * The Fock Build Tracer (US-058) shows how the HF Fock matrix F = H^core + G(P)
 * is constructed from its physical components, with J/K decomposition and
 * element-level drill-down. The ERI Browser (US-059) provides quartet-level
 * inspection of electron repulsion integrals.
 *
 * Code-split via React.lazy() to avoid loading integral-specific
 * components on other pages.
 *
 * @module pages/ModuleE
 * @see US-055 Integral Matrix Heatmap
 * @see US-056 Primitive Breakdown Panel
 * @see US-058 Fock Build Tracing
 * @see US-059 ERI Quartet Browser
 */

import { useMemo, useEffect, useRef, useCallback, useState } from 'react';
import {
  IntegralMatrixHeatmap,
  MatrixDetailPanel,
  MatrixTabBar,
  FockBuildTracer,
  EriBrowser,
} from '../components/integrals';
import { useIntegralStore } from '../stores/integralStore';
import type { MatrixTab, SelectedCell } from '../stores/integralStore';
import { useIntegralMatrices, PRESET_GEOMETRIES } from '../hooks/useIntegralMatrices';
import { useIntegralBreakdown } from '../hooks/useIntegralBreakdown';
import { useFockDecomposition } from '../hooks/useFockDecomposition';
import { useWorker } from '../hooks/useWorker';
import { getStateFromURL, updateURL, hasStateInURL } from '../lib/url';
import { APP_VERSION } from '../types/run-state';
import type { RunStateV1, IntegralParams } from '../types/run-state';
import type { GeometryInput, BasisSetName, KsScfResult } from '../worker/protocol';
import { ELEMENT_TO_Z } from '../lib/elements';

// ============================================================================
// Constants
// ============================================================================

/** Debounce delay for URL updates (ms) */
const URL_DEBOUNCE_MS = 300;

/** Molecule + basis set selectors shared with the SCF Sandbox (Module E, ModuleC.tsx) */
import { MOLECULE_PRESETS, BASIS_SET_OPTIONS, getPresetIdForMoleculeBasis } from '../stores/scfStore';

/** Available basis sets for custom geometry */
const BASIS_OPTIONS: { value: BasisSetName; label: string }[] = [
  { value: 'sto-3g', label: 'STO-3G' },
  { value: '3-21g', label: '3-21G' },
  { value: '6-31g', label: '6-31G' },
  { value: '6-31g*', label: '6-31G*' },
  { value: '6-31+g*', label: '6-31+G*' },
  { value: 'cc-pvdz', label: 'cc-pVDZ' },
];

// ============================================================================
// Deep Link Helpers
// ============================================================================

/**
 * Simple debounce utility for URL updates.
 */
function debounce<T extends (...args: Parameters<T>) => void>(
  fn: T,
  delay: number,
): T & { cancel: () => void } {
  let timeoutId: ReturnType<typeof setTimeout> | null = null;

  const debounced = ((...args: Parameters<T>) => {
    if (timeoutId) clearTimeout(timeoutId);
    timeoutId = setTimeout(() => {
      fn(...args);
      timeoutId = null;
    }, delay);
  }) as T & { cancel: () => void };

  debounced.cancel = () => {
    if (timeoutId) {
      clearTimeout(timeoutId);
      timeoutId = null;
    }
  };

  return debounced;
}

/**
 * Build a RunStateV1 object from the current Module B (Integral Inspector) store state.
 *
 * Serializes all user-visible state into the URL-encodable IntegralParams,
 * including Fock Build Tracer step/element and ERI Browser quartet selection.
 *
 * @see US-060 Module B Deep Links
 */
function buildRunState(
  systemId: string,
  activeTab: MatrixTab,
  selectedCell: SelectedCell | null,
  inputMode: 'preset' | 'custom',
  customGeometry: GeometryInput | null,
  customBasis: BasisSetName,
  fockBuildStep: 1 | 2 | 3,
  fockSelectedElement: [number, number] | null,
  eriSelectedQuartet: [number, number, number, number] | null,
): RunStateV1 {
  const integralParams: IntegralParams = {
    system_id: systemId,
    tab: activeTab,
  };

  if (selectedCell) {
    integralParams.cell = [selectedCell.row, selectedCell.col];
  }

  if (inputMode === 'custom' && customGeometry) {
    integralParams.custom = {
      atoms: customGeometry.atoms.map((a) => ({
        symbol: a.symbol,
        xyz: a.xyz,
      })),
      units: customGeometry.units,
      basis: customBasis,
    };
  }

  // Fock Build Tracing state (US-060)
  // Only include non-default fock_step to keep URLs compact
  if (fockBuildStep !== 1) {
    integralParams.fock_step = fockBuildStep;
  }
  if (fockSelectedElement) {
    integralParams.fock_element = fockSelectedElement;
  }

  // ERI Browser state (US-060)
  if (eriSelectedQuartet) {
    integralParams.eri_quartet = eriSelectedQuartet;
  }

  return {
    schema_version: 'run_state_v1',
    app_version: APP_VERSION,
    module: 'integrals',
    integrals: integralParams,
    ui: { mode: 'explain' },
  };
}

// ============================================================================
// Geometry Parser
// ============================================================================

/**
 * Parse XYZ text into a GeometryInput.
 * Returns an object with the parsed geometry (or null) and an optional
 * validation error message for unsupported elements.
 */
function parseGeometryText(
  text: string,
  units: 'bohr' | 'angstrom',
): { geometry: GeometryInput | null; error: string | null } {
  const lines = text.trim().split('\n').filter((line) => line.trim().length > 0);
  if (lines.length === 0) return { geometry: null, error: null };

  const atoms: GeometryInput['atoms'] = [];
  for (const line of lines) {
    const parts = line.trim().split(/\s+/);
    if (parts.length < 4) return { geometry: null, error: null };

    const symbol = parts[0];
    const x = parseFloat(parts[1]);
    const y = parseFloat(parts[2]);
    const z = parseFloat(parts[3]);

    if (!symbol || isNaN(x) || isNaN(y) || isNaN(z)) {
      return { geometry: null, error: null };
    }

    // Check if element is supported (H-Ar, Z=1-18)
    if (!(symbol in ELEMENT_TO_Z)) {
      return {
        geometry: null,
        error: `Unsupported element '${symbol}'. Only H-Ar (Z=1-18) are supported.`,
      };
    }

    atoms.push({ symbol, xyz: [x, y, z] });
  }

  return { geometry: { atoms, units }, error: null };
}

// ============================================================================
// Component
// ============================================================================

function ModuleE() {
  const { isReady, send } = useWorker();

  // Store selectors
  const systemId = useIntegralStore((s) => s.systemId);
  const inputMode = useIntegralStore((s) => s.inputMode);
  const customGeometryText = useIntegralStore((s) => s.customGeometryText);
  const customGeometry = useIntegralStore((s) => s.customGeometry);
  const customBasis = useIntegralStore((s) => s.customBasis);
  const activeTab = useIntegralStore((s) => s.activeTab);
  const selectedCell = useIntegralStore((s) => s.selectedCell);
  const compute = useIntegralStore((s) => s.compute);
  const eriCache = useIntegralStore((s) => s.eriCache);
  const showAllPrimitives = useIntegralStore((s) => s.showAllPrimitives);

  // Store actions
  const setSystemId = useIntegralStore((s) => s.setSystemId);
  const setInputMode = useIntegralStore((s) => s.setInputMode);
  const setCustomGeometryText = useIntegralStore((s) => s.setCustomGeometryText);
  const setCustomGeometry = useIntegralStore((s) => s.setCustomGeometry);
  const setCustomBasis = useIntegralStore((s) => s.setCustomBasis);
  const setActiveTab = useIntegralStore((s) => s.setActiveTab);
  const selectCell = useIntegralStore((s) => s.selectCell);
  const toggleShowAllPrimitives = useIntegralStore((s) => s.toggleShowAllPrimitives);

  // Fock Build Tracing store selectors (US-058)
  const fockBuildStep = useIntegralStore((s) => s.fockBuildStep);
  const fockStep2View = useIntegralStore((s) => s.fockStep2View);
  const fockSelectedElement = useIntegralStore((s) => s.fockSelectedElement);
  const fockDecomposition = useIntegralStore((s) => s.fockDecomposition);

  // Fock Build Tracing store actions (US-058)
  const setFockBuildStep = useIntegralStore((s) => s.setFockBuildStep);
  const setFockStep2View = useIntegralStore((s) => s.setFockStep2View);
  const setFockSelectedElement = useIntegralStore((s) => s.setFockSelectedElement);

  // ERI Browser store selectors (US-059, deep-linked via US-060)
  const eriSelectedQuartet = useIntegralStore((s) => s.eriSelectedQuartet);
  const setEriSelectedQuartet = useIntegralStore((s) => s.setEriSelectedQuartet);

  // -----------------------------------------------------------------------
  // UI Mode vs Store Mode
  // -----------------------------------------------------------------------
  // `uiMode` (local state) controls which UI elements are visible:
  //   - 'preset': shows molecule + basis set dropdowns
  //   - 'custom': shows XYZ textarea + units + basis selector
  //
  // `inputMode` (integral store) controls how useIntegralMatrices resolves
  // the geometry for WASM computation:
  //   - 'preset': uses PRESET_GEOMETRIES[systemId]
  //   - 'custom': uses customGeometry + customBasis from store
  //
  // They can diverge when a non-preset molecule+basis combo (e.g., CH4/3-21G)
  // is selected in preset UI mode. The store switches to 'custom' transparently
  // to enable on-the-fly computation, while the UI stays on the molecule/basis
  // selectors (not the textarea).
  // -----------------------------------------------------------------------
  const [uiMode, setUiMode] = useState<'preset' | 'custom'>(inputMode);

  // Local state for preset molecule + basis selection
  const [selectedMolecule, setSelectedMolecule] = useState('h2');
  const [selectedBasis, setSelectedBasis] = useState<BasisSetName>('sto-3g');

  // Local state for custom geometry units
  const [customUnits, setCustomUnits] = useState<'bohr' | 'angstrom'>('bohr');

  // Local state for geometry validation error (e.g., unsupported element)
  const [geometryError, setGeometryError] = useState<string | null>(null);

  // Local state for the density matrix obtained from a quick SCF run (US-058)
  const [scfDensity, setScfDensity] = useState<number[] | null>(null);
  const [scfRunning, setScfRunning] = useState(false);
  const [scfError, setScfError] = useState<string | null>(null);
  const scfRequestRef = useRef(0);

  // When molecule or basis changes in preset UI mode, update the integral store.
  // For pre-computed presets (5 molecules × STO-3G), use the fast preset path.
  // For all other combos, set up custom geometry transparently in the store.
  useEffect(() => {
    if (uiMode !== 'preset') return;
    const mol = MOLECULE_PRESETS.find((m) => m.id === selectedMolecule);
    if (!mol) return;

    const presetId = getPresetIdForMoleculeBasis(selectedMolecule, selectedBasis);

    if (presetId) {
      // Pre-computed preset exists — use fast preset path
      setInputMode('preset');
      setSystemId(presetId);
    } else {
      // No preset: set the integral store to 'custom' mode with the
      // molecule's geometry. The UI still shows molecule+basis selectors.
      const geometry: GeometryInput = {
        atoms: mol.atoms.map((a) => ({ symbol: a.symbol, xyz: a.xyz as [number, number, number] })),
        units: 'bohr',
      };
      setInputMode('custom');
      setCustomGeometry(geometry);
      setCustomBasis(selectedBasis);
      const text = mol.atoms
        .map((a) => `${a.symbol}  ${a.xyz[0].toFixed(6)}  ${a.xyz[1].toFixed(6)}  ${a.xyz[2].toFixed(6)}`)
        .join('\n');
      setCustomGeometryText(text);
    }
  }, [selectedMolecule, selectedBasis, uiMode, setSystemId, setInputMode, setCustomGeometry, setCustomBasis, setCustomGeometryText]);

  // Kick off integral computation via the hook
  useIntegralMatrices();

  // =========================================================================
  // Resolve current geometry and basis for breakdown requests (US-056)
  // =========================================================================

  const currentGeometry = useMemo((): GeometryInput | null => {
    if (inputMode === 'preset') {
      const preset = PRESET_GEOMETRIES[systemId];
      return preset?.geometry ?? null;
    }
    return customGeometry;
  }, [inputMode, systemId, customGeometry]);

  const currentBasisName = useMemo((): string => {
    if (inputMode === 'preset') {
      const preset = PRESET_GEOMETRIES[systemId];
      return preset?.basisName ?? 'sto-3g';
    }
    return customBasis;
  }, [inputMode, systemId, customBasis]);

  // =========================================================================
  // Fock decomposition hook (US-058)
  // =========================================================================

  const fockDecomp = useFockDecomposition({
    geometry: currentGeometry,
    basisName: currentBasisName,
    densityMatrix: scfDensity,
    isReady,
  });

  // Reset density when system changes
  const prevSystemRef = useRef<string | null>(null);
  useEffect(() => {
    const key = inputMode === 'preset' ? systemId : `custom:${currentBasisName}`;
    if (prevSystemRef.current !== null && prevSystemRef.current !== key) {
      setScfDensity(null);
      setScfError(null);
    }
    prevSystemRef.current = key;
  }, [inputMode, systemId, currentBasisName]);

  /**
   * Run SCF on the current system to obtain a converged density matrix,
   * then immediately compute the Fock decomposition.
   * Uses ks_scf with method='rhf' for both preset and custom systems
   * to avoid the "unknown system" error with scf_run for custom geometries.
   */
  const handleComputeFock = useCallback(async () => {
    if (!isReady || !currentGeometry) return;

    const requestId = ++scfRequestRef.current;
    setScfRunning(true);
    setScfError(null);
    setScfDensity(null);

    try {
      // Convert geometry to [Z, x, y, z] atom format for ks_scf
      const atoms: [number, number, number, number][] = currentGeometry.atoms.map((a) => [
        ELEMENT_TO_Z[a.symbol] ?? 1,
        ...(a.xyz as [number, number, number]),
      ]);

      // Use ks_scf (RHF mode) which computes integrals on-the-fly —
      // works for any molecule+basis combination without preset registration.
      const ksRequest = {
        type: 'ks_scf' as const,
        atoms,
        basisName: currentBasisName,
        method: 'rhf',
        convergenceProfile: 'tight',
        maxIterations: 200,
        useDiis: true,
      };

      const result = await send<KsScfResult>(ksRequest as never);

      if (scfRequestRef.current !== requestId) return; // stale

      if (!result.converged) {
        setScfError('SCF did not converge. Try a different system or basis set.');
        setScfRunning(false);
        return;
      }

      if (!result.densityMatrix || result.densityMatrix.length === 0) {
        setScfError('SCF converged but density matrix not available.');
        setScfRunning(false);
        return;
      }

      setScfDensity(result.densityMatrix);
      setScfRunning(false);
    } catch (err) {
      if (scfRequestRef.current !== requestId) return;
      const message = err instanceof Error ? err.message : 'SCF computation failed';
      setScfError(message);
      setScfRunning(false);
    }
  }, [isReady, currentGeometry, currentBasisName, send]);

  // After SCF density becomes available, auto-trigger Fock decomposition
  const prevDensityRef = useRef<number[] | null>(null);
  useEffect(() => {
    if (scfDensity && scfDensity !== prevDensityRef.current) {
      prevDensityRef.current = scfDensity;
      // Small delay to ensure state has settled
      const timer = setTimeout(() => {
        fockDecomp.compute();
      }, 50);
      return () => clearTimeout(timer);
    }
  }, [scfDensity, fockDecomp]);

  // =========================================================================
  // Primitive breakdown hook (US-056)
  // =========================================================================

  const breakdown = useIntegralBreakdown({
    geometry: currentGeometry,
    basisName: currentBasisName,
    integralType: activeTab,
    indexI: selectedCell?.row ?? 0,
    indexJ: selectedCell?.col ?? 0,
    isReady,
    enabled: selectedCell !== null && compute.status === 'success',
  });

  // =========================================================================
  // Deep Link: Read URL state on mount
  // =========================================================================

  const initializedRef = useRef(false);

  useEffect(() => {
    if (initializedRef.current) return;
    initializedRef.current = true;

    if (!hasStateInURL()) return;

    const urlState = getStateFromURL();
    if (!urlState || urlState.module !== 'integrals' || !urlState.integrals) return;

    const ip = urlState.integrals;

    if (ip.tab) setActiveTab(ip.tab);

    if (ip.custom) {
      setInputMode('custom');
      const geom: GeometryInput = {
        atoms: ip.custom.atoms,
        units: ip.custom.units,
      };
      setCustomGeometry(geom);
      setCustomBasis(ip.custom.basis as BasisSetName);
      setCustomUnits(ip.custom.units);
      // Build text representation
      const text = ip.custom.atoms
        .map((a) => `${a.symbol}  ${a.xyz[0].toFixed(6)}  ${a.xyz[1].toFixed(6)}  ${a.xyz[2].toFixed(6)}`)
        .join('\n');
      setCustomGeometryText(text);
    } else if (ip.system_id) {
      setSystemId(ip.system_id);
    }

    if (ip.cell) {
      selectCell({ row: ip.cell[0], col: ip.cell[1] });
    }

    // Fock Build Tracing state (US-060)
    if (ip.fock_step) {
      setFockBuildStep(ip.fock_step);
    }
    if (ip.fock_element) {
      setFockSelectedElement(ip.fock_element);
    }

    // ERI Browser state (US-060)
    if (ip.eri_quartet) {
      setEriSelectedQuartet(ip.eri_quartet);
    }
  }, [setActiveTab, setInputMode, setCustomGeometry, setCustomBasis, setCustomGeometryText, setSystemId, selectCell, setFockBuildStep, setFockSelectedElement, setEriSelectedQuartet]);

  // =========================================================================
  // Deep Link: Write URL state on changes (debounced)
  // =========================================================================

  const debouncedUpdateURL = useMemo(
    () => debounce((state: RunStateV1) => updateURL(state), URL_DEBOUNCE_MS),
    [],
  );

  useEffect(() => {
    if (!initializedRef.current) return;

    const state = buildRunState(
      systemId,
      activeTab,
      selectedCell,
      inputMode,
      customGeometry,
      customBasis,
      fockBuildStep,
      fockSelectedElement,
      eriSelectedQuartet,
    );
    debouncedUpdateURL(state);

    return () => debouncedUpdateURL.cancel();
  }, [systemId, activeTab, selectedCell, inputMode, customGeometry, customBasis, fockBuildStep, fockSelectedElement, eriSelectedQuartet, debouncedUpdateURL]);

  // =========================================================================
  // Custom geometry parsing
  // =========================================================================

  const handleGeometryTextChange = useCallback(
    (text: string) => {
      setCustomGeometryText(text);
      const { geometry, error } = parseGeometryText(text, customUnits);
      setCustomGeometry(geometry);
      setGeometryError(error);
    },
    [setCustomGeometryText, setCustomGeometry, customUnits],
  );

  const handleUnitsChange = useCallback(
    (units: 'bohr' | 'angstrom') => {
      setCustomUnits(units);
      // Re-parse with new units
      const { geometry, error } = parseGeometryText(customGeometryText, units);
      setCustomGeometry(geometry);
      setGeometryError(error);
    },
    [customGeometryText, setCustomGeometry],
  );

  // =========================================================================
  // Get the active matrix data
  // =========================================================================

  const activeMatrixData = useMemo(() => {
    if (compute.status !== 'success') return null;
    const data = compute.data;
    switch (activeTab) {
      case 'S': return data.sMatrix;
      case 'T': return data.tMatrix;
      case 'V': return data.vMatrix;
      case 'Hcore': return data.hCore;
    }
  }, [compute, activeTab]);

  // =========================================================================
  // Get detail panel data for selected cell
  // =========================================================================

  const detailData = useMemo(() => {
    if (!selectedCell || compute.status !== 'success' || !activeMatrixData) return null;
    const { row, col } = selectedCell;
    const data = compute.data;
    if (row >= data.nbf || col >= data.nbf) return null;

    return {
      row,
      col,
      rowLabel: data.labels[row] ?? `BF ${row + 1}`,
      colLabel: data.labels[col] ?? `BF ${col + 1}`,
      value: activeMatrixData[row * data.nbf + col],
    };
  }, [selectedCell, compute, activeMatrixData]);

  // =========================================================================
  // Advanced section state (Fock Build + ERI Browser)
  // =========================================================================

  const [advancedSection, setAdvancedSection] = useState<'none' | 'fock' | 'eri'>('none');

  const handleAdvancedToggle = useCallback((section: 'fock' | 'eri') => {
    setAdvancedSection((prev) => (prev === section ? 'none' : section));
  }, []);

  // Collapsible educational section state
  const [educationOpen, setEducationOpen] = useState(false);

  // =========================================================================
  // Matrix description for inline subtitle
  // =========================================================================

  const matrixDescription = useMemo((): React.ReactNode => {
    switch (activeTab) {
      case 'S': return 'Overlap matrix: measures spatial overlap between basis functions';
      case 'T': return 'Kinetic energy matrix: electron kinetic energy integrals';
      case 'V': return 'Nuclear attraction matrix: electron-nuclear Coulomb integrals';
      case 'Hcore': return <>Core Hamiltonian H<sup>core</sup> (T + V): one-electron part of the Fock matrix</>;
    }
  }, [activeTab]);

  // =========================================================================
  // Render
  // =========================================================================

  return (
    <div className="max-w-6xl mx-auto">
      {/* Header */}
      <div className="mb-6">
        <h1 className="text-3xl font-bold text-slate-800 mb-2">
          Module B: Integral Inspector
        </h1>
        <p className="text-slate-600">
          Visualize one-electron integral matrices (S, T, V, H<sup>core</sup>) as
          interactive heatmaps for 17 preset molecules with 6 basis sets. Click any
          matrix element to see its primitive-pair decomposition, trace the Fock matrix
          build step-by-step, and browse individual electron repulsion integrals.
        </p>
      </div>

      {/* Main content: sidebar + visualization */}
      <div className="grid grid-cols-1 lg:grid-cols-4 gap-6 mb-6">
        {/* ================================================================
            LEFT SIDEBAR: Compact controls
            ================================================================ */}
        <div className="lg:col-span-1 space-y-4">
          {/* System selector */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
            <h3 className="text-xs font-semibold text-slate-400 uppercase tracking-wider mb-3">
              System
            </h3>

            {/* Input mode toggle */}
            <div className="flex gap-2 mb-3">
              <button
                className={`flex-1 px-3 py-1.5 text-xs rounded-lg border transition-colors
                  ${uiMode === 'preset'
                    ? 'bg-primary-100 text-primary-800 border-primary-300'
                    : 'text-slate-600 border-slate-200 hover:bg-slate-50'
                  }`}
                onClick={() => { setUiMode('preset'); setInputMode('preset'); }}
              >
                Preset
              </button>
              <button
                className={`flex-1 px-3 py-1.5 text-xs rounded-lg border transition-colors
                  ${uiMode === 'custom'
                    ? 'bg-primary-100 text-primary-800 border-primary-300'
                    : 'text-slate-600 border-slate-200 hover:bg-slate-50'
                  }`}
                onClick={() => { setUiMode('custom'); setInputMode('custom'); }}
              >
                Custom
              </button>
            </div>

            {/* Preset molecule + basis selectors */}
            {uiMode === 'preset' && (
              <div className="space-y-3">
                <div>
                  <label htmlFor="mol-select" className="flex items-center gap-1 text-xs font-medium text-slate-600 mb-1">
                    Molecule
                    <span
                      className="inline-flex items-center justify-center w-3.5 h-3.5 rounded-full bg-slate-200 text-slate-500 text-[9px] font-bold cursor-help"
                      title="Pre-computed integrals are available for H2, HeH+, LiH, H2O, and NH3 with STO-3G. All other molecule/basis combinations compute integrals on-the-fly in your browser."
                    >
                      i
                    </span>
                  </label>
                  <select
                    id="mol-select"
                    value={selectedMolecule}
                    onChange={(e) => setSelectedMolecule(e.target.value)}
                    className="w-full px-3 py-2 text-sm border border-slate-300 rounded-lg bg-white focus:outline-none focus:ring-2 focus:ring-primary-500 focus:border-transparent"
                    aria-label="Select molecule"
                  >
                    {MOLECULE_PRESETS.map((mol) => (
                      <option key={mol.id} value={mol.id}>
                        {mol.label}
                      </option>
                    ))}
                  </select>
                </div>
                <div>
                  <label htmlFor="basis-select" className="block text-xs font-medium text-slate-600 mb-1">
                    Basis Set
                  </label>
                  <select
                    id="basis-select"
                    value={selectedBasis}
                    onChange={(e) => setSelectedBasis(e.target.value as BasisSetName)}
                    className="w-full px-3 py-2 text-sm border border-slate-300 rounded-lg bg-white focus:outline-none focus:ring-2 focus:ring-primary-500 focus:border-transparent"
                    aria-label="Select basis set"
                  >
                    {BASIS_SET_OPTIONS.map((b) => (
                      <option key={b.id} value={b.id}>
                        {b.label}
                      </option>
                    ))}
                  </select>
                </div>
                <p className="text-xs text-slate-500">
                  {MOLECULE_PRESETS.find((m) => m.id === selectedMolecule)?.description}
                  {' '}({MOLECULE_PRESETS.find((m) => m.id === selectedMolecule)?.nelec} electrons)
                  {!getPresetIdForMoleculeBasis(selectedMolecule, selectedBasis) && (
                    <span className="text-amber-600 block mt-0.5">On-the-fly integrals (may be slower)</span>
                  )}
                </p>
              </div>
            )}

            {/* Custom geometry input */}
            {uiMode === 'custom' && (
              <div className="space-y-3">
                <div>
                  <label
                    htmlFor="geometry-text"
                    className="block text-xs font-medium text-slate-600 mb-1"
                  >
                    Geometry (XYZ)
                  </label>
                  <textarea
                    id="geometry-text"
                    value={customGeometryText}
                    onChange={(e) => handleGeometryTextChange(e.target.value)}
                    rows={Math.min(Math.max(customGeometryText.split('\n').length + 1, 3), 12)}
                    className="w-full px-2 py-1.5 text-[11px] leading-tight font-mono border border-slate-300 rounded-lg bg-white focus:outline-none focus:ring-2 focus:ring-primary-500 focus:border-transparent resize-y"
                    placeholder={`H  0.000  0.000  0.000\nH  0.000  0.000  1.398`}
                    aria-label="Molecular geometry in XYZ format"
                  />
                  {customGeometryText.trim() && !customGeometry && (
                    <p className="text-xs text-red-600 mt-1">
                      {geometryError
                        ? geometryError
                        : 'Invalid geometry. Use: Symbol X Y Z (one atom per line)'}
                    </p>
                  )}
                </div>

                <div className="flex gap-3">
                  <div className="flex-1">
                    <label
                      htmlFor="units-select"
                      className="block text-xs font-medium text-slate-600 mb-1"
                    >
                      Units
                    </label>
                    <select
                      id="units-select"
                      value={customUnits}
                      onChange={(e) => handleUnitsChange(e.target.value as 'bohr' | 'angstrom')}
                      className="w-full px-2 py-1.5 text-xs border border-slate-300 rounded-lg bg-white focus:outline-none focus:ring-2 focus:ring-primary-500"
                    >
                      <option value="bohr">Bohr</option>
                      <option value="angstrom">Angstrom</option>
                    </select>
                  </div>

                  <div className="flex-1">
                    <label
                      htmlFor="basis-select"
                      className="block text-xs font-medium text-slate-600 mb-1"
                    >
                      Basis Set
                    </label>
                    <select
                      id="basis-select"
                      value={customBasis}
                      onChange={(e) => setCustomBasis(e.target.value as BasisSetName)}
                      className="w-full px-2 py-1.5 text-xs border border-slate-300 rounded-lg bg-white focus:outline-none focus:ring-2 focus:ring-primary-500"
                    >
                      {BASIS_OPTIONS.map((b) => (
                        <option key={b.value} value={b.value}>
                          {b.label}
                        </option>
                      ))}
                    </select>
                  </div>
                </div>
              </div>
            )}
          </div>

          {/* Basis / system info panel */}
          {compute.status === 'success' && (
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
              <h3 className="text-xs font-semibold text-slate-400 uppercase tracking-wider mb-2">
                Basis Info
              </h3>
              <dl className="text-xs space-y-1.5">
                <div className="flex justify-between">
                  <dt className="text-slate-500">Basis Functions</dt>
                  <dd className="text-slate-800 font-medium">{compute.data.nbf}</dd>
                </div>
                <div className="flex justify-between">
                  <dt className="text-slate-500">Matrix Size</dt>
                  <dd className="text-slate-800 font-medium">
                    {compute.data.nbf} x {compute.data.nbf}
                  </dd>
                </div>
                <div className="flex justify-between">
                  <dt className="text-slate-500">V_nn</dt>
                  <dd className="text-slate-800 font-mono text-[11px]">
                    {compute.data.nuclearRepulsion.toFixed(8)} Ha
                  </dd>
                </div>
                <div className="flex justify-between">
                  <dt className="text-slate-500">Compute</dt>
                  <dd className="text-slate-800 font-medium">
                    {compute.data.computeTimeMs.toFixed(1)} ms
                  </dd>
                </div>
                {/* ERI cache status (US-057 AC3) */}
                <div className="flex justify-between items-center">
                  <dt className="text-slate-500">ERIs</dt>
                  <dd className="text-slate-800 font-medium">
                    {eriCache.status === 'computing' && (
                      <span className="flex items-center gap-1.5 text-amber-600">
                        <span className="inline-block w-2 h-2 rounded-full bg-amber-500 animate-pulse" />
                        computing...
                      </span>
                    )}
                    {eriCache.status === 'success' && (
                      <span className="text-emerald-700">
                        {eriCache.eriCount} unique
                      </span>
                    )}
                    {eriCache.status === 'error' && (
                      <span className="text-red-600" title={eriCache.error}>failed</span>
                    )}
                    {eriCache.status === 'idle' && (
                      <span className="text-slate-400">pending</span>
                    )}
                  </dd>
                </div>
              </dl>
            </div>
          )}

          {/* Cell detail (when a heatmap cell is selected) */}
          {detailData && (
            <div className="bg-white rounded-xl shadow-sm border border-primary-200 p-4">
              <h3 className="text-xs font-semibold text-primary-600 uppercase tracking-wider mb-2">
                Selected Cell
              </h3>
              <dl className="text-xs space-y-1">
                <div className="flex justify-between">
                  <dt className="text-slate-500">Element</dt>
                  <dd className="text-slate-800 font-medium">
                    {activeTab}({detailData.rowLabel}, {detailData.colLabel})
                  </dd>
                </div>
                <div className="flex justify-between">
                  <dt className="text-slate-500">Index</dt>
                  <dd className="text-slate-800 font-mono text-[11px]">
                    ({detailData.row}, {detailData.col})
                  </dd>
                </div>
                <div className="flex justify-between">
                  <dt className="text-slate-500">Value</dt>
                  <dd className="text-slate-800 font-mono text-[11px]">
                    {detailData.value.toExponential(6)}
                  </dd>
                </div>
              </dl>
            </div>
          )}
        </div>

        {/* ================================================================
            MAIN CONTENT: Matrix tabs + heatmap + breakdown
            ================================================================ */}
        <div className="lg:col-span-3 space-y-4">
          {/* Prominent horizontal matrix tab bar */}
          <MatrixTabBar activeTab={activeTab} onTabChange={setActiveTab} />

          {/* Matrix description subtitle */}
          <p className="text-sm text-slate-500 -mt-2 px-1">
            {matrixDescription}
          </p>

          {/* Loading state */}
          {(compute.status === 'idle' || compute.status === 'computing') && (
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 flex items-center justify-center" style={{ minHeight: '400px' }}>
              <div className="flex items-center gap-3">
                {compute.status === 'computing' ? (
                  <>
                    <div className="animate-spin rounded-full h-6 w-6 border-b-2 border-primary-600" />
                    <span className="text-slate-600 text-sm">Computing integral matrices...</span>
                  </>
                ) : !isReady ? (
                  <>
                    <div className="animate-spin rounded-full h-6 w-6 border-b-2 border-slate-400" />
                    <span className="text-slate-500 text-sm">Initializing WASM module...</span>
                  </>
                ) : (
                  <span className="text-slate-500 text-sm">Select a molecule to compute integrals.</span>
                )}
              </div>
            </div>
          )}

          {/* Error state */}
          {compute.status === 'error' && (
            <div className="bg-white rounded-xl shadow-sm border border-red-200 p-6">
              <div className="flex items-start gap-3">
                <div className="flex-shrink-0 w-8 h-8 bg-red-100 rounded-full flex items-center justify-center">
                  <span className="text-red-600 text-sm font-bold">!</span>
                </div>
                <div>
                  <h3 className="text-sm font-semibold text-red-800">Computation Error</h3>
                  <p className="text-sm text-red-700 mt-1">{compute.error}</p>
                </div>
              </div>
            </div>
          )}

          {/* Heatmap */}
          {compute.status === 'success' && activeMatrixData && (
            <div
              role="tabpanel"
              id={`matrix-panel-${activeTab}`}
              aria-labelledby={`matrix-tab-${activeTab}`}
            >
              <IntegralMatrixHeatmap
                matrix={activeMatrixData}
                nbf={compute.data.nbf}
                labels={compute.data.labels}
                matrixType={activeTab}
                selectedCell={selectedCell}
                onCellSelect={selectCell}
                minHeight={420}
              />
            </div>
          )}

          {/* Detail panel with primitive breakdown (US-056) */}
          {detailData && (
            <MatrixDetailPanel
              row={detailData.row}
              col={detailData.col}
              rowLabel={detailData.rowLabel}
              colLabel={detailData.colLabel}
              value={detailData.value}
              matrixType={activeTab}
              breakdownData={breakdown.data}
              breakdownLoading={breakdown.loading}
              breakdownError={breakdown.error}
              showAllPrimitives={showAllPrimitives}
              onToggleShowAll={toggleShowAllPrimitives}
            />
          )}

          {/* =============================================================
              ADVANCED FEATURES: Fock Build + ERI Browser
              ============================================================= */}
          {compute.status === 'success' && (
            <div className="space-y-4">
              {/* Advanced toggle bar */}
              <div className="flex items-center gap-3">
                <span className="text-xs font-medium text-slate-400 uppercase tracking-wider">Advanced:</span>
                <button
                  type="button"
                  onClick={() => handleAdvancedToggle('fock')}
                  className={`
                    px-3 py-1.5 text-xs rounded-lg border transition-colors
                    focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
                    ${advancedSection === 'fock'
                      ? 'bg-primary-50 border-primary-300 text-primary-700 font-medium'
                      : 'bg-white border-slate-200 text-slate-600 hover:bg-slate-50 hover:border-slate-300'
                    }
                  `}
                  aria-pressed={advancedSection === 'fock'}
                >
                  Fock Build Tracer
                </button>
                <button
                  type="button"
                  onClick={() => handleAdvancedToggle('eri')}
                  disabled={eriCache.status !== 'success'}
                  className={`
                    px-3 py-1.5 text-xs rounded-lg border transition-colors
                    focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
                    ${advancedSection === 'eri'
                      ? 'bg-primary-50 border-primary-300 text-primary-700 font-medium'
                      : eriCache.status !== 'success'
                        ? 'bg-slate-50 border-slate-200 text-slate-400 cursor-not-allowed'
                        : 'bg-white border-slate-200 text-slate-600 hover:bg-slate-50 hover:border-slate-300'
                    }
                  `}
                  aria-pressed={advancedSection === 'eri'}
                  title={eriCache.status !== 'success' ? 'ERI computation must complete first' : undefined}
                >
                  ERI Browser
                </button>
              </div>

              {/* -------------------------------------------------------
                  Fock Build Section (expanded)
                  ------------------------------------------------------- */}
              {advancedSection === 'fock' && (
                <div className="space-y-4">
                  <p className="text-slate-600 text-sm">
                    Trace the step-by-step construction of the Hartree-Fock Fock matrix
                    F = H<sup>core</sup> + G(P), where G = J &minus; 0.5K (Coulomb minus
                    exchange). In Kohn-Sham DFT, the exact exchange K is replaced by
                    an exchange-correlation potential V<sub>xc</sub>, giving
                    F = H<sup>core</sup> + J + V<sub>xc</sub>. Requires an SCF
                    calculation to obtain the converged density matrix P.
                  </p>

                  {/* Action button: Compute Fock Decomposition */}
                  {fockDecomposition.status === 'idle' && !scfRunning && (
                    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5 text-center">
                      <p className="text-sm text-slate-600 mb-3">
                        Run a tight-convergence RHF calculation on the current system
                        to obtain the density matrix, then decompose the Fock matrix.
                      </p>
                      <button
                        onClick={handleComputeFock}
                        disabled={!isReady || !currentGeometry}
                        className="px-6 py-2.5 bg-primary-600 text-white text-sm font-medium rounded-lg
                          hover:bg-primary-700 disabled:bg-slate-300 disabled:cursor-not-allowed
                          transition-colors focus-visible:outline-none focus-visible:ring-2
                          focus-visible:ring-primary-500 focus-visible:ring-offset-2"
                      >
                        Run SCF + Decompose Fock Matrix
                      </button>
                      {inputMode === 'custom' && (
                        <p className="text-xs text-amber-600 mt-2">
                          Note: Custom geometries require the system to be registered first.
                          If SCF fails, try using a preset system.
                        </p>
                      )}
                    </div>
                  )}

                  {/* SCF running state */}
                  {scfRunning && (
                    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
                      <div className="flex items-center justify-center gap-3">
                        <div className="animate-spin rounded-full h-5 w-5 border-b-2 border-primary-600" />
                        <span className="text-sm text-slate-600">
                          Running SCF to obtain density matrix...
                        </span>
                      </div>
                    </div>
                  )}

                  {/* SCF error */}
                  {scfError && !scfRunning && (
                    <div className="bg-white rounded-xl shadow-sm border border-red-200 p-5">
                      <div className="flex items-start gap-3">
                        <div className="flex-shrink-0 w-7 h-7 bg-red-100 rounded-full flex items-center justify-center">
                          <span className="text-red-600 text-xs font-bold">!</span>
                        </div>
                        <div className="flex-1">
                          <h3 className="text-sm font-semibold text-red-800">SCF Error</h3>
                          <p className="text-sm text-red-700 mt-1">{scfError}</p>
                          <button
                            onClick={handleComputeFock}
                            className="mt-3 px-4 py-1.5 text-xs bg-red-100 text-red-800 rounded-lg
                              hover:bg-red-200 transition-colors"
                          >
                            Retry
                          </button>
                        </div>
                      </div>
                    </div>
                  )}

                  {/* Fock decomposition computing state */}
                  {fockDecomposition.status === 'computing' && (
                    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
                      <div className="flex items-center justify-center gap-3">
                        <div className="animate-spin rounded-full h-5 w-5 border-b-2 border-primary-600" />
                        <span className="text-sm text-slate-600">
                          Computing Fock decomposition (J, K, G, F)...
                        </span>
                      </div>
                    </div>
                  )}

                  {/* Fock decomposition error */}
                  {fockDecomposition.status === 'error' && (
                    <div className="bg-white rounded-xl shadow-sm border border-red-200 p-5">
                      <div className="flex items-start gap-3">
                        <div className="flex-shrink-0 w-7 h-7 bg-red-100 rounded-full flex items-center justify-center">
                          <span className="text-red-600 text-xs font-bold">!</span>
                        </div>
                        <div className="flex-1">
                          <h3 className="text-sm font-semibold text-red-800">Decomposition Error</h3>
                          <p className="text-sm text-red-700 mt-1">{fockDecomposition.error}</p>
                          <button
                            onClick={handleComputeFock}
                            className="mt-3 px-4 py-1.5 text-xs bg-red-100 text-red-800 rounded-lg
                              hover:bg-red-200 transition-colors"
                          >
                            Retry
                          </button>
                        </div>
                      </div>
                    </div>
                  )}

                  {/* Fock Build Tracer (shown when decomposition data is available) */}
                  {fockDecomposition.status === 'success' && (
                    <FockBuildTracer
                      decomposition={fockDecomposition.data}
                      step={fockBuildStep}
                      onStepChange={setFockBuildStep}
                      step2View={fockStep2View}
                      onStep2ViewChange={setFockStep2View}
                      selectedElement={fockSelectedElement}
                      onElementSelect={setFockSelectedElement}
                    />
                  )}
                </div>
              )}

              {/* -------------------------------------------------------
                  ERI Browser Section (expanded)
                  ------------------------------------------------------- */}
              {advancedSection === 'eri' && eriCache.status === 'success' && currentGeometry && (
                <EriBrowser
                  geometry={currentGeometry}
                  basisName={currentBasisName}
                  nbf={compute.data.nbf}
                  eriCache={eriCache}
                />
              )}
            </div>
          )}
        </div>
      </div>

      {/* ================================================================
          EDUCATIONAL BACKGROUND (collapsible, default collapsed)
          ================================================================ */}
      <div className="bg-slate-100 rounded-xl overflow-hidden">
        <button
          type="button"
          onClick={() => setEducationOpen(!educationOpen)}
          className="w-full flex items-center justify-between p-5 text-left focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-inset focus-visible:ring-primary-500"
          aria-expanded={educationOpen}
          aria-controls="education-content-e"
        >
          <h2 className="font-semibold text-slate-700">
            Understanding One-Electron Integrals
          </h2>
          <svg
            className={`w-5 h-5 text-slate-400 transition-transform duration-200 ${educationOpen ? 'rotate-180' : ''}`}
            fill="none"
            viewBox="0 0 24 24"
            strokeWidth={2}
            stroke="currentColor"
            aria-hidden="true"
          >
            <path strokeLinecap="round" strokeLinejoin="round" d="M19.5 8.25l-7.5 7.5-7.5-7.5" />
          </svg>
        </button>

        {educationOpen && (
          <div id="education-content-e" className="px-5 pb-5">
            <div className="text-slate-600 text-sm space-y-3">
              <p>
                In both Hartree-Fock and Kohn-Sham DFT, the electronic Hamiltonian is
                decomposed into one-electron and two-electron contributions. This module
                lets you explore the four one-electron integral matrices that form the
                foundation of every SCF calculation, regardless of the method used.
              </p>

              <div className="grid md:grid-cols-2 gap-4 mt-4">
                <div className="bg-white rounded-lg p-4 border border-slate-200">
                  <h3 className="font-semibold text-slate-700 mb-2">
                    Overlap Matrix (S)
                  </h3>
                  <p className="text-xs text-slate-600">
                    S_ij measures the spatial overlap between basis functions i and j.
                    Diagonal elements are always 1.0 for normalized functions.
                    Off-diagonal elements range from 0 (orthogonal) to near 1 (nearly
                    identical spatial extent). The overlap matrix is needed because
                    Gaussian basis functions are not orthogonal.
                  </p>
                </div>

                <div className="bg-white rounded-lg p-4 border border-slate-200">
                  <h3 className="font-semibold text-slate-700 mb-2">
                    Kinetic Energy Matrix (T)
                  </h3>
                  <p className="text-xs text-slate-600">
                    T_ij represents the kinetic energy contribution from basis functions
                    i and j. All diagonal elements are positive (kinetic energy is always
                    positive). Functions with higher angular momentum (p, d) have larger
                    kinetic energy because they vary more rapidly in space.
                  </p>
                </div>

                <div className="bg-white rounded-lg p-4 border border-slate-200">
                  <h3 className="font-semibold text-slate-700 mb-2">
                    Nuclear Attraction Matrix (V)
                  </h3>
                  <p className="text-xs text-slate-600">
                    V_ij quantifies the electrostatic attraction between electrons
                    (described by basis functions) and all atomic nuclei. All elements
                    are negative because the electron-nuclear interaction is attractive.
                    Core orbitals near high-Z nuclei have the largest magnitudes.
                    These integrals use the Boys function from Module C.
                  </p>
                </div>

                <div className="bg-white rounded-lg p-4 border border-slate-200">
                  <h3 className="font-semibold text-slate-700 mb-2">
                    Core Hamiltonian (H<sup>core</sup>)
                  </h3>
                  <p className="text-xs text-slate-600">
                    H<sup>core</sup> = T + V is the one-electron part of the Fock matrix. It
                    represents the combined effect of kinetic energy and nuclear
                    attraction without electron-electron interactions. This is the
                    initial Fock matrix guess (iteration 0) used in Module E's SCF
                    procedure. In HF: F = H<sup>core</sup> + J &minus; 0.5K. In DFT:
                    F = H<sup>core</sup> + J + V<sub>xc</sub>.
                  </p>
                </div>
              </div>

              <p className="text-xs text-slate-500 mt-4">
                All matrices are symmetric (M_ij = M_ji) because the underlying quantum
                mechanical operators are Hermitian. Click any cell to see its value and
                primitive-pair decomposition. Use the Fock Build Tracer to see how
                H<sup>core</sup> combines with two-electron integrals, and the ERI Browser
                to inspect individual (ij|kl) electron repulsion integrals.
              </p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

export default ModuleE;
