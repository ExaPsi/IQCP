/**
 * ModuleD - Module A: Basis Explorer page.
 *
 * Interactive exploration of 6 Gaussian basis sets (STO-3G, 3-21G, 6-31G,
 * 6-31G*, 6-31+G*, cc-pVDZ) for all 18 supported elements (H-Ar).
 * Students can select an element and basis set to inspect the shell structure,
 * understanding how atomic orbitals are represented by contracted Gaussians.
 *
 * When a shell is selected, the radial profile plot shows R(r) for that
 * shell, and the shell data table shows the primitive parameters (exponents,
 * coefficients, normalization, effective coefficients).
 *
 * US-052: Basis Set Comparison Mode adds overlaid radial profiles across
 * multiple basis sets, with color-coded traces, split-valence visualization,
 * and guided discussion prompts.
 *
 * US-053: Exponent Sliders allow real-time modification of Gaussian exponents
 * with live radial profile and overlap metric feedback.
 *
 * US-054: Overlap vs. Distance Plot shows how the overlap integral S_ab
 * between two basis functions decays with interatomic distance, connecting
 * basis function shape to molecular bonding.
 *
 * @module pages/ModuleD
 */

import { useMemo, useEffect, useRef, useCallback, useState } from 'react';
import {
  ElementSelector,
  BasisSelector,
  ShellList,
  RadialProfilePlot,
  ShellDataTable,
  deriveShellDisplayInfo,
  BasisComparisonControls,
  BasisComparisonPlot,
  DiscussionPrompt,
  getAvailableShellTypes,
  ExponentSliders,
  OverlapDistanceControls,
  OverlapDistancePlot,
} from '../components/basisExplorer';
import { useBasisStore } from '../stores/basisStore';
import type { OverlapTraceConfig } from '../stores/basisStore';
import { useBasisInfo } from '../hooks/useBasisInfo';
import { useRadialProfile } from '../hooks/useRadialProfile';
import { useRadialComparison } from '../hooks/useRadialComparison';
import { useExponentSliders } from '../hooks/useExponentSliders';
import { useOverlapDistance } from '../hooks/useOverlapDistance';
import type { OverlapDistanceResult } from '../hooks/useOverlapDistance';
import { useWorker } from '../hooks/useWorker';
import { getStateFromURL, updateURL, hasStateInURL } from '../lib/url';
import { APP_VERSION } from '../types/run-state';
import type { RunStateV1, BasisParams } from '../types/run-state';

// ============================================================================
// Deep Link Helpers
// ============================================================================

/**
 * Debounce delay for URL updates (ms).
 */
const URL_DEBOUNCE_MS = 300;

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
 * Build a RunStateV1 object from the current Module A (Basis Explorer) store state.
 */
function buildRunState(
  element: number,
  basis: string,
  shell: number | null,
  comparisonMode: boolean,
  comparisonBases: string[],
  comparisonShellType: string,
  exponentOverrides: Record<number, number[]>,
  overlapDistanceMode: boolean,
  overlapTraces: OverlapTraceConfig[],
): RunStateV1 {
  const basisParams: BasisParams = {
    element,
    basis,
    shell,
  };

  if (comparisonMode && comparisonBases.length > 0) {
    basisParams.compare = true;
    basisParams.compare_bases = comparisonBases;
    if (comparisonShellType) {
      basisParams.compare_shell_type = comparisonShellType;
    }
  }

  // Include exponent overrides in deep link (only non-empty entries)
  const overrideEntries = Object.entries(exponentOverrides);
  if (overrideEntries.length > 0) {
    const serialized: Record<string, number[]> = {};
    for (const [key, values] of overrideEntries) {
      serialized[key] = values;
    }
    basisParams.exponent_overrides = serialized;
  }

  // Include overlap distance state in deep link (US-054)
  if (overlapDistanceMode) {
    basisParams.overlap_mode = true;
    if (overlapTraces.length > 0) {
      basisParams.overlap_traces = overlapTraces.map((t) => ({
        element_b: t.elementB,
        basis_b: t.basisB,
        shell_index_b: t.shellIndexB,
      }));
    }
  }

  return {
    schema_version: 'run_state_v1',
    app_version: APP_VERSION,
    module: 'basis',
    basis: basisParams,
    ui: { mode: 'explain' },
  };
}

// ============================================================================
// Component
// ============================================================================

/**
 * Module A: Basis Explorer page.
 *
 * Provides an interactive interface for exploring the 6 Gaussian basis sets
 * (STO-3G, 3-21G, 6-31G, 6-31G*, 6-31+G*, cc-pVDZ) for all 18 supported
 * elements (H-Ar). Students select an element from a mini periodic table
 * and a basis set, then inspect the shell structure and view radial profile plots.
 *
 * US-052 adds comparison mode: when activated, the same shell type is shown
 * across 2-5 basis sets with color-coded overlaid traces.
 */
function ModuleD() {
  // Store selectors
  const selectedElement = useBasisStore((state) => state.selectedElement);
  const selectedBasis = useBasisStore((state) => state.selectedBasis);
  const selectedShell = useBasisStore((state) => state.selectedShell);
  const radialProfileScale = useBasisStore((state) => state.radialProfileScale);
  const comparisonMode = useBasisStore((state) => state.comparisonMode);
  const comparisonBases = useBasisStore((state) => state.comparisonBases);
  const exponentOverrides = useBasisStore((state) => state.exponentOverrides);
  const overlapDistanceMode = useBasisStore((state) => state.overlapDistanceMode);
  const overlapTraces = useBasisStore((state) => state.overlapTraces);
  const overlapPlotScale = useBasisStore((state) => state.overlapPlotScale);

  // Store actions
  const setElement = useBasisStore((state) => state.setElement);
  const setBasis = useBasisStore((state) => state.setBasis);
  const setShell = useBasisStore((state) => state.setShell);
  const setRadialProfileScale = useBasisStore((state) => state.setRadialProfileScale);
  const toggleComparisonMode = useBasisStore((state) => state.toggleComparisonMode);
  const setComparisonMode = useBasisStore((state) => state.setComparisonMode);
  const setComparisonBases = useBasisStore((state) => state.setComparisonBases);
  const addComparisonBasis = useBasisStore((state) => state.addComparisonBasis);
  const removeComparisonBasis = useBasisStore((state) => state.removeComparisonBasis);
  const toggleOverlapDistanceMode = useBasisStore((state) => state.toggleOverlapDistanceMode);
  const setOverlapDistanceMode = useBasisStore((state) => state.setOverlapDistanceMode);
  const addOverlapTrace = useBasisStore((state) => state.addOverlapTrace);
  const removeOverlapTrace = useBasisStore((state) => state.removeOverlapTrace);
  const clearOverlapTraces = useBasisStore((state) => state.clearOverlapTraces);
  const setOverlapPlotScale = useBasisStore((state) => state.setOverlapPlotScale);

  const { isReady } = useWorker();

  // Local state for comparison shell type (not in global store to avoid coupling)
  const [comparisonShellType, setComparisonShellType] = useState('1s');

  // Fetch basis shell data from WASM (for ShellDataTable and comparison)
  const { shells: wasmShells } = useBasisInfo(selectedElement, selectedBasis);

  // Derive display labels from WASM shell data
  const shellDisplayInfo = useMemo(() => {
    if (!wasmShells) return null;
    return deriveShellDisplayInfo(wasmShells, selectedElement);
  }, [wasmShells, selectedElement]);

  // Available shell types for comparison (derived from the primary basis)
  const availableShellTypes = useMemo(() => {
    if (!shellDisplayInfo) return [];
    return getAvailableShellTypes(shellDisplayInfo);
  }, [shellDisplayInfo]);

  // Auto-select the comparison shell type based on the currently selected shell
  useEffect(() => {
    if (selectedShell !== null && shellDisplayInfo && shellDisplayInfo[selectedShell]) {
      const shellLabel = shellDisplayInfo[selectedShell].label.replace("'", '');
      setComparisonShellType(shellLabel);
    }
  }, [selectedShell, shellDisplayInfo]);

  // Ensure the primary basis is always in comparison bases when comparison is active
  const effectiveComparisonBases = useMemo(() => {
    if (!comparisonMode) return [];
    const bases = comparisonBases.includes(selectedBasis)
      ? comparisonBases
      : [selectedBasis, ...comparisonBases];
    // Ensure at least 2 bases (add sto-3g or 3-21g as second if only primary)
    if (bases.length < 2) {
      const fallback = selectedBasis === 'sto-3g' ? '3-21g' : 'sto-3g';
      return [...bases, fallback];
    }
    return bases;
  }, [comparisonMode, comparisonBases, selectedBasis]);

  // Fetch radial profile data when a shell is selected (single-shell view)
  const { data: profileData, loading: profileLoading } = useRadialProfile({
    atomicNumber: selectedElement,
    basisName: selectedBasis,
    shellIndex: selectedShell,
    isReady,
  });

  // Fetch comparison profiles across multiple basis sets
  const {
    profiles: comparisonProfiles,
    loading: comparisonLoading,
    error: comparisonError,
  } = useRadialComparison({
    atomicNumber: selectedElement,
    shellTypeLabel: comparisonShellType,
    selectedBases: effectiveComparisonBases,
    isActive: comparisonMode && selectedShell !== null,
  });

  // Exponent sliders hook (US-053)
  const {
    defaultExponents,
    currentExponents,
    modifiedProfile,
    overlapMetric,
    isModified: exponentsModified,
    handleExponentChange,
    handleReset: handleExponentReset,
  } = useExponentSliders({
    shellIndex: selectedShell,
    profileData,
  });

  // ===========================================================================
  // Overlap vs. Distance (US-054)
  // ===========================================================================

  // Primary overlap pair: uses current element/basis/shell for both A and B
  // (homonuclear self-overlap as default)
  const primaryOverlap = useOverlapDistance({
    elementA: selectedElement,
    basisNameA: selectedBasis,
    shellIndexA: selectedShell ?? 0,
    elementB: selectedElement,
    basisNameB: selectedBasis,
    shellIndexB: selectedShell ?? 0,
    isReady: isReady && overlapDistanceMode && selectedShell !== null,
  });

  // Additional overlay traces (up to 4 total)
  const overlay0 = useOverlapDistance({
    elementA: selectedElement,
    basisNameA: selectedBasis,
    shellIndexA: selectedShell ?? 0,
    elementB: overlapTraces[0]?.elementB ?? selectedElement,
    basisNameB: overlapTraces[0]?.basisB ?? selectedBasis,
    shellIndexB: overlapTraces[0]?.shellIndexB ?? 0,
    isReady: isReady && overlapDistanceMode && selectedShell !== null && overlapTraces.length > 0,
  });

  const overlay1 = useOverlapDistance({
    elementA: selectedElement,
    basisNameA: selectedBasis,
    shellIndexA: selectedShell ?? 0,
    elementB: overlapTraces[1]?.elementB ?? selectedElement,
    basisNameB: overlapTraces[1]?.basisB ?? selectedBasis,
    shellIndexB: overlapTraces[1]?.shellIndexB ?? 0,
    isReady: isReady && overlapDistanceMode && selectedShell !== null && overlapTraces.length > 1,
  });

  const overlay2 = useOverlapDistance({
    elementA: selectedElement,
    basisNameA: selectedBasis,
    shellIndexA: selectedShell ?? 0,
    elementB: overlapTraces[2]?.elementB ?? selectedElement,
    basisNameB: overlapTraces[2]?.basisB ?? selectedBasis,
    shellIndexB: overlapTraces[2]?.shellIndexB ?? 0,
    isReady: isReady && overlapDistanceMode && selectedShell !== null && overlapTraces.length > 2,
  });

  const overlay3 = useOverlapDistance({
    elementA: selectedElement,
    basisNameA: selectedBasis,
    shellIndexA: selectedShell ?? 0,
    elementB: overlapTraces[3]?.elementB ?? selectedElement,
    basisNameB: overlapTraces[3]?.basisB ?? selectedBasis,
    shellIndexB: overlapTraces[3]?.shellIndexB ?? 0,
    isReady: isReady && overlapDistanceMode && selectedShell !== null && overlapTraces.length > 3,
  });

  // Extract data from overlay hooks for stable dependency tracking
  const overlay0Data = overlay0.data;
  const overlay1Data = overlay1.data;
  const overlay2Data = overlay2.data;
  const overlay3Data = overlay3.data;

  // Collect all active overlap traces for the plot
  const overlapPlotTraces = useMemo((): OverlapDistanceResult[] => {
    const results: OverlapDistanceResult[] = [];
    if (primaryOverlap.data) results.push(primaryOverlap.data);

    const overlayData = [overlay0Data, overlay1Data, overlay2Data, overlay3Data];
    for (let i = 0; i < overlapTraces.length; i++) {
      if (overlayData[i]) results.push(overlayData[i]!);
    }
    return results;
  }, [
    primaryOverlap.data,
    overlay0Data,
    overlay1Data,
    overlay2Data,
    overlay3Data,
    overlapTraces.length,
  ]);

  const overlapLoading =
    primaryOverlap.loading ||
    overlay0.loading ||
    overlay1.loading ||
    overlay2.loading ||
    overlay3.loading;

  // Determine the currently selected shell data for ShellDataTable
  const selectedShellData = useMemo(() => {
    if (selectedShell === null || !wasmShells || !shellDisplayInfo) return null;
    const ws = wasmShells[selectedShell];
    const di = shellDisplayInfo[selectedShell];
    if (!ws || !di) return null;
    return {
      exponents: ws.exponents,
      coefficients: ws.coefficients,
      angularMomentum: ws.angularMomentum,
      angularMomentumLabel: ws.angularMomentumLabel,
      shellLabel: di.label,
    };
  }, [selectedShell, wasmShells, shellDisplayInfo]);

  // ===========================================================================
  // Deep Link: Read URL state on mount
  // ===========================================================================

  const initializedRef = useRef(false);

  useEffect(() => {
    if (initializedRef.current) return;
    initializedRef.current = true;

    if (!hasStateInURL()) return;

    const urlState = getStateFromURL();
    if (!urlState || urlState.module !== 'basis' || !urlState.basis) return;

    const bp = urlState.basis;

    // Apply state from URL
    setElement(bp.element);
    setBasis(bp.basis);
    if (bp.shell !== null && bp.shell !== undefined) {
      setShell(bp.shell);
    }

    // Restore comparison state
    if (bp.compare && bp.compare_bases && bp.compare_bases.length > 0) {
      setComparisonBases(bp.compare_bases);
      setComparisonMode(true);
      if (bp.compare_shell_type) {
        setComparisonShellType(bp.compare_shell_type);
      }
    }

    // Restore exponent overrides from deep link (US-053)
    if (bp.exponent_overrides) {
      const overrides: Record<number, number[]> = {};
      for (const [key, values] of Object.entries(bp.exponent_overrides)) {
        const shellIdx = parseInt(key, 10);
        if (!isNaN(shellIdx) && Array.isArray(values) && values.every(v => v > 0)) {
          overrides[shellIdx] = values;
        }
      }
      if (Object.keys(overrides).length > 0) {
        useBasisStore.setState({ exponentOverrides: overrides });
      }
    }

    // Restore overlap distance state from deep link (US-054)
    if (bp.overlap_mode) {
      setOverlapDistanceMode(true);
      if (bp.overlap_traces && bp.overlap_traces.length > 0) {
        const traces: OverlapTraceConfig[] = bp.overlap_traces.map((t) => ({
          elementB: t.element_b,
          basisB: t.basis_b,
          shellIndexB: t.shell_index_b,
        }));
        useBasisStore.setState({ overlapTraces: traces });
      }
    }
  }, [setElement, setBasis, setShell, setComparisonBases, setComparisonMode, setOverlapDistanceMode]);

  // ===========================================================================
  // Deep Link: Write URL state on changes (debounced)
  // ===========================================================================

  const debouncedUpdateURL = useMemo(
    () => debounce((state: RunStateV1) => updateURL(state), URL_DEBOUNCE_MS),
    [],
  );

  useEffect(() => {
    // Don't update URL during initial load
    if (!initializedRef.current) return;

    const state = buildRunState(
      selectedElement,
      selectedBasis,
      selectedShell,
      comparisonMode,
      comparisonBases,
      comparisonShellType,
      exponentOverrides,
      overlapDistanceMode,
      overlapTraces,
    );
    debouncedUpdateURL(state);

    return () => debouncedUpdateURL.cancel();
  }, [
    selectedElement,
    selectedBasis,
    selectedShell,
    comparisonMode,
    comparisonBases,
    comparisonShellType,
    exponentOverrides,
    overlapDistanceMode,
    overlapTraces,
    debouncedUpdateURL,
  ]);

  // ===========================================================================
  // Comparison toggle handler
  // ===========================================================================

  const handleToggleComparison = useCallback(() => {
    if (!comparisonMode) {
      // Deactivate overlap mode when entering comparison mode
      if (overlapDistanceMode) {
        setOverlapDistanceMode(false);
      }
      // Entering comparison mode: ensure primary basis is in the list
      if (!comparisonBases.includes(selectedBasis)) {
        addComparisonBasis(selectedBasis);
      }
      // Add a second basis if list is empty or has only the primary
      const needsSecond = comparisonBases.length === 0 ||
        (comparisonBases.length === 1 && comparisonBases[0] === selectedBasis);
      if (needsSecond) {
        const fallback = selectedBasis === 'sto-3g' ? '3-21g' : 'sto-3g';
        addComparisonBasis(fallback);
      }
    }
    toggleComparisonMode();
  }, [comparisonMode, comparisonBases, selectedBasis, toggleComparisonMode, addComparisonBasis, overlapDistanceMode, setOverlapDistanceMode]);

  // ===========================================================================
  // Overlap distance toggle handler (US-054)
  // ===========================================================================

  const handleToggleOverlap = useCallback(() => {
    toggleOverlapDistanceMode();
  }, [toggleOverlapDistanceMode]);

  // ===========================================================================
  // Auto-select first shell when basis data loads and no shell is selected
  // ===========================================================================

  useEffect(() => {
    if (selectedShell === null && wasmShells && wasmShells.length > 0) {
      setShell(0);
    }
  }, [selectedShell, wasmShells, setShell]);

  // ===========================================================================
  // Collapsible educational section state
  // ===========================================================================

  const [educationOpen, setEducationOpen] = useState(false);

  // ===========================================================================
  // Advanced features section state
  // ===========================================================================

  const [advancedSection, setAdvancedSection] = useState<'none' | 'compare' | 'overlap'>(() => {
    if (comparisonMode) return 'compare';
    if (overlapDistanceMode) return 'overlap';
    return 'none';
  });

  // Keep advancedSection in sync with mode toggles from deep link restore
  useEffect(() => {
    if (comparisonMode) setAdvancedSection('compare');
    else if (overlapDistanceMode) setAdvancedSection('overlap');
  }, [comparisonMode, overlapDistanceMode]);

  const handleAdvancedToggle = useCallback((section: 'compare' | 'overlap') => {
    if (section === 'compare') {
      if (advancedSection === 'compare') {
        // Turning off comparison
        if (comparisonMode) toggleComparisonMode();
        setAdvancedSection('none');
      } else {
        // Turn off overlap if active
        if (overlapDistanceMode) setOverlapDistanceMode(false);
        // Turn on comparison
        handleToggleComparison();
        setAdvancedSection('compare');
      }
    } else if (section === 'overlap') {
      if (advancedSection === 'overlap') {
        // Turning off overlap
        if (overlapDistanceMode) toggleOverlapDistanceMode();
        setAdvancedSection('none');
      } else {
        // Turn off comparison if active
        if (comparisonMode) toggleComparisonMode();
        // Turn on overlap
        if (!overlapDistanceMode) toggleOverlapDistanceMode();
        setAdvancedSection('overlap');
      }
    }
  }, [advancedSection, comparisonMode, overlapDistanceMode, toggleComparisonMode, handleToggleComparison, toggleOverlapDistanceMode, setOverlapDistanceMode]);

  return (
    <div className="max-w-6xl mx-auto">
      {/* Header */}
      <div className="mb-6">
        <h1 className="text-3xl font-bold text-slate-800 mb-2">
          Module A: Basis Explorer
        </h1>
        <p className="text-slate-600">
          Explore 6 Gaussian basis sets (STO-3G through cc-pVDZ) for all 18
          supported elements (H-Ar). Select an element and basis set, then inspect
          the shell structure, radial profiles, and primitive parameters.
        </p>
      </div>

      {/* ================================================================
          PRIMARY AREA: Controls (left) + Radial Profile Plot (right)
          ================================================================ */}
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-6">
        {/* Left column: Element selector, Basis dropdown, Shell list */}
        <div className="lg:col-span-1 space-y-4">
          <ElementSelector />
          <BasisSelector />
          <ShellList />
        </div>

        {/* Right column: Data first, then Plot, then Advanced */}
        <div className="lg:col-span-2 space-y-4">
          {/* Primitive Data Table — numerical data FIRST for intuitive understanding */}
          {selectedShellData && (
            <ShellDataTable
              exponents={selectedShellData.exponents}
              coefficients={selectedShellData.coefficients}
              angularMomentum={selectedShellData.angularMomentum}
              angularMomentumLabel={selectedShellData.angularMomentumLabel}
              shellLabel={selectedShellData.shellLabel}
            />
          )}

          {/* Radial Profile Plot — visualizes the data above */}
          {comparisonMode && selectedShell !== null ? (
            <div className="space-y-4">
              <BasisComparisonPlot
                profiles={comparisonProfiles}
                shellTypeLabel={comparisonShellType}
                atomicNumber={selectedElement}
                isLoading={comparisonLoading}
                yAxisScale={radialProfileScale}
                onScaleChange={setRadialProfileScale}
                minHeight={420}
              />
              {comparisonError && (
                <div className="bg-red-50 rounded-lg px-4 py-3 border border-red-200">
                  <p className="text-xs text-red-800">{comparisonError}</p>
                </div>
              )}
              {comparisonProfiles.size >= 2 && (
                <DiscussionPrompt
                  comparedBases={effectiveComparisonBases}
                  shellType={comparisonShellType}
                  atomicNumber={selectedElement}
                />
              )}
            </div>
          ) : (
            <RadialProfilePlot
              profileData={profileData}
              isLoading={profileLoading}
              yAxisScale={radialProfileScale}
              onScaleChange={setRadialProfileScale}
              minHeight={420}
              modifiedProfileValues={modifiedProfile}
              exponentsModified={exponentsModified}
            />
          )}

          {/* Exponent sliders — modify primitives, see effect in plot above */}
          {selectedShellData && !comparisonMode && profileData && defaultExponents.length > 0 && (
            <ExponentSliders
              defaultExponents={defaultExponents}
              currentExponents={currentExponents}
              angularMomentumLabel={selectedShellData.angularMomentumLabel}
              shellLabel={selectedShellData.shellLabel}
              onExponentChange={handleExponentChange}
              onReset={handleExponentReset}
              isModified={exponentsModified}
              overlapMetric={overlapMetric}
            />
          )}

          {/* Advanced features: Compare Bases / Overlap — below the plot */}
          <div className="flex items-center gap-3">
            <span className="text-xs font-medium text-slate-400 uppercase tracking-wider">Advanced:</span>
            <button
              type="button"
              onClick={() => handleAdvancedToggle('compare')}
              className={`
                px-3 py-1.5 text-xs rounded-lg border transition-colors
                focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
                ${advancedSection === 'compare'
                  ? 'bg-primary-50 border-primary-300 text-primary-700 font-medium'
                  : 'bg-white border-slate-200 text-slate-600 hover:bg-slate-50 hover:border-slate-300'
                }
              `}
              aria-pressed={advancedSection === 'compare'}
            >
              Compare Basis Sets
            </button>
            <button
              type="button"
              onClick={() => handleAdvancedToggle('overlap')}
              disabled={selectedShell === null}
              className={`
                px-3 py-1.5 text-xs rounded-lg border transition-colors
                focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
                ${advancedSection === 'overlap'
                  ? 'bg-primary-50 border-primary-300 text-primary-700 font-medium'
                  : selectedShell === null
                    ? 'bg-slate-50 border-slate-200 text-slate-400 cursor-not-allowed'
                    : 'bg-white border-slate-200 text-slate-600 hover:bg-slate-50 hover:border-slate-300'
                }
              `}
              aria-pressed={advancedSection === 'overlap'}
            >
              Overlap vs. Distance
            </button>
          </div>

          {/* Comparison controls + plot */}
          {advancedSection === 'compare' && (
            <BasisComparisonControls
              isActive={comparisonMode}
              onToggle={handleToggleComparison}
              selectedBases={effectiveComparisonBases}
              onAddBasis={addComparisonBasis}
              onRemoveBasis={removeComparisonBasis}
              primaryBasis={selectedBasis}
              shellTypes={availableShellTypes}
              selectedShellType={comparisonShellType}
              onShellTypeChange={setComparisonShellType}
              hasShellSelected={selectedShell !== null}
            />
          )}

          {/* Overlap distance controls + plot */}
          {advancedSection === 'overlap' && (
            <div className="space-y-4">
              <OverlapDistanceControls
                isActive={overlapDistanceMode}
                onToggle={handleToggleOverlap}
                overlapTraces={overlapTraces}
                onAddTrace={addOverlapTrace}
                onRemoveTrace={removeOverlapTrace}
                onClearTraces={clearOverlapTraces}
                primaryElement={selectedElement}
                primaryBasis={selectedBasis}
                primaryShell={selectedShell}
                hasShellSelected={selectedShell !== null}
                comparisonModeActive={comparisonMode}
              />
              {overlapDistanceMode && selectedShell !== null && (
                <OverlapDistancePlot
                  traces={overlapPlotTraces}
                  isLoading={overlapLoading}
                  yAxisScale={overlapPlotScale}
                  onScaleChange={setOverlapPlotScale}
                  minHeight={350}
                />
              )}
            </div>
          )}

          {/* Primitive data + exponent sliders are now above the plot */}
        </div>
      </div>

      {/* Advanced features are now inside the right column above */}

      {/* ================================================================
          EDUCATIONAL BACKGROUND (collapsible)
          ================================================================ */}
      <div className="bg-slate-100 rounded-xl overflow-hidden">
        <button
          type="button"
          onClick={() => setEducationOpen(!educationOpen)}
          className="w-full flex items-center justify-between p-5 text-left focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-inset focus-visible:ring-primary-500"
          aria-expanded={educationOpen}
          aria-controls="education-content"
        >
          <h2 className="font-semibold text-slate-700">
            Understanding Basis Sets
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
          <div id="education-content" className="px-5 pb-5">
            <div className="text-slate-600 text-sm space-y-3">
              <p>
                In quantum chemistry, atomic orbitals are approximated using
                Gaussian-type functions. A <strong>basis set</strong> defines
                the collection of these functions used for each atom in a molecule.
              </p>

              <div className="grid md:grid-cols-2 gap-4 mt-4">
                <div className="bg-white rounded-lg p-4 border border-slate-200">
                  <h3 className="font-semibold text-slate-700 mb-2">Primitives vs. Contractions</h3>
                  <p className="text-xs text-slate-600">
                    A <strong>primitive</strong> Gaussian is a single exponential function.
                    A <strong>contracted</strong> shell combines several primitives with
                    fixed coefficients to better approximate atomic orbitals while keeping
                    computation tractable.
                  </p>
                </div>

                <div className="bg-white rounded-lg p-4 border border-slate-200">
                  <h3 className="font-semibold text-slate-700 mb-2">Angular Momentum</h3>
                  <p className="text-xs text-slate-600">
                    Shells are classified by angular momentum quantum number:
                    <strong> s</strong> (l=0, 1 component),
                    <strong> p</strong> (l=1, 3 components),
                    <strong> d</strong> (l=2, 5 or 6 components).
                    Higher angular momentum shells capture more directional bonding character.
                    For d-shells, <strong>Cartesian</strong> representation uses 6 components
                    (xx, yy, zz, xy, xz, yz) while <strong>spherical</strong> uses 5
                    (m = -2 to +2).
                  </p>
                </div>

                <div className="bg-white rounded-lg p-4 border border-slate-200">
                  <h3 className="font-semibold text-slate-700 mb-2">Split-Valence Basis Sets</h3>
                  <p className="text-xs text-slate-600">
                    In a split-valence basis (e.g., 3-21G, 6-31G), core orbitals use a
                    single contraction while valence orbitals are split into inner and outer
                    parts. This allows valence orbitals to adjust their size during
                    the SCF procedure.
                  </p>
                </div>

                <div className="bg-white rounded-lg p-4 border border-slate-200">
                  <h3 className="font-semibold text-slate-700 mb-2">Polarization and Diffuse Functions</h3>
                  <p className="text-xs text-slate-600">
                    The <strong>*</strong> (star) adds d-type polarization functions to
                    heavy atoms, improving descriptions of bond angles. The <strong>+</strong>
                    adds diffuse s and p functions, which improve descriptions of anions
                    and lone pairs extending far from the nucleus.
                  </p>
                </div>
              </div>

              <div className="mt-4 grid md:grid-cols-6 gap-2">
                <div className="bg-blue-50 rounded-lg p-3 border border-blue-100 text-center">
                  <p className="font-semibold text-blue-800 text-sm">STO-3G</p>
                  <p className="text-xs text-blue-600">Minimal</p>
                </div>
                <div className="bg-blue-50 rounded-lg p-3 border border-blue-100 text-center">
                  <p className="font-semibold text-blue-800 text-sm">3-21G</p>
                  <p className="text-xs text-blue-600">Split-valence</p>
                </div>
                <div className="bg-blue-50 rounded-lg p-3 border border-blue-100 text-center">
                  <p className="font-semibold text-blue-800 text-sm">6-31G</p>
                  <p className="text-xs text-blue-600">Split-valence</p>
                </div>
                <div className="bg-blue-50 rounded-lg p-3 border border-blue-100 text-center">
                  <p className="font-semibold text-blue-800 text-sm">6-31G*</p>
                  <p className="text-xs text-blue-600">+ polarization</p>
                </div>
                <div className="bg-blue-50 rounded-lg p-3 border border-blue-100 text-center">
                  <p className="font-semibold text-blue-800 text-sm">6-31+G*</p>
                  <p className="text-xs text-blue-600">+ diffuse</p>
                </div>
                <div className="bg-blue-50 rounded-lg p-3 border border-blue-100 text-center">
                  <p className="font-semibold text-blue-800 text-sm">cc-pVDZ</p>
                  <p className="text-xs text-blue-600">Correlation-consistent</p>
                </div>
              </div>

              <p className="text-xs text-slate-500 mt-4">
                All 6 basis sets support elements H through Ar (Z=1-18).
                Basis set data sourced from PySCF 2.11.0 and the Basis Set Exchange.
              </p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

export default ModuleD;
