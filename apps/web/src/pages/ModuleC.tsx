/**
 * ModuleC - SCF Sandbox page.
 *
 * Interactive self-consistent field (SCF) computation with
 * optional DIIS acceleration. Run Hartree-Fock calculations
 * on preset molecular systems and visualize convergence.
 *
 * Phase 2 (US-035): Integrates a 3D molecular viewer with
 * progressive disclosure, React.lazy() code splitting, and
 * WebGL fallback handling.
 *
 * Phase 4-5: Task-oriented workflow tabs separate the page into
 * five workflows: Single Point, Optimize, Frequency, PES Scan,
 * and Compare. Each tab has its own results display while sharing
 * molecule/basis selection state.
 *
 * @module pages/ModuleC
 */

import React, { useEffect, useRef, useMemo, useCallback, useState, Suspense } from 'react';
import { ScfControlsPanel, ScfResultDisplay, SystemInfoPanel, PesScanPanel, PesCurvePlot, CoordinateTrackingPanel, OrbitalSelector, IsovalueSlider, BasisComparisonPanel, BasisComparisonChart, BasisComparisonTable, DensityPanel, DensityCrossSection, DftInfoPanel, DftComparisonPanel, OptimizeButton, OptimizationProgressPlot, OptimizationResultPanel, FrequencyTab } from '../components/scf';
import type { DiatomicInfo } from '../components/scf';
import { ELEMENT_TO_Z, Z_TO_ELEMENT } from '../lib/elements';
import { useDftCompare } from '../hooks/useDftCompare';
import { useOptimizer } from '../hooks/useOptimizer';
import { useOrbitalViewer } from '../hooks/useOrbitalViewer';
import { useDensityViewer } from '../hooks/useDensityViewer';
import { CopyLinkButton, ExportButton, ImportButton } from '../components/common';
import { Math, MathBlock } from '../components/common/Math';
import { useScf } from '../hooks/useScf';
import { useScfStore, getMoleculePreset } from '../stores/scfStore';
import { getPreset } from '../worker/presets';
import { getStateFromURL, hasStateInURL, updateURL, clearURL } from '../lib/url';
import { createArtifact, downloadArtifact, restoreArtifact, getArtifactModule } from '../lib/artifact';
import { downloadPovRay } from '../lib/povrayExport';
import { isDftMethod } from '../types/dft';
import { DeepLinkError } from '../components/DeepLinkError';
import { DEFAULT_SCF_STATE, APP_VERSION } from '../types/run-state';
import type { RunStateV1 } from '../types/run-state';
import type { ScfArtifactResult, RunArtifactV1, PesScanArtifactData } from '../types/run-artifact';
import type { Atom3D } from '../components/viewer3d/AtomSpheres';

/**
 * Reverse mapping from atomic number to element symbol.
 * Used to reconstruct Atom3D from [Z, x, y, z] during trajectory animation.
 */
const ELEMENT_TO_Z_REVERSE = Z_TO_ELEMENT;

/**
 * Lazy-loaded MoleculeViewer component.
 *
 * Three.js and React Three Fiber are code-split into a separate chunk
 * that is only fetched when the 3D viewer panel is first rendered.
 * This keeps the initial Module E bundle size small.
 */
const LazyMoleculeViewer = React.lazy(() =>
  import('../components/viewer3d').then((m) => ({ default: m.MoleculeViewer }))
);

/**
 * Lazy-loaded OrbitalSurface component.
 *
 * Code-split with the viewer3d chunk so it's only loaded when Module E
 * renders the 3D viewer.
 */
const LazyOrbitalSurface = React.lazy(() =>
  import('../components/viewer3d').then((m) => ({ default: m.OrbitalSurface }))
);

/**
 * Lazy-loaded DensitySurface component.
 *
 * Code-split with the viewer3d chunk. Renders the electron density
 * isosurface as a single semi-transparent teal mesh.
 *
 * @see US-062 Density Isosurface & Cross-Sections
 */
const LazyDensitySurface = React.lazy(() =>
  import('../components/viewer3d').then((m) => ({ default: m.DensitySurface }))
);

/**
 * Lazy-loaded GhostAtoms component.
 *
 * Renders translucent gray spheres for the initial (pre-optimization)
 * geometry overlay in trajectory animation.
 *
 * @see US-075 Optimization UI + Trajectory Animation (AC4)
 */
const LazyGhostAtoms = React.lazy(() =>
  import('../components/viewer3d').then((m) => ({ default: m.GhostAtoms }))
);

/**
 * Workflow tab for Module E task-oriented navigation.
 *
 * - "single-point": Standard SCF calculation with full result inspection
 * - "optimize": Geometry optimization workflow with trajectory animation
 * - "pes-scan": Potential energy surface scanning along internal coordinates
 * - "compare": Basis set / method comparison across multiple configurations
 */
type WorkflowTab =
  | 'single-point'
  | 'optimize'
  | 'frequency'
  | 'pes-scan'
  | 'compare';

/**
 * Workflow tab definitions for rendering the tab bar.
 */
const WORKFLOW_TABS: { id: WorkflowTab; label: string; description: string }[] = [
  { id: 'single-point', label: 'Single Point', description: 'Run a single SCF calculation and inspect results' },
  { id: 'optimize', label: 'Optimize', description: 'Optimize molecular geometry' },
  { id: 'frequency', label: 'Frequency', description: 'Vibrational analysis: frequencies, IR/Raman spectra, thermochemistry' },
  { id: 'pes-scan', label: 'PES Scan', description: 'Scan potential energy surface along internal coordinates' },
  { id: 'compare', label: 'Compare', description: 'Compare methods and basis sets' },
];

/**
 * Simple debounce for URL updates.
 */
function debounce<T extends (...args: Parameters<T>) => void>(
  fn: T,
  delay: number
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
 * Debounce delay for URL updates (ms).
 *
 * Longer than typical to avoid excessive URL changes during parameter adjustment.
 */
const URL_DEBOUNCE_MS = 300;

/**
 * Extract Atom3D array from current system state.
 *
 * Handles both preset and custom geometry modes.
 * Returns null if no geometry data is available.
 *
 * @param systemId - Current system ID (for preset lookup)
 * @param inputMode - Current input mode (preset or custom)
 */
function useAtoms3D(
  systemId: string,
  inputMode: { mode: 'preset' } | { mode: 'custom'; geometry: { atoms: Array<{ symbol: string; xyz: [number, number, number] | number[] }>; units: string }; basisSet: string }
): Atom3D[] | null {
  return useMemo(() => {
    if (inputMode.mode === 'custom') {
      // Custom mode: use geometry from input mode
      const atoms = inputMode.geometry.atoms;
      if (!atoms || atoms.length === 0) return null;
      return atoms.map((a) => ({
        symbol: a.symbol,
        position: [a.xyz[0], a.xyz[1], a.xyz[2]] as [number, number, number],
      }));
    }

    // Preset mode: look up geometry from preset data
    const preset = getPreset(systemId);
    if (!preset?.geometry?.atoms || preset.geometry.atoms.length === 0) return null;

    return preset.geometry.atoms.map((a) => ({
      symbol: a.symbol,
      position: [a.xyz[0], a.xyz[1], a.xyz[2]] as [number, number, number],
    }));
  }, [systemId, inputMode]);
}

/**
 * Detect diatomic molecule and extract atom info for PES scan.
 *
 * Returns DiatomicInfo if the current system has exactly 2 atoms
 * (suitable for PES bond-length scanning), or null otherwise.
 *
 * @param systemId - Current system ID (for preset lookup)
 * @param inputMode - Current input mode (preset or custom)
 * @returns DiatomicInfo or null
 *
 * @see US-040 PES Scan UI
 */
function useDiatomicInfo(
  systemId: string,
  inputMode: { mode: 'preset' } | { mode: 'custom'; geometry: { atoms: Array<{ symbol: string; xyz: [number, number, number] | number[] }>; units: string }; basisSet: string }
): DiatomicInfo | null {
  return useMemo(() => {
    let atoms: Array<{ symbol: string }> | undefined;
    let basisName: string;

    if (inputMode.mode === 'custom') {
      atoms = inputMode.geometry.atoms;
      basisName = inputMode.basisSet;
    } else {
      const preset = getPreset(systemId);
      if (!preset?.geometry?.atoms) return null;
      atoms = preset.geometry.atoms;
      basisName = preset.basis_id;
    }

    if (!atoms || atoms.length !== 2) return null;

    const symbolA = atoms[0].symbol;
    const symbolB = atoms[1].symbol;
    const atomAZ = ELEMENT_TO_Z[symbolA];
    const atomBZ = ELEMENT_TO_Z[symbolB];

    // Both atoms must have known atomic numbers
    if (!atomAZ || !atomBZ) return null;

    return { symbolA, symbolB, atomAZ, atomBZ, basisName };
  }, [systemId, inputMode]);
}

/**
 * 3D Molecule Viewer panel with progressive disclosure.
 *
 * Renders a large panel containing the lazy-loaded MoleculeViewer
 * alongside orbital controls. Designed for the right-side main content
 * area with generous vertical space for the 3D scene.
 *
 * Only shown when atom geometry data is available (AC-3, AC-4).
 * Three.js is code-split via React.lazy() and loaded on demand (AC-2).
 */
function MoleculePanel({
  atoms,
  showLabels,
  onToggleLabels,
  orbitalControls,
  headerControls,
  orbitalMeshData,
  children,
}: {
  atoms: Atom3D[];
  showLabels: boolean;
  onToggleLabels: (show: boolean) => void;
  /** Orbital controls rendered beside the 3D viewer on large screens */
  orbitalControls?: React.ReactNode;
  /** Extra controls rendered in the header bar (e.g., charge overlay toggle) */
  headerControls?: React.ReactNode;
  /** Orbital mesh data for POV-Ray export (positive + negative lobes) */
  orbitalMeshData?: { positive: import('../components/viewer3d/OrbitalSurface').LobeMeshData | null; negative: import('../components/viewer3d/OrbitalSurface').LobeMeshData | null } | null;
  children?: React.ReactNode;
}) {
  // Track camera state for POV-Ray export (matches viewer orientation)
  const cameraStateRef = useRef<import('../components/viewer3d/MoleculeViewer').CameraState | null>(null);
  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      {/* Header */}
      <div className="px-4 py-3 bg-slate-50 border-b border-slate-200 flex items-center justify-between">
        <h3 className="text-sm font-semibold text-slate-700">
          3D Structure &amp; Orbitals
        </h3>
        <div className="flex items-center gap-3">
          {headerControls}
          <label className="flex items-center gap-1.5 text-xs text-slate-600 cursor-pointer select-none">
            <input
              type="checkbox"
              checked={showLabels}
              onChange={(e) => onToggleLabels(e.target.checked)}
              className="rounded border-slate-300 text-blue-600 focus:ring-blue-500 h-3.5 w-3.5"
            />
            Labels
          </label>
          <button
            type="button"
            onClick={(e) => {
              // Walk up to the panel container, then find the R3F canvas.
              // preserveDrawingBuffer is set on the Canvas so toDataURL works.
              const panel = (e.target as HTMLElement).closest('.bg-white');
              const canvas = panel?.querySelector('canvas');
              if (canvas instanceof HTMLCanvasElement) {
                // 2x supersampling for publication-quality export:
                // Create a high-res offscreen canvas, draw the WebGL canvas
                // scaled up, then export the high-res version.
                const scale = 2;
                const w = canvas.width * scale;
                const h = canvas.height * scale;
                const offscreen = document.createElement('canvas');
                offscreen.width = w;
                offscreen.height = h;
                const ctx = offscreen.getContext('2d');
                if (ctx) {
                  // White background (WebGL canvas may have transparency)
                  ctx.fillStyle = '#f8fafc';
                  ctx.fillRect(0, 0, w, h);
                  ctx.imageSmoothingEnabled = true;
                  ctx.imageSmoothingQuality = 'high';
                  ctx.drawImage(canvas, 0, 0, w, h);
                  const url = offscreen.toDataURL('image/png');
                  const a = document.createElement('a');
                  a.href = url;
                  a.download = 'iqcp-molecule-3d.png';
                  a.click();
                }
              }
            }}
            className="text-xs text-slate-500 hover:text-slate-700 flex items-center gap-1"
            title="Export 3D view as high-resolution PNG (2x)"
          >
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
            </svg>
            PNG
          </button>
          <button
            type="button"
            onClick={() => {
              downloadPovRay(
                atoms,
                'iqcp-molecule.pov',
                `${atoms.map(a => a.symbol).join('')} structure`,
                orbitalMeshData ?? undefined,
                cameraStateRef.current ?? undefined,
              );
            }}
            className="text-xs text-slate-500 hover:text-slate-700 flex items-center gap-1"
            title="Export POV-Ray scene file for ray-traced rendering"
          >
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
            </svg>
            POV-Ray
          </button>
        </div>
      </div>

      {/* Content: 3D viewer + orbital controls side by side on lg */}
      <div className="flex flex-col lg:flex-row">
        {/* 3D Viewer — large presentation area */}
        <div className="flex-1 h-[28rem]">
          <Suspense
            fallback={
              <div className="flex items-center justify-center h-full bg-slate-50">
                <div className="text-center">
                  <div className="animate-spin rounded-full h-6 w-6 border-b-2 border-blue-600 mx-auto mb-2" />
                  <span className="text-xs text-slate-500">Loading 3D viewer...</span>
                </div>
              </div>
            }
          >
            <LazyMoleculeViewer
              atoms={atoms}
              showLabels={showLabels}
              onCameraChange={(state) => { cameraStateRef.current = state; }}
            >
              {children}
            </LazyMoleculeViewer>
          </Suspense>
        </div>

        {/* Orbital controls sidebar (right of 3D viewer) */}
        {orbitalControls && (
          <div className="lg:w-56 border-t lg:border-t-0 lg:border-l border-slate-200 p-4 space-y-4 bg-slate-50/50">
            {orbitalControls}
          </div>
        )}
      </div>
    </div>
  );
}

/**
 * Module E: SCF Sandbox page.
 *
 * Provides interactive controls for running RHF SCF calculations
 * on preset molecular systems. Unlike Modules A/B, this uses
 * explicit Run/Cancel buttons due to longer computation times.
 *
 * Phase 2 (US-035): Includes a lazy-loaded 3D molecular viewer
 * that appears when a molecule is selected (progressive disclosure).
 */
function ModuleC() {
  const { isReady, workerError } = useScf();

  // Store state
  const systemId = useScfStore((state) => state.systemId);
  const convergenceProfile = useScfStore((state) => state.convergenceProfile);
  const maxIterations = useScfStore((state) => state.maxIterations);
  const useDiis = useScfStore((state) => state.useDiis);
  const useSpherical = useScfStore((state) => state.useSpherical);
  const mode = useScfStore((state) => state.mode);
  const compute = useScfStore((state) => state.compute);
  const history = useScfStore((state) => state.history);
  const urlInitialized = useScfStore((state) => state.urlInitialized);
  const initializeFromURL = useScfStore((state) => state.initializeFromURL);
  const reset = useScfStore((state) => state.reset);

  // Phase 3: DFT method selection (US-070)
  const method = useScfStore((state) => state.method);
  const ksResult = useScfStore((state) => state.ksResult);
  const isRunning = compute.status === 'running';

  // Molecule selection (needed for SystemInfoPanel data)
  const selectedMolecule = useScfStore((state) => state.selectedMolecule);

  // Phase 2: Viewer state (US-035)
  const inputMode = useScfStore((state) => state.inputMode);
  const viewerState = useScfStore((state) => state.viewerState);
  const setShowLabels = useScfStore((state) => state.setShowLabels);
  const setSelectedOrbital = useScfStore((state) => state.setSelectedOrbital);
  const setIsovalue = useScfStore((state) => state.setIsovalue);

  // Phase 4: Population analysis and charge overlay (S4)
  const populationResult = useScfStore((state) => state.populationResult);
  const showChargeOverlay = useScfStore((state) => state.showChargeOverlay);
  const setShowChargeOverlay = useScfStore((state) => state.setShowChargeOverlay);

  // Extract atoms for the 3D viewer from current system
  const atoms3D = useAtoms3D(systemId, inputMode);

  // Phase 4 (S4): Merge Mulliken charges into atoms when charge overlay is active
  const chargedAtoms3D = useMemo(() => {
    if (!showChargeOverlay || !populationResult || !atoms3D) return atoms3D;
    return atoms3D.map((atom, i) => ({
      ...atom,
      charge: populationResult.atoms[i]?.mullikenCharge ?? undefined,
    }));
  }, [atoms3D, populationResult, showChargeOverlay]);

  // Phase 2: Diatomic detection for PES scan (US-040)
  const diatomicInfo = useDiatomicInfo(systemId, inputMode);

  // Phase 2 (US-044): Orbital viewer — extract atoms in worker format [Z, x, y, z]
  const workerAtoms = useMemo(() => {
    if (!atoms3D) return null;
    return atoms3D.map((atom) => {
      const z = ELEMENT_TO_Z[atom.symbol] ?? 1;
      return [z, atom.position[0], atom.position[1], atom.position[2]] as [number, number, number, number];
    });
  }, [atoms3D]);

  // Phase 2 (US-044): Determine basis set name for orbital evaluation
  const basisName = useMemo(() => {
    if (inputMode.mode === 'custom' && inputMode.basisSet) {
      return inputMode.basisSet;
    }
    // Preset systems — extract basis from system ID
    const preset = getPreset(systemId);
    return preset?.basis_id ?? 'sto-3g';
  }, [inputMode, systemId]);

  // Phase 3 (US-070): DFT comparison hook
  const { compare: dftCompare, isComparing: isDftComparing } = useDftCompare({
    atoms: workerAtoms,
    basisName,
    isReady,
  });

  // Phase 3 (US-075): Geometry optimization hook
  const {
    optimize: runOptimization,
    cancel: cancelOptimization,
    loading: isOptimizing,
  } = useOptimizer({
    atoms: workerAtoms,
    basisName,
    method,
    isReady,
  });
  const optimizationState = useScfStore((state) => state.optimizationState);
  const setTrajectoryStep = useScfStore((state) => state.setTrajectoryStep);
  const setShowGhostOverlay = useScfStore((state) => state.setShowGhostOverlay);
  const applyOptimizedGeometry = useScfStore((state) => state.applyOptimizedGeometry);
  const [optimizedApplied, setOptimizedApplied] = useState(false);
  const [educationOpen, setEducationOpen] = useState(false);

  // Workflow tab state — UI-only, not persisted to URL
  const [workflowTab, setWorkflowTab] = useState<WorkflowTab>('single-point');

  // Reset the applied flag when optimization restarts or system changes
  useEffect(() => {
    setOptimizedApplied(false);
  }, [optimizationState.result, systemId]);

  // Compute atoms for the current trajectory step (for 3D viewer animation)
  const trajectoryAtoms3D = useMemo(() => {
    if (!optimizationState.result || !workerAtoms) return null;
    const { steps, totalSteps } = optimizationState.result;
    const stepIdx = globalThis.Math.min(optimizationState.trajectoryStep, totalSteps);
    const step = steps.find((s) => s.step === stepIdx);
    if (!step || step.geometry.length === 0) return null;

    // Use the geometry from the step (preserving original Z values)
    return workerAtoms.map((a, i) => ({
      symbol: ELEMENT_TO_Z_REVERSE[a[0]] ?? 'X',
      position: step.geometry[i] ?? [a[1], a[2], a[3]],
    })) as Atom3D[];
  }, [optimizationState.result, optimizationState.trajectoryStep, workerAtoms]);

  // Initial geometry for ghost overlay (step 0)
  const ghostAtoms3D = useMemo(() => {
    if (!optimizationState.result || !workerAtoms || !optimizationState.showGhostOverlay) return null;
    const step0 = optimizationState.result.steps.find((s) => s.step === 0);
    if (!step0 || step0.geometry.length === 0) {
      // Fall back to atoms3D (pre-optimization positions)
      return atoms3D;
    }
    return workerAtoms.map((a, i) => ({
      symbol: ELEMENT_TO_Z_REVERSE[a[0]] ?? 'X',
      position: step0.geometry[i] ?? [a[1], a[2], a[3]],
    })) as Atom3D[];
  }, [optimizationState.result, optimizationState.showGhostOverlay, workerAtoms, atoms3D]);

  // Callback to apply optimized geometry
  const handleApplyOptimizedGeometry = useCallback(() => {
    if (!optimizationState.result || !workerAtoms) return;
    const finalGeom = optimizationState.result.finalGeometry;
    const newAtoms: [number, number, number, number][] = workerAtoms.map((a, i) => [
      a[0],
      finalGeom[i][0],
      finalGeom[i][1],
      finalGeom[i][2],
    ]);
    applyOptimizedGeometry(newAtoms, basisName);
    setOptimizedApplied(true);
  }, [optimizationState.result, workerAtoms, basisName, applyOptimizedGeometry]);

  // Phase 2 (US-044): SCF result data for orbital viewer
  const scfResult = compute.status === 'success' ? compute.result : null;
  const orbitalEnergies = scfResult?.orbitalEnergies ?? null;
  const moCoefficients = scfResult?.matrices?.moCoefficients ?? null;
  const densityMatrix = scfResult?.matrices?.densityMatrix ?? null;
  const nbf = scfResult?.matrices?.nbf ?? 0;
  const nElectrons = orbitalEnergies ? orbitalEnergies.nOccupied * 2 : 0;

  // Derive correct system info for SystemInfoPanel.
  // The molecule preset has accurate nelec and geometry; nbf and e_nuc
  // come from the KS/SCF result (only known after computation).
  const moleculePreset = useMemo(() => getMoleculePreset(selectedMolecule), [selectedMolecule]);
  const sysInfoNelec = moleculePreset?.nelec;
  const sysInfoNbf = ksResult?.nBasis ?? (scfResult?.matrices?.nbf) ?? undefined;
  const sysInfoEnuc = ksResult?.energyNuc ?? undefined;
  const sysInfoGeometry = useMemo(() => {
    if (inputMode.mode === 'custom') return inputMode.geometry;
    if (moleculePreset) {
      return {
        atoms: moleculePreset.atoms.map((a) => ({ symbol: a.symbol, xyz: a.xyz as number[] })),
        units: 'bohr',
      };
    }
    return undefined;
  }, [inputMode, moleculePreset]);

  // Phase 2 (US-044): Orbital isosurface computation
  const orbitalViewer = useOrbitalViewer({
    selectedOrbital: viewerState.selectedOrbital,
    isovalue: viewerState.isovalue,
    moCoefficients,
    nbf,
    atoms: workerAtoms,
    basisName,
    isReady,
    useSpherical,
  });

  // Phase 3 (US-062, US-063): Density isosurface computation
  const densityViewer = useDensityViewer({
    active: viewerState.showDensity && scfResult !== null,
    densityMatrix,
    atoms: workerAtoms,
    basisName,
    nElectrons,
    useSpherical,
    isovalue: viewerState.densityIsovalue,
    densityMode: viewerState.densityMode,
    diffIsovalue: viewerState.diffIsovalue,
    isReady,
  });

  // Phase 3 (US-062): Compute grid bounds for DensityPanel slider ranges
  const densityGridBounds = useMemo(() => {
    const grid = densityViewer.gridData;
    if (!grid) return null;
    const [ox, oy, oz] = grid.gridOrigin;
    const [nx, ny, nz] = grid.gridDims;
    const s = grid.gridSpacing;
    return {
      min: [ox, oy, oz] as [number, number, number],
      max: [ox + (nx - 1) * s, oy + (ny - 1) * s, oz + (nz - 1) * s] as [number, number, number],
    };
  }, [densityViewer.gridData]);

  // Track if URL had invalid state
  const invalidURLRef = useRef(false);

  // Initialize state from URL on mount
  useEffect(() => {
    if (urlInitialized) return;

    if (hasStateInURL()) {
      const urlState = getStateFromURL();

      if (urlState && urlState.module === 'scf' && urlState.scf) {
        // Valid SCF state from URL
        initializeFromURL(urlState.scf);
      } else if (urlState && urlState.module !== 'scf') {
        // Valid state but wrong module - use defaults
        initializeFromURL(DEFAULT_SCF_STATE.scf!);
      } else {
        // Invalid URL state
        invalidURLRef.current = true;
        initializeFromURL(DEFAULT_SCF_STATE.scf!);
      }
    } else {
      // No URL state - use defaults
      initializeFromURL(DEFAULT_SCF_STATE.scf!);
    }
  }, [urlInitialized, initializeFromURL]);

  // Phase 2 (US-040): PES scan state for URL encoding
  const pesRMin = useScfStore((state) => state.pesState.rMin);
  const pesRMax = useScfStore((state) => state.pesState.rMax);
  const pesNPoints = useScfStore((state) => state.pesState.nPoints);

  // Phase 2 (US-044): Orbital state for URL encoding
  const selectedOrbital = viewerState.selectedOrbital;
  const isovalue = viewerState.isovalue;

  // Phase 2 (US-045): Basis comparison state for URL encoding
  const comparisonBases = useScfStore((state) => state.basisCompareState.selectedBases);

  // Phase 3 (US-062, US-063): Density state for URL encoding
  const densityIsovalue = viewerState.densityIsovalue;
  const diffIsovalue = viewerState.diffIsovalue;
  const densityMode = viewerState.densityMode;
  const planeAxis = viewerState.planeAxis;
  const planePosition = viewerState.planePosition;
  const colorScale = viewerState.colorScale;
  const showDensity = viewerState.showDensity;

  // Phase 5 (US-103): Frequency state for URL encoding
  const frequencyState = useScfStore((state) => state.frequencyState);
  const freqSelectedMode = frequencyState.selectedMode;
  const freqTemperatureK = frequencyState.temperatureK;
  const freqPressurePa = frequencyState.pressurePa;
  const freqBroadeningKind = frequencyState.broadeningKind;
  const freqFwhmCm1 = frequencyState.fwhmCm1;
  const freqSpectrumTab = frequencyState.spectrumTab;
  const freqShowArrows = frequencyState.showDisplacementArrows;
  const freqUnitsMode = frequencyState.unitsMode;

  // Create debounced URL update function
  const debouncedUpdateURL = useMemo(
    () =>
      debounce(
        (
          newSystemId: string,
          newConv: 'loose' | 'medium' | 'tight',
          newMaxIter: number,
          newDiis: boolean,
          newMethod: string,
          newPesRMin: number,
          newPesRMax: number,
          newPesNPoints: number,
          newOrbital: number | null,
          newIsovalue: number,
          newComparisonBases: string[],
          newDensityMode: string,
          newDensityIsovalue: number,
          newDiffIsovalue: number,
          newPlaneAxis: string,
          newPlanePosition: number,
          newColorScale: string,
          newShowDensity: boolean,
          // Phase 5 (US-103): Frequency state params
          newFreqSelectedMode: number | null,
          newFreqTempK: number,
          newFreqPressPa: number,
          newFreqBroadening: string,
          newFreqFwhm: number,
          newFreqTab: string,
          newFreqArrows: boolean,
          newFreqUnits: string
        ) => {
          // Only include PES config if non-default values
          const pesConfig =
            newPesRMin !== 0.5 || newPesRMax !== 5.0 || newPesNPoints !== 20
              ? { r_min: newPesRMin, r_max: newPesRMax, n_points: newPesNPoints }
              : undefined;

          // Only include orbital config if non-default values (US-044)
          const orbitalConfig =
            newOrbital !== null || newIsovalue !== 0.03
              ? { orbital: newOrbital ?? undefined, isovalue: newIsovalue !== 0.03 ? newIsovalue : undefined }
              : undefined;

          // Only include basis comparison if non-default (US-045)
          const defaultBases = ['sto-3g', '6-31g'];
          const basesChanged =
            newComparisonBases.length !== defaultBases.length ||
            newComparisonBases.some((b, i) => b !== defaultBases[i]);
          const comparisonConfig = basesChanged
            ? { comparison_bases: newComparisonBases }
            : undefined;

          // Only include method if non-default (US-070)
          const methodConfig = newMethod !== 'rhf'
            ? { method: newMethod as 'lda' | 'b3lyp' | 'b3lyp-d3bj' }
            : undefined;

          // Only include density config if non-default (US-062, US-063)
          const densityConfig: Record<string, unknown> = {};
          if (newShowDensity) {
            if (newDensityMode !== 'total') densityConfig.density_mode = newDensityMode;
            if (newDensityIsovalue !== 0.05) densityConfig.density_isovalue = newDensityIsovalue;
            if (newDiffIsovalue !== 0.005) densityConfig.diff_isovalue = newDiffIsovalue;
            if (newPlaneAxis !== 'xz') densityConfig.plane_axis = newPlaneAxis;
            if (newPlanePosition !== 0.0) densityConfig.plane_position = newPlanePosition;
            if (newColorScale !== 'linear') densityConfig.color_scale = newColorScale;
          }

          // Only include frequency config if non-default (US-103)
          const freqConfig: Record<string, unknown> = {};
          if (newFreqSelectedMode !== null) freqConfig.freq_mode = newFreqSelectedMode;
          if (newFreqTempK !== 298.15) freqConfig.freq_temp = newFreqTempK;
          if (newFreqPressPa !== 101325) freqConfig.freq_pres = newFreqPressPa;
          if (newFreqBroadening !== 'lorentzian') freqConfig.freq_broad = 'g';
          if (newFreqFwhm !== 8.0) freqConfig.freq_fwhm = newFreqFwhm;
          if (newFreqTab !== 'ir') freqConfig.freq_tab = newFreqTab;
          if (newFreqArrows) freqConfig.freq_arrows = true;
          if (newFreqUnits !== 'kcal_mol') {
            freqConfig.freq_units = newFreqUnits === 'hartree' ? 'ha' : 'kj';
          }

          const state: RunStateV1 = {
            schema_version: 'run_state_v1',
            app_version: APP_VERSION,
            module: 'scf',
            scf: {
              system_id: newSystemId,
              conv: newConv,
              max_iter: newMaxIter,
              diis: newDiis,
              pes: pesConfig,
              ...methodConfig,
              ...orbitalConfig,
              ...comparisonConfig,
              ...densityConfig,
              ...freqConfig,
            },
            ui: { mode: 'explain' },
          };
          updateURL(state);
        },
        URL_DEBOUNCE_MS
      ),
    []
  );

  // Update URL when parameters change (after initial load)
  useEffect(() => {
    if (!urlInitialized) return;

    debouncedUpdateURL(
      systemId, convergenceProfile, maxIterations, useDiis, method,
      pesRMin, pesRMax, pesNPoints,
      selectedOrbital, isovalue, comparisonBases,
      densityMode, densityIsovalue, diffIsovalue, planeAxis, planePosition, colorScale, showDensity,
      freqSelectedMode, freqTemperatureK, freqPressurePa, freqBroadeningKind, freqFwhmCm1, freqSpectrumTab, freqShowArrows, freqUnitsMode
    );

    return () => {
      debouncedUpdateURL.cancel();
    };
  }, [systemId, convergenceProfile, maxIterations, useDiis, method, pesRMin, pesRMax, pesNPoints, selectedOrbital, isovalue, comparisonBases, densityMode, densityIsovalue, diffIsovalue, planeAxis, planePosition, colorScale, showDensity, freqSelectedMode, freqTemperatureK, freqPressurePa, freqBroadeningKind, freqFwhmCm1, freqSpectrumTab, freqShowArrows, freqUnitsMode, urlInitialized, debouncedUpdateURL]);

  // Handle reset from invalid URL
  const handleReset = useCallback(() => {
    invalidURLRef.current = false;
    clearURL();
    reset();
  }, [reset]);

  // PES state for export and coordinate tracking (US-083)
  const pesState = useScfStore((state) => state.pesState);
  const pesInternalResult = useScfStore((state) => state.pesState.pesInternalResult);
  const pesCoordinateConfig = useScfStore((state) => state.pesState.pesCoordinateConfig);
  const trackedCoordinates = useScfStore((state) => state.pesState.trackedCoordinates);
  const setTrackedCoordinates = useScfStore((state) => state.setTrackedCoordinates);

  // PES scan geometry viewer: tracks which scan point is shown in the 3D viewer
  const [pesViewerPointIndex, setPesViewerPointIndex] = useState<number | null>(null);
  const [pesAnimating, setPesAnimating] = useState(false);
  const pesAnimationRef = useRef<ReturnType<typeof setInterval> | null>(null);

  // Reset PES viewer state when scan restarts or molecule changes
  useEffect(() => {
    setPesViewerPointIndex(null);
    setPesAnimating(false);
    if (pesAnimationRef.current) {
      clearInterval(pesAnimationRef.current);
      pesAnimationRef.current = null;
    }
  }, [systemId, pesState.scanning]);

  // Clean up animation interval on unmount
  useEffect(() => {
    return () => {
      if (pesAnimationRef.current) {
        clearInterval(pesAnimationRef.current);
      }
    };
  }, []);

  // Hover handler for PES curve: updates the viewer point on hover but does NOT
  // clear it on unhover (null). This keeps the last hovered/clicked point visible
  // in the 3D viewer even after the mouse leaves the plot area.
  const handlePesHover = useCallback((idx: number | null) => {
    if (idx !== null) {
      setPesViewerPointIndex(idx);
    }
    // Don't clear on unhover — keep the last clicked/hovered point
  }, []);

  // PES animation play/pause
  const togglePesAnimation = useCallback(() => {
    if (pesAnimating) {
      // Stop animation
      if (pesAnimationRef.current) {
        clearInterval(pesAnimationRef.current);
        pesAnimationRef.current = null;
      }
      setPesAnimating(false);
    } else {
      // Start animation
      const points = pesInternalResult?.points;
      if (!points || points.length === 0) return;
      setPesAnimating(true);
      setPesViewerPointIndex((prev) => (prev === null || prev >= points.length - 1) ? 0 : prev);
      pesAnimationRef.current = setInterval(() => {
        setPesViewerPointIndex((prev) => {
          const next = (prev ?? -1) + 1;
          if (next >= points.length) {
            // Loop back to start
            return 0;
          }
          return next;
        });
      }, 400); // 400ms per frame for smooth visualization
    }
  }, [pesAnimating, pesInternalResult]);

  // Compute 3D atoms for the selected PES scan point
  const pesAtoms3D = useMemo<Atom3D[] | null>(() => {
    if (!pesInternalResult || pesViewerPointIndex === null) return null;
    const point = pesInternalResult.points[pesViewerPointIndex];
    if (!point?.geometry || point.geometry.length === 0) return null;
    if (!workerAtoms) return null;

    // Combine Z numbers from workerAtoms with positions from scan point
    return point.geometry.map((pos, i) => ({
      symbol: ELEMENT_TO_Z_REVERSE[workerAtoms[i]?.[0] ?? 1] ?? 'X',
      position: pos as [number, number, number],
    }));
  }, [pesInternalResult, pesViewerPointIndex, workerAtoms]);

  // Handle export of results as artifact
  const handleExport = useCallback(() => {
    // Allow export if there's a successful SCF result OR PES scan results
    if (compute.status !== 'success' && pesState.results.length === 0) return;

    const result = compute.status === 'success' ? compute.result : null;

    // Create RunStateV1 for the artifact
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'scf',
      scf: {
        system_id: systemId,
        conv: convergenceProfile,
        max_iter: maxIterations,
        diis: useDiis,
      },
      ui: { mode },
    };

    // Build trace from history
    const trace = history.map((iter) => ({
      iteration: iter.iteration,
      energy: iter.energy,
      delta: iter.delta,
      diis_error: iter.diisError,
    }));

    // Include orbital energies if in internals mode and available
    const orbitalEnergies =
      mode === 'internals' && result?.orbitalEnergies
        ? result.orbitalEnergies.energies
        : undefined;

    const homoLumoGap =
      mode === 'internals' && result?.orbitalEnergies
        ? (() => {
            const energies = result.orbitalEnergies.energies;
            const nOcc = result.orbitalEnergies.nOccupied;
            if (nOcc > 0 && nOcc < energies.length) {
              return energies[nOcc] - energies[nOcc - 1];
            }
            return undefined;
          })()
        : undefined;

    // Build PES scan artifact data if results are available (US-041)
    let pesScanData: PesScanArtifactData | undefined;
    if (pesState.results.length > 0) {
      // Derive molecule name from diatomic info
      const moleculeName = diatomicInfo
        ? `${diatomicInfo.symbolA}${diatomicInfo.symbolB === diatomicInfo.symbolA ? '2' : diatomicInfo.symbolB}`
        : 'unknown';
      const basisName = diatomicInfo?.basisName ?? 'unknown';

      pesScanData = {
        molecule: moleculeName,
        basis: basisName,
        scan_parameter: 'bond_length',
        atom_indices: [0, 1],
        unit: 'bohr',
        points: pesState.results.map((p) => ({
          r_bohr: p.r,
          energy: p.energy,
          converged: p.converged,
          iterations: p.iterations,
        })),
        equilibrium: pesState.equilibrium
          ? {
              r_bohr: pesState.equilibrium.r_bohr,
              energy_hartree: pesState.equilibrium.energy_hartree,
            }
          : null,
        compute_time_ms: pesState.computeTimeMs ?? 0,
      };
    }

    // Create ScfArtifactResult from compute result
    const artifactResult: ScfArtifactResult = {
      type: 'scf',
      data: {
        energy: result?.energy ?? 0,
        converged: result?.converged ?? false,
        iterations: result?.iterations ?? 0,
        aborted: result?.aborted ?? false,
        trace,
        orbital_energies: orbitalEnergies,
        homo_lumo_gap: homoLumoGap,
        pes_scan: pesScanData,
      },
    };

    // Create and download the artifact
    const artifact = createArtifact(state, artifactResult);
    downloadArtifact(artifact);
  }, [compute, history, systemId, convergenceProfile, maxIterations, useDiis, mode, pesState, diatomicInfo]);

  // Handle import of artifact
  const handleImport = useCallback((artifact: RunArtifactV1) => {
    // Check if artifact is for the correct module
    const artifactModule = getArtifactModule(artifact);
    if (artifactModule !== 'scf') {
      console.warn(
        `Artifact is for module '${artifactModule}', but you're on Module E (SCF). ` +
        `Please navigate to the correct module to import this artifact.`
      );
      return;
    }

    // Restore the artifact
    const result = restoreArtifact(artifact);

    if (result.success) {
      // Log any warnings
      result.warnings.forEach((w) => console.warn('Import warning:', w));
    } else {
      console.error('Failed to restore artifact:', result.warnings);
    }
  }, []);

  // Show deep link error if URL was invalid
  if (invalidURLRef.current) {
    return (
      <div className="max-w-4xl mx-auto">
        <DeepLinkError onReset={handleReset} />
      </div>
    );
  }

  // Show worker error if present
  if (workerError) {
    return (
      <div className="max-w-4xl mx-auto">
        <div className="bg-red-50 border border-red-200 rounded-xl p-6">
          <h2 className="text-lg font-semibold text-red-800 mb-2">
            Worker Error
          </h2>
          <p className="text-red-600">{workerError.message}</p>
          <p className="text-sm text-red-500 mt-2">
            Please refresh the page to try again.
          </p>
        </div>
      </div>
    );
  }

  return (
    <div className="max-w-7xl mx-auto">
      {/* Header with compact toolbar */}
      <div className="mb-6">
        <div className="flex flex-wrap items-start justify-between gap-3">
          <h1 className="text-3xl font-bold text-slate-800">
            Module E: SCF Sandbox
          </h1>
          <div className="flex items-center gap-1.5">
            <ImportButton
              onImport={handleImport}
              disabled={!isReady}
            />
            <CopyLinkButton label="Share" />
            <ExportButton
              onExport={handleExport}
              disabled={compute.status !== 'success' && pesState.results.length === 0}
            />
          </div>
        </div>
        <p className="text-slate-600 mt-1">
          Interactive self-consistent field computation with optional DIIS acceleration.
          Explore how RHF and KS-DFT calculations converge for different molecular systems.
        </p>
      </div>

      {/* Worker status indicator */}
      {!isReady && (
        <div className="mb-6 bg-blue-50 border border-blue-200 rounded-lg p-4 flex items-center">
          <div className="animate-spin rounded-full h-5 w-5 border-b-2 border-blue-600 mr-3" />
          <span className="text-blue-800">Initializing compute worker...</span>
        </div>
      )}

      {/* Workflow tab bar */}
      <div className="mb-6">
        <div className="bg-slate-100 rounded-xl p-1.5 inline-flex gap-1" role="tablist" aria-label="SCF workflow">
          {WORKFLOW_TABS.map((tab) => {
            const isActive = workflowTab === tab.id;
            return (
              <button
                key={tab.id}
                type="button"
                role="tab"
                aria-selected={isActive}
                aria-controls={`panel-${tab.id}`}
                id={`tab-${tab.id}`}
                onClick={() => setWorkflowTab(tab.id)}
                title={tab.description}
                className={`px-5 py-2.5 rounded-lg text-sm font-semibold transition-all duration-150 ${
                  isActive
                    ? 'bg-white text-slate-800 shadow-sm ring-1 ring-slate-200'
                    : 'text-slate-500 hover:text-slate-700 hover:bg-white/50'
                }`}
              >
                {tab.label}
              </button>
            );
          })}
        </div>
      </div>

      {/* ================================================================
          TAB 1: SINGLE POINT — Standard SCF workflow
          ================================================================ */}
      {workflowTab === 'single-point' && (
        <div
          id="panel-single-point"
          role="tabpanel"
          aria-labelledby="tab-single-point"
          className="grid grid-cols-1 lg:grid-cols-4 gap-6 mb-6"
        >
          {/* Left sidebar: Controls and System Info */}
          <div className="lg:col-span-1">
            <div className="lg:sticky lg:top-4 lg:max-h-[calc(100vh-7rem)] lg:overflow-y-auto lg:overflow-x-hidden scrollbar-thin scrollbar-on-hover lg:pr-1 space-y-4">
              <ScfControlsPanel disabled={!isReady} workerAtoms={workerAtoms} basisName={basisName} />

              <SystemInfoPanel
                systemId={systemId}
                basisOverride={basisName}
                methodOverride={method}
                nbfOverride={sysInfoNbf}
                nelecOverride={sysInfoNelec}
                enucOverride={sysInfoEnuc}
                descriptionOverride={moleculePreset?.description}
                geometryOverride={sysInfoGeometry}
              />

              {/* Phase 3 (US-070): DFT info panel - shown when DFT method is selected */}
              {isDftMethod(method) && (
                <DftInfoPanel
                  method={method}
                  ksResult={ksResult}
                />
              )}

              {/* Phase 3 (US-062, US-063): Density visualization panel - after SCF convergence */}
              {scfResult && scfResult.converged && densityMatrix && (
                <DensityPanel
                  hasGrid={densityViewer.gridData !== null}
                  loading={densityViewer.loading}
                  error={densityViewer.error}
                  integratedDensity={densityViewer.gridData?.integratedDensity ?? null}
                  nElectronsExpected={densityViewer.gridData?.nElectronsExpected ?? null}
                  maxDensity={densityViewer.gridData?.maxDensity ?? null}
                  gridBounds={densityGridBounds}
                  diffGridData={densityViewer.diffGridData}
                />
              )}
            </div>
          </div>

          {/* Right main area: 3D viewer, results, PES */}
          <div className="lg:col-span-3 space-y-6">
            {/* Phase 2 (US-035): 3D molecular viewer - large presentation area */}
            {atoms3D && atoms3D.length > 0 && (
              <MoleculePanel
                atoms={pesAtoms3D ?? chargedAtoms3D ?? atoms3D}
                showLabels={viewerState.showLabels}
                onToggleLabels={setShowLabels}
                orbitalMeshData={orbitalViewer.meshData}
                headerControls={
                  populationResult ? (
                    <label className="flex items-center gap-1.5 text-xs text-slate-600 cursor-pointer select-none">
                      <input
                        type="checkbox"
                        checked={showChargeOverlay}
                        onChange={(e) => setShowChargeOverlay(e.target.checked)}
                        className="rounded border-slate-300 text-blue-600 focus:ring-blue-500 h-3.5 w-3.5"
                      />
                      Color by charge
                    </label>
                  ) : undefined
                }
                orbitalControls={
                  orbitalEnergies ? (
                    <>
                      <OrbitalSelector
                        energies={orbitalEnergies.energies}
                        nOccupied={orbitalEnergies.nOccupied}
                        selectedOrbital={viewerState.selectedOrbital}
                        onSelect={setSelectedOrbital}
                        disabled={orbitalViewer.loading}
                      />
                      <IsovalueSlider
                        isovalue={viewerState.isovalue}
                        onChange={setIsovalue}
                        disabled={viewerState.selectedOrbital === null}
                      />
                      {orbitalViewer.loading && (
                        <div className="flex items-center gap-2 text-xs text-slate-500">
                          <div className="animate-spin rounded-full h-3 w-3 border-b-2 border-blue-600" />
                          Computing...
                        </div>
                      )}
                      {orbitalViewer.error && (
                        <p className="text-xs text-red-600">{orbitalViewer.error}</p>
                      )}
                    </>
                  ) : atoms3D.length > 0 ? (
                    <div className="opacity-50">
                      <p className="text-xs font-semibold text-slate-400 mb-1">Orbitals</p>
                      <p className="text-xs text-slate-400">Run SCF to compute orbitals</p>
                    </div>
                  ) : undefined
                }
              >
                {/* Phase 2 (US-044): Orbital isosurface overlay — hidden when density is active */}
                {orbitalViewer.meshData && !viewerState.showDensity && (
                  <Suspense fallback={null}>
                    <LazyOrbitalSurface
                      positive={orbitalViewer.meshData.positive}
                      negative={orbitalViewer.meshData.negative}
                    />
                  </Suspense>
                )}
                {/* Phase 3 (US-062): Total density isosurface -- shown in total mode */}
                {viewerState.showDensity && viewerState.densityMode === 'total' && densityViewer.meshData && (
                  <Suspense fallback={null}>
                    <LazyDensitySurface
                      meshData={densityViewer.meshData}
                    />
                  </Suspense>
                )}
                {/* Phase 3 (US-063): Difference density dual isosurface -- shown in difference mode */}
                {viewerState.showDensity && viewerState.densityMode === 'difference' && densityViewer.diffMeshData && (
                  <Suspense fallback={null}>
                    <LazyOrbitalSurface
                      positive={densityViewer.diffMeshData.positive}
                      negative={densityViewer.diffMeshData.negative}
                      positiveColor="#22cc88"
                      positiveOpacity={0.65}
                      negativeColor="#cc4444"
                      negativeOpacity={0.35}
                    />
                  </Suspense>
                )}
              </MoleculePanel>
            )}

            {/* Phase 3 (US-062, US-063): Density cross-section plot */}
            {viewerState.showDensity && viewerState.densityMode === 'total' && densityViewer.gridData && (
              <DensityCrossSection
                gridValues={densityViewer.gridData.values}
                gridDims={densityViewer.gridData.gridDims}
                gridOrigin={densityViewer.gridData.gridOrigin}
                gridSpacing={densityViewer.gridData.gridSpacing}
                planeAxis={viewerState.planeAxis}
                planePosition={viewerState.planePosition}
                colorScale={viewerState.colorScale}
                loading={densityViewer.loading}
              />
            )}
            {/* Phase 3 (US-063): Difference density cross-section with diverging colorscale */}
            {viewerState.showDensity && viewerState.densityMode === 'difference' && densityViewer.diffGridData && (
              <DensityCrossSection
                gridValues={densityViewer.diffGridData.values}
                gridDims={densityViewer.diffGridData.gridDims}
                gridOrigin={densityViewer.diffGridData.gridOrigin}
                gridSpacing={densityViewer.diffGridData.gridSpacing}
                planeAxis={viewerState.planeAxis}
                planePosition={viewerState.planePosition}
                colorScale={viewerState.colorScale}
                loading={densityViewer.loading}
                diverging
              />
            )}

            {/* SCF Result display (includes population analysis table from US-076) */}
            <ScfResultDisplay />

            {/* Phase 3 (US-070): DFT vs HF comparison panel */}
            {isDftMethod(method) && compute.status !== 'running' && (
              <DftComparisonPanel
                onCompare={dftCompare}
                isReady={isReady}
                isRunning={isRunning || isDftComparing}
              />
            )}

          </div>
        </div>
      )}

      {/* ================================================================
          TAB 2: OPTIMIZE — Geometry optimization workflow
          ================================================================ */}
      {workflowTab === 'optimize' && (
        <div
          id="panel-optimize"
          role="tabpanel"
          aria-labelledby="tab-optimize"
          className="grid grid-cols-1 lg:grid-cols-4 gap-6 mb-6"
        >
          {/* Left sidebar: Controls + Optimize button */}
          <div className="lg:col-span-1">
            <div className="lg:sticky lg:top-4 lg:max-h-[calc(100vh-7rem)] lg:overflow-y-auto lg:overflow-x-hidden scrollbar-thin scrollbar-on-hover lg:pr-1 space-y-4">
              <ScfControlsPanel disabled={!isReady} workerAtoms={workerAtoms} basisName={basisName} hideRunButton />

              {/* Optimize Geometry button */}
              <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
                <OptimizeButton
                  onOptimize={runOptimization}
                  onCancel={cancelOptimization}
                  loading={isOptimizing}
                  disabled={!isReady || isRunning || !workerAtoms || workerAtoms.length === 0}
                />
                {optimizationState.error && (
                  <p className="mt-2 text-xs text-red-600">{optimizationState.error}</p>
                )}
              </div>

              <SystemInfoPanel
                systemId={systemId}
                basisOverride={basisName}
                methodOverride={method}
                nbfOverride={sysInfoNbf}
                nelecOverride={sysInfoNelec}
                enucOverride={sysInfoEnuc}
                descriptionOverride={moleculePreset?.description}
                geometryOverride={sysInfoGeometry}
              />

              {isDftMethod(method) && (
                <DftInfoPanel
                  method={method}
                  ksResult={ksResult}
                />
              )}
            </div>
          </div>

          {/* Right main area: 3D viewer with trajectory + optimization results */}
          <div className="lg:col-span-3 space-y-6">
            {/* 3D molecular viewer with trajectory animation */}
            {atoms3D && atoms3D.length > 0 && (
              <MoleculePanel
                atoms={trajectoryAtoms3D ?? atoms3D}
                showLabels={viewerState.showLabels}
                onToggleLabels={setShowLabels}
              >
                {/* Ghost overlay of initial geometry */}
                {ghostAtoms3D && (
                  <Suspense fallback={null}>
                    <LazyGhostAtoms atoms={ghostAtoms3D} />
                  </Suspense>
                )}
              </MoleculePanel>
            )}

            {/* Optimization progress plot */}
            {optimizationState.progress.length > 0 && (
              <OptimizationProgressPlot
                steps={optimizationState.progress}
                converged={optimizationState.result?.converged ?? null}
              />
            )}

            {/* Optimization result panel */}
            {optimizationState.result && workerAtoms && (
              <OptimizationResultPanel
                result={optimizationState.result}
                initialAtoms={workerAtoms}
                method={method}
                basisName={basisName}
                onApplyGeometry={handleApplyOptimizedGeometry}
                applied={optimizedApplied}
              />
            )}

            {/* Trajectory animation controls */}
            {optimizationState.result && optimizationState.result.steps.length > 1 && (
              <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
                <div className="flex items-center justify-between mb-2">
                  <h3 className="text-sm font-semibold text-slate-700">
                    Trajectory Animation
                  </h3>
                  <label className="flex items-center gap-1.5 text-xs text-slate-600 cursor-pointer select-none">
                    <input
                      type="checkbox"
                      checked={optimizationState.showGhostOverlay}
                      onChange={(e) => setShowGhostOverlay(e.target.checked)}
                      className="rounded border-slate-300 text-amber-600 focus:ring-amber-500 h-3.5 w-3.5"
                    />
                    Show Initial
                  </label>
                </div>
                <div className="flex items-center gap-3">
                  <span className="text-xs text-slate-500 w-10">Step {optimizationState.trajectoryStep}</span>
                  <input
                    type="range"
                    min={0}
                    max={optimizationState.result.totalSteps}
                    value={optimizationState.trajectoryStep}
                    onChange={(e) => setTrajectoryStep(Number(e.target.value))}
                    className="flex-1 h-1.5 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-amber-600"
                    aria-label="Trajectory step"
                    aria-valuemin={0}
                    aria-valuemax={optimizationState.result.totalSteps}
                    aria-valuenow={optimizationState.trajectoryStep}
                  />
                  <span className="text-xs font-mono text-slate-600 w-28 text-right">
                    {optimizationState.result.steps[
                      globalThis.Math.min(optimizationState.trajectoryStep, optimizationState.result.steps.length - 1)
                    ]?.energy.toFixed(6) ?? ''} Ha
                  </span>
                </div>
              </div>
            )}

            {/* Empty state when no optimization has been run */}
            {!isOptimizing && optimizationState.progress.length === 0 && !optimizationState.result && (
              <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 text-center">
                <div className="mx-auto w-12 h-12 rounded-full bg-amber-50 flex items-center justify-center mb-3">
                  <svg className="w-6 h-6 text-amber-500" fill="none" stroke="currentColor" viewBox="0 0 24 24" aria-hidden="true">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 7h8m0 0v8m0-8l-8 8-4-4-6 6" />
                  </svg>
                </div>
                <h3 className="text-sm font-semibold text-slate-700 mb-1">No Optimization Results</h3>
                <p className="text-xs text-slate-500 max-w-sm mx-auto">
                  Select a molecule and click &quot;Optimize Geometry&quot; in the left panel to find the minimum energy structure.
                </p>
              </div>
            )}
          </div>
        </div>
      )}

      {/* ================================================================
          TAB 3: FREQUENCY — Vibrational analysis (US-102)
          ================================================================ */}
      {workflowTab === 'frequency' && (
        <div
          id="panel-frequency"
          role="tabpanel"
          aria-labelledby="tab-frequency"
          className="grid grid-cols-1 lg:grid-cols-4 gap-6 mb-6"
        >
          <FrequencyTab
            isReady={isReady}
            workerAtoms={workerAtoms}
            atoms3D={atoms3D}
            basisName={basisName}
            method={method}
            onNavigateToOptimize={() => setWorkflowTab('optimize')}
            showLabels={viewerState.showLabels}
            onToggleLabels={setShowLabels}
          />
        </div>
      )}

      {/* ================================================================
          TAB 4: PES SCAN — Potential energy surface scanning
          ================================================================ */}
      {workflowTab === 'pes-scan' && (
        <div
          id="panel-pes-scan"
          role="tabpanel"
          aria-labelledby="tab-pes-scan"
          className="grid grid-cols-1 lg:grid-cols-4 gap-6 mb-6"
        >
          {/* Left sidebar: Controls + PES scan configuration */}
          <div className="lg:col-span-1">
            <div className="lg:sticky lg:top-4 lg:max-h-[calc(100vh-7rem)] lg:overflow-y-auto lg:overflow-x-hidden scrollbar-thin scrollbar-on-hover lg:pr-1 space-y-4">
              <ScfControlsPanel disabled={!isReady} workerAtoms={workerAtoms} basisName={basisName} hideRunButton />

              <SystemInfoPanel
                systemId={systemId}
                basisOverride={basisName}
                methodOverride={method}
                nbfOverride={sysInfoNbf}
                nelecOverride={sysInfoNelec}
                enucOverride={sysInfoEnuc}
                descriptionOverride={moleculePreset?.description}
                geometryOverride={sysInfoGeometry}
              />

              {isDftMethod(method) && (
                <DftInfoPanel
                  method={method}
                  ksResult={ksResult}
                />
              )}

              {/* PES scan controls: coordinate selector, scan mode, range, run button */}
              {workerAtoms && workerAtoms.length >= 2 && (
                <PesScanPanel
                  diatomic={diatomicInfo}
                  workerAtoms={workerAtoms}
                  basisName={basisName}
                  method={method}
                  useSpherical={useSpherical}
                  disabled={!isReady}
                />
              )}

              {/* Message when molecule has fewer than 2 atoms */}
              {(!workerAtoms || workerAtoms.length < 2) && (
                <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
                  <p className="text-xs text-slate-500">
                    Select a molecule with at least 2 atoms to enable PES scanning.
                  </p>
                </div>
              )}
            </div>
          </div>

          {/* Right main area: 3D viewer + PES curve side-by-side, then slider, then tracking */}
          <div className="lg:col-span-3 space-y-4">
            {/* Top row: 3D viewer and PES curve side-by-side */}
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
              {/* Left: 3D molecular viewer showing scan geometry */}
              {atoms3D && atoms3D.length > 0 && (
                <MoleculePanel
                  atoms={pesAtoms3D ?? atoms3D}
                  showLabels={viewerState.showLabels}
                  onToggleLabels={setShowLabels}
                />
              )}

              {/* Right: PES energy curve plot */}
              <div className="min-h-[300px]">
                {workerAtoms && workerAtoms.length >= 2 && (
                  <PesCurvePlot
                    coordinateType={(pesCoordinateConfig?.coordinateType ?? pesInternalResult?.coordinate_type ?? 'bond') as 'bond' | 'angle' | 'dihedral'}
                    onHoverPoint={handlePesHover}
                    onClickPoint={setPesViewerPointIndex}
                    selectedPointIndex={pesViewerPointIndex}
                    minHeight={300}
                  />
                )}
              </div>
            </div>

            {/* Scan Geometry Viewer: full-width slider below both panels */}
            {pesInternalResult && !pesState.scanning && pesInternalResult.points.length > 1 && (
              <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
                <div className="flex items-center justify-between mb-2">
                  <h3 className="text-sm font-semibold text-slate-700">
                    Scan Geometry Viewer
                  </h3>
                  <button
                    type="button"
                    onClick={togglePesAnimation}
                    className={`flex items-center gap-1 px-2.5 py-1 rounded-md text-xs font-medium transition-colors ${
                      pesAnimating
                        ? 'bg-amber-100 text-amber-700 hover:bg-amber-200'
                        : 'bg-blue-100 text-blue-700 hover:bg-blue-200'
                    }`}
                    aria-label={pesAnimating ? 'Pause scan animation' : 'Play scan animation'}
                  >
                    {pesAnimating ? (
                      <>
                        <svg className="w-3.5 h-3.5" fill="currentColor" viewBox="0 0 24 24" aria-hidden="true">
                          <path d="M6 4h4v16H6V4zm8 0h4v16h-4V4z" />
                        </svg>
                        Pause
                      </>
                    ) : (
                      <>
                        <svg className="w-3.5 h-3.5" fill="currentColor" viewBox="0 0 24 24" aria-hidden="true">
                          <path d="M8 5v14l11-7z" />
                        </svg>
                        Play
                      </>
                    )}
                  </button>
                </div>
                <div className="flex items-center gap-3">
                  <span className="text-xs text-slate-500 w-14 shrink-0">
                    Point {(pesViewerPointIndex ?? 0) + 1}/{pesInternalResult.points.length}
                  </span>
                  <input
                    type="range"
                    min={0}
                    max={pesInternalResult.points.length - 1}
                    value={pesViewerPointIndex ?? 0}
                    onChange={(e) => {
                      setPesViewerPointIndex(Number(e.target.value));
                      if (pesAnimating) {
                        if (pesAnimationRef.current) {
                          clearInterval(pesAnimationRef.current);
                          pesAnimationRef.current = null;
                        }
                        setPesAnimating(false);
                      }
                    }}
                    className="flex-1 h-1.5 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                    aria-label="PES scan point selector"
                    aria-valuemin={0}
                    aria-valuemax={pesInternalResult.points.length - 1}
                    aria-valuenow={pesViewerPointIndex ?? 0}
                  />
                  <span className="text-xs font-mono text-slate-600 min-w-[160px] text-right">
                    {(() => {
                      const idx = pesViewerPointIndex ?? 0;
                      const pt = pesInternalResult.points[idx];
                      if (!pt) return '';
                      const coordType = pesInternalResult.coordinate_type;
                      const coordVal = coordType === 'bond'
                        ? `R = ${pt.coordinate_value.toFixed(3)} bohr`
                        : `${(pt.coordinate_value * 180 / globalThis.Math.PI).toFixed(1)}\u00B0`;
                      return `${coordVal} | E = ${pt.energy.toFixed(6)} Ha`;
                    })()}
                  </span>
                </div>
              </div>
            )}

            {/* Coordinate tracking panel (US-083): shown after internal scan completes */}
            {pesInternalResult && !pesState.scanning && workerAtoms && (
              <CoordinateTrackingPanel
                scanPoints={pesInternalResult.points}
                coordinateType={pesInternalResult.coordinate_type as 'bond' | 'angle' | 'dihedral'}
                scannedAtomIndices={pesInternalResult.atom_indices}
                atoms={workerAtoms}
                selectedCoordinates={trackedCoordinates}
                onSelectedChange={setTrackedCoordinates}
              />
            )}

            {/* Empty state when no scan has been run and no results exist */}
            {(!pesInternalResult || pesInternalResult.points.length === 0) && !pesState.scanning && (
              <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 text-center">
                <div className="mx-auto w-12 h-12 rounded-full bg-emerald-50 flex items-center justify-center mb-3">
                  <svg className="w-6 h-6 text-emerald-500" fill="none" stroke="currentColor" viewBox="0 0 24 24" aria-hidden="true">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 7h8m0 0v8m0-8l-8 8-4-4-6 6" />
                  </svg>
                </div>
                <h3 className="text-sm font-semibold text-slate-700 mb-1">No PES Scan Results</h3>
                <p className="text-xs text-slate-500 max-w-sm mx-auto">
                  Select a molecule, choose a coordinate to scan, and click &quot;Start Scan&quot; in the left panel
                  to explore the potential energy surface.
                </p>
              </div>
            )}
          </div>
        </div>
      )}

      {/* ================================================================
          TAB 4: COMPARE — Basis set / method comparison
          ================================================================ */}
      {workflowTab === 'compare' && (
        <div
          id="panel-compare"
          role="tabpanel"
          aria-labelledby="tab-compare"
          className="grid grid-cols-1 lg:grid-cols-4 gap-6 mb-6"
        >
          {/* Left sidebar: Controls + Comparison configuration */}
          <div className="lg:col-span-1">
            <div className="lg:sticky lg:top-4 lg:max-h-[calc(100vh-7rem)] lg:overflow-y-auto lg:overflow-x-hidden scrollbar-thin scrollbar-on-hover lg:pr-1 space-y-4">
              <ScfControlsPanel disabled={!isReady} workerAtoms={workerAtoms} basisName={basisName} hideRunButton />

              {isDftMethod(method) && (
                <DftInfoPanel
                  method={method}
                  ksResult={ksResult}
                />
              )}
            </div>
          </div>

          {/* Right main area: Comparison results */}
          <div className="lg:col-span-3 space-y-6">
            {/* DFT vs HF comparison panel */}
            {isDftMethod(method) && compute.status !== 'running' && (
              <DftComparisonPanel
                onCompare={dftCompare}
                isReady={isReady}
                isRunning={isRunning || isDftComparing}
              />
            )}

            {/* Basis set comparison: selector + results */}
            {inputMode.mode === 'custom' && inputMode.geometry.atoms.length > 0 && (
              <BasisComparisonPanel
                geometry={inputMode.geometry}
                disabled={!isReady}
              />
            )}
            {inputMode.mode === 'custom' && (
              <>
                <BasisComparisonChart />
                <BasisComparisonTable />
              </>
            )}

            {/* Empty state when no comparisons configured */}
            {!isDftMethod(method) && inputMode.mode !== 'custom' && (
              <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 text-center">
                <div className="mx-auto w-12 h-12 rounded-full bg-indigo-50 flex items-center justify-center mb-3">
                  <svg className="w-6 h-6 text-indigo-500" fill="none" stroke="currentColor" viewBox="0 0 24 24" aria-hidden="true">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
                  </svg>
                </div>
                <h3 className="text-sm font-semibold text-slate-700 mb-1">No Comparisons Available</h3>
                <p className="text-xs text-slate-500 max-w-sm mx-auto">
                  To compare basis sets, switch to &quot;Custom&quot; input mode on the System tab.
                  To compare HF vs DFT, select a DFT method (LDA or B3LYP) on the System tab.
                </p>
              </div>
            )}
          </div>
        </div>
      )}

      {/* ================================================================
          MATHEMATICAL BACKGROUND (collapsible, default collapsed)
          ================================================================ */}
      <div className="bg-slate-100 rounded-xl overflow-hidden">
        <button
          type="button"
          onClick={() => setEducationOpen(!educationOpen)}
          className="w-full flex items-center justify-between p-5 text-left focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-inset focus-visible:ring-blue-500"
          aria-expanded={educationOpen}
          aria-controls="math-background-content"
        >
          <h2 className="font-semibold text-slate-700">
            Mathematical Background
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
          <div id="math-background-content" className="px-5 pb-5">
            <div className="text-slate-600 text-sm space-y-3">
              <p>
                The Hartree-Fock method seeks the ground-state electronic energy by
                optimizing a single Slater determinant wavefunction. The procedure
                is iterative because the Fock matrix depends on the density matrix,
                which depends on the molecular orbitals, which are eigenvectors of
                the Fock matrix.
              </p>

              {/* Roothaan-Hall Equation */}
              <MathBlock label="Roothaan-Hall Equation">
                {String.raw`\mathbf{F}\mathbf{C} = \mathbf{S}\mathbf{C}\boldsymbol{\varepsilon}`}
              </MathBlock>

              <p>
                The <Math>{String.raw`\mathbf{F}`}</Math> (Fock matrix) is the effective one-electron
                Hamiltonian that includes electron-electron interactions. The <Math>{String.raw`\mathbf{C}`}</Math> matrix
                contains the molecular orbital coefficients, <Math>{String.raw`\mathbf{S}`}</Math> is the overlap
                matrix, and <Math>{String.raw`\boldsymbol{\varepsilon}`}</Math> are the orbital energies.
              </p>

              {/* Fock Matrix */}
              <MathBlock label="Fock Matrix Construction">
                {String.raw`F_{\mu\nu} = H_{\mu\nu}^{\text{core}} + G_{\mu\nu}`}
              </MathBlock>

              <p>
                The two-electron contribution <Math>{String.raw`G_{\mu\nu}`}</Math> includes both
                Coulomb and exchange interactions:
              </p>

              <MathBlock label="Two-Electron Contribution">
                {String.raw`G_{\mu\nu} = \sum_{\lambda\sigma} P_{\lambda\sigma}\left[(\mu\nu|\lambda\sigma) - \frac{1}{2}(\mu\lambda|\nu\sigma)\right]`}
              </MathBlock>

              {/* Density Matrix */}
              <MathBlock label="Density Matrix">
                {String.raw`P_{\mu\nu} = 2\sum_i^{\text{occ}} C_{\mu i} C_{\nu i}`}
              </MathBlock>

              {/* Energy Expression */}
              <MathBlock label="Total Electronic Energy">
                {String.raw`E = \frac{1}{2}\text{Tr}[\mathbf{P}(\mathbf{H}^{\text{core}} + \mathbf{F})] + V_{nn}`}
              </MathBlock>

              <p>
                The Self-Consistent Field (SCF) procedure repeats until the density
                matrix (or energy) converges within a specified tolerance:
              </p>

              <ol className="list-decimal list-inside space-y-1 ml-2">
                <li>Form initial guess for density matrix <Math>{String.raw`\mathbf{P}`}</Math></li>
                <li>Build Fock matrix <Math>{String.raw`\mathbf{F}`}</Math> from <Math>{String.raw`\mathbf{P}`}</Math></li>
                <li>Solve generalized eigenvalue problem <Math>{String.raw`\mathbf{FC} = \mathbf{SC}\boldsymbol{\varepsilon}`}</Math></li>
                <li>Form new density matrix from occupied orbitals</li>
                <li>Check convergence; if not converged, go to step 2</li>
              </ol>

              {/* Convergence Criteria */}
              <MathBlock label="Convergence Criteria">
                {String.raw`\Delta E < \epsilon \quad \text{and} \quad \|\mathbf{P}_{n} - \mathbf{P}_{n-1}\| < \delta`}
              </MathBlock>

              <p>
                <strong>DIIS (Direct Inversion in the Iterative Subspace)</strong>{' '}
                accelerates convergence by extrapolating from a history of previous
                Fock matrices:
              </p>

              {/* DIIS Extrapolation */}
              <MathBlock label="DIIS Extrapolation">
                {String.raw`\mathbf{F}' = \sum_i c_i \mathbf{F}_i`}
              </MathBlock>

              <p>
                The DIIS error vector measures the deviation from self-consistency:
              </p>

              <MathBlock label="DIIS Error Vector">
                {String.raw`\mathbf{e}_i = \mathbf{F}_i \mathbf{P}_i \mathbf{S} - \mathbf{S} \mathbf{P}_i \mathbf{F}_i`}
              </MathBlock>

              <p>
                The coefficients <Math>{String.raw`c_i`}</Math> are determined by minimizing
                the norm of the extrapolated error vector subject to the
                constraint <Math>{String.raw`\sum_i c_i = 1`}</Math>.
              </p>

              {/* ---- Kohn-Sham DFT ---- */}
              <h3 className="font-semibold text-slate-700 mt-6 mb-2">Kohn-Sham DFT</h3>

              <p>
                Kohn-Sham density functional theory replaces exact exchange with an
                exchange-correlation (XC) functional of the electron density. The KS
                total energy is:
              </p>

              <MathBlock label="Kohn-Sham Energy">
                {String.raw`E_{\text{KS}} = T_s[\rho] + V_{ne}[\rho] + J[\rho] + E_{\text{xc}}[\rho] + V_{nn}`}
              </MathBlock>

              <p>
                The KS Fock matrix generalises the HF Fock matrix by replacing the
                exchange operator with the XC potential:
              </p>

              <MathBlock label="Kohn-Sham Fock Matrix">
                {String.raw`F_{\mu\nu}^{\text{KS}} = H_{\mu\nu}^{\text{core}} + J_{\mu\nu} + (1-a_x)\,V_{\mu\nu}^{\text{xc}} + a_x\,K_{\mu\nu}`}
              </MathBlock>

              <p>
                where <Math>{String.raw`a_x`}</Math> is the fraction of exact exchange
                (0 for pure DFT functionals, 0.20 for B3LYP). IQCP supports:
              </p>

              <ul className="list-disc list-inside space-y-1 ml-2">
                <li><strong>LDA</strong> — Slater exchange + VWN5 correlation (local density)</li>
                <li><strong>B3LYP</strong> — 20% exact exchange + Becke88 + LYP correlation</li>
                <li><strong>B3LYP-D3(BJ)</strong> — B3LYP with Grimme D3 dispersion (Becke-Johnson damping)</li>
              </ul>

              <p className="mt-2">
                The XC contribution is evaluated by numerical integration on a Becke
                grid (Mura-Knowles radial + Lebedev angular quadrature).
              </p>

              {/* ---- Geometry Optimization ---- */}
              <h3 className="font-semibold text-slate-700 mt-6 mb-2">Geometry Optimization</h3>

              <p>
                Given the analytical energy gradient, IQCP minimises the total energy
                with respect to nuclear positions using the L-BFGS algorithm:
              </p>

              <MathBlock label="Analytical Gradient">
                {String.raw`\frac{\partial E}{\partial R_{A,d}} = \frac{\partial V_{nn}}{\partial R_{A,d}} + \text{Tr}\!\left[\mathbf{P}\frac{\partial \mathbf{H}^{\text{core}}}{\partial R_{A,d}}\right] + \frac{1}{2}\text{Tr}\!\left[\mathbf{P}\,\mathbf{G}^{(1)}\right] - \text{Tr}\!\left[\mathbf{W}\frac{\partial \mathbf{S}}{\partial R_{A,d}}\right]`}
              </MathBlock>

              <p>
                The four terms are the nuclear repulsion gradient, one-electron
                (Hellmann-Feynman) contribution, two-electron derivative, and the
                Pulay force arising from the basis-set dependence on nuclear positions.
                L-BFGS builds an approximate inverse Hessian from gradient history,
                with trust-radius line search for stability.
              </p>

              {/* ---- Vibrational Analysis ---- */}
              <h3 className="font-semibold text-slate-700 mt-6 mb-2">Vibrational Analysis</h3>

              <p>
                The Frequency tab computes the analytical Hessian (second derivatives
                of the energy) and derives vibrational frequencies, IR/Raman spectra,
                and thermochemistry:
              </p>

              <MathBlock label="Hessian and Normal Modes">
                {String.raw`H_{A_d,B_e} = \frac{\partial^2 E}{\partial R_{A,d}\,\partial R_{B,e}} \quad\longrightarrow\quad \mathbf{H}^{\text{mw}} \mathbf{L} = \mathbf{L}\boldsymbol{\Lambda}`}
              </MathBlock>

              <p>
                The mass-weighted Hessian <Math>{String.raw`\mathbf{H}^{\text{mw}}`}</Math> is
                diagonalised after projecting out translational and rotational degrees
                of freedom. Eigenvalues <Math>{String.raw`\lambda_k`}</Math> give
                frequencies in cm<sup>-1</sup>; eigenvectors <Math>{String.raw`\mathbf{L}`}</Math> are
                the normal mode displacement vectors.
              </p>

              <p>
                The response of the wavefunction to nuclear displacements is computed
                via the Coupled-Perturbed Hartree-Fock (CPHF) equations, which also
                yield the ingredients for spectroscopic intensities:
              </p>

              <ul className="list-disc list-inside space-y-1 ml-2">
                <li>
                  <strong>IR intensity:</strong>{' '}
                  <Math>{String.raw`A_k \propto |\partial\boldsymbol{\mu}/\partial Q_k|^2`}</Math>{' '}
                  (dipole derivative, km/mol)
                </li>
                <li>
                  <strong>Raman activity:</strong>{' '}
                  <Math>{String.raw`S_k \propto \partial\boldsymbol{\alpha}/\partial Q_k`}</Math>{' '}
                  (polarizability derivative, A<sup>4</sup>/amu)
                </li>
                <li>
                  <strong>Thermochemistry:</strong> RRHO partition functions give
                  ZPE, enthalpy, entropy, and Gibbs free energy at user-specified <Math>T</Math> and <Math>P</Math>
                </li>
              </ul>

              <p className="text-xs text-slate-500 mt-6">
                References: Pulay, P. (1980). <em>Chem. Phys. Lett.</em> 73, 393-398.
                Pulay, P. (1982). <em>J. Comput. Chem.</em> 3, 556-560.
                Becke, A. D. (1993). <em>J. Chem. Phys.</em> 98, 5648-5652.
                Grimme, S. et al. (2011). <em>J. Comput. Chem.</em> 32, 1456-1465.
                See the{' '}
                <a href="/v1/reference" className="text-blue-600 hover:underline">Reference</a>{' '}
                page for full derivations and validation data.
              </p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

export default ModuleC;
