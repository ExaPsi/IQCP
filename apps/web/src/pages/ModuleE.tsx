/**
 * ModuleE - SCF Sandbox page.
 *
 * Interactive self-consistent field (SCF) computation with
 * optional DIIS acceleration. Run Hartree-Fock calculations
 * on preset molecular systems and visualize convergence.
 *
 * @module pages/ModuleE
 */

import { useEffect, useRef, useMemo, useCallback } from 'react';
import { ScfControlsPanel, ScfResultDisplay, SystemInfoPanel } from '../components/scf';
import { CopyLinkButton, ExportButton, ImportButton } from '../components/common';
import { Math, MathBlock } from '../components/common/Math';
import { useScf } from '../hooks/useScf';
import { useScfStore } from '../stores/scfStore';
import { getStateFromURL, hasStateInURL, updateURL, clearURL } from '../lib/url';
import { createArtifact, downloadArtifact, restoreArtifact, getArtifactModule } from '../lib/artifact';
import { DeepLinkError } from '../components/DeepLinkError';
import { DEFAULT_SCF_STATE, APP_VERSION } from '../types/run-state';
import type { RunStateV1 } from '../types/run-state';
import type { ScfArtifactResult, RunArtifactV1 } from '../types/run-artifact';

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
 * Module E: SCF Sandbox page.
 *
 * Provides interactive controls for running RHF SCF calculations
 * on preset molecular systems. Unlike Modules B/C, this uses
 * explicit Run/Cancel buttons due to longer computation times.
 */
function ModuleE() {
  const { isReady, workerError } = useScf();

  // Store state
  const systemId = useScfStore((state) => state.systemId);
  const convergenceProfile = useScfStore((state) => state.convergenceProfile);
  const maxIterations = useScfStore((state) => state.maxIterations);
  const useDiis = useScfStore((state) => state.useDiis);
  const mode = useScfStore((state) => state.mode);
  const compute = useScfStore((state) => state.compute);
  const history = useScfStore((state) => state.history);
  const urlInitialized = useScfStore((state) => state.urlInitialized);
  const initializeFromURL = useScfStore((state) => state.initializeFromURL);
  const reset = useScfStore((state) => state.reset);

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

  // Create debounced URL update function
  const debouncedUpdateURL = useMemo(
    () =>
      debounce(
        (
          newSystemId: string,
          newConv: 'loose' | 'medium' | 'tight',
          newMaxIter: number,
          newDiis: boolean
        ) => {
          const state: RunStateV1 = {
            schema_version: 'run_state_v1',
            app_version: APP_VERSION,
            module: 'scf',
            scf: {
              system_id: newSystemId,
              conv: newConv,
              max_iter: newMaxIter,
              diis: newDiis,
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

    debouncedUpdateURL(systemId, convergenceProfile, maxIterations, useDiis);

    return () => {
      debouncedUpdateURL.cancel();
    };
  }, [systemId, convergenceProfile, maxIterations, useDiis, urlInitialized, debouncedUpdateURL]);

  // Handle reset from invalid URL
  const handleReset = useCallback(() => {
    invalidURLRef.current = false;
    clearURL();
    reset();
  }, [reset]);

  // Handle export of results as artifact
  const handleExport = useCallback(() => {
    if (compute.status !== 'success') return;

    const result = compute.result;

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
      mode === 'internals' && result.orbitalEnergies
        ? result.orbitalEnergies.energies
        : undefined;

    const homoLumoGap =
      mode === 'internals' && result.orbitalEnergies
        ? (() => {
            const energies = result.orbitalEnergies.energies;
            const nOcc = result.orbitalEnergies.nOccupied;
            if (nOcc > 0 && nOcc < energies.length) {
              return energies[nOcc] - energies[nOcc - 1];
            }
            return undefined;
          })()
        : undefined;

    // Create ScfArtifactResult from compute result
    const artifactResult: ScfArtifactResult = {
      type: 'scf',
      data: {
        energy: result.energy,
        converged: result.converged,
        iterations: result.iterations,
        aborted: result.aborted,
        trace,
        orbital_energies: orbitalEnergies,
        homo_lumo_gap: homoLumoGap,
      },
    };

    // Create and download the artifact
    const artifact = createArtifact(state, artifactResult);
    downloadArtifact(artifact);
  }, [compute, history, systemId, convergenceProfile, maxIterations, useDiis, mode]);

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
      {/* Header */}
      <div className="mb-8 flex flex-wrap items-start justify-between gap-4">
        <div>
          <h1 className="text-3xl font-bold text-slate-800 mb-2">
            Module E: SCF Sandbox
          </h1>
          <p className="text-slate-600">
            Interactive self-consistent field computation with optional DIIS acceleration.
            Explore how RHF calculations converge for different molecular systems.
          </p>
        </div>
        <div className="flex gap-2">
          <ImportButton
            onImport={handleImport}
            disabled={!isReady}
          />
          <CopyLinkButton label="Share State" />
          <ExportButton
            onExport={handleExport}
            disabled={compute.status !== 'success'}
          />
        </div>
      </div>

      {/* Worker status indicator */}
      {!isReady && (
        <div className="mb-6 bg-blue-50 border border-blue-200 rounded-lg p-4 flex items-center">
          <div className="animate-spin rounded-full h-5 w-5 border-b-2 border-blue-600 mr-3" />
          <span className="text-blue-800">Initializing compute worker...</span>
        </div>
      )}

      {/* Main content grid */}
      <div className="grid grid-cols-1 lg:grid-cols-4 gap-6 mb-8">
        {/* Controls and System Info panel - independently scrollable */}
        <div className="lg:col-span-1">
          <div className="lg:sticky lg:top-4 lg:max-h-[calc(100vh-7rem)] lg:overflow-y-auto lg:overflow-x-hidden scrollbar-thin scrollbar-on-hover lg:pr-1 space-y-4">
            <ScfControlsPanel disabled={!isReady} />
            <SystemInfoPanel systemId={systemId} />
          </div>
        </div>

        {/* Result display */}
        <div className="lg:col-span-3">
          <ScfResultDisplay />
        </div>
      </div>

      {/* Mathematical background */}
      <div className="bg-slate-100 rounded-xl p-6">
        <h2 className="font-semibold text-slate-700 mb-3">
          Mathematical Background
        </h2>
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

          <p className="text-xs text-slate-500 mt-4">
            References: Pulay, P. (1980). Convergence acceleration of iterative sequences.{' '}
            <em>Chem. Phys. Lett.</em> 73, 393-398.
            Pulay, P. (1982). Improved SCF convergence acceleration.{' '}
            <em>J. Comput. Chem.</em> 3, 556-560.
          </p>
        </div>
      </div>
    </div>
  );
}

export default ModuleE;
