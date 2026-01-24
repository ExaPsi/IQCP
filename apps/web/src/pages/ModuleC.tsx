/**
 * ModuleC - Rys Quadrature Lab page.
 *
 * Interactive exploration of Rys quadrature for computing roots
 * and weights used in Gaussian integral evaluation.
 *
 * @module pages/ModuleC
 */

import { useEffect, useRef, useMemo, useCallback } from 'react';
import { RysControlsPanel, RysResultDisplay } from '../components/rys';
import { Math, MathBlock, MathDisplay } from '../components/common/Math';
import { CopyLinkButton, ExportButton, ImportButton } from '../components/common';
import { useRys } from '../hooks/useRys';
import { useRysStore } from '../stores/rysStore';
import { getStateFromURL, hasStateInURL, updateURL, clearURL } from '../lib/url';
import { createArtifact, downloadArtifact, restoreArtifact, getArtifactModule } from '../lib/artifact';
import { DeepLinkError } from '../components/DeepLinkError';
import { DEFAULT_RYS_STATE, APP_VERSION } from '../types/run-state';
import type { RunStateV1 } from '../types/run-state';
import type { RysArtifactResult, RunArtifactV1 } from '../types/run-artifact';

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
 * Longer than compute debounce to avoid excessive URL changes.
 */
const URL_DEBOUNCE_MS = 300;

/**
 * Module C: Rys Quadrature Lab page.
 *
 * Provides interactive controls for exploring Rys quadrature roots and weights.
 * Integrates with the Web Worker for computation and supports deep linking.
 */
function ModuleC() {
  const { isReady, workerError } = useRys();

  // Store state
  const n = useRysStore((state) => state.n);
  const T = useRysStore((state) => state.T);
  const target = useRysStore((state) => state.target);
  const mode = useRysStore((state) => state.mode);
  const compute = useRysStore((state) => state.compute);
  const errorCurve = useRysStore((state) => state.errorCurve);
  const urlInitialized = useRysStore((state) => state.urlInitialized);
  const initializeFromURL = useRysStore((state) => state.initializeFromURL);
  const reset = useRysStore((state) => state.reset);

  // Track if URL had invalid state
  const invalidURLRef = useRef(false);

  // Initialize state from URL on mount
  useEffect(() => {
    if (urlInitialized) return;

    if (hasStateInURL()) {
      const urlState = getStateFromURL();

      if (urlState && urlState.module === 'rys' && urlState.rys) {
        // Valid Rys state from URL
        initializeFromURL(urlState.rys);
      } else if (urlState && urlState.module !== 'rys') {
        // Valid state but wrong module - use defaults
        initializeFromURL(DEFAULT_RYS_STATE.rys!);
      } else {
        // Invalid URL state
        invalidURLRef.current = true;
        initializeFromURL(DEFAULT_RYS_STATE.rys!);
      }
    } else {
      // No URL state - use defaults
      initializeFromURL(DEFAULT_RYS_STATE.rys!);
    }
  }, [urlInitialized, initializeFromURL]);

  // Create debounced URL update function
  const debouncedUpdateURL = useMemo(
    () =>
      debounce(
        (
          newN: number,
          newT: number,
          newTarget: '1e-4' | '1e-6' | '1e-8'
        ) => {
          const state: RunStateV1 = {
            schema_version: 'run_state_v1',
            app_version: APP_VERSION,
            module: 'rys',
            rys: { n: newN, T: newT, target: newTarget },
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

    debouncedUpdateURL(n, T, target);

    return () => {
      debouncedUpdateURL.cancel();
    };
  }, [n, T, target, urlInitialized, debouncedUpdateURL]);

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
      module: 'rys',
      rys: { n, T, target },
      ui: { mode },
    };

    // Build error curve data if available
    const errorCurveData =
      errorCurve.status === 'success'
        ? errorCurve.result.points.map((point) => ({
            n: point.n,
            max_error: point.maxError,
          }))
        : undefined;

    // Create RysArtifactResult from compute result
    const artifactResult: RysArtifactResult = {
      type: 'rys',
      data: {
        roots: result.roots,
        weights: result.weights,
        error_curve: errorCurveData,
        reconstruction_error: 0, // Not computed in current implementation
      },
    };

    // Create and download the artifact
    const artifact = createArtifact(state, artifactResult);
    downloadArtifact(artifact);
  }, [compute, errorCurve, n, T, target, mode]);

  // Handle import of artifact
  const handleImport = useCallback((artifact: RunArtifactV1) => {
    // Check if artifact is for the correct module
    const artifactModule = getArtifactModule(artifact);
    if (artifactModule !== 'rys') {
      console.warn(
        `Artifact is for module '${artifactModule}', but you're on Module C (Rys). ` +
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
    <div className="max-w-6xl mx-auto">
      {/* Header */}
      <div className="mb-8 flex flex-wrap items-start justify-between gap-4">
        <div>
          <h1 className="text-3xl font-bold text-slate-800 mb-2">
            Module C: Rys Quadrature Lab
          </h1>
          <p className="text-slate-600">
            Interactive exploration of Rys quadrature roots and weights for
            Gaussian integral evaluation.
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
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-8">
        {/* Controls panel */}
        <div className="lg:col-span-1">
          <RysControlsPanel disabled={!isReady} />
        </div>

        {/* Result display */}
        <div className="lg:col-span-2">
          <RysResultDisplay />
        </div>
      </div>

      {/* Mathematical background */}
      <div className="bg-slate-100 rounded-xl p-6">
        <h2 className="font-semibold text-slate-700 mb-3">
          Mathematical Background
        </h2>
        <div className="text-slate-600 text-sm space-y-3">
          <p>
            Rys quadrature provides specialized Gaussian quadrature rules for
            integrals involving the Rys weight function:
          </p>
          <MathBlock label="Rys Weight Function">
            {`w_T(x) = x^{-1/2} e^{-Tx}, \\quad x \\in [0,1]`}
          </MathBlock>

          {/* Main Quadrature Formula */}
          <h3 className="font-semibold text-slate-700 mt-4 mb-2">
            Quadrature Approximation
          </h3>
          <p>
            The quadrature rule approximates integrals with the Rys weight function:
          </p>
          <MathBlock label="Quadrature Formula">
            {`\\int_0^1 f(t^2) e^{-Tt^2} dt \\approx \\sum_{i=1}^{n} w_i f(t_i^2)`}
          </MathBlock>

          {/* Boys-to-Moments Correspondence */}
          <h3 className="font-semibold text-slate-700 mt-4 mb-2">
            Boys Functions as Moments
          </h3>
          <p>
            The Boys function <Math>{'F_n(T)'}</Math> equals half the n-th moment of the Rys weight function:
          </p>
          <div className="bg-blue-50 rounded-lg p-4 border border-blue-200 mt-3">
            <MathDisplay>
              {`\\mu_k(T) = \\int_0^1 x^k w_T(x)\\, dx = \\mathbf{2F_k(T)}`}
            </MathDisplay>
            <p className="text-blue-700 text-xs">
              This correspondence is the foundation of Rys quadrature: computing Boys functions is equivalent to computing moments of <Math>{'w_T'}</Math>.
            </p>
          </div>

          {/* Hankel Matrix */}
          <h3 className="font-semibold text-slate-700 mt-4 mb-2">
            Hankel Matrix Structure
          </h3>
          <p>
            The moments form a symmetric Hankel matrix <Math>{'H'}</Math> where <Math>{'H_{ij} = \\mu_{i+j}'}</Math>:
          </p>
          <div className="bg-white rounded-lg p-4 border border-slate-200 mt-2 overflow-x-auto">
            <MathDisplay>
              {`H = \\begin{pmatrix} \\mu_0 & \\mu_1 & \\mu_2 \\\\ \\mu_1 & \\mu_2 & \\mu_3 \\\\ \\mu_2 & \\mu_3 & \\mu_4 \\end{pmatrix}`}
            </MathDisplay>
            <p className="text-xs text-slate-500 mt-2 text-center">Example: 3x3 Hankel matrix for <Math>{'n_r = 3'}</Math></p>
          </div>

          {/* Algorithm 5.1 Pipeline */}
          <h3 className="font-semibold text-slate-700 mt-4 mb-2">
            Algorithm 5.1: Moments to Nodes/Weights
          </h3>
          <div className="bg-gradient-to-r from-blue-50 via-green-50 to-purple-50 rounded-lg p-4 border border-slate-200">
            <div className="flex flex-wrap items-center justify-center gap-2 text-xs font-medium">
              <span className="bg-blue-100 text-blue-800 px-3 py-1.5 rounded-full border border-blue-200">
                1. <Math>{'\\mu_k = 2F_k(T)'}</Math>
              </span>
              <span className="text-slate-400">then</span>
              <span className="bg-blue-100 text-blue-800 px-3 py-1.5 rounded-full border border-blue-200">
                2. <Math>{'H, H^{(1)}'}</Math>
              </span>
              <span className="text-slate-400">then</span>
              <span className="bg-green-100 text-green-800 px-3 py-1.5 rounded-full border border-green-200">
                3. <Math>{'H = LL^T'}</Math>
              </span>
              <span className="text-slate-400">then</span>
              <span className="bg-green-100 text-green-800 px-3 py-1.5 rounded-full border border-green-200">
                4. <Math>{'J = L^{-1} H^{(1)} L^{-T}'}</Math>
              </span>
              <span className="text-slate-400">then</span>
              <span className="bg-purple-100 text-purple-800 px-3 py-1.5 rounded-full border border-purple-200">
                5. Eigendecompose <Math>{'J'}</Math>
              </span>
              <span className="text-slate-400">then</span>
              <span className="bg-purple-100 text-purple-800 px-3 py-1.5 rounded-full border border-purple-200">
                6. Nodes/Weights
              </span>
            </div>
            <p className="text-xs text-slate-600 mt-3 text-center">
              <Math>{'C = L^{-1}'}</Math> transforms monomials to orthonormal polynomials (Gram-Schmidt via Cholesky)
            </p>
          </div>

          {/* Golub-Welsch Theorem */}
          <h3 className="font-semibold text-slate-700 mt-4 mb-2">
            Golub-Welsch Theorem
          </h3>
          <div className="bg-green-50 rounded-lg p-4 border border-green-200">
            <MathDisplay>
              {`t_i = \\text{eigenvalues of } J, \\quad w_i = \\mu_0 \\cdot (V_{0i})^2`}
            </MathDisplay>
            <p className="text-green-700 text-xs">
              <Math>{'V_{0i}'}</Math> is the first component of the normalized eigenvector for eigenvalue <Math>{'t_i'}</Math>.
              The weights are always positive because <Math>{'\\mu_0 > 0'}</Math> and <Math>{'(V_{0i})^2 > 0'}</Math>.
            </p>
          </div>

          {/* Root Count Rule */}
          <h3 className="font-semibold text-slate-700 mt-4 mb-2">
            Root Count Rule for ERIs
          </h3>
          <div className="bg-violet-50 rounded-lg p-4 border border-violet-200">
            <MathDisplay>
              {`n_r = \\left\\lfloor \\frac{L}{2} \\right\\rfloor + 1, \\quad \\text{where } L = l_a + l_b + l_c + l_d`}
            </MathDisplay>
            <p className="text-violet-700 text-xs mb-2">
              For a shell quartet <Math>{'(l_a l_b | l_c l_d)'}</Math> with total angular momentum <Math>{'L'}</Math>, an <Math>{'n_r'}</Math>-point quadrature is exact for polynomials of degree <Math>{'\\leq 2n_r-1'}</Math>.
            </p>
            <div className="overflow-x-auto">
              <table className="w-full text-xs border-collapse">
                <thead>
                  <tr className="border-b border-violet-200">
                    <th className="text-left py-1 px-2 text-violet-800">Shell Quartet</th>
                    <th className="text-center py-1 px-2 text-violet-800"><Math>{'L'}</Math></th>
                    <th className="text-center py-1 px-2 text-violet-800"><Math>{'n_r'}</Math></th>
                    <th className="text-left py-1 px-2 text-violet-800">Boys orders</th>
                  </tr>
                </thead>
                <tbody className="text-violet-700">
                  <tr><td className="py-1 px-2"><Math>{'(ss|ss)'}</Math></td><td className="text-center py-1 px-2">0</td><td className="text-center py-1 px-2">1</td><td className="py-1 px-2"><Math>{'F_0'}</Math></td></tr>
                  <tr><td className="py-1 px-2"><Math>{'(ps|ss)'}</Math></td><td className="text-center py-1 px-2">1</td><td className="text-center py-1 px-2">1</td><td className="py-1 px-2"><Math>{'F_0, F_1'}</Math></td></tr>
                  <tr><td className="py-1 px-2"><Math>{'(pp|ss)'}</Math></td><td className="text-center py-1 px-2">2</td><td className="text-center py-1 px-2">2</td><td className="py-1 px-2"><Math>{'F_0 \\ldots F_3'}</Math></td></tr>
                  <tr><td className="py-1 px-2"><Math>{'(pp|pp)'}</Math></td><td className="text-center py-1 px-2">4</td><td className="text-center py-1 px-2">3</td><td className="py-1 px-2"><Math>{'F_0 \\ldots F_5'}</Math></td></tr>
                  <tr><td className="py-1 px-2"><Math>{'(dd|pp)'}</Math></td><td className="text-center py-1 px-2">6</td><td className="text-center py-1 px-2">4</td><td className="py-1 px-2"><Math>{'F_0 \\ldots F_7'}</Math></td></tr>
                </tbody>
              </table>
            </div>
          </div>

          <p className="text-xs text-slate-500 mt-4">
            Reference: Dupuis, M., Rys, J., & King, H. F. (1976).{' '}
            <em>J. Chem. Phys.</em>, 65, 111-116. See also lecture notes Chapter 5 for Algorithm 5.1.
          </p>
        </div>
      </div>
    </div>
  );
}

export default ModuleC;
