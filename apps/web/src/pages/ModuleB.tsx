/**
 * ModuleB - Rys Quadrature Lab page.
 *
 * Interactive exploration of Rys quadrature for computing roots
 * and weights used in Gaussian integral evaluation.
 *
 * @module pages/ModuleB
 */

import { useEffect, useRef, useMemo, useCallback, useState } from 'react';
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
 * Module D: Rys Quadrature Lab page.
 *
 * Provides interactive controls for exploring Rys quadrature roots and weights.
 * Integrates with the Web Worker for computation and supports deep linking.
 */
function ModuleB() {
  const { isReady, workerError } = useRys();

  // Store state
  const n = useRysStore((state) => state.n);
  const T = useRysStore((state) => state.T);
  const target = useRysStore((state) => state.target);
  const mode = useRysStore((state) => state.mode);
  const setMode = useRysStore((state) => state.setMode);
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
        // Restore display mode from URL if present
        if (urlState.ui?.mode === 'internals') {
          setMode('internals');
        }
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
          newTarget: '1e-4' | '1e-6' | '1e-8',
          newMode: string
        ) => {
          const state: RunStateV1 = {
            schema_version: 'run_state_v1',
            app_version: APP_VERSION,
            module: 'rys',
            rys: { n: newN, T: newT, target: newTarget },
            ui: { mode: newMode as 'explain' | 'internals' },
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

    debouncedUpdateURL(n, T, target, mode);

    return () => {
      debouncedUpdateURL.cancel();
    };
  }, [n, T, target, mode, urlInitialized, debouncedUpdateURL]);

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

    // Get reconstruction error from error curve data if available
    let reconstructionError = 0;
    if (errorCurve.status === 'success' && errorCurve.result) {
      const point = errorCurve.result.points.find((p) => p.n === n);
      if (point) reconstructionError = point.maxError;
    }

    // Create RysArtifactResult from compute result
    const artifactResult: RysArtifactResult = {
      type: 'rys',
      data: {
        roots: result.roots,
        weights: result.weights,
        error_curve: errorCurveData,
        reconstruction_error: reconstructionError,
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
        `Artifact is for module '${artifactModule}', but you're on Module D (Rys). ` +
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

  // Collapsible educational section state
  const [educationOpen, setEducationOpen] = useState(false);

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
      {/* Header with compact toolbar */}
      <div className="mb-6">
        <div className="flex flex-wrap items-start justify-between gap-3">
          <h1 className="text-3xl font-bold text-slate-800">
            Module D: Rys Quadrature Lab
          </h1>
          <div className="flex items-center gap-1.5">
            <ImportButton
              onImport={handleImport}
              disabled={!isReady}
            />
            <CopyLinkButton label="Share" />
            <ExportButton
              onExport={handleExport}
              disabled={compute.status !== 'success'}
            />
          </div>
        </div>
        <p className="text-slate-600 mt-1">
          Interactive exploration of Rys quadrature roots and weights for
          Gaussian integral evaluation.
        </p>
      </div>

      {/* Worker status indicator */}
      {!isReady && (
        <div className="mb-4 bg-blue-50 border border-blue-200 rounded-lg p-3 flex items-center">
          <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-600 mr-2.5" />
          <span className="text-blue-800 text-sm">Initializing compute worker...</span>
        </div>
      )}

      {/* Main content grid */}
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-6">
        {/* Controls panel */}
        <div className="lg:col-span-1">
          <RysControlsPanel disabled={!isReady} />
        </div>

        {/* Result display */}
        <div className="lg:col-span-2">
          <RysResultDisplay />
        </div>
      </div>

      {/* ================================================================
          MATHEMATICAL BACKGROUND (collapsible, default collapsed)
          ================================================================ */}
      <div className="bg-slate-100 rounded-xl overflow-hidden">
        <button
          type="button"
          onClick={() => setEducationOpen(!educationOpen)}
          className="w-full flex items-center justify-between p-5 text-left focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-inset focus-visible:ring-blue-500"
          aria-expanded={educationOpen}
          aria-controls="rys-math-background-content"
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
          <div id="rys-math-background-content" className="px-5 pb-5">
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
                The Boys function <Math>{'F_k(T)'}</Math> is proportional to the k-th moment of the Rys weight function:
              </p>
              <div className="bg-blue-50 rounded-lg p-4 border border-blue-200 mt-3">
                <MathDisplay>
                  {`\\mu_k(T) = \\int_0^1 x^k w_T(x)\\, dx = 2F_k(T)`}
                </MathDisplay>
                <p className="text-blue-700 text-xs">
                  This correspondence is the foundation of Rys quadrature: the moments of <Math>{'w_T'}</Math> are
                  directly determined by Boys function values. The implementation uses <Math>{'F_k(T)'}</Math> directly
                  as working moments (absorbing the factor of 2 into the normalization), which produces the same
                  quadrature nodes and equivalent weights.
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
                Algorithm 5.1: Moments to Nodes/Weights (Textbook Derivation)
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

              {/* Implementation Note */}
              <div className="bg-amber-50 rounded-lg p-4 border border-amber-200 mt-4">
                <h4 className="font-semibold text-amber-800 mb-2">
                  Implementation (libcint-based)
                </h4>
                <p className="text-amber-700 text-sm">
                  The textbook derivation above (Algorithm 5.1) shows the
                  Hankel/Cholesky/Jacobi pathway. The actual implementation follows
                  libcint's <strong>RDK algorithm</strong>, which is mathematically
                  equivalent but uses a different procedure:
                </p>
                <ol className="list-decimal list-inside space-y-1 ml-2 text-amber-700 text-sm mt-2">
                  <li>
                    Compute moments <Math>{'\\mu_k = F_k(T)'}</Math> via Boys function
                  </li>
                  <li>
                    <strong>Schmidt (Gram-Schmidt) orthogonalization</strong> to build
                    an orthonormal polynomial basis from the moments (equivalent to
                    implicit Cholesky factorization of the Hankel matrix)
                  </li>
                  <li>
                    Find polynomial roots via <strong>companion matrix QR</strong> decomposition
                  </li>
                  <li>
                    Compute weights from polynomial values at the roots:
                    {' '}<Math>{'w_k = 1 / \\sum_j P_j(r_k)^2'}</Math>
                  </li>
                </ol>
                <p className="text-amber-600 text-xs mt-2">
                  The Schmidt approach avoids explicitly forming the Hankel and Jacobi matrices,
                  achieving the same results with better numerical stability for the
                  orthogonalization step.
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
                        <th className="text-left py-1 px-2 text-violet-800">Moments needed</th>
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
                References: Dupuis, M., Rys, J., & King, H. F. (1976).{' '}
                <em>J. Chem. Phys.</em>, 65, 111-116.
                Implementation follows libcint (Sun, Q. <em>J. Comput. Chem.</em>, 2015, 36, 1664-1671).
              </p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

export default ModuleB;
