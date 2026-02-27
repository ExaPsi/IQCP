/**
 * ModuleA - Boys Function Lab page.
 *
 * Interactive exploration of the Boys function F_m(T) with
 * regime-based evaluation strategies.
 *
 * @module pages/ModuleA
 */

import { useEffect, useRef, useMemo, useCallback } from 'react';
import { BoysControlsPanel, BoysResultDisplay } from '../components/boys';
import { CopyLinkButton, ExportButton, ImportButton, Math as MathInline, MathBlock } from '../components/common';
import { useBoys } from '../hooks/useBoys';
import { useBoysStore } from '../stores/boysStore';
import { getStateFromURL, hasStateInURL, updateURL, clearURL } from '../lib/url';
import { createArtifact, downloadArtifact, restoreArtifact, getArtifactModule } from '../lib/artifact';
import { DeepLinkError } from '../components/DeepLinkError';
import { DEFAULT_BOYS_STATE, APP_VERSION } from '../types/run-state';
import type { RunStateV1 } from '../types/run-state';
import type { BoysArtifactResult, RunArtifactV1 } from '../types/run-artifact';

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
 * Module A: Boys Function Lab page.
 *
 * Provides interactive controls for exploring the Boys function F_m(T).
 * Integrates with the Web Worker for computation and supports deep linking.
 */
function ModuleA() {
  const { isReady, workerError } = useBoys();

  // Store state
  const m = useBoysStore((state) => state.m);
  const T = useBoysStore((state) => state.T);
  const view = useBoysStore((state) => state.view);
  const mode = useBoysStore((state) => state.mode);
  const compute = useBoysStore((state) => state.compute);
  const urlInitialized = useBoysStore((state) => state.urlInitialized);
  const initializeFromURL = useBoysStore((state) => state.initializeFromURL);
  const reset = useBoysStore((state) => state.reset);

  // Track if URL had invalid state
  const invalidURLRef = useRef(false);

  // Initialize state from URL on mount
  useEffect(() => {
    if (urlInitialized) return;

    if (hasStateInURL()) {
      const urlState = getStateFromURL();

      if (urlState && urlState.module === 'boys' && urlState.boys) {
        // Valid Boys state from URL
        initializeFromURL(urlState.boys);
      } else if (urlState && urlState.module !== 'boys') {
        // Valid state but wrong module - use defaults
        initializeFromURL(DEFAULT_BOYS_STATE.boys!);
      } else {
        // Invalid URL state
        invalidURLRef.current = true;
        initializeFromURL(DEFAULT_BOYS_STATE.boys!);
      }
    } else {
      // No URL state - use defaults
      initializeFromURL(DEFAULT_BOYS_STATE.boys!);
    }
  }, [urlInitialized, initializeFromURL]);

  // Create debounced URL update function
  const debouncedUpdateURL = useMemo(
    () =>
      debounce((newM: number, newT: number, newView: 'single' | 'sweep') => {
        const state: RunStateV1 = {
          schema_version: 'run_state_v1',
          app_version: APP_VERSION,
          module: 'boys',
          boys: { m: newM, T: newT, view: newView },
          ui: { mode: 'explain' },
        };
        updateURL(state);
      }, URL_DEBOUNCE_MS),
    []
  );

  // Update URL when parameters change (after initial load)
  useEffect(() => {
    if (!urlInitialized) return;

    debouncedUpdateURL(m, T, view);

    return () => {
      debouncedUpdateURL.cancel();
    };
  }, [m, T, view, urlInitialized, debouncedUpdateURL]);

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
      module: 'boys',
      boys: { m, T, view },
      ui: { mode },
    };

    // Create BoysArtifactResult from compute result
    const artifactResult: BoysArtifactResult = {
      type: 'boys',
      data: {
        value: result.value,
        method: result.method,
        terms_count: result.termsCount,
        estimated_error: result.estimatedError,
      },
    };

    // Create and download the artifact
    const artifact = createArtifact(state, artifactResult);
    downloadArtifact(artifact);
  }, [compute, m, T, view, mode]);

  // Handle import of artifact
  const handleImport = useCallback((artifact: RunArtifactV1) => {
    // Check if artifact is for the correct module
    const artifactModule = getArtifactModule(artifact);
    if (artifactModule !== 'boys') {
      console.warn(
        `Artifact is for module '${artifactModule}', but you're on Module A (Boys). ` +
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
            Module A: Boys Function Lab
          </h1>
          <p className="text-slate-600">
            Interactive exploration of the Boys function <MathInline>F_m(T)</MathInline> with
            regime-based evaluation strategies.
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
          <BoysControlsPanel disabled={!isReady} />
        </div>

        {/* Result display */}
        <div className="lg:col-span-2">
          <BoysResultDisplay />
        </div>
      </div>

      {/* Mathematical background */}
      <div className="bg-slate-100 rounded-xl p-6">
        <h2 className="font-semibold text-slate-700 mb-3">
          Mathematical Background
        </h2>
        <div className="text-slate-600 text-sm space-y-3">
          <p>
            The Boys function is defined as:
          </p>
          <MathBlock label="Boys Function Definition">
            {`F_m(T) = \\int_0^1 t^{2m} e^{-Tt^2} \\, dt`}
          </MathBlock>
          <p>
            This integral is fundamental to the evaluation of Gaussian integrals
            in quantum chemistry. The parameter <MathInline>{`T = \\rho R_{PQ}^2`}</MathInline> measures
            the effective squared distance between two Gaussian overlap distributions.
          </p>

          {/* 3 Theoretical Regimes */}
          <h3 className="font-semibold text-slate-700 mt-4 mb-2">
            Three Computational Regimes (Theory)
          </h3>
          <p>
            Different computational strategies are optimal depending on the value of <MathInline>T</MathInline>:
          </p>

          {/* Small T - Series */}
          <div className="bg-blue-50 rounded-lg p-4 border border-blue-200 mt-3">
            <h4 className="font-semibold text-blue-800 mb-2">
              1. Small T (<MathInline>{`T < 25`}</MathInline>): Series Expansion
            </h4>
            <MathBlock className="!bg-white !border-blue-100 !my-2">
              {`F_m(T) = \\sum_{k=0}^{\\infty} \\frac{(-T)^k}{k! \\cdot (2m+2k+1)}`}
            </MathBlock>
            <p className="text-blue-700 text-xs">
              The Taylor series expansion converges rapidly for small <MathInline>T</MathInline>.
              Truncating after ~15-25 terms typically achieves machine precision.
            </p>
          </div>

          {/* Moderate T - erf + Recurrence */}
          <div className="bg-green-50 rounded-lg p-4 border border-green-200 mt-3">
            <h4 className="font-semibold text-green-800 mb-2">
              2. Moderate T (<MathInline>{`25 \\leq T < 30+5m`}</MathInline>): erf + Upward Recurrence
            </h4>
            <MathBlock className="!bg-white !border-green-100 !my-2">
              {`F_0(T) = \\frac{1}{2}\\sqrt{\\frac{\\pi}{T}} \\cdot \\text{erf}(\\sqrt{T})`}
            </MathBlock>
            <MathBlock className="!bg-white !border-green-100 !my-2">
              {`F_{m+1}(T) = \\frac{(2m+1) \\cdot F_m(T) - e^{-T}}{2T}`}
            </MathBlock>
            <p className="text-green-700 text-xs">
              For moderate <MathInline>T</MathInline>, compute <MathInline>{`F_0(T)`}</MathInline> via the error function,
              then use upward recurrence. Cancellation is manageable in this regime.
            </p>
          </div>

          {/* Large T - Asymptotic */}
          <div className="bg-violet-50 rounded-lg p-4 border border-violet-200 mt-3">
            <h4 className="font-semibold text-violet-800 mb-2">
              3. Large T (<MathInline>{`T > 30+5m`}</MathInline>): Asymptotic Expansion
            </h4>
            <MathBlock className="!bg-white !border-violet-100 !my-2">
              {`F_m(T) \\sim \\frac{(2m-1)!!}{2^{m+1}} \\sqrt{\\frac{\\pi}{T^{2m+1}}}`}
            </MathBlock>
            <p className="text-violet-700 text-xs">
              For large <MathInline>T</MathInline>, the asymptotic expansion achieves machine precision.
              The threshold <MathInline>{`30+5m`}</MathInline> depends on the order <MathInline>m</MathInline> (e.g., <MathInline>{`T>30`}</MathInline> for <MathInline>{`m=0`}</MathInline>, <MathInline>{`T>55`}</MathInline> for <MathInline>{`m=5`}</MathInline>).
            </p>
          </div>

          {/* Implementation Note */}
          <div className="bg-amber-50 rounded-lg p-4 border border-amber-200 mt-4">
            <h4 className="font-semibold text-amber-800 mb-2">
              Implementation (libcint-based)
            </h4>
            <p className="text-amber-700 text-sm">
              The actual implementation uses <strong>2 methods</strong> with an
              <MathInline>m</MathInline>-dependent turnover point optimized for double precision:
            </p>
            <ul className="list-disc list-inside space-y-1 ml-2 text-amber-700 text-sm mt-2">
              <li>
                <strong>Series:</strong> <MathInline>{`T < \\text{turnover}(m)`}</MathInline> with downward recurrence
              </li>
              <li>
                <strong>Recurrence:</strong> <MathInline>{`T \\geq \\text{turnover}(m)`}</MathInline> using erf + upward recurrence
                (combines moderate and large T regimes)
              </li>
            </ul>
            <p className="text-amber-600 text-xs mt-2">
              Turnover examples: <MathInline>{`m=0,1 \\to 0`}</MathInline> (always recurrence), <MathInline>{`m=5 \\to {\\sim}2.1`}</MathInline>, <MathInline>{`m=10 \\to {\\sim}4.05`}</MathInline>
            </p>
          </div>

          <p className="text-xs text-slate-500 mt-4">
            Reference: Shavitt, I. (1963). <em>Methods in Computational Physics</em>, 2, 1-45.
            See also lecture notes ch04_s06 for detailed regime analysis.
          </p>
        </div>
      </div>
    </div>
  );
}

export default ModuleA;
