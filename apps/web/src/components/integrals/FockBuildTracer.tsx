/**
 * FockBuildTracer - Step-by-step visualization of Fock matrix construction.
 *
 * Provides a three-step guided exploration of how the Fock matrix is built:
 *   Step 1: H^core (one-electron core Hamiltonian)
 *   Step 2: G(P) = J - 0.5*K (two-electron contribution with Coulomb/Exchange decomposition)
 *   Step 3: F = H^core + G (complete Fock matrix)
 *
 * Includes step navigation, an animate mode that auto-advances through steps,
 * a self-consistency loop annotation, and element-level drill-down via
 * FockElementDetail when an element is clicked in step 3.
 *
 * @module components/integrals/FockBuildTracer
 * @see US-058 Fock Build Tracing
 * @see Szabo & Ostlund (1996), Eq. 3.154
 */

import React, { useCallback, useEffect, useRef, useState, memo } from 'react';
import type {
  FockDecompositionResult,
  FockBuildStep,
  FockStep2View,
} from '../../stores/integralStore';
import { FockStepView } from './FockStepView';
import { FockElementDetail } from './FockElementDetail';

// ============================================================================
// Types
// ============================================================================

export interface FockBuildTracerProps {
  /** Fock decomposition result data */
  decomposition: FockDecompositionResult;
  /** Current step (1, 2, or 3) */
  step: FockBuildStep;
  /** Step change handler */
  onStepChange: (step: FockBuildStep) => void;
  /** Sub-view within step 2 (G, J, or K) */
  step2View: FockStep2View;
  /** Step 2 sub-view change handler */
  onStep2ViewChange: (view: FockStep2View) => void;
  /** Selected element [row, col] for drill-down */
  selectedElement: [number, number] | null;
  /** Element selection handler */
  onElementSelect: (element: [number, number] | null) => void;
}

// ============================================================================
// Step Configuration
// ============================================================================

interface StepConfig {
  step: FockBuildStep;
  label: React.ReactNode;
  shortLabel: string;
  description: string;
}

const STEPS: StepConfig[] = [
  {
    step: 1,
    label: <span>H<sup>core</sup></span>,
    shortLabel: '1',
    description: 'Core Hamiltonian',
  },
  {
    step: 2,
    label: 'G(P)',
    shortLabel: '2',
    description: 'Two-Electron Part',
  },
  {
    step: 3,
    label: <span>F = H<sup>core</sup> + G</span>,
    shortLabel: '3',
    description: 'Fock Matrix',
  },
];

/** Delay between steps in animate mode (ms) */
const ANIMATE_DELAY_MS = 1500;

// ============================================================================
// Self-Consistency Annotation
// ============================================================================

/**
 * Visual annotation showing the SCF self-consistency loop.
 * F -> C -> P -> F circular dependency.
 */
const SelfConsistencyAnnotation = memo(function SelfConsistencyAnnotation() {
  return (
    <div className="bg-amber-50 rounded-lg px-4 py-3 border border-amber-200">
      <h4 className="text-xs font-semibold text-amber-800 uppercase tracking-wider mb-2">
        Self-Consistency Loop
      </h4>

      {/* Cycle diagram */}
      <div className="flex items-center justify-center gap-1 text-xs font-mono text-amber-900 mb-3">
        <span className="bg-amber-100 px-2 py-1 rounded font-semibold">F</span>
        <span className="text-amber-600">&rarr;</span>
        <span className="text-amber-700">diag(F)</span>
        <span className="text-amber-600">&rarr;</span>
        <span className="bg-amber-100 px-2 py-1 rounded font-semibold">C</span>
        <span className="text-amber-600">&rarr;</span>
        <span className="text-amber-700">P = 2CC&#x1D40;</span>
        <span className="text-amber-600">&rarr;</span>
        <span className="text-amber-700">G(P)</span>
        <span className="text-amber-600">&rarr;</span>
        <span className="bg-amber-100 px-2 py-1 rounded font-semibold">F</span>
      </div>

      <p className="text-[11px] text-amber-800 leading-relaxed">
        The Fock matrix <span className="font-semibold">F</span> depends on the
        density matrix <span className="font-semibold">P</span>, which depends on
        the MO coefficients <span className="font-semibold">C</span>, which are
        obtained by diagonalizing <span className="font-semibold">F</span>.
        This circular dependency is why the Hartree-Fock method requires
        iterative self-consistent field (SCF) convergence.
      </p>
    </div>
  );
});

// ============================================================================
// Component
// ============================================================================

/**
 * Interactive step-by-step Fock matrix construction visualization.
 *
 * Features:
 * - Three-step navigation (H^core -> G(P) -> F)
 * - Animate mode with auto-advancing steps
 * - J/K/G sub-view toggle in step 2
 * - Element click for detailed breakdown (step 3)
 * - Self-consistency loop annotation
 * - Keyboard navigation (left/right arrows for steps)
 */
export const FockBuildTracer = memo(function FockBuildTracer({
  decomposition,
  step,
  onStepChange,
  step2View,
  onStep2ViewChange,
  selectedElement,
  onElementSelect,
}: FockBuildTracerProps) {
  const [isAnimating, setIsAnimating] = useState(false);
  const animateTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  // Clean up animation timer on unmount
  useEffect(() => {
    return () => {
      if (animateTimerRef.current) {
        clearTimeout(animateTimerRef.current);
      }
    };
  }, []);

  // Animation logic: advance to next step after delay
  useEffect(() => {
    if (!isAnimating) return;

    if (step < 3) {
      animateTimerRef.current = setTimeout(() => {
        onStepChange((step + 1) as FockBuildStep);
      }, ANIMATE_DELAY_MS);
    } else {
      // Reached the end, stop animating
      setIsAnimating(false);
    }

    return () => {
      if (animateTimerRef.current) {
        clearTimeout(animateTimerRef.current);
      }
    };
  }, [isAnimating, step, onStepChange]);

  const handleAnimate = useCallback(() => {
    if (isAnimating) {
      setIsAnimating(false);
      return;
    }
    // Start from step 1
    onStepChange(1);
    onElementSelect(null);
    setIsAnimating(true);
  }, [isAnimating, onStepChange, onElementSelect]);

  const handleStepClick = useCallback(
    (newStep: FockBuildStep) => {
      setIsAnimating(false);
      onStepChange(newStep);
      onElementSelect(null);
    },
    [onStepChange, onElementSelect],
  );

  // Keyboard navigation
  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if (e.key === 'ArrowRight' && step < 3) {
        e.preventDefault();
        handleStepClick((step + 1) as FockBuildStep);
      } else if (e.key === 'ArrowLeft' && step > 1) {
        e.preventDefault();
        handleStepClick((step - 1) as FockBuildStep);
      }
    },
    [step, handleStepClick],
  );

  return (
    <div
      className="space-y-4"
      onKeyDown={handleKeyDown}
      tabIndex={0}
      role="region"
      aria-label="Fock matrix build tracer"
    >
      {/* Step navigation bar */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
        <div className="flex items-center justify-between mb-3">
          <h3 className="text-sm font-semibold text-slate-700">
            Fock Matrix Construction
          </h3>
          <button
            onClick={handleAnimate}
            className={`px-3 py-1.5 text-xs font-medium rounded-lg border transition-colors
              ${
                isAnimating
                  ? 'bg-amber-100 text-amber-800 border-amber-300'
                  : 'text-slate-600 border-slate-300 hover:bg-slate-50'
              }`}
            aria-label={isAnimating ? 'Stop animation' : 'Animate build steps'}
          >
            {isAnimating ? 'Stop' : 'Animate'}
          </button>
        </div>

        {/* Step buttons */}
        <div className="flex gap-2">
          {STEPS.map((s) => {
            const isActive = step === s.step;
            const isCompleted = step > s.step;
            return (
              <button
                key={s.step}
                onClick={() => handleStepClick(s.step)}
                className={`flex-1 px-3 py-2 text-sm rounded-lg border transition-colors
                  focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500 focus-visible:ring-offset-1
                  ${
                    isActive
                      ? 'bg-primary-100 text-primary-800 border-primary-300 font-semibold'
                      : isCompleted
                        ? 'bg-emerald-50 text-emerald-700 border-emerald-200'
                        : 'text-slate-500 border-slate-200 hover:bg-slate-50'
                  }`}
                aria-pressed={isActive}
                aria-label={`Step ${s.step}: ${s.description}`}
              >
                <span className="block text-xs text-slate-400 mb-0.5">
                  Step {s.shortLabel}
                </span>
                <span className="block">{s.label}</span>
              </button>
            );
          })}
        </div>
      </div>

      {/* Step 2 sub-view toggle */}
      {step === 2 && (
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-3">
          <div className="flex gap-1.5">
            {(['G', 'J', 'K'] as FockStep2View[]).map((view) => {
              const isActive = step2View === view;
              const viewLabels: Record<FockStep2View, { label: string; full: string }> = {
                G: { label: 'G = J - 0.5K', full: 'Two-electron contribution' },
                J: { label: 'J (Coulomb)', full: 'Coulomb matrix' },
                K: { label: 'K (Exchange)', full: 'Exchange matrix' },
              };
              const config = viewLabels[view];
              return (
                <button
                  key={view}
                  onClick={() => onStep2ViewChange(view)}
                  className={`flex-1 px-2 py-1.5 text-xs rounded-lg border transition-colors
                    ${
                      isActive
                        ? 'bg-primary-100 text-primary-800 border-primary-300 font-semibold'
                        : 'text-slate-500 border-slate-200 hover:bg-slate-50'
                    }`}
                  title={config.full}
                  aria-pressed={isActive}
                >
                  {config.label}
                </button>
              );
            })}
          </div>
        </div>
      )}

      {/* Matrix visualization for current step */}
      <FockStepView
        decomposition={decomposition}
        step={step}
        step2View={step2View}
        selectedElement={selectedElement}
        onElementSelect={onElementSelect}
      />

      {/* Element detail drill-down (only in step 3 when an element is selected) */}
      {step === 3 && selectedElement && (
        <FockElementDetail
          decomposition={decomposition}
          row={selectedElement[0]}
          col={selectedElement[1]}
          onClose={() => onElementSelect(null)}
        />
      )}

      {/* Self-consistency annotation */}
      <SelfConsistencyAnnotation />
    </div>
  );
});
