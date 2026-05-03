/**
 * DiscussionPrompt - Educational discussion prompt for basis set comparison.
 *
 * Addresses the common student misconception that "more basis functions always
 * means more accurate results." Uses Socratic guided questions to lead students
 * to discover the variational principle, diminishing returns, and computational
 * cost tradeoffs.
 *
 * Styled as an amber/yellow callout box consistent with educational notes
 * elsewhere in the application. The prompt is collapsible (starts expanded)
 * so students can expand it when ready to discuss.
 *
 * @module components/basisExplorer/DiscussionPrompt
 * @see US-052 AC6 - Guided discussion prompt
 */

import { useState, memo } from 'react';
import { ELEMENTS } from './basisData';

// ============================================================================
// Types
// ============================================================================

export interface DiscussionPromptProps {
  /** Basis sets being compared (for contextual references) */
  comparedBases: string[];
  /** Shell type being compared (e.g., "1s", "2p") */
  shellType: string;
  /** Atomic number of the element being examined */
  atomicNumber: number;
}

// ============================================================================
// Display Label Mapping
// ============================================================================

const BASIS_LABELS: Record<string, string> = {
  'sto-3g': 'STO-3G',
  '3-21g': '3-21G',
  '6-31g': '6-31G',
  '6-31g*': '6-31G*',
  '6-31+g*': '6-31+G*',
  'cc-pvdz': 'cc-pVDZ',
};

// ============================================================================
// Component
// ============================================================================

/**
 * Educational discussion prompt below the comparison plot.
 *
 * Content is responsive to the specific basis sets being compared,
 * mentioning them by name. Uses Socratic questioning (observation,
 * structure, critical thinking) rather than lecturing.
 */
export const DiscussionPrompt = memo(function DiscussionPrompt({
  comparedBases,
  shellType,
  atomicNumber,
}: DiscussionPromptProps) {
  const [isExpanded, setIsExpanded] = useState(false);

  const elementSymbol = ELEMENTS.find((e) => e.z === atomicNumber)?.symbol ?? `Z=${atomicNumber}`;
  const basisLabels = comparedBases.map((b) => BASIS_LABELS[b] ?? b.toUpperCase());
  const basisListText = formatBasisList(basisLabels);

  // Determine which guided questions are relevant
  const hasSplitValence = comparedBases.some(
    (b) => b === '3-21g' || b === '6-31g' || b === '6-31g*' || b === '6-31+g*' || b === 'cc-pvdz',
  );
  const hasMinimal = comparedBases.includes('sto-3g');
  const hasPolarization = comparedBases.some((b) => b === '6-31g*' || b === '6-31+g*' || b === 'cc-pvdz');

  return (
    <div className="bg-amber-50 rounded-xl border border-amber-200 overflow-hidden">
      {/* Header (always visible) */}
      <button
        type="button"
        onClick={() => setIsExpanded(!isExpanded)}
        className="w-full flex items-center justify-between px-5 py-3 text-left
          hover:bg-amber-100/50 transition-colors focus-visible:outline-none
          focus-visible:ring-2 focus-visible:ring-amber-500 focus-visible:ring-inset"
        aria-expanded={isExpanded}
      >
        <div className="flex items-center gap-2">
          <span className="text-amber-600 text-lg" aria-hidden="true">
            ?
          </span>
          <span className="font-semibold text-amber-900 text-sm">
            Discussion Questions
          </span>
        </div>
        <svg
          className={`w-4 h-4 text-amber-600 transition-transform ${isExpanded ? 'rotate-180' : ''}`}
          fill="none"
          viewBox="0 0 24 24"
          strokeWidth={2}
          stroke="currentColor"
          aria-hidden="true"
        >
          <path strokeLinecap="round" strokeLinejoin="round" d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {/* Expandable content */}
      {isExpanded && (
        <div className="px-5 pb-5 space-y-4">
          {/* Question 1: Observational */}
          <div className="space-y-1">
            <p className="text-sm font-medium text-amber-900">
              1. Observation
            </p>
            <p className="text-sm text-amber-800">
              How does the shape of the {shellType} radial profile for {elementSymbol} change
              as you go from {basisListText}? Pay attention to both the peak height
              and how far the function extends from the nucleus.
            </p>
          </div>

          {/* Question 2: Structural (if split-valence is present) */}
          {hasSplitValence && hasMinimal && (
            <div className="space-y-1">
              <p className="text-sm font-medium text-amber-900">
                2. Split-Valence Structure
              </p>
              <p className="text-sm text-amber-800">
                Notice that STO-3G has a single function for the valence shell, while
                split-valence bases ({basisLabels.filter((l) => l !== 'STO-3G').join(', ')})
                show <em>two separate functions</em> (inner and outer, solid and dashed lines).
                Why does splitting the valence description into two parts give the SCF
                procedure more variational freedom to minimize the energy?
              </p>
            </div>
          )}

          {/* Question 2 alt: if no minimal basis for comparison */}
          {hasSplitValence && !hasMinimal && (
            <div className="space-y-1">
              <p className="text-sm font-medium text-amber-900">
                2. Split-Valence Structure
              </p>
              <p className="text-sm text-amber-800">
                Split-valence bases represent each valence orbital with two separate
                contractions: a compact inner part and a more diffuse outer part.
                How does this allow the SCF procedure to adjust the effective
                size of the orbital during self-consistent optimization?
              </p>
            </div>
          )}

          {/* Question 3: Critical thinking */}
          <div className="space-y-1">
            <p className="text-sm font-medium text-amber-900">
              3. Cost vs. Accuracy
            </p>
            <p className="text-sm text-amber-800">
              The variational principle guarantees that a larger basis set gives an energy
              less than or equal to a smaller one. Does this mean the largest available
              basis is always the best choice for a practical calculation?
              {hasPolarization && (
                <> Consider when polarization functions (like d-orbitals in 6-31G*) matter
                most and when they add cost without significant improvement.</>
              )}
            </p>
          </div>

          {/* Misconception callout */}
          <div className="bg-amber-100 rounded-lg px-4 py-3 border border-amber-300">
            <p className="text-xs font-semibold text-amber-900 uppercase tracking-wide mb-1">
              Common Misconception
            </p>
            <p className="text-sm text-amber-800 italic mb-2">
              "More basis functions always means more accurate results."
            </p>
            <p className="text-xs font-semibold text-amber-900 uppercase tracking-wide mb-1">
              Clarification
            </p>
            <p className="text-sm text-amber-800">
              More basis functions do lower the energy (variational principle), but the
              improvement shows <strong>diminishing returns</strong>. Meanwhile, the
              computational cost of two-electron integrals scales as O(N<sup>4</sup>) with
              the number of basis functions N. Doubling the basis can increase computation
              time by up to 16x. The best basis set balances accuracy against cost for the
              specific problem at hand. Eventually, further improvement requires
              <em> correlation methods</em> (MP2, CCSD, ...), not more basis functions.
            </p>
          </div>
        </div>
      )}
    </div>
  );
});

// ============================================================================
// Utilities
// ============================================================================

/**
 * Format a list of basis set labels as a natural-language list.
 *
 * @example
 * formatBasisList(['STO-3G', '3-21G', '6-31G'])
 * // => "STO-3G, 3-21G, and 6-31G"
 */
function formatBasisList(labels: string[]): string {
  if (labels.length === 0) return '';
  if (labels.length === 1) return labels[0];
  if (labels.length === 2) return `${labels[0]} and ${labels[1]}`;
  return `${labels.slice(0, -1).join(', ')}, and ${labels[labels.length - 1]}`;
}
