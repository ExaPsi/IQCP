/**
 * BasisTypeToggle - Toggle between Cartesian and Spherical basis functions.
 *
 * Provides an educational toggle for selecting between Cartesian (6D, 10F)
 * and Spherical (5D, 7F) harmonic basis functions. Includes comprehensive
 * educational content explaining the differences and why spherical
 * harmonics are the standard in quantum chemistry.
 *
 * @module components/scf/custom/BasisTypeToggle
 */

import { useState } from 'react';
import { Math } from '../../common/Math';

/**
 * Basis set information for dynamic comparison display.
 */
export interface BasisInfo {
  /** Number of basis functions with current setting */
  nbf: number;
  /** Number of basis functions with Cartesian (for comparison) */
  nbfCartesian?: number;
  /** Number of basis functions with Spherical (for comparison) */
  nbfSpherical?: number;
  /** Whether the basis set contains d or higher angular momentum */
  hasHigherAngularMomentum?: boolean;
}

/**
 * Props for BasisTypeToggle component.
 */
export interface BasisTypeToggleProps {
  /** Whether spherical harmonics are enabled (true) or Cartesian (false) */
  useSpherical: boolean;
  /** Callback when the toggle changes */
  onChange: (useSpherical: boolean) => void;
  /** Disable the toggle */
  disabled?: boolean;
  /** Optional basis information for dynamic comparison display */
  basisInfo?: BasisInfo;
}

/**
 * Info icon SVG component.
 */
function InfoIcon({ className = '' }: { className?: string }) {
  return (
    <svg
      className={className}
      fill="none"
      viewBox="0 0 24 24"
      strokeWidth={1.5}
      stroke="currentColor"
      aria-hidden="true"
    >
      <path
        strokeLinecap="round"
        strokeLinejoin="round"
        d="M11.25 11.25l.041-.02a.75.75 0 011.063.852l-.708 2.836a.75.75 0 001.063.853l.041-.021M21 12a9 9 0 11-18 0 9 9 0 0118 0zm-9-3.75h.008v.008H12V8.25z"
      />
    </svg>
  );
}

/**
 * Chevron icon for collapsible sections.
 */
function ChevronIcon({ isOpen, className = '' }: { isOpen: boolean; className?: string }) {
  return (
    <svg
      className={`${className} transition-transform duration-200 ${isOpen ? 'rotate-180' : ''}`}
      fill="none"
      viewBox="0 0 24 24"
      strokeWidth={2}
      stroke="currentColor"
      aria-hidden="true"
    >
      <path strokeLinecap="round" strokeLinejoin="round" d="M19 9l-7 7-7-7" />
    </svg>
  );
}

/**
 * Educational content about basis function types.
 */
function BasisTypeEducation({ useSpherical }: { useSpherical: boolean }) {
  const [isExpanded, setIsExpanded] = useState(false);

  return (
    <div className="mt-3 bg-slate-50 rounded-lg border border-slate-200 overflow-hidden">
      {/* Collapsible header */}
      <button
        type="button"
        onClick={() => setIsExpanded(!isExpanded)}
        className="w-full px-3 py-2 flex items-center justify-between text-left hover:bg-slate-100 transition-colors"
        aria-expanded={isExpanded}
        aria-controls="basis-type-education-content"
      >
        <span className="flex items-center gap-2 text-xs font-medium text-slate-700">
          <InfoIcon className="w-4 h-4 text-blue-500" />
          Why does this matter?
        </span>
        <ChevronIcon isOpen={isExpanded} className="w-4 h-4 text-slate-500" />
      </button>

      {/* Expandable content */}
      {isExpanded && (
        <div id="basis-type-education-content" className="px-3 pb-3 text-xs text-slate-600 space-y-3">
          {/* Quick summary */}
          <div className="pt-1 border-t border-slate-200">
            <p>
              <strong>Key insight:</strong> Both basis types span the same physical space and give
              the <em>same energy</em>. The difference is in efficiency and avoiding redundancy.
            </p>
          </div>

          {/* Cartesian vs Spherical comparison table */}
          <div>
            <h4 className="font-medium text-slate-700 mb-1">Basis Function Count by Angular Momentum</h4>
            <table className="w-full text-xs border-collapse">
              <thead>
                <tr className="bg-slate-100">
                  <th className="px-2 py-1 text-left font-medium border border-slate-200">Shell</th>
                  <th className="px-2 py-1 text-center font-medium border border-slate-200">Cartesian</th>
                  <th className="px-2 py-1 text-center font-medium border border-slate-200">Spherical</th>
                  <th className="px-2 py-1 text-center font-medium border border-slate-200">Difference</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td className="px-2 py-1 border border-slate-200">s (l=0)</td>
                  <td className="px-2 py-1 text-center border border-slate-200">1</td>
                  <td className="px-2 py-1 text-center border border-slate-200">1</td>
                  <td className="px-2 py-1 text-center border border-slate-200 text-slate-400">0</td>
                </tr>
                <tr>
                  <td className="px-2 py-1 border border-slate-200">p (l=1)</td>
                  <td className="px-2 py-1 text-center border border-slate-200">3</td>
                  <td className="px-2 py-1 text-center border border-slate-200">3</td>
                  <td className="px-2 py-1 text-center border border-slate-200 text-slate-400">0</td>
                </tr>
                <tr className={!useSpherical ? 'bg-amber-50' : 'bg-green-50'}>
                  <td className="px-2 py-1 border border-slate-200 font-medium">d (l=2)</td>
                  <td className="px-2 py-1 text-center border border-slate-200">6</td>
                  <td className="px-2 py-1 text-center border border-slate-200">5</td>
                  <td className="px-2 py-1 text-center border border-slate-200 text-green-600 font-medium">-1</td>
                </tr>
                <tr className={!useSpherical ? 'bg-amber-50' : 'bg-green-50'}>
                  <td className="px-2 py-1 border border-slate-200 font-medium">f (l=3)</td>
                  <td className="px-2 py-1 text-center border border-slate-200">10</td>
                  <td className="px-2 py-1 text-center border border-slate-200">7</td>
                  <td className="px-2 py-1 text-center border border-slate-200 text-green-600 font-medium">-3</td>
                </tr>
                <tr>
                  <td className="px-2 py-1 border border-slate-200">g (l=4)</td>
                  <td className="px-2 py-1 text-center border border-slate-200">15</td>
                  <td className="px-2 py-1 text-center border border-slate-200">9</td>
                  <td className="px-2 py-1 text-center border border-slate-200 text-green-600 font-medium">-6</td>
                </tr>
              </tbody>
            </table>
            <p className="mt-1 text-slate-500">
              Formula: Cartesian = <Math>{String.raw`(l+1)(l+2)/2`}</Math>, Spherical = <Math>{String.raw`2l+1`}</Math>
            </p>
          </div>

          {/* Explanation of why the difference exists */}
          <div>
            <h4 className="font-medium text-slate-700 mb-1">Why the Difference?</h4>
            <p className="mb-2">
              <strong>Cartesian Gaussians</strong> use polynomial prefactors like <Math>{String.raw`x^i y^j z^k`}</Math> where <Math>{String.raw`i+j+k=l`}</Math>.
              For d-orbitals (l=2), this gives 6 functions:
            </p>
            <div className="bg-white rounded px-2 py-1 font-mono text-xs border border-slate-200 mb-2">
              <Math>{String.raw`d_{xx}, d_{xy}, d_{xz}, d_{yy}, d_{yz}, d_{zz}`}</Math>
            </div>
            <p className="mb-2">
              <strong>Spherical harmonics</strong> use angular functions <Math>{String.raw`Y_l^m(\theta, \phi)`}</Math> where <Math>{String.raw`m = -l, ..., +l`}</Math>.
              For d-orbitals, this gives 5 functions:
            </p>
            <div className="bg-white rounded px-2 py-1 font-mono text-xs border border-slate-200 mb-2">
              <Math>{String.raw`d_{-2}, d_{-1}, d_0, d_{+1}, d_{+2}`}</Math>
            </div>
            <p>
              The extra Cartesian function (<Math>{String.raw`d_{xx} + d_{yy} + d_{zz}`}</Math>) is proportional to <Math>{String.raw`r^2`}</Math>,
              which is an s-type function in disguise! This creates <strong>linear dependencies</strong> in larger basis sets.
            </p>
          </div>

          {/* Why spherical is preferred */}
          <div>
            <h4 className="font-medium text-slate-700 mb-1">Why Spherical is the Standard</h4>
            <ul className="list-disc list-inside space-y-1 ml-1">
              <li>Avoids linear dependencies that cause numerical instabilities</li>
              <li>Smaller matrices = faster computation</li>
              <li>Results in cleaner orbital analysis (pure angular momentum states)</li>
              <li>Used by default in PySCF, Gaussian, ORCA, and most QC programs</li>
            </ul>
          </div>

          {/* Important note */}
          <div className="bg-blue-50 rounded px-2 py-1.5 border border-blue-200">
            <p className="text-blue-800">
              <strong>Note:</strong> For s and p orbitals, there is no difference. The distinction only
              matters for basis sets with d, f, or higher angular momentum functions (like 6-31G*, cc-pVDZ, etc.).
            </p>
          </div>
        </div>
      )}
    </div>
  );
}

/**
 * BasisTypeToggle - Toggle between Cartesian and Spherical basis types.
 *
 * Provides a clear toggle with educational content explaining the difference
 * between Cartesian and Spherical harmonic basis functions.
 *
 * @example
 * ```tsx
 * const [useSpherical, setUseSpherical] = useState(false);
 * <BasisTypeToggle
 *   useSpherical={useSpherical}
 *   onChange={setUseSpherical}
 * />
 * ```
 */
export function BasisTypeToggle({
  useSpherical,
  onChange,
  disabled = false,
  basisInfo,
}: BasisTypeToggleProps) {
  // Determine if we can show meaningful comparison
  const showComparison = basisInfo?.nbfCartesian !== undefined &&
                         basisInfo?.nbfSpherical !== undefined &&
                         basisInfo.nbfCartesian !== basisInfo.nbfSpherical;

  return (
    <div className="space-y-2">
      {/* Label */}
      <label className="block text-sm font-medium text-slate-700">
        Basis Function Type
      </label>

      {/* Toggle buttons */}
      <div className="flex rounded-lg bg-slate-100 p-1">
        <button
          type="button"
          onClick={() => onChange(false)}
          disabled={disabled}
          className={`flex-1 px-3 py-1.5 rounded-md text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
            !useSpherical
              ? 'bg-white text-slate-800 shadow-sm'
              : 'text-slate-600 hover:text-slate-800'
          } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
          aria-pressed={!useSpherical}
        >
          Cartesian
          {basisInfo?.nbfCartesian !== undefined && (
            <span className="ml-1 text-xs opacity-70">({basisInfo.nbfCartesian})</span>
          )}
        </button>
        <button
          type="button"
          onClick={() => onChange(true)}
          disabled={disabled}
          className={`flex-1 px-3 py-1.5 rounded-md text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
            useSpherical
              ? 'bg-white text-slate-800 shadow-sm'
              : 'text-slate-600 hover:text-slate-800'
          } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
          aria-pressed={useSpherical}
        >
          Spherical
          {basisInfo?.nbfSpherical !== undefined && (
            <span className="ml-1 text-xs opacity-70">({basisInfo.nbfSpherical})</span>
          )}
        </button>
      </div>

      {/* Quick description */}
      <p className="text-xs text-slate-500">
        {useSpherical ? (
          <>
            Using <strong>spherical harmonics</strong> (5D, 7F) - the quantum chemistry standard
          </>
        ) : (
          <>
            Using <strong>Cartesian Gaussians</strong> (6D, 10F) - includes redundant components
          </>
        )}
      </p>

      {/* Dynamic comparison display when basis info is available */}
      {showComparison && (
        <div className="bg-green-50 rounded-lg px-3 py-2 border border-green-200 text-xs">
          <div className="flex items-center gap-2 text-green-800">
            <svg className="w-4 h-4 flex-shrink-0" fill="currentColor" viewBox="0 0 20 20">
              <path
                fillRule="evenodd"
                d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z"
                clipRule="evenodd"
              />
            </svg>
            <span>
              <strong>For this system:</strong>{' '}
              {basisInfo.nbfCartesian} Cartesian vs {basisInfo.nbfSpherical} Spherical basis functions
              {' '}
              ({basisInfo.nbfCartesian! - basisInfo.nbfSpherical!} fewer with spherical)
            </span>
          </div>
        </div>
      )}

      {/* Note when no difference */}
      {basisInfo?.hasHigherAngularMomentum === false && (
        <div className="bg-slate-50 rounded-lg px-3 py-2 border border-slate-200 text-xs">
          <div className="flex items-center gap-2 text-slate-600">
            <svg className="w-4 h-4 flex-shrink-0" fill="currentColor" viewBox="0 0 20 20">
              <path
                fillRule="evenodd"
                d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z"
                clipRule="evenodd"
              />
            </svg>
            <span>
              This basis set has only s and p orbitals - both types give identical results.
            </span>
          </div>
        </div>
      )}

      {/* Educational content */}
      <BasisTypeEducation useSpherical={useSpherical} />
    </div>
  );
}

export default BasisTypeToggle;
