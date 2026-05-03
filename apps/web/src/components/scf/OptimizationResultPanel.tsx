/**
 * OptimizationResultPanel - Display optimization results and export controls.
 *
 * Shows convergence status, final energy, gradient, step count, and
 * provides buttons to apply the optimized geometry and export as XYZ.
 *
 * @module components/scf/OptimizationResultPanel
 * @see US-075 Optimization UI + Trajectory Animation (AC4, AC5, AC6)
 */

import { useCallback } from 'react';
import type { OptimizationResult } from '../../worker/protocol';

/**
 * Bohr to Angstrom conversion factor.
 */
const BOHR_TO_ANGSTROM = 0.529177249;

/**
 * Map of atomic number to element symbol.
 */
const Z_TO_SYMBOL: Record<number, string> = {
  1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B',
  6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
};

/**
 * Props for OptimizationResultPanel.
 */
export interface OptimizationResultPanelProps {
  /** Optimization result */
  result: OptimizationResult;
  /** Initial atoms as [Z, x, y, z][] for the XYZ comment line */
  initialAtoms: [number, number, number, number][];
  /** Method used (for XYZ comment line) */
  method: string;
  /** Basis set (for XYZ comment line) */
  basisName: string;
  /** Callback to apply optimized geometry */
  onApplyGeometry: () => void;
  /** Whether geometry has already been applied */
  applied: boolean;
}

/**
 * Generate XYZ format text from optimized geometry.
 *
 * Converts coordinates from bohr to Angstrom (XYZ convention).
 */
function generateXyz(
  atoms: [number, number, number, number][],
  finalGeometry: [number, number, number][],
  energy: number,
  method: string,
  basisName: string
): string {
  const nAtoms = atoms.length;
  const lines: string[] = [];

  // Line 1: atom count
  lines.push(String(nAtoms));

  // Line 2: comment with method and energy
  lines.push(
    `Optimized (${method.toUpperCase()}/${basisName}, E=${energy.toFixed(6)} Ha) - IQCP`
  );

  // Lines 3+: atom coordinates in Angstrom
  for (let i = 0; i < nAtoms; i++) {
    const symbol = Z_TO_SYMBOL[atoms[i][0]] ?? 'X';
    const [x, y, z] = finalGeometry[i];
    const xA = x * BOHR_TO_ANGSTROM;
    const yA = y * BOHR_TO_ANGSTROM;
    const zA = z * BOHR_TO_ANGSTROM;
    lines.push(
      `${symbol.padEnd(4)}${xA.toFixed(6).padStart(12)}${yA.toFixed(6).padStart(12)}${zA.toFixed(6).padStart(12)}`
    );
  }

  return lines.join('\n') + '\n';
}

/**
 * Download a text string as a file.
 */
function downloadText(content: string, filename: string): void {
  const blob = new Blob([content], { type: 'text/plain' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

/**
 * Panel displaying geometry optimization results.
 *
 * Shown after optimization completes. Includes:
 * - Convergence status badge
 * - Final energy and gradient
 * - Step count and computation time
 * - "Apply Optimized Geometry" button
 * - "Export XYZ" button
 */
export function OptimizationResultPanel({
  result,
  initialAtoms,
  method,
  basisName,
  onApplyGeometry,
  applied,
}: OptimizationResultPanelProps) {
  const handleExportXyz = useCallback(() => {
    const xyz = generateXyz(
      initialAtoms,
      result.finalGeometry,
      result.finalEnergy,
      method,
      basisName
    );

    // Build molecule name from atom symbols
    const symbolCounts: Record<string, number> = {};
    for (const a of initialAtoms) {
      const sym = Z_TO_SYMBOL[a[0]] ?? 'X';
      symbolCounts[sym] = (symbolCounts[sym] ?? 0) + 1;
    }
    const moleculeName = Object.entries(symbolCounts)
      .map(([sym, count]) => (count > 1 ? `${sym}${count}` : sym))
      .join('');

    downloadText(xyz, `${moleculeName}_optimized.xyz`);
  }, [result, initialAtoms, method, basisName]);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      <div className="flex items-center justify-between mb-3">
        <h3 className="text-sm font-semibold text-slate-700">
          Optimization Result
        </h3>
        {/* Convergence badge */}
        {result.converged ? (
          <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs font-medium bg-green-100 text-green-800">
            <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20">
              <path
                fillRule="evenodd"
                d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z"
                clipRule="evenodd"
              />
            </svg>
            Converged
          </span>
        ) : (
          <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs font-medium bg-amber-100 text-amber-800">
            <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20">
              <path
                fillRule="evenodd"
                d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z"
                clipRule="evenodd"
              />
            </svg>
            Not Converged
          </span>
        )}
      </div>

      {/* Result details */}
      <div className="grid grid-cols-2 gap-3 text-sm mb-4">
        <div>
          <span className="text-slate-500">Final Energy</span>
          <p className="font-mono font-medium text-slate-800">
            {result.finalEnergy.toFixed(8)} Ha
          </p>
        </div>
        <div>
          <span className="text-slate-500">Max Gradient</span>
          <p className="font-mono font-medium text-slate-800">
            {result.steps.length > 0
              ? result.steps[result.steps.length - 1].maxGradient.toExponential(2)
              : 'N/A'}{' '}
            Ha/bohr
          </p>
        </div>
        <div>
          <span className="text-slate-500">Steps</span>
          <p className="font-medium text-slate-800">{result.totalSteps}</p>
        </div>
        <div>
          <span className="text-slate-500">Time</span>
          <p className="font-medium text-slate-800">
            {result.computeTimeMs < 1000
              ? `${result.computeTimeMs.toFixed(0)} ms`
              : `${(result.computeTimeMs / 1000).toFixed(1)} s`}
          </p>
        </div>
      </div>

      {/* Action buttons */}
      <div className="flex gap-2">
        <button
          type="button"
          onClick={onApplyGeometry}
          disabled={applied}
          className={`flex-1 px-3 py-2 rounded-lg text-sm font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-green-500 ${
            applied
              ? 'bg-green-100 text-green-700 cursor-default'
              : 'bg-green-600 text-white hover:bg-green-700'
          } disabled:opacity-70`}
          aria-label="Apply optimized geometry as active molecule"
        >
          {applied ? 'Geometry Applied' : 'Apply Optimized Geometry'}
        </button>
        <button
          type="button"
          onClick={handleExportXyz}
          className="px-3 py-2 rounded-lg text-sm font-medium transition-colors bg-slate-100 text-slate-700 hover:bg-slate-200 focus:outline-none focus-visible:ring-2 focus-visible:ring-slate-500"
          aria-label="Export optimized geometry as XYZ file"
        >
          Export XYZ
        </button>
      </div>
    </div>
  );
}
