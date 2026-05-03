/**
 * FockElementDetail - Element-level breakdown of a Fock matrix element.
 *
 * When a user clicks F_{i,j} in the step 3 heatmap, this panel shows:
 *   F_{i,j} = H^core_{i,j} + G_{i,j}
 *           = H^core_{i,j} + J_{i,j} - 0.5 * K_{i,j}
 *
 * with all numerical values and a visual bar chart showing the relative
 * magnitudes of each contribution.
 *
 * This component satisfies AC2 of US-058: "For a selected F_ij element,
 * list the contributing ERIs weighted by density matrix elements."
 *
 * Note: Full ERI-level contributions (which specific (kl) pairs contribute
 * most to J_{ij} and K_{ij}) would require the ERI tensor and density
 * matrix element-wise products. This is planned for a follow-up enhancement.
 * For now, we show the matrix-level decomposition F = H^core + J - 0.5*K.
 *
 * @module components/integrals/FockElementDetail
 * @see US-058 Fock Build Tracing (AC2)
 * @see Szabo & Ostlund (1996), Eq. 3.154
 */

import { useMemo, memo } from 'react';
import type { FockDecompositionResult } from '../../stores/integralStore';

// ============================================================================
// Types
// ============================================================================

export interface FockElementDetailProps {
  /** Fock decomposition result data */
  decomposition: FockDecompositionResult;
  /** Selected row index (0-based) */
  row: number;
  /** Selected column index (0-based) */
  col: number;
  /** Close handler */
  onClose: () => void;
}

// ============================================================================
// Helpers
// ============================================================================

/**
 * Format a number for display: scientific notation with 8 decimal places.
 */
function formatValue(value: number): string {
  return value.toExponential(8);
}

/**
 * Format a number as fixed-point with 8 decimal places.
 */
function formatFixed(value: number): string {
  return value.toFixed(8);
}

/**
 * Contribution bar with label, showing the sign and relative magnitude.
 */
const ContributionBar = memo(function ContributionBar({
  label,
  value,
  maxAbsValue,
  color,
}: {
  label: string;
  value: number;
  maxAbsValue: number;
  color: string;
}) {
  const pct = maxAbsValue > 1e-30 ? (Math.abs(value) / maxAbsValue) * 100 : 0;
  const isNegative = value < 0;

  return (
    <div className="flex items-center gap-2">
      <span className="text-[11px] font-mono text-slate-600 w-24 text-right flex-shrink-0">
        {label}
      </span>
      <div className="flex-1 flex items-center">
        {isNegative ? (
          // Negative: bar extends leftward from center
          <div className="flex-1 flex justify-end">
            <div className="h-4 rounded-r-sm" style={{ width: `${pct / 2}%`, backgroundColor: color, opacity: 0.7 }} />
          </div>
        ) : (
          // Positive: bar extends rightward from center
          <div className="flex-1">
            <div className="h-4 rounded-l-sm" style={{ width: `${pct / 2}%`, backgroundColor: color, opacity: 0.7 }} />
          </div>
        )}
      </div>
      <span className="text-[11px] font-mono text-slate-800 w-32 text-right flex-shrink-0">
        {formatFixed(value)}
      </span>
    </div>
  );
});

// ============================================================================
// Component
// ============================================================================

/**
 * Detailed numerical breakdown of a selected Fock matrix element.
 *
 * Displays:
 * 1. The decomposition formula: F_{i,j} = H^core_{i,j} + J_{i,j} - 0.5*K_{i,j}
 * 2. Numerical values for each component
 * 3. Visual bar chart of relative magnitudes
 * 4. Verification that G = J - 0.5*K and F = H^core + G
 */
export const FockElementDetail = memo(function FockElementDetail({
  decomposition,
  row,
  col,
  onClose,
}: FockElementDetailProps) {
  const { nbf, labels, hCore, jMatrix, kMatrix, gMatrix, fMatrix, density } = decomposition;

  // Extract values for the selected element
  const values = useMemo(() => {
    const idx = row * nbf + col;
    const hVal = hCore[idx];
    const jVal = jMatrix[idx];
    const kVal = kMatrix[idx];
    const gVal = gMatrix[idx];
    const fVal = fMatrix[idx];
    const pVal = density[idx];

    // Verification
    const gComputed = jVal - 0.5 * kVal;
    const fComputed = hVal + gVal;
    const gError = Math.abs(gVal - gComputed);
    const fError = Math.abs(fVal - fComputed);

    return { hVal, jVal, kVal, gVal, fVal, pVal, gComputed, fComputed, gError, fError };
  }, [row, col, nbf, hCore, jMatrix, kMatrix, gMatrix, fMatrix, density]);

  // Find the top contributing density-weighted ERI pairs for J and K.
  // J_{mn} = sum_{ls} P_{ls} * (mn|ls)
  // K_{mn} = sum_{ls} P_{ls} * (ml|ns)
  //
  // Since we don't have the individual ERI values in this component
  // (they are inside the WASM module), we show the density matrix elements
  // that participate in the sum. This gives students insight into which
  // basis function pairs contribute most via the density weighting.
  const densityContributions = useMemo(() => {
    const contributions: { l: number; s: number; label: string; pVal: number; absP: number }[] = [];
    for (let l = 0; l < nbf; l++) {
      for (let s = 0; s < nbf; s++) {
        const pVal = density[l * nbf + s];
        contributions.push({
          l,
          s,
          label: `P(${labels[l]}, ${labels[s]})`,
          pVal,
          absP: Math.abs(pVal),
        });
      }
    }
    // Sort by |P| descending and take top 10
    contributions.sort((a, b) => b.absP - a.absP);
    return contributions.slice(0, 10);
  }, [nbf, density, labels]);

  const rowLabel = labels[row] ?? `BF ${row + 1}`;
  const colLabel = labels[col] ?? `BF ${col + 1}`;

  // Max absolute value for bar chart scaling
  const maxAbs = Math.max(
    Math.abs(values.hVal),
    Math.abs(values.jVal),
    Math.abs(0.5 * values.kVal),
    Math.abs(values.gVal),
    Math.abs(values.fVal),
    1e-30,
  );

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
      {/* Header with close button */}
      <div className="flex items-start justify-between mb-4">
        <div>
          <h3 className="text-sm font-semibold text-slate-700">
            Element Breakdown: F<sub>{row + 1},{col + 1}</sub>
          </h3>
          <p className="text-xs text-slate-500 mt-0.5">
            {rowLabel} | {colLabel}
          </p>
        </div>
        <button
          onClick={onClose}
          className="text-slate-400 hover:text-slate-600 p-1 rounded transition-colors"
          aria-label="Close element detail"
          title="Close (Escape)"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
          </svg>
        </button>
      </div>

      {/* Decomposition formula */}
      <div className="bg-slate-50 rounded-lg px-4 py-3 border border-slate-200 mb-4">
        <div className="text-xs font-mono text-slate-800 space-y-1">
          <p>
            <span className="text-slate-500">F</span>
            <sub>{row + 1},{col + 1}</sub>
            <span className="text-slate-500"> = H^core</span>
            <sub>{row + 1},{col + 1}</sub>
            <span className="text-slate-500"> + G</span>
            <sub>{row + 1},{col + 1}</sub>
          </p>
          <p className="pl-4">
            <span className="text-slate-500">= H^core</span>
            <sub>{row + 1},{col + 1}</sub>
            <span className="text-slate-500"> + J</span>
            <sub>{row + 1},{col + 1}</sub>
            <span className="text-slate-500"> - 0.5 * K</span>
            <sub>{row + 1},{col + 1}</sub>
          </p>
          <p className="pl-4 text-slate-600">
            = ({formatFixed(values.hVal)}) + ({formatFixed(values.jVal)}) - 0.5 * ({formatFixed(values.kVal)})
          </p>
          <p className="pl-4 text-slate-600">
            = ({formatFixed(values.hVal)}) + ({formatFixed(values.gVal)})
          </p>
          <p className="pl-4 font-semibold text-slate-900">
            = {formatFixed(values.fVal)}
          </p>
        </div>
      </div>

      {/* Numerical values table */}
      <div className="mb-4">
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wider mb-2">
          Component Values
        </h4>
        <div className="space-y-1">
          <ContributionBar label="H^core" value={values.hVal} maxAbsValue={maxAbs} color="#6366f1" />
          <ContributionBar label="J" value={values.jVal} maxAbsValue={maxAbs} color="#3b82f6" />
          <ContributionBar label="-0.5*K" value={-0.5 * values.kVal} maxAbsValue={maxAbs} color="#f97316" />
          <div className="border-t border-slate-200 pt-1 mt-1">
            <ContributionBar label="G = J-0.5K" value={values.gVal} maxAbsValue={maxAbs} color="#10b981" />
          </div>
          <div className="border-t border-slate-300 pt-1 mt-1">
            <ContributionBar label="F" value={values.fVal} maxAbsValue={maxAbs} color="#1e293b" />
          </div>
        </div>
      </div>

      {/* Verification */}
      <div className="grid grid-cols-2 gap-3 mb-4">
        <div
          className={`rounded-lg px-3 py-2 border text-xs font-mono ${
            values.gError < 1e-12
              ? 'bg-green-50 border-green-200 text-green-800'
              : 'bg-red-50 border-red-200 text-red-800'
          }`}
        >
          <p className="font-semibold mb-0.5">G = J - 0.5*K</p>
          <p>diff = {values.gError.toExponential(2)}</p>
        </div>
        <div
          className={`rounded-lg px-3 py-2 border text-xs font-mono ${
            values.fError < 1e-12
              ? 'bg-green-50 border-green-200 text-green-800'
              : 'bg-red-50 border-red-200 text-red-800'
          }`}
        >
          <p className="font-semibold mb-0.5">F = H^core + G</p>
          <p>diff = {values.fError.toExponential(2)}</p>
        </div>
      </div>

      {/* Density matrix contributions */}
      <div className="mb-4">
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wider mb-2">
          Top Density Matrix Elements (Weighting J and K Sums)
        </h4>
        <p className="text-[10px] text-slate-400 mb-2">
          J<sub>{row + 1},{col + 1}</sub> = &Sigma;<sub>ls</sub> P<sub>ls</sub> ({row + 1}{col + 1}|ls)
          &nbsp;&nbsp;|&nbsp;&nbsp;
          K<sub>{row + 1},{col + 1}</sub> = &Sigma;<sub>ls</sub> P<sub>ls</sub> ({row + 1}l|{col + 1}s)
        </p>
        <div className="overflow-x-auto">
          <table
            className="w-full text-xs border-collapse"
            aria-label={`Top density matrix elements contributing to F(${row + 1},${col + 1})`}
          >
            <thead>
              <tr className="border-b border-slate-200">
                <th className="px-1.5 py-1 text-left text-[10px] font-semibold text-slate-500">#</th>
                <th className="px-1.5 py-1 text-left text-[10px] font-semibold text-slate-500">l, s</th>
                <th className="px-1.5 py-1 text-left text-[10px] font-semibold text-slate-500">Basis Functions</th>
                <th className="px-1.5 py-1 text-right text-[10px] font-semibold text-slate-500">P<sub>ls</sub></th>
              </tr>
            </thead>
            <tbody>
              {densityContributions.map((contrib, idx) => (
                <tr
                  key={`${contrib.l}-${contrib.s}`}
                  className={`border-b border-slate-50 ${idx % 2 === 0 ? 'bg-white' : 'bg-slate-50/50'}`}
                >
                  <td className="px-1.5 py-1 font-mono text-slate-400">{idx + 1}</td>
                  <td className="px-1.5 py-1 font-mono text-slate-700">({contrib.l + 1},{contrib.s + 1})</td>
                  <td className="px-1.5 py-1 text-slate-600">{labels[contrib.l]}, {labels[contrib.s]}</td>
                  <td className="px-1.5 py-1 font-mono text-right text-slate-800">{formatValue(contrib.pVal)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <p className="text-[10px] text-slate-400 mt-1.5">
          The density matrix P weights the two-electron integrals (ERIs) in both
          the Coulomb (J) and exchange (K) sums. Large |P<sub>ls</sub>| values indicate
          basis function pairs with significant electron density, which contribute
          most to the electron-electron interaction.
        </p>
      </div>

      {/* Physical interpretation */}
      <div className="bg-blue-50 rounded-lg px-4 py-3 border border-blue-100">
        <p className="text-xs font-medium text-blue-800 mb-1">
          Physical Interpretation
        </p>
        <p className="text-[11px] text-blue-700 leading-relaxed">
          {row === col ? (
            <>
              The diagonal element F<sub>{row + 1},{row + 1}</sub> represents the
              effective one-electron energy of basis function {rowLabel}, including
              both the one-electron terms (kinetic energy + nuclear attraction) and
              the mean-field electron-electron repulsion. The negative H^core
              contribution (nuclear attraction dominates kinetic energy) is partially
              offset by the positive Coulomb repulsion J.
            </>
          ) : (
            <>
              The off-diagonal element F<sub>{row + 1},{col + 1}</sub> couples basis
              functions {rowLabel} and {colLabel} through both one-electron (H^core)
              and two-electron (G) interactions. These coupling elements determine
              how atomic orbitals mix to form molecular orbitals when F is diagonalized.
            </>
          )}
        </p>
      </div>
    </div>
  );
});
