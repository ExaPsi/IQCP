/**
 * ShellDataTable - Displays primitive Gaussian parameters for a selected shell.
 *
 * Shows a table of exponents (alpha_k), raw contraction coefficients (c_k),
 * computed normalization constants (N_k), and effective coefficients (d_k = c_k * N_k)
 * for all primitives in a contracted shell.
 *
 * Features:
 * - Sort by exponent (ascending/descending toggle)
 * - Copy-to-clipboard as TSV for student reports
 * - Educational note for d-type shells explaining component-dependent normalization
 *
 * @module components/basisExplorer/ShellDataTable
 */

import { useState, useMemo, useCallback, memo } from 'react';
import { primitiveNormalization } from './normalization';

// ============================================================================
// Types
// ============================================================================

/**
 * Props for the ShellDataTable component.
 */
export interface ShellDataTableProps {
  /** Exponents of the primitive Gaussians in this shell */
  exponents: number[];
  /** Raw contraction coefficients c_k */
  coefficients: number[];
  /** Angular momentum quantum number (0=s, 1=p, 2=d) */
  angularMomentum: number;
  /** Angular momentum label for display ("s", "p", "d") */
  angularMomentumLabel: string;
  /** Shell display label (e.g., "1s", "2p", "3d") */
  shellLabel: string;
}

/**
 * A single row of computed primitive data.
 */
interface PrimitiveRow {
  /** 1-based index */
  k: number;
  /** Primitive exponent */
  alpha: number;
  /** Raw contraction coefficient */
  c: number;
  /** Normalization constant */
  N: number;
  /** Effective coefficient d = c * N */
  d: number;
}

/** Sort direction for the exponent column. */
type SortDirection = 'asc' | 'desc';

// ============================================================================
// Number formatting
// ============================================================================

/**
 * Format a number with 6 significant figures in scientific notation.
 *
 * Uses `toPrecision(6)` which produces scientific notation for very large
 * or very small values and fixed notation for moderate values, then
 * normalizes to always use exponential form for consistency.
 */
function formatSci(value: number): string {
  return value.toExponential(6);
}

// ============================================================================
// Component
// ============================================================================

/**
 * Shell data table displaying primitive Gaussian parameters.
 *
 * Computes normalization constants (N_k) and effective coefficients (d_k)
 * from the raw shell data using `primitiveNormalization`. All computation
 * is pure TypeScript in the render path -- no WASM calls needed.
 *
 * @example
 * ```tsx
 * <ShellDataTable
 *   exponents={[3.425, 0.624, 0.169]}
 *   coefficients={[0.154, 0.535, 0.445]}
 *   angularMomentum={0}
 *   angularMomentumLabel="s"
 *   shellLabel="1s"
 * />
 * ```
 */
export const ShellDataTable = memo(function ShellDataTable({
  exponents,
  coefficients,
  angularMomentum,
  angularMomentumLabel,
  shellLabel,
}: ShellDataTableProps) {
  // Sort state: default descending (largest exponent first)
  const [sortDirection, setSortDirection] = useState<SortDirection>('desc');
  // Copy feedback state: tracks which format was last copied
  const [copiedFormat, setCopiedFormat] = useState<'tsv' | 'json' | 'latex' | null>(null);

  // Compute rows with normalization and effective coefficients
  const rows = useMemo((): PrimitiveRow[] => {
    const unsorted: PrimitiveRow[] = exponents.map((alpha, i) => {
      const c = coefficients[i];
      const N = primitiveNormalization(alpha, angularMomentum);
      return {
        k: i + 1,
        alpha,
        c,
        N,
        d: c * N,
      };
    });

    // Sort by exponent
    const sorted = [...unsorted].sort((a, b) =>
      sortDirection === 'desc' ? b.alpha - a.alpha : a.alpha - b.alpha
    );

    return sorted;
  }, [exponents, coefficients, angularMomentum, sortDirection]);

  // Toggle sort direction
  const toggleSort = useCallback(() => {
    setSortDirection((prev) => (prev === 'desc' ? 'asc' : 'desc'));
  }, []);

  // Helper to copy text and show feedback
  const copyToClipboard = useCallback(async (text: string, format: 'tsv' | 'json' | 'latex') => {
    try {
      await navigator.clipboard.writeText(text);
      setCopiedFormat(format);
      setTimeout(() => setCopiedFormat(null), 2000);
    } catch {
      console.warn('Clipboard API not available');
    }
  }, []);

  // Copy table data as TSV
  const handleCopyTsv = useCallback(async () => {
    const header = `# Shell: ${shellLabel} (l=${angularMomentum}) | Type: ${angularMomentumLabel}`;
    const columnHeader = 'k\talpha_k\tc_k\tN_k\td_k';
    const dataLines = rows.map(
      (row) =>
        `${row.k}\t${formatSci(row.alpha)}\t${formatSci(row.c)}\t${formatSci(row.N)}\t${formatSci(row.d)}`
    );
    const tsv = [header, columnHeader, ...dataLines].join('\n');
    await copyToClipboard(tsv, 'tsv');
  }, [rows, shellLabel, angularMomentum, angularMomentumLabel, copyToClipboard]);

  // Copy table data as JSON
  const handleCopyJson = useCallback(async () => {
    const jsonData = rows.map((row) => ({
      k: row.k,
      alpha: row.alpha,
      c: row.c,
      N: row.N,
      d: row.d,
    }));
    const json = JSON.stringify(jsonData, null, 2);
    await copyToClipboard(json, 'json');
  }, [rows, copyToClipboard]);

  // Copy table data as LaTeX tabular
  const handleCopyLatex = useCallback(async () => {
    const lines: string[] = [];
    lines.push(`% Shell: ${shellLabel} (l=${angularMomentum}, ${angularMomentumLabel}-type)`);
    lines.push('\\begin{tabular}{ccccc}');
    lines.push('\\hline');
    lines.push('$k$ & $\\alpha_k$ & $c_k$ & $N_k$ & $d_k$ \\\\');
    lines.push('\\hline');
    for (const row of rows) {
      lines.push(
        `${row.k} & ${formatSci(row.alpha)} & ${formatSci(row.c)} & ${formatSci(row.N)} & ${formatSci(row.d)} \\\\`
      );
    }
    lines.push('\\hline');
    lines.push('\\end{tabular}');
    const latex = lines.join('\n');
    await copyToClipboard(latex, 'latex');
  }, [rows, shellLabel, angularMomentum, angularMomentumLabel, copyToClipboard]);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
      {/* Header with copy button */}
      <div className="flex items-center justify-between mb-4">
        <h3 className="text-sm font-semibold text-slate-700">
          Primitive Data: {shellLabel} shell
        </h3>
        <div className="inline-flex rounded-md border border-slate-300" role="group" aria-label="Copy table data">
          {([
            { format: 'tsv' as const, label: 'TSV', handler: handleCopyTsv },
            { format: 'json' as const, label: 'JSON', handler: handleCopyJson },
            { format: 'latex' as const, label: 'LaTeX', handler: handleCopyLatex },
          ]).map(({ format, label, handler }, idx) => (
            <button
              key={format}
              type="button"
              onClick={handler}
              aria-label={copiedFormat === format ? `${label} copied` : `Copy as ${label}`}
              className={`
                inline-flex items-center gap-1 px-2.5 py-1.5 text-xs font-medium
                transition-colors duration-150
                focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
                ${idx === 0 ? 'rounded-l-md' : ''}
                ${idx === 2 ? 'rounded-r-md' : ''}
                ${idx > 0 ? 'border-l border-slate-300' : ''}
                ${
                  copiedFormat === format
                    ? 'bg-green-50 text-green-700'
                    : 'bg-white text-slate-600 hover:bg-slate-50'
                }
              `}
            >
              {copiedFormat === format ? (
                <>
                  <svg className="w-3 h-3" fill="none" viewBox="0 0 24 24" strokeWidth={2} stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" d="M4.5 12.75l6 6 9-13.5" />
                  </svg>
                  {label}
                </>
              ) : (
                label
              )}
            </button>
          ))}
        </div>
      </div>

      {/* Data table */}
      <div className="overflow-x-auto">
        <table className="w-full text-sm font-mono" role="table">
          <thead>
            <tr className="border-b-2 border-slate-200">
              <th
                scope="col"
                className="px-3 py-2 text-left text-xs font-semibold text-slate-500 tracking-wider"
              >
                <span className="italic">k</span>
              </th>
              <th
                scope="col"
                role="columnheader"
                aria-sort={sortDirection === 'desc' ? 'descending' : 'ascending'}
                aria-label={`alpha_k, sort ${sortDirection === 'desc' ? 'descending' : 'ascending'}`}
                className="px-3 py-2 text-right text-xs font-semibold text-slate-500 tracking-wider cursor-pointer select-none hover:text-slate-700 transition-colors"
                onClick={toggleSort}
              >
                <span className="inline-flex items-center gap-1">
                  <span>&alpha;<sub>k</sub></span>
                  <span className="text-slate-400" aria-hidden="true">
                    {sortDirection === 'desc' ? '\u2193' : '\u2191'}
                  </span>
                </span>
              </th>
              <th
                scope="col"
                className="px-3 py-2 text-right text-xs font-semibold text-slate-500 tracking-wider"
              >
                <span className="italic">c</span><sub>k</sub>
              </th>
              <th
                scope="col"
                className="px-3 py-2 text-right text-xs font-semibold text-slate-500 tracking-wider"
              >
                N<sub>k</sub>
              </th>
              <th
                scope="col"
                className="px-3 py-2 text-right text-xs font-semibold text-slate-500 tracking-wider"
              >
                <span className="italic">d</span><sub>k</sub>
              </th>
            </tr>
          </thead>
          <tbody className="divide-y divide-slate-100">
            {rows.map((row, idx) => (
              <tr
                key={row.k}
                className={idx % 2 === 0 ? 'bg-white' : 'bg-slate-50/50'}
              >
                <td className="px-3 py-2 text-slate-500 text-left">
                  {row.k}
                </td>
                <td className="px-3 py-2 text-slate-800 text-right tabular-nums">
                  {formatSci(row.alpha)}
                </td>
                <td className="px-3 py-2 text-slate-800 text-right tabular-nums">
                  {formatSci(row.c)}
                </td>
                <td className="px-3 py-2 text-slate-800 text-right tabular-nums">
                  {formatSci(row.N)}
                </td>
                <td className="px-3 py-2 text-slate-800 text-right tabular-nums">
                  {formatSci(row.d)}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* d-type normalization note */}
      {angularMomentum === 2 && (
        <div className="mt-4 bg-blue-50 rounded-lg px-4 py-3 border border-blue-100">
          <p className="text-xs text-blue-800">
            <strong>Note:</strong> N<sub>k</sub> shown is for the d<sub>xy</sub> component
            (Cartesian powers 1,1,0). For d<sub>xx</sub>, d<sub>yy</sub>, d<sub>zz</sub> components,
            multiply N<sub>k</sub> by 1/sqrt(3) due to the (2n-1)!! double-factorial terms.
            The d<sub>xy</sub> normalization is: N = (2&alpha;/&pi;)<sup>3/4</sup> &times; 4&alpha;.
          </p>
        </div>
      )}

      {/* Educational footer */}
      <div className="mt-4 pt-3 border-t border-slate-100">
        <p className="text-xs text-slate-500">
          <strong>&alpha;<sub>k</sub></strong> = exponent,{' '}
          <strong>c<sub>k</sub></strong> = raw contraction coefficient,{' '}
          <strong>N<sub>k</sub></strong> = normalization constant,{' '}
          <strong>d<sub>k</sub></strong> = c<sub>k</sub> &times; N<sub>k</sub> (effective coefficient
          used in R(r) = r<sup>l</sup> &sum; d<sub>k</sub> exp(-&alpha;<sub>k</sub>r<sup>2</sup>)).
        </p>
      </div>
    </div>
  );
});

export default ShellDataTable;
