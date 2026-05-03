/**
 * RysInternalsPanel - Display computational internals for Rys quadrature.
 *
 * Shows detailed information about the quadrature computation including
 * algorithm steps (Schmidt orthogonalization + companion matrix QR),
 * sanity checks (roots in [0,1), sorted, weights non-negative),
 * weight function formula, three-term recurrence, moment reconstruction verification,
 * Hankel matrix structure (theoretical derivation), and numerical stability.
 *
 * @module components/rys/RysInternalsPanel
 */

import { useMemo, useState } from 'react';
import type { RysComputeResult, RysMethod, BoysEvalResult, BoysEvalRequest } from '../../worker/protocol';
import { Math as MathFormula, MathDisplay } from '../common/Math';
import { useWorker } from '../../hooks/useWorker';

/**
 * Method details for each computational approach.
 * Steps use LaTeX notation for rendering with KaTeX.
 */
const METHOD_DETAILS: Record<
  RysMethod,
  { label: string; description: string; steps: { text: string; latex?: string }[] }
> = {
  special: {
    label: 'Special Case',
    description:
      'Used for T=0 or n=1 where analytical formulas are available.',
    steps: [{ text: 'Direct analytical formula applied' }],
  },
  standard: {
    label: 'Standard Algorithm (RDK)',
    description:
      'Rys-Dupuis-King algorithm via Schmidt orthogonalization: builds orthonormal polynomials from Boys function moments, finds roots via companion matrix QR, and computes weights from polynomial evaluations at the roots.',
    steps: [
      { text: '1. Compute moments:', latex: '\\mu_k = F_k(T) \\text{ for } k = 0, \\ldots, 2n' },
      { text: '2. Schmidt orthogonalization:', latex: 'P_0, P_1, \\ldots, P_n \\text{ from } \\{\\mu_k\\}' },
      { text: '3. Companion matrix QR:', latex: '\\text{roots of } P_n' },
      { text: '4. Weights:', latex: 'w_k = 1 / \\sum_j P_j(r_k)^2' },
    ],
  },
};

/**
 * Approximate condition numbers for Hankel matrix at T~1.
 * Based on empirical observations from Appendix E of the lecture notes.
 */
const HANKEL_CONDITION_NUMBERS: Record<number, string> = {
  2: '10^2',
  3: '10^4',
  4: '10^6',
  5: '10^9',
  6: '10^{12}',
  7: '10^{14}',
  8: '10^{16}',
  9: '> 10^{16}',
  10: '> 10^{16}',
};

/**
 * Props for RysInternalsPanel.
 */
interface RysInternalsPanelProps {
  /** The computation result to display internals for */
  result: RysComputeResult;
}

/**
 * Check if all roots are in the open interval (0, 1).
 */
function checkRootsInRange(roots: number[]): { valid: boolean; count: number; total: number } {
  let validCount = 0;
  for (const root of roots) {
    if (root > 0 && root < 1) {
      validCount++;
    }
  }
  return {
    valid: validCount === roots.length,
    count: validCount,
    total: roots.length,
  };
}

/**
 * Check if roots are sorted in ascending order.
 */
function checkRootsSorted(roots: number[]): boolean {
  for (let i = 1; i < roots.length; i++) {
    if (roots[i] < roots[i - 1]) {
      return false;
    }
  }
  return true;
}

/**
 * Check if all weights are strictly positive.
 */
function checkWeightsPositive(weights: number[]): { valid: boolean; count: number; total: number } {
  let validCount = 0;
  for (const weight of weights) {
    if (weight > 0) {
      validCount++;
    }
  }
  return {
    valid: validCount === weights.length,
    count: validCount,
    total: weights.length,
  };
}

/**
 * Status indicator icon component.
 */
function StatusIcon({ valid }: { valid: boolean }) {
  if (valid) {
    return (
      <span className="text-green-600" aria-label="Check passed">
        <svg className="w-4 h-4 inline" fill="currentColor" viewBox="0 0 20 20">
          <path
            fillRule="evenodd"
            d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z"
            clipRule="evenodd"
          />
        </svg>
      </span>
    );
  }
  return (
    <span className="text-red-600" aria-label="Check failed">
      <svg className="w-4 h-4 inline" fill="currentColor" viewBox="0 0 20 20">
        <path
          fillRule="evenodd"
          d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z"
          clipRule="evenodd"
        />
      </svg>
    </span>
  );
}

/**
 * Collapsible section component.
 */
function CollapsibleSection({
  title,
  defaultOpen = false,
  className = '',
  children,
}: {
  title: string;
  defaultOpen?: boolean;
  className?: string;
  children: React.ReactNode;
}) {
  const [isOpen, setIsOpen] = useState(defaultOpen);

  return (
    <div className={`border border-slate-200 rounded-lg overflow-hidden ${className}`}>
      <button
        type="button"
        onClick={() => setIsOpen(!isOpen)}
        className="w-full px-3 py-2 flex items-center justify-between text-left bg-slate-100 hover:bg-slate-200 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-inset focus-visible:ring-blue-500"
        aria-expanded={isOpen}
      >
        <span className="text-xs font-semibold text-slate-600 uppercase tracking-wide">
          {title}
        </span>
        <svg
          className={`w-4 h-4 text-slate-500 transition-transform ${isOpen ? 'rotate-180' : ''}`}
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>
      {isOpen && <div className="p-3">{children}</div>}
    </div>
  );
}

/**
 * Component to compute and display moment reconstruction with actual Boys values.
 */
function MomentReconstructionTable({
  roots,
  weights,
  t,
  nroots,
}: {
  roots: number[];
  weights: number[];
  t: number;
  nroots: number;
}) {
  const { send, isReady } = useWorker();
  const [boysValues, setBoysValues] = useState<BoysEvalResult[] | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Compute the required number of moments: 0 to 2n-1
  const maxMomentOrder = 2 * nroots - 1;

  // Fetch Boys values when component mounts or parameters change
  // Use useMemo for triggering the async effect (to avoid stale closures)
  useMemo(() => {
    if (!isReady) return;

    const fetchBoysValues = async () => {
      setIsLoading(true);
      setError(null);
      try {
        // Make individual boys_eval calls for m = 0, 1, ..., maxMomentOrder
        // Use Promise.all for parallel execution
        const promises: Promise<BoysEvalResult>[] = [];
        for (let order = 0; order <= maxMomentOrder; order++) {
          promises.push(
            send<BoysEvalResult>({
              type: 'boys_eval',
              m: order,
              T: t,
            } as Omit<BoysEvalRequest, 'requestId'>)
          );
        }
        const results = await Promise.all(promises);
        setBoysValues(results);
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Failed to compute Boys values');
      } finally {
        setIsLoading(false);
      }
    };

    fetchBoysValues();
  }, [isReady, send, t, maxMomentOrder]);

  // Compute reconstructed moments from roots and weights
  const reconstructedMoments = useMemo(() => {
    const moments: number[] = [];
    const maxK = Math.min(maxMomentOrder, 5); // Show up to 6 moments (m=0 to m=5)

    for (let k = 0; k <= maxK; k++) {
      let reconstructed = 0;
      for (let i = 0; i < roots.length; i++) {
        reconstructed += weights[i] * Math.pow(roots[i], k);
      }
      moments.push(reconstructed);
    }
    return moments;
  }, [roots, weights, maxMomentOrder]);

  // Compute true moments from Boys values: mu_k = 2 * F_k(T)
  const trueMoments = useMemo(() => {
    if (!boysValues) return null;
    return boysValues.slice(0, reconstructedMoments.length).map((bv) => 2 * bv.value);
  }, [boysValues, reconstructedMoments.length]);

  // Compute errors
  const errors = useMemo(() => {
    if (!trueMoments) return null;
    return reconstructedMoments.map((recon, k) => {
      const trueVal = trueMoments[k];
      return Math.abs(recon - trueVal);
    });
  }, [reconstructedMoments, trueMoments]);

  const maxError = errors ? Math.max(...errors) : 0;

  if (isLoading) {
    return (
      <div className="text-xs text-purple-600 italic">
        Computing Boys function values...
      </div>
    );
  }

  if (error) {
    return (
      <div className="text-xs text-red-600">
        Error: {error}
      </div>
    );
  }

  return (
    <div>
      <p className="text-xs text-purple-700 mb-2">
        Verification: <MathFormula>{'\\sum_i w_i t_i^m'}</MathFormula> should equal{' '}
        <MathFormula>{'\\mu_m(T) = 2F_m(T)'}</MathFormula>
      </p>
      <div className="overflow-x-auto">
        <table className="w-full text-xs border-collapse">
          <thead>
            <tr className="border-b border-purple-200">
              <th className="text-left py-1 px-2 text-purple-800">
                <MathFormula>{'m'}</MathFormula>
              </th>
              <th className="text-right py-1 px-2 text-purple-800">
                <MathFormula>{'\\sum_i w_i t_i^m'}</MathFormula>
              </th>
              {trueMoments && (
                <>
                  <th className="text-right py-1 px-2 text-purple-800">
                    <MathFormula>{'2F_m(T)'}</MathFormula>
                  </th>
                  <th className="text-right py-1 px-2 text-purple-800">Error</th>
                </>
              )}
            </tr>
          </thead>
          <tbody className="text-purple-700 font-mono">
            {reconstructedMoments.map((recon, k) => (
              <tr key={k} className={k === 0 ? 'bg-purple-100' : ''}>
                <td className="py-0.5 px-2">{k}</td>
                <td className="text-right py-0.5 px-2">{recon.toExponential(8)}</td>
                {trueMoments && errors && (
                  <>
                    <td className="text-right py-0.5 px-2">{trueMoments[k].toExponential(8)}</td>
                    <td className="text-right py-0.5 px-2">
                      <span className={errors[k] > 1e-10 ? 'text-amber-600' : 'text-green-600'}>
                        {errors[k].toExponential(2)}
                      </span>
                    </td>
                  </>
                )}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      {errors && (
        <p className="text-xs text-purple-600 mt-2">
          Max reconstruction error:{' '}
          <span className={maxError > 1e-10 ? 'text-amber-600 font-medium' : 'text-green-600 font-medium'}>
            {maxError.toExponential(2)}
          </span>
          {maxError > 1e-10 && nroots > 6 && (
            <span className="text-amber-600"> (expected for n &gt; 6 due to Hankel ill-conditioning)</span>
          )}
        </p>
      )}
      <p className="text-xs text-purple-500 mt-1">
        Row <MathFormula>{'m=0'}</MathFormula> is highlighted: weight sum ={' '}
        <MathFormula>{'\\mu_0 = 2F_0(T)'}</MathFormula>
      </p>
    </div>
  );
}

/**
 * Component to display the Hankel matrix structure.
 */
function HankelMatrixDisplay({ nroots, t }: { nroots: number; t: number }) {
  const { send, isReady } = useWorker();
  const [moments, setMoments] = useState<number[] | null>(null);
  const [isLoading, setIsLoading] = useState(false);

  // Fetch moments when parameters change
  // Use useMemo for triggering the async effect
  useMemo(() => {
    if (!isReady) return;

    const fetchMoments = async () => {
      setIsLoading(true);
      try {
        // Make individual boys_eval calls for m = 0, 1, ..., 2*nroots-1
        // Use Promise.all for parallel execution
        const maxOrder = 2 * nroots - 1;
        const promises: Promise<BoysEvalResult>[] = [];
        for (let order = 0; order <= maxOrder; order++) {
          promises.push(
            send<BoysEvalResult>({
              type: 'boys_eval',
              m: order,
              T: t,
            } as Omit<BoysEvalRequest, 'requestId'>)
          );
        }
        const results = await Promise.all(promises);
        // Convert Boys values to moments: mu_k = 2 * F_k(T)
        setMoments(results.map((r) => 2 * r.value));
      } catch {
        // Silently fail - display symbolic form
      } finally {
        setIsLoading(false);
      }
    };

    fetchMoments();
  }, [isReady, send, t, nroots]);

  // Only show numerical values for small matrices
  const showNumerical = moments && nroots <= 4;
  const displaySize = Math.min(nroots, 4);

  return (
    <div className="space-y-3">
      <p className="text-xs text-slate-600">
        The Hankel matrix <MathFormula>{'H'}</MathFormula> is the Gram matrix of monomials{' '}
        <MathFormula>{'\\{1, x, x^2, \\ldots\\}'}</MathFormula> under the Rys inner product.
        Its entries <MathFormula>{'H_{ij} = \\mu_{i+j}'}</MathFormula> are constant along anti-diagonals.
      </p>

      {/* Symbolic form */}
      <div className="bg-white rounded p-2 border border-slate-200 overflow-x-auto">
        <MathDisplay>
          {`H = \\begin{pmatrix} \\mu_0 & \\mu_1 & \\cdots & \\mu_{${nroots - 1}} \\\\ \\mu_1 & \\mu_2 & \\cdots & \\mu_{${nroots}} \\\\ \\vdots & \\vdots & \\ddots & \\vdots \\\\ \\mu_{${nroots - 1}} & \\mu_{${nroots}} & \\cdots & \\mu_{${2 * nroots - 2}} \\end{pmatrix}`}
        </MathDisplay>
      </div>

      {/* Numerical values for small matrices */}
      {showNumerical && !isLoading && (
        <div className="bg-blue-50 rounded p-2 border border-blue-200">
          <p className="text-xs text-blue-700 mb-2 font-medium">
            Numerical values at T = {t.toFixed(4)}:
          </p>
          <div className="overflow-x-auto">
            <table className="text-xs font-mono">
              <tbody>
                {Array.from({ length: displaySize }, (_, i) => (
                  <tr key={i}>
                    {Array.from({ length: displaySize }, (_, j) => {
                      const moment = moments[i + j];
                      return (
                        <td key={j} className="px-2 py-0.5 text-right text-blue-800">
                          {moment !== undefined ? moment.toFixed(4) : '...'}
                        </td>
                      );
                    })}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Shifted Hankel matrix */}
      <div className="text-xs text-slate-600">
        <p>
          The shifted Hankel matrix <MathFormula>{'H^{(1)}'}</MathFormula> has entries{' '}
          <MathFormula>{'H^{(1)}_{ij} = \\mu_{i+j+1}'}</MathFormula>, representing
          multiplication by <MathFormula>{'x'}</MathFormula> in the monomial basis.
        </p>
      </div>

      {/* Properties */}
      <div className="bg-slate-50 rounded p-2 text-xs">
        <p className="font-medium text-slate-700 mb-1">Key Properties:</p>
        <ul className="list-disc list-inside text-slate-600 space-y-0.5">
          <li>
            <MathFormula>{'H'}</MathFormula> is symmetric positive definite
          </li>
          <li>
            Anti-diagonal structure: <MathFormula>{'H_{ij}'}</MathFormula> depends only on{' '}
            <MathFormula>{'i+j'}</MathFormula>
          </li>
          <li>
            Dimension: <MathFormula>{`${nroots} \\times ${nroots}`}</MathFormula>
          </li>
        </ul>
      </div>
    </div>
  );
}

/**
 * Component to display Cholesky factorization explanation.
 */
function CholeskySection() {
  return (
    <div className="space-y-3">
      <p className="text-xs text-slate-600">
        Cholesky factorization <MathFormula>{'H = LL^T'}</MathFormula> is equivalent to
        Gram-Schmidt orthonormalization of the monomial basis.
      </p>

      <MathDisplay>
        {`H = LL^T, \\quad C = L^{-1}`}
      </MathDisplay>

      <div className="bg-green-50 rounded p-3 border border-green-200">
        <p className="text-xs text-green-700">
          The transformation matrix <MathFormula>{'C = L^{-1}'}</MathFormula> converts
          monomials to orthonormal polynomials:
        </p>
        <MathDisplay>
          {`p_k(x) = \\sum_{j=0}^{k} C_{kj} x^j, \\quad \\langle p_i, p_j \\rangle_T = \\delta_{ij}`}
        </MathDisplay>
        <p className="text-xs text-green-600 mt-2">
          This is the same Gram-Schmidt procedure used for AO basis orthonormalization
          in SCF, but applied to polynomials instead of atomic orbitals.
        </p>
      </div>

      <div className="text-xs text-slate-600">
        <p className="font-medium mb-1">Why Cholesky?</p>
        <ul className="list-disc list-inside space-y-0.5">
          <li>Numerically stable for symmetric positive definite matrices</li>
          <li>Lower triangular structure preserves polynomial degree ordering</li>
          <li>
            Equivalent to Gram-Schmidt but encoded as matrix factorization
          </li>
        </ul>
      </div>
    </div>
  );
}

/**
 * Component to display Jacobi matrix structure.
 */
function JacobiMatrixDisplay({ nroots }: { nroots: number }) {
  return (
    <div className="space-y-3">
      <p className="text-xs text-slate-600">
        The Jacobi matrix <MathFormula>{'J'}</MathFormula> represents &quot;multiply by{' '}
        <MathFormula>{'x'}</MathFormula>&quot; in the orthonormal polynomial basis.
        It is symmetric tridiagonal:
      </p>

      <MathDisplay>
        {`J = \\begin{pmatrix} a_0 & b_1 & 0 & \\cdots & 0 \\\\ b_1 & a_1 & b_2 & \\cdots & 0 \\\\ 0 & b_2 & a_2 & \\ddots & \\vdots \\\\ \\vdots & \\ddots & \\ddots & \\ddots & b_{${nroots - 1}} \\\\ 0 & \\cdots & 0 & b_{${nroots - 1}} & a_{${nroots - 1}} \\end{pmatrix}`}
      </MathDisplay>

      <div className="bg-purple-50 rounded p-3 border border-purple-200">
        <p className="text-xs text-purple-700 font-medium mb-1">
          Three-Term Recurrence Coefficients:
        </p>
        <ul className="text-xs text-purple-600 space-y-0.5">
          <li>
            <MathFormula>{'a_k'}</MathFormula> (diagonal): &quot;centering&quot; of polynomial{' '}
            <MathFormula>{'p_k'}</MathFormula>
          </li>
          <li>
            <MathFormula>{'b_k > 0'}</MathFormula> (off-diagonal): coupling to neighbors
          </li>
        </ul>
        <MathDisplay>
          {`x\\,p_k(x) = b_{k+1}\\,p_{k+1}(x) + a_k\\,p_k(x) + b_k\\,p_{k-1}(x)`}
        </MathDisplay>
      </div>

      <div className="text-xs text-slate-600">
        <p className="font-medium mb-1">Physical Analogy:</p>
        <p>
          This structure is analogous to the tight-binding model in solid-state physics:
          each &quot;site&quot; (polynomial <MathFormula>{'p_k'}</MathFormula>) couples only
          to its nearest neighbors through the operator &quot;multiply by{' '}
          <MathFormula>{'x'}</MathFormula>&quot;.
        </p>
      </div>

      {/* Connection to nodes/weights */}
      <div className="bg-amber-50 rounded p-2 border border-amber-200 text-xs">
        <p className="text-amber-800">
          <span className="font-medium">Golub-Welsch (theory):</span> Eigendecomposition of{' '}
          <MathFormula>{'J'}</MathFormula> gives nodes and weights.
          The implementation achieves the same result via companion matrix QR of the
          orthogonal polynomials built by Schmidt orthogonalization.
        </p>
      </div>
    </div>
  );
}

/**
 * Component to display numerical stability information.
 */
function NumericalStabilitySection({ nroots }: { nroots: number }) {
  const conditionNumber = HANKEL_CONDITION_NUMBERS[nroots] || '> 10^{16}';

  return (
    <div className="space-y-3">
      <p className="text-xs text-slate-600">
        The Hankel matrix becomes increasingly ill-conditioned as{' '}
        <MathFormula>{'n'}</MathFormula> grows. The condition number{' '}
        <MathFormula>{'\\kappa(H)'}</MathFormula> grows exponentially.
      </p>

      {/* Condition number table */}
      <div className="overflow-x-auto">
        <table className="text-xs border-collapse w-full">
          <thead>
            <tr className="border-b border-slate-300">
              <th className="text-center py-1 px-2 text-slate-700">
                <MathFormula>{'n'}</MathFormula>
              </th>
              <th className="text-center py-1 px-2 text-slate-700">
                <MathFormula>{'\\kappa(H)'}</MathFormula> (approx)
              </th>
              <th className="text-center py-1 px-2 text-slate-700">Status</th>
            </tr>
          </thead>
          <tbody className="text-slate-600">
            {[2, 3, 4, 5, 6, 7, 8].map((n) => {
              const kappa = HANKEL_CONDITION_NUMBERS[n];
              const isCurrentN = n === nroots;
              const isBeyondPrecision = n > 6;
              return (
                <tr
                  key={n}
                  className={`border-b border-slate-200 ${isCurrentN ? 'bg-blue-100' : ''} ${isBeyondPrecision ? 'text-red-600' : ''}`}
                >
                  <td className="text-center py-0.5 px-2">{n}</td>
                  <td className="text-center py-0.5 px-2 font-mono">
                    <MathFormula>{kappa}</MathFormula>
                  </td>
                  <td className="text-center py-0.5 px-2">
                    {n <= 5 ? (
                      <span className="text-green-600">Safe</span>
                    ) : n === 6 ? (
                      <span className="text-amber-600">Borderline</span>
                    ) : (
                      <span className="text-red-600">Unreliable</span>
                    )}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      {/* Current status */}
      <div
        className={`rounded p-3 text-xs ${nroots > 6 ? 'bg-red-50 border border-red-200' : 'bg-green-50 border border-green-200'}`}
      >
        <p className={nroots > 6 ? 'text-red-700' : 'text-green-700'}>
          <span className="font-medium">Current:</span> n = {nroots},{' '}
          <MathFormula>{`\\kappa(H) \\approx ${conditionNumber}`}</MathFormula>
        </p>
        {nroots > 6 && (
          <p className="text-red-600 mt-1">
            Beyond double precision reliability. Results may contain significant numerical errors.
          </p>
        )}
      </div>

      {/* T dependence */}
      <div className="text-xs text-slate-600 space-y-1">
        <p className="font-medium">Dependence on T:</p>
        <ul className="list-disc list-inside space-y-0.5">
          <li>
            <span className="font-medium">Small T:</span> Moments vary slowly, rows of H are
            nearly proportional, poor conditioning
          </li>
          <li>
            <span className="font-medium">Large T:</span> Moments decay rapidly, better
            conditioning, but moments need careful evaluation
          </li>
          <li>
            <span className="font-medium">Moderate T:</span> Generally best-conditioned, but
            still limited to n &lt;= 6
          </li>
        </ul>
      </div>

      {/* Production note */}
      <div className="bg-slate-100 rounded p-2 text-xs text-slate-600">
        <p className="font-medium mb-1">In Production Codes:</p>
        <p>
          Libraries like libcint use precomputed tables, asymptotic expansions, or extended
          precision to handle high angular momentum (n &gt; 6). The moment-based algorithm is
          excellent for teaching but has precision limits.
        </p>
      </div>
    </div>
  );
}

/**
 * Panel displaying computational internals for Rys quadrature.
 *
 * Shows:
 * - Weight function formula
 * - Algorithm 5.1 steps (RDK algorithm)
 * - Three-term recurrence formula
 * - Sanity checks (roots in (0,1), sorted, weights positive)
 * - Moment reconstruction verification with actual Boys values
 * - Hankel matrix visualization
 * - Cholesky factorization explanation
 * - Jacobi matrix structure
 * - Numerical stability information
 *
 * @example
 * ```tsx
 * <RysInternalsPanel result={rysResult} />
 * ```
 */
export function RysInternalsPanel({ result }: RysInternalsPanelProps) {
  const methodDetails = METHOD_DETAILS[result.method];

  // Run sanity checks
  const rootsInRange = checkRootsInRange(result.roots);
  const rootsSorted = checkRootsSorted(result.roots);
  const weightsPositive = checkWeightsPositive(result.weights);
  const weightSum = result.weights.reduce((sum, w) => sum + w, 0);

  return (
    <div className="bg-slate-50 rounded-lg p-4 border border-slate-200 space-y-4">
      {/* Numerical Stability Warning for n > 6 */}
      {result.nroots > 6 && (
        <div className="bg-amber-50 border-2 border-amber-300 rounded-lg p-4">
          <div className="flex items-start gap-3">
            <svg
              className="w-6 h-6 text-amber-500 flex-shrink-0 mt-0.5"
              fill="currentColor"
              viewBox="0 0 20 20"
            >
              <path
                fillRule="evenodd"
                d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z"
                clipRule="evenodd"
              />
            </svg>
            <div>
              <h4 className="text-amber-800 font-semibold mb-1">
                High Quadrature Order Warning
              </h4>
              <p className="text-amber-700 text-sm">
                For <MathFormula>{`n > 6`}</MathFormula>, the underlying moment structure
                (Hankel matrix) becomes ill-conditioned in double precision (condition number{' '}
                <MathFormula>{`\\kappa > 10^{12}`}</MathFormula>). The Schmidt orthogonalization
                and root-finding may lose accuracy.
              </p>
              <p className="text-amber-600 text-sm mt-1">
                Consider using <MathFormula>{`n \\leq 6`}</MathFormula> for reliable results
                in double precision.
              </p>
            </div>
          </div>
        </div>
      )}

      {/* Weight Function */}
      <div className="bg-blue-50 rounded-lg p-3 border border-blue-200">
        <h4 className="text-xs font-semibold text-blue-800 uppercase tracking-wide mb-2">
          Rys Weight Function
        </h4>
        <MathDisplay>
          {`w_T(x) = x^{-1/2} e^{-Tx}, \\quad x \\in [0,1]`}
        </MathDisplay>
        <p className="text-xs text-blue-700 mt-2">
          The moments of this weight function are related to Boys functions:{' '}
          <MathFormula>{'\\mu_n(T) = 2F_n(T)'}</MathFormula>
        </p>
      </div>

      {/* Algorithm Steps */}
      <div>
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-2">
          Algorithm Steps (RDK / libcint)
        </h4>
        <p className="text-xs text-slate-600 mb-2">{methodDetails.description}</p>
        <div className="bg-white rounded p-2 border border-slate-200">
          <ul className="text-xs text-slate-700 space-y-1">
            {methodDetails.steps.map((step, i) => (
              <li
                key={i}
                className={`pl-2 border-l-2 ${
                  i === 0
                    ? 'border-blue-400'
                    : i <= 1
                      ? 'border-green-400'
                      : 'border-purple-400'
                }`}
              >
                {step.text}
                {step.latex && (
                  <>
                    {' '}
                    <MathFormula>{step.latex}</MathFormula>
                  </>
                )}
              </li>
            ))}
          </ul>
        </div>
      </div>

      {/* Hankel Matrix (Collapsible) */}
      <CollapsibleSection title="Hankel Matrix (Theory)" defaultOpen={false}>
        <HankelMatrixDisplay nroots={result.nroots} t={result.t} />
      </CollapsibleSection>

      {/* Cholesky / Gram-Schmidt (Collapsible) */}
      <CollapsibleSection title="Gram-Schmidt / Cholesky (Theory)" defaultOpen={false}>
        <CholeskySection />
      </CollapsibleSection>

      {/* Jacobi Matrix (Collapsible) */}
      <CollapsibleSection title="Jacobi / Companion Matrix (Theory)" defaultOpen={false}>
        <JacobiMatrixDisplay nroots={result.nroots} />
      </CollapsibleSection>

      {/* Three-Term Recurrence */}
      <div className="bg-green-50 rounded-lg p-3 border border-green-200">
        <h4 className="text-xs font-semibold text-green-800 uppercase tracking-wide mb-2">
          Three-Term Recurrence
        </h4>
        <MathDisplay>
          {`x\\, p_k(x) = b_{k+1}\\, p_{k+1}(x) + a_k\\, p_k(x) + b_k\\, p_{k-1}(x)`}
        </MathDisplay>
        <p className="text-xs text-green-700 mt-2">
          Orthonormal polynomials <MathFormula>{'p_k'}</MathFormula> satisfy this recurrence.
          The coefficients <MathFormula>{'a_k, b_k'}</MathFormula> form the Jacobi matrix{' '}
          <MathFormula>{'J'}</MathFormula>.
        </p>
      </div>

      {/* Golub-Welsch Weight Formula */}
      <div className="bg-purple-50 rounded-lg p-3 border border-purple-200">
        <h4 className="text-xs font-semibold text-purple-800 uppercase tracking-wide mb-2">
          Weight Formula
        </h4>
        <MathDisplay>{`w_k = \\frac{1}{\\sum_{j=0}^{n-1} P_j(r_k)^2}`}</MathDisplay>
        <div className="text-xs text-purple-700 mt-2 space-y-1">
          <p>
            Weights are computed from the orthonormal polynomials <MathFormula>{'P_j'}</MathFormula> evaluated
            at each root <MathFormula>{'r_k'}</MathFormula>. Since <MathFormula>{'P_0 = 1/\\sqrt{\\mu_0}'}</MathFormula>,
            this is equivalent to the Golub-Welsch formula{' '}
            <MathFormula>{'w_i = \\mu_0 \\cdot (V_{0i})^2'}</MathFormula>.
          </p>
          <p className="text-purple-600 italic">
            All squared terms are positive, so weights are{' '}
            <span className="font-medium">strictly positive</span> for non-degenerate cases.
          </p>
        </div>
      </div>

      {/* Sanity Checks Section */}
      <div>
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-2">
          Sanity Checks
        </h4>
        <div className="space-y-2">
          {/* Roots in (0, 1) */}
          <div className="flex items-center justify-between text-sm">
            <span className="text-slate-600">Roots in (0, 1)</span>
            <span className="flex items-center gap-2">
              <StatusIcon valid={rootsInRange.valid} />
              <span className={rootsInRange.valid ? 'text-green-700' : 'text-red-700'}>
                {rootsInRange.count}/{rootsInRange.total}
              </span>
            </span>
          </div>

          {/* Roots sorted */}
          <div className="flex items-center justify-between text-sm">
            <span className="text-slate-600">Roots sorted</span>
            <span className="flex items-center gap-2">
              <StatusIcon valid={rootsSorted} />
              <span className={rootsSorted ? 'text-green-700' : 'text-red-700'}>
                {rootsSorted ? 'Yes' : 'No'}
              </span>
            </span>
          </div>

          {/* Weights positive */}
          <div className="flex items-center justify-between text-sm">
            <span className="text-slate-600">Weights positive</span>
            <span className="flex items-center gap-2">
              <StatusIcon valid={weightsPositive.valid} />
              <span className={weightsPositive.valid ? 'text-green-700' : 'text-red-700'}>
                {weightsPositive.count}/{weightsPositive.total}
              </span>
            </span>
          </div>
        </div>
      </div>

      {/* Moment Reconstruction Verification */}
      <div className="bg-purple-50 rounded-lg p-3 border border-purple-200">
        <h4 className="text-xs font-semibold text-purple-800 uppercase tracking-wide mb-2">
          Moment Reconstruction Verification
        </h4>
        <MomentReconstructionTable
          roots={result.roots}
          weights={result.weights}
          t={result.t}
          nroots={result.nroots}
        />
      </div>

      {/* Weight Sum Verification */}
      <div>
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-1">
          Weight Sum Verification
        </h4>
        <div className="flex items-center gap-2">
          <p className="text-sm text-slate-700">
            <MathFormula>{`\\sum_i w_i = ${weightSum.toExponential(12)}`}</MathFormula>
          </p>
          <StatusIcon valid={weightSum > 0} />
        </div>
        <p className="text-xs text-slate-500 mt-1">
          Should equal{' '}
          <MathFormula>{`\\mu_0(T) = 2F_0(${result.t.toFixed(4)})`}</MathFormula> (zeroth moment)
        </p>
      </div>

      {/* Numerical Stability Section (Collapsible) */}
      <CollapsibleSection title="Numerical Stability" defaultOpen={result.nroots > 6}>
        <NumericalStabilitySection nroots={result.nroots} />
      </CollapsibleSection>

      {/* Quadrature Details */}
      <div className="pt-2 border-t border-slate-200">
        <h4 className="text-xs font-semibold text-slate-500 uppercase tracking-wide mb-2">
          Quadrature Details
        </h4>
        <div className="text-xs text-slate-600 space-y-1">
          <div className="flex justify-between">
            <span>Order:</span>
            <span>
              <MathFormula>{`n = ${result.nroots}`}</MathFormula>
            </span>
          </div>
          <div className="flex justify-between">
            <span>Parameter:</span>
            <span>
              <MathFormula>{`T = ${result.t.toFixed(6)}`}</MathFormula>
            </span>
          </div>
          <div className="flex justify-between">
            <span>Polynomial exactness:</span>
            <span>
              <MathFormula>{`\\deg \\leq ${2 * result.nroots - 1}`}</MathFormula>
            </span>
          </div>
          <div className="flex justify-between">
            <span>Boys orders needed:</span>
            <span>
              <MathFormula>{`F_0 \\ldots F_{${2 * result.nroots - 1}}`}</MathFormula>
            </span>
          </div>
        </div>
      </div>
    </div>
  );
}

export default RysInternalsPanel;
