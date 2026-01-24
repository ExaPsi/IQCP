/**
 * RysResultDisplay - Display area for Rys quadrature computation results.
 *
 * Shows the computed roots/weights table, internals panel (in internals mode),
 * and error curve section with order recommendation.
 * Handles all computation states (idle, computing, success, error).
 *
 * @module components/rys/RysResultDisplay
 */

import { useState } from 'react';
import { useRysStore, findRecommendedOrder } from '../../stores/rysStore';
import { useRysErrorCurve } from '../../hooks/useRysErrorCurve';
import { RysDataTable } from './RysDataTable';
import { RysErrorChart } from './RysErrorChart';
import { RysInternalsPanel } from './RysInternalsPanel';
import { Math, MathDisplay } from '../common/Math';

/**
 * Loading spinner component.
 */
function LoadingSpinner() {
  return (
    <div className="flex items-center justify-center py-8">
      <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600" />
      <span className="ml-3 text-slate-600">Computing...</span>
    </div>
  );
}

/**
 * Format target accuracy for display.
 * Converts '1e-6' to LaTeX '10^{-6}'.
 */
function formatTargetLatex(target: string): string {
  const match = target.match(/1e(-?\d+)/);
  if (match) {
    return `10^{${match[1]}}`;
  }
  return target;
}

/**
 * Explanatory box about reconstruction error.
 * Teaches students what the error metric means and why it matters.
 */
function ReconstructionErrorExplanation() {
  return (
    <div className="bg-blue-50 rounded-lg p-4 border border-blue-200">
      <h4 className="font-medium text-blue-800 mb-2 flex items-center gap-2">
        <svg className="w-5 h-5" fill="currentColor" viewBox="0 0 20 20">
          <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
        </svg>
        What is Reconstruction Error?
      </h4>
      <p className="text-sm text-blue-700 mb-3">
        Rys quadrature computes roots <Math>{'t_i'}</Math> and weights <Math>{'w_i'}</Math> such
        that weighted sums <em>exactly reproduce</em> the Boys function moments. The reconstruction
        error measures how well this moment-matching works:
      </p>
      <MathDisplay>
        {`\\varepsilon_n = \\max_{0 \\le k \\le 2n-1} \\left| \\underbrace{\\sum_{i=1}^n w_i t_i^k}_{\\text{computed moment } \\hat{\\mu}_k} - \\underbrace{F_k(T)}_{\\text{true moment } \\mu_k} \\right|`}
      </MathDisplay>
      <p className="text-xs text-blue-600 mt-2">
        The error takes the <strong>maximum</strong> absolute difference over all <Math>{'2n'}</Math> moments
        that an <Math>{'n'}</Math>-point quadrature should reproduce exactly.
      </p>
    </div>
  );
}

/**
 * Explanation of exponential convergence for Gaussian quadrature.
 */
function ConvergenceTheoryNote() {
  return (
    <div className="bg-violet-50 rounded-lg p-4 border border-violet-200">
      <h4 className="font-medium text-violet-800 mb-2 flex items-center gap-2">
        <svg className="w-5 h-5" fill="currentColor" viewBox="0 0 20 20">
          <path d="M9 4.804A7.968 7.968 0 005.5 4c-1.255 0-2.443.29-3.5.804v10A7.969 7.969 0 015.5 14c1.669 0 3.218.51 4.5 1.385A7.962 7.962 0 0114.5 14c1.255 0 2.443.29 3.5.804v-10A7.968 7.968 0 0014.5 4c-1.255 0-2.443.29-3.5.804V12a1 1 0 11-2 0V4.804z" />
        </svg>
        Convergence Theory
      </h4>
      <p className="text-sm text-violet-700 mb-2">
        For well-conditioned problems, Gaussian quadrature errors decrease <strong>exponentially</strong> with order:
      </p>
      <MathDisplay>
        {`\\varepsilon_n = O(c^{-n}) \\quad \\text{for some } c > 1`}
      </MathDisplay>
      <p className="text-sm text-violet-700">
        This is why you see a roughly straight line on the log-scale chart: each additional
        quadrature point multiplies accuracy by approximately the same factor.
      </p>
      <div className="bg-amber-50 rounded p-2 mt-3 border border-amber-200">
        <p className="text-xs text-amber-700">
          <strong>Numerical Warning:</strong> For very high orders (<Math>{'n > 6'}</Math>),
          numerical conditioning may limit achievable accuracy. The Jacobi matrix eigenvalue
          problem becomes increasingly sensitive, and roundoff errors may prevent reaching
          theoretical precision.
        </p>
      </div>
    </div>
  );
}

/**
 * Explanation of the threshold line and what it represents.
 */
function ThresholdExplanation({ target }: { target: string }) {
  const targetLatex = formatTargetLatex(target);
  return (
    <div className="bg-slate-50 rounded-lg p-3 border border-slate-200 text-sm">
      <div className="flex items-start gap-2">
        <div className="w-6 h-0.5 bg-red-500 mt-2 border-dashed border-t-2 border-red-500 flex-shrink-0" />
        <div className="text-slate-700">
          <strong>Red dashed line</strong> = Target accuracy <Math>{targetLatex}</Math>.
          Points below this line have sufficient precision for your application.
        </div>
      </div>
    </div>
  );
}

/**
 * Physical meaning explanation - why this error matters.
 */
function PhysicalMeaningNote() {
  return (
    <div className="bg-emerald-50 rounded-lg p-4 border border-emerald-200">
      <h4 className="font-medium text-emerald-800 mb-2 flex items-center gap-2">
        <svg className="w-5 h-5" fill="currentColor" viewBox="0 0 20 20">
          <path fillRule="evenodd" d="M5.05 4.05a7 7 0 119.9 9.9L10 18.9l-4.95-4.95a7 7 0 010-9.9zM10 11a2 2 0 100-4 2 2 0 000 4z" clipRule="evenodd" />
        </svg>
        Physical Significance
      </h4>
      <p className="text-sm text-emerald-700">
        This reconstruction error <strong>directly affects molecular integral accuracy</strong>.
        When computing electron repulsion integrals (ERIs), the quadrature error propagates into
        the final integral values. A reconstruction error of <Math>{'10^{-10}'}</Math> typically
        yields integral errors of similar magnitude, which is sufficient for most quantum
        chemistry applications requiring <Math>{'10^{-8}'}</Math> Hartree energy precision.
      </p>
    </div>
  );
}

/**
 * Order recommendation badge component with enhanced explanation.
 */
function RecommendationBadge({ order, target }: { order: number | null; target: string }) {
  const targetLatex = formatTargetLatex(target);

  if (order === null) {
    return (
      <div className="bg-amber-50 border border-amber-200 rounded-lg px-4 py-3">
        <div className="flex items-center gap-2">
          <svg
            className="w-5 h-5 text-amber-500"
            fill="currentColor"
            viewBox="0 0 20 20"
          >
            <path
              fillRule="evenodd"
              d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z"
              clipRule="evenodd"
            />
          </svg>
          <span className="text-amber-800 font-medium">
            No order <Math>{'n'}</Math> (1-10) achieves <Math>{targetLatex}</Math> accuracy
          </span>
        </div>
        <p className="text-sm text-amber-700 mt-1">
          Try relaxing the accuracy target. For very small <Math>{'T'}</Math> values,
          higher orders may be needed, or numerical conditioning limits precision.
        </p>
      </div>
    );
  }

  return (
    <div className="bg-green-50 border border-green-200 rounded-lg px-4 py-3">
      <div className="flex items-center gap-2">
        <svg
          className="w-5 h-5 text-green-600"
          fill="currentColor"
          viewBox="0 0 20 20"
        >
          <path
            fillRule="evenodd"
            d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z"
            clipRule="evenodd"
          />
        </svg>
        <span className="text-green-800 font-medium">
          Recommended order: <Math>{`n = ${order}`}</Math>
        </span>
      </div>
      <p className="text-sm text-green-700 mt-2">
        This is the <strong>minimum order</strong> that achieves <Math>{targetLatex}</Math> accuracy
        or better at the current <Math>{'T'}</Math> value.
      </p>
      <p className="text-xs text-green-600 mt-1">
        <strong>Cost-accuracy trade-off:</strong> Using <Math>{`n = ${order}`}</Math> minimizes
        computational cost (fewer roots to compute) while meeting your precision requirement.
        Higher orders would work but require more floating-point operations.
      </p>
    </div>
  );
}

/**
 * Error curve accordion section with full pedagogical content.
 *
 * This section teaches students:
 * - What reconstruction error means (moment matching)
 * - The mathematical definition with proper formulas
 * - Convergence theory (exponential decay)
 * - Physical significance for molecular integrals
 * - How to interpret the chart and recommendation
 */
function ErrorCurveSection() {
  const target = useRysStore((state) => state.target);
  const T = useRysStore((state) => state.T);
  const errorCurve = useRysStore((state) => state.errorCurve);
  const { isComputing, computeErrorCurve } = useRysErrorCurve();
  const [isExpanded, setIsExpanded] = useState(false);

  // Find recommended order if we have error curve data
  const recommendedOrder =
    errorCurve.status === 'success'
      ? findRecommendedOrder(errorCurve.result, target)
      : null;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      {/* Accordion header */}
      <button
        type="button"
        onClick={() => setIsExpanded(!isExpanded)}
        className="w-full px-6 py-4 flex items-center justify-between text-left hover:bg-slate-50 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-inset focus-visible:ring-blue-500"
        aria-expanded={isExpanded}
      >
        <div>
          <h3 className="text-lg font-semibold text-slate-800">
            Error Curve Analysis
          </h3>
          <p className="text-sm text-slate-500 mt-0.5">
            Explore how quadrature order affects moment reconstruction accuracy
          </p>
        </div>
        <svg
          className={`w-5 h-5 text-slate-500 transition-transform flex-shrink-0 ml-4 ${
            isExpanded ? 'transform rotate-180' : ''
          }`}
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
        >
          <path
            strokeLinecap="round"
            strokeLinejoin="round"
            strokeWidth={2}
            d="M19 9l-7 7-7-7"
          />
        </svg>
      </button>

      {/* Accordion content */}
      {isExpanded && (
        <div className="px-6 pb-6 space-y-4 border-t border-slate-100 pt-4">
          {/* Always show the explanation - this is the pedagogical core */}
          <ReconstructionErrorExplanation />

          {/* Compute button (shown when no data or to refresh) */}
          {(errorCurve.status === 'idle' || errorCurve.status === 'error') && (
            <div className="bg-slate-50 rounded-lg p-4 border border-slate-200">
              <p className="text-sm text-slate-600 mb-3">
                Click below to compute the error curve for <Math>{`T = ${T.toFixed(4)}`}</Math>,
                testing quadrature orders from <Math>{'n = 1'}</Math> to <Math>{'n = 10'}</Math>.
              </p>
              <button
                type="button"
                onClick={computeErrorCurve}
                disabled={isComputing}
                className="px-4 py-2 bg-blue-600 text-white rounded-lg font-medium hover:bg-blue-700 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
              >
                {isComputing ? 'Computing...' : 'Compute Error Curve'}
              </button>
            </div>
          )}

          {/* Computing state */}
          {errorCurve.status === 'computing' && <LoadingSpinner />}

          {/* Error state */}
          {errorCurve.status === 'error' && (
            <div className="bg-red-50 border border-red-200 rounded-lg p-4">
              <p className="text-red-800 font-medium">Computation Error</p>
              <p className="text-red-600 text-sm mt-1">{errorCurve.error}</p>
              <button
                type="button"
                onClick={computeErrorCurve}
                disabled={isComputing}
                className="mt-3 px-4 py-2 bg-red-600 text-white rounded-lg font-medium hover:bg-red-700 transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-red-500 disabled:opacity-50"
              >
                Try Again
              </button>
            </div>
          )}

          {/* Success state - show chart and all educational content */}
          {errorCurve.status === 'success' && (
            <>
              {/* Recommendation badge */}
              <RecommendationBadge order={recommendedOrder} target={target} />

              {/* The chart */}
              <div className="space-y-2">
                <RysErrorChart
                  errorCurve={errorCurve.result}
                  target={target}
                  loading={false}
                />
                {/* Chart legend explanation */}
                <ThresholdExplanation target={target} />
              </div>

              {/* Convergence theory explanation */}
              <ConvergenceTheoryNote />

              {/* Physical meaning */}
              <PhysicalMeaningNote />

              {/* Refresh button */}
              <div className="flex justify-between items-center pt-2 border-t border-slate-100">
                <p className="text-xs text-slate-500">
                  Computed at <Math>{`T = ${T.toFixed(4)}`}</Math>. Change <Math>{'T'}</Math> above and refresh to see how the curve changes.
                </p>
                <button
                  type="button"
                  onClick={computeErrorCurve}
                  disabled={isComputing}
                  className="text-sm text-blue-600 hover:text-blue-800 font-medium focus:outline-none focus-visible:underline flex items-center gap-1"
                >
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
                  </svg>
                  Refresh
                </button>
              </div>
            </>
          )}
        </div>
      )}
    </div>
  );
}

/**
 * Main result display component.
 */
function MainResultDisplay() {
  const compute = useRysStore((state) => state.compute);
  const mode = useRysStore((state) => state.mode);
  const n = useRysStore((state) => state.n);
  const T = useRysStore((state) => state.T);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 min-h-[200px]">
      <h2 className="text-lg font-semibold text-slate-800 mb-4">
        Roots & Weights
      </h2>

      {/* Idle state */}
      {compute.status === 'idle' && (
        <div className="text-center py-8 text-slate-500">
          <p className="text-lg">Waiting for computation...</p>
          <p className="text-sm mt-1">Adjust parameters to see results</p>
        </div>
      )}

      {/* Computing state */}
      {compute.status === 'computing' && <LoadingSpinner />}

      {/* Error state */}
      {compute.status === 'error' && (
        <div className="bg-red-50 border border-red-200 rounded-lg p-4">
          <p className="text-red-800 font-medium">Computation Error</p>
          <p className="text-red-600 text-sm mt-1">{compute.error}</p>
        </div>
      )}

      {/* Success state */}
      {compute.status === 'success' && (
        <div className="space-y-4">
          {/* Parameters summary */}
          <div className="bg-slate-50 rounded-lg p-3">
            <p className="text-sm text-slate-700">
              <Math>{`n = ${n}`}</Math> quadrature points at{' '}
              <Math>{`T = ${T.toFixed(4)}`}</Math>
            </p>
          </div>

          {/* Data table */}
          <RysDataTable result={compute.result} />

          {/* Conditional content based on display mode */}
          {mode === 'explain' ? (
            /* Explain mode - educational content */
            <div className="space-y-4">
              {/* Polynomial Exactness */}
              <div className="bg-blue-50 rounded-lg p-4 border border-blue-200">
                <h4 className="font-medium text-blue-800 mb-2">Polynomial Exactness</h4>
                <p className="text-sm text-blue-700">
                  With <Math>{`n = ${n}`}</Math> quadrature points, integrals of polynomials up to degree <Math>{`2n-1 = ${2*n - 1}`}</Math> are computed exactly:
                </p>
                <MathDisplay>
                  {`\\int p(x)\\, w_T(x)\\, dx = \\sum_i w_i\\, p(t_i), \\quad \\deg(p) \\leq ${2*n - 1}`}
                </MathDisplay>
              </div>

              {/* Root Count Rule Table */}
              <div className="bg-violet-50 rounded-lg p-4 border border-violet-200">
                <h4 className="font-medium text-violet-800 mb-2">Root Count Rule for Shell Quartets</h4>
                <p className="text-sm text-violet-700 mb-3">
                  The minimum number of Rys roots needed depends on the total angular momentum <Math>{'L = l_a + l_b + l_c + l_d'}</Math>:
                </p>
                <MathDisplay>
                  {`n_r = \\left\\lfloor \\frac{L}{2} \\right\\rfloor + 1`}
                </MathDisplay>
                <div className="overflow-x-auto">
                  <table className="w-full text-sm border-collapse">
                    <thead>
                      <tr className="border-b border-violet-200 text-violet-800">
                        <th className="text-left py-2 px-3">Shell Quartet</th>
                        <th className="text-center py-2 px-3"><Math>{'L'}</Math></th>
                        <th className="text-center py-2 px-3"><Math>{'n_r'}</Math></th>
                        <th className="text-left py-2 px-3">Boys Orders</th>
                      </tr>
                    </thead>
                    <tbody className="text-violet-700">
                      <tr className={n === 1 ? 'bg-violet-100' : ''}>
                        <td className="py-1.5 px-3"><Math>{'(ss|ss)'}</Math></td>
                        <td className="text-center py-1.5 px-3">0</td>
                        <td className="text-center py-1.5 px-3 font-medium">1</td>
                        <td className="py-1.5 px-3"><Math>{'F_0'}</Math></td>
                      </tr>
                      <tr className={n === 1 ? 'bg-violet-100' : ''}>
                        <td className="py-1.5 px-3"><Math>{'(ps|ss)'}</Math></td>
                        <td className="text-center py-1.5 px-3">1</td>
                        <td className="text-center py-1.5 px-3 font-medium">1</td>
                        <td className="py-1.5 px-3"><Math>{'F_0, F_1'}</Math></td>
                      </tr>
                      <tr className={n === 2 ? 'bg-violet-100' : ''}>
                        <td className="py-1.5 px-3"><Math>{'(pp|ss)'}</Math></td>
                        <td className="text-center py-1.5 px-3">2</td>
                        <td className="text-center py-1.5 px-3 font-medium">2</td>
                        <td className="py-1.5 px-3"><Math>{'F_0 \\ldots F_3'}</Math></td>
                      </tr>
                      <tr className={n === 3 ? 'bg-violet-100' : ''}>
                        <td className="py-1.5 px-3"><Math>{'(pp|pp)'}</Math></td>
                        <td className="text-center py-1.5 px-3">4</td>
                        <td className="text-center py-1.5 px-3 font-medium">3</td>
                        <td className="py-1.5 px-3"><Math>{'F_0 \\ldots F_5'}</Math></td>
                      </tr>
                      <tr className={n === 4 ? 'bg-violet-100' : ''}>
                        <td className="py-1.5 px-3"><Math>{'(dd|pp)'}</Math></td>
                        <td className="text-center py-1.5 px-3">6</td>
                        <td className="text-center py-1.5 px-3 font-medium">4</td>
                        <td className="py-1.5 px-3"><Math>{'F_0 \\ldots F_7'}</Math></td>
                      </tr>
                      <tr className={n === 5 ? 'bg-violet-100' : ''}>
                        <td className="py-1.5 px-3"><Math>{'(dd|dd)'}</Math></td>
                        <td className="text-center py-1.5 px-3">8</td>
                        <td className="text-center py-1.5 px-3 font-medium">5</td>
                        <td className="py-1.5 px-3"><Math>{'F_0 \\ldots F_9'}</Math></td>
                      </tr>
                      <tr className={n === 7 ? 'bg-violet-100' : ''}>
                        <td className="py-1.5 px-3"><Math>{'(ff|ff)'}</Math></td>
                        <td className="text-center py-1.5 px-3">12</td>
                        <td className="text-center py-1.5 px-3 font-medium">7</td>
                        <td className="py-1.5 px-3"><Math>{'F_0 \\ldots F_{13}'}</Math></td>
                      </tr>
                    </tbody>
                  </table>
                </div>
                <p className="text-xs text-violet-600 mt-2">
                  Row highlighted when current <Math>{'n'}</Math> matches that shell quartet's requirement.
                </p>
              </div>

              {/* About section */}
              <div className="text-sm text-slate-600 bg-slate-50 rounded-lg p-3">
                <span className="font-medium">Why This Matters: </span>
                These roots <Math>{'t_1, t_2, \\ldots, t_n'}</Math> and weights <Math>{'w_1, w_2, \\ldots, w_n'}</Math> allow exact integration of polynomials up to
                degree <Math>{`${2*n - 1}`}</Math> multiplied by the Rys weight function <Math>{'w_T(x) = x^{-1/2}e^{-Tx}'}</Math>.
                This is essential for computing electron-electron repulsion integrals (ERIs) in quantum chemistry,
                where the polynomial degree is determined by the angular momenta of the shell quartet.
              </div>
            </div>
          ) : (
            /* Internals mode - computational details */
            <RysInternalsPanel result={compute.result} />
          )}
        </div>
      )}
    </div>
  );
}

/**
 * Display area for Rys quadrature computation results.
 *
 * Shows:
 * - Loading state, error state, success with roots/weights table
 * - Display mode switching (explain vs internals)
 * - Error curve analysis section (collapsible)
 *
 * @example
 * ```tsx
 * <RysResultDisplay />
 * ```
 */
export function RysResultDisplay() {
  return (
    <div className="space-y-6">
      <MainResultDisplay />
      <ErrorCurveSection />
    </div>
  );
}

export default RysResultDisplay;
