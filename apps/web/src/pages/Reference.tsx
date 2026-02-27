/**
 * Reference page - Comprehensive technical documentation for quantum chemistry algorithms.
 * Content aligned with lecture notes covering Boys functions, Rys quadrature, and SCF/DIIS.
 */
import { useState } from 'react';
import { Math, MathBlock } from '@/components/common';

type SectionId = 'boys' | 'rys' | 'scf' | 'references';

interface CollapsibleSectionProps {
  title: string;
  id: SectionId;
  expandedSection: SectionId | null;
  onToggle: (id: SectionId) => void;
  children: React.ReactNode;
  color: 'blue' | 'green' | 'purple';
}

function CollapsibleSection({
  title,
  id,
  expandedSection,
  onToggle,
  children,
  color,
}: CollapsibleSectionProps) {
  const isExpanded = expandedSection === id;
  const colorClasses = {
    blue: {
      bg: 'bg-blue-50',
      border: 'border-blue-200',
      headerBg: 'bg-blue-100',
      headerText: 'text-blue-800',
      icon: 'text-blue-600',
    },
    green: {
      bg: 'bg-green-50',
      border: 'border-green-200',
      headerBg: 'bg-green-100',
      headerText: 'text-green-800',
      icon: 'text-green-600',
    },
    purple: {
      bg: 'bg-purple-50',
      border: 'border-purple-200',
      headerBg: 'bg-purple-100',
      headerText: 'text-purple-800',
      icon: 'text-purple-600',
    },
  };
  const colors = colorClasses[color];

  return (
    <div className={`rounded-xl border ${colors.border} overflow-hidden`}>
      <button
        onClick={() => onToggle(id)}
        className={`w-full px-6 py-4 flex items-center justify-between ${colors.headerBg} hover:opacity-90 transition-opacity`}
      >
        <h2 className={`text-xl font-semibold ${colors.headerText}`}>{title}</h2>
        <svg
          className={`w-6 h-6 ${colors.icon} transform transition-transform ${isExpanded ? 'rotate-180' : ''}`}
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>
      {isExpanded && <div className={`px-6 py-5 ${colors.bg}`}>{children}</div>}
    </div>
  );
}

interface SubsectionProps {
  title: string;
  children: React.ReactNode;
}

function Subsection({ title, children }: SubsectionProps) {
  return (
    <div className="mb-6 last:mb-0">
      <h3 className="text-lg font-semibold text-slate-700 mb-3">{title}</h3>
      {children}
    </div>
  );
}

function Reference() {
  const [expandedSection, setExpandedSection] = useState<SectionId | null>('boys');

  const handleToggle = (id: SectionId) => {
    setExpandedSection(expandedSection === id ? null : id);
  };

  return (
    <div className="max-w-4xl mx-auto">
      <h1 className="text-3xl font-bold text-slate-800 mb-4">Reference Documentation</h1>
      <p className="text-slate-600 mb-8">
        Technical reference for quantum chemistry algorithms implemented in IQCP. This
        documentation covers the mathematical foundations, computational methods, and key formulas
        for each interactive module.
      </p>

      <div className="space-y-4">
        {/* Module A: Boys Function */}
        <CollapsibleSection
          title="Module A: Boys Function F_m(T)"
          id="boys"
          expandedSection={expandedSection}
          onToggle={handleToggle}
          color="blue"
        >
          <Subsection title="Definition">
            <p className="text-slate-600 mb-3">
              The Boys function emerges as the central special function in the analytic evaluation
              of molecular integrals over Gaussian basis functions. For integer{' '}
              <Math>m</Math> and <Math>{'T \\geq 0'}</Math>:
            </p>
            <MathBlock label="Boys Function Definition">
              {'F_m(T) = \\int_0^1 t^{2m} e^{-Tt^2} \\, dt'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              The parameter <Math>{'T = \\rho \\cdot R_{PQ}^2'}</Math> measures the effective
              squared distance between two Gaussian overlap distributions. When{' '}
              <Math>{'T \\to 0'}</Math> (pair centers coincide),{' '}
              <Math>{'F_m(0) = \\frac{1}{2m+1}'}</Math>. When <Math>{'T \\to \\infty'}</Math>{' '}
              (pair centers far apart), <Math>{'F_m(T) \\to 0'}</Math>.
            </p>
          </Subsection>

          <Subsection title="Closed Form for F_0(T)">
            <MathBlock label="Error Function Form (T > 0)">
              {'F_0(T) = \\frac{\\sqrt{\\pi}}{2\\sqrt{T}} \\operatorname{erf}(\\sqrt{T})'}
            </MathBlock>
            <p className="text-slate-600">
              This provides the starting point for computing higher-order Boys functions via
              recurrence relations.
            </p>
          </Subsection>

          <Subsection title="Recurrence Relations">
            <p className="text-slate-600 mb-3">
              Integration by parts yields recurrence formulas for computing{' '}
              <Math>{'F_m'}</Math> from <Math>{'F_0'}</Math>:
            </p>
            <MathBlock label="Upward Recurrence (for T > 0)">
              {'F_{m+1}(T) = \\frac{(2m+1)F_m(T) - e^{-T}}{2T}'}
            </MathBlock>
            <MathBlock label="Downward Recurrence">
              {'F_m(T) = \\frac{2T \\cdot F_{m+1}(T) + e^{-T}}{2m+1}'}
            </MathBlock>
            <MathBlock label="Derivative Identity">
              {'\\frac{dF_m}{dT} = -F_{m+1}(T)'}
            </MathBlock>
          </Subsection>

          <Subsection title="Computational Regimes">
            <p className="text-slate-600 mb-4">
              The choice of evaluation method depends on the value of <Math>T</Math> to ensure
              numerical stability:
            </p>
            <div className="overflow-x-auto">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Regime
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Condition
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Method
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">Small T</td>
                    <td className="px-4 py-3 text-sm text-slate-700">
                      <Math>{'T < 25'}</Math>
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">Taylor Series Expansion</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">Moderate T</td>
                    <td className="px-4 py-3 text-sm text-slate-700">
                      <Math>{'25 \\leq T < 30 + 5m'}</Math>
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'\\operatorname{erf}(\\sqrt{T})'}</Math> + Upward Recurrence
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">Large T</td>
                    <td className="px-4 py-3 text-sm text-slate-700">
                      <Math>{'T \\geq 30 + 5m'}</Math>
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">Asymptotic Expansion</td>
                  </tr>
                </tbody>
              </table>
            </div>
          </Subsection>

          <Subsection title="Series and Asymptotic Expansions">
            <MathBlock label="Taylor Series (Small T)">
              {'F_m(T) = e^{-T} \\sum_{k=0}^{\\infty} \\frac{(2m-1)!! \\cdot T^k}{(2m+2k+1)!!}'}
            </MathBlock>
            <p className="text-slate-600 my-3">
              Truncating after approximately 15-25 terms typically achieves machine precision for{' '}
              <Math>{'T < 25'}</Math>.
            </p>
            <MathBlock label="Asymptotic Expansion (Large T)">
              {'F_m(T) \\approx \\frac{(2m-1)!!}{2^{m+1}} \\sqrt{\\frac{\\pi}{T^{2m+1}}}'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              where <Math>{'(2m-1)!! = 1 \\cdot 3 \\cdot 5 \\cdots (2m-1)'}</Math> is the double
              factorial with <Math>{'(-1)!! = 1'}</Math> by convention.
            </p>
          </Subsection>

          <Subsection title="Numerical Stability">
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
              <p className="text-amber-800 text-sm">
                <strong>Warning:</strong> The upward recurrence suffers from catastrophic
                cancellation when T is small. The numerator{' '}
                <Math>{'(2m+1)F_m(T) - e^{-T}'}</Math> subtracts nearly equal quantities. For
                small T, always use the series expansion or downward recurrence.
              </p>
            </div>
          </Subsection>

          <Subsection title="Connection to Incomplete Gamma Function">
            <MathBlock label="Boys-Gamma Relation">
              {'F_m(T) = \\frac{1}{2} T^{-(m+1/2)} \\gamma\\left(m + \\frac{1}{2}, T\\right)'}
            </MathBlock>
            <p className="text-slate-600">
              where <Math>{'\\gamma(a, x)'}</Math> is the lower incomplete gamma function. This
              relation is computationally valuable because many numerical libraries provide
              optimized routines for the incomplete gamma function.
            </p>
          </Subsection>
        </CollapsibleSection>

        {/* Module B: Rys Quadrature */}
        <CollapsibleSection
          title="Module B: Rys Quadrature"
          id="rys"
          expandedSection={expandedSection}
          onToggle={handleToggle}
          color="green"
        >
          <Subsection title="Purpose and Definition">
            <p className="text-slate-600 mb-3">
              Rys quadrature is a Gaussian quadrature scheme tailored to molecular integrals. It
              provides nodes <Math>{'{t_i}'}</Math> and weights <Math>{'{w_i}'}</Math> such that the
              weighted integral is exact for all polynomials up to degree{' '}
              <Math>{'2n_r - 1'}</Math>.
            </p>
            <MathBlock label="Rys Quadrature Formula">
              {'\\int_0^1 f(t^2) e^{-Tt^2} \\, dt \\approx \\sum_{i=1}^{n} w_i f(t_i^2)'}
            </MathBlock>
          </Subsection>

          <Subsection title="Moment-Boys Connection">
            <p className="text-slate-600 mb-3">
              The moments of the Rys weight function are directly related to Boys functions:
            </p>
            <MathBlock label="Moment Definition">
              {'\\mu_k(T) = \\int_0^1 t^{2k} e^{-Tt^2} \\, dt'}
            </MathBlock>
            <MathBlock label="Key Identity">{'\\mu_k(T) = 2 F_k(T)'}</MathBlock>
            <p className="text-slate-600 mt-3">
              This identity is the foundation of Rys quadrature: computing Boys functions is
              equivalent to computing moments of the Rys weight.
            </p>
          </Subsection>

          <Subsection title="Root Count Rule">
            <p className="text-slate-600 mb-3">
              For a shell quartet with total angular momentum{' '}
              <Math>{'L = l_A + l_B + l_C + l_D'}</Math>, the required number of Rys roots is:
            </p>
            <MathBlock label="Root Count Formula">
              {'n_r = \\left\\lfloor \\frac{L}{2} \\right\\rfloor + 1'}
            </MathBlock>
            <div className="overflow-x-auto mt-4">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Shell Quartet
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      <Math>L</Math>
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      <Math>n_r</Math>
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Boys Orders
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">(ss|ss)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">0</td>
                    <td className="px-4 py-3 text-sm text-slate-600">1</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>F_0</Math>
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">(pp|ss)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">2</td>
                    <td className="px-4 py-3 text-sm text-slate-600">2</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'F_0 \\text{ to } F_3'}</Math>
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">(pp|pp)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">4</td>
                    <td className="px-4 py-3 text-sm text-slate-600">3</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'F_0 \\text{ to } F_5'}</Math>
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">(dd|dd)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">8</td>
                    <td className="px-4 py-3 text-sm text-slate-600">5</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'F_0 \\text{ to } F_9'}</Math>
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">(ff|ff)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">12</td>
                    <td className="px-4 py-3 text-sm text-slate-600">7</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'F_0 \\text{ to } F_{13}'}</Math>
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
          </Subsection>

          <Subsection title="Quadrature Formula">
            <MathBlock label="Boys Function via Rys Quadrature">
              {'F_n(T) = \\frac{1}{2} \\sum_{i=1}^{n_r} w_i \\, t_i^{2n}, \\quad n = 0, 1, \\ldots, 2n_r - 1'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              This remarkable result shows that <Math>n_r</Math> Rys roots and weights exactly
              reproduce <Math>2n_r</Math> consecutive Boys function values.
            </p>
          </Subsection>

          <Subsection title="Golub-Welsch Algorithm">
            <p className="text-slate-600 mb-4">
              The algorithm computes quadrature nodes and weights from moments via a five-step
              pipeline:
            </p>
            <div className="bg-white rounded-lg border border-slate-200 p-4">
              <ol className="space-y-3 text-sm text-slate-700">
                <li className="flex items-start">
                  <span className="flex-shrink-0 w-6 h-6 bg-green-100 text-green-700 rounded-full flex items-center justify-center text-xs font-semibold mr-3">
                    1
                  </span>
                  <div>
                    <strong>Compute Moments:</strong>{' '}
                    <Math>{'\\mu_k = 2F_k(T)'}</Math> for{' '}
                    <Math>{'k = 0, \\ldots, 2n_r - 1'}</Math>
                  </div>
                </li>
                <li className="flex items-start">
                  <span className="flex-shrink-0 w-6 h-6 bg-green-100 text-green-700 rounded-full flex items-center justify-center text-xs font-semibold mr-3">
                    2
                  </span>
                  <div>
                    <strong>Form Hankel Matrices:</strong>{' '}
                    <Math>{'H_{ij} = \\mu_{i+j}'}</Math> and{' '}
                    <Math>{'H^{(1)}_{ij} = \\mu_{i+j+1}'}</Math>
                  </div>
                </li>
                <li className="flex items-start">
                  <span className="flex-shrink-0 w-6 h-6 bg-green-100 text-green-700 rounded-full flex items-center justify-center text-xs font-semibold mr-3">
                    3
                  </span>
                  <div>
                    <strong>Cholesky Factorization:</strong>{' '}
                    <Math>{'H = LL^T'}</Math>, then <Math>{'C = L^{-1}'}</Math>
                  </div>
                </li>
                <li className="flex items-start">
                  <span className="flex-shrink-0 w-6 h-6 bg-green-100 text-green-700 rounded-full flex items-center justify-center text-xs font-semibold mr-3">
                    4
                  </span>
                  <div>
                    <strong>Build Jacobi Matrix:</strong>{' '}
                    <Math>{'J = C H^{(1)} C^T'}</Math> (symmetric tridiagonal)
                  </div>
                </li>
                <li className="flex items-start">
                  <span className="flex-shrink-0 w-6 h-6 bg-green-100 text-green-700 rounded-full flex items-center justify-center text-xs font-semibold mr-3">
                    5
                  </span>
                  <div>
                    <strong>Eigendecomposition:</strong>{' '}
                    <Math>{'J = V \\Lambda V^T'}</Math>. Nodes:{' '}
                    <Math>{'t_i = \\Lambda_{ii}'}</Math>. Weights:{' '}
                    <Math>{'w_i = \\mu_0 (V_{0i})^2'}</Math>
                  </div>
                </li>
              </ol>
            </div>
          </Subsection>

          <Subsection title="Properties of Nodes and Weights">
            <div className="grid md:grid-cols-2 gap-4">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">
                  Nodes (<Math>t_i</Math>)
                </h4>
                <ul className="text-sm text-slate-600 space-y-1">
                  <li>
                    Must lie in the open interval <Math>{'(0, 1)'}</Math>
                  </li>
                  <li>Are eigenvalues of the Jacobi matrix</li>
                  <li>Are roots of orthogonal polynomials</li>
                  <li>Are distinct and well-separated</li>
                </ul>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">
                  Weights (<Math>w_i</Math>)
                </h4>
                <ul className="text-sm text-slate-600 space-y-1">
                  <li>Must be strictly positive</li>
                  <li>
                    Sum to <Math>{'\\mu_0 = 2F_0(T)'}</Math>
                  </li>
                  <li>Derived from eigenvector components</li>
                  <li>Ensure moment matching</li>
                </ul>
              </div>
            </div>
          </Subsection>

          <Subsection title="Verification: Moment Matching">
            <MathBlock label="Verification Condition">
              {'\\sum_{i=1}^{n_r} w_i \\, t_i^{2n} = \\mu_n, \\quad n = 0, 1, \\ldots, 2n_r - 1'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              This verification loop confirms that accumulated numerical errors remain within
              tolerance (typically <Math>{'10^{-12}'}</Math> or better for double precision with{' '}
              <Math>{'n_r \\leq 5'}</Math>).
            </p>
          </Subsection>

          <Subsection title="Numerical Limitations">
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
              <p className="text-amber-800 text-sm mb-2">
                <strong>Hankel Matrix Conditioning:</strong> The moment-based approach is limited to{' '}
                <Math>{'n_r \\leq 6'}</Math> in double precision due to ill-conditioning.
                Approximate condition numbers:
              </p>
              <ul className="text-amber-700 text-sm space-y-1 ml-4">
                <li>
                  <Math>{'n_r = 2'}</Math>: <Math>{'\\kappa(H) \\sim 10^2'}</Math>
                </li>
                <li>
                  <Math>{'n_r = 4'}</Math>: <Math>{'\\kappa(H) \\sim 10^6'}</Math>
                </li>
                <li>
                  <Math>{'n_r = 6'}</Math>: <Math>{'\\kappa(H) \\sim 10^{12}'}</Math>
                </li>
              </ul>
            </div>
          </Subsection>
        </CollapsibleSection>

        {/* Module C: SCF/DIIS */}
        <CollapsibleSection
          title="Module C: SCF and DIIS Acceleration"
          id="scf"
          expandedSection={expandedSection}
          onToggle={handleToggle}
          color="purple"
        >
          <Subsection title="Restricted Hartree-Fock (RHF) Overview">
            <p className="text-slate-600 mb-3">
              The Hartree-Fock method finds the best single Slater determinant wavefunction by
              solving a nonlinear eigenvalue problem iteratively. RHF applies to closed-shell
              systems where all electrons are paired.
            </p>
            <MathBlock label="Roothaan-Hall Equations">
              {'\\mathbf{F} \\mathbf{C} = \\mathbf{S} \\mathbf{C} \\boldsymbol{\\varepsilon}'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              where <Math>{'\\mathbf{F}'}</Math> is the Fock matrix,{' '}
              <Math>{'\\mathbf{C}'}</Math> contains MO coefficients,{' '}
              <Math>{'\\mathbf{S}'}</Math> is the overlap matrix, and{' '}
              <Math>{'\\boldsymbol{\\varepsilon}'}</Math> contains orbital energies. This
              generalized eigenvalue problem must be solved iteratively because{' '}
              <Math>{'\\mathbf{F}'}</Math> depends on <Math>{'\\mathbf{C}'}</Math> through the
              density matrix.
            </p>
          </Subsection>

          <Subsection title="Key Matrix Definitions">
            <div className="space-y-3">
              <MathBlock label="Density Matrix (closed-shell)">
                {'P_{\\mu\\nu} = 2 \\sum_i^{\\text{occ}} C_{\\mu i} C_{\\nu i}'}
              </MathBlock>
              <MathBlock label="Fock Matrix">
                {'F_{\\mu\\nu} = H_{\\mu\\nu}^{\\text{core}} + \\sum_{\\lambda\\sigma} P_{\\lambda\\sigma} \\left[ (\\mu\\nu|\\lambda\\sigma) - \\frac{1}{2}(\\mu\\lambda|\\nu\\sigma) \\right]'}
              </MathBlock>
              <MathBlock label="Coulomb Matrix">
                {'J_{\\mu\\nu} = \\sum_{\\lambda\\sigma} (\\mu\\nu|\\lambda\\sigma) P_{\\lambda\\sigma}'}
              </MathBlock>
              <MathBlock label="Exchange Matrix">
                {'K_{\\mu\\nu} = \\sum_{\\lambda\\sigma} (\\mu\\lambda|\\nu\\sigma) P_{\\lambda\\sigma}'}
              </MathBlock>
            </div>
          </Subsection>

          <Subsection title="SCF Iteration Algorithm">
            <div className="bg-white rounded-lg border border-slate-200 p-4">
              <ol className="space-y-2 text-sm text-slate-700">
                <li>
                  <strong>1. Initialize:</strong> Build <Math>{'\\mathbf{S}'}</Math>,{' '}
                  <Math>{'\\mathbf{H}^{\\text{core}}'}</Math>; compute orthogonalizer{' '}
                  <Math>{'\\mathbf{X} = \\mathbf{S}^{-1/2}'}</Math>; generate initial density{' '}
                  <Math>{'\\mathbf{P}^{(0)}'}</Math>
                </li>
                <li>
                  <strong>2. Build Fock matrix:</strong>{' '}
                  <Math>{'\\mathbf{F} = \\mathbf{H}^{\\text{core}} + \\mathbf{J}(\\mathbf{P}) - \\frac{1}{2}\\mathbf{K}(\\mathbf{P})'}</Math>
                </li>
                <li>
                  <strong>3. Transform:</strong>{' '}
                  <Math>{"\\mathbf{F}' = \\mathbf{X}^T \\mathbf{F} \\mathbf{X}"}</Math> (ordinary
                  eigenproblem)
                </li>
                <li>
                  <strong>4. Diagonalize:</strong>{' '}
                  <Math>{"\\mathbf{F}' \\mathbf{C}' = \\mathbf{C}' \\boldsymbol{\\varepsilon}"}</Math>
                </li>
                <li>
                  <strong>5. Back-transform:</strong>{' '}
                  <Math>{"\\mathbf{C} = \\mathbf{X} \\mathbf{C}'"}</Math>
                </li>
                <li>
                  <strong>6. Form new density:</strong>{' '}
                  <Math>{'\\mathbf{P}^{(\\text{new})}'}</Math> from occupied columns of{' '}
                  <Math>{'\\mathbf{C}'}</Math>
                </li>
                <li>
                  <strong>7. Check convergence:</strong> If{' '}
                  <Math>{'\\|\\mathbf{P}^{(\\text{new})} - \\mathbf{P}^{(\\text{old})}\\| < \\text{threshold}'}</Math>
                  , done
                </li>
                <li>
                  <strong>8. Update:</strong>{' '}
                  <Math>{'\\mathbf{P} = \\mathbf{P}^{(\\text{new})}'}</Math>, go to step 2
                </li>
              </ol>
            </div>
          </Subsection>

          <Subsection title="Energy Expression">
            <MathBlock label="Electronic Energy">
              {'E_{\\text{elec}} = \\frac{1}{2} \\text{Tr}[\\mathbf{P}(\\mathbf{H}^{\\text{core}} + \\mathbf{F})]'}
            </MathBlock>
            <MathBlock label="Total Energy">
              {'E_{\\text{total}} = E_{\\text{elec}} + V_{\\text{nn}}'}
            </MathBlock>
          </Subsection>

          <Subsection title="Orthogonalization">
            <p className="text-slate-600 mb-3">
              The Roothaan-Hall equations are a generalized eigenvalue problem because the AO basis
              is not orthonormal (<Math>{'\\mathbf{S} \\neq \\mathbf{I}'}</Math>). Symmetric
              (Lowdin) orthogonalization transforms to an ordinary eigenproblem:
            </p>
            <MathBlock label="Symmetric Orthogonalizer">
              {'\\mathbf{X} = \\mathbf{U} \\mathbf{s}^{-1/2} \\mathbf{U}^T = \\mathbf{S}^{-1/2}, \\quad \\text{where } \\mathbf{S} = \\mathbf{U} \\mathbf{s} \\mathbf{U}^T'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              The orthogonalizer <Math>{'\\mathbf{X}'}</Math> is computed once from{' '}
              <Math>{'\\mathbf{S}'}</Math> at the start of SCF and reused in every iteration.
            </p>
          </Subsection>

          <Subsection title="DIIS (Direct Inversion in the Iterative Subspace)">
            <p className="text-slate-600 mb-3">
              DIIS accelerates SCF convergence by constructing an improved Fock matrix as a linear
              combination of previous Fock matrices. Developed by Pulay (1980, 1982).
            </p>
            <MathBlock label="DIIS Extrapolation">
              {"\\mathbf{F}' = \\sum_{i=1}^{m} c_i \\mathbf{F}_i, \\quad \\text{with } \\sum_i c_i = 1"}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              The coefficients <Math>c_i</Math> are chosen to minimize the norm of the combined
              error vector.
            </p>
          </Subsection>

          <Subsection title="DIIS Error Vector">
            <MathBlock label="Commutator Error (Pulay)">
              {'\\mathbf{R} = \\mathbf{F} \\mathbf{P} \\mathbf{S} - \\mathbf{S} \\mathbf{P} \\mathbf{F}'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              The error matrix <Math>{'\\mathbf{R}'}</Math> vanishes when the Fock and density
              matrices commute under the metric defined by <Math>{'\\mathbf{S}'}</Math>, which
              occurs at the self-consistent solution.
            </p>
          </Subsection>

          <Subsection title="DIIS Linear System">
            <MathBlock label="Augmented System">
              {'\\begin{pmatrix} \\mathbf{B} & -\\mathbf{1} \\\\ -\\mathbf{1}^T & 0 \\end{pmatrix} \\begin{pmatrix} \\mathbf{c} \\\\ \\lambda \\end{pmatrix} = \\begin{pmatrix} \\mathbf{0} \\\\ -1 \\end{pmatrix}'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              where <Math>{'B_{ij} = \\mathbf{r}_i^T \\mathbf{r}_j'}</Math> (error inner products),
              and <Math>\\lambda</Math> is a Lagrange multiplier enforcing the constraint{' '}
              <Math>{'\\sum_i c_i = 1'}</Math>.
            </p>
          </Subsection>

          <Subsection title="Convergence Criteria">
            <div className="overflow-x-auto">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Profile
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Energy Threshold
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Density Threshold
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Loose</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-4}'}</Math> Ha
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-3}'}</Math>
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Medium</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-6}'}</Math> Ha
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-5}'}</Math>
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Tight</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-8}'}</Math> Ha
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-6}'}</Math>
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
          </Subsection>

          <Subsection title="DIIS Best Practices">
            <div className="grid md:grid-cols-2 gap-4">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">
                  Subspace Size (<Math>{'m_{\\max}'}</Math>)
                </h4>
                <ul className="text-sm text-slate-600 space-y-1">
                  <li>Typical range: 6-10 previous Fock matrices</li>
                  <li>Too small (&lt;4): Limited extrapolation power</li>
                  <li>Too large (&gt;15): B matrix becomes ill-conditioned</li>
                </ul>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Restart Conditions</h4>
                <ul className="text-sm text-slate-600 space-y-1">
                  <li>
                    Monitor condition number of <Math>{'\\mathbf{B}'}</Math> matrix
                  </li>
                  <li>
                    Clear history if <Math>{'\\kappa(\\mathbf{B})'}</Math> too large
                  </li>
                  <li>Residuals become collinear near convergence</li>
                </ul>
              </div>
            </div>
          </Subsection>

          <Subsection title="Additional Convergence Aids">
            <div className="space-y-3">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-1">Damping (Density Mixing)</h4>
                <p className="text-sm text-slate-600 mb-2">
                  <Math>
                    {'\\mathbf{P}^{(k+1)} = (1-\\alpha)\\mathbf{P}^{(k+1)}_{\\text{new}} + \\alpha\\mathbf{P}^{(k)}'}
                  </Math>
                  , where <Math>{'0 < \\alpha < 1'}</Math>
                </p>
                <p className="text-xs text-slate-500">
                  Prevents density from overshooting. Typical <Math>{'\\alpha = 0.3\\text{-}0.7'}</Math>.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-1">Level Shifting</h4>
                <p className="text-sm text-slate-600 mb-2">
                  Add constant <Math>b</Math> to virtual orbital energies:{' '}
                  <Math>{'\\varepsilon_a \\to \\varepsilon_a + b'}</Math>
                </p>
                <p className="text-xs text-slate-500">
                  Enlarges HOMO-LUMO gap. Typical <Math>{'b = 0.3\\text{-}1.0'}</Math> Ha. Useful for
                  near-metallic systems.
                </p>
              </div>
            </div>
          </Subsection>
        </CollapsibleSection>

        {/* References */}
        <CollapsibleSection
          title="Academic References"
          id="references"
          expandedSection={expandedSection}
          onToggle={handleToggle}
          color="blue"
        >
          <Subsection title="Boys Function and Molecular Integrals">
            <ul className="space-y-3 text-sm text-slate-600">
              <li className="border-l-4 border-blue-300 pl-4">
                <strong>Shavitt, I.</strong> (1963). The Gaussian Function in Calculations of
                Statistical Mechanics and Quantum Mechanics. In{' '}
                <em>Methods in Computational Physics</em>, Vol. 2, pp. 1-45.
                <p className="text-slate-500 text-xs mt-1">
                  Foundational reference for Boys function evaluation methods.
                </p>
              </li>
              <li className="border-l-4 border-blue-300 pl-4">
                <strong>Boys, S. F.</strong> (1950). Electronic Wave Functions. I. A General Method
                of Calculation for the Stationary States of Any Molecular System.{' '}
                <em>Proc. R. Soc. London A</em>, 200, 542-554.
                <p className="text-slate-500 text-xs mt-1">
                  Introduction of Gaussian-type orbitals in quantum chemistry.
                </p>
              </li>
            </ul>
          </Subsection>

          <Subsection title="Rys Quadrature">
            <ul className="space-y-3 text-sm text-slate-600">
              <li className="border-l-4 border-green-300 pl-4">
                <strong>Dupuis, M., Rys, J., & King, H. F.</strong> (1976). Evaluation of molecular
                integrals over Gaussian basis functions. <em>J. Chem. Phys.</em>, 65, 111-116.
                <p className="text-slate-500 text-xs mt-1">
                  Original paper introducing Rys quadrature for ERI evaluation.
                </p>
              </li>
              <li className="border-l-4 border-green-300 pl-4">
                <strong>Golub, G. H., & Welsch, J. H.</strong> (1969). Calculation of Gauss
                quadrature rules. <em>Math. Comp.</em>, 23, 221-230.
                <p className="text-slate-500 text-xs mt-1">
                  Algorithm for computing Gaussian quadrature from moments.
                </p>
              </li>
            </ul>
          </Subsection>

          <Subsection title="SCF and DIIS">
            <ul className="space-y-3 text-sm text-slate-600">
              <li className="border-l-4 border-purple-300 pl-4">
                <strong>Roothaan, C. C. J.</strong> (1951). New Developments in Molecular Orbital
                Theory. <em>Rev. Mod. Phys.</em>, 23, 69-89.
                <p className="text-slate-500 text-xs mt-1">
                  Foundation of the Roothaan-Hall equations for closed-shell HF.
                </p>
              </li>
              <li className="border-l-4 border-purple-300 pl-4">
                <strong>Pulay, P.</strong> (1980). Convergence acceleration of iterative sequences.
                The case of SCF iteration. <em>Chem. Phys. Lett.</em>, 73, 393-398.
                <p className="text-slate-500 text-xs mt-1">
                  Introduction of DIIS for SCF convergence acceleration.
                </p>
              </li>
              <li className="border-l-4 border-purple-300 pl-4">
                <strong>Pulay, P.</strong> (1982). Improved SCF convergence acceleration.{' '}
                <em>J. Comput. Chem.</em>, 3, 556-560.
                <p className="text-slate-500 text-xs mt-1">
                  Refinements to the original DIIS method.
                </p>
              </li>
            </ul>
          </Subsection>

          <Subsection title="Software References">
            <ul className="space-y-3 text-sm text-slate-600">
              <li className="border-l-4 border-slate-300 pl-4">
                <strong>Sun, Q.</strong> (2015). Libcint: An efficient general integral library for
                Gaussian basis functions. <em>J. Comput. Chem.</em>, 36, 1664-1671.
                <p className="text-slate-500 text-xs mt-1">
                  Reference for production-quality integral evaluation including Boys and Rys.
                </p>
              </li>
              <li className="border-l-4 border-slate-300 pl-4">
                <strong>Sun, Q. et al.</strong> (2018). PySCF: the Python-based simulations of
                chemistry framework. <em>WIREs Comput. Mol. Sci.</em>, 8, e1340.
                <p className="text-slate-500 text-xs mt-1">
                  Open-source quantum chemistry package used for validation.
                </p>
              </li>
            </ul>
          </Subsection>
        </CollapsibleSection>
      </div>

      {/* Quick Navigation */}
      <div className="mt-8 bg-slate-50 rounded-xl p-6 border border-slate-200">
        <h2 className="text-lg font-semibold text-slate-800 mb-4">Quick Navigation</h2>
        <div className="grid md:grid-cols-3 gap-4">
          <button
            onClick={() => setExpandedSection('boys')}
            className="text-left p-3 bg-white rounded-lg border border-slate-200 hover:border-blue-300 hover:bg-blue-50 transition-colors"
          >
            <div className="font-medium text-slate-700">Module A</div>
            <div className="text-sm text-slate-500">
              Boys Function <Math>F_m(T)</Math>
            </div>
          </button>
          <button
            onClick={() => setExpandedSection('rys')}
            className="text-left p-3 bg-white rounded-lg border border-slate-200 hover:border-green-300 hover:bg-green-50 transition-colors"
          >
            <div className="font-medium text-slate-700">Module B</div>
            <div className="text-sm text-slate-500">Rys Quadrature</div>
          </button>
          <button
            onClick={() => setExpandedSection('scf')}
            className="text-left p-3 bg-white rounded-lg border border-slate-200 hover:border-purple-300 hover:bg-purple-50 transition-colors"
          >
            <div className="font-medium text-slate-700">Module C</div>
            <div className="text-sm text-slate-500">SCF and DIIS</div>
          </button>
        </div>
      </div>
    </div>
  );
}

export default Reference;
