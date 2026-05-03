/**
 * Reference page - Comprehensive technical documentation for quantum chemistry algorithms.
 * Content aligned with lecture notes covering basis sets, integrals, Boys functions,
 * Rys quadrature, SCF/DIIS, DFT, geometry optimization, and population analysis.
 */
import { useState } from 'react';
import { Link } from 'react-router-dom';
import { Math, MathBlock } from '@/components/common';
import { ROUTES } from '@/lib/routes';

type SectionId = 'basis' | 'integrals' | 'boys' | 'rys' | 'scf' | 'frequency' | 'references';

interface CollapsibleSectionProps {
  title: React.ReactNode;
  id: SectionId;
  expandedSection: SectionId | null;
  onToggle: (id: SectionId) => void;
  children: React.ReactNode;
  color: 'blue' | 'green' | 'purple' | 'amber' | 'teal' | 'rose';
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
    amber: {
      bg: 'bg-amber-50',
      border: 'border-amber-200',
      headerBg: 'bg-amber-100',
      headerText: 'text-amber-800',
      icon: 'text-amber-600',
    },
    teal: {
      bg: 'bg-teal-50',
      border: 'border-teal-200',
      headerBg: 'bg-teal-100',
      headerText: 'text-teal-800',
      icon: 'text-teal-600',
    },
    rose: {
      bg: 'bg-rose-50',
      border: 'border-rose-200',
      headerBg: 'bg-rose-100',
      headerText: 'text-rose-800',
      icon: 'text-rose-600',
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
  const [expandedSection, setExpandedSection] = useState<SectionId | null>('basis');

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

      {/* Quick Navigation */}
      <div className="mb-8 bg-slate-50 rounded-xl p-6 border border-slate-200">
        <h2 className="text-lg font-semibold text-slate-800 mb-4">Quick Navigation</h2>
        <div className="grid md:grid-cols-3 gap-3">
          <div className="text-left p-3 bg-white rounded-lg border border-slate-200 hover:border-amber-300 hover:bg-amber-50 transition-colors">
            <button onClick={() => setExpandedSection('basis')} className="w-full text-left">
              <div className="font-medium text-slate-700">Module A</div>
              <div className="text-sm text-slate-500">Basis Explorer</div>
            </button>
            <Link to={ROUTES.BASIS} className="text-xs text-amber-600 hover:text-amber-800 mt-1 inline-block">
              Open Module &rarr;
            </Link>
          </div>
          <div className="text-left p-3 bg-white rounded-lg border border-slate-200 hover:border-teal-300 hover:bg-teal-50 transition-colors">
            <button onClick={() => setExpandedSection('integrals')} className="w-full text-left">
              <div className="font-medium text-slate-700">Module B</div>
              <div className="text-sm text-slate-500">Integral Inspector</div>
            </button>
            <Link to={ROUTES.INTEGRALS} className="text-xs text-teal-600 hover:text-teal-800 mt-1 inline-block">
              Open Module &rarr;
            </Link>
          </div>
          <div className="text-left p-3 bg-white rounded-lg border border-slate-200 hover:border-blue-300 hover:bg-blue-50 transition-colors">
            <button onClick={() => setExpandedSection('boys')} className="w-full text-left">
              <div className="font-medium text-slate-700">Module C</div>
              <div className="text-sm text-slate-500">
                Boys Function <Math>F_m(T)</Math>
              </div>
            </button>
            <Link to={ROUTES.BOYS} className="text-xs text-blue-600 hover:text-blue-800 mt-1 inline-block">
              Open Module &rarr;
            </Link>
          </div>
          <div className="text-left p-3 bg-white rounded-lg border border-slate-200 hover:border-green-300 hover:bg-green-50 transition-colors">
            <button onClick={() => setExpandedSection('rys')} className="w-full text-left">
              <div className="font-medium text-slate-700">Module D</div>
              <div className="text-sm text-slate-500">Rys Quadrature</div>
            </button>
            <Link to={ROUTES.RYS} className="text-xs text-green-600 hover:text-green-800 mt-1 inline-block">
              Open Module &rarr;
            </Link>
          </div>
          <div className="text-left p-3 bg-white rounded-lg border border-slate-200 hover:border-purple-300 hover:bg-purple-50 transition-colors">
            <button onClick={() => setExpandedSection('scf')} className="w-full text-left">
              <div className="font-medium text-slate-700">Module E</div>
              <div className="text-sm text-slate-500">SCF Sandbox</div>
            </button>
            <Link to={ROUTES.SCF} className="text-xs text-purple-600 hover:text-purple-800 mt-1 inline-block">
              Open Module &rarr;
            </Link>
          </div>
          <button
            onClick={() => setExpandedSection('frequency')}
            className="text-left p-3 bg-white rounded-lg border border-slate-200 hover:border-rose-300 hover:bg-rose-50 transition-colors"
          >
            <div className="font-medium text-slate-700">Frequency Analysis</div>
            <div className="text-sm text-slate-500">Vibrational Spectroscopy</div>
          </button>
          <button
            onClick={() => setExpandedSection('references')}
            className="text-left p-3 bg-white rounded-lg border border-slate-200 hover:border-blue-300 hover:bg-blue-50 transition-colors"
          >
            <div className="font-medium text-slate-700">References</div>
            <div className="text-sm text-slate-500">Academic Citations</div>
          </button>
        </div>
      </div>

      {/* Pipeline Overview */}
      <div className="bg-slate-50 rounded-xl p-6 border border-slate-200 mb-6">
        <h2 className="text-lg font-semibold text-slate-800 mb-3">How the Modules Connect</h2>
        <p className="text-slate-600 text-sm leading-relaxed">
          The IQCP modules form a natural progression through the computational chemistry
          pipeline: Module A explores how atoms are described by basis functions. Module B shows
          how these basis functions generate the integral matrices that encode molecular
          interactions. Modules C and D cover the special functions (Boys function and Rys
          quadrature) needed to evaluate these integrals efficiently. Module E brings everything
          together in the self-consistent field procedure, where the Hartree-Fock or Kohn-Sham
          equations are solved iteratively to find the electronic ground state. From the
          converged SCF wavefunction, Module E&apos;s Frequency tab computes the analytical
          Hessian via second derivatives and the CPHF equations, yielding normal mode
          frequencies, IR intensities, Raman activities, and thermochemical properties --
          connecting the electronic structure to experimentally observable spectra.
        </p>
      </div>

      <div className="space-y-4">
        {/* Module A: Basis Explorer */}
        <CollapsibleSection
          title="Module A: Basis Explorer"
          id="basis"
          expandedSection={expandedSection}
          onToggle={handleToggle}
          color="amber"
        >
          <Subsection title="Gaussian-Type Orbitals (GTOs)">
            <p className="text-slate-600 mb-3">
              Gaussian-type orbitals form the foundation of modern quantum chemistry calculations.
              Why Gaussians instead of the exact hydrogen-like orbitals (Slater-type orbitals)?
              The key advantage is that the product of two Gaussians centered on different atoms
              is itself a Gaussian centered at a point between them. This property makes all
              molecular integrals analytically tractable -- a feat impossible with Slater-type
              orbitals. A primitive Gaussian basis function is defined as:
            </p>
            <MathBlock label="Primitive Gaussian">
              {'g(\\mathbf{r}; \\alpha, \\mathbf{A}, l, m, n) = N \\cdot (x - A_x)^l (y - A_y)^m (z - A_z)^n \\exp(-\\alpha |\\mathbf{r} - \\mathbf{A}|^2)'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              where <Math>\alpha</Math> is the orbital exponent controlling the spatial extent,{' '}
              <Math>{'\\mathbf{A}'}</Math> is the center (usually a nucleus), <Math>N</Math> is a
              normalization constant ensuring <Math>{'\\int |g|^2 \\, d\\mathbf{r} = 1'}</Math>,
              and <Math>{'l + m + n'}</Math> determines the angular momentum type:
            </p>
            <div className="overflow-x-auto mt-3">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      <Math>{'l + m + n'}</Math>
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">Type</th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Components
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Count
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">0</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">s</td>
                    <td className="px-4 py-3 text-sm text-slate-600">1</td>
                    <td className="px-4 py-3 text-sm text-slate-600">1</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">1</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">p</td>
                    <td className="px-4 py-3 text-sm text-slate-600">x, y, z</td>
                    <td className="px-4 py-3 text-sm text-slate-600">3</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">2</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">d</td>
                    <td className="px-4 py-3 text-sm text-slate-600">xx, xy, xz, yy, yz, zz</td>
                    <td className="px-4 py-3 text-sm text-slate-600">6 (Cartesian)</td>
                  </tr>
                </tbody>
              </table>
            </div>
          </Subsection>

          <Subsection title="Contracted Gaussians">
            <p className="text-slate-600 mb-3">
              Individual primitive Gaussians are combined into contracted Gaussian functions to
              improve computational efficiency while retaining chemical accuracy:
            </p>
            <MathBlock label="Contracted Gaussian Function">
              {'\\chi_\\mu(\\mathbf{r}) = \\sum_{k=1}^{K} c_k \\, g_k(\\mathbf{r}; \\alpha_k)'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              where <Math>c_k</Math> are contraction coefficients and <Math>{`\\alpha_k`}</Math> are
              exponents, both predetermined and fixed during the SCF calculation.
              Contraction reduces the number of variational parameters while keeping the basis
              flexible enough to describe molecular wavefunctions accurately.
            </p>
          </Subsection>

          <Subsection title="Basis Set Families">
            <p className="text-slate-600 mb-3">
              IQCP supports six Gaussian basis sets spanning the Pople and Dunning families,
              covering elements H through Ar (Z = 1-18):
            </p>
            <div className="overflow-x-auto">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Basis Set
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Type
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Key Feature
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">STO-3G</td>
                    <td className="px-4 py-3 text-sm text-slate-600">Minimal</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      3 Gaussians approximate 1 Slater-type orbital; smallest possible basis
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">3-21G</td>
                    <td className="px-4 py-3 text-sm text-slate-600">Split-valence</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      Valence shell split into inner (2 primitives) and outer (1 primitive)
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">6-31G</td>
                    <td className="px-4 py-3 text-sm text-slate-600">Split-valence</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      Core: 6 primitives; valence: 3+1 split for better flexibility
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">6-31G*</td>
                    <td className="px-4 py-3 text-sm text-slate-600">Polarization</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      Adds d-functions on heavy atoms; essential for accurate bond angles
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">6-31+G*</td>
                    <td className="px-4 py-3 text-sm text-slate-600">Diffuse + Polarization</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      Adds diffuse sp-functions on heavy atoms; important for anions and lone pairs
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">cc-pVDZ</td>
                    <td className="px-4 py-3 text-sm text-slate-600">Correlation-consistent</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      Dunning double-zeta; designed for correlated methods; includes d-polarization on heavy atoms and p-polarization on H
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
          </Subsection>

          <Subsection title="Spherical vs Cartesian d-Orbitals">
            <p className="text-slate-600 mb-3">
              Cartesian d-functions have 6 components (<Math>{'xx, xy, xz, yy, yz, zz'}</Math>),
              while the true angular momentum <Math>{'l = 2'}</Math> subspace has only 5 linearly
              independent functions. The sixth Cartesian function is:
            </p>
            <MathBlock label="s-Type Contamination">
              {'x^2 + y^2 + z^2 = r^2 \\quad \\text{(transforms as } l = 0 \\text{)}'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              Spherical harmonics (<Math>{'d_{-2}, d_{-1}, d_0, d_{+1}, d_{+2}'}</Math>) are
              linear combinations of the 6 Cartesian components that form a pure{' '}
              <Math>{'l = 2'}</Math> set. Using spherical d-functions avoids s-type contamination
              and produces a smaller basis dimension.
            </p>

            <div className="overflow-x-auto mt-3 mb-3">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">Property</th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">Cartesian (6d)</th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">Spherical (5d)</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">d-functions per shell</td>
                    <td className="px-4 py-3 text-sm text-slate-600">6 (<Math>{'xx, xy, xz, yy, yz, zz'}</Math>)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">5 (<Math>{'d_{-2}, d_{-1}, d_0, d_{+1}, d_{+2}'}</Math>)</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">H₂O / 6-31G*</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">19 basis functions</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">18 basis functions</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">H₂O / cc-pVDZ</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">25 basis functions</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">24 basis functions</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Linear dependence</td>
                    <td className="px-4 py-3 text-sm text-slate-600">Contains s-type contamination</td>
                    <td className="px-4 py-3 text-sm text-slate-600">Pure angular momentum states</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Convention</td>
                    <td className="px-4 py-3 text-sm text-slate-600">PySCF with <code className="text-xs">mol.cart=True</code></td>
                    <td className="px-4 py-3 text-sm text-slate-600">PySCF default (<code className="text-xs">mol.cart=False</code>)</td>
                  </tr>
                </tbody>
              </table>
            </div>

            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
              <p className="text-blue-800 text-sm">
                <strong>IQCP implementation:</strong> IQCP supports both Cartesian (6d) and spherical (5d)
                harmonic basis functions, selectable via the Basis Function Type toggle in Module E.
                Integrals are evaluated in the Cartesian basis, then the Cartesian-to-spherical
                transformation (Schlegel &amp; Frisch, 1995) is applied to one-electron matrices,
                two-electron integrals, and DFT grid values. Both representations yield the same
                converged energy to sub-microhartree precision.
              </p>
            </div>
          </Subsection>

          <Subsection title="Basis Set Convergence">
            <p className="text-slate-600 mb-3">
              The variational principle guarantees that a larger basis set always gives an energy
              less than or equal to a smaller basis set. This monotonic convergence toward the
              Hartree-Fock limit is demonstrated by the following energies for{' '}
              <Math>{'\\text{H}_2\\text{O}'}</Math> (RHF):
            </p>
            <div className="overflow-x-auto">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Basis Set
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      N (Cartesian 6d)
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      <Math>{'E_{\\text{RHF}}'}</Math> (Ha)
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">STO-3G</td>
                    <td className="px-4 py-3 text-sm text-slate-600">7</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">-74.9630</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">3-21G</td>
                    <td className="px-4 py-3 text-sm text-slate-600">13</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">-75.5854</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">6-31G</td>
                    <td className="px-4 py-3 text-sm text-slate-600">13</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">-75.9840</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">6-31G*</td>
                    <td className="px-4 py-3 text-sm text-slate-600">19</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">-76.0105</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">6-31+G*</td>
                    <td className="px-4 py-3 text-sm text-slate-600">23</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">-76.0174</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">cc-pVDZ</td>
                    <td className="px-4 py-3 text-sm text-slate-600">25</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">-76.0271</td>
                  </tr>
                </tbody>
              </table>
            </div>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4 mt-3">
              <p className="text-amber-800 text-sm">
                <strong>Note:</strong> Energies computed at the IQCP H&#8322;O preset geometry
                (Cartesian 6d basis functions). Use Module E to compute values for your system
                of interest. The key observation is the monotonic decrease:{' '}
                <Math>{'E(\\text{STO-3G}) > E(\\text{3-21G}) > E(\\text{6-31G}) > E(\\text{6-31G*}) > E(\\text{6-31+G*}) > E(\\text{cc-pVDZ})'}</Math>.
              </p>
            </div>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4 mt-3">
              <p className="text-amber-800 text-sm">
                <strong>Common misconception:</strong> A larger basis set always gives better
                results. While the variational principle guarantees lower energies, other
                properties (like Mulliken charges) can become less physically meaningful with
                very large or diffuse basis sets.
              </p>
            </div>
          </Subsection>

          <Subsection title="Exponent Sliders">
            <p className="text-slate-600 mb-3">
              In the Basis Explorer, you can interactively adjust the exponent{' '}
              <Math>\alpha</Math> of a primitive Gaussian to see how the orbital shape changes in
              real time. The radial profile of a Gaussian basis function is:
            </p>
            <MathBlock label="Radial Profile">
              {'R(r) = N \\cdot r^l \\cdot \\exp(-\\alpha \\, r^2)'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              The exponent <Math>\alpha</Math> controls the spatial extent of the orbital.
              A large <Math>\alpha</Math> produces a tight, compact orbital concentrated
              near the nucleus, while a small <Math>\alpha</Math> produces a diffuse orbital
              that extends far from the center. This interactive control helps build intuition
              for why diffuse functions (as in 6-31+G*) are needed to describe anions and
              weakly bound electrons: their charge density extends further from the nucleus
              than standard valence functions can represent.
            </p>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4 mt-3">
              <p className="text-amber-800 text-sm">
                <strong>Try it:</strong> In <Link to={ROUTES.BASIS} className="text-amber-700 underline hover:text-amber-900">Module A</Link>, set a 1s orbital and sweep{' '}
                <Math>\alpha</Math> from 10 down to 0.1. Notice how the radial maximum
                moves outward and the orbital becomes broader. The integrated area under the
                curve always equals 1 (normalization), so spreading out reduces the peak height.
              </p>
            </div>
          </Subsection>

          <Subsection title="Overlap vs Distance">
            <p className="text-slate-600 mb-3">
              The overlap integral <Math>{'S_{ab}'}</Math> between two Gaussian basis functions
              measures how much they share the same region of space. For two s-type Gaussians
              centered at distance <Math>R</Math> apart:
            </p>
            <MathBlock label="Overlap Decay">
              {'S_{ab}(R) = \\left(\\frac{2\\sqrt{\\alpha_a \\alpha_b}}{\\alpha_a + \\alpha_b}\\right)^{3/2} \\exp\\!\\left(-\\frac{\\alpha_a \\alpha_b}{\\alpha_a + \\alpha_b} R^2\\right)'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              When <Math>{'R \\to 0'}</Math> (same center with identical exponents),{' '}
              <Math>{'S_{ab} \\to 1'}</Math>. When <Math>{'R \\to \\infty'}</Math>,{' '}
              <Math>{'S_{ab} \\to 0'}</Math> exponentially. The decay rate depends on the
              exponents: tighter functions (large <Math>\alpha</Math>) fall off faster with
              distance, while diffuse functions maintain appreciable overlap at longer range.
            </p>
            <p className="text-slate-600 mt-3">
              This connects basis function shape to molecular bonding: significant overlap
              between atomic orbitals on neighboring atoms is a prerequisite for strong covalent
              interaction. The overlap matrix <Math>{'\\mathbf{S}'}</Math> in an SCF calculation
              encodes all pairwise overlaps and is essential for orthogonalization.
            </p>
          </Subsection>

          <Subsection title="Basis Set Comparison Mode">
            <p className="text-slate-600 mb-3">
              The Basis Explorer allows overlaying radial profiles from multiple basis sets
              simultaneously, making it easy to see how different bases describe the same orbital.
              Each basis set is shown as a color-coded trace:
            </p>
            <ul className="text-sm text-slate-600 space-y-2 mb-3">
              <li>
                <strong>STO-3G:</strong> A single contracted Gaussian (3 primitives) approximating a
                Slater orbital. Shows only one peak per shell.
              </li>
              <li>
                <strong>3-21G / 6-31G:</strong> Split-valence bases show the hallmark two-component
                valence description: a tighter inner function and a more diffuse outer function.
                This splitting is directly visible as two distinct radial profiles for the valence
                shell.
              </li>
              <li>
                <strong>6-31G*:</strong> Adds d-type polarization functions on heavy atoms, visible
                as an additional set of radial curves with <Math>{'l = 2'}</Math> angular character.
              </li>
              <li>
                <strong>6-31+G*:</strong> Includes diffuse sp-functions with very small exponents,
                producing radial profiles that extend significantly further from the nucleus.
              </li>
              <li>
                <strong>cc-pVDZ:</strong> Dunning correlation-consistent design with generally
                contracted functions and p-polarization on hydrogen.
              </li>
            </ul>
            <p className="text-slate-600">
              Pedagogical discussion prompts in the module guide students through key observations:
              Why does split-valence splitting improve energies? How do polarization functions change
              bond angles? When are diffuse functions essential?
            </p>
          </Subsection>
        </CollapsibleSection>

        {/* Module B: Integral Inspector */}
        <CollapsibleSection
          title="Module B: Integral Inspector"
          id="integrals"
          expandedSection={expandedSection}
          onToggle={handleToggle}
          color="teal"
        >
          <Subsection title="One-Electron Integrals">
            <p className="text-slate-600 mb-3">
              One-electron integrals describe the interaction of an electron with the external
              potential (nuclei) and its kinetic energy. These matrices are computed once and
              remain constant throughout the SCF procedure.
            </p>
            <div className="space-y-3">
              <div>
                <MathBlock label="Overlap Integral">
                  {'S_{\\mu\\nu} = \\langle \\mu | \\nu \\rangle = \\int \\chi_\\mu(\\mathbf{r}) \\, \\chi_\\nu(\\mathbf{r}) \\, d\\mathbf{r}'}
                </MathBlock>
                <p className="text-slate-600 text-sm mt-2">
                  The overlap integral measures how much two basis functions occupy the same
                  region of space. If they are on distant atoms, <Math>{'S_{\\mu\\nu} \\approx 0'}</Math>.
                </p>
              </div>
              <div>
                <MathBlock label="Kinetic Energy Integral">
                  {'T_{\\mu\\nu} = \\left\\langle \\mu \\left| -\\frac{1}{2} \\nabla^2 \\right| \\nu \\right\\rangle'}
                </MathBlock>
                <p className="text-slate-600 text-sm mt-2">
                  The kinetic energy integral captures how much kinetic energy an electron has
                  when described by these basis functions.
                </p>
              </div>
              <div>
                <MathBlock label="Nuclear Attraction Integral">
                  {'V_{\\mu\\nu} = \\left\\langle \\mu \\left| -\\sum_A \\frac{Z_A}{|\\mathbf{r} - \\mathbf{R}_A|} \\right| \\nu \\right\\rangle'}
                </MathBlock>
                <p className="text-slate-600 text-sm mt-2">
                  The nuclear attraction integral measures the Coulomb attraction between an
                  electron in these orbitals and all nuclei.
                </p>
              </div>
              <div>
                <MathBlock label="Core Hamiltonian">
                  {'H^{\\text{core}}_{\\mu\\nu} = T_{\\mu\\nu} + V_{\\mu\\nu}'}
                </MathBlock>
                <p className="text-slate-600 text-sm mt-2">
                  The core Hamiltonian combines kinetic energy and nuclear attraction -- it
                  describes a single electron in the field of bare nuclei.
                </p>
              </div>
            </div>
            <p className="text-slate-600 mt-3">
              For an orthonormal basis,{' '}
              <Math>{'\\mathbf{S} = \\mathbf{I}'}</Math>, but Gaussian basis functions are
              generally non-orthogonal, making <Math>{'\\mathbf{S} \\neq \\mathbf{I}'}</Math>.
            </p>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 mt-3">
              <p className="text-blue-800 text-sm">
                <strong>Try it:</strong> Open <Link to={ROUTES.INTEGRALS} className="text-teal-700 underline hover:text-teal-900">Module B</Link> and select a molecule. Click on the overlap
                matrix heatmap to see how <Math>{'S_{\\mu\\nu}'}</Math> values relate to the
                spatial extent of each basis function pair.
              </p>
            </div>
          </Subsection>

          <Subsection title="Two-Electron Integrals (ERIs)">
            <p className="text-slate-600 mb-3">
              Electron repulsion integrals (ERIs) describe the Coulomb interaction between two
              electrons in the field of the basis functions. In chemist notation:
            </p>
            <MathBlock label="ERI in Chemist Notation">
              {'(\\mu\\nu|\\lambda\\sigma) = \\int\\!\\int \\chi_\\mu(\\mathbf{r}_1) \\chi_\\nu(\\mathbf{r}_1) \\frac{1}{r_{12}} \\chi_\\lambda(\\mathbf{r}_2) \\chi_\\sigma(\\mathbf{r}_2) \\, d\\mathbf{r}_1 \\, d\\mathbf{r}_2'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              For real basis functions, ERIs possess 8-fold permutational symmetry:
            </p>
            <MathBlock label="8-Fold Symmetry">
              {'(\\mu\\nu|\\lambda\\sigma) = (\\nu\\mu|\\lambda\\sigma) = (\\mu\\nu|\\sigma\\lambda) = (\\nu\\mu|\\sigma\\lambda) = (\\lambda\\sigma|\\mu\\nu) = \\ldots'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              For <Math>N</Math> basis functions, there are formally{' '}
              <Math>{'N^4'}</Math> ERIs, but 8-fold symmetry reduces the unique count to
              approximately <Math>{'N^4/8'}</Math>. This quartic scaling with basis size is the
              primary computational bottleneck in Hartree-Fock calculations.
            </p>
          </Subsection>

          <Subsection title="Schwarz Screening">
            <p className="text-slate-600 mb-3">
              Many ERIs are negligibly small. The Schwarz inequality provides an upper bound for
              efficient screening:
            </p>
            <MathBlock label="Schwarz Inequality">
              {'|(\\mu\\nu|\\lambda\\sigma)| \\leq \\sqrt{(\\mu\\nu|\\mu\\nu)} \\cdot \\sqrt{(\\lambda\\sigma|\\lambda\\sigma)}'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              The diagonal integrals <Math>{'(\\mu\\nu|\\mu\\nu)'}</Math> are computed once and
              stored. If the Schwarz bound falls below a threshold (typically{' '}
              <Math>{'10^{-12}'}</Math>), the full ERI is skipped. This reduces the effective
              scaling from <Math>{'O(N^4)'}</Math> toward <Math>{'O(N^2)'}</Math> for large
              molecules with localized basis functions.
            </p>
          </Subsection>

          <Subsection title="Fock Matrix Construction">
            <p className="text-slate-600 mb-3">
              The Fock matrix is assembled from the one-electron core Hamiltonian and two-electron
              contributions. The Coulomb (<Math>{'\\mathbf{J}'}</Math>) and exchange (
              <Math>{'\\mathbf{K}'}</Math>) matrices are built from ERIs and the density matrix:
            </p>
            <div className="space-y-3">
              <MathBlock label="Coulomb Matrix">
                {'J_{\\mu\\nu} = \\sum_{\\lambda\\sigma} P_{\\lambda\\sigma} \\, (\\mu\\nu|\\lambda\\sigma)'}
              </MathBlock>
              <MathBlock label="Exchange Matrix">
                {'K_{\\mu\\nu} = \\sum_{\\lambda\\sigma} P_{\\lambda\\sigma} \\, (\\mu\\lambda|\\nu\\sigma)'}
              </MathBlock>
              <MathBlock label="RHF Fock Matrix">
                {'F_{\\mu\\nu} = H^{\\text{core}}_{\\mu\\nu} + J_{\\mu\\nu} - \\frac{1}{2} K_{\\mu\\nu}'}
              </MathBlock>
            </div>
            <p className="text-slate-600 mt-3">
              The factor of <Math>{'\\frac{1}{2}'}</Math> before <Math>{'\\mathbf{K}'}</Math> in
              RHF arises from the closed-shell restriction. In the Integral Inspector module, you
              can trace exactly which ERIs contribute to each element of{' '}
              <Math>{'\\mathbf{J}'}</Math> and <Math>{'\\mathbf{K}'}</Math>.
            </p>
          </Subsection>

          <Subsection title="Interactive Heatmap Visualization">
            <p className="text-slate-600 mb-3">
              IQCP displays the overlap (<Math>{'\\mathbf{S}'}</Math>), kinetic energy (
              <Math>{'\\mathbf{T}'}</Math>), nuclear attraction (
              <Math>{'\\mathbf{V}'}</Math>), and core Hamiltonian (
              <Math>{'\\mathbf{H}^{\\text{core}}'}</Math>) matrices as color-coded heatmaps.
              Each cell is colored by magnitude, giving an immediate visual picture of the matrix
              structure.
            </p>
            <ul className="text-sm text-slate-600 space-y-2 mb-3">
              <li>
                Click any cell to see its numerical value and identify which pair of atomic
                orbitals (<Math>{'\\mu, \\nu'}</Math>) contribute to that element.
              </li>
              <li>
                Matrix symmetry (<Math>{'S_{\\mu\\nu} = S_{\\nu\\mu}'}</Math>) is immediately
                visible from the mirror symmetry of the heatmap pattern about the diagonal.
              </li>
              <li>
                Diagonal dominance of the overlap matrix (<Math>{'S_{\\mu\\mu} \\approx 1'}</Math>{' '}
                for normalized functions) appears as a bright band along the diagonal.
              </li>
              <li>
                Off-diagonal structure reveals which orbitals on different atoms have significant
                spatial overlap, directly connecting to bonding character.
              </li>
            </ul>
            <div className="bg-teal-50 border border-teal-200 rounded-lg p-4">
              <p className="text-teal-800 text-sm">
                <strong>Pedagogical note:</strong> Comparing the overlap heatmap across basis sets
                shows how split-valence and polarization functions add off-diagonal structure. A
                minimal basis like STO-3G produces a sparse overlap matrix; 6-31G* reveals richer
                inter-orbital interactions.
              </p>
            </div>
          </Subsection>

          <Subsection title="ERI Browser">
            <p className="text-slate-600 mb-3">
              A <em>shell</em> is a group of basis functions sharing the same exponents and
              contraction coefficients but differing in angular components (e.g., the three
              p-functions <Math>{'p_x, p_y, p_z'}</Math> form one p-shell). A <em>shell
              quartet</em> <Math>{'(ij|kl)'}</Math> indexes a block of ERIs between four
              shells.
            </p>
            <p className="text-slate-600 mb-3">
              The ERI Browser lets you explore electron repulsion integrals organized by shell
              quartet. This hierarchical view mirrors how ERIs are
              actually computed and stored in quantum chemistry codes.
            </p>
            <ul className="text-sm text-slate-600 space-y-2 mb-3">
              <li>
                <strong>Shell pair selection:</strong> Choose two shell pairs to see the full
                block of ERIs for that quartet, with magnitudes indicated by color intensity.
              </li>
              <li>
                <strong>Drill-down:</strong> Click any individual ERI value to see its
                primitive-pair decomposition — how each pair of contracted Gaussian primitives
                contributes to the total integral.
              </li>
              <li>
                <strong>Schwarz bounds:</strong> Each shell quartet displays its Schwarz screening
                bound <Math>{'\\sqrt{(\\mu\\nu|\\mu\\nu)} \\cdot \\sqrt{(\\lambda\\sigma|\\lambda\\sigma)}'}</Math>{' '}
                alongside the actual integral magnitude, showing which quartets would be screened
                out at a given threshold.
              </li>
              <li>
                <strong>Symmetry:</strong> The browser highlights the 8-fold permutational symmetry
                by grouping equivalent integrals and marking the unique representatives.
              </li>
            </ul>
          </Subsection>

          <Subsection title="Fock Matrix Decomposition">
            <p className="text-slate-600 mb-3">
              The Integral Inspector decomposes the Fock matrix into its constituent parts,
              allowing you to examine each contribution separately:
            </p>
            <div className="space-y-3 mb-3">
              <MathBlock label="HF Decomposition">
                {'\\mathbf{F} = \\underbrace{\\mathbf{H}^{\\text{core}}}_{\\text{1-electron}} + \\underbrace{\\mathbf{J}}_{\\text{Coulomb}} - \\underbrace{\\frac{1}{2}\\mathbf{K}}_{\\text{Exchange}}'}
              </MathBlock>
            </div>
            <p className="text-slate-600 mb-3">
              Each component is viewable as a separate heatmap. The core Hamiltonian{' '}
              <Math>{'\\mathbf{H}^{\\text{core}}'}</Math> captures single-electron physics
              (kinetic energy and nuclear attraction). The Coulomb matrix{' '}
              <Math>{'\\mathbf{J}'}</Math> represents the classical electron-electron repulsion,
              while the exchange matrix <Math>{'\\mathbf{K}'}</Math> encodes the quantum
              mechanical exchange interaction arising from the antisymmetry of the wavefunction.
              The exchange interaction has no classical analogue -- it arises from the quantum
              mechanical requirement that the electronic wavefunction be antisymmetric under
              electron exchange (the Pauli exclusion principle). Physically, electrons of the
              same spin tend to avoid each other, creating a &quot;Fermi hole&quot; that reduces
              their mutual repulsion and lowers the energy.
            </p>
            <p className="text-slate-600 mb-3">
              You can trace how specific ERIs contribute to particular{' '}
              <Math>{'\\mathbf{J}'}</Math> and <Math>{'\\mathbf{K}'}</Math> matrix elements,
              making the connection between four-index integrals and two-index matrices concrete.
            </p>
            <div className="bg-teal-50 border border-teal-200 rounded-lg p-4">
              <p className="text-teal-800 text-sm">
                <strong>DFT variant:</strong> For Kohn-Sham DFT calculations, the decomposition
                becomes <Math>{'\\mathbf{F}^{\\text{KS}} = \\mathbf{H}^{\\text{core}} + \\mathbf{J} + \\mathbf{V}^{\\text{XC}}'}</Math>{' '}
                (pure functionals like LDA) or includes a fraction of exact exchange for hybrid
                functionals like B3LYP. The XC potential matrix{' '}
                <Math>{'\\mathbf{V}^{\\text{XC}}'}</Math> replaces part or all of the exchange term.
              </p>
            </div>
          </Subsection>
        </CollapsibleSection>

        {/* Module C: Boys Function */}
        <CollapsibleSection
          title={<>Module C: Boys Function <Math>{'F_m(T)'}</Math></>}
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
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 mt-3">
              <p className="text-blue-800 text-sm">
                <strong>Try it:</strong> In <Link to={ROUTES.BOYS} className="text-blue-700 underline hover:text-blue-900">Module C</Link>, set <Math>{'m = 0'}</Math> and sweep{' '}
                <Math>T</Math> from 0 to 50 to see <Math>{'F_0(T)'}</Math> transition from 1 to
                0. Toggle the regime display to see which computational method is active at each{' '}
                <Math>T</Math> value.
              </p>
            </div>
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
            <p className="text-slate-600 text-sm mt-2 mb-3">
              This direction is numerically stable for large <Math>T</Math> because each
              successive <Math>{'F_m'}</Math> decreases in magnitude, so rounding errors do
              not amplify.
            </p>
            <MathBlock label="Downward Recurrence">
              {'F_m(T) = \\frac{2T \\cdot F_{m+1}(T) + e^{-T}}{2m+1}'}
            </MathBlock>
            <p className="text-slate-600 text-sm mt-2 mb-3">
              This direction is numerically stable for small <Math>T</Math> because the
              exponential term <Math>{'e^{-T}'}</Math> is close to 1 and the terms reinforce
              rather than cancel.
            </p>
            <MathBlock label="Derivative Identity">
              {'\\frac{dF_m}{dT} = -F_{m+1}(T)'}
            </MathBlock>
          </Subsection>

          <Subsection title="Computational Regimes">
            <p className="text-slate-600 mb-4">
              IQCP uses two computational regimes with an <Math>m</Math>-dependent turnover
              point derived from numerical error analysis (following libcint):
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
                    <td className="px-4 py-3 text-sm text-slate-600">Series</td>
                    <td className="px-4 py-3 text-sm text-slate-700">
                      <Math>{'T < \\text{turnover}(m)'}</Math>
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">Taylor series + downward recurrence</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">Recurrence</td>
                    <td className="px-4 py-3 text-sm text-slate-700">
                      <Math>{'T \\geq \\text{turnover}(m)'}</Math>
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'\\operatorname{erf}(\\sqrt{T})'}</Math> + upward recurrence
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
            <p className="text-slate-600 mt-3 mb-3">
              The turnover point is <Math>m</Math>-dependent and increases with order.
              Representative values:
            </p>
            <div className="overflow-x-auto">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      <Math>m</Math>
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Turnover <Math>T</Math>
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Notes
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">0, 1</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">0</td>
                    <td className="px-4 py-3 text-sm text-slate-600">Always use recurrence</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">2</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">~0.87</td>
                    <td className="px-4 py-3 text-sm text-slate-600"></td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">5</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">~2.1</td>
                    <td className="px-4 py-3 text-sm text-slate-600"></td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">10</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">~4.05</td>
                    <td className="px-4 py-3 text-sm text-slate-600"></td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-600">20</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">~7.84</td>
                    <td className="px-4 py-3 text-sm text-slate-600"></td>
                  </tr>
                </tbody>
              </table>
            </div>
          </Subsection>

          <Subsection title="Series Expansion and Limiting Form">
            <MathBlock label="Taylor Series (Small T)">
              {'F_m(T) = e^{-T} \\sum_{k=0}^{\\infty} \\frac{(2m-1)!! \\cdot T^k}{(2m+2k+1)!!}'}
            </MathBlock>
            <p className="text-slate-600 my-3">
              Truncating after approximately 15-25 terms typically achieves machine precision
              below the turnover point.
            </p>
            <MathBlock label="Large-T Limiting Form">
              {'F_m(T) \\sim \\frac{(2m-1)!!}{2^{m+1}} \\sqrt{\\frac{\\pi}{T^{2m+1}}}'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              where <Math>{'(2m-1)!! = 1 \\cdot 3 \\cdot 5 \\cdots (2m-1)'}</Math> is the double
              factorial with <Math>{'(-1)!! = 1'}</Math> by convention. This limiting form is
              useful for understanding the behavior of <Math>{'F_m(T)'}</Math> at large{' '}
              <Math>T</Math>, but is not used as a separate computational regime in IQCP. Instead,
              the erf-based recurrence formula handles all <Math>T</Math> values above the
              turnover point.
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

          <Subsection title="Implementation Details">
            <p className="text-slate-600 mb-3">
              IQCP uses <Math>m</Math>-dependent turnover points rather than fixed thresholds to
              switch between evaluation methods. These turnover points are derived from numerical
              error analysis, ensuring that each method operates in its region of optimal accuracy.
            </p>
            <div className="bg-white rounded-lg border border-slate-200 p-4 mb-3">
              <h4 className="font-semibold text-slate-700 mb-2">
                Direct Evaluation for <Math>{'m = 0'}</Math>
              </h4>
              <p className="text-sm text-slate-600">
                The zeroth-order Boys function is computed directly from the error function:
              </p>
              <MathBlock label="Direct Formula">
                {'F_0(T) = \\frac{\\sqrt{\\pi}}{2\\sqrt{T}} \\, \\operatorname{erf}(\\sqrt{T})'}
              </MathBlock>
              <p className="text-sm text-slate-600 mt-2">
                This closed-form expression serves as the anchor for computing higher orders.
              </p>
            </div>
            <div className="grid md:grid-cols-2 gap-4 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Below Turnover (Small T)</h4>
                <p className="text-sm text-slate-600">
                  Taylor series expansion with <strong>downward recurrence</strong>. Starting from
                  a high-order Boys function computed via the series, values are propagated downward
                  to lower orders. This direction is numerically stable because each step reinforces
                  rather than amplifies rounding errors.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Above Turnover (Large T)</h4>
                <p className="text-sm text-slate-600">
                  Starting from <Math>{'F_0(T)'}</Math> via the error function, then applying{' '}
                  <strong>upward recurrence</strong> to obtain higher orders. This direction is stable
                  for large <Math>T</Math> because <Math>{'F_m(T)'}</Math> decreases rapidly with{' '}
                  <Math>m</Math>, and the subtraction in the recurrence does not cause catastrophic
                  cancellation.
                </p>
              </div>
            </div>
            <p className="text-slate-600">
              The turnover point increases with <Math>m</Math>: higher-order Boys functions need
              larger <Math>T</Math> values before the recurrence approach becomes accurate.
              The turnover values in the regime table above are derived from error analysis
              of the two methods. IQCP visualizes the active regime for each computation, so
              students can see exactly which method is in use and why.
            </p>
          </Subsection>
        </CollapsibleSection>

        {/* Module D: Rys Quadrature */}
        <CollapsibleSection
          title="Module D: Rys Quadrature"
          id="rys"
          expandedSection={expandedSection}
          onToggle={handleToggle}
          color="green"
        >
          <Subsection title="Purpose and Definition">
            <p className="text-slate-600 mb-3">
              Evaluating molecular integrals requires computing integrals of the form{' '}
              <Math>{'\\int f(t^2) \\exp(-Tt^2) \\, dt'}</Math>. Rather than using
              general-purpose numerical integration, Rys quadrature exploits the specific
              structure of these integrals by finding optimal evaluation points (nodes) and
              weights tailored to the <Math>{'\\exp(-Tt^2)'}</Math> weight function.
            </p>
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
            <MathBlock label="Key Identity">{'\\mu_k(T) = F_k(T)'}</MathBlock>
            <p className="text-slate-600 mt-3">
              This identity is the foundation of Rys quadrature: computing Boys functions is
              equivalent to computing moments of the Rys weight function.
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
              {'F_n(T) = \\sum_{i=1}^{n_r} w_i \\, t_i^{2n}, \\quad n = 0, 1, \\ldots, 2n_r - 1'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              This remarkable result shows that <Math>n_r</Math> Rys roots and weights exactly
              reproduce <Math>2n_r</Math> consecutive Boys function values.
            </p>
          </Subsection>

          <Subsection title="Root-Finding Algorithm (Schmidt Orthogonalization)">
            <p className="text-slate-600 mb-4">
              IQCP uses the RDK algorithm (Rys, Dupuis, King) with Schmidt orthogonalization,
              following the <code className="text-xs">CINTrys_schmidt</code> implementation in
              libcint. The algorithm computes quadrature nodes and weights from moments via a
              four-step pipeline:
            </p>
            <div className="bg-white rounded-lg border border-slate-200 p-4">
              <ol className="space-y-3 text-sm text-slate-700">
                <li className="flex items-start">
                  <span className="flex-shrink-0 w-6 h-6 bg-green-100 text-green-700 rounded-full flex items-center justify-center text-xs font-semibold mr-3">
                    1
                  </span>
                  <div>
                    <strong>Compute Moments:</strong>{' '}
                    <Math>{'\\mu_k = F_k(T)'}</Math> for{' '}
                    <Math>{'k = 0, \\ldots, 2n_r'}</Math>
                  </div>
                </li>
                <li className="flex items-start">
                  <span className="flex-shrink-0 w-6 h-6 bg-green-100 text-green-700 rounded-full flex items-center justify-center text-xs font-semibold mr-3">
                    2
                  </span>
                  <div>
                    <strong>Schmidt Orthogonalization:</strong> Build orthonormal
                    polynomials <Math>{'P_0, P_1, \\ldots, P_{n_r}'}</Math> via Gram-Schmidt
                    using the moment inner product{' '}
                    <Math>{'\\langle f, g \\rangle = \\sum_k \\mu_k \\cdot (fg)_k'}</Math>.
                    Starting from <Math>{'P_0 = 1/\\sqrt{\\mu_0}'}</Math>, each successive
                    polynomial is orthogonalized against all previous ones and normalized.
                  </div>
                </li>
                <li className="flex items-start">
                  <span className="flex-shrink-0 w-6 h-6 bg-green-100 text-green-700 rounded-full flex items-center justify-center text-xs font-semibold mr-3">
                    3
                  </span>
                  <div>
                    <strong>Companion Matrix Eigenvalues:</strong> Build the companion
                    matrix from the coefficients of <Math>{'P_{n_r}'}</Math> and compute its
                    eigenvalues via QR iteration. The eigenvalues are the quadrature
                    nodes <Math>{'t_i \\in (0, 1)'}</Math>.
                  </div>
                </li>
                <li className="flex items-start">
                  <span className="flex-shrink-0 w-6 h-6 bg-green-100 text-green-700 rounded-full flex items-center justify-center text-xs font-semibold mr-3">
                    4
                  </span>
                  <div>
                    <strong>Weight Computation:</strong> For each root{' '}
                    <Math>t_i</Math>, the weight is{' '}
                    <Math>{'w_i = 1 \\big/ \\sum_{j=0}^{n_r - 1} P_j(t_i)^2'}</Math>.
                    This formula follows from the Christoffel-Darboux identity for
                    orthonormal polynomials.
                  </div>
                </li>
              </ol>
            </div>
            <div className="bg-green-50 border border-green-200 rounded-lg p-4 mt-3">
              <p className="text-green-800 text-sm">
                <strong>Alternative approach:</strong> The Golub-Welsch algorithm (1969) computes
                quadrature rules via Cholesky factorization of Hankel matrices and Jacobi matrix
                eigendecomposition. While mathematically equivalent, the Schmidt orthogonalization
                approach used in IQCP (and libcint) tends to be more numerically stable for the
                Rys weight function.
              </p>
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
                    Must lie in the interval <Math>{'[0, 1)'}</Math> (root = 0 can occur
                    in degenerate limits, e.g., <Math>{'T = 0'}</Math>)
                  </li>
                  <li>Are eigenvalues of the Jacobi matrix</li>
                  <li>Are roots of orthogonal polynomials</li>
                  <li>Are distinct and well-separated for non-degenerate cases</li>
                </ul>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">
                  Weights (<Math>w_i</Math>)
                </h4>
                <ul className="text-sm text-slate-600 space-y-1">
                  <li>Must be non-negative (strictly positive for non-degenerate cases)</li>
                  <li>
                    Sum to <Math>{'\\mu_0 = F_0(T)'}</Math>
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

          <Subsection title="Practical Limits and Conditioning">
            <p className="text-slate-600 mb-3">
              The Hankel matrix condition number <Math>{'\\kappa(H)'}</Math> grows exponentially
              with the quadrature order <Math>n</Math>. Because the moments{' '}
              <Math>{'\\mu_k = F_k(T)'}</Math> decay rapidly with <Math>k</Math>, the Hankel
              matrix becomes increasingly ill-conditioned as its dimension grows. In IEEE 754
              double precision (approximately 15-16 significant digits), the practical limit is:
            </p>
            <MathBlock label="Practical Double-Precision Limit">
              {'n_r \\leq 5\\text{-}6 \\quad (\\kappa(H) \\sim 10^{12} \\text{ at } n_r = 6)'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              Beyond this limit, the Cholesky factorization of <Math>{'\\mathbf{H}'}</Math> loses
              too many significant digits, producing roots and weights that fail the moment
              reconstruction test. For higher angular momentum integrals that would require{' '}
              <Math>{'n_r > 6'}</Math>, production codes combine multiple lower-order quadratures
              rather than attempting a single high-order rule.
            </p>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
              <p className="text-blue-800 text-sm">
                <strong>IQCP implementation:</strong> The Rys Quadrature Lab displays the Hankel
                matrix condition number <Math>{'\\kappa(H)'}</Math> in the internals panel,
                allowing students to observe the exponential growth first-hand. Try increasing{' '}
                <Math>n</Math> from 2 to 6 and watch the condition number climb by roughly two
                orders of magnitude per additional root.
              </p>
            </div>
          </Subsection>

          <Subsection title="Moment Reconstruction Verification">
            <p className="text-slate-600 mb-3">
              A correctly computed set of Rys roots and weights must exactly reproduce the moments
              of the weight function. This provides a self-consistency check:
            </p>
            <MathBlock label="Moment Reconstruction Condition">
              {'\\sum_{i=1}^{n_r} w_i \\, t_i^{2k} = F_k(T), \\quad k = 0, 1, \\ldots, 2n_r - 1'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              For well-conditioned cases (<Math>{'n_r \\leq 5'}</Math>), the reconstruction error
              for each moment should be at or near machine precision (approximately{' '}
              <Math>{'10^{-14}'}</Math> to <Math>{'10^{-15}'}</Math>). Larger errors indicate
              that the roots and weights have been corrupted by floating-point roundoff during
              the Cholesky or eigendecomposition steps.
            </p>
            <p className="text-slate-600 mb-3">
              IQCP verifies this reconstruction automatically and displays the per-moment error
              alongside the computed roots and weights. Perfect reconstruction to machine
              precision confirms that the quadrature rule is correct, serving as a built-in
              consistency check that students can inspect without needing external reference data.
            </p>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
              <p className="text-amber-800 text-sm">
                <strong>Diagnostic tip:</strong> When the reconstruction error for the highest
                moments exceeds <Math>{'10^{-10}'}</Math>, the quadrature order is approaching
                the conditioning limit. This is normal and expected -- it is not an implementation
                bug but a fundamental limitation of the moment-based algorithm in finite precision.
              </p>
            </div>
          </Subsection>
        </CollapsibleSection>

        {/* Module E: SCF/DIIS + DFT */}
        <CollapsibleSection
          title="Module E: SCF Sandbox"
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
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 mt-3">
              <p className="text-blue-800 text-sm">
                <strong>Initial guess:</strong> IQCP uses multiple initial guess strategies.
                For DFT calculations, a Superposition of Atomic Densities (SAD) guess is
                preferred. When SAD is unavailable (e.g., spherical harmonics), the MINAO
                (Minimal AO) projection guess is used, matching PySCF&apos;s default. The core
                Hamiltonian diagonalization serves as the final fallback. For PES scans, the
                converged density from the previous geometry point seeds the next point,
                with automatic fallback to MINAO if the seeded SCF takes too many iterations.
              </p>
            </div>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 mt-3">
              <p className="text-blue-800 text-sm">
                <strong>Try it:</strong> In <Link to={ROUTES.SCF} className="text-purple-700 underline hover:text-purple-900">Module E</Link>, run an RHF calculation with DIIS disabled
                to see how many iterations are needed. Then enable DIIS and compare -- you should
                see dramatically faster convergence.
              </p>
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
                      <Math>{'10^{-6}'}</Math> Ha
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-4}'}</Math>
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Medium</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-8}'}</Math> Ha
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-6}'}</Math>
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Tight</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-10}'}</Math> Ha
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      <Math>{'10^{-8}'}</Math>
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
                  <li>IQCP uses an 8-vector rolling window by default, with DIIS starting from iteration 1</li>
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
                  near-metallic systems. Available in Module E&apos;s Options tab as a dropdown
                  (Off / Mild 0.25 / Standard 0.5 / Strong 1.0).
                </p>
              </div>
            </div>
          </Subsection>

          <Subsection title="Kohn-Sham DFT">
            <p className="text-slate-600 mb-3">
              Kohn-Sham density functional theory replaces the exact exchange with an
              exchange-correlation (XC) functional of the electron density. The KS total energy is:
            </p>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4 mb-3">
              <p className="text-amber-800 text-sm">
                <strong>Note:</strong> DFT is not an approximation to Hartree-Fock. It is a
                fundamentally different theoretical framework based on the Hohenberg-Kohn
                theorems, which establish that the ground-state energy is a functional of the
                electron density.
              </p>
            </div>
            <MathBlock label="KS Total Energy">
              {'E_{\\text{KS}} = E_{\\text{1e}} + E_J + E_{\\text{XC}}[\\rho] + E_{\\text{nuc}}'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              IQCP supports three XC functional levels:
            </p>
            <div className="overflow-x-auto">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Functional
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Components
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Exact Exchange
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">LDA</td>
                    <td className="px-4 py-3 text-sm text-slate-600">Slater exchange + VWN5 correlation</td>
                    <td className="px-4 py-3 text-sm text-slate-600">0%</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">B3LYP</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      Three-parameter hybrid: <Math>{'(1\\!-\\!a_0)'}</Math>Slater +{' '}
                      <Math>{'a_0'}</Math>HF + <Math>{'a_x'}</Math>Becke88 +{' '}
                      <Math>{'(1\\!-\\!a_c)'}</Math>VWN5 + <Math>{'a_c'}</Math>LYP
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      20% (<Math>{'a_0 = 0.2'}</Math>)
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">B3LYP-D3(BJ)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">B3LYP + Grimme D3 dispersion (BJ damping)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      20% (<Math>{'a_0 = 0.2'}</Math>)
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
          </Subsection>

          <Subsection title="KS Fock Matrix">
            <p className="text-slate-600 mb-3">
              The KS Fock matrix differs from the HF Fock matrix by incorporating the XC potential
              and scaling exact exchange:
            </p>
            <MathBlock label="Kohn-Sham Fock Matrix">
              {'F^{\\text{KS}}_{\\mu\\nu} = H^{\\text{core}}_{\\mu\\nu} + J_{\\mu\\nu} - \\frac{a_0}{2} K_{\\mu\\nu} + V^{\\text{XC}}_{\\mu\\nu}'}
            </MathBlock>
            <p className="text-slate-600 mt-3">
              where <Math>{'a_0'}</Math> is the fraction of exact exchange (0 for pure DFT
              functionals like LDA, 0.2 for B3LYP), and <Math>{'V^{\\text{XC}}'}</Math> is the
              XC potential matrix obtained by numerical integration of the XC functional
              derivatives on a grid.
            </p>
          </Subsection>

          <Subsection title="Numerical Integration Grid">
            <p className="text-slate-600 mb-3">
              The XC energy and potential require numerical integration over 3D space. IQCP uses
              atom-centered grids with two components:
            </p>
            <div className="grid md:grid-cols-2 gap-4">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Radial Quadrature</h4>
                <p className="text-sm text-slate-600">
                  Mura-Knowles log-based radial grid. Maps <Math>{'[0, \\infty)'}</Math> to a
                  finite set of radial points with appropriate Jacobian weighting. Concentrates
                  points near the nucleus where the density varies rapidly.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Angular Quadrature</h4>
                <p className="text-sm text-slate-600">
                  Lebedev grids provide angular integration points on the unit sphere. These grids
                  are exact for spherical harmonics up to a given order and are highly efficient
                  due to octahedral symmetry.
                </p>
              </div>
            </div>
            <p className="text-slate-600 mt-3">
              Becke partitioning assigns each grid point a weight proportional to its proximity
              to the parent atom, ensuring smooth decomposition of the molecular integral into
              atomic contributions.
            </p>
          </Subsection>

          <Subsection title="Geometry Optimization">
            <p className="text-slate-600 mb-3">
              IQCP supports geometry optimization with automatic selection between two strategies
              based on the system size:
            </p>
            <div className="grid md:grid-cols-2 gap-4 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Polyatomic (3+ atoms)</h4>
                <p className="text-sm text-slate-600">
                  Trust-radius Newton-Raphson in redundant internal coordinates (bonds, angles,
                  dihedrals) with the Wilson B-matrix for Cartesian-to-internal transformation,
                  BFGS Hessian update with Wolfe curvature condition, Schlegel model Hessian
                  initial guess, and eigenvalue-floored Newton step. Adaptive trust radius
                  grows/shrinks based on quality ratio (actual vs predicted energy change).
                  Typically converges in 3-5 steps for molecules like H₂O.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Diatomic (2 atoms)</h4>
                <p className="text-sm text-slate-600">
                  L-BFGS quasi-Newton method in Cartesian coordinates. Used as a fallback
                  when internal coordinates are not beneficial (only one bond stretch degree
                  of freedom).
                </p>
              </div>
            </div>
            <p className="text-slate-600 mb-3">
              Both strategies use analytical nuclear gradients:
            </p>
            <MathBlock label="Nuclear Gradient (4-term decomposition)">
              {'\\frac{\\partial E}{\\partial R_A} = \\sum_{\\mu\\nu} D_{\\mu\\nu} \\frac{\\partial H^{\\text{core}}_{\\mu\\nu}}{\\partial R_A} + \\frac{1}{2} \\sum D_{\\mu\\nu} D_{\\lambda\\sigma} \\frac{\\partial}{\\partial R_A}\\!\\left[(\\mu\\nu|\\lambda\\sigma) - \\tfrac{1}{2}(\\mu\\lambda|\\nu\\sigma)\\right] - \\sum_{\\mu\\nu} W_{\\mu\\nu} \\frac{\\partial S_{\\mu\\nu}}{\\partial R_A} + \\frac{\\partial V_{\\text{nn}}}{\\partial R_A}'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              The four terms are: (1) one-electron integral derivatives, (2) two-electron Coulomb
              and exchange integral derivatives, (3) the Pulay force arising from basis function
              dependence on nuclear positions, and (4) the classical nuclear repulsion gradient.
              Here <Math>{'D_{\\mu\\nu}'}</Math> is the density matrix and{' '}
              <Math>{'W_{\\mu\\nu} = 2 \\sum_i^{\\text{occ}} \\varepsilon_i C_{\\mu i} C_{\\nu i}'}</Math>{' '}
              is the energy-weighted density matrix.
            </p>
            <p className="text-slate-600 mb-3">
              Convergence criteria for geometry optimization:
            </p>
            <div className="overflow-x-auto">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Criterion
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Threshold
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Max gradient component</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'4.5 \\times 10^{-4}'}</Math> Ha/Bohr
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">RMS gradient</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'3.0 \\times 10^{-4}'}</Math> Ha/Bohr
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Max displacement</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'1.8 \\times 10^{-3}'}</Math> Bohr
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">RMS displacement</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'1.2 \\times 10^{-3}'}</Math> Bohr
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 mt-3">
              <p className="text-blue-800 text-sm">
                <strong>IQCP implementation:</strong> IQCP checks two convergence criteria:
                maximum gradient component ({`<`} 4.5 × 10⁻⁴ Ha/bohr) and energy change
                ({`<`} 1.0 × 10⁻⁶ Ha). Both must be satisfied simultaneously. The RMS gradient
                and displacement criteria shown above are standard in geomeTRIC but are not
                used as convergence gates in IQCP. The max gradient is the primary indicator,
                matching the approach used by many production codes.
              </p>
            </div>
          </Subsection>

          <Subsection title="Population Analysis">
            <p className="text-slate-600 mb-3">
              Population analysis assigns partial atomic charges by partitioning the total electron
              density among atoms. IQCP supports two schemes:
            </p>
            <div className="space-y-3">
              <MathBlock label="Mulliken Charge">
                {'q_A = Z_A - \\sum_{\\mu \\in A} (\\mathbf{P} \\mathbf{S})_{\\mu\\mu}'}
              </MathBlock>
              <MathBlock label="Lowdin Charge">
                {'q_A = Z_A - \\sum_{\\mu \\in A} (\\mathbf{S}^{1/2} \\mathbf{P} \\mathbf{S}^{1/2})_{\\mu\\mu}'}
              </MathBlock>
            </div>
            <p className="text-slate-600 mt-3">
              Mulliken analysis equally splits overlap populations between atoms, making it
              basis-set dependent. Lowdin analysis uses symmetrically orthogonalized orbitals,
              which reduces (but does not eliminate) basis-set sensitivity. Both schemes sum to
              the total number of electrons:{' '}
              <Math>{'\\sum_A q_A = \\sum_A Z_A - N_{\\text{elec}} = Q_{\\text{total}}'}</Math>.
            </p>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4 mt-3">
              <p className="text-amber-800 text-sm">
                <strong>Warning:</strong> Mulliken charges can become unphysical (negative
                populations) with diffuse basis sets like 6-31+G*. Lowdin charges are generally
                more stable but still should be interpreted qualitatively rather than quantitatively.
              </p>
            </div>
          </Subsection>

          <Subsection title="Potential Energy Surface Scans">
            <p className="text-slate-600 mb-3">
              A potential energy surface (PES) scan maps the total electronic energy as a function
              of a geometric parameter. IQCP supports scanning along three types of internal
              coordinates:
            </p>
            <ul className="text-slate-600 text-sm space-y-1 mb-3 ml-4 list-disc">
              <li><strong>Bond length</strong> — distance between two atoms (bohr)</li>
              <li><strong>Bond angle</strong> — angle between three atoms (degrees)</li>
              <li><strong>Dihedral angle</strong> — torsion angle between four atoms (degrees)</li>
            </ul>
            <p className="text-slate-600 mb-3">
              Two scan modes are available:
            </p>
            <div className="grid md:grid-cols-2 gap-4 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Rigid Scan</h4>
                <p className="text-sm text-slate-600">
                  Only the scanned coordinate changes; all other atomic positions are frozen.
                  Fast but ignores geometric relaxation effects.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Relaxed Scan</h4>
                <p className="text-sm text-slate-600">
                  At each scan point, the remaining degrees of freedom are optimized while the
                  scanned coordinate is held fixed. More physically meaningful but computationally
                  more expensive.
                </p>
              </div>
            </div>
            <p className="text-slate-600 mb-3">
              To improve convergence, the converged density matrix from each geometry point
              seeds the initial guess for the next point. This density extrapolation exploits
              the continuity of electronic structure along the scan coordinate, reducing the
              number of SCF iterations required at each point.
            </p>
            <MathBlock label="Parabolic Equilibrium Interpolation">
              {'r_e \\approx \\frac{1}{2} \\frac{(r_2^2 - r_3^2) E_1 + (r_3^2 - r_1^2) E_2 + (r_1^2 - r_2^2) E_3}{(r_2 - r_3) E_1 + (r_3 - r_1) E_2 + (r_1 - r_2) E_3}'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              The equilibrium geometry is estimated by parabolic interpolation through the
              three lowest-energy grid points. This provides a sub-grid estimate without
              requiring additional SCF calculations.
            </p>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
              <p className="text-blue-800 text-sm">
                <strong>IQCP implementation:</strong> The PES Scan tab plots energy vs coordinate
                value as points accumulate, with geometry animation via a slider synchronized
                with the 3D viewer. A coordinate tracking panel shows how all internal coordinates
                (bonds, angles, dihedrals) respond to the scanned parameter. Students can observe
                how the PES shape changes with basis set (deeper well with larger basis) and
                method (RHF vs DFT), and compare rigid vs relaxed curves to see the effect of
                geometric relaxation on the energy profile.
              </p>
            </div>
          </Subsection>

          <Subsection title="3D Molecular Visualization">
            <p className="text-slate-600 mb-3">
              IQCP renders molecules as interactive 3D structures using Three.js and WebGL,
              providing spatial context that complements the numerical matrix data. The
              visualization employs a sphere-and-stick model:
            </p>
            <ul className="text-slate-600 text-sm space-y-2 mb-3 ml-4 list-disc">
              <li>
                <strong>Atoms</strong> are drawn as spheres with element-specific CPK colors
                (H white, C gray, N blue, O red, S yellow, F/Cl green) and radii proportional
                to covalent radii.
              </li>
              <li>
                <strong>Bonds</strong> are detected automatically using a distance criterion:
                two atoms are bonded when their separation is less than{' '}
                <Math>{'1.2 \\times (r_{\\text{cov}}(A) + r_{\\text{cov}}(B))'}</Math>, where{' '}
                <Math>{'r_{\\text{cov}}'}</Math> is the covalent radius.
              </li>
              <li>
                <strong>Interaction:</strong> left-drag rotates, scroll zooms, right-drag pans,
                with a reset button to restore the default viewpoint.
              </li>
            </ul>
            <p className="text-slate-600 mb-3">
              The visualization updates in real-time as geometry changes during optimization,
              allowing students to watch the molecule relax toward its equilibrium structure.
            </p>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
              <p className="text-amber-800 text-sm">
                <strong>Note:</strong> The 3D viewer requires WebGL support in the browser. If
                WebGL is unavailable, the viewer displays a fallback message while all other
                Module E features (matrices, plots, convergence data) remain fully functional.
              </p>
            </div>
          </Subsection>

          <Subsection title="Orbital Isosurface Visualization">
            <p className="text-slate-600 mb-3">
              After SCF convergence, molecular orbital wavefunctions can be visualized as 3D
              isosurfaces. The orbital value at any point in space is computed from the MO
              coefficients and basis functions:
            </p>
            <MathBlock label="Molecular Orbital Wavefunction">
              {'\\psi_i(\\mathbf{r}) = \\sum_{\\mu} C_{\\mu i} \\, \\chi_\\mu(\\mathbf{r})'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              IQCP evaluates <Math>{'\\psi_i(\\mathbf{r})'}</Math> on a uniform 3D grid
              encompassing the molecule, then extracts isosurfaces at a user-selected threshold
              using the marching cubes algorithm (Lorensen & Cline, 1987). Each grid cube is
              classified by the sign pattern of the orbital at its eight corners, and the
              appropriate triangle mesh is generated from a lookup table.
            </p>
            <div className="grid md:grid-cols-2 gap-4 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Positive Lobe</h4>
                <p className="text-sm text-slate-600">
                  Rendered as a solid blue mesh at the isovalue{' '}
                  <Math>{'+\\sigma'}</Math>. Indicates regions where the orbital
                  wavefunction is positive.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Negative Lobe</h4>
                <p className="text-sm text-slate-600">
                  Rendered as a translucent red mesh at the isovalue{' '}
                  <Math>{'-\\sigma'}</Math>. Indicates regions where the orbital
                  wavefunction is negative.
                </p>
              </div>
            </div>
            <p className="text-slate-600 mb-3">
              Students can switch between occupied and virtual orbitals using the orbital
              selector, with orbital energies <Math>{'\\varepsilon_i'}</Math> displayed
              alongside the visualization. The default isovalue of 0.03 a.u. can be adjusted
              via a slider (range 0.005 to 0.10) to explore how the orbital shape changes
              with the chosen contour level.
            </p>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
              <p className="text-blue-800 text-sm">
                <strong>IQCP implementation:</strong> The marching cubes computation runs in a
                Web Worker to keep the UI responsive. Positive and negative lobes use distinct
                rendering styles (solid vs translucent) rather than relying solely on color,
                ensuring accessibility for color-vision-deficient users.
              </p>
            </div>
          </Subsection>

          <Subsection title="Electron Density and Deformation Density">
            <p className="text-slate-600 mb-3">
              The total electron density is constructed from the converged density matrix and
              basis functions:
            </p>
            <MathBlock label="Total Electron Density">
              {'\\rho(\\mathbf{r}) = \\sum_{\\mu\\nu} P_{\\mu\\nu} \\, \\chi_\\mu(\\mathbf{r}) \\, \\chi_\\nu(\\mathbf{r})'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              While the total density shows where electrons are located in the molecule, the
              deformation density reveals how bonding redistributes charge relative to the
              non-interacting atoms:
            </p>
            <MathBlock label="Deformation Density">
              {'\\Delta\\rho(\\mathbf{r}) = \\rho_{\\text{mol}}(\\mathbf{r}) - \\rho_{\\text{promol}}(\\mathbf{r})'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              The promolecular density <Math>{'\\rho_{\\text{promol}}'}</Math> is the
              superposition of spherically averaged atomic densities placed at the nuclear
              positions. It represents the electron density of non-interacting atoms at the
              molecular geometry.
            </p>
            <div className="grid md:grid-cols-2 gap-4 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-blue-700 mb-2">
                  Positive <Math>{'\\Delta\\rho'}</Math> (blue)
                </h4>
                <p className="text-sm text-slate-600">
                  Charge accumulation -- regions where bonding concentrates electron density
                  beyond what isolated atoms would produce. Typically found between bonded atoms
                  and in lone pair regions.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-red-700 mb-2">
                  Negative <Math>{'\\Delta\\rho'}</Math> (red)
                </h4>
                <p className="text-sm text-slate-600">
                  Charge depletion -- regions where bonding removes electron density. Typically
                  found along the axis connecting bonded atoms, just outside the bonding region.
                </p>
              </div>
            </div>
            <p className="text-slate-600 mb-3">
              Cross-section plots show the density along user-selected planes through the
              molecule, providing a 2D view of the 3D density distribution. This is
              particularly useful for examining bonding patterns in planar or near-planar
              molecules.
            </p>
          </Subsection>

          <Subsection title="DFT Method Comparison">
            <p className="text-slate-600 mb-3">
              IQCP can run both RHF and a DFT method on the same geometry for direct comparison,
              highlighting how the treatment of electron-electron interactions differs between
              the two approaches. The energy decomposition makes this transparent:
            </p>
            <div className="overflow-x-auto mb-3">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Component
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      RHF
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      DFT (e.g., B3LYP)
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">One-electron energy</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'E_{\\text{1e}}'}</Math>
                    </td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'E_{\\text{1e}}'}</Math>
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Coulomb repulsion</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'E_J'}</Math>
                    </td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'E_J'}</Math>
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Exchange</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'E_K'}</Math> (exact, 100%)
                    </td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'a_0 E_K'}</Math> (fraction, e.g., 20%)
                    </td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">XC functional</td>
                    <td className="px-4 py-3 text-sm text-slate-600">--</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">
                      <Math>{'E_{\\text{XC}}[\\rho]'}</Math>
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
            <p className="text-slate-600 mb-3">
              Students can observe that DFT typically converges in fewer iterations than RHF
              for the same system, and that the inclusion of correlation via the XC functional
              generally lowers the total energy. The energy decomposition table makes clear
              exactly where HF and DFT differ: HF uses 100% exact exchange and no correlation,
              while B3LYP mixes 20% exact exchange with DFT exchange and adds correlation.
            </p>
          </Subsection>

          <Subsection title="Basis Set Comparison in SCF">
            <p className="text-slate-600 mb-3">
              IQCP can run the same molecule and method across multiple basis sets simultaneously,
              producing a side-by-side comparison that illustrates basis set convergence. Results
              are displayed in a comparison table showing:
            </p>
            <ul className="text-slate-600 text-sm space-y-2 mb-3 ml-4 list-disc">
              <li>
                <strong>Total energy</strong> in Hartree -- demonstrating the variational
                principle (larger basis sets always give lower or equal energies)
              </li>
              <li>
                <strong>Basis function count</strong> -- showing the computational cost of
                each basis set
              </li>
              <li>
                <strong>SCF iteration count</strong> -- revealing how basis set size affects
                convergence behavior
              </li>
              <li>
                <strong>Energy difference</strong> from the largest basis -- quantifying how
                far each result is from the basis set limit
              </li>
            </ul>
            <p className="text-slate-600 mb-3">
              A convergence chart plots the total energy against the number of basis functions,
              visually demonstrating how energy decreases as the basis grows. The variational
              principle guarantees this monotonic behavior:
            </p>
            <MathBlock label="Variational Principle">
              {'E(\\text{small basis}) \\geq E(\\text{large basis}) \\geq E_{\\text{HF limit}}'}
            </MathBlock>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 mt-3">
              <p className="text-blue-800 text-sm">
                <strong>IQCP implementation:</strong> The basis set comparison panel runs all
                selected basis sets in sequence and displays results as they complete. The
                convergence chart updates progressively, letting students watch the energy
                decrease as larger basis sets finish. This provides an intuitive understanding
                of diminishing returns -- the energy improvement from STO-3G to 6-31G* is much
                larger than from 6-31G* to cc-pVDZ, even though both steps add similar numbers
                of basis functions.
              </p>
            </div>
          </Subsection>
        </CollapsibleSection>

        {/* Frequency Analysis & Vibrational Spectroscopy */}
        <CollapsibleSection
          title="Frequency Analysis & Vibrational Spectroscopy"
          id="frequency"
          expandedSection={expandedSection}
          onToggle={handleToggle}
          color="rose"
        >
          <p className="text-slate-600 mb-6">
            Vibrational frequency analysis extends the SCF calculation to predict how molecules
            vibrate, absorb infrared light, scatter Raman photons, and behave thermochemically.
            Starting from the converged electronic wavefunction, IQCP computes the analytical
            Hessian (matrix of second derivatives of the energy with respect to nuclear
            coordinates), transforms it into normal mode frequencies, and derives all
            spectroscopic and thermodynamic observables. This entire pipeline runs in the
            browser via WebAssembly, accessible from Module E&apos;s Frequency tab.
          </p>

          <Subsection title="Analytical Hessian (Second Derivatives)">
            <p className="text-slate-600 mb-3">
              The Hessian matrix contains the second derivatives of the total energy with respect
              to nuclear coordinates. Each element describes how the force on one atom changes as
              another atom is displaced:
            </p>
            <MathBlock label="Hessian Element">
              {'H_{A_d, B_e} = \\frac{\\partial^2 E}{\\partial R_{A,d} \\, \\partial R_{B,e}}'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              For a molecule with <Math>N</Math> atoms, the Hessian is a{' '}
              <Math>{'3N \\times 3N'}</Math> symmetric matrix. IQCP computes it analytically
              (not by finite differences), which is more accurate and requires only one CPHF
              solve rather than <Math>{'6N'}</Math> displaced SCF calculations.
            </p>
            <p className="text-slate-600 mb-3">
              The analytical Hessian decomposes into four contributions:
            </p>
            <MathBlock label="Hessian Decomposition">
              {'\\mathbf{H} = \\mathbf{H}^{\\text{nuc}} + \\mathbf{H}^{\\text{core}} + \\mathbf{H}^{\\text{eri}} + \\mathbf{H}^{\\text{cphf}}'}
            </MathBlock>
            <div className="space-y-2 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-3">
                <p className="text-sm text-slate-700">
                  <strong><Math>{'\\mathbf{H}^{\\text{nuc}}'}</Math>:</strong> Nuclear repulsion
                  second derivatives -- the purely classical Coulomb contribution from nuclei.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-3">
                <p className="text-sm text-slate-700">
                  <strong><Math>{'\\mathbf{H}^{\\text{core}}'}</Math>:</strong> Second derivatives
                  of one-electron integrals (kinetic energy and nuclear attraction), contracted with
                  the density matrix.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-3">
                <p className="text-sm text-slate-700">
                  <strong><Math>{'\\mathbf{H}^{\\text{eri}}'}</Math>:</strong> Second derivatives
                  of two-electron repulsion integrals (ERIs), computed via Rys quadrature at
                  angular momentum <Math>{'L + 4'}</Math> to accommodate the two additional
                  derivative operators.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-3">
                <p className="text-sm text-slate-700">
                  <strong><Math>{'\\mathbf{H}^{\\text{cphf}}'}</Math>:</strong> The response
                  contribution from the Coupled-Perturbed Hartree-Fock equations. This accounts
                  for how the density matrix relaxes when nuclei move -- the most computationally
                  demanding piece.
                </p>
              </div>
            </div>
            <p className="text-slate-600 mb-3">
              For DFT methods (LDA, B3LYP), an additional exchange-correlation (XC) second
              derivative term <Math>{'\\mathbf{H}^{\\text{xc}}'}</Math> is included, built from
              the second functional derivative <Math>{'f_{\\text{xc}}'}</Math> evaluated on the
              numerical integration grid.
            </p>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
              <p className="text-blue-800 text-sm">
                <strong>IQCP validation:</strong> For H₂O/STO-3G/RHF, the IQCP Hessian matches
                PySCF to a maximum absolute difference of{' '}
                <Math>{'1.26 \\times 10^{-7}'}</Math> Ha/Bohr². For B3LYP, agreement reaches{' '}
                <Math>{'8.68 \\times 10^{-11}'}</Math> Ha/Bohr². These are analytical, not
                numerical, second derivatives.
              </p>
            </div>
          </Subsection>

          <Subsection title="Coupled-Perturbed Hartree-Fock (CPHF)">
            <p className="text-slate-600 mb-3">
              When a nucleus is displaced, the electronic wavefunction must relax in response.
              The CPHF equations determine this first-order response -- how the molecular orbital
              coefficients change with nuclear position. This is essential because the Hessian
              contains terms that depend on {' '}
              <Math>{'\\partial \\mathbf{C} / \\partial R_A'}</Math>.
            </p>
            <MathBlock label="CPHF Equations">
              {'(\\varepsilon_a - \\varepsilon_i) U^{(1)}_{ai} + \\sum_{bj} A_{ai,bj} U^{(1)}_{bj} = -h^{(1)}_{ai} + \\varepsilon_i s^{(1)}_{ai}'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              Here <Math>{'U^{(1)}_{ai}'}</Math> is the mixing of virtual orbital{' '}
              <Math>a</Math> into occupied orbital <Math>i</Math> due to the perturbation,{' '}
              <Math>{'A_{ai,bj}'}</Math> is the two-electron coupling kernel,{' '}
              <Math>{'h^{(1)}'}</Math> and <Math>{'s^{(1)}'}</Math> are the
              derivatives of the Fock and overlap matrices. The orbital energy denominator{' '}
              <Math>{'(\\varepsilon_a - \\varepsilon_i)'}</Math> appears naturally from the
              zeroth-order eigenvalue equation.
            </p>
            <p className="text-slate-600 mb-3">
              IQCP solves the CPHF equations using a Krylov subspace method that avoids
              constructing the full <Math>{'A'}</Math> tensor. The right-hand side is computed
              from overlap and Fock matrix derivatives (the &quot;skeleton&quot; derivatives),
              and the solver iterates until the residual falls below{' '}
              <Math>{'10^{-9}'}</Math>.
            </p>
            <p className="text-slate-600 mb-3">
              IQCP implements two CPHF paths:
            </p>
            <div className="grid md:grid-cols-2 gap-4 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Nuclear Perturbation</h4>
                <p className="text-sm text-slate-600">
                  Used for the Hessian and IR intensities. The basis functions depend on nuclear
                  positions, so <Math>{'\\mathbf{S}^{(1)} \\neq 0'}</Math>. The
                  occupied-occupied block is fixed at{' '}
                  <Math>{'-\\mathbf{S}^{(1)}/2'}</Math>.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Electric Field Perturbation</h4>
                <p className="text-sm text-slate-600">
                  Used for Raman activities. The basis functions are field-independent,
                  so <Math>{'\\mathbf{S}^{(1)} = 0'}</Math>. Only the
                  virtual-occupied block is solved.
                </p>
              </div>
            </div>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
              <p className="text-amber-800 text-sm">
                <strong>Key insight:</strong> The CPHF equations are a linear system, not an
                eigenvalue problem. They need to be solved only once for all perturbations
                simultaneously. The solution{' '}
                <Math>{'U^{(1)}'}</Math> is reused for the Hessian, IR intensities, and
                (in the field variant) Raman activities -- making the overall cost dominated
                by the CPHF solve rather than the property evaluation.
              </p>
            </div>
          </Subsection>

          <Subsection title="Normal Mode Analysis">
            <p className="text-slate-600 mb-3">
              Normal modes are the natural vibrational patterns of a molecule -- collective
              motions where all atoms oscillate at the same frequency. To find them, the Hessian
              must first be mass-weighted, then purified of translational and rotational degrees
              of freedom.
            </p>
            <MathBlock label="Mass-Weighted Hessian">
              {'H^{\\text{mw}}_{A_d, B_e} = \\frac{H_{A_d, B_e}}{\\sqrt{m_A \\cdot m_B}}'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              A nonlinear molecule with <Math>N</Math> atoms has{' '}
              <Math>{'3N - 6'}</Math> vibrational degrees of freedom (or{' '}
              <Math>{'3N - 5'}</Math> for a linear molecule). The remaining 6 (or 5) degrees
              correspond to translation and rotation, which must be projected out before
              diagonalization.
            </p>
            <p className="text-slate-600 mb-3">
              IQCP performs the projection in several steps:
            </p>
            <div className="bg-white rounded-lg border border-slate-200 p-4 mb-3">
              <ol className="space-y-2 text-sm text-slate-700">
                <li>
                  <strong>1. Center of mass:</strong> Translate to the center of mass frame.
                </li>
                <li>
                  <strong>2. Principal axes:</strong> Diagonalize the moment-of-inertia tensor to
                  find the principal axis frame and classify the rotor type (atom, linear,
                  spherical top, symmetric top, or asymmetric top).
                </li>
                <li>
                  <strong>3. Build TR subspace:</strong> Construct 6 (or 5) translation and
                  rotation vectors in mass-weighted coordinates.
                </li>
                <li>
                  <strong>4. QR orthogonalization:</strong> Apply QR decomposition to separate the
                  translational/rotational subspace from the vibrational subspace.
                </li>
                <li>
                  <strong>5. Project and diagonalize:</strong> Project the mass-weighted Hessian
                  into the vibrational subspace and diagonalize to obtain eigenvalues{' '}
                  <Math>{'\\lambda_k'}</Math> and eigenvectors (normal modes).
                </li>
              </ol>
            </div>
            <p className="text-slate-600 mb-3">
              The eigenvalues are converted to wavenumbers using the atomic unit conversion factor:
            </p>
            <MathBlock label="Frequency Conversion">
              {'\\tilde{\\nu}_k \\; [\\text{cm}^{-1}] = \\text{sign}(\\lambda_k) \\times \\sqrt{|\\lambda_k|} \\times \\text{AU\\_TO\\_CM1}'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              Negative eigenvalues produce imaginary frequencies (printed as negative cm⁻¹
              values), indicating that the geometry is a saddle point (transition state) rather
              than a true minimum. A true minimum has all real (positive) frequencies.
            </p>
            <p className="text-slate-600 mb-3">
              Each normal mode also yields a reduced mass and force constant:
            </p>
            <div className="grid md:grid-cols-2 gap-4 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Reduced Mass</h4>
                <p className="text-sm text-slate-600">
                  <Math>{'\\mu_k = \\left( \\sum_{A,d} |L_{A_d,k}|^2 / m_A \\right)^{-1}'}</Math>{' '}
                  in amu, where <Math>{'L_{A_d,k}'}</Math> is the normal mode displacement vector.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Force Constant</h4>
                <p className="text-sm text-slate-600">
                  <Math>{'k_k = \\lambda_k \\cdot \\mu_k'}</Math> converted to mDyne/A using
                  the appropriate unit factors. Measures the stiffness of the vibration.
                </p>
              </div>
            </div>
            <div className="overflow-x-auto mb-3">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Molecule
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Modes
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Frequencies (cm⁻¹, STO-3G/RHF)
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">H₂O</td>
                    <td className="px-4 py-3 text-sm text-slate-600">3</td>
                    <td className="px-4 py-3 text-sm text-slate-600">2043, 4488, 4790</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">CO₂</td>
                    <td className="px-4 py-3 text-sm text-slate-600">4</td>
                    <td className="px-4 py-3 text-sm text-slate-600">421, 421, 1571, 2829</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">CH₄</td>
                    <td className="px-4 py-3 text-sm text-slate-600">9</td>
                    <td className="px-4 py-3 text-sm text-slate-600">1688 (3x), 1908 (2x), 3476, 3733 (3x)</td>
                  </tr>
                </tbody>
              </table>
            </div>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4 mb-3">
              <p className="text-amber-800 text-sm">
                <strong>Note:</strong> STO-3G/RHF frequencies are substantially higher than
                experiment because the Hartree-Fock method neglects electron correlation and the
                minimal basis set is too inflexible. For example, the experimental O-H stretch
                in water is approximately 3657 cm⁻¹, while STO-3G/RHF gives 4488 cm⁻¹.
                Empirical scaling factors (typically 0.8-0.9) are commonly applied in practice
                to correct for this systematic overestimation.
              </p>
            </div>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
              <p className="text-blue-800 text-sm">
                <strong>Try it:</strong> In <Link to={ROUTES.SCF} className="text-purple-700 underline hover:text-purple-900">Module E</Link>, converge an SCF calculation, then open the
                Frequency tab and click &quot;Run.&quot; The sortable frequency table shows
                all modes with their wavenumber, symmetry classification, reduced mass, force
                constant, IR intensity, and Raman activity.
              </p>
            </div>
          </Subsection>

          <Subsection title="IR Spectroscopy">
            <p className="text-slate-600 mb-3">
              Infrared spectroscopy measures how strongly a molecule absorbs light at each
              vibrational frequency. The physical origin is the coupling between the oscillating
              electric field of light and the molecular dipole moment. A vibrational mode is
              IR-active if it changes the dipole moment of the molecule.
            </p>
            <MathBlock label="IR Selection Rule">
              {'\\text{Mode } k \\text{ is IR-active if } \\frac{\\partial \\boldsymbol{\\mu}}{\\partial Q_k} \\neq 0'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              To compute IR intensities, IQCP first evaluates the molecular dipole moment and
              its derivatives with respect to nuclear coordinates. The dipole moment is:
            </p>
            <MathBlock label="Dipole Moment">
              {'\\mu_d = \\sum_A Z_A R_{A,d} - \\text{Tr}(\\mathbf{D} \\cdot \\boldsymbol{\\mu}_d)'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              where the first term is the nuclear contribution and the second is the electronic
              contribution via dipole integrals{' '}
              <Math>{'(\\boldsymbol{\\mu}_d)_{\\mu\\nu} = \\langle \\mu | r_d | \\nu \\rangle'}</Math>.
              The dipole integrals are computed using the overlap recurrence relation at
              incremented angular momentum <Math>{'(a_d + 1, b)'}</Math>.
            </p>
            <p className="text-slate-600 mb-3">
              The dipole derivatives with respect to normal coordinates{' '}
              <Math>{'\\partial \\mu_d / \\partial Q_k'}</Math> are obtained analytically by
              transforming the Cartesian dipole derivative matrix (which uses the CPHF density
              response <Math>{'U^{(R)}'}</Math> from the Hessian calculation) into the
              normal mode basis. The IR intensity is then:
            </p>
            <MathBlock label="IR Intensity">
              {'A_k = \\frac{N_A \\pi}{3 c^2} \\left| \\frac{\\partial \\boldsymbol{\\mu}}{\\partial Q_k} \\right|^2 \\quad [\\text{km/mol}]'}
            </MathBlock>
            <div className="overflow-x-auto mb-3">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Molecule
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Mode (cm⁻¹)
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      IR Intensity (km/mol)
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Active?
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">H₂O</td>
                    <td className="px-4 py-3 text-sm text-slate-600">2043 (bend)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">11.4</td>
                    <td className="px-4 py-3 text-sm text-green-600 font-medium">Yes</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">H₂O</td>
                    <td className="px-4 py-3 text-sm text-slate-600">4488 (sym stretch)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">37.4</td>
                    <td className="px-4 py-3 text-sm text-green-600 font-medium">Yes</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">H₂O</td>
                    <td className="px-4 py-3 text-sm text-slate-600">4790 (asym stretch)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">20.0</td>
                    <td className="px-4 py-3 text-sm text-green-600 font-medium">Yes</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">CO₂</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      1571 (sym stretch, <Math>{'\\sigma_g'}</Math>)
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">~0</td>
                    <td className="px-4 py-3 text-sm text-red-600 font-medium">No</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">CO₂</td>
                    <td className="px-4 py-3 text-sm text-slate-600">2829 (asym stretch)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">342.7</td>
                    <td className="px-4 py-3 text-sm text-green-600 font-medium">Yes</td>
                  </tr>
                </tbody>
              </table>
            </div>
            <p className="text-slate-600 mb-3">
              All three modes of H₂O are IR-active because each vibration changes the dipole
              moment of this bent molecule. In contrast, the symmetric stretch of CO₂ is
              IR-inactive: the two C=O bonds stretch equally and in phase, so the dipole moment
              (which is zero by symmetry) does not change.
            </p>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
              <p className="text-blue-800 text-sm">
                <strong>IQCP implementation:</strong> The IR tab in the Frequency panel displays
                the simulated IR spectrum as both stick lines (at exact frequencies) and a
                broadened envelope. Click any peak to select the corresponding normal mode and
                see its displacement pattern.
              </p>
            </div>
          </Subsection>

          <Subsection title="Raman Spectroscopy">
            <p className="text-slate-600 mb-3">
              Raman spectroscopy measures inelastic light scattering, where the energy shift
              of scattered photons corresponds to vibrational frequencies. A mode is Raman-active
              if the vibration changes the molecular polarizability -- the ease with which the
              electron cloud is distorted by an external electric field.
            </p>
            <MathBlock label="Raman Selection Rule">
              {'\\text{Mode } k \\text{ is Raman-active if } \\frac{\\partial \\boldsymbol{\\alpha}}{\\partial Q_k} \\neq 0'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              IQCP computes the static polarizability tensor by solving the electric field CPHF
              equations (the <Math>{'\\mathbf{S}^{(1)} = 0'}</Math> variant):
            </p>
            <MathBlock label="Static Polarizability">
              {'\\alpha_{de} = -2 \\sum_{a,i} U^{(E_d)}_{ai} \\, (\\boldsymbol{\\mu}_e)_{ai}'}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              where <Math>{'U^{(E_d)}'}</Math> is the CPHF response to a field in direction{' '}
              <Math>d</Math>. The polarizability derivatives are computed semi-analytically
              by combining the field CPHF response with the nuclear CPHF response from the
              Hessian.
            </p>
            <p className="text-slate-600 mb-3">
              From the polarizability derivative tensor for each normal mode, two scalar
              invariants determine the Raman activity:
            </p>
            <div className="grid md:grid-cols-2 gap-4 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">
                  Mean Polarizability Derivative <Math>{"\\bar{\\alpha}'_k"}</Math>
                </h4>
                <p className="text-sm text-slate-600">
                  <Math>{"\\bar{\\alpha}'_k = \\frac{1}{3} \\text{Tr}(\\partial \\boldsymbol{\\alpha} / \\partial Q_k)"}</Math>{' '}
                  -- measures isotropic change in polarizability.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">
                  Anisotropy Derivative <Math>{"\\gamma'_k"}</Math>
                </h4>
                <p className="text-sm text-slate-600">
                  Measures the directional asymmetry of the polarizability change. Constructed
                  from off-diagonal and trace-removed diagonal elements of the derivative tensor.
                </p>
              </div>
            </div>
            <MathBlock label="Raman Activity">
              {"S_k = 45(\\bar{\\alpha}'_k)^2 + 7(\\gamma'_k)^2 \\quad [\\text{\\AA}^4/\\text{amu}]"}
            </MathBlock>
            <MathBlock label="Depolarization Ratio">
              {"\\rho_k = \\frac{3(\\gamma'_k)^2}{45(\\bar{\\alpha}'_k)^2 + 4(\\gamma'_k)^2}"}
            </MathBlock>
            <p className="text-slate-600 mt-3 mb-3">
              The depolarization ratio <Math>{'\\rho_k'}</Math> ranges from 0 (completely
              polarized, highly symmetric vibration) to 0.75 (depolarized). It provides
              experimental information about the symmetry of each vibrational mode.
            </p>
            <div className="overflow-x-auto mb-3">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Molecule
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Mode (cm⁻¹)
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Raman (A⁴/amu)
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      <Math>{'\\rho'}</Math>
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Active?
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">H₂O</td>
                    <td className="px-4 py-3 text-sm text-slate-600">2043 (bend)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">8.17</td>
                    <td className="px-4 py-3 text-sm text-slate-600">0.717</td>
                    <td className="px-4 py-3 text-sm text-green-600 font-medium">Yes</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">H₂O</td>
                    <td className="px-4 py-3 text-sm text-slate-600">4488 (sym stretch)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">43.81</td>
                    <td className="px-4 py-3 text-sm text-slate-600">0.185</td>
                    <td className="px-4 py-3 text-sm text-green-600 font-medium">Yes</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">H₂O</td>
                    <td className="px-4 py-3 text-sm text-slate-600">4790 (asym stretch)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">17.54</td>
                    <td className="px-4 py-3 text-sm text-slate-600">0.750</td>
                    <td className="px-4 py-3 text-sm text-green-600 font-medium">Yes</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">CO₂</td>
                    <td className="px-4 py-3 text-sm text-slate-600">
                      1571 (sym stretch, <Math>{'\\sigma_g'}</Math>)
                    </td>
                    <td className="px-4 py-3 text-sm text-slate-600">7.57</td>
                    <td className="px-4 py-3 text-sm text-slate-600">0.210</td>
                    <td className="px-4 py-3 text-sm text-green-600 font-medium">Yes</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm font-mono text-slate-700">CO₂</td>
                    <td className="px-4 py-3 text-sm text-slate-600">2829 (asym stretch)</td>
                    <td className="px-4 py-3 text-sm text-slate-600">~0</td>
                    <td className="px-4 py-3 text-sm text-slate-600">--</td>
                    <td className="px-4 py-3 text-sm text-red-600 font-medium">No</td>
                  </tr>
                </tbody>
              </table>
            </div>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4 mb-3">
              <p className="text-amber-800 text-sm">
                <strong>Mutual exclusion rule:</strong> For molecules with a center of inversion
                (like CO₂), vibrations that are IR-active are Raman-inactive, and vice versa.
                The symmetric stretch of CO₂ is Raman-active (polarizability changes) but
                IR-inactive (no dipole change). The asymmetric stretch is IR-active but
                Raman-inactive. IQCP computes both properties so students can verify this
                fundamental symmetry principle numerically.
              </p>
            </div>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
              <p className="text-blue-800 text-sm">
                <strong>IQCP implementation:</strong> The Raman tab in the Frequency panel shows
                the simulated Raman spectrum with activities as peak heights. IQCP is the first
                open-source tool to provide Raman activity calculations entirely in the browser
                via WebAssembly, enabling classroom demonstrations without any server
                infrastructure.
              </p>
            </div>
          </Subsection>

          <Subsection title="RRHO Thermochemistry">
            <p className="text-slate-600 mb-3">
              The Rigid-Rotor Harmonic-Oscillator (RRHO) model uses the computed vibrational
              frequencies and molecular geometry to estimate thermodynamic properties. It
              partitions the total energy and entropy into four independent contributions:
            </p>
            <div className="space-y-2 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-3">
                <h4 className="font-semibold text-slate-700 mb-1">Translational</h4>
                <p className="text-sm text-slate-600">
                  From the Sackur-Tetrode equation. Depends on molecular mass, temperature, and
                  pressure. Contributes <Math>{'C_v = \\frac{3}{2}R'}</Math>.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-3">
                <h4 className="font-semibold text-slate-700 mb-1">Rotational</h4>
                <p className="text-sm text-slate-600">
                  Classical rigid-rotor partition function using the principal moments of inertia.
                  Linear molecules contribute <Math>{'C_v = R'}</Math>, nonlinear contribute{' '}
                  <Math>{'C_v = \\frac{3}{2}R'}</Math>. The symmetry number{' '}
                  <Math>{'\\sigma'}</Math> prevents overcounting equivalent orientations
                  (e.g., <Math>{'\\sigma = 2'}</Math> for H₂O).
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-3">
                <h4 className="font-semibold text-slate-700 mb-1">Vibrational</h4>
                <p className="text-sm text-slate-600">
                  Quantum harmonic oscillator partition function for each normal mode. The
                  zero-point energy (ZPE) is the sum of half-quanta:{' '}
                  <Math>{'\\text{ZPE} = \\frac{1}{2} \\sum_k h \\nu_k'}</Math>. The thermal
                  vibrational contribution depends on the Boltzmann population of excited states.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-3">
                <h4 className="font-semibold text-slate-700 mb-1">Electronic</h4>
                <p className="text-sm text-slate-600">
                  <Math>{'S_{\\text{elec}} = R \\ln(2S + 1)'}</Math> where <Math>S</Math> is
                  the total spin. For closed-shell singlet states (<Math>{'2S + 1 = 1'}</Math>),
                  this contributes zero.
                </p>
              </div>
            </div>
            <p className="text-slate-600 mb-3">
              The Gibbs free energy assembles all contributions:
            </p>
            <MathBlock label="Gibbs Free Energy">
              {'G = E_{\\text{elec}} + \\text{ZPE} + E_{\\text{trans}} + E_{\\text{rot}} + E_{\\text{vib}} + RT - T S_{\\text{tot}}'}
            </MathBlock>
            <div className="overflow-x-auto mb-3">
              <table className="min-w-full bg-white rounded-lg overflow-hidden">
                <thead className="bg-slate-100">
                  <tr>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      Property
                    </th>
                    <th className="px-4 py-3 text-left text-sm font-semibold text-slate-700">
                      H₂O/STO-3G/RHF (298.15 K)
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-200">
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">SCF Energy</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">-74.96302 Ha</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">ZPE</td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">0.02579 Ha</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">
                      <Math>{'E_0'}</Math> (E + ZPE)
                    </td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">-74.93723 Ha</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Enthalpy <Math>H</Math></td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">-74.93345 Ha</td>
                  </tr>
                  <tr>
                    <td className="px-4 py-3 text-sm text-slate-700">Gibbs Free Energy <Math>G</Math></td>
                    <td className="px-4 py-3 text-sm font-mono text-slate-600">-74.95486 Ha</td>
                  </tr>
                </tbody>
              </table>
            </div>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 mb-3">
              <p className="text-blue-800 text-sm">
                <strong>IQCP implementation:</strong> The Thermochemistry panel in the Frequency
                tab provides a temperature slider (100-1000 K) that recomputes all thermodynamic
                quantities client-side as the slider moves. The translational, rotational,
                vibrational, and electronic contributions are displayed separately so students
                can see how each partition function contributes to the total. Only the vibrational
                partition function changes significantly with temperature for typical molecules
                at room temperature.
              </p>
            </div>
            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
              <p className="text-amber-800 text-sm">
                <strong>Limitation:</strong> The RRHO model assumes that vibrations are perfectly
                harmonic and that rotation-vibration coupling is negligible. For floppy molecules
                with low-frequency modes (below about 100 cm⁻¹), the harmonic approximation
                breaks down and more sophisticated treatments such as hindered rotors may be
                needed. IQCP flags low-frequency modes in the frequency table to alert students
                to this limitation.
              </p>
            </div>
          </Subsection>

          <Subsection title="Simulated IR and Raman Spectra">
            <p className="text-slate-600 mb-3">
              Experimental spectra show broadened peaks rather than infinitely sharp lines. IQCP
              simulates this broadening by convolving the stick spectrum (discrete
              frequency/intensity pairs) with a lineshape function. Two lineshapes are available:
            </p>
            <div className="grid md:grid-cols-2 gap-4 mb-3">
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Lorentzian</h4>
                <MathBlock label="">
                  {'L(\\tilde{\\nu}) = \\frac{A_k}{\\pi} \\cdot \\frac{\\gamma}{(\\tilde{\\nu} - \\tilde{\\nu}_k)^2 + \\gamma^2}'}
                </MathBlock>
                <p className="text-sm text-slate-600 mt-2">
                  where <Math>{'\\gamma = \\text{FWHM}/2'}</Math>. Models homogeneous broadening
                  from natural linewidth and pressure broadening. Has longer tails than Gaussian.
                </p>
              </div>
              <div className="bg-white rounded-lg border border-slate-200 p-4">
                <h4 className="font-semibold text-slate-700 mb-2">Gaussian</h4>
                <MathBlock label="">
                  {'G(\\tilde{\\nu}) = \\frac{A_k}{\\sigma\\sqrt{2\\pi}} \\exp\\!\\left(-\\frac{(\\tilde{\\nu} - \\tilde{\\nu}_k)^2}{2\\sigma^2}\\right)'}
                </MathBlock>
                <p className="text-sm text-slate-600 mt-2">
                  where <Math>{'\\sigma = \\text{FWHM} / (2\\sqrt{2 \\ln 2})'}</Math>. Models
                  inhomogeneous (Doppler) broadening. Peaks decay faster than Lorentzian.
                </p>
              </div>
            </div>
            <p className="text-slate-600 mb-3">
              Both kernels integrate to 1 over the full real line, so the integrated area of each
              peak in the broadened spectrum equals the original stick amplitude. The default
              full width at half maximum (FWHM) is 8 cm⁻¹, adjustable from 2 to 20 cm⁻¹ via a
              slider. The wavenumber grid spans 0 to 4500 cm⁻¹ at 1 cm⁻¹ resolution (4501
              points).
            </p>
            <p className="text-slate-600 mb-3">
              When the FWHM or broadening type is changed, the broadened spectrum is recomputed
              client-side in the browser without requiring a new WASM calculation -- only the
              stick frequencies and intensities come from the compute engine.
            </p>
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
              <p className="text-blue-800 text-sm">
                <strong>Try it:</strong> In the Frequency tab, toggle between Lorentzian and
                Gaussian broadening to see how the peak shapes differ. Increase the FWHM to
                simulate lower-resolution experimental conditions, or decrease it to resolve
                closely spaced peaks. Click any peak in the broadened spectrum to highlight the
                corresponding stick line and select that normal mode.
              </p>
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
          <Subsection title="Basis Sets and Gaussian Integrals">
            <ul className="space-y-3 text-sm text-slate-600">
              <li className="border-l-4 border-amber-300 pl-4">
                <strong>Boys, S. F.</strong> (1950). Electronic Wave Functions. I. A General Method
                of Calculation for the Stationary States of Any Molecular System.{' '}
                <em>Proc. R. Soc. London A</em>, 200, 542-554.
                <p className="text-slate-500 text-xs mt-1">
                  Introduction of Gaussian-type orbitals in quantum chemistry.
                </p>
              </li>
              <li className="border-l-4 border-amber-300 pl-4">
                <strong>Hehre, W. J., Stewart, R. F., & Pople, J. A.</strong> (1969). Self-consistent
                molecular-orbital methods. I. Use of Gaussian expansions of Slater-type atomic
                orbitals. <em>J. Chem. Phys.</em>, 51, 2657-2664.
                <p className="text-slate-500 text-xs mt-1">
                  Foundation of the STO-3G basis set and Pople-style basis set design.
                </p>
              </li>
              <li className="border-l-4 border-amber-300 pl-4">
                <strong>Ditchfield, R., Hehre, W. J., & Pople, J. A.</strong> (1971). Self-consistent
                molecular-orbital methods. IX. An extended Gaussian-type basis for molecular-orbital
                studies of organic molecules. <em>J. Chem. Phys.</em>, 54, 724-728.
                <p className="text-slate-500 text-xs mt-1">
                  Introduction of the 6-31G split-valence basis set.
                </p>
              </li>
              <li className="border-l-4 border-amber-300 pl-4">
                <strong>Hariharan, P. C., & Pople, J. A.</strong> (1973). The influence of polarization
                functions on molecular orbital hydrogenation energies. <em>Theor. Chim. Acta</em>, 28, 213-222.
                <p className="text-slate-500 text-xs mt-1">
                  Addition of d-type polarization functions to split-valence basis sets (6-31G*).
                </p>
              </li>
              <li className="border-l-4 border-amber-300 pl-4">
                <strong>Clark, T., Chandrasekhar, J., Spitznagel, G. W., & Schleyer, P. v. R.</strong> (1983).
                Efficient diffuse function-augmented basis sets for anion calculations. III. The
                3-21+G basis set for first-row elements, Li-F. <em>J. Comput. Chem.</em>, 4, 294-301.
                <p className="text-slate-500 text-xs mt-1">
                  Diffuse function augmentation for Pople basis sets (6-31+G*).
                </p>
              </li>
              <li className="border-l-4 border-amber-300 pl-4">
                <strong>Dunning, T. H. Jr.</strong> (1989). Gaussian basis sets for use in correlated
                molecular calculations. I. The atoms boron through neon and hydrogen.{' '}
                <em>J. Chem. Phys.</em>, 90, 1007-1023.
                <p className="text-slate-500 text-xs mt-1">
                  Correlation-consistent basis sets (cc-pVDZ, cc-pVTZ, cc-pVQZ).
                </p>
              </li>
            </ul>
          </Subsection>

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

          <Subsection title="SCF, DIIS, and DFT">
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
              <li className="border-l-4 border-purple-300 pl-4">
                <strong>Hohenberg, P., & Kohn, W.</strong> (1964). Inhomogeneous electron gas.{' '}
                <em>Phys. Rev.</em>, 136, B864-B871.
                <p className="text-slate-500 text-xs mt-1">
                  First Hohenberg-Kohn theorem — foundation of density functional theory.
                </p>
              </li>
              <li className="border-l-4 border-purple-300 pl-4">
                <strong>Kohn, W., & Sham, L. J.</strong> (1965). Self-consistent equations including
                exchange and correlation effects. <em>Phys. Rev.</em>, 140, A1133-A1138.
                <p className="text-slate-500 text-xs mt-1">
                  Kohn-Sham equations — practical DFT via auxiliary non-interacting system.
                </p>
              </li>
              <li className="border-l-4 border-purple-300 pl-4">
                <strong>Vosko, S. H., Wilk, L., & Nusair, M.</strong> (1980). Accurate spin-dependent
                electron liquid correlation energies for local spin density calculations: a critical
                analysis. <em>Can. J. Phys.</em>, 58, 1200-1211.
                <p className="text-slate-500 text-xs mt-1">
                  VWN correlation functional used in LDA (Slater + VWN5).
                </p>
              </li>
              <li className="border-l-4 border-purple-300 pl-4">
                <strong>Lee, C., Yang, W., & Parr, R. G.</strong> (1988). Development of the Colle-Salvetti
                correlation-energy formula into a functional of the electron density.{' '}
                <em>Phys. Rev. B</em>, 37, 785-789.
                <p className="text-slate-500 text-xs mt-1">
                  LYP correlation functional — the correlation component of B3LYP.
                </p>
              </li>
              <li className="border-l-4 border-purple-300 pl-4">
                <strong>Becke, A. D.</strong> (1993). Density-functional thermochemistry. III. The
                role of exact exchange. <em>J. Chem. Phys.</em>, 98, 5648-5652.
                <p className="text-slate-500 text-xs mt-1">
                  Three-parameter hybrid functional mixing exact and DFT exchange.
                </p>
              </li>
              <li className="border-l-4 border-purple-300 pl-4">
                <strong>Stephens, P. J., Devlin, F. J., Chabalowski, C. F., & Frisch, M. J.</strong> (1994).
                Ab initio calculation of vibrational absorption and circular dichroism spectra using
                density functional force fields. <em>J. Phys. Chem.</em>, 98, 11623-11627.
                <p className="text-slate-500 text-xs mt-1">
                  Defines the standard B3LYP parameterization used in modern codes.
                </p>
              </li>
              <li className="border-l-4 border-purple-300 pl-4">
                <strong>Grimme, S., Ehrlich, S., & Goerigk, L.</strong> (2011). Effect of the damping
                function in dispersion corrected density functional theory. <em>J. Comput. Chem.</em>,
                32, 1456-1465.
                <p className="text-slate-500 text-xs mt-1">
                  D3 dispersion correction with Becke-Johnson damping.
                </p>
              </li>
            </ul>
          </Subsection>

          <Subsection title="Population Analysis">
            <ul className="space-y-3 text-sm text-slate-600">
              <li className="border-l-4 border-teal-300 pl-4">
                <strong>Mulliken, R. S.</strong> (1955). Electronic population analysis on LCAO-MO
                molecular wave functions. I. <em>J. Chem. Phys.</em>, 23, 1833-1840.
                <p className="text-slate-500 text-xs mt-1">
                  Original formulation of Mulliken population analysis.
                </p>
              </li>
            </ul>
          </Subsection>

          <Subsection title="Vibrational Analysis, Spectroscopy, and Thermochemistry">
            <ul className="space-y-3 text-sm text-slate-600">
              <li className="border-l-4 border-rose-300 pl-4">
                <strong>Wilson, E. B., Decius, J. C., & Cross, P. C.</strong> (1955).{' '}
                <em>Molecular Vibrations: The Theory of Infrared and Raman Vibrational Spectra</em>.
                McGraw-Hill (reprinted Dover, 1980).
                <p className="text-slate-500 text-xs mt-1">
                  Foundational treatment of GF matrix method, normal mode analysis, and symmetry-based
                  selection rules for vibrational spectroscopy.
                </p>
              </li>
              <li className="border-l-4 border-rose-300 pl-4">
                <strong>Pople, J. A., Krishnan, R., Schlegel, H. B., & Binkley, J. S.</strong> (1979).
                Derivative studies in Hartree-Fock and Moller-Plesset theories.{' '}
                <em>Int. J. Quantum Chem. Symp.</em>, 13, 225-241.
                <p className="text-slate-500 text-xs mt-1">
                  CPHF theory and analytical second derivatives for Hartree-Fock wavefunctions.
                </p>
              </li>
              <li className="border-l-4 border-rose-300 pl-4">
                <strong>Amos, R. D.</strong> (1986). Calculation of polarizability derivatives using
                analytic gradient methods. <em>Chem. Phys. Lett.</em>, 124, 376-381.
                <p className="text-slate-500 text-xs mt-1">
                  Semi-analytical computation of polarizability derivatives for Raman activities.
                </p>
              </li>
              <li className="border-l-4 border-rose-300 pl-4">
                <strong>Herzberg, G.</strong> (1945).{' '}
                <em>Molecular Spectra and Molecular Structure. II. Infrared and Raman Spectra of
                Polyatomic Molecules</em>. Van Nostrand.
                <p className="text-slate-500 text-xs mt-1">
                  Classic reference for vibrational spectroscopy of polyatomic molecules.
                </p>
              </li>
              <li className="border-l-4 border-rose-300 pl-4">
                <strong>Long, D. A.</strong> (2002).{' '}
                <em>The Raman Effect: A Unified Treatment of the Theory of Raman Scattering by
                Molecules</em>. Wiley.
                <p className="text-slate-500 text-xs mt-1">
                  Comprehensive treatment of Raman scattering theory, polarizability tensors, and
                  depolarization ratios.
                </p>
              </li>
              <li className="border-l-4 border-rose-300 pl-4">
                <strong>McQuarrie, D. A.</strong> (1976).{' '}
                <em>Statistical Mechanics</em>. Harper & Row.
                <p className="text-slate-500 text-xs mt-1">
                  Standard reference for statistical thermodynamics, partition functions, and the
                  RRHO model used in thermochemistry calculations.
                </p>
              </li>
            </ul>
          </Subsection>

          <Subsection title="Visualization and Numerical Methods">
            <ul className="space-y-3 text-sm text-slate-600">
              <li className="border-l-4 border-slate-300 pl-4">
                <strong>Lorensen, W. E., & Cline, H. E.</strong> (1987). Marching cubes: A high
                resolution 3D surface construction algorithm. <em>ACM SIGGRAPH Computer
                Graphics</em>, 21(4), 163-169.
                <p className="text-slate-500 text-xs mt-1">
                  Isosurface extraction algorithm used for orbital and density visualization.
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
    </div>
  );
}

export default Reference;
