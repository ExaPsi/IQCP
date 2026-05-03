/**
 * About page - Comprehensive information about IQCP.
 *
 * Contains project overview, module descriptions, supported capabilities,
 * citation information, technical details, and version/license info.
 *
 * @module pages/About
 */

const GITHUB_URL = 'https://github.com/ExaPsi/IQCP';
const ZENODO_DOI = '10.5281/zenodo.18798310';
const IQCP_URL = 'https://iqcp.dev';

function About() {
  return (
    <div className="max-w-4xl mx-auto">
      <h1 className="text-3xl font-bold text-slate-800 mb-4">About IQCP</h1>
      <p className="text-slate-600 mb-6">
        Interactive Quantum Chemistry Playground — a browser-based Hartree-Fock
        and density functional theory engine for education and research, requiring
        no installation or server.
      </p>

      <div className="space-y-6">
        {/* Overview Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Overview</h2>
          <div className="space-y-4 text-slate-600">
            <p>
              IQCP is a validated computational chemistry implementation that runs
              entirely in the browser using WebAssembly compiled from Rust. It supports
              RHF, LDA (SVWN5), B3LYP, and B3LYP-D3(BJ) methods with six basis sets
              for elements H through Ar (Z = 1-18).
            </p>
            <p>
              Designed for both education and research, IQCP enables students and
              researchers to explore quantum chemistry algorithms through interactive
              visualization while providing sub-microhartree accuracy against PySCF
              for Hartree-Fock calculations.
            </p>
            <div className="grid md:grid-cols-2 gap-4 mt-4">
              <div className="bg-blue-50 rounded-lg p-4">
                <h3 className="font-semibold text-blue-800 text-sm mb-2">For Education</h3>
                <ul className="text-sm text-blue-700 space-y-1">
                  <li>Interactive exploration of Boys functions and Rys quadrature</li>
                  <li>Visualize basis sets, integrals, and SCF convergence</li>
                  <li>3D molecular viewer with orbital isosurfaces</li>
                  <li>Deep links for shareable, reproducible states</li>
                </ul>
              </div>
              <div className="bg-green-50 rounded-lg p-4">
                <h3 className="font-semibold text-green-800 text-sm mb-2">For Research</h3>
                <ul className="text-sm text-green-700 space-y-1">
                  <li>Validated HF and DFT with cross-code verification</li>
                  <li>Analytical gradients and geometry optimization</li>
                  <li>Vibrational frequencies, IR/Raman spectra, thermochemistry</li>
                  <li>PES scanning (rigid and relaxed) along internal coordinates</li>
                  <li>POV-Ray export for publication-quality figures</li>
                </ul>
              </div>
            </div>
          </div>
        </section>

        {/* Modules Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Modules</h2>
          <div className="space-y-4">
            <div className="border-l-4 border-indigo-400 pl-4">
              <h3 className="font-semibold text-slate-800">
                Module A: Basis Explorer
              </h3>
              <p className="text-slate-600 text-sm mt-1">
                Explore 6 Gaussian basis sets (STO-3G through cc-pVDZ) for all 18
                supported elements. Interactive radial profile plots show how contracted
                functions are built from primitives. Compare basis sets side-by-side,
                adjust exponents with real-time feedback, and visualize overlap integrals
                as a function of distance.
              </p>
            </div>
            <div className="border-l-4 border-teal-400 pl-4">
              <h3 className="font-semibold text-slate-800">
                Module B: Integral Inspector
              </h3>
              <p className="text-slate-600 text-sm mt-1">
                Visualize one-electron integral matrices (S, T, V, H<sup>core</sup>) as
                interactive heatmaps for 17 preset molecules. Click any matrix element
                to see its primitive-pair decomposition. Trace the Fock matrix build
                step-by-step (H<sup>core</sup> + J - 0.5K) and browse individual
                electron repulsion integrals.
              </p>
            </div>
            <div className="border-l-4 border-blue-400 pl-4">
              <h3 className="font-semibold text-slate-800">
                Module C: Boys Function Lab
              </h3>
              <p className="text-slate-600 text-sm mt-1">
                Explore the Boys function F<sub>m</sub>(T) with interactive regime
                visualization showing series expansion vs. erf-recurrence boundaries.
                Sweep across T values, compare multiple orders simultaneously, and
                inspect computational details including error estimation.
              </p>
            </div>
            <div className="border-l-4 border-violet-400 pl-4">
              <h3 className="font-semibold text-slate-800">
                Module D: Rys Quadrature Lab
              </h3>
              <p className="text-slate-600 text-sm mt-1">
                Compute Rys polynomial roots and weights for quadrature orders 1-10.
                Visualize order-error tradeoffs, inspect the Hankel/Jacobi matrix
                construction, and verify moment reconstruction against Boys function
                values.
              </p>
            </div>
            <div className="border-l-4 border-amber-400 pl-4">
              <h3 className="font-semibold text-slate-800">
                Module E: SCF Sandbox
              </h3>
              <p className="text-slate-600 text-sm mt-1">
                The primary computational module with five workflow tabs:
              </p>
              <ul className="text-slate-600 text-sm mt-2 space-y-1 ml-4 list-disc">
                <li>
                  <strong>Single Point:</strong> Run RHF or KS-DFT (LDA, B3LYP, B3LYP-D3BJ)
                  calculations. Inspect convergence plots, matrices (S, H, F, P), orbital
                  energies, 3D orbital/density isosurfaces, and Mulliken/Lowdin charges.
                </li>
                <li>
                  <strong>Optimize:</strong> L-BFGS geometry optimization with trust-radius
                  control in redundant internal coordinates. Trajectory animation with ghost
                  overlay of initial geometry.
                </li>
                <li>
                  <strong>Frequency:</strong> Analytical Hessian with CPHF for vibrational
                  frequencies, IR intensities, Raman activities, normal mode animation,
                  and RRHO thermochemistry with temperature/pressure controls.
                </li>
                <li>
                  <strong>PES Scan:</strong> Potential energy surface scanning along bonds,
                  angles, or dihedrals in rigid or relaxed mode. Side-by-side 3D viewer
                  with geometry animation and coordinate tracking.
                </li>
                <li>
                  <strong>Compare:</strong> Side-by-side RHF vs DFT comparison and multi-basis
                  set benchmarking.
                </li>
              </ul>
            </div>
          </div>
        </section>

        {/* Capabilities Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Capabilities</h2>
          <div className="grid md:grid-cols-3 gap-4 text-sm">
            <div>
              <h3 className="font-semibold text-slate-700 mb-2">Methods</h3>
              <ul className="text-slate-600 space-y-1">
                <li>Restricted Hartree-Fock (RHF)</li>
                <li>LDA (Slater + VWN5)</li>
                <li>B3LYP (Becke88 + LYP)</li>
                <li>B3LYP-D3(BJ) dispersion</li>
                <li>DIIS convergence acceleration</li>
                <li>Analytical gradients and Hessian</li>
                <li>Vibrational frequencies (IR/Raman)</li>
                <li>RRHO thermochemistry</li>
              </ul>
            </div>
            <div>
              <h3 className="font-semibold text-slate-700 mb-2">Basis Sets</h3>
              <ul className="text-slate-600 space-y-1">
                <li>STO-3G (minimal)</li>
                <li>3-21G (split-valence)</li>
                <li>6-31G (split-valence)</li>
                <li>6-31G* (polarized)</li>
                <li>6-31+G* (augmented)</li>
                <li>cc-pVDZ (correlation-consistent)</li>
              </ul>
            </div>
            <div>
              <h3 className="font-semibold text-slate-700 mb-2">Elements</h3>
              <ul className="text-slate-600 space-y-1">
                <li>H-He (period 1)</li>
                <li>Li-Ne (period 2)</li>
                <li>Na-Ar (period 3)</li>
                <li>Cartesian and spherical harmonics</li>
                <li>17 preset molecules</li>
                <li>Custom geometry input</li>
              </ul>
            </div>
          </div>
        </section>

        {/* Validation Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Validation</h2>
          <div className="text-slate-600 space-y-3">
            <p>
              IQCP is cross-validated against PySCF 2.11.0 across 120 benchmark
              combinations (6 molecules, 6 basis sets, 4 methods):
            </p>
            <div className="grid md:grid-cols-2 gap-4">
              <div className="bg-slate-50 rounded-lg p-3">
                <h4 className="font-semibold text-slate-700 text-sm">RHF Accuracy</h4>
                <p className="text-sm mt-1">
                  Sub-microhartree agreement ({`<`}10<sup>-6</sup> Ha) for total energies
                  across all supported molecules and basis sets.
                </p>
              </div>
              <div className="bg-slate-50 rounded-lg p-3">
                <h4 className="font-semibold text-slate-700 text-sm">DFT Accuracy</h4>
                <p className="text-sm mt-1">
                  Sub-nanoHartree agreement on the same numerical grid; systematic
                  grid-dependent deviations of ~0.03-0.05 mHa from PySCF defaults.
                </p>
              </div>
              <div className="bg-slate-50 rounded-lg p-3">
                <h4 className="font-semibold text-slate-700 text-sm">Gradients</h4>
                <p className="text-sm mt-1">
                  Analytical gradients validated against PySCF for RHF, LDA, and
                  B3LYP across multiple molecules and basis sets.
                </p>
              </div>
              <div className="bg-slate-50 rounded-lg p-3">
                <h4 className="font-semibold text-slate-700 text-sm">Frequencies</h4>
                <p className="text-sm mt-1">
                  Analytical Hessian with CPHF validated against PySCF for
                  vibrational frequencies, IR intensities, and Raman activities.
                </p>
              </div>
            </div>
          </div>
        </section>

        {/* Citation Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Citation</h2>
          <p className="text-slate-600 mb-4">
            If you use IQCP in your teaching or research, please cite:
          </p>
          <div className="bg-slate-50 p-4 rounded-lg text-sm text-slate-700 space-y-3">
            <div>
              <p className="font-medium text-slate-800 mb-1">Software:</p>
              <p>
                Sawatlon, B.; Paiboonvorachat, N.; Vchirawongkwin, V.
                IQCP: Interactive Quantum Chemistry Playground. 2026.{' '}
                <a
                  href={`https://doi.org/${ZENODO_DOI}`}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-blue-600 hover:underline"
                >
                  DOI: {ZENODO_DOI}
                </a>
              </p>
            </div>
          </div>
        </section>

        {/* Technical Details Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">
            Technical Details
          </h2>
          <div className="grid md:grid-cols-2 gap-4 text-sm">
            <div>
              <h3 className="font-semibold text-slate-700 mb-2">Frontend</h3>
              <ul className="text-slate-600 space-y-1">
                <li>React 18 with TypeScript (strict mode)</li>
                <li>Vite build system</li>
                <li>TailwindCSS styling</li>
                <li>Plotly.js interactive plots</li>
                <li>Three.js / React Three Fiber (3D viewer)</li>
                <li>Zustand state management</li>
              </ul>
            </div>
            <div>
              <h3 className="font-semibold text-slate-700 mb-2">Compute Core</h3>
              <ul className="text-slate-600 space-y-1">
                <li>Rust stable with wasm-bindgen</li>
                <li>Web Worker for non-blocking computation</li>
                <li>Deterministic floating-point operations</li>
                <li>~470 KB gzipped core WASM + ~300 KB spectra (lazy-loaded)</li>
                <li>Becke integration grid (Mura-Knowles + Lebedev)</li>
                <li>Schwarz screening for ERI evaluation</li>
              </ul>
            </div>
          </div>
          <div className="mt-4 text-sm text-slate-500">
            <p>
              Architecture: Browser SPA {'->'} Web Worker {'->'} WASM (qc-core) {'->'} Results.
              All computation runs client-side; no server required.
            </p>
          </div>
        </section>

        {/* Export Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Export</h2>
          <div className="text-slate-600 text-sm space-y-2">
            <p>IQCP supports multiple export formats for publication and sharing:</p>
            <ul className="list-disc list-inside space-y-1 ml-2">
              <li><strong>Deep Links:</strong> Every module state is URL-encodable for sharing reproducible results</li>
              <li><strong>SVG:</strong> All Plotly plots exportable as scalable vector graphics</li>
              <li><strong>PNG:</strong> 3D molecular viewer with 2x supersampled export</li>
              <li><strong>POV-Ray:</strong> Complete .pov scene files with orbital isosurfaces for ray-traced rendering</li>
              <li><strong>Clipboard:</strong> Matrix data in TSV, JSON, or LaTeX format</li>
            </ul>
          </div>
        </section>

        {/* Links Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Links</h2>
          <div className="flex flex-wrap gap-3">
            <a
              href={IQCP_URL}
              target="_blank"
              rel="noopener noreferrer"
              className="inline-flex items-center gap-2 px-4 py-2 rounded-lg bg-blue-50 text-blue-700 hover:bg-blue-100 text-sm font-medium transition-colors"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 12a9 9 0 01-9 9m9-9a9 9 0 00-9-9m9 9H3m9 9a9 9 0 01-9-9m9 9c1.657 0 3-4.03 3-9s-1.343-9-3-9m0 18c-1.657 0-3-4.03-3-9s1.343-9 3-9m-9 9a9 9 0 019-9" />
              </svg>
              iqcp.dev
            </a>
            <a
              href={GITHUB_URL}
              target="_blank"
              rel="noopener noreferrer"
              className="inline-flex items-center gap-2 px-4 py-2 rounded-lg bg-slate-100 text-slate-700 hover:bg-slate-200 text-sm font-medium transition-colors"
            >
              <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 24 24">
                <path d="M12 0C5.37 0 0 5.37 0 12c0 5.31 3.435 9.795 8.205 11.385.6.105.825-.255.825-.57 0-.285-.015-1.23-.015-2.235-3.015.555-3.795-.735-4.035-1.41-.135-.345-.72-1.41-1.23-1.695-.42-.225-1.02-.78-.015-.795.945-.015 1.62.87 1.845 1.23 1.08 1.815 2.805 1.305 3.495.99.105-.78.42-1.305.765-1.605-2.67-.3-5.46-1.335-5.46-5.925 0-1.305.465-2.385 1.23-3.225-.12-.3-.54-1.53.12-3.18 0 0 1.005-.315 3.3 1.23.96-.27 1.98-.405 3-.405s2.04.135 3 .405c2.295-1.56 3.3-1.23 3.3-1.23.66 1.65.24 2.88.12 3.18.765.84 1.23 1.905 1.23 3.225 0 4.605-2.805 5.625-5.475 5.925.435.375.81 1.095.81 2.22 0 1.605-.015 2.895-.015 3.3 0 .315.225.69.825.57A12.02 12.02 0 0024 12c0-6.63-5.37-12-12-12z" />
              </svg>
              GitHub
            </a>
            <a
              href={`https://doi.org/${ZENODO_DOI}`}
              target="_blank"
              rel="noopener noreferrer"
              className="inline-flex items-center gap-2 px-4 py-2 rounded-lg bg-green-50 text-green-700 hover:bg-green-100 text-sm font-medium transition-colors"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
              </svg>
              Zenodo DOI
            </a>
          </div>
        </section>

        {/* Version Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Version</h2>
          <div className="flex items-center space-x-4">
            <div className="inline-flex items-center px-3 py-1 rounded-full bg-primary-100 text-primary-800 text-sm font-medium">
              v1.0.0-study
            </div>
            <span className="text-slate-500 text-sm">Study Release</span>
          </div>
          <p className="text-slate-600 mt-3 text-sm">
            This version accompanies the submitted manuscripts (M1: J. Chem. Educ.,
            M2: J. Chem. Inf. Model.) and includes all 5 educational modules with
            validated RHF and DFT implementations.
          </p>
        </section>

        {/* Authors Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Authors</h2>
          <div className="text-slate-600 text-sm space-y-2">
            <p>
              <strong>Boodsarin Sawatlon</strong>,{' '}
              <strong>Nattapong Paiboonvorachat</strong>,{' '}
              <strong>Viwat Vchirawongkwin</strong>
            </p>
            <p>
              Department of Chemistry, Faculty of Science, Chulalongkorn University,
              Bangkok, Thailand
            </p>
            <p>
              Contact:{' '}
              <a href="mailto:viwat.v@chula.ac.th" className="text-blue-600 hover:underline">
                viwat.v@chula.ac.th
              </a>
            </p>
          </div>
        </section>

        {/* License Section */}
        <section className="bg-slate-50 rounded-xl p-6 border border-slate-200">
          <h2 className="text-lg font-semibold text-slate-800 mb-2">License</h2>
          <p className="text-slate-600 text-sm">
            IQCP is open source software released under the{' '}
            <strong>MIT License</strong>. See the{' '}
            <a
              href={GITHUB_URL}
              target="_blank"
              rel="noopener noreferrer"
              className="text-blue-600 hover:underline"
            >
              GitHub repository
            </a>
            {' '}for full license text and contribution guidelines.
          </p>
        </section>
      </div>
    </div>
  );
}

export default About;
