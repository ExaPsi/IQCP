/**
 * About page - Information about IQCP project.
 * Contains overview, citation information, and version details.
 */
function About() {
  return (
    <div className="max-w-4xl mx-auto">
      <h1 className="text-3xl font-bold text-slate-800 mb-4">About IQCP</h1>
      <p className="text-slate-600 mb-6">
        Interactive Quantum Chemistry Playground - An educational web
        application for teaching quantum chemistry concepts through transparent,
        reproducible computation.
      </p>

      <div className="space-y-6">
        {/* Overview Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Overview</h2>
          <div className="space-y-4 text-slate-600">
            <p>
              IQCP (Interactive Quantum Chemistry Playground) is designed for
              undergraduate and graduate chemistry education. It allows students
              to explore fundamental quantum chemistry algorithms through
              interactive visualizations, making abstract concepts tangible.
            </p>
            <p>
              All heavy computation runs in-browser using WebAssembly (WASM),
              compiled from Rust. This architecture ensures:
            </p>
            <ul className="list-disc list-inside space-y-1 ml-4">
              <li>Responsive interactions with latency under 200ms</li>
              <li>
                Cross-browser reproducibility through deterministic algorithms
              </li>
              <li>No server-side computation required</li>
              <li>
                Complete transparency - students can inspect intermediate
                results
              </li>
            </ul>
          </div>
        </section>

        {/* Modules Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Modules</h2>
          <div className="space-y-4">
            <div className="border-l-4 border-primary-500 pl-4">
              <h3 className="font-semibold text-slate-800">
                Module A: Boys Function Lab
              </h3>
              <p className="text-slate-600 text-sm mt-1">
                Explore the Boys function F_m(T) with interactive regime
                visualization. Understand how series expansion, erf-recurrence,
                and asymptotic methods each provide accurate evaluation in their
                respective domains.
              </p>
            </div>
            <div className="border-l-4 border-primary-500 pl-4">
              <h3 className="font-semibold text-slate-800">
                Module B: Rys Quadrature Lab
              </h3>
              <p className="text-slate-600 text-sm mt-1">
                Compute Rys polynomial roots and weights. Visualize order-error
                tradeoffs and understand how quadrature accuracy relates to
                integral evaluation.
              </p>
            </div>
            <div className="border-l-4 border-primary-500 pl-4">
              <h3 className="font-semibold text-slate-800">
                Module C: SCF Sandbox
              </h3>
              <p className="text-slate-600 text-sm mt-1">
                Run restricted Hartree-Fock calculations with optional DIIS
                acceleration. Inspect matrices (S, H_core, F, D, C) and watch
                energy converge iteration by iteration.
              </p>
            </div>
          </div>
        </section>

        {/* Citation Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Citation</h2>
          <p className="text-slate-600 mb-4">
            If you use IQCP in your teaching or research, please cite:
          </p>
          <pre className="bg-slate-100 p-4 rounded-lg text-sm text-slate-700 overflow-x-auto whitespace-pre-wrap">
            {`IQCP: Interactive Quantum Chemistry Playground.
https://iqcp.dev

Citation details will be updated upon publication in J. Chem. Educ.`}
          </pre>
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
                <li>React 18 with TypeScript</li>
                <li>Vite build system</li>
                <li>TailwindCSS styling</li>
                <li>Plotly.js visualizations</li>
              </ul>
            </div>
            <div>
              <h3 className="font-semibold text-slate-700 mb-2">Compute Core</h3>
              <ul className="text-slate-600 space-y-1">
                <li>Rust stable with wasm-bindgen</li>
                <li>Web Worker for non-blocking computation</li>
                <li>Deterministic floating-point operations</li>
                <li>Sub-500KB gzipped WASM bundle</li>
              </ul>
            </div>
          </div>
        </section>

        {/* Version Section */}
        <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h2 className="text-xl font-semibold text-slate-800 mb-3">Version</h2>
          <div className="flex items-center space-x-4">
            <div className="inline-flex items-center px-3 py-1 rounded-full bg-primary-100 text-primary-800 text-sm font-medium">
              v0.1.0
            </div>
            <span className="text-slate-500 text-sm">Development Release</span>
          </div>
          <p className="text-slate-600 mt-3 text-sm">
            IQCP is under active development. Features and APIs may change
            before the 1.0 release accompanying the J. Chem. Educ. publication.
          </p>
        </section>

        {/* License Section */}
        <section className="bg-slate-50 rounded-xl p-6 border border-slate-200">
          <h2 className="text-lg font-semibold text-slate-800 mb-2">License</h2>
          <p className="text-slate-600 text-sm">
            IQCP is open source software. See the GitHub repository for license
            details and contribution guidelines.
          </p>
        </section>
      </div>
    </div>
  );
}

export default About;
