/**
 * Home page — Landing page with module cards and navigation.
 *
 * Displays all 5 educational modules with color-coded cards,
 * plus links to Reference, Labs, and About pages.
 *
 * @module pages/Home
 */

import { Link } from 'react-router-dom';
import { ROUTES } from '@/lib/routes';

interface ModuleCardProps {
  to: string;
  title: string;
  subtitle: string;
  description: string;
  /** Accent color class for the left border */
  accentColor: string;
}

function ModuleCard({ to, title, subtitle, description, accentColor }: ModuleCardProps) {
  return (
    <Link
      to={to}
      className={`block bg-white rounded-xl shadow-sm border border-slate-200 p-6 hover:shadow-md hover:border-primary-300 transition-all focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500 focus-visible:ring-offset-2 border-l-4 ${accentColor}`}
    >
      <h3 className="text-lg font-semibold text-slate-800">{title}</h3>
      <p className="text-primary-600 font-medium mt-1">{subtitle}</p>
      <p className="text-slate-600 mt-3 text-sm">{description}</p>
    </Link>
  );
}

function Home() {
  return (
    <div className="max-w-5xl mx-auto">
      {/* Hero */}
      <div className="text-center mb-12">
        <h1 className="text-4xl font-bold text-slate-800 mb-4">
          Interactive Quantum Chemistry Playground
        </h1>
        <p className="text-lg text-slate-600 max-w-2xl mx-auto">
          Explore basis sets, molecular integrals, Boys functions, Rys quadrature,
          Hartree-Fock/DFT calculations, and vibrational spectroscopy — all running
          in your browser with sub-microhartree accuracy. No installation required.
        </p>
        <div className="flex flex-wrap justify-center gap-3 mt-6">
          <span className="inline-flex items-center px-3 py-1 rounded-full bg-blue-50 text-blue-700 text-xs font-medium">
            RHF / LDA / B3LYP / B3LYP-D3(BJ)
          </span>
          <span className="inline-flex items-center px-3 py-1 rounded-full bg-green-50 text-green-700 text-xs font-medium">
            6 Basis Sets (STO-3G — cc-pVDZ)
          </span>
          <span className="inline-flex items-center px-3 py-1 rounded-full bg-purple-50 text-purple-700 text-xs font-medium">
            Elements H — Ar
          </span>
          <span className="inline-flex items-center px-3 py-1 rounded-full bg-rose-50 text-rose-700 text-xs font-medium">
            IR / Raman Spectra
          </span>
          <span className="inline-flex items-center px-3 py-1 rounded-full bg-amber-50 text-amber-700 text-xs font-medium">
            100% Client-Side (WebAssembly)
          </span>
        </div>
      </div>

      {/* Module E — Featured primary module */}
      <ModuleCard
        to={ROUTES.SCF}
        title="Module E"
        subtitle="SCF Sandbox"
        accentColor="border-l-amber-400"
        description="The primary computational module. Single Point: run RHF or KS-DFT with convergence plots, matrix inspection, 3D orbital/density viewer, and population analysis. Optimize: L-BFGS geometry optimization with trajectory animation. Frequency: vibrational frequencies, IR/Raman spectra with broadening, normal mode animation, and RRHO thermochemistry. PES Scan: rigid and relaxed scanning along bonds, angles, and dihedrals. Compare: side-by-side RHF vs DFT and multi-basis benchmarking. Supports 17 preset molecules, custom geometry input, and POV-Ray export."
      />

      {/* Modules A-D — 2x2 grid */}
      <div className="grid md:grid-cols-2 gap-5 mt-5">
        <ModuleCard
          to={ROUTES.BASIS}
          title="Module A"
          subtitle="Basis Explorer"
          accentColor="border-l-indigo-400"
          description="Explore 6 Gaussian basis sets for 18 elements (H-Ar). Interactive radial profiles, multi-basis comparison, exponent adjustment with real-time feedback, and overlap distance visualization."
        />
        <ModuleCard
          to={ROUTES.INTEGRALS}
          title="Module B"
          subtitle="Integral Inspector"
          accentColor="border-l-teal-400"
          description="Visualize S, T, V, and Hcore matrices as interactive heatmaps for 17 molecules. Primitive-pair decomposition, step-by-step Fock matrix tracing, and ERI quartet browsing."
        />
        <ModuleCard
          to={ROUTES.BOYS}
          title="Module C"
          subtitle="Boys Function Lab"
          accentColor="border-l-blue-400"
          description="Explore the Boys function Fm(T) with regime visualization (series vs recurrence). Sweep plots, multi-order comparison, error estimation, and computational internals."
        />
        <ModuleCard
          to={ROUTES.RYS}
          title="Module D"
          subtitle="Rys Quadrature Lab"
          accentColor="border-l-violet-400"
          description="Compute Rys polynomial roots and weights (n=1-10). Error curves with recommended order, Hankel/Jacobi matrix internals, and moment reconstruction verification."
        />
      </div>

      {/* Secondary Navigation */}
      <div className="mt-8 grid md:grid-cols-3 gap-4">
        <Link
          to={ROUTES.REFERENCE}
          className="flex items-center gap-3 bg-white rounded-xl shadow-sm border border-slate-200 p-4 hover:shadow-md hover:border-primary-300 transition-all focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500"
        >
          <div className="w-10 h-10 rounded-lg bg-slate-100 flex items-center justify-center flex-shrink-0">
            <svg className="w-5 h-5 text-slate-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253" />
            </svg>
          </div>
          <div>
            <h3 className="text-sm font-semibold text-slate-800">Reference</h3>
            <p className="text-xs text-slate-500">Theory and equations</p>
          </div>
        </Link>
        <Link
          to={ROUTES.LABS}
          className="flex items-center gap-3 bg-white rounded-xl shadow-sm border border-slate-200 p-4 hover:shadow-md hover:border-primary-300 transition-all focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500"
        >
          <div className="w-10 h-10 rounded-lg bg-slate-100 flex items-center justify-center flex-shrink-0">
            <svg className="w-5 h-5 text-slate-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
            </svg>
          </div>
          <div>
            <h3 className="text-sm font-semibold text-slate-800">Lab Packs</h3>
            <p className="text-xs text-slate-500">Guided activities</p>
          </div>
        </Link>
        <Link
          to={ROUTES.ABOUT}
          className="flex items-center gap-3 bg-white rounded-xl shadow-sm border border-slate-200 p-4 hover:shadow-md hover:border-primary-300 transition-all focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500"
        >
          <div className="w-10 h-10 rounded-lg bg-slate-100 flex items-center justify-center flex-shrink-0">
            <svg className="w-5 h-5 text-slate-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
            </svg>
          </div>
          <div>
            <h3 className="text-sm font-semibold text-slate-800">About</h3>
            <p className="text-xs text-slate-500">Citation and details</p>
          </div>
        </Link>
      </div>

      {/* Footer note */}
      <div className="mt-10 text-center text-xs text-slate-400">
        <p>
          Validated against PySCF 2.11.0 | Rust + WebAssembly | MIT License |{' '}
          <a
            href="https://github.com/ExaPsi/IQCP"
            target="_blank"
            rel="noopener noreferrer"
            className="hover:text-slate-600 underline"
          >
            GitHub
          </a>
        </p>
      </div>
    </div>
  );
}

export default Home;
