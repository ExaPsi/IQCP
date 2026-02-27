import { Link } from 'react-router-dom';

interface ModuleCardProps {
  to: string;
  title: string;
  subtitle: string;
  description: string;
}

function ModuleCard({ to, title, subtitle, description }: ModuleCardProps) {
  return (
    <Link
      to={to}
      className="block bg-white rounded-xl shadow-sm border border-slate-200 p-6 hover:shadow-md hover:border-primary-300 transition-all focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500 focus-visible:ring-offset-2"
    >
      <h3 className="text-lg font-semibold text-slate-800">{title}</h3>
      <p className="text-primary-600 font-medium mt-1">{subtitle}</p>
      <p className="text-slate-600 mt-3">{description}</p>
    </Link>
  );
}

function Home() {
  return (
    <div className="max-w-4xl mx-auto">
      <div className="text-center mb-12">
        <h1 className="text-4xl font-bold text-slate-800 mb-4">
          Interactive Quantum Chemistry Playground
        </h1>
        <p className="text-xl text-slate-600 max-w-2xl mx-auto">
          Explore fundamental quantum chemistry algorithms through interactive
          visualizations. Understand Boys functions, Rys quadrature, and SCF
          convergence with transparent, reproducible computation.
        </p>
      </div>

      <div className="grid md:grid-cols-2 lg:grid-cols-3 gap-6">
        <ModuleCard
          to="/boys"
          title="Module A"
          subtitle="Boys Function Lab"
          description="Explore the Boys function F_m(T) with regime-based evaluation. Visualize series, recurrence, and asymptotic methods."
        />
        <ModuleCard
          to="/rys"
          title="Module B"
          subtitle="Rys Quadrature Lab"
          description="Compute Rys polynomial roots and weights. Understand order-error tradeoffs in Gaussian integral evaluation."
        />
        <ModuleCard
          to="/scf"
          title="Module C"
          subtitle="SCF Sandbox"
          description="Run restricted Hartree-Fock calculations with optional DIIS acceleration. Inspect iteration-by-iteration convergence."
        />
      </div>

      <div className="mt-12 bg-primary-50 rounded-xl p-6 border border-primary-100">
        <h2 className="text-lg font-semibold text-primary-800 mb-2">
          About IQCP
        </h2>
        <p className="text-primary-700">
          IQCP is an educational web application designed for teaching quantum
          chemistry concepts. All heavy computation runs in-browser using
          WebAssembly, ensuring responsive interactions and reproducible
          results. Designed for use in J. Chem. Educ. coursework.
        </p>
      </div>
    </div>
  );
}

export default Home;
