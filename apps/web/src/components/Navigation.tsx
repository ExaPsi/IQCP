import { NavLink } from 'react-router-dom';

/**
 * Navigation item component for module links.
 * Displays label and description with keyboard-accessible focus indicators.
 */
interface NavItemProps {
  to: string;
  label: string;
  description: string;
}

function NavItem({ to, label, description }: NavItemProps) {
  return (
    <li>
      <NavLink
        to={to}
        className={({ isActive }) =>
          `block px-4 py-3 rounded-lg transition-colors
          focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500 focus-visible:ring-offset-2
          ${
            isActive
              ? 'bg-primary-100 text-primary-800 border-l-4 border-primary-600'
              : 'text-slate-600 hover:bg-slate-100'
          }`
        }
      >
        <div className="font-medium">{label}</div>
        <div className="text-sm text-slate-500">{description}</div>
      </NavLink>
    </li>
  );
}

/**
 * Sidebar navigation component for module selection.
 * Provides keyboard-accessible navigation to Module A, B, and C.
 * Uses semantic HTML with proper ARIA labels.
 */
function Navigation() {
  return (
    <aside
      className="w-64 bg-white border-r border-slate-200 min-h-screen p-4"
      role="complementary"
      aria-label="Module navigation"
    >
      <nav aria-label="Module navigation">
        <section>
          <h2 className="text-xs font-semibold text-slate-400 uppercase tracking-wider px-4 py-2">
            Modules
          </h2>
          <ul className="space-y-1" role="list">
            <NavItem
              to="/boys"
              label="Module A: Boys Function"
              description="Explore F_m(T) evaluation"
            />
            <NavItem
              to="/rys"
              label="Module B: Rys Quadrature"
              description="Roots and weights"
            />
            <NavItem
              to="/scf"
              label="Module C: SCF Sandbox"
              description="RHF with DIIS"
            />
          </ul>
        </section>
      </nav>
    </aside>
  );
}

export default Navigation;
