import { Link, NavLink } from 'react-router-dom';
import { ROUTES } from '@/lib/routes';

/**
 * Top navigation link component with keyboard accessibility.
 * Uses focus-visible for clear focus indicators.
 */
interface TopNavLinkProps {
  to: string;
  children: React.ReactNode;
}

function TopNavLink({ to, children }: TopNavLinkProps) {
  return (
    <NavLink
      to={to}
      className={({ isActive }) =>
        `px-3 py-2 rounded-md text-sm font-medium transition-colors
        focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-white focus-visible:ring-offset-2 focus-visible:ring-offset-primary-800
        ${
          isActive
            ? 'bg-primary-700 text-white'
            : 'text-primary-200 hover:bg-primary-700 hover:text-white'
        }`
      }
    >
      {children}
    </NavLink>
  );
}

/**
 * Header component with IQCP branding and primary navigation.
 * Includes links to Labs, Reference, and About pages.
 * Fully keyboard accessible with focus-visible indicators.
 */
function Header() {
  return (
    <header className="bg-primary-800 text-white shadow-md" role="banner">
      <div className="container mx-auto px-6 py-4">
        <div className="flex items-center justify-between">
          {/* Branding */}
          <Link
            to="/"
            className="flex items-center space-x-3 rounded focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-white focus-visible:ring-offset-2 focus-visible:ring-offset-primary-800"
          >
            <span className="text-2xl font-bold tracking-tight">IQCP</span>
            <span className="hidden sm:inline-block text-primary-200 text-sm">
              Interactive Quantum Chemistry Playground
            </span>
          </Link>

          {/* Primary Navigation */}
          <nav
            className="flex items-center space-x-1"
            aria-label="Primary navigation"
          >
            <TopNavLink to={ROUTES.LABS}>Labs</TopNavLink>
            <TopNavLink to={ROUTES.REFERENCE}>Reference</TopNavLink>
            <TopNavLink to={ROUTES.ABOUT}>About</TopNavLink>
            <a
              href="https://github.com/ExaPsi/IQCP"
              target="_blank"
              rel="noopener noreferrer"
              className="ml-2 px-3 py-2 rounded-md text-sm font-medium text-primary-200 hover:bg-primary-700 hover:text-white transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-white focus-visible:ring-offset-2 focus-visible:ring-offset-primary-800"
              aria-label="GitHub repository (opens in new tab)"
            >
              GitHub
            </a>
          </nav>
        </div>
      </div>
    </header>
  );
}

export default Header;
