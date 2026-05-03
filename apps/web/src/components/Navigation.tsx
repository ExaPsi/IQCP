import { useState, useCallback } from 'react';
import { NavLink } from 'react-router-dom';
import { ROUTES } from '@/lib/routes';

/**
 * Navigation item component for module links.
 * Shows full label + description when expanded, icon-only when collapsed.
 */
interface NavItemProps {
  to: string;
  label: string;
  shortLabel: string;
  description: string;
  collapsed: boolean;
}

function NavItem({ to, label, shortLabel, description, collapsed }: NavItemProps) {
  return (
    <li>
      <NavLink
        to={to}
        title={collapsed ? `${label}: ${description}` : undefined}
        className={({ isActive }) =>
          `block rounded-lg transition-colors
          focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500 focus-visible:ring-offset-2
          ${collapsed ? 'px-2 py-2.5 text-center' : 'px-4 py-3'}
          ${
            isActive
              ? 'bg-primary-100 text-primary-800 border-l-4 border-primary-600'
              : 'text-slate-600 hover:bg-slate-100'
          }`
        }
      >
        {collapsed ? (
          <div className="font-bold text-sm">{shortLabel}</div>
        ) : (
          <>
            <div className="font-medium">{label}</div>
            <div className="text-sm text-slate-500">{description}</div>
          </>
        )}
      </NavLink>
    </li>
  );
}

/**
 * Hamburger / close icon for the sidebar toggle.
 */
function MenuIcon({ collapsed }: { collapsed: boolean }) {
  if (collapsed) {
    // Hamburger menu icon (3 lines)
    return (
      <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
        <path fillRule="evenodd" d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z" clipRule="evenodd" />
      </svg>
    );
  }
  // X close icon
  return (
    <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
      <path fillRule="evenodd" d="M4.293 4.293a1 1 0 011.414 0L10 8.586l4.293-4.293a1 1 0 111.414 1.414L11.414 10l4.293 4.293a1 1 0 01-1.414 1.414L10 11.414l-4.293 4.293a1 1 0 01-1.414-1.414L8.586 10 4.293 5.707a1 1 0 010-1.414z" clipRule="evenodd" />
    </svg>
  );
}

/**
 * Collapsible sidebar navigation for module selection.
 *
 * - Toggle button at top-right corner of sidebar header (always visible)
 * - Expanded: full module names with descriptions (w-64)
 * - Collapsed: compact letter labels (w-16)
 * - State persists in localStorage
 */
function Navigation() {
  const [collapsed, setCollapsed] = useState(() => {
    try {
      return localStorage.getItem('iqcp-nav-collapsed') === 'true';
    } catch {
      return false;
    }
  });

  const toggle = useCallback(() => {
    setCollapsed((prev) => {
      const next = !prev;
      try {
        localStorage.setItem('iqcp-nav-collapsed', String(next));
      } catch {
        // localStorage unavailable
      }
      return next;
    });
  }, []);

  return (
    <aside
      className={`${collapsed ? 'w-16' : 'w-64'} bg-white border-r border-slate-200 min-h-screen transition-all duration-200 flex flex-col`}
      role="complementary"
      aria-label="Module navigation"
    >
      {/* Sidebar header with toggle button */}
      <div className={`flex items-center border-b border-slate-200 ${collapsed ? 'justify-center px-2 py-3' : 'justify-between px-4 py-3'}`}>
        {!collapsed && (
          <h2 className="text-xs font-semibold text-slate-400 uppercase tracking-wider">
            Modules
          </h2>
        )}
        <button
          onClick={toggle}
          className="p-1.5 rounded-md text-slate-400 hover:text-slate-700 hover:bg-slate-100 transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500"
          aria-label={collapsed ? 'Open navigation panel' : 'Close navigation panel'}
          title={collapsed ? 'Open navigation panel' : 'Close navigation panel'}
        >
          <MenuIcon collapsed={collapsed} />
        </button>
      </div>

      {/* Navigation links */}
      <nav aria-label="Module navigation" className="flex-1 p-2">
        <ul className="space-y-1" role="list">
          <NavItem
            to={ROUTES.BASIS}
            label="Module A: Basis Explorer"
            shortLabel="A"
            description="Gaussian basis sets"
            collapsed={collapsed}
          />
          <NavItem
            to={ROUTES.INTEGRALS}
            label="Module B: Integral Inspector"
            shortLabel="B"
            description="Matrix heatmaps"
            collapsed={collapsed}
          />
          <NavItem
            to={ROUTES.BOYS}
            label="Module C: Boys Function"
            shortLabel="C"
            description="Explore Fₘ(T) evaluation"
            collapsed={collapsed}
          />
          <NavItem
            to={ROUTES.RYS}
            label="Module D: Rys Quadrature"
            shortLabel="D"
            description="Roots and weights"
            collapsed={collapsed}
          />
          <NavItem
            to={ROUTES.SCF}
            label="Module E: SCF Sandbox"
            shortLabel="E"
            description="RHF / DFT with DIIS"
            collapsed={collapsed}
          />
        </ul>
      </nav>
    </aside>
  );
}

export default Navigation;
