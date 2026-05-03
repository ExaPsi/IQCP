import { Link } from 'react-router-dom';
import { ROUTES, CURRENT_VERSION } from '@/lib/routes';

/**
 * Footer component displaying version and citation link.
 *
 * Shows the IQCP version derived from `package.json` (injected at build time
 * via Vite's `define` config) and the URL version prefix.
 *
 * @module components/Footer
 */
function Footer() {
  const appVersion = typeof __APP_VERSION__ !== 'undefined' ? __APP_VERSION__ : '0.1.0';

  return (
    <footer
      className="bg-white border-t border-slate-200 py-3 px-6"
      role="contentinfo"
    >
      <div className="container mx-auto flex items-center justify-between text-xs text-slate-500">
        <div className="flex items-center gap-4">
          <span>
            IQCP {CURRENT_VERSION} (app {appVersion})
          </span>
          <Link
            to={ROUTES.ABOUT}
            className="hover:text-slate-700 underline decoration-slate-300 hover:decoration-slate-500 transition-colors"
          >
            About & Citation
          </Link>
        </div>
        <a
          href="https://github.com/ExaPsi/IQCP"
          target="_blank"
          rel="noopener noreferrer"
          className="hover:text-slate-700 transition-colors"
        >
          GitHub
        </a>
      </div>
    </footer>
  );
}

export default Footer;
