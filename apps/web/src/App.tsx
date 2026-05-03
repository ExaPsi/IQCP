import { lazy, Suspense } from 'react';
import { Routes, Route, Navigate, useLocation } from 'react-router-dom';
import Header from '@/components/Header';
import Navigation from '@/components/Navigation';
import Footer from '@/components/Footer';
import { ErrorBoundary } from '@/components/ErrorBoundary';
import Home from '@/pages/Home';
import ModuleA from '@/pages/ModuleA';
import ModuleB from '@/pages/ModuleB';
import ModuleC from '@/pages/ModuleC';
import Labs from '@/pages/Labs';
import Reference from '@/pages/Reference';
import About from '@/pages/About';
import { ROUTES } from '@/lib/routes';

/** Module A: Basis Explorer (code-split -- only loaded when route is visited) */
const ModuleD = lazy(() => import('@/pages/ModuleD'));

/** Module B: Integral Inspector (code-split -- only loaded when route is visited) */
const ModuleE = lazy(() => import('@/pages/ModuleE'));

/**
 * Redirect component that preserves query parameters.
 *
 * React Router's Navigate component does not automatically carry
 * query parameters from the current URL to the target. This wrapper
 * appends the current search string to the redirect target so that
 * deep links like `/boys?run=X` correctly redirect to `/v1/boys?run=X`.
 */
function RedirectWithSearch({ to }: { to: string }) {
  const location = useLocation();
  return <Navigate to={`${to}${location.search}`} replace />;
}

/**
 * Root application component.
 * Provides the main layout structure with header, sidebar navigation, and content area.
 * All routes are configured here for the SPA shell.
 *
 * URL Versioning (Phase 2, US-036):
 * - All module and content routes use `/v1/` prefix for citation stability
 * - Unversioned paths redirect to versioned paths, preserving query parameters
 * - Deep links from manuscript and Lab Pack #1 continue to work
 */
function App() {
  return (
    <div className="min-h-screen bg-slate-50 flex flex-col">
      <Header />
      <div className="flex flex-1">
        <Navigation />
        <main className="flex-1 p-6" role="main">
          <ErrorBoundary>
            <Routes>
              {/* Home page (unversioned) */}
              <Route path="/" element={<Home />} />

              {/* === Versioned routes (v1) === */}
              <Route path="/v1/boys" element={<ModuleA />} />
              <Route path="/v1/rys" element={<ModuleB />} />
              <Route path="/v1/scf" element={<ModuleC />} />
              <Route
                path="/v1/basis"
                element={
                  <Suspense
                    fallback={
                      <div className="flex items-center justify-center py-12">
                        <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-primary-600" />
                        <span className="ml-3 text-slate-600">Loading Basis Explorer...</span>
                      </div>
                    }
                  >
                    <ModuleD />
                  </Suspense>
                }
              />
              <Route
                path="/v1/integrals"
                element={
                  <Suspense
                    fallback={
                      <div className="flex items-center justify-center py-12">
                        <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-primary-600" />
                        <span className="ml-3 text-slate-600">Loading Integral Inspector...</span>
                      </div>
                    }
                  >
                    <ModuleE />
                  </Suspense>
                }
              />
              <Route path="/v1/labs" element={<Labs />} />
              <Route path="/v1/reference" element={<Reference />} />
              <Route path="/v1/about" element={<About />} />

              {/* === Backward-compatible redirects === */}
              {/* Unversioned paths redirect to /v1/ equivalents.
                  Query parameters (deep link state) are preserved
                  by RedirectWithSearch. */}
              <Route path="/boys" element={<RedirectWithSearch to={ROUTES.BOYS} />} />
              <Route path="/rys" element={<RedirectWithSearch to={ROUTES.RYS} />} />
              <Route path="/scf" element={<RedirectWithSearch to={ROUTES.SCF} />} />
              <Route path="/basis" element={<RedirectWithSearch to={ROUTES.BASIS} />} />
              <Route path="/integrals" element={<RedirectWithSearch to={ROUTES.INTEGRALS} />} />
              <Route path="/labs" element={<RedirectWithSearch to={ROUTES.LABS} />} />
              <Route path="/reference" element={<RedirectWithSearch to={ROUTES.REFERENCE} />} />
              <Route path="/about" element={<RedirectWithSearch to={ROUTES.ABOUT} />} />
            </Routes>
          </ErrorBoundary>
        </main>
      </div>
      <Footer />
    </div>
  );
}

export default App;
