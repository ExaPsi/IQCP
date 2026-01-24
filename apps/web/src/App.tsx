import { Routes, Route } from 'react-router-dom';
import Header from '@/components/Header';
import Navigation from '@/components/Navigation';
import { ErrorBoundary } from '@/components/ErrorBoundary';
import Home from '@/pages/Home';
import ModuleB from '@/pages/ModuleB';
import ModuleC from '@/pages/ModuleC';
import ModuleE from '@/pages/ModuleE';
import Labs from '@/pages/Labs';
import Reference from '@/pages/Reference';
import About from '@/pages/About';

/**
 * Root application component.
 * Provides the main layout structure with header, sidebar navigation, and content area.
 * All routes are configured here for the SPA shell.
 */
function App() {
  return (
    <div className="min-h-screen bg-slate-50">
      <Header />
      <div className="flex">
        <Navigation />
        <main className="flex-1 p-6" role="main">
          <ErrorBoundary>
            <Routes>
              {/* Home page with module cards */}
              <Route path="/" element={<Home />} />

              {/* Module pages */}
              <Route path="/boys" element={<ModuleB />} />
              <Route path="/rys" element={<ModuleC />} />
              <Route path="/scf" element={<ModuleE />} />

              {/* Content pages */}
              <Route path="/labs" element={<Labs />} />
              <Route path="/reference" element={<Reference />} />
              <Route path="/about" element={<About />} />
            </Routes>
          </ErrorBoundary>
        </main>
      </div>
    </div>
  );
}

export default App;
