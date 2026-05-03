/**
 * Labs page - Displays Lab Pack worksheets for guided laboratory exercises.
 *
 * Supports multiple lab packs with a tab selector:
 * - Lab Pack #1: "From Boys Functions to SCF Convergence" (Modules C, D, E)
 * - Lab Pack #2: "3D Visualization, PES Scanning & Orbital Exploration" (Module E)
 * - Lab Pack #3: "Basis Functions, Integrals & Electron Density" (Modules A, B, E)
 *
 * Uses URL search param `?pack=1`, `?pack=2`, or `?pack=3` to select the active pack,
 * defaulting to the overview if no pack is specified.
 *
 * @module pages/Labs
 */

import { useSearchParams } from 'react-router-dom';
import { WorksheetViewer } from '@/components/labs';

// Import worksheet content as raw strings
import labpack1Content from '../../../../content/labpack1/worksheet-student.md?raw';
import labpack2Content from '../../../../content/labpack2/student-worksheet.md?raw';
import labpack3Content from '../../../../content/labpack3/student-worksheet.md?raw';

/**
 * Lab pack metadata for the selector.
 */
const LAB_PACKS = [
  {
    id: '1',
    labpack: 'labpack1',
    title: 'Lab Pack #1: From Boys Functions to SCF Convergence',
    subtitle: 'From Boys Functions to SCF Convergence',
    description:
      'Explore Boys function evaluation regimes, Rys quadrature accuracy vs order tradeoffs, and SCF convergence with DIIS acceleration.',
    modules: 'Modules C, D, E',
    content: labpack1Content,
    pdfFilename: 'labpack1-worksheet',
    phase: 'Phase 1',
    time: '60-90 min',
  },
  {
    id: '2',
    labpack: 'labpack2',
    title: 'Lab Pack #2: 3D Visualization, PES Scanning & Orbital Exploration',
    subtitle: '3D Visualization, PES Scanning & Orbital Exploration',
    description:
      'View molecules in 3D, scan potential energy surfaces to locate equilibria, visualize molecular orbital isosurfaces, and compare basis sets.',
    modules: 'Module E',
    content: labpack2Content,
    pdfFilename: 'labpack2-worksheet',
    phase: 'Phase 2',
    time: '~60 min',
    prereq: 'Lab Pack #1',
  },
  {
    id: '3',
    labpack: 'labpack3',
    title: 'Lab Pack #3: Basis Functions, Integrals & Electron Density',
    subtitle: 'Basis Functions, Integrals & Electron Density',
    description:
      'Examine contracted Gaussian basis functions, inspect one- and two-electron integral matrices, trace Fock matrix construction, and visualize electron density isosurfaces and difference maps.',
    modules: 'Modules A, B, E',
    content: labpack3Content,
    pdfFilename: 'labpack3-worksheet',
    phase: 'Phase 3',
    time: '~60 min',
    prereq: 'Lab Packs #1-2',
  },
] as const;

/**
 * Labs - Page component for Lab Pack materials.
 *
 * Shows a lab pack selector when no pack is chosen, or the
 * WorksheetViewer for the selected pack.
 */
function Labs() {
  const [searchParams, setSearchParams] = useSearchParams();
  const packId = searchParams.get('pack');

  const selectedPack = LAB_PACKS.find((p) => p.id === packId);

  // Show the selected lab pack
  if (selectedPack) {
    return (
      <div className="max-w-7xl mx-auto px-4">
        {/* Back to overview */}
        <button
          type="button"
          onClick={() => setSearchParams({})}
          className="mb-4 text-sm text-blue-600 hover:text-blue-800 flex items-center gap-1"
        >
          <span aria-hidden="true">&larr;</span> All Lab Packs
        </button>
        <WorksheetViewer
          labpack={selectedPack.labpack}
          title={selectedPack.title}
          content={selectedPack.content}
          pdfFilename={selectedPack.pdfFilename}
        />
      </div>
    );
  }

  // Lab pack overview / selector
  return (
    <div className="max-w-4xl mx-auto px-4">
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-slate-800 mb-2">
          Lab Packs
        </h1>
        <p className="text-slate-600">
          Guided laboratory exercises for teaching quantum chemistry concepts
          with IQCP. Each lab pack includes a student worksheet with interactive
          deep links to the computational modules.
        </p>
      </div>

      <div className="grid grid-cols-1 gap-5">
        {LAB_PACKS.map((pack) => (
          <button
            key={pack.id}
            type="button"
            onClick={() => setSearchParams({ pack: pack.id })}
            className="text-left bg-white rounded-xl shadow-sm border border-slate-200 p-6 hover:shadow-md hover:border-blue-300 transition-all group"
          >
            <div className="flex items-center gap-2 mb-2">
              <span className="text-xs font-medium text-blue-600 bg-blue-50 px-2 py-0.5 rounded-full">
                {pack.phase}
              </span>
              <span className="text-xs font-medium text-slate-500 bg-slate-100 px-2 py-0.5 rounded-full">
                {pack.modules}
              </span>
              <span className="text-xs text-slate-400">
                {pack.time}
              </span>
            </div>
            <h2 className="text-lg font-semibold text-slate-800 group-hover:text-blue-700 mb-1">
              Lab Pack #{pack.id}: {pack.subtitle}
            </h2>
            <p className="text-sm text-slate-600">
              {pack.description}
            </p>
            {'prereq' in pack && (
              <p className="text-xs text-slate-400 mt-2">
                Prerequisite: {pack.prereq}
              </p>
            )}
            <div className="mt-3 text-sm text-blue-600 font-medium group-hover:underline">
              Open worksheet &rarr;
            </div>
          </button>
        ))}
      </div>
    </div>
  );
}

export default Labs;
