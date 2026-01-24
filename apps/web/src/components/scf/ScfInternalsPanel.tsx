/**
 * ScfInternalsPanel - Container for SCF internals visualization.
 *
 * Combines matrix heatmaps, orbital energies table, and diagnostic panel
 * into a cohesive internals inspection view. Includes tabbed navigation
 * for switching between different SCF matrices.
 *
 * Layout:
 * ```
 * +------------------------------------------+
 * |  [S] [Hcore] [F] [D]    (matrix tabs)    |
 * +------------------------------------------+
 * |        Matrix Heatmap                    |
 * |        (selected matrix)                 |
 * +------------------------------------------+
 * |  Orbital Energies Table                  |
 * +------------------------------------------+
 * |  Diagnostic Panel                        |
 * +------------------------------------------+
 * ```
 *
 * @module components/scf/ScfInternalsPanel
 */

import { useState, useMemo } from 'react';
import { MatrixHeatmap, type MatrixColorscale } from './MatrixHeatmap';
import { OrbitalEnergiesTable } from './OrbitalEnergiesTable';
import { DiagnosticPanel } from './DiagnosticPanel';
import { Math } from '../common/Math';
import type { ScfMatrices, OrbitalEnergies, ConvergenceProfile } from '../../worker/protocol';

/**
 * Matrix tab configuration.
 */
interface MatrixTab {
  /** Tab identifier */
  id: 'S' | 'Hcore' | 'F' | 'P';
  /** Display label */
  label: string;
  /** Full title for heatmap */
  title: string;
  /** Appropriate colorscale for this matrix type */
  colorscale: MatrixColorscale;
  /** Brief description for tooltip */
  description: string;
}

/**
 * Configuration for all matrix tabs.
 * Note: descriptions use plain text here; KaTeX rendering happens in the component.
 */
const MATRIX_TABS: MatrixTab[] = [
  {
    id: 'S',
    label: 'S',
    title: 'Overlap Matrix (S)',
    colorscale: 'Blues',
    description: 'Overlap integrals between basis functions. Diagonal elements are 1.',
  },
  {
    id: 'Hcore',
    label: 'Hcore',
    title: 'Core Hamiltonian (Hcore)',
    colorscale: 'RdBu',
    description: 'One-electron integrals: kinetic energy + nuclear attraction.',
  },
  {
    id: 'F',
    label: 'F',
    title: 'Fock Matrix (F)',
    colorscale: 'RdBu',
    description: 'Effective one-electron operator including electron-electron repulsion.',
  },
  {
    id: 'P',
    label: 'P',
    title: 'Density Matrix (P)',
    colorscale: 'Greens',
    description: 'One-electron density matrix from occupied MO coefficients.',
  },
];

/**
 * Get math formula for each matrix type.
 */
function getMatrixFormula(tabId: MatrixTab['id']): string {
  switch (tabId) {
    case 'S':
      return String.raw`S_{\mu\nu} = \langle \phi_\mu | \phi_\nu \rangle`;
    case 'Hcore':
      return String.raw`H^{\text{core}}_{\mu\nu} = \langle \phi_\mu | \hat{T} + \hat{V}_{ne} | \phi_\nu \rangle`;
    case 'F':
      return String.raw`F_{\mu\nu} = H^{\text{core}}_{\mu\nu} + G_{\mu\nu}`;
    case 'P':
      return String.raw`P_{\mu\nu} = 2\sum_i^{\text{occ}} C_{\mu i} C_{\nu i}`;
  }
}

/**
 * Props for ScfInternalsPanel component.
 */
export interface ScfInternalsPanelProps {
  /** SCF matrices for visualization */
  matrices: ScfMatrices;
  /** Orbital energies data */
  orbitalEnergies: OrbitalEnergies;
  /** Whether SCF converged */
  converged: boolean;
  /** Number of iterations */
  iterations: number;
  /** Final total energy in Hartree */
  energy: number;
  /** Convergence profile used */
  convergenceProfile: ConvergenceProfile;
  /** Whether DIIS was enabled */
  diisEnabled: boolean;
  /** Whether computation was aborted */
  aborted?: boolean;
  /** Additional CSS classes */
  className?: string;
}

/**
 * Get matrix data array for the selected tab.
 */
function getMatrixData(matrices: ScfMatrices, tabId: MatrixTab['id']): number[] {
  switch (tabId) {
    case 'S':
      return matrices.sMatrix;
    case 'Hcore':
      return matrices.hCore;
    case 'F':
      return matrices.fockMatrix;
    case 'P':
      return matrices.densityMatrix;
  }
}

/**
 * Tab button component with accessibility support.
 */
function TabButton({
  tab,
  isSelected,
  onClick,
}: {
  tab: MatrixTab;
  isSelected: boolean;
  onClick: () => void;
}) {
  return (
    <button
      type="button"
      role="tab"
      aria-selected={isSelected}
      aria-controls={`matrix-panel-${tab.id}`}
      id={`matrix-tab-${tab.id}`}
      className={`px-4 py-2 text-sm font-medium transition-colors relative ${
        isSelected
          ? 'text-blue-700 bg-blue-50'
          : 'text-slate-600 hover:text-slate-800 hover:bg-slate-50'
      }`}
      onClick={onClick}
      title={tab.description}
    >
      {tab.label}
      {isSelected && (
        <span className="absolute bottom-0 left-0 right-0 h-0.5 bg-blue-600" />
      )}
    </button>
  );
}

/**
 * ScfInternalsPanel - Complete SCF internals visualization.
 *
 * Provides tabbed access to SCF matrices, orbital energies, and
 * computational diagnostics in a unified interface.
 *
 * @example
 * ```tsx
 * <ScfInternalsPanel
 *   matrices={result.matrices}
 *   orbitalEnergies={result.orbitalEnergies}
 *   converged={result.converged}
 *   iterations={result.iterations}
 *   energy={result.energy}
 *   convergenceProfile="medium"
 *   diisEnabled={true}
 * />
 * ```
 */
export function ScfInternalsPanel({
  matrices,
  orbitalEnergies,
  converged,
  iterations,
  energy,
  convergenceProfile,
  diisEnabled,
  aborted = false,
  className = '',
}: ScfInternalsPanelProps) {
  // Track selected matrix tab
  const [selectedTab, setSelectedTab] = useState<MatrixTab['id']>('S');

  // Get current tab configuration
  const currentTab = useMemo(
    () => MATRIX_TABS.find((t) => t.id === selectedTab) || MATRIX_TABS[0],
    [selectedTab]
  );

  // Get matrix data for current tab
  const currentMatrixData = useMemo(
    () => getMatrixData(matrices, selectedTab),
    [matrices, selectedTab]
  );

  return (
    <div className={`space-y-4 ${className}`}>
      {/* Matrix Visualization Section */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
        {/* Tab Bar */}
        <div
          className="border-b border-slate-200 flex"
          role="tablist"
          aria-label="SCF Matrix selection"
        >
          {MATRIX_TABS.map((tab) => (
            <TabButton
              key={tab.id}
              tab={tab}
              isSelected={selectedTab === tab.id}
              onClick={() => setSelectedTab(tab.id)}
            />
          ))}
        </div>

        {/* Matrix Description with Formula */}
        <div className="px-4 py-2 bg-slate-50 border-b border-slate-100 flex flex-wrap items-center gap-x-4 gap-y-1">
          <p className="text-xs text-slate-600">{currentTab.description}</p>
          <span className="text-xs text-slate-500">
            <Math>{getMatrixFormula(selectedTab)}</Math>
          </span>
        </div>

        {/* Matrix Heatmap Panel */}
        <div
          role="tabpanel"
          id={`matrix-panel-${selectedTab}`}
          aria-labelledby={`matrix-tab-${selectedTab}`}
          className="p-4"
        >
          <MatrixHeatmap
            key={selectedTab}
            data={currentMatrixData}
            nbf={matrices.nbf}
            title={currentTab.title}
            colorscale={currentTab.colorscale}
            minHeight={350}
          />
        </div>

        {/* Matrix Info Footer */}
        <div className="px-4 py-2 bg-slate-50 border-t border-slate-100 text-xs text-slate-500 flex justify-between">
          <span>
            Dimension: {matrices.nbf} x {matrices.nbf} ({matrices.nbf * matrices.nbf} elements)
          </span>
          <span>
            Data: Row-major storage
          </span>
        </div>
      </div>

      {/* Orbital Energies Table */}
      <OrbitalEnergiesTable
        energies={orbitalEnergies.energies}
        nOccupied={orbitalEnergies.nOccupied}
      />

      {/* Diagnostic Panel */}
      <DiagnosticPanel
        converged={converged}
        iterations={iterations}
        energy={energy}
        convergenceProfile={convergenceProfile}
        diisEnabled={diisEnabled}
        aborted={aborted}
      />
    </div>
  );
}

export default ScfInternalsPanel;
