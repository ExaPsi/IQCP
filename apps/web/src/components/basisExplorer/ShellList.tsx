/**
 * ShellList - Displays the shells for the selected element + basis combination.
 *
 * Shows each contracted shell with its angular momentum label,
 * number of primitives, and a brief description. Shells are clickable
 * to set the selectedShell in the store (for future use in US-050).
 *
 * Shell data is sourced from qc-core's builtin.rs via WASM (US-049 AC6),
 * ensuring a single source of truth for basis set information.
 *
 * @module components/basisExplorer/ShellList
 */

import { useMemo } from 'react';
import { useBasisStore } from '../../stores/basisStore';
import { useBasisInfo } from '../../hooks/useBasisInfo';
import { getElement, isElementSupported, BASIS_SETS } from './basisData';
import type { ShellInfo } from './basisData';
import { deriveShellDisplayInfo } from './shellDisplayInfo';

// ============================================================================
// Components
// ============================================================================

/**
 * Color class for angular momentum type badges.
 */
function getAngularMomentumBadge(am: number): { bg: string; text: string; label: string } {
  switch (am) {
    case 0:
      return { bg: 'bg-blue-100', text: 'text-blue-800', label: 's' };
    case 1:
      return { bg: 'bg-green-100', text: 'text-green-800', label: 'p' };
    case 2:
      return { bg: 'bg-purple-100', text: 'text-purple-800', label: 'd' };
    default:
      return { bg: 'bg-slate-100', text: 'text-slate-800', label: '?' };
  }
}

/**
 * Individual shell row component.
 */
interface ShellRowProps {
  shell: ShellInfo;
  index: number;
  isSelected: boolean;
  onSelect: (index: number) => void;
}

function ShellRow({ shell, index, isSelected, onSelect }: ShellRowProps) {
  const badge = getAngularMomentumBadge(shell.angularMomentum);

  return (
    <button
      type="button"
      onClick={() => onSelect(index)}
      aria-pressed={isSelected}
      aria-label={`Shell ${shell.label}: ${shell.nPrimitives} primitive${shell.nPrimitives !== 1 ? 's' : ''}, ${shell.description}`}
      className={`
        w-full flex items-center gap-3 px-4 py-3 rounded-lg border text-left
        transition-all duration-150
        focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500 focus-visible:ring-offset-1
        ${
          isSelected
            ? 'bg-primary-50 border-primary-300 shadow-sm'
            : 'bg-white border-slate-200 hover:bg-slate-50 hover:border-slate-300'
        }
      `}
    >
      {/* Shell index */}
      <span className="text-xs text-slate-400 font-mono w-5 text-right shrink-0">
        {index + 1}
      </span>

      {/* Angular momentum badge */}
      <span
        className={`inline-flex items-center justify-center w-8 h-8 rounded-full text-sm font-bold shrink-0 ${badge.bg} ${badge.text}`}
      >
        {badge.label}
      </span>

      {/* Shell info */}
      <div className="flex-1 min-w-0">
        <div className="flex items-center gap-2">
          <span className="font-semibold text-slate-800">{shell.label}</span>
          <span className="text-xs text-slate-500">
            {shell.nPrimitives} primitive{shell.nPrimitives !== 1 ? 's' : ''}
          </span>
        </div>
        <p className="text-xs text-slate-500 truncate">{shell.description}</p>
      </div>

      {/* Selection indicator */}
      {isSelected && (
        <span className="text-primary-600 shrink-0" aria-hidden="true">
          <svg className="w-5 h-5" fill="currentColor" viewBox="0 0 20 20">
            <path
              fillRule="evenodd"
              d="M16.704 4.153a.75.75 0 01.143 1.052l-8 10.5a.75.75 0 01-1.127.075l-4.5-4.5a.75.75 0 011.06-1.06l3.894 3.893 7.48-9.817a.75.75 0 011.05-.143z"
              clipRule="evenodd"
            />
          </svg>
        </span>
      )}
    </button>
  );
}

/**
 * Shell list panel.
 *
 * Displays all shells for the currently selected element and basis set.
 * Shell data is fetched from WASM (qc-core builtin.rs) via the Web Worker.
 */
export function ShellList() {
  const selectedElement = useBasisStore((state) => state.selectedElement);
  const selectedBasis = useBasisStore((state) => state.selectedBasis);
  const selectedShell = useBasisStore((state) => state.selectedShell);
  const setShell = useBasisStore((state) => state.setShell);

  const element = getElement(selectedElement);
  const basisLabel = BASIS_SETS.find((b) => b.value === selectedBasis)?.label ?? selectedBasis;

  // Fetch shell data from WASM via the Web Worker (US-049 AC6)
  const { shells: wasmShells, loading, error } = useBasisInfo(selectedElement, selectedBasis);

  // Derive display-friendly shell info from WASM data
  const shells = useMemo(() => {
    if (!wasmShells) return null;
    return deriveShellDisplayInfo(wasmShells, selectedElement);
  }, [wasmShells, selectedElement]);

  // Check if element is supported
  if (!isElementSupported(selectedElement)) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
        <h2 className="text-lg font-semibold text-slate-800 mb-4">
          Shell Structure
        </h2>
        <div className="text-center py-8">
          <div className="inline-flex items-center justify-center w-12 h-12 rounded-full bg-slate-100 mb-3">
            <svg
              className="w-6 h-6 text-slate-400"
              fill="none"
              viewBox="0 0 24 24"
              strokeWidth={1.5}
              stroke="currentColor"
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                d="M12 9v3.75m9-.75a9 9 0 11-18 0 9 9 0 0118 0zm-9 3.75h.008v.008H12v-.008z"
              />
            </svg>
          </div>
          <p className="text-slate-600 font-medium">
            {element?.name ?? `Element Z=${selectedElement}`} is not yet supported
          </p>
          <p className="text-sm text-slate-500 mt-1">
            Basis data for elements beyond Ar will be added in a future update.
          </p>
        </div>
      </div>
    );
  }

  // Loading state
  if (loading) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
        <h2 className="text-lg font-semibold text-slate-800 mb-4">
          Shell Structure
        </h2>
        <div className="text-center py-8">
          <div className="inline-flex items-center justify-center w-12 h-12 rounded-full bg-slate-100 mb-3 animate-pulse">
            <svg
              className="w-6 h-6 text-slate-400"
              fill="none"
              viewBox="0 0 24 24"
              strokeWidth={1.5}
              stroke="currentColor"
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                d="M9.75 3.104v5.714a2.25 2.25 0 01-.659 1.591L5 14.5M9.75 3.104c-.251.023-.501.05-.75.082m.75-.082a24.301 24.301 0 014.5 0m0 0v5.714c0 .597.237 1.17.659 1.591L19.8 15.3M14.25 3.104c.251.023.501.05.75.082M19.8 15.3l-1.57.393A9.065 9.065 0 0112 15a9.065 9.065 0 00-6.23.693L5 14.5m14.8.8l1.402 1.402c1.232 1.232.65 3.318-1.067 3.611A48.309 48.309 0 0112 21c-2.773 0-5.491-.235-8.135-.687-1.718-.293-2.3-2.379-1.067-3.61L5 14.5"
              />
            </svg>
          </div>
          <p className="text-slate-600">
            Loading shell data...
          </p>
        </div>
      </div>
    );
  }

  // Error state
  if (error) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
        <h2 className="text-lg font-semibold text-slate-800 mb-4">
          Shell Structure
        </h2>
        <div className="text-center py-8">
          <div className="inline-flex items-center justify-center w-12 h-12 rounded-full bg-red-100 mb-3">
            <svg
              className="w-6 h-6 text-red-500"
              fill="none"
              viewBox="0 0 24 24"
              strokeWidth={1.5}
              stroke="currentColor"
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                d="M12 9v3.75m-9.303 3.376c-.866 1.5.217 3.374 1.948 3.374h14.71c1.73 0 2.813-1.874 1.948-3.374L13.949 3.378c-.866-1.5-3.032-1.5-3.898 0L2.697 16.126zM12 15.75h.007v.008H12v-.008z"
              />
            </svg>
          </div>
          <p className="text-red-600 font-medium">
            Failed to load shell data
          </p>
          <p className="text-sm text-slate-500 mt-1">
            {error}
          </p>
        </div>
      </div>
    );
  }

  // No data (shouldn't happen for supported elements, but handle gracefully)
  if (!shells) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
        <h2 className="text-lg font-semibold text-slate-800 mb-4">
          Shell Structure
        </h2>
        <div className="text-center py-8">
          <p className="text-slate-600">
            No shell data available for {element?.symbol ?? '?'} with {basisLabel}.
          </p>
        </div>
      </div>
    );
  }

  // Count d-shells for spherical vs Cartesian distinction
  const nDShells = shells.filter((s) => s.angularMomentum === 2).length;

  // Count total basis functions (using spherical 5d by default):
  // s -> 1, p -> 3 (px, py, pz), d -> 5 (spherical)
  const totalBasisFunctions = shells.reduce((sum, shell) => {
    switch (shell.angularMomentum) {
      case 0: return sum + 1;  // s
      case 1: return sum + 3;  // px, py, pz
      case 2: return sum + 5;  // 5 spherical d
      default: return sum + (2 * shell.angularMomentum + 1);
    }
  }, 0);

  // Cartesian count: d-shells contribute 6 instead of 5
  const totalBasisFunctionsCartesian = totalBasisFunctions + nDShells;

  const totalPrimitives = shells.reduce((sum, shell) => {
    const multiplier = shell.angularMomentum === 0 ? 1 : shell.angularMomentum === 1 ? 3 : 5;
    return sum + shell.nPrimitives * multiplier;
  }, 0);

  const handleSelect = (index: number) => {
    // Toggle selection
    setShell(selectedShell === index ? null : index);
  };

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
      <div className="flex items-center justify-between mb-4">
        <h2 className="text-lg font-semibold text-slate-800">
          Shell Structure
        </h2>
        <span className="text-sm text-slate-500">
          {element?.symbol} / {basisLabel}
        </span>
      </div>

      {/* Summary stats */}
      <div className="grid grid-cols-3 gap-3 mb-4">
        <div className="bg-slate-50 rounded-lg p-3 text-center">
          <p className="text-2xl font-bold text-slate-800">{shells.length}</p>
          <p className="text-xs text-slate-500">Shells</p>
        </div>
        <div className="bg-slate-50 rounded-lg p-3 text-center">
          <p className="text-2xl font-bold text-slate-800">{totalBasisFunctions}</p>
          <p className="text-xs text-slate-500">Basis Functions</p>
        </div>
        <div className="bg-slate-50 rounded-lg p-3 text-center">
          <p className="text-2xl font-bold text-slate-800">{totalPrimitives}</p>
          <p className="text-xs text-slate-500">Primitives</p>
        </div>
      </div>

      {/* Spherical vs Cartesian note for bases with d-shells */}
      {nDShells > 0 && (
        <div className="mb-4 px-3 py-2 bg-purple-50 rounded-lg border border-purple-100">
          <p className="text-xs text-purple-800">
            <strong>Cartesian:</strong> {totalBasisFunctionsCartesian} basis functions (6d) &nbsp;|&nbsp;{' '}
            <strong>Spherical:</strong> {totalBasisFunctions} basis functions (5d)
          </p>
          <p className="text-[10px] text-purple-600 mt-0.5">
            {nDShells} d-shell{nDShells !== 1 ? 's' : ''}: Cartesian uses 6 components (xx, yy, zz, xy, xz, yz), spherical uses 5 (m = -2..+2).
          </p>
        </div>
      )}

      {/* Shell list */}
      <div className="space-y-2" role="listbox" aria-label={`Shells for ${element?.name} in ${basisLabel}`}>
        {shells.map((shell, index) => (
          <ShellRow
            key={`${index}-${shell.label}`}
            shell={shell}
            index={index}
            isSelected={selectedShell === index}
            onSelect={handleSelect}
          />
        ))}
      </div>

      {/* Educational note */}
      <div className="mt-4 pt-3 border-t border-slate-100">
        <p className="text-xs text-slate-500">
          Each shell is a contracted Gaussian function built from one or more primitive
          Gaussians. Click a shell to select it for detailed inspection.
        </p>
      </div>
    </div>
  );
}
