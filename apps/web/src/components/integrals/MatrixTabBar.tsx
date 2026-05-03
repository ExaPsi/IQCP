/**
 * MatrixTabBar - Tab navigation for switching between integral matrices.
 *
 * Provides keyboard-accessible tab buttons for S, T, V, and H^core
 * matrices with ARIA roles and proper focus management.
 *
 * @module components/integrals/MatrixTabBar
 * @see US-055 Integral Matrix Heatmap
 */

import React, { useCallback, useRef } from 'react';
import type { MatrixTab } from '../../stores/integralStore';

// ============================================================================
// Types
// ============================================================================

export interface MatrixTabBarProps {
  /** Currently active tab */
  activeTab: MatrixTab;
  /** Callback when a tab is selected */
  onTabChange: (tab: MatrixTab) => void;
}

// ============================================================================
// Tab Configuration
// ============================================================================

interface TabConfig {
  id: MatrixTab;
  label: React.ReactNode;
  fullName: string;
  description: string;
}

const TABS: TabConfig[] = [
  {
    id: 'S',
    label: 'S',
    fullName: 'Overlap',
    description: 'Overlap integrals between basis functions',
  },
  {
    id: 'T',
    label: 'T',
    fullName: 'Kinetic',
    description: 'Kinetic energy integrals',
  },
  {
    id: 'V',
    label: 'V',
    fullName: 'Nuclear',
    description: 'Nuclear attraction integrals',
  },
  {
    id: 'Hcore',
    label: <span>H<sup>core</sup></span>,
    fullName: 'Core Hamiltonian',
    description: 'T + V: one-electron part of the Fock matrix',
  },
];

// ============================================================================
// Component
// ============================================================================

/**
 * Tab bar for switching between S, T, V, and H^core matrices.
 *
 * Implements the WAI-ARIA Tabs Pattern:
 * - Tab buttons have role="tab" and aria-selected
 * - Left/Right arrow keys navigate between tabs
 * - Home/End keys jump to first/last tab
 */
export function MatrixTabBar({ activeTab, onTabChange }: MatrixTabBarProps) {
  const tabRefs = useRef<(HTMLButtonElement | null)[]>([]);

  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent, currentIndex: number) => {
      let newIndex = currentIndex;

      switch (e.key) {
        case 'ArrowRight':
          e.preventDefault();
          newIndex = (currentIndex + 1) % TABS.length;
          break;
        case 'ArrowLeft':
          e.preventDefault();
          newIndex = (currentIndex - 1 + TABS.length) % TABS.length;
          break;
        case 'Home':
          e.preventDefault();
          newIndex = 0;
          break;
        case 'End':
          e.preventDefault();
          newIndex = TABS.length - 1;
          break;
        default:
          return;
      }

      onTabChange(TABS[newIndex].id);
      tabRefs.current[newIndex]?.focus();
    },
    [onTabChange],
  );

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-2">
      <h3 className="text-xs font-semibold text-slate-400 uppercase tracking-wider px-2 pb-2">
        Matrix
      </h3>
      <div
        role="tablist"
        aria-label="Integral matrix selection"
        className="flex gap-1"
      >
        {TABS.map((tab, index) => {
          const isActive = activeTab === tab.id;
          return (
            <button
              key={tab.id}
              ref={(el) => { tabRefs.current[index] = el; }}
              role="tab"
              id={`matrix-tab-${tab.id}`}
              aria-selected={isActive}
              aria-controls={`matrix-panel-${tab.id}`}
              tabIndex={isActive ? 0 : -1}
              title={`${tab.fullName}: ${tab.description}`}
              className={`flex-1 px-3 py-2 text-sm font-medium rounded-lg transition-colors
                focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500 focus-visible:ring-offset-1
                ${
                  isActive
                    ? 'bg-primary-100 text-primary-800 border border-primary-300'
                    : 'text-slate-600 hover:bg-slate-100 border border-transparent'
                }`}
              onClick={() => onTabChange(tab.id)}
              onKeyDown={(e) => handleKeyDown(e, index)}
            >
              <span className="block font-semibold">{tab.label}</span>
              <span className="block text-xs opacity-75 mt-0.5">{tab.fullName}</span>
            </button>
          );
        })}
      </div>
    </div>
  );
}
