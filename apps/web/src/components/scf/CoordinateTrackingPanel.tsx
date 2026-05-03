/**
 * CoordinateTrackingPanel - Visualize non-scanned coordinate response during PES scans.
 *
 * Shows how internal coordinates (bonds, angles, dihedrals) change across a
 * PES scan. For rigid scans, tracked coordinates appear as flat horizontal
 * lines (confirming rigidity). For relaxed scans, tracked coordinates show
 * the physical response of the molecular geometry to the constraint.
 *
 * Features:
 * - Multi-select dropdown grouped by coordinate category (bonds, angles, dihedrals)
 * - Plotly traces with Okabe-Ito colorblind-safe palette
 * - Dynamic y-axis label based on tracked coordinate types
 * - Shared x-axis with PesCurvePlot for visual synchronization
 * - Excludes the scanned coordinate from the selection list
 *
 * @module components/scf/CoordinateTrackingPanel
 * @see US-083 Extended PES Plot + Coordinate Tracking Panel
 */

import { useMemo, useState, useCallback } from 'react';
import { PlotPanel } from '../common/PlotPanel';
import type { Data, Layout } from 'plotly.js';
import type { PesInternalPoint, InternalCoordinateSnapshot } from '../../worker/protocol';
import type { CoordinateType } from './CoordinateSelector';

// ============================================================================
// Constants
// ============================================================================

/**
 * Okabe-Ito colorblind-safe palette for tracked coordinate traces.
 *
 * These colors are distinguishable by individuals with the most common
 * forms of color vision deficiency (deuteranopia and protanopia).
 *
 * @see https://jfly.uni-koeln.de/color/
 */
const TRACKING_COLORS = [
  '#0072B2', // blue
  '#D55E00', // vermillion
  '#009E73', // bluish green
  '#CC79A7', // reddish purple
  '#F0E442', // yellow
  '#56B4E9', // sky blue
];

/**
 * Dynamic x-axis labels for the tracking plot, matching PesCurvePlot.
 */
const X_AXIS_LABELS: Record<CoordinateType, string> = {
  bond: 'Bond Length (bohr)',
  angle: 'Angle (degrees)',
  dihedral: 'Dihedral Angle (degrees)',
};

/**
 * Map atomic number Z to element symbol for labels.
 */
const ELEMENT_SYMBOLS: Record<number, string> = {
  1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
  11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
};

// ============================================================================
// Types
// ============================================================================

/**
 * Coordinate key format: "type:idx1-idx2[-idx3[-idx4]]"
 * Examples: "bond:0-1", "angle:1-0-2", "dihedral:0-1-2-3"
 */
type CoordinateKey = string;

/**
 * A selectable internal coordinate with label and category.
 */
interface SelectableCoordinate {
  /** Unique key identifying this coordinate */
  key: CoordinateKey;
  /** Human-readable label, e.g., "O0-H1 bond" */
  label: string;
  /** Category for grouped display */
  category: 'bond' | 'angle' | 'dihedral';
  /** Atom indices defining this coordinate */
  indices: number[];
}

/**
 * Props for CoordinateTrackingPanel.
 */
export interface CoordinateTrackingPanelProps {
  /** Scan points with internal coordinate snapshots */
  scanPoints: PesInternalPoint[];
  /** Type of the scanned coordinate (for x-axis labeling) */
  coordinateType: CoordinateType;
  /** Atom indices of the scanned coordinate (excluded from dropdown) */
  scannedAtomIndices: number[];
  /** Atom data as [Z, x, y, z] for generating labels */
  atoms: [number, number, number, number][];
  /** Currently selected coordinate keys (multi-select, from store) */
  selectedCoordinates: string[];
  /** Handler for selection changes (updates store) */
  onSelectedChange: (keys: string[]) => void;
  /** Minimum height of the plot in pixels (default: 280) */
  minHeight?: number;
}

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Build the key for a scanned coordinate to exclude it from the dropdown.
 */
function buildScannedKey(coordinateType: CoordinateType, atomIndices: number[]): string {
  return `${coordinateType}:${atomIndices.join('-')}`;
}

/**
 * Extract available coordinates from the first scan point's snapshot.
 *
 * All scan points have the same set of internal coordinates (bonds, angles,
 * dihedrals), differing only in their values. We extract the coordinate
 * definitions from the first point.
 *
 * @param points - Scan points with internal coordinate snapshots
 * @param atoms - Atom data [Z, x, y, z] for generating element labels
 * @param excludeKey - Key of the scanned coordinate to exclude from the list
 * @returns Array of selectable coordinates grouped by category
 */
export function extractAvailableCoordinates(
  points: PesInternalPoint[],
  atoms: [number, number, number, number][],
  excludeKey: string,
): SelectableCoordinate[] {
  if (points.length === 0) return [];

  const snapshot = points[0].internal_coordinates;
  if (!snapshot) return [];

  const result: SelectableCoordinate[] = [];

  // Helper to get element label for an atom index
  const atomLabel = (idx: number): string => {
    const z = atoms[idx]?.[0] ?? 0;
    return `${ELEMENT_SYMBOLS[z] ?? 'X'}${idx}`;
  };

  for (const [i, j] of snapshot.bonds) {
    const key = `bond:${i}-${j}`;
    if (key === excludeKey) continue;
    result.push({
      key,
      label: `${atomLabel(i)}-${atomLabel(j)} bond`,
      category: 'bond',
      indices: [i, j],
    });
  }

  for (const [i, j, k] of snapshot.angles) {
    const key = `angle:${i}-${j}-${k}`;
    if (key === excludeKey) continue;
    result.push({
      key,
      label: `${atomLabel(i)}-${atomLabel(j)}-${atomLabel(k)} angle`,
      category: 'angle',
      indices: [i, j, k],
    });
  }

  for (const [i, j, k, l] of snapshot.dihedrals) {
    const key = `dihedral:${i}-${j}-${k}-${l}`;
    if (key === excludeKey) continue;
    result.push({
      key,
      label: `${atomLabel(i)}-${atomLabel(j)}-${atomLabel(k)}-${atomLabel(l)} dihedral`,
      category: 'dihedral',
      indices: [i, j, k, l],
    });
  }

  return result;
}

/**
 * Extract a specific coordinate's value from each scan point's snapshot.
 *
 * @param scanPoints - All scan points from the PES scan
 * @param coordKey - Coordinate key identifying which coordinate to extract
 * @returns Object with x (scanned coordinate values) and y (tracked coordinate values)
 */
export function extractCoordinateValues(
  scanPoints: PesInternalPoint[],
  coordKey: CoordinateKey,
): { x: number[]; y: number[] } {
  const sepIdx = coordKey.indexOf(':');
  const type = coordKey.slice(0, sepIdx);
  const indices = coordKey.slice(sepIdx + 1).split('-').map(Number);

  const x: number[] = [];
  const y: number[] = [];

  for (const point of scanPoints) {
    if (!point.internal_coordinates) continue;

    x.push(point.coordinate_value);

    const snapshot: InternalCoordinateSnapshot = point.internal_coordinates;
    let value: number | null = null;

    if (type === 'bond') {
      const entry = snapshot.bonds.find(
        ([a, b]) => a === indices[0] && b === indices[1],
      );
      value = entry ? entry[2] : null;
    } else if (type === 'angle') {
      const entry = snapshot.angles.find(
        ([a, b, c]) => a === indices[0] && b === indices[1] && c === indices[2],
      );
      // Convert radians to degrees for display
      value = entry ? entry[3] * (180 / Math.PI) : null;
    } else if (type === 'dihedral') {
      const entry = snapshot.dihedrals.find(
        ([a, b, c, d]) =>
          a === indices[0] && b === indices[1] &&
          c === indices[2] && d === indices[3],
      );
      // Convert radians to degrees for display
      value = entry ? entry[4] * (180 / Math.PI) : null;
    }

    if (value !== null) {
      y.push(value);
    }
  }

  return { x, y };
}

/**
 * Convert scanned coordinate x values from internal units to display units.
 *
 * Bond lengths remain in bohr. Angles and dihedrals are converted
 * from radians to degrees.
 */
function convertXValues(
  values: number[],
  coordinateType: CoordinateType,
): number[] {
  if (coordinateType === 'bond') return values;
  return values.map((v) => v * (180 / Math.PI));
}

/**
 * Determine the y-axis label based on the categories of tracked coordinates.
 *
 * - All bonds: "Distance (bohr)"
 * - All angles/dihedrals: "Angle (degrees)"
 * - Mixed: "Coordinate Value (mixed units)"
 */
function getYAxisLabel(
  selectedCoords: SelectableCoordinate[],
): string {
  if (selectedCoords.length === 0) return 'Coordinate Value';

  const categories = new Set(selectedCoords.map((c) => c.category));

  if (categories.size === 1) {
    const cat = categories.values().next().value;
    if (cat === 'bond') return 'Distance (bohr)';
    return 'Angle (degrees)';
  }

  return 'Coordinate Value (mixed units)';
}

// ============================================================================
// Sub-Components
// ============================================================================

/**
 * Multi-select dropdown for choosing internal coordinates to track.
 *
 * Coordinates are grouped by category (Bonds, Angles, Dihedrals)
 * with group headers.
 */
function CoordinateMultiSelect({
  available,
  selected,
  onSelectedChange,
}: {
  available: SelectableCoordinate[];
  selected: string[];
  onSelectedChange: (keys: string[]) => void;
}) {
  const [isOpen, setIsOpen] = useState(false);

  const toggleOpen = useCallback(() => setIsOpen((prev) => !prev), []);

  const handleToggle = useCallback(
    (key: string) => {
      if (selected.includes(key)) {
        onSelectedChange(selected.filter((k) => k !== key));
      } else {
        onSelectedChange([...selected, key]);
      }
    },
    [selected, onSelectedChange],
  );

  const handleClearAll = useCallback(() => {
    onSelectedChange([]);
  }, [onSelectedChange]);

  const handleSelectAll = useCallback(() => {
    onSelectedChange(available.map((c) => c.key));
  }, [available, onSelectedChange]);

  // Group by category
  const grouped = useMemo(() => {
    const groups: Record<string, SelectableCoordinate[]> = {};
    for (const coord of available) {
      if (!groups[coord.category]) {
        groups[coord.category] = [];
      }
      groups[coord.category].push(coord);
    }
    return groups;
  }, [available]);

  const categoryLabels: Record<string, string> = {
    bond: 'Bonds',
    angle: 'Angles',
    dihedral: 'Dihedrals',
  };

  const selectedCount = selected.length;

  return (
    <div className="relative">
      {/* Dropdown trigger */}
      <button
        type="button"
        onClick={toggleOpen}
        className="w-full px-3 py-2 border border-slate-300 rounded-lg text-sm text-left bg-white hover:bg-slate-50 transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 flex items-center justify-between"
        aria-expanded={isOpen}
        aria-haspopup="listbox"
        aria-label="Select coordinates to track"
      >
        <span className="text-slate-700">
          {selectedCount === 0
            ? 'Select coordinates to track...'
            : `${selectedCount} coordinate${selectedCount !== 1 ? 's' : ''} selected`}
        </span>
        <svg
          className={`w-4 h-4 text-slate-400 transition-transform ${isOpen ? 'rotate-180' : ''}`}
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {/* Dropdown menu */}
      {isOpen && (
        <div
          className="absolute z-10 mt-1 w-full bg-white border border-slate-200 rounded-lg shadow-lg max-h-60 overflow-y-auto"
          role="listbox"
          aria-multiselectable="true"
          aria-label="Available internal coordinates"
        >
          {/* Bulk actions */}
          <div className="flex items-center justify-between px-3 py-1.5 border-b border-slate-100 bg-slate-50">
            <button
              type="button"
              onClick={handleSelectAll}
              className="text-xs text-blue-600 hover:text-blue-800 font-medium"
            >
              Select all
            </button>
            <button
              type="button"
              onClick={handleClearAll}
              className="text-xs text-slate-500 hover:text-slate-700 font-medium"
            >
              Clear all
            </button>
          </div>

          {/* Grouped checkboxes */}
          {Object.entries(grouped).map(([category, coords]) => (
            <div key={category}>
              <div className="px-3 py-1.5 text-xs font-semibold text-slate-500 uppercase tracking-wider bg-slate-50 border-b border-slate-100">
                {categoryLabels[category] ?? category}
              </div>
              {coords.map((coord) => {
                const isChecked = selected.includes(coord.key);
                const colorIdx = selected.indexOf(coord.key);
                return (
                  <label
                    key={coord.key}
                    className="flex items-center gap-2 px-3 py-1.5 hover:bg-blue-50 cursor-pointer transition-colors"
                    role="option"
                    aria-selected={isChecked}
                  >
                    <input
                      type="checkbox"
                      checked={isChecked}
                      onChange={() => handleToggle(coord.key)}
                      className="h-3.5 w-3.5 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
                    />
                    {isChecked && colorIdx >= 0 && (
                      <span
                        className="inline-block w-2.5 h-2.5 rounded-full flex-shrink-0"
                        style={{
                          backgroundColor: TRACKING_COLORS[colorIdx % TRACKING_COLORS.length],
                        }}
                      />
                    )}
                    <span className="text-sm text-slate-700">{coord.label}</span>
                  </label>
                );
              })}
            </div>
          ))}

          {available.length === 0 && (
            <div className="px-3 py-4 text-sm text-slate-400 text-center">
              No coordinates available
            </div>
          )}
        </div>
      )}
    </div>
  );
}

// ============================================================================
// Main Component
// ============================================================================

/**
 * Panel for tracking how non-scanned internal coordinates respond during a
 * PES scan.
 *
 * Features:
 * - Multi-select dropdown with coordinates grouped by category (AC4, AC6)
 * - Plotly line traces with Okabe-Ito palette for each selected coordinate (AC5)
 * - Dynamic y-axis label based on tracked coordinate types
 * - Shared x-axis domain with PesCurvePlot for visual synchronization (AC5)
 * - Flat lines for rigid scan tracked coordinates (AC7)
 *
 * @example
 * ```tsx
 * <CoordinateTrackingPanel
 *   scanPoints={pesInternalResult.points}
 *   coordinateType="angle"
 *   scannedAtomIndices={[1, 0, 2]}
 *   atoms={workerAtoms}
 *   selectedCoordinates={trackedKeys}
 *   onSelectedChange={setTrackedCoordinates}
 * />
 * ```
 */
export function CoordinateTrackingPanel({
  scanPoints,
  coordinateType,
  scannedAtomIndices,
  atoms,
  selectedCoordinates,
  onSelectedChange,
  minHeight = 280,
}: CoordinateTrackingPanelProps) {
  // Build the key for the scanned coordinate to exclude it
  const scannedKey = useMemo(
    () => buildScannedKey(coordinateType, scannedAtomIndices),
    [coordinateType, scannedAtomIndices],
  );

  // Extract available coordinates from the first scan point
  const availableCoordinates = useMemo(
    () => extractAvailableCoordinates(scanPoints, atoms, scannedKey),
    [scanPoints, atoms, scannedKey],
  );

  // Map selected keys to their SelectableCoordinate objects
  const selectedCoords = useMemo(
    () => availableCoordinates.filter((c) => selectedCoordinates.includes(c.key)),
    [availableCoordinates, selectedCoordinates],
  );

  // Y-axis label
  const yAxisLabel = useMemo(
    () => getYAxisLabel(selectedCoords),
    [selectedCoords],
  );

  // Build Plotly traces for selected coordinates
  const plotData = useMemo((): Data[] => {
    if (selectedCoords.length === 0) return [];

    return selectedCoords.map((coord, idx): Data => {
      const { x, y } = extractCoordinateValues(scanPoints, coord.key);
      const xDisplay = convertXValues(x, coordinateType);
      const color = TRACKING_COLORS[idx % TRACKING_COLORS.length];

      // Determine value unit for hover
      const valueUnit = coord.category === 'bond' ? 'bohr' : 'deg';

      return {
        x: xDisplay,
        y,
        type: 'scatter',
        mode: 'lines+markers',
        name: coord.label,
        line: {
          color,
          width: 2,
        },
        marker: {
          color,
          size: 6,
        },
        hovertemplate: `${coord.label}: %{y:.4f} ${valueUnit}<extra></extra>`,
      };
    });
  }, [selectedCoords, scanPoints, coordinateType]);

  // Build layout
  const layout = useMemo((): Partial<Layout> => ({
    title: {
      text: 'Coordinate Response',
      font: { size: 14 },
    },
    xaxis: {
      title: { text: X_AXIS_LABELS[coordinateType], font: { size: 12 } },
    },
    yaxis: {
      title: { text: yAxisLabel, font: { size: 12 } },
    },
    margin: { t: 50, r: 30, b: 50, l: 80 },
    showlegend: true,
    legend: {
      x: 1,
      y: 1,
      xanchor: 'right',
      yanchor: 'top',
      font: { size: 10 },
      bgcolor: 'rgba(255,255,255,0.8)',
      bordercolor: '#e2e8f0',
      borderwidth: 1,
    },
    hovermode: 'x unified' as const,
  }), [coordinateType, yAxisLabel]);

  // No scan data yet
  if (scanPoints.length === 0) return null;

  // No internal coordinates available
  if (availableCoordinates.length === 0) return null;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
      {/* Header */}
      <div className="px-4 py-3 bg-slate-50 border-b border-slate-200">
        <h3 className="text-sm font-semibold text-slate-700">
          Track Coordinates
        </h3>
        <p className="text-xs text-slate-500 mt-0.5">
          Select internal coordinates to see how they respond during the scan
        </p>
      </div>

      {/* Coordinate selector */}
      <div className="p-4 space-y-4">
        <CoordinateMultiSelect
          available={availableCoordinates}
          selected={selectedCoordinates}
          onSelectedChange={onSelectedChange}
        />

        {/* Plot (only when coordinates are selected) */}
        {selectedCoords.length > 0 && (
          <PlotPanel
            data={plotData}
            layout={layout}
            minHeight={minHeight}
            ariaLabel={`Coordinate response plot showing ${selectedCoords.length} tracked coordinate${selectedCoords.length !== 1 ? 's' : ''} during PES scan`}
          />
        )}

        {/* Empty state when nothing selected */}
        {selectedCoords.length === 0 && (
          <div className="py-8 text-center text-sm text-slate-400">
            Select one or more coordinates above to see how they change during the scan
          </div>
        )}
      </div>
    </div>
  );
}

export type { CoordinateKey, SelectableCoordinate };
