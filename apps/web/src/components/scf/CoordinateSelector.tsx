/**
 * CoordinateSelector - Select internal coordinate type and atoms for PES scan.
 *
 * Provides a controlled component for choosing coordinate type
 * (Bond/Angle/Dihedral), selecting atom indices via dropdowns populated
 * from the current geometry, and displaying the current coordinate value.
 *
 * Also includes scan range inputs (min/max/nPoints) with sensible defaults
 * per coordinate type.
 *
 * @module components/scf/CoordinateSelector
 * @see US-082 Coordinate Selector + Scan Mode UI
 */

import { useMemo } from 'react';

// ============================================================================
// Types
// ============================================================================

/**
 * Type of internal coordinate for PES scanning.
 * Matches the coordinateType field on PesScanInternalRequest.
 */
export type CoordinateType = 'bond' | 'angle' | 'dihedral';

/**
 * Number of atom indices required per coordinate type.
 */
const COORDINATE_ATOM_COUNT: Record<CoordinateType, number> = {
  bond: 2,
  angle: 3,
  dihedral: 4,
};

/**
 * Human-readable labels for coordinate types.
 */
const COORDINATE_LABELS: Record<CoordinateType, string> = {
  bond: 'Bond Length',
  angle: 'Bond Angle',
  dihedral: 'Dihedral Angle',
};

/**
 * Unit labels for coordinate types.
 */
const COORDINATE_UNITS: Record<CoordinateType, string> = {
  bond: 'bohr',
  angle: 'degrees',
  dihedral: 'degrees',
};

/**
 * Atom label descriptions for each position in the coordinate definition.
 *
 * For angles the central atom is the second (index 1) position.
 * For dihedrals the bond is between atoms j-k (positions 1-2).
 */
const ATOM_POSITION_LABELS: Record<CoordinateType, string[]> = {
  bond: ['Atom i', 'Atom j'],
  angle: ['Atom i', 'Central atom j', 'Atom k'],
  dihedral: ['Atom i', 'Atom j', 'Atom k', 'Atom l'],
};

/**
 * Map atomic number Z to element symbol for labels.
 */
const ELEMENT_SYMBOLS: Record<number, string> = {
  1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
  11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
};

// ============================================================================
// Geometry Computations
// ============================================================================

/**
 * Compute the current value of an internal coordinate from atom positions.
 *
 * @param atoms - Array of [Z, x, y, z] per atom (coordinates in bohr)
 * @param coordinateType - Type of coordinate to compute
 * @param atomIndices - Atom indices defining the coordinate
 * @returns The coordinate value, or null if the selection is invalid.
 *          Bond: distance in bohr. Angle/Dihedral: value in degrees.
 */
export function computeCoordinateValue(
  atoms: [number, number, number, number][],
  coordinateType: CoordinateType,
  atomIndices: number[],
): number | null {
  const requiredCount = COORDINATE_ATOM_COUNT[coordinateType];
  if (atomIndices.length < requiredCount) return null;

  // Check all indices are valid and distinct
  for (let k = 0; k < requiredCount; k++) {
    const idx = atomIndices[k];
    if (idx < 0 || idx >= atoms.length) return null;
    for (let j = k + 1; j < requiredCount; j++) {
      if (idx === atomIndices[j]) return null;
    }
  }

  if (coordinateType === 'bond') {
    const [, x1, y1, z1] = atoms[atomIndices[0]];
    const [, x2, y2, z2] = atoms[atomIndices[1]];
    const dx = x2 - x1;
    const dy = y2 - y1;
    const dz = z2 - z1;
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  if (coordinateType === 'angle') {
    // Angle at central atom j between atoms i-j-k
    const [, xi, yi, zi] = atoms[atomIndices[0]];
    const [, xj, yj, zj] = atoms[atomIndices[1]];
    const [, xk, yk, zk] = atoms[atomIndices[2]];

    const v1x = xi - xj, v1y = yi - yj, v1z = zi - zj;
    const v2x = xk - xj, v2y = yk - yj, v2z = zk - zj;

    const dot = v1x * v2x + v1y * v2y + v1z * v2z;
    const len1 = Math.sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
    const len2 = Math.sqrt(v2x * v2x + v2y * v2y + v2z * v2z);

    if (len1 < 1e-10 || len2 < 1e-10) return null;

    // Clamp to avoid NaN from floating-point drift
    const cosAngle = Math.max(-1, Math.min(1, dot / (len1 * len2)));
    return Math.acos(cosAngle) * (180 / Math.PI);
  }

  if (coordinateType === 'dihedral') {
    // Dihedral angle for chain i-j-k-l
    const [, xi, yi, zi] = atoms[atomIndices[0]];
    const [, xj, yj, zj] = atoms[atomIndices[1]];
    const [, xk, yk, zk] = atoms[atomIndices[2]];
    const [, xl, yl, zl] = atoms[atomIndices[3]];

    // Vectors along the chain
    const b1x = xj - xi, b1y = yj - yi, b1z = zj - zi;
    const b2x = xk - xj, b2y = yk - yj, b2z = zk - zj;
    const b3x = xl - xk, b3y = yl - yk, b3z = zl - zk;

    // Normal vectors: n1 = b1 x b2, n2 = b2 x b3
    const n1x = b1y * b2z - b1z * b2y;
    const n1y = b1z * b2x - b1x * b2z;
    const n1z = b1x * b2y - b1y * b2x;

    const n2x = b2y * b3z - b2z * b3y;
    const n2y = b2z * b3x - b2x * b3z;
    const n2z = b2x * b3y - b2y * b3x;

    // Unit vector along b2
    const b2len = Math.sqrt(b2x * b2x + b2y * b2y + b2z * b2z);
    if (b2len < 1e-10) return null;
    const e2x = b2x / b2len, e2y = b2y / b2len, e2z = b2z / b2len;

    // m = n1 x e2
    const mx = n1y * e2z - n1z * e2y;
    const my = n1z * e2x - n1x * e2z;
    const mz = n1x * e2y - n1y * e2x;

    const x = n1x * n2x + n1y * n2y + n1z * n2z;
    const y = mx * n2x + my * n2y + mz * n2z;

    return Math.atan2(y, x) * (180 / Math.PI);
  }

  return null;
}

// ============================================================================
// Props
// ============================================================================

/**
 * Props for CoordinateSelector.
 *
 * This is a controlled component -- all state is managed by the parent
 * (PesScanPanel) which passes values and change handlers.
 */
export interface CoordinateSelectorProps {
  /** Currently selected coordinate type */
  coordinateType: CoordinateType;
  /** Handler for coordinate type changes */
  onCoordinateTypeChange: (type: CoordinateType) => void;
  /** Currently selected atom indices */
  atomIndices: number[];
  /** Handler for atom index changes */
  onAtomIndicesChange: (indices: number[]) => void;
  /** Available atoms from the current geometry as [Z, x, y, z] */
  atoms: [number, number, number, number][];
  /** Min value for scan range */
  valueMin: number;
  /** Handler for min value changes */
  onValueMinChange: (v: number) => void;
  /** Max value for scan range */
  valueMax: number;
  /** Handler for max value changes */
  onValueMaxChange: (v: number) => void;
  /** Number of scan points */
  nPoints: number;
  /** Handler for nPoints changes */
  onNPointsChange: (n: number) => void;
  /** Whether controls are disabled (scanning in progress, worker not ready) */
  disabled?: boolean;
}

// ============================================================================
// Component
// ============================================================================

/**
 * CoordinateSelector - Select internal coordinate and scan range for PES scan.
 *
 * Renders:
 * 1. Three radio/tab buttons for coordinate type (Bond / Angle / Dihedral)
 * 2. Atom index dropdowns (2, 3, or 4 depending on coordinate type)
 * 3. Current coordinate value display
 * 4. Scan range inputs (min, max, nPoints)
 *
 * Dihedral is only available when the molecule has >= 4 atoms.
 *
 * @example
 * ```tsx
 * <CoordinateSelector
 *   coordinateType={coordType}
 *   onCoordinateTypeChange={setCoordType}
 *   atomIndices={indices}
 *   onAtomIndicesChange={setIndices}
 *   atoms={workerAtoms}
 *   valueMin={min}
 *   onValueMinChange={setMin}
 *   valueMax={max}
 *   onValueMaxChange={setMax}
 *   nPoints={nPts}
 *   onNPointsChange={setNPts}
 * />
 * ```
 */
export function CoordinateSelector({
  coordinateType,
  onCoordinateTypeChange,
  atomIndices,
  onAtomIndicesChange,
  atoms,
  valueMin,
  onValueMinChange,
  valueMax,
  onValueMaxChange,
  nPoints,
  onNPointsChange,
  disabled = false,
}: CoordinateSelectorProps) {
  const nAtoms = atoms.length;

  // Atom labels for dropdown options: "O1", "H2", etc.
  // Matches the 3D viewer's atom label format (1-indexed, no parentheses).
  const atomOptions = useMemo(() => {
    return atoms.map((a, i) => {
      const symbol = ELEMENT_SYMBOLS[a[0]] ?? 'X';
      return { index: i, label: `${symbol}${i + 1}` };
    });
  }, [atoms]);

  // Current coordinate value
  const currentValue = useMemo(
    () => computeCoordinateValue(atoms, coordinateType, atomIndices),
    [atoms, coordinateType, atomIndices],
  );

  // Available coordinate types based on atom count
  const availableTypes = useMemo(() => {
    const types: CoordinateType[] = ['bond'];
    if (nAtoms >= 3) types.push('angle');
    if (nAtoms >= 4) types.push('dihedral');
    return types;
  }, [nAtoms]);

  // Required number of atom indices for current type
  const requiredAtomCount = COORDINATE_ATOM_COUNT[coordinateType];
  const positionLabels = ATOM_POSITION_LABELS[coordinateType];
  const unit = COORDINATE_UNITS[coordinateType];

  // Handle coordinate type change -- reset atom indices to first N atoms
  const handleTypeChange = (newType: CoordinateType) => {
    onCoordinateTypeChange(newType);
    const count = COORDINATE_ATOM_COUNT[newType];
    const newIndices = Array.from({ length: count }, (_, i) => Math.min(i, nAtoms - 1));
    onAtomIndicesChange(newIndices);
  };

  // Handle individual atom dropdown change
  const handleAtomChange = (position: number, newIdx: number) => {
    const updated = [...atomIndices];
    updated[position] = newIdx;
    onAtomIndicesChange(updated);
  };

  // Determine step size for range inputs
  const step = coordinateType === 'bond' ? 0.1 : 1;

  return (
    <div className="space-y-3">
      {/* Coordinate type selector (toggle buttons) */}
      <div>
        <label className="block text-xs font-medium text-slate-700 mb-1.5">
          Coordinate Type
        </label>
        <div className="flex rounded-lg bg-slate-100 p-1">
          {availableTypes.map((type) => (
            <button
              key={type}
              type="button"
              onClick={() => handleTypeChange(type)}
              disabled={disabled}
              className={`flex-1 px-2 py-1.5 rounded-md text-xs font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
                coordinateType === type
                  ? 'bg-white text-slate-800 shadow-sm'
                  : 'text-slate-600 hover:text-slate-800'
              } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
              aria-pressed={coordinateType === type}
            >
              {COORDINATE_LABELS[type]}
            </button>
          ))}
        </div>
      </div>

      {/* Atom index dropdowns */}
      <div>
        <label className="block text-xs font-medium text-slate-700 mb-1.5">
          Atoms
        </label>
        <div className={`grid gap-2 ${requiredAtomCount <= 2 ? 'grid-cols-2' : requiredAtomCount === 3 ? 'grid-cols-3' : 'grid-cols-2 sm:grid-cols-4'}`}>
          {Array.from({ length: requiredAtomCount }, (_, pos) => (
            <div key={pos}>
              <label
                htmlFor={`coord-atom-${pos}`}
                className="block text-[10px] text-slate-500 mb-0.5"
              >
                {positionLabels[pos]}
              </label>
              <select
                id={`coord-atom-${pos}`}
                value={atomIndices[pos] ?? 0}
                onChange={(e) => handleAtomChange(pos, parseInt(e.target.value, 10))}
                disabled={disabled}
                className="w-full px-2 py-1.5 border border-slate-300 rounded-lg text-sm text-slate-800 bg-white focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
                aria-label={`${positionLabels[pos]} for ${COORDINATE_LABELS[coordinateType]}`}
              >
                {atomOptions.map((opt) => (
                  <option key={opt.index} value={opt.index}>
                    {opt.label}
                  </option>
                ))}
              </select>
            </div>
          ))}
        </div>
      </div>

      {/* Current value display */}
      <div className="bg-slate-50 rounded-lg px-3 py-2 border border-slate-200">
        <div className="flex items-center justify-between">
          <span className="text-xs text-slate-500">
            Current {COORDINATE_LABELS[coordinateType]}
          </span>
          <span className="text-sm font-mono font-medium text-slate-800">
            {currentValue !== null ? (
              <>
                {coordinateType === 'bond'
                  ? currentValue.toFixed(4)
                  : currentValue.toFixed(2)
                }
                {' '}
                <span className="text-xs text-slate-500">{unit}</span>
              </>
            ) : (
              <span className="text-xs text-amber-600">
                Select distinct atoms
              </span>
            )}
          </span>
        </div>
      </div>

      {/* Scan range inputs */}
      <div className="grid grid-cols-2 gap-3">
        <div>
          <label
            htmlFor="coord-range-min"
            className="block text-xs font-medium text-slate-700 mb-1"
          >
            Min ({unit})
          </label>
          <input
            id="coord-range-min"
            type="number"
            step={step}
            value={valueMin}
            onChange={(e) => {
              const v = parseFloat(e.target.value);
              if (!Number.isNaN(v)) onValueMinChange(v);
            }}
            disabled={disabled}
            className="w-full px-2 py-1.5 border border-slate-300 rounded-lg text-sm text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
            aria-label={`Minimum ${COORDINATE_LABELS[coordinateType]} in ${unit}`}
          />
        </div>
        <div>
          <label
            htmlFor="coord-range-max"
            className="block text-xs font-medium text-slate-700 mb-1"
          >
            Max ({unit})
          </label>
          <input
            id="coord-range-max"
            type="number"
            step={step}
            value={valueMax}
            onChange={(e) => {
              const v = parseFloat(e.target.value);
              if (!Number.isNaN(v)) onValueMaxChange(v);
            }}
            disabled={disabled}
            className="w-full px-2 py-1.5 border border-slate-300 rounded-lg text-sm text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
            aria-label={`Maximum ${COORDINATE_LABELS[coordinateType]} in ${unit}`}
          />
        </div>
      </div>

      <div>
        <label
          htmlFor="coord-npoints"
          className="block text-xs font-medium text-slate-700 mb-1"
        >
          Number of Points
        </label>
        <input
          id="coord-npoints"
          type="number"
          min={2}
          max={100}
          step={1}
          value={nPoints}
          onChange={(e) => {
            const v = parseInt(e.target.value, 10);
            if (!Number.isNaN(v)) onNPointsChange(v);
          }}
          disabled={disabled}
          className="w-full px-2 py-1.5 border border-slate-300 rounded-lg text-sm text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
          aria-label="Number of PES scan points"
        />
        <p className="mt-0.5 text-xs text-slate-400">Range: 2-100</p>
      </div>
    </div>
  );
}

// ============================================================================
// Exports
// ============================================================================

export { COORDINATE_ATOM_COUNT, COORDINATE_LABELS, COORDINATE_UNITS, ELEMENT_SYMBOLS };
