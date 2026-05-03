/**
 * ElementSelector - Compact periodic table for element selection (H-Ar).
 *
 * Displays elements H through Ar in an authentic periodic-table layout
 * sized to fit within the Module D sidebar panel.
 *
 * @module components/basisExplorer/ElementSelector
 */

import { useBasisStore } from '../../stores/basisStore';
import { ELEMENTS } from './basisData';

// ============================================================================
// Element Category Styling
// ============================================================================

type ElementCategory = 'alkali' | 'alkaline' | 'metalloid' | 'nonmetal' | 'halogen' | 'noble' | 'post-transition';

function getCategory(z: number): ElementCategory {
  if (z === 3 || z === 11) return 'alkali';
  if (z === 4 || z === 12) return 'alkaline';
  if (z === 5 || z === 14) return 'metalloid';
  if (z === 13) return 'post-transition';
  if (z === 1 || z === 6 || z === 7 || z === 8 || z === 15 || z === 16) return 'nonmetal';
  if (z === 9 || z === 17) return 'halogen';
  if (z === 2 || z === 10 || z === 18) return 'noble';
  return 'nonmetal';
}

const CAT_STYLE: Record<ElementCategory, string> = {
  alkali:            'bg-red-50 border-red-200 text-red-800 hover:bg-red-100',
  alkaline:          'bg-orange-50 border-orange-200 text-orange-800 hover:bg-orange-100',
  metalloid:         'bg-teal-50 border-teal-200 text-teal-800 hover:bg-teal-100',
  'post-transition': 'bg-cyan-50 border-cyan-200 text-cyan-800 hover:bg-cyan-100',
  nonmetal:          'bg-emerald-50 border-emerald-200 text-emerald-800 hover:bg-emerald-100',
  halogen:           'bg-sky-50 border-sky-200 text-sky-800 hover:bg-sky-100',
  noble:             'bg-violet-50 border-violet-200 text-violet-800 hover:bg-violet-100',
};

const CAT_DOT: Record<ElementCategory, string> = {
  alkali: 'bg-red-200', alkaline: 'bg-orange-200', metalloid: 'bg-teal-200',
  'post-transition': 'bg-cyan-200', nonmetal: 'bg-emerald-200', halogen: 'bg-sky-200', noble: 'bg-violet-200',
};

const CAT_LABEL: Record<ElementCategory, string> = {
  alkali: 'Alkali', alkaline: 'Alk. earth', metalloid: 'Metalloid',
  'post-transition': 'Post-trans.', nonmetal: 'Nonmetal', halogen: 'Halogen', noble: 'Noble gas',
};

// ============================================================================
// Compact periodic table layout for H-Ar only.
// Instead of 18 columns (mostly empty), use a condensed 8-column layout:
//
//   Row 1:  H                              He
//   Row 2:  Li  Be      B   C   N   O   F  Ne
//   Row 3:  Na  Mg      Al  Si  P   S   Cl Ar
//
// 8 visible columns with a spacer gap between col 2 and col 3.
// ============================================================================

interface CellInfo {
  z: number;
  symbol: string;
  name: string;
  gridCol: number; // 1-based column in our compact 9-col grid (col 3 is spacer)
  gridRow: number;
}

const LAYOUT: CellInfo[] = [
  // Row 1
  { z: 1,  symbol: 'H',  name: 'Hydrogen',   gridRow: 1, gridCol: 1 },
  { z: 2,  symbol: 'He', name: 'Helium',     gridRow: 1, gridCol: 9 },
  // Row 2 s-block
  { z: 3,  symbol: 'Li', name: 'Lithium',    gridRow: 2, gridCol: 1 },
  { z: 4,  symbol: 'Be', name: 'Beryllium',  gridRow: 2, gridCol: 2 },
  // Row 2 p-block
  { z: 5,  symbol: 'B',  name: 'Boron',      gridRow: 2, gridCol: 4 },
  { z: 6,  symbol: 'C',  name: 'Carbon',     gridRow: 2, gridCol: 5 },
  { z: 7,  symbol: 'N',  name: 'Nitrogen',   gridRow: 2, gridCol: 6 },
  { z: 8,  symbol: 'O',  name: 'Oxygen',     gridRow: 2, gridCol: 7 },
  { z: 9,  symbol: 'F',  name: 'Fluorine',   gridRow: 2, gridCol: 8 },
  { z: 10, symbol: 'Ne', name: 'Neon',       gridRow: 2, gridCol: 9 },
  // Row 3 s-block
  { z: 11, symbol: 'Na', name: 'Sodium',     gridRow: 3, gridCol: 1 },
  { z: 12, symbol: 'Mg', name: 'Magnesium',  gridRow: 3, gridCol: 2 },
  // Row 3 p-block
  { z: 13, symbol: 'Al', name: 'Aluminium',  gridRow: 3, gridCol: 4 },
  { z: 14, symbol: 'Si', name: 'Silicon',    gridRow: 3, gridCol: 5 },
  { z: 15, symbol: 'P',  name: 'Phosphorus', gridRow: 3, gridCol: 6 },
  { z: 16, symbol: 'S',  name: 'Sulfur',     gridRow: 3, gridCol: 7 },
  { z: 17, symbol: 'Cl', name: 'Chlorine',   gridRow: 3, gridCol: 8 },
  { z: 18, symbol: 'Ar', name: 'Argon',      gridRow: 3, gridCol: 9 },
];

// ============================================================================
// Component
// ============================================================================

export function ElementSelector() {
  const selectedElement = useBasisStore((state) => state.selectedElement);
  const setElement = useBasisStore((state) => state.setElement);
  const selectedInfo = ELEMENTS.find((el) => el.z === selectedElement);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      {/* Header with selected element badge */}
      <div className="flex items-center justify-between mb-3">
        <h2 className="text-base font-semibold text-slate-800">Select Element</h2>
        {selectedInfo && (
          <span className="inline-flex items-center gap-1.5 px-2.5 py-0.5 bg-primary-50 text-primary-700 rounded-full text-xs font-semibold">
            {selectedInfo.symbol}
            <span className="font-normal text-primary-500">{selectedInfo.name}</span>
          </span>
        )}
      </div>

      {/* Compact Periodic Table Grid: 9 columns (col 3 is narrow spacer) */}
      <div
        className="grid gap-1"
        style={{
          gridTemplateColumns: 'repeat(2, 1fr) 0.4fr repeat(6, 1fr)',
          gridTemplateRows: 'repeat(3, auto)',
        }}
        role="radiogroup"
        aria-label="Select an element"
      >
        {LAYOUT.map((cell) => {
          const isSelected = selectedElement === cell.z;
          const cat = getCategory(cell.z);
          const baseStyle = CAT_STYLE[cat];

          return (
            <button
              key={cell.z}
              type="button"
              role="radio"
              aria-checked={isSelected}
              aria-label={`${cell.name} (${cell.symbol}, Z=${cell.z})`}
              title={`${cell.name} (Z=${cell.z})`}
              onClick={() => setElement(cell.z)}
              className={`
                flex flex-col items-center justify-center rounded border
                transition-all duration-100 cursor-pointer
                focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-500
                h-9
                ${isSelected
                  ? 'bg-primary-600 text-white border-primary-700 shadow-md ring-2 ring-primary-300 font-bold'
                  : baseStyle
                }
              `}
              style={{ gridRow: cell.gridRow, gridColumn: cell.gridCol }}
            >
              <span className={`text-[0.55rem] leading-none ${isSelected ? 'text-primary-200' : 'opacity-40'}`}>
                {cell.z}
              </span>
              <span className="text-xs font-bold leading-tight">{cell.symbol}</span>
            </button>
          );
        })}
      </div>

      {/* Compact legend */}
      <div className="mt-2.5 pt-2 border-t border-slate-100 flex flex-wrap gap-x-3 gap-y-1">
        {(['nonmetal', 'halogen', 'noble', 'alkali', 'alkaline', 'metalloid'] as ElementCategory[]).map((cat) => (
          <div key={cat} className="flex items-center gap-1">
            <div className={`w-2 h-2 rounded-sm ${CAT_DOT[cat]}`} />
            <span className="text-[0.6rem] text-slate-400">{CAT_LABEL[cat]}</span>
          </div>
        ))}
      </div>
    </div>
  );
}
