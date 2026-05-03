/**
 * Constants for 3D molecular visualization.
 *
 * CPK colors and covalent radii for elements H through Ar (Z = 1..18).
 * These cover all elements present in the current preset systems
 * and supported basis sets.
 *
 * @module viewer3d/constants
 */

/**
 * CPK (Corey-Pauling-Koltun) color palette for elements H-Ar.
 *
 * Standard CPK colors used in molecular visualization software
 * (Jmol, PyMOL, VMD). Values from TDD Section 9.1 (Phase 2).
 * Keys are element symbols (uppercase first letter).
 */
export const CPK_COLORS: Record<string, string> = {
  H: '#FFFFFF',
  He: '#D9FFFF',
  Li: '#CC80FF',
  Be: '#C2FF00',
  B: '#FFB5B5',
  C: '#909090',
  N: '#3050F8',
  O: '#FF0D0D',
  F: '#90E050',
  Ne: '#B3E3F5',
  Na: '#AB5CF2',
  Mg: '#8AFF00',
  Al: '#BFA6A6',
  Si: '#F0C8A0',
  P: '#FF8000',
  S: '#FFFF30',
  Cl: '#1FF01F',
  Ar: '#80D1E3',
};

/** Fallback color for unrecognized elements (hot pink -- visually obvious). */
export const DEFAULT_ATOM_COLOR = '#FF69B4';

/**
 * Covalent radii in bohr for elements H-Ar.
 *
 * Used for:
 * 1. Visual sphere sizing (scaled by SPHERE_RADIUS_SCALE)
 * 2. Bond detection in US-033 (distance < BOND_TOLERANCE * (r_A + r_B))
 *
 * Values from TDD Section 8.4 (Phase 2) and Pyykko & Atsumi (2009).
 * Conversion: r_bohr = r_angstrom / 0.529177
 */
export const COVALENT_RADII_BOHR: Record<string, number> = {
  H: 0.586,
  // He radius increased from 0.529 (standard covalent estimate) to 0.700 bohr
  // to ensure HeH+ bond detection succeeds. The standard value of 0.28 Angstrom
  // is an estimate for hypothetical covalent He; HeH+ is ionic/charge-transfer
  // with a bond length (1.4632 bohr) exceeding the threshold at 0.529.
  // 0.700 bohr (~0.37 Angstrom) is within the range of published estimates
  // (Pyykko 2009 gives 0.46 Angstrom = 0.869 bohr).
  He: 0.700,
  Li: 2.419,
  Be: 1.814,
  B: 1.587,
  C: 1.455,
  N: 1.361,
  O: 1.304,
  F: 1.285,
  Ne: 1.247,
  Na: 3.023,
  Mg: 2.627,
  Al: 2.362,
  Si: 2.192,
  P: 2.098,
  S: 1.984,
  Cl: 1.928,
  Ar: 1.852,
};

/** Default covalent radius for unrecognized elements (bohr). */
export const DEFAULT_COVALENT_RADIUS_BOHR = 1.0;

/**
 * Visual scaling factor for sphere radii.
 *
 * The raw covalent radii in bohr produce spheres that overlap heavily
 * in a ball-and-stick model. This factor is multiplied by the covalent
 * radius to produce the rendered sphere radius.
 *
 * Tuned so that:
 * - H atoms are visible but smaller than heavy atoms
 * - Bonds remain visible between spheres (not obscured)
 * - Size differences between elements are clearly apparent
 *
 * Example rendered radii:
 *   H:  0.586 * 0.35 = 0.205 bohr
 *   O:  1.304 * 0.35 = 0.456 bohr
 *   Li: 2.419 * 0.35 = 0.847 bohr
 */
export const SPHERE_RADIUS_SCALE = 0.35;

/**
 * Bond detection tolerance factor.
 *
 * Two atoms A, B are considered bonded if:
 *   distance(A, B) < BOND_TOLERANCE * (r_cov(A) + r_cov(B))
 *
 * The 1.2 factor accounts for bonds slightly longer than ideal
 * covalent distances while avoiding false positives.
 *
 * TDD Section 8.4.
 */
export const BOND_TOLERANCE = 1.2;

/**
 * Bohr to Angstrom conversion factor.
 *
 * 1 bohr = 0.529177 Angstrom (CODATA 2018).
 */
export const BOHR_TO_ANGSTROM = 0.529177;
