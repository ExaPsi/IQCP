/**
 * AtomSpheres - Sphere rendering for molecular atoms.
 *
 * Renders atoms as CPK-colored spheres using individual Three.js meshes.
 * Each atom's color follows the CPK convention and its radius is
 * proportional to the element's covalent radius.
 *
 * Uses individual <mesh> components (same pattern as BondCylinders)
 * rather than instancedMesh, which has unreliable per-instance color
 * support in Three.js r183 + React Three Fiber 8.x due to shader
 * compilation timing. For IQCP's small molecules (2-10 atoms), individual
 * meshes have negligible performance overhead.
 *
 * @module viewer3d/AtomSpheres
 */

import { useMemo } from 'react';
import {
  CPK_COLORS,
  DEFAULT_ATOM_COLOR,
  COVALENT_RADII_BOHR,
  DEFAULT_COVALENT_RADIUS_BOHR,
  SPHERE_RADIUS_SCALE,
} from './constants';

/**
 * Atom data for 3D rendering.
 *
 * A simplified representation containing only the information
 * needed for visualization: element identity and position.
 * Positions are always in bohr (atomic units).
 */
export interface Atom3D {
  /** Element symbol (e.g., "H", "O", "Li") */
  symbol: string;
  /** Position [x, y, z] in bohr */
  position: [number, number, number];
  /**
   * Mulliken partial charge on this atom (in electron units).
   *
   * When present, overrides the default CPK color with a diverging
   * red (negative/electron-rich) to white (neutral) to blue (positive/electron-poor)
   * color scale. Charges are clamped to [-0.5, 0.5] for color mapping.
   */
  charge?: number;
}

/**
 * Props for the AtomSpheres component.
 */
export interface AtomSpheresProps {
  /** Array of atoms to render as spheres */
  atoms: Atom3D[];
}

/**
 * Convert a Mulliken partial charge to a diverging color string.
 *
 * Maps charge to a red-white-blue diverging scale:
 * - Negative charge (electron-rich): red
 * - Zero charge (neutral): white
 * - Positive charge (electron-poor): blue
 *
 * Charges are clamped to [-0.5, 0.5] for a consistent visual range.
 *
 * @param charge - Mulliken partial charge in electron units
 * @returns CSS color string (e.g., "rgb(255,128,128)")
 */
function chargeToColor(charge: number): string {
  const clamped = Math.max(-0.5, Math.min(0.5, charge));
  // t: 0 = most negative (red), 0.5 = neutral (white), 1 = most positive (blue)
  const t = (clamped + 0.5);
  if (t < 0.5) {
    // Red to white: R stays at 255, G and B increase from 0 to 255
    const f = t * 2; // 0 to 1
    const r = 1.0;
    const g = f;
    const b = f;
    return `rgb(${Math.round(r * 255)},${Math.round(g * 255)},${Math.round(b * 255)})`;
  } else {
    // White to blue: B stays at 255, R and G decrease from 255 to 0
    const f = (1 - t) * 2; // 1 to 0
    const r = f;
    const g = f;
    const b = 1.0;
    return `rgb(${Math.round(r * 255)},${Math.round(g * 255)},${Math.round(b * 255)})`;
  }
}

/** Sphere geometry segments for smooth appearance. */
const SPHERE_WIDTH_SEGMENTS = 32;
const SPHERE_HEIGHT_SEGMENTS = 24;

/**
 * Renders atoms as CPK-colored spheres.
 *
 * Each atom is rendered as an individual <mesh> with a shared sphere
 * geometry (via R3F caching) and per-atom meshStandardMaterial color.
 * Sphere radius is proportional to the element's covalent radius
 * scaled by SPHERE_RADIUS_SCALE.
 */
export function AtomSpheres({ atoms }: AtomSpheresProps) {
  const atomData = useMemo(
    () =>
      atoms.map((atom) => ({
        key: `${atom.symbol}-${atom.position[0]}-${atom.position[1]}-${atom.position[2]}`,
        position: atom.position as [number, number, number],
        radius:
          (COVALENT_RADII_BOHR[atom.symbol] ?? DEFAULT_COVALENT_RADIUS_BOHR) *
          SPHERE_RADIUS_SCALE,
        color: atom.charge !== undefined
          ? chargeToColor(atom.charge)
          : (CPK_COLORS[atom.symbol] ?? DEFAULT_ATOM_COLOR),
      })),
    [atoms]
  );

  if (atoms.length === 0) return null;

  return (
    <group>
      {atomData.map((atom) => (
        <mesh key={atom.key} position={atom.position} scale={atom.radius}>
          <sphereGeometry args={[1, SPHERE_WIDTH_SEGMENTS, SPHERE_HEIGHT_SEGMENTS]} />
          <meshStandardMaterial color={atom.color} />
        </mesh>
      ))}
    </group>
  );
}
