/**
 * GhostAtoms - Translucent sphere rendering for before/after geometry comparison.
 *
 * Renders atoms as translucent gray spheres to show the initial (pre-optimization)
 * geometry overlaid on the current/optimized geometry. Used for the ghost overlay
 * in trajectory animation (US-075, AC4).
 *
 * @module viewer3d/GhostAtoms
 */

import { useMemo } from 'react';
import {
  COVALENT_RADII_BOHR,
  DEFAULT_COVALENT_RADIUS_BOHR,
  SPHERE_RADIUS_SCALE,
} from './constants';
import type { Atom3D } from './AtomSpheres';

/**
 * Props for GhostAtoms.
 */
export interface GhostAtomsProps {
  /** Atoms to render as ghost spheres */
  atoms: Atom3D[];
  /** Opacity for the ghost spheres (default: 0.25) */
  opacity?: number;
}

/** Sphere geometry segments. */
const SPHERE_WIDTH_SEGMENTS = 24;
const SPHERE_HEIGHT_SEGMENTS = 16;

/**
 * Renders atoms as translucent gray spheres for geometry comparison.
 *
 * Used to show the initial geometry as a "ghost" overlay while viewing
 * the optimized geometry at full opacity.
 */
export function GhostAtoms({ atoms, opacity = 0.25 }: GhostAtomsProps) {
  const atomData = useMemo(
    () =>
      atoms.map((atom, i) => ({
        key: `ghost-${i}-${atom.symbol}`,
        position: atom.position as [number, number, number],
        radius:
          (COVALENT_RADII_BOHR[atom.symbol] ?? DEFAULT_COVALENT_RADIUS_BOHR) *
          SPHERE_RADIUS_SCALE,
      })),
    [atoms]
  );

  if (atoms.length === 0) return null;

  return (
    <group>
      {atomData.map((atom) => (
        <mesh key={atom.key} position={atom.position} scale={atom.radius}>
          <sphereGeometry args={[1, SPHERE_WIDTH_SEGMENTS, SPHERE_HEIGHT_SEGMENTS]} />
          <meshStandardMaterial
            color="#94a3b8"
            transparent
            opacity={opacity}
            depthWrite={false}
          />
        </mesh>
      ))}
    </group>
  );
}
