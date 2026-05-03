/**
 * BondCylinders - Cylinder rendering for molecular bonds.
 *
 * Renders detected bonds as gray cylinders connecting bonded atom pairs.
 * Bond detection is performed internally via detectBonds(), so this
 * component accepts the same atoms prop as AtomSpheres.
 *
 * Each bond is rendered as an individual <mesh> with a CylinderGeometry.
 * This is simpler than instancedMesh (which requires per-instance quaternion
 * setup) and is efficient for our molecule sizes (1-3 bonds max in presets).
 *
 * @module viewer3d/BondCylinders
 */

import { useMemo } from 'react';
import * as THREE from 'three';
import { detectBonds } from './bondDetection';
import type { Atom3D } from './AtomSpheres';

/**
 * Props for the BondCylinders component.
 */
export interface BondCylindersProps {
  /** Atoms to detect and render bonds for (positions in bohr) */
  atoms: Atom3D[];
}

/** Cylinder radius in bohr -- thin enough to not obscure atoms, thick enough to be visible. */
const BOND_RADIUS = 0.08;

/** Neutral gray -- does not compete with CPK atom colors for visual attention. */
const BOND_COLOR = '#888888';

/** Radial segments for cylinder smoothness (8 produces a visually round cylinder). */
const BOND_RADIAL_SEGMENTS = 8;

/** The default cylinder axis in Three.js is the Y-axis. */
const DEFAULT_AXIS = new THREE.Vector3(0, 1, 0);

/**
 * Renders bonds between atoms as cylinder meshes.
 *
 * Internally calls detectBonds() on the atoms array and renders
 * a cylinder for each detected bond. Cylinders are positioned at
 * bond midpoints and rotated to align with bond direction vectors.
 *
 * The quaternion rotation approach (setFromUnitVectors) is robust for
 * all bond orientations including bonds along Y-axis (identity) and
 * bonds along negative Y-axis (180-degree rotation).
 */
export function BondCylinders({ atoms }: BondCylindersProps) {
  const bonds = useMemo(() => detectBonds(atoms), [atoms]);

  // Precompute quaternions for each bond (avoids recomputation on every render frame)
  const bondQuaternions = useMemo(() => {
    return bonds.map((bond) => {
      const bondDir = new THREE.Vector3(...bond.direction);
      const quaternion = new THREE.Quaternion();
      quaternion.setFromUnitVectors(DEFAULT_AXIS, bondDir);
      return quaternion;
    });
  }, [bonds]);

  if (bonds.length === 0) return null;

  return (
    <group>
      {bonds.map((bond, i) => (
        <mesh
          key={`bond-${bond.atomIndexA}-${bond.atomIndexB}`}
          position={bond.midpoint}
          quaternion={bondQuaternions[i]}
        >
          <cylinderGeometry
            args={[BOND_RADIUS, BOND_RADIUS, bond.length, BOND_RADIAL_SEGMENTS]}
          />
          <meshStandardMaterial color={BOND_COLOR} />
        </mesh>
      ))}
    </group>
  );
}
