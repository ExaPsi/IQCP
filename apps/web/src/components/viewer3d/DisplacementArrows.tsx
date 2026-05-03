/**
 * DisplacementArrows — R3F overlay for animated vibrational displacement arrows.
 *
 * Rendered as a child of `<MoleculeViewer>` (i.e., inside an R3F `<Canvas>`).
 * Draws one orange arrow per atom showing the direction of the selected
 * normal mode's cartesian displacement. When `isAnimating` is true, the
 * endpoint oscillates via `useFrame` with phase `sin(2π · speed · t)`, so
 * each atom sweeps back and forth along its displacement axis. When
 * `isAnimating` is false, the arrow is frozen at `phase = 1` (maximum
 * displacement).
 *
 * Bohr → Ångström conversion: qc-core returns atom positions in bohr and
 * normal mode coefficients in mass-weighted bohr · amu^(-1/2). This component
 * converts atom positions to Ångstroms for display (matching the coordinate
 * system the existing viewer3d components use) and treats the normalized
 * displacement as dimensionless unit vectors scaled by `amplitudeBohr`.
 *
 * @module components/viewer3d/DisplacementArrows
 * @see US-102 Frequency Tab UI, AC3
 */

import { useMemo, useRef } from 'react';
import { useFrame } from '@react-three/fiber';
import { Line } from '@react-three/drei';
import * as THREE from 'three';
import type { Atom3D } from './AtomSpheres';
import { BOHR_TO_ANGSTROM } from './constants';
import {
  normalizeDisplacement,
  animationPhase,
  displacementEndpoint,
  type DisplacementVec,
} from './displacementMath';

// ============================================================================
// Props
// ============================================================================

export interface DisplacementArrowsProps {
  /**
   * Atoms in the scene. Positions are in bohr (standard viewer3d convention),
   * and are converted to Å internally for rendering.
   */
  atoms: Atom3D[];
  /**
   * Cartesian displacement vectors, one per atom. May be any normalization
   * (qc-core returns mass-weighted normals); the component rescales to
   * unit L2 norm internally.
   */
  displacement: DisplacementVec[];
  /** Displacement amplitude in Ångströms (user slider). */
  amplitudeAng: number;
  /** Animation speed multiplier (user slider). */
  animationSpeed: number;
  /** Whether animation is playing. */
  isAnimating: boolean;
}

// ============================================================================
// Component
// ============================================================================

/**
 * R3F component that renders one orange displacement arrow per atom.
 *
 * Must be used inside an R3F `<Canvas>`.
 *
 * @see US-102 Frequency Tab UI, AC3
 */
export function DisplacementArrows({
  atoms,
  displacement,
  amplitudeAng,
  animationSpeed,
  isAnimating,
}: DisplacementArrowsProps): JSX.Element {
  // One ref per atom so we can update endpoint geometry without reallocation.
  // (drei `<Line>` forwards to a Line2 instance with a geometry.setFromPoints
  // method.)
  const lineRefs = useRef<Array<THREE.Object3D | null>>([]);

  // Normalize displacement to unit L2 norm so slider amplitude = physical Å.
  const normalized = useMemo(
    () => normalizeDisplacement(displacement),
    [displacement]
  );

  // Convert atom positions to Angstroms once per render.
  const startsAng = useMemo(
    () =>
      atoms.map((a) => [
        a.position[0] * BOHR_TO_ANGSTROM,
        a.position[1] * BOHR_TO_ANGSTROM,
        a.position[2] * BOHR_TO_ANGSTROM,
      ]) as Array<[number, number, number]>,
    [atoms]
  );

  // Initial endpoints (phase = 1 when static, phase = 0 at t=0 when animating).
  const initialEndpoints = useMemo(() => {
    const staticPhase = isAnimating ? 0 : 1;
    return startsAng.map((start, i) => {
      const d = normalized[i] ?? [0, 0, 0];
      return displacementEndpoint(start, d, amplitudeAng, staticPhase);
    });
  }, [startsAng, normalized, amplitudeAng, isAnimating]);

  // Animation loop: update each line's geometry from (start, end) without
  // creating new buffers. Runs at the R3F frame rate.
  useFrame(({ clock }) => {
    if (!isAnimating) return;
    const phase = animationPhase(clock.getElapsedTime(), animationSpeed);
    for (let i = 0; i < startsAng.length; i++) {
      const ref = lineRefs.current[i];
      if (!ref) continue;
      // drei Line2 wraps a LineGeometry with a setPositions method.
      const geom =
        (ref as unknown as { geometry?: { setPositions?: (arr: number[]) => void } })
          .geometry;
      if (!geom || typeof geom.setPositions !== 'function') continue;
      const start = startsAng[i];
      const d = normalized[i] ?? [0, 0, 0];
      const end = displacementEndpoint(start, d, amplitudeAng, phase);
      geom.setPositions([start[0], start[1], start[2], end[0], end[1], end[2]]);
    }
  });

  return (
    <group>
      {startsAng.map((start, i) => {
        const end = initialEndpoints[i] ?? start;
        return (
          <Line
            key={i}
            ref={(el) => {
              lineRefs.current[i] = el as unknown as THREE.Object3D | null;
            }}
            points={[start, end]}
            color="#f97316"
            lineWidth={3}
          />
        );
      })}
    </group>
  );
}

export default DisplacementArrows;
