/**
 * AtomLabels - Billboard text labels for molecular atoms.
 *
 * Renders element symbol + 1-based index (e.g., "O1", "H2", "H3") at each
 * atom position using Drei's `<Html>` component for automatic billboard
 * behavior (labels always face the camera).
 *
 * Labels are designed to be togglable via the `visible` prop without
 * causing re-renders of sibling components (AtomSpheres, BondCylinders).
 *
 * @module viewer3d/AtomLabels
 */

import { useMemo } from 'react';
import { Html } from '@react-three/drei';
import type { Atom3D } from './AtomSpheres';

/**
 * Props for the AtomLabels component.
 */
export interface AtomLabelsProps {
  /** Array of atoms to label */
  atoms: Atom3D[];
  /** Whether labels are visible. When false, nothing is rendered. */
  visible: boolean;
}

/**
 * Precomputed label data for a single atom.
 */
interface LabelData {
  /** Display text (e.g., "O1", "H2") */
  text: string;
  /** 3D position [x, y, z] in bohr */
  position: [number, number, number];
  /** Unique key for React reconciliation */
  key: string;
}

/**
 * Renders billboard text labels for each atom in the molecule.
 *
 * Uses Drei's `<Html>` component which projects 3D positions onto the
 * screen and renders standard HTML elements that always face the camera.
 * This approach provides:
 * - Automatic billboard behavior (no quaternion math needed)
 * - Standard CSS styling for text
 * - Size attenuation via `distanceFactor` (labels shrink with distance)
 * - Proper depth occlusion ordering via `zIndexRange`
 *
 * Performance: Label data is memoized on the atoms array. When `visible`
 * is false, the component returns null immediately, producing no DOM
 * elements and no Three.js objects.
 */
export function AtomLabels({ atoms, visible }: AtomLabelsProps) {
  // Precompute label data to avoid recalculation on re-renders
  const labels: LabelData[] = useMemo(() => {
    return atoms.map((atom, index) => ({
      text: `${atom.symbol}${index + 1}`,
      position: atom.position,
      key: `label-${index}-${atom.symbol}`,
    }));
  }, [atoms]);

  // When not visible, render nothing -- no DOM elements, no overhead
  if (!visible) return null;
  if (labels.length === 0) return null;

  return (
    <group>
      {labels.map((label) => (
        <Html
          key={label.key}
          position={label.position}
          center
          distanceFactor={8}
          zIndexRange={[1, 0]}
          sprite
          pointerEvents="none"
        >
          <div
            style={{
              background: 'rgba(30, 41, 59, 0.75)',
              color: '#f1f5f9',
              padding: '1px 4px',
              borderRadius: '3px',
              fontSize: '11px',
              fontFamily:
                'ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace',
              fontWeight: 600,
              letterSpacing: '0.02em',
              whiteSpace: 'nowrap',
              userSelect: 'none',
              lineHeight: '1.3',
            }}
          >
            {label.text}
          </div>
        </Html>
      ))}
    </group>
  );
}
