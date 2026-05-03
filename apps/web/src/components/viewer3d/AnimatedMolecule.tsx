/**
 * AnimatedMolecule — 3D molecular viewer with normal-mode atom animation.
 *
 * Self-contained Canvas wrapper that renders sphere-and-stick atoms whose
 * positions oscillate along a normal mode coordinate. Used by the
 * Frequency tab's `NormalModeViewer` instead of the static `MoleculeViewer`
 * + `DisplacementArrows` overlay pattern (US-102 UX refinement).
 *
 * For each atom A and selected mode k, the position at time t is:
 *
 *   r_A(t) = R_A + amplitude · sin(2π · speed · t) · q^(k)_A
 *
 * where R_A is the equilibrium position (in bohr) and q^(k)_A is the
 * cartesian displacement vector for atom A in mode k. All atoms share
 * the same `sin(2π · speed · t)` factor — that is what gives a normal
 * mode its coordinated bending/stretching appearance.
 *
 * Bonds are recomputed every frame from the displaced atom positions, so
 * cylinders smoothly follow the atoms (bend, stretch, contract). Bond
 * connectivity is detected once at equilibrium and reused throughout the
 * animation, since at the small amplitudes in scope no bond can break.
 *
 * **Performance**: 60 fps target. Uses ref-based direct Three.js mutation
 * inside `useFrame`. No React re-renders happen during animation. The
 * scratch `Float32Array` for displaced positions is allocated once.
 *
 * **Accessibility**: When `prefers-reduced-motion` is set system-wide, the
 * caller (`NormalModeViewer`) is responsible for forcing `isAnimating=false`.
 * Static rendering at equilibrium is the default behavior in that case.
 *
 * **Code split**: This module imports Three.js and `@react-three/fiber`,
 * so it should be loaded via `React.lazy()` to keep the initial bundle
 * small. The existing `viewer3d/index.ts` barrel re-exports it for that
 * purpose.
 *
 * @module components/viewer3d/AnimatedMolecule
 * @see US-102 Frequency Tab UI (UX refinement)
 */

import { Suspense, useCallback, useEffect, useMemo, useRef } from 'react';
import { Canvas, useFrame, useThree } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import type { OrbitControls as OrbitControlsImpl } from 'three-stdlib';
import * as THREE from 'three';
import type { Atom3D } from './AtomSpheres';
import { AtomLabels } from './AtomLabels';
import { DisplacementArrows } from './DisplacementArrows';
import { detectBonds, type Bond } from './bondDetection';
import {
  CPK_COLORS,
  DEFAULT_ATOM_COLOR,
  COVALENT_RADII_BOHR,
  DEFAULT_COVALENT_RADIUS_BOHR,
  SPHERE_RADIUS_SCALE,
} from './constants';
import {
  computeDisplacedPositions,
  computeBondGeometry,
  type AnimAtom,
  type AnimDisplacement,
} from './animatedMoleculePositions';

// ============================================================================
// Constants — mirror AtomSpheres / BondCylinders for visual consistency
// ============================================================================

const SPHERE_WIDTH_SEGMENTS = 32;
const SPHERE_HEIGHT_SEGMENTS = 24;

/** Cylinder radius in bohr (matches BondCylinders). */
const BOND_RADIUS = 0.08;
/** Neutral gray (matches BondCylinders). */
const BOND_COLOR = '#888888';
const BOND_RADIAL_SEGMENTS = 8;

/** Three.js cylinders default to the +Y axis. */
const DEFAULT_CYLINDER_AXIS = new THREE.Vector3(0, 1, 0);

// ============================================================================
// Props
// ============================================================================

/**
 * Props for the AnimatedMolecule component.
 */
export interface AnimatedMoleculeProps {
  /**
   * Equilibrium atoms (in bohr) — same shape as `MoleculeViewer.atoms`.
   * Positions are the unperturbed nuclear coordinates from the geometry
   * optimization that produced the FrequencyResult.
   */
  atoms: Atom3D[];
  /**
   * Cartesian displacement vectors for the selected normal mode, one per
   * atom, in bohr. Length should equal `atoms.length`. If empty / shorter,
   * missing atoms are treated as having zero displacement (static).
   */
  displacement: ReadonlyArray<readonly [number, number, number]>;
  /** Peak oscillation amplitude in bohr (slider; 0.1–2.0 typical). */
  amplitude: number;
  /** Animation frequency in cycles per second (slider; 0.5–3.0 typical). */
  speed: number;
  /**
   * Whether the animation loop is running. When false, atoms remain at
   * equilibrium (the displacement is multiplied by phase=0).
   *
   * The caller is responsible for honoring `prefers-reduced-motion` —
   * typically by forcing this to false when the media query matches.
   */
  isAnimating: boolean;
  /**
   * Whether to overlay static displacement-vector arrows at the equilibrium
   * positions. The arrows are NOT pulsed; they show the direction and
   * magnitude of the mode's displacement field as a static visualization.
   */
  showArrows: boolean;
  /** Whether to render `AtomLabels` (element symbol + 1-based index). */
  showLabels: boolean;
  /** Optional ARIA label for the canvas wrapper. */
  ariaLabel?: string;
  /** Optional className for the outer container. */
  className?: string;
}

// ============================================================================
// Internal helpers
// ============================================================================

/**
 * Per-atom static descriptor (color + sphere radius).
 *
 * Computed once from the equilibrium atoms array and reused across frames
 * — only the world position is mutated by `useFrame`.
 */
interface AtomDescriptor {
  color: string;
  radius: number;
}

/** Compute the per-atom static descriptors. */
function makeAtomDescriptors(atoms: Atom3D[]): AtomDescriptor[] {
  return atoms.map((atom) => ({
    color: CPK_COLORS[atom.symbol] ?? DEFAULT_ATOM_COLOR,
    radius:
      (COVALENT_RADII_BOHR[atom.symbol] ?? DEFAULT_COVALENT_RADIUS_BOHR) *
      SPHERE_RADIUS_SCALE,
  }));
}

/**
 * Convert an `Atom3D[]` (the existing viewer3d shape, with bohr positions)
 * into the flat `[Z, x, y, z]` form expected by `computeDisplacedPositions`.
 *
 * The atomic number is set to 0 (unused; the helper only reads the position).
 */
function atomsToAnimAtoms(atoms: Atom3D[]): AnimAtom[] {
  return atoms.map(
    (a) => [0, a.position[0], a.position[1], a.position[2]] as AnimAtom
  );
}

/** Compute camera framing for an `Atom3D[]` (mirrors MoleculeViewer). */
function computeCameraFraming(atoms: Atom3D[]): {
  position: [number, number, number];
  target: [number, number, number];
} {
  if (atoms.length === 0) {
    return { position: [0, 0, 10], target: [0, 0, 0] };
  }
  let cx = 0;
  let cy = 0;
  let cz = 0;
  for (const a of atoms) {
    cx += a.position[0];
    cy += a.position[1];
    cz += a.position[2];
  }
  cx /= atoms.length;
  cy /= atoms.length;
  cz /= atoms.length;
  let maxDist = 0;
  for (const a of atoms) {
    const dx = a.position[0] - cx;
    const dy = a.position[1] - cy;
    const dz = a.position[2] - cz;
    const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
    if (d > maxDist) maxDist = d;
  }
  const cameraDistance = Math.max(maxDist * 2.5 + 2.0, 5.0);
  return {
    position: [cx, cy, cz + cameraDistance],
    target: [cx, cy, cz],
  };
}

function isWebGLAvailable(): boolean {
  try {
    const canvas = document.createElement('canvas');
    return !!(
      window.WebGLRenderingContext &&
      (canvas.getContext('webgl') || canvas.getContext('webgl2'))
    );
  } catch {
    return false;
  }
}

function checkPrefersReducedMotion(): boolean {
  if (typeof window === 'undefined') return false;
  return window.matchMedia('(prefers-reduced-motion: reduce)').matches;
}

function WebGLFallback() {
  return (
    <div className="flex items-center justify-center h-full min-h-[16rem] bg-slate-50 rounded-xl border border-slate-200">
      <div className="text-center p-6">
        <p className="text-slate-600 font-medium">
          3D visualization requires WebGL
        </p>
        <p className="text-sm text-slate-500 mt-1">
          Please use a WebGL-capable browser to view the molecular structure.
        </p>
      </div>
    </div>
  );
}

// ============================================================================
// Inner scene — runs inside the Canvas (uses R3F hooks)
// ============================================================================

interface AnimatedSceneProps {
  atoms: Atom3D[];
  descriptors: AtomDescriptor[];
  bonds: Bond[];
  animAtoms: AnimAtom[];
  displacement: ReadonlyArray<readonly [number, number, number]>;
  amplitude: number;
  speed: number;
  isAnimating: boolean;
  showArrows: boolean;
  showLabels: boolean;
  prefersReducedMotion: boolean;
  defaultCameraPosition: [number, number, number];
  defaultCameraTarget: [number, number, number];
}

/**
 * The animated scene contents.
 *
 * Mounts:
 *  - Lights (ambient + 2 directional)
 *  - One mesh per atom (with refs for per-frame position updates)
 *  - One mesh per bond (with refs for per-frame position/quaternion/scale)
 *  - AtomLabels (re-rendered when atoms move via React because Drei `<Html>`
 *    is positioned in scene coordinates)
 *  - Optional static DisplacementArrows overlay
 *  - OrbitControls (rotate / zoom / pan)
 *
 * The `useFrame` callback computes displaced positions and bond geometry
 * each frame and writes them directly into the mesh refs — no React state
 * updates, no re-renders. This is the standard R3F pattern for 60 fps
 * mutating animation.
 */
function AnimatedScene({
  atoms,
  descriptors,
  bonds,
  animAtoms,
  displacement,
  amplitude,
  speed,
  isAnimating,
  showArrows,
  showLabels,
  prefersReducedMotion,
  defaultCameraPosition,
  defaultCameraTarget,
}: AnimatedSceneProps): JSX.Element {
  const controlsRef = useRef<OrbitControlsImpl>(null);
  const { camera } = useThree();

  // ---- Refs to each atom mesh and each bond mesh -------------------------
  // Pre-sized arrays so the per-frame loop can index into them without
  // bounds checks.
  const atomMeshRefs = useRef<Array<THREE.Mesh | null>>([]);
  const bondMeshRefs = useRef<Array<THREE.Mesh | null>>([]);
  // Resize the ref arrays whenever the molecule shape changes.
  if (atomMeshRefs.current.length !== atoms.length) {
    atomMeshRefs.current = new Array(atoms.length).fill(null);
  }
  if (bondMeshRefs.current.length !== bonds.length) {
    bondMeshRefs.current = new Array(bonds.length).fill(null);
  }

  // ---- Scratch buffers, allocated once per molecule ----------------------
  const scratchPositions = useRef<Float32Array>(
    new Float32Array(atoms.length * 3)
  );
  if (scratchPositions.current.length !== atoms.length * 3) {
    scratchPositions.current = new Float32Array(atoms.length * 3);
  }
  // Reusable Three.js objects for the per-frame loop. Allocating these in
  // the closure (not inside useFrame) keeps GC pressure to zero.
  const tmpQuaternion = useRef(new THREE.Quaternion());
  const tmpDir = useRef(new THREE.Vector3());

  // ---- Bond index pairs (stable across the animation) -------------------
  const bondPairs = useMemo<Array<readonly [number, number]>>(
    () => bonds.map((b) => [b.atomIndexA, b.atomIndexB] as const),
    [bonds]
  );

  // ---- Per-frame mutation -----------------------------------------------
  useFrame(({ clock }) => {
    const phase = isAnimating
      ? Math.sin(2 * Math.PI * speed * clock.elapsedTime)
      : 0;

    // 1. Compute displaced positions directly into the pre-allocated scratch
    //    buffer — zero per-frame allocation for the position calculation.
    computeDisplacedPositions(
      animAtoms,
      displacement,
      amplitude,
      phase,
      scratchPositions.current
    );

    // 2. Update each atom mesh's position.
    for (let i = 0; i < atoms.length; i++) {
      const mesh = atomMeshRefs.current[i];
      if (!mesh) continue;
      mesh.position.set(
        scratchPositions.current[i * 3 + 0],
        scratchPositions.current[i * 3 + 1],
        scratchPositions.current[i * 3 + 2]
      );
    }

    // 3. Update each bond mesh's transform from the displaced positions.
    //    `computeBondGeometry` allocates its result array, but the array is
    //    short-lived (one record per bond, 1–3 bonds for typical small
    //    molecules in scope) and well within the per-frame GC budget.
    const geoms = computeBondGeometry(scratchPositions.current, bondPairs);
    for (let i = 0; i < geoms.length; i++) {
      const mesh = bondMeshRefs.current[i];
      if (!mesh) continue;
      const g = geoms[i];
      mesh.position.set(g.midpoint[0], g.midpoint[1], g.midpoint[2]);
      // Orient the cylinder along the bond direction.
      tmpDir.current.set(g.direction[0], g.direction[1], g.direction[2]);
      tmpQuaternion.current.setFromUnitVectors(
        DEFAULT_CYLINDER_AXIS,
        tmpDir.current
      );
      mesh.quaternion.copy(tmpQuaternion.current);
      // Scale the cylinder along its local Y axis (its long axis) so that
      // the rendered length tracks the actual bond length. The base
      // cylinder geometry has length 1; multiplying by `g.length` gives the
      // correct stretched/compressed visual.
      mesh.scale.set(1, g.length, 1);
    }
  });

  // ---- Reset camera helper (exposed via window for the HTML overlay button)
  const resetCamera = useCallback(() => {
    if (!controlsRef.current) return;
    camera.position.set(...defaultCameraPosition);
    controlsRef.current.target.set(...defaultCameraTarget);
    controlsRef.current.update();
  }, [camera, defaultCameraPosition, defaultCameraTarget]);

  // Reset camera when the molecule changes (e.g., user picks a new system).
  const isInitialMountRef = useRef(true);
  useEffect(() => {
    if (isInitialMountRef.current) {
      isInitialMountRef.current = false;
      return;
    }
    if (controlsRef.current) {
      camera.position.set(...defaultCameraPosition);
      controlsRef.current.target.set(...defaultCameraTarget);
      controlsRef.current.update();
    }
  }, [defaultCameraPosition, defaultCameraTarget, camera]);

  // Bridge to the HTML overlay reset button. Same pragmatic mechanism used by
  // MoleculeViewer (see SceneContent in MoleculeViewer.tsx).
  useMemo(() => {
    (window as unknown as Record<string, unknown>).__iqcp_reset_camera =
      resetCamera;
    return resetCamera;
  }, [resetCamera]);

  return (
    <>
      <ambientLight intensity={0.6} />
      <directionalLight position={[10, 10, 10]} intensity={0.8} />
      <directionalLight position={[-5, -5, 5]} intensity={0.3} />

      {/* Animated atom spheres -- one mesh per atom, position updated per
          frame via ref mutation. Initial position is the equilibrium so
          the first frame matches static rendering. */}
      <group>
        {atoms.map((atom, i) => (
          <mesh
            key={`atom-${i}`}
            ref={(el) => {
              atomMeshRefs.current[i] = el;
            }}
            position={atom.position}
            scale={descriptors[i]?.radius ?? 1}
          >
            <sphereGeometry
              args={[1, SPHERE_WIDTH_SEGMENTS, SPHERE_HEIGHT_SEGMENTS]}
            />
            <meshStandardMaterial color={descriptors[i]?.color ?? DEFAULT_ATOM_COLOR} />
          </mesh>
        ))}
      </group>

      {/* Animated bond cylinders -- one mesh per bond, transform updated per
          frame via ref mutation. Base cylinder has unit height; the per-frame
          scale.y mutation gives the rendered length. */}
      <group>
        {bonds.map((bond, i) => (
          <mesh
            key={`bond-${bond.atomIndexA}-${bond.atomIndexB}`}
            ref={(el) => {
              bondMeshRefs.current[i] = el;
            }}
            position={bond.midpoint}
            scale={[1, bond.length, 1]}
          >
            <cylinderGeometry
              args={[BOND_RADIUS, BOND_RADIUS, 1, BOND_RADIAL_SEGMENTS]}
            />
            <meshStandardMaterial color={BOND_COLOR} />
          </mesh>
        ))}
      </group>

      {/* Atom labels (no per-frame update; positioned at equilibrium). For the
          UX-refined view we accept that labels stay anchored to the equilibrium
          position rather than tracking the moving atom — that's the convention
          in GaussView and Avogadro and avoids extra HTML reflow per frame. */}
      <AtomLabels atoms={atoms} visible={showLabels} />

      {/* Optional static displacement-vector arrows. We pass `isAnimating=false`
          so the arrows do NOT pulse; they show the static direction field of
          the selected mode. */}
      {showArrows && displacement.length > 0 && (
        <DisplacementArrows
          atoms={atoms}
          displacement={displacement.map(
            (d) => [d[0], d[1], d[2]] as [number, number, number]
          )}
          amplitudeAng={amplitude}
          animationSpeed={speed}
          isAnimating={false}
        />
      )}

      <OrbitControls
        ref={controlsRef}
        makeDefault
        enableDamping={!prefersReducedMotion}
        dampingFactor={0.1}
      />
    </>
  );
}

// ============================================================================
// Public component
// ============================================================================

/**
 * AnimatedMolecule — top-level Canvas wrapper for normal-mode animation.
 *
 * Renders a self-contained 3D scene with atoms that oscillate along the
 * selected normal mode. Drop-in replacement for `MoleculeViewer +
 * DisplacementArrows` overlay in the Frequency tab.
 *
 * **Lazy load this component** — it pulls in Three.js and React Three Fiber.
 *
 * @see US-102 Frequency Tab UI (UX refinement)
 */
export function AnimatedMolecule({
  atoms,
  displacement,
  amplitude,
  speed,
  isAnimating,
  showArrows,
  showLabels,
  ariaLabel,
  className = '',
}: AnimatedMoleculeProps): JSX.Element {
  const webGLAvailable = useMemo(() => isWebGLAvailable(), []);
  const prefersReducedMotion = useMemo(() => checkPrefersReducedMotion(), []);

  // Camera framing (recomputed when the molecule changes shape).
  const { position: cameraPosition, target: cameraTarget } = useMemo(
    () => computeCameraFraming(atoms),
    [atoms]
  );

  // Stable per-atom descriptors (color + radius).
  const descriptors = useMemo(() => makeAtomDescriptors(atoms), [atoms]);

  // Bonds detected once at equilibrium. Connectivity is reused across the
  // animation (see Bond docstring on AnimatedMolecule).
  const bonds = useMemo(() => detectBonds(atoms), [atoms]);

  // Reshape atoms into the AnimAtom format expected by the math helpers.
  const animAtoms = useMemo(() => atomsToAnimAtoms(atoms), [atoms]);

  // Coerce the displacement to a typed array form once per render. We do
  // not normalize: the magnitude of the user-supplied amplitude is the
  // single source of truth for visual scale. (Normalization can be added
  // here if needed but for the typical FrequencyResult the cartesian
  // displacements are already O(0.01–0.1) bohr and respond well to a
  // 0.5 bohr default amplitude.)
  const animDisp = useMemo<AnimDisplacement[]>(
    () =>
      displacement.map(
        (d) => [d[0], d[1], d[2]] as AnimDisplacement
      ),
    [displacement]
  );

  // ARIA label
  const effectiveAriaLabel = useMemo(() => {
    if (ariaLabel) return ariaLabel;
    if (atoms.length === 0)
      return '3D molecular viewer (no molecule loaded)';
    const symbolCounts: Record<string, number> = {};
    for (const a of atoms) {
      symbolCounts[a.symbol] = (symbolCounts[a.symbol] ?? 0) + 1;
    }
    const formula = Object.entries(symbolCounts)
      .map(([sym, count]) => (count > 1 ? `${sym}${count}` : sym))
      .join('');
    return `3D animated view of molecule ${formula}`;
  }, [atoms, ariaLabel]);

  // Reset camera button (HTML overlay).
  const handleResetCamera = useCallback(() => {
    const fn = (window as unknown as Record<string, unknown>)
      .__iqcp_reset_camera;
    if (typeof fn === 'function') (fn as () => void)();
  }, []);

  if (!webGLAvailable) {
    return <WebGLFallback />;
  }

  return (
    <div
      className={`relative w-full h-full min-h-[16rem] ${className}`}
      role="img"
      aria-label={effectiveAriaLabel}
    >
      <Canvas
        camera={{
          position: cameraPosition,
          fov: 50,
          near: 0.1,
          far: 1000,
        }}
        gl={{ preserveDrawingBuffer: true }}
        style={{ background: '#f8fafc' }}
      >
        <Suspense fallback={null}>
          <AnimatedScene
            atoms={atoms}
            descriptors={descriptors}
            bonds={bonds}
            animAtoms={animAtoms}
            displacement={animDisp}
            amplitude={amplitude}
            speed={speed}
            isAnimating={isAnimating}
            showArrows={showArrows}
            showLabels={showLabels}
            prefersReducedMotion={prefersReducedMotion}
            defaultCameraPosition={cameraPosition}
            defaultCameraTarget={cameraTarget}
          />
        </Suspense>
      </Canvas>

      {/* Reset camera button (HTML overlay) */}
      <button
        type="button"
        onClick={handleResetCamera}
        className="absolute top-2 right-2 px-2 py-1 text-xs font-medium text-slate-600 bg-white/80 hover:bg-white border border-slate-200 rounded shadow-sm backdrop-blur-sm transition-colors focus:outline-none focus:ring-2 focus:ring-blue-400 focus:ring-offset-1"
        aria-label="Reset camera to default view"
        title="Reset camera to default view"
      >
        Reset View
      </button>
    </div>
  );
}

export default AnimatedMolecule;
