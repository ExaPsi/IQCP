/**
 * MoleculeViewer - R3F Canvas wrapper for 3D molecular visualization.
 *
 * Provides the top-level Three.js Canvas with camera, lighting, OrbitControls,
 * and child components for rendering atoms, bonds, and labels.
 *
 * This component is designed to be loaded via React.lazy() to ensure
 * Three.js is code-split into a separate chunk that is only fetched
 * when the 3D viewer is actually rendered.
 *
 * Features (US-032, US-033, US-034):
 * - Instanced mesh rendering for performance (AtomSpheres)
 * - Bond detection and cylinder rendering (BondCylinders)
 * - Togglable billboard atom labels (AtomLabels)
 * - OrbitControls: click-drag rotate, scroll zoom, right-drag pan
 * - Reset camera button to restore default view
 * - Camera state tracking via onCameraChange callback
 * - ARIA accessibility (role="img" + descriptive aria-label)
 * - Respects prefers-reduced-motion media query
 *
 * @module viewer3d/MoleculeViewer
 */

import { useCallback, useEffect, useMemo, useRef } from 'react';
import { Canvas, useThree } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import type { OrbitControls as OrbitControlsImpl } from 'three-stdlib';
import { AtomSpheres } from './AtomSpheres';
import { BondCylinders } from './BondCylinders';
import { AtomLabels } from './AtomLabels';
import type { Atom3D } from './AtomSpheres';

/**
 * Camera state for deep link serialization.
 *
 * Represents the camera orientation and framing of the 3D scene.
 * These values are approximate -- not pixel-identical across displays
 * or aspect ratios, but sufficient for restoring a similar view.
 */
export interface CameraState {
  /** Point the camera looks at [x, y, z] in scene coordinates */
  target: [number, number, number];
  /** Camera position [x, y, z] in scene coordinates */
  position: [number, number, number];
  /** Distance from camera to target */
  distance: number;
  /** Azimuthal angle in radians (rotation around Y axis) */
  azimuth: number;
  /** Polar angle in radians (rotation from Y axis, 0 = top, PI = bottom) */
  polar: number;
}

/**
 * Props for the MoleculeViewer component.
 */
export interface MoleculeViewerProps {
  /** Atoms to render (positions in bohr) */
  atoms: Atom3D[];
  /** Descriptive label for ARIA accessibility */
  ariaLabel?: string;
  /** Additional CSS classes for the container */
  className?: string;
  /** Whether to show atom labels (billboard text at each atom position) */
  showLabels?: boolean;
  /** Callback fired when camera position/orientation changes */
  onCameraChange?: (state: CameraState) => void;
  /** Additional R3F children to render inside the 3D scene (e.g., OrbitalSurface) */
  children?: React.ReactNode;
}

/**
 * Compute camera position to frame all atoms.
 *
 * Calculates the centroid and bounding sphere of the atom positions,
 * then places the camera at a distance that ensures all atoms are
 * visible with some padding.
 */
function computeCameraPosition(atoms: Atom3D[]): {
  position: [number, number, number];
  target: [number, number, number];
} {
  if (atoms.length === 0) {
    return { position: [0, 0, 10], target: [0, 0, 0] };
  }

  // Compute centroid
  let cx = 0,
    cy = 0,
    cz = 0;
  for (const atom of atoms) {
    cx += atom.position[0];
    cy += atom.position[1];
    cz += atom.position[2];
  }
  cx /= atoms.length;
  cy /= atoms.length;
  cz /= atoms.length;

  // Compute bounding sphere radius from centroid
  let maxDist = 0;
  for (const atom of atoms) {
    const dx = atom.position[0] - cx;
    const dy = atom.position[1] - cy;
    const dz = atom.position[2] - cz;
    const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
    if (dist > maxDist) maxDist = dist;
  }

  // Camera distance: bounding sphere radius + padding
  // Minimum distance of 5 bohr to handle single-atom or linear molecules
  const cameraDistance = Math.max(maxDist * 2.5 + 2.0, 5.0);

  return {
    position: [cx, cy, cz + cameraDistance],
    target: [cx, cy, cz],
  };
}

/**
 * Check if WebGL is available in the current browser.
 */
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

/**
 * Detect whether the user prefers reduced motion.
 *
 * When true, OrbitControls damping (smooth deceleration) is disabled
 * to respect the user's accessibility preference.
 */
function checkPrefersReducedMotion(): boolean {
  if (typeof window === 'undefined') return false;
  return window.matchMedia('(prefers-reduced-motion: reduce)').matches;
}

/**
 * Fallback component shown when WebGL is not available.
 */
function WebGLFallback() {
  return (
    <div className="flex items-center justify-center h-full min-h-[16rem] bg-slate-50 rounded-xl border border-slate-200">
      <div className="text-center p-6">
        <svg
          className="w-12 h-12 text-slate-400 mx-auto mb-3"
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
          aria-hidden="true"
        >
          <path
            strokeLinecap="round"
            strokeLinejoin="round"
            strokeWidth={1.5}
            d="M20 7l-8-4-8 4m16 0l-8 4m8-4v10l-8 4m0-10L4 7m8 4v10M4 7v10l8 4"
          />
        </svg>
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

/**
 * Inner scene component that has access to the R3F context (useThree).
 *
 * This is separated from the outer MoleculeViewer because R3F hooks
 * (useThree, useFrame) can only be called inside a <Canvas> tree.
 * The OrbitControls ref and camera reset logic live here.
 */
function SceneContent({
  atoms,
  showLabels,
  defaultCameraPosition,
  defaultCameraTarget,
  prefersReducedMotion,
  onCameraChange,
  children,
}: {
  atoms: Atom3D[];
  showLabels: boolean;
  defaultCameraPosition: [number, number, number];
  defaultCameraTarget: [number, number, number];
  prefersReducedMotion: boolean;
  onCameraChange?: (state: CameraState) => void;
  children?: React.ReactNode;
}) {
  const controlsRef = useRef<OrbitControlsImpl>(null);
  const { camera } = useThree();

  // Store the onCameraChange in a ref to avoid re-creating the onChange handler
  const onCameraChangeRef = useRef(onCameraChange);
  onCameraChangeRef.current = onCameraChange;

  // Store default positions for reset
  const defaultsRef = useRef({
    position: defaultCameraPosition,
    target: defaultCameraTarget,
  });
  defaultsRef.current = {
    position: defaultCameraPosition,
    target: defaultCameraTarget,
  };

  /**
   * Extract camera state from OrbitControls on change.
   *
   * Computes azimuthal and polar angles from the camera position
   * relative to the controls target, matching Three.js spherical
   * coordinate conventions.
   */
  const handleCameraChange = useCallback(() => {
    const cb = onCameraChangeRef.current;
    if (!cb || !controlsRef.current) return;

    const controls = controlsRef.current;
    const target = controls.target;

    // Distance from camera to target
    const distance = camera.position.distanceTo(target);

    // Compute spherical coordinates (relative to target)
    const dx = camera.position.x - target.x;
    const dy = camera.position.y - target.y;
    const dz = camera.position.z - target.z;

    // Polar angle: angle from Y axis (0 = top, PI = bottom)
    const polar = distance > 0 ? Math.acos(Math.min(1, Math.max(-1, dy / distance))) : 0;

    // Azimuthal angle: angle in XZ plane from Z axis
    const azimuth = Math.atan2(dx, dz);

    cb({
      target: [target.x, target.y, target.z],
      position: [camera.position.x, camera.position.y, camera.position.z],
      distance,
      azimuth,
      polar,
    });
  }, [camera]);

  /**
   * Reset camera to the default auto-computed position.
   *
   * This is called from the parent via the exposed resetCamera ref.
   * It restores both the camera position and the OrbitControls target.
   */
  const resetCamera = useCallback(() => {
    if (!controlsRef.current) return;

    const { position, target } = defaultsRef.current;
    camera.position.set(...position);
    controlsRef.current.target.set(...target);
    controlsRef.current.update();
  }, [camera]);

  // Auto-reset camera when the molecule changes (system switch).
  // Track whether this is the initial mount to avoid resetting on first render
  // (the Canvas camera prop handles the initial position).
  const isInitialMountRef = useRef(true);
  useEffect(() => {
    if (isInitialMountRef.current) {
      isInitialMountRef.current = false;
      return;
    }
    // Molecule changed -- reset camera to frame the new geometry
    if (controlsRef.current) {
      camera.position.set(...defaultCameraPosition);
      controlsRef.current.target.set(...defaultCameraTarget);
      controlsRef.current.update();
    }
  }, [defaultCameraPosition, defaultCameraTarget, camera]);

  // Expose resetCamera to parent via window.
  // This is a pragmatic bridge between the R3F scene context (inside Canvas)
  // and the HTML overlay button (outside Canvas). The alternative -- complex
  // ref forwarding through Canvas -- adds significant complexity for no benefit.
  useMemo(() => {
    (window as unknown as Record<string, unknown>).__iqcp_reset_camera = resetCamera;
    return resetCamera;
  }, [resetCamera]);

  return (
    <>
      {/* Ambient light for base illumination */}
      <ambientLight intensity={0.6} />

      {/* Directional lights for depth and shading */}
      <directionalLight position={[10, 10, 10]} intensity={0.8} />
      <directionalLight position={[-5, -5, 5]} intensity={0.3} />

      {/* Render atom spheres and bond cylinders */}
      <AtomSpheres atoms={atoms} />
      <BondCylinders atoms={atoms} />

      {/* Render atom labels (conditionally) */}
      <AtomLabels atoms={atoms} visible={showLabels} />

      {/* Additional scene children (e.g., OrbitalSurface for isosurface rendering) */}
      {children}

      {/* OrbitControls: click-drag rotate, scroll zoom, right-drag pan */}
      <OrbitControls
        ref={controlsRef}
        makeDefault
        enableDamping={!prefersReducedMotion}
        dampingFactor={0.1}
        onChange={handleCameraChange}
      />
    </>
  );
}

/**
 * 3D molecular viewer using React Three Fiber.
 *
 * Renders atoms as CPK-colored spheres with bond cylinders and optional
 * billboard atom labels in a 3D scene with interactive OrbitControls.
 *
 * The camera is automatically positioned to frame all atoms in the molecule.
 * Users can rotate (left-drag), zoom (scroll), and pan (right-drag) the view.
 * A reset button restores the default camera position.
 *
 * Features:
 * - Instanced mesh rendering for performance (AtomSpheres)
 * - Bond detection and cylinder rendering (BondCylinders)
 * - Togglable billboard atom labels (AtomLabels, US-034)
 * - OrbitControls with camera state tracking (US-034)
 * - Reset camera button (US-034)
 * - WebGL fallback for unsupported browsers
 * - ARIA accessibility (role="img" + descriptive aria-label)
 * - Respects prefers-reduced-motion media query (US-034)
 *
 * US-044: Accepts optional children for rendering additional scene
 * content (e.g., OrbitalSurface for MO isosurface visualization).
 */
export function MoleculeViewer({
  atoms,
  ariaLabel,
  className = '',
  showLabels = false,
  onCameraChange,
  children,
}: MoleculeViewerProps) {
  // Check WebGL support (once)
  const webGLAvailable = useMemo(() => isWebGLAvailable(), []);

  // Check reduced motion preference (once)
  const prefersReducedMotion = useMemo(() => checkPrefersReducedMotion(), []);

  // Compute camera position to frame the molecule
  const { position: cameraPosition, target: cameraTarget } = useMemo(
    () => computeCameraPosition(atoms),
    [atoms]
  );

  // Generate a descriptive ARIA label if none provided
  const effectiveAriaLabel = useMemo(() => {
    if (ariaLabel) return ariaLabel;
    if (atoms.length === 0) return '3D molecular viewer (no molecule loaded)';

    // Build a simple description from atom symbols
    const symbolCounts: Record<string, number> = {};
    for (const atom of atoms) {
      symbolCounts[atom.symbol] = (symbolCounts[atom.symbol] ?? 0) + 1;
    }
    const formula = Object.entries(symbolCounts)
      .map(([sym, count]) => (count > 1 ? `${sym}${count}` : sym))
      .join('');
    return `3D view of molecule ${formula}`;
  }, [atoms, ariaLabel]);

  /**
   * Reset camera to the default auto-computed position.
   *
   * Calls the reset function stored by SceneContent via window.
   * This is a pragmatic approach to communicate between the HTML overlay
   * button and the R3F scene context without complex ref forwarding.
   */
  const handleResetCamera = useCallback(() => {
    const resetFn = (window as unknown as Record<string, unknown>).__iqcp_reset_camera;
    if (typeof resetFn === 'function') {
      (resetFn as () => void)();
    }
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
        <SceneContent
          atoms={atoms}
          showLabels={showLabels}
          defaultCameraPosition={cameraPosition}
          defaultCameraTarget={cameraTarget}
          prefersReducedMotion={prefersReducedMotion}
          onCameraChange={onCameraChange}
        >
          {children}
        </SceneContent>
      </Canvas>

      {/* Reset camera button -- HTML overlay in top-right corner */}
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
