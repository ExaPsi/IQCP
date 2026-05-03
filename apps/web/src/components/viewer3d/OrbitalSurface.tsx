/**
 * OrbitalSurface - R3F component for rendering orbital isosurfaces.
 *
 * Renders positive and negative orbital lobes as separate meshes inside
 * the Three.js Canvas. Positive lobe uses solid blue material; negative
 * lobe uses translucent red/orange material (not relying on color alone
 * for accessibility -- solid vs translucent provides shape distinction).
 *
 * Geometry is created from flat arrays (vertices, normals, indices)
 * produced by the marching cubes worker. Resources are properly cleaned
 * up when props change or the component unmounts.
 *
 * @module viewer3d/OrbitalSurface
 * @see US-044 Orbital Viewer UI
 */

import { useMemo, useEffect, useRef } from 'react';
import * as THREE from 'three';

// ============================================================================
// Types
// ============================================================================

/**
 * Mesh data for a single isosurface lobe.
 *
 * Contains interleaved flat arrays directly from marching cubes output,
 * ready to be consumed by Three.js BufferGeometry.
 */
export interface LobeMeshData {
  /** Interleaved vertex positions [x0,y0,z0, x1,y1,z1, ...] */
  vertices: number[];
  /** Triangle indices (3 per triangle) */
  indices: number[];
  /** Interleaved vertex normals [nx0,ny0,nz0, nx1,ny1,nz1, ...] */
  normals: number[];
}

/**
 * Props for the OrbitalSurface component.
 */
export interface OrbitalSurfaceProps {
  /** Positive lobe mesh data (field > +isovalue), or null if none */
  positive: LobeMeshData | null;
  /** Negative lobe mesh data (field < -isovalue), or null if none */
  negative: LobeMeshData | null;
  /**
   * Color for the positive lobe mesh.
   * Defaults to blue (#4488ff) for orbital visualization.
   * Use green (#22cc88) for difference density accumulation.
   * @see US-063 Difference Density
   */
  positiveColor?: string;
  /**
   * Opacity for the positive lobe mesh.
   * Defaults to 0.6 for orbital visualization.
   * Use 0.65 for difference density accumulation (solid).
   * @see US-063 Difference Density
   */
  positiveOpacity?: number;
  /**
   * Color for the negative lobe mesh.
   * Defaults to orange (#ff8844) for orbital visualization.
   * Use red (#cc4444) for difference density depletion.
   * @see US-063 Difference Density
   */
  negativeColor?: string;
  /**
   * Opacity for the negative lobe mesh.
   * Defaults to 0.4 for orbital visualization.
   * Use 0.35 for difference density depletion (translucent).
   * @see US-063 Difference Density
   */
  negativeOpacity?: number;
}

// ============================================================================
// Helpers
// ============================================================================

/**
 * Create a BufferGeometry from flat mesh data arrays.
 *
 * Converts the flat number arrays from marching cubes into typed arrays
 * and sets them as attributes on a Three.js BufferGeometry.
 *
 * @param data - Mesh data with vertices, normals, and indices
 * @returns BufferGeometry ready for rendering, or null if data is empty
 */
function createGeometry(data: LobeMeshData): THREE.BufferGeometry | null {
  if (data.vertices.length === 0 || data.indices.length === 0) {
    return null;
  }

  const geometry = new THREE.BufferGeometry();

  // Set vertex positions (itemSize = 3 for x,y,z)
  const positionArray = new Float32Array(data.vertices);
  geometry.setAttribute('position', new THREE.BufferAttribute(positionArray, 3));

  // Set vertex normals (itemSize = 3 for nx,ny,nz)
  if (data.normals.length > 0) {
    const normalArray = new Float32Array(data.normals);
    geometry.setAttribute('normal', new THREE.BufferAttribute(normalArray, 3));
  } else {
    // Compute normals from faces if not provided
    geometry.computeVertexNormals();
  }

  // Set triangle indices
  const indexArray = new Uint32Array(data.indices);
  geometry.setIndex(new THREE.BufferAttribute(indexArray, 1));

  return geometry;
}

// ============================================================================
// Component
// ============================================================================

/**
 * OrbitalSurface - renders positive and negative orbital isosurface lobes.
 *
 * This component lives inside a React Three Fiber <Canvas> and renders
 * two meshes: one for the positive lobe (solid blue) and one for the
 * negative lobe (translucent orange). Each mesh uses BufferGeometry
 * constructed from flat arrays produced by marching cubes.
 *
 * Accessibility: positive and negative lobes are distinguished by both
 * color (blue vs orange) and opacity (solid vs translucent), ensuring
 * the distinction does not rely on color alone.
 *
 * @example
 * ```tsx
 * <Canvas>
 *   <OrbitalSurface
 *     positive={{ vertices: [...], indices: [...], normals: [...] }}
 *     negative={{ vertices: [...], indices: [...], normals: [...] }}
 *   />
 * </Canvas>
 * ```
 */
export function OrbitalSurface({
  positive,
  negative,
  positiveColor = '#4488ff',
  positiveOpacity = 0.6,
  negativeColor = '#ff8844',
  negativeOpacity = 0.4,
}: OrbitalSurfaceProps) {
  const positiveMeshRef = useRef<THREE.Mesh>(null);
  const negativeMeshRef = useRef<THREE.Mesh>(null);

  // Create geometries from mesh data, memoized to avoid recreation on every render.
  // When the data changes, old geometries are disposed by the cleanup effect below.
  const positiveGeometry = useMemo(
    () => (positive ? createGeometry(positive) : null),
    [positive]
  );

  const negativeGeometry = useMemo(
    () => (negative ? createGeometry(negative) : null),
    [negative]
  );

  // Dispose old geometries on change or unmount to prevent GPU memory leaks.
  useEffect(() => {
    return () => {
      positiveGeometry?.dispose();
    };
  }, [positiveGeometry]);

  useEffect(() => {
    return () => {
      negativeGeometry?.dispose();
    };
  }, [negativeGeometry]);

  return (
    <>
      {/* Positive lobe (default: translucent blue for orbitals, solid green for diff density) */}
      {positiveGeometry && (
        <mesh ref={positiveMeshRef} geometry={positiveGeometry}>
          <meshStandardMaterial
            color={positiveColor}
            transparent
            opacity={positiveOpacity}
            side={THREE.DoubleSide}
            depthWrite={false}
          />
        </mesh>
      )}

      {/* Negative lobe (default: translucent orange for orbitals, translucent red for diff density) */}
      {negativeGeometry && (
        <mesh ref={negativeMeshRef} geometry={negativeGeometry}>
          <meshStandardMaterial
            color={negativeColor}
            transparent
            opacity={negativeOpacity}
            side={THREE.DoubleSide}
            depthWrite={false}
          />
        </mesh>
      )}
    </>
  );
}
