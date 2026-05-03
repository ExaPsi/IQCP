/**
 * DensitySurface - R3F component for rendering electron density isosurfaces.
 *
 * Renders the total electron density as a single semi-transparent mesh
 * inside the Three.js Canvas. Unlike OrbitalSurface (which renders two
 * lobes for positive/negative MO values), DensitySurface renders only
 * one surface because the electron density rho(r) is always non-negative.
 *
 * The mesh uses a teal/cyan color (#33ccaa) at 35% opacity so that
 * the underlying molecular structure (atoms and bonds) remains visible.
 *
 * Geometry is created from flat arrays (vertices, normals, indices)
 * produced by the marching cubes worker. Resources are properly cleaned
 * up when props change or the component unmounts.
 *
 * @module viewer3d/DensitySurface
 * @see US-062 Density Isosurface & Cross-Sections
 */

import { useMemo, useEffect } from 'react';
import * as THREE from 'three';
import type { LobeMeshData } from './OrbitalSurface';

// ============================================================================
// Types
// ============================================================================

/**
 * Props for the DensitySurface component.
 *
 * Unlike OrbitalSurface (which accepts positive + negative lobes),
 * DensitySurface accepts a single mesh because electron density
 * rho(r) >= 0 everywhere -- there are no "negative lobes."
 */
export interface DensitySurfaceProps {
  /** Mesh data from marching cubes (vertices, indices, normals), or null if none */
  meshData: LobeMeshData | null;
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
 * DensitySurface - renders electron density isosurface.
 *
 * This component lives inside a React Three Fiber <Canvas> and renders
 * a single semi-transparent teal mesh representing the electron density
 * at a given isovalue threshold.
 *
 * The surface is rendered with:
 * - Teal color (#33ccaa) distinct from orbital blue/orange
 * - 35% opacity for semi-transparency (AC3)
 * - Double-sided rendering for correct appearance from all angles
 * - depthWrite disabled for proper transparency compositing
 *
 * @example
 * ```tsx
 * <Canvas>
 *   <DensitySurface
 *     meshData={{ vertices: [...], indices: [...], normals: [...] }}
 *   />
 * </Canvas>
 * ```
 */
export function DensitySurface({ meshData }: DensitySurfaceProps) {
  // Create geometry from mesh data, memoized to avoid recreation on every render.
  const geometry = useMemo(
    () => (meshData ? createGeometry(meshData) : null),
    [meshData]
  );

  // Dispose old geometry on change or unmount to prevent GPU memory leaks.
  useEffect(() => {
    return () => {
      geometry?.dispose();
    };
  }, [geometry]);

  if (!geometry) return null;

  return (
    <mesh geometry={geometry}>
      <meshStandardMaterial
        color="#33ccaa"
        transparent
        opacity={0.35}
        side={THREE.DoubleSide}
        depthWrite={false}
      />
    </mesh>
  );
}
