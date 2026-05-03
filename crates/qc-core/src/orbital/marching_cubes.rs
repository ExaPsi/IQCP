//! Marching cubes isosurface extraction
//!
//! Extracts a triangle mesh isosurface from a 3D scalar field.
//! Returns vertex positions, triangle indices, and smooth normals
//! directly consumable by Three.js `BufferGeometry`.
//!
//! # Algorithm
//!
//! The standard marching cubes algorithm (Lorensen & Cline, 1987):
//!
//! 1. Iterate over all cells in the 3D grid (each cell has 8 corner vertices).
//! 2. For each cell, determine the cube index from the 8 corner values
//!    relative to the isovalue (256 possible configurations).
//! 3. Use lookup tables to determine intersected edges and triangle topology.
//! 4. Interpolate vertex positions along intersected edges.
//! 5. Compute smooth normals via central-difference gradient estimation.
//!
//! # Output Format
//!
//! - `vertices: Vec<f32>` — interleaved [x0,y0,z0, x1,y1,z1, ...] (GPU-native f32)
//! - `indices: Vec<u32>` — triangle indices (3 per triangle, u32 for >65536 vertices)
//! - `normals: Vec<f32>` — interleaved [nx0,ny0,nz0, ...] (unit-length smooth normals)
//!
//! # Determinism
//!
//! The output is fully deterministic: same input always produces identical output.
//!
//! # References
//!
//! - Lorensen, W.E. & Cline, H.E. (1987). "Marching Cubes: A High Resolution
//!   3D Surface Construction Algorithm." Computer Graphics, 21(4), 163-169.
//! - Phase 2 TDD Section 8.3

use serde::{Deserialize, Serialize};

// =============================================================================
// Result Type
// =============================================================================

/// Result from marching cubes isosurface extraction.
///
/// Vertex/index/normal arrays are directly consumable by Three.js BufferGeometry.
/// Vertices and normals use f32 for GPU efficiency; indices use u32 to support
/// meshes with >65536 vertices.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MarchingCubesResult {
    /// Interleaved vertex positions [x0, y0, z0, x1, y1, z1, ...]
    pub vertices: Vec<f32>,
    /// Triangle indices (3 per triangle)
    pub indices: Vec<u32>,
    /// Interleaved vertex normals [nx0, ny0, nz0, nx1, ny1, nz1, ...]
    /// Smooth normals computed via central-difference gradient estimation
    pub normals: Vec<f32>,
}

// =============================================================================
// Public API
// =============================================================================

/// Extract an isosurface from a 3D scalar field using marching cubes.
///
/// # Arguments
///
/// * `grid` - Scalar field values (C-order: x-slowest, z-fastest).
///   Length must equal `dims[0] * dims[1] * dims[2]`.
/// * `dims` - Grid dimensions [nx, ny, nz], each >= 2.
/// * `origin` - Grid origin [x, y, z] in world coordinates.
/// * `spacing` - Uniform grid spacing.
/// * `isovalue` - Isosurface threshold value.
///
/// # Returns
///
/// A `MarchingCubesResult` with triangle mesh data. Returns an empty mesh
/// if no cells intersect the isosurface.
///
/// # Panics
///
/// Panics if `grid.len() != dims[0] * dims[1] * dims[2]` or if any
/// dimension is less than 2.
pub fn marching_cubes(
    grid: &[f64],
    dims: [usize; 3],
    origin: [f64; 3],
    spacing: f64,
    isovalue: f64,
) -> MarchingCubesResult {
    let [nx, ny, nz] = dims;
    assert_eq!(grid.len(), nx * ny * nz, "Grid size mismatch");
    assert!(nx >= 2 && ny >= 2 && nz >= 2, "Grid dims must be >= 2");

    let mut vertices: Vec<f32> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    let mut normals: Vec<f32> = Vec::new();

    // Precompute gradients at all grid points (for smooth normals)
    let gradients = compute_all_gradients(grid, dims, spacing);

    // Iterate over all cells (a cell spans from [ix, ix+1] x [iy, iy+1] x [iz, iz+1])
    for ix in 0..nx - 1 {
        for iy in 0..ny - 1 {
            for iz in 0..nz - 1 {
                process_cell(
                    grid,
                    &gradients,
                    dims,
                    origin,
                    spacing,
                    isovalue,
                    ix,
                    iy,
                    iz,
                    &mut vertices,
                    &mut indices,
                    &mut normals,
                );
            }
        }
    }

    MarchingCubesResult {
        vertices,
        indices,
        normals,
    }
}

/// Extract dual isosurfaces for orbital visualization.
///
/// Orbital wavefunctions have positive and negative lobes. This function
/// extracts isosurfaces at `+isovalue` and `-isovalue` in a single call.
///
/// # Returns
///
/// A tuple of `(positive_lobe, negative_lobe)` meshes.
pub fn dual_marching_cubes(
    grid: &[f64],
    dims: [usize; 3],
    origin: [f64; 3],
    spacing: f64,
    isovalue: f64,
) -> (MarchingCubesResult, MarchingCubesResult) {
    let positive = marching_cubes(grid, dims, origin, spacing, isovalue);
    let negative = marching_cubes(grid, dims, origin, spacing, -isovalue);
    (positive, negative)
}

// =============================================================================
// Cell Processing
// =============================================================================

/// Process a single cell and emit triangles.
#[allow(clippy::too_many_arguments)]
fn process_cell(
    grid: &[f64],
    gradients: &[[f64; 3]],
    dims: [usize; 3],
    origin: [f64; 3],
    spacing: f64,
    isovalue: f64,
    ix: usize,
    iy: usize,
    iz: usize,
    vertices: &mut Vec<f32>,
    indices: &mut Vec<u32>,
    normals: &mut Vec<f32>,
) {
    let [_nx, ny, nz] = dims;

    // Corner indices in the grid (standard cube vertex numbering)
    //   4----5
    //  /|   /|
    // 7----6 |   y up, x right, z toward viewer
    // | 0--|-1
    // |/   |/
    // 3----2
    let corner_indices = [
        flat_idx(ix, iy, iz, ny, nz),             // 0
        flat_idx(ix + 1, iy, iz, ny, nz),         // 1
        flat_idx(ix + 1, iy, iz + 1, ny, nz),     // 2
        flat_idx(ix, iy, iz + 1, ny, nz),         // 3
        flat_idx(ix, iy + 1, iz, ny, nz),         // 4
        flat_idx(ix + 1, iy + 1, iz, ny, nz),     // 5
        flat_idx(ix + 1, iy + 1, iz + 1, ny, nz), // 6
        flat_idx(ix, iy + 1, iz + 1, ny, nz),     // 7
    ];

    // Read corner values
    let vals: [f64; 8] = [
        grid[corner_indices[0]],
        grid[corner_indices[1]],
        grid[corner_indices[2]],
        grid[corner_indices[3]],
        grid[corner_indices[4]],
        grid[corner_indices[5]],
        grid[corner_indices[6]],
        grid[corner_indices[7]],
    ];

    // Compute cube index (bit i set if corner i is above isovalue)
    let mut cube_index: u8 = 0;
    for (i, &val) in vals.iter().enumerate() {
        if val > isovalue {
            cube_index |= 1 << i;
        }
    }

    // Look up which edges are intersected
    let edge_bits = EDGE_TABLE[cube_index as usize];
    if edge_bits == 0 {
        return; // No intersection
    }

    // Corner positions in world coordinates
    let corner_pos: [[f64; 3]; 8] = [
        pos3(origin, spacing, ix, iy, iz),
        pos3(origin, spacing, ix + 1, iy, iz),
        pos3(origin, spacing, ix + 1, iy, iz + 1),
        pos3(origin, spacing, ix, iy, iz + 1),
        pos3(origin, spacing, ix, iy + 1, iz),
        pos3(origin, spacing, ix + 1, iy + 1, iz),
        pos3(origin, spacing, ix + 1, iy + 1, iz + 1),
        pos3(origin, spacing, ix, iy + 1, iz + 1),
    ];

    // Corner gradients
    let corner_grads: [[f64; 3]; 8] = [
        gradients[corner_indices[0]],
        gradients[corner_indices[1]],
        gradients[corner_indices[2]],
        gradients[corner_indices[3]],
        gradients[corner_indices[4]],
        gradients[corner_indices[5]],
        gradients[corner_indices[6]],
        gradients[corner_indices[7]],
    ];

    // Edge endpoint pairs (vertex indices for each of the 12 edges)
    const EDGE_VERTICES: [[usize; 2]; 12] = [
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 0], // bottom face
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 4], // top face
        [0, 4],
        [1, 5],
        [2, 6],
        [3, 7], // vertical edges
    ];

    // Interpolate vertex positions and normals for each intersected edge
    let mut edge_verts: [[f32; 3]; 12] = [[0.0; 3]; 12];
    let mut edge_norms: [[f32; 3]; 12] = [[0.0; 3]; 12];

    for edge in 0..12u32 {
        if edge_bits & (1 << edge) != 0 {
            let [a, b] = EDGE_VERTICES[edge as usize];
            let (vert, norm) = interpolate_edge(
                &corner_pos[a],
                &corner_pos[b],
                vals[a],
                vals[b],
                &corner_grads[a],
                &corner_grads[b],
                isovalue,
            );
            edge_verts[edge as usize] = vert;
            edge_norms[edge as usize] = norm;
        }
    }

    // Emit triangles from the triangle table
    let tri_row = &TRI_TABLE[cube_index as usize];
    let mut i = 0;
    while i < 16 && tri_row[i] != -1 {
        let e0 = tri_row[i] as usize;
        let e1 = tri_row[i + 1] as usize;
        let e2 = tri_row[i + 2] as usize;

        let base_idx = (vertices.len() / 3) as u32;

        // Add three vertices
        vertices.extend_from_slice(&edge_verts[e0]);
        vertices.extend_from_slice(&edge_verts[e1]);
        vertices.extend_from_slice(&edge_verts[e2]);

        // Add three normals
        normals.extend_from_slice(&edge_norms[e0]);
        normals.extend_from_slice(&edge_norms[e1]);
        normals.extend_from_slice(&edge_norms[e2]);

        // Add triangle indices
        indices.push(base_idx);
        indices.push(base_idx + 1);
        indices.push(base_idx + 2);

        i += 3;
    }
}

// =============================================================================
// Vertex Interpolation
// =============================================================================

/// Linearly interpolate a vertex position and normal along an edge.
///
/// Given two edge endpoints with positions, scalar values, and gradients,
/// find the point where the scalar field equals `isovalue` and interpolate
/// the gradient-based normal at that point.
#[inline]
fn interpolate_edge(
    p1: &[f64; 3],
    p2: &[f64; 3],
    v1: f64,
    v2: f64,
    g1: &[f64; 3],
    g2: &[f64; 3],
    iso: f64,
) -> ([f32; 3], [f32; 3]) {
    let dv = v2 - v1;
    let t = if dv.abs() < 1e-30 {
        0.5
    } else {
        (iso - v1) / dv
    };
    // Clamp to avoid floating point overshoot
    let t = t.clamp(0.0, 1.0);

    let pos = [
        (p1[0] + t * (p2[0] - p1[0])) as f32,
        (p1[1] + t * (p2[1] - p1[1])) as f32,
        (p1[2] + t * (p2[2] - p1[2])) as f32,
    ];

    // Interpolate gradient
    let gx = g1[0] + t * (g2[0] - g1[0]);
    let gy = g1[1] + t * (g2[1] - g1[1]);
    let gz = g1[2] + t * (g2[2] - g1[2]);

    // Normalize and negate (normals point away from higher values)
    let len = (gx * gx + gy * gy + gz * gz).sqrt();
    let norm = if len > 1e-10 {
        [(-gx / len) as f32, (-gy / len) as f32, (-gz / len) as f32]
    } else {
        [0.0, 0.0, 1.0] // Default normal for zero gradient
    };

    (pos, norm)
}

// =============================================================================
// Gradient Computation
// =============================================================================

/// Compute gradients at all grid points using central differences.
///
/// At interior points: grad_x = (f[ix+1] - f[ix-1]) / (2 * spacing)
/// At boundary points: forward or backward differences.
fn compute_all_gradients(grid: &[f64], dims: [usize; 3], spacing: f64) -> Vec<[f64; 3]> {
    let [nx, ny, nz] = dims;
    let total = nx * ny * nz;
    let mut gradients = vec![[0.0f64; 3]; total];
    let inv_2h = 1.0 / (2.0 * spacing);
    let inv_h = 1.0 / spacing;

    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let idx = flat_idx(ix, iy, iz, ny, nz);

                // X gradient
                let gx = if ix == 0 {
                    (grid[flat_idx(1, iy, iz, ny, nz)] - grid[idx]) * inv_h
                } else if ix == nx - 1 {
                    (grid[idx] - grid[flat_idx(nx - 2, iy, iz, ny, nz)]) * inv_h
                } else {
                    (grid[flat_idx(ix + 1, iy, iz, ny, nz)]
                        - grid[flat_idx(ix - 1, iy, iz, ny, nz)])
                        * inv_2h
                };

                // Y gradient
                let gy = if iy == 0 {
                    (grid[flat_idx(ix, 1, iz, ny, nz)] - grid[idx]) * inv_h
                } else if iy == ny - 1 {
                    (grid[idx] - grid[flat_idx(ix, ny - 2, iz, ny, nz)]) * inv_h
                } else {
                    (grid[flat_idx(ix, iy + 1, iz, ny, nz)]
                        - grid[flat_idx(ix, iy - 1, iz, ny, nz)])
                        * inv_2h
                };

                // Z gradient
                let gz = if iz == 0 {
                    (grid[flat_idx(ix, iy, 1, ny, nz)] - grid[idx]) * inv_h
                } else if iz == nz - 1 {
                    (grid[idx] - grid[flat_idx(ix, iy, nz - 2, ny, nz)]) * inv_h
                } else {
                    (grid[flat_idx(ix, iy, iz + 1, ny, nz)]
                        - grid[flat_idx(ix, iy, iz - 1, ny, nz)])
                        * inv_2h
                };

                gradients[idx] = [gx, gy, gz];
            }
        }
    }

    gradients
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Flat index from 3D coordinates (C-order: x-slowest, z-fastest).
#[inline]
fn flat_idx(ix: usize, iy: usize, iz: usize, ny: usize, nz: usize) -> usize {
    ix * ny * nz + iy * nz + iz
}

/// World position from grid indices.
#[inline]
fn pos3(origin: [f64; 3], spacing: f64, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
    [
        origin[0] + ix as f64 * spacing,
        origin[1] + iy as f64 * spacing,
        origin[2] + iz as f64 * spacing,
    ]
}

// =============================================================================
// Lookup Tables
// =============================================================================

/// Edge table: for each of the 256 cube configurations, which edges are
/// intersected by the isosurface (bitmask of 12 edges).
///
/// Reference: Lorensen & Cline (1987), standard marching cubes table.
#[rustfmt::skip]
const EDGE_TABLE: [u16; 256] = [
    0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000,
];

/// Triangle table: for each of the 256 cube configurations, lists the
/// edge indices forming triangles. Terminated by -1.
///
/// Each row has at most 5 triangles (15 edge indices) plus a terminator.
/// Edge indices refer to the 12 edges of the cube (0-11).
///
/// Reference: Lorensen & Cline (1987), standard marching cubes table.
#[rustfmt::skip]
const TRI_TABLE: [[i8; 16]; 256] = [
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 1, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 8, 3, 9, 8, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3, 1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 2,10, 0, 2, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 8, 3, 2,10, 8,10, 9, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0,11, 2, 8,11, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 9, 0, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1,11, 2, 1, 9,11, 9, 8,11,-1,-1,-1,-1,-1,-1,-1],
    [ 3,10, 1,11,10, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0,10, 1, 0, 8,10, 8,11,10,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 9, 0, 3,11, 9,11,10, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 3, 0, 7, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 1, 9, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 1, 9, 4, 7, 1, 7, 3, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 4, 7, 3, 0, 4, 1, 2,10,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 2,10, 9, 0, 2, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 2,10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4,-1,-1,-1,-1],
    [ 8, 4, 7, 3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11, 4, 7,11, 2, 4, 2, 0, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 0, 1, 8, 4, 7, 2, 3,11,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 7,11, 9, 4,11, 9,11, 2, 9, 2, 1,-1,-1,-1,-1],
    [ 3,10, 1, 3,11,10, 7, 8, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 1,11,10, 1, 4,11, 1, 0, 4, 7,11, 4,-1,-1,-1,-1],
    [ 4, 7, 8, 9, 0,11, 9,11,10,11, 0, 3,-1,-1,-1,-1],
    [ 4, 7,11, 4,11, 9, 9,11,10,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 5, 4, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 5, 4, 1, 5, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 5, 4, 8, 3, 5, 3, 1, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 0, 8, 1, 2,10, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 5, 2,10, 5, 4, 2, 4, 0, 2,-1,-1,-1,-1,-1,-1,-1],
    [ 2,10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8,-1,-1,-1,-1],
    [ 9, 5, 4, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0,11, 2, 0, 8,11, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 5, 4, 0, 1, 5, 2, 3,11,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 1, 5, 2, 5, 8, 2, 8,11, 4, 8, 5,-1,-1,-1,-1],
    [10, 3,11,10, 1, 3, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 9, 5, 0, 8, 1, 8,10, 1, 8,11,10,-1,-1,-1,-1],
    [ 5, 4, 0, 5, 0,11, 5,11,10,11, 0, 3,-1,-1,-1,-1],
    [ 5, 4, 8, 5, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 7, 8, 5, 7, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 3, 0, 9, 5, 3, 5, 7, 3,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 7, 8, 0, 1, 7, 1, 5, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 7, 8, 9, 5, 7,10, 1, 2,-1,-1,-1,-1,-1,-1,-1],
    [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3,-1,-1,-1,-1],
    [ 8, 0, 2, 8, 2, 5, 8, 5, 7,10, 5, 2,-1,-1,-1,-1],
    [ 2,10, 5, 2, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 7, 9, 5, 7, 8, 9, 3,11, 2,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7,11,-1,-1,-1,-1],
    [ 2, 3,11, 0, 1, 8, 1, 7, 8, 1, 5, 7,-1,-1,-1,-1],
    [11, 2, 1,11, 1, 7, 7, 1, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 5, 8, 8, 5, 7,10, 1, 3,10, 3,11,-1,-1,-1,-1],
    [ 5, 7, 0, 5, 0, 9, 7,11, 0, 1, 0,10,11,10, 0,-1],
    [11,10, 0,11, 0, 3,10, 5, 0, 8, 0, 7, 5, 7, 0,-1],
    [11,10, 5, 7,11, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 0, 1, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 8, 3, 1, 9, 8, 5,10, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 6, 5, 2, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 6, 5, 1, 2, 6, 3, 0, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 6, 5, 9, 0, 6, 0, 2, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8,-1,-1,-1,-1],
    [ 2, 3,11,10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11, 0, 8,11, 2, 0,10, 6, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 1, 9, 2, 3,11, 5,10, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 5,10, 6, 1, 9, 2, 9,11, 2, 9, 8,11,-1,-1,-1,-1],
    [ 6, 3,11, 6, 5, 3, 5, 1, 3,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8,11, 0,11, 5, 0, 5, 1, 5,11, 6,-1,-1,-1,-1],
    [ 3,11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9,-1,-1,-1,-1],
    [ 6, 5, 9, 6, 9,11,11, 9, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 5,10, 6, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 3, 0, 4, 7, 3, 6, 5,10,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 9, 0, 5,10, 6, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1],
    [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4,-1,-1,-1,-1],
    [ 6, 1, 2, 6, 5, 1, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7,-1,-1,-1,-1],
    [ 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6,-1,-1,-1,-1],
    [ 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9,-1],
    [ 3,11, 2, 7, 8, 4,10, 6, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 5,10, 6, 4, 7, 2, 4, 2, 0, 2, 7,11,-1,-1,-1,-1],
    [ 0, 1, 9, 4, 7, 8, 2, 3,11, 5,10, 6,-1,-1,-1,-1],
    [ 9, 2, 1, 9,11, 2, 9, 4,11, 7,11, 4, 5,10, 6,-1],
    [ 8, 4, 7, 3,11, 5, 3, 5, 1, 5,11, 6,-1,-1,-1,-1],
    [ 5, 1,11, 5,11, 6, 1, 0,11, 7,11, 4, 0, 4,11,-1],
    [ 0, 5, 9, 0, 6, 5, 0, 3, 6,11, 6, 3, 8, 4, 7,-1],
    [ 6, 5, 9, 6, 9,11, 4, 7, 9, 7,11, 9,-1,-1,-1,-1],
    [10, 4, 9, 6, 4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4,10, 6, 4, 9,10, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1],
    [10, 0, 1,10, 6, 0, 6, 4, 0,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1,10,-1,-1,-1,-1],
    [ 1, 4, 9, 1, 2, 4, 2, 6, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4,-1,-1,-1,-1],
    [ 0, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 3, 2, 8, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1],
    [10, 4, 9,10, 6, 4,11, 2, 3,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 2, 2, 8,11, 4, 9,10, 4,10, 6,-1,-1,-1,-1],
    [ 3,11, 2, 0, 1, 6, 0, 6, 4, 6, 1,10,-1,-1,-1,-1],
    [ 6, 4, 1, 6, 1,10, 4, 8, 1, 2, 1,11, 8,11, 1,-1],
    [ 9, 6, 4, 9, 3, 6, 9, 1, 3,11, 6, 3,-1,-1,-1,-1],
    [ 8,11, 1, 8, 1, 0,11, 6, 1, 9, 1, 4, 6, 4, 1,-1],
    [ 3,11, 6, 3, 6, 0, 0, 6, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 6, 4, 8,11, 6, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 7,10, 6, 7, 8,10, 8, 9,10,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 7, 3, 0,10, 7, 0, 9,10, 6, 7,10,-1,-1,-1,-1],
    [10, 6, 7, 1,10, 7, 1, 7, 8, 1, 8, 0,-1,-1,-1,-1],
    [10, 6, 7,10, 7, 1, 1, 7, 3,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7,-1,-1,-1,-1],
    [ 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9,-1],
    [ 7, 8, 0, 7, 0, 6, 6, 0, 2,-1,-1,-1,-1,-1,-1,-1],
    [ 7, 3, 2, 6, 7, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 3,11,10, 6, 8,10, 8, 9, 8, 6, 7,-1,-1,-1,-1],
    [ 2, 0, 7, 2, 7,11, 0, 9, 7, 6, 7,10, 9,10, 7,-1],
    [ 1, 8, 0, 1, 7, 8, 1,10, 7, 6, 7,10, 2, 3,11,-1],
    [11, 2, 1,11, 1, 7,10, 6, 1, 6, 7, 1,-1,-1,-1,-1],
    [ 8, 9, 6, 8, 6, 7, 9, 1, 6,11, 6, 3, 1, 3, 6,-1],
    [ 0, 9, 1,11, 6, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 7, 8, 0, 7, 0, 6, 3,11, 0,11, 6, 0,-1,-1,-1,-1],
    [ 7,11, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 0, 8,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 1, 9,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 1, 9, 8, 3, 1,11, 7, 6,-1,-1,-1,-1,-1,-1,-1],
    [10, 1, 2, 6,11, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10, 3, 0, 8, 6,11, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 9, 0, 2,10, 9, 6,11, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 6,11, 7, 2,10, 3,10, 8, 3,10, 9, 8,-1,-1,-1,-1],
    [ 7, 2, 3, 6, 2, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 7, 0, 8, 7, 6, 0, 6, 2, 0,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 7, 6, 2, 3, 7, 0, 1, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6,-1,-1,-1,-1],
    [10, 7, 6,10, 1, 7, 1, 3, 7,-1,-1,-1,-1,-1,-1,-1],
    [10, 7, 6, 1, 7,10, 1, 8, 7, 1, 0, 8,-1,-1,-1,-1],
    [ 0, 3, 7, 0, 7,10, 0,10, 9, 6,10, 7,-1,-1,-1,-1],
    [ 7, 6,10, 7,10, 8, 8,10, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 6, 8, 4,11, 8, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 6,11, 3, 0, 6, 0, 4, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 6,11, 8, 4, 6, 9, 0, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 4, 6, 9, 6, 3, 9, 3, 1,11, 3, 6,-1,-1,-1,-1],
    [ 6, 8, 4, 6,11, 8, 2,10, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10, 3, 0,11, 0, 6,11, 0, 4, 6,-1,-1,-1,-1],
    [ 4,11, 8, 4, 6,11, 0, 2, 9, 2,10, 9,-1,-1,-1,-1],
    [10, 9, 3,10, 3, 2, 9, 4, 3,11, 3, 6, 4, 6, 3,-1],
    [ 8, 2, 3, 8, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8,-1,-1,-1,-1],
    [ 1, 9, 4, 1, 4, 2, 2, 4, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 1, 3, 8, 6, 1, 8, 4, 6, 6,10, 1,-1,-1,-1,-1],
    [10, 1, 0,10, 0, 6, 6, 0, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 6, 3, 4, 3, 8, 6,10, 3, 0, 3, 9,10, 9, 3,-1],
    [10, 9, 4, 6,10, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 9, 5, 7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3, 4, 9, 5,11, 7, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 5, 0, 1, 5, 4, 0, 7, 6,11,-1,-1,-1,-1,-1,-1,-1],
    [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5,-1,-1,-1,-1],
    [ 9, 5, 4,10, 1, 2, 7, 6,11,-1,-1,-1,-1,-1,-1,-1],
    [ 6,11, 7, 1, 2,10, 0, 8, 3, 4, 9, 5,-1,-1,-1,-1],
    [ 7, 6,11, 5, 4,10, 4, 2,10, 4, 0, 2,-1,-1,-1,-1],
    [ 3, 4, 8, 3, 5, 4, 3, 2, 5,10, 5, 2,11, 7, 6,-1],
    [ 7, 2, 3, 7, 6, 2, 5, 4, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7,-1,-1,-1,-1],
    [ 3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0,-1,-1,-1,-1],
    [ 6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8,-1],
    [ 9, 5, 4,10, 1, 6, 1, 7, 6, 1, 3, 7,-1,-1,-1,-1],
    [ 1, 6,10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4,-1],
    [ 4, 0,10, 4,10, 5, 0, 3,10, 6,10, 7, 3, 7,10,-1],
    [ 7, 6,10, 7,10, 8, 5, 4,10, 4, 8,10,-1,-1,-1,-1],
    [ 6, 9, 5, 6,11, 9,11, 8, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 6,11, 0, 6, 3, 0, 5, 6, 0, 9, 5,-1,-1,-1,-1],
    [ 0,11, 8, 0, 5,11, 0, 1, 5, 5, 6,11,-1,-1,-1,-1],
    [ 6,11, 3, 6, 3, 5, 5, 3, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10, 9, 5,11, 9,11, 8,11, 5, 6,-1,-1,-1,-1],
    [ 0,11, 3, 0, 6,11, 0, 9, 6, 5, 6, 9, 1, 2,10,-1],
    [11, 8, 5,11, 5, 6, 8, 0, 5,10, 5, 2, 0, 2, 5,-1],
    [ 6,11, 3, 6, 3, 5, 2,10, 3,10, 5, 3,-1,-1,-1,-1],
    [ 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2,-1,-1,-1,-1],
    [ 9, 5, 6, 9, 6, 0, 0, 6, 2,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8,-1],
    [ 1, 5, 6, 2, 1, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 3, 6, 1, 6,10, 3, 8, 6, 5, 6, 9, 8, 9, 6,-1],
    [10, 1, 0,10, 0, 6, 9, 5, 0, 5, 6, 0,-1,-1,-1,-1],
    [ 0, 3, 8, 5, 6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [10, 5, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11, 5,10, 7, 5,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11, 5,10,11, 7, 5, 8, 3, 0,-1,-1,-1,-1,-1,-1,-1],
    [ 5,11, 7, 5,10,11, 1, 9, 0,-1,-1,-1,-1,-1,-1,-1],
    [10, 7, 5,10,11, 7, 9, 8, 1, 8, 3, 1,-1,-1,-1,-1],
    [11, 1, 2,11, 7, 1, 7, 5, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2,11,-1,-1,-1,-1],
    [ 9, 7, 5, 9, 2, 7, 9, 0, 2, 2,11, 7,-1,-1,-1,-1],
    [ 7, 5, 2, 7, 2,11, 5, 9, 2, 3, 2, 8, 9, 8, 2,-1],
    [ 2, 5,10, 2, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 2, 0, 8, 5, 2, 8, 7, 5,10, 2, 5,-1,-1,-1,-1],
    [ 9, 0, 1, 5,10, 3, 5, 3, 7, 3,10, 2,-1,-1,-1,-1],
    [ 9, 8, 2, 9, 2, 1, 8, 7, 2,10, 2, 5, 7, 5, 2,-1],
    [ 1, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 7, 0, 7, 1, 1, 7, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 0, 3, 9, 3, 5, 5, 3, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 8, 7, 5, 9, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 5, 8, 4, 5,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 5, 0, 4, 5,11, 0, 5,10,11,11, 3, 0,-1,-1,-1,-1],
    [ 0, 1, 9, 8, 4,10, 8,10,11,10, 4, 5,-1,-1,-1,-1],
    [10,11, 4,10, 4, 5,11, 3, 4, 9, 4, 1, 3, 1, 4,-1],
    [ 2, 5, 1, 2, 8, 5, 2,11, 8, 4, 5, 8,-1,-1,-1,-1],
    [ 0, 4,11, 0,11, 3, 4, 5,11, 2,11, 1, 5, 1,11,-1],
    [ 0, 2, 5, 0, 5, 9, 2,11, 5, 4, 5, 8,11, 8, 5,-1],
    [ 9, 4, 5, 2,11, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 5,10, 3, 5, 2, 3, 4, 5, 3, 8, 4,-1,-1,-1,-1],
    [ 5,10, 2, 5, 2, 4, 4, 2, 0,-1,-1,-1,-1,-1,-1,-1],
    [ 3,10, 2, 3, 5,10, 3, 8, 5, 4, 5, 8, 0, 1, 9,-1],
    [ 5,10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2,-1,-1,-1,-1],
    [ 8, 4, 5, 8, 5, 3, 3, 5, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 4, 5, 1, 0, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5,-1,-1,-1,-1],
    [ 9, 4, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4,11, 7, 4, 9,11, 9,10,11,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3, 4, 9, 7, 9,11, 7, 9,10,11,-1,-1,-1,-1],
    [ 1,10,11, 1,11, 4, 1, 4, 0, 7, 4,11,-1,-1,-1,-1],
    [ 3, 1, 4, 3, 4, 8, 1,10, 4, 7, 4,11,10,11, 4,-1],
    [ 4,11, 7, 9,11, 4, 9, 2,11, 9, 1, 2,-1,-1,-1,-1],
    [ 9, 7, 4, 9,11, 7, 9, 1,11, 2,11, 1, 0, 8, 3,-1],
    [11, 7, 4,11, 4, 2, 2, 4, 0,-1,-1,-1,-1,-1,-1,-1],
    [11, 7, 4,11, 4, 2, 8, 3, 4, 3, 2, 4,-1,-1,-1,-1],
    [ 2, 9,10, 2, 7, 9, 2, 3, 7, 7, 4, 9,-1,-1,-1,-1],
    [ 9,10, 7, 9, 7, 4,10, 2, 7, 8, 7, 0, 2, 0, 7,-1],
    [ 3, 7,10, 3,10, 2, 7, 4,10, 1,10, 0, 4, 0,10,-1],
    [ 1,10, 2, 8, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 9, 1, 4, 1, 7, 7, 1, 3,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1,-1,-1,-1,-1],
    [ 4, 0, 3, 7, 4, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 8, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 0, 9, 3, 9,11,11, 9,10,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 1,10, 0,10, 8, 8,10,11,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 1,10,11, 3,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,11, 1,11, 9, 9,11, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 0, 9, 3, 9,11, 1, 2, 9, 2,11, 9,-1,-1,-1,-1],
    [ 0, 2,11, 8, 0,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 3, 8, 2, 8,10,10, 8, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 9,10, 2, 0, 9, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 3, 8, 2, 8,10, 0, 1, 8, 1,10, 8,-1,-1,-1,-1],
    [ 1,10, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 3, 8, 9, 1, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 9, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 3, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
];

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------------------------------------------------------
    // Determinism test
    // -------------------------------------------------------------------------

    #[test]
    fn test_determinism() {
        // Create a Gaussian blob scalar field
        let n = 15;
        let spacing = 0.5;
        let origin = [-3.5, -3.5, -3.5];
        let mut grid = vec![0.0f64; n * n * n];

        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let x = origin[0] + ix as f64 * spacing;
                    let y = origin[1] + iy as f64 * spacing;
                    let z = origin[2] + iz as f64 * spacing;
                    let r2 = x * x + y * y + z * z;
                    grid[ix * n * n + iy * n + iz] = (-r2).exp();
                }
            }
        }

        let result1 = marching_cubes(&grid, [n, n, n], origin, spacing, 0.3);
        let result2 = marching_cubes(&grid, [n, n, n], origin, spacing, 0.3);

        assert_eq!(
            result1.vertices, result2.vertices,
            "Vertices must be identical"
        );
        assert_eq!(
            result1.indices, result2.indices,
            "Indices must be identical"
        );
        assert_eq!(
            result1.normals, result2.normals,
            "Normals must be identical"
        );
    }

    // -------------------------------------------------------------------------
    // Empty isosurface tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_empty_all_above() {
        // All values above isovalue -> no triangles (entire grid is "inside")
        let grid = vec![1.0f64; 8];
        let result = marching_cubes(&grid, [2, 2, 2], [0.0, 0.0, 0.0], 1.0, 0.5);
        assert!(result.vertices.is_empty());
        assert!(result.indices.is_empty());
        assert!(result.normals.is_empty());
    }

    #[test]
    fn test_empty_all_below() {
        // All values below isovalue -> no triangles
        let grid = vec![0.0f64; 8];
        let result = marching_cubes(&grid, [2, 2, 2], [0.0, 0.0, 0.0], 1.0, 0.5);
        assert!(result.vertices.is_empty());
        assert!(result.indices.is_empty());
        assert!(result.normals.is_empty());
    }

    // -------------------------------------------------------------------------
    // Single cell with one corner above
    // -------------------------------------------------------------------------

    #[test]
    fn test_single_corner() {
        // 2x2x2 grid: one corner above isovalue
        // Corner 0 (0,0,0) = 1.0, all others = 0.0
        // cube_index = 0b00000001 = 1
        let mut grid = vec![0.0f64; 8];
        grid[0] = 1.0; // corner 0 at (0,0,0)

        let result = marching_cubes(&grid, [2, 2, 2], [0.0, 0.0, 0.0], 1.0, 0.5);

        // cube_index=1 -> TRI_TABLE[1] = [0, 8, 3, -1, ...]
        // One triangle from edges 0, 8, 3
        assert_eq!(result.indices.len(), 3, "One triangle = 3 indices");
        assert_eq!(result.vertices.len(), 9, "3 vertices * 3 coords = 9");
        assert_eq!(result.normals.len(), 9, "3 normals * 3 coords = 9");
    }

    // -------------------------------------------------------------------------
    // Sphere isosurface
    // -------------------------------------------------------------------------

    #[test]
    fn test_sphere_isosurface() {
        // Create a Gaussian field centered at origin
        let n = 20;
        let spacing = 0.5;
        let offset = -(n as f64 - 1.0) * spacing / 2.0;
        let origin = [offset, offset, offset];
        let mut grid = vec![0.0f64; n * n * n];

        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let x = origin[0] + ix as f64 * spacing;
                    let y = origin[1] + iy as f64 * spacing;
                    let z = origin[2] + iz as f64 * spacing;
                    let r2 = x * x + y * y + z * z;
                    grid[ix * n * n + iy * n + iz] = (-r2 / 4.0).exp();
                }
            }
        }

        let isovalue = 0.3;
        let result = marching_cubes(&grid, [n, n, n], origin, spacing, isovalue);

        // Should produce a closed surface (non-empty)
        assert!(
            result.vertices.len() > 30,
            "Sphere should have many vertices, got {}",
            result.vertices.len() / 3
        );
        assert!(
            result.indices.len() > 30,
            "Sphere should have many triangles"
        );

        // Verify normals are approximately unit length
        let n_verts = result.vertices.len() / 3;
        for i in 0..n_verts {
            let nx = result.normals[i * 3] as f64;
            let ny = result.normals[i * 3 + 1] as f64;
            let nz = result.normals[i * 3 + 2] as f64;
            let len = (nx * nx + ny * ny + nz * nz).sqrt();
            assert!(
                (len - 1.0).abs() < 0.01,
                "Normal at vertex {} has length {}, expected ~1.0",
                i,
                len
            );
        }

        // For a Gaussian blob, the isosurface should be approximately spherical.
        // Check that all vertices are at roughly the same distance from origin.
        let expected_radius = (-4.0 * isovalue.ln()).sqrt(); // from exp(-r^2/4) = isovalue
        let mut max_deviation = 0.0f64;
        for i in 0..n_verts {
            let vx = result.vertices[i * 3] as f64;
            let vy = result.vertices[i * 3 + 1] as f64;
            let vz = result.vertices[i * 3 + 2] as f64;
            let r = (vx * vx + vy * vy + vz * vz).sqrt();
            let deviation = (r - expected_radius).abs();
            max_deviation = max_deviation.max(deviation);
        }
        // Grid spacing is 0.5, so vertex positions can deviate by up to ~spacing
        assert!(
            max_deviation < spacing * 1.5,
            "Max radius deviation {:.4} exceeds tolerance (expected radius {:.4})",
            max_deviation,
            expected_radius
        );
    }

    // -------------------------------------------------------------------------
    // Dual isosurface
    // -------------------------------------------------------------------------

    #[test]
    fn test_dual_isosurface() {
        // Create a field with both positive and negative values
        // (like an orbital with two lobes)
        let n = 15;
        let spacing = 0.4;
        let offset = -(n as f64 - 1.0) * spacing / 2.0;
        let origin = [offset, offset, offset];
        let mut grid = vec![0.0f64; n * n * n];

        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let x = origin[0] + ix as f64 * spacing;
                    let y = origin[1] + iy as f64 * spacing;
                    let z = origin[2] + iz as f64 * spacing;
                    let r2 = x * x + y * y + z * z;
                    // p-type orbital shape: x * exp(-r^2)
                    grid[ix * n * n + iy * n + iz] = x * (-r2).exp();
                }
            }
        }

        let (positive, negative) = dual_marching_cubes(&grid, [n, n, n], origin, spacing, 0.05);

        // Both lobes should have non-empty meshes
        assert!(
            !positive.vertices.is_empty(),
            "Positive lobe should have vertices"
        );
        assert!(
            !negative.vertices.is_empty(),
            "Negative lobe should have vertices"
        );

        // Positive lobe vertices should have x > 0 (mostly)
        let pos_verts = positive.vertices.len() / 3;
        let pos_x_positive = (0..pos_verts)
            .filter(|&i| positive.vertices[i * 3] > 0.0)
            .count();
        assert!(
            pos_x_positive > pos_verts / 2,
            "Most positive lobe vertices should have x > 0"
        );

        // Negative lobe vertices should have x < 0 (mostly)
        let neg_verts = negative.vertices.len() / 3;
        let neg_x_negative = (0..neg_verts)
            .filter(|&i| negative.vertices[i * 3] < 0.0)
            .count();
        assert!(
            neg_x_negative > neg_verts / 2,
            "Most negative lobe vertices should have x < 0"
        );
    }

    // -------------------------------------------------------------------------
    // Normal direction test
    // -------------------------------------------------------------------------

    #[test]
    fn test_normals_point_outward_for_gaussian() {
        // For a Gaussian blob at origin, normals should point outward
        // (away from the center / away from higher values)
        let n = 12;
        let spacing = 0.6;
        let offset = -(n as f64 - 1.0) * spacing / 2.0;
        let origin = [offset, offset, offset];
        let mut grid = vec![0.0f64; n * n * n];

        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let x = origin[0] + ix as f64 * spacing;
                    let y = origin[1] + iy as f64 * spacing;
                    let z = origin[2] + iz as f64 * spacing;
                    let r2 = x * x + y * y + z * z;
                    grid[ix * n * n + iy * n + iz] = (-r2 / 2.0).exp();
                }
            }
        }

        let result = marching_cubes(&grid, [n, n, n], origin, spacing, 0.3);
        let n_verts = result.vertices.len() / 3;

        // For a centered Gaussian, normals should point roughly radially outward
        // (dot product of normal with position should be positive)
        let mut outward_count = 0;
        for i in 0..n_verts {
            let vx = result.vertices[i * 3] as f64;
            let vy = result.vertices[i * 3 + 1] as f64;
            let vz = result.vertices[i * 3 + 2] as f64;
            let nx = result.normals[i * 3] as f64;
            let ny = result.normals[i * 3 + 1] as f64;
            let nz = result.normals[i * 3 + 2] as f64;
            let dot = vx * nx + vy * ny + vz * nz;
            if dot > 0.0 {
                outward_count += 1;
            }
        }

        // At least 90% of normals should point outward
        let fraction = outward_count as f64 / n_verts as f64;
        assert!(
            fraction > 0.9,
            "Expected >90% outward normals, got {:.1}%",
            fraction * 100.0
        );
    }

    // -------------------------------------------------------------------------
    // Index validity
    // -------------------------------------------------------------------------

    #[test]
    fn test_index_validity() {
        let n = 10;
        let spacing = 0.5;
        let origin = [-2.25, -2.25, -2.25];
        let mut grid = vec![0.0f64; n * n * n];

        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let x = origin[0] + ix as f64 * spacing;
                    let y = origin[1] + iy as f64 * spacing;
                    let z = origin[2] + iz as f64 * spacing;
                    let r2 = x * x + y * y + z * z;
                    grid[ix * n * n + iy * n + iz] = (-r2).exp();
                }
            }
        }

        let result = marching_cubes(&grid, [n, n, n], origin, spacing, 0.2);
        let n_verts = (result.vertices.len() / 3) as u32;

        // All indices must be valid
        for &idx in &result.indices {
            assert!(
                idx < n_verts,
                "Index {} out of bounds (n_verts = {})",
                idx,
                n_verts
            );
        }

        // Number of indices must be divisible by 3 (triangles)
        assert_eq!(
            result.indices.len() % 3,
            0,
            "Indices must form complete triangles"
        );

        // Vertices and normals must have the same length
        assert_eq!(
            result.vertices.len(),
            result.normals.len(),
            "Vertices and normals must have same length"
        );
    }

    // -------------------------------------------------------------------------
    // Result struct serialization
    // -------------------------------------------------------------------------

    #[test]
    fn test_result_serializable() {
        let result = MarchingCubesResult {
            vertices: vec![0.0, 1.0, 2.0],
            indices: vec![0],
            normals: vec![0.0, 0.0, 1.0],
        };
        let json = serde_json::to_string(&result).unwrap();
        let deserialized: MarchingCubesResult = serde_json::from_str(&json).unwrap();
        assert_eq!(result.vertices, deserialized.vertices);
        assert_eq!(result.indices, deserialized.indices);
        assert_eq!(result.normals, deserialized.normals);
    }
}
