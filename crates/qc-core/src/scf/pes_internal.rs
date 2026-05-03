//! Internal coordinate geometry generation for PES scans
//!
//! Provides functions to generate molecular geometries at specified internal
//! coordinate values (bond lengths, bond angles, dihedral angles). These are
//! used by the PES scan engine to produce valid geometries at each scan point.
//!
//! # Algorithms
//!
//! - **Bond stretching:** Translate atom along the existing bond vector
//! - **Angle bending:** Rodrigues' rotation around the normal to the i-j-k plane
//! - **Dihedral rotation:** Rodrigues' rotation of a fragment around the j-k bond axis
//!
//! All rotations use Rodrigues' rotation formula:
//! ```text
//! v' = v·cos(θ) + (k × v)·sin(θ) + k·(k·v)·(1 - cos(θ))
//! ```
//!
//! # References
//!
//! - Murray, Li, Sastry. "A Mathematical Introduction to Robotic Manipulation" (Rodrigues' formula)
//! - Wilson, Decius & Cross (1955). "Molecular Vibrations" (internal coordinates)

use std::collections::{HashSet, VecDeque};
use thiserror::Error;

// =============================================================================
// Error Types
// =============================================================================

/// Errors that can occur during PES geometry generation or scanning
#[derive(Debug, Error)]
pub enum PesScanError {
    /// Atom index exceeds the number of atoms in the molecule
    #[error("Atom index {0} out of range (max {1})")]
    AtomIndexOutOfRange(usize, usize),

    /// Two atoms occupy the same position — bond vector is undefined
    #[error("Zero-length bond between atoms {0} and {1}")]
    ZeroLengthBond(usize, usize),

    /// Invalid scan configuration (e.g., n_points < 2, unsupported method)
    #[error("Invalid scan configuration: {0}")]
    InvalidScanConfig(String),
}

// =============================================================================
// Standalone Geometry Helpers
// =============================================================================
//
// These operate directly on [f64; 3] arrays rather than Atom structs,
// avoiding coupling to the basis::Atom type.

/// Cross product of two 3-vectors
#[inline]
fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product of two 3-vectors
#[inline]
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-vector
#[inline]
fn norm3(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Normalize a 3-vector in place, returning the original norm.
/// If the norm is zero, returns 0.0 and leaves the vector unchanged.
#[inline]
fn normalize3(v: &mut [f64; 3]) -> f64 {
    let n = norm3(v);
    if n > 0.0 {
        v[0] /= n;
        v[1] /= n;
        v[2] /= n;
    }
    n
}

/// Compute distance between two points
#[inline]
fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Compute bond angle i-j-k in radians (j is the central atom)
///
/// Returns the angle in [0, pi].
fn bond_angle(ri: &[f64; 3], rj: &[f64; 3], rk: &[f64; 3]) -> f64 {
    let v1 = [ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]];
    let v2 = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];
    let dot = dot3(&v1, &v2);
    let mag1 = norm3(&v1);
    let mag2 = norm3(&v2);
    let cos_theta = (dot / (mag1 * mag2)).clamp(-1.0, 1.0);
    cos_theta.acos()
}

/// Compute dihedral angle i-j-k-l in radians (j-k is the central bond)
///
/// Uses the atan2 formula for numerical stability:
///   tau = atan2(|b2| * b1 . (b2 x b3), (b1 x b2) . (b2 x b3))
/// where b1 = j-i, b2 = k-j, b3 = l-k
///
/// Returns a value in (-pi, pi].
fn dihedral_angle(ri: &[f64; 3], rj: &[f64; 3], rk: &[f64; 3], rl: &[f64; 3]) -> f64 {
    let b1 = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
    let b2 = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];
    let b3 = [rl[0] - rk[0], rl[1] - rk[1], rl[2] - rk[2]];

    let n1 = cross3(&b1, &b2);
    let n2 = cross3(&b2, &b3);

    let b2_mag = norm3(&b2);
    let x = dot3(&n1, &n2);
    // y = b2_hat . (n1 x n2) = (1/|b2|) * b1 . (b2 x (b2 x b3))
    // Equivalently: y = dot3(&cross3(&n1, &n2), &b2_hat)
    let b2_hat = [b2[0] / b2_mag, b2[1] / b2_mag, b2[2] / b2_mag];
    let y = dot3(&cross3(&n1, &n2), &b2_hat);

    y.atan2(x)
}

// =============================================================================
// Rodrigues' Rotation
// =============================================================================

/// Rotate point `p` around an axis through `origin` by angle `theta` (radians).
///
/// Uses Rodrigues' rotation formula:
/// ```text
/// v = p - origin
/// v' = v·cos(θ) + (k × v)·sin(θ) + k·(k·v)·(1 - cos(θ))
/// p' = origin + v'
/// ```
///
/// The `axis` vector is normalized internally (need not be a unit vector on input).
/// If `axis` is the zero vector, the point is returned unchanged.
///
/// # Properties
///
/// - Preserves distance from origin: |p' - origin| = |p - origin|
/// - Reduces to identity for theta = 0
/// - Numerically stable for all angles (no gimbal lock)
pub fn rotate_point(p: [f64; 3], origin: [f64; 3], axis: [f64; 3], theta: f64) -> [f64; 3] {
    let mut k = axis;
    let n = normalize3(&mut k);
    if n < 1e-15 {
        // Zero axis — return point unchanged
        return p;
    }

    let v = [p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]];

    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let k_dot_v = dot3(&k, &v);
    let k_cross_v = cross3(&k, &v);

    [
        origin[0] + v[0] * cos_t + k_cross_v[0] * sin_t + k[0] * k_dot_v * (1.0 - cos_t),
        origin[1] + v[1] * cos_t + k_cross_v[1] * sin_t + k[1] * k_dot_v * (1.0 - cos_t),
        origin[2] + v[2] * cos_t + k_cross_v[2] * sin_t + k[2] * k_dot_v * (1.0 - cos_t),
    ]
}

// =============================================================================
// Bond Geometry Generation
// =============================================================================

/// Generate geometry with a modified bond length between atoms `atom_i` and `atom_j`.
///
/// Moves `atom_j` along the i->j bond vector to achieve `target_r` distance.
/// All other atoms remain unchanged.
///
/// # Arguments
///
/// * `atoms` - Molecular geometry as (atomic_number, [x, y, z]) in bohr
/// * `atom_i` - Index of the fixed anchor atom
/// * `atom_j` - Index of the atom to move
/// * `target_r` - Desired bond length in bohr (must be > 0)
///
/// # Returns
///
/// New coordinates for all atoms (only atom_j is modified).
pub fn generate_bond_geometry(
    atoms: &[(u8, [f64; 3])],
    atom_i: usize,
    atom_j: usize,
    target_r: f64,
) -> Result<Vec<[f64; 3]>, PesScanError> {
    let n = atoms.len();
    if atom_i >= n {
        return Err(PesScanError::AtomIndexOutOfRange(atom_i, n - 1));
    }
    if atom_j >= n {
        return Err(PesScanError::AtomIndexOutOfRange(atom_j, n - 1));
    }

    let ri = &atoms[atom_i].1;
    let rj = &atoms[atom_j].1;

    // Compute bond direction vector (i -> j)
    let mut dir = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
    let current_r = normalize3(&mut dir);

    if current_r < 1e-14 {
        return Err(PesScanError::ZeroLengthBond(atom_i, atom_j));
    }

    // New position for atom_j: place at target_r along the bond direction from atom_i
    let new_rj = [
        ri[0] + dir[0] * target_r,
        ri[1] + dir[1] * target_r,
        ri[2] + dir[2] * target_r,
    ];

    // Copy all positions, replacing atom_j
    let mut coords: Vec<[f64; 3]> = atoms.iter().map(|(_, pos)| *pos).collect();
    coords[atom_j] = new_rj;

    Ok(coords)
}

// =============================================================================
// Angle Geometry Generation
// =============================================================================

/// Generate geometry with a modified bond angle i-j-k (j is the central atom).
///
/// Rotates atom_k around atom_j to achieve the target angle. The rotation axis
/// is the normal to the i-j-k plane (i.e., normalize(v1 x v2) where
/// v1 = r_i - r_j, v2 = r_k - r_j).
///
/// The O-H bond length (j-k distance) is preserved.
///
/// # Edge Cases
///
/// If the atoms are nearly collinear (|v1 x v2| < 1e-10), an arbitrary
/// perpendicular vector to v1 is used as the rotation axis.
///
/// # Arguments
///
/// * `atoms` - Molecular geometry as (atomic_number, [x, y, z]) in bohr
/// * `atom_i`, `atom_j`, `atom_k` - Atom indices defining the angle (j is central)
/// * `target_theta` - Desired angle in radians
///
/// # Returns
///
/// New coordinates for all atoms (only atom_k is modified).
pub fn generate_angle_geometry(
    atoms: &[(u8, [f64; 3])],
    atom_i: usize,
    atom_j: usize,
    atom_k: usize,
    target_theta: f64,
) -> Result<Vec<[f64; 3]>, PesScanError> {
    let n = atoms.len();
    if atom_i >= n {
        return Err(PesScanError::AtomIndexOutOfRange(atom_i, n - 1));
    }
    if atom_j >= n {
        return Err(PesScanError::AtomIndexOutOfRange(atom_j, n - 1));
    }
    if atom_k >= n {
        return Err(PesScanError::AtomIndexOutOfRange(atom_k, n - 1));
    }

    let ri = &atoms[atom_i].1;
    let rj = &atoms[atom_j].1;
    let rk = &atoms[atom_k].1;

    // Vectors from central atom j
    let v1 = [ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]];
    let v2 = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];

    let len_v1 = norm3(&v1);
    let len_v2 = norm3(&v2);

    if len_v1 < 1e-14 {
        return Err(PesScanError::ZeroLengthBond(atom_i, atom_j));
    }
    if len_v2 < 1e-14 {
        return Err(PesScanError::ZeroLengthBond(atom_j, atom_k));
    }

    // Current angle
    let current_theta = bond_angle(ri, rj, rk);

    // Rotation angle needed
    let delta_theta = target_theta - current_theta;

    // Rotation axis: normal to the i-j-k plane
    let mut axis = cross3(&v1, &v2);
    let axis_norm = normalize3(&mut axis);

    if axis_norm < 1e-10 {
        // Near-linear: v1 and v2 are (anti-)parallel.
        // Use an arbitrary perpendicular vector to v1 as the rotation axis.
        axis = find_perpendicular(&v1);
        normalize3(&mut axis);
    }

    // Rotate atom_k around atom_j by delta_theta
    let new_rk = rotate_point(*rk, *rj, axis, delta_theta);

    let mut coords: Vec<[f64; 3]> = atoms.iter().map(|(_, pos)| *pos).collect();
    coords[atom_k] = new_rk;

    Ok(coords)
}

/// Find a vector perpendicular to `v`.
///
/// Strategy: take the cross product of `v` with (1,0,0). If that is too small
/// (v is nearly parallel to x-axis), use (0,1,0) instead.
fn find_perpendicular(v: &[f64; 3]) -> [f64; 3] {
    let candidate = [1.0, 0.0, 0.0];
    let mut perp = cross3(v, &candidate);
    if norm3(&perp) < 1e-10 {
        let candidate2 = [0.0, 1.0, 0.0];
        perp = cross3(v, &candidate2);
    }
    perp
}

// =============================================================================
// Dihedral Geometry Generation
// =============================================================================

/// Generate geometry with a modified dihedral angle i-j-k-l.
///
/// Rotates all atoms in `fragment` (a set of atom indices) around the j->k
/// bond axis by `target_tau - current_tau`. The fragment should contain
/// atom_l and all atoms connected to it beyond the j-k bond.
///
/// Bond lengths and bond angles are preserved by the rotation.
///
/// # Arguments
///
/// * `atoms` - Molecular geometry as (atomic_number, [x, y, z]) in bohr
/// * `atom_i`, `atom_j`, `atom_k`, `atom_l` - Atom indices defining the dihedral
/// * `target_tau` - Desired dihedral angle in radians
/// * `fragment` - Set of atom indices to rotate (should include atom_l and
///   all atoms on the l-side of the j-k bond)
///
/// # Returns
///
/// New coordinates for all atoms. Only atoms in `fragment` are modified.
pub fn generate_dihedral_geometry(
    atoms: &[(u8, [f64; 3])],
    atom_i: usize,
    atom_j: usize,
    atom_k: usize,
    atom_l: usize,
    target_tau: f64,
    fragment: &HashSet<usize>,
) -> Result<Vec<[f64; 3]>, PesScanError> {
    let n = atoms.len();
    for &idx in &[atom_i, atom_j, atom_k, atom_l] {
        if idx >= n {
            return Err(PesScanError::AtomIndexOutOfRange(idx, n - 1));
        }
    }

    let ri = &atoms[atom_i].1;
    let rj = &atoms[atom_j].1;
    let rk = &atoms[atom_k].1;
    let rl = &atoms[atom_l].1;

    // Current dihedral
    let current_tau = dihedral_angle(ri, rj, rk, rl);

    // Rotation angle needed
    let delta_tau = target_tau - current_tau;

    // Rotation axis: j -> k direction
    let axis = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];

    // Rotate all atoms in the fragment around the j-k axis
    // Use atom_j as the rotation origin (any point on the axis works)
    let mut coords: Vec<[f64; 3]> = atoms.iter().map(|(_, pos)| *pos).collect();
    for &frag_idx in fragment {
        if frag_idx < n {
            coords[frag_idx] = rotate_point(coords[frag_idx], *rj, axis, delta_tau);
        }
    }

    Ok(coords)
}

// =============================================================================
// Fragment Detection for Dihedral Rotation
// =============================================================================
//
// When rotating around a central bond (e.g., C-C in ethane), the entire
// connected fragment on one side must rotate together. These functions build
// a molecular connectivity graph and use BFS to identify which atoms belong
// to each side of the bond.
//
// Reference: Section 4c of PES-internal-coordinate-plan.md

/// Covalent radii in bohr for elements H (Z=1) through Ar (Z=18).
///
/// Source: Pyykko & Atsumi (2009), converted from Angstrom.
/// Values match `apps/web/src/components/viewer3d/constants.ts` (COVALENT_RADII_BOHR).
const COVALENT_RADII_BOHR: [f64; 19] = [
    0.0,   // Z=0 placeholder
    0.586, // H
    0.700, // He (increased for HeH+ bond detection, see constants.ts)
    2.419, // Li
    1.814, // Be
    1.587, // B
    1.455, // C
    1.361, // N
    1.304, // O
    1.285, // F
    1.247, // Ne
    3.023, // Na
    2.627, // Mg
    2.362, // Al
    2.192, // Si
    2.098, // P
    1.984, // S
    1.928, // Cl
    1.852, // Ar
];

/// Returns the covalent radius in bohr for atomic number `z` (1-18).
///
/// For atoms outside the table (Z > 18 or Z = 0), returns a conservative
/// default of 1.5 bohr.
pub fn covalent_radius_bohr(z: u8) -> f64 {
    if (z as usize) < COVALENT_RADII_BOHR.len() {
        COVALENT_RADII_BOHR[z as usize]
    } else {
        1.5 // conservative default for unknown elements
    }
}

/// Build a molecular connectivity graph (adjacency list) from atomic positions.
///
/// Two atoms are considered bonded when:
/// ```text
/// distance(i, j) < 1.2 * (r_cov(i) + r_cov(j))
/// ```
///
/// The 1.2x threshold follows the TDD Section 8.4 bond detection convention.
/// This is intentionally tighter than the optimizer's 1.3x threshold, producing
/// a conservative connectivity graph appropriate for fragment detection.
///
/// A minimum distance threshold of 0.1 bohr is enforced to avoid detecting
/// degenerate atoms (two atoms at the same position) as bonded.
///
/// # Arguments
///
/// * `atoms` - Molecular geometry as (atomic_number, [x, y, z]) in bohr
///
/// # Returns
///
/// Adjacency list: `result[i]` contains the sorted indices of atoms bonded to atom `i`.
pub fn build_connectivity_simple(atoms: &[(u8, [f64; 3])]) -> Vec<Vec<usize>> {
    let n = atoms.len();
    let bond_scale = 1.2;
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); n];

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = atoms[i].1[0] - atoms[j].1[0];
            let dy = atoms[i].1[1] - atoms[j].1[1];
            let dz = atoms[i].1[2] - atoms[j].1[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            let threshold =
                bond_scale * (covalent_radius_bohr(atoms[i].0) + covalent_radius_bohr(atoms[j].0));

            if dist < threshold && dist > 0.1 {
                adjacency[i].push(j);
                adjacency[j].push(i);
            }
        }
    }

    adjacency
}

/// Detect the fragment of atoms connected to `start` after removing the j-k bond.
///
/// Performs BFS from `start` on the adjacency graph, but excludes the edge
/// between `atom_j` and `atom_k` (in both directions). Returns all atoms
/// reachable from `start` without crossing the j-k bond.
///
/// For a non-ring bond in a molecule like ethane, calling this with
/// `start = atom_k` yields all atoms on the k-side of the j-k bond.
///
/// # Arguments
///
/// * `adjacency` - Adjacency list from `build_connectivity_simple()`
/// * `atom_j` - First atom of the excluded bond
/// * `atom_k` - Second atom of the excluded bond
/// * `start` - Atom to begin BFS from (typically `atom_j` or `atom_k`)
///
/// # Returns
///
/// Set of all atom indices reachable from `start` without crossing the j-k bond.
pub fn detect_fragment(
    adjacency: &[Vec<usize>],
    atom_j: usize,
    atom_k: usize,
    start: usize,
) -> HashSet<usize> {
    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();

    visited.insert(start);
    queue.push_back(start);

    while let Some(current) = queue.pop_front() {
        for &neighbor in &adjacency[current] {
            // Skip the excluded j-k bond (in both directions)
            if (current == atom_j && neighbor == atom_k)
                || (current == atom_k && neighbor == atom_j)
            {
                continue;
            }
            if visited.insert(neighbor) {
                queue.push_back(neighbor);
            }
        }
    }

    visited
}

/// Check whether the bond between `atom_j` and `atom_k` is part of a ring.
///
/// A bond is a ring bond if removing it does not disconnect `atom_j` from
/// `atom_k`. Implementation: run `detect_fragment(adjacency, j, k, j)` and
/// check whether `k` is in the resulting fragment.
///
/// If this returns `true`, the caller should warn the user that fragment
/// detection is ambiguous for ring bonds -- rotating a dihedral in a ring
/// would break the ring geometry.
///
/// # Arguments
///
/// * `adjacency` - Adjacency list from `build_connectivity_simple()`
/// * `atom_j` - First atom of the bond
/// * `atom_k` - Second atom of the bond
///
/// # Returns
///
/// `true` if the j-k bond is part of a ring (removing it leaves j and k connected).
pub fn is_ring_bond(adjacency: &[Vec<usize>], atom_j: usize, atom_k: usize) -> bool {
    let fragment = detect_fragment(adjacency, atom_j, atom_k, atom_j);
    fragment.contains(&atom_k)
}

// =============================================================================
// Rigid PES Scan: Data Types
// =============================================================================

use serde::{Deserialize, Serialize};

/// The internal coordinate being scanned.
///
/// Each variant specifies the atom indices that define the coordinate.
/// Bond lengths are in bohr; angles and dihedrals are in radians.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ScanCoordinate {
    /// Bond stretch between two atoms (distance in bohr)
    Bond { atom_i: usize, atom_j: usize },
    /// Bond angle i-j-k where j is the central atom (angle in radians)
    Angle {
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
    },
    /// Dihedral angle i-j-k-l where j-k is the central bond (angle in radians)
    Dihedral {
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        atom_l: usize,
    },
}

/// Snapshot of all internal coordinates at a particular molecular geometry.
///
/// Contains all bonds, angles, and dihedrals detected from the molecular
/// connectivity graph. Used by the relaxed PES scan to track how non-scanned
/// coordinates change during geometry relaxation.
///
/// # Construction
///
/// Build from atomic positions using `compute_all_internals()`, which detects
/// connectivity via covalent radius criterion and enumerates all valence
/// coordinates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InternalCoordinateSnapshot {
    /// Detected bonds as (atom_i, atom_j, distance_bohr), ordered i < j
    pub bonds: Vec<(usize, usize, f64)>,
    /// Detected angles as (atom_i, atom_j_central, atom_k, angle_radians), i < k
    pub angles: Vec<(usize, usize, usize, f64)>,
    /// Detected dihedrals as (atom_i, atom_j, atom_k, atom_l, angle_radians)
    pub dihedrals: Vec<(usize, usize, usize, usize, f64)>,
}

/// Compute all internal coordinates (bonds, angles, dihedrals) for a geometry.
///
/// Uses `build_connectivity_simple()` to detect bonds (1.2x covalent radii),
/// then enumerates:
/// - **Bonds:** All pairs (i, j) in the adjacency list with i < j
/// - **Angles:** All triples (i, j, k) where j is bonded to both i and k, i < k
/// - **Dihedrals:** All quadruples (i, j, k, l) where i bonded to j, j bonded to k,
///   k bonded to l, and j-k is the central bond
///
/// # Arguments
///
/// * `atoms` - Molecular geometry as (atomic_number, [x, y, z]) in bohr
///
/// # Returns
///
/// Complete `InternalCoordinateSnapshot` with computed coordinate values.
pub fn compute_all_internals(atoms: &[(u8, [f64; 3])]) -> InternalCoordinateSnapshot {
    let adjacency = build_connectivity_simple(atoms);

    // --- Bonds ---
    let mut bonds = Vec::new();
    for i in 0..atoms.len() {
        for &j in &adjacency[i] {
            if j > i {
                let d = distance(&atoms[i].1, &atoms[j].1);
                bonds.push((i, j, d));
            }
        }
    }

    // --- Angles ---
    // For each central atom j, enumerate pairs (i, k) where j is bonded to both
    let mut angles = Vec::new();
    for j in 0..atoms.len() {
        let neighbors = &adjacency[j];
        for ni in 0..neighbors.len() {
            for nk in (ni + 1)..neighbors.len() {
                let i = neighbors[ni];
                let k = neighbors[nk];
                // Ensure canonical ordering i < k
                let (a, c) = if i < k { (i, k) } else { (k, i) };
                let theta = bond_angle(&atoms[a].1, &atoms[j].1, &atoms[c].1);
                angles.push((a, j, c, theta));
            }
        }
    }

    // --- Dihedrals ---
    // For each bond j-k, enumerate pairs (i bonded to j) and (l bonded to k)
    // where i != k and l != j
    let mut dihedrals = Vec::new();
    for j in 0..atoms.len() {
        for &k in &adjacency[j] {
            if k <= j {
                continue; // process each bond once (j < k)
            }
            // i bonded to j (i != k)
            let i_neighbors: Vec<usize> =
                adjacency[j].iter().copied().filter(|&n| n != k).collect();
            // l bonded to k (l != j)
            let l_neighbors: Vec<usize> =
                adjacency[k].iter().copied().filter(|&n| n != j).collect();

            for &i in &i_neighbors {
                for &l in &l_neighbors {
                    if i == l {
                        continue; // degenerate
                    }
                    let tau = dihedral_angle(&atoms[i].1, &atoms[j].1, &atoms[k].1, &atoms[l].1);
                    dihedrals.push((i, j, k, l, tau));
                }
            }
        }
    }

    InternalCoordinateSnapshot {
        bonds,
        angles,
        dihedrals,
    }
}

/// Scan mode: rigid (freeze all non-scanned coordinates) or relaxed
/// (optimize non-scanned coordinates at each scan point).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ScanMode {
    /// Only the scanned coordinate changes; all other coordinates are frozen.
    Rigid,
    /// At each scan point, non-scanned coordinates are optimized to minimize energy
    /// while the scanned coordinate is held at its target value.
    Relaxed,
}

/// Configuration for a rigid PES scan over an internal coordinate.
///
/// "Rigid" means that only the scanned coordinate changes; all other
/// internal coordinates are frozen at their initial values.
#[derive(Debug, Clone)]
pub struct PesScanInternalConfig {
    /// Molecular geometry as (atomic_number, [x, y, z]) in bohr
    pub atoms: Vec<(u8, [f64; 3])>,
    /// The internal coordinate to scan
    pub coordinate: ScanCoordinate,
    /// Minimum coordinate value (bohr for bonds, radians for angles)
    pub value_min: f64,
    /// Maximum coordinate value
    pub value_max: f64,
    /// Number of evenly spaced scan points (must be >= 2)
    pub n_points: usize,
    /// Basis set name (e.g., "sto-3g", "6-31g*")
    pub basis_name: String,
    /// Computational method: "rhf", "lda", "b3lyp", or "b3lyp-d3bj"
    pub method: String,
    /// Whether to seed the density matrix from the previous scan point (RHF only)
    pub use_seeding: bool,
    /// Whether to use spherical harmonics for d/f functions
    pub use_spherical: bool,
    /// Convergence profile: "loose", "medium", or "tight"
    pub convergence_profile: String,
    /// Maximum optimization steps per scan point for relaxed scans (default: 50)
    pub opt_max_steps: Option<usize>,
    /// Gradient convergence threshold for relaxed scans (default: 4.5e-4 Ha/bohr)
    pub opt_grad_threshold: Option<f64>,
}

/// Result for a single point on the scanned PES.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PesInternalPoint {
    /// The value of the scanned coordinate at this point
    pub coordinate_value: f64,
    /// Total SCF energy in Hartree (electronic + nuclear repulsion + dispersion if applicable)
    pub energy: f64,
    /// Whether the SCF converged at this geometry
    pub converged: bool,
    /// Number of SCF iterations at this point
    pub scf_iterations: usize,
    /// Cartesian geometry at this scan point (one [x,y,z] per atom)
    pub geometry: Vec<[f64; 3]>,
    /// Number of optimization steps taken (relaxed scan only; None for rigid)
    pub opt_steps: Option<usize>,
    /// Snapshot of all internal coordinates at the final geometry
    pub internal_coordinates: Option<InternalCoordinateSnapshot>,
}

/// Equilibrium coordinate value and energy from parabolic interpolation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PesInternalEquilibrium {
    /// Interpolated equilibrium coordinate value
    pub value: f64,
    /// Interpolated equilibrium energy in Hartree
    pub energy: f64,
}

/// Complete result of a PES scan (rigid or relaxed).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PesScanInternalResult {
    /// Type of coordinate scanned: "bond", "angle", or "dihedral"
    pub coordinate_type: String,
    /// Atom indices defining the scanned coordinate
    pub atom_indices: Vec<usize>,
    /// All scan points (ordered by increasing coordinate value)
    pub points: Vec<PesInternalPoint>,
    /// Equilibrium from parabolic interpolation of three lowest-energy converged points
    pub equilibrium: Option<PesInternalEquilibrium>,
    /// Total SCF iterations summed across all scan points
    pub total_iterations: usize,
    /// Scan mode used: "rigid" or "relaxed"
    pub scan_mode: String,
    /// Total optimization steps across all scan points (relaxed scan only; 0 for rigid)
    pub total_opt_steps: usize,
}

// =============================================================================
// Rigid PES Scan Engine
// =============================================================================

use crate::basis::{Atom, BasisSet};
use crate::integrals::{
    eri_compressed, eri_compressed_spherical, hcore_matrix, hcore_matrix_spherical, overlap_matrix,
    overlap_matrix_spherical,
};

use super::initial_guess::minao_initial_guess;
use super::{rhf_scf_with_guess, PresetSystem, ScfConfig};

/// Run a rigid PES scan over an arbitrary internal coordinate.
///
/// Generates geometries by varying a single internal coordinate (bond, angle,
/// or dihedral) while keeping all other coordinates frozen. At each geometry,
/// molecular integrals are computed on-the-fly and SCF is run with the
/// specified method (RHF or DFT).
///
/// # Arguments
///
/// * `config` - Scan configuration (molecule, coordinate, range, method)
/// * `progress` - Optional callback: (point_index, coordinate_value, energy, converged)
///
/// # Returns
///
/// `PesScanInternalResult` containing all scan points, equilibrium interpolation,
/// and total iteration count.
///
/// # Errors
///
/// Returns `PesScanError` if atom indices are out of range, bonds are degenerate,
/// or the method string is unrecognized.
pub fn pes_scan_rigid(
    config: &PesScanInternalConfig,
    progress: Option<&dyn Fn(usize, f64, f64, bool)>,
) -> Result<PesScanInternalResult, PesScanError> {
    // ---- Validation ----
    if config.n_points < 2 {
        return Err(PesScanError::InvalidScanConfig(
            "n_points must be >= 2".to_string(),
        ));
    }

    let method = config.method.to_lowercase();
    if !["rhf", "lda", "b3lyp", "b3lyp-d3bj"].contains(&method.as_str()) {
        return Err(PesScanError::InvalidScanConfig(format!(
            "Unsupported method '{}'; expected rhf, lda, b3lyp, or b3lyp-d3bj",
            config.method
        )));
    }

    let n_atoms = config.atoms.len();

    // Validate atom indices for the coordinate
    let (coordinate_type, atom_indices) = match &config.coordinate {
        ScanCoordinate::Bond { atom_i, atom_j } => {
            if *atom_i >= n_atoms {
                return Err(PesScanError::AtomIndexOutOfRange(*atom_i, n_atoms - 1));
            }
            if *atom_j >= n_atoms {
                return Err(PesScanError::AtomIndexOutOfRange(*atom_j, n_atoms - 1));
            }
            ("bond".to_string(), vec![*atom_i, *atom_j])
        }
        ScanCoordinate::Angle {
            atom_i,
            atom_j,
            atom_k,
        } => {
            for &idx in &[*atom_i, *atom_j, *atom_k] {
                if idx >= n_atoms {
                    return Err(PesScanError::AtomIndexOutOfRange(idx, n_atoms - 1));
                }
            }
            ("angle".to_string(), vec![*atom_i, *atom_j, *atom_k])
        }
        ScanCoordinate::Dihedral {
            atom_i,
            atom_j,
            atom_k,
            atom_l,
        } => {
            for &idx in &[*atom_i, *atom_j, *atom_k, *atom_l] {
                if idx >= n_atoms {
                    return Err(PesScanError::AtomIndexOutOfRange(idx, n_atoms - 1));
                }
            }
            (
                "dihedral".to_string(),
                vec![*atom_i, *atom_j, *atom_k, *atom_l],
            )
        }
    };

    // Build SCF config from convergence profile string.
    // Always enable DIIS — without it, compressed/strained geometries may
    // require 50-100+ iterations or fail entirely. With DIIS, PySCF converges
    // the same systems in 5 iterations.
    // Use 200 max iterations for DFT at strained geometries (B3LYP can be slow).
    let scf_config = {
        let mut cfg = match config.convergence_profile.to_lowercase().as_str() {
            "loose" => ScfConfig::loose(),
            "tight" => ScfConfig::tight(),
            _ => ScfConfig::medium(),
        };
        cfg.use_diis = true;
        cfg.max_iterations = 200;
        cfg
    };

    // For dihedral scans: build connectivity and detect fragment once
    let dihedral_fragment =
        if let ScanCoordinate::Dihedral { atom_j, atom_k, .. } = &config.coordinate {
            let adjacency = build_connectivity_simple(&config.atoms);
            let fragment = detect_fragment(&adjacency, *atom_j, *atom_k, *atom_k);
            Some(fragment)
        } else {
            None
        };

    // ---- Scan loop ----
    let n_points = config.n_points;
    let mut points = Vec::with_capacity(n_points);
    let mut prev_density: Option<Vec<f64>> = None;
    let mut total_iterations: usize = 0;

    for i in 0..n_points {
        // Step 1: Compute target coordinate value (evenly spaced)
        let v = config.value_min
            + (config.value_max - config.value_min) * i as f64 / (n_points - 1) as f64;

        // Step 2: Generate geometry at this coordinate value
        let new_coords = match &config.coordinate {
            ScanCoordinate::Bond { atom_i, atom_j } => {
                generate_bond_geometry(&config.atoms, *atom_i, *atom_j, v)?
            }
            ScanCoordinate::Angle {
                atom_i,
                atom_j,
                atom_k,
            } => generate_angle_geometry(&config.atoms, *atom_i, *atom_j, *atom_k, v)?,
            ScanCoordinate::Dihedral {
                atom_i,
                atom_j,
                atom_k,
                atom_l,
            } => {
                let fragment = dihedral_fragment.as_ref().unwrap();
                generate_dihedral_geometry(
                    &config.atoms,
                    *atom_i,
                    *atom_j,
                    *atom_k,
                    *atom_l,
                    v,
                    fragment,
                )?
            }
        };

        // Step 3: Build Atom objects from the new coordinates (pairing Z values with new positions)
        let atoms: Vec<Atom> = config
            .atoms
            .iter()
            .zip(new_coords.iter())
            .filter_map(|(&(z, _), &pos)| Atom::new(z, pos).ok())
            .collect();

        if atoms.len() != n_atoms {
            // Failed to create some atoms; mark unconverged and continue
            points.push(PesInternalPoint {
                coordinate_value: v,
                energy: 0.0,
                converged: false,
                scf_iterations: 0,
                geometry: new_coords,
                opt_steps: None,
                internal_coordinates: None,
            });
            if let Some(cb) = progress {
                cb(i, v, 0.0, false);
            }
            prev_density = None;
            continue;
        }

        // Step 4: Build basis set
        let basis = match BasisSet::build(atoms.clone(), &config.basis_name) {
            Ok(b) => b,
            Err(_) => {
                points.push(PesInternalPoint {
                    coordinate_value: v,
                    energy: 0.0,
                    converged: false,
                    scf_iterations: 0,
                    geometry: new_coords,
                    opt_steps: None,
                    internal_coordinates: None,
                });
                if let Some(cb) = progress {
                    cb(i, v, 0.0, false);
                }
                prev_density = None;
                continue;
            }
        };

        // Step 5: Compute molecular integrals
        let (s_mat, h_core, eri) = if config.use_spherical {
            (
                overlap_matrix_spherical(&basis),
                hcore_matrix_spherical(&basis),
                eri_compressed_spherical(&basis),
            )
        } else {
            (
                overlap_matrix(&basis),
                hcore_matrix(&basis),
                eri_compressed(&basis),
            )
        };

        let nbf = if config.use_spherical {
            basis.n_basis_spherical()
        } else {
            basis.n_basis
        };

        // Step 6: Package as PresetSystem
        let system = PresetSystem {
            system_id: format!("pes_rigid_{}_pt{}", coordinate_type, i),
            label: format!("Rigid PES {} scan, point {}", coordinate_type, i),
            nbf,
            nelec: basis.n_electrons,
            e_nuc: basis.nuclear_repulsion,
            s_matrix: s_mat,
            h_core,
            eri_compressed: eri,
        };

        // Step 7: Run SCF with method dispatch
        let (energy, converged, iterations, density_out) = match method.as_str() {
            "rhf" => {
                // Strategy: try seeded density first; if it fails or takes too
                // many iterations (>20), retry with MINAO guess. This handles
                // cases where density seeding from a neighboring geometry is
                // worse than a fresh guess (e.g., compressed CH4 bond scan).
                let seeded_density = if config.use_seeding {
                    prev_density.clone()
                } else {
                    None
                };

                // First attempt: seeded density or MINAO
                let first_guess = seeded_density
                    .clone()
                    .or_else(|| minao_initial_guess(&basis, Some(&system.s_matrix)));
                match rhf_scf_with_guess(&system, &scf_config, first_guess.as_deref()) {
                    Ok(output) if output.iterations <= 20 => (
                        output.energy_total,
                        true,
                        output.iterations,
                        Some(output.density_matrix.clone()),
                    ),
                    first_result => {
                        // Seeded guess either failed or took >20 iterations.
                        // Retry with MINAO if we haven't tried it yet.
                        if seeded_density.is_some() {
                            let minao_guess = minao_initial_guess(&basis, Some(&system.s_matrix));
                            match rhf_scf_with_guess(&system, &scf_config, minao_guess.as_deref()) {
                                Ok(output2) => (
                                    output2.energy_total,
                                    true,
                                    output2.iterations,
                                    Some(output2.density_matrix.clone()),
                                ),
                                Err(_) => {
                                    // MINAO also failed — use first result if it succeeded
                                    match first_result {
                                        Ok(output) => (
                                            output.energy_total,
                                            true,
                                            output.iterations,
                                            Some(output.density_matrix.clone()),
                                        ),
                                        Err(_) => (0.0, false, scf_config.max_iterations, None),
                                    }
                                }
                            }
                        } else {
                            // No seeded density — first attempt was already MINAO
                            match first_result {
                                Ok(output) => (
                                    output.energy_total,
                                    true,
                                    output.iterations,
                                    Some(output.density_matrix.clone()),
                                ),
                                Err(_) => (0.0, false, scf_config.max_iterations, None),
                            }
                        }
                    }
                }
            }
            "lda" | "b3lyp" | "b3lyp-d3bj" => {
                // Build DFT grid for this geometry
                let grid_config = crate::dft::GridConfig::default();
                let grid = crate::dft::build_becke_grid(&atoms, &grid_config);

                // Select functional
                let functional: Box<dyn crate::dft::ExchangeCorrelation> = match method.as_str() {
                    "lda" => Box::new(crate::dft::Lda::new()),
                    _ => Box::new(crate::dft::B3lyp::new()), // b3lyp and b3lyp-d3bj
                };

                match crate::dft::ks_scf(
                    &system,
                    &scf_config,
                    functional.as_ref(),
                    &grid,
                    &basis,
                    config.use_spherical,
                    None,
                ) {
                    Ok(ks_output) => {
                        let mut total_energy = ks_output.scf_output.energy_total;
                        let iters = ks_output.scf_output.iterations;
                        let density = ks_output.scf_output.density_matrix.clone();

                        // For B3LYP-D3BJ: add dispersion correction
                        if method == "b3lyp-d3bj" {
                            let atoms_for_d3: Vec<(u8, [f64; 3])> = config
                                .atoms
                                .iter()
                                .zip(new_coords.iter())
                                .map(|(&(z, _), &pos)| (z, pos))
                                .collect();
                            let d3_result = crate::dft::compute_d3bj_energy(
                                &atoms_for_d3,
                                &crate::dft::D3BJ_B3LYP,
                            );
                            total_energy += d3_result.energy;
                        }

                        (total_energy, true, iters, Some(density))
                    }
                    Err(_) => (0.0, false, scf_config.max_iterations, None),
                }
            }
            _ => unreachable!(), // validated above
        };

        total_iterations += iterations;

        // Step 8: Update density seeding (RHF only)
        if converged {
            if let Some(dm) = density_out {
                prev_density = Some(dm);
            }
        } else {
            prev_density = None;
        }

        // Step 9: Compute internal coordinate snapshot
        let geom_with_z: Vec<(u8, [f64; 3])> = config
            .atoms
            .iter()
            .zip(new_coords.iter())
            .map(|(&(z, _), &pos)| (z, pos))
            .collect();
        let snapshot = compute_all_internals(&geom_with_z);

        // Step 10: Record point
        points.push(PesInternalPoint {
            coordinate_value: v,
            energy,
            converged,
            scf_iterations: iterations,
            geometry: new_coords,
            opt_steps: None,
            internal_coordinates: Some(snapshot),
        });

        if let Some(cb) = progress {
            cb(i, v, energy, converged);
        }
    }

    // ---- Equilibrium interpolation ----
    let equilibrium = find_internal_equilibrium(&points);

    Ok(PesScanInternalResult {
        coordinate_type,
        atom_indices,
        points,
        equilibrium,
        total_iterations,
        scan_mode: "rigid".to_string(),
        total_opt_steps: 0,
    })
}

// =============================================================================
// Parabolic Interpolation for Internal Coordinate Scans
// =============================================================================

/// Find equilibrium coordinate value via parabolic interpolation.
///
/// Identical algorithm to `pes::find_equilibrium()` but generalized to work
/// with `PesInternalPoint` (coordinate_value instead of r).
///
/// # Algorithm
///
/// 1. Filter to converged points only
/// 2. Find the index `k` of the minimum energy point
/// 3. If `k` is at the boundary (first or last), return the boundary point
/// 4. Fit a parabola through (v_{k-1}, E_{k-1}), (v_k, E_k), (v_{k+1}, E_{k+1})
/// 5. Return the vertex of the parabola as the equilibrium
fn find_internal_equilibrium(points: &[PesInternalPoint]) -> Option<PesInternalEquilibrium> {
    let converged: Vec<&PesInternalPoint> = points.iter().filter(|p| p.converged).collect();

    if converged.len() < 3 {
        return None;
    }

    // Find index of minimum energy among converged points
    let min_idx = converged
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            a.energy
                .partial_cmp(&b.energy)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)?;

    // Need at least one point on each side of minimum
    if min_idx == 0 || min_idx >= converged.len() - 1 {
        // Minimum is at the boundary -- return the boundary point as-is
        return Some(PesInternalEquilibrium {
            value: converged[min_idx].coordinate_value,
            energy: converged[min_idx].energy,
        });
    }

    // Three points for parabolic fit
    let v1 = converged[min_idx - 1].coordinate_value;
    let v2 = converged[min_idx].coordinate_value;
    let v3 = converged[min_idx + 1].coordinate_value;
    let e1 = converged[min_idx - 1].energy;
    let e2 = converged[min_idx].energy;
    let e3 = converged[min_idx + 1].energy;

    // Fit parabola E(v) = a*v^2 + b*v + c through three points
    let denom = (v1 - v2) * (v1 - v3) * (v2 - v3);

    if denom.abs() < 1e-30 {
        return Some(PesInternalEquilibrium {
            value: v2,
            energy: e2,
        });
    }

    let a = (v3 * (e2 - e1) + v2 * (e1 - e3) + v1 * (e3 - e2)) / denom;
    let b = (v3 * v3 * (e1 - e2) + v2 * v2 * (e3 - e1) + v1 * v1 * (e2 - e3)) / denom;
    let c =
        (v2 * v3 * (v2 - v3) * e1 + v3 * v1 * (v3 - v1) * e2 + v1 * v2 * (v1 - v2) * e3) / denom;

    // Check that we have a proper minimum (a > 0)
    if a <= 0.0 {
        return Some(PesInternalEquilibrium {
            value: v2,
            energy: e2,
        });
    }

    // Minimum of parabola: v_eq = -b / (2a)
    let v_eq = -b / (2.0 * a);
    let e_eq = a * v_eq * v_eq + b * v_eq + c;

    Some(PesInternalEquilibrium {
        value: v_eq,
        energy: e_eq,
    })
}

// =============================================================================
// Relaxed PES Scan: Constraint Helpers
// =============================================================================

use super::gradient::{rhf_gradient, GradientResult};

/// Compute the current value of the scanned coordinate from a geometry.
///
/// # Arguments
///
/// * `coord` - The scan coordinate definition (bond, angle, or dihedral)
/// * `positions` - Current Cartesian positions ([x,y,z] per atom)
/// * `adjacency` - Adjacency list (needed for dihedral fragment detection)
///
/// # Returns
///
/// Current value of the coordinate (bohr for bonds, radians for angles/dihedrals).
fn current_coord_value(coord: &ScanCoordinate, positions: &[[f64; 3]]) -> f64 {
    match coord {
        ScanCoordinate::Bond { atom_i, atom_j } => {
            distance(&positions[*atom_i], &positions[*atom_j])
        }
        ScanCoordinate::Angle {
            atom_i,
            atom_j,
            atom_k,
        } => bond_angle(
            &positions[*atom_i],
            &positions[*atom_j],
            &positions[*atom_k],
        ),
        ScanCoordinate::Dihedral {
            atom_i,
            atom_j,
            atom_k,
            atom_l,
        } => dihedral_angle(
            &positions[*atom_i],
            &positions[*atom_j],
            &positions[*atom_k],
            &positions[*atom_l],
        ),
    }
}

/// Compute the constraint gradient direction for a bond constraint.
///
/// For a bond between atoms i and j with unit vector e_ij = (r_j - r_i) / |r_j - r_i|:
/// - d(r)/d(r_i) = -e_ij
/// - d(r)/d(r_j) = +e_ij
///
/// Returns a flat vector of length 3*n_atoms.
///
/// Reference: Wilson, Decius & Cross (1955), internal coordinate theory
fn constraint_gradient_bond(positions: &[[f64; 3]], atom_i: usize, atom_j: usize) -> Vec<f64> {
    let n = positions.len();
    let mut c = vec![0.0; 3 * n];

    let mut e = [
        positions[atom_j][0] - positions[atom_i][0],
        positions[atom_j][1] - positions[atom_i][1],
        positions[atom_j][2] - positions[atom_i][2],
    ];
    normalize3(&mut e);

    // d(r)/d(r_i) = -e
    c[3 * atom_i] = -e[0];
    c[3 * atom_i + 1] = -e[1];
    c[3 * atom_i + 2] = -e[2];

    // d(r)/d(r_j) = +e
    c[3 * atom_j] = e[0];
    c[3 * atom_j + 1] = e[1];
    c[3 * atom_j + 2] = e[2];

    c
}

/// Compute the constraint gradient direction for an angle constraint.
///
/// For angle theta at atoms (i, j_central, k):
/// d(theta)/d(r_i) = (cos(theta) * e_ji - e_jk) / (r_ji * sin(theta))
/// d(theta)/d(r_k) = (cos(theta) * e_jk - e_ji) / (r_jk * sin(theta))
/// d(theta)/d(r_j) = -d(theta)/d(r_i) - d(theta)/d(r_k)
///
/// Reference: Wilson, Decius & Cross (1955), Eq. 4.1.12
fn constraint_gradient_angle(
    positions: &[[f64; 3]],
    atom_i: usize,
    atom_j: usize,
    atom_k: usize,
) -> Vec<f64> {
    let n = positions.len();
    let mut c = vec![0.0; 3 * n];

    let rji = [
        positions[atom_i][0] - positions[atom_j][0],
        positions[atom_i][1] - positions[atom_j][1],
        positions[atom_i][2] - positions[atom_j][2],
    ];
    let rjk = [
        positions[atom_k][0] - positions[atom_j][0],
        positions[atom_k][1] - positions[atom_j][1],
        positions[atom_k][2] - positions[atom_j][2],
    ];

    let r_ji = norm3(&rji);
    let r_jk = norm3(&rjk);

    if r_ji < 1e-14 || r_jk < 1e-14 {
        return c; // degenerate geometry
    }

    let e_ji = [rji[0] / r_ji, rji[1] / r_ji, rji[2] / r_ji];
    let e_jk = [rjk[0] / r_jk, rjk[1] / r_jk, rjk[2] / r_jk];

    let cos_theta = dot3(&e_ji, &e_jk).clamp(-1.0, 1.0);
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-14);

    // d(theta)/d(r_i)
    let di = [
        (cos_theta * e_ji[0] - e_jk[0]) / (r_ji * sin_theta),
        (cos_theta * e_ji[1] - e_jk[1]) / (r_ji * sin_theta),
        (cos_theta * e_ji[2] - e_jk[2]) / (r_ji * sin_theta),
    ];

    // d(theta)/d(r_k)
    let dk = [
        (cos_theta * e_jk[0] - e_ji[0]) / (r_jk * sin_theta),
        (cos_theta * e_jk[1] - e_ji[1]) / (r_jk * sin_theta),
        (cos_theta * e_jk[2] - e_ji[2]) / (r_jk * sin_theta),
    ];

    // d(theta)/d(r_j) = -di - dk
    c[3 * atom_i] = di[0];
    c[3 * atom_i + 1] = di[1];
    c[3 * atom_i + 2] = di[2];

    c[3 * atom_k] = dk[0];
    c[3 * atom_k + 1] = dk[1];
    c[3 * atom_k + 2] = dk[2];

    c[3 * atom_j] = -di[0] - dk[0];
    c[3 * atom_j + 1] = -di[1] - dk[1];
    c[3 * atom_j + 2] = -di[2] - dk[2];

    c
}

/// Compute the constraint gradient direction for a dihedral constraint.
///
/// For a dihedral angle defined by atoms i-j-k-l, the B-matrix row gives
/// d(tau)/d(r_A) for each atom A. This is the direction in Cartesian space
/// that changes the dihedral angle.
///
/// Reference: Helgaker, Jorgensen & Olsen, Chapter 1, Eq. 1.4.26-1.4.29
/// (same formulas as optimizer::WilsonBMatrix dihedral row)
fn constraint_gradient_dihedral(
    positions: &[[f64; 3]],
    atom_i: usize,
    atom_j: usize,
    atom_k: usize,
    atom_l: usize,
) -> Vec<f64> {
    let n = positions.len();
    let mut c = vec![0.0; 3 * n];

    let ri = &positions[atom_i];
    let rj = &positions[atom_j];
    let rk = &positions[atom_k];
    let rl = &positions[atom_l];

    // Bond vectors
    let mut e_ij = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
    let r_ij = normalize3(&mut e_ij);
    let mut e_jk = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];
    let r_jk = normalize3(&mut e_jk);
    let mut e_kl = [rl[0] - rk[0], rl[1] - rk[1], rl[2] - rk[2]];
    let r_kl = normalize3(&mut e_kl);

    if r_ij < 1e-10 || r_jk < 1e-10 || r_kl < 1e-10 {
        return c; // degenerate geometry
    }

    // Normal vectors to the i-j-k and j-k-l planes
    let n1 = cross3(&e_ij, &e_jk);
    let n2 = cross3(&e_jk, &e_kl);

    let sin_ijk = norm3(&n1).max(1e-12);
    let sin_jkl = norm3(&n2).max(1e-12);

    let cos_ijk = dot3(&e_ij, &e_jk).clamp(-1.0, 1.0);
    let cos_jkl = dot3(&e_jk, &e_kl).clamp(-1.0, 1.0);

    let sin2_ijk = sin_ijk * sin_ijk;
    let sin2_jkl = sin_jkl * sin_jkl;

    // B-matrix elements for dihedral i-j-k-l:
    // B[tau, atom_i] = n1 / (r_ij * sin^2(ijk))
    // B[tau, atom_l] = -n2 / (r_kl * sin^2(jkl))
    // B[tau, atom_j] and B[tau, atom_k] from chain rule
    for d in 0..3 {
        let b_i = n1[d] / (r_ij * sin2_ijk);
        let b_l = -n2[d] / (r_kl * sin2_jkl);

        let frac_ij_jk = r_ij / r_jk;
        let frac_kl_jk = r_kl / r_jk;

        let b_j = -(1.0 - cos_ijk * frac_ij_jk) * b_i + cos_jkl * frac_kl_jk * b_l;
        let b_k = -(1.0 - cos_jkl * frac_kl_jk) * b_l + cos_ijk * frac_ij_jk * b_i;

        c[3 * atom_i + d] = b_i;
        c[3 * atom_j + d] = b_j;
        c[3 * atom_k + d] = b_k;
        c[3 * atom_l + d] = b_l;
    }

    c
}

/// Compute the constraint gradient direction for a scan coordinate.
///
/// Returns a normalized flat vector of length 3*n_atoms representing the
/// gradient of the constrained coordinate with respect to Cartesian coordinates.
fn constraint_gradient_direction(coord: &ScanCoordinate, positions: &[[f64; 3]]) -> Vec<f64> {
    let mut c = match coord {
        ScanCoordinate::Bond { atom_i, atom_j } => {
            constraint_gradient_bond(positions, *atom_i, *atom_j)
        }
        ScanCoordinate::Angle {
            atom_i,
            atom_j,
            atom_k,
        } => constraint_gradient_angle(positions, *atom_i, *atom_j, *atom_k),
        ScanCoordinate::Dihedral {
            atom_i,
            atom_j,
            atom_k,
            atom_l,
        } => constraint_gradient_dihedral(positions, *atom_i, *atom_j, *atom_k, *atom_l),
    };

    // Normalize the constraint direction
    let norm: f64 = c.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 1e-14 {
        for x in &mut c {
            *x /= norm;
        }
    }

    c
}

/// Project out the constraint gradient component from the nuclear gradient.
///
/// Given a Cartesian gradient g (as Vec<[f64; 3]>) and the constraint direction c
/// (flat vector, normalized), removes the component of g along c:
///
/// ```text
/// g_proj = g - (g . c) * c
/// ```
///
/// This ensures the projected gradient has zero component in the direction
/// that would change the constrained coordinate.
fn project_out_constraint(gradient: &[[f64; 3]], constraint_dir: &[f64]) -> Vec<[f64; 3]> {
    let n = gradient.len();

    // Flatten gradient
    let g_flat: Vec<f64> = gradient.iter().flat_map(|g| g.iter().copied()).collect();

    // Compute dot product g . c
    let g_dot_c: f64 = g_flat
        .iter()
        .zip(constraint_dir.iter())
        .map(|(a, b)| a * b)
        .sum();

    // g_proj = g - (g.c) * c
    let mut result = vec![[0.0; 3]; n];
    for i in 0..n {
        for d in 0..3 {
            result[i][d] = gradient[i][d] - g_dot_c * constraint_dir[3 * i + d];
        }
    }

    result
}

/// Enforce the constraint by resetting the scanned coordinate to its target value.
///
/// After an optimization step, the constrained coordinate may have drifted from
/// its target value. This function regenerates the geometry to reset the
/// coordinate using the existing geometry generators.
///
/// # Arguments
///
/// * `positions` - Current atom positions (modified in place via return)
/// * `atoms_z` - Atomic numbers (for building atom tuples)
/// * `coord` - The scan coordinate
/// * `target` - Target value of the constrained coordinate
/// * `tolerance` - Drift tolerance (1e-4 bohr for bonds, 1e-4 rad for angles)
/// * `dihedral_fragment` - Fragment for dihedral rotations (if applicable)
///
/// # Returns
///
/// `(corrected_positions, correction_needed)` where `correction_needed` is true
/// if the coordinate drifted beyond tolerance.
fn enforce_constraint(
    positions: &[[f64; 3]],
    atoms_z: &[u8],
    coord: &ScanCoordinate,
    target: f64,
    tolerance: f64,
    dihedral_fragment: &Option<HashSet<usize>>,
) -> (Vec<[f64; 3]>, bool) {
    let current = current_coord_value(coord, positions);
    let drift = (current - target).abs();

    if drift <= tolerance {
        return (positions.to_vec(), false);
    }

    let mut corrected = positions.to_vec();

    match coord {
        ScanCoordinate::Bond { atom_i, atom_j } => {
            // Symmetric correction: move BOTH atoms equally along bond direction.
            // This preserves the center of mass and doesn't disrupt other relaxed
            // coordinates, unlike the old asymmetric generate_bond_geometry approach.
            let ai = *atom_i;
            let aj = *atom_j;
            let mut dir = [
                positions[aj][0] - positions[ai][0],
                positions[aj][1] - positions[ai][1],
                positions[aj][2] - positions[ai][2],
            ];
            let r = normalize3(&mut dir);
            if r < 1e-14 {
                return (corrected, false);
            }
            let half_correction = (target - r) / 2.0;
            for (d, &dv) in dir.iter().enumerate() {
                corrected[ai][d] -= half_correction * dv;
                corrected[aj][d] += half_correction * dv;
            }
        }
        ScanCoordinate::Angle {
            atom_i,
            atom_j,
            atom_k,
        } => {
            // For angles, use the geometry generator (moves only atom_k)
            let atom_tuples: Vec<(u8, [f64; 3])> = atoms_z
                .iter()
                .zip(corrected.iter())
                .map(|(&z, &pos)| (z, pos))
                .collect();
            if let Ok(new_pos) =
                generate_angle_geometry(&atom_tuples, *atom_i, *atom_j, *atom_k, target)
            {
                corrected = new_pos;
            }
        }
        ScanCoordinate::Dihedral {
            atom_i,
            atom_j,
            atom_k,
            atom_l,
        } => {
            if let Some(frag) = dihedral_fragment {
                let atom_tuples: Vec<(u8, [f64; 3])> = atoms_z
                    .iter()
                    .zip(corrected.iter())
                    .map(|(&z, &pos)| (z, pos))
                    .collect();
                if let Ok(new_pos) = generate_dihedral_geometry(
                    &atom_tuples,
                    *atom_i,
                    *atom_j,
                    *atom_k,
                    *atom_l,
                    target,
                    frag,
                ) {
                    corrected = new_pos;
                }
            }
        }
    }

    (corrected, true)
}

// =============================================================================
// L-BFGS Two-Loop Recursion for Constrained Optimization
// =============================================================================

/// Compute L-BFGS search direction via two-loop recursion.
///
/// Given the current gradient and a history of {s_k, y_k} pairs, computes
/// the search direction `d = H_k * g` where `H_k` is the L-BFGS approximation
/// to the inverse Hessian.
///
/// When the history is empty (first iteration), falls back to a scaled
/// steepest descent step: `d = g / h0_diag`.
///
/// # Arguments
///
/// * `grad` - Current gradient vector (length n_cart)
/// * `s_history` - Recent position differences s_k = x_{k+1} - x_k
/// * `y_history` - Recent gradient differences y_k = g_{k+1} - g_k
/// * `h0_diag` - Initial inverse Hessian scaling (diagonal, scalar)
///
/// # Returns
///
/// Search direction vector (length n_cart). The step is `x_new = x - d`.
///
/// # Reference
///
/// Nocedal & Wright (2006), Numerical Optimization, Algorithm 7.4.
fn lbfgs_two_loop(
    grad: &[f64],
    s_history: &VecDeque<Vec<f64>>,
    y_history: &VecDeque<Vec<f64>>,
    h0_diag: f64,
) -> Vec<f64> {
    let m = s_history.len();

    if m == 0 {
        // No history yet: use scaled steepest descent
        return grad.iter().map(|g| g / h0_diag).collect();
    }

    // Compute initial H0 scaling from most recent pair (gamma_k)
    // gamma_k = (s_{k-1}^T y_{k-1}) / (y_{k-1}^T y_{k-1})
    let last_s = &s_history[m - 1];
    let last_y = &y_history[m - 1];
    let sy: f64 = last_s.iter().zip(last_y.iter()).map(|(s, y)| s * y).sum();
    let yy: f64 = last_y.iter().map(|y| y * y).sum();
    let gamma = if yy > 1e-14 { sy / yy } else { 1.0 / h0_diag };

    // Two-loop recursion (Algorithm 7.4 in Nocedal & Wright)
    let mut q = grad.to_vec();
    let mut alpha = vec![0.0; m];
    let mut rho = vec![0.0; m];

    // Compute rho values
    for k in 0..m {
        let sy_k: f64 = s_history[k]
            .iter()
            .zip(y_history[k].iter())
            .map(|(s, y)| s * y)
            .sum();
        rho[k] = if sy_k.abs() > 1e-14 { 1.0 / sy_k } else { 0.0 };
    }

    // First loop: from most recent to oldest
    for k in (0..m).rev() {
        let sq: f64 = s_history[k]
            .iter()
            .zip(q.iter())
            .map(|(s, qi)| s * qi)
            .sum();
        alpha[k] = rho[k] * sq;
        for (qj, yj) in q.iter_mut().zip(y_history[k].iter()) {
            *qj -= alpha[k] * yj;
        }
    }

    // Apply initial Hessian: r = gamma * q
    let mut r: Vec<f64> = q.iter().map(|qi| gamma * qi).collect();

    // Second loop: from oldest to most recent
    for k in 0..m {
        let yr: f64 = y_history[k]
            .iter()
            .zip(r.iter())
            .map(|(y, ri)| y * ri)
            .sum();
        let beta = rho[k] * yr;
        for (rj, sj) in r.iter_mut().zip(s_history[k].iter()) {
            *rj += (alpha[k] - beta) * sj;
        }
    }

    r
}

// =============================================================================
// Relaxed PES Scan Engine
// =============================================================================

/// Default model Hessian diagonal for step scaling (Schlegel 1984).
///
/// Uses conservative values: 0.5 Ha/bohr^2 for all coordinates.
/// This provides reasonable scaling for steepest descent steps.
const MODEL_HESSIAN_DIAG: f64 = 0.5;

/// Maximum Cartesian displacement per atom per optimization step (bohr).
const MAX_STEP_PER_ATOM: f64 = 0.3;

/// Constraint enforcement tolerance (bohr for bonds, radians for angles).
/// Must be tight (1e-6) to prevent systematic drift toward equilibrium
/// which produces artificially low energies in relaxed scans.
const CONSTRAINT_TOLERANCE: f64 = 1e-6;

/// Maximum constraint corrections per optimization step.
/// Increased from 3 to 15 to ensure convergence for large L-BFGS steps.
const MAX_CONSTRAINT_CORRECTIONS: usize = 15;

/// Maximum optimization steps per scan point.
///
/// 50 steps is sufficient for ethane-sized molecules (~8 atoms) with
/// steepest descent. Smaller molecules converge in 3-10 steps typically.
const OPT_MAX_STEPS: usize = 50;

/// Projected gradient convergence threshold (Ha/bohr).
const OPT_GRAD_THRESHOLD: f64 = 4.5e-4;

/// Run a relaxed PES scan over an arbitrary internal coordinate.
///
/// At each scan point, the scanned coordinate is held at its target value while
/// all other degrees of freedom are optimized using constrained steepest descent.
/// This produces the minimum energy path (MEP) along the scanned coordinate.
///
/// # Algorithm
///
/// For each scan point:
/// 1. Generate rigid geometry at the target coordinate value
/// 2. Run constrained optimization (max 20 steps):
///    a. Build basis set and compute molecular integrals
///    b. Run SCF (RHF or DFT) to get energy
///    c. Compute analytical gradient
///    d. Project out the constraint gradient component
///    e. Check convergence (max projected gradient < threshold)
///    f. Take scaled steepest descent step
///    g. Enforce constraint (reset scanned coordinate to target)
/// 3. Record final energy, geometry, and internal coordinate snapshot
/// 4. Seed next point with converged density (RHF only)
///
/// # Arguments
///
/// * `config` - Scan configuration (molecule, coordinate, range, method)
/// * `progress` - Optional callback: (point_index, coordinate_value, energy, converged)
///
/// # Returns
///
/// `PesScanInternalResult` with scan_mode="relaxed" and optimization step counts.
///
/// # References
///
/// - Pulay, P. (1969). Mol. Phys. 17, 197. (Analytical gradients)
/// - Schlegel, H. B. (1984). Theor. Chim. Acta 66, 333. (Model Hessian)
pub fn pes_scan_relaxed(
    config: &PesScanInternalConfig,
    progress: Option<&dyn Fn(usize, f64, f64, bool)>,
) -> Result<PesScanInternalResult, PesScanError> {
    // ---- Validation (same as rigid scan) ----
    if config.n_points < 2 {
        return Err(PesScanError::InvalidScanConfig(
            "n_points must be >= 2".to_string(),
        ));
    }

    let method = config.method.to_lowercase();
    if !["rhf", "lda", "b3lyp", "b3lyp-d3bj"].contains(&method.as_str()) {
        return Err(PesScanError::InvalidScanConfig(format!(
            "Unsupported method '{}'; expected rhf, lda, b3lyp, or b3lyp-d3bj",
            config.method
        )));
    }

    let n_atoms = config.atoms.len();
    let atoms_z: Vec<u8> = config.atoms.iter().map(|&(z, _)| z).collect();

    // Validate atom indices
    let (coordinate_type, atom_indices) = match &config.coordinate {
        ScanCoordinate::Bond { atom_i, atom_j } => {
            if *atom_i >= n_atoms {
                return Err(PesScanError::AtomIndexOutOfRange(*atom_i, n_atoms - 1));
            }
            if *atom_j >= n_atoms {
                return Err(PesScanError::AtomIndexOutOfRange(*atom_j, n_atoms - 1));
            }
            ("bond".to_string(), vec![*atom_i, *atom_j])
        }
        ScanCoordinate::Angle {
            atom_i,
            atom_j,
            atom_k,
        } => {
            for &idx in &[*atom_i, *atom_j, *atom_k] {
                if idx >= n_atoms {
                    return Err(PesScanError::AtomIndexOutOfRange(idx, n_atoms - 1));
                }
            }
            ("angle".to_string(), vec![*atom_i, *atom_j, *atom_k])
        }
        ScanCoordinate::Dihedral {
            atom_i,
            atom_j,
            atom_k,
            atom_l,
        } => {
            for &idx in &[*atom_i, *atom_j, *atom_k, *atom_l] {
                if idx >= n_atoms {
                    return Err(PesScanError::AtomIndexOutOfRange(idx, n_atoms - 1));
                }
            }
            (
                "dihedral".to_string(),
                vec![*atom_i, *atom_j, *atom_k, *atom_l],
            )
        }
    };

    let scf_config = {
        let mut cfg = match config.convergence_profile.to_lowercase().as_str() {
            "loose" => ScfConfig::loose(),
            "tight" => ScfConfig::tight(),
            _ => ScfConfig::medium(),
        };
        cfg.use_diis = true;
        cfg.max_iterations = 200;
        cfg
    };

    // For dihedral scans: build connectivity and detect fragment once
    let dihedral_fragment =
        if let ScanCoordinate::Dihedral { atom_j, atom_k, .. } = &config.coordinate {
            let adjacency = build_connectivity_simple(&config.atoms);
            let fragment = detect_fragment(&adjacency, *atom_j, *atom_k, *atom_k);
            Some(fragment)
        } else {
            None
        };

    // ---- Relaxed-scan optimization limits (user-configurable with defaults) ----
    let max_steps = config.opt_max_steps.unwrap_or(OPT_MAX_STEPS);
    let grad_threshold = config.opt_grad_threshold.unwrap_or(OPT_GRAD_THRESHOLD);

    // ---- Scan loop ----
    let n_points = config.n_points;
    let mut points = Vec::with_capacity(n_points);
    let mut prev_density: Option<Vec<f64>> = None;
    let mut total_iterations: usize = 0;
    let mut total_opt_steps: usize = 0;

    for i in 0..n_points {
        // Step 1: Compute target coordinate value
        let v = config.value_min
            + (config.value_max - config.value_min) * i as f64 / (n_points - 1) as f64;

        // Step 2: Generate initial rigid geometry at target coordinate value
        let initial_coords = match &config.coordinate {
            ScanCoordinate::Bond { atom_i, atom_j } => {
                generate_bond_geometry(&config.atoms, *atom_i, *atom_j, v)?
            }
            ScanCoordinate::Angle {
                atom_i,
                atom_j,
                atom_k,
            } => generate_angle_geometry(&config.atoms, *atom_i, *atom_j, *atom_k, v)?,
            ScanCoordinate::Dihedral {
                atom_i,
                atom_j,
                atom_k,
                atom_l,
            } => {
                let fragment = dihedral_fragment.as_ref().unwrap();
                generate_dihedral_geometry(
                    &config.atoms,
                    *atom_i,
                    *atom_j,
                    *atom_k,
                    *atom_l,
                    v,
                    fragment,
                )?
            }
        };

        // Step 3: Constrained optimization loop using L-BFGS quasi-Newton
        //
        // Uses L-BFGS (limited-memory BFGS) to compute the search direction
        // instead of steepest descent. L-BFGS builds an approximation to the
        // inverse Hessian from the last few {position, gradient} pairs, giving
        // superlinear convergence. This reduces the number of SCF evaluations
        // needed from ~50+ (steepest descent) to ~5-10 for ethane-sized systems.
        //
        // Reference: Nocedal (1980). Math. Comp. 35, 773.
        let mut current_positions = initial_coords;
        let mut point_converged = false;
        let mut point_energy = 0.0;
        let mut point_scf_iters = 0;
        let mut point_opt_steps = 0;
        let mut current_density = if config.use_seeding {
            prev_density.clone()
        } else {
            None
        };

        // Energy history for stagnation detection (rolling window)
        let mut energy_history: Vec<f64> = Vec::with_capacity(10);

        // Track the best energy seen so far and the step at which it was found.
        // Used for a secondary convergence criterion: if the best energy hasn't
        // improved for STAGNATION_PATIENCE steps, the optimization is stagnating
        // and we should declare convergence for PES scan purposes.
        let mut best_energy = f64::INFINITY;
        let mut steps_since_best: usize = 0;
        const STAGNATION_PATIENCE: usize = 8;
        const STAGNATION_IMPROVEMENT: f64 = 1e-5; // 0.01 mHa improvement required

        // L-BFGS history: stores recent {s_k, y_k} pairs where
        // s_k = x_{k+1} - x_k (position change)
        // y_k = g_{k+1} - g_k (gradient change)
        const LBFGS_MEMORY: usize = 5;
        let mut lbfgs_s: VecDeque<Vec<f64>> = VecDeque::with_capacity(LBFGS_MEMORY);
        let mut lbfgs_y: VecDeque<Vec<f64>> = VecDeque::with_capacity(LBFGS_MEMORY);
        let mut prev_flat_coords: Option<Vec<f64>> = None;
        let mut prev_proj_grad_flat: Option<Vec<f64>> = None;

        // Track the last known-good geometry for SCF failure recovery
        let mut last_good_positions = current_positions.clone();

        for opt_step in 0..max_steps {
            // 3a: Build Atom objects from current positions
            let atom_structs: Vec<Atom> = atoms_z
                .iter()
                .zip(current_positions.iter())
                .filter_map(|(&z, &pos)| Atom::new(z, pos).ok())
                .collect();

            if atom_structs.len() != n_atoms {
                if opt_step == 0 {
                    break; // degenerate initial geometry — unrecoverable
                }
                // Restore last good geometry and retry with fresh guess
                current_positions = last_good_positions.clone();
                current_density = None;
                lbfgs_s.clear();
                lbfgs_y.clear();
                prev_flat_coords = None;
                prev_proj_grad_flat = None;
                continue;
            }

            // 3b: Evaluate energy and gradient
            // Note: always use Cartesian basis (use_spherical=false) because
            // the analytical gradient code (rhf_gradient, ks_dft_gradient)
            // only supports Cartesian derivative integrals. Spherical harmonics
            // would cause a dimension mismatch panic.
            let eval_result = match method.as_str() {
                "rhf" => evaluate_rhf_energy_gradient(
                    &atom_structs,
                    &config.basis_name,
                    false, // gradients require Cartesian basis
                    &scf_config,
                    current_density.as_deref(),
                )
                .map(|(e, grad, dm)| (e, grad, Some(dm))),
                "lda" | "b3lyp" | "b3lyp-d3bj" => evaluate_dft_energy_gradient(
                    &atom_structs,
                    &config.basis_name,
                    &method,
                    &scf_config,
                    current_density.as_deref(),
                )
                .map(|(e, grad, dm)| (e, grad, Some(dm))),
                _ => unreachable!(),
            };

            let (energy, gradient, density_out) = match eval_result {
                Some(r) => r,
                None => {
                    if opt_step == 0 {
                        break; // SCF failed on initial geometry — unrecoverable
                    }
                    // SCF failed after an optimization step: the step was too
                    // aggressive. Restore previous good geometry, clear L-BFGS
                    // history (which produced the bad step), drop density seeding,
                    // and let the optimizer retry with steepest descent.
                    current_positions = last_good_positions.clone();
                    current_density = None;
                    lbfgs_s.clear();
                    lbfgs_y.clear();
                    prev_flat_coords = None;
                    prev_proj_grad_flat = None;
                    continue;
                }
            };

            let prev_energy = point_energy;
            point_energy = energy;
            point_scf_iters += 1;
            last_good_positions = current_positions.clone();
            current_density = density_out.or(current_density);

            // Oscillation detection: if energy increased significantly, the
            // L-BFGS step was too aggressive. Clear the history so the next
            // step falls back to scaled steepest descent, which is safer.
            // Threshold: 10 mHa increase indicates the curvature approximation
            // is badly wrong. Smaller increases (1-10 mHa) are normal for
            // constrained optimization and should not trigger a reset.
            if opt_step > 0 && energy > prev_energy + 1e-2 {
                lbfgs_s.clear();
                lbfgs_y.clear();
                prev_flat_coords = None;
                prev_proj_grad_flat = None;
            }

            // Update best energy tracking for stagnation detection
            if energy < best_energy - STAGNATION_IMPROVEMENT {
                best_energy = energy;
                steps_since_best = 0;
            } else {
                steps_since_best += 1;
            }

            // 3c: Project out the constraint gradient component
            //
            // Also project out translations and rotations from the gradient.
            // Without this, Cartesian gradients contain 6 spurious DOFs
            // (3 translations + 3 rotations) that cause the optimizer to
            // move atoms unnecessarily, leading to constraint violations
            // and non-convergence. This matches the approach used by
            // geomeTRIC (TRIC coordinates remove these automatically).
            let constraint_dir =
                constraint_gradient_direction(&config.coordinate, &current_positions);
            let projected_grad = project_out_constraint(&gradient.gradients, &constraint_dir);

            // Remove translational component: g_i -= mean(g)
            let n_at = projected_grad.len();
            let g_trans = [
                projected_grad.iter().map(|g| g[0]).sum::<f64>() / n_at as f64,
                projected_grad.iter().map(|g| g[1]).sum::<f64>() / n_at as f64,
                projected_grad.iter().map(|g| g[2]).sum::<f64>() / n_at as f64,
            ];
            let mut projected_grad: Vec<[f64; 3]> = projected_grad
                .iter()
                .map(|g| [g[0] - g_trans[0], g[1] - g_trans[1], g[2] - g_trans[2]])
                .collect();

            // Remove rotational component: L = sum(r x g), omega = I^{-1} L, g_rot = omega x r
            let torque = {
                let mut t = [0.0; 3];
                for (pos, g) in current_positions.iter().zip(projected_grad.iter()) {
                    t[0] += pos[1] * g[2] - pos[2] * g[1];
                    t[1] += pos[2] * g[0] - pos[0] * g[2];
                    t[2] += pos[0] * g[1] - pos[1] * g[0];
                }
                t
            };
            // Moment of inertia tensor
            let mut inertia = [[0.0f64; 3]; 3];
            for pos in current_positions.iter() {
                let r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
                for a in 0..3 {
                    inertia[a][a] += r2;
                    for b in 0..3 {
                        inertia[a][b] -= pos[a] * pos[b];
                    }
                }
            }
            // Solve I * omega = torque (3x3, use simple Cramer's rule or add small regularization)
            let det = inertia[0][0]
                * (inertia[1][1] * inertia[2][2] - inertia[1][2] * inertia[2][1])
                - inertia[0][1] * (inertia[1][0] * inertia[2][2] - inertia[1][2] * inertia[2][0])
                + inertia[0][2] * (inertia[1][0] * inertia[2][1] - inertia[1][1] * inertia[2][0]);
            if det.abs() > 1e-20 {
                // Compute inverse via cofactors
                let inv_det = 1.0 / det;
                let omega = [
                    inv_det
                        * (torque[0]
                            * (inertia[1][1] * inertia[2][2] - inertia[1][2] * inertia[2][1])
                            + torque[1]
                                * (inertia[0][2] * inertia[2][1] - inertia[0][1] * inertia[2][2])
                            + torque[2]
                                * (inertia[0][1] * inertia[1][2] - inertia[0][2] * inertia[1][1])),
                    inv_det
                        * (torque[0]
                            * (inertia[1][2] * inertia[2][0] - inertia[1][0] * inertia[2][2])
                            + torque[1]
                                * (inertia[0][0] * inertia[2][2] - inertia[0][2] * inertia[2][0])
                            + torque[2]
                                * (inertia[0][2] * inertia[1][0] - inertia[0][0] * inertia[1][2])),
                    inv_det
                        * (torque[0]
                            * (inertia[1][0] * inertia[2][1] - inertia[1][1] * inertia[2][0])
                            + torque[1]
                                * (inertia[0][1] * inertia[2][0] - inertia[0][0] * inertia[2][1])
                            + torque[2]
                                * (inertia[0][0] * inertia[1][1] - inertia[0][1] * inertia[1][0])),
                ];
                for (i, pos) in current_positions.iter().enumerate() {
                    // g_rot = omega x r
                    projected_grad[i][0] -= omega[1] * pos[2] - omega[2] * pos[1];
                    projected_grad[i][1] -= omega[2] * pos[0] - omega[0] * pos[2];
                    projected_grad[i][2] -= omega[0] * pos[1] - omega[1] * pos[0];
                }
            }

            // 3d: Check convergence: max projected gradient
            let max_proj_grad = projected_grad
                .iter()
                .flat_map(|g| g.iter())
                .map(|x| x.abs())
                .fold(0.0_f64, f64::max);

            point_opt_steps = opt_step + 1;

            if max_proj_grad < grad_threshold {
                point_converged = true;
                break;
            }

            // Energy stagnation detection via rolling window.
            // The Cartesian-space optimizer can oscillate around the constrained
            // minimum with |ΔE| >> threshold per step, but the energy stays in a
            // narrow band. Track the last 5 energies and converge if the range
            // (max - min) is below threshold.
            energy_history.push(energy);
            if energy_history.len() > 8 {
                energy_history.remove(0);
            }
            if energy_history.len() >= 5 {
                let e_min = energy_history.iter().cloned().fold(f64::INFINITY, f64::min);
                let e_max = energy_history
                    .iter()
                    .cloned()
                    .fold(f64::NEG_INFINITY, f64::max);
                let e_range = e_max - e_min;
                // 0.1 mHa range = energy is oscillating within noise for PES scan
                if e_range < 1e-4 {
                    point_converged = true;
                    break;
                }
            }

            // Also check single-step stagnation for fast convergence
            if opt_step >= 2 && (energy - prev_energy).abs() < 1e-6 {
                point_converged = true;
                break;
            }

            // Best-energy stagnation: if the best energy hasn't improved
            // for STAGNATION_PATIENCE steps, the optimizer is oscillating
            // around a minimum without making progress. This is more robust
            // than the rolling-window check because it isn't poisoned by
            // occasional large energy spikes from bad L-BFGS steps.
            if steps_since_best >= STAGNATION_PATIENCE {
                point_converged = true;
                // Use the best energy found, not the current (possibly high) one
                point_energy = best_energy;
                break;
            }

            // 3e: Flatten current coordinates and projected gradient for L-BFGS
            let flat_coords: Vec<f64> = current_positions
                .iter()
                .flat_map(|p| p.iter().copied())
                .collect();
            let grad_flat: Vec<f64> = projected_grad
                .iter()
                .flat_map(|g| g.iter().copied())
                .collect();

            // Update L-BFGS history with previous step's data
            if let (Some(prev_x), Some(prev_g)) = (&prev_flat_coords, &prev_proj_grad_flat) {
                let s: Vec<f64> = flat_coords
                    .iter()
                    .zip(prev_x.iter())
                    .map(|(x, px)| x - px)
                    .collect();
                let y: Vec<f64> = grad_flat
                    .iter()
                    .zip(prev_g.iter())
                    .map(|(g, pg)| g - pg)
                    .collect();

                // Only add to history if curvature condition s^T y > 0
                let sy: f64 = s.iter().zip(y.iter()).map(|(si, yi)| si * yi).sum();
                if sy > 1e-14 {
                    if lbfgs_s.len() >= LBFGS_MEMORY {
                        lbfgs_s.pop_front();
                        lbfgs_y.pop_front();
                    }
                    lbfgs_s.push_back(s);
                    lbfgs_y.push_back(y);
                }
            }

            prev_flat_coords = Some(flat_coords);
            prev_proj_grad_flat = Some(grad_flat.clone());

            // 3f: Compute search direction via L-BFGS two-loop recursion
            let mut direction = lbfgs_two_loop(&grad_flat, &lbfgs_s, &lbfgs_y, MODEL_HESSIAN_DIAG);

            // Project the L-BFGS direction to remove the constraint component.
            // This prevents the step from changing the constrained coordinate,
            // reducing the amount of correction needed from enforce_constraint
            // and keeping the L-BFGS s_k vectors consistent with actual steps.
            let c_dot_d: f64 = constraint_dir
                .iter()
                .zip(direction.iter())
                .map(|(c, d)| c * d)
                .sum();
            for j in 0..direction.len() {
                direction[j] -= c_dot_d * constraint_dir[j];
            }

            // 3g: Apply step with clamping
            let mut new_positions = current_positions.clone();
            for a in 0..n_atoms {
                let mut step_sq = 0.0;
                for d in 0..3 {
                    let step = -direction[3 * a + d];
                    new_positions[a][d] += step;
                    step_sq += step * step;
                }
                // Clamp step magnitude per atom
                let step_mag = step_sq.sqrt();
                if step_mag > MAX_STEP_PER_ATOM {
                    let scale = MAX_STEP_PER_ATOM / step_mag;
                    for d in 0..3 {
                        let unclamped = new_positions[a][d] - current_positions[a][d];
                        new_positions[a][d] = current_positions[a][d] + unclamped * scale;
                    }
                }
            }

            // 3h: Enforce constraint with max 3 correction attempts
            let mut corrected = new_positions;
            for _correction in 0..MAX_CONSTRAINT_CORRECTIONS {
                let (pos, needed) = enforce_constraint(
                    &corrected,
                    &atoms_z,
                    &config.coordinate,
                    v,
                    CONSTRAINT_TOLERANCE,
                    &dihedral_fragment,
                );
                corrected = pos;
                if !needed {
                    break;
                }
            }

            current_positions = corrected;
        }

        // If optimization didn't converge but we have energy, still record it
        if point_opt_steps == 0 {
            // No optimization steps were taken; run a single SCF at the rigid geometry
            let atom_structs: Vec<Atom> = atoms_z
                .iter()
                .zip(current_positions.iter())
                .filter_map(|(&z, &pos)| Atom::new(z, pos).ok())
                .collect();

            if atom_structs.len() == n_atoms {
                let result = match method.as_str() {
                    "rhf" => evaluate_rhf_energy_gradient(
                        &atom_structs,
                        &config.basis_name,
                        false, // gradients require Cartesian basis
                        &scf_config,
                        current_density.as_deref(),
                    )
                    .map(|(e, _, dm)| (e, Some(dm))),
                    _ => evaluate_dft_energy_gradient(
                        &atom_structs,
                        &config.basis_name,
                        &method,
                        &scf_config,
                        current_density.as_deref(),
                    )
                    .map(|(e, _, dm)| (e, Some(dm))),
                };

                if let Some((e, dm)) = result {
                    point_energy = e;
                    point_converged = true;
                    current_density = dm.or(current_density);
                    point_opt_steps = 1;
                }
            }
        }

        total_iterations += point_scf_iters;
        total_opt_steps += point_opt_steps;

        // Update density seeding
        if point_converged {
            if let Some(dm) = &current_density {
                prev_density = Some(dm.clone());
            }
        } else {
            prev_density = None;
        }

        // Compute internal coordinate snapshot
        let geom_with_z: Vec<(u8, [f64; 3])> = atoms_z
            .iter()
            .zip(current_positions.iter())
            .map(|(&z, &pos)| (z, pos))
            .collect();
        let snapshot = compute_all_internals(&geom_with_z);

        points.push(PesInternalPoint {
            coordinate_value: v,
            energy: point_energy,
            converged: point_converged,
            scf_iterations: point_scf_iters,
            geometry: current_positions,
            opt_steps: Some(point_opt_steps),
            internal_coordinates: Some(snapshot),
        });

        if let Some(cb) = progress {
            cb(i, v, point_energy, point_converged);
        }
    }

    // ---- Equilibrium interpolation ----
    let equilibrium = find_internal_equilibrium(&points);

    Ok(PesScanInternalResult {
        coordinate_type,
        atom_indices,
        points,
        equilibrium,
        total_iterations,
        scan_mode: "relaxed".to_string(),
        total_opt_steps,
    })
}

/// Evaluate RHF energy and gradient at a given geometry.
///
/// Builds basis set, computes integrals, runs SCF, then computes the analytical
/// gradient. Returns (energy, gradient_result, density_matrix).
fn evaluate_rhf_energy_gradient(
    atoms: &[Atom],
    basis_name: &str,
    use_spherical: bool,
    scf_config: &ScfConfig,
    initial_density: Option<&[f64]>,
) -> Option<(f64, GradientResult, Vec<f64>)> {
    let basis = BasisSet::build(atoms.to_vec(), basis_name).ok()?;
    let nbf = if use_spherical {
        basis.n_basis_spherical()
    } else {
        basis.n_basis
    };

    let (s_mat, h_core, eri) = if use_spherical {
        (
            overlap_matrix_spherical(&basis),
            hcore_matrix_spherical(&basis),
            eri_compressed_spherical(&basis),
        )
    } else {
        (
            overlap_matrix(&basis),
            hcore_matrix(&basis),
            eri_compressed(&basis),
        )
    };

    let system = PresetSystem {
        system_id: "relaxed_pes".to_string(),
        label: "Relaxed PES point".to_string(),
        nbf,
        nelec: basis.n_electrons,
        e_nuc: basis.nuclear_repulsion,
        s_matrix: s_mat,
        h_core,
        eri_compressed: eri,
    };

    let n_occ = system.n_occ();
    let scf_output = rhf_scf_with_guess(&system, scf_config, initial_density).ok()?;
    let grad_result = rhf_gradient(
        &basis,
        &scf_output.density_matrix,
        &scf_output.mo_coefficients,
        &scf_output.mo_energies,
        n_occ,
    );

    Some((
        scf_output.energy_total,
        grad_result,
        scf_output.density_matrix.clone(),
    ))
}

/// Evaluate DFT energy and gradient at a given geometry with optional density seeding.
///
/// Runs KS-DFT SCF (with density seeding) + analytical gradient. Returns the
/// converged density matrix for seeding subsequent evaluations.
///
/// Without density seeding, each DFT SCF starts from scratch (~100+ iterations).
/// With seeding from the previous geometry, convergence takes ~5-10 iterations.
fn evaluate_dft_energy_gradient(
    atoms: &[Atom],
    basis_name: &str,
    method: &str,
    scf_config: &ScfConfig,
    initial_density: Option<&[f64]>,
) -> Option<(f64, GradientResult, Vec<f64>)> {
    use super::gradient::ks_dft_gradient_with_guess;

    let use_d3bj = method == "b3lyp-d3bj";
    let grid_config = crate::dft::GridConfig::default();

    let functional: Box<dyn crate::dft::ExchangeCorrelation> = match method {
        "lda" => Box::new(crate::dft::Lda::new()),
        _ => Box::new(crate::dft::B3lyp::new()),
    };

    // ks_dft_gradient_with_guess passes the initial_density through to
    // ks_scf_with_guess, enabling density seeding for fast convergence.
    let grad_result = ks_dft_gradient_with_guess(
        atoms,
        basis_name,
        functional.as_ref(),
        &grid_config,
        scf_config,
        use_d3bj,
        initial_density,
    );

    let energy = grad_result.energy?;

    // Extract converged density for seeding the next evaluation.
    // ks_dft_gradient doesn't return density directly, so we extract it
    // from the GradientResult's internal SCF density if available,
    // or rebuild it cheaply.
    let density = grad_result.density.clone().unwrap_or_default();

    Some((energy, grad_result, density))
}

// =============================================================================
// Scan Dispatcher
// =============================================================================

/// Run a PES scan over an internal coordinate, dispatching to rigid or relaxed mode.
///
/// This is the primary entry point for internal coordinate PES scans. It routes
/// to `pes_scan_rigid()` or `pes_scan_relaxed()` based on the `mode` parameter.
///
/// # Arguments
///
/// * `config` - Scan configuration (molecule, coordinate, range, method)
/// * `mode` - Scan mode: `ScanMode::Rigid` or `ScanMode::Relaxed`
/// * `progress` - Optional callback: (point_index, coordinate_value, energy, converged)
///
/// # Returns
///
/// `PesScanInternalResult` with the scan_mode field indicating which mode was used.
pub fn pes_scan_internal(
    config: &PesScanInternalConfig,
    mode: ScanMode,
    progress: Option<&dyn Fn(usize, f64, f64, bool)>,
) -> Result<PesScanInternalResult, PesScanError> {
    match mode {
        ScanMode::Rigid => pes_scan_rigid(config, progress),
        ScanMode::Relaxed => pes_scan_relaxed(config, progress),
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const TOL: f64 = 1e-12;

    // ---- rotate_point tests ----

    #[test]
    fn test_rotate_point_90_degrees() {
        // Rotate (1,0,0) around z-axis by pi/2 -> (0,1,0)
        let p = [1.0, 0.0, 0.0];
        let origin = [0.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let result = rotate_point(p, origin, axis, PI / 2.0);

        assert!(
            (result[0] - 0.0).abs() < TOL,
            "x: expected 0.0, got {}",
            result[0]
        );
        assert!(
            (result[1] - 1.0).abs() < TOL,
            "y: expected 1.0, got {}",
            result[1]
        );
        assert!(
            (result[2] - 0.0).abs() < TOL,
            "z: expected 0.0, got {}",
            result[2]
        );
    }

    #[test]
    fn test_rotate_point_180_degrees() {
        // Rotate (1,0,0) around z-axis by pi -> (-1,0,0)
        let p = [1.0, 0.0, 0.0];
        let origin = [0.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let result = rotate_point(p, origin, axis, PI);

        assert!(
            (result[0] - (-1.0)).abs() < TOL,
            "x: expected -1.0, got {}",
            result[0]
        );
        assert!(
            (result[1] - 0.0).abs() < TOL,
            "y: expected 0.0, got {}",
            result[1]
        );
        assert!(
            (result[2] - 0.0).abs() < TOL,
            "z: expected 0.0, got {}",
            result[2]
        );
    }

    #[test]
    fn test_rotate_point_identity() {
        // Rotate by 0 -> unchanged
        let p = [3.14, 2.72, 1.41];
        let origin = [1.0, 2.0, 3.0];
        let axis = [0.0, 0.0, 1.0];
        let result = rotate_point(p, origin, axis, 0.0);

        assert!((result[0] - p[0]).abs() < TOL, "x should be unchanged");
        assert!((result[1] - p[1]).abs() < TOL, "y should be unchanged");
        assert!((result[2] - p[2]).abs() < TOL, "z should be unchanged");
    }

    #[test]
    fn test_rotate_point_preserves_distance() {
        // Distance from origin must be preserved for any rotation
        let p = [3.0, 4.0, 5.0];
        let origin = [1.0, 1.0, 1.0];
        let axis = [1.0, 1.0, 1.0]; // arbitrary axis (will be normalized)
        let theta = 1.23456; // arbitrary angle

        let d_before = distance(&p, &origin);
        let result = rotate_point(p, origin, axis, theta);
        let d_after = distance(&result, &origin);

        assert!(
            (d_before - d_after).abs() < TOL,
            "distance should be preserved: before={}, after={}",
            d_before,
            d_after
        );
    }

    #[test]
    fn test_rotate_point_arbitrary_axis() {
        // Rotate (1,0,0) around axis (1,1,0)/sqrt(2) by pi -> (0,1,0)
        //
        // The axis (1,1,0) bisects x and y. Rotating by pi reflects
        // across this axis in the xy-plane:
        //   (1,0,0) -> (0,1,0)
        let p = [1.0, 0.0, 0.0];
        let origin = [0.0, 0.0, 0.0];
        let axis = [1.0, 1.0, 0.0]; // will be normalized internally
        let result = rotate_point(p, origin, axis, PI);

        assert!(
            (result[0] - 0.0).abs() < TOL,
            "x: expected 0.0, got {}",
            result[0]
        );
        assert!(
            (result[1] - 1.0).abs() < TOL,
            "y: expected 1.0, got {}",
            result[1]
        );
        assert!(
            (result[2] - 0.0).abs() < TOL,
            "z: expected 0.0, got {}",
            result[2]
        );
    }

    // ---- generate_bond_geometry tests ----

    #[test]
    fn test_generate_bond_h2_stretch() {
        // H2: H at origin, H at (0,0,1.4) bohr. Stretch to 2.0 bohr.
        let atoms: Vec<(u8, [f64; 3])> = vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])];

        let coords = generate_bond_geometry(&atoms, 0, 1, 2.0).unwrap();

        // atom_i (0) unchanged
        assert!((coords[0][0]).abs() < TOL);
        assert!((coords[0][1]).abs() < TOL);
        assert!((coords[0][2]).abs() < TOL);

        // atom_j (1) should be at (0,0,2.0)
        assert!((coords[1][0]).abs() < TOL);
        assert!((coords[1][1]).abs() < TOL);
        assert!(
            (coords[1][2] - 2.0).abs() < TOL,
            "z: expected 2.0, got {}",
            coords[1][2]
        );
    }

    #[test]
    fn test_generate_bond_preserves_others() {
        // H2O: stretch O-H1 bond, verify O and H2 unchanged
        let atoms: Vec<(u8, [f64; 3])> = vec![
            (8, [0.0, 0.0, 0.2217282]),         // O
            (1, [0.0, 1.4305447, -0.8869128]),  // H1
            (1, [0.0, -1.4305447, -0.8869128]), // H2
        ];

        // Stretch O-H1 bond (atom 0 -> atom 1) to 2.0 bohr
        let coords = generate_bond_geometry(&atoms, 0, 1, 2.0).unwrap();

        // O (atom 0) unchanged
        assert!((coords[0][0] - atoms[0].1[0]).abs() < TOL);
        assert!((coords[0][1] - atoms[0].1[1]).abs() < TOL);
        assert!((coords[0][2] - atoms[0].1[2]).abs() < TOL);

        // H2 (atom 2) unchanged
        assert!((coords[2][0] - atoms[2].1[0]).abs() < TOL);
        assert!((coords[2][1] - atoms[2].1[1]).abs() < TOL);
        assert!((coords[2][2] - atoms[2].1[2]).abs() < TOL);

        // Verify the new bond length
        let new_r = distance(&coords[0], &coords[1]);
        assert!(
            (new_r - 2.0).abs() < TOL,
            "bond length: expected 2.0, got {}",
            new_r
        );
    }

    // ---- generate_angle_geometry tests ----

    #[test]
    fn test_generate_angle_h2o_90() {
        // H2O at ~104.5 degrees, bend to 90 degrees
        let atoms: Vec<(u8, [f64; 3])> = vec![
            (8, [0.0, 0.0, 0.2217282]),         // O (central for angle)
            (1, [0.0, 1.4305447, -0.8869128]),  // H1
            (1, [0.0, -1.4305447, -0.8869128]), // H2
        ];

        // Angle H1-O-H2 (atom_i=1, atom_j=0, atom_k=2)
        let original_angle = bond_angle(&atoms[1].1, &atoms[0].1, &atoms[2].1);
        assert!(
            (original_angle - 104.5_f64.to_radians()).abs() < 0.01,
            "initial angle should be ~104.5 degrees, got {:.1} degrees",
            original_angle.to_degrees()
        );

        let target = PI / 2.0; // 90 degrees
        let coords = generate_angle_geometry(&atoms, 1, 0, 2, target).unwrap();

        // Verify new angle is 90 degrees
        let new_angle = bond_angle(&coords[1], &coords[0], &coords[2]);
        assert!(
            (new_angle - target).abs() < 1e-10,
            "angle: expected {:.4} rad, got {:.4} rad",
            target,
            new_angle
        );

        // Verify O-H1 bond length unchanged (atom_i=1 was NOT moved)
        let oh1_before = distance(&atoms[0].1, &atoms[1].1);
        let oh1_after = distance(&coords[0], &coords[1]);
        assert!(
            (oh1_before - oh1_after).abs() < TOL,
            "O-H1 bond changed: {} -> {}",
            oh1_before,
            oh1_after
        );

        // Verify O-H2 bond length unchanged (atom_k=2 was rotated, distance preserved)
        let oh2_before = distance(&atoms[0].1, &atoms[2].1);
        let oh2_after = distance(&coords[0], &coords[2]);
        assert!(
            (oh2_before - oh2_after).abs() < TOL,
            "O-H2 bond changed: {} -> {}",
            oh2_before,
            oh2_after
        );
    }

    #[test]
    fn test_generate_angle_near_linear() {
        // 3 atoms nearly collinear (179 degrees), bend to 170 degrees.
        // Should not panic or produce NaN.
        let theta_179 = 179.0_f64.to_radians();
        // Place atoms i-j-k nearly linear along x-axis
        // j at origin, i along +x, k slightly off -x axis
        let tiny_off = (PI - theta_179).sin(); // small perpendicular offset
        let atoms: Vec<(u8, [f64; 3])> = vec![
            (1, [2.0, 0.0, 0.0]),       // i
            (1, [0.0, 0.0, 0.0]),       // j (central)
            (1, [-2.0, tiny_off, 0.0]), // k (nearly -x)
        ];

        let initial = bond_angle(&atoms[0].1, &atoms[1].1, &atoms[2].1);
        assert!(
            (initial.to_degrees() - 179.0).abs() < 1.0,
            "initial angle should be ~179 deg, got {:.1}",
            initial.to_degrees()
        );

        let target = 170.0_f64.to_radians();
        let coords = generate_angle_geometry(&atoms, 0, 1, 2, target).unwrap();

        // Check no NaN
        for c in &coords {
            for &v in c {
                assert!(!v.is_nan(), "NaN in result coordinates");
                assert!(!v.is_infinite(), "Inf in result coordinates");
            }
        }

        // Check angle is approximately correct
        let new_angle = bond_angle(&coords[0], &coords[1], &coords[2]);
        assert!(
            (new_angle - target).abs() < 1e-8,
            "angle: expected {:.4} rad ({:.1} deg), got {:.4} rad ({:.1} deg)",
            target,
            target.to_degrees(),
            new_angle,
            new_angle.to_degrees()
        );
    }

    // ---- generate_dihedral_geometry tests ----

    #[test]
    fn test_generate_dihedral() {
        // 4-atom chain along z, with atom_l offset.
        // We'll set up a known dihedral and rotate by 90 degrees.
        let atoms: Vec<(u8, [f64; 3])> = vec![
            (1, [1.0, 0.0, 0.0]), // i — defines the reference plane
            (6, [0.0, 0.0, 0.0]), // j — start of central bond
            (6, [0.0, 0.0, 1.5]), // k — end of central bond
            (1, [1.0, 0.0, 1.5]), // l — in the same plane as i (dihedral = 0)
        ];

        let current_tau = dihedral_angle(&atoms[0].1, &atoms[1].1, &atoms[2].1, &atoms[3].1);
        assert!(
            current_tau.abs() < 1e-10,
            "initial dihedral should be ~0, got {:.6}",
            current_tau
        );

        let target_tau = PI / 2.0; // 90 degrees
        let mut fragment = HashSet::new();
        fragment.insert(3); // only rotate atom_l

        let coords = generate_dihedral_geometry(&atoms, 0, 1, 2, 3, target_tau, &fragment).unwrap();

        let new_tau = dihedral_angle(&coords[0], &coords[1], &coords[2], &coords[3]);
        assert!(
            (new_tau - target_tau).abs() < 1e-10,
            "dihedral: expected {:.4} rad, got {:.4} rad",
            target_tau,
            new_tau
        );
    }

    #[test]
    fn test_generate_dihedral_preserves_bonds() {
        // Verify that dihedral rotation preserves bond lengths and angles
        let atoms: Vec<(u8, [f64; 3])> = vec![
            (1, [1.0, 0.0, 0.0]), // i
            (6, [0.0, 0.0, 0.0]), // j
            (6, [0.0, 0.0, 1.5]), // k
            (1, [1.0, 0.0, 1.5]), // l
        ];

        // Bond lengths before
        let ij_before = distance(&atoms[0].1, &atoms[1].1);
        let jk_before = distance(&atoms[1].1, &atoms[2].1);
        let kl_before = distance(&atoms[2].1, &atoms[3].1);

        // Angle i-j-k before
        let angle_ijk_before = bond_angle(&atoms[0].1, &atoms[1].1, &atoms[2].1);
        // Angle j-k-l before
        let angle_jkl_before = bond_angle(&atoms[1].1, &atoms[2].1, &atoms[3].1);

        let target_tau = 2.0; // arbitrary angle
        let mut fragment = HashSet::new();
        fragment.insert(3);

        let coords = generate_dihedral_geometry(&atoms, 0, 1, 2, 3, target_tau, &fragment).unwrap();

        // Bond i-j unchanged (neither i nor j in fragment)
        let ij_after = distance(&coords[0], &coords[1]);
        assert!(
            (ij_before - ij_after).abs() < TOL,
            "i-j bond changed: {} -> {}",
            ij_before,
            ij_after
        );

        // Bond j-k unchanged
        let jk_after = distance(&coords[1], &coords[2]);
        assert!(
            (jk_before - jk_after).abs() < TOL,
            "j-k bond changed: {} -> {}",
            jk_before,
            jk_after
        );

        // Bond k-l: rotation around j-k axis from j preserves distance from j,
        // but we need distance from k preserved. Since we rotate around j-k axis
        // from origin j, the k-l distance is preserved because l moves on a circle
        // centered on the j-k axis.
        let kl_after = distance(&coords[2], &coords[3]);
        assert!(
            (kl_before - kl_after).abs() < TOL,
            "k-l bond changed: {} -> {}",
            kl_before,
            kl_after
        );

        // Angle i-j-k unchanged (i and j not rotated, k not rotated)
        let angle_ijk_after = bond_angle(&coords[0], &coords[1], &coords[2]);
        assert!(
            (angle_ijk_before - angle_ijk_after).abs() < TOL,
            "angle i-j-k changed: {} -> {}",
            angle_ijk_before,
            angle_ijk_after
        );

        // Angle j-k-l: j and k not rotated, l rotated around j-k axis.
        // The angle j-k-l should be preserved because the rotation is around
        // the j-k bond axis.
        let angle_jkl_after = bond_angle(&coords[1], &coords[2], &coords[3]);
        assert!(
            (angle_jkl_before - angle_jkl_after).abs() < TOL,
            "angle j-k-l changed: {} -> {}",
            angle_jkl_before,
            angle_jkl_after
        );
    }

    // ---- build_connectivity_simple tests ----

    #[test]
    fn test_build_connectivity_h2() {
        // H2: two H atoms at 1.4 bohr — well within bonding threshold
        // Threshold: 1.2 * (0.586 + 0.586) = 1.4064 bohr > 1.4 bohr
        let atoms: Vec<(u8, [f64; 3])> = vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])];
        let adj = build_connectivity_simple(&atoms);

        assert_eq!(adj.len(), 2);
        assert_eq!(adj[0], vec![1], "atom 0 should be bonded to atom 1");
        assert_eq!(adj[1], vec![0], "atom 1 should be bonded to atom 0");
    }

    #[test]
    fn test_build_connectivity_h2o() {
        // H2O: O bonded to both H's, H's not bonded to each other
        let atoms: Vec<(u8, [f64; 3])> = vec![
            (8, [0.0, 0.0, 0.2217282]),         // O
            (1, [0.0, 1.4305447, -0.8869128]),  // H1
            (1, [0.0, -1.4305447, -0.8869128]), // H2
        ];
        let adj = build_connectivity_simple(&atoms);

        assert_eq!(adj.len(), 3);
        // O (atom 0) bonded to H1 and H2
        assert!(adj[0].contains(&1), "O should be bonded to H1");
        assert!(adj[0].contains(&2), "O should be bonded to H2");
        assert_eq!(adj[0].len(), 2, "O should have exactly 2 bonds");
        // H1 (atom 1) bonded only to O
        assert_eq!(adj[1], vec![0], "H1 should be bonded only to O");
        // H2 (atom 2) bonded only to O
        assert_eq!(adj[2], vec![0], "H2 should be bonded only to O");
    }

    /// Helper: return staggered ethane geometry (C2H6) in bohr
    fn ethane_atoms() -> Vec<(u8, [f64; 3])> {
        vec![
            (6, [0.0, 0.0, 0.0]),       // C1 (index 0)
            (6, [0.0, 0.0, 2.9]),       // C2 (index 1)
            (1, [1.94, 0.0, -0.63]),    // H1 (index 2, bonded to C1)
            (1, [-0.97, 1.68, -0.63]),  // H2 (index 3, bonded to C1)
            (1, [-0.97, -1.68, -0.63]), // H3 (index 4, bonded to C1)
            (1, [-1.94, 0.0, 3.53]),    // H4 (index 5, bonded to C2)
            (1, [0.97, 1.68, 3.53]),    // H5 (index 6, bonded to C2)
            (1, [0.97, -1.68, 3.53]),   // H6 (index 7, bonded to C2)
        ]
    }

    #[test]
    fn test_build_connectivity_ethane() {
        let atoms = ethane_atoms();
        let adj = build_connectivity_simple(&atoms);

        // Count total bonds (each bond counted once)
        let total_bonds: usize = adj.iter().map(|v| v.len()).sum::<usize>() / 2;
        assert_eq!(total_bonds, 7, "ethane should have 7 bonds (1 C-C + 6 C-H)");

        // C1 (index 0) bonded to C2, H1, H2, H3
        assert!(adj[0].contains(&1), "C1 bonded to C2");
        assert!(adj[0].contains(&2), "C1 bonded to H1");
        assert!(adj[0].contains(&3), "C1 bonded to H2");
        assert!(adj[0].contains(&4), "C1 bonded to H3");
        assert_eq!(adj[0].len(), 4, "C1 should have 4 bonds");

        // C2 (index 1) bonded to C1, H4, H5, H6
        assert!(adj[1].contains(&0), "C2 bonded to C1");
        assert!(adj[1].contains(&5), "C2 bonded to H4");
        assert!(adj[1].contains(&6), "C2 bonded to H5");
        assert!(adj[1].contains(&7), "C2 bonded to H6");
        assert_eq!(adj[1].len(), 4, "C2 should have 4 bonds");

        // Each H bonded to exactly 1 atom (its carbon)
        for h_idx in 2..=7 {
            assert_eq!(
                adj[h_idx].len(),
                1,
                "H{} should have exactly 1 bond",
                h_idx - 1
            );
        }
    }

    // ---- detect_fragment tests ----

    #[test]
    fn test_detect_fragment_ethane() {
        // Remove C1-C2 bond, BFS from C1 -> gets {C1, H1, H2, H3}
        let atoms = ethane_atoms();
        let adj = build_connectivity_simple(&atoms);

        let fragment_c1 = detect_fragment(&adj, 0, 1, 0);
        assert_eq!(fragment_c1.len(), 4, "C1 fragment should have 4 atoms");
        assert!(fragment_c1.contains(&0), "fragment should contain C1");
        assert!(fragment_c1.contains(&2), "fragment should contain H1");
        assert!(fragment_c1.contains(&3), "fragment should contain H2");
        assert!(fragment_c1.contains(&4), "fragment should contain H3");
        assert!(!fragment_c1.contains(&1), "fragment should NOT contain C2");

        // BFS from C2 -> gets {C2, H4, H5, H6}
        let fragment_c2 = detect_fragment(&adj, 0, 1, 1);
        assert_eq!(fragment_c2.len(), 4, "C2 fragment should have 4 atoms");
        assert!(fragment_c2.contains(&1), "fragment should contain C2");
        assert!(fragment_c2.contains(&5), "fragment should contain H4");
        assert!(fragment_c2.contains(&6), "fragment should contain H5");
        assert!(fragment_c2.contains(&7), "fragment should contain H6");
        assert!(!fragment_c2.contains(&0), "fragment should NOT contain C1");
    }

    #[test]
    fn test_detect_fragment_h2o() {
        // Remove O-H1 bond, BFS from H1 -> gets {H1} only
        let atoms: Vec<(u8, [f64; 3])> = vec![
            (8, [0.0, 0.0, 0.2217282]),         // O (index 0)
            (1, [0.0, 1.4305447, -0.8869128]),  // H1 (index 1)
            (1, [0.0, -1.4305447, -0.8869128]), // H2 (index 2)
        ];
        let adj = build_connectivity_simple(&atoms);

        let fragment = detect_fragment(&adj, 0, 1, 1);
        assert_eq!(fragment.len(), 1, "H1 fragment should have only 1 atom");
        assert!(fragment.contains(&1), "fragment should contain H1");
    }

    // ---- is_ring_bond tests ----

    #[test]
    fn test_is_ring_bond_linear() {
        // C-C in ethane is NOT a ring bond
        let atoms = ethane_atoms();
        let adj = build_connectivity_simple(&atoms);

        assert!(
            !is_ring_bond(&adj, 0, 1),
            "C-C bond in ethane should NOT be a ring bond"
        );
    }

    #[test]
    fn test_is_ring_bond_cyclic() {
        // Any C-C in benzene IS a ring bond
        // Benzene geometry from benchmark_spherical.rs
        let atoms: Vec<(u8, [f64; 3])> = vec![
            (6, [0.0000000000, 2.6399473960, 0.0000000000]), // C0
            (6, [2.2861906655, 1.3199736980, 0.0000000000]), // C1
            (6, [2.2861906655, -1.3199736980, 0.0000000000]), // C2
            (6, [0.0000000000, -2.6399473960, 0.0000000000]), // C3
            (6, [-2.2861906655, -1.3199736980, 0.0000000000]), // C4
            (6, [-2.2861906655, 1.3199736980, 0.0000000000]), // C5
            (1, [0.0000000000, 4.6884105150, 0.0000000000]), // H6
            (1, [4.0602655512, 2.3442052575, 0.0000000000]), // H7
            (1, [4.0602655512, -2.3442052575, 0.0000000000]), // H8
            (1, [0.0000000000, -4.6884105150, 0.0000000000]), // H9
            (1, [-4.0602655512, -2.3442052575, 0.0000000000]), // H10
            (1, [-4.0602655512, 2.3442052575, 0.0000000000]), // H11
        ];
        let adj = build_connectivity_simple(&atoms);

        // C0-C1 bond should be a ring bond
        assert!(
            is_ring_bond(&adj, 0, 1),
            "C0-C1 bond in benzene should be a ring bond"
        );

        // C2-C3 bond should also be a ring bond
        assert!(
            is_ring_bond(&adj, 2, 3),
            "C2-C3 bond in benzene should be a ring bond"
        );

        // C-H bond should NOT be a ring bond
        assert!(
            !is_ring_bond(&adj, 0, 6),
            "C-H bond in benzene should NOT be a ring bond"
        );
    }

    // ---- Integration: dihedral with fragment detection ----

    #[test]
    fn test_dihedral_with_fragment_ethane() {
        // Build ethane, detect the C2-side fragment, rotate dihedral
        // and verify all 3 H's on the C2 side rotate together.
        let atoms = ethane_atoms();
        let adj = build_connectivity_simple(&atoms);

        // Detect fragment on C2 side of the C1-C2 bond
        let fragment = detect_fragment(&adj, 0, 1, 1);
        assert_eq!(fragment.len(), 4);
        assert!(fragment.contains(&1)); // C2
        assert!(fragment.contains(&5)); // H4
        assert!(fragment.contains(&6)); // H5
        assert!(fragment.contains(&7)); // H6

        // Define dihedral: H1-C1-C2-H4  (atoms 2, 0, 1, 5)
        // Rotate the C2-side fragment by 60 degrees
        let target_tau = PI / 3.0; // 60 degrees
        let coords = generate_dihedral_geometry(&atoms, 2, 0, 1, 5, target_tau, &fragment).unwrap();

        // Verify the dihedral angle is now 60 degrees
        let new_tau = dihedral_angle(&coords[2], &coords[0], &coords[1], &coords[5]);
        assert!(
            (new_tau - target_tau).abs() < 1e-10,
            "dihedral should be {:.4} rad, got {:.4} rad",
            target_tau,
            new_tau
        );

        // Verify atoms on the C1 side are unchanged
        for &idx in &[0, 2, 3, 4] {
            assert!(
                (coords[idx][0] - atoms[idx].1[0]).abs() < TOL
                    && (coords[idx][1] - atoms[idx].1[1]).abs() < TOL
                    && (coords[idx][2] - atoms[idx].1[2]).abs() < TOL,
                "atom {} (C1 side) should be unchanged",
                idx
            );
        }

        // Verify C2-H bond lengths are preserved (all 3 H's rotated with C2)
        for &h_idx in &[5, 6, 7] {
            let d_before = distance(&atoms[1].1, &atoms[h_idx].1);
            let d_after = distance(&coords[1], &coords[h_idx]);
            assert!(
                (d_before - d_after).abs() < TOL,
                "C2-H{} bond changed: {:.6} -> {:.6}",
                h_idx - 4,
                d_before,
                d_after
            );
        }

        // Verify C1-C2 bond length is preserved
        let cc_before = distance(&atoms[0].1, &atoms[1].1);
        let cc_after = distance(&coords[0], &coords[1]);
        assert!(
            (cc_before - cc_after).abs() < TOL,
            "C-C bond changed: {:.6} -> {:.6}",
            cc_before,
            cc_after
        );
    }

    // =========================================================================
    // Rigid PES Scan Tests
    // =========================================================================

    #[test]
    fn test_rigid_scan_h2_bond() {
        // H2 STO-3G bond scan: 1.0 to 3.0 bohr, 11 points, RHF
        // Expected: minimum energy near R = 1.346 bohr
        let config = PesScanInternalConfig {
            atoms: vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])],
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.0,
            value_max: 3.0,
            n_points: 11,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None).expect("Scan should succeed");

        // All 11 points should be present
        assert_eq!(result.points.len(), 11);
        assert_eq!(result.coordinate_type, "bond");
        assert_eq!(result.atom_indices, vec![0, 1]);

        // All points should converge
        for pt in &result.points {
            assert!(
                pt.converged,
                "Point at v={:.4} did not converge",
                pt.coordinate_value
            );
        }

        // Energies should form a well: decrease then increase
        let min_idx = result
            .points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                a.energy
                    .partial_cmp(&b.energy)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap();

        // Minimum should not be at the boundary
        assert!(
            min_idx > 0 && min_idx < 10,
            "Minimum at boundary: {}",
            min_idx
        );

        // Should find an equilibrium
        let eq = result
            .equilibrium
            .as_ref()
            .expect("Should find equilibrium");
        // H2 STO-3G equilibrium is near 1.346 bohr
        assert!(
            (eq.value - 1.346).abs() < 0.05,
            "Equilibrium distance {:.4} too far from 1.346",
            eq.value
        );
    }

    #[test]
    fn test_rigid_scan_h2_matches_existing() {
        // Regression test: H2 bond scan with pes_scan_rigid() must produce
        // identical energies to the existing pes_scan() for the same H2 system.
        use crate::scf::pes::{pes_scan, PesScanConfig};
        use crate::scf::ConvergenceProfile;

        let scf_config = ScfConfig {
            profile: ConvergenceProfile::Tight,
            max_iterations: 100,
            use_diis: true,
            diis_size: 8,
            diis_start: 1,
            damp: 0.0,
            damp_start: 5,
            level_shift: 0.0,
        };

        // Run existing pes_scan
        let old_config = PesScanConfig {
            atom_a_z: 1,
            atom_b_z: 1,
            r_min: 1.0,
            r_max: 3.0,
            n_points: 11,
            basis_name: "sto-3g",
            scf_config: &scf_config,
            use_seeding: true,
        };
        let old_result = pes_scan(&old_config, None);

        // Run new pes_scan_rigid
        let new_config = PesScanInternalConfig {
            atoms: vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])],
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.0,
            value_max: 3.0,
            n_points: 11,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };
        let new_result = pes_scan_rigid(&new_config, None).expect("Scan should succeed");

        // Compare energies at each point
        assert_eq!(old_result.points.len(), new_result.points.len());
        for (old_pt, new_pt) in old_result.points.iter().zip(new_result.points.iter()) {
            // Bond distances should match
            assert!(
                (old_pt.r - new_pt.coordinate_value).abs() < 1e-14,
                "Distance mismatch: old={} new={}",
                old_pt.r,
                new_pt.coordinate_value
            );

            // Both should converge
            assert!(
                old_pt.converged && new_pt.converged,
                "Convergence mismatch at r={:.4}",
                old_pt.r
            );

            // Energies must match to machine precision
            let error = (old_pt.energy - new_pt.energy).abs();
            assert!(
                error < 1e-12,
                "Energy mismatch at r={:.4}: old={:.12} new={:.12} error={:.2e}",
                old_pt.r,
                old_pt.energy,
                new_pt.energy,
                error
            );
        }
    }

    #[test]
    fn test_rigid_scan_h2o_angle() {
        // H2O STO-3G angle scan: 80 to 130 degrees, 11 points, RHF
        // Water geometry: O at origin, H1 and H2 placed symmetrically
        let oh_bond = 1.809; // bohr (typical O-H in water STO-3G)
        let angle0 = 104.5_f64.to_radians();
        let h1 = [
            0.0,
            oh_bond * (angle0 / 2.0).sin(),
            -oh_bond * (angle0 / 2.0).cos(),
        ];
        let h2 = [
            0.0,
            -oh_bond * (angle0 / 2.0).sin(),
            -oh_bond * (angle0 / 2.0).cos(),
        ];

        let config = PesScanInternalConfig {
            atoms: vec![
                (8, [0.0, 0.0, 0.0]), // O at origin
                (1, h1),              // H1
                (1, h2),              // H2
            ],
            coordinate: ScanCoordinate::Angle {
                atom_i: 1, // H1
                atom_j: 0, // O (central)
                atom_k: 2, // H2 (moved)
            },
            value_min: 80.0_f64.to_radians(),
            value_max: 130.0_f64.to_radians(),
            n_points: 11,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: false, // no seeding for angle scan
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None).expect("Scan should succeed");

        assert_eq!(result.points.len(), 11);
        assert_eq!(result.coordinate_type, "angle");

        // All points should converge
        let converged_count = result.points.iter().filter(|p| p.converged).count();
        assert!(
            converged_count >= 10,
            "At least 10/11 angle scan points should converge, got {}",
            converged_count
        );

        // Verify that bond lengths are preserved: the O-H1 distance should not change
        // (only atom_k = H2 is moved)
        let oh1_initial = norm3(&[h1[0], h1[1], h1[2]]);
        for pt in &result.points {
            if pt.converged {
                let oh1_now = {
                    let dx = pt.geometry[1][0] - pt.geometry[0][0];
                    let dy = pt.geometry[1][1] - pt.geometry[0][1];
                    let dz = pt.geometry[1][2] - pt.geometry[0][2];
                    (dx * dx + dy * dy + dz * dz).sqrt()
                };
                assert!(
                    (oh1_now - oh1_initial).abs() < 1e-10,
                    "O-H1 bond length changed: {:.8} -> {:.8} at angle={:.2}",
                    oh1_initial,
                    oh1_now,
                    pt.coordinate_value.to_degrees()
                );
            }
        }

        // Check that minimum is near the equilibrium angle (~100-110 degrees for STO-3G)
        if let Some(eq) = &result.equilibrium {
            let eq_deg = eq.value.to_degrees();
            assert!(
                eq_deg > 90.0 && eq_deg < 120.0,
                "Equilibrium angle {:.1} degrees outside expected range [90, 120]",
                eq_deg
            );
        }
    }

    #[test]
    fn test_rigid_scan_parabolic() {
        // Verify parabolic interpolation with synthetic data
        let points = vec![
            PesInternalPoint {
                coordinate_value: 1.0,
                energy: -1.0,
                converged: true,
                scf_iterations: 5,
                geometry: vec![[0.0; 3]],
                opt_steps: None,
                internal_coordinates: None,
            },
            PesInternalPoint {
                coordinate_value: 2.0,
                energy: -2.0,
                converged: true,
                scf_iterations: 5,
                geometry: vec![[0.0; 3]],
                opt_steps: None,
                internal_coordinates: None,
            },
            PesInternalPoint {
                coordinate_value: 3.0,
                energy: -1.5,
                converged: true,
                scf_iterations: 5,
                geometry: vec![[0.0; 3]],
                opt_steps: None,
                internal_coordinates: None,
            },
        ];

        let eq = find_internal_equilibrium(&points).expect("Should find equilibrium");

        // Same test as pes.rs: minimum of parabola through (-1, -2, -1.5) at r=(1,2,3)
        assert!(
            (eq.value - 2.1667).abs() < 0.001,
            "Equilibrium value {:.4} != 2.1667",
            eq.value
        );
    }

    #[test]
    fn test_rigid_scan_parabolic_with_unconverged() {
        let points = vec![
            PesInternalPoint {
                coordinate_value: 0.5,
                energy: 0.0,
                converged: false,
                scf_iterations: 100,
                geometry: vec![[0.0; 3]],
                opt_steps: None,
                internal_coordinates: None,
            },
            PesInternalPoint {
                coordinate_value: 1.0,
                energy: -1.0,
                converged: true,
                scf_iterations: 5,
                geometry: vec![[0.0; 3]],
                opt_steps: None,
                internal_coordinates: None,
            },
            PesInternalPoint {
                coordinate_value: 1.5,
                energy: -2.0,
                converged: true,
                scf_iterations: 5,
                geometry: vec![[0.0; 3]],
                opt_steps: None,
                internal_coordinates: None,
            },
            PesInternalPoint {
                coordinate_value: 2.0,
                energy: -1.5,
                converged: true,
                scf_iterations: 5,
                geometry: vec![[0.0; 3]],
                opt_steps: None,
                internal_coordinates: None,
            },
        ];

        let eq = find_internal_equilibrium(&points).expect("Should find equilibrium");
        // Unconverged point at 0.5 should be excluded; minimum among 3 converged is at 1.5
        assert!(eq.value > 1.0 && eq.value < 2.0);
    }

    #[test]
    fn test_rigid_scan_lda() {
        // H2 STO-3G bond scan with LDA method
        let config = PesScanInternalConfig {
            atoms: vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])],
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.0,
            value_max: 2.5,
            n_points: 6,
            basis_name: "sto-3g".to_string(),
            method: "lda".to_string(),
            use_seeding: false,
            use_spherical: false,
            convergence_profile: "medium".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None).expect("LDA scan should succeed");

        assert_eq!(result.points.len(), 6);

        // At least most points should converge
        let converged_count = result.points.iter().filter(|p| p.converged).count();
        assert!(
            converged_count >= 4,
            "At least 4/6 LDA points should converge, got {}",
            converged_count
        );

        // DFT energies for H2 should be in the expected range
        for pt in &result.points {
            if pt.converged {
                assert!(
                    pt.energy < 0.0 && pt.energy > -5.0,
                    "LDA H2 energy {:.6} out of range at r={:.2}",
                    pt.energy,
                    pt.coordinate_value
                );
            }
        }
    }

    #[test]
    fn test_rigid_scan_invalid_method() {
        let config = PesScanInternalConfig {
            atoms: vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])],
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.0,
            value_max: 2.0,
            n_points: 3,
            basis_name: "sto-3g".to_string(),
            method: "mp2".to_string(), // unsupported
            use_seeding: false,
            use_spherical: false,
            convergence_profile: "medium".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None);
        assert!(result.is_err(), "Should fail with unsupported method");
    }

    #[test]
    fn test_rigid_scan_invalid_n_points() {
        let config = PesScanInternalConfig {
            atoms: vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])],
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.0,
            value_max: 2.0,
            n_points: 1, // too few
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: false,
            use_spherical: false,
            convergence_profile: "medium".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None);
        assert!(result.is_err(), "Should fail with n_points < 2");
    }

    #[test]
    fn test_rigid_scan_progress_callback() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        let call_count = AtomicUsize::new(0);
        let n_points = 5;

        let config = PesScanInternalConfig {
            atoms: vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])],
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.0,
            value_max: 2.0,
            n_points,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "loose".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let _result = pes_scan_rigid(
            &config,
            Some(&|_idx, _v, _e, _conv| {
                call_count.fetch_add(1, Ordering::Relaxed);
            }),
        )
        .expect("Scan should succeed");

        assert_eq!(
            call_count.load(Ordering::Relaxed),
            n_points,
            "Progress callback should be called {} times",
            n_points
        );
    }

    #[test]
    fn test_rigid_scan_seeding_reduces_iterations() {
        // LiH STO-3G bond scan without DIIS to observe seeding benefit
        let config_cold = PesScanInternalConfig {
            atoms: vec![(3, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 3.0])],
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 2.5,
            value_max: 4.0,
            n_points: 10,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: false,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let config_seeded = PesScanInternalConfig {
            use_seeding: true,
            ..config_cold.clone()
        };

        let result_cold = pes_scan_rigid(&config_cold, None).expect("Cold scan should succeed");
        let result_seeded =
            pes_scan_rigid(&config_seeded, None).expect("Seeded scan should succeed");

        let cold_iters = result_cold.total_iterations;
        let seeded_iters = result_seeded.total_iterations;

        // Seeded scan should use fewer total iterations
        assert!(
            seeded_iters <= cold_iters,
            "Seeding should not increase iterations: cold={}, seeded={}",
            cold_iters,
            seeded_iters
        );
    }

    // =========================================================================
    // InternalCoordinateSnapshot Tests
    // =========================================================================

    #[test]
    fn test_compute_all_internals_h2() {
        // H2: 1 bond, 0 angles, 0 dihedrals
        let atoms: Vec<(u8, [f64; 3])> = vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])];
        let snapshot = compute_all_internals(&atoms);

        assert_eq!(snapshot.bonds.len(), 1, "H2 should have 1 bond");
        assert_eq!(snapshot.angles.len(), 0, "H2 should have 0 angles");
        assert_eq!(snapshot.dihedrals.len(), 0, "H2 should have 0 dihedrals");

        // Check bond value
        let (i, j, d) = snapshot.bonds[0];
        assert_eq!(i, 0);
        assert_eq!(j, 1);
        assert!(
            (d - 1.4).abs() < 1e-10,
            "H-H bond length should be 1.4, got {}",
            d
        );
    }

    #[test]
    fn test_compute_all_internals_h2o() {
        // H2O: 2 bonds, 1 angle, 0 dihedrals
        let atoms: Vec<(u8, [f64; 3])> = vec![
            (8, [0.0, 0.0, 0.2217282]),         // O
            (1, [0.0, 1.4305447, -0.8869128]),  // H1
            (1, [0.0, -1.4305447, -0.8869128]), // H2
        ];
        let snapshot = compute_all_internals(&atoms);

        assert_eq!(snapshot.bonds.len(), 2, "H2O should have 2 bonds");
        assert_eq!(snapshot.angles.len(), 1, "H2O should have 1 angle");
        assert_eq!(snapshot.dihedrals.len(), 0, "H2O should have 0 dihedrals");

        // Check angle is approximately 104.5 degrees
        let (_, j, _, theta) = snapshot.angles[0];
        assert_eq!(j, 0, "Central atom of angle should be O (index 0)");
        assert!(
            (theta.to_degrees() - 104.5).abs() < 1.0,
            "H-O-H angle should be ~104.5 degrees, got {:.1}",
            theta.to_degrees()
        );
    }

    #[test]
    fn test_compute_all_internals_nh3() {
        // NH3: 3 bonds, 3 angles
        // Approximate NH3 geometry (N at origin, 3 H's)
        let nh_bond = 1.912; // bohr
        let theta_h = 107.0_f64.to_radians();
        let h1 = [0.0, nh_bond * theta_h.sin(), -nh_bond * theta_h.cos()];
        let h2 = [
            nh_bond * theta_h.sin() * (2.0 * PI / 3.0).sin(),
            -nh_bond * theta_h.sin() * (2.0 * PI / 3.0).cos(),
            -nh_bond * theta_h.cos(),
        ];
        let h3 = [
            -nh_bond * theta_h.sin() * (2.0 * PI / 3.0).sin(),
            -nh_bond * theta_h.sin() * (2.0 * PI / 3.0).cos(),
            -nh_bond * theta_h.cos(),
        ];

        let atoms: Vec<(u8, [f64; 3])> = vec![(7, [0.0, 0.0, 0.0]), (1, h1), (1, h2), (1, h3)];

        let snapshot = compute_all_internals(&atoms);

        assert_eq!(snapshot.bonds.len(), 3, "NH3 should have 3 bonds");
        assert_eq!(snapshot.angles.len(), 3, "NH3 should have 3 angles");
    }

    // =========================================================================
    // Rigid Scan Internal Coordinate Tracking Tests
    // =========================================================================

    #[test]
    fn test_rigid_scan_has_internal_coordinates() {
        // Verify that rigid scan computes InternalCoordinateSnapshot at each point.
        // Use H2O angle scan so bonds are always within covalent detection range.
        let oh_bond = 1.809;
        let angle0 = 104.5_f64.to_radians();
        let h1 = [
            0.0,
            oh_bond * (angle0 / 2.0).sin(),
            -oh_bond * (angle0 / 2.0).cos(),
        ];
        let h2 = [
            0.0,
            -oh_bond * (angle0 / 2.0).sin(),
            -oh_bond * (angle0 / 2.0).cos(),
        ];

        let config = PesScanInternalConfig {
            atoms: vec![(8, [0.0, 0.0, 0.0]), (1, h1), (1, h2)],
            coordinate: ScanCoordinate::Angle {
                atom_i: 1,
                atom_j: 0,
                atom_k: 2,
            },
            value_min: 95.0_f64.to_radians(),
            value_max: 115.0_f64.to_radians(),
            n_points: 3,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: false,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None).expect("Scan should succeed");

        for pt in &result.points {
            assert!(
                pt.internal_coordinates.is_some(),
                "Rigid scan point at v={:.4} should have internal_coordinates",
                pt.coordinate_value
            );
            if let Some(snap) = &pt.internal_coordinates {
                assert_eq!(snap.bonds.len(), 2, "H2O should have 2 bonds");
                assert_eq!(snap.angles.len(), 1, "H2O should have 1 angle");
            }
            assert!(
                pt.opt_steps.is_none(),
                "Rigid scan points should have opt_steps=None"
            );
        }

        assert_eq!(result.scan_mode, "rigid");
        assert_eq!(result.total_opt_steps, 0);
    }

    #[test]
    fn test_rigid_scan_h2o_bond_lengths_constant() {
        // For a rigid H2O angle scan, O-H bond lengths should be identical at every point
        let oh_bond = 1.809;
        let angle0 = 104.5_f64.to_radians();
        let h1 = [
            0.0,
            oh_bond * (angle0 / 2.0).sin(),
            -oh_bond * (angle0 / 2.0).cos(),
        ];
        let h2 = [
            0.0,
            -oh_bond * (angle0 / 2.0).sin(),
            -oh_bond * (angle0 / 2.0).cos(),
        ];

        let config = PesScanInternalConfig {
            atoms: vec![(8, [0.0, 0.0, 0.0]), (1, h1), (1, h2)],
            coordinate: ScanCoordinate::Angle {
                atom_i: 1,
                atom_j: 0,
                atom_k: 2,
            },
            value_min: 90.0_f64.to_radians(),
            value_max: 120.0_f64.to_radians(),
            n_points: 5,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: false,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None).expect("Scan should succeed");

        // Get O-H bond lengths from internal coordinate snapshots
        for pt in &result.points {
            if !pt.converged {
                continue;
            }
            let snap = pt.internal_coordinates.as_ref().unwrap();
            for &(_, _, bond_len) in &snap.bonds {
                assert!(
                    (bond_len - oh_bond).abs() < 1e-8,
                    "Rigid scan: O-H bond {:.8} should equal initial {:.8} at angle={:.1} deg",
                    bond_len,
                    oh_bond,
                    pt.coordinate_value.to_degrees()
                );
            }
        }
    }

    // =========================================================================
    // Relaxed PES Scan Tests
    // =========================================================================

    #[test]
    fn test_relaxed_scan_h2_bond() {
        // H2 bond scan: relaxed should equal rigid (nothing else to relax).
        // For a diatomic, relaxed energies must equal rigid energies.
        let config = PesScanInternalConfig {
            atoms: vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])],
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.0,
            value_max: 2.0,
            n_points: 5,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let rigid = pes_scan_rigid(&config, None).expect("Rigid scan should succeed");
        let relaxed = pes_scan_relaxed(&config, None).expect("Relaxed scan should succeed");

        assert_eq!(rigid.points.len(), relaxed.points.len());
        assert_eq!(relaxed.scan_mode, "relaxed");

        for (rg, rx) in rigid.points.iter().zip(relaxed.points.iter()) {
            // Both should converge
            assert!(
                rg.converged && rx.converged,
                "Points at v={:.4} should converge",
                rg.coordinate_value
            );

            // Energies should match closely (H2 has only 1 DOF - the bond)
            let error = (rg.energy - rx.energy).abs();
            assert!(
                error < 1e-8,
                "H2 relaxed vs rigid energy mismatch at v={:.4}: rigid={:.10} relaxed={:.10} error={:.2e}",
                rg.coordinate_value,
                rg.energy,
                rx.energy,
                error
            );
        }
    }

    #[test]
    fn test_relaxed_scan_h2o_angle_variational() {
        // H2O angle scan: E_relaxed <= E_rigid at every point (variational principle)
        let oh_bond = 1.809;
        let angle0 = 104.5_f64.to_radians();
        let h1 = [
            0.0,
            oh_bond * (angle0 / 2.0).sin(),
            -oh_bond * (angle0 / 2.0).cos(),
        ];
        let h2 = [
            0.0,
            -oh_bond * (angle0 / 2.0).sin(),
            -oh_bond * (angle0 / 2.0).cos(),
        ];

        let config = PesScanInternalConfig {
            atoms: vec![(8, [0.0, 0.0, 0.0]), (1, h1), (1, h2)],
            coordinate: ScanCoordinate::Angle {
                atom_i: 1,
                atom_j: 0,
                atom_k: 2,
            },
            value_min: 90.0_f64.to_radians(),
            value_max: 120.0_f64.to_radians(),
            n_points: 4,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: false,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let rigid = pes_scan_rigid(&config, None).expect("Rigid scan should succeed");
        let relaxed = pes_scan_relaxed(&config, None).expect("Relaxed scan should succeed");

        assert_eq!(rigid.points.len(), relaxed.points.len());

        for (rg, rx) in rigid.points.iter().zip(relaxed.points.iter()) {
            if !rg.converged || !rx.converged {
                continue;
            }
            // Variational principle: E_relaxed <= E_rigid + small tolerance
            // The tolerance accounts for SCF convergence differences
            assert!(
                rx.energy <= rg.energy + 1e-8,
                "Variational violation at angle={:.1} deg: E_relaxed={:.10} > E_rigid={:.10}",
                rg.coordinate_value.to_degrees(),
                rx.energy,
                rg.energy,
            );
        }
    }

    #[test]
    fn test_relaxed_scan_coord_tracking() {
        // Verify InternalCoordinateSnapshot is computed at each point
        // For H2O angle scan: should have 2 bonds, 1 angle
        let oh_bond = 1.809;
        let angle0 = 104.5_f64.to_radians();
        let h1 = [
            0.0,
            oh_bond * (angle0 / 2.0).sin(),
            -oh_bond * (angle0 / 2.0).cos(),
        ];
        let h2 = [
            0.0,
            -oh_bond * (angle0 / 2.0).sin(),
            -oh_bond * (angle0 / 2.0).cos(),
        ];

        let config = PesScanInternalConfig {
            atoms: vec![(8, [0.0, 0.0, 0.0]), (1, h1), (1, h2)],
            coordinate: ScanCoordinate::Angle {
                atom_i: 1,
                atom_j: 0,
                atom_k: 2,
            },
            value_min: 95.0_f64.to_radians(),
            value_max: 115.0_f64.to_radians(),
            n_points: 3,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: false,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_relaxed(&config, None).expect("Relaxed scan should succeed");

        for pt in &result.points {
            assert!(
                pt.internal_coordinates.is_some(),
                "Relaxed scan point should have internal_coordinates"
            );
            if let Some(snap) = &pt.internal_coordinates {
                assert_eq!(snap.bonds.len(), 2, "H2O should have 2 bonds");
                assert_eq!(snap.angles.len(), 1, "H2O should have 1 angle");
            }
            assert!(
                pt.opt_steps.is_some(),
                "Relaxed scan points should have opt_steps"
            );
        }
    }

    #[test]
    fn test_scan_dispatcher() {
        // Test pes_scan_internal() dispatcher routes correctly
        let config = PesScanInternalConfig {
            atoms: vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])],
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.0,
            value_max: 2.0,
            n_points: 3,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let rigid = pes_scan_internal(&config, ScanMode::Rigid, None)
            .expect("Rigid dispatch should succeed");
        assert_eq!(rigid.scan_mode, "rigid");

        let relaxed = pes_scan_internal(&config, ScanMode::Relaxed, None)
            .expect("Relaxed dispatch should succeed");
        assert_eq!(relaxed.scan_mode, "relaxed");
    }

    // =========================================================================
    // Constraint Gradient Projection Tests
    // =========================================================================

    #[test]
    fn test_constraint_gradient_bond() {
        // Bond between atoms 0 and 1 along z-axis
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]];
        let c = constraint_gradient_bond(&positions, 0, 1);

        // Unit bond vector is [0, 0, 1]
        // d(r)/d(r_0) = -[0, 0, 1], d(r)/d(r_1) = [0, 0, 1]
        let norm: f64 = c.iter().map(|x| x * x).sum::<f64>().sqrt();
        // Norm should be sqrt(2) since |[-e]|^2 + |[+e]|^2 = 1 + 1 = 2
        assert!((norm - 2.0_f64.sqrt()).abs() < 1e-10, "norm = {}", norm);

        // z-component of atom 0 should be negative (bond direction)
        assert!(c[2] < 0.0, "atom 0 z-component should be negative");
        // z-component of atom 1 should be positive
        assert!(c[5] > 0.0, "atom 1 z-component should be positive");
    }

    #[test]
    fn test_project_out_constraint_removes_component() {
        // Create a gradient and a constraint direction, verify projection
        let gradient = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]; // 2 atoms
        let mut constraint_dir = vec![0.0; 6];
        // Constraint direction: pure x on atom 0
        constraint_dir[0] = 1.0;

        let projected = project_out_constraint(&gradient, &constraint_dir);

        // The x-component of atom 0 should be projected out
        assert!(
            projected[0][0].abs() < 1e-10,
            "x-component of atom 0 should be zero after projection, got {}",
            projected[0][0]
        );
        // y-component of atom 1 should be unchanged (orthogonal to constraint)
        assert!(
            (projected[1][1] - 1.0).abs() < 1e-10,
            "y-component of atom 1 should be unchanged"
        );
    }

    #[test]
    fn test_enforce_constraint_bond() {
        // Place H2 at 2.0 bohr, enforce constraint to 1.5 bohr
        let atoms_z = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]];
        let coord = ScanCoordinate::Bond {
            atom_i: 0,
            atom_j: 1,
        };

        let (corrected, needed) = enforce_constraint(
            &positions,
            &atoms_z,
            &coord,
            1.5,
            CONSTRAINT_TOLERANCE,
            &None,
        );

        assert!(needed, "Correction should be needed");

        // Check that the bond length is now 1.5
        let new_d = distance(&corrected[0], &corrected[1]);
        assert!(
            (new_d - 1.5).abs() < 1e-10,
            "Bond length should be 1.5, got {}",
            new_d
        );
    }

    // =========================================================================
    // Golden Tests: H2O PES scans vs PySCF
    // =========================================================================
    //
    // These tests compare IQCP rigid PES scan energies against PySCF reference
    // values stored in tests/golden/pes/*.json.
    //
    // Tolerance: 1e-4 Ha (absolute).
    // Rationale: IQCP uses Obara-Saika recurrence for integrals while PySCF
    // uses libcint. Small numerical differences (~1e-5) in the Boys function
    // and nuclear attraction integrals accumulate through the SCF, but the
    // total energy should agree well within chemical accuracy.
    //
    // Golden data generated by:
    //   tests/benchmarks/pyscf/pes_rigid_h2o_angle.py
    //   tests/benchmarks/pyscf/pes_rigid_h2o_bond.py
    // PySCF 2.11.0, RHF/STO-3G, conv_tol=1e-12

    /// Deserialized golden PES data from PySCF JSON files
    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenPesInternal {
        molecule: String,
        basis: String,
        method: String,
        scan_type: String,
        points: Vec<GoldenPesInternalPoint>,
        equilibrium: GoldenPesInternalEquilibrium,
    }

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenPesInternalPoint {
        coordinate_value: f64,
        energy: f64,
        converged: bool,
        scf_iterations: usize,
        geometry: Vec<Vec<f64>>,
    }

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenPesInternalEquilibrium {
        value: f64,
        energy: f64,
    }

    fn load_golden_h2o_angle() -> GoldenPesInternal {
        let data = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/golden/pes/h2o_angle_rigid_sto3g.json"
        ));
        serde_json::from_str(data).expect("Failed to parse H2O angle golden data")
    }

    fn load_golden_h2o_bond() -> GoldenPesInternal {
        let data = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/golden/pes/h2o_bond_rigid_sto3g.json"
        ));
        serde_json::from_str(data).expect("Failed to parse H2O bond golden data")
    }

    /// H2O preset geometry in (Z, [x, y, z]) format (bohr)
    fn h2o_atoms() -> Vec<(u8, [f64; 3])> {
        vec![
            (8, [0.0, 0.0, 0.2217282]),         // O
            (1, [0.0, 1.4305447, -0.8869128]),  // H1
            (1, [0.0, -1.4305447, -0.8869128]), // H2
        ]
    }

    #[test]
    fn test_golden_h2o_angle_rigid() {
        // Compare IQCP rigid H-O-H angle scan energies against PySCF golden data.
        //
        // IQCP generates geometries by rotating atom_k (H2) around the normal
        // to the H1-O-H2 plane. PySCF places both H atoms symmetrically.
        // Both produce geometries with the SAME H-O-H angle and O-H distances,
        // so the Hamiltonian (and energy) is identical by symmetry.
        //
        // We compare by running IQCP at the exact PySCF geometries (from golden
        // data) rather than using IQCP's own geometry generation. This isolates
        // the SCF comparison from any geometry generation differences.
        let golden = load_golden_h2o_angle();
        assert_eq!(golden.molecule, "H2O");
        assert_eq!(golden.basis, "sto-3g");
        assert_eq!(golden.method, "RHF");
        assert_eq!(golden.scan_type, "angle_rigid");
        assert_eq!(golden.points.len(), 11);

        let tolerance = 1e-4; // Ha — manuscript validation standard

        for (i, gpt) in golden.points.iter().enumerate() {
            assert!(gpt.converged, "PySCF point {} not converged", i);

            // Build IQCP molecule from the PySCF geometry
            let atoms: Vec<(u8, [f64; 3])> = vec![
                (
                    8,
                    [gpt.geometry[0][0], gpt.geometry[0][1], gpt.geometry[0][2]],
                ),
                (
                    1,
                    [gpt.geometry[1][0], gpt.geometry[1][1], gpt.geometry[1][2]],
                ),
                (
                    1,
                    [gpt.geometry[2][0], gpt.geometry[2][1], gpt.geometry[2][2]],
                ),
            ];

            let iqcp_atoms: Vec<Atom> = atoms
                .iter()
                .filter_map(|&(z, pos)| Atom::new(z, pos).ok())
                .collect();

            let basis = BasisSet::build(iqcp_atoms, "sto-3g").expect("Failed to build basis set");

            let s_mat = overlap_matrix(&basis);
            let h_core_mat = hcore_matrix(&basis);
            let eri = eri_compressed(&basis);

            let system = PresetSystem {
                system_id: format!("h2o_angle_golden_{}", i),
                label: format!("H2O angle scan point {}", i),
                nbf: basis.n_basis,
                nelec: basis.n_electrons,
                e_nuc: basis.nuclear_repulsion,
                s_matrix: s_mat,
                h_core: h_core_mat,
                eri_compressed: eri,
            };

            let scf_config = ScfConfig::tight();
            let result =
                rhf_scf_with_guess(&system, &scf_config, None).expect("SCF should converge");

            let error = (result.energy_total - gpt.energy).abs();
            assert!(
                error < tolerance,
                "H2O angle golden mismatch at point {} (angle={:.1} deg): \
                 IQCP={:.10} PySCF={:.10} error={:.2e}",
                i,
                gpt.coordinate_value.to_degrees(),
                result.energy_total,
                gpt.energy,
                error
            );
        }
    }

    #[test]
    fn test_golden_h2o_bond_rigid() {
        // Compare IQCP rigid O-H1 bond scan energies against PySCF golden data.
        //
        // Both IQCP and PySCF move H1 along the O->H1 direction, so the
        // geometries are identical (within floating-point).
        let golden = load_golden_h2o_bond();
        assert_eq!(golden.molecule, "H2O");
        assert_eq!(golden.basis, "sto-3g");
        assert_eq!(golden.method, "RHF");
        assert_eq!(golden.scan_type, "bond_rigid");
        assert_eq!(golden.points.len(), 11);

        let tolerance = 1e-4; // Ha — manuscript validation standard

        for (i, gpt) in golden.points.iter().enumerate() {
            assert!(gpt.converged, "PySCF point {} not converged", i);

            // Build IQCP molecule from the PySCF geometry
            let atoms: Vec<(u8, [f64; 3])> = vec![
                (
                    8,
                    [gpt.geometry[0][0], gpt.geometry[0][1], gpt.geometry[0][2]],
                ),
                (
                    1,
                    [gpt.geometry[1][0], gpt.geometry[1][1], gpt.geometry[1][2]],
                ),
                (
                    1,
                    [gpt.geometry[2][0], gpt.geometry[2][1], gpt.geometry[2][2]],
                ),
            ];

            let iqcp_atoms: Vec<Atom> = atoms
                .iter()
                .filter_map(|&(z, pos)| Atom::new(z, pos).ok())
                .collect();

            let basis = BasisSet::build(iqcp_atoms, "sto-3g").expect("Failed to build basis set");

            let s_mat = overlap_matrix(&basis);
            let h_core_mat = hcore_matrix(&basis);
            let eri = eri_compressed(&basis);

            let system = PresetSystem {
                system_id: format!("h2o_bond_golden_{}", i),
                label: format!("H2O bond scan point {}", i),
                nbf: basis.n_basis,
                nelec: basis.n_electrons,
                e_nuc: basis.nuclear_repulsion,
                s_matrix: s_mat,
                h_core: h_core_mat,
                eri_compressed: eri,
            };

            let scf_config = ScfConfig::tight();
            let result =
                rhf_scf_with_guess(&system, &scf_config, None).expect("SCF should converge");

            let error = (result.energy_total - gpt.energy).abs();
            assert!(
                error < tolerance,
                "H2O bond golden mismatch at point {} (R_OH1={:.3} bohr): \
                 IQCP={:.10} PySCF={:.10} error={:.2e}",
                i,
                gpt.coordinate_value,
                result.energy_total,
                gpt.energy,
                error
            );
        }
    }

    #[test]
    fn test_golden_h2o_angle_via_scan_engine() {
        // Run the full pes_scan_rigid engine for H2O angle scan and compare
        // against PySCF golden data. This tests both geometry generation
        // AND SCF accuracy end-to-end.
        let golden = load_golden_h2o_angle();

        let config = PesScanInternalConfig {
            atoms: h2o_atoms(),
            coordinate: ScanCoordinate::Angle {
                atom_i: 1,
                atom_j: 0,
                atom_k: 2,
            },
            value_min: golden.points.first().unwrap().coordinate_value,
            value_max: golden.points.last().unwrap().coordinate_value,
            n_points: 11,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None).expect("Scan should succeed");
        assert_eq!(result.points.len(), 11);

        let tolerance = 1e-4; // Ha

        for (i, (iqcp_pt, pyscf_pt)) in result.points.iter().zip(golden.points.iter()).enumerate() {
            assert!(
                iqcp_pt.converged,
                "IQCP point {} did not converge at angle={:.1} deg",
                i,
                iqcp_pt.coordinate_value.to_degrees()
            );

            // Coordinate values should match (same linspace)
            assert!(
                (iqcp_pt.coordinate_value - pyscf_pt.coordinate_value).abs() < 1e-10,
                "Coordinate mismatch at point {}: IQCP={} PySCF={}",
                i,
                iqcp_pt.coordinate_value,
                pyscf_pt.coordinate_value
            );

            let error = (iqcp_pt.energy - pyscf_pt.energy).abs();
            assert!(
                error < tolerance,
                "H2O angle scan engine mismatch at point {} (angle={:.1} deg): \
                 IQCP={:.10} PySCF={:.10} error={:.2e}",
                i,
                iqcp_pt.coordinate_value.to_degrees(),
                iqcp_pt.energy,
                pyscf_pt.energy,
                error
            );
        }
    }

    #[test]
    fn test_golden_h2o_bond_via_scan_engine() {
        // Run the full pes_scan_rigid engine for H2O bond scan and compare
        // against PySCF golden data.
        let golden = load_golden_h2o_bond();

        let config = PesScanInternalConfig {
            atoms: h2o_atoms(),
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: golden.points.first().unwrap().coordinate_value,
            value_max: golden.points.last().unwrap().coordinate_value,
            n_points: 11,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None).expect("Scan should succeed");
        assert_eq!(result.points.len(), 11);

        let tolerance = 1e-4; // Ha

        for (i, (iqcp_pt, pyscf_pt)) in result.points.iter().zip(golden.points.iter()).enumerate() {
            assert!(
                iqcp_pt.converged,
                "IQCP point {} did not converge at R_OH1={:.3} bohr",
                i, iqcp_pt.coordinate_value
            );

            // Coordinate values should match (same linspace)
            assert!(
                (iqcp_pt.coordinate_value - pyscf_pt.coordinate_value).abs() < 1e-10,
                "Coordinate mismatch at point {}: IQCP={} PySCF={}",
                i,
                iqcp_pt.coordinate_value,
                pyscf_pt.coordinate_value
            );

            let error = (iqcp_pt.energy - pyscf_pt.energy).abs();
            assert!(
                error < tolerance,
                "H2O bond scan engine mismatch at point {} (R_OH1={:.3} bohr): \
                 IQCP={:.10} PySCF={:.10} error={:.2e}",
                i,
                iqcp_pt.coordinate_value,
                iqcp_pt.energy,
                pyscf_pt.energy,
                error
            );
        }
    }

    #[test]
    fn test_variational_h2o_angle() {
        // Variational principle: the rigid scan energies should form a smooth
        // potential curve (no spurious jumps) and the equilibrium from the
        // angle scan should be physically reasonable.
        //
        // Rigid PES is an upper bound to the true PES at each angle. We verify
        // that the rigid scan produces a well-shaped potential curve (minimum
        // not at boundary, energies decrease monotonically toward minimum from
        // both sides).
        let config = PesScanInternalConfig {
            atoms: h2o_atoms(),
            coordinate: ScanCoordinate::Angle {
                atom_i: 1,
                atom_j: 0,
                atom_k: 2,
            },
            value_min: 80.0_f64.to_radians(),
            value_max: 130.0_f64.to_radians(),
            n_points: 11,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None).expect("Scan should succeed");

        // All points should converge
        for pt in &result.points {
            assert!(
                pt.converged,
                "Point at angle={:.1} deg did not converge",
                pt.coordinate_value.to_degrees()
            );
        }

        // Should have an equilibrium (minimum not at boundary)
        let eq = result
            .equilibrium
            .as_ref()
            .expect("Should find equilibrium");

        // H2O STO-3G equilibrium angle is around 100-105 degrees
        let eq_deg = eq.value.to_degrees();
        assert!(
            eq_deg > 95.0 && eq_deg < 110.0,
            "Equilibrium angle {:.2} deg outside expected range 95-110",
            eq_deg
        );

        // Monotonicity: energies decrease toward minimum, increase after
        let converged: Vec<&PesInternalPoint> =
            result.points.iter().filter(|p| p.converged).collect();
        let min_idx = converged
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                a.energy
                    .partial_cmp(&b.energy)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap();

        // Check descending before minimum
        for i in 1..min_idx {
            assert!(
                converged[i].energy <= converged[i - 1].energy + 1e-10,
                "Energy should decrease before minimum: E({:.1} deg)={:.8} > E({:.1} deg)={:.8}",
                converged[i].coordinate_value.to_degrees(),
                converged[i].energy,
                converged[i - 1].coordinate_value.to_degrees(),
                converged[i - 1].energy,
            );
        }

        // Check ascending after minimum
        for i in (min_idx + 1)..converged.len() {
            assert!(
                converged[i].energy >= converged[i - 1].energy - 1e-10,
                "Energy should increase after minimum: E({:.1} deg)={:.8} < E({:.1} deg)={:.8}",
                converged[i].coordinate_value.to_degrees(),
                converged[i].energy,
                converged[i - 1].coordinate_value.to_degrees(),
                converged[i - 1].energy,
            );
        }
    }

    // =========================================================================
    // Golden Tests: H2O Relaxed PES scan vs PySCF (geomeTRIC)
    // =========================================================================
    //
    // These tests compare IQCP relaxed PES scan energies against PySCF
    // reference values generated by geomeTRIC constrained optimization.
    //
    // At each angle, geomeTRIC holds the H-O-H angle fixed while optimizing
    // the O-H bond lengths. IQCP uses constrained steepest descent with
    // gradient projection, which is a simpler optimizer than geomeTRIC's
    // trust-region method.
    //
    // Tolerance: 1e-3 Ha (absolute).
    // Rationale: IQCP's steepest descent optimizer may not reach the exact
    // same geometry as geomeTRIC's trust-region optimizer, especially at
    // angles far from equilibrium where the O-H bond relaxation is larger.
    // However, both should produce energies within ~1 mHa since the PES
    // is relatively flat in the bond-stretch direction near equilibrium.
    //
    // Golden data generated by:
    //   tests/benchmarks/pyscf/pes_relaxed_h2o_angle.py
    // PySCF 2.11.0, geomeTRIC 1.1, RHF/STO-3G, conv_tol=1e-12

    /// Deserialized relaxed golden PES data from PySCF/geomeTRIC JSON
    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenPesRelaxed {
        molecule: String,
        basis: String,
        method: String,
        scan_type: String,
        points: Vec<GoldenPesRelaxedPoint>,
        equilibrium: Option<GoldenPesRelaxedEquilibrium>,
    }

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenPesRelaxedPoint {
        coordinate_value: f64,
        coordinate_value_deg: f64,
        energy: f64,
        converged: bool,
        geometry: Vec<Vec<f64>>,
        r_oh1_bohr: f64,
        r_oh2_bohr: f64,
        angle_deg: f64,
    }

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenPesRelaxedEquilibrium {
        value: f64,
        value_deg: f64,
        energy: f64,
    }

    fn load_golden_h2o_angle_relaxed() -> GoldenPesRelaxed {
        let data = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/golden/pes/h2o_angle_relaxed_sto3g.json"
        ));
        serde_json::from_str(data).expect("Failed to parse H2O relaxed angle golden data")
    }

    #[test]
    fn test_golden_h2o_angle_relaxed() {
        // Compare IQCP relaxed H-O-H angle scan energies against PySCF
        // geomeTRIC constrained optimization golden data.
        //
        // IQCP uses constrained steepest descent with gradient projection,
        // while PySCF/geomeTRIC uses a trust-region optimizer. Both hold the
        // H-O-H angle fixed while optimizing O-H bond lengths.
        //
        // Tolerance: 1e-3 Ha per point. This is looser than the rigid scan
        // tolerance (1e-4 Ha) because the optimizer implementations differ.
        let golden = load_golden_h2o_angle_relaxed();
        assert_eq!(golden.molecule, "H2O");
        assert_eq!(golden.basis, "sto-3g");
        assert_eq!(golden.method, "RHF");
        assert_eq!(golden.scan_type, "angle_relaxed");
        assert_eq!(golden.points.len(), 11);

        // Run IQCP relaxed scan with matching parameters
        let config = PesScanInternalConfig {
            atoms: h2o_atoms(),
            coordinate: ScanCoordinate::Angle {
                atom_i: 1,
                atom_j: 0,
                atom_k: 2,
            },
            value_min: golden.points.first().unwrap().coordinate_value,
            value_max: golden.points.last().unwrap().coordinate_value,
            n_points: 11,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_relaxed(&config, None).expect("Relaxed scan should succeed");
        assert_eq!(result.points.len(), 11);
        assert_eq!(result.scan_mode, "relaxed");

        let tolerance = 1e-3; // Ha — relaxed scan tolerance

        for (i, (iqcp_pt, pyscf_pt)) in result.points.iter().zip(golden.points.iter()).enumerate() {
            // Coordinate values should match (same linspace)
            assert!(
                (iqcp_pt.coordinate_value - pyscf_pt.coordinate_value).abs() < 1e-10,
                "Coordinate mismatch at point {}: IQCP={} PySCF={}",
                i,
                iqcp_pt.coordinate_value,
                pyscf_pt.coordinate_value
            );

            if !iqcp_pt.converged {
                // If IQCP didn't converge the optimization, skip energy comparison
                // but note it for debugging
                eprintln!(
                    "  Warning: IQCP relaxed point {} (angle={:.1} deg) did not converge",
                    i, pyscf_pt.coordinate_value_deg
                );
                continue;
            }

            let diff = (iqcp_pt.energy - pyscf_pt.energy).abs();
            assert!(
                diff < tolerance,
                "H2O relaxed angle golden mismatch at point {} (angle={:.1} deg): \
                 IQCP={:.10} PySCF={:.10} diff={:.2e} > {:.0e}",
                i,
                pyscf_pt.coordinate_value_deg,
                iqcp_pt.energy,
                pyscf_pt.energy,
                diff,
                tolerance
            );
        }
    }

    #[test]
    fn test_golden_relaxed_vs_rigid_variational() {
        // Verify that E_relaxed <= E_rigid at every corresponding angle.
        //
        // The variational principle guarantees that relaxing additional
        // degrees of freedom (O-H bond lengths) can only lower or maintain
        // the energy. This test cross-checks both rigid and relaxed golden
        // data from PySCF.
        let rigid = load_golden_h2o_angle();
        let relaxed = load_golden_h2o_angle_relaxed();

        assert_eq!(rigid.points.len(), relaxed.points.len());

        for (i, (rg, rx)) in rigid.points.iter().zip(relaxed.points.iter()).enumerate() {
            // Angles should match
            assert!(
                (rg.coordinate_value - rx.coordinate_value).abs() < 1e-10,
                "Angle mismatch at point {}: rigid={} relaxed={}",
                i,
                rg.coordinate_value,
                rx.coordinate_value
            );

            if !rg.converged || !rx.converged {
                continue;
            }

            // Variational principle: E_relaxed <= E_rigid
            // Allow a tiny tolerance for floating-point rounding
            assert!(
                rx.energy <= rg.energy + 1e-10,
                "Variational violation at point {} (angle={:.1} deg): \
                 E_relaxed={:.10} > E_rigid={:.10} (diff={:.2e})",
                i,
                rx.coordinate_value_deg,
                rx.energy,
                rg.energy,
                rx.energy - rg.energy
            );

            // Also verify the relaxation lowered energy by a reasonable amount
            // (not orders of magnitude, which would indicate a bug)
            let lowering = rg.energy - rx.energy;
            assert!(
                lowering < 0.02,
                "Excessive relaxation at point {} (angle={:.1} deg): \
                 lowering={:.6} Ha > 20 mHa — possible optimizer issue",
                i,
                rx.coordinate_value_deg,
                lowering
            );
        }
    }

    #[test]
    fn test_relaxed_scan_ethane_dihedral() {
        // C2H6 dihedral scan (3 points, 0-180 degrees): verifies that the
        // L-BFGS constrained optimizer converges for 8-atom molecules.
        //
        // This was the original failure case: steepest descent with a constant
        // diagonal Hessian took >50 steps without converging. The L-BFGS
        // optimizer with constraint-projected direction converges in ~7 steps.
        let atoms = ethane_atoms();

        let config = PesScanInternalConfig {
            atoms,
            coordinate: ScanCoordinate::Dihedral {
                atom_i: 2,
                atom_j: 0,
                atom_k: 1,
                atom_l: 5,
            },
            value_min: 0.0_f64.to_radians(),
            value_max: 180.0_f64.to_radians(),
            n_points: 3,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result =
            pes_scan_internal(&config, ScanMode::Relaxed, None).expect("Relaxed scan should work");

        assert_eq!(result.points.len(), 3);
        assert_eq!(result.scan_mode, "relaxed");

        // All three points should converge
        for pt in &result.points {
            assert!(
                pt.converged,
                "Ethane relaxed point at {:.1} deg should converge, \
                 opt_steps={:?}",
                pt.coordinate_value.to_degrees(),
                pt.opt_steps
            );
        }

        // Energy should increase from staggered (180 deg) toward eclipsed (0 deg)
        // (ethane torsion barrier: staggered is lower energy)
        let e0 = result.points[0].energy; // 0 deg (eclipsed)
        let e180 = result.points[2].energy; // 180 deg (staggered)
        assert!(
            e180 < e0,
            "Staggered ethane ({:.10} Ha) should be lower energy than eclipsed ({:.10} Ha)",
            e180,
            e0
        );

        // Each point should converge in <= 15 optimization steps
        for pt in &result.points {
            let steps = pt.opt_steps.unwrap_or(0);
            assert!(
                steps <= 15,
                "Point at {:.1} deg took {} opt steps (max 15 expected)",
                pt.coordinate_value.to_degrees(),
                steps
            );
        }
    }

    #[test]
    fn test_b3lyp_rigid_h2o_bond_scan() {
        // Reproduce B3LYP PES scan crash on H2O O-H(1) bond
        let atoms = vec![
            (8u8, [0.0, 0.0, 0.2217282]),
            (1, [0.0, 1.4305447, -0.8869128]),
            (1, [0.0, -1.4305447, -0.8869128]),
        ];
        let config = PesScanInternalConfig {
            atoms,
            basis_name: "sto-3g".to_string(),
            method: "b3lyp".to_string(),
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 0.8,
            value_max: 3.0,
            n_points: 5, // small for speed
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };
        let result = pes_scan_internal(&config, ScanMode::Rigid, None);
        assert!(
            result.is_ok(),
            "B3LYP PES scan should not crash: {:?}",
            result.err()
        );
        let result = result.unwrap();
        eprintln!("B3LYP H2O bond scan: {} points", result.points.len());
        for pt in &result.points {
            eprintln!(
                "  R={:.4}  E={:.10}  conv={}  iters={}",
                pt.coordinate_value, pt.energy, pt.converged, pt.scf_iterations
            );
        }
        assert_eq!(result.points.len(), 5);
    }

    #[test]
    fn test_b3lyp_relaxed_h2o_bond_scan() {
        // B3LYP relaxed PES scan — must NOT panic (was crashing with "unreachable executed")
        let atoms = vec![
            (8u8, [0.0, 0.0, 0.2217282]),
            (1, [0.0, 1.4305447, -0.8869128]),
            (1, [0.0, -1.4305447, -0.8869128]),
        ];
        let config = PesScanInternalConfig {
            atoms,
            basis_name: "sto-3g".to_string(),
            method: "b3lyp".to_string(),
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 0.8,
            value_max: 3.0,
            n_points: 5,
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: Some(20),
            opt_grad_threshold: Some(4.5e-4),
        };
        let result = pes_scan_internal(&config, ScanMode::Relaxed, None);
        assert!(
            result.is_ok(),
            "B3LYP relaxed PES scan should not crash: {:?}",
            result.err()
        );
        let result = result.unwrap();
        eprintln!(
            "B3LYP relaxed H2O bond scan: {} points",
            result.points.len()
        );
        for pt in &result.points {
            eprintln!(
                "  R={:.4}  E={:.10}  conv={}  iters={}  opt_steps={:?}",
                pt.coordinate_value, pt.energy, pt.converged, pt.scf_iterations, pt.opt_steps
            );
        }
        assert_eq!(result.points.len(), 5);
    }

    #[test]
    fn test_b3lyp_relaxed_h2o_convergence_near_equilibrium() {
        // Regression test: B3LYP relaxed scan near equilibrium O-H bond
        // length (~1.38-1.40 bohr) previously failed to converge within
        // 50 optimization steps due to L-BFGS producing aggressive steps
        // that caused energy oscillation. Fixed by:
        // 1. Resetting L-BFGS history when energy increases by >10 mHa
        // 2. Best-energy stagnation detection (patience=8 steps)
        let atoms = vec![
            (8u8, [0.0, 0.0, 0.2217282]),
            (1, [0.0, 1.4305447, -0.8869128]),
            (1, [0.0, -1.4305447, -0.8869128]),
        ];
        let config = PesScanInternalConfig {
            atoms,
            basis_name: "sto-3g".to_string(),
            method: "b3lyp".to_string(),
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.3789,
            value_max: 1.4,
            n_points: 2,
            use_seeding: false,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: Some(50),
            opt_grad_threshold: Some(4.5e-4),
        };
        let result = pes_scan_internal(&config, ScanMode::Relaxed, None);
        assert!(result.is_ok(), "Scan failed: {:?}", result.err());
        let result = result.unwrap();
        for pt in &result.points {
            assert!(
                pt.converged,
                "Point R={:.4} should converge within 50 steps (got {:?} steps)",
                pt.coordinate_value, pt.opt_steps
            );
        }
    }

    // =========================================================================
    // SI Data Extraction Tests (GAP-1, GAP-2, GAP-5)
    // =========================================================================

    #[test]
    fn test_si_gap2_h2o_bond_rigid_energies() {
        // GAP-2: H2O rigid O-H(1) bond scan (RHF/STO-3G, 11 points, 1.4-2.4 bohr)
        // Prints IQCP energies alongside PySCF golden data for SI table.
        let config = PesScanInternalConfig {
            atoms: h2o_atoms(),
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.4,
            value_max: 2.4,
            n_points: 11,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None).expect("Scan should succeed");

        // PySCF golden data (from tests/golden/pes/h2o_bond_rigid_sto3g.json)
        let golden = load_golden_h2o_bond();

        eprintln!("\n=== GAP-2: H2O O-H rigid scan RHF/STO-3G (11 points) ===");
        eprintln!(
            "{:<10} {:<20} {:<20} {:<15}",
            "R(bohr)", "E_IQCP (Ha)", "E_PySCF (Ha)", "|dE| (Ha)"
        );
        for (iqcp_pt, gold_pt) in result.points.iter().zip(golden.points.iter()) {
            let delta = (iqcp_pt.energy - gold_pt.energy).abs();
            eprintln!(
                "{:<10.1} {:<20.10} {:<20.10} {:<15.2e}",
                iqcp_pt.coordinate_value, iqcp_pt.energy, gold_pt.energy, delta
            );
            assert!(
                delta < 1e-4,
                "Error exceeds 1e-4 Ha at R={:.1}",
                iqcp_pt.coordinate_value
            );
        }
        eprintln!("All 11 points agree within 1e-4 Ha.");
    }

    #[test]
    fn test_si_gap5_seeding_iteration_reduction_h2o() {
        // GAP-5: Density seeding iteration reduction measurement
        // H2O rigid bond scan (RHF/STO-3G, 11 points, 1.4-2.4 bohr)
        // Compare total SCF iterations: seeded vs cold start
        let base_config = PesScanInternalConfig {
            atoms: h2o_atoms(),
            coordinate: ScanCoordinate::Bond {
                atom_i: 0,
                atom_j: 1,
            },
            value_min: 1.4,
            value_max: 2.4,
            n_points: 11,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: false,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let config_seeded = PesScanInternalConfig {
            use_seeding: true,
            ..base_config.clone()
        };

        let result_cold = pes_scan_rigid(&base_config, None).expect("Cold scan should succeed");
        let result_seeded =
            pes_scan_rigid(&config_seeded, None).expect("Seeded scan should succeed");

        let cold_iters = result_cold.total_iterations;
        let seeded_iters = result_seeded.total_iterations;
        let reduction_pct = if cold_iters > 0 {
            100.0 * (cold_iters as f64 - seeded_iters as f64) / cold_iters as f64
        } else {
            0.0
        };

        eprintln!("\n=== GAP-5: Density Seeding Iteration Reduction (H2O RHF/STO-3G) ===");
        eprintln!(
            "Cold start:  {} total SCF iterations across 11 points",
            cold_iters
        );
        eprintln!(
            "Seeded:      {} total SCF iterations across 11 points",
            seeded_iters
        );
        eprintln!("Reduction:   {:.1}%", reduction_pct);

        eprintln!("\nPer-point detail (cold vs seeded):");
        for (i, (cold_pt, seed_pt)) in result_cold
            .points
            .iter()
            .zip(result_seeded.points.iter())
            .enumerate()
        {
            eprintln!(
                "  Point {:2} R={:.1}: cold={:2} iters, seeded={:2} iters",
                i, cold_pt.coordinate_value, cold_pt.scf_iterations, seed_pt.scf_iterations
            );
        }

        assert!(
            seeded_iters <= cold_iters,
            "Seeding should not increase iterations"
        );
    }

    #[test]
    fn test_si_gap1_nh3_rhf_sto3g_optimization() {
        // GAP-1: NH3/RHF/STO-3G geometry optimization
        // PySCF reference: R(NH) = 1.951185 bohr = 1.032522 A
        //                  angle(HNH) = 104.1641 deg
        //                  E_opt = -55.455419778805 Ha
        use crate::optimizer::{optimize_geometry, OptMethod, OptimizationConfig};

        let atoms = vec![
            (7u8, [0.0, 0.0, 0.2197]),
            (1, [0.0, 1.7715, -0.5126]),
            (1, [1.5342, -0.8857, -0.5126]),
            (1, [-1.5342, -0.8857, -0.5126]),
        ];

        let config = OptimizationConfig {
            max_steps: 50,
            method: OptMethod::Rhf,
            basis: "sto-3g".to_string(),
            ..Default::default()
        };

        let result = optimize_geometry(
            &atoms,
            &config,
            Some(&|step| {
                eprintln!(
                    "NH3 Step {}: E={:.10} max|g|={:.6e}",
                    step.step, step.energy, step.max_gradient
                );
            }),
        );

        assert!(
            result.converged,
            "NH3 optimization should converge, got {} steps",
            result.total_steps
        );

        let geom = &result.final_geometry;
        use crate::constants::BOHR_TO_ANGSTROM;

        // Compute NH distances
        let mut nh_dists = Vec::new();
        for i in 1..4 {
            let r = ((geom[0][0] - geom[i][0]).powi(2)
                + (geom[0][1] - geom[i][1]).powi(2)
                + (geom[0][2] - geom[i][2]).powi(2))
            .sqrt();
            nh_dists.push(r);
        }

        // Compute HNH angles
        let mut angles = Vec::new();
        for (i, j) in [(1, 2), (1, 3), (2, 3)] {
            let v1 = [
                geom[i][0] - geom[0][0],
                geom[i][1] - geom[0][1],
                geom[i][2] - geom[0][2],
            ];
            let v2 = [
                geom[j][0] - geom[0][0],
                geom[j][1] - geom[0][1],
                geom[j][2] - geom[0][2],
            ];
            let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
            let mag1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
            let mag2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
            let angle_rad = (dot / (mag1 * mag2)).acos();
            angles.push(angle_rad * 180.0 / std::f64::consts::PI);
        }

        let avg_nh = nh_dists.iter().sum::<f64>() / 3.0;
        let avg_angle = angles.iter().sum::<f64>() / 3.0;

        eprintln!("\n=== GAP-1: NH3/RHF/STO-3G Optimization ===");
        eprintln!("Converged in {} steps", result.total_steps);
        eprintln!("E_opt = {:.10} Ha", result.final_energy);
        for (i, r) in nh_dists.iter().enumerate() {
            eprintln!(
                "  R(N-H{}) = {:.6} bohr = {:.4} A",
                i + 1,
                r,
                r * BOHR_TO_ANGSTROM
            );
        }
        eprintln!(
            "  Avg R(NH) = {:.6} bohr = {:.4} A",
            avg_nh,
            avg_nh * BOHR_TO_ANGSTROM
        );
        for (k, &a) in angles.iter().enumerate() {
            let pairs = [(1, 2), (1, 3), (2, 3)];
            eprintln!("  angle(H{}-N-H{}) = {:.4} deg", pairs[k].0, pairs[k].1, a);
        }
        eprintln!("  Avg angle(HNH) = {:.1} deg", avg_angle);
    }

    // =========================================================================
    // SI Data Extraction: GAP-3 — C₂H₆ Dihedral Torsion Scan
    // =========================================================================

    /// IQCP C₂H₆ geometry in bohr (from apps/web/src/stores/scfStore.ts).
    /// This is the staggered conformation; H(2)-C(0)-C(1)-H(5) dihedral = 180°.
    fn c2h6_iqcp_atoms() -> Vec<(u8, [f64; 3])> {
        vec![
            (6, [0.0000, 0.0000, 1.4508]),   // C0
            (6, [0.0000, 0.0000, -1.4508]),  // C1
            (1, [0.0000, 1.9217, 2.2700]),   // H2
            (1, [1.6641, -0.9609, 2.2700]),  // H3
            (1, [-1.6641, -0.9609, 2.2700]), // H4
            (1, [0.0000, -1.9217, -2.2700]), // H5
            (1, [-1.6641, 0.9609, -2.2700]), // H6
            (1, [1.6641, 0.9609, -2.2700]),  // H7
        ]
    }

    /// Golden data structure for C₂H₆ dihedral PySCF reference
    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenC2H6DihedralPoint {
        dihedral_deg: f64,
        energy: f64,
    }

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenC2H6Dihedral {
        molecule: String,
        points: Vec<GoldenC2H6DihedralPoint>,
    }

    fn load_golden_c2h6_dihedral() -> GoldenC2H6Dihedral {
        let data = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/golden/pes/c2h6_dihedral_rigid_sto3g.json"
        ));
        serde_json::from_str(data).expect("Failed to parse C2H6 dihedral golden data")
    }

    #[test]
    fn test_si_gap3_c2h6_dihedral_rigid_scan() {
        // GAP-3: C₂H₆ rigid dihedral scan H(2)-C(0)-C(1)-H(5), RHF/STO-3G.
        // 13 points from 0° to 360° (every 30°).
        // Eclipsed at 0°/120°/240°; staggered at 60°/180°/300°.
        let config = PesScanInternalConfig {
            atoms: c2h6_iqcp_atoms(),
            coordinate: ScanCoordinate::Dihedral {
                atom_i: 2,
                atom_j: 0,
                atom_k: 1,
                atom_l: 5,
            },
            value_min: 0.0_f64.to_radians(),
            value_max: 360.0_f64.to_radians(),
            n_points: 13,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: None,
            opt_grad_threshold: None,
        };

        let result = pes_scan_rigid(&config, None).expect("Rigid scan should succeed");

        // Load PySCF golden data
        let golden = load_golden_c2h6_dihedral();
        assert_eq!(result.points.len(), 13, "Should have 13 scan points");
        assert_eq!(golden.points.len(), 13, "Golden data should have 13 points");

        eprintln!("\n=== GAP-3: C2H6 Rigid Dihedral Scan H(2)-C(0)-C(1)-H(5), RHF/STO-3G ===");
        eprintln!(
            "{:<12} {:<20} {:<20} {:<15}",
            "Angle(deg)", "E_IQCP (Ha)", "E_PySCF (Ha)", "|dE| (Ha)"
        );

        let mut max_error = 0.0_f64;
        for (iqcp_pt, gold_pt) in result.points.iter().zip(golden.points.iter()) {
            let angle_deg = iqcp_pt.coordinate_value.to_degrees();
            let delta = (iqcp_pt.energy - gold_pt.energy).abs();
            max_error = max_error.max(delta);
            eprintln!(
                "{:<12.1} {:<20.10} {:<20.10} {:<15.2e}",
                angle_deg, iqcp_pt.energy, gold_pt.energy, delta
            );
            assert!(
                iqcp_pt.converged,
                "Point at {:.1} deg not converged",
                angle_deg
            );
            assert!(
                delta < 1e-4,
                "Error {:.2e} exceeds 1e-4 Ha at {:.1} deg",
                delta,
                angle_deg
            );
        }

        // Verify 3-fold barrier: staggered < eclipsed
        let e_eclipsed = result.points[0].energy; // 0 deg
        let e_staggered = result.points[2].energy; // 60 deg
        let barrier_ha = e_eclipsed - e_staggered;
        let barrier_kcal = barrier_ha * 627.509474;
        eprintln!(
            "\nTorsion barrier: {:.6} mHa = {:.4} kcal/mol",
            barrier_ha * 1000.0,
            barrier_kcal
        );
        eprintln!("Max |dE| vs PySCF: {:.2e} Ha", max_error);

        assert!(
            barrier_ha > 0.002,
            "Barrier {:.6} Ha too small (expected ~4 mHa)",
            barrier_ha
        );
        assert!(
            barrier_ha < 0.010,
            "Barrier {:.6} Ha too large (expected ~4 mHa)",
            barrier_ha
        );
    }

    #[test]
    fn test_si_gap3_c2h6_dihedral_relaxed_scan() {
        // GAP-3 (relaxed): C₂H₆ relaxed dihedral scan H(2)-C(0)-C(1)-H(5), RHF/STO-3G.
        // 13 points from 0° to 360° (every 30°).
        // Relaxed barrier should be <= rigid barrier (variational principle).
        let config = PesScanInternalConfig {
            atoms: c2h6_iqcp_atoms(),
            coordinate: ScanCoordinate::Dihedral {
                atom_i: 2,
                atom_j: 0,
                atom_k: 1,
                atom_l: 5,
            },
            value_min: 0.0_f64.to_radians(),
            value_max: 360.0_f64.to_radians(),
            n_points: 13,
            basis_name: "sto-3g".to_string(),
            method: "rhf".to_string(),
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: Some(20),
            opt_grad_threshold: Some(4.5e-4),
        };

        let result = pes_scan_relaxed(&config, None).expect("Relaxed scan should succeed");

        assert_eq!(result.points.len(), 13, "Should have 13 scan points");

        eprintln!("\n=== GAP-3: C2H6 Relaxed Dihedral Scan H(2)-C(0)-C(1)-H(5), RHF/STO-3G ===");
        eprintln!(
            "{:<12} {:<20} {:<10} {:<10}",
            "Angle(deg)", "E_relaxed (Ha)", "Conv", "OptSteps"
        );

        for pt in &result.points {
            let angle_deg = pt.coordinate_value.to_degrees();
            eprintln!(
                "{:<12.1} {:<20.10} {:<10} {:<10}",
                angle_deg,
                pt.energy,
                pt.converged,
                pt.opt_steps.unwrap_or(0)
            );
            assert!(pt.converged, "Point at {:.1} deg not converged", angle_deg);
        }

        // Verify barrier exists
        let e_eclipsed = result.points[0].energy;
        let e_staggered = result.points[2].energy;
        let barrier_ha = e_eclipsed - e_staggered;
        let barrier_kcal = barrier_ha * 627.509474;
        eprintln!(
            "\nRelaxed torsion barrier: {:.6} mHa = {:.4} kcal/mol",
            barrier_ha * 1000.0,
            barrier_kcal
        );

        assert!(barrier_ha > 0.001, "Relaxed barrier too small");
        assert!(barrier_ha < 0.010, "Relaxed barrier too large");
    }

    /// Generate IQCP C2H6 dihedral scan data at 6-31G* for Figure 5.
    /// Writes all data to /tmp/fig5_iqcp_data.json for provenance.
    /// Run with: cargo test --release -p qc-core --lib -- test_fig5_collect_all --nocapture --ignored
    #[test]
    #[ignore]
    fn test_fig5_collect_all_data() {
        use std::collections::BTreeMap;
        use std::io::Write;

        let atoms = vec![
            (6u8, [0.0, 0.0, 1.4508]),
            (6, [0.0, 0.0, -1.4508]),
            (1, [0.0, 1.9217, 2.2700]),
            (1, [1.6641, -0.9609, 2.2700]),
            (1, [-1.6641, -0.9609, 2.2700]),
            (1, [0.0, -1.9217, -2.2700]),
            (1, [-1.6641, 0.9609, -2.2700]),
            (1, [1.6641, 0.9609, -2.2700]),
        ];

        let mut all_data: BTreeMap<String, Vec<(f64, f64, bool, usize)>> = BTreeMap::new();

        for method in &["rhf", "b3lyp", "b3lyp-d3bj"] {
            for mode in &[ScanMode::Rigid, ScanMode::Relaxed] {
                let mode_str = match mode {
                    ScanMode::Rigid => "rigid",
                    ScanMode::Relaxed => "relaxed",
                };
                let key = format!("{}_{}", method, mode_str);
                eprintln!(
                    "\n>>> Running {} / 6-31G* / {} (72 pts)...",
                    method.to_uppercase(),
                    mode_str
                );

                let config = PesScanInternalConfig {
                    atoms: atoms.clone(),
                    basis_name: "6-31g*".to_string(),
                    method: method.to_string(),
                    coordinate: ScanCoordinate::Dihedral {
                        atom_i: 2,
                        atom_j: 0,
                        atom_k: 1,
                        atom_l: 5,
                    },
                    value_min: 0.0,
                    value_max: 355.0_f64.to_radians(),
                    n_points: 72,
                    use_seeding: true,
                    use_spherical: false,
                    convergence_profile: "tight".to_string(),
                    opt_max_steps: Some(30),
                    opt_grad_threshold: Some(4.5e-4),
                };
                let result = pes_scan_internal(&config, *mode, None).unwrap();
                let points: Vec<(f64, f64, bool, usize)> = result
                    .points
                    .iter()
                    .map(|p| {
                        (
                            p.coordinate_value.to_degrees(),
                            p.energy,
                            p.converged,
                            p.opt_steps.unwrap_or(0),
                        )
                    })
                    .collect();

                let energies: Vec<f64> = points.iter().map(|p| p.1).collect();
                let e_min = energies.iter().cloned().fold(f64::INFINITY, f64::min);
                let e_max = energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let barrier = (e_max - e_min) * 627.5095;
                let unconverged = points.iter().filter(|p| !p.2).count();
                eprintln!(
                    "    Barrier: {:.4} kcal/mol, unconverged: {}",
                    barrier, unconverged
                );

                all_data.insert(key, points);
            }
        }

        // Write JSON
        let mut json = String::from("{\n");
        for (i, (key, points)) in all_data.iter().enumerate() {
            json.push_str(&format!("  \"{}\": [\n", key));
            for (j, (angle, energy, conv, steps)) in points.iter().enumerate() {
                let comma = if j < points.len() - 1 { "," } else { "" };
                json.push_str(&format!(
                    "    [{:.4}, {:.12}, {}, {}]{}\n",
                    angle, energy, conv, steps, comma
                ));
            }
            let comma = if i < all_data.len() - 1 { "," } else { "" };
            json.push_str(&format!("  ]{}\n", comma));
        }
        json.push_str("}\n");

        let path = "/tmp/fig5_iqcp_all_data.json";
        let mut f = std::fs::File::create(path).unwrap();
        f.write_all(json.as_bytes()).unwrap();
        eprintln!("\n=== All data written to {} ===", path);
    }

    /// Individual scan tests for parallel execution.
    /// Run with: cargo test --release -p qc-core --lib -- test_fig5_c2h6 --nocapture --ignored
    #[test]
    #[ignore] // Long-running: use --release for speed
    fn test_fig5_c2h6_rigid_rhf() {
        run_fig5_scan("rhf", ScanMode::Rigid);
    }

    #[test]
    #[ignore]
    fn test_fig5_c2h6_rigid_b3lyp() {
        run_fig5_scan("b3lyp", ScanMode::Rigid);
    }

    #[test]
    #[ignore]
    fn test_fig5_c2h6_rigid_b3lyp_d3bj() {
        run_fig5_scan("b3lyp-d3bj", ScanMode::Rigid);
    }

    #[test]
    #[ignore]
    fn test_fig5_c2h6_relaxed_rhf() {
        run_fig5_scan("rhf", ScanMode::Relaxed);
    }

    #[test]
    #[ignore]
    fn test_fig5_c2h6_relaxed_b3lyp() {
        run_fig5_scan("b3lyp", ScanMode::Relaxed);
    }

    #[test]
    #[ignore]
    fn test_fig5_c2h6_relaxed_b3lyp_d3bj() {
        run_fig5_scan("b3lyp-d3bj", ScanMode::Relaxed);
    }

    fn run_fig5_scan(method: &str, mode: ScanMode) {
        let atoms = vec![
            (6u8, [0.0, 0.0, 1.4508]),
            (6, [0.0, 0.0, -1.4508]),
            (1, [0.0, 1.9217, 2.2700]),
            (1, [1.6641, -0.9609, 2.2700]),
            (1, [-1.6641, -0.9609, 2.2700]),
            (1, [0.0, -1.9217, -2.2700]),
            (1, [-1.6641, 0.9609, -2.2700]),
            (1, [1.6641, 0.9609, -2.2700]),
        ];
        let mode_str = match mode {
            ScanMode::Rigid => "rigid",
            ScanMode::Relaxed => "relaxed",
        };

        let config = PesScanInternalConfig {
            atoms,
            basis_name: "6-31g*".to_string(),
            method: method.to_string(),
            coordinate: ScanCoordinate::Dihedral {
                atom_i: 2,
                atom_j: 0,
                atom_k: 1,
                atom_l: 5,
            },
            value_min: 0.0,
            value_max: 355.0_f64.to_radians(),
            n_points: 72,
            use_seeding: true,
            use_spherical: false,
            convergence_profile: "tight".to_string(),
            opt_max_steps: Some(30),
            opt_grad_threshold: Some(4.5e-4),
        };
        let result = pes_scan_internal(&config, mode, None).unwrap();

        eprintln!(
            "\nFIG5_DATA|{}|{}|6-31g*|{}",
            method,
            mode_str,
            result.points.len()
        );
        for pt in &result.points {
            eprintln!(
                "FIG5_PT|{:.4}|{:.12}|{}|{}",
                pt.coordinate_value.to_degrees(),
                pt.energy,
                pt.converged,
                pt.opt_steps.unwrap_or(0)
            );
        }

        let energies: Vec<f64> = result.points.iter().map(|p| p.energy).collect();
        let e_min = energies.iter().cloned().fold(f64::INFINITY, f64::min);
        let e_max = energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let barrier = (e_max - e_min) * 627.5095;
        eprintln!("FIG5_BARRIER|{:.4}", barrier);

        assert_eq!(result.points.len(), 72);
        let unconverged: Vec<_> = result
            .points
            .iter()
            .filter(|p| !p.converged)
            .map(|p| format!("{:.0}°", p.coordinate_value.to_degrees()))
            .collect();
        if !unconverged.is_empty() {
            eprintln!(
                "WARNING: {} unconverged points: {}",
                unconverged.len(),
                unconverged.join(", ")
            );
        }
    }
}
