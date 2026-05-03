//! Superposition of Atomic Densities (SAD) initial guess
//!
//! Constructs an initial density matrix as the sum of free-atom densities,
//! providing a much better starting point for SCF convergence than the
//! core Hamiltonian guess.
//!
//! # Theory
//!
//! The SAD guess constructs the initial density matrix as:
//!
//! ```text
//! D_SAD[mu,nu] = sum_A D_A[mu,nu]
//! ```
//!
//! where D_A is the density matrix of free atom A in the molecular basis,
//! computed via spherically-averaged atomic SCF calculations.
//!
//! # Atomic SCF with Spherical Averaging
//!
//! For each unique element, we run a mini atomic HF calculation following
//! PySCF's `AtomSphAverageRHF` approach (pyscf/scf/atom_hf.py, lines 89-188):
//!
//! 1. Build a single-atom BasisSet with only that atom's shells
//! 2. Compute atomic overlap, core Hamiltonian, and ERIs
//! 3. Run SCF iterations with:
//!    - **Spherical averaging** of the Fock matrix: for each l-subshell,
//!      average F and S blocks over the (2l+1) angular components, solve
//!      the reduced (radial-only) eigenproblem, expand back
//!    - **Fractional occupation** for open-shell atoms: distribute electrons
//!      equally among degenerate orbitals (Aufbau with spherical symmetry)
//!
//! This produces off-diagonal density matrix elements (inter-shell coupling
//! within each atom) that capture proper inner/outer valence shell ratios
//! in split-valence bases -- a significant improvement over simple diagonal
//! occupation for molecules like CH4 with 6-31G*.
//!
//! # References
//!
//! - Van Lenthe et al. J. Comput. Chem. 2006, 27, 926-932.
//! - PySCF: `pyscf/scf/atom_hf.py` `AtomSphAverageRHF` class
//! - PySCF: `pyscf/scf/hf.py` `init_guess_by_atom`

use std::collections::HashMap;

use nalgebra::{DMatrix, SymmetricEigen};

use crate::basis::{BasisSet, ContractedShell};
use crate::integrals::{eri_compressed, hcore_matrix, overlap_matrix};

/// Maximum number of atomic SCF iterations.
const ATOMIC_SCF_MAX_ITER: usize = 40;

/// Convergence threshold for atomic SCF energy (Hartree).
const ATOMIC_SCF_ENERGY_TOL: f64 = 1e-10;

// =============================================================================
// Atomic electron configuration
// =============================================================================

/// Ground-state electron configuration per angular momentum subshell.
///
/// Returns (n_electrons_s, n_electrons_p, n_electrons_d) for a given element.
/// These are the total electron counts per l-subshell, used for fractional
/// occupation in the atomic SCF.
///
/// Reference: PySCF `pyscf/data/elements.py` NRSRHF_CONFIGURATION
fn atomic_electron_config(z: u8) -> [usize; 3] {
    match z {
        1 => [1, 0, 0],
        2 => [2, 0, 0],
        3 => [3, 0, 0],   // Li: 1s^2 2s^1
        4 => [4, 0, 0],   // Be: 1s^2 2s^2
        5 => [4, 1, 0],   // B
        6 => [4, 2, 0],   // C
        7 => [4, 3, 0],   // N
        8 => [4, 4, 0],   // O
        9 => [4, 5, 0],   // F
        10 => [4, 6, 0],  // Ne
        11 => [5, 6, 0],  // Na
        12 => [6, 6, 0],  // Mg
        13 => [6, 7, 0],  // Al
        14 => [6, 8, 0],  // Si
        15 => [6, 9, 0],  // P
        16 => [6, 10, 0], // S
        17 => [6, 11, 0], // Cl
        18 => [6, 12, 0], // Ar
        _ => [0, 0, 0],
    }
}

// =============================================================================
// Mini basis set construction
// =============================================================================

/// Build a single-atom BasisSet for atomic SCF calculation.
///
/// Creates a BasisSet containing only one atom (placed at the origin) and its
/// shells, with all shell atom_idx set to 0. Nuclear repulsion is zero.
///
/// Reference: PySCF atom_hf.py line 44-51 (single-atom molecule construction)
fn build_mini_basis(atom_z: u8, atom_position: [f64; 3], shells: &[&ContractedShell]) -> BasisSet {
    let atom = crate::basis::Atom::new(atom_z, atom_position).unwrap();
    let mini_shells: Vec<ContractedShell> = shells
        .iter()
        .map(|s| ContractedShell::new(s.angular_momentum, s.primitives.clone(), atom_position, 0))
        .collect();
    let n_basis = mini_shells.iter().map(|s| s.n_basis_functions()).sum();

    BasisSet {
        atoms: vec![atom],
        shells: mini_shells,
        name: String::new(),
        n_basis,
        n_electrons: atom_z as usize,
        nuclear_repulsion: 0.0,
    }
}

// =============================================================================
// Spherical-average eigensolver
// =============================================================================

/// Information about shell groups for spherical averaging.
///
/// For each angular momentum l present in the basis, we track which basis
/// function indices belong to it and how many radial shells there are.
#[derive(Debug)]
struct LGroupInfo {
    /// Angular momentum quantum number
    l: u32,
    /// Number of angular components: 1 for S, 3 for P, 6 for D (Cartesian)
    n_comp: usize,
    /// Number of radial shells of this l type
    n_radial: usize,
    /// Basis function indices belonging to this l-group, ordered as:
    /// [shell0_comp0, shell0_comp1, ..., shell1_comp0, shell1_comp1, ...]
    bf_indices: Vec<usize>,
}

/// Analyze the shell structure of a mini basis to extract l-group information.
///
/// This identifies which basis functions belong to each angular momentum type
/// and how they're organized into radial shells, mirroring PySCF's
/// `_angular_momentum_for_each_ao` (atom_hf.py line 208-214).
fn analyze_shell_structure(shells: &[ContractedShell]) -> Vec<LGroupInfo> {
    // Group shells by l value, preserving order
    let mut l_groups: Vec<(u32, Vec<(usize, usize)>)> = Vec::new(); // (l, [(bf_start, n_bf)])

    let mut bf_offset = 0;
    for shell in shells {
        let l = shell.l_value();
        let n_bf = shell.n_basis_functions();

        // Find or create group for this l
        if let Some(group) = l_groups.iter_mut().find(|(gl, _)| *gl == l) {
            group.1.push((bf_offset, n_bf));
        } else {
            l_groups.push((l, vec![(bf_offset, n_bf)]));
        }

        bf_offset += n_bf;
    }

    // Convert to LGroupInfo
    l_groups
        .into_iter()
        .map(|(l, shell_ranges)| {
            let n_comp = match l {
                0 => 1,
                1 => 3,
                2 => 6, // Cartesian d
                _ => (2 * l + 1) as usize,
            };
            let n_radial = shell_ranges.len();

            // Build bf_indices in the order PySCF expects:
            // For each radial shell, list its angular components in order
            let mut bf_indices = Vec::with_capacity(n_radial * n_comp);
            for &(start, _n_bf) in &shell_ranges {
                for c in 0..n_comp {
                    bf_indices.push(start + c);
                }
            }

            LGroupInfo {
                l,
                n_comp,
                n_radial,
                bf_indices,
            }
        })
        .collect()
}

/// Spherically-averaged eigensolver following PySCF's `AtomSphAverageRHF.eig`.
///
/// For each angular momentum l:
/// 1. Extract the F and S subblocks for basis functions of that l
/// 2. Reshape into (n_radial, n_comp, n_radial, n_comp)
/// 3. Average over angular components: F_avg[i,j] = sum_m F[i,m,j,m] / n_comp
/// 4. Solve the reduced eigenproblem: F_avg * c = S_avg * c * e
/// 5. Expand eigenvectors back: each angular component gets the same radial coefficients
///
/// Reference: PySCF atom_hf.py lines 109-140
///
/// Returns (eigenvalues, eigenvectors) in the full AO basis, sorted by energy.
fn spherical_average_eig(
    fock: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
    l_groups: &[LGroupInfo],
    nbf: usize,
) -> (Vec<f64>, DMatrix<f64>) {
    let mut all_energies: Vec<f64> = Vec::with_capacity(nbf);
    let mut mo_coeff = DMatrix::zeros(nbf, nbf);
    let mut col_offset = 0;

    for group in l_groups {
        let n_rad = group.n_radial;
        let n_comp = group.n_comp;

        if n_rad == 0 {
            continue;
        }

        // Extract F and S subblocks and average over angular components.
        // PySCF: f_l = f[idx[:,None],idx].reshape(nsh, degen, nsh, degen)
        //        f_l = einsum('piqi->pq', f_l) / degen
        let mut f_avg = DMatrix::zeros(n_rad, n_rad);
        let mut s_avg = DMatrix::zeros(n_rad, n_rad);

        for i_rad in 0..n_rad {
            for j_rad in 0..n_rad {
                let mut f_sum = 0.0;
                let mut s_sum = 0.0;
                for m in 0..n_comp {
                    let row = group.bf_indices[i_rad * n_comp + m];
                    let col = group.bf_indices[j_rad * n_comp + m];
                    f_sum += fock[(row, col)];
                    s_sum += overlap[(row, col)];
                }
                f_avg[(i_rad, j_rad)] = f_sum / n_comp as f64;
                s_avg[(i_rad, j_rad)] = s_sum / n_comp as f64;
            }
        }

        // Solve generalized eigenproblem F_avg * c = S_avg * c * e
        // Use S^{-1/2} transformation
        let (energies, coeffs) = solve_gen_eig(&f_avg, &s_avg);

        // Expand back: each angular component gets the same radial coefficients
        // PySCF: mo[idx[i::degen],:,i] = c  (atom_hf.py line 137)
        for i_rad in 0..n_rad {
            // Eigenvalue repeated for each angular component
            for _m in 0..n_comp {
                all_energies.push(energies[i_rad]);
            }

            // Eigenvector: set coefficients for each angular component
            for m in 0..n_comp {
                let mo_col = col_offset + i_rad * n_comp + m;
                for j_rad in 0..n_rad {
                    let ao_row = group.bf_indices[j_rad * n_comp + m];
                    mo_coeff[(ao_row, mo_col)] = coeffs[(j_rad, i_rad)];
                }
            }
        }

        col_offset += n_rad * n_comp;
    }

    // Sort by eigenvalue (ascending energy)
    let mut indices: Vec<usize> = (0..all_energies.len()).collect();
    indices.sort_by(|&a, &b| all_energies[a].partial_cmp(&all_energies[b]).unwrap());

    let sorted_energies: Vec<f64> = indices.iter().map(|&i| all_energies[i]).collect();
    let mut sorted_mo = DMatrix::zeros(nbf, nbf);
    for (new_col, &old_col) in indices.iter().enumerate() {
        for row in 0..nbf {
            sorted_mo[(row, new_col)] = mo_coeff[(row, old_col)];
        }
    }

    (sorted_energies, sorted_mo)
}

/// Solve generalized eigenvalue problem F*c = S*c*e via symmetric orthogonalization.
///
/// For small matrices (atomic blocks, typically 1-6), this is efficient.
fn solve_gen_eig(f: &DMatrix<f64>, s: &DMatrix<f64>) -> (Vec<f64>, DMatrix<f64>) {
    let n = f.nrows();
    if n == 0 {
        return (vec![], DMatrix::zeros(0, 0));
    }
    if n == 1 {
        // Trivial case: single function
        let e = f[(0, 0)] / s[(0, 0)];
        let c = 1.0 / s[(0, 0)].sqrt();
        return (vec![e], DMatrix::from_element(1, 1, c));
    }

    // S^{-1/2} via eigendecomposition
    let s_eig = SymmetricEigen::new(s.clone());
    let mut x = DMatrix::zeros(n, n);
    for k in 0..n {
        let val = s_eig.eigenvalues[k];
        if val < 1e-10 {
            continue; // Skip near-zero eigenvalues (linear dependence)
        }
        let inv_sqrt = 1.0 / val.sqrt();
        for i in 0..n {
            for j in 0..n {
                x[(i, j)] += s_eig.eigenvectors[(i, k)] * inv_sqrt * s_eig.eigenvectors[(j, k)];
            }
        }
    }

    // Transform: F' = X^T F X
    let f_prime = x.transpose() * f * &x;

    // Diagonalize F'
    let eig = SymmetricEigen::new(f_prime);

    // Sort eigenvalues
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| eig.eigenvalues[a].partial_cmp(&eig.eigenvalues[b]).unwrap());

    let energies: Vec<f64> = indices.iter().map(|&i| eig.eigenvalues[i]).collect();

    // Back-transform eigenvectors: C = X * C'
    let c_prime_sorted = {
        let mut m = DMatrix::zeros(n, n);
        for (new_col, &old_col) in indices.iter().enumerate() {
            for row in 0..n {
                m[(row, new_col)] = eig.eigenvectors[(row, old_col)];
            }
        }
        m
    };
    let c = &x * c_prime_sorted;

    (energies, c)
}

// =============================================================================
// Fractional occupation
// =============================================================================

/// Build density matrix with fractional occupation for spherically-averaged atom.
///
/// For each l-subshell, occupies orbitals Aufbau-style:
/// - Fully doubly-occupied radial shells get occupation 2.0 per angular component
/// - The last partially-filled shell gets fractional occupation spread equally
///   over all (2l+1) angular components
///
/// Reference: PySCF atom_hf.py lines 142-171 (get_occ method)
///
/// # Arguments
/// * `mo_coeff` - MO coefficients from spherical-average eigensolver
/// * `mo_energies` - Orbital energies (sorted ascending)
/// * `l_groups` - Shell structure info
/// * `config` - Atomic electron configuration [n_s, n_p, n_d]
fn build_density_fractional(
    mo_coeff: &DMatrix<f64>,
    l_groups: &[LGroupInfo],
    config: &[usize; 3],
) -> DMatrix<f64> {
    let nbf = mo_coeff.nrows();
    // Assign occupation by l-group
    // MO columns are ordered by l-group (from spherical_average_eig), then sorted
    // But after sorting by energy, the l-group ordering is mixed.
    // We need to assign occupations per l-subshell.
    //
    // Strategy: Build occupation vector in unsorted MO order (by l-group),
    // since density D = sum_i occ_i * c_i * c_i^T uses MO coefficients directly.
    //
    // Actually, the sorted eigenvectors from spherical_average_eig already have
    // degenerate orbitals grouped together. We need to identify which MO belongs
    // to which l by looking at the eigenvalue degeneracy pattern.

    // Simpler approach: directly compute D from occupation per l-group.
    // For each l-group, we know exactly which MO columns correspond to it
    // (since spherical averaging produces block-diagonal structure).
    // After sorting, we need a different approach.

    // Most robust: use PySCF's approach of assigning occupation by l-subshell
    // independently. The MO coefficients from spherical_average_eig have the
    // property that each MO is purely of one l-character. So we can identify
    // the l-character of each MO by checking which l-group has nonzero coefficients.

    // Actually, for density matrix construction, we can work directly in the
    // l-group framework without worrying about MO ordering:
    //   D = sum_l sum_i^{n_rad_l} occ_li * C_li * C_li^T
    // where C_li has entries only in the bf_indices for l-group l.

    // Let's use the l-group-based approach directly.
    let mut density = DMatrix::zeros(nbf, nbf);

    // For each l, figure out occupation
    for group in l_groups {
        let l = group.l as usize;
        if l >= 3 || l >= config.len() {
            continue;
        }

        let n_elec = config[l]; // Total electrons in this l-subshell
        let n_comp = group.n_comp;
        let n_rad = group.n_radial;

        if n_elec == 0 || n_rad == 0 {
            continue;
        }

        // Max capacity of this l: 2 * n_comp electrons per radial shell
        // PySCF: ndocc = ne // (degen*2), frac = (ne/(degen*2) - ndocc) * 2
        let degen = n_comp; // For Cartesian: S=1, P=3, D=6
        let max_per_radial = 2 * degen; // max electrons per radial shell (accounting for spin)
        let n_doubly_occ = n_elec / max_per_radial;
        let remaining = n_elec - n_doubly_occ * max_per_radial;
        let frac_per_component = if remaining > 0 {
            remaining as f64 / degen as f64
        } else {
            0.0
        };

        // Build occupation array for radial shells of this l
        // Each radial shell has n_comp MOs with the SAME occupation
        let mut radial_occ = vec![0.0f64; n_rad];
        for occ_slot in radial_occ.iter_mut().take(n_rad.min(n_doubly_occ)) {
            *occ_slot = 2.0; // Fully occupied
        }
        if n_doubly_occ < n_rad && frac_per_component > 0.0 {
            radial_occ[n_doubly_occ] = frac_per_component;
        }

        // We need MO coefficients for this l-group.
        // From spherical_average_eig, each radial shell's MO has nonzero entries
        // only in the bf_indices of this l-group, and each angular component
        // shares the same radial coefficients.
        //
        // We can extract these from mo_coeff, but it's easier to just re-solve
        // the reduced eigenproblem for this l-group.
        //
        // Actually, mo_coeff already has the correct structure after
        // spherical_average_eig. We just need to find which columns correspond
        // to this l-group.

        // Find MO columns that have nonzero entries in this l-group's bf_indices
        let mut l_mo_cols: Vec<usize> = Vec::new();
        for col in 0..nbf {
            let mut has_l = false;
            for &bf in &group.bf_indices {
                if mo_coeff[(bf, col)].abs() > 1e-15 {
                    has_l = true;
                    break;
                }
            }
            if has_l {
                l_mo_cols.push(col);
            }
        }

        // l_mo_cols should have n_rad * n_comp entries
        // Group them into radial shells (every n_comp consecutive columns
        // should be degenerate)
        // But after energy sorting, they might be interleaved with other l-groups.
        // We need to identify the radial shell each MO belongs to.

        // Actually, for the density construction, we just need:
        // D += occ * C_col * C_col^T
        // where occ is the same for all angular components of the same radial shell.

        // Group l_mo_cols by their radial shell.
        // All n_comp MOs of the same radial shell have the same eigenvalue.
        // We can use eigenvalue-based grouping, but we don't have eigenvalues here.

        // Simpler: since spherical_average_eig assigns the same radial coefficients
        // to all angular components of one radial shell, we can identify radial
        // shells by looking at the first n_comp entries in sorted order.

        // Actually, the simplest correct approach: just iterate over l_mo_cols
        // in order and assign them to radial shells in groups of n_comp.
        if l_mo_cols.len() != n_rad * n_comp {
            // Fallback: this shouldn't happen, but if it does, skip this l-group
            continue;
        }

        for i_rad in 0..n_rad {
            let occupation = radial_occ[i_rad];
            if occupation.abs() < 1e-15 {
                continue;
            }
            for m in 0..n_comp {
                let col = l_mo_cols[i_rad * n_comp + m];
                // D += occ * |C_col><C_col|
                for mu in 0..nbf {
                    let c_mu = mo_coeff[(mu, col)];
                    if c_mu.abs() < 1e-15 {
                        continue;
                    }
                    for nu in 0..nbf {
                        density[(mu, nu)] += occupation * c_mu * mo_coeff[(nu, col)];
                    }
                }
            }
        }
    }

    density
}

// =============================================================================
// Atomic SCF
// =============================================================================

/// Run a spherically-averaged atomic SCF for a single atom.
///
/// Follows PySCF's `AtomSphAverageRHF` (atom_hf.py lines 89-188):
/// 1. Build overlap, core Hamiltonian, and ERI for the single-atom basis
/// 2. Iteratively solve the Roothaan-Hall equations with spherical averaging
/// 3. Use fractional occupation for open-shell subshells
///
/// Returns the atomic density matrix (nbf_atom x nbf_atom) in the atomic basis.
///
/// For atoms with 0-1 electrons, returns immediately without iterative SCF.
fn atomic_scf(atom_z: u8, atom_position: [f64; 3], shells: &[&ContractedShell]) -> DMatrix<f64> {
    let nbf_atom: usize = shells.iter().map(|s| s.n_basis_functions()).sum();
    if nbf_atom == 0 {
        return DMatrix::zeros(0, 0);
    }

    let n_electrons = atom_z as usize;
    let config = atomic_electron_config(atom_z);

    // Build mini basis and compute one-electron integrals
    let mini_basis = build_mini_basis(atom_z, atom_position, shells);
    let s_flat = overlap_matrix(&mini_basis);
    let h_flat = hcore_matrix(&mini_basis);

    // Convert to nalgebra matrices (row-major from integral routines)
    let s_mat = DMatrix::from_row_slice(nbf_atom, nbf_atom, &s_flat);
    let h_core = DMatrix::from_row_slice(nbf_atom, nbf_atom, &h_flat);

    // Analyze shell structure for spherical averaging
    let l_groups = analyze_shell_structure(&mini_basis.shells);

    // For single-electron atoms (H): just diagonalize H_core, no ERIs needed.
    // The eigenvector of H_core gives the optimal atomic orbital for the
    // free atom, which properly distributes density across inner/outer shells
    // in split-valence bases.
    if n_electrons <= 1 {
        let (_energies, mo_coeff) = spherical_average_eig(&h_core, &s_mat, &l_groups, nbf_atom);
        return build_density_fractional(&mo_coeff, &l_groups, &config);
    }

    // Multi-electron atoms: need ERIs for Fock matrix
    let eri = eri_compressed(&mini_basis);

    // Initial guess: diagonalize H_core with spherical averaging
    let (_energies, mut mo_coeff) = spherical_average_eig(&h_core, &s_mat, &l_groups, nbf_atom);
    let mut density = build_density_fractional(&mo_coeff, &l_groups, &config);

    // Atomic SCF iteration loop
    let mut e_old = f64::MAX;

    for _iter in 0..ATOMIC_SCF_MAX_ITER {
        // Build Fock matrix: F = H_core + G(D)
        let fock = build_atomic_fock(&h_core, &density, &eri, nbf_atom);

        // Compute energy for convergence check
        let e_elec = compute_atomic_energy(&density, &h_core, &fock);

        // Check convergence
        if (e_elec - e_old).abs() < ATOMIC_SCF_ENERGY_TOL {
            break;
        }
        e_old = e_elec;

        // Spherical-average eigensolver
        let (_energies, new_mo) = spherical_average_eig(&fock, &s_mat, &l_groups, nbf_atom);
        mo_coeff = new_mo;

        // Build new density with fractional occupation
        density = build_density_fractional(&mo_coeff, &l_groups, &config);
    }

    density
}

/// Build Fock matrix for atomic SCF: F = H_core + G(D)
///
/// G_mn = sum_{ls} D_{ls} * [(mn|ls) - 0.5*(ml|ns)]
fn build_atomic_fock(
    h_core: &DMatrix<f64>,
    density: &DMatrix<f64>,
    eri: &[f64],
    nbf: usize,
) -> DMatrix<f64> {
    let mut fock = h_core.clone();

    for mu in 0..nbf {
        for nu in 0..=mu {
            let mut g_mn = 0.0;
            for lambda in 0..nbf {
                for sigma in 0..nbf {
                    let d_ls = density[(lambda, sigma)];
                    if d_ls.abs() < 1e-15 {
                        continue;
                    }

                    let j_int = eri_get_local(eri, mu, nu, lambda, sigma);
                    let k_int = eri_get_local(eri, mu, lambda, nu, sigma);
                    g_mn += d_ls * (j_int - 0.5 * k_int);
                }
            }
            fock[(mu, nu)] += g_mn;
            if mu != nu {
                fock[(nu, mu)] += g_mn;
            }
        }
    }

    fock
}

/// Compute electronic energy: E = 0.5 * Tr(D * (H + F))
fn compute_atomic_energy(
    density: &DMatrix<f64>,
    h_core: &DMatrix<f64>,
    fock: &DMatrix<f64>,
) -> f64 {
    let hf = h_core + fock;
    let mut energy = 0.0;
    let n = density.nrows();
    for i in 0..n {
        for j in 0..n {
            energy += density[(i, j)] * hf[(j, i)];
        }
    }
    0.5 * energy
}

/// Get ERI value from compressed storage (local version for atomic SCF)
#[inline]
fn eri_get_local(eri: &[f64], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let p = if i >= j {
        i * (i + 1) / 2 + j
    } else {
        j * (j + 1) / 2 + i
    };
    let q = if k >= l {
        k * (k + 1) / 2 + l
    } else {
        l * (l + 1) / 2 + k
    };
    let (p, q) = if p >= q { (p, q) } else { (q, p) };
    eri[p * (p + 1) / 2 + q]
}

// =============================================================================
// Public API
// =============================================================================

/// Build a SAD (Superposition of Atomic Densities) initial density matrix.
///
/// For each unique element in the molecule, runs a spherically-averaged
/// atomic SCF calculation (following PySCF's `AtomSphAverageRHF`) to produce
/// an atomic density matrix. These atomic densities are then block-diagonally
/// assembled into the molecular density matrix.
///
/// The atomic SCF produces off-diagonal elements within each atom's block,
/// capturing inter-shell coupling (e.g., 1s-2s interaction, inner/outer valence
/// mixing in split-valence bases). This is a significant improvement over
/// simple diagonal occupation guesses, especially for split-valence bases
/// like 6-31G* where the inner/outer valence ratio matters.
///
/// # Arguments
///
/// * `basis` - The molecular basis set containing atom and shell information
///
/// # Returns
///
/// A flat column-major density matrix of size `nbf * nbf`, suitable for
/// passing to `rhf_scf_with_guess` or use as the initial density in `ks_scf`.
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{Atom, BasisSet};
/// use qc_core::scf::sad::build_sad_density;
///
/// let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
/// let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
/// let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
///
/// let d_sad = build_sad_density(&basis);
/// assert_eq!(d_sad.len(), 4); // 2x2 matrix
/// ```
pub fn build_sad_density(basis: &BasisSet) -> Vec<f64> {
    let nbf = basis.n_basis;
    let mut density = DMatrix::<f64>::zeros(nbf, nbf);

    // Cache atomic density matrices by element AND basis name combination.
    // For most molecules, the same element appears multiple times (e.g., H in H2O)
    // and all instances at the same geometry produce identical atomic densities
    // (since the atomic integrals depend on the basis exponents/coefficients, not position,
    // EXCEPT for nuclear attraction integrals which depend on the atom position).
    //
    // Actually, for a single atom in vacuum, the nuclear attraction integral
    // V_mu,nu = -Z * <mu|1/r_A|nu> only involves the atom's own nucleus,
    // and since mu and nu are centered on the same atom, and V depends on
    // |r - R_A| where R_A is the atom center = center of mu and nu,
    // the result is INDEPENDENT of the atom's absolute position.
    //
    // So we CAN cache by (element, basis_set_name).
    // But to be safe and simple, we cache by atomic_number only,
    // since all atoms of the same element use the same basis shells.
    let mut atom_density_cache: HashMap<u8, DMatrix<f64>> = HashMap::new();

    for atom_idx in 0..basis.atoms.len() {
        let atom = &basis.atoms[atom_idx];
        let z = atom.atomic_number;

        // Get cached or compute atomic density
        let atom_density = atom_density_cache.entry(z).or_insert_with(|| {
            let atom_shells: Vec<&ContractedShell> = basis
                .shells
                .iter()
                .filter(|s| s.atom_idx == atom_idx)
                .collect();
            atomic_scf(z, atom.position, &atom_shells)
        });

        // Find basis function offset for this atom
        let bf_offset: usize = basis
            .shells
            .iter()
            .filter(|s| s.atom_idx < atom_idx)
            .map(|s| s.n_basis_functions())
            .sum();

        let n_atom_bf = atom_density.nrows();

        // Copy atomic density block into molecular density
        for mu in 0..n_atom_bf {
            for nu in 0..n_atom_bf {
                density[(bf_offset + mu, bf_offset + nu)] = atom_density[(mu, nu)];
            }
        }
    }

    // Return as column-major flat array (nalgebra default storage order)
    density.as_slice().to_vec()
}

// =============================================================================
// Legacy API (kept for reference/testing)
// =============================================================================

/// Atomic subshell occupation for simple diagonal SAD guess.
///
/// This was the original diagonal-only approach. Kept for testing purposes.
#[cfg(test)]
#[derive(Debug, Clone, Copy)]
struct SubshellOccupation {
    l: u32,
    electrons: f64,
}

/// Get the ground-state subshell occupations for a free atom.
///
/// Returns the occupations in shell order (matching the basis set ordering):
/// 1s, 2s, 2p, 3s, 3p, 3d, ...
#[cfg(test)]
fn atomic_subshell_occupations(z: u8) -> Vec<SubshellOccupation> {
    match z {
        // Period 1
        1 => vec![SubshellOccupation {
            l: 0,
            electrons: 1.0,
        }],
        2 => vec![SubshellOccupation {
            l: 0,
            electrons: 2.0,
        }],
        // Period 2
        3 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 1.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 0.0,
            },
        ],
        4 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 0.0,
            },
        ],
        5 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 1.0,
            },
        ],
        6 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 2.0,
            },
        ],
        7 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 3.0,
            },
        ],
        8 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 4.0,
            },
        ],
        9 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 5.0,
            },
        ],
        10 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 6.0,
            },
        ],
        // Period 3
        11 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 6.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 1.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 0.0,
            },
        ],
        12 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 6.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 0.0,
            },
        ],
        13 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 6.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 1.0,
            },
        ],
        14 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 6.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 2.0,
            },
        ],
        15 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 6.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 3.0,
            },
        ],
        16 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 6.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 4.0,
            },
        ],
        17 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 6.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 5.0,
            },
        ],
        18 => vec![
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 6.0,
            },
            SubshellOccupation {
                l: 0,
                electrons: 2.0,
            },
            SubshellOccupation {
                l: 1,
                electrons: 6.0,
            },
        ],
        _ => vec![],
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::Atom;
    use approx::assert_abs_diff_eq;

    // -------------------------------------------------------------------------
    // Atomic configuration tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_atomic_config_electron_count() {
        // Total electrons from config must match Z
        for z in 1..=18u8 {
            let config = atomic_electron_config(z);
            let total: usize = config.iter().sum();
            assert_eq!(
                total, z as usize,
                "Z={}: config {:?} sums to {}",
                z, config, total
            );
        }
    }

    #[test]
    fn test_hydrogen_config() {
        assert_eq!(atomic_electron_config(1), [1, 0, 0]);
    }

    #[test]
    fn test_carbon_config() {
        assert_eq!(atomic_electron_config(6), [4, 2, 0]);
    }

    #[test]
    fn test_oxygen_config() {
        assert_eq!(atomic_electron_config(8), [4, 4, 0]);
    }

    // -------------------------------------------------------------------------
    // Shell structure analysis tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_shell_structure_sto3g_oxygen() {
        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let basis = BasisSet::build(vec![o], "sto-3g").unwrap();
        let l_groups = analyze_shell_structure(&basis.shells);

        // O STO-3G: 1s, 2s (2 S-shells), 2p (1 P-shell)
        assert_eq!(l_groups.len(), 2); // S and P groups
        assert_eq!(l_groups[0].l, 0);
        assert_eq!(l_groups[0].n_radial, 2);
        assert_eq!(l_groups[0].n_comp, 1);
        assert_eq!(l_groups[1].l, 1);
        assert_eq!(l_groups[1].n_radial, 1);
        assert_eq!(l_groups[1].n_comp, 3);
    }

    #[test]
    fn test_shell_structure_631gs_oxygen() {
        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let basis = BasisSet::build(vec![o], "6-31g*").unwrap();
        let l_groups = analyze_shell_structure(&basis.shells);

        // O 6-31G*: 3 S-shells, 2 P-shells, 1 D-shell
        assert_eq!(l_groups.len(), 3); // S, P, D
        assert_eq!(l_groups[0].l, 0);
        assert_eq!(l_groups[0].n_radial, 3);
        assert_eq!(l_groups[1].l, 1);
        assert_eq!(l_groups[1].n_radial, 2);
        assert_eq!(l_groups[2].l, 2);
        assert_eq!(l_groups[2].n_radial, 1);
        assert_eq!(l_groups[2].n_comp, 6); // Cartesian d
    }

    // -------------------------------------------------------------------------
    // SAD density matrix tests -- basic properties
    // -------------------------------------------------------------------------

    #[test]
    fn test_sad_h2_sto3g() {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let d = build_sad_density(&basis);
        assert_eq!(d.len(), 4); // 2x2

        // Column-major: D[i,j] = d[i + j*nbf]
        let nbf = 2;
        let get = |i: usize, j: usize| d[i + j * nbf];

        // H1: 1s occupation = 1.0 (single s function, atomic SCF trivial)
        assert_abs_diff_eq!(get(0, 0), 1.0, epsilon = 1e-8);
        // H2: 1s occupation = 1.0
        assert_abs_diff_eq!(get(1, 1), 1.0, epsilon = 1e-8);

        // Off-diagonal should be zero (different atoms)
        assert_abs_diff_eq!(get(0, 1), 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(get(1, 0), 0.0, epsilon = 1e-10);

        // Total trace should equal total atomic electrons (2)
        let trace = get(0, 0) + get(1, 1);
        assert_abs_diff_eq!(trace, 2.0, epsilon = 1e-8);
    }

    #[test]
    fn test_sad_electron_count_via_overlap() {
        // The correct electron count from the density matrix is Tr(D*S),
        // not Tr(D). In a non-orthogonal basis, only Tr(D*S) = N_electrons.
        let systems: Vec<(&str, Vec<Atom>, f64)> = vec![
            (
                "H2",
                vec![
                    Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
                    Atom::new(1, [0.0, 0.0, 1.4]).unwrap(),
                ],
                2.0,
            ),
            (
                "H2O",
                vec![
                    Atom::new(8, [0.0, 0.0, 0.2216]).unwrap(),
                    Atom::new(1, [0.0, 1.4309, -0.8864]).unwrap(),
                    Atom::new(1, [0.0, -1.4309, -0.8864]).unwrap(),
                ],
                10.0,
            ),
            (
                "LiH",
                vec![
                    Atom::new(3, [0.0, 0.0, 0.0]).unwrap(),
                    Atom::new(1, [0.0, 0.0, 3.0139]).unwrap(),
                ],
                4.0,
            ),
            (
                "CH4",
                vec![
                    Atom::new(6, [0.0, 0.0, 0.0]).unwrap(),
                    Atom::new(1, [1.1851, 1.1851, 1.1851]).unwrap(),
                    Atom::new(1, [-1.1851, -1.1851, 1.1851]).unwrap(),
                    Atom::new(1, [-1.1851, 1.1851, -1.1851]).unwrap(),
                    Atom::new(1, [1.1851, -1.1851, -1.1851]).unwrap(),
                ],
                10.0,
            ),
        ];

        for (name, atoms, expected_n_elec) in systems {
            for basis_name in &["sto-3g", "6-31g*"] {
                let basis = BasisSet::build(atoms.clone(), basis_name).unwrap();
                let d_flat = build_sad_density(&basis);
                let nbf = basis.n_basis;
                let s_flat = overlap_matrix(&basis);

                // Tr(D*S) = sum_ij D_ij * S_ji = sum_ij D_ij * S_ij (S symmetric)
                // D is column-major, S is row-major
                let mut tr_ds = 0.0;
                for i in 0..nbf {
                    for j in 0..nbf {
                        let d_ij = d_flat[i + j * nbf]; // column-major
                        let s_ij = s_flat[i * nbf + j]; // row-major
                        tr_ds += d_ij * s_ij;
                    }
                }

                assert_abs_diff_eq!(
                    tr_ds,
                    expected_n_elec,
                    epsilon = 0.1, // SAD is approximate; Tr(D*S) ~ N within ~0.1
                );
            }
        }
    }

    #[test]
    fn test_sad_density_is_symmetric() {
        // The SAD density matrix should be symmetric: D[i,j] = D[j,i]
        let o = Atom::new(8, [0.0, 0.0, 0.2216]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.4309, -0.8864]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.4309, -0.8864]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();

        let d = build_sad_density(&basis);
        let nbf = basis.n_basis;

        for i in 0..nbf {
            for j in 0..nbf {
                let dij = d[i + j * nbf];
                let dji = d[j + i * nbf];
                assert_abs_diff_eq!(
                    dij,
                    dji,
                    epsilon = 1e-12,
                    // "D[{},{}] = {} != D[{},{}] = {}", i, j, dij, j, i, dji
                );
            }
        }
    }

    #[test]
    fn test_sad_inter_atom_blocks_zero() {
        // Off-diagonal blocks between different atoms should be zero
        // (SAD is block-diagonal in atom blocks)
        let o = Atom::new(8, [0.0, 0.0, 0.2216]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.4309, -0.8864]).unwrap();
        let basis = BasisSet::build(vec![o, h1], "sto-3g").unwrap();

        let d = build_sad_density(&basis);
        let nbf = basis.n_basis; // O: 5 bf, H: 1 bf = 6 total

        // Oxygen block: 0..5, Hydrogen: 5
        for i in 0..5 {
            let val = d[i + 5 * nbf]; // D[i, 5]
            assert_abs_diff_eq!(
                val,
                0.0,
                epsilon = 1e-15,
                // "D[{}, 5] = {} should be zero", i, val
            );
            let val = d[5 + i * nbf]; // D[5, i]
            assert_abs_diff_eq!(val, 0.0, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_sad_has_off_diagonal_within_atom() {
        // For multi-shell atoms with split-valence basis, the atomic SCF should
        // produce nonzero off-diagonal elements within the atom's block.
        // This is the key improvement over the simple diagonal approach.
        let c = Atom::new(6, [0.0, 0.0, 0.0]).unwrap();
        let basis = BasisSet::build(vec![c], "6-31g*").unwrap();

        let d = build_sad_density(&basis);
        let nbf = basis.n_basis;

        // Carbon 6-31G* has 3 S-shells, 2 P-shells, 1 D-shell = 15 BFs
        // The s-block should have off-diagonal elements (1s-2s coupling)
        // D[0,1] and D[1,0] should be nonzero
        let d01 = d[0 + 1 * nbf]; // D[0,1]
                                  // With atomic SCF, at least the s-s off-diagonal should be nonzero
                                  // (1s-2s mixing from the Fock matrix diagonalization)
        assert!(
            d01.abs() > 1e-6,
            "D[0,1] = {} should be nonzero for C 6-31G*",
            d01
        );
    }

    // -------------------------------------------------------------------------
    // Helper for SCF convergence tests
    // -------------------------------------------------------------------------

    fn build_test_system(basis: &BasisSet) -> crate::scf::PresetSystem {
        crate::scf::PresetSystem {
            system_id: "sad_test".to_string(),
            label: "SAD test system".to_string(),
            nbf: basis.n_basis,
            nelec: basis.n_electrons,
            e_nuc: basis.nuclear_repulsion,
            s_matrix: overlap_matrix(basis),
            h_core: hcore_matrix(basis),
            eri_compressed: eri_compressed(basis),
        }
    }

    // -------------------------------------------------------------------------
    // SCF convergence tests -- same final energy regardless of guess
    // -------------------------------------------------------------------------

    #[test]
    fn test_sad_guess_matches_core_hamiltonian_energy() {
        use crate::scf::{rhf_scf, rhf_scf_with_guess, ConvergenceProfile, ScfConfig};

        let o = Atom::new(8, [0.0, 0.0, 0.2216]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.4309, -0.8864]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.4309, -0.8864]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();
        let system = build_test_system(&basis);

        let config = ScfConfig {
            profile: ConvergenceProfile::Tight,
            max_iterations: 100,
            use_diis: true,
            ..ScfConfig::default()
        };

        let result_core = rhf_scf(&system, &config).expect("Core Hamiltonian SCF failed");

        let d_sad = build_sad_density(&basis);
        let result_sad =
            rhf_scf_with_guess(&system, &config, Some(&d_sad)).expect("SAD SCF failed");

        // Converged energies must agree
        assert_abs_diff_eq!(
            result_core.energy_total,
            result_sad.energy_total,
            epsilon = 1e-10
        );

        // SAD should converge in fewer or equal iterations
        assert!(
            result_sad.iterations <= result_core.iterations,
            "SAD ({} iters) should use <= core Hamiltonian ({} iters)",
            result_sad.iterations,
            result_core.iterations
        );
    }

    #[test]
    fn test_sad_reduces_iterations_h2() {
        use crate::scf::{rhf_scf, rhf_scf_with_guess, ConvergenceProfile, ScfConfig};

        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
        let system = build_test_system(&basis);

        let config = ScfConfig {
            profile: ConvergenceProfile::Tight,
            max_iterations: 100,
            use_diis: true,
            ..ScfConfig::default()
        };

        let result_core = rhf_scf(&system, &config).unwrap();
        let d_sad = build_sad_density(&basis);
        let result_sad = rhf_scf_with_guess(&system, &config, Some(&d_sad)).unwrap();

        assert_abs_diff_eq!(
            result_core.energy_total,
            result_sad.energy_total,
            epsilon = 1e-10
        );

        assert!(result_sad.iterations <= result_core.iterations);
    }

    #[test]
    fn test_sad_iteration_reduction_summary() {
        use crate::scf::{rhf_scf, rhf_scf_with_guess, ConvergenceProfile, ScfConfig};

        let config = ScfConfig {
            profile: ConvergenceProfile::Tight,
            max_iterations: 100,
            use_diis: true,
            ..ScfConfig::default()
        };

        let systems: Vec<(&str, Vec<Atom>)> = vec![
            (
                "H2",
                vec![
                    Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
                    Atom::new(1, [0.0, 0.0, 1.4]).unwrap(),
                ],
            ),
            (
                "LiH",
                vec![
                    Atom::new(3, [0.0, 0.0, 0.0]).unwrap(),
                    Atom::new(1, [0.0, 0.0, 3.0139]).unwrap(),
                ],
            ),
            (
                "H2O",
                vec![
                    Atom::new(8, [0.0, 0.0, 0.2216]).unwrap(),
                    Atom::new(1, [0.0, 1.4309, -0.8864]).unwrap(),
                    Atom::new(1, [0.0, -1.4309, -0.8864]).unwrap(),
                ],
            ),
        ];

        for (name, atoms) in &systems {
            let basis = BasisSet::build(atoms.clone(), "sto-3g").unwrap();
            let system = build_test_system(&basis);

            let result_core = rhf_scf(&system, &config).unwrap();
            let d_sad = build_sad_density(&basis);
            let result_sad = rhf_scf_with_guess(&system, &config, Some(&d_sad)).unwrap();

            assert_abs_diff_eq!(
                result_core.energy_total,
                result_sad.energy_total,
                epsilon = 1e-10
            );

            assert!(
                result_sad.iterations <= result_core.iterations,
                "{}: SAD ({}) > core Hamiltonian ({})",
                name,
                result_sad.iterations,
                result_core.iterations
            );
        }
    }

    // -------------------------------------------------------------------------
    // Legacy API tests (atomic_subshell_occupations)
    // -------------------------------------------------------------------------

    #[test]
    fn test_hydrogen_occupations() {
        let occ = atomic_subshell_occupations(1);
        assert_eq!(occ.len(), 1);
        assert_eq!(occ[0].l, 0);
        assert_abs_diff_eq!(occ[0].electrons, 1.0);
    }

    #[test]
    fn test_helium_occupations() {
        let occ = atomic_subshell_occupations(2);
        assert_eq!(occ.len(), 1);
        assert_eq!(occ[0].l, 0);
        assert_abs_diff_eq!(occ[0].electrons, 2.0);
    }

    #[test]
    fn test_oxygen_occupations() {
        let occ = atomic_subshell_occupations(8);
        assert_eq!(occ.len(), 3);
        assert_eq!(occ[0].l, 0);
        assert_abs_diff_eq!(occ[0].electrons, 2.0);
        assert_eq!(occ[1].l, 0);
        assert_abs_diff_eq!(occ[1].electrons, 2.0);
        assert_eq!(occ[2].l, 1);
        assert_abs_diff_eq!(occ[2].electrons, 4.0);
    }

    #[test]
    fn test_chlorine_occupations() {
        let occ = atomic_subshell_occupations(17);
        assert_eq!(occ.len(), 5);
        let total: f64 = occ.iter().map(|o| o.electrons).sum();
        assert_abs_diff_eq!(total, 17.0);
    }

    #[test]
    fn test_total_electrons_match_z() {
        for z in 1..=18u8 {
            let occ = atomic_subshell_occupations(z);
            let total: f64 = occ.iter().map(|o| o.electrons).sum();
            assert_abs_diff_eq!(total, z as f64, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_sad_matches_pyscf_h_631gs() {
        // Verify our H atomic density in 6-31G* matches PySCF's AtomSphAverageRHF
        // PySCF reference (from atom_hf.get_atm_nrhf):
        //   D[0,0] = 0.18269670, D[0,1] = 0.28443329, D[1,1] = 0.44282296
        //   Tr(D*S) = 1.0
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "6-31g*").unwrap();
        let nbf = basis.n_basis;
        let d = build_sad_density(&basis);

        // H1 block (first 2x2)
        let d00 = d[0];
        let d01 = d[0 + 1 * nbf];
        let d11 = d[1 + 1 * nbf];

        assert_abs_diff_eq!(d00, 0.18269670, epsilon = 1e-6);
        assert_abs_diff_eq!(d01, 0.28443329, epsilon = 1e-6);
        assert_abs_diff_eq!(d11, 0.44282296, epsilon = 1e-6);
    }
}
