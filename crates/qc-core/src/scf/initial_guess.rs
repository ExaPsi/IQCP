//! MINAO (Minimal Basis Atomic Orbitals) initial guess for SCF calculations.
//!
//! Projects a minimal-basis (STO-3G) density matrix onto the target basis set,
//! providing a much better initial guess than the core Hamiltonian for most
//! molecules, and a more robust alternative to SAD for certain systems
//! (e.g., CH4/6-31G* where SAD produces poor starting densities).
//!
//! # Algorithm (following PySCF `hf.py` `init_guess_by_minao`)
//!
//! 1. Build STO-3G basis for the same molecular geometry
//! 2. Compute cross-overlap `S_cross` between target basis and STO-3G basis
//! 3. Compute target self-overlap `S_target`
//! 4. Solve `S_target @ C_proj = S_cross` for `C_proj` (projection coefficients)
//! 5. Assign atomic occupations to STO-3G basis functions
//! 6. Build density: `D = C_proj @ diag(occ) @ C_proj^T`
//!
//! # Key insight
//!
//! The STO-3G basis functions are atom-centered and naturally ordered by energy
//! (1s, 2s, 2p, 3s, 3p), so we can assign occupations directly without
//! diagonalizing any Hamiltonian. The projection onto the larger basis captures
//! the spatial character of the minimal-basis orbitals in the target basis.
//!
//! # References
//!
//! - PySCF: `scf/hf.py` `init_guess_by_minao` (lines 347-473)
//! - PySCF: `scf/addons.py` `project_mo_nr2nr` (lines 580-603)

use nalgebra::DMatrix;

use crate::basis::BasisSet;
use crate::integrals::shell_overlap;

// =============================================================================
// Atomic occupation for MINAO guess
// =============================================================================

/// Electron configuration per angular momentum subshell [s, p, d, f].
///
/// These are total electron counts per l-subshell for neutral atoms,
/// matching PySCF's `NRSRHF_CONFIGURATION` from `pyscf/data/elements.py`.
///
/// The values represent the number of electrons in all shells of that l type:
/// - s: includes 1s, 2s, 3s electrons
/// - p: includes 2p, 3p electrons
/// - d: includes 3d electrons (not relevant for STO-3G/H-Ar)
fn nrsrhf_configuration(z: u8) -> [usize; 4] {
    match z {
        0 => [0, 0, 0, 0],   // Ghost
        1 => [1, 0, 0, 0],   // H
        2 => [2, 0, 0, 0],   // He
        3 => [3, 0, 0, 0],   // Li
        4 => [4, 0, 0, 0],   // Be
        5 => [4, 1, 0, 0],   // B
        6 => [4, 2, 0, 0],   // C
        7 => [4, 3, 0, 0],   // N
        8 => [4, 4, 0, 0],   // O
        9 => [4, 5, 0, 0],   // F
        10 => [4, 6, 0, 0],  // Ne
        11 => [5, 6, 0, 0],  // Na
        12 => [6, 6, 0, 0],  // Mg
        13 => [6, 7, 0, 0],  // Al
        14 => [6, 8, 0, 0],  // Si
        15 => [6, 9, 0, 0],  // P
        16 => [6, 10, 0, 0], // S
        17 => [6, 11, 0, 0], // Cl
        18 => [6, 12, 0, 0], // Ar
        _ => [0, 0, 0, 0],   // Unsupported
    }
}

/// Compute fractional occupation for shells of angular momentum `l` for element `z`.
///
/// Returns `(ndocc, frac)` where:
/// - `ndocc`: number of fully doubly-occupied shells of this `l`
/// - `frac`: fractional occupation per orbital in the partially-filled shell (0.0..2.0)
///
/// Following PySCF `atom_hf.frac_occ` (atom_hf.py lines 196-205).
///
/// Example: Carbon (Z=6), l=1 (p):
///   ne=2 electrons in p-orbitals, nd=6 (3 orbitals * 2 electrons max)
///   ndocc = 2 // 6 = 0
///   frac = (2.0/6.0 - 0) * 2 = 0.6667
///   So each p orbital gets 0.6667 electrons (fractional)
fn frac_occ(z: u8, l: u32) -> (usize, f64) {
    let config = nrsrhf_configuration(z);
    if (l as usize) < 4 && config[l as usize] > 0 {
        let ne = config[l as usize];
        let nd = (l as usize * 2 + 1) * 2; // max electrons per l-subshell
        let ndocc = ne / nd;
        let frac = (ne as f64 / nd as f64 - ndocc as f64) * 2.0;
        (ndocc, frac)
    } else {
        (0, 0.0)
    }
}

/// Compute occupation vector for an atom's STO-3G basis functions.
///
/// For each angular momentum shell in STO-3G:
/// - Doubly-occupied shells get 2.0 per orbital
/// - Partially-occupied shell gets fractional occupation per orbital
/// - Empty shells get 0.0
///
/// The occupation is distributed equally across the (2l+1) degenerate components
/// (spherically-averaged occupation), matching PySCF's MINAO approach.
///
/// Returns: vector of occupations, one per Cartesian basis function in STO-3G
/// (Cartesian because STO-3G has no d-functions, so Cartesian == spherical).
fn atom_occupation(z: u8, shells: &[&crate::basis::ContractedShell]) -> Vec<f64> {
    let mut occ = Vec::new();

    // Group shells by angular momentum to determine shell ordering
    // For STO-3G: H has [1s], Li-Ne have [1s, 2s, 2p], Na-Ar have [1s, 2s, 2p, 3s, 3p]
    for l in 0..=2u32 {
        let (ndocc, frac) = frac_occ(z, l);
        let degen = (2 * l + 1) as usize; // 1 for s, 3 for p, 5 for d

        // Count how many shells of this l the atom has in the basis
        let n_shells_l: usize = shells.iter().filter(|s| s.l_value() == l).count();

        if n_shells_l == 0 {
            continue;
        }

        // Assign occupation: first ndocc shells are doubly occupied,
        // next shell (if any) gets fractional occupation
        for shell_idx in 0..n_shells_l {
            if shell_idx < ndocc {
                // Fully doubly-occupied shell
                occ.extend(std::iter::repeat_n(2.0, degen));
            } else if shell_idx == ndocc && frac > 1e-15 {
                // Partially occupied shell
                occ.extend(std::iter::repeat_n(frac, degen));
            } else {
                // Empty shell
                occ.extend(std::iter::repeat_n(0.0, degen));
            }
        }
    }

    occ
}

// =============================================================================
// Cross-overlap matrix computation
// =============================================================================

/// Compute the rectangular overlap matrix between two different basis sets.
///
/// Returns an `n_a x n_b` matrix (flattened row-major) where:
/// - `n_a` is the number of Cartesian basis functions in `basis_a`
/// - `n_b` is the number of Cartesian basis functions in `basis_b`
///
/// `S_cross[mu, nu] = <chi_mu^a | chi_nu^b>` where `chi_mu^a` is from `basis_a`
/// and `chi_nu^b` is from `basis_b`.
///
/// This reuses the existing `shell_overlap` function which already computes
/// overlap between any two shells regardless of their basis set origin.
///
/// Reference: PySCF `gto.mole.intor_cross('int1e_ovlp', mol_a, mol_b)`
pub fn cross_overlap_matrix(basis_a: &BasisSet, basis_b: &BasisSet) -> Vec<f64> {
    let n_a = basis_a.n_basis;
    let n_b = basis_b.n_basis;
    let mut s_cross = vec![0.0; n_a * n_b];

    let mut mu = 0; // basis function index in basis_a
    for shell_a in &basis_a.shells {
        let n_a_shell = shell_a.n_basis_functions();

        let mut nu = 0; // basis function index in basis_b
        for shell_b in &basis_b.shells {
            let n_b_shell = shell_b.n_basis_functions();

            // Compute shell block of overlap integrals
            let block = shell_overlap(shell_a, shell_b);

            // Copy into the rectangular matrix
            for ia in 0..n_a_shell {
                for ib in 0..n_b_shell {
                    s_cross[(mu + ia) * n_b + (nu + ib)] = block[ia * n_b_shell + ib];
                }
            }

            nu += n_b_shell;
        }
        mu += n_a_shell;
    }

    s_cross
}

/// Compute the cross-overlap matrix with spherical transform on the "a" side.
///
/// When the target basis (`basis_a`) has d-functions, the cross-overlap is first
/// computed in the Cartesian basis and then the rows (target side) are transformed
/// to spherical harmonics. The STO-3G minimal basis (columns) has only s and p
/// shells, so no column transformation is needed.
///
/// Returns an `n_a_sph x n_b` matrix (flattened row-major).
pub fn cross_overlap_matrix_spherical(basis_a: &BasisSet, basis_b: &BasisSet) -> Vec<f64> {
    use crate::integrals::spherical::CART2SPH_D;

    let n_a_cart = basis_a.n_basis;
    let n_a_sph = basis_a.n_basis_spherical();
    let n_b = basis_b.n_basis; // STO-3G: always Cartesian = spherical

    // If no difference, just return the Cartesian cross-overlap
    if n_a_cart == n_a_sph {
        return cross_overlap_matrix(basis_a, basis_b);
    }

    // Compute Cartesian cross-overlap: n_a_cart x n_b
    let cart_cross = cross_overlap_matrix(basis_a, basis_b);

    // Transform rows from Cartesian to spherical for basis_a
    // This applies the cart->sph transform to each shell block of rows
    let mut sph_cross = vec![0.0; n_a_sph * n_b];

    let mut mu_cart = 0;
    let mut mu_sph = 0;
    for shell_a in &basis_a.shells {
        let n_cart = shell_a.n_basis_functions();
        let n_sph = shell_a.n_basis_functions_spherical();

        if n_cart == n_sph {
            // S or P shell: copy directly (no transform needed)
            for ic in 0..n_cart {
                for j in 0..n_b {
                    sph_cross[(mu_sph + ic) * n_b + j] = cart_cross[(mu_cart + ic) * n_b + j];
                }
            }
        } else {
            // D shell: transform 6 Cartesian -> 5 spherical
            // CART2SPH_D is indexed as [cart_idx][sph_idx] (shape 6x5)
            // S_sph[s, j] = sum_c CART2SPH_D[c][s] * S_cart[c, j]
            for s in 0..n_sph {
                for j in 0..n_b {
                    let mut val = 0.0;
                    for c in 0..n_cart {
                        val += CART2SPH_D[c][s] * cart_cross[(mu_cart + c) * n_b + j];
                    }
                    sph_cross[(mu_sph + s) * n_b + j] = val;
                }
            }
        }

        mu_cart += n_cart;
        mu_sph += n_sph;
    }

    sph_cross
}

// =============================================================================
// MINAO initial guess
// =============================================================================

/// Compute MINAO initial guess density matrix.
///
/// Projects a minimal basis (STO-3G) density onto the target basis,
/// providing a better initial guess than core Hamiltonian or SAD for most molecules.
///
/// # Algorithm (following PySCF `hf.py` `init_guess_by_minao`)
///
/// 1. Build STO-3G basis for the same molecular geometry
/// 2. Compute cross-overlap `S_cross` between target and STO-3G: (n_target x n_sto)
/// 3. Compute target self-overlap `S_target`: (n_target x n_target)
/// 4. Solve `S_target @ C_proj = S_cross` for projection coefficients
/// 5. Set occupation vector from atomic electron configurations
/// 6. Build density: `D = C_proj @ diag(occ) @ C_proj^T`
///
/// # Arguments
///
/// * `basis` - Target basis set (e.g., 6-31G*)
/// * `s_target` - Pre-computed overlap matrix of target basis (row-major, n_target x n_target).
///   If `None`, it will be computed internally. When calling from SCF context where
///   S is already available, pass it to avoid redundant computation.
///
/// # Returns
///
/// Density matrix as a flattened column-major `Vec<f64>` of length n_target^2,
/// matching the format used by `ScfOutput.density_matrix`.
///
/// Returns `None` if the STO-3G basis cannot be built (unsupported elements)
/// or if the linear solve fails.
///
/// # References
///
/// - PySCF: `scf/hf.py` `init_guess_by_minao` (lines 347-473)
/// - PySCF: `scf/addons.py` `project_mo_nr2nr` (lines 580-603)
pub fn minao_initial_guess(basis: &BasisSet, s_target: Option<&[f64]>) -> Option<Vec<f64>> {
    // Step 1: Build STO-3G minimal basis for same geometry
    let atoms = basis.atoms.clone();
    let sto3g = match BasisSet::build(atoms, "sto-3g") {
        Ok(b) => b,
        Err(_) => return None, // STO-3G not available for some element
    };

    let n_cart = basis.n_basis;
    let n_sph = basis.n_basis_spherical();
    let n_sto = sto3g.n_basis; // STO-3G is always Cartesian (no d-functions)

    if n_sto == 0 || n_cart == 0 {
        return None;
    }

    // Determine whether to work in Cartesian or spherical representation.
    // If the caller provides an overlap matrix, infer from its size.
    // If the overlap is Cartesian (n_cart x n_cart), produce Cartesian density.
    // If the overlap is spherical (n_sph x n_sph), produce spherical density.
    // If no overlap is provided, use spherical if the basis has d-functions.
    let use_spherical = if let Some(s_flat) = s_target {
        let expected_cart = n_cart * n_cart;
        let expected_sph = n_sph * n_sph;
        if s_flat.len() == expected_sph && n_cart != n_sph {
            true
        } else if s_flat.len() == expected_cart {
            false
        } else {
            // Unexpected size -- this shouldn't happen
            return None;
        }
    } else {
        basis.has_spherical_difference()
    };

    let n_target = if use_spherical { n_sph } else { n_cart };

    // Step 2: Compute cross-overlap S_cross between target and STO-3G
    // S_cross[mu, nu] = <target_mu | sto3g_nu>
    // This is (n_target x n_sto) -- may need spherical transform on target side
    let s_cross_flat = if use_spherical {
        cross_overlap_matrix_spherical(basis, &sto3g)
    } else {
        cross_overlap_matrix(basis, &sto3g)
    };

    let s_cross = DMatrix::from_row_slice(n_target, n_sto, &s_cross_flat);

    // Step 3: Get target overlap matrix S_target (n_target x n_target)
    let s_target_mat = if let Some(s_flat) = s_target {
        // Use pre-computed overlap matrix from caller
        DMatrix::from_row_slice(n_target, n_target, s_flat)
    } else {
        // Compute from scratch
        let s_flat = if use_spherical {
            crate::integrals::overlap_matrix_spherical(basis)
        } else {
            crate::integrals::overlap_matrix(basis)
        };
        DMatrix::from_row_slice(n_target, n_target, &s_flat)
    };

    // Step 4: Solve S_target @ C_proj = S_cross for C_proj
    // C_proj = S_target^{-1} @ S_cross (n_target x n_sto)
    //
    // PySCF uses Cholesky decomposition: cho_solve(s22, s21 @ mo1)
    // We use LU decomposition via nalgebra for robustness (handles near-singular S)
    let c_proj = match s_target_mat.clone().lu().solve(&s_cross) {
        Some(c) => c,
        None => {
            // LU solve failed, try with pseudoinverse approach
            // This shouldn't happen for well-conditioned basis sets
            return None;
        }
    };

    // Step 5: Build occupation vector
    let mut occ = Vec::with_capacity(n_sto);
    for atom_idx in 0..sto3g.atoms.len() {
        let z = sto3g.atoms[atom_idx].atomic_number;
        let atom_shells: Vec<&crate::basis::ContractedShell> = sto3g
            .shells
            .iter()
            .filter(|s| s.atom_idx == atom_idx)
            .collect();
        let atom_occ = atom_occupation(z, &atom_shells);
        occ.extend_from_slice(&atom_occ);
    }

    // Verify occupation vector length matches STO-3G basis size
    if occ.len() != n_sto {
        #[cfg(not(target_arch = "wasm32"))]
        eprintln!(
            "MINAO: occupation vector length {} != STO-3G basis size {}",
            occ.len(),
            n_sto
        );
        return None;
    }

    // Step 6: Build density matrix D = C_proj @ diag(occ) @ C_proj^T
    //
    // PySCF: dm = lib.dot(mo*occ, mo.conj().T)
    // This is: D = (C_proj * occ) @ C_proj^T
    //
    // Efficiently: D[mu, nu] = sum_i occ[i] * C_proj[mu, i] * C_proj[nu, i]
    let mut density = DMatrix::zeros(n_target, n_target);
    for mu in 0..n_target {
        for nu in 0..=mu {
            let mut val = 0.0;
            for i in 0..n_sto {
                val += occ[i] * c_proj[(mu, i)] * c_proj[(nu, i)];
            }
            density[(mu, nu)] = val;
            density[(nu, mu)] = val;
        }
    }

    // Return density as column-major flat array (matching nalgebra convention)
    Some(density.as_slice().to_vec())
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::Atom;
    use approx::assert_abs_diff_eq;

    /// Build H2 at R=1.4 bohr
    fn h2_atoms() -> Vec<Atom> {
        vec![
            Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
            Atom::new(1, [0.0, 0.0, 1.4]).unwrap(),
        ]
    }

    /// Build H2O geometry (same as used in golden tests)
    fn h2o_atoms() -> Vec<Atom> {
        vec![
            Atom::new(8, [0.0, 0.0, 0.2173202]).unwrap(),
            Atom::new(1, [0.0, 1.4309847, -0.8692808]).unwrap(),
            Atom::new(1, [0.0, -1.4309847, -0.8692808]).unwrap(),
        ]
    }

    /// Build CH4 geometry (tetrahedral)
    fn ch4_atoms() -> Vec<Atom> {
        // Tetrahedral CH4 with C-H = 2.05 bohr
        let d = 2.05 / 3.0_f64.sqrt();
        vec![
            Atom::new(6, [0.0, 0.0, 0.0]).unwrap(),
            Atom::new(1, [d, d, d]).unwrap(),
            Atom::new(1, [-d, -d, d]).unwrap(),
            Atom::new(1, [-d, d, -d]).unwrap(),
            Atom::new(1, [d, -d, -d]).unwrap(),
        ]
    }

    // =========================================================================
    // Occupation tests
    // =========================================================================

    #[test]
    fn test_frac_occ_hydrogen() {
        let (ndocc, frac) = frac_occ(1, 0);
        assert_eq!(ndocc, 0);
        assert_abs_diff_eq!(frac, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_frac_occ_carbon_s() {
        // C (Z=6): config [4, 2, 0, 0]
        // s: 4 electrons, max per subshell = 2, ndocc = 4/2 = 2, frac = 0
        let (ndocc, frac) = frac_occ(6, 0);
        assert_eq!(ndocc, 2);
        assert_abs_diff_eq!(frac, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_frac_occ_carbon_p() {
        // C (Z=6): config [4, 2, 0, 0]
        // p: 2 electrons, max per subshell = 6, ndocc = 2/6 = 0
        // frac = (2/6 - 0) * 2 = 0.6667
        let (ndocc, frac) = frac_occ(6, 1);
        assert_eq!(ndocc, 0);
        assert_abs_diff_eq!(frac, 2.0 / 3.0, epsilon = 1e-12);
    }

    #[test]
    fn test_frac_occ_oxygen_p() {
        // O (Z=8): config [4, 4, 0, 0]
        // p: 4 electrons, max per subshell = 6, ndocc = 4/6 = 0
        // frac = (4/6 - 0) * 2 = 1.3333
        let (ndocc, frac) = frac_occ(8, 1);
        assert_eq!(ndocc, 0);
        assert_abs_diff_eq!(frac, 4.0 / 3.0, epsilon = 1e-12);
    }

    #[test]
    fn test_atom_occupation_hydrogen() {
        // H in STO-3G: 1 s-shell
        let h = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let basis = BasisSet::build(vec![h], "sto-3g").unwrap();
        let shells: Vec<&_> = basis.shells.iter().collect();
        let occ = atom_occupation(1, &shells);
        assert_eq!(occ.len(), 1);
        assert_abs_diff_eq!(occ[0], 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_atom_occupation_carbon() {
        // C in STO-3G: 1s, 2s, 2p (5 basis functions)
        let c = Atom::new(6, [0.0, 0.0, 0.0]).unwrap();
        let basis = BasisSet::build(vec![c], "sto-3g").unwrap();
        let shells: Vec<&_> = basis.shells.iter().collect();
        let occ = atom_occupation(6, &shells);
        // s-shells: 2 (1s) + 2 (2s) = 4 electrons, ndocc=2
        // p-shells: 2/3 per orbital (3 orbitals) = 2 electrons total
        assert_eq!(occ.len(), 5);
        assert_abs_diff_eq!(occ[0], 2.0, epsilon = 1e-12); // 1s
        assert_abs_diff_eq!(occ[1], 2.0, epsilon = 1e-12); // 2s
        assert_abs_diff_eq!(occ[2], 2.0 / 3.0, epsilon = 1e-12); // 2px
        assert_abs_diff_eq!(occ[3], 2.0 / 3.0, epsilon = 1e-12); // 2py
        assert_abs_diff_eq!(occ[4], 2.0 / 3.0, epsilon = 1e-12); // 2pz
                                                                 // Total: 4 + 2 = 6 electrons (correct for carbon)
        let total: f64 = occ.iter().sum();
        assert_abs_diff_eq!(total, 6.0, epsilon = 1e-12);
    }

    #[test]
    fn test_atom_occupation_oxygen() {
        // O in STO-3G: 1s, 2s, 2p (5 basis functions)
        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let basis = BasisSet::build(vec![o], "sto-3g").unwrap();
        let shells: Vec<&_> = basis.shells.iter().collect();
        let occ = atom_occupation(8, &shells);
        assert_eq!(occ.len(), 5);
        assert_abs_diff_eq!(occ[0], 2.0, epsilon = 1e-12); // 1s
        assert_abs_diff_eq!(occ[1], 2.0, epsilon = 1e-12); // 2s
        assert_abs_diff_eq!(occ[2], 4.0 / 3.0, epsilon = 1e-12); // 2px
        assert_abs_diff_eq!(occ[3], 4.0 / 3.0, epsilon = 1e-12); // 2py
        assert_abs_diff_eq!(occ[4], 4.0 / 3.0, epsilon = 1e-12); // 2pz
        let total: f64 = occ.iter().sum();
        assert_abs_diff_eq!(total, 8.0, epsilon = 1e-12);
    }

    // =========================================================================
    // Cross-overlap tests
    // =========================================================================

    #[test]
    fn test_cross_overlap_same_basis() {
        // Cross-overlap of a basis with itself should equal the self-overlap
        let atoms = h2_atoms();
        let basis = BasisSet::build(atoms, "sto-3g").unwrap();

        let s_self = crate::integrals::overlap_matrix(&basis);
        let s_cross = cross_overlap_matrix(&basis, &basis);

        assert_eq!(s_self.len(), s_cross.len());
        for (a, b) in s_self.iter().zip(s_cross.iter()) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_cross_overlap_rectangular() {
        // H2O: STO-3G (7 functions) vs 6-31G (13 functions, no d)
        let atoms = h2o_atoms();
        let basis_sto = BasisSet::build(atoms.clone(), "sto-3g").unwrap();
        let basis_631g = BasisSet::build(atoms, "6-31g").unwrap();

        let s_cross = cross_overlap_matrix(&basis_631g, &basis_sto);
        assert_eq!(s_cross.len(), basis_631g.n_basis * basis_sto.n_basis);

        // Check that diagonal-ish elements are non-zero (same-center overlap)
        // The first basis function in both should have large overlap
        assert!(s_cross[0].abs() > 0.1);
    }

    // =========================================================================
    // MINAO density matrix tests
    // =========================================================================

    #[test]
    fn test_minao_h2_sto3g() {
        // When target basis IS STO-3G, MINAO should recover
        // something close to the correct density
        let atoms = h2_atoms();
        let basis = BasisSet::build(atoms, "sto-3g").unwrap();
        let s = crate::integrals::overlap_matrix(&basis);

        let density =
            minao_initial_guess(&basis, Some(&s)).expect("MINAO should succeed for H2 STO-3G");

        let nbf = basis.n_basis;
        assert_eq!(density.len(), nbf * nbf);

        // Check symmetry
        let d_mat = DMatrix::from_column_slice(nbf, nbf, &density);
        for i in 0..nbf {
            for j in 0..nbf {
                assert_abs_diff_eq!(d_mat[(i, j)], d_mat[(j, i)], epsilon = 1e-12);
            }
        }

        // Check trace: Tr(D @ S) should equal number of electrons
        let s_mat = DMatrix::from_row_slice(nbf, nbf, &s);
        let ds = &d_mat * &s_mat;
        let trace: f64 = (0..nbf).map(|i| ds[(i, i)]).sum();
        assert_abs_diff_eq!(trace, basis.n_electrons as f64, epsilon = 0.5);
    }

    #[test]
    fn test_minao_h2o_631g_star() {
        // H2O with 6-31G* -- MINAO should produce valid density
        let atoms = h2o_atoms();
        let basis = BasisSet::build(atoms, "6-31g*").unwrap();

        let density =
            minao_initial_guess(&basis, None).expect("MINAO should succeed for H2O 6-31G*");

        let n_sph = basis.n_basis_spherical();
        assert_eq!(density.len(), n_sph * n_sph);

        // Check symmetry
        let d_mat = DMatrix::from_column_slice(n_sph, n_sph, &density);
        for i in 0..n_sph {
            for j in 0..n_sph {
                assert_abs_diff_eq!(d_mat[(i, j)], d_mat[(j, i)], epsilon = 1e-12);
            }
        }

        // Check electron count via Tr(D @ S) ~ n_electrons
        let s = crate::integrals::overlap_matrix_spherical(&basis);
        let s_mat = DMatrix::from_row_slice(n_sph, n_sph, &s);
        let ds = &d_mat * &s_mat;
        let trace: f64 = (0..n_sph).map(|i| ds[(i, i)]).sum();
        // MINAO trace won't be exact but should be close
        assert!(
            (trace - basis.n_electrons as f64).abs() < 2.0,
            "Tr(D*S)={} should be near {} electrons",
            trace,
            basis.n_electrons
        );
    }

    #[test]
    fn test_minao_ch4_631g_star() {
        // CH4 with 6-31G* -- the key test case
        let atoms = ch4_atoms();
        let basis = BasisSet::build(atoms, "6-31g*").unwrap();

        let density =
            minao_initial_guess(&basis, None).expect("MINAO should succeed for CH4 6-31G*");

        let n_sph = basis.n_basis_spherical();
        assert_eq!(density.len(), n_sph * n_sph);

        // Check symmetry
        let d_mat = DMatrix::from_column_slice(n_sph, n_sph, &density);
        for i in 0..n_sph {
            for j in 0..n_sph {
                assert_abs_diff_eq!(d_mat[(i, j)], d_mat[(j, i)], epsilon = 1e-12);
            }
        }

        // Check electron count
        let s = crate::integrals::overlap_matrix_spherical(&basis);
        let s_mat = DMatrix::from_row_slice(n_sph, n_sph, &s);
        let ds = &d_mat * &s_mat;
        let trace: f64 = (0..n_sph).map(|i| ds[(i, i)]).sum();
        assert!(
            (trace - 10.0).abs() < 2.0,
            "Tr(D*S)={} should be near 10 electrons for CH4",
            trace,
        );
    }

    #[test]
    fn test_minao_density_all_positive_diagonal() {
        // Density matrix diagonal should be non-negative
        let atoms = h2o_atoms();
        let basis = BasisSet::build(atoms, "sto-3g").unwrap();

        let density = minao_initial_guess(&basis, None).expect("MINAO should succeed");

        let nbf = basis.n_basis;
        let d_mat = DMatrix::from_column_slice(nbf, nbf, &density);
        for i in 0..nbf {
            assert!(
                d_mat[(i, i)] >= -1e-10,
                "D[{},{}] = {} should be non-negative",
                i,
                i,
                d_mat[(i, i)]
            );
        }
    }

    // =========================================================================
    // SCF convergence tests with MINAO initial guess
    // =========================================================================

    #[test]
    fn test_minao_improves_ch4_631g_star_rhf() {
        // CH4/6-31G* RHF converges in ~44 iterations with core Hamiltonian + DIIS,
        // or 53 iterations without DIIS. With MINAO + DIIS, it should converge
        // in ~26 iterations -- a ~40% reduction.
        //
        // Note: PySCF achieves 8 iterations because it uses the ANO (Atomic
        // Natural Orbital) basis for MINAO, which provides much better atomic
        // orbital shapes than STO-3G. Our STO-3G projection is still a major
        // improvement over core Hamiltonian.
        let atoms = ch4_atoms();
        let basis = BasisSet::build(atoms, "6-31g*").unwrap();

        // Build system with Cartesian integrals (matching test convention)
        let s_matrix = crate::integrals::overlap_matrix(&basis);
        let h_core = crate::integrals::hcore_matrix(&basis);
        let eri = crate::integrals::eri_compressed(&basis);

        let system = crate::scf::PresetSystem {
            system_id: "ch4_631gs".to_string(),
            label: "CH4 6-31G*".to_string(),
            nbf: basis.n_basis,
            nelec: basis.n_electrons,
            e_nuc: basis.nuclear_repulsion,
            s_matrix: s_matrix.clone(),
            h_core,
            eri_compressed: eri,
        };

        let config = crate::scf::ScfConfig {
            use_diis: true,
            ..crate::scf::ScfConfig::tight()
        };

        // Get MINAO density
        let minao_density = minao_initial_guess(&basis, Some(&s_matrix))
            .expect("MINAO should succeed for CH4 6-31G*");

        // Run RHF with MINAO initial guess
        let result = crate::scf::rhf_scf_with_guess(&system, &config, Some(&minao_density))
            .expect("RHF with MINAO should converge for CH4 6-31G*");

        assert!(result.converged);
        // MINAO + DIIS should converge in <=35 iterations
        assert!(
            result.iterations <= 35,
            "CH4/6-31G* RHF with MINAO should converge in <=35 iterations, got {}",
            result.iterations
        );

        // Also run without MINAO for comparison
        let result_core = crate::scf::rhf_scf(&system, &config)
            .expect("RHF with core Hamiltonian should converge for CH4 6-31G*");
        assert!(result_core.converged);

        // MINAO should use fewer iterations than core Hamiltonian
        eprintln!(
            "CH4/6-31G* RHF iterations: MINAO={}, core={}",
            result.iterations, result_core.iterations
        );
        assert!(
            result.iterations < result_core.iterations,
            "MINAO ({}) should converge faster than core Hamiltonian ({})",
            result.iterations,
            result_core.iterations
        );

        // Both should give the same final energy
        assert_abs_diff_eq!(
            result.energy_total,
            result_core.energy_total,
            epsilon = 1e-8
        );
    }

    #[test]
    fn test_minao_h2o_631g_star_rhf() {
        // H2O/6-31G* RHF with MINAO should also converge faster
        let atoms = h2o_atoms();
        let basis = BasisSet::build(atoms, "6-31g*").unwrap();

        let s_matrix = crate::integrals::overlap_matrix(&basis);
        let h_core = crate::integrals::hcore_matrix(&basis);
        let eri = crate::integrals::eri_compressed(&basis);

        let system = crate::scf::PresetSystem {
            system_id: "h2o_631gs".to_string(),
            label: "H2O 6-31G*".to_string(),
            nbf: basis.n_basis,
            nelec: basis.n_electrons,
            e_nuc: basis.nuclear_repulsion,
            s_matrix: s_matrix.clone(),
            h_core,
            eri_compressed: eri,
        };

        let config = crate::scf::ScfConfig {
            use_diis: true,
            ..crate::scf::ScfConfig::tight()
        };

        let minao_density = minao_initial_guess(&basis, Some(&s_matrix))
            .expect("MINAO should succeed for H2O 6-31G*");

        let result = crate::scf::rhf_scf_with_guess(&system, &config, Some(&minao_density))
            .expect("RHF with MINAO should converge for H2O 6-31G*");

        assert!(result.converged);
        // Should converge quickly
        assert!(
            result.iterations <= 15,
            "H2O/6-31G* RHF with MINAO should converge in <=15 iterations, got {}",
            result.iterations
        );
    }
}
