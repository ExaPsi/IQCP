//! Kinetic energy integral computation via Obara-Saika recurrence relations
//!
//! This module computes kinetic energy integrals `T_ij = ⟨i|-1/2 nabla^2|j⟩` between Gaussian
//! basis functions using the Obara-Saika (OS) scheme.
//!
//! # Algorithm
//!
//! The kinetic energy integral is computed by expressing it in terms of overlap integrals.
//! For each Cartesian direction, the 1D kinetic contribution is:
//!
//! ```text
//! T_1D(a, b) = -2*beta^2 * S_1D(a, b+2)
//!            + beta * (2b + 1) * S_1D(a, b)
//!            - (1/2) * b * (b - 1) * S_1D(a, b-2)
//! ```
//!
//! where S_1D is the 1D overlap integral and beta is the exponent on center B.
//!
//! The total 3D kinetic integral is:
//! ```text
//! T_3D = T_x * S_y * S_z + S_x * T_y * S_z + S_x * S_y * T_z
//! ```
//!
//! # References
//!
//! - Obara & Saika (1986), J. Chem. Phys. 84, 3963, Eq. 3.22
//! - Helgaker, Jorgensen & Olsen (2000), Molecular Electronic-Structure Theory, Ch. 9.3
//! - libcint implementation: `references/libcint/src/g1e.c`
//! - PySCF: `references/pyscf/pyscf/gto/moleintor.py`
//!
//! # Example
//!
//! ```rust
//! use qc_core::basis::{Atom, BasisSet};
//! use qc_core::integrals::kinetic_matrix;
//!
//! // Build H2 molecule
//! let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
//! let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
//! let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
//!
//! // Compute kinetic matrix
//! let t = kinetic_matrix(&basis);
//!
//! // T is a 2x2 symmetric matrix (flattened row-major)
//! assert_eq!(t.len(), 4);
//! // Diagonal elements are positive (kinetic energy is positive)
//! assert!(t[0] > 0.0);  // T[0,0]
//! assert!(t[3] > 0.0);  // T[1,1]
//! ```

use super::cartesian::{cartesian_components, CartesianPower};
use super::overlap::{cartesian_gaussian_normalization, overlap_1d};
use super::GaussianProduct;
use crate::basis::{BasisSet, ContractedShell};

// =============================================================================
// 1D Kinetic Integral Component
// =============================================================================

/// Compute 1D kinetic integral component using overlap integrals
///
/// This computes the 1D kinetic energy contribution:
/// ```text
/// T_1D(a, b) = -2*beta^2 * S_1D(a, b+2)
///            + beta * (2b + 1) * S_1D(a, b)
///            - (1/2) * b * (b - 1) * S_1D(a, b-2)
/// ```
///
/// where S_1D is the 1D overlap integral from overlap.rs.
///
/// # Arguments
///
/// * `pa` - P - A (product center minus bra center)
/// * `ab` - B - A (ket center minus bra center)
/// * `one_over_2p` - 1/(2p) where p = alpha + beta
/// * `a` - Angular momentum on center A
/// * `b` - Angular momentum on center B
/// * `beta` - Exponent on center B
///
/// # Reference
///
/// Obara & Saika (1986), Eq. 3.22
#[inline]
pub(crate) fn kinetic_1d(pa: f64, ab: f64, one_over_2p: f64, a: i32, b: i32, beta: f64) -> f64 {
    // Term 1: -2*beta^2 * S_1D(a, b+2)
    let term1 = -2.0 * beta * beta * overlap_1d(pa, ab, one_over_2p, a, b + 2);

    // Term 2: beta * (2b + 1) * S_1D(a, b)
    let term2 = beta * (2.0 * b as f64 + 1.0) * overlap_1d(pa, ab, one_over_2p, a, b);

    // Term 3: -1/2 * b * (b - 1) * S_1D(a, b-2)
    // This term is zero if b < 2
    let term3 = if b >= 2 {
        -0.5 * (b as f64) * ((b - 1) as f64) * overlap_1d(pa, ab, one_over_2p, a, b - 2)
    } else {
        0.0
    };

    term1 + term2 + term3
}

// =============================================================================
// Primitive Kinetic Integral
// =============================================================================

/// Compute the kinetic energy integral between two primitive Cartesian Gaussians
///
/// This implements the Obara-Saika scheme to compute:
/// ```text
/// T = ⟨G_a|-1/2 nabla^2|G_b⟩
/// ```
/// where `G_a` and `G_b` are primitive Gaussians with angular momentum
/// specified by `a_powers` and `b_powers`.
///
/// # Arguments
///
/// * `gp` - Pre-computed Gaussian product data
/// * `a_powers` - Cartesian powers (i, j, k) for bra Gaussian
/// * `b_powers` - Cartesian powers (i, j, k) for ket Gaussian
/// * `beta` - Exponent of the ket Gaussian (center B)
///
/// # Returns
///
/// The unnormalized kinetic energy integral value
///
/// # Algorithm
///
/// The 3D kinetic integral separates as:
/// ```text
/// T_3D = T_x * S_y * S_z + S_x * T_y * S_z + S_x * S_y * T_z
/// ```
/// where T_i is the 1D kinetic contribution and S_i is the 1D overlap.
///
/// # Reference
///
/// Obara & Saika (1986), Eq. 3.22; libcint g1e.c
pub fn primitive_kinetic(
    gp: &GaussianProduct,
    a_powers: &CartesianPower,
    b_powers: &CartesianPower,
    beta: f64,
) -> f64 {
    let a_x = a_powers.i as i32;
    let a_y = a_powers.j as i32;
    let a_z = a_powers.k as i32;

    let b_x = b_powers.i as i32;
    let b_y = b_powers.j as i32;
    let b_z = b_powers.k as i32;

    // Compute 1D overlap integrals for each direction
    let s_x = overlap_1d(gp.pa[0], gp.ab[0], gp.one_over_2p, a_x, b_x);
    let s_y = overlap_1d(gp.pa[1], gp.ab[1], gp.one_over_2p, a_y, b_y);
    let s_z = overlap_1d(gp.pa[2], gp.ab[2], gp.one_over_2p, a_z, b_z);

    // Compute 1D kinetic integrals for each direction
    let t_x = kinetic_1d(gp.pa[0], gp.ab[0], gp.one_over_2p, a_x, b_x, beta);
    let t_y = kinetic_1d(gp.pa[1], gp.ab[1], gp.one_over_2p, a_y, b_y, beta);
    let t_z = kinetic_1d(gp.pa[2], gp.ab[2], gp.one_over_2p, a_z, b_z, beta);

    // Total kinetic integral: T_x*S_y*S_z + S_x*T_y*S_z + S_x*S_y*T_z
    // Multiply by the [s|s] prefactor
    gp.ss_integral * (t_x * s_y * s_z + s_x * t_y * s_z + s_x * s_y * t_z)
}

// =============================================================================
// Shell Kinetic Integral
// =============================================================================

/// Compute kinetic energy integrals between two contracted shells
///
/// This computes all `n_a * n_b` kinetic energy integrals between the Cartesian
/// components of shells A and B:
/// ```text
/// T_ij = sum_{p in A} sum_{q in B} c_p * N_p * c_q * N_q * T[p|q]
/// ```
/// where the sum runs over primitive pairs and N are normalization factors.
///
/// # Arguments
///
/// * `shell_a` - First contracted shell (bra)
/// * `shell_b` - Second contracted shell (ket)
///
/// # Returns
///
/// Vector of kinetic energy integrals in row-major order:
/// `[T(a0,b0), T(a0,b1), ..., T(a0,bn), T(a1,b0), ...]`
/// where a0, a1, ... are Cartesian components of shell A.
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{ContractedShell, AngularMomentum, GaussianPrimitive};
/// use qc_core::integrals::shell_kinetic;
///
/// // Create two s-shells at different centers
/// let h_prims = vec![
///     GaussianPrimitive::new(3.425251, 0.154329),
///     GaussianPrimitive::new(0.623914, 0.535328),
///     GaussianPrimitive::new(0.168855, 0.444635),
/// ];
/// let shell_a = ContractedShell::new(AngularMomentum::S, h_prims.clone(), [0.0, 0.0, 0.0], 0);
/// let shell_b = ContractedShell::new(AngularMomentum::S, h_prims, [0.0, 0.0, 1.4], 0);
///
/// let t = shell_kinetic(&shell_a, &shell_b);
/// assert_eq!(t.len(), 1);  // 1x1 for s-s kinetic integral
/// ```
pub fn shell_kinetic(shell_a: &ContractedShell, shell_b: &ContractedShell) -> Vec<f64> {
    let l_a = shell_a.l_value();
    let l_b = shell_b.l_value();

    // Get Cartesian components for each shell
    let comps_a = cartesian_components(l_a).expect("Angular momentum within supported range");
    let comps_b = cartesian_components(l_b).expect("Angular momentum within supported range");

    let n_a = comps_a.len();
    let n_b = comps_b.len();

    // Output array: n_a x n_b integrals
    let mut integrals = vec![0.0; n_a * n_b];

    // Loop over primitive pairs
    for prim_a in &shell_a.primitives {
        for prim_b in &shell_b.primitives {
            // Compute Gaussian product data
            let gp = GaussianProduct::new(
                prim_a.exponent,
                &shell_a.center,
                prim_b.exponent,
                &shell_b.center,
            );

            // Raw contraction coefficient product
            let coef = prim_a.coefficient * prim_b.coefficient;

            // Loop over Cartesian components
            for (i, pow_a) in comps_a.iter().enumerate() {
                for (j, pow_b) in comps_b.iter().enumerate() {
                    // Compute primitive integral (unnormalized)
                    let prim_integral = primitive_kinetic(&gp, pow_a, pow_b, prim_b.exponent);

                    // Apply Cartesian Gaussian normalization for each primitive
                    let norm_a = cartesian_gaussian_normalization(prim_a.exponent, pow_a);
                    let norm_b = cartesian_gaussian_normalization(prim_b.exponent, pow_b);

                    // Add contribution to contracted integral
                    integrals[i * n_b + j] += coef * norm_a * norm_b * prim_integral;
                }
            }
        }
    }

    integrals
}

// =============================================================================
// Kinetic Matrix
// =============================================================================

/// Compute the full kinetic energy matrix for a basis set
///
/// Returns the N x N kinetic energy matrix T where N is the total number of
/// basis functions. The matrix is symmetric: T_ij = T_ji.
///
/// # Arguments
///
/// * `basis` - The molecular basis set
///
/// # Returns
///
/// Vector of length N*N containing the flattened row-major kinetic matrix
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{Atom, BasisSet};
/// use qc_core::integrals::kinetic_matrix;
///
/// // H2 molecule
/// let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
/// let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
/// let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
///
/// let t = kinetic_matrix(&basis);
/// assert_eq!(t.len(), 4);  // 2x2 matrix
///
/// // Check symmetry
/// assert!((t[1] - t[2]).abs() < 1e-15);
///
/// // Check positive diagonal (kinetic energy is positive)
/// assert!(t[0] > 0.0);
/// assert!(t[3] > 0.0);
/// ```
pub fn kinetic_matrix(basis: &BasisSet) -> Vec<f64> {
    let n = basis.n_basis;
    let mut t_matrix = vec![0.0; n * n];

    // Iterate over shell pairs
    let mut mu = 0; // Basis function index for shell A
    for (i, shell_a) in basis.shells.iter().enumerate() {
        let n_a = shell_a.n_basis_functions();

        let mut nu = 0; // Basis function index for shell B
        for (j, shell_b) in basis.shells.iter().enumerate() {
            let n_b = shell_b.n_basis_functions();

            // Only compute upper triangle (i <= j) and symmetrize
            if i <= j {
                // Compute shell block
                let block = shell_kinetic(shell_a, shell_b);

                // Copy to matrix
                for ia in 0..n_a {
                    for ib in 0..n_b {
                        let val = block[ia * n_b + ib];
                        // Upper triangle
                        t_matrix[(mu + ia) * n + (nu + ib)] = val;
                        // Lower triangle (symmetry)
                        if i != j {
                            t_matrix[(nu + ib) * n + (mu + ia)] = val;
                        }
                    }
                }
            }

            nu += n_b;
        }
        mu += n_a;
    }

    t_matrix
}

// =============================================================================
// Spherical Kinetic Matrix
// =============================================================================

/// Compute the kinetic energy matrix in spherical harmonic basis
///
/// This computes the kinetic matrix and transforms it from Cartesian to
/// spherical harmonic basis. For basis sets without d-orbitals or higher,
/// this is identical to `kinetic_matrix()`.
///
/// # Arguments
///
/// * `basis` - The molecular basis set
///
/// # Returns
///
/// Vector of length N_sph * N_sph containing the flattened row-major
/// spherical kinetic matrix, where N_sph is the number of spherical
/// basis functions.
pub fn kinetic_matrix_spherical(basis: &BasisSet) -> Vec<f64> {
    use super::spherical::transform_one_electron_matrix;

    // If no spherical/Cartesian difference, just return the Cartesian matrix
    if !basis.has_spherical_difference() {
        return kinetic_matrix(basis);
    }

    // Compute Cartesian matrix and transform
    let cart_matrix = kinetic_matrix(basis);
    transform_one_electron_matrix(&cart_matrix, &basis.shells)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;
    use crate::basis::{AngularMomentum, Atom, GaussianPrimitive};
    use approx::assert_abs_diff_eq;

    // Tolerance for golden test comparisons
    // Similar to overlap.rs, this accounts for accumulated floating-point errors
    const TOL: f64 = 1e-9;

    // -------------------------------------------------------------------------
    // Utility functions
    // -------------------------------------------------------------------------

    fn h_sto3g_primitives() -> Vec<GaussianPrimitive> {
        // Hydrogen STO-3G from PySCF
        vec![
            GaussianPrimitive::new(3.4252509100, 0.1543289707),
            GaussianPrimitive::new(0.6239137300, 0.5353281424),
            GaussianPrimitive::new(0.1688554000, 0.4446345420),
        ]
    }

    // -------------------------------------------------------------------------
    // 1D kinetic integral tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_kinetic_1d_ss() {
        // For s-s (a=0, b=0), the kinetic contribution is:
        // T_1D(0, 0) = -2*beta^2 * S(0,2) + beta * 1 * S(0,0)
        // (the b-2 term is zero since b < 2)
        let pa = 0.0;
        let ab = 0.0;
        let beta = 1.0;
        let p = 2.0 * beta; // assuming alpha = beta
        let one_over_2p = 0.5 / p;

        let t = kinetic_1d(pa, ab, one_over_2p, 0, 0, beta);

        // S(0,0) = 1.0
        // S(0,2) = 1/(2p) (from VRR)
        // T = -2*1*1/(2*2) + 1*1*1 = -0.25 + 1.0 = 0.75
        let s_02 = one_over_2p; // [2|0] for same center = 1/(2p) from VRR
        let expected = -2.0 * beta * beta * s_02 + beta * 1.0;

        assert_abs_diff_eq!(t, expected, epsilon = TOL);
    }

    #[test]
    fn test_kinetic_1d_with_angular_momentum() {
        // Test with b >= 2 to engage all three terms
        let pa = 0.3;
        let ab = 1.0;
        let beta = 1.5;
        let p = 2.5;
        let one_over_2p = 0.5 / p;

        let a = 1;
        let b = 2;

        let t = kinetic_1d(pa, ab, one_over_2p, a, b, beta);

        // Verify it's a reasonable value (not NaN or infinite)
        assert!(t.is_finite());
    }

    // -------------------------------------------------------------------------
    // Primitive kinetic integral tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_primitive_kinetic_ss_same_center() {
        // Two s-type primitives at the same center
        let alpha = 1.0;
        let beta = 1.0;
        let center = [0.0, 0.0, 0.0];
        let gp = GaussianProduct::new(alpha, &center, beta, &center);

        let s = CartesianPower::new(0, 0, 0);
        let t = primitive_kinetic(&gp, &s, &s, beta);

        // For same center s-s, the kinetic integral has a known form:
        // T = 3 * alpha * beta / (alpha + beta) * (pi / (alpha + beta))^(3/2)
        // For alpha = beta = 1: T = 3 * 1 / 2 * (pi/2)^(3/2)
        let p = alpha + beta;
        let expected = 3.0 * alpha * beta / p * (std::f64::consts::PI / p).powf(1.5);

        assert_abs_diff_eq!(t, expected, epsilon = TOL);
    }

    #[test]
    fn test_primitive_kinetic_positive_definite() {
        // Kinetic energy should be positive for any non-zero function
        let alpha = 2.0;
        let beta = 1.5;
        let center_a = [0.0, 0.0, 0.0];
        let center_b = [0.5, 0.0, 0.0];
        let gp = GaussianProduct::new(alpha, &center_a, beta, &center_b);

        let s = CartesianPower::new(0, 0, 0);
        let t = primitive_kinetic(&gp, &s, &s, beta);

        // Kinetic matrix should be positive semi-definite
        // For s-s, diagonal elements should be positive
        assert!(t > 0.0, "Kinetic integral should be positive, got {}", t);
    }

    // -------------------------------------------------------------------------
    // Shell kinetic integral tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_shell_kinetic_same_center() {
        // Single H 1s shell kinetic energy with itself
        let prims = h_sto3g_primitives();
        let shell = ContractedShell::new(AngularMomentum::S, prims, [0.0, 0.0, 0.0], 0);

        let t = shell_kinetic(&shell, &shell);

        assert_eq!(t.len(), 1);
        // Kinetic energy should be positive
        assert!(
            t[0] > 0.0,
            "Kinetic energy should be positive, got {}",
            t[0]
        );
    }

    // -------------------------------------------------------------------------
    // Golden tests: H2 STO-3G
    // -------------------------------------------------------------------------

    #[test]
    fn test_golden_h2_sto3g_kinetic() {
        // Reference: PySCF 2.11.0
        //
        // H2 at 1.3984 Bohr separation (same geometry as overlap tests)
        // Expected kinetic matrix (2x2):
        // [[0.76003188, 0.23695942],
        //  [0.23695942, 0.76003188]]
        //
        // Generated with:
        // ```python
        // from pyscf import gto
        // mol = gto.Mole()
        // mol.atom = 'H 0 0 0; H 0 0 1.3984'
        // mol.basis = 'sto-3g'
        // mol.unit = 'B'
        // mol.build()
        // print(mol.intor('int1e_kin'))
        // ```

        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let t = kinetic_matrix(&basis);

        // Reference values from PySCF
        let t_00 = 7.600318835666091e-01;
        let t_01 = 2.369594228630360e-01;

        assert_abs_diff_eq!(t[0], t_00, epsilon = TOL); // T[0,0]
        assert_abs_diff_eq!(t[1], t_01, epsilon = TOL); // T[0,1]
        assert_abs_diff_eq!(t[2], t_01, epsilon = TOL); // T[1,0] (symmetry)
        assert_abs_diff_eq!(t[3], t_00, epsilon = TOL); // T[1,1]
    }

    // -------------------------------------------------------------------------
    // Golden tests: H2O STO-3G
    // -------------------------------------------------------------------------

    #[test]
    fn test_golden_h2o_sto3g_kinetic() {
        // Reference: PySCF 2.11.0
        //
        // H2O geometry (Bohr):
        // O: [0, 0, 0.2216656]
        // H1: [0, 1.4309295, -0.8866625]
        // H2: [0, -1.4309295, -0.8866625]
        //
        // 7x7 kinetic matrix
        //
        // Generated with:
        // ```python
        // from pyscf import gto
        // mol = gto.Mole()
        // mol.atom = '''
        // O   0.0000000   0.0000000   0.2216656303019316
        // H   0.0000000   1.4309295343303710  -0.8866625212077263
        // H   0.0000000  -1.4309295343303710  -0.8866625212077263
        // '''
        // mol.basis = 'sto-3g'
        // mol.unit = 'B'
        // mol.build()
        // print(mol.intor('int1e_kin').flatten())
        // ```

        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        let t = kinetic_matrix(&basis);
        assert_eq!(t.len(), 49); // 7x7

        // Reference values from PySCF 2.11.0 (flattened row-major)
        // Generated with:
        // mol.atom = 'O 0 0 0.2216656303019316; H 0 1.430929534330371 -0.8866625212077263; H 0 -1.430929534330371 -0.8866625212077263'
        // mol.basis = 'sto-3g'; mol.unit = 'B'; mol.build()
        // mol.intor('int1e_kin').flatten()
        let reference = [
            2.900319994553957e+01,
            -1.680109393164921e-01,
            0.0,
            0.0,
            0.0,
            -2.552810828575374e-03,
            -2.552810828575374e-03,
            -1.680109393164919e-01,
            8.081279549303472e-01,
            0.0,
            0.0,
            0.0,
            1.283275865537280e-01,
            1.283275865537280e-01,
            0.0,
            0.0,
            2.528731198194763e+00,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2.528731198194763e+00,
            0.0,
            2.243872717590316e-01,
            -2.243872717590316e-01,
            0.0,
            0.0,
            0.0,
            0.0,
            2.528731198194763e+00,
            -1.737994248943669e-01,
            -1.737994248943669e-01,
            -2.552810828575348e-03,
            1.283275865537280e-01,
            0.0,
            2.243872717590316e-01,
            -1.737994248943670e-01,
            7.600318835666091e-01,
            8.441933805978971e-03,
            -2.552810828575348e-03,
            1.283275865537280e-01,
            0.0,
            -2.243872717590316e-01,
            -1.737994248943670e-01,
            8.441933805978971e-03,
            7.600318835666091e-01,
        ];

        // Check all elements
        // Note: The H2O matrix has larger elements (up to ~29 Ha for O 1s kinetic energy)
        // which can accumulate more floating-point error. We use a slightly relaxed
        // absolute tolerance that still represents excellent numerical agreement
        // (relative error < 1e-10 for all elements).
        let h2o_tol = 2e-9;
        for (i, &ref_val) in reference.iter().enumerate() {
            assert!(
                (t[i] - ref_val).abs() < h2o_tol,
                "Mismatch at index {}: computed {} vs reference {} (error: {:.2e})",
                i,
                t[i],
                ref_val,
                (t[i] - ref_val).abs()
            );
        }
    }

    // -------------------------------------------------------------------------
    // Matrix properties tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_kinetic_matrix_symmetry() {
        // Kinetic matrix should be symmetric: T_ij = T_ji
        let o = Atom::new(8, [0.0, 0.0, 0.0]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.43, -0.47]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.43, -0.47]).unwrap();
        let basis = BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap();

        let t = kinetic_matrix(&basis);
        let n = basis.n_basis;

        for i in 0..n {
            for j in 0..n {
                assert_abs_diff_eq!(t[i * n + j], t[j * n + i], epsilon = 1e-14);
            }
        }
    }

    #[test]
    fn test_kinetic_matrix_positive_diagonal() {
        // Kinetic energy diagonal elements should be positive
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();

        let t = kinetic_matrix(&basis);
        let n = basis.n_basis;

        for i in 0..n {
            assert!(
                t[i * n + i] > 0.0,
                "Diagonal element {} should be positive, got {}",
                i,
                t[i * n + i]
            );
        }
    }

    #[test]
    fn test_kinetic_integral_scaling() {
        // Kinetic energy scales with exponent squared
        // For tighter Gaussians (higher exponent), kinetic energy should be larger
        let center = [0.0, 0.0, 0.0];

        // Tight Gaussian
        let tight_prims = vec![GaussianPrimitive::new(5.0, 1.0)];
        let tight_shell = ContractedShell::new(AngularMomentum::S, tight_prims, center, 0);

        // Diffuse Gaussian
        let diffuse_prims = vec![GaussianPrimitive::new(0.5, 1.0)];
        let diffuse_shell = ContractedShell::new(AngularMomentum::S, diffuse_prims, center, 0);

        let t_tight = shell_kinetic(&tight_shell, &tight_shell)[0];
        let t_diffuse = shell_kinetic(&diffuse_shell, &diffuse_shell)[0];

        // Kinetic energy is proportional to exponent for Gaussians
        // T ~ alpha for s-type Gaussians
        assert!(
            t_tight > t_diffuse,
            "Tight Gaussian ({}) should have larger kinetic energy than diffuse ({})",
            t_tight,
            t_diffuse
        );
    }
}
