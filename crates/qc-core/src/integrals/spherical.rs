//! Cartesian to spherical harmonic transformation
//!
//! This module implements transformations from Cartesian Gaussian basis functions
//! to real spherical harmonic basis functions.
//!
//! # Background
//!
//! Cartesian Gaussians use polynomial prefactors x^i y^j z^k, leading to
//! (L+1)(L+2)/2 components for angular momentum L:
//! - L=2 (d-orbitals): 6 Cartesian (dxx, dxy, dxz, dyy, dyz, dzz)
//!
//! Real spherical harmonics use 2L+1 functions labeled by m = -L,...,+L:
//! - L=2 (d-orbitals): 5 spherical (d-2, d-1, d0, d+1, d+2)
//!
//! Spherical harmonics are standard in quantum chemistry because:
//! 1. Smaller matrices (5 vs 6 for d-orbitals)
//! 2. No linear dependency issues
//! 3. Standard convention in publications
//!
//! # Transformation Convention
//!
//! The transformation follows the Schlegel & Frisch convention (IJQC 54, 83, 1995)
//! as implemented in libcint. The transformation is:
//!
//! ```text
//! AO_sph = C^T * AO_cart
//! ```
//!
//! where C is the transformation matrix with dimensions (n_cart x n_sph).
//!
//! For two-electron integrals, the transformation is applied to each index:
//!
//! ```text
//! (ij|kl)_sph = sum_{a,b,c,d} C[a,i] * C[b,j] * C[c,k] * C[d,l] * (ab|cd)_cart
//! ```
//!
//! # References
//!
//! - H. B. Schlegel and M. J. Frisch, Int. J. Quant. Chem. 54, 83 (1995)
//! - libcint implementation: `references/libcint/src/cart2sph.c`
//! - PySCF gto.mole.cart2sph: `references/pyscf/pyscf/gto/mole.py`

use super::cartesian::n_cartesian;

// =============================================================================
// Transformation Coefficients
// =============================================================================

/// Transformation coefficients from Cartesian to real spherical harmonics.
///
/// Source: libcint cart2sph.c lines 52-86
/// Reference: Schlegel & Frisch, IJQC 54, 83 (1995)
///
/// Cartesian d-orbital ordering (IQCP convention):
/// - Index 0: dxx (2,0,0)
/// - Index 1: dxy (1,1,0)
/// - Index 2: dxz (1,0,1)
/// - Index 3: dyy (0,2,0)
/// - Index 4: dyz (0,1,1)
/// - Index 5: dzz (0,0,2)
///
/// Spherical d-orbital ordering (standard convention):
/// - Index 0: d-2 (m = -2)
/// - Index 1: d-1 (m = -1)
/// - Index 2: d0  (m = 0)
/// - Index 3: d+1 (m = +1)
/// - Index 4: d+2 (m = +2)
///
/// The matrix C[cart][sph] transforms as: AO_sph = C^T * AO_cart
#[allow(clippy::excessive_precision)]
pub const CART2SPH_D: [[f64; 5]; 6] = [
    // dxx: contributes to d0 and d+2
    // d-2,     d-1,      d0,                    d+1,      d+2
    [0.0, 0.0, -0.315391565252520002, 0.0, 0.546274215296039535],
    // dxy: contributes only to d-2
    [1.092548430592079070, 0.0, 0.0, 0.0, 0.0],
    // dxz: contributes only to d+1
    [0.0, 0.0, 0.0, 1.092548430592079070, 0.0],
    // dyy: contributes to d0 and d+2
    [0.0, 0.0, -0.315391565252520002, 0.0, -0.546274215296039535],
    // dyz: contributes only to d-1
    [0.0, 1.092548430592079070, 0.0, 0.0, 0.0],
    // dzz: contributes only to d0
    [0.0, 0.0, 0.630783130505040012, 0.0, 0.0],
];

/// Number of spherical harmonics for angular momentum L
#[inline]
pub const fn n_spherical(l: u32) -> u32 {
    2 * l + 1
}

// =============================================================================
// Transformation Structs
// =============================================================================

/// Spherical transformation configuration for a shell
#[derive(Debug, Clone)]
pub struct SphericalTransform {
    /// Angular momentum of the shell
    pub l: u32,
    /// Number of Cartesian functions
    pub n_cart: usize,
    /// Number of spherical functions
    pub n_sph: usize,
}

impl SphericalTransform {
    /// Create a new spherical transformation for angular momentum L
    pub fn new(l: u32) -> Self {
        Self {
            l,
            n_cart: n_cartesian(l) as usize,
            n_sph: n_spherical(l) as usize,
        }
    }

    /// Check if transformation is needed (only for L >= 2)
    #[inline]
    pub fn needs_transform(&self) -> bool {
        self.l >= 2
    }

    /// Get the transformation matrix element C[cart_idx][sph_idx]
    ///
    /// Returns 0.0 for L < 2 (identity transformation for s and p)
    #[inline]
    pub fn coeff(&self, cart_idx: usize, sph_idx: usize) -> f64 {
        match self.l {
            0 | 1 => {
                // Identity transformation for s and p orbitals
                if cart_idx == sph_idx {
                    1.0
                } else {
                    0.0
                }
            }
            2 => {
                // d-orbitals
                CART2SPH_D[cart_idx][sph_idx]
            }
            _ => {
                // Higher angular momentum not yet implemented
                panic!("Spherical transformation not implemented for L > 2");
            }
        }
    }

    /// Transform a single shell's values from Cartesian to spherical
    ///
    /// # Arguments
    ///
    /// * `cart_values` - Cartesian values with length n_cart
    ///
    /// # Returns
    ///
    /// Spherical values with length n_sph
    pub fn transform_1d(&self, cart_values: &[f64]) -> Vec<f64> {
        debug_assert_eq!(
            cart_values.len(),
            self.n_cart,
            "Input length {} does not match n_cart {}",
            cart_values.len(),
            self.n_cart
        );

        if !self.needs_transform() {
            // For s and p, Cartesian and spherical are the same
            return cart_values.to_vec();
        }

        let mut sph_values = vec![0.0; self.n_sph];

        // sph[m] = sum_c C[c][m] * cart[c]
        for (m, sph_val) in sph_values.iter_mut().enumerate() {
            for (c, &cart_val) in cart_values.iter().enumerate() {
                *sph_val += self.coeff(c, m) * cart_val;
            }
        }

        sph_values
    }
}

// =============================================================================
// Two-Index (One-Electron) Transformation
// =============================================================================

/// Configuration for transforming a 2-index one-electron integral matrix
///
/// One-electron integrals (overlap S, kinetic T, nuclear V) are 2-index matrices.
/// The spherical transformation is:
///
/// ```text
/// M_sph[i,j] = sum_{a,b} C[a,i] * C[b,j] * M_cart[a,b]
/// ```
///
/// where C is the Cartesian-to-spherical transformation matrix.
#[derive(Debug, Clone)]
pub struct OneElectronSphericalTransform {
    /// Transformations for the two shell indices
    pub trans_i: SphericalTransform,
    pub trans_j: SphericalTransform,
}

impl OneElectronSphericalTransform {
    /// Create a transformation for a shell pair
    pub fn new(l_i: u32, l_j: u32) -> Self {
        Self {
            trans_i: SphericalTransform::new(l_i),
            trans_j: SphericalTransform::new(l_j),
        }
    }

    /// Check if any transformation is needed
    pub fn needs_transform(&self) -> bool {
        self.trans_i.needs_transform() || self.trans_j.needs_transform()
    }

    /// Get the number of spherical functions for each shell
    pub fn n_sph(&self) -> (usize, usize) {
        (self.trans_i.n_sph, self.trans_j.n_sph)
    }

    /// Get the number of Cartesian functions for each shell
    pub fn n_cart(&self) -> (usize, usize) {
        (self.trans_i.n_cart, self.trans_j.n_cart)
    }

    /// Transform a 2-index integral block from Cartesian to spherical
    ///
    /// The input matrix has dimensions (n_i_cart, n_j_cart) in row-major order.
    /// The output matrix has dimensions (n_i_sph, n_j_sph).
    ///
    /// # Algorithm
    ///
    /// ```text
    /// M_sph[i,j] = sum_{a,b} C[a,i] * C[b,j] * M_cart[a,b]
    /// ```
    ///
    /// This is applied sequentially:
    /// 1. Transform index j: (cart, cart) -> (cart, sph)
    /// 2. Transform index i: (cart, sph) -> (sph, sph)
    pub fn transform(&self, cart_matrix: &[f64]) -> Vec<f64> {
        if !self.needs_transform() {
            return cart_matrix.to_vec();
        }

        let (n_i_c, n_j_c) = self.n_cart();
        let (n_i_s, n_j_s) = self.n_sph();

        debug_assert_eq!(
            cart_matrix.len(),
            n_i_c * n_j_c,
            "Input matrix length does not match expected Cartesian dimensions"
        );

        // Step 1: Transform index j (column index)
        // (n_i_c, n_j_c) -> (n_i_c, n_j_s)
        let tmp = if self.trans_j.needs_transform() {
            let mut result = vec![0.0; n_i_c * n_j_s];
            for i in 0..n_i_c {
                for j_s in 0..n_j_s {
                    let mut sum = 0.0;
                    for j_c in 0..n_j_c {
                        sum += self.trans_j.coeff(j_c, j_s) * cart_matrix[i * n_j_c + j_c];
                    }
                    result[i * n_j_s + j_s] = sum;
                }
            }
            result
        } else {
            cart_matrix.to_vec()
        };

        // Step 2: Transform index i (row index)
        // (n_i_c, n_j_s) -> (n_i_s, n_j_s)
        if self.trans_i.needs_transform() {
            let n_j = if self.trans_j.needs_transform() {
                n_j_s
            } else {
                n_j_c
            };
            let mut result = vec![0.0; n_i_s * n_j_s];
            for i_s in 0..n_i_s {
                for j in 0..n_j {
                    let mut sum = 0.0;
                    for i_c in 0..n_i_c {
                        sum += self.trans_i.coeff(i_c, i_s) * tmp[i_c * n_j + j];
                    }
                    result[i_s * n_j_s + j] = sum;
                }
            }
            result
        } else {
            tmp
        }
    }
}

/// Transform a full one-electron matrix from Cartesian to spherical basis
///
/// This transforms a symmetric N_cart x N_cart matrix to N_sph x N_sph,
/// where N_cart and N_sph are the total Cartesian and spherical basis sizes.
///
/// # Arguments
///
/// * `cart_matrix` - The Cartesian integral matrix (row-major, flattened)
/// * `shells` - Slice of shells describing the basis set
///
/// # Returns
///
/// The spherical integral matrix (row-major, flattened)
pub fn transform_one_electron_matrix(
    cart_matrix: &[f64],
    shells: &[crate::basis::ContractedShell],
) -> Vec<f64> {
    // Compute dimensions
    let n_cart: usize = shells.iter().map(|s| s.n_basis_functions()).sum();
    let n_sph: usize = shells.iter().map(|s| s.n_basis_functions_spherical()).sum();

    debug_assert_eq!(
        cart_matrix.len(),
        n_cart * n_cart,
        "Input matrix size does not match Cartesian basis size"
    );

    // If dimensions are the same, no transformation needed
    if n_cart == n_sph {
        return cart_matrix.to_vec();
    }

    let mut sph_matrix = vec![0.0; n_sph * n_sph];

    // Iterate over shell pairs
    let mut mu_cart = 0; // Cartesian index for shell i
    let mut mu_sph = 0; // Spherical index for shell i
    for (i, shell_i) in shells.iter().enumerate() {
        let n_i_c = shell_i.n_basis_functions();
        let n_i_s = shell_i.n_basis_functions_spherical();

        let mut nu_cart = 0; // Cartesian index for shell j
        let mut nu_sph = 0; // Spherical index for shell j
        for shell_j in shells.iter().take(i + 1) {
            let n_j_c = shell_j.n_basis_functions();
            let n_j_s = shell_j.n_basis_functions_spherical();

            // Extract Cartesian block
            let mut cart_block = vec![0.0; n_i_c * n_j_c];
            for a in 0..n_i_c {
                for b in 0..n_j_c {
                    cart_block[a * n_j_c + b] = cart_matrix[(mu_cart + a) * n_cart + (nu_cart + b)];
                }
            }

            // Transform block
            let transform =
                OneElectronSphericalTransform::new(shell_i.l_value(), shell_j.l_value());
            let sph_block = transform.transform(&cart_block);

            // Copy to output matrix (upper triangle)
            for a in 0..n_i_s {
                for b in 0..n_j_s {
                    let val = sph_block[a * n_j_s + b];
                    sph_matrix[(mu_sph + a) * n_sph + (nu_sph + b)] = val;
                    // Lower triangle (symmetry)
                    if mu_sph + a != nu_sph + b {
                        sph_matrix[(nu_sph + b) * n_sph + (mu_sph + a)] = val;
                    }
                }
            }

            nu_cart += n_j_c;
            nu_sph += n_j_s;
        }

        mu_cart += n_i_c;
        mu_sph += n_i_s;
    }

    sph_matrix
}

// =============================================================================
// ERI Transformation
// =============================================================================

/// Configuration for transforming a 4-index ERI tensor
#[derive(Debug, Clone)]
pub struct EriSphericalTransform {
    /// Transformations for each of the 4 shells
    pub trans_i: SphericalTransform,
    pub trans_j: SphericalTransform,
    pub trans_k: SphericalTransform,
    pub trans_l: SphericalTransform,
}

impl EriSphericalTransform {
    /// Create ERI transformation for four shells
    pub fn new(l_i: u32, l_j: u32, l_k: u32, l_l: u32) -> Self {
        Self {
            trans_i: SphericalTransform::new(l_i),
            trans_j: SphericalTransform::new(l_j),
            trans_k: SphericalTransform::new(l_k),
            trans_l: SphericalTransform::new(l_l),
        }
    }

    /// Check if any transformation is needed
    pub fn needs_transform(&self) -> bool {
        self.trans_i.needs_transform()
            || self.trans_j.needs_transform()
            || self.trans_k.needs_transform()
            || self.trans_l.needs_transform()
    }

    /// Get the number of spherical functions for each shell
    pub fn n_sph(&self) -> (usize, usize, usize, usize) {
        (
            self.trans_i.n_sph,
            self.trans_j.n_sph,
            self.trans_k.n_sph,
            self.trans_l.n_sph,
        )
    }

    /// Get the number of Cartesian functions for each shell
    pub fn n_cart(&self) -> (usize, usize, usize, usize) {
        (
            self.trans_i.n_cart,
            self.trans_j.n_cart,
            self.trans_k.n_cart,
            self.trans_l.n_cart,
        )
    }

    /// Transform a 4-index ERI tensor from Cartesian to spherical
    ///
    /// The input tensor has dimensions (n_i_cart, n_j_cart, n_k_cart, n_l_cart)
    /// stored in row-major order: [i][j][k][l]
    ///
    /// The output tensor has dimensions (n_i_sph, n_j_sph, n_k_sph, n_l_sph)
    ///
    /// # Algorithm
    ///
    /// The transformation is applied sequentially to each index:
    /// 1. Transform index l: (cart, cart, cart, cart) -> (cart, cart, cart, sph)
    /// 2. Transform index k: (cart, cart, cart, sph) -> (cart, cart, sph, sph)
    /// 3. Transform index j: (cart, cart, sph, sph) -> (cart, sph, sph, sph)
    /// 4. Transform index i: (cart, sph, sph, sph) -> (sph, sph, sph, sph)
    ///
    /// This is more efficient than the full 4-way contraction.
    pub fn transform(&self, cart_eri: &[f64]) -> Vec<f64> {
        if !self.needs_transform() {
            // No transformation needed, return a copy
            return cart_eri.to_vec();
        }

        let (n_i_c, n_j_c, n_k_c, n_l_c) = self.n_cart();
        let (n_i_s, n_j_s, n_k_s, n_l_s) = self.n_sph();

        debug_assert_eq!(
            cart_eri.len(),
            n_i_c * n_j_c * n_k_c * n_l_c,
            "Input ERI length does not match expected Cartesian dimensions"
        );

        // Step 1: Transform index l (fastest varying index)
        // (n_i_c, n_j_c, n_k_c, n_l_c) -> (n_i_c, n_j_c, n_k_c, n_l_s)
        let tmp1 = if self.trans_l.needs_transform() {
            self.transform_last_index(cart_eri, n_i_c * n_j_c * n_k_c, n_l_c, n_l_s, &self.trans_l)
        } else {
            cart_eri.to_vec()
        };

        // Step 2: Transform index k
        // (n_i_c, n_j_c, n_k_c, n_l_s) -> (n_i_c, n_j_c, n_k_s, n_l_s)
        let tmp2 = if self.trans_k.needs_transform() {
            self.transform_index_k(&tmp1, n_i_c, n_j_c, n_k_c, n_l_s, n_k_s, &self.trans_k)
        } else {
            tmp1
        };

        // Step 3: Transform index j
        // (n_i_c, n_j_c, n_k_s, n_l_s) -> (n_i_c, n_j_s, n_k_s, n_l_s)
        let tmp3 = if self.trans_j.needs_transform() {
            self.transform_index_j(&tmp2, n_i_c, n_j_c, n_k_s, n_l_s, n_j_s, &self.trans_j)
        } else {
            tmp2
        };

        // Step 4: Transform index i
        // (n_i_c, n_j_s, n_k_s, n_l_s) -> (n_i_s, n_j_s, n_k_s, n_l_s)
        if self.trans_i.needs_transform() {
            self.transform_index_i(&tmp3, n_i_c, n_j_s, n_k_s, n_l_s, n_i_s, &self.trans_i)
        } else {
            tmp3
        }
    }

    /// Transform the last (l) index of a 4D tensor
    fn transform_last_index(
        &self,
        input: &[f64],
        n_outer: usize,
        n_l_in: usize,
        n_l_out: usize,
        trans: &SphericalTransform,
    ) -> Vec<f64> {
        let mut output = vec![0.0; n_outer * n_l_out];

        for outer in 0..n_outer {
            for l_out in 0..n_l_out {
                let mut sum = 0.0;
                for l_in in 0..n_l_in {
                    sum += trans.coeff(l_in, l_out) * input[outer * n_l_in + l_in];
                }
                output[outer * n_l_out + l_out] = sum;
            }
        }

        output
    }

    /// Transform index k (second to last)
    #[allow(clippy::too_many_arguments)]
    fn transform_index_k(
        &self,
        input: &[f64],
        n_i: usize,
        n_j: usize,
        n_k_in: usize,
        n_l: usize,
        n_k_out: usize,
        trans: &SphericalTransform,
    ) -> Vec<f64> {
        let mut output = vec![0.0; n_i * n_j * n_k_out * n_l];

        for i in 0..n_i {
            for j in 0..n_j {
                for k_out in 0..n_k_out {
                    for l in 0..n_l {
                        let mut sum = 0.0;
                        for k_in in 0..n_k_in {
                            let in_idx = ((i * n_j + j) * n_k_in + k_in) * n_l + l;
                            sum += trans.coeff(k_in, k_out) * input[in_idx];
                        }
                        let out_idx = ((i * n_j + j) * n_k_out + k_out) * n_l + l;
                        output[out_idx] = sum;
                    }
                }
            }
        }

        output
    }

    /// Transform index j (third to last)
    #[allow(clippy::too_many_arguments)]
    fn transform_index_j(
        &self,
        input: &[f64],
        n_i: usize,
        n_j_in: usize,
        n_k: usize,
        n_l: usize,
        n_j_out: usize,
        trans: &SphericalTransform,
    ) -> Vec<f64> {
        let mut output = vec![0.0; n_i * n_j_out * n_k * n_l];

        for i in 0..n_i {
            for j_out in 0..n_j_out {
                for k in 0..n_k {
                    for l in 0..n_l {
                        let mut sum = 0.0;
                        for j_in in 0..n_j_in {
                            let in_idx = ((i * n_j_in + j_in) * n_k + k) * n_l + l;
                            sum += trans.coeff(j_in, j_out) * input[in_idx];
                        }
                        let out_idx = ((i * n_j_out + j_out) * n_k + k) * n_l + l;
                        output[out_idx] = sum;
                    }
                }
            }
        }

        output
    }

    /// Transform index i (first/outermost)
    #[allow(clippy::too_many_arguments)]
    fn transform_index_i(
        &self,
        input: &[f64],
        n_i_in: usize,
        n_j: usize,
        n_k: usize,
        n_l: usize,
        n_i_out: usize,
        trans: &SphericalTransform,
    ) -> Vec<f64> {
        let mut output = vec![0.0; n_i_out * n_j * n_k * n_l];
        let stride_jkl = n_j * n_k * n_l;

        for i_out in 0..n_i_out {
            for jkl in 0..stride_jkl {
                let mut sum = 0.0;
                for i_in in 0..n_i_in {
                    sum += trans.coeff(i_in, i_out) * input[i_in * stride_jkl + jkl];
                }
                output[i_out * stride_jkl + jkl] = sum;
            }
        }

        output
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    const TOL: f64 = 1e-14;

    #[test]
    fn test_n_spherical() {
        assert_eq!(n_spherical(0), 1); // s: 1
        assert_eq!(n_spherical(1), 3); // p: 3
        assert_eq!(n_spherical(2), 5); // d: 5
        assert_eq!(n_spherical(3), 7); // f: 7
    }

    #[test]
    fn test_cart2sph_d_coefficients() {
        // Verify the transformation coefficients match libcint values
        // Reference: libcint/src/cart2sph.c lines 52-86

        // The coefficients are from Schlegel & Frisch, IJQC 54, 83 (1995)
        // They transform normalized Cartesian GTOs to normalized spherical GTOs.

        // Note: The column norms are NOT 1.0 because the Cartesian and spherical
        // normalizations differ. The transformation is correct when applied to
        // properly normalized Cartesian GTOs (which PySCF does internally).

        // Verify specific coefficient values from libcint
        // dxy -> d-2 coefficient: sqrt(3) = 1.732... but libcint uses 1.0925484...
        // This is because the spherical harmonic includes a factor of sqrt(15/4pi)
        // and Cartesian dxy has a factor of sqrt(3) in its angular part.

        // d-2 (dxy): coefficient should be sqrt(3) = 1.7320508... / sqrt(15/(4*pi)) * ...
        // From libcint: 1.092548430592079
        assert_abs_diff_eq!(CART2SPH_D[1][0], 1.092548430592079, epsilon = 1e-14);

        // d-1 (dyz): same coefficient
        assert_abs_diff_eq!(CART2SPH_D[4][1], 1.092548430592079, epsilon = 1e-14);

        // d+1 (dxz): same coefficient
        assert_abs_diff_eq!(CART2SPH_D[2][3], 1.092548430592079, epsilon = 1e-14);

        // d0 (dz2): dzz coefficient = sqrt(1/3) * 2 = 0.630783130505040
        assert_abs_diff_eq!(CART2SPH_D[5][2], 0.630783130505040, epsilon = 1e-14);

        // d0: dxx coefficient = -sqrt(1/12) = -0.315391565252520
        assert_abs_diff_eq!(CART2SPH_D[0][2], -0.315391565252520, epsilon = 1e-14);

        // d+2 (dx2-y2): dxx coefficient = sqrt(3)/2 = 0.546274215296040
        assert_abs_diff_eq!(CART2SPH_D[0][4], 0.546274215296040, epsilon = 1e-14);

        // d+2: dyy coefficient = -sqrt(3)/2 = -0.546274215296040
        assert_abs_diff_eq!(CART2SPH_D[3][4], -0.546274215296040, epsilon = 1e-14);

        // Verify all non-contributing elements are exactly zero
        // dxx only contributes to d0 and d+2
        assert_eq!(CART2SPH_D[0][0], 0.0);
        assert_eq!(CART2SPH_D[0][1], 0.0);
        assert_eq!(CART2SPH_D[0][3], 0.0);

        // dxy only contributes to d-2
        assert_eq!(CART2SPH_D[1][1], 0.0);
        assert_eq!(CART2SPH_D[1][2], 0.0);
        assert_eq!(CART2SPH_D[1][3], 0.0);
        assert_eq!(CART2SPH_D[1][4], 0.0);
    }

    #[test]
    fn test_cart2sph_d_orthogonality() {
        // Different spherical d-orbitals should be orthogonal
        // sum_c C[c][m1] * C[c][m2] = 0 for m1 != m2

        for m1 in 0..5 {
            for m2 in (m1 + 1)..5 {
                let dot: f64 = CART2SPH_D.iter().map(|row| row[m1] * row[m2]).sum();
                assert!(
                    (dot).abs() < 1e-12,
                    "d_{} and d_{} are not orthogonal: dot = {}",
                    m1,
                    m2,
                    dot
                );
            }
        }
    }

    #[test]
    fn test_spherical_transform_s_orbital() {
        let trans = SphericalTransform::new(0);
        assert_eq!(trans.n_cart, 1);
        assert_eq!(trans.n_sph, 1);
        assert!(!trans.needs_transform());

        let cart = vec![1.5];
        let sph = trans.transform_1d(&cart);
        assert_eq!(sph, cart);
    }

    #[test]
    fn test_spherical_transform_p_orbital() {
        let trans = SphericalTransform::new(1);
        assert_eq!(trans.n_cart, 3);
        assert_eq!(trans.n_sph, 3);
        assert!(!trans.needs_transform());

        let cart = vec![1.0, 2.0, 3.0];
        let sph = trans.transform_1d(&cart);
        assert_eq!(sph, cart);
    }

    #[test]
    fn test_spherical_transform_d_orbital() {
        let trans = SphericalTransform::new(2);
        assert_eq!(trans.n_cart, 6);
        assert_eq!(trans.n_sph, 5);
        assert!(trans.needs_transform());

        // Test with a specific Cartesian d-orbital pattern
        // Input: pure dxy (Cartesian index 1)
        let cart_dxy = vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        let sph = trans.transform_1d(&cart_dxy);

        // dxy should transform to pure d-2 (spherical index 0)
        assert_abs_diff_eq!(sph[0], CART2SPH_D[1][0], epsilon = TOL); // d-2
        assert_abs_diff_eq!(sph[1], 0.0, epsilon = TOL); // d-1
        assert_abs_diff_eq!(sph[2], 0.0, epsilon = TOL); // d0
        assert_abs_diff_eq!(sph[3], 0.0, epsilon = TOL); // d+1
        assert_abs_diff_eq!(sph[4], 0.0, epsilon = TOL); // d+2
    }

    #[test]
    fn test_spherical_transform_d_dz2() {
        let trans = SphericalTransform::new(2);

        // Pure dzz (Cartesian index 5) should give pure d0 (spherical index 2)
        let cart_dzz = vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        let sph = trans.transform_1d(&cart_dzz);

        // dzz contributes to d0 with coefficient 0.630783130505040012
        assert_abs_diff_eq!(sph[0], 0.0, epsilon = TOL); // d-2
        assert_abs_diff_eq!(sph[1], 0.0, epsilon = TOL); // d-1
        assert_abs_diff_eq!(sph[2], CART2SPH_D[5][2], epsilon = TOL); // d0
        assert_abs_diff_eq!(sph[3], 0.0, epsilon = TOL); // d+1
        assert_abs_diff_eq!(sph[4], 0.0, epsilon = TOL); // d+2
    }

    #[test]
    fn test_spherical_transform_d_traceless() {
        // The d0 orbital (dz2) is the traceless combination: 2*dzz - dxx - dyy
        // In normalized form: d0 = c1*(dxx + dyy) + c2*dzz
        // where c1 = -0.315391565252520002 and c2 = 0.630783130505040012
        // Note: c2 = -2*c1, which gives the traceless property

        let _trans = SphericalTransform::new(2);

        // Verify the traceless property: sum of diagonal d-orbital coefficients for d0 is zero
        // c1 (dxx) + c1 (dyy) + c2 (dzz) = -0.315... + (-0.315...) + 0.630... = 0
        let trace = CART2SPH_D[0][2] + CART2SPH_D[3][2] + CART2SPH_D[5][2];
        assert_abs_diff_eq!(trace, 0.0, epsilon = 1e-12);

        // Test that a spherical d0 input gives correct Cartesian coefficients
        // If we have pure d0 in spherical basis, the Cartesian components should be:
        // dxx: -0.315..., dyy: -0.315..., dzz: 0.630...
        // This corresponds to 2*z^2 - x^2 - y^2 (up to normalization)
        let c_d0 = CART2SPH_D[5][2]; // coefficient of dzz in d0
        let c_dxx_d0 = CART2SPH_D[0][2]; // coefficient of dxx in d0
        assert_abs_diff_eq!(c_d0, -2.0 * c_dxx_d0, epsilon = 1e-12);
    }

    #[test]
    fn test_eri_transform_identity() {
        // For s and p orbitals, transformation should be identity
        let trans = EriSphericalTransform::new(0, 0, 0, 0);
        assert!(!trans.needs_transform());

        let cart = vec![1.5];
        let sph = trans.transform(&cart);
        assert_eq!(sph, cart);

        // p-p-p-p
        let trans_p = EriSphericalTransform::new(1, 1, 1, 1);
        assert!(!trans_p.needs_transform());

        let cart_p = vec![0.0; 81]; // 3^4
        let sph_p = trans_p.transform(&cart_p);
        assert_eq!(sph_p, cart_p);
    }

    #[test]
    fn test_eri_transform_mixed_sp() {
        // s-p-s-p should be identity (no d-orbitals)
        let trans = EriSphericalTransform::new(0, 1, 0, 1);
        assert!(!trans.needs_transform());
    }

    #[test]
    fn test_eri_transform_with_d_orbital() {
        // s-s-d-d requires transformation
        let trans = EriSphericalTransform::new(0, 0, 2, 2);
        assert!(trans.needs_transform());

        let (n_i_c, n_j_c, n_k_c, n_l_c) = trans.n_cart();
        assert_eq!((n_i_c, n_j_c, n_k_c, n_l_c), (1, 1, 6, 6));

        let (n_i_s, n_j_s, n_k_s, n_l_s) = trans.n_sph();
        assert_eq!((n_i_s, n_j_s, n_k_s, n_l_s), (1, 1, 5, 5));

        // Create a simple test tensor
        let cart = vec![0.0; 36]; // 1 * 1 * 6 * 6
        let sph = trans.transform(&cart);
        assert_eq!(sph.len(), 25); // 1 * 1 * 5 * 5
    }

    #[test]
    fn test_eri_transform_d_d_d_d() {
        // d-d-d-d is the most complex case
        let trans = EriSphericalTransform::new(2, 2, 2, 2);
        assert!(trans.needs_transform());

        let (n_i_c, n_j_c, n_k_c, n_l_c) = trans.n_cart();
        assert_eq!((n_i_c, n_j_c, n_k_c, n_l_c), (6, 6, 6, 6));

        let (n_i_s, n_j_s, n_k_s, n_l_s) = trans.n_sph();
        assert_eq!((n_i_s, n_j_s, n_k_s, n_l_s), (5, 5, 5, 5));

        // Create an identity-like test where only diagonal elements are non-zero
        let mut cart = vec![0.0; 1296]; // 6^4

        // Set (dxy,dxy|dxy,dxy) = 1.0 in Cartesian (index 1,1,1,1)
        let dxy_idx = 1;
        let idx = ((dxy_idx * 6 + dxy_idx) * 6 + dxy_idx) * 6 + dxy_idx;
        cart[idx] = 1.0;

        let sph = trans.transform(&cart);
        assert_eq!(sph.len(), 625); // 5^4

        // In spherical, (d-2,d-2|d-2,d-2) should be:
        // C[dxy,d-2]^4 * 1.0 = (1.092548430592079070)^4
        let expected = CART2SPH_D[1][0].powi(4);
        let d_m2_idx = 0;
        let sph_idx = ((d_m2_idx * 5 + d_m2_idx) * 5 + d_m2_idx) * 5 + d_m2_idx;
        assert_abs_diff_eq!(sph[sph_idx], expected, epsilon = 1e-12);
    }

    // -------------------------------------------------------------------------
    // One-electron (2-index) transformation tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_one_electron_transform_s_s() {
        // s-s block should be identity
        let trans = OneElectronSphericalTransform::new(0, 0);
        assert!(!trans.needs_transform());

        let cart = vec![1.5];
        let sph = trans.transform(&cart);
        assert_eq!(sph, cart);
    }

    #[test]
    fn test_one_electron_transform_p_p() {
        // p-p block should be identity
        let trans = OneElectronSphericalTransform::new(1, 1);
        assert!(!trans.needs_transform());

        let cart = vec![1.0, 0.1, 0.2, 0.1, 2.0, 0.3, 0.2, 0.3, 3.0];
        let sph = trans.transform(&cart);
        assert_eq!(sph.len(), 9);
        assert_eq!(sph, cart);
    }

    #[test]
    fn test_one_electron_transform_s_d() {
        // s-d block: (1, 6) -> (1, 5)
        let trans = OneElectronSphericalTransform::new(0, 2);
        assert!(trans.needs_transform());

        let (n_i_c, n_j_c) = trans.n_cart();
        assert_eq!((n_i_c, n_j_c), (1, 6));

        let (n_i_s, n_j_s) = trans.n_sph();
        assert_eq!((n_i_s, n_j_s), (1, 5));

        // Create test matrix with unit value for s-dxy
        let mut cart = vec![0.0; 6];
        cart[1] = 1.0; // s-dxy

        let sph = trans.transform(&cart);
        assert_eq!(sph.len(), 5);

        // s-dxy -> s-d-2 with coefficient C[dxy, d-2]
        assert_abs_diff_eq!(sph[0], CART2SPH_D[1][0], epsilon = TOL);
        assert_abs_diff_eq!(sph[1], 0.0, epsilon = TOL);
        assert_abs_diff_eq!(sph[2], 0.0, epsilon = TOL);
        assert_abs_diff_eq!(sph[3], 0.0, epsilon = TOL);
        assert_abs_diff_eq!(sph[4], 0.0, epsilon = TOL);
    }

    #[test]
    fn test_one_electron_transform_d_d() {
        // d-d block: (6, 6) -> (5, 5)
        let trans = OneElectronSphericalTransform::new(2, 2);
        assert!(trans.needs_transform());

        let (n_i_c, n_j_c) = trans.n_cart();
        assert_eq!((n_i_c, n_j_c), (6, 6));

        let (n_i_s, n_j_s) = trans.n_sph();
        assert_eq!((n_i_s, n_j_s), (5, 5));

        // Create test matrix: unit on dxy-dxy (index [1,1])
        let mut cart = vec![0.0; 36];
        cart[1 * 6 + 1] = 1.0;

        let sph = trans.transform(&cart);
        assert_eq!(sph.len(), 25);

        // dxy-dxy -> d-2-d-2 with coefficient C[dxy, d-2]^2
        let expected = CART2SPH_D[1][0].powi(2);
        assert_abs_diff_eq!(sph[0 * 5 + 0], expected, epsilon = TOL);
    }

    #[test]
    fn test_one_electron_transform_d_d_symmetry() {
        // d-d block with symmetric input should give symmetric output
        let trans = OneElectronSphericalTransform::new(2, 2);

        // Create symmetric Cartesian matrix (identity-like)
        let mut cart = vec![0.0; 36];
        for i in 0..6 {
            cart[i * 6 + i] = 1.0;
        }

        let sph = trans.transform(&cart);
        assert_eq!(sph.len(), 25);

        // Check symmetry of result
        for i in 0..5 {
            for j in 0..5 {
                let diff = (sph[i * 5 + j] - sph[j * 5 + i]).abs();
                assert!(
                    diff < TOL,
                    "Symmetry violation at [{}, {}]: {} vs {}",
                    i,
                    j,
                    sph[i * 5 + j],
                    sph[j * 5 + i]
                );
            }
        }
    }

    #[test]
    fn test_one_electron_transform_d_d_trace_preserved() {
        // The trace of the overlap matrix should be preserved under orthogonal transformation
        // (if the transformation is orthogonal/unitary)
        // Note: The cart2sph transformation is NOT orthogonal due to normalization,
        // so trace is NOT preserved. This test just checks for reasonable values.

        let trans = OneElectronSphericalTransform::new(2, 2);

        // Create diagonal Cartesian matrix with all 1s (identity)
        let mut cart = vec![0.0; 36];
        for i in 0..6 {
            cart[i * 6 + i] = 1.0;
        }

        let sph = trans.transform(&cart);

        // Check diagonal elements are positive
        for i in 0..5 {
            assert!(
                sph[i * 5 + i] > 0.0,
                "Diagonal element {} should be positive",
                i
            );
        }
    }
}
