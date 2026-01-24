//! Gaussian Product Theorem implementation
//!
//! When two Gaussian functions centered at A and B are multiplied, the result
//! is proportional to a Gaussian centered at P (the "product center"):
//!
//! ```text
//! exp(-α|r-A|²) · exp(-β|r-B|²) = K_AB · exp(-p|r-P|²)
//! ```
//!
//! where:
//! - `p = α + β` (combined exponent)
//! - `P = (αA + βB) / p` (product center)
//! - `μ = αβ / p` (reduced exponent)
//! - `K_AB = exp(-μ|A-B|²)` (overlap prefactor)
//!
//! # Reference
//!
//! - Helgaker et al. (2000), Chapter 9.2
//! - libcint: `references/libcint/src/g1e.c`, line 125-184

use std::f64::consts::PI;

/// Pre-computed data for a pair of Gaussian primitives
///
/// Contains all quantities needed for Obara-Saika recurrence relations.
#[derive(Debug, Clone)]
pub struct GaussianProduct {
    /// Combined exponent: p = α + β
    pub p: f64,

    /// Reduced exponent: μ = αβ / p
    pub mu: f64,

    /// Product center: P = (αA + βB) / p
    pub center_p: [f64; 3],

    /// Vector from center A to product center P: PA = P - A
    pub pa: [f64; 3],

    /// Vector from center B to product center P: PB = P - B
    pub pb: [f64; 3],

    /// Vector from A to B: AB = B - A
    pub ab: [f64; 3],

    /// Squared distance |A - B|²
    pub ab_squared: f64,

    /// Overlap prefactor: K_AB = exp(-μ|A-B|²)
    pub k_ab: f64,

    /// Base integral coefficient: (π/p)^(3/2) · K_AB
    /// This is [s|s] = ⟨0|0⟩ for the primitive pair
    pub ss_integral: f64,

    /// Coefficient 1/(2p) used in recurrence relations
    pub one_over_2p: f64,
}

impl GaussianProduct {
    /// Create a GaussianProduct from two primitive Gaussians
    ///
    /// # Arguments
    ///
    /// * `alpha` - Exponent of first Gaussian
    /// * `center_a` - Center of first Gaussian [x, y, z]
    /// * `beta` - Exponent of second Gaussian
    /// * `center_b` - Center of second Gaussian [x, y, z]
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::integrals::GaussianProduct;
    ///
    /// let a = [0.0, 0.0, 0.0];
    /// let b = [0.0, 0.0, 1.4];
    /// let gp = GaussianProduct::new(3.425250914, &a, 3.425250914, &b);
    ///
    /// // Product center is midway for equal exponents
    /// assert!((gp.center_p[2] - 0.7).abs() < 1e-10);
    /// ```
    pub fn new(alpha: f64, center_a: &[f64; 3], beta: f64, center_b: &[f64; 3]) -> Self {
        let p = alpha + beta;
        let mu = alpha * beta / p;
        let one_over_2p = 0.5 / p;

        // Product center P = (αA + βB) / p
        let center_p = [
            (alpha * center_a[0] + beta * center_b[0]) / p,
            (alpha * center_a[1] + beta * center_b[1]) / p,
            (alpha * center_a[2] + beta * center_b[2]) / p,
        ];

        // PA = P - A
        let pa = [
            center_p[0] - center_a[0],
            center_p[1] - center_a[1],
            center_p[2] - center_a[2],
        ];

        // PB = P - B
        let pb = [
            center_p[0] - center_b[0],
            center_p[1] - center_b[1],
            center_p[2] - center_b[2],
        ];

        // AB = B - A
        let ab = [
            center_b[0] - center_a[0],
            center_b[1] - center_a[1],
            center_b[2] - center_a[2],
        ];

        // |A - B|²
        let ab_squared = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];

        // K_AB = exp(-μ|A-B|²)
        let k_ab = (-mu * ab_squared).exp();

        // Base integral: (π/p)^(3/2) · K_AB
        let pi_over_p = PI / p;
        let ss_integral = pi_over_p.sqrt() * pi_over_p * k_ab;

        Self {
            p,
            mu,
            center_p,
            pa,
            pb,
            ab,
            ab_squared,
            k_ab,
            ss_integral,
            one_over_2p,
        }
    }

    /// Check if the primitives are on the same center (within tolerance)
    #[inline]
    pub fn is_same_center(&self) -> bool {
        self.ab_squared < 1e-20
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    #[test]
    fn gaussian_product_same_center() {
        let center = [0.0, 0.0, 0.0];
        let gp = GaussianProduct::new(1.0, &center, 2.0, &center);

        assert!((gp.p - 3.0).abs() < TOL);
        assert!((gp.mu - 2.0 / 3.0).abs() < TOL);
        assert!((gp.k_ab - 1.0).abs() < TOL);
        assert!(gp.is_same_center());

        // Product center is at origin
        for i in 0..3 {
            assert!(gp.center_p[i].abs() < TOL);
        }
    }

    #[test]
    fn gaussian_product_equal_exponents() {
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 2.0];
        let alpha = 1.0;

        let gp = GaussianProduct::new(alpha, &a, alpha, &b);

        // With equal exponents, P is the midpoint
        assert!((gp.center_p[2] - 1.0).abs() < TOL);

        // p = 2α
        assert!((gp.p - 2.0).abs() < TOL);

        // μ = α²/(2α) = α/2
        assert!((gp.mu - 0.5).abs() < TOL);

        // K_AB = exp(-0.5 * 4) = exp(-2)
        let expected_k = (-2.0_f64).exp();
        assert!((gp.k_ab - expected_k).abs() < TOL);
    }

    #[test]
    fn gaussian_product_unequal_exponents() {
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 1.0];
        let alpha = 1.0;
        let beta = 3.0;

        let gp = GaussianProduct::new(alpha, &a, beta, &b);

        // p = 4
        assert!((gp.p - 4.0).abs() < TOL);

        // Product center: P_z = (1*0 + 3*1)/4 = 0.75
        assert!((gp.center_p[2] - 0.75).abs() < TOL);

        // PA_z = 0.75 - 0 = 0.75
        assert!((gp.pa[2] - 0.75).abs() < TOL);

        // PB_z = 0.75 - 1 = -0.25
        assert!((gp.pb[2] - (-0.25)).abs() < TOL);
    }

    #[test]
    fn gaussian_product_ss_integral_same_center() {
        let center = [0.0, 0.0, 0.0];
        let alpha = 1.0;
        let gp = GaussianProduct::new(alpha, &center, alpha, &center);

        // [s|s] = (π/p)^(3/2) for same center
        let p = 2.0 * alpha;
        let expected = (PI / p).powf(1.5);
        assert!((gp.ss_integral - expected).abs() < TOL);
    }

    #[test]
    fn gaussian_product_ab_vector() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];

        let gp = GaussianProduct::new(1.0, &a, 1.0, &b);

        // AB = B - A
        assert!((gp.ab[0] - 3.0).abs() < TOL);
        assert!((gp.ab[1] - 3.0).abs() < TOL);
        assert!((gp.ab[2] - 3.0).abs() < TOL);

        // |AB|² = 27
        assert!((gp.ab_squared - 27.0).abs() < TOL);
    }
}
