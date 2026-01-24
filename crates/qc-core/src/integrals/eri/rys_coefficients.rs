//! Rys quadrature recurrence coefficients
//!
//! After obtaining Rys roots u_n, we compute the recurrence coefficients
//! needed for the 2D vertical recurrence relations (VRR).
//!
//! # Coefficient Formulas
//!
//! For a given Rys root u (where u is in [0, 1)):
//!
//! ```text
//! a0 = p * q / (p + q)              // Same as rho
//! u2 = a0 * u                       // Scaled root
//! tmp4 = 0.5 / (u2 * (p + q) + p * q)
//! tmp5 = u2 * tmp4
//!
//! b00 = tmp5                        // Cross-term coupling
//! b10 = tmp5 + tmp4 * q             // VRR coefficient for bra
//! b01 = tmp5 + tmp4 * p             // VRR coefficient for ket
//!
//! c00[i] = PA[i] - 2 * tmp5 * q * PQ[i]    // Center offset for bra
//! c0p[i] = QC[i] + 2 * tmp5 * p * PQ[i]    // Center offset for ket
//! ```
//!
//! # Reference
//!
//! - libcint g2e.c lines 4542-4564 (CINTg0_2e coefficient computation)
//! - Dupuis, Rys & King (1976), J. Chem. Phys. 65, 111

use super::GaussianProduct2e;

/// Rys quadrature recurrence coefficients for a single root
///
/// These coefficients are used in the 2D VRR to build integrals g[n][m].
#[derive(Debug, Clone, Copy)]
pub struct RysCoefficients {
    /// Cross-term coupling coefficient: b00
    ///
    /// Used in both VRR directions when both n > 0 and m > 0.
    pub b00: f64,

    /// Bra VRR coefficient: b10
    ///
    /// Coefficient of [n-1, m] term in the bra recurrence.
    pub b10: f64,

    /// Ket VRR coefficient: b01
    ///
    /// Coefficient of [n, m-1] term in the ket recurrence.
    pub b01: f64,

    /// Center offset for bra (x, y, z components): c00
    ///
    /// The linear term coefficient in the bra recurrence:
    /// [n+1, m] = c00 * [n, m] + ...
    pub c00: [f64; 3],

    /// Center offset for ket (x, y, z components): c0p
    ///
    /// The linear term coefficient in the ket recurrence:
    /// [n, m+1] = c0p * [n, m] + ...
    pub c0p: [f64; 3],
}

impl RysCoefficients {
    /// Compute Rys coefficients for a given root
    ///
    /// # Arguments
    ///
    /// * `gp2e` - Pre-computed two-electron Gaussian product data
    /// * `u` - The Rys quadrature root (in [0, 1))
    ///
    /// # Returns
    ///
    /// RysCoefficients containing b00, b10, b01, c00[3], c0p[3]
    ///
    /// # Reference
    ///
    /// libcint g2e.c lines 4542-4564
    ///
    /// ```text
    /// u2 = a0 * u;
    /// tmp4 = .5 / (u2 * (aij + akl) + a1);
    /// tmp5 = u2 * tmp4;
    /// b00 = tmp5;
    /// b10 = tmp5 + tmp4 * akl;
    /// b01 = tmp5 + tmp4 * aij;
    /// c00x = rijrx[0] - tmp2 * rijrkl[0];
    /// c00y = rijrx[1] - tmp2 * rijrkl[1];
    /// c00z = rijrx[2] - tmp2 * rijrkl[2];
    /// c0px = rklrx[0] + tmp3 * rijrkl[0];
    /// c0py = rklrx[1] + tmp3 * rijrkl[1];
    /// c0pz = rklrx[2] + tmp3 * rijrkl[2];
    /// ```
    pub fn compute(gp2e: &GaussianProduct2e, u: f64) -> Self {
        let p = gp2e.p;
        let q = gp2e.q;

        // a0 = p * q / (p + q) = rho
        let a0 = gp2e.rho;

        // a1 = p * q
        let a1 = p * q;

        // CRITICAL: Transform root from IQCP convention to libcint convention
        //
        // IQCP rys_roots returns: u = t^2 where t is the actual quadrature point in (0, 1)
        // libcint formula expects: u_libcint = t^2 / (1 - t^2)
        //
        // Transformation: u_libcint = u / (1 - u)
        //
        // Reference: libcint g2e.c lines 4543-4546:
        //   "u(irys) = t2/(1-t2)
        //    t2 = u(irys)/(1+u(irys))"
        //
        // For numerical stability, handle u approaching 1 (root near boundary)
        let u_transformed = if u < 1.0 - 1e-15 {
            u / (1.0 - u)
        } else {
            // When u is very close to 1, the root is at the boundary
            // This can happen for very large T; return a large value
            1e15
        };

        // u2 = a0 * u_libcint
        let u2 = a0 * u_transformed;

        // tmp4 = 0.5 / (u2 * (p + q) + p * q)
        let tmp4 = 0.5 / (u2 * (p + q) + a1);

        // tmp5 = u2 * tmp4
        let tmp5 = u2 * tmp4;

        // b00 = tmp5 (cross-term coupling)
        let b00 = tmp5;

        // b10 = tmp5 + tmp4 * q (bra coefficient)
        let b10 = tmp5 + tmp4 * q;

        // b01 = tmp5 + tmp4 * p (ket coefficient)
        let b01 = tmp5 + tmp4 * p;

        // tmp2 = 2 * tmp5 * q (for c00)
        let tmp2 = 2.0 * tmp5 * q;

        // tmp3 = 2 * tmp5 * p (for c0p)
        let tmp3 = 2.0 * tmp5 * p;

        // c00[i] = PA[i] - tmp2 * PQ[i]
        // c0p[i] = QC[i] + tmp3 * PQ[i]
        let c00 = [
            gp2e.pa[0] - tmp2 * gp2e.pq[0],
            gp2e.pa[1] - tmp2 * gp2e.pq[1],
            gp2e.pa[2] - tmp2 * gp2e.pq[2],
        ];

        let c0p = [
            gp2e.qc[0] + tmp3 * gp2e.pq[0],
            gp2e.qc[1] + tmp3 * gp2e.pq[1],
            gp2e.qc[2] + tmp3 * gp2e.pq[2],
        ];

        Self {
            b00,
            b10,
            b01,
            c00,
            c0p,
        }
    }

    /// Compute Rys coefficients for T = 0 (special case)
    ///
    /// When T = 0, the problem simplifies and we use root = 0, weight = 1.
    /// This method computes the coefficients for this limiting case.
    ///
    /// # Arguments
    ///
    /// * `gp2e` - Pre-computed two-electron Gaussian product data
    ///
    /// # Returns
    ///
    /// RysCoefficients for the T = 0 case
    pub fn compute_t_zero(gp2e: &GaussianProduct2e) -> Self {
        let p = gp2e.p;
        let q = gp2e.q;

        // At u = 0:
        // u2 = 0
        // tmp4 = 0.5 / (0 + p * q) = 0.5 / (p * q)
        // tmp5 = 0
        let a1 = p * q;
        let tmp4 = 0.5 / a1;

        // b00 = 0
        let b00 = 0.0;

        // b10 = tmp4 * q = q / (2 * p * q) = 1 / (2 * p)
        let b10 = tmp4 * q;

        // b01 = tmp4 * p = p / (2 * p * q) = 1 / (2 * q)
        let b01 = tmp4 * p;

        // c00[i] = PA[i] (since tmp2 = 0)
        let c00 = gp2e.pa;

        // c0p[i] = QC[i] (since tmp3 = 0)
        let c0p = gp2e.qc;

        Self {
            b00,
            b10,
            b01,
            c00,
            c0p,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    #[test]
    fn test_coefficients_t_zero() {
        // At T = 0, u = 0, coefficients simplify
        let center = [0.0, 0.0, 0.0];
        let gp = GaussianProduct2e::new(1.0, &center, 1.0, &center, 1.0, &center, 1.0, &center);

        let coeffs = RysCoefficients::compute(&gp, 0.0);

        // b00 = 0 at u = 0
        assert!((coeffs.b00).abs() < TOL);

        // b10 = 1/(2p) = 1/4
        assert!((coeffs.b10 - 0.25).abs() < TOL);

        // b01 = 1/(2q) = 1/4
        assert!((coeffs.b01 - 0.25).abs() < TOL);

        // c00 = PA = 0 (same center)
        for i in 0..3 {
            assert!((coeffs.c00[i]).abs() < TOL);
            assert!((coeffs.c0p[i]).abs() < TOL);
        }
    }

    #[test]
    fn test_coefficients_nonzero_root() {
        let a = [0.0, 0.0, 0.0];
        let c = [0.0, 0.0, 2.0];
        let gp = GaussianProduct2e::new(1.0, &a, 1.0, &a, 1.0, &c, 1.0, &c);

        // Use a typical root value
        let u = 0.3;
        let coeffs = RysCoefficients::compute(&gp, u);

        // All coefficients should be finite
        assert!(coeffs.b00.is_finite());
        assert!(coeffs.b10.is_finite());
        assert!(coeffs.b01.is_finite());

        for i in 0..3 {
            assert!(coeffs.c00[i].is_finite());
            assert!(coeffs.c0p[i].is_finite());
        }

        // b10, b01 should be positive
        assert!(coeffs.b10 > 0.0);
        assert!(coeffs.b01 > 0.0);

        // b00 should be non-negative
        assert!(coeffs.b00 >= 0.0);
    }

    #[test]
    fn test_coefficients_symmetry() {
        // When p = q, b10 and b01 should be equal
        let center = [0.0, 0.0, 0.0];
        let gp = GaussianProduct2e::new(1.0, &center, 1.0, &center, 1.0, &center, 1.0, &center);

        let u = 0.5;
        let coeffs = RysCoefficients::compute(&gp, u);

        // p = q = 2, so b10 = b01
        assert!((coeffs.b10 - coeffs.b01).abs() < TOL);
    }

    #[test]
    fn test_compute_t_zero_same_as_u_zero() {
        // compute_t_zero should give same result as compute with u=0
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let gp = GaussianProduct2e::new(1.0, &a, 1.0, &b, 1.0, &a, 1.0, &b);

        let coeffs_u0 = RysCoefficients::compute(&gp, 0.0);
        let coeffs_t0 = RysCoefficients::compute_t_zero(&gp);

        assert!((coeffs_u0.b00 - coeffs_t0.b00).abs() < TOL);
        assert!((coeffs_u0.b10 - coeffs_t0.b10).abs() < TOL);
        assert!((coeffs_u0.b01 - coeffs_t0.b01).abs() < TOL);

        for i in 0..3 {
            assert!((coeffs_u0.c00[i] - coeffs_t0.c00[i]).abs() < TOL);
            assert!((coeffs_u0.c0p[i] - coeffs_t0.c0p[i]).abs() < TOL);
        }
    }

    #[test]
    fn test_coefficients_with_different_exponents() {
        let a = [0.0, 0.0, 0.0];
        let b = [0.5, 0.0, 0.0];
        let c = [0.0, 0.5, 0.0];
        let d = [0.0, 0.0, 0.5];

        // Different exponents
        let gp = GaussianProduct2e::new(1.0, &a, 2.0, &b, 1.5, &c, 0.5, &d);

        // u is the Rys root in IQCP convention (t^2 where t is quadrature point in [0,1])
        let u = 0.4;
        let coeffs = RysCoefficients::compute(&gp, u);

        // p = 3, q = 2
        // b10 should have q contribution, b01 should have p contribution
        // So with q < p, we expect b10 < b01
        assert!(coeffs.b10 < coeffs.b01);

        // Verify formula manually (with root transformation)
        // The transformation is: u_libcint = u / (1 - u)
        let p = gp.p;
        let q = gp.q;
        let a0 = gp.rho;
        let u_transformed = u / (1.0 - u); // Transform to libcint convention
        let u2 = a0 * u_transformed;
        let tmp4 = 0.5 / (u2 * (p + q) + p * q);
        let tmp5 = u2 * tmp4;

        let expected_b10 = tmp5 + tmp4 * q;
        let expected_b01 = tmp5 + tmp4 * p;

        assert!((coeffs.b10 - expected_b10).abs() < TOL);
        assert!((coeffs.b01 - expected_b01).abs() < TOL);
    }
}
