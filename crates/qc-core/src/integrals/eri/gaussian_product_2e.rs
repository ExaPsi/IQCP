//! Two-electron Gaussian product data structure
//!
//! Extends the Gaussian Product Theorem to four-center integrals by
//! combining the bra pair (A, B) -> P and ket pair (C, D) -> Q.
//!
//! # Gaussian Product Theorem (Two Pairs)
//!
//! For bra primitives at A and B with exponents alpha, beta:
//! ```text
//! p = alpha + beta
//! P = (alpha * A + beta * B) / p
//! K_ij = exp(-alpha * beta / p * |A - B|^2)
//! ```
//!
//! For ket primitives at C and D with exponents gamma, delta:
//! ```text
//! q = gamma + delta
//! Q = (gamma * C + delta * D) / q
//! K_kl = exp(-gamma * delta / q * |C - D|^2)
//! ```
//!
//! The combined quantities for the Rys quadrature:
//! ```text
//! rho = p * q / (p + q)
//! T = rho * |P - Q|^2
//! prefactor = 2 * pi^(5/2) / (p * q * sqrt(p + q)) * K_ij * K_kl
//! ```
//!
//! # Reference
//!
//! - Dupuis, Rys & King (1976), J. Chem. Phys. 65, 111
//! - libcint g2e.c lines 4425-4470

use std::f64::consts::PI;

/// Pre-computed data for a two-electron primitive quartet
///
/// Contains all quantities needed for Rys quadrature ERI evaluation.
#[derive(Debug, Clone)]
pub struct GaussianProduct2e {
    // === Bra pair (i, j) centered at A, B ===
    /// Combined exponent: p = alpha_i + alpha_j
    pub p: f64,

    /// Product center: P = (alpha_i * A + alpha_j * B) / p
    pub center_p: [f64; 3],

    /// P - A vector
    pub pa: [f64; 3],

    /// A - B vector (for HTR)
    pub ab: [f64; 3],

    /// Overlap prefactor: K_ij = exp(-mu_ij * |A-B|^2)
    pub k_ij: f64,

    // === Ket pair (k, l) centered at C, D ===
    /// Combined exponent: q = alpha_k + alpha_l
    pub q: f64,

    /// Product center: Q = (alpha_k * C + alpha_l * D) / q
    pub center_q: [f64; 3],

    /// Q - C vector
    pub qc: [f64; 3],

    /// C - D vector (for HTR)
    pub cd: [f64; 3],

    /// Overlap prefactor: K_kl = exp(-mu_kl * |C-D|^2)
    pub k_kl: f64,

    // === Combined quantities ===
    /// P - Q vector
    pub pq: [f64; 3],

    /// |P - Q|^2
    pub pq_squared: f64,

    /// Reduced exponent: rho = p * q / (p + q)
    pub rho: f64,

    /// Rys quadrature argument: T = rho * |P - Q|^2
    pub t: f64,

    /// ERI prefactor: 2 * pi^(5/2) / (p * q * sqrt(p + q)) * K_ij * K_kl
    pub prefactor: f64,

    /// 1 / (2 * p) - coefficient for VRR
    pub one_over_2p: f64,

    /// 1 / (2 * q) - coefficient for VRR
    pub one_over_2q: f64,
}

impl GaussianProduct2e {
    /// Create a GaussianProduct2e from four primitive Gaussians
    ///
    /// # Arguments
    ///
    /// * `alpha_i` - Exponent of first Gaussian (bra, center A)
    /// * `center_a` - Center of first Gaussian [x, y, z]
    /// * `alpha_j` - Exponent of second Gaussian (bra, center B)
    /// * `center_b` - Center of second Gaussian [x, y, z]
    /// * `alpha_k` - Exponent of third Gaussian (ket, center C)
    /// * `center_c` - Center of third Gaussian [x, y, z]
    /// * `alpha_l` - Exponent of fourth Gaussian (ket, center D)
    /// * `center_d` - Center of fourth Gaussian [x, y, z]
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_core::integrals::GaussianProduct2e;
    ///
    /// let a = [0.0, 0.0, 0.0];
    /// let b = [0.0, 0.0, 0.0];
    /// let c = [0.0, 0.0, 1.4];
    /// let d = [0.0, 0.0, 1.4];
    ///
    /// let gp = GaussianProduct2e::new(1.0, &a, 1.0, &b, 1.0, &c, 1.0, &d);
    /// assert!((gp.p - 2.0).abs() < 1e-12);
    /// assert!((gp.q - 2.0).abs() < 1e-12);
    /// ```
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        alpha_i: f64,
        center_a: &[f64; 3],
        alpha_j: f64,
        center_b: &[f64; 3],
        alpha_k: f64,
        center_c: &[f64; 3],
        alpha_l: f64,
        center_d: &[f64; 3],
    ) -> Self {
        // === Bra pair (i, j) ===
        let p = alpha_i + alpha_j;
        let one_over_2p = 0.5 / p;
        let mu_ij = alpha_i * alpha_j / p;

        // Product center P
        let center_p = [
            (alpha_i * center_a[0] + alpha_j * center_b[0]) / p,
            (alpha_i * center_a[1] + alpha_j * center_b[1]) / p,
            (alpha_i * center_a[2] + alpha_j * center_b[2]) / p,
        ];

        // PA = P - A
        let pa = [
            center_p[0] - center_a[0],
            center_p[1] - center_a[1],
            center_p[2] - center_a[2],
        ];

        // AB = A - B (for HTR)
        let ab = [
            center_a[0] - center_b[0],
            center_a[1] - center_b[1],
            center_a[2] - center_b[2],
        ];

        // |A - B|^2
        let ab_squared = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];

        // K_ij = exp(-mu_ij * |A-B|^2)
        let k_ij = (-mu_ij * ab_squared).exp();

        // === Ket pair (k, l) ===
        let q = alpha_k + alpha_l;
        let one_over_2q = 0.5 / q;
        let mu_kl = alpha_k * alpha_l / q;

        // Product center Q
        let center_q = [
            (alpha_k * center_c[0] + alpha_l * center_d[0]) / q,
            (alpha_k * center_c[1] + alpha_l * center_d[1]) / q,
            (alpha_k * center_c[2] + alpha_l * center_d[2]) / q,
        ];

        // QC = Q - C
        let qc = [
            center_q[0] - center_c[0],
            center_q[1] - center_c[1],
            center_q[2] - center_c[2],
        ];

        // CD = C - D (for HTR)
        let cd = [
            center_c[0] - center_d[0],
            center_c[1] - center_d[1],
            center_c[2] - center_d[2],
        ];

        // |C - D|^2
        let cd_squared = cd[0] * cd[0] + cd[1] * cd[1] + cd[2] * cd[2];

        // K_kl = exp(-mu_kl * |C-D|^2)
        let k_kl = (-mu_kl * cd_squared).exp();

        // === Combined quantities ===
        // PQ = P - Q
        let pq = [
            center_p[0] - center_q[0],
            center_p[1] - center_q[1],
            center_p[2] - center_q[2],
        ];

        // |P - Q|^2
        let pq_squared = pq[0] * pq[0] + pq[1] * pq[1] + pq[2] * pq[2];

        // rho = p * q / (p + q)
        let rho = p * q / (p + q);

        // T = rho * |P - Q|^2 (Rys quadrature argument)
        let t = rho * pq_squared;

        // Prefactor: 2 * pi^(5/2) / (p * q * sqrt(p + q)) * K_ij * K_kl
        // Reference: libcint g2e.c line 4441
        // pi^(5/2) = pi^2 * sqrt(pi), avoids powf()
        let prefactor = 2.0 * PI * PI * PI.sqrt() / (p * q * (p + q).sqrt()) * k_ij * k_kl;

        Self {
            p,
            center_p,
            pa,
            ab,
            k_ij,
            q,
            center_q,
            qc,
            cd,
            k_kl,
            pq,
            pq_squared,
            rho,
            t,
            prefactor,
            one_over_2p,
            one_over_2q,
        }
    }

    /// Fast constructor that accepts pre-computed shell pair data.
    ///
    /// Avoids recomputing |A-B|², |C-D|², AB, CD vectors, and exp() calls
    /// for K_ij and K_kl, which are constant across all primitive quartets
    /// within a shell quartet (or depend only on the bra/ket pair).
    ///
    /// # Arguments
    ///
    /// * `alpha_i`, `alpha_j` - Bra primitive exponents
    /// * `center_a`, `center_b` - Bra primitive centers
    /// * `alpha_k`, `alpha_l` - Ket primitive exponents
    /// * `center_c`, `center_d` - Ket primitive centers
    /// * `ab` - Pre-computed A - B vector
    /// * `cd` - Pre-computed C - D vector
    /// * `k_ij` - Pre-computed exp(-mu_ij * |A-B|²)
    /// * `k_kl` - Pre-computed exp(-mu_kl * |C-D|²)
    #[allow(clippy::too_many_arguments)]
    #[inline]
    pub fn new_prescreened(
        alpha_i: f64,
        center_a: &[f64; 3],
        alpha_j: f64,
        center_b: &[f64; 3],
        alpha_k: f64,
        center_c: &[f64; 3],
        alpha_l: f64,
        center_d: &[f64; 3],
        ab: [f64; 3],
        cd: [f64; 3],
        k_ij: f64,
        k_kl: f64,
    ) -> Self {
        // === Bra pair (i, j) ===
        let p = alpha_i + alpha_j;
        let one_over_2p = 0.5 / p;

        // Product center P
        let center_p = [
            (alpha_i * center_a[0] + alpha_j * center_b[0]) / p,
            (alpha_i * center_a[1] + alpha_j * center_b[1]) / p,
            (alpha_i * center_a[2] + alpha_j * center_b[2]) / p,
        ];

        // PA = P - A
        let pa = [
            center_p[0] - center_a[0],
            center_p[1] - center_a[1],
            center_p[2] - center_a[2],
        ];

        // === Ket pair (k, l) ===
        let q = alpha_k + alpha_l;
        let one_over_2q = 0.5 / q;

        // Product center Q
        let center_q = [
            (alpha_k * center_c[0] + alpha_l * center_d[0]) / q,
            (alpha_k * center_c[1] + alpha_l * center_d[1]) / q,
            (alpha_k * center_c[2] + alpha_l * center_d[2]) / q,
        ];

        // QC = Q - C
        let qc = [
            center_q[0] - center_c[0],
            center_q[1] - center_c[1],
            center_q[2] - center_c[2],
        ];

        // === Combined quantities ===
        let pq = [
            center_p[0] - center_q[0],
            center_p[1] - center_q[1],
            center_p[2] - center_q[2],
        ];
        let pq_squared = pq[0] * pq[0] + pq[1] * pq[1] + pq[2] * pq[2];
        let rho = p * q / (p + q);
        let t = rho * pq_squared;

        // pi^(5/2) = pi^2 * sqrt(pi)
        let prefactor = 2.0 * PI * PI * PI.sqrt() / (p * q * (p + q).sqrt()) * k_ij * k_kl;

        Self {
            p,
            center_p,
            pa,
            ab,
            k_ij,
            q,
            center_q,
            qc,
            cd,
            k_kl,
            pq,
            pq_squared,
            rho,
            t,
            prefactor,
            one_over_2p,
            one_over_2q,
        }
    }

    /// Check if all four primitives are at the same center (within tolerance)
    #[inline]
    pub fn is_same_center(&self) -> bool {
        self.pq_squared < 1e-20
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    #[test]
    fn test_same_center() {
        let center = [0.0, 0.0, 0.0];
        let gp = GaussianProduct2e::new(1.0, &center, 1.0, &center, 1.0, &center, 1.0, &center);

        assert!((gp.p - 2.0).abs() < TOL);
        assert!((gp.q - 2.0).abs() < TOL);
        assert!(gp.t.abs() < TOL); // T = 0 for same center
        assert!(gp.is_same_center());

        // Check prefactor is finite and positive
        assert!(gp.prefactor > 0.0);
        assert!(gp.prefactor.is_finite());
    }

    #[test]
    fn test_two_centers() {
        // Two H atoms along z-axis
        let a = [0.0, 0.0, 0.0];
        let c = [0.0, 0.0, 2.0];

        let gp = GaussianProduct2e::new(1.0, &a, 1.0, &a, 1.0, &c, 1.0, &c);

        assert!((gp.p - 2.0).abs() < TOL);
        assert!((gp.q - 2.0).abs() < TOL);

        // P = A = (0,0,0), Q = C = (0,0,2)
        assert!((gp.center_p[2]).abs() < TOL);
        assert!((gp.center_q[2] - 2.0).abs() < TOL);

        // |P - Q|^2 = 4
        assert!((gp.pq_squared - 4.0).abs() < TOL);

        // rho = 2*2/(2+2) = 1
        assert!((gp.rho - 1.0).abs() < TOL);

        // T = rho * |P-Q|^2 = 4
        assert!((gp.t - 4.0).abs() < TOL);

        assert!(!gp.is_same_center());
    }

    #[test]
    fn test_four_centers() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];

        let gp = GaussianProduct2e::new(1.0, &a, 2.0, &b, 1.5, &c, 0.5, &d);

        // p = 1 + 2 = 3
        assert!((gp.p - 3.0).abs() < TOL);

        // q = 1.5 + 0.5 = 2
        assert!((gp.q - 2.0).abs() < TOL);

        // P = (1*0 + 2*1, 0, 0) / 3 = (2/3, 0, 0)
        assert!((gp.center_p[0] - 2.0 / 3.0).abs() < TOL);
        assert!((gp.center_p[1]).abs() < TOL);
        assert!((gp.center_p[2]).abs() < TOL);

        // Q = (0, 1.5*1 + 0.5*0, 0.5*1) / 2 = (0, 0.75, 0.25)
        assert!((gp.center_q[0]).abs() < TOL);
        assert!((gp.center_q[1] - 0.75).abs() < TOL);
        assert!((gp.center_q[2] - 0.25).abs() < TOL);

        // Prefactor should be positive and finite
        assert!(gp.prefactor > 0.0);
        assert!(gp.prefactor.is_finite());
    }

    #[test]
    fn test_k_factors() {
        // K factors should be 1.0 when centers are same, < 1.0 otherwise
        let a = [0.0, 0.0, 0.0];
        let gp_same = GaussianProduct2e::new(1.0, &a, 1.0, &a, 1.0, &a, 1.0, &a);
        assert!((gp_same.k_ij - 1.0).abs() < TOL);
        assert!((gp_same.k_kl - 1.0).abs() < TOL);

        let b = [1.0, 0.0, 0.0];
        let gp_diff = GaussianProduct2e::new(1.0, &a, 1.0, &b, 1.0, &a, 1.0, &b);
        assert!(gp_diff.k_ij < 1.0);
        assert!(gp_diff.k_kl < 1.0);
    }

    #[test]
    fn test_pa_qc_vectors() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [0.0, 0.0, 0.0];
        let d = [0.0, 2.0, 0.0];

        // With equal exponents, P = midpoint of A and B
        let gp = GaussianProduct2e::new(1.0, &a, 1.0, &b, 1.0, &c, 1.0, &d);

        // P = (1, 0, 0)
        assert!((gp.center_p[0] - 1.0).abs() < TOL);

        // PA = P - A = (1, 0, 0)
        assert!((gp.pa[0] - 1.0).abs() < TOL);

        // AB = A - B = (-2, 0, 0)
        assert!((gp.ab[0] - (-2.0)).abs() < TOL);

        // Q = (0, 1, 0)
        assert!((gp.center_q[1] - 1.0).abs() < TOL);

        // QC = Q - C = (0, 1, 0)
        assert!((gp.qc[1] - 1.0).abs() < TOL);

        // CD = C - D = (0, -2, 0)
        assert!((gp.cd[1] - (-2.0)).abs() < TOL);
    }
}
