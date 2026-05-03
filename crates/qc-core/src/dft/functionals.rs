//! Exchange-correlation functionals for Kohn-Sham DFT.
//!
//! This module implements LDA-level exchange-correlation functionals:
//!
//! - **Slater exchange**: Dirac (1930), Slater (1951) local exchange
//! - **VWN5 correlation**: Vosko, Wilk & Nusair parametrization V (1980)
//! - **LDA composite**: Slater + VWN5
//!
//! # Convention
//!
//! - `epsilon_xc(rho)`: energy density per electron
//! - `v_xc(rho) = d(rho * epsilon_xc) / d(rho)`: the xc potential
//! - Total xc energy: `E_xc = sum_g w_g * rho(r_g) * epsilon_xc(rho(r_g))`
//!
//! # References
//!
//! - Dirac, P. A. M. (1930). Proc. Camb. Phil. Soc. 26, 376.
//! - Slater, J. C. (1951). Phys. Rev. 81, 385.
//! - Vosko, S. H., Wilk, L. & Nusair, M. (1980). Can. J. Phys. 58, 1200.

use std::f64::consts::PI;

// =============================================================================
// Constants
// =============================================================================

/// Density threshold below which xc contributions are set to zero.
///
/// Prevents numerical issues from `rho^{1/3}` at near-zero densities.
/// Matches PySCF/libxc behavior which returns zero for `rho < 1e-20`.
pub const DENSITY_THRESHOLD: f64 = 1e-20;

// =============================================================================
// ExchangeCorrelation Trait
// =============================================================================

/// Trait for exchange-correlation functionals.
///
/// Implementations provide the energy density per electron (`epsilon_xc`)
/// and the exchange-correlation potential (`v_xc = d(rho*epsilon_xc)/drho`)
/// evaluated at an array of grid-point densities.
///
/// # Convention
///
/// - `exc[i]` = `epsilon_xc(rho[i])`: energy density per electron
/// - `vxc[i]` = `v_xc(rho[i])` = `d(rho * epsilon_xc) / d(rho)`: the potential
///
/// The total xc energy is: `E_xc = sum_g w_g * rho(r_g) * exc(rho(r_g))`
///
/// # GGA Extension
///
/// GGA functionals additionally depend on `sigma = |grad rho|^2`.
/// Override `eval_xc_gga()` to provide gradient-dependent evaluation.
/// The default implementation falls back to `eval_xc()` (LDA behavior).
pub trait ExchangeCorrelation {
    /// Evaluate xc energy density and potential at given density values (LDA).
    ///
    /// For each grid point with density `rho[i]`:
    /// - `exc[i]` = exchange-correlation energy density per electron
    /// - `vxc[i]` = exchange-correlation potential
    ///
    /// Points with `rho <= DENSITY_THRESHOLD` are set to zero.
    fn eval_xc(&self, rho: &[f64], exc: &mut [f64], vxc: &mut [f64]);

    /// Evaluate xc energy density and potential for GGA functionals.
    ///
    /// For each grid point with density `rho[i]` and gradient squared `sigma[i]`:
    /// - `exc[i]` = exchange-correlation energy density per electron
    /// - `vrho[i]` = d(rho * epsilon_xc) / d(rho)
    /// - `vsigma[i]` = d(rho * epsilon_xc) / d(sigma)
    ///
    /// The GGA V_xc matrix includes both rho and sigma contributions:
    /// ```text
    /// V_xc[mu,nu] += w_g * vrho(g) * chi_mu(g) * chi_nu(g)
    ///             + 2 * w_g * vsigma(g) * (grad_chi_mu . grad_rho) * chi_nu
    ///             + 2 * w_g * vsigma(g) * chi_mu * (grad_chi_nu . grad_rho)
    /// ```
    ///
    /// Default implementation calls `eval_xc()` and sets `vsigma = 0`.
    fn eval_xc_gga(
        &self,
        rho: &[f64],
        _sigma: &[f64],
        exc: &mut [f64],
        vrho: &mut [f64],
        vsigma: &mut [f64],
    ) {
        self.eval_xc(rho, exc, vrho);
        for v in vsigma.iter_mut() {
            *v = 0.0;
        }
    }

    /// Evaluate second functional derivatives (XC kernel fxc).
    ///
    /// For LDA functionals, only `v2rho2` is meaningful; `v2rhosigma` and
    /// `v2sigma2` are set to zero by the default implementation.
    ///
    /// For GGA functionals, all three components are needed:
    /// - `v2rho2[i]`     = d²(rho * eps_xc) / d(rho)²
    /// - `v2rhosigma[i]` = d²(rho * eps_xc) / d(rho) d(sigma)
    /// - `v2sigma2[i]`   = d²(rho * eps_xc) / d(sigma)²
    ///
    /// The 0.5 GGA scale factor for vsigma terms is NOT applied here —
    /// the caller (Hessian assembler) handles it.
    ///
    /// # References
    ///
    /// - Phase 5 plan, Section 4f: DFT XC Second Derivatives
    /// - PySCF `hessian/rks.py`: uses fxc in `_get_vxc_deriv2`
    fn eval_xc_second_deriv(
        &self,
        rho: &[f64],
        sigma: &[f64],
        v2rho2: &mut [f64],
        v2rhosigma: &mut [f64],
        v2sigma2: &mut [f64],
    );

    /// Whether this functional requires the density gradient (GGA).
    ///
    /// Returns `false` for LDA, `true` for GGA/hybrid.
    fn needs_gradient(&self) -> bool;

    /// Fraction of Hartree-Fock exchange to include (0.0 for pure DFT, 0.20 for B3LYP).
    fn hf_exchange_fraction(&self) -> f64;

    /// Human-readable name of the functional.
    fn name(&self) -> &str;
}

// =============================================================================
// Slater Exchange
// =============================================================================

/// Slater (Dirac) local exchange functional.
///
/// The Slater exchange energy for a uniform electron gas is:
///
/// ```text
/// E_x[rho] = -C_x * integral rho(r)^{4/3} dr
/// ```
///
/// where `C_x = (3/4)(3/pi)^{1/3}`.
///
/// The energy density per electron and potential are:
///
/// ```text
/// epsilon_x(rho) = -C_x * rho^{1/3}
/// v_x(rho)       = -(4/3) * C_x * rho^{1/3}
/// ```
///
/// Note that `v_x = (4/3) * epsilon_x`. The 4/3 factor arises from the
/// functional derivative chain rule and is essential for correct
/// self-consistent KS-DFT.
///
/// # References
///
/// - Dirac, P. A. M. (1930). Proc. Camb. Phil. Soc. 26, 376.
/// - Slater, J. C. (1951). Phys. Rev. 81, 385.
pub struct SlaterExchange;

impl SlaterExchange {
    /// Exchange constant `C_x = (3/4)(3/pi)^{1/3}`.
    pub const C_X: f64 = 0.7385587663820224;

    /// `(4/3) * C_x`, used in the exchange potential.
    const FOUR_THIRDS_C_X: f64 = 0.9847450218426965;

    /// Evaluate Slater exchange energy density and potential at a single density.
    ///
    /// Returns `(epsilon_x, v_x)`.
    ///
    /// ```text
    /// epsilon_x(rho) = -C_x * rho^{1/3}
    /// v_x(rho) = -(4/3) * C_x * rho^{1/3}
    /// ```
    #[inline]
    pub fn eval(&self, rho: f64) -> (f64, f64) {
        if rho <= DENSITY_THRESHOLD {
            return (0.0, 0.0);
        }
        let rho_third = rho.cbrt();
        let eps_x = -Self::C_X * rho_third;
        let v_x = -Self::FOUR_THIRDS_C_X * rho_third;
        (eps_x, v_x)
    }
}

impl SlaterExchange {
    /// Evaluate Slater exchange second derivative at a single density.
    ///
    /// Returns `d²(rho * epsilon_x) / d(rho)²`.
    ///
    /// ```text
    /// E_x = -C_x * rho^{4/3}
    /// dE_x/drho = -(4/3) * C_x * rho^{1/3}
    /// d²E_x/drho² = -(4/9) * C_x * rho^{-2/3}
    /// ```
    ///
    /// Low-density guard: if rho < DENSITY_THRESHOLD, returns 0.
    #[inline]
    pub fn eval_fxc(&self, rho: f64) -> f64 {
        if rho <= DENSITY_THRESHOLD {
            return 0.0;
        }
        let rho_13 = rho.cbrt();
        // d²(rho * eps_x)/drho² = -(4/9) * C_x * rho^{-2/3}
        -4.0 / 9.0 * Self::C_X / (rho_13 * rho_13)
    }
}

impl ExchangeCorrelation for SlaterExchange {
    fn eval_xc(&self, rho: &[f64], exc: &mut [f64], vxc: &mut [f64]) {
        for i in 0..rho.len() {
            let (e, v) = self.eval(rho[i]);
            exc[i] = e;
            vxc[i] = v;
        }
    }

    fn eval_xc_second_deriv(
        &self,
        rho: &[f64],
        _sigma: &[f64],
        v2rho2: &mut [f64],
        v2rhosigma: &mut [f64],
        v2sigma2: &mut [f64],
    ) {
        for i in 0..rho.len() {
            v2rho2[i] = self.eval_fxc(rho[i]);
            v2rhosigma[i] = 0.0;
            v2sigma2[i] = 0.0;
        }
    }

    fn needs_gradient(&self) -> bool {
        false
    }

    fn hf_exchange_fraction(&self) -> f64 {
        0.0
    }

    fn name(&self) -> &str {
        "Slater Exchange"
    }
}

// =============================================================================
// VWN5 Correlation
// =============================================================================

/// Vosko-Wilk-Nusair parametrization V (VWN5) correlation functional.
///
/// Parametrizes the correlation energy of the uniform electron gas as a
/// function of the Wigner-Seitz radius `r_s = (3/(4*pi*rho))^{1/3}`.
///
/// The functional form is (VWN 1980, Eq. 4.4):
///
/// ```text
/// epsilon_c(r_s) = A * [ ln(x^2/X(x))
///                      + 2*b/Q * arctan(Q/(2*x+b))
///                      - b*x0/X(x0) * ( ln((x-x0)^2/X(x))
///                                      + 2*(b+2*x0)/Q * arctan(Q/(2*x+b)) ) ]
/// ```
///
/// where `x = sqrt(r_s)`, `X(x) = x^2 + b*x + c`, and `Q = sqrt(4c - b^2)`.
///
/// The correlation potential is:
///
/// ```text
/// v_c = epsilon_c - (r_s/3) * d(epsilon_c)/d(r_s)
/// ```
///
/// # A Parameter Convention
///
/// The original VWN paper defines `A = 0.0621814` (Hartree) for the
/// spin-pair correlation energy. PySCF and libxc use `A/2 = 0.0310907`
/// (per-electron convention). We follow PySCF/libxc.
///
/// # Reference
///
/// Vosko, S. H., Wilk, L. & Nusair, M. (1980). Can. J. Phys. 58, 1200,
/// Eq. 4.4, Table 4.4 parametrization V.
pub struct Vwn5Correlation;

impl Vwn5Correlation {
    // VWN5 parameters (spin-unpolarized, parametrization V)
    // Table 4.4 of VWN (1980).

    /// `A = 0.0621814 / 2` (per-electron convention matching PySCF/libxc).
    pub const A: f64 = 0.0310907;

    /// Fitting parameter `x0`.
    pub const X0: f64 = -0.10498;

    /// Fitting parameter `b`.
    pub const B: f64 = 3.72744;

    /// Fitting parameter `c`.
    pub const C: f64 = 12.9352;

    // Derived constants (precomputed for efficiency).

    /// `Q = sqrt(4c - b^2)`.
    pub const Q: f64 = 6.15199081975908;

    /// `X(x0) = x0^2 + b*x0 + c`.
    pub const X_X0: f64 = 12.5549141492;

    /// Evaluate VWN5 correlation energy density and potential at a single density.
    ///
    /// Returns `(epsilon_c, v_c)`.
    ///
    /// # Algorithm
    ///
    /// 1. Compute `r_s` (Wigner-Seitz radius) and `x = sqrt(r_s)`
    /// 2. Evaluate `epsilon_c` from VWN Eq. 4.4
    /// 3. Compute `d(epsilon_c)/dx` analytically
    /// 4. Chain rule: `d(epsilon_c)/d(r_s) = d(epsilon_c)/dx / (2*x)`
    /// 5. Potential: `v_c = epsilon_c - (r_s/3) * d(epsilon_c)/d(r_s)`
    #[inline]
    pub fn eval(&self, rho: f64) -> (f64, f64) {
        if rho <= DENSITY_THRESHOLD {
            return (0.0, 0.0);
        }

        // Wigner-Seitz radius: r_s = (3/(4*pi*rho))^{1/3}
        let r_s = (3.0 / (4.0 * PI * rho)).cbrt();
        let x = r_s.sqrt();

        // X(x) = x^2 + b*x + c
        let x_sq = x * x; // = r_s
        let x_val = x_sq + Self::B * x + Self::C;

        // Common arctan argument: arctan(Q / (2*x + b))
        let two_x_plus_b = 2.0 * x + Self::B;
        let atan_arg = (Self::Q / two_x_plus_b).atan();

        // =====================================================================
        // epsilon_c: VWN Eq. 4.4
        // Three terms: ln term, arctan term, correction term involving x0
        // =====================================================================

        // Term 1: ln(x^2 / X(x))
        let term1 = (x_sq / x_val).ln();

        // Term 2: 2*b/Q * arctan(Q/(2*x+b))
        let term2 = 2.0 * Self::B / Self::Q * atan_arg;

        // Term 3: -b*x0/X(x0) * [ ln((x-x0)^2/X(x)) + 2*(b+2*x0)/Q * arctan(Q/(2*x+b)) ]
        let x_minus_x0 = x - Self::X0;
        let term3_ln = (x_minus_x0 * x_minus_x0 / x_val).ln();
        let term3_atan = 2.0 * (Self::B + 2.0 * Self::X0) / Self::Q * atan_arg;
        let term3 = -Self::B * Self::X0 / Self::X_X0 * (term3_ln + term3_atan);

        let eps_c = Self::A * (term1 + term2 + term3);

        // =====================================================================
        // d(epsilon_c)/dx: analytical derivative
        //
        // d/dx[ln(x^2/X(x))] = 2/x - (2x+b)/X(x)
        // d/dx[arctan(Q/(2x+b))] = -2Q / ((2x+b)^2 + Q^2)
        // =====================================================================

        // Derivative of arctan term w.r.t. x
        let denom = two_x_plus_b * two_x_plus_b + Self::Q * Self::Q;
        let d_atan = -2.0 * Self::Q / denom;

        // d(term1)/dx = 2/x - (2x+b)/X(x)
        let d_term1 = 2.0 / x - two_x_plus_b / x_val;

        // d(term2)/dx = 2*b/Q * d_atan
        let d_term2 = 2.0 * Self::B / Self::Q * d_atan;

        // d(term3)/dx involves:
        //   d/dx[ln((x-x0)^2/X(x))] = 2/(x-x0) - (2x+b)/X(x)
        //   d/dx[arctan term] = 2*(b+2*x0)/Q * d_atan
        let d_term3 = -Self::B * Self::X0 / Self::X_X0
            * (2.0 / x_minus_x0 - two_x_plus_b / x_val
                + 2.0 * (Self::B + 2.0 * Self::X0) / Self::Q * d_atan);

        let deps_dx = Self::A * (d_term1 + d_term2 + d_term3);

        // Chain rule: d(eps_c)/d(r_s) = d(eps_c)/dx * dx/d(r_s)
        // dx/d(r_s) = 1/(2*x)
        let deps_drs = deps_dx / (2.0 * x);

        // Correlation potential: v_c = eps_c - (r_s/3) * d(eps_c)/d(r_s)
        let v_c = eps_c - (r_s / 3.0) * deps_drs;

        (eps_c, v_c)
    }
}

impl Vwn5Correlation {
    /// Evaluate VWN5 second derivative d²(rho * eps_c)/drho² at a single density.
    ///
    /// # Algorithm
    ///
    /// Chain rule through r_s and x = sqrt(r_s):
    ///
    /// ```text
    /// v_c = eps_c - (r_s/3) * d(eps_c)/d(r_s)
    /// d(v_c)/d(rho) = d(r_s)/d(rho) * [(2/3) * d(eps_c)/d(r_s) - (r_s/3) * d²(eps_c)/d(r_s)²]
    /// ```
    ///
    /// where d(r_s)/d(rho) = -r_s / (3*rho).
    ///
    /// # Reference
    ///
    /// VWN (1980) Eq. 4.4, differentiated twice.
    #[inline]
    pub fn eval_fxc(&self, rho: f64) -> f64 {
        if rho <= DENSITY_THRESHOLD {
            return 0.0;
        }

        let r_s = (3.0 / (4.0 * PI * rho)).cbrt();
        let x = r_s.sqrt();

        let x_sq = x * x; // = r_s
        let x_val = x_sq + Self::B * x + Self::C;

        let two_x_plus_b = 2.0 * x + Self::B;
        let x_minus_x0 = x - Self::X0;
        let b_plus_2x0 = Self::B + 2.0 * Self::X0;

        // Pre-compute common sub-expressions
        let q_sq = Self::Q * Self::Q;
        let denom = two_x_plus_b * two_x_plus_b + q_sq; // (2x+b)^2 + Q^2

        // =====================================================================
        // First derivative deps_c/dx (same as in eval(), kept for chain rule)
        // =====================================================================

        // d(arctan(Q/(2x+b)))/dx = -2Q / denom
        let d_atan = -2.0 * Self::Q / denom;

        // d(term1)/dx = 2/x - (2x+b)/X(x)
        let d_term1 = 2.0 / x - two_x_plus_b / x_val;

        // d(term2)/dx = 2*b/Q * d_atan
        let d_term2 = 2.0 * Self::B / Self::Q * d_atan;

        // d(term3)/dx:
        // d/dx[ln((x-x0)^2/X(x))] = 2/(x-x0) - (2x+b)/X(x)
        // d/dx[arctan_term] = 2*(b+2*x0)/Q * d_atan
        let bx0_over_xx0 = Self::B * Self::X0 / Self::X_X0;
        let d_term3 = -bx0_over_xx0
            * (2.0 / x_minus_x0 - two_x_plus_b / x_val + 2.0 * b_plus_2x0 / Self::Q * d_atan);

        let deps_dx = Self::A * (d_term1 + d_term2 + d_term3);

        // =====================================================================
        // Second derivative d²eps_c/dx² — differentiate d_term{1,2,3} w.r.t. x
        // =====================================================================

        // d²(arctan(Q/(2x+b)))/dx² = 8*Q*(2x+b) / denom²
        let d2_atan = 8.0 * Self::Q * two_x_plus_b / (denom * denom);

        // dX/dx = 2x + b = two_x_plus_b
        // d²X/dx² = 2

        // d²(term1)/dx²:
        // d/dx[2/x] = -2/x²
        // d/dx[-(2x+b)/X(x)] = [-2*X(x) - (2x+b)*dX/dx] / X(x)² + correction
        //   = [-2*X(x) + (2x+b)*(2x+b)] / X(x)²  ... but let me be precise:
        // d/dx[(2x+b)/X(x)] = [2*X(x) - (2x+b)²] / X(x)²
        let d2_term1 =
            -2.0 / (x * x) - (2.0 * x_val - two_x_plus_b * two_x_plus_b) / (x_val * x_val);

        // d²(term2)/dx² = 2*b/Q * d2_atan
        let d2_term2 = 2.0 * Self::B / Self::Q * d2_atan;

        // d²(term3)/dx²:
        // d/dx[2/(x-x0)] = -2/(x-x0)²
        // d/dx[-(2x+b)/X(x)] already computed above
        // d/dx[2*(b+2*x0)/Q * d_atan] = 2*(b+2*x0)/Q * d2_atan
        let d2_term3 = -bx0_over_xx0
            * (-2.0 / (x_minus_x0 * x_minus_x0)
                - (2.0 * x_val - two_x_plus_b * two_x_plus_b) / (x_val * x_val)
                + 2.0 * b_plus_2x0 / Self::Q * d2_atan);

        let d2eps_dx2 = Self::A * (d2_term1 + d2_term2 + d2_term3);

        // =====================================================================
        // Chain rule to get d²(rho*eps_c)/drho²
        // =====================================================================
        //
        // deps_c/drs = deps_dx / (2x)
        // d²eps_c/drs² = [d²eps_dx2 - deps_dx/x] / (4*x²)
        //
        // dv_c/drho = drs/drho * [(2/3)*deps_c/drs - (r_s/3)*d²eps_c/drs²]
        // where drs/drho = -r_s / (3*rho)

        let deps_drs = deps_dx / (2.0 * x);
        let d2eps_drs2 = (d2eps_dx2 - deps_dx / x) / (4.0 * x_sq);

        let drs_drho = -r_s / (3.0 * rho);

        drs_drho * (2.0 / 3.0 * deps_drs - r_s / 3.0 * d2eps_drs2)
    }
}

impl ExchangeCorrelation for Vwn5Correlation {
    fn eval_xc(&self, rho: &[f64], exc: &mut [f64], vxc: &mut [f64]) {
        for i in 0..rho.len() {
            let (e, v) = self.eval(rho[i]);
            exc[i] = e;
            vxc[i] = v;
        }
    }

    fn eval_xc_second_deriv(
        &self,
        rho: &[f64],
        _sigma: &[f64],
        v2rho2: &mut [f64],
        v2rhosigma: &mut [f64],
        v2sigma2: &mut [f64],
    ) {
        for i in 0..rho.len() {
            v2rho2[i] = self.eval_fxc(rho[i]);
            v2rhosigma[i] = 0.0;
            v2sigma2[i] = 0.0;
        }
    }

    fn needs_gradient(&self) -> bool {
        false
    }

    fn hf_exchange_fraction(&self) -> f64 {
        0.0
    }

    fn name(&self) -> &str {
        "VWN5 Correlation"
    }
}

// =============================================================================
// LDA Composite Functional
// =============================================================================

/// Local Density Approximation: Slater exchange + VWN5 correlation.
///
/// Combines Slater local exchange with VWN parametrization V correlation:
///
/// ```text
/// E_xc[rho] = E_x[rho] + E_c[rho]
///           = -C_x * integral rho^{4/3} dr + integral rho * eps_c(r_s) dr
/// ```
///
/// # References
///
/// - Slater, J. C. (1951). Phys. Rev. 81, 385.
/// - Vosko, S. H., Wilk, L. & Nusair, M. (1980). Can. J. Phys. 58, 1200.
pub struct Lda {
    exchange: SlaterExchange,
    correlation: Vwn5Correlation,
}

impl Lda {
    /// Create a new LDA functional (Slater exchange + VWN5 correlation).
    pub fn new() -> Self {
        Self {
            exchange: SlaterExchange,
            correlation: Vwn5Correlation,
        }
    }
}

impl Default for Lda {
    fn default() -> Self {
        Self::new()
    }
}

impl ExchangeCorrelation for Lda {
    fn eval_xc(&self, rho: &[f64], exc: &mut [f64], vxc: &mut [f64]) {
        for i in 0..rho.len() {
            let (ex, vx) = self.exchange.eval(rho[i]);
            let (ec, vc) = self.correlation.eval(rho[i]);
            exc[i] = ex + ec;
            vxc[i] = vx + vc;
        }
    }

    fn eval_xc_second_deriv(
        &self,
        rho: &[f64],
        sigma: &[f64],
        v2rho2: &mut [f64],
        v2rhosigma: &mut [f64],
        v2sigma2: &mut [f64],
    ) {
        // LDA: sum of Slater and VWN5 second derivatives
        self.exchange
            .eval_xc_second_deriv(rho, sigma, v2rho2, v2rhosigma, v2sigma2);
        let mut v2rho2_c = vec![0.0; rho.len()];
        let mut v2rhosigma_c = vec![0.0; rho.len()];
        let mut v2sigma2_c = vec![0.0; rho.len()];
        self.correlation.eval_xc_second_deriv(
            rho,
            sigma,
            &mut v2rho2_c,
            &mut v2rhosigma_c,
            &mut v2sigma2_c,
        );
        for i in 0..rho.len() {
            v2rho2[i] += v2rho2_c[i];
            v2rhosigma[i] += v2rhosigma_c[i];
            v2sigma2[i] += v2sigma2_c[i];
        }
    }

    fn needs_gradient(&self) -> bool {
        false
    }

    fn hf_exchange_fraction(&self) -> f64 {
        0.0
    }

    fn name(&self) -> &str {
        "LDA (Slater + VWN5)"
    }
}

// =============================================================================
// Becke88 Exchange (GGA)
// =============================================================================

/// Becke88 gradient-corrected exchange functional.
///
/// The Becke88 exchange functional adds a gradient correction to the
/// Slater (LDA) exchange. The full B88 exchange energy density per electron:
///
/// ```text
/// eps_x^B88(rho, sigma) = eps_x^Slater(rho) + dEps_x^B88(rho, sigma)
/// ```
///
/// where the gradient correction is:
///
/// ```text
/// dEps_x^B88 = -beta * rho^{1/3} * x^2 / (1 + 6*beta*x*arcsinh(x))
/// x = |grad rho| / rho^{4/3}   (reduced density gradient)
/// sigma = |grad rho|^2
/// beta = 0.0042
/// ```
///
/// Note: Like libxc's `gga_x_b88`, this INCLUDES the Slater exchange.
/// At sigma=0, it exactly reduces to Slater exchange.
///
/// # References
///
/// - Becke, A. D. (1988). Phys. Rev. A, 38, 3098.
/// - PySCF: `dft.libxc.eval_xc('gga_x_b88', ...)`
pub struct Becke88Exchange;

impl Becke88Exchange {
    /// Becke88 empirical parameter beta = 0.0042.
    const BETA: f64 = 0.0042;

    /// Evaluate Becke88 exchange at a single grid point (closed-shell).
    ///
    /// Returns `(eps_x, vrho, vsigma)` where:
    /// - `eps_x` = exchange energy density per electron (includes Slater)
    /// - `vrho` = d(rho * eps_x) / d(rho)
    /// - `vsigma` = d(rho * eps_x) / d(sigma)
    ///
    /// # Spin-density formulation
    ///
    /// Following libxc convention, exchange functionals are evaluated using
    /// spin-density variables for the closed-shell case:
    ///
    /// ```text
    /// rho_a = rho/2, sigma_a = sigma/4
    /// E_x = 2 * integral(rho_a * eps_x^spin(rho_a, sigma_a)) dr
    /// eps_x(rho, sigma) = eps_x^spin(rho/2, sigma/4)
    /// ```
    ///
    /// The spin-polarized Slater exchange uses the 2^{1/3} scaling factor:
    /// ```text
    /// eps_x^spin(rho_a) = -2^{1/3} * C_x * rho_a^{1/3}
    /// ```
    #[inline]
    pub fn eval_gga(&self, rho: f64, sigma: f64) -> (f64, f64, f64) {
        if rho <= DENSITY_THRESHOLD {
            return (0.0, 0.0, 0.0);
        }

        // Spin-density variables for closed shell
        let rho_a = rho / 2.0;
        let sigma_a = sigma.max(0.0) / 4.0;

        let rho_a_13 = rho_a.cbrt();
        let rho_a_43 = rho_a * rho_a_13;

        // Spin-polarized Slater: eps_x^a = -2^{1/3} * C_x * rho_a^{1/3}
        let two_13: f64 = 2.0_f64.cbrt(); // 2^{1/3}
        let eps_slater_a = -two_13 * SlaterExchange::C_X * rho_a_13;

        let grad_rho_a = sigma_a.sqrt();

        if grad_rho_a < 1e-30 {
            // No gradient: pure Slater
            // vrho = d(rho * eps_x)/d(rho) where eps_x = eps_slater_a
            // d(rho * eps)/drho = eps + rho * deps/drho
            // eps = -2^{1/3}*C_x*(rho/2)^{1/3} = -C_x * rho^{1/3}  (by algebra)
            // This is exactly Slater's result.
            let slater = SlaterExchange;
            let (eps_s, v_s) = slater.eval(rho);
            return (eps_s, v_s, 0.0);
        }

        // Reduced density gradient: x_a = |grad rho_a| / rho_a^{4/3}
        let x = grad_rho_a / rho_a_43;
        let x2 = x * x;

        let asinh_x = x.asinh();
        let denom = 1.0 + 6.0 * Self::BETA * x * asinh_x;

        // B88 gradient correction per spin-electron:
        // dEps_a = -beta * rho_a^{1/3} * x^2 / denom
        let d_eps_a = -Self::BETA * rho_a_13 * x2 / denom;

        // Total exchange energy per electron:
        // eps_x = eps_slater_a + d_eps_a
        let eps_x = eps_slater_a + d_eps_a;

        // =====================================================================
        // Derivatives
        //
        // F(rho, sigma) = rho * eps_x
        //               = rho * [eps_slater_a(rho) + d_eps_a(rho, sigma)]
        //
        // where rho_a = rho/2, sigma_a = sigma/4
        //
        // d(rho_a)/d(rho) = 1/2
        // d(sigma_a)/d(sigma) = 1/4
        //
        // F = rho * eps_x^spin(rho/2, sigma/4)
        //
        // vrho = dF/d(rho) = eps_x + rho * d(eps_x)/d(rho)
        //      = eps_x + rho * d(eps_x)/d(rho_a) * (1/2)
        //
        // vsigma = dF/d(sigma) = rho * d(eps_x)/d(sigma_a) * (1/4)
        // =====================================================================

        // d(arcsinh(x))/dx = 1/sqrt(x^2+1)
        let sqrt_x2p1 = (x2 + 1.0).sqrt();
        let d_asinh = 1.0 / sqrt_x2p1;

        let d_denom_dx = 6.0 * Self::BETA * (asinh_x + x * d_asinh);

        // g = x^2/denom
        let dg_dx = (2.0 * x * denom - x2 * d_denom_dx) / (denom * denom);

        // d(eps_x)/d(rho_a) at fixed sigma_a:
        // eps_slater_a = -2^{1/3}*C_x*rho_a^{1/3}
        // d(eps_slater_a)/d(rho_a) = -2^{1/3}*C_x/(3*rho_a^{2/3})
        let deps_slater_drho_a = -two_13 * SlaterExchange::C_X / (3.0 * rho_a_13 * rho_a_13);

        // d(d_eps_a)/d(rho_a) at fixed sigma_a:
        // d_eps_a = -beta * rho_a^{1/3} * g(x)
        // x = sqrt(sigma_a)/rho_a^{4/3}
        // dx/d(rho_a) = -4/3 * sqrt(sigma_a) / rho_a^{7/3} = -4/3 * x / rho_a
        let dx_drho_a = -4.0 / 3.0 * x / rho_a;

        // d(d_eps_a)/d(rho_a) = -beta * [x2/(denom*3*rho_a^{2/3}) + rho_a^{1/3} * dg/dx * dx/drho_a]
        let ddeps_drho_a =
            -Self::BETA * (x2 / denom / (3.0 * rho_a_13 * rho_a_13) + rho_a_13 * dg_dx * dx_drho_a);

        let deps_x_drho_a = deps_slater_drho_a + ddeps_drho_a;

        // vrho = eps_x + rho * deps_x/drho_a * (1/2)
        let vrho = eps_x + rho * deps_x_drho_a * 0.5;

        // d(eps_x)/d(sigma_a) at fixed rho_a:
        // Only d_eps_a depends on sigma_a via x
        // dx/d(sigma_a) = 1/(2*sqrt(sigma_a)*rho_a^{4/3})
        let dx_dsigma_a = 1.0 / (2.0 * grad_rho_a * rho_a_43);

        // d(d_eps_a)/d(sigma_a) = -beta * rho_a^{1/3} * dg/dx * dx/dsigma_a
        let ddeps_dsigma_a = -Self::BETA * rho_a_13 * dg_dx * dx_dsigma_a;

        // vsigma = rho * ddeps_dsigma_a * (1/4)
        let vsigma = rho * ddeps_dsigma_a * 0.25;

        (eps_x, vrho, vsigma)
    }
}

impl Becke88Exchange {
    /// Evaluate Becke88 exchange second derivatives at a single grid point.
    ///
    /// Returns `(v2rho2, v2rhosigma, v2sigma2)` = second partial derivatives
    /// of F(rho, sigma) = rho * eps_x(rho, sigma).
    ///
    /// Follows the same spin-density formulation as `eval_gga()`:
    /// rho_a = rho/2, sigma_a = sigma/4, x = |grad rho_a| / rho_a^{4/3}.
    ///
    /// # Reference
    ///
    /// Becke, A. D. (1988). Phys. Rev. A, 38, 3098. Second derivatives
    /// obtained by analytical differentiation of the first-derivative
    /// expressions in `eval_gga()`.
    #[inline]
    pub fn eval_fxc_gga(&self, rho: f64, sigma: f64) -> (f64, f64, f64) {
        if rho <= DENSITY_THRESHOLD {
            return (0.0, 0.0, 0.0);
        }

        let rho_a = rho / 2.0;
        let sigma_a = sigma.max(0.0) / 4.0;
        let rho_a_13 = rho_a.cbrt();
        let rho_a_23 = rho_a_13 * rho_a_13;
        let rho_a_43 = rho_a * rho_a_13;

        let two_13: f64 = 2.0_f64.cbrt();

        let grad_rho_a = sigma_a.sqrt();

        if grad_rho_a < 1e-30 {
            // At sigma=0: B88 = Slater, second derivatives:
            // v2rho2 = Slater fxc = -(4/9)*C_x*rho^{-2/3}
            // v2rhosigma = 0, v2sigma2 = 0 (evaluated at sigma→0 limit below)
            let slater = SlaterExchange;
            let v2rho2 = slater.eval_fxc(rho);

            // vsigma at sigma=0: need limit of d(vsigma)/d(rho) and d(vsigma)/d(sigma)
            // vsigma = rho * d(d_eps_a)/d(sigma_a) * 1/4
            // d_eps_a = -beta * rho_a^{1/3} * x^2/denom
            // At sigma=0 (x=0): dg/dx at x=0 needs care
            // g(x) = x^2/denom, dg/dx|_{x=0} = 0 (since x^2 term), d²g/dx²|_{x=0} = 2
            // d(d_eps_a)/d(sigma_a) = -beta * rho_a^{1/3} * dg/dx * dx/dsigma_a
            // But at x=0, dg/dx = 0, so vsigma = 0.

            // For v2rhosigma and v2sigma2 at sigma=0, we need limits:
            // Since vsigma is zero at sigma=0 and continuous, d(vsigma)/drho|_{sigma=0} = ?
            // We need the limit as sigma→0. In practice these components are finite but
            // require careful L'Hopital analysis. Use the full expression with sigma=0
            // guard to avoid division by zero and just return the Slater v2rho2.
            //
            // At exactly sigma=0: vsigma depends on sigma through x, so:
            //   v2sigma2 = d(vsigma)/d(sigma) involves terms like d²g/dx² * (dx/dsigma)²
            //   which go as 1/(rho_a^{8/3} * 4 * sigma_a) → diverges
            //
            // But PySCF returns finite values. The clean approach: evaluate at a tiny
            // positive sigma to get the correct limits for the gradient-dependent terms.
            // However, this function is called FROM the eval_xc_second_deriv() which
            // passes the actual grid sigma values. If sigma is truly zero on the grid,
            // these terms are physically irrelevant (multiplied by zero in the contraction).
            return (v2rho2, 0.0, 0.0);
        }

        // Reduced density gradient
        let x = grad_rho_a / rho_a_43;
        let x2 = x * x;

        let asinh_x = x.asinh();
        let sqrt_x2p1 = (x2 + 1.0).sqrt();
        let d_asinh = 1.0 / sqrt_x2p1; // d(asinh(x))/dx

        let denom = 1.0 + 6.0 * Self::BETA * x * asinh_x;
        let d_denom_dx = 6.0 * Self::BETA * (asinh_x + x * d_asinh);
        let denom2 = denom * denom;

        // g(x) = x^2 / denom
        let g = x2 / denom;
        let dg_dx = (2.0 * x * denom - x2 * d_denom_dx) / denom2;

        // =====================================================================
        // Second derivative of g(x) = x^2/denom w.r.t. x
        // =====================================================================

        // d²(asinh(x))/dx² = -x / (x²+1)^{3/2}
        let d2_asinh = -x / (sqrt_x2p1 * sqrt_x2p1 * sqrt_x2p1);

        let d2_denom_dx2 = 6.0 * Self::BETA * (2.0 * d_asinh + x * d2_asinh);

        // d²g/dx² = d/dx[(2x*D - x²*D')/D²]
        // = [(2D + 2x*D' - 2x*D' - x²*D'')*D² - (2x*D - x²*D')*2D*D'] / D⁴
        // = [(2D - x²*D'')*D - 2*(2x*D - x²*D')*D'] / D³
        let numerator_d2g = (2.0 * denom - x2 * d2_denom_dx2) * denom
            - 2.0 * (2.0 * x * denom - x2 * d_denom_dx) * d_denom_dx;
        let d2g_dx2 = numerator_d2g / (denom2 * denom);

        // =====================================================================
        // Key intermediate derivatives in spin-density variables
        // =====================================================================

        // dx/drho_a = -(4/3) * x / rho_a
        let dx_drho_a = -4.0 / 3.0 * x / rho_a;

        // dx/dsigma_a = 1 / (2 * sqrt(sigma_a) * rho_a^{4/3})
        let dx_dsigma_a = 1.0 / (2.0 * grad_rho_a * rho_a_43);

        // eps_x = eps_slater_a + d_eps_a where d_eps_a = -beta * rho_a^{1/3} * g(x)

        // =====================================================================
        // vrho and vsigma (already computed in eval_gga, reproduced here)
        // =====================================================================

        let deps_slater_drho_a = -two_13 * SlaterExchange::C_X / (3.0 * rho_a_23);

        // d(d_eps_a)/d(rho_a) at fixed sigma_a
        let ddeps_drho_a = -Self::BETA * (g / (3.0 * rho_a_23) + rho_a_13 * dg_dx * dx_drho_a);

        let deps_x_drho_a = deps_slater_drho_a + ddeps_drho_a;

        // d(d_eps_a)/d(sigma_a) at fixed rho_a
        let ddeps_dsigma_a = -Self::BETA * rho_a_13 * dg_dx * dx_dsigma_a;

        // eps_x = eps_slater_a + d_eps_a
        let eps_slater_a = -two_13 * SlaterExchange::C_X * rho_a_13;
        let d_eps_a = -Self::BETA * rho_a_13 * g;
        let eps_x = eps_slater_a + d_eps_a;

        // vrho = eps_x + rho * deps_x_drho_a * 0.5
        let _vrho = eps_x + rho * deps_x_drho_a * 0.5;

        // vsigma = rho * ddeps_dsigma_a * 0.25
        let _vsigma = rho * ddeps_dsigma_a * 0.25;

        // =====================================================================
        // v2rho2 = d(vrho)/d(rho)
        // vrho = eps_x(rho_a, sigma_a) + rho * d(eps_x)/d(rho_a) * (1/2)
        // d(vrho)/d(rho) = d(eps_x)/d(rho_a)*1/2 + d(eps_x)/d(rho_a)*1/2
        //                + rho * d²(eps_x)/d(rho_a)² * (1/2) * (1/2)
        //                = d(eps_x)/d(rho_a) + rho/4 * d²(eps_x)/d(rho_a)²
        // =====================================================================

        // d²(eps_slater_a)/d(rho_a)² = 2^{1/3}*C_x/(9*rho_a^{5/3})  ... let me compute:
        // deps_slater/drho_a = -2^{1/3}*C_x / (3*rho_a^{2/3})
        // d²eps_slater/drho_a² = 2^{1/3}*C_x*2 / (9*rho_a^{5/3})
        let rho_a_53 = rho_a_43 * rho_a_13; // rho_a^{5/3}
        let d2eps_slater_drho_a2 = two_13 * SlaterExchange::C_X * 2.0 / (9.0 * rho_a_53);

        // d²(d_eps_a)/d(rho_a)²:
        // d_eps_a_drho_a = -beta * [g/(3*rho_a^{2/3}) + rho_a^{1/3} * dg/dx * dx/drho_a]
        //
        // Need d/drho_a of each term (at fixed sigma_a):
        //
        // Term A: -beta * g / (3*rho_a^{2/3})
        // d(A)/drho_a = -beta * [dg/dx*dx/drho_a / (3*rho_a^{2/3})
        //                       + g * (-2/3) / (3*rho_a^{5/3})]
        //             = -beta * [dg/dx * dx/drho_a / (3*rho_a^{2/3})
        //                       - 2*g / (9*rho_a^{5/3})]
        //
        // Term B: -beta * rho_a^{1/3} * dg/dx * dx/drho_a
        //       = -beta * rho_a^{1/3} * dg/dx * (-4/3*x/rho_a)
        //       = beta * (4/3) * rho_a^{1/3} * dg/dx * x / rho_a
        //       = beta * (4/3) * dg/dx * x / rho_a^{2/3}
        //
        // d(B)/drho_a = -beta * d/drho_a[rho_a^{1/3} * dg/dx * dx/drho_a]
        //   = -beta * [rho_a^{1/3}/(3*rho_a) * dg/dx * dx/drho_a      (from d(rho_a^{1/3})/drho_a)
        //            + rho_a^{1/3} * d²g/dx² * dx/drho_a * dx/drho_a   (chain rule on dg/dx)
        //            + rho_a^{1/3} * dg/dx * d²x/drho_a²]              (second deriv of x)
        //
        // d²x/drho_a² = (-4/3) * d/drho_a[x/rho_a]
        //             = (-4/3) * [dx/drho_a * rho_a - x] / rho_a²
        //             = (-4/3) * [(-4/3*x/rho_a)*rho_a - x] / rho_a²
        //             = (-4/3) * [-4x/3 - x] / rho_a²
        //             = (-4/3) * (-7x/3) / rho_a²
        //             = 28*x / (9*rho_a²)
        let d2x_drho_a2 = 28.0 * x / (9.0 * rho_a * rho_a);

        let d_term_a =
            -Self::BETA * (dg_dx * dx_drho_a / (3.0 * rho_a_23) - 2.0 * g / (9.0 * rho_a_53));

        let d_term_b = -Self::BETA
            * (rho_a_13 / (3.0 * rho_a) * dg_dx * dx_drho_a
                + rho_a_13 * d2g_dx2 * dx_drho_a * dx_drho_a
                + rho_a_13 * dg_dx * d2x_drho_a2);

        let d2eps_drho_a2 = d2eps_slater_drho_a2 + d_term_a + d_term_b;

        let v2rho2 = deps_x_drho_a + rho / 4.0 * d2eps_drho_a2;

        // =====================================================================
        // v2rhosigma = d(vrho)/d(sigma) = d(vsigma)/d(rho) (by symmetry of mixed partials)
        //
        // vrho = eps_x + rho * deps_x_drho_a / 2
        // d(vrho)/d(sigma) = d(eps_x)/d(sigma_a)*1/4 + rho/2 * d²(eps_x)/d(rho_a)d(sigma_a) * 1/4
        //                  = [d(eps_x)/d(sigma_a) + rho/2 * d²(eps_x)/(d(rho_a) d(sigma_a))] / 4
        // =====================================================================

        // d²(d_eps_a)/(drho_a dsigma_a):
        // d(d_eps_a)/d(sigma_a) = -beta * rho_a^{1/3} * dg/dx * dx/dsigma_a
        // d/drho_a of above at fixed sigma_a:
        //   = -beta * [rho_a^{1/3}/(3*rho_a) * dg/dx * dx/dsigma_a
        //            + rho_a^{1/3} * d²g/dx² * dx/drho_a * dx/dsigma_a
        //            + rho_a^{1/3} * dg/dx * d²x/(drho_a dsigma_a)]
        //
        // d²x/(drho_a dsigma_a) = d/drho_a[1/(2*sqrt(sigma_a)*rho_a^{4/3})]
        //                       = -4/(3*2*sqrt(sigma_a)*rho_a^{7/3})
        //                       = -4/(6*sqrt(sigma_a)*rho_a^{7/3})
        //                       = -(4/3) * dx_dsigma_a / rho_a
        let d2x_drho_a_dsigma_a = -4.0 / 3.0 * dx_dsigma_a / rho_a;

        let d2eps_drho_a_dsigma_a = -Self::BETA
            * (rho_a_13 / (3.0 * rho_a) * dg_dx * dx_dsigma_a
                + rho_a_13 * d2g_dx2 * dx_drho_a * dx_dsigma_a
                + rho_a_13 * dg_dx * d2x_drho_a_dsigma_a);

        // deps_x/dsigma_a = ddeps_dsigma_a (Slater has no sigma dependence)
        let v2rhosigma = (ddeps_dsigma_a + rho / 2.0 * d2eps_drho_a_dsigma_a) / 4.0;

        // =====================================================================
        // v2sigma2 = d(vsigma)/d(sigma)
        // vsigma = rho * d(eps_x)/d(sigma_a) * 1/4
        //        = rho * (-beta * rho_a^{1/3} * dg/dx * dx/dsigma_a) * 1/4
        // d(vsigma)/d(sigma) = rho/4 * d²(eps_x)/d(sigma_a)² * 1/4
        //                    = rho/16 * d²(eps_x)/d(sigma_a)²
        // =====================================================================

        // d²(d_eps_a)/d(sigma_a)² = -beta * rho_a^{1/3} * [d²g/dx² * (dx/dsigma_a)²
        //                                                   + dg/dx * d²x/dsigma_a²]
        //
        // d²x/dsigma_a² = d/dsigma_a[1/(2*sqrt(sigma_a)*rho_a^{4/3})]
        //               = -1/(4*sigma_a^{3/2}*rho_a^{4/3})
        let d2x_dsigma_a2 = -1.0 / (4.0 * sigma_a * grad_rho_a * rho_a_43);

        let d2eps_dsigma_a2 =
            -Self::BETA * rho_a_13 * (d2g_dx2 * dx_dsigma_a * dx_dsigma_a + dg_dx * d2x_dsigma_a2);

        let v2sigma2 = rho / 16.0 * d2eps_dsigma_a2;

        (v2rho2, v2rhosigma, v2sigma2)
    }
}

impl ExchangeCorrelation for Becke88Exchange {
    fn eval_xc(&self, rho: &[f64], exc: &mut [f64], vxc: &mut [f64]) {
        // LDA fallback: evaluate at zero gradient (= Slater)
        for i in 0..rho.len() {
            let (e, v, _) = self.eval_gga(rho[i], 0.0);
            exc[i] = e;
            vxc[i] = v;
        }
    }

    fn eval_xc_gga(
        &self,
        rho: &[f64],
        sigma: &[f64],
        exc: &mut [f64],
        vrho: &mut [f64],
        vsigma: &mut [f64],
    ) {
        for i in 0..rho.len() {
            let (e, vr, vs) = self.eval_gga(rho[i], sigma[i]);
            exc[i] = e;
            vrho[i] = vr;
            vsigma[i] = vs;
        }
    }

    fn eval_xc_second_deriv(
        &self,
        rho: &[f64],
        sigma: &[f64],
        v2rho2: &mut [f64],
        v2rhosigma: &mut [f64],
        v2sigma2: &mut [f64],
    ) {
        for i in 0..rho.len() {
            let (vr2, vrs, vs2) = self.eval_fxc_gga(rho[i], sigma[i]);
            v2rho2[i] = vr2;
            v2rhosigma[i] = vrs;
            v2sigma2[i] = vs2;
        }
    }

    fn needs_gradient(&self) -> bool {
        true
    }

    fn hf_exchange_fraction(&self) -> f64 {
        0.0
    }

    fn name(&self) -> &str {
        "Becke88 Exchange"
    }
}

// =============================================================================
// LYP Correlation (GGA)
// =============================================================================

/// Lee-Yang-Parr correlation functional.
///
/// The LYP functional provides a gradient-corrected correlation energy.
/// For closed-shell systems (rho_alpha = rho_beta = rho/2), the functional
/// form simplifies significantly.
///
/// # Parameters (Lee, Yang, Parr 1988)
///
/// - a = 0.04918
/// - b = 0.132
/// - c = 0.2533
/// - d = 0.349
///
/// # Closed-Shell Formula
///
/// ```text
/// E_c^LYP = -a * integral { rho / (1 + d*rho^{-1/3})
///           + a*b * omega(rho) * [rho^2 * C_F * 2^{2/3} * rho^{-5/3}
///             - (11/24) * sigma
///             + (delta/8 + gamma/4) * sigma ] } dr
/// ```
///
/// where:
/// ```text
/// omega = exp(-c * rho^{-1/3}) * rho^{-11/3} / (1 + d*rho^{-1/3})
/// delta = c*rho^{-1/3} + d*rho^{-1/3}/(1 + d*rho^{-1/3})
/// C_F = (3/10)(3*pi^2)^{2/3}
/// ```
///
/// For the closed-shell case, gamma_ab = 0 (same-spin kinetic energy
/// correction vanishes for paired electrons), and the sigma terms
/// combine as shown.
///
/// # References
///
/// - Lee, C., Yang, W., Parr, R. G. (1988). Phys. Rev. B, 37, 785.
/// - Miehlich, B., Savin, A., Stoll, H., Preuss, H. (1989). CPL 157, 200.
///   (Analytic derivatives, closed-shell formula)
/// - PySCF: `dft.libxc.eval_xc('gga_c_lyp', ...)`
pub struct LypCorrelation;

impl LypCorrelation {
    // LYP parameters from Lee, Yang, Parr (1988)
    const A: f64 = 0.04918;
    const B: f64 = 0.132;
    const C: f64 = 0.2533;
    const D: f64 = 0.349;

    /// Fermi constant: C_F = (3/10)(3*pi^2)^{2/3}
    #[allow(clippy::excessive_precision)]
    const C_F: f64 = 2.8712340001881918;

    /// Evaluate LYP correlation at a single grid point (closed-shell).
    ///
    /// Returns `(eps_c, vrho, vsigma)`.
    ///
    /// # Formula
    ///
    /// The energy per unit volume F(rho, sigma) = rho * eps_c decomposes as:
    ///
    /// ```text
    /// F = -a * rho / (1+d*rho^{-1/3})                          [term 1]
    ///     - a*b * CF * omega * rho^{14/3}                       [term 2: TF kinetic]
    ///     + a*b * omega * rho^2 * (7*delta+3)/72 * sigma        [term 3: gradient]
    /// ```
    ///
    /// Verified against PySCF/libxc `gga_c_lyp` to machine precision.
    ///
    /// # References
    ///
    /// - Lee, Yang, Parr (1988). Phys. Rev. B, 37, 785.
    /// - Miehlich, Savin, Stoll, Preuss (1989). CPL 157, 200.
    #[inline]
    pub fn eval_gga(&self, rho: f64, sigma: f64) -> (f64, f64, f64) {
        if rho <= DENSITY_THRESHOLD {
            return (0.0, 0.0, 0.0);
        }

        let sigma_safe = sigma.max(0.0);

        // Build all needed powers of rho from a single cbrt() call.
        // This eliminates 3 powf() calls (rho^{-11/3}, rho^{14/3}, rho^{11/3})
        // which each cost ~50ns (log + multiply + exp), replaced by ~1ns multiplies.
        let rho_13 = rho.cbrt(); // rho^{1/3}
        let rho_23 = rho_13 * rho_13; // rho^{2/3}
        let rho_m13 = 1.0 / rho_13; // rho^{-1/3}
        let rho2 = rho * rho; // rho^2
        let rho3 = rho2 * rho; // rho^3
        let rho_m113 = 1.0 / (rho3 * rho_23); // rho^{-11/3} = 1/(rho^3 * rho^{2/3})

        let d_rho_m13 = Self::D * rho_m13;
        let inv_1pd = 1.0 / (1.0 + d_rho_m13);

        let c_rho_m13 = Self::C * rho_m13;
        let exp_term = (-c_rho_m13).exp();
        let omega = exp_term * inv_1pd * rho_m113;

        let delta = c_rho_m13 + d_rho_m13 * inv_1pd;

        // =====================================================================
        // F(rho, sigma) = term1 + term2 + term3
        // =====================================================================

        // Term 1: base correlation
        let term1 = -Self::A * rho * inv_1pd;

        // Term 2: Thomas-Fermi kinetic energy contribution
        // rho^{14/3} = rho^4 * rho^{2/3}
        let rho4 = rho2 * rho2;
        let rho_143 = rho4 * rho_23;
        let term2 = -Self::A * Self::B * Self::C_F * omega * rho_143;

        // Term 3: gradient-dependent term (linear in sigma)
        // B(rho) = a*b*omega*rho^2*(7*delta+3)/72
        let rho2 = rho * rho;
        let b_coeff = Self::A * Self::B * omega * rho2 * (7.0 * delta + 3.0) / 72.0;
        let term3 = b_coeff * sigma_safe;

        let f_total = term1 + term2 + term3;
        let eps_c = f_total / rho;

        // =====================================================================
        // Derivatives: vrho = dF/drho, vsigma = dF/dsigma
        // =====================================================================

        // vsigma = B(rho)
        let vsigma = b_coeff;

        // --- d(term1)/drho ---
        let h = 1.0 + d_rho_m13;
        let dterm1 = -Self::A * (1.0 + 4.0 / 3.0 * d_rho_m13) / (h * h);

        // --- d(omega)/drho ---
        // omega = exp(-c*rho^{-1/3}) * inv_1pd * rho^{-11/3}
        // d(omega)/drho = omega * [rho^{-1/3}/(3*rho)*(c + d*inv_1pd) - 11/(3*rho)]
        let domega =
            omega * (rho_m13 / (3.0 * rho) * (Self::C + Self::D * inv_1pd) - 11.0 / (3.0 * rho));

        // --- d(term2)/drho ---
        // term2 = -a*b*CF * omega * rho^{14/3}
        let rho_113 = rho3 * rho_23; // rho^{11/3} = rho^3 * rho^{2/3}
        let dterm2 =
            -Self::A * Self::B * Self::C_F * (domega * rho_143 + omega * 14.0 / 3.0 * rho_113);

        // --- d(delta)/drho ---
        let rho_43 = rho * rho_13;
        let ddelta = -Self::C * rho_m13 / (3.0 * rho)
            + Self::D * inv_1pd / (3.0 * rho_43) * (-1.0 + d_rho_m13 * inv_1pd);

        // --- d(term3)/drho ---
        // term3 = a*b * Q * sigma, where Q = omega*rho^2*(7*delta+3)/72
        // dQ/drho = [(domega*rho^2 + omega*2*rho)*(7*delta+3) + omega*rho^2*7*ddelta] / 72
        let dq = ((domega * rho2 + omega * 2.0 * rho) * (7.0 * delta + 3.0)
            + omega * rho2 * 7.0 * ddelta)
            / 72.0;
        let dterm3 = Self::A * Self::B * dq * sigma_safe;

        let vrho = dterm1 + dterm2 + dterm3;

        (eps_c, vrho, vsigma)
    }
}

impl LypCorrelation {
    /// Evaluate LYP correlation second derivatives at a single grid point.
    ///
    /// Returns `(v2rho2, v2rhosigma, v2sigma2)`.
    ///
    /// For closed-shell LYP, the energy per unit volume F(rho, sigma) is
    /// LINEAR in sigma, therefore `v2sigma2 = 0` identically.
    ///
    /// ```text
    /// F = term1(rho) + term2(rho) + B(rho) * sigma
    /// ```
    ///
    /// where B(rho) = a*b*omega*rho^2*(7*delta+3)/72.
    ///
    /// So:
    /// - vrho = d(term1)/drho + d(term2)/drho + dB/drho * sigma
    /// - vsigma = B(rho)                        [linear in sigma]
    /// - v2rho2 = d²F/drho² = d²(term1+term2)/drho² + d²B/drho² * sigma
    /// - v2rhosigma = dB/drho                    [cross derivative]
    /// - v2sigma2 = 0                            [linear in sigma]
    ///
    /// # Reference
    ///
    /// Lee, Yang, Parr (1988). Phys. Rev. B, 37, 785.
    /// Second derivatives by analytical differentiation.
    #[inline]
    pub fn eval_fxc_gga(&self, rho: f64, sigma: f64) -> (f64, f64, f64) {
        if rho <= DENSITY_THRESHOLD {
            return (0.0, 0.0, 0.0);
        }

        let sigma_safe = sigma.max(0.0);

        // Build powers of rho (same as eval_gga)
        let rho_13 = rho.cbrt();
        let rho_23 = rho_13 * rho_13;
        let rho_m13 = 1.0 / rho_13;
        let rho2 = rho * rho;
        let rho3 = rho2 * rho;
        let rho_m113 = 1.0 / (rho3 * rho_23); // rho^{-11/3}

        let d_rho_m13 = Self::D * rho_m13;
        let inv_1pd = 1.0 / (1.0 + d_rho_m13);

        let c_rho_m13 = Self::C * rho_m13;
        let exp_term = (-c_rho_m13).exp();
        let omega = exp_term * inv_1pd * rho_m113;

        let delta = c_rho_m13 + d_rho_m13 * inv_1pd;

        let rho4 = rho2 * rho2;
        let rho_143 = rho4 * rho_23; // rho^{14/3}
        let rho_113 = rho3 * rho_23; // rho^{11/3}
        let rho_43 = rho * rho_13; // rho^{4/3}

        let h = 1.0 + d_rho_m13; // = 1/inv_1pd

        // =====================================================================
        // d(omega)/d(rho) — same as eval_gga
        // =====================================================================
        let domega =
            omega * (rho_m13 / (3.0 * rho) * (Self::C + Self::D * inv_1pd) - 11.0 / (3.0 * rho));

        // =====================================================================
        // d(delta)/d(rho) — same as eval_gga
        // =====================================================================
        let ddelta = -Self::C * rho_m13 / (3.0 * rho)
            + Self::D * inv_1pd / (3.0 * rho_43) * (-1.0 + d_rho_m13 * inv_1pd);

        // =====================================================================
        // Term1: F1 = -a * rho / (1+d*rho^{-1/3})
        // dF1/drho and d²F1/drho²
        // =====================================================================
        // dterm1 is computed implicitly through d2term1 below

        // d²(term1)/drho²:
        // dterm1 = -A * (1 + 4/3*d*rho^{-1/3}) / h²
        // d/drho[(1 + 4/3*d*rho^{-1/3})] = -4/9 * d * rho^{-4/3}
        // d/drho[1/h²] = -2/h³ * dh/drho where dh/drho = -d/(3*rho^{4/3})
        //              = 2d/(3*rho^{4/3}*h³)
        let dh_drho = -Self::D / (3.0 * rho_43);
        let numer1 = 1.0 + 4.0 / 3.0 * d_rho_m13;
        let d_numer1 = -4.0 / 9.0 * Self::D / (rho_43);
        let d2term1 = -Self::A * (d_numer1 / (h * h) + numer1 * (-2.0) * dh_drho / (h * h * h));

        // =====================================================================
        // Term2: F2 = -a*b*CF*omega*rho^{14/3}
        // dF2/drho and d²F2/drho²
        // =====================================================================
        // dterm2 is computed implicitly through d2term2 below

        // d²(omega)/drho²:
        // domega = omega * f(rho) where f = (rho^{-1/3}/(3rho))*(C+D*inv_1pd) - 11/(3rho)
        // d(domega)/drho = domega*f + omega*df
        // Let's denote f_omega = domega/omega for clarity:
        let c_plus_d_inv = Self::C + Self::D * inv_1pd;
        let f_omega = rho_m13 / (3.0 * rho) * c_plus_d_inv - 11.0 / (3.0 * rho);

        // df/drho:
        // f = rho^{-4/3}/(3) * c_plus_d_inv - 11/(3*rho)
        // d/drho[rho^{-4/3}/3 * c_plus_d_inv]:
        //   = (-4/3)*rho^{-7/3}/3 * c_plus_d_inv + rho^{-4/3}/3 * d(c_plus_d_inv)/drho
        //   = -4/(9*rho^{7/3}) * c_plus_d_inv + rho^{-4/3}/3 * D*(-1)*inv_1pd²*dh_drho
        //   (recall d(inv_1pd)/drho = -inv_1pd² * dh_drho = -inv_1pd² * (-D/(3*rho^{4/3}))
        //                           = D*inv_1pd² / (3*rho^{4/3}))
        let rho_73 = rho_43 * rho; // rho^{7/3}
                                   // d(inv_1pd)/drho = -inv_1pd² * dh/drho = -inv_1pd² * (-D/(3*rho^{4/3}))
                                   //                 = D*inv_1pd²/(3*rho^{4/3})
        let d_inv_1pd_drho = Self::D * inv_1pd * inv_1pd / (3.0 * rho_43);
        let d_c_plus_d_inv_drho = Self::D * d_inv_1pd_drho;

        let rho_m43 = 1.0 / rho_43;
        let df_omega = -4.0 / (9.0 * rho_73) * c_plus_d_inv
            + rho_m43 / 3.0 * d_c_plus_d_inv_drho
            + 11.0 / (3.0 * rho * rho);

        let d2omega = domega * f_omega + omega * df_omega;

        // d²(term2)/drho²:
        // dterm2 = -A*B*CF * [domega * rho^{14/3} + omega * 14/3 * rho^{11/3}]
        // d²term2 = -A*B*CF * [d2omega * rho^{14/3} + domega * 14/3 * rho^{11/3}
        //                     + domega * 14/3 * rho^{11/3} + omega * 14/3 * 11/3 * rho^{8/3}]
        //         = -A*B*CF * [d2omega * rho^{14/3} + 2*domega * 14/3 * rho^{11/3}
        //                     + omega * 154/9 * rho^{8/3}]
        let rho_83 = rho2 * rho_23;
        let d2term2 = -Self::A
            * Self::B
            * Self::C_F
            * (d2omega * rho_143
                + 2.0 * domega * 14.0 / 3.0 * rho_113
                + omega * 154.0 / 9.0 * rho_83);

        // =====================================================================
        // Term3: F3 = A*B * Q(rho) * sigma where Q = omega*rho^2*(7*delta+3)/72
        // B_coeff = A*B*Q in the notation of eval_gga
        //
        // vsigma = B_coeff = A*B*Q
        // vrho_from_term3 = dB_coeff/drho * sigma = A*B * dQ/drho * sigma
        // v2rho2_from_term3 = A*B * d²Q/drho² * sigma
        // v2rhosigma = A*B * dQ/drho
        // v2sigma2 = 0 (linear in sigma)
        // =====================================================================

        // Q = omega * rho^2 * (7*delta+3) / 72
        let seven_delta_3 = 7.0 * delta + 3.0;

        // dQ/drho (already computed in eval_gga):
        let dq = ((domega * rho2 + omega * 2.0 * rho) * seven_delta_3
            + omega * rho2 * 7.0 * ddelta)
            / 72.0;

        // d²Q/drho²:
        // Q = omega*rho²*(7delta+3)/72
        // dQ = [(domega*rho² + omega*2rho)*(7delta+3) + omega*rho²*7*ddelta] / 72
        //
        // Need d²omega, d(ddelta), and chain-rule products.

        // d²(delta)/drho²:
        // delta = c*rho^{-1/3} + d*rho^{-1/3}*inv_1pd
        // ddelta = -c*rho^{-1/3}/(3*rho) + D*inv_1pd/(3*rho^{4/3})*(-1 + d*rho^{-1/3}*inv_1pd)
        //
        // Let me compute d²delta by differentiating ddelta:
        // ddelta = -C/(3*rho^{4/3}) + D*inv_1pd*(d*rho^{-1/3}*inv_1pd - 1) / (3*rho^{4/3})
        //
        // Split: ddelta = P1 + P2
        // P1 = -C/(3*rho^{4/3})
        // dP1/drho = C*4/(9*rho^{7/3})
        let dp1 = Self::C * 4.0 / (9.0 * rho_73);

        // P2 = D*inv_1pd/(3*rho^{4/3}) * (d*rho^{-1/3}*inv_1pd - 1)
        //    = D*inv_1pd * (d_rho_m13*inv_1pd - 1) / (3*rho^{4/3})
        //
        // Let u = D*inv_1pd, v = d_rho_m13*inv_1pd - 1, w = 1/(3*rho^{4/3})
        // P2 = u * v * w
        // dP2/drho = du*v*w + u*dv*w + u*v*dw
        let u = Self::D * inv_1pd;
        let v = d_rho_m13 * inv_1pd - 1.0;
        let w = 1.0 / (3.0 * rho_43);

        let du = Self::D * d_inv_1pd_drho;
        // dv = d(d_rho_m13*inv_1pd)/drho - 0
        //    = D*(-1/3)*rho^{-4/3} * inv_1pd + d_rho_m13 * d_inv_1pd_drho
        let dv = Self::D * (-1.0 / 3.0) * rho_m43 * inv_1pd + d_rho_m13 * d_inv_1pd_drho;
        let dw = -4.0 / (9.0 * rho_73);

        let dp2 = du * v * w + u * dv * w + u * v * dw;

        let d2delta = dp1 + dp2;

        // Now d²Q/drho²:
        // dQ = [(domega*rho² + omega*2rho)*(7delta+3) + omega*rho²*7*ddelta] / 72
        //    = [A_term * B_term + C_term] / 72
        // where A_term = domega*rho² + omega*2rho, B_term = 7delta+3, C_term = omega*rho²*7*ddelta
        //
        // d²Q = d/drho[A_term * B_term + C_term] / 72
        // d(A_term) = d2omega*rho² + domega*2rho + domega*2rho + omega*2
        //           = d2omega*rho² + 4*domega*rho + 2*omega
        let a_term = domega * rho2 + omega * 2.0 * rho;
        let b_term = seven_delta_3;

        let da_term = d2omega * rho2 + 4.0 * domega * rho + 2.0 * omega;
        let db_term = 7.0 * ddelta;
        let dc_term = domega * rho2 * 7.0 * ddelta
            + omega * 2.0 * rho * 7.0 * ddelta
            + omega * rho2 * 7.0 * d2delta;

        let d2q = (da_term * b_term + a_term * db_term + dc_term) / 72.0;

        // =====================================================================
        // Assemble second derivatives
        // =====================================================================

        let v2rho2 = d2term1 + d2term2 + Self::A * Self::B * d2q * sigma_safe;
        let v2rhosigma = Self::A * Self::B * dq;
        // v2sigma2 = 0 (LYP is linear in sigma for closed-shell)

        (v2rho2, v2rhosigma, 0.0)
    }
}

impl ExchangeCorrelation for LypCorrelation {
    fn eval_xc(&self, rho: &[f64], exc: &mut [f64], vxc: &mut [f64]) {
        for i in 0..rho.len() {
            let (e, v, _) = self.eval_gga(rho[i], 0.0);
            exc[i] = e;
            vxc[i] = v;
        }
    }

    fn eval_xc_gga(
        &self,
        rho: &[f64],
        sigma: &[f64],
        exc: &mut [f64],
        vrho: &mut [f64],
        vsigma: &mut [f64],
    ) {
        for i in 0..rho.len() {
            let (e, vr, vs) = self.eval_gga(rho[i], sigma[i]);
            exc[i] = e;
            vrho[i] = vr;
            vsigma[i] = vs;
        }
    }

    fn eval_xc_second_deriv(
        &self,
        rho: &[f64],
        sigma: &[f64],
        v2rho2: &mut [f64],
        v2rhosigma: &mut [f64],
        v2sigma2: &mut [f64],
    ) {
        for i in 0..rho.len() {
            let (vr2, vrs, vs2) = self.eval_fxc_gga(rho[i], sigma[i]);
            v2rho2[i] = vr2;
            v2rhosigma[i] = vrs;
            v2sigma2[i] = vs2;
        }
    }

    fn needs_gradient(&self) -> bool {
        true
    }

    fn hf_exchange_fraction(&self) -> f64 {
        0.0
    }

    fn name(&self) -> &str {
        "LYP Correlation"
    }
}

// =============================================================================
// B3LYP Composite Functional
// =============================================================================

/// B3LYP hybrid GGA functional (VWN5 variant).
///
/// The B3LYP functional uses three empirical parameters to mix:
/// - Slater local exchange (LDA)
/// - Becke88 gradient-corrected exchange (GGA)
/// - VWN5 local correlation (LDA)
/// - LYP gradient-corrected correlation (GGA)
/// - Hartree-Fock exact exchange
///
/// The total xc energy is:
///
/// ```text
/// E_xc = (1-a)*E_x^Slater + a*E_x^HF + b*dE_x^B88 + (1-c)*E_c^VWN5 + c*E_c^LYP
/// ```
///
/// with a = 0.20, b = 0.72, c = 0.81.
///
/// Since `gga_x_b88` includes Slater exchange, the DFT-only part is:
///
/// ```text
/// eps_xc^DFT = (1-a-b)*Slater + b*B88_full + (1-c)*VWN5 + c*LYP
///            = 0.08*Slater + 0.72*B88_full + 0.19*VWN5 + 0.81*LYP
/// ```
///
/// which is equivalent to:
///
/// ```text
/// eps_xc^DFT = (1-a)*Slater + b*(B88_full - Slater) + (1-c)*VWN5 + c*LYP
///            = 0.80*Slater + 0.72*dB88 + 0.19*VWN5 + 0.81*LYP
/// ```
///
/// The HF exchange fraction (a = 0.20) is applied in the SCF loop
/// using the existing K matrix infrastructure.
///
/// # VWN5 Variant
///
/// This uses VWN5 (parametrization V from Vosko, Wilk, Nusair 1980),
/// NOT VWN3 (RPA). In PySCF this corresponds to `xc='b3lyp5'`.
///
/// # References
///
/// - Stephens, P. J., et al. (1994). J. Phys. Chem. 98, 11623.
/// - Becke, A. D. (1993). J. Chem. Phys. 98, 5648.
pub struct B3lyp {
    slater: SlaterExchange,
    becke88: Becke88Exchange,
    vwn5: Vwn5Correlation,
    lyp: LypCorrelation,
}

impl B3lyp {
    /// B3LYP mixing parameter a (HF exchange fraction).
    pub const A: f64 = 0.20;
    /// B3LYP mixing parameter b (Becke88 gradient correction weight).
    pub const B_PARAM: f64 = 0.72;
    /// B3LYP mixing parameter c (LYP correlation weight).
    pub const C_PARAM: f64 = 0.81;

    /// Slater coefficient: (1-a-b) = 0.08
    const SLATER_COEFF: f64 = 1.0 - Self::A - Self::B_PARAM; // 0.08

    /// Becke88 (full, including Slater) coefficient: b = 0.72
    const B88_COEFF: f64 = Self::B_PARAM; // 0.72

    /// VWN5 coefficient: (1-c) = 0.19
    const VWN5_COEFF: f64 = 1.0 - Self::C_PARAM; // 0.19

    /// LYP coefficient: c = 0.81
    const LYP_COEFF: f64 = Self::C_PARAM; // 0.81

    /// Create a new B3LYP functional.
    pub fn new() -> Self {
        Self {
            slater: SlaterExchange,
            becke88: Becke88Exchange,
            vwn5: Vwn5Correlation,
            lyp: LypCorrelation,
        }
    }
}

impl Default for B3lyp {
    fn default() -> Self {
        Self::new()
    }
}

impl ExchangeCorrelation for B3lyp {
    fn eval_xc(&self, rho: &[f64], exc: &mut [f64], vxc: &mut [f64]) {
        // LDA evaluation: B88 at sigma=0 = Slater, LYP at sigma=0
        for i in 0..rho.len() {
            let (eps_s, v_s) = self.slater.eval(rho[i]);
            let (eps_v, v_v) = self.vwn5.eval(rho[i]);
            let (eps_l, v_l, _) = self.lyp.eval_gga(rho[i], 0.0);

            // At sigma=0: B88 = Slater, so:
            // 0.08*Slater + 0.72*Slater + 0.19*VWN5 + 0.81*LYP
            // = 0.80*Slater + 0.19*VWN5 + 0.81*LYP
            exc[i] = (Self::SLATER_COEFF + Self::B88_COEFF) * eps_s
                + Self::VWN5_COEFF * eps_v
                + Self::LYP_COEFF * eps_l;
            vxc[i] = (Self::SLATER_COEFF + Self::B88_COEFF) * v_s
                + Self::VWN5_COEFF * v_v
                + Self::LYP_COEFF * v_l;
        }
    }

    fn eval_xc_gga(
        &self,
        rho: &[f64],
        sigma: &[f64],
        exc: &mut [f64],
        vrho: &mut [f64],
        vsigma: &mut [f64],
    ) {
        // B3LYP DFT part: 0.08*Slater + 0.72*B88_full + 0.19*VWN5 + 0.81*LYP
        for i in 0..rho.len() {
            let (eps_s, v_s) = self.slater.eval(rho[i]);
            let (eps_v, v_v) = self.vwn5.eval(rho[i]);
            let (eps_b88, vr_b88, vs_b88) = self.becke88.eval_gga(rho[i], sigma[i]);
            let (eps_lyp, vr_lyp, vs_lyp) = self.lyp.eval_gga(rho[i], sigma[i]);

            exc[i] = Self::SLATER_COEFF * eps_s
                + Self::B88_COEFF * eps_b88
                + Self::VWN5_COEFF * eps_v
                + Self::LYP_COEFF * eps_lyp;

            vrho[i] = Self::SLATER_COEFF * v_s
                + Self::B88_COEFF * vr_b88
                + Self::VWN5_COEFF * v_v
                + Self::LYP_COEFF * vr_lyp;

            vsigma[i] = Self::B88_COEFF * vs_b88 + Self::LYP_COEFF * vs_lyp;
        }
    }

    fn eval_xc_second_deriv(
        &self,
        rho: &[f64],
        sigma: &[f64],
        v2rho2: &mut [f64],
        v2rhosigma: &mut [f64],
        v2sigma2: &mut [f64],
    ) {
        // B3LYP second derivatives:
        // 0.08*Slater_fxc + 0.72*B88_fxc + 0.19*VWN5_fxc + 0.81*LYP_fxc
        for i in 0..rho.len() {
            let slater_fxc = self.slater.eval_fxc(rho[i]);
            let vwn5_fxc = self.vwn5.eval_fxc(rho[i]);
            let (b88_v2rho2, b88_v2rhosigma, b88_v2sigma2) =
                self.becke88.eval_fxc_gga(rho[i], sigma[i]);
            let (lyp_v2rho2, lyp_v2rhosigma, lyp_v2sigma2) =
                self.lyp.eval_fxc_gga(rho[i], sigma[i]);

            v2rho2[i] = Self::SLATER_COEFF * slater_fxc
                + Self::B88_COEFF * b88_v2rho2
                + Self::VWN5_COEFF * vwn5_fxc
                + Self::LYP_COEFF * lyp_v2rho2;

            v2rhosigma[i] = Self::B88_COEFF * b88_v2rhosigma + Self::LYP_COEFF * lyp_v2rhosigma;

            v2sigma2[i] = Self::B88_COEFF * b88_v2sigma2 + Self::LYP_COEFF * lyp_v2sigma2;
        }
    }

    fn needs_gradient(&self) -> bool {
        true
    }

    fn hf_exchange_fraction(&self) -> f64 {
        Self::A // 0.20
    }

    fn name(&self) -> &str {
        "B3LYP"
    }
}

// =============================================================================
// DftFunctional Enum
// =============================================================================

/// Available DFT functionals.
///
/// Used for method selection in the KS-SCF loop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DftFunctional {
    /// Local Density Approximation: Slater exchange + VWN5 correlation.
    Lda,
    /// B3LYP hybrid GGA: Becke88 + LYP + VWN5 + 20% HF exchange.
    B3lyp,
}

impl DftFunctional {
    /// Human-readable name of the functional.
    pub fn name(&self) -> &str {
        match self {
            Self::Lda => "LDA (Slater + VWN5)",
            Self::B3lyp => "B3LYP",
        }
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // PySCF Golden Reference Data (PySCF 2.11.0, dft.libxc.eval_xc)
    // =========================================================================

    const TEST_DENSITIES: [f64; 6] = [0.01, 0.10, 0.50, 1.00, 5.00, 10.00];

    // Slater exchange: dft.libxc.eval_xc('lda_x', rho)
    const PYSCF_EPS_X: [f64; 6] = [
        -1.591176626920582e-01,
        -3.428086123005624e-01,
        -5.861944813475790e-01,
        -7.385587663820223e-01,
        -1.262917725720472e+00,
        -1.591176626920582e+00,
    ];
    const PYSCF_V_X: [f64; 6] = [
        -2.121568835894110e-01,
        -4.570781497340832e-01,
        -7.815926417967720e-01,
        -9.847450218426964e-01,
        -1.683890300960629e+00,
        -2.121568835894110e+00,
    ];

    // VWN5 correlation: dft.libxc.eval_xc('lda_c_vwn', rho)
    const PYSCF_EPS_C: [f64; 6] = [
        -3.764519026217142e-02,
        -5.339728918594981e-02,
        -6.589401996674182e-02,
        -7.159261230679065e-02,
        -8.544337816136141e-02,
        -9.163970578243799e-02,
    ];
    const PYSCF_V_C: [f64; 6] = [
        -4.387265644739365e-02,
        -6.081203033126155e-02,
        -7.398704765154965e-02,
        -7.993838317598562e-02,
        -9.429029059231650e-02,
        -1.006684090462797e-01,
    ];

    // LDA composite: dft.libxc.eval_xc('lda,vwn5', rho)
    const PYSCF_EPS_XC: [f64; 6] = [
        -1.967628529542297e-01,
        -3.962059014865122e-01,
        -6.520885013143208e-01,
        -8.101513786888130e-01,
        -1.348361103881834e+00,
        -1.682816332703021e+00,
    ];
    const PYSCF_V_XC: [f64; 6] = [
        -2.560295400368047e-01,
        -5.178901800653447e-01,
        -8.555796894483216e-01,
        -1.064683405018682e+00,
        -1.778180591552946e+00,
        -2.222237244940390e+00,
    ];

    // =========================================================================
    // Slater Exchange Tests
    // =========================================================================

    #[test]
    fn slater_cx_constant() {
        // C_x = (3/4)(3/pi)^{1/3}
        let cx_computed = 0.75 * (3.0 / PI).cbrt();
        assert!(
            (SlaterExchange::C_X - cx_computed).abs() < 1e-15,
            "C_x mismatch: {} vs {}",
            SlaterExchange::C_X,
            cx_computed
        );
    }

    #[test]
    fn slater_golden_values() {
        // SX-6: All 6 test densities match PySCF to 1e-14
        let slater = SlaterExchange;
        for i in 0..TEST_DENSITIES.len() {
            let (eps_x, v_x) = slater.eval(TEST_DENSITIES[i]);
            assert!(
                (eps_x - PYSCF_EPS_X[i]).abs() < 1e-14,
                "eps_x mismatch at rho={}: {} vs {} (diff={:.2e})",
                TEST_DENSITIES[i],
                eps_x,
                PYSCF_EPS_X[i],
                (eps_x - PYSCF_EPS_X[i]).abs()
            );
            assert!(
                (v_x - PYSCF_V_X[i]).abs() < 1e-14,
                "v_x mismatch at rho={}: {} vs {} (diff={:.2e})",
                TEST_DENSITIES[i],
                v_x,
                PYSCF_V_X[i],
                (v_x - PYSCF_V_X[i]).abs()
            );
        }
    }

    #[test]
    fn slater_four_thirds_relation() {
        // SX-4: v_x = (4/3) * eps_x for all test densities
        let slater = SlaterExchange;
        for &rho in &TEST_DENSITIES {
            let (eps_x, v_x) = slater.eval(rho);
            let ratio = v_x / eps_x;
            assert!(
                (ratio - 4.0 / 3.0).abs() < 1e-14,
                "v_x/eps_x ratio mismatch at rho={}: {} (expected 4/3)",
                rho,
                ratio
            );
        }
    }

    #[test]
    fn slater_zero_density() {
        // SX-5: rho=0 returns (0, 0)
        let slater = SlaterExchange;
        let (eps_x, v_x) = slater.eval(0.0);
        assert_eq!(eps_x, 0.0);
        assert_eq!(v_x, 0.0);
    }

    #[test]
    fn slater_below_threshold() {
        // rho below DENSITY_THRESHOLD returns (0, 0)
        let slater = SlaterExchange;
        let (eps_x, v_x) = slater.eval(1e-30);
        assert_eq!(eps_x, 0.0);
        assert_eq!(v_x, 0.0);
    }

    #[test]
    fn slater_eps_x_always_negative() {
        // PROP-1: eps_x negative for all rho > 0
        let slater = SlaterExchange;
        let test_rhos = [1e-10, 1e-5, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0];
        for &rho in &test_rhos {
            let (eps_x, _) = slater.eval(rho);
            assert!(
                eps_x < 0.0,
                "eps_x should be negative at rho={}, got {}",
                rho,
                eps_x
            );
        }
    }

    #[test]
    fn slater_v_x_always_negative() {
        let slater = SlaterExchange;
        let test_rhos = [1e-10, 1e-5, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0];
        for &rho in &test_rhos {
            let (_, v_x) = slater.eval(rho);
            assert!(
                v_x < 0.0,
                "v_x should be negative at rho={}, got {}",
                rho,
                v_x
            );
        }
    }

    #[test]
    fn slater_monotonically_decreasing() {
        // PROP-3: eps_x decreases with increasing rho
        // (becomes more negative for higher density)
        let slater = SlaterExchange;
        let rhos = [0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 100.0];
        for w in rhos.windows(2) {
            let (eps1, _) = slater.eval(w[0]);
            let (eps2, _) = slater.eval(w[1]);
            assert!(
                eps2 < eps1,
                "eps_x should decrease: eps_x({})={} > eps_x({})={}",
                w[0],
                eps1,
                w[1],
                eps2
            );
        }
    }

    // =========================================================================
    // VWN5 Correlation Tests
    // =========================================================================

    #[test]
    fn vwn5_a_parameter() {
        // VWN-1: A = 0.0621814 / 2 = 0.0310907
        assert!(
            (Vwn5Correlation::A - 0.0621814 / 2.0).abs() < 1e-10,
            "A parameter mismatch"
        );
    }

    #[test]
    fn vwn5_q_constant() {
        // VWN-2: Q = sqrt(4c - b^2)
        let q_computed =
            (4.0 * Vwn5Correlation::C - Vwn5Correlation::B * Vwn5Correlation::B).sqrt();
        assert!(
            (Vwn5Correlation::Q - q_computed).abs() < 1e-12,
            "Q mismatch: {} vs {}",
            Vwn5Correlation::Q,
            q_computed
        );
    }

    #[test]
    fn vwn5_x_x0_constant() {
        // X(x0) = x0^2 + b*x0 + c
        let x_x0_computed = Vwn5Correlation::X0 * Vwn5Correlation::X0
            + Vwn5Correlation::B * Vwn5Correlation::X0
            + Vwn5Correlation::C;
        assert!(
            (Vwn5Correlation::X_X0 - x_x0_computed).abs() < 1e-10,
            "X(x0) mismatch: {} vs {}",
            Vwn5Correlation::X_X0,
            x_x0_computed
        );
    }

    #[test]
    fn vwn5_golden_values() {
        // VWN-6: All 6 test densities match PySCF to 1e-14
        let vwn5 = Vwn5Correlation;
        for i in 0..TEST_DENSITIES.len() {
            let (eps_c, v_c) = vwn5.eval(TEST_DENSITIES[i]);
            assert!(
                (eps_c - PYSCF_EPS_C[i]).abs() < 1e-14,
                "eps_c mismatch at rho={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                eps_c,
                PYSCF_EPS_C[i],
                (eps_c - PYSCF_EPS_C[i]).abs()
            );
            assert!(
                (v_c - PYSCF_V_C[i]).abs() < 1e-14,
                "v_c mismatch at rho={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                v_c,
                PYSCF_V_C[i],
                (v_c - PYSCF_V_C[i]).abs()
            );
        }
    }

    #[test]
    fn vwn5_zero_density() {
        // VWN-5: rho=0 returns (0, 0)
        let vwn5 = Vwn5Correlation;
        let (eps_c, v_c) = vwn5.eval(0.0);
        assert_eq!(eps_c, 0.0);
        assert_eq!(v_c, 0.0);
    }

    #[test]
    fn vwn5_below_threshold() {
        let vwn5 = Vwn5Correlation;
        let (eps_c, v_c) = vwn5.eval(1e-30);
        assert_eq!(eps_c, 0.0);
        assert_eq!(v_c, 0.0);
    }

    #[test]
    fn vwn5_eps_c_always_negative() {
        // PROP-2: eps_c negative for all rho > 0
        let vwn5 = Vwn5Correlation;
        let test_rhos = [1e-10, 1e-5, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0];
        for &rho in &test_rhos {
            let (eps_c, _) = vwn5.eval(rho);
            assert!(
                eps_c < 0.0,
                "eps_c should be negative at rho={}, got {}",
                rho,
                eps_c
            );
        }
    }

    #[test]
    fn vwn5_derivative_consistency() {
        // Numerical derivative check: v_c = d(rho*eps_c)/d(rho)
        // Use finite differences to verify the analytical potential.
        let vwn5 = Vwn5Correlation;
        let test_rhos = [0.01, 0.1, 0.5, 1.0, 5.0, 10.0];
        let h = 1e-7;

        for &rho in &test_rhos {
            let (eps_c, v_c) = vwn5.eval(rho);

            // Finite difference: d(rho*eps_c)/d(rho) ~ [(rho+h)*eps_c(rho+h) - (rho-h)*eps_c(rho-h)] / (2h)
            let (eps_plus, _) = vwn5.eval(rho + h);
            let (eps_minus, _) = vwn5.eval(rho - h);
            let v_c_numerical = ((rho + h) * eps_plus - (rho - h) * eps_minus) / (2.0 * h);

            let diff = (v_c - v_c_numerical).abs();
            assert!(
                diff < 1e-6,
                "VWN5 potential mismatch at rho={}: analytical={:.12e}, numerical={:.12e}, diff={:.2e}",
                rho,
                v_c,
                v_c_numerical,
                diff
            );

            // Also verify eps_c is used correctly
            let _ = eps_c; // suppress unused warning in the check above
        }
    }

    // =========================================================================
    // LDA Composite Tests
    // =========================================================================

    #[test]
    fn lda_golden_values() {
        // LDA-3: All 6 test densities match PySCF to 1e-14
        let lda = Lda::new();
        let mut exc = vec![0.0; TEST_DENSITIES.len()];
        let mut vxc = vec![0.0; TEST_DENSITIES.len()];
        lda.eval_xc(&TEST_DENSITIES, &mut exc, &mut vxc);

        for i in 0..TEST_DENSITIES.len() {
            assert!(
                (exc[i] - PYSCF_EPS_XC[i]).abs() < 1e-14,
                "eps_xc mismatch at rho={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                exc[i],
                PYSCF_EPS_XC[i],
                (exc[i] - PYSCF_EPS_XC[i]).abs()
            );
            assert!(
                (vxc[i] - PYSCF_V_XC[i]).abs() < 1e-14,
                "v_xc mismatch at rho={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                vxc[i],
                PYSCF_V_XC[i],
                (vxc[i] - PYSCF_V_XC[i]).abs()
            );
        }
    }

    #[test]
    fn lda_additivity() {
        // LDA-1 & LDA-2: eps_xc = eps_x + eps_c, v_xc = v_x + v_c
        let slater = SlaterExchange;
        let vwn5 = Vwn5Correlation;
        let lda = Lda::new();

        let mut exc = vec![0.0; TEST_DENSITIES.len()];
        let mut vxc = vec![0.0; TEST_DENSITIES.len()];
        lda.eval_xc(&TEST_DENSITIES, &mut exc, &mut vxc);

        for i in 0..TEST_DENSITIES.len() {
            let (ex, vx) = slater.eval(TEST_DENSITIES[i]);
            let (ec, vc) = vwn5.eval(TEST_DENSITIES[i]);

            assert!(
                (exc[i] - (ex + ec)).abs() < 1e-15,
                "eps_xc additivity failed at rho={}",
                TEST_DENSITIES[i]
            );
            assert!(
                (vxc[i] - (vx + vc)).abs() < 1e-15,
                "v_xc additivity failed at rho={}",
                TEST_DENSITIES[i]
            );
        }
    }

    #[test]
    fn lda_edge_cases() {
        // LDA-4: Edge case densities
        let lda = Lda::new();

        // Zero density
        let mut exc = [0.0];
        let mut vxc = [0.0];
        lda.eval_xc(&[0.0], &mut exc, &mut vxc);
        assert_eq!(exc[0], 0.0);
        assert_eq!(vxc[0], 0.0);

        // Below threshold
        lda.eval_xc(&[1e-30], &mut exc, &mut vxc);
        assert_eq!(exc[0], 0.0);
        assert_eq!(vxc[0], 0.0);

        // At threshold boundary (1e-20 should be zero per PySCF)
        lda.eval_xc(&[1e-20], &mut exc, &mut vxc);
        assert_eq!(exc[0], 0.0);
        assert_eq!(vxc[0], 0.0);

        // Just above threshold: rho=1e-10 (PySCF gives nonzero)
        let rho_small = [1e-10];
        let pyscf_eps_xc_small = -6.319450843668625e-04;
        let pyscf_v_xc_small = -8.391659194511679e-04;
        lda.eval_xc(&rho_small, &mut exc, &mut vxc);
        assert!(
            (exc[0] - pyscf_eps_xc_small).abs() < 1e-12,
            "eps_xc at rho=1e-10: {:.15e} vs {:.15e}",
            exc[0],
            pyscf_eps_xc_small
        );
        assert!(
            (vxc[0] - pyscf_v_xc_small).abs() < 1e-12,
            "v_xc at rho=1e-10: {:.15e} vs {:.15e}",
            vxc[0],
            pyscf_v_xc_small
        );

        // High densities
        let high_density_cases = [
            (1e-5, -2.472222743994018e-02, -3.226218538943055e-02),
            (100.0, -3.541100547469025e+00, -4.693303254169443e+00),
            (1000.0, -7.520891784779448e+00, -9.992585690290449e+00),
        ];
        for (rho, ref_eps, ref_v) in &high_density_cases {
            lda.eval_xc(&[*rho], &mut exc, &mut vxc);
            assert!(
                (exc[0] - ref_eps).abs() < 1e-12,
                "eps_xc at rho={}: {:.15e} vs {:.15e} (diff={:.2e})",
                rho,
                exc[0],
                ref_eps,
                (exc[0] - ref_eps).abs()
            );
            assert!(
                (vxc[0] - ref_v).abs() < 1e-12,
                "v_xc at rho={}: {:.15e} vs {:.15e} (diff={:.2e})",
                rho,
                vxc[0],
                ref_v,
                (vxc[0] - ref_v).abs()
            );
        }
    }

    // =========================================================================
    // Trait Interface Tests
    // =========================================================================

    #[test]
    fn trait_needs_gradient() {
        let slater = SlaterExchange;
        let vwn5 = Vwn5Correlation;
        let lda = Lda::new();

        assert!(!slater.needs_gradient());
        assert!(!vwn5.needs_gradient());
        assert!(!lda.needs_gradient());
    }

    #[test]
    fn trait_hf_exchange_fraction() {
        let slater = SlaterExchange;
        let vwn5 = Vwn5Correlation;
        let lda = Lda::new();

        assert_eq!(slater.hf_exchange_fraction(), 0.0);
        assert_eq!(vwn5.hf_exchange_fraction(), 0.0);
        assert_eq!(lda.hf_exchange_fraction(), 0.0);
    }

    #[test]
    fn trait_name() {
        let slater = SlaterExchange;
        let vwn5 = Vwn5Correlation;
        let lda = Lda::new();

        assert_eq!(slater.name(), "Slater Exchange");
        assert_eq!(vwn5.name(), "VWN5 Correlation");
        assert_eq!(lda.name(), "LDA (Slater + VWN5)");
    }

    #[test]
    fn dft_functional_enum() {
        assert_eq!(DftFunctional::Lda.name(), "LDA (Slater + VWN5)");
        assert_eq!(DftFunctional::B3lyp.name(), "B3LYP");
    }

    // =========================================================================
    // Slice Consistency Test
    // =========================================================================

    #[test]
    fn slice_matches_pointwise() {
        // PROP-4: eval_xc on slice == pointwise evaluation
        let lda = Lda::new();
        let rhos: Vec<f64> = (1..=100).map(|i| i as f64 * 0.1).collect();

        // Slice evaluation
        let mut exc_slice = vec![0.0; rhos.len()];
        let mut vxc_slice = vec![0.0; rhos.len()];
        lda.eval_xc(&rhos, &mut exc_slice, &mut vxc_slice);

        // Pointwise evaluation
        let slater = SlaterExchange;
        let vwn5 = Vwn5Correlation;
        for (i, &rho) in rhos.iter().enumerate() {
            let (ex, vx) = slater.eval(rho);
            let (ec, vc) = vwn5.eval(rho);
            assert!(
                (exc_slice[i] - (ex + ec)).abs() < 1e-15,
                "Slice/pointwise mismatch at index {}, rho={}",
                i,
                rho
            );
            assert!(
                (vxc_slice[i] - (vx + vc)).abs() < 1e-15,
                "Slice/pointwise vxc mismatch at index {}, rho={}",
                i,
                rho
            );
        }
    }

    // =========================================================================
    // Slater Scaling Test
    // =========================================================================

    #[test]
    fn slater_scaling_relationship() {
        // Verify eps_x(rho) = -C_x * rho^{1/3}
        let slater = SlaterExchange;
        let test_rhos = [0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 100.0];
        for &rho in &test_rhos {
            let (eps_x, _) = slater.eval(rho);
            let expected = -SlaterExchange::C_X * rho.cbrt();
            assert!(
                (eps_x - expected).abs() < 1e-15,
                "Slater scaling mismatch at rho={}: {} vs {}",
                rho,
                eps_x,
                expected
            );
        }
    }

    // =========================================================================
    // Becke88 Exchange Tests (GGA)
    // =========================================================================

    // PySCF reference: dft.libxc.eval_xc('gga_x_b88', rho_gga)
    // Test pairs: (rho, sigma) = (0.01,0.001), (0.10,0.01), (0.50,0.1),
    //                            (1.00,0.5), (5.00,5.0), (10.00,10.0)
    const TEST_SIGMAS: [f64; 6] = [0.001, 0.01, 0.1, 0.5, 5.0, 10.0];

    // Becke88 exchange: includes Slater exchange
    const PYSCF_EPS_B88: [f64; 6] = [
        -2.506652126879046e-01,
        -3.530065209596222e-01,
        -5.888029649933622e-01,
        -7.411577988590337e-01,
        -1.263534983272655e+00,
        -1.591422034099187e+00,
    ];
    const PYSCF_VRHO_B88: [f64; 6] = [
        -1.878323091027048e-01,
        -4.456959951430178e-01,
        -7.782517401608836e-01,
        -9.813917809027930e-01,
        -1.683071630954492e+00,
        -2.121242183556790e+00,
    ];
    const PYSCF_VSIGMA_B88: [f64; 6] = [
        -5.489549043043804e-01,
        -9.367262301179430e-02,
        -1.278539968174857e-02,
        -5.113963181938868e-03,
        -6.156300283933023e-04,
        -2.451982157970300e-04,
    ];

    #[test]
    fn becke88_golden_values() {
        let b88 = Becke88Exchange;
        for i in 0..TEST_DENSITIES.len() {
            let (eps, vrho, vsigma) = b88.eval_gga(TEST_DENSITIES[i], TEST_SIGMAS[i]);
            assert!(
                (eps - PYSCF_EPS_B88[i]).abs() < 1e-10,
                "B88 eps mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                eps,
                PYSCF_EPS_B88[i],
                (eps - PYSCF_EPS_B88[i]).abs()
            );
            assert!(
                (vrho - PYSCF_VRHO_B88[i]).abs() < 1e-10,
                "B88 vrho mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                vrho,
                PYSCF_VRHO_B88[i],
                (vrho - PYSCF_VRHO_B88[i]).abs()
            );
            assert!(
                (vsigma - PYSCF_VSIGMA_B88[i]).abs() < 1e-10,
                "B88 vsigma mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                vsigma,
                PYSCF_VSIGMA_B88[i],
                (vsigma - PYSCF_VSIGMA_B88[i]).abs()
            );
        }
    }

    #[test]
    fn becke88_reduces_to_slater_at_zero_gradient() {
        // At sigma=0, B88 should give exactly Slater exchange
        let b88 = Becke88Exchange;
        let slater = SlaterExchange;
        for &rho in &TEST_DENSITIES {
            let (eps_b88, vrho_b88, vsigma_b88) = b88.eval_gga(rho, 0.0);
            let (eps_slater, v_slater) = slater.eval(rho);
            assert!(
                (eps_b88 - eps_slater).abs() < 1e-14,
                "B88 should equal Slater at sigma=0: rho={}, B88={:.15e}, Slater={:.15e}",
                rho,
                eps_b88,
                eps_slater
            );
            assert!(
                (vrho_b88 - v_slater).abs() < 1e-14,
                "B88 vrho should equal Slater v_x at sigma=0: rho={}",
                rho
            );
            assert_eq!(vsigma_b88, 0.0, "B88 vsigma should be 0 at sigma=0");
        }
    }

    #[test]
    fn becke88_zero_density() {
        let b88 = Becke88Exchange;
        let (eps, vrho, vsigma) = b88.eval_gga(0.0, 1.0);
        assert_eq!(eps, 0.0);
        assert_eq!(vrho, 0.0);
        assert_eq!(vsigma, 0.0);
    }

    #[test]
    fn becke88_trait_interface() {
        let b88 = Becke88Exchange;
        assert!(b88.needs_gradient());
        assert_eq!(b88.hf_exchange_fraction(), 0.0);
        assert_eq!(b88.name(), "Becke88 Exchange");
    }

    #[test]
    fn becke88_eps_always_negative() {
        let b88 = Becke88Exchange;
        let test_rhos = [0.01, 0.1, 1.0, 10.0, 100.0];
        let test_sigmas = [0.0, 0.001, 0.1, 1.0, 10.0];
        for &rho in &test_rhos {
            for &sigma in &test_sigmas {
                let (eps, _, _) = b88.eval_gga(rho, sigma);
                assert!(
                    eps < 0.0,
                    "B88 eps should be negative at rho={}, sigma={}, got {}",
                    rho,
                    sigma,
                    eps
                );
            }
        }
    }

    // =========================================================================
    // LYP Correlation Tests (GGA)
    // =========================================================================

    // PySCF reference: dft.libxc.eval_xc('gga_c_lyp', rho_gga)
    const PYSCF_EPS_LYP: [f64; 6] = [
        1.463107561391230e-02,
        -3.287738143188335e-02,
        -4.335580886547560e-02,
        -4.701120977382874e-02,
        -5.416948681855666e-02,
        -5.657833371646331e-02,
    ];
    const PYSCF_VRHO_LYP: [f64; 6] = [
        -7.172821668465290e-02,
        -4.233661639110011e-02,
        -4.935191701685578e-02,
        -5.243885400951806e-02,
        -5.795256673206874e-02,
        -5.978598211515272e-02,
    ];
    const PYSCF_VSIGMA_LYP: [f64; 6] = [
        3.559820169927601e-01,
        1.359846061050416e-02,
        1.065244936754465e-03,
        3.415904564702420e-04,
        2.306969617864156e-05,
        7.132433916288656e-06,
    ];

    #[test]
    fn lyp_golden_values() {
        let lyp = LypCorrelation;
        for i in 0..TEST_DENSITIES.len() {
            let (eps, vrho, vsigma) = lyp.eval_gga(TEST_DENSITIES[i], TEST_SIGMAS[i]);
            assert!(
                (eps - PYSCF_EPS_LYP[i]).abs() < 1e-10,
                "LYP eps mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                eps,
                PYSCF_EPS_LYP[i],
                (eps - PYSCF_EPS_LYP[i]).abs()
            );
            assert!(
                (vrho - PYSCF_VRHO_LYP[i]).abs() < 1e-10,
                "LYP vrho mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                vrho,
                PYSCF_VRHO_LYP[i],
                (vrho - PYSCF_VRHO_LYP[i]).abs()
            );
            assert!(
                (vsigma - PYSCF_VSIGMA_LYP[i]).abs() < 1e-10,
                "LYP vsigma mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                vsigma,
                PYSCF_VSIGMA_LYP[i],
                (vsigma - PYSCF_VSIGMA_LYP[i]).abs()
            );
        }
    }

    #[test]
    fn lyp_zero_density() {
        let lyp = LypCorrelation;
        let (eps, vrho, vsigma) = lyp.eval_gga(0.0, 1.0);
        assert_eq!(eps, 0.0);
        assert_eq!(vrho, 0.0);
        assert_eq!(vsigma, 0.0);
    }

    #[test]
    fn lyp_trait_interface() {
        let lyp = LypCorrelation;
        assert!(lyp.needs_gradient());
        assert_eq!(lyp.hf_exchange_fraction(), 0.0);
        assert_eq!(lyp.name(), "LYP Correlation");
    }

    // =========================================================================
    // B3LYP Composite Tests
    // =========================================================================

    // PySCF reference: dft.libxc.eval_xc('b3lyp5', rho_gga)
    // B3LYP5 = VWN5 variant
    const PYSCF_EPS_B3LYP5: [f64; 6] = [
        -1.885097810531995e-01,
        -3.183655479801289e-01,
        -5.184717622777433e-01,
        -6.443999927441574e-01,
        -1.070890132187639e+00,
        -1.336357989114059e+00,
    ];
    const PYSCF_VRHO_B3LYP5: [f64; 6] = [
        -2.186474734806739e-01,
        -4.033143135214302e-01,
        -6.769012560970256e-01,
        -8.430454485485736e-01,
        -1.411379532629600e+00,
        -1.764573522264485e+00,
    ];
    const PYSCF_VSIGMA_B3LYP5: [f64; 6] = [
        -1.069020973350182e-01,
        -5.642953547398352e-02,
        -8.342639372087853e-03,
        -3.405365221255089e-03,
        -4.245671665384780e-04,
        -1.707654439016678e-04,
    ];

    #[test]
    fn b3lyp_golden_values() {
        let b3lyp = B3lyp::new();
        let n = TEST_DENSITIES.len();
        let mut exc = vec![0.0; n];
        let mut vrho = vec![0.0; n];
        let mut vsigma = vec![0.0; n];

        b3lyp.eval_xc_gga(
            &TEST_DENSITIES,
            &TEST_SIGMAS,
            &mut exc,
            &mut vrho,
            &mut vsigma,
        );

        for i in 0..n {
            assert!(
                (exc[i] - PYSCF_EPS_B3LYP5[i]).abs() < 1e-10,
                "B3LYP exc mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                exc[i],
                PYSCF_EPS_B3LYP5[i],
                (exc[i] - PYSCF_EPS_B3LYP5[i]).abs()
            );
            assert!(
                (vrho[i] - PYSCF_VRHO_B3LYP5[i]).abs() < 1e-10,
                "B3LYP vrho mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                vrho[i],
                PYSCF_VRHO_B3LYP5[i],
                (vrho[i] - PYSCF_VRHO_B3LYP5[i]).abs()
            );
            assert!(
                (vsigma[i] - PYSCF_VSIGMA_B3LYP5[i]).abs() < 1e-10,
                "B3LYP vsigma mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                vsigma[i],
                PYSCF_VSIGMA_B3LYP5[i],
                (vsigma[i] - PYSCF_VSIGMA_B3LYP5[i]).abs()
            );
        }
    }

    #[test]
    fn b3lyp_trait_interface() {
        let b3lyp = B3lyp::new();
        assert!(b3lyp.needs_gradient());
        assert_eq!(b3lyp.hf_exchange_fraction(), 0.20);
        assert_eq!(b3lyp.name(), "B3LYP");
    }

    #[test]
    fn b3lyp_uses_vwn5() {
        // Verify B3LYP uses VWN5 (not VWN3/RPA) by comparing against
        // PySCF b3lyp5 reference values. If VWN3 were used, the values
        // would differ by ~1e-3 (the known VWN3/VWN5 discrepancy).
        let b3lyp = B3lyp::new();
        let mut exc = [0.0];
        let mut vrho = [0.0];
        let mut vsigma = [0.0];
        b3lyp.eval_xc_gga(&[1.0], &[0.5], &mut exc, &mut vrho, &mut vsigma);

        // This value is from b3lyp5 (VWN5). If VWN3 were used, exc would be
        // -6.482394766934618e-01 (diff ~3.8e-3 from VWN5 value).
        assert!(
            (exc[0] - PYSCF_EPS_B3LYP5[3]).abs() < 1e-10,
            "B3LYP should use VWN5: got {:.15e}, expected {:.15e} (b3lyp5)",
            exc[0],
            PYSCF_EPS_B3LYP5[3]
        );
    }

    #[test]
    fn b3lyp_mixing_coefficients() {
        // Verify the mixing coefficients are correct
        assert_eq!(B3lyp::A, 0.20);
        assert_eq!(B3lyp::B_PARAM, 0.72);
        assert_eq!(B3lyp::C_PARAM, 0.81);

        // Derived
        assert!((B3lyp::SLATER_COEFF - 0.08).abs() < 1e-15);
        assert!((B3lyp::VWN5_COEFF - 0.19).abs() < 1e-15);
    }

    // =========================================================================
    // Second Derivative (fxc) Tests — US-094
    // PySCF 2.11.0, dft.libxc.eval_xc(..., deriv=2)
    // =========================================================================

    // --- Slater fxc reference: libxc.eval_xc('lda_x', rho, deriv=2) ---
    const PYSCF_V2RHO2_SLATER: [f64; 6] = [
        -7.071896119647036e+00,
        -1.523593832446944e+00,
        -5.210617611978479e-01,
        -3.282483406142321e-01,
        -1.122593533973753e-01,
        -7.071896119647036e-02,
    ];

    // --- VWN5 fxc reference: libxc.eval_xc('lda_c_vwn', rho, deriv=2) ---
    const PYSCF_V2RHO2_VWN5: [f64; 6] = [
        -6.792284745355519e-01,
        -7.876341792538874e-02,
        -1.694895124061815e-02,
        -8.693785423436672e-03,
        -1.824645219162099e-03,
        -9.277261503265210e-04,
    ];

    // --- Becke88 fxc reference: libxc.eval_xc('gga_x_b88', rho_gga, deriv=2) ---
    const PYSCF_V2RHO2_B88: [f64; 6] = [
        -4.722623270100552e+00,
        -1.728078709928415e+00,
        -5.357067809128221e-01,
        -3.356776447733500e-01,
        -1.126380150097347e-01,
        -7.079495907164067e-02,
    ];
    const PYSCF_V2RHOSIGMA_B88: [f64; 6] = [
        -5.769201374961042e+00,
        9.090952229438343e-01,
        3.163553901043711e-02,
        6.410288354314328e-03,
        1.624648547881943e-04,
        3.258235740535974e-05,
    ];
    const PYSCF_V2SIGMA2_B88: [f64; 6] = [
        2.961119573082944e+02,
        1.274524064550337e+00,
        4.610362764173303e-03,
        3.062469162031238e-04,
        6.386822937574376e-07,
        4.152676284160585e-08,
    ];

    // --- LYP fxc reference: libxc.eval_xc('gga_c_lyp', rho_gga, deriv=2) ---
    const PYSCF_V2RHO2_LYP: [f64; 6] = [
        9.282375637359824e+00,
        -7.832598282198205e-03,
        -7.874099681665563e-03,
        -3.506854587775689e-03,
        -5.658288706070812e-04,
        -2.411261979204161e-04,
    ];
    const PYSCF_V2RHOSIGMA_LYP: [f64; 6] = [
        -4.557927359955423e+01,
        -2.079875296827534e-01,
        -3.465539825496245e-03,
        -5.648978842763125e-04,
        -7.796833398042325e-06,
        -1.210022762293925e-06,
    ];
    // LYP v2sigma2 = 0 for all densities (closed-shell LYP is linear in sigma)

    // --- LDA composite fxc reference ---
    const PYSCF_V2RHO2_LDA: [f64; 6] = [
        -7.751124594182588e+00,
        -1.602357250372333e+00,
        -5.380107124384660e-01,
        -3.369421260376688e-01,
        -1.140839986165374e-01,
        -7.164668734679688e-02,
    ];

    // --- B3LYP5 composite fxc reference ---
    const PYSCF_V2RHO2_B3LYP5: [f64; 6] = [
        3.423630412055543e+00,
        -1.387413631758618e+00,
        -4.369921446309263e-01,
        -2.724401429325019e-01,
        -9.088512305563154e-02,
        -5.700146761617648e-02,
    ];
    const PYSCF_V2RHOSIGMA_B3LYP5: [f64; 6] = [
        -4.107303660561088e+01,
        4.860786614765304e-01,
        1.997050082886277e-02,
        4.157840328842503e-03,
        1.106592603950856e-04,
        2.247917889440093e-05,
    ];
    const PYSCF_V2SIGMA2_B3LYP5: [f64; 6] = [
        2.132006092619719e+02,
        9.176573264762427e-01,
        3.319461190204778e-03,
        2.204977796662492e-04,
        4.598512515053550e-07,
        2.989926924595621e-08,
    ];

    // --- Test 1: Slater fxc (AC2) ---
    #[test]
    fn test_slater_fxc() {
        let slater = SlaterExchange;
        let n = TEST_DENSITIES.len();
        let mut v2rho2 = vec![0.0; n];
        let mut v2rhosigma = vec![0.0; n];
        let mut v2sigma2 = vec![0.0; n];
        let sigma = vec![0.0; n];

        slater.eval_xc_second_deriv(
            &TEST_DENSITIES,
            &sigma,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );

        for i in 0..n {
            assert!(
                (v2rho2[i] - PYSCF_V2RHO2_SLATER[i]).abs() < 1e-10,
                "Slater v2rho2 mismatch at rho={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                v2rho2[i],
                PYSCF_V2RHO2_SLATER[i],
                (v2rho2[i] - PYSCF_V2RHO2_SLATER[i]).abs()
            );
            assert_eq!(v2rhosigma[i], 0.0, "Slater v2rhosigma should be 0");
            assert_eq!(v2sigma2[i], 0.0, "Slater v2sigma2 should be 0");
        }
    }

    // --- Test 2: VWN5 fxc (AC3) ---
    #[test]
    fn test_vwn5_fxc() {
        let vwn5 = Vwn5Correlation;
        let n = TEST_DENSITIES.len();
        let mut v2rho2 = vec![0.0; n];
        let mut v2rhosigma = vec![0.0; n];
        let mut v2sigma2 = vec![0.0; n];
        let sigma = vec![0.0; n];

        vwn5.eval_xc_second_deriv(
            &TEST_DENSITIES,
            &sigma,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );

        for i in 0..n {
            assert!(
                (v2rho2[i] - PYSCF_V2RHO2_VWN5[i]).abs() < 1e-10,
                "VWN5 v2rho2 mismatch at rho={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                v2rho2[i],
                PYSCF_V2RHO2_VWN5[i],
                (v2rho2[i] - PYSCF_V2RHO2_VWN5[i]).abs()
            );
            assert_eq!(v2rhosigma[i], 0.0, "VWN5 v2rhosigma should be 0");
            assert_eq!(v2sigma2[i], 0.0, "VWN5 v2sigma2 should be 0");
        }
    }

    // --- Test 3: Becke88 fxc (AC4) ---
    #[test]
    fn test_becke88_fxc() {
        let b88 = Becke88Exchange;
        let n = TEST_DENSITIES.len();
        let mut v2rho2 = vec![0.0; n];
        let mut v2rhosigma = vec![0.0; n];
        let mut v2sigma2 = vec![0.0; n];

        b88.eval_xc_second_deriv(
            &TEST_DENSITIES,
            &TEST_SIGMAS,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );

        for i in 0..n {
            assert!(
                (v2rho2[i] - PYSCF_V2RHO2_B88[i]).abs() < 1e-8,
                "B88 v2rho2 mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                v2rho2[i],
                PYSCF_V2RHO2_B88[i],
                (v2rho2[i] - PYSCF_V2RHO2_B88[i]).abs()
            );
            assert!(
                (v2rhosigma[i] - PYSCF_V2RHOSIGMA_B88[i]).abs() < 1e-8,
                "B88 v2rhosigma mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                v2rhosigma[i],
                PYSCF_V2RHOSIGMA_B88[i],
                (v2rhosigma[i] - PYSCF_V2RHOSIGMA_B88[i]).abs()
            );
            assert!(
                (v2sigma2[i] - PYSCF_V2SIGMA2_B88[i]).abs() < 1e-8,
                "B88 v2sigma2 mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                v2sigma2[i],
                PYSCF_V2SIGMA2_B88[i],
                (v2sigma2[i] - PYSCF_V2SIGMA2_B88[i]).abs()
            );
        }
    }

    // --- Test 4: LYP fxc (AC5) ---
    #[test]
    fn test_lyp_fxc() {
        let lyp = LypCorrelation;
        let n = TEST_DENSITIES.len();
        let mut v2rho2 = vec![0.0; n];
        let mut v2rhosigma = vec![0.0; n];
        let mut v2sigma2 = vec![0.0; n];

        lyp.eval_xc_second_deriv(
            &TEST_DENSITIES,
            &TEST_SIGMAS,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );

        for i in 0..n {
            assert!(
                (v2rho2[i] - PYSCF_V2RHO2_LYP[i]).abs() < 1e-8,
                "LYP v2rho2 mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                v2rho2[i],
                PYSCF_V2RHO2_LYP[i],
                (v2rho2[i] - PYSCF_V2RHO2_LYP[i]).abs()
            );
            assert!(
                (v2rhosigma[i] - PYSCF_V2RHOSIGMA_LYP[i]).abs() < 1e-8,
                "LYP v2rhosigma mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                v2rhosigma[i],
                PYSCF_V2RHOSIGMA_LYP[i],
                (v2rhosigma[i] - PYSCF_V2RHOSIGMA_LYP[i]).abs()
            );
            assert_eq!(
                v2sigma2[i], 0.0,
                "LYP v2sigma2 should be 0 (linear in sigma for closed-shell)"
            );
        }
    }

    // --- Test 5: LDA composite fxc ---
    #[test]
    fn test_lda_fxc() {
        let lda = Lda::new();
        let n = TEST_DENSITIES.len();
        let mut v2rho2 = vec![0.0; n];
        let mut v2rhosigma = vec![0.0; n];
        let mut v2sigma2 = vec![0.0; n];
        let sigma = vec![0.0; n];

        lda.eval_xc_second_deriv(
            &TEST_DENSITIES,
            &sigma,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );

        for i in 0..n {
            assert!(
                (v2rho2[i] - PYSCF_V2RHO2_LDA[i]).abs() < 1e-10,
                "LDA v2rho2 mismatch at rho={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                v2rho2[i],
                PYSCF_V2RHO2_LDA[i],
                (v2rho2[i] - PYSCF_V2RHO2_LDA[i]).abs()
            );
        }

        // Verify additivity: LDA fxc = Slater fxc + VWN5 fxc
        for i in 0..n {
            let expected = PYSCF_V2RHO2_SLATER[i] + PYSCF_V2RHO2_VWN5[i];
            assert!(
                (v2rho2[i] - expected).abs() < 1e-10,
                "LDA additivity failed at rho={}: {} vs {}",
                TEST_DENSITIES[i],
                v2rho2[i],
                expected
            );
        }
    }

    // --- Test 6: B3LYP composite fxc ---
    #[test]
    fn test_b3lyp_fxc() {
        let b3lyp = B3lyp::new();
        let n = TEST_DENSITIES.len();
        let mut v2rho2 = vec![0.0; n];
        let mut v2rhosigma = vec![0.0; n];
        let mut v2sigma2 = vec![0.0; n];

        b3lyp.eval_xc_second_deriv(
            &TEST_DENSITIES,
            &TEST_SIGMAS,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );

        for i in 0..n {
            assert!(
                (v2rho2[i] - PYSCF_V2RHO2_B3LYP5[i]).abs() < 1e-8,
                "B3LYP v2rho2 mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                v2rho2[i],
                PYSCF_V2RHO2_B3LYP5[i],
                (v2rho2[i] - PYSCF_V2RHO2_B3LYP5[i]).abs()
            );
            assert!(
                (v2rhosigma[i] - PYSCF_V2RHOSIGMA_B3LYP5[i]).abs() < 1e-8,
                "B3LYP v2rhosigma mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                v2rhosigma[i],
                PYSCF_V2RHOSIGMA_B3LYP5[i],
                (v2rhosigma[i] - PYSCF_V2RHOSIGMA_B3LYP5[i]).abs()
            );
            assert!(
                (v2sigma2[i] - PYSCF_V2SIGMA2_B3LYP5[i]).abs() < 1e-8,
                "B3LYP v2sigma2 mismatch at rho={}, sigma={}: {:.15e} vs {:.15e} (diff={:.2e})",
                TEST_DENSITIES[i],
                TEST_SIGMAS[i],
                v2sigma2[i],
                PYSCF_V2SIGMA2_B3LYP5[i],
                (v2sigma2[i] - PYSCF_V2SIGMA2_B3LYP5[i]).abs()
            );
        }
    }

    // --- Test 7: Low density guard ---
    #[test]
    fn test_fxc_low_density_guard() {
        // Test at DENSITY_THRESHOLD boundary: rho <= 1e-20 should return 0.
        let rho_at_threshold = [DENSITY_THRESHOLD];
        let rho_below = [1e-25];
        let sigma_tiny = [1e-30];
        let mut v2rho2 = [0.0];
        let mut v2rhosigma = [0.0];
        let mut v2sigma2 = [0.0];

        // All functionals should return 0 at the threshold
        for rho_test in [&rho_at_threshold[..], &rho_below[..]] {
            // Slater
            SlaterExchange.eval_xc_second_deriv(
                rho_test,
                &sigma_tiny,
                &mut v2rho2,
                &mut v2rhosigma,
                &mut v2sigma2,
            );
            assert_eq!(v2rho2[0], 0.0, "Slater fxc should be 0 at/below threshold");

            // VWN5
            Vwn5Correlation.eval_xc_second_deriv(
                rho_test,
                &sigma_tiny,
                &mut v2rho2,
                &mut v2rhosigma,
                &mut v2sigma2,
            );
            assert_eq!(v2rho2[0], 0.0, "VWN5 fxc should be 0 at/below threshold");

            // Becke88
            Becke88Exchange.eval_xc_second_deriv(
                rho_test,
                &sigma_tiny,
                &mut v2rho2,
                &mut v2rhosigma,
                &mut v2sigma2,
            );
            assert_eq!(v2rho2[0], 0.0, "B88 fxc should be 0 at/below threshold");

            // LYP
            LypCorrelation.eval_xc_second_deriv(
                rho_test,
                &sigma_tiny,
                &mut v2rho2,
                &mut v2rhosigma,
                &mut v2sigma2,
            );
            assert_eq!(v2rho2[0], 0.0, "LYP fxc should be 0 at/below threshold");

            // LDA
            Lda::new().eval_xc_second_deriv(
                rho_test,
                &sigma_tiny,
                &mut v2rho2,
                &mut v2rhosigma,
                &mut v2sigma2,
            );
            assert_eq!(v2rho2[0], 0.0, "LDA fxc should be 0 at/below threshold");

            // B3LYP
            B3lyp::new().eval_xc_second_deriv(
                rho_test,
                &sigma_tiny,
                &mut v2rho2,
                &mut v2rhosigma,
                &mut v2sigma2,
            );
            assert_eq!(v2rho2[0], 0.0, "B3LYP fxc should be 0 at/below threshold");
        }

        // Above threshold: finite, no NaN/Inf
        let rho_small = [1e-10];
        let sigma_small = [1e-20];

        SlaterExchange.eval_xc_second_deriv(
            &rho_small,
            &sigma_small,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );
        assert!(
            v2rho2[0].is_finite(),
            "Slater fxc should be finite above threshold"
        );
        assert!(
            v2rho2[0] != 0.0,
            "Slater fxc should be nonzero above threshold"
        );

        Vwn5Correlation.eval_xc_second_deriv(
            &rho_small,
            &sigma_small,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );
        assert!(
            v2rho2[0].is_finite(),
            "VWN5 fxc should be finite above threshold"
        );
        assert!(
            v2rho2[0] != 0.0,
            "VWN5 fxc should be nonzero above threshold"
        );

        Becke88Exchange.eval_xc_second_deriv(
            &rho_small,
            &sigma_small,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );
        assert!(
            v2rho2[0].is_finite(),
            "B88 fxc should be finite above threshold"
        );

        LypCorrelation.eval_xc_second_deriv(
            &rho_small,
            &sigma_small,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );
        assert!(
            v2rho2[0].is_finite(),
            "LYP fxc should be finite above threshold"
        );

        B3lyp::new().eval_xc_second_deriv(
            &rho_small,
            &sigma_small,
            &mut v2rho2,
            &mut v2rhosigma,
            &mut v2sigma2,
        );
        assert!(
            v2rho2[0].is_finite(),
            "B3LYP fxc should be finite above threshold"
        );
        assert!(
            v2rhosigma[0].is_finite(),
            "B3LYP fxc v2rhosigma should be finite"
        );
        assert!(
            v2sigma2[0].is_finite(),
            "B3LYP fxc v2sigma2 should be finite"
        );
    }

    // --- Test 8: Finite difference validation ---
    #[test]
    fn test_fxc_finite_diff() {
        // Compare analytical fxc against FD of first derivative:
        // v2rho2 ≈ [vrho(rho+h) - vrho(rho-h)] / (2h)
        //
        // Use relative tolerance: |analytical - fd| / max(|analytical|, 1) < tol
        // Central FD with h=1e-5 has O(h²) truncation error ~ 1e-10 for smooth
        // functions, but the ratio can be larger for stiff functions at low density.
        let h = 1e-6;
        let fd_rtol = 1e-5; // Relative tolerance for FD comparison
        let test_rhos = [0.01, 0.1, 1.0, 5.0, 10.0];

        // Helper: relative error check
        let check_rel = |label: &str, analytical: f64, fd: f64, rho: f64, sigma: f64| {
            let scale = analytical.abs().max(1.0);
            let rel_err = (analytical - fd).abs() / scale;
            assert!(
                rel_err < fd_rtol,
                "{} FD mismatch at rho={}, sigma={}: analytical={:.10e}, fd={:.10e}, rel_err={:.2e}",
                label, rho, sigma, analytical, fd, rel_err
            );
        };

        // --- Slater ---
        let slater = SlaterExchange;
        for &rho in &test_rhos {
            let fxc = slater.eval_fxc(rho);
            let (_, v_plus) = slater.eval(rho + h);
            let (_, v_minus) = slater.eval(rho - h);
            let fd = (v_plus - v_minus) / (2.0 * h);
            check_rel("Slater v2rho2", fxc, fd, rho, 0.0);
        }

        // --- VWN5 ---
        let vwn5 = Vwn5Correlation;
        for &rho in &test_rhos {
            let fxc = vwn5.eval_fxc(rho);
            let (_, v_plus) = vwn5.eval(rho + h);
            let (_, v_minus) = vwn5.eval(rho - h);
            let fd = (v_plus - v_minus) / (2.0 * h);
            check_rel("VWN5 v2rho2", fxc, fd, rho, 0.0);
        }

        // --- Becke88 ---
        let b88 = Becke88Exchange;
        let test_sigmas_fd = [0.01, 0.1, 0.5, 5.0, 10.0];
        for (&rho, &sigma) in test_rhos.iter().zip(test_sigmas_fd.iter()) {
            let (v2rho2, v2rhosigma, v2sigma2) = b88.eval_fxc_gga(rho, sigma);

            // FD for v2rho2: d(vrho)/drho at fixed sigma
            let (_, vrho_plus, _) = b88.eval_gga(rho + h, sigma);
            let (_, vrho_minus, _) = b88.eval_gga(rho - h, sigma);
            let fd_v2rho2 = (vrho_plus - vrho_minus) / (2.0 * h);
            check_rel("B88 v2rho2", v2rho2, fd_v2rho2, rho, sigma);

            // FD for v2rhosigma: d(vrho)/dsigma at fixed rho
            let h_sig = h;
            let (_, vrho_sigp, _) = b88.eval_gga(rho, sigma + h_sig);
            let (_, vrho_sigm, _) = b88.eval_gga(rho, sigma - h_sig);
            let fd_v2rhosigma = (vrho_sigp - vrho_sigm) / (2.0 * h_sig);
            check_rel("B88 v2rhosigma", v2rhosigma, fd_v2rhosigma, rho, sigma);

            // FD for v2sigma2: d(vsigma)/dsigma at fixed rho
            let (_, _, vsig_plus) = b88.eval_gga(rho, sigma + h_sig);
            let (_, _, vsig_minus) = b88.eval_gga(rho, sigma - h_sig);
            let fd_v2sigma2 = (vsig_plus - vsig_minus) / (2.0 * h_sig);
            check_rel("B88 v2sigma2", v2sigma2, fd_v2sigma2, rho, sigma);
        }

        // --- LYP ---
        let lyp = LypCorrelation;
        for (&rho, &sigma) in test_rhos.iter().zip(test_sigmas_fd.iter()) {
            let (v2rho2, v2rhosigma, _v2sigma2) = lyp.eval_fxc_gga(rho, sigma);

            // FD for v2rho2
            let (_, vrho_plus, _) = lyp.eval_gga(rho + h, sigma);
            let (_, vrho_minus, _) = lyp.eval_gga(rho - h, sigma);
            let fd_v2rho2 = (vrho_plus - vrho_minus) / (2.0 * h);
            check_rel("LYP v2rho2", v2rho2, fd_v2rho2, rho, sigma);

            // FD for v2rhosigma
            let h_sig = h;
            let (_, vrho_sigp, _) = lyp.eval_gga(rho, sigma + h_sig);
            let (_, vrho_sigm, _) = lyp.eval_gga(rho, sigma - h_sig);
            let fd_v2rhosigma = (vrho_sigp - vrho_sigm) / (2.0 * h_sig);
            check_rel("LYP v2rhosigma", v2rhosigma, fd_v2rhosigma, rho, sigma);
        }
    }
}
