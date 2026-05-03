//! Property-based tests for qc-core algorithms using proptest.
//!
//! These tests verify mathematical invariants that must hold across the full
//! domain of each algorithm, going beyond the fixed test points used in unit tests.
//!
//! # Properties tested
//!
//! 1. Boys function non-negativity: F_m(T) >= 0 for all valid m, T
//! 2. Boys function monotonic decrease in T: F_m(T1) >= F_m(T2) when T1 < T2
//! 3. Boys function recurrence consistency: upward recurrence matches direct evaluation
//! 4. Rys quadrature roots bounded in [0, 1)
//! 5. Rys quadrature weights non-negative
//! 6. Rys moment reconstruction: sum_i w_i * r_i^(2k) ~ F_k(T)
//! 7. Overlap matrix positive semi-definite (all eigenvalues >= 0)
//! 8. SCF energy below nuclear repulsion (electrons lower the energy)
//!
//! # References
//!
//! - Shavitt, I. (1963). Methods in Computational Physics, Vol. 2.
//! - Dupuis, M., Rys, J., & King, H. F. (1976). J. Chem. Phys. 65, 111.

use proptest::prelude::*;

// =============================================================================
// Boys Function Properties
//
// The Boys function F_m(T) = integral_0^1 t^{2m} exp(-T t^2) dt satisfies:
//   (i)   F_m(T) >= 0 for all m >= 0, T >= 0
//   (ii)  F_m is monotonically decreasing in T (since exp(-T t^2) decreases)
//   (iii) The upward recurrence F_{m+1}(T) = [(2m+1) F_m(T) - exp(-T)] / (2T)
//         connects adjacent orders
// =============================================================================

proptest! {
    /// **Property 1: Non-negativity**
    ///
    /// F_m(T) is defined as an integral of a non-negative integrand (t^{2m} >= 0
    /// and exp(-T t^2) >= 0 on [0,1]), so the result must be non-negative.
    ///
    /// Reference: Shavitt (1963), Eq. 1 -- integral definition.
    #[test]
    fn boys_nonnegative(m in 0u32..=20, t in 0.0f64..120.0) {
        let result = qc_core::boys::boys_eval(m, t).unwrap();
        prop_assert!(
            result.value >= 0.0,
            "F_{}({}) = {} is negative",
            m, t, result.value,
        );
    }

    /// **Property 2: Monotonic decrease in T**
    ///
    /// For fixed m, increasing T makes exp(-T t^2) smaller at every interior point
    /// of [0,1], so the integral F_m(T) must decrease. Formally:
    ///   T1 < T2  =>  F_m(T1) >= F_m(T2)
    ///
    /// We pick two random T values and compare. The tolerance accounts for
    /// floating-point evaluation error (1e-12 absolute, per spec).
    #[test]
    fn boys_monotone_decreasing_in_t(m in 0u32..=15, t1 in 0.0f64..60.0, delta in 0.01f64..40.0) {
        let t2 = t1 + delta;
        let v1 = qc_core::boys::boys_eval(m, t1).unwrap().value;
        let v2 = qc_core::boys::boys_eval(m, t2).unwrap().value;
        prop_assert!(
            v1 >= v2 - 1e-12,
            "F_{}({}) = {} < F_{}({}) = {} -- monotonicity violated",
            m, t1, v1, m, t2, v2,
        );
    }

    /// **Property 3: Upward recurrence consistency**
    ///
    /// The Boys function satisfies the recurrence relation (Shavitt 1963, Eq. 7-8):
    ///
    ///   F_{m+1}(T) = [(2m+1) F_m(T) - exp(-T)] / (2T)
    ///
    /// We verify that the directly computed F_{m+1}(T) agrees with the value
    /// obtained by applying this recurrence to the directly computed F_m(T).
    ///
    /// Skip T < 0.01 to avoid division by near-zero (2T) in the recurrence.
    #[test]
    fn boys_recurrence_consistency(m in 0u32..=14, t in 0.01f64..80.0) {
        let fm = qc_core::boys::boys_eval(m, t).unwrap().value;
        let fm1_direct = qc_core::boys::boys_eval(m + 1, t).unwrap().value;

        // Upward recurrence: F_{m+1} = [(2m+1) F_m - exp(-T)] / (2T)
        let exp_neg_t = (-t).exp();
        let fm1_recurrence = ((2 * m + 1) as f64 * fm - exp_neg_t) / (2.0 * t);

        let abs_err = (fm1_direct - fm1_recurrence).abs();

        // Tolerance: allow accumulated rounding at 1e-10 level.
        // Both sides have ~1e-12 absolute error; the subtraction and division
        // can amplify this slightly, especially when F_m is small.
        prop_assert!(
            abs_err < 1e-10,
            "Recurrence mismatch for m={}, T={}: direct={:.15e}, recurrence={:.15e}, err={:.3e}",
            m, t, fm1_direct, fm1_recurrence, abs_err,
        );
    }
}

// =============================================================================
// Rys Quadrature Properties
//
// Rys quadrature computes n-point {root_i, weight_i} sets such that:
//   integral_0^1 P(t^2) exp(-T t^2) dt = sum_i w_i P(root_i)
//   for any polynomial P of degree < 2n.
//
// Required invariants:
//   (i)   All roots in [0, 1)
//   (ii)  All weights >= 0 (strictly positive for non-degenerate T > 0)
//   (iii) Moment reconstruction: sum_i w_i root_i^k ~ F_k(T)
// =============================================================================

proptest! {
    /// **Property 4: Rys roots bounded in [0, 1)**
    ///
    /// All quadrature roots must lie in [0, 1). The upper bound is strict because
    /// the integration variable t is in [0, 1] and roots are t^2.
    /// Root = 0 is allowed for degenerate cases (T = 0).
    ///
    /// Reference: Dupuis, Rys & King (1976); libcint rys_roots.c validate().
    #[test]
    fn rys_roots_in_unit_interval(n in 1usize..=5, t in 0.01f64..50.0) {
        let result = qc_core::rys::rys_roots(n, t).unwrap();
        for (i, &root) in result.roots.iter().enumerate() {
            prop_assert!(
                (0.0..1.0).contains(&root),
                "Root {} = {} out of [0, 1) for n={}, T={}",
                i, root, n, t,
            );
        }
    }

    /// **Property 5: Rys weights non-negative**
    ///
    /// Quadrature weights must be non-negative. For non-degenerate cases (T > 0),
    /// they should in fact be strictly positive.
    ///
    /// Reference: Dupuis, Rys & King (1976).
    #[test]
    fn rys_weights_nonnegative(n in 1usize..=5, t in 0.01f64..50.0) {
        let result = qc_core::rys::rys_roots(n, t).unwrap();
        for (i, &weight) in result.weights.iter().enumerate() {
            prop_assert!(
                weight >= 0.0,
                "Weight {} = {} is negative for n={}, T={}",
                i, weight, n, t,
            );
        }
        // For T > 0 (non-degenerate), weights should be strictly positive
        for (i, &weight) in result.weights.iter().enumerate() {
            prop_assert!(
                weight > 0.0,
                "Weight {} = {} is zero (expected positive) for non-degenerate n={}, T={}",
                i, weight, n, t,
            );
        }
    }

    /// **Property 6: Rys moment reconstruction**
    ///
    /// An n-point Rys quadrature must exactly reproduce the first 2n moments.
    /// The k-th moment is F_k(T), and the quadrature gives:
    ///
    ///   sum_i w_i * root_i^k  ~  F_k(T)   for k = 0, 1, ..., 2n-1
    ///
    /// We test the zeroth moment (k=0) as the most robust check:
    ///   sum_i w_i = F_0(T)
    ///
    /// Tolerance: 1e-10 (Rys moment reconstruction spec from TDD).
    #[test]
    fn rys_zeroth_moment_reconstruction(n in 1usize..=5, t in 0.01f64..50.0) {
        let result = qc_core::rys::rys_roots(n, t).unwrap();
        let f0 = qc_core::boys::boys_eval(0, t).unwrap().value;

        let weight_sum: f64 = result.weights.iter().sum();
        let abs_err = (weight_sum - f0).abs();

        prop_assert!(
            abs_err < 1e-10,
            "Moment reconstruction failed for n={}, T={}: sum(w)={:.15e}, F_0={:.15e}, err={:.3e}",
            n, t, weight_sum, f0, abs_err,
        );
    }
}

// =============================================================================
// Overlap Matrix Properties
//
// The overlap matrix S_uv = <phi_u | phi_v> for a real basis must be:
//   (i)  Symmetric: S = S^T
//   (ii) Positive semi-definite: all eigenvalues >= 0
//   (iii) Diagonal elements = 1 for normalized basis functions
//
// We test with H2 at random bond lengths in STO-3G.
// =============================================================================

proptest! {
    /// **Property 7: Overlap matrix positive semi-definite**
    ///
    /// The overlap matrix S of any real Gaussian basis must have all eigenvalues >= 0,
    /// since S_uv = <phi_u | phi_v> forms a Gram matrix. The diagonal must be 1.0
    /// for normalized basis functions.
    ///
    /// We construct H2 in STO-3G at random bond lengths and verify.
    #[test]
    fn overlap_positive_semidefinite(r in 0.5f64..8.0) {
        let h1 = qc_core::basis::Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = qc_core::basis::Atom::new(1, [0.0, 0.0, r]).unwrap();
        let basis = qc_core::basis::BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
        let s = qc_core::integrals::overlap_matrix(&basis);

        let n = basis.n_basis;

        // Check diagonal elements are ~1.0 (normalized basis).
        // Tolerance 1e-8: contraction normalization introduces small errors
        // at extreme geometries where Gaussian product overlap is very large.
        for i in 0..n {
            prop_assert!(
                (s[i * n + i] - 1.0).abs() < 1e-8,
                "S[{},{}] = {} != 1.0 at R={} bohr",
                i, i, s[i * n + i], r,
            );
        }

        // Check symmetry
        for i in 0..n {
            for j in (i + 1)..n {
                prop_assert!(
                    (s[i * n + j] - s[j * n + i]).abs() < 1e-14,
                    "S not symmetric: S[{},{}]={}, S[{},{}]={} at R={}",
                    i, j, s[i * n + j], j, i, s[j * n + i], r,
                );
            }
        }

        // Positive semi-definite check via eigenvalues.
        // For a 2x2 symmetric matrix [[a, b], [b, a]], eigenvalues are a+b and a-b.
        // With normalized basis: eigenvalues = 1 + S_01 and 1 - S_01.
        // Both must be >= 0 (they will be since |S_01| < 1 for distinct centers).
        let s01 = s[1]; // S[0, 1] at index 0*n + 1 = 1
        let eig1 = 1.0 + s01;
        let eig2 = 1.0 - s01;
        prop_assert!(
            eig1 >= -1e-14,
            "Negative eigenvalue {} at R={} bohr (S_01={})",
            eig1, r, s01,
        );
        prop_assert!(
            eig2 >= -1e-14,
            "Negative eigenvalue {} at R={} bohr (S_01={})",
            eig2, r, s01,
        );

        // Additional: off-diagonal element should be in (-1, 1) for distinct atoms
        prop_assert!(
            s01.abs() < 1.0,
            "Off-diagonal overlap |S_01| = {} >= 1.0 at R={} bohr",
            s01.abs(), r,
        );
    }
}

// =============================================================================
// SCF Energy Properties
//
// For a converged SCF calculation, the total energy must satisfy:
//   E_total = E_electronic + E_nuclear
//   E_total < E_nuclear (electrons lower the energy)
//   E_electronic < 0 (bound state)
//
// We test with H2 at various bond lengths, computing all integrals on-the-fly.
// =============================================================================

proptest! {
    #![proptest_config(ProptestConfig::with_cases(10))]  // Fewer cases: SCF is expensive

    /// **Property 8: SCF energy below nuclear repulsion**
    ///
    /// In a converged SCF calculation, the electrons must lower the total energy
    /// below the bare nuclear repulsion. Equivalently, E_electronic < 0 (the
    /// electronic Hamiltonian expectation value is negative for a bound state).
    ///
    /// We build H2 (STO-3G) from scratch at random bond lengths and run SCF.
    /// The range 0.8..6.0 bohr avoids extreme compression where convergence
    /// may be difficult, and extreme extension where the basis is insufficient.
    #[test]
    fn scf_energy_below_nuclear_repulsion(r in 0.8f64..6.0) {
        let h1 = qc_core::basis::Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2_atom = qc_core::basis::Atom::new(1, [0.0, 0.0, r]).unwrap();
        let basis = qc_core::basis::BasisSet::build(vec![h1, h2_atom], "sto-3g").unwrap();

        let s = qc_core::integrals::overlap_matrix(&basis);
        let hc = qc_core::integrals::hcore_matrix(&basis);
        let eri = qc_core::integrals::eri_compressed(&basis);

        let system = qc_core::scf::PresetSystem {
            system_id: format!("h2_proptest_r{:.4}", r),
            label: format!("H2 (STO-3G, R={:.4} bohr)", r),
            nbf: basis.n_basis,
            nelec: basis.n_electrons,
            e_nuc: basis.nuclear_repulsion,
            s_matrix: s,
            h_core: hc,
            eri_compressed: eri,
        };

        let config = qc_core::scf::ScfConfig {
            profile: qc_core::scf::ConvergenceProfile::Medium,
            max_iterations: 100,
            use_diis: true,
            diis_size: 6,
            diis_start: 2,
            damp: 0.0,
            damp_start: 5,
            level_shift: 0.0,
        };

        let result = qc_core::scf::rhf_scf(&system, &config)
            .expect("SCF should not error for H2 at reasonable bond length");

        // SCF must converge at these geometries
        prop_assert!(
            result.converged,
            "SCF did not converge for H2 at R={:.4} bohr",
            r,
        );

        // Total energy must be below nuclear repulsion
        // (electrons always lower the total energy for a bound system)
        prop_assert!(
            result.energy_total < system.e_nuc,
            "E_total ({:.10}) >= E_nuc ({:.10}) at R={:.4} bohr -- electrons did not lower energy",
            result.energy_total, system.e_nuc, r,
        );

        // Electronic energy must be negative
        prop_assert!(
            result.energy_electronic < 0.0,
            "E_electronic ({:.10}) >= 0 at R={:.4} bohr -- not a bound state",
            result.energy_electronic, r,
        );

        // Energy partitioning consistency: E_total = E_electronic + E_nuclear
        let partition_err = (result.energy_total - (result.energy_electronic + result.energy_nuclear)).abs();
        prop_assert!(
            partition_err < 1e-12,
            "Energy partition inconsistent: E_tot={:.12}, E_el+E_nuc={:.12}, err={:.3e}",
            result.energy_total,
            result.energy_electronic + result.energy_nuclear,
            partition_err,
        );
    }
}
