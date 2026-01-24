// Debug example for p-orbital ERI computation
// Traces through each step to identify the source of discrepancy
#![allow(dead_code, unused_variables, clippy::too_many_arguments)]

use qc_core::integrals::{GaussianProduct2e, RysCoefficients};
use qc_core::rys::rys_roots;
use std::f64::consts::PI;

// Define CartesianPower locally since it's not public
#[derive(Debug, Clone, Copy)]
struct CartesianPower {
    i: u32,
    j: u32,
    k: u32,
}

impl CartesianPower {
    fn new(i: u32, j: u32, k: u32) -> Self {
        Self { i, j, k }
    }
}

// 2D VRR building (copied from module for inspection)
fn build_2d(
    n_bra: usize,
    n_ket: usize,
    c00: f64,
    c0p: f64,
    b00: f64,
    b10: f64,
    b01: f64,
) -> Vec<f64> {
    let rows = n_bra + 1;
    let cols = n_ket + 1;
    let mut g = vec![0.0; rows * cols];
    g[0] = 1.0;

    if n_bra == 0 && n_ket == 0 {
        return g;
    }

    if n_bra > 0 {
        g[cols] = c00;
        for n in 1..n_bra {
            g[(n + 1) * cols] = c00 * g[n * cols] + (n as f64) * b10 * g[(n - 1) * cols];
        }
    }

    for m in 0..n_ket {
        if m == 0 {
            g[m + 1] = c0p * g[0];
        } else {
            g[m + 1] = c0p * g[m] + (m as f64) * b01 * g[m - 1];
        }

        for n in 1..=n_bra {
            let g_n_m = g[n * cols + m];
            let g_n_m_minus_1 = if m > 0 { g[n * cols + m - 1] } else { 0.0 };
            let g_n_minus_1_m = g[(n - 1) * cols + m];
            g[n * cols + m + 1] =
                c0p * g_n_m + (m as f64) * b01 * g_n_m_minus_1 + (n as f64) * b00 * g_n_minus_1_m;
        }
    }

    g
}

fn get_2d(g: &[f64], n_ket: usize, n: usize, m: usize) -> f64 {
    g[n * (n_ket + 1) + m]
}

// HTR for 1D
fn htr_recursive(
    g: &[f64],
    n_ket: usize,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    ab: f64,
    cd: f64,
) -> f64 {
    if j == 0 && l == 0 {
        return get_2d(g, n_ket, i, k);
    }
    if l > 0 {
        let term1 = htr_recursive(g, n_ket, i, j, k + 1, l - 1, ab, cd);
        let term2 = htr_recursive(g, n_ket, i, j, k, l - 1, ab, cd);
        return term1 + cd * term2;
    }
    if j > 0 {
        let term1 = htr_recursive(g, n_ket, i + 1, j - 1, k, l, ab, cd);
        let term2 = htr_recursive(g, n_ket, i, j - 1, k, l, ab, cd);
        return term1 + ab * term2;
    }
    0.0
}

fn main() {
    println!("=== Debug p-orbital ERI computation ===\n");

    // All four centers at origin, all exponents = 1.0
    let center = [0.0, 0.0, 0.0];
    let alpha = 1.0;

    println!("Setup:");
    println!("  All centers at origin");
    println!("  All exponents = {}", alpha);
    println!();

    // Compute GaussianProduct2e
    let gp2e = GaussianProduct2e::new(
        alpha, &center, alpha, &center, alpha, &center, alpha, &center,
    );

    println!("GaussianProduct2e:");
    println!("  p = {}", gp2e.p);
    println!("  q = {}", gp2e.q);
    println!("  rho = {}", gp2e.rho);
    println!("  T = {}", gp2e.t);
    println!("  prefactor = {}", gp2e.prefactor);
    println!("  PA = {:?}", gp2e.pa);
    println!("  QC = {:?}", gp2e.qc);
    println!("  AB = {:?}", gp2e.ab);
    println!("  CD = {:?}", gp2e.cd);
    println!();

    // Verify prefactor formula
    let expected_prefactor = 2.0 * PI.powf(2.5) / (gp2e.p * gp2e.q * (gp2e.p + gp2e.q).sqrt());
    println!("Prefactor verification:");
    println!(
        "  Expected: 2 * pi^2.5 / (p * q * sqrt(p+q)) = {}",
        expected_prefactor
    );
    println!("  Actual: {}", gp2e.prefactor);
    println!();

    // For (pppp), nroots = (1+1+1+1)/2 + 1 = 3
    let nroots = 3;
    let n_bra = 2; // li + lj = 1 + 1
    let n_ket = 2; // lk + ll = 1 + 1

    println!("Rys quadrature:");
    println!("  nroots = {}", nroots);
    println!("  n_bra = {}", n_bra);
    println!("  n_ket = {}", n_ket);

    // Get Rys roots and weights
    let rys_result = rys_roots(nroots, gp2e.t).expect("Rys roots failed");
    println!("  Rys roots: {:?}", rys_result.roots);
    println!("  Rys weights: {:?}", rys_result.weights);
    println!(
        "  Sum of weights = {}",
        rys_result.weights.iter().sum::<f64>()
    );
    println!();

    // For (px,px|px,px), powers are (1,0,0) for all four
    let px = CartesianPower::new(1, 0, 0);

    println!("Computing (px,px|px,px):");
    println!("  Powers: i=(1,0,0), j=(1,0,0), k=(1,0,0), l=(1,0,0)");
    println!();

    let mut total_sum = 0.0;

    for (root_idx, &root) in rys_result.roots.iter().enumerate() {
        let weight = rys_result.weights[root_idx];

        println!("  Root {}: u = {}, w = {}", root_idx, root, weight);

        // Compute Rys coefficients
        let coeffs = RysCoefficients::compute(&gp2e, root);

        println!("    b00 = {}", coeffs.b00);
        println!("    b10 = {}", coeffs.b10);
        println!("    b01 = {}", coeffs.b01);
        println!("    c00 = {:?}", coeffs.c00);
        println!("    c0p = {:?}", coeffs.c0p);

        // Build 2D VRR tables for each direction
        let g_x = build_2d(
            n_bra,
            n_ket,
            coeffs.c00[0],
            coeffs.c0p[0],
            coeffs.b00,
            coeffs.b10,
            coeffs.b01,
        );
        let g_y = build_2d(
            n_bra,
            n_ket,
            coeffs.c00[1],
            coeffs.c0p[1],
            coeffs.b00,
            coeffs.b10,
            coeffs.b01,
        );
        let g_z = build_2d(
            n_bra,
            n_ket,
            coeffs.c00[2],
            coeffs.c0p[2],
            coeffs.b00,
            coeffs.b10,
            coeffs.b01,
        );

        println!("    g_x: {:?}", g_x);
        println!("    g_y: {:?}", g_y);
        println!("    g_z: {:?}", g_z);

        // Apply HTR to get [i,j,k,l] = [1,1,1,1] for x, [0,0,0,0] for y,z
        // For x: we need [1,1,1,1] from powers (1,0,0) on each center
        let i_x = htr_recursive(&g_x, n_ket, 1, 1, 1, 1, gp2e.ab[0], gp2e.cd[0]);
        let i_y = htr_recursive(&g_y, n_ket, 0, 0, 0, 0, gp2e.ab[1], gp2e.cd[1]);
        let i_z = htr_recursive(&g_z, n_ket, 0, 0, 0, 0, gp2e.ab[2], gp2e.cd[2]);

        println!("    I_x[1,1,1,1] = {}", i_x);
        println!("    I_y[0,0,0,0] = {}", i_y);
        println!("    I_z[0,0,0,0] = {}", i_z);

        let contribution = weight * i_x * i_y * i_z;
        println!("    Contribution = w * Ix * Iy * Iz = {}", contribution);

        total_sum += contribution;
        println!();
    }

    println!("Total sum over roots = {}", total_sum);

    // Now apply prefactor and normalization
    let raw_integral = gp2e.prefactor * total_sum;
    println!("Raw integral (prefactor * sum) = {}", raw_integral);

    // Cartesian normalization for p_x
    let cart_norm = (2.0 * alpha / PI).powf(0.75) * (4.0 * alpha).powf(0.5);
    println!("\nCartesian normalization:");
    println!("  cart_norm(px) = {}", cart_norm);
    println!("  cart_norm^4 = {}", cart_norm.powi(4));

    let normalized = raw_integral * cart_norm.powi(4);
    println!("\nFinal IQCP value = raw * norm^4 = {}", normalized);

    println!("\nPySCF reference = 0.9215096531280023");
    println!("Ratio = {}", normalized / 0.9215096531280023);
}
