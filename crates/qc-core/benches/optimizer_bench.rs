//! Criterion benchmarks for geometry optimization
//!
//! Measures full optimization trajectory timings:
//! - H2O RHF/STO-3G optimization (L-BFGS)
//! - H2O B3LYP/6-31G* optimization (L-BFGS)
//!
//! These are expensive benchmarks (multiple SCF + gradient evaluations per
//! optimization step), so sample sizes are reduced accordingly.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use qc_core::optimizer::{optimize_geometry, OptMethod, OptimizationConfig};

// ============================================================================
// Helpers: molecule geometries as (atomic_number, [x, y, z]) tuples
// ============================================================================

fn h2o_atoms() -> Vec<(u8, [f64; 3])> {
    vec![
        (8, [0.0, 0.0, 0.2217282]),
        (1, [0.0, 1.4305447, -0.8869128]),
        (1, [0.0, -1.4305447, -0.8869128]),
    ]
}

// ============================================================================
// Optimization benchmarks
// ============================================================================

/// Geometry optimization: H2O RHF/STO-3G
fn bench_optimize_h2o_rhf_sto3g(c: &mut Criterion) {
    let atoms = h2o_atoms();
    let config = OptimizationConfig {
        method: OptMethod::Rhf,
        basis: "sto-3g".to_string(),
        max_steps: 50,
        ..Default::default()
    };

    let mut group = c.benchmark_group("optimizer");
    group.sample_size(10);
    group.bench_function("optimize_h2o_rhf_sto3g", |b| {
        b.iter(|| optimize_geometry(black_box(&atoms), black_box(&config), None))
    });
    group.finish();
}

/// Geometry optimization: H2O B3LYP/6-31G*
fn bench_optimize_h2o_b3lyp_631gs(c: &mut Criterion) {
    let atoms = h2o_atoms();
    let config = OptimizationConfig {
        method: OptMethod::B3lyp,
        basis: "6-31g*".to_string(),
        max_steps: 50,
        ..Default::default()
    };

    let mut group = c.benchmark_group("optimizer");
    group.sample_size(10);
    group.bench_function("optimize_h2o_b3lyp_631gs", |b| {
        b.iter(|| optimize_geometry(black_box(&atoms), black_box(&config), None))
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_optimize_h2o_rhf_sto3g,
    bench_optimize_h2o_b3lyp_631gs,
);
criterion_main!(benches);
