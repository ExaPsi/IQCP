//! Criterion benchmarks for the Boys function F_m(T)
//!
//! Measures evaluation time across computational regimes:
//! - Series expansion (T < turnover)
//! - erf + upward recurrence (T >= turnover)
//! - High-order evaluation
//! - Sweep over many T values

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use qc_core::boys::{boys_eval, boys_eval_many};

/// F_0(0.5): series expansion regime (small T)
fn bench_boys_single_series(c: &mut Criterion) {
    c.bench_function("boys_single_series [F_0(0.5)]", |b| {
        b.iter(|| boys_eval(black_box(0), black_box(0.5)).unwrap())
    });
}

/// F_0(20.0): erf + upward recurrence regime (moderate T)
fn bench_boys_single_recurrence(c: &mut Criterion) {
    c.bench_function("boys_single_recurrence [F_0(20.0)]", |b| {
        b.iter(|| boys_eval(black_box(0), black_box(20.0)).unwrap())
    });
}

/// F_0(50.0): asymptotic regime (large T)
fn bench_boys_single_asymptotic(c: &mut Criterion) {
    c.bench_function("boys_single_asymptotic [F_0(50.0)]", |b| {
        b.iter(|| boys_eval(black_box(0), black_box(50.0)).unwrap())
    });
}

/// F_10(5.0): high-order evaluation
fn bench_boys_high_order(c: &mut Criterion) {
    c.bench_function("boys_high_order [F_10(5.0)]", |b| {
        b.iter(|| boys_eval(black_box(10), black_box(5.0)).unwrap())
    });
}

/// boys_eval_many(5, &ts) where ts = 100 values spanning 0.0 to 100.0
fn bench_boys_sweep(c: &mut Criterion) {
    let ts: Vec<f64> = (0..100).map(|i| i as f64 * 1.0).collect();
    c.bench_function("boys_sweep [m=5, 100 pts, T=0..100]", |b| {
        b.iter(|| boys_eval_many(black_box(5), black_box(&ts)).unwrap())
    });
}

criterion_group!(
    benches,
    bench_boys_single_series,
    bench_boys_single_recurrence,
    bench_boys_single_asymptotic,
    bench_boys_high_order,
    bench_boys_sweep,
);
criterion_main!(benches);
