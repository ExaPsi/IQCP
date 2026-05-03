//! Criterion benchmarks for two-electron repulsion integral (ERI) computation
//!
//! Measures ERI evaluation for molecules of increasing size:
//! - H2/STO-3G (2 basis functions)
//! - H2O/STO-3G (7 basis functions)
//! - H2O/6-31G* (19 basis functions, Cartesian d)
//! - CH4/6-31G* (22 basis functions, Cartesian d)

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use qc_core::basis::{Atom, BasisSet};
use qc_core::integrals::eri_compressed;

/// Build H2 molecule (STO-3G, R=1.4 bohr)
fn h2_basis_sto3g() -> BasisSet {
    let atoms = vec![
        Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
        Atom::new(1, [0.0, 0.0, 1.4]).unwrap(),
    ];
    BasisSet::build(atoms, "sto-3g").unwrap()
}

/// Build H2O molecule (STO-3G)
fn h2o_basis_sto3g() -> BasisSet {
    let atoms = vec![
        Atom::new(8, [0.0, 0.0, 0.2217282]).unwrap(),
        Atom::new(1, [0.0, 1.4305447, -0.8869128]).unwrap(),
        Atom::new(1, [0.0, -1.4305447, -0.8869128]).unwrap(),
    ];
    BasisSet::build(atoms, "sto-3g").unwrap()
}

/// Build H2O molecule (6-31G*)
fn h2o_basis_631gs() -> BasisSet {
    let atoms = vec![
        Atom::new(8, [0.0, 0.0, 0.2217282]).unwrap(),
        Atom::new(1, [0.0, 1.4305447, -0.8869128]).unwrap(),
        Atom::new(1, [0.0, -1.4305447, -0.8869128]).unwrap(),
    ];
    BasisSet::build(atoms, "6-31g*").unwrap()
}

/// Build CH4 molecule (6-31G*)
fn ch4_basis_631gs() -> BasisSet {
    let atoms = vec![
        Atom::new(6, [0.0, 0.0, 0.0]).unwrap(),
        Atom::new(1, [1.1851, 1.1851, 1.1851]).unwrap(),
        Atom::new(1, [-1.1851, -1.1851, 1.1851]).unwrap(),
        Atom::new(1, [-1.1851, 1.1851, -1.1851]).unwrap(),
        Atom::new(1, [1.1851, -1.1851, -1.1851]).unwrap(),
    ];
    BasisSet::build(atoms, "6-31g*").unwrap()
}

/// ERI for H2/STO-3G (2 BFs, 6 unique integrals)
fn bench_eri_h2_sto3g(c: &mut Criterion) {
    let basis = h2_basis_sto3g();
    c.bench_function("eri_h2_sto3g [2 BFs]", |b| {
        b.iter(|| eri_compressed(black_box(&basis)))
    });
}

/// ERI for H2O/STO-3G (7 BFs)
fn bench_eri_h2o_sto3g(c: &mut Criterion) {
    let basis = h2o_basis_sto3g();
    c.bench_function("eri_h2o_sto3g [7 BFs]", |b| {
        b.iter(|| eri_compressed(black_box(&basis)))
    });
}

/// ERI for H2O/6-31G* (Cartesian d-functions)
fn bench_eri_h2o_631gs(c: &mut Criterion) {
    let basis = h2o_basis_631gs();
    c.bench_function(&format!("eri_h2o_631gs [{} BFs]", basis.n_basis), |b| {
        b.iter(|| eri_compressed(black_box(&basis)))
    });
}

/// ERI for CH4/6-31G* (Cartesian d-functions)
fn bench_eri_ch4_631gs(c: &mut Criterion) {
    let basis = ch4_basis_631gs();
    c.bench_function(&format!("eri_ch4_631gs [{} BFs]", basis.n_basis), |b| {
        b.iter(|| eri_compressed(black_box(&basis)))
    });
}

criterion_group!(
    benches,
    bench_eri_h2_sto3g,
    bench_eri_h2o_sto3g,
    bench_eri_h2o_631gs,
    bench_eri_ch4_631gs,
);
criterion_main!(benches);
