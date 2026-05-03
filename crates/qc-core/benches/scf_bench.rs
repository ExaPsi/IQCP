//! Criterion benchmarks for SCF (Self-Consistent Field) calculations
//!
//! Benchmarks all 108 combinations: 6 molecules × 6 basis sets × 3 methods
//! All systems compute integrals on-the-fly for fair comparison with PySCF.
//!
//! Molecules: H2, HF, H2O, NH3, CH4, C6H6
//! Basis sets: STO-3G, 3-21G, 6-31G, 6-31G*, 6-31+G*, cc-pVDZ
//! Methods: RHF, LDA, B3LYP

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use qc_core::basis::{Atom, BasisSet};
use qc_core::dft::{build_becke_grid, ks_scf, B3lyp, GridConfig, GridQuality, Lda};
use qc_core::integrals;
use qc_core::scf::{rhf_scf, PresetSystem, ScfConfig};
use std::time::Duration;

// ============================================================================
// Molecule definitions (coordinates in bohr, matching validation.json)
// ============================================================================

fn h2_atoms() -> Vec<Atom> {
    vec![
        Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
        Atom::new(1, [0.0, 0.0, 1.4]).unwrap(),
    ]
}

fn hf_atoms() -> Vec<Atom> {
    vec![
        Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
        Atom::new(9, [0.0, 0.0, 1.7328]).unwrap(),
    ]
}

fn h2o_atoms() -> Vec<Atom> {
    vec![
        Atom::new(8, [0.0, 0.0, 0.2217282]).unwrap(),
        Atom::new(1, [0.0, 1.4305447, -0.8869128]).unwrap(),
        Atom::new(1, [0.0, -1.4305447, -0.8869128]).unwrap(),
    ]
}

fn nh3_atoms() -> Vec<Atom> {
    vec![
        Atom::new(7, [0.0, 0.0, 0.219705]).unwrap(),
        Atom::new(1, [0.0, 1.7714918, -0.512645]).unwrap(),
        Atom::new(1, [1.5342036, -0.8857459, -0.512645]).unwrap(),
        Atom::new(1, [-1.5342036, -0.8857459, -0.512645]).unwrap(),
    ]
}

fn ch4_atoms() -> Vec<Atom> {
    vec![
        Atom::new(6, [0.0, 0.0, 0.0]).unwrap(),
        Atom::new(1, [1.1851, 1.1851, 1.1851]).unwrap(),
        Atom::new(1, [-1.1851, -1.1851, 1.1851]).unwrap(),
        Atom::new(1, [-1.1851, 1.1851, -1.1851]).unwrap(),
        Atom::new(1, [1.1851, -1.1851, -1.1851]).unwrap(),
    ]
}

// ============================================================================
// System builder
// ============================================================================

fn build_system(atoms: &[Atom], basis_name: &str) -> (PresetSystem, BasisSet) {
    let basis = BasisSet::build(atoms.to_vec(), basis_name).unwrap();
    let nbf = basis.n_basis;
    let nelec = basis.n_electrons;

    let s_matrix = integrals::overlap_matrix(&basis);
    let h_core = integrals::hcore_matrix(&basis);
    let eri = integrals::eri_compressed(&basis);

    let system = PresetSystem {
        system_id: "bench".to_string(),
        label: "Benchmark system".to_string(),
        nbf,
        nelec,
        e_nuc: basis.nuclear_repulsion,
        s_matrix,
        h_core,
        eri_compressed: eri,
    };

    (system, basis)
}

// ============================================================================
// Benchmark definitions
// ============================================================================

struct MolDef {
    name: &'static str,
    atoms_fn: fn() -> Vec<Atom>,
}

fn c6h6_atoms() -> Vec<Atom> {
    vec![
        Atom::new(6, [0.0000000000, 2.6399473960, 0.0000000000]).unwrap(),
        Atom::new(6, [2.2861906655, 1.3199736980, 0.0000000000]).unwrap(),
        Atom::new(6, [2.2861906655, -1.3199736980, 0.0000000000]).unwrap(),
        Atom::new(6, [0.0000000000, -2.6399473960, 0.0000000000]).unwrap(),
        Atom::new(6, [-2.2861906655, -1.3199736980, 0.0000000000]).unwrap(),
        Atom::new(6, [-2.2861906655, 1.3199736980, 0.0000000000]).unwrap(),
        Atom::new(1, [0.0000000000, 4.6884105150, 0.0000000000]).unwrap(),
        Atom::new(1, [4.0602655512, 2.3442052575, 0.0000000000]).unwrap(),
        Atom::new(1, [4.0602655512, -2.3442052575, 0.0000000000]).unwrap(),
        Atom::new(1, [0.0000000000, -4.6884105150, 0.0000000000]).unwrap(),
        Atom::new(1, [-4.0602655512, -2.3442052575, 0.0000000000]).unwrap(),
        Atom::new(1, [-4.0602655512, 2.3442052575, 0.0000000000]).unwrap(),
    ]
}

const MOLECULES: [MolDef; 6] = [
    MolDef {
        name: "H2",
        atoms_fn: h2_atoms,
    },
    MolDef {
        name: "HF",
        atoms_fn: hf_atoms,
    },
    MolDef {
        name: "H2O",
        atoms_fn: h2o_atoms,
    },
    MolDef {
        name: "NH3",
        atoms_fn: nh3_atoms,
    },
    MolDef {
        name: "CH4",
        atoms_fn: ch4_atoms,
    },
    MolDef {
        name: "C6H6",
        atoms_fn: c6h6_atoms,
    },
];

const BASIS_SETS: [&str; 6] = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];

// ============================================================================
// RHF benchmarks: 5 molecules × 2 basis sets = 10
// ============================================================================

fn bench_rhf(c: &mut Criterion) {
    let config = ScfConfig {
        use_diis: true,
        ..ScfConfig::tight()
    };

    let mut group = c.benchmark_group("rhf");
    group.warm_up_time(Duration::from_secs(3));

    for mol in &MOLECULES {
        for basis_name in &BASIS_SETS {
            let atoms = (mol.atoms_fn)();
            let (system, _basis) = build_system(&atoms, basis_name);
            let id = BenchmarkId::new(format!("{}/{}", mol.name, basis_name), system.nbf);

            group.bench_with_input(id, &(system, config.clone()), |b, (sys, cfg)| {
                b.iter(|| rhf_scf(black_box(sys), black_box(cfg)).unwrap())
            });
        }
    }
    group.finish();
}

// ============================================================================
// LDA benchmarks: 5 molecules × 2 basis sets = 10
// ============================================================================

fn bench_lda(c: &mut Criterion) {
    let config = ScfConfig {
        use_diis: true,
        ..ScfConfig::tight()
    };
    // Match PySCF default grid: 75 radial, Standard (Lebedev-302), SG-1 pruning
    let grid_config = GridConfig {
        n_radial: 75,
        quality: GridQuality::Standard,
        pruning: true,
    };
    let lda = Lda::new();

    let mut group = c.benchmark_group("lda");
    group.warm_up_time(Duration::from_secs(3));
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(30));

    for mol in &MOLECULES {
        for basis_name in &BASIS_SETS {
            let atoms = (mol.atoms_fn)();
            let (system, basis) = build_system(&atoms, basis_name);
            let grid = build_becke_grid(&basis.atoms, &grid_config);
            let id = BenchmarkId::new(format!("{}/{}", mol.name, basis_name), system.nbf);

            group.bench_with_input(id, &(), |b, _| {
                b.iter(|| {
                    ks_scf(
                        black_box(&system),
                        black_box(&config),
                        black_box(&lda),
                        black_box(&grid),
                        black_box(&basis),
                        false,
                        None,
                    )
                    .unwrap()
                })
            });
        }
    }
    group.finish();
}

// ============================================================================
// B3LYP benchmarks: 5 molecules × 2 basis sets = 10
// ============================================================================

fn bench_b3lyp(c: &mut Criterion) {
    let config = ScfConfig {
        use_diis: true,
        ..ScfConfig::tight()
    };
    // Match PySCF default grid: 75 radial, Standard (Lebedev-302), SG-1 pruning
    let grid_config = GridConfig {
        n_radial: 75,
        quality: GridQuality::Standard,
        pruning: true,
    };
    let b3lyp = B3lyp::new();

    let mut group = c.benchmark_group("b3lyp");
    group.warm_up_time(Duration::from_secs(3));
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(30));

    for mol in &MOLECULES {
        for basis_name in &BASIS_SETS {
            let atoms = (mol.atoms_fn)();
            let (system, basis) = build_system(&atoms, basis_name);
            let grid = build_becke_grid(&basis.atoms, &grid_config);
            let id = BenchmarkId::new(format!("{}/{}", mol.name, basis_name), system.nbf);

            group.bench_with_input(id, &(), |b, _| {
                b.iter(|| {
                    ks_scf(
                        black_box(&system),
                        black_box(&config),
                        black_box(&b3lyp),
                        black_box(&grid),
                        black_box(&basis),
                        false,
                        None,
                    )
                    .unwrap()
                })
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_rhf, bench_lda, bench_b3lyp);
criterion_main!(benches);
