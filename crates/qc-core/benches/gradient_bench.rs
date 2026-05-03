//! Criterion benchmarks for gradient computation (SCF + gradient)
//!
//! Measures full gradient time: SCF convergence + analytical gradient.
//! 6 molecules × B3LYP/6-31G* (matching PySCF comparison).

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use qc_core::basis::Atom;
use qc_core::dft::{B3lyp, GridConfig, GridQuality};
use qc_core::scf::gradient::ks_dft_gradient;
use qc_core::scf::ScfConfig;
use std::time::Duration;

// ============================================================================
// Molecule definitions (coordinates in bohr)
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

fn c6h6_atoms() -> Vec<Atom> {
    vec![
        Atom::new(6, [0.0, 2.6399473960, 0.0]).unwrap(),
        Atom::new(6, [2.2861906655, 1.3199736980, 0.0]).unwrap(),
        Atom::new(6, [2.2861906655, -1.3199736980, 0.0]).unwrap(),
        Atom::new(6, [0.0, -2.6399473960, 0.0]).unwrap(),
        Atom::new(6, [-2.2861906655, -1.3199736980, 0.0]).unwrap(),
        Atom::new(6, [-2.2861906655, 1.3199736980, 0.0]).unwrap(),
        Atom::new(1, [0.0, 4.6884105150, 0.0]).unwrap(),
        Atom::new(1, [4.0602655512, 2.3442052575, 0.0]).unwrap(),
        Atom::new(1, [4.0602655512, -2.3442052575, 0.0]).unwrap(),
        Atom::new(1, [0.0, -4.6884105150, 0.0]).unwrap(),
        Atom::new(1, [-4.0602655512, -2.3442052575, 0.0]).unwrap(),
        Atom::new(1, [-4.0602655512, 2.3442052575, 0.0]).unwrap(),
    ]
}

struct MolDef {
    name: &'static str,
    atoms_fn: fn() -> Vec<Atom>,
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

// ============================================================================
// B3LYP/6-31G* gradient benchmarks (SCF + gradient)
// ============================================================================

fn bench_b3lyp_gradient(c: &mut Criterion) {
    let b3lyp = B3lyp::new();
    // Match PySCF default grid
    let grid_config = GridConfig {
        n_radial: 75,
        quality: GridQuality::Standard,
        pruning: true,
    };
    let scf_config = ScfConfig {
        use_diis: true,
        ..ScfConfig::tight()
    };

    let mut group = c.benchmark_group("b3lyp_gradient");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(3));
    group.measurement_time(Duration::from_secs(30));

    for mol in &MOLECULES {
        let atoms = (mol.atoms_fn)();
        let basis = qc_core::basis::BasisSet::build(atoms.clone(), "6-31g*").unwrap();
        let nbf = basis.n_basis;

        let id = BenchmarkId::new(mol.name, nbf);

        group.bench_with_input(id, &atoms, |b, atoms| {
            b.iter(|| {
                ks_dft_gradient(
                    black_box(atoms),
                    black_box("6-31g*"),
                    black_box(&b3lyp),
                    black_box(&grid_config),
                    black_box(&scf_config),
                    black_box(false),
                )
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_b3lyp_gradient);
criterion_main!(benches);
