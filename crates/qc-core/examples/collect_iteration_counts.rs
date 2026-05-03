//! Collect SCF iteration counts for all 108 benchmark systems.
//! Outputs JSON to stdout for use in benchmark_comparison.json.
//!
//! Run: cargo run --release --example collect_iteration_counts

use qc_core::basis::{Atom, BasisSet};
use qc_core::dft::{build_becke_grid, ks_scf, B3lyp, GridConfig, GridQuality, Lda};
use qc_core::integrals;
use qc_core::scf::{rhf_scf, PresetSystem, ScfConfig};

fn build_system(atoms: &[Atom], basis_name: &str) -> (PresetSystem, BasisSet) {
    let basis = BasisSet::build(atoms.to_vec(), basis_name).unwrap();
    let nbf = basis.n_basis;
    let nelec = basis.n_electrons;
    let s_matrix = integrals::overlap_matrix(&basis);
    let h_core = integrals::hcore_matrix(&basis);
    let eri = integrals::eri_compressed(&basis);
    let system = PresetSystem {
        system_id: "bench".to_string(),
        label: "bench".to_string(),
        nbf,
        nelec,
        e_nuc: basis.nuclear_repulsion,
        s_matrix,
        h_core,
        eri_compressed: eri,
    };
    (system, basis)
}

struct MolDef {
    name: &'static str,
    atoms: Vec<Atom>,
}

fn main() {
    let molecules = vec![
        MolDef {
            name: "H2",
            atoms: vec![
                Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
                Atom::new(1, [0.0, 0.0, 1.4]).unwrap(),
            ],
        },
        MolDef {
            name: "HF",
            atoms: vec![
                Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
                Atom::new(9, [0.0, 0.0, 1.7328]).unwrap(),
            ],
        },
        MolDef {
            name: "H2O",
            atoms: vec![
                Atom::new(8, [0.0, 0.0, 0.2217282]).unwrap(),
                Atom::new(1, [0.0, 1.4305447, -0.8869128]).unwrap(),
                Atom::new(1, [0.0, -1.4305447, -0.8869128]).unwrap(),
            ],
        },
        MolDef {
            name: "NH3",
            atoms: vec![
                Atom::new(7, [0.0, 0.0, 0.219705]).unwrap(),
                Atom::new(1, [0.0, 1.7714918, -0.512645]).unwrap(),
                Atom::new(1, [1.5342036, -0.8857459, -0.512645]).unwrap(),
                Atom::new(1, [-1.5342036, -0.8857459, -0.512645]).unwrap(),
            ],
        },
        MolDef {
            name: "CH4",
            atoms: vec![
                Atom::new(6, [0.0, 0.0, 0.0]).unwrap(),
                Atom::new(1, [1.1851, 1.1851, 1.1851]).unwrap(),
                Atom::new(1, [-1.1851, -1.1851, 1.1851]).unwrap(),
                Atom::new(1, [-1.1851, 1.1851, -1.1851]).unwrap(),
                Atom::new(1, [1.1851, -1.1851, -1.1851]).unwrap(),
            ],
        },
        MolDef {
            name: "C6H6",
            atoms: vec![
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
            ],
        },
    ];

    let bases = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];
    let config = ScfConfig {
        use_diis: true,
        ..ScfConfig::tight()
    };
    let grid_config = GridConfig {
        n_radial: 75,
        quality: GridQuality::Standard,
        pruning: true,
    };

    println!("{{");
    println!("  \"description\": \"IQCP iteration counts for all 108 benchmark systems\",");
    println!("  \"settings\": \"ScfConfig::tight(), DIIS enabled, grid: 75/Standard/SG-1\",");
    println!("  \"results\": [");

    let mut first = true;
    for mol in &molecules {
        for basis_name in &bases {
            let (system, basis) = build_system(&mol.atoms, basis_name);
            let n_occ = system.n_occ();

            // RHF
            let rhf_result = rhf_scf(&system, &config).unwrap();
            if !first {
                print!(",\n");
            }
            first = false;
            print!(
                "    {{\"molecule\":\"{}\",\"basis\":\"{}\",\"method\":\"RHF\",\"n_basis\":{},\"iterations\":{},\"energy\":{:.12},\"converged\":{}}}",
                mol.name, basis_name, system.nbf, rhf_result.iterations, rhf_result.energy_total, rhf_result.converged
            );

            // LDA
            let lda = Lda::new();
            let grid = build_becke_grid(&basis.atoms, &grid_config);
            let lda_result = ks_scf(&system, &config, &lda, &grid, &basis, false, None).unwrap();
            print!(
                ",\n    {{\"molecule\":\"{}\",\"basis\":\"{}\",\"method\":\"LDA\",\"n_basis\":{},\"iterations\":{},\"energy\":{:.12},\"converged\":{}}}",
                mol.name, basis_name, system.nbf, lda_result.scf_output.iterations, lda_result.scf_output.energy_total, lda_result.scf_output.converged
            );

            // B3LYP
            let b3lyp = B3lyp::new();
            let b3lyp_result =
                ks_scf(&system, &config, &b3lyp, &grid, &basis, false, None).unwrap();
            print!(
                ",\n    {{\"molecule\":\"{}\",\"basis\":\"{}\",\"method\":\"B3LYP\",\"n_basis\":{},\"iterations\":{},\"energy\":{:.12},\"converged\":{}}}",
                mol.name, basis_name, system.nbf, b3lyp_result.scf_output.iterations, b3lyp_result.scf_output.energy_total, b3lyp_result.scf_output.converged
            );
        }
    }

    println!("\n  ]");
    println!("}}");
}
