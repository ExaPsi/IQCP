//! Generate IQCP spherical harmonic benchmark data for all 108 systems.
//!
//! Runs 6 molecules × 6 basis sets × 3 methods (RHF, LDA, B3LYP) with
//! spherical harmonic basis functions (use_spherical=true).
//!
//! Output: JSON array to stdout for comparison with PySCF.
//! Progress: printed to stderr.
//!
//! Run: cargo run --release --example benchmark_spherical

use qc_core::basis::{Atom, BasisSet};
use qc_core::dft::{build_becke_grid, ks_scf, B3lyp, GridConfig, GridQuality, Lda};
use qc_core::integrals;
use qc_core::scf::{rhf_scf, PresetSystem, ScfConfig};

// ============================================================================
// Molecule definitions (coordinates in bohr, matching scf_bench.rs)
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

// ============================================================================
// System builder using spherical harmonic integrals
// ============================================================================

fn build_system_spherical(atoms: &[Atom], basis_name: &str) -> (PresetSystem, BasisSet) {
    let basis = BasisSet::build(atoms.to_vec(), basis_name).unwrap();
    let nbf_sph = basis.n_basis_spherical();
    let nelec = basis.n_electrons;

    // Compute integrals in spherical harmonic basis
    let s_matrix = integrals::overlap_matrix_spherical(&basis);
    let h_core = integrals::hcore_matrix_spherical(&basis);
    let eri = integrals::eri_compressed_spherical(&basis);

    let system = PresetSystem {
        system_id: "bench_sph".to_string(),
        label: "Benchmark (spherical)".to_string(),
        nbf: nbf_sph,
        nelec,
        e_nuc: basis.nuclear_repulsion,
        s_matrix,
        h_core,
        eri_compressed: eri,
    };

    (system, basis)
}

// ============================================================================
// Molecule/basis/method definitions
// ============================================================================

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

const BASIS_SETS: [&str; 6] = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];

fn main() {
    let config = ScfConfig {
        use_diis: true,
        ..ScfConfig::tight()
    };
    let grid_config = GridConfig {
        n_radial: 75,
        quality: GridQuality::Standard,
        pruning: true,
    };

    let total = MOLECULES.len() * BASIS_SETS.len() * 3;
    let mut count = 0;

    // Collect all results into a vector, then serialize at the end
    let mut results: Vec<String> = Vec::with_capacity(total);

    for mol in &MOLECULES {
        let atoms = (mol.atoms_fn)();
        for basis_name in &BASIS_SETS {
            eprintln!(
                "[{}/{}] Building {} / {} (spherical integrals)...",
                count + 1,
                total,
                mol.name,
                basis_name
            );

            let (system, basis) = build_system_spherical(&atoms, basis_name);
            let nbf_sph = system.nbf;
            let nbf_cart = basis.n_basis;

            eprintln!(
                "  nbf_cart={}, nbf_sph={}, nelec={}, has_sph_diff={}",
                nbf_cart,
                nbf_sph,
                system.nelec,
                basis.has_spherical_difference()
            );

            // --- RHF ---
            count += 1;
            eprintln!("  [{}/{}] RHF...", count, total);
            let rhf_result = rhf_scf(&system, &config).unwrap();
            results.push(format!(
                "  {{\"molecule\":\"{}\",\"basis\":\"{}\",\"method\":\"RHF\",\"nbf_sph\":{},\"nbf_cart\":{},\"energy\":{:.12},\"converged\":{},\"iterations\":{}}}",
                mol.name, basis_name, nbf_sph, nbf_cart,
                rhf_result.energy_total, rhf_result.converged, rhf_result.iterations
            ));
            eprintln!(
                "    E={:.10} conv={} iter={}",
                rhf_result.energy_total, rhf_result.converged, rhf_result.iterations
            );

            // --- LDA ---
            count += 1;
            eprintln!("  [{}/{}] LDA...", count, total);
            let lda = Lda::new();
            let grid = build_becke_grid(&basis.atoms, &grid_config);
            let lda_result = ks_scf(&system, &config, &lda, &grid, &basis, true, None).unwrap();
            results.push(format!(
                "  {{\"molecule\":\"{}\",\"basis\":\"{}\",\"method\":\"LDA\",\"nbf_sph\":{},\"nbf_cart\":{},\"energy\":{:.12},\"converged\":{},\"iterations\":{}}}",
                mol.name, basis_name, nbf_sph, nbf_cart,
                lda_result.scf_output.energy_total, lda_result.scf_output.converged,
                lda_result.scf_output.iterations
            ));
            eprintln!(
                "    E={:.10} conv={} iter={}",
                lda_result.scf_output.energy_total,
                lda_result.scf_output.converged,
                lda_result.scf_output.iterations
            );

            // --- B3LYP ---
            count += 1;
            eprintln!("  [{}/{}] B3LYP...", count, total);
            let b3lyp = B3lyp::new();
            let b3lyp_result = ks_scf(&system, &config, &b3lyp, &grid, &basis, true, None).unwrap();
            results.push(format!(
                "  {{\"molecule\":\"{}\",\"basis\":\"{}\",\"method\":\"B3LYP\",\"nbf_sph\":{},\"nbf_cart\":{},\"energy\":{:.12},\"converged\":{},\"iterations\":{}}}",
                mol.name, basis_name, nbf_sph, nbf_cart,
                b3lyp_result.scf_output.energy_total, b3lyp_result.scf_output.converged,
                b3lyp_result.scf_output.iterations
            ));
            eprintln!(
                "    E={:.10} conv={} iter={}",
                b3lyp_result.scf_output.energy_total,
                b3lyp_result.scf_output.converged,
                b3lyp_result.scf_output.iterations
            );
        }
    }

    // Output JSON to stdout
    println!("[");
    for (i, entry) in results.iter().enumerate() {
        if i < results.len() - 1 {
            println!("{},", entry);
        } else {
            println!("{}", entry);
        }
    }
    println!("]");

    eprintln!("\nDone! {} systems completed.", count);
}
