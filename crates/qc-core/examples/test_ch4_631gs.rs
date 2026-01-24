// Test CH4 6-31G* integrals and SCF (with d-polarization)
use qc_core::basis::{Atom, BasisSet};
use qc_core::integrals::{eri_compressed, hcore_matrix, overlap_matrix};
use qc_core::scf::{rhf_scf, ConvergenceProfile, PresetSystem, ScfConfig};

fn main() {
    // CH4 geometry in Bohr (same as IQCP exampleGeometries.ts)
    let c = Atom::new(6, [0.0, 0.0, 0.0]).unwrap();
    let h1 = Atom::new(1, [1.19, 1.19, 1.19]).unwrap();
    let h2 = Atom::new(1, [-1.19, -1.19, 1.19]).unwrap();
    let h3 = Atom::new(1, [-1.19, 1.19, -1.19]).unwrap();
    let h4 = Atom::new(1, [1.19, -1.19, -1.19]).unwrap();

    let atoms = vec![c, h1, h2, h3, h4];

    // Build 6-31G* basis (with d polarization on Carbon)
    let basis = BasisSet::build(atoms, "6-31g*").expect("Failed to build basis");

    println!("=== CH4 6-31G* Basis Info ===");
    println!("Number of basis functions: {}", basis.n_basis);
    println!("Number of electrons: {}", basis.n_electrons);
    println!("Number of shells: {}", basis.n_shells());
    println!("Nuclear repulsion: {:.12e}", basis.nuclear_repulsion);

    // Compute integrals
    println!("\n=== Computing Integrals ===");
    let s_matrix = overlap_matrix(&basis);
    let h_core = hcore_matrix(&basis);
    let eri = eri_compressed(&basis);

    println!(
        "S matrix size: {} (expected {})",
        s_matrix.len(),
        basis.n_basis * basis.n_basis
    );
    println!(
        "H_core size: {} (expected {})",
        h_core.len(),
        basis.n_basis * basis.n_basis
    );

    let n_pairs = basis.n_basis * (basis.n_basis + 1) / 2;
    let expected_eri = n_pairs * (n_pairs + 1) / 2;
    println!("ERI size: {} (expected {})", eri.len(), expected_eri);

    // Print traces
    let s_trace: f64 = (0..basis.n_basis)
        .map(|i| s_matrix[i * basis.n_basis + i])
        .sum();
    let h_trace: f64 = (0..basis.n_basis)
        .map(|i| h_core[i * basis.n_basis + i])
        .sum();
    println!("\nS trace: {:.6}", s_trace);
    println!("H_core trace: {:.6}", h_trace);

    // Check some matrix elements
    println!("\n=== Overlap Matrix (first 5x5) ===");
    let nbf = basis.n_basis;
    for i in 0..5.min(nbf) {
        for j in 0..5.min(nbf) {
            print!("{:12.6}", s_matrix[i * nbf + j]);
        }
        println!();
    }

    println!("\n=== H_core Matrix (first 5x5) ===");
    for i in 0..5.min(nbf) {
        for j in 0..5.min(nbf) {
            print!("{:12.6}", h_core[i * nbf + j]);
        }
        println!();
    }

    // Create preset system for SCF
    let system = PresetSystem {
        system_id: "ch4_631gs".to_string(),
        label: "CH4 (6-31G*)".to_string(),
        nbf: basis.n_basis,
        nelec: basis.n_electrons,
        e_nuc: basis.nuclear_repulsion,
        s_matrix,
        h_core,
        eri_compressed: eri,
    };

    // Run SCF
    println!("\n=== Running SCF ===");

    let mut config = ScfConfig::new(ConvergenceProfile::Tight);
    config.max_iterations = 200;
    config.use_diis = true;
    config.diis_start = 3;
    config.diis_size = 8;

    match rhf_scf(&system, &config) {
        Ok(result) => {
            println!("Converged: {}", result.converged);
            println!("Iterations: {}", result.iterations);
            println!("Final energy: {:.12} Ha", result.energy_total);
            println!("Electronic energy: {:.12} Ha", result.energy_electronic);

            println!("\nPySCF reference: -40.19484951 Ha");
            println!(
                "Difference: {:.8e} Ha",
                result.energy_total - (-40.19484951)
            );
        }
        Err(e) => {
            println!("SCF failed: {}", e);
        }
    }
}
