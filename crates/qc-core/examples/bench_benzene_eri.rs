//! Benchmark ERI computation for benzene (C6H6) with 6-31G* basis
//!
//! Run with:
//!   cargo run --release --example bench_benzene_eri --features parallel
//!   cargo run --release --example bench_benzene_eri  # without parallel

use qc_core::basis::{Atom, BasisSet, ANGSTROM_TO_BOHR};
use std::time::Instant;

fn main() {
    println!("=== Benzene (C6H6) 6-31G* ERI Benchmark ===\n");

    // Benzene geometry in Angstroms
    let atoms: Vec<Atom> = vec![
        Atom::new(1, [1.2194, -0.1652, 2.1600].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(6, [0.6825, -0.0924, 1.2087].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(6, [-0.7075, -0.0352, 1.1973].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(1, [-1.2644, -0.0630, 2.1393].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(6, [-1.3898, 0.0572, -0.0114].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(1, [-2.4836, 0.1021, -0.0204].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(6, [-0.6824, 0.0925, -1.2088].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(1, [-1.2194, 0.1652, -2.1599].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(6, [0.7075, 0.0352, -1.1973].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(1, [1.2641, 0.0628, -2.1395].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(6, [1.3899, -0.0572, 0.0114].map(|x| x * ANGSTROM_TO_BOHR)),
        Atom::new(1, [2.4836, -0.1022, 0.0205].map(|x| x * ANGSTROM_TO_BOHR)),
    ]
    .into_iter()
    .collect::<Result<Vec<_>, _>>()
    .expect("Failed to create atoms");

    // Build basis set
    let basis = BasisSet::build(atoms, "6-31g*").expect("Failed to build basis set");

    println!("Basis functions: {}", basis.n_basis);
    println!("Number of shells: {}", basis.shells.len());

    let n = basis.n_basis;
    let n_pairs = n * (n + 1) / 2;
    let n_eri = n_pairs * (n_pairs + 1) / 2;
    println!("Unique ERIs (8-fold symmetry): {}", n_eri);

    // Count shell quartets
    let mut n_quartets = 0usize;
    for si in 0..basis.shells.len() {
        for sj in 0..=si {
            for sk in 0..basis.shells.len() {
                for sl in 0..=sk {
                    n_quartets += 1;
                }
            }
        }
    }
    println!("Shell quartets to compute: {}", n_quartets);

    #[cfg(feature = "parallel")]
    {
        println!("\n--- Parallel ERI (Rayon) ---");
        use qc_core::integrals::eri_compressed_parallel;

        // Warmup
        let _ = eri_compressed_parallel(&basis);

        // Timed run
        let start = Instant::now();
        let eri = eri_compressed_parallel(&basis);
        let elapsed = start.elapsed();

        println!("Time: {:.3} seconds", elapsed.as_secs_f64());
        println!("ERIs computed: {}", eri.len());
        println!(
            "Throughput: {:.2} shell quartets/second",
            n_quartets as f64 / elapsed.as_secs_f64()
        );
        println!(
            "Throughput: {:.2} million ERIs/second",
            eri.len() as f64 / elapsed.as_secs_f64() / 1_000_000.0
        );
    }

    #[cfg(not(feature = "parallel"))]
    {
        println!("\n--- Sequential ERI ---");
        use qc_core::integrals::eri_compressed;

        // Warmup
        let _ = eri_compressed(&basis);

        // Timed run
        let start = Instant::now();
        let eri = eri_compressed(&basis);
        let elapsed = start.elapsed();

        println!("Time: {:.3} seconds", elapsed.as_secs_f64());
        println!("ERIs computed: {}", eri.len());
        println!(
            "Throughput: {:.2} shell quartets/second",
            n_quartets as f64 / elapsed.as_secs_f64()
        );
        println!(
            "Throughput: {:.2} million ERIs/second",
            eri.len() as f64 / elapsed.as_secs_f64() / 1_000_000.0
        );
    }

    println!("\n=== Benchmark Complete ===");
}
