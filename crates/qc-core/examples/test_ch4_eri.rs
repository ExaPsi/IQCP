// Test CH4 6-31G ERI values
use qc_core::basis::{Atom, BasisSet};
use qc_core::integrals::eri_compressed;

fn main() {
    // CH4 geometry in Bohr
    let c = Atom::new(6, [0.0, 0.0, 0.0]).unwrap();
    let h1 = Atom::new(1, [1.19, 1.19, 1.19]).unwrap();
    let h2 = Atom::new(1, [-1.19, -1.19, 1.19]).unwrap();
    let h3 = Atom::new(1, [-1.19, 1.19, -1.19]).unwrap();
    let h4 = Atom::new(1, [1.19, -1.19, -1.19]).unwrap();

    let atoms = vec![c, h1, h2, h3, h4];
    let basis = BasisSet::build(atoms, "6-31g").expect("Failed to build basis");

    println!("=== CH4 6-31G ERI Check ===");
    println!("Number of basis functions: {}", basis.n_basis);

    let eri = eri_compressed(&basis);
    println!("Number of ERIs: {}", eri.len());

    println!("\nFirst 10 compressed ERIs:");
    for (i, val) in eri.iter().take(10).enumerate() {
        println!("  eri[{}] = {:.12e}", i, val);
    }

    println!("\nERIs at indices 0-5:");
    for i in 0..6 {
        println!("  eri[{}] = {:.12e}", i, eri[i]);
    }

    // Reference from PySCF:
    // eri[0] = (0,0|0,0) = 3.534810920238e+00
    // eri[1] = (0,0|0,1) = 4.003061234815e-01
    // eri[2] = (0,1|0,1) = 7.145087334305e-02
    println!("\n=== Comparison with PySCF ===");
    let pyscf_eri = [
        3.53481092e+00,
        4.00306123e-01,
        7.14508733e-02,
        9.50464728e-01,
        2.00655930e-01,
        7.49198320e-01,
    ];

    for i in 0..6 {
        let diff = (eri[i] - pyscf_eri[i]).abs();
        println!(
            "  eri[{}]: IQCP={:.10e}, PySCF={:.10e}, diff={:.2e}",
            i, eri[i], pyscf_eri[i], diff
        );
    }

    // The most important one: (0,0|0,0) should be ~3.53
    println!("\n=== Critical check: (0,0|0,0) ===");
    println!("IQCP:  {:.12e}", eri[0]);
    println!("PySCF: 3.534810920238e+00");
    println!("This ERI is for the C 1s orbital.");
}
