// Test C2H4 6-31G* ERI values
use qc_core::basis::{Atom, BasisSet};
use qc_core::integrals::{eri_compressed, eri_get};

fn main() {
    // C2H4 geometry in Bohr (same as test)
    let c1 = Atom::new(6, [0.0, 0.0, 1.2654]).unwrap();
    let c2 = Atom::new(6, [0.0, 0.0, -1.2654]).unwrap();
    let h1 = Atom::new(1, [0.0, 1.7453, 2.3280]).unwrap();
    let h2 = Atom::new(1, [0.0, -1.7453, 2.3280]).unwrap();
    let h3 = Atom::new(1, [0.0, 1.7453, -2.3280]).unwrap();
    let h4 = Atom::new(1, [0.0, -1.7453, -2.3280]).unwrap();

    let atoms = vec![c1, c2, h1, h2, h3, h4];
    let basis = BasisSet::build(atoms, "6-31g*").unwrap();

    let eri = eri_compressed(&basis);
    let n = basis.n_basis;

    println!("C2H4 6-31G* IQCP ERI values:");
    println!("Basis functions: {}", n);

    println!("\nAll d-orbital related ERIs from the test:");
    println!("  (0,0|9,9) = {:.12e}", eri_get(&eri, n, 0, 0, 9, 9));
    println!("  (9,9|12,12) = {:.12e}", eri_get(&eri, n, 9, 9, 12, 12));
    println!(
        "  (10,10|11,11) = {:.12e}",
        eri_get(&eri, n, 10, 10, 11, 11)
    );
    println!("  (9,9|24,24) = {:.12e}", eri_get(&eri, n, 9, 9, 24, 24));
    println!(
        "  (10,10|25,25) = {:.12e}",
        eri_get(&eri, n, 10, 10, 25, 25)
    );
    println!("  (9,9|10,10) = {:.12e}", eri_get(&eri, n, 9, 9, 10, 10));
    println!("  (9,9|30,30) = {:.12e}", eri_get(&eri, n, 9, 9, 30, 30));
    println!(
        "  (14,14|31,31) = {:.12e}",
        eri_get(&eri, n, 14, 14, 31, 31)
    );

    // Diagonal d-orbitals for verification
    println!("\nDiagonal d-orbital ERIs:");
    println!("  (9,9|9,9) = {:.12e}", eri_get(&eri, n, 9, 9, 9, 9));
    println!(
        "  (10,10|10,10) = {:.12e}",
        eri_get(&eri, n, 10, 10, 10, 10)
    );
    println!(
        "  (24,24|24,24) = {:.12e}",
        eri_get(&eri, n, 24, 24, 24, 24)
    );
}
