// Test H2O 6-31G* ERI values
use qc_core::basis::{Atom, BasisSet};
use qc_core::integrals::{eri_compressed, eri_get};

fn main() {
    let o = Atom::new(8, [0.0, 0.0, 0.1173]).unwrap();
    let h1 = Atom::new(1, [0.0, 0.7572, -0.4692]).unwrap();
    let h2 = Atom::new(1, [0.0, -0.7572, -0.4692]).unwrap();
    let basis = BasisSet::build(vec![o, h1, h2], "6-31g*").unwrap();

    let eri = eri_compressed(&basis);
    let n = basis.n_basis;

    println!("H2O 6-31G* IQCP ERI values:");
    println!("Basis functions: {}", n);

    println!("\nDiagonal ERIs:");
    for i in 0..n {
        let val = eri_get(&eri, n, i, i, i, i);
        println!("  ({i},{i}|{i},{i}) = {val:.12e}");
    }

    println!("\nOff-diagonal ERIs:");
    println!("  (0,0|1,1) = {:.12e}", eri_get(&eri, n, 0, 0, 1, 1));
    println!("  (0,0|2,2) = {:.12e}", eri_get(&eri, n, 0, 0, 2, 2));
    println!("  (0,1|0,1) = {:.12e}", eri_get(&eri, n, 0, 1, 0, 1));
    println!("  (1,2|1,2) = {:.12e}", eri_get(&eri, n, 1, 2, 1, 2));
    println!("  (0,0|9,9) = {:.12e}", eri_get(&eri, n, 0, 0, 9, 9));
    println!(
        "  (15,15|17,17) = {:.12e}",
        eri_get(&eri, n, 15, 15, 17, 17)
    );
}
