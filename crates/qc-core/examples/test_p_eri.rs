// Quick test for p-orbital ERI
use qc_core::basis::{AngularMomentum, ContractedShell, GaussianPrimitive};
use qc_core::integrals::shell_eri;

fn main() {
    // Single primitive p-shell at origin with exp=1.0, coef=1.0
    // This should match PySCF's value of 0.9215096531280023
    let shell_p = ContractedShell::new(
        AngularMomentum::P,
        vec![GaussianPrimitive::new(1.0, 1.0)],
        [0.0, 0.0, 0.0],
        0,
    );

    let result = shell_eri(&shell_p, &shell_p, &shell_p, &shell_p);

    println!("Single primitive p-shell (exp=1.0, coef=1.0):");
    println!("  (px,px|px,px) = {}", result.get(0, 0, 0, 0));
    println!("  (py,py|py,py) = {}", result.get(1, 1, 1, 1));
    println!("  (pz,pz|pz,pz) = {}", result.get(2, 2, 2, 2));
    println!("  (px,py|px,py) = {}", result.get(0, 1, 0, 1));

    println!("\nPySCF reference:");
    println!("  (px,px|px,px) = 0.9215096531280023");
    println!("  (px,py|px,py) = 0.056418958354775645");

    let iqcp = result.get(0, 0, 0, 0);
    let pyscf = 0.9215096531280023;
    println!("\nRatio IQCP/PySCF = {}", iqcp / pyscf);
}
