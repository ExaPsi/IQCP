// Debug example for d-orbital ERI computation
// Compares IQCP computation with PySCF reference values

use qc_core::basis::{AngularMomentum, ContractedShell, GaussianPrimitive};
use qc_core::integrals::shell_eri;

fn main() {
    println!("=== Debug d-orbital ERI computation ===\n");

    // Single primitive d-shell at origin with 6-31G* exponent for O
    // From 6-31G*: O d exponent = 0.8 (polarization)
    let d_exp = 0.8;
    let shell_d = ContractedShell::new(
        AngularMomentum::D,
        vec![GaussianPrimitive::new(d_exp, 1.0)],
        [0.0, 0.0, 0.0],
        0,
    );

    println!("D-shell configuration:");
    println!("  Angular momentum: D (L=2)");
    println!("  Exponent: {}", d_exp);
    println!("  Center: origin");
    println!();

    let result = shell_eri(&shell_d, &shell_d, &shell_d, &shell_d);

    println!("IQCP computed ERIs:");
    // Cartesian ordering: dxx, dxy, dxz, dyy, dyz, dzz
    println!("  (dxx,dxx|dxx,dxx) = {}", result.get(0, 0, 0, 0));
    println!("  (dxy,dxy|dxy,dxy) = {}", result.get(1, 1, 1, 1));
    println!("  (dxz,dxz|dxz,dxz) = {}", result.get(2, 2, 2, 2));
    println!("  (dyy,dyy|dyy,dyy) = {}", result.get(3, 3, 3, 3));
    println!("  (dyz,dyz|dyz,dyz) = {}", result.get(4, 4, 4, 4));
    println!("  (dzz,dzz|dzz,dzz) = {}", result.get(5, 5, 5, 5));

    println!("\nPySCF reference (d shell at origin, exp=0.8):");
    // Need to compute these with PySCF
    // For now, let's compute ratio with expected values based on test

    // Also test with exp=1.0 for simpler comparison
    println!("\n\n=== With exp=1.0 ===");
    let shell_d_1 = ContractedShell::new(
        AngularMomentum::D,
        vec![GaussianPrimitive::new(1.0, 1.0)],
        [0.0, 0.0, 0.0],
        0,
    );

    let result_1 = shell_eri(&shell_d_1, &shell_d_1, &shell_d_1, &shell_d_1);

    println!("IQCP computed ERIs (exp=1.0):");
    println!("  (dxx,dxx|dxx,dxx) = {}", result_1.get(0, 0, 0, 0));
    println!("  (dxy,dxy|dxy,dxy) = {}", result_1.get(1, 1, 1, 1));
    println!("  (dxz,dxz|dxz,dxz) = {}", result_1.get(2, 2, 2, 2));
    println!("  (dyy,dyy|dyy,dyy) = {}", result_1.get(3, 3, 3, 3));
    println!("  (dyz,dyz|dyz,dyz) = {}", result_1.get(4, 4, 4, 4));
    println!("  (dzz,dzz|dzz,dzz) = {}", result_1.get(5, 5, 5, 5));

    println!("\nPySCF reference (exp=1.0, Cartesian):");
    println!("  (dxx,dxx|dxx,dxx) = 5.396976697204335");
    println!("  (dxy,dxy|dxy,dxy) = 0.5397448090050628");
    println!("  (dxz,dxz|dxz,dxz) = 0.5397448090050628");
    println!("  (dyy,dyy|dyy,dyy) = 5.396976697204335");
    println!("  (dyz,dyz|dyz,dyz) = 0.5397448090050628");
    println!("  (dzz,dzz|dzz,dzz) = 5.396976697204336");

    let iqcp_dxx = result_1.get(0, 0, 0, 0);
    let pyscf_dxx = 5.396976697204335;
    println!(
        "\nRatio IQCP/PySCF for (dxx,dxx|dxx,dxx): {}",
        iqcp_dxx / pyscf_dxx
    );

    let iqcp_dxy = result_1.get(1, 1, 1, 1);
    let pyscf_dxy = 0.5397448090050628;
    println!(
        "Ratio IQCP/PySCF for (dxy,dxy|dxy,dxy): {}",
        iqcp_dxy / pyscf_dxy
    );

    // Check normalization factor for d orbitals
    use std::f64::consts::PI;
    let alpha = 1.0;

    // dxx: powers (2,0,0), L=2
    // N = (2*alpha/pi)^(3/4) * (4*alpha)^(L/2) / sqrt((2i-1)!! * (2j-1)!! * (2k-1)!!)
    // For dxx: i=2, j=0, k=0
    // (2*2-1)!! = 3!! = 3, (2*0-1)!! = (-1)!! = 1, (2*0-1)!! = 1
    let base = (2.0 * alpha / PI).powf(0.75);
    let ang = (4.0 * alpha).powf(1.0); // L/2 = 1
    let denom_dxx = 3.0_f64.sqrt(); // sqrt(3 * 1 * 1)
    let norm_dxx = base * ang / denom_dxx;
    println!("\nNormalization constants:");
    println!("  Base: {}", base);
    println!("  Angular: {}", ang);
    println!("  Denom for dxx: {}", denom_dxx);
    println!("  Full norm for dxx: {}", norm_dxx);
    println!("  norm_dxx^4 = {}", norm_dxx.powi(4));

    // For dxy: i=1, j=1, k=0
    // (2*1-1)!! = 1, (2*1-1)!! = 1, (-1)!! = 1
    let denom_dxy = 1.0; // sqrt(1 * 1 * 1)
    let norm_dxy = base * ang / denom_dxy;
    println!("  Full norm for dxy: {}", norm_dxy);
    println!("  norm_dxy^4 = {}", norm_dxy.powi(4));

    // The raw integral before normalization
    // For spherically symmetric (dddd) at same center, the raw integral
    // should be related to Boys function F_k(0)
    println!("\nRaw integral (before normalization):");
    println!("  dxx/norm^4 = {}", iqcp_dxx / norm_dxx.powi(4));
    println!("  dxy/norm^4 = {}", iqcp_dxy / norm_dxy.powi(4));
}
