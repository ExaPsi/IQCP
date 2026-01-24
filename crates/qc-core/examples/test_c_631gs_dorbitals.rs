// Test d-orbital integrals for Carbon 6-31G* basis
use qc_core::basis::{Atom, BasisSet};
use qc_core::integrals::{hcore_matrix, kinetic_matrix, nuclear_matrix, overlap_matrix};

fn main() {
    // Single Carbon atom at origin with 6-31G* basis
    let c = Atom::new(6, [0.0, 0.0, 0.0]).unwrap();
    let atoms = vec![c];

    // Build 6-31G* basis (includes d polarization)
    let basis = BasisSet::build(atoms, "6-31g*").expect("Failed to build basis");

    println!("=== Carbon 6-31G* Basis Info ===");
    println!("Number of basis functions: {}", basis.n_basis);
    println!("Number of electrons: {}", basis.n_electrons);
    println!("Number of shells: {}", basis.n_shells());

    // Print shell information
    println!("\nShell info:");
    for (i, shell) in basis.shells.iter().enumerate() {
        println!(
            "  Shell {}: L={}, center={:?}, n_prims={}",
            i,
            shell.l_value(),
            shell.center,
            shell.primitives.len()
        );
    }

    // Compute matrices
    let s_matrix = overlap_matrix(&basis);
    let t_matrix = kinetic_matrix(&basis);
    let v_matrix = nuclear_matrix(&basis);
    let h_core = hcore_matrix(&basis);

    let nbf = basis.n_basis;

    // Print diagonal elements to compare with PySCF
    println!("\n=== Diagonal Elements Comparison ===");
    println!("                                 S           T           V         H_core");
    for i in 0..nbf {
        let label = match i {
            0 => "1s    ",
            1 => "2s    ",
            2 => "3s    ",
            3 => "2px   ",
            4 => "2py   ",
            5 => "2pz   ",
            6 => "3px   ",
            7 => "3py   ",
            8 => "3pz   ",
            9 => "3dxx  ",
            10 => "3dxy  ",
            11 => "3dxz  ",
            12 => "3dyy  ",
            13 => "3dyz  ",
            14 => "3dzz  ",
            _ => "??    ",
        };
        println!(
            "  {} ({}): {:12.6} {:12.6} {:12.6} {:12.6}",
            label,
            i,
            s_matrix[i * nbf + i],
            t_matrix[i * nbf + i],
            v_matrix[i * nbf + i],
            h_core[i * nbf + i]
        );
    }

    // Print traces
    let s_trace: f64 = (0..nbf).map(|i| s_matrix[i * nbf + i]).sum();
    let t_trace: f64 = (0..nbf).map(|i| t_matrix[i * nbf + i]).sum();
    let v_trace: f64 = (0..nbf).map(|i| v_matrix[i * nbf + i]).sum();
    let h_trace: f64 = (0..nbf).map(|i| h_core[i * nbf + i]).sum();

    println!("\n=== Traces ===");
    println!("  S trace:      {:.8}", s_trace);
    println!("  T trace:      {:.8}", t_trace);
    println!("  V trace:      {:.8}", v_trace);
    println!("  H_core trace: {:.8}", h_trace);

    // PySCF reference values
    println!("\n=== PySCF Reference ===");
    println!("  T trace (PySCF): 45.29942880");
    println!("  V trace (PySCF): -116.32314761");
    println!("  H trace (PySCF): -71.02371881");

    // D-orbital specific comparison
    println!("\n=== D-orbital Diagonals (PySCF Reference) ===");
    println!("  dxx V: -11.479022 (PySCF)");
    println!("  dxy V: -3.826341 (PySCF)");
    println!("  dxz V: -3.826341 (PySCF)");
    println!("  dyy V: -11.479022 (PySCF)");
    println!("  dyz V: -3.826341 (PySCF)");
    println!("  dzz V: -11.479022 (PySCF)");

    // Check if all d-orbitals have the same value (which would be wrong)
    let d_indices = [9, 10, 11, 12, 13, 14];
    let d_v_values: Vec<f64> = d_indices.iter().map(|&i| v_matrix[i * nbf + i]).collect();

    println!("\n=== IQCP D-orbital V diagonals ===");
    for (idx, &i) in d_indices.iter().enumerate() {
        let label = match idx {
            0 => "dxx",
            1 => "dxy",
            2 => "dxz",
            3 => "dyy",
            4 => "dyz",
            5 => "dzz",
            _ => "???",
        };
        println!("  {} V[{},{}]: {:.6}", label, i, i, d_v_values[idx]);
    }

    // Check if values are incorrectly identical
    let all_same = d_v_values.windows(2).all(|w| (w[0] - w[1]).abs() < 1e-10);
    if all_same {
        println!("\n  *** BUG DETECTED: All d-orbital V diagonals are identical! ***");
        println!("  This is WRONG - 'squared' types (dxx,dyy,dzz) should differ from 'product' types (dxy,dxz,dyz)");
    }

    // Also check T
    let d_t_values: Vec<f64> = d_indices.iter().map(|&i| t_matrix[i * nbf + i]).collect();
    println!("\n=== IQCP D-orbital T diagonals ===");
    for (idx, &i) in d_indices.iter().enumerate() {
        let label = match idx {
            0 => "dxx",
            1 => "dxy",
            2 => "dxz",
            3 => "dyy",
            4 => "dyz",
            5 => "dzz",
            _ => "???",
        };
        println!("  {} T[{},{}]: {:.6}", label, i, i, d_t_values[idx]);
    }

    let t_all_same = d_t_values.windows(2).all(|w| (w[0] - w[1]).abs() < 1e-10);
    if t_all_same {
        println!("\n  *** BUG DETECTED: All d-orbital T diagonals are identical! ***");
        println!("  This is WRONG - 'squared' types should differ from 'product' types");
    }
}
