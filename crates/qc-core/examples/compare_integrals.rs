//! Compare IQCP integrals against PySCF reference for HF/6-31G* and NH3/STO-3G.
//!
//! Run: cargo run --release --example compare_integrals

use qc_core::basis::{Atom, BasisSet};
use qc_core::integrals;
use serde::Deserialize;

#[derive(Deserialize)]
struct IntegralRef {
    n_basis: usize,
    nuclear_repulsion: f64,
    overlap: Vec<f64>,
    hcore: Vec<f64>,
    eri_compressed: Vec<f64>,
}

fn compare_system(name: &str, atoms: &[Atom], basis_name: &str, ref_path: &str) {
    println!("\n{}", "=".repeat(70));
    println!("  {} / {}", name, basis_name);
    println!("{}", "=".repeat(70));

    // Build IQCP integrals
    let basis = BasisSet::build(atoms.to_vec(), basis_name).unwrap();
    let n = basis.n_basis;
    println!("IQCP basis functions: {}", n);

    let s_iqcp = integrals::overlap_matrix(&basis);
    let h_iqcp = integrals::hcore_matrix(&basis);
    let eri_iqcp = integrals::eri_compressed(&basis);

    // Load PySCF reference
    let ref_json = std::fs::read_to_string(ref_path).unwrap();
    let pyscf: IntegralRef = serde_json::from_str(&ref_json).unwrap();
    assert_eq!(n, pyscf.n_basis, "Basis size mismatch!");

    // Compare nuclear repulsion
    let enuc_diff = (basis.nuclear_repulsion - pyscf.nuclear_repulsion).abs();
    println!("\nNuclear repulsion:");
    println!("  IQCP:  {:.15e}", basis.nuclear_repulsion);
    println!("  PySCF: {:.15e}", pyscf.nuclear_repulsion);
    println!("  Diff:  {:.2e}", enuc_diff);

    // Compare overlap matrix
    let mut max_s_err = 0.0_f64;
    let mut max_s_idx = (0, 0);
    for i in 0..n {
        for j in 0..n {
            let iqcp_val = s_iqcp[i * n + j];
            let pyscf_val = pyscf.overlap[i * n + j];
            let err = (iqcp_val - pyscf_val).abs();
            if err > max_s_err {
                max_s_err = err;
                max_s_idx = (i, j);
            }
        }
    }
    println!("\nOverlap matrix (S):");
    println!(
        "  Max |error|: {:.2e} at ({}, {})",
        max_s_err, max_s_idx.0, max_s_idx.1
    );
    if max_s_err > 1e-12 {
        let i = max_s_idx.0;
        let j = max_s_idx.1;
        println!("  IQCP:  {:.15e}", s_iqcp[i * n + j]);
        println!("  PySCF: {:.15e}", pyscf.overlap[i * n + j]);
    }

    // Print all S differences > 1e-13
    println!("\n  S elements with |error| > 1e-13:");
    for i in 0..n {
        for j in i..n {
            let iqcp_val = s_iqcp[i * n + j];
            let pyscf_val = pyscf.overlap[i * n + j];
            let err = (iqcp_val - pyscf_val).abs();
            if err > 1e-13 {
                println!(
                    "    S[{},{}]: IQCP={:+.15e}, PySCF={:+.15e}, diff={:+.2e}",
                    i,
                    j,
                    iqcp_val,
                    pyscf_val,
                    iqcp_val - pyscf_val
                );
            }
        }
    }

    // Compare core Hamiltonian
    let mut max_h_err = 0.0_f64;
    let mut max_h_idx = (0, 0);
    for i in 0..n {
        for j in 0..n {
            let iqcp_val = h_iqcp[i * n + j];
            let pyscf_val = pyscf.hcore[i * n + j];
            let err = (iqcp_val - pyscf_val).abs();
            if err > max_h_err {
                max_h_err = err;
                max_h_idx = (i, j);
            }
        }
    }
    println!("\nCore Hamiltonian (H):");
    println!(
        "  Max |error|: {:.2e} at ({}, {})",
        max_h_err, max_h_idx.0, max_h_idx.1
    );
    if max_h_err > 1e-12 {
        let i = max_h_idx.0;
        let j = max_h_idx.1;
        println!("  IQCP:  {:.15e}", h_iqcp[i * n + j]);
        println!("  PySCF: {:.15e}", pyscf.hcore[i * n + j]);
    }

    // Print all H differences > 1e-10
    println!("\n  H elements with |error| > 1e-10:");
    for i in 0..n {
        for j in i..n {
            let iqcp_val = h_iqcp[i * n + j];
            let pyscf_val = pyscf.hcore[i * n + j];
            let err = (iqcp_val - pyscf_val).abs();
            if err > 1e-10 {
                println!(
                    "    H[{},{}]: IQCP={:+.15e}, PySCF={:+.15e}, diff={:+.2e}",
                    i,
                    j,
                    iqcp_val,
                    pyscf_val,
                    iqcp_val - pyscf_val
                );
            }
        }
    }

    // Compare ERI
    let eri_len = eri_iqcp.len();
    assert_eq!(eri_len, pyscf.eri_compressed.len(), "ERI size mismatch!");
    let mut max_eri_err = 0.0_f64;
    let mut max_eri_idx = 0;
    let mut eri_errors_above_1e10 = 0;
    let mut eri_errors_above_1e8 = 0;
    let mut sum_sq_err = 0.0_f64;

    for k in 0..eri_len {
        let iqcp_val = eri_iqcp[k];
        let pyscf_val = pyscf.eri_compressed[k];
        let err = (iqcp_val - pyscf_val).abs();
        sum_sq_err += err * err;
        if err > 1e-10 {
            eri_errors_above_1e10 += 1;
        }
        if err > 1e-8 {
            eri_errors_above_1e8 += 1;
        }
        if err > max_eri_err {
            max_eri_err = err;
            max_eri_idx = k;
        }
    }
    let rms_eri_err = (sum_sq_err / eri_len as f64).sqrt();

    println!("\nERI (compressed, {} elements):", eri_len);
    println!(
        "  Max |error|:         {:.2e} at index {}",
        max_eri_err, max_eri_idx
    );
    println!("  RMS error:           {:.2e}", rms_eri_err);
    println!("  Elements > 1e-10:    {}", eri_errors_above_1e10);
    println!("  Elements > 1e-8:     {}", eri_errors_above_1e8);
    println!("  IQCP:  {:.15e}", eri_iqcp[max_eri_idx]);
    println!("  PySCF: {:.15e}", pyscf.eri_compressed[max_eri_idx]);

    // Decode the compressed index to (i,j,k,l)
    // ij >= kl, i >= j, k >= l
    let mut idx = 0;
    'outer: for i in 0..n {
        for j in 0..=i {
            let ij = i * (i + 1) / 2 + j;
            for k in 0..n {
                for l in 0..=k {
                    let kl = k * (k + 1) / 2 + l;
                    if ij >= kl {
                        if idx == max_eri_idx {
                            println!(
                                "  Indices: ({},{} | {},{})  [compressed idx={}]",
                                i, j, k, l, max_eri_idx
                            );
                            break 'outer;
                        }
                        idx += 1;
                    }
                }
            }
        }
    }

    // Print top 20 largest ERI errors
    let mut errors: Vec<(usize, f64, f64, f64)> = (0..eri_len)
        .map(|k| {
            (
                k,
                eri_iqcp[k],
                pyscf.eri_compressed[k],
                (eri_iqcp[k] - pyscf.eri_compressed[k]).abs(),
            )
        })
        .filter(|(_, _, _, e)| *e > 1e-12)
        .collect();
    errors.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap());
    if !errors.is_empty() {
        println!("\n  Top 20 ERI errors:");
        for (k, iqcp_val, pyscf_val, err) in errors.iter().take(20) {
            // Decode index
            let mut idx = 0;
            let mut found = (0, 0, 0, 0);
            'decode: for i in 0..n {
                for j in 0..=i {
                    let ij = i * (i + 1) / 2 + j;
                    for kk in 0..n {
                        for ll in 0..=kk {
                            let kl = kk * (kk + 1) / 2 + ll;
                            if ij >= kl {
                                if idx == *k {
                                    found = (i, j, kk, ll);
                                    break 'decode;
                                }
                                idx += 1;
                            }
                        }
                    }
                }
            }
            println!(
                "    ({:2},{:2}|{:2},{:2}): IQCP={:+.12e}, PySCF={:+.12e}, diff={:+.2e}",
                found.0,
                found.1,
                found.2,
                found.3,
                iqcp_val,
                pyscf_val,
                iqcp_val - pyscf_val
            );
        }
    }
}

fn main() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let base = std::path::Path::new(manifest_dir).join("../../tests/benchmarks/pyscf/results");

    // System 1: HF/6-31G* (65 nHa error)
    let hf_atoms = vec![
        Atom::new(1, [0.0, 0.0, 0.0]).unwrap(),
        Atom::new(9, [0.0, 0.0, 1.7328]).unwrap(),
    ];
    compare_system(
        "HF",
        &hf_atoms,
        "6-31g*",
        base.join("hf_631gs_integrals.json").to_str().unwrap(),
    );

    // System 2: NH3/STO-3G (28 nHa error)
    let nh3_atoms = vec![
        Atom::new(7, [0.0, 0.0, 0.219705]).unwrap(),
        Atom::new(1, [0.0, 1.7714918, -0.512645]).unwrap(),
        Atom::new(1, [1.5342036, -0.8857459, -0.512645]).unwrap(),
        Atom::new(1, [-1.5342036, -0.8857459, -0.512645]).unwrap(),
    ];
    compare_system(
        "NH3",
        &nh3_atoms,
        "sto-3g",
        base.join("nh3_sto3g_integrals.json").to_str().unwrap(),
    );
}
