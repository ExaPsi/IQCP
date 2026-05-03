//! Validate IQCP against PySCF for cc-pVDZ basis set.
//!
//! Run: cargo run --release --example validate_ccpvdz
//!
//! Computes RHF, LDA, and B3LYP energies for 6 molecules x 1 basis set = 18 systems
//! and compares against PySCF 2.11.0 reference values (Cartesian d-functions).

use qc_core::basis::{Atom, BasisSet, ANGSTROM_TO_BOHR};
use qc_core::dft::ks_scf::ks_scf;
use qc_core::dft::{build_becke_grid, B3lyp, GridConfig, GridQuality, Lda};
use qc_core::integrals;
use qc_core::scf::{rhf_scf, PresetSystem, ScfConfig};

// ============================================================================
// Molecule builders (coordinates in Bohr, matching validation.json)
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
    // C6H6 coordinates in Angstrom, converted to Bohr
    let coords_ang: Vec<(u8, [f64; 3])> = vec![
        (6, [0.0, 1.397, 0.0]),
        (6, [1.2098, 0.6985, 0.0]),
        (6, [1.2098, -0.6985, 0.0]),
        (6, [0.0, -1.397, 0.0]),
        (6, [-1.2098, -0.6985, 0.0]),
        (6, [-1.2098, 0.6985, 0.0]),
        (1, [0.0, 2.481, 0.0]),
        (1, [2.1486, 1.2405, 0.0]),
        (1, [2.1486, -1.2405, 0.0]),
        (1, [0.0, -2.481, 0.0]),
        (1, [-2.1486, -1.2405, 0.0]),
        (1, [-2.1486, 1.2405, 0.0]),
    ];

    coords_ang
        .into_iter()
        .map(|(z, [x, y, z_coord])| {
            Atom::new(
                z,
                [
                    x * ANGSTROM_TO_BOHR,
                    y * ANGSTROM_TO_BOHR,
                    z_coord * ANGSTROM_TO_BOHR,
                ],
            )
            .unwrap()
        })
        .collect()
}

// ============================================================================
// System builder
// ============================================================================

fn build_system(atoms: &[Atom], basis_name: &str) -> (PresetSystem, BasisSet) {
    let basis = BasisSet::build(atoms.to_vec(), basis_name).unwrap();
    let nbf = basis.n_basis;
    let nelec = basis.n_electrons;

    let s_matrix = integrals::overlap_matrix(&basis);
    let h_core = integrals::hcore_matrix(&basis);
    let eri = integrals::eri_compressed(&basis);

    let system = PresetSystem {
        system_id: "validate".to_string(),
        label: "Validation system".to_string(),
        nbf,
        nelec,
        e_nuc: basis.nuclear_repulsion,
        s_matrix,
        h_core,
        eri_compressed: eri,
    };

    (system, basis)
}

// ============================================================================
// Computation runners
// ============================================================================

fn run_rhf(atoms: &[Atom], basis_name: &str) -> (f64, bool) {
    let (system, _basis) = build_system(atoms, basis_name);
    let config = ScfConfig {
        use_diis: true,
        max_iterations: 200,
        ..ScfConfig::tight()
    };
    match rhf_scf(&system, &config) {
        Ok(result) => (result.energy_total, result.converged),
        Err(e) => {
            eprintln!("    RHF error: {}", e);
            (f64::NAN, false)
        }
    }
}

fn run_lda(atoms: &[Atom], basis_name: &str) -> (f64, bool) {
    let (system, basis) = build_system(atoms, basis_name);
    let grid_config = GridConfig::default();
    let grid = build_becke_grid(&basis.atoms, &grid_config);
    let lda = Lda::new();
    let config = ScfConfig {
        use_diis: true,
        max_iterations: 200,
        ..ScfConfig::tight()
    };
    match ks_scf(&system, &config, &lda, &grid, &basis, false, None) {
        Ok(result) => (result.scf_output.energy_total, result.scf_output.converged),
        Err(e) => {
            eprintln!("    LDA error: {}", e);
            (f64::NAN, false)
        }
    }
}

fn run_b3lyp(atoms: &[Atom], basis_name: &str) -> (f64, bool) {
    let (system, basis) = build_system(atoms, basis_name);
    let grid_config = GridConfig {
        n_radial: 99,
        quality: GridQuality::Fine,
        pruning: false,
    };
    let grid = build_becke_grid(&basis.atoms, &grid_config);
    let b3lyp = B3lyp::new();
    let config = ScfConfig {
        use_diis: true,
        max_iterations: 200,
        ..ScfConfig::tight()
    };
    match ks_scf(&system, &config, &b3lyp, &grid, &basis, false, None) {
        Ok(result) => (result.scf_output.energy_total, result.scf_output.converged),
        Err(e) => {
            eprintln!("    B3LYP error: {}", e);
            (f64::NAN, false)
        }
    }
}

// ============================================================================
// PySCF reference values (PySCF 2.11.0, conv_tol=1e-10, mol.cart=True)
// ============================================================================

struct RefEntry {
    molecule: &'static str,
    basis: &'static str,
    method: &'static str,
    energy: f64,
    n_basis: usize,
}

fn pyscf_references() -> Vec<RefEntry> {
    vec![
        // H2: 10 Cartesian basis functions
        RefEntry {
            molecule: "H2",
            basis: "cc-pvdz",
            method: "RHF",
            energy: -1.128709448979890e+00,
            n_basis: 10,
        },
        RefEntry {
            molecule: "H2",
            basis: "cc-pvdz",
            method: "LDA",
            energy: -1.131458096785010e+00,
            n_basis: 10,
        },
        RefEntry {
            molecule: "H2",
            basis: "cc-pvdz",
            method: "B3LYP",
            energy: -1.166534775939096e+00,
            n_basis: 10,
        },
        // HF: 20 Cartesian basis functions
        RefEntry {
            molecule: "HF",
            basis: "cc-pvdz",
            method: "RHF",
            energy: -1.000198311869585e+02,
            n_basis: 20,
        },
        RefEntry {
            molecule: "HF",
            basis: "cc-pvdz",
            method: "LDA",
            energy: -9.978413784925793e+01,
            n_basis: 20,
        },
        RefEntry {
            molecule: "HF",
            basis: "cc-pvdz",
            method: "B3LYP",
            energy: -1.003990982893055e+02,
            n_basis: 20,
        },
        // H2O: 25 Cartesian basis functions
        RefEntry {
            molecule: "H2O",
            basis: "cc-pvdz",
            method: "RHF",
            energy: -7.602711576564408e+01,
            n_basis: 25,
        },
        RefEntry {
            molecule: "H2O",
            basis: "cc-pvdz",
            method: "LDA",
            energy: -7.585751315191276e+01,
            n_basis: 25,
        },
        RefEntry {
            molecule: "H2O",
            basis: "cc-pvdz",
            method: "B3LYP",
            energy: -7.638441293702004e+01,
            n_basis: 25,
        },
        // NH3: 30 Cartesian basis functions
        RefEntry {
            molecule: "NH3",
            basis: "cc-pvdz",
            method: "RHF",
            energy: -5.619570893849166e+01,
            n_basis: 30,
        },
        RefEntry {
            molecule: "NH3",
            basis: "cc-pvdz",
            method: "LDA",
            energy: -5.606876302616544e+01,
            n_basis: 30,
        },
        RefEntry {
            molecule: "NH3",
            basis: "cc-pvdz",
            method: "B3LYP",
            energy: -5.651825018412060e+01,
            n_basis: 30,
        },
        // CH4: 35 Cartesian basis functions
        RefEntry {
            molecule: "CH4",
            basis: "cc-pvdz",
            method: "RHF",
            energy: -4.019872492896991e+01,
            n_basis: 35,
        },
        RefEntry {
            molecule: "CH4",
            basis: "cc-pvdz",
            method: "LDA",
            energy: -4.009614548862446e+01,
            n_basis: 35,
        },
        RefEntry {
            molecule: "CH4",
            basis: "cc-pvdz",
            method: "B3LYP",
            energy: -4.048069420110117e+01,
            n_basis: 35,
        },
        // C6H6: 120 Cartesian basis functions
        RefEntry {
            molecule: "C6H6",
            basis: "cc-pvdz",
            method: "RHF",
            energy: -2.307226386431662e+02,
            n_basis: 120,
        },
        RefEntry {
            molecule: "C6H6",
            basis: "cc-pvdz",
            method: "LDA",
            energy: -2.301079958936073e+02,
            n_basis: 120,
        },
        RefEntry {
            molecule: "C6H6",
            basis: "cc-pvdz",
            method: "B3LYP",
            energy: -2.321128375305050e+02,
            n_basis: 120,
        },
    ]
}

// ============================================================================
// Main
// ============================================================================

fn main() {
    println!("==========================================================================");
    println!("  IQCP cc-pVDZ Validation vs PySCF 2.11.0");
    println!("  Basis set: cc-pVDZ (Cartesian d-functions, 6d)");
    println!("  Methods: RHF, LDA (Slater+VWN5), B3LYP (b3lyp5)");
    println!("==========================================================================\n");

    let refs = pyscf_references();

    // Collect results
    struct Result {
        molecule: String,
        basis: String,
        method: String,
        pyscf_energy: f64,
        iqcp_energy: f64,
        delta_e: f64,
        converged: bool,
        pass: bool,
        n_basis_expected: usize,
        n_basis_actual: usize,
    }

    let mut results: Vec<Result> = Vec::new();

    let molecules: Vec<(&str, Vec<Atom>)> = vec![
        ("H2", h2_atoms()),
        ("HF", hf_atoms()),
        ("H2O", h2o_atoms()),
        ("NH3", nh3_atoms()),
        ("CH4", ch4_atoms()),
        ("C6H6", c6h6_atoms()),
    ];

    let methods = ["RHF", "LDA", "B3LYP"];

    for (mol_name, atoms) in &molecules {
        // Check basis function count
        let basis = BasisSet::build(atoms.clone(), "cc-pvdz").unwrap();
        let nbf_actual = basis.n_basis;

        for method in &methods {
            // Find PySCF reference
            let pyscf_ref = refs
                .iter()
                .find(|r| r.molecule == *mol_name && r.method == *method)
                .expect("Missing PySCF reference");

            print!(
                "  {:<5} / cc-pVDZ / {:<5} (N_bf={:>3}) ... ",
                mol_name, method, nbf_actual
            );

            let (iqcp_energy, converged) = match *method {
                "RHF" => run_rhf(atoms, "cc-pvdz"),
                "LDA" => run_lda(atoms, "cc-pvdz"),
                "B3LYP" => run_b3lyp(atoms, "cc-pvdz"),
                _ => unreachable!(),
            };

            let delta_e = (iqcp_energy - pyscf_ref.energy).abs();

            // Tolerance: 100 nHa for RHF, 200 uHa for DFT (grid differences)
            let tolerance = if *method == "RHF" { 1e-7 } else { 2e-4 };
            let pass = converged && delta_e < tolerance;

            let status = if !converged {
                "FAIL (not converged)"
            } else if pass {
                "PASS"
            } else {
                "FAIL"
            };

            println!("dE = {:.2e}  [{}]", delta_e, status);

            results.push(Result {
                molecule: mol_name.to_string(),
                basis: "cc-pvdz".to_string(),
                method: method.to_string(),
                pyscf_energy: pyscf_ref.energy,
                iqcp_energy,
                delta_e,
                converged,
                pass,
                n_basis_expected: pyscf_ref.n_basis,
                n_basis_actual: nbf_actual,
            });
        }
    }

    // Summary
    println!("\n==========================================================================");
    println!("  SUMMARY");
    println!("==========================================================================\n");

    let total = results.len();
    let passed = results.iter().filter(|r| r.pass).count();
    let failed = total - passed;

    println!("Total: {}  Passed: {}  Failed: {}\n", total, passed, failed);

    // Check basis function counts
    println!("Basis function counts (IQCP vs PySCF):");
    for r in &results {
        if r.method == "RHF" {
            let bf_match = if r.n_basis_actual == r.n_basis_expected {
                "OK"
            } else {
                "MISMATCH"
            };
            println!(
                "  {:<5}: IQCP={}, PySCF={} [{}]",
                r.molecule, r.n_basis_actual, r.n_basis_expected, bf_match
            );
        }
    }
    println!();

    if failed > 0 {
        println!("FAILURES:");
        for r in &results {
            if !r.pass {
                println!(
                    "  {} / {} / {}: PySCF={:.12}, IQCP={:.12}, dE={:.2e}, converged={}",
                    r.molecule,
                    r.basis,
                    r.method,
                    r.pyscf_energy,
                    r.iqcp_energy,
                    r.delta_e,
                    r.converged
                );
            }
        }
        println!();
    }

    // Print deviation table for documentation
    println!("\n==========================================================================");
    println!("  DEVIATION TABLE (nHa for RHF, uHa for DFT)");
    println!("==========================================================================\n");
    println!(
        "{:<5} {:<8} {:<5} {:>5} {:>20} {:>20} {:>12} {:>6}",
        "Mol", "Basis", "Meth", "N_bf", "PySCF (Ha)", "IQCP (Ha)", "|dE|", "Pass?"
    );
    println!("{}", "-".repeat(85));

    for r in &results {
        let pass_str = if r.pass { "OK" } else { "FAIL" };
        let (value, unit) = if r.method == "RHF" {
            (r.delta_e * 1e9, "nHa")
        } else {
            (r.delta_e * 1e6, "uHa")
        };
        println!(
            "{:<5} {:<8} {:<5} {:>5} {:>20.12} {:>20.12} {:>8.2} {:<3} {:>6}",
            r.molecule,
            r.basis,
            r.method,
            r.n_basis_actual,
            r.pyscf_energy,
            r.iqcp_energy,
            value,
            unit,
            pass_str
        );
    }

    if failed > 0 {
        std::process::exit(1);
    }
}
