//! Validate IQCP against PySCF for 3-21G, 6-31G, and 6-31+G* basis sets.
//!
//! Run: cargo run --release --example validate_additional_bases
//!
//! Computes RHF, LDA, and B3LYP energies for 6 molecules x 3 basis sets = 54 systems
//! and compares against PySCF 2.11.0 reference values.

use qc_core::basis::{Atom, BasisSet, ANGSTROM_TO_BOHR};
use qc_core::dft::ks_scf::ks_scf;
use qc_core::dft::{build_becke_grid, B3lyp, GridConfig, GridQuality, Lda};
use qc_core::integrals;
use qc_core::scf::{rhf_scf, PresetSystem, ScfConfig};

// ============================================================================
// Molecule builders (coordinates in Bohr, matching tests/golden/dft/validation.json)
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
        label: format!("Validation system"),
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
// PySCF reference values (from generate_additional_basis_refs.py, PySCF 2.11.0)
// ============================================================================

struct RefEntry {
    molecule: &'static str,
    basis: &'static str,
    method: &'static str,
    energy: f64,
}

fn pyscf_references() -> Vec<RefEntry> {
    vec![
        // PySCF 2.11.0, conv_tol=1e-12, single-threaded
        // ===== H2 =====
        RefEntry {
            molecule: "H2",
            basis: "3-21g",
            method: "RHF",
            energy: -1.1229333636171903,
        },
        RefEntry {
            molecule: "H2",
            basis: "3-21g",
            method: "LDA",
            energy: -1.1281068730763897,
        },
        RefEntry {
            molecule: "H2",
            basis: "3-21g",
            method: "B3LYP",
            energy: -1.163816800691356,
        },
        RefEntry {
            molecule: "H2",
            basis: "6-31g",
            method: "RHF",
            energy: -1.1267427044518272,
        },
        RefEntry {
            molecule: "H2",
            basis: "6-31g",
            method: "LDA",
            energy: -1.1326691268172882,
        },
        RefEntry {
            molecule: "H2",
            basis: "6-31g",
            method: "B3LYP",
            energy: -1.1687129897605577,
        },
        RefEntry {
            molecule: "H2",
            basis: "6-31+g*",
            method: "RHF",
            energy: -1.1267427044518272,
        },
        RefEntry {
            molecule: "H2",
            basis: "6-31+g*",
            method: "LDA",
            energy: -1.1326691268172882,
        },
        RefEntry {
            molecule: "H2",
            basis: "6-31+g*",
            method: "B3LYP",
            energy: -1.1687129897605577,
        },
        // ===== HF =====
        RefEntry {
            molecule: "HF",
            basis: "3-21g",
            method: "RHF",
            energy: -99.45974981319155,
        },
        RefEntry {
            molecule: "HF",
            basis: "3-21g",
            method: "LDA",
            energy: -99.20896816812308,
        },
        RefEntry {
            molecule: "HF",
            basis: "3-21g",
            method: "B3LYP",
            energy: -99.82260737485956,
        },
        RefEntry {
            molecule: "HF",
            basis: "6-31g",
            method: "RHF",
            energy: -99.98340857072898,
        },
        RefEntry {
            molecule: "HF",
            basis: "6-31g",
            method: "LDA",
            energy: -99.74715390420434,
        },
        RefEntry {
            molecule: "HF",
            basis: "6-31g",
            method: "B3LYP",
            energy: -100.36475793613916,
        },
        RefEntry {
            molecule: "HF",
            basis: "6-31+g*",
            method: "RHF",
            energy: -100.01485609801351,
        },
        RefEntry {
            molecule: "HF",
            basis: "6-31+g*",
            method: "LDA",
            energy: -99.79217408401911,
        },
        RefEntry {
            molecule: "HF",
            basis: "6-31+g*",
            method: "B3LYP",
            energy: -100.40526778906425,
        },
        // ===== H2O =====
        RefEntry {
            molecule: "H2O",
            basis: "3-21g",
            method: "RHF",
            energy: -75.58540012457607,
        },
        RefEntry {
            molecule: "H2O",
            basis: "3-21g",
            method: "LDA",
            energy: -75.40742660806734,
        },
        RefEntry {
            molecule: "H2O",
            basis: "3-21g",
            method: "B3LYP",
            energy: -75.93431793211136,
        },
        RefEntry {
            molecule: "H2O",
            basis: "6-31g",
            method: "RHF",
            energy: -75.98396375334872,
        },
        RefEntry {
            molecule: "H2O",
            basis: "6-31g",
            method: "LDA",
            energy: -75.8179168183294,
        },
        RefEntry {
            molecule: "H2O",
            basis: "6-31g",
            method: "B3LYP",
            energy: -76.34780377303699,
        },
        RefEntry {
            molecule: "H2O",
            basis: "6-31+g*",
            method: "RHF",
            energy: -76.01744020557382,
        },
        RefEntry {
            molecule: "H2O",
            basis: "6-31+g*",
            method: "LDA",
            energy: -75.85968105803343,
        },
        RefEntry {
            molecule: "H2O",
            basis: "6-31+g*",
            method: "B3LYP",
            energy: -76.3852346530983,
        },
        // ===== NH3 =====
        RefEntry {
            molecule: "NH3",
            basis: "3-21g",
            method: "RHF",
            energy: -55.87022574746423,
        },
        RefEntry {
            molecule: "NH3",
            basis: "3-21g",
            method: "LDA",
            energy: -55.74189049833342,
        },
        RefEntry {
            molecule: "NH3",
            basis: "3-21g",
            method: "B3LYP",
            energy: -56.191338633184614,
        },
        RefEntry {
            molecule: "NH3",
            basis: "6-31g",
            method: "RHF",
            energy: -56.16061034823184,
        },
        RefEntry {
            molecule: "NH3",
            basis: "6-31g",
            method: "LDA",
            energy: -56.039965930025446,
        },
        RefEntry {
            molecule: "NH3",
            basis: "6-31g",
            method: "B3LYP",
            energy: -56.49206630151258,
        },
        RefEntry {
            molecule: "NH3",
            basis: "6-31+g*",
            method: "RHF",
            energy: -56.188994544888104,
        },
        RefEntry {
            molecule: "NH3",
            basis: "6-31+g*",
            method: "LDA",
            energy: -56.07153688376367,
        },
        RefEntry {
            molecule: "NH3",
            basis: "6-31+g*",
            method: "B3LYP",
            energy: -56.5203443531623,
        },
        // ===== CH4 =====
        RefEntry {
            molecule: "CH4",
            basis: "3-21g",
            method: "RHF",
            energy: -39.976847325231226,
        },
        RefEntry {
            molecule: "CH4",
            basis: "3-21g",
            method: "LDA",
            energy: -39.88279709708183,
        },
        RefEntry {
            molecule: "CH4",
            basis: "3-21g",
            method: "B3LYP",
            energy: -40.26547078101928,
        },
        RefEntry {
            molecule: "CH4",
            basis: "6-31g",
            method: "RHF",
            energy: -40.180507904,
        },
        RefEntry {
            molecule: "CH4",
            basis: "6-31g",
            method: "LDA",
            energy: -40.08939576149906,
        },
        RefEntry {
            molecule: "CH4",
            basis: "6-31g",
            method: "B3LYP",
            energy: -40.47452571835594,
        },
        RefEntry {
            molecule: "CH4",
            basis: "6-31+g*",
            method: "RHF",
            energy: -40.19565408921565,
        },
        RefEntry {
            molecule: "CH4",
            basis: "6-31+g*",
            method: "LDA",
            energy: -40.09971326050533,
        },
        RefEntry {
            molecule: "CH4",
            basis: "6-31+g*",
            method: "B3LYP",
            energy: -40.484487022022414,
        },
        // ===== C6H6 =====
        RefEntry {
            molecule: "C6H6",
            basis: "3-21g",
            method: "RHF",
            energy: -229.41789108926517,
        },
        RefEntry {
            molecule: "C6H6",
            basis: "3-21g",
            method: "LDA",
            energy: -228.83079659739045,
        },
        RefEntry {
            molecule: "C6H6",
            basis: "3-21g",
            method: "B3LYP",
            energy: -230.82146740765066,
        },
        RefEntry {
            molecule: "C6H6",
            basis: "6-31g",
            method: "RHF",
            energy: -230.623508002016,
        },
        RefEntry {
            molecule: "C6H6",
            basis: "6-31g",
            method: "LDA",
            energy: -230.03682236143194,
        },
        RefEntry {
            molecule: "C6H6",
            basis: "6-31g",
            method: "B3LYP",
            energy: -232.0442874355734,
        },
        RefEntry {
            molecule: "C6H6",
            basis: "6-31+g*",
            method: "RHF",
            energy: -230.71030489827106,
        },
        RefEntry {
            molecule: "C6H6",
            basis: "6-31+g*",
            method: "LDA",
            energy: -230.10009752381293,
        },
        RefEntry {
            molecule: "C6H6",
            basis: "6-31+g*",
            method: "B3LYP",
            energy: -232.1046321682426,
        },
    ]
}

// ============================================================================
// Main
// ============================================================================

fn main() {
    println!("==========================================================================");
    println!("  IQCP Additional Basis Set Validation vs PySCF 2.11.0");
    println!("  Basis sets: 3-21G, 6-31G, 6-31+G*");
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

    let basis_sets = ["3-21g", "6-31g", "6-31+g*"];
    let methods = ["RHF", "LDA", "B3LYP"];

    for (mol_name, atoms) in &molecules {
        for basis_name in &basis_sets {
            for method in &methods {
                // Find PySCF reference
                let pyscf_ref = refs
                    .iter()
                    .find(|r| {
                        r.molecule == *mol_name && r.basis == *basis_name && r.method == *method
                    })
                    .expect("Missing PySCF reference");

                print!("  {:<5} / {:<8} / {:<5} ... ", mol_name, basis_name, method);

                let (iqcp_energy, converged) = match *method {
                    "RHF" => run_rhf(atoms, basis_name),
                    "LDA" => run_lda(atoms, basis_name),
                    "B3LYP" => run_b3lyp(atoms, basis_name),
                    _ => unreachable!(),
                };

                let delta_e = (iqcp_energy - pyscf_ref.energy).abs();

                // Tolerance: 1e-5 for RHF, 1e-4 for DFT (grid differences)
                let tolerance = if *method == "RHF" { 1e-5 } else { 1e-4 };
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
                    basis: basis_name.to_string(),
                    method: method.to_string(),
                    pyscf_energy: pyscf_ref.energy,
                    iqcp_energy,
                    delta_e,
                    converged,
                    pass,
                });
            }
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

    // Print full table for documentation
    println!("\n==========================================================================");
    println!("  FULL RESULTS TABLE");
    println!("==========================================================================\n");
    println!(
        "{:<5} {:<8} {:<5} {:>20} {:>20} {:>12} {:>6}",
        "Mol", "Basis", "Meth", "PySCF (Ha)", "IQCP (Ha)", "|dE| (Ha)", "Pass?"
    );
    println!("{}", "-".repeat(82));

    for r in &results {
        let pass_str = if r.pass { "OK" } else { "FAIL" };
        println!(
            "{:<5} {:<8} {:<5} {:>20.12} {:>20.12} {:>12.2e} {:>6}",
            r.molecule, r.basis, r.method, r.pyscf_energy, r.iqcp_energy, r.delta_e, pass_str
        );
    }

    // Print nHa/uHa table for documentation
    println!("\n==========================================================================");
    println!("  DEVIATION TABLE (nHa for RHF, uHa for DFT)");
    println!("==========================================================================\n");
    println!(
        "{:<5} {:<8} {:<5} {:>12} {:>6}",
        "Mol", "Basis", "Meth", "|dE|", "Pass?"
    );
    println!("{}", "-".repeat(42));

    for r in &results {
        let pass_str = if r.pass { "OK" } else { "FAIL" };
        let (value, unit) = if r.method == "RHF" {
            (r.delta_e * 1e9, "nHa")
        } else {
            (r.delta_e * 1e6, "uHa")
        };
        println!(
            "{:<5} {:<8} {:<5} {:>8.2} {:<3} {:>6}",
            r.molecule, r.basis, r.method, value, unit, pass_str
        );
    }

    if failed > 0 {
        std::process::exit(1);
    }
}
