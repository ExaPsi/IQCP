//! Generate IQCP analytical gradient values and compare against PySCF references.
//!
//! Computes RHF, LDA, and B3LYP analytical gradients for 5 molecules x 2 basis sets,
//! compares against PySCF reference values from gradient_validation.json, and outputs
//! a JSON comparison table to stdout.
//!
//! Run: cargo run --release --example gradient_comparison

use qc_core::basis::{Atom, BasisSet};
use qc_core::dft::{B3lyp, GridConfig, GridQuality, Lda};
use qc_core::integrals;
use qc_core::scf::gradient::{ks_dft_gradient, rhf_gradient, GradientResult};
use qc_core::scf::{rhf_scf, PresetSystem, ScfConfig};

use std::fs;

// ============================================================================
// PySCF reference data structures
// ============================================================================

#[derive(serde::Deserialize)]
struct PyscfReferenceFile {
    results: Vec<PyscfEntry>,
}

#[derive(serde::Deserialize)]
#[allow(dead_code)]
struct PyscfEntry {
    molecule: String,
    basis: String,
    method: String,
    n_basis: usize,
    energy: f64,
    gradient: Vec<Vec<f64>>,
    max_gradient: f64,
    rms_gradient: f64,
}

// ============================================================================
// Output structures
// ============================================================================

#[derive(serde::Serialize)]
struct ComparisonEntry {
    molecule: String,
    basis: String,
    method: String,
    n_basis: usize,
    iqcp_energy: f64,
    pyscf_energy: Option<f64>,
    energy_error: Option<f64>,
    iqcp_gradient: Vec<Vec<f64>>,
    pyscf_gradient: Option<Vec<Vec<f64>>>,
    max_abs_error: Option<f64>,
    rms_abs_error: Option<f64>,
    iqcp_max_gradient: f64,
    pyscf_max_gradient: Option<f64>,
}

// ============================================================================
// Molecule definitions (coordinates in bohr)
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

// ============================================================================
// Helpers
// ============================================================================

fn build_system(atoms: &[Atom], basis_name: &str) -> (PresetSystem, BasisSet) {
    let basis = BasisSet::build(atoms.to_vec(), basis_name).unwrap();
    let nbf = basis.n_basis;
    let nelec = basis.n_electrons;
    let s_matrix = integrals::overlap_matrix(&basis);
    let h_core = integrals::hcore_matrix(&basis);
    let eri = integrals::eri_compressed(&basis);
    let system = PresetSystem {
        system_id: "grad_cmp".to_string(),
        label: "Gradient comparison".to_string(),
        nbf,
        nelec,
        e_nuc: basis.nuclear_repulsion,
        s_matrix,
        h_core,
        eri_compressed: eri,
    };
    (system, basis)
}

/// Convert gradient result to Vec<Vec<f64>> for JSON output
fn grad_to_vecs(result: &GradientResult) -> Vec<Vec<f64>> {
    result
        .gradients
        .iter()
        .map(|g| vec![g[0], g[1], g[2]])
        .collect()
}

/// Compute max absolute error and RMS error between two gradient arrays
fn gradient_errors(iqcp: &[[f64; 3]], pyscf: &[Vec<f64>]) -> (f64, f64) {
    let mut max_err = 0.0_f64;
    let mut sum_sq = 0.0;
    let mut n = 0;
    for (ig, pg) in iqcp.iter().zip(pyscf.iter()) {
        for dir in 0..3 {
            let err = (ig[dir] - pg[dir]).abs();
            max_err = max_err.max(err);
            sum_sq += err * err;
            n += 1;
        }
    }
    let rms = (sum_sq / n as f64).sqrt();
    (max_err, rms)
}

/// Find matching PySCF reference entry
fn find_pyscf_ref<'a>(
    refs: &'a [PyscfEntry],
    molecule: &str,
    basis: &str,
    method: &str,
) -> Option<&'a PyscfEntry> {
    refs.iter()
        .find(|e| e.molecule == molecule && e.basis == basis && e.method == method)
}

// ============================================================================
// Main
// ============================================================================

fn main() {
    // Load PySCF reference data (RHF + B3LYP)
    let ref_path = "tests/golden/dft/gradient_validation.json";
    let ref_data: PyscfReferenceFile = {
        let contents = fs::read_to_string(ref_path).expect("Read gradient_validation.json");
        serde_json::from_str(&contents).expect("Parse gradient_validation.json")
    };
    let mut pyscf_refs = ref_data.results;

    // Load PySCF LDA reference data
    let lda_ref_path = "tests/golden/dft/gradient_validation_lda.json";
    if let Ok(contents) = fs::read_to_string(lda_ref_path) {
        let lda_data: PyscfReferenceFile =
            serde_json::from_str(&contents).expect("Parse gradient_validation_lda.json");
        eprintln!(
            "Loaded {} LDA reference entries from gradient_validation_lda.json",
            lda_data.results.len()
        );
        pyscf_refs.extend(lda_data.results);
    }

    eprintln!("Loaded {} total PySCF reference entries", pyscf_refs.len());

    let scf_config = ScfConfig {
        use_diis: true,
        ..ScfConfig::tight()
    };
    let grid_config = GridConfig {
        n_radial: 75,
        quality: GridQuality::Standard,
        pruning: true,
    };

    struct MolDef {
        name: &'static str,
        atoms_fn: fn() -> Vec<Atom>,
    }

    let molecules = [
        MolDef {
            name: "H2",
            atoms_fn: h2_atoms,
        },
        MolDef {
            name: "HF",
            atoms_fn: hf_atoms,
        },
        MolDef {
            name: "H2O",
            atoms_fn: h2o_atoms,
        },
        MolDef {
            name: "NH3",
            atoms_fn: nh3_atoms,
        },
        MolDef {
            name: "CH4",
            atoms_fn: ch4_atoms,
        },
    ];

    let basis_sets = ["sto-3g", "6-31g*"];
    let total = molecules.len() * basis_sets.len() * 3; // RHF + LDA + B3LYP
    let mut count = 0;
    let mut results: Vec<ComparisonEntry> = Vec::with_capacity(total);

    for mol in &molecules {
        let atoms = (mol.atoms_fn)();
        for &basis_name in &basis_sets {
            // Build system and integrals (shared by RHF)
            eprintln!(
                "\n[{}/{}] {} / {} — building integrals...",
                count + 1,
                total,
                mol.name,
                basis_name
            );
            let (system, basis) = build_system(&atoms, basis_name);
            let n_occ = system.nelec / 2;
            let nbf = system.nbf;

            // ============================================================
            // RHF gradient
            // ============================================================
            count += 1;
            eprintln!("  [{}/{}] RHF gradient...", count, total);

            let rhf_output = rhf_scf(&system, &scf_config).expect("RHF SCF convergence");
            let rhf_grad = rhf_gradient(
                &basis,
                &rhf_output.density_matrix,
                &rhf_output.mo_coefficients,
                &rhf_output.mo_energies,
                n_occ,
            );

            let pyscf_rhf = find_pyscf_ref(&pyscf_refs, mol.name, basis_name, "RHF");
            let (max_err, rms_err) = pyscf_rhf
                .map(|r| gradient_errors(&rhf_grad.gradients, &r.gradient))
                .unzip();

            eprintln!(
                "    E={:.10}  max_grad={:.6e}  max_err={:.2e}",
                rhf_output.energy_total,
                rhf_grad.max_gradient,
                max_err.unwrap_or(f64::NAN)
            );

            results.push(ComparisonEntry {
                molecule: mol.name.to_string(),
                basis: basis_name.to_string(),
                method: "RHF".to_string(),
                n_basis: nbf,
                iqcp_energy: rhf_output.energy_total,
                pyscf_energy: pyscf_rhf.map(|r| r.energy),
                energy_error: pyscf_rhf.map(|r| (rhf_output.energy_total - r.energy).abs()),
                iqcp_gradient: grad_to_vecs(&rhf_grad),
                pyscf_gradient: pyscf_rhf.map(|r| r.gradient.clone()),
                max_abs_error: max_err,
                rms_abs_error: rms_err,
                iqcp_max_gradient: rhf_grad.max_gradient,
                pyscf_max_gradient: pyscf_rhf.map(|r| r.max_gradient),
            });

            // ============================================================
            // LDA gradient
            // ============================================================
            count += 1;
            eprintln!("  [{}/{}] LDA gradient...", count, total);

            let lda = Lda::new();
            let lda_grad = ks_dft_gradient(
                &atoms,
                basis_name,
                &lda,
                &grid_config,
                &scf_config,
                false, // no D3-BJ
            );

            // No PySCF LDA reference in gradient_validation.json,
            // so we report IQCP values only
            let pyscf_lda = find_pyscf_ref(&pyscf_refs, mol.name, basis_name, "LDA");
            let (lda_max_err, lda_rms_err) = pyscf_lda
                .map(|r| gradient_errors(&lda_grad.gradients, &r.gradient))
                .unzip();

            eprintln!(
                "    E={:.10}  max_grad={:.6e}  max_err={:.2e}",
                lda_grad.energy.unwrap_or(f64::NAN),
                lda_grad.max_gradient,
                lda_max_err.unwrap_or(f64::NAN)
            );

            results.push(ComparisonEntry {
                molecule: mol.name.to_string(),
                basis: basis_name.to_string(),
                method: "LDA".to_string(),
                n_basis: nbf,
                iqcp_energy: lda_grad.energy.unwrap_or(f64::NAN),
                pyscf_energy: pyscf_lda.map(|r| r.energy),
                energy_error: pyscf_lda.and_then(|r| lda_grad.energy.map(|e| (e - r.energy).abs())),
                iqcp_gradient: grad_to_vecs(&lda_grad),
                pyscf_gradient: pyscf_lda.map(|r| r.gradient.clone()),
                max_abs_error: lda_max_err,
                rms_abs_error: lda_rms_err,
                iqcp_max_gradient: lda_grad.max_gradient,
                pyscf_max_gradient: pyscf_lda.map(|r| r.max_gradient),
            });

            // ============================================================
            // B3LYP gradient
            // ============================================================
            count += 1;
            eprintln!("  [{}/{}] B3LYP gradient...", count, total);

            let b3lyp = B3lyp::new();
            let b3lyp_grad = ks_dft_gradient(
                &atoms,
                basis_name,
                &b3lyp,
                &grid_config,
                &scf_config,
                false, // no D3-BJ
            );

            let pyscf_b3lyp = find_pyscf_ref(&pyscf_refs, mol.name, basis_name, "B3LYP");
            let (b3_max_err, b3_rms_err) = pyscf_b3lyp
                .map(|r| gradient_errors(&b3lyp_grad.gradients, &r.gradient))
                .unzip();

            eprintln!(
                "    E={:.10}  max_grad={:.6e}  max_err={:.2e}",
                b3lyp_grad.energy.unwrap_or(f64::NAN),
                b3lyp_grad.max_gradient,
                b3_max_err.unwrap_or(f64::NAN)
            );

            results.push(ComparisonEntry {
                molecule: mol.name.to_string(),
                basis: basis_name.to_string(),
                method: "B3LYP".to_string(),
                n_basis: nbf,
                iqcp_energy: b3lyp_grad.energy.unwrap_or(f64::NAN),
                pyscf_energy: pyscf_b3lyp.map(|r| r.energy),
                energy_error: pyscf_b3lyp
                    .and_then(|r| b3lyp_grad.energy.map(|e| (e - r.energy).abs())),
                iqcp_gradient: grad_to_vecs(&b3lyp_grad),
                pyscf_gradient: pyscf_b3lyp.map(|r| r.gradient.clone()),
                max_abs_error: b3_max_err,
                rms_abs_error: b3_rms_err,
                iqcp_max_gradient: b3lyp_grad.max_gradient,
                pyscf_max_gradient: pyscf_b3lyp.map(|r| r.max_gradient),
            });
        }
    }

    // ========================================================================
    // Summary statistics to stderr
    // ========================================================================
    eprintln!("\n========================================");
    eprintln!("Gradient Comparison Summary");
    eprintln!("========================================");
    eprintln!(
        "{:<6} {:<8} {:<8} {:>12} {:>12} {:>12}",
        "Mol", "Basis", "Method", "Max|dE/dR|", "MaxErr", "RMS Err"
    );
    eprintln!("{}", "-".repeat(70));

    for r in &results {
        eprintln!(
            "{:<6} {:<8} {:<8} {:>12.6e} {:>12.2e} {:>12.2e}",
            r.molecule,
            r.basis,
            r.method,
            r.iqcp_max_gradient,
            r.max_abs_error.unwrap_or(f64::NAN),
            r.rms_abs_error.unwrap_or(f64::NAN),
        );
    }

    eprintln!("\nTotal systems: {}", results.len());

    // ========================================================================
    // JSON output to stdout
    // ========================================================================
    let json = serde_json::to_string_pretty(&results).expect("Serialize results");
    println!("{}", json);
}
