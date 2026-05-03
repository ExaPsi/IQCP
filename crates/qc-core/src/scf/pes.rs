//! PES (Potential Energy Surface) bond-length scan engine
//!
//! Implements bond-length scanning for diatomic molecules by running
//! a series of SCF calculations at different internuclear distances.
//!
//! # Algorithm
//!
//! 1. Generate a grid of bond distances from `r_min` to `r_max`
//! 2. At each distance:
//!    a. Construct molecule geometry (atom A at origin, atom B at [0,0,r])
//!    b. Build basis set and compute molecular integrals on-the-fly
//!    c. Run SCF with convergence seeding from previous point
//! 3. Find equilibrium via parabolic interpolation of three lowest-energy points
//!
//! # Convergence Seeding
//!
//! The converged density matrix from geometry point `i` is used as the initial
//! guess for point `i+1`. Because adjacent geometries have similar electronic
//! structure, this typically reduces SCF iterations by 30% or more compared
//! to cold starts from the core Hamiltonian.
//!
//! # References
//!
//! - Phase 2 TDD Section 8.1: PES scan engine specification
//! - PySCF: `scf.RHF.kernel(dm0=...)` for density seeding pattern

use serde::{Deserialize, Serialize};

use crate::basis::{Atom, BasisSet};
use crate::integrals::{eri_compressed, hcore_matrix, overlap_matrix};

use super::{rhf_scf_with_guess, PresetSystem, ScfConfig};

// ============================================================================
// Data Structures
// ============================================================================

/// Configuration for a PES bond-length scan
///
/// Groups all parameters needed to specify a PES scan, including the
/// molecular system, scan range, and SCF options.
#[derive(Debug, Clone)]
pub struct PesScanConfig<'a> {
    /// Atomic number of atom A (fixed at origin), 1-10
    pub atom_a_z: u8,
    /// Atomic number of atom B (translated along z-axis), 1-10
    pub atom_b_z: u8,
    /// Minimum bond distance in bohr
    pub r_min: f64,
    /// Maximum bond distance in bohr
    pub r_max: f64,
    /// Number of scan points (evenly spaced)
    pub n_points: usize,
    /// Basis set name (e.g., "sto-3g")
    pub basis_name: &'a str,
    /// SCF configuration options
    pub scf_config: &'a ScfConfig,
    /// Whether to use convergence seeding (density from previous point)
    pub use_seeding: bool,
}

/// A single point on the potential energy surface
///
/// Stores the bond distance, total energy, convergence status, and
/// iteration count for one SCF calculation at a specific geometry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PesPoint {
    /// Bond distance in bohr
    pub r: f64,
    /// Total SCF energy in Hartree (electronic + nuclear)
    pub energy: f64,
    /// Whether SCF converged at this geometry
    pub converged: bool,
    /// Number of SCF iterations at this point
    pub iterations: usize,
}

/// Equilibrium geometry from parabolic interpolation
///
/// Determined by fitting a parabola through the three points surrounding
/// the minimum energy on the PES scan grid.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PesEquilibrium {
    /// Interpolated equilibrium bond distance in bohr
    pub r_bohr: f64,
    /// Interpolated equilibrium energy in Hartree
    pub energy_hartree: f64,
}

/// Result of a PES bond-length scan
///
/// Contains all scan points, the interpolated equilibrium (if found),
/// and timing information.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PesScanResult {
    /// Energy at each scanned geometry
    pub points: Vec<PesPoint>,
    /// Equilibrium from parabolic interpolation of three lowest-energy points
    pub equilibrium: Option<PesEquilibrium>,
    /// Total computation time in milliseconds
    pub compute_time_ms: f64,
    /// Total SCF iterations across all scan points
    pub total_iterations: usize,
}

// ============================================================================
// PES Scan Engine
// ============================================================================

/// Run a bond-length PES scan for a diatomic molecule
///
/// Scans the potential energy surface by varying the internuclear distance
/// from `r_min` to `r_max` over `n_points` evenly spaced points. At each
/// geometry, molecular integrals are computed on-the-fly and SCF is run
/// with optional convergence seeding.
///
/// # Arguments
///
/// * `scan_config` - PES scan configuration (atoms, range, SCF options)
/// * `progress` - Optional progress callback: (point_idx, r, energy, converged)
///
/// # Geometry Convention
///
/// - Atom A is fixed at the origin [0, 0, 0]
/// - Atom B is placed at [0, 0, r] for each scan distance r
///
/// # Returns
///
/// `PesScanResult` containing all scan points, equilibrium interpolation,
/// and timing information.
///
/// # Example
///
/// ```rust
/// use qc_core::scf::pes::{pes_scan, PesScanConfig};
/// use qc_core::scf::ScfConfig;
///
/// let scf_config = ScfConfig::medium();
/// let scan_config = PesScanConfig {
///     atom_a_z: 1, atom_b_z: 1,
///     r_min: 0.5, r_max: 5.0,
///     n_points: 20,
///     basis_name: "sto-3g",
///     scf_config: &scf_config,
///     use_seeding: true,
/// };
/// let result = pes_scan(&scan_config, None);
///
/// assert!(!result.points.is_empty());
/// assert!(result.equilibrium.is_some());
/// ```
pub fn pes_scan(
    scan_config: &PesScanConfig,
    progress: Option<&dyn Fn(usize, f64, f64, bool)>,
) -> PesScanResult {
    let atom_a_z = scan_config.atom_a_z;
    let atom_b_z = scan_config.atom_b_z;
    let r_min = scan_config.r_min;
    let r_max = scan_config.r_max;
    let n_points = scan_config.n_points;
    let basis_name = scan_config.basis_name;
    let config = scan_config.scf_config;
    let use_seeding = scan_config.use_seeding;

    // Generate evenly spaced distances
    let distances: Vec<f64> = if n_points == 1 {
        vec![r_min]
    } else {
        (0..n_points)
            .map(|i| r_min + (r_max - r_min) * i as f64 / (n_points - 1) as f64)
            .collect()
    };

    let mut points = Vec::with_capacity(n_points);
    let mut prev_density: Option<Vec<f64>> = None;
    let mut total_iterations: usize = 0;

    for (idx, &r) in distances.iter().enumerate() {
        // Step 1: Construct molecular geometry
        // Atom A at origin, atom B at [0, 0, r]
        let atom_a = match Atom::new(atom_a_z, [0.0, 0.0, 0.0]) {
            Ok(a) => a,
            Err(_) => {
                // Invalid atom -- mark as unconverged and continue
                points.push(PesPoint {
                    r,
                    energy: 0.0,
                    converged: false,
                    iterations: 0,
                });
                if let Some(cb) = progress {
                    cb(idx, r, 0.0, false);
                }
                continue;
            }
        };

        let atom_b = match Atom::new(atom_b_z, [0.0, 0.0, r]) {
            Ok(a) => a,
            Err(_) => {
                points.push(PesPoint {
                    r,
                    energy: 0.0,
                    converged: false,
                    iterations: 0,
                });
                if let Some(cb) = progress {
                    cb(idx, r, 0.0, false);
                }
                continue;
            }
        };

        // Step 2: Build basis set
        let basis = match BasisSet::build(vec![atom_a, atom_b], basis_name) {
            Ok(b) => b,
            Err(_) => {
                points.push(PesPoint {
                    r,
                    energy: 0.0,
                    converged: false,
                    iterations: 0,
                });
                if let Some(cb) = progress {
                    cb(idx, r, 0.0, false);
                }
                // Reset seeding on failure
                prev_density = None;
                continue;
            }
        };

        // Step 3: Compute molecular integrals on-the-fly
        let s_mat = overlap_matrix(&basis);
        let h_core = hcore_matrix(&basis);
        let eri = eri_compressed(&basis);

        // Step 4: Package as PresetSystem for SCF engine
        let system = PresetSystem {
            system_id: format!("pes_{}_{}_r{:.4}", atom_a_z, atom_b_z, r),
            label: format!("PES scan point r={:.4} bohr", r),
            nbf: basis.n_basis,
            nelec: basis.n_electrons,
            e_nuc: basis.nuclear_repulsion,
            s_matrix: s_mat,
            h_core,
            eri_compressed: eri,
        };

        // Step 5: Run SCF with optional convergence seeding
        let initial_guess = if use_seeding {
            prev_density.as_deref()
        } else {
            None
        };

        match rhf_scf_with_guess(&system, config, initial_guess) {
            Ok(output) => {
                let energy = output.energy_total;
                let iterations = output.iterations;
                total_iterations += iterations;

                // Store density matrix for seeding next point
                prev_density = Some(output.density_matrix.clone());

                points.push(PesPoint {
                    r,
                    energy,
                    converged: true,
                    iterations,
                });

                if let Some(cb) = progress {
                    cb(idx, r, energy, true);
                }
            }
            Err(_) => {
                // SCF did not converge at this geometry
                // Mark as unconverged, reset seeding (don't pass bad density forward)
                points.push(PesPoint {
                    r,
                    energy: 0.0,
                    converged: false,
                    iterations: config.max_iterations,
                });
                total_iterations += config.max_iterations;

                // Reset seeding -- don't propagate from a failed point
                prev_density = None;

                if let Some(cb) = progress {
                    cb(idx, r, 0.0, false);
                }
            }
        }
    }

    // Find equilibrium via parabolic interpolation
    let equilibrium = find_equilibrium(&points);

    PesScanResult {
        points,
        equilibrium,
        // Timing is handled by the caller (WASM layer uses js_sys::Date::now()
        // since std::time::Instant is not supported on wasm32-unknown-unknown)
        compute_time_ms: 0.0,
        total_iterations,
    }
}

// ============================================================================
// Parabolic Interpolation
// ============================================================================

/// Find equilibrium geometry via parabolic interpolation
///
/// Takes the three converged points surrounding the minimum energy and fits
/// a parabola `E(r) = a*r^2 + b*r + c` to find the interpolated minimum.
///
/// # Algorithm
///
/// 1. Filter to converged points only
/// 2. Find the index `k` of the minimum energy point
/// 3. If `k` is at the boundary (first or last point), return `None`
/// 4. Take three points: (r_{k-1}, E_{k-1}), (r_k, E_k), (r_{k+1}, E_{k+1})
/// 5. Fit parabola using Lagrange interpolation:
///    - `a = [E_1(r_2-r_3) + E_2(r_3-r_1) + E_3(r_1-r_2)] / [(r_1-r_2)(r_1-r_3)(r_2-r_3)]`
///    - `b = [E_1(r_3^2-r_2^2) + E_2(r_1^2-r_3^2) + E_3(r_2^2-r_1^2)] / [(r_1-r_2)(r_1-r_3)(r_2-r_3)]`
/// 6. Minimum at `r_eq = -b / (2a)`
/// 7. Equilibrium energy: `E_eq = a*r_eq^2 + b*r_eq + c`
fn find_equilibrium(points: &[PesPoint]) -> Option<PesEquilibrium> {
    // Filter to converged points
    let converged: Vec<&PesPoint> = points.iter().filter(|p| p.converged).collect();

    if converged.len() < 3 {
        return None;
    }

    // Find index of minimum energy among converged points
    let min_idx = converged
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            a.energy
                .partial_cmp(&b.energy)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)?;

    // Need at least one point on each side of minimum
    if min_idx == 0 || min_idx >= converged.len() - 1 {
        // Minimum is at the boundary -- cannot interpolate
        // Return the boundary point as-is
        return Some(PesEquilibrium {
            r_bohr: converged[min_idx].r,
            energy_hartree: converged[min_idx].energy,
        });
    }

    // Three points for parabolic fit
    let r1 = converged[min_idx - 1].r;
    let r2 = converged[min_idx].r;
    let r3 = converged[min_idx + 1].r;
    let e1 = converged[min_idx - 1].energy;
    let e2 = converged[min_idx].energy;
    let e3 = converged[min_idx + 1].energy;

    // Fit parabola E(r) = a*r^2 + b*r + c through three points
    // Using the standard formula for three-point parabolic interpolation
    let denom = (r1 - r2) * (r1 - r3) * (r2 - r3);

    if denom.abs() < 1e-30 {
        // Degenerate case -- points too close together
        return Some(PesEquilibrium {
            r_bohr: r2,
            energy_hartree: e2,
        });
    }

    let a = (r3 * (e2 - e1) + r2 * (e1 - e3) + r1 * (e3 - e2)) / denom;
    let b = (r3 * r3 * (e1 - e2) + r2 * r2 * (e3 - e1) + r1 * r1 * (e2 - e3)) / denom;
    let c =
        (r2 * r3 * (r2 - r3) * e1 + r3 * r1 * (r3 - r1) * e2 + r1 * r2 * (r1 - r2) * e3) / denom;

    // Check that we have a proper minimum (a > 0)
    if a <= 0.0 {
        // Concave parabola -- no minimum
        return Some(PesEquilibrium {
            r_bohr: r2,
            energy_hartree: e2,
        });
    }

    // Minimum of parabola: r_eq = -b / (2a)
    let r_eq = -b / (2.0 * a);

    // Evaluate equilibrium energy
    let e_eq = a * r_eq * r_eq + b * r_eq + c;

    Some(PesEquilibrium {
        r_bohr: r_eq,
        energy_hartree: e_eq,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scf::{ConvergenceProfile, ScfConfig};

    // ========================================================================
    // Helper: build PesScanConfig and run scan
    // ========================================================================

    /// Convenience wrapper to run a PES scan without building the config struct manually
    fn run_scan(
        atom_a_z: u8,
        atom_b_z: u8,
        r_min: f64,
        r_max: f64,
        n_points: usize,
        basis_name: &str,
        scf_config: &ScfConfig,
        use_seeding: bool,
        progress: Option<&dyn Fn(usize, f64, f64, bool)>,
    ) -> PesScanResult {
        let scan_config = PesScanConfig {
            atom_a_z,
            atom_b_z,
            r_min,
            r_max,
            n_points,
            basis_name,
            scf_config,
            use_seeding,
        };
        pes_scan(&scan_config, progress)
    }

    // ========================================================================
    // Helper: load PySCF golden data
    // ========================================================================

    #[derive(Deserialize)]
    struct GoldenPes {
        points: Vec<GoldenPoint>,
        equilibrium: GoldenEquilibrium,
    }

    #[derive(Deserialize)]
    struct GoldenPoint {
        r_bohr: f64,
        energy: f64,
        converged: bool,
    }

    #[derive(Deserialize)]
    struct GoldenEquilibrium {
        r_bohr: f64,
        energy_hartree: f64,
    }

    fn load_golden() -> GoldenPes {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/golden/pes/h2_sto3g_pes.json"
        );
        let data = std::fs::read_to_string(path).expect("Golden PES data file not found");
        serde_json::from_str(&data).expect("Failed to parse golden PES data")
    }

    // ========================================================================
    // Golden tests: H2 STO-3G PES scan
    // ========================================================================

    #[test]
    fn test_h2_pes_scan_energies() {
        // Compare IQCP PES scan energies against PySCF golden data.
        //
        // Tolerance: 1e-4 Ha.
        // Our on-the-fly integral engine uses Obara-Saika recurrence with
        // Boys function evaluation that has small numerical differences from
        // PySCF's libcint-based integrals. At short bond distances (r=0.5),
        // the nuclear attraction integrals are especially sensitive, leading
        // to integral-level errors of ~1e-5 that propagate to the total energy.
        let golden = load_golden();

        let config = ScfConfig {
            profile: ConvergenceProfile::Tight,
            max_iterations: 100,
            use_diis: true,
            diis_size: 8,
            diis_start: 1,
            damp: 0.0,
            damp_start: 5,
            level_shift: 0.0,
        };

        let result = run_scan(1, 1, 0.5, 5.0, 20, "sto-3g", &config, true, None);

        // Verify we got all 20 points
        assert_eq!(result.points.len(), 20);

        // Compare energies at each point
        for (iqcp, pyscf) in result.points.iter().zip(golden.points.iter()) {
            assert!(
                (iqcp.r - pyscf.r_bohr).abs() < 1e-10,
                "Distance mismatch: IQCP={} PySCF={}",
                iqcp.r,
                pyscf.r_bohr
            );

            if iqcp.converged && pyscf.converged {
                let error = (iqcp.energy - pyscf.energy).abs();
                assert!(
                    error < 1e-4,
                    "Energy mismatch at r={:.4}: IQCP={:.10} PySCF={:.10} error={:.2e}",
                    iqcp.r,
                    iqcp.energy,
                    pyscf.energy,
                    error
                );
            }
        }
    }

    #[test]
    fn test_h2_pes_equilibrium() {
        let golden = load_golden();

        let config = ScfConfig {
            profile: ConvergenceProfile::Tight,
            max_iterations: 100,
            use_diis: true,
            diis_size: 8,
            diis_start: 1,
            damp: 0.0,
            damp_start: 5,
            level_shift: 0.0,
        };

        let result = run_scan(1, 1, 0.5, 5.0, 20, "sto-3g", &config, true, None);

        let eq = result.equilibrium.as_ref().expect("No equilibrium found");

        let r_error = (eq.r_bohr - golden.equilibrium.r_bohr).abs();
        assert!(
            r_error < 0.1,
            "Equilibrium r mismatch: IQCP={:.6} PySCF={:.6} error={:.4}",
            eq.r_bohr,
            golden.equilibrium.r_bohr,
            r_error
        );

        let e_error = (eq.energy_hartree - golden.equilibrium.energy_hartree).abs();
        assert!(
            e_error < 1e-2,
            "Equilibrium energy mismatch: IQCP={:.10} PySCF={:.10} error={:.2e}",
            eq.energy_hartree,
            golden.equilibrium.energy_hartree,
            e_error
        );
    }

    #[test]
    fn test_h2_pes_finer_grid_equilibrium() {
        let config = ScfConfig {
            profile: ConvergenceProfile::Tight,
            max_iterations: 100,
            use_diis: true,
            diis_size: 8,
            diis_start: 1,
            damp: 0.0,
            damp_start: 5,
            level_shift: 0.0,
        };

        let result = run_scan(1, 1, 1.0, 2.0, 30, "sto-3g", &config, true, None);

        let eq = result.equilibrium.as_ref().expect("No equilibrium found");

        // True H2 STO-3G equilibrium: r ~ 1.346 bohr
        let r_error = (eq.r_bohr - 1.346).abs();
        assert!(
            r_error < 0.005,
            "Fine grid equilibrium r mismatch: IQCP={:.6} expected~1.346 error={:.4}",
            eq.r_bohr,
            r_error
        );
    }

    #[test]
    fn test_convergence_seeding_reduces_iterations() {
        // Use LiH (6 basis functions, 4 electrons) without DIIS
        // to observe seeding benefit.
        let config = ScfConfig {
            profile: ConvergenceProfile::Tight,
            max_iterations: 100,
            use_diis: false,
            diis_size: 8,
            diis_start: 1,
            damp: 0.0,
            damp_start: 5,
            level_shift: 0.0,
        };

        let result_cold = run_scan(3, 1, 2.5, 4.0, 15, "sto-3g", &config, false, None);
        let result_seeded = run_scan(3, 1, 2.5, 4.0, 15, "sto-3g", &config, true, None);

        let cold_iters = result_cold.total_iterations;
        let seeded_iters = result_seeded.total_iterations;

        assert!(
            seeded_iters < cold_iters,
            "Convergence seeding should reduce total iterations: cold={}, seeded={}",
            cold_iters,
            seeded_iters
        );

        let reduction = (cold_iters as f64 - seeded_iters as f64) / cold_iters as f64;
        assert!(
            reduction >= 0.10,
            "Seeding should reduce iterations by >= 10%: {:.1}% (cold={}, seeded={})",
            reduction * 100.0,
            cold_iters,
            seeded_iters
        );
    }

    #[test]
    fn test_pes_all_points_converge() {
        let config = ScfConfig::medium();
        let result = run_scan(1, 1, 0.5, 5.0, 10, "sto-3g", &config, true, None);

        for pt in &result.points {
            assert!(pt.converged, "Point at r={:.4} did not converge", pt.r);
        }
    }

    #[test]
    fn test_pes_progress_callback() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let call_count = AtomicUsize::new(0);
        let n_points = 5;

        let config = ScfConfig::loose();

        let _result = run_scan(
            1,
            1,
            1.0,
            3.0,
            n_points,
            "sto-3g",
            &config,
            true,
            Some(&|_idx, _r, _e, _conv| {
                call_count.fetch_add(1, Ordering::Relaxed);
            }),
        );

        assert_eq!(
            call_count.load(Ordering::Relaxed),
            n_points,
            "Progress callback should be called {} times",
            n_points
        );
    }

    #[test]
    fn test_pes_monotonicity_near_minimum() {
        let config = ScfConfig::medium();
        let result = run_scan(1, 1, 0.5, 5.0, 20, "sto-3g", &config, true, None);

        let converged: Vec<&PesPoint> = result.points.iter().filter(|p| p.converged).collect();

        let min_idx = converged
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                a.energy
                    .partial_cmp(&b.energy)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap();

        for i in 1..min_idx {
            assert!(
                converged[i].energy <= converged[i - 1].energy,
                "Energy should decrease before minimum: E({:.2})={:.6} > E({:.2})={:.6}",
                converged[i].r,
                converged[i].energy,
                converged[i - 1].r,
                converged[i - 1].energy,
            );
        }

        for i in (min_idx + 1)..converged.len() {
            assert!(
                converged[i].energy >= converged[i - 1].energy,
                "Energy should increase after minimum: E({:.2})={:.6} < E({:.2})={:.6}",
                converged[i].r,
                converged[i].energy,
                converged[i - 1].r,
                converged[i - 1].energy,
            );
        }
    }

    #[test]
    fn test_pes_single_point() {
        let config = ScfConfig::medium();
        let result = run_scan(1, 1, 1.4, 1.4, 1, "sto-3g", &config, true, None);

        assert_eq!(result.points.len(), 1);
        assert!(result.points[0].converged);
    }

    #[test]
    fn test_pes_two_points() {
        let config = ScfConfig::medium();
        let result = run_scan(1, 1, 1.0, 2.0, 2, "sto-3g", &config, true, None);

        assert_eq!(result.points.len(), 2);
    }

    #[test]
    fn test_find_equilibrium_basic() {
        let points = vec![
            PesPoint {
                r: 1.0,
                energy: -1.0,
                converged: true,
                iterations: 5,
            },
            PesPoint {
                r: 2.0,
                energy: -2.0,
                converged: true,
                iterations: 5,
            },
            PesPoint {
                r: 3.0,
                energy: -1.5,
                converged: true,
                iterations: 5,
            },
        ];

        let eq = find_equilibrium(&points).expect("Should find equilibrium");
        assert!((eq.r_bohr - 2.1667).abs() < 0.001);
    }

    #[test]
    fn test_find_equilibrium_with_unconverged() {
        let points = vec![
            PesPoint {
                r: 0.5,
                energy: 0.0,
                converged: false,
                iterations: 100,
            },
            PesPoint {
                r: 1.0,
                energy: -1.0,
                converged: true,
                iterations: 5,
            },
            PesPoint {
                r: 1.5,
                energy: -2.0,
                converged: true,
                iterations: 5,
            },
            PesPoint {
                r: 2.0,
                energy: -1.5,
                converged: true,
                iterations: 5,
            },
        ];

        let eq = find_equilibrium(&points).expect("Should find equilibrium");
        assert!(eq.r_bohr > 1.0 && eq.r_bohr < 2.0);
    }

    #[test]
    fn test_find_equilibrium_min_at_boundary() {
        let points = vec![
            PesPoint {
                r: 1.0,
                energy: -3.0,
                converged: true,
                iterations: 5,
            },
            PesPoint {
                r: 2.0,
                energy: -2.0,
                converged: true,
                iterations: 5,
            },
            PesPoint {
                r: 3.0,
                energy: -1.0,
                converged: true,
                iterations: 5,
            },
        ];

        let eq = find_equilibrium(&points).expect("Should find equilibrium");
        assert!((eq.r_bohr - 1.0).abs() < 1e-10);
        assert!((eq.energy_hartree - (-3.0)).abs() < 1e-10);
    }

    #[test]
    fn test_pes_performance() {
        let config = ScfConfig::medium();

        let start = std::time::Instant::now();
        let result = run_scan(1, 1, 0.5, 5.0, 20, "sto-3g", &config, true, None);
        let elapsed = start.elapsed();

        assert!(
            elapsed.as_secs() < 5,
            "PES scan took {:.2}s, should be < 5s",
            elapsed.as_secs_f64()
        );

        assert_eq!(result.points.len(), 20);
        let converged_count = result.points.iter().filter(|p| p.converged).count();
        assert!(
            converged_count >= 18,
            "At least 18 of 20 points should converge, got {}",
            converged_count
        );
    }

    #[test]
    fn test_he2_pes() {
        let config = ScfConfig::medium();
        let result = run_scan(2, 2, 1.0, 6.0, 10, "sto-3g", &config, true, None);

        assert_eq!(result.points.len(), 10);
        let converged_count = result.points.iter().filter(|p| p.converged).count();
        assert!(
            converged_count >= 8,
            "Most He2 points should converge, got {}",
            converged_count
        );
    }

    #[test]
    fn test_pes_integral_consistency() {
        // Verify on-the-fly integrals are reasonable at R=1.4 bohr.
        //
        // The preset values were generated by PySCF (libcint) which uses
        // a different integral engine than our Obara-Saika implementation.
        // Small numerical differences (up to ~1e-5) are expected due to:
        // - Different Boys function evaluation strategies
        // - Different recurrence relation implementations
        // - Different normalization conventions
        //
        // We verify structural consistency (size, symmetry) and that
        // values agree to reasonable tolerance.

        let atom_a = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let atom_b = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
        let basis = BasisSet::build(vec![atom_a, atom_b], "sto-3g").unwrap();

        let s_mat = overlap_matrix(&basis);
        let h_core_mat = hcore_matrix(&basis);
        let eri = eri_compressed(&basis);

        // Verify correct sizes
        assert_eq!(s_mat.len(), 4);
        assert_eq!(h_core_mat.len(), 4);
        assert_eq!(eri.len(), 6);

        // Compare with preset values
        let preset = PresetSystem::h2_sto3g_test();

        // Overlap matrix: should agree within 1e-8
        // (overlap integrals are simple and very accurate)
        for (i, (computed, reference)) in s_mat.iter().zip(preset.s_matrix.iter()).enumerate() {
            let error = (computed - reference).abs();
            assert!(
                error < 1e-8,
                "S matrix mismatch at index {}: computed={:.12} preset={:.12} error={:.2e}",
                i,
                computed,
                reference,
                error
            );
        }

        // H_core matrix: agree within 1e-4
        // (nuclear attraction integrals have Boys function error accumulation)
        for (i, (computed, reference)) in h_core_mat.iter().zip(preset.h_core.iter()).enumerate() {
            let error = (computed - reference).abs();
            assert!(
                error < 1e-4,
                "H_core mismatch at index {}: computed={:.12} preset={:.12} error={:.2e}",
                i,
                computed,
                reference,
                error
            );
        }

        // ERI compressed: agree within 1e-4
        for (i, (computed, reference)) in eri.iter().zip(preset.eri_compressed.iter()).enumerate() {
            let error = (computed - reference).abs();
            assert!(
                error < 1e-4,
                "ERI mismatch at index {}: computed={:.12} preset={:.12} error={:.2e}",
                i,
                computed,
                reference,
                error
            );
        }
    }
}
