//! Density Functional Theory (DFT) numerical integration infrastructure.
//!
//! This module provides the Becke atom-centered numerical integration grid
//! for evaluating exchange-correlation functionals in KS-DFT calculations.
//!
//! # Components
//!
//! - **Grid construction** (`grid`): Becke partitioning, Mura-Knowles radial
//!   quadrature, SG-1 pruning, and full grid assembly.
//! - **Lebedev quadrature** (`lebedev`): Angular quadrature data tables for
//!   6, 38, 86, 194, 302, and 590-point grids on the unit sphere.
//!
//! # Algorithm
//!
//! The grid construction follows Becke's atom-centered scheme (1988):
//!
//! 1. For each atom, generate a radial grid using the Mura-Knowles log3
//!    mapping (1996).
//! 2. At each radial point, select an angular Lebedev grid using SG-1
//!    pruning (Gill, Johnson & Pople, 1993).
//! 3. Combine radial and angular grids into 3D Cartesian coordinates.
//! 4. Apply Becke's smooth partitioning with 3-iteration hardening to
//!    distribute weights among atoms.
//!
//! # References
//!
//! - Becke, A. D. (1988). J. Chem. Phys. 88, 2547. (partitioning)
//! - Mura, M. E. & Knowles, P. J. (1996). J. Chem. Phys. 104, 9848. (radial)
//! - Lebedev, V. I. & Laikov, D. N. (1999). Doklady Math. 59, 477. (angular)
//! - Gill, P. M. W., Johnson, B. G. & Pople, J. A. (1993). CPL 209, 506. (SG-1)
//!
//! # Example
//!
//! ```rust
//! use qc_core::basis::{Atom, BasisSet};
//! use qc_core::dft::{build_becke_grid, GridConfig, GridQuality};
//!
//! let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
//! let h2 = Atom::new(1, [0.0, 0.0, 1.4]).unwrap();
//! let basis = BasisSet::build(vec![h1, h2], "sto-3g").unwrap();
//!
//! let config = GridConfig::default();
//! let grid = build_becke_grid(&basis.atoms, &config);
//!
//! assert!(grid.n_points > 0);
//! assert_eq!(grid.n_atoms, 2);
//! ```

pub mod dispersion;
pub mod functionals;
pub mod grid;
pub mod ks_scf;
#[allow(clippy::excessive_precision, clippy::unreadable_literal)]
pub mod lebedev;

use serde::{Deserialize, Serialize};

// Re-export primary public API
pub use dispersion::{
    compute_d3bj_energy, compute_d3bj_gradient, D3bjParams, D3bjResult, D3BJ_B3LYP,
};
pub use functionals::{
    B3lyp, Becke88Exchange, DftFunctional, ExchangeCorrelation, Lda, LypCorrelation,
    SlaterExchange, Vwn5Correlation, DENSITY_THRESHOLD,
};
pub use grid::build_becke_grid;
pub use ks_scf::{ks_scf, ks_scf_with_guess, KsScfOutput};

/// Module version (matches crate version).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// =============================================================================
// Data Structures
// =============================================================================

/// Becke atom-centered numerical integration grid.
///
/// Contains quadrature points and weights for DFT numerical integration,
/// constructed using Mura-Knowles radial quadrature, Lebedev angular
/// quadrature, Becke partitioning, and SG-1 pruning.
///
/// # Weight composition
///
/// Each weight incorporates three factors:
/// - Radial weight: `4*pi * r^2 * dr` (spherical volume element)
/// - Angular weight: Lebedev quadrature weight (sums to 4*pi)
/// - Becke partition weight: smooth space-partitioning factor in \[0, 1\]
///
/// The final weight for grid point g belonging to atom A is:
/// ```text
/// w_g = w_radial * w_angular / (4*pi) * P_A(r_g)
/// ```
/// where `P_A` is the normalized Becke partition weight.
///
/// # Reference
///
/// Becke (1988) JCP 88, 2547; Gill, Johnson & Pople (1993) CPL 209, 506.
#[derive(Debug, Clone)]
pub struct BeckeGrid {
    /// Grid point coordinates \[x, y, z\] in bohr.
    pub points: Vec<[f64; 3]>,
    /// Quadrature weights (includes Becke partitioning and radial/angular weights).
    pub weights: Vec<f64>,
    /// Which atom each grid point belongs to (index into atoms array).
    pub atom_indices: Vec<usize>,
    /// Number of atoms in the molecule.
    pub n_atoms: usize,
    /// Total number of grid points.
    pub n_points: usize,
    /// Grid quality level used.
    pub quality: GridQuality,
}

/// Grid quality presets controlling angular quadrature resolution.
///
/// The quality setting determines the maximum Lebedev angular grid used
/// in the SG-1 pruning scheme's "standard region" (region 4).
///
/// | Quality | Region-4 Angular | Approx pts/atom |
/// |---------|------------------|-----------------|
/// | Standard | Lebedev-302 | ~5000 |
/// | Fine | Lebedev-590 | ~8000 |
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
pub enum GridQuality {
    /// Standard quality: SG-1 pruning with Lebedev-302 in the main region.
    Standard,
    /// Fine quality: SG-1 pruning with Lebedev-590 in the main region.
    Fine,
}

/// Configuration for Becke grid generation.
#[derive(Debug, Clone)]
pub struct GridConfig {
    /// Number of radial points per atom (default: 75).
    pub n_radial: usize,
    /// Angular grid quality.
    pub quality: GridQuality,
    /// Whether to apply SG-1 pruning (default: true).
    pub pruning: bool,
}

impl Default for GridConfig {
    fn default() -> Self {
        Self {
            n_radial: 75,
            quality: GridQuality::Standard,
            pruning: true,
        }
    }
}
