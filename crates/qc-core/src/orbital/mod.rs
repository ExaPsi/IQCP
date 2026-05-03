//! Molecular orbital grid evaluation and isosurface extraction
//!
//! This module provides tools for evaluating molecular orbital wavefunctions
//! on 3D grids and extracting isosurfaces for visualization.
//!
//! # Overview
//!
//! After an SCF calculation produces MO coefficients C, the orbital
//! wavefunction at any point r is:
//!
//! ```text
//! psi_i(r) = sum_mu { C_{mu,i} * chi_mu(r) }
//! ```
//!
//! where `chi_mu(r)` is the contracted Gaussian basis function evaluated at r.
//!
//! # Modules
//!
//! - [`grid`]: MO grid evaluator with Gaussian screening and shell-batched evaluation
//! - [`marching_cubes`]: Marching cubes isosurface extraction for 3D visualization
//! - [`promolecule`]: Promolecule density and difference density computation
//!
//! # References
//!
//! - Szabo & Ostlund (1996). "Modern Quantum Chemistry". Dover, Chapter 3.
//! - PySCF cubegen implementation: `references/pyscf/pyscf/tools/cubegen.py`
//! - Lorensen & Cline (1987). "Marching Cubes." Computer Graphics, 21(4), 163-169.

pub mod density;
pub mod grid;
pub mod marching_cubes;
pub mod promolecule;

// Re-export public types
pub use density::{evaluate_density_on_grid, DensityGridResult};
pub use grid::{
    auto_grid_bounds, evaluate_mo_on_grid, GridConfig, MoGridInput, MoGridResult, OrbitalError,
};
pub use marching_cubes::{
    dual_marching_cubes, marching_cubes as extract_isosurface, MarchingCubesResult,
};
pub use promolecule::{
    compute_difference_density, evaluate_promolecule_on_grid, get_atomic_density,
    unsupported_elements, AtomicDensityProfile, DifferenceDensityResult,
};

/// Module version (matches crate version)
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
