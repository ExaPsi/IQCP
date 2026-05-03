//! Core quantum chemistry algorithms for IQCP
//!
//! This crate provides implementations of fundamental quantum chemistry
//! algorithms for the Interactive Quantum Chemistry Playground:
//!
//! - **Boys function** `F_m(T)` with regime-based evaluation
//! - **Rys quadrature** roots and weights computation
//! - **RHF SCF** with optional DIIS acceleration
//!
//! # Design Philosophy
//!
//! All algorithms are implemented in pure Rust with no WASM dependencies,
//! enabling use in both native and WebAssembly contexts. The implementations
//! prioritize:
//!
//! 1. **Numerical accuracy**: Tolerances documented and tested
//! 2. **Pedagogical clarity**: Code structure mirrors textbook presentation
//! 3. **Performance**: Optimized for interactive use (<200ms latency)
//!
//! # References
//!
//! - Boys function: Shavitt (1963), Methods in Computational Physics, Vol. 2
//! - Rys quadrature: Dupuis, Rys & King (1976), J. Chem. Phys. 65, 111
//! - SCF/DIIS: Pulay (1980, 1982), Chem. Phys. Lett. 73, 393
//!
//! # Example
//!
//! ```rust
//! use qc_core::VERSION;
//!
//! // VERSION is set at compile time from Cargo.toml
//! assert!(VERSION.contains('.'));
//! ```

pub mod basis;
pub mod boys;
pub mod constants;
pub mod dft;
pub mod integrals;
pub mod ir;
pub mod optimizer;
pub mod orbital;
pub mod population;
pub mod raman;
pub mod rys;
pub mod scf;
pub mod simd;
pub mod spectra;
pub mod thermo;
pub mod thermochemistry;

#[cfg(test)]
mod golden_phase5;

/// Library version from Cargo.toml
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::const_is_empty)]
    fn version_is_valid() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn modules_are_accessible() {
        // Verify all modules compile and are accessible
        let _ = basis::VERSION;
        let _ = boys::VERSION;
        let _ = dft::VERSION;
        let _ = integrals::VERSION;
        let _ = optimizer::VERSION;
        let _ = orbital::VERSION;
        let _ = population::VERSION;
        let _ = rys::VERSION;
        let _ = scf::VERSION;
    }
}
