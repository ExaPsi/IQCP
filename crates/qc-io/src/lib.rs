//! Shared data schemas for IQCP
//!
//! Defines RunState, RunArtifact, and serialization formats for the
//! Interactive Quantum Chemistry Playground.
//!
//! # Overview
//!
//! This crate provides the data structures used throughout IQCP:
//! - `RunState`: URL-encodable state for deep linking
//! - `RunArtifact`: Exportable computation results with metadata
//! - Module-specific parameter and result types
//!
//! # Example
//!
//! ```rust
//! use qc_io::VERSION;
//! use qc_io::run_state::{RunStateV1, BoysView};
//!
//! // VERSION is set at compile time from Cargo.toml
//! assert!(VERSION.contains('.'));
//!
//! let state = RunStateV1::boys(3, 5.0, BoysView::Single);
//! assert!(state.validate().is_ok());
//! ```

pub mod run_state;

pub use run_state::*;

/// Library version from Cargo.toml
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::const_is_empty)]
    fn version_is_valid() {
        // VERSION is set at compile time from Cargo.toml
        assert!(!VERSION.is_empty());
        // Version should follow semver format
        let parts: Vec<&str> = VERSION.split('.').collect();
        assert!(parts.len() >= 2, "Version should have at least major.minor");
    }
}
