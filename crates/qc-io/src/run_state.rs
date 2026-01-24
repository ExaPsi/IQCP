//! RunStateV1 - URL-encodable application state
//!
//! This module defines the complete state schema needed to reproduce a user's view
//! of any module in the IQCP application. It includes schema versioning for
//! future migration support.
//!
//! # Overview
//!
//! The state is designed to be:
//! - JSON-serializable with deterministic output via serde
//! - Compressible to under 500 characters for URL embedding
//! - Version-aware for future schema migrations
//!
//! # Example
//!
//! ```rust
//! use qc_io::run_state::{RunStateV1, ModuleId, BoysParams, BoysView, UiSettings, UiMode};
//!
//! let state = RunStateV1::boys(3, 5.0, BoysView::Single);
//! let json = serde_json::to_string(&state).unwrap();
//! let parsed: RunStateV1 = serde_json::from_str(&json).unwrap();
//! assert_eq!(state, parsed);
//! ```

use serde::{Deserialize, Serialize};

/// Module identifier for routing and type discrimination
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ModuleId {
    /// Boys function exploration module
    Boys,
    /// Rys quadrature visualization module
    Rys,
    /// SCF convergence sandbox module
    Scf,
}

/// Boys module view mode
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum BoysView {
    /// Single point evaluation
    Single,
    /// Sweep over T range
    Sweep,
}

/// Boys module parameters
///
/// Controls the visualization of the Boys function F_m(T) with
/// support for single-point evaluation or parameter sweeps.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoysParams {
    /// Order m of the Boys function F_m(T), must be >= 0
    pub m: u32,
    /// Argument T (must be >= 0)
    #[serde(rename = "T")]
    pub t: f64,
    /// View mode: single point evaluation or sweep over T range
    pub view: BoysView,
}

/// Rys quadrature target accuracy
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RysTarget {
    /// Loose accuracy (1e-4)
    #[serde(rename = "1e-4")]
    Loose,
    /// Medium accuracy (1e-6)
    #[serde(rename = "1e-6")]
    Medium,
    /// Tight accuracy (1e-8)
    #[serde(rename = "1e-8")]
    Tight,
}

/// Rys module parameters
///
/// Controls the Rys quadrature visualization showing roots, weights,
/// and accuracy analysis for different orders.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RysParams {
    /// Number of quadrature points (order), typically 1-10
    pub n: u32,
    /// Parameter T for the Boys function moments
    #[serde(rename = "T")]
    pub t: f64,
    /// Target accuracy for order recommendation
    pub target: RysTarget,
}

/// SCF convergence profile
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ConvergenceProfile {
    /// Loose convergence (1e-4 energy threshold)
    Loose,
    /// Medium convergence (1e-6 energy threshold)
    Medium,
    /// Tight convergence (1e-8 energy threshold)
    Tight,
}

/// SCF module parameters
///
/// Controls the Self-Consistent Field computation including
/// convergence settings and DIIS acceleration.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScfParams {
    /// Identifier for the curated molecular system
    pub system_id: String,
    /// Maximum number of SCF iterations
    pub max_iter: u32,
    /// Convergence profile determining energy/density thresholds
    pub conv: ConvergenceProfile,
    /// Enable DIIS (Direct Inversion in the Iterative Subspace) acceleration
    pub diis: bool,
}

/// Preset reference for curated systems
///
/// Optional references to pre-configured systems used for
/// demonstration and educational purposes.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct PresetRef {
    /// SCF curated system ID (e.g., "h2_sto3g_r1.4")
    #[serde(skip_serializing_if = "Option::is_none")]
    pub system_id: Option<String>,
    /// Molecule ID for display in Boys/Rys modules
    #[serde(skip_serializing_if = "Option::is_none")]
    pub molecule_id: Option<String>,
    /// Basis set ID for display in Boys/Rys modules
    #[serde(skip_serializing_if = "Option::is_none")]
    pub basis_id: Option<String>,
}

/// UI display mode
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum UiMode {
    /// Educational explanations mode
    Explain,
    /// Raw computational internals mode
    Internals,
}

/// UI color theme
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum UiTheme {
    /// Light color theme
    Light,
    /// Dark color theme
    Dark,
}

/// UI settings
///
/// Controls the display mode and visual theme of the application.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UiSettings {
    /// Display mode: educational explanations or raw computational internals
    pub mode: UiMode,
    /// Color theme (optional, defaults to system preference)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub theme: Option<UiTheme>,
}

/// RunStateV1 - URL-encodable application state
///
/// This is the complete state needed to reproduce a user's view
/// of any module in the application. It includes schema versioning
/// for future migration support.
///
/// # JSON Format
///
/// The struct serializes to JSON with the following structure:
///
/// ```json
/// {
///   "schema_version": "run_state_v1",
///   "app_version": "0.1.0",
///   "module": "boys",
///   "boys": { "m": 3, "T": 5.0, "view": "single" },
///   "ui": { "mode": "explain" }
/// }
/// ```
///
/// # Example
///
/// ```rust
/// use qc_io::run_state::{RunStateV1, BoysView};
///
/// let state = RunStateV1::boys(3, 5.0, BoysView::Single);
/// assert_eq!(state.module, qc_io::run_state::ModuleId::Boys);
/// assert!(state.boys.is_some());
/// assert!(state.validate().is_ok());
/// ```
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RunStateV1 {
    /// Schema version for migration support (always "run_state_v1" for this version)
    pub schema_version: String,
    /// Application version that created this state
    pub app_version: String,
    /// Active module determining which view is displayed
    pub module: ModuleId,
    /// Preset references for curated system configurations
    #[serde(skip_serializing_if = "Option::is_none")]
    pub preset: Option<PresetRef>,
    /// Boys module parameters (present when module === Boys)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub boys: Option<BoysParams>,
    /// Rys module parameters (present when module === Rys)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rys: Option<RysParams>,
    /// SCF module parameters (present when module === Scf)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scf: Option<ScfParams>,
    /// UI settings controlling display mode and theme
    pub ui: UiSettings,
}

impl RunStateV1 {
    /// Schema version constant
    pub const SCHEMA_VERSION: &'static str = "run_state_v1";

    /// Create a new RunStateV1 for the Boys module
    ///
    /// # Arguments
    ///
    /// * `m` - Order of the Boys function
    /// * `t` - Argument T (must be >= 0)
    /// * `view` - View mode (single or sweep)
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_io::run_state::{RunStateV1, BoysView};
    ///
    /// let state = RunStateV1::boys(3, 5.0, BoysView::Single);
    /// assert_eq!(state.boys.unwrap().m, 3);
    /// ```
    pub fn boys(m: u32, t: f64, view: BoysView) -> Self {
        Self {
            schema_version: Self::SCHEMA_VERSION.to_string(),
            app_version: crate::VERSION.to_string(),
            module: ModuleId::Boys,
            preset: None,
            boys: Some(BoysParams { m, t, view }),
            rys: None,
            scf: None,
            ui: UiSettings {
                mode: UiMode::Explain,
                theme: None,
            },
        }
    }

    /// Create a new RunStateV1 for the Rys module
    ///
    /// # Arguments
    ///
    /// * `n` - Number of quadrature points
    /// * `t` - Parameter T for Boys function moments
    /// * `target` - Target accuracy level
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_io::run_state::{RunStateV1, RysTarget};
    ///
    /// let state = RunStateV1::rys(3, 1.0, RysTarget::Medium);
    /// assert_eq!(state.rys.unwrap().n, 3);
    /// ```
    pub fn rys(n: u32, t: f64, target: RysTarget) -> Self {
        Self {
            schema_version: Self::SCHEMA_VERSION.to_string(),
            app_version: crate::VERSION.to_string(),
            module: ModuleId::Rys,
            preset: None,
            boys: None,
            rys: Some(RysParams { n, t, target }),
            scf: None,
            ui: UiSettings {
                mode: UiMode::Explain,
                theme: None,
            },
        }
    }

    /// Create a new RunStateV1 for the SCF module
    ///
    /// # Arguments
    ///
    /// * `system_id` - Identifier for the curated molecular system
    /// * `max_iter` - Maximum number of SCF iterations
    /// * `conv` - Convergence profile
    /// * `diis` - Enable DIIS acceleration
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_io::run_state::{RunStateV1, ConvergenceProfile};
    ///
    /// let state = RunStateV1::scf("h2_sto3g_r1.4", 50, ConvergenceProfile::Medium, true);
    /// assert!(state.scf.unwrap().diis);
    /// ```
    pub fn scf(system_id: &str, max_iter: u32, conv: ConvergenceProfile, diis: bool) -> Self {
        Self {
            schema_version: Self::SCHEMA_VERSION.to_string(),
            app_version: crate::VERSION.to_string(),
            module: ModuleId::Scf,
            preset: None,
            boys: None,
            rys: None,
            scf: Some(ScfParams {
                system_id: system_id.to_string(),
                max_iter,
                conv,
                diis,
            }),
            ui: UiSettings {
                mode: UiMode::Explain,
                theme: None,
            },
        }
    }

    /// Validate that the state is internally consistent
    ///
    /// Checks that the active module has its corresponding parameters defined.
    ///
    /// # Returns
    ///
    /// * `Ok(())` if the state is valid
    /// * `Err(&'static str)` with error message if invalid
    ///
    /// # Example
    ///
    /// ```rust
    /// use qc_io::run_state::{RunStateV1, BoysView};
    ///
    /// let mut state = RunStateV1::boys(0, 1.0, BoysView::Single);
    /// assert!(state.validate().is_ok());
    ///
    /// state.boys = None;
    /// assert!(state.validate().is_err());
    /// ```
    pub fn validate(&self) -> Result<(), &'static str> {
        if self.schema_version != Self::SCHEMA_VERSION {
            return Err("Invalid schema version");
        }

        match self.module {
            ModuleId::Boys if self.boys.is_none() => Err("Boys module requires boys parameters"),
            ModuleId::Rys if self.rys.is_none() => Err("Rys module requires rys parameters"),
            ModuleId::Scf if self.scf.is_none() => Err("SCF module requires scf parameters"),
            _ => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boys_serialization() {
        let state = RunStateV1::boys(3, 5.0, BoysView::Single);
        let json = serde_json::to_string(&state).unwrap();

        // Verify JSON contains expected fields
        assert!(json.contains("\"schema_version\":\"run_state_v1\""));
        assert!(json.contains("\"module\":\"boys\""));
        assert!(json.contains("\"m\":3"));
        assert!(json.contains("\"T\":5"));
        assert!(json.contains("\"view\":\"single\""));

        // Verify round-trip
        let parsed: RunStateV1 = serde_json::from_str(&json).unwrap();
        assert_eq!(state, parsed);
    }

    #[test]
    fn test_rys_serialization() {
        let state = RunStateV1::rys(5, 2.5, RysTarget::Tight);
        let json = serde_json::to_string(&state).unwrap();

        assert!(json.contains("\"module\":\"rys\""));
        assert!(json.contains("\"n\":5"));
        assert!(json.contains("\"target\":\"1e-8\""));

        let parsed: RunStateV1 = serde_json::from_str(&json).unwrap();
        assert_eq!(state, parsed);
    }

    #[test]
    fn test_scf_serialization() {
        let state = RunStateV1::scf("h2_sto3g_r1.4", 100, ConvergenceProfile::Loose, false);
        let json = serde_json::to_string(&state).unwrap();

        assert!(json.contains("\"module\":\"scf\""));
        assert!(json.contains("\"system_id\":\"h2_sto3g_r1.4\""));
        assert!(json.contains("\"max_iter\":100"));
        assert!(json.contains("\"conv\":\"loose\""));
        assert!(json.contains("\"diis\":false"));

        let parsed: RunStateV1 = serde_json::from_str(&json).unwrap();
        assert_eq!(state, parsed);
    }

    #[test]
    fn test_validation_boys() {
        let mut state = RunStateV1::boys(0, 1.0, BoysView::Single);
        assert!(state.validate().is_ok());

        state.boys = None;
        assert!(state.validate().is_err());
    }

    #[test]
    fn test_validation_rys() {
        let mut state = RunStateV1::rys(3, 1.0, RysTarget::Medium);
        assert!(state.validate().is_ok());

        state.rys = None;
        assert!(state.validate().is_err());
    }

    #[test]
    fn test_validation_scf() {
        let mut state = RunStateV1::scf("test", 50, ConvergenceProfile::Medium, true);
        assert!(state.validate().is_ok());

        state.scf = None;
        assert!(state.validate().is_err());
    }

    #[test]
    fn test_validation_schema_version() {
        let mut state = RunStateV1::boys(0, 1.0, BoysView::Single);
        state.schema_version = "invalid".to_string();
        assert!(state.validate().is_err());
    }

    #[test]
    fn test_skip_serializing_none_fields() {
        let state = RunStateV1::boys(0, 1.0, BoysView::Single);
        let json = serde_json::to_string(&state).unwrap();

        // Should not contain optional None fields
        assert!(!json.contains("\"preset\""));
        assert!(!json.contains("\"rys\""));
        assert!(!json.contains("\"scf\""));
        assert!(!json.contains("\"theme\""));
    }

    #[test]
    fn test_ui_settings_with_theme() {
        let mut state = RunStateV1::boys(0, 1.0, BoysView::Single);
        state.ui.theme = Some(UiTheme::Dark);

        let json = serde_json::to_string(&state).unwrap();
        assert!(json.contains("\"theme\":\"dark\""));

        let parsed: RunStateV1 = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.ui.theme, Some(UiTheme::Dark));
    }

    #[test]
    fn test_preset_ref_serialization() {
        let mut state = RunStateV1::scf("h2_sto3g_r1.4", 50, ConvergenceProfile::Medium, true);
        state.preset = Some(PresetRef {
            system_id: Some("h2_sto3g_r1.4".to_string()),
            molecule_id: None,
            basis_id: Some("sto-3g".to_string()),
        });

        let json = serde_json::to_string(&state).unwrap();
        assert!(json.contains("\"preset\""));
        assert!(json.contains("\"system_id\":\"h2_sto3g_r1.4\""));
        assert!(json.contains("\"basis_id\":\"sto-3g\""));
        // molecule_id should be skipped
        assert!(!json.contains("\"molecule_id\""));

        let parsed: RunStateV1 = serde_json::from_str(&json).unwrap();
        assert_eq!(state.preset, parsed.preset);
    }

    #[test]
    fn test_json_size_estimate() {
        // Verify typical states are reasonably compact
        let boys_state = RunStateV1::boys(3, 5.0, BoysView::Single);
        let boys_json = serde_json::to_string(&boys_state).unwrap();
        assert!(
            boys_json.len() < 200,
            "Boys JSON too large: {} bytes",
            boys_json.len()
        );

        let rys_state = RunStateV1::rys(5, 10.0, RysTarget::Medium);
        let rys_json = serde_json::to_string(&rys_state).unwrap();
        assert!(
            rys_json.len() < 200,
            "Rys JSON too large: {} bytes",
            rys_json.len()
        );

        let scf_state = RunStateV1::scf("h2_sto3g_r1.4", 50, ConvergenceProfile::Medium, true);
        let scf_json = serde_json::to_string(&scf_state).unwrap();
        assert!(
            scf_json.len() < 300,
            "SCF JSON too large: {} bytes",
            scf_json.len()
        );
    }
}
