//! Simulated IR and Raman vibrational spectra via Lorentzian or Gaussian
//! broadening.
//!
//! This module implements US-100 (Sprint S41, final story of Milestone M19):
//! closed-form Lorentzian and Gaussian line-shape kernels that convert the
//! stick spectra produced by US-097 (IR intensities in km/mol) and US-098
//! (Raman activities in Å⁴/amu) into continuous broadened curves on a
//! wavenumber grid. The broadening module is shared between IR and Raman —
//! the kernel does not interpret the units of the intensity input and is
//! therefore reusable for any per-mode integrated intensity.
//!
//! # Algorithm
//!
//! A **line-shape function** $f(\tilde\nu; \tilde\nu_k, A_k, \text{width})$
//! describes how a single sharp transition at $\tilde\nu_k$ with integrated
//! intensity $A_k$ contributes to the observed spectrum at wavenumber
//! $\tilde\nu$. The standard normalization is
//!
//! $$\int_{-\infty}^{\infty} f(\tilde\nu;\tilde\nu_k,A_k,\text{width})\,d\tilde\nu = A_k$$
//!
//! so the **integrated intensity** (area under the curve) equals the physical
//! oscillator strength regardless of shape. The total spectrum is the sum
//! over all peaks
//!
//! $$S(\tilde\nu) = \sum_{k: \tilde\nu_k > 0,\, \tilde\nu_k \in [\tilde\nu_\text{min},\tilde\nu_\text{max}]} f(\tilde\nu; \tilde\nu_k, A_k, \text{width})$$
//!
//! with imaginary modes ($\tilde\nu_k < 0$) and out-of-grid peaks filtered
//! out silently. Two kernels are supported:
//!
//! **Lorentzian** (Cauchy-Lorentz, physical origin: damped harmonic
//! oscillator — natural line width via the uncertainty principle):
//!
//! $$L(\tilde\nu) = \frac{A_k}{\pi} \cdot \frac{\gamma}{(\tilde\nu - \tilde\nu_k)^2 + \gamma^2}, \qquad \gamma = \frac{\text{FWHM}}{2}$$
//!
//! **Gaussian** (normal distribution, physical origin: Doppler broadening
//! via the Maxwell-Boltzmann velocity distribution and/or instrumental
//! broadening via a finite-slit spectrometer):
//!
//! $$G(\tilde\nu) = \frac{A_k}{\sigma\sqrt{2\pi}} \exp\!\left(-\frac{(\tilde\nu - \tilde\nu_k)^2}{2\sigma^2}\right), \qquad \sigma = \frac{\text{FWHM}}{2\sqrt{2\ln 2}}$$
//!
//! Both kernels achieve their peak value at $\tilde\nu = \tilde\nu_k$
//! (respectively $A_k/(\pi\gamma)$ and $A_k/(\sigma\sqrt{2\pi})$) and drop
//! to half their peak value at exactly $\tilde\nu_k \pm \text{FWHM}/2$.
//! Both are normalized to the same integrated intensity $A_k$.
//!
//! # Default parameters
//!
//! The `SpectrumParams::default()` values match the IQCP educational
//! defaults: Lorentzian broadening, FWHM = 8 cm⁻¹, wavenumber grid from
//! 0 to 4500 cm⁻¹ at 1 cm⁻¹ resolution (4501 points). This choice
//! matches the typical experimental resolution of a benchtop FT-IR
//! spectrometer or a high-resolution dispersive Raman instrument and is
//! consistent with the default used in Gaussian 16 and ORCA 5 for simulated
//! spectra visualization.
//!
//! The ~8 points per FWHM gives <1% relative error in the trapezoidal
//! integrated-intensity reconstruction, which is the central correctness
//! test for the discretized kernel.
//!
//! # Consumers and producers
//!
//! **Upstream (providers of `frequencies_cm1`):**
//! - [`crate::thermo::FrequencyInfo`] (US-096) — `freq_wavenumber` field
//!
//! **Upstream (providers of `intensities`):**
//! - [`crate::ir::IrSpectrumResult`] (US-097) — `intensities_km_per_mol` field
//! - [`crate::raman::RamanResult`] (US-098) — `raman_activities_a4_amu` field
//!
//! **Downstream consumers:**
//! - `qc-wasm` (US-101, Sprint S42) — serializes [`Spectrum`] over postMessage
//! - Frontend Module C Spectrum tab (US-103, Sprint S42) — Plotly rendering
//!
//! # References
//!
//! - **Griffiths, D.J. (2013).** *Introduction to Electrodynamics* (4th ed.).
//!   Pearson. Chapter 11.2 (Radiation damping) derives the Lorentzian line
//!   shape as the frequency-space response of a damped charged oscillator.
//! - **Herzberg, G. (1945).** *Molecular Spectra and Molecular Structure*
//!   Vol II: Infrared and Raman Spectra of Polyatomic Molecules. Van
//!   Nostrand. Chapter II.B covers vibrational spectra and derives the
//!   Gaussian Doppler width from the Maxwell-Boltzmann velocity distribution.
//! - **Demtröder, W. (2008).** *Laser Spectroscopy, Volume 1: Basic
//!   Principles* (4th ed.). Springer. Chapter 3 is a comprehensive review
//!   of line-shape theory including Lorentzian, Gaussian, and Voigt profiles.
//!
//! # Rationale for closed-form kernels
//!
//! The Lorentzian and Gaussian kernels are evaluated pointwise via their
//! closed-form probability density functions; there is no finite-difference
//! or numerical-integration step at evaluation time. This makes the
//! discretization error entirely attributable to the grid spacing, which
//! is under the user's control. The trapezoidal integrated-intensity
//! reconstruction test at ~8 points per FWHM is the worst-case check for
//! this discretization.

use std::f64::consts::{LN_2, PI};
use thiserror::Error;

/// Module version (matches crate version).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// =============================================================================
// Types
// =============================================================================

/// Broadening kernel for simulated spectra.
///
/// - [`BroadeningKind::Lorentzian`]: Cauchy-Lorentz line shape, arising
///   physically from a damped harmonic oscillator (natural line width).
///   Has slowly-decaying $1/\tilde\nu^2$ tails.
/// - [`BroadeningKind::Gaussian`]: Gauss (normal) line shape, arising
///   physically from Doppler broadening or instrumental (slit-width)
///   broadening. Has exponentially-decaying tails.
///
/// Both kernels are normalized to the same integrated intensity, so the
/// choice only affects the peak shape — not the total area.
///
/// # References
///
/// - Griffiths (2013) *Electrodynamics* 4e, Ch 11.2
/// - Demtröder (2008) *Laser Spectroscopy* Vol 1, Ch 3.1-3.2
/// - Herzberg (1945) *Molecular Spectra* Vol II, Ch II.B
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BroadeningKind {
    /// Cauchy-Lorentz line shape:
    /// `L(ν) = (A/π) · γ / ((ν - ν_k)² + γ²)`, with `γ = FWHM/2`.
    Lorentzian,
    /// Gauss line shape:
    /// `G(ν) = (A / (σ·sqrt(2π))) · exp(-(ν - ν_k)² / (2σ²))`,
    /// with `σ = FWHM / (2·sqrt(2·ln 2))`.
    Gaussian,
}

/// Parameters for spectrum simulation.
///
/// Use [`SpectrumParams::default()`] to get the standard IQCP defaults:
/// Lorentzian broadening, FWHM = 8 cm⁻¹, grid from 0 to 4500 cm⁻¹ at
/// 1 cm⁻¹ resolution (4501 grid points total). These defaults match AC4
/// of US-100 and the conventions used by Gaussian 16 and ORCA 5.
///
/// # Fields
///
/// - `kind`: Broadening kernel shape ([`BroadeningKind::Lorentzian`] or
///   [`BroadeningKind::Gaussian`])
/// - `fwhm_cm1`: Full width at half maximum in cm⁻¹ (not HWHM, not σ).
///   Must be > 0 and finite. Default 8.0.
/// - `grid_min_cm1`: Lower bound of the wavenumber grid in cm⁻¹
///   (inclusive). Must be finite. Default 0.0.
/// - `grid_max_cm1`: Upper bound of the wavenumber grid in cm⁻¹
///   (inclusive). Must be finite and > `grid_min_cm1`. Default 4500.0.
/// - `grid_step_cm1`: Grid spacing in cm⁻¹. Must be > 0 and finite.
///   Default 1.0 gives approximately 8 points per FWHM for the default
///   FWHM, which keeps the trapezoidal integrated-intensity
///   reconstruction within ~1% relative error.
///
/// # Example
///
/// ```rust,ignore
/// use qc_core::spectra::{SpectrumParams, BroadeningKind};
///
/// // Default: Lorentzian, FWHM=8, grid 0-4500 at 1 cm⁻¹
/// let default_params = SpectrumParams::default();
///
/// // Gaussian with FWHM=16, same grid
/// let wide_gauss = SpectrumParams {
///     kind: BroadeningKind::Gaussian,
///     fwhm_cm1: 16.0,
///     ..SpectrumParams::default()
/// };
///
/// // Fingerprint region at high resolution
/// let fingerprint = SpectrumParams {
///     grid_min_cm1: 500.0,
///     grid_max_cm1: 1800.0,
///     grid_step_cm1: 0.25,
///     fwhm_cm1: 4.0,
///     ..SpectrumParams::default()
/// };
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectrumParams {
    pub kind: BroadeningKind,
    pub fwhm_cm1: f64,
    pub grid_min_cm1: f64,
    pub grid_max_cm1: f64,
    pub grid_step_cm1: f64,
}

impl Default for SpectrumParams {
    /// Standard IQCP defaults matching AC4 of US-100: Lorentzian
    /// broadening, FWHM = 8 cm⁻¹, grid from 0 to 4500 cm⁻¹ at 1 cm⁻¹
    /// resolution (4501 points).
    fn default() -> Self {
        Self {
            kind: BroadeningKind::Lorentzian,
            fwhm_cm1: 8.0,
            grid_min_cm1: 0.0,
            grid_max_cm1: 4500.0,
            grid_step_cm1: 1.0,
        }
    }
}

/// Simulated IR or Raman spectrum with broadened peaks on a wavenumber
/// grid.
///
/// Fields are ordered for convenient consumption by visualization layers
/// (Plotly, Three.js): `wavenumbers_cm1` as the x-axis, `intensity` as
/// the y-axis, `kind` and `fwhm_cm1` as plot annotations.
///
/// # Fields
///
/// - `wavenumbers_cm1`: Wavenumber grid in cm⁻¹, monotonically increasing
///   from `params.grid_min_cm1` in `params.grid_step_cm1` increments,
///   with `((grid_max - grid_min) / step).round() + 1` points total.
/// - `intensity`: Broadened intensity at each grid point, same length as
///   `wavenumbers_cm1`. Units depend on input: km/mol for IR,
///   Å⁴/amu for Raman, or arbitrary user-defined units.
/// - `kind`: Broadening kernel used
///   ([`BroadeningKind::Lorentzian`] or [`BroadeningKind::Gaussian`]).
/// - `fwhm_cm1`: FWHM used for the broadening, in cm⁻¹. Copied directly
///   from `SpectrumParams.fwhm_cm1`.
#[derive(Debug, Clone, PartialEq)]
pub struct Spectrum {
    pub wavenumbers_cm1: Vec<f64>,
    pub intensity: Vec<f64>,
    pub kind: BroadeningKind,
    pub fwhm_cm1: f64,
}

/// Errors returned by [`simulate_spectrum`] and its wrappers.
///
/// All errors are structural (invalid input shapes or parameters) and
/// the function never panics on valid data.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum SpectrumError {
    /// Input slices `frequencies_cm1` and `intensities` have different
    /// lengths.
    #[error("mismatched lengths: frequencies = {frequencies}, intensities = {intensities}")]
    MismatchedLengths {
        frequencies: usize,
        intensities: usize,
    },

    /// An input value (frequency, intensity, or a `SpectrumParams` field)
    /// is NaN or infinity.
    #[error("non-finite input: {field} = {value}")]
    NonFiniteInput { field: &'static str, value: f64 },

    /// `fwhm_cm1` is not strictly positive (or is non-finite).
    #[error("FWHM must be > 0 cm⁻¹, got {0}")]
    InvalidFwhm(f64),

    /// `grid_max_cm1` is not strictly greater than `grid_min_cm1`.
    #[error("grid range invalid: grid_max_cm1 ({max}) must be > grid_min_cm1 ({min})")]
    InvalidGridRange { min: f64, max: f64 },

    /// `grid_step_cm1` is not strictly positive (or is non-finite).
    #[error("grid step must be > 0 cm⁻¹, got {0}")]
    InvalidGridStep(f64),
}

// =============================================================================
// Private helpers
// =============================================================================

/// Convert FWHM to Lorentzian HWHM (gamma).
///
/// For the Cauchy-Lorentz distribution the full width at half maximum
/// (FWHM) is exactly twice the half-width at half maximum (HWHM), so
/// $\gamma = \text{FWHM}/2$.
///
/// Reference: Griffiths (2013) *Electrodynamics* 4e, Ch 11.2.
#[inline]
fn fwhm_to_gamma(fwhm: f64) -> f64 {
    fwhm / 2.0
}

/// Convert FWHM to Gaussian standard deviation (sigma).
///
/// For the Gauss (normal) distribution the relationship between FWHM and
/// standard deviation is derived by setting
/// $\exp(-(\tilde\nu - \tilde\nu_k)^2 / (2\sigma^2)) = 1/2$ and solving
/// for $\tilde\nu - \tilde\nu_k$:
///
/// $$\text{FWHM} = 2\sigma\sqrt{2\ln 2}, \qquad \sigma = \frac{\text{FWHM}}{2\sqrt{2\ln 2}}$$
///
/// The factor $2\sqrt{2\ln 2} \approx 2.3548200450309493$ uses
/// `std::f64::consts::LN_2` for bit-identical cross-platform behavior.
///
/// Reference: Herzberg (1945) *Molecular Spectra* Vol II, Ch II.B.
#[inline]
fn fwhm_to_sigma(fwhm: f64) -> f64 {
    fwhm / (2.0 * (2.0 * LN_2).sqrt())
}

/// Lorentzian kernel value at a single grid point.
///
/// Formula:
///
/// $$L(\tilde\nu; \tilde\nu_k, A_k, \gamma) = \frac{A_k}{\pi} \cdot \frac{\gamma}{(\tilde\nu - \tilde\nu_k)^2 + \gamma^2}$$
///
/// At the peak center ($\tilde\nu = \tilde\nu_k$) this equals
/// $A_k/(\pi\gamma)$; at the half-max points ($\tilde\nu = \tilde\nu_k \pm
/// \gamma$) this equals $A_k/(2\pi\gamma)$.
///
/// # Arguments
/// * `nu` — Grid wavenumber
/// * `nu_k` — Peak center wavenumber
/// * `amplitude` — Peak integrated intensity `A_k`
/// * `gamma` — Half-width at half maximum (HWHM); not the FWHM directly.
///   Use [`fwhm_to_gamma`] to convert from FWHM.
///
/// Reference: Griffiths (2013) *Electrodynamics* 4e, Ch 11.2
/// (damped harmonic oscillator response function).
#[inline]
fn lorentzian_value(nu: f64, nu_k: f64, amplitude: f64, gamma: f64) -> f64 {
    let dx = nu - nu_k;
    (amplitude / PI) * gamma / (dx * dx + gamma * gamma)
}

/// Gaussian kernel value at a single grid point.
///
/// Formula:
///
/// $$G(\tilde\nu; \tilde\nu_k, A_k, \sigma) = \frac{A_k}{\sigma\sqrt{2\pi}} \exp\!\left(-\frac{(\tilde\nu - \tilde\nu_k)^2}{2\sigma^2}\right)$$
///
/// At the peak center this equals $A_k/(\sigma\sqrt{2\pi})$; at the
/// half-max points ($\tilde\nu = \tilde\nu_k \pm \sigma\sqrt{2\ln 2}$)
/// this equals half the peak value.
///
/// For $|u| = |\tilde\nu - \tilde\nu_k|/\sigma > 40$, the exponent
/// $-u^2/2$ underflows to $-\infty$ and `exp` returns exactly 0, which
/// is the correct behavior (the Gaussian tail is effectively zero
/// at that distance).
///
/// # Arguments
/// * `nu` — Grid wavenumber
/// * `nu_k` — Peak center wavenumber
/// * `amplitude` — Peak integrated intensity `A_k`
/// * `sigma` — Standard deviation; not the FWHM directly. Use
///   [`fwhm_to_sigma`] to convert from FWHM.
///
/// Reference: Herzberg (1945) *Molecular Spectra* Vol II, Ch II.B
/// (Doppler broadening from Maxwell-Boltzmann velocity distribution).
#[inline]
fn gaussian_value(nu: f64, nu_k: f64, amplitude: f64, sigma: f64) -> f64 {
    let dx = nu - nu_k;
    let norm = amplitude / (sigma * (2.0 * PI).sqrt());
    let exponent = -(dx * dx) / (2.0 * sigma * sigma);
    norm * exponent.exp()
}

/// Build a wavenumber grid from `min` to `max` in `step` increments,
/// inclusive of both endpoints.
///
/// The number of grid points is computed as
/// `((max - min) / step).round() as usize + 1`, which correctly gives
/// 4501 points for the default `(0.0, 4500.0, 1.0)`, 2601 points for
/// `(500.0, 1800.0, 0.5)`, etc. The rounding handles the common case
/// where `(max - min) / step` is mathematically an integer but
/// floating-point arithmetic would give something like `15000.0000...01`.
///
/// Grid point $i$ is at `min + (i as f64) * step` for
/// `i = 0, 1, ..., n_points - 1`. The last point is thus
/// `min + (n_points - 1) * step`, which equals `max` when
/// `(max - min) / step` is integer; otherwise, the last point may be
/// slightly less than `max` (the grid uniformly spans
/// `[min, min + (n_points-1)*step]`).
///
/// # Preconditions
///
/// - `max > min`
/// - `step > 0`
/// - All inputs finite
///
/// These are validated by [`validate_params`] before this helper is
/// called, so this function does not return a `Result`.
fn build_wavenumber_grid(min: f64, max: f64, step: f64) -> Vec<f64> {
    let n_points = ((max - min) / step).round() as usize + 1;
    let mut grid = Vec::with_capacity(n_points);
    for i in 0..n_points {
        grid.push(min + (i as f64) * step);
    }
    grid
}

/// Validate the spectrum parameters: FWHM, grid bounds, grid step,
/// and that all values are finite.
///
/// Validation order (fail-fast):
/// 1. `params.fwhm_cm1 > 0` and finite
/// 2. `params.grid_step_cm1 > 0` and finite
/// 3. `params.grid_min_cm1` and `params.grid_max_cm1` both finite
/// 4. `params.grid_max_cm1 > params.grid_min_cm1`
fn validate_params(params: &SpectrumParams) -> Result<(), SpectrumError> {
    // FWHM: must be finite and > 0. NaN fails both `is_finite` and `> 0.0`.
    if !params.fwhm_cm1.is_finite() || params.fwhm_cm1 <= 0.0 {
        return Err(SpectrumError::InvalidFwhm(params.fwhm_cm1));
    }

    // Grid step: must be finite and > 0.
    if !params.grid_step_cm1.is_finite() || params.grid_step_cm1 <= 0.0 {
        return Err(SpectrumError::InvalidGridStep(params.grid_step_cm1));
    }

    // Grid bounds: must be finite.
    if !params.grid_min_cm1.is_finite() {
        return Err(SpectrumError::NonFiniteInput {
            field: "grid_min_cm1",
            value: params.grid_min_cm1,
        });
    }
    if !params.grid_max_cm1.is_finite() {
        return Err(SpectrumError::NonFiniteInput {
            field: "grid_max_cm1",
            value: params.grid_max_cm1,
        });
    }

    // Grid range: max must be strictly greater than min.
    if params.grid_max_cm1 <= params.grid_min_cm1 {
        return Err(SpectrumError::InvalidGridRange {
            min: params.grid_min_cm1,
            max: params.grid_max_cm1,
        });
    }

    Ok(())
}

// =============================================================================
// Public API
// =============================================================================

/// Simulate a broadened spectrum by summing per-peak kernels on a
/// wavenumber grid.
///
/// This is the core broadening function shared between IR and Raman
/// (AC1 of US-100). Given a list of vibrational frequencies in cm⁻¹ and
/// their associated peak intensities (in arbitrary units), the function
/// returns a continuous spectrum on a wavenumber grid, computed by
/// summing a Lorentzian or Gaussian kernel at each peak center.
///
/// # Arguments
///
/// * `frequencies_cm1` — Vibrational frequencies in cm⁻¹. May include
///   imaginary modes as negative values (filtered out).
/// * `intensities` — Peak intensities, same length as `frequencies_cm1`.
///   Units are not interpreted; IR is typically km/mol, Raman is
///   Å⁴/amu. Negative intensities are allowed (difference spectra).
/// * `params` — Broadening kernel, FWHM, and grid bounds.
///
/// # Returns
///
/// [`Spectrum`] with the wavenumber grid, summed intensity, kernel
/// type, and FWHM. The grid always has
/// `((grid_max - grid_min) / grid_step).round() + 1` points; the
/// intensity vector has the same length.
///
/// # Filtering
///
/// - **Imaginary modes** (`freq < 0`) are silently skipped. This matches
///   the PySCF convention (`freq.real > 0`) and the US-099
///   `vibrational_partition` filter.
/// - **Out-of-range peaks** (`freq < grid_min_cm1` or
///   `freq > grid_max_cm1`) are silently skipped. This matches the
///   matplotlib `xlim` convention (peaks outside the view are simply
///   not displayed).
/// - **Empty input** (zero-length frequency list) returns an all-zero
///   spectrum of the correct grid length (not an error).
///
/// # Errors
///
/// - [`SpectrumError::MismatchedLengths`] if
///   `frequencies_cm1.len() != intensities.len()`
/// - [`SpectrumError::InvalidFwhm`] if `params.fwhm_cm1 <= 0` or NaN
/// - [`SpectrumError::InvalidGridRange`] if
///   `params.grid_max_cm1 <= params.grid_min_cm1`
/// - [`SpectrumError::InvalidGridStep`] if `params.grid_step_cm1 <= 0`
///   or NaN
/// - [`SpectrumError::NonFiniteInput`] if any input (frequency,
///   intensity, grid bound) is NaN or Inf
///
/// # Loop order
///
/// The accumulation loop is outer-peaks, inner-grid: for each peak,
/// compute the kernel values at all grid points and add them to the
/// corresponding `intensity[i]`. This is cache-friendly (grid data
/// stays in cache) and deterministic (fixed accumulation order).
///
/// # References
///
/// - Griffiths (2013) *Electrodynamics* 4e, Ch 11.2
///   (damped harmonic oscillator → Lorentzian)
/// - Demtröder (2008) *Laser Spectroscopy* Vol 1, Ch 3.1-3.2
///   (natural, Doppler, instrumental broadening)
/// - Herzberg (1945) *Molecular Spectra* Vol II, Ch II.B
///   (vibrational spectra)
///
/// # Example
///
/// ```rust,ignore
/// use qc_core::spectra::{simulate_spectrum, SpectrumParams, BroadeningKind};
///
/// let freqs = vec![2043.1, 4488.1, 4790.3]; // H2O / STO-3G frequencies
/// let intensities = vec![11.36, 37.42, 20.01]; // km/mol
///
/// let spec = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default())?;
/// assert_eq!(spec.wavenumbers_cm1.len(), 4501);
/// assert_eq!(spec.intensity.len(), 4501);
/// assert_eq!(spec.kind, BroadeningKind::Lorentzian);
/// assert_eq!(spec.fwhm_cm1, 8.0);
/// # Ok::<(), qc_core::spectra::SpectrumError>(())
/// ```
pub fn simulate_spectrum(
    frequencies_cm1: &[f64],
    intensities: &[f64],
    params: &SpectrumParams,
) -> Result<Spectrum, SpectrumError> {
    // Step 1: Validate the parameter struct (FWHM, grid step, grid bounds).
    validate_params(params)?;

    // Step 2: Check that the two input slices have matching lengths.
    if frequencies_cm1.len() != intensities.len() {
        return Err(SpectrumError::MismatchedLengths {
            frequencies: frequencies_cm1.len(),
            intensities: intensities.len(),
        });
    }

    // Step 3: Check that all input values are finite. We allow negative
    // frequencies (imaginary modes — filtered out below) and negative
    // intensities (difference spectra), but not NaN or infinity.
    for &f in frequencies_cm1 {
        if !f.is_finite() {
            return Err(SpectrumError::NonFiniteInput {
                field: "frequencies_cm1",
                value: f,
            });
        }
    }
    for &a in intensities {
        if !a.is_finite() {
            return Err(SpectrumError::NonFiniteInput {
                field: "intensities",
                value: a,
            });
        }
    }

    // Step 4: Build the wavenumber grid (inclusive of both endpoints).
    let grid = build_wavenumber_grid(
        params.grid_min_cm1,
        params.grid_max_cm1,
        params.grid_step_cm1,
    );
    let n_points = grid.len();

    // Step 5: Convert FWHM to the shape-specific width parameter.
    let width = match params.kind {
        BroadeningKind::Lorentzian => fwhm_to_gamma(params.fwhm_cm1),
        BroadeningKind::Gaussian => fwhm_to_sigma(params.fwhm_cm1),
    };

    // Step 6: Allocate the output intensity vector, initialized to zero.
    let mut intensity = vec![0.0_f64; n_points];

    // Step 7: Accumulate per-peak contributions.
    // Loop order: outer over peaks (few), inner over grid (many).
    // This keeps the peak parameters (nu_k, A_k, width) in registers
    // while sweeping the grid, and keeps the grid cache-hot.
    for (&nu_k, &a_k) in frequencies_cm1.iter().zip(intensities.iter()) {
        // Filter: skip imaginary modes (negative frequencies) silently.
        if nu_k < 0.0 {
            continue;
        }
        // Filter: skip out-of-range peaks silently.
        if nu_k < params.grid_min_cm1 || nu_k > params.grid_max_cm1 {
            continue;
        }

        match params.kind {
            BroadeningKind::Lorentzian => {
                for (i, &nu) in grid.iter().enumerate() {
                    intensity[i] += lorentzian_value(nu, nu_k, a_k, width);
                }
            }
            BroadeningKind::Gaussian => {
                for (i, &nu) in grid.iter().enumerate() {
                    intensity[i] += gaussian_value(nu, nu_k, a_k, width);
                }
            }
        }
    }

    Ok(Spectrum {
        wavenumbers_cm1: grid,
        intensity,
        kind: params.kind,
        fwhm_cm1: params.fwhm_cm1,
    })
}

/// Convenience wrapper: simulate an IR spectrum from IR intensities.
///
/// This is a thin wrapper over [`simulate_spectrum`] using
/// [`SpectrumParams::default()`] but overriding `kind` and `fwhm_cm1`
/// from the function arguments. The resulting `SpectrumParams` uses
/// the IQCP default wavenumber grid `(0.0, 4500.0, 1.0)` (4501 points).
///
/// # Arguments
///
/// * `frequencies_cm1` — Vibrational frequencies in cm⁻¹ (e.g., from
///   [`crate::thermo::FrequencyInfo::freq_wavenumber`])
/// * `ir_intensities_km_per_mol` — IR absorption intensities in km/mol
///   (e.g., from [`crate::ir::IrSpectrumResult::intensities_km_per_mol`])
/// * `fwhm_cm1` — FWHM for the broadening in cm⁻¹ (typical: 8.0)
/// * `kind` — [`BroadeningKind::Lorentzian`] or [`BroadeningKind::Gaussian`]
///
/// # Errors
///
/// Same as [`simulate_spectrum`].
pub fn simulate_ir_spectrum(
    frequencies_cm1: &[f64],
    ir_intensities_km_per_mol: &[f64],
    fwhm_cm1: f64,
    kind: BroadeningKind,
) -> Result<Spectrum, SpectrumError> {
    let params = SpectrumParams {
        kind,
        fwhm_cm1,
        ..SpectrumParams::default()
    };
    simulate_spectrum(frequencies_cm1, ir_intensities_km_per_mol, &params)
}

/// Convenience wrapper: simulate a Raman spectrum from Raman activities.
///
/// This is a thin wrapper over [`simulate_spectrum`] using
/// [`SpectrumParams::default()`] but overriding `kind` and `fwhm_cm1`
/// from the function arguments. The resulting `SpectrumParams` uses
/// the IQCP default wavenumber grid `(0.0, 4500.0, 1.0)` (4501 points).
///
/// # Arguments
///
/// * `frequencies_cm1` — Vibrational frequencies in cm⁻¹
/// * `raman_activities_a4_amu` — Raman scattering activities in Å⁴/amu
///   (e.g., from [`crate::raman::RamanResult::raman_activities_a4_amu`])
/// * `fwhm_cm1` — FWHM for the broadening in cm⁻¹
/// * `kind` — [`BroadeningKind::Lorentzian`] or [`BroadeningKind::Gaussian`]
///
/// # Errors
///
/// Same as [`simulate_spectrum`].
pub fn simulate_raman_spectrum(
    frequencies_cm1: &[f64],
    raman_activities_a4_amu: &[f64],
    fwhm_cm1: f64,
    kind: BroadeningKind,
) -> Result<Spectrum, SpectrumError> {
    let params = SpectrumParams {
        kind,
        fwhm_cm1,
        ..SpectrumParams::default()
    };
    simulate_spectrum(frequencies_cm1, raman_activities_a4_amu, &params)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;

    // -------------------------------------------------------------------------
    // Helper unit tests
    // -------------------------------------------------------------------------

    /// `fwhm_to_gamma(8.0) == 4.0` exactly (HWHM is half of FWHM).
    #[test]
    fn test_fwhm_to_gamma() {
        assert_eq!(fwhm_to_gamma(8.0), 4.0);
        assert_eq!(fwhm_to_gamma(1.0), 0.5);
        assert_eq!(fwhm_to_gamma(16.0), 8.0);
    }

    /// `fwhm_to_sigma(8.0)` matches the analytical value
    /// $8 / (2\sqrt{2\ln 2}) = 8 / 2.3548200450309493 \approx 3.3972872011520763$.
    #[test]
    fn test_fwhm_to_sigma() {
        let sigma = fwhm_to_sigma(8.0);
        let expected = 8.0 / (2.0 * (2.0 * LN_2).sqrt());
        assert!(
            (sigma - expected).abs() < 1e-14,
            "fwhm_to_sigma(8.0) = {sigma}, expected {expected}"
        );
        // Also cross-check against the hardcoded numerical value.
        assert!(
            (sigma - 3.3972872011520763).abs() < 1e-14,
            "fwhm_to_sigma(8.0) = {sigma}, numerical expected ~3.3972872"
        );
    }

    /// `fwhm_to_sigma` inverse check: at `fwhm = 2·sqrt(2·ln 2)`, sigma
    /// is exactly 1 by construction.
    #[test]
    fn test_fwhm_to_sigma_unit() {
        let fwhm = 2.0 * (2.0 * LN_2).sqrt();
        let sigma = fwhm_to_sigma(fwhm);
        assert!(
            (sigma - 1.0).abs() < 1e-15,
            "fwhm_to_sigma({fwhm}) = {sigma}, expected 1.0"
        );
    }

    /// `lorentzian_value` at the peak center equals $A/(\pi\gamma)$.
    #[test]
    fn test_lorentzian_peak_value() {
        let value = lorentzian_value(1600.0, 1600.0, 50.0, 4.0);
        let expected = 50.0 / (PI * 4.0);
        assert!(
            (value - expected).abs() < 1e-12,
            "lorentzian_value at peak = {value}, expected {expected}"
        );
        // Cross-check against hardcoded strategy value.
        assert!((value - 3.978873577297384).abs() < 1e-12);
    }

    /// `lorentzian_value` at $\pm\gamma$ from the peak center is exactly
    /// half the peak value (FWHM recovery).
    #[test]
    fn test_lorentzian_half_max_point() {
        let peak = lorentzian_value(1600.0, 1600.0, 50.0, 4.0);
        let right = lorentzian_value(1604.0, 1600.0, 50.0, 4.0);
        let left = lorentzian_value(1596.0, 1600.0, 50.0, 4.0);
        // Both sides should be exactly peak/2.
        assert!(
            (right - peak / 2.0).abs() < 1e-12,
            "lorentzian at +gamma = {right}, expected {}",
            peak / 2.0
        );
        assert!(
            (left - peak / 2.0).abs() < 1e-12,
            "lorentzian at -gamma = {left}, expected {}",
            peak / 2.0
        );
        // Cross-check against hardcoded strategy value
        // (50 / (π · 8) = 1.989436788648692...).
        assert!((right - 1.989436788648692).abs() < 1e-12);
    }

    /// `lorentzian_value` is symmetric about the peak center.
    #[test]
    fn test_lorentzian_symmetry() {
        let nu_k = 1500.0;
        let a = 10.0;
        let gamma = 3.0;
        for delta in [0.1, 1.0, 5.0, 10.0, 100.0] {
            let right = lorentzian_value(nu_k + delta, nu_k, a, gamma);
            let left = lorentzian_value(nu_k - delta, nu_k, a, gamma);
            assert!(
                (right - left).abs() < 1e-14,
                "Lorentzian asymmetry at delta={delta}: right={right}, left={left}"
            );
        }
    }

    /// `gaussian_value` at the peak center equals $A/(\sigma\sqrt{2\pi})$.
    #[test]
    fn test_gaussian_peak_value() {
        let sigma = fwhm_to_sigma(8.0); // ≈ 3.3972872011520763
        let value = gaussian_value(1600.0, 1600.0, 50.0, sigma);
        let expected = 50.0 / (sigma * (2.0 * PI).sqrt());
        assert!(
            (value - expected).abs() < 1e-12,
            "gaussian_value at peak = {value}, expected {expected}"
        );
        // Cross-check against the numerical value
        // 50 / (3.3972872011520763 · sqrt(2π)) ≈ 5.871482991872821.
        assert!(
            (value - 5.871482991872821).abs() < 1e-10,
            "gaussian_value at peak = {value}, numerical expected ~5.871482"
        );
    }

    /// `gaussian_value` at $\pm\sigma\sqrt{2\ln 2}$ from center is
    /// exactly half the peak value (FWHM recovery).
    #[test]
    fn test_gaussian_half_max_point() {
        let sigma = fwhm_to_sigma(8.0);
        let peak = gaussian_value(1600.0, 1600.0, 50.0, sigma);
        // sigma · sqrt(2·ln 2) = FWHM/2 = 4.0
        let right = gaussian_value(1604.0, 1600.0, 50.0, sigma);
        let left = gaussian_value(1596.0, 1600.0, 50.0, sigma);
        assert!(
            (right - peak / 2.0).abs() < 1e-12,
            "Gaussian at +FWHM/2 = {right}, expected {}",
            peak / 2.0
        );
        assert!(
            (left - peak / 2.0).abs() < 1e-12,
            "Gaussian at -FWHM/2 = {left}, expected {}",
            peak / 2.0
        );
    }

    /// `gaussian_value` is symmetric about the peak center.
    #[test]
    fn test_gaussian_symmetry() {
        let nu_k = 2500.0;
        let a = 7.5;
        let sigma = 2.0;
        for delta in [0.1, 1.0, 3.0, 5.0, 20.0] {
            let right = gaussian_value(nu_k + delta, nu_k, a, sigma);
            let left = gaussian_value(nu_k - delta, nu_k, a, sigma);
            assert!(
                (right - left).abs() < 1e-14,
                "Gaussian asymmetry at delta={delta}: right={right}, left={left}"
            );
        }
    }

    /// Default grid: 0 to 4500 at 1 cm⁻¹ gives 4501 points with correct
    /// endpoint values.
    #[test]
    fn test_build_wavenumber_grid_default() {
        let grid = build_wavenumber_grid(0.0, 4500.0, 1.0);
        assert_eq!(grid.len(), 4501);
        assert_eq!(grid[0], 0.0);
        assert_eq!(grid[4500], 4500.0);
        assert_eq!(grid[2250], 2250.0);
        // Uniform spacing check at several indices.
        for i in [0, 100, 1000, 2000, 4500] {
            assert_eq!(grid[i], i as f64);
        }
    }

    /// Custom grid: (2000, 4000, 0.5) gives 4001 points.
    #[test]
    fn test_build_wavenumber_grid_custom() {
        let grid = build_wavenumber_grid(2000.0, 4000.0, 0.5);
        assert_eq!(grid.len(), 4001);
        assert_eq!(grid[0], 2000.0);
        assert_eq!(grid[4000], 4000.0);
        assert_eq!(grid[2000], 3000.0);
        // Check spacing at i=1
        assert_eq!(grid[1], 2000.5);
    }

    /// Coarser grid with integer step size: (500, 1000, 10) gives 51
    /// points, monotonic increasing.
    #[test]
    fn test_build_wavenumber_grid_step_size() {
        let grid = build_wavenumber_grid(500.0, 1000.0, 10.0);
        assert_eq!(grid.len(), 51);
        assert_eq!(grid[0], 500.0);
        assert_eq!(grid[50], 1000.0);
        for i in 0..50 {
            assert_eq!(grid[i + 1] - grid[i], 10.0);
        }
    }

    /// `SpectrumParams::default()` matches AC4 exactly.
    #[test]
    fn test_spectrum_params_default() {
        let params = SpectrumParams::default();
        assert_eq!(params.kind, BroadeningKind::Lorentzian);
        assert_eq!(params.fwhm_cm1, 8.0);
        assert_eq!(params.grid_min_cm1, 0.0);
        assert_eq!(params.grid_max_cm1, 4500.0);
        assert_eq!(params.grid_step_cm1, 1.0);
    }

    // -------------------------------------------------------------------------
    // Spectrum tests (synthetic + H2O hand inputs)
    // -------------------------------------------------------------------------

    /// Helper: find the grid index of the maximum intensity in a
    /// spectrum (for assertion convenience in tests).
    fn argmax(v: &[f64]) -> usize {
        let mut best = 0usize;
        let mut best_val = f64::NEG_INFINITY;
        for (i, &x) in v.iter().enumerate() {
            if x > best_val {
                best_val = x;
                best = i;
            }
        }
        best
    }

    /// Trapezoidal integration of a spectrum on a uniform grid.
    fn trapezoidal(intensity: &[f64], step: f64) -> f64 {
        if intensity.len() < 2 {
            return 0.0;
        }
        let n = intensity.len();
        let mut sum = 0.5 * (intensity[0] + intensity[n - 1]);
        for x in &intensity[1..n - 1] {
            sum += *x;
        }
        sum * step
    }

    /// Single Lorentzian peak: max at grid[1600] is $50/(\pi\cdot 4)$,
    /// half-max at grid[1604] and grid[1596].
    #[test]
    fn test_simulate_spectrum_single_peak_lorentzian() {
        let freqs = [1600.0];
        let intensities = [50.0];
        let spec = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        assert_eq!(spec.wavenumbers_cm1.len(), 4501);
        assert_eq!(spec.intensity.len(), 4501);
        assert_eq!(spec.kind, BroadeningKind::Lorentzian);
        assert_eq!(spec.fwhm_cm1, 8.0);
        // Max is at grid index 1600 (grid point 1600.0 cm⁻¹).
        let imax = argmax(&spec.intensity);
        assert_eq!(imax, 1600);
        let peak_value = spec.intensity[1600];
        let expected_peak = 50.0 / (PI * 4.0);
        assert!(
            (peak_value - expected_peak).abs() < 1e-12,
            "peak value = {peak_value}, expected {expected_peak}"
        );
        // Half-max at ±FWHM/2 = ±4 cm⁻¹.
        let half = peak_value / 2.0;
        assert!((spec.intensity[1604] - half).abs() < 1e-12);
        assert!((spec.intensity[1596] - half).abs() < 1e-12);
    }

    /// Single Gaussian peak: max at grid[1600] is $50/(\sigma\sqrt{2\pi})$.
    #[test]
    fn test_simulate_spectrum_single_peak_gaussian() {
        let freqs = [1600.0];
        let intensities = [50.0];
        let params = SpectrumParams {
            kind: BroadeningKind::Gaussian,
            ..SpectrumParams::default()
        };
        let spec = simulate_spectrum(&freqs, &intensities, &params).unwrap();
        assert_eq!(spec.kind, BroadeningKind::Gaussian);
        let imax = argmax(&spec.intensity);
        assert_eq!(imax, 1600);
        let peak_value = spec.intensity[1600];
        let sigma = fwhm_to_sigma(8.0);
        let expected_peak = 50.0 / (sigma * (2.0 * PI).sqrt());
        assert!(
            (peak_value - expected_peak).abs() < 1e-12,
            "peak value = {peak_value}, expected {expected_peak}"
        );
        // Half-max at ±FWHM/2 = ±4 cm⁻¹ (FWHM conversion).
        let half = peak_value / 2.0;
        assert!((spec.intensity[1604] - half).abs() < 1e-12);
        assert!((spec.intensity[1596] - half).abs() < 1e-12);
    }

    /// Three synthetic H₂O-like peaks at [1600, 3700, 3850] cm⁻¹
    /// with intensities [60, 10, 50] produce three local maxima at
    /// the expected heights.
    #[test]
    fn test_simulate_spectrum_three_peaks_h2o() {
        let freqs = [1600.0, 3700.0, 3850.0];
        let intensities = [60.0, 10.0, 50.0];
        let spec = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        // Height at each peak center equals A_k / (π · γ).
        let gamma = 4.0;
        let peak_1600 = spec.intensity[1600];
        let peak_3700 = spec.intensity[3700];
        let peak_3850 = spec.intensity[3850];
        // Note: the two 3700/3850 peaks are 150 cm⁻¹ apart, which is
        // ~19x the FWHM, so their tails contribute only ~1.4e-3 to each
        // other — not significant to the first 3 digits. Still, we have
        // to allow for some cross-contamination.
        let expected_1600 = 60.0 / (PI * gamma);
        let expected_3700 = 10.0 / (PI * gamma);
        let expected_3850 = 50.0 / (PI * gamma);
        // At the 1600 cm⁻¹ peak, the other two peaks are >2000 cm⁻¹ away,
        // so contamination is <1e-6 — we can use tight tolerance.
        assert!(
            (peak_1600 - expected_1600).abs() < 1e-4,
            "peak at 1600: got {peak_1600}, expected {expected_1600}"
        );
        // At the 3700/3850 peaks, cross-contamination is larger (~0.0014)
        // but manageable.
        assert!(
            (peak_3700 - expected_3700).abs() < 1e-2,
            "peak at 3700: got {peak_3700}, expected {expected_3700}"
        );
        assert!(
            (peak_3850 - expected_3850).abs() < 1e-2,
            "peak at 3850: got {peak_3850}, expected {expected_3850}"
        );
    }

    /// Integrated intensity of a single Lorentzian peak on the default
    /// grid is within 1% of the input `A_k`.
    #[test]
    fn test_simulate_spectrum_integrated_intensity_lorentzian() {
        let freqs = [2000.0];
        let intensities = [100.0];
        let spec = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        let integrated = trapezoidal(&spec.intensity, 1.0);
        let relative_error = (integrated - 100.0).abs() / 100.0;
        assert!(
            relative_error < 0.01,
            "Lorentzian integrated intensity = {integrated}, expected 100, relative error {relative_error}"
        );
    }

    /// Integrated intensity of a single Gaussian peak on the default
    /// grid is within 0.1% of the input `A_k`.
    #[test]
    fn test_simulate_spectrum_integrated_intensity_gaussian() {
        let freqs = [2000.0];
        let intensities = [100.0];
        let params = SpectrumParams {
            kind: BroadeningKind::Gaussian,
            ..SpectrumParams::default()
        };
        let spec = simulate_spectrum(&freqs, &intensities, &params).unwrap();
        let integrated = trapezoidal(&spec.intensity, 1.0);
        let relative_error = (integrated - 100.0).abs() / 100.0;
        assert!(
            relative_error < 0.001,
            "Gaussian integrated intensity = {integrated}, expected 100, relative error {relative_error}"
        );
    }

    /// Integrated intensity of three H₂O-like Lorentzian peaks is
    /// within 1% of the sum of input intensities.
    #[test]
    fn test_simulate_spectrum_integrated_intensity_three_peaks() {
        let freqs = [1600.0, 3700.0, 3850.0];
        let intensities = [60.0, 10.0, 50.0];
        let spec = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        let integrated = trapezoidal(&spec.intensity, 1.0);
        let expected = 120.0;
        let relative_error = (integrated - expected).abs() / expected;
        assert!(
            relative_error < 0.01,
            "3-peak integrated intensity = {integrated}, expected {expected}, relative error {relative_error}"
        );
    }

    /// Imaginary modes (negative frequency) are skipped — the spectrum
    /// is identical whether or not they are present in the input.
    #[test]
    fn test_simulate_spectrum_imaginary_skipped() {
        // With an imaginary mode included.
        let freqs_with_imag = [100.0, -200.0, 1500.0, 3000.0];
        let intensities_with_imag = [1.0, 1.0, 1.0, 1.0];
        let spec_with_imag = simulate_spectrum(
            &freqs_with_imag,
            &intensities_with_imag,
            &SpectrumParams::default(),
        )
        .unwrap();

        // Without the imaginary mode.
        let freqs = [100.0, 1500.0, 3000.0];
        let intensities = [1.0, 1.0, 1.0];
        let spec = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();

        // The two spectra should be bit-identical.
        assert_eq!(spec_with_imag.intensity.len(), spec.intensity.len());
        for (i, (a, b)) in spec_with_imag
            .intensity
            .iter()
            .zip(spec.intensity.iter())
            .enumerate()
        {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "differ at index {i}: with_imag={a}, without={b}"
            );
        }
    }

    /// Out-of-range peaks (above `grid_max_cm1`) are silently skipped.
    #[test]
    fn test_simulate_spectrum_out_of_range_skipped() {
        // 5000 cm⁻¹ is above the default grid_max = 4500.
        let freqs_with_oor = [100.0, 5000.0, 3000.0];
        let intensities_with_oor = [1.0, 1.0, 1.0];
        let spec_with_oor = simulate_spectrum(
            &freqs_with_oor,
            &intensities_with_oor,
            &SpectrumParams::default(),
        )
        .unwrap();

        let freqs = [100.0, 3000.0];
        let intensities = [1.0, 1.0];
        let spec = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();

        // The two spectra should be bit-identical.
        for (a, b) in spec_with_oor.intensity.iter().zip(spec.intensity.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    /// Out-of-range peaks below `grid_min_cm1` are also silently
    /// skipped (symmetric case).
    #[test]
    fn test_simulate_spectrum_below_grid_min_skipped() {
        // With a custom grid starting at 500, peak at 100 is below range.
        let params = SpectrumParams {
            grid_min_cm1: 500.0,
            grid_max_cm1: 2000.0,
            ..SpectrumParams::default()
        };
        let freqs = [100.0, 1200.0];
        let intensities = [1.0, 1.0];
        let spec = simulate_spectrum(&freqs, &intensities, &params).unwrap();
        // Spectrum should be equivalent to just [1200] (peak at 100 skipped).
        let freqs_single = [1200.0];
        let intensities_single = [1.0];
        let spec_single = simulate_spectrum(&freqs_single, &intensities_single, &params).unwrap();
        for (a, b) in spec.intensity.iter().zip(spec_single.intensity.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    /// Empty input returns an all-zero spectrum of the correct length
    /// (not an error).
    #[test]
    fn test_simulate_spectrum_empty_modes() {
        let spec = simulate_spectrum(&[], &[], &SpectrumParams::default()).unwrap();
        assert_eq!(spec.wavenumbers_cm1.len(), 4501);
        assert_eq!(spec.intensity.len(), 4501);
        assert!(spec.intensity.iter().all(|&x| x == 0.0));
        assert_eq!(spec.kind, BroadeningKind::Lorentzian);
        assert_eq!(spec.fwhm_cm1, 8.0);
    }

    /// Lorentzian tail at 2000 cm⁻¹ away from peak matches analytical
    /// formula `(A/π)·γ/(dx² + γ²)`.
    #[test]
    fn test_simulate_spectrum_baseline_far_from_peak() {
        let freqs = [2000.0];
        let intensities = [50.0];
        let spec = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        let tail_at_0 = spec.intensity[0];
        let gamma = 4.0;
        let dx: f64 = 2000.0;
        let expected = (50.0 / PI) * gamma / (dx * dx + gamma * gamma);
        assert!(
            (tail_at_0 - expected).abs() < 1e-14,
            "tail at 0 = {tail_at_0}, expected {expected}"
        );
        // Should be very small (~1.59e-5).
        assert!(tail_at_0 > 0.0 && tail_at_0 < 1e-4);
    }

    /// Lorentzian vs Gaussian peak value ratio at the same FWHM
    /// matches analytical ratio $1/\sqrt{\pi\ln 2} \approx 0.6777$.
    #[test]
    fn test_simulate_spectrum_lorentzian_vs_gaussian_same_fwhm() {
        let freqs = [1600.0];
        let intensities = [50.0];
        let lor = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        let gauss_params = SpectrumParams {
            kind: BroadeningKind::Gaussian,
            ..SpectrumParams::default()
        };
        let gauss = simulate_spectrum(&freqs, &intensities, &gauss_params).unwrap();
        // Ratio of peak values.
        let ratio = lor.intensity[1600] / gauss.intensity[1600];
        let expected_ratio = 1.0 / (PI * LN_2).sqrt();
        assert!(
            (ratio - expected_ratio).abs() < 1e-12,
            "L/G peak ratio = {ratio}, expected {expected_ratio}"
        );
    }

    /// Custom grid (500, 1800, 0.5) gives 2601 points and handles peaks
    /// correctly.
    #[test]
    fn test_simulate_spectrum_custom_grid() {
        let params = SpectrumParams {
            grid_min_cm1: 500.0,
            grid_max_cm1: 1800.0,
            grid_step_cm1: 0.5,
            ..SpectrumParams::default()
        };
        let freqs = [800.0, 1200.0, 1600.0];
        let intensities = [1.0, 1.0, 1.0];
        let spec = simulate_spectrum(&freqs, &intensities, &params).unwrap();
        // (1800 - 500) / 0.5 + 1 = 2601
        assert_eq!(spec.wavenumbers_cm1.len(), 2601);
        assert_eq!(spec.intensity.len(), 2601);
        // First grid point is 500, last is 1800.
        assert_eq!(spec.wavenumbers_cm1[0], 500.0);
        assert_eq!(spec.wavenumbers_cm1[2600], 1800.0);
        // Peaks should have non-trivial intensity.
        assert!(spec.intensity.iter().any(|&x| x > 0.01));
    }

    /// `Spectrum` struct invariants: lengths match, metadata populated.
    #[test]
    fn test_simulate_spectrum_struct_fields() {
        let freqs = [1000.0];
        let intensities = [1.0];
        let spec = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        assert_eq!(spec.kind, BroadeningKind::Lorentzian);
        assert_eq!(spec.fwhm_cm1, 8.0);
        assert_eq!(spec.wavenumbers_cm1.len(), spec.intensity.len());
        assert_eq!(spec.wavenumbers_cm1.len(), 4501);
    }

    /// Determinism: repeated calls with the same inputs produce
    /// bit-identical output.
    #[test]
    fn test_simulate_spectrum_determinism() {
        let freqs = [1600.0, 3700.0, 3850.0];
        let intensities = [60.0, 10.0, 50.0];
        let a = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        let b = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        for (x, y) in a.intensity.iter().zip(b.intensity.iter()) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
    }

    // -------------------------------------------------------------------------
    // Error handling tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_simulate_spectrum_mismatched_lengths_error() {
        let err =
            simulate_spectrum(&[1000.0], &[1.0, 2.0], &SpectrumParams::default()).unwrap_err();
        assert_eq!(
            err,
            SpectrumError::MismatchedLengths {
                frequencies: 1,
                intensities: 2
            }
        );
    }

    #[test]
    fn test_simulate_spectrum_zero_fwhm_error() {
        let params = SpectrumParams {
            fwhm_cm1: 0.0,
            ..SpectrumParams::default()
        };
        let err = simulate_spectrum(&[1000.0], &[1.0], &params).unwrap_err();
        assert_eq!(err, SpectrumError::InvalidFwhm(0.0));
    }

    #[test]
    fn test_simulate_spectrum_negative_fwhm_error() {
        let params = SpectrumParams {
            fwhm_cm1: -1.0,
            ..SpectrumParams::default()
        };
        let err = simulate_spectrum(&[1000.0], &[1.0], &params).unwrap_err();
        assert_eq!(err, SpectrumError::InvalidFwhm(-1.0));
    }

    #[test]
    fn test_simulate_spectrum_nan_fwhm_error() {
        let params = SpectrumParams {
            fwhm_cm1: f64::NAN,
            ..SpectrumParams::default()
        };
        let err = simulate_spectrum(&[1000.0], &[1.0], &params).unwrap_err();
        // NaN fails the `is_finite` and `> 0` check, mapping to InvalidFwhm.
        match err {
            SpectrumError::InvalidFwhm(v) => assert!(v.is_nan()),
            other => panic!("expected InvalidFwhm(NaN), got {other:?}"),
        }
    }

    #[test]
    fn test_simulate_spectrum_inf_fwhm_error() {
        let params = SpectrumParams {
            fwhm_cm1: f64::INFINITY,
            ..SpectrumParams::default()
        };
        let err = simulate_spectrum(&[1000.0], &[1.0], &params).unwrap_err();
        // Inf fails `is_finite`, mapping to InvalidFwhm.
        assert_eq!(err, SpectrumError::InvalidFwhm(f64::INFINITY));
    }

    #[test]
    fn test_simulate_spectrum_nan_frequency_error() {
        let err = simulate_spectrum(&[f64::NAN], &[1.0], &SpectrumParams::default()).unwrap_err();
        match err {
            SpectrumError::NonFiniteInput { field, value } => {
                assert_eq!(field, "frequencies_cm1");
                assert!(value.is_nan());
            }
            other => panic!("expected NonFiniteInput, got {other:?}"),
        }
    }

    #[test]
    fn test_simulate_spectrum_inf_frequency_error() {
        let err =
            simulate_spectrum(&[f64::INFINITY], &[1.0], &SpectrumParams::default()).unwrap_err();
        assert_eq!(
            err,
            SpectrumError::NonFiniteInput {
                field: "frequencies_cm1",
                value: f64::INFINITY,
            }
        );
    }

    #[test]
    fn test_simulate_spectrum_nan_intensity_error() {
        let err =
            simulate_spectrum(&[1000.0], &[f64::NAN], &SpectrumParams::default()).unwrap_err();
        match err {
            SpectrumError::NonFiniteInput { field, value } => {
                assert_eq!(field, "intensities");
                assert!(value.is_nan());
            }
            other => panic!("expected NonFiniteInput, got {other:?}"),
        }
    }

    #[test]
    fn test_simulate_spectrum_invalid_grid_range_error() {
        let params = SpectrumParams {
            grid_min_cm1: 4500.0,
            grid_max_cm1: 1000.0,
            ..SpectrumParams::default()
        };
        let err = simulate_spectrum(&[1000.0], &[1.0], &params).unwrap_err();
        assert_eq!(
            err,
            SpectrumError::InvalidGridRange {
                min: 4500.0,
                max: 1000.0,
            }
        );
    }

    #[test]
    fn test_simulate_spectrum_equal_grid_bounds_error() {
        let params = SpectrumParams {
            grid_min_cm1: 1000.0,
            grid_max_cm1: 1000.0,
            ..SpectrumParams::default()
        };
        let err = simulate_spectrum(&[1000.0], &[1.0], &params).unwrap_err();
        assert_eq!(
            err,
            SpectrumError::InvalidGridRange {
                min: 1000.0,
                max: 1000.0,
            }
        );
    }

    #[test]
    fn test_simulate_spectrum_zero_grid_step_error() {
        let params = SpectrumParams {
            grid_step_cm1: 0.0,
            ..SpectrumParams::default()
        };
        let err = simulate_spectrum(&[1000.0], &[1.0], &params).unwrap_err();
        assert_eq!(err, SpectrumError::InvalidGridStep(0.0));
    }

    #[test]
    fn test_simulate_spectrum_negative_grid_step_error() {
        let params = SpectrumParams {
            grid_step_cm1: -1.0,
            ..SpectrumParams::default()
        };
        let err = simulate_spectrum(&[1000.0], &[1.0], &params).unwrap_err();
        assert_eq!(err, SpectrumError::InvalidGridStep(-1.0));
    }

    // -------------------------------------------------------------------------
    // Property-based tests (mathematical invariants)
    // -------------------------------------------------------------------------

    /// Lorentzian values are non-negative for any positive amplitude.
    #[test]
    fn test_property_lorentzian_positivity() {
        let cases: &[(f64, f64, f64, f64)] = &[
            (1600.0, 1600.0, 50.0, 4.0),
            (1500.0, 1600.0, 50.0, 4.0),
            (1700.0, 1600.0, 50.0, 4.0),
            (0.0, 1600.0, 50.0, 4.0),
            (5000.0, 1600.0, 50.0, 4.0),
            (500.0, 2000.0, 10.0, 1.0),
            (2500.0, 2000.0, 10.0, 1.0),
            (2000.0, 2000.0, 100.0, 8.0),
            (3700.0, 3850.0, 20.0, 4.0),
            (3850.0, 3700.0, 20.0, 4.0),
        ];
        for &(nu, nu_k, a, gamma) in cases {
            let v = lorentzian_value(nu, nu_k, a, gamma);
            assert!(
                v >= 0.0 && v.is_finite(),
                "L({nu},{nu_k},{a},{gamma}) = {v} not non-negative/finite"
            );
        }
    }

    /// Gaussian values are non-negative for any positive amplitude.
    #[test]
    fn test_property_gaussian_positivity() {
        let cases: &[(f64, f64, f64, f64)] = &[
            (1600.0, 1600.0, 50.0, 3.4),
            (1500.0, 1600.0, 50.0, 3.4),
            (1700.0, 1600.0, 50.0, 3.4),
            (0.0, 1600.0, 50.0, 3.4),
            (5000.0, 1600.0, 50.0, 3.4),
            (500.0, 2000.0, 10.0, 1.0),
            (2500.0, 2000.0, 10.0, 1.0),
            (2000.0, 2000.0, 100.0, 8.0),
            (3700.0, 3850.0, 20.0, 2.0),
            (3850.0, 3700.0, 20.0, 2.0),
        ];
        for &(nu, nu_k, a, sigma) in cases {
            let v = gaussian_value(nu, nu_k, a, sigma);
            assert!(
                v >= 0.0 && v.is_finite(),
                "G({nu},{nu_k},{a},{sigma}) = {v} not non-negative/finite"
            );
        }
    }

    /// Total spectrum intensity is non-negative when all inputs are positive.
    #[test]
    fn test_simulate_spectrum_non_negative_for_positive_intensities() {
        let freqs = [1000.0, 2000.0, 3000.0];
        let intensities = [50.0, 75.0, 25.0];
        let spec_l = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        assert!(spec_l.intensity.iter().all(|&x| x >= 0.0));
        // Same for Gaussian.
        let params_g = SpectrumParams {
            kind: BroadeningKind::Gaussian,
            ..SpectrumParams::default()
        };
        let spec_g = simulate_spectrum(&freqs, &intensities, &params_g).unwrap();
        assert!(spec_g.intensity.iter().all(|&x| x >= 0.0));
    }

    /// Linearity: the spectrum for `[A1 + A2]` at frequency `nu` equals
    /// the sum of the spectrum for `A1` and the spectrum for `A2`
    /// (same frequency).
    #[test]
    fn test_simulate_spectrum_linear_in_intensities() {
        let freqs = [1600.0];
        let a1 = [30.0];
        let a2 = [20.0];
        let a_sum = [50.0];
        let s1 = simulate_spectrum(&freqs, &a1, &SpectrumParams::default()).unwrap();
        let s2 = simulate_spectrum(&freqs, &a2, &SpectrumParams::default()).unwrap();
        let s_sum = simulate_spectrum(&freqs, &a_sum, &SpectrumParams::default()).unwrap();
        for ((&x1, &x2), &xs) in s1
            .intensity
            .iter()
            .zip(s2.intensity.iter())
            .zip(s_sum.intensity.iter())
        {
            let combined = x1 + x2;
            assert!(
                (combined - xs).abs() < 1e-14 * (1.0 + xs.abs()),
                "linearity: x1+x2={combined}, xs={xs}"
            );
        }
    }

    /// Symmetry at N grid points: for a Lorentzian peak at mid-grid,
    /// `intensity[k_peak + delta] == intensity[k_peak - delta]` for any
    /// delta (within grid bounds).
    #[test]
    fn test_property_lorentzian_symmetry_at_n_points() {
        let freqs = [2000.0];
        let intensities = [50.0];
        let spec = simulate_spectrum(&freqs, &intensities, &SpectrumParams::default()).unwrap();
        // Peak at grid index 2000.
        for delta in [1, 5, 10, 50, 100, 500, 1000] {
            let right = spec.intensity[2000 + delta];
            let left = spec.intensity[2000 - delta];
            assert!(
                (right - left).abs() < 1e-14,
                "Lorentzian spec asymmetry at delta={delta}: right={right}, left={left}"
            );
        }
    }

    /// Symmetry at N grid points for Gaussian.
    #[test]
    fn test_property_gaussian_symmetry_at_n_points() {
        let freqs = [2000.0];
        let intensities = [50.0];
        let params = SpectrumParams {
            kind: BroadeningKind::Gaussian,
            ..SpectrumParams::default()
        };
        let spec = simulate_spectrum(&freqs, &intensities, &params).unwrap();
        for delta in [1, 3, 5, 10, 15] {
            let right = spec.intensity[2000 + delta];
            let left = spec.intensity[2000 - delta];
            assert!(
                (right - left).abs() < 1e-14,
                "Gaussian spec asymmetry at delta={delta}: right={right}, left={left}"
            );
        }
    }

    // -------------------------------------------------------------------------
    // Integration tests (US-097/US-098 goldens)
    // -------------------------------------------------------------------------

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct IrGolden {
        name: String,
        freq_wavenumber: Vec<f64>,
        ir_intensities_km_mol: Vec<f64>,
    }

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct RamanGolden {
        name: String,
        freq_wavenumber: Vec<f64>,
        raman_activities_a4_amu: Vec<f64>,
    }

    fn load_ir_golden_h2o() -> IrGolden {
        serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/golden/ir/h2o_sto3g_rhf.json"
        )))
        .expect("parse h2o IR golden")
    }

    fn load_raman_golden_h2o() -> RamanGolden {
        serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/golden/raman/h2o_sto3g_rhf.json"
        )))
        .expect("parse h2o Raman golden")
    }

    /// End-to-end integration test: load the US-097 H2O IR golden,
    /// broaden the stick spectrum through `simulate_ir_spectrum`, and
    /// verify structural properties of the output.
    ///
    /// Note: H2O STO-3G has an O-H stretch at ~4790 cm⁻¹ which is
    /// above the default `grid_max_cm1=4500`. This test uses a custom
    /// grid extending to 5000 cm⁻¹ so all 3 peaks fit, via
    /// `simulate_spectrum` directly.
    #[test]
    fn test_simulate_ir_spectrum_h2o_end_to_end() {
        let golden = load_ir_golden_h2o();
        assert_eq!(golden.name, "h2o");
        assert_eq!(golden.freq_wavenumber.len(), 3);
        assert_eq!(golden.ir_intensities_km_mol.len(), 3);

        // Extend grid_max to 5000 cm⁻¹ to include the H2O STO-3G
        // stretches that exceed the default 4500 cm⁻¹ limit.
        let params = SpectrumParams {
            kind: BroadeningKind::Lorentzian,
            fwhm_cm1: 8.0,
            grid_min_cm1: 0.0,
            grid_max_cm1: 5000.0,
            grid_step_cm1: 1.0,
        };
        let spec = simulate_spectrum(
            &golden.freq_wavenumber,
            &golden.ir_intensities_km_mol,
            &params,
        )
        .unwrap();

        // Grid: (5000 - 0) / 1 + 1 = 5001 points.
        assert_eq!(spec.wavenumbers_cm1.len(), 5001);
        assert_eq!(spec.intensity.len(), 5001);
        assert_eq!(spec.kind, BroadeningKind::Lorentzian);
        assert_eq!(spec.fwhm_cm1, 8.0);

        // The H2O/STO-3G frequencies are approximately 2043, 4488, 4790 cm⁻¹.
        // Verify that intensity at each peak center (rounded to nearest
        // grid point) is close to the analytical Lorentzian peak height
        // A/(π·γ) with γ=4. Account for inter-peak contamination: the
        // 4488 and 4790 peaks are ~302 cm⁻¹ apart (~38× FWHM), so their
        // tails contribute ~5e-4 to each other, negligible to 3 digits.
        let gamma = 4.0;
        for (k, (&nu_k, &a_k)) in golden
            .freq_wavenumber
            .iter()
            .zip(golden.ir_intensities_km_mol.iter())
            .enumerate()
        {
            let idx = nu_k.round() as usize;
            assert!(idx < 5001, "peak {k} index {idx} out of grid");
            let intensity_at_peak = spec.intensity[idx];
            let expected = a_k / (PI * gamma);
            // The peak center is rounded to the nearest cm⁻¹, so we
            // could be up to 0.5 cm⁻¹ off the center. At FWHM=8 the
            // Lorentzian value at 0.5 cm⁻¹ off-center is
            // L · γ² / (0.25 + γ²) ≈ L · 16/16.25 ≈ 0.985 · L.
            // So at most ~1.5% underestimate; allow 3% tolerance.
            let relative_error = (intensity_at_peak - expected).abs() / expected;
            assert!(
                relative_error < 0.03,
                "H2O IR peak {k} at {nu_k} cm⁻¹: intensity={intensity_at_peak}, expected {expected}, relative error {relative_error}"
            );
        }

        // All intensities should be finite and non-negative (all input
        // intensities are positive).
        for &x in &spec.intensity {
            assert!(x.is_finite() && x >= 0.0);
        }

        // Baseline far from the peaks: at grid index 0 (far from all
        // H2O peaks), the intensity should be very small (Lorentzian
        // tail from ~2000 cm⁻¹ away).
        assert!(spec.intensity[0] < 0.01);
    }

    /// Same as above but using the default grid — only the peaks in
    /// range [0, 4500] contribute, so the out-of-range 4790 cm⁻¹ peak
    /// is silently skipped. Verifies the default-grid pathway for
    /// a realistic H2O input from US-097.
    #[test]
    fn test_simulate_ir_spectrum_h2o_default_grid() {
        let golden = load_ir_golden_h2o();
        let spec = simulate_ir_spectrum(
            &golden.freq_wavenumber,
            &golden.ir_intensities_km_mol,
            8.0,
            BroadeningKind::Lorentzian,
        )
        .unwrap();
        assert_eq!(spec.wavenumbers_cm1.len(), 4501);
        // The 4790 peak is filtered because it exceeds grid_max=4500.
        // Only the 2043 and 4488 peaks contribute to the spectrum.
        for &x in &spec.intensity {
            assert!(x.is_finite() && x >= 0.0);
        }
        // Non-trivial intensity at the in-range peaks.
        let gamma = 4.0;
        let idx_2043 = golden.freq_wavenumber[0].round() as usize;
        let expected_2043 = golden.ir_intensities_km_mol[0] / (PI * gamma);
        assert!((spec.intensity[idx_2043] - expected_2043).abs() / expected_2043 < 0.03);
        let idx_4488 = golden.freq_wavenumber[1].round() as usize;
        let expected_4488 = golden.ir_intensities_km_mol[1] / (PI * gamma);
        assert!((spec.intensity[idx_4488] - expected_4488).abs() / expected_4488 < 0.03);
    }

    /// End-to-end integration test: load the US-098 H2O Raman golden,
    /// broaden the stick spectrum through `simulate_spectrum` with an
    /// extended grid to include all 3 H2O STO-3G peaks.
    ///
    /// Note: H2O STO-3G has an O-H stretch at ~4790 cm⁻¹ which is
    /// above the default `grid_max_cm1=4500`. This test uses a custom
    /// grid extending to 5000 cm⁻¹.
    #[test]
    fn test_simulate_raman_spectrum_h2o_end_to_end() {
        let golden = load_raman_golden_h2o();
        assert_eq!(golden.name, "h2o");
        assert_eq!(golden.freq_wavenumber.len(), 3);
        assert_eq!(golden.raman_activities_a4_amu.len(), 3);

        let params = SpectrumParams {
            kind: BroadeningKind::Lorentzian,
            fwhm_cm1: 8.0,
            grid_min_cm1: 0.0,
            grid_max_cm1: 5000.0,
            grid_step_cm1: 1.0,
        };
        let spec = simulate_spectrum(
            &golden.freq_wavenumber,
            &golden.raman_activities_a4_amu,
            &params,
        )
        .unwrap();

        // Grid: 5001 points.
        assert_eq!(spec.wavenumbers_cm1.len(), 5001);
        assert_eq!(spec.intensity.len(), 5001);
        assert_eq!(spec.kind, BroadeningKind::Lorentzian);
        assert_eq!(spec.fwhm_cm1, 8.0);

        // Verify each peak center intensity matches the analytical
        // Lorentzian peak height. All 3 peaks fit in the extended grid.
        let gamma = 4.0;
        for (k, (&nu_k, &a_k)) in golden
            .freq_wavenumber
            .iter()
            .zip(golden.raman_activities_a4_amu.iter())
            .enumerate()
        {
            let idx = nu_k.round() as usize;
            assert!(idx < 5001, "Raman peak {k} index {idx} out of grid");
            let intensity_at_peak = spec.intensity[idx];
            let expected = a_k / (PI * gamma);
            let relative_error = (intensity_at_peak - expected).abs() / expected;
            assert!(
                relative_error < 0.03,
                "H2O Raman peak {k} at {nu_k} cm⁻¹: intensity={intensity_at_peak}, expected {expected}, relative error {relative_error}"
            );
        }

        // All intensities should be finite and non-negative.
        for &x in &spec.intensity {
            assert!(x.is_finite() && x >= 0.0);
        }
    }

    /// Default-grid pathway for Raman: the 4790 cm⁻¹ peak is silently
    /// filtered, leaving only the 2043 and 4488 peaks.
    #[test]
    fn test_simulate_raman_spectrum_h2o_default_grid() {
        let golden = load_raman_golden_h2o();
        let spec = simulate_raman_spectrum(
            &golden.freq_wavenumber,
            &golden.raman_activities_a4_amu,
            8.0,
            BroadeningKind::Lorentzian,
        )
        .unwrap();
        assert_eq!(spec.wavenumbers_cm1.len(), 4501);
        for &x in &spec.intensity {
            assert!(x.is_finite() && x >= 0.0);
        }
    }

    /// Also run the Gaussian kernel through the IR golden to ensure
    /// both kernels work end-to-end.
    #[test]
    fn test_simulate_ir_spectrum_h2o_gaussian() {
        let golden = load_ir_golden_h2o();
        let spec = simulate_ir_spectrum(
            &golden.freq_wavenumber,
            &golden.ir_intensities_km_mol,
            8.0,
            BroadeningKind::Gaussian,
        )
        .unwrap();
        assert_eq!(spec.kind, BroadeningKind::Gaussian);
        assert_eq!(spec.wavenumbers_cm1.len(), 4501);
        for &x in &spec.intensity {
            assert!(x.is_finite() && x >= 0.0);
        }
    }

    // -------------------------------------------------------------------------
    // Convenience wrapper tests
    // -------------------------------------------------------------------------

    /// `simulate_ir_spectrum` delegates correctly to `simulate_spectrum`
    /// and matches the output of an equivalent direct call.
    #[test]
    fn test_simulate_ir_spectrum_wrapper() {
        let freqs = [1500.0, 3000.0];
        let intensities = [10.0, 20.0];
        let direct = simulate_spectrum(
            &freqs,
            &intensities,
            &SpectrumParams {
                kind: BroadeningKind::Lorentzian,
                fwhm_cm1: 8.0,
                ..SpectrumParams::default()
            },
        )
        .unwrap();
        let wrapper =
            simulate_ir_spectrum(&freqs, &intensities, 8.0, BroadeningKind::Lorentzian).unwrap();
        for (a, b) in direct.intensity.iter().zip(wrapper.intensity.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    /// `simulate_raman_spectrum` delegates correctly to `simulate_spectrum`.
    #[test]
    fn test_simulate_raman_spectrum_wrapper() {
        let freqs = [1500.0, 3000.0];
        let intensities = [5.0, 15.0];
        let direct = simulate_spectrum(
            &freqs,
            &intensities,
            &SpectrumParams {
                kind: BroadeningKind::Gaussian,
                fwhm_cm1: 4.0,
                ..SpectrumParams::default()
            },
        )
        .unwrap();
        let wrapper =
            simulate_raman_spectrum(&freqs, &intensities, 4.0, BroadeningKind::Gaussian).unwrap();
        for (a, b) in direct.intensity.iter().zip(wrapper.intensity.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
}
