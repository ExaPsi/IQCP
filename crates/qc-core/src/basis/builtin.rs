// Allow excessive precision for basis set coefficients from reference implementations.
// These values are taken directly from PySCF 2.11.0 for numerical accuracy.
#![allow(clippy::excessive_precision)]

//! Built-in basis set data
//!
//! This module provides pre-defined basis set parameters for elements H through Ne.
//! The data is sourced from PySCF 2.11.0 and uses pre-normalized coefficients.
//!
//! # Supported Basis Sets
//!
//! | Basis | Description | Elements |
//! |-------|-------------|----------|
//! | STO-3G | Minimal basis (3 Gaussians fit to Slater) | H-Ne |
//! | 3-21G | Split-valence | H-Ne |
//! | 6-31G | Split-valence | H-Ne |
//! | 6-31G* | Split-valence + d polarization | H-Ne |
//! | 6-31+G* | Split-valence + diffuse sp + d polarization | H-Ne |
//!
//! # Data Format
//!
//! For each element and basis set, the data is returned as a vector of shells,
//! where each shell is a tuple of (AngularMomentum, Vec<(exponent, coefficient)>).
//!
//! # References
//!
//! - PySCF basis set library: `references/pyscf/pyscf/gto/basis/`
//! - Basis Set Exchange: https://www.basissetexchange.org/

use super::primitives::AngularMomentum;
use thiserror::Error;

// =============================================================================
// Error Types
// =============================================================================

/// Errors that can occur when retrieving basis set data
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum BasisError {
    /// Unknown basis set name
    #[error("Unknown basis set: '{0}' (supported: sto-3g, 3-21g, 6-31g, 6-31g*, 6-31+g*)")]
    UnknownBasis(String),

    /// Element not available in basis set
    #[error("Element Z={0} not available in basis '{1}'")]
    ElementNotInBasis(u8, String),

    /// Element not supported by IQCP
    #[error("Unsupported element Z={0} (IQCP supports H-Ne, Z=1-10)")]
    UnsupportedElement(u8),

    /// Invalid molecular data
    #[error("Invalid molecular data: {0}")]
    InvalidMolecule(String),
}

/// Shell data type: (angular_momentum, [(exponent, coefficient), ...])
pub type ShellData = (AngularMomentum, Vec<(f64, f64)>);

// =============================================================================
// Main API
// =============================================================================

/// Get basis set data for an element
///
/// Returns a vector of shells for the specified element and basis set.
/// Each shell is a tuple of (AngularMomentum, Vec<(exponent, coefficient)>).
///
/// # Arguments
///
/// * `z` - Atomic number (1-10)
/// * `basis_name` - Basis set name (case insensitive): "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*"
///
/// # Returns
///
/// - `Ok(Vec<ShellData>)` - Vector of shells for the element
/// - `Err(BasisError)` - If basis or element is not available
///
/// # Example
///
/// ```rust
/// use qc_core::basis::{get_element_basis, AngularMomentum};
///
/// // Get hydrogen STO-3G basis
/// let h_basis = get_element_basis(1, "sto-3g").unwrap();
/// assert_eq!(h_basis.len(), 1);  // H has 1 shell (1s)
/// assert_eq!(h_basis[0].0, AngularMomentum::S);
/// assert_eq!(h_basis[0].1.len(), 3);  // 3 primitives
///
/// // Get oxygen STO-3G basis
/// let o_basis = get_element_basis(8, "sto-3g").unwrap();
/// assert_eq!(o_basis.len(), 3);  // O has 3 shells (1s, 2s, 2p)
/// ```
pub fn get_element_basis(z: u8, basis_name: &str) -> Result<Vec<ShellData>, BasisError> {
    // Validate element
    if !(1..=10).contains(&z) {
        return Err(BasisError::UnsupportedElement(z));
    }

    // Normalize basis name
    let basis_lower = basis_name.to_lowercase().replace(" ", "").replace("_", "");

    match basis_lower.as_str() {
        "sto-3g" | "sto3g" => get_sto3g(z),
        "3-21g" | "321g" => get_321g(z),
        "6-31g" | "631g" => get_631g(z),
        "6-31g*" | "631g*" | "6-31gs" | "631gs" | "6-31g(d)" | "631g(d)" => get_631gs(z),
        "6-31+g*" | "631+g*" | "6-31+gs" | "631+gs" | "6-31+g(d)" | "631+g(d)" => get_631pgs(z),
        _ => Err(BasisError::UnknownBasis(basis_name.to_string())),
    }
}

/// Get a list of supported basis set names
pub fn supported_basis_sets() -> &'static [&'static str] {
    &["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*"]
}

/// Check if a basis set is supported
pub fn is_supported_basis(basis_name: &str) -> bool {
    let basis_lower = basis_name.to_lowercase().replace(" ", "").replace("_", "");
    matches!(
        basis_lower.as_str(),
        "sto-3g"
            | "sto3g"
            | "3-21g"
            | "321g"
            | "6-31g"
            | "631g"
            | "6-31g*"
            | "631g*"
            | "6-31gs"
            | "631gs"
            | "6-31g(d)"
            | "631g(d)"
            | "6-31+g*"
            | "631+g*"
            | "6-31+gs"
            | "631+gs"
            | "6-31+g(d)"
            | "631+g(d)"
    )
}

// =============================================================================
// STO-3G Basis Set Data
// Reference: PySCF 2.11.0 gto/basis/sto-3g.dat
// =============================================================================

fn get_sto3g(z: u8) -> Result<Vec<ShellData>, BasisError> {
    match z {
        1 => Ok(sto3g_h()),
        2 => Ok(sto3g_he()),
        3 => Ok(sto3g_li()),
        4 => Ok(sto3g_be()),
        5 => Ok(sto3g_b()),
        6 => Ok(sto3g_c()),
        7 => Ok(sto3g_n()),
        8 => Ok(sto3g_o()),
        9 => Ok(sto3g_f()),
        10 => Ok(sto3g_ne()),
        _ => Err(BasisError::ElementNotInBasis(z, "sto-3g".to_string())),
    }
}

/// Hydrogen STO-3G basis
fn sto3g_h() -> Vec<ShellData> {
    vec![(
        AngularMomentum::S,
        vec![
            (3.4252509100, 0.1543289707),
            (0.6239137300, 0.5353281424),
            (0.1688554000, 0.4446345420),
        ],
    )]
}

/// Helium STO-3G basis
/// Coefficients matched to PySCF 2.11.0 for consistency
fn sto3g_he() -> Vec<ShellData> {
    vec![(
        AngularMomentum::S,
        vec![
            (6.36242139, 0.15432897),
            (1.15892300, 0.53532814),
            (0.31364979, 0.44463454),
        ],
    )]
}

/// Lithium STO-3G basis
fn sto3g_li() -> Vec<ShellData> {
    vec![
        // 1s shell
        (
            AngularMomentum::S,
            vec![
                (16.1195750000, 0.1543289703),
                (2.9362007000, 0.5353281412),
                (0.7946505000, 0.4446345410),
            ],
        ),
        // 2s shell (valence)
        (
            AngularMomentum::S,
            vec![
                (0.6362897000, -0.0999672298),
                (0.1478601000, 0.3995128292),
                (0.0480887000, 0.7001154686),
            ],
        ),
        // 2p shell (valence)
        (
            AngularMomentum::P,
            vec![
                (0.6362897000, 0.1559162650),
                (0.1478601000, 0.6076837007),
                (0.0480887000, 0.3919573775),
            ],
        ),
    ]
}

/// Beryllium STO-3G basis
fn sto3g_be() -> Vec<ShellData> {
    vec![
        // 1s shell
        (
            AngularMomentum::S,
            vec![
                (30.1678710000, 0.1543289700),
                (5.4951153000, 0.5353281400),
                (1.4871927000, 0.4446345400),
            ],
        ),
        // 2s shell (valence)
        (
            AngularMomentum::S,
            vec![
                (1.3148331000, -0.0999672315),
                (0.3055389000, 0.3995128361),
                (0.0993707000, 0.7001154807),
            ],
        ),
        // 2p shell (valence)
        (
            AngularMomentum::P,
            vec![
                (1.3148331000, 0.1559162776),
                (0.3055389000, 0.6076837498),
                (0.0993707000, 0.3919574092),
            ],
        ),
    ]
}

/// Boron STO-3G basis
fn sto3g_b() -> Vec<ShellData> {
    vec![
        // 1s shell
        (
            AngularMomentum::S,
            vec![
                (48.7911130000, 0.1543289687),
                (8.8873622000, 0.5353281356),
                (2.4052670000, 0.4446345363),
            ],
        ),
        // 2s shell
        (
            AngularMomentum::S,
            vec![
                (2.2369561000, -0.0999672287),
                (0.5198205000, 0.3995128246),
                (0.1690618000, 0.7001154606),
            ],
        ),
        // 2p shell
        (
            AngularMomentum::P,
            vec![
                (2.2369561000, 0.1559162685),
                (0.5198205000, 0.6076837141),
                (0.1690618000, 0.3919573862),
            ],
        ),
    ]
}

/// Carbon STO-3G basis
fn sto3g_c() -> Vec<ShellData> {
    vec![
        // 1s shell
        (
            AngularMomentum::S,
            vec![
                (71.6168370000, 0.1543289700),
                (13.0450960000, 0.5353281400),
                (3.5305122000, 0.4446345400),
            ],
        ),
        // 2s shell
        (
            AngularMomentum::S,
            vec![
                (2.9412494000, -0.0999672301),
                (0.6834831000, 0.3995128303),
                (0.2222899000, 0.7001154705),
            ],
        ),
        // 2p shell
        (
            AngularMomentum::P,
            vec![
                (2.9412494000, 0.1559162721),
                (0.6834831000, 0.6076837282),
                (0.2222899000, 0.3919573953),
            ],
        ),
    ]
}

/// Nitrogen STO-3G basis
fn sto3g_n() -> Vec<ShellData> {
    vec![
        // 1s shell
        (
            AngularMomentum::S,
            vec![
                (99.1061690000, 0.1543289687),
                (18.0523120000, 0.5353281356),
                (4.8856602000, 0.4446345363),
            ],
        ),
        // 2s shell
        (
            AngularMomentum::S,
            vec![
                (3.7804559000, -0.0999672287),
                (0.8784966000, 0.3995128246),
                (0.2857143600, 0.7001154606),
            ],
        ),
        // 2p shell
        (
            AngularMomentum::P,
            vec![
                (3.7804559000, 0.1559162685),
                (0.8784966000, 0.6076837141),
                (0.2857143600, 0.3919573862),
            ],
        ),
    ]
}

/// Oxygen STO-3G basis
fn sto3g_o() -> Vec<ShellData> {
    vec![
        // 1s shell
        (
            AngularMomentum::S,
            vec![
                (130.7093200000, 0.1543289687),
                (23.8088610000, 0.5353281356),
                (6.4436083000, 0.4446345363),
            ],
        ),
        // 2s shell
        (
            AngularMomentum::S,
            vec![
                (5.0331513000, -0.0999672287),
                (1.1695961000, 0.3995128246),
                (0.3803890000, 0.7001154606),
            ],
        ),
        // 2p shell
        (
            AngularMomentum::P,
            vec![
                (5.0331513000, 0.1559162685),
                (1.1695961000, 0.6076837141),
                (0.3803890000, 0.3919573862),
            ],
        ),
    ]
}

/// Fluorine STO-3G basis
fn sto3g_f() -> Vec<ShellData> {
    vec![
        // 1s shell
        (
            AngularMomentum::S,
            vec![
                (166.6791300000, 0.1543289687),
                (30.3608120000, 0.5353281356),
                (8.2168207000, 0.4446345363),
            ],
        ),
        // 2s shell
        (
            AngularMomentum::S,
            vec![
                (6.4648032000, -0.0999672287),
                (1.5022812000, 0.3995128246),
                (0.4885885000, 0.7001154606),
            ],
        ),
        // 2p shell
        (
            AngularMomentum::P,
            vec![
                (6.4648032000, 0.1559162685),
                (1.5022812000, 0.6076837141),
                (0.4885885000, 0.3919573862),
            ],
        ),
    ]
}

/// Neon STO-3G basis
fn sto3g_ne() -> Vec<ShellData> {
    vec![
        // 1s shell
        (
            AngularMomentum::S,
            vec![
                (207.0156100000, 0.1543289687),
                (37.7081510000, 0.5353281356),
                (10.2052970000, 0.4446345363),
            ],
        ),
        // 2s shell
        (
            AngularMomentum::S,
            vec![
                (8.2463151000, -0.0999672287),
                (1.9162662000, 0.3995128246),
                (0.6232293000, 0.7001154606),
            ],
        ),
        // 2p shell
        (
            AngularMomentum::P,
            vec![
                (8.2463151000, 0.1559162685),
                (1.9162662000, 0.6076837141),
                (0.6232293000, 0.3919573862),
            ],
        ),
    ]
}

// =============================================================================
// 3-21G Basis Set Data
// Reference: PySCF 2.11.0 gto/basis/3-21g.dat
// =============================================================================

fn get_321g(z: u8) -> Result<Vec<ShellData>, BasisError> {
    match z {
        1 => Ok(b321g_h()),
        2 => Ok(b321g_he()),
        3 => Ok(b321g_li()),
        4 => Ok(b321g_be()),
        5 => Ok(b321g_b()),
        6 => Ok(b321g_c()),
        7 => Ok(b321g_n()),
        8 => Ok(b321g_o()),
        9 => Ok(b321g_f()),
        10 => Ok(b321g_ne()),
        _ => Err(BasisError::ElementNotInBasis(z, "3-21g".to_string())),
    }
}

fn b321g_h() -> Vec<ShellData> {
    vec![
        // Inner s
        (
            AngularMomentum::S,
            vec![(5.4471780000, 0.1562849787), (0.8245472400, 0.9046908767)],
        ),
        // Outer s
        (AngularMomentum::S, vec![(0.1831915800, 1.0000000000)]),
    ]
}

fn b321g_he() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![(13.6267000000, 0.1752300), (1.9993500000, 0.8934830)],
        ),
        (AngularMomentum::S, vec![(0.3829930000, 1.0000000)]),
    ]
}

fn b321g_li() -> Vec<ShellData> {
    vec![
        // 1s
        (
            AngularMomentum::S,
            vec![
                (36.8382000000, 0.0696686),
                (5.4817200000, 0.3813460),
                (1.1132700000, 0.6817020),
            ],
        ),
        // 2s inner
        (
            AngularMomentum::S,
            vec![(0.5402050000, -0.2631270), (0.1022550000, 1.1433900)],
        ),
        // 2s outer
        (AngularMomentum::S, vec![(0.0285650000, 1.0000000)]),
    ]
}

fn b321g_be() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (71.8876000000, 0.0644263),
                (10.7289000000, 0.3660960),
                (2.2220500000, 0.6959340),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.2954800000, -0.4210640), (0.2688810000, 1.2240700)],
        ),
        (AngularMomentum::S, vec![(0.0773500000, 1.0000000)]),
    ]
}

fn b321g_b() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (116.4340000000, 0.0629605),
                (17.4314000000, 0.3633040),
                (3.6801600000, 0.6972550),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(2.2818700000, -0.3687570), (0.4652480000, 1.1994600)],
        ),
        (AngularMomentum::S, vec![(0.1243280000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![(2.2818700000, 0.2311390), (0.4652480000, 0.8668700)],
        ),
        (AngularMomentum::P, vec![(0.1243280000, 1.0000000)]),
    ]
}

fn b321g_c() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (172.2560000000, 0.0617669),
                (25.9109000000, 0.3587940),
                (5.5333500000, 0.7007130),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(3.6649800000, -0.3958970), (0.7705450000, 1.2158400)],
        ),
        (AngularMomentum::S, vec![(0.1958570000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![(3.6649800000, 0.2364600), (0.7705450000, 0.8606190)],
        ),
        (AngularMomentum::P, vec![(0.1958570000, 1.0000000)]),
    ]
}

fn b321g_n() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (242.7660000000, 0.0598657),
                (36.4851000000, 0.3529550),
                (7.8144900000, 0.7065130),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(5.4252200000, -0.4133010), (1.1491500000, 1.2244200)],
        ),
        (AngularMomentum::S, vec![(0.2832050000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![(5.4252200000, 0.2379720), (1.1491500000, 0.8589530)],
        ),
        (AngularMomentum::P, vec![(0.2832050000, 1.0000000)]),
    ]
}

fn b321g_o() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (322.0370000000, 0.0592394),
                (48.4308000000, 0.3515000),
                (10.4206000000, 0.7076580),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(7.4029400000, -0.4044530), (1.5762000000, 1.2215600)],
        ),
        (AngularMomentum::S, vec![(0.3736840000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![(7.4029400000, 0.2445860), (1.5762000000, 0.8539550)],
        ),
        (AngularMomentum::P, vec![(0.3736840000, 1.0000000)]),
    ]
}

fn b321g_f() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (413.8010000000, 0.0585483),
                (62.2446000000, 0.3493080),
                (13.4340000000, 0.7096320),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(9.7775900000, -0.4073270), (2.0861700000, 1.2231400)],
        ),
        (AngularMomentum::S, vec![(0.4823830000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![(9.7775900000, 0.2466800), (2.0861700000, 0.8523210)],
        ),
        (AngularMomentum::P, vec![(0.4823830000, 1.0000000)]),
    ]
}

fn b321g_ne() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (515.7240000000, 0.0581430),
                (77.6538000000, 0.3479510),
                (16.8136000000, 0.7107140),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(12.4830000000, -0.4099300), (2.6645100000, 1.2243800)],
        ),
        (AngularMomentum::S, vec![(0.6062500000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![(12.4830000000, 0.2472050), (2.6645100000, 0.8518030)],
        ),
        (AngularMomentum::P, vec![(0.6062500000, 1.0000000)]),
    ]
}

// =============================================================================
// 6-31G Basis Set Data
// Reference: PySCF 2.11.0 gto/basis/6-31g.dat
// =============================================================================

fn get_631g(z: u8) -> Result<Vec<ShellData>, BasisError> {
    match z {
        1 => Ok(b631g_h()),
        2 => Ok(b631g_he()),
        3 => Ok(b631g_li()),
        4 => Ok(b631g_be()),
        5 => Ok(b631g_b()),
        6 => Ok(b631g_c()),
        7 => Ok(b631g_n()),
        8 => Ok(b631g_o()),
        9 => Ok(b631g_f()),
        10 => Ok(b631g_ne()),
        _ => Err(BasisError::ElementNotInBasis(z, "6-31g".to_string())),
    }
}

fn b631g_h() -> Vec<ShellData> {
    // Coefficients matched to PySCF 2.11.0 for consistency
    // Reference: pyscf/gto/basis/6-31g.dat
    vec![
        (
            AngularMomentum::S,
            vec![
                (18.7311370000, 0.0334946),
                (2.8253937000, 0.23472695),
                (0.6401217000, 0.81375733),
            ],
        ),
        (AngularMomentum::S, vec![(0.1612778000, 1.0000000000)]),
    ]
}

fn b631g_he() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (38.4216340000, 0.02376600),
                (5.7780300000, 0.15467900),
                (1.2417740000, 0.46963000),
            ],
        ),
        (AngularMomentum::S, vec![(0.2979640000, 1.0000000)]),
    ]
}

fn b631g_li() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (642.4189200000, 0.00214260),
                (96.7985150000, 0.01620890),
                (22.0911210000, 0.07731560),
                (6.2010703000, 0.24578600),
                (1.9351177000, 0.47018900),
                (0.6367358000, 0.34547080),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.3249184000, -0.03509170),
                (0.6324306000, -0.19123280),
                (0.0790534000, 1.08398620),
            ],
        ),
        (AngularMomentum::S, vec![(0.0359620000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![
                (2.3249184000, 0.00893560),
                (0.6324306000, 0.14100950),
                (0.0790534000, 0.94536370),
            ],
        ),
        (AngularMomentum::P, vec![(0.0359620000, 1.0000000)]),
    ]
}

fn b631g_be() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1264.5857000000, 0.00194475),
                (189.9368100000, 0.01483510),
                (43.1590890000, 0.07209060),
                (12.0986630000, 0.23715420),
                (3.8063232000, 0.46919870),
                (1.2728903000, 0.35652020),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.1964631000, -0.11264870),
                (0.7478133000, -0.22950640),
                (0.2199663000, 1.18691670),
            ],
        ),
        (AngularMomentum::S, vec![(0.0823099000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![
                (3.1964631000, 0.05598020),
                (0.7478133000, 0.26155060),
                (0.2199663000, 0.79397230),
            ],
        ),
        (AngularMomentum::P, vec![(0.0823099000, 1.0000000)]),
    ]
}

fn b631g_b() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2068.8823000000, 0.00186627),
                (310.6495700000, 0.01425150),
                (70.6830330000, 0.06955160),
                (19.8610800000, 0.23257290),
                (6.2993048000, 0.46707870),
                (2.1270270000, 0.36343140),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.7279710000, -0.10392900),
                (1.1903377000, -0.19835100),
                (0.3594117000, 1.15996100),
            ],
        ),
        (AngularMomentum::S, vec![(0.1267512000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![
                (4.7279710000, 0.07459300),
                (1.1903377000, 0.30784700),
                (0.3594117000, 0.74346500),
            ],
        ),
        (AngularMomentum::P, vec![(0.1267512000, 1.0000000)]),
    ]
}

fn b631g_c() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (3047.5249000000, 0.00183474),
                (457.3695100000, 0.01403730),
                (103.9486900000, 0.06884260),
                (29.2101550000, 0.23218440),
                (9.2866630000, 0.46794130),
                (3.1639270000, 0.36231200),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (7.8682724000, -0.11933240),
                (1.8812885000, -0.16085420),
                (0.5442493000, 1.14345640),
            ],
        ),
        (AngularMomentum::S, vec![(0.1687144000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![
                (7.8682724000, 0.06899910),
                (1.8812885000, 0.31642400),
                (0.5442493000, 0.74430830),
            ],
        ),
        (AngularMomentum::P, vec![(0.1687144000, 1.0000000)]),
    ]
}

fn b631g_n() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (4173.5110000000, 0.00183480),
                (627.4579000000, 0.01399500),
                (142.9021000000, 0.06858700),
                (40.2343300000, 0.23224100),
                (12.8202100000, 0.46907000),
                (4.3904370000, 0.36045500),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (11.6263580000, -0.11496100),
                (2.7162800000, -0.16911800),
                (0.7722180000, 1.14585200),
            ],
        ),
        (AngularMomentum::S, vec![(0.2120313000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![
                (11.6263580000, 0.06758000),
                (2.7162800000, 0.32390700),
                (0.7722180000, 0.74089500),
            ],
        ),
        (AngularMomentum::P, vec![(0.2120313000, 1.0000000)]),
    ]
}

fn b631g_o() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (5484.6717000000, 0.00183110),
                (825.2349500000, 0.01395010),
                (188.0469600000, 0.06844510),
                (52.9645000000, 0.23271430),
                (16.8975700000, 0.47019300),
                (5.7996353000, 0.35852090),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (15.5396160000, -0.11077750),
                (3.5999336000, -0.14802630),
                (1.0137618000, 1.13076700),
            ],
        ),
        (AngularMomentum::S, vec![(0.2700058000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![
                (15.5396160000, 0.07087430),
                (3.5999336000, 0.33975280),
                (1.0137618000, 0.72715860),
            ],
        ),
        (AngularMomentum::P, vec![(0.2700058000, 1.0000000)]),
    ]
}

fn b631g_f() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (7001.7130000000, 0.00181960),
                (1051.3660000000, 0.01391600),
                (239.2857000000, 0.06840500),
                (67.3974500000, 0.23318600),
                (21.5199600000, 0.47126700),
                (7.4031010000, 0.35661900),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (20.8479400000, -0.10850700),
                (4.8083080000, -0.14645200),
                (1.3440700000, 1.12868900),
            ],
        ),
        (AngularMomentum::S, vec![(0.3581514000, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![
                (20.8479400000, 0.07162900),
                (4.8083080000, 0.34591200),
                (1.3440700000, 0.72247000),
            ],
        ),
        (AngularMomentum::P, vec![(0.3581514000, 1.0000000)]),
    ]
}

fn b631g_ne() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (8425.8515300, 0.00188435),
                (1268.5194000, 0.01433690),
                (289.6214140, 0.07010960),
                (81.8590040, 0.23737280),
                (26.2515079, 0.47300710),
                (9.0947205, 0.34840400),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (26.5321310, -0.10711880),
                (6.1017550, -0.14616380),
                (1.6962715, 1.12777350),
            ],
        ),
        (AngularMomentum::S, vec![(0.4458187, 1.0000000)]),
        (
            AngularMomentum::P,
            vec![
                (26.5321310, 0.07190960),
                (6.1017550, 0.34951340),
                (1.6962715, 0.71994050),
            ],
        ),
        (AngularMomentum::P, vec![(0.4458187, 1.0000000)]),
    ]
}

// =============================================================================
// 6-31G* Basis Set Data (6-31G with d polarization on heavy atoms)
// Reference: PySCF 2.11.0 gto/basis/6-31gs.dat
// =============================================================================

fn get_631gs(z: u8) -> Result<Vec<ShellData>, BasisError> {
    // 6-31G* = 6-31G + d polarization on Li-Ne
    let mut shells = get_631g(z)?;

    // Add d polarization functions for Li-Ne (Z >= 3)
    if z >= 3 {
        let d_exp = match z {
            3 => 0.2000000,  // Li
            4 => 0.4000000,  // Be
            5 => 0.6000000,  // B
            6 => 0.8000000,  // C
            7 => 0.8000000,  // N
            8 => 0.8000000,  // O
            9 => 0.8000000,  // F
            10 => 0.8000000, // Ne
            _ => return Err(BasisError::ElementNotInBasis(z, "6-31g*".to_string())),
        };
        shells.push((AngularMomentum::D, vec![(d_exp, 1.0000000)]));
    }

    Ok(shells)
}

// =============================================================================
// 6-31+G* Basis Set Data (6-31G with diffuse sp + d polarization on heavy atoms)
// Reference: Basis Set Exchange (https://www.basissetexchange.org/)
// 6-31+G* = 6-31G + diffuse sp functions on Li-Ne + d polarization on Li-Ne
// Note: H and He do NOT get diffuse functions in this basis set
// =============================================================================

fn get_631pgs(z: u8) -> Result<Vec<ShellData>, BasisError> {
    match z {
        1 => Ok(b631pgs_h()),
        2 => Ok(b631pgs_he()),
        3 => Ok(b631pgs_li()),
        4 => Ok(b631pgs_be()),
        5 => Ok(b631pgs_b()),
        6 => Ok(b631pgs_c()),
        7 => Ok(b631pgs_n()),
        8 => Ok(b631pgs_o()),
        9 => Ok(b631pgs_f()),
        10 => Ok(b631pgs_ne()),
        _ => Err(BasisError::ElementNotInBasis(z, "6-31+g*".to_string())),
    }
}

/// H 6-31+G* basis (same as 6-31G - no diffuse functions on hydrogen)
fn b631pgs_h() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (18.7311369600, 0.0334946043),
                (2.8253943650, 0.2347269535),
                (0.6401216923, 0.8137573261),
            ],
        ),
        (AngularMomentum::S, vec![(0.1612777588, 1.0000000000)]),
    ]
}

/// He 6-31+G* basis (same as 6-31G - no diffuse functions on helium)
fn b631pgs_he() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (38.4216340000, 0.0401397394),
                (5.7780300000, 0.2612460970),
                (1.2417740000, 0.7931846246),
            ],
        ),
        (AngularMomentum::S, vec![(0.2979640000, 1.0000000000)]),
    ]
}

/// Li 6-31+G* basis
fn b631pgs_li() -> Vec<ShellData> {
    vec![
        // 1s shell (6 primitives)
        (
            AngularMomentum::S,
            vec![
                (642.4189150000, 0.0021426078),
                (96.7985153000, 0.0162088715),
                (22.0911212000, 0.0773155725),
                (6.2010702500, 0.2457860520),
                (1.9351176800, 0.4701890040),
                (0.6367357890, 0.3454708450),
            ],
        ),
        // SP shell - S component (inner valence)
        (
            AngularMomentum::S,
            vec![
                (2.3249184080, -0.0350917457),
                (0.6324303556, -0.1912328431),
                (0.0790534347, 1.0839877950),
            ],
        ),
        // SP shell - P component (inner valence)
        (
            AngularMomentum::P,
            vec![
                (2.3249184080, 0.0089415080),
                (0.6324303556, 0.1410094640),
                (0.0790534347, 0.9453636953),
            ],
        ),
        // SP shell - S component (outer valence)
        (AngularMomentum::S, vec![(0.0359619718, 1.0000000000)]),
        // SP shell - P component (outer valence)
        (AngularMomentum::P, vec![(0.0359619718, 1.0000000000)]),
        // d polarization
        (AngularMomentum::D, vec![(0.2000000000, 1.0000000000)]),
        // Diffuse SP shell - S component
        (AngularMomentum::S, vec![(0.0074000000, 1.0000000000)]),
        // Diffuse SP shell - P component
        (AngularMomentum::P, vec![(0.0074000000, 1.0000000000)]),
    ]
}

/// Be 6-31+G* basis
fn b631pgs_be() -> Vec<ShellData> {
    vec![
        // 1s shell (6 primitives)
        (
            AngularMomentum::S,
            vec![
                (1264.5856900000, 0.0019447576),
                (189.9368060000, 0.0148350520),
                (43.1590890000, 0.0720905463),
                (12.0986627000, 0.2371541500),
                (3.8063232200, 0.4691986519),
                (1.2728903000, 0.3565202279),
            ],
        ),
        // SP shell - S component (inner valence)
        (
            AngularMomentum::S,
            vec![
                (3.1964630980, -0.1126487285),
                (0.7478133038, -0.2295064079),
                (0.2199663302, 1.1869167640),
            ],
        ),
        // SP shell - P component (inner valence)
        (
            AngularMomentum::P,
            vec![
                (3.1964630980, 0.0559801998),
                (0.7478133038, 0.2615506110),
                (0.2199663302, 0.7939723389),
            ],
        ),
        // SP shell - S component (outer valence)
        (AngularMomentum::S, vec![(0.0823099007, 1.0000000000)]),
        // SP shell - P component (outer valence)
        (AngularMomentum::P, vec![(0.0823099007, 1.0000000000)]),
        // d polarization
        (AngularMomentum::D, vec![(0.4000000000, 1.0000000000)]),
        // Diffuse SP shell - S component
        (AngularMomentum::S, vec![(0.0207000000, 1.0000000000)]),
        // Diffuse SP shell - P component
        (AngularMomentum::P, vec![(0.0207000000, 1.0000000000)]),
    ]
}

/// B 6-31+G* basis
fn b631pgs_b() -> Vec<ShellData> {
    vec![
        // 1s shell (6 primitives)
        (
            AngularMomentum::S,
            vec![
                (2068.8822500000, 0.0018662746),
                (310.6495700000, 0.0142514817),
                (70.6830330000, 0.0695516185),
                (19.8610803000, 0.2325729330),
                (6.2993048400, 0.4670787120),
                (2.1270269700, 0.3634314400),
            ],
        ),
        // SP shell - S component (inner valence)
        (
            AngularMomentum::S,
            vec![
                (4.7279710710, -0.1303937974),
                (1.1903377360, -0.1307889514),
                (0.3594116829, 1.1309444840),
            ],
        ),
        // SP shell - P component (inner valence)
        (
            AngularMomentum::P,
            vec![
                (4.7279710710, 0.0745975799),
                (1.1903377360, 0.3078466771),
                (0.3594116829, 0.7434568342),
            ],
        ),
        // SP shell - S component (outer valence)
        (AngularMomentum::S, vec![(0.1267512469, 1.0000000000)]),
        // SP shell - P component (outer valence)
        (AngularMomentum::P, vec![(0.1267512469, 1.0000000000)]),
        // d polarization
        (AngularMomentum::D, vec![(0.6000000000, 1.0000000000)]),
        // Diffuse SP shell - S component
        (AngularMomentum::S, vec![(0.0315000000, 1.0000000000)]),
        // Diffuse SP shell - P component
        (AngularMomentum::P, vec![(0.0315000000, 1.0000000000)]),
    ]
}

/// C 6-31+G* basis
fn b631pgs_c() -> Vec<ShellData> {
    vec![
        // 1s shell (6 primitives)
        (
            AngularMomentum::S,
            vec![
                (3047.5248800000, 0.0018347371),
                (457.3695180000, 0.0140373228),
                (103.9486850000, 0.0688426223),
                (29.2101553000, 0.2321844432),
                (9.2866629600, 0.4679413484),
                (3.1639269600, 0.3623119853),
            ],
        ),
        // SP shell - S component (inner valence)
        (
            AngularMomentum::S,
            vec![
                (7.8682723500, -0.1193324198),
                (1.8812885400, -0.1608541517),
                (0.5442492580, 1.1434564380),
            ],
        ),
        // SP shell - P component (inner valence)
        (
            AngularMomentum::P,
            vec![
                (7.8682723500, 0.0689990666),
                (1.8812885400, 0.3164239610),
                (0.5442492580, 0.7443082909),
            ],
        ),
        // SP shell - S component (outer valence)
        (AngularMomentum::S, vec![(0.1687144782, 1.0000000000)]),
        // SP shell - P component (outer valence)
        (AngularMomentum::P, vec![(0.1687144782, 1.0000000000)]),
        // d polarization
        (AngularMomentum::D, vec![(0.8000000000, 1.0000000000)]),
        // Diffuse SP shell - S component
        (AngularMomentum::S, vec![(0.0438000000, 1.0000000000)]),
        // Diffuse SP shell - P component
        (AngularMomentum::P, vec![(0.0438000000, 1.0000000000)]),
    ]
}

/// N 6-31+G* basis
fn b631pgs_n() -> Vec<ShellData> {
    vec![
        // 1s shell (6 primitives)
        (
            AngularMomentum::S,
            vec![
                (4173.5114600000, 0.0018347722),
                (627.4579110000, 0.0139946270),
                (142.9020930000, 0.0685865518),
                (40.2343293000, 0.2322408730),
                (12.8202129000, 0.4690699481),
                (4.3904370100, 0.3604551991),
            ],
        ),
        // SP shell - S component (inner valence)
        (
            AngularMomentum::S,
            vec![
                (11.6263618600, -0.1149611817),
                (2.7162798070, -0.1691174786),
                (0.7722183966, 1.1458519470),
            ],
        ),
        // SP shell - P component (inner valence)
        (
            AngularMomentum::P,
            vec![
                (11.6263618600, 0.0675797439),
                (2.7162798070, 0.3239072959),
                (0.7722183966, 0.7408951398),
            ],
        ),
        // SP shell - S component (outer valence)
        (AngularMomentum::S, vec![(0.2120314975, 1.0000000000)]),
        // SP shell - P component (outer valence)
        (AngularMomentum::P, vec![(0.2120314975, 1.0000000000)]),
        // d polarization
        (AngularMomentum::D, vec![(0.8000000000, 1.0000000000)]),
        // Diffuse SP shell - S component
        (AngularMomentum::S, vec![(0.0639000000, 1.0000000000)]),
        // Diffuse SP shell - P component
        (AngularMomentum::P, vec![(0.0639000000, 1.0000000000)]),
    ]
}

/// O 6-31+G* basis
fn b631pgs_o() -> Vec<ShellData> {
    vec![
        // 1s shell (6 primitives)
        (
            AngularMomentum::S,
            vec![
                (5484.6716600000, 0.0018310744),
                (825.2349460000, 0.0139501722),
                (188.0469580000, 0.0684450781),
                (52.9645000000, 0.2327143360),
                (16.8975704000, 0.4701928980),
                (5.7996353400, 0.3585208530),
            ],
        ),
        // SP shell - S component (inner valence)
        (
            AngularMomentum::S,
            vec![
                (15.5396162500, -0.1107775495),
                (3.5999335860, -0.1480262627),
                (1.0137617500, 1.1307670150),
            ],
        ),
        // SP shell - P component (inner valence)
        (
            AngularMomentum::P,
            vec![
                (15.5396162500, 0.0708742682),
                (3.5999335860, 0.3397528391),
                (1.0137617500, 0.7271585773),
            ],
        ),
        // SP shell - S component (outer valence)
        (AngularMomentum::S, vec![(0.2700058226, 1.0000000000)]),
        // SP shell - P component (outer valence)
        (AngularMomentum::P, vec![(0.2700058226, 1.0000000000)]),
        // d polarization
        (AngularMomentum::D, vec![(0.8000000000, 1.0000000000)]),
        // Diffuse SP shell - S component
        (AngularMomentum::S, vec![(0.0845000000, 1.0000000000)]),
        // Diffuse SP shell - P component
        (AngularMomentum::P, vec![(0.0845000000, 1.0000000000)]),
    ]
}

/// F 6-31+G* basis
fn b631pgs_f() -> Vec<ShellData> {
    vec![
        // 1s shell (6 primitives)
        (
            AngularMomentum::S,
            vec![
                (7001.7130900000, 0.0018196169),
                (1051.3660900000, 0.0139160796),
                (239.2856900000, 0.0684053245),
                (67.3974453000, 0.2331857601),
                (21.5199573000, 0.4712674392),
                (7.4031013000, 0.3566185462),
            ],
        ),
        // SP shell - S component (inner valence)
        (
            AngularMomentum::S,
            vec![
                (20.8479528000, -0.1085069751),
                (4.8083083400, -0.1464516581),
                (1.3440698600, 1.1286885810),
            ],
        ),
        // SP shell - P component (inner valence)
        (
            AngularMomentum::P,
            vec![
                (20.8479528000, 0.0716287242),
                (4.8083083400, 0.3459121027),
                (1.3440698600, 0.7224699564),
            ],
        ),
        // SP shell - S component (outer valence)
        (AngularMomentum::S, vec![(0.3581513930, 1.0000000000)]),
        // SP shell - P component (outer valence)
        (AngularMomentum::P, vec![(0.3581513930, 1.0000000000)]),
        // d polarization
        (AngularMomentum::D, vec![(0.8000000000, 1.0000000000)]),
        // Diffuse SP shell - S component
        (AngularMomentum::S, vec![(0.1076000000, 1.0000000000)]),
        // Diffuse SP shell - P component
        (AngularMomentum::P, vec![(0.1076000000, 1.0000000000)]),
    ]
}

/// Ne 6-31+G* basis
fn b631pgs_ne() -> Vec<ShellData> {
    vec![
        // 1s shell (6 primitives)
        (
            AngularMomentum::S,
            vec![
                (8425.8515300000, 0.0018843480),
                (1268.5194000000, 0.0143368994),
                (289.6214140000, 0.0701096233),
                (81.8590040000, 0.2373732660),
                (26.2515079000, 0.4730071261),
                (9.0947205100, 0.3484012410),
            ],
        ),
        // SP shell - S component (inner valence)
        (
            AngularMomentum::S,
            vec![
                (26.5321310000, -0.1071182872),
                (6.1017550100, -0.1461638213),
                (1.6962715300, 1.1277735030),
            ],
        ),
        // SP shell - P component (inner valence)
        (
            AngularMomentum::P,
            vec![
                (26.5321310000, 0.0719095885),
                (6.1017550100, 0.3495133720),
                (1.6962715300, 0.7199405121),
            ],
        ),
        // SP shell - S component (outer valence)
        (AngularMomentum::S, vec![(0.4458187000, 1.0000000000)]),
        // SP shell - P component (outer valence)
        (AngularMomentum::P, vec![(0.4458187000, 1.0000000000)]),
        // d polarization
        (AngularMomentum::D, vec![(0.8000000000, 1.0000000000)]),
        // Diffuse SP shell - S component
        (AngularMomentum::S, vec![(0.1300000000, 1.0000000000)]),
        // Diffuse SP shell - P component
        (AngularMomentum::P, vec![(0.1300000000, 1.0000000000)]),
    ]
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------------------------------------------------------
    // API tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_get_element_basis_h_sto3g() {
        let basis = get_element_basis(1, "sto-3g").unwrap();

        // H has 1 shell (1s)
        assert_eq!(basis.len(), 1);

        // Check angular momentum
        assert_eq!(basis[0].0, AngularMomentum::S);

        // Check number of primitives (STO-3G = 3 primitives)
        assert_eq!(basis[0].1.len(), 3);
    }

    #[test]
    fn test_get_element_basis_c_sto3g() {
        let basis = get_element_basis(6, "sto-3g").unwrap();

        // C has 3 shells (1s, 2s, 2p)
        assert_eq!(basis.len(), 3);

        // Check angular momenta
        assert_eq!(basis[0].0, AngularMomentum::S); // 1s
        assert_eq!(basis[1].0, AngularMomentum::S); // 2s
        assert_eq!(basis[2].0, AngularMomentum::P); // 2p

        // Each shell has 3 primitives in STO-3G
        for shell in &basis {
            assert_eq!(shell.1.len(), 3);
        }
    }

    #[test]
    fn test_get_element_basis_o_sto3g() {
        let basis = get_element_basis(8, "sto-3g").unwrap();

        // O has 3 shells (1s, 2s, 2p)
        assert_eq!(basis.len(), 3);
        assert_eq!(basis[0].0, AngularMomentum::S);
        assert_eq!(basis[1].0, AngularMomentum::S);
        assert_eq!(basis[2].0, AngularMomentum::P);
    }

    #[test]
    fn test_get_element_basis_case_insensitive() {
        // All these should work
        assert!(get_element_basis(1, "STO-3G").is_ok());
        assert!(get_element_basis(1, "sto-3g").is_ok());
        assert!(get_element_basis(1, "Sto-3g").is_ok());
        assert!(get_element_basis(1, "sto3g").is_ok());
    }

    #[test]
    fn test_get_element_basis_unknown() {
        let result = get_element_basis(1, "unknown-basis");
        assert!(matches!(result, Err(BasisError::UnknownBasis(_))));
    }

    #[test]
    fn test_get_element_basis_unsupported_element() {
        // Z = 11 (Sodium) not supported
        let result = get_element_basis(11, "sto-3g");
        assert!(matches!(result, Err(BasisError::UnsupportedElement(11))));

        // Z = 0 invalid
        let result = get_element_basis(0, "sto-3g");
        assert!(matches!(result, Err(BasisError::UnsupportedElement(0))));
    }

    #[test]
    fn test_supported_basis_sets() {
        let basis_sets = supported_basis_sets();
        assert!(basis_sets.contains(&"sto-3g"));
        assert!(basis_sets.contains(&"3-21g"));
        assert!(basis_sets.contains(&"6-31g"));
        assert!(basis_sets.contains(&"6-31g*"));
        assert!(basis_sets.contains(&"6-31+g*"));
    }

    #[test]
    fn test_is_supported_basis() {
        assert!(is_supported_basis("sto-3g"));
        assert!(is_supported_basis("STO-3G"));
        assert!(is_supported_basis("3-21g"));
        assert!(is_supported_basis("6-31g"));
        assert!(is_supported_basis("6-31g*"));
        assert!(is_supported_basis("6-31g(d)"));
        assert!(is_supported_basis("6-31+g*"));
        assert!(is_supported_basis("6-31+G*"));
        assert!(is_supported_basis("6-31+g(d)"));
        assert!(!is_supported_basis("cc-pVDZ"));
        assert!(!is_supported_basis("unknown"));
    }

    // -------------------------------------------------------------------------
    // 3-21G tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_321g_h() {
        let basis = get_element_basis(1, "3-21g").unwrap();

        // H 3-21G has 2 S shells (split valence)
        assert_eq!(basis.len(), 2);
        assert_eq!(basis[0].0, AngularMomentum::S);
        assert_eq!(basis[1].0, AngularMomentum::S);

        // Inner: 2 primitives, Outer: 1 primitive
        assert_eq!(basis[0].1.len(), 2);
        assert_eq!(basis[1].1.len(), 1);
    }

    #[test]
    fn test_321g_c() {
        let basis = get_element_basis(6, "3-21g").unwrap();

        // C 3-21G: 1s(3), 2s(2), 2s(1), 2p(2), 2p(1) = 5 shells
        assert_eq!(basis.len(), 5);
    }

    // -------------------------------------------------------------------------
    // 6-31G tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_631g_h() {
        let basis = get_element_basis(1, "6-31g").unwrap();

        // H 6-31G has 2 S shells
        assert_eq!(basis.len(), 2);

        // Inner: 3 primitives, Outer: 1 primitive
        assert_eq!(basis[0].1.len(), 3);
        assert_eq!(basis[1].1.len(), 1);
    }

    #[test]
    fn test_631g_c() {
        let basis = get_element_basis(6, "6-31g").unwrap();

        // C 6-31G: 1s(6), 2s(3), 2s(1), 2p(3), 2p(1) = 5 shells
        assert_eq!(basis.len(), 5);

        // Check shell types
        assert_eq!(basis[0].0, AngularMomentum::S); // 1s
        assert_eq!(basis[1].0, AngularMomentum::S); // 2s inner
        assert_eq!(basis[2].0, AngularMomentum::S); // 2s outer
        assert_eq!(basis[3].0, AngularMomentum::P); // 2p inner
        assert_eq!(basis[4].0, AngularMomentum::P); // 2p outer
    }

    // -------------------------------------------------------------------------
    // 6-31G* tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_631gs_h() {
        let basis = get_element_basis(1, "6-31g*").unwrap();

        // H has no d functions (same as 6-31G)
        assert_eq!(basis.len(), 2);
        assert_eq!(basis[0].0, AngularMomentum::S);
        assert_eq!(basis[1].0, AngularMomentum::S);
    }

    #[test]
    fn test_631gs_c() {
        let basis = get_element_basis(6, "6-31g*").unwrap();

        // C 6-31G* = 6-31G + 1 d shell = 6 shells
        assert_eq!(basis.len(), 6);

        // Last shell should be D
        assert_eq!(basis[5].0, AngularMomentum::D);
        assert_eq!(basis[5].1.len(), 1); // Single d exponent
    }

    #[test]
    fn test_631gs_o() {
        let basis = get_element_basis(8, "6-31g*").unwrap();

        // O 6-31G* = 6-31G + 1 d shell
        let last_shell = basis.last().unwrap();
        assert_eq!(last_shell.0, AngularMomentum::D);

        // Check d exponent is 0.8 for oxygen
        assert!((last_shell.1[0].0 - 0.8).abs() < 1e-6);
    }

    // -------------------------------------------------------------------------
    // 6-31+G* tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_631pgs_h() {
        let basis = get_element_basis(1, "6-31+g*").unwrap();

        // H has no diffuse functions in 6-31+G* (same as 6-31G)
        assert_eq!(basis.len(), 2);
        assert_eq!(basis[0].0, AngularMomentum::S);
        assert_eq!(basis[1].0, AngularMomentum::S);
    }

    #[test]
    fn test_631pgs_c() {
        let basis = get_element_basis(6, "6-31+g*").unwrap();

        // C 6-31+G* = 6-31G + d polarization + diffuse sp = 8 shells
        // Shells: 1s(S), 2s inner(S), 2p inner(P), 2s outer(S), 2p outer(P), d(D), diffuse S, diffuse P
        assert_eq!(basis.len(), 8);

        // Check that we have the expected shell types
        let shell_types: Vec<_> = basis.iter().map(|(am, _)| *am).collect();
        assert!(shell_types.contains(&AngularMomentum::S));
        assert!(shell_types.contains(&AngularMomentum::P));
        assert!(shell_types.contains(&AngularMomentum::D));

        // Check that diffuse exponent is present (0.0438 for C)
        let diffuse_s = basis
            .iter()
            .find(|(am, prims)| *am == AngularMomentum::S && prims.len() == 1 && prims[0].0 < 0.05);
        assert!(diffuse_s.is_some(), "Should have diffuse S shell");
    }

    #[test]
    fn test_631pgs_o() {
        let basis = get_element_basis(8, "6-31+g*").unwrap();

        // O 6-31+G* = 6-31G + d polarization + diffuse sp = 8 shells
        assert_eq!(basis.len(), 8);

        // Check d polarization exponent is 0.8
        let d_shell = basis.iter().find(|(am, _)| *am == AngularMomentum::D);
        assert!(d_shell.is_some());
        assert!((d_shell.unwrap().1[0].0 - 0.8).abs() < 1e-6);

        // Check diffuse exponent (0.0845 for O)
        let diffuse_s = basis.iter().find(|(am, prims)| {
            *am == AngularMomentum::S && prims.len() == 1 && (prims[0].0 - 0.0845).abs() < 1e-6
        });
        assert!(
            diffuse_s.is_some(),
            "Should have diffuse S shell with exp ~0.0845"
        );
    }

    #[test]
    fn test_631pgs_basis_count_h2o() {
        // H2O 6-31+G*: 1 O + 2 H
        // O 6-31+G*: 4 S + 3 P + 1 D = 4 + 9 + 6 = 19 basis functions (Cartesian)
        // H 6-31+G*: 2 S = 2 basis functions (same as 6-31G)
        // Total: 19 + 2 + 2 = 23 (Cartesian)

        let o_basis = get_element_basis(8, "6-31+g*").unwrap();
        let h_basis = get_element_basis(1, "6-31+g*").unwrap();

        let o_nao: usize = o_basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        let h_nao: usize = h_basis.iter().map(|(am, _)| am.n_cartesian()).sum();

        // O: 4 S shells (1 each) + 3 P shells (3 each) + 1 D shell (6) = 4 + 9 + 6 = 19
        assert_eq!(o_nao, 19);
        // H: 2 S shells = 2
        assert_eq!(h_nao, 2);
        // Total for H2O
        assert_eq!(o_nao + 2 * h_nao, 23);
    }

    // -------------------------------------------------------------------------
    // Basis function count tests (from golden reference)
    // -------------------------------------------------------------------------

    #[test]
    fn test_h2_sto3g_basis_count() {
        // H2 STO-3G: 2 H atoms, each with 1 S shell
        // nao = 2, nbas = 2
        let h_basis = get_element_basis(1, "sto-3g").unwrap();

        let h_nao: usize = h_basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(h_nao, 1); // H has 1 basis function

        // Total for H2: 2 * 1 = 2
        assert_eq!(h_nao * 2, 2);
    }

    #[test]
    fn test_h2o_sto3g_basis_count() {
        // H2O STO-3G: 1 O + 2 H
        // O: 1s + 2s + 2p = 1 + 1 + 3 = 5 basis functions
        // H: 1s = 1 basis function
        // Total: 5 + 1 + 1 = 7

        let o_basis = get_element_basis(8, "sto-3g").unwrap();
        let h_basis = get_element_basis(1, "sto-3g").unwrap();

        let o_nao: usize = o_basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        let h_nao: usize = h_basis.iter().map(|(am, _)| am.n_cartesian()).sum();

        assert_eq!(o_nao, 5);
        assert_eq!(h_nao, 1);
        assert_eq!(o_nao + 2 * h_nao, 7);
    }

    #[test]
    fn test_h2o_631gs_basis_count() {
        // H2O 6-31G*: 1 O + 2 H
        // O 6-31G*: 2 S + 1 S + 2 P + 1 P + 1 D = 1 + 1 + 3 + 3 + 6 = 14? No wait...
        // Actually: (1s) 1 + (2s inner) 1 + (2s outer) 1 + (2p inner) 3 + (2p outer) 3 + (d) 6 = 15
        // But that's not right either. Let me recalculate:
        // O 6-31G = 5 shells: 1s(6prims, 1 func), 2s(3prims, 1 func), 2s(1prim, 1 func),
        //                     2p(3prims, 3 funcs), 2p(1prim, 3 funcs) = 1+1+1+3+3 = 9
        // O 6-31G* = + d(1prim, 6 funcs) = 9 + 6 = 15
        // H 6-31G = 2 shells: 1s(3prims, 1 func), 1s(1prim, 1 func) = 2
        // Total: 15 + 2 + 2 = 19 -- but reference says 18!
        //
        // Let me check: the reference says nao=18 for H2O 6-31G*
        // O: 1s(1) + 2s(1) + 2s(1) + 2p(3) + 2p(3) + d(6) = 15 -- only if Cartesian d
        // But wait, with spherical d, d has 5 functions: 15 - 1 = 14
        // Actually, looking at the reference: nao = 18
        // Let me reconsider: maybe the reference uses spherical harmonics
        //
        // Actually the problem statement says: H2O 6-31G* nao=18
        // With Cartesian d (6 functions): 9 + 6 + 2 + 2 = 19
        // With spherical d (5 functions): 9 + 5 + 2 + 2 = 18 -- this matches!
        //
        // But our implementation uses Cartesian, so we should get 19.
        // The reference must be using spherical harmonics.

        let o_basis = get_element_basis(8, "6-31g*").unwrap();
        let h_basis = get_element_basis(1, "6-31g*").unwrap();

        let o_nao: usize = o_basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        let h_nao: usize = h_basis.iter().map(|(am, _)| am.n_cartesian()).sum();

        // With Cartesian d functions
        assert_eq!(o_nao, 15); // 1 + 1 + 1 + 3 + 3 + 6
        assert_eq!(h_nao, 2); // 1 + 1
                              // Total with Cartesian: 15 + 2 + 2 = 19
                              // (Reference with spherical d: 18)
    }
}
