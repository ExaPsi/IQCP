// Allow excessive precision for basis set coefficients from reference implementations.
// These values are taken directly from PySCF 2.11.0 for numerical accuracy.
#![allow(clippy::excessive_precision)]

//! Built-in basis set data
//!
//! This module provides pre-defined basis set parameters for elements H through Ar.
//! The data is sourced from PySCF 2.11.0 and uses pre-normalized coefficients.
//!
//! # Supported Basis Sets
//!
//! | Basis | Description | Elements |
//! |-------|-------------|----------|
//! | STO-3G | Minimal basis (3 Gaussians fit to Slater) | H-Ar |
//! | 3-21G | Split-valence | H-Ar |
//! | 6-31G | Split-valence | H-Ar |
//! | 6-31G* | Split-valence + d polarization | H-Ar |
//! | 6-31+G* | Split-valence + diffuse sp + d polarization | H-Ar |
//! | cc-pVDZ | Correlation-consistent polarized valence double-zeta | H-Ar |
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
    #[error(
        "Unknown basis set: '{0}' (supported: sto-3g, 3-21g, 6-31g, 6-31g*, 6-31+g*, cc-pvdz)"
    )]
    UnknownBasis(String),

    /// Element not available in basis set
    #[error("Element Z={0} not available in basis '{1}'")]
    ElementNotInBasis(u8, String),

    /// Element not supported by IQCP
    #[error("Unsupported element Z={0} (IQCP supports H-Ar, Z=1-18)")]
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
/// * `z` - Atomic number (1-18)
/// * `basis_name` - Basis set name (case insensitive): "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"
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
    if !(1..=18).contains(&z) {
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
        "cc-pvdz" | "ccpvdz" => get_ccpvdz(z),
        _ => Err(BasisError::UnknownBasis(basis_name.to_string())),
    }
}

/// Get a list of supported basis set names
pub fn supported_basis_sets() -> &'static [&'static str] {
    &["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"]
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
            | "cc-pvdz"
            | "ccpvdz"
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
        11 => Ok(sto3g_na()),
        12 => Ok(sto3g_mg()),
        13 => Ok(sto3g_al()),
        14 => Ok(sto3g_si()),
        15 => Ok(sto3g_p()),
        16 => Ok(sto3g_s()),
        17 => Ok(sto3g_cl()),
        18 => Ok(sto3g_ar()),
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
        11 => Ok(b321g_na()),
        12 => Ok(b321g_mg()),
        13 => Ok(b321g_al()),
        14 => Ok(b321g_si()),
        15 => Ok(b321g_p()),
        16 => Ok(b321g_s()),
        17 => Ok(b321g_cl()),
        18 => Ok(b321g_ar()),
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
        11 => Ok(b631g_na()),
        12 => Ok(b631g_mg()),
        13 => Ok(b631g_al()),
        14 => Ok(b631g_si()),
        15 => Ok(b631g_p()),
        16 => Ok(b631g_s()),
        17 => Ok(b631g_cl()),
        18 => Ok(b631g_ar()),
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
    // 6-31G* = 6-31G + d polarization on heavy atoms (Z >= 3)
    // For Z=11-18, use dedicated functions from PySCF extraction (include D shell)
    match z {
        1 | 2 => get_631g(z), // H, He: same as 6-31G (no d polarization)
        3..=10 => {
            let mut shells = get_631g(z)?;
            let d_exp = match z {
                3 => 0.2000000,  // Li
                4 => 0.4000000,  // Be
                5 => 0.6000000,  // B
                6 => 0.8000000,  // C
                7 => 0.8000000,  // N
                8 => 0.8000000,  // O
                9 => 0.8000000,  // F
                10 => 0.8000000, // Ne
                _ => unreachable!(),
            };
            shells.push((AngularMomentum::D, vec![(d_exp, 1.0000000)]));
            Ok(shells)
        }
        11 => Ok(b631gs_na()),
        12 => Ok(b631gs_mg()),
        13 => Ok(b631gs_al()),
        14 => Ok(b631gs_si()),
        15 => Ok(b631gs_p()),
        16 => Ok(b631gs_s()),
        17 => Ok(b631gs_cl()),
        18 => Ok(b631gs_ar()),
        _ => Err(BasisError::ElementNotInBasis(z, "6-31g*".to_string())),
    }
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
        11 => Ok(b631pgs_na()),
        12 => Ok(b631pgs_mg()),
        13 => Ok(b631pgs_al()),
        14 => Ok(b631pgs_si()),
        15 => Ok(b631pgs_p()),
        16 => Ok(b631pgs_s()),
        17 => Ok(b631pgs_cl()),
        18 => Ok(b631pgs_ar()),
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
// Third-Row Element Basis Data (Na-Ar, Z=11-18)
// All data extracted from PySCF 2.11.0 gto.basis.load()
// SP-combined shells split into separate S and P shells
// =============================================================================

// ====== STO-3G Na-Ar ======

/// Na (Z=11) STO-3G basis
fn sto3g_na() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.5077243000e+02, 1.5432896730e-01),
                (4.5678511000e+01, 5.3532814230e-01),
                (1.2362388000e+01, 4.4463454220e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.2040193000e+01, -9.9967229190e-02),
                (2.7978819000e+00, 3.9951282610e-01),
                (9.0995800000e-01, 7.0011546890e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.4787406000e+00, -2.1962036900e-01),
                (4.1256490000e-01, 2.2559543360e-01),
                (1.6147510000e-01, 9.0039842600e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.2040193000e+01, 1.5591627500e-01),
                (2.7978819000e+00, 6.0768371860e-01),
                (9.0995800000e-01, 3.9195739310e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.4787406000e+00, 1.0587604290e-02),
                (4.1256490000e-01, 5.9516700530e-01),
                (1.6147510000e-01, 4.6200101200e-01),
            ],
        ),
    ]
}

/// Mg (Z=12) STO-3G basis
fn sto3g_mg() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.9923740000e+02, 1.5432896730e-01),
                (5.4506470000e+01, 5.3532814230e-01),
                (1.4751580000e+01, 4.4463454220e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.5121820000e+01, -9.9967229190e-02),
                (3.5139870000e+00, 3.9951282610e-01),
                (1.1428570000e+00, 7.0011546890e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.3954480000e+00, -2.1962036900e-01),
                (3.8932600000e-01, 2.2559543360e-01),
                (1.5238000000e-01, 9.0039842600e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.5121820000e+01, 1.5591627500e-01),
                (3.5139870000e+00, 6.0768371860e-01),
                (1.1428570000e+00, 3.9195739310e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.3954480000e+00, 1.0587604290e-02),
                (3.8932600000e-01, 5.9516700530e-01),
                (1.5238000000e-01, 4.6200101200e-01),
            ],
        ),
    ]
}

/// Al (Z=13) STO-3G basis
fn sto3g_al() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (3.5142147670e+02, 1.5432896730e-01),
                (6.4011860670e+01, 5.3532814230e-01),
                (1.7324107610e+01, 4.4463454220e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.8899396210e+01, -9.9967229190e-02),
                (4.3918132330e+00, 3.9951282610e-01),
                (1.4283539700e+00, 7.0011546890e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.3954482930e+00, -2.1962036900e-01),
                (3.8932653180e-01, 2.2559543360e-01),
                (1.5237976590e-01, 9.0039842600e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.8899396210e+01, 1.5591627500e-01),
                (4.3918132330e+00, 6.0768371860e-01),
                (1.4283539700e+00, 3.9195739310e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.3954482930e+00, 1.0587604290e-02),
                (3.8932653180e-01, 5.9516700530e-01),
                (1.5237976590e-01, 4.6200101200e-01),
            ],
        ),
    ]
}

/// Si (Z=14) STO-3G basis
fn sto3g_si() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (4.0779755140e+02, 1.5432896730e-01),
                (7.4280833050e+01, 5.3532814230e-01),
                (2.0103292290e+01, 4.4463454220e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.3193656060e+01, -9.9967229190e-02),
                (5.3897068710e+00, 3.9951282610e-01),
                (1.7528999520e+00, 7.0011546890e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.4787406220e+00, -2.1962036900e-01),
                (4.1256488010e-01, 2.2559543360e-01),
                (1.6147509790e-01, 9.0039842600e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.3193656060e+01, 1.5591627500e-01),
                (5.3897068710e+00, 6.0768371860e-01),
                (1.7528999520e+00, 3.9195739310e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.4787406220e+00, 1.0587604290e-02),
                (4.1256488010e-01, 5.9516700530e-01),
                (1.6147509790e-01, 4.6200101200e-01),
            ],
        ),
    ]
}

/// P (Z=15) STO-3G basis
fn sto3g_p() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (4.6836563780e+02, 1.5432896730e-01),
                (8.5313385590e+01, 5.3532814230e-01),
                (2.3089131560e+01, 4.4463454220e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.8032639580e+01, -9.9967229190e-02),
                (6.5141825770e+00, 3.9951282610e-01),
                (2.1186143520e+00, 7.0011546890e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.7431032310e+00, -2.1962036900e-01),
                (4.8632137710e-01, 2.2559543360e-01),
                (1.9034289090e-01, 9.0039842600e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.8032639580e+01, 1.5591627500e-01),
                (6.5141825770e+00, 6.0768371860e-01),
                (2.1186143520e+00, 3.9195739310e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.7431032310e+00, 1.0587604290e-02),
                (4.8632137710e-01, 5.9516700530e-01),
                (1.9034289090e-01, 4.6200101200e-01),
            ],
        ),
    ]
}

/// S (Z=16) STO-3G basis
fn sto3g_s() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (5.3312573590e+02, 1.5432896730e-01),
                (9.7109518300e+01, 5.3532814230e-01),
                (2.6281625420e+01, 4.4463454220e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.3329751730e+01, -9.9967229190e-02),
                (7.7451175210e+00, 3.9951282610e-01),
                (2.5189525990e+00, 7.0011546890e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.0291942740e+00, -2.1962036900e-01),
                (5.6614005180e-01, 2.2559543360e-01),
                (2.2158337920e-01, 9.0039842600e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.3329751730e+01, 1.5591627500e-01),
                (7.7451175210e+00, 6.0768371860e-01),
                (2.5189525990e+00, 3.9195739310e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.0291942740e+00, 1.0587604290e-02),
                (5.6614005180e-01, 5.9516700530e-01),
                (2.2158337920e-01, 4.6200101200e-01),
            ],
        ),
    ]
}

/// Cl (Z=17) STO-3G basis
fn sto3g_cl() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (6.0134561360e+02, 1.5432896730e-01),
                (1.0953585420e+02, 5.3532814230e-01),
                (2.9644676860e+01, 4.4463454220e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.8960418890e+01, -9.9967229190e-02),
                (9.0535634770e+00, 3.9951282610e-01),
                (2.9444998340e+00, 7.0011546890e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.1293864950e+00, -2.1962036900e-01),
                (5.9409342740e-01, 2.2559543360e-01),
                (2.3252414100e-01, 9.0039842600e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.8960418890e+01, 1.5591627500e-01),
                (9.0535634770e+00, 6.0768371860e-01),
                (2.9444998340e+00, 3.9195739310e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.1293864950e+00, 1.0587604290e-02),
                (5.9409342740e-01, 5.9516700530e-01),
                (2.3252414100e-01, 4.6200101200e-01),
            ],
        ),
    ]
}

/// Ar (Z=18) STO-3G basis
fn sto3g_ar() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (6.7444651840e+02, 1.5432896730e-01),
                (1.2285127530e+02, 5.3532814230e-01),
                (3.3248349450e+01, 4.4463454220e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.5164243920e+01, -9.9967229190e-02),
                (1.0495199000e+01, 3.9951282610e-01),
                (3.4133644480e+00, 7.0011546890e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.6213665180e+00, -2.1962036900e-01),
                (7.3135460500e-01, 2.2559543360e-01),
                (2.8624723560e-01, 9.0039842600e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.5164243920e+01, 1.5591627500e-01),
                (1.0495199000e+01, 6.0768371860e-01),
                (3.4133644480e+00, 3.9195739310e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.6213665180e+00, 1.0587604290e-02),
                (7.3135460500e-01, 5.9516700530e-01),
                (2.8624723560e-01, 4.6200101200e-01),
            ],
        ),
    ]
}

// ====== 3-21G Na-Ar ======

/// Na (Z=11) 3-21G basis
fn b321g_na() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (5.4761300000e+02, 6.7491100000e-02),
                (8.2067800000e+01, 3.9350500000e-01),
                (1.7691700000e+01, 6.6560500000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.7540700000e+01, -1.1193700000e-01),
                (3.7939800000e+00, 2.5465400000e-01),
                (9.0644100000e-01, 8.4441700000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (5.0182400000e-01, -2.1966000000e-01),
                (6.0945800000e-02, 1.0891200000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(2.4434900000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.7540700000e+01, 1.2823300000e-01),
                (3.7939800000e+00, 4.7153300000e-01),
                (9.0644100000e-01, 6.0427300000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (5.0182400000e-01, 9.0665000000e-03),
                (6.0945800000e-02, 9.9720200000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(2.4434900000e-02, 1.0000000000e+00)],
        ),
    ]
}

/// Mg (Z=12) 3-21G basis
fn b321g_mg() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (6.5284100000e+02, 6.7598200000e-02),
                (9.8380500000e+01, 3.9177800000e-01),
                (2.1299600000e+01, 6.6666100000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.3372700000e+01, -1.1024600000e-01),
                (5.1995300000e+00, 1.8411900000e-01),
                (1.3150800000e+00, 8.9639900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (6.1134900000e-01, -3.6110100000e-01),
                (1.4184100000e-01, 1.2150500000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(4.6401100000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.3372700000e+01, 1.2101400000e-01),
                (5.1995300000e+00, 4.6281000000e-01),
                (1.3150800000e+00, 6.0690700000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (6.1134900000e-01, 2.4263300000e-02),
                (1.4184100000e-01, 9.8667300000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(4.6401100000e-02, 1.0000000000e+00)],
        ),
    ]
}

/// Al (Z=13) 3-21G basis
fn b321g_al() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (7.7573700000e+02, 6.6834700000e-02),
                (1.1695200000e+02, 3.8906100000e-01),
                (2.5332600000e+01, 6.6946800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.9479600000e+01, -1.0790200000e-01),
                (6.6331400000e+00, 1.4624500000e-01),
                (1.7267500000e+00, 9.2373000000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (9.4616000000e-01, -3.2032700000e-01),
                (2.0250600000e-01, 1.1841200000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(6.3908800000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.9479600000e+01, 1.1757400000e-01),
                (6.6331400000e+00, 4.6117400000e-01),
                (1.7267500000e+00, 6.0553500000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (9.4616000000e-01, 5.1938300000e-02),
                (2.0250600000e-01, 9.7266000000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(6.3908800000e-02, 1.0000000000e+00)],
        ),
    ]
}

/// Si (Z=14) 3-21G basis
fn b321g_si() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (9.1065500000e+02, 6.6082300000e-02),
                (1.3733600000e+02, 3.8622900000e-01),
                (2.9760100000e+01, 6.7238000000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.6671600000e+01, -1.0451100000e-01),
                (8.3172900000e+00, 1.0741000000e-01),
                (2.2164500000e+00, 9.5144600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.0791300000e+00, -3.7610800000e-01),
                (3.0242200000e-01, 1.2516500000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(9.3339200000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.6671600000e+01, 1.1335500000e-01),
                (8.3172900000e+00, 4.5757800000e-01),
                (2.2164500000e+00, 6.0742700000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.0791300000e+00, 6.7103000000e-02),
                (3.0242200000e-01, 9.5688300000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(9.3339200000e-02, 1.0000000000e+00)],
        ),
    ]
}

/// P (Z=15) 3-21G basis
fn b321g_p() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.0549000000e+03, 6.5541000000e-02),
                (1.5919500000e+02, 3.8403600000e-01),
                (3.4530400000e+01, 6.7454100000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.4286600000e+01, -1.0213000000e-01),
                (1.0101900000e+01, 8.1592000000e-02),
                (2.7399700000e+00, 9.6978800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.2186500000e+00, -3.7149500000e-01),
                (3.9554600000e-01, 1.2709900000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.2281100000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.4286600000e+01, 1.1085100000e-01),
                (1.0101900000e+01, 4.5649500000e-01),
                (2.7399700000e+00, 6.0693600000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.2186500000e+00, 9.1582000000e-02),
                (3.9554600000e-01, 9.3492400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.2281100000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// S (Z=16) 3-21G basis
fn b321g_s() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.2106200000e+03, 6.5007000000e-02),
                (1.8274700000e+02, 3.8204000000e-01),
                (3.9667300000e+01, 6.7654500000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (5.2223600000e+01, -1.0031000000e-01),
                (1.1962900000e+01, 6.5088000000e-02),
                (3.2891100000e+00, 9.8145500000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.2238400000e+00, -2.8608900000e-01),
                (4.5730300000e-01, 1.2280600000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.4226900000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (5.2223600000e+01, 1.0964600000e-01),
                (1.1962900000e+01, 4.5764900000e-01),
                (3.2891100000e+00, 6.0426100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.2238400000e+00, 1.6477700000e-01),
                (4.5730300000e-01, 8.7085500000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.4226900000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Cl (Z=17) 3-21G basis
fn b321g_cl() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.3764000000e+03, 6.4582700000e-02),
                (2.0785700000e+02, 3.8036300000e-01),
                (4.5155400000e+01, 6.7819000000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (6.0801400000e+01, -9.8763900000e-02),
                (1.3976500000e+01, 5.1133800000e-02),
                (3.8871000000e+00, 9.9133700000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.3529900000e+00, -2.2240100000e-01),
                (5.2695500000e-01, 1.1825200000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.6671400000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (6.0801400000e+01, 1.0859800000e-01),
                (1.3976500000e+01, 4.5868200000e-01),
                (3.8871000000e+00, 6.0196200000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.3529900000e+00, 2.1921600000e-01),
                (5.2695500000e-01, 8.2232100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.6671400000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Ar (Z=18) 3-21G basis
fn b321g_ar() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.5537100000e+03, 6.4170700000e-02),
                (2.3467800000e+02, 3.7879700000e-01),
                (5.1012100000e+01, 6.7975200000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (7.0045300000e+01, -9.7466100000e-02),
                (1.6147300000e+01, 3.9056900000e-02),
                (4.5349200000e+00, 9.9991600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.5420900000e+00, -1.7686600000e-01),
                (6.0726700000e-01, 1.1469000000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.9537300000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (7.0045300000e+01, 1.0761900000e-01),
                (1.6147300000e+01, 4.5957600000e-01),
                (4.5349200000e+00, 6.0004100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.5420900000e+00, 2.5568700000e-01),
                (6.0726700000e-01, 7.8984200000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.9537300000e-01, 1.0000000000e+00)],
        ),
    ]
}

// ====== 6-31G Na-Ar ======

/// Na (Z=11) 6-31G basis
fn b631g_na() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (9.9932000000e+03, 1.9377000000e-03),
                (1.4998900000e+03, 1.4807000000e-02),
                (3.4195100000e+02, 7.2706000000e-02),
                (9.4679700000e+01, 2.5262900000e-01),
                (2.9734500000e+01, 4.9324200000e-01),
                (1.0006300000e+01, 3.1316900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.5096300000e+02, -3.5421000000e-03),
                (3.5587800000e+01, -4.3959000000e-02),
                (1.1168300000e+01, -1.0975210000e-01),
                (3.9020100000e+00, 1.8739800000e-01),
                (1.3817700000e+00, 6.4669900000e-01),
                (4.6638200000e-01, 3.0605800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.9796600000e-01, -2.4850300000e-01),
                (8.4353000000e-02, -1.3170400000e-01),
                (6.6635000000e-02, 1.2335200000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(2.5954400000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.5096300000e+02, 5.0017000000e-03),
                (3.5587800000e+01, 3.5511000000e-02),
                (1.1168300000e+01, 1.4282500000e-01),
                (3.9020100000e+00, 3.3862000000e-01),
                (1.3817700000e+00, 4.5157900000e-01),
                (4.6638200000e-01, 2.7327100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.9796600000e-01, -2.3023000000e-02),
                (8.4353000000e-02, 9.5035900000e-01),
                (6.6635000000e-02, 5.9858000000e-02),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(2.5954400000e-02, 1.0000000000e+00)],
        ),
    ]
}

/// Mg (Z=12) 6-31G basis
fn b631g_mg() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.1722800000e+04, 1.9778000000e-03),
                (1.7599300000e+03, 1.5114000000e-02),
                (4.0084600000e+02, 7.3911000000e-02),
                (1.1280700000e+02, 2.4919100000e-01),
                (3.5999700000e+01, 4.8792800000e-01),
                (1.2182800000e+01, 3.1966200000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.8918000000e+02, -3.2372000000e-03),
                (4.5211900000e+01, -4.1008000000e-02),
                (1.4356300000e+01, -1.1260000000e-01),
                (5.1388600000e+00, 1.4863300000e-01),
                (1.9065200000e+00, 6.1649700000e-01),
                (7.0588700000e-01, 3.6482900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (9.2934000000e-01, -2.1229000000e-01),
                (2.6903500000e-01, -1.0798500000e-01),
                (1.1737900000e-01, 1.1758400000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(4.2106100000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.8918000000e+02, 4.9281000000e-03),
                (4.5211900000e+01, 3.4989000000e-02),
                (1.4356300000e+01, 1.4072500000e-01),
                (5.1388600000e+00, 3.3364200000e-01),
                (1.9065200000e+00, 4.4494000000e-01),
                (7.0588700000e-01, 2.6925400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (9.2934000000e-01, -2.2419000000e-02),
                (2.6903500000e-01, 1.9227000000e-01),
                (1.1737900000e-01, 8.4618100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(4.2106100000e-02, 1.0000000000e+00)],
        ),
    ]
}

/// Al (Z=13) 6-31G basis
fn b631g_al() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.3983100000e+04, 1.9426700000e-03),
                (2.0987500000e+03, 1.4859900000e-02),
                (4.7770500000e+02, 7.2849400000e-02),
                (1.3436000000e+02, 2.4683000000e-01),
                (4.2870900000e+01, 4.8725800000e-01),
                (1.4518900000e+01, 3.2349600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.3966800000e+02, -2.9261900000e-03),
                (5.7441900000e+01, -3.7408000000e-02),
                (1.8285900000e+01, -1.1448700000e-01),
                (6.5991400000e+00, 1.1563500000e-01),
                (2.4904900000e+00, 6.1259500000e-01),
                (9.4454000000e-01, 3.9379900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.2779000000e+00, -2.2760600000e-01),
                (3.9759000000e-01, 1.4458300000e-03),
                (1.6009500000e-01, 1.0927900000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(5.5657700000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.3966800000e+02, 4.6028500000e-03),
                (5.7441900000e+01, 3.3199000000e-02),
                (1.8285900000e+01, 1.3628200000e-01),
                (6.5991400000e+00, 3.3047600000e-01),
                (2.4904900000e+00, 4.4914600000e-01),
                (9.4454000000e-01, 2.6570400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.2779000000e+00, -1.7513000000e-02),
                (3.9759000000e-01, 2.4453300000e-01),
                (1.6009500000e-01, 8.0493400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(5.5657700000e-02, 1.0000000000e+00)],
        ),
    ]
}

/// Si (Z=14) 6-31G basis
fn b631g_si() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.6115900000e+04, 1.9594800000e-03),
                (2.4255800000e+03, 1.4928800000e-02),
                (5.5386700000e+02, 7.2847800000e-02),
                (1.5634000000e+02, 2.4613000000e-01),
                (5.0068300000e+01, 4.8591400000e-01),
                (1.7017800000e+01, 3.2500200000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.9271800000e+02, -2.7809400000e-03),
                (6.9873100000e+01, -3.5714600000e-02),
                (2.2336300000e+01, -1.1498500000e-01),
                (8.1503900000e+00, 9.3563400000e-02),
                (3.1345800000e+00, 6.0301700000e-01),
                (1.2254300000e+00, 4.1895900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.7273800000e+00, -2.4463000000e-01),
                (5.7292200000e-01, 4.3157200000e-03),
                (2.2219200000e-01, 1.0981800000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(7.7836900000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.9271800000e+02, 4.4382600000e-03),
                (6.9873100000e+01, 3.2667900000e-02),
                (2.2336300000e+01, 1.3472100000e-01),
                (8.1503900000e+00, 3.2867800000e-01),
                (3.1345800000e+00, 4.4964000000e-01),
                (1.2254300000e+00, 2.6137200000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.7273800000e+00, -1.7795100000e-02),
                (5.7292200000e-01, 2.5353900000e-01),
                (2.2219200000e-01, 8.0066900000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(7.7836900000e-02, 1.0000000000e+00)],
        ),
    ]
}

/// P (Z=15) 6-31G basis
fn b631g_p() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.9413300000e+04, 1.8516000000e-03),
                (2.9094200000e+03, 1.4206200000e-02),
                (6.6136400000e+02, 6.9999500000e-02),
                (1.8575900000e+02, 2.4007900000e-01),
                (5.9194300000e+01, 4.8476200000e-01),
                (2.0031000000e+01, 3.3520000000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.3947800000e+02, -2.7821700000e-03),
                (8.1010100000e+01, -3.6049900000e-02),
                (2.5878000000e+01, -1.1663100000e-01),
                (9.4522100000e+00, 9.6832800000e-02),
                (3.6656600000e+00, 6.1441800000e-01),
                (1.4674600000e+00, 4.0379800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.1562300000e+00, -2.5292300000e-01),
                (7.4899700000e-01, 3.2851700000e-02),
                (2.8314500000e-01, 1.0812500000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(9.9831700000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.3947800000e+02, 4.5646200000e-03),
                (8.1010100000e+01, 3.3693600000e-02),
                (2.5878000000e+01, 1.3975500000e-01),
                (9.4522100000e+00, 3.3936200000e-01),
                (3.6656600000e+00, 4.5092100000e-01),
                (1.4674600000e+00, 2.3858600000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.1562300000e+00, -1.7765300000e-02),
                (7.4899700000e-01, 2.7405800000e-01),
                (2.8314500000e-01, 7.8542100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(9.9831700000e-02, 1.0000000000e+00)],
        ),
    ]
}

/// S (Z=16) 6-31G basis
fn b631g_s() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.1917100000e+04, 1.8690000000e-03),
                (3.3014900000e+03, 1.4230000000e-02),
                (7.5414600000e+02, 6.9696000000e-02),
                (2.1271100000e+02, 2.3848700000e-01),
                (6.7989600000e+01, 4.8330700000e-01),
                (2.3051500000e+01, 3.3807400000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.2373500000e+02, -2.3767000000e-03),
                (1.0071000000e+02, -3.1693000000e-02),
                (3.2159900000e+01, -1.1331700000e-01),
                (1.1807900000e+01, 5.6090000000e-02),
                (4.6311000000e+00, 5.9225500000e-01),
                (1.8702500000e+00, 4.5500600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.6158400000e+00, -2.5037400000e-01),
                (9.2216700000e-01, 6.6957000000e-02),
                (3.4128700000e-01, 1.0545100000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.1716700000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.2373500000e+02, 4.0610000000e-03),
                (1.0071000000e+02, 3.0681000000e-02),
                (3.2159900000e+01, 1.3045200000e-01),
                (1.1807900000e+01, 3.2720500000e-01),
                (4.6311000000e+00, 4.5285100000e-01),
                (1.8702500000e+00, 2.5604200000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.6158400000e+00, -1.4511000000e-02),
                (9.2216700000e-01, 3.1026300000e-01),
                (3.4128700000e-01, 7.5448300000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.1716700000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Cl (Z=17) 6-31G basis
fn b631g_cl() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.5180100000e+04, 1.8330000000e-03),
                (3.7803500000e+03, 1.4034000000e-02),
                (8.6047400000e+02, 6.9097000000e-02),
                (2.4214500000e+02, 2.3745200000e-01),
                (7.7334900000e+01, 4.8303400000e-01),
                (2.6247000000e+01, 3.3985600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.9176500000e+02, -2.2974000000e-03),
                (1.1698400000e+02, -3.0714000000e-02),
                (3.7415300000e+01, -1.1252800000e-01),
                (1.3783400000e+01, 4.5016000000e-02),
                (5.4521500000e+00, 5.8935300000e-01),
                (2.2258800000e+00, 4.6520600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.1864900000e+00, -2.5183000000e-01),
                (1.1442700000e+00, 6.1589000000e-02),
                (4.2037700000e-01, 1.0601800000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.4265700000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.9176500000e+02, 3.9894000000e-03),
                (1.1698400000e+02, 3.0318000000e-02),
                (3.7415300000e+01, 1.2988000000e-01),
                (1.3783400000e+01, 3.2795100000e-01),
                (5.4521500000e+00, 4.5352700000e-01),
                (2.2258800000e+00, 2.5215400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.1864900000e+00, -1.4299000000e-02),
                (1.1442700000e+00, 3.2357200000e-01),
                (4.2037700000e-01, 7.4350700000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.4265700000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Ar (Z=18) 6-31G basis
fn b631g_ar() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.8348300000e+04, 1.8252600000e-03),
                (4.2576200000e+03, 1.3968600000e-02),
                (9.6985700000e+02, 6.8707300000e-02),
                (2.7326300000e+02, 2.3620400000e-01),
                (8.7369500000e+01, 4.8221400000e-01),
                (2.9686700000e+01, 3.4204300000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (5.7589100000e+02, -2.1597200000e-03),
                (1.3681600000e+02, -2.9077500000e-02),
                (4.3809800000e+01, -1.1082700000e-01),
                (1.6209400000e+01, 2.7699900000e-02),
                (6.4608400000e+00, 5.7761300000e-01),
                (2.6511400000e+00, 4.8868800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.8602800000e+00, -2.5559200000e-01),
                (1.4137300000e+00, 3.7806600000e-02),
                (5.1664600000e-01, 1.0805600000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.7388800000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (5.7589100000e+02, 3.8066500000e-03),
                (1.3681600000e+02, 2.9230500000e-02),
                (4.3809800000e+01, 1.2646700000e-01),
                (1.6209400000e+01, 3.2351000000e-01),
                (6.4608400000e+00, 4.5489600000e-01),
                (2.6511400000e+00, 2.5663000000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.8602800000e+00, -1.5919700000e-02),
                (1.4137300000e+00, 3.2464600000e-01),
                (5.1664600000e-01, 7.4399000000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.7388800000e-01, 1.0000000000e+00)],
        ),
    ]
}

// ====== 6-31G* Na-Ar ======

/// Na (Z=11) 6-31G* basis
fn b631gs_na() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (9.9932000000e+03, 1.9377000000e-03),
                (1.4998900000e+03, 1.4807000000e-02),
                (3.4195100000e+02, 7.2706000000e-02),
                (9.4679700000e+01, 2.5262900000e-01),
                (2.9734500000e+01, 4.9324200000e-01),
                (1.0006300000e+01, 3.1316900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.5096300000e+02, -3.5421000000e-03),
                (3.5587800000e+01, -4.3959000000e-02),
                (1.1168300000e+01, -1.0975210000e-01),
                (3.9020100000e+00, 1.8739800000e-01),
                (1.3817700000e+00, 6.4669900000e-01),
                (4.6638200000e-01, 3.0605800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.9796600000e-01, -2.4850300000e-01),
                (8.4353000000e-02, -1.3170400000e-01),
                (6.6635000000e-02, 1.2335200000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(2.5954400000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.5096300000e+02, 5.0017000000e-03),
                (3.5587800000e+01, 3.5511000000e-02),
                (1.1168300000e+01, 1.4282500000e-01),
                (3.9020100000e+00, 3.3862000000e-01),
                (1.3817700000e+00, 4.5157900000e-01),
                (4.6638200000e-01, 2.7327100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.9796600000e-01, -2.3023000000e-02),
                (8.4353000000e-02, 9.5035900000e-01),
                (6.6635000000e-02, 5.9858000000e-02),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(2.5954400000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(1.7500000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Mg (Z=12) 6-31G* basis
fn b631gs_mg() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.1722800000e+04, 1.9778000000e-03),
                (1.7599300000e+03, 1.5114000000e-02),
                (4.0084600000e+02, 7.3911000000e-02),
                (1.1280700000e+02, 2.4919100000e-01),
                (3.5999700000e+01, 4.8792800000e-01),
                (1.2182800000e+01, 3.1966200000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.8918000000e+02, -3.2372000000e-03),
                (4.5211900000e+01, -4.1008000000e-02),
                (1.4356300000e+01, -1.1260000000e-01),
                (5.1388600000e+00, 1.4863300000e-01),
                (1.9065200000e+00, 6.1649700000e-01),
                (7.0588700000e-01, 3.6482900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (9.2934000000e-01, -2.1229000000e-01),
                (2.6903500000e-01, -1.0798500000e-01),
                (1.1737900000e-01, 1.1758400000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(4.2106100000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.8918000000e+02, 4.9281000000e-03),
                (4.5211900000e+01, 3.4989000000e-02),
                (1.4356300000e+01, 1.4072500000e-01),
                (5.1388600000e+00, 3.3364200000e-01),
                (1.9065200000e+00, 4.4494000000e-01),
                (7.0588700000e-01, 2.6925400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (9.2934000000e-01, -2.2419000000e-02),
                (2.6903500000e-01, 1.9227000000e-01),
                (1.1737900000e-01, 8.4618100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(4.2106100000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(1.7500000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Al (Z=13) 6-31G* basis
fn b631gs_al() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.3983100000e+04, 1.9426700000e-03),
                (2.0987500000e+03, 1.4859900000e-02),
                (4.7770500000e+02, 7.2849400000e-02),
                (1.3436000000e+02, 2.4683000000e-01),
                (4.2870900000e+01, 4.8725800000e-01),
                (1.4518900000e+01, 3.2349600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.3966800000e+02, -2.9261900000e-03),
                (5.7441900000e+01, -3.7408000000e-02),
                (1.8285900000e+01, -1.1448700000e-01),
                (6.5991400000e+00, 1.1563500000e-01),
                (2.4904900000e+00, 6.1259500000e-01),
                (9.4454000000e-01, 3.9379900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.2779000000e+00, -2.2760600000e-01),
                (3.9759000000e-01, 1.4458300000e-03),
                (1.6009500000e-01, 1.0927900000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(5.5657700000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.3966800000e+02, 4.6028500000e-03),
                (5.7441900000e+01, 3.3199000000e-02),
                (1.8285900000e+01, 1.3628200000e-01),
                (6.5991400000e+00, 3.3047600000e-01),
                (2.4904900000e+00, 4.4914600000e-01),
                (9.4454000000e-01, 2.6570400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.2779000000e+00, -1.7513000000e-02),
                (3.9759000000e-01, 2.4453300000e-01),
                (1.6009500000e-01, 8.0493400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(5.5657700000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(3.2500000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Si (Z=14) 6-31G* basis
fn b631gs_si() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.6115900000e+04, 1.9594800000e-03),
                (2.4255800000e+03, 1.4928800000e-02),
                (5.5386700000e+02, 7.2847800000e-02),
                (1.5634000000e+02, 2.4613000000e-01),
                (5.0068300000e+01, 4.8591400000e-01),
                (1.7017800000e+01, 3.2500200000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.9271800000e+02, -2.7809400000e-03),
                (6.9873100000e+01, -3.5714600000e-02),
                (2.2336300000e+01, -1.1498500000e-01),
                (8.1503900000e+00, 9.3563400000e-02),
                (3.1345800000e+00, 6.0301700000e-01),
                (1.2254300000e+00, 4.1895900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.7273800000e+00, -2.4463000000e-01),
                (5.7292200000e-01, 4.3157200000e-03),
                (2.2219200000e-01, 1.0981800000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(7.7836900000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.9271800000e+02, 4.4382600000e-03),
                (6.9873100000e+01, 3.2667900000e-02),
                (2.2336300000e+01, 1.3472100000e-01),
                (8.1503900000e+00, 3.2867800000e-01),
                (3.1345800000e+00, 4.4964000000e-01),
                (1.2254300000e+00, 2.6137200000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.7273800000e+00, -1.7795100000e-02),
                (5.7292200000e-01, 2.5353900000e-01),
                (2.2219200000e-01, 8.0066900000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(7.7836900000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(4.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// P (Z=15) 6-31G* basis
fn b631gs_p() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.9413300000e+04, 1.8516000000e-03),
                (2.9094200000e+03, 1.4206200000e-02),
                (6.6136400000e+02, 6.9999500000e-02),
                (1.8575900000e+02, 2.4007900000e-01),
                (5.9194300000e+01, 4.8476200000e-01),
                (2.0031000000e+01, 3.3520000000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.3947800000e+02, -2.7821700000e-03),
                (8.1010100000e+01, -3.6049900000e-02),
                (2.5878000000e+01, -1.1663100000e-01),
                (9.4522100000e+00, 9.6832800000e-02),
                (3.6656600000e+00, 6.1441800000e-01),
                (1.4674600000e+00, 4.0379800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.1562300000e+00, -2.5292300000e-01),
                (7.4899700000e-01, 3.2851700000e-02),
                (2.8314500000e-01, 1.0812500000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(9.9831700000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.3947800000e+02, 4.5646200000e-03),
                (8.1010100000e+01, 3.3693600000e-02),
                (2.5878000000e+01, 1.3975500000e-01),
                (9.4522100000e+00, 3.3936200000e-01),
                (3.6656600000e+00, 4.5092100000e-01),
                (1.4674600000e+00, 2.3858600000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.1562300000e+00, -1.7765300000e-02),
                (7.4899700000e-01, 2.7405800000e-01),
                (2.8314500000e-01, 7.8542100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(9.9831700000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(5.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// S (Z=16) 6-31G* basis
fn b631gs_s() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.1917100000e+04, 1.8690000000e-03),
                (3.3014900000e+03, 1.4230000000e-02),
                (7.5414600000e+02, 6.9696000000e-02),
                (2.1271100000e+02, 2.3848700000e-01),
                (6.7989600000e+01, 4.8330700000e-01),
                (2.3051500000e+01, 3.3807400000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.2373500000e+02, -2.3767000000e-03),
                (1.0071000000e+02, -3.1693000000e-02),
                (3.2159900000e+01, -1.1331700000e-01),
                (1.1807900000e+01, 5.6090000000e-02),
                (4.6311000000e+00, 5.9225500000e-01),
                (1.8702500000e+00, 4.5500600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.6158400000e+00, -2.5037400000e-01),
                (9.2216700000e-01, 6.6957000000e-02),
                (3.4128700000e-01, 1.0545100000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.1716700000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.2373500000e+02, 4.0610000000e-03),
                (1.0071000000e+02, 3.0681000000e-02),
                (3.2159900000e+01, 1.3045200000e-01),
                (1.1807900000e+01, 3.2720500000e-01),
                (4.6311000000e+00, 4.5285100000e-01),
                (1.8702500000e+00, 2.5604200000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.6158400000e+00, -1.4511000000e-02),
                (9.2216700000e-01, 3.1026300000e-01),
                (3.4128700000e-01, 7.5448300000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.1716700000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(6.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Cl (Z=17) 6-31G* basis
fn b631gs_cl() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.5180100000e+04, 1.8330000000e-03),
                (3.7803500000e+03, 1.4034000000e-02),
                (8.6047400000e+02, 6.9097000000e-02),
                (2.4214500000e+02, 2.3745200000e-01),
                (7.7334900000e+01, 4.8303400000e-01),
                (2.6247000000e+01, 3.3985600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.9176500000e+02, -2.2974000000e-03),
                (1.1698400000e+02, -3.0714000000e-02),
                (3.7415300000e+01, -1.1252800000e-01),
                (1.3783400000e+01, 4.5016000000e-02),
                (5.4521500000e+00, 5.8935300000e-01),
                (2.2258800000e+00, 4.6520600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.1864900000e+00, -2.5183000000e-01),
                (1.1442700000e+00, 6.1589000000e-02),
                (4.2037700000e-01, 1.0601800000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.4265700000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.9176500000e+02, 3.9894000000e-03),
                (1.1698400000e+02, 3.0318000000e-02),
                (3.7415300000e+01, 1.2988000000e-01),
                (1.3783400000e+01, 3.2795100000e-01),
                (5.4521500000e+00, 4.5352700000e-01),
                (2.2258800000e+00, 2.5215400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.1864900000e+00, -1.4299000000e-02),
                (1.1442700000e+00, 3.2357200000e-01),
                (4.2037700000e-01, 7.4350700000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.4265700000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(7.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Ar (Z=18) 6-31G* basis
fn b631gs_ar() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.8348300000e+04, 1.8252600000e-03),
                (4.2576200000e+03, 1.3968600000e-02),
                (9.6985700000e+02, 6.8707300000e-02),
                (2.7326300000e+02, 2.3620400000e-01),
                (8.7369500000e+01, 4.8221400000e-01),
                (2.9686700000e+01, 3.4204300000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (5.7589100000e+02, -2.1597200000e-03),
                (1.3681600000e+02, -2.9077500000e-02),
                (4.3809800000e+01, -1.1082700000e-01),
                (1.6209400000e+01, 2.7699900000e-02),
                (6.4608400000e+00, 5.7761300000e-01),
                (2.6511400000e+00, 4.8868800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.8602800000e+00, -2.5559200000e-01),
                (1.4137300000e+00, 3.7806600000e-02),
                (5.1664600000e-01, 1.0805600000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.7388800000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (5.7589100000e+02, 3.8066500000e-03),
                (1.3681600000e+02, 2.9230500000e-02),
                (4.3809800000e+01, 1.2646700000e-01),
                (1.6209400000e+01, 3.2351000000e-01),
                (6.4608400000e+00, 4.5489600000e-01),
                (2.6511400000e+00, 2.5663000000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.8602800000e+00, -1.5919700000e-02),
                (1.4137300000e+00, 3.2464600000e-01),
                (5.1664600000e-01, 7.4399000000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.7388800000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(8.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

// ====== 6-31+G* Na-Ar ======

/// Na (Z=11) 6-31+G* basis
fn b631pgs_na() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (9.9932000000e+03, 1.9377000000e-03),
                (1.4998900000e+03, 1.4807000000e-02),
                (3.4195100000e+02, 7.2706000000e-02),
                (9.4679700000e+01, 2.5262900000e-01),
                (2.9734500000e+01, 4.9324200000e-01),
                (1.0006300000e+01, 3.1316900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.5096300000e+02, -3.5421000000e-03),
                (3.5587800000e+01, -4.3959000000e-02),
                (1.1168300000e+01, -1.0975210000e-01),
                (3.9020100000e+00, 1.8739800000e-01),
                (1.3817700000e+00, 6.4669900000e-01),
                (4.6638200000e-01, 3.0605800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.9796600000e-01, -2.4850300000e-01),
                (8.4353000000e-02, -1.3170400000e-01),
                (6.6635000000e-02, 1.2335200000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(2.5954400000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::S,
            vec![(7.6000000000e-03, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.5096300000e+02, 5.0017000000e-03),
                (3.5587800000e+01, 3.5511000000e-02),
                (1.1168300000e+01, 1.4282500000e-01),
                (3.9020100000e+00, 3.3862000000e-01),
                (1.3817700000e+00, 4.5157900000e-01),
                (4.6638200000e-01, 2.7327100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.9796600000e-01, -2.3023000000e-02),
                (8.4353000000e-02, 9.5035900000e-01),
                (6.6635000000e-02, 5.9858000000e-02),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(2.5954400000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![(7.6000000000e-03, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(1.7500000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Mg (Z=12) 6-31+G* basis
fn b631pgs_mg() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.1722800000e+04, 1.9778000000e-03),
                (1.7599300000e+03, 1.5114000000e-02),
                (4.0084600000e+02, 7.3911000000e-02),
                (1.1280700000e+02, 2.4919100000e-01),
                (3.5999700000e+01, 4.8792800000e-01),
                (1.2182800000e+01, 3.1966200000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.8918000000e+02, -3.2372000000e-03),
                (4.5211900000e+01, -4.1008000000e-02),
                (1.4356300000e+01, -1.1260000000e-01),
                (5.1388600000e+00, 1.4863300000e-01),
                (1.9065200000e+00, 6.1649700000e-01),
                (7.0588700000e-01, 3.6482900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (9.2934000000e-01, -2.1229000000e-01),
                (2.6903500000e-01, -1.0798500000e-01),
                (1.1737900000e-01, 1.1758400000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(4.2106100000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::S,
            vec![(1.4600000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.8918000000e+02, 4.9281000000e-03),
                (4.5211900000e+01, 3.4989000000e-02),
                (1.4356300000e+01, 1.4072500000e-01),
                (5.1388600000e+00, 3.3364200000e-01),
                (1.9065200000e+00, 4.4494000000e-01),
                (7.0588700000e-01, 2.6925400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (9.2934000000e-01, -2.2419000000e-02),
                (2.6903500000e-01, 1.9227000000e-01),
                (1.1737900000e-01, 8.4618100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(4.2106100000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![(1.4600000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(1.7500000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Al (Z=13) 6-31+G* basis
fn b631pgs_al() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.3983100000e+04, 1.9426700000e-03),
                (2.0987500000e+03, 1.4859900000e-02),
                (4.7770500000e+02, 7.2849400000e-02),
                (1.3436000000e+02, 2.4683000000e-01),
                (4.2870900000e+01, 4.8725800000e-01),
                (1.4518900000e+01, 3.2349600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.3966800000e+02, -2.9261900000e-03),
                (5.7441900000e+01, -3.7408000000e-02),
                (1.8285900000e+01, -1.1448700000e-01),
                (6.5991400000e+00, 1.1563500000e-01),
                (2.4904900000e+00, 6.1259500000e-01),
                (9.4454000000e-01, 3.9379900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.2779000000e+00, -2.2760600000e-01),
                (3.9759000000e-01, 1.4458300000e-03),
                (1.6009500000e-01, 1.0927900000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(5.5657700000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::S,
            vec![(3.1800000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.3966800000e+02, 4.6028500000e-03),
                (5.7441900000e+01, 3.3199000000e-02),
                (1.8285900000e+01, 1.3628200000e-01),
                (6.5991400000e+00, 3.3047600000e-01),
                (2.4904900000e+00, 4.4914600000e-01),
                (9.4454000000e-01, 2.6570400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.2779000000e+00, -1.7513000000e-02),
                (3.9759000000e-01, 2.4453300000e-01),
                (1.6009500000e-01, 8.0493400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(5.5657700000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![(3.1800000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(3.2500000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Si (Z=14) 6-31+G* basis
fn b631pgs_si() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.6115900000e+04, 1.9594800000e-03),
                (2.4255800000e+03, 1.4928800000e-02),
                (5.5386700000e+02, 7.2847800000e-02),
                (1.5634000000e+02, 2.4613000000e-01),
                (5.0068300000e+01, 4.8591400000e-01),
                (1.7017800000e+01, 3.2500200000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.9271800000e+02, -2.7809400000e-03),
                (6.9873100000e+01, -3.5714600000e-02),
                (2.2336300000e+01, -1.1498500000e-01),
                (8.1503900000e+00, 9.3563400000e-02),
                (3.1345800000e+00, 6.0301700000e-01),
                (1.2254300000e+00, 4.1895900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.7273800000e+00, -2.4463000000e-01),
                (5.7292200000e-01, 4.3157200000e-03),
                (2.2219200000e-01, 1.0981800000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(7.7836900000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::S,
            vec![(3.3100000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.9271800000e+02, 4.4382600000e-03),
                (6.9873100000e+01, 3.2667900000e-02),
                (2.2336300000e+01, 1.3472100000e-01),
                (8.1503900000e+00, 3.2867800000e-01),
                (3.1345800000e+00, 4.4964000000e-01),
                (1.2254300000e+00, 2.6137200000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.7273800000e+00, -1.7795100000e-02),
                (5.7292200000e-01, 2.5353900000e-01),
                (2.2219200000e-01, 8.0066900000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(7.7836900000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![(3.3100000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(4.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// P (Z=15) 6-31+G* basis
fn b631pgs_p() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.9413300000e+04, 1.8516000000e-03),
                (2.9094200000e+03, 1.4206200000e-02),
                (6.6136400000e+02, 6.9999500000e-02),
                (1.8575900000e+02, 2.4007900000e-01),
                (5.9194300000e+01, 4.8476200000e-01),
                (2.0031000000e+01, 3.3520000000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.3947800000e+02, -2.7821700000e-03),
                (8.1010100000e+01, -3.6049900000e-02),
                (2.5878000000e+01, -1.1663100000e-01),
                (9.4522100000e+00, 9.6832800000e-02),
                (3.6656600000e+00, 6.1441800000e-01),
                (1.4674600000e+00, 4.0379800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.1562300000e+00, -2.5292300000e-01),
                (7.4899700000e-01, 3.2851700000e-02),
                (2.8314500000e-01, 1.0812500000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(9.9831700000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::S,
            vec![(3.4800000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.3947800000e+02, 4.5646200000e-03),
                (8.1010100000e+01, 3.3693600000e-02),
                (2.5878000000e+01, 1.3975500000e-01),
                (9.4522100000e+00, 3.3936200000e-01),
                (3.6656600000e+00, 4.5092100000e-01),
                (1.4674600000e+00, 2.3858600000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.1562300000e+00, -1.7765300000e-02),
                (7.4899700000e-01, 2.7405800000e-01),
                (2.8314500000e-01, 7.8542100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(9.9831700000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![(3.4800000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(5.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// S (Z=16) 6-31+G* basis
fn b631pgs_s() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.1917100000e+04, 1.8690000000e-03),
                (3.3014900000e+03, 1.4230000000e-02),
                (7.5414600000e+02, 6.9696000000e-02),
                (2.1271100000e+02, 2.3848700000e-01),
                (6.7989600000e+01, 4.8330700000e-01),
                (2.3051500000e+01, 3.3807400000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.2373500000e+02, -2.3767000000e-03),
                (1.0071000000e+02, -3.1693000000e-02),
                (3.2159900000e+01, -1.1331700000e-01),
                (1.1807900000e+01, 5.6090000000e-02),
                (4.6311000000e+00, 5.9225500000e-01),
                (1.8702500000e+00, 4.5500600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.6158400000e+00, -2.5037400000e-01),
                (9.2216700000e-01, 6.6957000000e-02),
                (3.4128700000e-01, 1.0545100000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.1716700000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::S,
            vec![(4.0500000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.2373500000e+02, 4.0610000000e-03),
                (1.0071000000e+02, 3.0681000000e-02),
                (3.2159900000e+01, 1.3045200000e-01),
                (1.1807900000e+01, 3.2720500000e-01),
                (4.6311000000e+00, 4.5285100000e-01),
                (1.8702500000e+00, 2.5604200000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.6158400000e+00, -1.4511000000e-02),
                (9.2216700000e-01, 3.1026300000e-01),
                (3.4128700000e-01, 7.5448300000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.1716700000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![(4.0500000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(6.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Cl (Z=17) 6-31+G* basis
fn b631pgs_cl() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.5180100000e+04, 1.8330000000e-03),
                (3.7803500000e+03, 1.4034000000e-02),
                (8.6047400000e+02, 6.9097000000e-02),
                (2.4214500000e+02, 2.3745200000e-01),
                (7.7334900000e+01, 4.8303400000e-01),
                (2.6247000000e+01, 3.3985600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.9176500000e+02, -2.2974000000e-03),
                (1.1698400000e+02, -3.0714000000e-02),
                (3.7415300000e+01, -1.1252800000e-01),
                (1.3783400000e+01, 4.5016000000e-02),
                (5.4521500000e+00, 5.8935300000e-01),
                (2.2258800000e+00, 4.6520600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.1864900000e+00, -2.5183000000e-01),
                (1.1442700000e+00, 6.1589000000e-02),
                (4.2037700000e-01, 1.0601800000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.4265700000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::S,
            vec![(4.8300000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.9176500000e+02, 3.9894000000e-03),
                (1.1698400000e+02, 3.0318000000e-02),
                (3.7415300000e+01, 1.2988000000e-01),
                (1.3783400000e+01, 3.2795100000e-01),
                (5.4521500000e+00, 4.5352700000e-01),
                (2.2258800000e+00, 2.5215400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.1864900000e+00, -1.4299000000e-02),
                (1.1442700000e+00, 3.2357200000e-01),
                (4.2037700000e-01, 7.4350700000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.4265700000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![(4.8300000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(7.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Ar (Z=18) 6-31+G* basis
fn b631pgs_ar() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.8348300000e+04, 1.8252600000e-03),
                (4.2576200000e+03, 1.3968600000e-02),
                (9.6985700000e+02, 6.8707300000e-02),
                (2.7326300000e+02, 2.3620400000e-01),
                (8.7369500000e+01, 4.8221400000e-01),
                (2.9686700000e+01, 3.4204300000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (5.7589100000e+02, -2.1597200000e-03),
                (1.3681600000e+02, -2.9077500000e-02),
                (4.3809800000e+01, -1.1082700000e-01),
                (1.6209400000e+01, 2.7699900000e-02),
                (6.4608400000e+00, 5.7761300000e-01),
                (2.6511400000e+00, 4.8868800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.8602800000e+00, -2.5559200000e-01),
                (1.4137300000e+00, 3.7806600000e-02),
                (5.1664600000e-01, 1.0805600000e+00),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.7388800000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::S,
            vec![(6.0000000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (5.7589100000e+02, 3.8066500000e-03),
                (1.3681600000e+02, 2.9230500000e-02),
                (4.3809800000e+01, 1.2646700000e-01),
                (1.6209400000e+01, 3.2351000000e-01),
                (6.4608400000e+00, 4.5489600000e-01),
                (2.6511400000e+00, 2.5663000000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.8602800000e+00, -1.5919700000e-02),
                (1.4137300000e+00, 3.2464600000e-01),
                (5.1664600000e-01, 7.4399000000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.7388800000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![(6.0000000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(8.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

// =============================================================================
// cc-pVDZ Basis Set Data (Dunning correlation-consistent polarized valence double-zeta)
// Reference: Dunning, T.H. Jr., J. Chem. Phys. 90, 1007 (1989)
// Data extracted from PySCF 2.11.0
//
// IMPORTANT: cc-pVDZ is conventionally used with spherical d-functions (5d).
// IQCP currently uses Cartesian d-functions (6d) for consistency with the
// existing Pople basis set implementations. This means basis function counts
// differ from standard cc-pVDZ tabulations:
//   - H, He: 5 BFs (no d-functions, so no difference)
//   - Li-Ne: 15 BFs (Cartesian) vs 14 BFs (spherical)
//   - Na-Ar: 19 BFs (Cartesian) vs 18 BFs (spherical)
//
// General contractions from the original cc-pVDZ definition are "unrolled"
// into separate segmented contractions sharing the same primitive exponents.
// For example, Carbon's (9s4p1d) -> [3s2p1d] general contraction becomes:
//   - S shell 1: 8 primitives with contraction coefficients from column 1
//   - S shell 2: 8 primitives with contraction coefficients from column 2
//   - S shell 3: 1 primitive (uncontracted)
//   - P shell 1: 3 primitives (contracted)
//   - P shell 2: 1 primitive (uncontracted)
//   - D shell: 1 primitive (uncontracted)
// =============================================================================

fn get_ccpvdz(z: u8) -> Result<Vec<ShellData>, BasisError> {
    match z {
        1 => Ok(ccpvdz_h()),
        2 => Ok(ccpvdz_he()),
        3 => Ok(ccpvdz_li()),
        4 => Ok(ccpvdz_be()),
        5 => Ok(ccpvdz_b()),
        6 => Ok(ccpvdz_c()),
        7 => Ok(ccpvdz_n()),
        8 => Ok(ccpvdz_o()),
        9 => Ok(ccpvdz_f()),
        10 => Ok(ccpvdz_ne()),
        11 => Ok(ccpvdz_na()),
        12 => Ok(ccpvdz_mg()),
        13 => Ok(ccpvdz_al()),
        14 => Ok(ccpvdz_si()),
        15 => Ok(ccpvdz_p()),
        16 => Ok(ccpvdz_s()),
        17 => Ok(ccpvdz_cl()),
        18 => Ok(ccpvdz_ar()),
        _ => Err(BasisError::ElementNotInBasis(z, "cc-pvdz".to_string())),
    }
}

/// Hydrogen cc-pVDZ basis
/// (4s1p) -> [2s1p]: 5 Cartesian basis functions
fn ccpvdz_h() -> Vec<ShellData> {
    vec![
        // 3-primitive contracted s shell
        (
            AngularMomentum::S,
            vec![
                (1.3010000000e+01, 1.9685000000e-02),
                (1.9620000000e+00, 1.3797700000e-01),
                (4.4460000000e-01, 4.7814800000e-01),
            ],
        ),
        // uncontracted s shell
        (
            AngularMomentum::S,
            vec![(1.2200000000e-01, 1.0000000000e+00)],
        ),
        // polarization p shell
        (
            AngularMomentum::P,
            vec![(7.2700000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Helium cc-pVDZ basis
/// (4s1p) -> [2s1p]: 5 Cartesian basis functions
fn ccpvdz_he() -> Vec<ShellData> {
    vec![
        // 3-primitive contracted s shell
        (
            AngularMomentum::S,
            vec![
                (3.8360000000e+01, 2.3809000000e-02),
                (5.7700000000e+00, 1.5489100000e-01),
                (1.2400000000e+00, 4.6998700000e-01),
            ],
        ),
        // uncontracted s shell
        (
            AngularMomentum::S,
            vec![(2.9760000000e-01, 1.0000000000e+00)],
        ),
        // polarization p shell
        (
            AngularMomentum::P,
            vec![(1.2750000000e+00, 1.0000000000e+00)],
        ),
    ]
}

/// Lithium cc-pVDZ basis
/// (9s4p1d) -> [3s2p1d]: 15 Cartesian basis functions
/// General contraction: 8-primitive s-block produces 2 contracted s functions
fn ccpvdz_li() -> Vec<ShellData> {
    vec![
        // 8-primitive contracted s shell (contraction 1)
        (
            AngularMomentum::S,
            vec![
                (1.4690000000e+03, 7.6600000000e-04),
                (2.2050000000e+02, 5.8920000000e-03),
                (5.0260000000e+01, 2.9671000000e-02),
                (1.4240000000e+01, 1.0918000000e-01),
                (4.5810000000e+00, 2.8278900000e-01),
                (1.5800000000e+00, 4.5312300000e-01),
                (5.6400000000e-01, 2.7477400000e-01),
                (7.3450000000e-02, 9.7510000000e-03),
            ],
        ),
        // 8-primitive contracted s shell (contraction 2)
        (
            AngularMomentum::S,
            vec![
                (1.4690000000e+03, -1.2000000000e-04),
                (2.2050000000e+02, -9.2300000000e-04),
                (5.0260000000e+01, -4.6890000000e-03),
                (1.4240000000e+01, -1.7682000000e-02),
                (4.5810000000e+00, -4.8902000000e-02),
                (1.5800000000e+00, -9.6009000000e-02),
                (5.6400000000e-01, -1.3638000000e-01),
                (7.3450000000e-02, 5.7510200000e-01),
            ],
        ),
        // uncontracted s shell
        (
            AngularMomentum::S,
            vec![(2.8050000000e-02, 1.0000000000e+00)],
        ),
        // 3-primitive contracted p shell
        (
            AngularMomentum::P,
            vec![
                (1.5340000000e+00, 2.2784000000e-02),
                (2.7490000000e-01, 1.3910700000e-01),
                (7.3620000000e-02, 5.0037500000e-01),
            ],
        ),
        // uncontracted p shell
        (
            AngularMomentum::P,
            vec![(2.4030000000e-02, 1.0000000000e+00)],
        ),
        // d polarization shell
        (
            AngularMomentum::D,
            vec![(1.2390000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Beryllium cc-pVDZ basis
/// (9s4p1d) -> [3s2p1d]: 15 Cartesian basis functions
fn ccpvdz_be() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (2.9400000000e+03, 6.8000000000e-04),
                (4.4120000000e+02, 5.2360000000e-03),
                (1.0050000000e+02, 2.6606000000e-02),
                (2.8430000000e+01, 9.9993000000e-02),
                (9.1690000000e+00, 2.6970200000e-01),
                (3.1960000000e+00, 4.5146900000e-01),
                (1.1590000000e+00, 2.9507400000e-01),
                (1.8110000000e-01, 1.2587000000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (2.9400000000e+03, -1.2300000000e-04),
                (4.4120000000e+02, -9.6600000000e-04),
                (1.0050000000e+02, -4.8310000000e-03),
                (2.8430000000e+01, -1.9314000000e-02),
                (9.1690000000e+00, -5.3280000000e-02),
                (3.1960000000e+00, -1.2072300000e-01),
                (1.1590000000e+00, -1.3343500000e-01),
                (1.8110000000e-01, 5.3076700000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(5.8900000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.6190000000e+00, 2.9111000000e-02),
                (7.1100000000e-01, 1.6936500000e-01),
                (1.9510000000e-01, 5.1345800000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(6.0180000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(2.3800000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Boron cc-pVDZ basis
/// (9s4p1d) -> [3s2p1d]: 15 Cartesian basis functions
fn ccpvdz_b() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (4.5700000000e+03, 6.9600000000e-04),
                (6.8590000000e+02, 5.3530000000e-03),
                (1.5650000000e+02, 2.7134000000e-02),
                (4.4470000000e+01, 1.0138000000e-01),
                (1.4480000000e+01, 2.7205500000e-01),
                (5.1310000000e+00, 4.4840300000e-01),
                (1.8980000000e+00, 2.9012300000e-01),
                (3.3290000000e-01, 1.4322000000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.5700000000e+03, -1.3900000000e-04),
                (6.8590000000e+02, -1.0970000000e-03),
                (1.5650000000e+02, -5.4440000000e-03),
                (4.4470000000e+01, -2.1916000000e-02),
                (1.4480000000e+01, -5.9751000000e-02),
                (5.1310000000e+00, -1.3873200000e-01),
                (1.8980000000e+00, -1.3148200000e-01),
                (3.3290000000e-01, 5.3952600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.0430000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (6.0010000000e+00, 3.5481000000e-02),
                (1.2410000000e+00, 1.9807200000e-01),
                (3.3640000000e-01, 5.0523000000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(9.5380000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(3.4300000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Carbon cc-pVDZ basis
/// (9s4p1d) -> [3s2p1d]: 15 Cartesian basis functions
fn ccpvdz_c() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (6.6650000000e+03, 6.9200000000e-04),
                (1.0000000000e+03, 5.3290000000e-03),
                (2.2800000000e+02, 2.7077000000e-02),
                (6.4710000000e+01, 1.0171800000e-01),
                (2.1060000000e+01, 2.7474000000e-01),
                (7.4950000000e+00, 4.4856400000e-01),
                (2.7970000000e+00, 2.8507400000e-01),
                (5.2150000000e-01, 1.5204000000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (6.6650000000e+03, -1.4600000000e-04),
                (1.0000000000e+03, -1.1540000000e-03),
                (2.2800000000e+02, -5.7250000000e-03),
                (6.4710000000e+01, -2.3312000000e-02),
                (2.1060000000e+01, -6.3955000000e-02),
                (7.4950000000e+00, -1.4998100000e-01),
                (2.7970000000e+00, -1.2726200000e-01),
                (5.2150000000e-01, 5.4452900000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.5960000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (9.4390000000e+00, 3.8109000000e-02),
                (2.0020000000e+00, 2.0948000000e-01),
                (5.4560000000e-01, 5.0855700000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.5170000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(5.5000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Nitrogen cc-pVDZ basis
/// (9s4p1d) -> [3s2p1d]: 15 Cartesian basis functions
fn ccpvdz_n() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (9.0460000000e+03, 7.0000000000e-04),
                (1.3570000000e+03, 5.3890000000e-03),
                (3.0930000000e+02, 2.7406000000e-02),
                (8.7730000000e+01, 1.0320700000e-01),
                (2.8560000000e+01, 2.7872300000e-01),
                (1.0210000000e+01, 4.4854000000e-01),
                (3.8380000000e+00, 2.7823800000e-01),
                (7.4660000000e-01, 1.5440000000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (9.0460000000e+03, -1.5300000000e-04),
                (1.3570000000e+03, -1.2080000000e-03),
                (3.0930000000e+02, -5.9920000000e-03),
                (8.7730000000e+01, -2.4544000000e-02),
                (2.8560000000e+01, -6.7459000000e-02),
                (1.0210000000e+01, -1.5807800000e-01),
                (3.8380000000e+00, -1.2183100000e-01),
                (7.4660000000e-01, 5.4900300000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(2.2480000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.3550000000e+01, 3.9919000000e-02),
                (2.9170000000e+00, 2.1716900000e-01),
                (7.9730000000e-01, 5.1031900000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(2.1850000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(8.1700000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Oxygen cc-pVDZ basis
/// (9s4p1d) -> [3s2p1d]: 15 Cartesian basis functions
fn ccpvdz_o() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.1720000000e+04, 7.1000000000e-04),
                (1.7590000000e+03, 5.4700000000e-03),
                (4.0080000000e+02, 2.7837000000e-02),
                (1.1370000000e+02, 1.0480000000e-01),
                (3.7030000000e+01, 2.8306200000e-01),
                (1.3270000000e+01, 4.4871900000e-01),
                (5.0250000000e+00, 2.7095200000e-01),
                (1.0130000000e+00, 1.5458000000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.1720000000e+04, -1.6000000000e-04),
                (1.7590000000e+03, -1.2630000000e-03),
                (4.0080000000e+02, -6.2670000000e-03),
                (1.1370000000e+02, -2.5716000000e-02),
                (3.7030000000e+01, -7.0924000000e-02),
                (1.3270000000e+01, -1.6541100000e-01),
                (5.0250000000e+00, -1.1695500000e-01),
                (1.0130000000e+00, 5.5736800000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(3.0230000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.7700000000e+01, 4.3018000000e-02),
                (3.8540000000e+00, 2.2891300000e-01),
                (1.0460000000e+00, 5.0872800000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(2.7530000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(1.1850000000e+00, 1.0000000000e+00)],
        ),
    ]
}

/// Fluorine cc-pVDZ basis
/// (9s4p1d) -> [3s2p1d]: 15 Cartesian basis functions
fn ccpvdz_f() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.4710000000e+04, 7.2100000000e-04),
                (2.2070000000e+03, 5.5530000000e-03),
                (5.0280000000e+02, 2.8267000000e-02),
                (1.4260000000e+02, 1.0644400000e-01),
                (4.6470000000e+01, 2.8681400000e-01),
                (1.6700000000e+01, 4.4864100000e-01),
                (6.3560000000e+00, 2.6476100000e-01),
                (1.3160000000e+00, 1.5333000000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.4710000000e+04, -1.6500000000e-04),
                (2.2070000000e+03, -1.3080000000e-03),
                (5.0280000000e+02, -6.4950000000e-03),
                (1.4260000000e+02, -2.6691000000e-02),
                (4.6470000000e+01, -7.3690000000e-02),
                (1.6700000000e+01, -1.7077600000e-01),
                (6.3560000000e+00, -1.1232700000e-01),
                (1.3160000000e+00, 5.6281400000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(3.8970000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.2670000000e+01, 4.4878000000e-02),
                (4.9770000000e+00, 2.3571800000e-01),
                (1.3470000000e+00, 5.0852100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(3.4710000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(1.6400000000e+00, 1.0000000000e+00)],
        ),
    ]
}

/// Neon cc-pVDZ basis
/// (9s4p1d) -> [3s2p1d]: 15 Cartesian basis functions
fn ccpvdz_ne() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.7880000000e+04, 7.3800000000e-04),
                (2.6830000000e+03, 5.6770000000e-03),
                (6.1150000000e+02, 2.8883000000e-02),
                (1.7350000000e+02, 1.0854000000e-01),
                (5.6640000000e+01, 2.9090700000e-01),
                (2.0420000000e+01, 4.4832400000e-01),
                (7.8100000000e+00, 2.5802600000e-01),
                (1.6530000000e+00, 1.5063000000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.7880000000e+04, -1.7200000000e-04),
                (2.6830000000e+03, -1.3570000000e-03),
                (6.1150000000e+02, -6.7370000000e-03),
                (1.7350000000e+02, -2.7663000000e-02),
                (5.6640000000e+01, -7.6208000000e-02),
                (2.0420000000e+01, -1.7522700000e-01),
                (7.8100000000e+00, -1.0703800000e-01),
                (1.6530000000e+00, 5.6705000000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(4.8690000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.8390000000e+01, 4.6087000000e-02),
                (6.2700000000e+00, 2.4018100000e-01),
                (1.6950000000e+00, 5.0874400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(4.3170000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(2.2020000000e+00, 1.0000000000e+00)],
        ),
    ]
}

/// Sodium cc-pVDZ basis
/// (12s8p1d) -> [4s3p1d]: 19 Cartesian basis functions
/// General contraction: 11-primitive s-block produces 3 contracted s functions
fn ccpvdz_na() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (3.1700000000e+04, 4.5887800000e-04),
                (4.7550000000e+03, 3.5507000000e-03),
                (1.0820000000e+03, 1.8261800000e-02),
                (3.0640000000e+02, 7.1665000000e-02),
                (9.9530000000e+01, 2.1234600000e-01),
                (3.5420000000e+01, 4.1620300000e-01),
                (1.3300000000e+01, 3.7302000000e-01),
                (4.3920000000e+00, 6.2505400000e-02),
                (1.6760000000e+00, -6.2453200000e-03),
                (5.8890000000e-01, 2.4337400000e-03),
                (5.6400000000e-02, -4.4238100000e-04),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.1700000000e+04, -1.1216200000e-04),
                (4.7550000000e+03, -8.6851200000e-04),
                (1.0820000000e+03, -4.5133000000e-03),
                (3.0640000000e+02, -1.8143600000e-02),
                (9.9530000000e+01, -5.8079900000e-02),
                (3.5420000000e+01, -1.3765300000e-01),
                (1.3300000000e+01, -1.9390800000e-01),
                (4.3920000000e+00, 8.5800900000e-02),
                (1.6760000000e+00, 6.0441900000e-01),
                (5.8890000000e-01, 4.4171900000e-01),
                (5.6400000000e-02, 1.3054700000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (3.1700000000e+04, 1.7016000000e-05),
                (4.7550000000e+03, 1.3069300000e-04),
                (1.0820000000e+03, 6.8778400000e-04),
                (3.0640000000e+02, 2.7235900000e-03),
                (9.9530000000e+01, 8.9552900000e-03),
                (3.5420000000e+01, 2.0783200000e-02),
                (1.3300000000e+01, 3.1938000000e-02),
                (4.3920000000e+00, -1.9136800000e-02),
                (1.6760000000e+00, -1.0259500000e-01),
                (5.8890000000e-01, -1.9894500000e-01),
                (5.6400000000e-02, 6.5595200000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(2.3070000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.3810000000e+02, 5.7964100000e-03),
                (3.2240000000e+01, 4.1575600000e-02),
                (9.9850000000e+00, 1.6287300000e-01),
                (3.4840000000e+00, 3.5940100000e-01),
                (1.2310000000e+00, 4.4998800000e-01),
                (4.1770000000e-01, 2.2750700000e-01),
                (6.5130000000e-02, 8.0824700000e-03),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.3810000000e+02, -5.8153100000e-04),
                (3.2240000000e+01, -4.0730600000e-03),
                (9.9850000000e+00, -1.6793700000e-02),
                (3.4840000000e+00, -3.5326800000e-02),
                (1.2310000000e+00, -5.2197100000e-02),
                (4.1770000000e-01, -1.6835900000e-02),
                (6.5130000000e-02, 4.3461300000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(2.0530000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(9.7300000000e-02, 1.0000000000e+00)],
        ),
    ]
}

/// Magnesium cc-pVDZ basis
/// (12s8p1d) -> [4s3p1d]: 19 Cartesian basis functions
fn ccpvdz_mg() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (4.7390000000e+04, 3.4602300000e-04),
                (7.1080000000e+03, 2.6807700000e-03),
                (1.6180000000e+03, 1.3836700000e-02),
                (4.5840000000e+02, 5.5176700000e-02),
                (1.4930000000e+02, 1.6966000000e-01),
                (5.3590000000e+01, 3.6470300000e-01),
                (2.0700000000e+01, 4.0685600000e-01),
                (8.3840000000e+00, 1.3508900000e-01),
                (2.5420000000e+00, 4.9088400000e-03),
                (8.7870000000e-01, 2.8646000000e-04),
                (1.0770000000e-01, 2.6459000000e-05),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.7390000000e+04, -8.7783900000e-05),
                (7.1080000000e+03, -6.7472500000e-04),
                (1.6180000000e+03, -3.5560300000e-03),
                (4.5840000000e+02, -1.4215400000e-02),
                (1.4930000000e+02, -4.7674800000e-02),
                (5.3590000000e+01, -1.1489200000e-01),
                (2.0700000000e+01, -2.0067600000e-01),
                (8.3840000000e+00, -3.4122400000e-02),
                (2.5420000000e+00, 5.7045400000e-01),
                (8.7870000000e-01, 5.4230900000e-01),
                (1.0770000000e-01, 2.1812800000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (4.7390000000e+04, 1.6962800000e-05),
                (7.1080000000e+03, 1.2986500000e-04),
                (1.6180000000e+03, 6.8883100000e-04),
                (4.5840000000e+02, 2.7353300000e-03),
                (1.4930000000e+02, 9.3122400000e-03),
                (5.3590000000e+01, 2.2326500000e-02),
                (2.0700000000e+01, 4.1119500000e-02),
                (8.3840000000e+00, 5.4564200000e-03),
                (2.5420000000e+00, -1.3401200000e-01),
                (8.7870000000e-01, -2.5617600000e-01),
                (1.0770000000e-01, 6.0585600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(3.9990000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.7990000000e+02, 5.3816100000e-03),
                (4.2140000000e+01, 3.9241800000e-02),
                (1.3130000000e+01, 1.5744500000e-01),
                (4.6280000000e+00, 3.5853500000e-01),
                (1.6700000000e+00, 4.5722600000e-01),
                (5.8570000000e-01, 2.1591800000e-01),
                (1.3110000000e-01, 6.6494800000e-03),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (1.7990000000e+02, -8.6594800000e-04),
                (4.2140000000e+01, -6.1597800000e-03),
                (1.3130000000e+01, -2.6151900000e-02),
                (4.6280000000e+00, -5.7064700000e-02),
                (1.6700000000e+00, -8.7390600000e-02),
                (5.8570000000e-01, -1.2299000000e-02),
                (1.3110000000e-01, 5.0208500000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(4.1120000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(1.8700000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Aluminum cc-pVDZ basis
/// (12s8p1d) -> [4s3p1d]: 19 Cartesian basis functions
fn ccpvdz_al() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (6.4150000000e+04, 2.9025000000e-04),
                (9.6170000000e+03, 2.2506400000e-03),
                (2.1890000000e+03, 1.1645900000e-02),
                (6.2050000000e+02, 4.6737700000e-02),
                (2.0270000000e+02, 1.4629900000e-01),
                (7.3150000000e+01, 3.3028300000e-01),
                (2.8550000000e+01, 4.1586100000e-01),
                (1.1770000000e+01, 1.8925300000e-01),
                (3.3000000000e+00, 1.1588900000e-02),
                (1.1730000000e+00, -1.2838500000e-03),
                (1.7520000000e-01, 4.2588300000e-04),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (6.4150000000e+04, -7.5804800000e-05),
                (9.6170000000e+03, -5.8179100000e-04),
                (2.1890000000e+03, -3.0811300000e-03),
                (6.2050000000e+02, -1.2311200000e-02),
                (2.0270000000e+02, -4.1978100000e-02),
                (7.3150000000e+01, -1.0337100000e-01),
                (2.8550000000e+01, -1.9630800000e-01),
                (1.1770000000e+01, -8.3000200000e-02),
                (3.3000000000e+00, 5.4104000000e-01),
                (1.1730000000e+00, 5.7879600000e-01),
                (1.7520000000e-01, 2.8814700000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (6.4150000000e+04, 1.7507800000e-05),
                (9.6170000000e+03, 1.3420800000e-04),
                (2.1890000000e+03, 7.1244200000e-04),
                (6.2050000000e+02, 2.8433000000e-03),
                (2.0270000000e+02, 9.7684200000e-03),
                (7.3150000000e+01, 2.4185000000e-02),
                (2.8550000000e+01, 4.7499300000e-02),
                (1.1770000000e+01, 2.0362100000e-02),
                (3.3000000000e+00, -1.5878800000e-01),
                (1.1730000000e+00, -3.1169400000e-01),
                (1.7520000000e-01, 6.2014700000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(6.4730000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.5880000000e+02, 4.0684700000e-03),
                (6.0890000000e+01, 3.0681500000e-02),
                (1.9140000000e+01, 1.2914900000e-01),
                (6.8810000000e+00, 3.2083100000e-01),
                (2.5740000000e+00, 4.5381500000e-01),
                (9.5720000000e-01, 2.7506600000e-01),
                (2.0990000000e-01, 1.9080700000e-02),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (2.5880000000e+02, -7.4805300000e-04),
                (6.0890000000e+01, -5.4579600000e-03),
                (1.9140000000e+01, -2.4537100000e-02),
                (6.8810000000e+00, -5.8213800000e-02),
                (2.5740000000e+00, -9.8375600000e-02),
                (9.5720000000e-01, -2.6006400000e-02),
                (2.0990000000e-01, 4.6402000000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(5.9860000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(1.8900000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Silicon cc-pVDZ basis
/// (12s8p1d) -> [4s3p1d]: 19 Cartesian basis functions
fn ccpvdz_si() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (7.8860000000e+04, 2.7044300000e-04),
                (1.1820000000e+04, 2.0971700000e-03),
                (2.6920000000e+03, 1.0850600000e-02),
                (7.6340000000e+02, 4.3675400000e-02),
                (2.4960000000e+02, 1.3765300000e-01),
                (9.0280000000e+01, 3.1664400000e-01),
                (3.5290000000e+01, 4.1858100000e-01),
                (1.4510000000e+01, 2.1021200000e-01),
                (4.0530000000e+00, 1.4495200000e-02),
                (1.4820000000e+00, -2.0359000000e-03),
                (2.5170000000e-01, 6.2418600000e-04),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (7.8860000000e+04, -7.2317700000e-05),
                (1.1820000000e+04, -5.5511600000e-04),
                (2.6920000000e+03, -2.9380500000e-03),
                (7.6340000000e+02, -1.1768700000e-02),
                (2.4960000000e+02, -4.0290700000e-02),
                (9.0280000000e+01, -1.0060900000e-01),
                (3.5290000000e+01, -1.9652800000e-01),
                (1.4510000000e+01, -1.0238200000e-01),
                (4.0530000000e+00, 5.2719000000e-01),
                (1.4820000000e+00, 5.9325100000e-01),
                (2.5170000000e-01, 3.3265200000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (7.8860000000e+04, 1.8511300000e-05),
                (1.1820000000e+04, 1.4223600000e-04),
                (2.6920000000e+03, 7.5218500000e-04),
                (7.6340000000e+02, 3.0227900000e-03),
                (2.4960000000e+02, 1.0367700000e-02),
                (9.0280000000e+01, 2.6256300000e-02),
                (3.5290000000e+01, 5.2398900000e-02),
                (1.4510000000e+01, 2.9095900000e-02),
                (4.0530000000e+00, -1.7800300000e-01),
                (1.4820000000e+00, -3.4687400000e-01),
                (2.5170000000e-01, 6.2302000000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(9.2430000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.1590000000e+02, 3.9265600000e-03),
                (7.4420000000e+01, 2.9881100000e-02),
                (2.3480000000e+01, 1.2721200000e-01),
                (8.4880000000e+00, 3.2094300000e-01),
                (3.2170000000e+00, 4.5542900000e-01),
                (1.2290000000e+00, 2.6856300000e-01),
                (2.9640000000e-01, 1.8833600000e-02),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.1590000000e+02, -8.5830200000e-04),
                (7.4420000000e+01, -6.3032800000e-03),
                (2.3480000000e+01, -2.8825500000e-02),
                (8.4880000000e+00, -6.9456000000e-02),
                (3.2170000000e+00, -1.1949300000e-01),
                (1.2290000000e+00, -1.9958100000e-02),
                (2.9640000000e-01, 5.1026800000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(8.7680000000e-02, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(2.7500000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Phosphorus cc-pVDZ basis
/// (12s8p1d) -> [4s3p1d]: 19 Cartesian basis functions
fn ccpvdz_p() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (9.4840000000e+04, 2.5550900000e-04),
                (1.4220000000e+04, 1.9819300000e-03),
                (3.2360000000e+03, 1.0276000000e-02),
                (9.1710000000e+02, 4.1482300000e-02),
                (2.9950000000e+02, 1.3198400000e-01),
                (1.0810000000e+02, 3.0866200000e-01),
                (4.2180000000e+01, 4.2064700000e-01),
                (1.7280000000e+01, 2.2287800000e-01),
                (4.8580000000e+00, 1.6403500000e-02),
                (1.8180000000e+00, -2.5425500000e-03),
                (3.3720000000e-01, 7.4805000000e-04),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (9.4840000000e+04, -6.9693900000e-05),
                (1.4220000000e+04, -5.3526600000e-04),
                (3.2360000000e+03, -2.8370900000e-03),
                (9.1710000000e+02, -1.1398300000e-02),
                (2.9950000000e+02, -3.9292900000e-02),
                (1.0810000000e+02, -9.9636400000e-02),
                (4.2180000000e+01, -1.9798300000e-01),
                (1.7280000000e+01, -1.1486000000e-01),
                (4.8580000000e+00, 5.1859500000e-01),
                (1.8180000000e+00, 6.0184700000e-01),
                (3.3720000000e-01, 3.6861200000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (9.4840000000e+04, 1.9119900000e-05),
                (1.4220000000e+04, 1.4722300000e-04),
                (3.2360000000e+03, 7.7791200000e-04),
                (9.1710000000e+02, 3.1454600000e-03),
                (2.9950000000e+02, 1.0820000000e-02),
                (1.0810000000e+02, 2.7995700000e-02),
                (4.2180000000e+01, 5.6397800000e-02),
                (1.7280000000e+01, 3.5819000000e-02),
                (4.8580000000e+00, -1.9338700000e-01),
                (1.8180000000e+00, -3.7209700000e-01),
                (3.3720000000e-01, 6.2424600000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.2320000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.7050000000e+02, 3.9500500000e-03),
                (8.7330000000e+01, 3.0249200000e-02),
                (2.7590000000e+01, 1.2955400000e-01),
                (1.0000000000e+01, 3.2759400000e-01),
                (3.8250000000e+00, 4.5699200000e-01),
                (1.4940000000e+00, 2.5308600000e-01),
                (3.9210000000e-01, 1.6879800000e-02),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.7050000000e+02, -9.5983200000e-04),
                (8.7330000000e+01, -7.1117700000e-03),
                (2.7590000000e+01, -3.2712200000e-02),
                (1.0000000000e+01, -7.9578400000e-02),
                (3.8250000000e+00, -1.3501600000e-01),
                (1.4940000000e+00, -9.1058500000e-03),
                (3.9210000000e-01, 5.3780200000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.1860000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(3.7300000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Sulfur cc-pVDZ basis
/// (12s8p1d) -> [4s3p1d]: 19 Cartesian basis functions
fn ccpvdz_s() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.1080000000e+05, 2.4763500000e-04),
                (1.6610000000e+04, 1.9202600000e-03),
                (3.7810000000e+03, 9.9619200000e-03),
                (1.0710000000e+03, 4.0297500000e-02),
                (3.4980000000e+02, 1.2860400000e-01),
                (1.2630000000e+02, 3.0348000000e-01),
                (4.9260000000e+01, 4.2143200000e-01),
                (2.0160000000e+01, 2.3078100000e-01),
                (5.7200000000e+00, 1.7897100000e-02),
                (2.1820000000e+00, -2.9751600000e-03),
                (4.3270000000e-01, 8.4952200000e-04),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.1080000000e+05, -6.8703900000e-05),
                (1.6610000000e+04, -5.2768100000e-04),
                (3.7810000000e+03, -2.7967100000e-03),
                (1.0710000000e+03, -1.1265100000e-02),
                (3.4980000000e+02, -3.8883400000e-02),
                (1.2630000000e+02, -9.9502500000e-02),
                (4.9260000000e+01, -1.9974000000e-01),
                (2.0160000000e+01, -1.2336000000e-01),
                (5.7200000000e+00, 5.1319400000e-01),
                (2.1820000000e+00, 6.0712000000e-01),
                (4.3270000000e-01, 3.9675300000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.1080000000e+05, 1.9907700000e-05),
                (1.6610000000e+04, 1.5348300000e-04),
                (3.7810000000e+03, 8.0950300000e-04),
                (1.0710000000e+03, 3.2897400000e-03),
                (3.4980000000e+02, 1.1296700000e-02),
                (1.2630000000e+02, 2.9638500000e-02),
                (4.9260000000e+01, 5.9985100000e-02),
                (2.0160000000e+01, 4.1324800000e-02),
                (5.7200000000e+00, -2.0747400000e-01),
                (2.1820000000e+00, -3.9288900000e-01),
                (4.3270000000e-01, 6.3284000000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.5700000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.9970000000e+02, 4.4754100000e-03),
                (9.4190000000e+01, 3.4170800000e-02),
                (2.9750000000e+01, 1.4425000000e-01),
                (1.0770000000e+01, 3.5392800000e-01),
                (4.1190000000e+00, 4.5908500000e-01),
                (1.6250000000e+00, 2.0638300000e-01),
                (4.7260000000e-01, 1.0214100000e-02),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (3.9970000000e+02, -1.1625100000e-03),
                (9.4190000000e+01, -8.6566400000e-03),
                (2.9750000000e+01, -3.9088600000e-02),
                (1.0770000000e+01, -9.3462500000e-02),
                (4.1190000000e+00, -1.4799400000e-01),
                (1.6250000000e+00, 3.0190400000e-02),
                (4.7260000000e-01, 5.6157300000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.4070000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(4.7900000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Chlorine cc-pVDZ basis
/// (12s8p1d) -> [4s3p1d]: 19 Cartesian basis functions
fn ccpvdz_cl() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.2790000000e+05, 2.4115300000e-04),
                (1.9170000000e+04, 1.8709500000e-03),
                (4.3630000000e+03, 9.7082700000e-03),
                (1.2360000000e+03, 3.9315300000e-02),
                (4.0360000000e+02, 1.2593200000e-01),
                (1.4570000000e+02, 2.9934100000e-01),
                (5.6810000000e+01, 4.2188600000e-01),
                (2.3230000000e+01, 2.3720100000e-01),
                (6.6440000000e+00, 1.9153100000e-02),
                (2.5750000000e+00, -3.3479200000e-03),
                (5.3710000000e-01, 9.2988300000e-04),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.2790000000e+05, -6.7892200000e-05),
                (1.9170000000e+04, -5.2183600000e-04),
                (4.3630000000e+03, -2.7651300000e-03),
                (1.2360000000e+03, -1.1153700000e-02),
                (4.0360000000e+02, -3.8591900000e-02),
                (1.4570000000e+02, -9.9484800000e-02),
                (5.6810000000e+01, -2.0139200000e-01),
                (2.3230000000e+01, -1.3031300000e-01),
                (6.6440000000e+00, 5.0944300000e-01),
                (2.5750000000e+00, 6.1072500000e-01),
                (5.3710000000e-01, 4.2154900000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.2790000000e+05, 2.0498600000e-05),
                (1.9170000000e+04, 1.5829800000e-04),
                (4.3630000000e+03, 8.3363900000e-04),
                (1.2360000000e+03, 3.3988000000e-03),
                (4.0360000000e+02, 1.1673800000e-02),
                (1.4570000000e+02, 3.0962200000e-02),
                (5.6810000000e+01, 6.2953300000e-02),
                (2.3230000000e+01, 4.6025700000e-02),
                (6.6440000000e+00, -2.1931200000e-01),
                (2.5750000000e+00, -4.0877300000e-01),
                (5.3710000000e-01, 6.3846500000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(1.9380000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.1760000000e+02, 5.2598200000e-03),
                (9.8330000000e+01, 3.9833200000e-02),
                (3.1040000000e+01, 1.6465500000e-01),
                (1.1190000000e+01, 3.8732200000e-01),
                (4.2490000000e+00, 4.5707200000e-01),
                (1.6240000000e+00, 1.5163600000e-01),
                (5.3220000000e-01, 1.8161500000e-03),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.1760000000e+02, -1.4357000000e-03),
                (9.8330000000e+01, -1.0779600000e-02),
                (3.1040000000e+01, -4.7007500000e-02),
                (1.1190000000e+01, -1.1103000000e-01),
                (4.2490000000e+00, -1.5327500000e-01),
                (1.6240000000e+00, 8.9460900000e-02),
                (5.3220000000e-01, 5.7944400000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.6200000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(6.0000000000e-01, 1.0000000000e+00)],
        ),
    ]
}

/// Argon cc-pVDZ basis
/// (12s8p1d) -> [4s3p1d]: 19 Cartesian basis functions
fn ccpvdz_ar() -> Vec<ShellData> {
    vec![
        (
            AngularMomentum::S,
            vec![
                (1.4570000000e+05, 2.3670000000e-04),
                (2.1840000000e+04, 1.8352300000e-03),
                (4.9720000000e+03, 9.5286000000e-03),
                (1.4080000000e+03, 3.8628300000e-02),
                (4.5970000000e+02, 1.2408100000e-01),
                (1.6590000000e+02, 2.9647100000e-01),
                (6.4690000000e+01, 4.2206800000e-01),
                (2.6440000000e+01, 2.4171100000e-01),
                (7.6280000000e+00, 2.0050900000e-02),
                (2.9960000000e+00, -3.6100000000e-03),
                (6.5040000000e-01, 9.7560700000e-04),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.4570000000e+05, -6.7491000000e-05),
                (2.1840000000e+04, -5.1852200000e-04),
                (4.9720000000e+03, -2.7482500000e-03),
                (1.4080000000e+03, -1.1100700000e-02),
                (4.5970000000e+02, -3.8482000000e-02),
                (1.6590000000e+02, -9.9759900000e-02),
                (6.4690000000e+01, -2.0308800000e-01),
                (2.6440000000e+01, -1.3560800000e-01),
                (7.6280000000e+00, 5.0719500000e-01),
                (2.9960000000e+00, 6.1289800000e-01),
                (6.5040000000e-01, 4.4296800000e-02),
            ],
        ),
        (
            AngularMomentum::S,
            vec![
                (1.4570000000e+05, 2.1045700000e-05),
                (2.1840000000e+04, 1.6256500000e-04),
                (4.9720000000e+03, 8.5546300000e-04),
                (1.4080000000e+03, 3.4974500000e-03),
                (4.5970000000e+02, 1.2015600000e-02),
                (1.6590000000e+02, 3.2136800000e-02),
                (6.4690000000e+01, 6.5527900000e-02),
                (2.6440000000e+01, 4.9937000000e-02),
                (7.6280000000e+00, -2.2976900000e-01),
                (2.9960000000e+00, -4.2100600000e-01),
                (6.5040000000e-01, 6.4233100000e-01),
            ],
        ),
        (
            AngularMomentum::S,
            vec![(2.3370000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.5370000000e+02, 5.7055500000e-03),
                (1.0680000000e+02, 4.3046000000e-02),
                (3.3730000000e+01, 1.7659100000e-01),
                (1.2130000000e+01, 4.0686300000e-01),
                (4.5940000000e+00, 4.5254900000e-01),
                (1.6780000000e+00, 1.2280100000e-01),
                (5.9090000000e-01, -4.4599600000e-03),
            ],
        ),
        (
            AngularMomentum::P,
            vec![
                (4.5370000000e+02, -1.6065500000e-03),
                (1.0680000000e+02, -1.2171400000e-02),
                (3.3730000000e+01, -5.2078900000e-02),
                (1.2130000000e+01, -1.2373700000e-01),
                (4.5940000000e+00, -1.5161900000e-01),
                (1.6780000000e+00, 1.4242500000e-01),
                (5.9090000000e-01, 5.8450100000e-01),
            ],
        ),
        (
            AngularMomentum::P,
            vec![(1.8520000000e-01, 1.0000000000e+00)],
        ),
        (
            AngularMomentum::D,
            vec![(7.3800000000e-01, 1.0000000000e+00)],
        ),
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
        // Z = 19 (Potassium) not supported
        let result = get_element_basis(19, "sto-3g");
        assert!(matches!(result, Err(BasisError::UnsupportedElement(19))));

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
        assert!(is_supported_basis("cc-pvdz"));
        assert!(is_supported_basis("cc-pVDZ"));
        assert!(is_supported_basis("ccpvdz"));
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

    // -------------------------------------------------------------------------
    // Third-row element tests (Na-Ar, Z=11-18)
    // Shell counts verified against PySCF 2.11.0 mol.nbas
    // -------------------------------------------------------------------------

    #[test]
    fn test_sto3g_na() {
        let basis = get_element_basis(11, "sto-3g").unwrap();
        assert_eq!(basis.len(), 5); // 3 S + 2 P
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 9); // 3*1 + 2*3 = 9
    }

    #[test]
    fn test_sto3g_si() {
        let basis = get_element_basis(14, "sto-3g").unwrap();
        assert_eq!(basis.len(), 5);
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 9);
    }

    #[test]
    fn test_sto3g_cl() {
        let basis = get_element_basis(17, "sto-3g").unwrap();
        assert_eq!(basis.len(), 5);
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 9);
    }

    #[test]
    fn test_sto3g_ar() {
        let basis = get_element_basis(18, "sto-3g").unwrap();
        assert_eq!(basis.len(), 5);
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 9);
    }

    #[test]
    fn test_321g_na() {
        let basis = get_element_basis(11, "3-21g").unwrap();
        assert_eq!(basis.len(), 7); // 4 S + 3 P
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 13); // 4*1 + 3*3 = 13
    }

    #[test]
    fn test_321g_si() {
        let basis = get_element_basis(14, "3-21g").unwrap();
        assert_eq!(basis.len(), 7);
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 13);
    }

    #[test]
    fn test_321g_cl() {
        let basis = get_element_basis(17, "3-21g").unwrap();
        assert_eq!(basis.len(), 7);
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 13);
    }

    #[test]
    fn test_631g_na() {
        let basis = get_element_basis(11, "6-31g").unwrap();
        assert_eq!(basis.len(), 7); // 4 S + 3 P
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 13);
    }

    #[test]
    fn test_631g_p_element() {
        let basis = get_element_basis(15, "6-31g").unwrap();
        assert_eq!(basis.len(), 7);
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 13);
    }

    #[test]
    fn test_631gs_si() {
        let basis = get_element_basis(14, "6-31g*").unwrap();
        assert_eq!(basis.len(), 8); // 4 S + 3 P + 1 D
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 19); // 4*1 + 3*3 + 1*6 = 19
    }

    #[test]
    fn test_631gs_cl() {
        let basis = get_element_basis(17, "6-31g*").unwrap();
        assert_eq!(basis.len(), 8);
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 19);
    }

    #[test]
    fn test_631pgs_al() {
        let basis = get_element_basis(13, "6-31+g*").unwrap();
        assert_eq!(basis.len(), 10); // 5 S + 4 P + 1 D
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        // PySCF 2.11.0 confirms: nao_cart=23 = 5*1 + 4*3 + 1*6
        assert_eq!(nao, 23);
    }

    #[test]
    fn test_631pgs_ar() {
        let basis = get_element_basis(18, "6-31+g*").unwrap();
        assert_eq!(basis.len(), 10);
        let nao: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
        assert_eq!(nao, 23);
    }

    // -------------------------------------------------------------------------
    // Exponent spot-check tests (vs PySCF 2.11.0)
    // -------------------------------------------------------------------------

    #[test]
    fn test_na_sto3g_1s_exponent() {
        let basis = get_element_basis(11, "sto-3g").unwrap();
        // First shell (1s), first exponent
        assert!((basis[0].1[0].0 - 250.77243).abs() < 1e-4);
    }

    #[test]
    fn test_cl_sto3g_1s_exponent() {
        let basis = get_element_basis(17, "sto-3g").unwrap();
        assert!((basis[0].1[0].0 - 601.3456136).abs() < 1e-4);
    }

    #[test]
    fn test_si_631gs_d_exponent() {
        let basis = get_element_basis(14, "6-31g*").unwrap();
        let d_shell = basis
            .iter()
            .find(|(am, _)| *am == AngularMomentum::D)
            .unwrap();
        assert!((d_shell.1[0].0 - 0.45).abs() < 1e-6);
    }

    #[test]
    fn test_s_631pgs_diffuse_exponent() {
        let basis = get_element_basis(16, "6-31+g*").unwrap();
        // Find diffuse S shell (single primitive with small exponent)
        let diffuse_s = basis
            .iter()
            .find(|(am, prims)| *am == AngularMomentum::S && prims.len() == 1 && prims[0].0 < 0.1);
        assert!(diffuse_s.is_some(), "Should have diffuse S shell");
        assert!((diffuse_s.unwrap().1[0].0 - 0.0405).abs() < 1e-3);
    }

    // -------------------------------------------------------------------------
    // Shell type tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_na_sto3g_shell_types() {
        let basis = get_element_basis(11, "sto-3g").unwrap();
        // First 3 shells are S, last 2 are P
        assert_eq!(basis[0].0, AngularMomentum::S);
        assert_eq!(basis[1].0, AngularMomentum::S);
        assert_eq!(basis[2].0, AngularMomentum::S);
        assert_eq!(basis[3].0, AngularMomentum::P);
        assert_eq!(basis[4].0, AngularMomentum::P);
    }

    #[test]
    fn test_si_631gs_has_d() {
        let basis = get_element_basis(14, "6-31g*").unwrap();
        let last = basis.last().unwrap();
        assert_eq!(last.0, AngularMomentum::D);
    }

    #[test]
    fn test_cl_631pgs_has_diffuse() {
        let basis = get_element_basis(17, "6-31+g*").unwrap();
        let has_diffuse_s = basis
            .iter()
            .any(|(am, prims)| *am == AngularMomentum::S && prims.len() == 1 && prims[0].0 < 0.1);
        assert!(has_diffuse_s, "6-31+G* should have diffuse S on Cl");
    }

    // -------------------------------------------------------------------------
    // All 6 basis sets work for all Na-Ar
    // -------------------------------------------------------------------------

    #[test]
    fn test_all_basis_sets_na_ar() {
        let basis_names = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"];
        let expected_shells = [5, 7, 7, 8, 10, 8];
        for z in 11..=18 {
            for (basis_name, expected) in basis_names.iter().zip(expected_shells.iter()) {
                let result = get_element_basis(z, basis_name);
                assert!(
                    result.is_ok(),
                    "Z={} basis={} should succeed",
                    z,
                    basis_name
                );
                assert_eq!(
                    result.unwrap().len(),
                    *expected,
                    "Z={} basis={} should have {} shells",
                    z,
                    basis_name,
                    expected,
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // cc-pVDZ specific tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_ccpvdz_loads_all_elements() {
        // Verify cc-pVDZ loads for all elements H-Ar (Z=1-18)
        for z in 1..=18u8 {
            let result = get_element_basis(z, "cc-pvdz");
            assert!(result.is_ok(), "cc-pVDZ should load for Z={}", z);
        }
    }

    #[test]
    fn test_ccpvdz_alternate_names() {
        // All alternate name forms should work
        let h_basis1 = get_element_basis(1, "cc-pvdz").unwrap();
        let h_basis2 = get_element_basis(1, "ccpvdz").unwrap();
        let h_basis3 = get_element_basis(1, "CC-PVDZ").unwrap();

        assert_eq!(h_basis1.len(), h_basis2.len());
        assert_eq!(h_basis1.len(), h_basis3.len());
    }

    #[test]
    fn test_ccpvdz_h_shell_structure() {
        // H cc-pVDZ: (4s1p) -> [2s1p]
        // Shells: 2 S shells + 1 P shell = 3 shells
        // Cartesian BFs: 2(s) + 3(p) = 5
        let basis = get_element_basis(1, "cc-pvdz").unwrap();
        assert_eq!(basis.len(), 3); // 3 shells

        // Shell types
        assert_eq!(basis[0].0, AngularMomentum::S);
        assert_eq!(basis[1].0, AngularMomentum::S);
        assert_eq!(basis[2].0, AngularMomentum::P);

        // Shell sizes
        assert_eq!(basis[0].1.len(), 3); // 3-primitive contracted
        assert_eq!(basis[1].1.len(), 1); // uncontracted
        assert_eq!(basis[2].1.len(), 1); // uncontracted p
    }

    #[test]
    fn test_ccpvdz_c_shell_structure() {
        // C cc-pVDZ: (9s4p1d) -> [3s2p1d]
        // Shells: 3 S + 2 P + 1 D = 6 shells
        // Cartesian BFs: 3(s) + 6(p) + 6(d) = 15
        let basis = get_element_basis(6, "cc-pvdz").unwrap();
        assert_eq!(basis.len(), 6);

        // Two 8-primitive contracted S shells (general contraction unrolled)
        assert_eq!(basis[0].0, AngularMomentum::S);
        assert_eq!(basis[0].1.len(), 8);
        assert_eq!(basis[1].0, AngularMomentum::S);
        assert_eq!(basis[1].1.len(), 8);

        // One uncontracted S shell
        assert_eq!(basis[2].0, AngularMomentum::S);
        assert_eq!(basis[2].1.len(), 1);

        // One 3-primitive contracted P shell
        assert_eq!(basis[3].0, AngularMomentum::P);
        assert_eq!(basis[3].1.len(), 3);

        // One uncontracted P shell
        assert_eq!(basis[4].0, AngularMomentum::P);
        assert_eq!(basis[4].1.len(), 1);

        // One d polarization shell
        assert_eq!(basis[5].0, AngularMomentum::D);
        assert_eq!(basis[5].1.len(), 1);
    }

    #[test]
    fn test_ccpvdz_na_shell_structure() {
        // Na cc-pVDZ: (12s8p1d) -> [4s3p1d]
        // Shells: 4 S + 3 P + 1 D = 8 shells
        // Cartesian BFs: 4(s) + 9(p) + 6(d) = 19
        let basis = get_element_basis(11, "cc-pvdz").unwrap();
        assert_eq!(basis.len(), 8);

        // Three 11-primitive contracted S shells
        assert_eq!(basis[0].0, AngularMomentum::S);
        assert_eq!(basis[0].1.len(), 11);
        assert_eq!(basis[1].0, AngularMomentum::S);
        assert_eq!(basis[1].1.len(), 11);
        assert_eq!(basis[2].0, AngularMomentum::S);
        assert_eq!(basis[2].1.len(), 11);

        // One uncontracted S shell
        assert_eq!(basis[3].0, AngularMomentum::S);
        assert_eq!(basis[3].1.len(), 1);

        // Two 7-primitive contracted P shells
        assert_eq!(basis[4].0, AngularMomentum::P);
        assert_eq!(basis[4].1.len(), 7);
        assert_eq!(basis[5].0, AngularMomentum::P);
        assert_eq!(basis[5].1.len(), 7);

        // One uncontracted P shell
        assert_eq!(basis[6].0, AngularMomentum::P);
        assert_eq!(basis[6].1.len(), 1);

        // One d polarization shell
        assert_eq!(basis[7].0, AngularMomentum::D);
        assert_eq!(basis[7].1.len(), 1);
    }

    #[test]
    fn test_ccpvdz_basis_function_counts() {
        // Verify Cartesian basis function counts match PySCF with cart=True
        // Expected counts from PySCF 2.11.0:
        // H, He: 5 BFs each (2s + 1p = 2 + 3 = 5)
        // Li-Ne: 15 BFs each (3s + 2p + 1d = 3 + 6 + 6 = 15)
        // Na-Ar: 19 BFs each (4s + 3p + 1d = 4 + 9 + 6 = 19)

        // Expected Cartesian BF count per element
        let expected_counts: [(u8, usize); 18] = [
            (1, 5),
            (2, 5),
            (3, 15),
            (4, 15),
            (5, 15),
            (6, 15),
            (7, 15),
            (8, 15),
            (9, 15),
            (10, 15),
            (11, 19),
            (12, 19),
            (13, 19),
            (14, 19),
            (15, 19),
            (16, 19),
            (17, 19),
            (18, 19),
        ];

        for (z, expected_nao) in expected_counts.iter() {
            let basis = get_element_basis(*z, "cc-pvdz").unwrap();
            let total_bf: usize = basis.iter().map(|(am, _)| am.n_cartesian()).sum();
            assert_eq!(
                total_bf, *expected_nao,
                "Z={} cc-pVDZ should have {} Cartesian BFs, got {}",
                z, expected_nao, total_bf
            );
        }
    }

    #[test]
    fn test_ccpvdz_has_d_polarization() {
        // All elements Li-Ar should have d polarization in cc-pVDZ
        for z in 3..=18u8 {
            let basis = get_element_basis(z, "cc-pvdz").unwrap();
            let has_d = basis.iter().any(|(am, _)| *am == AngularMomentum::D);
            assert!(has_d, "Z={} cc-pVDZ should have D shell", z);
        }
    }

    #[test]
    fn test_ccpvdz_h_has_p_polarization() {
        // Even H and He have p polarization in cc-pVDZ
        for z in 1..=2u8 {
            let basis = get_element_basis(z, "cc-pvdz").unwrap();
            let has_p = basis.iter().any(|(am, _)| *am == AngularMomentum::P);
            assert!(has_p, "Z={} cc-pVDZ should have P shell", z);
        }
    }

    #[test]
    fn test_ccpvdz_general_contraction_exponents_match() {
        // For elements with general contractions (Li-Ne), the first two S shells
        // should share the same exponents (both come from the same general contraction)
        for z in 3..=10u8 {
            let basis = get_element_basis(z, "cc-pvdz").unwrap();
            let s_shells: Vec<_> = basis
                .iter()
                .filter(|(am, _)| *am == AngularMomentum::S)
                .collect();
            assert!(s_shells.len() >= 2, "Z={} should have >= 2 S shells", z);

            // First two S shells should have same number of primitives
            assert_eq!(
                s_shells[0].1.len(),
                s_shells[1].1.len(),
                "Z={} first two S shells should have same number of primitives",
                z
            );

            // Exponents should match
            for (i, (p0, p1)) in s_shells[0].1.iter().zip(s_shells[1].1.iter()).enumerate() {
                assert!(
                    (p0.0 - p1.0).abs() < 1e-8,
                    "Z={} S shell exponent {} mismatch: {} vs {}",
                    z,
                    i,
                    p0.0,
                    p1.0
                );
            }
        }
    }

    #[test]
    fn test_ccpvdz_in_supported_list() {
        let supported = supported_basis_sets();
        assert!(
            supported.contains(&"cc-pvdz"),
            "cc-pvdz should be in supported_basis_sets()"
        );
    }

    #[test]
    fn test_ccpvdz_is_supported() {
        assert!(is_supported_basis("cc-pvdz"));
        assert!(is_supported_basis("ccpvdz"));
        assert!(is_supported_basis("CC-PVDZ"));
        assert!(is_supported_basis("cc_pvdz"));
    }

    #[test]
    fn test_ccpvdz_c_d_exponent() {
        // Carbon cc-pVDZ d-polarization exponent should be 0.55
        let basis = get_element_basis(6, "cc-pvdz").unwrap();
        let d_shell = basis
            .iter()
            .find(|(am, _)| *am == AngularMomentum::D)
            .unwrap();
        assert!((d_shell.1[0].0 - 0.55).abs() < 1e-6);
    }

    #[test]
    fn test_ccpvdz_unsupported_element() {
        // Z=19 (K) is not supported
        let result = get_element_basis(19, "cc-pvdz");
        assert!(result.is_err());
    }
}
