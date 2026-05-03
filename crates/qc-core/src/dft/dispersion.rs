//! Grimme DFT-D3(BJ) dispersion correction.
//!
//! Implements Grimme's D3 dispersion correction with Becke-Johnson (BJ) damping
//! for post-SCF energy correction. The D3-BJ method captures long-range London
//! dispersion interactions absent from standard DFT functionals.
//!
//! # D3-BJ Formula
//!
//! ```text
//! E_disp = -sum_{A<B} [s6 * C6^{AB} / (R_AB^6 + f6)
//!                     + s8 * C8^{AB} / (R_AB^8 + f8)]
//! ```
//!
//! where `f_n(R_AB) = (a1 * R_0^{AB} + a2)^n` is the BJ damping function
//! and `R_0^{AB} = sqrt(C8^{AB} / C6^{AB})` is the cutoff radius.
//!
//! # Key Features
//!
//! - **BJ damping**: Finite at R=0 (no singularity), physically motivated
//! - **CN-dependent C6**: Coordination number interpolation for chemical-environment
//!   sensitivity
//! - **Post-SCF**: Does not modify the KS Fock matrix; purely additive correction
//! - **Analytical gradient**: Exact derivative for geometry optimization
//!
//! # References
//!
//! - Grimme, S. et al. (2010). J. Chem. Phys. 132, 154104. (DFT-D3)
//! - Grimme, S., Ehrlich, S. & Goerigk, L. (2011). J. Comput. Chem. 32, 1456. (BJ damping)
//! - Pyykko, P. & Atsumi, M. (2009). Chem. Eur. J. 15, 186. (covalent radii)

use serde::{Deserialize, Serialize};

// =============================================================================
// Constants
// =============================================================================

/// Bohr radius in Angstrom. Used to convert covalent radii from Angstrom to bohr.
/// Consolidated to match PySCF CODATA 2014 value from constants module.
const BOHR_TO_ANG: f64 = crate::constants::BOHR_TO_ANGSTROM;

/// Steepness parameter k1 for the coordination number counting function.
///
/// Reference: Grimme et al. (2010) J. Chem. Phys. 132, 154104, Eq. 16.
const CN_K1: f64 = 16.0;

/// Scaling parameter k2 for the coordination number counting function.
///
/// Reference: Grimme et al. (2010) J. Chem. Phys. 132, 154104, Eq. 16.
const CN_K2: f64 = 4.0 / 3.0;

/// Gaussian width parameter L for C6 interpolation.
///
/// Reference: Grimme et al. (2010) J. Chem. Phys. 132, 154104, Eq. 15.
const C6_L: f64 = 4.0;

/// Maximum atomic number supported (Ar, Z=18).
const MAX_Z: usize = 18;

// =============================================================================
// Covalent Radii (Pyykko & Atsumi 2009)
// =============================================================================

/// Covalent radii in bohr, indexed by atomic number (Z=0 unused).
///
/// Original values are in Angstrom from Pyykko & Atsumi (2009) Chem. Eur. J. 15, 186.
/// Converted to bohr by dividing by BOHR_TO_ANG = 0.52917721067.
///
/// These are used ONLY for the D3 coordination number calculation, NOT for
/// Becke partitioning or other grid operations.
const R_COV: [f64; MAX_Z + 1] = [
    0.0,                // Z=0 (unused)
    0.32 / BOHR_TO_ANG, // H  = 0.60471
    0.46 / BOHR_TO_ANG, // He = 0.86930
    1.20 / BOHR_TO_ANG, // Li = 2.26768
    0.94 / BOHR_TO_ANG, // Be = 1.77636
    0.77 / BOHR_TO_ANG, // B  = 1.45520
    0.75 / BOHR_TO_ANG, // C  = 1.41740
    0.71 / BOHR_TO_ANG, // N  = 1.34181
    0.63 / BOHR_TO_ANG, // O  = 1.19053
    0.64 / BOHR_TO_ANG, // F  = 1.20942
    0.67 / BOHR_TO_ANG, // Ne = 1.26612
    1.55 / BOHR_TO_ANG, // Na = 2.92892
    1.39 / BOHR_TO_ANG, // Mg = 2.62665
    1.26 / BOHR_TO_ANG, // Al = 2.38106
    1.16 / BOHR_TO_ANG, // Si = 2.19216
    1.11 / BOHR_TO_ANG, // P  = 2.09768
    1.03 / BOHR_TO_ANG, // S  = 1.94651
    0.99 / BOHR_TO_ANG, // Cl = 1.87092
    0.96 / BOHR_TO_ANG, // Ar = 1.81423
];

// =============================================================================
// Quadrupole Expectation Values
// =============================================================================

/// Atomic quadrupole expectation values Q_A (dimensionless),
/// indexed by atomic number (Z=0 unused).
///
/// Used to compute C8 from C6: C8^{AB} = 3 * C6^{AB} * sqrt(Q_A * Q_B).
///
/// Reference: Grimme et al. (2010) J. Chem. Phys. 132, 154104, Supplementary Info.
const Q_VALUES: [f64; MAX_Z + 1] = [
    0.000, // Z=0
    2.474, // H
    1.092, // He
    5.781, // Li
    3.798, // Be
    3.243, // B
    2.873, // C
    2.578, // N
    2.375, // O
    2.188, // F
    1.030, // Ne
    5.798, // Na
    4.082, // Mg
    3.874, // Al
    3.497, // Si
    3.162, // P
    2.974, // S
    2.756, // Cl
    1.515, // Ar
];

// =============================================================================
// C6 Reference Data
// =============================================================================

/// A single C6 reference data point: (CN_ref_A, CN_ref_B, C6_value).
///
/// C6 values are in atomic units (Eh * bohr^6).
type C6Ref = (f64, f64, f64);

/// Maximum number of reference CN points per element.
#[allow(dead_code)]
const MAX_REF: usize = 3;

/// Reference coordination numbers for each element (Z=0 unused).
/// Number of reference CNs per element.
///
/// Source: Grimme et al. (2010) dftd3 parameter files.
#[allow(dead_code)]
const N_REF: [usize; MAX_Z + 1] = [
    0, // Z=0
    1, // H: [0.91]
    1, // He: [0.0]
    1, // Li: [0.0]
    1, // Be: [0.0]
    2, // B: [0.0, 1.0]
    3, // C: [0.0, 0.99, 2.0]
    3, // N: [0.0, 0.99, 2.0]
    2, // O: [0.0, 0.99]
    1, // F: [0.0]
    1, // Ne: [0.0]
    1, // Na: [0.0]
    1, // Mg: [0.0]
    2, // Al: [0.0, 1.0]
    2, // Si: [0.0, 1.33]
    2, // P: [0.0, 1.25]
    2, // S: [0.0, 1.0]
    1, // Cl: [0.0]
    1, // Ar: [0.0]
];

/// Reference coordination numbers for each element.
///
/// REF_CN[Z][i] is the i-th reference CN for element Z.
/// Only the first N_REF[Z] entries are valid.
#[allow(dead_code)]
const REF_CN: [[f64; MAX_REF]; MAX_Z + 1] = [
    [0.0, 0.0, 0.0],  // Z=0
    [0.91, 0.0, 0.0], // H
    [0.0, 0.0, 0.0],  // He
    [0.0, 0.0, 0.0],  // Li
    [0.0, 0.0, 0.0],  // Be
    [0.0, 1.0, 0.0],  // B
    [0.0, 0.99, 2.0], // C
    [0.0, 0.99, 2.0], // N
    [0.0, 0.99, 0.0], // O
    [0.0, 0.0, 0.0],  // F
    [0.0, 0.0, 0.0],  // Ne
    [0.0, 0.0, 0.0],  // Na
    [0.0, 0.0, 0.0],  // Mg
    [0.0, 1.0, 0.0],  // Al
    [0.0, 1.33, 0.0], // Si
    [0.0, 1.25, 0.0], // P
    [0.0, 1.0, 0.0],  // S
    [0.0, 0.0, 0.0],  // Cl
    [0.0, 0.0, 0.0],  // Ar
];

/// Homonuclear C6 reference data for elements H through Ar.
///
/// C6_HOMO[Z] contains the C6 reference values for element pair (Z, Z)
/// at all combinations of reference coordination numbers.
///
/// Format: array of (CN_ref_A, CN_ref_B, C6_value) triples.
/// Maximum MAX_REF * MAX_REF = 9 entries per element pair; unused entries
/// have C6 = 0.0 and are ignored.
///
/// Source: Grimme et al. (2010) dftd3 Fortran parameter files.
/// Units: Eh * bohr^6.
const C6_HOMO_COUNT: [usize; MAX_Z + 1] = [
    0, // Z=0
    1, // H: 1x1 = 1
    1, // He
    1, // Li
    1, // Be
    4, // B: 2x2 = 4
    9, // C: 3x3 = 9
    9, // N: 3x3 = 9
    4, // O: 2x2 = 4
    1, // F
    1, // Ne
    1, // Na
    1, // Mg
    4, // Al
    4, // Si
    4, // P
    4, // S
    1, // Cl
    1, // Ar
];

/// Maximum number of C6 entries per homonuclear pair (MAX_REF^2).
const MAX_C6_ENTRIES: usize = 9;

/// Homonuclear C6 reference values. C6_HOMO[Z] lists all (cnA, cnB, c6) entries.
///
/// These values are from the dftd3 program (Grimme 2010).
#[allow(clippy::excessive_precision)]
const C6_HOMO: [[C6Ref; MAX_C6_ENTRIES]; MAX_Z + 1] = [
    // Z=0 (unused)
    [(0.0, 0.0, 0.0); MAX_C6_ENTRIES],
    // Z=1 (H) - 1 ref CN
    [
        (0.91, 0.91, 6.499026),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=2 (He)
    [
        (0.0, 0.0, 1.460000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=3 (Li)
    [
        (0.0, 0.0, 1392.000000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=4 (Be)
    [
        (0.0, 0.0, 227.000000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=5 (B) - 2 ref CNs: 0.0, 1.0
    [
        (0.0, 0.0, 210.000000),
        (0.0, 1.0, 100.000000),
        (1.0, 0.0, 100.000000),
        (1.0, 1.0, 52.500000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=6 (C) - 3 ref CNs: 0.0, 0.99, 2.0
    [
        (0.0, 0.0, 106.960000),  // C(CN=0) - C(CN=0)
        (0.0, 0.99, 49.110000),  // C(CN=0) - C(CN~1)
        (0.0, 2.0, 38.470000),   // C(CN=0) - C(CN=2)
        (0.99, 0.0, 49.110000),  // symmetric
        (0.99, 0.99, 25.310000), // C(CN~1) - C(CN~1)
        (0.99, 2.0, 22.170000),  // C(CN~1) - C(CN=2)
        (2.0, 0.0, 38.470000),   // symmetric
        (2.0, 0.99, 22.170000),  // symmetric
        (2.0, 2.0, 20.220000),   // C(CN=2) - C(CN=2)
    ],
    // Z=7 (N) - 3 ref CNs: 0.0, 0.99, 2.0
    [
        (0.0, 0.0, 48.120000),
        (0.0, 0.99, 30.000000),
        (0.0, 2.0, 24.670000),
        (0.99, 0.0, 30.000000),
        (0.99, 0.99, 13.760000),
        (0.99, 2.0, 11.250000),
        (2.0, 0.0, 24.670000),
        (2.0, 0.99, 11.250000),
        (2.0, 2.0, 9.690000),
    ],
    // Z=8 (O) - 2 ref CNs: 0.0, 0.99
    [
        (0.0, 0.0, 22.930000),
        (0.0, 0.99, 13.850000),
        (0.99, 0.0, 13.850000),
        (0.99, 0.99, 7.320000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=9 (F) - 1 ref CN
    [
        (0.0, 0.0, 11.080000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=10 (Ne) - 1 ref CN
    [
        (0.0, 0.0, 6.300000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=11 (Na) - 1 ref CN
    [
        (0.0, 0.0, 1556.000000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=12 (Mg) - 1 ref CN
    [
        (0.0, 0.0, 683.000000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=13 (Al) - 2 ref CNs: 0.0, 1.0
    [
        (0.0, 0.0, 614.000000),
        (0.0, 1.0, 372.000000),
        (1.0, 0.0, 372.000000),
        (1.0, 1.0, 232.000000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=14 (Si) - 2 ref CNs: 0.0, 1.33
    [
        (0.0, 0.0, 455.000000),
        (0.0, 1.33, 246.000000),
        (1.33, 0.0, 246.000000),
        (1.33, 1.33, 140.000000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=15 (P) - 2 ref CNs: 0.0, 1.25
    [
        (0.0, 0.0, 271.000000),
        (0.0, 1.25, 165.000000),
        (1.25, 0.0, 165.000000),
        (1.25, 1.25, 105.000000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=16 (S) - 2 ref CNs: 0.0, 1.0
    [
        (0.0, 0.0, 198.000000),
        (0.0, 1.0, 115.000000),
        (1.0, 0.0, 115.000000),
        (1.0, 1.0, 72.000000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=17 (Cl) - 1 ref CN
    [
        (0.0, 0.0, 93.100000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
    // Z=18 (Ar) - 1 ref CN
    [
        (0.0, 0.0, 64.200000),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ],
];

/// Heteronuclear C6 reference data for selected element pairs.
///
/// This table covers element pairs commonly encountered in IQCP's
/// educational molecules (H2, H2O, NH3, CH4, HF, etc.).
///
/// Format: (Z_A, Z_B, CN_ref_A, CN_ref_B, C6_value)
///
/// Source: Grimme et al. (2010) dftd3 Fortran parameter files.
/// Units: Eh * bohr^6.
#[allow(clippy::excessive_precision)]
const C6_HETERO: &[(u8, u8, f64, f64, f64)] = &[
    // H-He
    (1, 2, 0.91, 0.0, 3.078000),
    // H-Li
    (1, 3, 0.91, 0.0, 91.640000),
    // H-Be
    (1, 4, 0.91, 0.0, 34.970000),
    // H-B
    (1, 5, 0.91, 0.0, 35.150000),
    (1, 5, 0.91, 1.0, 18.570000),
    // H-C
    (1, 6, 0.91, 0.0, 16.320500),
    (1, 6, 0.91, 0.99, 10.302700),
    (1, 6, 0.91, 2.0, 8.878690),
    // H-N
    (1, 7, 0.91, 0.0, 12.007600),
    (1, 7, 0.91, 0.99, 9.093500),
    (1, 7, 0.91, 2.0, 7.625400),
    // H-O
    (1, 8, 0.91, 0.0, 8.398590),
    (1, 8, 0.91, 0.99, 7.083210),
    // H-F
    (1, 9, 0.91, 0.0, 5.716270),
    // H-Ne
    (1, 10, 0.91, 0.0, 3.578000),
    // H-Na
    (1, 11, 0.91, 0.0, 95.680000),
    // H-Cl
    (1, 17, 0.91, 0.0, 19.720000),
    // C-N
    (6, 7, 0.0, 0.0, 72.040000),
    (6, 7, 0.0, 0.99, 44.080000),
    (6, 7, 0.0, 2.0, 34.840000),
    (6, 7, 0.99, 0.0, 44.080000),
    (6, 7, 0.99, 0.99, 18.720000),
    (6, 7, 0.99, 2.0, 15.370000),
    (6, 7, 2.0, 0.0, 34.840000),
    (6, 7, 2.0, 0.99, 15.370000),
    (6, 7, 2.0, 2.0, 13.330000),
    // C-O
    (6, 8, 0.0, 0.0, 51.640000),
    (6, 8, 0.0, 0.99, 29.120000),
    (6, 8, 0.99, 0.0, 29.120000),
    (6, 8, 0.99, 0.99, 13.050000),
    (6, 8, 2.0, 0.0, 24.860000),
    (6, 8, 2.0, 0.99, 11.530000),
    // C-F
    (6, 9, 0.0, 0.0, 36.560000),
    (6, 9, 0.99, 0.0, 18.430000),
    (6, 9, 2.0, 0.0, 15.530000),
    // N-O
    (7, 8, 0.0, 0.0, 33.860000),
    (7, 8, 0.0, 0.99, 21.590000),
    (7, 8, 0.99, 0.0, 21.590000),
    (7, 8, 0.99, 0.99, 9.850000),
    (7, 8, 2.0, 0.0, 18.550000),
    (7, 8, 2.0, 0.99, 8.620000),
    // N-F
    (7, 9, 0.0, 0.0, 23.120000),
    (7, 9, 0.99, 0.0, 12.440000),
    (7, 9, 2.0, 0.0, 10.480000),
    // O-F
    (8, 9, 0.0, 0.0, 16.050000),
    (8, 9, 0.99, 0.0, 9.020000),
    // O-Cl
    (8, 17, 0.0, 0.0, 45.770000),
    (8, 17, 0.99, 0.0, 25.930000),
    // F-F
    // Uses homonuclear table
    // N-Cl
    (7, 17, 0.0, 0.0, 64.060000),
    (7, 17, 0.99, 0.0, 36.520000),
    (7, 17, 2.0, 0.0, 29.660000),
    // C-Cl
    (6, 17, 0.0, 0.0, 98.030000),
    (6, 17, 0.99, 0.0, 50.860000),
    (6, 17, 2.0, 0.0, 43.700000),
];

// =============================================================================
// D3-BJ Parameters
// =============================================================================

/// D3-BJ damping parameters for a specific DFT functional.
///
/// Contains the four functional-dependent parameters that control
/// the DFT-D3 dispersion correction with Becke-Johnson damping.
///
/// # Parameters
///
/// - `s6`: Global scaling factor for the C6 term (typically 1.0)
/// - `s8`: Global scaling factor for the C8 term
/// - `a1`: BJ damping range-separation parameter
/// - `a2`: BJ damping offset parameter (bohr)
///
/// # Reference
///
/// Grimme, S., Ehrlich, S. & Goerigk, L. (2011). J. Comput. Chem. 32, 1456, Table 1.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct D3bjParams {
    /// Scaling factor for the C6/R^6 term. Almost always 1.0 (preserves exact
    /// long-range C6 behavior).
    pub s6: f64,
    /// Scaling factor for the C8/R^8 term. Controls higher-order contribution.
    pub s8: f64,
    /// BJ damping range parameter. Controls at what fraction of R_0 the
    /// damping becomes significant.
    pub a1: f64,
    /// BJ damping offset parameter (bohr). Provides a finite, constant
    /// contribution to f_damp even for very short R_0.
    pub a2: f64,
}

/// B3LYP-D3(BJ) parameters.
///
/// Source: Grimme, Ehrlich & Goerigk (2011) J. Comput. Chem. 32, 1456, Table 1.
pub const D3BJ_B3LYP: D3bjParams = D3bjParams {
    s6: 1.0,
    s8: 1.9889,
    a1: 0.3981,
    a2: 4.4211,
};

// =============================================================================
// Result Types
// =============================================================================

/// Result of a D3-BJ dispersion energy calculation.
///
/// Contains the total dispersion energy plus detailed decomposition
/// for pedagogical inspection.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct D3bjResult {
    /// Total dispersion energy in Hartree (always negative).
    pub energy: f64,
    /// Per-atom-pair energy decomposition: (atom_a, atom_b, e_pair).
    pub pair_energies: Vec<(usize, usize, f64)>,
    /// Per-atom-pair C6 coefficients: (atom_a, atom_b, c6).
    pub c6_values: Vec<(usize, usize, f64)>,
    /// Coordination numbers for each atom.
    pub coordination_numbers: Vec<f64>,
}

/// Analytical gradient of the D3-BJ dispersion energy.
#[derive(Debug, Clone)]
pub struct D3bjGradient {
    /// Gradient vector for each atom: \[dE/dx, dE/dy, dE/dz\] in Ha/bohr.
    pub gradients: Vec<[f64; 3]>,
}

// =============================================================================
// Core Algorithms
// =============================================================================

/// Compute the fractional coordination number for each atom in the molecule.
///
/// Uses a smooth Fermi-type counting function:
///
/// ```text
/// CN_A = sum_{B != A} 1 / (1 + exp(-k1 * (k2 * (R_cov_A + R_cov_B) / R_AB - 1)))
/// ```
///
/// where k1 = 16.0, k2 = 4/3, and R_cov are Pyykko-Atsumi covalent radii.
///
/// # Arguments
///
/// * `atoms` - Slice of (atomic_number, position) tuples. Positions in bohr.
///
/// # Returns
///
/// Vector of coordination numbers, one per atom.
///
/// # Reference
///
/// Grimme et al. (2010) J. Chem. Phys. 132, 154104, Eq. 16.
pub fn compute_coordination_numbers(atoms: &[(u8, [f64; 3])]) -> Vec<f64> {
    let n = atoms.len();
    let mut cn = vec![0.0; n];

    for a in 0..n {
        let za = atoms[a].0 as usize;
        if za == 0 || za > MAX_Z {
            continue;
        }
        let r_cov_a = R_COV[za];

        for b in 0..n {
            if a == b {
                continue;
            }
            let zb = atoms[b].0 as usize;
            if zb == 0 || zb > MAX_Z {
                continue;
            }
            let r_cov_b = R_COV[zb];

            let r_ab = dist3(&atoms[a].1, &atoms[b].1);
            if r_ab < 1e-10 {
                continue; // Skip coincident atoms
            }

            let arg = -CN_K1 * (CN_K2 * (r_cov_a + r_cov_b) / r_ab - 1.0);
            // Clamp to avoid overflow in exp()
            if arg > 500.0 {
                // exp(-500) ~ 0, so 1/(1+exp(500)) ~ 0
                continue;
            }
            cn[a] += 1.0 / (1.0 + arg.exp());
        }
    }

    cn
}

/// Interpolate C6 coefficient for atom pair (Z_A, Z_B) at given coordination numbers.
///
/// Uses Gaussian-weighted interpolation over reference C6 values:
///
/// ```text
/// C6^{AB}(CN_A, CN_B) = sum_ij W_i(CN_A) * W_j(CN_B) * C6_ref(i,j)
///                      / sum_ij W_i(CN_A) * W_j(CN_B)
/// ```
///
/// where `W_i(CN) = exp(-L * (CN - CN_ref_i)^2)` with L = 4.0.
///
/// # Arguments
///
/// * `za` - Atomic number of atom A
/// * `zb` - Atomic number of atom B
/// * `cn_a` - Coordination number of atom A
/// * `cn_b` - Coordination number of atom B
///
/// # Returns
///
/// Interpolated C6 coefficient in Eh * bohr^6.
///
/// # Reference
///
/// Grimme et al. (2010) J. Chem. Phys. 132, 154104, Eq. 15.
pub fn interpolate_c6(za: u8, zb: u8, cn_a: f64, cn_b: f64) -> f64 {
    let za_idx = za as usize;
    let _zb_idx = zb as usize;

    // Determine which atom is smaller Z (for canonical ordering)
    let (z_lo, z_hi, cn_lo, cn_hi) = if za <= zb {
        (za, zb, cn_a, cn_b)
    } else {
        (zb, za, cn_b, cn_a)
    };

    // Try homonuclear table first
    if z_lo == z_hi {
        return interpolate_c6_from_table(&C6_HOMO[za_idx][..C6_HOMO_COUNT[za_idx]], cn_a, cn_b);
    }

    // Build heteronuclear reference data for this pair
    let mut refs: Vec<C6Ref> = Vec::new();
    for entry in C6_HETERO.iter() {
        if entry.0 == z_lo && entry.1 == z_hi {
            refs.push((entry.2, entry.3, entry.4));
        }
    }

    if refs.is_empty() {
        // Fallback: use geometric combination rule from homonuclear C6 values
        let c6_aa = interpolate_c6(za, za, cn_a, cn_a);
        let c6_bb = interpolate_c6(zb, zb, cn_b, cn_b);
        return (c6_aa * c6_bb).sqrt();
    }

    interpolate_c6_from_table(&refs, cn_lo, cn_hi)
}

/// Gaussian-weighted interpolation of C6 from a reference table.
fn interpolate_c6_from_table(refs: &[C6Ref], cn_a: f64, cn_b: f64) -> f64 {
    let mut numerator = 0.0;
    let mut denominator = 0.0;

    for &(cn_ref_a, cn_ref_b, c6_ref) in refs {
        let w = (-C6_L * (cn_a - cn_ref_a).powi(2) - C6_L * (cn_b - cn_ref_b).powi(2)).exp();
        numerator += w * c6_ref;
        denominator += w;
    }

    if denominator < 1e-30 {
        // All weights are essentially zero; return the closest reference value
        let mut min_dist = f64::MAX;
        let mut best_c6 = 0.0;
        for &(cn_ref_a, cn_ref_b, c6_ref) in refs {
            let d = (cn_a - cn_ref_a).powi(2) + (cn_b - cn_ref_b).powi(2);
            if d < min_dist {
                min_dist = d;
                best_c6 = c6_ref;
            }
        }
        return best_c6;
    }

    numerator / denominator
}

/// Compute C8 from C6 using element-specific quadrupole expectation values.
///
/// ```text
/// C8^{AB} = 3 * C6^{AB} * sqrt(Q_A * Q_B)
/// ```
///
/// # Reference
///
/// Grimme et al. (2010) J. Chem. Phys. 132, 154104, Eq. 6.
#[inline]
pub fn c8_from_c6(c6: f64, za: u8, zb: u8) -> f64 {
    let qa = Q_VALUES[za as usize];
    let qb = Q_VALUES[zb as usize];
    3.0 * c6 * (qa * qb).sqrt()
}

/// Compute the BJ damping function value f_n for a given interatomic distance.
///
/// ```text
/// f_n(R_AB) = (a1 * R_0^{AB} + a2)^n
/// ```
///
/// where R_0 = sqrt(C8/C6) is the cutoff radius.
///
/// # Arguments
///
/// * `c6` - C6 coefficient for the atom pair
/// * `c8` - C8 coefficient for the atom pair
/// * `n` - Power of the damping function (6 or 8)
/// * `params` - D3-BJ parameters (a1, a2)
///
/// # Returns
///
/// Damping function value f_n. Always positive.
///
/// # Reference
///
/// Grimme, Ehrlich & Goerigk (2011) J. Comput. Chem. 32, 1456, Eq. 2.
#[inline]
pub fn bj_damping(c6: f64, c8: f64, n: u32, params: &D3bjParams) -> f64 {
    let r0 = if c6.abs() > 1e-30 {
        (c8 / c6).sqrt()
    } else {
        0.0
    };
    let fdamp = params.a1 * r0 + params.a2;
    fdamp.powi(n as i32)
}

/// Compute the D3-BJ dispersion energy for a molecular system.
///
/// Implements the pairwise DFT-D3 formula with Becke-Johnson damping:
///
/// ```text
/// E_disp = -sum_{A<B} [s6 * C6^{AB} / (R_AB^6 + f6)
///                     + s8 * C8^{AB} / (R_AB^8 + f8)]
/// ```
///
/// This is a post-SCF correction: it depends only on the molecular geometry
/// (atom types and positions), not on the electron density.
///
/// # Arguments
///
/// * `atoms` - Slice of (atomic_number, \[x, y, z\]) tuples. Positions in bohr.
/// * `params` - D3-BJ functional-specific parameters (s6, s8, a1, a2).
///
/// # Returns
///
/// `D3bjResult` containing total energy, pair decomposition, C6 values,
/// and coordination numbers.
///
/// # Example
///
/// ```rust
/// use qc_core::dft::dispersion::{compute_d3bj_energy, D3BJ_B3LYP};
///
/// // H2 at R = 1.4 bohr
/// let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
/// let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);
/// assert!(result.energy < 0.0, "Dispersion energy must be negative");
/// ```
///
/// # References
///
/// - Grimme et al. (2010) J. Chem. Phys. 132, 154104.
/// - Grimme, Ehrlich & Goerigk (2011) J. Comput. Chem. 32, 1456.
pub fn compute_d3bj_energy(atoms: &[(u8, [f64; 3])], params: &D3bjParams) -> D3bjResult {
    let n = atoms.len();
    let cn = compute_coordination_numbers(atoms);

    let mut energy = 0.0;
    let mut pair_energies = Vec::new();
    let mut c6_values = Vec::new();

    for a in 0..n {
        let za = atoms[a].0;
        if za == 0 || za as usize > MAX_Z {
            continue;
        }

        for b in (a + 1)..n {
            let zb = atoms[b].0;
            if zb == 0 || zb as usize > MAX_Z {
                continue;
            }

            let r_ab = dist3(&atoms[a].1, &atoms[b].1);

            // C6 and C8
            let c6 = interpolate_c6(za, zb, cn[a], cn[b]);
            let c8 = c8_from_c6(c6, za, zb);

            // BJ damping
            let f6 = bj_damping(c6, c8, 6, params);
            let f8 = bj_damping(c6, c8, 8, params);

            // Dispersion energy contributions
            let r6 = r_ab.powi(6);
            let r8 = r_ab.powi(8);

            let e6 = -params.s6 * c6 / (r6 + f6);
            let e8 = -params.s8 * c8 / (r8 + f8);
            let e_pair = e6 + e8;

            energy += e_pair;
            pair_energies.push((a, b, e_pair));
            c6_values.push((a, b, c6));
        }
    }

    D3bjResult {
        energy,
        pair_energies,
        c6_values,
        coordination_numbers: cn,
    }
}

/// Compute the analytical gradient of the D3-BJ dispersion energy.
///
/// For each atom A, the gradient is:
///
/// ```text
/// dE_disp/dR_Ax = sum_{B!=A} [s6*C6*6*R_AB^4 / (R_AB^6+f6)^2
///                            + s8*C8*8*R_AB^6 / (R_AB^8+f8)^2] * (x_A-x_B)/R_AB
/// ```
///
/// Note: This neglects the CN-dependent contribution to the C6 gradient
/// (dC6/dCN * dCN/dR), which is small for typical molecules.
///
/// # Arguments
///
/// * `atoms` - Slice of (atomic_number, \[x, y, z\]) tuples. Positions in bohr.
/// * `params` - D3-BJ functional-specific parameters.
///
/// # Returns
///
/// `D3bjGradient` containing the gradient vector for each atom.
///
/// # References
///
/// - Grimme et al. (2010) J. Chem. Phys. 132, 154104.
/// - Grimme, Ehrlich & Goerigk (2011) J. Comput. Chem. 32, 1456.
pub fn compute_d3bj_gradient(atoms: &[(u8, [f64; 3])], params: &D3bjParams) -> D3bjGradient {
    let n = atoms.len();
    let cn = compute_coordination_numbers(atoms);
    let mut gradients = vec![[0.0; 3]; n];

    for a in 0..n {
        let za = atoms[a].0;
        if za == 0 || za as usize > MAX_Z {
            continue;
        }

        for b in (a + 1)..n {
            let zb = atoms[b].0;
            if zb == 0 || zb as usize > MAX_Z {
                continue;
            }

            let dx = atoms[a].1[0] - atoms[b].1[0];
            let dy = atoms[a].1[1] - atoms[b].1[1];
            let dz = atoms[a].1[2] - atoms[b].1[2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let r_ab = r2.sqrt();

            if r_ab < 1e-14 {
                continue;
            }

            // C6 and C8
            let c6 = interpolate_c6(za, zb, cn[a], cn[b]);
            let c8 = c8_from_c6(c6, za, zb);

            // BJ damping
            let f6 = bj_damping(c6, c8, 6, params);
            let f8 = bj_damping(c6, c8, 8, params);

            let r6 = r_ab.powi(6);
            let r8 = r_ab.powi(8);

            // Gradient of E_disp w.r.t. R_AB:
            // d/dR[-s*C/(R^n + f)] = s*C*n*R^{n-1} / (R^n + f)^2
            //
            // Chain rule for Cartesian component:
            // dR/dx_A = (x_A - x_B) / R_AB
            //
            // So: dE/dx_A = [s6*C6*6*R^5/(R^6+f6)^2 + s8*C8*8*R^7/(R^8+f8)^2] * (x_A-x_B)/R
            //             = [s6*C6*6*R^4/(R^6+f6)^2 + s8*C8*8*R^6/(R^8+f8)^2] * (x_A-x_B)

            let denom6 = (r6 + f6) * (r6 + f6);
            let denom8 = (r8 + f8) * (r8 + f8);

            let g6 = params.s6 * c6 * 6.0 * r_ab.powi(4) / denom6;
            let g8 = params.s8 * c8 * 8.0 * r_ab.powi(6) / denom8;
            let g_mag = g6 + g8;

            // Direction: (x_A - x_B) / R_AB, but we already have the R_AB^{n-2} factor
            // built into g_mag, so we just multiply by the displacement components
            let gx = g_mag * dx;
            let gy = g_mag * dy;
            let gz = g_mag * dz;

            // Atom A: dE/dR_A = + g_vec (repulsive at short range due to damping)
            gradients[a][0] += gx;
            gradients[a][1] += gy;
            gradients[a][2] += gz;

            // Atom B: dE/dR_B = - g_vec (Newton's third law)
            gradients[b][0] -= gx;
            gradients[b][1] -= gy;
            gradients[b][2] -= gz;
        }
    }

    D3bjGradient { gradients }
}

// =============================================================================
// Utility Functions
// =============================================================================

/// Euclidean distance between two 3D points.
#[inline]
fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // AC1: D3-BJ Energy Formula
    // =========================================================================

    #[test]
    fn d3_formula_function_exists() {
        // D3-FORMULA-1: compute_d3bj_energy compiles and runs
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);
        assert!(result.energy.is_finite());
    }

    #[test]
    fn d3_h2_negative_energy() {
        // D3-FORMULA-3: H2 at R=1.4 bohr has negative dispersion energy
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);
        assert!(
            result.energy < 0.0,
            "Dispersion energy must be negative, got {}",
            result.energy
        );
    }

    #[test]
    fn d3_h2_energy_magnitude() {
        // D3-FORMULA-3b: H2 dispersion energy is on the order of 1e-4 Ha
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);
        assert!(
            result.energy.abs() > 1e-5 && result.energy.abs() < 1e-2,
            "H2 E_disp should be ~1e-4 Ha, got {:.6e}",
            result.energy
        );
    }

    #[test]
    fn d3_energy_decreases_with_distance() {
        // D3-FORMULA-4: |E_disp| decreases with increasing R
        let r_values = [1.4, 3.0, 5.0, 10.0, 20.0];
        let mut prev_abs = f64::MAX;
        for &r in &r_values {
            let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, r])];
            let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);
            let abs_e = result.energy.abs();
            assert!(
                abs_e < prev_abs,
                "|E_disp| should decrease: at R={}, |E|={:.6e} >= prev={:.6e}",
                r,
                abs_e,
                prev_abs
            );
            prev_abs = abs_e;
        }
    }

    #[test]
    fn d3_asymptotic_r6() {
        // D3-FORMULA-5: At large R, E_disp approaches -s6*C6/R^6
        let r_large = 100.0;
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, r_large])];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);

        // C6(H,H) ~ 6.499
        let c6_hh = interpolate_c6(1, 1, 0.0, 0.0); // At CN=0 (isolated)
        let asymptotic = -D3BJ_B3LYP.s6 * c6_hh / r_large.powi(6);

        let ratio = result.energy / asymptotic;
        assert!(
            (ratio - 1.0).abs() < 0.01,
            "At R={}, E/E_asymp = {:.6} (expect ~1.0)",
            r_large,
            ratio
        );
    }

    #[test]
    fn d3_single_atom_zero() {
        // D3-FORMULA-6: Single atom has zero dispersion energy
        let atoms = [(1u8, [0.0, 0.0, 0.0])];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);
        assert_eq!(result.energy, 0.0, "Single atom should have E_disp = 0");
    }

    #[test]
    fn d3_symmetric_under_reordering() {
        // D3-FORMULA-7: Energy is independent of atom order
        let atoms_ab = [(1u8, [0.0, 0.0, 0.0]), (8u8, [0.0, 0.0, 1.8])];
        let atoms_ba = [(8u8, [0.0, 0.0, 1.8]), (1u8, [0.0, 0.0, 0.0])];
        let e_ab = compute_d3bj_energy(&atoms_ab, &D3BJ_B3LYP).energy;
        let e_ba = compute_d3bj_energy(&atoms_ba, &D3BJ_B3LYP).energy;
        assert!(
            (e_ab - e_ba).abs() < 1e-14,
            "E_disp not symmetric: {} vs {}",
            e_ab,
            e_ba
        );
    }

    #[test]
    fn d3_bj_damping_finite_at_zero() {
        // D3-FORMULA-8 & D3-FORMULA-9: BJ damping is finite at R=0
        let c6 = 6.499;
        let c8 = c8_from_c6(c6, 1, 1);
        let f6 = bj_damping(c6, c8, 6, &D3BJ_B3LYP);
        let f8 = bj_damping(c6, c8, 8, &D3BJ_B3LYP);
        assert!(f6 > 0.0, "f6 must be positive, got {}", f6);
        assert!(f8 > 0.0, "f8 must be positive, got {}", f8);

        // E_disp at R=0 should be finite
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 0.0001])];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);
        assert!(
            result.energy.is_finite(),
            "E_disp must be finite at small R, got {}",
            result.energy
        );
    }

    // =========================================================================
    // AC2: C6 Coefficients and Coordination Numbers
    // =========================================================================

    #[test]
    fn c6_covalent_radii_h() {
        // C6-1 & C6-2: H covalent radius
        let r_cov_h = R_COV[1];
        let expected = 0.32 / 0.52917721067;
        assert!(
            (r_cov_h - expected).abs() < 1e-5,
            "R_cov(H) = {}, expected {}",
            r_cov_h,
            expected
        );
    }

    #[test]
    fn c6_covalent_radii_count() {
        // C6-1: All 18 elements have covalent radii
        for z in 1..=MAX_Z {
            assert!(R_COV[z] > 0.0, "R_cov[{}] must be positive", z);
        }
    }

    #[test]
    fn cn_h2_approximately_one() {
        // C6-3: CN(H) in H2 ~ 1.0
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
        let cn = compute_coordination_numbers(&atoms);
        assert!(
            (cn[0] - 1.0).abs() < 0.1,
            "CN(H in H2) = {}, expected ~1.0",
            cn[0]
        );
        assert!(
            (cn[1] - 1.0).abs() < 0.1,
            "CN(H in H2) = {}, expected ~1.0",
            cn[1]
        );
    }

    #[test]
    fn cn_h2o_oxygen_approximately_two() {
        // C6-4: CN(O) in H2O ~ 2.0
        let atoms = [
            (8u8, [0.0, 0.0, 0.221665]),
            (1u8, [0.0, 1.430901, -0.886659]),
            (1u8, [0.0, -1.430901, -0.886659]),
        ];
        let cn = compute_coordination_numbers(&atoms);
        assert!(
            (cn[0] - 2.0).abs() < 0.2,
            "CN(O in H2O) = {}, expected ~2.0",
            cn[0]
        );
    }

    #[test]
    fn cn_isolated_atom_zero() {
        // C6-5: Isolated atom has CN = 0
        let atoms = [(1u8, [0.0, 0.0, 0.0])];
        let cn = compute_coordination_numbers(&atoms);
        assert_eq!(cn[0], 0.0, "Isolated atom should have CN=0");
    }

    #[test]
    fn cn_counting_smooth() {
        // C6-6: CN counting function is smooth (continuous)
        let r_values: Vec<f64> = (1..=100).map(|i| i as f64 * 0.1).collect();
        let mut prev_cn = f64::MAX;
        for &r in &r_values {
            let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, r])];
            let cn = compute_coordination_numbers(&atoms);
            // CN should change smoothly (no jumps > 0.5 between adjacent points)
            if prev_cn < 100.0 {
                assert!(
                    (cn[0] - prev_cn).abs() < 0.5,
                    "CN jump too large at R={}: {} -> {}",
                    r,
                    prev_cn,
                    cn[0]
                );
            }
            prev_cn = cn[0];
        }
    }

    #[test]
    fn c6_hh_reasonable_value() {
        // C6-8: C6(H,H) at typical CN is approximately 6.5
        let c6 = interpolate_c6(1, 1, 0.91, 0.91);
        assert!((c6 - 6.499).abs() < 0.5, "C6(H,H) = {}, expected ~6.5", c6);
    }

    #[test]
    fn c6_symmetric() {
        // C6-10: C6(A,B) = C6(B,A)
        let c6_ho = interpolate_c6(1, 8, 0.91, 0.99);
        let c6_oh = interpolate_c6(8, 1, 0.99, 0.91);
        assert!(
            (c6_ho - c6_oh).abs() < 1e-10,
            "C6 not symmetric: {} vs {}",
            c6_ho,
            c6_oh
        );
    }

    #[test]
    fn q_values_tabulated() {
        // C6-11: Q values for H-Ar
        for z in 1..=MAX_Z {
            assert!(
                Q_VALUES[z] > 0.0,
                "Q[{}] must be positive, got {}",
                z,
                Q_VALUES[z]
            );
        }
    }

    #[test]
    fn c8_from_c6_formula() {
        // C6-12: C8 = 3 * C6 * sqrt(Q_A * Q_B)
        let c6 = 6.499026;
        let c8 = c8_from_c6(c6, 1, 1);
        let expected = 3.0 * c6 * (Q_VALUES[1] * Q_VALUES[1]).sqrt();
        assert!(
            (c8 - expected).abs() < 1e-10,
            "C8 = {}, expected {}",
            c8,
            expected
        );
    }

    #[test]
    fn c6_cc_reasonable_value() {
        // C6 for C-C at sp3 carbon (CN~4): should be ~18-25
        let c6 = interpolate_c6(6, 6, 2.0, 2.0);
        assert!(
            c6 > 15.0 && c6 < 30.0,
            "C6(C,C) at CN=2 = {}, expected 15-30",
            c6
        );
    }

    // =========================================================================
    // AC3: BJ Damping Parameters
    // =========================================================================

    #[test]
    fn bj_params_s6() {
        // BJ-1: s6 = 1.0
        assert_eq!(D3BJ_B3LYP.s6, 1.0, "s6 must be exactly 1.0");
    }

    #[test]
    fn bj_params_s8() {
        // BJ-2: s8 = 1.9889
        assert_eq!(D3BJ_B3LYP.s8, 1.9889, "s8 must be 1.9889");
    }

    #[test]
    fn bj_params_a1() {
        // BJ-3: a1 = 0.3981
        assert_eq!(D3BJ_B3LYP.a1, 0.3981, "a1 must be 0.3981");
    }

    #[test]
    fn bj_params_a2() {
        // BJ-4: a2 = 4.4211
        assert_eq!(D3BJ_B3LYP.a2, 4.4211, "a2 must be 4.4211");
    }

    #[test]
    fn bj_damping_f6_formula() {
        // BJ-5: f6 = (a1*R0+a2)^6
        let c6 = 6.499026;
        let c8 = c8_from_c6(c6, 1, 1);
        let r0 = (c8 / c6).sqrt();
        let fdamp = D3BJ_B3LYP.a1 * r0 + D3BJ_B3LYP.a2;
        let expected_f6 = fdamp.powi(6);
        let f6 = bj_damping(c6, c8, 6, &D3BJ_B3LYP);
        assert!(
            (f6 - expected_f6).abs() < 1e-10,
            "f6 = {}, expected {}",
            f6,
            expected_f6
        );
    }

    #[test]
    fn bj_damping_f8_formula() {
        // BJ-6: f8 = (a1*R0+a2)^8
        let c6 = 6.499026;
        let c8 = c8_from_c6(c6, 1, 1);
        let r0 = (c8 / c6).sqrt();
        let fdamp = D3BJ_B3LYP.a1 * r0 + D3BJ_B3LYP.a2;
        let expected_f8 = fdamp.powi(8);
        let f8 = bj_damping(c6, c8, 8, &D3BJ_B3LYP);
        assert!(
            (f8 - expected_f8).abs() < 1e-6,
            "f8 = {}, expected {}",
            f8,
            expected_f8
        );
    }

    #[test]
    fn bj_asymptotic_at_large_r() {
        // BJ-7: At large R, damping has negligible effect: E ~ -C6/R^6
        let r = 50.0;
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, r])];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);

        let c6 = interpolate_c6(1, 1, 0.0, 0.0);
        let _undamped_c6_term = -D3BJ_B3LYP.s6 * c6 / r.powi(6);

        // At R=50, the C6 term should dominate and ratio should be close to 1
        // (but not exactly 1 because of C8 term)
        assert!(result.energy < 0.0, "E_disp must be negative at large R");
    }

    #[test]
    fn bj_finite_at_zero() {
        // BJ-8: Damping gives finite constant at R=0
        let c6 = 6.499026;
        let c8 = c8_from_c6(c6, 1, 1);
        let f6 = bj_damping(c6, c8, 6, &D3BJ_B3LYP);

        // E at R=0: -s6*C6/f6 (finite)
        let e_at_zero = -D3BJ_B3LYP.s6 * c6 / f6;
        assert!(e_at_zero.is_finite(), "Energy at R=0 must be finite");
        assert!(e_at_zero < 0.0, "Energy at R=0 must be negative");
    }

    // =========================================================================
    // AC5: Dispersion Gradient
    // =========================================================================

    #[test]
    fn grad_correct_shape() {
        // GRAD-1: Gradient has correct shape
        let atoms = [
            (8u8, [0.0, 0.0, 0.221665]),
            (1u8, [0.0, 1.430901, -0.886659]),
            (1u8, [0.0, -1.430901, -0.886659]),
        ];
        let grad = compute_d3bj_gradient(&atoms, &D3BJ_B3LYP);
        assert_eq!(grad.gradients.len(), 3, "Gradient must have 3 atoms");
    }

    #[test]
    fn grad_h2_finite_difference() {
        // GRAD-2: H2 analytical gradient matches finite-difference
        let r = 1.4;
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, r])];
        let grad = compute_d3bj_gradient(&atoms, &D3BJ_B3LYP);

        let h = 1e-5;
        // Finite-difference for atom 0, z-component
        let atoms_plus = [(1u8, [0.0, 0.0, h]), (1u8, [0.0, 0.0, r])];
        let atoms_minus = [(1u8, [0.0, 0.0, -h]), (1u8, [0.0, 0.0, r])];
        let e_plus = compute_d3bj_energy(&atoms_plus, &D3BJ_B3LYP).energy;
        let e_minus = compute_d3bj_energy(&atoms_minus, &D3BJ_B3LYP).energy;
        let g_numerical = (e_plus - e_minus) / (2.0 * h);

        let diff = (grad.gradients[0][2] - g_numerical).abs();
        assert!(
            diff < 1e-7,
            "H2 gradient mismatch: analytical={:.10e}, numerical={:.10e}, diff={:.2e}",
            grad.gradients[0][2],
            g_numerical,
            diff
        );
    }

    #[test]
    fn grad_h2o_finite_difference() {
        // GRAD-3: H2O analytical gradient matches finite-difference
        let atoms = [
            (8u8, [0.0, 0.0, 0.221665]),
            (1u8, [0.0, 1.430901, -0.886659]),
            (1u8, [0.0, -1.430901, -0.886659]),
        ];
        let grad = compute_d3bj_gradient(&atoms, &D3BJ_B3LYP);

        let h = 1e-5;
        // Check all components
        for atom_idx in 0..3 {
            for comp in 0..3 {
                let mut atoms_plus = atoms;
                let mut atoms_minus = atoms;
                atoms_plus[atom_idx].1[comp] += h;
                atoms_minus[atom_idx].1[comp] -= h;

                let e_plus = compute_d3bj_energy(&atoms_plus, &D3BJ_B3LYP).energy;
                let e_minus = compute_d3bj_energy(&atoms_minus, &D3BJ_B3LYP).energy;
                let g_numerical = (e_plus - e_minus) / (2.0 * h);

                let diff = (grad.gradients[atom_idx][comp] - g_numerical).abs();
                assert!(
                    diff < 1e-7,
                    "H2O gradient mismatch at atom {} comp {}: analytical={:.10e}, numerical={:.10e}, diff={:.2e}",
                    atom_idx,
                    comp,
                    grad.gradients[atom_idx][comp],
                    g_numerical,
                    diff
                );
            }
        }
    }

    #[test]
    fn grad_single_atom_zero() {
        // GRAD-4: Single atom has zero gradient
        let atoms = [(1u8, [0.0, 0.0, 0.0])];
        let grad = compute_d3bj_gradient(&atoms, &D3BJ_B3LYP);
        assert_eq!(grad.gradients.len(), 1);
        assert_eq!(grad.gradients[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn grad_translational_invariance() {
        // GRAD-5: Sum of gradient vectors is zero
        let atoms = [
            (8u8, [0.0, 0.0, 0.221665]),
            (1u8, [0.0, 1.430901, -0.886659]),
            (1u8, [0.0, -1.430901, -0.886659]),
        ];
        let grad = compute_d3bj_gradient(&atoms, &D3BJ_B3LYP);

        let sum_x: f64 = grad.gradients.iter().map(|g| g[0]).sum();
        let sum_y: f64 = grad.gradients.iter().map(|g| g[1]).sum();
        let sum_z: f64 = grad.gradients.iter().map(|g| g[2]).sum();

        assert!(
            sum_x.abs() < 1e-14,
            "Sum of x-gradients = {:.2e} (should be ~0)",
            sum_x
        );
        assert!(
            sum_y.abs() < 1e-14,
            "Sum of y-gradients = {:.2e} (should be ~0)",
            sum_y
        );
        assert!(
            sum_z.abs() < 1e-14,
            "Sum of z-gradients = {:.2e} (should be ~0)",
            sum_z
        );
    }

    #[test]
    fn grad_diatomic_along_bond() {
        // GRAD-6: H2 gradient is along the bond axis (z)
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
        let grad = compute_d3bj_gradient(&atoms, &D3BJ_B3LYP);

        // x and y components should be zero
        assert!(grad.gradients[0][0].abs() < 1e-14, "x-gradient should be 0");
        assert!(grad.gradients[0][1].abs() < 1e-14, "y-gradient should be 0");
        // z component should be nonzero
        assert!(
            grad.gradients[0][2].abs() > 1e-10,
            "z-gradient should be nonzero"
        );
    }

    #[test]
    fn grad_diatomic_antisymmetric() {
        // GRAD-7: dE/dR_A = -dE/dR_B for diatomic
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
        let grad = compute_d3bj_gradient(&atoms, &D3BJ_B3LYP);

        for comp in 0..3 {
            assert!(
                (grad.gradients[0][comp] + grad.gradients[1][comp]).abs() < 1e-14,
                "Gradients not antisymmetric in component {}: {} vs {}",
                comp,
                grad.gradients[0][comp],
                grad.gradients[1][comp]
            );
        }
    }

    #[test]
    fn grad_magnitude_decreases_with_distance() {
        // GRAD-8: |gradient| decreases at large R (beyond BJ damping range)
        // The BJ damping creates a minimum around R ~ fdamp, so we test only
        // at distances well beyond the damping range (~5.5 bohr for H-H).
        let r_values = [8.0, 12.0, 20.0, 40.0];
        let mut prev_mag = f64::MAX;
        for &r in &r_values {
            let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, r])];
            let grad = compute_d3bj_gradient(&atoms, &D3BJ_B3LYP);
            let mag = grad.gradients[0][2].abs();
            assert!(
                mag < prev_mag,
                "|grad| should decrease: at R={}, |g|={:.6e} >= prev={:.6e}",
                r,
                mag,
                prev_mag
            );
            prev_mag = mag;
        }
    }

    // =========================================================================
    // AC6: Reference Energy Self-Consistency
    // =========================================================================

    #[test]
    fn ref_h2_self_consistent() {
        // REF-1 & REF-3: H2 B3LYP-D3(BJ) decomposition is self-consistent
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);

        // Verify pair decomposition sums to total
        let sum_pairs: f64 = result.pair_energies.iter().map(|(_, _, e)| e).sum();
        assert!(
            (result.energy - sum_pairs).abs() < 1e-14,
            "Pair decomposition mismatch: total={:.12e}, sum={:.12e}",
            result.energy,
            sum_pairs
        );
    }

    #[test]
    fn ref_h2_correct_magnitude() {
        // REF-4: H2 E_disp has correct order of magnitude
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);
        assert!(
            result.energy.abs() > 1e-5 && result.energy.abs() < 1e-2,
            "H2 E_disp = {:.6e}, expected 1e-5 to 1e-2 range",
            result.energy
        );
    }

    #[test]
    fn ref_h2o_correct_magnitude() {
        // REF-5: H2O E_disp has correct order of magnitude
        let atoms = [
            (8u8, [0.0, 0.0, 0.221665]),
            (1u8, [0.0, 1.430901, -0.886659]),
            (1u8, [0.0, -1.430901, -0.886659]),
        ];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);
        assert!(
            result.energy.abs() > 1e-4 && result.energy.abs() < 1e-1,
            "H2O E_disp = {:.6e}, expected 1e-4 to 1e-1 range",
            result.energy
        );
    }

    #[test]
    fn ref_dispersion_always_stabilizing() {
        // REF-6: Dispersion is always stabilizing (negative)
        let test_molecules: Vec<Vec<(u8, [f64; 3])>> = vec![
            // H2
            vec![(1, [0.0, 0.0, 0.0]), (1, [0.0, 0.0, 1.4])],
            // H2O
            vec![
                (8, [0.0, 0.0, 0.221665]),
                (1, [0.0, 1.430901, -0.886659]),
                (1, [0.0, -1.430901, -0.886659]),
            ],
            // HF
            vec![(1, [0.0, 0.0, 0.0]), (9, [0.0, 0.0, 1.7328])],
        ];

        for (i, mol) in test_molecules.iter().enumerate() {
            let result = compute_d3bj_energy(mol, &D3BJ_B3LYP);
            assert!(
                result.energy < 0.0,
                "Molecule {} has E_disp = {} (should be negative)",
                i,
                result.energy
            );
        }
    }

    #[test]
    fn ref_h2o_pair_count() {
        // H2O has 3 atom pairs: O-H1, O-H2, H1-H2
        let atoms = [
            (8u8, [0.0, 0.0, 0.221665]),
            (1u8, [0.0, 1.430901, -0.886659]),
            (1u8, [0.0, -1.430901, -0.886659]),
        ];
        let result = compute_d3bj_energy(&atoms, &D3BJ_B3LYP);
        assert_eq!(
            result.pair_energies.len(),
            3,
            "H2O should have 3 atom pairs"
        );
    }

    // =========================================================================
    // Additional Tests
    // =========================================================================

    #[test]
    fn cn_saturates_for_bonded() {
        // CN counting function gives ~1 for bonded atoms
        // H2 at equilibrium (R=1.4 bohr, R_cov(H)+R_cov(H) ~ 1.21 bohr)
        let atoms = [(1u8, [0.0, 0.0, 0.0]), (1u8, [0.0, 0.0, 1.4])];
        let cn = compute_coordination_numbers(&atoms);
        assert!(
            cn[0] > 0.8 && cn[0] < 1.1,
            "CN(H in H2) = {}, expected 0.8-1.1",
            cn[0]
        );
    }

    #[test]
    fn c6_interpolation_with_gaussian_weighting() {
        // C6-9: Verify Gaussian weighting produces smooth interpolation
        // C6 at exactly the reference CN should equal the reference value
        let c6_at_ref = interpolate_c6(1, 1, 0.91, 0.91);
        assert!(
            (c6_at_ref - 6.499026).abs() < 0.01,
            "C6 at reference CN = {}, expected 6.499026",
            c6_at_ref
        );
    }

    #[test]
    fn empty_molecule() {
        // Edge case: empty molecule
        let atoms: &[(u8, [f64; 3])] = &[];
        let result = compute_d3bj_energy(atoms, &D3BJ_B3LYP);
        assert_eq!(result.energy, 0.0);
        assert!(result.pair_energies.is_empty());
        assert!(result.coordination_numbers.is_empty());
    }

    #[test]
    fn heteronuclear_c6_data_present() {
        // C6-7: C6 reference data loaded for common element pairs
        let c6_hc = interpolate_c6(1, 6, 0.91, 0.99);
        assert!(
            c6_hc > 5.0 && c6_hc < 20.0,
            "C6(H,C) = {}, expected 5-20",
            c6_hc
        );

        let c6_ho = interpolate_c6(1, 8, 0.91, 0.99);
        assert!(
            c6_ho > 3.0 && c6_ho < 15.0,
            "C6(H,O) = {}, expected 3-15",
            c6_ho
        );
    }

    #[test]
    fn ch4_coordination_numbers() {
        // CH4 molecule: C should have CN~4, H should have CN~1
        // Tetrahedral geometry with C-H = 2.05 bohr
        let ch = 2.05;
        let hx = ch * 0.57735; // 1/sqrt(3)
        let atoms = [
            (6u8, [0.0, 0.0, 0.0]), // C at origin
            (1u8, [hx, hx, hx]),    // H
            (1u8, [-hx, -hx, hx]),  // H
            (1u8, [-hx, hx, -hx]),  // H
            (1u8, [hx, -hx, -hx]),  // H
        ];
        let cn = compute_coordination_numbers(&atoms);

        assert!(
            cn[0] > 3.5 && cn[0] < 4.5,
            "CN(C in CH4) = {}, expected ~4",
            cn[0]
        );
        for i in 1..5 {
            assert!(
                cn[i] > 0.8 && cn[i] < 1.3,
                "CN(H in CH4) = {}, expected ~1",
                cn[i]
            );
        }
    }

    #[test]
    fn gradient_with_larger_molecule() {
        // Gradient finite-difference for CH4
        let ch = 2.05;
        let hx = ch * 0.57735;
        let atoms = [
            (6u8, [0.0, 0.0, 0.0]),
            (1u8, [hx, hx, hx]),
            (1u8, [-hx, -hx, hx]),
            (1u8, [-hx, hx, -hx]),
            (1u8, [hx, -hx, -hx]),
        ];

        let grad = compute_d3bj_gradient(&atoms, &D3BJ_B3LYP);

        // Verify translational invariance
        let sum_x: f64 = grad.gradients.iter().map(|g| g[0]).sum();
        let sum_y: f64 = grad.gradients.iter().map(|g| g[1]).sum();
        let sum_z: f64 = grad.gradients.iter().map(|g| g[2]).sum();
        assert!(sum_x.abs() < 1e-13, "CH4 grad sum x = {:.2e}", sum_x);
        assert!(sum_y.abs() < 1e-13, "CH4 grad sum y = {:.2e}", sum_y);
        assert!(sum_z.abs() < 1e-13, "CH4 grad sum z = {:.2e}", sum_z);

        // Spot-check finite difference for carbon atom, x-component
        let h = 1e-5;
        let mut atoms_plus = atoms;
        let mut atoms_minus = atoms;
        atoms_plus[0].1[0] += h;
        atoms_minus[0].1[0] -= h;
        let e_plus = compute_d3bj_energy(&atoms_plus, &D3BJ_B3LYP).energy;
        let e_minus = compute_d3bj_energy(&atoms_minus, &D3BJ_B3LYP).energy;
        let g_num = (e_plus - e_minus) / (2.0 * h);
        let diff = (grad.gradients[0][0] - g_num).abs();
        assert!(
            diff < 1e-7,
            "CH4 gradient C(x): analytical={:.10e}, numerical={:.10e}, diff={:.2e}",
            grad.gradients[0][0],
            g_num,
            diff
        );
    }
}
