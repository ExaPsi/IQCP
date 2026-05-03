//! Harmonic normal mode analysis.
//!
//! Given a molecular geometry and an analytical Hessian in Hartree/bohr²,
//! this module computes vibrational frequencies, normal modes, reduced
//! masses, force constants (in both atomic units and mDyne/Å), rotor type
//! classification, and rotational constants — all consistent with PySCF's
//! `pyscf/hessian/thermo.py::harmonic_analysis`.
//!
//! # Algorithm outline
//!
//! 1. **Mass-weighted Hessian.** For each atom pair (A, B) and Cartesian
//!    directions (d, e):
//!    $$\tilde H_{3A+d,\,3B+e} = H_{3A+d,\,3B+e} / \sqrt{m_A m_B}.$$
//!    Masses are in atomic mass units (amu) from
//!    [`crate::constants::ATOMIC_MASSES`]. Matches PySCF `thermo.py:55`.
//!
//! 2. **Center of mass + principal axis frame.** Compute the COM and the
//!    moment-of-inertia tensor, then diagonalize to get principal moments
//!    (ascending, then reversed to match PySCF `thermo.py:246-249`).
//!
//! 3. **Rotor classification.** From the rotational constants in GHz,
//!    classify as one of
//!    [`RotorType::Atom`], [`RotorType::Linear`],
//!    [`RotorType::SphericalTop`], [`RotorType::SymmetricTop`],
//!    [`RotorType::AsymmetricTop`].  Uses PySCF's thresholds
//!    (`thermo.py:263-270`) and adds internal refinement of the generic
//!    `REGULAR` category into spherical/symmetric/asymmetric for pedagogy.
//!
//! 4. **Translation + rotation basis.** Build 6 (or 5 for linear, 3 for
//!    single atom) basis vectors that span the translational and
//!    rotational modes in the mass-weighted space. The three translational
//!    vectors are `T_d[3A+d] = √m_A`; the rotational vectors come from the
//!    cross product of the COM-shifted coordinate (in principal-axis frame)
//!    with each principal axis — see `build_tr_basis` below. Matches PySCF
//!    `thermo.py:233-260`.
//!
//! 5. **Projection.** Use QR decomposition on the T/R basis transposed
//!    (shape `3N × N_TR`), form the projection operator
//!    $P = I - Q Q^\top$, then diagonalize $P$ and keep the eigenvectors
//!    with eigenvalue $> 10^{-7}$. These form an orthonormal basis `bvec`
//!    for the internal subspace. Project the mass-weighted Hessian as
//!    $H_\text{proj} = \mathrm{bvec}^\top \tilde H\, \mathrm{bvec}$.
//!    Matches PySCF `thermo.py:73-79`.
//!
//! 6. **Eigendecomposition.** Diagonalize the projected Hessian with
//!    [`nalgebra::SymmetricEigen`], sort eigenvalues ascending, and
//!    back-transform: `mode_3N = bvec @ mode_internal`.
//!
//! 7. **Frequency conversion.** For each eigenvalue $\lambda$:
//!    - if $\lambda \ge 0$: $\nu = \sqrt{\lambda}$ (real vibration),
//!    - if $\lambda < 0$: $\nu = -\sqrt{-\lambda}$ (imaginary, stored
//!      as negative real per Gaussian/PySCF convention).
//!
//!    Then convert to cm⁻¹ using [`crate::constants::AU_TO_CM1`].
//!
//! 8. **Derived quantities.** Reduced mass (amu), force constant
//!    (mDyne/Å), vibrational temperature (K), Cartesian and mass-weighted
//!    normal modes, rotational constants (GHz), and principal moments
//!    (amu·bohr²).
//!
//! # References
//!
//! - Wilson, Decius & Cross (1955), *Molecular Vibrations*. McGraw-Hill.
//! - Herzberg (1945), *Molecular Spectra and Molecular Structure, Vol. II*.
//! - Gaussian vibrational analysis manual:
//!   <https://gaussian.com/vib/>
//! - PySCF `pyscf/hessian/thermo.py:40-270` (function-level traceability).
//! - Dupuis, Rys & King (1976) — context (not used directly here).

use nalgebra::{DMatrix, Matrix3, SymmetricEigen, Vector3};
use thiserror::Error;

use crate::basis::Atom;
use crate::constants::{
    AMU_TO_KG, ATOMIC_MASSES, AU_TO_CM1, AU_TO_HZ, BOHR_SI, BOLTZMANN, HARTREE_TO_JOULE, HBAR,
    PLANCK,
};

/// Module version (matches crate version).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Linear-dependency threshold for the T/R projection step.
///
/// Matches PySCF `pyscf/hessian/thermo.py:37` (`LINDEP_THRESHOLD = 1e-7`).
const LINDEP_THRESHOLD: f64 = 1e-7;

// ============================================================================
// Public types
// ============================================================================

/// Classification of the molecular rotor type based on principal moments
/// of inertia / rotational constants.
///
/// Matches PySCF `_get_rotor_type` (`thermo.py:263-270`) with IQCP
/// refinement of the generic `REGULAR` category into sub-types for
/// pedagogical display.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RotorType {
    /// Single atom: all three principal moments are ~zero, so all three
    /// rotational constants are infinite.
    Atom,
    /// Linear molecule: one zero principal moment, two equal non-zero
    /// moments. PySCF test: `rot_const[0] > 1e8 && abs(rot_const[1] -
    /// rot_const[2]) < 1e-3`.
    Linear,
    /// Spherical top: all three rotational constants equal
    /// (e.g., CH₄, SF₆).
    SphericalTop,
    /// Symmetric top: two rotational constants equal, one different
    /// (e.g., NH₃, CH₃F). Prolate vs oblate is encoded in the relative
    /// ordering of the moments.
    SymmetricTop,
    /// Asymmetric top: all three rotational constants distinct
    /// (e.g., H₂O, NO₂).
    AsymmetricTop,
}

/// Errors returned by [`harmonic_analysis`] and its helpers.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum ThermoError {
    /// The Hessian is not square.
    #[error("Hessian is not square: rows={rows}, cols={cols}")]
    NonSquareHessian { rows: usize, cols: usize },

    /// The Hessian dimensions do not match `3 * n_atoms`.
    #[error(
        "Hessian size mismatch: expected {expected}x{expected} for {n_atoms} atoms, got {actual}x{actual}"
    )]
    SizeMismatch {
        expected: usize,
        actual: usize,
        n_atoms: usize,
    },

    /// A Hessian element is not finite (NaN or infinity).
    #[error("Non-finite Hessian element at ({row}, {col}): {value}")]
    NonFiniteHessian { row: usize, col: usize, value: f64 },

    /// Atomic number not supported by IQCP (must be 1–18).
    #[error("Unsupported atomic number: {0} (must be 1-18)")]
    UnsupportedElement(u8),
}

/// Full output of a normal mode analysis: frequencies, modes, and all
/// derived quantities.
///
/// # Layout notes
///
/// - `freq_wavenumber[i]` is the frequency of mode `i` in cm⁻¹. Imaginary
///   frequencies (from a saddle point or non-stationary geometry) are
///   stored as negative real values, matching PySCF's `imaginary_freq=False`
///   convention.
/// - `norm_mode[i][A][d]` is the Cartesian displacement of atom `A` in
///   direction `d` (0=x, 1=y, 2=z) for mode `i`. Modes are normalized
///   such that the mass-weighted mode vector has unit length
///   (`norm_mode_mw`).
/// - `norm_mode_mw[i][A][d]` is the mass-weighted form of the same mode;
///   downstream US-097/098 use this to project derivatives.
/// - `rotational_constants_ghz` are sorted from largest to smallest
///   (A ≥ B ≥ C), which corresponds to principal moments sorted
///   ascending (smallest moment ↔ largest constant).
/// - `principal_moments_amu_bohr2` are sorted ascending in amu·bohr².
/// - `vib_temperature[i] = freq_au[i] * AU_TO_HZ * h / k_B` in Kelvin.
#[derive(Debug, Clone)]
pub struct FrequencyInfo {
    /// Vibrational frequencies in cm⁻¹ (negative = imaginary).
    pub freq_wavenumber: Vec<f64>,
    /// Vibrational frequencies in atomic units (sign matches `freq_wavenumber`).
    pub freq_au: Vec<f64>,
    /// Reduced masses in amu.
    pub reduced_mass: Vec<f64>,
    /// Force constant eigenvalues in Ha/bohr² (sorted ascending).
    pub force_const_au: Vec<f64>,
    /// Force constants in mDyne/Å.
    pub force_const_dyne: Vec<f64>,
    /// Cartesian normal modes: `[mode][atom][xyz]`.
    pub norm_mode: Vec<Vec<[f64; 3]>>,
    /// Mass-weighted normal modes: `[mode][atom][xyz]`.
    pub norm_mode_mw: Vec<Vec<[f64; 3]>>,
    /// Rotor type classification.
    pub rot_type: RotorType,
    /// Rotational constants in GHz (A ≥ B ≥ C; ∞ for zero moments).
    pub rotational_constants_ghz: [f64; 3],
    /// Principal moments of inertia in amu·bohr² (sorted ascending).
    pub principal_moments_amu_bohr2: [f64; 3],
    /// Vibrational temperatures in Kelvin.
    pub vib_temperature: Vec<f64>,
    /// Number of vibrational modes (3N-6 / 3N-5 / 0).
    pub n_modes: usize,
    /// Number of atoms.
    pub n_atoms: usize,
}

// ============================================================================
// Helpers: mass lookup + Hessian validation
// ============================================================================

/// Look up the mass of an atom in amu, mapping IQCP element errors to
/// [`ThermoError::UnsupportedElement`].
#[inline]
fn mass_of(atom: &Atom) -> Result<f64, ThermoError> {
    let z = atom.atomic_number;
    if z == 0 || (z as usize) >= ATOMIC_MASSES.len() {
        return Err(ThermoError::UnsupportedElement(z));
    }
    Ok(ATOMIC_MASSES[z as usize])
}

/// Build a `Vec<f64>` of atomic masses in amu, one per atom.
fn masses_of(atoms: &[Atom]) -> Result<Vec<f64>, ThermoError> {
    atoms.iter().map(mass_of).collect()
}

/// Validate that the Hessian has shape `(3N, 3N)` with all finite entries.
fn validate_hessian(hessian: &DMatrix<f64>, n_atoms: usize) -> Result<(), ThermoError> {
    let rows = hessian.nrows();
    let cols = hessian.ncols();
    if rows != cols {
        return Err(ThermoError::NonSquareHessian { rows, cols });
    }
    let expected = 3 * n_atoms;
    if rows != expected {
        return Err(ThermoError::SizeMismatch {
            expected,
            actual: rows,
            n_atoms,
        });
    }
    for i in 0..rows {
        for j in 0..cols {
            let v = hessian[(i, j)];
            if !v.is_finite() {
                return Err(ThermoError::NonFiniteHessian {
                    row: i,
                    col: j,
                    value: v,
                });
            }
        }
    }
    Ok(())
}

// ============================================================================
// Phase 1: Mass-weighted Hessian
// ============================================================================

/// Construct the mass-weighted Hessian.
///
/// For each atom pair (A, B) and Cartesian directions (d, e):
/// ```text
/// H_mw[3A+d, 3B+e] = H[3A+d, 3B+e] / sqrt(m_A * m_B)
/// ```
/// Masses are looked up by atomic number from [`ATOMIC_MASSES`] (amu).
///
/// # Reference
/// - PySCF `pyscf/hessian/thermo.py:55`:
///   `mass_hess = numpy.einsum('pqxy,p,q->pqxy', hess, mass**-.5, mass**-.5)`
pub fn build_mass_weighted_hessian(
    atoms: &[Atom],
    hessian: &DMatrix<f64>,
) -> Result<DMatrix<f64>, ThermoError> {
    validate_hessian(hessian, atoms.len())?;
    let masses = masses_of(atoms)?;
    let n3 = 3 * atoms.len();
    let mut out = DMatrix::zeros(n3, n3);
    for (a, &m_a) in masses.iter().enumerate() {
        for (b, &m_b) in masses.iter().enumerate() {
            let scale = 1.0 / (m_a * m_b).sqrt();
            for d in 0..3 {
                for e in 0..3 {
                    let i = 3 * a + d;
                    let j = 3 * b + e;
                    out[(i, j)] = hessian[(i, j)] * scale;
                }
            }
        }
    }
    Ok(out)
}

// ============================================================================
// Phase 2: Center of mass, inertia tensor, principal axes, rotor type
// ============================================================================

/// Compute the center of mass of a molecule in bohr.
///
/// `R_com = (Σ m_A * R_A) / (Σ m_A)` with masses in amu from
/// [`ATOMIC_MASSES`].
///
/// # Reference
/// - PySCF `pyscf/hessian/thermo.py:51`:
///   `mass_center = einsum('z,zx->x', mass, atom_coords) / mass.sum()`
pub fn center_of_mass(atoms: &[Atom]) -> [f64; 3] {
    let mut num = [0.0f64; 3];
    let mut den = 0.0f64;
    for atom in atoms {
        // Safe: mass_of errors are propagated in public API; here we only
        // use COM internally, and the public harmonic_analysis validates
        // atoms up-front. Fall back to 0 to avoid panics.
        let m = ATOMIC_MASSES
            .get(atom.atomic_number as usize)
            .copied()
            .unwrap_or(0.0);
        num[0] += m * atom.position[0];
        num[1] += m * atom.position[1];
        num[2] += m * atom.position[2];
        den += m;
    }
    if den == 0.0 {
        return [0.0, 0.0, 0.0];
    }
    [num[0] / den, num[1] / den, num[2] / den]
}

/// Compute the moment-of-inertia tensor for COM-shifted coordinates
/// (in bohr) and masses (in amu).
///
/// ```text
/// I_{de} = δ_{de} Σ_A m_A |r_A|² - Σ_A m_A r_{A,d} r_{A,e}
/// ```
///
/// Returns a symmetric 3×3 matrix in amu·bohr².
///
/// # Reference
/// - PySCF `pyscf/hessian/thermo.py:119-120, 244-245`:
///   `im = einsum('m,mx,my->xy', mass, coords, coords);
///    im = eye(3)*im.trace() - im`
pub fn moment_of_inertia_tensor(coords_com: &[[f64; 3]], masses: &[f64]) -> Matrix3<f64> {
    assert_eq!(coords_com.len(), masses.len());
    let mut im = Matrix3::<f64>::zeros();
    for (pos, &m) in coords_com.iter().zip(masses.iter()) {
        for d in 0..3 {
            for e in 0..3 {
                im[(d, e)] += m * pos[d] * pos[e];
            }
        }
    }
    // I = tr(M) * I3 - M  with M = Σ m r_d r_e
    let tr = im.trace();
    let mut out = Matrix3::<f64>::identity() * tr;
    out -= im;
    out
}

/// Diagonalize the moment-of-inertia tensor, returning principal moments
/// (ascending, in amu·bohr²) and the principal axes as a 3×3 matrix
/// whose columns are the axes.
///
/// # PySCF compatibility
///
/// PySCF reverses the eigenvector ordering after diagonalization
/// (`thermo.py:246-249`):
/// ```python
/// w, paxes = eigh(im)
/// w = w[::-1]
/// paxes = paxes[:, ::-1]
/// ```
/// so that the **smallest moment ends up last**. This function matches
/// PySCF exactly: the returned `moments` array is sorted ascending (i.e.,
/// reversed from the PySCF-internal descending order), and the `paxes`
/// column order matches PySCF's `paxes[:, ::-1]` — i.e., `paxes[:, 0]`
/// corresponds to the **largest** moment and `paxes[:, 2]` to the
/// **smallest**.
///
/// The T/R basis builder [`build_tr_basis`] uses this convention when
/// forming `coords_in_rot_frame = coords @ paxes`.
pub fn principal_axes(inertia: Matrix3<f64>) -> ([f64; 3], Matrix3<f64>) {
    // nalgebra SymmetricEigen for fixed-size matrices returns eigenvalues
    // and eigenvectors. Order is not guaranteed, so we sort ascending then
    // reverse (to match PySCF's `w[::-1], paxes[:,::-1]`).
    let eig = inertia.symmetric_eigen();
    let mut pairs: Vec<(f64, Vector3<f64>)> = (0..3)
        .map(|i| (eig.eigenvalues[i], eig.eigenvectors.column(i).into_owned()))
        .collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let moments_asc = [pairs[0].0, pairs[1].0, pairs[2].0];

    // PySCF reversal: reversed so column 0 = largest moment, column 2 = smallest
    let mut paxes = Matrix3::<f64>::zeros();
    paxes.set_column(0, &pairs[2].1);
    paxes.set_column(1, &pairs[1].1);
    paxes.set_column(2, &pairs[0].1);

    (moments_asc, paxes)
}

/// Compute rotational constants in GHz from principal moments
/// (amu·bohr², sorted ascending).
///
/// Uses the PySCF formula (`thermo.py:123-127`):
/// ```text
/// unit_im = ATOMIC_MASS * BOHR_SI²        [kg·m²]
/// unit_hz = ℏ / (4π * unit_im)            [Hz]
/// B_i [GHz] = unit_hz / I_i * 1e-9
/// ```
///
/// Zero moments (e.g., for a single atom or the unique axis of a linear
/// molecule) produce `f64::INFINITY` for the corresponding constant.
///
/// # Ordering
///
/// Because moments are sorted ascending, the returned `[B_A, B_B, B_C]`
/// are sorted **descending**: `B_A ≥ B_B ≥ B_C` (with ∞ ≥ ∞).
pub fn rotational_constants_ghz(moments_amu_bohr2: [f64; 3]) -> [f64; 3] {
    let unit_im = AMU_TO_KG * BOHR_SI * BOHR_SI;
    let unit_hz = HBAR / (4.0 * std::f64::consts::PI * unit_im);
    let mut out = [f64::INFINITY; 3];
    // moments are ascending: smallest moment → largest B
    // PySCF rot_const sorts the eigenvalues ascending then does unit_hz/e
    // which gives descending rotational constants (largest first).
    for (i, &i_val) in moments_amu_bohr2.iter().enumerate() {
        if i_val <= 0.0 {
            out[i] = f64::INFINITY;
        } else {
            out[i] = unit_hz / i_val * 1e-9;
        }
    }
    out
}

/// Classify the rotor type from rotational constants in GHz.
///
/// # Thresholds (PySCF-compatible, from `thermo.py:263-270`)
///
/// - All three `> 1e8 GHz` → [`RotorType::Atom`]
/// - `rot_const[0] > 1e8 && |rot_const[1] - rot_const[2]| < 1e-3` → [`RotorType::Linear`]
/// - Otherwise PySCF returns `REGULAR`, which IQCP further refines:
///   - All three equal (within 1e-3 relative) → [`RotorType::SphericalTop`]
///   - Two equal, one different → [`RotorType::SymmetricTop`]
///   - All three different → [`RotorType::AsymmetricTop`]
pub fn classify_rotor(rot_const_ghz: [f64; 3]) -> RotorType {
    // Expect rot_const_ghz in PySCF convention: sorted descending
    // (because moments are ascending). rot_const_ghz[0] is the largest.
    let ra = rot_const_ghz[0];
    let rb = rot_const_ghz[1];
    let rc = rot_const_ghz[2];

    // Atom: all three infinite (or above 1e8 GHz threshold)
    if ra > 1e8 && rb > 1e8 && rc > 1e8 {
        return RotorType::Atom;
    }

    // Linear: largest is infinite (zero moment along unique axis), other two equal
    // PySCF: rot_const[0] > 1e8 && abs(rot_const[1] - rot_const[2]) < 1e-3
    if ra > 1e8 && (rb - rc).abs() < 1e-3 {
        return RotorType::Linear;
    }

    // Regular sub-classification. Use relative tolerance for comparing
    // finite rotational constants.
    let rel_tol = 1e-3;
    let approx_eq = |x: f64, y: f64| {
        let denom = x.abs().max(y.abs()).max(1e-30);
        (x - y).abs() / denom < rel_tol
    };

    let ab_eq = approx_eq(ra, rb);
    let bc_eq = approx_eq(rb, rc);
    let ac_eq = approx_eq(ra, rc);

    if ab_eq && bc_eq && ac_eq {
        RotorType::SphericalTop
    } else if ab_eq || bc_eq {
        // Two rotational constants equal, one different → symmetric top
        RotorType::SymmetricTop
    } else {
        RotorType::AsymmetricTop
    }
}

// ============================================================================
// Phase 3: T/R basis construction and projection
// ============================================================================

/// Build the translation + rotation basis vectors for the T/R projection.
///
/// Returns an `(N_TR, 3N)` matrix where each **row** is a T/R basis
/// vector of length `3N`. `N_TR` is:
///
/// - `3` for a single atom (translations only)
/// - `5` for a linear molecule (3 translations + 2 rotations)
/// - `6` otherwise (3 translations + 3 rotations)
///
/// # Translations
///
/// ```text
/// T_d[3A+d] = √m_A  (else 0)       for d ∈ {x, y, z}
/// ```
///
/// # Rotations (in the principal axis frame)
///
/// ```text
/// coords_in_rot_frame = coords_com @ paxes         # (natm, 3) @ (3, 3)
/// (cx, cy, cz) = coords_in_rot_frame[A]
/// Rx[3A:3A+3] = √m_A · (cy · ez - cz · ey)
/// Ry[3A:3A+3] = √m_A · (cz · ex - cx · ez)
/// Rz[3A:3A+3] = √m_A · (cx · ey - cy · ex)
/// ```
/// where `ex, ey, ez` are the three columns of `paxes`.
///
/// For a linear molecule, PySCF takes `TR[3:5]` = (Rx, Ry) — the rotation
/// about the unique axis has zero mass-weighted norm and is filtered out.
///
/// # Reference
///
/// PySCF `pyscf/hessian/thermo.py:233-260` (`_get_TR`).
pub fn build_tr_basis(
    masses: &[f64],
    coords_com: &[[f64; 3]],
    paxes: &Matrix3<f64>,
    rotor: RotorType,
) -> DMatrix<f64> {
    let natm = masses.len();
    let n3 = 3 * natm;
    let massp: Vec<f64> = masses.iter().map(|m| m.sqrt()).collect();

    // Translation rows are always present (3 of them)
    // Layout: row d, columns [3A+d] = √m_A
    let n_tr_rows = match rotor {
        RotorType::Atom => 3,
        RotorType::Linear => 5,
        _ => 6,
    };

    let mut tr = DMatrix::<f64>::zeros(n_tr_rows, n3);

    // T_x, T_y, T_z rows
    for a in 0..natm {
        for d in 0..3 {
            tr[(d, 3 * a + d)] = massp[a];
        }
    }

    if matches!(rotor, RotorType::Atom) {
        return tr;
    }

    // Principal axes (columns of paxes)
    let ex = paxes.column(0).into_owned();
    let ey = paxes.column(1).into_owned();
    let ez = paxes.column(2).into_owned();

    // Compute (cx, cy, cz) per atom = coords_com[a] dot (ex, ey, ez)
    // (this is coords_in_rot_frame = coords_com @ paxes)
    // Build all three rotation rows; we'll select the right subset for Linear.
    let mut rx = vec![0.0f64; n3];
    let mut ry = vec![0.0f64; n3];
    let mut rz = vec![0.0f64; n3];
    for a in 0..natm {
        let r = Vector3::new(coords_com[a][0], coords_com[a][1], coords_com[a][2]);
        let cx = r.dot(&ex);
        let cy = r.dot(&ey);
        let cz = r.dot(&ez);
        // Rx vector for atom a: √m_A * (cy * ez - cz * ey), a 3-vector
        let rxa = massp[a] * (ez * cy - ey * cz);
        let rya = massp[a] * (ex * cz - ez * cx);
        let rza = massp[a] * (ey * cx - ex * cy);
        for d in 0..3 {
            rx[3 * a + d] = rxa[d];
            ry[3 * a + d] = rya[d];
            rz[3 * a + d] = rza[d];
        }
    }

    match rotor {
        RotorType::Linear => {
            // PySCF: TR[3:5] = Rx, Ry
            for k in 0..n3 {
                tr[(3, k)] = rx[k];
                tr[(4, k)] = ry[k];
            }
        }
        _ => {
            // Regular: all 6 rows
            for k in 0..n3 {
                tr[(3, k)] = rx[k];
                tr[(4, k)] = ry[k];
                tr[(5, k)] = rz[k];
            }
        }
    }

    tr
}

/// Project out the T/R basis from the mass-weighted Hessian.
///
/// Returns `(h_proj, bvec)` where:
/// - `bvec` is an orthonormal `3N × n_internal` matrix whose columns
///   span the internal (vibrational) subspace.
/// - `h_proj = bvec^T · mass_weighted_hessian · bvec` is the projected
///   Hessian in the internal subspace, shape `(n_internal, n_internal)`.
///
/// # Algorithm (matches PySCF `thermo.py:73-79`)
///
/// 1. QR decomposition of `tr_basis.transpose()` (shape `3N × N_TR`).
///    The `Q` factor has orthonormal columns spanning the T/R subspace.
/// 2. Projector `P = I - Q Q^T` (shape `3N × 3N`).
/// 3. Eigendecompose `P`: eigenvectors with eigenvalue `> LINDEP_THRESHOLD`
///    form an orthonormal basis of the internal subspace (the eigenvalues
///    of `P` are 0 on the T/R subspace and 1 on its orthogonal
///    complement).
/// 4. Collect those eigenvectors as columns of `bvec`.
/// 5. Compute `h_proj = bvec^T @ H_mw @ bvec`.
///
/// # PySCF correspondence
///
/// ```python
/// q, r = numpy.linalg.qr(TRspace.T)            # thermo.py:75
/// P = numpy.eye(n3) - q.dot(q.T)                # thermo.py:76
/// w, v = numpy.linalg.eigh(P)                   # thermo.py:77
/// bvec = v[:, w > LINDEP_THRESHOLD]             # thermo.py:78
/// h = bvec.T @ h_mw @ bvec                      # thermo.py:79
/// ```
pub fn project_out_tr(
    mass_weighted_hessian: &DMatrix<f64>,
    tr_basis: &DMatrix<f64>,
) -> (DMatrix<f64>, DMatrix<f64>) {
    let n3 = mass_weighted_hessian.nrows();
    let n_tr = tr_basis.nrows();

    // QR of TRspace.T (shape 3N × N_TR)
    let tr_t = tr_basis.transpose();
    let qr = tr_t.qr();
    let q_full = qr.q(); // shape 3N × min(3N, N_TR) = 3N × N_TR
                         // Keep only the first n_tr columns (those are the orthonormal basis
                         // of the column space). nalgebra's .q() returns a matrix of size
                         // (3N, min(3N, N_TR)).
    let k = q_full.ncols().min(n_tr);
    let q = q_full.columns(0, k).into_owned();

    // Projection operator P = I - Q Q^T  (3N × 3N)
    let qqt = &q * q.transpose();
    let identity = DMatrix::<f64>::identity(n3, n3);
    let p = &identity - qqt;

    // Eigendecompose P; its eigenvalues are 0 on the TR subspace and
    // 1 on the internal subspace. We keep eigenvectors where eigenvalue
    // > LINDEP_THRESHOLD.
    let eig = SymmetricEigen::new(p);
    let mut keep_indices = Vec::new();
    for i in 0..eig.eigenvalues.len() {
        if eig.eigenvalues[i] > LINDEP_THRESHOLD {
            keep_indices.push(i);
        }
    }

    let n_internal = keep_indices.len();
    let mut bvec = DMatrix::<f64>::zeros(n3, n_internal);
    for (new_col, &old_col) in keep_indices.iter().enumerate() {
        for row in 0..n3 {
            bvec[(row, new_col)] = eig.eigenvectors[(row, old_col)];
        }
    }

    // h_proj = bvec^T @ H_mw @ bvec
    let h_proj = bvec.transpose() * mass_weighted_hessian * &bvec;

    (h_proj, bvec)
}

// ============================================================================
// Phase 4: Frequency conversion + derived quantities
// ============================================================================

/// Convert a set of force-constant eigenvalues (atomic units of
/// mass-weighted Hessian) to vibrational frequencies in atomic units and
/// cm⁻¹.
///
/// For each eigenvalue $\lambda$:
/// - $\lambda \ge 0$: real frequency $\nu = \sqrt\lambda$
/// - $\lambda < 0$: imaginary frequency stored as $-\sqrt{-\lambda}$
///   (Gaussian / PySCF `imaginary_freq=False` convention)
///
/// Returns `(freq_au, freq_cm1)`.
///
/// # Reference
/// - PySCF `pyscf/hessian/thermo.py:85-93`
pub fn force_const_to_freq(force_const_au: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut freq_au = Vec::with_capacity(force_const_au.len());
    let mut freq_cm1 = Vec::with_capacity(force_const_au.len());
    for &lambda in force_const_au {
        let fau = if lambda >= 0.0 {
            lambda.sqrt()
        } else {
            -((-lambda).sqrt())
        };
        freq_au.push(fau);
        freq_cm1.push(fau * AU_TO_CM1);
    }
    (freq_au, freq_cm1)
}

// ============================================================================
// Phase 5: Top-level harmonic_analysis
// ============================================================================

/// Perform harmonic normal mode analysis.
///
/// Given a molecular geometry (`atoms`, positions in bohr) and a 3N×3N
/// Hessian in Ha/bohr² with layout `hessian[3A+d, 3B+e]`, compute
/// vibrational frequencies, normal modes, reduced masses, force
/// constants, rotor type, and rotational constants.
///
/// This is the main entry point for normal mode analysis in IQCP.
///
/// # Errors
/// - [`ThermoError::NonSquareHessian`] if the Hessian is not square
/// - [`ThermoError::SizeMismatch`] if Hessian dimensions ≠ `3 * n_atoms`
/// - [`ThermoError::NonFiniteHessian`] if any element is NaN or inf
/// - [`ThermoError::UnsupportedElement`] if any atomic number is not in 1..=18
///
/// # Reference
/// - PySCF `pyscf/hessian/thermo.py:40-109` (`harmonic_analysis`)
pub fn harmonic_analysis(
    atoms: &[Atom],
    hessian: &DMatrix<f64>,
) -> Result<FrequencyInfo, ThermoError> {
    let natm = atoms.len();

    // Validate Hessian shape and values
    validate_hessian(hessian, natm)?;
    let masses = masses_of(atoms)?;

    // ------------------------------------------------------------------
    // Edge case: single atom — no vibrations, no rotations.
    // ------------------------------------------------------------------
    if natm == 1 {
        return Ok(FrequencyInfo {
            freq_wavenumber: Vec::new(),
            freq_au: Vec::new(),
            reduced_mass: Vec::new(),
            force_const_au: Vec::new(),
            force_const_dyne: Vec::new(),
            norm_mode: Vec::new(),
            norm_mode_mw: Vec::new(),
            rot_type: RotorType::Atom,
            rotational_constants_ghz: [f64::INFINITY; 3],
            principal_moments_amu_bohr2: [0.0; 3],
            vib_temperature: Vec::new(),
            n_modes: 0,
            n_atoms: 1,
        });
    }

    // ------------------------------------------------------------------
    // Step 1: Mass-weighted Hessian (PySCF thermo.py:55)
    // ------------------------------------------------------------------
    let h_mw = build_mass_weighted_hessian(atoms, hessian)?;

    // ------------------------------------------------------------------
    // Step 2: Center of mass, COM-shifted coords, inertia tensor,
    // principal axes (PySCF thermo.py:51-53, 244-249)
    // ------------------------------------------------------------------
    let com = center_of_mass(atoms);
    let coords_com: Vec<[f64; 3]> = atoms
        .iter()
        .map(|a| {
            [
                a.position[0] - com[0],
                a.position[1] - com[1],
                a.position[2] - com[2],
            ]
        })
        .collect();
    let inertia = moment_of_inertia_tensor(&coords_com, &masses);
    let (moments, paxes) = principal_axes(inertia);

    // ------------------------------------------------------------------
    // Step 3: Rotational constants + rotor classification
    // (PySCF thermo.py:111-132, 263-270)
    // ------------------------------------------------------------------
    let rot_const_ghz = rotational_constants_ghz(moments);
    let rot_type = classify_rotor(rot_const_ghz);

    // ------------------------------------------------------------------
    // Step 4: Build T/R basis in principal axis frame and project it out
    // of the mass-weighted Hessian (PySCF thermo.py:58-79)
    // ------------------------------------------------------------------
    let tr_basis = build_tr_basis(&masses, &coords_com, &paxes, rot_type);
    let (h_proj, bvec) = project_out_tr(&h_mw, &tr_basis);

    // ------------------------------------------------------------------
    // Step 5: Eigendecompose the projected Hessian and back-transform
    // modes to 3N space (PySCF thermo.py:80-81)
    // ------------------------------------------------------------------
    let (force_const_au, mode_3n) = if h_proj.nrows() == 0 {
        (Vec::<f64>::new(), DMatrix::<f64>::zeros(3 * natm, 0))
    } else {
        let eig = SymmetricEigen::new(h_proj);
        let m = eig.eigenvalues.len();
        // Sort ascending
        let mut indices: Vec<usize> = (0..m).collect();
        indices.sort_by(|&a, &b| eig.eigenvalues[a].partial_cmp(&eig.eigenvalues[b]).unwrap());
        let fc: Vec<f64> = indices.iter().map(|&i| eig.eigenvalues[i]).collect();
        // Reorder eigenvectors (columns)
        let n_internal = m;
        let mut mode_internal = DMatrix::<f64>::zeros(n_internal, n_internal);
        for (new_col, &old_col) in indices.iter().enumerate() {
            for row in 0..n_internal {
                mode_internal[(row, new_col)] = eig.eigenvectors[(row, old_col)];
            }
        }
        // mode_3N = bvec @ mode_internal  (shape (3N, n_internal))
        let modes = &bvec * mode_internal;
        (fc, modes)
    };

    // ------------------------------------------------------------------
    // Step 6: Frequency conversion (PySCF thermo.py:85-93)
    // ------------------------------------------------------------------
    let (freq_au, freq_wavenumber) = force_const_to_freq(&force_const_au);

    // ------------------------------------------------------------------
    // Step 7: Cartesian modes + reduced masses + force constants in mDyne/Å
    // + vibrational temperatures (PySCF thermo.py:95-106)
    // ------------------------------------------------------------------
    let n_modes = force_const_au.len();
    let mut norm_mode = Vec::with_capacity(n_modes);
    let mut norm_mode_mw = Vec::with_capacity(n_modes);
    for i in 0..n_modes {
        let mut cart = Vec::with_capacity(natm);
        let mut mw = Vec::with_capacity(natm);
        for a in 0..natm {
            let qx_mw = mode_3n[(3 * a, i)];
            let qy_mw = mode_3n[(3 * a + 1, i)];
            let qz_mw = mode_3n[(3 * a + 2, i)];
            // Cartesian = mass-weighted / √m_A  (PySCF thermo.py:95)
            let inv_sqrt_m = 1.0 / masses[a].sqrt();
            cart.push([qx_mw * inv_sqrt_m, qy_mw * inv_sqrt_m, qz_mw * inv_sqrt_m]);
            mw.push([qx_mw, qy_mw, qz_mw]);
        }
        norm_mode.push(cart);
        norm_mode_mw.push(mw);
    }

    // Reduced mass (amu): 1 / Σ_A |norm_mode[i][A]|²  (PySCF thermo.py:97)
    let reduced_mass: Vec<f64> = norm_mode
        .iter()
        .map(|mode| {
            let s: f64 = mode
                .iter()
                .map(|v| v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
                .sum();
            if s > 0.0 {
                1.0 / s
            } else {
                0.0
            }
        })
        .collect();

    // Force constant in mDyne/Å
    // PySCF `dyne = 1e-2 * HARTREE2J / BOHR_SI**2`  (thermo.py:104)
    //   force_const_dyne = reduced_mass * force_const_au * dyne
    // The numerical value corresponds to 15.569… per Hartree/bohr² per amu,
    // i.e. mDyne/Å when reduced_mass is in amu and force_const_au in Ha/bohr².
    let dyne_factor = 1e-2 * HARTREE_TO_JOULE / (BOHR_SI * BOHR_SI);
    let force_const_dyne: Vec<f64> = reduced_mass
        .iter()
        .zip(force_const_au.iter())
        .map(|(mu, fc)| mu * fc * dyne_factor)
        .collect();

    // Vibrational temperature (K): freq_au * AU_TO_HZ * h / k_B
    // PySCF thermo.py:100-101
    let vib_temperature: Vec<f64> = freq_au
        .iter()
        .map(|fau| fau * AU_TO_HZ * PLANCK / BOLTZMANN)
        .collect();

    Ok(FrequencyInfo {
        freq_wavenumber,
        freq_au,
        reduced_mass,
        force_const_au,
        force_const_dyne,
        norm_mode,
        norm_mode_mw,
        rot_type,
        rotational_constants_ghz: rot_const_ghz,
        principal_moments_amu_bohr2: moments,
        vib_temperature,
        n_modes,
        n_atoms: natm,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn h2o_atoms_bohr() -> Vec<Atom> {
        // H2O bohr coords (arbitrary — used for geometry tests that don't
        // require SCF). Same coords as VP test skeleton.
        vec![
            Atom::new(8, [0.0, 0.0, 0.0]).unwrap(),
            Atom::new(1, [0.0, 1.43, 1.11]).unwrap(),
            Atom::new(1, [0.0, -1.43, 1.11]).unwrap(),
        ]
    }

    // -------------------------------------------------------------
    // Unit tests — mass-weighted Hessian
    // -------------------------------------------------------------

    #[test]
    fn test_build_mass_weighted_hessian_h2o() {
        let atoms = h2o_atoms_bohr();
        let mut h = DMatrix::zeros(9, 9);
        h[(0, 0)] = 1.0;
        h[(1, 1)] = 2.0;
        h[(3, 3)] = 3.0;

        let mwh = build_mass_weighted_hessian(&atoms, &h).unwrap();

        let m_o = ATOMIC_MASSES[8];
        let m_h = ATOMIC_MASSES[1];
        assert_relative_eq!(mwh[(0, 0)], 1.0 / m_o, epsilon = 1e-12);
        assert_relative_eq!(mwh[(1, 1)], 2.0 / m_o, epsilon = 1e-12);
        assert_relative_eq!(mwh[(3, 3)], 3.0 / m_h, epsilon = 1e-12);
    }

    #[test]
    fn test_mass_weighted_hessian_symmetry() {
        let atoms = h2o_atoms_bohr();
        let mut h = DMatrix::<f64>::zeros(9, 9);
        for i in 0..9 {
            for j in i..9 {
                let v = 0.1 * ((i + 1) * (j + 2)) as f64;
                h[(i, j)] = v;
                h[(j, i)] = v;
            }
        }
        let mwh = build_mass_weighted_hessian(&atoms, &h).unwrap();
        let asym = (&mwh - mwh.transpose()).abs().max();
        assert!(asym < 1e-14, "mass-weighted Hessian not symmetric: {asym}");
    }

    #[test]
    fn test_size_mismatch_error() {
        let atoms = h2o_atoms_bohr(); // 3 atoms → expects 9x9
        let h = DMatrix::zeros(6, 6);
        let err = build_mass_weighted_hessian(&atoms, &h).unwrap_err();
        assert!(matches!(err, ThermoError::SizeMismatch { .. }));
    }

    #[test]
    fn test_non_finite_error() {
        let atoms = h2o_atoms_bohr();
        let mut h = DMatrix::<f64>::zeros(9, 9);
        h[(0, 0)] = f64::NAN;
        let err = build_mass_weighted_hessian(&atoms, &h).unwrap_err();
        assert!(matches!(err, ThermoError::NonFiniteHessian { .. }));
    }

    // -------------------------------------------------------------
    // Unit tests — COM, inertia, principal axes
    // -------------------------------------------------------------

    #[test]
    fn test_center_of_mass_h2o() {
        let atoms = h2o_atoms_bohr();
        let com = center_of_mass(&atoms);
        // y components cancel
        assert!(com[1].abs() < 1e-14);
        // O heavy → COM has small z shift (hydrogens pull it toward them)
        assert!(com[2] > 0.0 && com[2] < 0.5);
        // x is zero by geometry
        assert!(com[0].abs() < 1e-14);
    }

    #[test]
    fn test_moment_of_inertia_tensor_symmetry() {
        let atoms = h2o_atoms_bohr();
        let com = center_of_mass(&atoms);
        let coords_com: Vec<[f64; 3]> = atoms
            .iter()
            .map(|a| {
                [
                    a.position[0] - com[0],
                    a.position[1] - com[1],
                    a.position[2] - com[2],
                ]
            })
            .collect();
        let masses: Vec<f64> = atoms
            .iter()
            .map(|a| ATOMIC_MASSES[a.atomic_number as usize])
            .collect();
        let i = moment_of_inertia_tensor(&coords_com, &masses);
        for d in 0..3 {
            for e in 0..3 {
                assert_relative_eq!(i[(d, e)], i[(e, d)], epsilon = 1e-14);
            }
        }
        // All diagonal elements should be positive
        for d in 0..3 {
            assert!(i[(d, d)] > 0.0);
        }
    }

    #[test]
    fn test_principal_axes_reversal_consistency() {
        // Build a trivial diagonal inertia tensor with known moments
        let mut i = Matrix3::zeros();
        i[(0, 0)] = 5.0;
        i[(1, 1)] = 3.0;
        i[(2, 2)] = 7.0;
        let (moments, paxes) = principal_axes(i);
        // Moments must be sorted ascending
        assert!(moments[0] <= moments[1] && moments[1] <= moments[2]);
        assert_relative_eq!(moments[0], 3.0, epsilon = 1e-12);
        assert_relative_eq!(moments[1], 5.0, epsilon = 1e-12);
        assert_relative_eq!(moments[2], 7.0, epsilon = 1e-12);
        // PySCF reversal: paxes col 0 ↔ largest moment = 7 ↔ axis z (2,2)
        // col 2 ↔ smallest moment = 3 ↔ axis y (1,1)
        // paxes columns should be |eigenvector| == 1
        for c in 0..3 {
            let col = paxes.column(c);
            let norm = (col[0] * col[0] + col[1] * col[1] + col[2] * col[2]).sqrt();
            assert_relative_eq!(norm, 1.0, epsilon = 1e-12);
        }
    }

    // -------------------------------------------------------------
    // Unit tests — rotor classification
    // -------------------------------------------------------------

    #[test]
    fn test_classify_rotor_atom() {
        let rot_const = [f64::INFINITY, f64::INFINITY, f64::INFINITY];
        assert_eq!(classify_rotor(rot_const), RotorType::Atom);
    }

    #[test]
    fn test_classify_rotor_linear_co2() {
        // CO2: I_a = 0 → rot_const[0] = inf, rot_const[1] = rot_const[2]
        let rot_const = [
            f64::INFINITY,
            11.737_575_172_773_878,
            11.737_575_172_773_878,
        ];
        assert_eq!(classify_rotor(rot_const), RotorType::Linear);
    }

    #[test]
    fn test_classify_rotor_spherical_top_ch4() {
        let rot_const = [
            158.554_933_354_607_93,
            158.554_933_354_607_93,
            158.554_933_354_607_93,
        ];
        assert_eq!(classify_rotor(rot_const), RotorType::SphericalTop);
    }

    #[test]
    fn test_classify_rotor_asymmetric_top_h2o() {
        let rot_const = [
            820.601_046_730_622_8,
            437.225_477_643_510_77,
            285.244_171_321_706_06,
        ];
        assert_eq!(classify_rotor(rot_const), RotorType::AsymmetricTop);
    }

    #[test]
    fn test_classify_rotor_symmetric_top_nh3_like() {
        // Symmetric top: two rotational constants equal, one different
        // Like NH3 in prolate form: B_A > B_B = B_C
        let rot_const = [298.0, 298.0, 189.0];
        assert_eq!(classify_rotor(rot_const), RotorType::SymmetricTop);
    }

    // -------------------------------------------------------------
    // Unit tests — rotational constants
    // -------------------------------------------------------------

    #[test]
    fn test_rotational_constants_h2o_reasonable_magnitude() {
        // From PySCF golden: ~[820.60, 437.23, 285.24] GHz
        // We don't use the exact hand-H2O here, just check order of magnitude
        // for the actual H2O geometry from the golden reference.
        let atoms = vec![
            Atom::new(8, [0.0, 0.0, 0.221_666_677_2]).unwrap(),
            Atom::new(1, [0.0, 1.431_044_7, -0.886_666_709]).unwrap(),
            Atom::new(1, [0.0, -1.431_044_7, -0.886_666_709]).unwrap(),
        ];
        let masses: Vec<f64> = atoms
            .iter()
            .map(|a| ATOMIC_MASSES[a.atomic_number as usize])
            .collect();
        let com = center_of_mass(&atoms);
        let coords_com: Vec<[f64; 3]> = atoms
            .iter()
            .map(|a| {
                [
                    a.position[0] - com[0],
                    a.position[1] - com[1],
                    a.position[2] - com[2],
                ]
            })
            .collect();
        let i = moment_of_inertia_tensor(&coords_com, &masses);
        let (moments, _) = principal_axes(i);
        let rot = rotational_constants_ghz(moments);
        // A > B > C for asymmetric top
        assert!(rot[0] > rot[1] && rot[1] > rot[2]);
        // H2O rotational constants are in the hundreds of GHz range
        assert!(rot[0] > 100.0 && rot[0] < 2000.0);
        assert!(rot[2] > 50.0 && rot[2] < 1000.0);
    }

    #[test]
    fn test_rotational_constants_zero_moment_infinite() {
        let rot = rotational_constants_ghz([0.0, 1.0, 1.0]);
        assert_eq!(rot[0], f64::INFINITY);
        assert!(rot[1].is_finite());
        assert!(rot[2].is_finite());
    }

    // -------------------------------------------------------------
    // Unit tests — T/R basis
    // -------------------------------------------------------------

    #[test]
    fn test_build_tr_basis_atom_shape() {
        let masses = vec![1.0];
        let coords = vec![[0.0, 0.0, 0.0]];
        let paxes = Matrix3::identity();
        let tr = build_tr_basis(&masses, &coords, &paxes, RotorType::Atom);
        assert_eq!(tr.nrows(), 3);
        assert_eq!(tr.ncols(), 3);
    }

    #[test]
    fn test_build_tr_basis_linear_shape() {
        // 2 atoms, linear → 5 rows, 6 columns
        let masses = vec![1.0, 1.0];
        let coords = vec![[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]];
        let paxes = Matrix3::identity();
        let tr = build_tr_basis(&masses, &coords, &paxes, RotorType::Linear);
        assert_eq!(tr.nrows(), 5);
        assert_eq!(tr.ncols(), 6);
    }

    #[test]
    fn test_build_tr_basis_regular_shape() {
        // 3 atoms, asymmetric → 6 rows, 9 columns
        let masses = vec![16.0, 1.0, 1.0];
        let coords = vec![[0.0, 0.0, 0.0], [0.0, 1.43, 1.11], [0.0, -1.43, 1.11]];
        let paxes = Matrix3::identity();
        let tr = build_tr_basis(&masses, &coords, &paxes, RotorType::AsymmetricTop);
        assert_eq!(tr.nrows(), 6);
        assert_eq!(tr.ncols(), 9);
    }

    #[test]
    fn test_tr_basis_translation_rows() {
        // For any molecule, the first 3 rows are always translations with
        // √m_A scaling. Verify this.
        let masses = vec![16.0, 1.0, 1.0];
        let coords = vec![[0.0, 0.0, 0.0], [0.0, 1.43, 1.11], [0.0, -1.43, 1.11]];
        let paxes = Matrix3::identity();
        let tr = build_tr_basis(&masses, &coords, &paxes, RotorType::AsymmetricTop);
        let sqrt_16 = 4.0;
        let sqrt_1 = 1.0;
        // Row 0 = T_x: [√m_O, 0, 0, √m_H, 0, 0, √m_H, 0, 0]
        assert_relative_eq!(tr[(0, 0)], sqrt_16, epsilon = 1e-12);
        assert_relative_eq!(tr[(0, 3)], sqrt_1, epsilon = 1e-12);
        assert_relative_eq!(tr[(0, 6)], sqrt_1, epsilon = 1e-12);
        assert_relative_eq!(tr[(0, 1)], 0.0, epsilon = 1e-12);
        assert_relative_eq!(tr[(1, 1)], sqrt_16, epsilon = 1e-12);
        assert_relative_eq!(tr[(2, 2)], sqrt_16, epsilon = 1e-12);
    }

    // -------------------------------------------------------------
    // Unit tests — frequency conversion
    // -------------------------------------------------------------

    #[test]
    fn test_force_const_to_freq_positive() {
        let (freq_au, freq_cm1) = force_const_to_freq(&[0.09]);
        assert_relative_eq!(freq_au[0], 0.3, epsilon = 1e-12);
        assert!(freq_cm1[0] > 1500.0 && freq_cm1[0] < 1700.0);
        // 0.3 au * AU_TO_CM1 ≈ 0.3 * 5140.5 ≈ 1542 cm^-1
    }

    #[test]
    fn test_force_const_to_freq_negative_imaginary() {
        let (freq_au, freq_cm1) = force_const_to_freq(&[-0.04]);
        assert!(freq_au[0] < 0.0);
        assert_relative_eq!(freq_au[0], -0.2, epsilon = 1e-12);
        assert!(freq_cm1[0] < 0.0);
    }

    #[test]
    fn test_force_const_to_freq_mixed() {
        let (_, freq_cm1) = force_const_to_freq(&[-0.01, 0.0, 0.01, 0.1]);
        assert!(freq_cm1[0] < 0.0);
        assert_eq!(freq_cm1[1], 0.0);
        assert!(freq_cm1[2] > 0.0);
        assert!(freq_cm1[3] > 0.0);
        // Magnitude check: freq_cm1[2] and freq_cm1[0] should have equal abs
        assert_relative_eq!(freq_cm1[0].abs(), freq_cm1[2].abs(), epsilon = 1e-12);
    }

    // -------------------------------------------------------------
    // Integration test: single atom
    // -------------------------------------------------------------

    #[test]
    fn test_harmonic_analysis_single_atom() {
        let atoms = vec![Atom::new(1, [0.0, 0.0, 0.0]).unwrap()];
        let hessian = DMatrix::zeros(3, 3);
        let result = harmonic_analysis(&atoms, &hessian).unwrap();
        assert_eq!(result.rot_type, RotorType::Atom);
        assert_eq!(result.n_modes, 0);
        assert_eq!(result.n_atoms, 1);
        assert!(result.freq_wavenumber.is_empty());
        assert!(result.norm_mode.is_empty());
        assert!(result.norm_mode_mw.is_empty());
        assert!(result.reduced_mass.is_empty());
    }

    // -------------------------------------------------------------
    // Golden tests: PySCF reference data
    //
    // These tests load the Hessian from PySCF golden JSON files
    // and verify that IQCP's `harmonic_analysis` reproduces the
    // PySCF reference frequencies (within ±2 cm⁻¹), reduced masses,
    // force constants, and rotor classification exactly as PySCF
    // computes them.
    //
    // Generated by: `scripts/phase5/generate_thermo_golden.py`
    // PySCF 2.11.0, RHF/STO-3G, conv_tol=1e-12, cart=True
    // -------------------------------------------------------------

    use serde::Deserialize;

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenAtom {
        #[serde(rename = "Z")]
        z: u8,
        symbol: String,
        pos_bohr: [f64; 3],
    }

    #[derive(Deserialize)]
    #[allow(dead_code)]
    struct GoldenThermo {
        name: String,
        atoms: Vec<GoldenAtom>,
        hessian: Vec<Vec<f64>>,
        energy: f64,
        rot_type: String,
        rotational_constants_ghz: Vec<f64>,
        principal_moments_amu_bohr2: Vec<f64>,
        freq_wavenumber: Vec<f64>,
        freq_au: Vec<f64>,
        reduced_mass: Vec<f64>,
        force_const_dyne: Vec<f64>,
        force_const_au: Vec<f64>,
        norm_mode: Vec<Vec<[f64; 3]>>,
    }

    impl GoldenThermo {
        fn atoms(&self) -> Vec<Atom> {
            self.atoms
                .iter()
                .map(|a| Atom::new(a.z, a.pos_bohr).unwrap())
                .collect()
        }

        fn hessian_matrix(&self) -> DMatrix<f64> {
            let n3 = self.hessian.len();
            let mut h = DMatrix::<f64>::zeros(n3, n3);
            for i in 0..n3 {
                for j in 0..n3 {
                    h[(i, j)] = self.hessian[i][j];
                }
            }
            h
        }
    }

    fn load_golden(name: &str) -> GoldenThermo {
        // include_str! requires a literal path; delegate to per-molecule functions.
        match name {
            "h2o" => {
                let data = include_str!(concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../tests/golden/thermo/h2o_sto3g_rhf.json"
                ));
                serde_json::from_str(data).expect("parse h2o golden")
            }
            "ch4" => {
                let data = include_str!(concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../tests/golden/thermo/ch4_sto3g_rhf.json"
                ));
                serde_json::from_str(data).expect("parse ch4 golden")
            }
            "co2" => {
                let data = include_str!(concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../tests/golden/thermo/co2_sto3g_rhf.json"
                ));
                serde_json::from_str(data).expect("parse co2 golden")
            }
            "h2" => {
                let data = include_str!(concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../tests/golden/thermo/h2_sto3g_rhf.json"
                ));
                serde_json::from_str(data).expect("parse h2 golden")
            }
            _ => panic!("unknown golden name: {}", name),
        }
    }

    /// Compare two sets of sorted frequencies and return the maximum
    /// absolute difference in cm⁻¹.
    fn max_abs_freq_diff(a: &[f64], b: &[f64]) -> f64 {
        assert_eq!(a.len(), b.len());
        let mut asorted: Vec<f64> = a.to_vec();
        let mut bsorted: Vec<f64> = b.to_vec();
        asorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
        bsorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
        asorted
            .iter()
            .zip(bsorted.iter())
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max)
    }

    #[test]
    fn test_harmonic_analysis_h2o_golden_hessian() {
        // AC9: H2O/STO-3G/RHF 3 vibrational frequencies within ±2 cm⁻¹ of PySCF.
        // Uses the PySCF Hessian directly to isolate thermo correctness.
        let golden = load_golden("h2o");
        let atoms = golden.atoms();
        let hessian = golden.hessian_matrix();
        let freq_info = harmonic_analysis(&atoms, &hessian).unwrap();

        assert_eq!(freq_info.n_atoms, 3);
        assert_eq!(freq_info.n_modes, 3);
        assert_eq!(freq_info.rot_type, RotorType::AsymmetricTop);

        // Frequencies within ±2 cm⁻¹ (AC9)
        let mut freq_sorted = freq_info.freq_wavenumber.clone();
        freq_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut ref_sorted = golden.freq_wavenumber.clone();
        ref_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for i in 0..3 {
            let err = (freq_sorted[i] - ref_sorted[i]).abs();
            assert!(
                err < 2.0,
                "H2O mode {}: IQCP {:.4} vs PySCF {:.4}, err {:.4} cm^-1",
                i,
                freq_sorted[i],
                ref_sorted[i],
                err
            );
        }

        // Reduced masses within 0.01 amu. Sort both sets by frequency
        // (ascending) so that mode i corresponds to mode i.
        let iqcp_sorted: Vec<(f64, f64)> = {
            let mut pairs: Vec<(f64, f64)> = freq_info
                .freq_wavenumber
                .iter()
                .zip(freq_info.reduced_mass.iter())
                .map(|(f, m)| (*f, *m))
                .collect();
            pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            pairs
        };
        let pyscf_sorted: Vec<(f64, f64)> = {
            let mut pairs: Vec<(f64, f64)> = golden
                .freq_wavenumber
                .iter()
                .zip(golden.reduced_mass.iter())
                .map(|(f, m)| (*f, *m))
                .collect();
            pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            pairs
        };
        for i in 0..3 {
            let rm = iqcp_sorted[i].1;
            let rm_ref = pyscf_sorted[i].1;
            let err = (rm - rm_ref).abs();
            assert!(
                err < 0.01,
                "H2O reduced mass {}: IQCP {:.5} vs PySCF {:.5}, err {:.5} amu",
                i,
                rm,
                rm_ref,
                err
            );
        }

        // Rotational constants within 0.1%
        for i in 0..3 {
            let err =
                (freq_info.rotational_constants_ghz[i] - golden.rotational_constants_ghz[i]).abs();
            let rel = err / golden.rotational_constants_ghz[i].abs().max(1e-30);
            assert!(
                rel < 1e-3,
                "H2O rot const {}: IQCP {:.4} vs PySCF {:.4} GHz",
                i,
                freq_info.rotational_constants_ghz[i],
                golden.rotational_constants_ghz[i]
            );
        }
    }

    #[test]
    fn test_harmonic_analysis_ch4_golden_hessian() {
        let golden = load_golden("ch4");
        let atoms = golden.atoms();
        let hessian = golden.hessian_matrix();
        let freq_info = harmonic_analysis(&atoms, &hessian).unwrap();

        assert_eq!(freq_info.n_atoms, 5);
        assert_eq!(freq_info.n_modes, 9);
        assert_eq!(freq_info.rot_type, RotorType::SphericalTop);

        let max_err = max_abs_freq_diff(&freq_info.freq_wavenumber, &golden.freq_wavenumber);
        assert!(max_err < 2.0, "CH4 max freq err = {:.4} cm^-1", max_err);

        // Check that 2 triply-degenerate modes exist (2 T2 modes: stretching & bending)
        let mut freq = freq_info.freq_wavenumber.clone();
        freq.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut degenerate_triplets = 0;
        let mut i = 0;
        while i + 2 < freq.len() {
            if (freq[i + 2] - freq[i]).abs() < 1.0 {
                degenerate_triplets += 1;
                i += 3;
            } else {
                i += 1;
            }
        }
        assert!(
            degenerate_triplets >= 2,
            "Expected ≥2 triply-degenerate groups in CH4, found {}",
            degenerate_triplets
        );
    }

    #[test]
    fn test_harmonic_analysis_co2_golden_hessian() {
        let golden = load_golden("co2");
        let atoms = golden.atoms();
        let hessian = golden.hessian_matrix();
        let freq_info = harmonic_analysis(&atoms, &hessian).unwrap();

        assert_eq!(freq_info.n_atoms, 3);
        assert_eq!(freq_info.rot_type, RotorType::Linear);
        assert_eq!(freq_info.n_modes, 4); // 3N-5 = 4

        let max_err = max_abs_freq_diff(&freq_info.freq_wavenumber, &golden.freq_wavenumber);
        assert!(max_err < 2.0, "CO2 max freq err = {:.4} cm^-1", max_err);

        // Bending modes are doubly degenerate
        let mut freq = freq_info.freq_wavenumber.clone();
        freq.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!(
            (freq[1] - freq[0]).abs() < 1.0,
            "CO2 bending degeneracy broken: {} vs {}",
            freq[0],
            freq[1]
        );
    }

    #[test]
    fn test_harmonic_analysis_h2_golden_hessian() {
        let golden = load_golden("h2");
        let atoms = golden.atoms();
        let hessian = golden.hessian_matrix();
        let freq_info = harmonic_analysis(&atoms, &hessian).unwrap();

        assert_eq!(freq_info.n_atoms, 2);
        assert_eq!(freq_info.rot_type, RotorType::Linear);
        assert_eq!(freq_info.n_modes, 1); // 3N-5 = 1

        let max_err = max_abs_freq_diff(&freq_info.freq_wavenumber, &golden.freq_wavenumber);
        assert!(max_err < 2.0, "H2 freq err = {:.4} cm^-1", max_err);
    }

    // -------------------------------------------------------------
    // End-to-end test: IQCP Hessian + IQCP thermo vs PySCF
    //
    // These tests exercise the full pipeline from SCF → analytical
    // Hessian → harmonic analysis, not just the thermo.rs module
    // in isolation. Tolerances are looser to account for IQCP's own
    // Hessian numerical error (~1e-6 Ha/bohr² vs PySCF).
    // -------------------------------------------------------------

    #[test]
    fn test_harmonic_analysis_h2o_rhf_end_to_end() {
        use crate::scf::hessian::rhf_hessian;
        use crate::scf::ScfConfig;

        let golden = load_golden("h2o");
        let atoms = golden.atoms();
        let atom_tuples: Vec<(u8, [f64; 3])> = atoms
            .iter()
            .map(|a| (a.atomic_number, a.position))
            .collect();

        let config = ScfConfig::default();
        let hess_result =
            rhf_hessian(&atom_tuples, "sto-3g", &config).expect("rhf_hessian failed for h2o");

        let freq_info =
            harmonic_analysis(&atoms, &hess_result.hessian).expect("harmonic_analysis failed");

        assert_eq!(freq_info.n_modes, 3);
        assert_eq!(freq_info.rot_type, RotorType::AsymmetricTop);

        // AC9: within ±2 cm⁻¹
        let max_err = max_abs_freq_diff(&freq_info.freq_wavenumber, &golden.freq_wavenumber);
        assert!(
            max_err < 2.0,
            "H2O end-to-end max freq err = {:.4} cm^-1 (AC9 tolerance 2.0)",
            max_err
        );
    }

    #[test]
    fn test_harmonic_analysis_co2_rhf_end_to_end() {
        use crate::scf::hessian::rhf_hessian;
        use crate::scf::ScfConfig;

        let golden = load_golden("co2");
        let atoms = golden.atoms();
        let atom_tuples: Vec<(u8, [f64; 3])> = atoms
            .iter()
            .map(|a| (a.atomic_number, a.position))
            .collect();

        let config = ScfConfig::default();
        let hess_result =
            rhf_hessian(&atom_tuples, "sto-3g", &config).expect("rhf_hessian failed for co2");

        let freq_info =
            harmonic_analysis(&atoms, &hess_result.hessian).expect("harmonic_analysis failed");

        assert_eq!(freq_info.rot_type, RotorType::Linear);
        assert_eq!(freq_info.n_modes, 4);

        let max_err = max_abs_freq_diff(&freq_info.freq_wavenumber, &golden.freq_wavenumber);
        assert!(
            max_err < 2.0,
            "CO2 end-to-end max freq err = {:.4} cm^-1",
            max_err
        );
    }

    #[test]
    fn test_harmonic_analysis_h2_rhf_end_to_end() {
        use crate::scf::hessian::rhf_hessian;
        use crate::scf::ScfConfig;

        let golden = load_golden("h2");
        let atoms = golden.atoms();
        let atom_tuples: Vec<(u8, [f64; 3])> = atoms
            .iter()
            .map(|a| (a.atomic_number, a.position))
            .collect();

        let config = ScfConfig::default();
        let hess_result =
            rhf_hessian(&atom_tuples, "sto-3g", &config).expect("rhf_hessian failed for h2");

        let freq_info =
            harmonic_analysis(&atoms, &hess_result.hessian).expect("harmonic_analysis failed");

        assert_eq!(freq_info.rot_type, RotorType::Linear);
        assert_eq!(freq_info.n_modes, 1);

        let max_err = max_abs_freq_diff(&freq_info.freq_wavenumber, &golden.freq_wavenumber);
        assert!(
            max_err < 2.0,
            "H2 end-to-end freq err = {:.4} cm^-1",
            max_err
        );
    }
}
