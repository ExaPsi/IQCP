//! Primitive-pair breakdown of contracted one-electron integrals
//!
//! This module decomposes a contracted integral I_{ij} = <chi_i | O | chi_j>
//! into its constituent primitive-pair contributions:
//!
//! ```text
//! I_{ij} = sum_{p=1}^{n_p} sum_{q=1}^{n_q} c_p * N_p * c_q * N_q * <g_p | O | g_q>
//! ```
//!
//! where c_p is the contraction coefficient, N_p is the primitive Gaussian
//! normalization factor, and <g_p | O | g_q> is the bare primitive integral.
//!
//! # Pedagogical Purpose
//!
//! This breakdown allows students to see that each matrix element is a
//! structured summation over primitive Gaussian pairs, each with clear
//! physical content. The Gaussian product theorem (Helgaker et al., Ch. 9.2)
//! determines which pairs dominate through the overlap prefactor K_AB.
//!
//! # Invariant
//!
//! `sum(contributions.weighted_contribution) == contracted_value` within 1e-12.
//!
//! # Reference
//!
//! Phase 3 TDD Section 8.2; Phase 3 PRD FR-INT-03; Helgaker et al. (2000) Ch. 9

use serde::{Deserialize, Serialize};

use super::cartesian::cartesian_components;
use super::eri::{primitive_eri, GaussianProduct2e};
use super::kinetic::primitive_kinetic;
use super::nuclear::primitive_nuclear;
use super::overlap::{cartesian_gaussian_normalization, primitive_overlap};
use super::GaussianProduct;
use super::IntegralError;
use crate::basis::BasisSet;
use crate::rys::rys_roots;

// =============================================================================
// Data Types
// =============================================================================

/// A single primitive-pair contribution to a contracted integral.
///
/// Represents one (p, q) pair in the double sum:
/// `I_ij = sum_p sum_q c_p * N_p * c_q * N_q * <g_p | O | g_q>`
///
/// # Fields
///
/// - `prim_indices`: Which primitives within their shells (0-based)
/// - `exponents`: [alpha_p, alpha_q]
/// - `coefficients`: [c_p, c_q] (raw contraction coefficients, without normalization)
/// - `norm_coefficients`: [c_p * N_p, c_q * N_q] (effective coefficients with normalization)
/// - `primitive_value`: The bare primitive integral <g_p | O | g_q> (unnormalized)
/// - `weighted_contribution`: c_p * N_p * c_q * N_q * <g_p | O | g_q>
///
/// # Reference
///
/// Phase 3 TDD Section 8.2, PrimitiveContribution struct
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PrimitiveContribution {
    /// Primitive indices within their respective shells [p, q] (0-based)
    pub prim_indices: [usize; 2],
    /// Primitive exponents [alpha_p, alpha_q]
    pub exponents: [f64; 2],
    /// Raw contraction coefficients [c_p, c_q] (without normalization)
    pub coefficients: [f64; 2],
    /// Contraction coefficients including normalization [c_p * N_p, c_q * N_q]
    pub norm_coefficients: [f64; 2],
    /// The bare primitive integral value (before contraction weighting)
    pub primitive_value: f64,
    /// The weighted contribution: c_p * N_p * c_q * N_q * primitive_value
    pub weighted_contribution: f64,
}

/// Decomposition of a contracted integral into primitive-pair contributions.
///
/// This struct is the return type of `integral_with_breakdown`, providing
/// both the contracted integral value and the individual primitive-pair
/// contributions that sum to it.
///
/// # Invariant
///
/// `sum(pc.weighted_contribution for pc in primitive_contributions) == contracted_value`
/// within 1e-12 (floating-point arithmetic).
///
/// # Reference
///
/// Phase 3 TDD Section 8.2; Phase 3 PRD FR-INT-03
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct IntegralBreakdown {
    /// The contracted integral value (sum of all primitive contributions)
    pub contracted_value: f64,
    /// Integral type: "overlap", "kinetic", "nuclear", "hcore"
    pub integral_type: String,
    /// Basis function indices [i, j] (0-based)
    pub indices: [usize; 2],
    /// Basis function labels for the two functions (e.g., ["O 1s", "H1 1s"])
    pub labels: [String; 2],
    /// Individual primitive-pair contributions, sorted by |weighted_contribution| descending
    pub primitive_contributions: Vec<PrimitiveContribution>,
    /// Number of primitives in shell i
    pub n_prim_i: usize,
    /// Number of primitives in shell j
    pub n_prim_j: usize,
}

// =============================================================================
// ERI Breakdown Data Types (US-059)
// =============================================================================

/// A single primitive-quartet contribution to a contracted ERI.
///
/// Represents one (p, q, r, s) quartet in the quadruple sum:
/// `(ij|kl) = sum_p sum_q sum_r sum_s c_p*N_p * c_q*N_q * c_r*N_r * c_s*N_s * (g_p g_q|g_r g_s)`
///
/// # Reference
///
/// Phase 3 PRD FR-INT-04; Szabo & Ostlund Eq. 3.155
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct EriPrimitiveContribution {
    /// Primitive indices within their respective shells [p, q, r, s] (0-based)
    pub prim_indices: [usize; 4],
    /// Primitive exponents [alpha_p, alpha_q, alpha_r, alpha_s]
    pub exponents: [f64; 4],
    /// Raw contraction coefficients [c_p, c_q, c_r, c_s] (without normalization)
    pub coefficients: [f64; 4],
    /// Contraction coefficients including normalization [c_p*N_p, c_q*N_q, c_r*N_r, c_s*N_s]
    pub norm_coefficients: [f64; 4],
    /// The bare primitive ERI value (before contraction weighting)
    pub primitive_value: f64,
    /// The weighted contribution: product(norm_coefficients) * primitive_value
    pub weighted_contribution: f64,
    /// T parameter for this primitive quartet: rho * |P-Q|^2
    pub t_parameter: f64,
}

/// Method used to compute the ERI.
///
/// Pedagogical indicator connecting ERI computation to Module B (Rys) concepts.
/// For (ss|ss) quartets, the integral reduces to a single Boys function F_0(T),
/// which is mathematically equivalent to single-root Rys quadrature. For higher
/// angular momentum, multi-root Rys quadrature is required.
///
/// # Reference
///
/// Dupuis, Rys & King (1976), J. Chem. Phys. 65, 111.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(tag = "type", rename_all = "camelCase")]
pub enum EriMethod {
    /// All s-type shells: (ss|ss). Equivalent to F_0(T).
    #[serde(rename_all = "camelCase")]
    BoysFunction {
        /// T parameter for the representative primitive quartet (most diffuse)
        t_parameter: f64,
    },
    /// Higher angular momentum: uses Rys quadrature with nroots >= 2
    #[serde(rename_all = "camelCase")]
    RysQuadrature {
        /// Number of Rys quadrature roots
        nroots: usize,
        /// T parameter for the representative primitive quartet (most diffuse)
        t_parameter: f64,
        /// Rys roots for the representative primitive quartet
        roots: Vec<f64>,
        /// Rys weights for the representative primitive quartet
        weights: Vec<f64>,
    },
}

/// Decomposition of a contracted ERI into primitive-quartet contributions.
///
/// Analogous to `IntegralBreakdown` for one-electron integrals, but with
/// 4 indices and primitive-quartet (not pair) contributions.
///
/// # Invariant
///
/// `sum(pc.weighted_contribution for pc in contributions) == contracted_value`
/// within 1e-12 (floating-point arithmetic).
///
/// # Reference
///
/// Phase 3 PRD FR-INT-04; Dupuis, Rys & King (1976)
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct EriBreakdown {
    /// The contracted ERI value (sum of all primitive contributions)
    pub contracted_value: f64,
    /// Basis function indices [i, j, k, l] (0-based)
    pub indices: [usize; 4],
    /// Basis function labels (e.g., ["O 1s", "O 2s", "H1 1s", "H1 1s"])
    pub labels: [String; 4],
    /// Computational method used (Boys function or Rys quadrature)
    pub method: EriMethod,
    /// Individual primitive-quartet contributions, sorted by |weighted_contribution| descending
    pub contributions: Vec<EriPrimitiveContribution>,
    /// Number of primitives in each shell [n_i, n_j, n_k, n_l]
    pub n_primitives: [usize; 4],
    /// Total angular momentum: L_i + L_j + L_k + L_l
    pub total_angular_momentum: u32,
    /// Number of Rys quadrature roots: L_total/2 + 1
    pub nroots: usize,
}

// =============================================================================
// Helper: Find shell and component for a basis function index
// =============================================================================

/// Map a linear basis function index to (shell_index, component_index_within_shell).
///
/// Iterates through basis.shells, accumulating `n_basis_functions()`, to find
/// which shell contains the given basis function index and what Cartesian
/// component it corresponds to within that shell.
///
/// # Returns
///
/// `(shell_idx, component_idx)` or `Err` if index is out of range.
fn find_shell_and_component(
    basis: &BasisSet,
    bf_index: usize,
) -> Result<(usize, usize), IntegralError> {
    if bf_index >= basis.n_basis {
        return Err(IntegralError::InvalidBasis(format!(
            "Basis function index {} out of range (n_basis = {})",
            bf_index, basis.n_basis
        )));
    }

    let mut offset = 0;
    for (shell_idx, shell) in basis.shells.iter().enumerate() {
        let n_funcs = shell.n_basis_functions();
        if bf_index < offset + n_funcs {
            return Ok((shell_idx, bf_index - offset));
        }
        offset += n_funcs;
    }

    // Should be unreachable if n_basis is consistent with shells
    Err(IntegralError::InvalidBasis(format!(
        "Basis function index {} not found in any shell (internal inconsistency)",
        bf_index
    )))
}

/// Generate a label for a basis function at a given index.
///
/// Uses the atom symbol and shell type to produce labels like
/// "H1 1s", "O 2px", etc.
fn generate_bf_label(basis: &BasisSet, bf_index: usize) -> String {
    use crate::basis::AngularMomentum;
    use std::collections::HashMap;

    // Count atoms of each element for numbering
    let mut element_counts: HashMap<&str, usize> = HashMap::new();
    for atom in &basis.atoms {
        *element_counts.entry(atom.symbol.as_str()).or_insert(0) += 1;
    }

    // Build atom labels
    let mut element_instance: HashMap<&str, usize> = HashMap::new();
    let mut atom_labels: Vec<String> = Vec::with_capacity(basis.atoms.len());
    for atom in &basis.atoms {
        let sym = atom.symbol.as_str();
        let count = element_counts[sym];
        if count > 1 {
            let instance = element_instance.entry(sym).or_insert(0);
            *instance += 1;
            atom_labels.push(format!("{}{}", sym, *instance));
        } else {
            atom_labels.push(sym.to_string());
        }
    }

    // Track shell counts per atom for principal quantum number
    let mut atom_s_count: Vec<usize> = vec![0; basis.atoms.len()];
    let mut atom_p_count: Vec<usize> = vec![0; basis.atoms.len()];
    let mut atom_d_count: Vec<usize> = vec![0; basis.atoms.len()];

    let cart_p = ["px", "py", "pz"];
    let cart_d = ["dxx", "dxy", "dxz", "dyy", "dyz", "dzz"];

    let mut offset = 0;
    for shell in &basis.shells {
        let atom_idx = shell.atom_idx;
        let n_funcs = shell.n_basis_functions();

        let (n_label, comp_labels): (usize, &[&str]) = match shell.angular_momentum {
            AngularMomentum::S => {
                atom_s_count[atom_idx] += 1;
                (atom_s_count[atom_idx], &["s"])
            }
            AngularMomentum::P => {
                atom_p_count[atom_idx] += 1;
                (atom_p_count[atom_idx] + 1, &cart_p)
            }
            AngularMomentum::D => {
                atom_d_count[atom_idx] += 1;
                (atom_d_count[atom_idx] + 2, &cart_d)
            }
        };

        if bf_index >= offset && bf_index < offset + n_funcs {
            let comp_idx = bf_index - offset;
            return format!(
                "{} {}{}",
                atom_labels[atom_idx], n_label, comp_labels[comp_idx]
            );
        }

        offset += n_funcs;
    }

    format!("bf{}", bf_index)
}

// =============================================================================
// Core Function: integral_with_breakdown
// =============================================================================

/// Compute a contracted one-electron integral and decompose it into primitive-pair contributions.
///
/// For a contracted integral I_ij = <chi_i | O | chi_j> where shell i has n_p primitives
/// and shell j has n_q primitives, this function returns n_p * n_q primitive contributions.
///
/// # Arguments
///
/// * `basis` - The molecular basis set
/// * `integral_type` - Type of integral: "overlap", "kinetic", "nuclear", or "hcore"
/// * `i` - Basis function index (row in the integral matrix)
/// * `j` - Basis function index (column in the integral matrix)
///
/// # Returns
///
/// An `IntegralBreakdown` containing the contracted value and all primitive-pair contributions
/// sorted by |weighted_contribution| descending.
///
/// # Errors
///
/// Returns `IntegralError` if indices are out of range or integral type is invalid.
///
/// # Invariant
///
/// `sum(contributions.weighted_contribution) == contracted_value` within 1e-12.
///
/// # Reference
///
/// Phase 3 TDD Section 8.2; Phase 3 PRD FR-INT-03
pub fn integral_with_breakdown(
    basis: &BasisSet,
    integral_type: &str,
    i: usize,
    j: usize,
) -> Result<IntegralBreakdown, IntegralError> {
    // Validate integral type
    match integral_type {
        "overlap" | "kinetic" | "nuclear" | "hcore" => {}
        other => {
            return Err(IntegralError::InvalidBasis(format!(
                "Invalid integral type '{}'. Must be 'overlap', 'kinetic', 'nuclear', or 'hcore'.",
                other
            )));
        }
    }

    // Map basis function indices to shells and Cartesian components
    let (shell_idx_i, comp_idx_i) = find_shell_and_component(basis, i)?;
    let (shell_idx_j, comp_idx_j) = find_shell_and_component(basis, j)?;

    let shell_i = &basis.shells[shell_idx_i];
    let shell_j = &basis.shells[shell_idx_j];

    // Get Cartesian powers for the specific components
    let comps_i = cartesian_components(shell_i.l_value())?;
    let comps_j = cartesian_components(shell_j.l_value())?;

    let pow_i = &comps_i[comp_idx_i];
    let pow_j = &comps_j[comp_idx_j];

    let n_prim_i = shell_i.n_primitives();
    let n_prim_j = shell_j.n_primitives();

    // Generate labels
    let label_i = generate_bf_label(basis, i);
    let label_j = generate_bf_label(basis, j);

    // Iterate over all primitive pairs
    let mut contributions = Vec::with_capacity(n_prim_i * n_prim_j);

    for (p_idx, prim_a) in shell_i.primitives.iter().enumerate() {
        for (q_idx, prim_b) in shell_j.primitives.iter().enumerate() {
            // Compute Gaussian product
            let gp = GaussianProduct::new(
                prim_a.exponent,
                &shell_i.center,
                prim_b.exponent,
                &shell_j.center,
            );

            // Compute Cartesian Gaussian normalization factors
            let norm_a = cartesian_gaussian_normalization(prim_a.exponent, pow_i);
            let norm_b = cartesian_gaussian_normalization(prim_b.exponent, pow_j);

            // Compute the bare primitive integral depending on type
            let prim_value = match integral_type {
                "overlap" => primitive_overlap(&gp, pow_i, pow_j),

                "kinetic" => primitive_kinetic(&gp, pow_i, pow_j, prim_b.exponent),

                "nuclear" => {
                    // Sum over all nuclei, matching nuclear_matrix logic
                    let mut v_total = 0.0;
                    for atom in &basis.atoms {
                        let v_prim = primitive_nuclear(&gp, pow_i, pow_j, &atom.position);
                        let z_factor = atom.atomic_number as f64;
                        v_total += v_prim * z_factor;
                    }
                    v_total
                }

                "hcore" => {
                    // Kinetic + nuclear attraction
                    let t_prim = primitive_kinetic(&gp, pow_i, pow_j, prim_b.exponent);

                    let mut v_total = 0.0;
                    for atom in &basis.atoms {
                        let v_prim = primitive_nuclear(&gp, pow_i, pow_j, &atom.position);
                        let z_factor = atom.atomic_number as f64;
                        v_total += v_prim * z_factor;
                    }

                    t_prim + v_total
                }

                _ => unreachable!(), // Already validated above
            };

            // Weighted contribution = c_p * N_p * c_q * N_q * <g_p | O | g_q>
            let c_p = prim_a.coefficient;
            let c_q = prim_b.coefficient;
            let cn_p = c_p * norm_a;
            let cn_q = c_q * norm_b;
            let weighted = cn_p * cn_q * prim_value;

            contributions.push(PrimitiveContribution {
                prim_indices: [p_idx, q_idx],
                exponents: [prim_a.exponent, prim_b.exponent],
                coefficients: [c_p, c_q],
                norm_coefficients: [cn_p, cn_q],
                primitive_value: prim_value,
                weighted_contribution: weighted,
            });
        }
    }

    // Sort by |weighted_contribution| descending
    contributions.sort_by(|a, b| {
        b.weighted_contribution
            .abs()
            .partial_cmp(&a.weighted_contribution.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Sum all weighted contributions to get the contracted value
    let contracted_value: f64 = contributions.iter().map(|c| c.weighted_contribution).sum();

    Ok(IntegralBreakdown {
        contracted_value,
        integral_type: integral_type.to_string(),
        indices: [i, j],
        labels: [label_i, label_j],
        primitive_contributions: contributions,
        n_prim_i,
        n_prim_j,
    })
}

// =============================================================================
// ERI Breakdown Function (US-059)
// =============================================================================

/// Compute a contracted ERI and decompose it into primitive-quartet contributions.
///
/// For a contracted ERI (ij|kl), this function returns all primitive-quartet
/// contributions sorted by |weighted_contribution| descending. It also determines
/// the computational method (Boys function for (ss|ss) or Rys quadrature for
/// higher angular momentum) and provides representative Rys roots/weights for
/// the most diffuse primitive quartet.
///
/// # Arguments
///
/// * `basis` - The molecular basis set
/// * `i`, `j`, `k`, `l` - Basis function indices (0-based)
///
/// # Returns
///
/// An `EriBreakdown` containing the contracted value, all primitive-quartet
/// contributions sorted by |weighted_contribution| descending, and method
/// metadata (Boys/Rys with roots, weights, T parameter).
///
/// # Errors
///
/// Returns `IntegralError` if any index is out of range.
///
/// # Invariant
///
/// `sum(contributions.weighted_contribution) == contracted_value` within 1e-12.
///
/// # Reference
///
/// Phase 3 PRD FR-INT-04; Szabo & Ostlund Eq. 3.155; Dupuis, Rys & King (1976)
pub fn eri_with_breakdown(
    basis: &BasisSet,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Result<EriBreakdown, IntegralError> {
    // 1. Map each basis function index to its shell and Cartesian component
    let (shell_idx_i, comp_idx_i) = find_shell_and_component(basis, i)?;
    let (shell_idx_j, comp_idx_j) = find_shell_and_component(basis, j)?;
    let (shell_idx_k, comp_idx_k) = find_shell_and_component(basis, k)?;
    let (shell_idx_l, comp_idx_l) = find_shell_and_component(basis, l)?;

    let shell_i = &basis.shells[shell_idx_i];
    let shell_j = &basis.shells[shell_idx_j];
    let shell_k = &basis.shells[shell_idx_k];
    let shell_l = &basis.shells[shell_idx_l];

    // Get Cartesian powers for the specific components
    let comps_i = cartesian_components(shell_i.l_value())?;
    let comps_j = cartesian_components(shell_j.l_value())?;
    let comps_k = cartesian_components(shell_k.l_value())?;
    let comps_l = cartesian_components(shell_l.l_value())?;

    let pow_i = &comps_i[comp_idx_i];
    let pow_j = &comps_j[comp_idx_j];
    let pow_k = &comps_k[comp_idx_k];
    let pow_l = &comps_l[comp_idx_l];

    let n_prim_i = shell_i.n_primitives();
    let n_prim_j = shell_j.n_primitives();
    let n_prim_k = shell_k.n_primitives();
    let n_prim_l = shell_l.n_primitives();

    // 2. Angular momentum and Rys root count
    let l_i = shell_i.l_value();
    let l_j = shell_j.l_value();
    let l_k = shell_k.l_value();
    let l_l = shell_l.l_value();
    let l_total = l_i + l_j + l_k + l_l;
    let nroots = (l_total / 2 + 1) as usize;

    // Generate labels
    let label_i = generate_bf_label(basis, i);
    let label_j = generate_bf_label(basis, j);
    let label_k = generate_bf_label(basis, k);
    let label_l = generate_bf_label(basis, l);

    // 3. Iterate over all primitive quartets
    let mut contributions = Vec::with_capacity(n_prim_i * n_prim_j * n_prim_k * n_prim_l);

    // Track the most diffuse quartet (smallest rho => smallest T) for representative
    // Rys roots/weights display
    let mut min_rho: f64 = f64::MAX;
    let mut representative_t: f64 = 0.0;

    for (p_idx, prim_a) in shell_i.primitives.iter().enumerate() {
        let norm_a = cartesian_gaussian_normalization(prim_a.exponent, pow_i);

        for (q_idx, prim_b) in shell_j.primitives.iter().enumerate() {
            let norm_b = cartesian_gaussian_normalization(prim_b.exponent, pow_j);

            for (r_idx, prim_c) in shell_k.primitives.iter().enumerate() {
                let norm_c = cartesian_gaussian_normalization(prim_c.exponent, pow_k);

                for (s_idx, prim_d) in shell_l.primitives.iter().enumerate() {
                    let norm_d = cartesian_gaussian_normalization(prim_d.exponent, pow_l);

                    // Compute GaussianProduct2e for this primitive quartet
                    let gp2e = GaussianProduct2e::new(
                        prim_a.exponent,
                        &shell_i.center,
                        prim_b.exponent,
                        &shell_j.center,
                        prim_c.exponent,
                        &shell_k.center,
                        prim_d.exponent,
                        &shell_l.center,
                    );

                    // Track representative quartet (most diffuse = smallest rho)
                    if gp2e.rho < min_rho {
                        min_rho = gp2e.rho;
                        representative_t = gp2e.t;
                    }

                    // Compute the bare primitive ERI
                    let prim_val = primitive_eri(&gp2e, pow_i, pow_j, pow_k, pow_l);

                    // Weighted contribution = c_p*N_p * c_q*N_q * c_r*N_r * c_s*N_s * (g_p g_q|g_r g_s)
                    let c_p = prim_a.coefficient;
                    let c_q = prim_b.coefficient;
                    let c_r = prim_c.coefficient;
                    let c_s = prim_d.coefficient;
                    let cn_p = c_p * norm_a;
                    let cn_q = c_q * norm_b;
                    let cn_r = c_r * norm_c;
                    let cn_s = c_s * norm_d;
                    let weighted = cn_p * cn_q * cn_r * cn_s * prim_val;

                    contributions.push(EriPrimitiveContribution {
                        prim_indices: [p_idx, q_idx, r_idx, s_idx],
                        exponents: [
                            prim_a.exponent,
                            prim_b.exponent,
                            prim_c.exponent,
                            prim_d.exponent,
                        ],
                        coefficients: [c_p, c_q, c_r, c_s],
                        norm_coefficients: [cn_p, cn_q, cn_r, cn_s],
                        primitive_value: prim_val,
                        weighted_contribution: weighted,
                        t_parameter: gp2e.t,
                    });
                }
            }
        }
    }

    // 4. Sort contributions by |weighted_contribution| descending
    contributions.sort_by(|a, b| {
        b.weighted_contribution
            .abs()
            .partial_cmp(&a.weighted_contribution.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // 5. Sum all weighted contributions to get the contracted value
    let contracted_value: f64 = contributions.iter().map(|c| c.weighted_contribution).sum();

    // 6. Determine method metadata
    // For the representative Rys roots/weights, use the most diffuse primitive quartet
    // (smallest rho, giving the smallest T and largest spatial extent).
    let method = if l_total == 0 {
        // (ss|ss): pedagogically labeled as Boys function F_0(T)
        EriMethod::BoysFunction {
            t_parameter: representative_t,
        }
    } else {
        // Higher angular momentum: compute Rys roots/weights for representative T
        let (roots, weights) = match rys_roots(nroots, representative_t) {
            Ok(rys_result) => (rys_result.roots.to_vec(), rys_result.weights.to_vec()),
            Err(_) => {
                // Fallback: report empty roots/weights if Rys fails at this T
                (vec![], vec![])
            }
        };

        EriMethod::RysQuadrature {
            nroots,
            t_parameter: representative_t,
            roots,
            weights,
        }
    };

    Ok(EriBreakdown {
        contracted_value,
        indices: [i, j, k, l],
        labels: [label_i, label_j, label_k, label_l],
        method,
        contributions,
        n_primitives: [n_prim_i, n_prim_j, n_prim_k, n_prim_l],
        total_angular_momentum: l_total,
        nroots,
    })
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;
    use crate::basis::Atom;
    use crate::integrals::{kinetic_matrix, nuclear_matrix, overlap_matrix};
    use approx::assert_abs_diff_eq;

    // -------------------------------------------------------------------------
    // Test helpers
    // -------------------------------------------------------------------------

    fn h2_basis() -> BasisSet {
        let h1 = Atom::new(1, [0.0, 0.0, 0.0]).unwrap();
        let h2 = Atom::new(1, [0.0, 0.0, 1.3984]).unwrap();
        BasisSet::build(vec![h1, h2], "sto-3g").unwrap()
    }

    fn h2o_basis() -> BasisSet {
        let o = Atom::new(8, [0.0, 0.0, 0.2216656303019316]).unwrap();
        let h1 = Atom::new(1, [0.0, 1.430929534330371, -0.8866625212077263]).unwrap();
        let h2 = Atom::new(1, [0.0, -1.430929534330371, -0.8866625212077263]).unwrap();
        BasisSet::build(vec![o, h1, h2], "sto-3g").unwrap()
    }

    fn lih_basis() -> BasisSet {
        let li = Atom::new(3, [0.0, 0.0, 0.0]).unwrap();
        let h = Atom::new(1, [0.0, 0.0, 3.015437]).unwrap();
        BasisSet::build(vec![li, h], "sto-3g").unwrap()
    }

    // =========================================================================
    // PB-R1: Sum of weighted contributions equals contracted value (overlap)
    // =========================================================================

    #[test]
    fn pb_r1_sum_equals_contracted_overlap_h2() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();

        let sum: f64 = bd
            .primitive_contributions
            .iter()
            .map(|c| c.weighted_contribution)
            .sum();

        assert_abs_diff_eq!(sum, bd.contracted_value, epsilon = 1e-15);

        // Also verify against full matrix computation
        let s = overlap_matrix(&basis);
        assert_abs_diff_eq!(bd.contracted_value, s[0 * 2 + 1], epsilon = 1e-12);
    }

    // =========================================================================
    // PB-R2: Sum invariant for kinetic integral
    // =========================================================================

    #[test]
    fn pb_r2_sum_equals_contracted_kinetic_h2() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "kinetic", 0, 1).unwrap();

        let sum: f64 = bd
            .primitive_contributions
            .iter()
            .map(|c| c.weighted_contribution)
            .sum();

        assert_abs_diff_eq!(sum, bd.contracted_value, epsilon = 1e-15);

        let t = kinetic_matrix(&basis);
        assert_abs_diff_eq!(bd.contracted_value, t[0 * 2 + 1], epsilon = 1e-12);
    }

    // =========================================================================
    // PB-R3: Sum invariant for nuclear integral
    // =========================================================================

    #[test]
    fn pb_r3_sum_equals_contracted_nuclear_h2() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "nuclear", 0, 1).unwrap();

        let sum: f64 = bd
            .primitive_contributions
            .iter()
            .map(|c| c.weighted_contribution)
            .sum();

        assert_abs_diff_eq!(sum, bd.contracted_value, epsilon = 1e-15);

        let v = nuclear_matrix(&basis);
        assert_abs_diff_eq!(bd.contracted_value, v[0 * 2 + 1], epsilon = 1e-5);
    }

    // =========================================================================
    // PB-R4: Sum invariant for Hcore integral
    // =========================================================================

    #[test]
    fn pb_r4_sum_equals_contracted_hcore_h2() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "hcore", 0, 1).unwrap();

        let sum: f64 = bd
            .primitive_contributions
            .iter()
            .map(|c| c.weighted_contribution)
            .sum();

        assert_abs_diff_eq!(sum, bd.contracted_value, epsilon = 1e-15);

        // Hcore = T + V
        let t = kinetic_matrix(&basis);
        let v = nuclear_matrix(&basis);
        let hcore_ref = t[0 * 2 + 1] + v[0 * 2 + 1];
        assert_abs_diff_eq!(bd.contracted_value, hcore_ref, epsilon = 1e-5);
    }

    // =========================================================================
    // PB-R5: Sum invariant for diagonal element (S[0,0] = 1)
    // =========================================================================

    #[test]
    fn pb_r5_sum_equals_contracted_diagonal_overlap() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 0).unwrap();

        // S[0,0] should be ~1.0 for normalized basis functions.
        // The tolerance accounts for floating-point accumulation in the primitive-pair
        // summation (9 terms for STO-3G). The overlap_matrix function produces
        // the same result via the same shell_overlap code path.
        let s = overlap_matrix(&basis);
        assert_abs_diff_eq!(bd.contracted_value, s[0], epsilon = 1e-12);
        assert_abs_diff_eq!(bd.contracted_value, 1.0, epsilon = 1e-9);
    }

    // =========================================================================
    // PB-R6: Sum invariant for H2O system
    // =========================================================================

    #[test]
    fn pb_r6_sum_invariant_h2o() {
        let basis = h2o_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 5).unwrap();

        let sum: f64 = bd
            .primitive_contributions
            .iter()
            .map(|c| c.weighted_contribution)
            .sum();

        assert_abs_diff_eq!(sum, bd.contracted_value, epsilon = 1e-15);

        let s = overlap_matrix(&basis);
        let n = basis.n_basis;
        assert_abs_diff_eq!(bd.contracted_value, s[0 * n + 5], epsilon = 1e-12);
    }

    // =========================================================================
    // PB-R7: Sum invariant for p-type shell (H2O: O 2px - H1 1s kinetic)
    // =========================================================================

    #[test]
    fn pb_r7_sum_invariant_p_type() {
        let basis = h2o_basis();
        // Basis function 2 is O 2px, function 5 is H1 1s
        let bd = integral_with_breakdown(&basis, "kinetic", 2, 5).unwrap();

        let sum: f64 = bd
            .primitive_contributions
            .iter()
            .map(|c| c.weighted_contribution)
            .sum();

        assert_abs_diff_eq!(sum, bd.contracted_value, epsilon = 1e-15);

        let t = kinetic_matrix(&basis);
        let n = basis.n_basis;
        assert_abs_diff_eq!(bd.contracted_value, t[2 * n + 5], epsilon = 1e-9);
    }

    // =========================================================================
    // PB-R8 through PB-R11: Primitive count tests
    // =========================================================================

    #[test]
    fn pb_r8_sto3g_ss_pair_gives_9_contributions() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();
        assert_eq!(bd.primitive_contributions.len(), 9); // 3 x 3
    }

    #[test]
    fn pb_r9_sto3g_ss_same_shell_gives_9() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 0).unwrap();
        assert_eq!(bd.primitive_contributions.len(), 9);
    }

    #[test]
    fn pb_r10_h2o_o1s_h1s_gives_9() {
        let basis = h2o_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 5).unwrap();
        assert_eq!(bd.primitive_contributions.len(), 9); // O 1s (3 prims) x H1 1s (3 prims)
    }

    #[test]
    fn pb_r11_lih_inner_outer_ss_gives_9() {
        let basis = lih_basis();
        // Li 1s is bf 0, Li 2s is bf 1
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();
        assert_eq!(bd.primitive_contributions.len(), 9); // Li 1s (3 prims) x Li 2s (3 prims)
    }

    // =========================================================================
    // PB-R12: Sorting test
    // =========================================================================

    #[test]
    fn pb_r12_contributions_sorted_by_magnitude_descending() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();

        for window in bd.primitive_contributions.windows(2) {
            assert!(
                window[0].weighted_contribution.abs() >= window[1].weighted_contribution.abs(),
                "Contributions not sorted: {} >= {} failed",
                window[0].weighted_contribution.abs(),
                window[1].weighted_contribution.abs()
            );
        }
    }

    // =========================================================================
    // PB-R13 through PB-R16: Metadata tests
    // =========================================================================

    #[test]
    fn pb_r13_integral_type_matches_request() {
        let basis = h2_basis();

        let bd_s = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();
        assert_eq!(bd_s.integral_type, "overlap");

        let bd_t = integral_with_breakdown(&basis, "kinetic", 0, 1).unwrap();
        assert_eq!(bd_t.integral_type, "kinetic");

        let bd_v = integral_with_breakdown(&basis, "nuclear", 0, 1).unwrap();
        assert_eq!(bd_v.integral_type, "nuclear");

        let bd_h = integral_with_breakdown(&basis, "hcore", 0, 1).unwrap();
        assert_eq!(bd_h.integral_type, "hcore");
    }

    #[test]
    fn pb_r14_indices_match_request() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();
        assert_eq!(bd.indices, [0, 1]);
    }

    #[test]
    fn pb_r15_n_prim_correct() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();
        assert_eq!(bd.n_prim_i, 3);
        assert_eq!(bd.n_prim_j, 3);
    }

    #[test]
    fn pb_r16_labels_populated() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();
        assert!(bd.labels[0].contains("H"));
        assert!(bd.labels[0].contains("1s"));
        assert!(bd.labels[1].contains("H"));
        assert!(bd.labels[1].contains("1s"));
    }

    // =========================================================================
    // PB-R17 through PB-R19: Error handling tests
    // =========================================================================

    #[test]
    fn pb_r17_out_of_range_index_i() {
        let basis = h2_basis();
        let result = integral_with_breakdown(&basis, "overlap", 5, 0);
        assert!(result.is_err());
    }

    #[test]
    fn pb_r18_out_of_range_index_j() {
        let basis = h2_basis();
        let result = integral_with_breakdown(&basis, "overlap", 0, 5);
        assert!(result.is_err());
    }

    #[test]
    fn pb_r19_invalid_integral_type() {
        let basis = h2_basis();
        let result = integral_with_breakdown(&basis, "invalid", 0, 1);
        assert!(result.is_err());
    }

    // =========================================================================
    // PB-R20 through PB-R21: Exponent and coefficient consistency
    // =========================================================================

    #[test]
    fn pb_r20_exponents_match_basis_set() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();

        let expected_exponents = [3.4252509100, 0.6239137300, 0.1688554000];

        for pc in &bd.primitive_contributions {
            assert!(
                expected_exponents
                    .iter()
                    .any(|&e| (pc.exponents[0] - e).abs() < 1e-6),
                "Unexpected exponent_i: {}",
                pc.exponents[0]
            );
            assert!(
                expected_exponents
                    .iter()
                    .any(|&e| (pc.exponents[1] - e).abs() < 1e-6),
                "Unexpected exponent_j: {}",
                pc.exponents[1]
            );
        }
    }

    #[test]
    fn pb_r21_coefficients_match_basis_set() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();

        let expected_coefficients = [0.1543289707, 0.5353281424, 0.4446345420];

        for pc in &bd.primitive_contributions {
            assert!(
                expected_coefficients
                    .iter()
                    .any(|&c| (pc.coefficients[0] - c).abs() < 1e-6),
                "Unexpected coefficient_i: {}",
                pc.coefficients[0]
            );
            assert!(
                expected_coefficients
                    .iter()
                    .any(|&c| (pc.coefficients[1] - c).abs() < 1e-6),
                "Unexpected coefficient_j: {}",
                pc.coefficients[1]
            );
        }
    }

    // =========================================================================
    // PB-G1 through PB-G7: PySCF golden value tests
    // =========================================================================

    #[test]
    fn pb_g1_h2_overlap_01_contracted_matches_pyscf() {
        // PySCF 2.11.0: S[0,1] = 6.598721980070731e-01
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();
        let s = overlap_matrix(&basis);
        assert_abs_diff_eq!(bd.contracted_value, s[0 * 2 + 1], epsilon = 1e-12);
    }

    #[test]
    fn pb_g5_h2_overlap_00_diagonal_equals_1() {
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 0).unwrap();
        // Tolerance matches the overlap_matrix golden test (1e-9)
        assert_abs_diff_eq!(bd.contracted_value, 1.0, epsilon = 1e-9);
    }

    #[test]
    fn pb_g6_h2_kinetic_01_contracted_matches_pyscf() {
        // PySCF 2.11.0: T[0,1] = 2.369594228630360e-01
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "kinetic", 0, 1).unwrap();
        let t = kinetic_matrix(&basis);
        assert_abs_diff_eq!(bd.contracted_value, t[0 * 2 + 1], epsilon = 1e-12);
    }

    #[test]
    fn pb_g7_h2_nuclear_01_contracted_matches_pyscf() {
        // PySCF 2.11.0: V[0,1] = -1.196333536602375e+00
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "nuclear", 0, 1).unwrap();
        let v = nuclear_matrix(&basis);
        assert_abs_diff_eq!(bd.contracted_value, v[0 * 2 + 1], epsilon = 1e-5);
    }

    // =========================================================================
    // Weighted contribution verification
    // =========================================================================

    #[test]
    fn weighted_equals_norm_coefficients_times_primitive() {
        // Verify the relationship: weighted = cn_p * cn_q * primitive_value
        let basis = h2_basis();
        let bd = integral_with_breakdown(&basis, "overlap", 0, 1).unwrap();

        for pc in &bd.primitive_contributions {
            let expected = pc.norm_coefficients[0] * pc.norm_coefficients[1] * pc.primitive_value;
            assert_abs_diff_eq!(pc.weighted_contribution, expected, epsilon = 1e-15);
        }
    }

    // =========================================================================
    // H2O multi-type test
    // =========================================================================

    #[test]
    fn all_integral_types_consistent_h2o() {
        let basis = h2o_basis();
        let n = basis.n_basis;

        let s = overlap_matrix(&basis);
        let t = kinetic_matrix(&basis);
        let v = nuclear_matrix(&basis);

        // Test a few elements across all types
        for &(ii, jj) in &[(0, 0), (0, 1), (0, 5), (2, 5), (3, 4)] {
            let bd_s = integral_with_breakdown(&basis, "overlap", ii, jj).unwrap();
            assert_abs_diff_eq!(bd_s.contracted_value, s[ii * n + jj], epsilon = 1e-12);

            let bd_t = integral_with_breakdown(&basis, "kinetic", ii, jj).unwrap();
            assert_abs_diff_eq!(bd_t.contracted_value, t[ii * n + jj], epsilon = 1e-9);

            let bd_v = integral_with_breakdown(&basis, "nuclear", ii, jj).unwrap();
            assert_abs_diff_eq!(bd_v.contracted_value, v[ii * n + jj], epsilon = 1e-5);

            let bd_h = integral_with_breakdown(&basis, "hcore", ii, jj).unwrap();
            let hcore_ref = t[ii * n + jj] + v[ii * n + jj];
            assert_abs_diff_eq!(bd_h.contracted_value, hcore_ref, epsilon = 1e-5);
        }
    }

    // =========================================================================
    // ERI Breakdown Tests (US-059)
    // =========================================================================

    // =========================================================================
    // EB-R1: Sum of weighted contributions equals contracted value
    // =========================================================================

    #[test]
    fn eb_r1_sum_equals_contracted_h2_0000() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();

        let sum: f64 = bd
            .contributions
            .iter()
            .map(|c| c.weighted_contribution)
            .sum();
        assert_abs_diff_eq!(sum, bd.contracted_value, epsilon = 1e-12);
    }

    // =========================================================================
    // EB-R2: Contracted value matches eri_get for (00|00)
    // =========================================================================

    #[test]
    fn eb_r2_contracted_matches_eri_get_0000() {
        let basis = h2_basis();
        let eri = crate::integrals::eri_compressed(&basis);
        let n = basis.n_basis;

        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();
        let ref_val = crate::integrals::eri_get(&eri, n, 0, 0, 0, 0);

        assert_abs_diff_eq!(bd.contracted_value, ref_val, epsilon = 1e-10);
    }

    // =========================================================================
    // EB-R3: Contracted value matches eri_get for (00|11)
    // =========================================================================

    #[test]
    fn eb_r3_contracted_matches_eri_get_0011() {
        let basis = h2_basis();
        let eri = crate::integrals::eri_compressed(&basis);
        let n = basis.n_basis;

        let bd = eri_with_breakdown(&basis, 0, 0, 1, 1).unwrap();
        let ref_val = crate::integrals::eri_get(&eri, n, 0, 0, 1, 1);

        assert_abs_diff_eq!(bd.contracted_value, ref_val, epsilon = 1e-10);
    }

    // =========================================================================
    // EB-R4: Contracted value matches eri_get for (01|01)
    // =========================================================================

    #[test]
    fn eb_r4_contracted_matches_eri_get_0101() {
        let basis = h2_basis();
        let eri = crate::integrals::eri_compressed(&basis);
        let n = basis.n_basis;

        let bd = eri_with_breakdown(&basis, 0, 1, 0, 1).unwrap();
        let ref_val = crate::integrals::eri_get(&eri, n, 0, 1, 0, 1);

        assert_abs_diff_eq!(bd.contracted_value, ref_val, epsilon = 1e-10);
    }

    // =========================================================================
    // EB-R5: Correct number of primitive quartets
    // =========================================================================

    #[test]
    fn eb_r5_sto3g_ssss_gives_81_contributions() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();
        // STO-3G: 3 primitives per s-shell, 3^4 = 81 quartets
        assert_eq!(bd.contributions.len(), 81);
    }

    // =========================================================================
    // EB-R6: Indices match request
    // =========================================================================

    #[test]
    fn eb_r6_indices_match_request() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 1, 1).unwrap();
        assert_eq!(bd.indices, [0, 0, 1, 1]);
    }

    // =========================================================================
    // EB-R7: Labels populated
    // =========================================================================

    #[test]
    fn eb_r7_labels_populated() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 1, 1).unwrap();
        for label in &bd.labels {
            assert!(label.contains("H"), "Label should contain 'H': {}", label);
            assert!(label.contains("1s"), "Label should contain '1s': {}", label);
        }
    }

    // =========================================================================
    // EB-R8: Contributions sorted by magnitude descending
    // =========================================================================

    #[test]
    fn eb_r8_contributions_sorted_descending() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();

        for window in bd.contributions.windows(2) {
            assert!(
                window[0].weighted_contribution.abs() >= window[1].weighted_contribution.abs(),
                "ERI contributions not sorted: {} >= {} failed",
                window[0].weighted_contribution.abs(),
                window[1].weighted_contribution.abs()
            );
        }
    }

    // =========================================================================
    // EB-R9: n_primitives correct
    // =========================================================================

    #[test]
    fn eb_r9_n_primitives_correct() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();
        assert_eq!(bd.n_primitives, [3, 3, 3, 3]);
    }

    // =========================================================================
    // EB-R10: nroots = 1 for (ss|ss)
    // =========================================================================

    #[test]
    fn eb_r10_nroots_1_for_ssss() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();
        assert_eq!(bd.nroots, 1);
        assert_eq!(bd.total_angular_momentum, 0);
    }

    // =========================================================================
    // EB-R11: (ss|ss) uses BoysFunction method
    // =========================================================================

    #[test]
    fn eb_r11_ssss_uses_boys_method() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();
        assert!(
            matches!(bd.method, EriMethod::BoysFunction { .. }),
            "Expected BoysFunction method for (ss|ss)"
        );
    }

    // =========================================================================
    // EB-R12: Method T parameter positive for non-zero separation
    // =========================================================================

    #[test]
    fn eb_r12_t_parameter_positive_for_separation() {
        let basis = h2_basis();
        // (00|11) has non-zero separation between bra and ket centers
        let bd = eri_with_breakdown(&basis, 0, 0, 1, 1).unwrap();
        match &bd.method {
            EriMethod::BoysFunction { t_parameter } => {
                assert!(
                    *t_parameter > 0.0,
                    "T should be positive for non-zero separation"
                );
            }
            _ => panic!("Expected BoysFunction for (ss|ss)"),
        }
    }

    // =========================================================================
    // EB-R13: p-type quartet uses RysQuadrature method (H2O)
    // =========================================================================

    #[test]
    fn eb_r13_p_type_uses_rys_method() {
        let basis = h2o_basis();
        // (00|22) is (ss|pp) with L_total = 2, nroots = 2
        let bd = eri_with_breakdown(&basis, 0, 0, 2, 2).unwrap();
        assert!(
            matches!(bd.method, EriMethod::RysQuadrature { .. }),
            "Expected RysQuadrature method for (ss|pp)"
        );
    }

    // =========================================================================
    // EB-R14: RysQuadrature has correct nroots for (ss|pp)
    // =========================================================================

    #[test]
    fn eb_r14_rys_correct_nroots() {
        let basis = h2o_basis();
        // (00|22) is (ss|pp), L_total = 0+0+1+1 = 2, nroots = 2
        let bd = eri_with_breakdown(&basis, 0, 0, 2, 2).unwrap();
        assert_eq!(bd.nroots, 2);
        if let EriMethod::RysQuadrature { nroots, .. } = &bd.method {
            assert_eq!(*nroots, 2);
        }
    }

    // =========================================================================
    // EB-R15: Rys roots in [0, 1)
    // =========================================================================

    #[test]
    fn eb_r15_rys_roots_bounded() {
        let basis = h2o_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 2, 2).unwrap();
        if let EriMethod::RysQuadrature { roots, .. } = &bd.method {
            for &r in roots {
                assert!(r >= 0.0 && r < 1.0, "Root {} out of [0, 1)", r);
            }
        }
    }

    // =========================================================================
    // EB-R16: Rys weights positive
    // =========================================================================

    #[test]
    fn eb_r16_rys_weights_positive() {
        let basis = h2o_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 2, 2).unwrap();
        if let EriMethod::RysQuadrature { weights, .. } = &bd.method {
            for &w in weights {
                assert!(w > 0.0, "Weight {} should be positive", w);
            }
        }
    }

    // =========================================================================
    // EB-R17 through EB-R21: Cross-system validation
    // =========================================================================

    #[test]
    fn eb_r17_heh_plus_0000_matches_eri_get() {
        let he = Atom::new(2, [0.0, 0.0, 0.0]).unwrap();
        let h = Atom::new(1, [0.0, 0.0, 1.4632]).unwrap();
        let basis = BasisSet::build(vec![he, h], "sto-3g").unwrap();
        let eri = crate::integrals::eri_compressed(&basis);
        let n = basis.n_basis;

        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();
        let ref_val = crate::integrals::eri_get(&eri, n, 0, 0, 0, 0);
        assert_abs_diff_eq!(bd.contracted_value, ref_val, epsilon = 1e-10);
    }

    #[test]
    fn eb_r18_lih_0000_matches_eri_get() {
        let basis = lih_basis();
        let eri = crate::integrals::eri_compressed(&basis);
        let n = basis.n_basis;

        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();
        let ref_val = crate::integrals::eri_get(&eri, n, 0, 0, 0, 0);
        assert_abs_diff_eq!(bd.contracted_value, ref_val, epsilon = 1e-10);
    }

    #[test]
    fn eb_r19_h2o_0055_matches_eri_get() {
        let basis = h2o_basis();
        let eri = crate::integrals::eri_compressed(&basis);
        let n = basis.n_basis;

        let bd = eri_with_breakdown(&basis, 0, 0, 5, 5).unwrap();
        let ref_val = crate::integrals::eri_get(&eri, n, 0, 0, 5, 5);
        assert_abs_diff_eq!(bd.contracted_value, ref_val, epsilon = 1e-10);
    }

    #[test]
    fn eb_r20_h2o_0022_p_type_matches_eri_get() {
        let basis = h2o_basis();
        let eri = crate::integrals::eri_compressed(&basis);
        let n = basis.n_basis;

        // (00|22) is (ss|pp) - involves p-type shell
        let bd = eri_with_breakdown(&basis, 0, 0, 2, 2).unwrap();
        let ref_val = crate::integrals::eri_get(&eri, n, 0, 0, 2, 2);
        assert_abs_diff_eq!(bd.contracted_value, ref_val, epsilon = 1e-10);
    }

    #[test]
    fn eb_r21_h2o_2233_pp_pp_matches_eri_get() {
        let basis = h2o_basis();
        let eri = crate::integrals::eri_compressed(&basis);
        let n = basis.n_basis;

        // (22|33) is (pp|pp) - all p-type
        let bd = eri_with_breakdown(&basis, 2, 2, 3, 3).unwrap();
        let ref_val = crate::integrals::eri_get(&eri, n, 2, 2, 3, 3);
        assert_abs_diff_eq!(bd.contracted_value, ref_val, epsilon = 1e-10);
    }

    // =========================================================================
    // EB-P1 through EB-P5: Property-based tests
    // =========================================================================

    #[test]
    fn eb_p1_all_contributions_finite() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 1, 1).unwrap();
        for c in &bd.contributions {
            assert!(c.weighted_contribution.is_finite());
            assert!(c.primitive_value.is_finite());
        }
    }

    #[test]
    fn eb_p2_all_exponents_positive() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 1, 1).unwrap();
        for c in &bd.contributions {
            for &e in &c.exponents {
                assert!(e > 0.0, "Exponent should be positive: {}", e);
            }
        }
    }

    #[test]
    fn eb_p3_all_t_parameters_non_negative() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 1, 1).unwrap();
        for c in &bd.contributions {
            assert!(
                c.t_parameter >= 0.0,
                "T should be non-negative: {}",
                c.t_parameter
            );
        }
    }

    #[test]
    fn eb_p4_permutation_symmetry_ij_swap() {
        let basis = h2_basis();
        let bd_01_00 = eri_with_breakdown(&basis, 0, 1, 0, 0).unwrap();
        let bd_10_00 = eri_with_breakdown(&basis, 1, 0, 0, 0).unwrap();
        assert_abs_diff_eq!(
            bd_01_00.contracted_value,
            bd_10_00.contracted_value,
            epsilon = 1e-12
        );
    }

    #[test]
    fn eb_p5_permutation_symmetry_bra_ket_swap() {
        let basis = h2_basis();
        let bd_0011 = eri_with_breakdown(&basis, 0, 0, 1, 1).unwrap();
        let bd_1100 = eri_with_breakdown(&basis, 1, 1, 0, 0).unwrap();
        assert_abs_diff_eq!(
            bd_0011.contracted_value,
            bd_1100.contracted_value,
            epsilon = 1e-12
        );
    }

    // =========================================================================
    // EB-E1 through EB-E2: Error handling
    // =========================================================================

    #[test]
    fn eb_e1_out_of_range_index_i() {
        let basis = h2_basis();
        let result = eri_with_breakdown(&basis, 5, 0, 0, 0);
        assert!(result.is_err());
    }

    #[test]
    fn eb_e2_out_of_range_index_l() {
        let basis = h2_basis();
        let result = eri_with_breakdown(&basis, 0, 0, 0, 5);
        assert!(result.is_err());
    }

    // =========================================================================
    // EB-G1 through EB-G5: PySCF golden value tests
    // PySCF 2.11.0, H2 STO-3G at R = 1.3984 bohr
    // =========================================================================

    #[test]
    fn eb_g1_h2_0000_matches_pyscf() {
        // PySCF 2.11.0: (00|00) = 7.746059439198978e-01
        // Note: IQCP Rys quadrature matches PySCF within ~2e-9 (well within
        // the 1e-8 SCF energy convergence threshold). The small difference
        // arises from primitive summation order in the breakdown function.
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();
        assert_abs_diff_eq!(bd.contracted_value, 7.746059439198978e-01, epsilon = 5e-9);
    }

    #[test]
    fn eb_g2_h2_0011_matches_pyscf() {
        // PySCF 2.11.0: (00|11) = 5.699943512349178e-01
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 1, 1).unwrap();
        assert_abs_diff_eq!(bd.contracted_value, 5.699943512349178e-01, epsilon = 5e-9);
    }

    #[test]
    fn eb_g3_h2_0101_matches_pyscf() {
        // PySCF 2.11.0: (01|01) = 2.975896148546511e-01
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 1, 0, 1).unwrap();
        assert_abs_diff_eq!(bd.contracted_value, 2.975896148546511e-01, epsilon = 5e-9);
    }

    #[test]
    fn eb_g4_h2o_0000_matches_pyscf() {
        // PySCF 2.11.0: (00|00) = 4.785065404705503e+00
        let basis = h2o_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 0, 0).unwrap();
        assert_abs_diff_eq!(bd.contracted_value, 4.785065404705503e+00, epsilon = 5e-9);
    }

    #[test]
    fn eb_g5_h2o_0022_matches_pyscf() {
        // PySCF 2.11.0: (00|22) = 1.115813812152427e+00
        let basis = h2o_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 2, 2).unwrap();
        assert_abs_diff_eq!(bd.contracted_value, 1.115813812152427e+00, epsilon = 5e-9);
    }

    // =========================================================================
    // Weighted contribution consistency for ERI
    // =========================================================================

    #[test]
    fn eb_weighted_equals_norm_coeff_product_times_primitive() {
        let basis = h2_basis();
        let bd = eri_with_breakdown(&basis, 0, 0, 1, 1).unwrap();

        for pc in &bd.contributions {
            let expected = pc.norm_coefficients[0]
                * pc.norm_coefficients[1]
                * pc.norm_coefficients[2]
                * pc.norm_coefficients[3]
                * pc.primitive_value;
            assert_abs_diff_eq!(pc.weighted_contribution, expected, epsilon = 1e-15);
        }
    }

    // =========================================================================
    // Sum invariant across all H2 unique quartets
    // =========================================================================

    #[test]
    fn eb_sum_invariant_all_h2_quartets() {
        let basis = h2_basis();
        let eri = crate::integrals::eri_compressed(&basis);
        let n = basis.n_basis;

        // All unique quartets for H2 (2 basis functions)
        let quartets = [
            (0, 0, 0, 0),
            (0, 0, 0, 1),
            (0, 0, 1, 1),
            (0, 1, 0, 1),
            (0, 1, 1, 1),
            (1, 1, 1, 1),
        ];

        for (i, j, k, l) in quartets {
            let bd = eri_with_breakdown(&basis, i, j, k, l).unwrap();
            let sum: f64 = bd
                .contributions
                .iter()
                .map(|c| c.weighted_contribution)
                .sum();
            let ref_val = crate::integrals::eri_get(&eri, n, i, j, k, l);

            assert_abs_diff_eq!(sum, bd.contracted_value, epsilon = 1e-12);
            assert_abs_diff_eq!(bd.contracted_value, ref_val, epsilon = 1e-10);
        }
    }
}
