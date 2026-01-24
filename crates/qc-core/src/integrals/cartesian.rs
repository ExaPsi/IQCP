//! Cartesian angular momentum components
//!
//! Manages the Cartesian powers (i, j, k) for Gaussian basis functions:
//! ```text
//! χ(r) = x^i · y^j · z^k · exp(-α|r-A|²)
//! ```
//!
//! # Angular Momentum Ordering
//!
//! Cartesian components are ordered lexicographically for each L:
//!
//! | L | Components | Powers (i,j,k) |
//! |---|------------|----------------|
//! | 0 | s          | (0,0,0)        |
//! | 1 | px,py,pz   | (1,0,0), (0,1,0), (0,0,1) |
//! | 2 | dxx,dxy,dxz,dyy,dyz,dzz | (2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2) |

use super::{IntegralError, MAX_ANGULAR_MOMENTUM};

/// Cartesian angular momentum powers (i, j, k) where L = i + j + k
///
/// Represents the polynomial prefactor x^i · y^j · z^k in a Cartesian Gaussian.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CartesianPower {
    /// Power of x component
    pub i: u32,
    /// Power of y component
    pub j: u32,
    /// Power of z component
    pub k: u32,
}

impl CartesianPower {
    /// Create a new CartesianPower
    pub const fn new(i: u32, j: u32, k: u32) -> Self {
        Self { i, j, k }
    }

    /// Total angular momentum L = i + j + k
    #[inline]
    pub const fn angular_momentum(&self) -> u32 {
        self.i + self.j + self.k
    }

    /// Get power for a specific axis (0=x, 1=y, 2=z)
    #[inline]
    pub const fn power(&self, axis: usize) -> u32 {
        match axis {
            0 => self.i,
            1 => self.j,
            2 => self.k,
            _ => 0,
        }
    }

    /// Create a new CartesianPower with one unit added to the specified axis
    #[inline]
    pub const fn increment(&self, axis: usize) -> Self {
        match axis {
            0 => Self::new(self.i + 1, self.j, self.k),
            1 => Self::new(self.i, self.j + 1, self.k),
            2 => Self::new(self.i, self.j, self.k + 1),
            _ => *self,
        }
    }

    /// Create a new CartesianPower with one unit removed from the specified axis
    /// Returns None if the power would become negative
    #[inline]
    pub const fn decrement(&self, axis: usize) -> Option<Self> {
        match axis {
            0 if self.i > 0 => Some(Self::new(self.i - 1, self.j, self.k)),
            1 if self.j > 0 => Some(Self::new(self.i, self.j - 1, self.k)),
            2 if self.k > 0 => Some(Self::new(self.i, self.j, self.k - 1)),
            _ => None,
        }
    }
}

/// Generate all Cartesian components for angular momentum L
///
/// Returns components in lexicographic order of (i, j, k) powers.
///
/// # Arguments
///
/// * `l` - Total angular momentum quantum number
///
/// # Returns
///
/// Vector of CartesianPower in standard ordering
///
/// # Errors
///
/// Returns error if `l > MAX_ANGULAR_MOMENTUM`
///
/// # Example
///
/// ```rust
/// use qc_core::integrals::cartesian_components;
///
/// // S orbital: one component
/// let s = cartesian_components(0).unwrap();
/// assert_eq!(s.len(), 1);
///
/// // P orbitals: three components (px, py, pz)
/// let p = cartesian_components(1).unwrap();
/// assert_eq!(p.len(), 3);
///
/// // D orbitals: six components
/// let d = cartesian_components(2).unwrap();
/// assert_eq!(d.len(), 6);
/// ```
pub fn cartesian_components(l: u32) -> Result<Vec<CartesianPower>, IntegralError> {
    if l > MAX_ANGULAR_MOMENTUM {
        return Err(IntegralError::AngularMomentumTooHigh(
            l,
            MAX_ANGULAR_MOMENTUM,
        ));
    }

    let mut components = Vec::with_capacity(n_cartesian(l) as usize);

    // Generate in lexicographic order: i decreases, j+k increases
    // For each i from l down to 0, j goes from l-i down to 0
    for i in (0..=l).rev() {
        for j in (0..=(l - i)).rev() {
            let k = l - i - j;
            components.push(CartesianPower::new(i, j, k));
        }
    }

    Ok(components)
}

/// Number of Cartesian components for angular momentum L
///
/// Formula: (L+1)(L+2)/2
#[inline]
pub const fn n_cartesian(l: u32) -> u32 {
    (l + 1) * (l + 2) / 2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cartesian_power_new() {
        let p = CartesianPower::new(1, 2, 3);
        assert_eq!(p.i, 1);
        assert_eq!(p.j, 2);
        assert_eq!(p.k, 3);
    }

    #[test]
    fn cartesian_power_angular_momentum() {
        assert_eq!(CartesianPower::new(0, 0, 0).angular_momentum(), 0); // s
        assert_eq!(CartesianPower::new(1, 0, 0).angular_momentum(), 1); // px
        assert_eq!(CartesianPower::new(2, 0, 0).angular_momentum(), 2); // dxx
        assert_eq!(CartesianPower::new(1, 1, 0).angular_momentum(), 2); // dxy
    }

    #[test]
    fn cartesian_power_increment() {
        let s = CartesianPower::new(0, 0, 0);
        assert_eq!(s.increment(0), CartesianPower::new(1, 0, 0)); // px
        assert_eq!(s.increment(1), CartesianPower::new(0, 1, 0)); // py
        assert_eq!(s.increment(2), CartesianPower::new(0, 0, 1)); // pz
    }

    #[test]
    fn cartesian_power_decrement() {
        let p = CartesianPower::new(1, 1, 1);
        assert_eq!(p.decrement(0), Some(CartesianPower::new(0, 1, 1)));
        assert_eq!(p.decrement(1), Some(CartesianPower::new(1, 0, 1)));
        assert_eq!(p.decrement(2), Some(CartesianPower::new(1, 1, 0)));

        let s = CartesianPower::new(0, 0, 0);
        assert_eq!(s.decrement(0), None);
        assert_eq!(s.decrement(1), None);
        assert_eq!(s.decrement(2), None);
    }

    #[test]
    fn n_cartesian_values() {
        assert_eq!(n_cartesian(0), 1); // s
        assert_eq!(n_cartesian(1), 3); // p
        assert_eq!(n_cartesian(2), 6); // d
        assert_eq!(n_cartesian(3), 10); // f (not currently supported but formula is general)
    }

    #[test]
    fn cartesian_components_s_orbital() {
        let comps = cartesian_components(0).unwrap();
        assert_eq!(comps.len(), 1);
        assert_eq!(comps[0], CartesianPower::new(0, 0, 0));
    }

    #[test]
    fn cartesian_components_p_orbitals() {
        let comps = cartesian_components(1).unwrap();
        assert_eq!(comps.len(), 3);
        // Lexicographic order: px, py, pz
        assert_eq!(comps[0], CartesianPower::new(1, 0, 0)); // px
        assert_eq!(comps[1], CartesianPower::new(0, 1, 0)); // py
        assert_eq!(comps[2], CartesianPower::new(0, 0, 1)); // pz
    }

    #[test]
    fn cartesian_components_d_orbitals() {
        let comps = cartesian_components(2).unwrap();
        assert_eq!(comps.len(), 6);
        // Lexicographic order
        assert_eq!(comps[0], CartesianPower::new(2, 0, 0)); // dxx
        assert_eq!(comps[1], CartesianPower::new(1, 1, 0)); // dxy
        assert_eq!(comps[2], CartesianPower::new(1, 0, 1)); // dxz
        assert_eq!(comps[3], CartesianPower::new(0, 2, 0)); // dyy
        assert_eq!(comps[4], CartesianPower::new(0, 1, 1)); // dyz
        assert_eq!(comps[5], CartesianPower::new(0, 0, 2)); // dzz
    }

    #[test]
    fn cartesian_components_angular_momentum_too_high() {
        let result = cartesian_components(3);
        assert!(matches!(
            result,
            Err(IntegralError::AngularMomentumTooHigh(3, 2))
        ));
    }

    #[test]
    fn cartesian_components_all_have_correct_angular_momentum() {
        for l in 0..=MAX_ANGULAR_MOMENTUM {
            let comps = cartesian_components(l).unwrap();
            for comp in comps {
                assert_eq!(comp.angular_momentum(), l);
            }
        }
    }
}
