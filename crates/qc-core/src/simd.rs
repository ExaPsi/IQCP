//! SIMD-optimized numerical operations
//!
//! This module provides SIMD-accelerated versions of common numerical operations
//! used throughout the quantum chemistry algorithms. The implementations use the
//! `wide` crate which provides portable SIMD that works on both native and WASM.
//!
//! # Feature Gating
//!
//! All SIMD functionality is gated behind the `simd` feature. When disabled,
//! scalar fallbacks are used automatically.
//!
//! # Performance Notes
//!
//! SIMD provides significant speedups for:
//! - Dot products and vector operations (2-4x)
//! - Matrix element-wise operations (2-4x)
//! - Reduction operations (sums, norms)
//!
//! The `wide` crate uses `f64x4` (4 doubles per vector) on most platforms.
//!
//! # Example
//!
//! ```ignore
//! use qc_core::simd::dot_product;
//!
//! let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
//! let b = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
//! let result = dot_product(&a, &b);
//! assert!((result - 36.0).abs() < 1e-12);
//! ```

#[cfg(feature = "simd")]
use wide::f64x4;

/// SIMD lane width for f64 (4 elements per vector)
#[cfg(feature = "simd")]
pub const SIMD_WIDTH: usize = 4;

/// Compute dot product of two slices using SIMD when available.
///
/// This is a key primitive used in:
/// - Density matrix trace operations
/// - Fock matrix accumulation
/// - Error vector computation for DIIS
///
/// # Arguments
///
/// * `a` - First input slice
/// * `b` - Second input slice (must have same length as `a`)
///
/// # Returns
///
/// The dot product sum(a[i] * b[i])
///
/// # Panics
///
/// Panics if slices have different lengths.
#[cfg(feature = "simd")]
#[inline]
pub fn dot_product(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "Slices must have equal length");

    let n = a.len();
    let chunks = n / SIMD_WIDTH;

    // Process full SIMD chunks
    let mut sum_vec = f64x4::ZERO;
    for i in 0..chunks {
        let offset = i * SIMD_WIDTH;
        let va = f64x4::new([a[offset], a[offset + 1], a[offset + 2], a[offset + 3]]);
        let vb = f64x4::new([b[offset], b[offset + 1], b[offset + 2], b[offset + 3]]);
        sum_vec += va * vb;
    }

    // Reduce SIMD vector to scalar
    let arr: [f64; 4] = sum_vec.into();
    let mut result = arr[0] + arr[1] + arr[2] + arr[3];

    // Handle remainder elements
    for i in (chunks * SIMD_WIDTH)..n {
        result += a[i] * b[i];
    }

    result
}

/// Scalar fallback for dot product
#[cfg(not(feature = "simd"))]
#[inline]
pub fn dot_product(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "Slices must have equal length");
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Compute sum of squares of elements using SIMD when available.
///
/// Used for:
/// - RMS calculations
/// - Norm computations
/// - Convergence checks
///
/// # Arguments
///
/// * `a` - Input slice
///
/// # Returns
///
/// The sum of squared elements sum(a[i]^2)
#[cfg(feature = "simd")]
#[inline]
pub fn sum_of_squares(a: &[f64]) -> f64 {
    let n = a.len();
    let chunks = n / SIMD_WIDTH;

    let mut sum_vec = f64x4::ZERO;
    for i in 0..chunks {
        let offset = i * SIMD_WIDTH;
        let va = f64x4::new([a[offset], a[offset + 1], a[offset + 2], a[offset + 3]]);
        sum_vec += va * va;
    }

    let arr: [f64; 4] = sum_vec.into();
    let mut result = arr[0] + arr[1] + arr[2] + arr[3];

    // Handle remainder using iterator
    for &ai in a.iter().skip(chunks * SIMD_WIDTH) {
        result += ai * ai;
    }

    result
}

/// Scalar fallback for sum of squares
#[cfg(not(feature = "simd"))]
#[inline]
pub fn sum_of_squares(a: &[f64]) -> f64 {
    a.iter().map(|x| x * x).sum()
}

/// Add scaled vector to another: y += alpha * x
///
/// This is the AXPY operation (alpha * x plus y), a BLAS Level 1 routine.
///
/// Used for:
/// - Fock matrix accumulation
/// - DIIS extrapolation
/// - Density matrix updates
///
/// # Arguments
///
/// * `y` - Output slice (modified in place)
/// * `alpha` - Scalar multiplier
/// * `x` - Input slice (must have same length as `y`)
///
/// # Panics
///
/// Panics if slices have different lengths.
#[cfg(feature = "simd")]
#[inline]
pub fn axpy(y: &mut [f64], alpha: f64, x: &[f64]) {
    assert_eq!(y.len(), x.len(), "Slices must have equal length");

    let n = y.len();
    let chunks = n / SIMD_WIDTH;
    let alpha_vec = f64x4::splat(alpha);

    for i in 0..chunks {
        let offset = i * SIMD_WIDTH;
        let vx = f64x4::new([x[offset], x[offset + 1], x[offset + 2], x[offset + 3]]);
        let vy = f64x4::new([y[offset], y[offset + 1], y[offset + 2], y[offset + 3]]);
        let result = vy + alpha_vec * vx;
        let arr: [f64; 4] = result.into();
        y[offset] = arr[0];
        y[offset + 1] = arr[1];
        y[offset + 2] = arr[2];
        y[offset + 3] = arr[3];
    }

    // Handle remainder
    for i in (chunks * SIMD_WIDTH)..n {
        y[i] += alpha * x[i];
    }
}

/// Scalar fallback for axpy
#[cfg(not(feature = "simd"))]
#[inline]
pub fn axpy(y: &mut [f64], alpha: f64, x: &[f64]) {
    assert_eq!(y.len(), x.len(), "Slices must have equal length");
    for i in 0..y.len() {
        y[i] += alpha * x[i];
    }
}

/// Compute element-wise sum of a slice.
///
/// # Arguments
///
/// * `a` - Input slice
///
/// # Returns
///
/// The sum of all elements
#[cfg(feature = "simd")]
#[inline]
pub fn simd_sum(a: &[f64]) -> f64 {
    let n = a.len();
    let chunks = n / SIMD_WIDTH;

    let mut sum_vec = f64x4::ZERO;
    for i in 0..chunks {
        let offset = i * SIMD_WIDTH;
        let va = f64x4::new([a[offset], a[offset + 1], a[offset + 2], a[offset + 3]]);
        sum_vec += va;
    }

    let arr: [f64; 4] = sum_vec.into();
    let mut result = arr[0] + arr[1] + arr[2] + arr[3];

    // Handle remainder using iterator
    for &ai in a.iter().skip(chunks * SIMD_WIDTH) {
        result += ai;
    }

    result
}

/// Scalar fallback for simd_sum
#[cfg(not(feature = "simd"))]
#[inline]
pub fn simd_sum(a: &[f64]) -> f64 {
    a.iter().sum()
}

/// Scale a vector in place: x *= alpha
///
/// # Arguments
///
/// * `x` - Vector to scale (modified in place)
/// * `alpha` - Scalar multiplier
#[cfg(feature = "simd")]
#[inline]
pub fn scale(x: &mut [f64], alpha: f64) {
    let n = x.len();
    let chunks = n / SIMD_WIDTH;
    let alpha_vec = f64x4::splat(alpha);

    for i in 0..chunks {
        let offset = i * SIMD_WIDTH;
        let vx = f64x4::new([x[offset], x[offset + 1], x[offset + 2], x[offset + 3]]);
        let result = vx * alpha_vec;
        let arr: [f64; 4] = result.into();
        x[offset] = arr[0];
        x[offset + 1] = arr[1];
        x[offset + 2] = arr[2];
        x[offset + 3] = arr[3];
    }

    // Handle remainder using iterator
    for xi in x.iter_mut().skip(chunks * SIMD_WIDTH) {
        *xi *= alpha;
    }
}

/// Scalar fallback for scale
#[cfg(not(feature = "simd"))]
#[inline]
pub fn scale(x: &mut [f64], alpha: f64) {
    for xi in x.iter_mut() {
        *xi *= alpha;
    }
}

/// Compute the difference of two vectors and return sum of squared differences.
///
/// This is more efficient than computing the difference and then sum of squares
/// separately.
///
/// Used for:
/// - Density matrix RMS change
/// - Convergence criteria
///
/// # Arguments
///
/// * `a` - First vector
/// * `b` - Second vector (must have same length as `a`)
///
/// # Returns
///
/// sum((a[i] - b[i])^2)
#[cfg(feature = "simd")]
#[inline]
pub fn sum_squared_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "Slices must have equal length");

    let n = a.len();
    let chunks = n / SIMD_WIDTH;

    let mut sum_vec = f64x4::ZERO;
    for i in 0..chunks {
        let offset = i * SIMD_WIDTH;
        let va = f64x4::new([a[offset], a[offset + 1], a[offset + 2], a[offset + 3]]);
        let vb = f64x4::new([b[offset], b[offset + 1], b[offset + 2], b[offset + 3]]);
        let diff = va - vb;
        sum_vec += diff * diff;
    }

    let arr: [f64; 4] = sum_vec.into();
    let mut result = arr[0] + arr[1] + arr[2] + arr[3];

    for i in (chunks * SIMD_WIDTH)..n {
        let diff = a[i] - b[i];
        result += diff * diff;
    }

    result
}

/// Scalar fallback for sum_squared_diff
#[cfg(not(feature = "simd"))]
#[inline]
pub fn sum_squared_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "Slices must have equal length");
    a.iter()
        .zip(b.iter())
        .map(|(ai, bi)| {
            let d = ai - bi;
            d * d
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-14;

    #[test]
    fn test_dot_product_exact() {
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![1.0, 1.0, 1.0, 1.0];
        assert!((dot_product(&a, &b) - 10.0).abs() < TOL);
    }

    #[test]
    fn test_dot_product_with_remainder() {
        // 9 elements to test remainder handling
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let b = vec![1.0; 9];
        assert!((dot_product(&a, &b) - 45.0).abs() < TOL);
    }

    #[test]
    fn test_dot_product_empty() {
        let a: Vec<f64> = vec![];
        let b: Vec<f64> = vec![];
        assert!((dot_product(&a, &b) - 0.0).abs() < TOL);
    }

    #[test]
    fn test_sum_of_squares() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((sum_of_squares(&a) - 55.0).abs() < TOL);
    }

    #[test]
    fn test_axpy() {
        let mut y = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let x = vec![1.0; 8];
        axpy(&mut y, 2.0, &x);
        let expected = vec![3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        for (yi, ei) in y.iter().zip(expected.iter()) {
            assert!((yi - ei).abs() < TOL);
        }
    }

    #[test]
    fn test_simd_sum() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        assert!((simd_sum(&a) - 28.0).abs() < TOL);
    }

    #[test]
    fn test_scale() {
        let mut x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        scale(&mut x, 2.0);
        let expected = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        for (xi, ei) in x.iter().zip(expected.iter()) {
            assert!((xi - ei).abs() < TOL);
        }
    }

    #[test]
    fn test_sum_squared_diff() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let b = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        // Each diff is 1.0, so sum of squares is 5.0
        assert!((sum_squared_diff(&a, &b) - 5.0).abs() < TOL);
    }

    #[test]
    fn test_sum_squared_diff_zero() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 2.0, 3.0];
        assert!((sum_squared_diff(&a, &b) - 0.0).abs() < TOL);
    }
}
