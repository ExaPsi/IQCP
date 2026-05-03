/**
 * Gaussian primitive normalization utility.
 *
 * Computes the normalization constant N for a Cartesian Gaussian primitive
 * using a representative component for each angular momentum type:
 *
 *   s (l=0): N = (2*alpha/pi)^(3/4)
 *   p (l=1): N = (2*alpha/pi)^(3/4) * (4*alpha)^(1/2)
 *   d (l=2): N = (2*alpha/pi)^(3/4) * 4*alpha        (d_xy representative)
 *
 * This is a pure TypeScript port of `primitive_normalization` from
 * `crates/qc-core/src/basis/radial.rs`, using identical formulas.
 *
 * Reference: Szabo & Ostlund (1996), Eq. 3.203.
 *
 * @module components/basisExplorer/normalization
 */

/**
 * Compute the normalization constant for a Gaussian primitive.
 *
 * Uses a representative Cartesian component for each angular momentum:
 * - s (l=0): (0,0,0) -> N = (2*alpha/pi)^(3/4)
 * - p (l=1): (0,0,1) -> N = (2*alpha/pi)^(3/4) * (4*alpha)^(1/2)
 * - d (l=2): (1,1,0) -> N = (2*alpha/pi)^(3/4) * 4*alpha
 *
 * Consistent with `primitive_normalization` in `crates/qc-core/src/basis/radial.rs`.
 *
 * @param alpha - Primitive Gaussian exponent (must be positive)
 * @param l - Angular momentum quantum number (0=s, 1=p, 2=d)
 * @returns Normalization constant N_k
 */
export function primitiveNormalization(alpha: number, l: number): number {
  // Base normalization: (2*alpha/pi)^(3/4)
  const base = Math.pow((2.0 * alpha) / Math.PI, 0.75);

  switch (l) {
    case 0:
      // s-type: N = (2*alpha/pi)^(3/4)
      // Representative Cartesian: (0,0,0), no angular factor
      return base;
    case 1:
      // p-type: N = (2*alpha/pi)^(3/4) * (4*alpha)^(1/2)
      // Representative Cartesian: (0,0,1)
      // Double factorials: (-1)!!=1, (-1)!!=1, (1)!!=1, product=1
      return base * Math.sqrt(4.0 * alpha);
    case 2:
      // d-type: N = (2*alpha/pi)^(3/4) * 4*alpha
      // Representative Cartesian: (1,1,0) (d_xy)
      // Powers (1,1,0): (2*1-1)!!=1, (2*1-1)!!=1, (-1)!!=1
      // Angular factor: (4*alpha)^(2/2) / sqrt(1*1*1) = 4*alpha
      return base * 4.0 * alpha;
    default:
      // Unsupported angular momentum; fall back to base
      return base;
  }
}
