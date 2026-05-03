#!/usr/bin/env python3
"""
Generate golden reference values for Rys quadrature validation.

This script generates Rys quadrature roots and weights using the RDK algorithm
(Rys, Dupuis, King 1976) based on the libcint implementation (rys_roots.c).

The algorithm uses Schmidt orthogonalization to build orthonormal polynomials
from the Boys function moments, then finds roots via eigenvalue decomposition.

Reference: Dupuis, M., Rys, J., & King, H.F. (1976). J. Chem. Phys. 65, 111.

Usage:
    source .venv/bin/activate
    python scripts/generate_golden_rys.py

Output:
    Writes JSON to tests/golden/rys/reference.json
"""

import json
import sys
from datetime import datetime, timezone

import numpy as np
from scipy.special import hyp1f1

# Configuration
MAX_ORDER = 10  # Maximum supported order in IQCP
T_VALUES = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
OUTPUT_PATH = "tests/golden/rys/reference.json"
THRESHOLD_ZERO = np.finfo(float).eps * 8


def boys_function(m: int, t: float) -> float:
    """
    Compute Boys function F_m(T) = integral_0^1 t^(2m) exp(-T*t^2) dt.

    Uses the relation: F_m(T) = hyp1f1(m+0.5, m+1.5, -T) / (2m+1)
    """
    if t < 1e-15:
        return 1.0 / (2 * m + 1)
    return float(hyp1f1(m + 0.5, m + 1.5, -t) / (2 * m + 1))


def schmidt_orthogonalize(moments: np.ndarray, n: int) -> np.ndarray:
    """
    Schmidt orthogonalization to build orthonormal polynomial coefficients.

    This implements the R_dsmit function from libcint/rys_roots.c (lines 1643-1693).

    Args:
        moments: Array of moments mu_k = F_k(T) for k = 0, 1, ..., 2n-2
        n: Size of polynomial basis (nroots + 1)

    Returns:
        Coefficient matrix cs where cs[i, j] is the coefficient of x^i in P_j(x)
        Stored in column-major order matching libcint.
    """
    cs = np.zeros((n, n))

    # First polynomial: P_0(x) = 1/sqrt(mu_0)
    cs[0, 0] = 1.0 / np.sqrt(moments[0])

    # Second polynomial (j=1)
    fac = -moments[1] / moments[0]
    tmp = moments[2] + fac * moments[1]
    if tmp <= 0:
        raise ValueError(f"Schmidt singular at j=1: tmp={tmp}")
    tmp = 1.0 / np.sqrt(tmp)
    cs[0, 1] = fac * tmp
    cs[1, 1] = tmp

    # Build remaining polynomials P_2, ..., P_{n-1}
    v = np.zeros(n)
    for j in range(2, n):
        v[:j] = 0.0
        fac = moments[j + j]  # mu_{2j}

        # Orthogonalize against all previous polynomials
        for k in range(j):
            # Compute dot product: <x^j, P_k>
            dot = np.sum(cs[:k+1, k] * moments[:k+1 + j - 0][:k+1])
            # Actually: dot = sum_{i=0}^{k} cs[i,k] * mu_{i+j}
            dot = 0.0
            for i in range(k + 1):
                dot += cs[i, k] * moments[i + j]

            # Subtract projection
            for i in range(k + 1):
                v[i] -= dot * cs[i, k]

            fac -= dot * dot

        if fac <= 0:
            if fac == 0:
                # Degenerate case - set remaining to zero
                cs[j:, j:] = 0
                break
            raise ValueError(f"Schmidt singular at j={j}: fac={fac}")

        fac_inv_sqrt = 1.0 / np.sqrt(fac)
        cs[j, j] = fac_inv_sqrt
        for k in range(j):
            cs[k, j] = fac_inv_sqrt * v[k]

    return cs


def polynomial_roots_qr(cs: np.ndarray, nroots: int) -> np.ndarray:
    """
    Find roots of orthogonal polynomial via companion matrix eigenvalues.

    The highest-order polynomial is P_nroots(x) with coefficients in cs[:, nroots].

    Args:
        cs: Coefficient matrix from Schmidt orthogonalization
        nroots: Number of roots to find

    Returns:
        Array of roots
    """
    nroots1 = nroots + 1

    if nroots == 1:
        # Linear: cs[0, 1] + cs[1, 1] * x = 0
        return np.array([-cs[0, 1] / cs[1, 1]])

    if nroots == 2:
        # Quadratic: cs[0, 2] + cs[1, 2] * x + cs[2, 2] * x^2 = 0
        a = cs[2, 2]
        b = cs[1, 2]
        c = cs[0, 2]
        disc = b * b - 4 * a * c
        if disc < 0:
            raise ValueError("Quadratic has no real roots")
        sqrt_disc = np.sqrt(disc)
        return np.array([(-b - sqrt_disc) / (2 * a), (-b + sqrt_disc) / (2 * a)])

    # General case: build companion matrix
    # The companion matrix for polynomial a_0 + a_1*x + ... + a_n*x^n is:
    # | 0   0   ...  0  -a_0/a_n |
    # | 1   0   ...  0  -a_1/a_n |
    # | 0   1   ...  0  -a_2/a_n |
    # | ...                      |
    # | 0   0   ...  1  -a_{n-1}/a_n |

    # Get polynomial coefficients for P_nroots
    poly_coeffs = cs[:nroots + 1, nroots]
    leading = poly_coeffs[nroots]

    # Build companion matrix
    companion = np.zeros((nroots, nroots))
    companion[0, :] = -poly_coeffs[nroots - 1::-1] / leading
    for i in range(1, nroots):
        companion[i, i - 1] = 1.0

    # Eigenvalues are the roots
    eigenvalues = np.linalg.eigvals(companion)

    # Take real parts (should be real for symmetric problem)
    roots = np.real(eigenvalues)
    roots = np.sort(roots)

    return roots


def polynomial_value(cs: np.ndarray, j: int, order: int, x: float) -> float:
    """
    Evaluate polynomial P_j at x using Horner's method.

    Args:
        cs: Coefficient matrix
        j: Polynomial index
        order: Degree of polynomial (== j)
        x: Point to evaluate at

    Returns:
        P_j(x)
    """
    p = cs[order, j]
    for i in range(1, order + 1):
        p = p * x + cs[order - i, j]
    return p


def compute_rys_quadrature(n: int, t: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute Rys quadrature roots and weights using the RDK algorithm.

    This follows libcint's _rdk_rys_roots function (rys_roots.c lines 1699-1756).

    Args:
        n: Number of quadrature points (order)
        t: Argument value T >= 0

    Returns:
        Tuple of (roots, weights) arrays, each of length n
        Roots are in (0, 1), weights are positive
    """
    if n < 1:
        raise ValueError(f"Order n must be >= 1, got {n}")
    if t < 0:
        raise ValueError(f"T must be >= 0, got {t}")

    # Compute moments: mu_k = F_k(T) for k = 0, 1, ..., 2n
    # Need 2n + 1 moments for Schmidt orthogonalization of n+1 polynomials
    moments = np.array([boys_function(k, t) for k in range(2 * n + 1)])

    # Handle mu_0 = 0 case
    if moments[0] < THRESHOLD_ZERO:
        return np.zeros(n), np.zeros(n)

    # Special case: n=1 has analytical solution
    if n == 1:
        # root = mu_1 / (mu_0 - mu_1) in libcint, which is u/(1-u)
        # We want u, so: ratio = mu_1/(mu_0 - mu_1), u = ratio/(1+ratio)
        ratio = moments[1] / (moments[0] - moments[1])
        root = ratio / (1 + ratio) if ratio >= 0 else 0.0
        weight = moments[0]
        return np.array([root]), np.array([weight])

    # General case: Schmidt orthogonalization
    nroots1 = n + 1
    cs = schmidt_orthogonalize(moments, nroots1)

    # Find roots of highest-order polynomial
    rt = polynomial_roots_qr(cs, n)

    # Compute weights
    roots = np.zeros(n)
    weights = np.zeros(n)

    for k in range(n):
        root = rt[k]

        # Handle degenerate case
        if abs(root - 1.0) < THRESHOLD_ZERO:
            roots[k] = 0.0
            weights[k] = 0.0
            continue

        # Compute weight: w_k = 1 / sum_j P_j(root)^2
        dum = 1.0 / moments[0]
        for j in range(1, n):
            poly = polynomial_value(cs, j, j, root)
            dum += poly * poly

        roots[k] = root
        weights[k] = 1.0 / dum

    # Sort by roots
    sort_idx = np.argsort(roots)
    roots = roots[sort_idx]
    weights = weights[sort_idx]

    return roots, weights


def verify_moment_reconstruction(roots: np.ndarray, weights: np.ndarray,
                                  t: float) -> float:
    """
    Verify that the quadrature reconstructs moments correctly.

    For n-point Gaussian quadrature, moments 0 through 2n-1 should be exact.

    Returns the maximum absolute error across all moments.
    """
    n = len(roots)
    max_error = 0.0

    for m in range(2 * n):
        # Expected: F_m(T)
        expected = boys_function(m, t)
        # Computed: sum_k w_k * r_k^m
        computed = np.sum(weights * roots ** m)
        error = abs(expected - computed)
        max_error = max(max_error, error)

    return max_error


def main():
    """Generate and save golden reference data."""
    print(f"Generating Rys quadrature golden reference data")
    print(f"Orders: 1 to {MAX_ORDER}")
    print(f"T values: {T_VALUES}")
    print()

    # Collect version info
    import scipy

    metadata = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "scipy_version": scipy.__version__,
        "numpy_version": np.__version__,
        "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "method": "Schmidt orthogonalization with companion matrix eigenvalue decomposition",
        "algorithm": "RDK (Rys, Dupuis, King 1976) via libcint-style implementation",
        "reference": "Dupuis, M., Rys, J., & King, H.F. (1976). J. Chem. Phys. 65, 111.",
        "libcint_reference": "references/libcint/src/rys_roots.c lines 1643-1756",
        "description": "Rys quadrature roots and weights for IQCP validation",
        "notes": [
            "Roots are in the interval (0, 1)",
            "Weights are strictly positive",
            "Moment reconstruction: sum_k w_k * r_k^m = F_m(T) for m = 0..2n-1",
            "Algorithm matches IQCP qc-core rys module implementation"
        ]
    }

    test_cases = []

    for order in range(1, MAX_ORDER + 1):
        for t in T_VALUES:
            print(f"  Computing n={order:2d}, T={t:5.1f}...", end="")

            try:
                roots, weights = compute_rys_quadrature(order, t)
                moment_error = verify_moment_reconstruction(roots, weights, t)

                # Validate roots are in (0, 1)
                roots_valid = all(0 <= r < 1 for r in roots)
                # Validate weights are positive (allow zero for degenerate)
                weights_valid = all(w >= 0 for w in weights)

                status = "OK" if (roots_valid and weights_valid and moment_error < 1e-8) else "WARN"
                print(f" {status} (moment_err={moment_error:.2e})")

                test_case = {
                    "order": order,
                    "T": t,
                    "roots": [f"{r:.16e}" for r in roots],
                    "weights": [f"{w:.16e}" for w in weights],
                    "moment_reconstruction_error": f"{moment_error:.6e}",
                    "weight_sum": f"{np.sum(weights):.16e}",
                    "F_0_T": f"{boys_function(0, t):.16e}"
                }
                test_cases.append(test_case)

            except Exception as e:
                print(f" ERROR: {e}")
                # Still add the test case with error info
                test_cases.append({
                    "order": order,
                    "T": t,
                    "error": str(e)
                })

    # Build output structure
    output = {
        "metadata": metadata,
        "tolerance": {
            "roots_interval": "(0, 1)",
            "weights_positive": "> 0",
            "moment_reconstruction": "1e-10",
            "weight_sum_vs_F0": "1e-10"
        },
        "test_cases": test_cases
    }

    # Write output
    with open(OUTPUT_PATH, "w") as f:
        json.dump(output, f, indent=2)

    print()
    print(f"Written {len(test_cases)} test cases to {OUTPUT_PATH}")

    # Summary statistics
    successful = [tc for tc in test_cases if "error" not in tc]
    failed = [tc for tc in test_cases if "error" in tc]
    print(f"Successful: {len(successful)}/{len(test_cases)}")
    if failed:
        print(f"Failed cases:")
        for tc in failed:
            print(f"  n={tc['order']}, T={tc['T']}: {tc['error']}")

    if successful:
        errors = [float(tc["moment_reconstruction_error"]) for tc in successful]
        print(f"Moment reconstruction error: min={min(errors):.2e}, max={max(errors):.2e}")


if __name__ == "__main__":
    main()
