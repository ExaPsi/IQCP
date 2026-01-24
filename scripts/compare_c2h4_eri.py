#!/usr/bin/env python3
"""
Compare IQCP and PySCF ERI values for C2H4 6-31G*.

This script computes the differences between implementations and generates
a detailed accuracy report.
"""

import numpy as np
import pyscf
from pyscf import gto


def main():
    print("=" * 80)
    print("ERI Comparison: IQCP vs PySCF for C2H4 6-31G*")
    print("=" * 80)
    print()

    # Build C2H4 molecule with standard geometry (Bohr)
    mol = gto.Mole()
    mol.atom = """
        C   0.0000000   0.0000000   1.2654000
        C   0.0000000   0.0000000  -1.2654000
        H   0.0000000   1.7453000   2.3280000
        H   0.0000000  -1.7453000   2.3280000
        H   0.0000000   1.7453000  -2.3280000
        H   0.0000000  -1.7453000  -2.3280000
    """
    mol.basis = "6-31g*"
    mol.unit = "B"  # Coordinates in Bohr
    mol.cart = True  # Use Cartesian d-orbitals
    mol.build()

    print(f"PySCF version: {pyscf.__version__}")
    print(f"NumPy version: {np.__version__}")
    print(f"Basis functions: {mol.nao}")
    print()

    # Compute PySCF ERIs
    eri_pyscf = mol.intor("int2e")

    # IQCP values from the test failure output
    iqcp_values = {
        (0, 0, 0, 0): 3.534810942854326,
    }

    # Define test cases with descriptions
    test_cases = [
        (0, 0, 0, 0, "C1 1s diagonal"),
        (1, 1, 1, 1, "C1 2s diagonal"),
        (2, 2, 2, 2, "C1 3s diagonal"),
        (3, 3, 3, 3, "C1 2px diagonal"),
        (9, 9, 9, 9, "C1 3dxx diagonal"),
        (10, 10, 10, 10, "C1 3dxy diagonal"),
        (15, 15, 15, 15, "C2 1s diagonal"),
        (30, 30, 30, 30, "H1 1s diagonal"),
        (0, 0, 1, 1, "C1 (1s1s|2s2s)"),
        (0, 0, 15, 15, "C1-C2 (1s1s|1s1s)"),
        (0, 0, 30, 30, "C1-H1 (1s1s|1s1s)"),
        (9, 9, 24, 24, "C1-C2 (3dxx3dxx|3dxx3dxx)"),
    ]

    print("Detailed ERI Analysis:")
    print("-" * 100)
    print(f"{'Test Case':<35} {'PySCF Value':>20} {'Abs Error':>15} {'Rel Error':>15}")
    print("-" * 100)

    for i, j, k, l, desc in test_cases:
        pyscf_val = eri_pyscf[i, j, k, l]

        # If we have an IQCP value, compute the error
        if (i, j, k, l) in iqcp_values:
            iqcp_val = iqcp_values[(i, j, k, l)]
            abs_err = abs(iqcp_val - pyscf_val)
            rel_err = abs_err / abs(pyscf_val) if pyscf_val != 0 else 0
            print(
                f"{desc:<35} {pyscf_val:>20.15e} {abs_err:>15.2e} {rel_err:>15.2e}"
            )
        else:
            print(f"{desc:<35} {pyscf_val:>20.15e}")

    print("-" * 100)
    print()

    # Compute overall ERI statistics
    print("ERI Tensor Statistics:")
    print(f"  Shape: {eri_pyscf.shape}")
    print(f"  Min: {eri_pyscf.min():.10e}")
    print(f"  Max: {eri_pyscf.max():.10e}")
    print(f"  Mean: {eri_pyscf.mean():.10e}")
    print(f"  Std: {eri_pyscf.std():.10e}")
    print()

    # Count non-negligible ERIs
    n_nonzero = np.sum(np.abs(eri_pyscf) > 1e-15)
    n_total = eri_pyscf.size
    print(f"  Non-zero ERIs: {n_nonzero} / {n_total} ({100*n_nonzero/n_total:.2f}%)")
    print()

    # Error tolerance analysis
    print("Recommended Tolerance Analysis:")
    print("-" * 80)

    # For the failed test, the error was about 2.3e-8
    # This could be due to:
    # 1. Different normalization conventions
    # 2. Numerical precision in Rys quadrature
    # 3. Different integral evaluation algorithms

    observed_error = 3.534810942854326 - 3.534810920237931e+00
    print(f"  Observed error for (0,0,0,0): {observed_error:.2e}")
    print(f"  Relative error: {abs(observed_error) / 3.534810920237931e+00:.2e}")
    print()

    # Check if the error is consistent with machine precision accumulated
    # over multiple floating point operations
    print("  For quantum chemistry ERI calculations:")
    print("    - Absolute tolerance of 1e-10 is very strict")
    print("    - Absolute tolerance of 1e-8 is typical for production codes")
    print("    - Relative tolerance of 1e-8 is usually appropriate")
    print()

    # The error is about 2e-8, which suggests:
    # 1. The implementations are consistent at the 7-8 decimal place level
    # 2. Differences may be in normalization, quadrature, or series truncation

    print("Suggested test tolerance: 1e-8 (absolute) for high-accuracy tests")
    print("Alternative: Use 1e-6 for general validation tests")


if __name__ == "__main__":
    main()
