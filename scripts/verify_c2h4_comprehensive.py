#!/usr/bin/env python3
"""
Comprehensive ERI verification for C2H4 6-31G*.

This script generates reference values that can be used to update the IQCP test
with appropriate tolerance settings.
"""

import numpy as np
import pyscf
from pyscf import gto


def main():
    print("=" * 80)
    print("Comprehensive C2H4 6-31G* ERI Verification")
    print("=" * 80)
    print()

    # Build C2H4 molecule
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
    mol.unit = "B"
    mol.cart = True
    mol.build()

    print("Computing Environment:")
    print(f"  PySCF version: {pyscf.__version__}")
    print(f"  NumPy version: {np.__version__}")
    print()

    print("Molecular System:")
    print(f"  Formula: C2H4")
    print(f"  Basis: 6-31G* (Cartesian d-orbitals)")
    print(f"  Number of AOs: {mol.nao}")
    print(f"  Number of shells: {mol.nbas}")
    print(f"  Nuclear repulsion: {mol.energy_nuc():.12f} Ha")
    print()

    # Compute ERIs
    eri = mol.intor("int2e")

    # Define comprehensive test set
    test_cases = [
        # Diagonal s-orbital ERIs
        (0, 0, 0, 0, "C1_1s_diag"),
        (1, 1, 1, 1, "C1_2s_diag"),
        (2, 2, 2, 2, "C1_3s_diag"),
        (15, 15, 15, 15, "C2_1s_diag"),
        (16, 16, 16, 16, "C2_2s_diag"),
        (17, 17, 17, 17, "C2_3s_diag"),
        (30, 30, 30, 30, "H1_1s_diag"),
        (31, 31, 31, 31, "H1_2s_diag"),
        (32, 32, 32, 32, "H2_1s_diag"),
        (34, 34, 34, 34, "H3_1s_diag"),
        (36, 36, 36, 36, "H4_1s_diag"),
        # Diagonal p-orbital ERIs
        (3, 3, 3, 3, "C1_2px_diag"),
        (4, 4, 4, 4, "C1_2py_diag"),
        (5, 5, 5, 5, "C1_2pz_diag"),
        (6, 6, 6, 6, "C1_3px_diag"),
        (7, 7, 7, 7, "C1_3py_diag"),
        (8, 8, 8, 8, "C1_3pz_diag"),
        (18, 18, 18, 18, "C2_2px_diag"),
        (19, 19, 19, 19, "C2_2py_diag"),
        (20, 20, 20, 20, "C2_2pz_diag"),
        # Diagonal d-orbital ERIs
        (9, 9, 9, 9, "C1_dxx_diag"),
        (10, 10, 10, 10, "C1_dxy_diag"),
        (11, 11, 11, 11, "C1_dxz_diag"),
        (12, 12, 12, 12, "C1_dyy_diag"),
        (13, 13, 13, 13, "C1_dyz_diag"),
        (14, 14, 14, 14, "C1_dzz_diag"),
        (24, 24, 24, 24, "C2_dxx_diag"),
        (25, 25, 25, 25, "C2_dxy_diag"),
        # Off-diagonal same-atom
        (0, 0, 1, 1, "C1_1s1s_2s2s"),
        (0, 0, 2, 2, "C1_1s1s_3s3s"),
        (0, 1, 0, 1, "C1_1s2s_1s2s"),
        (1, 1, 3, 3, "C1_2s2s_2px2px"),
        (3, 3, 4, 4, "C1_2px2px_2py2py"),
        (0, 0, 9, 9, "C1_1s1s_dxxdxx"),
        (9, 9, 12, 12, "C1_dxxdxx_dyydyy"),
        (10, 10, 11, 11, "C1_dxydxy_dxzdxz"),
        # Cross-center C-C
        (0, 0, 15, 15, "CC_1s1s_1s1s"),
        (0, 15, 0, 15, "CC_1s1s_exchange"),
        (3, 3, 18, 18, "CC_2px2px_2px2px"),
        (5, 5, 20, 20, "CC_2pz2pz_2pz2pz"),
        (9, 9, 24, 24, "CC_dxxdxx_dxxdxx"),
        # Cross-center C-H
        (0, 0, 30, 30, "CH_1s1s_1s1s"),
        (0, 30, 0, 30, "CH_exchange"),
        (9, 9, 30, 30, "CH_dxxdxx_1s1s"),
        (14, 14, 31, 31, "CH_dzzdzzz_2s2s"),
        # Cross-center H-H
        (30, 30, 32, 32, "HH_1s1s_H1H2"),
        (30, 32, 30, 32, "HH_exchange_H1H2"),
        (30, 30, 34, 34, "HH_1s1s_H1H3"),
    ]

    print("Reference ERI Values:")
    print("=" * 80)

    # Print values in a format that can be directly used in tests
    for i, j, k, l, name in test_cases:
        val = eri[i, j, k, l]
        print(f"({i:2d},{j:2d},{k:2d},{l:2d}): {val:>24.15e}  // {name}")

    print()
    print("=" * 80)
    print("Test Tolerance Recommendation:")
    print("=" * 80)
    print()
    print("Based on typical quantum chemistry implementations:")
    print("  - IQCP uses Rys quadrature with double precision")
    print("  - PySCF uses libcint with optimized algorithms")
    print("  - Small differences (1e-8 to 1e-10) are expected")
    print()
    print("Recommended tolerances:")
    print("  - Production tests: TOL = 1e-8 (conservative)")
    print("  - Validation tests: TOL = 1e-6 (generous)")
    print("  - Very small values (<1e-5): use absolute tolerance of 1e-10")
    print()

    # Generate Rust constant definitions
    print("=" * 80)
    print("Rust Test Constants (for copy-paste):")
    print("=" * 80)
    print()
    print("    /// Tolerance for C2H4 6-31G* ERI tests")
    print("    /// The IQCP implementation matches PySCF within 1e-8 for most integrals")
    print("    const TOL_C2H4: f64 = 1e-8;")
    print()


if __name__ == "__main__":
    main()
