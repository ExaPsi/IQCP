#!/usr/bin/env python3
"""
Generate PySCF reference ERI values for C2H4 (ethylene) with 6-31G* basis set.

This script computes two-electron repulsion integrals using PySCF as the
authoritative reference for validating IQCP's ERI implementation.

C2H4 standard geometry (Bohr):
  C   0.0000   0.0000   1.2654
  C   0.0000   0.0000  -1.2654
  H   0.0000   1.7453   2.3280
  H   0.0000  -1.7453   2.3280
  H   0.0000   1.7453  -2.3280
  H   0.0000  -1.7453  -2.3280

Expected basis functions for C2H4 6-31G* (Cartesian d-orbitals):
- Each C: 1s, 2s, 3s, 2px, 2py, 2pz, 3px, 3py, 3pz, 3dxx, 3dxy, 3dxz, 3dyy, 3dyz, 3dzz (15 functions)
- Each H: 1s, 2s (2 functions)
- Total: 2*15 + 4*2 = 38 basis functions
"""

import numpy as np
import pyscf
from pyscf import gto


def main():
    print("=" * 80)
    print("PySCF Reference ERI Generation for C2H4 6-31G*")
    print("=" * 80)
    print()

    # Print environment info
    print("Computing Environment:")
    print(f"  PySCF version: {pyscf.__version__}")
    print(f"  NumPy version: {np.__version__}")
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
    mol.cart = True  # Use Cartesian d-orbitals (6 components: xx, xy, xz, yy, yz, zz)
    mol.build()

    print("Molecular System:")
    print(f"  Formula: C2H4 (ethylene)")
    print(f"  Basis set: 6-31G* (Cartesian d-orbitals)")
    print(f"  Number of atoms: {mol.natm}")
    print(f"  Number of electrons: {mol.nelectron}")
    print(f"  Number of basis functions (nao): {mol.nao}")
    print(f"  Number of shells (nbas): {mol.nbas}")
    print(f"  Nuclear repulsion energy: {mol.energy_nuc():.12f} Ha")
    print()

    # Print shell information
    print("Shell Information:")
    print("-" * 60)
    ao_labels = mol.ao_labels()
    shell_offset = 0
    for ish in range(mol.nbas):
        atom_idx = mol.bas_atom(ish)
        atom_symbol = mol.atom_symbol(atom_idx)
        ang_mom = mol.bas_angular(ish)
        n_cart = (ang_mom + 1) * (ang_mom + 2) // 2
        ang_label = ["s", "p", "d", "f", "g"][ang_mom]
        print(
            f"  Shell {ish:2d}: Atom {atom_idx} ({atom_symbol:2s}) "
            f"l={ang_mom} ({ang_label:1s})  "
            f"AO range: [{shell_offset}, {shell_offset + n_cart})"
        )
        shell_offset += n_cart
    print()

    # Compute ERIs
    print("Computing ERIs...")
    eri = mol.intor("int2e")
    print(f"  ERI tensor shape: {eri.shape}")
    print()

    # Select specific ERIs for validation
    n = mol.nao

    # Define test cases: (i, j, k, l, description)
    # Categories:
    # 1. Diagonal ERIs (iiii) - various orbital types
    # 2. Off-diagonal same-atom ERIs
    # 3. Cross-center ERIs (C-C, C-H)
    # 4. ERIs involving d-orbitals

    test_cases = []

    # --- Diagonal ERIs (iiii) ---
    # C1 orbitals (indices 0-14)
    test_cases.append((0, 0, 0, 0, "C1 1s diagonal"))
    test_cases.append((1, 1, 1, 1, "C1 2s diagonal"))
    test_cases.append((2, 2, 2, 2, "C1 3s diagonal"))
    test_cases.append((3, 3, 3, 3, "C1 2px diagonal"))
    test_cases.append((4, 4, 4, 4, "C1 2py diagonal"))
    test_cases.append((5, 5, 5, 5, "C1 2pz diagonal"))
    test_cases.append((6, 6, 6, 6, "C1 3px diagonal"))
    test_cases.append((7, 7, 7, 7, "C1 3py diagonal"))
    test_cases.append((8, 8, 8, 8, "C1 3pz diagonal"))
    test_cases.append((9, 9, 9, 9, "C1 3dxx diagonal"))
    test_cases.append((10, 10, 10, 10, "C1 3dxy diagonal"))
    test_cases.append((11, 11, 11, 11, "C1 3dxz diagonal"))
    test_cases.append((12, 12, 12, 12, "C1 3dyy diagonal"))
    test_cases.append((13, 13, 13, 13, "C1 3dyz diagonal"))
    test_cases.append((14, 14, 14, 14, "C1 3dzz diagonal"))

    # C2 orbitals (indices 15-29)
    test_cases.append((15, 15, 15, 15, "C2 1s diagonal"))
    test_cases.append((24, 24, 24, 24, "C2 3dxx diagonal"))

    # H atom orbitals (indices 30-37)
    test_cases.append((30, 30, 30, 30, "H1 1s diagonal"))
    test_cases.append((31, 31, 31, 31, "H1 2s diagonal"))
    test_cases.append((32, 32, 32, 32, "H2 1s diagonal"))
    test_cases.append((34, 34, 34, 34, "H3 1s diagonal"))
    test_cases.append((36, 36, 36, 36, "H4 1s diagonal"))

    # --- Off-diagonal same-atom ERIs ---
    test_cases.append((0, 0, 1, 1, "C1 (1s1s|2s2s)"))
    test_cases.append((0, 0, 2, 2, "C1 (1s1s|3s3s)"))
    test_cases.append((0, 1, 0, 1, "C1 (1s2s|1s2s)"))
    test_cases.append((1, 1, 3, 3, "C1 (2s2s|2px2px)"))
    test_cases.append((3, 3, 4, 4, "C1 (2px2px|2py2py)"))
    test_cases.append((0, 0, 9, 9, "C1 (1s1s|3dxx3dxx)"))
    test_cases.append((9, 9, 12, 12, "C1 (3dxx3dxx|3dyy3dyy)"))
    test_cases.append((10, 10, 11, 11, "C1 (3dxy3dxy|3dxz3dxz)"))

    # --- Cross-center ERIs (C-C) ---
    test_cases.append((0, 0, 15, 15, "C1-C2 (1s1s|1s1s)"))
    test_cases.append((0, 15, 0, 15, "C1-C2 (1s1s|1s1s) exchange"))
    test_cases.append((3, 3, 18, 18, "C1-C2 (2px2px|2px2px)"))
    test_cases.append((5, 5, 20, 20, "C1-C2 (2pz2pz|2pz2pz)"))
    test_cases.append((9, 9, 24, 24, "C1-C2 (3dxx3dxx|3dxx3dxx)"))

    # --- Cross-center ERIs (C-H) ---
    test_cases.append((0, 0, 30, 30, "C1-H1 (1s1s|1s1s)"))
    test_cases.append((0, 30, 0, 30, "C1-H1 (1s1s|1s1s) exchange"))
    test_cases.append((9, 9, 30, 30, "C1-H1 (3dxx3dxx|1s1s)"))
    test_cases.append((14, 14, 31, 31, "C1-H1 (3dzz3dzz|2s2s)"))

    # --- Cross-center ERIs (H-H) ---
    test_cases.append((30, 30, 32, 32, "H1-H2 (1s1s|1s1s)"))
    test_cases.append((30, 32, 30, 32, "H1-H2 (1s1s|1s1s) exchange"))
    test_cases.append((30, 30, 34, 34, "H1-H3 (1s1s|1s1s)"))

    # --- Mixed ERIs ---
    test_cases.append((0, 1, 2, 3, "C1 (1s2s|3s2px)"))
    test_cases.append((5, 6, 7, 8, "C1 (2pz3px|3py3pz)"))
    test_cases.append((0, 15, 30, 32, "C1-C2-H1-H2 mixed"))

    print("Reference ERI Values:")
    print("-" * 80)
    print(f"{'Test Case':<45} {'(i,j,k,l)':<15} {'Value':>20}")
    print("-" * 80)

    # Store values for Rust test generation
    rust_test_data = []

    for i, j, k, l, desc in test_cases:
        if i < n and j < n and k < n and l < n:
            val = eri[i, j, k, l]
            print(f"{desc:<45} ({i:2d},{j:2d},{k:2d},{l:2d})   {val:>20.15e}")
            rust_test_data.append((i, j, k, l, val, desc))
        else:
            print(f"{desc:<45} ({i:2d},{j:2d},{k:2d},{l:2d})   SKIPPED (out of range)")

    print("-" * 80)
    print()

    # Generate Rust test code
    print("=" * 80)
    print("Rust Test Code for IQCP:")
    print("=" * 80)
    print()
    print(
        """#[test]
fn test_golden_c2h4_631gs_eri() {
    // Reference values from PySCF 2.11.0 with CARTESIAN d-orbitals
    //
    // C2H4 (ethylene) with 6-31G* basis
    // Geometry (Bohr):
    //   C   0.0000   0.0000   1.2654
    //   C   0.0000   0.0000  -1.2654
    //   H   0.0000   1.7453   2.3280
    //   H   0.0000  -1.7453   2.3280
    //   H   0.0000   1.7453  -2.3280
    //   H   0.0000  -1.7453  -2.3280
    //
    // Generated with:
    // ```python
    // from pyscf import gto
    // mol = gto.Mole()
    // mol.atom = '''
    //     C   0.0000000   0.0000000   1.2654000
    //     C   0.0000000   0.0000000  -1.2654000
    //     H   0.0000000   1.7453000   2.3280000
    //     H   0.0000000  -1.7453000   2.3280000
    //     H   0.0000000   1.7453000  -2.3280000
    //     H   0.0000000  -1.7453000  -2.3280000
    // '''
    // mol.basis = '6-31g*'
    // mol.unit = 'B'
    // mol.cart = True
    // mol.build()
    // eri = mol.intor('int2e')
    // ```
    //
    // Number of basis functions: 38
    // - C1: indices 0-14 (15 functions)
    // - C2: indices 15-29 (15 functions)
    // - H1: indices 30-31 (2 functions)
    // - H2: indices 32-33 (2 functions)
    // - H3: indices 34-35 (2 functions)
    // - H4: indices 36-37 (2 functions)

    let c1 = Atom::new(6, [0.0, 0.0, 1.2654]).unwrap();
    let c2 = Atom::new(6, [0.0, 0.0, -1.2654]).unwrap();
    let h1 = Atom::new(1, [0.0, 1.7453, 2.3280]).unwrap();
    let h2 = Atom::new(1, [0.0, -1.7453, 2.3280]).unwrap();
    let h3 = Atom::new(1, [0.0, 1.7453, -2.3280]).unwrap();
    let h4 = Atom::new(1, [0.0, -1.7453, -2.3280]).unwrap();
    let basis = BasisSet::build(vec![c1, c2, h1, h2, h3, h4], "6-31g*").unwrap();

    assert_eq!(
        basis.n_basis, 38,
        "C2H4 6-31G* should have 38 basis functions (Cartesian d-orbitals)"
    );

    let eri = eri_compressed(&basis);
    let n = basis.n_basis;

    // Test tolerance
    const TOL: f64 = 1e-10;
"""
    )

    # Group tests by category
    print("    // === Diagonal ERIs (iiii) ===")
    for i, j, k, l, val, desc in rust_test_data:
        if i == j == k == l:
            print(
                f'    assert_abs_diff_eq!(eri_get(&eri, n, {i}, {j}, {k}, {l}), {val:.15e}, epsilon = TOL);  // {desc}'
            )

    print()
    print("    // === Off-diagonal same-atom ERIs ===")
    for i, j, k, l, val, desc in rust_test_data:
        if i != j or i != k or i != l:
            if "C1-" not in desc and "H1-" not in desc and "mixed" not in desc:
                if desc.startswith("C1 "):
                    print(
                        f'    assert_abs_diff_eq!(eri_get(&eri, n, {i}, {j}, {k}, {l}), {val:.15e}, epsilon = TOL);  // {desc}'
                    )

    print()
    print("    // === Cross-center ERIs (C-C) ===")
    for i, j, k, l, val, desc in rust_test_data:
        if "C1-C2" in desc:
            print(
                f'    assert_abs_diff_eq!(eri_get(&eri, n, {i}, {j}, {k}, {l}), {val:.15e}, epsilon = TOL);  // {desc}'
            )

    print()
    print("    // === Cross-center ERIs (C-H) ===")
    for i, j, k, l, val, desc in rust_test_data:
        if "C1-H1" in desc:
            print(
                f'    assert_abs_diff_eq!(eri_get(&eri, n, {i}, {j}, {k}, {l}), {val:.15e}, epsilon = TOL);  // {desc}'
            )

    print()
    print("    // === Cross-center ERIs (H-H) ===")
    for i, j, k, l, val, desc in rust_test_data:
        if "H1-H" in desc:
            print(
                f'    assert_abs_diff_eq!(eri_get(&eri, n, {i}, {j}, {k}, {l}), {val:.15e}, epsilon = TOL);  // {desc}'
            )

    print()
    print("    // === Mixed ERIs ===")
    for i, j, k, l, val, desc in rust_test_data:
        if "mixed" in desc or (
            desc.startswith("C1 ") and "diagonal" not in desc and "(" in desc
        ):
            if any(
                x in desc for x in ["(1s2s", "(2pz3px", "C1-C2-H1-H2"]
            ):  # Filter specific mixed cases
                print(
                    f'    assert_abs_diff_eq!(eri_get(&eri, n, {i}, {j}, {k}, {l}), {val:.15e}, epsilon = TOL);  // {desc}'
                )

    print("}")
    print()

    # Also print summary statistics
    print("=" * 80)
    print("Summary Statistics:")
    print("=" * 80)

    all_diag = [eri[i, i, i, i] for i in range(n)]
    print(f"  Diagonal ERIs (iiii):")
    print(f"    Min: {min(all_diag):.10e}")
    print(f"    Max: {max(all_diag):.10e}")
    print(f"    Mean: {np.mean(all_diag):.10e}")

    # Check symmetry
    print()
    print("  Symmetry verification:")
    max_asym = 0.0
    for i in range(min(10, n)):
        for j in range(min(10, n)):
            for k in range(min(10, n)):
                for l in range(min(10, n)):
                    val1 = eri[i, j, k, l]
                    val2 = eri[j, i, k, l]
                    val3 = eri[i, j, l, k]
                    val4 = eri[k, l, i, j]
                    max_asym = max(
                        max_asym,
                        abs(val1 - val2),
                        abs(val1 - val3),
                        abs(val1 - val4),
                    )
    print(f"    Max asymmetry (first 10x10x10x10): {max_asym:.2e}")


if __name__ == "__main__":
    main()
