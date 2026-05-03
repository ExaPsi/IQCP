#!/usr/bin/env python3
"""Investigate RHF energy discrepancies between IQCP and PySCF.

Systems with >10 nHa errors:
  - HF/6-31G*:   65 nHa
  - NH3/STO-3G:  28 nHa
  - C6H6/6-31G*: 25 nHa

Strategy: Compare individual integral elements (S, H, ERI) between PySCF
and IQCP to identify which integral classes contribute to the error.
"""

import numpy as np
from pyscf import gto, scf, lib
import json

lib.num_threads(1)


def analyze_system(name, mol_spec, basis_name, cart=False, unit="bohr"):
    """Run PySCF and print integral diagnostics for a system."""
    print(f"\n{'='*70}")
    print(f"  {name} / {basis_name}")
    print(f"{'='*70}")

    mol = gto.Mole()
    mol.atom = mol_spec
    mol.unit = unit
    mol.basis = basis_name
    mol.cart = cart
    mol.build()

    n = mol.nao
    print(f"Basis functions: {n}")
    print(f"Cartesian d-functions: {cart}")
    print(f"Shells: {mol.nbas}")

    # Print shell info
    print(f"\nShell structure:")
    for i in range(mol.nbas):
        l = mol.bas_angular(i)
        atom_id = mol.bas_atom(i)
        nprim = mol.bas_nprim(i)
        nctr = mol.bas_nctr(i)
        print(f"  Shell {i}: atom={atom_id}, l={l}, nprim={nprim}, nctr={nctr}")

    # Compute integrals
    s = mol.intor("int1e_ovlp")
    t = mol.intor("int1e_kin")
    v = mol.intor("int1e_nuc")
    h = t + v

    print(f"\nOverlap matrix (S) - selected elements:")
    print(f"  S[0,0] = {s[0,0]:.15e}")
    if n > 1:
        print(f"  S[0,1] = {s[0,1]:.15e}")
        print(f"  S[1,1] = {s[1,1]:.15e}")
    if n > 5:
        print(f"  S[0,5] = {s[0,5]:.15e}")
        print(f"  S[5,5] = {s[5,5]:.15e}")

    print(f"\nCore Hamiltonian (H) - selected elements:")
    print(f"  H[0,0] = {h[0,0]:.15e}")
    if n > 1:
        print(f"  H[0,1] = {h[0,1]:.15e}")
        print(f"  H[1,1] = {h[1,1]:.15e}")
    if n > 5:
        print(f"  H[0,5] = {h[0,5]:.15e}")
        print(f"  H[5,5] = {h[5,5]:.15e}")

    # ERI computation
    eri = mol.intor("int2e")
    print(f"\nERI tensor shape: {eri.shape}")
    print(f"ERI - selected elements:")
    print(f"  (0,0|0,0) = {eri[0,0,0,0]:.15e}")
    if n > 1:
        print(f"  (0,1|0,1) = {eri[0,1,0,1]:.15e}")
        print(f"  (1,1|1,1) = {eri[1,1,1,1]:.15e}")
        print(f"  (0,0|1,1) = {eri[0,0,1,1]:.15e}")
    if n > 5:
        print(f"  (0,0|5,5) = {eri[0,0,5,5]:.15e}")
        print(f"  (5,5|5,5) = {eri[5,5,5,5]:.15e}")

    # Run SCF
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    e = mf.kernel()
    print(f"\nRHF energy: {e:.15f}")
    print(f"Converged: {mf.converged}")
    print(f"Iterations: {mf.cycles}")

    # For d-function analysis: find d-type shells and print their ERIs
    d_shells = []
    d_bf_start = []
    bf_offset = 0
    for i in range(mol.nbas):
        l = mol.bas_angular(i)
        if cart:
            n_bf = (l + 1) * (l + 2) // 2
        else:
            n_bf = 2 * l + 1
        if l == 2:
            d_shells.append(i)
            d_bf_start.append(bf_offset)
        bf_offset += n_bf

    if d_shells:
        print(f"\nd-type shells found: {d_shells}")
        print(f"d-function basis function ranges:")
        for shell_idx, bf_start in zip(d_shells, d_bf_start):
            l = mol.bas_angular(shell_idx)
            n_bf = 6 if cart else 5  # Cartesian d has 6 components
            print(f"  Shell {shell_idx}: BFs [{bf_start}..{bf_start+n_bf})")
            # Print some ERI elements involving d-functions
            for d_bf in range(bf_start, min(bf_start + n_bf, n)):
                print(f"    (0,0|{d_bf},{d_bf}) = {eri[0,0,d_bf,d_bf]:.15e}")
                print(f"    ({d_bf},{d_bf}|{d_bf},{d_bf}) = {eri[d_bf,d_bf,d_bf,d_bf]:.15e}")

        # Print max ERI involving at least one d-function
        max_d_eri = 0.0
        max_d_idx = (0, 0, 0, 0)
        for ds, bf_s in zip(d_shells, d_bf_start):
            n_bf_d = 6 if cart else 5
            for d_bf in range(bf_s, bf_s + n_bf_d):
                for j in range(n):
                    for k in range(n):
                        for l_idx in range(n):
                            val = abs(eri[d_bf, j, k, l_idx])
                            if val > max_d_eri:
                                max_d_eri = val
                                max_d_idx = (d_bf, j, k, l_idx)
        print(f"\n  Max |ERI| involving d-function: {max_d_eri:.15e}")
        print(f"  At indices: {max_d_idx}")
    else:
        print(f"\nNo d-type shells in this basis set.")

    # Save full integral arrays for detailed comparison
    result = {
        "name": name,
        "basis": basis_name,
        "n_basis": n,
        "energy": float(e),
        "converged": bool(mf.converged),
        "S_diag": s.diagonal().tolist(),
        "H_diag": h.diagonal().tolist(),
        "S_flat": s.flatten().tolist(),
        "H_flat": h.flatten().tolist(),
    }

    # For small systems, also save full ERI
    if n <= 21:
        result["ERI_flat"] = eri.flatten().tolist()

    return result, s, h, eri, mf


def main():
    print("=" * 70)
    print("  RHF Energy Discrepancy Investigation")
    print("  PySCF reference integrals for comparison with IQCP")
    print("=" * 70)

    results = {}

    # System 1: HF/6-31G* (65 nHa error - LARGEST)
    r1, s1, h1, eri1, mf1 = analyze_system(
        "HF",
        "H 0 0 0; F 0 0 1.7328",
        "6-31g*",
        cart=True,
    )
    results["HF_631gs"] = r1

    # System 2: NH3/STO-3G (28 nHa error - NO d-functions)
    r2, s2, h2, eri2, mf2 = analyze_system(
        "NH3",
        "N 0.0 0.0 0.219705; H 0.0 1.7714918 -0.512645; "
        "H 1.5342036 -0.8857459 -0.512645; H -1.5342036 -0.8857459 -0.512645",
        "sto-3g",
        cart=False,
    )
    results["NH3_sto3g"] = r2

    # System 3: C6H6/6-31G* (25 nHa error)
    # Use Angstrom coordinates matching the validation.json geometry
    r3, s3, h3, eri3, mf3 = analyze_system(
        "C6H6",
        "C 0 1.397 0; C 1.2098 0.6985 0; C 1.2098 -0.6985 0; "
        "C 0 -1.397 0; C -1.2098 -0.6985 0; C -1.2098 0.6985 0; "
        "H 0 2.481 0; H 2.1486 1.2405 0; H 2.1486 -1.2405 0; "
        "H 0 -2.481 0; H -2.1486 -1.2405 0; H -2.1486 1.2405 0",
        "6-31g*",
        cart=True,
        unit="angstrom",
    )
    results["C6H6_631gs"] = r3

    # Also check a few "good" systems for comparison
    r4, s4, h4, eri4, mf4 = analyze_system(
        "H2O",
        "O 0.0 0.0 0.2217282; H 0.0 1.4305447 -0.8869128; H 0.0 -1.4305447 -0.8869128",
        "6-31g*",
        cart=True,
    )
    results["H2O_631gs"] = r4

    # Summary
    print(f"\n{'='*70}")
    print(f"  SUMMARY: RHF Energy Comparison")
    print(f"{'='*70}")

    iqcp_energies = {
        "HF_631gs": -100.002863192006,
        "NH3_sto3g": -55.454361686086,
        "C6H6_631gs": -230.702101338829,
        "H2O_631gs": -76.010505696475,
    }

    for key, r in results.items():
        pyscf_e = r["energy"]
        iqcp_e = iqcp_energies.get(key, None)
        if iqcp_e is not None:
            diff_nha = (iqcp_e - pyscf_e) * 1e9
            print(f"  {key:15s}: PySCF={pyscf_e:.12f}, IQCP={iqcp_e:.12f}, diff={diff_nha:+.1f} nHa")
        else:
            print(f"  {key:15s}: PySCF={pyscf_e:.12f}")

    # Save reference data
    output_path = "tests/benchmarks/pyscf/results/rhf_discrepancy_integrals.json"
    # Only save small-system data (skip C6H6 full ERI)
    save_results = {k: v for k, v in results.items() if k != "C6H6_631gs"}
    # For C6H6, save only diagonal elements
    save_results["C6H6_631gs"] = {
        "name": results["C6H6_631gs"]["name"],
        "basis": results["C6H6_631gs"]["basis"],
        "n_basis": results["C6H6_631gs"]["n_basis"],
        "energy": results["C6H6_631gs"]["energy"],
        "converged": results["C6H6_631gs"]["converged"],
        "S_diag": results["C6H6_631gs"]["S_diag"],
        "H_diag": results["C6H6_631gs"]["H_diag"],
    }

    with open(output_path, "w") as f:
        json.dump(save_results, f, indent=2)
    print(f"\nReference data saved to {output_path}")


if __name__ == "__main__":
    main()
