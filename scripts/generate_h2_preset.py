#!/usr/bin/env python3
"""
Generate H2 STO-3G preset integrals for IQCP.

This script generates the preset JSON file for the H2 molecule
at R=1.4 bohr with STO-3G basis set.

Usage:
    source .venv/bin/activate
    python scripts/generate_h2_preset.py

Output:
    content/presets/systems/h2_sto3g_r1.4.json

Reference:
    - PySCF 2.11.0
    - US-015a: Preset Format & H2
"""
import json
import sys
from pathlib import Path

import numpy as np

try:
    import pyscf
    from pyscf import gto, scf
except ImportError:
    print("Error: PySCF not found. Please activate the virtual environment:")
    print("  source .venv/bin/activate")
    sys.exit(1)


def compress_eri(eri_4d: np.ndarray, nbf: int) -> list[float]:
    """
    Compress 4D ERI tensor to 1D with 8-fold symmetry.

    Uses the indexing scheme from TDD Section 8.3.3:
    - pair(i,j) = i*(i+1)/2 + j (for i >= j)
    - idx(P,Q) = P*(P+1)/2 + Q (for P >= Q)

    Args:
        eri_4d: Full 4D ERI tensor (nbf, nbf, nbf, nbf)
        nbf: Number of basis functions

    Returns:
        Compressed 1D array of unique ERIs
    """
    def pair_index(i: int, j: int) -> int:
        if i < j:
            i, j = j, i
        return i * (i + 1) // 2 + j

    n_pairs = nbf * (nbf + 1) // 2
    n_unique = n_pairs * (n_pairs + 1) // 2
    eri_compressed = np.zeros(n_unique)

    for i in range(nbf):
        for j in range(nbf):
            p = pair_index(i, j)
            for k in range(nbf):
                for l in range(nbf):
                    q = pair_index(k, l)
                    if p >= q:
                        idx = p * (p + 1) // 2 + q
                        eri_compressed[idx] = eri_4d[i, j, k, l]

    return eri_compressed.tolist()


def main():
    print(f"PySCF version: {pyscf.__version__}")
    print()

    # Build molecule
    mol = gto.Mole()
    mol.atom = "H 0 0 0; H 0 0 1.4"  # 1.4 bohr bond length
    mol.basis = "sto-3g"
    mol.unit = "bohr"
    mol.build()

    print(f"System: H2 (STO-3G, R=1.4 bohr)")
    print(f"Number of basis functions: {mol.nao}")
    print(f"Number of electrons: {mol.nelectron}")

    # Get one-electron integrals
    S = mol.intor("int1e_ovlp")  # Overlap
    T = mol.intor("int1e_kin")   # Kinetic
    V = mol.intor("int1e_nuc")   # Nuclear attraction
    H = T + V                     # Core Hamiltonian

    # Get two-electron integrals
    eri_4d = mol.intor("int2e")  # (nbf, nbf, nbf, nbf)

    # Compress ERIs
    eri_compressed = compress_eri(eri_4d, mol.nao)

    # Nuclear repulsion energy
    e_nuc = mol.energy_nuc()

    # Run SCF for reference energy
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    e_tot = mf.kernel()

    print()
    print("=" * 50)
    print("Integral Summary")
    print("=" * 50)
    print(f"Nuclear repulsion: {e_nuc:.15f} Ha")
    print()
    print("Overlap matrix S:")
    print(f"  S[0,0] = {S[0,0]:.15f}")
    print(f"  S[0,1] = {S[0,1]:.15f}")
    print(f"  S[1,0] = {S[1,0]:.15f}")
    print(f"  S[1,1] = {S[1,1]:.15f}")
    print()
    print("Core Hamiltonian H = T + V:")
    print(f"  H[0,0] = {H[0,0]:.15f}")
    print(f"  H[0,1] = {H[0,1]:.15f}")
    print(f"  H[1,0] = {H[1,0]:.15f}")
    print(f"  H[1,1] = {H[1,1]:.15f}")
    print()
    print("Two-electron integrals (compressed):")
    print(f"  (00|00) = {eri_compressed[0]:.15f}")
    print(f"  (01|00) = {eri_compressed[1]:.15f}")
    print(f"  (01|01) = {eri_compressed[2]:.15f}")
    print(f"  (11|00) = {eri_compressed[3]:.15f}")
    print(f"  (11|01) = {eri_compressed[4]:.15f}")
    print(f"  (11|11) = {eri_compressed[5]:.15f}")
    print()
    print("=" * 50)
    print("SCF Results")
    print("=" * 50)
    print(f"Converged: {mf.converged}")
    print(f"Iterations: {mf.cycles}")
    print(f"Total energy: {e_tot:.15f} Ha")
    print(f"Electronic energy: {mf.e_tot - e_nuc:.15f} Ha")
    print()
    print("MO energies:")
    print(f"  HOMO (occupied): {mf.mo_energy[0]:.15f} Ha")
    print(f"  LUMO (virtual):  {mf.mo_energy[1]:.15f} Ha")
    print()

    # Build JSON preset
    preset = {
        "format_version": 1,
        "system_id": "h2_sto3g_r1.4",
        "label": "H2 (STO-3G, R=1.4 bohr)",
        "description": "Hydrogen molecule at equilibrium bond length. Reference system for IQCP Module E.",
        "geometry": {
            "atoms": [
                {"symbol": "H", "xyz": [0.0, 0.0, 0.0]},
                {"symbol": "H", "xyz": [0.0, 0.0, 1.4]}
            ],
            "units": "bohr"
        },
        "basis_id": "sto-3g",
        "nbf": int(mol.nao),
        "nelec": int(mol.nelectron),
        "e_nuc": float(e_nuc),
        "s_matrix": S.flatten().tolist(),
        "h_core": H.flatten().tolist(),
        "eri_compressed": eri_compressed,
        "eri_indexing": "8-fold symmetry: pair(i,j)=i*(i+1)/2+j, idx(P,Q)=P*(P+1)/2+Q",
        "reference": {
            "software": "PySCF",
            "version": pyscf.__version__,
            "energy": float(e_tot)
        }
    }

    # Save to file
    output_path = Path(__file__).parent.parent / "content" / "presets" / "systems" / "h2_sto3g_r1.4.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(preset, f, indent=2)

    print(f"Preset saved to: {output_path}")
    print()
    print("Verification command:")
    print(f"  cargo test --package qc-core preset_json")


if __name__ == "__main__":
    main()
