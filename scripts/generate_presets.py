#!/usr/bin/env python3
"""
Generate preset integrals for IQCP Module E.

This script generates JSON preset files for molecular systems
using PySCF to compute integrals and reference SCF energies.

Usage:
    source .venv/bin/activate
    python scripts/generate_presets.py                    # Generate all systems
    python scripts/generate_presets.py heh_plus          # Generate specific system
    python scripts/generate_presets.py --list            # List available systems

Output:
    content/presets/systems/<system_id>.json

Reference:
    - PySCF 2.11.0
    - US-015b: Additional Presets
"""
import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

try:
    import pyscf
    from pyscf import gto, scf
except ImportError:
    print("Error: PySCF not found. Please activate the virtual environment:")
    print("  source .venv/bin/activate")
    sys.exit(1)


# System definitions with geometries in bohr
SYSTEMS: dict[str, dict[str, Any]] = {
    "heh_plus_sto3g": {
        "atom": "He 0 0 0; H 0 0 1.4632",
        "charge": 1,  # HeH+ is a cation
        "label": "HeH+ (STO-3G)",
        "description": "Helium hydride ion - heteronuclear diatomic with 2 electrons. "
                       "Good comparison to H2 as same basis size but different nuclear charges.",
        "output_file": "heh_plus_sto3g.json",
    },
    "lih_sto3g": {
        "atom": "Li 0 0 0; H 0 0 3.0139",
        "charge": 0,
        "label": "LiH (STO-3G)",
        "description": "Lithium hydride - polar diatomic with 4 electrons. "
                       "First system with core and valence orbitals.",
        "output_file": "lih_sto3g.json",
    },
    "h2o_sto3g": {
        "atom": "O 0 0 0.2217282; H 0 1.4305447 -0.8869128; H 0 -1.4305447 -0.8869128",
        "charge": 0,
        "label": "H2O (STO-3G)",
        "description": "Water molecule - bent triatomic with 10 electrons. "
                       "Familiar molecule demonstrating non-linear geometry.",
        "output_file": "h2o_sto3g.json",
    },
    "nh3_sto3g": {
        "atom": "N 0 0 0.2197050; H 0 1.7714918 -0.5126450; "
                "H 1.5342036 -0.8857459 -0.5126450; H -1.5342036 -0.8857459 -0.5126450",
        "charge": 0,
        "label": "NH3 (STO-3G)",
        "description": "Ammonia - pyramidal molecule with 10 electrons. "
                       "Demonstrates non-planar molecular geometry.",
        "output_file": "nh3_sto3g.json",
    },
}


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


def parse_geometry(atom_str: str) -> list[dict[str, Any]]:
    """Parse PySCF atom string into geometry structure."""
    atoms = []
    for atom_spec in atom_str.split(";"):
        parts = atom_spec.strip().split()
        if len(parts) >= 4:
            symbol = parts[0]
            xyz = [float(parts[1]), float(parts[2]), float(parts[3])]
            atoms.append({"symbol": symbol, "xyz": xyz})
    return atoms


def generate_preset(system_id: str, system_def: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """
    Generate a preset JSON file for a molecular system.

    Args:
        system_id: Unique identifier for the system
        system_def: System definition dictionary
        output_dir: Directory to write the JSON file

    Returns:
        Dictionary with generation results
    """
    print(f"\n{'='*60}")
    print(f"Generating: {system_def['label']}")
    print(f"{'='*60}")

    # Build molecule
    mol = gto.Mole()
    mol.atom = system_def["atom"]
    mol.basis = "sto-3g"
    mol.unit = "bohr"
    mol.charge = system_def.get("charge", 0)
    mol.build()

    print(f"System ID: {system_id}")
    print(f"Number of basis functions: {mol.nao}")
    print(f"Number of electrons: {mol.nelectron}")
    print(f"Charge: {mol.charge}")

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
    print("Integrals Summary:")
    print(f"  Nuclear repulsion: {e_nuc:.15f} Ha")
    print(f"  Matrix dimensions: {mol.nao} x {mol.nao}")
    print(f"  Unique ERIs: {len(eri_compressed)}")
    print()
    print("SCF Results:")
    print(f"  Converged: {mf.converged}")
    print(f"  Iterations: {mf.cycles}")
    print(f"  Total energy: {e_tot:.15f} Ha")
    print(f"  Electronic energy: {mf.e_tot - e_nuc:.15f} Ha")
    print()
    print("MO energies:")
    n_occ = mol.nelectron // 2
    for i, e in enumerate(mf.mo_energy):
        occ_label = "(occ)" if i < n_occ else "(vir)"
        homo_lumo = ""
        if i == n_occ - 1:
            homo_lumo = " <- HOMO"
        elif i == n_occ:
            homo_lumo = " <- LUMO"
        print(f"  MO {i:2d} {occ_label}: {e:15.12f} Ha{homo_lumo}")

    # Parse geometry for JSON
    geometry_atoms = parse_geometry(system_def["atom"])

    # Build JSON preset
    preset = {
        "format_version": 1,
        "system_id": system_id,
        "label": system_def["label"],
        "description": system_def["description"],
        "geometry": {
            "atoms": geometry_atoms,
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
    output_file = output_dir / system_def["output_file"]
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        json.dump(preset, f, indent=2)

    print(f"Preset saved to: {output_file}")

    return {
        "system_id": system_id,
        "nbf": mol.nao,
        "nelec": mol.nelectron,
        "n_eri": len(eri_compressed),
        "e_nuc": e_nuc,
        "e_tot": e_tot,
        "converged": mf.converged,
        "iterations": mf.cycles,
        "output_file": str(output_file),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate IQCP preset integrals from PySCF"
    )
    parser.add_argument(
        "systems",
        nargs="*",
        help="System IDs to generate (default: all)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available systems",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent.parent / "content" / "presets" / "systems",
        help="Output directory for JSON files",
    )

    args = parser.parse_args()

    if args.list:
        print("Available systems:")
        for sys_id, sys_def in SYSTEMS.items():
            print(f"  {sys_id}: {sys_def['label']}")
        return

    # Determine which systems to generate
    if args.systems:
        systems_to_generate = {}
        for sys_id in args.systems:
            if sys_id not in SYSTEMS:
                print(f"Error: Unknown system '{sys_id}'")
                print(f"Available: {', '.join(SYSTEMS.keys())}")
                sys.exit(1)
            systems_to_generate[sys_id] = SYSTEMS[sys_id]
    else:
        systems_to_generate = SYSTEMS

    print(f"PySCF version: {pyscf.__version__}")
    print(f"NumPy version: {np.__version__}")
    print(f"Generation time: {datetime.now(timezone.utc).isoformat()}")
    print(f"Output directory: {args.output_dir}")

    # Generate each system
    results = []
    for sys_id, sys_def in systems_to_generate.items():
        result = generate_preset(sys_id, sys_def, args.output_dir)
        results.append(result)

    # Print summary
    print()
    print("=" * 60)
    print("GENERATION SUMMARY")
    print("=" * 60)
    print()
    print(f"{'System':<20} {'nbf':>4} {'nelec':>6} {'n_eri':>7} {'E_total (Ha)':>20} {'Conv':>5}")
    print("-" * 70)
    for r in results:
        conv_str = "Yes" if r["converged"] else "NO!"
        print(f"{r['system_id']:<20} {r['nbf']:>4} {r['nelec']:>6} {r['n_eri']:>7} {r['e_tot']:>20.12f} {conv_str:>5}")
    print()
    print("Verification commands:")
    print("  cargo test --package qc-core preset")


if __name__ == "__main__":
    main()
