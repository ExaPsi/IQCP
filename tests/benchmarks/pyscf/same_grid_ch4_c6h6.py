#!/usr/bin/env python3
"""
Same-Grid DFT Validation for CH4 and C6H6
==========================================

Extension of same_grid_validation.py to CH4 and C6H6 molecules.
Demonstrates whether the systematic ~50 uHa CH4/LDA offset disappears
under matched grids, confirming grid vs algorithmic origin.

Usage:
    source .venv/bin/activate
    python tests/benchmarks/pyscf/same_grid_ch4_c6h6.py

Output:
    - Summary table to stdout
    - JSON results to docs/Manuscripts/M2/Data/same_grid_ch4_c6h6.json
"""

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from pyscf import __version__ as pyscf_version
from pyscf import dft, gto, lib

# Force single-threaded execution for deterministic results
lib.num_threads(1)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent.parent
OUTPUT_DIR = PROJECT_ROOT / "docs" / "Manuscripts" / "M2" / "Data"
IQCP_DATA = PROJECT_ROOT / "tests" / "benchmarks" / "criterion" / "iqcp_iteration_counts.json"

# -------------------------------------------------------------------------
# IQCP molecule geometries (coordinates in Bohr)
# From crates/qc-core/src/dft/ks_scf.rs: ch4_atoms(), c6h6_atoms()
# -------------------------------------------------------------------------

# CH4: tetrahedral geometry, C at origin
CH4_GEOM = """
C  0.0       0.0       0.0
H  1.1851    1.1851    1.1851
H -1.1851   -1.1851    1.1851
H -1.1851    1.1851   -1.1851
H  1.1851   -1.1851   -1.1851
"""

# C6H6: planar D6h benzene
# From crates/qc-core/src/scf/pes_internal.rs (benchmark_spherical.rs)
C6H6_GEOM = """
C  0.0000000000  2.6399473960  0.0000000000
C  2.2861906655  1.3199736980  0.0000000000
C  2.2861906655 -1.3199736980  0.0000000000
C  0.0000000000 -2.6399473960  0.0000000000
C -2.2861906655 -1.3199736980  0.0000000000
C -2.2861906655  1.3199736980  0.0000000000
H  0.0000000000  4.6884105150  0.0000000000
H  4.0602655512  2.3442052575  0.0000000000
H  4.0602655512 -2.3442052575  0.0000000000
H  0.0000000000 -4.6884105150  0.0000000000
H -4.0602655512 -2.3442052575  0.0000000000
H -4.0602655512  2.3442052575  0.0000000000
"""

BASIS_SETS = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"]
DFT_METHODS = [
    ("LDA", "lda,vwn5"),
    ("B3LYP", "b3lyp5"),
]
CONV_TOL = 1e-12


def load_iqcp_energies():
    """Load IQCP computed energies from iteration counts JSON."""
    with open(IQCP_DATA) as f:
        data = json.load(f)
    iqcp = {}
    for entry in data["results"]:
        key = (entry["molecule"], entry["basis"], entry["method"])
        iqcp[key] = entry["energy"]
    return iqcp


def build_mol(geom_str, basis_name):
    """Build PySCF Mole object."""
    mol = gto.Mole()
    mol.atom = geom_str
    mol.unit = "bohr"
    mol.basis = basis_name
    mol.cart = True  # Cartesian d-functions (6d, matching IQCP)
    mol.verbose = 0
    mol.build()
    return mol


def run_dft_same_grid(mol, xc_string):
    """Run DFT with IQCP-equivalent SG-1 grid."""
    mf = dft.RKS(mol)
    mf.xc = xc_string
    mf.conv_tol = CONV_TOL
    mf.max_cycle = 200

    # Match IQCP's SG-1 grid exactly:
    # 1. Mura-Knowles Log3 radial quadrature (75 points)
    # 2. SG-1 pruning with zones [6, 38, 86, 194, 86]
    # 3. Original Becke partition (NO atomic-size adjustment)
    mf.grids.radi_method = dft.radi.mura_knowles
    mf.grids.prune = dft.gen_grid.sg1_prune
    mf.grids.becke_scheme = dft.gen_grid.original_becke
    mf.grids.radii_adjust = None
    mf.grids.atom_grid = (75, 194)
    mf.grids.build()

    energy = mf.kernel()
    n_grid = mf.grids.weights.shape[0]
    converged = mf.converged
    return energy, n_grid, converged


def run_dft_default_grid(mol, xc_string):
    """Run DFT with PySCF default grid."""
    mf = dft.RKS(mol)
    mf.xc = xc_string
    mf.conv_tol = CONV_TOL
    mf.max_cycle = 200
    energy = mf.kernel()
    n_grid = mf.grids.weights.shape[0]
    converged = mf.converged
    return energy, n_grid, converged


def main():
    print("=" * 110)
    print("Same-Grid DFT Validation: CH4 and C6H6")
    print("=" * 110)
    print()
    print(f"PySCF version: {pyscf_version}")
    print(f"NumPy version: {np.__version__}")
    print(f"Date: {datetime.now(timezone.utc).isoformat()}")
    print()

    iqcp_energies = load_iqcp_energies()
    results = []

    molecules = [
        ("CH4", CH4_GEOM),
        ("C6H6", C6H6_GEOM),
    ]

    print("-" * 110)
    print(f"{'Mol':<6} {'Basis':<10} {'Method':<8} {'IQCP Energy':>18} {'PySCF(same)':>18} "
          f"{'|dE| same':>12} {'PySCF(default)':>18} {'|dE| default':>12}")
    print(f"{'':6} {'':10} {'':8} {'(Ha)':>18} {'(Ha)':>18} "
          f"{'(uHa)':>12} {'(Ha)':>18} {'(uHa)':>12}")
    print("-" * 110)

    for mol_name, geom in molecules:
        for basis_name in BASIS_SETS:
            mol = build_mol(geom, basis_name)

            for method_label, xc_string in DFT_METHODS:
                iqcp_key = (mol_name, basis_name, method_label)
                iqcp_energy = iqcp_energies.get(iqcp_key)

                if iqcp_energy is None:
                    print(f"{mol_name:<6} {basis_name:<10} {method_label:<8} {'N/A':>18}")
                    continue

                e_same, n_grid_same, conv_same = run_dft_same_grid(mol, xc_string)
                e_default, n_grid_default, conv_default = run_dft_default_grid(mol, xc_string)

                delta_same = abs(e_same - iqcp_energy)
                delta_default = abs(e_default - iqcp_energy)

                print(f"{mol_name:<6} {basis_name:<10} {method_label:<8} {iqcp_energy:>18.12f} "
                      f"{e_same:>18.12f} {delta_same*1e6:>12.4f} "
                      f"{e_default:>18.12f} {delta_default*1e6:>12.4f}")

                results.append({
                    "molecule": mol_name,
                    "basis": basis_name,
                    "method": method_label,
                    "xc_string": xc_string,
                    "iqcp_energy": iqcp_energy,
                    "pyscf_same_grid_energy": e_same,
                    "pyscf_default_grid_energy": e_default,
                    "delta_same_grid_Ha": delta_same,
                    "delta_same_grid_uHa": delta_same * 1e6,
                    "delta_same_grid_nHa": delta_same * 1e9,
                    "delta_default_grid_Ha": delta_default,
                    "delta_default_grid_uHa": delta_default * 1e6,
                    "n_grid_same": n_grid_same,
                    "n_grid_default": n_grid_default,
                    "converged_same": conv_same,
                    "converged_default": conv_default,
                })

    print("-" * 110)

    # Split summary by molecule
    for mol_name in ["CH4", "C6H6"]:
        mol_results = [r for r in results if r["molecule"] == mol_name]
        if not mol_results:
            continue
        same_deltas = [r["delta_same_grid_uHa"] for r in mol_results]
        default_deltas = [r["delta_default_grid_uHa"] for r in mol_results]
        print()
        print(f"  {mol_name} SUMMARY:")
        print(f"    Same-grid max:    {max(same_deltas):.4f} uHa = {max(same_deltas)*1e3:.2f} nHa")
        print(f"    Same-grid mean:   {np.mean(same_deltas):.4f} uHa = {np.mean(same_deltas)*1e3:.2f} nHa")
        print(f"    Default-grid max: {max(default_deltas):.4f} uHa")
        print(f"    Default-grid mean:{np.mean(default_deltas):.4f} uHa")

    # Save JSON
    output = {
        "metadata": {
            "description": "Same-grid DFT validation for CH4 and C6H6",
            "generated": datetime.now(timezone.utc).isoformat(),
            "pyscf_version": pyscf_version,
            "numpy_version": np.__version__,
            "convergence_threshold": CONV_TOL,
            "molecules": {
                "CH4": {"geometry_bohr": CH4_GEOM.strip(), "cartesian_d": True},
                "C6H6": {"geometry_bohr": C6H6_GEOM.strip(), "cartesian_d": True},
            },
        },
        "grid_configurations": {
            "same_grid": {
                "description": "IQCP-equivalent SG-1 grid",
                "radial_method": "Mura-Knowles Log3",
                "n_radial": 75,
                "max_angular": 194,
                "pruning": "SG-1 [6, 38, 86, 194, 86]",
                "partition": "Original Becke (1988)",
                "atomic_size_adjustment": False,
            },
            "default_grid": {
                "description": "PySCF default grid (level=3)",
                "radial_method": "Treutler-Ahlrichs",
                "pruning": "nwchem_prune",
                "partition": "Original Becke with Treutler atomic-size adjustment",
                "atomic_size_adjustment": True,
            },
        },
        "results": results,
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "same_grid_ch4_c6h6.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
