#!/usr/bin/env python3
"""
Generate PySCF analytical gradients for 3-21G and 6-31+G* basis sets.
=====================================================================

Extends gradient validation (SI Table S3) to cover basis sets with the
largest energy deviations, as requested by Reviewer 3.

Systems: H2, HF, H2O, NH3, CH4
Methods: RHF, LDA (lda,vwn5), B3LYP (b3lyp5)
Basis sets: 3-21g, 6-31+g*

All geometries use IQCP coordinates (Bohr), Cartesian d-functions (cart=True),
and conv_tol=1e-12.

Output: tests/golden/dft/gradient_validation_extended.json

Run:
    source .venv/bin/activate
    python tests/benchmarks/pyscf/generate_extended_gradients.py
"""

import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from pyscf import dft, gto, lib, scf

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CONV_TOL = 1e-12
lib.num_threads(1)

OUTPUT_PATH = Path(__file__).resolve().parents[2] / "golden" / "dft" / "gradient_validation_extended.json"

# ---------------------------------------------------------------------------
# Molecule definitions (coordinates in Bohr, matching IQCP exactly)
# ---------------------------------------------------------------------------

MOLECULES = {
    "H2": {
        "atoms": "H 0.0 0.0 0.0; H 0.0 0.0 1.4",
        "charge": 0,
        "spin": 0,
    },
    "HF": {
        "atoms": "H 0.0 0.0 0.0; F 0.0 0.0 1.7328",
        "charge": 0,
        "spin": 0,
    },
    "H2O": {
        "atoms": "O 0.0 0.0 0.2217282; H 0.0 1.4305447 -0.8869128; H 0.0 -1.4305447 -0.8869128",
        "charge": 0,
        "spin": 0,
    },
    "NH3": {
        "atoms": "N 0.0 0.0 0.219705; H 0.0 1.7714918 -0.512645; H 1.5342036 -0.8857459 -0.512645; H -1.5342036 -0.8857459 -0.512645",
        "charge": 0,
        "spin": 0,
    },
    "CH4": {
        "atoms": "C 0 0 0; H 1.1851 1.1851 1.1851; H -1.1851 -1.1851 1.1851; H -1.1851 1.1851 -1.1851; H 1.1851 -1.1851 -1.1851",
        "charge": 0,
        "spin": 0,
    },
}

BASIS_SETS = ["3-21g", "6-31+g*"]
METHODS = ["RHF", "LDA", "B3LYP"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def build_mol(mol_name: str, basis: str) -> gto.Mole:
    """Build a PySCF Mole with Cartesian d-functions."""
    mol_def = MOLECULES[mol_name]
    mol = gto.Mole()
    mol.atom = mol_def["atoms"]
    mol.basis = basis
    mol.charge = mol_def["charge"]
    mol.spin = mol_def["spin"]
    mol.unit = "Bohr"
    mol.cart = True
    mol.verbose = 0
    mol.build()
    return mol


def run_scf_and_gradient(mol: gto.Mole, method: str):
    """Run SCF + analytical gradient, return (energy, gradient, converged, n_iter)."""
    if method == "RHF":
        mf = scf.RHF(mol)
    elif method == "LDA":
        mf = dft.RKS(mol)
        mf.xc = "lda,vwn5"
    elif method == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "b3lyp5"
    else:
        raise ValueError(f"Unknown method: {method}")

    mf.conv_tol = CONV_TOL
    energy = mf.kernel()

    grad_obj = mf.nuc_grad_method()
    grad_obj.verbose = 0
    g = grad_obj.kernel()

    return energy, g, mf.converged, getattr(mf, "cycles", None)


def gradient_stats(grad: np.ndarray):
    """Compute max and RMS gradient magnitudes."""
    flat = np.abs(grad).flatten()
    return float(np.max(flat)), float(np.sqrt(np.mean(flat**2)))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    import pyscf
    import scipy

    print("PySCF Extended Gradient Generation (3-21G, 6-31+G*)")
    print("=" * 65)
    print(f"PySCF {pyscf.__version__}  |  NumPy {np.__version__}  |  SciPy {scipy.__version__}")
    print(f"Python {platform.python_version()}")
    print(f"Cartesian d-functions: True")
    print(f"conv_tol: {CONV_TOL}")
    print()

    metadata = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "description": "PySCF analytical gradients for 3-21G and 6-31+G* (cart=True, b3lyp5)",
        "pyscf_version": pyscf.__version__,
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "python_version": platform.python_version(),
        "conv_tol": CONV_TOL,
        "cartesian_d": True,
        "basis_sets": BASIS_SETS,
        "methods": METHODS,
    }

    results = []
    total = len(MOLECULES) * len(BASIS_SETS) * len(METHODS)
    count = 0

    for mol_name in MOLECULES:
        for basis in BASIS_SETS:
            mol = build_mol(mol_name, basis)
            n_basis = mol.nao_nr()

            for method in METHODS:
                count += 1
                label = f"{mol_name}/{basis}/{method}"
                print(f"  [{count:2d}/{total}] {label:<35s}", end="", flush=True)

                try:
                    energy, grad, converged, n_iter = run_scf_and_gradient(mol, method)
                    max_grad, rms_grad = gradient_stats(grad)

                    result = {
                        "molecule": mol_name,
                        "basis": basis,
                        "method": method,
                        "n_basis": n_basis,
                        "energy": energy,
                        "gradient": grad.tolist(),
                        "max_gradient": max_grad,
                        "rms_gradient": rms_grad,
                        "converged": bool(converged),
                    }
                    results.append(result)

                    status = "OK" if converged else "NOT CONVERGED"
                    print(f"  E={energy:18.12f}  |grad|_max={max_grad:.6e}  [{status}]")

                except Exception as e:
                    print(f"  FAILED: {e}")
                    results.append({
                        "molecule": mol_name,
                        "basis": basis,
                        "method": method,
                        "n_basis": n_basis,
                        "error": str(e),
                    })

    # Summary table
    print()
    print("=" * 100)
    print(f"{'Molecule':<8s} {'Basis':<10s} {'Method':<8s} {'N_bas':>5s} "
          f"{'Energy (Ha)':>18s} {'Max|grad|':>12s} {'RMS|grad|':>12s}")
    print("-" * 100)

    for r in results:
        if "error" in r:
            print(f"{r['molecule']:<8s} {r['basis']:<10s} {r['method']:<8s} "
                  f"{r['n_basis']:5d}  ERROR: {r['error']}")
        else:
            print(f"{r['molecule']:<8s} {r['basis']:<10s} {r['method']:<8s} "
                  f"{r['n_basis']:5d} {r['energy']:18.12f} "
                  f"{r['max_gradient']:12.6e} {r['rms_gradient']:12.6e}")
    print("=" * 100)

    # Save output
    output_data = {**metadata, "results": results}
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(output_data, f, indent=2)

    print(f"\nResults saved to: {OUTPUT_PATH}")
    print(f"Total entries: {len(results)} ({len([r for r in results if 'error' not in r])} successful)")


if __name__ == "__main__":
    main()
