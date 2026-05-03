#!/usr/bin/env python3
"""
PySCF Gradient Timing Benchmarks for JCIM M2 Manuscript
=========================================================

Times analytical nuclear gradient calculations for selected molecule/method
combinations. Each includes the full SCF solve + gradient evaluation.

Systems
-------
- H2O RHF/STO-3G
- H2O RHF/6-31G*
- H2O B3LYP/6-31G*
- CH4 B3LYP/6-31G*

Methodology: 10 runs, single-threaded, mean +/- std.
"""

import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from pyscf import dft, gto, lib, scf

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

N_RUNS = 10
CONV_TOL = 1e-12
SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"

lib.num_threads(1)

# ---------------------------------------------------------------------------
# Molecule definitions (bohr)
# ---------------------------------------------------------------------------

MOLECULES = {
    "H2O": {
        "atoms": "O 0.0 0.0 0.2217282; H 0.0 1.4305447 -0.8869128; H 0.0 -1.4305447 -0.8869128",
        "charge": 0,
        "spin": 0,
    },
    "CH4": {
        "atoms": "C 0 0 0; H 1.1851 1.1851 1.1851; H -1.1851 -1.1851 1.1851; H -1.1851 1.1851 -1.1851; H 1.1851 -1.1851 -1.1851",
        "charge": 0,
        "spin": 0,
    },
}

# (molecule, basis, method) combinations to benchmark
GRADIENT_CASES = [
    ("H2O", "sto-3g", "RHF"),
    ("H2O", "6-31g*", "RHF"),
    ("H2O", "6-31g*", "B3LYP"),
    ("CH4", "6-31g*", "B3LYP"),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def build_mol(mol_name: str, basis: str) -> gto.Mole:
    mol_def = MOLECULES[mol_name]
    mol = gto.Mole()
    mol.atom = mol_def["atoms"]
    mol.basis = basis
    mol.charge = mol_def["charge"]
    mol.spin = mol_def["spin"]
    mol.unit = "Bohr"
    mol.verbose = 0
    # NOTE: Do NOT set mol.cart = True for gradient benchmarks.
    # PySCF's analytical gradient implementation does not support
    # Cartesian d-functions (raises NotImplementedError). We use
    # spherical harmonics here. The purpose of these benchmarks is
    # timing comparison, not exact numerical agreement with IQCP.
    # Basis function counts will differ slightly for 6-31G* (e.g.,
    # 18 spherical vs 19 Cartesian for H2O).
    mol.build()
    return mol


def run_scf_and_gradient(mol: gto.Mole, method: str):
    """Run SCF + gradient, return (energy, gradient_array, converged)."""
    if method == "RHF":
        mf = scf.RHF(mol)
    elif method == "LDA":
        mf = dft.RKS(mol)
        mf.xc = "LDA,VWN5"
    elif method == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "b3lyp5"  # VWN5 variant, matching IQCP implementation
    else:
        raise ValueError(f"Unknown method: {method}")

    mf.conv_tol = CONV_TOL
    energy = mf.kernel()

    grad = mf.nuc_grad_method()
    grad.verbose = 0
    g = grad.kernel()

    return energy, g, mf.converged


def get_cpu_info() -> str:
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":")[1].strip()
    except (FileNotFoundError, PermissionError):
        pass
    return platform.processor() or "unknown"


def collect_metadata() -> dict:
    import pyscf
    import scipy

    return {
        "pyscf_version": pyscf.__version__,
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "cpu": get_cpu_info(),
        "date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "n_runs": N_RUNS,
        "conv_tol": CONV_TOL,
        "threading": "single-threaded (lib.num_threads(1))",
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def run_benchmarks() -> list[dict]:
    results = []
    total = len(GRADIENT_CASES)

    for idx, (mol_name, basis, method) in enumerate(GRADIENT_CASES, 1):
        label = f"{mol_name}/{basis}/{method}"
        print(f"  [{idx}/{total}] {label:<30s}", end="", flush=True)

        mol = build_mol(mol_name, basis)

        timings = []
        energy = None
        converged = None
        grad_norm = None

        for i in range(N_RUNS):
            t0 = time.perf_counter()
            e, g, conv = run_scf_and_gradient(mol, method)
            t1 = time.perf_counter()
            timings.append(t1 - t0)
            if i == 0:
                energy = e
                converged = conv
                grad_norm = float(np.linalg.norm(g))

        mean_s = float(np.mean(timings))
        std_s = float(np.std(timings, ddof=1))
        min_s = float(np.min(timings))
        max_s = float(np.max(timings))

        result = {
            "molecule": mol_name,
            "basis": basis,
            "method": method,
            "n_basis": mol.nao_nr(),
            "n_atoms": mol.natm,
            "energy_ha": energy,
            "grad_norm": grad_norm,
            "mean_s": round(mean_s, 6),
            "std_s": round(std_s, 6),
            "min_s": round(min_s, 6),
            "max_s": round(max_s, 6),
            "n_runs": N_RUNS,
            "converged": converged,
        }
        results.append(result)
        print(f"  {mean_s:8.4f} +/- {std_s:.4f} s   |grad| = {grad_norm:.6e}")

    return results


def print_table(results: list[dict]) -> None:
    print()
    print("=" * 110)
    print(
        f"{'Molecule':<8s} {'Basis':<10s} {'Method':<8s} {'N_bas':>5s} "
        f"{'Mean (s)':>10s} {'Std (s)':>10s} {'Min (s)':>10s} {'Max (s)':>10s} "
        f"{'|grad|':>12s} {'Energy (Ha)':>18s}"
    )
    print("-" * 110)
    for r in results:
        print(
            f"{r['molecule']:<8s} {r['basis']:<10s} {r['method']:<8s} "
            f"{r['n_basis']:5d} "
            f"{r['mean_s']:10.4f} {r['std_s']:10.4f} "
            f"{r['min_s']:10.4f} {r['max_s']:10.4f} "
            f"{r['grad_norm']:12.6e} "
            f"{r['energy_ha']:18.10f}"
        )
    print("=" * 110)


def save_results(metadata: dict, results: list[dict]) -> Path:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RESULTS_DIR / "gradient_timing.json"
    data = {"metadata": metadata, "benchmarks": results}
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    return output_path


if __name__ == "__main__":
    print("PySCF Gradient Timing Benchmarks")
    print("=" * 60)

    metadata = collect_metadata()
    print(f"PySCF {metadata['pyscf_version']}  |  "
          f"NumPy {metadata['numpy_version']}  |  "
          f"Python {metadata['python_version']}")
    print(f"CPU: {metadata['cpu']}")
    print(f"Threading: {metadata['threading']}")
    print(f"Runs per combination: {N_RUNS}")
    print()

    print("Running gradient benchmarks...")
    results = run_benchmarks()

    print_table(results)

    output_path = save_results(metadata, results)
    print(f"\nResults saved to: {output_path}")
