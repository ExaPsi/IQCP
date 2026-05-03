#!/usr/bin/env python3
"""
PySCF Geometry Optimization Timing Benchmarks for JCIM M2 Manuscript
=====================================================================

Times full geometry optimizations for selected molecule/method combinations.
Each optimization starts from the IQCP preset geometry and runs to convergence.

Systems
-------
- H2O RHF/STO-3G
- H2O B3LYP/6-31G*

Methodology: 5 runs (optimization is expensive), single-threaded, mean +/- std.

Note: Uses PySCF's built-in geometric optimizer (berny_solver or geomopt).
"""

import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from pyscf import dft, gto, lib, scf

lib.num_threads(1)

# Try geometric optimizer -- PySCF uses geomeTRIC or berny
try:
    from pyscf.geomopt.geometric_solver import optimize as geometric_optimize
    HAS_GEOMETRIC = True
except ImportError:
    HAS_GEOMETRIC = False

try:
    from pyscf.geomopt.berny_solver import optimize as berny_optimize
    HAS_BERNY = True
except ImportError:
    HAS_BERNY = False

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

N_RUNS = 5
CONV_TOL = 1e-12
SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"

# ---------------------------------------------------------------------------
# Molecules (bohr)
# ---------------------------------------------------------------------------

MOLECULES = {
    "H2O": {
        "atoms": "O 0.0 0.0 0.2217282; H 0.0 1.4305447 -0.8869128; H 0.0 -1.4305447 -0.8869128",
        "charge": 0,
        "spin": 0,
    },
}

OPTIMIZATION_CASES = [
    ("H2O", "sto-3g", "RHF"),
    ("H2O", "6-31g*", "B3LYP"),
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
    if basis == "6-31g*":
        mol.cart = True
    mol.build()
    return mol


def run_optimization(mol: gto.Mole, method: str):
    """Run geometry optimization, return (energy, converged, opt_mol)."""
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
    mf.verbose = 0

    # Use geometric if available, fall back to berny, then raise
    if HAS_GEOMETRIC:
        opt_mol = geometric_optimize(mf, verbose=0)
    elif HAS_BERNY:
        opt_mol = berny_optimize(mf, verbose=0)
    else:
        raise RuntimeError(
            "No geometry optimizer available. Install geometric or pyberny: "
            "pip install geometric pyberny"
        )

    # Recompute final energy at optimized geometry
    mf_final = mf.__class__(opt_mol)
    if hasattr(mf, "xc"):
        mf_final.xc = mf.xc
    mf_final.conv_tol = CONV_TOL
    mf_final.verbose = 0
    energy = mf_final.kernel()

    return energy, mf_final.converged, opt_mol


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

    optimizer = "geometric" if HAS_GEOMETRIC else ("berny" if HAS_BERNY else "none")
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
        "optimizer": optimizer,
        "threading": "single-threaded (lib.num_threads(1))",
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def run_benchmarks() -> list[dict]:
    results = []
    total = len(OPTIMIZATION_CASES)

    for idx, (mol_name, basis, method) in enumerate(OPTIMIZATION_CASES, 1):
        label = f"{mol_name}/{basis}/{method}"
        print(f"  [{idx}/{total}] {label:<30s}", end="", flush=True)

        timings = []
        energy = None
        converged = None

        for i in range(N_RUNS):
            mol = build_mol(mol_name, basis)  # fresh mol each run
            t0 = time.perf_counter()
            e, conv, opt_mol = run_optimization(mol, method)
            t1 = time.perf_counter()
            timings.append(t1 - t0)
            if i == 0:
                energy = e
                converged = conv

        mean_s = float(np.mean(timings))
        std_s = float(np.std(timings, ddof=1))
        min_s = float(np.min(timings))
        max_s = float(np.max(timings))

        mol_ref = build_mol(mol_name, basis)
        result = {
            "molecule": mol_name,
            "basis": basis,
            "method": method,
            "n_basis": mol_ref.nao_nr(),
            "energy_ha": energy,
            "mean_s": round(mean_s, 6),
            "std_s": round(std_s, 6),
            "min_s": round(min_s, 6),
            "max_s": round(max_s, 6),
            "n_runs": N_RUNS,
            "converged": converged,
        }
        results.append(result)
        print(f"  {mean_s:8.4f} +/- {std_s:.4f} s   E = {energy:.10f} Ha")

    return results


def print_table(results: list[dict]) -> None:
    print()
    print("=" * 100)
    print(
        f"{'Molecule':<8s} {'Basis':<10s} {'Method':<8s} {'N_bas':>5s} "
        f"{'Mean (s)':>10s} {'Std (s)':>10s} {'Min (s)':>10s} {'Max (s)':>10s} "
        f"{'Energy (Ha)':>18s}"
    )
    print("-" * 100)
    for r in results:
        print(
            f"{r['molecule']:<8s} {r['basis']:<10s} {r['method']:<8s} "
            f"{r['n_basis']:5d} "
            f"{r['mean_s']:10.4f} {r['std_s']:10.4f} "
            f"{r['min_s']:10.4f} {r['max_s']:10.4f} "
            f"{r['energy_ha']:18.10f}"
        )
    print("=" * 100)


def save_results(metadata: dict, results: list[dict]) -> Path:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RESULTS_DIR / "optimization_timing.json"
    data = {"metadata": metadata, "benchmarks": results}
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    return output_path


if __name__ == "__main__":
    print("PySCF Geometry Optimization Timing Benchmarks")
    print("=" * 60)

    if not HAS_GEOMETRIC and not HAS_BERNY:
        print("ERROR: No geometry optimizer available.")
        print("Install one: pip install geometric  OR  pip install pyberny")
        sys.exit(1)

    metadata = collect_metadata()
    print(f"PySCF {metadata['pyscf_version']}  |  "
          f"NumPy {metadata['numpy_version']}  |  "
          f"Python {metadata['python_version']}")
    print(f"CPU: {metadata['cpu']}")
    print(f"Optimizer: {metadata['optimizer']}")
    print(f"Threading: {metadata['threading']}")
    print(f"Runs per combination: {N_RUNS}")
    print()

    print("Running optimization benchmarks...")
    results = run_benchmarks()

    print_table(results)

    output_path = save_results(metadata, results)
    print(f"\nResults saved to: {output_path}")
