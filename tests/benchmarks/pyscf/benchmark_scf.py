#!/usr/bin/env python3
"""
PySCF SCF Timing Benchmarks for JCIM M2 Manuscript
====================================================

Times SCF calculations for 6 molecules x 6 basis sets x 3 methods = 108 combinations.
Each combination is run 10 times; mean, std, min, max wall-clock times are reported.

Methodology
-----------
- Single-threaded execution (lib.num_threads(1)) for fair comparison with
  single-threaded WASM execution in IQCP.
- Geometries, basis sets, and convergence thresholds match IQCP validation data
  (tests/golden/dft/validation.json).
- Coordinates are in bohr; 6-31G* uses Cartesian d-functions (mol.cart = True).
- conv_tol = 1e-10 throughout (matching IQCP's ScfConfig::tight() energy threshold).

Output
------
- Human-readable table to stdout
- JSON results to tests/benchmarks/pyscf/results/scf_timing.json
"""

import json
import os
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
# Match IQCP's ScfConfig::tight() energy threshold for fair comparison
CONV_TOL = 1e-10
SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"

# Force single-threaded execution for fair WASM comparison
lib.num_threads(1)

# ---------------------------------------------------------------------------
# Molecule definitions (coordinates in bohr, matching validation.json)
# ---------------------------------------------------------------------------

MOLECULES = {
    "H2": {
        "atoms": "H 0 0 0; H 0 0 1.4",
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
    "HF": {
        "atoms": "H 0 0 0; F 0 0 1.7328",
        "charge": 0,
        "spin": 0,
    },
    "CH4": {
        "atoms": "C 0 0 0; H 1.1851 1.1851 1.1851; H -1.1851 -1.1851 1.1851; H -1.1851 1.1851 -1.1851; H 1.1851 -1.1851 -1.1851",
        "charge": 0,
        "spin": 0,
    },
    "C6H6": {
        "atoms": "C 0 1.397 0; C 1.2098 0.6985 0; C 1.2098 -0.6985 0; C 0 -1.397 0; C -1.2098 -0.6985 0; C -1.2098 0.6985 0; H 0 2.481 0; H 2.1486 1.2405 0; H 2.1486 -1.2405 0; H 0 -2.481 0; H -2.1486 -1.2405 0; H -2.1486 1.2405 0",
        "charge": 0,
        "spin": 0,
        "unit": "Angstrom",
    },
}

BASIS_SETS = ["sto-3g", "3-21g", "6-31g", "6-31g*", "6-31+g*", "cc-pvdz"]
METHODS = ["RHF", "LDA", "B3LYP"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def build_mol(mol_name: str, basis: str) -> gto.Mole:
    """Build a PySCF Mole object with IQCP-consistent parameters."""
    mol_def = MOLECULES[mol_name]
    mol = gto.Mole()
    mol.atom = mol_def["atoms"]
    mol.basis = basis
    mol.charge = mol_def["charge"]
    mol.spin = mol_def["spin"]
    mol.unit = mol_def.get("unit", "Bohr")
    mol.verbose = 0  # suppress output during benchmarks
    if basis in ("6-31g*", "6-31+g*", "cc-pvdz"):
        mol.cart = True  # Cartesian d-functions to match IQCP
    mol.build()
    return mol


def run_scf(mol: gto.Mole, method: str, precomputed_eri=None):
    """Run a single SCF calculation and return (energy, converged, mf, cycles)."""
    if method == "RHF":
        mf = scf.RHF(mol)
        if precomputed_eri is not None:
            mf._eri = precomputed_eri
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
    return energy, mf.converged, mf, getattr(mf, 'cycles', None)


def get_cpu_info() -> str:
    """Best-effort CPU identification."""
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":")[1].strip()
    except (FileNotFoundError, PermissionError):
        pass
    return platform.processor() or "unknown"


def collect_metadata() -> dict:
    """Collect environment metadata for reproducibility."""
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
# Main benchmark loop
# ---------------------------------------------------------------------------


def run_benchmarks() -> list[dict]:
    """Run all 30 benchmark combinations and return results."""
    results = []
    total = len(MOLECULES) * len(BASIS_SETS) * len(METHODS)
    count = 0

    for mol_name in MOLECULES:
        for basis in BASIS_SETS:
            mol = build_mol(mol_name, basis)
            n_basis = mol.nao_nr()
            n_electrons = mol.nelectron

            for method in METHODS:
                count += 1
                label = f"{mol_name}/{basis}/{method}"
                print(f"  [{count:2d}/{total}] {label:<30s}", end="", flush=True)

                # Pre-compute ERIs for RHF (matching IQCP benchmark which uses pre-stored integrals)
                precomputed_eri = mol.intor('int2e') if method == "RHF" else None

                timings = []
                energy = None
                converged = None
                cycles = None

                for i in range(N_RUNS):
                    t0 = time.perf_counter()
                    e, conv, _, cyc = run_scf(mol, method, precomputed_eri)
                    t1 = time.perf_counter()
                    timings.append(t1 - t0)
                    if i == 0:
                        energy = e
                        converged = conv
                        cycles = cyc

                mean_s = float(np.mean(timings))
                std_s = float(np.std(timings, ddof=1))
                min_s = float(np.min(timings))
                max_s = float(np.max(timings))

                result = {
                    "molecule": mol_name,
                    "basis": basis,
                    "method": method,
                    "n_basis": n_basis,
                    "n_electrons": n_electrons,
                    "energy_ha": energy,
                    "iterations": cycles,
                    "mean_s": round(mean_s, 6),
                    "std_s": round(std_s, 6),
                    "min_s": round(min_s, 6),
                    "max_s": round(max_s, 6),
                    "n_runs": N_RUNS,
                    "converged": converged,
                }
                results.append(result)
                print(
                    f"  {mean_s:8.4f} +/- {std_s:.4f} s   E = {energy:.10f} Ha"
                )

    return results


def print_table(results: list[dict]) -> None:
    """Print a human-readable summary table."""
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
    """Save results to JSON."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RESULTS_DIR / "scf_timing.json"
    data = {
        "metadata": metadata,
        "benchmarks": results,
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    return output_path


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("PySCF SCF Timing Benchmarks")
    print("=" * 60)

    metadata = collect_metadata()
    print(f"PySCF {metadata['pyscf_version']}  |  "
          f"NumPy {metadata['numpy_version']}  |  "
          f"Python {metadata['python_version']}")
    print(f"CPU: {metadata['cpu']}")
    print(f"Threading: {metadata['threading']}")
    print(f"Runs per combination: {N_RUNS}")
    print(f"Convergence tolerance: {CONV_TOL}")
    print()

    print("Running benchmarks...")
    results = run_benchmarks()

    print_table(results)

    output_path = save_results(metadata, results)
    print(f"\nResults saved to: {output_path}")
