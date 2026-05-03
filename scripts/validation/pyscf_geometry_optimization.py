#!/usr/bin/env python3
"""
PySCF/geomeTRIC geometry optimization for the four molecules cited in
M2 manuscript §4.5 (Geometry Optimization Validation).

Usage:
    source .venv/bin/activate
    python scripts/validation/pyscf_geometry_optimization.py

Reports (per molecule, RHF/STO-3G):
- final energy
- step count (number of evaluations after the initial point)
- canonical bond length

Reference: PySCF 2.11.0, geomeTRIC 1.1, conv_tol = 1e-12.
geomeTRIC default convergence: |dE| < 1e-6, RMS-Grad < 3e-4,
Max-Grad < 4.5e-4, RMS-Disp < 1.2e-3, Max-Disp < 1.8e-3.
"""
import json
from pathlib import Path

import numpy as np
from pyscf import gto, scf
from pyscf.geomopt import geometric_solver

SYSTEMS = {
    "H2O": (
        "O 0.0 0.0  0.117370; "
        "H 0.0 0.756950 -0.469483; "
        "H 0.0 -0.756950 -0.469483",
        "Angstrom",
        ("O", "H"),
    ),
    "NH3": (
        "N 0.0 0.0 0.1173; "
        "H 0.0 0.9377 -0.2737; "
        "H 0.8120 -0.4689 -0.2737; "
        "H -0.8120 -0.4689 -0.2737",
        "Angstrom",
        ("N", "H"),
    ),
    "CH4": (
        "C 0.0 0.0 0.0; "
        "H 0.6276 0.6276 0.6276; "
        "H -0.6276 -0.6276 0.6276; "
        "H -0.6276 0.6276 -0.6276; "
        "H 0.6276 -0.6276 -0.6276",
        "Angstrom",
        ("C", "H"),
    ),
    "HF": (
        "H 0.0 0.0 0.0; F 0.0 0.0 0.917",
        "Angstrom",
        ("H", "F"),
    ),
}


def first_pair_distance(coords, atoms, a, b):
    ia = atoms.index(a)
    ib = next(i for i, sym in enumerate(atoms) if sym == b and i != ia)
    return float(np.linalg.norm(coords[ia] - coords[ib]))


def main() -> None:
    results = {}
    for name, (atom_str, unit, bond) in SYSTEMS.items():
        mol = gto.Mole()
        mol.atom = atom_str
        mol.unit = unit
        mol.basis = "sto-3g"
        mol.verbose = 0
        mol.build()

        mf = scf.RHF(mol)
        mf.conv_tol = 1e-12

        # geomeTRIC writes step trajectory to its own log; we count points
        # by collecting energies from the converged callback.
        history = []

        def _record(envs):
            history.append(envs["energy"])

        mol_opt = geometric_solver.optimize(
            mf, callback=_record, maxsteps=100
        )

        coords = mol_opt.atom_coords()
        atoms = [mol_opt.atom_symbol(i) for i in range(mol_opt.natm)]
        bond_len_bohr = first_pair_distance(coords, atoms, *bond)

        results[name] = {
            "method": "RHF",
            "basis": "STO-3G",
            "n_points_evaluated": len(history),
            "n_optimization_steps": max(0, len(history) - 1),
            "converged": mf.converged,
            "energy_initial": history[0] if history else None,
            "energy_final": history[-1] if history else None,
            f"d({bond[0]}-{bond[1]})_bohr": bond_len_bohr,
            f"d({bond[0]}-{bond[1]})_Angstrom": bond_len_bohr * 0.529177,
        }

    out = Path(__file__).resolve().parents[2] / "docs/Manuscripts/M2/Data/optimization_validation_pyscf.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2))
    print(f"Wrote {out}")
    for name, r in results.items():
        print(
            f"{name:5s}  steps={r['n_optimization_steps']:2d}  "
            f"E={r['energy_final']:.10f}  "
            f"d={r[next(k for k in r if k.endswith('_Angstrom'))]:.4f} A"
        )


if __name__ == "__main__":
    main()
