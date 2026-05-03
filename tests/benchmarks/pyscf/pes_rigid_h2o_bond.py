#!/usr/bin/env python3
"""Generate PySCF golden data for H2O rigid O-H bond scan (RHF/STO-3G).

H2O preset geometry (bohr):
  O  [0.0,  0.0,        0.2217282]
  H1 [0.0,  1.4305447, -0.8869128]
  H2 [0.0, -1.4305447, -0.8869128]

The scan varies the O-H1 bond length from 1.4 to 2.4 bohr (11 points).
At each bond length:
  - O remains at its position
  - H1 is moved along the O->H1 direction to the target distance
  - H2 and the H-O-H angle remain fixed

This matches IQCP's generate_bond_geometry(atoms, 0, 1, target_r):
  - atom_i = O (fixed anchor at index 0)
  - atom_j = H1 (moved to target distance, index 1)

Output: tests/golden/pes/h2o_bond_rigid_sto3g.json
"""

import json
import sys
import numpy as np
from datetime import datetime, timezone

import pyscf
from pyscf import gto, scf

print(f"PySCF version: {pyscf.__version__}")
print(f"NumPy version: {np.__version__}")

# ---------------------------------------------------------------------------
# Reference geometry (bohr)
# ---------------------------------------------------------------------------
O_pos  = np.array([0.0,  0.0,        0.2217282])
H1_pos = np.array([0.0,  1.4305447, -0.8869128])
H2_pos = np.array([0.0, -1.4305447, -0.8869128])

# O-H1 bond direction (unit vector)
oh1_vec = H1_pos - O_pos
R_OH1_init = np.linalg.norm(oh1_vec)
oh1_dir = oh1_vec / R_OH1_init

print(f"R_OH1 initial = {R_OH1_init:.10f} bohr")
print(f"O->H1 direction = {oh1_dir}")

# Verify initial O-H2 distance
R_OH2 = np.linalg.norm(H2_pos - O_pos)
print(f"R_OH2 = {R_OH2:.10f} bohr (unchanged throughout scan)")

# ---------------------------------------------------------------------------
# Scan parameters
# ---------------------------------------------------------------------------
r_min = 1.4
r_max = 2.4
n_points = 11

distances = np.linspace(r_min, r_max, n_points)

print(f"\nScanning O-H1 bond from {r_min} to {r_max} bohr, {n_points} points")
print("-" * 70)

# ---------------------------------------------------------------------------
# Run scan
# ---------------------------------------------------------------------------
points = []

for i, r in enumerate(distances):
    # Generate geometry: move H1 along the O->H1 direction to distance r
    H1_new = O_pos + oh1_dir * r

    # Verify bond length
    r_check = np.linalg.norm(H1_new - O_pos)

    # Build PySCF molecule (input in bohr)
    mol = gto.Mole()
    mol.unit = 'bohr'
    mol.atom = [
        ['O',  O_pos.tolist()],
        ['H',  H1_new.tolist()],
        ['H',  H2_pos.tolist()],
    ]
    mol.basis = 'sto-3g'
    mol.verbose = 0
    mol.build()

    # Run RHF
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 100
    energy = mf.kernel()

    print(f"  [{i+1:2d}/{n_points}]  R_OH1 = {r:6.3f} bohr  "
          f"E = {energy:.12f} Ha  converged = {mf.converged}  "
          f"cycles = {mf.cycles}  r_check = {r_check:.6f}")

    points.append({
        "coordinate_value": float(r),
        "energy": float(energy),
        "converged": bool(mf.converged),
        "scf_iterations": int(mf.cycles),
        "geometry": [
            O_pos.tolist(),
            H1_new.tolist(),
            H2_pos.tolist(),
        ],
    })

# ---------------------------------------------------------------------------
# Find equilibrium via parabolic interpolation
# ---------------------------------------------------------------------------
energies = np.array([p["energy"] for p in points])
coords = np.array([p["coordinate_value"] for p in points])

min_idx = np.argmin(energies)
if 0 < min_idx < len(energies) - 1:
    v1, v2, v3 = coords[min_idx-1], coords[min_idx], coords[min_idx+1]
    e1, e2, e3 = energies[min_idx-1], energies[min_idx], energies[min_idx+1]
    denom = (v1 - v2) * (v1 - v3) * (v2 - v3)
    a = (v3*(e2-e1) + v2*(e1-e3) + v1*(e3-e2)) / denom
    b = (v3**2*(e1-e2) + v2**2*(e3-e1) + v1**2*(e2-e3)) / denom
    c = (v2*v3*(v2-v3)*e1 + v3*v1*(v3-v1)*e2 + v1*v2*(v1-v2)*e3) / denom
    v_eq = -b / (2*a)
    e_eq = a*v_eq**2 + b*v_eq + c
    equilibrium = {
        "value": float(v_eq),
        "energy": float(e_eq),
    }
    print(f"\nEquilibrium: R_OH1 = {v_eq:.6f} bohr, E = {e_eq:.12f} Ha")
else:
    equilibrium = {
        "value": float(coords[min_idx]),
        "energy": float(energies[min_idx]),
    }
    print(f"\nMinimum at boundary: R_OH1 = {coords[min_idx]:.6f} bohr")

# ---------------------------------------------------------------------------
# Save JSON
# ---------------------------------------------------------------------------
result = {
    "molecule": "H2O",
    "basis": "sto-3g",
    "method": "RHF",
    "scan_type": "bond_rigid",
    "scan_coordinate": {
        "type": "bond",
        "atom_indices": [0, 1],
        "description": "O-H1 bond length (O is anchor, H1 is moved)",
    },
    "reference_geometry_bohr": {
        "O": O_pos.tolist(),
        "H1": H1_pos.tolist(),
        "H2": H2_pos.tolist(),
        "R_OH1_init_bohr": float(R_OH1_init),
        "R_OH2_bohr": float(R_OH2),
    },
    "scan_range": {
        "min_bohr": float(r_min),
        "max_bohr": float(r_max),
        "n_points": n_points,
    },
    "points": points,
    "equilibrium": equilibrium,
    "provenance": {
        "generated": datetime.now(timezone.utc).isoformat(),
        "pyscf_version": pyscf.__version__,
        "numpy_version": np.__version__,
        "conv_tol": 1e-12,
        "max_cycle": 100,
    },
}

output_path = "tests/golden/pes/h2o_bond_rigid_sto3g.json"
with open(output_path, "w") as f:
    json.dump(result, f, indent=2)

print(f"\nSaved to {output_path}")
