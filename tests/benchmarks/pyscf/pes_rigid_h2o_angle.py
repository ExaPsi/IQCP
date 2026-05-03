#!/usr/bin/env python3
"""Generate PySCF golden data for H2O rigid angle scan (RHF/STO-3G).

H2O preset geometry (bohr):
  O  [0.0,  0.0,        0.2217282]
  H1 [0.0,  1.4305447, -0.8869128]
  H2 [0.0, -1.4305447, -0.8869128]

The scan varies the H-O-H angle from 80 to 130 degrees (11 points).
At each angle:
  - O remains at its position
  - H1 at: O + R_OH * [0,  sin(theta/2), -cos(theta/2)]
  - H2 at: O + R_OH * [0, -sin(theta/2), -cos(theta/2)]
where R_OH is the initial O-H distance and theta is the full H-O-H angle.

This geometry generation must match IQCP's generate_angle_geometry().
IQCP rotates atom_k (H2) around the normal to the H1-O-H2 plane by
delta_theta = target - current. Both approaches produce geometries with
the SAME H-O-H angle and O-H distances. The absolute Cartesian coords
may differ, but the internal coordinates (and thus energy) are identical.

Output: tests/golden/pes/h2o_angle_rigid_sto3g.json
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

R_OH = np.linalg.norm(H1_pos - O_pos)
print(f"R_OH = {R_OH:.10f} bohr")

# Verify initial angle
v1 = H1_pos - O_pos
v2 = H2_pos - O_pos
cos_init = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
theta_init = np.degrees(np.arccos(np.clip(cos_init, -1, 1)))
print(f"Initial H-O-H angle = {theta_init:.6f} degrees")

# ---------------------------------------------------------------------------
# Scan parameters
# ---------------------------------------------------------------------------
angle_min_deg = 80.0
angle_max_deg = 130.0
n_points = 11

angles_deg = np.linspace(angle_min_deg, angle_max_deg, n_points)
angles_rad = np.radians(angles_deg)

print(f"\nScanning H-O-H angle from {angle_min_deg} to {angle_max_deg} degrees, {n_points} points")
print("-" * 70)

# ---------------------------------------------------------------------------
# Run scan
# ---------------------------------------------------------------------------
points = []

for i, (theta_deg, theta_rad) in enumerate(zip(angles_deg, angles_rad)):
    # Generate geometry: place H atoms symmetrically
    half_theta = theta_rad / 2.0
    H1_new = O_pos + R_OH * np.array([0.0,  np.sin(half_theta), -np.cos(half_theta)])
    H2_new = O_pos + R_OH * np.array([0.0, -np.sin(half_theta), -np.cos(half_theta)])

    # Verify angle
    vv1 = H1_new - O_pos
    vv2 = H2_new - O_pos
    cos_check = np.dot(vv1, vv2) / (np.linalg.norm(vv1) * np.linalg.norm(vv2))
    angle_check = np.degrees(np.arccos(np.clip(cos_check, -1, 1)))

    # Build PySCF molecule (input in bohr)
    mol = gto.Mole()
    mol.unit = 'bohr'
    mol.atom = [
        ['O',  O_pos.tolist()],
        ['H',  H1_new.tolist()],
        ['H',  H2_new.tolist()],
    ]
    mol.basis = 'sto-3g'
    mol.verbose = 0
    mol.build()

    # Run RHF
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 100
    energy = mf.kernel()

    print(f"  [{i+1:2d}/{n_points}]  angle = {theta_deg:7.2f} deg  "
          f"E = {energy:.12f} Ha  converged = {mf.converged}  "
          f"cycles = {mf.cycles}  check_angle = {angle_check:.4f} deg")

    points.append({
        "coordinate_value": float(theta_rad),
        "coordinate_value_deg": float(theta_deg),
        "energy": float(energy),
        "converged": bool(mf.converged),
        "scf_iterations": int(mf.cycles),
        "geometry": [
            O_pos.tolist(),
            H1_new.tolist(),
            H2_new.tolist(),
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
        "value_deg": float(np.degrees(v_eq)),
        "energy": float(e_eq),
    }
    print(f"\nEquilibrium: angle = {np.degrees(v_eq):.4f} deg, E = {e_eq:.12f} Ha")
else:
    equilibrium = {
        "value": float(coords[min_idx]),
        "value_deg": float(np.degrees(coords[min_idx])),
        "energy": float(energies[min_idx]),
    }
    print(f"\nMinimum at boundary: angle = {np.degrees(coords[min_idx]):.4f} deg")

# ---------------------------------------------------------------------------
# Save JSON
# ---------------------------------------------------------------------------
result = {
    "molecule": "H2O",
    "basis": "sto-3g",
    "method": "RHF",
    "scan_type": "angle_rigid",
    "scan_coordinate": {
        "type": "angle",
        "atom_indices": [1, 0, 2],
        "description": "H1-O-H2 angle (O is central atom)",
    },
    "reference_geometry_bohr": {
        "O": O_pos.tolist(),
        "H1": H1_pos.tolist(),
        "H2": H2_pos.tolist(),
        "R_OH_bohr": float(R_OH),
    },
    "scan_range": {
        "min_rad": float(angles_rad[0]),
        "max_rad": float(angles_rad[-1]),
        "min_deg": float(angle_min_deg),
        "max_deg": float(angle_max_deg),
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

output_path = "tests/golden/pes/h2o_angle_rigid_sto3g.json"
with open(output_path, "w") as f:
    json.dump(result, f, indent=2)

print(f"\nSaved to {output_path}")
