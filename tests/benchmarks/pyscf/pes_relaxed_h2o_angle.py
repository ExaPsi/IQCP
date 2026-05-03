#!/usr/bin/env python3
"""Generate PySCF relaxed PES scan for H2O angle (geomeTRIC constrained optimization).

Generates golden reference data for the relaxed H-O-H angle PES scan by running
geomeTRIC constrained optimizations at each target angle. At each scan point,
the H-O-H angle is held fixed while the O-H bond lengths are optimized.

This produces the minimum energy path (MEP) along the angle coordinate, which
serves as a reference for IQCP's relaxed PES scan implementation.

Usage:
    source .venv/bin/activate
    python tests/benchmarks/pyscf/pes_relaxed_h2o_angle.py

Environment:
    PySCF 2.11.0
    geomeTRIC 1.1
    NumPy 2.4.1
    SciPy 1.17.0

Reference geometry:
    IQCP H2O preset (bohr):
      O  [0.0, 0.0, 0.2217282]
      H1 [0.0, 1.4305447, -0.8869128]
      H2 [0.0, -1.4305447, -0.8869128]
    Initial angle: ~104.5 degrees
    Initial R_OH: ~1.8098 bohr

Output:
    tests/golden/pes/h2o_angle_relaxed_sto3g.json
"""
import json
import os
import sys
import tempfile
from datetime import datetime, timezone

import numpy as np
import pyscf
from pyscf import gto, scf
from pyscf.geomopt.geometric_solver import optimize

# =============================================================================
# Configuration
# =============================================================================

# IQCP H2O preset geometry (bohr)
O  = [0.0, 0.0, 0.2217282]
H1 = [0.0, 1.4305447, -0.8869128]
H2 = [0.0, -1.4305447, -0.8869128]

# Scan range: 80 to 130 degrees, 11 points (matches rigid scan)
ANGLES_DEG = np.linspace(80.0, 130.0, 11)

# Convergence settings
SCF_CONV_TOL = 1e-12
SCF_MAX_CYCLE = 100
OPT_MAX_STEPS = 50

# Output path
OUTPUT_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    '..', '..', 'golden', 'pes', 'h2o_angle_relaxed_sto3g.json'
)


def compute_angle(coords, i, j, k):
    """Compute bond angle i-j-k in degrees (j is central)."""
    v1 = coords[i] - coords[j]
    v2 = coords[k] - coords[j]
    cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0)))


def compute_bond(coords, i, j):
    """Compute distance between atoms i and j."""
    return np.linalg.norm(coords[i] - coords[j])


def build_starting_geometry(target_angle_deg):
    """Build a starting geometry for H2O at the target H-O-H angle.

    Places both H atoms symmetrically around the O atom to create
    the specified angle while preserving the original O-H bond lengths.
    This is the same approach used in the rigid scan PySCF script.
    """
    # Compute original O-H bond length
    r_oh = np.linalg.norm(np.array(H1) - np.array(O))
    theta = np.radians(target_angle_deg)

    # Place H atoms symmetrically around O
    h1_new = [
        O[0],
        O[1] + r_oh * np.sin(theta / 2.0),
        O[2] - r_oh * np.cos(theta / 2.0),
    ]
    h2_new = [
        O[0],
        O[1] - r_oh * np.sin(theta / 2.0),
        O[2] - r_oh * np.cos(theta / 2.0),
    ]

    return O, h1_new, h2_new


def run_constrained_optimization(target_angle_deg):
    """Run geomeTRIC constrained optimization at the target angle.

    Returns:
        dict with keys: energy, converged, geometry, opt_steps,
                        r_oh1, r_oh2, angle_deg
    """
    o_pos, h1_pos, h2_pos = build_starting_geometry(target_angle_deg)

    # Build PySCF molecule
    mol = gto.Mole()
    mol.unit = 'bohr'
    mol.atom = (
        f'O  {o_pos[0]:.10f} {o_pos[1]:.10f} {o_pos[2]:.10f}; '
        f'H  {h1_pos[0]:.10f} {h1_pos[1]:.10f} {h1_pos[2]:.10f}; '
        f'H  {h2_pos[0]:.10f} {h2_pos[1]:.10f} {h2_pos[2]:.10f}'
    )
    mol.basis = 'sto-3g'
    mol.verbose = 0
    mol.build()

    # Run initial SCF to get a starting wavefunction
    mf = scf.RHF(mol)
    mf.conv_tol = SCF_CONV_TOL
    mf.max_cycle = SCF_MAX_CYCLE
    mf.kernel()

    # Write geomeTRIC constraint file
    # NOTE: geomeTRIC uses 1-based atom indices
    constraint_text = f'$set\nangle 2 1 3 {target_angle_deg:.6f}\n'

    constraint_fd, constraint_path = tempfile.mkstemp(suffix='.txt')
    try:
        with os.fdopen(constraint_fd, 'w') as f:
            f.write(constraint_text)

        # Run constrained optimization
        mol_opt = optimize(
            mf,
            maxsteps=OPT_MAX_STEPS,
            constraints=constraint_path,
        )

        # Extract optimized geometry
        coords = mol_opt.atom_coords(unit='Bohr')

        # Run a final high-precision SCF at the optimized geometry
        mf_opt = scf.RHF(mol_opt)
        mf_opt.conv_tol = SCF_CONV_TOL
        mf_opt.max_cycle = SCF_MAX_CYCLE
        e_relaxed = mf_opt.kernel()

        # Compute geometric parameters
        angle_final = compute_angle(coords, 1, 0, 2)
        r_oh1 = compute_bond(coords, 0, 1)
        r_oh2 = compute_bond(coords, 0, 2)

        return {
            'energy': float(e_relaxed),
            'converged': True,
            'geometry': coords.tolist(),
            'r_oh1_bohr': float(r_oh1),
            'r_oh2_bohr': float(r_oh2),
            'angle_deg': float(angle_final),
        }

    except Exception as e:
        print(f'    geomeTRIC failed at {target_angle_deg:.1f} deg: {e}',
              file=sys.stderr)
        # Fall back to rigid geometry energy
        return {
            'energy': float(mf.e_tot),
            'converged': False,
            'geometry': [list(o_pos), list(h1_pos), list(h2_pos)],
            'r_oh1_bohr': float(np.linalg.norm(np.array(h1_pos) - np.array(o_pos))),
            'r_oh2_bohr': float(np.linalg.norm(np.array(h2_pos) - np.array(o_pos))),
            'angle_deg': float(target_angle_deg),
        }
    finally:
        os.unlink(constraint_path)


def main():
    print(f'PySCF version: {pyscf.__version__}')
    print(f'NumPy version: {np.__version__}')

    import geometric
    print(f'geomeTRIC version: {geometric.__version__}')
    print()

    # Original geometry info
    r_oh_orig = np.linalg.norm(np.array(H1) - np.array(O))
    angle_orig = compute_angle(
        np.array([O, H1, H2]), 1, 0, 2
    )
    print(f'Reference geometry: R_OH = {r_oh_orig:.8f} bohr, angle = {angle_orig:.4f} deg')
    print(f'Scan range: {ANGLES_DEG[0]:.1f} to {ANGLES_DEG[-1]:.1f} degrees, {len(ANGLES_DEG)} points')
    print(f'Method: RHF/STO-3G, conv_tol={SCF_CONV_TOL}')
    print()

    results = []
    for i, angle_deg in enumerate(ANGLES_DEG):
        print(f'Point {i+1}/{len(ANGLES_DEG)}: angle = {angle_deg:.1f} deg ...')
        result = run_constrained_optimization(angle_deg)
        result['coordinate_value'] = float(np.radians(angle_deg))
        result['coordinate_value_deg'] = float(angle_deg)
        results.append(result)

        if result['converged']:
            print(f'  E = {result["energy"]:.10f} Ha, '
                  f'R_OH1 = {result["r_oh1_bohr"]:.6f}, '
                  f'R_OH2 = {result["r_oh2_bohr"]:.6f}, '
                  f'angle = {result["angle_deg"]:.4f} deg')
        else:
            print(f'  FAILED: using rigid energy {result["energy"]:.10f} Ha')

    # Find equilibrium via parabolic interpolation
    converged_pts = [(r['coordinate_value'], r['energy'])
                     for r in results if r['converged']]
    equilibrium = None
    if len(converged_pts) >= 3:
        vals = [v for v, e in converged_pts]
        energies = [e for v, e in converged_pts]
        min_idx = np.argmin(energies)
        if 0 < min_idx < len(energies) - 1:
            v1, v2, v3 = vals[min_idx-1], vals[min_idx], vals[min_idx+1]
            e1, e2, e3 = energies[min_idx-1], energies[min_idx], energies[min_idx+1]
            denom = (v1 - v2) * (v1 - v3) * (v2 - v3)
            if abs(denom) > 1e-30:
                a = (v3 * (e2 - e1) + v2 * (e1 - e3) + v1 * (e3 - e2)) / denom
                b = (v3**2 * (e1 - e2) + v2**2 * (e3 - e1) + v1**2 * (e2 - e3)) / denom
                c = (v2*v3*(v2-v3)*e1 + v3*v1*(v3-v1)*e2 + v1*v2*(v1-v2)*e3) / denom
                if a > 0:
                    v_eq = -b / (2.0 * a)
                    e_eq = a * v_eq**2 + b * v_eq + c
                    equilibrium = {
                        'value': float(v_eq),
                        'value_deg': float(np.degrees(v_eq)),
                        'energy': float(e_eq),
                    }

    # Build output JSON
    output = {
        'molecule': 'H2O',
        'basis': 'sto-3g',
        'method': 'RHF',
        'scan_type': 'angle_relaxed',
        'scan_coordinate': {
            'type': 'angle',
            'atom_indices': [1, 0, 2],
            'description': 'H1-O-H2 angle (O is central atom), relaxed O-H bonds',
        },
        'reference_geometry_bohr': {
            'O': O,
            'H1': H1,
            'H2': H2,
            'R_OH_bohr': float(r_oh_orig),
        },
        'scan_range': {
            'min_rad': float(np.radians(ANGLES_DEG[0])),
            'max_rad': float(np.radians(ANGLES_DEG[-1])),
            'min_deg': float(ANGLES_DEG[0]),
            'max_deg': float(ANGLES_DEG[-1]),
            'n_points': len(ANGLES_DEG),
        },
        'provenance': {
            'generated': datetime.now(timezone.utc).isoformat(),
            'pyscf_version': pyscf.__version__,
            'numpy_version': np.__version__,
            'geometric_version': geometric.__version__,
            'scf_conv_tol': SCF_CONV_TOL,
            'scf_max_cycle': SCF_MAX_CYCLE,
            'opt_max_steps': OPT_MAX_STEPS,
            'script': 'tests/benchmarks/pyscf/pes_relaxed_h2o_angle.py',
        },
        'points': results,
        'equilibrium': equilibrium,
    }

    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(OUTPUT_PATH)), exist_ok=True)

    with open(OUTPUT_PATH, 'w') as f:
        json.dump(output, f, indent=2)

    print()
    print(f'Golden data written to: {OUTPUT_PATH}')
    print(f'  {len(results)} points, {sum(1 for r in results if r["converged"])} converged')
    if equilibrium:
        print(f'  Equilibrium: {equilibrium["value_deg"]:.2f} deg, '
              f'E = {equilibrium["energy"]:.10f} Ha')

    # Summary table
    print()
    print('Summary:')
    print(f'{"Point":>5} {"Angle (deg)":>12} {"Energy (Ha)":>18} {"R_OH1":>10} {"R_OH2":>10} {"Conv":>6}')
    print('-' * 70)
    for i, r in enumerate(results):
        conv_str = 'YES' if r['converged'] else 'NO'
        print(f'{i+1:5d} {r["coordinate_value_deg"]:12.1f} {r["energy"]:18.10f} '
              f'{r["r_oh1_bohr"]:10.6f} {r["r_oh2_bohr"]:10.6f} {conv_str:>6}')


if __name__ == '__main__':
    main()
