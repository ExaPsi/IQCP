#!/usr/bin/env python3
"""Generate PySCF golden reference data for US-096 (Normal Mode Analysis).

Runs PySCF's RHF/Hessian/thermo pipeline on H2O, CH4, CO2, and H2 at
STO-3G with mol.cart=True, then dumps:

  - atoms (symbols, coords in bohr)
  - hessian (3N x 3N in Hartree/bohr^2)
  - freq_wavenumber (cm^-1)
  - freq_au (atomic units)
  - reduced_mass (amu)
  - force_const_dyne (mDyne/A)
  - force_const_au (Ha/bohr^2)
  - norm_mode (Cartesian, [mode][atom][dir])
  - rot_type (string as per PySCF _get_rotor_type)
  - rotational_constants_ghz
  - principal_moments_amu_bohr2

to JSON files under tests/golden/thermo/.

Also prints all values to stdout so they can be captured in the
verification report.

Usage (from ):
    source .venv/bin/activate
    python scripts/phase5/generate_thermo_golden.py
"""

import json
import os
from pathlib import Path

import numpy as np
import pyscf
from pyscf import gto, scf
from pyscf import hessian as hess_mod
from pyscf.hessian import thermo as pyscf_thermo


OUT_DIR = Path(__file__).resolve().parents[2] / "tests" / "golden" / "thermo"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_molecule(name, atom_str, basis="sto-3g"):
    print()
    print("=" * 78)
    print(f"=== {name} / {basis} / RHF ===")
    print("=" * 78)

    mol = gto.Mole()
    mol.atom = atom_str
    mol.basis = basis
    mol.cart = True
    mol.unit = "angstrom"
    mol.build()

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    e_tot = mf.kernel()
    print(f"SCF energy: {e_tot:.12f} Ha")
    assert mf.converged, f"SCF did not converge for {name}"

    # Analytical Hessian
    hobj = hess_mod.RHF(mf)
    hessian = hobj.kernel()  # shape (natm, natm, 3, 3)
    natm = mol.natm
    # Flatten to (3N, 3N) in layout H[3A+d, 3B+e]
    # PySCF returns shape (natm, natm, 3, 3) where [A, B, d, e] = d^2 E / dR_Ad dR_Be
    h_flat = hessian.transpose(0, 2, 1, 3).reshape(3 * natm, 3 * natm)

    # Normal mode analysis
    freq_info = pyscf_thermo.harmonic_analysis(mol, hessian)
    # Sort frequencies ascending (real part)
    freq_cm1 = np.asarray(freq_info["freq_wavenumber"]).real
    freq_au = np.asarray(freq_info["freq_au"]).real
    reduced_mass = np.asarray(freq_info["reduced_mass"]).real
    force_const_dyne = np.asarray(freq_info["force_const_dyne"]).real
    force_const_au = np.asarray(freq_info["force_const_au"]).real
    # norm_mode is shape (n_modes, n_atoms, 3)
    norm_mode = np.asarray(freq_info["norm_mode"]).real

    sort_idx = np.argsort(freq_cm1)
    freq_cm1 = freq_cm1[sort_idx]
    freq_au = freq_au[sort_idx]
    reduced_mass = reduced_mass[sort_idx]
    force_const_dyne = force_const_dyne[sort_idx]
    force_const_au = force_const_au[sort_idx]
    norm_mode = norm_mode[sort_idx]

    # Rotor type + rotational constants
    atom_mass = mol.atom_mass_list(isotope_avg=True)
    atom_coords = mol.atom_coords()  # bohr
    rot_const = pyscf_thermo.rotation_const(atom_mass, atom_coords, "GHz")
    rotor_type = pyscf_thermo._get_rotor_type(rot_const)

    # Principal moments (amu * bohr^2)
    mass_center = np.einsum("z,zx->x", atom_mass, atom_coords) / atom_mass.sum()
    coords_com = atom_coords - mass_center
    im = np.einsum("m,mx,my->xy", atom_mass, coords_com, coords_com)
    im = np.eye(3) * im.trace() - im
    moments = np.sort(np.linalg.eigvalsh(im))

    print(f"n_atoms = {natm}")
    print(f"n_modes = {len(freq_cm1)}")
    print(f"rot_type = {rotor_type}")
    print(f"rot_const (GHz) = {rot_const.tolist()}")
    print(f"principal_moments (amu*bohr^2) = {moments.tolist()}")
    print(f"freq_wavenumber (cm^-1, sorted asc) = {freq_cm1.tolist()}")
    print(f"freq_au (sorted) = {freq_au.tolist()}")
    print(f"reduced_mass (amu) = {reduced_mass.tolist()}")
    print(f"force_const_dyne (mDyne/A) = {force_const_dyne.tolist()}")
    print(f"force_const_au (Ha/bohr^2) = {force_const_au.tolist()}")

    # Collect atoms as (Z, [x,y,z] in bohr)
    atoms_out = []
    for i in range(natm):
        z = int(mol.atom_charge(i))
        pos = atom_coords[i].tolist()
        atoms_out.append({"Z": z, "symbol": mol.atom_symbol(i), "pos_bohr": pos})

    # Replace Infinity with a very large finite value for JSON compatibility.
    # IQCP's classify_rotor only checks `> 1e8`, so any large finite value works.
    def _sanitize(arr):
        return [1e300 if (not np.isfinite(x)) and x > 0 else (-1e300 if (not np.isfinite(x)) else float(x)) for x in arr]

    rot_const_json = _sanitize(rot_const.tolist())

    data = {
        "_provenance": {
            "generator": "scripts/phase5/generate_thermo_golden.py",
            "pyscf_version": pyscf.__version__,
            "numpy_version": np.__version__,
            "basis": basis,
            "method": "RHF",
            "conv_tol": 1e-12,
            "cart": True,
            "unit_input": "angstrom",
            "note_infinity": "rot_const entries of ±1e300 represent ±Infinity (linear molecules / single atoms)",
        },
        "name": name,
        "atoms": atoms_out,
        "hessian": h_flat.tolist(),
        "energy": float(e_tot),
        "rot_type": rotor_type,
        "rotational_constants_ghz": rot_const_json,
        "principal_moments_amu_bohr2": moments.tolist(),
        "freq_wavenumber": freq_cm1.tolist(),
        "freq_au": freq_au.tolist(),
        "reduced_mass": reduced_mass.tolist(),
        "force_const_dyne": force_const_dyne.tolist(),
        "force_const_au": force_const_au.tolist(),
        "norm_mode": norm_mode.tolist(),
    }

    out_path = OUT_DIR / f"{name}_sto3g_rhf.json"
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Wrote {out_path}")


def main():
    print(f"PySCF version: {pyscf.__version__}")
    print(f"NumPy version: {np.__version__}")
    print(f"Output dir:    {OUT_DIR}")

    # --- H2O (asymmetric top, 3N-6 = 3 modes) ---
    run_molecule(
        "h2o",
        """
        O  0.0000000000  0.0000000000  0.1173000000
        H  0.0000000000  0.7572000000 -0.4692000000
        H  0.0000000000 -0.7572000000 -0.4692000000
        """,
    )

    # --- CH4 (spherical top, 3N-6 = 9 modes) ---
    run_molecule(
        "ch4",
        """
        C  0.0000000000  0.0000000000  0.0000000000
        H  0.6287000000  0.6287000000  0.6287000000
        H -0.6287000000 -0.6287000000  0.6287000000
        H -0.6287000000  0.6287000000 -0.6287000000
        H  0.6287000000 -0.6287000000 -0.6287000000
        """,
    )

    # --- CO2 (linear, 3N-5 = 4 modes) ---
    run_molecule(
        "co2",
        """
        C  0.0000000000  0.0000000000  0.0000000000
        O  0.0000000000  0.0000000000  1.1600000000
        O  0.0000000000  0.0000000000 -1.1600000000
        """,
    )

    # --- H2 (linear, 3N-5 = 1 mode) ---
    run_molecule(
        "h2",
        """
        H  0.0000000000  0.0000000000  0.0000000000
        H  0.0000000000  0.0000000000  0.7400000000
        """,
    )

    print()
    print("Done. Wrote golden data to", OUT_DIR)


if __name__ == "__main__":
    main()
