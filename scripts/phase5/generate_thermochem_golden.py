#!/usr/bin/env python3
"""Generate PySCF RRHO thermochemistry golden reference data for US-099.

Runs PySCF's RHF/Hessian/thermo pipeline on H2O, CH4, CO2, H2 at
STO-3G with ``mol.cart=True``, then dumps for each molecule:

  - atoms (Z, symbol, pos_bohr)               ← geometry in bohr
  - scf_energy_ha                              ← electronic energy
  - rot_type                                   ← "ATOM" / "LINEAR" / "REGULAR"
  - rotational_constants_ghz                   ← [A, B, C] (GHz)
  - freq_wavenumber_cm1                        ← vibrational frequencies (cm^-1)
  - freq_au                                    ← same in atomic units
  - symmetry_number                            ← sigma (PySCF-detected)
  - multiplicity                               ← 2S+1
  - temperature_k, pressure_pa                 ← reference conditions
  - thermo                                     ← full pyscf.hessian.thermo.thermo() dict

to JSON files under ``tests/golden/thermochem/{h2o,ch4,co2,h2}_sto3g_rhf.json``.

The JSONs are the single source of truth for IQCP's RRHO integration tests in
``crates/qc-core/src/thermochemistry.rs``. IQCP consumes these via
``include_str!(...)`` + ``serde_json::from_str(...)`` and runs
``compute_thermochemistry()`` on its own Hessian-derived ``FrequencyInfo``.

Usage (from ````):

    source .venv/bin/activate
    python scripts/phase5/generate_thermochem_golden.py

This script is independent of any IQCP Rust code — the reference is
generated purely from PySCF 2.11.0.

PySCF line-level references:

  - ``pyscf/hessian/thermo.py::thermo()``  (lines 135-230)
  - ``pyscf/hessian/thermo.py::harmonic_analysis()``
  - ``pyscf/hessian/thermo.py::rotational_symmetry_number()``

See also: ``docs/stories/US-099_rrho_thermochemistry.md`` Section 8.4.
"""

import json
from pathlib import Path

import numpy as np
import pyscf
from pyscf import gto, scf
from pyscf import hessian as hess_mod
from pyscf.hessian import thermo as pyscf_thermo


OUT_DIR = Path(__file__).resolve().parents[2] / "tests" / "golden" / "thermochem"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TEMPERATURE_K = 298.15
PRESSURE_PA = 101325.0


def _sanitize(arr):
    """Replace +/- infinity with +/- 1e300 so the JSON is a valid number.

    IQCP's rotor classification and linear-rotor branch both accept 1e300 as
    "infinity-like"; the existing thermo golden JSONs use the same convention.
    """
    out = []
    for x in arr:
        if np.isfinite(x):
            out.append(float(x))
        elif x > 0:
            out.append(1e300)
        else:
            out.append(-1e300)
    return out


def run_molecule(name, atom_str, basis="sto-3g"):
    print()
    print("=" * 78)
    print(f"=== {name} / {basis} / RHF — RRHO thermochemistry ===")
    print("=" * 78)

    mol = gto.Mole()
    mol.atom = atom_str
    mol.basis = basis
    mol.cart = True
    mol.unit = "angstrom"
    mol.build()

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.verbose = 0
    e_tot = mf.kernel()
    print(f"SCF energy: {e_tot:.12f} Ha")
    assert mf.converged, f"SCF did not converge for {name}"

    # Analytical Hessian: shape (natm, natm, 3, 3)
    hobj = hess_mod.RHF(mf)
    hess = hobj.kernel()
    natm = mol.natm

    # Harmonic analysis → freq_wavenumber, freq_au, rot_const, rot_type
    freq_info = pyscf_thermo.harmonic_analysis(mol, hess)
    freq_cm1 = np.asarray(freq_info["freq_wavenumber"]).real
    freq_au = np.asarray(freq_info["freq_au"]).real

    # Sort ascending (match US-096 convention)
    sort_idx = np.argsort(freq_cm1)
    freq_cm1 = freq_cm1[sort_idx]
    freq_au = freq_au[sort_idx]

    # Rotational constants (GHz) and rotor type from PySCF internals
    atom_mass = mol.atom_mass_list(isotope_avg=True)
    atom_coords = mol.atom_coords()  # bohr
    rot_const = pyscf_thermo.rotation_const(atom_mass, atom_coords, "GHz")
    rotor_type = pyscf_thermo._get_rotor_type(rot_const)
    sigma = int(pyscf_thermo.rotational_symmetry_number(mol))

    # Run the full PySCF thermo() pipeline
    therm = pyscf_thermo.thermo(mf, freq_info["freq_au"], TEMPERATURE_K, PRESSURE_PA)

    # ========================================================================
    # Pretty-print all key quantities
    # ========================================================================
    print(f"n_atoms         = {natm}")
    print(f"n_modes         = {len(freq_cm1)}")
    print(f"rot_type        = {rotor_type}")
    print(f"rotational_constants_ghz = {rot_const.tolist()}")
    print(f"freq_wavenumber = {freq_cm1.tolist()}")
    print(f"freq_au         = {freq_au.tolist()}")
    print(f"symmetry_number = {sigma}")
    print(f"multiplicity    = {mol.multiplicity}")
    print()
    print(f"{'Quantity':20s}  {'Value':>22s}  Unit")
    print("-" * 60)
    for key in ("ZPE", "E_0K", "E_tot", "H_tot", "G_tot",
                "S_tot", "Cv_tot", "Cp_tot",
                "S_trans", "S_rot", "S_vib", "S_elec",
                "E_trans", "E_rot", "E_vib", "E_elec",
                "Cv_trans", "Cv_rot", "Cv_vib",
                "H_trans", "H_rot", "H_vib", "H_elec",
                "G_trans", "G_rot", "G_vib", "G_elec"):
        val, unit = therm[key]
        print(f"{key:20s}  {float(val):>22.14e}  {unit}")
    print()

    # ========================================================================
    # Build the JSON payload
    # ========================================================================
    atoms_out = []
    for i in range(natm):
        atoms_out.append({
            "Z": int(mol.atom_charge(i)),
            "symbol": mol.atom_symbol(i),
            "pos_bohr": atom_coords[i].tolist(),
        })

    # Convert the (value, unit) tuples in `therm` to a plain dict of scalars
    # Skipping the rot_const entry (already stored separately) and temperature/pressure
    # (echo-only).
    thermo_dict = {}
    for key, val in therm.items():
        raw = val[0]
        if isinstance(raw, np.ndarray) or isinstance(raw, (list, tuple)):
            continue  # skip array-valued entries (e.g. rot_const)
        if not np.isfinite(raw):
            continue  # shouldn't happen, but be defensive
        thermo_dict[key] = {"value": float(raw), "unit": val[1]}

    data = {
        "_provenance": {
            "generator": "scripts/phase5/generate_thermochem_golden.py",
            "pyscf_version": pyscf.__version__,
            "numpy_version": np.__version__,
            "basis": basis,
            "method": "RHF",
            "conv_tol": 1e-12,
            "cart": True,
            "unit_input": "angstrom",
            "temperature_k": TEMPERATURE_K,
            "pressure_pa": PRESSURE_PA,
            "note_infinity": "rot_const entries of +/-1e300 represent +/-Infinity (linear molecules / single atoms)",
        },
        "name": name,
        "atoms": atoms_out,
        "scf_energy_ha": float(e_tot),
        "rot_type": rotor_type,
        "rotational_constants_ghz": _sanitize(rot_const.tolist()),
        "freq_wavenumber_cm1": freq_cm1.tolist(),
        "freq_au": freq_au.tolist(),
        "symmetry_number": sigma,
        "multiplicity": int(mol.multiplicity),
        "temperature_k": TEMPERATURE_K,
        "pressure_pa": PRESSURE_PA,
        "thermo": thermo_dict,
    }

    out_path = OUT_DIR / f"{name}_sto3g_rhf.json"
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Wrote {out_path}")
    print(f"  ZPE = {thermo_dict['ZPE']['value']:.10e} Ha")
    print(f"  G   = {thermo_dict['G_tot']['value']:.10e} Ha")


def main():
    print(f"PySCF version: {pyscf.__version__}")
    print(f"NumPy version: {np.__version__}")
    print(f"Output dir:    {OUT_DIR}")
    print(f"Temperature:   {TEMPERATURE_K} K")
    print(f"Pressure:      {PRESSURE_PA} Pa")

    # --- H2O (asymmetric top, 3N-6 = 3 modes, C2v, sigma=2) ---
    run_molecule(
        "h2o",
        """
        O  0.0000000000  0.0000000000  0.1173000000
        H  0.0000000000  0.7572000000 -0.4692000000
        H  0.0000000000 -0.7572000000 -0.4692000000
        """,
    )

    # --- CH4 (spherical top, 3N-6 = 9 modes, T_d, sigma=12) ---
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

    # --- CO2 (linear, 3N-5 = 4 modes, Dinfh, sigma=2) ---
    run_molecule(
        "co2",
        """
        C  0.0000000000  0.0000000000  0.0000000000
        O  0.0000000000  0.0000000000  1.1600000000
        O  0.0000000000  0.0000000000 -1.1600000000
        """,
    )

    # --- H2 (linear, 3N-5 = 1 mode, Dinfh, sigma=2) ---
    run_molecule(
        "h2",
        """
        H  0.0000000000  0.0000000000  0.0000000000
        H  0.0000000000  0.0000000000  0.7400000000
        """,
    )

    print()
    print("Done. 4 thermochemistry golden JSON files written.")


if __name__ == "__main__":
    main()
