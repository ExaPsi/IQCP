#!/usr/bin/env python3
"""Generate PySCF golden reference data for US-098 (Field CPHF + Raman Activity).

Runs PySCF's RHF/Hessian pipeline on H2, H2O, CH4, CO2 at STO-3G with
mol.cart=True, then computes:

  - static polarizability via pyscf-forge Polarizability.polarizability()
  - polarizability derivatives via finite difference of the analytical
    polarizability at displaced geometries (central FD, step 1e-3 bohr).
    The CPHF is analytical at each sample point — only the outer
    differentiation is numerical. This matches exactly what IQCP's Rust
    implementation does in compute_polarizability_derivatives.
  - Raman invariants + scattering activities + depolarization ratios
    using the Long (2002) Chapter 5 formulas.

The output JSON files in tests/golden/raman/ are consumed by the Rust
golden tests in crates/qc-core/src/raman.rs.

Note on "analytical" polarizability derivatives:
  pyscf-forge does not ship an analytical polarizability-derivative
  module. Deriving the correct interchange-theorem formula from scratch
  is non-trivial (see US-098 derive notes). Instead, both IQCP and this
  script compute ∂α/∂R by central finite differences of the analytical
  polarizability at displaced geometries — the CPHF is fully analytical
  at each sample point, giving a correct and auditable result. Cost is
  6N SCF+CPHF solves per molecule.

Usage (from ):
    source .venv/bin/activate
    python scripts/phase5/generate_raman_golden.py
"""

import json
from pathlib import Path

import numpy as np
import pyscf
from pyscf import gto, scf
from pyscf import hessian as hess_mod
from pyscf.hessian import thermo as pyscf_thermo

# pyscf-forge Polarizability module (from pyscf-properties pip package).
from pyscf.prop.polarizability.rhf import Polarizability

OUT_DIR = Path(__file__).resolve().parents[2] / "tests" / "golden" / "raman"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BOHR_TO_ANG = 0.52917721092
BOHR4_TO_ANG4 = BOHR_TO_ANG ** 4  # ≈ 0.07842
AU_POLARIZABILITY_TO_ANG3 = BOHR_TO_ANG ** 3  # ≈ 0.14818


def compute_alpha_at_geometry(mol_template, coords_bohr, basis="sto-3g"):
    """Compute the analytical polarizability at the given bohr coordinates.

    Uses pyscf-forge Polarizability for analytical CPHF.
    """
    mol = gto.Mole()
    mol.atom = [
        (mol_template.atom_symbol(i), tuple(coords_bohr[i]))
        for i in range(mol_template.natm)
    ]
    mol.basis = basis
    mol.cart = True
    mol.unit = "bohr"
    mol.build()
    mf = scf.RHF(mol).run(conv_tol=1e-13, verbose=0)
    polobj = Polarizability(mf)
    polobj.conv_tol = 1e-11
    return polobj.polarizability()  # 3×3 in bohr³


def compute_polar_derivs_fd(mol, step=1e-3, basis="sto-3g"):
    """Compute ∂α/∂R via central finite difference of the analytical
    polarizability at displaced geometries.

    Returns
    -------
    polar_deriv : (3, 3, 3*natm) array
        polar_deriv[d, e, 3*A + γ] = ∂α_{de}/∂R_{A,γ} in au/bohr
    """
    natm = mol.natm
    coords0 = mol.atom_coords().copy()  # bohr
    polar_deriv = np.zeros((3, 3, 3 * natm))
    for A in range(natm):
        for gamma in range(3):
            cp = coords0.copy()
            cp[A, gamma] += step
            alpha_p = compute_alpha_at_geometry(mol, cp, basis)
            cm = coords0.copy()
            cm[A, gamma] -= step
            alpha_m = compute_alpha_at_geometry(mol, cm, basis)
            polar_deriv[:, :, 3 * A + gamma] = (alpha_p - alpha_m) / (2.0 * step)
    # Symmetrize d<->e
    polar_deriv = 0.5 * (polar_deriv + polar_deriv.transpose(1, 0, 2))
    return polar_deriv


def compute_raman_from_polar_deriv(polar_deriv, norm_mode):
    """Compute Raman invariants and activities from the Cartesian polar
    derivative and normal modes.

    Parameters
    ----------
    polar_deriv : (3, 3, 3*natm) array, au/bohr
    norm_mode : (n_modes, natm, 3) array, Cartesian displacements

    Returns
    -------
    dict with activities_a4_amu, depolarization_ratios,
    isotropic_invariants_au, anisotropy_invariants_au,
    polar_derivs_per_mode
    """
    n_modes = norm_mode.shape[0]
    natm = norm_mode.shape[1]
    activities = np.zeros(n_modes)
    depols = np.zeros(n_modes)
    bar_sq_list = np.zeros(n_modes)
    g_sq_list = np.zeros(n_modes)
    per_mode = np.zeros((n_modes, 3, 3))
    for k in range(n_modes):
        alpha_prime = np.zeros((3, 3))
        for d in range(3):
            for e in range(3):
                s = 0.0
                for A in range(natm):
                    for gamma in range(3):
                        s += polar_deriv[d, e, 3 * A + gamma] * norm_mode[k, A, gamma]
                alpha_prime[d, e] = s
        per_mode[k] = alpha_prime
        xx, yy, zz = alpha_prime[0, 0], alpha_prime[1, 1], alpha_prime[2, 2]
        xy, yz, xz = alpha_prime[0, 1], alpha_prime[1, 2], alpha_prime[0, 2]
        bar = (xx + yy + zz) / 3.0
        bar_sq = bar * bar
        g_sq = (
            0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2)
            + 3.0 * (xy ** 2 + yz ** 2 + xz ** 2)
        )
        act = (45.0 * bar_sq + 7.0 * g_sq) * BOHR4_TO_ANG4
        denom = 45.0 * bar_sq + 4.0 * g_sq
        rho = 3.0 * g_sq / denom if abs(denom) > 1e-30 else 0.0
        rho = max(0.0, min(rho, 0.75))
        activities[k] = act
        depols[k] = rho
        bar_sq_list[k] = bar_sq
        g_sq_list[k] = g_sq
    return dict(
        activities_a4_amu=activities,
        depolarization_ratios=depols,
        isotropic_invariants_au=bar_sq_list,
        anisotropy_invariants_au=g_sq_list,
        polar_derivs_per_mode=per_mode,
    )


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
    assert mf.converged, f"SCF did not converge for {name}"
    print(f"SCF energy: {e_tot:.12f} Ha")
    natm = mol.natm

    # Polarizability via pyscf-forge
    polobj = Polarizability(mf)
    polobj.conv_tol = 1e-11
    alpha = polobj.polarizability()
    print(f"polarizability (trace): {np.trace(alpha):.6f} bohr³")
    print(f"  diag: {np.diag(alpha)}")
    alpha_ang3 = alpha * AU_POLARIZABILITY_TO_ANG3

    # Nuclear Hessian for normal modes
    h_obj = hess_mod.RHF(mf)
    h_obj.verbose = 0
    hessian = h_obj.kernel()  # (natm, natm, 3, 3)

    freq_info = pyscf_thermo.harmonic_analysis(mol, hessian)
    freq_cm1 = np.asarray(freq_info["freq_wavenumber"]).real
    norm_mode = np.asarray(freq_info["norm_mode"]).real
    n_modes = len(freq_cm1)

    # Sort by frequency
    sort_idx = np.argsort(freq_cm1)
    freq_cm1 = freq_cm1[sort_idx]
    norm_mode = norm_mode[sort_idx]

    # Polarizability derivatives (central FD of analytical alpha)
    polar_deriv = compute_polar_derivs_fd(mol, step=1e-3, basis=basis)
    print(f"polar_deriv max: {np.abs(polar_deriv).max():.6f}")

    # Raman invariants and activities
    raman_data = compute_raman_from_polar_deriv(polar_deriv, norm_mode)
    print(f"frequencies (cm⁻¹, sorted): {freq_cm1.tolist()}")
    print(f"Raman activities (Å⁴/amu): {raman_data['activities_a4_amu'].tolist()}")
    print(f"Depolarization ratios: {raman_data['depolarization_ratios'].tolist()}")

    # Atom metadata
    atom_coords = mol.atom_coords()
    atomic_numbers = [int(mol.atom_charge(i)) for i in range(natm)]
    atom_symbols = [mol.atom_symbol(i) for i in range(natm)]
    atoms_out = []
    for i in range(natm):
        atoms_out.append(
            {
                "Z": atomic_numbers[i],
                "symbol": atom_symbols[i],
                "pos_bohr": atom_coords[i].tolist(),
            }
        )

    data = {
        "_provenance": {
            "generator": "scripts/phase5/generate_raman_golden.py",
            "pyscf_version": pyscf.__version__,
            "numpy_version": np.__version__,
            "basis": basis,
            "method": "RHF",
            "conv_tol_scf": 1e-12,
            "conv_tol_cphf": 1e-11,
            "cart": True,
            "unit_input": "angstrom",
            "gauge_origin": "charge_center",
            "bohr4_to_ang4_factor": BOHR4_TO_ANG4,
            "au_polarizability_to_ang3_factor": AU_POLARIZABILITY_TO_ANG3,
            "fd_step_bohr": 1e-3,
            "description": "US-098 field CPHF + polarizability derivatives + Raman activities golden reference",
        },
        "name": name,
        "atoms": atoms_out,
        "n_atoms": natm,
        "n_modes": int(n_modes),
        "energy": float(e_tot),
        "polarizability_au": alpha.tolist(),
        "polarizability_ang3": alpha_ang3.tolist(),
        "polarizability_derivs_cartesian": polar_deriv.tolist(),
        "freq_wavenumber": freq_cm1.tolist(),
        "norm_mode": norm_mode.tolist(),
        "raman_activities_a4_amu": raman_data["activities_a4_amu"].tolist(),
        "depolarization_ratios": raman_data["depolarization_ratios"].tolist(),
        "isotropic_invariants_au": raman_data["isotropic_invariants_au"].tolist(),
        "anisotropy_invariants_au": raman_data["anisotropy_invariants_au"].tolist(),
        "polarizability_derivs_normal_mode": raman_data["polar_derivs_per_mode"].tolist(),
    }

    out_path = OUT_DIR / f"{name}_sto3g_rhf.json"
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Wrote {out_path}")
    return data


def main():
    print(f"PySCF version: {pyscf.__version__}")
    print(f"NumPy version: {np.__version__}")
    print(f"Output dir:    {OUT_DIR}")
    print(f"BOHR4_TO_ANG4: {BOHR4_TO_ANG4:.6e}")

    # --- H2O (asymmetric top, 3 modes) ---
    run_molecule(
        "h2o",
        """
        O  0.0000000000  0.0000000000  0.1173000000
        H  0.0000000000  0.7572000000 -0.4692000000
        H  0.0000000000 -0.7572000000 -0.4692000000
        """,
    )

    # --- CH4 (tetrahedral Td, 9 modes with A1 symmetric stretch) ---
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

    # --- CO2 (linear D∞h, 4 modes with mutual exclusion) ---
    run_molecule(
        "co2",
        """
        C  0.0000000000  0.0000000000  0.0000000000
        O  0.0000000000  0.0000000000  1.1600000000
        O  0.0000000000  0.0000000000 -1.1600000000
        """,
    )

    # --- H2 (homonuclear diatomic, 1 mode, σ_g is Raman-active) ---
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
