#!/usr/bin/env python3
"""Generate comprehensive golden reference data for US-104 Phase 5 golden tests.

Generates PySCF reference data for molecules NOT already covered by existing
golden files (US-096 through US-100):
  - HF / STO-3G / RHF  (new)
  - NH3 / STO-3G / RHF  (new, with geometry optimization)
  - H2O / 6-31G* / RHF  (new)
  - H2O / STO-3G / LDA   (new)
  - H2O / STO-3G / B3LYP (new)

Usage:
    source .venv/bin/activate
    python scripts/phase5/generate_phase5_golden.py

Output: tests/golden/phase5/*.json
"""

import json
import os
import datetime
import numpy as np

import pyscf
from pyscf import gto, scf, dft
from pyscf.hessian import rhf as rhf_hess
from pyscf.hessian import rks as rks_hess
from pyscf.hessian.thermo import harmonic_analysis

AU_TO_DEBYE = 2.541746473

# H2O geometry (same as US-096, bohr)
H2O_BOHR = [
    ("O", [0.0, 0.0, 0.22166487441148175]),
    ("H", [0.0, 1.4309006215206648, -0.886659497645927]),
    ("H", [0.0, -1.4309006215206648, -0.886659497645927]),
]

# HF geometry (0.917 A = 1.7328 bohr)
HF_BOHR = [
    ("F", [0.0, 0.0, 0.0]),
    ("H", [0.0, 0.0, 1.7328]),
]

# NH3 initial geometry (will be optimized)
NH3_ANGSTROM = [
    ("N", [0.0, 0.0, 0.1169]),
    ("H", [0.0, 0.9377, -0.2727]),
    ("H", [0.8121, -0.4689, -0.2727]),
    ("H", [-0.8121, -0.4689, -0.2727]),
]

OUTPUT_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "../../tests/golden/phase5",
)


def make_mol(atoms, basis="sto-3g", unit="bohr", cart=True):
    mol = gto.Mole()
    atom_str = "; ".join(f"{sym} {x} {y} {z}" for sym, (x, y, z) in atoms)
    mol.atom = atom_str
    mol.basis = basis
    mol.unit = unit
    mol.cart = cart
    mol.build()
    return mol


def run_scf(mol, method="rhf"):
    if method.lower() == "rhf":
        mf = scf.RHF(mol)
    elif method.lower() in ("lda", "svwn"):
        mf = dft.RKS(mol)
        mf.xc = "svwn"
    elif method.lower() in ("b3lyp",):
        mf = dft.RKS(mol)
        mf.xc = "b3lyp5"
    else:
        raise ValueError(f"Unknown method: {method}")
    mf.conv_tol = 1e-12
    mf.kernel()
    assert mf.converged, f"SCF did not converge"
    return mf


def compute_hessian(mf):
    if isinstance(mf, dft.rks.RKS):
        hess_obj = rks_hess.Hessian(mf)
    else:
        hess_obj = rhf_hess.Hessian(mf)
    return hess_obj.kernel()


def hessian_4d_to_2d(hess_4d, natm):
    """Convert PySCF 4D Hessian (natm, natm, 3, 3) to 2D (3N, 3N).

    PySCF convention: hess[A, B, i, j] = d^2E / (dR_Ai dR_Bj)
    We want H[3A+i, 3B+j].
    """
    n3 = 3 * natm
    return hess_4d.transpose(0, 2, 1, 3).reshape(n3, n3)


def compute_dipole_derivs_fd(mf, dx=1e-4):
    """Finite difference of dipole moment w.r.t. nuclear coordinates."""
    mol = mf.mol
    natm = mol.natm
    n3 = 3 * natm
    coords_orig = mol.atom_coords().copy()

    derivs = np.zeros((3, n3))
    for a in range(natm):
        for g in range(3):
            coords_p = coords_orig.copy()
            coords_p[a, g] += dx
            mol.set_geom_(coords_p, unit="Bohr")
            mol.build()
            mf_p = mf.__class__(mol)
            if hasattr(mf, 'xc'):
                mf_p.xc = mf.xc
            mf_p.conv_tol = 1e-12
            mf_p.kernel()
            mu_p = np.array(mf_p.dip_moment(unit="AU"))

            coords_m = coords_orig.copy()
            coords_m[a, g] -= dx
            mol.set_geom_(coords_m, unit="Bohr")
            mol.build()
            mf_m = mf.__class__(mol)
            if hasattr(mf, 'xc'):
                mf_m.xc = mf.xc
            mf_m.conv_tol = 1e-12
            mf_m.kernel()
            mu_m = np.array(mf_m.dip_moment(unit="AU"))

            derivs[:, 3*a + g] = (mu_p - mu_m) / (2 * dx)

    mol.set_geom_(coords_orig, unit="Bohr")
    mol.build()
    return derivs


def compute_polarizability_fd(mf, dx=1e-4):
    """Compute static polarizability and its derivatives via FD."""
    from pyscf.prop.polarizability import rhf as rhf_polar
    from pyscf.prop.polarizability import rks as rks_polar

    mol = mf.mol
    natm = mol.natm
    n3 = 3 * natm

    if isinstance(mf, dft.rks.RKS):
        polar = rks_polar.Polarizability(mf)
    else:
        polar = rhf_polar.Polarizability(mf)
    alpha_eq = polar.polarizability()

    dalpha = np.zeros((3, 3, n3))
    coords_orig = mol.atom_coords().copy()

    for a in range(natm):
        for g in range(3):
            coords_p = coords_orig.copy()
            coords_p[a, g] += dx
            mol.set_geom_(coords_p, unit="Bohr")
            mol.build()
            mf_p = mf.__class__(mol)
            if hasattr(mf, 'xc'):
                mf_p.xc = mf.xc
            mf_p.conv_tol = 1e-12
            mf_p.kernel()
            if isinstance(mf_p, dft.rks.RKS):
                polar_p = rks_polar.Polarizability(mf_p)
            else:
                polar_p = rhf_polar.Polarizability(mf_p)
            alpha_p = polar_p.polarizability()

            coords_m = coords_orig.copy()
            coords_m[a, g] -= dx
            mol.set_geom_(coords_m, unit="Bohr")
            mol.build()
            mf_m = mf.__class__(mol)
            if hasattr(mf, 'xc'):
                mf_m.xc = mf.xc
            mf_m.conv_tol = 1e-12
            mf_m.kernel()
            if isinstance(mf_m, dft.rks.RKS):
                polar_m = rks_polar.Polarizability(mf_m)
            else:
                polar_m = rhf_polar.Polarizability(mf_m)
            alpha_m = polar_m.polarizability()

            dalpha[:, :, 3*a + g] = (alpha_p - alpha_m) / (2 * dx)

    mol.set_geom_(coords_orig, unit="Bohr")
    mol.build()
    return alpha_eq, dalpha


def compute_ir(dipole_derivs, norm_modes):
    """IR intensities from dipole derivs and normal modes."""
    n_modes = len(norm_modes)
    AU_TO_KMMOL = 974.8801193081486
    ir = []
    dmu_dQ_all = []
    for k in range(n_modes):
        dmu_dQ = np.zeros(3)
        for d in range(3):
            dmu_dQ[d] = np.dot(dipole_derivs[d], norm_modes[k])
        dsq = np.sum(dmu_dQ**2)
        ir.append(float(dsq * AU_TO_KMMOL))
        dmu_dQ_all.append(dmu_dQ.tolist())
    return ir, dmu_dQ_all


def compute_raman(dalpha, norm_modes):
    """Raman activities from polarizability derivs and normal modes."""
    n_modes = len(norm_modes)
    n3 = len(norm_modes[0])
    raman = []
    depol = []
    for k in range(n_modes):
        dalpha_dQ = np.zeros((3, 3))
        for i in range(n3):
            dalpha_dQ += dalpha[:, :, i] * norm_modes[k][i]

        bar_alpha = np.trace(dalpha_dQ) / 3.0
        gamma2 = 0.0
        for i in range(3):
            for j in range(3):
                if i == j:
                    gamma2 += 0.5 * (dalpha_dQ[i, i] - dalpha_dQ[(i+1) % 3, (i+1) % 3])**2
                else:
                    gamma2 += 0.75 * dalpha_dQ[i, j]**2

        raman.append(float(45.0 * bar_alpha**2 + 7.0 * gamma2))
        denom = 45.0 * bar_alpha**2 + 4.0 * gamma2
        depol.append(float(3.0 * gamma2 / denom) if denom > 1e-20 else 0.75)

    return raman, depol


def compute_rot_info(mol):
    """Rotational constants and rotor type."""
    natm = mol.natm
    coords = mol.atom_coords()
    masses = mol.atom_mass_list(isotope_avg=True)
    com = np.average(coords, weights=masses, axis=0)
    coords_com = coords - com

    I = np.zeros((3, 3))
    for i in range(natm):
        r = coords_com[i]
        m = masses[i]
        I[0, 0] += m * (r[1]**2 + r[2]**2)
        I[1, 1] += m * (r[0]**2 + r[2]**2)
        I[2, 2] += m * (r[0]**2 + r[1]**2)
        I[0, 1] -= m * r[0] * r[1]
        I[0, 2] -= m * r[0] * r[2]
        I[1, 2] -= m * r[1] * r[2]
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]

    # Moments in amu*bohr^2
    moments = np.sort(np.linalg.eigvalsh(I))[::-1]

    AMU_TO_KG = 1.66053906660e-27
    BOHR_M = 5.29177210903e-11
    HBAR = 1.0545718e-34
    PI = np.pi

    rot_ghz = []
    for mom in moments:
        mom_si = mom * AMU_TO_KG * BOHR_M**2
        if mom_si > 1e-60:
            rot_ghz.append(float(HBAR / (4 * PI * mom_si) / 1e9))
        else:
            rot_ghz.append(float("inf"))

    if natm == 1:
        rot_type = "atom"
    elif natm == 2:
        rot_type = "linear"
    else:
        tol = 1e-4
        if moments[2] < tol * moments[0]:
            rot_type = "linear"
        elif (abs(moments[0] - moments[1]) < tol * moments[0] and
              abs(moments[1] - moments[2]) < tol * moments[1]):
            rot_type = "spherical_top"
        elif (abs(moments[0] - moments[1]) < tol * moments[0] or
              abs(moments[1] - moments[2]) < tol * moments[1]):
            rot_type = "symmetric_top"
        else:
            rot_type = "asymmetric_top"

    return rot_ghz, rot_type, moments


def compute_thermochem(mol, mf, freq_cm1, temperature=298.15, pressure=101325.0):
    """RRHO thermochemistry with pre-SI-revision constants (matching IQCP)."""
    kB = 1.38064852e-23
    h = 6.62607004e-34
    NA = 6.022140857e23
    R = kB * NA
    c = 2.99792458e10
    AMU_TO_KG = 1.66053906660e-27
    BOHR_M = 5.29177210903e-11
    HARTREE_J = 4.3597447222071e-18
    PI = np.pi

    T = temperature
    P = pressure
    E_scf = mf.e_tot
    natm = mol.natm

    real_freqs = [f for f in freq_cm1 if f > 0]
    zpe_J = 0.5 * h * c * sum(real_freqs) * NA
    zpe_Ha = zpe_J / (HARTREE_J * NA)

    mult = mol.spin + 1
    s_elec = R * np.log(mult) if mult > 1 else 0.0

    masses = mol.atom_mass_list(isotope_avg=True)
    total_mass_kg = sum(masses) * AMU_TO_KG
    lambda_th = h / np.sqrt(2 * PI * total_mass_kg * kB * T)
    q_trans = (kB * T / P) / lambda_th**3
    s_trans = R * (np.log(q_trans) + 2.5)
    e_trans = 1.5 * R * T
    cv_trans = 1.5 * R

    if natm == 1:
        s_rot = e_rot = cv_rot = 0.0
    else:
        coords = mol.atom_coords()
        com = np.average(coords, weights=masses, axis=0)
        I_tensor = np.zeros((3, 3))
        for i in range(natm):
            r_m = (coords[i] - com) * BOHR_M
            m = masses[i] * AMU_TO_KG
            I_tensor[0, 0] += m * (r_m[1]**2 + r_m[2]**2)
            I_tensor[1, 1] += m * (r_m[0]**2 + r_m[2]**2)
            I_tensor[2, 2] += m * (r_m[0]**2 + r_m[1]**2)
            I_tensor[0, 1] -= m * r_m[0] * r_m[1]
            I_tensor[0, 2] -= m * r_m[0] * r_m[2]
            I_tensor[1, 2] -= m * r_m[1] * r_m[2]
        I_tensor[1, 0] = I_tensor[0, 1]
        I_tensor[2, 0] = I_tensor[0, 2]
        I_tensor[2, 1] = I_tensor[1, 2]
        evals = np.sort(np.linalg.eigvalsh(I_tensor))[::-1]
        is_linear = evals[2] < 1e-40
        sigma = 1

        if is_linear:
            theta = h**2 / (8 * PI**2 * evals[0] * kB)
            q_rot = T / (sigma * theta)
            s_rot = R * (np.log(q_rot) + 1.0)
            e_rot = R * T
            cv_rot = R
        else:
            q_rot = (np.sqrt(PI) / sigma) * (T**1.5) * np.sqrt(
                (8 * PI**2 * kB)**3 * np.prod(evals)
            ) / h**3
            s_rot = R * (np.log(q_rot) + 1.5)
            e_rot = 1.5 * R * T
            cv_rot = 1.5 * R

    e_vib = s_vib = cv_vib = 0.0
    for f in real_freqs:
        x = h * c * f / (kB * T)
        if x < 500:
            e_vib += R * T * x / (np.exp(x) - 1)
            s_vib += R * (x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)))
            cv_vib += R * x**2 * np.exp(x) / (np.exp(x) - 1)**2

    conv = 1.0 / (HARTREE_J * NA)
    e_total = E_scf + (zpe_J + e_trans + e_rot + e_vib) * conv
    pV = R * T * conv
    h_total = e_total + pV
    s_total = (s_elec + s_trans + s_rot + s_vib) * conv
    g_total = h_total - T * s_total
    cv_total = (cv_trans + cv_rot + cv_vib) * conv
    cp_total = cv_total + R * conv

    return {
        "temperature": temperature,
        "pressure": pressure,
        "zpe_hartree": float(zpe_Ha),
        "e_total": float(e_total),
        "h_total": float(h_total),
        "s_total": float(s_total),
        "g_total": float(g_total),
        "cv_total": float(cv_total),
        "cp_total": float(cp_total),
    }


def generate_golden(mol_name, atoms, basis, method, unit="bohr", optimize=False):
    """Generate full golden reference for one molecule/method/basis."""
    print(f"\n{'='*60}")
    print(f"Generating: {mol_name} / {basis} / {method}")
    print(f"{'='*60}")

    mol = make_mol(atoms, basis=basis, unit=unit, cart=True)

    if optimize:
        print("  Optimizing geometry...")
        from pyscf.geomopt import geometric_solver
        mf_opt = run_scf(mol, method)
        mol = geometric_solver.optimize(mf_opt, maxsteps=100)
        print(f"  Optimization converged")
        for i in range(mol.natm):
            c = mol.atom_coords()[i]
            print(f"    {mol.atom_symbol(i)} {c[0]:.10f} {c[1]:.10f} {c[2]:.10f}")

    natm = mol.natm
    print(f"  Atoms: {natm}, AOs: {mol.nao_nr()}")

    # SCF
    mf = run_scf(mol, method)
    energy = float(mf.e_tot)
    print(f"  Energy: {energy:.12f} Ha")

    # Dipole
    mu_au = np.array(mf.dip_moment(unit="AU"))
    mu_debye = mu_au * AU_TO_DEBYE
    print(f"  Dipole: {np.linalg.norm(mu_debye):.6f} D")

    # Hessian
    hess_4d = compute_hessian(mf)
    n3 = 3 * natm
    hess_2d = hessian_4d_to_2d(hess_4d, natm)
    print(f"  Hessian: {n3}x{n3}")

    # PySCF harmonic analysis
    ha = harmonic_analysis(mol, hess_4d)
    freq_cm1 = ha["freq_wavenumber"].tolist()
    freq_au = ha["freq_au"].tolist()
    reduced_mass = ha["reduced_mass"].tolist()
    force_const = ha["force_const_dyne"].tolist()
    n_modes = len(freq_cm1)

    # Normal modes: ha['norm_mode'] is (n_modes, natm, 3)
    # Flatten to (n_modes, 3N) for our format
    norm_mode_raw = ha["norm_mode"]
    norm_modes_flat = []
    for k in range(n_modes):
        mode = norm_mode_raw[k].flatten().tolist()
        norm_modes_flat.append(mode)

    print(f"  Frequencies: {[f'{f:.1f}' for f in freq_cm1]}")

    # Rotational info
    rot_ghz, rot_type, _ = compute_rot_info(mol)
    print(f"  Rotor: {rot_type}")

    # Dipole derivatives + IR
    print("  Computing dipole derivs (FD)...")
    dipole_derivs = compute_dipole_derivs_fd(mf)
    ir_intensities, dmu_dQ = compute_ir(dipole_derivs, norm_modes_flat)
    print(f"  IR: {[f'{i:.2f}' for i in ir_intensities]}")

    # Polarizability + Raman
    print("  Computing polar derivs (FD)...")
    alpha_eq, dalpha = compute_polarizability_fd(mf)
    raman_activities, depol_ratios = compute_raman(dalpha, norm_modes_flat)
    print(f"  Raman: {[f'{r:.2f}' for r in raman_activities]}")

    # Thermochemistry
    thermo = compute_thermochem(mol, mf, freq_cm1)
    print(f"  ZPE: {thermo['zpe_hartree']:.8f} Ha")
    print(f"  G: {thermo['g_total']:.8f} Ha")

    # Assemble
    atoms_list = []
    coords_bohr = mol.atom_coords()
    for i in range(natm):
        atoms_list.append({
            "Z": int(mol.atom_charge(i)),
            "symbol": mol.atom_symbol(i),
            "pos_bohr": coords_bohr[i].tolist(),
        })

    golden = {
        "_provenance": {
            "generator": "scripts/phase5/generate_phase5_golden.py",
            "pyscf_version": pyscf.__version__,
            "numpy_version": np.__version__,
            "date": datetime.datetime.now().isoformat(),
            "basis": basis,
            "method": method.upper(),
            "conv_tol": 1e-12,
            "cart": True,
            "unit": "bohr",
            "description": f"US-104 Phase 5 golden reference for {mol_name}",
        },
        "name": mol_name,
        "atoms": atoms_list,
        "n_atoms": natm,
        "n_modes": n_modes,
        "energy": energy,
        "hessian": hess_2d.flatten().tolist(),
        "freq_wavenumber": freq_cm1,
        "freq_au": freq_au,
        "reduced_mass": reduced_mass,
        "force_const_dyne": force_const,
        "dipole_au": mu_au.tolist(),
        "dipole_debye": mu_debye.tolist(),
        "dipole_derivs_cartesian": dipole_derivs.tolist(),
        "dmu_dQ": dmu_dQ,
        "ir_intensities_km_mol": ir_intensities,
        "polarizability_au": alpha_eq.tolist(),
        "raman_activities_a4_amu": raman_activities,
        "depolarization_ratios": depol_ratios,
        "rot_type": rot_type,
        "rotational_constants_ghz": rot_ghz,
        "norm_mode": norm_modes_flat,
        "thermochemistry": thermo,
    }

    return golden


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output: {OUTPUT_DIR}")
    print(f"PySCF {pyscf.__version__}, NumPy {np.__version__}")

    # HF / STO-3G / RHF
    g = generate_golden("hf", HF_BOHR, "sto-3g", "rhf", unit="bohr")
    with open(os.path.join(OUTPUT_DIR, "hf_sto3g_rhf.json"), "w") as f:
        json.dump(g, f, indent=2)

    # NH3 / STO-3G / RHF (optimize)
    g = generate_golden("nh3", NH3_ANGSTROM, "sto-3g", "rhf",
                        unit="angstrom", optimize=True)
    with open(os.path.join(OUTPUT_DIR, "nh3_sto3g_rhf.json"), "w") as f:
        json.dump(g, f, indent=2)

    # H2O / 6-31G* / RHF
    g = generate_golden("h2o", H2O_BOHR, "6-31g*", "rhf", unit="bohr")
    with open(os.path.join(OUTPUT_DIR, "h2o_631gs_rhf.json"), "w") as f:
        json.dump(g, f, indent=2)

    # H2O / STO-3G / LDA
    g = generate_golden("h2o", H2O_BOHR, "sto-3g", "lda", unit="bohr")
    with open(os.path.join(OUTPUT_DIR, "h2o_sto3g_lda.json"), "w") as f:
        json.dump(g, f, indent=2)

    # H2O / STO-3G / B3LYP
    g = generate_golden("h2o", H2O_BOHR, "sto-3g", "b3lyp", unit="bohr")
    with open(os.path.join(OUTPUT_DIR, "h2o_sto3g_b3lyp.json"), "w") as f:
        json.dump(g, f, indent=2)

    print(f"\nDone! All golden files in {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
