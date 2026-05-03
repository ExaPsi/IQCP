#!/usr/bin/env python3
"""
PySCF population analysis (Mulliken and true Löwdin).

Regenerates `tests/golden/dft/population_validation.json` with correct
Löwdin **charges** (Z - Löwdin population), not raw populations.

The previous version of the file labeled raw Löwdin populations
(e.g., 0.637 on each H of H2) as "lowdin_charges", which is incorrect:
charges on a neutral molecule must sum to zero.
"""
import json
from pathlib import Path

import numpy as np
import scipy.linalg
from pyscf import gto, scf, dft

SYSTEMS = [
    ("H2", "sto-3g", "RHF", "H 0.0 0.0 0.0; H 0.0 0.0 0.74", "Angstrom"),
    ("H2", "sto-3g", "B3LYP", "H 0.0 0.0 0.0; H 0.0 0.0 0.74", "Angstrom"),
    ("H2", "6-31g*", "RHF", "H 0.0 0.0 0.0; H 0.0 0.0 0.74", "Angstrom"),
    ("H2O", "sto-3g", "RHF",
     "O 0.0 0.0 0.117370; H 0.0 0.756950 -0.469483; H 0.0 -0.756950 -0.469483",
     "Angstrom"),
    ("H2O", "sto-3g", "B3LYP",
     "O 0.0 0.0 0.117370; H 0.0 0.756950 -0.469483; H 0.0 -0.756950 -0.469483",
     "Angstrom"),
    ("H2O", "6-31g*", "RHF",
     "O 0.0 0.0 0.117370; H 0.0 0.756950 -0.469483; H 0.0 -0.756950 -0.469483",
     "Angstrom"),
]


def lowdin_populations(mol, dm):
    """True Löwdin populations: q_A = sum_{m in A} (S^(1/2) D S^(1/2))_{mm}.

    Returns array of populations per AO, summed by atom -> charges.
    """
    s = mol.intor_symmetric("int1e_ovlp")
    s_half = scipy.linalg.sqrtm(s).real
    dm_lowdin = s_half @ dm @ s_half
    pop_per_ao = np.diag(dm_lowdin)
    # Sum populations by atom
    n_atom = mol.natm
    atom_pop = np.zeros(n_atom)
    ao_atom = mol.ao_labels(fmt=None)  # list of (atom_idx, atom_sym, ao_label)
    for i, lbl in enumerate(ao_atom):
        atom_pop[lbl[0]] += pop_per_ao[i]
    # Charge = nuclear_charge - electron_population
    z = np.array([mol.atom_charge(i) for i in range(n_atom)], dtype=float)
    return z - atom_pop, atom_pop


def main() -> None:
    results = []
    for molname, basis, method, atom_str, unit in SYSTEMS:
        mol = gto.Mole()
        mol.atom = atom_str
        mol.unit = unit
        mol.basis = basis
        mol.cart = True
        mol.verbose = 0
        mol.build()

        if method == "RHF":
            mf = scf.RHF(mol)
        elif method == "B3LYP":
            mf = dft.RKS(mol)
            mf.xc = "b3lyp5"
        mf.conv_tol = 1e-12
        mf.kernel()

        dm = mf.make_rdm1()
        atoms = [mol.atom_symbol(i) for i in range(mol.natm)]

        # Mulliken (PySCF built-in)
        mull_pop, mull_chg = mf.mulliken_pop(verbose=0)

        # True Löwdin
        low_chg, low_pop = lowdin_populations(mol, dm)

        results.append({
            "molecule": molname,
            "basis": basis,
            "method": method,
            "atoms": atoms,
            "mulliken_charges": [float(c) for c in mull_chg],
            "lowdin_charges": [float(c) for c in low_chg],
            "lowdin_populations_per_atom": [float(p) for p in low_pop],
            "charge_sum_check": float(sum(low_chg)),
        })

    out = {
        "generated": "2026-04-15T11:30:00",
        "method_note": (
            "Löwdin charges = Z_nuc - sum(diag(S^(1/2) D S^(1/2))) per atom; "
            "Mulliken charges via PySCF mulliken_pop()."
        ),
        "results": results,
    }
    path = Path(__file__).resolve().parents[2] / "tests/golden/dft/population_validation.json"
    path.write_text(json.dumps(out, indent=2))
    print(f"Wrote {path}")
    for r in results:
        sys = f"{r['molecule']}/{r['basis']}/{r['method']}"
        print(f"  {sys:25s}  Mull: {r['mulliken_charges']}  Löwdin: {r['lowdin_charges']}")


if __name__ == "__main__":
    main()
