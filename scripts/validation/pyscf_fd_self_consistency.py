#!/usr/bin/env python3
"""
PySCF analytical-vs-finite-difference gradient self-consistency.

Reproduces the values cited in M2 manuscript §4.4 (Gradient Validation).

Usage:
    source .venv/bin/activate
    python scripts/validation/pyscf_fd_self_consistency.py

Reference geometry (H2O, validation):
    O  0.0  0.0  0.22166487441148175  bohr
    H  0.0  +1.4309006215206648  -0.886659497645927  bohr
    H  0.0  -1.4309006215206648  -0.886659497645927  bohr

PySCF version: 2.11.0
NumPy: 2.4.1, SciPy: 1.17.0
Step size: h = 1e-5 bohr (central difference)
SCF tolerance: conv_tol = 1e-12
"""
import numpy as np
from pyscf import gto, scf, dft

H_STEP = 1e-5  # bohr
GEOM = """
O  0.0  0.0   0.22166487441148175
H  0.0  1.4309006215206648  -0.886659497645927
H  0.0 -1.4309006215206648  -0.886659497645927
"""


def make_mf(method: str, mol):
    if method == "RHF":
        mf = scf.RHF(mol)
    elif method == "LDA":
        mf = dft.RKS(mol)
        mf.xc = "lda,vwn"
    elif method == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "b3lyp5"
    else:
        raise ValueError(method)
    mf.conv_tol = 1e-12
    mf.verbose = 0
    return mf


def fd_gradient(method: str, mol):
    """Central-difference gradient of the SCF energy."""
    natm = mol.natm
    g = np.zeros((natm, 3))
    coords = mol.atom_coords().copy()
    for ia in range(natm):
        for ix in range(3):
            for sign in (+1, -1):
                disp = coords.copy()
                disp[ia, ix] += sign * H_STEP
                mp = mol.copy()
                mp.set_geom_(disp, unit="Bohr")
                mp.build()
                e = make_mf(method, mp).kernel()
                g[ia, ix] += sign * e
            g[ia, ix] /= 2 * H_STEP
    return g


def main() -> None:
    print(f"# PySCF FD self-consistency (h = {H_STEP} bohr)\n")
    mol = gto.Mole()
    mol.atom = GEOM
    mol.unit = "Bohr"
    mol.basis = "sto-3g"
    mol.cart = True
    mol.verbose = 0
    mol.build()

    print(f"{'method':<8s}  {'max|analytic|':>15s}  {'max|FD|':>13s}  {'max|delta|':>14s}")
    for method in ("RHF", "LDA", "B3LYP"):
        mf = make_mf(method, mol)
        mf.kernel()
        g_ana = mf.nuc_grad_method().kernel()
        g_fd = fd_gradient(method, mol)
        delta = np.abs(g_ana - g_fd)
        print(
            f"{method:<8s}  {np.max(np.abs(g_ana)):15.6e}  "
            f"{np.max(np.abs(g_fd)):13.6e}  {np.max(delta):14.6e}"
        )


if __name__ == "__main__":
    main()
