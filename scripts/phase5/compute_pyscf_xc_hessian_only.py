#!/usr/bin/env python3
"""
Extract just the XC contribution to the PySCF Hessian on IQCP's exact grid.

Computes both `_get_vxc_diag` (Part 1 = veff_diag) and `_get_vxc_deriv2`
(Parts 2+3 = ipip + fxc bilinear), and prints them as Rust literals for
comparison against IQCP's `xc_hessian_gga()` output piece-by-piece.
"""
import json
import sys

import numpy as np

import pyscf
from pyscf import gto, dft
from pyscf.hessian import rks as rks_hess


def dump_array(name, arr):
    n3 = arr.shape[0]
    print(f"#[rustfmt::skip]")
    print(f"let {name}: [f64; 81] = [")
    for i in range(n3):
        row = "    "
        for j in range(n3):
            row += f"{arr[i,j]: 22.15e},"
            if j < n3 - 1:
                row += " "
        print(row)
    print("];")


def main():
    functional = sys.argv[1] if len(sys.argv) > 1 else "b3lyp5"

    print(f"# PySCF version: {pyscf.__version__}", file=sys.stderr)
    print(f"# Functional: {functional}", file=sys.stderr)

    with open("/tmp/iqcp_grid_h2o.json") as f:
        grid_data = json.load(f)

    coords = np.array(grid_data["points"], dtype=np.float64)
    weights = np.array(grid_data["weights"], dtype=np.float64)

    mol = gto.Mole()
    mol.atom = [
        ["O", (0.0, 0.0, 0.0)],
        ["H", (0.0, 1.43, 1.11)],
        ["H", (0.0, -1.43, 1.11)],
    ]
    mol.basis = "sto-3g"
    mol.unit = "bohr"
    mol.cart = True
    mol.verbose = 0
    mol.build()

    mf = dft.RKS(mol)
    mf.xc = functional
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.max_cycle = 200

    # Assign IQCP's grid
    mf.grids.coords = coords
    mf.grids.weights = weights
    mf.grids.non0tab = None

    def noop_build(*args, **kwargs):
        return mf.grids
    mf.grids.build = noop_build

    mf.kernel()
    print(f"# E_scf = {mf.e_tot:.12f} Ha", file=sys.stderr)

    hobj = mf.Hessian()
    hobj.grids = mf.grids

    # Part 1: veff_diag
    veff_diag = rks_hess._get_vxc_diag(hobj, mf.mo_coeff, mf.mo_occ, 4000)
    # Shape (3,3,nao,nao)
    nao = mol.nao
    dm0 = mf.make_rdm1()
    n_atoms = mol.natm
    n3 = 3 * n_atoms
    aoslices = mol.aoslice_by_atom()

    # Assemble de2_part1 from veff_diag. Only diagonal atom blocks.
    # PySCF line 106: `de2[i0,i0] += einsum('xypq,pq->xy', veff_diag[:,:,p0:p1], dm0[p0:p1])*2`
    # The slice veff_diag[:,:,p0:p1] slices the 3rd axis (index 2 of (3,3,nao,nao)),
    # which is the "p" axis corresponding to mu. So mu ∈ atom ia.
    de2_part1 = np.zeros((n3, n3))
    for ia in range(n_atoms):
        p0, p1 = aoslices[ia][2:]
        block = np.einsum("xypq,pq->xy", veff_diag[:, :, p0:p1, :], dm0[p0:p1, :]) * 2
        for d1 in range(3):
            for d2 in range(3):
                de2_part1[3 * ia + d1, 3 * ia + d2] = block[d1, d2]

    # Part 2+3: _get_vxc_deriv2
    vxc = rks_hess._get_vxc_deriv2(hobj, mf.mo_coeff, mf.mo_occ, 4000)
    # Shape (natm, 3, 3, nao, nao)
    # PySCF indexes `vxc[ia,:,:,:,q0:q1]` with p as first AO axis; when
    # contracted with dm0[q0:q1] (mu ∈ ja as bra, nu as ket), the einsum is
    #   einsum('xypq,pq->xy', vxc[ia,:,:,:,q0:q1], dm0[q0:q1])
    # Note: PySCF's `vxc[ia]` has shape (3,3,nao,nao), where the LAST axis
    # is what `q0:q1` is applied to — but then dm0 is indexed as dm0[q0:q1]
    # = [mu in ja, nu]. So p=mu, q=nu, with the einsum contracting (mu,nu).
    # Let me replicate line 109 of PySCF's partial_hess_elec exactly:
    #   de2[i0,j0] += einsum('xypq,pq->xy', veff[:,:,q0:q1], dm0[q0:q1])*2
    # Here `veff[:,:,q0:q1]` slices the THIRD axis (index 2 of (3,3,nao,nao)),
    # i.e. the "p" axis, so p ∈ [q0,q1) = atom ja. Then p is mu ∈ ja, and the
    # full dm0 row gives us all q=nu. Einsum contracts (p,q).
    de2_part23 = np.zeros((n3, n3))
    for ia in range(n_atoms):
        for ja in range(n_atoms):
            q0, q1 = aoslices[ja][2:]
            # vxc[ia,:,:,p_in_ja,:] has shape (3,3,p1-p0,nao)
            # dm0[p_in_ja,:] has shape (p1-p0,nao)
            block = np.einsum("xypq,pq->xy", vxc[ia, :, :, q0:q1, :], dm0[q0:q1, :]) * 2
            for d1 in range(3):
                for d2 in range(3):
                    de2_part23[3 * ia + d1, 3 * ja + d2] = block[d1, d2]

    # Symmetrize Part 1 (it's diagonal blocks only, but PySCF relies on the
    # j<=i assembly for symmetry; here we compute both triangles). Part 23
    # is asymmetric and gets symmetrized in de2 = 0.5*(de2 + de2.T).
    # For direct comparison against IQCP's assembled H_xc, PySCF's final
    # symmetrization is handled in partial_hess_elec. Replicate here:
    total = de2_part1 + de2_part23
    # PySCF symmetrizes via `for j0 in range(i0): de2[j0,i0] = de2[i0,j0].T`
    # i.e., it copies the lower triangle to the upper. Then 0.5*(de2+de2.T).
    # The net effect is that the block for (i0, j0) and (j0, i0) are the
    # transpose of each other.

    dump_array("pyscf_xc_part1", de2_part1)
    print()
    dump_array("pyscf_xc_parts23", de2_part23)
    print()
    dump_array("pyscf_xc_total", total)

    print(f"\n# Max |Part 1|: {np.abs(de2_part1).max():.6e}", file=sys.stderr)
    print(f"# Max |Parts 2+3|: {np.abs(de2_part23).max():.6e}", file=sys.stderr)
    print(f"# Max |Total XC|: {np.abs(total).max():.6e}", file=sys.stderr)


if __name__ == "__main__":
    main()
