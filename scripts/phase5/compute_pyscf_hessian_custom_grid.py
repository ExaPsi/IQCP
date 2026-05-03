#!/usr/bin/env python3
"""
Compute the PySCF DFT Hessian for H2O/STO-3G on IQCP's exact custom grid.

This script loads IQCP's exported Becke grid (points + weights) and assigns
it to PySCF's Hessian object via `mf.grids.coords` / `mf.grids.weights`. With
identical grids, remaining differences between IQCP and PySCF are limited to
SCF convergence + XC kernel evaluation precision (well below 1e-6).

Usage:
    source .venv/bin/activate
    python scripts/phase5/compute_pyscf_hessian_custom_grid.py <functional>

Where <functional> is 'lda', 'b3lyp', or 'b3lyp5'. Defaults to b3lyp5.

Input:  /tmp/iqcp_grid_h2o.json  (produced by `cargo test export_h2o_lda_grid_for_pyscf`)
Output: prints 9x9 Hessian (row-major) to stdout as a Rust `[f64; 81]` literal
"""
import json
import sys

import numpy as np

import pyscf
from pyscf import gto, dft


def main():
    functional = sys.argv[1] if len(sys.argv) > 1 else "b3lyp5"
    # IQCP 'lda' = Slater exchange + VWN5 correlation. PySCF's 'lda' is Slater only.
    # Translate to the PySCF name that matches IQCP's implementation.
    if functional.lower() == "lda":
        functional = "lda,vwn5"

    print(f"# PySCF version: {pyscf.__version__}", file=sys.stderr)
    print(f"# Functional: {functional}", file=sys.stderr)

    with open("/tmp/iqcp_grid_h2o.json") as f:
        grid_data = json.load(f)

    coords = np.array(grid_data["points"], dtype=np.float64)
    weights = np.array(grid_data["weights"], dtype=np.float64)
    print(
        f"# Loaded IQCP grid: {grid_data['n_points']} points "
        f"(coords.shape={coords.shape}, weights.shape={weights.shape})",
        file=sys.stderr,
    )

    # Build H2O molecule - must match IQCP test geometry EXACTLY
    mol = gto.Mole()
    mol.atom = [
        ["O", (0.0, 0.0, 0.0)],
        ["H", (0.0, 1.43, 1.11)],
        ["H", (0.0, -1.43, 1.11)],
    ]
    mol.basis = "sto-3g"
    mol.unit = "bohr"
    mol.cart = True  # Cartesian basis (6d), matches IQCP
    mol.verbose = 0
    mol.build()

    # DFT setup
    mf = dft.RKS(mol)
    mf.xc = functional
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.max_cycle = 200

    # Override the DFT grid with IQCP's grid. PySCF's grid object is
    # mf.grids which has .coords and .weights attributes.
    mf.grids.coords = coords
    mf.grids.weights = weights
    # Disable PySCF's own grid build by setting level to -1 semantics:
    # Actually, once coords/weights are set, PySCF will use them as long as
    # we don't call grids.build() again. We can mark it "already built" by
    # setting non0tab = None and atom_grids_tab to empty.
    mf.grids.non0tab = None

    # Monkey-patch build() so it's a no-op (prevents regeneration on reset).
    original_build = mf.grids.build

    def noop_build(*args, **kwargs):
        return mf.grids

    mf.grids.build = noop_build

    print("# Running SCF...", file=sys.stderr)
    e_scf = mf.kernel()
    print(f"# SCF converged: {mf.converged}", file=sys.stderr)
    print(f"# E_scf = {e_scf:.12f} Ha", file=sys.stderr)
    print(f"# Number of SCF cycles: {mf.cycles}", file=sys.stderr)

    # Compute analytical Hessian with the custom grid
    print("# Computing analytical Hessian...", file=sys.stderr)
    hessobj = mf.Hessian()
    # Make sure the Hessian uses the same custom grid
    hessobj.grids = mf.grids
    hess = hessobj.kernel()

    # hess has shape (natm, natm, 3, 3); we need 9x9
    natm = mol.natm
    n3 = 3 * natm
    h9 = np.zeros((n3, n3))
    for i in range(natm):
        for j in range(natm):
            for a in range(3):
                for b in range(3):
                    h9[3 * i + a, 3 * j + b] = hess[i, j, a, b]

    # Print as Rust literal
    print(f"// PySCF reference Hessian for H2O/STO-3G, functional={functional}")
    print(f"// PySCF version: {pyscf.__version__}")
    print(f"// Grid: IQCP Becke grid ({grid_data['n_points']} points, loaded from /tmp/iqcp_grid_h2o.json)")
    print(f"// SCF energy: {e_scf:.12f} Ha")
    print("#[rustfmt::skip]")
    print(f"let pyscf_{functional.replace('5','')}: [f64; 81] = [")
    for i in range(n3):
        row = "    "
        for j in range(n3):
            row += f"{h9[i,j]: 18.12f},"
            if j < n3 - 1:
                row += " "
        print(row)
    print("];")

    print(f"\n# Max |Hessian|: {np.abs(h9).max():.6e}", file=sys.stderr)
    print(f"# Hessian symmetry check: {np.abs(h9 - h9.T).max():.2e}", file=sys.stderr)


if __name__ == "__main__":
    main()
