# Gaussian 09 PolarDeriv Reference Data

Captured 2026-04-11 from Gaussian 09 Rev A.02 (Gaussian 09 (local installation)).

## Purpose

These files contain the **analytical polarizability derivatives** computed by Gaussian's
`DPolar` subroutine (`l1002.F` line 6240+) for the four Raman test molecules. They serve
as a **machine-precision oracle** for cross-validating the IQCP semi-analytical
polarizability derivative implementation.

## Files

| Input | Output | Molecule | Geometry | Basis | Method |
|-------|--------|----------|----------|-------|--------|
| `h2_test.gjf` | `h2_test.out` | H₂ | (0,0,±0.37 Å) | STO-3G | RHF Freq=Raman |
| `h2o_test.gjf` | `h2o_test.out` | H₂O | matches IQCP standard geometry | STO-3G | RHF Freq=Raman |
| `ch4_test.gjf` | `ch4_test.out` | CH₄ | tetrahedral, R(C-H)=1.089 Å | STO-3G | RHF Freq=Raman |
| `co2_test.gjf` | `co2_test.out` | CO₂ | linear, R(C=O)=1.162 Å | STO-3G | RHF Freq=Raman |

## How to extract `PolarDeriv`

The `PolarDeriv` array is in the **archive line** at the bottom of the Gaussian output:

```python
import re
def parse_gaussian_polarderiv(path):
    """Returns numpy array of shape (n_atoms, 3, 3, 3) indexed [A, gamma, d, e]."""
    import numpy as np
    with open(path) as f:
        text = f.read()
    m = re.search(r'1\\1\\.*?\\\\@', text, re.DOTALL)
    archive = m.group(0).replace('\n ', '').replace('\n', '')
    m = re.search(r'PolarDeriv=([^\\]+)', archive)
    values = [float(x) for x in m.group(1).split(',')]
    n_atoms = len(values) // 18
    result = np.zeros((n_atoms, 3, 3, 3))
    pair_idx = [(0,0),(0,1),(1,1),(0,2),(1,2),(2,2)]
    idx = 0
    for A in range(n_atoms):
        for gamma in range(3):
            for d, e in pair_idx:
                v = values[idx]; idx += 1
                result[A, gamma, d, e] = v
                result[A, gamma, e, d] = v
    return result
```

## Key reference values

### H₂/STO-3G

- `∂α_zz/∂R_{H1,z} = +3.4917361 Ha/bohr` (the only nonzero element)
- IQCP semi-analytical FD-of-CPHF result: `+3.4917363` (agreement to 2e-7)

### H₂O/STO-3G

| Element | Value (Ha/bohr) |
|---------|-----------------|
| `∂α_zz/∂R_{O,z}` | +4.5298984 |
| `∂α_yy/∂R_{O,z}` | +2.9577730 |
| `∂α_yz/∂R_{O,y}` | +2.5607138 |
| `∂α_xz/∂R_{O,x}` | +1.8816255 |
| `∂α_yy/∂R_{H1,y}` | +4.2895547 |
| `∂α_zz/∂R_{H1,z}` | -2.2649492 |
| `∂α_yz/∂R_{H1,z}` | +2.0688887 |

(All 54 H₂O values are in `h2o_test.out` archive line.)

### CO₂/STO-3G

- `∂α_zz/∂R_{O1,z} = +7.6621102 Ha/bohr`
- `∂α_xz/∂R_{O1,x} = +2.7592642`
- `∂α_yz/∂R_{O1,y} = +2.7592642`

## Translational invariance check

For all four molecules: `Σ_atoms ∂α/∂R_atom = 0` to ~1e-7. This is a sanity check
for the analytical implementation.

## Status

These reference values are **NOT currently used** by the IQCP test suite. The
existing golden data in `tests/golden/raman/{h2,h2o,ch4,co2}_sto3g_rhf.json` is
generated from PySCF FD-of-analytical-polarizability (see
`scripts/phase5/generate_raman_golden.py`), which agrees with Gaussian to ~1e-7.

The Gaussian reference is preserved here for **future use** when an analytical
formula for `∂α/∂R` is implemented (see Appendix H of `docs/stories/US-098`).

## How Gaussian computes these values (`l1002.F:DPolar`)

Gaussian's `DPolar` subroutine combines TWO contributions:

1. **`DPolDs` (line 6286)**: Trace of field-perturbed density `D^(Ed)` with the
   nuclear-derivative dipole integral matrix. Equivalent to Amos eq (16) first
   term.

2. **`CPTrFx` + `CPTrac` (lines 6304-6345)**: Contraction of the **second-order
   density `P²A` and "W²" matrix** with bare dipole / overlap matrices. The P²
   and W² matrices come from Gaussian's **second-order CPHF solver** (`Force`
   link's "2nd order cphf"; see line 154).

The IQCP analytical formula attempts (documented in US-098 Appendix G/H) failed
because they tried to use the Wigner 2n+1 symmetric formula (Amos eq 12) instead
of explicit second-order CPHF. Future analytical implementation should follow
Gaussian's two-part approach.
