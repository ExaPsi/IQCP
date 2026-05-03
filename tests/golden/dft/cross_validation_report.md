# DFT Cross-Validation Report: IQCP vs PySCF

**Date:** 2026-03-19 (revised after integral engine fix + geometry alignment)
**IQCP Version:** qc-core 0.1.0 (feature/phase3)
**PySCF Version:** 2.11.0

---

## Integral Engine Fix

Two bugs fixed in the on-the-fly integral engine:
1. **Inaccurate `erf()`** — Custom Cody approximation (~7 digits) replaced with `libm::erf()` (machine precision)
2. **Nuclear VRR cancellation** — Auxiliary-index VRR replaced with Rys-quadrature approach (following libcint)

Result: Max H_core error improved from **7.6e-6** to **1.6e-9** (4600x improvement).

## Geometry Alignment

All PySCF reference values generated using the **exact same coordinates** as IQCP Rust test helpers. Previous discrepancies (up to 1e-3 Ha) were caused by coordinate mismatches.

## Conventions

- **STO-3G:** No d-functions (Cartesian = Spherical)
- **6-31G*:** IQCP uses **Cartesian** d-functions. PySCF references use `mol.cart=True`.

## Cross-Validation: 30/30 PASS at 1e-5 Ha

### STO-3G

| Molecule | Method | PySCF (Ha) | Tolerance | Status |
|----------|--------|-----------|-----------|--------|
| H2 | RHF | -1.11671432506255 | <1e-5 | PASS |
| H2 | LDA | -1.12120070415970 | <1e-5 | PASS |
| H2 | B3LYP | -1.15860014805365 | <1e-5 | PASS |
| H2O | RHF | -74.96302571754660 | <1e-5 | PASS |
| H2O | LDA | -74.73203834615985 | <1e-5 | PASS |
| H2O | B3LYP | -75.27523821491260 | <1e-5 | PASS |
| NH3 | RHF | -55.45436165848794 | <1e-5 | PASS |
| NH3 | LDA | -55.29075298644806 | <1e-5 | PASS |
| NH3 | B3LYP | -55.74932495255848 | <1e-5 | PASS |
| HF | RHF | -98.57077532424054 | <1e-5 | PASS |
| HF | LDA | -98.24587315661505 | <1e-5 | PASS |
| HF | B3LYP | -98.88126529002626 | <1e-5 | PASS |
| CH4 | RHF | -39.72682902478177 | <1e-5 | PASS |
| CH4 | LDA | -39.61679835214327 | <1e-5 | PASS |
| CH4 | B3LYP | -40.00251270098887 | <1e-5 | PASS |

### 6-31G* (Cartesian d-functions)

| Molecule | Method | PySCF (Ha) | Tolerance | Status |
|----------|--------|-----------|-----------|--------|
| H2 | RHF | -1.12674270445184 | <1e-5 | PASS |
| H2 | LDA | -1.13266912681730 | <1e-5 | PASS |
| H2 | B3LYP | -1.16871298976056 | <1e-5 | PASS |
| H2O | RHF | -76.01050569647504 | <1e-5 | PASS |
| H2O | LDA | -75.84438175789660 | <1e-5 | PASS |
| H2O | B3LYP | -76.37159375341544 | <1e-5 | PASS |
| NH3 | RHF | -56.18399100969380 | <1e-5 | PASS |
| NH3 | LDA | -56.06153732874716 | <1e-5 | PASS |
| NH3 | B3LYP | -56.51134161643392 | <1e-5 | PASS |
| HF | RHF | -100.00286325730460 | <1e-5 | PASS |
| HF | LDA | -99.76633158746822 | <1e-5 | PASS |
| HF | B3LYP | -100.38211096579278 | <1e-5 | PASS |
| CH4 | RHF | -40.19515375820445 | <1e-5 | PASS |
| CH4 | LDA | -40.09718983409379 | <1e-5 | PASS |
| CH4 | B3LYP | -40.48225959744568 | <1e-5 | PASS |

## Error Budget

| Source | Contribution |
|--------|-------------|
| One-electron integrals | ~1e-9 Ha (after libm::erf fix) |
| Two-electron integrals | ~1e-10 Ha |
| DFT grid quadrature | ~1e-6 to 5e-6 Ha (dominant) |
| **Total** | **~5e-6 Ha** |

The residual ~5e-6 error is from the Becke grid quadrature (Mura-Knowles 75-point radial + Lebedev-302 angular + SG-1 pruning). This is the expected accuracy for a standard-quality DFT grid.
