# Reproducibility Guide

This document is the reviewer-facing recipe for reproducing every numerical
result reported in the IQCP JCIM Full Article. Every golden reference in
`tests/golden/` was generated from the canonical reference implementations
listed below, and every assertion in the Rust test suite is a direct
comparison against those goldens.

If you can run `cargo test --workspace` and `pytest`-equivalent PySCF scripts,
you can reproduce the entire validation suite.

---

## 1. Computing Environment

### Required toolchains

| Tool         | Minimum version | Purpose                                            |
|--------------|-----------------|----------------------------------------------------|
| Rust         | 1.94 stable     | `qc-core`, `qc-wasm`, `qc-io` workspace            |
| Cargo        | bundled         | Build / test / bench driver                        |
| `wasm-pack`  | 0.13+           | Builds `qc-wasm` for the browser                   |
| `wasm32-unknown-unknown` target | via `rustup target add` | WASM compilation         |
| Node         | 18 LTS or newer | Vite, frontend build, integration tests            |
| npm          | 9+              | JS dependency resolution                           |
| Python       | 3.12.x          | PySCF reference & golden regeneration scripts      |
| `gzip`       | system          | Bundle-size measurement                            |

### Reference (PySCF) environment

All `tests/golden/` data and every cross-code validation script in
`scripts/validation/`, `scripts/phase5/`, and `tests/benchmarks/pyscf/` was
produced under the following pinned stack. Reproducing within these versions
is required for bit-identical reference matches; small package-version drift
will produce equivalent science but may shift trailing digits.

| Package | Version |
|---------|---------|
| Python  | 3.12.3  |
| PySCF   | 2.11.0  |
| NumPy   | 2.4.1 (some Phase 5 goldens were regenerated with NumPy 1.26.4 — see per-file `_provenance`) |
| SciPy   | 1.17.0  |
| h5py    | 3.15.1  |

The exact toolchain used to generate each golden is recorded in the file's
`metadata` / `_provenance` block (see Section 5).

---

## 2. PySCF Environment Setup

The PySCF environment is *not* checked in. Recreate it from a clean shell:

```bash
cd /path/to/IQCP

# Create and activate the virtual environment
python3.12 -m venv .venv
source .venv/bin/activate

# Pin to the exact reference stack
pip install --upgrade pip
pip install \
    pyscf==2.11.0 \
    numpy==2.4.1 \
    scipy==1.17.0 \
    h5py==3.15.1
```

Verify:

```bash
python -c "import pyscf, numpy, scipy, h5py; \
print('PySCF', pyscf.__version__); \
print('NumPy', numpy.__version__); \
print('SciPy', scipy.__version__); \
print('h5py', h5py.__version__)"
```

Expected output:

```
PySCF 2.11.0
NumPy 2.4.1
SciPy 1.17.0
h5py 3.15.1
```

Single-threaded execution is enforced inside the benchmark scripts via
`pyscf.lib.num_threads(1)` — there is no need to set OMP/BLAS env vars
manually.

---

## 3. Quick Validation (~5 minutes)

This pipeline confirms the build is healthy and the smallest reference case
(H₂ / STO-3G / RHF) matches PySCF to sub-microhartree.

```bash
# 1. Build the Rust workspace
cargo build --workspace --release

# 2. Run the Boys + Rys + smallest SCF golden tests
cargo test -p qc-core --release boys
cargo test -p qc-core --release rys
cargo test -p qc-core --release scf::tests::h2

# 3. Cross-check against PySCF (requires .venv from Section 2)
source .venv/bin/activate
python -c "
from pyscf import gto, scf, lib
lib.num_threads(1)
mol = gto.M(atom='H 0 0 0; H 0 0 1.4', basis='sto-3g', unit='bohr')
mf = scf.RHF(mol); mf.conv_tol = 1e-12
print(f'PySCF H2/STO-3G: {mf.kernel():.12f} Ha')
"
```

The PySCF energy must agree with the IQCP value cached in
`tests/golden/orbital/h2_sigma_g.json` (`-1.116714325062551 Ha`) to all 12
digits.

---

## 4. Full Validation (~30 minutes)

The full Rust suite covers every category listed in Section 5:

```bash
# Full Rust workspace (golden + property + integration tests)
cargo test --workspace --release

# WASM module build + browser bindings sanity
wasm-pack build crates/qc-wasm --release --target web \
    --out-dir ../../apps/web/src/wasm

# Frontend type check + lint
cd apps/web
npm install
npm run lint
npm run typecheck
cd -
```

### Per-category recipe table

Every entry below is a self-contained validation: `Goldens path` is the
checked-in reference, `Generator script` regenerates it from PySCF / SciPy,
and `Tolerance` is the assertion used in the Rust test suite. **Do not
regenerate goldens** unless an algorithm has materially changed and the
change has been independently audited.

| Category    | Goldens path                                  | Generator script                                            | Reference                                  | Tolerance                |
|-------------|-----------------------------------------------|-------------------------------------------------------------|--------------------------------------------|--------------------------|
| Boys        | `tests/golden/boys/reference.json`            | uses `scipy.special.hyp1f1`                                 | Shavitt (1963); SciPy 1.17.0               | 1e-12 absolute           |
| Rys         | `tests/golden/rys/reference.json`             | `scripts/generate_golden_rys.py`                            | Dupuis–Rys–King (1976); libcint            | roots ∈ (0,1); reconstruct 1e-10 |
| ERI (primitive) | `tests/golden/eri/c2h4_631gs_verification_report.md` | `scripts/generate_c2h4_eri_reference.py`, `scripts/verify_c2h4_comprehensive.py` | PySCF 2.11.0                  | 1e-10 per ERI            |
| Integrals (1e) | `tests/golden/integrals/h2*_sto3g_overlap.json` | PySCF `mol.intor('int1e_ovlp')`                          | PySCF 2.11.0                               | 1e-12                    |
| SCF / DFT (energy) | `tests/golden/dft/validation.json`, `validation_all_bases.json` | covered by `scripts/validation/` and `tests/benchmarks/pyscf/benchmark_scf.py` | PySCF 2.11.0 (b3lyp5, cart=True, conv_tol=1e-12) | 1e-7 Ha (1e-5 Ha for cc-pVDZ — see notes) |
| Gradients   | `tests/golden/dft/gradient_validation*.json`  | regenerated via PySCF `mf.nuc_grad_method().kernel()`       | PySCF 2.11.0                               | 1e-7 Ha/bohr (RHF/B3LYP) |
| Geometry    | `tests/golden/dft/geometry_validation.json`   | PySCF `geometric` + `mol.RHF()`                             | PySCF 2.11.0                               | 1e-4 bohr / 1e-6 Ha      |
| Population  | `tests/golden/dft/population_validation.json` | PySCF `mulliken_pop()` + Löwdin from S^(1/2)·D·S^(1/2)      | PySCF 2.11.0                               | 1e-6 e                   |
| Orbital     | `tests/golden/orbital/*.json`                 | PySCF MO coeffs + on-grid evaluation                        | PySCF 2.11.0                               | 1e-8 at grid points      |
| PES         | `tests/golden/pes/*.json`                     | `tests/benchmarks/pyscf/pes_*_h2o_*.py`                     | PySCF 2.11.0                               | 1e-8 Ha per point        |
| IR          | `tests/golden/ir/*.json`                      | `scripts/phase5/generate_ir_golden.py`                      | PySCF 2.11.0 + analytic dipole derivs      | 1e-6 km/mol intensity    |
| Raman       | `tests/golden/raman/*.json`                   | `scripts/phase5/generate_raman_golden.py`                   | PySCF 2.11.0 CPHF + FD polarizability      | semi-analytical match (see §7) |
| Thermo      | `tests/golden/thermo/*.json`                  | `scripts/phase5/generate_thermo_golden.py`                  | PySCF 2.11.0 RRHO partition functions      | 1e-7 Ha / 1e-3 J·K⁻¹·mol⁻¹ |
| Thermochem  | `tests/golden/thermochem/*.json`              | `scripts/phase5/generate_thermochem_golden.py`              | PySCF 2.11.0 + RRHO @ 298.15 K, 101325 Pa  | 1e-6 Ha (ZPE, H, S, G)   |
| Phase 5 full | `tests/golden/phase5/*.json`                 | `scripts/phase5/generate_phase5_golden.py`                  | PySCF 2.11.0 (RHF/LDA/B3LYP, conv_tol=1e-12) | per-field (see file)   |

Note: `tests/golden/scf/` is reserved for raw SCF iteration traces and is
currently empty in the public release; SCF energy convergence is exercised
through the `dft/`, `phase5/`, and `pes/` goldens above.

---

## 5. Cross-Code Validation against PySCF

The cross-code scripts demonstrate that IQCP and PySCF agree on the same
science from independent code paths.

### `scripts/validation/`

| Script | What it verifies | How to run |
|--------|------------------|------------|
| `pyscf_fd_self_consistency.py` | PySCF analytic gradients vs PySCF finite-difference (sanity of the reference itself) | `python scripts/validation/pyscf_fd_self_consistency.py` |
| `pyscf_geometry_optimization.py` | PySCF-side reference optimized geometries that IQCP must reproduce | `python scripts/validation/pyscf_geometry_optimization.py` |
| `pyscf_population_analysis.py` | Mulliken + Löwdin charges from a PySCF SCF density | `python scripts/validation/pyscf_population_analysis.py` |

Expected output: each script prints a tabular comparison with per-atom /
per-element residuals and exits 0 if all residuals are within the documented
tolerance.

### `scripts/phase5/`

Re-runs the Phase 5 reference data generation (IR, Raman, thermo, full
Phase 5 RHF/LDA/B3LYP). Each generator writes back into
`tests/golden/<category>/`. For reviewers wanting to compare without
overwriting checked-in files, redirect output via the `--out` flag if the
script exposes one; otherwise diff against `git`.

```bash
source .venv/bin/activate
python scripts/phase5/generate_ir_golden.py
python scripts/phase5/generate_raman_golden.py
python scripts/phase5/generate_thermo_golden.py
python scripts/phase5/generate_thermochem_golden.py
python scripts/phase5/generate_phase5_golden.py

# Hessian helpers (used by IR/Raman intermediate verification)
python scripts/phase5/compute_pyscf_hessian_custom_grid.py
python scripts/phase5/compute_pyscf_xc_hessian_only.py
```

### `tests/benchmarks/pyscf/`

These scripts produce PySCF timing references *and* energy/gradient cross-
checks for the JCIM Performance section (see `BENCHMARKS.md`). To run only
the validation portion (energies, not timings), each script's `--validate`
mode (where exposed) skips the timing loop.

```bash
source .venv/bin/activate
bash tests/benchmarks/pyscf/run_all.sh          # full suite
python tests/benchmarks/pyscf/same_grid_validation.py   # DFT same-grid energy match
python tests/benchmarks/pyscf/same_grid_ch4_c6h6.py     # extended same-grid set
python tests/benchmarks/pyscf/investigate_rhf_discrepancy.py  # diagnostic
python tests/benchmarks/pyscf/pes_rigid_h2o_bond.py           # PES rigid scan
python tests/benchmarks/pyscf/pes_rigid_h2o_angle.py
python tests/benchmarks/pyscf/pes_relaxed_h2o_angle.py
```

`tests/benchmarks/pyscf/results/` holds the cached output of the most recent
run; the JCIM tables are populated from these files.

---

## 6. Determinism

IQCP is designed to produce **bit-identical** outputs across all WASM-capable
browsers and across native Rust on x86_64 / aarch64.

| Layer | Determinism mechanism |
|-------|------------------------|
| Boys / Rys | Pure functions; series and recurrences are evaluated in fixed order. No reductions over hash maps. |
| ERI / SCF  | Iteration order over shell pairs and AO indices is deterministic; no parallel reductions in the hot path. |
| Linear algebra | All matrix operations use deterministic in-house implementations (no LAPACK BLAS dispatch); double precision throughout. |
| WASM | Single-threaded; IEEE-754 double precision per the WebAssembly spec; no SIMD-induced fused-multiply-add (FMA) reordering. |

For cross-browser regression testing, the `apps/web/` integration suite
captures `Float64Array` payloads from each WASM call and compares them
byte-for-byte against the Chromium reference output.

---

## 7. Known Differences

Not every reference is matched to machine precision. The following
deviations are **expected and load-bearing** — they reflect either
methodological differences between IQCP and PySCF or the intrinsic limit of
finite-difference cross-validation.

| Quantity | Expected residual | Reason |
|----------|-------------------|--------|
| DFT (different XC grids) | up to ~1e-5 Ha | IQCP and PySCF use different default radial/angular grids. The `same_grid_*.py` scripts in `tests/benchmarks/pyscf/` perform the apples-to-apples comparison and recover sub-µHa agreement. |
| DFT (cc-pVDZ basis) | up to ~5e-6 Ha | PySCF reference uses `conv_tol=1e-10` (looser than `1e-12` for smaller bases) because cc-pVDZ has a larger condition number on these systems. Documented in `tests/golden/dft/validation.json` `settings.convergence`. |
| Raman activities (H₂O OH-stretch) | semi-analytical agreement to ~2.5e-7 | After 60+ purely analytical attempts on H₂O failed the 1e-8 bar, the production path uses semi-analytical CPHF + FD polarizability derivatives. Matches Gaussian 09 to 2.5e-7. See `scripts/phase5/generate_raman_golden.py` and the `_provenance.fd_step_bohr` field. |
| B3LYP convention | exact match to `b3lyp5` | IQCP uses **VWN5 correlation** (matching PySCF's `b3lyp5`, not the default `b3lyp` which uses VWN3). All DFT goldens carry `B3LYP_variant: "b3lyp5"` in metadata. |
| 6-31G* d-functions | exact match to Cartesian (6d) | IQCP uses 6 Cartesian d-functions for 6-31G*; PySCF reference data was generated with `mol.cart = True` to match. |

Additional regime / edge-case tolerances live alongside the relevant golden
in a sibling `tolerance.md` (per-category) where applicable.

---

## 8. Provenance Quick Reference

Every golden carries a `metadata` or `_provenance` block recording at minimum:

- **`generated`**: ISO-8601 timestamp
- **`pyscf_version`** (and SciPy / NumPy where relevant)
- **`generator`**: relative path to the script that produced the file
- **Method-specific knobs**: `conv_tol`, `cart`, `B3LYP_variant`,
  `gauge_origin`, `fd_step_bohr`, `temperature_k`, etc.

For files lacking explicit version metadata (a small number of pre-Phase-5
goldens), the upstream PySCF version is `2.11.0` and the SciPy version is
`1.17.0` per the project-wide reference stack documented in Section 1.

---

## 9. Where to file reproduction issues

If you encounter a residual outside the tolerances above:

1. Check the per-file `metadata` / `_provenance` block — most discrepancies
   trace to a version mismatch (PySCF, NumPy, basis convention).
2. Re-run the relevant generator script and `git diff` against the
   checked-in golden. Any non-zero diff that survives a clean
   `pip install`-pinned environment is a real regression.
3. Open an issue at https://github.com/ExaPsi/IQCP with the diff and your
   `pip freeze`.
