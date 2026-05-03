# Performance Benchmarks

This document is the reproduction guide for every performance number quoted in
the IQCP JCIM Full Article. Three execution environments are measured side
by side: native Rust (Criterion.rs), PySCF reference (`time.perf_counter()`),
and browser WASM (`performance.now()`).

---

## 1. Performance Methodology

### What is measured

| Layer            | Tool                | Harness path                              | What it captures                                  |
|------------------|---------------------|-------------------------------------------|---------------------------------------------------|
| Native Rust      | Criterion.rs        | `crates/qc-core/benches/`                 | Pure-CPU algorithm time, no I/O                   |
| Browser WASM     | `performance.now()` | `tests/benchmarks/wasm/`                  | End-to-end WASM compute time (no Worker overhead) |
| PySCF reference  | `time.perf_counter()` | `tests/benchmarks/pyscf/`               | Native-Python single-threaded SCF / gradient time |
| Bundle size      | `gzip` + `du`       | `tests/benchmarks/bundle/`                | Raw + gzipped delivered bytes                     |

### Common methodology

- **Single-threaded everywhere.** PySCF runs are pinned with
  `pyscf.lib.num_threads(1)`; the WASM harness runs on the main thread (no
  Web Worker, no SharedArrayBuffer); Criterion benches are inherently
  single-threaded per iteration.
- **Same molecules, same basis sets, same convergence.** All three layers
  use the geometries cached in `content/presets/systems.json`, the basis
  sets `STO-3G` and `6-31G*` (Cartesian d-functions, `mol.cart = True`),
  and `conv_tol = 1e-12` for the SCF.
- **Warmup before timing.** Each WASM case runs N=2 warmup iterations that
  are discarded (JIT + cache warmup), then N=10 timed iterations. Criterion
  uses its own statistical sampling; PySCF runs N=10 timed iterations.
- **Wall-clock only.** No CPU-time, no perf counters. Mean, sample standard
  deviation (`ddof=1`), min, and max are reported.
- **Geometries in bohr** to avoid implicit Å↔bohr conversion noise.

### Hardware reference

JCIM Performance numbers were collected on a single workstation:

> Intel x86_64, 16 cores, 32 GB RAM, Linux 6.x, Rust 1.94 stable,
> Node 20 LTS, Chrome 130 (WASM measurements).

Reproductions on similar hardware (any modern x86_64 / aarch64 with ≥16 GB
RAM and a current browser) reproduce the same *ordering* and *order of
magnitude* of timings. Absolute milliseconds vary with CPU clock, cache,
and browser version.

---

## 2. Running the Benchmarks

### Native Rust (Criterion)

```bash
# All benchmark groups
cargo bench --workspace

# Individual groups
cargo bench --package qc-core -- boys
cargo bench --package qc-core -- eri
cargo bench --package qc-core -- scf
cargo bench --package qc-core -- gradient
cargo bench --package qc-core -- optimizer

# Compile-check without running (CI smoke)
cargo bench --package qc-core -- --test
```

Bench targets (one Rust file per target):

| File                                 | What it measures                                  |
|--------------------------------------|---------------------------------------------------|
| `crates/qc-core/benches/boys_bench.rs`      | F_m(T) across the three regimes + sweep   |
| `crates/qc-core/benches/eri_bench.rs`       | ERI for H₂, H₂O, CH₄ (STO-3G, 6-31G*)     |
| `crates/qc-core/benches/scf_bench.rs`       | RHF + B3LYP SCF for H₂, H₂O, NH₃, CH₄     |
| `crates/qc-core/benches/gradient_bench.rs`  | RHF + DFT analytical gradients on H₂O     |
| `crates/qc-core/benches/optimizer_bench.rs` | Full geometry optimization trajectories   |

Criterion writes raw timing data + HTML reports to `target/criterion/`.
For a manuscript-ready JSON summary:

```bash
bash tests/benchmarks/criterion/run_and_collect.sh
# → tests/benchmarks/criterion/results_native.json
```

Pre-collected reference outputs live in
`tests/benchmarks/criterion/results_native.json`,
`results_scf_all.json`, `results_spherical.json`, and
`iqcp_iteration_counts.json`.

### Browser WASM

```bash
# 1. Build the WASM module first
wasm-pack build crates/qc-wasm --release --target web \
    --out-dir ../../apps/web/src/wasm

# 2. Serve the project root over HTTP
cd /path/to/IQCP
python -m http.server 8080

# 3. Open the benchmark page
# http://localhost:8080/tests/benchmarks/wasm/benchmark.html
```

The page exposes the SCF (`scf_h2_sto3g_rhf` … `scf_ch4_631gs_b3lyp`) and
geometry-optimization (`opt_h2o_sto3g_rhf`, `opt_h2o_631gs_b3lyp`) cases.
JSON results can be downloaded directly from the page UI.

### PySCF reference

```bash
source .venv/bin/activate
bash tests/benchmarks/pyscf/run_all.sh

# Individual scripts
python tests/benchmarks/pyscf/benchmark_scf.py
python tests/benchmarks/pyscf/benchmark_gradients.py
python tests/benchmarks/pyscf/benchmark_optimization.py

# Extended gradients (3-21G, 6-31+G*)
python tests/benchmarks/pyscf/generate_extended_gradients.py
```

Results are written to `tests/benchmarks/pyscf/results/` as
`scf_timing.json`, `gradient_timing.json`, `optimization_timing.json`.

### Bundle size

```bash
# Full measurement (production build + gzip)
bash tests/benchmarks/bundle/measure_bundle.sh

# Or measure an existing build without rebuilding
bash tests/benchmarks/bundle/measure_bundle.sh --skip-build
```

Output: `tests/benchmarks/bundle/results/bundle_sizes.json` with raw +
gzipped sizes per chunk plus the budget table.

---

## 3. Cross-Code Performance Comparison

### What to expect

The JCIM Performance section reports IQCP-WASM vs PySCF-native side by
side. The headline finding: **WASM is slower than native Python+C, by a
factor that varies with system size**.

| System (basis / method)      | Order-of-magnitude PySCF (native, single-thread) | Order-of-magnitude WASM (browser) |
|------------------------------|--------------------------------------------------|-----------------------------------|
| H₂ STO-3G RHF                | ~10 ms                                           | ~10–50 ms                         |
| H₂O STO-3G RHF               | ~30 ms                                           | ~50–200 ms                        |
| H₂O 6-31G* RHF               | ~100 ms                                          | ~300 ms – 1 s                     |
| H₂O 6-31G* B3LYP             | ~300 ms                                          | ~1–3 s                            |
| CH₄ 6-31G* B3LYP             | ~1 s                                             | ~3–10 s                           |
| H₂O 6-31G* B3LYP optimization (full)  | ~5 s                                    | ~30 s – 2 min                     |

These ranges are **deliberately loose** — they cover modern desktop
hardware from a recent x86_64 CPU through to mid-range laptops. The exact
numbers in the manuscript table come from the workstation described in
Section 1.

### The value proposition

IQCP is intentionally not competing with PySCF on raw speed. The value the
manuscript argues for is:

1. **Zero install** — every measurement above is reproducible by opening
   a URL.
2. **Zero server compute** — the entire SCF + DFT runs client-side in WASM,
   making per-user cost flat.
3. **Bit-identical determinism** across browsers, which native Python +
   threaded BLAS does not provide out of the box.
4. **Sub-µHa agreement** with PySCF on the 108-system benchmark
   (6 molecules × 6 basis sets × 3 methods); see `REPRODUCIBILITY.md`.

The manuscript frames the WASM↔native gap as the cost of (1)–(3).

---

## 4. Frontend Performance

### Bundle-size budgets (gzipped)

Enforced by `tests/benchmarks/bundle/measure_bundle.sh` and the CI helper
`apps/web/scripts/check-bundle-size.sh`.

| Asset            | Budget  | Enforcement | Notes                                          |
|------------------|---------|-------------|------------------------------------------------|
| Initial JS       | 600 KB  | Soft        | Plotly.js is the dominant contributor          |
| WASM module      | 250 KB  | Hard        | `qc_wasm_bg.wasm` after `wasm-pack --release`  |
| Three.js chunk   | 300 KB  | Hard        | Lazy-loaded (`viewer3d-*.js`); never on Module B/E |
| CSS              | n/a     | Tracked     | TailwindCSS purged + component styles          |
| Web Worker       | n/a     | Tracked     | `compute.worker-*.js`                          |
| html2pdf         | n/a     | Tracked     | Lazy-loaded for PDF export only                |
| Module chunks    | n/a     | Tracked     | Route-level code split (ModuleD, ModuleE, …)   |
| KaTeX fonts      | n/a     | Tracked     | Aggregate                                      |

### Code-split assertions

The bundle measurement script also verifies:

- Three.js is **not** present in the initial JS bundle (loaded only when
  the 3D viewer route is opened).
- The WASM module is **not** inlined into JS (delivered as a separate
  `.wasm` file so the browser can stream-compile).
- The compute worker is in its own chunk (no main-thread blocking).

### Lighthouse-style metrics

Bundle size is the load-bearing front-end metric in the manuscript. Time-
to-interactive and Largest Contentful Paint are dominated by the initial
JS + WASM download, so reproducing the bundle table reproduces the
front-end performance story. We do not currently check in Lighthouse
HTML reports.

---

## 5. Result Interpretation Guide

What "good" looks like for each category, and when to investigate.

| Benchmark | Healthy range | Investigate if |
|-----------|---------------|----------------|
| Boys F_m(T) (Criterion) | ~10–100 ns per evaluation across all regimes | >1 µs per evaluation, or large jump at a regime boundary (T ≈ 12 or T ≈ 30) — likely a missing regime branch |
| ERI (Criterion) | Linear in number of significant shell quartets | Super-linear scaling vs system size — missing screening |
| SCF (Criterion / WASM) | Iteration count within ±1 of PySCF for the same convergence; per-iteration time monotonic in basis size | Iteration count diverges by >2 — DIIS regression. Per-iteration time scales worse than O(N⁴) — Fock build regression |
| Gradient (Criterion) | ~2–5× SCF time for analytical RHF gradients; ~5–10× for DFT | Gradient time dominates SCF by >20× — analytical formula regression |
| Optimizer (Criterion) | Steps to converge within ±1 of PySCF for the same `gtol` | Step count diverges by >3 — L-BFGS recovery regression (see `f40731f`) |
| WASM SCF | 2–10× PySCF native | >50× — likely a debug build, missing `--release`, or the WASM module not rebuilt |
| Bundle: Initial JS | ≤600 KB gzipped | >700 KB — Plotly variant regression or accidental Three.js inclusion |
| Bundle: WASM | ≤250 KB gzipped | >250 KB — `wasm-opt` not run, or symbols not stripped |
| Bundle: Three.js chunk | ≤300 KB gzipped | Three.js appears in initial JS — code-split regression |

### Quick triage flow

If a Criterion bench regresses:

```bash
# 1. Confirm release build
cargo bench --package qc-core -- --test         # smoke
cargo bench --package qc-core -- <group>        # full

# 2. Re-run the cross-code reference and diff
source .venv/bin/activate
python tests/benchmarks/pyscf/benchmark_scf.py
diff tests/benchmarks/pyscf/results/scf_timing.json <previous>
```

If the WASM bench regresses:

```bash
# Always rebuild WASM with --release before timing
wasm-pack build crates/qc-wasm --release --target web \
    --out-dir ../../apps/web/src/wasm
```

If the bundle exceeds budget:

```bash
bash tests/benchmarks/bundle/measure_bundle.sh
# Inspect the per-chunk breakdown in
# tests/benchmarks/bundle/results/bundle_sizes.json
```
