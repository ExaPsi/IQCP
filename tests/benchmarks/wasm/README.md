# WASM Timing Benchmarks

Wall-clock timing measurements for IQCP WASM calculations, used in the JCIM
manuscript (M2) to compare against native Rust and PySCF performance.

## Purpose

This benchmark harness measures the end-to-end computation time of quantum
chemistry calculations running through the WASM module. Results feed directly
into Table X of the JCIM manuscript where we compare three execution
environments:

| Environment | Tool | Location |
|-------------|------|----------|
| Native Rust | `cargo bench` (Criterion.rs) | `crates/qc-core/benches/` |
| PySCF | `benchmark_scf.py` | `tests/benchmarks/pyscf/` |
| **WASM** | **`benchmark.html`** | **`tests/benchmarks/wasm/`** |

## Benchmark Cases

### SCF Calculations

| Case | Molecule | Method | Basis |
|------|----------|--------|-------|
| `scf_h2_sto3g_rhf` | H2 | RHF | STO-3G |
| `scf_h2o_sto3g_rhf` | H2O | RHF | STO-3G |
| `scf_h2o_631gs_rhf` | H2O | RHF | 6-31G* |
| `scf_h2o_631gs_b3lyp` | H2O | B3LYP | 6-31G* |
| `scf_ch4_631gs_b3lyp` | CH4 | B3LYP | 6-31G* |

### Geometry Optimizations (includes gradient computation)

| Case | Molecule | Method | Basis |
|------|----------|--------|-------|
| `opt_h2o_sto3g_rhf` | H2O | RHF | STO-3G |
| `opt_h2o_631gs_b3lyp` | H2O | B3LYP | 6-31G* |

All geometries are specified in bohr and match the PySCF benchmarks exactly.

## How to Run

### Option 1: Serve from project root

Serve the project root directory with any static file server and open the
benchmark page in your browser:

```bash
cd /path/to/IQCP
python -m http.server 8080

# Open in browser:
# http://localhost:8080/tests/benchmarks/wasm/benchmark.html
```

### Option 2: Via Vite dev server

Start the Vite dev server and navigate to the benchmark page. Note that
Vite must be configured to serve files outside `apps/web/`:

```bash
cd apps/web
npm run dev

# Open in browser:
# http://localhost:5173/../../tests/benchmarks/wasm/benchmark.html
```

**Note:** For the simplest experience, use Option 1. The HTML page loads the
WASM module directly from `apps/web/src/wasm/` using relative paths, bypassing
the React app entirely.

### Option 3: Import into React app

The `benchmark.ts` module exports types and functions that can be imported:

```typescript
import { runBenchmarkSuite, BENCHMARK_CASES } from '../../../tests/benchmarks/wasm/benchmark';
```

## Prerequisites

1. **WASM module must be built** before running benchmarks:

```bash
wasm-pack build crates/qc-wasm --release --target web --out-dir ../../apps/web/src/wasm
```

2. **Browser requirements:**
   - Chrome 90+, Firefox 89+, Edge 90+, or Safari 15+
   - WebAssembly support (all modern browsers)
   - `performance.now()` high-resolution timer support

## Methodology

1. **WASM Initialization:** The WASM module is loaded and compiled once before
   any benchmarks run. Compilation time is not included in measurements.

2. **Warmup Phase:** Each benchmark case runs `N_WARMUP` (default: 2) iterations
   that are discarded. This ensures JIT compilation of JS glue code and cache
   warming.

3. **Timing Phase:** Each case then runs `N_RUNS` (default: 10) timed iterations.
   Each iteration is timed independently using `performance.now()`, which
   provides sub-millisecond precision in modern browsers.

4. **Statistics:** Mean, sample standard deviation (ddof=1), min, and max are
   computed from the timed runs.

5. **Synchronous Execution:** All WASM functions are called synchronously from
   the main thread (no Web Worker). This measures pure WASM computation time
   without message-passing overhead.

6. **UI Yielding:** Between runs, `setTimeout(0)` is used to keep the browser
   responsive. This adds negligible overhead (not included in timings).

## Output Format

Results are saved as JSON with the following structure:

```json
{
  "metadata": {
    "browser": "Chrome 130.0.6723.91",
    "platform": "Linux x86_64",
    "date": "2026-03-25",
    "user_agent": "Mozilla/5.0 ...",
    "wasm_module_size_bytes": 1153434,
    "wasm_version": "0.1.0",
    "n_warmup": 2,
    "n_runs": 10
  },
  "benchmarks": [
    {
      "name": "scf_h2o_631gs_b3lyp",
      "description": "H2O B3LYP/6-31G* full SCF",
      "category": "scf",
      "n_warmup": 2,
      "n_runs": 10,
      "mean_ms": 234.56,
      "std_ms": 12.34,
      "min_ms": 220.00,
      "max_ms": 260.00,
      "times_ms": [220.00, 225.12, ...],
      "valid": true,
      "validation_info": "E = -76.3891234567 Ha, converged = true"
    }
  ]
}
```

## Comparing with Other Benchmarks

### PySCF Comparison

PySCF benchmarks use single-threaded execution (`lib.num_threads(1)`) to match
the single-threaded WASM environment. Both use identical geometries (in bohr),
basis sets, and convergence thresholds (`conv_tol = 1e-12`).

Run PySCF benchmarks:

```bash
source .venv/bin/activate
python tests/benchmarks/pyscf/benchmark_scf.py
```

### Native Rust Comparison

Native benchmarks use Criterion.rs with the same molecular systems and
convergence parameters. The speedup ratio (Native / WASM) characterizes
the WASM overhead.

### Interpretation Notes

- WASM times include serde serialization overhead (JS -> Rust -> JS)
- WASM runs single-threaded (no SIMD, no threads unless parallel feature enabled)
- Browser JIT optimization may cause the first few timed runs to be slower
  than subsequent ones (hence the warmup phase)
- `performance.now()` resolution may be reduced in cross-origin contexts;
  run from the same origin as the WASM module for best precision

## File Structure

```
tests/benchmarks/wasm/
  README.md             # This file
  benchmark.html        # Standalone HTML benchmark page
  benchmark.ts          # TypeScript benchmark logic (importable)
  results/              # Saved JSON results
    .gitkeep
```

## Results Directory

Benchmark results are saved in `results/` with filenames like
`wasm_benchmark_2026-03-25_chrome-130.json`. These JSON files are committed
to the repository for reproducibility and manuscript reference.
