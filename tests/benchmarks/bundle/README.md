# Bundle Size Measurement

Reproducible bundle size measurement infrastructure for the IQCP JCIM manuscript's Performance Benchmarks section.

## Purpose

Generates precise, machine-readable bundle size data for:

- **JCIM Application Note (M2):** Performance Benchmarks section reporting total JS, WASM, Three.js chunk, and CSS sizes
- **Budget enforcement:** Verifies gzipped sizes stay within defined limits
- **Regression tracking:** JSON output enables historical comparison across commits

## Quick Start

```bash
# Full measurement (builds first, then measures)
bash tests/benchmarks/bundle/measure_bundle.sh

# Measure existing build (skip npm run build)
bash tests/benchmarks/bundle/measure_bundle.sh --skip-build
```

## What It Measures

All sizes reported in both raw bytes and gzipped bytes:

| Asset | Description |
|-------|-------------|
| **Initial JS** | Main application bundle loaded on every page |
| **WASM module** | `qc_wasm_bg.wasm` -- Rust compute core |
| **Three.js chunk** | `viewer3d-*.js` -- lazy-loaded 3D rendering |
| **CSS** | Compiled TailwindCSS + component styles |
| **Web Worker** | `compute.worker-*.js` -- background compute |
| **html2pdf** | Lazy-loaded PDF export |
| **Module chunks** | Route-level code-split chunks (ModuleD, ModuleE, etc.) |
| **Fonts** | KaTeX math rendering fonts (aggregate) |

## Budget Limits (gzipped)

| Asset | Budget | Enforcement |
|-------|--------|-------------|
| Initial JS | 600 KB | Soft (Plotly.js known contributor) |
| WASM module | 250 KB | Hard |
| Three.js chunk | 300 KB | Hard |

See also: `apps/web/scripts/check-bundle-size.sh` for CI-oriented budget checks.

## Output

### Human-readable table

Printed to stdout with color-coded budget pass/fail status.

### JSON file

Written to `tests/benchmarks/bundle/results/bundle_sizes.json` with:

- **metadata:** date, toolchain versions (Node, Vite, Rust, wasm-pack)
- **sizes:** per-category raw and gzipped byte counts
- **budgets:** defined limits for comparison
- **chunks:** detailed per-file breakdown with category tags

## Dependencies

- `gzip` (system)
- `node`, `npm` (for production build)
- `wasm-pack`, `rustc` (for WASM build, if module not pre-built)

## Directory Structure

```
tests/benchmarks/bundle/
  measure_bundle.sh      # Measurement script
  README.md              # This file
  results/
    .gitkeep
    bundle_sizes.json    # Generated measurement output
```
