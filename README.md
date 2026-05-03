# IQCP: Interactive Quantum Chemistry Playground

A client-side quantum chemistry engine that compiles a Rust implementation of Hartree–Fock and Kohn–Sham DFT to WebAssembly, executing all numerical work directly in the browser with no server.

[![Live Demo](https://img.shields.io/badge/demo-iqcp.dev-blue.svg)](https://iqcp.dev)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18798309.svg)](https://doi.org/10.5281/zenodo.18798309)

**Live Site:** [https://iqcp.dev](https://iqcp.dev)

## Overview

IQCP supports closed-shell **RHF**, **LDA**, **B3LYP**, and **B3LYP-D3(BJ)** with six basis sets (STO-3G, 3-21G, 6-31G, 6-31G\*, 6-31+G\*, cc-pVDZ) for elements H–Ar. The engine implements the full computational pipeline:

- SCF convergence with DIIS acceleration
- Analytical energy gradients and geometry optimization
- Potential-energy-surface (PES) scans with rigid and relaxed modes
- Analytical Hessians via Coupled-Perturbed Hartree–Fock (CPHF)
- Vibrational frequencies and normal-mode analysis
- IR intensities from dipole derivatives
- Semi-analytical Raman activities from polarizability derivatives
- RRHO thermochemistry (ZPE, enthalpy, entropy, Gibbs free energy)
- Electron density and molecular orbital visualization
- Mulliken and Löwdin population analysis

All computation runs entirely client-side via two WebAssembly modules (470 KB core + 289 KB spectra, gzipped).

## Validation

Validated against PySCF 2.11.0 across **108 single-point benchmark systems** (six molecules × six basis sets × three methods):

| Quantity                                     | Maximum deviation   | Median  |
| -------------------------------------------- | ------------------- | ------- |
| RHF energy                                   | 77 nHa              | 0.4 nHa |
| DFT energy                                   | 148 µHa             | 21 µHa  |
| RHF gradient                                 | 10⁻⁵ Ha/bohr        | —       |
| H₂O frequency                                | <0.01 cm⁻¹          | —       |
| H₂O ZPE                                      | 7.0 × 10⁻⁷ Ha       | —       |
| H₂O Gibbs (298 K)                            | 7.6 × 10⁻⁷ Ha       | —       |
| Raman polarizability deriv. (vs Gaussian 09) | 2.5 × 10⁻⁷ Ha/bohr² | —       |

The 1,418 Rust tests + 361 TypeScript tests in this repository cover all reported numerical claims. See **Reproducing Manuscript Results** below.

## Quick Start

### Try It Now

Visit [https://iqcp.dev](https://iqcp.dev) — no installation required.

### Local Build

**Prerequisites:** Rust 1.94+, wasm-pack, Node.js 18+, npm

```bash
git clone https://github.com/ExaPsi/IQCP.git
cd IQCP

# Build Rust workspace (4 crates)
cargo build --workspace

# Build both WASM modules
wasm-pack build crates/qc-wasm --release --target web --out-dir ../../apps/web/src/wasm
wasm-pack build crates/qc-wasm-spectra --release --target web --out-dir ../../apps/web/src/wasm-spectra

# Frontend
cd apps/web
npm install
npm run dev
```

Application available at `http://localhost:5173`.

## Architecture

```
Browser SPA (React + TypeScript + Vite)
    ↓ postMessage
Web Worker (off-main-thread compute)
    ↓ wasm-bindgen (lazy import for spectra)
qc-wasm (470 KB) + qc-wasm-spectra (289 KB)
    ↓ Rust FFI
qc-core (pure Rust algorithms, native + WASM)
```

**Crates:**

- `qc-core` — pure Rust algorithms (SCF, DFT, gradients, Hessians, IR, Raman, thermo)
- `qc-wasm` — eager-loaded WASM bindings (SCF/DFT/gradients/optimization)
- `qc-wasm-spectra` — lazy-loaded WASM bindings (Hessian, CPHF, IR, Raman, thermochemistry)
- `qc-io` — shared serde schemas (RunState, RunArtifact)

## Six Interactive Modules

1. **Module A — Basis Set Explorer:** Radial profiles, contraction coefficients, basis comparison
2. **Module B — Integral Inspector:** Overlap, kinetic, nuclear, ERI matrices with primitive breakdown
3. **Module C — Boys Function Lab:** Series and recurrence regimes, regime visualization
4. **Module D — Rys Quadrature Lab:** Roots, weights, polynomial reconstruction, order-error trade-offs
5. **Module E — SCF Sandbox:** RHF/LDA/B3LYP/B3LYP-D3(BJ) with DIIS, 3D molecular viewer, PES scans, geometry optimization, electron-density and molecular-orbital isosurface visualization, Mulliken/Löwdin population analysis
6. **Module F — Frequency & Vibrational Spectra:** Analytical Hessians via CPHF, normal-mode analysis, IR intensities, semi-analytical Raman activities, RRHO thermochemistry (ZPE, H, S, G), simulated IR/Raman spectra with Lorentzian/Gaussian broadening (lazy-loaded via `qc-wasm-spectra`)

Every module state is URL-encodable via deep links for reproducibility.

## Reproducing Manuscript Results

The numerical claims in the JCIM manuscript trace to data and scripts committed in this repository.

### Validation Scripts (PySCF reference generators)

```bash
# Activate Python environment with PySCF 2.11.0
python -m venv .venv && source .venv/bin/activate
pip install pyscf==2.11.0 numpy scipy

# Reproduce gradient FD self-consistency (manuscript §4.4)
python scripts/validation/pyscf_fd_self_consistency.py

# Reproduce geometry optimization (manuscript §4.5)
python scripts/validation/pyscf_geometry_optimization.py

# Regenerate population analysis golden file
python scripts/validation/pyscf_population_analysis.py

# Phase 5 golden generators (Hessian, IR, Raman, thermo)
python scripts/phase5/generate_phase5_golden.py
python scripts/phase5/generate_ir_golden.py
python scripts/phase5/generate_raman_golden.py
python scripts/phase5/generate_thermochem_golden.py
```

### Golden Test Data (`tests/golden/`)

Reference values used by `cargo test`:

- `boys/`, `rys/`, `eri/`, `integrals/` — integral evaluation
- `scf/`, `dft/` — SCF energies, DFT functionals, gradients, populations
- `pes/` — potential energy surface scans
- `phase5/` — Hessians, frequencies, optimized geometries
- `ir/`, `raman/` — vibrational spectroscopy intensities
- `thermo/`, `thermochem/` — thermodynamic properties

### Performance Benchmarks (`tests/benchmarks/`)

```bash
# Run native Rust benchmarks (Criterion)
cargo bench -p qc-core

# PySCF timing reference (regenerate)
python tests/benchmarks/pyscf/run_scf_timing.py
```

Browser benchmarks: see WASM benchmark JSONs in `tests/benchmarks/wasm/results/`.

### Run the Full Test Suite

```bash
# Workspace tests (1,418 Rust)
cargo test --workspace

# Frontend tests (361 TypeScript)
cd apps/web
npm test
```

## Project Structure

```
IQCP/
├── apps/
│   └── web/                 # React 18 + TypeScript + Vite
│       ├── src/
│       │   ├── components/  # Per-module UI (basis, boys, rys, scf, ...)
│       │   ├── hooks/       # useWorker, useFrequency, useOptimizer, ...
│       │   ├── stores/      # Zustand stores
│       │   ├── worker/      # Web Worker handlers
│       │   └── wasm*/       # Built WASM modules
├── crates/
│   ├── qc-core/             # Pure Rust algorithms
│   ├── qc-wasm/             # Core WASM bindings (eager)
│   ├── qc-wasm-spectra/     # Spectroscopy WASM bindings (lazy)
│   └── qc-io/               # Shared schemas
├── content/
│   └── presets/             # Curated systems and pre-computed integrals
├── scripts/
│   ├── validation/          # PySCF reference scripts (FD, opt, population)
│   ├── phase5/              # Phase 5 golden data generators
│   └── *.py                 # Other utility scripts
└── tests/
    ├── golden/              # Reference values for cargo test
    └── benchmarks/          # Criterion (Rust), PySCF, WASM browser timings
```

## Development

### Common Commands

```bash
# Rust
cargo build --workspace
cargo test --workspace
cargo fmt --all
cargo clippy --workspace

# WASM (both modules)
wasm-pack build crates/qc-wasm --release --target web --out-dir ../../apps/web/src/wasm
wasm-pack build crates/qc-wasm-spectra --release --target web --out-dir ../../apps/web/src/wasm-spectra

# Web
cd apps/web
npm run dev
npm run build
npm run lint
npm run typecheck
npm test
```

### Pre-Commit Checklist

```bash
cargo fmt --all -- --check
cargo clippy --workspace -- -D warnings
cargo test --workspace
cd apps/web && npm run lint && npm run typecheck && npm run build
```

## Numerical Settings

| Setting                       | Value                                                                             |
| ----------------------------- | --------------------------------------------------------------------------------- |
| SCF energy convergence        | 1 × 10⁻¹⁰ Ha                                                                      |
| DIIS subspace size            | 6 vectors                                                                         |
| ERI storage                   | 8-fold symmetry, stored once before SCF                                           |
| DFT grid                      | SG-1 pruning, 75 Mura–Knowles radial × Lebedev (max 194)                          |
| B3LYP variant                 | B3LYP5 (VWN5 correlation)                                                         |
| Cartesian d-functions         | Default (`mol.cart=True` analog)                                                  |
| Cross-browser reproducibility | Bit-identical IEEE 754 Float64 (verified across V8, SpiderMonkey, JavaScriptCore) |

## Performance

| Workflow                      | Time (Apple M2 Max, Chrome) |
| ----------------------------- | --------------------------- |
| WASM init + first calculation | ~50 ms                      |
| H₂O B3LYP/6-31G\* SCF         | 183 ms                      |
| C₆H₆ B3LYP/6-31G\* SCF        | ~28 s                       |
| H₂O frequency analysis        | ~5 s                        |

Native Rust outperforms single-threaded PySCF for 95 of 108 benchmark systems under end-to-end conditions; the crossover with PySCF for B3LYP falls between 102 and 126 basis functions.

## Citation

If you use IQCP in your research, please cite both the software (Zenodo) and the JCIM manuscript:

```bibtex
@software{iqcp2026software,
  title   = {IQCP: Client-Side Quantum Chemistry from SCF to Vibrational
             Spectra via WebAssembly},
  author  = {Sawatlon, Boodsarin and Paiboonvorachat, Nattapong and
             Vchirawongkwin, Viwat},
  year    = {2026},
  version = {2.0.0},
  doi     = {10.5281/zenodo.19996612},
  url     = {https://github.com/ExaPsi/IQCP}
}
```

A companion JCIM Article is in preparation; the citation will be updated upon publication.

## References

Key algorithmic references implemented in this work:

- **Boys function:** Shavitt, I. _Methods in Computational Physics_ **1963**, _2_, 1–45.
- **Rys quadrature:** Dupuis, M.; Rys, J.; King, H. F. _J. Chem. Phys._ **1976**, _65_, 111–116.
- **DIIS:** Pulay, P. _Chem. Phys. Lett._ **1980**, _73_, 393–398; _J. Comput. Chem._ **1982**, _3_, 556–560.
- **B3LYP:** Becke, A. D. _J. Chem. Phys._ **1993**, _98_, 5648–5652; Stephens et al. _J. Phys. Chem._ **1994**, _98_, 11623–11627.
- **D3(BJ) dispersion:** Grimme, S. et al. _J. Chem. Phys._ **2010**, _132_, 154104; _J. Comput. Chem._ **2011**, _32_, 1456–1465.
- **CPHF for Raman:** Amos, R. D. _Chem. Phys. Lett._ **1986**, _124_, 376–381.
- **PySCF (validation reference):** Sun, Q. et al. _J. Chem. Phys._ **2020**, _153_, 024109.
- **WebAssembly:** Haas, A. et al. _PLDI_ **2017**.

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgments

We thank the developers of PySCF and Gaussian 09 for the reference implementations used in cross-code validation, and the Rust and WebAssembly communities for the toolchains that make browser-native scientific computing possible.
