# Interactive Quantum Chemistry Playground (IQCP)

A browser-based educational web application for teaching quantum chemistry concepts through interactive exploration.

[![Live Demo](https://img.shields.io/badge/demo-iqcp.dev-blue.svg)](https://iqcp.dev)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Live Site:** [https://iqcp.dev](https://iqcp.dev)

## Overview

IQCP enables students to interactively explore **Boys functions**, **Rys quadrature**, and **SCF convergence** through transparent, reproducible computation - all within a web browser with no installation required.

**Publication:** This software accompanies a Technology Report submitted to the *Journal of Chemical Education*.

### Key Features

- **Zero Installation:** Runs entirely in the browser using WebAssembly
- **Interactive Exploration:** Real-time parameter adjustment with instant visual feedback
- **Transparent Computation:** Inspect intermediate values and understand algorithm internals
- **Reproducible Results:** Deep links and exportable artifacts for every calculation
- **Custom Molecules:** Enter arbitrary geometries and compute integrals on-the-fly (5 basis sets)
- **Educational Focus:** Designed for classroom use with guided lab activities

## Architecture

```
Browser SPA (React + TypeScript + Vite)
    | postMessage
Web Worker (non-blocking computation)
    | wasm-bindgen
qc-wasm (Rust WASM module)
    |
qc-core (pure Rust algorithms)
```

All heavy computation runs in a Web Worker to keep the UI responsive.

## Core Modules

### Module B: Boys Function Lab
Interactive exploration of F_m(T) with regime visualization:
- **Series expansion** (T < 12): Taylor expansion with double factorial
- **Recurrence** (12 <= T < 30): erf(sqrt(T)) + upward recurrence
- **Asymptotic** (T >= 30): Asymptotic series expansion

### Module C: Rys Quadrature Lab
Roots and weights computation for Gaussian integrals:
- Moment computation via Boys functions
- Modified Chebyshev algorithm for recurrence coefficients
- Eigenvalue decomposition for roots and weights
- Order-error tradeoff visualization

### Module E: SCF Sandbox (RHF)
Restricted Hartree-Fock self-consistent field calculations:
- Interactive convergence visualization
- DIIS acceleration toggle
- Matrix inspection (Fock, density, coefficients)
- Convergence profiles: loose (1e-4), medium (1e-6), tight (1e-8)

## Getting Started

### Try It Now

Visit [https://iqcp.dev](https://iqcp.dev) to use IQCP directly in your browser - no installation required.

### Building from Source

#### Prerequisites

- **Rust** (stable, 1.70+)
- **wasm-pack** (`cargo install wasm-pack`)
- **Node.js** (18+)
- **npm**

#### Installation

```bash
# Clone the repository
git clone https://github.com/ExaPsi/IQCP.git
cd IQCP

# Build the Rust workspace
cargo build --workspace

# Build the WASM module
wasm-pack build crates/qc-wasm --release --target web --out-dir ../../apps/web/src/wasm

# Install frontend dependencies
cd apps/web
npm install

# Start the development server
npm run dev
```

The application will be available at `http://localhost:5173`.

## Project Structure

```
IQCP/
├── apps/
│   └── web/                 # React SPA (Vite + TypeScript)
│       ├── src/
│       │   ├── components/  # React components
│       │   ├── hooks/       # Custom hooks
│       │   ├── lib/         # Utilities
│       │   ├── stores/      # Zustand state management
│       │   ├── types/       # TypeScript types
│       │   └── worker/      # Web Worker
│       └── ...
├── crates/
│   ├── qc-core/             # Pure Rust algorithms
│   │   └── src/
│   │       ├── boys/        # Boys function F_m(T)
│   │       ├── rys/         # Rys quadrature
│   │       └── scf/         # SCF engine with DIIS
│   ├── qc-wasm/             # wasm-bindgen exports
│   └── qc-io/               # Shared schemas
├── content/
│   ├── labpack1/            # Lab Pack #1 materials
│   └── presets/             # Pre-computed molecular systems
└── tests/
    └── golden/              # Golden test reference data
```

## Technology Stack

### Frontend
- React 18+ with TypeScript (strict mode)
- Vite for bundling
- TailwindCSS for styling
- Zustand for state management
- Plotly.js for interactive plots

### Compute
- Rust (stable)
- wasm-bindgen + wasm-pack
- serde + serde-wasm-bindgen

## Numerical Specifications

| Algorithm | Tolerance | Notes |
|-----------|-----------|-------|
| Boys F_m(T) | 1e-12 | Absolute, all regimes |
| Rys roots | (0, 1) | Strict open interval |
| Rys weights | > 0 | Strictly positive |
| Rys reconstruction | 1e-10 | Moment accuracy |
| SCF energy | 1e-8 Ha | Final converged value |
| SCF density | 1e-6 | Frobenius norm |

## Deep Links and Artifacts

Every module state is URL-encodable for sharing and reproducibility:

```
https://iqcp.dev/boys?s=<base64-encoded-state>
https://iqcp.dev/rys?s=<base64-encoded-state>
https://iqcp.dev/scf?s=<base64-encoded-state>
```

Exportable artifacts include computation results and metadata for assessment and reproducibility.

## Lab Materials

The `content/labpack1/` directory contains guided lab activities:
- Student worksheets with deep links to pre-configured calculations
- Instructor answer keys and grading rubrics
- Performance assessment rubrics aligned with learning objectives

## References

### Primary Literature

- **Boys function:** Shavitt, I. (1963). *Methods in Computational Physics*, 2, 1-45.
- **Rys quadrature:** Dupuis, M., Rys, J., & King, H. F. (1976). *J. Chem. Phys.*, 65, 111-116.
- **DIIS acceleration:** Pulay, P. (1980). *Chem. Phys. Lett.*, 73, 393-398; (1982). *J. Comput. Chem.*, 3, 556-560.

## Citation

If you use IQCP in your teaching or research, please cite:

```bibtex
@software{iqcp2026,
  author = {{IQCP Contributors}},
  title = {Interactive Quantum Chemistry Playground (IQCP)},
  year = {2026},
  url = {https://github.com/ExaPsi/IQCP},
  note = {Supporting Information for J. Chem. Educ.}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

IQCP is developed for educational use in quantum chemistry courses. We thank the computational chemistry community for the foundational algorithms and the educators who provided feedback on pedagogical approaches.
