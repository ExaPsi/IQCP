# Changelog

All notable changes to IQCP (Interactive Quantum Chemistry Playground) will be
documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Releases are archived on Zenodo. The concept DOI
[10.5281/zenodo.18798309](https://doi.org/10.5281/zenodo.18798309) always
resolves to the latest version; individual versions have their own DOIs
(e.g. v2.0.0 → [10.5281/zenodo.19996612](https://doi.org/10.5281/zenodo.19996612)).

## [Unreleased]

## [2.0.0] - 2026-05-03

JCIM Full Article release. IQCP graduates from an educational SCF playground
into a complete browser-native quantum chemistry engine validated against
PySCF 2.11.0 across 108 single-point benchmark systems with sub-microhartree
accuracy for HF, LDA, and B3LYP.

### Added

- **Density Functional Theory pipeline:** closed-shell LDA (S-VWN5), B3LYP,
  and B3LYP-D3(BJ) with Becke-Johnson damping for elements H-Ar.
- **Six basis sets:** STO-3G, 3-21G, 6-31G, 6-31G*, 6-31+G*, cc-pVDZ
  (Cartesian d-functions by default).
- **DFT integration grid:** SG-1 pruning, 75 Mura-Knowles radial points
  combined with Lebedev angular grids (max 194 points).
- **Analytical energy gradients** for RHF and DFT, with
  finite-difference self-consistency cross-checks against PySCF.
- **Geometry optimization** with redundant internal coordinates and an
  L-BFGS optimizer hardened with consecutive-rejection recovery.
- **Potential-energy-surface (PES) scans** in rigid and relaxed modes.
- **Analytical Hessians** via Coupled-Perturbed Hartree-Fock (CPHF).
- **Vibrational frequency analysis** with normal-mode visualization.
- **IR intensities** from analytical dipole derivatives.
- **Semi-analytical Raman activities** from polarizability derivatives
  (validated against Gaussian 09 to 2.5e-7 Ha/bohr^2).
- **RRHO thermochemistry:** zero-point energy, enthalpy, entropy, and
  Gibbs free energy at user-selected temperature and pressure.
- **3D molecular viewer** (Three.js) with electron density and molecular
  orbital isosurface rendering.
- **Lazy-loaded spectroscopy WASM module** (`qc-wasm-spectra`, 289 KB
  gzipped) keeping the eager core (`qc-wasm`, 470 KB gzipped) lean.
- **Mulliken and Lowdin population analysis** with golden cross-validation
  against PySCF.
- **108-system PySCF validation suite** (six molecules x six basis sets x
  three methods) with all reference data committed in `tests/golden/`.
- **Performance benchmarks** (Criterion native + WASM browser timings) in
  `tests/benchmarks/`.

### Changed

- **Module renumbering** for pedagogical flow: A=Basis Set Explorer,
  B=Integral Inspector, C=Boys Function Lab, D=Rys Quadrature Lab,
  E=SCF Sandbox.
- **SCF tightened** to 1e-10 Ha energy convergence with DIIS subspace
  size of 6.
- **B3LYP variant standardized** to B3LYP5 (VWN5 correlation) to match
  PySCF defaults exactly.
- **ERI storage** consolidated to 8-fold symmetry, computed once before
  SCF iteration begins.
- **Build pipeline split** to two WASM modules to enable lazy loading of
  spectroscopy code paths on demand.

### Validated

- RHF energy: max deviation 77 nHa, median 0.4 nHa (vs PySCF 2.11.0).
- DFT energy: max deviation 148 uHa, median 21 uHa (vs PySCF 2.11.0).
- RHF gradients: 1e-5 Ha/bohr agreement (vs PySCF analytical).
- H2O frequencies: <0.01 cm^-1 deviation (vs PySCF).
- H2O ZPE: 7.0e-7 Ha agreement.
- H2O Gibbs energy at 298 K: 7.6e-7 Ha agreement.
- Cross-browser bit-identical IEEE 754 Float64 reproducibility verified
  across V8, SpiderMonkey, and JavaScriptCore.

### Test Coverage

- 1,418 Rust unit and integration tests.
- 361 TypeScript component and integration tests.

## [1.0.0] - 2026-01-24

Initial public release.

### Added

- **Module A - Boys Function Lab:** interactive F_m(T) exploration with
  series, recurrence, and asymptotic regime visualization.
- **Module B - Rys Quadrature Lab:** roots, weights, and polynomial
  reconstruction with order-error trade-off plots.
- **Module C - SCF Sandbox (RHF):** closed-shell Hartree-Fock with DIIS
  acceleration toggle, iteration trace, and Fock/density/MO inspection.
- **Three basis sets:** STO-3G, 3-21G, 6-31G.
- **WebAssembly compute core** (`qc-wasm`) with Web Worker offloading
  to keep the UI responsive (slider updates <200 ms).
- **URL-encodable deep links** for every module state, enabling
  reproducible classroom worksheets and sharing.
- **Lab Pack #1** materials in `content/labpack1/`.
- **Golden test suite** for Boys, Rys, and SCF against SciPy and PySCF
  reference values.
- **MIT license**, Zenodo DOI minted, citation metadata
  (`CITATION.cff`, `.zenodo.json`).

[Unreleased]: https://github.com/ExaPsi/IQCP/compare/v2.0.0...HEAD
[2.0.0]: https://github.com/ExaPsi/IQCP/releases/tag/v2.0.0
[1.0.0]: https://github.com/ExaPsi/IQCP/releases/tag/v1.0.0
