# PySCF Timing Benchmarks

Timing comparison data for the JCIM M2 manuscript. Measures PySCF performance on the same molecules, basis sets, and methods used by IQCP, providing a reference for browser-based WASM vs. native Python execution.

## How to Run

```bash
# From project root
bash tests/benchmarks/pyscf/run_all.sh

# Or individual scripts
source .venv/bin/activate
python tests/benchmarks/pyscf/benchmark_scf.py
python tests/benchmarks/pyscf/benchmark_gradients.py
python tests/benchmarks/pyscf/benchmark_optimization.py
```

## Environment

- **PySCF** 2.11.0
- **Python** 3.12.x
- **NumPy** 1.26.4
- **SciPy** 1.17.0
- Virtual environment: `.venv/`

## Methodology

- **Single-threaded** execution (`lib.num_threads(1)`) for fair comparison with single-threaded WASM in the browser.
- Geometries and convergence thresholds match `tests/golden/dft/validation.json` exactly.
- Coordinates in bohr; 6-31G\* uses Cartesian d-functions (`mol.cart = True`) to match IQCP convention.
- `conv_tol = 1e-12` for all calculations.
- Wall-clock timing via `time.perf_counter()`.

### Benchmark Scripts

| Script | What it measures | Runs | Combinations |
|--------|------------------|------|-------------|
| `benchmark_scf.py` | SCF energy convergence | 10 | 5 molecules x 2 basis sets x 3 methods = 30 |
| `benchmark_gradients.py` | SCF + analytical gradients | 10 | 4 selected cases |
| `benchmark_optimization.py` | Full geometry optimization | 5 | 2 selected cases |

### Molecules

H2, H2O, NH3, HF, CH4 (geometries from IQCP presets).

### Methods

RHF, LDA (VWN5), B3LYP.

### Basis Sets

STO-3G, 6-31G\*.

## Output

Results are saved as JSON in `results/`:

- `results/scf_timing.json`
- `results/gradient_timing.json`
- `results/optimization_timing.json`

Each JSON file includes full environment metadata (versions, CPU, date, threading mode) for reproducibility.

## Output Format

```json
{
  "metadata": { ... },
  "benchmarks": [
    {
      "molecule": "H2",
      "basis": "sto-3g",
      "method": "RHF",
      "n_basis": 2,
      "n_electrons": 2,
      "energy_ha": -1.11671432506255,
      "mean_s": 0.012,
      "std_s": 0.001,
      "min_s": 0.011,
      "max_s": 0.014,
      "n_runs": 10,
      "converged": true
    }
  ]
}
```
