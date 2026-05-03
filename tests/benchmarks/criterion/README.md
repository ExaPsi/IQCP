# Criterion.rs Benchmark Infrastructure

Performance benchmarks for qc-core quantum chemistry algorithms, used to generate
timing data for the JCIM Application Note (M2 manuscript).

## Running Benchmarks

### All benchmarks

```bash
cargo bench --package qc-core
```

### Individual benchmark groups

```bash
cargo bench --package qc-core -- boys       # Boys function F_m(T)
cargo bench --package qc-core -- eri         # Electron repulsion integrals
cargo bench --package qc-core -- scf         # SCF convergence
cargo bench --package qc-core -- gradient    # Energy gradients
cargo bench --package qc-core -- optimizer   # Geometry optimization
```

### Compile-check only (no execution)

```bash
cargo bench --package qc-core -- --test
```

## Benchmark Targets

| File | Benchmarks | Description |
|------|-----------|-------------|
| `boys_bench.rs` | 5 | Boys function across regimes + sweep |
| `eri_bench.rs` | 4 | ERI for H2, H2O, CH4 with STO-3G and 6-31G* |
| `scf_bench.rs` | 5 | RHF and B3LYP SCF for various systems |
| `gradient_bench.rs` | 2 | RHF and DFT gradients for H2O |
| `optimizer_bench.rs` | 2 | Full geometry optimization trajectories |

## Results

Criterion stores raw results in `target/criterion/` with HTML reports.
The `run_and_collect.sh` script extracts timing data and system information
into `results_native.json` for manuscript inclusion.

```bash
bash tests/benchmarks/criterion/run_and_collect.sh
```

## Output Format

`results_native.json` contains:
- System information (CPU, OS, Rust version)
- Per-benchmark timing statistics (mean, std_dev, median)
- Benchmark metadata (date, criterion version)

## Notes

- Expensive benchmarks (DFT gradient, optimizer) use reduced sample sizes (10)
  to keep total runtime reasonable.
- Results are deterministic for a given hardware/compiler combination.
- Always run benchmarks on the same machine with minimal background load
  for comparable results.
