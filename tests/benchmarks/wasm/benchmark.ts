/**
 * WASM Timing Benchmark Harness
 *
 * Measures wall-clock time for IQCP WASM calculations to compare against
 * native Rust (Criterion.rs) and PySCF benchmarks for the JCIM manuscript.
 *
 * This module defines benchmark cases, timing utilities, and the main
 * benchmark runner. It can be used standalone in benchmark.html or
 * imported into the React app.
 *
 * Methodology:
 * - Each benchmark case runs N_WARMUP warmup iterations (discarded)
 * - Then N_RUNS timed iterations using performance.now()
 * - Statistics computed: mean, std, min, max
 * - All times in milliseconds
 *
 * @module tests/benchmarks/wasm/benchmark
 */

// ============================================================================
// Types
// ============================================================================

/**
 * A single benchmark case definition.
 */
export interface BenchmarkCase {
  /** Short identifier (e.g., "scf_h2o_sto3g_rhf") */
  name: string;
  /** Human-readable description */
  description: string;
  /** Category for grouping in output */
  category: 'scf' | 'gradient' | 'optimization';
  /** The function to benchmark. Returns a result (used to verify correctness). */
  run: (wasm: WasmModule) => unknown;
}

/**
 * Timing result for a single benchmark case.
 */
export interface BenchmarkResult {
  /** Benchmark case name */
  name: string;
  /** Human-readable description */
  description: string;
  /** Category */
  category: string;
  /** Number of warmup runs (discarded) */
  n_warmup: number;
  /** Number of timed runs */
  n_runs: number;
  /** Mean wall-clock time in milliseconds */
  mean_ms: number;
  /** Standard deviation in milliseconds */
  std_ms: number;
  /** Minimum time in milliseconds */
  min_ms: number;
  /** Maximum time in milliseconds */
  max_ms: number;
  /** All individual timed run durations in milliseconds */
  times_ms: number[];
  /** Whether the computation produced a valid result */
  valid: boolean;
  /** Optional validation info (e.g., energy value) */
  validation_info?: string;
}

/**
 * Complete benchmark suite output.
 */
export interface BenchmarkSuiteResult {
  metadata: {
    browser: string;
    platform: string;
    date: string;
    user_agent: string;
    wasm_module_size_bytes: number | null;
    wasm_version: string;
    n_warmup: number;
    n_runs: number;
  };
  benchmarks: BenchmarkResult[];
}

/**
 * Progress callback for UI updates during benchmark execution.
 */
export type ProgressCallback = (
  caseIndex: number,
  totalCases: number,
  caseName: string,
  phase: 'warmup' | 'timing',
  runIndex: number,
  totalRuns: number
) => void;

/**
 * Minimal interface for the WASM module functions used in benchmarks.
 *
 * The benchmark loads the WASM module and passes it to each case's run function.
 */
export interface WasmModule {
  /** ks_scf: unified entry point for RHF, LDA, B3LYP */
  ks_scf: (input: unknown, progressCallback?: ((p: unknown) => void) | null) => unknown;
  /** optimize_geometry: L-BFGS geometry optimization */
  optimize_geometry: (input: unknown, progressCallback: (p: unknown) => void) => unknown;
  /** version string */
  version: () => string;
}

// ============================================================================
// Configuration
// ============================================================================

/** Number of warmup iterations (discarded from timing) */
export const N_WARMUP = 2;

/** Number of timed iterations */
export const N_RUNS = 10;

// ============================================================================
// Molecule Geometries (coordinates in bohr, matching PySCF benchmarks)
// ============================================================================

/**
 * Atom format for ks_scf: [Z, x, y, z] in bohr.
 */
type AtomSpec = [number, number, number, number];

const MOLECULES: Record<string, { atoms: AtomSpec[]; label: string }> = {
  H2: {
    label: 'H2',
    atoms: [
      [1, 0, 0, 0],
      [1, 0, 0, 1.4],
    ],
  },
  H2O: {
    label: 'H2O',
    atoms: [
      [8, 0.0, 0.0, 0.2217282],
      [1, 0.0, 1.4305447, -0.8869128],
      [1, 0.0, -1.4305447, -0.8869128],
    ],
  },
  CH4: {
    label: 'CH4',
    atoms: [
      [6, 0, 0, 0],
      [1, 1.1851, 1.1851, 1.1851],
      [1, -1.1851, -1.1851, 1.1851],
      [1, -1.1851, 1.1851, -1.1851],
      [1, 1.1851, -1.1851, -1.1851],
    ],
  },
};

// ============================================================================
// Benchmark Case Definitions
// ============================================================================

/**
 * Build a ks_scf input object.
 */
function ksInput(
  atoms: AtomSpec[],
  basisName: string,
  method: string,
  opts?: { useSpherical?: boolean; gridQuality?: string }
): Record<string, unknown> {
  return {
    atoms,
    basisName,
    method,
    convergenceProfile: 'tight',
    maxIterations: 100,
    useDiis: true,
    gridQuality: opts?.gridQuality ?? 'standard',
    useSpherical: opts?.useSpherical ?? false,
  };
}

/**
 * Build an optimize_geometry input object.
 */
function optInput(
  atoms: AtomSpec[],
  basisName: string,
  method: string
): Record<string, unknown> {
  return {
    atoms,
    basisName,
    method,
    maxSteps: 50,
    gradThreshold: 4.5e-4,
    energyThreshold: 1.0e-6,
    memorySize: 7,
  };
}

/**
 * Extract energy from a ks_scf result for validation.
 */
function extractEnergy(result: unknown): number | null {
  if (result && typeof result === 'object' && 'energy' in result) {
    return (result as { energy: number }).energy;
  }
  return null;
}

/**
 * Extract convergence status from a result.
 */
function isConverged(result: unknown): boolean {
  if (result && typeof result === 'object' && 'converged' in result) {
    return (result as { converged: boolean }).converged;
  }
  return false;
}

/**
 * All benchmark cases.
 *
 * These match the cases in:
 * - tests/benchmarks/pyscf/benchmark_scf.py (PySCF comparison)
 * - Native Criterion benchmarks (when available)
 */
export const BENCHMARK_CASES: BenchmarkCase[] = [
  // =========================================================================
  // SCF benchmarks
  // =========================================================================
  {
    name: 'scf_h2_sto3g_rhf',
    description: 'H2 RHF/STO-3G full SCF',
    category: 'scf',
    run: (wasm) => wasm.ks_scf(ksInput(MOLECULES.H2.atoms, 'sto-3g', 'rhf')),
  },
  {
    name: 'scf_h2o_sto3g_rhf',
    description: 'H2O RHF/STO-3G full SCF',
    category: 'scf',
    run: (wasm) => wasm.ks_scf(ksInput(MOLECULES.H2O.atoms, 'sto-3g', 'rhf')),
  },
  {
    name: 'scf_h2o_631gs_rhf',
    description: 'H2O RHF/6-31G* full SCF',
    category: 'scf',
    run: (wasm) => wasm.ks_scf(ksInput(MOLECULES.H2O.atoms, '6-31g*', 'rhf')),
  },
  {
    name: 'scf_h2o_631gs_b3lyp',
    description: 'H2O B3LYP/6-31G* full SCF',
    category: 'scf',
    run: (wasm) =>
      wasm.ks_scf(ksInput(MOLECULES.H2O.atoms, '6-31g*', 'b3lyp')),
  },
  {
    name: 'scf_ch4_631gs_b3lyp',
    description: 'CH4 B3LYP/6-31G* full SCF',
    category: 'scf',
    run: (wasm) =>
      wasm.ks_scf(ksInput(MOLECULES.CH4.atoms, '6-31g*', 'b3lyp')),
  },

  // =========================================================================
  // Optimization benchmarks (includes gradient computation internally)
  // =========================================================================
  {
    name: 'opt_h2o_sto3g_rhf',
    description: 'H2O RHF/STO-3G geometry optimization',
    category: 'optimization',
    run: (wasm) => {
      const noop = () => {};
      return wasm.optimize_geometry(
        optInput(MOLECULES.H2O.atoms, 'sto-3g', 'rhf'),
        noop
      );
    },
  },
  {
    name: 'opt_h2o_631gs_b3lyp',
    description: 'H2O B3LYP/6-31G* geometry optimization',
    category: 'optimization',
    run: (wasm) => {
      const noop = () => {};
      return wasm.optimize_geometry(
        optInput(MOLECULES.H2O.atoms, '6-31g*', 'b3lyp'),
        noop
      );
    },
  },
];

// ============================================================================
// Statistics Helpers
// ============================================================================

/**
 * Compute mean of a numeric array.
 */
export function mean(values: number[]): number {
  if (values.length === 0) return 0;
  return values.reduce((sum, v) => sum + v, 0) / values.length;
}

/**
 * Compute sample standard deviation (ddof=1).
 */
export function std(values: number[]): number {
  if (values.length <= 1) return 0;
  const m = mean(values);
  const squaredDiffs = values.map((v) => (v - m) ** 2);
  return Math.sqrt(
    squaredDiffs.reduce((sum, v) => sum + v, 0) / (values.length - 1)
  );
}

/**
 * Compute min of a numeric array.
 */
export function min(values: number[]): number {
  return Math.min(...values);
}

/**
 * Compute max of a numeric array.
 */
export function max(values: number[]): number {
  return Math.max(...values);
}

// ============================================================================
// Benchmark Runner
// ============================================================================

/**
 * Run a single benchmark case with warmup and timing.
 *
 * @param benchCase - The benchmark case to run
 * @param wasm - The initialized WASM module
 * @param nWarmup - Number of warmup iterations
 * @param nRuns - Number of timed iterations
 * @param onProgress - Optional progress callback
 * @param caseIndex - Index of this case in the suite (for progress)
 * @param totalCases - Total number of cases (for progress)
 * @returns Timing result for this benchmark case
 */
export function runBenchmarkCase(
  benchCase: BenchmarkCase,
  wasm: WasmModule,
  nWarmup: number = N_WARMUP,
  nRuns: number = N_RUNS,
  onProgress?: ProgressCallback,
  caseIndex: number = 0,
  totalCases: number = 1
): BenchmarkResult {
  // Warmup phase
  let lastResult: unknown = null;
  for (let i = 0; i < nWarmup; i++) {
    onProgress?.(caseIndex, totalCases, benchCase.name, 'warmup', i, nWarmup);
    try {
      lastResult = benchCase.run(wasm);
    } catch (e) {
      console.error(`[Benchmark] Warmup failed for ${benchCase.name}:`, e);
      return {
        name: benchCase.name,
        description: benchCase.description,
        category: benchCase.category,
        n_warmup: nWarmup,
        n_runs: nRuns,
        mean_ms: 0,
        std_ms: 0,
        min_ms: 0,
        max_ms: 0,
        times_ms: [],
        valid: false,
        validation_info: `Error during warmup: ${e instanceof Error ? e.message : String(e)}`,
      };
    }
  }

  // Timing phase
  const times: number[] = [];
  for (let i = 0; i < nRuns; i++) {
    onProgress?.(caseIndex, totalCases, benchCase.name, 'timing', i, nRuns);
    const t0 = performance.now();
    try {
      lastResult = benchCase.run(wasm);
    } catch (e) {
      console.error(`[Benchmark] Run ${i} failed for ${benchCase.name}:`, e);
      return {
        name: benchCase.name,
        description: benchCase.description,
        category: benchCase.category,
        n_warmup: nWarmup,
        n_runs: nRuns,
        mean_ms: 0,
        std_ms: 0,
        min_ms: 0,
        max_ms: 0,
        times_ms: times,
        valid: false,
        validation_info: `Error on run ${i}: ${e instanceof Error ? e.message : String(e)}`,
      };
    }
    const t1 = performance.now();
    times.push(t1 - t0);
  }

  // Validation
  const energy = extractEnergy(lastResult);
  const converged = isConverged(lastResult);
  const validationInfo =
    energy !== null
      ? `E = ${energy.toFixed(10)} Ha, converged = ${converged}`
      : 'no energy in result';

  return {
    name: benchCase.name,
    description: benchCase.description,
    category: benchCase.category,
    n_warmup: nWarmup,
    n_runs: nRuns,
    mean_ms: mean(times),
    std_ms: std(times),
    min_ms: min(times),
    max_ms: max(times),
    times_ms: times,
    valid: true,
    validation_info: validationInfo,
  };
}

/**
 * Run the complete benchmark suite.
 *
 * Runs all benchmark cases sequentially with warmup + timing,
 * collects metadata, and returns a complete results object.
 *
 * @param wasm - The initialized WASM module
 * @param wasmSizeBytes - Size of the WASM binary in bytes (null if unknown)
 * @param onProgress - Optional progress callback
 * @param cases - Specific cases to run (defaults to all)
 * @param nWarmup - Number of warmup iterations (defaults to N_WARMUP)
 * @param nRuns - Number of timed iterations (defaults to N_RUNS)
 * @returns Complete benchmark suite results
 */
export function runBenchmarkSuite(
  wasm: WasmModule,
  wasmSizeBytes: number | null = null,
  onProgress?: ProgressCallback,
  cases: BenchmarkCase[] = BENCHMARK_CASES,
  nWarmup: number = N_WARMUP,
  nRuns: number = N_RUNS
): BenchmarkSuiteResult {
  // Detect browser
  const ua = navigator.userAgent;
  let browser = 'Unknown';
  if (ua.includes('Firefox/')) {
    const match = ua.match(/Firefox\/([\d.]+)/);
    browser = match ? `Firefox ${match[1]}` : 'Firefox';
  } else if (ua.includes('Chrome/')) {
    const match = ua.match(/Chrome\/([\d.]+)/);
    browser = match ? `Chrome ${match[1]}` : 'Chrome';
  } else if (ua.includes('Safari/') && !ua.includes('Chrome/')) {
    const match = ua.match(/Version\/([\d.]+)/);
    browser = match ? `Safari ${match[1]}` : 'Safari';
  } else if (ua.includes('Edg/')) {
    const match = ua.match(/Edg\/([\d.]+)/);
    browser = match ? `Edge ${match[1]}` : 'Edge';
  }

  // Get WASM version
  let wasmVersion = 'unknown';
  try {
    wasmVersion = wasm.version();
  } catch {
    // version() might not be available
  }

  const metadata = {
    browser,
    platform: navigator.platform || 'Unknown',
    date: new Date().toISOString().split('T')[0],
    user_agent: ua,
    wasm_module_size_bytes: wasmSizeBytes,
    wasm_version: wasmVersion,
    n_warmup: nWarmup,
    n_runs: nRuns,
  };

  const benchmarks: BenchmarkResult[] = [];
  for (let i = 0; i < cases.length; i++) {
    const result = runBenchmarkCase(
      cases[i],
      wasm,
      nWarmup,
      nRuns,
      onProgress,
      i,
      cases.length
    );
    benchmarks.push(result);
  }

  return { metadata, benchmarks };
}

/**
 * Format a BenchmarkSuiteResult as a human-readable table string.
 */
export function formatResultsTable(results: BenchmarkSuiteResult): string {
  const lines: string[] = [];
  lines.push('IQCP WASM Benchmark Results');
  lines.push('='.repeat(80));
  lines.push(`Browser:     ${results.metadata.browser}`);
  lines.push(`Platform:    ${results.metadata.platform}`);
  lines.push(`Date:        ${results.metadata.date}`);
  lines.push(`WASM v:      ${results.metadata.wasm_version}`);
  if (results.metadata.wasm_module_size_bytes !== null) {
    const kb = (results.metadata.wasm_module_size_bytes / 1024).toFixed(1);
    lines.push(`WASM size:   ${kb} KB`);
  }
  lines.push(`Warmup/Runs: ${results.metadata.n_warmup} / ${results.metadata.n_runs}`);
  lines.push('');
  lines.push(
    'Case'.padEnd(30) +
      'Mean (ms)'.padStart(12) +
      'Std (ms)'.padStart(12) +
      'Min (ms)'.padStart(12) +
      'Max (ms)'.padStart(12) +
      '  Valid'
  );
  lines.push('-'.repeat(80));

  for (const b of results.benchmarks) {
    lines.push(
      b.name.padEnd(30) +
        b.mean_ms.toFixed(2).padStart(12) +
        b.std_ms.toFixed(2).padStart(12) +
        b.min_ms.toFixed(2).padStart(12) +
        b.max_ms.toFixed(2).padStart(12) +
        (b.valid ? '  YES' : '  NO')
    );
    if (b.validation_info) {
      lines.push('  ' + b.validation_info);
    }
  }

  lines.push('='.repeat(80));
  return lines.join('\n');
}
