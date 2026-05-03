/**
 * Pure math helpers for SpectrumPlot.tsx.
 *
 * Extracted into its own module (no React or Plotly imports) so that the
 * helpers can be unit-tested in a Node environment without pulling in the
 * browser-only Plotly bundle.
 *
 * Algorithms match `crates/qc-core/src/spectra.rs` bit-for-bit (Lorentzian
 * HWHM = FWHM/2; Gaussian σ = FWHM/(2√(2·ln 2))).
 *
 * @module components/scf/spectrumMath
 * @see US-102 Frequency Tab UI, AC2
 */

import type { BroadeningKind } from '../../worker/protocol';

/** Grid used for the broadened trace. Matches qc-core: 0..=4500 step 1. */
export const SPECTRUM_GRID_POINTS = 4501;

/**
 * Build the default wavenumber grid used by the client-side broadening
 * recompute and by the qc-core `simulate_ir_spectrum` function. Length 4501.
 */
export function buildWavenumberGrid(): number[] {
  const grid = new Array<number>(SPECTRUM_GRID_POINTS);
  for (let i = 0; i < SPECTRUM_GRID_POINTS; i++) grid[i] = i;
  return grid;
}

/**
 * Compute a broadened spectrum from stick frequencies + amplitudes on the
 * default 0..=4500 cm⁻¹ grid.
 *
 * Imaginary frequencies (ν ≤ 0) are skipped. Kernel normalization matches
 * `crates/qc-core/src/spectra.rs`:
 *
 *   Lorentzian:  L(Δω) = (γ/π) / (Δω² + γ²)    with γ = FWHM/2
 *   Gaussian:    G(Δω) = 1/(σ√(2π)) · exp(−Δω²/(2σ²))  with σ = FWHM/(2√(2·ln 2))
 *
 * Both kernels integrate to 1 over the full real line, so the integrated
 * area of each peak in the output equals the original stick amplitude.
 *
 * @param frequenciesCm1 - Raw frequency values in cm⁻¹ (imaginary allowed, skipped)
 * @param amplitudes     - Intensity or activity values (same length as frequencies)
 * @param kind           - Broadening kernel
 * @param fwhmCm1        - FWHM in cm⁻¹
 * @returns              - { grid, intensity } with grid length 4501
 */
export function computeBroadenedSpectrum(
  frequenciesCm1: readonly number[],
  amplitudes: readonly number[],
  kind: BroadeningKind,
  fwhmCm1: number
): { grid: number[]; intensity: number[] } {
  const grid = buildWavenumberGrid();
  const intensity = new Array<number>(grid.length).fill(0);

  // Filter real modes
  const realFreqs: number[] = [];
  const realAmps: number[] = [];
  for (let i = 0; i < frequenciesCm1.length; i++) {
    if (frequenciesCm1[i] > 0) {
      realFreqs.push(frequenciesCm1[i]);
      realAmps.push(amplitudes[i] ?? 0);
    }
  }

  if (realFreqs.length === 0 || fwhmCm1 <= 0) {
    return { grid, intensity };
  }

  // Lorentzian HWHM
  const halfFwhm = fwhmCm1 / 2;

  // Gaussian σ and prefactor. Using the area-normalized form
  //   G(Δω) = 1/(σ·√(2π)) · exp(-Δω²/(2σ²))
  // which is mathematically identical to the "4 ln 2 / FWHM²" form but
  // matches the Rust reference in `lorentzian_value`/`gaussian_value`.
  const sigma = fwhmCm1 / (2 * Math.sqrt(2 * Math.log(2)));
  const gaussPrefac = 1 / (sigma * Math.sqrt(2 * Math.PI));
  const gaussExpDenom = 2 * sigma * sigma;

  for (let k = 0; k < grid.length; k++) {
    const w = grid[k];
    let sum = 0;
    for (let i = 0; i < realFreqs.length; i++) {
      const dw = w - realFreqs[i];
      if (kind === 'lorentzian') {
        sum +=
          (realAmps[i] * (halfFwhm / Math.PI)) /
          (dw * dw + halfFwhm * halfFwhm);
      } else {
        sum += realAmps[i] * gaussPrefac * Math.exp(-(dw * dw) / gaussExpDenom);
      }
    }
    intensity[k] = sum;
  }

  return { grid, intensity };
}

/**
 * Find the index of the mode whose frequency is closest to `clickedX`
 * within `thresholdCm1`. Returns -1 if no mode qualifies.
 *
 * Imaginary modes (frequency ≤ 0) are excluded because they are not part
 * of the displayed spectrum on a [0, 4500] cm⁻¹ axis.
 *
 * Uses linear search (N ≤ ~30 for typical molecules).
 */
export function findNearestMode(
  frequenciesCm1: readonly number[],
  clickedX: number,
  thresholdCm1: number = 50
): number {
  if (frequenciesCm1.length === 0) return -1;
  let bestIdx = -1;
  let bestDist = Infinity;
  for (let i = 0; i < frequenciesCm1.length; i++) {
    if (frequenciesCm1[i] <= 0) continue;
    const d = Math.abs(frequenciesCm1[i] - clickedX);
    if (d < bestDist) {
      bestDist = d;
      bestIdx = i;
    }
  }
  if (bestIdx === -1) return -1;
  return bestDist <= thresholdCm1 ? bestIdx : -1;
}
