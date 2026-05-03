/**
 * Unit tests for SpectrumPlot helper functions (US-102).
 *
 * Pure-logic tests for the exported helpers in `SpectrumPlot.tsx`
 * (`buildWavenumberGrid`, `computeBroadenedSpectrum`, `findNearestMode`).
 * The component itself is a thin JSX/Plotly wrapper and is exercised via
 * the manual smoke test documented in VP-US-102 Section 17.
 *
 * Kernel normalizations are validated both against analytic values
 * (peak center, half-width symmetry) and integrated area (trapezoidal
 * integration of the broadened curve returns the original stick amplitude
 * within <1% error at 8 points per FWHM).
 *
 * @module components/scf/__tests__/SpectrumPlot.test
 */

import { describe, it, expect } from 'vitest';
import {
  buildWavenumberGrid,
  computeBroadenedSpectrum,
  findNearestMode,
  SPECTRUM_GRID_POINTS,
} from '../spectrumMath';

// ============================================================================
// Helper: trapezoidal integration
// ============================================================================

function trapz(x: number[], y: number[]): number {
  let sum = 0;
  for (let i = 0; i < x.length - 1; i++) {
    sum += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i]);
  }
  return sum;
}

// ============================================================================
// buildWavenumberGrid
// ============================================================================

describe('buildWavenumberGrid', () => {
  it('has exactly 4501 points', () => {
    const grid = buildWavenumberGrid();
    expect(grid).toHaveLength(SPECTRUM_GRID_POINTS);
    expect(grid).toHaveLength(4501);
  });

  it('spans 0..=4500 in steps of 1', () => {
    const grid = buildWavenumberGrid();
    expect(grid[0]).toBe(0);
    expect(grid[grid.length - 1]).toBe(4500);
    expect(grid[1] - grid[0]).toBe(1);
  });
});

// ============================================================================
// computeBroadenedSpectrum — Lorentzian
// ============================================================================

describe('computeBroadenedSpectrum — Lorentzian', () => {
  it('returns an all-zero intensity array for empty input', () => {
    const { grid, intensity } = computeBroadenedSpectrum(
      [],
      [],
      'lorentzian',
      8
    );
    expect(grid).toHaveLength(4501);
    expect(intensity.every((v) => v === 0)).toBe(true);
  });

  it('returns zeros for zero FWHM', () => {
    const { intensity } = computeBroadenedSpectrum(
      [1000],
      [1],
      'lorentzian',
      0
    );
    expect(intensity.every((v) => v === 0)).toBe(true);
  });

  it('peak center value for Lorentzian is A/(π·γ) with γ = FWHM/2', () => {
    const fwhm = 8;
    const { intensity } = computeBroadenedSpectrum(
      [1000],
      [1],
      'lorentzian',
      fwhm
    );
    const gamma = fwhm / 2;
    const peakExpected = 1 / (Math.PI * gamma);
    // Grid point exactly at 1000
    expect(intensity[1000]).toBeCloseTo(peakExpected, 12);
  });

  it('Lorentzian falls to exactly half peak at ±FWHM/2', () => {
    const fwhm = 10;
    const { intensity } = computeBroadenedSpectrum(
      [500],
      [1],
      'lorentzian',
      fwhm
    );
    const peak = intensity[500];
    const half = intensity[495]; // 500 - 5
    expect(half).toBeCloseTo(peak / 2, 12);
  });

  it('Lorentzian area ≈ stick amplitude (trapezoidal integration)', () => {
    const stickAmp = 42;
    const fwhm = 8;
    const { grid, intensity } = computeBroadenedSpectrum(
      [1500],
      [stickAmp],
      'lorentzian',
      fwhm
    );
    const area = trapz(grid, intensity);
    // Lorentzian tails are heavy — integrate only the stick region, but
    // the [0, 4500] range at 1 cm⁻¹ resolution with FWHM = 8 captures
    // well over 99% of the area. Accept 1% relative error (matches qc-core
    // reconstruction accuracy at 8 points per FWHM).
    expect(area).toBeGreaterThan(stickAmp * 0.98);
    expect(area).toBeLessThan(stickAmp * 1.02);
  });

  it('skips imaginary (negative) frequencies', () => {
    const { intensity } = computeBroadenedSpectrum(
      [-500, 1000, -200, 2000],
      [1, 1, 1, 1],
      'lorentzian',
      8
    );
    // Area should be ~2 (two real modes), not 4
    const grid = buildWavenumberGrid();
    const area = trapz(grid, intensity);
    expect(area).toBeGreaterThan(1.95);
    expect(area).toBeLessThan(2.05);
  });

  it('linear in stick amplitude: 2·A produces 2× intensity', () => {
    const { intensity: i1 } = computeBroadenedSpectrum(
      [1500],
      [1],
      'lorentzian',
      8
    );
    const { intensity: i2 } = computeBroadenedSpectrum(
      [1500],
      [2],
      'lorentzian',
      8
    );
    expect(i2[1500]).toBeCloseTo(2 * i1[1500], 12);
  });
});

// ============================================================================
// computeBroadenedSpectrum — Gaussian
// ============================================================================

describe('computeBroadenedSpectrum — Gaussian', () => {
  it('peak center value is A / (σ√(2π)) with σ = FWHM/(2√(2 ln 2))', () => {
    const fwhm = 10;
    const { intensity } = computeBroadenedSpectrum(
      [2000],
      [1],
      'gaussian',
      fwhm
    );
    const sigma = fwhm / (2 * Math.sqrt(2 * Math.log(2)));
    const peakExpected = 1 / (sigma * Math.sqrt(2 * Math.PI));
    expect(intensity[2000]).toBeCloseTo(peakExpected, 12);
  });

  it('Gaussian falls to half peak at ±FWHM/2 (by definition of FWHM)', () => {
    const fwhm = 20;
    const { intensity } = computeBroadenedSpectrum(
      [3000],
      [1],
      'gaussian',
      fwhm
    );
    const peak = intensity[3000];
    const half = intensity[2990]; // 3000 - 10
    expect(half).toBeCloseTo(peak / 2, 12);
  });

  it('Gaussian area ≈ stick amplitude', () => {
    const stickAmp = 7;
    const fwhm = 8;
    const { grid, intensity } = computeBroadenedSpectrum(
      [2500],
      [stickAmp],
      'gaussian',
      fwhm
    );
    const area = trapz(grid, intensity);
    // Gaussian tails decay fast; at 8 points per FWHM, area error < 0.5%.
    expect(area).toBeGreaterThan(stickAmp * 0.995);
    expect(area).toBeLessThan(stickAmp * 1.005);
  });

  it('multiple modes → curves add linearly', () => {
    const single1 = computeBroadenedSpectrum([1000], [1], 'gaussian', 10);
    const single2 = computeBroadenedSpectrum([2000], [1], 'gaussian', 10);
    const both = computeBroadenedSpectrum(
      [1000, 2000],
      [1, 1],
      'gaussian',
      10
    );
    for (const k of [1000, 1500, 2000]) {
      expect(both.intensity[k]).toBeCloseTo(
        single1.intensity[k] + single2.intensity[k],
        12
      );
    }
  });
});

// ============================================================================
// findNearestMode
// ============================================================================

describe('findNearestMode', () => {
  it('returns -1 for empty frequency list', () => {
    expect(findNearestMode([], 1000)).toBe(-1);
  });

  it('returns the only mode if it is within threshold', () => {
    expect(findNearestMode([1500], 1510)).toBe(0);
  });

  it('returns -1 if the clicked x is > threshold from any mode', () => {
    expect(findNearestMode([1500], 2000)).toBe(-1); // 500 cm⁻¹ away
  });

  it('picks the closest real mode when within threshold', () => {
    const freqs = [500, 1500, 2500, 3500];
    expect(findNearestMode(freqs, 1520)).toBe(1); // 20 cm⁻¹ from 1500
    expect(findNearestMode(freqs, 2490)).toBe(2); // 10 cm⁻¹ from 2500
    // Click at exactly 2000 is 500 cm⁻¹ from both 1500 and 2500 → outside threshold
    expect(findNearestMode(freqs, 2000)).toBe(-1);
  });

  it('skips imaginary (negative) frequencies', () => {
    const freqs = [-1000, 1500, 2500];
    // Click at -900 is close to the imaginary -1000 but imaginary should be skipped
    expect(findNearestMode(freqs, -900, 200)).toBe(-1);
  });

  it('accepts a custom threshold', () => {
    const freqs = [1000];
    expect(findNearestMode(freqs, 1200, 50)).toBe(-1);
    expect(findNearestMode(freqs, 1200, 250)).toBe(0);
  });

  it('returns the first tied mode (linear search, stable)', () => {
    const freqs = [1000, 1200];
    // Clicked at 1100 — both are 100 away
    expect(findNearestMode(freqs, 1100, 150)).toBe(0);
  });
});
