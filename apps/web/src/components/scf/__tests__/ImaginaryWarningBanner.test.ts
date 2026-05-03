/**
 * Unit tests for the imaginary-warning-banner decision logic (US-102).
 *
 * The banner component itself (JSX + Tailwind classes) is exercised via the
 * manual smoke test documented in VP-US-102 Section 17. These tests verify
 * the two pure helpers the banner depends on (`countImaginary` and
 * `firstImaginaryIndex`) under scenarios specific to the banner's rendering
 * decisions:
 *   - Banner hidden when no imaginary modes
 *   - Banner hidden when result is null
 *   - Banner shown with count = 1 for a single imaginary
 *   - Banner shown with count ≥ 2 for multiple imaginary
 *   - "Show mode" button target index = firstImaginaryIndex
 *
 * @module components/scf/__tests__/ImaginaryWarningBanner.test
 */

import { describe, it, expect } from 'vitest';
import {
  countImaginary,
  firstImaginaryIndex,
} from '../frequencyPanelLogic';
import type { FrequencyResult } from '../../../worker/protocol';

function makeResultWithFreqs(freqs: number[]): FrequencyResult {
  return {
    nAtoms: 3,
    nModes: freqs.length,
    rotorType: 'asymmetric_top',
    electronicEnergyHa: -76.0,
    dipoleAu: [0, 0, 0],
    dipoleDebye: [0, 0, 0],
    polarizabilityAu: [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ],
    polarizabilityAng3: [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ],
    frequenciesCm1: freqs,
    reducedMassesAmu: freqs.map(() => 1),
    forceConstantsMdyne: freqs.map(() => 1),
    normalModesCartesian: [],
    rotationalConstantsGhz: [835, 435, 286],
    irIntensitiesKmPerMol: freqs.map(() => 0),
    ramanActivitiesA4Amu: freqs.map(() => 0),
    depolarizationRatios: freqs.map(() => 0),
    thermochemistry: {
      temperatureK: 298.15,
      pressurePa: 101325,
      symmetryNumber: 1,
      multiplicity: 1,
      totalMassAmu: 18,
      nVibModesUsed: 0,
      nImag: 0,
      zpeHa: 0,
      e0kHa: 0,
      internalEnergyHa: 0,
      enthalpyHa: 0,
      entropyHaPerK: 0,
      gibbsHa: 0,
      cvHaPerK: 0,
      cpHaPerK: 0,
      eTransHa: 0,
      hTransHa: 0,
      sTransHaPerK: 0,
      cvTransHaPerK: 0,
      cpTransHaPerK: 0,
      eRotHa: 0,
      hRotHa: 0,
      sRotHaPerK: 0,
      cvRotHaPerK: 0,
      cpRotHaPerK: 0,
      eVibThermalHa: 0,
      hVibHa: 0,
      sVibHaPerK: 0,
      cvVibHaPerK: 0,
      cpVibHaPerK: 0,
      sElecHaPerK: 0,
    },
    irSpectrum: {
      wavenumbersCm1: [],
      intensity: [],
      kind: 'lorentzian',
      fwhmCm1: 8,
    },
    ramanSpectrum: {
      wavenumbersCm1: [],
      intensity: [],
      kind: 'lorentzian',
      fwhmCm1: 8,
    },
    timingMs: {
      integralsMs: 0,
      nuclearCphfMs: 0,
      fieldCphfMs: 0,
      assemblyMs: 0,
      modesMs: 0,
      totalMs: 0,
    },
    aborted: false,
  };
}

// ============================================================================
// Banner visibility logic
// ============================================================================

describe('ImaginaryWarningBanner visibility logic', () => {
  it('hidden when result is null', () => {
    // The banner renders null when countImaginary(result) === 0.
    expect(countImaginary(null)).toBe(0);
  });

  it('hidden when all frequencies are positive', () => {
    const r = makeResultWithFreqs([100, 200, 300]);
    expect(countImaginary(r)).toBe(0);
  });

  it('hidden when the frequency list is empty', () => {
    const r = makeResultWithFreqs([]);
    expect(countImaginary(r)).toBe(0);
  });

  it('shown with count=1 for a single imaginary', () => {
    const r = makeResultWithFreqs([-200, 1000, 2000]);
    expect(countImaginary(r)).toBe(1);
  });

  it('shown with count=2 for two imaginary', () => {
    const r = makeResultWithFreqs([-500, -200, 1000, 2000]);
    expect(countImaginary(r)).toBe(2);
  });

  it('shown with count=3 for an entirely imaginary spectrum', () => {
    const r = makeResultWithFreqs([-100, -200, -300]);
    expect(countImaginary(r)).toBe(3);
  });
});

// ============================================================================
// "Show mode" button target
// ============================================================================

describe('ImaginaryWarningBanner "Show mode" target', () => {
  it('returns -1 when no imaginary (banner would be hidden anyway)', () => {
    const r = makeResultWithFreqs([100, 200, 300]);
    expect(firstImaginaryIndex(r)).toBe(-1);
  });

  it('returns 0 when the first mode is imaginary', () => {
    const r = makeResultWithFreqs([-500, 1000, 2000]);
    expect(firstImaginaryIndex(r)).toBe(0);
  });

  it('returns the first imaginary index when the imaginary is in the middle', () => {
    const r = makeResultWithFreqs([1000, -400, 2000, 3000]);
    expect(firstImaginaryIndex(r)).toBe(1);
  });

  it('returns the first imaginary index when there are multiple imaginaries', () => {
    const r = makeResultWithFreqs([1000, -400, -300, 2000]);
    expect(firstImaginaryIndex(r)).toBe(1);
  });

  it('returns the last index when only the last mode is imaginary', () => {
    const r = makeResultWithFreqs([1000, 2000, 3000, -100]);
    expect(firstImaginaryIndex(r)).toBe(3);
  });
});

// ============================================================================
// Pluralization — banner text shows singular/plural based on count
// ============================================================================

describe('ImaginaryWarningBanner plural form decision', () => {
  it('single → plural flag is false', () => {
    const count = countImaginary(makeResultWithFreqs([-100, 1000]));
    expect(count !== 1).toBe(false);
  });

  it('multiple → plural flag is true', () => {
    const count = countImaginary(makeResultWithFreqs([-100, -200, 1000]));
    expect(count !== 1).toBe(true);
  });
});
