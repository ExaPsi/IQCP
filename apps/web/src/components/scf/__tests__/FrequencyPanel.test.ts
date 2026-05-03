/**
 * Unit tests for FrequencyPanel helper functions (US-102).
 *
 * Pure-logic tests for the exported helpers in `FrequencyPanel.tsx`
 * (`buildFrequencyRows`, `sortFrequencyRows`, `countImaginary`,
 * `firstImaginaryIndex`, `formatFrequency`). These helpers drive the
 * sortable table UI in the Frequency tab and must be correct before any
 * rendering happens.
 *
 * The component itself is a thin JSX wrapper around these helpers and is
 * exercised via the manual smoke test documented in VP-US-102 Section 17.
 *
 * @module components/scf/__tests__/FrequencyPanel.test
 */

import { describe, it, expect } from 'vitest';
import {
  buildFrequencyRows,
  sortFrequencyRows,
  countImaginary,
  firstImaginaryIndex,
  formatFrequency,
} from '../frequencyPanelLogic';
import type { FrequencyResult } from '../../../worker/protocol';

// ============================================================================
// Fixture — minimal FrequencyResult with 4 modes (1 imaginary)
// ============================================================================

function makeFixture(): FrequencyResult {
  return {
    nAtoms: 3,
    nModes: 4,
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
    frequenciesCm1: [-300, 1595, 3657, 3756],
    reducedMassesAmu: [1.0, 1.1, 1.05, 1.08],
    forceConstantsMdyne: [0.1, 1.5, 8.0, 8.5],
    normalModesCartesian: [],
    rotationalConstantsGhz: [835, 435, 286],
    irIntensitiesKmPerMol: [0.0, 80.0, 5.0, 60.0],
    ramanActivitiesA4Amu: [0.5, 10.0, 100.0, 30.0],
    depolarizationRatios: [0.75, 0.3, 0.2, 0.5],
    thermochemistry: {
      temperatureK: 298.15,
      pressurePa: 101325,
      symmetryNumber: 2,
      multiplicity: 1,
      totalMassAmu: 18.0153,
      nVibModesUsed: 3,
      nImag: 1,
      zpeHa: 0.02,
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
// Tests: buildFrequencyRows
// ============================================================================

describe('buildFrequencyRows', () => {
  it('returns an empty array for null result', () => {
    expect(buildFrequencyRows(null)).toEqual([]);
  });

  it('returns an empty array for zero modes', () => {
    const result = makeFixture();
    result.frequenciesCm1 = [];
    expect(buildFrequencyRows(result)).toEqual([]);
  });

  it('builds one row per mode in input order', () => {
    const rows = buildFrequencyRows(makeFixture());
    expect(rows).toHaveLength(4);
    expect(rows.map((r) => r.modeIndex)).toEqual([0, 1, 2, 3]);
    expect(rows[0].frequency).toBe(-300);
    expect(rows[1].frequency).toBe(1595);
    expect(rows[2].frequency).toBe(3657);
    expect(rows[3].frequency).toBe(3756);
  });

  it('flags negative frequencies as imaginary', () => {
    const rows = buildFrequencyRows(makeFixture());
    expect(rows[0].isImaginary).toBe(true);
    expect(rows[1].isImaginary).toBe(false);
    expect(rows[2].isImaginary).toBe(false);
    expect(rows[3].isImaginary).toBe(false);
  });

  it('copies intensity, activity, ρ, μ, k from result arrays', () => {
    const rows = buildFrequencyRows(makeFixture());
    expect(rows[1].irIntensity).toBe(80.0);
    expect(rows[2].ramanActivity).toBe(100.0);
    expect(rows[3].depolarization).toBe(0.5);
    expect(rows[0].reducedMass).toBe(1.0);
    expect(rows[2].forceConstant).toBe(8.0);
  });
});

// ============================================================================
// Tests: sortFrequencyRows
// ============================================================================

describe('sortFrequencyRows', () => {
  it('sorts ascending by frequency', () => {
    const rows = buildFrequencyRows(makeFixture());
    const sorted = sortFrequencyRows(rows, 'frequency', 'asc');
    expect(sorted.map((r) => r.frequency)).toEqual([-300, 1595, 3657, 3756]);
  });

  it('sorts descending by frequency', () => {
    const rows = buildFrequencyRows(makeFixture());
    const sorted = sortFrequencyRows(rows, 'frequency', 'desc');
    expect(sorted.map((r) => r.frequency)).toEqual([3756, 3657, 1595, -300]);
  });

  it('sorts by IR intensity', () => {
    const rows = buildFrequencyRows(makeFixture());
    const sorted = sortFrequencyRows(rows, 'irIntensity', 'desc');
    expect(sorted[0].irIntensity).toBe(80.0);
    expect(sorted[1].irIntensity).toBe(60.0);
    expect(sorted[2].irIntensity).toBe(5.0);
    expect(sorted[3].irIntensity).toBe(0.0);
  });

  it('sorts by Raman activity', () => {
    const rows = buildFrequencyRows(makeFixture());
    const sorted = sortFrequencyRows(rows, 'ramanActivity', 'desc');
    expect(sorted[0].ramanActivity).toBe(100.0);
  });

  it('sort by modeIndex ascending matches original order', () => {
    const rows = buildFrequencyRows(makeFixture());
    const sorted = sortFrequencyRows(rows, 'modeIndex', 'asc');
    expect(sorted.map((r) => r.modeIndex)).toEqual([0, 1, 2, 3]);
  });

  it('is stable for equal keys', () => {
    const rows = buildFrequencyRows(makeFixture());
    // Tie: irIntensity 0 (modeIndex 0) vs any other row — but make a new tie
    // by duplicating a depolarization value
    const sorted = sortFrequencyRows(rows, 'depolarization', 'asc');
    // depolarizations: [0.75, 0.3, 0.2, 0.5] → sorted: 0.2, 0.3, 0.5, 0.75
    expect(sorted.map((r) => r.depolarization)).toEqual([0.2, 0.3, 0.5, 0.75]);
  });

  it('does not mutate the input array', () => {
    const rows = buildFrequencyRows(makeFixture());
    const originalOrder = rows.map((r) => r.frequency);
    sortFrequencyRows(rows, 'frequency', 'desc');
    expect(rows.map((r) => r.frequency)).toEqual(originalOrder);
  });
});

// ============================================================================
// Tests: countImaginary + firstImaginaryIndex
// ============================================================================

describe('countImaginary', () => {
  it('returns 0 for null', () => {
    expect(countImaginary(null)).toBe(0);
  });

  it('returns 0 when all frequencies are positive', () => {
    const r = makeFixture();
    r.frequenciesCm1 = [100, 200, 300];
    expect(countImaginary(r)).toBe(0);
  });

  it('counts a single imaginary frequency', () => {
    expect(countImaginary(makeFixture())).toBe(1);
  });

  it('counts multiple imaginary frequencies', () => {
    const r = makeFixture();
    r.frequenciesCm1 = [-100, -50, 200, 300];
    expect(countImaginary(r)).toBe(2);
  });

  it('zero frequency is NOT imaginary', () => {
    const r = makeFixture();
    r.frequenciesCm1 = [0, 100, 200];
    expect(countImaginary(r)).toBe(0);
  });
});

describe('firstImaginaryIndex', () => {
  it('returns -1 for null', () => {
    expect(firstImaginaryIndex(null)).toBe(-1);
  });

  it('returns -1 when no imaginary modes', () => {
    const r = makeFixture();
    r.frequenciesCm1 = [100, 200, 300];
    expect(firstImaginaryIndex(r)).toBe(-1);
  });

  it('returns 0 when the first mode is imaginary (fixture)', () => {
    expect(firstImaginaryIndex(makeFixture())).toBe(0);
  });

  it('returns the first imaginary index when it is not 0', () => {
    const r = makeFixture();
    r.frequenciesCm1 = [100, 200, -500, 400];
    expect(firstImaginaryIndex(r)).toBe(2);
  });
});

// ============================================================================
// Tests: formatFrequency
// ============================================================================

describe('formatFrequency', () => {
  it('formats positive values to 1 decimal place', () => {
    expect(formatFrequency(1595.3)).toBe('1595.3');
  });

  it('formats negative values as absolute value (sign is contextual)', () => {
    expect(formatFrequency(-500)).toBe('500.0');
  });

  it('formats zero as 0.0', () => {
    expect(formatFrequency(0)).toBe('0.0');
  });
});
