/**
 * Unit tests for ThermochemistryPanel helper functions (US-102).
 *
 * Pure-logic tests for the exported format helpers + the supporting
 * `lib/units.ts` helpers. The component itself is exercised via the manual
 * smoke test documented in VP-US-102 Section 17.
 *
 * @module components/scf/__tests__/ThermochemistryPanel.test
 */

import { describe, it, expect } from 'vitest';
import {
  formatEnergyValue,
  formatEntropyValue,
} from '../thermochemistryFormat';
import {
  convertEnergy,
  convertEntropyPerK,
  energyUnitsLabel,
  entropyUnitsLabel,
  HA_TO_KCAL_MOL,
  HA_TO_KJ_MOL,
} from '../../../lib/units';

// ============================================================================
// convertEnergy
// ============================================================================

describe('convertEnergy', () => {
  it('returns the same value for hartree', () => {
    expect(convertEnergy(1.0, 'hartree')).toBe(1.0);
    expect(convertEnergy(-76.0, 'hartree')).toBe(-76.0);
  });

  it('converts to kcal/mol using 627.5094740631', () => {
    expect(convertEnergy(1.0, 'kcal_mol')).toBeCloseTo(HA_TO_KCAL_MOL, 12);
    expect(convertEnergy(2.0, 'kcal_mol')).toBeCloseTo(2 * HA_TO_KCAL_MOL, 12);
  });

  it('converts to kJ/mol using 2625.4996394799', () => {
    expect(convertEnergy(1.0, 'kj_mol')).toBeCloseTo(HA_TO_KJ_MOL, 12);
  });

  it('is linear in the input', () => {
    const a = convertEnergy(0.5, 'kcal_mol');
    const b = convertEnergy(1.0, 'kcal_mol');
    expect(b).toBeCloseTo(2 * a, 12);
  });

  it('preserves sign', () => {
    expect(convertEnergy(-1.0, 'kj_mol')).toBeLessThan(0);
    expect(convertEnergy(-1.0, 'kj_mol')).toBeCloseTo(-HA_TO_KJ_MOL, 12);
  });
});

// ============================================================================
// convertEntropyPerK
// ============================================================================

describe('convertEntropyPerK', () => {
  it('uses the same conversion factors as convertEnergy', () => {
    expect(convertEntropyPerK(3.166811e-6, 'hartree')).toBe(3.166811e-6);
    expect(convertEntropyPerK(3.166811e-6, 'kcal_mol')).toBeCloseTo(
      3.166811e-6 * HA_TO_KCAL_MOL,
      12
    );
  });
});

// ============================================================================
// Unit labels
// ============================================================================

describe('energyUnitsLabel', () => {
  it('maps hartree → "Ha"', () => {
    expect(energyUnitsLabel('hartree')).toBe('Ha');
  });
  it('maps kcal_mol → "kcal/mol"', () => {
    expect(energyUnitsLabel('kcal_mol')).toBe('kcal/mol');
  });
  it('maps kj_mol → "kJ/mol"', () => {
    expect(energyUnitsLabel('kj_mol')).toBe('kJ/mol');
  });
});

describe('entropyUnitsLabel', () => {
  it('includes the per-mol·K denominator', () => {
    expect(entropyUnitsLabel('hartree')).toBe('Ha/(mol·K)');
    expect(entropyUnitsLabel('kcal_mol')).toBe('kcal/(mol·K)');
    expect(entropyUnitsLabel('kj_mol')).toBe('kJ/(mol·K)');
  });
});

// ============================================================================
// formatEnergyValue
// ============================================================================

describe('formatEnergyValue', () => {
  it('formats Hartree values with 6 decimal places', () => {
    expect(formatEnergyValue(1.234567890123, 'hartree')).toBe('1.234568');
  });

  it('formats kcal/mol values with 3 decimal places', () => {
    expect(formatEnergyValue(0.025, 'kcal_mol')).toBe(
      (0.025 * HA_TO_KCAL_MOL).toFixed(3)
    );
  });

  it('formats kJ/mol values with 3 decimal places', () => {
    expect(formatEnergyValue(0.025, 'kj_mol')).toBe(
      (0.025 * HA_TO_KJ_MOL).toFixed(3)
    );
  });

  it('returns "—" for non-finite input', () => {
    expect(formatEnergyValue(NaN, 'hartree')).toBe('—');
    expect(formatEnergyValue(Infinity, 'hartree')).toBe('—');
  });

  it('handles zero', () => {
    expect(formatEnergyValue(0, 'hartree')).toBe('0.000000');
    expect(formatEnergyValue(0, 'kcal_mol')).toBe('0.000');
  });

  it('handles negative values', () => {
    expect(formatEnergyValue(-76.0, 'hartree')).toBe('-76.000000');
  });
});

// ============================================================================
// formatEntropyValue
// ============================================================================

describe('formatEntropyValue', () => {
  it('formats Hartree entropy in scientific notation', () => {
    const s = 7.179e-5;
    const out = formatEntropyValue(s, 'hartree');
    expect(out).toMatch(/e[+-]?\d/);
  });

  it('formats kcal/(mol·K) in fixed notation', () => {
    const s = 7.179e-5;
    const out = formatEntropyValue(s, 'kcal_mol');
    expect(out).not.toMatch(/e[+-]?\d/);
    // 7.179e-5 × 627.509 ≈ 0.04504890... The panel rounds to 5 decimals
    // → "0.04505"; verify within ±5e-6 (half a unit of the last decimal).
    expect(Math.abs(Number(out) - s * HA_TO_KCAL_MOL)).toBeLessThan(5e-6);
  });

  it('returns "—" for non-finite values', () => {
    expect(formatEntropyValue(NaN, 'hartree')).toBe('—');
    expect(formatEntropyValue(Infinity, 'kcal_mol')).toBe('—');
  });

  it('handles zero', () => {
    expect(formatEntropyValue(0, 'hartree')).toMatch(/0/);
    expect(formatEntropyValue(0, 'kj_mol')).toBe('0.00000');
  });
});

// ============================================================================
// Cross-consistency: round-trip
// ============================================================================

describe('convertEnergy round-trip', () => {
  it('hartree → kcal/mol → hartree returns the original', () => {
    const original = -76.02341;
    const kcal = convertEnergy(original, 'kcal_mol');
    const back = kcal / HA_TO_KCAL_MOL;
    expect(back).toBeCloseTo(original, 12);
  });

  it('HA_TO_KJ_MOL ≈ HA_TO_KCAL_MOL × 4.184 (thermochemical cal/kcal)', () => {
    // This relationship is fundamental: 1 kcal (thermochem) = 4.184 kJ exactly
    expect(HA_TO_KJ_MOL / HA_TO_KCAL_MOL).toBeCloseTo(4.184, 8);
  });
});
