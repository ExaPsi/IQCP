/**
 * Unit tests for FrequencySummaryPanel helper functions (US-102).
 *
 * Pure-logic tests for the exported helpers (`dipoleMagnitude`,
 * `isotropicPolarizability`, `formatRotorType`, `formatRotConst`).
 *
 * @module components/scf/__tests__/FrequencySummaryPanel.test
 */

import { describe, it, expect } from 'vitest';
import {
  dipoleMagnitude,
  isotropicPolarizability,
  formatRotorType,
  formatRotConst,
} from '../frequencySummaryLogic';

// ============================================================================
// dipoleMagnitude
// ============================================================================

describe('dipoleMagnitude', () => {
  it('computes the L2 norm of a 3-vector', () => {
    expect(dipoleMagnitude([1, 0, 0])).toBeCloseTo(1, 15);
    expect(dipoleMagnitude([0, 2, 0])).toBeCloseTo(2, 15);
    expect(dipoleMagnitude([3, 4, 0])).toBeCloseTo(5, 15);
    expect(dipoleMagnitude([1, 2, 2])).toBeCloseTo(3, 15);
  });

  it('returns 0 for the zero vector', () => {
    expect(dipoleMagnitude([0, 0, 0])).toBe(0);
  });

  it('handles negative components', () => {
    expect(dipoleMagnitude([-3, -4, 0])).toBeCloseTo(5, 15);
  });

  it('H₂O-like dipole ~1.85 Debye', () => {
    // Experimental H₂O dipole: 1.8546 D, directed along the C₂v axis (y-ish)
    const mu = dipoleMagnitude([0.1, 1.85, 0.2]);
    expect(mu).toBeGreaterThan(1.85);
    expect(mu).toBeLessThan(2.0);
  });
});

// ============================================================================
// isotropicPolarizability
// ============================================================================

describe('isotropicPolarizability', () => {
  it('returns (αxx + αyy + αzz) / 3 for a diagonal tensor', () => {
    const tensor: [
      [number, number, number],
      [number, number, number],
      [number, number, number],
    ] = [
      [9, 0, 0],
      [0, 6, 0],
      [0, 0, 3],
    ];
    expect(isotropicPolarizability(tensor)).toBeCloseTo(6, 15);
  });

  it('ignores off-diagonal elements', () => {
    const tensor: [
      [number, number, number],
      [number, number, number],
      [number, number, number],
    ] = [
      [1, 99, 99],
      [99, 2, 99],
      [99, 99, 3],
    ];
    expect(isotropicPolarizability(tensor)).toBeCloseTo(2, 15);
  });

  it('returns 0 for the zero tensor', () => {
    const tensor: [
      [number, number, number],
      [number, number, number],
      [number, number, number],
    ] = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];
    expect(isotropicPolarizability(tensor)).toBe(0);
  });
});

// ============================================================================
// formatRotorType
// ============================================================================

describe('formatRotorType', () => {
  it('maps every RotorType variant to a human-readable label', () => {
    expect(formatRotorType('atom')).toBe('Atom');
    expect(formatRotorType('linear')).toBe('Linear');
    expect(formatRotorType('spherical_top')).toBe('Spherical top');
    expect(formatRotorType('symmetric_top')).toBe('Symmetric top');
    expect(formatRotorType('asymmetric_top')).toBe('Asymmetric top');
  });
});

// ============================================================================
// formatRotConst
// ============================================================================

describe('formatRotConst', () => {
  it('formats a finite value to 3 decimal places', () => {
    expect(formatRotConst(12.3456)).toBe('12.346');
    expect(formatRotConst(1.0)).toBe('1.000');
    expect(formatRotConst(0.123)).toBe('0.123');
  });

  it('returns "—" for Infinity', () => {
    expect(formatRotConst(Infinity)).toBe('—');
    expect(formatRotConst(-Infinity)).toBe('—');
  });

  it('returns "—" for NaN', () => {
    expect(formatRotConst(NaN)).toBe('—');
  });

  it('handles zero', () => {
    expect(formatRotConst(0)).toBe('0.000');
  });
});
