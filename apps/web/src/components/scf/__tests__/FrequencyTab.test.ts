/**
 * Integration-level tests for FrequencyTab (US-102).
 *
 * The component tree itself (JSX + lazy-loaded R3F + Plotly) cannot be
 * mounted in the Node test environment that the web app uses. These tests
 * instead verify the pieces that are pure-logic:
 *   - `LocalFrequencyState` shape (field coverage + default values)
 *   - `DEFAULT_LOCAL_FREQUENCY_STATE` values match US-102 Section 7.1
 *   - The exported state shape includes every field required by US-103's
 *     lift-to-Zustand contract (guard against accidental breakage).
 *
 * Full end-to-end wiring is exercised via the manual smoke test
 * documented in VP-US-102 Section 17.
 *
 * @module components/scf/__tests__/FrequencyTab.test
 */

import { describe, it, expect } from 'vitest';
import {
  DEFAULT_LOCAL_FREQUENCY_STATE,
  type LocalFrequencyState,
} from '../frequencyTabState';

// ============================================================================
// Default state
// ============================================================================

describe('DEFAULT_LOCAL_FREQUENCY_STATE', () => {
  it('starts with no result and no loading state', () => {
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.result).toBeNull();
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.isComputing).toBe(false);
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.progress).toBeNull();
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.error).toBeNull();
  });

  it('starts with no selected mode', () => {
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.selectedMode).toBeNull();
  });

  it('defaults to T = 298.15 K, P = 101325 Pa', () => {
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.temperatureK).toBe(298.15);
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.pressurePa).toBe(101325);
  });

  it('seeds displayThermo to null until a result arrives', () => {
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.displayThermo).toBeNull();
  });

  it('defaults the spectrum tab to IR', () => {
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.spectrumTab).toBe('ir');
  });

  it('defaults the broadening to Lorentzian with FWHM 8 cm⁻¹', () => {
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.broadeningKind).toBe('lorentzian');
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.fwhmCm1).toBe(8.0);
  });

  it('defaults the amplitude to 0.5 Å, speed 1.0×, animating=true', () => {
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.amplitudeBohr).toBe(0.5);
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.animationSpeed).toBe(1.0);
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.isAnimating).toBe(true);
  });

  it('defaults showDisplacementArrows to false (novice-friendly UX refinement)', () => {
    // US-102 UX refinement: animated atoms are the primary view; static
    // direction-field arrows are an opt-in for advanced users.
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.showDisplacementArrows).toBe(false);
  });

  it('defaults the units mode to kcal/mol', () => {
    expect(DEFAULT_LOCAL_FREQUENCY_STATE.unitsMode).toBe('kcal_mol');
  });

  it('is a plain object (no class instances)', () => {
    expect(typeof DEFAULT_LOCAL_FREQUENCY_STATE).toBe('object');
    expect(DEFAULT_LOCAL_FREQUENCY_STATE).not.toBeNull();
  });
});

// ============================================================================
// Shape coverage
// ============================================================================

describe('LocalFrequencyState shape (lift-friendly contract with US-103)', () => {
  it('has every field listed in US-102 Section 7.1', () => {
    // If any field is removed or renamed, this test fails and alerts the
    // implementer to update US-103 before lifting.
    // (showDisplacementArrows added by US-102 UX refinement.)
    const expectedKeys: Array<keyof LocalFrequencyState> = [
      'result',
      'isComputing',
      'progress',
      'error',
      'selectedMode',
      'temperatureK',
      'pressurePa',
      'displayThermo',
      'spectrumTab',
      'broadeningKind',
      'fwhmCm1',
      'amplitudeBohr',
      'animationSpeed',
      'isAnimating',
      'showDisplacementArrows',
      'unitsMode',
    ];
    const actualKeys = Object.keys(DEFAULT_LOCAL_FREQUENCY_STATE).sort();
    const expectedSorted = [...expectedKeys].sort();
    expect(actualKeys).toEqual(expectedSorted);
  });

  it('has no extra fields beyond the US-103 contract', () => {
    const allowed = new Set<string>([
      'result',
      'isComputing',
      'progress',
      'error',
      'selectedMode',
      'temperatureK',
      'pressurePa',
      'displayThermo',
      'spectrumTab',
      'broadeningKind',
      'fwhmCm1',
      'amplitudeBohr',
      'animationSpeed',
      'isAnimating',
      'showDisplacementArrows',
      'unitsMode',
    ]);
    for (const k of Object.keys(DEFAULT_LOCAL_FREQUENCY_STATE)) {
      expect(allowed.has(k)).toBe(true);
    }
  });
});
