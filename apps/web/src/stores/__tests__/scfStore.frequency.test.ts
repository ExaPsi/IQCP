/**
 * Unit tests for the FrequencyState slice in scfStore (US-103).
 *
 * Tests:
 *   - DEFAULT_FREQUENCY_STATE shape and values
 *   - All 17 store actions (16 setters + 1 reset)
 *   - Tab switch persistence (state survives after simulated unmount)
 *   - Global reset includes frequency state
 *
 * @module stores/__tests__/scfStore.frequency.test
 * @see US-103 Frequency State + Deep Links
 */

import { describe, it, expect, beforeEach } from 'vitest';
import { useScfStore, DEFAULT_FREQUENCY_STATE } from '../scfStore';

// ============================================================================
// Setup
// ============================================================================

beforeEach(() => {
  // Reset store to default state before each test
  useScfStore.getState().reset();
});

// ============================================================================
// DEFAULT_FREQUENCY_STATE shape and values
// ============================================================================

describe('DEFAULT_FREQUENCY_STATE', () => {
  it('has exactly 16 keys', () => {
    expect(Object.keys(DEFAULT_FREQUENCY_STATE).length).toBe(16);
  });

  it('has all 16 expected field names', () => {
    const expectedKeys = [
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
    ].sort();
    const actualKeys = Object.keys(DEFAULT_FREQUENCY_STATE).sort();
    expect(actualKeys).toEqual(expectedKeys);
  });

  it('has amplitudeBohr (not amplitudeAngstrom)', () => {
    expect('amplitudeBohr' in DEFAULT_FREQUENCY_STATE).toBe(true);
    expect('amplitudeAngstrom' in DEFAULT_FREQUENCY_STATE).toBe(false);
  });

  it('result is null', () => {
    expect(DEFAULT_FREQUENCY_STATE.result).toBeNull();
  });

  it('isComputing is false', () => {
    expect(DEFAULT_FREQUENCY_STATE.isComputing).toBe(false);
  });

  it('progress is null', () => {
    expect(DEFAULT_FREQUENCY_STATE.progress).toBeNull();
  });

  it('error is null', () => {
    expect(DEFAULT_FREQUENCY_STATE.error).toBeNull();
  });

  it('selectedMode is null', () => {
    expect(DEFAULT_FREQUENCY_STATE.selectedMode).toBeNull();
  });

  it('temperatureK is 298.15', () => {
    expect(DEFAULT_FREQUENCY_STATE.temperatureK).toBe(298.15);
  });

  it('pressurePa is 101325', () => {
    expect(DEFAULT_FREQUENCY_STATE.pressurePa).toBe(101325);
  });

  it('displayThermo is null', () => {
    expect(DEFAULT_FREQUENCY_STATE.displayThermo).toBeNull();
  });

  it('spectrumTab is ir', () => {
    expect(DEFAULT_FREQUENCY_STATE.spectrumTab).toBe('ir');
  });

  it('broadeningKind is lorentzian', () => {
    expect(DEFAULT_FREQUENCY_STATE.broadeningKind).toBe('lorentzian');
  });

  it('fwhmCm1 is 8.0', () => {
    expect(DEFAULT_FREQUENCY_STATE.fwhmCm1).toBe(8.0);
  });

  it('amplitudeBohr is 0.5', () => {
    expect(DEFAULT_FREQUENCY_STATE.amplitudeBohr).toBe(0.5);
  });

  it('animationSpeed is 1.0', () => {
    expect(DEFAULT_FREQUENCY_STATE.animationSpeed).toBe(1.0);
  });

  it('isAnimating is true', () => {
    expect(DEFAULT_FREQUENCY_STATE.isAnimating).toBe(true);
  });

  it('showDisplacementArrows is false', () => {
    expect(DEFAULT_FREQUENCY_STATE.showDisplacementArrows).toBe(false);
  });

  it('unitsMode is kcal_mol', () => {
    expect(DEFAULT_FREQUENCY_STATE.unitsMode).toBe('kcal_mol');
  });
});

// ============================================================================
// Store integration
// ============================================================================

describe('Store integration', () => {
  it('frequencyState exists in the store', () => {
    expect(useScfStore.getState().frequencyState).toBeDefined();
  });

  it('initial frequencyState matches DEFAULT_FREQUENCY_STATE', () => {
    const state = useScfStore.getState().frequencyState;
    expect(state).toEqual(DEFAULT_FREQUENCY_STATE);
  });
});

// ============================================================================
// Action unit tests
// ============================================================================

describe('Frequency actions', () => {
  it('setFrequencyResult sets result', () => {
    const mockResult = { frequenciesCm1: [100, 200] } as never;
    useScfStore.getState().setFrequencyResult(mockResult);
    expect(useScfStore.getState().frequencyState.result).toBe(mockResult);
  });

  it('setFrequencyResult(null) clears result', () => {
    useScfStore.getState().setFrequencyResult({ frequenciesCm1: [] } as never);
    useScfStore.getState().setFrequencyResult(null);
    expect(useScfStore.getState().frequencyState.result).toBeNull();
  });

  it('setFrequencyIsComputing(true)', () => {
    useScfStore.getState().setFrequencyIsComputing(true);
    expect(useScfStore.getState().frequencyState.isComputing).toBe(true);
  });

  it('setFrequencyIsComputing(false)', () => {
    useScfStore.getState().setFrequencyIsComputing(true);
    useScfStore.getState().setFrequencyIsComputing(false);
    expect(useScfStore.getState().frequencyState.isComputing).toBe(false);
  });

  it('setFrequencyProgress sets progress', () => {
    const mockProgress = { module: 'frequency', phase: 'integrals', percent: 0.5 } as never;
    useScfStore.getState().setFrequencyProgress(mockProgress);
    expect(useScfStore.getState().frequencyState.progress).toBe(mockProgress);
  });

  it('setFrequencyProgress(null) clears progress', () => {
    useScfStore.getState().setFrequencyProgress(null);
    expect(useScfStore.getState().frequencyState.progress).toBeNull();
  });

  it('setFrequencyError sets error message', () => {
    useScfStore.getState().setFrequencyError('test error');
    expect(useScfStore.getState().frequencyState.error).toBe('test error');
  });

  it('setFrequencyError(null) clears error', () => {
    useScfStore.getState().setFrequencyError('test error');
    useScfStore.getState().setFrequencyError(null);
    expect(useScfStore.getState().frequencyState.error).toBeNull();
  });

  it('setFrequencySelectedMode(3)', () => {
    useScfStore.getState().setFrequencySelectedMode(3);
    expect(useScfStore.getState().frequencyState.selectedMode).toBe(3);
  });

  it('setFrequencySelectedMode(null)', () => {
    useScfStore.getState().setFrequencySelectedMode(3);
    useScfStore.getState().setFrequencySelectedMode(null);
    expect(useScfStore.getState().frequencyState.selectedMode).toBeNull();
  });

  it('setFrequencyTemperature(500)', () => {
    useScfStore.getState().setFrequencyTemperature(500);
    expect(useScfStore.getState().frequencyState.temperatureK).toBe(500);
  });

  it('setFrequencyPressure(200000)', () => {
    useScfStore.getState().setFrequencyPressure(200000);
    expect(useScfStore.getState().frequencyState.pressurePa).toBe(200000);
  });

  it('setFrequencyDisplayThermo sets thermo', () => {
    const mockThermo = { temperatureK: 300, pressurePa: 101325 } as never;
    useScfStore.getState().setFrequencyDisplayThermo(mockThermo);
    expect(useScfStore.getState().frequencyState.displayThermo).toBe(mockThermo);
  });

  it('setFrequencyDisplayThermo(null) clears thermo', () => {
    useScfStore.getState().setFrequencyDisplayThermo(null);
    expect(useScfStore.getState().frequencyState.displayThermo).toBeNull();
  });

  it('setFrequencySpectrumTab(raman)', () => {
    useScfStore.getState().setFrequencySpectrumTab('raman');
    expect(useScfStore.getState().frequencyState.spectrumTab).toBe('raman');
  });

  it('setFrequencySpectrumTab(ir)', () => {
    useScfStore.getState().setFrequencySpectrumTab('raman');
    useScfStore.getState().setFrequencySpectrumTab('ir');
    expect(useScfStore.getState().frequencyState.spectrumTab).toBe('ir');
  });

  it('setFrequencyBroadeningKind(gaussian)', () => {
    useScfStore.getState().setFrequencyBroadeningKind('gaussian');
    expect(useScfStore.getState().frequencyState.broadeningKind).toBe('gaussian');
  });

  it('setFrequencyBroadeningKind(lorentzian)', () => {
    useScfStore.getState().setFrequencyBroadeningKind('gaussian');
    useScfStore.getState().setFrequencyBroadeningKind('lorentzian');
    expect(useScfStore.getState().frequencyState.broadeningKind).toBe('lorentzian');
  });

  it('setFrequencyFwhmCm1(12)', () => {
    useScfStore.getState().setFrequencyFwhmCm1(12);
    expect(useScfStore.getState().frequencyState.fwhmCm1).toBe(12);
  });

  it('setFrequencyAmplitudeBohr(1.5)', () => {
    useScfStore.getState().setFrequencyAmplitudeBohr(1.5);
    expect(useScfStore.getState().frequencyState.amplitudeBohr).toBe(1.5);
  });

  it('setFrequencyAnimationSpeed(2.0)', () => {
    useScfStore.getState().setFrequencyAnimationSpeed(2.0);
    expect(useScfStore.getState().frequencyState.animationSpeed).toBe(2.0);
  });

  it('setFrequencyIsAnimating(false)', () => {
    useScfStore.getState().setFrequencyIsAnimating(false);
    expect(useScfStore.getState().frequencyState.isAnimating).toBe(false);
  });

  it('setFrequencyShowDisplacementArrows(true)', () => {
    useScfStore.getState().setFrequencyShowDisplacementArrows(true);
    expect(useScfStore.getState().frequencyState.showDisplacementArrows).toBe(true);
  });

  it('setFrequencyUnitsMode(hartree)', () => {
    useScfStore.getState().setFrequencyUnitsMode('hartree');
    expect(useScfStore.getState().frequencyState.unitsMode).toBe('hartree');
  });

  it('setFrequencyUnitsMode(kj_mol)', () => {
    useScfStore.getState().setFrequencyUnitsMode('kj_mol');
    expect(useScfStore.getState().frequencyState.unitsMode).toBe('kj_mol');
  });
});

// ============================================================================
// Reset behavior
// ============================================================================

describe('Reset behavior', () => {
  it('resetFrequencyState restores defaults after modifications', () => {
    const store = useScfStore.getState();
    store.setFrequencySelectedMode(5);
    store.setFrequencyTemperature(500);
    store.setFrequencySpectrumTab('raman');
    store.setFrequencyBroadeningKind('gaussian');
    store.setFrequencyFwhmCm1(15);
    store.setFrequencyAmplitudeBohr(1.5);
    store.setFrequencyAnimationSpeed(2.0);
    store.setFrequencyIsAnimating(false);
    store.setFrequencyShowDisplacementArrows(true);
    store.setFrequencyUnitsMode('hartree');

    useScfStore.getState().resetFrequencyState();

    expect(useScfStore.getState().frequencyState).toEqual(DEFAULT_FREQUENCY_STATE);
  });

  it('global reset() also resets frequency state', () => {
    useScfStore.getState().setFrequencySelectedMode(7);
    useScfStore.getState().setFrequencyTemperature(400);

    useScfStore.getState().reset();

    expect(useScfStore.getState().frequencyState).toEqual(DEFAULT_FREQUENCY_STATE);
  });

  it('resetFrequencyState does NOT affect other store slices', () => {
    // Modify optimization state
    useScfStore.getState().setTrajectoryStep(5);
    // Modify frequency state
    useScfStore.getState().setFrequencySelectedMode(3);

    // Reset only frequency
    useScfStore.getState().resetFrequencyState();

    // Frequency should be reset
    expect(useScfStore.getState().frequencyState).toEqual(DEFAULT_FREQUENCY_STATE);
    // Optimization should remain modified
    expect(useScfStore.getState().optimizationState.trajectoryStep).toBe(5);
  });
});

// ============================================================================
// Tab switch persistence (AC4)
// ============================================================================

describe('Tab switch persistence', () => {
  it('result persists after simulated tab switch', () => {
    const mockResult = { frequenciesCm1: [100, 200, 300] } as never;
    useScfStore.getState().setFrequencyResult(mockResult);

    // "Unmount" tab — store state persists because it is global
    // Re-read store state (simulates re-mount)
    expect(useScfStore.getState().frequencyState.result).toBe(mockResult);
  });

  it('selectedMode persists after simulated tab switch', () => {
    useScfStore.getState().setFrequencySelectedMode(5);
    expect(useScfStore.getState().frequencyState.selectedMode).toBe(5);
  });

  it('temperature persists after simulated tab switch', () => {
    useScfStore.getState().setFrequencyTemperature(500);
    expect(useScfStore.getState().frequencyState.temperatureK).toBe(500);
  });

  it('spectrum settings persist after simulated tab switch', () => {
    useScfStore.getState().setFrequencySpectrumTab('raman');
    useScfStore.getState().setFrequencyBroadeningKind('gaussian');
    useScfStore.getState().setFrequencyFwhmCm1(12);

    expect(useScfStore.getState().frequencyState.spectrumTab).toBe('raman');
    expect(useScfStore.getState().frequencyState.broadeningKind).toBe('gaussian');
    expect(useScfStore.getState().frequencyState.fwhmCm1).toBe(12);
  });

  it('units mode persists after simulated tab switch', () => {
    useScfStore.getState().setFrequencyUnitsMode('hartree');
    expect(useScfStore.getState().frequencyState.unitsMode).toBe('hartree');
  });

  it('viewer settings persist after simulated tab switch', () => {
    useScfStore.getState().setFrequencyAmplitudeBohr(1.5);
    useScfStore.getState().setFrequencyAnimationSpeed(2.0);
    useScfStore.getState().setFrequencyIsAnimating(false);
    useScfStore.getState().setFrequencyShowDisplacementArrows(true);

    expect(useScfStore.getState().frequencyState.amplitudeBohr).toBe(1.5);
    expect(useScfStore.getState().frequencyState.animationSpeed).toBe(2.0);
    expect(useScfStore.getState().frequencyState.isAnimating).toBe(false);
    expect(useScfStore.getState().frequencyState.showDisplacementArrows).toBe(true);
  });
});

// ============================================================================
// URL initialization (initializeFromURL)
// ============================================================================

describe('initializeFromURL with frequency params', () => {
  it('decodes freq_mode from URL', () => {
    useScfStore.getState().initializeFromURL({
      system_id: 'h2_sto3g_r1.4',
      conv: 'medium',
      max_iter: 50,
      diis: true,
      freq_mode: 3,
    });
    expect(useScfStore.getState().frequencyState.selectedMode).toBe(3);
  });

  it('decodes freq_temp and freq_broad from URL', () => {
    useScfStore.getState().initializeFromURL({
      system_id: 'h2_sto3g_r1.4',
      conv: 'medium',
      max_iter: 50,
      diis: true,
      freq_temp: 400,
      freq_broad: 'g',
    });
    expect(useScfStore.getState().frequencyState.temperatureK).toBe(400);
    expect(useScfStore.getState().frequencyState.broadeningKind).toBe('gaussian');
  });

  it('defaults frequency state when no freq_* params present', () => {
    useScfStore.getState().initializeFromURL({
      system_id: 'h2_sto3g_r1.4',
      conv: 'medium',
      max_iter: 50,
      diis: true,
    });
    expect(useScfStore.getState().frequencyState).toEqual(DEFAULT_FREQUENCY_STATE);
  });

  it('freq_units=ha maps to hartree', () => {
    useScfStore.getState().initializeFromURL({
      system_id: 'h2_sto3g_r1.4',
      conv: 'medium',
      max_iter: 50,
      diis: true,
      freq_units: 'ha',
    });
    expect(useScfStore.getState().frequencyState.unitsMode).toBe('hartree');
  });

  it('freq_units=kj maps to kj_mol', () => {
    useScfStore.getState().initializeFromURL({
      system_id: 'h2_sto3g_r1.4',
      conv: 'medium',
      max_iter: 50,
      diis: true,
      freq_units: 'kj',
    });
    expect(useScfStore.getState().frequencyState.unitsMode).toBe('kj_mol');
  });

  it('freq params do not affect PES or optimization state', () => {
    useScfStore.getState().initializeFromURL({
      system_id: 'h2_sto3g_r1.4',
      conv: 'medium',
      max_iter: 50,
      diis: true,
      freq_mode: 5,
      freq_temp: 400,
    });
    expect(useScfStore.getState().pesState.scanning).toBe(false);
    expect(useScfStore.getState().optimizationState.running).toBe(false);
  });
});
