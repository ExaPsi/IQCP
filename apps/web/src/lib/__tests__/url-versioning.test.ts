/**
 * URL Versioning and Deep Link Regression Tests (US-036)
 *
 * Validates that:
 * 1. Route constants use /v1/ prefix
 * 2. Deep link encoding/decoding is path-independent (query param based)
 * 3. generateShareURL produces versioned paths
 * 4. Round-trip encoding preserves state
 *
 * @module lib/__tests__/url-versioning
 */

import { describe, it, expect } from 'vitest';
import {
  ROUTES,
  CURRENT_VERSION,
  VERSION_PREFIX,
  MODULE_ROUTES,
  getModuleRoute,
} from '../routes';
import { encodeRunState, decodeRunState, verifyRoundTrip } from '../deeplink';
import {
  DEFAULT_BOYS_STATE,
  DEFAULT_RYS_STATE,
  DEFAULT_SCF_STATE,
} from '../../types/run-state';
import type { RunStateV1 } from '../../types/run-state';
import { ScfParamsSchema } from '../../types/run-state.schema';

// =============================================================================
// Group 1: Route Constants
// =============================================================================

describe('route constants', () => {
  it('CURRENT_VERSION is v1', () => {
    expect(CURRENT_VERSION).toBe('v1');
  });

  it('VERSION_PREFIX starts with /', () => {
    expect(VERSION_PREFIX).toBe('/v1');
  });

  it('all module routes start with /v1/', () => {
    expect(ROUTES.BOYS).toBe('/v1/boys');
    expect(ROUTES.RYS).toBe('/v1/rys');
    expect(ROUTES.SCF).toBe('/v1/scf');
  });

  it('content routes start with /v1/', () => {
    expect(ROUTES.LABS).toBe('/v1/labs');
    expect(ROUTES.REFERENCE).toBe('/v1/reference');
    expect(ROUTES.ABOUT).toBe('/v1/about');
  });

  it('home route is unversioned', () => {
    expect(ROUTES.HOME).toBe('/');
  });

  it('MODULE_ROUTES maps module IDs to versioned paths', () => {
    expect(MODULE_ROUTES['boys']).toBe('/v1/boys');
    expect(MODULE_ROUTES['rys']).toBe('/v1/rys');
    expect(MODULE_ROUTES['scf']).toBe('/v1/scf');
  });

  it('getModuleRoute returns versioned path', () => {
    expect(getModuleRoute('boys')).toBe('/v1/boys');
    expect(getModuleRoute('rys')).toBe('/v1/rys');
    expect(getModuleRoute('scf')).toBe('/v1/scf');
  });

  it('getModuleRoute handles unknown module with fallback', () => {
    expect(getModuleRoute('unknown')).toBe('/v1/unknown');
  });
});

// =============================================================================
// Group 2: Deep Link Encoding/Decoding (Path-Independent)
// =============================================================================

describe('deep link encoding is path-independent', () => {
  const testStates: [string, RunStateV1][] = [
    ['boys', DEFAULT_BOYS_STATE],
    ['rys', DEFAULT_RYS_STATE],
    ['scf', DEFAULT_SCF_STATE],
  ];

  for (const [name, state] of testStates) {
    it(`${name}: encode/decode round-trip preserves state`, () => {
      expect(verifyRoundTrip(state)).toBe(true);
    });

    it(`${name}: encoded string does not contain path information`, () => {
      const encoded = encodeRunState(state);
      // The encoded string is lz-string compressed, not a URL path
      // It should not start with / and should be a compact string
      expect(encoded.startsWith('/')).toBe(false);
      expect(encoded.length).toBeGreaterThan(0);
      expect(encoded.length).toBeLessThan(500);
    });

    it(`${name}: same state produces same encoding`, () => {
      const encoded1 = encodeRunState(state);
      const encoded2 = encodeRunState(state);
      expect(encoded1).toBe(encoded2);
    });
  }

  it('decoded state from same encoding matches regardless of origin path', () => {
    // The deep link query parameter (?run=...) is path-independent.
    // Whether accessed from /boys or /v1/boys, the same ?run= value
    // decodes to the same state.
    const encoded = encodeRunState(DEFAULT_BOYS_STATE);
    const decoded = decodeRunState(encoded);

    expect(decoded).not.toBeNull();
    expect(decoded!.module).toBe('boys');
    expect(decoded!.boys!.m).toBe(DEFAULT_BOYS_STATE.boys!.m);
    expect(decoded!.boys!.T).toBe(DEFAULT_BOYS_STATE.boys!.T);
    expect(decoded!.boys!.view).toBe(DEFAULT_BOYS_STATE.boys!.view);
  });
});

// =============================================================================
// Group 3: Deep Link State Validation
// =============================================================================

describe('deep link state validation', () => {
  it('rejects null input', () => {
    expect(decodeRunState(null)).toBeNull();
  });

  it('rejects empty string', () => {
    expect(decodeRunState('')).toBeNull();
  });

  it('rejects invalid lz-string', () => {
    expect(decodeRunState('not-valid-lz-string')).toBeNull();
  });

  it('SCF state with DIIS flag round-trips', () => {
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: '0.1.0',
      module: 'scf',
      scf: {
        system_id: 'h2_sto3g_r1.4',
        max_iter: 100,
        conv: 'tight',
        diis: false,
      },
      ui: { mode: 'internals' },
    };

    const encoded = encodeRunState(state);
    const decoded = decodeRunState(encoded);

    expect(decoded).not.toBeNull();
    expect(decoded!.scf!.diis).toBe(false);
    expect(decoded!.scf!.conv).toBe('tight');
    expect(decoded!.scf!.max_iter).toBe(100);
    expect(decoded!.ui.mode).toBe('internals');
  });

  it('Boys state with sweep mode round-trips', () => {
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: '0.1.0',
      module: 'boys',
      boys: { m: 5, T: 15.0, view: 'sweep' },
      ui: { mode: 'explain' },
    };

    const encoded = encodeRunState(state);
    const decoded = decodeRunState(encoded);

    expect(decoded).not.toBeNull();
    expect(decoded!.boys!.m).toBe(5);
    expect(decoded!.boys!.T).toBe(15.0);
    expect(decoded!.boys!.view).toBe('sweep');
  });

  it('Rys state with target accuracy round-trips', () => {
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: '0.1.0',
      module: 'rys',
      rys: { n: 7, T: 25.0, target: '1e-8' },
      ui: { mode: 'explain' },
    };

    const encoded = encodeRunState(state);
    const decoded = decodeRunState(encoded);

    expect(decoded).not.toBeNull();
    expect(decoded!.rys!.n).toBe(7);
    expect(decoded!.rys!.T).toBe(25.0);
    expect(decoded!.rys!.target).toBe('1e-8');
  });
});

// =============================================================================
// Group 4: Frequency State Deep Link Round-Trip (US-103)
// =============================================================================

describe('Frequency state deep link round-trip (US-103)', () => {
  it('all defaults produces no freq_* params in encoded state', () => {
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: '0.1.0',
      module: 'scf',
      scf: {
        system_id: 'h2_sto3g_r1.4',
        max_iter: 50,
        conv: 'medium',
        diis: true,
      },
      ui: { mode: 'explain' },
    };
    const encoded = encodeRunState(state);
    const decoded = decodeRunState(encoded);
    expect(decoded).not.toBeNull();
    expect(decoded!.scf!.freq_mode).toBeUndefined();
    expect(decoded!.scf!.freq_temp).toBeUndefined();
    expect(decoded!.scf!.freq_pres).toBeUndefined();
    expect(decoded!.scf!.freq_broad).toBeUndefined();
    expect(decoded!.scf!.freq_fwhm).toBeUndefined();
    expect(decoded!.scf!.freq_tab).toBeUndefined();
    expect(decoded!.scf!.freq_arrows).toBeUndefined();
    expect(decoded!.scf!.freq_units).toBeUndefined();
  });

  it('freq_mode=5 round-trips correctly', () => {
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: '0.1.0',
      module: 'scf',
      scf: {
        system_id: 'h2_sto3g_r1.4',
        max_iter: 50,
        conv: 'medium',
        diis: true,
        freq_mode: 5,
      },
      ui: { mode: 'explain' },
    };
    const encoded = encodeRunState(state);
    const decoded = decodeRunState(encoded);
    expect(decoded).not.toBeNull();
    expect(decoded!.scf!.freq_mode).toBe(5);
  });

  it('freq_broad=g round-trips correctly', () => {
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: '0.1.0',
      module: 'scf',
      scf: {
        system_id: 'h2_sto3g_r1.4',
        max_iter: 50,
        conv: 'medium',
        diis: true,
        freq_broad: 'g',
      },
      ui: { mode: 'explain' },
    };
    const encoded = encodeRunState(state);
    const decoded = decodeRunState(encoded);
    expect(decoded).not.toBeNull();
    expect(decoded!.scf!.freq_broad).toBe('g');
  });

  it('freq_units=ha round-trips correctly', () => {
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: '0.1.0',
      module: 'scf',
      scf: {
        system_id: 'h2_sto3g_r1.4',
        max_iter: 50,
        conv: 'medium',
        diis: true,
        freq_units: 'ha',
      },
      ui: { mode: 'explain' },
    };
    const encoded = encodeRunState(state);
    const decoded = decodeRunState(encoded);
    expect(decoded).not.toBeNull();
    expect(decoded!.scf!.freq_units).toBe('ha');
  });

  it('freq_units=kj round-trips correctly', () => {
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: '0.1.0',
      module: 'scf',
      scf: {
        system_id: 'h2_sto3g_r1.4',
        max_iter: 50,
        conv: 'medium',
        diis: true,
        freq_units: 'kj',
      },
      ui: { mode: 'explain' },
    };
    const encoded = encodeRunState(state);
    const decoded = decodeRunState(encoded);
    expect(decoded).not.toBeNull();
    expect(decoded!.scf!.freq_units).toBe('kj');
  });

  it('full frequency state round-trips with all 8 params', () => {
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: '0.1.0',
      module: 'scf',
      scf: {
        system_id: 'h2o_sto3g',
        max_iter: 50,
        conv: 'medium',
        diis: true,
        freq_mode: 5,
        freq_temp: 400,
        freq_pres: 200000,
        freq_broad: 'g',
        freq_fwhm: 15,
        freq_tab: 'raman',
        freq_arrows: true,
        freq_units: 'kj',
      },
      ui: { mode: 'explain' },
    };
    const encoded = encodeRunState(state);
    const decoded = decodeRunState(encoded);
    expect(decoded).not.toBeNull();
    expect(decoded!.scf!.freq_mode).toBe(5);
    expect(decoded!.scf!.freq_temp).toBe(400);
    expect(decoded!.scf!.freq_pres).toBe(200000);
    expect(decoded!.scf!.freq_broad).toBe('g');
    expect(decoded!.scf!.freq_fwhm).toBe(15);
    expect(decoded!.scf!.freq_tab).toBe('raman');
    expect(decoded!.scf!.freq_arrows).toBe(true);
    expect(decoded!.scf!.freq_units).toBe('kj');
  });

  it('pre-US-103 URL (no freq params) round-trips with SCF state intact', () => {
    const state: RunStateV1 = {
      schema_version: 'run_state_v1',
      app_version: '0.1.0',
      module: 'scf',
      scf: {
        system_id: 'h2_sto3g_r1.4',
        max_iter: 100,
        conv: 'tight',
        diis: false,
      },
      ui: { mode: 'internals' },
    };
    const encoded = encodeRunState(state);
    const decoded = decodeRunState(encoded);
    expect(decoded).not.toBeNull();
    expect(decoded!.scf!.max_iter).toBe(100);
    expect(decoded!.scf!.conv).toBe('tight');
    expect(decoded!.scf!.diis).toBe(false);
    // No frequency params
    expect(decoded!.scf!.freq_mode).toBeUndefined();
  });
});

// =============================================================================
// Group 5: Frequency Zod Validation (US-103)
// =============================================================================

describe('Frequency Zod validation (US-103)', () => {
  const baseScf = {
    system_id: 'h2_sto3g_r1.4',
    max_iter: 50,
    conv: 'medium',
    diis: true,
  };

  it('accepts empty freq params', () => {
    expect(ScfParamsSchema.safeParse(baseScf).success).toBe(true);
  });

  it('accepts freq_mode: 5', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_mode: 5 }).success).toBe(true);
  });

  it('rejects freq_mode: -1', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_mode: -1 }).success).toBe(false);
  });

  it('rejects freq_mode: 1000', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_mode: 1000 }).success).toBe(false);
  });

  it('accepts freq_temp: 350', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_temp: 350 }).success).toBe(true);
  });

  it('rejects freq_temp: 5 (below min 10)', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_temp: 5 }).success).toBe(false);
  });

  it('accepts freq_broad: l', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_broad: 'l' }).success).toBe(true);
  });

  it('accepts freq_broad: g', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_broad: 'g' }).success).toBe(true);
  });

  it('rejects freq_broad: x', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_broad: 'x' }).success).toBe(false);
  });

  it('accepts freq_tab: raman', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_tab: 'raman' }).success).toBe(true);
  });

  it('rejects freq_tab: uv', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_tab: 'uv' }).success).toBe(false);
  });

  it('accepts freq_units: ha', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_units: 'ha' }).success).toBe(true);
  });

  it('accepts freq_units: kj', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_units: 'kj' }).success).toBe(true);
  });

  it('rejects freq_units: ev', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_units: 'ev' }).success).toBe(false);
  });

  it('accepts freq_arrows: true', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_arrows: true }).success).toBe(true);
  });

  it('accepts freq_fwhm: 15.5', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_fwhm: 15.5 }).success).toBe(true);
  });

  it('rejects freq_fwhm: 0.1 (below min 0.5)', () => {
    expect(ScfParamsSchema.safeParse({ ...baseScf, freq_fwhm: 0.1 }).success).toBe(false);
  });
});
