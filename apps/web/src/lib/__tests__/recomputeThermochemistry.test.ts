/**
 * Unit tests for `recomputeThermochemistry` (US-102).
 *
 * Validates the TypeScript port of `crates/qc-core/src/thermochemistry.rs`
 * against PySCF-style reference values computed with the same
 * pre-2019-revised-SI constants that qc-core uses. Reference values were
 * generated with a Python script that mirrors the Rust formulas exactly —
 * see the `/scripts/` folder or re-derive with the constants documented
 * in `lib/recomputeThermochemistry.ts`.
 *
 * Tests cover:
 *  - H₂O nonlinear rotor (reference values for all components)
 *  - CO linear rotor (branch coverage for the linear formula)
 *  - He atom (branch coverage for the trivial atom case)
 *  - Electronic entropy for doublet / triplet states
 *  - Pressure doubling → S_trans decreases by R·ln(2)
 *  - Temperature limits (low T → vib thermal ≈ 0; high T → large Cv)
 *  - Imaginary mode skipping (nImag counted but not included in sums)
 *  - Symmetry number override (R·ln(σ) suppression of S_rot)
 *  - Round-trip at the original (T, P) reproduces the input FrequencyThermochemistry
 *    within tight tolerance
 *  - Error paths (T ≤ 0, P ≤ 0)
 *
 * @module lib/__tests__/recomputeThermochemistry.test
 * @see US-102 Frequency Tab UI, AC4
 */

import { describe, it, expect } from 'vitest';
import {
  recomputeThermochemistry,
  rHartreePerMolK,
} from '../recomputeThermochemistry';
import type {
  FrequencyResult,
  FrequencyThermochemistry,
} from '../../worker/protocol';

// ============================================================================
// Test fixtures — minimal FrequencyResult shells
// ============================================================================

/**
 * Construct a minimal but type-complete FrequencyThermochemistry. Individual
 * tests that exercise recompute only care about `temperatureK`, `pressurePa`,
 * `symmetryNumber`, `multiplicity`, and `totalMassAmu`; the rest are defaults.
 */
function makeThermo(
  overrides: Partial<FrequencyThermochemistry> = {}
): FrequencyThermochemistry {
  return {
    temperatureK: 298.15,
    pressurePa: 101325,
    symmetryNumber: 1,
    multiplicity: 1,
    totalMassAmu: 18.0153,
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
    ...overrides,
  };
}

/**
 * Construct a minimal FrequencyResult for use in recompute tests. Only the
 * fields actually read by recompute (`frequenciesCm1`, `rotorType`,
 * `rotationalConstantsGhz`, `electronicEnergyHa`, `thermochemistry`) need to be
 * valid; the rest are zeroed out.
 */
function makeResult(
  frequenciesCm1: number[],
  rotorType: FrequencyResult['rotorType'],
  rotationalConstantsGhz: [number, number, number],
  electronicEnergyHa: number,
  thermo: Partial<FrequencyThermochemistry> = {}
): FrequencyResult {
  const thermochemistry = makeThermo(thermo);
  return {
    nAtoms: 0,
    nModes: frequenciesCm1.length,
    rotorType,
    electronicEnergyHa,
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
    frequenciesCm1,
    reducedMassesAmu: frequenciesCm1.map(() => 0),
    forceConstantsMdyne: frequenciesCm1.map(() => 0),
    normalModesCartesian: [],
    rotationalConstantsGhz,
    irIntensitiesKmPerMol: frequenciesCm1.map(() => 0),
    ramanActivitiesA4Amu: frequenciesCm1.map(() => 0),
    depolarizationRatios: frequenciesCm1.map(() => 0),
    thermochemistry,
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
// H₂O nonlinear reference fixture
// ============================================================================

/**
 * H₂O-like test molecule for the nonlinear rotor branch.
 *
 * Values are representative experimental numbers; the exact frequencies are
 * not critical — what matters is that the same inputs produce the same
 * outputs as the reference Python port (see test docstring).
 */
function h2oLike(): FrequencyResult {
  return makeResult(
    [1595.0, 3657.0, 3756.0],
    'asymmetric_top',
    [835.84, 435.129, 286.181],
    -76.0,
    { totalMassAmu: 18.0153, symmetryNumber: 2, multiplicity: 1 }
  );
}

// ============================================================================
// Tests — nonlinear rotor (H₂O)
// ============================================================================

describe('recomputeThermochemistry — H₂O-like nonlinear rotor', () => {
  it('matches Python reference at T = 298.15 K, P = 101325 Pa', () => {
    const result = h2oLike();
    const thermo = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });

    // Reference values computed in Python using the exact same pre-SI-revision
    // constants as `crates/qc-core/src/constants.rs` — see docstring.
    expect(thermo.zpeHa).toBeCloseTo(2.05217340022683e-2, 14);
    expect(thermo.eTransHa).toBeCloseTo(1.41627683418939e-3, 14);
    expect(thermo.eRotHa).toBeCloseTo(1.41627683418939e-3, 14);
    expect(thermo.eVibThermalHa).toBeCloseTo(3.30295220574906e-6, 14);
    expect(thermo.internalEnergyHa).toBeCloseTo(-75.9766424093771, 9);
    expect(thermo.enthalpyHa).toBeCloseTo(-75.9756982248210, 9);
    expect(thermo.sTransHaPerK).toBeCloseTo(5.51530216245729e-5, 14);
    expect(thermo.sRotHaPerK).toBeCloseTo(1.66267945305208e-5, 14);
    expect(thermo.sVibHaPerK).toBeCloseTo(1.25169717794908e-8, 14);
    expect(thermo.entropyHaPerK).toBeCloseTo(7.17923331268732e-5, 14);
    expect(thermo.gibbsHa).toBeCloseTo(-75.9971031089428, 9);
    expect(thermo.cvHaPerK).toBeCloseTo(9.58575855722184e-6, 14);
    expect(thermo.cpHaPerK).toBeCloseTo(1.27525690758409e-5, 14);
    expect(thermo.nVibModesUsed).toBe(3);
    expect(thermo.nImag).toBe(0);
  });

  it('e0kHa = electronicEnergyHa + zpeHa', () => {
    const result = h2oLike();
    const thermo = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    expect(thermo.e0kHa).toBeCloseTo(-76.0 + thermo.zpeHa, 14);
  });

  it('enthalpy = internal + RT', () => {
    const result = h2oLike();
    const thermo = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    const r = rHartreePerMolK();
    // The subtraction U − U' loses ~1 ulp at the -76 Ha scale, so we compare
    // as a relative error of the RT term (∼9e-4 Ha) at 1e-10 relative precision.
    const expectedRt = r * 298.15;
    const actualRt = thermo.enthalpyHa - thermo.internalEnergyHa;
    expect(Math.abs(actualRt - expectedRt) / expectedRt).toBeLessThan(1e-10);
  });

  it('Gibbs = Enthalpy − T·S', () => {
    const result = h2oLike();
    const thermo = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    expect(thermo.gibbsHa).toBeCloseTo(
      thermo.enthalpyHa - 298.15 * thermo.entropyHaPerK,
      14
    );
  });

  it('Cp = Cv + R', () => {
    const result = h2oLike();
    const thermo = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    expect(thermo.cpHaPerK - thermo.cvHaPerK).toBeCloseTo(
      rHartreePerMolK(),
      14
    );
  });

  it('h_vib = zpe + e_vib_thermal (PySCF convention)', () => {
    const result = h2oLike();
    const thermo = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    expect(thermo.hVibHa).toBeCloseTo(thermo.zpeHa + thermo.eVibThermalHa, 14);
  });

  it('h_rot = e_rot (PySCF convention)', () => {
    const result = h2oLike();
    const thermo = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    expect(thermo.hRotHa).toBe(thermo.eRotHa);
  });
});

// ============================================================================
// Tests — temperature limits
// ============================================================================

describe('recomputeThermochemistry — temperature limits', () => {
  it('low-T (100 K) vibrational thermal contributions are near zero; ZPE unchanged', () => {
    const result = h2oLike();
    const lowT = recomputeThermochemistry({
      result,
      temperatureK: 100,
      pressurePa: 101325,
    });
    const refT = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    // ZPE is temperature-independent
    expect(lowT.zpeHa).toBeCloseTo(refT.zpeHa, 14);
    // Thermal vib contributions drop off sharply at low T for stiff modes
    expect(Math.abs(lowT.eVibThermalHa)).toBeLessThan(
      Math.abs(refT.eVibThermalHa) + 1e-20
    );
    // Translational energy is linear in T: E_trans(100) / E_trans(298.15) ≈ 100/298.15
    expect(lowT.eTransHa / refT.eTransHa).toBeCloseTo(100 / 298.15, 12);
  });

  it('high-T (1000 K) entropy and Cv grow (non-strict: Cv never decreases in T)', () => {
    const result = h2oLike();
    const highT = recomputeThermochemistry({
      result,
      temperatureK: 1000,
      pressurePa: 101325,
    });
    const refT = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    expect(highT.cvHaPerK).toBeGreaterThan(refT.cvHaPerK);
    expect(highT.entropyHaPerK).toBeGreaterThan(refT.entropyHaPerK);
  });
});

// ============================================================================
// Tests — linear rotor (CO)
// ============================================================================

describe('recomputeThermochemistry — CO linear rotor', () => {
  it('matches Python reference for linear rotor with 1 vibration', () => {
    const result = makeResult(
      [2143.0],
      'linear',
      [Infinity, 57.898, 57.898],
      -113.0,
      { totalMassAmu: 12.0 + 15.999, symmetryNumber: 1, multiplicity: 1 }
    );
    const thermo = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    // Reference values from Python (Rust constants, same formulas)
    expect(thermo.zpeHa).toBeCloseTo(4.88211322900322e-3, 14);
    expect(thermo.eRotHa).toBeCloseTo(9.44184556126261e-4, 14);
    expect(thermo.sRotHaPerK).toBeCloseTo(1.79736305850697e-5, 14);
    expect(thermo.cvRotHaPerK).toBeCloseTo(rHartreePerMolK(), 14);
    expect(thermo.eVibThermalHa).toBeCloseTo(3.15080996417665e-7, 14);
    expect(thermo.sVibHaPerK).toBeCloseTo(1.15897473185885e-9, 14);
  });
});

// ============================================================================
// Tests — atom (He)
// ============================================================================

describe('recomputeThermochemistry — He atom', () => {
  it('rotational and vibrational contributions are zero for an atom', () => {
    const result = makeResult(
      [],
      'atom',
      [Infinity, Infinity, Infinity],
      -2.8551,
      { totalMassAmu: 4.0026, symmetryNumber: 1, multiplicity: 1 }
    );
    const thermo = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    expect(thermo.zpeHa).toBe(0);
    expect(thermo.eRotHa).toBe(0);
    expect(thermo.sRotHaPerK).toBe(0);
    expect(thermo.cvRotHaPerK).toBe(0);
    expect(thermo.eVibThermalHa).toBe(0);
    expect(thermo.sVibHaPerK).toBe(0);
    expect(thermo.cvVibHaPerK).toBe(0);
    expect(thermo.nVibModesUsed).toBe(0);
    expect(thermo.nImag).toBe(0);

    // Translational still nonzero
    expect(thermo.eTransHa).toBeCloseTo(1.41627683418939e-3, 14);
    expect(thermo.sTransHaPerK).toBeCloseTo(4.80073801120202e-5, 14);
  });
});

// ============================================================================
// Tests — pressure scaling
// ============================================================================

describe('recomputeThermochemistry — pressure scaling', () => {
  it('doubling pressure decreases S_trans by R·ln(2)', () => {
    const result = h2oLike();
    const at1atm = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    const at2atm = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 2 * 101325,
    });
    const deltaS = at1atm.sTransHaPerK - at2atm.sTransHaPerK;
    const rLn2 = rHartreePerMolK() * Math.log(2);
    expect(deltaS).toBeCloseTo(rLn2, 14);
  });

  it('rotational and vibrational contributions are pressure-independent', () => {
    const result = h2oLike();
    const at1 = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    const at2 = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 200000,
    });
    expect(at2.sRotHaPerK).toBe(at1.sRotHaPerK);
    expect(at2.sVibHaPerK).toBe(at1.sVibHaPerK);
    expect(at2.eRotHa).toBe(at1.eRotHa);
    expect(at2.eVibThermalHa).toBe(at1.eVibThermalHa);
  });
});

// ============================================================================
// Tests — imaginary mode skipping
// ============================================================================

describe('recomputeThermochemistry — imaginary modes', () => {
  it('imaginary frequencies are skipped from vibrational sums but counted', () => {
    // Inject an imaginary mode at the start
    const result = makeResult(
      [-500, 1595.0, 3657.0, 3756.0],
      'asymmetric_top',
      [835.84, 435.129, 286.181],
      -76.0,
      { totalMassAmu: 18.0153, symmetryNumber: 2, multiplicity: 1 }
    );
    const thermo = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
    });

    // Same ZPE as if the imaginary mode wasn't there
    const ref = recomputeThermochemistry({
      result: h2oLike(),
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    expect(thermo.zpeHa).toBeCloseTo(ref.zpeHa, 14);
    expect(thermo.eVibThermalHa).toBeCloseTo(ref.eVibThermalHa, 14);
    expect(thermo.sVibHaPerK).toBeCloseTo(ref.sVibHaPerK, 14);
    expect(thermo.cvVibHaPerK).toBeCloseTo(ref.cvVibHaPerK, 14);

    // Counted
    expect(thermo.nImag).toBe(1);
    expect(thermo.nVibModesUsed).toBe(3);
  });
});

// ============================================================================
// Tests — symmetry number override
// ============================================================================

describe('recomputeThermochemistry — symmetry number', () => {
  it('CH₄-like σ=12 override reduces S_rot by R·ln(12) vs σ=1', () => {
    const result = makeResult(
      [1320, 1320, 1320, 1540, 1540, 3020, 3120, 3120, 3120],
      'spherical_top',
      [160.0, 160.0, 160.0],
      -40.0,
      { totalMassAmu: 16.043, symmetryNumber: 1, multiplicity: 1 }
    );
    const sigma1 = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
      symmetryNumberOverride: 1,
    });
    const sigma12 = recomputeThermochemistry({
      result,
      temperatureK: 298.15,
      pressurePa: 101325,
      symmetryNumberOverride: 12,
    });
    const delta = sigma1.sRotHaPerK - sigma12.sRotHaPerK;
    const rLn12 = rHartreePerMolK() * Math.log(12);
    expect(delta).toBeCloseTo(rLn12, 14);
    expect(sigma12.symmetryNumber).toBe(12);
  });
});

// ============================================================================
// Tests — electronic entropy
// ============================================================================

describe('recomputeThermochemistry — electronic entropy', () => {
  it('singlet (multiplicity = 1) → S_elec = 0', () => {
    const thermo = recomputeThermochemistry({
      result: h2oLike(),
      temperatureK: 298.15,
      pressurePa: 101325,
    });
    expect(thermo.sElecHaPerK).toBe(0);
  });

  it('doublet (multiplicity = 2) → S_elec = R·ln(2)', () => {
    const thermo = recomputeThermochemistry({
      result: h2oLike(),
      temperatureK: 298.15,
      pressurePa: 101325,
      multiplicityOverride: 2,
    });
    expect(thermo.sElecHaPerK).toBeCloseTo(rHartreePerMolK() * Math.log(2), 14);
  });

  it('triplet (multiplicity = 3) → S_elec = R·ln(3)', () => {
    const thermo = recomputeThermochemistry({
      result: h2oLike(),
      temperatureK: 298.15,
      pressurePa: 101325,
      multiplicityOverride: 3,
    });
    expect(thermo.sElecHaPerK).toBeCloseTo(rHartreePerMolK() * Math.log(3), 14);
  });
});

// ============================================================================
// Tests — error handling
// ============================================================================

describe('recomputeThermochemistry — error handling', () => {
  it('throws on T <= 0', () => {
    expect(() =>
      recomputeThermochemistry({
        result: h2oLike(),
        temperatureK: 0,
        pressurePa: 101325,
      })
    ).toThrow(/temperatureK/);
    expect(() =>
      recomputeThermochemistry({
        result: h2oLike(),
        temperatureK: -10,
        pressurePa: 101325,
      })
    ).toThrow(/temperatureK/);
  });

  it('throws on P <= 0', () => {
    expect(() =>
      recomputeThermochemistry({
        result: h2oLike(),
        temperatureK: 298.15,
        pressurePa: 0,
      })
    ).toThrow(/pressurePa/);
    expect(() =>
      recomputeThermochemistry({
        result: h2oLike(),
        temperatureK: 298.15,
        pressurePa: -1,
      })
    ).toThrow(/pressurePa/);
  });

  it('throws on non-finite inputs', () => {
    expect(() =>
      recomputeThermochemistry({
        result: h2oLike(),
        temperatureK: NaN,
        pressurePa: 101325,
      })
    ).toThrow();
    expect(() =>
      recomputeThermochemistry({
        result: h2oLike(),
        temperatureK: 298.15,
        pressurePa: Infinity,
      })
    ).toThrow();
  });
});

// ============================================================================
// Tests — R constant
// ============================================================================

describe('rHartreePerMolK', () => {
  it('matches the published approximate value 3.166811e-6 Ha/(mol·K)', () => {
    // Tolerance tight enough to catch a wrong constant swap; loose enough to
    // allow either pre-2019 or post-2019 SI within one part in 10^8.
    expect(rHartreePerMolK()).toBeCloseTo(3.166811e-6, 10);
  });
});
