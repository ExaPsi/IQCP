/**
 * TypeScript port of the `crates/qc-core/src/thermochemistry.rs` RRHO
 * thermochemistry routine for client-side recomputation in the Frequency tab.
 *
 * When the user drags the temperature slider, re-running the full WASM
 * `compute_frequencies` pipeline would take seconds. Instead, the frequencies,
 * rotational constants, total mass, and electronic energy from a prior WASM
 * result are reused, and only the (cheap) RRHO partition functions are
 * evaluated in TypeScript.
 *
 * This implementation is a faithful line-by-line port of the Rust reference
 * at `crates/qc-core/src/thermochemistry.rs:373-776`. Physical constants are
 * matched to the Rust values in `crates/qc-core/src/constants.rs` (NOT the
 * CODATA 2018 values) so that the output agrees bit-for-bit with the
 * WASM-produced `FrequencyResult.thermochemistry` for the call's native T and P.
 *
 * @module lib/recomputeThermochemistry
 * @see crates/qc-core/src/thermochemistry.rs — authoritative Rust source
 * @see US-099 RRHO Thermochemistry (original Rust implementation)
 * @see US-102 Frequency Tab UI
 */

import type {
  FrequencyResult,
  FrequencyThermochemistry,
  RotorType,
} from '../worker/protocol';

// ============================================================================
// Physical constants — matched to `crates/qc-core/src/constants.rs`
// ============================================================================

/**
 * Boltzmann constant in J/K.
 * Matches `qc-core::constants::BOLTZMANN` (pre-2019 revised SI value).
 */
const BOLTZMANN = 1.38064852e-23;

/**
 * Planck constant in J·s.
 * Matches `qc-core::constants::PLANCK` (pre-2019 revised SI value).
 */
const PLANCK = 6.62607004e-34;

/**
 * Avogadro constant (per mol).
 * Matches `qc-core::constants::AVOGADRO` (pre-2019 revised SI value).
 */
const AVOGADRO = 6.022140857e23;

/**
 * Hartree → Joule conversion.
 * Matches `qc-core::constants::HARTREE_TO_JOULE`.
 */
const HARTREE_TO_JOULE = 4.359744644911914e-18;

/**
 * Atomic mass unit → kilogram conversion.
 * Matches `qc-core::constants::AMU_TO_KG = 1e-3 / AVOGADRO`.
 */
const AMU_TO_KG = 1e-3 / AVOGADRO;

/**
 * Speed of light in m/s (exact in SI).
 * Matches `qc-core::constants::LIGHT_SPEED_SI`.
 */
const LIGHT_SPEED_SI = 299792458.0;

const PI = Math.PI;

/**
 * Gas constant R in Hartree per mole per Kelvin.
 *
 * `R_Eh = kB / HARTREE_TO_J` — the Avogadro factor cancels because Hartree is
 * already a molar-scale unit in the quantum chemistry convention. Matches
 * `qc-core::thermochemistry::r_hartree_per_mol_k`.
 *
 * @returns ≈ 3.166811e-6 Ha/(mol·K)
 */
export function rHartreePerMolK(): number {
  return BOLTZMANN / HARTREE_TO_JOULE;
}

// ============================================================================
// Public input shape
// ============================================================================

/**
 * Inputs to `recomputeThermochemistry`.
 *
 * Everything that changes when the user drags the temperature or pressure
 * slider is here; everything else (frequencies, rotor type, rotational
 * constants, electronic energy, total mass) is pulled from the prior
 * `FrequencyResult`.
 */
export interface RecomputeThermochemistryInput {
  /**
   * Original WASM result. Frequencies, rotational constants, rotor type,
   * total mass, and electronic energy are read from here.
   */
  result: FrequencyResult;
  /** New temperature in Kelvin. Must be > 0. */
  temperatureK: number;
  /** New pressure in Pascals. Must be > 0. */
  pressurePa: number;
  /**
   * Override the rotational symmetry number σ.
   * Defaults to `result.thermochemistry.symmetryNumber`.
   */
  symmetryNumberOverride?: number;
  /**
   * Override the spin multiplicity (2S+1).
   * Defaults to `result.thermochemistry.multiplicity`.
   */
  multiplicityOverride?: number;
}

// ============================================================================
// Internal: translational partition
// ============================================================================

/**
 * Translational partition function contributions (Sackur-Tetrode).
 *
 * Rust reference: `thermochemistry.rs:373-391`.
 */
function translationalPartition(
  totalMassAmu: number,
  temperatureK: number,
  pressurePa: number
): {
  eTransHa: number;
  hTransHa: number;
  sTransHaPerK: number;
  cvTransHaPerK: number;
  cpTransHaPerK: number;
} {
  const r = rHartreePerMolK();
  const massKg = totalMassAmu * AMU_TO_KG;

  // q_trans = (2π m kB T / h²)^(3/2) · kB T / P
  const factor = (2 * PI * massKg * BOLTZMANN * temperatureK) / (PLANCK * PLANCK);
  const qTrans =
    Math.pow(factor, 1.5) * ((BOLTZMANN * temperatureK) / pressurePa);

  const eTransHa = 1.5 * r * temperatureK;
  const hTransHa = 2.5 * r * temperatureK;
  const sTransHaPerK = r * (Math.log(qTrans) + 2.5);
  const cvTransHaPerK = 1.5 * r;
  const cpTransHaPerK = 2.5 * r;

  return { eTransHa, hTransHa, sTransHaPerK, cvTransHaPerK, cpTransHaPerK };
}

// ============================================================================
// Internal: rotational partition
// ============================================================================

/**
 * Rotational partition function contributions.
 *
 * Rust reference: `thermochemistry.rs:417-457`.
 *
 * - `atom` → all zero (no rotational degrees of freedom).
 * - `linear` → `B = rotationalConstantsGhz[1] * 1e9 Hz` (index 0 is ∞).
 * - nonlinear (spherical / symmetric / asymmetric top) → geometric mean of the
 *   three finite constants.
 */
function rotationalPartition(
  rotType: RotorType,
  rotationalConstantsGhz: [number, number, number],
  symmetryNumber: number,
  temperatureK: number
): { eRotHa: number; sRotHaPerK: number; cvRotHaPerK: number } {
  const r = rHartreePerMolK();
  const sigma = symmetryNumber;
  const t = temperatureK;

  if (rotType === 'atom') {
    return { eRotHa: 0, sRotHaPerK: 0, cvRotHaPerK: 0 };
  }

  if (rotType === 'linear') {
    // PySCF convention: B = rotational_constants_ghz[1] * 1e9 (index 0 is ∞)
    const bHz = rotationalConstantsGhz[1] * 1e9;
    const qRot = (BOLTZMANN * t) / (sigma * PLANCK * bHz);
    return {
      eRotHa: r * t,
      sRotHaPerK: r * (1 + Math.log(qRot)),
      cvRotHaPerK: r,
    };
  }

  // Nonlinear rotor (spherical / symmetric / asymmetric top)
  const abcHz: [number, number, number] = [
    rotationalConstantsGhz[0] * 1e9,
    rotationalConstantsGhz[1] * 1e9,
    rotationalConstantsGhz[2] * 1e9,
  ];
  const prodAbc = abcHz[0] * abcHz[1] * abcHz[2];
  const ktOverH = (BOLTZMANN * t) / PLANCK;
  const qRot =
    (Math.pow(ktOverH, 1.5) * Math.sqrt(PI)) / (sigma * Math.sqrt(prodAbc));
  return {
    eRotHa: 1.5 * r * t,
    sRotHaPerK: r * (1.5 + Math.log(qRot)),
    cvRotHaPerK: 1.5 * r,
  };
}

// ============================================================================
// Internal: vibrational partition
// ============================================================================

/**
 * Vibrational partition function contributions (quantum harmonic oscillator).
 *
 * Imaginary frequencies (ν ≤ 0) are skipped from every sum and counted.
 *
 * Rust reference: `thermochemistry.rs:504-564`.
 */
function vibrationalPartition(
  frequencyWavenumberCm1: number[],
  temperatureK: number
): {
  zpeHa: number;
  eVibThermalHa: number;
  sVibHaPerK: number;
  cvVibHaPerK: number;
  nUsed: number;
  nImag: number;
} {
  const r = rHartreePerMolK();
  // cm/s conversion for c (LIGHT_SPEED_SI is in m/s)
  const cCmPerS = LIGHT_SPEED_SI * 100;

  let zpe = 0;
  let eThermal = 0;
  let sVib = 0;
  let cvVib = 0;
  let nUsed = 0;
  let nImag = 0;

  for (const nu of frequencyWavenumberCm1) {
    // Skip imaginary modes (negative-real by convention)
    if (nu <= 0) {
      nImag += 1;
      continue;
    }
    nUsed += 1;

    // Vibrational temperature Θ_k = h·c·ν / k_B  (ν in cm⁻¹, c in cm/s)
    const thetaK = (PLANCK * cCmPerS * nu) / BOLTZMANN;

    // ZPE contribution is temperature-independent.
    zpe += 0.5 * r * thetaK;

    // Reduced temperature u_k = Θ_k / T
    const uK = thetaK / temperatureK;
    // exp(-u_k) ∈ [0, 1]; underflows to 0 for u_k ≳ 700
    const eNeg = Math.exp(-uK);
    if (eNeg <= 0) {
      // Low-T limit: thermal contributions → 0 for this mode.
      continue;
    }
    const oneMinus = 1 - eNeg;
    if (oneMinus <= 0) {
      // Guard against numerical underflow (u_k near 0 in double precision).
      continue;
    }

    eThermal += (r * temperatureK * uK * eNeg) / oneMinus;
    sVib += r * ((uK * eNeg) / oneMinus - Math.log(oneMinus));
    cvVib += (r * uK * uK * eNeg) / (oneMinus * oneMinus);
  }

  return {
    zpeHa: zpe,
    eVibThermalHa: eThermal,
    sVibHaPerK: sVib,
    cvVibHaPerK: cvVib,
    nUsed,
    nImag,
  };
}

// ============================================================================
// Internal: electronic entropy
// ============================================================================

/**
 * Electronic entropy for a ground-state multiplet.
 *
 * `S_elec = R · ln(multiplicity)`. For a closed-shell singlet
 * (multiplicity = 1) the result is exactly 0.
 *
 * Rust reference: `thermochemistry.rs:580-582`.
 */
function electronicEntropy(multiplicity: number): number {
  return rHartreePerMolK() * Math.log(multiplicity);
}

// ============================================================================
// Public: top-level recompute
// ============================================================================

/**
 * Recompute RRHO thermochemistry at a new (T, P) using the frozen frequencies,
 * rotor type, rotational constants, mass, and electronic energy from an
 * earlier WASM `FrequencyResult`.
 *
 * Pure function — does NOT call the worker, does NOT mutate the input.
 *
 * @throws Error if `temperatureK <= 0` or `pressurePa <= 0`.
 *
 * @see crates/qc-core/src/thermochemistry.rs (authoritative source)
 */
export function recomputeThermochemistry(
  input: RecomputeThermochemistryInput
): FrequencyThermochemistry {
  const { result, temperatureK, pressurePa } = input;

  if (!Number.isFinite(temperatureK) || temperatureK <= 0) {
    throw new Error(
      `recomputeThermochemistry: temperatureK must be > 0 (got ${temperatureK})`
    );
  }
  if (!Number.isFinite(pressurePa) || pressurePa <= 0) {
    throw new Error(
      `recomputeThermochemistry: pressurePa must be > 0 (got ${pressurePa})`
    );
  }

  const symmetryNumber =
    input.symmetryNumberOverride ?? result.thermochemistry.symmetryNumber;
  const multiplicity =
    input.multiplicityOverride ?? result.thermochemistry.multiplicity;
  const totalMassAmu = result.thermochemistry.totalMassAmu;

  const trans = translationalPartition(totalMassAmu, temperatureK, pressurePa);
  const rot = rotationalPartition(
    result.rotorType,
    result.rotationalConstantsGhz,
    symmetryNumber,
    temperatureK
  );
  const vib = vibrationalPartition(result.frequenciesCm1, temperatureK);
  const sElec = electronicEntropy(multiplicity);

  const r = rHartreePerMolK();

  // Totals — matches `thermochemistry.rs:725-746`.
  const internalEnergyHa =
    result.electronicEnergyHa +
    vib.zpeHa +
    trans.eTransHa +
    rot.eRotHa +
    vib.eVibThermalHa;
  const enthalpyHa = internalEnergyHa + r * temperatureK; // H = U + RT
  const entropyHaPerK =
    trans.sTransHaPerK + rot.sRotHaPerK + vib.sVibHaPerK + sElec;
  const gibbsHa = enthalpyHa - temperatureK * entropyHaPerK;
  const cvHaPerK = trans.cvTransHaPerK + rot.cvRotHaPerK + vib.cvVibHaPerK;
  const cpHaPerK = cvHaPerK + r; // Cp = Cv + R
  const e0kHa = result.electronicEnergyHa + vib.zpeHa;

  // PySCF convention for H_rot / H_vib:
  //   H_rot = E_rot, H_vib = ZPE + E_vib_thermal (matches Rust lines 708, 715).
  const hRotHa = rot.eRotHa;
  const hVibHa = vib.zpeHa + vib.eVibThermalHa;

  return {
    // Input echoes
    temperatureK,
    pressurePa,
    symmetryNumber,
    multiplicity,
    totalMassAmu,
    nVibModesUsed: vib.nUsed,
    nImag: vib.nImag,

    // Zero-point + 0 K
    zpeHa: vib.zpeHa,
    e0kHa,

    // Totals
    internalEnergyHa,
    enthalpyHa,
    entropyHaPerK,
    gibbsHa,
    cvHaPerK,
    cpHaPerK,

    // Translational
    eTransHa: trans.eTransHa,
    hTransHa: trans.hTransHa,
    sTransHaPerK: trans.sTransHaPerK,
    cvTransHaPerK: trans.cvTransHaPerK,
    cpTransHaPerK: trans.cpTransHaPerK,

    // Rotational (H_rot = E_rot, Cp_rot = Cv_rot per PySCF convention)
    eRotHa: rot.eRotHa,
    hRotHa,
    sRotHaPerK: rot.sRotHaPerK,
    cvRotHaPerK: rot.cvRotHaPerK,
    cpRotHaPerK: rot.cvRotHaPerK,

    // Vibrational (H_vib = ZPE + thermal, Cp_vib = Cv_vib per PySCF convention)
    eVibThermalHa: vib.eVibThermalHa,
    hVibHa,
    sVibHaPerK: vib.sVibHaPerK,
    cvVibHaPerK: vib.cvVibHaPerK,
    cpVibHaPerK: vib.cvVibHaPerK,

    // Electronic (only entropy nonzero in RRHO)
    sElecHaPerK: sElec,
  };
}
