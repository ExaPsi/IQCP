/**
 * Shell display info derivation utility.
 *
 * Generates human-readable shell labels (e.g., "1s", "2p'") and descriptions
 * from raw WASM shell data. Extracted from ShellList.tsx to avoid the
 * react-refresh/only-export-components lint warning when exporting non-component
 * functions alongside React components.
 *
 * @module components/basisExplorer/shellDisplayInfo
 */

import type { ShellInfo } from './basisData';
import type { BasisShellInfo } from '../../worker/protocol';

/**
 * Generate display labels and descriptions from WASM shell data.
 *
 * Since the WASM layer provides raw shell info (angular momentum + primitives),
 * we derive human-readable labels (e.g., "1s", "2p") and descriptions from
 * the shell sequence and element context.
 *
 * @param wasmShells - Shell data from WASM (exponents, coefficients, etc.)
 * @param z - Atomic number of the element
 * @returns Array of display-friendly shell info objects
 */
export function deriveShellDisplayInfo(
  wasmShells: BasisShellInfo[],
  z: number
): ShellInfo[] {
  // Track principal quantum numbers per angular momentum type
  const amCounters: Record<string, number> = { s: 0, p: 0, d: 0 };

  // Determine shell role based on element and position
  const isHeavyAtom = z >= 3; // Li and above have core electrons
  const isThirdRow = z >= 11; // Na-Ar

  return wasmShells.map((ws) => {
    const amLabel = ws.angularMomentumLabel;
    amCounters[amLabel] = (amCounters[amLabel] ?? 0) + 1;
    const count = amCounters[amLabel];

    // Derive principal quantum number from shell ordering
    let principalN: number;
    if (amLabel === 's') {
      principalN = count;
    } else if (amLabel === 'p') {
      // p orbitals start from n=2
      principalN = count + 1;
    } else {
      // d orbitals start from n=3
      principalN = count + 2;
    }

    // For split-valence bases (3-21G, 6-31G, etc.), multiple shells of the
    // same angular momentum appear consecutively. Use prime notation for
    // the outer valence shell.
    let label = `${principalN}${amLabel}`;

    // Check if this is a duplicate principal quantum number for split-valence
    // (The prior shell with same n and am already exists)
    const priorSameShells = wasmShells
      .slice(0, wasmShells.indexOf(ws))
      .filter((s) => s.angularMomentumLabel === amLabel);
    const priorCount = priorSameShells.length;

    // In split-valence bases, if there's already a shell with the same derived
    // principal number, mark this one as primed (outer)
    if (priorCount > 0) {
      // Recalculate: if prior shells produced the same principalN, add prime
      const priorPrincipalNs = priorSameShells.map((_, idx) => {
        const c = idx + 1;
        if (amLabel === 's') return c;
        if (amLabel === 'p') return c + 1;
        return c + 2;
      });
      if (priorPrincipalNs.includes(principalN)) {
        label = `${principalN}${amLabel}'`;
      }
    }

    // Generate description based on shell role
    let description: string;
    if (amLabel === 'd') {
      description = 'Polarization d orbital';
    } else if (!isHeavyAtom) {
      // H, He: all s shells are core/valence
      if (label.includes("'")) {
        description = `Outer valence ${amLabel}`;
      } else if (amCounters[amLabel] > 1) {
        description = `Inner valence ${amLabel}`;
      } else {
        description = `Core/valence ${amLabel} orbital`;
      }
    } else if (principalN === 1) {
      description = `Core ${amLabel} orbital`;
    } else if (label.includes("'")) {
      description = `Outer valence ${amLabel}`;
    } else if (isThirdRow && principalN === 2) {
      description = `Inner ${amLabel} orbital`;
    } else {
      description = `Valence ${amLabel} orbital`;
    }

    return {
      label,
      angularMomentum: ws.angularMomentum,
      angularMomentumLabel: amLabel as 's' | 'p' | 'd',
      nPrimitives: ws.nPrimitives,
      description,
    };
  });
}
