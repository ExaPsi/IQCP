/**
 * Shell type matching utility for cross-basis comparison.
 *
 * When comparing the "1s" shell across STO-3G, 3-21G, and 6-31G, each basis
 * set may have a different number of shells for that angular momentum type.
 * STO-3G has one 1s shell (3 primitives), while 3-21G has two (inner: 2
 * primitives, outer: 1 primitive). This utility matches shells across basis
 * sets by angular momentum and principal quantum number, labeling each with
 * its structural role (core, inner valence, outer valence).
 *
 * @module components/basisExplorer/shellMatching
 */

import type { BasisShellInfo } from '../../worker/protocol';
import type { ShellInfo } from './basisData';

/**
 * A matched shell with its index and structural role.
 */
export interface MatchedShell {
  /** Shell index in the WASM shell array */
  shellIndex: number;
  /** Structural role for labeling: "core", "inner valence", "outer valence", "valence", "polarization" */
  role: string;
}

/**
 * Available shell types for a given element across all basis sets.
 *
 * Used to populate the shell type selector in comparison mode.
 */
export interface ShellTypeOption {
  /** Display label, e.g., "1s", "2s", "2p" */
  label: string;
  /** Angular momentum (0=s, 1=p, 2=d) */
  angularMomentum: number;
  /** Principal quantum number */
  principalN: number;
}

/**
 * Find all shells in a basis set that match a target shell type.
 *
 * For a given shell type label (e.g., "1s" or "2p"), identifies all shells
 * in the basis set with that angular momentum in the appropriate principal
 * quantum number range. For split-valence bases, this returns multiple shells
 * (inner and outer valence) for the same angular momentum type.
 *
 * @param wasmShells - Shell data from WASM for one element+basis combination
 * @param shellDisplayInfo - Derived display labels for each shell
 * @param targetLabel - Target shell type label (e.g., "1s", "2p")
 * @returns Array of matched shell indices with role labels
 *
 * @example
 * ```typescript
 * // For H in 3-21G, matching "1s":
 * // Returns [
 * //   { shellIndex: 0, role: "inner valence" },  // 1s (2 primitives)
 * //   { shellIndex: 1, role: "outer valence" },   // 1s' (1 primitive)
 * // ]
 * ```
 */
export function findMatchingShells(
  _wasmShells: BasisShellInfo[],
  shellDisplayInfo: ShellInfo[],
  targetLabel: string,
): MatchedShell[] {
  const matches: MatchedShell[] = [];

  // Parse the target label to extract the base label without prime
  // "1s" -> base "1s", "2p" -> base "2p", "1s'" -> base "1s"
  const baseLabel = targetLabel.replace("'", '');

  for (let i = 0; i < shellDisplayInfo.length; i++) {
    const di = shellDisplayInfo[i];
    // Match if the base label (without prime) matches the target
    const shellBase = di.label.replace("'", '');

    if (shellBase === baseLabel) {
      // Determine role from the display info description
      const role = categorizeRole(di);
      matches.push({ shellIndex: i, role });
    }
  }

  return matches;
}

/**
 * Categorize the structural role of a shell from its display info.
 *
 * @param di - Shell display info
 * @returns Human-readable role string
 */
function categorizeRole(di: ShellInfo): string {
  const desc = di.description.toLowerCase();
  if (desc.includes('polarization')) return 'polarization';
  if (desc.includes('outer valence')) return 'outer valence';
  if (desc.includes('inner valence')) return 'inner valence';
  if (desc.includes('core/valence')) return 'valence';
  if (desc.includes('core')) return 'core';
  if (desc.includes('valence')) return 'valence';
  if (desc.includes('inner ')) return 'inner';
  return 'valence';
}

/**
 * Get the unique shell types available for an element in a given basis set.
 *
 * Returns deduplicated shell type labels (e.g., ["1s", "2s", "2p"]) with
 * angular momentum info. The prime notation is collapsed: "1s" and "1s'" both
 * map to the "1s" entry (since they represent inner/outer parts of the same
 * orbital type).
 *
 * @param shellDisplayInfo - Derived display labels for each shell
 * @returns Array of unique shell type options
 */
export function getAvailableShellTypes(
  shellDisplayInfo: ShellInfo[],
): ShellTypeOption[] {
  const seen = new Set<string>();
  const options: ShellTypeOption[] = [];

  for (const di of shellDisplayInfo) {
    // Collapse prime notation: "1s'" -> "1s"
    const baseLabel = di.label.replace("'", '');

    if (!seen.has(baseLabel)) {
      seen.add(baseLabel);

      // Extract principal quantum number from label
      const match = baseLabel.match(/^(\d+)/);
      const principalN = match ? parseInt(match[1], 10) : 1;

      options.push({
        label: baseLabel,
        angularMomentum: di.angularMomentum,
        principalN,
      });
    }
  }

  return options;
}
