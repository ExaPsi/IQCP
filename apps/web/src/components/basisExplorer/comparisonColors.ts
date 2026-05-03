/**
 * Color palette for basis set comparison mode.
 *
 * Five distinct, colorblind-accessible colors for identifying basis sets
 * in overlaid comparison plots. Colors are high-contrast on white backgrounds
 * and supplemented by dash patterns for accessibility.
 *
 * @module components/basisExplorer/comparisonColors
 * @see US-052 Basis Set Comparison Mode
 */

/**
 * Distinct colors for each basis set in comparison mode.
 *
 * Keyed by internal basis set name (lowercase).
 * Selected for accessibility: distinguishable by people with
 * common color vision deficiencies (protanopia, deuteranopia).
 */
export const COMPARISON_COLORS: Record<string, string> = {
  'sto-3g': '#2563eb',  // blue-600
  '3-21g': '#dc2626',   // red-600
  '6-31g': '#16a34a',   // green-600
  '6-31g*': '#9333ea',  // purple-600
  '6-31+g*': '#ea580c', // orange-600
  'cc-pvdz': '#0891b2', // cyan-600
};

/**
 * Dash patterns for distinguishing inner/outer valence shells.
 *
 * When a basis set has multiple shells for the same angular momentum
 * (split-valence), the first (inner) shell uses a solid line and
 * subsequent (outer) shells use dash patterns.
 *
 * Plotly dash values: 'solid', 'dash', 'dot', 'dashdot', 'longdash'
 */
export const SHELL_DASH_PATTERNS: string[] = [
  'solid',    // First matching shell (core or inner valence)
  'dash',     // Second matching shell (outer valence)
  'dot',      // Third (rare, but for completeness)
  'dashdot',  // Fourth
  'longdash', // Fifth
];
