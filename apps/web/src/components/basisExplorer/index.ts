/**
 * Basis Explorer module components.
 *
 * @module components/basisExplorer
 */

export { ElementSelector } from './ElementSelector';
export { BasisSelector } from './BasisSelector';
export { ShellList } from './ShellList';
export { deriveShellDisplayInfo } from './shellDisplayInfo';
export { RadialProfilePlot } from './RadialProfilePlot';
export type { RadialProfilePlotProps } from './RadialProfilePlot';
export { ShellDataTable } from './ShellDataTable';
export type { ShellDataTableProps } from './ShellDataTable';
export { primitiveNormalization } from './normalization';
export type { ShellInfo, ElementInfo, BasisSetOption } from './basisData';
export { ELEMENTS, BASIS_SETS, getElement, isElementSupported } from './basisData';

// US-053: Exponent Sliders
export { ExponentSliders } from './ExponentSliders';
export type { ExponentSlidersProps } from './ExponentSliders';
export { evaluateRadialProfileTS, computeOverlapMetric } from './radialMath';

// US-052: Basis Set Comparison Mode
export { BasisComparisonControls } from './BasisComparisonControls';
export type { BasisComparisonControlsProps } from './BasisComparisonControls';
export { BasisComparisonPlot } from './BasisComparisonPlot';
export type { BasisComparisonPlotProps } from './BasisComparisonPlot';
export { DiscussionPrompt } from './DiscussionPrompt';
export type { DiscussionPromptProps } from './DiscussionPrompt';
export { COMPARISON_COLORS, SHELL_DASH_PATTERNS } from './comparisonColors';
export { findMatchingShells, getAvailableShellTypes } from './shellMatching';
export type { MatchedShell, ShellTypeOption } from './shellMatching';

// US-054: Overlap vs. Distance Plot
export { OverlapDistancePlot } from './OverlapDistancePlot';
export type { OverlapDistancePlotProps } from './OverlapDistancePlot';
export { OverlapDistanceControls } from './OverlapDistanceControls';
export type { OverlapDistanceControlsProps, OverlapTraceConfig } from './OverlapDistanceControls';
