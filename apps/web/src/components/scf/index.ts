/**
 * SCF Components Barrel Export
 *
 * Re-exports all SCF module components for convenient importing.
 *
 * @module components/scf
 */

// Main components
export { ScfControlsPanel } from './ScfControlsPanel';
export { ScfIterationTable } from './ScfIterationTable';
export { ScfEnergyPlot } from './ScfEnergyPlot';
export { ScfResidualPlot } from './ScfResidualPlot';
export { ScfResultDisplay } from './ScfResultDisplay';

// Internals mode components (US-018)
export { MatrixHeatmap } from './MatrixHeatmap';
export type { MatrixHeatmapProps, MatrixColorscale } from './MatrixHeatmap';

export { OrbitalEnergiesTable } from './OrbitalEnergiesTable';
export type { OrbitalEnergiesTableProps } from './OrbitalEnergiesTable';

export { DiagnosticPanel } from './DiagnosticPanel';
export type { DiagnosticPanelProps } from './DiagnosticPanel';

export { ScfInternalsPanel } from './ScfInternalsPanel';
export type { ScfInternalsPanelProps } from './ScfInternalsPanel';

export { SystemInfoPanel } from './SystemInfoPanel';
export type { SystemInfoPanelProps } from './SystemInfoPanel';

// Phase 2: PES scan components (US-040)
export { PesScanPanel } from './PesScanPanel';
export type { DiatomicInfo } from './PesScanPanel';
export { validatePesConfig } from './pesValidation';

// Phase 4: Internal coordinate PES scan (US-082)
export { CoordinateSelector, computeCoordinateValue, COORDINATE_UNITS, COORDINATE_LABELS, ELEMENT_SYMBOLS } from './CoordinateSelector';
export type { CoordinateType, CoordinateSelectorProps } from './CoordinateSelector';
export { ScanModeToggle } from './ScanModeToggle';
export type { ScanMode, ScanModeToggleProps } from './ScanModeToggle';

// Phase 2: PES curve visualization (US-041)
export { PesCurvePlot } from './PesCurvePlot';

// Phase 4: Coordinate tracking panel (US-083)
export { CoordinateTrackingPanel, extractAvailableCoordinates, extractCoordinateValues } from './CoordinateTrackingPanel';
export type { CoordinateTrackingPanelProps, CoordinateKey, SelectableCoordinate } from './CoordinateTrackingPanel';

// Phase 2: Orbital viewer UI (US-044)
export { OrbitalSelector } from './OrbitalSelector';
export type { OrbitalSelectorProps } from './OrbitalSelector';
export { IsovalueSlider } from './IsovalueSlider';
export type { IsovalueSliderProps } from './IsovalueSlider';

// Phase 2: Basis set comparison (US-045)
export { BasisComparisonPanel } from './BasisComparisonPanel';
export type { BasisComparisonPanelProps } from './BasisComparisonPanel';
export { BasisComparisonChart } from './BasisComparisonChart';
export { BasisComparisonTable } from './BasisComparisonTable';

// Phase 3: Density visualization (US-062)
export { DensityPanel } from './DensityPanel';
export type { DensityPanelProps } from './DensityPanel';
export { DensityCrossSection } from './DensityCrossSection';
export type { DensityCrossSectionProps } from './DensityCrossSection';

// Phase 3: DFT method selection (US-070)
export { MethodSelector } from './MethodSelector';
export type { MethodSelectorProps } from './MethodSelector';
export { DftInfoPanel } from './DftInfoPanel';
export type { DftInfoPanelProps } from './DftInfoPanel';
export { DftComparisonPanel } from './DftComparisonPanel';
export type { DftComparisonPanelProps } from './DftComparisonPanel';

// Phase 3: Population analysis (US-076)
export { PopulationTable } from './PopulationTable';
export type { PopulationTableProps } from './PopulationTable';

// Phase 3: Geometry optimization (US-075)
export { OptimizeButton } from './OptimizeButton';
export type { OptimizeButtonProps } from './OptimizeButton';
export { OptimizationProgressPlot } from './OptimizationProgressPlot';
export type { OptimizationProgressPlotProps } from './OptimizationProgressPlot';
export { OptimizationResultPanel } from './OptimizationResultPanel';
export type { OptimizationResultPanelProps } from './OptimizationResultPanel';

// Phase 5: Frequency Tab UI (US-102) — components
export { FrequencyTab } from './FrequencyTab';
export type { FrequencyTabProps } from './FrequencyTab';
export {
  DEFAULT_LOCAL_FREQUENCY_STATE,
} from './frequencyTabState';
export type { LocalFrequencyState } from './frequencyTabState';

export { FrequencyPanel } from './FrequencyPanel';
export type { FrequencyPanelProps } from './FrequencyPanel';
export {
  buildFrequencyRows,
  sortFrequencyRows,
  countImaginary,
  firstImaginaryIndex,
  formatFrequency,
} from './frequencyPanelLogic';
export type {
  FrequencyRow,
  FrequencyColumnKey,
  SortDirection,
} from './frequencyPanelLogic';

export { FrequencySummaryPanel } from './FrequencySummaryPanel';
export type { FrequencySummaryPanelProps } from './FrequencySummaryPanel';
export {
  dipoleMagnitude,
  isotropicPolarizability,
  formatRotorType,
  formatRotConst,
} from './frequencySummaryLogic';

export { SpectrumPlot } from './SpectrumPlot';
export type { SpectrumPlotProps } from './SpectrumPlot';
export {
  buildWavenumberGrid,
  computeBroadenedSpectrum,
  findNearestMode,
  SPECTRUM_GRID_POINTS,
} from './spectrumMath';

export { NormalModeViewer } from './NormalModeViewer';
export type { NormalModeViewerProps } from './NormalModeViewer';

export { ThermochemistryPanel } from './ThermochemistryPanel';
export type { ThermochemistryPanelProps } from './ThermochemistryPanel';
export {
  formatEnergyValue,
  formatEntropyValue,
} from './thermochemistryFormat';

export { ImaginaryWarningBanner } from './ImaginaryWarningBanner';
export type { ImaginaryWarningBannerProps } from './ImaginaryWarningBanner';
