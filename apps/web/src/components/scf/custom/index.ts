/**
 * Custom Geometry UI Components
 *
 * Components for the custom geometry input mode in the SCF Sandbox.
 * These components enable users to enter molecular geometries manually
 * and compute integrals on-the-fly.
 *
 * @module components/scf/custom
 */

// Input mode toggle
export { InputModeToggle } from './InputModeToggle';
export type { InputModeToggleProps, InputMode } from './InputModeToggle';

// Geometry editor
export { GeometryEditor } from './GeometryEditor';
export type { GeometryEditorProps } from './GeometryEditor';

// Basis set selector
export { BasisSetSelector } from './BasisSetSelector';
export type { BasisSetSelectorProps } from './BasisSetSelector';

// Basis type toggle (Cartesian vs Spherical)
export { BasisTypeToggle } from './BasisTypeToggle';
export type { BasisTypeToggleProps, BasisInfo } from './BasisTypeToggle';

// Integral compute button
export { IntegralComputeButton } from './IntegralComputeButton';
export type { IntegralComputeButtonProps } from './IntegralComputeButton';

// Example geometry selector
export { ExampleGeometrySelector } from './ExampleGeometrySelector';
export type { ExampleGeometrySelectorProps } from './ExampleGeometrySelector';
