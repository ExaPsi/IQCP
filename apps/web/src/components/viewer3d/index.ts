/**
 * 3D Molecular Viewer components.
 *
 * This barrel export serves as the code-split entry point for Three.js
 * dependencies. Import this module via React.lazy() to ensure Three.js
 * is loaded on demand and NOT included in the initial bundle.
 *
 * Usage (in US-035 Module E integration):
 *   const MoleculeViewer = React.lazy(() =>
 *     import('./components/viewer3d').then(m => ({ default: m.MoleculeViewer }))
 *   );
 *
 * @module viewer3d
 */

export { MoleculeViewer } from './MoleculeViewer';
export type { MoleculeViewerProps, CameraState } from './MoleculeViewer';

export { AtomSpheres } from './AtomSpheres';
export type { Atom3D, AtomSpheresProps } from './AtomSpheres';

export { BondCylinders } from './BondCylinders';
export type { BondCylindersProps } from './BondCylinders';

export { AtomLabels } from './AtomLabels';
export type { AtomLabelsProps } from './AtomLabels';

export { OrbitalSurface } from './OrbitalSurface';
export type { OrbitalSurfaceProps, LobeMeshData } from './OrbitalSurface';

export { DensitySurface } from './DensitySurface';
export type { DensitySurfaceProps } from './DensitySurface';

export { GhostAtoms } from './GhostAtoms';
export type { GhostAtomsProps } from './GhostAtoms';

// Phase 5 (US-102): Normal-mode displacement arrows
export { DisplacementArrows } from './DisplacementArrows';
export type { DisplacementArrowsProps } from './DisplacementArrows';

// Phase 5 (US-102 UX refinement): Animated atoms following normal-mode coordinate
export { AnimatedMolecule } from './AnimatedMolecule';
export type { AnimatedMoleculeProps } from './AnimatedMolecule';

export { detectBonds } from './bondDetection';
export type { Bond } from './bondDetection';
