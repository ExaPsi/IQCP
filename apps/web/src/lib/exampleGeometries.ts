/**
 * Example Geometries for Custom Geometry Mode
 *
 * Pre-defined molecular geometries that users can load to explore
 * the custom geometry functionality without typing coordinates.
 *
 * @module lib/exampleGeometries
 */

import type { GeometryInput, BasisSetName } from '../worker/protocol';

/**
 * An example geometry with metadata.
 */
export interface ExampleGeometry {
  /** Unique identifier */
  id: string;
  /** Display label with Unicode subscripts */
  label: string;
  /** Brief description */
  description: string;
  /** The geometry data */
  geometry: GeometryInput;
  /** Recommended basis set for this geometry */
  recommendedBasis: BasisSetName;
}

/**
 * Available example geometries.
 *
 * Geometries are in Bohr by default.
 * Coordinates sourced from optimized structures at specified basis.
 */
export const EXAMPLE_GEOMETRIES: ExampleGeometry[] = [
  {
    id: 'h2_eq',
    label: 'H\u2082 (equilibrium)',
    description: 'Hydrogen molecule at equilibrium bond length (R=1.4 bohr)',
    geometry: {
      atoms: [
        { symbol: 'H', xyz: [0.0, 0.0, 0.0] },
        { symbol: 'H', xyz: [0.0, 0.0, 1.4] },
      ],
      units: 'bohr',
    },
    recommendedBasis: 'sto-3g',
  },
  {
    id: 'heh_plus',
    label: 'HeH\u207A',
    description: 'Helium hydride cation (equilibrium geometry)',
    geometry: {
      atoms: [
        { symbol: 'He', xyz: [0.0, 0.0, 0.0] },
        { symbol: 'H', xyz: [0.0, 0.0, 1.46] },
      ],
      units: 'bohr',
    },
    recommendedBasis: 'sto-3g',
  },
  {
    id: 'lih',
    label: 'LiH',
    description: 'Lithium hydride (equilibrium geometry)',
    geometry: {
      atoms: [
        { symbol: 'Li', xyz: [0.0, 0.0, 0.0] },
        { symbol: 'H', xyz: [0.0, 0.0, 3.015] },
      ],
      units: 'bohr',
    },
    recommendedBasis: 'sto-3g',
  },
  {
    id: 'h2o',
    label: 'H\u2082O',
    description: 'Water molecule (experimental geometry)',
    geometry: {
      atoms: [
        { symbol: 'O', xyz: [0.0, 0.0, 0.1173] },
        { symbol: 'H', xyz: [0.0, 1.4299, -0.4692] },
        { symbol: 'H', xyz: [0.0, -1.4299, -0.4692] },
      ],
      units: 'bohr',
    },
    recommendedBasis: 'sto-3g',
  },
  {
    id: 'nh3',
    label: 'NH\u2083',
    description: 'Ammonia molecule (experimental geometry)',
    geometry: {
      atoms: [
        { symbol: 'N', xyz: [0.0, 0.0, 0.2201] },
        { symbol: 'H', xyz: [0.0, 1.7717, -0.5135] },
        { symbol: 'H', xyz: [1.5344, -0.8858, -0.5135] },
        { symbol: 'H', xyz: [-1.5344, -0.8858, -0.5135] },
      ],
      units: 'bohr',
    },
    recommendedBasis: 'sto-3g',
  },
  {
    id: 'ch4',
    label: 'CH\u2084',
    description: 'Methane molecule (tetrahedral geometry)',
    geometry: {
      atoms: [
        { symbol: 'C', xyz: [0.0, 0.0, 0.0] },
        { symbol: 'H', xyz: [1.19, 1.19, 1.19] },
        { symbol: 'H', xyz: [-1.19, -1.19, 1.19] },
        { symbol: 'H', xyz: [-1.19, 1.19, -1.19] },
        { symbol: 'H', xyz: [1.19, -1.19, -1.19] },
      ],
      units: 'bohr',
    },
    recommendedBasis: 'sto-3g',
  },
  {
    id: 'hf',
    label: 'HF',
    description: 'Hydrogen fluoride (equilibrium geometry)',
    geometry: {
      atoms: [
        { symbol: 'H', xyz: [0.0, 0.0, 0.0] },
        { symbol: 'F', xyz: [0.0, 0.0, 1.7328] },
      ],
      units: 'bohr',
    },
    recommendedBasis: 'sto-3g',
  },
];

/**
 * Get an example geometry by ID.
 *
 * @param id - Example ID
 * @returns ExampleGeometry or undefined
 */
export function getExampleGeometry(id: string): ExampleGeometry | undefined {
  return EXAMPLE_GEOMETRIES.find((eg) => eg.id === id);
}

/**
 * Get all example geometries for a given basis set.
 *
 * Currently returns all examples since they work with all basis sets.
 * In the future, this could filter based on recommended basis.
 *
 * @param basis - Basis set name (for future filtering)
 * @returns Array of example geometries that work well with this basis
 */
export function getExamplesForBasis(basis: BasisSetName): ExampleGeometry[] {
  // All our examples work with all basis sets, but some may recommend specific ones
  // Use basis to filter in the future if needed
  void basis;
  return EXAMPLE_GEOMETRIES;
}
