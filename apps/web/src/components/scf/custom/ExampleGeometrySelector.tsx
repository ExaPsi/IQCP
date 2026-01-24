/**
 * ExampleGeometrySelector - Dropdown to select example geometries.
 *
 * Provides a dropdown with pre-defined molecular geometries that users
 * can load to explore custom geometry functionality without typing.
 *
 * @module components/scf/custom/ExampleGeometrySelector
 */

import {
  EXAMPLE_GEOMETRIES,
  type ExampleGeometry,
} from '../../../lib/exampleGeometries';

/**
 * Props for ExampleGeometrySelector component.
 */
export interface ExampleGeometrySelectorProps {
  /** Callback when an example is selected */
  onSelect: (example: ExampleGeometry) => void;
  /** Disable the selector */
  disabled?: boolean;
}

/**
 * Dropdown to select and load example geometries.
 *
 * When an example is selected, calls `onSelect` with the full
 * ExampleGeometry object, allowing the parent to populate the
 * geometry editor and optionally set the recommended basis.
 *
 * @example
 * ```tsx
 * const handleExample = (example: ExampleGeometry) => {
 *   setXyzText(geometryToXyzText(example.geometry));
 *   setUnits(example.geometry.units);
 *   setBasis(example.recommendedBasis);
 * };
 *
 * <ExampleGeometrySelector onSelect={handleExample} />
 * ```
 */
export function ExampleGeometrySelector({
  onSelect,
  disabled = false,
}: ExampleGeometrySelectorProps) {
  const handleChange = (event: React.ChangeEvent<HTMLSelectElement>) => {
    const selectedId = event.target.value;
    if (selectedId === '') return;

    const example = EXAMPLE_GEOMETRIES.find((eg) => eg.id === selectedId);
    if (example) {
      onSelect(example);
      // Reset select to placeholder after selection
      event.target.value = '';
    }
  };

  return (
    <div>
      <label
        htmlFor="example-geometry-select"
        className="block text-sm font-medium text-slate-700 mb-2"
      >
        Load Example
      </label>
      <select
        id="example-geometry-select"
        onChange={handleChange}
        disabled={disabled}
        defaultValue=""
        className="w-full px-3 py-2 border border-slate-300 rounded-lg bg-white text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
        aria-label="Select an example geometry to load"
      >
        <option value="" disabled>
          Select an example geometry...
        </option>
        {EXAMPLE_GEOMETRIES.map((example) => (
          <option key={example.id} value={example.id}>
            {example.label} - {example.description}
          </option>
        ))}
      </select>
      <p className="mt-1 text-xs text-slate-500">
        Load a pre-defined geometry to see the format and try computations
      </p>
    </div>
  );
}

export default ExampleGeometrySelector;
