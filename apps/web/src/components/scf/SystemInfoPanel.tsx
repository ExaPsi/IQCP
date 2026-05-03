/**
 * SystemInfoPanel - Display molecular system information.
 *
 * Shows the selected preset system's geometry (atomic coordinates),
 * basis set, and other properties for educational transparency.
 *
 * Supports override props so the panel reflects the user's current
 * selections (e.g. a different basis set) rather than only the
 * hardcoded preset values.
 *
 * @module components/scf/SystemInfoPanel
 */

import { getPreset } from '../../worker/presets';
import { Math } from '../common/Math';

/**
 * Human-readable basis set label mapping.
 *
 * Converts lowercase internal identifiers to the conventional
 * typographic form used in chemistry literature.
 */
const BASIS_LABELS: Record<string, string> = {
  'sto-3g': 'STO-3G',
  '3-21g': '3-21G',
  '6-31g': '6-31G',
  '6-31g*': '6-31G*',
  '6-31+g*': '6-31+G*',
  'cc-pvdz': 'cc-pVDZ',
};

/**
 * Format a basis set identifier for display.
 *
 * Looks up the canonical label first; falls back to uppercase.
 */
function formatBasisLabel(basis: string): string {
  return BASIS_LABELS[basis.toLowerCase()] ?? basis.toUpperCase();
}

/**
 * Human-readable method name mapping.
 */
const METHOD_LABELS: Record<string, string> = {
  rhf: 'RHF',
  lda: 'LDA',
  b3lyp: 'B3LYP',
  'b3lyp-d3bj': 'B3LYP-D3(BJ)',
};

/**
 * Format a method identifier for display.
 */
function formatMethodLabel(method: string): string {
  return METHOD_LABELS[method.toLowerCase()] ?? method.toUpperCase();
}

/**
 * Props for SystemInfoPanel component.
 */
export interface SystemInfoPanelProps {
  /** Selected system ID */
  systemId: string;
  /** Additional CSS classes */
  className?: string;
  /** Override basis set name (when user selects a different basis than the preset) */
  basisOverride?: string;
  /** Override number of basis functions */
  nbfOverride?: number;
  /** Override method name for display */
  methodOverride?: string;
  /** Override electron count */
  nelecOverride?: number;
  /** Override nuclear repulsion energy (from SCF/KS result) */
  enucOverride?: number;
  /** Override description (e.g., molecule name from preset list) */
  descriptionOverride?: string;
  /** Override geometry (for custom input mode) */
  geometryOverride?: { atoms: Array<{ symbol: string; xyz: number[] }>; units: string };
}

/**
 * Format a number for display (6 decimal places).
 */
function formatCoord(value: number): string {
  return value.toFixed(6);
}

/**
 * SystemInfoPanel - Displays molecular geometry and system properties.
 *
 * Shows atomic coordinates, basis set info, and computed properties
 * like nuclear repulsion energy. Helps students understand what
 * input data goes into an SCF calculation.
 *
 * @example
 * ```tsx
 * <SystemInfoPanel systemId="h2o_sto3g" />
 * ```
 */
export function SystemInfoPanel({
  systemId,
  className = '',
  basisOverride,
  nbfOverride,
  methodOverride,
  nelecOverride,
  enucOverride,
  descriptionOverride,
  geometryOverride,
}: SystemInfoPanelProps) {
  const preset = getPreset(systemId);

  // In custom geometry mode a preset may not exist, but we can still
  // show geometry information from the override.
  if (!preset && !geometryOverride) {
    return (
      <div className={`bg-slate-50 rounded-xl p-4 ${className}`}>
        <p className="text-slate-500 text-sm">Unknown system: {systemId}</p>
      </div>
    );
  }

  // Resolve displayed values: overrides take priority over preset data.
  const presetBasis = preset?.basis_id ?? '';
  const displayBasis = basisOverride ?? presetBasis;
  const displayNbf = nbfOverride ?? preset?.nbf;
  const displayNelec = nelecOverride ?? preset?.nelec;
  const displayGeometry = geometryOverride ?? preset?.geometry;
  const description = descriptionOverride ?? preset?.description;
  const e_nuc = enucOverride ?? preset?.e_nuc;

  // Determine whether the basis differs from the preset default
  const basisDiffers =
    basisOverride !== undefined &&
    presetBasis !== '' &&
    basisOverride.toLowerCase() !== presetBasis.toLowerCase();

  return (
    <div className={`bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden ${className}`}>
      {/* Header */}
      <div className="px-4 py-3 bg-slate-50 border-b border-slate-200">
        <h3 className="text-sm font-semibold text-slate-700">System Information</h3>
      </div>

      <div className="p-4 space-y-4">
        {/* Description */}
        {description && (
          <p className="text-xs text-slate-600">{description}</p>
        )}

        {/* System Properties */}
        <div className="grid grid-cols-2 gap-x-4 gap-y-2 text-sm">
          {/* Method (only shown when an override is provided) */}
          {methodOverride && (
            <>
              <div className="text-slate-500">Method:</div>
              <div className="font-mono text-slate-800">{formatMethodLabel(methodOverride)}</div>
            </>
          )}

          <div className="text-slate-500">Basis Set:</div>
          <div className="font-mono text-slate-800">
            {formatBasisLabel(displayBasis)}
            {basisDiffers && (
              <span
                className="ml-1.5 inline-flex items-center rounded px-1.5 py-0.5 text-[10px] font-medium bg-amber-100 text-amber-800"
                title={`Preset default: ${formatBasisLabel(presetBasis)}`}
              >
                custom
              </span>
            )}
          </div>

          {displayNbf !== undefined && (
            <>
              <div className="text-slate-500">Basis Functions:</div>
              <div className="font-mono text-slate-800">{displayNbf}</div>
            </>
          )}

          {displayNelec !== undefined && (
            <>
              <div className="text-slate-500">Electrons:</div>
              <div className="font-mono text-slate-800">{displayNelec}</div>
            </>
          )}

          {e_nuc !== undefined && (
            <>
              <div className="text-slate-500">
                <Math>{String.raw`V_{nn}`}</Math>:
              </div>
              <div className="font-mono text-slate-800">{e_nuc.toFixed(8)} Ha</div>
            </>
          )}
        </div>

        {/* Geometry Table */}
        {displayGeometry && displayGeometry.atoms.length > 0 && (
          <div>
            <h4 className="text-xs font-medium text-slate-600 mb-2">
              Atomic Coordinates ({displayGeometry.units})
            </h4>
            <div className="overflow-x-auto">
              <table className="w-full text-xs">
                <thead>
                  <tr className="bg-slate-50">
                    <th className="px-2 py-1.5 text-left font-medium text-slate-600">Atom</th>
                    <th className="px-2 py-1.5 text-right font-medium text-slate-600">X</th>
                    <th className="px-2 py-1.5 text-right font-medium text-slate-600">Y</th>
                    <th className="px-2 py-1.5 text-right font-medium text-slate-600">Z</th>
                  </tr>
                </thead>
                <tbody className="font-mono">
                  {displayGeometry.atoms.map((atom, index) => (
                    <tr
                      key={index}
                      className={index % 2 === 0 ? 'bg-white' : 'bg-slate-50/50'}
                    >
                      <td className="px-2 py-1.5 text-slate-800 font-medium">
                        {atom.symbol}
                      </td>
                      <td className="px-2 py-1.5 text-right text-slate-700">
                        {formatCoord(atom.xyz[0])}
                      </td>
                      <td className="px-2 py-1.5 text-right text-slate-700">
                        {formatCoord(atom.xyz[1])}
                      </td>
                      <td className="px-2 py-1.5 text-right text-slate-700">
                        {formatCoord(atom.xyz[2])}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

export default SystemInfoPanel;
