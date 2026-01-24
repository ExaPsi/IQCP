/**
 * SystemInfoPanel - Display molecular system information.
 *
 * Shows the selected preset system's geometry (atomic coordinates),
 * basis set, and other properties for educational transparency.
 *
 * @module components/scf/SystemInfoPanel
 */

import { getPreset } from '../../worker/presets';
import { Math } from '../common/Math';

/**
 * Props for SystemInfoPanel component.
 */
export interface SystemInfoPanelProps {
  /** Selected system ID */
  systemId: string;
  /** Additional CSS classes */
  className?: string;
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
export function SystemInfoPanel({ systemId, className = '' }: SystemInfoPanelProps) {
  const preset = getPreset(systemId);

  if (!preset) {
    return (
      <div className={`bg-slate-50 rounded-xl p-4 ${className}`}>
        <p className="text-slate-500 text-sm">Unknown system: {systemId}</p>
      </div>
    );
  }

  const { geometry, basis_id, nbf, nelec, e_nuc, description } = preset;

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
          <div className="text-slate-500">Basis Set:</div>
          <div className="font-mono text-slate-800">{basis_id.toUpperCase()}</div>

          <div className="text-slate-500">Basis Functions:</div>
          <div className="font-mono text-slate-800">{nbf}</div>

          <div className="text-slate-500">Electrons:</div>
          <div className="font-mono text-slate-800">{nelec}</div>

          <div className="text-slate-500">
            <Math>{String.raw`V_{nn}`}</Math>:
          </div>
          <div className="font-mono text-slate-800">{e_nuc.toFixed(8)} Ha</div>
        </div>

        {/* Geometry Table */}
        {geometry && geometry.atoms.length > 0 && (
          <div>
            <h4 className="text-xs font-medium text-slate-600 mb-2">
              Atomic Coordinates ({geometry.units})
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
                  {geometry.atoms.map((atom, index) => (
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
