/**
 * ScanModeToggle - Toggle between Rigid and Relaxed PES scan modes.
 *
 * Provides a radio-button toggle with brief educational explanations
 * of each scan mode. Follows the same styling pattern as BasisTypeToggle.
 *
 * @module components/scf/ScanModeToggle
 * @see US-082 Coordinate Selector + Scan Mode UI
 */

// ============================================================================
// Types
// ============================================================================

/**
 * Scan mode for PES scanning.
 * Matches the scanMode field on PesScanInternalRequest.
 */
export type ScanMode = 'rigid' | 'relaxed';

/**
 * Props for ScanModeToggle component.
 */
export interface ScanModeToggleProps {
  /** Currently selected scan mode */
  scanMode: ScanMode;
  /** Handler for scan mode changes */
  onScanModeChange: (mode: ScanMode) => void;
  /** Disable the toggle */
  disabled?: boolean;
}

// ============================================================================
// Component
// ============================================================================

/**
 * ScanModeToggle - Toggle between Rigid and Relaxed PES scan modes.
 *
 * Rigid scans fix all other coordinates while varying the scanned one.
 * Relaxed scans optimize all other coordinates at each scan point.
 *
 * @example
 * ```tsx
 * <ScanModeToggle
 *   scanMode={mode}
 *   onScanModeChange={setMode}
 * />
 * ```
 */
export function ScanModeToggle({
  scanMode,
  onScanModeChange,
  disabled = false,
}: ScanModeToggleProps) {
  return (
    <div className="space-y-2">
      {/* Label */}
      <label className="block text-xs font-medium text-slate-700">
        Scan Mode
      </label>

      {/* Toggle buttons */}
      <div className="flex rounded-lg bg-slate-100 p-1">
        <button
          type="button"
          onClick={() => onScanModeChange('rigid')}
          disabled={disabled}
          className={`flex-1 px-3 py-1.5 rounded-md text-xs font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
            scanMode === 'rigid'
              ? 'bg-white text-slate-800 shadow-sm'
              : 'text-slate-600 hover:text-slate-800'
          } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
          aria-pressed={scanMode === 'rigid'}
        >
          Rigid
        </button>
        <button
          type="button"
          onClick={() => onScanModeChange('relaxed')}
          disabled={disabled}
          className={`flex-1 px-3 py-1.5 rounded-md text-xs font-medium transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 ${
            scanMode === 'relaxed'
              ? 'bg-white text-slate-800 shadow-sm'
              : 'text-slate-600 hover:text-slate-800'
          } ${disabled ? 'opacity-50 cursor-not-allowed' : ''}`}
          aria-pressed={scanMode === 'relaxed'}
        >
          Relaxed
        </button>
      </div>

      {/* Description */}
      <p className="text-xs text-slate-500">
        {scanMode === 'rigid' ? (
          <>
            <strong>Rigid scan:</strong> Only the selected coordinate changes; all other degrees of freedom are frozen.
          </>
        ) : (
          <>
            <strong>Relaxed scan:</strong> All other degrees of freedom are optimized at each scan point. Slower but physically more realistic.
          </>
        )}
      </p>
    </div>
  );
}
