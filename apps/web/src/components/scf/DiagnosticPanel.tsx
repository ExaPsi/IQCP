/**
 * DiagnosticPanel - SCF computation diagnostic display.
 *
 * Shows convergence status, iteration count, final energy, and
 * educational hints to help students understand their results.
 *
 * Features:
 * - Status badge (Converged/Not Converged/Aborted)
 * - Key metrics display
 * - Context-aware educational hints
 * - Clean, informative layout
 *
 * @module components/scf/DiagnosticPanel
 */

import { Math } from '../common/Math';
import type { ConvergenceProfile } from '../../worker/protocol';

/**
 * Props for DiagnosticPanel component.
 */
export interface DiagnosticPanelProps {
  /** Whether SCF converged */
  converged: boolean;
  /** Number of iterations performed */
  iterations: number;
  /** Final total energy in Hartree */
  energy: number;
  /** Convergence profile used */
  convergenceProfile: ConvergenceProfile;
  /** Whether DIIS was enabled */
  diisEnabled: boolean;
  /** Whether computation was aborted */
  aborted?: boolean;
  /** Additional CSS classes */
  className?: string;
}

/**
 * Convergence status type for styling and display.
 */
type StatusType = 'converged' | 'not_converged' | 'aborted';

/**
 * Determine the status type from props.
 */
function getStatusType(converged: boolean, aborted?: boolean): StatusType {
  if (aborted) return 'aborted';
  return converged ? 'converged' : 'not_converged';
}

/**
 * Get status badge styling.
 */
function getStatusStyles(statusType: StatusType): {
  bgClass: string;
  textClass: string;
  borderClass: string;
  label: string;
} {
  switch (statusType) {
    case 'converged':
      return {
        bgClass: 'bg-green-100',
        textClass: 'text-green-800',
        borderClass: 'border-green-200',
        label: 'Converged',
      };
    case 'not_converged':
      return {
        bgClass: 'bg-red-100',
        textClass: 'text-red-800',
        borderClass: 'border-red-200',
        label: 'Not Converged',
      };
    case 'aborted':
      return {
        bgClass: 'bg-amber-100',
        textClass: 'text-amber-800',
        borderClass: 'border-amber-200',
        label: 'Aborted',
      };
  }
}

/**
 * Generate educational hints based on the computation outcome.
 *
 * Provides context-aware suggestions for students to improve their
 * understanding or troubleshoot convergence issues.
 */
function getEducationalHints(
  statusType: StatusType,
  iterations: number,
  diisEnabled: boolean,
  convergenceProfile: ConvergenceProfile
): string[] {
  const hints: string[] = [];

  if (statusType === 'converged') {
    // Success hints
    if (diisEnabled && iterations <= 10) {
      hints.push(
        'DIIS significantly accelerated convergence. Try disabling it to see how many more iterations would be needed.'
      );
    }
    if (!diisEnabled && iterations > 15) {
      hints.push(
        'SCF converged without DIIS but required many iterations. DIIS could speed this up considerably.'
      );
    }
    if (convergenceProfile === 'tight') {
      hints.push(
        'Tight convergence achieved. This precision is appropriate for comparing energies between similar systems.'
      );
    }
  } else if (statusType === 'not_converged') {
    // Failure hints
    if (!diisEnabled) {
      hints.push(
        'Try enabling DIIS acceleration. It often helps challenging systems converge.'
      );
    }
    if (convergenceProfile === 'tight') {
      hints.push(
        'Try a looser convergence profile (medium or loose) to see if convergence is achievable at lower precision.'
      );
    }
    if (iterations >= 50) {
      hints.push(
        'Maximum iterations reached. The system may have convergence difficulties. Consider checking initial guess or trying different options.'
      );
    }
  } else if (statusType === 'aborted') {
    hints.push(
      'The calculation was stopped before completion. The displayed data shows partial results.'
    );
  }

  return hints;
}

/**
 * Status badge component.
 */
function StatusBadge({ statusType }: { statusType: StatusType }) {
  const styles = getStatusStyles(statusType);

  return (
    <span
      className={`inline-flex items-center px-3 py-1 rounded-full text-sm font-medium ${styles.bgClass} ${styles.textClass} border ${styles.borderClass}`}
      role="status"
      aria-label={`SCF status: ${styles.label}`}
    >
      {statusType === 'converged' && (
        <svg
          className="w-4 h-4 mr-1.5"
          fill="currentColor"
          viewBox="0 0 20 20"
          aria-hidden="true"
        >
          <path
            fillRule="evenodd"
            d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z"
            clipRule="evenodd"
          />
        </svg>
      )}
      {statusType === 'not_converged' && (
        <svg
          className="w-4 h-4 mr-1.5"
          fill="currentColor"
          viewBox="0 0 20 20"
          aria-hidden="true"
        >
          <path
            fillRule="evenodd"
            d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z"
            clipRule="evenodd"
          />
        </svg>
      )}
      {statusType === 'aborted' && (
        <svg
          className="w-4 h-4 mr-1.5"
          fill="currentColor"
          viewBox="0 0 20 20"
          aria-hidden="true"
        >
          <path
            fillRule="evenodd"
            d="M10 18a8 8 0 100-16 8 8 0 000 16zM9 9V5a1 1 0 012 0v4a1 1 0 01-2 0zm1 6a1 1 0 100-2 1 1 0 000 2z"
            clipRule="evenodd"
          />
        </svg>
      )}
      {styles.label}
    </span>
  );
}

/**
 * Metric item component.
 */
function MetricItem({
  label,
  value,
  className = '',
}: {
  label: string;
  value: string | number;
  className?: string;
}) {
  return (
    <div className={className}>
      <dt className="text-xs text-slate-500 uppercase tracking-wider">{label}</dt>
      <dd className="mt-1 font-mono text-sm text-slate-800">{value}</dd>
    </div>
  );
}

/**
 * DiagnosticPanel - SCF computation diagnostic display.
 *
 * Provides clear visual feedback on computation status and educational
 * hints to help students understand their results.
 *
 * @example
 * ```tsx
 * <DiagnosticPanel
 *   converged={true}
 *   iterations={12}
 *   energy={-1.116714}
 *   convergenceProfile="medium"
 *   diisEnabled={true}
 * />
 * ```
 */
export function DiagnosticPanel({
  converged,
  iterations,
  energy,
  convergenceProfile,
  diisEnabled,
  aborted = false,
  className = '',
}: DiagnosticPanelProps) {
  const statusType = getStatusType(converged, aborted);
  const hints = getEducationalHints(
    statusType,
    iterations,
    diisEnabled,
    convergenceProfile
  );

  // Format convergence profile for display
  const profileDisplay =
    convergenceProfile.charAt(0).toUpperCase() + convergenceProfile.slice(1);

  return (
    <div
      className={`bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden ${className}`}
      role="region"
      aria-label="SCF computation diagnostics"
    >
      {/* Header with status badge */}
      <div className="px-4 py-3 border-b border-slate-200 flex items-center justify-between">
        <h3 className="text-sm font-semibold text-slate-700">Diagnostics</h3>
        <StatusBadge statusType={statusType} />
      </div>

      {/* Metrics grid */}
      <div className="px-4 py-4">
        <dl className="grid grid-cols-2 md:grid-cols-4 gap-4">
          <div>
            <dt className="text-xs text-slate-500 uppercase tracking-wider flex items-center gap-1">
              Final Energy <Math>{String.raw`E`}</Math>
            </dt>
            <dd className="mt-1 font-mono text-sm text-slate-800">{energy.toFixed(10)} Ha</dd>
          </div>
          <MetricItem label="Iterations" value={iterations} />
          <MetricItem label="Profile" value={profileDisplay} />
          <MetricItem label="DIIS" value={diisEnabled ? 'Enabled' : 'Disabled'} />
        </dl>
      </div>

      {/* Educational hints */}
      {hints.length > 0 && (
        <div className="px-4 py-3 bg-blue-50 border-t border-blue-100">
          <h4 className="text-xs font-semibold text-blue-800 uppercase tracking-wider mb-2">
            Insights
          </h4>
          <ul className="space-y-1">
            {hints.map((hint, index) => (
              <li key={index} className="text-sm text-blue-700 flex items-start gap-2">
                <svg
                  className="w-4 h-4 mt-0.5 flex-shrink-0"
                  fill="currentColor"
                  viewBox="0 0 20 20"
                  aria-hidden="true"
                >
                  <path
                    fillRule="evenodd"
                    d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z"
                    clipRule="evenodd"
                  />
                </svg>
                <span>{hint}</span>
              </li>
            ))}
          </ul>
        </div>
      )}
    </div>
  );
}

export default DiagnosticPanel;
