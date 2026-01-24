/**
 * ImportButton - Button to import computation results from JSON artifact.
 *
 * Used for loading run artifacts containing state, results, and metadata
 * to restore previous sessions without re-computation.
 * Shows "Imported!" feedback after successful import.
 *
 * @module components/common/ImportButton
 */

import { useState, useCallback, useRef } from 'react';
import { validateRunArtifact, type RunArtifactV1 } from '../../types/run-artifact';

/**
 * Import error types for error handling.
 */
export type ImportErrorType = 'file_read' | 'json_parse' | 'validation' | 'version_mismatch';

/**
 * Import error object for error callbacks.
 */
export interface ImportError {
  type: ImportErrorType;
  message: string;
  details?: unknown;
}

/**
 * Props for ImportButton component.
 */
export interface ImportButtonProps {
  /** Callback when artifact is successfully imported */
  onImport: (artifact: RunArtifactV1) => void;
  /** Callback when import fails */
  onError?: (error: ImportError) => void;
  /** Whether the button should be disabled */
  disabled?: boolean;
  /** Additional CSS classes */
  className?: string;
  /** Button label (default: "Import Results") */
  label?: string;
}

/**
 * Duration to show feedback state (ms).
 */
const FEEDBACK_DURATION_MS = 2000;

/**
 * Read a file as text using the native FileReader API.
 *
 * @param file - File to read
 * @returns Promise resolving to file contents as string
 */
function readFileAsText(file: File): Promise<string> {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(reader.result as string);
    reader.onerror = () => reject(reader.error);
    reader.readAsText(file);
  });
}

/**
 * ImportButton - Loads artifact from JSON file with visual feedback.
 *
 * Features:
 * - Upload icon (SVG)
 * - Hidden file input triggered by button click
 * - Shows "Imported!" feedback for 2 seconds after success
 * - Shows "Failed" feedback on error
 * - Validates artifact against Zod schema
 * - Warns on version mismatch but attempts import
 * - Uses same button style as ExportButton for consistency
 * - Accessible with proper ARIA attributes
 *
 * @example
 * ```tsx
 * // Basic usage
 * <ImportButton
 *   onImport={handleImport}
 *   disabled={!isReady}
 * />
 *
 * // With error callback
 * <ImportButton
 *   onImport={handleImport}
 *   onError={(err) => console.error(err.message)}
 *   label="Load Session"
 * />
 * ```
 */
export function ImportButton({
  onImport,
  onError,
  disabled = false,
  className = '',
  label = 'Import Results',
}: ImportButtonProps) {
  const [imported, setImported] = useState(false);
  const [error, setError] = useState(false);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleError = useCallback(
    (type: ImportErrorType, message: string, details?: unknown) => {
      console.error(`Import failed (${type}):`, message, details);
      setError(true);
      setErrorMessage(message);

      // Call error callback if provided
      onError?.({ type, message, details });

      // Reset error state after feedback duration
      setTimeout(() => {
        setError(false);
        setErrorMessage(null);
      }, FEEDBACK_DURATION_MS);
    },
    [onError]
  );

  const handleFileChange = useCallback(
    async (event: React.ChangeEvent<HTMLInputElement>) => {
      const file = event.target.files?.[0];
      if (!file) return;

      // Reset input value so the same file can be re-selected
      event.target.value = '';

      try {
        // 1. Read file as text
        let contents: string;
        try {
          contents = await readFileAsText(file);
        } catch (readError) {
          handleError('file_read', 'Failed to read file', readError);
          return;
        }

        // 2. Parse JSON
        let data: unknown;
        try {
          data = JSON.parse(contents);
        } catch (parseError) {
          handleError('json_parse', 'Invalid JSON format', parseError);
          return;
        }

        // 3. Check version (warn but continue)
        const rawData = data as Record<string, unknown>;
        if (rawData.schema_version && rawData.schema_version !== 'run_artifact_v1') {
          console.warn(
            `Artifact version mismatch: expected 'run_artifact_v1', got '${rawData.schema_version}'. Attempting import anyway.`
          );
        }

        // 4. Validate against Zod schema
        const validationResult = validateRunArtifact(data);
        if (!validationResult.success) {
          const formattedError = validationResult.error.format();
          handleError(
            'validation',
            'Invalid artifact structure',
            formattedError
          );
          return;
        }

        // 5. Success - call import callback
        const artifact = validationResult.data;
        onImport(artifact);

        // Show success feedback
        setImported(true);
        setError(false);

        // Reset after feedback duration
        setTimeout(() => {
          setImported(false);
        }, FEEDBACK_DURATION_MS);
      } catch (err) {
        // Catch-all for unexpected errors
        handleError('validation', 'Unexpected error during import', err);
      }
    },
    [onImport, handleError]
  );

  const handleClick = useCallback(() => {
    fileInputRef.current?.click();
  }, []);

  // Determine button text and style based on state
  const buttonText = imported ? 'Imported!' : error ? 'Failed' : label;
  const buttonStyle = imported
    ? 'bg-green-600 hover:bg-green-700 text-white'
    : error
      ? 'bg-red-600 hover:bg-red-700 text-white'
      : 'bg-slate-100 hover:bg-slate-200 text-slate-700';

  return (
    <>
      {/* Hidden file input */}
      <input
        ref={fileInputRef}
        type="file"
        accept=".json,application/json"
        onChange={handleFileChange}
        className="hidden"
        aria-hidden="true"
      />

      {/* Visible button */}
      <button
        type="button"
        onClick={handleClick}
        disabled={disabled || imported}
        className={`
          inline-flex items-center gap-2 px-4 py-2 rounded-lg
          text-sm font-medium transition-colors
          focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500
          disabled:opacity-50 disabled:cursor-not-allowed
          ${buttonStyle}
          ${className}
        `}
        aria-label={
          imported
            ? 'Results imported successfully'
            : error && errorMessage
              ? `Import failed: ${errorMessage}`
              : 'Import results from JSON file'
        }
        title={error && errorMessage ? errorMessage : undefined}
      >
        {/* Icon */}
        {imported ? (
          // Checkmark icon
          <svg
            className="w-4 h-4"
            fill="none"
            stroke="currentColor"
            viewBox="0 0 24 24"
            aria-hidden="true"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M5 13l4 4L19 7"
            />
          </svg>
        ) : error ? (
          // X icon
          <svg
            className="w-4 h-4"
            fill="none"
            stroke="currentColor"
            viewBox="0 0 24 24"
            aria-hidden="true"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M6 18L18 6M6 6l12 12"
            />
          </svg>
        ) : (
          // Upload icon
          <svg
            className="w-4 h-4"
            fill="none"
            stroke="currentColor"
            viewBox="0 0 24 24"
            aria-hidden="true"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-8l-4-4m0 0L8 8m4-4v12"
            />
          </svg>
        )}
        {buttonText}
      </button>
    </>
  );
}

export default ImportButton;
