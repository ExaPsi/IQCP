/**
 * DownloadPDF - Button to download worksheet as PDF.
 *
 * Uses html2pdf.js to convert the rendered worksheet HTML to PDF.
 * Provides visual feedback during PDF generation.
 *
 * @module components/labs/DownloadPDF
 */

import { useState, useCallback } from 'react';

/**
 * Props for DownloadPDF component.
 */
export interface DownloadPDFProps {
  /** ID of the element to convert to PDF */
  contentId: string;
  /** Filename for the downloaded PDF (without .pdf extension) */
  filename?: string;
  /** Document title for PDF metadata */
  title?: string;
  /** Additional CSS classes */
  className?: string;
}

/**
 * DownloadPDF - Converts worksheet to PDF for offline use.
 *
 * Features:
 * - Converts rendered HTML to PDF using html2pdf.js
 * - Shows loading state during generation
 * - Configurable filename and title
 * - Print-optimized styling
 *
 * @example
 * ```tsx
 * <div id="worksheet-content">
 *   <MarkdownRenderer content={worksheetContent} />
 * </div>
 *
 * <DownloadPDF
 *   contentId="worksheet-content"
 *   filename="labpack1-worksheet"
 *   title="Lab Pack #1: From Boys Functions to SCF Convergence"
 * />
 * ```
 */
export function DownloadPDF({
  contentId,
  filename = 'worksheet',
  title: _title = 'IQCP Worksheet',
  className = '',
}: DownloadPDFProps) {
  // Note: title is reserved for future PDF metadata support
  void _title;
  const [isGenerating, setIsGenerating] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleDownload = useCallback(async () => {
    const element = document.getElementById(contentId);
    if (!element) {
      setError('Content not found');
      return;
    }

    setIsGenerating(true);
    setError(null);

    try {
      // Dynamically import html2pdf.js to avoid SSR issues
      const html2pdf = (await import('html2pdf.js')).default;

      // Configure PDF options
      const opt = {
        margin: [10, 10, 10, 10] as [number, number, number, number], // mm
        filename: `${filename}.pdf`,
        image: { type: 'jpeg' as const, quality: 0.98 },
        html2canvas: {
          scale: 2,
          useCORS: true,
          letterRendering: true,
        },
        jsPDF: {
          unit: 'mm' as const,
          format: 'a4' as const,
          orientation: 'portrait' as const,
        },
        pagebreak: { mode: ['avoid-all', 'css', 'legacy'] as const },
      };

      // Add print-specific styles temporarily
      const style = document.createElement('style');
      style.id = 'pdf-print-styles';
      style.textContent = `
        @media print {
          .no-print { display: none !important; }
          a { text-decoration: underline; color: #1e40af !important; }
          pre { white-space: pre-wrap !important; word-break: break-word; }
          table { page-break-inside: avoid; }
          h2, h3 { page-break-after: avoid; }
        }
      `;
      document.head.appendChild(style);

      // Generate PDF
      await html2pdf().set(opt).from(element).save();

      // Remove print styles
      style.remove();
    } catch (err) {
      console.error('PDF generation failed:', err);
      setError('Failed to generate PDF');
    } finally {
      setIsGenerating(false);
    }
  }, [contentId, filename]);

  return (
    <div className={`${className}`}>
      <button
        type="button"
        onClick={handleDownload}
        disabled={isGenerating}
        className="
          inline-flex items-center gap-2 px-4 py-2 rounded-lg
          text-sm font-medium transition-colors
          bg-slate-100 hover:bg-slate-200 text-slate-700
          disabled:opacity-50 disabled:cursor-not-allowed
          focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500
        "
        aria-label={isGenerating ? 'Generating PDF...' : 'Download as PDF'}
      >
        {isGenerating ? (
          // Loading spinner
          <svg
            className="w-4 h-4 animate-spin"
            fill="none"
            viewBox="0 0 24 24"
            aria-hidden="true"
          >
            <circle
              className="opacity-25"
              cx="12"
              cy="12"
              r="10"
              stroke="currentColor"
              strokeWidth="4"
            />
            <path
              className="opacity-75"
              fill="currentColor"
              d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
            />
          </svg>
        ) : (
          // Download icon
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
              d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z"
            />
          </svg>
        )}
        {isGenerating ? 'Generating...' : 'Download PDF'}
      </button>

      {error && (
        <p className="mt-1 text-xs text-red-600" role="alert">
          {error}
        </p>
      )}
    </div>
  );
}

export default DownloadPDF;
