/**
 * Labs page - Displays Lab Pack worksheets for guided laboratory exercises.
 *
 * Features Lab Pack #1: "From Boys Functions to SCF Convergence" with
 * interactive deep links, progress tracking, and PDF export.
 *
 * @module pages/Labs
 */

import { WorksheetViewer } from '@/components/labs';

// Import worksheet content as raw string
import worksheetContent from '../../../../content/labpack1/worksheet-student.md?raw';

/**
 * Labs - Page component for Lab Pack materials.
 *
 * Displays the student worksheet with:
 * - Table of contents navigation
 * - Deep links to module states
 * - Progress tracking (localStorage)
 * - PDF download option
 *
 * @example
 * ```tsx
 * // In router configuration:
 * { path: '/labs', element: <Labs /> }
 * ```
 */
function Labs() {
  return (
    <div className="max-w-7xl mx-auto px-4">
      <WorksheetViewer
        labpack="labpack1"
        title="Lab Pack #1: From Boys Functions to SCF Convergence"
        content={worksheetContent}
        pdfFilename="labpack1-worksheet"
      />
    </div>
  );
}

export default Labs;
