/**
 * WorksheetViewer - Container component for displaying Lab Pack worksheets.
 *
 * Combines the Markdown renderer, table of contents, progress tracking,
 * and PDF download into a responsive layout with independently scrollable
 * sidebar navigation.
 *
 * @module components/labs/WorksheetViewer
 */

import { useState, useMemo, useCallback } from 'react';
import { MarkdownRenderer } from './MarkdownRenderer';
import { TableOfContents, useScrollSpy } from './TableOfContents';
import { ProgressTracker } from './ProgressTracker';
import { DownloadPDF } from './DownloadPDF';
import { useWorksheetProgress } from '@/hooks/useWorksheetProgress';
import { buildToc, flattenToc, getSectionIds } from '@/lib/tocUtils';

/**
 * Props for WorksheetViewer component.
 */
export interface WorksheetViewerProps {
  /** Lab pack identifier (e.g., 'labpack1') */
  labpack: string;
  /** Worksheet title */
  title: string;
  /** Raw Markdown content */
  content: string;
  /** PDF filename (without extension) */
  pdfFilename?: string;
}

/**
 * WorksheetViewer - Full-featured worksheet display.
 *
 * Features:
 * - Responsive layout (sidebar TOC on desktop, collapsed on mobile)
 * - Markdown rendering with deep link support
 * - Table of contents with scroll spy and collapsible sections
 * - Progress tracking with localStorage persistence
 * - PDF download option
 * - Independently scrollable TOC sidebar
 *
 * @example
 * ```tsx
 * import worksheetContent from '../../../content/labpack1/worksheet-student.md?raw';
 *
 * <WorksheetViewer
 *   labpack="labpack1"
 *   title="Lab Pack #1: From Boys Functions to SCF Convergence"
 *   content={worksheetContent}
 *   pdfFilename="labpack1-worksheet"
 * />
 * ```
 */
export function WorksheetViewer({
  labpack,
  title,
  content,
  pdfFilename = 'worksheet',
}: WorksheetViewerProps) {
  // Mobile TOC state
  const [isTocOpen, setIsTocOpen] = useState(false);

  // Build TOC from content
  const toc = useMemo(() => buildToc(content, 3), [content]);
  const flatToc = useMemo(() => flattenToc(toc), [toc]);
  const sectionIds = useMemo(() => flatToc.map((item) => item.id), [flatToc]);

  // Get main section IDs (level 2 headings) for progress tracking
  const mainSectionIds = useMemo(() => getSectionIds(content, 2), [content]);

  // Scroll spy for active section
  const activeId = useScrollSpy(sectionIds);

  // Progress tracking
  const {
    completedSections,
    toggleCompleted,
    clearProgress,
    progressPercent,
    setLastVisited,
  } = useWorksheetProgress(labpack, mainSectionIds.length);

  // Handle TOC navigation
  const handleNavigate = useCallback(
    (id: string) => {
      setLastVisited(id);
      setIsTocOpen(false); // Close mobile TOC after navigation
    },
    [setLastVisited]
  );

  // Toggle mobile TOC
  const toggleToc = useCallback(() => {
    setIsTocOpen((prev) => !prev);
  }, []);

  return (
    <div className="relative">
      {/* Header with title and actions */}
      <div className="sticky top-0 z-20 bg-white border-b border-slate-200 px-4 py-3 -mx-4 mb-6">
        <div className="max-w-7xl mx-auto flex items-center justify-between gap-4">
          {/* Mobile TOC toggle */}
          <button
            type="button"
            onClick={toggleToc}
            className="lg:hidden p-2 -ml-2 rounded-lg hover:bg-slate-100 transition-colors"
            aria-label={isTocOpen ? 'Close table of contents' : 'Open table of contents'}
            aria-expanded={isTocOpen}
          >
            <svg
              className="w-5 h-5 text-slate-600"
              fill="none"
              stroke="currentColor"
              viewBox="0 0 24 24"
              aria-hidden="true"
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeWidth={2}
                d="M4 6h16M4 12h16M4 18h16"
              />
            </svg>
          </button>

          {/* Title (hidden on mobile) */}
          <h1 className="hidden sm:block text-lg font-semibold text-slate-800 truncate">
            {title}
          </h1>

          {/* Actions */}
          <div className="flex items-center gap-2">
            <DownloadPDF
              contentId="worksheet-content"
              filename={pdfFilename}
              title={title}
              className="no-print"
            />
          </div>
        </div>
      </div>

      {/* Main layout */}
      <div className="flex gap-8">
        {/* Desktop TOC sidebar - sticky with independent scrolling */}
        <aside className="hidden lg:block w-64 flex-shrink-0">
          <div
            className="sticky top-20"
            style={{
              maxHeight: 'calc(100vh - 6rem)',
              overflowY: 'auto',
              overflowX: 'hidden',
            }}
          >
            {/* Progress tracker */}
            <ProgressTracker
              percent={progressPercent}
              completedCount={completedSections.size}
              totalCount={mainSectionIds.length}
              onClearProgress={clearProgress}
              className="mb-4 p-4 bg-slate-50 rounded-lg"
            />

            {/* Table of contents with internal scrolling */}
            <TableOfContents
              items={toc}
              activeId={activeId}
              onNavigate={handleNavigate}
              completedSections={completedSections}
              className="p-4 bg-slate-50 rounded-lg"
              maxHeight="calc(100vh - 22rem)"
            />

            {/* Section completion checkboxes */}
            <div className="mt-4 p-4 bg-slate-50 rounded-lg">
              <h3 className="text-sm font-semibold text-slate-700 mb-3">
                Mark Sections Complete
              </h3>
              <ul className="space-y-2 max-h-48 overflow-y-auto">
                {mainSectionIds.map((id) => {
                  const item = flatToc.find((t) => t.id === id);
                  return (
                    <li key={id} className="flex items-center gap-2">
                      <input
                        type="checkbox"
                        id={`complete-${id}`}
                        checked={completedSections.has(id)}
                        onChange={() => toggleCompleted(id)}
                        className="
                          h-4 w-4 rounded border-slate-300
                          text-blue-600 focus:ring-blue-500
                          cursor-pointer
                        "
                      />
                      <label
                        htmlFor={`complete-${id}`}
                        className="text-sm text-slate-600 cursor-pointer truncate"
                      >
                        {item?.text ?? id}
                      </label>
                    </li>
                  );
                })}
              </ul>
            </div>
          </div>
        </aside>

        {/* Mobile TOC overlay */}
        {isTocOpen && (
          <>
            {/* Backdrop */}
            <div
              className="fixed inset-0 z-30 bg-black/50 lg:hidden"
              onClick={() => setIsTocOpen(false)}
              aria-hidden="true"
            />

            {/* Slide-out panel */}
            <aside
              className="
                fixed top-0 left-0 z-40 w-80 h-full
                bg-white shadow-xl overflow-y-auto
                lg:hidden
              "
            >
              <div className="p-4">
                <div className="flex items-center justify-between mb-4">
                  <h2 className="text-lg font-semibold text-slate-800">
                    Contents
                  </h2>
                  <button
                    type="button"
                    onClick={() => setIsTocOpen(false)}
                    className="p-2 -mr-2 rounded-lg hover:bg-slate-100 transition-colors"
                    aria-label="Close table of contents"
                  >
                    <svg
                      className="w-5 h-5 text-slate-600"
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
                  </button>
                </div>

                {/* Progress tracker */}
                <ProgressTracker
                  percent={progressPercent}
                  completedCount={completedSections.size}
                  totalCount={mainSectionIds.length}
                  onClearProgress={clearProgress}
                  className="mb-6 p-3 bg-slate-50 rounded-lg"
                />

                {/* Table of contents */}
                <TableOfContents
                  items={toc}
                  activeId={activeId}
                  onNavigate={handleNavigate}
                  completedSections={completedSections}
                  maxHeight="calc(100vh - 16rem)"
                />

                {/* Section completion checkboxes (mobile) */}
                <div className="mt-6 p-3 bg-slate-50 rounded-lg">
                  <h3 className="text-sm font-semibold text-slate-700 mb-3">
                    Mark Sections Complete
                  </h3>
                  <ul className="space-y-2">
                    {mainSectionIds.map((id) => {
                      const item = flatToc.find((t) => t.id === id);
                      return (
                        <li key={id} className="flex items-center gap-2">
                          <input
                            type="checkbox"
                            id={`mobile-complete-${id}`}
                            checked={completedSections.has(id)}
                            onChange={() => toggleCompleted(id)}
                            className="
                              h-4 w-4 rounded border-slate-300
                              text-blue-600 focus:ring-blue-500
                              cursor-pointer
                            "
                          />
                          <label
                            htmlFor={`mobile-complete-${id}`}
                            className="text-sm text-slate-600 cursor-pointer truncate"
                          >
                            {item?.text ?? id}
                          </label>
                        </li>
                      );
                    })}
                  </ul>
                </div>
              </div>
            </aside>
          </>
        )}

        {/* Main content */}
        <main className="flex-1 min-w-0">
          <article id="worksheet-content">
            <MarkdownRenderer content={content} />
          </article>
        </main>
      </div>
    </div>
  );
}

export default WorksheetViewer;
