/**
 * TableOfContents - Navigation sidebar for worksheet sections.
 *
 * Displays a hierarchical list of sections extracted from Markdown headings.
 * Supports scroll-to-section navigation, scroll spy to highlight the current
 * visible section, and collapsible subsections.
 *
 * @module components/labs/TableOfContents
 */

import { useState, useEffect, useCallback, useMemo } from 'react';
import { type TocItem } from '@/lib/tocUtils';

/**
 * Props for TableOfContents component.
 */
export interface TableOfContentsProps {
  /** Hierarchical TOC items */
  items: TocItem[];
  /** Currently active section ID (from scroll spy) */
  activeId?: string;
  /** Callback when user clicks a TOC item */
  onNavigate?: (id: string) => void;
  /** Progress data: which sections are completed */
  completedSections?: Set<string>;
  /** Additional CSS classes */
  className?: string;
  /** Maximum height for scrollable container (default: calc(100vh - 12rem)) */
  maxHeight?: string;
}

/**
 * Props for individual TOC item.
 */
interface TocItemProps {
  item: TocItem;
  activeId?: string;
  onNavigate?: (id: string) => void;
  isCompleted?: boolean;
  depth: number;
  expandedSections: Set<string>;
  onToggleExpand: (id: string) => void;
}

/**
 * Chevron icon component for expand/collapse indicator.
 */
function ChevronIcon({ expanded }: { expanded: boolean }) {
  return (
    <svg
      className={`
        w-4 h-4 text-slate-400 transition-transform duration-200 flex-shrink-0
        ${expanded ? 'rotate-90' : 'rotate-0'}
      `}
      fill="none"
      stroke="currentColor"
      viewBox="0 0 24 24"
      aria-hidden="true"
    >
      <path
        strokeLinecap="round"
        strokeLinejoin="round"
        strokeWidth={2}
        d="M9 5l7 7-7 7"
      />
    </svg>
  );
}

/**
 * Checkmark icon for completed sections.
 */
function CheckmarkIcon() {
  return (
    <svg
      className="w-3.5 h-3.5 text-green-600 flex-shrink-0"
      fill="currentColor"
      viewBox="0 0 20 20"
      aria-hidden="true"
    >
      <path
        fillRule="evenodd"
        d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
        clipRule="evenodd"
      />
    </svg>
  );
}

/**
 * Individual TOC item with optional nesting and collapsible behavior.
 */
function TocItemComponent({
  item,
  activeId,
  onNavigate,
  isCompleted,
  depth,
  expandedSections,
  onToggleExpand,
}: TocItemProps) {
  const isActive = item.id === activeId;
  const hasChildren = item.children.length > 0;
  const isExpanded = expandedSections.has(item.id);

  // Check if any child or descendant is active
  const isChildActive = useMemo(() => {
    function checkActive(children: TocItem[]): boolean {
      for (const child of children) {
        if (child.id === activeId) return true;
        if (checkActive(child.children)) return true;
      }
      return false;
    }
    return checkActive(item.children);
  }, [item.children, activeId]);

  const handleClick = useCallback(
    (e: React.MouseEvent) => {
      e.preventDefault();
      onNavigate?.(item.id);
    },
    [item.id, onNavigate]
  );

  const handleToggle = useCallback(
    (e: React.MouseEvent) => {
      e.preventDefault();
      e.stopPropagation();
      onToggleExpand(item.id);
    },
    [item.id, onToggleExpand]
  );

  return (
    <li className={depth > 0 ? 'ml-3' : ''}>
      <div className="flex items-center gap-1">
        {/* Expand/collapse button for items with children */}
        {hasChildren ? (
          <button
            type="button"
            onClick={handleToggle}
            className="p-0.5 -ml-1 rounded hover:bg-slate-200 transition-colors"
            aria-label={isExpanded ? 'Collapse section' : 'Expand section'}
            aria-expanded={isExpanded}
          >
            <ChevronIcon expanded={isExpanded} />
          </button>
        ) : (
          <span className="w-5 flex-shrink-0" aria-hidden="true" />
        )}

        {/* TOC link */}
        <a
          href={`#${item.id}`}
          onClick={handleClick}
          className={`
            flex-1 block py-1.5 px-2 rounded text-sm transition-colors
            ${
              isActive
                ? 'bg-blue-100 text-blue-800 font-medium'
                : isChildActive
                ? 'text-blue-600'
                : 'text-slate-600 hover:text-slate-900 hover:bg-slate-100'
            }
            ${isCompleted ? 'border-l-2 border-green-500 -ml-0.5 pl-2.5' : ''}
          `}
          aria-current={isActive ? 'location' : undefined}
        >
          <span className="flex items-center gap-2">
            {isCompleted && <CheckmarkIcon />}
            <span className="truncate">{item.text}</span>
          </span>
        </a>
      </div>

      {/* Children (collapsible) */}
      {hasChildren && (
        <ul
          className={`
            mt-1 space-y-0.5 overflow-hidden transition-all duration-200
            ${isExpanded ? 'max-h-[1000px] opacity-100' : 'max-h-0 opacity-0'}
          `}
          aria-hidden={!isExpanded}
        >
          {item.children.map((child) => (
            <TocItemComponent
              key={child.id}
              item={child}
              activeId={activeId}
              onNavigate={onNavigate}
              isCompleted={false}
              depth={depth + 1}
              expandedSections={expandedSections}
              onToggleExpand={onToggleExpand}
            />
          ))}
        </ul>
      )}
    </li>
  );
}

/**
 * useScrollSpy - Hook to track which section is currently visible.
 *
 * Uses IntersectionObserver to detect when headings enter/leave viewport.
 *
 * @param sectionIds - Array of section IDs to observe
 * @returns Currently active section ID
 */
export function useScrollSpy(sectionIds: string[]): string | undefined {
  const [activeId, setActiveId] = useState<string | undefined>();

  useEffect(() => {
    if (sectionIds.length === 0) return;

    const observer = new IntersectionObserver(
      (entries) => {
        // Find the first visible section
        for (const entry of entries) {
          if (entry.isIntersecting) {
            setActiveId(entry.target.id);
            break;
          }
        }
      },
      {
        rootMargin: '-80px 0px -70% 0px', // Trigger when heading is near top
        threshold: 0,
      }
    );

    // Observe all sections
    for (const id of sectionIds) {
      const element = document.getElementById(id);
      if (element) {
        observer.observe(element);
      }
    }

    return () => observer.disconnect();
  }, [sectionIds]);

  return activeId;
}

/**
 * Scroll position indicator showing current position in document.
 */
function ScrollIndicator({ percent }: { percent: number }) {
  return (
    <div className="mb-3">
      <div className="flex items-center justify-between text-xs text-slate-500 mb-1">
        <span>Reading progress</span>
        <span>{Math.round(percent)}%</span>
      </div>
      <div className="h-1 bg-slate-200 rounded-full overflow-hidden">
        <div
          className="h-full bg-blue-500 transition-all duration-150"
          style={{ width: `${percent}%` }}
        />
      </div>
    </div>
  );
}

/**
 * useScrollProgress - Hook to track scroll progress through the document.
 *
 * @returns Current scroll progress as percentage (0-100)
 */
export function useScrollProgress(): number {
  const [progress, setProgress] = useState(0);

  useEffect(() => {
    const updateProgress = () => {
      const scrollTop = window.scrollY;
      const docHeight = document.documentElement.scrollHeight - window.innerHeight;
      const percent = docHeight > 0 ? (scrollTop / docHeight) * 100 : 0;
      setProgress(Math.min(100, Math.max(0, percent)));
    };

    // Initial calculation
    updateProgress();

    // Throttled scroll handler
    let ticking = false;
    const handleScroll = () => {
      if (!ticking) {
        window.requestAnimationFrame(() => {
          updateProgress();
          ticking = false;
        });
        ticking = true;
      }
    };

    window.addEventListener('scroll', handleScroll, { passive: true });
    return () => window.removeEventListener('scroll', handleScroll);
  }, []);

  return progress;
}

/**
 * TableOfContents - Sidebar navigation for worksheet sections.
 *
 * Features:
 * - Hierarchical section display with collapsible subsections
 * - Smooth scroll-to-section navigation
 * - Scroll spy highlighting current section
 * - Progress indicators for completed sections
 * - Independent scrolling within the TOC container
 * - Scroll position indicator
 * - Collapsible on mobile (controlled externally)
 *
 * @example
 * ```tsx
 * const toc = buildToc(worksheetContent);
 * const completedSections = new Set(['introduction', 'section-1']);
 *
 * <TableOfContents
 *   items={toc}
 *   completedSections={completedSections}
 *   onNavigate={(id) => {
 *     document.getElementById(id)?.scrollIntoView({ behavior: 'smooth' });
 *   }}
 * />
 * ```
 */
export function TableOfContents({
  items,
  activeId,
  onNavigate,
  completedSections = new Set(),
  className = '',
  maxHeight = 'calc(100vh - 12rem)',
}: TableOfContentsProps) {
  // Track which sections are expanded (default: all expanded)
  const [expandedSections, setExpandedSections] = useState<Set<string>>(() => {
    // Initialize with all parent sections expanded
    const expanded = new Set<string>();
    function addParents(tocItems: TocItem[]) {
      for (const item of tocItems) {
        if (item.children.length > 0) {
          expanded.add(item.id);
          addParents(item.children);
        }
      }
    }
    addParents(items);
    return expanded;
  });

  // Scroll progress for indicator
  const scrollProgress = useScrollProgress();

  // Toggle section expansion
  const handleToggleExpand = useCallback((id: string) => {
    setExpandedSections((prev) => {
      const next = new Set(prev);
      if (next.has(id)) {
        next.delete(id);
      } else {
        next.add(id);
      }
      return next;
    });
  }, []);

  // Auto-expand section when active item is inside a collapsed section
  useEffect(() => {
    if (!activeId) return;

    // Find parent sections that contain the active item
    function findParents(
      tocItems: TocItem[],
      targetId: string,
      path: string[] = []
    ): string[] | null {
      for (const item of tocItems) {
        if (item.id === targetId) {
          return path;
        }
        if (item.children.length > 0) {
          const result = findParents(item.children, targetId, [...path, item.id]);
          if (result) return result;
        }
      }
      return null;
    }

    const parents = findParents(items, activeId);
    if (parents && parents.length > 0) {
      setExpandedSections((prev) => {
        const next = new Set(prev);
        for (const parentId of parents) {
          next.add(parentId);
        }
        return next;
      });
    }
  }, [activeId, items]);

  const handleNavigate = useCallback(
    (id: string) => {
      // Smooth scroll to section
      const element = document.getElementById(id);
      if (element) {
        element.scrollIntoView({ behavior: 'smooth', block: 'start' });
      }
      onNavigate?.(id);
    },
    [onNavigate]
  );

  if (items.length === 0) {
    return null;
  }

  return (
    <nav
      className={`${className}`}
      aria-label="Table of contents"
    >
      {/* Scroll progress indicator */}
      <ScrollIndicator percent={scrollProgress} />

      <h2 className="text-sm font-semibold text-slate-900 uppercase tracking-wide mb-3">
        Contents
      </h2>

      {/* Scrollable TOC container */}
      <div
        className="overflow-y-auto overflow-x-hidden scrollbar-thin"
        style={{ maxHeight }}
      >
        <ul className="space-y-0.5 pr-2">
          {items.map((item) => (
            <TocItemComponent
              key={item.id}
              item={item}
              activeId={activeId}
              onNavigate={handleNavigate}
              isCompleted={completedSections.has(item.id)}
              depth={0}
              expandedSections={expandedSections}
              onToggleExpand={handleToggleExpand}
            />
          ))}
        </ul>
      </div>
    </nav>
  );
}

export default TableOfContents;
