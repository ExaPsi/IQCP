/**
 * TOC Utilities - Extract and manage table of contents from Markdown.
 *
 * Provides functions to parse headings from Markdown content
 * and build a hierarchical table of contents structure.
 *
 * @module lib/tocUtils
 */

import { slugify } from '@/components/labs/MarkdownRenderer';

/**
 * Represents a single TOC item (heading).
 */
export interface TocItem {
  /** Unique ID for the heading (slugified text) */
  id: string;
  /** Original heading text */
  text: string;
  /** Heading level (1-6) */
  level: number;
  /** Child headings (nested under this heading) */
  children: TocItem[];
}

/**
 * Pattern to match Markdown headings.
 * Captures: level (# count), text (heading content)
 */
const HEADING_PATTERN = /^(#{1,6})\s+(.+)$/gm;

/**
 * Extract headings from Markdown content.
 *
 * @param content - Raw Markdown string
 * @returns Flat array of heading objects
 */
export function extractHeadings(
  content: string
): { level: number; text: string; id: string }[] {
  const headings: { level: number; text: string; id: string }[] = [];
  let match;

  // Reset regex lastIndex
  HEADING_PATTERN.lastIndex = 0;

  while ((match = HEADING_PATTERN.exec(content)) !== null) {
    const level = match[1].length;
    const text = match[2].trim();
    const id = slugify(text);

    headings.push({ level, text, id });
  }

  return headings;
}

/**
 * Build a hierarchical TOC from flat headings.
 *
 * Groups headings based on their levels, creating a nested structure.
 * Only includes headings at level 2 and 3 for cleaner navigation.
 *
 * @param content - Raw Markdown string
 * @param maxLevel - Maximum heading level to include (default: 3)
 * @returns Array of top-level TOC items with nested children
 */
export function buildToc(content: string, maxLevel: number = 3): TocItem[] {
  const headings = extractHeadings(content).filter(
    (h) => h.level >= 2 && h.level <= maxLevel
  );

  const toc: TocItem[] = [];
  const stack: TocItem[] = [];

  for (const heading of headings) {
    const item: TocItem = {
      id: heading.id,
      text: heading.text,
      level: heading.level,
      children: [],
    };

    // Find parent for this heading
    while (stack.length > 0 && stack[stack.length - 1].level >= heading.level) {
      stack.pop();
    }

    if (stack.length === 0) {
      // Top-level heading
      toc.push(item);
    } else {
      // Nested under parent
      stack[stack.length - 1].children.push(item);
    }

    stack.push(item);
  }

  return toc;
}

/**
 * Flatten a hierarchical TOC back to a linear list.
 *
 * Useful for progress tracking and indexing.
 *
 * @param toc - Hierarchical TOC items
 * @returns Flat array of all items
 */
export function flattenToc(toc: TocItem[]): TocItem[] {
  const flat: TocItem[] = [];

  function traverse(items: TocItem[]) {
    for (const item of items) {
      flat.push(item);
      if (item.children.length > 0) {
        traverse(item.children);
      }
    }
  }

  traverse(toc);
  return flat;
}

/**
 * Get section IDs at a specific level.
 *
 * Used for progress tracking where we track completion
 * of major sections (level 2 headings).
 *
 * @param content - Raw Markdown string
 * @param level - Heading level to extract (default: 2)
 * @returns Array of section IDs
 */
export function getSectionIds(content: string, level: number = 2): string[] {
  return extractHeadings(content)
    .filter((h) => h.level === level)
    .map((h) => h.id);
}
