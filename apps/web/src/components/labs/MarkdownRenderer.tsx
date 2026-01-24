/**
 * MarkdownRenderer - Renders Markdown content with custom components.
 *
 * Uses react-markdown with remark-gfm for GitHub-flavored Markdown support
 * (tables, checkboxes, strikethrough). Includes custom heading renderer
 * for TOC anchors and DeepLinkAnchor for IQCP deep links.
 *
 * Math rendering is supported via KaTeX:
 * - Inline math: $F_m(T)$ or \(F_m(T)\)
 * - Display math: $$F_m(T) = \int_0^1 t^{2m} e^{-Tt^2} dt$$ or \[...\]
 *
 * @module components/labs/MarkdownRenderer
 */

import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import remarkMath from 'remark-math';
import rehypeRaw from 'rehype-raw';
import rehypeKatex from 'rehype-katex';
import { type ReactNode, useCallback } from 'react';
import { DeepLinkAnchor } from './DeepLinkAnchor';

// Import KaTeX CSS for proper math rendering
import 'katex/dist/katex.min.css';

/**
 * Generate a slug from heading text for anchor links.
 */
export function slugify(text: string): string {
  return text
    .toLowerCase()
    .trim()
    .replace(/[^\w\s-]/g, '') // Remove special characters
    .replace(/\s+/g, '-') // Replace spaces with hyphens
    .replace(/-+/g, '-'); // Remove consecutive hyphens
}

/**
 * Props for MarkdownRenderer component.
 */
export interface MarkdownRendererProps {
  /** Markdown content to render */
  content: string;
  /** Additional CSS classes for the container */
  className?: string;
  /** Callback when a checkbox is toggled */
  onCheckboxToggle?: (index: number, checked: boolean) => void;
}

/**
 * Custom heading component with anchor IDs for TOC navigation.
 */
function createHeadingComponent(level: 1 | 2 | 3 | 4 | 5 | 6) {
  const baseStyles: Record<number, string> = {
    1: 'text-3xl font-bold text-slate-900 mt-8 mb-4 pb-2 border-b border-slate-200',
    2: 'text-2xl font-bold text-slate-800 mt-8 mb-3 pb-2 border-b border-slate-200',
    3: 'text-xl font-semibold text-slate-800 mt-6 mb-2',
    4: 'text-lg font-semibold text-slate-700 mt-4 mb-2',
    5: 'text-base font-semibold text-slate-700 mt-3 mb-1',
    6: 'text-sm font-semibold text-slate-600 mt-2 mb-1',
  };

  return function HeadingComponent({
    children,
  }: {
    children?: ReactNode;
  }) {
    // Extract text content for slug generation
    const text =
      typeof children === 'string'
        ? children
        : Array.isArray(children)
          ? children
              .map((child) => (typeof child === 'string' ? child : ''))
              .join('')
          : '';
    const id = slugify(text);

    // Use explicit JSX elements for proper typing
    const className = baseStyles[level];
    switch (level) {
      case 1:
        return <h1 id={id} className={className}>{children}</h1>;
      case 2:
        return <h2 id={id} className={className}>{children}</h2>;
      case 3:
        return <h3 id={id} className={className}>{children}</h3>;
      case 4:
        return <h4 id={id} className={className}>{children}</h4>;
      case 5:
        return <h5 id={id} className={className}>{children}</h5>;
      case 6:
        return <h6 id={id} className={className}>{children}</h6>;
    }
  };
}

/**
 * MarkdownRenderer - Renders Markdown with custom styling and components.
 *
 * Features:
 * - GitHub-flavored Markdown (tables, checkboxes, strikethrough)
 * - Custom heading components with anchor IDs
 * - Deep link handling via DeepLinkAnchor
 * - Styled code blocks and tables
 * - Accessible and responsive design
 *
 * @example
 * ```tsx
 * const worksheetContent = `# Lab Pack #1\n\nWelcome to the lab...`;
 *
 * <MarkdownRenderer content={worksheetContent} />
 * ```
 */
export function MarkdownRenderer({
  content,
  className = '',
  onCheckboxToggle,
}: MarkdownRendererProps) {
  // Track checkbox indices for toggle callbacks
  let checkboxIndex = 0;

  const handleCheckboxChange = useCallback(
    (index: number) => (e: React.ChangeEvent<HTMLInputElement>) => {
      onCheckboxToggle?.(index, e.target.checked);
    },
    [onCheckboxToggle]
  );

  return (
    <div
      className={`
        prose prose-slate max-w-none
        prose-headings:scroll-mt-20
        prose-a:no-underline
        prose-code:before:content-none prose-code:after:content-none
        prose-pre:bg-slate-900 prose-pre:text-slate-100
        ${className}
      `}
    >
      <ReactMarkdown
        remarkPlugins={[remarkGfm, remarkMath]}
        rehypePlugins={[rehypeRaw, rehypeKatex]}
        components={{
          // Custom headings with anchor IDs
          h1: createHeadingComponent(1),
          h2: createHeadingComponent(2),
          h3: createHeadingComponent(3),
          h4: createHeadingComponent(4),
          h5: createHeadingComponent(5),
          h6: createHeadingComponent(6),

          // Custom anchor with deep link support
          a: ({ href, children, ...props }) => (
            <DeepLinkAnchor href={href} {...props}>
              {children}
            </DeepLinkAnchor>
          ),

          // Styled paragraphs
          p: ({ children }) => (
            <p className="text-slate-700 leading-relaxed mb-4">{children}</p>
          ),

          // Styled lists
          ul: ({ children }) => (
            <ul className="list-disc list-outside ml-6 mb-4 space-y-1 text-slate-700">
              {children}
            </ul>
          ),
          ol: ({ children }) => (
            <ol className="list-decimal list-outside ml-6 mb-4 space-y-1 text-slate-700">
              {children}
            </ol>
          ),
          li: ({ children, ...props }) => {
            // Check if this is a task list item (has checkbox)
            const childArray = Array.isArray(children) ? children : [children];
            const hasCheckbox = childArray.some(
              (child) =>
                child &&
                typeof child === 'object' &&
                'type' in child &&
                child.type === 'input'
            );

            if (hasCheckbox) {
              return (
                <li className="list-none -ml-6 flex items-start gap-2" {...props}>
                  {children}
                </li>
              );
            }

            return (
              <li className="text-slate-700" {...props}>
                {children}
              </li>
            );
          },

          // Task list checkboxes
          input: ({ type, checked, ...props }) => {
            if (type === 'checkbox') {
              const currentIndex = checkboxIndex++;
              return (
                <input
                  type="checkbox"
                  checked={checked}
                  onChange={handleCheckboxChange(currentIndex)}
                  className="
                    mt-1 h-4 w-4 rounded border-slate-300
                    text-blue-600 focus:ring-blue-500
                    cursor-pointer
                  "
                  {...props}
                />
              );
            }
            return <input type={type} {...props} />;
          },

          // Styled blockquotes
          blockquote: ({ children }) => (
            <blockquote
              className="
                border-l-4 border-blue-500 bg-blue-50
                pl-4 py-2 my-4 italic text-slate-600
              "
            >
              {children}
            </blockquote>
          ),

          // Styled code blocks
          pre: ({ children }) => (
            <pre
              className="
                bg-slate-900 text-slate-100
                rounded-lg p-4 overflow-x-auto
                text-sm font-mono my-4
              "
            >
              {children}
            </pre>
          ),
          code: ({ children, className }) => {
            // Check if this is inline code or a code block
            const isBlock = className?.includes('language-');

            if (isBlock) {
              return <code className={className}>{children}</code>;
            }

            return (
              <code
                className="
                  bg-slate-100 text-slate-800
                  px-1.5 py-0.5 rounded text-sm
                  font-mono
                "
              >
                {children}
              </code>
            );
          },

          // Styled tables
          table: ({ children }) => (
            <div className="overflow-x-auto my-4">
              <table className="min-w-full divide-y divide-slate-200 border border-slate-200 rounded-lg">
                {children}
              </table>
            </div>
          ),
          thead: ({ children }) => (
            <thead className="bg-slate-50">{children}</thead>
          ),
          tbody: ({ children }) => (
            <tbody className="divide-y divide-slate-200 bg-white">
              {children}
            </tbody>
          ),
          tr: ({ children }) => <tr>{children}</tr>,
          th: ({ children }) => (
            <th
              className="
                px-4 py-2 text-left text-sm font-semibold
                text-slate-700 border-b border-slate-200
              "
            >
              {children}
            </th>
          ),
          td: ({ children }) => (
            <td className="px-4 py-2 text-sm text-slate-700">{children}</td>
          ),

          // Horizontal rules
          hr: () => <hr className="my-8 border-slate-300" />,

          // Strong and emphasis
          strong: ({ children }) => (
            <strong className="font-semibold text-slate-900">{children}</strong>
          ),
          em: ({ children }) => (
            <em className="italic text-slate-700">{children}</em>
          ),
        }}
      >
        {content}
      </ReactMarkdown>
    </div>
  );
}

export default MarkdownRenderer;
