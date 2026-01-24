/**
 * Math rendering components using KaTeX.
 * Provides inline and block math display with graceful error handling.
 */
import 'katex/dist/katex.min.css';
import katex from 'katex';
import { useMemo } from 'react';

interface MathProps {
  /** LaTeX expression to render */
  children: string;
  /** Additional CSS classes */
  className?: string;
  /** Error handler callback */
  onError?: (error: Error) => void;
}

interface MathBlockProps extends MathProps {
  /** Optional label displayed above the equation */
  label?: string;
}

/**
 * Inline math component for use within text.
 * Renders LaTeX expressions inline with surrounding content.
 *
 * @example
 * <p>The function <Math>F_m(T)</Math> is defined as...</p>
 */
export function Math({ children, className = '', onError }: MathProps) {
  const html = useMemo(() => {
    try {
      return katex.renderToString(children, {
        displayMode: false,
        throwOnError: false,
        errorColor: '#dc2626',
        strict: false,
        trust: true,
      });
    } catch (error) {
      if (onError && error instanceof Error) {
        onError(error);
      }
      // Return fallback display for parse errors
      return `<span class="text-red-600 font-mono text-sm">${children}</span>`;
    }
  }, [children, onError]);

  return (
    <span
      className={`math-inline ${className}`}
      dangerouslySetInnerHTML={{ __html: html }}
      role="math"
      aria-label={children}
    />
  );
}

/**
 * Block math component for display equations.
 * Centers the equation and renders it in display mode.
 *
 * @example
 * <MathBlock label="Boys Function Definition">
 *   F_m(T) = \int_0^1 t^{2m} e^{-Tt^2} dt
 * </MathBlock>
 */
export function MathBlock({ children, label, className = '', onError }: MathBlockProps) {
  const html = useMemo(() => {
    try {
      return katex.renderToString(children, {
        displayMode: true,
        throwOnError: false,
        errorColor: '#dc2626',
        strict: false,
        trust: true,
      });
    } catch (error) {
      if (onError && error instanceof Error) {
        onError(error);
      }
      // Return fallback display for parse errors
      return `<span class="text-red-600 font-mono text-sm block text-center">${children}</span>`;
    }
  }, [children, onError]);

  return (
    <div
      className={`bg-white rounded-lg border border-slate-200 p-4 my-4 ${className}`}
      role="math"
      aria-label={children}
    >
      {label && (
        <div className="text-xs font-medium text-slate-500 uppercase tracking-wide mb-2">
          {label}
        </div>
      )}
      <div
        className="overflow-x-auto"
        dangerouslySetInnerHTML={{ __html: html }}
      />
    </div>
  );
}

/**
 * Simple display math without the box styling.
 * Useful for equations that need to flow with content.
 */
export function MathDisplay({ children, className = '', onError }: MathProps) {
  const html = useMemo(() => {
    try {
      return katex.renderToString(children, {
        displayMode: true,
        throwOnError: false,
        errorColor: '#dc2626',
        strict: false,
        trust: true,
      });
    } catch (error) {
      if (onError && error instanceof Error) {
        onError(error);
      }
      return `<span class="text-red-600 font-mono text-sm block text-center">${children}</span>`;
    }
  }, [children, onError]);

  return (
    <div
      className={`my-3 overflow-x-auto ${className}`}
      dangerouslySetInnerHTML={{ __html: html }}
      role="math"
      aria-label={children}
    />
  );
}

export default Math;
