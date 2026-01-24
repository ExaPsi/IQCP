/**
 * DeepLinkAnchor - Custom anchor component for handling IQCP deep links.
 *
 * Intercepts links to iqcp.dev and uses client-side routing instead of
 * full page navigation. External links open in new tabs.
 *
 * @module components/labs/DeepLinkAnchor
 */

import { useCallback, type ReactNode, type AnchorHTMLAttributes } from 'react';
import { useNavigate } from 'react-router-dom';

/**
 * Production domain for IQCP deep links.
 */
const IQCP_DOMAIN = 'iqcp.dev';

/**
 * Props for DeepLinkAnchor component.
 */
export interface DeepLinkAnchorProps
  extends AnchorHTMLAttributes<HTMLAnchorElement> {
  /** Link href */
  href?: string;
  /** Link content */
  children?: ReactNode;
}

/**
 * Check if a URL is an IQCP deep link (points to iqcp.dev).
 */
function isIqcpDeepLink(href: string): boolean {
  try {
    const url = new URL(href);
    return url.hostname === IQCP_DOMAIN || url.hostname === `www.${IQCP_DOMAIN}`;
  } catch {
    return false;
  }
}

/**
 * Extract the path and query string from a URL.
 */
function extractPath(href: string): string {
  try {
    const url = new URL(href);
    return url.pathname + url.search;
  } catch {
    return href;
  }
}

/**
 * DeepLinkAnchor - Handles deep links with client-side routing.
 *
 * For IQCP deep links (iqcp.dev/*):
 * - Extracts path and query params
 * - Uses react-router navigate() for client-side navigation
 * - No page reload required
 *
 * For external links:
 * - Opens in new tab with proper security attributes
 *
 * @example
 * ```tsx
 * // IQCP deep link - navigates client-side
 * <DeepLinkAnchor href="https://iqcp.dev/boys?s=...">
 *   Open Boys Module
 * </DeepLinkAnchor>
 *
 * // External link - opens in new tab
 * <DeepLinkAnchor href="https://example.com">
 *   External Resource
 * </DeepLinkAnchor>
 * ```
 */
export function DeepLinkAnchor({
  href,
  children,
  className,
  ...rest
}: DeepLinkAnchorProps) {
  const navigate = useNavigate();

  const handleClick = useCallback(
    (e: React.MouseEvent<HTMLAnchorElement>) => {
      if (!href) return;

      if (isIqcpDeepLink(href)) {
        e.preventDefault();
        const path = extractPath(href);
        navigate(path);
      }
      // For external links, let the default behavior handle it
    },
    [href, navigate]
  );

  // IQCP deep links use client-side navigation
  if (href && isIqcpDeepLink(href)) {
    return (
      <a
        href={href}
        onClick={handleClick}
        className={`
          text-blue-600 hover:text-blue-800
          underline decoration-blue-400 hover:decoration-blue-600
          font-medium transition-colors
          ${className ?? ''}
        `}
        {...rest}
      >
        {children}
      </a>
    );
  }

  // External links open in new tab
  return (
    <a
      href={href}
      target="_blank"
      rel="noopener noreferrer"
      className={`
        text-blue-600 hover:text-blue-800
        underline decoration-blue-400 hover:decoration-blue-600
        transition-colors
        ${className ?? ''}
      `}
      {...rest}
    >
      {children}
    </a>
  );
}

export default DeepLinkAnchor;
