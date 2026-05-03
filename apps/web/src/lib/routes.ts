/**
 * Centralized Route Constants
 *
 * All application routes are defined here with path-based versioning.
 * This ensures citation-stable URLs: manuscript and Lab Pack #1 deep links
 * reference `/v1/boys`, `/v1/scf`, etc. and must work indefinitely.
 *
 * To add a new version in the future:
 * 1. Add `/v2/` routes alongside `/v1/` routes
 * 2. Update CURRENT_VERSION
 * 3. Keep `/v1/` routes pointing to the v1 implementation
 *
 * @module routes
 */

/**
 * Current API/URL version prefix.
 */
export const CURRENT_VERSION = 'v1';

/**
 * Version prefix for URL paths.
 */
export const VERSION_PREFIX = `/${CURRENT_VERSION}`;

/**
 * All application routes with versioned paths.
 *
 * Module routes use the `/v1/` prefix for citation stability.
 * Content pages (labs, reference, about) also use `/v1/` for consistency.
 */
export const ROUTES = {
  /** Home page (unversioned) */
  HOME: '/',

  /** Module C: Boys Function Lab */
  BOYS: `${VERSION_PREFIX}/boys`,

  /** Module D: Rys Quadrature Lab */
  RYS: `${VERSION_PREFIX}/rys`,

  /** Module E: SCF Sandbox */
  SCF: `${VERSION_PREFIX}/scf`,

  /** Module A: Basis Explorer */
  BASIS: `${VERSION_PREFIX}/basis`,

  /** Module B: Integral Inspector */
  INTEGRALS: `${VERSION_PREFIX}/integrals`,

  /** Lab Pack materials */
  LABS: `${VERSION_PREFIX}/labs`,

  /** Reference documentation */
  REFERENCE: `${VERSION_PREFIX}/reference`,

  /** About page with citation info */
  ABOUT: `${VERSION_PREFIX}/about`,
} as const;

/**
 * Map from module ID to versioned route path.
 *
 * Used by `generateShareURL()` and deep link construction.
 */
export const MODULE_ROUTES: Record<string, string> = {
  boys: ROUTES.BOYS,
  rys: ROUTES.RYS,
  scf: ROUTES.SCF,
  basis: ROUTES.BASIS,
  integrals: ROUTES.INTEGRALS,
};

/**
 * Get the versioned route path for a module ID.
 *
 * @param moduleId - Module identifier ('boys', 'rys', 'scf')
 * @returns Versioned path (e.g., '/v1/boys')
 *
 * @example
 * ```typescript
 * getModuleRoute('boys') // '/v1/boys'
 * getModuleRoute('scf')  // '/v1/scf'
 * ```
 */
export function getModuleRoute(moduleId: string): string {
  return MODULE_ROUTES[moduleId] ?? `${VERSION_PREFIX}/${moduleId}`;
}
