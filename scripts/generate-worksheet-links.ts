#!/usr/bin/env npx ts-node

/**
 * Lab Pack #1 Deep Link Generator
 *
 * This script generates and validates all deep link URLs for the student worksheet.
 * It uses the same encoding as the IQCP web application to ensure compatibility.
 *
 * Usage:
 *   npx ts-node scripts/generate-worksheet-links.ts
 *   npx ts-node scripts/generate-worksheet-links.ts --verify
 *   npx ts-node scripts/generate-worksheet-links.ts --output markdown
 *
 * @module generate-worksheet-links
 */

import LZString from 'lz-string';

// ============================================================================
// Type Definitions (mirroring apps/web/src/types/run-state.ts)
// ============================================================================

type ModuleId = 'boys' | 'rys' | 'scf';

interface BoysParams {
  m: number;
  T: number;
  view: 'single' | 'sweep';
}

interface RysParams {
  n: number;
  T: number;
  target: '1e-4' | '1e-6' | '1e-8';
}

interface ScfParams {
  system_id: string;
  max_iter: number;
  conv: 'loose' | 'medium' | 'tight';
  diis: boolean;
}

interface UiSettings {
  mode: 'explain' | 'internals';
  theme?: 'light' | 'dark';
}

interface RunStateV1 {
  schema_version: 'run_state_v1';
  app_version: string;
  module: ModuleId;
  boys?: BoysParams;
  rys?: RysParams;
  scf?: ScfParams;
  ui: UiSettings;
}

interface DeepLinkEntry {
  id: string;
  label: string;
  description: string;
  section: string;
  state: RunStateV1;
  url?: string;
}

// ============================================================================
// Configuration
// ============================================================================

const APP_VERSION = '0.1.0';
const BASE_URL_PROD = 'https://iqcp.dev';
const BASE_URL_LOCAL = 'http://localhost:5173';

// ============================================================================
// Link Encoding (mirrors apps/web/src/lib/deeplink.ts)
// ============================================================================

function encodeRunState(state: RunStateV1): string {
  const json = JSON.stringify(state);
  return LZString.compressToEncodedURIComponent(json);
}

function decodeRunState(encoded: string): RunStateV1 | null {
  try {
    const json = LZString.decompressFromEncodedURIComponent(encoded);
    if (!json) return null;
    return JSON.parse(json) as RunStateV1;
  } catch {
    return null;
  }
}

function verifyRoundTrip(state: RunStateV1): boolean {
  const encoded = encodeRunState(state);
  const decoded = decodeRunState(encoded);
  if (!decoded) return false;
  return JSON.stringify(state) === JSON.stringify(decoded);
}

// ============================================================================
// Link Definitions
// ============================================================================

const worksheetLinks: DeepLinkEntry[] = [
  // Section 2: Boys Function Exploration
  {
    id: 'boys-2.1',
    label: 'Default Boys view',
    description: 'F_0(T=1.0) with explain mode - starting point for exploration',
    section: 'boys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'boys',
      boys: { m: 0, T: 1.0, view: 'single' },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'boys-2.2a',
    label: 'Small T (series regime)',
    description: 'F_0(T=0.5) showing series regime behavior with internals',
    section: 'boys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'boys',
      boys: { m: 0, T: 0.5, view: 'single' },
      ui: { mode: 'internals' },
    },
  },
  {
    id: 'boys-2.2b',
    label: 'Very small T',
    description: 'F_0(T=0.01) approaching the limiting value 1/(2m+1)',
    section: 'boys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'boys',
      boys: { m: 0, T: 0.01, view: 'single' },
      ui: { mode: 'internals' },
    },
  },
  {
    id: 'boys-2.3a',
    label: 'Moderate T (recurrence regime)',
    description: 'F_0(T=15.0) showing recurrence regime behavior',
    section: 'boys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'boys',
      boys: { m: 0, T: 15.0, view: 'single' },
      ui: { mode: 'internals' },
    },
  },
  {
    id: 'boys-2.3b',
    label: 'Large T (asymptotic regime)',
    description: 'F_0(T=35.0) showing asymptotic behavior',
    section: 'boys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'boys',
      boys: { m: 0, T: 35.0, view: 'single' },
      ui: { mode: 'internals' },
    },
  },
  {
    id: 'boys-2.4a',
    label: 'Higher order m',
    description: 'F_5(T=10.0) comparing higher order behavior',
    section: 'boys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'boys',
      boys: { m: 5, T: 10.0, view: 'single' },
      ui: { mode: 'internals' },
    },
  },
  {
    id: 'boys-2.4b',
    label: 'Regime boundary',
    description: 'F_0(T=12.0) near series/recurrence transition',
    section: 'boys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'boys',
      boys: { m: 0, T: 12.0, view: 'single' },
      ui: { mode: 'internals' },
    },
  },
  {
    id: 'boys-checkpoint',
    label: 'Checkpoint artifact state',
    description: 'F_3(T=8.0) for artifact export checkpoint',
    section: 'boys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'boys',
      boys: { m: 3, T: 8.0, view: 'single' },
      ui: { mode: 'internals' },
    },
  },

  // Section 3: Rys Quadrature Exploration
  {
    id: 'rys-3.1',
    label: 'Default Rys view',
    description: 'n=3 quadrature at T=5.0 - starting point',
    section: 'rys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'rys',
      rys: { n: 3, T: 5.0, target: '1e-6' },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'rys-3.2a',
    label: 'Inspect roots and weights',
    description: 'n=5 at T=10.0 - examine quadrature points',
    section: 'rys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'rys',
      rys: { n: 5, T: 10.0, target: '1e-6' },
      ui: { mode: 'internals' },
    },
  },
  {
    id: 'rys-3.2b',
    label: 'High quadrature order',
    description: 'n=10 at T=10.0 - compare with lower orders',
    section: 'rys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'rys',
      rys: { n: 10, T: 10.0, target: '1e-6' },
      ui: { mode: 'internals' },
    },
  },
  {
    id: 'rys-3.3a',
    label: 'Error curve at moderate T',
    description: 'n=5 at T=10.0 - examine error vs order relationship',
    section: 'rys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'rys',
      rys: { n: 5, T: 10.0, target: '1e-6' },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'rys-3.3b',
    label: 'Error curve at larger T',
    description: 'n=5 at T=25.0 - observe T dependence',
    section: 'rys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'rys',
      rys: { n: 5, T: 25.0, target: '1e-6' },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'rys-3.4a',
    label: 'Target 1e-6 accuracy',
    description: 'Finding order needed for moderate accuracy',
    section: 'rys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'rys',
      rys: { n: 5, T: 10.0, target: '1e-6' },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'rys-3.4b',
    label: 'Target 1e-8 accuracy',
    description: 'Finding order needed for high accuracy',
    section: 'rys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'rys',
      rys: { n: 5, T: 10.0, target: '1e-8' },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'rys-checkpoint',
    label: 'Checkpoint artifact state',
    description: 'T=15.0 with 1e-8 target for artifact export',
    section: 'rys',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'rys',
      rys: { n: 5, T: 15.0, target: '1e-8' },
      ui: { mode: 'internals' },
    },
  },

  // Section 4: SCF Convergence Exploration
  {
    id: 'scf-4.1',
    label: 'H2 default run',
    description: 'H2 molecule with medium convergence and DIIS',
    section: 'scf',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'scf',
      scf: { system_id: 'h2_sto3g_r1.4', max_iter: 50, conv: 'medium', diis: true },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'scf-4.2a',
    label: 'H2 tight convergence (no DIIS)',
    description: 'H2 with tight tolerance, DIIS disabled - observe slow convergence',
    section: 'scf',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'scf',
      scf: { system_id: 'h2_sto3g_r1.4', max_iter: 50, conv: 'tight', diis: false },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'scf-4.2b',
    label: 'H2 tight convergence (with DIIS)',
    description: 'H2 with tight tolerance, DIIS enabled - observe acceleration',
    section: 'scf',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'scf',
      scf: { system_id: 'h2_sto3g_r1.4', max_iter: 50, conv: 'tight', diis: true },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'scf-4.3a',
    label: 'H2O without DIIS',
    description: 'Water molecule with medium convergence, DIIS disabled',
    section: 'scf',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'scf',
      scf: { system_id: 'h2o_sto3g', max_iter: 50, conv: 'medium', diis: false },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'scf-4.3b',
    label: 'H2O with DIIS',
    description: 'Water molecule with medium convergence, DIIS enabled',
    section: 'scf',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'scf',
      scf: { system_id: 'h2o_sto3g', max_iter: 50, conv: 'medium', diis: true },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'scf-4.4',
    label: 'H2 matrix inspection',
    description: 'H2 with internals mode - view Fock, density matrices',
    section: 'scf',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'scf',
      scf: { system_id: 'h2_sto3g_r1.4', max_iter: 50, conv: 'medium', diis: true },
      ui: { mode: 'internals' },
    },
  },
  {
    id: 'scf-4.5',
    label: 'LiH system',
    description: 'Lithium hydride - more complex system with DIIS',
    section: 'scf',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'scf',
      scf: { system_id: 'lih_sto3g', max_iter: 50, conv: 'tight', diis: true },
      ui: { mode: 'explain' },
    },
  },
  {
    id: 'scf-checkpoint',
    label: 'Checkpoint artifact state',
    description: 'H2O with DIIS for artifact export',
    section: 'scf',
    state: {
      schema_version: 'run_state_v1',
      app_version: APP_VERSION,
      module: 'scf',
      scf: { system_id: 'h2o_sto3g', max_iter: 50, conv: 'medium', diis: true },
      ui: { mode: 'internals' },
    },
  },
];

// ============================================================================
// URL Generation
// ============================================================================

function generateUrl(entry: DeepLinkEntry, baseUrl: string = BASE_URL_PROD): string {
  const encoded = encodeRunState(entry.state);
  return `${baseUrl}/${entry.state.module}?s=${encoded}`;
}

function generateAllUrls(
  links: DeepLinkEntry[],
  baseUrl: string = BASE_URL_PROD
): DeepLinkEntry[] {
  return links.map((link) => ({
    ...link,
    url: generateUrl(link, baseUrl),
  }));
}

// ============================================================================
// Verification
// ============================================================================

interface VerificationResult {
  id: string;
  label: string;
  roundTripOk: boolean;
  urlLength: number;
  urlSafe: boolean;
  issues: string[];
}

function verifyLink(entry: DeepLinkEntry): VerificationResult {
  const issues: string[] = [];
  const roundTripOk = verifyRoundTrip(entry.state);
  const encoded = encodeRunState(entry.state);
  const urlLength = encoded.length;
  const urlSafe = urlLength < 500;

  if (!roundTripOk) {
    issues.push('Round-trip encoding/decoding failed');
  }
  if (!urlSafe) {
    issues.push(`URL encoding too long (${urlLength} chars, max 500)`);
  }

  // Validate state structure
  const state = entry.state;
  if (state.module === 'boys' && !state.boys) {
    issues.push('Boys module missing boys params');
  }
  if (state.module === 'rys' && !state.rys) {
    issues.push('Rys module missing rys params');
  }
  if (state.module === 'scf' && !state.scf) {
    issues.push('SCF module missing scf params');
  }

  return {
    id: entry.id,
    label: entry.label,
    roundTripOk,
    urlLength,
    urlSafe,
    issues,
  };
}

function verifyAllLinks(links: DeepLinkEntry[]): VerificationResult[] {
  return links.map(verifyLink);
}

// ============================================================================
// Output Formatters
// ============================================================================

function outputJson(links: DeepLinkEntry[]): void {
  const output = {
    version: '1.0',
    generated: new Date().toISOString(),
    base_url: BASE_URL_PROD,
    local_url: BASE_URL_LOCAL,
    links: generateAllUrls(links),
  };
  console.log(JSON.stringify(output, null, 2));
}

function outputMarkdown(links: DeepLinkEntry[]): void {
  const linksWithUrls = generateAllUrls(links);

  console.log('# Lab Pack #1 Deep Links\n');
  console.log(`Generated: ${new Date().toISOString()}\n`);

  // Group by section
  const sections = new Map<string, DeepLinkEntry[]>();
  for (const link of linksWithUrls) {
    const existing = sections.get(link.section) || [];
    existing.push(link);
    sections.set(link.section, existing);
  }

  for (const [section, sectionLinks] of sections) {
    console.log(`## ${section.charAt(0).toUpperCase() + section.slice(1)} Module\n`);
    console.log('| ID | Label | Description |');
    console.log('|----|-------|-------------|');
    for (const link of sectionLinks) {
      console.log(`| [${link.id}](${link.url}) | ${link.label} | ${link.description} |`);
    }
    console.log('');
  }
}

function outputVerification(links: DeepLinkEntry[]): void {
  const results = verifyAllLinks(links);
  let allPassed = true;

  console.log('# Deep Link Verification Report\n');
  console.log(`Verified: ${new Date().toISOString()}`);
  console.log(`Total links: ${results.length}\n`);

  console.log('| ID | Round-trip | URL Length | URL Safe | Issues |');
  console.log('|----|------------|------------|----------|--------|');

  for (const result of results) {
    const rtIcon = result.roundTripOk ? 'OK' : 'FAIL';
    const safeIcon = result.urlSafe ? 'OK' : 'WARN';
    const issueText = result.issues.length > 0 ? result.issues.join('; ') : '-';

    if (result.issues.length > 0) {
      allPassed = false;
    }

    console.log(
      `| ${result.id} | ${rtIcon} | ${result.urlLength} | ${safeIcon} | ${issueText} |`
    );
  }

  console.log('');
  console.log(`## Summary\n`);
  console.log(`- Total links: ${results.length}`);
  console.log(`- All round-trips OK: ${results.every((r) => r.roundTripOk)}`);
  console.log(`- All URLs safe: ${results.every((r) => r.urlSafe)}`);
  console.log(`- Overall: ${allPassed ? 'PASS' : 'FAIL'}`);
}

// ============================================================================
// CLI Interface
// ============================================================================

function main(): void {
  const args = process.argv.slice(2);

  if (args.includes('--help') || args.includes('-h')) {
    console.log(`
Lab Pack #1 Deep Link Generator

Usage:
  npx ts-node scripts/generate-worksheet-links.ts [options]

Options:
  --output json       Output all links as JSON (default)
  --output markdown   Output all links as Markdown table
  --verify            Verify all links and output report
  --local             Use localhost URLs instead of production
  --help, -h          Show this help message

Examples:
  npx ts-node scripts/generate-worksheet-links.ts
  npx ts-node scripts/generate-worksheet-links.ts --verify
  npx ts-node scripts/generate-worksheet-links.ts --output markdown
  npx ts-node scripts/generate-worksheet-links.ts --local --output json
`);
    process.exit(0);
  }

  if (args.includes('--verify')) {
    outputVerification(worksheetLinks);
  } else if (args.includes('--output') && args.includes('markdown')) {
    outputMarkdown(worksheetLinks);
  } else {
    outputJson(worksheetLinks);
  }
}

// Run if executed directly
main();
