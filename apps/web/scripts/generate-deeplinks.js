// Generate deep links for IQCP lab pack
import LZString from 'lz-string';

function encodeState(state) {
  const json = JSON.stringify(state);
  return LZString.compressToEncodedURIComponent(json);
}

function makeBoysLink(m, T, view = 'single', mode = 'explain') {
  const state = {
    schema_version: 'run_state_v1',
    app_version: '0.1.0',
    module: 'boys',
    boys: { m, T, view },
    ui: { mode }
  };
  return `https://iqcp.dev/boys?s=${encodeState(state)}`;
}

function makeRysLink(n, T, target = '1e-6', mode = 'explain') {
  const state = {
    schema_version: 'run_state_v1',
    app_version: '0.1.0',
    module: 'rys',
    rys: { n, T, target },
    ui: { mode }
  };
  return `https://iqcp.dev/rys?s=${encodeState(state)}`;
}

function makeScfLink(system_id, diis, conv = 'medium', max_iter = 50, mode = 'explain') {
  const state = {
    schema_version: 'run_state_v1',
    app_version: '0.1.0',
    module: 'scf',
    scf: { system_id, max_iter, conv, diis },
    ui: { mode }
  };
  return `https://iqcp.dev/scf?s=${encodeState(state)}`;
}

// Use localhost for testing
const BASE = process.argv[2] === '--prod' ? 'https://iqcp.dev' : 'http://localhost:5173';

function makeBoysLinkLocal(m, T, view = 'single', mode = 'explain') {
  const state = {
    schema_version: 'run_state_v1',
    app_version: '0.1.0',
    module: 'boys',
    boys: { m, T, view },
    ui: { mode }
  };
  return `${BASE}/boys?s=${encodeState(state)}`;
}

function makeRysLinkLocal(n, T, target = '1e-6', mode = 'explain') {
  const state = {
    schema_version: 'run_state_v1',
    app_version: '0.1.0',
    module: 'rys',
    rys: { n, T, target },
    ui: { mode }
  };
  return `${BASE}/rys?s=${encodeState(state)}`;
}

function makeScfLinkLocal(system_id, diis, conv = 'medium', max_iter = 50, mode = 'explain') {
  const state = {
    schema_version: 'run_state_v1',
    app_version: '0.1.0',
    module: 'scf',
    scf: { system_id, max_iter, conv, diis },
    ui: { mode }
  };
  return `${BASE}/scf?s=${encodeState(state)}`;
}

console.log(`Using base URL: ${BASE}\n`);
console.log('=== BOYS MODULE DEEP LINKS ===\n');

// Boys links for worksheet
console.log('Very Small T (m=0, T=0.01):');
console.log(makeBoysLinkLocal(0, 0.01, 'single'));
console.log('\nSmall T (m=0, T=0.5):');
console.log(makeBoysLinkLocal(0, 0.5, 'single'));
console.log('\nMedium T (m=0, T=10.0):');
console.log(makeBoysLinkLocal(0, 10.0, 'single'));
console.log('\nLarge T (m=0, T=35.0):');
console.log(makeBoysLinkLocal(0, 35.0, 'single'));
console.log('\nVery Large T (m=0, T=100.0):');
console.log(makeBoysLinkLocal(0, 100.0, 'single'));
console.log('\nm=5, T=10.0:');
console.log(makeBoysLinkLocal(5, 10.0, 'single'));
console.log('\nm=0, T=15.0:');
console.log(makeBoysLinkLocal(0, 15.0, 'single'));
console.log('\nm=5, T=45.0 (Q2.6):');
console.log(makeBoysLinkLocal(5, 45.0, 'single'));
console.log('\nm=0, T=35.0 (Q2.6):');
console.log(makeBoysLinkLocal(0, 35.0, 'single'));
console.log('\nSweep view (m=0):');
console.log(makeBoysLinkLocal(0, 10.0, 'sweep'));
console.log('\nInternals mode (m=0, T=0.5):');
console.log(makeBoysLinkLocal(0, 0.5, 'single', 'internals'));
console.log('\nInternals mode (m=0, T=15.0):');
console.log(makeBoysLinkLocal(0, 15.0, 'single', 'internals'));
console.log('\nInternals mode (m=0, T=35.0):');
console.log(makeBoysLinkLocal(0, 35.0, 'single', 'internals'));

// Performance task Boys
console.log('\n--- PT-Boys Scenarios ---');
console.log('Scenario A (T=9.5, m=3):');
console.log(makeBoysLinkLocal(3, 9.5, 'single', 'internals'));
console.log('\nScenario B (T=28.0, m=0):');
console.log(makeBoysLinkLocal(0, 28.0, 'single', 'internals'));
console.log('\nScenario C (T=14.0, m=8):');
console.log(makeBoysLinkLocal(8, 14.0, 'single', 'internals'));

console.log('\n=== RYS MODULE DEEP LINKS ===\n');

// Rys links for worksheet
console.log('n=3, T=10.0:');
console.log(makeRysLinkLocal(3, 10.0));
console.log('\nn=5, T=10.0:');
console.log(makeRysLinkLocal(5, 10.0));
console.log('\nn=7, T=10.0:');
console.log(makeRysLinkLocal(7, 10.0));
console.log('\nn=3, T=25.0:');
console.log(makeRysLinkLocal(3, 25.0));
console.log('\nn=5, T=15.0, target 1e-8:');
console.log(makeRysLinkLocal(5, 15.0, '1e-8'));
console.log('\nn=3, T=10.0 (shell quartet pp|pp):');
console.log(makeRysLinkLocal(3, 10.0));
console.log('\nn=4, T=10.0 (shell quartet dd|pp):');
console.log(makeRysLinkLocal(4, 10.0));

// Performance task Rys
console.log('\n--- PT-Rys Scenarios ---');
console.log('Scenario A (T=8.0, target=1e-6):');
console.log(makeRysLinkLocal(3, 8.0, '1e-6'));
console.log('\nScenario B (T=18.0, target=1e-8):');
console.log(makeRysLinkLocal(5, 18.0, '1e-8'));
console.log('\nScenario C (T=5.0, target=1e-8):');
console.log(makeRysLinkLocal(5, 5.0, '1e-8'));

console.log('\n=== SCF MODULE DEEP LINKS ===\n');

// SCF links for worksheet
console.log('H2 with DIIS:');
console.log(makeScfLinkLocal('h2_sto3g_r1.4', true));
console.log('\nH2 without DIIS:');
console.log(makeScfLinkLocal('h2_sto3g_r1.4', false));
console.log('\nH2O with DIIS:');
console.log(makeScfLinkLocal('h2o_sto3g', true));
console.log('\nH2O without DIIS:');
console.log(makeScfLinkLocal('h2o_sto3g', false));
console.log('\nH2 tight convergence with DIIS:');
console.log(makeScfLinkLocal('h2_sto3g_r1.4', true, 'tight'));
console.log('\nH2 tight convergence without DIIS:');
console.log(makeScfLinkLocal('h2_sto3g_r1.4', false, 'tight'));
console.log('\nH2 internals mode:');
console.log(makeScfLinkLocal('h2_sto3g_r1.4', true, 'medium', 50, 'internals'));

// Performance task SCF
console.log('\n--- PT-SCF Scenarios ---');
console.log('LiH without DIIS:');
console.log(makeScfLinkLocal('lih_sto3g', false));
console.log('\nLiH with DIIS:');
console.log(makeScfLinkLocal('lih_sto3g', true));
console.log('\nNH3 without DIIS:');
console.log(makeScfLinkLocal('nh3_sto3g', false));
console.log('\nNH3 with DIIS:');
console.log(makeScfLinkLocal('nh3_sto3g', true));
