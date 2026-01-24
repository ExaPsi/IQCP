// Generate all deep links for IQCP lab pack with verified encoding
import LZString from 'lz-string';

const BASE = process.argv[2] === '--prod' ? 'https://iqcp.dev' : 'http://localhost:5173';

function encode(state) {
  const json = JSON.stringify(state);
  const encoded = LZString.compressToEncodedURIComponent(json);
  // Verify round-trip
  const decoded = LZString.decompressFromEncodedURIComponent(encoded);
  if (decoded !== json) {
    throw new Error('Round-trip verification failed!');
  }
  return encoded;
}

function boysLink(m, T, view = 'single', mode = 'explain') {
  const state = {
    schema_version: 'run_state_v1',
    app_version: '0.1.0',
    module: 'boys',
    boys: { m, T, view },
    ui: { mode }
  };
  return `${BASE}/boys?run=${encode(state)}`;
}

function rysLink(n, T, target = '1e-6', mode = 'explain') {
  const state = {
    schema_version: 'run_state_v1',
    app_version: '0.1.0',
    module: 'rys',
    rys: { n, T, target },
    ui: { mode }
  };
  return `${BASE}/rys?run=${encode(state)}`;
}

function scfLink(system_id, diis, conv = 'medium', max_iter = 50, mode = 'explain') {
  const state = {
    schema_version: 'run_state_v1',
    app_version: '0.1.0',
    module: 'scf',
    scf: { system_id, max_iter, conv, diis },
    ui: { mode }
  };
  return `${BASE}/scf?run=${encode(state)}`;
}

console.log(`Base URL: ${BASE}\n`);

// Output as JSON for easy parsing
const links = {
  boys: {
    // Worksheet links
    "default_m0_T1": boysLink(0, 1.0),
    "very_small_T_m0_T0.01": boysLink(0, 0.01),
    "small_T_m0_T0.5": boysLink(0, 0.5),
    "medium_T_m0_T10": boysLink(0, 10.0),
    "large_T_m0_T35": boysLink(0, 35.0),
    "very_large_T_m0_T100": boysLink(0, 100.0),
    "m5_T10": boysLink(5, 10.0),
    "m0_T15": boysLink(0, 15.0),
    "m5_T45_Q2.6": boysLink(5, 45.0),
    "sweep_m0": boysLink(0, 10.0, 'sweep'),
    // Internals mode
    "internals_m0_T0.5": boysLink(0, 0.5, 'single', 'internals'),
    "internals_m0_T15": boysLink(0, 15.0, 'single', 'internals'),
    "internals_m0_T35": boysLink(0, 35.0, 'single', 'internals'),
    // PT-Boys scenarios
    "pt_scenario_A_m3_T9.5": boysLink(3, 9.5, 'single', 'internals'),
    "pt_scenario_B_m0_T28": boysLink(0, 28.0, 'single', 'internals'),
    "pt_scenario_C_m8_T14": boysLink(8, 14.0, 'single', 'internals'),
  },
  rys: {
    // Worksheet links
    "default_n3_T10": rysLink(3, 10.0),
    "n5_T10": rysLink(5, 10.0),
    "n7_T10": rysLink(7, 10.0),
    "n3_T25": rysLink(3, 25.0),
    "n5_T15_1e-8": rysLink(5, 15.0, '1e-8'),
    "pp_pp_n3_T10": rysLink(3, 10.0),  // (pp|pp) L=4, n=3
    "dd_pp_n4_T10": rysLink(4, 10.0),  // (dd|pp) L=6, n=4
    // PT-Rys scenarios
    "pt_scenario_A_T8_1e-6": rysLink(3, 8.0, '1e-6'),
    "pt_scenario_B_T18_1e-8": rysLink(5, 18.0, '1e-8'),
    "pt_scenario_C_T5_1e-8": rysLink(5, 5.0, '1e-8'),
  },
  scf: {
    // Worksheet links
    "h2_diis": scfLink('h2_sto3g_r1.4', true),
    "h2_no_diis": scfLink('h2_sto3g_r1.4', false),
    "h2o_diis": scfLink('h2o_sto3g', true),
    "h2o_no_diis": scfLink('h2o_sto3g', false),
    "h2_tight_diis": scfLink('h2_sto3g_r1.4', true, 'tight'),
    "h2_tight_no_diis": scfLink('h2_sto3g_r1.4', false, 'tight'),
    "h2_internals": scfLink('h2_sto3g_r1.4', true, 'medium', 50, 'internals'),
    // PT-SCF scenarios
    "lih_no_diis": scfLink('lih_sto3g', false),
    "lih_diis": scfLink('lih_sto3g', true),
    "nh3_no_diis": scfLink('nh3_sto3g', false),
    "nh3_diis": scfLink('nh3_sto3g', true),
  }
};

console.log(JSON.stringify(links, null, 2));
