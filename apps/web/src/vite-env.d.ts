/// <reference types="vite/client" />

// Type declaration for plotly.js-basic-dist (uses same types as plotly.js)
declare module 'plotly.js-basic-dist' {
  import Plotly from 'plotly.js';
  export = Plotly;
}

// Type declaration for plotly.js-cartesian-dist (includes heatmap, contour, etc.)
declare module 'plotly.js-cartesian-dist' {
  import Plotly from 'plotly.js';
  export = Plotly;
}

// Type declarations for preset system JSON files
// These are bundled by Vite at build time
declare module '../../../../content/presets/systems/h2_sto3g_r1.4.json' {
  const value: {
    system_id: string;
    label: string;
    description?: string;
    geometry?: { atoms: Array<{ symbol: string; xyz: number[] }>; units: string };
    basis_id: string;
    nbf: number;
    nelec: number;
    e_nuc: number;
    s_matrix: number[];
    h_core: number[];
    eri_compressed: number[];
  };
  export default value;
}

declare module '../../../../content/presets/systems/heh_plus_sto3g.json' {
  const value: {
    system_id: string;
    label: string;
    description?: string;
    geometry?: { atoms: Array<{ symbol: string; xyz: number[] }>; units: string };
    basis_id: string;
    nbf: number;
    nelec: number;
    e_nuc: number;
    s_matrix: number[];
    h_core: number[];
    eri_compressed: number[];
  };
  export default value;
}

declare module '../../../../content/presets/systems/lih_sto3g.json' {
  const value: {
    system_id: string;
    label: string;
    description?: string;
    geometry?: { atoms: Array<{ symbol: string; xyz: number[] }>; units: string };
    basis_id: string;
    nbf: number;
    nelec: number;
    e_nuc: number;
    s_matrix: number[];
    h_core: number[];
    eri_compressed: number[];
  };
  export default value;
}

declare module '../../../../content/presets/systems/h2o_sto3g.json' {
  const value: {
    system_id: string;
    label: string;
    description?: string;
    geometry?: { atoms: Array<{ symbol: string; xyz: number[] }>; units: string };
    basis_id: string;
    nbf: number;
    nelec: number;
    e_nuc: number;
    s_matrix: number[];
    h_core: number[];
    eri_compressed: number[];
  };
  export default value;
}

declare module '../../../../content/presets/systems/nh3_sto3g.json' {
  const value: {
    system_id: string;
    label: string;
    description?: string;
    geometry?: { atoms: Array<{ symbol: string; xyz: number[] }>; units: string };
    basis_id: string;
    nbf: number;
    nelec: number;
    e_nuc: number;
    s_matrix: number[];
    h_core: number[];
    eri_compressed: number[];
  };
  export default value;
}
