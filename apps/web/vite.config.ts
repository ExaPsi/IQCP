import { defineConfig, type Plugin } from 'vite';
import react from '@vitejs/plugin-react-swc';
import wasm from 'vite-plugin-wasm';
import topLevelAwait from 'vite-plugin-top-level-await';
import path from 'path';
import { readFileSync } from 'fs';

// Custom plugin for COOP/COEP headers - enables SharedArrayBuffer for WASM threading
// NOTE: These headers break browser extensions (MetaMask, etc.)
// Test in Incognito mode or browser profile without extensions
const coopCoepPlugin = (): Plugin => ({
  name: 'coop-coep-headers',
  configureServer(server) {
    console.log('[COOP/COEP] Installing cross-origin isolation headers middleware');
    server.middlewares.use((_req, res, next) => {
      res.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
      res.setHeader('Cross-Origin-Embedder-Policy', 'credentialless');
      next();
    });
  },
  configurePreviewServer(server) {
    server.middlewares.use((_req, res, next) => {
      res.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
      res.setHeader('Cross-Origin-Embedder-Policy', 'credentialless');
      next();
    });
  },
});

// Read version from package.json for build-time injection
const pkg = JSON.parse(readFileSync(path.resolve(__dirname, 'package.json'), 'utf-8'));

// https://vite.dev/config/
export default defineConfig({
  plugins: [
    coopCoepPlugin(),
    react(),
    wasm(),
    topLevelAwait(),
  ],
  define: {
    __APP_VERSION__: JSON.stringify(pkg.version),
  },
  assetsInclude: ['**/*.md'],
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
  build: {
    target: 'esnext',
    rollupOptions: {
      // `initThreadPool` is only present in the parallel-feature WASM build;
      // the worker probes for it at runtime via dynamic property access. The
      // single-threaded build legitimately omits the export, so silence the
      // expected static-analysis warning to keep the build log clean.
      onwarn(warning, warn) {
        if (
          warning.code === 'MISSING_EXPORT' &&
          typeof warning.message === 'string' &&
          warning.message.includes('"initThreadPool" is not exported')
        ) {
          return;
        }
        warn(warning);
      },
      output: {
        manualChunks: {
          // Phase 2 (US-035): Code-split Three.js and R3F into a separate chunk
          // that is only loaded when the 3D molecular viewer is activated.
          // This keeps the initial Module C bundle small.
          'viewer3d': [
            'three',
            '@react-three/fiber',
            '@react-three/drei',
          ],
        },
      },
    },
  },
  worker: {
    format: 'es',
    plugins: () => [wasm(), topLevelAwait()],
  },
  optimizeDeps: {
    exclude: ['qc-wasm', 'qc-wasm-spectra'],
  },
});
