import type { Plugin } from 'vite';

/**
 * Vite plugin to add COOP/COEP headers for SharedArrayBuffer support.
 *
 * These headers are required for WASM threading via wasm-bindgen-rayon.
 *
 * NOTE: These headers break browser extensions (MetaMask, etc.)
 * Test in Incognito mode or browser profile without extensions.
 */
export function coopCoepPlugin(): Plugin {
  return {
    name: 'coop-coep-headers',
    configureServer(server) {
      server.middlewares.use((_req, res, next) => {
        res.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
        res.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
        next();
      });
    },
    configurePreviewServer(server) {
      server.middlewares.use((_req, res, next) => {
        res.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
        res.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
        next();
      });
    },
  };
}

export default coopCoepPlugin;
