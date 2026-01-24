import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react-swc';
import wasm from 'vite-plugin-wasm';
import topLevelAwait from 'vite-plugin-top-level-await';
import path from 'path';
// Custom plugin for COOP/COEP headers - enables SharedArrayBuffer for WASM threading
// NOTE: These headers break browser extensions (MetaMask, etc.)
// Test in Incognito mode or browser profile without extensions
const coopCoepPlugin = () => ({
    name: 'coop-coep-headers',
    configureServer(server) {
        console.log('[COOP/COEP] Installing cross-origin isolation headers middleware');
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
});
// https://vite.dev/config/
export default defineConfig({
    plugins: [
        coopCoepPlugin(),
        react(),
        wasm(),
        topLevelAwait(),
    ],
    assetsInclude: ['**/*.md'],
    resolve: {
        alias: {
            '@': path.resolve(__dirname, './src'),
        },
    },
    build: {
        target: 'esnext',
    },
    worker: {
        format: 'es',
        plugins: () => [wasm(), topLevelAwait()],
    },
    optimizeDeps: {
        exclude: ['qc-wasm'],
    },
});
