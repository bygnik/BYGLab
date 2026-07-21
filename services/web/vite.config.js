import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import path from "node:path";
import { fileURLToPath } from "node:url";

const repoRoot = path.resolve(fileURLToPath(new URL(".", import.meta.url)), "../..");

export default defineConfig({
  plugins: [react()],
  server: {
    host: true,
    port: 5173,
    fs: {
      // wasm-worker.js imports the wasm-pack output from rust/crates/byglab-wasm/pkg/,
      // which lives outside this project's root (services/web) - allow the whole repo
      // so Vite serves it instead of silently falling back to the SPA's index.html.
      allow: [repoRoot],
    },
  },
});
