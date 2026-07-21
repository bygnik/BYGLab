# BYGLab

Engine modeling suite. Core solver: a 1-D finite-volume gas-dynamics + 0-D cylinder model, written in Rust, compiled to WebAssembly, and run **entirely client-side in the browser** — no server-side compute. BMW S54B32 is the reference engine.

## Architecture

| Layer | Stack |
|-------|--------|
| Solver | `rust/crates/byglab-core` — pure Rust, no wasm-specific deps, `cargo test`-able natively |
| Browser binding | `rust/crates/byglab-wasm` — thin `wasm-bindgen` wrapper, no physics of its own |
| Native CLI | `rust/crates/byglab-cli` — same core, JSON in/out, for fast local iteration and reference-data comparison |
| Web | React + Vite + Recharts (`services/web`), solver runs in a Web Worker via wasm |
| Local run | Docker Compose (`web` on :8080) — static SPA, single container, no backend |
| Solver validation | `benchmarks/openwam/` — legacy OpenWAM reference cases (exact Riemann solutions, physically-validated single- and 6-cylinder S54B32 engine models) used as ground truth for the Rust port |

Future modules (Flowbench, Dyno, Turbo, Camshaft, Exhaust, Intake, Combustion) appear in the nav as Coming Soon.

## Quick start (Docker)

```bash
docker compose up --build
```

- App: http://localhost:8080

## Local development (without Docker)

### Rust solver

Requires `rustup` (native install, not Docker — Rust cross-compiles to `wasm32-unknown-unknown` fine on any host) and `wasm-pack`:

```bash
rustup target add wasm32-unknown-unknown
cargo install wasm-pack
```

Fast core-logic iteration (no wasm/browser involved — this is where nearly all solver work happens):

```bash
cd rust
cargo test
```

Rebuild the wasm package after touching `byglab-wasm`'s binding surface:

```bash
wasm-pack build rust/crates/byglab-wasm --target web --dev   # fast, for iteration
wasm-pack build rust/crates/byglab-wasm --target web         # release, for production builds
```

### Web

```bash
cd services/web
npm install
npm run dev
```

The Worker (`src/wasm-worker.js`) imports the wasm-pack output via a relative path into `rust/crates/byglab-wasm/pkg/` — rebuild wasm (above) after Rust changes; Vite picks up the new `pkg/` files on refresh.

### CLI (fast native iteration / reference comparison)

```bash
cd rust
cargo run -p byglab-cli -- crates/byglab-cli/example-config.json
```

## S54B32 reference specs

Validated against OpenWAM in `benchmarks/openwam/` (see that directory's README for the full derivation and physical-plausibility checks):

| Param | Value |
|-------|-------|
| Bore / stroke / rod | 87 / 91 / 139 mm |
| Compression ratio | 11.5:1 |
| Firing order | 1-5-3-6-2-4 |
| Valves | Ø35 mm intake ×2, Ø30.5 mm exhaust ×2 |
| Throat | 85% |
| Cd | 0.75 |
| Cam | 260° @ 0.15 mm, 12 mm lift |

## Solver porting roadmap

`byglab-core` currently contains only a placeholder (`SimConfig`/`SimResult`/`run()`) proving the Rust → wasm → Worker → React pipeline end to end — see the "WASM Smoke Test" nav item. Real physics is ported next, in phases mirroring the OpenWAM validation ladder already built and checked in `benchmarks/openwam/`:

1. 1-D Euler finite-volume pipe solver (quiescent pipe → linear acoustic resonance → Sod shock tube), checked against `benchmarks/openwam/analysis/exact_riemann.py` and the captured OpenWAM `.DAT` reference output.
2. Pipe network topology (junctions, branches).
3. 0-D cylinder + Wiebe combustion + valve models, checked against the single-cylinder S54B32 OpenWAM case.
4. Multi-cylinder firing order + full 6-2-1 exhaust manifold, checked against the 6-cylinder S54B32 OpenWAM case.

`WavePage.jsx` (the pre-Rust prototype UI) is not yet wired to the new solver and will show a broken/empty state until phase 3-4 lands — it's kept as a reference for the eventual real wave-dynamics UI, not currently functional.

## Future cloud deployment

The web container is a single static SPA with no backend — deployable to any static host (Azure Static Web Apps, Netlify, GitHub Pages, S3+CloudFront, etc.) in addition to the Docker/nginx path above. No server-side compute, no COOP/COEP headers required (the wasm build is single-threaded for now — see `rust/crates/byglab-wasm`).

## Legacy files

Root `IntakeMachCalculator.jsx` and the incomplete `.crdownload` Python twin are superseded by the Rust solver + `services/web`. The prior Python prototype (`packages/byglab_engine`, `services/api`) has been retired — OpenWAM (`benchmarks/openwam/`) is the reference implementation the Rust port is validated against instead.
