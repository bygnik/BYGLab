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

**Phase 1 — done.** `byglab-core` is a real 1D finite-volume gas-dynamics solver: second-order MUSCL-Hancock scheme (piecewise-linear reconstruction with a minmod slope limiter, HLLC approximate Riemann solver, explicit CFL-limited time-stepping — see `reconstruction.rs`), supporting networks of pipes joined by same-area junctions with fully symmetric second-order accuracy at the junction faces too. A clean-room Rust reimplementation informed by OpenWAM's validated physics, not a line-by-line C++ translation. Started as a simpler first-order scheme, then upgraded once validated — both stages checked against the same three analytically-checked OpenWAM cases:

- `tests/quiescent.rs` — uniform at-rest closed pipe stays at rest (no discretization error possible to hide behind).
- `tests/acoustic_resonance.rs` — closed-closed resonance period matches the exact `T = 2L/c0` to **0.09%** (OpenWAM's own higher-order scheme: 0.07% — comparable; a first-order-only version of this solver measured 0.002% here, tighter still, since numerical dissipation for a small linear disturbance mostly damps amplitude rather than shifting phase — the slope limiter changes that balance slightly, still well within margin either way).
- `tests/sod_shock_tube.rs` — compared pointwise against a from-scratch Rust port of Toro's exact Riemann solver (`tests/support/exact_riemann.rs`, self-tested against the textbook Sod case). Star region matches to 4 significant figures (32 Pa / 0.07 m/s); whole-profile RMS error (0.075 bar pressure, 8.2 m/s velocity) is concentrated at the shock/rarefaction fronts and now runs at essential parity with OpenWAM's own TVD-scheme result (0.068 bar / 6.8 m/s) — a first-order-only version of this solver measured roughly 2.5x worse on both (0.186 bar / 14.3 m/s), confirming the MUSCL-Hancock upgrade's benefit directly on this case. Whole-profile RMS in the low tens of Pa isn't a reachable target for any practical shock-capturing scheme — a genuine discontinuity is always smeared over a few cells no matter the scheme order, and refining the mesh only helps at roughly `sqrt(cell count)`; the *smooth* star-region result above (32 Pa) is where sub-100 Pa accuracy is actually achievable.

Junctions (phase 2 of the original plan) landed as part of phase 1, since both multi-pipe validation cases needed them. Remaining phases:

2. 0-D cylinder + Wiebe combustion + valve models, checked against the single-cylinder S54B32 OpenWAM case.
3. Multi-cylinder firing order + full 6-2-1 exhaust manifold (N-way branch junctions), checked against the 6-cylinder S54B32 OpenWAM case.
4. Wire the real solver into the web UI — `services/web/src/pages/SmokeTestPage.jsx` and `WavePage.jsx` still call the old placeholder API and show a broken/stale demo until this lands (chart design for wave profiles over time/crank-angle is a meaningfully different UI problem than a single sine curve, and deserves its own pass).

## Future cloud deployment

The web container is a single static SPA with no backend — deployable to any static host (Azure Static Web Apps, Netlify, GitHub Pages, S3+CloudFront, etc.) in addition to the Docker/nginx path above. No server-side compute, no COOP/COEP headers required (the wasm build is single-threaded for now — see `rust/crates/byglab-wasm`).

## Legacy files

Root `IntakeMachCalculator.jsx` and the incomplete `.crdownload` Python twin are superseded by the Rust solver + `services/web`. The prior Python prototype (`packages/byglab_engine`, `services/api`) has been retired — OpenWAM (`benchmarks/openwam/`) is the reference implementation the Rust port is validated against instead.
