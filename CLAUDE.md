# BYGLab

Open-source internal combustion engine dynamics simulator: a 1D/0D finite-volume + lumped-parameter gas-dynamics solver written in Rust, compiled to WebAssembly, running client-side in a browser (no server-side compute). Long-term goal: a full engine-modeling suite for intake manifold, exhaust header/muffler, camshaft, cylinder head porting, valve sizing/seat angle, and throttle body optimization studies.

Read `README.md` first — architecture overview, the full roadmap (phases 2-6 not yet built), and the S54B32 reference engine spec table. Read `rust/crates/byglab-core/src/lib.rs`'s module doc comment for the current solver's scope. Read `benchmarks/openwam/README.md` for how the OpenWAM reference cases were built/validated (input format quirks, known OpenWAM bugs worked around, exact measured numbers for every case).

`BLAIR-Gordon-P-Design-and-Simulation-of-Four-Stroke-Engine.pdf` (repo root) is the primary textbook reference for the physics in the upcoming phases (Wiebe combustion, valve discharge coefficients, wave action, muffler design) — check it before inventing a correlation from scratch.

## Rust toolchain environment (Windows, GNU toolchain — not MSVC)

C: drive has no space and no MSVC linker; the whole toolchain lives on D:. Every `cargo`/`rustc` invocation in this project needs:

```bash
export RUSTUP_HOME="/d/rust-toolchain/rustup"
export CARGO_HOME="/d/rust-toolchain/cargo"
export TEMP="/d/temp"
export TMP="/d/temp"
export PATH="/d/rust-toolchain/rustup/toolchains/stable-x86_64-pc-windows-gnu/lib/rustlib/x86_64-pc-windows-gnu/bin/self-contained:/d/rust-toolchain/cargo/bin:$PATH"
```

- Toolchain: `stable-x86_64-pc-windows-gnu` (an MSVC toolchain is also installed but unused — no linker for it).
- `wasm-pack` is NOT installed via cargo (compiling it from source hit a `dlltool`/`as.exe` gap in the self-contained GNU toolchain) — a prebuilt binary lives at `/d/rust-toolchain/wasm-pack-bin/wasm-pack.exe`; add its directory to `PATH` to use it.
- Run cargo commands from `rust/` (the workspace root: `crates/byglab-core`, `crates/byglab-cli`, `crates/byglab-wasm`).

## Project structure

- `benchmarks/openwam/` — legacy C++ OpenWAM (git submodule), built via Docker, used ONLY to generate reference numerical output (exact solutions, validated cases) to check the Rust port against. Not linked into the Rust code in any way — a clean-room reimplementation, not a port of the C++.
- `rust/` — the real solver. `byglab-core` (pure logic, no I/O, wasm-portable), `byglab-wasm` (thin wasm-bindgen wrapper), `byglab-cli` (native binary for manual debugging; `example-config.json` is a working sample config).
- `services/web/` — Vite/React frontend. Currently STALE/disconnected: `SmokeTestPage.jsx`/`WavePage.jsx`/`wasm-worker.js` still call an old placeholder sine-wave API, not the real pipe solver. This is intentionally deferred (see README roadmap phase 6), not a bug to fix opportunistically.
- `packages/`, root `IntakeMachCalculator.jsx`, `.venv` — retired Python prototype, superseded by the Rust solver. Don't build on these.

## Established conventions (demonstrated throughout this project, not just stated)

- **Every new physics piece gets validated against real ground truth** before being considered done: an exact analytical solution, a closed-form relation, or an OpenWAM reference case — never "looks plausible." Real measured numbers go in doc comments/README, never aspirational ones (if a target isn't reachable, say so and explain why — e.g. the whole-profile RMS-error discussion in `sod_shock_tube.rs`).
- **No backwards-compatibility shims.** Pre-1.0, breaking API renames (e.g. `PipeSpec.diameter` → `diameter_left`/`diameter_right`) are done as clean breaks, not aliased/deprecated.
- **Delete dead code** rather than leave it with `#[allow(dead_code)]` — e.g. `ConservedState::advance` was deleted outright once `apply_face_fluxes` stopped calling it.
- **No comments unless the WHY is non-obvious** (a hidden constraint, a subtle invariant, a workaround). Doc comments on public items explain what/why, not restate the signature.
- New solver capabilities get a design-review pass (a Plan/Explore subagent stress-testing the numerical approach against the actual current source) before implementation, when the physics is nontrivial — this caught a real Reynolds=0 NaN bug and a units double-division bug in one session, both before they shipped.

## Current status

Phase 1 (quasi-1D pipe finite-volume solver: MUSCL-Hancock + HLLC, taper, wall friction, wall heat transfer) is done and validated — see `README.md`'s "Solver porting roadmap" section for the full test list and measured numbers. Phase 2 (0D cylinder + combustion + valve flow + camshaft) is far along: `crank_mechanism.rs` (exact finite-rod-length slider-crank kinematics), `cylinder.rs` (cylinder volume + motoring/fired-cycle/breathing energy balance — three independent integration modes), `combustion.rs` (Wiebe heat release + Woschni wall heat transfer, validated against the real OpenWAM S54 2500rpm case — systematic ~30%/40% high on pressure/temperature, attributed to an unmodeled residual-gas fraction, not a physics/integration bug), `camshaft.rs` (versine lift profile), `valve.rs` (compressible orifice mass flow, correct choked/subsonic and reverse-flow handling), and now `valve_port.rs` (binding a valve to a real 1D `Pipe` end) are all done. The cylinder can "breathe" mass/energy either to/from a fixed reservoir (`integrate_breathing`, validated: exact linear growth under guaranteed-choked conditions to 3.3e-12 relative error, correct reverse-flow direction, independent trapezoidal cross-check to 5.5e-9) or now to/from a real, wave-propagating pipe (`valve_port::step_pipe_cylinder`/`run_pipe_cylinder_to_time`, validated: exact — 0.0 relative error — combined mass/energy conservation on a rigid cylinder over a short choked-flow window, a first-step analytic-mdot cross-check to ~1e-9, and soft sanity agreement with the fixed-reservoir path). Phases 3-6 (multi-cylinder + branched manifolds; muffler elements; performance metrics/sweeps; web UI wiring) are scoped in the README but not started.

Three real bugs were caught and fixed while building this: a curtain-area formula used `sin(seat_angle)` instead of `cos(seat_angle)` (invisible at exactly 45°, ~40% wrong elsewhere — caught by a pre-implementation design review); `crank_mechanism.rs`'s TDC-finding Newton solve divided `0.0/0.0` for a zero-stroke ("rigid") mechanism, poisoning everything downstream with NaN (caught by the first breathing test run, fixed with an explicit degenerate-case guard); and an early version of the pipe-coupling "exact conservation" test flagged a real (if honest) modeling confound rather than a code bug — a moving-piston cylinder's own compression/expansion work is a genuine external energy source unrelated to valve/pipe exchange conservation, so that check needs a rigid cylinder to isolate the property it's actually meant to test.

`camshaft_presets.rs` now holds 10 REAL BMW S54 Schrick camshaft grinds (duration/lift/installed timing for every published spread option, low/medium/high street, intake and exhaust) ingested directly from official Schrick PDF data sheets — not estimated. Important subtlety if you extend this: Schrick's data sheets reference crank angle "0°" to *gas-exchange* TDC, not the firing TDC `combustion.rs` uses — a ±360° shift is needed before combining this data with the combustion model in one simulation. Checking the versine approximation against the manufacturer's own `@1mm` timing across all 10 grinds shows a consistent, honest 7-12° discrepancy (real cam lobes ramp up faster than a symmetric raised-cosine) — documented, not hidden; motivates an eventual direct-lift-table `CamProfile` variant if higher fidelity is ever needed.

See the plan file (`agile-wondering-grove.md` in the plans directory) for the design writeup of the pipe-coupling piece just completed. Remaining work, roughly in order: (1) unifying `integrate_motoring`/`integrate_fired_cycle`/`integrate_breathing` (plus the new pipe-coupled path) into one full-cycle ODE, (2) reconciling the crank-angle convention mismatch between `combustion.rs` (firing TDC = 0) and the ingested Schrick cam data (gas-exchange TDC = 0) — item 1 naturally wants item 2 solved at the same time. Residual-gas combustion calibration and CAT Cams data ingestion are lower-priority, independent refinements, not blockers. Key scoping decision already made: valve discharge coefficients are a parametric INPUT to the model, not predicted from 3D port geometry — this is what makes porting/seat-angle/valve-size studies tractable, matching how professional 1D tools (GT-Power, WAVE) work.

## Running tests

```bash
cd rust
cargo test -p byglab-core --release -- --nocapture   # release mode matters — the heat-transfer relaxation test runs many thousands of explicit timesteps
cargo build --workspace && cargo test --workspace
cargo run -p byglab-cli -- crates/byglab-cli/example-config.json
```
