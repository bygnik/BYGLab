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

Phase 1 (quasi-1D pipe finite-volume solver: MUSCL-Hancock + HLLC, taper, wall friction, wall heat transfer) is done and validated — see `README.md`'s "Solver porting roadmap" section for the full test list and measured numbers. Phase 2 (0D cylinder + combustion + valve flow + camshaft) is underway: `crank_mechanism.rs` (exact finite-rod-length slider-crank kinematics, piston pin offset support) and `cylinder.rs` (cylinder volume + motoring energy balance, RK4-integrated, validated against the exact isentropic relation to ~1e-11 relative error) are both done. No combustion, wall heat transfer, or valve flow yet. Phases 3-6 (multi-cylinder + branched manifolds; muffler elements; performance metrics/sweeps; web UI wiring) are scoped in the README but not started.

Suggested next step (discussed, not started): add Wiebe combustion heat release and Woschni wall heat transfer to `cylinder.rs`'s energy-balance ODE (same RK4 integration, new source terms), then valve flow (quasi-steady compressible orifice flow, `Cd(lift/diameter)` curves) to connect the cylinder to the intake/exhaust runners. Key scoping decision already made: valve discharge coefficients are a parametric INPUT to the model, not predicted from 3D port geometry — this is what makes porting/seat-angle/valve-size studies tractable, matching how professional 1D tools (GT-Power, WAVE) work.

## Running tests

```bash
cd rust
cargo test -p byglab-core --release -- --nocapture   # release mode matters — the heat-transfer relaxation test runs many thousands of explicit timesteps
cargo build --workspace && cargo test --workspace
cargo run -p byglab-cli -- crates/byglab-cli/example-config.json
```
