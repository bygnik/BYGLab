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

**Phase 1 — done.** `byglab-core` is a real quasi-1D finite-volume gas-dynamics solver: second-order MUSCL-Hancock scheme (piecewise-linear reconstruction with a minmod slope limiter, HLLC approximate Riemann solver, explicit CFL-limited time-stepping — see `reconstruction.rs`), supporting networks of pipes joined by same-area junctions with fully symmetric second-order accuracy at the junction faces too, variable cross-sectional area (taper — velocity stacks, diffusers, header collector cones), and optional wall friction/heat transfer (`source_terms.rs`). A clean-room Rust reimplementation informed by OpenWAM's validated physics, not a line-by-line C++ translation. Started as a simpler first-order, constant-area, frictionless/adiabatic scheme, then upgraded twice once each stage was validated — all stages checked against real ground truth:

- `tests/quiescent.rs` — uniform at-rest closed pipe stays at rest (no discretization error possible to hide behind).
- `tests/acoustic_resonance.rs` — closed-closed resonance period matches the exact `T = 2L/c0` to **0.09%** (OpenWAM's own higher-order scheme: 0.07% — comparable; a first-order-only version of this solver measured 0.002% here, tighter still, since numerical dissipation for a small linear disturbance mostly damps amplitude rather than shifting phase — the slope limiter changes that balance slightly, still well within margin either way).
- `tests/sod_shock_tube.rs` — compared pointwise against a from-scratch Rust port of Toro's exact Riemann solver (`tests/support/exact_riemann.rs`, self-tested against the textbook Sod case). Star region matches to 4 significant figures (32 Pa / 0.07 m/s); whole-profile RMS error (0.075 bar pressure, 8.2 m/s velocity) is concentrated at the shock/rarefaction fronts and now runs at essential parity with OpenWAM's own TVD-scheme result (0.068 bar / 6.8 m/s) — a first-order-only version of this solver measured roughly 2.5x worse on both (0.186 bar / 14.3 m/s), confirming the MUSCL-Hancock upgrade's benefit directly on this case. Whole-profile RMS in the low tens of Pa isn't a reachable target for any practical shock-capturing scheme — a genuine discontinuity is always smeared over a few cells no matter the scheme order, and refining the mesh only helps at roughly `sqrt(cell count)`; the *smooth* star-region result above (32 Pa) is where sub-100 Pa accuracy is actually achievable.
- `tests/mesh_convergence.rs` — doubling mesh resolution on the Sod case shrinks error at close to the scheme's design rate: smooth star-region error shrinks at order 3.9 (pressure) / 2.6 (velocity) — even *exceeding* the formal 2nd-order rate, since that probe's only error source is decaying numerical diffusion bleeding in from the distant shock, not classical truncation error — while whole-profile RMS error (which includes the shock/rarefaction fronts) shrinks at a reduced but real order 0.55 / 0.58, the expected signature of discontinuity smearing capping the rate below the smooth-region behavior. A meaningfully stronger correctness signal than a fixed-tolerance check alone — a subtly broken scheme could still look "small enough" at one resolution while failing to converge under refinement.
- `tests/tapered_pipe_stays_at_rest.rs` — the well-balanced check for the taper geometric source term: a uniform, at-rest gas column in a pipe tapering from ⌀0.05m to ⌀0.10m generates **zero** spurious velocity (< 1e-6 m/s), proven algebraically exact (not just empirically small) by construction — the source term reuses the same MUSCL-reconstructed face pressures already computed for the flux, rather than a separately-averaged value, so the two cancel exactly for any uniform state regardless of how the area varies.
- `tests/isentropic_nozzle.rs` — steady subsonic flow through a converging duct, calibrated from one interior station's own state and checked against every other interior station via the exact isentropic area-Mach relation (`tests/support/isentropic_nozzle.rs`, a from-scratch bisection solver, self-tested against a textbook A/A* table value). Interior pressure matches to **0.11%** of the driving pressure range, comfortably inside the 0.2% target — validates that the taper source term correctly reproduces quasi-1D isentropic flow acceleration under genuine (non-zero-velocity) flow, isolated from the separately-documented approximation in the `Reservoir` boundary condition (see `boundary.rs`).
- `tests/shu_osher.rs` — the classic Shu-Osher shock/entropy-wave interaction problem (a Mach-3 shock advancing into a sinusoidal density field, `gamma=1.4`, standard non-dimensional problem values, run to `t=1.8`). Checks that pressure stays at its initial uniform value in the region the shock has not yet reached — measured deviation was **exactly `0.0`**, bit-for-bit, not just "small." A stronger property than it first looks: since velocity and pressure are both spatially uniform ahead of the shock while density varies sinusoidally, and this solver's MUSCL reconstruction limits each primitive variable independently (no cross-variable coupling), every face flux ahead of the shock is *exactly* constant regardless of the density gradient, so nothing should change there at all until the shock's actual domain of dependence reaches it — a real correctness signature, not just "small enough," and the measured result confirms the scheme has no hidden cross-variable coupling.
- `tests/heat_transfer_relaxation.rs` — a quiescent gas column with a fixed differing wall temperature relaxes toward it exactly along the closed-form exponential time constant `τ = ρ·cv·D/(4h)`, matching to **0.0003 K** out of a ~107 K total temperature rise. Also regression coverage for a zero-Reynolds-number singularity caught during design review (naive friction/heat-transfer formulas are singular at u=0 — routine in pulsating engine flow — which would otherwise produce `NaN` instead of the physically correct zero).
- `source_terms.rs`'s own unit tests check the Haaland friction factor and Colburn Nusselt-number correlations against Moody-chart/textbook reference points.
- Mass conservation in a closed pipe with a moving pressure pulse (mass has no source term — flux differencing conserves it exactly regardless of scheme details, so this is a bookkeeping check, not a physics one): **0.001%** for a constant-diameter pipe (`solver.rs`'s own unit test) and **0.000145%** for a tapered pipe (`tests/mass_conservation_taper.rs`, checking the area-weighted flux differencing specifically — a genuinely different code path added for taper support), both far inside the 0.1% target.

Wall friction (Haaland correlation) and heat transfer (Colburn correlation, fixed wall temperature) are simpler standard choices deliberately *not* ported from OpenWAM's own piecewise-polynomial Colebrook cascade, regime-switched Nusselt correlations, and full multi-layer transient wall-conduction model (`benchmarks/openwam/OpenWAM/Source/1DPipes/TTubo.cpp`) — same governing physics, independently re-derived and cross-checked against OpenWAM's formulas, simpler code. A junction between two mismatched-area pipe ends (a genuine sudden area change, needing a different loss-coefficient-based treatment) is explicitly out of scope and rejected at runtime with a clear panic.

Junctions (phase 2 of the original plan) landed as part of phase 1, since both multi-pipe validation cases needed them. Remaining phases, scoped against the end goal of a full engine-modeling suite (intake manifold, exhaust header/muffler, camshaft, cylinder head porting, valve sizing/seat angle, and throttle body optimization studies):

2. **0-D cylinder + Wiebe combustion + valve flow + camshaft**, checked against the single-cylinder S54B32 OpenWAM case. Cylinder: crank-slider `V(θ)`, single/double Wiebe heat release, Woschni-type wall heat transfer. Valve flow: quasi-steady compressible orifice flow through the curtain area (the same subsonic/choked-flow physics already validated by `tests/isentropic_nozzle.rs`), scaled by a discharge coefficient `Cd(lift/diameter)`. **Scoping decision: `Cd` curves are a parametric *input*, not predicted from 3D port geometry** — this is how professional 1D tools (GT-Power, WAVE) handle it too, and it's what makes "cylinder head porting," "valve seat angle," and "valve size" studies tractable: you compare different `Cd` curves (measured on a flow bench, or estimated), not run port-level CFD inside the cycle solver. Camshaft: a lift-vs-crank-angle profile (lift/duration/timing/ramps or a direct table) per valve, feeding the curtain area; must support reverse flow (blowdown, overlap reversion), since that's the actual physics of exhaust pulse tuning. Throttle body reuses this same orifice-flow submodel with a `Cd(throttle angle)` curve instead.
   - **Crank mechanism — done.** `crank_mechanism.rs`: exact (finite connecting-rod-length, not the infinite-rod/simple-harmonic-motion approximation) slider-crank kinematics — piston position, velocity, and acceleration as functions of crank angle, with piston pin offset ("desaxé") support. Offset shifts true TDC away from the "crank pin collinear with the cylinder axis" reference — a real second-order effect, solved exactly via a one-time Newton iteration at construction (not a small-angle series approximation) so every public method's crank angle is measured from *true* TDC, matching how a physical TDC sensor is calibrated. Validated against the standard textbook slider-crank position/velocity formulas (exact match, no offset case), central finite differences of the solver's own position/velocity (velocity and acceleration respectively), the simple-harmonic-motion limit as rod length → ∞, and an exact identity (piston velocity is exactly zero at true TDC and true BDC, with or without offset — confirmed to < 1e-9 m/s). One genuine finding along the way: with offset, total TDC-to-BDC piston travel is *not* exactly the nominal stroke (TDC and BDC each get an independent, generally unequal phase correction) — measured ~0.003-0.007% deviation for realistic 1-1.5mm offsets, a real physical effect now documented rather than assumed away.
   - **Motoring cylinder model — done.** `cylinder.rs`: cylinder volume `V(θ)` built on `CrankMechanism`, plus a first-law energy balance (`dU = -p dV`) for the trapped charge, integrated in crank-angle space (not real time — a reversible adiabatic process is exactly rate-independent) via classic 4th-order Runge-Kutta. Deliberately scoped to *motoring only* (no combustion, no valve mass flow, no wall heat transfer yet) because that has an exact closed-form answer to validate against before adding physics that don't: the isentropic relation `p·V^γ = const`. Using real S54B32 geometry (CR 11.5:1), compressing from BDC to TDC matches the exact isentropic pressure and temperature to **~1e-11 relative error** — essentially machine precision, not just "small." Mass conservation is exact (no valve flow exists in this model yet, so there's nothing to leak). A full BDC→TDC→BDC round trip — compression then expansion retracing the same reversible path — returns to the exact initial energy and pressure to **~6.6e-14 relative error**, confirming the numerical integration doesn't manufacture or leak energy over a closed cycle. Halving the RK4 step size shrinks the error at **observed order 4.26**, confirming genuine 4th-order convergence rather than a scheme that merely looks accurate at one step count. Combustion (Wiebe heat release) and wall heat transfer (Woschni correlation) are the natural next additions to the same energy-balance ODE.
3. **Multi-cylinder firing order + branched intake/exhaust manifolds** (N-way junctions, generalizing today's 2-pipe `Junction`), checked against the 6-cylinder S54B32 OpenWAM case. Firing order sets the relative phase of each cylinder's pulses arriving at a shared plenum/collector — that phase relationship is the entire physics of intake plenum and header/collector design.
4. **Muffler/silencer elements** — a genuinely different category from a tapered pipe: expansion chambers (sudden area change — the case explicitly deferred in phase 1's taper work, see above), perforated-tube sections, and side-branch resonators (Helmholtz/quarter-wave), each needing its own loss-coefficient/perforate-impedance model. Bundle in a revisit of `BoundaryCondition::Reservoir`'s simplified open-end treatment (see its doc comment) for proper reflection-coefficient accuracy — matters for acoustic muffler tuning in a way it didn't for the runner-only validation so far.
5. **Performance metrics + operating-point sweeps** — none of the above studies are judged from raw pressure/velocity/temperature traces directly; they're judged by volumetric efficiency, torque/power vs. RPM, trapped mass, pumping losses, and exhaust backpressure, derived automatically rather than read off a chart by hand. Needs support for running a single geometry across many RPM/load points (a torque curve is many solver runs, not one) — shapes `case.rs`'s API and the browser-side WASM compute budget for an optimization loop.
6. Wire the real solver into the web UI — `services/web/src/pages/SmokeTestPage.jsx` and `WavePage.jsx` still call the old placeholder API and show a broken/stale demo until this lands (chart design for wave profiles over time/crank-angle is a meaningfully different UI problem than a single sine curve, and deserves its own pass).

## Future cloud deployment

The web container is a single static SPA with no backend — deployable to any static host (Azure Static Web Apps, Netlify, GitHub Pages, S3+CloudFront, etc.) in addition to the Docker/nginx path above. No server-side compute, no COOP/COEP headers required (the wasm build is single-threaded for now — see `rust/crates/byglab-wasm`).

## Legacy files

Root `IntakeMachCalculator.jsx` and the incomplete `.crdownload` Python twin are superseded by the Rust solver + `services/web`. The prior Python prototype (`packages/byglab_engine`, `services/api`) has been retired — OpenWAM (`benchmarks/openwam/`) is the reference implementation the Rust port is validated against instead.
