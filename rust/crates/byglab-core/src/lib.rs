//! `byglab-core` — a 1D finite-volume gas-dynamics solver for engine intake
//! and exhaust runners, ported (as a clean-room reimplementation, not a
//! line-by-line translation) from the physics validated in this project's
//! OpenWAM reference suite (`benchmarks/openwam/`).
//!
//! # What's here
//!
//! A second-order MUSCL-Hancock finite-volume scheme (HLLC approximate
//! Riemann solver fed slope-limited, half-timestep-evolved reconstructed
//! states) for the quasi-1D compressible Euler equations — variable
//! cross-sectional area (taper), wall friction, and wall heat transfer are
//! all supported — over networks of pipes joined by same-area junctions.
//! Phase 2 (0D cylinder + combustion + valve flow + camshaft) is underway
//! — [`crank_mechanism`] (piston kinematics), [`cylinder`] (volume +
//! motoring/fired-cycle/breathing energy balance), [`combustion`] (Wiebe
//! heat release + Woschni wall heat transfer), [`camshaft`] (valve lift
//! profile), [`valve`] (compressible orifice mass flow), and
//! [`valve_port`] (binding a valve to a real 1D [`pipe::Pipe`] end, so a
//! cylinder can breathe against genuine wave-propagating flow instead of
//! only a fixed reservoir) are done. [`branch_junction`] generalizes
//! [`network::Junction`]'s 2-pipe, same-area connection to any number of
//! pipes with any areas (a single intake/exhaust runner splitting into
//! several, or several merging into one). Multi-cylinder firing order is
//! not implemented yet (see the root `README.md`'s roadmap section).
//!
//! # Module map
//!
//! - [`gas`] — ideal-gas thermodynamics: the state types ([`gas::PrimitiveState`],
//!   [`gas::ConservedState`], [`gas::Flux`]) everything else is built on.
//! - [`mesh`] — pipe discretization into finite-volume cells, including taper.
//! - [`reconstruction`] — MUSCL slope-limited reconstruction and the
//!   MUSCL-Hancock half-timestep predictor.
//! - [`riemann`] — the HLLC numerical flux function.
//! - [`source_terms`] — wall friction and wall heat transfer
//!   ([`source_terms::WallProperties`]); the taper geometric source term
//!   lives in `pipe.rs` itself (it needs each cell's face areas and
//!   reconstructed pressures, both already local to that module).
//! - [`boundary`] — how a pipe's end is terminated ([`boundary::BoundaryCondition`]).
//! - [`pipe`] — a single pipe: mesh + gas state + boundary conditions + wall properties.
//! - [`network`] — multiple pipes joined by [`network::Junction`]s.
//! - [`branch_junction`] — an N-way [`branch_junction::BranchJunction`]
//!   generalizing `Junction` to any number of pipes/areas, via a
//!   method-of-characteristics closure (each leg's own Riemann invariant
//!   plus a shared junction pressure, solved by bisection). Validated:
//!   exact structural mass conservation, a symmetric N-branch split/merge
//!   matching the closed-form linear-acoustics junction-pressure
//!   prediction (error shrinking ~linearly with disturbance amplitude), a
//!   quantified (not just sanity-checked) N=2 regression against
//!   `Junction`'s own exact HLLC solve, and a measured (not assumed-zero)
//!   energy-conservation residual (~0.14% small-amplitude, ~15.8% at a
//!   realistic exhaust-pulse ratio) — two independently-wrong closures were
//!   tried and disproved by direct numerical scans before landing on this
//!   one; see the module's own doc comment for what failed and why.
//! - [`lossy_branch_junction`] — closes that energy gap for the specific
//!   case Gordon Blair worked out in closed form: a 3-way (only) branch
//!   with a real angle-dependent pressure loss and a genuinely separate,
//!   jointly-solved stagnation-enthalpy energy equation
//!   ([`lossy_branch_junction::LossyBranchJunction`]/
//!   [`lossy_branch_junction::resolve_lossy_branch_junction`]). The
//!   one-supplier/two-supplied case matches all 4 of Blair's own published
//!   test cases; the two-supplier/one-supplied case has no published
//!   reference and is validated only for internal consistency. Rerunning
//!   `branch_junction`'s own "realistic exhaust pulse" scenario through
//!   this model drops the energy residual from ~15.8% to ~3e-7.
//! - [`solver`] — the explicit time-stepping driver.
//! - [`case`] — the serializable public API ([`case::PipeCaseConfig`]/
//!   [`case::PipeCaseResult`]/[`case::run_pipe_case`]) that `byglab-cli` and
//!   `byglab-wasm` bind to.
//! - [`crank_mechanism`] — slider-crank kinematics ([`crank_mechanism::CrankMechanism`]):
//!   exact (finite rod length, optional piston pin offset) piston
//!   position/velocity/acceleration as a function of crank angle.
//! - [`cylinder`] — 0D cylinder volume ([`cylinder::Cylinder`]) and
//!   thermodynamic state ([`cylinder::CylinderState`]), with both a
//!   motoring-only energy balance (validated against the exact isentropic
//!   relation) and a fired-cycle energy balance (combustion + wall heat
//!   transfer, validated against a real OpenWAM reference case).
//! - [`combustion`] — Wiebe mass-fraction-burned heat release and Woschni
//!   in-cylinder wall heat transfer; the source terms `cylinder`'s
//!   fired-cycle integration adds to the motoring energy balance.
//! - [`camshaft`] — valve lift as a function of crank angle
//!   ([`camshaft::CamProfile`]).
//! - [`camshaft_presets`] — real BMW S54 camshaft grind data ingested
//!   directly from official Schrick data sheets (duration, lift, and
//!   installed timing for 10 grinds spanning low/medium/high street
//!   intake and exhaust), not estimated.
//! - [`valve`] — poppet valve curtain area and quasi-steady compressible
//!   orifice mass flow rate ([`valve::mass_flow_rate`]); the source term
//!   `cylinder`'s breathing integration adds to the motoring energy
//!   balance (mass exchange with an external reservoir/pipe, not yet
//!   combined with combustion/wall heat transfer in the same integration).
//! - [`valve_port`] — binds a [`valve::ValveGeometry`]/[`camshaft::CamProfile`]
//!   pair to one end of a real [`pipe::Pipe`] ([`valve_port::ValvePort`]),
//!   and drives a [`network::PipeNetwork`] and a [`cylinder::Cylinder`]
//!   together one shared timestep at a time
//!   ([`valve_port::step_pipe_cylinder`]/[`valve_port::run_pipe_cylinder_to_time`]).
//!   Validated: exact (to floating-point precision) combined mass/energy
//!   conservation over a short choked-flow window, an exact first-step
//!   cross-check against the fixed-reservoir path's own analytic mass flow
//!   rate, and the Left/Right sign convention.
//!
//! # Design for WebAssembly
//!
//! This crate has no file or network I/O, no threads, and no `wasm-bindgen`
//! dependency of its own — every public type is plain data (`Vec<f64>`-
//! backed, `serde`-serializable), and [`case::run_pipe_case`] is a pure
//! function from config to result. `byglab-wasm` wraps this crate with
//! nothing but marshalling; `byglab-core` itself compiles and tests
//! identically on native and `wasm32-unknown-unknown` targets.
//!
//! # Validation
//!
//! `tests/` checks this solver against real ground truth throughout: the
//! same analytically-validated cases used to validate OpenWAM itself
//! (`benchmarks/openwam/cases/{quiescent,acoustic_resonance,sod_shock_tube}/`
//! — a trivial rest-state consistency check, a closed-form acoustic
//! resonance period, and a nonlinear Riemann problem with a known exact
//! solution), an empirical mesh-refinement convergence-order check, and —
//! for the quasi-1D extensions — a well-balanced at-rest check for the
//! taper source term, a closed-form isentropic area-Mach-relation check
//! for tapered flow, and a closed-form exponential relaxation check for
//! wall heat transfer. See the root `README.md`'s roadmap section for the
//! actual measured numbers.

pub mod boundary;
pub mod branch_junction;
pub mod camshaft;
pub mod camshaft_presets;
pub mod case;
pub mod combustion;
pub mod crank_mechanism;
pub mod cylinder;
pub mod gas;
pub mod lossy_branch_junction;
pub mod mesh;
pub mod network;
pub mod pipe;
pub mod reconstruction;
pub mod riemann;
pub mod solver;
pub mod source_terms;
pub mod valve;
pub mod valve_port;

pub use boundary::BoundaryCondition;
pub use branch_junction::{resolve_branch_junction, BranchJunction};
pub use camshaft::CamProfile;
pub use case::{run_pipe_case, PipeCaseConfig, PipeCaseResult, PipeResult, PipeSpec};
pub use crank_mechanism::CrankMechanism;
pub use cylinder::{Cylinder, CylinderState};
pub use gas::{ConservedState, Flux, GasProperties, PrimitiveState};
pub use lossy_branch_junction::{resolve_lossy_branch_junction, LossyBranchJunction};
pub use mesh::{Cell, Mesh};
pub use network::{ExternalPortFlux, Junction, PipeEnd, PipeEndRef, PipeNetwork};
pub use pipe::Pipe;
pub use source_terms::WallProperties;
pub use valve::ValveGeometry;
pub use valve_port::{run_pipe_cylinder_to_time, step_pipe_cylinder, ValvePort};
