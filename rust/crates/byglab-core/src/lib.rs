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
//! This is the foundation phase of a larger effort — 0D cylinder +
//! combustion + valve models, multi-cylinder firing order, and branched
//! exhaust manifolds are not implemented yet (see the root `README.md`'s
//! roadmap section).
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
//! - [`solver`] — the explicit time-stepping driver.
//! - [`case`] — the serializable public API ([`case::PipeCaseConfig`]/
//!   [`case::PipeCaseResult`]/[`case::run_pipe_case`]) that `byglab-cli` and
//!   `byglab-wasm` bind to.
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
pub mod case;
pub mod gas;
pub mod mesh;
pub mod network;
pub mod pipe;
pub mod reconstruction;
pub mod riemann;
pub mod solver;
pub mod source_terms;

pub use boundary::BoundaryCondition;
pub use case::{run_pipe_case, PipeCaseConfig, PipeCaseResult, PipeResult, PipeSpec};
pub use gas::{ConservedState, Flux, GasProperties, PrimitiveState};
pub use mesh::{Cell, Mesh};
pub use network::{Junction, PipeEnd, PipeEndRef, PipeNetwork};
pub use pipe::Pipe;
pub use source_terms::WallProperties;
