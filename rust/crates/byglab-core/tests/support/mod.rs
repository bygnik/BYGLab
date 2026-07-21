//! Shared test-only support code, used by more than one integration test
//! file. Not itself a test binary (Cargo only auto-discovers files placed
//! directly in `tests/`, not in subdirectories) — each test file that
//! needs this declares `mod support;` to pull it in.

pub mod exact_riemann;
pub mod isentropic_nozzle;
pub mod sod_shock_tube_case;
