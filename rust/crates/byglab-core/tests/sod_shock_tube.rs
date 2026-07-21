//! The flagship validation case: two 1m pipes joined at a shared node, a
//! strong pressure ratio (10 bar vs 1 bar) generating a genuine shock,
//! contact discontinuity, and rarefaction fan — compared pointwise against
//! Toro's exact Riemann solution (`tests/support/exact_riemann.rs`, itself
//! validated against the textbook Sod test).
//!
//! Matches `benchmarks/openwam/cases/sod_shock_tube/`. OpenWAM's
//! higher-order TVD scheme achieved RMS pressure error 0.068 bar and RMS
//! velocity error 6.8 m/s on this same case. This solver uses the same
//! class of scheme — MUSCL-Hancock reconstruction with a minmod limiter
//! (see `reconstruction.rs`) — and measures essentially at parity: RMS
//! pressure error 0.075 bar (max 0.55 bar), RMS velocity error 8.2 m/s
//! (max 91.6 m/s). (An earlier first-order-only version of this solver,
//! before the MUSCL-Hancock upgrade, measured roughly 2.5x worse on both —
//! 0.186 bar / 14.3 m/s RMS — confirming the upgrade's benefit directly on
//! this exact case.) The star region itself (p* = 2.8482 bar,
//! u* = 281.84 m/s) matches to 4 significant figures either way — remaining
//! error is concentrated at the shock/rarefaction fronts, which no
//! practical shock-capturing scheme at reasonable mesh resolution
//! eliminates entirely — a genuine mathematical discontinuity is always
//! smeared over at least a couple of cells, and the local error in those
//! cells is a real fraction of the actual jump (hundreds of thousands of
//! Pa here) regardless of scheme order. Getting whole-profile RMS down to,
//! say, 10 Pa isn't a reachable target for any practical FV scheme on a
//! problem with a genuine shock — even a much finer mesh only improves
//! this at roughly `sqrt(cell count)`, not linearly.
//!
//! See `tests/mesh_convergence.rs` for a check that error actually shrinks
//! at the rate a correctly-implemented second-order scheme should under
//! mesh refinement — a stronger correctness signal than the fixed-
//! tolerance checks here, which only confirm the error is small enough at
//! one specific resolution.

mod support;
use support::sod_shock_tube_case;

#[test]
fn matches_the_exact_riemann_solution_before_reflections_arrive() {
    let run = sod_shock_tube_case::run(0.01);

    println!("exact star region: p* = {:.4} bar, u* = {:.2} m/s", run.exact.p_star / 1e5, run.exact.u_star);

    // A sharper, more direct check than the whole-profile RMS below: the
    // solver's own state well inside the star region (constant p*/u* on
    // both sides of the contact, away from both the rarefaction and the
    // shock) should match the exact star state closely - most of the
    // whole-profile error is concentrated at the two moving fronts, not
    // spread through the plateau in between.
    let star_pressure_error = run.star_region_pressure_error();
    let star_velocity_error = run.star_region_velocity_error();
    println!("star-region errors: {star_pressure_error:.0} Pa, {star_velocity_error:.2} m/s");

    // Measured at 32 Pa / 0.07 m/s with MUSCL-Hancock (was 68 Pa / 0.24 m/s
    // with the earlier first-order-only version) - tight margins below
    // since this is the smooth, unsmeared part of the flow where accuracy
    // genuinely can be very high, unlike the whole-profile bound below.
    assert!(star_pressure_error < 500.0, "star-region pressure off by {star_pressure_error:.0} Pa");
    assert!(star_velocity_error < 1.0, "star-region velocity off by {star_velocity_error:.2} m/s");

    let pressure_error = run.whole_profile_pressure_error();
    let velocity_error = run.whole_profile_velocity_error();
    println!(
        "pressure error: RMS {:.0} Pa, max {:.0} Pa (state range 1-10 bar)",
        pressure_error.rms, pressure_error.max
    );
    println!(
        "velocity error: RMS {:.2} m/s, max {:.2} m/s (state range 0-282 m/s)",
        velocity_error.rms, velocity_error.max
    );

    // MUSCL-Hancock still smears the shock/rarefaction fronts over a
    // handful of cells (no shock-capturing scheme avoids this entirely);
    // these bounds leave comfortable margin above the measured values
    // (RMS 7,468 Pa / 8.19 m/s) while still catching a real regression.
    assert!(pressure_error.rms < 1.5e4, "RMS pressure error too high: {:.0} Pa", pressure_error.rms);
    assert!(velocity_error.rms < 15.0, "RMS velocity error too high: {:.2} m/s", velocity_error.rms);
}
