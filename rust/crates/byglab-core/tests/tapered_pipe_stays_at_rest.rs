//! The well-balanced correctness check for the taper (variable
//! cross-sectional area) geometric source term: a uniform, at-rest gas
//! column in a *tapered* pipe must stay exactly as it started, exactly
//! like the constant-diameter case in `quiescent.rs`.
//!
//! This is the critical test for the source term added in `pipe.rs`'s
//! `apply_face_fluxes` (`p̄ * (A_right_face - A_left_face)`, momentum only):
//! a naive or inconsistent discretization of the `p * dA/dx` term would
//! generate spurious velocity here even though nothing is physically
//! driving the gas — reducing the momentum equation for a steady, u≡0 duct
//! with no body force gives `dp/dx = 0` unconditionally, so uniform
//! pressure is the only at-rest equilibrium, independent of how the area
//! varies along the pipe.
//!
//! Mechanistically: for a uniform state, MUSCL reconstruction gives zero
//! slope everywhere, so every face's two approaching states are identical
//! and HLLC returns the exact analytic flux there (already covered by
//! `riemann.rs`'s own flux-consistency tests) — so the area-weighted
//! momentum-flux divergence across each cell is exactly
//! `p0 * (A_right_face - A_left_face)`, which exactly cancels the
//! geometric source term as long as it uses the same `p0` (which it does,
//! by construction: `p̄` is the average of the cell's own two
//! MUSCL-reconstructed face pressures, both exactly `p0` here).

use byglab_core::boundary::BoundaryCondition;
use byglab_core::gas::{GasProperties, PrimitiveState};
use byglab_core::mesh::Mesh;
use byglab_core::network::PipeNetwork;
use byglab_core::pipe::Pipe;
use byglab_core::solver::run_to_time;

#[test]
fn uniform_at_rest_tapered_pipe_stays_at_rest() {
    let gas = GasProperties::AIR;
    // A 2x taper over the pipe's length (0.05 m to 0.10 m diameter, a 4x
    // area ratio) - aggressive enough that a well-balancedness bug would
    // produce an obvious, not marginal, spurious velocity.
    let mesh = Mesh::tapered(1.0, 0.05, 0.10, 0.01);
    let initial_state = PrimitiveState::from_pressure_temperature(1.5e5, 293.15, 0.0, &gas);
    let pipe =
        Pipe::uniform_initial_state(mesh, initial_state, &gas, BoundaryCondition::ClosedEnd, BoundaryCondition::ClosedEnd);

    let mut network = PipeNetwork::single_pipe(pipe);
    run_to_time(&mut network, &gas, 0.5, 0.02);

    let mut max_pressure_deviation = 0.0_f64;
    let mut max_speed = 0.0_f64;
    for state in &network.pipes[0].state {
        let primitive = state.to_primitive(&gas);
        max_pressure_deviation = max_pressure_deviation.max((primitive.pressure - 1.5e5).abs());
        max_speed = max_speed.max(primitive.velocity.abs());
    }

    assert!(max_pressure_deviation < 1e-3, "pressure drifted by {max_pressure_deviation} Pa");
    assert!(max_speed < 1e-6, "spurious velocity of {max_speed} m/s appeared from taper alone - the geometric source term is not well-balanced");
}
