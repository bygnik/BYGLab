//! Trivial consistency check: a uniform, at-rest, closed-closed pipe must
//! stay exactly as it started, for as long as it runs. No discretization
//! error is possible to hide behind here — any spurious motion or pressure
//! drift would indicate a real bug in the flux/boundary-condition
//! implementation. Matches `benchmarks/openwam/cases/quiescent/`, where
//! OpenWAM itself achieved pressure exact to the last printed digit and
//! velocity at ~1e-14 m/s (floating-point noise).

use byglab_core::boundary::BoundaryCondition;
use byglab_core::gas::{GasProperties, PrimitiveState};
use byglab_core::mesh::Mesh;
use byglab_core::network::PipeNetwork;
use byglab_core::pipe::Pipe;
use byglab_core::solver::run_to_time;

#[test]
fn uniform_at_rest_pipe_stays_at_rest() {
    let gas = GasProperties::AIR;
    let mesh = Mesh::uniform(1.0, 0.05, 0.01);
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
    assert!(max_speed < 1e-6, "spurious velocity of {max_speed} m/s appeared from nothing");
}
