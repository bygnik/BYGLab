//! Closed-form validation of the wall heat-transfer source term
//! (`source_terms::wall_sources`'s energy component, applied by `pipe.rs`).
//!
//! A quiescent gas column in a closed pipe with a fixed wall temperature
//! has u≡0 for all time: heating is spatially uniform (every cell has the
//! same state, so the same source), so it raises pressure uniformly with
//! no pressure gradient to drive flow (the same mechanism that keeps
//! `quiescent.rs`/`tapered_pipe_stays_at_rest.rs` at rest — a uniform state
//! produces zero net flux at every face). With u≡0 for all time, Reynolds
//! number is exactly 0 throughout, so the Nusselt number stays pinned at
//! its laminar/Re=0 floor (3.66, see `source_terms.rs`) rather than
//! varying — giving a genuinely linear relaxation ODE,
//! `rho*cv*dT/dt = (4h/D)*(T_wall - T)`, with a closed-form exponential
//! solution `T(t) = T_wall - (T_wall - T0)*exp(-t/tau)`,
//! `tau = rho*cv*D/(4h)`. This doubles as regression coverage for the
//! Reynolds-number-zero friction/heat-transfer guard (`source_terms.rs`'s
//! `MIN_REYNOLDS_NUMBER_FOR_FRICTION`) — this scenario is u≡0 for its
//! entire duration.
//!
//! (Applying an internal-pipe-flow Nusselt correlation to genuinely
//! stagnant gas is a modeling simplification, not a claim that this is the
//! most physically realistic treatment of pure conduction to a stagnant
//! gas column — this test validates that the *implementation* correctly
//! integrates the ODE its own source term implies, not the choice of
//! correlation itself.)

use byglab_core::boundary::BoundaryCondition;
use byglab_core::gas::{GasProperties, PrimitiveState};
use byglab_core::mesh::Mesh;
use byglab_core::network::PipeNetwork;
use byglab_core::pipe::Pipe;
use byglab_core::solver::run_to_time;
use byglab_core::WallProperties;

#[test]
fn quiescent_gas_relaxes_toward_wall_temperature_with_the_closed_form_time_constant() {
    let gas = GasProperties::AIR;
    let diameter = 0.05;
    let initial_pressure = 1.5e5;
    let initial_temperature_kelvin = 293.15;
    let wall_temperature_kelvin = 400.0;

    // A coarse mesh is fine here: the state is spatially uniform for the
    // whole run (no wave or gradient to resolve), and coarser cells give a
    // larger CFL timestep, so far fewer steps are needed to reach a
    // meaningful fraction of the (multi-second) thermal time constant.
    let mesh = Mesh::uniform(1.0, diameter, 0.1);
    let initial_state =
        PrimitiveState::from_pressure_temperature(initial_pressure, initial_temperature_kelvin, 0.0, &gas);
    let wall = WallProperties { roughness: 0.0, wall_temperature_kelvin: Some(wall_temperature_kelvin) };
    let pipe =
        Pipe::uniform_initial_state(mesh, initial_state, &gas, BoundaryCondition::ClosedEnd, BoundaryCondition::ClosedEnd)
            .with_wall(wall);

    let mut network = PipeNetwork::single_pipe(pipe);

    let specific_heat_cv = gas.gas_constant / (gas.gamma - 1.0);
    let specific_heat_cp = gas.gamma * gas.gas_constant / (gas.gamma - 1.0);
    let thermal_conductivity = gas.dynamic_viscosity * specific_heat_cp / gas.prandtl_number;
    let nusselt_laminar_floor = 3.66;
    let convective_coefficient = nusselt_laminar_floor * thermal_conductivity / diameter;
    let density = initial_pressure / (gas.gas_constant * initial_temperature_kelvin);
    let tau = density * specific_heat_cv * diameter / (4.0 * convective_coefficient);

    let run_duration = 0.5 * tau;
    let elapsed = run_to_time(&mut network, &gas, 0.5, run_duration);

    let expected_temperature =
        wall_temperature_kelvin - (wall_temperature_kelvin - initial_temperature_kelvin) * (-elapsed / tau).exp();

    let mut max_temperature_deviation = 0.0_f64;
    let mut max_speed = 0.0_f64;
    for state in &network.pipes[0].state {
        let primitive = state.to_primitive(&gas);
        let temperature = primitive.temperature_kelvin(&gas);
        max_temperature_deviation = max_temperature_deviation.max((temperature - expected_temperature).abs());
        max_speed = max_speed.max(primitive.velocity.abs());
    }

    println!(
        "tau = {tau:.3} s, ran to {elapsed:.3} s, expected T = {expected_temperature:.3} K, max deviation = {max_temperature_deviation:.5} K"
    );

    assert!(max_speed < 1e-6, "heat transfer alone should not generate velocity, got {max_speed} m/s");
    // Forward-Euler integration of a stable relaxation ODE with
    // dt/tau ~ 1e-5 here - expect the accumulated integration error to be
    // tiny relative to the ~100 K total temperature rise; see the doc
    // comment printed above for the actual measured value.
    assert!(
        max_temperature_deviation < 0.5,
        "temperature deviated from the closed-form relaxation by {max_temperature_deviation} K"
    );
}
