//! End-to-end validation for `valve_port.rs`: a cylinder breathing against
//! a real, wave-propagating 1D pipe through a `ValvePort`, instead of
//! `cylinder::BreathingParameters`'s fixed reservoir.
//!
//! `valve_port.rs`'s own unit tests already validate the mdot->Flux
//! conversion formula in isolation (matching the direct analytic
//! `valve::mass_flow_rate` call, the Left/Right sign convention, and the
//! at-rest/moving-interior `ClosedEnd` identity). These tests instead
//! validate the *coupled driver* (`step_pipe_cylinder`/
//! `run_pipe_cylinder_to_time`) against real ground truth.

use byglab_core::boundary::BoundaryCondition;
use byglab_core::camshaft::{self, CamProfile};
use byglab_core::crank_mechanism::CrankMechanism;
use byglab_core::cylinder::{BreathingParameters, Cylinder, CylinderState};
use byglab_core::gas::{GasProperties, PrimitiveState};
use byglab_core::mesh::Mesh;
use byglab_core::network::{PipeEnd, PipeEndRef, PipeNetwork};
use byglab_core::pipe::Pipe;
use byglab_core::valve::{self, ValveGeometry};
use byglab_core::valve_port::{run_pipe_cylinder_to_time, step_pipe_cylinder, ValvePort};

/// S54B32 reference geometry (bore 87mm, stroke 91mm, rod 139mm, CR
/// 11.5:1), the same real numbers used throughout `cylinder.rs`'s own
/// tests (see its `s54b32_cylinder` helper).
fn s54b32_cylinder() -> Cylinder {
    let bore = 0.087;
    let stroke = 0.091;
    let compression_ratio = 11.5;
    let piston_area = std::f64::consts::PI / 4.0 * bore * bore;
    let displaced_volume = piston_area * stroke;
    let clearance_volume = displaced_volume / (compression_ratio - 1.0);
    Cylinder { crank_mechanism: CrankMechanism::new(stroke, 0.139, 0.0), bore, clearance_volume }
}

fn total_pipe_mass(pipe: &Pipe) -> f64 {
    pipe.state.iter().zip(&pipe.mesh.cells).map(|(s, cell)| s.mass * cell.area * cell.width).sum()
}

fn total_pipe_energy(pipe: &Pipe) -> f64 {
    pipe.state.iter().zip(&pipe.mesh.cells).map(|(s, cell)| s.energy * cell.area * cell.width).sum()
}

#[test]
fn total_mass_and_energy_are_conserved_across_a_short_choked_flow_window() {
    // A *rigid* cylinder (zero crank radius, same technique
    // `cylinder::tests::breathing_choked_mass_flow_produces_exactly_linear_growth_in_a_rigid_cylinder`
    // uses) over a tiny crank-angle window, reusing that same test's exact
    // pressures/temperatures. Two deliberate choices isolate exactly the
    // property this coupling exists to deliver, without confounding it
    // with unrelated physics:
    // - Rigid cylinder => zero piston motion => zero motoring work, so
    //   "cylinder + pipe" energy is expected to conserve exactly; with a
    //   real moving piston, piston work is a genuine external energy
    //   source/sink that has nothing to do with valve/pipe conservation
    //   (confirmed - an earlier version of this test used the real S54B32
    //   geometry and measured a ~15% energy "error" that was entirely this
    //   motoring-work term, not a conservation bug).
    // - A tiny window (matching the original choked-flow test's) means the
    //   wave launched into the pipe by the valve never reaches the far
    //   closed end within the test - avoiding wave reflections, which
    //   (even though they conserve exactly in the exact HLLC algebra for a
    //   perfectly mirrored boundary state) pick up a small extra
    //   MUSCL-reconstruction-driven asymmetry once real gradients build up
    //   at a boundary cell after reflecting.
    let gas = GasProperties::AIR;
    let clearance_volume = 5.0e-5;
    let cylinder = Cylinder { crank_mechanism: CrankMechanism::new(0.0, 0.139, 0.0), bore: 0.087, clearance_volume };

    let cam = CamProfile { max_lift: 0.010, opening_angle_radians: -1000.0, duration_radians: 2000.0 };
    let valve_geometry = ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 45.0_f64.to_radians() };
    let port = ValvePort {
        pipe_end: PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
        cam,
        valve: valve_geometry,
        discharge_coefficient: 0.7,
    };
    let angular_velocity = 500.0;

    let pipe_pressure = 10.0e5;
    let pipe_temperature_kelvin = 400.0;
    let mesh = Mesh::uniform(0.3, 0.04, 0.01);
    let pipe = Pipe::uniform_initial_state(
        mesh,
        PrimitiveState::from_pressure_temperature(pipe_pressure, pipe_temperature_kelvin, 0.0, &gas),
        &gas,
        BoundaryCondition::ClosedEnd,
        BoundaryCondition::ClosedEnd, // overridden every step by the ValvePort above
    );
    let mut network = PipeNetwork::single_pipe(pipe);

    let theta_start = -0.002_f64;
    let theta_end = 0.002_f64;
    let initial_cylinder_state =
        CylinderState::from_pressure_temperature(1.0e5, 350.0, cylinder.volume(theta_start), &gas);

    let initial_total_mass = total_pipe_mass(&network.pipes[0]) + initial_cylinder_state.mass;
    let initial_total_energy = total_pipe_energy(&network.pipes[0]) + initial_cylinder_state.internal_energy;

    let cfl = 0.5;
    let final_cylinder_state = run_pipe_cylinder_to_time(
        &mut network,
        &cylinder,
        initial_cylinder_state,
        theta_start,
        theta_end,
        &port,
        angular_velocity,
        &gas,
        cfl,
    );

    let final_total_mass = total_pipe_mass(&network.pipes[0]) + final_cylinder_state.mass;
    let final_total_energy = total_pipe_energy(&network.pipes[0]) + final_cylinder_state.internal_energy;

    let mass_relative_error = (final_total_mass - initial_total_mass).abs() / initial_total_mass;
    let energy_relative_error = (final_total_energy - initial_total_energy).abs() / initial_total_energy;
    println!(
        "pipe+cylinder conservation: mass relative error {mass_relative_error:e} ({:.8}%), energy relative error {energy_relative_error:e} ({:.8}%)",
        mass_relative_error * 100.0,
        energy_relative_error * 100.0
    );

    // Matches this crate's existing precedent tolerance for MUSCL
    // round-trip/floating-point rounding noise in a closed system (see
    // `tests/mass_conservation_taper.rs`), not a "roughly conserved" bound.
    assert!(mass_relative_error < 1e-4, "mass not conserved: relative error {mass_relative_error:e}");
    assert!(energy_relative_error < 1e-4, "energy not conserved: relative error {energy_relative_error:e}");
}

#[test]
fn coupled_path_agrees_with_the_fixed_reservoir_path_over_a_short_window() {
    // A softer, corroborating sanity check (explicitly not claimed exact -
    // these are two different integrators: fixed-reservoir RK4 vs. the
    // coupled path's pipe-CFL-substepped scheme). Over a short-enough
    // window that the pipe's boundary-cell state hasn't meaningfully
    // evolved from its initial uniform value, the coupled path should
    // track the already-validated fixed-reservoir `integrate_breathing`
    // path closely.
    let gas = GasProperties::AIR;
    let cylinder = s54b32_cylinder();
    let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();

    let cam = CamProfile { max_lift: 0.009, opening_angle_radians: bdc_angle - 1.2, duration_radians: 2.4 };
    let valve_geometry = ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 45.0_f64.to_radians() };
    let discharge_coefficient = 0.7;
    let angular_velocity = 500.0;
    let reservoir_pressure = 1.0e5;
    let reservoir_temperature_kelvin = 300.0;

    // Narrow enough that a wave launched into the pipe by the valve
    // travels only a small fraction of the pipe's length within this
    // window (at omega=500 rad/s this is ~4e-5s of real time, versus
    // ~3.4e-4s for sound to cross the 0.3m pipe once) - the regime where
    // "the pipe's boundary-cell state hasn't meaningfully evolved from its
    // initial uniform value" is actually true, not just assumed.
    let theta_start = bdc_angle - 0.02;
    let theta_end = bdc_angle + 0.02;
    let initial_cylinder_state = CylinderState::from_pressure_temperature(0.5e5, 350.0, cylinder.volume(theta_start), &gas);

    let breathing_params = BreathingParameters {
        cam,
        valve: valve_geometry,
        discharge_coefficient,
        reservoir_pressure,
        reservoir_temperature_kelvin,
        angular_velocity_radians_per_second: angular_velocity,
    };
    let reservoir_final_state =
        cylinder.integrate_breathing(&gas, initial_cylinder_state, theta_start, theta_end, 400, &breathing_params);

    let port = ValvePort {
        pipe_end: PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
        cam,
        valve: valve_geometry,
        discharge_coefficient,
    };
    let mesh = Mesh::uniform(0.3, 0.04, 0.01);
    let pipe = Pipe::uniform_initial_state(
        mesh,
        PrimitiveState::from_pressure_temperature(reservoir_pressure, reservoir_temperature_kelvin, 0.0, &gas),
        &gas,
        BoundaryCondition::ClosedEnd,
        BoundaryCondition::ClosedEnd,
    );
    let mut network = PipeNetwork::single_pipe(pipe);
    let coupled_final_state = run_pipe_cylinder_to_time(
        &mut network,
        &cylinder,
        initial_cylinder_state,
        theta_start,
        theta_end,
        &port,
        angular_velocity,
        &gas,
        0.5,
    );

    let volume_at_end = cylinder.volume(theta_end);
    let reservoir_pressure_at_end = reservoir_final_state.pressure(volume_at_end, &gas);
    let coupled_pressure_at_end = coupled_final_state.pressure(volume_at_end, &gas);
    let relative_difference = (coupled_pressure_at_end - reservoir_pressure_at_end).abs() / reservoir_pressure_at_end;
    println!(
        "coupled vs fixed-reservoir over a short window: reservoir path={reservoir_pressure_at_end:e} Pa, coupled path={coupled_pressure_at_end:e} Pa, relative difference={relative_difference:e}"
    );
    assert!(
        relative_difference < 0.05,
        "expected the coupled path to track the fixed-reservoir path within a few percent over a short window, got {relative_difference:e}"
    );
}

#[test]
fn closed_valve_coupled_trajectory_loosely_tracks_pure_motoring() {
    // Another soft sanity check: with the valve held fully closed, the
    // pipe exchanges no mass/energy with the cylinder at all, so the
    // coupled path's cylinder trajectory should be governed by motoring
    // alone (mass exactly constant, energy following `-p dV/dtheta`) - not
    // claimed to match `integrate_motoring`'s RK4 trajectory to machine
    // precision (different integrators/step counts), just loosely.
    let gas = GasProperties::AIR;
    let cylinder = s54b32_cylinder();
    let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();

    let cam = CamProfile { max_lift: 0.009, opening_angle_radians: bdc_angle - 1.2, duration_radians: 2.4 };
    let valve_geometry = ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 45.0_f64.to_radians() };
    let port = ValvePort {
        pipe_end: PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
        cam,
        valve: valve_geometry,
        discharge_coefficient: 0.0, // valve held fully closed
    };
    let angular_velocity = 500.0;

    // Window far from the cam's own open event, well away from BDC, so the
    // real cam lift would also be zero here regardless (belt-and-braces -
    // discharge_coefficient=0.0 above is what actually forces mdot=0).
    let theta_start = -1.0;
    let theta_end = 1.0;
    let initial_cylinder_state = CylinderState::from_pressure_temperature(20.0e5, 900.0, cylinder.volume(theta_start), &gas);

    let motoring_final_state = cylinder.integrate_motoring(&gas, initial_cylinder_state, theta_start, theta_end, 200);

    let mesh = Mesh::uniform(0.3, 0.04, 0.01);
    let pipe = Pipe::uniform_initial_state(
        mesh,
        PrimitiveState::from_pressure_temperature(1.0e5, 300.0, 0.0, &gas),
        &gas,
        BoundaryCondition::ClosedEnd,
        BoundaryCondition::ClosedEnd,
    );
    let mut network = PipeNetwork::single_pipe(pipe);
    let coupled_final_state = run_pipe_cylinder_to_time(
        &mut network,
        &cylinder,
        initial_cylinder_state,
        theta_start,
        theta_end,
        &port,
        angular_velocity,
        &gas,
        0.5,
    );

    assert_eq!(coupled_final_state.mass, initial_cylinder_state.mass, "a fully closed valve must not change cylinder mass at all");

    let volume_at_end = cylinder.volume(theta_end);
    let motoring_pressure = motoring_final_state.pressure(volume_at_end, &gas);
    let coupled_pressure = coupled_final_state.pressure(volume_at_end, &gas);
    let relative_difference = (coupled_pressure - motoring_pressure).abs() / motoring_pressure;
    println!("closed-valve coupled path vs. pure motoring: relative difference {relative_difference:e}");
    assert!(relative_difference < 0.01, "expected close agreement with pure motoring, got {relative_difference:e}");
}

/// Sanity check that `step_pipe_cylinder`'s single-step contract behaves as
/// documented (nonzero `dt`, a `CylinderState` with finite fields) before
/// relying on it inside the loops above.
#[test]
fn a_single_coupled_step_returns_a_positive_dt_and_a_finite_cylinder_state() {
    let gas = GasProperties::AIR;
    let cylinder = s54b32_cylinder();
    let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();
    let cam = CamProfile { max_lift: 0.009, opening_angle_radians: bdc_angle - 1.2, duration_radians: 2.4 };
    let valve_geometry = ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 45.0_f64.to_radians() };
    let port = ValvePort {
        pipe_end: PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
        cam,
        valve: valve_geometry,
        discharge_coefficient: 0.7,
    };

    let mesh = Mesh::uniform(0.3, 0.04, 0.01);
    let pipe = Pipe::uniform_initial_state(
        mesh,
        PrimitiveState::from_pressure_temperature(1.0e5, 300.0, 0.0, &gas),
        &gas,
        BoundaryCondition::ClosedEnd,
        BoundaryCondition::ClosedEnd,
    );
    let mut network = PipeNetwork::single_pipe(pipe);
    let initial_cylinder_state = CylinderState::from_pressure_temperature(0.5e5, 350.0, cylinder.volume(bdc_angle), &gas);

    let (dt, new_state) =
        step_pipe_cylinder(&mut network, &cylinder, initial_cylinder_state, bdc_angle, &port, 500.0, &gas, 0.5);

    assert!(dt > 0.0);
    assert!(new_state.mass.is_finite() && new_state.mass > 0.0);
    assert!(new_state.internal_energy.is_finite() && new_state.internal_energy > 0.0);

    // Independent sanity: the valve is wide open here (lift near its max
    // at the cam's midpoint), and the pipe (1.0e5 Pa) is higher pressure
    // than the cylinder (0.5e5 Pa) in this scenario, so mass should flow
    // INTO the cylinder.
    let lift = camshaft::lift_at(&cam, bdc_angle);
    let effective_area = 0.7 * valve::curtain_area(&valve_geometry, lift);
    assert!(effective_area > 0.0, "test setup should have a genuinely open valve at this crank angle");
    assert!(new_state.mass > initial_cylinder_state.mass, "expected net inflow from the higher-pressure pipe");
}
