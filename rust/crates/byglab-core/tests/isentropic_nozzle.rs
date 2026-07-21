//! Steady subsonic flow through a converging (tapered) duct, compared
//! against the exact isentropic area-Mach relation
//! (`tests/support/isentropic_nozzle.rs`) — validates the taper source
//! term under genuine flow, complementing the at-rest well-balanced check
//! in `tests/tapered_pipe_stays_at_rest.rs`.
//!
//! Driven by `Reservoir` boundary conditions at both ends of a converging
//! duct, run until the pressure profile stops changing between two check
//! points in time (rather than an arbitrary fixed duration).
//!
//! `BoundaryCondition::Reservoir` is a simplified, non-characteristic
//! far-field state (see its own doc comment) - not a real stagnation-plenum
//! or static-pressure-outlet condition - so the *absolute* stagnation
//! conditions it settles to are not expected to exactly match the nominal
//! `stagnation_pressure`/`exit_static_pressure` values used to drive it
//! (measured deviation there was large: tens of percent). Comparing
//! against those nominal values would conflate that already-documented
//! boundary-condition approximation with the thing this test actually
//! wants to check: whether the *interior* flow correctly follows quasi-1D
//! isentropic theory as area varies, which is what the taper source term
//! is responsible for. So instead: calibrate the exact solution from one
//! interior reference station's own simulated state (trivially exact
//! *there*, by construction), then check that every *other* interior
//! station's pressure follows the same isentropic curve - this isolates
//! taper-term correctness from boundary-condition fidelity.

mod support;
use byglab_core::boundary::BoundaryCondition;
use byglab_core::gas::GasProperties;
use byglab_core::mesh::Mesh;
use byglab_core::network::PipeNetwork;
use byglab_core::pipe::Pipe;
use byglab_core::solver::run_to_time;
use support::isentropic_nozzle::ExactNozzleFlow;

#[test]
fn steady_flow_through_a_converging_duct_matches_the_isentropic_area_mach_relation() {
    let gas = GasProperties::AIR;

    let stagnation_pressure = 1.5e5;
    let stagnation_temperature_kelvin = 293.15;
    let exit_static_pressure = 1.3e5;
    let diameter_left = 0.08;
    let diameter_right = 0.05;

    let mesh = Mesh::tapered(1.0, diameter_left, diameter_right, 0.01);

    // Initial condition is an arbitrary at-rest guess between the two
    // driving pressures - the point of running to a self-consistency check
    // below is to not depend on this being close to the eventual steady
    // state.
    let initial_state = byglab_core::gas::PrimitiveState::from_pressure_temperature(
        0.5 * (stagnation_pressure + exit_static_pressure),
        stagnation_temperature_kelvin,
        0.0,
        &gas,
    );
    let pipe = Pipe::uniform_initial_state(
        mesh,
        initial_state,
        &gas,
        BoundaryCondition::Reservoir { pressure: stagnation_pressure, temperature_kelvin: stagnation_temperature_kelvin },
        BoundaryCondition::Reservoir { pressure: exit_static_pressure, temperature_kelvin: stagnation_temperature_kelvin },
    );
    let mut network = PipeNetwork::single_pipe(pipe);

    // Run until the pressure profile stops changing between successive
    // checks, rather than assuming an arbitrary fixed duration reaches
    // steady state.
    let check_interval = 0.005;
    let mut elapsed = run_to_time(&mut network, &gas, 0.5, 0.02);
    let mut previous_pressures: Vec<f64> = pressures(&network, &gas);
    let mut converged = false;
    for _ in 0..20 {
        // `run_to_time` runs for approximately `check_interval` more
        // seconds of new simulated time each call (its own `elapsed`
        // return value is local to that call, not a target on an
        // absolute clock) - accumulate the running total ourselves.
        elapsed += run_to_time(&mut network, &gas, 0.5, check_interval);
        let current_pressures = pressures(&network, &gas);
        let max_relative_change = previous_pressures
            .iter()
            .zip(&current_pressures)
            .map(|(prev, curr)| (curr - prev).abs() / prev)
            .fold(0.0_f64, f64::max);
        previous_pressures = current_pressures;
        if max_relative_change < 1e-6 {
            converged = true;
            break;
        }
    }
    assert!(converged, "pressure profile did not settle to a steady state within the iteration budget");

    // Skip the first/last few cells, closest to the approximate boundary
    // conditions, where the largest deviation from ideal quasi-1D
    // isentropic behavior is expected (see the module doc comment).
    let margin = 10;
    let pipe = &network.pipes[0];

    // Calibrate the exact solution from one interior reference station's
    // own simulated state - trivially exact *there* by construction (its
    // local stagnation pressure/temperature, derived from its own Mach
    // number, are used as the exact solution's nominal stagnation
    // conditions, and its own area/pressure as the "exit" station).
    let reference_index = margin;
    let reference_cell = pipe.mesh.cells[reference_index];
    let reference_primitive = pipe.state[reference_index].to_primitive(&gas);
    let reference_mach = reference_primitive.velocity / reference_primitive.sound_speed(&gas);
    let isentropic_factor = (1.0 + 0.5 * (gas.gamma - 1.0) * reference_mach * reference_mach).powf(gas.gamma / (gas.gamma - 1.0));
    let local_stagnation_pressure = reference_primitive.pressure * isentropic_factor;
    let local_stagnation_temperature_kelvin =
        reference_primitive.temperature_kelvin(&gas) * (1.0 + 0.5 * (gas.gamma - 1.0) * reference_mach * reference_mach);

    let exact = ExactNozzleFlow::new(
        local_stagnation_pressure,
        local_stagnation_temperature_kelvin,
        reference_cell.area,
        reference_primitive.pressure,
        gas.gamma,
    );

    let mut max_relative_pressure_error = 0.0_f64;
    let mut max_relative_pressure_error_x = 0.0_f64;
    for i in margin..pipe.cell_count() - margin {
        let cell = pipe.mesh.cells[i];
        let primitive = pipe.state[i].to_primitive(&gas);
        let (_, exact_pressure, _) = exact.state_at_area(cell.area);
        let relative_error = (primitive.pressure - exact_pressure).abs() / (stagnation_pressure - exit_static_pressure);
        if relative_error > max_relative_pressure_error {
            max_relative_pressure_error = relative_error;
            max_relative_pressure_error_x = cell.center;
        }
    }

    println!(
        "settled after t={elapsed:.4} s; reference Mach {reference_mach:.4} at x={:.3} m; max interior pressure error {:.2}% of driving range at x={max_relative_pressure_error_x:.3} m",
        reference_cell.center,
        max_relative_pressure_error * 100.0
    );

    assert!(
        max_relative_pressure_error < 0.02,
        "max interior pressure error {:.2}% of driving range exceeds 2%",
        max_relative_pressure_error * 100.0
    );
}

fn pressures(network: &PipeNetwork, gas: &GasProperties) -> Vec<f64> {
    network.pipes[0].state.iter().map(|state| state.to_primitive(gas).pressure).collect()
}
