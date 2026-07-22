//! Total mass conservation in a closed, *tapered* pipe under real (moving,
//! non-trivial) flow.
//!
//! `solver::tests::total_mass_is_conserved_in_a_closed_pipe_with_a_traveling_pressure_pulse`
//! (in `src/solver.rs`) already checks this for a constant-diameter pipe,
//! but `apply_face_fluxes`'s area-weighted flux differencing (added for
//! taper support — see `pipe.rs`) is a genuinely different code path, and
//! deserves its own check rather than assuming the constant-area result
//! carries over. The mass equation has no source term at all (only
//! momentum gets the geometric taper source, and momentum/energy get the
//! optional wall sources) — mass changes only via flux differencing, and a
//! numerical flux applied with opposite sign to both neighboring cells
//! conserves exactly regardless of how the flux itself was computed
//! (Lax-Wendroff conservation theorem), so exact conservation is expected
//! here just as it was for the constant-area case, modulo the same
//! floating-point rounding from repeated `PrimitiveState`<->`ConservedState`
//! conversions already documented on the constant-area test.

use byglab_core::boundary::BoundaryCondition;
use byglab_core::gas::{GasProperties, PrimitiveState};
use byglab_core::mesh::Mesh;
use byglab_core::network::PipeNetwork;
use byglab_core::pipe::Pipe;
use byglab_core::solver::run_to_time;

/// Total mass in the pipe, properly area-weighted per cell (`rho * area *
/// width`) - not just `rho * width`, since a tapered pipe's cells don't
/// all have the same cross-sectional area.
fn total_mass(pipe: &Pipe, gas: &GasProperties) -> f64 {
    pipe.state
        .iter()
        .zip(pipe.mesh.cells.iter())
        .map(|(state, cell)| state.to_primitive(gas).density * cell.area * cell.width)
        .sum()
}

#[test]
fn total_mass_is_conserved_in_a_closed_tapered_pipe_with_a_traveling_pressure_pulse() {
    let gas = GasProperties::AIR;
    let mesh = Mesh::tapered(1.0, 0.05, 0.10, 0.01);
    let uniform_state = PrimitiveState::from_pressure_temperature(150_000.0, 293.15, 0.0, &gas);
    let mut pipe =
        Pipe::uniform_initial_state(mesh, uniform_state, &gas, BoundaryCondition::ClosedEnd, BoundaryCondition::ClosedEnd);

    // Perturb a handful of cells near the middle to create a non-trivial,
    // moving disturbance instead of a trivial at-rest state (same
    // technique as the constant-area test in `src/solver.rs`).
    for cell_state in pipe.state.iter_mut().skip(45).take(10) {
        let mut primitive = cell_state.to_primitive(&gas);
        primitive.pressure *= 1.5;
        *cell_state = primitive.to_conserved(&gas);
    }

    let initial_mass = total_mass(&pipe, &gas);

    let mut network = PipeNetwork::single_pipe(pipe);
    run_to_time(&mut network, &gas, 0.5, 0.005);

    let final_mass = total_mass(&network.pipes[0], &gas);

    let relative_error = (final_mass - initial_mass).abs() / initial_mass;
    println!("closed pipe, tapered (0.05m to 0.10m): mass conservation relative error {relative_error:e} ({:.6}%)", relative_error * 100.0);
    // Measured 0.000145% - two orders of magnitude inside the 0.1% target,
    // matching the same order of floating-point rounding measured for the
    // constant-area case (see `solver.rs`). 1e-4 (0.01%) leaves real
    // margin below the measured value while still comfortably beating the
    // 0.1% target.
    assert!(relative_error < 1e-4, "mass not conserved in a tapered pipe: relative error {relative_error:e} exceeds the 0.1% target");
}
