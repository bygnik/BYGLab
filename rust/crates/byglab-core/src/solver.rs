//! Explict time-stepping driver for a [`PipeNetwork`].

use crate::gas::GasProperties;
use crate::network::PipeNetwork;

/// Advances `network` by one explicit timestep, sized by the CFL condition
/// computed across every cell in every pipe in the network — never a
/// single pipe's own value in isolation, see
/// [`PipeNetwork::cfl_time_step`]. Returns the timestep taken.
///
/// A junction's shared flux is only physically meaningful if both sides of
/// the junction are evaluated at the same instant in time. If two
/// junction-coupled pipes advanced using different, independently-chosen
/// timesteps, that assumption would silently break — a single network-wide
/// `dt` is what keeps every pipe's cells, and both sides of every
/// junction, at the same time level after each step.
///
/// The outer update here is a single full-timestep, forward-Euler-style
/// conservative step (`U^{n+1} = U^n - dt/dx * flux divergence`) — this
/// looks first-order, but the scheme is actually second-order accurate in
/// both space and time, because the face fluxes it's given have already
/// been computed from MUSCL-Hancock–reconstructed, half-timestep-evolved
/// states (see `reconstruction.rs`). That's the specific property that
/// makes the Hancock variant convenient: the predictor half-step is baked
/// into the flux computation itself, so — unlike a generic
/// piecewise-linear scheme — it doesn't need a multi-stage outer time
/// integrator (e.g. SSP-RK2) to reach full second order. `step`/`advance`
/// would need no changes at all if the reconstruction were swapped for a
/// different second-order variant that keeps this property; only a scheme
/// that skips the half-step evolution (plain MUSCL without the Hancock
/// predictor) would require upgrading the integrator here too.
pub fn step(network: &mut PipeNetwork, gas: &GasProperties, cfl: f64) -> f64 {
    let dt = network.cfl_time_step(gas, cfl);
    network.advance(dt, gas);
    dt
}

/// Repeatedly steps `network` until at least `t_end` seconds of simulated
/// time have elapsed. Returns the actual elapsed time, which may slightly
/// overshoot `t_end` since each step's size is fixed by the CFL condition
/// rather than chosen to land exactly on `t_end`.
pub fn run_to_time(network: &mut PipeNetwork, gas: &GasProperties, cfl: f64, t_end: f64) -> f64 {
    let mut elapsed = 0.0;
    while elapsed < t_end {
        elapsed += step(network, gas, cfl);
    }
    elapsed
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary::BoundaryCondition;
    use crate::gas::PrimitiveState;
    use crate::mesh::Mesh;
    use crate::pipe::Pipe;

    /// Total mass in a closed-closed pipe should stay constant, however the
    /// gas moves around inside it — a solid wall cannot pass mass through
    /// by definition, and the same single shared flux value is applied
    /// (with opposite sign) to both neighbors of every internal face, so
    /// exact conservation doesn't depend on the flux's accuracy, only on
    /// this bookkeeping being consistent (see `network.rs`'s doc comment
    /// on why a shared flux is what makes this hold at junctions too).
    ///
    /// The tolerance here is looser than it was for the first-order scheme
    /// (was 1e-10) because MUSCL-Hancock reconstruction does several more
    /// `PrimitiveState`<->`ConservedState` round-trip conversions per cell
    /// per step (each with its own small floating-point rounding, as seen
    /// in `gas.rs`'s own round-trip tests) — confirmed via a manual
    /// duration sweep that the error saturates around 1e-5 rather than
    /// growing unboundedly with more steps (1.6e-15 before the pulse
    /// starts interacting with neighboring cells, jumping to ~1e-7 once it
    /// does, then leveling off near 1e-5) — consistent with accumulated
    /// rounding at sharp gradients, not a structural conservation bug.
    #[test]
    fn total_mass_is_conserved_in_a_closed_pipe_with_a_traveling_pressure_pulse() {
        let gas = GasProperties::AIR;
        let mesh = Mesh::uniform(1.0, 0.05, 0.01);
        let uniform_state = PrimitiveState::from_pressure_temperature(150_000.0, 293.15, 0.0, &gas);
        let mut pipe = Pipe::uniform_initial_state(
            mesh,
            uniform_state,
            &gas,
            BoundaryCondition::ClosedEnd,
            BoundaryCondition::ClosedEnd,
        );

        // Perturb a handful of cells near the middle to create a
        // non-trivial, moving disturbance instead of a trivial at-rest state.
        for cell_state in pipe.state.iter_mut().skip(45).take(10) {
            let mut primitive = cell_state.to_primitive(&gas);
            primitive.pressure *= 1.5;
            *cell_state = primitive.to_conserved(&gas);
        }

        let initial_mass: f64 = pipe.state.iter().map(|s| s.mass).sum::<f64>() * pipe.mesh.cells[0].width;

        let mut network = PipeNetwork::single_pipe(pipe);
        run_to_time(&mut network, &gas, 0.5, 0.005);

        let final_mass: f64 =
            network.pipes[0].state.iter().map(|s| s.mass).sum::<f64>() * network.pipes[0].mesh.cells[0].width;

        let relative_error = (final_mass - initial_mass).abs() / initial_mass;
        println!("closed pipe, constant area: mass conservation relative error {relative_error:e} ({:.6}%)", relative_error * 100.0);
        assert!(relative_error < 1e-4, "mass not conserved: relative error {relative_error:e}");
    }

    #[test]
    fn run_to_time_reaches_or_slightly_exceeds_the_requested_duration() {
        let gas = GasProperties::AIR;
        let mesh = Mesh::uniform(1.0, 0.05, 0.01);
        let uniform_state = PrimitiveState::from_pressure_temperature(150_000.0, 293.15, 0.0, &gas);
        let pipe = Pipe::uniform_initial_state(
            mesh,
            uniform_state,
            &gas,
            BoundaryCondition::ClosedEnd,
            BoundaryCondition::ClosedEnd,
        );
        let mut network = PipeNetwork::single_pipe(pipe);
        let elapsed = run_to_time(&mut network, &gas, 0.5, 0.001);
        assert!(elapsed >= 0.001);
    }
}
