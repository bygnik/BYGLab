//! The classic Shu-Osher shock/entropy-wave interaction problem (Shu &
//! Osher, "Efficient Implementation of Essentially Non-oscillatory
//! Shock-Capturing Schemes, II", 1989) — a standard benchmark for whether a
//! shock-capturing scheme spuriously couples a smooth density (entropy)
//! wave into the pressure field.
//!
//! Standard non-dimensional problem (`gamma=1.4`, domain `x in [-5, 5]`):
//! a Mach-3 shock, initially at `x=-4`, advances rightward into a region
//! at rest with sinusoidally-varying density but *uniform* pressure. Most
//! published treatments compare the post-shock density profile at `t=1.8`
//! against a very fine reference-grid solution (there's no closed-form
//! exact solution once the shock has interacted with the density wave) —
//! that comparison is not attempted here.
//!
//! What *is* checked, and does have an exact, checkable answer: pressure
//! in the region the shock has not yet reached must stay at its initial
//! uniform value (1.0), to floating-point precision. This is a stronger
//! property than it looks. Ahead of the shock, velocity is uniformly zero
//! and pressure is uniformly 1.0, while density varies sinusoidally —
//! since this solver's MUSCL reconstruction (`reconstruction.rs`) limits
//! each primitive variable (density, velocity, pressure) *independently*,
//! a spatially-varying density has no way to leak into the reconstructed
//! pressure or velocity at any face: both stay exactly at their uniform
//! values regardless of the density gradient, so momentum flux
//! (`rho*u^2 + p`) is exactly `1.0` at every face ahead of the shock, mass
//! flux (`rho*u`) is exactly `0.0` there regardless of density, and energy
//! flux is exactly `0.0` too — every face flux ahead of the shock is
//! exactly constant, so nothing in that region should change at all until
//! the shock's actual domain of dependence reaches it. A scheme that
//! cross-couples variables during reconstruction or flux computation (a
//! real failure mode - see the module's opening paragraph) would fail
//! this cleanly.

use byglab_core::boundary::BoundaryCondition;
use byglab_core::gas::{ConservedState, GasProperties, PrimitiveState};
use byglab_core::mesh::Mesh;
use byglab_core::network::PipeNetwork;
use byglab_core::pipe::Pipe;
use byglab_core::solver::run_to_time;

const DOMAIN_LENGTH: f64 = 10.0;
const DOMAIN_LEFT_X: f64 = -5.0;
const SHOCK_POSITION_X: f64 = -4.0;
const RUN_DURATION: f64 = 1.8;

#[test]
fn pressure_ahead_of_the_shock_stays_exactly_at_its_initial_uniform_value() {
    let gas = GasProperties::AIR; // gamma=1.4 matches the standard problem; gas_constant is unused (states are built directly, not via temperature)

    let mesh = Mesh::uniform(DOMAIN_LENGTH, 1.0, DOMAIN_LENGTH / 400.0);
    let state: Vec<ConservedState> = mesh
        .cells
        .iter()
        .map(|cell| {
            let x = DOMAIN_LEFT_X + cell.center;
            let primitive = if x < SHOCK_POSITION_X {
                PrimitiveState { density: 3.857143, velocity: 2.629369, pressure: 10.333333 }
            } else {
                PrimitiveState { density: 1.0 + 0.2 * (5.0 * x).sin(), velocity: 0.0, pressure: 1.0 }
            };
            primitive.to_conserved(&gas)
        })
        .collect();

    let pipe = Pipe {
        mesh,
        state,
        left_boundary: BoundaryCondition::Outflow,
        right_boundary: BoundaryCondition::Outflow,
        wall: None,
    };
    let mut network = PipeNetwork::single_pipe(pipe);
    let elapsed = run_to_time(&mut network, &gas, 0.5, RUN_DURATION);

    // Shock speed is close to Mach 3 relative to the undisturbed ahead gas
    // (sound speed ~1.18 there), so by t=1.8 the shock/post-shock
    // oscillatory region has advanced to roughly x=-4 + 3.5*1.8 =~ 2.3.
    // x > 3.5 gives comfortable margin against that estimate while still
    // comparing well inside the domain (right edge at x=5).
    let unshocked_region_start_x = 3.5;

    let pipe = &network.pipes[0];
    let mut max_pressure_deviation = 0.0_f64;
    let mut checked_cells = 0;
    for (cell, state) in pipe.mesh.cells.iter().zip(pipe.state.iter()) {
        let x = DOMAIN_LEFT_X + cell.center;
        if x < unshocked_region_start_x {
            continue;
        }
        checked_cells += 1;
        let primitive = state.to_primitive(&gas);
        max_pressure_deviation = max_pressure_deviation.max((primitive.pressure - 1.0).abs());
    }

    println!(
        "ran to t={elapsed:.3}; checked {checked_cells} cells with x > {unshocked_region_start_x}; max pressure deviation from the exact value of 1.0: {max_pressure_deviation:e}"
    );

    assert!(checked_cells > 0, "the unshocked check region was empty - adjust unshocked_region_start_x");
    assert!(
        max_pressure_deviation < 1e-9,
        "pressure ahead of the shock deviated from its exact initial value of 1.0 by {max_pressure_deviation:e} - density variation is leaking into pressure"
    );
}
