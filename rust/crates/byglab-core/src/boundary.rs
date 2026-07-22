//! How a pipe's end is terminated.
//!
//! A junction-coupled end (two pipes joined together, see
//! [`crate::network::Junction`]) is deliberately NOT a variant here — a
//! junction supplies an externally-computed face flux instead, bypassing
//! `BoundaryCondition` for that end entirely (see [`crate::pipe::Pipe`]'s
//! doc comment for how that seam works).

use crate::gas::{GasProperties, PrimitiveState};
use serde::{Deserialize, Serialize};

/// The three ways (besides a junction) a pipe's end can be terminated.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum BoundaryCondition {
    /// A solid, reflective wall: flow cannot cross it, so the outward
    /// velocity at the wall is zero.
    ClosedEnd,
    /// A large reservoir held at a fixed pressure and temperature — e.g.
    /// atmosphere at the open end of an intake or exhaust runner.
    ///
    /// This is a simplification: the ghost cell is simply the reservoir's
    /// own state, regardless of the interior flow direction or speed. A
    /// physically complete open-end model (matching OpenWAM's own
    /// quasi-steady discharge-coefficient treatment) would distinguish
    /// subsonic/supersonic and inflow/outflow cases — deliberately not
    /// attempted here. Not validated against an exact solution for this
    /// reason (see `benchmarks/openwam/cases/single_pipe/` in the OpenWAM
    /// reference suite, which carries the same caveat).
    Reservoir { pressure: f64, temperature_kelvin: f64 },
    /// Zero-gradient (transmissive/non-reflecting) outflow: the ghost cell
    /// is simply a copy of the adjacent interior cell, so a wave crossing
    /// this boundary sees no impedance change and generates no reflection.
    /// Used for numerical test domains meant to be "large enough that the
    /// boundary doesn't matter" (e.g. `tests/shu_osher.rs`) rather than a
    /// physical pipe termination.
    Outflow,
}

impl BoundaryCondition {
    /// Computes the ghost-cell state to pair with `interior` (the
    /// adjacent real cell's state) when solving the Riemann problem at
    /// this boundary.
    pub fn ghost_state(&self, interior: PrimitiveState, gas: &GasProperties) -> PrimitiveState {
        match self {
            BoundaryCondition::ClosedEnd => PrimitiveState {
                density: interior.density,
                velocity: -interior.velocity,
                pressure: interior.pressure,
            },
            BoundaryCondition::Reservoir { pressure, temperature_kelvin } => {
                PrimitiveState::from_pressure_temperature(*pressure, *temperature_kelvin, 0.0, gas)
            }
            BoundaryCondition::Outflow => interior,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closed_end_mirrors_velocity_and_keeps_density_and_pressure() {
        let gas = GasProperties::AIR;
        let interior = PrimitiveState { density: 1.2, velocity: 30.0, pressure: 150_000.0 };
        let ghost = BoundaryCondition::ClosedEnd.ghost_state(interior, &gas);
        assert_eq!(ghost.velocity, -30.0);
        assert_eq!(ghost.density, interior.density);
        assert_eq!(ghost.pressure, interior.pressure);
    }

    #[test]
    fn reservoir_ghost_state_ignores_interior_and_uses_reservoir_conditions() {
        let gas = GasProperties::AIR;
        let interior = PrimitiveState { density: 5.0, velocity: 200.0, pressure: 500_000.0 };
        let boundary = BoundaryCondition::Reservoir { pressure: 101_325.0, temperature_kelvin: 293.15 };
        let ghost = boundary.ghost_state(interior, &gas);
        assert_eq!(ghost.pressure, 101_325.0);
        assert_eq!(ghost.velocity, 0.0);
    }

    #[test]
    fn outflow_ghost_state_exactly_copies_the_interior_state() {
        let gas = GasProperties::AIR;
        let interior = PrimitiveState { density: 3.857143, velocity: 2.629369, pressure: 10.333333 };
        let ghost = BoundaryCondition::Outflow.ghost_state(interior, &gas);
        assert_eq!(ghost, interior);
    }
}
