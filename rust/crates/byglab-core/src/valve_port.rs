//! Binds a poppet valve ([`camshaft`]/[`valve`]) to one end of a real
//! [`crate::pipe::Pipe`], letting a [`Cylinder`] breathe against genuine
//! wave-propagating 1D flow instead of [`crate::cylinder::BreathingParameters`]'s
//! fixed reservoir.
//!
//! Deliberately kept separate from both `network.rs` (which stays
//! cylinder/valve-agnostic — Phase 1 must not learn about Phase 2 concepts)
//! and `cylinder.rs` (which stays pipe-agnostic). This module is the driver
//! that sits above both.
//!
//! # The mdot -> Flux conversion
//!
//! [`valve::mass_flow_rate`] returns a scalar `kg/s` — there are two
//! *different* areas involved, and conflating them would silently create or
//! destroy mass at every coupled step:
//!
//! - `effective_area = Cd * curtain_area(...)` — the valve's restricted
//!   throat area. Feeds [`valve::mass_flow_rate`] only.
//! - `pipe_face_area` — the pipe's own full bore cross-section at the
//!   coupled end ([`crate::mesh::Mesh::face_areas`]). Used only to convert
//!   the scalar `mdot` (kg/s, already area-integrated) into a flux density
//!   (kg/(m^2*s)) via continuity, since [`crate::pipe::Pipe::apply_face_fluxes`]
//!   always multiplies a face's flux density back by the pipe's own face
//!   area to recover a total rate.
//!
//! The face state feeding [`crate::gas::ConservedState::physical_flux`] is
//! built the same way every other flux in this solver is: the pipe's own
//! local static pressure, whichever side is upstream's temperature (mirrors
//! `cylinder.rs`'s existing `breathing_derivatives`'s upstream-temperature
//! branch verbatim, just relabeled reservoir -> pipe), and a velocity
//! derived from continuity. This means the resulting energy flux carries a
//! kinetic-energy term (`0.5*rho*u_face^3` per unit area) that the OLD
//! fixed-reservoir path's enthalpy-only formula never had (that path never
//! modeled a through-valve velocity) — this coupled path's cylinder-side
//! energy update must read its rate directly off this same [`gas::Flux`]
//! rather than reusing that older enthalpy-only formula, or exact
//! conservation would silently break by that missing term as flow rate
//! grows.
//!
//! # Sign convention
//!
//! A [`gas::Flux`] is always "quantity moving in the +x direction through
//! this face" (`riemann::hllc_flux`'s convention). So a positive flux at a
//! pipe's `Right` end is outflow; at a `Left` end, inflow. This module
//! always computes `mdot_out_of_pipe` (positive = net flow *out of* the
//! pipe, into the cylinder, regardless of which end is coupled) via
//! `valve::mass_flow_rate(pipe_side, cylinder_side, ...)`, then applies
//! `end_sign` (`Right` = `+1`, `Left` = `-1`) to get the correctly-signed
//! face velocity — verified algebraically: `flux.mass * pipe_face_area ==
//! end_sign * mdot_out_of_pipe` exactly (checked by a `debug_assert` in
//! [`valve_port_face`] and by a dedicated test below).
//!
//! # A real, honest limitation found while deriving this (not hidden)
//!
//! At zero lift (`mdot=0`), this construction gives `u_face=0` exactly, so
//! `flux = {mass: 0, momentum: pipe_pressure, energy: 0}`. It's tempting to
//! assume this always exactly matches what `riemann::hllc_flux` would give
//! against `BoundaryCondition::ClosedEnd`'s own mirrored ghost state — but
//! working through HLLC's actual approximate wave speeds (`riemann.rs`)
//! shows that's only exactly true when the pipe's boundary-cell velocity is
//! zero. For nonzero interior velocity, a real closed-end reflection adds a
//! genuine dynamic-pressure term (`rho*v*c*q`, `q` the shock/rarefaction
//! factor) that a naive "just use the local static pressure" construction
//! doesn't capture — a real wall reflection compresses (or expands) gas
//! moving toward (or away from) it. `tests` below demonstrate both the
//! at-rest identity (holds exactly) and the moving-interior divergence
//! (does not), rather than silently assuming the stronger claim.

use crate::camshaft::{self, CamProfile};
use crate::cylinder::{Cylinder, CylinderState};
use crate::gas::{Flux, GasProperties, PrimitiveState};
use crate::network::{ExternalPortFlux, PipeEnd, PipeEndRef, PipeNetwork};
use crate::valve::{self, ValveGeometry};

/// A poppet valve bound to one specific end of one specific pipe within a
/// [`PipeNetwork`]. Angular velocity is deliberately excluded (unlike
/// [`crate::cylinder::BreathingParameters`], which bundles it in) — it's a
/// crankshaft-wide operating point, not a property of one valve's geometry,
/// so it's passed separately to [`step_pipe_cylinder`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ValvePort {
    pub pipe_end: PipeEndRef,
    pub cam: CamProfile,
    pub valve: ValveGeometry,
    pub discharge_coefficient: f64,
}

fn end_sign(end: PipeEnd) -> f64 {
    match end {
        PipeEnd::Right => 1.0,
        PipeEnd::Left => -1.0,
    }
}

/// Computes `(mdot_out_of_pipe, face_state, flux)` for one `ValvePort` at
/// the given crank angle, given the pipe's current boundary-cell state (and
/// face area) and the cylinder's current pressure/temperature. See the
/// module doc comment for the full derivation. `face_state` doubles as the
/// [`ExternalPortFlux::neighbor_state`] the pipe should use for MUSCL
/// reconstruction at this end, the same way [`crate::boundary::BoundaryCondition::Reservoir`]'s
/// ghost state serves both roles.
fn valve_port_face(
    port: &ValvePort,
    crank_angle_from_tdc_radians: f64,
    pipe_boundary_state: PrimitiveState,
    pipe_face_area: f64,
    cylinder_pressure: f64,
    cylinder_temperature_kelvin: f64,
    gas: &GasProperties,
) -> (f64, PrimitiveState, Flux) {
    let lift = camshaft::lift_at(&port.cam, crank_angle_from_tdc_radians);
    let curtain_area = valve::curtain_area(&port.valve, lift);
    let effective_area = port.discharge_coefficient * curtain_area;

    let pipe_pressure = pipe_boundary_state.pressure;
    let pipe_temperature_kelvin = pipe_boundary_state.temperature_kelvin(gas);

    // Positive = net flow pipe -> cylinder, matching
    // `cylinder::breathing_derivatives`'s existing reservoir->cylinder sign
    // convention verbatim, just relabeled reservoir -> pipe.
    let mdot_out_of_pipe = valve::mass_flow_rate(
        pipe_pressure,
        pipe_temperature_kelvin,
        cylinder_pressure,
        cylinder_temperature_kelvin,
        effective_area,
        gas,
    );

    let upstream_temperature_kelvin =
        if mdot_out_of_pipe >= 0.0 { pipe_temperature_kelvin } else { cylinder_temperature_kelvin };

    let sign = end_sign(port.pipe_end.end);
    let density_face = pipe_pressure / (gas.gas_constant * upstream_temperature_kelvin);
    let velocity_face =
        if pipe_face_area > 0.0 { sign * mdot_out_of_pipe / (density_face * pipe_face_area) } else { 0.0 };

    let face_state =
        PrimitiveState::from_pressure_temperature(pipe_pressure, upstream_temperature_kelvin, velocity_face, gas);
    let flux = face_state.to_conserved(gas).physical_flux(gas);

    debug_assert!(
        (flux.mass * pipe_face_area - sign * mdot_out_of_pipe).abs() < 1e-6 * mdot_out_of_pipe.abs().max(1e-9),
        "valve-port flux/mdot invariant violated: flux.mass*area={:e}, sign*mdot={:e}",
        flux.mass * pipe_face_area,
        sign * mdot_out_of_pipe
    );

    (mdot_out_of_pipe, face_state, flux)
}

/// The pipe's boundary-cell state and face area at `end_ref`.
fn pipe_end_state_and_area(network: &PipeNetwork, end_ref: PipeEndRef, gas: &GasProperties) -> (PrimitiveState, f64) {
    let pipe = &network.pipes[end_ref.pipe_index];
    match end_ref.end {
        PipeEnd::Left => (pipe.left_boundary_cell_state(gas), pipe.mesh.face_areas[0]),
        PipeEnd::Right => {
            let face_areas = &pipe.mesh.face_areas;
            (pipe.right_boundary_cell_state(gas), face_areas[face_areas.len() - 1])
        }
    }
}

/// Advances a pipe network and a cylinder together by one shared,
/// CFL-limited timestep, exchanging mass/energy through `port`.
///
/// Evaluates `mdot`/`flux` **once**, at the pre-step pipe and cylinder
/// state, and applies that single value to both sides — mirroring how
/// [`crate::network::Junction`] guarantees exact conservation by using one
/// shared flux value on both sides of a pipe-to-pipe connection. RK4
/// sub-stepping within one network `dt` is deliberately not used here:
/// re-evaluating the flux at 3 additional perturbed states would only
/// apply one of those four evaluations to the pipe (the last-computed
/// `dt` already advanced it once), breaking that same-flux invariant. This
/// is a non-issue for accuracy in practice — CFL-limited `dt` for an
/// intake/exhaust-scale pipe is orders of magnitude finer than a valve
/// event's timescale (see the root README's roadmap notes), so a
/// breathing event is resolved by many thousands of these steps.
///
/// Returns the `dt` actually taken and the cylinder's updated state.
pub fn step_pipe_cylinder(
    network: &mut PipeNetwork,
    cylinder: &Cylinder,
    cylinder_state: CylinderState,
    crank_angle_from_tdc_radians: f64,
    port: &ValvePort,
    angular_velocity_radians_per_second: f64,
    gas: &GasProperties,
    cfl: f64,
) -> (f64, CylinderState) {
    let dt = network.cfl_time_step(gas, cfl);
    let dtheta = dt * angular_velocity_radians_per_second;

    let (pipe_boundary_state, pipe_face_area) = pipe_end_state_and_area(network, port.pipe_end, gas);
    let volume = cylinder.volume(crank_angle_from_tdc_radians);
    let cylinder_pressure = cylinder_state.pressure(volume, gas);
    let cylinder_temperature_kelvin = cylinder_state.temperature_kelvin(gas);

    let (_mdot_out_of_pipe, face_state, flux) = valve_port_face(
        port,
        crank_angle_from_tdc_radians,
        pipe_boundary_state,
        pipe_face_area,
        cylinder_pressure,
        cylinder_temperature_kelvin,
        gas,
    );

    network.advance_with_external_fluxes(
        dt,
        gas,
        &[ExternalPortFlux { end: port.pipe_end, neighbor_state: face_state, flux }],
    );

    // Read mass/energy rates directly off the SAME flux just applied to the
    // pipe (not a re-derived formula) - this is what makes the combined
    // (pipe + cylinder) system's conservation structural rather than
    // approximate. See the module doc comment for why the energy rate here
    // is not simply an enthalpy-flux term.
    let sign = end_sign(port.pipe_end.end);
    let mass_rate_out_of_pipe = sign * flux.mass * pipe_face_area;
    let energy_rate_out_of_pipe = sign * flux.energy * pipe_face_area;

    let motoring_term = cylinder.motoring_energy_derivative(gas, &cylinder_state, crank_angle_from_tdc_radians);

    let new_state = CylinderState {
        mass: cylinder_state.mass + mass_rate_out_of_pipe * dt,
        internal_energy: cylinder_state.internal_energy + motoring_term * dtheta + energy_rate_out_of_pipe * dt,
    };

    (dt, new_state)
}

/// Repeatedly calls [`step_pipe_cylinder`] until the cylinder's crank angle
/// has advanced from `crank_angle_start` to at least `crank_angle_end`
/// (`angular_velocity_radians_per_second` must be positive). Mirrors
/// [`crate::solver::run_to_time`]'s pairing with
/// [`crate::solver::step`].
pub fn run_pipe_cylinder_to_time(
    network: &mut PipeNetwork,
    cylinder: &Cylinder,
    initial_state: CylinderState,
    crank_angle_start: f64,
    crank_angle_end: f64,
    port: &ValvePort,
    angular_velocity_radians_per_second: f64,
    gas: &GasProperties,
    cfl: f64,
) -> CylinderState {
    assert!(angular_velocity_radians_per_second > 0.0, "angular velocity must be positive");
    let mut state = initial_state;
    let mut crank_angle = crank_angle_start;
    while crank_angle < crank_angle_end {
        let (dt, next_state) = step_pipe_cylinder(
            network,
            cylinder,
            state,
            crank_angle,
            port,
            angular_velocity_radians_per_second,
            gas,
            cfl,
        );
        state = next_state;
        crank_angle += dt * angular_velocity_radians_per_second;
    }
    state
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary::BoundaryCondition;
    use crate::riemann::hllc_flux;

    fn s54_cam() -> CamProfile {
        CamProfile { max_lift: 0.010, opening_angle_radians: -1000.0, duration_radians: 2000.0 }
    }

    fn s54_valve() -> ValveGeometry {
        ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 45.0_f64.to_radians() }
    }

    fn right_port() -> ValvePort {
        ValvePort {
            pipe_end: PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
            cam: s54_cam(),
            valve: s54_valve(),
            discharge_coefficient: 0.7,
        }
    }

    #[test]
    fn valve_port_face_matches_the_direct_analytic_mass_flow_rate() {
        // Reuses the exact numbers from
        // `cylinder::tests::breathing_choked_mass_flow_produces_exactly_linear_growth_in_a_rigid_cylinder`,
        // with the old "reservoir" now playing the role of the pipe's
        // boundary-cell state. A pure plumbing check: does `valve_port_face`
        // reproduce the same analytic `valve::mass_flow_rate` call it wraps.
        let gas = GasProperties::AIR;
        let port = right_port();

        let pipe_pressure = 10.0e5;
        let pipe_temperature_kelvin = 400.0;
        let pipe_boundary_state = PrimitiveState::from_pressure_temperature(pipe_pressure, pipe_temperature_kelvin, 0.0, &gas);
        let pipe_face_area = 1.0e-3;

        let cylinder_pressure = 1.0e5;
        let cylinder_temperature_kelvin = 350.0;

        let (mdot_out_of_pipe, _face_state, _flux) = valve_port_face(
            &port, 0.0, pipe_boundary_state, pipe_face_area, cylinder_pressure, cylinder_temperature_kelvin, &gas,
        );

        let lift = camshaft::lift_at(&port.cam, 0.0);
        let effective_area = port.discharge_coefficient * valve::curtain_area(&port.valve, lift);
        let expected_mdot = valve::mass_flow_rate(
            pipe_pressure, pipe_temperature_kelvin, cylinder_pressure, cylinder_temperature_kelvin, effective_area, &gas,
        );

        let relative_error = (mdot_out_of_pipe - expected_mdot).abs() / expected_mdot.abs();
        assert!(relative_error < 1e-9, "expected mdot {expected_mdot:e}, got {mdot_out_of_pipe:e} (rel. error {relative_error:e})");
    }

    #[test]
    fn valve_port_flux_mass_times_area_equals_signed_mdot() {
        let gas = GasProperties::AIR;
        let pipe_boundary_state = PrimitiveState::from_pressure_temperature(3.0e5, 380.0, 0.0, &gas);
        let pipe_face_area = 2.0e-3;

        for end in [PipeEnd::Right, PipeEnd::Left] {
            let port = ValvePort {
                pipe_end: PipeEndRef { pipe_index: 0, end },
                cam: s54_cam(),
                valve: s54_valve(),
                discharge_coefficient: 0.7,
            };
            let (mdot_out_of_pipe, _face_state, flux) =
                valve_port_face(&port, 0.0, pipe_boundary_state, pipe_face_area, 1.0e5, 320.0, &gas);
            let expected = end_sign(end) * mdot_out_of_pipe;
            let actual = flux.mass * pipe_face_area;
            let relative_error = (actual - expected).abs() / expected.abs();
            assert!(relative_error < 1e-9, "{end:?}: expected flux.mass*area={expected:e}, got {actual:e}");
        }
    }

    #[test]
    fn valve_port_flux_sign_flips_between_left_and_right_for_the_same_physical_scenario() {
        let gas = GasProperties::AIR;
        let pipe_boundary_state = PrimitiveState::from_pressure_temperature(3.0e5, 380.0, 0.0, &gas);
        let pipe_face_area = 2.0e-3;

        let right_port = ValvePort {
            pipe_end: PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
            cam: s54_cam(),
            valve: s54_valve(),
            discharge_coefficient: 0.7,
        };
        let left_port = ValvePort { pipe_end: PipeEndRef { pipe_index: 0, end: PipeEnd::Left }, ..right_port };

        let (_, _, right_flux) = valve_port_face(&right_port, 0.0, pipe_boundary_state, pipe_face_area, 1.0e5, 320.0, &gas);
        let (_, _, left_flux) = valve_port_face(&left_port, 0.0, pipe_boundary_state, pipe_face_area, 1.0e5, 320.0, &gas);

        assert!(right_flux.mass > 0.0, "pipe is upstream (outflow) - expected a positive flux at the Right end");
        assert!(left_flux.mass < 0.0, "same outflow at the Left end should be a negative (leftward) flux");
        assert!((right_flux.mass + left_flux.mass).abs() < 1e-12 * right_flux.mass.abs());
    }

    #[test]
    fn zero_lift_gives_exactly_zero_face_velocity_and_mass_flow() {
        let gas = GasProperties::AIR;
        let port = ValvePort { discharge_coefficient: 0.0, ..right_port() }; // forces effective_area = 0 regardless of lift
        let pipe_boundary_state = PrimitiveState::from_pressure_temperature(2.0e5, 350.0, 15.0, &gas);
        let (mdot_out_of_pipe, face_state, flux) =
            valve_port_face(&port, 0.0, pipe_boundary_state, 1.0e-3, 1.0e5, 300.0, &gas);
        assert_eq!(mdot_out_of_pipe, 0.0);
        assert_eq!(face_state.velocity, 0.0);
        assert_eq!(flux.mass, 0.0);
        assert_eq!(flux.energy, 0.0);
        assert!((flux.momentum - pipe_boundary_state.pressure).abs() < 1e-6);
    }

    #[test]
    fn closed_valve_flux_matches_closed_end_hllc_flux_when_the_pipe_is_at_rest() {
        // See the module doc comment: this identity holds exactly only
        // when the pipe's boundary-cell velocity is zero.
        let gas = GasProperties::AIR;
        let port = ValvePort { discharge_coefficient: 0.0, ..right_port() };
        let interior = PrimitiveState::from_pressure_temperature(2.0e5, 350.0, 0.0, &gas);

        let (_, _, valve_port_flux) = valve_port_face(&port, 0.0, interior, 1.0e-3, 1.0e5, 300.0, &gas);
        let closed_end_ghost = BoundaryCondition::ClosedEnd.ghost_state(interior, &gas);
        let closed_end_flux = hllc_flux(interior, closed_end_ghost, &gas);

        assert!((valve_port_flux.mass - closed_end_flux.mass).abs() < 1e-9);
        assert!((valve_port_flux.momentum - closed_end_flux.momentum).abs() < 1e-6);
        assert!((valve_port_flux.energy - closed_end_flux.energy).abs() < 1e-6);
    }

    #[test]
    fn closed_valve_flux_diverges_from_closed_end_hllc_flux_when_the_pipe_has_nonzero_velocity() {
        // The real, documented limitation from the module doc comment: a
        // moving interior reflecting off a closed valve generates a real
        // dynamic-pressure term HLLC captures and this construction (built
        // from the pipe's own static pressure alone) does not. Demonstrated
        // numerically rather than assumed - the two should meaningfully
        // disagree here, not silently match by coincidence.
        let gas = GasProperties::AIR;
        let port = ValvePort { discharge_coefficient: 0.0, ..right_port() };
        let interior = PrimitiveState::from_pressure_temperature(2.0e5, 350.0, 60.0, &gas);

        let (_, _, valve_port_flux) = valve_port_face(&port, 0.0, interior, 1.0e-3, 1.0e5, 300.0, &gas);
        let closed_end_ghost = BoundaryCondition::ClosedEnd.ghost_state(interior, &gas);
        let closed_end_flux = hllc_flux(interior, closed_end_ghost, &gas);

        let relative_difference =
            (valve_port_flux.momentum - closed_end_flux.momentum).abs() / closed_end_flux.momentum.abs();
        assert!(
            relative_difference > 1e-3,
            "expected a meaningfully different momentum flux for a moving interior, got relative difference {relative_difference:e}"
        );
    }
}
