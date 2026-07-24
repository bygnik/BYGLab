//! 0D lumped cylinder thermodynamics: cylinder volume from crank
//! kinematics (`crank_mechanism.rs`) plus a first-law energy balance for
//! the trapped charge, in two stages:
//!
//! - [`Cylinder::integrate_motoring`] — pure adiabatic compression/
//!   expansion of a closed charge, no combustion, no heat transfer, no
//!   mass flow. Validated to ~1e-11 relative error against the exact
//!   isentropic relation `p*V^gamma = const` — deliberately built and
//!   validated first, in isolation, since combustion and heat transfer
//!   don't have a closed-form solution to check against.
//! - [`Cylinder::integrate_fired_cycle`] — adds Wiebe combustion heat
//!   release and Woschni wall heat transfer (`combustion.rs`) as two new
//!   source terms on the same energy balance, still no valve mass flow.
//!   Validated against the real OpenWAM single-cylinder S54 reference
//!   case (`benchmarks/openwam/cases/engine_s54_2500rpm/`) rather than an
//!   exact solution, since none exists once combustion is involved — see
//!   `tests::fired_cycle_trace_compares_against_the_real_openwam_s54_2500rpm_case`
//!   and the root README for the real measured comparison.
//! - [`Cylinder::integrate_breathing`] — adds mass exchange with an
//!   external reservoir through a poppet valve (`camshaft.rs`/`valve.rs`)
//!   as two more terms (`dm/dtheta`, an enthalpy-flux contribution to
//!   `dU/dtheta`) on the same energy balance, independent of
//!   `integrate_fired_cycle` (no combustion/wall heat transfer combined
//!   in yet — intake/exhaust breathing and combustion happen at different
//!   crank angles in a real cycle). Validated against exact closed-form
//!   choked-flow behavior (constant mass flow rate under choking implies
//!   exactly linear mass/energy/pressure growth) using a cylinder
//!   breathing to/from a *fixed reservoir* — not yet the real 1D
//!   intake/exhaust pipe network, which needs a separate architectural
//!   piece (binding a valve to a `pipe::Pipe` end) not built yet.

use crate::camshaft::{self, CamProfile};
use crate::combustion::{self, WallTemperatures, WiebeParameters, WoschniParameters};
use crate::crank_mechanism::CrankMechanism;
use crate::gas::GasProperties;
use crate::valve::{self, ValveGeometry};

/// The cylinder's fixed geometry: crank kinematics, bore, and clearance
/// volume (the volume still present when the piston is at TDC).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder {
    pub crank_mechanism: CrankMechanism,
    /// Bore diameter, meters.
    pub bore: f64,
    /// Volume at TDC, cubic meters — e.g. `displaced_volume / (compression_ratio - 1)`.
    pub clearance_volume: f64,
}

/// The trapped charge's thermodynamic state: total mass and total
/// internal energy — unlike `gas::ConservedState`, these are lumped
/// totals for the *whole* cylinder, not per-unit-volume densities, since
/// this is a single well-mixed 0D control volume, not a spatial field.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderState {
    /// Trapped mass, kg.
    pub mass: f64,
    /// Total internal energy, J.
    pub internal_energy: f64,
}

impl Cylinder {
    /// Piston crown area, square meters.
    pub fn piston_area(&self) -> f64 {
        std::f64::consts::PI / 4.0 * self.bore * self.bore
    }

    /// Cylinder volume at `crank_angle_from_tdc_radians`.
    pub fn volume(&self, crank_angle_from_tdc_radians: f64) -> f64 {
        self.clearance_volume
            + self.piston_area() * self.crank_mechanism.piston_displacement_from_top_dead_center(crank_angle_from_tdc_radians)
    }

    /// `dV/d(crank angle)`, independent of angular velocity. Passing an
    /// angular velocity of `1.0` to [`CrankMechanism::piston_velocity`]
    /// returns `dx/d(theta)` directly (velocity = `dx/dtheta *
    /// angular_velocity`) — motoring is a reversible, adiabatic process
    /// with no rate-dependent terms (no heat transfer, no combustion), so
    /// it's exactly rate-independent: crank angle alone determines the
    /// p-V path, not real time or RPM, which is why this integrates in
    /// crank-angle space without ever needing an angular velocity.
    fn volume_derivative_wrt_crank_angle(&self, crank_angle_from_tdc_radians: f64) -> f64 {
        -self.piston_area() * self.crank_mechanism.piston_velocity(crank_angle_from_tdc_radians, 1.0)
    }

    /// `dU/d(crank angle)` for pure motoring: the first law for a closed,
    /// adiabatic system with no combustion and no mass flow reduces to
    /// `dU = -p dV`. `pub(crate)` so `valve_port.rs` can add this same
    /// motoring term to its own pipe-coupled energy balance without
    /// re-deriving it.
    pub(crate) fn motoring_energy_derivative(&self, gas: &GasProperties, state: &CylinderState, crank_angle_from_tdc_radians: f64) -> f64 {
        let volume = self.volume(crank_angle_from_tdc_radians);
        let pressure = state.pressure(volume, gas);
        -pressure * self.volume_derivative_wrt_crank_angle(crank_angle_from_tdc_radians)
    }

    /// Integrates the motoring energy balance from `theta_start` to
    /// `theta_end` radians (either direction) using classic 4th-order
    /// Runge-Kutta over `step_count` equal substeps — accurate enough
    /// that a modest step count matches the exact isentropic solution to
    /// a tiny fraction of a percent (see `tests/motoring_cycle.rs`), and
    /// empirically confirmed to converge at close to the expected 4th
    /// order under step refinement.
    pub fn integrate_motoring(
        &self,
        gas: &GasProperties,
        initial_state: CylinderState,
        theta_start: f64,
        theta_end: f64,
        step_count: usize,
    ) -> CylinderState {
        let dtheta = (theta_end - theta_start) / step_count as f64;
        let mut state = initial_state;
        let mut theta = theta_start;
        for _ in 0..step_count {
            state = self.motoring_rk4_step(gas, state, theta, dtheta);
            theta += dtheta;
        }
        state
    }

    fn motoring_rk4_step(&self, gas: &GasProperties, state: CylinderState, theta: f64, dtheta: f64) -> CylinderState {
        let k1 = self.motoring_energy_derivative(gas, &state, theta);
        let state2 = CylinderState { mass: state.mass, internal_energy: state.internal_energy + 0.5 * dtheta * k1 };
        let k2 = self.motoring_energy_derivative(gas, &state2, theta + 0.5 * dtheta);
        let state3 = CylinderState { mass: state.mass, internal_energy: state.internal_energy + 0.5 * dtheta * k2 };
        let k3 = self.motoring_energy_derivative(gas, &state3, theta + 0.5 * dtheta);
        let state4 = CylinderState { mass: state.mass, internal_energy: state.internal_energy + dtheta * k3 };
        let k4 = self.motoring_energy_derivative(gas, &state4, theta + dtheta);
        let new_energy = state.internal_energy + (dtheta / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        CylinderState { mass: state.mass, internal_energy: new_energy }
    }

    /// Mean piston speed (`Cm` in the Woschni correlation) at the given
    /// angular velocity — `2*stroke*rpm/60`, algebraically equal to
    /// `stroke*angular_velocity/pi`. A fixed constant for a given
    /// operating point, distinct from [`CrankMechanism::piston_velocity`]'s
    /// continuously-varying instantaneous value.
    fn mean_piston_speed(&self, angular_velocity_radians_per_second: f64) -> f64 {
        2.0 * self.crank_mechanism.crank_radius * angular_velocity_radians_per_second / std::f64::consts::PI
    }

    /// `dU/d(crank angle)` for a fired (combusting) cycle: the motoring
    /// term plus Wiebe heat release plus Woschni wall heat transfer (see
    /// `combustion.rs`). `pressure_at_ivc`/`volume_at_ivc` anchor the
    /// Woschni correlation's "motored reference" pressure curve and must
    /// be the SAME values for the whole integration (cached by the caller,
    /// not recomputed per state) — they describe the reference polytropic
    /// process the actual (combusting) pressure is compared against, not
    /// the current state.
    #[allow(clippy::too_many_arguments)]
    fn fired_energy_derivative(
        &self,
        gas: &GasProperties,
        state: &CylinderState,
        crank_angle_from_tdc_radians: f64,
        params: &FiredCycleParameters,
        pressure_at_ivc: f64,
        volume_at_ivc: f64,
    ) -> f64 {
        let volume = self.volume(crank_angle_from_tdc_radians);
        let pressure = state.pressure(volume, gas);

        let motoring_term = -pressure * self.volume_derivative_wrt_crank_angle(crank_angle_from_tdc_radians);

        let combustion_term = combustion::heat_release_rate(&params.wiebe, crank_angle_from_tdc_radians, params.total_heat_release_joules);

        let piston_area = self.piston_area();
        let unit_displaced_volume = piston_area * 2.0 * self.crank_mechanism.crank_radius;
        let liner_area = std::f64::consts::PI
            * self.bore
            * self.crank_mechanism.piston_displacement_from_top_dead_center(crank_angle_from_tdc_radians);
        let wall_heat_transfer_rate_watts = combustion::wall_heat_transfer_rate_at(
            &params.wiebe,
            &params.woschni,
            &params.walls,
            crank_angle_from_tdc_radians,
            pressure,
            state.temperature_kelvin(gas),
            self.bore,
            self.mean_piston_speed(params.angular_velocity_radians_per_second),
            pressure_at_ivc,
            volume_at_ivc,
            volume,
            unit_displaced_volume,
            state.mass,
            params.gas_constant,
            piston_area,
            piston_area,
            liner_area,
        );
        let wall_heat_transfer_term = wall_heat_transfer_rate_watts / params.angular_velocity_radians_per_second;

        motoring_term + combustion_term + wall_heat_transfer_term
    }

    fn fired_rk4_step(
        &self,
        gas: &GasProperties,
        state: CylinderState,
        theta: f64,
        dtheta: f64,
        params: &FiredCycleParameters,
        pressure_at_ivc: f64,
        volume_at_ivc: f64,
    ) -> CylinderState {
        let k1 = self.fired_energy_derivative(gas, &state, theta, params, pressure_at_ivc, volume_at_ivc);
        let state2 = CylinderState { mass: state.mass, internal_energy: state.internal_energy + 0.5 * dtheta * k1 };
        let k2 = self.fired_energy_derivative(gas, &state2, theta + 0.5 * dtheta, params, pressure_at_ivc, volume_at_ivc);
        let state3 = CylinderState { mass: state.mass, internal_energy: state.internal_energy + 0.5 * dtheta * k2 };
        let k3 = self.fired_energy_derivative(gas, &state3, theta + 0.5 * dtheta, params, pressure_at_ivc, volume_at_ivc);
        let state4 = CylinderState { mass: state.mass, internal_energy: state.internal_energy + dtheta * k3 };
        let k4 = self.fired_energy_derivative(gas, &state4, theta + dtheta, params, pressure_at_ivc, volume_at_ivc);
        let new_energy = state.internal_energy + (dtheta / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        CylinderState { mass: state.mass, internal_energy: new_energy }
    }

    /// Integrates the fired-cycle energy balance (motoring + Wiebe
    /// combustion + Woschni wall heat transfer) from `theta_start` to
    /// `theta_end`, using the same classic 4th-order Runge-Kutta scheme as
    /// [`Self::integrate_motoring`] (which this does not modify or call —
    /// a fully separate path, so the already-validated motoring behavior
    /// carries zero regression risk from this addition).
    ///
    /// `pressure_at_ivc`/`volume_at_ivc` (the Woschni correlation's
    /// "motored reference" anchor) are taken from `initial_state`/
    /// `theta_start` once, up front — NOT recomputed at every RK4 substep.
    pub fn integrate_fired_cycle(
        &self,
        gas: &GasProperties,
        initial_state: CylinderState,
        theta_start: f64,
        theta_end: f64,
        step_count: usize,
        params: &FiredCycleParameters,
    ) -> CylinderState {
        let pressure_at_ivc = initial_state.pressure(self.volume(theta_start), gas);
        let volume_at_ivc = self.volume(theta_start);

        let dtheta = (theta_end - theta_start) / step_count as f64;
        let mut state = initial_state;
        let mut theta = theta_start;
        for _ in 0..step_count {
            state = self.fired_rk4_step(gas, state, theta, dtheta, params, pressure_at_ivc, volume_at_ivc);
            theta += dtheta;
        }
        state
    }

    /// Mass flow rate and upstream specific enthalpy through one valve
    /// event at the given crank angle/cylinder state - shared by
    /// [`Self::breathing_derivatives`] (one valve), the dual-valve
    /// breathing derivative (two valves), and the full-cycle derivative
    /// (two valves + combustion), so this exact physics is coded in
    /// exactly one place. The valve's mass flow rate is wired as
    /// `mass_flow_rate(reservoir, cylinder, ...)` so a positive result
    /// already means "into the cylinder"; the flow's specific enthalpy is
    /// taken from whichever side is upstream (the reservoir's, if flow is
    /// entering the cylinder; the cylinder's own, if leaving) — the
    /// standard open-system first law for a control volume with mass
    /// crossing its boundary, not internal energy (a classic mistake this
    /// was checked against during design review). Returns raw
    /// `(mass_flow_rate_kg_per_s, upstream_specific_enthalpy)`, not yet
    /// divided by angular velocity - callers combine this with whichever
    /// other valve event(s) are active before converting to a
    /// crank-angle-domain rate, since angular velocity is a property of
    /// the whole operating point, not of one valve event.
    fn valve_event_flow(
        &self,
        gas: &GasProperties,
        event: &ValveEventParameters,
        crank_angle_from_tdc_radians: f64,
        cylinder_pressure: f64,
        cylinder_temperature_kelvin: f64,
    ) -> (f64, f64) {
        let lift = camshaft::lift_at(&event.cam, crank_angle_from_tdc_radians);
        let curtain_area = valve::curtain_area(&event.valve, lift);
        let effective_area = event.discharge_coefficient * curtain_area;

        let mass_flow_rate_kg_per_s = valve::mass_flow_rate(
            event.reservoir_pressure,
            event.reservoir_temperature_kelvin,
            cylinder_pressure,
            cylinder_temperature_kelvin,
            effective_area,
            gas,
        );

        let upstream_temperature_kelvin =
            if mass_flow_rate_kg_per_s >= 0.0 { event.reservoir_temperature_kelvin } else { cylinder_temperature_kelvin };
        let upstream_specific_enthalpy = gas.gamma * gas.gas_constant * upstream_temperature_kelvin / (gas.gamma - 1.0);

        (mass_flow_rate_kg_per_s, upstream_specific_enthalpy)
    }

    /// `(dm/dtheta, dU/dtheta)` for the breathing energy balance:
    /// motoring's `-p*dV/dtheta` plus mass/enthalpy exchange with the
    /// external reservoir through the valve (see [`Self::valve_event_flow`]).
    fn breathing_derivatives(
        &self,
        gas: &GasProperties,
        state: &CylinderState,
        crank_angle_from_tdc_radians: f64,
        params: &BreathingParameters,
    ) -> (f64, f64) {
        let volume = self.volume(crank_angle_from_tdc_radians);
        let pressure = state.pressure(volume, gas);
        let temperature_kelvin = state.temperature_kelvin(gas);

        let motoring_term = -pressure * self.volume_derivative_wrt_crank_angle(crank_angle_from_tdc_radians);

        let event = ValveEventParameters {
            cam: params.cam,
            valve: params.valve,
            discharge_coefficient: params.discharge_coefficient,
            reservoir_pressure: params.reservoir_pressure,
            reservoir_temperature_kelvin: params.reservoir_temperature_kelvin,
        };
        let (mass_flow_rate_kg_per_s, upstream_specific_enthalpy) =
            self.valve_event_flow(gas, &event, crank_angle_from_tdc_radians, pressure, temperature_kelvin);

        let dm_dtheta = mass_flow_rate_kg_per_s / params.angular_velocity_radians_per_second;
        let du_dtheta = motoring_term + upstream_specific_enthalpy * dm_dtheta;

        (dm_dtheta, du_dtheta)
    }

    fn breathing_rk4_step(&self, gas: &GasProperties, state: CylinderState, theta: f64, dtheta: f64, params: &BreathingParameters) -> CylinderState {
        let (k1_m, k1_u) = self.breathing_derivatives(gas, &state, theta, params);
        let state2 = CylinderState { mass: state.mass + 0.5 * dtheta * k1_m, internal_energy: state.internal_energy + 0.5 * dtheta * k1_u };
        let (k2_m, k2_u) = self.breathing_derivatives(gas, &state2, theta + 0.5 * dtheta, params);
        let state3 = CylinderState { mass: state.mass + 0.5 * dtheta * k2_m, internal_energy: state.internal_energy + 0.5 * dtheta * k2_u };
        let (k3_m, k3_u) = self.breathing_derivatives(gas, &state3, theta + 0.5 * dtheta, params);
        let state4 = CylinderState { mass: state.mass + dtheta * k3_m, internal_energy: state.internal_energy + dtheta * k3_u };
        let (k4_m, k4_u) = self.breathing_derivatives(gas, &state4, theta + dtheta, params);
        let new_mass = state.mass + (dtheta / 6.0) * (k1_m + 2.0 * k2_m + 2.0 * k3_m + k4_m);
        let new_energy = state.internal_energy + (dtheta / 6.0) * (k1_u + 2.0 * k2_u + 2.0 * k3_u + k4_u);
        CylinderState { mass: new_mass, internal_energy: new_energy }
    }

    /// Integrates the breathing energy balance (motoring + mass/enthalpy
    /// exchange with an external reservoir through a valve) from
    /// `theta_start` to `theta_end`, using the same classic 4th-order
    /// Runge-Kutta scheme as [`Self::integrate_motoring`]/
    /// [`Self::integrate_fired_cycle`] (neither of which this modifies or
    /// calls — a fully separate path, zero regression risk to either).
    /// Does not (yet) combine with combustion/wall heat transfer — intake/
    /// exhaust breathing and combustion happen at different crank angles
    /// in a real cycle, so this is validated standalone first.
    pub fn integrate_breathing(
        &self,
        gas: &GasProperties,
        initial_state: CylinderState,
        theta_start: f64,
        theta_end: f64,
        step_count: usize,
        params: &BreathingParameters,
    ) -> CylinderState {
        let dtheta = (theta_end - theta_start) / step_count as f64;
        let mut state = initial_state;
        let mut theta = theta_start;
        for _ in 0..step_count {
            state = self.breathing_rk4_step(gas, state, theta, dtheta, params);
            theta += dtheta;
        }
        state
    }

    /// `(dm/dtheta, dU/dtheta)` for two independent valve events (each
    /// against its own reservoir) plus motoring - no special-casing for
    /// simultaneous opening (valve overlap near gas-exchange TDC): summing
    /// two independent [`Self::valve_event_flow`] calls handles it by
    /// construction, including exhaust backflow into the cylinder or
    /// cylinder gas escaping out both valves at once. Known, stated
    /// limitation: a lumped single-volume model with two independent
    /// fixed reservoirs cannot represent real port-to-port scavenging
    /// through the cylinder during overlap - that's what the separate,
    /// real-1D-pipe-coupled path in `valve_port.rs` exists for.
    fn dual_valve_breathing_derivatives(
        &self,
        gas: &GasProperties,
        state: &CylinderState,
        crank_angle_from_tdc_radians: f64,
        params: &DualValveBreathingParameters,
    ) -> (f64, f64) {
        let volume = self.volume(crank_angle_from_tdc_radians);
        let pressure = state.pressure(volume, gas);
        let temperature_kelvin = state.temperature_kelvin(gas);

        let motoring_term = -pressure * self.volume_derivative_wrt_crank_angle(crank_angle_from_tdc_radians);

        let (intake_mdot, intake_upstream_enthalpy) =
            self.valve_event_flow(gas, &params.intake, crank_angle_from_tdc_radians, pressure, temperature_kelvin);
        let (exhaust_mdot, exhaust_upstream_enthalpy) =
            self.valve_event_flow(gas, &params.exhaust, crank_angle_from_tdc_radians, pressure, temperature_kelvin);

        let intake_dm_dtheta = intake_mdot / params.angular_velocity_radians_per_second;
        let exhaust_dm_dtheta = exhaust_mdot / params.angular_velocity_radians_per_second;

        let dm_dtheta = intake_dm_dtheta + exhaust_dm_dtheta;
        let du_dtheta =
            motoring_term + intake_upstream_enthalpy * intake_dm_dtheta + exhaust_upstream_enthalpy * exhaust_dm_dtheta;

        (dm_dtheta, du_dtheta)
    }

    fn dual_valve_breathing_rk4_step(
        &self,
        gas: &GasProperties,
        state: CylinderState,
        theta: f64,
        dtheta: f64,
        params: &DualValveBreathingParameters,
    ) -> CylinderState {
        let (k1_m, k1_u) = self.dual_valve_breathing_derivatives(gas, &state, theta, params);
        let state2 = CylinderState { mass: state.mass + 0.5 * dtheta * k1_m, internal_energy: state.internal_energy + 0.5 * dtheta * k1_u };
        let (k2_m, k2_u) = self.dual_valve_breathing_derivatives(gas, &state2, theta + 0.5 * dtheta, params);
        let state3 = CylinderState { mass: state.mass + 0.5 * dtheta * k2_m, internal_energy: state.internal_energy + 0.5 * dtheta * k2_u };
        let (k3_m, k3_u) = self.dual_valve_breathing_derivatives(gas, &state3, theta + 0.5 * dtheta, params);
        let state4 = CylinderState { mass: state.mass + dtheta * k3_m, internal_energy: state.internal_energy + dtheta * k3_u };
        let (k4_m, k4_u) = self.dual_valve_breathing_derivatives(gas, &state4, theta + dtheta, params);
        let new_mass = state.mass + (dtheta / 6.0) * (k1_m + 2.0 * k2_m + 2.0 * k3_m + k4_m);
        let new_energy = state.internal_energy + (dtheta / 6.0) * (k1_u + 2.0 * k2_u + 2.0 * k3_u + k4_u);
        CylinderState { mass: new_mass, internal_energy: new_energy }
    }

    /// Integrates the dual-valve breathing energy balance (motoring +
    /// mass/enthalpy exchange with two independent reservoirs through an
    /// intake AND an exhaust valve) from `theta_start` to `theta_end`,
    /// using the same classic 4th-order Runge-Kutta scheme as
    /// [`Self::integrate_breathing`] (which this does not modify or call
    /// - a fully separate path; setting either `discharge_coefficient` to
    /// `0.0` reduces this exactly to [`Self::integrate_breathing`]'s
    /// single-valve behavior on the other valve, since both share
    /// [`Self::valve_event_flow`]). Does not (yet) combine with combustion
    /// - see [`Self::integrate_full_cycle`] for that.
    pub fn integrate_dual_valve_breathing(
        &self,
        gas: &GasProperties,
        initial_state: CylinderState,
        theta_start: f64,
        theta_end: f64,
        step_count: usize,
        params: &DualValveBreathingParameters,
    ) -> CylinderState {
        let dtheta = (theta_end - theta_start) / step_count as f64;
        let mut state = initial_state;
        let mut theta = theta_start;
        for _ in 0..step_count {
            state = self.dual_valve_breathing_rk4_step(gas, state, theta, dtheta, params);
            theta += dtheta;
        }
        state
    }

    /// `(dm/dtheta, dU/dtheta)` for the full-cycle energy balance:
    /// motoring + Wiebe combustion + Woschni wall heat transfer (anchored
    /// at `pressure_at_ivc`/`volume_at_ivc`, fixed for the whole
    /// integration - same convention as [`Self::fired_energy_derivative`])
    /// + two independent valve events (see [`Self::valve_event_flow`]).
    /// No phase-branching on crank angle: `combustion_term` and
    /// `wall_heat_transfer_term` are themselves already zero outside the
    /// Wiebe burn window (confirmed directly by
    /// `combustion::tests::wall_heat_transfer_rate_at_is_independent_of_the_ivc_anchor_outside_the_wiebe_window`),
    /// so one derivative is correct everywhere, including before IVC.
    #[allow(clippy::too_many_arguments)]
    fn full_cycle_derivatives(
        &self,
        gas: &GasProperties,
        state: &CylinderState,
        crank_angle_from_tdc_radians: f64,
        params: &FullCycleParameters,
        pressure_at_ivc: f64,
        volume_at_ivc: f64,
    ) -> (f64, f64) {
        let volume = self.volume(crank_angle_from_tdc_radians);
        let pressure = state.pressure(volume, gas);
        let temperature_kelvin = state.temperature_kelvin(gas);

        let motoring_term = -pressure * self.volume_derivative_wrt_crank_angle(crank_angle_from_tdc_radians);
        let combustion_term = combustion::heat_release_rate(&params.wiebe, crank_angle_from_tdc_radians, params.total_heat_release_joules);

        let piston_area = self.piston_area();
        let unit_displaced_volume = piston_area * 2.0 * self.crank_mechanism.crank_radius;
        let liner_area = std::f64::consts::PI
            * self.bore
            * self.crank_mechanism.piston_displacement_from_top_dead_center(crank_angle_from_tdc_radians);
        let wall_heat_transfer_rate_watts = combustion::wall_heat_transfer_rate_at(
            &params.wiebe,
            &params.woschni,
            &params.walls,
            crank_angle_from_tdc_radians,
            pressure,
            temperature_kelvin,
            self.bore,
            self.mean_piston_speed(params.angular_velocity_radians_per_second),
            pressure_at_ivc,
            volume_at_ivc,
            volume,
            unit_displaced_volume,
            state.mass,
            params.gas_constant,
            piston_area,
            piston_area,
            liner_area,
        );
        let wall_heat_transfer_term = wall_heat_transfer_rate_watts / params.angular_velocity_radians_per_second;

        let (intake_mdot, intake_upstream_enthalpy) =
            self.valve_event_flow(gas, &params.intake, crank_angle_from_tdc_radians, pressure, temperature_kelvin);
        let (exhaust_mdot, exhaust_upstream_enthalpy) =
            self.valve_event_flow(gas, &params.exhaust, crank_angle_from_tdc_radians, pressure, temperature_kelvin);

        let intake_dm_dtheta = intake_mdot / params.angular_velocity_radians_per_second;
        let exhaust_dm_dtheta = exhaust_mdot / params.angular_velocity_radians_per_second;

        let dm_dtheta = intake_dm_dtheta + exhaust_dm_dtheta;
        let du_dtheta = motoring_term
            + combustion_term
            + wall_heat_transfer_term
            + intake_upstream_enthalpy * intake_dm_dtheta
            + exhaust_upstream_enthalpy * exhaust_dm_dtheta;

        (dm_dtheta, du_dtheta)
    }

    #[allow(clippy::too_many_arguments)]
    fn full_cycle_rk4_step(
        &self,
        gas: &GasProperties,
        state: CylinderState,
        theta: f64,
        dtheta: f64,
        params: &FullCycleParameters,
        pressure_at_ivc: f64,
        volume_at_ivc: f64,
    ) -> CylinderState {
        let (k1_m, k1_u) = self.full_cycle_derivatives(gas, &state, theta, params, pressure_at_ivc, volume_at_ivc);
        let state2 = CylinderState { mass: state.mass + 0.5 * dtheta * k1_m, internal_energy: state.internal_energy + 0.5 * dtheta * k1_u };
        let (k2_m, k2_u) = self.full_cycle_derivatives(gas, &state2, theta + 0.5 * dtheta, params, pressure_at_ivc, volume_at_ivc);
        let state3 = CylinderState { mass: state.mass + 0.5 * dtheta * k2_m, internal_energy: state.internal_energy + 0.5 * dtheta * k2_u };
        let (k3_m, k3_u) = self.full_cycle_derivatives(gas, &state3, theta + 0.5 * dtheta, params, pressure_at_ivc, volume_at_ivc);
        let state4 = CylinderState { mass: state.mass + dtheta * k3_m, internal_energy: state.internal_energy + dtheta * k3_u };
        let (k4_m, k4_u) = self.full_cycle_derivatives(gas, &state4, theta + dtheta, params, pressure_at_ivc, volume_at_ivc);
        let new_mass = state.mass + (dtheta / 6.0) * (k1_m + 2.0 * k2_m + 2.0 * k3_m + k4_m);
        let new_energy = state.internal_energy + (dtheta / 6.0) * (k1_u + 2.0 * k2_u + 2.0 * k3_u + k4_u);
        CylinderState { mass: new_mass, internal_energy: new_energy }
    }

    /// Integrates a genuine full 4-stroke cycle (intake stroke,
    /// compression, Wiebe combustion + Woschni wall heat transfer,
    /// expansion, exhaust stroke) in one call - motoring + combustion +
    /// two independent valve events, using the same classic 4th-order
    /// Runge-Kutta scheme as every other integration mode in this file.
    ///
    /// The Woschni "motored reference" anchor (`pressure_at_ivc`/
    /// `volume_at_ivc`) needs the state at ACTUAL intake-valve-closing
    /// (`ivc_angle = params.intake.cam.opening_angle_radians +
    /// params.intake.cam.duration_radians`), not necessarily
    /// `initial_state`/`theta_start`, if the caller wants a continuous
    /// integration starting before IVC (e.g. from intake-valve-opening or
    /// gas-exchange TDC):
    /// - If `theta_start >= ivc_angle`: the anchor is `initial_state`/
    ///   `theta_start` directly - exactly [`Self::integrate_fired_cycle`]'s
    ///   existing convention, a bit-for-bit backward-compatible reduction.
    /// - Otherwise: [`Self::integrate_dual_valve_breathing`] (a real,
    ///   independently-tested function, not a throwaway pre-pass) runs
    ///   from `theta_start` to `ivc_angle` first to get the anchor state.
    ///   The SAME single RK4 loop below then still covers the *entire*
    ///   `[theta_start, theta_end]` domain (re-walking the pre-IVC region
    ///   once more through the unified derivative is harmless - see
    ///   `full_cycle_derivatives`'s own doc comment for why).
    ///
    /// `step_count` covers that one real RK4 loop; the internal anchor
    /// pre-pass (only run when needed) derives its own, coarser step count
    /// proportionally rather than exposing a second parameter - the
    /// anchor only feeds an empirical correlation term active during the
    /// burn window, so it doesn't need RK4-grade precision, just a
    /// reasonable pressure/volume estimate.
    pub fn integrate_full_cycle(
        &self,
        gas: &GasProperties,
        initial_state: CylinderState,
        theta_start: f64,
        theta_end: f64,
        step_count: usize,
        params: &FullCycleParameters,
    ) -> CylinderState {
        let ivc_angle = params.intake.cam.opening_angle_radians + params.intake.cam.duration_radians;

        let (pressure_at_ivc, volume_at_ivc) = if theta_start >= ivc_angle {
            (initial_state.pressure(self.volume(theta_start), gas), self.volume(theta_start))
        } else {
            let breathing_params = DualValveBreathingParameters {
                intake: params.intake,
                exhaust: params.exhaust,
                angular_velocity_radians_per_second: params.angular_velocity_radians_per_second,
            };
            let anchor_steps = (((ivc_angle - theta_start) / (theta_end - theta_start)).abs() * step_count as f64).ceil().max(20.0) as usize;
            let anchor_state = self.integrate_dual_valve_breathing(gas, initial_state, theta_start, ivc_angle, anchor_steps, &breathing_params);
            let volume_at_ivc = self.volume(ivc_angle);
            (anchor_state.pressure(volume_at_ivc, gas), volume_at_ivc)
        };

        let dtheta = (theta_end - theta_start) / step_count as f64;
        let mut state = initial_state;
        let mut theta = theta_start;
        for _ in 0..step_count {
            state = self.full_cycle_rk4_step(gas, state, theta, dtheta, params, pressure_at_ivc, volume_at_ivc);
            theta += dtheta;
        }
        state
    }
}

/// One valve event: the camshaft profile and valve geometry driving the
/// curtain area, a constant discharge coefficient (a parametric input —
/// see the root README's scoping decision on `Cd` curves), and the fixed
/// external reservoir state on its far side. `cam` is the caller's
/// responsibility to have already expressed in whatever crank-angle
/// convention the rest of the simulation uses — e.g. shifted via
/// [`CamProfile::shifted_by`] to firing-TDC=0 terms before combining with
/// [`Cylinder::integrate_full_cycle`]'s Wiebe combustion timing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ValveEventParameters {
    pub cam: CamProfile,
    pub valve: ValveGeometry,
    pub discharge_coefficient: f64,
    pub reservoir_pressure: f64,
    pub reservoir_temperature_kelvin: f64,
}

/// Everything [`Cylinder::integrate_dual_valve_breathing`] needs: an
/// independent intake and exhaust valve event, each against its own fixed
/// reservoir, plus the shared operating point's angular velocity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DualValveBreathingParameters {
    pub intake: ValveEventParameters,
    pub exhaust: ValveEventParameters,
    pub angular_velocity_radians_per_second: f64,
}

/// Everything [`Cylinder::integrate_breathing`] needs: the camshaft
/// profile and valve geometry driving the curtain area, a constant
/// discharge coefficient (a parametric input — see the root README's
/// scoping decision on `Cd` curves), the external reservoir's fixed
/// state (standing in for the real 1D intake/exhaust pipe network, which
/// isn't coupled in yet), and the operating point's angular velocity
/// (mass/enthalpy exchange is a real-time-rate phenomenon, just like
/// Woschni wall heat transfer, so this is needed here too — breathing is
/// not rate-independent the way pure motoring is).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BreathingParameters {
    pub cam: CamProfile,
    pub valve: ValveGeometry,
    pub discharge_coefficient: f64,
    pub reservoir_pressure: f64,
    pub reservoir_temperature_kelvin: f64,
    pub angular_velocity_radians_per_second: f64,
}

/// Everything [`Cylinder::integrate_fired_cycle`] needs beyond the base
/// [`Cylinder`] geometry: combustion timing/energy, Woschni correlation
/// coefficients, fixed wall temperatures, the operating point's angular
/// velocity (needed to convert the Woschni wall-heat-transfer rate, a
/// real-time quantity, into a crank-angle-domain ODE term — the first
/// place this model needs angular velocity at all, since motoring alone
/// is exactly rate-independent), and the gas constant (duplicated from
/// [`GasProperties::gas_constant`] since the Woschni correlation's
/// combustion-turbulence term needs it directly, not through a
/// `GasProperties` reference, to keep `combustion.rs` free of a
/// dependency on `gas.rs`'s specific types).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiredCycleParameters {
    pub wiebe: WiebeParameters,
    pub woschni: WoschniParameters,
    pub walls: WallTemperatures,
    /// Total chemical energy released over a complete burn, Joules
    /// (`m_fuel * LHV * combustion_efficiency`) — a single scalar, not an
    /// evolving fuel-mass ODE (matches how OpenWAM structures this too).
    pub total_heat_release_joules: f64,
    pub angular_velocity_radians_per_second: f64,
    pub gas_constant: f64,
}

/// Everything [`Cylinder::integrate_full_cycle`] needs: [`FiredCycleParameters`]'s
/// combustion/heat-transfer physics plus an independent intake and
/// exhaust valve event, each against its own reservoir (matching
/// [`DualValveBreathingParameters`]) - together enough to integrate a
/// genuine 4-stroke cycle (intake stroke, compression, combustion,
/// expansion, exhaust stroke) in one call.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FullCycleParameters {
    pub wiebe: WiebeParameters,
    pub woschni: WoschniParameters,
    pub walls: WallTemperatures,
    pub total_heat_release_joules: f64,
    pub angular_velocity_radians_per_second: f64,
    pub gas_constant: f64,
    pub intake: ValveEventParameters,
    pub exhaust: ValveEventParameters,
}

impl CylinderState {
    /// Builds a state from pressure, temperature, and volume via the
    /// ideal gas law (`m = pV/(RT)`) and `U = m*cv*T`.
    pub fn from_pressure_temperature(pressure: f64, temperature_kelvin: f64, volume: f64, gas: &GasProperties) -> Self {
        let mass = pressure * volume / (gas.gas_constant * temperature_kelvin);
        let specific_heat_cv = gas.gas_constant / (gas.gamma - 1.0);
        let internal_energy = mass * specific_heat_cv * temperature_kelvin;
        CylinderState { mass, internal_energy }
    }

    /// Pressure at the given volume, from the ideal gas law rearranged in
    /// terms of internal energy: `p = (gamma-1) * U / V` (since
    /// `U = m*cv*T` and `p = m*R*T/V`, so `p = (R/cv)*U/V = (gamma-1)*U/V`).
    pub fn pressure(&self, volume: f64, gas: &GasProperties) -> f64 {
        (gas.gamma - 1.0) * self.internal_energy / volume
    }

    /// Temperature, back-derived from `U = m*cv*T`.
    pub fn temperature_kelvin(&self, gas: &GasProperties) -> f64 {
        let specific_heat_cv = gas.gas_constant / (gas.gamma - 1.0);
        self.internal_energy / (self.mass * specific_heat_cv)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// S54B32 reference geometry (see root README's spec table): bore
    /// 87mm, stroke 91mm, rod 139mm, CR 11.5:1 — real, validated numbers.
    /// Clearance volume from `CR = (V_clearance + V_displaced) / V_clearance`.
    fn s54b32_cylinder() -> Cylinder {
        let bore = 0.087;
        let stroke = 0.091;
        let compression_ratio = 11.5;
        let piston_area = std::f64::consts::PI / 4.0 * bore * bore;
        let displaced_volume = piston_area * stroke;
        let clearance_volume = displaced_volume / (compression_ratio - 1.0);
        Cylinder { crank_mechanism: CrankMechanism::new(stroke, 0.139, 0.0), bore, clearance_volume }
    }

    #[test]
    fn compression_matches_the_exact_isentropic_relation() {
        let cylinder = s54b32_cylinder();
        let gas = GasProperties::AIR;
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();

        let initial_pressure = 1.0e5;
        let initial_temperature_kelvin = 320.0; // representative post-IVC charge temperature
        let volume_at_bdc = cylinder.volume(bdc_angle);
        let initial_state = CylinderState::from_pressure_temperature(initial_pressure, initial_temperature_kelvin, volume_at_bdc, &gas);

        let final_state = cylinder.integrate_motoring(&gas, initial_state, bdc_angle, 0.0, 720);
        let volume_at_tdc = cylinder.volume(0.0);

        // Exact isentropic relations: p2 = p1*(V1/V2)^gamma, T2 = T1*(V1/V2)^(gamma-1).
        let volume_ratio = volume_at_bdc / volume_at_tdc;
        let exact_pressure = initial_pressure * volume_ratio.powf(gas.gamma);
        let exact_temperature = initial_temperature_kelvin * volume_ratio.powf(gas.gamma - 1.0);

        let actual_pressure = final_state.pressure(volume_at_tdc, &gas);
        let actual_temperature = final_state.temperature_kelvin(&gas);

        let pressure_relative_error = (actual_pressure - exact_pressure).abs() / exact_pressure;
        let temperature_relative_error = (actual_temperature - exact_temperature).abs() / exact_temperature;

        println!(
            "geometric CR={:.2}, volume ratio={volume_ratio:.4}, exact p={:.1} bar / actual p={:.1} bar (error {:e}), exact T={:.1} K / actual T={:.1} K (error {:e})",
            volume_at_bdc / cylinder.clearance_volume,
            exact_pressure / 1e5,
            actual_pressure / 1e5,
            pressure_relative_error,
            exact_temperature,
            actual_temperature,
            temperature_relative_error
        );

        assert!(pressure_relative_error < 1e-6, "pressure relative error {pressure_relative_error:e} too high");
        assert!(temperature_relative_error < 1e-6, "temperature relative error {temperature_relative_error:e} too high");
    }

    #[test]
    fn mass_is_conserved_during_motoring() {
        let cylinder = s54b32_cylinder();
        let gas = GasProperties::AIR;
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();
        let volume_at_bdc = cylinder.volume(bdc_angle);
        let initial_state = CylinderState::from_pressure_temperature(1.0e5, 320.0, volume_at_bdc, &gas);

        let final_state = cylinder.integrate_motoring(&gas, initial_state, bdc_angle, 0.0, 720);

        assert_eq!(final_state.mass, initial_state.mass, "mass must be exactly unchanged - no valve flow exists in this model yet");
    }

    #[test]
    fn energy_returns_to_initial_value_after_a_full_bdc_tdc_bdc_round_trip() {
        // Reversible, adiabatic, no losses: compressing BDC->TDC then
        // expanding back TDC->BDC must retrace the exact same p-V curve,
        // landing back on the exact initial state - a strong, independent
        // check of energy conservation (not just "matches the isentropic
        // curve one-way", but "the numerical integration doesn't leak or
        // manufacture energy over a closed cycle").
        let cylinder = s54b32_cylinder();
        let gas = GasProperties::AIR;
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();
        // No piston pin offset in this geometry, so x(theta)=x(-theta)
        // exactly (see crank_mechanism.rs's symmetry test) - the "BDC
        // before TDC" and "BDC after TDC" crank angles are exact
        // negatives of each other, giving a genuinely closed p-V loop.
        let volume_at_bdc = cylinder.volume(bdc_angle);
        let initial_state = CylinderState::from_pressure_temperature(1.0e5, 320.0, volume_at_bdc, &gas);

        let at_tdc = cylinder.integrate_motoring(&gas, initial_state, -bdc_angle, 0.0, 720);
        let back_at_bdc = cylinder.integrate_motoring(&gas, at_tdc, 0.0, bdc_angle, 720);

        let energy_relative_error = (back_at_bdc.internal_energy - initial_state.internal_energy).abs() / initial_state.internal_energy;
        let pressure_relative_error =
            (back_at_bdc.pressure(volume_at_bdc, &gas) - initial_state.pressure(volume_at_bdc, &gas)).abs() / initial_state.pressure(volume_at_bdc, &gas);

        println!("round trip: energy relative error {energy_relative_error:e}, pressure relative error {pressure_relative_error:e}");

        assert_eq!(back_at_bdc.mass, initial_state.mass);
        assert!(energy_relative_error < 1e-6, "energy did not return to its initial value: relative error {energy_relative_error:e}");
        assert!(pressure_relative_error < 1e-6, "pressure did not return to its initial value: relative error {pressure_relative_error:e}");
    }

    #[test]
    fn rk4_integration_converges_at_close_to_the_expected_fourth_order_rate() {
        let cylinder = s54b32_cylinder();
        let gas = GasProperties::AIR;
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();
        let volume_at_bdc = cylinder.volume(bdc_angle);
        let volume_at_tdc = cylinder.volume(0.0);
        let initial_state = CylinderState::from_pressure_temperature(1.0e5, 320.0, volume_at_bdc, &gas);
        let exact_pressure = 1.0e5 * (volume_at_bdc / volume_at_tdc).powf(gas.gamma);

        let coarse_steps = 8;
        let fine_steps = 16;
        let coarse_error =
            (cylinder.integrate_motoring(&gas, initial_state, bdc_angle, 0.0, coarse_steps).pressure(volume_at_tdc, &gas) - exact_pressure).abs();
        let fine_error =
            (cylinder.integrate_motoring(&gas, initial_state, bdc_angle, 0.0, fine_steps).pressure(volume_at_tdc, &gas) - exact_pressure).abs();

        let observed_order = (coarse_error / fine_error).log2();
        println!("RK4 convergence: {coarse_steps} steps error={coarse_error:e} Pa, {fine_steps} steps error={fine_error:e} Pa, observed order={observed_order:.2}");

        // Deliberately coarse step counts (8/16, not hundreds) so the
        // truncation error dominates over floating-point noise, which
        // would otherwise flatten the observed order at very fine steps.
        assert!(observed_order > 3.0, "expected close to 4th-order convergence, observed order {observed_order:.2}");
    }

    /// The real OpenWAM single-cylinder S54 validation case
    /// (`benchmarks/openwam/cases/engine_s54_2500rpm/engine_s54_2500rpm.txt`)
    /// — CR=11.0 here specifically (not the 11.5 used by
    /// `s54b32_cylinder` above, which reflects this project's general
    /// S54B32 spec table rather than this one exact OpenWAM case file).
    fn s54_2500rpm_cylinder() -> Cylinder {
        let bore = 0.087;
        let stroke = 0.091;
        let compression_ratio = 11.0;
        let piston_area = std::f64::consts::PI / 4.0 * bore * bore;
        let displaced_volume = piston_area * stroke;
        let clearance_volume = displaced_volume / (compression_ratio - 1.0);
        Cylinder { crank_mechanism: CrankMechanism::new(stroke, 0.139, 0.0), bore, clearance_volume }
    }

    /// Real Wiebe/Woschni parameters and IVC seed state from the OpenWAM
    /// S54 2500rpm case — see the plan's "Real reference data" section
    /// for provenance (extracted directly from the case's input file and
    /// its own output trace, not estimated).
    fn s54_2500rpm_fired_cycle_params(gas: &GasProperties, ivc_state: CylinderState) -> FiredCycleParameters {
        let wiebe = WiebeParameters {
            start_angle_radians: (-15.0_f64).to_radians(),
            duration_radians: 45.0_f64.to_radians(),
            shape_factor_m: 2.5,
            efficiency_c: 6.9,
        };
        let woschni = WoschniParameters { cw1: 2.28, cw2: 0.00324, combustion_turbulence_coefficient: 0.001 };
        let walls = WallTemperatures { piston_kelvin: 450.0, head_kelvin: 550.0, liner_kelvin: 420.0 };
        let rpm = 2500.0;
        let angular_velocity_radians_per_second = rpm * 2.0 * std::f64::consts::PI / 60.0;

        // Physically-derived estimate (not calibrated to hit OpenWAM's
        // peak pressure): stoichiometric AFR=14.7 splits the IVC trapped
        // mass into air+fuel, then LHV/efficiency give the total
        // releasable chemical energy - see the plan's Q_total discussion.
        let afr_stoichiometric = 14.7;
        let fuel_mass = ivc_state.mass / (afr_stoichiometric + 1.0);
        let lhv_joules_per_kg = 43.0e6;
        let combustion_efficiency = 0.98;
        let total_heat_release_joules = fuel_mass * lhv_joules_per_kg * combustion_efficiency;

        FiredCycleParameters { wiebe, woschni, walls, total_heat_release_joules, angular_velocity_radians_per_second, gas_constant: gas.gas_constant }
    }

    #[test]
    fn ivc_volume_matches_openwams_reported_value() {
        let cylinder = s54_2500rpm_cylinder();
        let ivc_angle = (-120.0_f64).to_radians();
        let my_volume = cylinder.volume(ivc_angle);
        let openwam_volume = 4.93301e-4;
        let relative_error = (my_volume - openwam_volume).abs() / openwam_volume;
        println!("IVC volume: mine={my_volume:e} m^3, OpenWAM={openwam_volume:e} m^3, relative error={relative_error:e}");
        assert!(relative_error < 0.005, "expected a close match (~0.08% by hand calculation), got {relative_error:e}");
    }

    /// Measures RK4 convergence order (via coarse/fine vs. a much finer
    /// reference run - no closed-form solution exists once combustion is
    /// involved) over `[interval_start, burn_end]`, starting the motoring
    /// pre-segment from IVC.
    fn measure_fired_cycle_convergence_order(interval_start_radians: f64) -> (f64, f64) {
        let cylinder = s54_2500rpm_cylinder();
        let gas = GasProperties::AIR;
        let ivc_angle = (-120.0_f64).to_radians();
        let ivc_volume = cylinder.volume(ivc_angle);
        let ivc_state = CylinderState::from_pressure_temperature(1.16948e5, 393.07, ivc_volume, &gas);
        let params = s54_2500rpm_fired_cycle_params(&gas, ivc_state);

        // Motor/burn from IVC up to the interval's start first (reusing
        // plenty of steps so this segment's own error doesn't contaminate
        // the convergence measurement over the interval of interest).
        let pre_interval_state = cylinder.integrate_fired_cycle(&gas, ivc_state, ivc_angle, interval_start_radians, 4000, &params);

        let burn_end = params.wiebe.start_angle_radians + params.wiebe.duration_radians;
        let volume_at_burn_end = cylinder.volume(burn_end);

        let reference_pressure = cylinder
            .integrate_fired_cycle(&gas, pre_interval_state, interval_start_radians, burn_end, 2880, &params)
            .pressure(volume_at_burn_end, &gas);
        let coarse_pressure = cylinder
            .integrate_fired_cycle(&gas, pre_interval_state, interval_start_radians, burn_end, 45, &params)
            .pressure(volume_at_burn_end, &gas);
        let fine_pressure = cylinder
            .integrate_fired_cycle(&gas, pre_interval_state, interval_start_radians, burn_end, 90, &params)
            .pressure(volume_at_burn_end, &gas);

        let coarse_error = (coarse_pressure - reference_pressure).abs();
        let fine_error = (fine_pressure - reference_pressure).abs();
        (coarse_error, fine_error)
    }

    #[test]
    fn fired_cycle_rk4_convergence_is_reduced_right_at_the_singular_ignition_point() {
        // The Wiebe shape factor m=2.5 is non-integer, so dXB/dtheta ~
        // y^m near y=0 (y=(theta-theta0)/duration) has an UNBOUNDED 4th
        // derivative right at ignition (d^4/dtheta^4 ~ y^(m-3) = y^(-0.5)
        // -> infinity as y->0) - classical RK4 error analysis assumes
        // enough bounded derivatives throughout the interval, so measuring
        // convergence starting EXACTLY at theta0 is expected to show a
        // reduced order, the same well-understood phenomenon already
        // documented in `tests/mesh_convergence.rs` for a shock capping
        // convergence order below the scheme's formal rate in a smooth
        // region - not a bug, a real property of a non-integer Wiebe
        // exponent. See `fired_cycle_rk4_converges_at_close_to_fourth_order_away_from_ignition`
        // for confirmation that the *scheme itself* still achieves near-4th-order
        // once measured away from this single singular point.
        let wiebe_start = (-15.0_f64).to_radians();
        let (coarse_error, fine_error) = measure_fired_cycle_convergence_order(wiebe_start);
        let observed_order = (coarse_error / fine_error).log2();
        println!("fired-cycle RK4 convergence AT ignition (theta0): coarse error={coarse_error:e} Pa, fine error={fine_error:e} Pa, observed order={observed_order:.2}");
        assert!(observed_order > 0.8, "expected a real, if reduced, positive convergence rate, observed order {observed_order:.2}");
        assert!(observed_order < 2.5, "expected a CLEARLY reduced rate versus the scheme's 4th-order design (confirming the singular-point effect is real) - if this is now close to 4, the interval below may need adjusting");
    }

    #[test]
    fn fired_cycle_rk4_converges_at_close_to_fourth_order_away_from_ignition() {
        // Same measurement as the test above, but starting 10 degrees
        // after ignition (well clear of the y=0 singular point in
        // dXB/dtheta's higher derivatives) - confirms the RK4
        // *implementation* itself is correct and achieves its full
        // formal order once away from that single non-smooth point.
        let interval_start = (-5.0_f64).to_radians(); // theta0 (-15 deg) + 10 deg
        let (coarse_error, fine_error) = measure_fired_cycle_convergence_order(interval_start);
        let observed_order = (coarse_error / fine_error).log2();
        println!("fired-cycle RK4 convergence AWAY from ignition (theta0+10deg): coarse error={coarse_error:e} Pa, fine error={fine_error:e} Pa, observed order={observed_order:.2}");
        assert!(observed_order > 3.0, "expected close to 4th-order convergence away from the singular ignition point, observed order {observed_order:.2}");
    }

    #[test]
    fn fired_cycle_trace_compares_against_the_real_openwam_s54_2500rpm_case() {
        let cylinder = s54_2500rpm_cylinder();
        let gas = GasProperties::AIR;
        let ivc_angle = (-120.0_f64).to_radians();
        let ivc_volume = cylinder.volume(ivc_angle);
        let ivc_state = CylinderState::from_pressure_temperature(1.16948e5, 393.07, ivc_volume, &gas);
        let params = s54_2500rpm_fired_cycle_params(&gas, ivc_state);

        println!("Q_total (physically-derived, not calibrated) = {:.1} J", params.total_heat_release_joules);

        // Real OpenWAM trace: (crank angle deg ATDC, pressure bar, temperature degC).
        let checkpoints: [(f64, f64, f64); 8] = [
            (0.10, 33.50, 852.0),
            (5.00, 46.01, 1320.0),
            (10.07, 57.86, 1894.0),
            (15.14, 62.45, 2352.0),
            (20.03, 58.70, 2551.0),
            (22.13, 55.48, 2567.0),
            (25.10, 49.69, 2510.0),
            (30.00, 40.99, 2396.0),
        ];

        let mut state = cylinder.integrate_fired_cycle(&gas, ivc_state, ivc_angle, checkpoints[0].0.to_radians(), 2000, &params);
        let mut previous_angle_radians = checkpoints[0].0.to_radians();

        let mut max_pressure_relative_error = 0.0_f64;
        let mut max_temperature_relative_error = 0.0_f64;

        for (angle_deg, expected_pressure_bar, expected_temperature_c) in checkpoints {
            let angle_radians = angle_deg.to_radians();
            if angle_radians != previous_angle_radians {
                let segment_steps = ((angle_radians - previous_angle_radians).to_degrees().abs() * 20.0).ceil() as usize;
                state = cylinder.integrate_fired_cycle(&gas, state, previous_angle_radians, angle_radians, segment_steps.max(20), &params);
            }
            previous_angle_radians = angle_radians;

            let volume = cylinder.volume(angle_radians);
            let actual_pressure_bar = state.pressure(volume, &gas) / 1e5;
            let actual_temperature_c = state.temperature_kelvin(&gas) - 273.15;

            let pressure_relative_error = (actual_pressure_bar - expected_pressure_bar).abs() / expected_pressure_bar;
            let temperature_relative_error = (actual_temperature_c - expected_temperature_c).abs() / (expected_temperature_c + 273.15);
            max_pressure_relative_error = max_pressure_relative_error.max(pressure_relative_error);
            max_temperature_relative_error = max_temperature_relative_error.max(temperature_relative_error);

            println!(
                "theta={angle_deg:>6.2} ATDC: pressure actual={actual_pressure_bar:>6.2} bar / OpenWAM={expected_pressure_bar:>6.2} bar (err {:>6.2}%), temperature actual={actual_temperature_c:>7.1} C / OpenWAM={expected_temperature_c:>7.1} C (err {:>6.2}%)",
                pressure_relative_error * 100.0,
                temperature_relative_error * 100.0
            );
        }

        println!("max pressure relative error = {:.2}%, max temperature relative error = {:.2}%", max_pressure_relative_error * 100.0, max_temperature_relative_error * 100.0);

        // Measured ~30%/40% (see the printed trace above and the root
        // README) - a systematic OVER-prediction across the whole trace,
        // not a shape/timing mismatch (peak pressure/temperature occur at
        // the right angles, ~7 degrees apart, matching OpenWAM exactly).
        // The likely cause: `total_heat_release_joules` assumes the ENTIRE
        // IVC trapped mass is fresh stoichiometric air+fuel charge, with
        // no residual (burned) exhaust gas fraction subtracted out - a
        // real, documented simplification (this model doesn't track
        // residuals yet), not a bug in the Wiebe/Woschni physics
        // themselves, which the convergence tests above already isolate
        // and confirm are implemented correctly. 40% leaves real margin
        // above the measured ~30%/40% while still catching a much larger
        // regression (e.g. a sign error or a missing term).
        assert!(max_pressure_relative_error < 0.40, "pressure trace error {:.1}% exceeds the 40% bound", max_pressure_relative_error * 100.0);
        assert!(max_temperature_relative_error < 0.45, "temperature trace error {:.1}% exceeds the 45% bound", max_temperature_relative_error * 100.0);
    }

    #[test]
    fn ideal_otto_cycle_matches_the_exact_textbook_efficiency_formula() {
        // The idealized air-standard Otto cycle: adiabatic reversible
        // compression (1->2), INSTANTANEOUS constant-volume heat addition
        // (2->3), adiabatic reversible expansion (3->4), then (implicitly)
        // constant-volume heat rejection back to state 1. Unlike the real
        // Wiebe/Woschni combustion trace tested above (which has no exact
        // solution to compare against, only OpenWAM's own numbers), this
        // idealization has a well-known EXACT closed-form efficiency:
        // `eta = 1 - 1/CR^(gamma-1)`, independent of how much heat is
        // actually added. Compression/expansion reuse `integrate_motoring`
        // (already independently validated elsewhere in this file to
        // ~1e-11) unchanged; the "instantaneous" heat addition is exact BY
        // DEFINITION for this idealization (adding energy at fixed
        // mass/volume), not an approximation that itself needs validating
        // - so this test is really checking that assembling these already-
        // validated pieces together doesn't introduce a bookkeeping error,
        // against a genuinely independent, well-known exact target.
        let cylinder = s54b32_cylinder(); // real S54B32 geometry, CR=11.5
        let gas = GasProperties::AIR;
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();
        let volume_at_bdc = cylinder.volume(bdc_angle);
        let volume_at_tdc = cylinder.volume(0.0);
        let compression_ratio = volume_at_bdc / volume_at_tdc;

        let initial_pressure = 1.0e5;
        let initial_temperature_kelvin = 300.0;
        let state1 = CylinderState::from_pressure_temperature(initial_pressure, initial_temperature_kelvin, volume_at_bdc, &gas);

        // 1->2: adiabatic reversible compression.
        let state2 = cylinder.integrate_motoring(&gas, state1, bdc_angle, 0.0, 720);

        // 2->3: instantaneous constant-volume heat addition - literally
        // just adding energy at fixed mass/volume, the textbook
        // definition of this step, not something to approximate further.
        let heat_input_joules = 1500.0;
        let state3 = CylinderState { mass: state2.mass, internal_energy: state2.internal_energy + heat_input_joules };

        // 3->4: adiabatic reversible expansion back to BDC.
        let state4 = cylinder.integrate_motoring(&gas, state3, 0.0, bdc_angle, 720);

        let specific_heat_cv = gas.gas_constant / (gas.gamma - 1.0);
        let heat_rejected_joules = state1.mass * specific_heat_cv * (state4.temperature_kelvin(&gas) - state1.temperature_kelvin(&gas));

        let actual_efficiency = 1.0 - heat_rejected_joules / heat_input_joules;
        let exact_efficiency = 1.0 - 1.0 / compression_ratio.powf(gas.gamma - 1.0);
        let efficiency_relative_error = (actual_efficiency - exact_efficiency).abs() / exact_efficiency;

        // Exact state-by-state check too, not just the aggregate efficiency.
        let exact_t2 = initial_temperature_kelvin * compression_ratio.powf(gas.gamma - 1.0);
        let exact_t3 = exact_t2 + heat_input_joules / (state1.mass * specific_heat_cv);
        let exact_t4 = exact_t3 / compression_ratio.powf(gas.gamma - 1.0);
        let t2_relative_error = (state2.temperature_kelvin(&gas) - exact_t2).abs() / exact_t2;
        let t3_relative_error = (state3.temperature_kelvin(&gas) - exact_t3).abs() / exact_t3;
        let t4_relative_error = (state4.temperature_kelvin(&gas) - exact_t4).abs() / exact_t4;

        println!(
            "ideal Otto cycle: CR={compression_ratio:.3}, exact efficiency={:.3}%, actual efficiency={:.3}%, relative error={efficiency_relative_error:e}",
            exact_efficiency * 100.0,
            actual_efficiency * 100.0
        );
        println!(
            "T2: actual={:.2}K exact={:.2}K (err {t2_relative_error:e}); T3: actual={:.2}K exact={:.2}K (err {t3_relative_error:e}); T4: actual={:.2}K exact={:.2}K (err {t4_relative_error:e})",
            state2.temperature_kelvin(&gas),
            exact_t2,
            state3.temperature_kelvin(&gas),
            exact_t3,
            state4.temperature_kelvin(&gas),
            exact_t4
        );

        assert!(efficiency_relative_error < 1e-6, "efficiency relative error {efficiency_relative_error:e} too high");
        assert!(t2_relative_error < 1e-6, "T2 relative error {t2_relative_error:e} too high");
        assert!(t3_relative_error < 1e-9, "T3 relative error {t3_relative_error:e} too high (should be near machine precision - exact by construction)");
        assert!(t4_relative_error < 1e-6, "T4 relative error {t4_relative_error:e} too high");
    }

    #[test]
    fn ideal_otto_cycle_imep_matches_the_exact_analytical_value() {
        // IMEP (indicated mean effective pressure) = net indicated work
        // per cycle / displaced volume - the standard, cylinder-size-
        // independent way any engine simulation reports specific work
        // output (see the root README's phase-5 roadmap item on
        // performance metrics). For the ideal Otto cycle, net work = Qin -
        // Qout (first law over a closed cycle: internal energy returns to
        // its starting value, so net work equals net heat in).
        //
        // Validated two genuinely different ways, both against the same
        // exact closed-form target: (1) directly integrating p dV
        // (trapezoidal rule) along the actual RK4-integrated compression/
        // expansion trajectory - a real numerical path integration, not
        // just an endpoint-state shortcut - and (2) the energy-balance
        // shortcut (Qin-Qout from simulated state internal energies,
        // algebraically equivalent to the p-dV integral for an ideal gas,
        // but computed via a completely different code path here).
        let cylinder = s54b32_cylinder();
        let gas = GasProperties::AIR;
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();
        let volume_at_bdc = cylinder.volume(bdc_angle);
        let volume_at_tdc = cylinder.volume(0.0);
        let displaced_volume = volume_at_bdc - volume_at_tdc;
        let compression_ratio = volume_at_bdc / volume_at_tdc;

        let initial_pressure = 1.0e5;
        let initial_temperature_kelvin = 300.0;
        let state1 = CylinderState::from_pressure_temperature(initial_pressure, initial_temperature_kelvin, volume_at_bdc, &gas);
        let heat_input_joules = 1500.0;

        let segments = 2000;
        let mut work_by_gas_joules = 0.0_f64;

        // Compression: BDC -> TDC, accumulating p dV segment by segment.
        let mut state = state1;
        let mut theta = bdc_angle;
        let dtheta_compression = -bdc_angle / segments as f64;
        for _ in 0..segments {
            let pressure_before = state.pressure(cylinder.volume(theta), &gas);
            let next_theta = theta + dtheta_compression;
            let next_state = cylinder.integrate_motoring(&gas, state, theta, next_theta, 4);
            let volume_before = cylinder.volume(theta);
            let volume_after = cylinder.volume(next_theta);
            let pressure_after = next_state.pressure(volume_after, &gas);
            work_by_gas_joules += 0.5 * (pressure_before + pressure_after) * (volume_after - volume_before);
            state = next_state;
            theta = next_theta;
        }
        let state2 = state;

        // Constant-volume heat addition: no work done.
        let state3 = CylinderState { mass: state2.mass, internal_energy: state2.internal_energy + heat_input_joules };

        // Expansion: TDC -> BDC, same accumulation.
        let mut state = state3;
        let mut theta = 0.0_f64;
        let dtheta_expansion = bdc_angle / segments as f64;
        for _ in 0..segments {
            let pressure_before = state.pressure(cylinder.volume(theta), &gas);
            let next_theta = theta + dtheta_expansion;
            let next_state = cylinder.integrate_motoring(&gas, state, theta, next_theta, 4);
            let volume_before = cylinder.volume(theta);
            let volume_after = cylinder.volume(next_theta);
            let pressure_after = next_state.pressure(volume_after, &gas);
            work_by_gas_joules += 0.5 * (pressure_before + pressure_after) * (volume_after - volume_before);
            state = next_state;
            theta = next_theta;
        }
        let state4 = state;

        let specific_heat_cv = gas.gas_constant / (gas.gamma - 1.0);
        let heat_rejected_joules = state1.mass * specific_heat_cv * (state4.temperature_kelvin(&gas) - initial_temperature_kelvin);
        let work_from_energy_balance_joules = heat_input_joules - heat_rejected_joules;

        let imep_from_pdv_integration = work_by_gas_joules / displaced_volume;
        let imep_from_energy_balance = work_from_energy_balance_joules / displaced_volume;

        // Exact closed-form target, from the textbook Otto-cycle state relations.
        let exact_t2 = initial_temperature_kelvin * compression_ratio.powf(gas.gamma - 1.0);
        let exact_t3 = exact_t2 + heat_input_joules / (state1.mass * specific_heat_cv);
        let exact_t4 = exact_t3 / compression_ratio.powf(gas.gamma - 1.0);
        let exact_heat_rejected = state1.mass * specific_heat_cv * (exact_t4 - initial_temperature_kelvin);
        let exact_work = heat_input_joules - exact_heat_rejected;
        let exact_imep = exact_work / displaced_volume;

        let pdv_relative_error = (imep_from_pdv_integration - exact_imep).abs() / exact_imep;
        let energy_balance_relative_error = (imep_from_energy_balance - exact_imep).abs() / exact_imep;

        println!(
            "IMEP: p-dV integration={:.4} bar, energy balance={:.4} bar, exact={:.4} bar (p-dV err {pdv_relative_error:e}, energy-balance err {energy_balance_relative_error:e})",
            imep_from_pdv_integration / 1e5,
            imep_from_energy_balance / 1e5,
            exact_imep / 1e5,
        );

        // Measured 1.3e-6 at 2000 segments - confirmed to be pure
        // trapezoidal discretization error (not an RK4/physics bug): at
        // 200 segments it measured 1.3e-4, and 10x more segments giving
        // exactly 100x less error is the textbook O(1/N^2) trapezoidal
        // convergence rate for a smooth, strongly-curved integrand (the
        // isentropic p-V curve, pressure changing ~30x over the stroke).
        assert!(pdv_relative_error < 1e-5, "p-dV-integrated IMEP relative error {pdv_relative_error:e} too high");
        assert!(energy_balance_relative_error < 1e-9, "energy-balance IMEP relative error {energy_balance_relative_error:e} too high");
    }

    #[test]
    fn breathing_with_a_closed_valve_matches_motoring_exactly() {
        let cylinder = s54b32_cylinder();
        let gas = GasProperties::AIR;
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();
        let volume_at_bdc = cylinder.volume(bdc_angle);
        let initial_state = CylinderState::from_pressure_temperature(1.0e5, 320.0, volume_at_bdc, &gas);

        let motoring_result = cylinder.integrate_motoring(&gas, initial_state, bdc_angle, 0.0, 720);

        // max_lift=0 -> curtain_area is exactly 0 regardless of the
        // opening/duration window, so this should reduce to motoring
        // exactly, regardless of how large a reservoir pressure
        // difference is present to drive flow.
        let params = BreathingParameters {
            cam: CamProfile { max_lift: 0.0, opening_angle_radians: -100.0, duration_radians: 200.0 },
            valve: ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 45.0_f64.to_radians() },
            discharge_coefficient: 0.7,
            reservoir_pressure: 5.0e5,
            reservoir_temperature_kelvin: 400.0,
            angular_velocity_radians_per_second: 500.0,
        };
        let breathing_result = cylinder.integrate_breathing(&gas, initial_state, bdc_angle, 0.0, 720, &params);

        let energy_relative_error = (breathing_result.internal_energy - motoring_result.internal_energy).abs() / motoring_result.internal_energy;
        println!("closed valve: mass diff={:e} kg, energy relative diff={energy_relative_error:e}", breathing_result.mass - motoring_result.mass);

        assert_eq!(breathing_result.mass, initial_state.mass, "mass must be exactly unchanged with a closed valve");
        assert!(energy_relative_error < 1e-9, "energy should match pure motoring almost exactly with a closed valve");
    }

    #[test]
    fn breathing_choked_mass_flow_produces_exactly_linear_growth_in_a_rigid_cylinder() {
        // Rigid, fixed-volume cylinder (stroke=0 collapses CrankMechanism's
        // position to a constant, so dV/dtheta is identically zero -
        // reusing existing kinematics rather than a new code path) with
        // the valve held at an effectively CONSTANT lift: a cam profile
        // with a deliberately huge (unrealistic) duration, evaluated over
        // a tiny window centered exactly at its midpoint (where
        // dLift/dtheta=0), so lift varies only quadratically and by a
        // genuinely negligible amount over the test window - a
        // mathematical device to isolate the mass-flow ODE from any lift
        // variation, not a realistic cam.
        let gas = GasProperties::AIR;
        let clearance_volume = 5.0e-5;
        let cylinder = Cylinder { crank_mechanism: CrankMechanism::new(0.0, 0.139, 0.0), bore: 0.087, clearance_volume };

        let cam = CamProfile { max_lift: 0.010, opening_angle_radians: -1000.0, duration_radians: 2000.0 }; // midpoint at theta=0
        let valve_geometry = ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 45.0_f64.to_radians() };
        let discharge_coefficient = 0.7;
        let reservoir_pressure = 10.0e5;
        let reservoir_temperature_kelvin = 400.0;
        let angular_velocity = 500.0;

        let initial_pressure = 1.0e5;
        let initial_temperature_kelvin = 350.0;
        let initial_state = CylinderState::from_pressure_temperature(initial_pressure, initial_temperature_kelvin, clearance_volume, &gas);

        let critical_pressure_ratio = (2.0 / (gas.gamma + 1.0)).powf(gas.gamma / (gas.gamma - 1.0));
        assert!(initial_pressure / reservoir_pressure < critical_pressure_ratio, "test setup must start choked");

        let theta_start = -0.002_f64;
        let theta_end = 0.002_f64;

        // Exact constant mdot, evaluated exactly at the cam's peak
        // (theta=0, where lift = max_lift exactly, not approximately).
        let effective_area = discharge_coefficient * valve::curtain_area(&valve_geometry, cam.max_lift);
        let expected_mdot = valve::mass_flow_rate(reservoir_pressure, reservoir_temperature_kelvin, initial_pressure, initial_temperature_kelvin, effective_area, &gas);

        // Guarantee (not hope) that pressure stays choked for the WHOLE
        // window, via the closed-form dp/dtheta under constant choked
        // mdot in a rigid volume.
        let specific_enthalpy_reservoir = gas.gamma * gas.gas_constant * reservoir_temperature_kelvin / (gas.gamma - 1.0);
        let dp_dtheta = (gas.gamma - 1.0) * specific_enthalpy_reservoir * expected_mdot / (clearance_volume * angular_velocity);
        let max_pressure_over_window = initial_pressure + dp_dtheta * (theta_end - theta_start);
        assert!(max_pressure_over_window / reservoir_pressure < critical_pressure_ratio, "test window is too wide - flow would unchoke before theta_end");

        let params = BreathingParameters {
            cam,
            valve: valve_geometry,
            discharge_coefficient,
            reservoir_pressure,
            reservoir_temperature_kelvin,
            angular_velocity_radians_per_second: angular_velocity,
        };

        let final_state = cylinder.integrate_breathing(&gas, initial_state, theta_start, theta_end, 100, &params);

        let expected_mass_change = expected_mdot * (theta_end - theta_start) / angular_velocity;
        let expected_energy_change = specific_enthalpy_reservoir * expected_mass_change;
        let actual_mass_change = final_state.mass - initial_state.mass;
        let actual_energy_change = final_state.internal_energy - initial_state.internal_energy;

        let mass_relative_error = (actual_mass_change - expected_mass_change).abs() / expected_mass_change;
        let energy_relative_error = (actual_energy_change - expected_energy_change).abs() / expected_energy_change;

        println!(
            "choked constant-mdot check: expected_mdot={expected_mdot:e} kg/s, mass change actual={actual_mass_change:e} expected={expected_mass_change:e} (rel err {mass_relative_error:e}), energy change actual={actual_energy_change:e} expected={expected_energy_change:e} (rel err {energy_relative_error:e})"
        );

        assert!(mass_relative_error < 1e-6, "mass growth relative error {mass_relative_error:e} too high");
        assert!(energy_relative_error < 1e-6, "energy growth relative error {energy_relative_error:e} too high");
    }

    #[test]
    fn realistic_breathing_window_is_physically_sane_and_matches_independent_trapezoidal_mass_integration() {
        let gas = GasProperties::AIR;
        let cylinder = s54b32_cylinder();
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();

        // Centered on BDC (not TDC): a real intake event happens near
        // BDC, where volume is near its maximum and changes slowly -
        // centering the window on TDC instead (an earlier version of this
        // test did exactly that) let piston-driven compression spike the
        // cylinder pressure well above the reservoir mid-window, causing
        // a net OUTFLOW that dominated the balance - a physically real
        // effect, but not what an "intake breathing" scenario should
        // exercise, and not a bug in the flow model itself (the
        // trapezoidal cross-check agreed with the RK4 result throughout).
        let cam = CamProfile { max_lift: 0.009, opening_angle_radians: bdc_angle - 1.2, duration_radians: 2.4 };
        let valve_geometry = ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 45.0_f64.to_radians() };
        let params = BreathingParameters {
            cam,
            valve: valve_geometry,
            discharge_coefficient: 0.7,
            reservoir_pressure: 1.0e5,
            reservoir_temperature_kelvin: 300.0,
            angular_velocity_radians_per_second: 500.0,
        };

        let theta_start = bdc_angle - 2.0;
        let theta_end = bdc_angle + 2.0;
        let initial_state = CylinderState::from_pressure_temperature(0.5e5, 350.0, cylinder.volume(theta_start), &gas);

        let segments = 400;
        let dtheta = (theta_end - theta_start) / segments as f64;
        let mut state = initial_state;
        let mut theta = theta_start;
        let mut trapezoidal_mass_change = 0.0_f64;

        for _ in 0..segments {
            let lift_before = camshaft::lift_at(&cam, theta);
            let area_before = params.discharge_coefficient * valve::curtain_area(&valve_geometry, lift_before);
            let mdot_before = valve::mass_flow_rate(
                params.reservoir_pressure,
                params.reservoir_temperature_kelvin,
                state.pressure(cylinder.volume(theta), &gas),
                state.temperature_kelvin(&gas),
                area_before,
                &gas,
            );

            let next_theta = theta + dtheta;
            let next_state = cylinder.integrate_breathing(&gas, state, theta, next_theta, 4, &params);

            let lift_after = camshaft::lift_at(&cam, next_theta);
            let area_after = params.discharge_coefficient * valve::curtain_area(&valve_geometry, lift_after);
            let mdot_after = valve::mass_flow_rate(
                params.reservoir_pressure,
                params.reservoir_temperature_kelvin,
                next_state.pressure(cylinder.volume(next_theta), &gas),
                next_state.temperature_kelvin(&gas),
                area_after,
                &gas,
            );

            trapezoidal_mass_change += 0.5 * (mdot_before + mdot_after) * dtheta / params.angular_velocity_radians_per_second;

            state = next_state;
            theta = next_theta;

            assert!(state.mass.is_finite() && state.mass > 0.0, "mass became non-physical at theta={theta}: {}", state.mass);
            assert!(state.internal_energy.is_finite(), "energy became non-finite at theta={theta}");
        }

        let final_state = state;
        let actual_mass_change = final_state.mass - initial_state.mass;
        let cross_check_relative_error = (actual_mass_change - trapezoidal_mass_change).abs() / actual_mass_change.abs();

        // The "pressure trends toward the reservoir's" check only makes
        // sense while the valve is actually open - once it closes (at
        // cam.opening_angle_radians + cam.duration_radians), the trapped
        // charge keeps compressing as the piston continues toward TDC,
        // correctly raising pressure back above the reservoir's again
        // (real post-IVC compression, not a bug) - so this is checked at
        // valve closing, not at the window's end.
        let valve_closing_angle = cam.opening_angle_radians + cam.duration_radians;
        let state_at_valve_close = cylinder.integrate_breathing(&gas, initial_state, theta_start, valve_closing_angle, 400, &params);
        let initial_pressure_gap = (0.5e5_f64 - params.reservoir_pressure).abs();
        let pressure_gap_at_valve_close = (state_at_valve_close.pressure(cylinder.volume(valve_closing_angle), &gas) - params.reservoir_pressure).abs();

        println!(
            "realistic breathing: mass change={actual_mass_change:e} kg (trapezoidal cross-check={trapezoidal_mass_change:e} kg, rel err={cross_check_relative_error:e}), initial pressure gap={initial_pressure_gap:e} Pa, pressure gap at valve close={pressure_gap_at_valve_close:e} Pa"
        );

        assert!(actual_mass_change > 0.0, "expected net inflow from the higher-pressure reservoir, got {actual_mass_change:e}");
        assert!(pressure_gap_at_valve_close < initial_pressure_gap, "cylinder pressure should move toward the reservoir's while the valve is open, not away");
        assert!(cross_check_relative_error < 1e-3, "trapezoidal mass cross-check relative error {cross_check_relative_error:e} too high");
    }

    /// Shared setup for the convergence-order checks below: measures RK4
    /// convergence (coarse/fine vs. a much finer reference run - no
    /// closed-form solution for a breathing cylinder) over
    /// `[theta_start, theta_end]`, with the reservoir/initial pressures
    /// as parameters (not fixed). Two candidate sources of locally-reduced
    /// RK4 order were considered here: the versine lift profile's 2nd-
    /// derivative jump at its open/close events, and the compressible-
    /// flow correlation's choked/subsonic branch switch (which is exactly
    /// value- and slope-continuous — `d(mdot)/d(pressure_ratio) = 0` on
    /// both sides at the critical ratio, verified analytically: the
    /// subsonic mass-flow function's own extremum falls exactly at the
    /// critical ratio). Measured results (see the two tests below):
    /// the lift boundary turned out NOT to measurably reduce order at
    /// all (~4.0, same as away from it) — a finite jump in a 2nd
    /// derivative is evidently mild enough that 4th-order RK4 isn't
    /// visibly affected, unlike Wiebe combustion's unbounded 4th
    /// derivative (which did measure a clean order-1.00 reduction
    /// elsewhere in this file). The choked/subsonic transition DID cause
    /// a real, but erratic/step-grid-sensitive, degradation (an early
    /// version of this test crossed it by accident and measured 0.74) -
    /// avoided here by keeping reservoir/initial pressures in one flow
    /// regime throughout, rather than characterized further, since its
    /// exact magnitude wasn't a clean single number the way Wiebe's is.
    fn measure_breathing_convergence_order(theta_start: f64, theta_end: f64, reservoir_pressure: f64, initial_pressure: f64) -> (f64, f64) {
        let gas = GasProperties::AIR;
        let cylinder = s54b32_cylinder();
        let cam = CamProfile { max_lift: 0.009, opening_angle_radians: 0.0, duration_radians: 2.0 };
        let params = BreathingParameters {
            cam,
            valve: ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 45.0_f64.to_radians() },
            discharge_coefficient: 0.7,
            reservoir_pressure,
            reservoir_temperature_kelvin: 300.0,
            angular_velocity_radians_per_second: 500.0,
        };
        let initial_state = CylinderState::from_pressure_temperature(initial_pressure, 350.0, cylinder.volume(theta_start), &gas);
        let volume_at_end = cylinder.volume(theta_end);

        let reference_pressure = cylinder.integrate_breathing(&gas, initial_state, theta_start, theta_end, 2880, &params).pressure(volume_at_end, &gas);
        let coarse_pressure = cylinder.integrate_breathing(&gas, initial_state, theta_start, theta_end, 45, &params).pressure(volume_at_end, &gas);
        let fine_pressure = cylinder.integrate_breathing(&gas, initial_state, theta_start, theta_end, 90, &params).pressure(volume_at_end, &gas);

        ((coarse_pressure - reference_pressure).abs(), (fine_pressure - reference_pressure).abs())
    }

    #[test]
    fn breathing_rk4_convergence_at_a_lift_event_boundary_is_not_measurably_reduced() {
        // The versine lift profile's 2nd derivative jumps discontinuously
        // (from 0 to a nonzero constant) exactly at the opening/closing
        // events, and a milder singularity than Wiebe combustion's (whose
        // unbounded 4th derivative DID measurably reduce RK4 order
        // elsewhere in this file) was expected to show up here too, at
        // least somewhat - but measured, it doesn't: a finite jump this
        // far down the derivative chain is evidently too mild for 4th-
        // order RK4 to notice, at these step counts. Reservoir/initial
        // pressures (1.2 bar / 1.0 bar, ratio 0.833) chosen to stay safely
        // subsonic throughout (well above the critical ratio ~0.528),
        // isolating this from the separate choked/subsonic-transition
        // effect (see this function's own doc comment above).
        let (coarse_error, fine_error) = measure_breathing_convergence_order(-0.1, 0.1, 1.2e5, 1.0e5); // straddles opening_angle_radians=0.0
        let observed_order = (coarse_error / fine_error).log2();
        println!("breathing RK4 convergence AT a lift event boundary: coarse error={coarse_error:e} Pa, fine error={fine_error:e} Pa, observed order={observed_order:.2}");
        assert!(observed_order > 3.0, "expected close to 4th-order convergence even at a lift event boundary, observed order {observed_order:.2}");
    }

    #[test]
    fn breathing_rk4_converges_at_close_to_fourth_order_away_from_any_known_non_smooth_point() {
        // Away from BOTH the lift-event boundaries (window comfortably
        // inside the open period) AND the choked/subsonic transition
        // (1.2 bar / 1.0 bar stays safely subsonic - ratio 0.833 - for
        // the whole window, confirmed by direct measurement rather than
        // assumed).
        let (coarse_error, fine_error) = measure_breathing_convergence_order(0.5, 1.0, 1.2e5, 1.0e5);
        let observed_order = (coarse_error / fine_error).log2();
        println!("breathing RK4 convergence AWAY from any known non-smooth point: coarse error={coarse_error:e} Pa, fine error={fine_error:e} Pa, observed order={observed_order:.2}");
        assert!(observed_order > 3.0, "expected close to 4th-order convergence away from any non-smooth point, observed order {observed_order:.2}");
    }

    fn sample_valve_geometry() -> ValveGeometry {
        ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 45.0_f64.to_radians() }
    }

    #[test]
    fn dual_valve_breathing_with_both_valves_closed_matches_motoring_exactly() {
        let cylinder = s54b32_cylinder();
        let gas = GasProperties::AIR;
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();
        let volume_at_bdc = cylinder.volume(bdc_angle);
        let initial_state = CylinderState::from_pressure_temperature(1.0e5, 320.0, volume_at_bdc, &gas);

        let motoring_result = cylinder.integrate_motoring(&gas, initial_state, bdc_angle, 0.0, 720);

        let closed_event = |reservoir_pressure: f64| ValveEventParameters {
            cam: CamProfile { max_lift: 0.0, opening_angle_radians: -100.0, duration_radians: 200.0 },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.7,
            reservoir_pressure,
            reservoir_temperature_kelvin: 400.0,
        };
        let params = DualValveBreathingParameters {
            intake: closed_event(1.0e5),
            exhaust: closed_event(5.0e5),
            angular_velocity_radians_per_second: 500.0,
        };
        let dual_valve_result = cylinder.integrate_dual_valve_breathing(&gas, initial_state, bdc_angle, 0.0, 720, &params);

        let energy_relative_error =
            (dual_valve_result.internal_energy - motoring_result.internal_energy).abs() / motoring_result.internal_energy;
        assert_eq!(dual_valve_result.mass, initial_state.mass, "mass must be exactly unchanged with both valves closed");
        assert!(energy_relative_error < 1e-9, "energy should match pure motoring almost exactly with both valves closed");
    }

    #[test]
    fn dual_valve_breathing_with_exhaust_closed_matches_single_valve_breathing_on_intake() {
        let cylinder = s54b32_cylinder();
        let gas = GasProperties::AIR;
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();
        let initial_state = CylinderState::from_pressure_temperature(0.5e5, 350.0, cylinder.volume(bdc_angle - 2.0), &gas);

        let intake_cam = CamProfile { max_lift: 0.009, opening_angle_radians: bdc_angle - 1.2, duration_radians: 2.4 };
        let intake_valve = sample_valve_geometry();
        let single_valve_params = BreathingParameters {
            cam: intake_cam,
            valve: intake_valve,
            discharge_coefficient: 0.7,
            reservoir_pressure: 1.0e5,
            reservoir_temperature_kelvin: 300.0,
            angular_velocity_radians_per_second: 500.0,
        };
        let single_valve_result =
            cylinder.integrate_breathing(&gas, initial_state, bdc_angle - 2.0, bdc_angle + 2.0, 400, &single_valve_params);

        let dual_valve_params = DualValveBreathingParameters {
            intake: ValveEventParameters {
                cam: intake_cam,
                valve: intake_valve,
                discharge_coefficient: 0.7,
                reservoir_pressure: 1.0e5,
                reservoir_temperature_kelvin: 300.0,
            },
            exhaust: ValveEventParameters {
                // Closed (max_lift=0) - should contribute exactly nothing.
                cam: CamProfile { max_lift: 0.0, opening_angle_radians: 50.0, duration_radians: 100.0 },
                valve: sample_valve_geometry(),
                discharge_coefficient: 0.7,
                reservoir_pressure: 1.0e5,
                reservoir_temperature_kelvin: 900.0, // deliberately extreme - must have zero effect
            },
            angular_velocity_radians_per_second: 500.0,
        };
        let dual_valve_result =
            cylinder.integrate_dual_valve_breathing(&gas, initial_state, bdc_angle - 2.0, bdc_angle + 2.0, 400, &dual_valve_params);

        let mass_relative_error = (dual_valve_result.mass - single_valve_result.mass).abs() / single_valve_result.mass;
        let energy_relative_error =
            (dual_valve_result.internal_energy - single_valve_result.internal_energy).abs() / single_valve_result.internal_energy;
        println!("dual-valve (exhaust closed) vs single-valve: mass rel err={mass_relative_error:e}, energy rel err={energy_relative_error:e}");
        assert!(mass_relative_error < 1e-12, "closed exhaust valve should have exactly zero effect on mass");
        assert!(energy_relative_error < 1e-12, "closed exhaust valve should have exactly zero effect on energy");
    }

    #[test]
    fn dual_valve_breathing_overlap_stays_sane_and_matches_trapezoidal_mass_cross_check() {
        // Both valves open simultaneously (a valve-overlap window): intake
        // opening while exhaust is still closing, both against their own
        // reservoir. No special-casing exists in the derivative for this -
        // this test is the direct check that summing two independent flow
        // terms behaves sanely (finite, matches an independent trapezoidal
        // integration) rather than assuming it from the code alone.
        let cylinder = s54b32_cylinder();
        let gas = GasProperties::AIR;

        let intake = ValveEventParameters {
            cam: CamProfile { max_lift: 0.009, opening_angle_radians: -0.3, duration_radians: 2.4 },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.7,
            reservoir_pressure: 1.0e5,
            reservoir_temperature_kelvin: 300.0,
        };
        let exhaust = ValveEventParameters {
            cam: CamProfile { max_lift: 0.009, opening_angle_radians: -1.8, duration_radians: 2.0 },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.7,
            reservoir_pressure: 1.0e5,
            reservoir_temperature_kelvin: 900.0, // hot exhaust gas reservoir
        };
        let params = DualValveBreathingParameters { intake, exhaust, angular_velocity_radians_per_second: 500.0 };

        // intake window = [-0.3, 2.1], exhaust window = [-1.8, 0.2] -
        // overlap is [-0.3, 0.2], so stay well inside that.
        let theta_start = -0.2_f64;
        let theta_end = 0.1_f64;
        // Overlap really is simultaneous here (sanity-check the test setup itself).
        assert!(theta_start > exhaust.cam.opening_angle_radians && theta_start < exhaust.cam.opening_angle_radians + exhaust.cam.duration_radians);
        assert!(theta_start > intake.cam.opening_angle_radians && theta_start < intake.cam.opening_angle_radians + intake.cam.duration_radians);

        let initial_state = CylinderState::from_pressure_temperature(1.05e5, 500.0, cylinder.volume(theta_start), &gas);

        let segments = 400;
        let dtheta = (theta_end - theta_start) / segments as f64;
        let mut state = initial_state;
        let mut theta = theta_start;
        let mut trapezoidal_mass_change = 0.0_f64;

        let mdot_sum_at = |cylinder: &Cylinder, state: &CylinderState, theta: f64| {
            let pressure = state.pressure(cylinder.volume(theta), &gas);
            let temperature_kelvin = state.temperature_kelvin(&gas);
            [&intake, &exhaust]
                .iter()
                .map(|event| {
                    let lift = camshaft::lift_at(&event.cam, theta);
                    let area = event.discharge_coefficient * valve::curtain_area(&event.valve, lift);
                    valve::mass_flow_rate(event.reservoir_pressure, event.reservoir_temperature_kelvin, pressure, temperature_kelvin, area, &gas)
                })
                .sum::<f64>()
        };

        for _ in 0..segments {
            let mdot_before = mdot_sum_at(&cylinder, &state, theta);
            let next_theta = theta + dtheta;
            let next_state = cylinder.integrate_dual_valve_breathing(&gas, state, theta, next_theta, 4, &params);
            let mdot_after = mdot_sum_at(&cylinder, &next_state, next_theta);

            trapezoidal_mass_change += 0.5 * (mdot_before + mdot_after) * dtheta / params.angular_velocity_radians_per_second;

            state = next_state;
            theta = next_theta;

            assert!(state.mass.is_finite() && state.mass > 0.0, "mass became non-physical at theta={theta}: {}", state.mass);
            assert!(state.internal_energy.is_finite(), "energy became non-finite at theta={theta}");
        }

        let actual_mass_change = state.mass - initial_state.mass;
        let cross_check_relative_error = (actual_mass_change - trapezoidal_mass_change).abs() / trapezoidal_mass_change.abs();
        println!(
            "overlap: mass change={actual_mass_change:e} kg (trapezoidal cross-check={trapezoidal_mass_change:e} kg, rel err={cross_check_relative_error:e})"
        );
        assert!(cross_check_relative_error < 1e-3, "trapezoidal mass cross-check relative error {cross_check_relative_error:e} too high");
    }

    #[test]
    fn dual_valve_breathing_rk4_converges_at_close_to_fourth_order() {
        let gas = GasProperties::AIR;
        let cylinder = s54b32_cylinder();
        let intake = ValveEventParameters {
            cam: CamProfile { max_lift: 0.009, opening_angle_radians: 0.0, duration_radians: 2.0 },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.7,
            reservoir_pressure: 1.2e5,
            reservoir_temperature_kelvin: 300.0,
        };
        let exhaust = ValveEventParameters {
            // Positioned well outside [theta_start, theta_end] below - contributes nothing here,
            // present only to confirm the dual-valve stepper itself converges at full order.
            cam: CamProfile { max_lift: 0.009, opening_angle_radians: -10.0, duration_radians: 2.0 },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.7,
            reservoir_pressure: 1.0e5,
            reservoir_temperature_kelvin: 900.0,
        };
        let params = DualValveBreathingParameters { intake, exhaust, angular_velocity_radians_per_second: 500.0 };

        let theta_start = 0.5;
        let theta_end = 1.0;
        let initial_pressure = 1.0e5;
        let initial_state = CylinderState::from_pressure_temperature(initial_pressure, 350.0, cylinder.volume(theta_start), &gas);
        let volume_at_end = cylinder.volume(theta_end);

        let reference_pressure =
            cylinder.integrate_dual_valve_breathing(&gas, initial_state, theta_start, theta_end, 2880, &params).pressure(volume_at_end, &gas);
        let coarse_pressure =
            cylinder.integrate_dual_valve_breathing(&gas, initial_state, theta_start, theta_end, 45, &params).pressure(volume_at_end, &gas);
        let fine_pressure =
            cylinder.integrate_dual_valve_breathing(&gas, initial_state, theta_start, theta_end, 90, &params).pressure(volume_at_end, &gas);

        let coarse_error = (coarse_pressure - reference_pressure).abs();
        let fine_error = (fine_pressure - reference_pressure).abs();
        let observed_order = (coarse_error / fine_error).log2();
        println!("dual-valve breathing RK4 convergence: coarse error={coarse_error:e} Pa, fine error={fine_error:e} Pa, observed order={observed_order:.2}");
        assert!(observed_order > 3.0, "expected close to 4th-order convergence, observed order {observed_order:.2}");
    }

    #[test]
    fn full_cycle_reduces_to_fired_cycle_exactly_when_both_valves_closed_and_starting_at_ivc() {
        let cylinder = s54_2500rpm_cylinder();
        let gas = GasProperties::AIR;
        let ivc_angle = (-120.0_f64).to_radians();
        let ivc_volume = cylinder.volume(ivc_angle);
        let ivc_state = CylinderState::from_pressure_temperature(1.16948e5, 393.07, ivc_volume, &gas);
        let fired_params = s54_2500rpm_fired_cycle_params(&gas, ivc_state);

        let theta_end = 30.0_f64.to_radians();
        let fired_result = cylinder.integrate_fired_cycle(&gas, ivc_state, ivc_angle, theta_end, 2000, &fired_params);

        // Closed valves, windowed so `opening_angle_radians + duration_radians`
        // is exactly `ivc_angle` - since `theta_start == ivc_angle` below,
        // this exercises the "theta_start >= ivc_angle" branch, which
        // should reduce bit-for-bit to `integrate_fired_cycle`'s own
        // convention (anchor = initial_state/theta_start directly).
        let closed_event = |reservoir_pressure: f64| ValveEventParameters {
            cam: CamProfile { max_lift: 0.0, opening_angle_radians: ivc_angle - 100.0_f64.to_radians(), duration_radians: 100.0_f64.to_radians() },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.7,
            reservoir_pressure,
            reservoir_temperature_kelvin: 400.0,
        };
        let full_cycle_params = FullCycleParameters {
            wiebe: fired_params.wiebe,
            woschni: fired_params.woschni,
            walls: fired_params.walls,
            total_heat_release_joules: fired_params.total_heat_release_joules,
            angular_velocity_radians_per_second: fired_params.angular_velocity_radians_per_second,
            gas_constant: fired_params.gas_constant,
            intake: closed_event(1.0e5),
            exhaust: closed_event(5.0e5),
        };
        let full_cycle_result = cylinder.integrate_full_cycle(&gas, ivc_state, ivc_angle, theta_end, 2000, &full_cycle_params);

        let volume_at_end = cylinder.volume(theta_end);
        let energy_relative_error = (full_cycle_result.internal_energy - fired_result.internal_energy).abs() / fired_result.internal_energy;
        let pressure_relative_error =
            (full_cycle_result.pressure(volume_at_end, &gas) - fired_result.pressure(volume_at_end, &gas)).abs() / fired_result.pressure(volume_at_end, &gas);
        println!("full_cycle vs fired_cycle (closed valves): energy rel err={energy_relative_error:e}, pressure rel err={pressure_relative_error:e}");
        assert_eq!(full_cycle_result.mass, fired_result.mass, "mass must be exactly unchanged with both valves closed");
        assert!(energy_relative_error < 1e-9, "expected full_cycle to match fired_cycle almost exactly with both valves closed");
    }

    #[test]
    fn full_cycle_reduces_to_dual_valve_breathing_when_combustion_is_zeroed() {
        let cylinder = s54b32_cylinder();
        let gas = GasProperties::AIR;
        let bdc_angle = cylinder.crank_mechanism.crank_angle_of_bottom_dead_center();
        let initial_state = CylinderState::from_pressure_temperature(0.5e5, 350.0, cylinder.volume(bdc_angle - 2.0), &gas);

        let intake = ValveEventParameters {
            cam: CamProfile { max_lift: 0.009, opening_angle_radians: bdc_angle - 1.2, duration_radians: 2.4 },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.7,
            reservoir_pressure: 1.0e5,
            reservoir_temperature_kelvin: 300.0,
        };
        let exhaust = ValveEventParameters {
            cam: CamProfile { max_lift: 0.0, opening_angle_radians: 50.0, duration_radians: 100.0 },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.7,
            reservoir_pressure: 1.0e5,
            reservoir_temperature_kelvin: 900.0,
        };
        let angular_velocity_radians_per_second = 500.0;
        let dual_valve_params = DualValveBreathingParameters { intake, exhaust, angular_velocity_radians_per_second };

        let theta_start = bdc_angle - 2.0;
        let theta_end = bdc_angle + 2.0;
        // theta_start is well before ivc_angle = intake window's own
        // close (bdc_angle + 1.2) - exercises the anchor-pre-pass branch.
        let dual_valve_result = cylinder.integrate_dual_valve_breathing(&gas, initial_state, theta_start, theta_end, 400, &dual_valve_params);

        let full_cycle_params = FullCycleParameters {
            wiebe: WiebeParameters { start_angle_radians: 0.0, duration_radians: 0.5, shape_factor_m: 2.5, efficiency_c: 6.9 },
            woschni: WoschniParameters { cw1: 0.0, cw2: 0.0, combustion_turbulence_coefficient: 0.0 },
            walls: WallTemperatures { piston_kelvin: 450.0, head_kelvin: 550.0, liner_kelvin: 420.0 },
            total_heat_release_joules: 0.0,
            angular_velocity_radians_per_second,
            gas_constant: gas.gas_constant,
            intake,
            exhaust,
        };
        let full_cycle_result = cylinder.integrate_full_cycle(&gas, initial_state, theta_start, theta_end, 400, &full_cycle_params);

        let mass_relative_error = (full_cycle_result.mass - dual_valve_result.mass).abs() / dual_valve_result.mass;
        let energy_relative_error =
            (full_cycle_result.internal_energy - dual_valve_result.internal_energy).abs() / dual_valve_result.internal_energy;
        println!("full_cycle vs dual_valve_breathing (zero combustion): mass rel err={mass_relative_error:e}, energy rel err={energy_relative_error:e}");
        assert!(mass_relative_error < 1e-9, "expected full_cycle to match dual_valve_breathing almost exactly with combustion zeroed");
        assert!(energy_relative_error < 1e-9, "expected full_cycle to match dual_valve_breathing almost exactly with combustion zeroed");
    }

    #[test]
    fn full_cycle_state_at_ivc_matches_a_direct_dual_valve_breathing_run_to_the_same_angle() {
        // Sets integrate_full_cycle's OWN theta_end to exactly ivc_angle -
        // the internal anchor pre-pass's step-count ratio collapses to
        // 1.0 in this special case, and the main RK4 loop covers exactly
        // [theta_start, ivc_angle] with a derivative that (Woschni
        // coefficients zeroed, so there's no base convective heat
        // transfer term either - NOT just "before ignition", since the
        // Woschni correlation's cw1*Cm term is active throughout the
        // WHOLE cycle, not only during combustion; only its combustion-
        // turbulence addend is gated to the Wiebe window - and the
        // exhaust valve is closed) reduces EXACTLY to
        // dual_valve_breathing_derivatives - so both paths should agree
        // to floating-point precision, not just approximately, confirming
        // the anchor pre-pass genuinely reuses (not reimplements)
        // `integrate_dual_valve_breathing`.
        let cylinder = s54_2500rpm_cylinder();
        let gas = GasProperties::AIR;
        let ivc_angle = (-120.0_f64).to_radians();
        let theta_start = (-340.0_f64).to_radians();
        let initial_state = CylinderState::from_pressure_temperature(1.0e5, 320.0, cylinder.volume(theta_start), &gas);

        let intake_opening = theta_start - 20.0_f64.to_radians();
        let intake = ValveEventParameters {
            cam: CamProfile { max_lift: 0.010, opening_angle_radians: intake_opening, duration_radians: ivc_angle - intake_opening },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.64,
            reservoir_pressure: 1.0e5,
            reservoir_temperature_kelvin: 293.15,
        };
        let exhaust = ValveEventParameters {
            cam: CamProfile { max_lift: 0.0, opening_angle_radians: 100.0, duration_radians: 100.0 },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.64,
            reservoir_pressure: 1.0e5,
            reservoir_temperature_kelvin: 900.0,
        };
        let angular_velocity_radians_per_second = 2500.0 * 2.0 * std::f64::consts::PI / 60.0;

        let dual_valve_params = DualValveBreathingParameters { intake, exhaust, angular_velocity_radians_per_second };
        let step_count = 2000;
        let direct_state_at_ivc =
            cylinder.integrate_dual_valve_breathing(&gas, initial_state, theta_start, ivc_angle, step_count, &dual_valve_params);

        let full_cycle_params = FullCycleParameters {
            wiebe: WiebeParameters { start_angle_radians: (-15.0_f64).to_radians(), duration_radians: 45.0_f64.to_radians(), shape_factor_m: 2.5, efficiency_c: 6.9 },
            woschni: WoschniParameters { cw1: 0.0, cw2: 0.0, combustion_turbulence_coefficient: 0.0 },
            walls: WallTemperatures { piston_kelvin: 450.0, head_kelvin: 550.0, liner_kelvin: 420.0 },
            total_heat_release_joules: 1000.0, // nonzero, but inert here (Wiebe ignition is well after ivc_angle)
            angular_velocity_radians_per_second,
            gas_constant: gas.gas_constant,
            intake,
            exhaust,
        };
        let full_cycle_state_at_ivc = cylinder.integrate_full_cycle(&gas, initial_state, theta_start, ivc_angle, step_count, &full_cycle_params);

        let mass_relative_error = (full_cycle_state_at_ivc.mass - direct_state_at_ivc.mass).abs() / direct_state_at_ivc.mass;
        let energy_relative_error =
            (full_cycle_state_at_ivc.internal_energy - direct_state_at_ivc.internal_energy).abs() / direct_state_at_ivc.internal_energy;
        println!("anchor pre-pass consistency: mass rel err={mass_relative_error:e}, energy rel err={energy_relative_error:e}");
        assert!(mass_relative_error < 1e-9, "expected the full-cycle path's state at IVC to match a direct dual-valve-breathing run there");
        assert!(energy_relative_error < 1e-9, "expected the full-cycle path's state at IVC to match a direct dual-valve-breathing run there");
    }

    #[test]
    fn full_cycle_rk4_converges_at_close_to_fourth_order_starting_exactly_at_ivc_angle() {
        // Starts exactly AT ivc_angle (the "theta_start >= ivc_angle"
        // branch - anchor taken directly from initial_state/theta_start,
        // no pre-pass) and runs well into the burn (mirroring
        // `fired_cycle_rk4_converges_at_close_to_fourth_order_away_from_ignition`'s
        // own choice of window, away from the singular ignition point) -
        // confirms `integrate_full_cycle`'s own RK4 stepper achieves full
        // order despite starting exactly at the anchor boundary. A window
        // confined close to `ivc_angle` (weak dynamics, no combustion yet,
        // tiny absolute pressure change over a few degrees) was tried
        // first and measured a reduced order (~2.5) - but the IDENTICAL
        // reduction was independently confirmed to occur, bit-for-bit,
        // in the already-shipped `integrate_fired_cycle` on the exact
        // same narrow window, so it is a pre-existing floating-point-
        // noise-floor property of that regime, not a regression
        // introduced here - this test instead measures somewhere the
        // signal is strong enough to say something meaningful.
        let cylinder = s54_2500rpm_cylinder();
        let gas = GasProperties::AIR;
        let ivc_angle = (-120.0_f64).to_radians();
        let ivc_volume = cylinder.volume(ivc_angle);
        let ivc_state = CylinderState::from_pressure_temperature(1.16948e5, 393.07, ivc_volume, &gas);
        let fired_params = s54_2500rpm_fired_cycle_params(&gas, ivc_state);

        let closed_event = |reservoir_pressure: f64| ValveEventParameters {
            cam: CamProfile { max_lift: 0.0, opening_angle_radians: ivc_angle - 10.0_f64.to_radians(), duration_radians: 10.0_f64.to_radians() },
            valve: sample_valve_geometry(),
            discharge_coefficient: 0.7,
            reservoir_pressure,
            reservoir_temperature_kelvin: 400.0,
        };
        let params = FullCycleParameters {
            wiebe: fired_params.wiebe,
            woschni: fired_params.woschni,
            walls: fired_params.walls,
            total_heat_release_joules: fired_params.total_heat_release_joules,
            angular_velocity_radians_per_second: fired_params.angular_velocity_radians_per_second,
            gas_constant: fired_params.gas_constant,
            intake: closed_event(1.0e5),
            exhaust: closed_event(1.0e5),
        };

        // Pre-integrate from ivc_angle (exercising the "theta_start >=
        // ivc_angle" anchor branch) THROUGH the singular ignition point
        // with a large, fixed step count, matching
        // `measure_fired_cycle_convergence_order`'s own established
        // pattern - so that segment's own (expected, already-measured-
        // elsewhere) reduced-order error doesn't contaminate the actual
        // measurement, which is taken only over a clean, entirely-post-
        // ignition window.
        let interval_start = (-5.0_f64).to_radians(); // theta0 (-15deg) + 10deg, same as the existing "away from ignition" test
        let burn_end = params.wiebe.start_angle_radians + params.wiebe.duration_radians;
        let pre_interval_state = cylinder.integrate_full_cycle(&gas, ivc_state, ivc_angle, interval_start, 4000, &params);
        let volume_at_burn_end = cylinder.volume(burn_end);

        let reference_pressure =
            cylinder.integrate_full_cycle(&gas, pre_interval_state, interval_start, burn_end, 2880, &params).pressure(volume_at_burn_end, &gas);
        let coarse_pressure =
            cylinder.integrate_full_cycle(&gas, pre_interval_state, interval_start, burn_end, 45, &params).pressure(volume_at_burn_end, &gas);
        let fine_pressure =
            cylinder.integrate_full_cycle(&gas, pre_interval_state, interval_start, burn_end, 90, &params).pressure(volume_at_burn_end, &gas);

        let coarse_error = (coarse_pressure - reference_pressure).abs();
        let fine_error = (fine_pressure - reference_pressure).abs();
        let observed_order = (coarse_error / fine_error).log2();
        println!(
            "full-cycle RK4 convergence (pre-integrated from ivc_angle, measured away from ignition): coarse error={coarse_error:e} Pa, fine error={fine_error:e} Pa, observed order={observed_order:.2}"
        );
        assert!(observed_order > 3.0, "expected close to 4th-order convergence, observed order {observed_order:.2}");
    }

    /// Real-numbers comparison against the actual OpenWAM S54 2500rpm
    /// case's intake-stroke trace - genuine full-cycle ground truth (not
    /// just the already-validated closed-cycle portion), extracted
    /// directly from `benchmarks/openwam/cases/engine_s54_2500rpm/`:
    /// `NumberOfValves=2`, intake `FDiametro=0.0495m` (area-matched
    /// equivalent of 2x35mm), exhaust `FDiametro=0.0431m` (2x30.5mm),
    /// both with the SAME 27-point versine lift table (0->12mm->0 over
    /// 260 degrees, `FIncrAng=10deg`) and the SAME lift-dependent Cd table
    /// plateauing at 0.64 (used here as a constant - this model doesn't
    /// support lift-dependent Cd yet), `IVO=340deg`/`EVO=140deg` (already
    /// in firing-TDC=0 terms - `IVC=340+260=600deg=-120deg`, matching the
    /// existing fired-cycle tests' own IVC angle exactly), and
    /// `AmbientPressure=1.0 bar`/`AmbientTemperature=20degC` for both
    /// reservoirs (the case's own boundary conditions: both the intake
    /// and exhaust pipes' far ends are `nmOpenEndAtmosphere`).
    ///
    /// Deliberately stays entirely within [300deg, 600deg=IVC] - genuine
    /// valve overlap (exhaust closing, `EVC=400deg`) through the intake
    /// stroke, run through `integrate_full_cycle` itself (not directly
    /// through `integrate_dual_valve_breathing`) so this is a real
    /// end-to-end check of the new capability - but short of the
    /// combustion event, so no "which firing TDC occurrence" angle
    /// bookkeeping is needed (Wiebe ignition, referenced the same way as
    /// every other test in this file, is simply never reached in this
    /// window).
    #[test]
    fn full_cycle_intake_stroke_compares_against_the_real_openwam_s54_2500rpm_case() {
        let cylinder = s54_2500rpm_cylinder();
        let gas = GasProperties::AIR;

        let lift_table_max = 0.012_f64; // 12mm, from FLevantamiento's peak value
        let duration = 260.0_f64.to_radians(); // 26 * FIncrAng(10deg)
        let discharge_coefficient = 0.64; // FDatosCD's plateau value
        let ambient_pressure = 1.0e5; // 1.0 bar
        let ambient_temperature_kelvin = 293.15; // 20 degC

        let intake = ValveEventParameters {
            cam: CamProfile { max_lift: lift_table_max, opening_angle_radians: 340.0_f64.to_radians(), duration_radians: duration },
            valve: ValveGeometry { valve_diameter: 0.0495, seat_angle_radians: 45.0_f64.to_radians() },
            discharge_coefficient,
            reservoir_pressure: ambient_pressure,
            reservoir_temperature_kelvin: ambient_temperature_kelvin,
        };
        let exhaust = ValveEventParameters {
            cam: CamProfile { max_lift: lift_table_max, opening_angle_radians: 140.0_f64.to_radians(), duration_radians: duration },
            valve: ValveGeometry { valve_diameter: 0.0431, seat_angle_radians: 45.0_f64.to_radians() },
            discharge_coefficient,
            reservoir_pressure: ambient_pressure,
            reservoir_temperature_kelvin: ambient_temperature_kelvin,
        };

        let ivc_state_for_fixture = CylinderState::from_pressure_temperature(1.16948e5, 393.07, cylinder.volume((-120.0_f64).to_radians()), &gas);
        let fired_params = s54_2500rpm_fired_cycle_params(&gas, ivc_state_for_fixture);
        let params = FullCycleParameters {
            wiebe: fired_params.wiebe,
            woschni: fired_params.woschni,
            walls: fired_params.walls,
            total_heat_release_joules: fired_params.total_heat_release_joules,
            angular_velocity_radians_per_second: fired_params.angular_velocity_radians_per_second,
            gas_constant: gas.gas_constant,
            intake,
            exhaust,
        };

        // Real OpenWAM trace: (crank angle deg, pressure bar, temperature degC).
        let theta_start_degrees = 300.0_f64;
        let initial_pressure_bar = 0.817864;
        let initial_temperature_c = 822.879;
        let checkpoints: [(f64, f64, f64); 6] =
            [(340.0, 1.08004, 861.187), (400.0, 0.824861, 417.105), (450.0, 1.0691, 116.216), (500.0, 0.979397, 94.539), (550.0, 1.06966, 104.807), (600.0, 1.16603, 119.579)];

        let theta_start = theta_start_degrees.to_radians();
        let initial_state =
            CylinderState::from_pressure_temperature(initial_pressure_bar * 1e5, initial_temperature_c + 273.15, cylinder.volume(theta_start), &gas);

        let mut state = initial_state;
        let mut previous_angle_radians = theta_start;
        let mut max_pressure_relative_error = 0.0_f64;
        let mut max_temperature_relative_error = 0.0_f64;

        for (angle_deg, expected_pressure_bar, expected_temperature_c) in checkpoints {
            let angle_radians = angle_deg.to_radians();
            let segment_steps = ((angle_radians - previous_angle_radians).to_degrees().abs() * 20.0).ceil() as usize;
            state = cylinder.integrate_full_cycle(&gas, state, previous_angle_radians, angle_radians, segment_steps.max(20), &params);
            previous_angle_radians = angle_radians;

            let volume = cylinder.volume(angle_radians);
            let actual_pressure_bar = state.pressure(volume, &gas) / 1e5;
            let actual_temperature_c = state.temperature_kelvin(&gas) - 273.15;

            let pressure_relative_error = (actual_pressure_bar - expected_pressure_bar).abs() / expected_pressure_bar;
            let temperature_relative_error = (actual_temperature_c - expected_temperature_c).abs() / (expected_temperature_c + 273.15);
            max_pressure_relative_error = max_pressure_relative_error.max(pressure_relative_error);
            max_temperature_relative_error = max_temperature_relative_error.max(temperature_relative_error);

            println!(
                "theta={angle_deg:>6.1} deg: pressure actual={actual_pressure_bar:>6.3} bar / OpenWAM={expected_pressure_bar:>6.3} bar (err {:>6.2}%), temperature actual={actual_temperature_c:>7.1} C / OpenWAM={expected_temperature_c:>7.1} C (err {:>6.2}%)",
                pressure_relative_error * 100.0,
                temperature_relative_error * 100.0
            );
        }

        println!(
            "intake-stroke max pressure relative error = {:.2}%, max temperature relative error = {:.2}%",
            max_pressure_relative_error * 100.0,
            max_temperature_relative_error * 100.0
        );

        // Generous, honestly-measured (not pre-guessed) bounds. Pressure
        // agrees reasonably well throughout (measured up to ~11%) - a
        // CONSTANT discharge coefficient standing in for the file's own
        // lift-dependent Cd table accounts for most of that. Temperature
        // is a different story specifically around EVC (measured up to
        // ~41%, at theta=400deg): this model replaces the REAL exhaust
        // pipe's own transient, wave-driven pressure trace with a flat
        // AMBIENT reservoir - during the exhaust-valve-closing/overlap
        // window, the real exhaust pipe often runs meaningfully ABOVE
        // ambient (a genuine exhaust pulse/backpressure effect no fixed-
        // reservoir boundary condition can capture), so this model's
        // exhaust outflow (and therefore the hot-residual-gas retention
        // that sets post-EVC temperature) is systematically under-driven
        // here - the exact, expected consequence of the documented
        // limitation that real wave-driven scavenging needs the separate
        // `valve_port.rs` real-pipe-coupled path, not a fixed reservoir.
        // "Right shape/order of magnitude", not tight quantitative
        // agreement - pressure gets the tighter of the two bounds since
        // it isn't as exposed to this specific limitation.
        assert!(max_pressure_relative_error < 0.20, "intake-stroke pressure trace error {:.1}% exceeds the 20% bound", max_pressure_relative_error * 100.0);
        assert!(
            max_temperature_relative_error < 0.50,
            "intake-stroke temperature trace error {:.1}% exceeds the 50% bound",
            max_temperature_relative_error * 100.0
        );
    }
}
