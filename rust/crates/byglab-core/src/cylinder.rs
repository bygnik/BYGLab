//! 0D lumped cylinder thermodynamics: cylinder volume from crank
//! kinematics (`crank_mechanism.rs`) plus a first-law energy balance for
//! the trapped charge.
//!
//! This pass implements *motoring* only — no combustion, no valve mass
//! flow, no wall heat transfer, just the pure compression/expansion of a
//! closed, adiabatic charge. Deliberately scoped this way: motoring has an
//! exact analytical solution (the isentropic relation `p*V^gamma = const`)
//! to validate the energy-balance integration against, *before* combustion
//! heat release and valve mass flow — each with their own, much harder to
//! isolate, sources of error — get added on top of it.

use crate::crank_mechanism::CrankMechanism;
use crate::gas::GasProperties;

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
    /// `dU = -p dV`.
    fn motoring_energy_derivative(&self, gas: &GasProperties, state: &CylinderState, crank_angle_from_tdc_radians: f64) -> f64 {
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
}
