//! HLLC approximate Riemann solver: computes the numerical flux across the
//! interface between two different gas states (the "Riemann problem" that
//! arises at every cell face in a finite-volume scheme).
//!
//! HLLC extends the simpler HLL solver with an explicit middle ("star")
//! wave, which is what lets it resolve a contact discontinuity (a jump in
//! density/temperature with *no* jump in pressure or velocity) instead of
//! smearing it into the surrounding waves. That matters directly for this
//! solver's `sod_shock_tube` validation case, which has a genuine contact
//! discontinuity between the shock and the rarefaction fan.
//!
//! Reference: Toro, "Riemann Solvers and Numerical Methods for Fluid
//! Dynamics", 3rd ed., Chapter 10.

use crate::gas::{ConservedState, Flux, GasProperties, PrimitiveState};

/// A floor on the estimated star-region pressure. Strong shocks or
/// near-vacuum states can otherwise push the pressure-based wave-speed
/// estimate below zero, which would propagate NaN through the sqrt below —
/// clamping here is cheap insurance against that failure mode.
const MIN_PRESSURE: f64 = 1e-6;

/// Computes the HLLC numerical flux across the interface between `left` and
/// `right` states of the same gas.
///
/// For identical `left == right`, this returns exactly the analytic
/// physical flux of that state (see the `flux_of_identical_states_...` unit
/// test below) — a basic consistency property every valid numerical flux
/// function must satisfy, and the reason a uniform, at-rest pipe stays
/// exactly at rest under this scheme (see `tests/quiescent.rs`).
pub fn hllc_flux(left: PrimitiveState, right: PrimitiveState, gas: &GasProperties) -> Flux {
    let left_sound_speed = left.sound_speed(gas);
    let right_sound_speed = right.sound_speed(gas);

    // Pressure-based wave speed estimate (Toro's "PVRS", eq. 10.61),
    // clamped to stay positive.
    let density_avg = 0.5 * (left.density + right.density);
    let sound_speed_avg = 0.5 * (left_sound_speed + right_sound_speed);
    let pressure_star_estimate = (0.5 * (left.pressure + right.pressure)
        - 0.5 * (right.velocity - left.velocity) * density_avg * sound_speed_avg)
        .max(MIN_PRESSURE);

    let wave_speed_left =
        left.velocity - left_sound_speed * shock_or_rarefaction_factor(pressure_star_estimate, left.pressure, gas.gamma);
    let wave_speed_right =
        right.velocity + right_sound_speed * shock_or_rarefaction_factor(pressure_star_estimate, right.pressure, gas.gamma);

    let left_conserved = left.to_conserved(gas);
    let right_conserved = right.to_conserved(gas);

    // Beyond either outer wave, the flow at the interface is just the
    // undisturbed state on that side.
    if wave_speed_left >= 0.0 {
        return left_conserved.physical_flux(gas);
    }
    if wave_speed_right <= 0.0 {
        return right_conserved.physical_flux(gas);
    }

    // Contact/entropy wave speed (Toro eq. 10.37) - the middle wave that
    // distinguishes HLLC from plain HLL.
    let wave_speed_star = {
        let numerator = right.pressure - left.pressure
            + left.density * left.velocity * (wave_speed_left - left.velocity)
            - right.density * right.velocity * (wave_speed_right - right.velocity);
        let denominator =
            left.density * (wave_speed_left - left.velocity) - right.density * (wave_speed_right - right.velocity);
        numerator / denominator
    };

    if wave_speed_star >= 0.0 {
        star_region_flux(left, left_conserved, wave_speed_left, wave_speed_star, gas)
    } else {
        star_region_flux(right, right_conserved, wave_speed_right, wave_speed_star, gas)
    }
}

/// The factor `q_K` distinguishing a shock (compression, `p_star > p_K`)
/// from a rarefaction (expansion, `p_star <= p_K`) on one side of the wave
/// (Toro eq. 10.59-10.60), used to estimate that side's outer wave speed.
fn shock_or_rarefaction_factor(pressure_star: f64, pressure_side: f64, gamma: f64) -> f64 {
    if pressure_star <= pressure_side {
        1.0
    } else {
        (1.0 + (gamma + 1.0) / (2.0 * gamma) * (pressure_star / pressure_side - 1.0)).sqrt()
    }
}

/// Flux in the "star region" between one outer wave (`side_wave_speed`) and
/// the contact wave (`star_wave_speed`), for whichever side (`side`) that
/// region borders (Toro eq. 10.38-10.39).
fn star_region_flux(
    side: PrimitiveState,
    side_conserved: ConservedState,
    side_wave_speed: f64,
    star_wave_speed: f64,
    gas: &GasProperties,
) -> Flux {
    let density_scale = side.density * (side_wave_speed - side.velocity) / (side_wave_speed - star_wave_speed);
    let specific_energy = side_conserved.energy / side.density;

    let star_state = ConservedState {
        mass: density_scale,
        momentum: density_scale * star_wave_speed,
        energy: density_scale
            * (specific_energy
                + (star_wave_speed - side.velocity)
                    * (star_wave_speed + side.pressure / (side.density * (side_wave_speed - side.velocity)))),
    };

    let side_flux = side_conserved.physical_flux(gas);
    Flux {
        mass: side_flux.mass + side_wave_speed * (star_state.mass - side_conserved.mass),
        momentum: side_flux.momentum + side_wave_speed * (star_state.momentum - side_conserved.momentum),
        energy: side_flux.energy + side_wave_speed * (star_state.energy - side_conserved.energy),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn air_at_rest(pressure: f64, temperature_kelvin: f64) -> PrimitiveState {
        PrimitiveState::from_pressure_temperature(pressure, temperature_kelvin, 0.0, &GasProperties::AIR)
    }

    #[test]
    fn flux_of_identical_states_matches_the_analytic_physical_flux() {
        let gas = GasProperties::AIR;
        let state = air_at_rest(150_000.0, 293.15);
        let numerical = hllc_flux(state, state, &gas);
        let analytic = state.to_conserved(&gas).physical_flux(&gas);
        assert!((numerical.mass - analytic.mass).abs() < 1e-9);
        assert!((numerical.momentum - analytic.momentum).abs() < 1e-6);
        assert!((numerical.energy - analytic.energy).abs() < 1e-6);
    }

    #[test]
    fn flux_of_identical_moving_states_matches_the_analytic_physical_flux() {
        let gas = GasProperties::AIR;
        let state = PrimitiveState::from_pressure_temperature(300_000.0, 350.0, 120.0, &gas);
        let numerical = hllc_flux(state, state, &gas);
        let analytic = state.to_conserved(&gas).physical_flux(&gas);
        assert!((numerical.mass - analytic.mass).abs() < 1e-6);
        assert!((numerical.momentum - analytic.momentum).abs() < 1e-3);
        assert!((numerical.energy - analytic.energy).abs() < 1.0);
    }

    #[test]
    fn strong_pressure_ratio_does_not_produce_nan_or_infinite_flux() {
        let gas = GasProperties::AIR;
        let left = air_at_rest(1000e5, 293.15); // 1000 bar
        let right = air_at_rest(1.0, 293.15); // ~vacuum
        let flux = hllc_flux(left, right, &gas);
        assert!(flux.mass.is_finite());
        assert!(flux.momentum.is_finite());
        assert!(flux.energy.is_finite());
    }

    #[test]
    fn higher_pressure_on_the_left_produces_a_rightward_mass_flux() {
        let gas = GasProperties::AIR;
        let left = air_at_rest(2e5, 293.15);
        let right = air_at_rest(1e5, 293.15);
        let flux = hllc_flux(left, right, &gas);
        assert!(flux.mass > 0.0, "expected mass to flow from high to low pressure");
    }
}
