//! Poppet valve flow: curtain area from lift, and quasi-steady
//! compressible orifice mass flow rate through it.
//!
//! Discharge coefficient is a constant scalar for this pass — this
//! matches the project's established scoping decision (see the root
//! README's roadmap): `Cd(lift/diameter)` curves are a parametric INPUT
//! in general (measured on a flow bench, or estimated), not predicted
//! from 3D port geometry, since that's what makes cylinder-head-porting
//! and valve-sizing studies tractable without running port-level CFD
//! inside the cycle solver. A real `Cd` curve is a straightforward
//! substitution later; a constant is enough to validate the underlying
//! compressible-flow physics now.

use crate::gas::GasProperties;

/// A poppet valve's seat geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ValveGeometry {
    /// Nominal seat diameter, meters.
    pub valve_diameter: f64,
    /// Seat angle, radians, measured from the seat plane (perpendicular
    /// to the valve stem) — NOT from the stem axis. A flat seat
    /// (`seat_angle_radians = 0`) should give a curtain gap approximately
    /// equal to the lift itself; a common real seat angle is 45 degrees.
    pub seat_angle_radians: f64,
}

/// Curtain area at the given lift: `pi * valve_diameter * lift *
/// cos(seat_angle_radians)`. No port-area cap for large lift (a real
/// valve's effective area eventually becomes limited by the port's own
/// inner diameter rather than growing linearly with lift forever) — a
/// known, documented simplification; stay within the physically-linear
/// regime when using this for validation.
pub fn curtain_area(geometry: &ValveGeometry, lift: f64) -> f64 {
    std::f64::consts::PI * geometry.valve_diameter * lift * geometry.seat_angle_radians.cos()
}

/// Ratio of downstream to upstream pressure below which flow through the
/// orifice is choked (sonic at the throat) — `(2/(gamma+1))^(gamma/(gamma-1))`,
/// ~0.528 for air.
fn critical_pressure_ratio(gas: &GasProperties) -> f64 {
    (2.0 / (gas.gamma + 1.0)).powf(gas.gamma / (gas.gamma - 1.0))
}

/// Quasi-steady compressible orifice mass flow rate (kg/s), given the
/// already-identified upstream state and downstream pressure. Standard
/// choked/subsonic isentropic-orifice relations (independently re-derived
/// from the isentropic area-Mach/pressure relations during design review,
/// not assumed from memory) — the same physics already validated for duct
/// flow in `tests/isentropic_nozzle.rs`, here in mass-flow-rate form.
/// Always returns a non-negative value; direction/sign is handled by the
/// caller, [`mass_flow_rate`].
fn choked_or_subsonic_mass_flow(
    upstream_pressure: f64,
    upstream_temperature_kelvin: f64,
    downstream_pressure: f64,
    area: f64,
    gas: &GasProperties,
) -> f64 {
    let gamma = gas.gamma;
    let r = gas.gas_constant;
    let p0 = upstream_pressure;
    let t0 = upstream_temperature_kelvin;

    // Always <= 1 by construction (the caller picks upstream as whichever
    // side has the higher pressure) - clamped defensively anyway, since
    // floating-point roundoff at near-equal pressures could otherwise
    // push this fractionally above 1.0 and send the subsonic branch's
    // sqrt() argument negative (-> NaN).
    let pressure_ratio = (downstream_pressure / p0).min(1.0);

    if pressure_ratio <= critical_pressure_ratio(gas) {
        area * p0 * (gamma / (r * t0)).sqrt() * (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    } else {
        area * p0 * (2.0 * gamma / (r * t0 * (gamma - 1.0))).sqrt()
            * (pressure_ratio.powf(2.0 / gamma) - pressure_ratio.powf((gamma + 1.0) / gamma)).sqrt()
    }
}

/// Signed mass flow rate (kg/s) through an orifice of `effective_area`
/// (`Cd * curtain_area`) between two lumped states — positive means net
/// flow from `side_a` to `side_b`, negative means the reverse (this
/// correctly handles back-flow, e.g. valve overlap reversion, since the
/// higher-pressure side is always treated as upstream regardless of which
/// argument position it's passed in).
///
/// Each side's *static* temperature is used directly (no stagnation
/// correction) — the standard simplification for a 0D lumped-volume-to-
/// valve interface, the same one `boundary::BoundaryCondition::Reservoir`
/// already uses at the 1D pipe-network boundary.
pub fn mass_flow_rate(
    side_a_pressure: f64,
    side_a_temperature_kelvin: f64,
    side_b_pressure: f64,
    side_b_temperature_kelvin: f64,
    effective_area: f64,
    gas: &GasProperties,
) -> f64 {
    if side_a_pressure >= side_b_pressure {
        choked_or_subsonic_mass_flow(side_a_pressure, side_a_temperature_kelvin, side_b_pressure, effective_area, gas)
    } else {
        -choked_or_subsonic_mass_flow(side_b_pressure, side_b_temperature_kelvin, side_a_pressure, effective_area, gas)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curtain_area_at_a_flat_seat_approximately_equals_lift_times_circumference() {
        // beta -> 0 (flat seat): the curtain gap should approach the
        // lift itself, so area -> pi*D*lift - this is what cos(0)=1
        // gives; sin(0)=0 would (wrongly) give zero area for a flat seat.
        let geometry = ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 0.0 };
        let lift = 0.008;
        let area = curtain_area(&geometry, lift);
        let expected = std::f64::consts::PI * geometry.valve_diameter * lift;
        assert!((area - expected).abs() < 1e-12, "expected {expected}, got {area}");
    }

    #[test]
    fn curtain_area_at_a_non_45_degree_seat_uses_cosine_not_sine() {
        // At exactly 45 degrees sin and cos are equal, hiding a sin/cos
        // mix-up entirely - a 30 degree seat distinguishes them clearly
        // (cos(30deg)=0.866, sin(30deg)=0.5, a ~42% difference).
        let geometry = ValveGeometry { valve_diameter: 0.035, seat_angle_radians: 30.0_f64.to_radians() };
        let lift = 0.008;
        let area = curtain_area(&geometry, lift);
        let expected_with_cosine = std::f64::consts::PI * geometry.valve_diameter * lift * 30.0_f64.to_radians().cos();
        let wrong_with_sine = std::f64::consts::PI * geometry.valve_diameter * lift * 30.0_f64.to_radians().sin();
        assert!((area - expected_with_cosine).abs() < 1e-12, "expected the cosine-based area {expected_with_cosine}, got {area}");
        assert!((area - wrong_with_sine).abs() > 1e-6, "area should NOT match the sine-based (wrong) formula");
    }

    #[test]
    fn choked_and_subsonic_branches_agree_exactly_at_the_critical_pressure_ratio() {
        let gas = GasProperties::AIR;
        let upstream_pressure = 10.0e5;
        let upstream_temperature = 900.0;
        let area = 1.0e-4;
        let critical_ratio = critical_pressure_ratio(&gas);
        let downstream_at_critical = upstream_pressure * critical_ratio;

        let just_below = choked_or_subsonic_mass_flow(upstream_pressure, upstream_temperature, downstream_at_critical * 0.999999, area, &gas);
        let just_above = choked_or_subsonic_mass_flow(upstream_pressure, upstream_temperature, downstream_at_critical * 1.000001, area, &gas);

        let relative_difference = (just_below - just_above).abs() / just_below;
        println!("choked/subsonic continuity: just_below={just_below:e}, just_above={just_above:e}, relative difference={relative_difference:e}");
        assert!(relative_difference < 1e-5, "expected the two branches to agree closely at the critical pressure ratio, got {relative_difference:e}");
    }

    #[test]
    fn choked_mass_flow_is_independent_of_downstream_pressure() {
        let gas = GasProperties::AIR;
        let upstream_pressure = 10.0e5;
        let upstream_temperature = 900.0;
        let area = 1.0e-4;
        let critical_ratio = critical_pressure_ratio(&gas);

        let flow_at_half_critical = choked_or_subsonic_mass_flow(upstream_pressure, upstream_temperature, upstream_pressure * critical_ratio * 0.5, area, &gas);
        let flow_at_tenth_critical = choked_or_subsonic_mass_flow(upstream_pressure, upstream_temperature, upstream_pressure * critical_ratio * 0.1, area, &gas);
        let flow_near_vacuum = choked_or_subsonic_mass_flow(upstream_pressure, upstream_temperature, 1.0, area, &gas);

        assert!((flow_at_half_critical - flow_at_tenth_critical).abs() / flow_at_half_critical < 1e-9);
        assert!((flow_at_half_critical - flow_near_vacuum).abs() / flow_at_half_critical < 1e-9);
    }

    #[test]
    fn mass_flow_rate_is_zero_at_equal_pressures() {
        let gas = GasProperties::AIR;
        let rate = mass_flow_rate(2.0e5, 400.0, 2.0e5, 300.0, 1.0e-4, &gas);
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn mass_flow_rate_is_zero_at_zero_effective_area() {
        let gas = GasProperties::AIR;
        let rate = mass_flow_rate(5.0e5, 900.0, 1.0e5, 300.0, 0.0, &gas);
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn mass_flow_rate_sign_matches_the_higher_pressure_side() {
        let gas = GasProperties::AIR;
        let area = 1.0e-4;

        // side_a (reservoir-like) higher pressure -> positive (a->b).
        let into_b = mass_flow_rate(3.0e5, 350.0, 1.0e5, 900.0, area, &gas);
        assert!(into_b > 0.0, "expected positive (a->b) flow, got {into_b}");

        // side_b higher pressure -> negative (net flow is b->a).
        let into_a = mass_flow_rate(1.0e5, 900.0, 3.0e5, 350.0, area, &gas);
        assert!(into_a < 0.0, "expected negative (b->a) flow, got {into_a}");

        // Swapping which argument holds the higher pressure should flip
        // the sign but not the magnitude - this is the exact sign-
        // convention wiring the cylinder breathing integration depends on.
        assert!((into_b + into_a).abs() < 1e-12, "magnitudes should match with opposite sign: {into_b} vs {into_a}");
    }
}
