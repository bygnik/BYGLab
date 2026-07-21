//! Wall friction and wall heat transfer: momentum/energy source terms for
//! flow through a real (rough, non-adiabatic) pipe, as opposed to the
//! inviscid, adiabatic flow `riemann.rs`/`reconstruction.rs` model.
//!
//! Applied by `pipe.rs` as an explicit, operator-split source term added to
//! the same-step conservative update — evaluated from the pre-step cell
//! state, alongside the taper geometric source term (see `pipe.rs`'s
//! `apply_face_fluxes`). Opt-in per pipe via [`WallProperties`]; a pipe with
//! no `WallProperties` is exactly frictionless and adiabatic, matching
//! every validation case that predates this module.
//!
//! Correlations are standard textbook choices (Haaland for the Darcy
//! friction factor, Colburn for the Nusselt number) rather than a port of
//! OpenWAM's own piecewise-polynomial-approximated Colebrook cascade and
//! regime-switched Nusselt correlations — simpler code, same governing
//! physics, independently re-derived and cross-checked against OpenWAM's
//! formulas (see `benchmarks/openwam/OpenWAM/Source/1DPipes/TTubo.cpp`,
//! `Colebrook`/`TransmisionCalor`).

use crate::gas::{GasProperties, PrimitiveState};

/// Below this Reynolds number, both the laminar (`64/Re`) and Haaland
/// (`6.9/Re` inside a log) friction-factor formulas are singular. The
/// physically correct momentum source at (near-)zero velocity is exactly
/// zero — not `NaN` from `0 * Inf` — so friction is skipped entirely below
/// this floor rather than evaluated. Flow reversal through zero velocity is
/// routine in pulsating engine flow, so this is a normal code path, not a
/// rare edge case.
const MIN_REYNOLDS_NUMBER_FOR_FRICTION: f64 = 1.0;

/// Wall properties of a pipe, controlling friction and heat transfer with
/// the gas inside it. `None` on [`crate::pipe::Pipe::wall`] means
/// frictionless and adiabatic — the source terms in this module are never
/// evaluated for such a pipe.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct WallProperties {
    /// Absolute wall roughness height, meters. `0.0` = hydraulically smooth
    /// (friction still applies via the smooth-pipe limit of the Haaland
    /// correlation, it just has no roughness contribution).
    pub roughness: f64,
    /// Fixed wall temperature, Kelvin. `None` means adiabatic (no heat
    /// transfer) even though friction may still be active.
    pub wall_temperature_kelvin: Option<f64>,
}

/// Reynolds number `rho * |u| * D / mu` for flow through a pipe of the
/// given diameter.
fn reynolds_number(primitive: &PrimitiveState, diameter: f64, gas: &GasProperties) -> f64 {
    primitive.density * primitive.velocity.abs() * diameter / gas.dynamic_viscosity
}

/// Darcy friction factor via the Haaland approximation to the implicit
/// Colebrook-White equation (explicit, within ~1-2% of Colebrook — the
/// standard modern substitute), with a laminar (`f = 64/Re`) branch below
/// the transitional Reynolds number 2300. Returns `None` below
/// [`MIN_REYNOLDS_NUMBER_FOR_FRICTION`], where both formulas are singular
/// but the physically correct answer is "no friction" (see that constant's
/// doc comment).
fn darcy_friction_factor(reynolds: f64, relative_roughness: f64) -> Option<f64> {
    if reynolds < MIN_REYNOLDS_NUMBER_FOR_FRICTION {
        return None;
    }
    if reynolds < 2300.0 {
        return Some(64.0 / reynolds);
    }
    let inner = (relative_roughness / 3.7).powf(1.11) + 6.9 / reynolds;
    let sqrt_inv_f = -1.8 * inner.log10();
    Some(1.0 / (sqrt_inv_f * sqrt_inv_f))
}

/// Nusselt number for internal pipe flow: the laminar, constant-wall-
/// temperature value (`3.66`, textbook exact) below the transitional
/// Reynolds number 2300, and the Colburn correlation
/// (`Nu = 0.023 * Re^0.8 * Pr^(1/3)`) above it.
fn nusselt_number(reynolds: f64, prandtl: f64) -> f64 {
    if reynolds < 2300.0 {
        3.66
    } else {
        0.023 * reynolds.powf(0.8) * prandtl.powf(1.0 / 3.0)
    }
}

/// Momentum and energy source terms (per unit volume) from wall friction
/// and wall heat transfer, evaluated at one cell's current state.
///
/// - Friction (Darcy-Weisbach): `-f * rho * u * |u| / (2 * D)`, opposing
///   the flow direction. Zero at zero velocity (both because `u*|u| = 0`
///   and because [`darcy_friction_factor`] returns `None` there).
/// - Heat transfer: `(4 * h / D) * (T_wall - T_gas)`, where `h` is the
///   convective coefficient from [`nusselt_number`] (`h = Nu * k / D`,
///   `k` from `Pr = mu * cp / k`) and `4/D` is the surface-to-volume ratio
///   of a circular duct. Zero if `wall.wall_temperature_kelvin` is `None`.
pub fn wall_sources(primitive: PrimitiveState, wall: &WallProperties, diameter: f64, gas: &GasProperties) -> (f64, f64) {
    let reynolds = reynolds_number(&primitive, diameter, gas);

    let momentum_source = match darcy_friction_factor(reynolds, wall.roughness / diameter) {
        Some(friction_factor) => {
            -friction_factor * primitive.density * primitive.velocity * primitive.velocity.abs() / (2.0 * diameter)
        }
        None => 0.0,
    };

    let energy_source = match wall.wall_temperature_kelvin {
        Some(wall_temperature_kelvin) => {
            let specific_heat_cp = gas.gamma * gas.gas_constant / (gas.gamma - 1.0);
            let thermal_conductivity = gas.dynamic_viscosity * specific_heat_cp / gas.prandtl_number;
            let nusselt = nusselt_number(reynolds, gas.prandtl_number);
            let convective_coefficient = nusselt * thermal_conductivity / diameter;
            let gas_temperature_kelvin = primitive.temperature_kelvin(gas);
            4.0 * convective_coefficient * (wall_temperature_kelvin - gas_temperature_kelvin) / diameter
        }
        None => 0.0,
    };

    (momentum_source, energy_source)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn friction_factor_matches_known_laminar_value() {
        // f = 64/Re is exact in the laminar regime, independent of roughness.
        let f = darcy_friction_factor(1000.0, 0.0).unwrap();
        assert!((f - 0.064).abs() < 1e-9, "expected 0.064, got {f}");
    }

    #[test]
    fn friction_factor_matches_moody_chart_reference_point() {
        // Smooth pipe, Re = 1e5: Moody chart / Colebrook gives f ~ 0.0180.
        let f = darcy_friction_factor(1.0e5, 0.0).unwrap();
        assert!((f - 0.0180).abs() < 0.001, "expected ~0.0180, got {f}");
    }

    #[test]
    fn friction_factor_increases_with_roughness_at_fixed_reynolds_number() {
        let smooth = darcy_friction_factor(1.0e5, 0.0).unwrap();
        let rough = darcy_friction_factor(1.0e5, 0.001).unwrap();
        assert!(rough > smooth);
    }

    #[test]
    fn friction_factor_is_none_below_the_reynolds_floor() {
        assert_eq!(darcy_friction_factor(0.1, 0.0), None);
    }

    #[test]
    fn nusselt_number_matches_laminar_constant_wall_temperature_value() {
        assert!((nusselt_number(1000.0, 0.71) - 3.66).abs() < 1e-9);
    }

    #[test]
    fn nusselt_number_matches_colburn_reference_point() {
        // Re = 1e4, Pr = 0.71: Nu = 0.023 * 10000^0.8 * 0.71^(1/3) ~= 32.5.
        let nu = nusselt_number(1.0e4, 0.71);
        assert!((nu - 32.5).abs() < 1.0, "expected ~32.5, got {nu}");
    }

    #[test]
    fn friction_source_is_zero_at_zero_velocity_not_nan() {
        let gas = GasProperties::AIR;
        let at_rest = PrimitiveState::from_pressure_temperature(150_000.0, 293.15, 0.0, &gas);
        let wall = WallProperties { roughness: 0.0001, wall_temperature_kelvin: None };
        let (momentum_source, _) = wall_sources(at_rest, &wall, 0.05, &gas);
        assert_eq!(momentum_source, 0.0, "expected exactly zero, not NaN, at u=0");
    }

    #[test]
    fn friction_source_opposes_the_flow_direction() {
        let gas = GasProperties::AIR;
        let forward = PrimitiveState::from_pressure_temperature(150_000.0, 293.15, 50.0, &gas);
        let backward = PrimitiveState::from_pressure_temperature(150_000.0, 293.15, -50.0, &gas);
        let wall = WallProperties { roughness: 0.0, wall_temperature_kelvin: None };
        let (forward_source, _) = wall_sources(forward, &wall, 0.05, &gas);
        let (backward_source, _) = wall_sources(backward, &wall, 0.05, &gas);
        assert!(forward_source < 0.0, "friction should decelerate forward flow");
        assert!(backward_source > 0.0, "friction should decelerate backward flow");
    }

    #[test]
    fn heat_transfer_source_is_zero_with_no_wall_temperature() {
        let gas = GasProperties::AIR;
        let flowing = PrimitiveState::from_pressure_temperature(150_000.0, 293.15, 20.0, &gas);
        let wall = WallProperties { roughness: 0.0, wall_temperature_kelvin: None };
        let (_, energy_source) = wall_sources(flowing, &wall, 0.05, &gas);
        assert_eq!(energy_source, 0.0);
    }

    #[test]
    fn heat_transfer_source_warms_gas_toward_a_hotter_wall() {
        let gas = GasProperties::AIR;
        let cool_gas = PrimitiveState::from_pressure_temperature(150_000.0, 293.15, 20.0, &gas);
        let wall = WallProperties { roughness: 0.0, wall_temperature_kelvin: Some(400.0) };
        let (_, energy_source) = wall_sources(cool_gas, &wall, 0.05, &gas);
        assert!(energy_source > 0.0, "a hotter wall should add energy to the gas");
    }
}
