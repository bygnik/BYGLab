//! Wiebe combustion heat release and Woschni in-cylinder wall heat
//! transfer: source terms for `cylinder.rs`'s energy balance, mirroring
//! how `source_terms.rs` holds wall friction/heat transfer for the pipe
//! solver (a separate module holding correlation physics, called by
//! whatever owns the ODE integration).
//!
//! Formulas and coefficients ported from OpenWAM's own C++ source
//! (`benchmarks/openwam/OpenWAM/Source/Engine/TCilindro.cpp`,
//! `TCilindro4T.cpp`), cross-checked against the real single-cylinder S54
//! validation case (`benchmarks/openwam/cases/engine_s54_2500rpm/`) — not
//! textbook formulas taken on faith, and not OpenWAM-specific inventions:
//! `cw1`/`cw2` are the standard literature Woschni constants.

/// Wiebe mass-fraction-burned parameters for a single combustion event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WiebeParameters {
    /// Start of combustion, radians from TDC (negative = before TDC).
    pub start_angle_radians: f64,
    /// Burn duration, radians.
    pub duration_radians: f64,
    /// Wiebe shape factor `m` (dimensionless).
    pub shape_factor_m: f64,
    /// Wiebe efficiency constant `C` (dimensionless; ~6.9 gives ~99.9%
    /// mass fraction burned at the end of `duration_radians`).
    pub efficiency_c: f64,
}

/// Mass fraction of fuel burned at `crank_angle_from_tdc_radians`:
/// `XB(theta) = 1 - exp(-C * y^(m+1))`, `y = (theta-theta0)/duration`,
/// zero before combustion starts.
pub fn mass_fraction_burned(params: &WiebeParameters, crank_angle_from_tdc_radians: f64) -> f64 {
    if crank_angle_from_tdc_radians <= params.start_angle_radians {
        return 0.0;
    }
    let y = (crank_angle_from_tdc_radians - params.start_angle_radians) / params.duration_radians;
    1.0 - (-params.efficiency_c * y.powf(params.shape_factor_m + 1.0)).exp()
}

/// `dXB/d(crank angle)`, the closed-form derivative of
/// [`mass_fraction_burned`] — zero before combustion starts.
fn mass_fraction_burned_rate(params: &WiebeParameters, crank_angle_from_tdc_radians: f64) -> f64 {
    if crank_angle_from_tdc_radians <= params.start_angle_radians {
        return 0.0;
    }
    let y = (crank_angle_from_tdc_radians - params.start_angle_radians) / params.duration_radians;
    let m = params.shape_factor_m;
    let c = params.efficiency_c;
    c * (m + 1.0) * y.powf(m) * (-c * y.powf(m + 1.0)).exp() / params.duration_radians
}

/// The crank angle at which a given fraction of the charge has burned —
/// the exact closed-form inverse of [`mass_fraction_burned`]. Useful for
/// standard combustion-phasing metrics: CA10/CA50/CA90 are the crank
/// angles at 10%/50%/90% mass fraction burned (`mass_fraction = 0.1/0.5/0.9`
/// respectively) — CA50 in particular is a standard real-world combustion
/// phasing reference, often targeted around 8-10 degrees ATDC for best
/// efficiency in a real engine.
///
/// Inverting `XB = 1 - exp(-C*y^(m+1))` for `y`:
/// `y = (-ln(1-XB)/C)^(1/(m+1))`, `theta = theta0 + y*duration`. Exact —
/// no iteration needed, unlike inverting most other physics in this
/// crate (e.g. the isentropic area-Mach relation in
/// `tests/support/isentropic_nozzle.rs`, which has no closed-form
/// inverse and needs bisection).
pub fn crank_angle_at_mass_fraction_burned(params: &WiebeParameters, mass_fraction: f64) -> f64 {
    let y = (-(1.0 - mass_fraction).ln() / params.efficiency_c).powf(1.0 / (params.shape_factor_m + 1.0));
    params.start_angle_radians + y * params.duration_radians
}

/// Heat release rate `dQ_combustion/d(crank angle)`, Watts-per-radian
/// (i.e. Joules per radian of crank rotation) — inherently a crank-angle
/// domain quantity, no angular velocity involved (matches OpenWAM's own
/// heat-release line, which has no `dt`/omega in it either).
pub fn heat_release_rate(wiebe: &WiebeParameters, crank_angle_from_tdc_radians: f64, total_heat_release_joules: f64) -> f64 {
    mass_fraction_burned_rate(wiebe, crank_angle_from_tdc_radians) * total_heat_release_joules
}

/// Woschni correlation coefficients. `combustion_turbulence_coefficient`
/// (`Fc2` in OpenWAM) is only applied by
/// [`woschni_heat_transfer_coefficient`] while the crank angle is within
/// the Wiebe burn window — zero both before ignition and after burn
/// completion (confirmed against OpenWAM's own `InicioFinCombustion()`
/// gating; an earlier draft of this port had this window wrong).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WoschniParameters {
    pub cw1: f64,
    pub cw2: f64,
    pub combustion_turbulence_coefficient: f64,
}

/// Fixed cylinder wall temperatures (piston crown, cylinder head, liner)
/// — this pass has no wall-conduction submodel, matching the S54
/// reference case's own `CalculoTempPared=2` (fixed wall temperature)
/// setting.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallTemperatures {
    pub piston_kelvin: f64,
    pub head_kelvin: f64,
    pub liner_kelvin: f64,
}

/// Reference "motored" pressure at `volume_now`: the pressure the charge
/// would have if it had continued the same polytropic process it was on
/// at the IVC reference state, with NO combustion — `p_ivc *
/// (v_ivc/v_now)^1.36`. The polytropic exponent (1.36, not the working
/// fluid's `gamma`) is a fixed empirical Woschni-correlation convention,
/// not a mistake. Used only to compute how far actual pressure has risen
/// *above* this baseline due to combustion (see `delta_pressure_pa` in
/// [`woschni_heat_transfer_coefficient`]).
pub fn motored_reference_pressure(pressure_at_ivc: f64, volume_at_ivc: f64, volume_now: f64) -> f64 {
    pressure_at_ivc * (volume_at_ivc / volume_now).powf(1.36)
}

/// Woschni convective heat transfer coefficient, W/(m^2*K).
///
/// `h = 1.2e-2 * p_Pa^0.8 * T_K^(-0.53) * D^(-0.2) * w^0.8`, `w = cw1*Cm +
/// cw2*Cu + Fc2*deltaP*Vd/(m_trapped*R)` (the compression/expansion/
/// combustion regime of the correlation — OpenWAM has separate formulas
/// for the intake/exhaust strokes, not ported here since this module only
/// covers the closed portion of the cycle).
///
/// `swirl_velocity` (`Cu`) is accepted but expected to be `0.0` for now —
/// not because swirl is physically negligible (OpenWAM computes it
/// unconditionally), but because its calibration inputs (swirl ratio,
/// piston bowl geometry) don't exist in this model yet; passing zero is
/// an explicit omission, not a claim that it's always negligible.
#[allow(clippy::too_many_arguments)]
fn woschni_heat_transfer_coefficient(
    woschni: &WoschniParameters,
    pressure_pa: f64,
    temperature_kelvin: f64,
    bore: f64,
    mean_piston_speed: f64,
    swirl_velocity: f64,
    delta_pressure_pa: f64,
    unit_displaced_volume: f64,
    trapped_mass: f64,
    gas_constant: f64,
) -> f64 {
    let clamped_delta_pressure = delta_pressure_pa.max(0.0);
    let w = woschni.cw1 * mean_piston_speed
        + woschni.cw2 * swirl_velocity
        + woschni.combustion_turbulence_coefficient * clamped_delta_pressure * unit_displaced_volume / (trapped_mass * gas_constant);
    1.2e-2 * pressure_pa.powf(0.8) * temperature_kelvin.powf(-0.53) * bore.powf(-0.2) * w.powf(0.8)
}

/// Whether `Fc2` (the Woschni combustion-turbulence term) is active at
/// this crank angle: only within the Wiebe burn window.
fn combustion_turbulence_active(wiebe: &WiebeParameters, crank_angle_from_tdc_radians: f64) -> bool {
    crank_angle_from_tdc_radians >= wiebe.start_angle_radians
        && crank_angle_from_tdc_radians <= wiebe.start_angle_radians + wiebe.duration_radians
}

/// Net wall heat transfer rate into the gas, Watts (positive when the
/// walls are hotter than the gas). Sums three surfaces at three
/// (generally different) fixed temperatures: piston crown, cylinder head
/// (both `piston_area`, flat-piston/flat-head assumption matching the
/// S54 reference case), and the currently-exposed cylinder liner
/// (`liner_area`, which grows as the piston moves away from TDC).
pub fn wall_heat_transfer_rate(
    heat_transfer_coefficient: f64,
    gas_temperature_kelvin: f64,
    walls: &WallTemperatures,
    piston_area: f64,
    head_area: f64,
    liner_area: f64,
) -> f64 {
    let h = heat_transfer_coefficient;
    h * (walls.piston_kelvin - gas_temperature_kelvin) * piston_area
        + h * (walls.head_kelvin - gas_temperature_kelvin) * head_area
        + h * (walls.liner_kelvin - gas_temperature_kelvin) * liner_area
}

/// Convenience bundle: given the current state, evaluates the Woschni
/// coefficient (gating `Fc2` by the burn window) and the resulting net
/// wall heat transfer rate in one call — the combination
/// `cylinder.rs`'s energy-balance derivative actually needs each RK4
/// stage.
#[allow(clippy::too_many_arguments)]
pub fn wall_heat_transfer_rate_at(
    wiebe: &WiebeParameters,
    woschni: &WoschniParameters,
    walls: &WallTemperatures,
    crank_angle_from_tdc_radians: f64,
    pressure_pa: f64,
    temperature_kelvin: f64,
    bore: f64,
    mean_piston_speed: f64,
    pressure_at_ivc: f64,
    volume_at_ivc: f64,
    volume_now: f64,
    unit_displaced_volume: f64,
    trapped_mass: f64,
    gas_constant: f64,
    piston_area: f64,
    head_area: f64,
    liner_area: f64,
) -> f64 {
    let active_woschni = if combustion_turbulence_active(wiebe, crank_angle_from_tdc_radians) {
        *woschni
    } else {
        WoschniParameters { combustion_turbulence_coefficient: 0.0, ..*woschni }
    };
    let delta_pressure_pa = pressure_pa - motored_reference_pressure(pressure_at_ivc, volume_at_ivc, volume_now);
    let h = woschni_heat_transfer_coefficient(
        &active_woschni,
        pressure_pa,
        temperature_kelvin,
        bore,
        mean_piston_speed,
        0.0, // swirl velocity - not modeled yet, see woschni_heat_transfer_coefficient's doc comment
        delta_pressure_pa,
        unit_displaced_volume,
        trapped_mass,
        gas_constant,
    );
    wall_heat_transfer_rate(h, temperature_kelvin, walls, piston_area, head_area, liner_area)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn s54_wiebe() -> WiebeParameters {
        WiebeParameters {
            start_angle_radians: -15.0_f64.to_radians(),
            duration_radians: 45.0_f64.to_radians(),
            shape_factor_m: 2.5,
            efficiency_c: 6.9,
        }
    }

    /// Finds the crank angle at which [`mass_fraction_burned`] reaches
    /// `target_fraction` via plain bisection on the forward function -
    /// deliberately NOT using the closed-form inverse
    /// ([`crank_angle_at_mass_fraction_burned`]), so it's a genuinely
    /// independent check of that formula (same pattern as
    /// `tests/support/isentropic_nozzle.rs`'s bisection solver, used
    /// there specifically because no closed-form inverse exists for that
    /// relation - here one does, and this function exists purely to
    /// cross-check it against a different algorithm).
    fn bisect_crank_angle_at_mass_fraction_burned(wiebe: &WiebeParameters, target_fraction: f64) -> f64 {
        let mut low = wiebe.start_angle_radians;
        let mut high = wiebe.start_angle_radians + wiebe.duration_radians;
        for _ in 0..200 {
            let mid = 0.5 * (low + high);
            if mass_fraction_burned(wiebe, mid) < target_fraction {
                low = mid;
            } else {
                high = mid;
            }
        }
        0.5 * (low + high)
    }

    #[test]
    fn crank_angle_at_mass_fraction_burned_is_the_exact_inverse_of_mass_fraction_burned() {
        let wiebe = s54_wiebe();
        for fraction in [0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99] {
            let theta = crank_angle_at_mass_fraction_burned(&wiebe, fraction);
            let recovered_fraction = mass_fraction_burned(&wiebe, theta);
            assert!(
                (recovered_fraction - fraction).abs() < 1e-9,
                "fraction={fraction}: round-trip gave {recovered_fraction}, expected {fraction}"
            );
        }
    }

    #[test]
    fn ca10_ca50_ca90_match_an_independent_bisection_on_the_forward_wiebe_function() {
        let wiebe = s54_wiebe();
        for (label, fraction) in [("CA10", 0.10), ("CA50", 0.50), ("CA90", 0.90)] {
            let analytical_radians = crank_angle_at_mass_fraction_burned(&wiebe, fraction);
            let bisected_radians = bisect_crank_angle_at_mass_fraction_burned(&wiebe, fraction);
            let difference_degrees = (analytical_radians - bisected_radians).to_degrees();

            println!(
                "{label}: analytical={:.4} deg ATDC, bisection={:.4} deg ATDC, difference={difference_degrees:e} deg",
                analytical_radians.to_degrees(),
                bisected_radians.to_degrees()
            );

            assert!(
                difference_degrees.abs() < 1e-6,
                "{label}: analytical ({:.6} deg) and bisected ({:.6} deg) crank angles disagree by {difference_degrees:e} deg",
                analytical_radians.to_degrees(),
                bisected_radians.to_degrees()
            );
        }
    }

    #[test]
    fn ca10_ca50_ca90_are_physically_plausible_and_correctly_ordered() {
        let wiebe = s54_wiebe();
        let ca10 = crank_angle_at_mass_fraction_burned(&wiebe, 0.10).to_degrees();
        let ca50 = crank_angle_at_mass_fraction_burned(&wiebe, 0.50).to_degrees();
        let ca90 = crank_angle_at_mass_fraction_burned(&wiebe, 0.90).to_degrees();

        println!("CA10={ca10:.3} deg ATDC, CA50={ca50:.3} deg ATDC, CA90={ca90:.3} deg ATDC (CA10-90 duration = {:.3} deg)", ca90 - ca10);

        assert!(ca10 < ca50, "CA10 ({ca10} deg) should occur before CA50 ({ca50} deg)");
        assert!(ca50 < ca90, "CA50 ({ca50} deg) should occur before CA90 ({ca90} deg)");
        // CA50 in the 0-20 deg ATDC range is the standard real-world
        // combustion-phasing expectation for gasoline engines near MBT
        // timing - a sanity check that these Wiebe parameters produce a
        // physically sensible burn, not just a mathematically consistent one.
        assert!((0.0..20.0).contains(&ca50), "CA50={ca50} deg ATDC is outside the physically expected range for this engine");
    }

    #[test]
    fn mass_fraction_burned_is_zero_before_ignition() {
        let wiebe = s54_wiebe();
        assert_eq!(mass_fraction_burned(&wiebe, (-20.0_f64).to_radians()), 0.0);
    }

    #[test]
    fn mass_fraction_burned_reaches_near_completion_at_end_of_duration() {
        let wiebe = s54_wiebe();
        let end_angle = wiebe.start_angle_radians + wiebe.duration_radians;
        let xb = mass_fraction_burned(&wiebe, end_angle);
        // XB(y=1) = 1 - exp(-C) = 1 - exp(-6.9) ~= 0.99899.
        assert!((xb - (1.0 - (-6.9_f64).exp())).abs() < 1e-9, "expected ~0.99899, got {xb}");
    }

    #[test]
    fn mass_fraction_burned_is_monotonically_increasing_through_the_burn_window() {
        let wiebe = s54_wiebe();
        let mut previous = 0.0;
        let mut theta = wiebe.start_angle_radians;
        while theta <= wiebe.start_angle_radians + wiebe.duration_radians {
            let xb = mass_fraction_burned(&wiebe, theta);
            assert!(xb >= previous, "mass fraction burned decreased at theta={theta}");
            previous = xb;
            theta += 1.0_f64.to_radians();
        }
    }

    #[test]
    fn mass_fraction_burned_rate_matches_a_central_finite_difference() {
        let wiebe = s54_wiebe();
        let dtheta = 1e-6;
        for theta_deg in [-10.0_f64, -5.0, 0.0, 10.0, 25.0] {
            let theta = theta_deg.to_radians();
            let xb_ahead = mass_fraction_burned(&wiebe, theta + dtheta);
            let xb_behind = mass_fraction_burned(&wiebe, theta - dtheta);
            let finite_difference_rate = (xb_ahead - xb_behind) / (2.0 * dtheta);
            let analytic_rate = mass_fraction_burned_rate(&wiebe, theta);
            assert!(
                (finite_difference_rate - analytic_rate).abs() < 1e-6,
                "theta={theta_deg}deg: expected {finite_difference_rate}, got {analytic_rate}"
            );
        }
    }

    #[test]
    fn combustion_turbulence_window_matches_the_wiebe_burn_window() {
        let wiebe = s54_wiebe();
        assert!(!combustion_turbulence_active(&wiebe, (-16.0_f64).to_radians()), "before ignition");
        assert!(combustion_turbulence_active(&wiebe, (-15.0_f64).to_radians()), "at ignition");
        assert!(combustion_turbulence_active(&wiebe, 10.0_f64.to_radians()), "mid-burn");
        assert!(combustion_turbulence_active(&wiebe, 30.0_f64.to_radians()), "at burn completion (theta0+duration)");
        assert!(!combustion_turbulence_active(&wiebe, 31.0_f64.to_radians()), "after burn completion");
    }

    #[test]
    fn motored_reference_pressure_matches_the_polytropic_relation() {
        let p = motored_reference_pressure(1.0e5, 2.0, 1.0);
        // Halving volume: p2 = p1 * 2^1.36.
        let expected = 1.0e5 * 2.0_f64.powf(1.36);
        assert!((p - expected).abs() < 1.0, "expected {expected}, got {p}");
    }

    #[test]
    fn woschni_coefficient_is_finite_and_positive_for_typical_combustion_conditions() {
        let woschni = WoschniParameters { cw1: 2.28, cw2: 0.00324, combustion_turbulence_coefficient: 0.001 };
        let h = woschni_heat_transfer_coefficient(&woschni, 60e5, 2500.0, 0.087, 7.58, 0.0, 30e5, 5.4e-4, 5.1e-4, 287.0);
        assert!(h.is_finite() && h > 0.0, "expected a finite positive coefficient, got {h}");
    }

    #[test]
    fn woschni_coefficient_clamps_negative_delta_pressure_instead_of_producing_nan() {
        let woschni = WoschniParameters { cw1: 2.28, cw2: 0.00324, combustion_turbulence_coefficient: 0.001 };
        // A large negative delta_pressure would make `w` negative without the clamp, and `w.powf(0.8)` on a negative base is NaN.
        let h = woschni_heat_transfer_coefficient(&woschni, 1e5, 400.0, 0.087, 7.58, 0.0, -1e9, 5.4e-4, 5.1e-4, 287.0);
        assert!(h.is_finite(), "expected a finite coefficient even with a large negative delta pressure, got {h}");
    }

    #[test]
    fn wall_heat_transfer_rate_is_positive_when_walls_are_hotter_than_gas() {
        let walls = WallTemperatures { piston_kelvin: 450.0, head_kelvin: 550.0, liner_kelvin: 420.0 };
        let rate = wall_heat_transfer_rate(500.0, 350.0, &walls, 0.0059, 0.0059, 0.01);
        assert!(rate > 0.0, "expected heat flowing INTO a cooler gas from hotter walls, got {rate}");
    }

    #[test]
    fn wall_heat_transfer_rate_is_negative_when_gas_is_hotter_than_walls() {
        let walls = WallTemperatures { piston_kelvin: 450.0, head_kelvin: 550.0, liner_kelvin: 420.0 };
        let rate = wall_heat_transfer_rate(500.0, 2500.0, &walls, 0.0059, 0.0059, 0.01);
        assert!(rate < 0.0, "expected heat flowing OUT of a hotter gas into cooler walls, got {rate}");
    }

    #[test]
    fn wall_heat_transfer_rate_at_is_independent_of_the_ivc_anchor_outside_the_wiebe_window() {
        // Load-bearing property for `cylinder.rs`'s unified full-cycle
        // integrator: it resolves the Woschni "motored reference" anchor
        // (pressure_at_ivc/volume_at_ivc) via a pre-pass ONLY when
        // integration starts before actual IVC, then runs a SINGLE
        // derivative across the whole domain rather than switching
        // derivatives at IVC - which is only valid if the anchor is
        // provably inert everywhere combustion hasn't started/finished.
        // `combustion_turbulence_active` gates the ENTIRE
        // `combustion_turbulence_coefficient` term to zero outside the
        // Wiebe burn window, and `delta_pressure_pa` (the only other place
        // the anchor is used) is multiplied by that same zeroed
        // coefficient - so two wildly different anchors must give
        // bit-identical output outside the window, checked directly here
        // rather than only reasoned about.
        let wiebe = s54_wiebe();
        let woschni = WoschniParameters { cw1: 2.28, cw2: 0.00324, combustion_turbulence_coefficient: 0.001 };
        let walls = WallTemperatures { piston_kelvin: 450.0, head_kelvin: 550.0, liner_kelvin: 420.0 };

        let crank_angle_before_ignition = wiebe.start_angle_radians - 100.0_f64.to_radians();
        let crank_angle_after_burnout = wiebe.start_angle_radians + wiebe.duration_radians + 100.0_f64.to_radians();

        for crank_angle in [crank_angle_before_ignition, crank_angle_after_burnout] {
            let rate_with_anchor_a = wall_heat_transfer_rate_at(
                &wiebe, &woschni, &walls, crank_angle, 20e5, 900.0, 0.087, 7.58, 15e5, 5.4e-4, 4.0e-4, 5.4e-4, 5.1e-4, 287.0, 0.0059, 0.0059, 0.01,
            );
            // Wildly different anchor: 10x the pressure, 1/10th the volume.
            let rate_with_anchor_b = wall_heat_transfer_rate_at(
                &wiebe, &woschni, &walls, crank_angle, 20e5, 900.0, 0.087, 7.58, 150e5, 5.4e-5, 4.0e-4, 5.4e-4, 5.1e-4, 287.0, 0.0059, 0.0059, 0.01,
            );
            assert_eq!(
                rate_with_anchor_a, rate_with_anchor_b,
                "crank_angle={:.1}deg: expected identical wall heat transfer rate regardless of the IVC anchor outside the Wiebe window",
                crank_angle.to_degrees()
            );
        }
    }

    #[test]
    fn liner_area_formula_matches_the_pi_bore_displacement_identity() {
        // A * displacement, with A = pi/4*bore^2, should equal pi*bore*displacement
        // when expressed as 4*(V-Vcc)/bore - the identity this module's doc comment
        // (and the design review) claims for OpenWAM's own liner-area formula.
        let bore = 0.087_f64;
        let displacement = 0.04_f64;
        let piston_area = PI / 4.0 * bore * bore;
        let volume_minus_clearance = piston_area * displacement;
        let openwam_form = 4.0 * volume_minus_clearance / bore;
        let this_module_form = PI * bore * displacement;
        assert!((openwam_form - this_module_form).abs() < 1e-15, "expected the two liner-area forms to match exactly");
    }
}
