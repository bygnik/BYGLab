//! Exact steady, isentropic, quasi-1D compressible flow through a duct of
//! varying area — used to validate the solver's taper source term under
//! genuine (non-zero-velocity) flow, complementing
//! `tests/tapered_pipe_stays_at_rest.rs`'s at-rest well-balanced check.
//!
//! For steady isentropic flow, the local Mach number M and the local
//! area-to-reference-("throat")-area ratio A/A* are related by (Anderson,
//! "Modern Compressible Flow", eq. 5.20):
//!
//!   (A/A*)^2 = (1/M^2) * [ (2/(gamma+1)) * (1 + (gamma-1)/2 * M^2) ]^((gamma+1)/(gamma-1))
//!
//! `A*` is a reference area (where M would equal 1) - it need not exist
//! physically in a purely subsonic duct, it's just a normalization. Given
//! stagnation conditions (p0, T0) and the static pressure at one station,
//! the Mach number there follows directly from the isentropic pressure
//! relation (closed form, no iteration needed); from that Mach number and
//! that station's known physical area, `A*` follows; from `A*` and any
//! other station's physical area, that station's Mach number follows by
//! inverting the area-Mach relation above (via bisection - the relation is
//! monotonic on the subsonic branch, decreasing from +infinity at M->0 to
//! 1 at M=1, so bisection is simpler and more robust here than Newton).

/// A/A* for a given Mach number and ratio of specific heats, subsonic or
/// supersonic branch alike (this module only ever calls it with M<1).
fn area_ratio(mach: f64, gamma: f64) -> f64 {
    let exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    let bracket = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * mach * mach);
    (1.0 / mach) * bracket.powf(exponent)
}

/// Solves for the subsonic Mach number (0 < M < 1) giving the requested
/// A/A*, by bisection.
fn subsonic_mach_for_area_ratio(target_area_ratio: f64, gamma: f64) -> f64 {
    assert!(target_area_ratio >= 1.0, "A/A* must be >= 1, got {target_area_ratio}");
    let mut low = 1e-6_f64;
    let mut high = 1.0 - 1e-9_f64;
    for _ in 0..200 {
        let mid = 0.5 * (low + high);
        if area_ratio(mid, gamma) > target_area_ratio {
            low = mid;
        } else {
            high = mid;
        }
    }
    0.5 * (low + high)
}

/// Mach number from the ratio of static to stagnation pressure (closed
/// form, isentropic relation).
fn mach_from_pressure_ratio(static_pressure: f64, stagnation_pressure: f64, gamma: f64) -> f64 {
    let pressure_ratio = stagnation_pressure / static_pressure;
    (2.0 / (gamma - 1.0) * (pressure_ratio.powf((gamma - 1.0) / gamma) - 1.0)).sqrt()
}

/// Static pressure at the given Mach number, from stagnation pressure `p0`.
fn static_pressure(mach: f64, stagnation_pressure: f64, gamma: f64) -> f64 {
    stagnation_pressure * (1.0 + 0.5 * (gamma - 1.0) * mach * mach).powf(-gamma / (gamma - 1.0))
}

/// Static temperature at the given Mach number, from stagnation
/// temperature `T0`.
fn static_temperature_kelvin(mach: f64, stagnation_temperature_kelvin: f64, gamma: f64) -> f64 {
    stagnation_temperature_kelvin / (1.0 + 0.5 * (gamma - 1.0) * mach * mach)
}

/// The exact steady solution for a converging (or diverging), purely
/// subsonic duct, given stagnation conditions at the (upstream) entrance
/// and a prescribed static back-pressure at the (downstream) exit.
pub struct ExactNozzleFlow {
    stagnation_pressure: f64,
    stagnation_temperature_kelvin: f64,
    gamma: f64,
    /// Reference ("throat") area - a normalization, not necessarily a
    /// physical station in this duct.
    reference_area: f64,
}

impl ExactNozzleFlow {
    /// Solves for the flow uniquely determined by stagnation conditions
    /// `(stagnation_pressure, stagnation_temperature_kelvin)`, the known
    /// physical area `exit_area` at the exit station, and the prescribed
    /// static pressure `exit_static_pressure` there.
    pub fn new(
        stagnation_pressure: f64,
        stagnation_temperature_kelvin: f64,
        exit_area: f64,
        exit_static_pressure: f64,
        gamma: f64,
    ) -> Self {
        let exit_mach = mach_from_pressure_ratio(exit_static_pressure, stagnation_pressure, gamma);
        let reference_area = exit_area / area_ratio(exit_mach, gamma);
        ExactNozzleFlow { stagnation_pressure, stagnation_temperature_kelvin, gamma, reference_area }
    }

    /// The exact (Mach number, static pressure, static temperature Kelvin)
    /// at a station with the given physical cross-sectional area.
    pub fn state_at_area(&self, area: f64) -> (f64, f64, f64) {
        let mach = subsonic_mach_for_area_ratio(area / self.reference_area, self.gamma);
        let pressure = static_pressure(mach, self.stagnation_pressure, self.gamma);
        let temperature_kelvin = static_temperature_kelvin(mach, self.stagnation_temperature_kelvin, self.gamma);
        (mach, pressure, temperature_kelvin)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn area_ratio_matches_a_known_isentropic_table_value() {
        // A well-known reference point (gamma=1.4): M=0.5 <-> A/A*=1.33984.
        let ratio = area_ratio(0.5, 1.4);
        assert!((ratio - 1.33984).abs() < 1e-4, "expected 1.33984, got {ratio}");
    }

    #[test]
    fn subsonic_mach_for_area_ratio_inverts_area_ratio() {
        let gamma = 1.4;
        for mach in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let ratio = area_ratio(mach, gamma);
            let recovered = subsonic_mach_for_area_ratio(ratio, gamma);
            assert!((recovered - mach).abs() < 1e-6, "expected {mach}, got {recovered}");
        }
    }

    #[test]
    fn mach_number_is_zero_at_the_stagnation_pressure() {
        let mach = mach_from_pressure_ratio(1.5e5, 1.5e5, 1.4);
        assert!(mach.abs() < 1e-9);
    }

    #[test]
    fn exact_nozzle_flow_is_self_consistent_at_the_exit_station() {
        // By construction, evaluating at the exit area should recover the
        // exit static pressure used to define the flow.
        let gamma = 1.4;
        let exit_area = 1.9635e-3;
        let flow = ExactNozzleFlow::new(1.5e5, 293.15, exit_area, 1.3e5, gamma);
        let (_, pressure, _) = flow.state_at_area(exit_area);
        assert!((pressure - 1.3e5).abs() < 1.0, "expected 1.3e5 Pa, got {pressure}");
    }
}
