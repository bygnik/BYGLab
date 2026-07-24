//! Camshaft lift profile: valve lift as a function of crank angle.
//!
//! A single valve event is modeled as a raised-cosine ("versine") lift
//! profile — a standard first-approximation cam shape, smooth (zero lift
//! AND zero lift-rate at both the opening and closing events, so a valve
//! seats/unseats without an instantaneous velocity jump). Real production
//! cams use more sophisticated jerk-limited profiles with dedicated ramp
//! sections for quiet seating — a direct measured-lift-table variant is
//! the natural future addition for matching a specific real cam, not
//! needed for this pass.

use std::f64::consts::PI;

/// A single valve event: maximum lift, the crank angle (from TDC) at
/// which the valve begins to open, and the total duration it stays open.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CamProfile {
    /// Maximum valve lift, meters.
    pub max_lift: f64,
    /// Crank angle from TDC at which the valve begins to lift, radians.
    pub opening_angle_radians: f64,
    /// Total open duration, radians.
    pub duration_radians: f64,
}

impl CamProfile {
    /// Returns a copy shifted by a constant crank-angle offset - the whole
    /// event moves together, `duration_radians`/`max_lift` unchanged.
    ///
    /// Needed because real ingested cam data (`camshaft_presets.rs`) is
    /// referenced to *gas-exchange* TDC = 0, while [`crate::combustion`]'s
    /// Wiebe timing is referenced to *firing* TDC = 0 - combining the two
    /// in one simulation requires re-expressing one convention in terms
    /// of the other first. The needed shift is `+-2*pi`, NOT the same
    /// sign for intake and exhaust: exhaust events happen *after*
    /// combustion (their own gas-exchange TDC occurrence is ahead, so
    /// `+2*pi`), intake events happen *before* combustion (their relevant
    /// occurrence is behind, so `-2*pi`) - confirmed against real Schrick
    /// numbers and independently cross-checked against the OpenWAM S54
    /// case file (which already encodes valve timing directly in
    /// firing-TDC terms) in `camshaft_presets.rs`'s own tests.
    pub fn shifted_by(&self, offset_radians: f64) -> CamProfile {
        CamProfile { opening_angle_radians: self.opening_angle_radians + offset_radians, ..*self }
    }
}

/// Valve lift at `crank_angle_from_tdc_radians`: zero outside
/// `[opening, opening+duration]`, otherwise
/// `0.5 * max_lift * (1 - cos(2*pi*(theta-opening)/duration))` — a
/// half-cosine hump reaching exactly `max_lift` at the midpoint and
/// returning smoothly (zero slope) to zero at both ends.
pub fn lift_at(profile: &CamProfile, crank_angle_from_tdc_radians: f64) -> f64 {
    let closing_angle_radians = profile.opening_angle_radians + profile.duration_radians;
    if crank_angle_from_tdc_radians < profile.opening_angle_radians || crank_angle_from_tdc_radians > closing_angle_radians {
        return 0.0;
    }
    let fraction = (crank_angle_from_tdc_radians - profile.opening_angle_radians) / profile.duration_radians;
    0.5 * profile.max_lift * (1.0 - (2.0 * PI * fraction).cos())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_profile() -> CamProfile {
        CamProfile { max_lift: 0.010, opening_angle_radians: (-10.0_f64).to_radians(), duration_radians: 220.0_f64.to_radians() }
    }

    #[test]
    fn lift_is_zero_before_opening_and_after_closing() {
        let profile = sample_profile();
        assert_eq!(lift_at(&profile, profile.opening_angle_radians - 0.01), 0.0);
        assert_eq!(lift_at(&profile, profile.opening_angle_radians + profile.duration_radians + 0.01), 0.0);
    }

    #[test]
    fn lift_is_exactly_zero_at_the_opening_and_closing_angles() {
        let profile = sample_profile();
        assert_eq!(lift_at(&profile, profile.opening_angle_radians), 0.0);
        assert_eq!(lift_at(&profile, profile.opening_angle_radians + profile.duration_radians), 0.0);
    }

    #[test]
    fn lift_reaches_exactly_max_lift_at_the_midpoint() {
        let profile = sample_profile();
        let midpoint = profile.opening_angle_radians + 0.5 * profile.duration_radians;
        assert!((lift_at(&profile, midpoint) - profile.max_lift).abs() < 1e-12);
    }

    #[test]
    fn lift_is_symmetric_about_the_midpoint() {
        let profile = sample_profile();
        let midpoint = profile.opening_angle_radians + 0.5 * profile.duration_radians;
        for offset_deg in [5.0_f64, 30.0, 60.0, 100.0] {
            let offset = offset_deg.to_radians();
            let before = lift_at(&profile, midpoint - offset);
            let after = lift_at(&profile, midpoint + offset);
            assert!((before - after).abs() < 1e-12, "offset={offset_deg}deg: expected symmetric lift, got {before} vs {after}");
        }
    }

    #[test]
    fn lift_rate_is_zero_at_the_opening_and_closing_angles() {
        // Finite-difference check that d(lift)/dtheta -> 0 at both event
        // boundaries (the "smooth seating" property) - a real, checkable
        // consequence of the cosine profile's construction, not assumed.
        // Each boundary's OUTSIDE direction is exactly zero by the early
        // return (opening: backward; closing: forward) - the other
        // (inside) direction is checked via a near-zero finite-difference
        // rate instead, since it's not exactly zero there, just close.
        let profile = sample_profile();
        let dtheta = 1e-6;
        let opening = profile.opening_angle_radians;
        let closing = profile.opening_angle_radians + profile.duration_radians;

        assert_eq!(lift_at(&profile, opening - dtheta), 0.0, "just before opening should be exactly zero");
        assert_eq!(lift_at(&profile, closing + dtheta), 0.0, "just after closing should be exactly zero");

        let rate_at_opening = (lift_at(&profile, opening + dtheta) - lift_at(&profile, opening)) / dtheta;
        let rate_at_closing = (lift_at(&profile, closing) - lift_at(&profile, closing - dtheta)) / dtheta;
        assert!(rate_at_opening.abs() < 1e-3, "expected near-zero lift rate at opening, got {rate_at_opening}");
        assert!(rate_at_closing.abs() < 1e-3, "expected near-zero lift rate at closing, got {rate_at_closing}");
    }

    #[test]
    fn lift_is_never_negative_across_the_whole_event() {
        let profile = sample_profile();
        let mut theta = profile.opening_angle_radians;
        while theta <= profile.opening_angle_radians + profile.duration_radians {
            assert!(lift_at(&profile, theta) >= 0.0, "negative lift at theta={theta}");
            theta += 1.0_f64.to_radians();
        }
    }

    #[test]
    fn shifted_by_moves_the_whole_event_and_leaves_shape_unchanged() {
        let profile = sample_profile();
        let offset = 2.0 * PI;
        let shifted = profile.shifted_by(offset);

        assert_eq!(shifted.max_lift, profile.max_lift);
        assert_eq!(shifted.duration_radians, profile.duration_radians);
        assert!((shifted.opening_angle_radians - (profile.opening_angle_radians + offset)).abs() < 1e-12);

        // The lift curve is identical, just relabeled at a shifted angle.
        for offset_deg in [0.0_f64, 30.0, 110.0, 219.9] {
            let theta = profile.opening_angle_radians + offset_deg.to_radians();
            let original_lift = lift_at(&profile, theta);
            let shifted_lift = lift_at(&shifted, theta + offset);
            assert!(
                (original_lift - shifted_lift).abs() < 1e-12,
                "theta_offset={offset_deg}deg: expected {original_lift}, got {shifted_lift}"
            );
        }
    }
}
