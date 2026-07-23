//! Slider-crank mechanism kinematics: piston position, velocity, and
//! acceleration as exact functions of crank angle — the geometric
//! foundation the 0D cylinder model (phase 2) builds on.
//!
//! Uses the *finite* connecting-rod-length formula throughout (not the
//! infinite-rod/simple-harmonic-motion approximation, and not a truncated
//! series expansion). Supports piston pin offset ("desaxé" mechanisms,
//! used in production engines to reduce piston slap): with a nonzero
//! offset, true TDC no longer occurs exactly when the crank pin is
//! collinear with the cylinder axis — a genuine second-order effect this
//! module accounts for exactly (via a one-time Newton solve at
//! construction), not approximated away.
//!
//! All angles are radians (SI convention, matching the rest of this
//! crate); a config/UI layer converts to/from degrees at its boundary.

/// A slider-crank mechanism: crank radius, connecting rod length, and
/// (optional) piston pin offset.
///
/// Deliberately not `Serialize`/`Deserialize` — `reference_angle_at_top_dead_center`
/// is a cached result of a Newton solve done in [`Self::new`], not user
/// input; constructing this struct any other way (e.g. deserializing a
/// stale/mismatched value for it) would silently break every method below.
/// A future JSON-facing "spec" struct (mirroring how `case::PipeSpec`
/// relates to `Pipe`) should hold the three physical inputs and build this
/// via `new()`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrankMechanism {
    /// Crank radius (half the stroke), meters.
    pub crank_radius: f64,
    /// Connecting rod length, center-to-center, meters.
    pub rod_length: f64,
    /// Piston pin offset from the cylinder axis ("desaxé" distance),
    /// meters — the perpendicular distance between the cylinder axis and
    /// the crank's center of rotation. `0.0` for a conventional
    /// (non-offset) mechanism. Sign follows the same convention as the
    /// crank pin's perpendicular displacement in
    /// [`Self::position_at_reference_angle`]'s derivation.
    pub piston_pin_offset: f64,
    /// The crank's rotation angle, measured from the "crank pin collinear
    /// with the cylinder axis" reference (called `psi` in this module's
    /// internal methods), at which the piston actually reaches TDC. Zero
    /// when `piston_pin_offset` is zero (by symmetry); otherwise a small
    /// angle found by Newton's method at construction time (see
    /// [`Self::solve_for_true_top_dead_center`]), so every public method's
    /// `crank_angle_from_tdc` parameter is measured from *true* TDC — the
    /// same convention a physical TDC sensor uses — rather than an
    /// offset mechanism's slightly-shifted geometric reference.
    reference_angle_at_top_dead_center: f64,
}

impl CrankMechanism {
    /// Builds a mechanism from stroke (not crank radius directly — stroke
    /// is the number usually quoted for an engine), connecting rod
    /// length, and piston pin offset (`0.0` for a conventional mechanism).
    pub fn new(stroke: f64, rod_length: f64, piston_pin_offset: f64) -> Self {
        let mut mechanism =
            CrankMechanism { crank_radius: stroke / 2.0, rod_length, piston_pin_offset, reference_angle_at_top_dead_center: 0.0 };
        mechanism.reference_angle_at_top_dead_center = mechanism.solve_for_true_extremum(0.0);
        mechanism
    }

    /// Converts a crank angle measured from true TDC to this mechanism's
    /// internal reference frame (see
    /// [`Self::reference_angle_at_top_dead_center`]'s doc comment).
    fn to_reference_angle(&self, crank_angle_from_tdc_radians: f64) -> f64 {
        crank_angle_from_tdc_radians + self.reference_angle_at_top_dead_center
    }

    /// Piston pin position along the cylinder axis, measured from the
    /// crank's center of rotation, at reference angle `psi` (radians).
    ///
    /// Derived from the rigid connecting-rod constraint: the crank pin
    /// traces a circle of radius `r` about the crank center,
    /// `(r sin(psi), r cos(psi))`; the piston pin is constrained to the
    /// (offset) cylinder axis, `(offset, x)`; the rod length between them
    /// is fixed at `L`. Substituting into
    /// `(offset - r sin(psi))^2 + (x - r cos(psi))^2 = L^2` and solving
    /// for `x` (taking the positive root, since the piston pin sits
    /// "above" the crank pin toward TDC):
    ///
    /// `x(psi) = r cos(psi) + sqrt(L^2 - (offset - r sin(psi))^2)`
    ///
    /// Exact for any rod length and offset — not a small-angle,
    /// infinite-rod, or small-offset approximation. Reduces to the
    /// standard textbook slider-crank formula when `offset = 0`.
    fn position_at_reference_angle(&self, psi: f64) -> f64 {
        let r = self.crank_radius;
        let g = self.piston_pin_offset - r * psi.sin();
        r * psi.cos() + (self.rod_length * self.rod_length - g * g).sqrt()
    }

    /// `dx/d(psi)`, in closed form (analytic derivative of
    /// [`Self::position_at_reference_angle`] — not a finite difference).
    fn position_derivative_at_reference_angle(&self, psi: f64) -> f64 {
        let r = self.crank_radius;
        let g = self.piston_pin_offset - r * psi.sin();
        let s = (self.rod_length * self.rod_length - g * g).sqrt();
        -r * psi.sin() + r * psi.cos() * g / s
    }

    /// `d^2x/d(psi)^2`, in closed form (analytic second derivative — not a
    /// finite difference). Built from the same intermediate quantities
    /// (`g`, `s`) as [`Self::position_derivative_at_reference_angle`],
    /// differentiated once more by the product/quotient rule.
    fn position_second_derivative_at_reference_angle(&self, psi: f64) -> f64 {
        let r = self.crank_radius;
        let g = self.piston_pin_offset - r * psi.sin();
        let g_prime = -r * psi.cos();
        let s = (self.rod_length * self.rod_length - g * g).sqrt();
        let s_prime = r * g * psi.cos() / s;

        // The second term of dx/dpsi is n/s, with n = r*g*cos(psi).
        let n = r * g * psi.cos();
        let n_prime = r * (g_prime * psi.cos() - g * psi.sin());

        -r * psi.cos() + (n_prime * s - n * s_prime) / (s * s)
    }

    /// Finds the reference-frame angle `psi` nearest `initial_guess` at
    /// which the piston reaches a true extremum of position (`dx/dpsi =
    /// 0`) — TDC near `psi=0`, BDC near `psi=pi`. Exactly `initial_guess`
    /// when `piston_pin_offset` is zero (by symmetry); Newton's method
    /// converges in a handful of iterations otherwise, since the
    /// offset-induced shift is small for any physically realistic offset.
    ///
    /// Degenerate case: a zero crank radius (`crank_radius == 0.0`) is a
    /// legitimate way to model a rigid, non-moving 0D element (nothing
    /// ever displaces, e.g. a fixed-volume control volume built by reusing
    /// this same kinematics code rather than a separate code path) — but
    /// there both `dx/dpsi` and `d^2x/dpsi^2` are identically zero
    /// (position is constant everywhere), so the Newton step would
    /// otherwise divide `0.0/0.0 = NaN` and poison every subsequent
    /// calculation. Every angle is trivially an "extremum" of a constant
    /// function, so `initial_guess` itself is as valid a reference as any
    /// other — returned directly, skipping the Newton iteration entirely.
    fn solve_for_true_extremum(&self, initial_guess: f64) -> f64 {
        if self.crank_radius == 0.0 {
            return initial_guess;
        }
        let mut psi = initial_guess;
        for _ in 0..50 {
            let velocity_term = self.position_derivative_at_reference_angle(psi);
            let acceleration_term = self.position_second_derivative_at_reference_angle(psi);
            let step = velocity_term / acceleration_term;
            psi -= step;
            if step.abs() < 1e-14 {
                break;
            }
        }
        psi
    }

    /// Piston pin position along the cylinder axis, measured from the
    /// crank's center of rotation, at `crank_angle_from_tdc_radians`
    /// radians past true TDC.
    pub fn piston_position_from_crank_center(&self, crank_angle_from_tdc_radians: f64) -> f64 {
        self.position_at_reference_angle(self.to_reference_angle(crank_angle_from_tdc_radians))
    }

    /// Piston displacement from true TDC (`>= 0`, increasing toward BDC)
    /// — the quantity that maps directly to swept cylinder volume:
    /// `V(theta) = V_clearance + piston_area * displacement_from_tdc(theta)`.
    pub fn piston_displacement_from_top_dead_center(&self, crank_angle_from_tdc_radians: f64) -> f64 {
        self.piston_position_from_crank_center(0.0) - self.piston_position_from_crank_center(crank_angle_from_tdc_radians)
    }

    /// Instantaneous piston velocity (m/s) at `crank_angle_from_tdc_radians`
    /// past TDC, for a crank rotating at angular velocity
    /// `angular_velocity_radians_per_second` (rad/s) — `d(position)/dt`,
    /// using the same sign convention as
    /// [`Self::piston_position_from_crank_center`] (i.e. *not* flipped to
    /// read positive-away-from-TDC): negative while the piston moves away
    /// from TDC toward BDC (for positive angular velocity), positive on
    /// the return stroke.
    pub fn piston_velocity(&self, crank_angle_from_tdc_radians: f64, angular_velocity_radians_per_second: f64) -> f64 {
        let psi = self.to_reference_angle(crank_angle_from_tdc_radians);
        self.position_derivative_at_reference_angle(psi) * angular_velocity_radians_per_second
    }

    /// Instantaneous piston acceleration (m/s^2). Assumes *constant*
    /// angular velocity (no angular acceleration) — the standard,
    /// appropriate assumption for a steady-state engine cycle evaluated
    /// at one fixed operating point (RPM), which is how this crate's
    /// simulations are structured (see the root README's roadmap, phase
    /// 5's "operating-point sweeps"). If angular velocity itself varies
    /// with time, a `(dx/dpsi) * angular_acceleration` term would need
    /// adding — not implemented, since it's not needed by that use case.
    pub fn piston_acceleration(&self, crank_angle_from_tdc_radians: f64, angular_velocity_radians_per_second: f64) -> f64 {
        let psi = self.to_reference_angle(crank_angle_from_tdc_radians);
        self.position_second_derivative_at_reference_angle(psi) * angular_velocity_radians_per_second * angular_velocity_radians_per_second
    }

    /// The true BDC crank angle, measured from true TDC (radians) — `pi`
    /// when `piston_pin_offset` is zero; otherwise slightly different,
    /// since TDC and BDC each get an independent (generally unequal)
    /// phase correction under offset.
    pub fn crank_angle_of_bottom_dead_center(&self) -> f64 {
        // `solve_for_true_extremum` returns an absolute reference-frame
        // angle (psi), not an offset from `initial_guess` - converting to
        // the "crank angle from TDC" convention (see `to_reference_angle`)
        // means subtracting `reference_angle_at_top_dead_center` once, not
        // adding another `pi` on top of the already-absolute result.
        let reference_angle_at_bdc = self.solve_for_true_extremum(std::f64::consts::PI);
        reference_angle_at_bdc - self.reference_angle_at_top_dead_center
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// S54B32 reference geometry (see root README's spec table): bore
    /// isn't needed for pure kinematics, but stroke/rod length are real,
    /// validated numbers, not invented ones. No documented pin offset for
    /// this engine, so 0.0 here; offset-specific behavior is covered by
    /// the dedicated offset tests below with a deliberately exercised
    /// nonzero value.
    fn s54b32_mechanism() -> CrankMechanism {
        CrankMechanism::new(0.091, 0.139, 0.0)
    }

    #[test]
    fn no_offset_position_matches_the_standard_slider_crank_formula() {
        // Heywood, "Internal Combustion Engine Fundamentals", eq. 2.3:
        // s/a = (1 - cos(theta)) + R - sqrt(R^2 - sin^2(theta)), R = L/a.
        let mechanism = s54b32_mechanism();
        let r = mechanism.crank_radius;
        let big_r = mechanism.rod_length / r;
        for theta_deg in [0.0_f64, 30.0, 90.0, 150.0, 180.0, 270.0, 350.0] {
            let theta = theta_deg.to_radians();
            let expected_displacement = r * ((1.0 - theta.cos()) + big_r - (big_r * big_r - theta.sin().powi(2)).sqrt());
            let actual_displacement = mechanism.piston_displacement_from_top_dead_center(theta);
            assert!(
                (actual_displacement - expected_displacement).abs() < 1e-12,
                "theta={theta_deg}: expected {expected_displacement}, got {actual_displacement}"
            );
        }
    }

    #[test]
    fn no_offset_velocity_matches_the_standard_closed_form_formula() {
        // dx/dtheta = -r*[sin(theta) + sin(theta)*cos(theta)/sqrt(R^2-sin^2(theta))],
        // independently re-derived from this module's general (offset-aware)
        // formula and cross-checked by hand against Heywood eq. 2.4 before
        // writing this test.
        let mechanism = s54b32_mechanism();
        let r = mechanism.crank_radius;
        let big_r = mechanism.rod_length / r;
        let omega = 500.0; // rad/s, arbitrary nonzero value
        for theta_deg in [10.0_f64, 45.0, 90.0, 135.0, 225.0, 300.0] {
            let theta = theta_deg.to_radians();
            let expected_dx_dtheta =
                -r * (theta.sin() + theta.sin() * theta.cos() / (big_r * big_r - theta.sin().powi(2)).sqrt());
            let expected_velocity = expected_dx_dtheta * omega;
            let actual_velocity = mechanism.piston_velocity(theta, omega);
            assert!(
                (actual_velocity - expected_velocity).abs() < 1e-9,
                "theta={theta_deg}: expected {expected_velocity}, got {actual_velocity}"
            );
        }
    }

    #[test]
    fn velocity_is_exactly_zero_at_true_top_dead_center_with_and_without_offset() {
        for offset in [0.0, 0.001, -0.0015] {
            let mechanism = CrankMechanism::new(0.091, 0.139, offset);
            let velocity = mechanism.piston_velocity(0.0, 500.0);
            assert!(velocity.abs() < 1e-9, "offset={offset}: expected ~0, got {velocity}");
        }
    }

    #[test]
    fn velocity_is_exactly_zero_at_true_bottom_dead_center_with_and_without_offset() {
        for offset in [0.0, 0.001, -0.0015] {
            let mechanism = CrankMechanism::new(0.091, 0.139, offset);
            let bdc_angle = mechanism.crank_angle_of_bottom_dead_center();
            let velocity = mechanism.piston_velocity(bdc_angle, 500.0);
            assert!(velocity.abs() < 1e-9, "offset={offset}: expected ~0 at bdc_angle={bdc_angle}, got {velocity}");
        }
    }

    #[test]
    fn velocity_matches_a_central_finite_difference_of_position() {
        let mechanism = CrankMechanism::new(0.091, 0.139, 0.0015);
        let omega = 500.0;
        let dtheta = 1e-6;
        for theta_deg in [15.0_f64, 75.0, 100.0, 200.0, 340.0] {
            let theta = theta_deg.to_radians();
            let position_ahead = mechanism.piston_position_from_crank_center(theta + dtheta);
            let position_behind = mechanism.piston_position_from_crank_center(theta - dtheta);
            let finite_difference_dx_dtheta = (position_ahead - position_behind) / (2.0 * dtheta);
            let expected_velocity = finite_difference_dx_dtheta * omega;
            let actual_velocity = mechanism.piston_velocity(theta, omega);
            let relative_error = (actual_velocity - expected_velocity).abs() / actual_velocity.abs();
            assert!(relative_error < 1e-6, "theta={theta_deg}: expected {expected_velocity}, got {actual_velocity}");
        }
    }

    #[test]
    fn acceleration_matches_a_central_finite_difference_of_velocity() {
        let mechanism = CrankMechanism::new(0.091, 0.139, 0.0015);
        let omega = 500.0;
        let dtheta = 1e-6;
        for theta_deg in [15.0_f64, 75.0, 100.0, 200.0, 340.0] {
            let theta = theta_deg.to_radians();
            let velocity_ahead = mechanism.piston_velocity(theta + dtheta, omega);
            let velocity_behind = mechanism.piston_velocity(theta - dtheta, omega);
            let finite_difference_dv_dtheta = (velocity_ahead - velocity_behind) / (2.0 * dtheta);
            let expected_acceleration = finite_difference_dv_dtheta * omega;
            let actual_acceleration = mechanism.piston_acceleration(theta, omega);
            let relative_error = (actual_acceleration - expected_acceleration).abs() / actual_acceleration.abs();
            assert!(relative_error < 1e-5, "theta={theta_deg}: expected {expected_acceleration}, got {actual_acceleration}");
        }
    }

    #[test]
    fn infinite_rod_limit_approaches_simple_harmonic_motion() {
        // As L/r -> infinity, the connecting rod's obliquity vanishes and
        // the finite-rod formula must approach pure sinusoidal motion,
        // x(theta) -> r*cos(theta) + const - the textbook infinite-rod
        // approximation many simplified engine models use instead of the
        // exact formula this module implements.
        let very_long_rod_mechanism = CrankMechanism::new(0.091, 1000.0, 0.0);
        let r = very_long_rod_mechanism.crank_radius;
        for theta_deg in [30.0_f64, 90.0, 150.0] {
            let theta = theta_deg.to_radians();
            let expected_shm_displacement = r * (1.0 - theta.cos());
            let actual_displacement = very_long_rod_mechanism.piston_displacement_from_top_dead_center(theta);
            let relative_error = (actual_displacement - expected_shm_displacement).abs() / (2.0 * r);
            assert!(relative_error < 1e-4, "theta={theta_deg}: SHM limit not approached closely enough");
        }
    }

    #[test]
    fn nonzero_offset_converges_to_the_no_offset_case_as_offset_shrinks() {
        let baseline = s54b32_mechanism(); // offset = 0.0
        let tiny_offset = CrankMechanism::new(0.091, 0.139, 1e-8);
        for theta_deg in [45.0_f64, 135.0, 225.0] {
            let theta = theta_deg.to_radians();
            let baseline_position = baseline.piston_position_from_crank_center(theta);
            let tiny_offset_position = tiny_offset.piston_position_from_crank_center(theta);
            assert!(
                (baseline_position - tiny_offset_position).abs() < 1e-6,
                "theta={theta_deg}: a 1e-8 m offset should barely perturb position"
            );
        }
    }

    #[test]
    fn realistic_offset_produces_a_small_top_dead_center_phase_correction() {
        // A physically realistic piston pin offset (1 mm, on a 45.5 mm
        // crank radius / 139 mm rod) should shift true TDC by a small
        // fraction of a degree, not several degrees - a sanity check on
        // the Newton solve's result, not just "some number came out."
        let mechanism = CrankMechanism::new(0.091, 0.139, 0.001);
        let phase_correction_degrees = mechanism.reference_angle_at_top_dead_center.to_degrees();
        println!("1mm pin offset -> TDC phase correction {phase_correction_degrees:.4} deg");
        assert!(phase_correction_degrees.abs() < 1.0, "expected a sub-degree correction, got {phase_correction_degrees} deg");
        assert!(phase_correction_degrees.abs() > 1e-4, "expected a nonzero correction for a nonzero offset");
    }

    #[test]
    fn total_piston_travel_from_true_tdc_to_true_bdc_is_close_to_the_nominal_stroke() {
        // With no offset this is an exact identity (stroke = 2 * crank
        // radius, by definition). With offset, TDC and BDC each get an
        // independent phase correction, so this is measured rather than
        // assumed to be exact - printed for visibility, asserted against
        // a generous but real bound.
        let stroke = 0.091;
        for offset in [0.0, 0.001, -0.0015] {
            let mechanism = CrankMechanism::new(stroke, 0.139, offset);
            let bdc_angle = mechanism.crank_angle_of_bottom_dead_center();
            let travel = mechanism.piston_displacement_from_top_dead_center(bdc_angle);
            let relative_deviation = (travel - stroke).abs() / stroke;
            println!("offset={offset}: TDC-to-BDC travel = {travel:.9} m (nominal stroke {stroke} m), relative deviation {relative_deviation:e}");
            if offset == 0.0 {
                assert!(relative_deviation < 1e-12, "expected exact equality at zero offset");
            } else {
                // Measured ~2.9e-5 (0.0029%) for a 1mm offset on this
                // geometry - a real, physically genuine effect (TDC and
                // BDC each get an independent phase correction under
                // offset, so there's no reason their travel distance
                // should exactly equal 2*crank_radius), not numerical
                // error. 1e-3 leaves real margin above the measured value
                // while still catching a much larger regression.
                assert!(relative_deviation < 1e-3, "offset={offset}: travel deviated from nominal stroke by {relative_deviation:e}");
            }
        }
    }

    #[test]
    fn position_is_symmetric_and_velocity_antisymmetric_about_tdc_with_no_offset() {
        let mechanism = s54b32_mechanism();
        for theta_deg in [20.0_f64, 60.0, 130.0] {
            let theta = theta_deg.to_radians();
            let position_positive = mechanism.piston_position_from_crank_center(theta);
            let position_negative = mechanism.piston_position_from_crank_center(-theta);
            assert!((position_positive - position_negative).abs() < 1e-12, "theta={theta_deg}: position should be symmetric");

            let velocity_positive = mechanism.piston_velocity(theta, 500.0);
            let velocity_negative = mechanism.piston_velocity(-theta, 500.0);
            assert!(
                (velocity_positive + velocity_negative).abs() < 1e-9,
                "theta={theta_deg}: velocity should be antisymmetric"
            );
        }
    }

    #[test]
    fn bottom_dead_center_angle_is_exactly_pi_with_no_offset() {
        let mechanism = s54b32_mechanism();
        let bdc_angle = mechanism.crank_angle_of_bottom_dead_center();
        assert!((bdc_angle - PI).abs() < 1e-12, "expected exactly pi, got {bdc_angle}");
    }

    #[test]
    fn zero_crank_radius_does_not_produce_nan() {
        // A zero-stroke mechanism is a legitimate way to model a rigid,
        // non-moving 0D element by reusing this kinematics code (position
        // is constant everywhere) - the Newton solve for true TDC must not
        // divide 0.0/0.0 in this degenerate case (both derivatives are
        // identically zero when nothing moves), which would otherwise
        // poison every subsequent position/velocity/acceleration query
        // with NaN.
        let mechanism = CrankMechanism::new(0.0, 0.139, 0.0);
        assert!(!mechanism.piston_position_from_crank_center(0.5).is_nan());
        assert!(!mechanism.piston_displacement_from_top_dead_center(0.5).is_nan());
        assert_eq!(mechanism.piston_displacement_from_top_dead_center(0.5), 0.0, "a rigid mechanism should show exactly zero displacement everywhere");
        assert!(!mechanism.piston_velocity(0.5, 500.0).is_nan());
        assert_eq!(mechanism.piston_velocity(0.5, 500.0), 0.0);
        assert!(!mechanism.crank_angle_of_bottom_dead_center().is_nan());
    }
}
