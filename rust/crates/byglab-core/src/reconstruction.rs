//! Second-order MUSCL (piecewise-linear) reconstruction with a slope
//! limiter, paired with the MUSCL-Hancock half-timestep evolution that
//! makes the resulting scheme second-order accurate in time as well as
//! space.
//!
//! [`crate::riemann::hllc_flux`] itself does not change at all to get this
//! upgrade — it only ever sees two states on either side of a face and has
//! no notion of whether those came straight from a cell center
//! (first-order) or from the more careful reconstruction in this module
//! (second-order). Upgrading accuracy is entirely a matter of computing
//! better face states to hand it.

use crate::gas::{Flux, GasProperties, PrimitiveState};

/// The slope-limited left and right face values for one cell, evolved
/// forward by half a timestep (Toro, "Riemann Solvers and Numerical
/// Methods for Fluid Dynamics", 3rd ed., section 14.4 — the MUSCL-Hancock
/// "predictor" step).
#[derive(Debug, Clone, Copy)]
pub struct ReconstructedCell {
    pub left_face: PrimitiveState,
    pub right_face: PrimitiveState,
}

/// Reconstructs one cell, given the (unreconstructed) state just before it
/// (`before`), its own cell-center state (`center`), and the state just
/// after it (`after`).
///
/// Three steps: (1) compute a limited slope from the neighboring
/// differences, (2) linearly extrapolate to the cell's two faces, (3)
/// evolve both face values forward by half a timestep using the flux
/// difference across this (still-unevolved) reconstructed profile. That
/// half-step evolution is what makes this genuinely second-order in time,
/// not just in space — the reason a generic "piecewise-linear + plain
/// forward Euler" scheme would NOT be fully second-order, but
/// MUSCL-Hancock specifically is, without needing a multi-stage outer
/// time integrator (see `solver.rs`'s doc comment).
pub fn reconstruct_cell(
    before: PrimitiveState,
    center: PrimitiveState,
    after: PrimitiveState,
    cell_width: f64,
    dt: f64,
    gas: &GasProperties,
) -> ReconstructedCell {
    let slope = limited_slope(before, center, after);
    let left_face = extrapolate_to_face(center, slope, -0.5);
    let right_face = extrapolate_to_face(center, slope, 0.5);

    let flux_at_left_face = left_face.to_conserved(gas).physical_flux(gas);
    let flux_at_right_face = right_face.to_conserved(gas).physical_flux(gas);
    let half_dt_over_dx = 0.5 * dt / cell_width;

    ReconstructedCell {
        left_face: evolve_half_step(left_face, flux_at_left_face, flux_at_right_face, half_dt_over_dx, gas),
        right_face: evolve_half_step(right_face, flux_at_left_face, flux_at_right_face, half_dt_over_dx, gas),
    }
}

/// The minmod slope limiter, applied to one scalar quantity: takes
/// whichever of the backward/forward differences has the smaller
/// magnitude, or zero if they disagree in sign (a local minimum or
/// maximum, where any nonzero slope would overshoot the neighboring
/// values).
///
/// This is the most diffusive (least sharp, but most robust) of the
/// standard TVD limiters — a deliberately conservative choice for a first
/// upgrade beyond first-order. Sharper limiters (van Leer, superbee) are a
/// natural future tuning step once this one is validated.
fn minmod(backward_difference: f64, forward_difference: f64) -> f64 {
    if backward_difference * forward_difference <= 0.0 {
        0.0
    } else if backward_difference.abs() < forward_difference.abs() {
        backward_difference
    } else {
        forward_difference
    }
}

/// Component-wise minmod slope for a full primitive state.
fn limited_slope(before: PrimitiveState, center: PrimitiveState, after: PrimitiveState) -> PrimitiveState {
    PrimitiveState {
        density: minmod(center.density - before.density, after.density - center.density),
        velocity: minmod(center.velocity - before.velocity, after.velocity - center.velocity),
        pressure: minmod(center.pressure - before.pressure, after.pressure - center.pressure),
    }
}

/// Linearly extrapolates from a cell center to one of its faces:
/// `center + fraction * slope`, where `fraction` is -0.5 for the left face
/// and +0.5 for the right face (the face sits half a cell width away from
/// the center).
fn extrapolate_to_face(center: PrimitiveState, slope: PrimitiveState, fraction: f64) -> PrimitiveState {
    PrimitiveState {
        density: center.density + fraction * slope.density,
        velocity: center.velocity + fraction * slope.velocity,
        pressure: center.pressure + fraction * slope.pressure,
    }
}

/// Advances one reconstructed face value forward by half a timestep,
/// using the flux difference across the cell's *unevolved* reconstructed
/// profile (`flux_at_right_face - flux_at_left_face`) — the same
/// correction is applied to both the left and right face (Toro eq. 14.20),
/// which is what keeps the reconstructed profile's slope unchanged while
/// shifting both faces to their values half a timestep later.
fn evolve_half_step(
    face: PrimitiveState,
    flux_at_left_face: Flux,
    flux_at_right_face: Flux,
    half_dt_over_dx: f64,
    gas: &GasProperties,
) -> PrimitiveState {
    let mut conserved = face.to_conserved(gas);
    conserved.mass -= half_dt_over_dx * (flux_at_right_face.mass - flux_at_left_face.mass);
    conserved.momentum -= half_dt_over_dx * (flux_at_right_face.momentum - flux_at_left_face.momentum);
    conserved.energy -= half_dt_over_dx * (flux_at_right_face.energy - flux_at_left_face.energy);
    conserved.to_primitive(gas)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minmod_picks_the_smaller_magnitude_when_signs_agree() {
        assert_eq!(minmod(2.0, 4.0), 2.0);
        assert_eq!(minmod(4.0, 2.0), 2.0);
        assert_eq!(minmod(-2.0, -4.0), -2.0);
    }

    #[test]
    fn minmod_is_zero_at_a_local_extremum() {
        assert_eq!(minmod(2.0, -1.0), 0.0);
        assert_eq!(minmod(-2.0, 1.0), 0.0);
        assert_eq!(minmod(0.0, 5.0), 0.0);
    }

    #[test]
    fn uniform_state_reconstructs_to_itself_with_zero_slope() {
        let gas = GasProperties::AIR;
        let state = PrimitiveState { density: 1.2, velocity: 10.0, pressure: 150_000.0 };
        let reconstructed = reconstruct_cell(state, state, state, 0.01, 1e-6, &gas);
        // Flux is identical everywhere for a uniform state, so the
        // half-step correction is also zero - both faces should land
        // exactly on the original state.
        assert!((reconstructed.left_face.pressure - state.pressure).abs() < 1e-6);
        assert!((reconstructed.right_face.pressure - state.pressure).abs() < 1e-6);
        assert!((reconstructed.left_face.velocity - state.velocity).abs() < 1e-9);
        assert!((reconstructed.right_face.velocity - state.velocity).abs() < 1e-9);
    }

    #[test]
    fn reconstruction_does_not_overshoot_beyond_neighboring_values() {
        // TVD property: the reconstructed face values must stay within the
        // range spanned by the cell and its neighbors - this is what keeps
        // MUSCL-Hancock from producing non-physical (e.g. negative
        // density/pressure) states even next to a strong discontinuity.
        let gas = GasProperties::AIR;
        let before = PrimitiveState { density: 1.0, velocity: 0.0, pressure: 100_000.0 };
        let center = PrimitiveState { density: 5.0, velocity: 0.0, pressure: 1_000_000.0 };
        let after = PrimitiveState { density: 1.0, velocity: 0.0, pressure: 100_000.0 };
        let reconstructed = reconstruct_cell(before, center, after, 0.01, 1e-7, &gas);
        assert!(reconstructed.left_face.pressure <= center.pressure);
        assert!(reconstructed.right_face.pressure <= center.pressure);
        assert!(reconstructed.left_face.pressure >= before.pressure.min(after.pressure));
    }
}
