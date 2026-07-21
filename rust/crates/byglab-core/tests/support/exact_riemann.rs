//! Exact solution of the 1D Euler Riemann problem (Toro, "Riemann Solvers
//! and Numerical Methods for Fluid Dynamics", 3rd ed., Ch. 4).
//!
//! Test-only support code: a direct Rust port of the already-validated
//! `benchmarks/openwam/analysis/exact_riemann.py`, used to check
//! `byglab-core`'s `sod_shock_tube` result against a genuine closed-form
//! answer (not just OpenWAM's own numerical approximation). Ported
//! line-for-line rather than restructured, so it stays easy to compare
//! against the original if either ever needs a second look.
//!
//! Only handles the non-vacuum case, matching the Python original.

/// The exact solution of a 1D Riemann problem: a discontinuity at `x0`
/// between a left state `(rho_l, u_l, p_l)` and a right state
/// `(rho_r, u_r, p_r)`, both in the same ideal gas (`gamma`).
pub struct RiemannProblem {
    rho_l: f64,
    u_l: f64,
    p_l: f64,
    c_l: f64,
    rho_r: f64,
    u_r: f64,
    p_r: f64,
    c_r: f64,
    gamma: f64,
    pub p_star: f64,
    pub u_star: f64,
    rho_star_l: f64,
    rho_star_r: f64,
}

impl RiemannProblem {
    /// Solves the Riemann problem. Panics if the result would be a vacuum
    /// (two rarefactions that don't overlap) — this simplified exact
    /// solver doesn't handle that case, matching the Python original's
    /// behavior of raising rather than returning a bogus answer.
    #[allow(clippy::too_many_arguments)]
    pub fn new(rho_l: f64, u_l: f64, p_l: f64, rho_r: f64, u_r: f64, p_r: f64, gamma: f64) -> Self {
        let c_l = (gamma * p_l / rho_l).sqrt();
        let c_r = (gamma * p_r / rho_r).sqrt();

        let du_crit = 2.0 * c_l / (gamma - 1.0) + 2.0 * c_r / (gamma - 1.0);
        assert!(
            du_crit > (u_r - u_l),
            "vacuum generated: this simplified exact solver does not handle the vacuum case"
        );

        let (p_star, u_star) = Self::solve_star_region(rho_l, u_l, p_l, c_l, rho_r, u_r, p_r, c_r, gamma);
        let rho_star_l = Self::star_density(rho_l, p_l, p_star, gamma);
        let rho_star_r = Self::star_density(rho_r, p_r, p_star, gamma);

        RiemannProblem { rho_l, u_l, p_l, c_l, rho_r, u_r, p_r, c_r, gamma, p_star, u_star, rho_star_l, rho_star_r }
    }

    /// Sampling function for the Newton iteration on the star-region
    /// pressure (Toro eq. 4.6/4.7): the velocity jump across one side's
    /// wave, as a function of a trial star pressure `p`.
    fn f_k(p: f64, rho_k: f64, p_k: f64, c_k: f64, gamma: f64) -> f64 {
        if p > p_k {
            // shock
            let a_k = 2.0 / ((gamma + 1.0) * rho_k);
            let b_k = (gamma - 1.0) / (gamma + 1.0) * p_k;
            (p - p_k) * (a_k / (p + b_k)).sqrt()
        } else {
            // rarefaction
            2.0 * c_k / (gamma - 1.0) * ((p / p_k).powf((gamma - 1.0) / (2.0 * gamma)) - 1.0)
        }
    }

    /// Derivative of [`Self::f_k`] with respect to `p`, used by the Newton
    /// iteration.
    fn f_k_prime(p: f64, rho_k: f64, p_k: f64, c_k: f64, gamma: f64) -> f64 {
        if p > p_k {
            let a_k = 2.0 / ((gamma + 1.0) * rho_k);
            let b_k = (gamma - 1.0) / (gamma + 1.0) * p_k;
            (a_k / (b_k + p)).sqrt() * (1.0 - (p - p_k) / (2.0 * (b_k + p)))
        } else {
            1.0 / (rho_k * c_k) * (p / p_k).powf(-(gamma + 1.0) / (2.0 * gamma))
        }
    }

    /// Newton-Raphson iteration for the star-region pressure and velocity
    /// (Toro section 4.3), starting from the "primitive variable Riemann
    /// solver" (PVRS) estimate.
    #[allow(clippy::too_many_arguments)]
    fn solve_star_region(rho_l: f64, u_l: f64, p_l: f64, c_l: f64, rho_r: f64, u_r: f64, p_r: f64, c_r: f64, gamma: f64) -> (f64, f64) {
        let tolerance = 1e-12;
        let max_iterations = 100;

        let mut p = (0.5 * (p_l + p_r) - 0.125 * (u_r - u_l) * (rho_l + rho_r) * (c_l + c_r)).max(tolerance);

        for _ in 0..max_iterations {
            let f_l = Self::f_k(p, rho_l, p_l, c_l, gamma);
            let f_r = Self::f_k(p, rho_r, p_r, c_r, gamma);
            let f = f_l + f_r + (u_r - u_l);
            let f_l_prime = Self::f_k_prime(p, rho_l, p_l, c_l, gamma);
            let f_r_prime = Self::f_k_prime(p, rho_r, p_r, c_r, gamma);

            let mut p_new = p - f / (f_l_prime + f_r_prime);
            if p_new < tolerance {
                p_new = tolerance;
            }
            let converged = (p_new - p).abs() / (0.5 * (p_new + p)) < tolerance;
            p = p_new;
            if converged {
                break;
            }
        }

        let u_star = 0.5 * (u_l + u_r) + 0.5 * (Self::f_k(p, rho_r, p_r, c_r, gamma) - Self::f_k(p, rho_l, p_l, c_l, gamma));
        (p, u_star)
    }

    /// Density in the star region on one side, given that side's original
    /// state and the solved star pressure.
    fn star_density(rho_k: f64, p_k: f64, p_star: f64, gamma: f64) -> f64 {
        if p_star > p_k {
            // shock
            let ratio = p_star / p_k;
            rho_k * ((ratio + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * ratio + 1.0))
        } else {
            // rarefaction (isentropic)
            rho_k * (p_star / p_k).powf(1.0 / gamma)
        }
    }

    /// Samples the exact solution `(density, velocity, pressure)` at
    /// position `x` and time `t > 0`, given the initial discontinuity was
    /// at `x0`. Uses the self-similar structure of the Riemann problem:
    /// the whole solution only depends on `x` and `t` through the ratio
    /// `s = (x - x0) / t`.
    pub fn sample(&self, x: f64, t: f64, x0: f64) -> (f64, f64, f64) {
        let g = self.gamma;
        let s = (x - x0) / t;

        if s <= self.u_star {
            self.sample_left_of_contact(s)
        } else {
            self.sample_right_of_contact(s, g)
        }
    }

    fn sample_left_of_contact(&self, s: f64) -> (f64, f64, f64) {
        let g = self.gamma;
        if self.p_star > self.p_l {
            // left-moving shock
            let shock_speed =
                self.u_l - self.c_l * ((g + 1.0) / (2.0 * g) * (self.p_star / self.p_l) + (g - 1.0) / (2.0 * g)).sqrt();
            if s < shock_speed {
                (self.rho_l, self.u_l, self.p_l)
            } else {
                (self.rho_star_l, self.u_star, self.p_star)
            }
        } else {
            // left-moving rarefaction fan
            let c_star_l = self.c_l * (self.p_star / self.p_l).powf((g - 1.0) / (2.0 * g));
            let head = self.u_l - self.c_l;
            let tail = self.u_star - c_star_l;
            if s < head {
                (self.rho_l, self.u_l, self.p_l)
            } else if s < tail {
                let velocity = 2.0 / (g + 1.0) * (self.c_l + (g - 1.0) / 2.0 * self.u_l + s);
                let sound_speed = 2.0 / (g + 1.0) * (self.c_l + (g - 1.0) / 2.0 * (self.u_l - s));
                let density = self.rho_l * (sound_speed / self.c_l).powf(2.0 / (g - 1.0));
                let pressure = self.p_l * (sound_speed / self.c_l).powf(2.0 * g / (g - 1.0));
                (density, velocity, pressure)
            } else {
                (self.rho_star_l, self.u_star, self.p_star)
            }
        }
    }

    fn sample_right_of_contact(&self, s: f64, g: f64) -> (f64, f64, f64) {
        if self.p_star > self.p_r {
            // right-moving shock
            let shock_speed =
                self.u_r + self.c_r * ((g + 1.0) / (2.0 * g) * (self.p_star / self.p_r) + (g - 1.0) / (2.0 * g)).sqrt();
            if s > shock_speed {
                (self.rho_r, self.u_r, self.p_r)
            } else {
                (self.rho_star_r, self.u_star, self.p_star)
            }
        } else {
            // right-moving rarefaction fan
            let c_star_r = self.c_r * (self.p_star / self.p_r).powf((g - 1.0) / (2.0 * g));
            let head = self.u_r + self.c_r;
            let tail = self.u_star + c_star_r;
            if s > head {
                (self.rho_r, self.u_r, self.p_r)
            } else if s > tail {
                let velocity = 2.0 / (g + 1.0) * (-self.c_r + (g - 1.0) / 2.0 * self.u_r + s);
                let sound_speed = 2.0 / (g + 1.0) * (self.c_r - (g - 1.0) / 2.0 * (self.u_r - s));
                let density = self.rho_r * (sound_speed / self.c_r).powf(2.0 / (g - 1.0));
                let pressure = self.p_r * (sound_speed / self.c_r).powf(2.0 * g / (g - 1.0));
                (density, velocity, pressure)
            } else {
                (self.rho_star_r, self.u_star, self.p_star)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Toro's classic Sod test (Ch. 4, Test 1): exact star pressure
    /// p* = 0.30313, u* = 0.92745. Mirrors the self-test already passing
    /// in `benchmarks/openwam/analysis/exact_riemann.py` — confirms this
    /// port is correct before trusting it to validate the real solver.
    #[test]
    fn matches_toros_classic_sod_test() {
        let solution = RiemannProblem::new(1.0, 0.0, 1.0, 0.125, 0.0, 0.1, 1.4);
        assert!((solution.p_star - 0.30313).abs() < 1e-4, "p* = {}", solution.p_star);
        assert!((solution.u_star - 0.92745).abs() < 1e-4, "u* = {}", solution.u_star);
    }
}
