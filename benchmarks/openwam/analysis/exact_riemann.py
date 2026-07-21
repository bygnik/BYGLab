"""
Exact solution of the 1D Euler Riemann problem (Toro, "Riemann Solvers and
Numerical Methods for Fluid Dynamics", 3rd ed., Ch. 4).

Used to generate closed-form reference profiles for the OpenWAM two-pipe
Riemann-problem benchmark cases (benchmarks/openwam/cases/sod_shock_tube/,
acoustic_resonance/), and will be reused later to validate the Rust FV
solver directly against the same exact solution.

Usage as a library:
    sol = RiemannProblem(rho_l, u_l, p_l, rho_r, u_r, p_r, gamma=1.4)
    rho, u, p = sol.sample(x, t, x0=0.0)   # x, t, x0 in consistent SI units
"""
from __future__ import annotations
import numpy as np


class RiemannProblem:
    def __init__(self, rho_l, u_l, p_l, rho_r, u_r, p_r, gamma=1.4):
        self.rho_l, self.u_l, self.p_l = rho_l, u_l, p_l
        self.rho_r, self.u_r, self.p_r = rho_r, u_r, p_r
        self.gamma = gamma
        self.c_l = np.sqrt(gamma * p_l / rho_l)
        self.c_r = np.sqrt(gamma * p_r / rho_r)
        self._check_vacuum()
        self.p_star, self.u_star = self._solve_star_region()
        self.rho_star_l = self._star_density_left()
        self.rho_star_r = self._star_density_right()

    def _check_vacuum(self):
        g = self.gamma
        du_crit = 2 * self.c_l / (g - 1) + 2 * self.c_r / (g - 1)
        if du_crit <= (self.u_r - self.u_l):
            raise ValueError("Vacuum is generated: exact solver (this simplified "
                              "version) does not handle the vacuum case.")

    def _f_k(self, p, rho_k, p_k, c_k):
        g = self.gamma
        if p > p_k:  # shock
            A_k = 2.0 / ((g + 1) * rho_k)
            B_k = (g - 1) / (g + 1) * p_k
            return (p - p_k) * np.sqrt(A_k / (p + B_k))
        else:  # rarefaction
            return 2 * c_k / (g - 1) * ((p / p_k) ** ((g - 1) / (2 * g)) - 1)

    def _f_k_prime(self, p, rho_k, p_k, c_k):
        g = self.gamma
        if p > p_k:
            A_k = 2.0 / ((g + 1) * rho_k)
            B_k = (g - 1) / (g + 1) * p_k
            return np.sqrt(A_k / (B_k + p)) * (1 - (p - p_k) / (2 * (B_k + p)))
        else:
            return 1.0 / (rho_k * c_k) * (p / p_k) ** (-(g + 1) / (2 * g))

    def _solve_star_region(self, tol=1e-12, max_iter=100):
        g = self.gamma
        p_pv = max(tol, 0.5 * (self.p_l + self.p_r) - 0.125 * (self.u_r - self.u_l) *
                   (self.rho_l + self.rho_r) * (self.c_l + self.c_r))
        p = p_pv
        for _ in range(max_iter):
            f_l = self._f_k(p, self.rho_l, self.p_l, self.c_l)
            f_r = self._f_k(p, self.rho_r, self.p_r, self.c_r)
            f = f_l + f_r + (self.u_r - self.u_l)
            f_l_p = self._f_k_prime(p, self.rho_l, self.p_l, self.c_l)
            f_r_p = self._f_k_prime(p, self.rho_r, self.p_r, self.c_r)
            d = f_l_p + f_r_p
            p_new = p - f / d
            if p_new < tol:
                p_new = tol
            if abs(p_new - p) / (0.5 * (p_new + p)) < tol:
                p = p_new
                break
            p = p_new
        u = 0.5 * (self.u_l + self.u_r) + 0.5 * (
            self._f_k(p, self.rho_r, self.p_r, self.c_r) -
            self._f_k(p, self.rho_l, self.p_l, self.c_l))
        return p, u

    def _star_density_left(self):
        g = self.gamma
        if self.p_star > self.p_l:  # left shock
            ratio = self.p_star / self.p_l
            return self.rho_l * ((ratio + (g - 1) / (g + 1)) /
                                  ((g - 1) / (g + 1) * ratio + 1))
        else:  # left rarefaction
            return self.rho_l * (self.p_star / self.p_l) ** (1 / g)

    def _star_density_right(self):
        g = self.gamma
        if self.p_star > self.p_r:  # right shock
            ratio = self.p_star / self.p_r
            return self.rho_r * ((ratio + (g - 1) / (g + 1)) /
                                  ((g - 1) / (g + 1) * ratio + 1))
        else:  # right rarefaction
            return self.rho_r * (self.p_star / self.p_r) ** (1 / g)

    def wave_speeds(self):
        """Return a dict describing each wave's type and propagation speed(s)."""
        g = self.gamma
        out = {"p_star": self.p_star, "u_star": self.u_star}
        # left wave
        if self.p_star > self.p_l:
            q_l = np.sqrt((g + 1) / (2 * g) * (self.p_star / self.p_l) + (g - 1) / (2 * g))
            out["left"] = {"type": "shock", "speed": self.u_l - self.c_l * q_l}
        else:
            c_star_l = self.c_l * (self.p_star / self.p_l) ** ((g - 1) / (2 * g))
            out["left"] = {"type": "rarefaction",
                            "head_speed": self.u_l - self.c_l,
                            "tail_speed": self.u_star - c_star_l}
        # right wave
        if self.p_star > self.p_r:
            q_r = np.sqrt((g + 1) / (2 * g) * (self.p_star / self.p_r) + (g - 1) / (2 * g))
            out["right"] = {"type": "shock", "speed": self.u_r + self.c_r * q_r}
        else:
            c_star_r = self.c_r * (self.p_star / self.p_r) ** ((g - 1) / (2 * g))
            out["right"] = {"type": "rarefaction",
                             "head_speed": self.u_r + self.c_r,
                             "tail_speed": self.u_star + c_star_r}
        out["contact_speed"] = self.u_star
        return out

    def sample(self, x, t, x0=0.0):
        """Sample (rho, u, p) at position(s) x and time t > 0. x may be a numpy array."""
        g = self.gamma
        x = np.atleast_1d(np.asarray(x, dtype=float))
        s = (x - x0) / t  # similarity variable

        rho = np.empty_like(x)
        u = np.empty_like(x)
        p = np.empty_like(x)

        # left of contact
        left_of_contact = s <= self.u_star

        # --- left side ---
        if self.p_star > self.p_l:  # left shock
            s_l = self.u_l - self.c_l * np.sqrt((g + 1) / (2 * g) * (self.p_star / self.p_l) + (g - 1) / (2 * g))
            region_l_state = s < s_l
            region_star_l = left_of_contact & ~region_l_state
        else:  # left rarefaction
            c_star_l = self.c_l * (self.p_star / self.p_l) ** ((g - 1) / (2 * g))
            head = self.u_l - self.c_l
            tail = self.u_star - c_star_l
            region_l_state = s < head
            region_fan_l = left_of_contact & (s >= head) & (s < tail)
            region_star_l = left_of_contact & (s >= tail)

        # --- right side ---
        if self.p_star > self.p_r:  # right shock
            s_r = self.u_r + self.c_r * np.sqrt((g + 1) / (2 * g) * (self.p_star / self.p_r) + (g - 1) / (2 * g))
            region_r_state = s > s_r
            region_star_r = ~left_of_contact & ~region_r_state
        else:  # right rarefaction
            c_star_r = self.c_r * (self.p_star / self.p_r) ** ((g - 1) / (2 * g))
            head = self.u_r + self.c_r
            tail = self.u_star + c_star_r
            region_r_state = s > head
            region_fan_r = (~left_of_contact) & (s <= head) & (s > tail)
            region_star_r = (~left_of_contact) & (s <= tail)

        # fill left original state
        rho[region_l_state] = self.rho_l
        u[region_l_state] = self.u_l
        p[region_l_state] = self.p_l

        # fill left star state
        rho[region_star_l] = self.rho_star_l
        u[region_star_l] = self.u_star
        p[region_star_l] = self.p_star

        # fill left fan (if rarefaction)
        if self.p_star <= self.p_l:
            idx = region_fan_l
            u_fan = 2 / (g + 1) * (self.c_l + (g - 1) / 2 * self.u_l + s[idx])
            c_fan = 2 / (g + 1) * (self.c_l + (g - 1) / 2 * (self.u_l - s[idx]))
            rho[idx] = self.rho_l * (c_fan / self.c_l) ** (2 / (g - 1))
            u[idx] = u_fan
            p[idx] = self.p_l * (c_fan / self.c_l) ** (2 * g / (g - 1))

        # fill right original state
        rho[region_r_state] = self.rho_r
        u[region_r_state] = self.u_r
        p[region_r_state] = self.p_r

        # fill right star state
        rho[region_star_r] = self.rho_star_r
        u[region_star_r] = self.u_star
        p[region_star_r] = self.p_star

        # fill right fan (if rarefaction)
        if self.p_star <= self.p_r:
            idx = region_fan_r
            u_fan = 2 / (g + 1) * (-self.c_r + (g - 1) / 2 * self.u_r + s[idx])
            c_fan = 2 / (g + 1) * (self.c_r - (g - 1) / 2 * (self.u_r - s[idx]))
            rho[idx] = self.rho_r * (c_fan / self.c_r) ** (2 / (g - 1))
            u[idx] = u_fan
            p[idx] = self.p_r * (c_fan / self.c_r) ** (2 * g / (g - 1))

        return rho, u, p


if __name__ == "__main__":
    # Quick self-test against Toro's classic Sod test (Ch. 4, Test 1):
    # exact star pressure p* = 0.30313, u* = 0.92745
    sol = RiemannProblem(rho_l=1.0, u_l=0.0, p_l=1.0, rho_r=0.125, u_r=0.0, p_r=0.1, gamma=1.4)
    print("Sod test star region: p* =", sol.p_star, " u* =", sol.u_star)
    assert abs(sol.p_star - 0.30313) < 1e-4
    assert abs(sol.u_star - 0.92745) < 1e-4
    print("Self-test passed.")
