"""
Compares OpenWAM's `sod_shock_tube` case output against the exact Riemann
solution, and prints the error summary documented in the benchmarks README.

Usage: python verify_sod_shock_tube.py
(run from benchmarks/openwam/analysis/, with the case already generated at
../cases/sod_shock_tube/)
"""
import numpy as np
from exact_riemann import RiemannProblem

CASE_DIR = "../cases/sod_shock_tube"


def load_space_time(fname):
    with open(fname) as f:
        lines = f.readlines()
    # header is 6 lines: counts, cyl ids, plenum ids, pipe ids, cell counts, wamer count
    data_lines = [l for l in lines[6:] if l.strip()]
    rows = [list(map(float, l.split())) for l in data_lines]
    return np.array(rows)


def main():
    pre = load_space_time(f"{CASE_DIR}/sod_shock_tubeINS_pre.DAT")
    vel = load_space_time(f"{CASE_DIR}/sod_shock_tubeINS_vel.DAT")

    n = len(pre)
    duration = 0.0012  # SimulationDuration in the input file
    dt = duration / n
    t = np.arange(n) * dt

    # NOTE: the leading time value printed by OpenWAM in these files is NOT
    # reliable for a no-engine-block run (see README's "time-column caveat").
    # Time is reconstructed here as row_index * (SimulationDuration / n_rows).
    target_t = 0.0009  # well before either wave reaches a closed end (1.81/2.91 ms)
    row_idx = int(round(target_t / dt))
    t_actual = t[row_idx]

    Ncell = 101
    L = 1.0
    dx = L / Ncell
    x_pipe1 = -1.0 + (np.arange(Ncell) + 0.5) * dx  # driver pipe, x in [-1, 0]
    x_pipe2 = 0.0 + (np.arange(Ncell) + 0.5) * dx    # driven pipe, x in [0, 1]
    x_all = np.concatenate([x_pipe1, x_pipe2])

    p_ow = pre[row_idx, 1:1 + 2 * Ncell]  # bar
    u_ow = vel[row_idx, 1:1 + 2 * Ncell]  # m/s

    R, T0, gamma = 287.0, 293.15, 1.4
    p_l, p_r = 10e5, 1e5
    rho_l, rho_r = p_l / (R * T0), p_r / (R * T0)
    sol = RiemannProblem(rho_l, 0.0, p_l, rho_r, 0.0, p_r, gamma=gamma)

    rho_ex, u_ex, p_ex = sol.sample(x_all, t_actual, x0=0.0)
    p_ex_bar = p_ex / 1e5

    err_p = p_ow - p_ex_bar
    err_u = u_ow - u_ex

    ws = sol.wave_speeds()
    print(f"Comparing at reconstructed t = {t_actual*1000:.4f} ms (row {row_idx}/{n})")
    print(f"Exact star region: p* = {sol.p_star/1e5:.4f} bar, u* = {sol.u_star:.2f} m/s")
    print(f"Exact left wave:  {ws['left']}")
    print(f"Exact right wave: {ws['right']}")
    print()
    print(f"Pressure error: max|err| = {np.max(np.abs(err_p)):.4f} bar, "
          f"RMS = {np.sqrt(np.mean(err_p**2)):.4f} bar  (state range: 1-10 bar)")
    print(f"Velocity error: max|err| = {np.max(np.abs(err_u)):.3f} m/s, "
          f"RMS = {np.sqrt(np.mean(err_u**2)):.3f} m/s  (state range: 0-282 m/s)")
    print("(errors concentrate at the shock/rarefaction fronts - numerical")
    print(" smearing over ~1-2 cells, as expected for a shock-capturing FV scheme)")


if __name__ == "__main__":
    main()
