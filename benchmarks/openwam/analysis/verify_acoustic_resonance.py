"""
Verifies OpenWAM's `acoustic_resonance` case against the linear-acoustics
prediction for a closed-closed resonant tube: fundamental period T = 2L/c0,
independent of the excitation shape.

Usage: python verify_acoustic_resonance.py
"""
import numpy as np

CASE_DIR = "../cases/acoustic_resonance"


def load_space_time(fname):
    with open(fname) as f:
        lines = f.readlines()
    data_lines = [l for l in lines[6:] if l.strip()]
    rows = [list(map(float, l.split())) for l in data_lines]
    return np.array(rows)


def threshold_crossings(sig, t, mid):
    down, up = [], []
    for i in range(1, len(sig)):
        if sig[i - 1] >= mid and sig[i] < mid:
            frac = (mid - sig[i - 1]) / (sig[i] - sig[i - 1])
            down.append(t[i - 1] + frac * (t[i] - t[i - 1]))
        if sig[i - 1] < mid and sig[i] >= mid:
            frac = (mid - sig[i - 1]) / (sig[i] - sig[i - 1])
            up.append(t[i - 1] + frac * (t[i] - t[i - 1]))
    return down, up


def main():
    pre = load_space_time(f"{CASE_DIR}/acoustic_resonanceINS_pre.DAT")
    n = len(pre)
    duration = 0.03
    dt = duration / n
    t = np.arange(n) * dt

    p_left = pre[:, 1]  # pipe1 cell0 = closed end at x=0 (BC#1)

    p_hi, p_lo = 1.05, 0.99983  # the two plateau pressures (bar)
    mid = (p_hi + p_lo) / 2
    down, up = threshold_crossings(p_left, t, mid)

    gamma, R, T0 = 1.4, 287.0, 293.15
    c0 = np.sqrt(gamma * R * T0)
    L = 2.0  # total closed-closed length (two 1m pipes)
    T_theory = 2 * L / c0

    print(f"c0 = {c0:.3f} m/s, theoretical fundamental period T = 2L/c0 = {T_theory*1000:.4f} ms")
    print(f"downward crossings (ms): {[round(x*1000,4) for x in down]}")
    print(f"upward crossings (ms):   {[round(x*1000,4) for x in up]}")

    for label, xs in [("downward", down), ("upward", up)]:
        if len(xs) > 1:
            diffs = np.diff(xs)
            mean_T = np.mean(diffs)
            err_pct = 100 * (mean_T - T_theory) / T_theory
            print(f"{label}: mean period = {mean_T*1000:.4f} ms  (error vs theory: {err_pct:+.3f}%)")

    L1 = 1.0
    t_transit_theory = L1 / c0
    print(f"\nfirst disturbance arrival at left wall: observed {down[0]*1000:.4f} ms, "
          f"theory L1/c0 = {t_transit_theory*1000:.4f} ms "
          f"(error: {100*(down[0]-t_transit_theory)/t_transit_theory:+.2f}%)")


if __name__ == "__main__":
    main()
