"""
Physical-plausibility check for the engine_s54_*rpm cases: no closed-form
solution exists for a real combustion cycle, so this checks internal
consistency instead (compression ratio recovered from the volume trace,
peak pressure/temperature in a sane range and at a sane crank angle,
consistent behavior across RPM points) per benchmarks/openwam/README.md.

Usage: python verify_engine_plausibility.py
"""
import numpy as np

RPMS = [2500, 5000, 7000]


def load_cylinder_trace(rpm):
    fname = f"../cases/engine_s54_{rpm}rpm/engine_s54_{rpm}rpmINS.DAT"
    with open(fname) as f:
        lines = f.readlines()
    data = [list(map(float, l.split())) for l in lines[1:] if l.strip()]
    return np.array(data)  # columns: Time, Angle(deg), Pressure(bar), Temperature(degC), Volume(m^3)


def main():
    print(f"{'RPM':>6} {'valid/total':>12} {'angle range':>18} {'peak P(bar)':>12} "
          f"{'@deg':>8} {'peak T(degC)':>13} {'CR check':>10}")
    for rpm in RPMS:
        arr = load_cylinder_trace(rpm)
        t, ang, p, T, V = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4]
        nan_idx = np.where(np.isnan(p))[0]
        n_valid = nan_idx[0] if len(nan_idx) else len(arr)
        p_v, T_v, ang_v, V_v = p[:n_valid], T[:n_valid], ang[:n_valid], V[:n_valid]

        imax = np.nanargmax(p_v)
        cr_check = V_v.max() / V_v.min()

        print(f"{rpm:>6} {n_valid:>6}/{len(arr):<5} {ang_v[0]:>7.2f}-{ang_v[-1]:<7.2f} "
              f"{p_v[imax]:>12.2f} {ang_v[imax]:>8.2f} {np.nanmax(T_v):>13.1f} {cr_check:>10.3f}")

    print("\nExpected: CR check ~= 11.000 (set via FGeom.RelaCompresion), peak pressure in the")
    print("60-90 bar range typical of an 11:1 NA gasoline engine, peak pressure at roughly")
    print("10-20 deg ATDC (textbook SI combustion phasing for the Wiebe parameters used).")
    print("\nKnown limitation: all three cases go NaN deterministically at ~600 deg crank angle")
    print("(exactly the intake-valve-closing angle, FAnguloApertura=340 + duration 260 = 600),")
    print("independent of Courant number and valve-lift-table endpoint epsilon - a genuine")
    print("numerical singularity in the 'closed cycle' initialization triggered at IVC, not yet")
    print("root-caused beyond localizing it to TCilindro4T.cpp's FCicloCerrado block. Valid data")
    print("covers combustion, expansion, exhaust blowdown, and intake stroke (~600 of 720 deg).")


if __name__ == "__main__":
    main()
