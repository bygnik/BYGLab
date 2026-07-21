"""Calibration probe for S54 power (run: python scripts/probe_s54.py)."""
from byglab_engine import S54B32, run_full_cycle

def run(**kw):
    d = {**S54B32.to_dict(), "n_cycles": 4, "n_cells": 28, "hist_stride": 40, **kw}
    r = run_full_cycle(d)
    print(
        f"hp={r.brake_hp_n:6.1f}  IMEP={r.imep_bar:5.2f}  BMEP={r.bmep_bar:5.2f}  "
        f"FMEP={r.fmep_bar:4.2f}  air={r.air_mg:6.1f}mg  Ppeak={r.Ppeak_bar:5.1f}  "
        f"mech={r.mech_eff:4.1f}%  kw={kw}"
    )
    return r

if __name__ == "__main__":
    run()
