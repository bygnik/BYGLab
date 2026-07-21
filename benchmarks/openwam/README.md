# OpenWAM reference benchmarks

This directory builds [OpenWAM](https://github.com/CMT-UPV/OpenWAM) (CMT-Motores
Térmicos / Universitat Politècnica de València's open-source 1D gas-dynamic
engine code) and runs it to produce **reference numerical output** that the
Rust finite-volume solver will be checked against as it's developed.

OpenWAM (`OpenWAM/`, a git submodule) is GPLv3. We only ever consume its
*numerical output* here (pressure/temperature/velocity fields) — nothing from
its source is copied into our own code, so this has no license implications
for the Rust project.

## Building

```bash
docker build -t openwam:2.2.0 -f benchmarks/openwam/Dockerfile benchmarks/openwam
```

Builds cleanly on Ubuntu 22.04 with stock `cmake`/`g++` (C++11, no external
deps). Only warnings are unchecked `fscanf` return values — no errors.

## Running a case

```bash
docker run --rm -v "<absolute path to a case dir>:/work" openwam:2.2.0 /work/<case>.txt
```

The container's entrypoint is the `OpenWAM` binary; it takes the input file
path as `argv[1]` and writes output files into the current directory (i.e.
into the mounted case dir).

**Windows/Git Bash note:** MSYS mangles `-v` and the `/work/...` argument as
if they were local paths. Prefix the command with `MSYS_NO_PATHCONV=1` and
use a real Windows-style host path (`D:/...`) on the left side of the `-v`
mount, e.g.:

```bash
MSYS_NO_PATHCONV=1 docker run --rm -v "D:/Desktop/BYGLab/benchmarks/openwam/cases/single_pipe:/work" openwam:2.2.0 /work/single_pipe.txt
```

## The input file format

There is no XML/JSON input format here — OpenWAM's CLI solver
(`TOpenWAM::ReadInputData` in `Source/TOpenWAM.cpp`) reads a legacy
whitespace-delimited ASCII format via a long fixed sequence of `fscanf`
calls, one field/line at a time, in an exact hardcoded order (general data →
engine block → pipes → valves → plenums → compressors → connections →
turbo axis → sensors → controllers → output selectors → DLL flag).

This format has **no example files anywhere** (upstream only ships
`.PCS` binary project files consumed by **WAMer**, a separate closed,
Windows-only Delphi GUI that isn't part of this repo — see "Unusable
examples" below). Every case here was built by tracing the `fscanf` call
sequence directly in the C++ source.

**Comment syntax exists but is fragile.** On Linux, `CleanLabelsX()`
(`TOpenWAM.cpp:349`) strips `<...>` before parsing, by toggling on `<` and
back off at the *next* `>` — character by character, no nesting awareness.
**Any stray `>` inside a comment (e.g. from `->`) closes the comment early**
and leaks trailing text into the token stream, silently misaligning every
`fscanf` call after it (this manifested as a `std::bad_alloc` crash deep in
`ReadGeneralData`, from a corrupted array-size read many fields later —
nasty to debug blind). Rule: **never put a literal `>` inside a `<...>`
comment.**

See [`cases/single_pipe/single_pipe.txt`](cases/single_pipe/single_pipe.txt)
for a fully-annotated minimal case (one straight pipe, closed at one end,
open to atmosphere at the other, no engine/valves/turbo) as a format
reference for building further cases.

**Pipe-to-pipe junctions (`TipoCC=6`, `nmPipesConnection`) auto-detect their
pipes but still read extra fields.** Two pipes sharing the same boundary-condition
node number are automatically linked (`TCCUnionEntreTubos::ReadBoundaryData`
scans all pipes' `FNodoIzq`/`FNodoDer` for a match — no pipe IDs needed in the
file). But the class *also* reads two more doubles right after the `TipoCC=6`
token: wall thickness and conductivity (`FEspesor`, `FConductividad`) —
easy to miss since they're not read in the main `ReadConnections` switch
statement but deep in the boundary class's own read method (same pattern as
the general/pipes/output sections: **always check whether a boundary class
has its own follow-on reads beyond the dispatch switch**). See
[`cases/acoustic_resonance/acoustic_resonance.txt`](cases/acoustic_resonance/acoustic_resonance.txt)
for a working two-pipe example.

## Output format and a time-column caveat

Two different output channels exist, on two different clocks:

- **`<case>INS.DAT` / `<case>AVG.DAT`** — global instantaneous/average
  results, gated by the "results selector" block at the end of the input
  file. Their `Time` column is `AcumulatedTime` in real seconds, output
  every `agincr` seconds (confirmed: with `agincr=0.0005`,
  `SimulationDuration=0.05`, `single_pipeINS.DAT` has exactly 100 rows).
  Reliable, but this run requests no fields in it (all selectors set to 0),
  so it's just a bare time column.

- **`<case>INS_pre.DAT` / `_tem.DAT` / `_vel.DAT`** — full per-cell
  space-time profiles for a chosen pipe (requested via the
  `FNumMagnitudesEspTemp` block), one row per *physics timestep* (not per
  `agincr`), one column per mesh cell. This is the data that actually
  matters for benchmarking a finite-volume solver cell-by-cell.

  **The leading time value on each row of these files is unreliable.** It's
  computed from crank-angle bookkeeping (`Theta`, `RegimenFicticio` — a
  "fictitious" engine speed) that OpenWAM's architecture assumes even for a
  no-engine, pipe-only run (`TOutputResults.cpp:1720-1725`). With
  `EngineBlock=false` this produces a garbage, non-monotonic sequence (it
  visibly wraps, e.g. `..., 719.509, 719.802, 0.0965288, ...` at the tail of
  `single_pipeINS_pre.DAT`) — do not use it.

  **Validated workaround:** reconstruct time as
  `row_index * (SimulationDuration / total_data_rows)`. This assumes
  ~uniform timestep, which held for this case (explicit CFL-limited
  scheme, mild transient). Sanity check performed: the closed-end cell
  (column 1) first visibly departs from its initial 2 bar at data-row 163;
  with 2648 total rows that reconstructs to ≈3.08 ms, versus the analytic
  rarefaction-wave transit time for a 1 m tube at 20°C
  (`L / sqrt(γRT) = 1.0 / 343.2 m/s ≈ 2.91 ms`) — a ~6% match, consistent
  with the coarse 101-cell mesh. Treat this reconstruction as good enough
  to unblock benchmarking, not as exact ground truth; if a future case
  needs high-precision timing, prefer adding a minimal dummy engine block
  so the crank-angle clock is meaningful, or patch `DoSpaceTimeFiles` to
  print `AcumulatedTime` directly.

## Cases

Cases are grouped by what validates them. `analysis/` holds Python scripts
that compute the closed-form reference and compare it against each case's
OpenWAM output — run them after regenerating a case to reproduce the numbers
quoted below. They need `numpy` (already available via `.venv`).

### Analytically-validated cases

These have exact or closed-form solutions and are checked against them
directly — the strongest form of validation available, and the primary
target once the Rust solver exists (same scripts, swap the data source).

#### `cases/quiescent/` — trivial consistency check

1 m pipe, **closed at both ends**, uniform initial state (1.5 bar/20°C, at
rest), no perturbation anywhere. Exact solution: nothing changes, ever —
this checks that the FV scheme doesn't spuriously generate motion from a
rest state (a basic well-posedness/consistency requirement any correct
Euler solver must satisfy exactly, no discretization error possible).
**Result: pressure stayed exactly 1.5000 bar in every cell for the whole
run; velocity stayed at ~7.6e-14 m/s (floating-point noise) instead of
exactly 0.** Passes.

```bash
MSYS_NO_PATHCONV=1 docker run --rm -v "$(pwd)/cases/quiescent:/work" openwam:2.2.0 /work/quiescent.txt
```

#### `cases/acoustic_resonance/` — linear acoustics, closed-form period

Two 1 m pipes joined at a shared node (`TipoCC=6`), closed at both outer
ends, small pressure step (1.05 bar vs 1.00 bar, both 20°C) — small enough
to stay in the linear-acoustics regime (no shock forms). A closed-closed
tube's fundamental resonance period `T = 2L/c₀` is exact regardless of the
excitation's shape (all harmonics are integer multiples of the same
fundamental), so this checks the solver's phase/dispersion accuracy over
many wave transits, not just one.

**Result** (`analysis/verify_acoustic_resonance.py`): theoretical period
`T = 2×2m / 343.2 m/s = 11.6549 ms`; measured period from pressure-plateau
threshold crossings at the closed end = **11.647–11.656 ms (≤0.07% error)**.
First wall-arrival time of the initial disturbance: 2.9349 ms observed vs
2.9137 ms theoretical (0.73% error, consistent with the 101-cell mesh
resolution).

```bash
MSYS_NO_PATHCONV=1 docker run --rm -v "$(pwd)/cases/acoustic_resonance:/work" openwam:2.2.0 /work/acoustic_resonance.txt
cd analysis && python verify_acoustic_resonance.py
```

#### `cases/sod_shock_tube/` — nonlinear Riemann problem, exact solution

Same two-pipe junction topology, strong pressure ratio (10 bar vs 1 bar,
both 20°C air) — a Sod-shock-tube-style problem, generating a genuine
shock, contact discontinuity, and rarefaction fan. Compared pointwise
against `analysis/exact_riemann.py`, a from-scratch implementation of
Toro's exact Riemann solver (self-tested against the textbook Sod test:
`p* = 0.30313`, `u* = 0.92745`, matches to 1e-4).

**Result** (`analysis/verify_sod_shock_tube.py`, compared at t≈0.893 ms,
safely before either wave reaches a closed end at 1.81/2.91 ms): exact
star region `p* = 2.848 bar`, `u* = 281.84 m/s` — OpenWAM's plateau matches
to 4 significant figures. Full-profile pressure error: RMS 0.068 bar,
max 0.62 bar (on a 1–10 bar range); velocity error: RMS 6.8 m/s, max
69.6 m/s (on a 0–282 m/s range). Errors are concentrated exactly at the
shock and rarefaction fronts — ~1–2 cells of numerical smearing, exactly
what's expected of a shock-capturing FV scheme resolving a true
discontinuity. This is the flagship validation case: it exercises shocks,
rarefactions, and a contact discontinuity all in one run with an exact
reference.

```bash
MSYS_NO_PATHCONV=1 docker run --rm -v "$(pwd)/cases/sod_shock_tube:/work" openwam:2.2.0 /work/sod_shock_tube.txt
cd analysis && python verify_sod_shock_tube.py
```

### `cases/single_pipe/` — done, physically sanity-checked (not exact)

1 m × Ø50 mm pipe, 101-cell mesh, TVD scheme, Courant=0.8, closed at the
left end, **open to atmosphere** (1 bar) at the right end, initialized at
rest at 2 bar/20°C. Superseded for rigorous validation by
`sod_shock_tube/` (an internal junction is a true Riemann problem; this
open-end BC involves a discharge-coefficient/quasi-steady model that isn't
exactly the ideal Riemann solution) but kept as the first working case and
because it exercises the `nmOpenEndAtmosphere` boundary instead of a
junction. Committed outputs: `run.log`, `single_pipeAVG.DAT`,
`single_pipeINS.DAT`, `single_pipeINS_pre/tem/vel.DAT`.

```bash
MSYS_NO_PATHCONV=1 docker run --rm -v "$(pwd)/cases/single_pipe:/work" openwam:2.2.0 /work/single_pipe.txt
```

### Real engine models — single-cylinder proxy — `cases/engine_s54_{2500,5000,7000}rpm/` — done, physically plausible

Single-cylinder 4-stroke NA engine using BMW S54B32-derived parameters
(bore 87mm, stroke 91mm, rod 139mm, CR 11:1 — an assumed placeholder at the
time this case was built, superseded by the accurate 11.5:1 spec used in
the full 6-cylinder case below; intake ⌀49.5mm/exhaust ⌀43.1mm,
area-matched single-valve equivalents of the real 2× intake ⌀35mm / 2×
exhaust ⌀30.5mm ports; Wiebe combustion, gasoline/stoichiometric, IVO 340°/
IVC 600°, EVO 140°/EVC 400°, all referenced to firing TDC=0°), run at three
steady RPM points identical except for `FRegimen`. No closed-form solution
exists for a real combustion cycle, so these are validated by internal
physical consistency instead — see `analysis/verify_engine_plausibility.py`.

**Results, all three RPM points:**

| RPM | Peak cylinder pressure | @ crank angle | Peak temperature | CR recovered from V(θ) |
|-----|------------------------|----------------|-------------------|--------------------------|
| 2500 | 62.45 bar | 15.1° ATDC | 2566.9 °C | 11.000 |
| 5000 | 63.57 bar | 15.3° ATDC | 2549.2 °C | 11.000 |
| 7000 | 63.91 bar | 15.1° ATDC | 2545.1 °C | 11.000 |

All physically sane for an 11:1 NA gasoline engine (peak pressure in the
60-90 bar range, combustion phasing at 10-20° ATDC is textbook), and the
small monotonic pressure/temperature rise with RPM is a real, expected
trend (less time for heat loss per cycle at higher speed), not noise. The
compression ratio is recovered exactly from the geometric volume trace at
all three RPM points, confirming `CalculaVolumen` is being driven correctly
by the crank kinematics regardless of engine speed. The cylinder trace
(`<case>INS.DAT`, using the reliable `AcumulatedTime`-based channel, not
the buggy per-cell space-time one) shows a complete, qualitatively correct
combustion → expansion → exhaust blowdown → intake stroke sequence.

**Known limitation:** all three cases go NaN deterministically at crank
angle ≈600° — exactly the computed intake-valve-closing angle
(IVO 340° + 260° duration). Confirmed independent of Courant number (tried
0.8 and 0.4, same failure point) and of the valve lift table's endpoint
value (tried exact 0.0 and a 1e-6 m epsilon, no change) — this rules out a
CFL/timestep instability and points to a genuine numerical singularity.
Localized (not fully root-caused) to `TCilindro4T.cpp`'s `FCicloCerrado`
("closed cycle") initialization block, which runs once per cycle exactly
when crank angle crosses the intake valve's closing angle — the analogous
exhaust-valve-closing event (EVC≈400°) causes no issue, so it's specific to
whatever "start of closed cycle" bookkeeping triggers on IVC, not valve
closure in general. Valid data covers ~600 of 720 degrees (everything
except the tail of compression and the next combustion event) — enough for
the physical-plausibility validation above, but not a clean multi-cycle
run. Worth a focused follow-up before relying on this for anything beyond
sanity-checking.

**Other format notes learned building this case** (see the input files for
the fully-annotated field-by-field sequence):
- Section order matters and does NOT match visual/logical grouping:
  `ReadGeneralData → ReadEngine → ReadPipes → ReadValves → ...` — the
  engine block must come **before** the pipes block in the file, easy to
  get backwards.
- For engine-block runs, `SimulationDuration` (in `ReadGeneralData`) is
  reinterpreted as **number of cycles**, not seconds (`thmax =
  SimulationDuration * 720`) — using a small seconds-like value (e.g.
  0.048) terminates the run after a fraction of one crank degree.
- `TOutputResults::OutputAverageResults` dereferences `AvgEngine`
  unconditionally near the top of the function, before its own later
  null-check — if `NumEnginesAvg=0` in the AVG output selectors,
  `AvgEngine` stays `NULL` and this crashes. Must request at least one
  engine-average variable to avoid it, even if you don't care about the
  values.
- The per-cylinder instantaneous results channel (`NumCylindersIns`
  in `ReadInstantaneousResults` → `TCilindro::ReadInstantaneousResultsCilindro`)
  writes into the same reliable `AcumulatedTime`-based `<case>INS.DAT` file
  used by pipe-only cases — this is the channel to use for a P-V-T trace,
  not the per-cell space-time dump (which, in these engine-block runs,
  produced no rows at all — not investigated further since `INS.DAT` fully
  covers the validation need).
- A pipe-to-cylinder valve boundary condition (`TipoCC=7`/`8`) binds to
  exactly one pipe, so modeling a real engine's 2 intake + 2 exhaust valves
  requires either a manifold/junction network upstream of a single
  equivalent valve (what we did — area-matched single valve per port) or
  multiple separate intake/exhaust pipes merged at a plenum.

### Real engine models — full 6-cylinder — `cases/engine_s54b32_6cyl_4900rpm/` — done, physically plausible

The complete BMW S54B32 straight-6, not the single-cylinder proxy above:
6 cylinders, correct firing order, individual throttle-body-style intake
runners (matching the real S54's ITB design), and a proper 6-2-1 exhaust
header — built with accurate specs confirmed via web search rather than
the single-cylinder proxy's placeholder values: bore 87mm, stroke 91mm,
rod 139mm, **CR 11.5:1** (real spec, not the 11.0 guess used above),
firing order **1-5-3-6-2-4**, displacement checks out to 3246cc exactly
matching the real engine. Run at 4900 RPM (the engine's real torque-peak
operating point). Input file is generated programmatically by
[`generate_input.py`](cases/engine_s54b32_6cyl_4900rpm/generate_input.py)
(15 pipes, 22 connections, 12 valve definitions — hand-typing this
reliably wasn't realistic).

**Topology:**
- **Intake**: 6 independent runners (⌀45mm × 0.3m), each straight to
  atmosphere — real ITBs still share an airbox upstream of the throttles,
  but atmosphere is already an effectively-infinite constant-pressure
  reservoir, so modeling the airbox as a separate plenum was judged not
  worth the added complexity for this pass.
- **Exhaust**: a genuine 6-2-1 header using `TCCRamificacion`
  (`TipoCC=12`, `nmBranch`) — an uncapped N-way junction discovered while
  researching this case (distinct from the 2-pipe-only `TipoCC=6` used in
  the analytical cases). Cylinders grouped by firing-order adjacency for
  even pulse spacing within each group (a real 6-2-1 tuning heuristic, not
  a solver constraint): group A = cyl {1,3,2} (fire at 0°/240°/480°, 240°
  apart), group B = cyl {5,6,4} (120°/360°/600°, also 240° apart). Each
  group's 3 runners (⌀40mm × 0.35m) merge at one branch into a secondary
  pipe (⌀55mm × 0.3m); the two secondaries merge at a third branch into a
  ⌀65mm × 0.5m tailpipe to atmosphere. `TCCRamificacion` needs only the
  bare `12` token per connection — no extra fields, and no valve model
  required (unlike the plenum alternative, which would need one per
  attached pipe) — see the README's format notes below.

**Results:** compression ratio recovers to exactly **11.500** independently
for all 6 cylinders (confirms per-cylinder crank kinematics is correct
across the whole engine, not just cylinder 1). Cylinders 1, 3, and 5 (the
first three in the firing order) show complete, clean combustion events
within the captured window, each peaking at essentially the same crank
angle after their own TDC (14.88–15.01° ATDC — consistent with the
single-cylinder result and with each other), and — this is the real
validation of the firing-order/phase-offset implementation — **spaced
119.92–119.95° apart**, matching the expected 120° firing interval for a
6-cylinder 4-stroke engine to within 0.08°. See
`analysis/verify_6cyl_plausibility.py`.

```bash
MSYS_NO_PATHCONV=1 docker run --rm -v "$(pwd)/cases/engine_s54b32_6cyl_4900rpm:/work" openwam:2.2.0 /work/engine_s54b32_6cyl_4900rpm.txt
cd analysis && python verify_6cyl_plausibility.py
```

**Known limitation — same IVC bug, now independently confirmed 5 times
over:** every cylinder's data goes NaN at exactly **its own local crank
angle ≈600°** (the shared intake-valve-closing angle, each cylinder just
reaching it at a different *global* angle because of its firing-order
phase offset — cyl1 at global 600.2°, cyl6 at global 240.5°, cyl4 at
global 480.2°, cyl2 at global 360.4°, cyl5 at global 0.2° (wrapped from
720.2°) — every single one lands within 0.5° of local angle 600 once the
phase offset is subtracted out). This is strong independent confirmation
that the single-cylinder case's finding is a genuine, general OpenWAM bug
tied to IVC bookkeeping, not an artifact of that simpler topology.
Cylinders 2, 4, and 6 (the second half of the firing order) hit their own
IVC — and go NaN — before their combustion peak falls inside this
particular run's captured output window, so their peak-pressure numbers
in the raw data are not meaningful (the verification script filters these
out).

Also notable: cylinders 1, 3, and 5 — despite identical geometry, valve
timing, and combustion parameters — show meaningfully different peak
pressures (65.8, 81.2, 76.8 bar) and temperatures (2547, 2599, 2599 °C).
This is very unlikely to be numerical noise (their ATDC timing all agrees
to within 0.13° of each other) — plausible real cause is gas-dynamic
cross-talk through the shared exhaust manifold: each cylinder's exhaust
runner receives reflected pressure waves from its manifold-mates at a
different phase relative to its own cycle, changing effective trapped
mass/scavenging cylinder-to-cylinder. This is a genuine, well-known real
phenomenon in multi-cylinder engines (part of why exhaust manifold design
is a real tuning discipline) — a decoupled single-cylinder-times-6 model
would not reproduce it, so seeing it here is a point in favor of the
6-2-1 topology actually being modeled correctly, not a red flag.

**New format pieces learned building this case** (beyond the
single-cylinder engine notes above):
- The multi-cylinder firing-order block in `LeeMotor` (skipped entirely
  when `NCilin=1`) offers two modes: `tipodesfa=0` (enter each cylinder's
  phase offset directly) or `tipodesfa=1` (enter the firing order as a
  permutation of cylinder numbers; offsets are auto-computed as
  `position * 720/NCilin`). We used mode 1. There is also a **global**
  engine-level controllers-count field (distinct from the per-cylinder
  one) sandwiched between the firing-order block and the per-cylinder
  controllers loop — easy to miss since it looks like it should be part
  of one or the other.
- `TCCRamificacion` (`TipoCC=12`) vs `TCCUnionEntreTubos` (`TipoCC=6`):
  both auto-detect their pipes by matching node numbers, but the union
  type is hardcoded to exactly 2 pipes (a straight coupling, functionally
  one continuous duct) while the branch type accepts however many pipes
  reference its node — a real N-way junction with symmetric,
  direction-agnostic physics (any leg can be inflow or outflow depending
  on instantaneous conditions). Use branch for anything that's actually a
  merge/split, not union.
- `FMasaInicial`, wall temperatures, Wiebe combustion parameters, and
  essentially everything else in the engine block are **single shared
  values applied identically to all cylinders** — only the phase offset
  is genuinely per-cylinder. A more detailed model (different combustion
  quality per cylinder, individual wall temps, etc.) isn't representable
  without patching the solver.

### `cases/raw/*.zip` — unusable as-is

`2StrokeEngine_v2.1.zip`, `4S_HDSI_Turbo_v2.1.zip`, `6CylinderTwinT_v2.1.zip`
downloaded from the [SourceForge OpenWAM v2.1 Examples](https://sourceforge.net/projects/openwam/files/OpenWAM%20v2.1/Examples/)
page. Each unzips to a single `.PCS` file — WAMer's binary GUI project
format, not the plain-text format the CLI solver reads (confirmed by
inspecting the raw bytes: compressed/binary floats, not ASCII). WAMer
itself is a separate, Windows-only, closed-source(ish) Delphi GUI not
included in this repository, so these can't currently be converted.
Kept for reference in case a WAMer install or a `.PCS` parser becomes
available later; not part of the automated benchmark suite.

## Next steps

- Root-cause the IVC NaN (see "Known limitation" in the 6-cylinder
  section — now confirmed 8 times total across two independent
  topologies) to get clean full-run data for every cylinder, not just
  whichever ones happen to fire early enough in the captured window.
- A plenum/volume case and a standalone valve boundary condition case
  (short of a full engine block) would round out the "building blocks"
  library.
- Consider a small Python/Rust script to parse `_pre/_tem/_vel.DAT` and
  `INS.DAT` into a more convenient array format (e.g. `.npy`/`.csv`, with
  the row-index time reconstruction baked in where needed) for the future
  Rust test harness to load directly.
- This is the reference model the Rust port's 0D cylinder + valve +
  combustion + manifold-network implementation should ultimately be
  checked against.
