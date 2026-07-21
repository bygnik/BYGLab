import { useCallback, useEffect, useMemo, useState } from "react";
import {
  LineChart, Line, Scatter, LabelList, XAxis, YAxis, CartesianGrid, Tooltip,
  ReferenceLine, ResponsiveContainer, Legend,
} from "recharts";
import { fetchPresetS54, simulateWave } from "../api";

const T = {
  void: "#0A0D10", panel: "#12161B", panel2: "#171C22", line: "#232A32", hair: "#2E363F",
  cyan: "#4FD8C4", amber: "#F2A93B", red: "#E5644B", green: "#7FD858", hi: "#E7ECEF", lo: "#7C8791",
};

function Field({ label, unit, value, onChange, min, max, step }) {
  return (
    <div className="mb-4">
      <div className="flex items-baseline justify-between mb-1">
        <label className="text-xs tracking-wide uppercase text-lo">{label}</label>
        <span className="text-sm font-mono text-cyan">{value}{unit}</span>
      </div>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        className="w-full"
      />
    </div>
  );
}

function Stat({ label, value, sub, color }) {
  return (
    <div className="px-4 py-3 rounded-md bg-panel2 border border-hair">
      <div className="text-xs md:text-sm uppercase tracking-wider mb-1 text-lo">{label}</div>
      <div className="text-2xl font-mono font-semibold" style={{ color: color || T.hi }}>{value}</div>
      {sub && <div className="text-xs md:text-sm mt-0.5 font-mono text-lo">{sub}</div>}
    </div>
  );
}

function Tab({ active, onClick, children }) {
  return (
    <button
      type="button"
      onClick={onClick}
      className="px-3 py-1.5 text-xs uppercase tracking-wider rounded-t-md border border-hair"
      style={{
        color: active ? T.void : T.lo,
        background: active ? T.cyan : "transparent",
        borderBottom: active ? "none" : `1px solid ${T.hair}`,
      }}
    >
      {children}
    </button>
  );
}

const DEFAULTS = {
  bore: 87, stroke: 91, rod: 139, rpm: 8000, CR: 11.5, intake_std: "j1349", n_cylinders: 6,
  l_int: 0.35, d_int: 46, l_exh: 0.50, d_exh: 42,
  afr: 12.8, eta_comb: 0.97, wall_heat_frac: 0.12, soc: 340, burn_dur: 55, n_cycles: 5,
  fmep_a: 0.28, fmep_b: 0.004, fmep_c: 0.060, fmep_d: 0.0006,
  i_valve_dia: 35, i_throat_pct: 85, i_stem_dia: 6, i_n_valves: 2, i_seat: 45, i_cd: 0.75,
  i_dur006: 260, i_dur040: 218, i_max_lift: 12, i_centerline: 125,
  e_valve_dia: 30.5, e_throat_pct: 85, e_stem_dia: 6, e_n_valves: 2, e_seat: 45, e_cd: 0.75,
  e_dur006: 260, e_dur040: 218, e_max_lift: 12, e_centerline: -130,
  n_cells: 36, hist_stride: 10, cfl: 0.45,
};

function computeCylinderVolumeL(boreMm, strokeMm, rodMm, compressionRatio, thetaDeg) {
  const boreM = boreMm / 1000;
  const strokeM = strokeMm / 1000;
  const rodM = rodMm / 1000;
  const areaM2 = (Math.PI / 4) * boreM * boreM;
  const sweptVolM3 = areaM2 * strokeM;
  const clearanceVolM3 = sweptVolM3 / (compressionRatio - 1);
  const thetaRad = ((thetaDeg % 360) * Math.PI) / 180;
  const a = strokeM / 2;
  const pistonTravelM = a * (1 - Math.cos(thetaRad)) + rodM - Math.sqrt((rodM * rodM) - (a * a * Math.sin(thetaRad) * Math.sin(thetaRad)));
  const volumeM3 = clearanceVolM3 + areaM2 * Math.max(0, pistonTravelM);
  return volumeM3 * 1000;
}

export default function WavePage() {
  const [s, setS] = useState(DEFAULTS);
  const [tab, setTab] = useState("runners");
  const [result, setResult] = useState(null);
  const [running, setRunning] = useState(false);
  const [error, setError] = useState(null);

  const set = (key) => (v) => setS((prev) => ({ ...prev, [key]: v }));

  useEffect(() => {
    fetchPresetS54()
      .then((preset) => setS((prev) => ({ ...prev, ...preset })))
      .catch(() => { /* keep local defaults if API not up yet */ });
  }, []);

  const camValid =
    s.i_dur040 < s.i_dur006 && s.i_max_lift > 1 &&
    s.e_dur040 < s.e_dur006 && s.e_max_lift > 1;

  const pvData = useMemo(() => {
    if (!result?.hist) return { workingLoop: [], pumpingLoop: [] };

    const points = result.hist.map((point) => ({
      ...point,
      volume_l: computeCylinderVolumeL(s.bore, s.stroke, s.rod, s.CR, point.theta),
    }));

    const byTheta = [...points].sort((a, b) => a.theta - b.theta);
    const workingLoop = byTheta.filter((p) => p.theta >= 180 && p.theta <= 540);
    const pumpingLoop = [...byTheta.filter((p) => p.theta >= 540), ...byTheta.filter((p) => p.theta <= 180)];

    const pressureOffset = 0.12;
    const closeTheLoop = (series) => {
      if (!series.length) return series;
      const first = series[0];
      return series.concat({
        ...first,
        theta: series[series.length - 1].theta + 0.0001,
        volume_l: first.volume_l,
        P: first.P,
      });
    };
    const offsetLoop = (series, offset) => series.map((point) => ({
      ...point,
      P: point.P + offset,
    }));

    const smoothPressure = (series) => {
      if (series.length < 3) return series;
      const smoothed = series.map((point, index) => {
        const prev = series[Math.max(0, index - 1)]?.P ?? point.P;
        const next = series[Math.min(series.length - 1, index + 1)]?.P ?? point.P;
        return {
          ...point,
          P: (point.P * 0.6) + (prev * 0.2) + (next * 0.2),
        };
      });
      return smoothed;
    };

    const working = smoothPressure(workingLoop).map((point) => ({
      ...point,
      volume_l: Number(point.volume_l.toFixed(3)),
    }));
    const pumping = offsetLoop(smoothPressure(pumpingLoop), pressureOffset).map((point) => ({
      ...point,
      volume_l: Number(point.volume_l.toFixed(3)),
    }));

    const pvPressureMin = Math.max(0.05, Math.min(...working.map((point) => point.P)));
    const pvPressureMax = Math.max(...working.map((point) => point.P));

    return {
      workingLoop: working,
      pumpingLoop: pumping,
      tdcVolume: Math.min(...working.map((point) => point.volume_l)),
      bdcVolume: Math.max(...working.map((point) => point.volume_l)),
      pvPressureMin,
      pvPressureMax,
    };
  }, [result, s.bore, s.stroke, s.rod, s.CR]);

  const eventMarkers = useMemo(() => {
    const allPoints = [...pvData.workingLoop, ...pvData.pumpingLoop];
    if (!allPoints.length) return [];

    const normalize = (value) => ((value % 360) + 360) % 360;
    const findPoint = (theta, preferredLoop = "working") => {
      const target = normalize(theta);
      const candidates = preferredLoop === "working"
        ? pvData.workingLoop
        : pvData.pumpingLoop;
      let best = candidates[0] || allPoints[0];
      let bestDiff = 360;
      for (const point of candidates) {
        const diff = Math.abs(normalize(point.theta) - target);
        if (diff < bestDiff) {
          bestDiff = diff;
          best = point;
        }
      }
      return best;
    };

    const peakPoint = allPoints.reduce((best, point) => (point.P > best.P ? point : best), allPoints[0]);
    const ignitionPoint = findPoint(s.soc, "working");

    const makeLabel = (key, point) => `${key} ${point.P.toFixed(1)} bar @ ${Math.round(point.theta)}°`;

    const ivoPoint = findPoint(s.i_centerline - s.i_dur006 / 2, "pumping");
    const ivcPoint = findPoint(s.i_centerline + s.i_dur006 / 2, "pumping");
    const evoPoint = findPoint(s.e_centerline + s.e_dur006 / 2, "pumping");
    const evcPoint = findPoint(s.e_centerline - s.e_dur006 / 2, "pumping");

    return [
      {
        key: "IVO",
        label: makeLabel("IVO", ivoPoint),
        point: ivoPoint,
        color: T.cyan,
      },
      {
        key: "IVC",
        label: makeLabel("IVC", ivcPoint),
        point: ivcPoint,
        color: T.cyan,
      },
      {
        key: "EVO",
        label: makeLabel("EVO", evoPoint),
        point: evoPoint,
        color: T.red,
      },
      {
        key: "EVC",
        label: makeLabel("EVC", evcPoint),
        point: evcPoint,
        color: T.red,
      },
      {
        key: "IGN",
        label: makeLabel("IGN", ignitionPoint),
        point: ignitionPoint,
        color: T.amber,
      },
      {
        key: "Pmax",
        label: makeLabel("Pmax", peakPoint),
        point: peakPoint,
        color: T.green,
      },
    ];
  }, [pvData, s.i_centerline, s.i_dur006, s.e_centerline, s.e_dur006, s.soc]);


  const runSim = useCallback(async () => {
    if (!camValid) return;
    setRunning(true);
    setError(null);
    try {
      const res = await simulateWave(s);
      setResult(res);
    } catch (err) {
      setError(err.message || String(err));
    } finally {
      setRunning(false);
    }
  }, [s, camValid]);

  return (
    <div>
      <div className="flex items-baseline justify-between mb-6 pb-4 border-b border-hair">
        <div>
          <div className="text-xs md:text-sm tracking-[0.2em] uppercase mb-1 text-cyan">
            Full four-stroke cycle · 1-D wave dynamics · Blair-inspired FV
          </div>
          <h1 className="font-display text-2xl tracking-wide uppercase">720° Closed-Loop Wave Solver</h1>
          <div className="text-xs md:text-sm font-mono text-lo mt-1">BMW S54B32 defaults · Cd 0.75 · 85% throat</div>
        </div>
        <div className="text-xs md:text-sm font-mono text-right text-lo hidden md:block">
          compression · Wiebe combustion · expansion · gas exchange<br />
          Solved on the Python API — P@EVO, peak P, IMEP, power emerge
        </div>
      </div>

      <div className="grid grid-cols-12 gap-6">
        <div className="col-span-12 md:col-span-4">
          <div className="rounded-lg p-5 bg-panel border border-hair">
            <div className="text-xs uppercase tracking-wider mb-4 text-lo">Engine</div>
            <Field label="Bore" unit="mm" value={s.bore} onChange={set("bore")} min={60} max={120} step={0.5} />
            <Field label="Stroke" unit="mm" value={s.stroke} onChange={set("stroke")} min={50} max={120} step={0.5} />
            <Field label="Rod length" unit="mm" value={s.rod} onChange={set("rod")} min={90} max={200} step={0.5} />
            <Field label="Engine speed" unit=" rpm" value={s.rpm} onChange={set("rpm")} min={1000} max={14000} step={100} />
            <Field label="Compression ratio" unit=":1" value={s.CR} onChange={set("CR")} min={8} max={14} step={0.1} />
            <div className="flex gap-2 mt-2">
              {["j1349", "j607"].map((k) => (
                <button
                  key={k}
                  type="button"
                  onClick={() => set("intake_std")(k)}
                  className="px-2 py-1 text-xs md:text-sm rounded border border-hair"
                  style={{
                    background: s.intake_std === k ? T.cyan : T.panel2,
                    color: s.intake_std === k ? T.void : T.lo,
                  }}
                >
                  {k.toUpperCase()}
                </button>
              ))}
            </div>

            <div className="mt-5">
              <div className="flex gap-1 flex-wrap">
                <Tab active={tab === "runners"} onClick={() => setTab("runners")}>Runners</Tab>
                <Tab active={tab === "comb"} onClick={() => setTab("comb")}>Combustion</Tab>
                <Tab active={tab === "fric"} onClick={() => setTab("fric")}>Friction</Tab>
                <Tab active={tab === "intake"} onClick={() => setTab("intake")}>Intake</Tab>
                <Tab active={tab === "exhaust"} onClick={() => setTab("exhaust")}>Exhaust</Tab>
              </div>
              <div className="rounded-b-lg rounded-tr-lg p-5 bg-panel border border-hair border-t-0">
                {tab === "runners" && (
                  <>
                    <Field label="Intake runner length" unit="m" value={s.l_int} onChange={set("l_int")} min={0.1} max={0.8} step={0.01} />
                    <Field label="Intake runner dia" unit="mm" value={s.d_int} onChange={set("d_int")} min={25} max={60} step={1} />
                    <Field label="Exhaust runner length" unit="m" value={s.l_exh} onChange={set("l_exh")} min={0.15} max={1.2} step={0.01} />
                    <Field label="Exhaust runner dia" unit="mm" value={s.d_exh} onChange={set("d_exh")} min={25} max={60} step={1} />
                  </>
                )}
                {tab === "comb" && (
                  <>
                    <Field label="AFR" unit=":1" value={s.afr} onChange={set("afr")} min={11} max={16} step={0.1} />
                    <Field label="Combustion efficiency" unit="" value={s.eta_comb} onChange={set("eta_comb")} min={0.85} max={1} step={0.01} />
                    <Field label="Wall heat loss frac" unit="" value={s.wall_heat_frac} onChange={set("wall_heat_frac")} min={0} max={0.35} step={0.01} />
                    <Field label="Start of combustion" unit="° (360=TDC)" value={s.soc} onChange={set("soc")} min={320} max={365} step={1} />
                    <Field label="Burn duration" unit="°" value={s.burn_dur} onChange={set("burn_dur")} min={30} max={90} step={1} />
                    <Field label="Convergence cycles" unit="" value={s.n_cycles} onChange={set("n_cycles")} min={3} max={8} step={1} />
                  </>
                )}
                {tab === "fric" && (
                  <>
                    <Field label="Chen-Flynn A (const)" unit=" bar" value={s.fmep_a} onChange={set("fmep_a")} min={0.1} max={0.8} step={0.01} />
                    <Field label="B (×Ppeak)" unit="" value={s.fmep_b} onChange={set("fmep_b")} min={0.002} max={0.01} step={0.0005} />
                    <Field label="C (×Sp)" unit="" value={s.fmep_c} onChange={set("fmep_c")} min={0.04} max={0.16} step={0.005} />
                    <Field label="D (×Sp²)" unit="" value={s.fmep_d} onChange={set("fmep_d")} min={0.0003} max={0.002} step={0.0001} />
                  </>
                )}
                {tab === "intake" && (
                  <>
                    <Field label="Valve diameter" unit="mm" value={s.i_valve_dia} onChange={set("i_valve_dia")} min={20} max={55} step={0.5} />
                    <Field label="Throat %" unit="%" value={s.i_throat_pct} onChange={set("i_throat_pct")} min={70} max={100} step={1} />
                    <Field label="Stem diameter" unit="mm" value={s.i_stem_dia} onChange={set("i_stem_dia")} min={4} max={10} step={0.5} />
                    <Field label="Seat angle" unit="°" value={s.i_seat} onChange={set("i_seat")} min={30} max={50} step={1} />
                    <Field label="Cd" unit="" value={s.i_cd} onChange={set("i_cd")} min={0.3} max={1} step={0.01} />
                    <Field label="Dur @0.15mm" unit="°" value={s.i_dur006} onChange={set("i_dur006")} min={180} max={320} step={1} />
                    <Field label="Dur @1mm" unit="°" value={s.i_dur040} onChange={set("i_dur040")} min={140} max={300} step={1} />
                    <Field label="Max lift" unit="mm" value={s.i_max_lift} onChange={set("i_max_lift")} min={2} max={16} step={0.1} />
                    <Field label="Centerline" unit="°" value={s.i_centerline} onChange={set("i_centerline")} min={90} max={150} step={1} />
                  </>
                )}
                {tab === "exhaust" && (
                  <>
                    <Field label="Valve diameter" unit="mm" value={s.e_valve_dia} onChange={set("e_valve_dia")} min={18} max={50} step={0.5} />
                    <Field label="Throat %" unit="%" value={s.e_throat_pct} onChange={set("e_throat_pct")} min={70} max={100} step={1} />
                    <Field label="Stem diameter" unit="mm" value={s.e_stem_dia} onChange={set("e_stem_dia")} min={4} max={10} step={0.5} />
                    <Field label="Seat angle" unit="°" value={s.e_seat} onChange={set("e_seat")} min={30} max={50} step={1} />
                    <Field label="Cd" unit="" value={s.e_cd} onChange={set("e_cd")} min={0.3} max={1} step={0.01} />
                    <Field label="Dur @0.15mm" unit="°" value={s.e_dur006} onChange={set("e_dur006")} min={180} max={320} step={1} />
                    <Field label="Dur @1mm" unit="°" value={s.e_dur040} onChange={set("e_dur040")} min={140} max={300} step={1} />
                    <Field label="Max lift" unit="mm" value={s.e_max_lift} onChange={set("e_max_lift")} min={2} max={16} step={0.1} />
                    <Field label="Centerline (signed)" unit="°" value={s.e_centerline} onChange={set("e_centerline")} min={-150} max={-70} step={1} />
                  </>
                )}
              </div>
            </div>

            <button
              type="button"
              onClick={runSim}
              disabled={running || !camValid}
              className="mt-5 w-full py-3 rounded-lg text-sm uppercase tracking-widest font-semibold border border-hair"
              style={{
                background: running ? T.panel2 : T.cyan,
                color: running ? T.lo : T.void,
              }}
            >
              {running ? "Solving cycles…" : "Run full-cycle simulation"}
            </button>
            {error && (
              <div className="mt-3 text-sm md:text-base font-mono text-danger border border-danger/40 rounded p-2">
                {error}
              </div>
            )}
          </div>
        </div>

        <div className="col-span-12 md:col-span-8">
          {!result && !error && (
            <div className="rounded-lg p-10 text-center bg-panel border border-hair text-lo">
              <div className="text-sm">Configure the engine, then run the simulation.</div>
              <div className="text-xs md:text-sm mt-2 font-mono">
                Solves the full 720° cycle over multiple iterations until periodic —
                everything (P@EVO, peak pressure, IMEP, power) emerges from first principles.
              </div>
            </div>
          )}

          {result && (
            <>
              <div className="grid grid-cols-2 md:grid-cols-4 gap-3 mb-5">
                <Stat
                  label={`Brake power ×${result.n_cylinders}`}
                  value={`${result.brake_hp_n.toFixed(0)}hp`}
                  sub={`indicated ${(result.n_cylinders * result.hp_per_cyl).toFixed(0)}hp − friction`}
                  color={T.cyan}
                />
                <Stat
                  label={`Brake torque ×${result.n_cylinders}`}
                  value={`${result.brake_torque_n.toFixed(0)}Nm`}
                  sub={`BMEP ${result.bmep_bar.toFixed(2)}bar`}
                  color={T.green}
                />
                <Stat
                  label="FMEP (Chen-Flynn)"
                  value={`${result.fmep_bar.toFixed(2)}bar`}
                  sub={`mech eff ${result.mech_eff.toFixed(1)}% · Sp ${result.Sp_mean.toFixed(1)}m/s`}
                  color={T.amber}
                />
                <Stat
                  label="Peak P / P@EVO"
                  value={`${result.Ppeak_bar.toFixed(1)}bar`}
                  sub={`EVO ${result.PatEVO.toFixed(2)}bar · air ${result.air_mg.toFixed(0)}mg · VE ${result.ve_pct?.toFixed(0) ?? "—"}%`}
                />
              </div>

              <div className="rounded-lg p-5 mb-5 bg-panel border border-hair">
                <div className="text-xs uppercase tracking-wider mb-3 text-lo">
                  Cylinder pressure, full cycle (log scale) — 360° = firing TDC
                </div>
                <ResponsiveContainer width="100%" height={240}>
                  <LineChart data={result.hist} margin={{ top: 5, right: 15, bottom: 5, left: 0 }}>
                    <CartesianGrid stroke={T.line} strokeDasharray="2 4" />
                    <XAxis
                      dataKey="theta"
                      stroke={T.lo}
                      tick={{ fontSize: 13, fill: T.lo }}
                      ticks={[0, 360, 720]}
                      label={{ value: "Crank angle (°) — 360° = firing TDC", position: "insideBottom", offset: -5, fill: T.lo }}
                    />
                    <YAxis
                      dataKey="logP"
                      stroke={T.lo}
                      tick={{ fontSize: 13, fill: T.lo }}
                      tickFormatter={(v) => `${Math.pow(10, v).toFixed(1)}`}
                    />
                    <Tooltip
                      contentStyle={{ background: T.panel2, border: `1px solid ${T.hair}`, fontSize: 14 }}
                      labelStyle={{ color: T.lo }}
                      formatter={(v, n, entry) => {
                        if (n === "logP") {
                          const pressure = Math.pow(10, v);
                          return [`${pressure.toFixed(2)} bar`, "Cylinder P"];
                        }
                        return [v, n];
                      }}
                      labelFormatter={(value) => `θ = ${value}°`}
                    />
                    <ReferenceLine x={360} stroke={T.amber} strokeDasharray="4 4" label={{ value: "TDC fire", fill: T.amber, fontSize: 10 }} />
                    {eventMarkers.map((event) => (
                      <Scatter
                        key={event.key}
                        data={[{
                          theta: event.point.theta,
                          logP: Math.log10(event.point.P),
                          label: event.label,
                        }]}
                        fill={event.color}
                        shape="circle"
                        line={false}
                        r={6}
                        label={{
                          position: "top",
                          value: event.label,
                          fill: event.color,
                          fontSize: 12,
                          fontWeight: 700,
                        }}
                      />
                    ))}
                    <Line type="linear" dataKey="logP" stroke={T.green} strokeWidth={2} dot={false} name="Cylinder P" />
                  </LineChart>
                </ResponsiveContainer>
              </div>

              <div className="rounded-lg p-5 mb-5 bg-panel border border-hair">
                <div className="text-xs uppercase tracking-wider mb-3 text-lo">Pressure–volume diagram (cylinder)</div>
                <ResponsiveContainer width="100%" height={280}>
                  <LineChart data={pvData.workingLoop} margin={{ top: 20, right: 15, bottom: 15, left: 0 }}>
                    <CartesianGrid stroke={T.line} strokeDasharray="2 4" />
                    <XAxis
                      dataKey="volume_l"
                      stroke={T.lo}
                      tick={{ fontSize: 11, fill: T.lo }}
                      domain={["dataMin", "dataMax"]}
                      type="number"
                      ticks={[pvData.tdcVolume, pvData.bdcVolume]}
                      tickFormatter={(value) => {
                        if (Math.abs(value - (pvData.tdcVolume ?? 0)) < 0.01) return "TDC";
                        if (Math.abs(value - (pvData.bdcVolume ?? 0)) < 0.01) return "BDC";
                        return "";
                      }}
                      label={{ value: "TDC ← volume → BDC", position: "insideBottom", offset: -5, fill: T.lo }}
                    />
                    <YAxis dataKey="P" stroke={T.lo} tick={{ fontSize: 13, fill: T.lo }} scale="log" domain={[pvData.pvPressureMin, pvData.pvPressureMax]} />
                    <Tooltip
                      contentStyle={{ background: T.panel2, border: `1px solid ${T.hair}`, fontSize: 14 }}
                      labelStyle={{ color: T.lo }}
                      formatter={(value, name, entry) => {
                        const payload = entry?.payload;
                        if (!payload) return [`${Number(value).toFixed(2)} bar`, "Pressure"];
                        return [
                          `${Number(payload.P).toFixed(2)} bar · V = ${Number(payload.volume_l).toFixed(2)} L`,
                          "PV point",
                        ];
                      }}
                      labelFormatter={(value) => `Volume = ${Number(value).toFixed(2)} L`}
                    />
                    <ReferenceLine y={1} stroke={T.hair} strokeDasharray="2 2" />
                    {eventMarkers.map((event) => (
                      <Scatter
                        key={event.key}
                        data={[{
                          volume_l: event.point.volume_l,
                          P: event.point.P,
                          label: event.label,
                        }]}
                        fill={event.color}
                        shape="circle"
                        line={false}
                        r={6}
                        label={{
                          position: "top",
                          value: event.label,
                          fill: event.color,
                          fontSize: 12,
                          fontWeight: 700,
                        }}
                      />
                    ))}
                    <Line type="linear" dataKey="P" stroke={T.cyan} strokeWidth={2.2} dot={false} name="Working loop" />
                    <Line type="linear" dataKey="P" data={pvData.pumpingLoop} stroke={T.red} strokeWidth={1.8} dot={false} name="Pumping loop" />
                  </LineChart>
                </ResponsiveContainer>
                <div className="mt-3 text-xs md:text-sm font-mono text-lo leading-relaxed">
                  Working loop: compression → ignition → pressure rise → peak pressure after TDC → blowdown → expansion. Pumping loop: intake and exhaust gas exchange with IVO/IVC/EVO/EVC markers.
                </div>
              </div>

              <div className="rounded-lg p-5 mb-5 bg-panel border border-hair">
                <div className="text-xs uppercase tracking-wider mb-3 text-lo">Wave-affected port pressures at the valves</div>
                <ResponsiveContainer width="100%" height={200}>
                  <LineChart data={result.hist} margin={{ top: 5, right: 15, bottom: 5, left: 0 }}>
                    <CartesianGrid stroke={T.line} strokeDasharray="2 4" />
                    <XAxis dataKey="theta" stroke={T.lo} tick={{ fontSize: 13, fill: T.lo }} />
                    <YAxis stroke={T.lo} tick={{ fontSize: 13, fill: T.lo }} />
                    <Tooltip contentStyle={{ background: T.panel2, border: `1px solid ${T.hair}`, fontSize: 14 }} labelStyle={{ color: T.lo }} />
                    <Legend wrapperStyle={{ fontSize: 13 }} />
                    <Line type="monotone" dataKey="Pi" stroke={T.cyan} strokeWidth={1.6} dot={false} name="Intake port @ valve" />
                    <Line type="monotone" dataKey="Pe" stroke={T.red} strokeWidth={1.6} dot={false} name="Exhaust port @ valve" />
                  </LineChart>
                </ResponsiveContainer>
              </div>

              <div className="rounded-lg p-5 bg-panel border border-hair">
                <div className="text-xs uppercase tracking-wider mb-3 text-lo">True local Mach at each valve (wave-coupled)</div>
                <ResponsiveContainer width="100%" height={200}>
                  <LineChart data={result.hist} margin={{ top: 5, right: 15, bottom: 5, left: 0 }}>
                    <CartesianGrid stroke={T.line} strokeDasharray="2 4" />
                    <XAxis dataKey="theta" stroke={T.lo} tick={{ fontSize: 13, fill: T.lo }} />
                    <YAxis stroke={T.lo} tick={{ fontSize: 13, fill: T.lo }} domain={[-1.05, 1.05]} />
                    <Tooltip contentStyle={{ background: T.panel2, border: `1px solid ${T.hair}`, fontSize: 14 }} labelStyle={{ color: T.lo }} />
                    <Legend wrapperStyle={{ fontSize: 13 }} />
                    <ReferenceLine y={0} stroke={T.hair} />
                    <ReferenceLine y={1} stroke={T.red} strokeDasharray="4 4" />
                    <ReferenceLine y={-1} stroke={T.red} strokeDasharray="4 4" />
                    <Line type="monotone" dataKey="eM" stroke={T.red} strokeWidth={2} dot={false} name="Exhaust Mach" />
                    <Line type="monotone" dataKey="iM" stroke={T.cyan} strokeWidth={2} dot={false} name="Intake Mach" />
                  </LineChart>
                </ResponsiveContainer>
              </div>

              <div className="mt-4 text-xs md:text-sm font-mono leading-relaxed text-lo">
                Complete 720° closed-loop cycle: 1-D finite-volume Euler runners (HLL, CFL-stepped)
                run continuously through all four strokes. Cylinder: 0-D mass/energy with piston work
                and single-zone Wiebe heat release. Simplifications: γ=1.38, adiabatic walls,
                no pipe friction, straight runners, infinite plenum/collector.
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  );
}
