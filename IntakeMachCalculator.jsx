import React, { useMemo, useState, useCallback } from "react";
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip,
  ReferenceLine, ResponsiveContainer, Legend
} from "recharts";

// ---- design tokens -------------------------------------------------------
const T = {
  void: "#0A0D10", panel: "#12161B", panel2: "#171C22", line: "#232A32", hair: "#2E363F",
  cyan: "#4FD8C4", amber: "#F2A93B", red: "#E5644B", green: "#7FD858", hi: "#E7ECEF", lo: "#7C8791",
};

// ---- physics constants ------------------------------------------------------
const R_GAS = 287.05, GAMMA = 1.38;
const CV = R_GAS / (GAMMA - 1), CP = (GAMMA * R_GAS) / (GAMMA - 1);
const LHV = 44.0e6;
const CRIT = Math.pow((GAMMA + 1) / 2, GAMMA / (GAMMA - 1));
const LIFT_THRESHOLD_MM = 0.15;
const SAE_STANDARDS = {
  j1349: { P: 99000.0, T: 298.15, label: "SAE J1349 (99 kPa, 25°C)" },
  j607: { P: 101325.0, T: 288.71, label: "SAE J607 (29.92 inHg, 60°F)" },
};

function pistonS(thRad, r, L) {
  const s = Math.sin(thRad), c = Math.cos(thRad);
  return r * (1 - c) + L - Math.sqrt(Math.max(L * L - r * r * s * s, 1e-12));
}
function solveCam(D006, D040, Lmax) {
  const x1 = (1 - D040 / D006) / 2;
  return { Dt: D006, p: Math.log(1.0 / Lmax) / Math.log(Math.sin(Math.PI * x1)) };
}
function liftAt(th, cl, Dt, p, Lmax) {
  let x = ((th - (cl - Dt / 2)) % 720 + 720) % 720;
  x = x / Dt;
  if (x < 0 || x > 1) return 0;
  return Lmax * Math.pow(Math.max(Math.sin(Math.PI * x), 0), p);
}
function curtain(l, Dv, seat, n) {
  return n * Math.PI * (Dv / 1000) * (l / 1000) * Math.cos((seat * Math.PI) / 180);
}
function throat(Dt, Ds, n) {
  return n * (Math.PI / 4) * (Math.pow(Dt / 1000, 2) - Math.pow(Ds / 1000, 2));
}
function orifice(Pu, Tu, Pd, A) {
  if (A <= 0 || Pu <= Pd) return [0, 0];
  const ratio = Pu / Pd;
  const M = ratio >= CRIT ? 1
    : Math.sqrt(Math.max(0, (2 / (GAMMA - 1)) * (Math.pow(ratio, (GAMMA - 1) / GAMMA) - 1)));
  const t = 1 + ((GAMMA - 1) / 2) * M * M;
  return [A * Pu * Math.sqrt(GAMMA / (R_GAS * Tu)) * M * Math.pow(t, -((GAMMA + 1) / (2 * (GAMMA - 1)))), M];
}
function wiebe(th, soc, dur) {
  const x = (th - soc) / dur;
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  return 1 - Math.exp(-5.0 * x * x);
}
function primOf(rho, mom, E) {
  const r = Math.max(rho, 1e-6), u = mom / r;
  const p = Math.max((GAMMA - 1) * (E - 0.5 * r * u * u), 1e2);
  return [r, u, p];
}
function hllInto(out, rL0, mL, EL, rR0, mR, ER) {
  const [rL, uL, pL] = primOf(rL0, mL, EL);
  const [rR, uR, pR] = primOf(rR0, mR, ER);
  const aL = Math.sqrt(GAMMA * pL / rL), aR = Math.sqrt(GAMMA * pR / rR);
  const SL = Math.min(uL - aL, uR - aR), SR = Math.max(uL + aL, uR + aR);
  const FL0 = rL * uL, FL1 = rL * uL * uL + pL, FL2 = uL * (EL + pL);
  const FR0 = rR * uR, FR1 = rR * uR * uR + pR, FR2 = uR * (ER + pR);
  if (SL >= 0) { out[0] = FL0; out[1] = FL1; out[2] = FL2; return; }
  if (SR <= 0) { out[0] = FR0; out[1] = FR1; out[2] = FR2; return; }
  const inv = 1 / (SR - SL);
  out[0] = (SR * FL0 - SL * FR0 + SL * SR * (rR0 - rL0)) * inv;
  out[1] = (SR * FL1 - SL * FR1 + SL * SR * (mR - mL)) * inv;
  out[2] = (SR * FL2 - SL * FR2 + SL * SR * (ER - EL)) * inv;
}

function runFullCycle(s) {
  const bore = s.bore / 1000, stroke = s.stroke / 1000, rod = s.rod / 1000;
  const r = stroke / 2, Ap = (Math.PI / 4) * bore * bore;
  const Vc = (Ap * stroke) / (s.CR - 1), Vd = Ap * stroke;
  const omega = (2 * Math.PI * s.rpm) / 60;

  const std = SAE_STANDARDS[s.intakeStd];
  const Pamb = std.P, Tamb = std.T;
  const Pcoll = 1.05e5, Tcoll = 1173;

  const N = 36;
  const dxI = s.lInt / N, dxE = s.lExh / N;
  const Aint = (Math.PI / 4) * Math.pow(s.dInt / 1000, 2);
  const Aexh = (Math.PI / 4) * Math.pow(s.dExh / 1000, 2);

  const iCam = solveCam(s.iDur006, s.iDur040, s.iMaxLift);
  const eCam = solveCam(s.eDur006, s.eDur040, s.eMaxLift);
  const AiT = throat(s.iValveDia * s.iThroatPct / 100, s.iStemDia, s.iNValves);
  const AeT = throat(s.eValveDia * s.eThroatPct / 100, s.eStemDia, s.eNValves);

  const rhoI = new Float64Array(N), momI = new Float64Array(N), EI = new Float64Array(N);
  const rhoE = new Float64Array(N), momE = new Float64Array(N), EE = new Float64Array(N);
  const rho0i = Pamb / (R_GAS * Tamb), rho0e = Pcoll / (R_GAS * Tcoll);
  for (let k = 0; k < N; k++) {
    rhoI[k] = rho0i; momI[k] = 0; EI[k] = Pamb / (GAMMA - 1);
    rhoE[k] = rho0e; momE[k] = 0; EE[k] = Pcoll / (GAMMA - 1);
  }

  let V = Vc + Ap * pistonS(0, r, rod);
  let Tcyl = 800, mCyl = (1.2e5 * V) / (R_GAS * Tcyl);

  const f = [0, 0, 0];
  const FI = new Float64Array(3 * (N + 1)), FE = new Float64Array(3 * (N + 1));
  const cycles = [];
  let lastHist = null;

  for (let cyc = 0; cyc < s.nCycles; cyc++) {
    let theta = 0, mAir = 0, work = 0, Qt = null, bprev = 0;
    const hist = [];
    let Ppeak = 0, PatEVO = 0;
    let step = 0;
    let guard = 0;

    while (theta < 720 && guard < 120000) {
      guard++;
      let amax = 1;
      for (let k = 0; k < N; k++) {
        let [rr, u, p] = primOf(rhoI[k], momI[k], EI[k]);
        let a = Math.abs(u) + Math.sqrt(GAMMA * p / rr);
        if (a > amax) amax = a;
        [rr, u, p] = primOf(rhoE[k], momE[k], EE[k]);
        a = Math.abs(u) + Math.sqrt(GAMMA * p / rr);
        if (a > amax) amax = a;
      }
      let dt = 0.45 * Math.min(dxI, dxE) / amax;
      let dth = (omega * dt * 180) / Math.PI;
      if (theta + dth > 720) { dth = 720 - theta; dt = (dth * Math.PI) / 180 / omega; }

      const Pcyl = (mCyl * R_GAS * Tcyl) / V;
      if (Pcyl > Ppeak) Ppeak = Pcyl;
      if (Math.abs(theta - 470) < 1.0 && PatEVO === 0) PatEVO = Pcyl;

      const li = liftAt(theta, s.iCenterline, iCam.Dt, iCam.p, s.iMaxLift);
      const le = liftAt(theta, s.eCenterline, eCam.Dt, eCam.p, s.eMaxLift);
      const AeffI = s.iCd * (li > LIFT_THRESHOLD_MM ? Math.min(curtain(li, s.iValveDia, s.iSeat, s.iNValves), AiT) : 0);
      const AeffE = s.eCd * (le > LIFT_THRESHOLD_MM ? Math.min(curtain(le, s.eValveDia, s.eSeat, s.eNValves), AeT) : 0);

      const [rIv, , pIv] = primOf(rhoI[N - 1], momI[N - 1], EI[N - 1]);
      const TIv = pIv / (rIv * R_GAS);
      const [rEv, , pEv] = primOf(rhoE[0], momE[0], EE[0]);
      const TEv = pEv / (rEv * R_GAS);

      const [miIn, Mif] = orifice(pIv, TIv, Pcyl, AeffI);
      const [miOut, Mir] = orifice(Pcyl, Tcyl, pIv, AeffI);
      const [meOut, Mef] = orifice(Pcyl, Tcyl, pEv, AeffE);
      const [meIn, Mer] = orifice(pEv, TEv, Pcyl, AeffE);
      const mdotI = miIn - miOut, mdotE = meOut - meIn;
      const MiS = mdotI >= 0 ? Mif : -Mir;
      const MeS = mdotE >= 0 ? Mef : -Mer;

      // intake pipe (valve at right end)
      hllInto(f, rho0i, 0, Pamb / (GAMMA - 1), rhoI[0], momI[0], EI[0]);
      FI[0] = f[0]; FI[1] = f[1]; FI[2] = f[2];
      for (let k = 1; k < N; k++) {
        hllInto(f, rhoI[k - 1], momI[k - 1], EI[k - 1], rhoI[k], momI[k], EI[k]);
        FI[3 * k] = f[0]; FI[3 * k + 1] = f[1]; FI[3 * k + 2] = f[2];
      }
      hllInto(f, rhoI[N - 1], momI[N - 1], EI[N - 1], rhoI[N - 1], -momI[N - 1], EI[N - 1]);
      FI[3 * N] = f[0]; FI[3 * N + 1] = f[1]; FI[3 * N + 2] = f[2];
      for (let k = 0; k < N; k++) {
        rhoI[k] -= (dt / dxI) * (FI[3 * (k + 1)] - FI[3 * k]);
        momI[k] -= (dt / dxI) * (FI[3 * (k + 1) + 1] - FI[3 * k + 1]);
        EI[k] -= (dt / dxI) * (FI[3 * (k + 1) + 2] - FI[3 * k + 2]);
      }
      // exhaust pipe (valve at left end)
      hllInto(f, rhoE[0], -momE[0], EE[0], rhoE[0], momE[0], EE[0]);
      FE[0] = f[0]; FE[1] = f[1]; FE[2] = f[2];
      for (let k = 1; k < N; k++) {
        hllInto(f, rhoE[k - 1], momE[k - 1], EE[k - 1], rhoE[k], momE[k], EE[k]);
        FE[3 * k] = f[0]; FE[3 * k + 1] = f[1]; FE[3 * k + 2] = f[2];
      }
      hllInto(f, rhoE[N - 1], momE[N - 1], EE[N - 1], rho0e, 0, Pcoll / (GAMMA - 1));
      FE[3 * N] = f[0]; FE[3 * N + 1] = f[1]; FE[3 * N + 2] = f[2];
      for (let k = 0; k < N; k++) {
        rhoE[k] -= (dt / dxE) * (FE[3 * (k + 1)] - FE[3 * k]);
        momE[k] -= (dt / dxE) * (FE[3 * (k + 1) + 1] - FE[3 * k + 1]);
        EE[k] -= (dt / dxE) * (FE[3 * (k + 1) + 2] - FE[3 * k + 2]);
      }

      // valve exchanges
      const cVi = Aint * dxI, cVe = Aexh * dxE;
      const dm = mdotI * dt;
      if (dm > 0) {
        const h = CP * TIv;
        rhoI[N - 1] -= dm / cVi; EI[N - 1] -= (dm * h) / cVi;
        const mn = mCyl + dm;
        Tcyl = (mCyl * CV * Tcyl + dm * h) / (mn * CV);
        mCyl = mn; mAir += dm;
      } else if (dm < 0) {
        const h = CP * Tcyl;
        rhoI[N - 1] += -dm / cVi; EI[N - 1] += (-dm * h) / cVi;
        mCyl = Math.max(mCyl + dm, 1e-9); mAir = Math.max(mAir + dm, 0);
      }
      const dme = mdotE * dt;
      if (dme > 0) {
        const h = CP * Tcyl;
        rhoE[0] += dme / cVe; EE[0] += (dme * h) / cVe;
        mCyl = Math.max(mCyl - dme, 1e-9);
      } else if (dme < 0) {
        const h = CP * TEv;
        rhoE[0] -= -dme / cVe; EE[0] -= (-dme * h) / cVe;
        const mn = mCyl - dme;
        Tcyl = (mCyl * CV * Tcyl + (-dme) * h) / (mn * CV);
        mCyl = mn;
      }

      // combustion
      if (Qt === null && theta >= s.soc) Qt = s.etaComb * (mAir / s.afr) * LHV;
      if (Qt !== null) {
        const b = wiebe(theta, s.soc, s.burnDur);
        Tcyl += (Qt * (b - bprev)) / (mCyl * 1021.0); // burned-gas cv (gamma_burn~1.28)
        bprev = b;
      }

      // piston work
      const Vnext = Vc + Ap * pistonS(((theta + dth) * Math.PI) / 180, r, rod);
      const dV = Vnext - V;
      work += Pcyl * dV;
      Tcyl = (mCyl * CV * Tcyl - Pcyl * dV) / (mCyl * CV);
      V = Vnext; theta += dth;

      if (step % 10 === 0 && cyc === s.nCycles - 1) {
        hist.push({
          theta: Math.round(theta * 10) / 10,
          P: Number((Pcyl / 1e5).toFixed(3)),
          logP: Number(Math.log10(Pcyl / 1e5).toFixed(4)),
          Tc: Math.round(Tcyl),
          Pi: Number((pIv / 1e5).toFixed(3)),
          Pe: Number((pEv / 1e5).toFixed(3)),
          iM: Number(MiS.toFixed(3)),
          eM: Number(MeS.toFixed(3)),
        });
      }
      step++;
    }

    const imep = work / Vd;
    cycles.push({ cycle: cyc + 1, airMg: mAir * 1e6, imepBar: imep / 1e5, PatEVO: PatEVO / 1e5, PpeakBar: Ppeak / 1e5 });
    if (cyc === s.nCycles - 1) {
      const torque = (imep * Vd) / (4 * Math.PI);
      const power = torque * omega;
      // Chen-Flynn friction: FMEP = A + B*Ppeak + C*Sp + D*Sp^2 (bar)
      const Sp = (2 * stroke * s.rpm) / 60;
      const PpeakBar = Ppeak / 1e5;
      const fmepBar = s.fmepA + s.fmepB * PpeakBar + s.fmepC * Sp + s.fmepD * Sp * Sp;
      const bmep = imep - fmepBar * 1e5;
      const torqueB = (bmep * Vd) / (4 * Math.PI);
      const powerB = torqueB * omega;
      lastHist = {
        hist, imepBar: imep / 1e5, torqueNm: torque, hpPerCyl: power / 745.7,
        airMg: mAir * 1e6, PpeakBar, PatEVO: PatEVO / 1e5,
        fmepBar, bmepBar: bmep / 1e5, SpMean: Sp,
        brakeTorque6: 6 * torqueB, brakeHp6: (6 * powerB) / 745.7,
        mechEff: (100 * bmep) / imep,
      };
    }
    Qt = null; bprev = 0;
  }

  return { cycles, ...lastHist };
}

// ---- UI --------------------------------------------------------------------
function Field({ label, unit, value, onChange, min, max, step }) {
  return (
    <div className="mb-4">
      <div className="flex items-baseline justify-between mb-1">
        <label className="text-xs tracking-wide uppercase" style={{ color: T.lo }}>{label}</label>
        <span className="text-sm font-mono" style={{ color: T.cyan }}>{value}{unit}</span>
      </div>
      <input type="range" min={min} max={max} step={step} value={value}
             onChange={(e) => onChange(parseFloat(e.target.value))}
             className="w-full" style={{ accentColor: T.cyan }} />
    </div>
  );
}
function Stat({ label, value, sub, color }) {
  return (
    <div className="px-4 py-3 rounded-md" style={{ background: T.panel2, border: `1px solid ${T.hair}` }}>
      <div className="text-[10px] uppercase tracking-wider mb-1" style={{ color: T.lo }}>{label}</div>
      <div className="text-2xl font-mono font-semibold" style={{ color: color || T.hi }}>{value}</div>
      {sub && <div className="text-[11px] mt-0.5 font-mono" style={{ color: T.lo }}>{sub}</div>}
    </div>
  );
}
function Tab({ active, onClick, children }) {
  return (
    <button onClick={onClick} className="px-3 py-1.5 text-xs uppercase tracking-wider rounded-t-md"
      style={{ color: active ? T.void : T.lo, background: active ? T.cyan : "transparent",
        border: `1px solid ${T.hair}`, borderBottom: active ? "none" : `1px solid ${T.hair}` }}>
      {children}
    </button>
  );
}

export default function FullCycleSimulator() {
  const [bore, setBore] = useState(87);
  const [stroke, setStroke] = useState(91);
  const [rod, setRod] = useState(139);
  const [rpm, setRpm] = useState(8000);
  const [CR, setCR] = useState(11.5);
  const [intakeStd, setIntakeStd] = useState("j1349");

  const [lInt, setLInt] = useState(0.30);
  const [dInt, setDInt] = useState(40);
  const [lExh, setLExh] = useState(0.45);
  const [dExh, setDExh] = useState(38);

  const [afr, setAfr] = useState(12.8);
  const [etaComb, setEtaComb] = useState(0.97);
  const [soc, setSoc] = useState(340);
  const [burnDur, setBurnDur] = useState(55);
  const [nCycles, setNCycles] = useState(5);
  const [fmepA, setFmepA] = useState(0.35);
  const [fmepB, setFmepB] = useState(0.005);
  const [fmepC, setFmepC] = useState(0.09);
  const [fmepD, setFmepD] = useState(0.0009);

  const [iValveDia, setIValveDia] = useState(35);
  const [iThroatPct, setIThroatPct] = useState(85);
  const [iStemDia, setIStemDia] = useState(6);
  const [iNValves] = useState(2);
  const [iSeat, setISeat] = useState(45);
  const [iCd, setICd] = useState(0.85);
  const [iDur006, setIDur006] = useState(260);
  const [iDur040, setIDur040] = useState(218);
  const [iMaxLift, setIMaxLift] = useState(12.0);
  const [iCenterline, setICenterline] = useState(130);

  const [eValveDia, setEValveDia] = useState(30.5);
  const [eThroatPct, setEThroatPct] = useState(83);
  const [eStemDia, setEStemDia] = useState(6);
  const [eNValves] = useState(2);
  const [eSeat, setESeat] = useState(45);
  const [eCd, setECd] = useState(0.85);
  const [eDur006, setEDur006] = useState(260);
  const [eDur040, setEDur040] = useState(218);
  const [eMaxLift, setEMaxLift] = useState(12.0);
  const [eCenterline, setECenterline] = useState(-120);

  const [tab, setTab] = useState("runners");
  const [result, setResult] = useState(null);
  const [running, setRunning] = useState(false);

  const camValid = iDur040 < iDur006 && iMaxLift > 1 && eDur040 < eDur006 && eMaxLift > 1;

  const runSim = useCallback(() => {
    if (!camValid) return;
    setRunning(true);
    // let the UI paint the "running" state before the heavy compute
    setTimeout(() => {
      const res = runFullCycle({
        bore, stroke, rod, rpm, CR, intakeStd, lInt, dInt, lExh, dExh,
        afr, etaComb, soc, burnDur, nCycles, fmepA, fmepB, fmepC, fmepD,
        iValveDia, iThroatPct, iStemDia, iNValves, iSeat, iCd, iDur006, iDur040, iMaxLift, iCenterline,
        eValveDia, eThroatPct, eStemDia, eNValves, eSeat, eCd, eDur006, eDur040, eMaxLift, eCenterline,
      });
      setResult(res);
      setRunning(false);
    }, 30);
  }, [bore, stroke, rod, rpm, CR, intakeStd, lInt, dInt, lExh, dExh, afr, etaComb, soc, burnDur, nCycles, fmepA, fmepB, fmepC, fmepD,
      iValveDia, iThroatPct, iStemDia, iNValves, iSeat, iCd, iDur006, iDur040, iMaxLift, iCenterline,
      eValveDia, eThroatPct, eStemDia, eNValves, eSeat, eCd, eDur006, eDur040, eMaxLift, eCenterline, camValid]);

  return (
    <div className="min-h-screen w-full font-sans" style={{ background: T.void, color: T.hi }}>
      <div className="max-w-6xl mx-auto px-6 py-8">
        <div className="flex items-baseline justify-between mb-6 pb-4" style={{ borderBottom: `1px solid ${T.hair}` }}>
          <div>
            <div className="text-[11px] tracking-[0.2em] uppercase mb-1" style={{ color: T.cyan }}>
              Full four-stroke cycle · 1-D wave dynamics · Blair-inspired
            </div>
            <h1 className="text-xl font-semibold tracking-tight">720° Closed-Loop Engine Simulator</h1>
          </div>
          <div className="text-[11px] font-mono text-right" style={{ color: T.lo }}>
            compression · Wiebe combustion · expansion · gas exchange<br />
            P@EVO, peak P, IMEP, torque &amp; power all EMERGE — nothing prescribed
          </div>
        </div>

        <div className="grid grid-cols-12 gap-6">
          {/* left: inputs */}
          <div className="col-span-12 md:col-span-4">
            <div className="rounded-lg p-5" style={{ background: T.panel, border: `1px solid ${T.hair}` }}>
              <div className="text-xs uppercase tracking-wider mb-4" style={{ color: T.lo }}>Engine</div>
              <Field label="Bore" unit="mm" value={bore} onChange={setBore} min={60} max={120} step={0.5} />
              <Field label="Stroke" unit="mm" value={stroke} onChange={setStroke} min={50} max={120} step={0.5} />
              <Field label="Rod length" unit="mm" value={rod} onChange={setRod} min={90} max={200} step={0.5} />
              <Field label="Engine speed" unit=" rpm" value={rpm} onChange={setRpm} min={1000} max={14000} step={100} />
              <Field label="Compression ratio" unit=":1" value={CR} onChange={setCR} min={8} max={14} step={0.1} />
              <div className="flex gap-2 mt-2">
                {Object.keys(SAE_STANDARDS).map((k) => (
                  <button key={k} onClick={() => setIntakeStd(k)} className="px-2 py-1 text-[11px] rounded"
                    style={{ background: intakeStd === k ? T.cyan : T.panel2, color: intakeStd === k ? T.void : T.lo, border: `1px solid ${T.hair}` }}>
                    {k.toUpperCase()}
                  </button>
                ))}
              </div>
            </div>

            <div className="mt-5">
              <div className="flex gap-1 flex-wrap">
                <Tab active={tab === "runners"} onClick={() => setTab("runners")}>Runners</Tab>
                <Tab active={tab === "comb"} onClick={() => setTab("comb")}>Combustion</Tab>
                <Tab active={tab === "fric"} onClick={() => setTab("fric")}>Friction</Tab>
                <Tab active={tab === "intake"} onClick={() => setTab("intake")}>Intake</Tab>
                <Tab active={tab === "exhaust"} onClick={() => setTab("exhaust")}>Exhaust</Tab>
              </div>
              <div className="rounded-b-lg rounded-tr-lg p-5" style={{ background: T.panel, border: `1px solid ${T.hair}` }}>
                {tab === "runners" && (
                  <>
                    <Field label="Intake runner length" unit="m" value={lInt} onChange={setLInt} min={0.10} max={0.80} step={0.01} />
                    <Field label="Intake runner dia" unit="mm" value={dInt} onChange={setDInt} min={25} max={60} step={1} />
                    <Field label="Exhaust runner length" unit="m" value={lExh} onChange={setLExh} min={0.15} max={1.20} step={0.01} />
                    <Field label="Exhaust runner dia" unit="mm" value={dExh} onChange={setDExh} min={25} max={60} step={1} />
                  </>
                )}
                {tab === "comb" && (
                  <>
                    <Field label="AFR" unit=":1" value={afr} onChange={setAfr} min={11} max={16} step={0.1} />
                    <Field label="Combustion efficiency" unit="" value={etaComb} onChange={setEtaComb} min={0.85} max={1.0} step={0.01} />
                    <Field label="Start of combustion" unit="° (360=TDC)" value={soc} onChange={setSoc} min={320} max={365} step={1} />
                    <Field label="Burn duration" unit="°" value={burnDur} onChange={setBurnDur} min={30} max={90} step={1} />
                    <Field label="Convergence cycles" unit="" value={nCycles} onChange={setNCycles} min={3} max={8} step={1} />
                    <div className="text-[11px] leading-relaxed mt-2" style={{ color: T.lo }}>
                      Single-zone Wiebe heat release (a=5, m=2). Fuel mass = inducted air / AFR.
                      Simulation loops full cycles until periodic — watch the cycle table converge.
                    </div>
                  </>
                )}
                {tab === "fric" && (
                  <>
                    <Field label="Chen-Flynn A (const)" unit=" bar" value={fmepA} onChange={setFmepA} min={0.1} max={0.8} step={0.01} />
                    <Field label="B (×Ppeak)" unit="" value={fmepB} onChange={setFmepB} min={0.002} max={0.01} step={0.0005} />
                    <Field label="C (×Sp)" unit="" value={fmepC} onChange={setFmepC} min={0.04} max={0.16} step={0.005} />
                    <Field label="D (×Sp²)" unit="" value={fmepD} onChange={setFmepD} min={0.0003} max={0.002} step={0.0001} />
                    <div className="text-[11px] leading-relaxed mt-2" style={{ color: T.lo }}>
                      FMEP = A + B·P_peak + C·S_p + D·S_p². Brake = indicated − friction.
                      The Sp and Sp² terms dominate at redline — friction eats
                      disproportionately more at high rpm.
                    </div>
                  </>
                )}
                {tab === "intake" && (
                  <>
                    <Field label="Valve diameter" unit="mm" value={iValveDia} onChange={setIValveDia} min={20} max={55} step={0.5} />
                    <Field label="Throat %" unit="%" value={iThroatPct} onChange={setIThroatPct} min={70} max={100} step={1} />
                    <Field label="Stem diameter" unit="mm" value={iStemDia} onChange={setIStemDia} min={4} max={10} step={0.5} />
                    <Field label="Seat angle" unit="°" value={iSeat} onChange={setISeat} min={30} max={50} step={1} />
                    <Field label="Cd" unit="" value={iCd} onChange={setICd} min={0.3} max={1.0} step={0.01} />
                    <Field label="Dur @0.15mm" unit="°" value={iDur006} onChange={setIDur006} min={180} max={320} step={1} />
                    <Field label="Dur @1mm" unit="°" value={iDur040} onChange={setIDur040} min={140} max={300} step={1} />
                    <Field label="Max lift" unit="mm" value={iMaxLift} onChange={setIMaxLift} min={2} max={16} step={0.1} />
                    <Field label="Centerline" unit="°" value={iCenterline} onChange={setICenterline} min={90} max={150} step={1} />
                  </>
                )}
                {tab === "exhaust" && (
                  <>
                    <Field label="Valve diameter" unit="mm" value={eValveDia} onChange={setEValveDia} min={18} max={50} step={0.5} />
                    <Field label="Throat %" unit="%" value={eThroatPct} onChange={setEThroatPct} min={70} max={100} step={1} />
                    <Field label="Stem diameter" unit="mm" value={eStemDia} onChange={setEStemDia} min={4} max={10} step={0.5} />
                    <Field label="Seat angle" unit="°" value={eSeat} onChange={setESeat} min={30} max={50} step={1} />
                    <Field label="Cd" unit="" value={eCd} onChange={setECd} min={0.3} max={1.0} step={0.01} />
                    <Field label="Dur @0.15mm" unit="°" value={eDur006} onChange={setEDur006} min={180} max={320} step={1} />
                    <Field label="Dur @1mm" unit="°" value={eDur040} onChange={setEDur040} min={140} max={300} step={1} />
                    <Field label="Max lift" unit="mm" value={eMaxLift} onChange={setEMaxLift} min={2} max={16} step={0.1} />
                    <Field label="Centerline (signed)" unit="°" value={eCenterline} onChange={setECenterline} min={-150} max={-70} step={1} />
                  </>
                )}
              </div>
            </div>

            <button onClick={runSim} disabled={running || !camValid}
              className="mt-5 w-full py-3 rounded-lg text-sm uppercase tracking-widest font-semibold"
              style={{ background: running ? T.panel2 : T.cyan, color: running ? T.lo : T.void, border: `1px solid ${T.hair}` }}>
              {running ? "Solving cycles…" : "Run full-cycle simulation"}
            </button>
          </div>

          {/* right: results */}
          <div className="col-span-12 md:col-span-8">
            {!result && (
              <div className="rounded-lg p-10 text-center" style={{ background: T.panel, border: `1px solid ${T.hair}`, color: T.lo }}>
                <div className="text-sm">Configure the engine, then run the simulation.</div>
                <div className="text-[11px] mt-2 font-mono">
                  Solves the full 720° cycle over multiple iterations until periodic —
                  everything (P@EVO, peak pressure, IMEP, power) emerges from first principles.
                </div>
              </div>
            )}
            {result && (
              <>
                <div className="grid grid-cols-2 md:grid-cols-4 gap-3 mb-5">
                  <Stat label="Brake power ×6" value={`${result.brakeHp6.toFixed(0)}hp`}
                        sub={`indicated ${(6 * result.hpPerCyl).toFixed(0)}hp − friction`} color={T.cyan} />
                  <Stat label="Brake torque ×6" value={`${result.brakeTorque6.toFixed(0)}Nm`}
                        sub={`BMEP ${result.bmepBar.toFixed(2)}bar`} color={T.green} />
                  <Stat label="FMEP (Chen-Flynn)" value={`${result.fmepBar.toFixed(2)}bar`}
                        sub={`mech eff ${result.mechEff.toFixed(1)}% · Sp ${result.SpMean.toFixed(1)}m/s`} color={T.amber} />
                  <Stat label="Peak P / P@EVO" value={`${result.PpeakBar.toFixed(1)}bar`}
                        sub={`EVO ${result.PatEVO.toFixed(2)}bar · air ${result.airMg.toFixed(0)}mg`} />
                </div>

                <div className="rounded-lg p-4 mb-5 overflow-x-auto" style={{ background: T.panel, border: `1px solid ${T.hair}` }}>
                  <div className="text-xs uppercase tracking-wider mb-2" style={{ color: T.lo }}>Cycle convergence</div>
                  <table className="w-full text-[12px] font-mono" style={{ color: T.hi }}>
                    <thead>
                      <tr style={{ color: T.lo }}>
                        <th className="text-left pr-4">cyc</th><th className="text-right pr-4">air (mg)</th>
                        <th className="text-right pr-4">IMEP (bar)</th><th className="text-right pr-4">P@EVO</th>
                        <th className="text-right">Ppeak</th>
                      </tr>
                    </thead>
                    <tbody>
                      {result.cycles.map((c) => (
                        <tr key={c.cycle}>
                          <td className="pr-4">{c.cycle}</td>
                          <td className="text-right pr-4">{c.airMg.toFixed(1)}</td>
                          <td className="text-right pr-4">{c.imepBar.toFixed(2)}</td>
                          <td className="text-right pr-4">{c.PatEVO.toFixed(2)}</td>
                          <td className="text-right">{c.PpeakBar.toFixed(1)}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>

                <div className="rounded-lg p-5 mb-5" style={{ background: T.panel, border: `1px solid ${T.hair}` }}>
                  <div className="text-xs uppercase tracking-wider mb-3" style={{ color: T.lo }}>
                    Cylinder pressure, full cycle (log scale) — 360° = firing TDC
                  </div>
                  <ResponsiveContainer width="100%" height={240}>
                    <LineChart data={result.hist} margin={{ top: 5, right: 15, bottom: 5, left: 0 }}>
                      <CartesianGrid stroke={T.line} strokeDasharray="2 4" />
                      <XAxis dataKey="theta" stroke={T.lo} tick={{ fontSize: 11, fill: T.lo }} />
                      <YAxis dataKey="logP" stroke={T.lo} tick={{ fontSize: 11, fill: T.lo }}
                             tickFormatter={(v) => `${Math.pow(10, v).toFixed(1)}`}
                             label={{ value: "bar (log)", angle: -90, position: "insideLeft", fill: T.lo, fontSize: 11 }} />
                      <Tooltip contentStyle={{ background: T.panel2, border: `1px solid ${T.hair}`, fontSize: 12 }}
                               labelStyle={{ color: T.lo }}
                               formatter={(v, n) => n === "logP" ? [`${Math.pow(10, v).toFixed(2)} bar`, "P"] : [v, n]} />
                      <ReferenceLine x={360} stroke={T.amber} strokeDasharray="4 4"
                                      label={{ value: "TDC fire", fill: T.amber, fontSize: 10 }} />
                      <Line type="monotone" dataKey="logP" stroke={T.green} strokeWidth={2} dot={false} name="Cylinder P" />
                    </LineChart>
                  </ResponsiveContainer>
                </div>

                <div className="rounded-lg p-5 mb-5" style={{ background: T.panel, border: `1px solid ${T.hair}` }}>
                  <div className="text-xs uppercase tracking-wider mb-3" style={{ color: T.lo }}>
                    Wave-affected port pressures at the valves
                  </div>
                  <ResponsiveContainer width="100%" height={200}>
                    <LineChart data={result.hist} margin={{ top: 5, right: 15, bottom: 5, left: 0 }}>
                      <CartesianGrid stroke={T.line} strokeDasharray="2 4" />
                      <XAxis dataKey="theta" stroke={T.lo} tick={{ fontSize: 11, fill: T.lo }} />
                      <YAxis stroke={T.lo} tick={{ fontSize: 11, fill: T.lo }}
                             label={{ value: "bar", angle: -90, position: "insideLeft", fill: T.lo, fontSize: 11 }} />
                      <Tooltip contentStyle={{ background: T.panel2, border: `1px solid ${T.hair}`, fontSize: 12 }}
                               labelStyle={{ color: T.lo }} />
                      <Legend wrapperStyle={{ fontSize: 11 }} />
                      <Line type="monotone" dataKey="Pi" stroke={T.cyan} strokeWidth={1.6} dot={false} name="Intake port @ valve" />
                      <Line type="monotone" dataKey="Pe" stroke={T.red} strokeWidth={1.6} dot={false} name="Exhaust port @ valve" />
                    </LineChart>
                  </ResponsiveContainer>
                </div>

                <div className="rounded-lg p-5" style={{ background: T.panel, border: `1px solid ${T.hair}` }}>
                  <div className="text-xs uppercase tracking-wider mb-3" style={{ color: T.lo }}>
                    True local Mach at each valve (wave-coupled)
                  </div>
                  <ResponsiveContainer width="100%" height={200}>
                    <LineChart data={result.hist} margin={{ top: 5, right: 15, bottom: 5, left: 0 }}>
                      <CartesianGrid stroke={T.line} strokeDasharray="2 4" />
                      <XAxis dataKey="theta" stroke={T.lo} tick={{ fontSize: 11, fill: T.lo }}
                             label={{ value: "Crank angle (°, 360=firing TDC)", position: "insideBottom", offset: -3, fill: T.lo, fontSize: 11 }} />
                      <YAxis stroke={T.lo} tick={{ fontSize: 11, fill: T.lo }} domain={[-1.05, 1.05]} />
                      <Tooltip contentStyle={{ background: T.panel2, border: `1px solid ${T.hair}`, fontSize: 12 }}
                               labelStyle={{ color: T.lo }} />
                      <Legend wrapperStyle={{ fontSize: 11 }} />
                      <ReferenceLine y={0} stroke={T.hair} />
                      <ReferenceLine y={1} stroke={T.red} strokeDasharray="4 4" />
                      <ReferenceLine y={-1} stroke={T.red} strokeDasharray="4 4" />
                      <Line type="monotone" dataKey="eM" stroke={T.red} strokeWidth={2} dot={false} name="Exhaust Mach" />
                      <Line type="monotone" dataKey="iM" stroke={T.cyan} strokeWidth={2} dot={false} name="Intake Mach" />
                    </LineChart>
                  </ResponsiveContainer>
                </div>

                <div className="mt-4 text-[11px] font-mono leading-relaxed" style={{ color: T.lo }}>
                  Complete 720° closed-loop cycle: 1-D finite-volume Euler runners (HLL, CFL-stepped,
                  Sod-validated) run continuously through all four strokes — waves keep ringing while
                  valves are closed, which is what makes runner tuning real. Cylinder: 0-D mass/energy
                  with piston work and single-zone Wiebe heat release (fuel = inducted air / AFR).
                  Iterated cycle-over-cycle until periodic; the residual gas, P@EVO, blowdown strength,
                  and induction all feed back on each other and settle together. Simplifications:
                  single species (γ=1.38), adiabatic walls (indicated figures run slightly high),
                  no pipe friction, straight runners, infinite plenum/collector.
                </div>
              </>
            )}
          </div>
        </div>
      </div>
    </div>
  );
}
