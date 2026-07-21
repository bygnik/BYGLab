import { useEffect, useRef, useState } from "react";
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer,
} from "recharts";

const T = {
  panel: "#101623", panel2: "#152138", line: "#2C3A5C", hair: "#273149",
  cyan: "#3A7CF6", danger: "#E5644B", hi: "#EBF1FF", lo: "#98A6C0",
};

function Field({ label, value, onChange, min, max, step }) {
  return (
    <div className="mb-4">
      <div className="flex items-baseline justify-between mb-1">
        <label className="text-xs tracking-wide uppercase text-lo">{label}</label>
        <span className="text-sm font-mono text-cyan">{value}</span>
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

// Proves the Rust -> wasm -> Web Worker -> React pipeline end to end before
// any real solver physics is ported. See rust/crates/byglab-{core,wasm,cli}
// and services/web/src/wasm-worker.js.
export default function SmokeTestPage() {
  const [nPoints, setNPoints] = useState(200);
  const [amplitude, setAmplitude] = useState(1.0);
  const [series, setSeries] = useState(null);
  const [running, setRunning] = useState(false);
  const [error, setError] = useState(null);
  const workerRef = useRef(null);
  const requestIdRef = useRef(0);

  useEffect(() => {
    const worker = new Worker(new URL("../wasm-worker.js", import.meta.url), { type: "module" });
    worker.onmessage = (event) => {
      const { type, requestId, series: resultSeries, message } = event.data;
      if (requestId !== requestIdRef.current) return; // stale response, ignore
      setRunning(false);
      if (type === "result") {
        setError(null);
        const chartData = Array.from(resultSeries.theta).map((theta, i) => ({
          theta,
          value: resultSeries.value[i],
        }));
        setSeries(chartData);
      } else if (type === "error") {
        setError(message);
      }
    };
    workerRef.current = worker;
    return () => worker.terminate();
  }, []);

  const runSimulation = () => {
    requestIdRef.current += 1;
    setRunning(true);
    setError(null);
    workerRef.current.postMessage({
      type: "run",
      requestId: requestIdRef.current,
      config: { label: "smoke", n_points: Math.round(nPoints), amplitude },
    });
  };

  return (
    <div className="rounded-lg border border-hair bg-panel p-6">
      <div className="font-display text-2xl tracking-wide uppercase text-hi">WASM Smoke Test</div>
      <p className="mt-2 text-sm text-lo max-w-xl">
        Exercises the full pipeline — Rust (<code className="text-cyan">byglab-core</code>) compiled
        to WebAssembly, running in a Web Worker, driven from React — before any real solver physics
        is ported. See <code className="text-cyan">rust/</code>.
      </p>

      <div className="mt-6 grid grid-cols-1 md:grid-cols-3 gap-6">
        <div className="rounded-md border border-hair bg-panel2 p-4">
          <Field label="Points" value={nPoints} onChange={setNPoints} min={2} max={2000} step={1} />
          <Field label="Amplitude" value={amplitude} onChange={setAmplitude} min={0.1} max={5} step={0.1} />
          <button
            type="button"
            onClick={runSimulation}
            disabled={running}
            className="w-full px-3 py-2 text-xs uppercase tracking-wider rounded border border-hair bg-cyan text-void font-semibold disabled:opacity-50"
          >
            {running ? "Running..." : "Run"}
          </button>
          {error && (
            <div className="mt-3 text-xs text-danger break-words">{error}</div>
          )}
        </div>

        <div className="md:col-span-2 rounded-md border border-hair bg-panel2 p-4">
          {series ? (
            <ResponsiveContainer width="100%" height={280}>
              <LineChart data={series} margin={{ top: 10, right: 15, bottom: 5, left: 0 }}>
                <CartesianGrid stroke={T.line} strokeDasharray="2 4" />
                <XAxis dataKey="theta" stroke={T.lo} tick={{ fontSize: 13, fill: T.lo }} />
                <YAxis stroke={T.lo} tick={{ fontSize: 13, fill: T.lo }} />
                <Tooltip contentStyle={{ background: T.panel, border: `1px solid ${T.hair}` }} />
                <Line type="monotone" dataKey="value" stroke={T.cyan} dot={false} strokeWidth={2} />
              </LineChart>
            </ResponsiveContainer>
          ) : (
            <div className="h-[280px] flex items-center justify-center text-sm text-lo">
              Run the simulation to see wasm-sourced output rendered here.
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
