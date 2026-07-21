// Web Worker hosting the byglab-wasm module so simulation runs never block
// the UI thread. Protocol:
//   main -> worker: { type: 'run', requestId, config }
//   worker -> main: { type: 'result', requestId, series: { theta, value } }
//                 | { type: 'error', requestId, message }
//
// Single-threaded wasm-bindgen build (no SharedArrayBuffer, no COOP/COEP
// headers required) — see rust/crates/byglab-wasm.

import init, { run as wasmRun } from "../../../rust/crates/byglab-wasm/pkg/byglab_wasm.js";

const ready = init();

self.onmessage = async (event) => {
  const { type, requestId, config } = event.data;
  if (type !== "run") return;

  try {
    await ready;
    const result = wasmRun(JSON.stringify(config));
    self.postMessage({
      type: "result",
      requestId,
      series: { theta: result.theta, value: result.value },
    });
  } catch (err) {
    self.postMessage({
      type: "error",
      requestId,
      message: err?.message || String(err),
    });
  }
};
