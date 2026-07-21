const API_BASE = import.meta.env.VITE_API_BASE || "";

export async function fetchPresetS54() {
  const res = await fetch(`${API_BASE}/api/v1/presets/s54b32`);
  if (!res.ok) throw new Error(`Failed to load preset (${res.status})`);
  return res.json();
}

export async function simulateWave(inputs) {
  const res = await fetch(`${API_BASE}/api/v1/wave/simulate`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(inputs),
  });
  if (!res.ok) {
    let detail = `Simulation failed (${res.status})`;
    try {
      const body = await res.json();
      detail = body.detail || detail;
    } catch {
      /* ignore */
    }
    throw new Error(typeof detail === "string" ? detail : JSON.stringify(detail));
  }
  return res.json();
}

export async function fetchModules() {
  const res = await fetch(`${API_BASE}/api/v1/modules`);
  if (!res.ok) throw new Error(`Failed to load modules (${res.status})`);
  return res.json();
}
