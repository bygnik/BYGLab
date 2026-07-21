import { NavLink, Outlet } from "react-router-dom";

const MODULES = [
  { to: "/wave", label: "Wave Dynamics", ready: true },
  { to: "/smoke-test", label: "WASM Smoke Test", ready: true },
  { to: "/flowbench", label: "Flowbench", ready: false },
  { to: "/dyno", label: "Dyno", ready: false },
  { to: "/turbo", label: "Turbo", ready: false },
  { to: "/camshaft", label: "Camshaft", ready: false },
  { to: "/exhaust", label: "Exhaust Header", ready: false },
  { to: "/intake", label: "Intake", ready: false },
  { to: "/combustion", label: "Combustion", ready: false },
];

export default function Layout() {
  return (
    <div className="min-h-screen font-sans text-base md:text-lg text-hi">
      <header className="border-b border-hair/80 bg-panel/80 backdrop-blur-sm sticky top-0 z-20">
        <div className="max-w-7xl mx-auto px-5 py-4 flex flex-col gap-3 md:flex-row md:items-end md:justify-between">
          <div className="flex items-center gap-4">
            <img
              src="/byg_racing_logo1.JPEG"
              alt="BYG Racing logo"
              className="h-12 w-auto rounded-md border border-hair bg-panel2 p-1 object-contain"
            />
            <div>
              <div className="font-display text-3xl md:text-4xl tracking-[0.08em] uppercase text-cyan leading-none">
                BYGLab
              </div>
              <div className="text-xs md:text-sm tracking-[0.18em] uppercase text-lo mt-1">
                Engine modeling suite · 1-D wave dynamics
              </div>
            </div>
          </div>
          <nav className="flex flex-wrap gap-1.5">
            {MODULES.map((m) =>
              m.ready ? (
                <NavLink
                  key={m.to}
                  to={m.to}
                  className={({ isActive }) =>
                    `px-2.5 py-1.5 text-xs md:text-sm uppercase tracking-wider rounded border transition-colors ${
                      isActive
                        ? "bg-cyan text-void border-cyan"
                        : "border-hair text-xs md:text-sm text-lo hover:text-hi hover:border-lo"
                    }`
                  }
                >
                  {m.label}
                </NavLink>
              ) : (
                <span
                  key={m.to}
                  title="Coming soon"
                  className="px-2.5 py-1.5 text-xs md:text-sm uppercase tracking-wider rounded border border-hair/50 text-lo/50 cursor-not-allowed"
                >
                  {m.label}
                </span>
              )
            )}
          </nav>
        </div>
      </header>
      <main className="max-w-7xl mx-auto px-5 py-6">
        <Outlet />
      </main>
    </div>
  );
}
