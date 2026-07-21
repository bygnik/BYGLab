import { useParams } from "react-router-dom";

const TITLES = {
  flowbench: "Flowbench",
  dyno: "Dyno",
  turbo: "Turbo",
  camshaft: "Camshaft",
  exhaust: "Exhaust Header",
  intake: "Intake",
  combustion: "Combustion",
};

export default function ComingSoon() {
  const { moduleId } = useParams();
  const title = TITLES[moduleId] || "Module";

  return (
    <div className="rounded-lg border border-hair bg-panel p-10 text-center">
      <div className="font-display text-2xl tracking-wide uppercase text-hi">{title}</div>
      <p className="mt-3 text-sm text-lo max-w-md mx-auto">
        This module is scaffolded for the BYGLab suite and will ship in a later release.
        Use <span className="text-cyan">Wave Dynamics</span> for the 1-D closed-loop cycle solver.
      </p>
    </div>
  );
}
