import { Navigate, Route, Routes } from "react-router-dom";
import Layout from "./components/Layout";
import WavePage from "./pages/WavePage";
import SmokeTestPage from "./pages/SmokeTestPage";
import ComingSoon from "./pages/ComingSoon";

export default function App() {
  return (
    <Routes>
      <Route element={<Layout />}>
        <Route index element={<Navigate to="/wave" replace />} />
        <Route path="wave" element={<WavePage />} />
        <Route path="smoke-test" element={<SmokeTestPage />} />
        <Route path=":moduleId" element={<ComingSoon />} />
      </Route>
    </Routes>
  );
}
