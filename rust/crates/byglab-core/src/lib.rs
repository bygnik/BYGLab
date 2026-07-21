//! Placeholder solver core for the Rust/wasm integration smoke test.
//!
//! This crate is the single source of truth shared by the CLI and the wasm
//! bindings — see `byglab-wasm` and `byglab-cli`. `run()` here is a stand-in
//! for the eventual 1D finite-volume solver; it exists only to prove the
//! config-in/time-series-out shape end to end before any real physics lands.

use serde::{Deserialize, Serialize};
use std::f64::consts::TAU;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimConfig {
    pub label: String,
    pub n_points: u32,
    pub amplitude: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimResult {
    pub theta: Vec<f64>,
    pub value: Vec<f64>,
}

/// Samples `amplitude * sin(theta)` over `n_points` evenly-spaced points
/// covering `[0, 2*pi]` inclusive.
pub fn run(config: &SimConfig) -> SimResult {
    let n = config.n_points as usize;
    if n == 0 {
        return SimResult { theta: vec![], value: vec![] };
    }
    if n == 1 {
        return SimResult { theta: vec![0.0], value: vec![0.0] };
    }

    let step = TAU / (n - 1) as f64;
    let theta: Vec<f64> = (0..n).map(|i| i as f64 * step).collect();
    let value: Vec<f64> = theta.iter().map(|t| config.amplitude * t.sin()).collect();

    SimResult { theta, value }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg(n_points: u32, amplitude: f64) -> SimConfig {
        SimConfig { label: "test".into(), n_points, amplitude }
    }

    #[test]
    fn length_matches_n_points() {
        let result = run(&cfg(200, 1.0));
        assert_eq!(result.theta.len(), 200);
        assert_eq!(result.value.len(), 200);
    }

    #[test]
    fn starts_and_ends_at_zero_crossing() {
        let result = run(&cfg(200, 2.5));
        assert!((result.theta[0] - 0.0).abs() < 1e-12);
        assert!((result.value[0] - 0.0).abs() < 1e-9);
        assert!((result.theta[199] - TAU).abs() < 1e-9);
        // sin(2*pi) == 0
        assert!((result.value[199] - 0.0).abs() < 1e-6);
    }

    #[test]
    fn amplitude_scales_linearly() {
        let a = run(&cfg(50, 1.0));
        let b = run(&cfg(50, 3.0));
        for (va, vb) in a.value.iter().zip(b.value.iter()) {
            assert!((vb - 3.0 * va).abs() < 1e-9);
        }
    }

    #[test]
    fn zero_and_one_point_edge_cases_do_not_panic() {
        assert_eq!(run(&cfg(0, 1.0)).theta.len(), 0);
        assert_eq!(run(&cfg(1, 1.0)).theta.len(), 1);
    }

    #[test]
    fn config_round_trips_through_json() {
        let config = cfg(10, 1.5);
        let json = serde_json::to_string(&config).unwrap();
        let back: SimConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n_points, 10);
        assert!((back.amplitude - 1.5).abs() < 1e-12);
    }
}
