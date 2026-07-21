//! Two 1m pipes joined at a shared node, closed at both outer ends, a small
//! pressure step (1.05 bar vs 1.00 bar) — small enough to stay in the
//! linear-acoustics regime (no shock forms). A closed-closed tube's
//! fundamental resonance period `T = 2L/c0` is exact regardless of the
//! excitation's shape (every harmonic is an integer multiple of the same
//! fundamental), so this checks the solver's phase/dispersion accuracy
//! over many wave transits, not just one.
//!
//! Matches `benchmarks/openwam/cases/acoustic_resonance/`, where OpenWAM's
//! higher-order TVD scheme measured the period to within 0.07% of exact.
//! This solver is first-order (see `solver.rs`'s doc comment), so more
//! numerical dispersion accumulating over many wave transits is expected —
//! the tolerance here is looser accordingly.

use byglab_core::boundary::BoundaryCondition;
use byglab_core::gas::{GasProperties, PrimitiveState};
use byglab_core::mesh::Mesh;
use byglab_core::network::{Junction, PipeEnd, PipeEndRef, PipeNetwork};
use byglab_core::pipe::Pipe;
use byglab_core::solver::step;

#[test]
fn resonance_period_matches_the_closed_form_prediction() {
    let gas = GasProperties::AIR;

    let pipe_a = Pipe::uniform_initial_state(
        Mesh::uniform(1.0, 0.05, 0.01),
        PrimitiveState::from_pressure_temperature(1.05e5, 293.15, 0.0, &gas),
        &gas,
        BoundaryCondition::ClosedEnd,
        BoundaryCondition::ClosedEnd, // right end overridden by the junction below
    );
    let pipe_b = Pipe::uniform_initial_state(
        Mesh::uniform(1.0, 0.05, 0.01),
        PrimitiveState::from_pressure_temperature(1.00e5, 293.15, 0.0, &gas),
        &gas,
        BoundaryCondition::ClosedEnd, // left end overridden by the junction below
        BoundaryCondition::ClosedEnd,
    );

    let mut network = PipeNetwork {
        pipes: vec![pipe_a, pipe_b],
        junctions: vec![Junction {
            a: PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
            b: PipeEndRef { pipe_index: 1, end: PipeEnd::Left },
        }],
    };

    let cfl = 0.5;
    let duration = 0.030; // 30 ms - several resonance periods
    let mut elapsed = 0.0;
    let mut pressure_history_at_left_wall: Vec<(f64, f64)> = Vec::new();

    while elapsed < duration {
        elapsed += step(&mut network, &gas, cfl);
        let pressure = network.pipes[0].state[0].to_primitive(&gas).pressure;
        pressure_history_at_left_wall.push((elapsed, pressure));
    }

    let observed_period = mean_period_from_threshold_crossings(&pressure_history_at_left_wall);

    let sound_speed = (gas.gamma * gas.gas_constant * 293.15_f64).sqrt();
    let total_closed_length = 2.0; // two 1m pipes
    let exact_period = 2.0 * total_closed_length / sound_speed;

    let relative_error = (observed_period - exact_period).abs() / exact_period;
    println!(
        "observed period {:.4} ms vs exact {:.4} ms ({:+.3}% error)",
        observed_period * 1000.0,
        exact_period * 1000.0,
        relative_error * 100.0 * (observed_period - exact_period).signum()
    );
    // Measured at ~0.09% with the MUSCL-Hancock (2nd-order) scheme -
    // still comfortably tighter than OpenWAM's own TVD-scheme result
    // (0.07%)... actually comparable, not tighter (0.002% with this
    // solver's earlier first-order version, before the MUSCL-Hancock
    // upgrade - the slope limiter's min/max clipping shifts numerical
    // dispersion slightly differently than plain upwinding for a small
    // linear disturbance; both are well within this margin either way).
    // 1% leaves comfortable room above the measured value while still
    // catching a real regression.
    assert!(
        relative_error < 0.01,
        "resonance period off by {:.2}%: observed {:.4} ms vs exact {:.4} ms",
        relative_error * 100.0,
        observed_period * 1000.0,
        exact_period * 1000.0
    );
}

/// Finds the mean time between consecutive downward crossings of the
/// pressure series through its own midpoint — the same threshold-crossing
/// technique `benchmarks/openwam/analysis/verify_acoustic_resonance.py`
/// uses to measure OpenWAM's resonance period, so the two are
/// apples-to-apples comparable.
fn mean_period_from_threshold_crossings(samples: &[(f64, f64)]) -> f64 {
    let max_pressure = samples.iter().map(|(_, p)| *p).fold(f64::MIN, f64::max);
    let min_pressure = samples.iter().map(|(_, p)| *p).fold(f64::MAX, f64::min);
    let midpoint = 0.5 * (max_pressure + min_pressure);

    let mut crossing_times = Vec::new();
    for pair in samples.windows(2) {
        let (t0, p0) = pair[0];
        let (t1, p1) = pair[1];
        if p0 >= midpoint && p1 < midpoint {
            let fraction = (midpoint - p0) / (p1 - p0);
            crossing_times.push(t0 + fraction * (t1 - t0));
        }
    }
    assert!(crossing_times.len() >= 2, "expected at least 2 downward crossings, found {}", crossing_times.len());

    let periods: Vec<f64> = crossing_times.windows(2).map(|pair| pair[1] - pair[0]).collect();
    periods.iter().sum::<f64>() / periods.len() as f64
}
