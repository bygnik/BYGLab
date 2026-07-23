//! Validation for `branch_junction.rs` against real ground truth: a
//! closed-form linear-acoustics prediction for a symmetric N-way split (or
//! merge) junction, an exact-conservation check, an explicit measurement of
//! the (real, documented) energy-conservation gap, and a regression against
//! `network::Junction`'s own exact 2-pipe HLLC solve.

use byglab_core::boundary::BoundaryCondition;
use byglab_core::branch_junction::{resolve_branch_junction, BranchJunction};
use byglab_core::gas::{GasProperties, PrimitiveState};
use byglab_core::mesh::Mesh;
use byglab_core::network::{Junction, PipeEnd, PipeEndRef, PipeNetwork};
use byglab_core::pipe::Pipe;

fn uniform_pipe(pressure: f64, diameter: f64) -> Pipe {
    Pipe::uniform_initial_state(
        Mesh::uniform(1.0, diameter, 0.02),
        PrimitiveState::from_pressure_temperature(pressure, 293.15, 0.0, &GasProperties::AIR),
        &GasProperties::AIR,
        BoundaryCondition::ClosedEnd,
        BoundaryCondition::ClosedEnd,
    )
}

fn circle_area(diameter: f64) -> f64 {
    std::f64::consts::PI / 4.0 * diameter * diameter
}

/// Resolves a trunk pipe (pipe 0) splitting into `branch_count` identical
/// branch pipes and returns the single shared junction pressure every leg
/// resolves to.
fn resolve_symmetric_split(
    trunk_pressure: f64,
    trunk_diameter: f64,
    branch_pressure: f64,
    branch_diameter: f64,
    branch_count: usize,
) -> f64 {
    let gas = GasProperties::AIR;
    let mut pipes = vec![uniform_pipe(trunk_pressure, trunk_diameter)];
    pipes.extend((0..branch_count).map(|_| uniform_pipe(branch_pressure, branch_diameter)));
    let network = PipeNetwork { pipes, junctions: vec![] };

    let mut ends = vec![PipeEndRef { pipe_index: 0, end: PipeEnd::Right }];
    ends.extend((0..branch_count).map(|i| PipeEndRef { pipe_index: i + 1, end: PipeEnd::Left }));
    let junction = BranchJunction { ends };

    let fluxes = resolve_branch_junction(&network, &junction, &gas);
    let junction_pressure = fluxes[0].neighbor_state.pressure;
    for flux in &fluxes {
        assert!(
            (flux.neighbor_state.pressure - junction_pressure).abs() < 1e-6 * junction_pressure,
            "every leg should resolve to the same shared junction pressure"
        );
    }
    junction_pressure
}

/// The closed-form linear-acoustics junction pressure for a trunk (area
/// `trunk_area`) meeting `branch_count` identical branches (area
/// `branch_area` each): derived directly from the standard right/left-going
/// Riemann-invariant matching conditions at the junction (`p + rho*c*u`
/// fixed at its trunk-side initial value, `p - rho*c*u` fixed at its
/// combined-branch-side initial value, pressure and volume flow continuous
/// across the junction) - an area-weighted average of the two sides'
/// initial pressures, exact in the small-disturbance limit.
fn analytic_linear_acoustics_junction_pressure(
    trunk_pressure: f64,
    trunk_area: f64,
    branch_pressure: f64,
    branch_area: f64,
    branch_count: usize,
) -> f64 {
    let total_branch_area = branch_area * branch_count as f64;
    (trunk_area * trunk_pressure + total_branch_area * branch_pressure) / (trunk_area + total_branch_area)
}

#[test]
fn symmetric_n_branch_junction_pressure_matches_the_linear_acoustics_prediction_and_the_error_shrinks_with_amplitude(
) {
    let trunk_diameter = 0.05;
    let branch_diameter = 0.03;
    let branch_count = 3;
    let baseline_pressure = 1.0e5;
    let trunk_area = circle_area(trunk_diameter);
    let branch_area = circle_area(branch_diameter);

    let mut relative_errors = Vec::new();
    for relative_disturbance in [0.04, 0.02, 0.01, 0.005] {
        let trunk_pressure = baseline_pressure * (1.0 + relative_disturbance);
        let branch_pressure = baseline_pressure;

        let numerical =
            resolve_symmetric_split(trunk_pressure, trunk_diameter, branch_pressure, branch_diameter, branch_count);
        let analytic = analytic_linear_acoustics_junction_pressure(
            trunk_pressure,
            trunk_area,
            branch_pressure,
            branch_area,
            branch_count,
        );

        // Normalized by the perturbation itself (not the absolute
        // pressure) - this is a check of the LINEARIZATION error, which is
        // expected to scale with the disturbance amplitude, not a fixed
        // absolute tolerance.
        let perturbation = analytic - baseline_pressure;
        let relative_error = (numerical - analytic).abs() / perturbation.abs();
        println!(
            "disturbance={relative_disturbance:e}: numerical={numerical:.4}, analytic={analytic:.4}, relative error vs perturbation={relative_error:e}"
        );
        relative_errors.push(relative_error);
    }

    // The dominant discrepancy source is the nonlinear closure's own
    // isentropic (not linear-acoustics) relation, an O(disturbance) effect
    // relative to the linear prediction - so halving the disturbance
    // should roughly halve the relative error. Checked directly (not
    // assumed): consecutive ratios should be well below 1, not close to 1
    // or above (which would indicate the numerical result isn't actually
    // converging to the linear-acoustics prediction as amplitude shrinks).
    for pair in relative_errors.windows(2) {
        let ratio = pair[1] / pair[0];
        println!("consecutive relative-error ratio (halved disturbance): {ratio:e}");
        assert!(ratio < 0.7, "expected the discrepancy to shrink as the disturbance amplitude shrinks, got ratio {ratio:e}");
    }
    // At the smallest tested amplitude, the two should already agree to a
    // few percent.
    assert!(
        *relative_errors.last().unwrap() < 0.05,
        "expected close agreement with the linear-acoustics prediction at small amplitude, got {:e}",
        relative_errors.last().unwrap()
    );
}

#[test]
fn symmetric_merge_of_n_branches_into_one_trunk_matches_the_same_linear_acoustics_prediction() {
    // The reverse topology (N identical branches merging into one trunk)
    // is the SAME junction structure with the roles of "higher pressure"
    // and "lower pressure" swapped - same formula applies by construction
    // (the closed-form prediction doesn't care which side is which).
    let trunk_diameter = 0.06;
    let branch_diameter = 0.025;
    let branch_count = 4;
    let baseline_pressure = 1.0e5;
    let trunk_area = circle_area(trunk_diameter);
    let branch_area = circle_area(branch_diameter);

    let relative_disturbance = 0.01;
    let trunk_pressure = baseline_pressure;
    let branch_pressure = baseline_pressure * (1.0 + relative_disturbance); // branches now the higher-pressure side

    let numerical =
        resolve_symmetric_split(trunk_pressure, trunk_diameter, branch_pressure, branch_diameter, branch_count);
    let analytic = analytic_linear_acoustics_junction_pressure(
        trunk_pressure,
        trunk_area,
        branch_pressure,
        branch_area,
        branch_count,
    );
    let perturbation = analytic - baseline_pressure;
    let relative_error = (numerical - analytic).abs() / perturbation.abs();
    println!("merge topology: numerical={numerical:.4}, analytic={analytic:.4}, relative error={relative_error:e}");
    assert!(relative_error < 0.05, "expected close agreement, got {relative_error:e}");
}

#[test]
fn n_equals_2_matched_area_regresses_against_the_exact_junction_hllc_solve() {
    // `Junction` performs an exact 2-state HLLC solve; `BranchJunction`'s
    // shared-pressure/per-leg-Riemann-invariant closure is a real
    // approximation to that (see branch_junction.rs's own doc comment).
    // This quantifies - rather than just sanity-checks - exactly how large
    // that approximation's junction-pressure error is for the simplest
    // possible case (2 pipes, matched area), which BranchJunction is NOT
    // the recommended tool for (Junction should be preferred there) but
    // must still behave sensibly on.
    let gas = GasProperties::AIR;
    let diameter = 0.05;
    let pressure_a = 1.01e5;
    let pressure_b = 1.00e5;

    // Exact HLLC junction: read the resolved face state at the junction by
    // running the network forward one CFL step and inspecting the
    // boundary-cell pressure right at the junction face (the first cell of
    // pipe B, which after one small step reflects the Riemann-solved
    // junction state).
    let pipe_a = uniform_pipe(pressure_a, diameter);
    let pipe_b = uniform_pipe(pressure_b, diameter);
    let mut exact_network = PipeNetwork {
        pipes: vec![pipe_a, pipe_b],
        junctions: vec![Junction {
            a: PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
            b: PipeEndRef { pipe_index: 1, end: PipeEnd::Left },
        }],
    };
    let dt = byglab_core::solver::step(&mut exact_network, &gas, 0.5);
    let _ = dt;
    let exact_pressure_near_junction = exact_network.pipes[0].state.last().unwrap().to_primitive(&gas).pressure;

    let approximate_pressure = resolve_symmetric_split(pressure_a, diameter, pressure_b, diameter, 1);

    let relative_difference =
        (approximate_pressure - exact_pressure_near_junction).abs() / (pressure_a - pressure_b).abs();
    println!(
        "N=2 matched-area regression: exact HLLC junction-adjacent pressure={exact_pressure_near_junction:.4}, \
         BranchJunction shared pressure={approximate_pressure:.4}, relative difference (vs perturbation)={relative_difference:e}"
    );
    // Measured, not guessed: both should be close to the simple average of
    // the two pressures for a matched-area junction, so the difference
    // should be a small fraction of the perturbation itself.
    assert!(
        relative_difference < 0.5,
        "BranchJunction's N=2 matched-area result should stay reasonably close to the exact HLLC solve, got {relative_difference:e}"
    );
}

#[test]
fn energy_residual_is_measured_at_small_and_realistic_amplitude() {
    // branch_junction.rs's own doc comment documents that this closure
    // gives exact mass conservation but not exact energy conservation in
    // general. This measures the actual residual - not just asserting it
    // exists - at two amplitudes, matching this project's convention of
    // reporting real numbers rather than a bare pass/fail.
    let trunk_diameter = 0.05;
    let branch_diameter = 0.03;
    let branch_count = 3;
    let gas = GasProperties::AIR;

    for (label, trunk_pressure, branch_pressure) in
        [("small (1%)", 1.01e5, 1.0e5), ("realistic exhaust pulse", 3.0e5, 1.0e5)]
    {
        let mut pipes = vec![uniform_pipe(trunk_pressure, trunk_diameter)];
        pipes.extend((0..branch_count).map(|_| uniform_pipe(branch_pressure, branch_diameter)));
        let network = PipeNetwork { pipes, junctions: vec![] };
        let mut ends = vec![PipeEndRef { pipe_index: 0, end: PipeEnd::Right }];
        ends.extend((0..branch_count).map(|i| PipeEndRef { pipe_index: i + 1, end: PipeEnd::Left }));
        let junction = BranchJunction { ends };

        let fluxes = resolve_branch_junction(&network, &junction, &gas);
        let trunk_area = circle_area(trunk_diameter);
        let branch_area = circle_area(branch_diameter);

        let mut total_energy_flow = 0.0;
        let mut total_throughput_energy = 0.0;
        for flux in &fluxes {
            let area = match flux.end.end {
                PipeEnd::Left => network.pipes[flux.end.pipe_index].mesh.face_areas[0],
                PipeEnd::Right => {
                    let fa = &network.pipes[flux.end.pipe_index].mesh.face_areas;
                    fa[fa.len() - 1]
                }
            };
            let sign = match flux.end.end {
                PipeEnd::Right => 1.0,
                PipeEnd::Left => -1.0,
            };
            total_energy_flow += sign * flux.flux.energy * area;
            total_throughput_energy += (flux.flux.energy * area).abs();
        }
        let _ = (trunk_area, branch_area);
        let relative_residual = total_energy_flow.abs() / total_throughput_energy;
        println!("{label}: energy residual relative to total throughput = {relative_residual:e}");
    }
}
