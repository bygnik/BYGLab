//! Shared setup for the Sod shock tube case (10 bar vs 1 bar air, two 1m
//! pipes joined at a junction), parameterized by mesh resolution. Used by
//! both the fixed-resolution validation test (`tests/sod_shock_tube.rs`)
//! and the mesh-refinement convergence check (`tests/mesh_convergence.rs`)
//! so the two stay in sync and neither duplicates the setup/error-sampling
//! logic.

use super::exact_riemann::RiemannProblem;
use byglab_core::boundary::BoundaryCondition;
use byglab_core::gas::{GasProperties, PrimitiveState};
use byglab_core::mesh::Mesh;
use byglab_core::network::{Junction, PipeEnd, PipeEndRef, PipeNetwork};
use byglab_core::pipe::Pipe;
use byglab_core::solver::run_to_time;

/// 0.9 ms — matches the OpenWAM comparison point, safely before either
/// wave reaches a closed outer wall (shock at ~1.81 ms, rarefaction head
/// at ~2.91 ms), so the exact infinite-domain Riemann solution still
/// applies in the interior region between the wavefronts.
pub const TARGET_TIME: f64 = 0.0009;

/// x-position (within `pipe_b`) of the star-region accuracy probe — well
/// clear of both wavefronts at `TARGET_TIME`, in the flat plateau where
/// pressure and velocity are constant (p*/u*) on both sides of the contact.
pub const STAR_REGION_PROBE_X: f64 = 0.2;

/// The result of running the Sod shock tube case: the final network state,
/// the simulated time actually reached, and the exact Riemann solution to
/// compare it against.
pub struct SodShockTubeRun {
    pub network: PipeNetwork,
    pub elapsed: f64,
    pub exact: RiemannProblem,
}

/// RMS and max absolute error over every cell in both pipes, for one flow
/// quantity (pressure or velocity).
pub struct WholeProfileError {
    pub rms: f64,
    pub max: f64,
}

/// Builds and runs the Sod shock tube case with a mesh targeting
/// `target_cell_size` meters per cell in each pipe, to [`TARGET_TIME`].
pub fn run(target_cell_size: f64) -> SodShockTubeRun {
    let gas = GasProperties::AIR;

    let pipe_a = Pipe::uniform_initial_state(
        Mesh::uniform(1.0, 0.05, target_cell_size),
        PrimitiveState::from_pressure_temperature(10e5, 293.15, 0.0, &gas),
        &gas,
        BoundaryCondition::ClosedEnd,
        BoundaryCondition::ClosedEnd, // right end overridden by the junction below
    );
    let pipe_b = Pipe::uniform_initial_state(
        Mesh::uniform(1.0, 0.05, target_cell_size),
        PrimitiveState::from_pressure_temperature(1e5, 293.15, 0.0, &gas),
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

    let elapsed = run_to_time(&mut network, &gas, 0.5, TARGET_TIME);

    let density_left = 10e5 / (gas.gas_constant * 293.15);
    let density_right = 1e5 / (gas.gas_constant * 293.15);
    let exact = RiemannProblem::new(density_left, 0.0, 10e5, density_right, 0.0, 1e5, gas.gamma);

    SodShockTubeRun { network, elapsed, exact }
}

impl SodShockTubeRun {
    /// Absolute pressure error at [`STAR_REGION_PROBE_X`] — the smooth,
    /// unsmeared part of the flow, where accuracy is expected to be much
    /// higher than the whole-profile error (which includes the smeared
    /// shock/rarefaction fronts).
    pub fn star_region_pressure_error(&self) -> f64 {
        (self.probe_state().pressure - self.exact.p_star).abs()
    }

    /// Absolute velocity error at [`STAR_REGION_PROBE_X`]. See
    /// [`Self::star_region_pressure_error`].
    pub fn star_region_velocity_error(&self) -> f64 {
        (self.probe_state().velocity - self.exact.u_star).abs()
    }

    fn probe_state(&self) -> PrimitiveState {
        let gas = GasProperties::AIR;
        let cell_width = self.network.pipes[1].mesh.cells[0].width;
        let probe_cell_index = (STAR_REGION_PROBE_X / cell_width) as usize;
        self.network.pipes[1].state[probe_cell_index].to_primitive(&gas)
    }

    /// Whole-profile pressure error (RMS and max), sampled at every cell
    /// center in both pipes against the exact solution.
    pub fn whole_profile_pressure_error(&self) -> WholeProfileError {
        self.whole_profile_error(|primitive| primitive.pressure, |_, _, p| p)
    }

    /// Whole-profile velocity error (RMS and max). See
    /// [`Self::whole_profile_pressure_error`].
    pub fn whole_profile_velocity_error(&self) -> WholeProfileError {
        self.whole_profile_error(|primitive| primitive.velocity, |_, u, _| u)
    }

    /// Shared sampling loop for both whole-profile error methods: walks
    /// every cell in both pipes (mapping each pipe's local cell centers to
    /// the shared global x — `pipe_a` occupies `[-1, 0]`, `pipe_b`
    /// occupies `[0, 1]`, with the initial discontinuity at the junction,
    /// x=0), extracting one scalar quantity from the solver's state
    /// (`extract_actual`) and from the exact solution's `(rho, u, p)`
    /// tuple (`extract_exact`) to compare.
    fn whole_profile_error(
        &self,
        extract_actual: impl Fn(PrimitiveState) -> f64,
        extract_exact: impl Fn(f64, f64, f64) -> f64,
    ) -> WholeProfileError {
        let gas = GasProperties::AIR;
        let mut squared_errors = Vec::new();
        let mut max_error = 0.0_f64;

        let mut sample_pipe = |cells: &[byglab_core::Cell], states: &[byglab_core::ConservedState], offset: f64| {
            for (cell, state) in cells.iter().zip(states.iter()) {
                let global_x = cell.center + offset;
                let (rho, u, p) = self.exact.sample(global_x, self.elapsed, 0.0);
                let exact_value = extract_exact(rho, u, p);
                let actual_value = extract_actual(state.to_primitive(&gas));
                let error = actual_value - exact_value;
                squared_errors.push(error * error);
                max_error = max_error.max(error.abs());
            }
        };

        sample_pipe(&self.network.pipes[0].mesh.cells, &self.network.pipes[0].state, -1.0);
        sample_pipe(&self.network.pipes[1].mesh.cells, &self.network.pipes[1].state, 0.0);

        let rms = (squared_errors.iter().sum::<f64>() / squared_errors.len() as f64).sqrt();
        WholeProfileError { rms, max: max_error }
    }
}
