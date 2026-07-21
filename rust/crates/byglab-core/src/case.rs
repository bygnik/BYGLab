//! Serializable, JSON-friendly description of a pipe-network simulation and
//! its result. This is the stable public API surface `byglab-cli` and
//! `byglab-wasm` bind to — everything in here is plain data (no file I/O,
//! no threads), which is exactly what makes it straightforward to expose
//! from WebAssembly: build a [`PipeCaseConfig`] from JSON, call
//! [`run_pipe_case`], serialize the [`PipeCaseResult`] back to JSON.

use crate::boundary::BoundaryCondition;
use crate::gas::{GasProperties, PrimitiveState};
use crate::mesh::Mesh;
use crate::network::{Junction, PipeNetwork};
use crate::pipe::Pipe;
use crate::solver::run_to_time;
use serde::{Deserialize, Serialize};

/// One pipe's geometry, initial condition, and end boundary conditions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipeSpec {
    /// Pipe length, meters.
    pub length: f64,
    /// Pipe diameter, meters (constant along its length).
    pub diameter: f64,
    /// Target finite-volume cell size, meters — see [`Mesh::uniform`].
    pub target_cell_size: f64,
    /// Initial pressure, Pa, uniform along the whole pipe.
    pub initial_pressure: f64,
    /// Initial temperature, Kelvin, uniform along the whole pipe.
    pub initial_temperature_kelvin: f64,
    /// Initial velocity, m/s, uniform along the whole pipe.
    pub initial_velocity: f64,
    pub left_boundary: BoundaryCondition,
    pub right_boundary: BoundaryCondition,
}

/// A full pipe-network case: the working gas, every pipe's setup, how
/// pipes are joined, and how long/how carefully to run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipeCaseConfig {
    pub gas: GasProperties,
    /// Referenced by index (position in this list) from `junctions`.
    pub pipes: Vec<PipeSpec>,
    pub junctions: Vec<Junction>,
    /// How long to simulate, seconds.
    pub duration: f64,
    /// CFL number for the timestep (see [`crate::solver::step`]) — 1.0 is
    /// the stability limit for this scheme; use a smaller value (e.g. 0.5)
    /// for a comfortable safety margin.
    pub cfl: f64,
}

/// The final state of one pipe, sampled at every cell center.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipeResult {
    /// Position of each cell center along the pipe, meters. Same length
    /// and order as `pressure`/`velocity`/`temperature_kelvin`.
    pub cell_centers: Vec<f64>,
    pub pressure: Vec<f64>,
    pub velocity: Vec<f64>,
    pub temperature_kelvin: Vec<f64>,
}

/// The result of running a [`PipeCaseConfig`] to completion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipeCaseResult {
    /// Simulated time actually reached, seconds (may slightly exceed
    /// `PipeCaseConfig::duration` — see [`run_to_time`]).
    pub elapsed_time: f64,
    /// Same order as `PipeCaseConfig::pipes`.
    pub pipes: Vec<PipeResult>,
}

/// Builds the pipe network described by `config`, runs it for
/// `config.duration` seconds, and returns the final state of every pipe.
pub fn run_pipe_case(config: &PipeCaseConfig) -> PipeCaseResult {
    let pipes = config
        .pipes
        .iter()
        .map(|spec| {
            let mesh = Mesh::uniform(spec.length, spec.diameter, spec.target_cell_size);
            let initial_state = PrimitiveState::from_pressure_temperature(
                spec.initial_pressure,
                spec.initial_temperature_kelvin,
                spec.initial_velocity,
                &config.gas,
            );
            Pipe::uniform_initial_state(mesh, initial_state, &config.gas, spec.left_boundary, spec.right_boundary)
        })
        .collect();

    let mut network = PipeNetwork { pipes, junctions: config.junctions.clone() };
    let elapsed_time = run_to_time(&mut network, &config.gas, config.cfl, config.duration);

    let pipes = network
        .pipes
        .iter()
        .map(|pipe| {
            let cell_centers = pipe.mesh.cells.iter().map(|cell| cell.center).collect();
            let mut pressure = Vec::with_capacity(pipe.cell_count());
            let mut velocity = Vec::with_capacity(pipe.cell_count());
            let mut temperature_kelvin = Vec::with_capacity(pipe.cell_count());
            for state in &pipe.state {
                let primitive = state.to_primitive(&config.gas);
                pressure.push(primitive.pressure);
                velocity.push(primitive.velocity);
                temperature_kelvin.push(primitive.temperature_kelvin(&config.gas));
            }
            PipeResult { cell_centers, pressure, velocity, temperature_kelvin }
        })
        .collect();

    PipeCaseResult { elapsed_time, pipes }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_closed_pipe_at_rest_stays_at_rest() {
        let config = PipeCaseConfig {
            gas: GasProperties::AIR,
            pipes: vec![PipeSpec {
                length: 1.0,
                diameter: 0.05,
                target_cell_size: 0.01,
                initial_pressure: 150_000.0,
                initial_temperature_kelvin: 293.15,
                initial_velocity: 0.0,
                left_boundary: BoundaryCondition::ClosedEnd,
                right_boundary: BoundaryCondition::ClosedEnd,
            }],
            junctions: vec![],
            duration: 0.01,
            cfl: 0.5,
        };

        let result = run_pipe_case(&config);
        assert_eq!(result.pipes.len(), 1);
        for pressure in &result.pipes[0].pressure {
            assert!((pressure - 150_000.0).abs() < 1e-3);
        }
        for velocity in &result.pipes[0].velocity {
            assert!(velocity.abs() < 1e-6);
        }
    }

    #[test]
    fn config_round_trips_through_json() {
        let config = PipeCaseConfig {
            gas: GasProperties::AIR,
            pipes: vec![PipeSpec {
                length: 1.0,
                diameter: 0.05,
                target_cell_size: 0.01,
                initial_pressure: 101_325.0,
                initial_temperature_kelvin: 293.15,
                initial_velocity: 0.0,
                left_boundary: BoundaryCondition::ClosedEnd,
                right_boundary: BoundaryCondition::Reservoir { pressure: 101_325.0, temperature_kelvin: 293.15 },
            }],
            junctions: vec![],
            duration: 0.001,
            cfl: 0.5,
        };
        let json = serde_json::to_string(&config).expect("config should serialize");
        let recovered: PipeCaseConfig = serde_json::from_str(&json).expect("config should deserialize");
        assert_eq!(recovered.pipes.len(), 1);
        assert_eq!(recovered.pipes[0].length, 1.0);
    }
}
