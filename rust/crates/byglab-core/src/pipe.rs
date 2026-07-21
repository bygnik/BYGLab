//! A single 1D pipe: its mesh, the gas state in every cell, and how its two
//! ends are terminated.

use crate::boundary::BoundaryCondition;
use crate::gas::{ConservedState, Flux, GasProperties, PrimitiveState};
use crate::mesh::Mesh;
use crate::reconstruction::{reconstruct_cell, ReconstructedCell};
use crate::riemann::hllc_flux;
use crate::source_terms::{self, WallProperties};
use serde::{Deserialize, Serialize};

/// A pipe's gas state discretized over its [`Mesh`], with a
/// [`BoundaryCondition`] at each end.
///
/// A pipe's own boundary conditions only apply when that end is *not*
/// joined to another pipe by a [`crate::network::Junction`]. The
/// neighbor/reconstruction methods below are deliberately separate from
/// the assembled-flux convenience method — this is the seam
/// [`crate::network::PipeNetwork`] uses to substitute the joined pipe's
/// real boundary cell for either end's neighbor (and to build a properly
/// symmetric two-sided junction flux) without `Pipe` needing to know a
/// junction exists.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pipe {
    pub mesh: Mesh,
    /// One conserved gas state per cell, same length and order as `mesh.cells`.
    pub state: Vec<ConservedState>,
    pub left_boundary: BoundaryCondition,
    pub right_boundary: BoundaryCondition,
    /// Wall friction/heat-transfer properties, or `None` for a
    /// frictionless, adiabatic pipe (the default, and the only behavior
    /// that existed before `source_terms.rs` — every pre-existing
    /// validation case relies on this being `None`).
    pub wall: Option<WallProperties>,
}

impl Pipe {
    /// Builds a pipe filled with the same gas state everywhere — the
    /// common case for setting up a validation/initial-condition case.
    /// Frictionless and adiabatic by default; see [`Self::with_wall`].
    pub fn uniform_initial_state(
        mesh: Mesh,
        initial_state: PrimitiveState,
        gas: &GasProperties,
        left_boundary: BoundaryCondition,
        right_boundary: BoundaryCondition,
    ) -> Self {
        let state = vec![initial_state.to_conserved(gas); mesh.cell_count()];
        Pipe { mesh, state, left_boundary, right_boundary, wall: None }
    }

    /// Opts this pipe into wall friction and/or heat transfer.
    pub fn with_wall(mut self, wall: WallProperties) -> Self {
        self.wall = Some(wall);
        self
    }

    pub fn cell_count(&self) -> usize {
        self.mesh.cell_count()
    }

    /// The (unreconstructed) state to treat as this pipe's neighbor just
    /// past its left end, derived from its own boundary condition.
    /// [`crate::network::PipeNetwork`] substitutes the joined pipe's real
    /// boundary cell instead, for a junction-coupled end.
    pub fn own_left_neighbor(&self, gas: &GasProperties) -> PrimitiveState {
        self.left_boundary.ghost_state(self.state[0].to_primitive(gas), gas)
    }

    /// The (unreconstructed) state to treat as this pipe's neighbor just
    /// past its right end. See [`Self::own_left_neighbor`].
    pub fn own_right_neighbor(&self, gas: &GasProperties) -> PrimitiveState {
        self.right_boundary.ghost_state(self.state[self.state.len() - 1].to_primitive(gas), gas)
    }

    /// This pipe's actual leftmost cell state — used by
    /// [`crate::network::PipeNetwork`] as the neighbor another pipe should
    /// use for reconstruction across a junction at this end.
    pub fn left_boundary_cell_state(&self, gas: &GasProperties) -> PrimitiveState {
        self.state[0].to_primitive(gas)
    }

    /// This pipe's actual rightmost cell state. See [`Self::left_boundary_cell_state`].
    pub fn right_boundary_cell_state(&self, gas: &GasProperties) -> PrimitiveState {
        self.state[self.state.len() - 1].to_primitive(gas)
    }

    /// Slope-limits and half-timestep-evolves every cell (the
    /// MUSCL-Hancock predictor step, see `reconstruction.rs`), given the
    /// neighbor states to use just past each end. Used directly by
    /// [`crate::network::PipeNetwork`] (with the correct cross-pipe
    /// neighbor for junction-coupled ends) as well as by
    /// [`Self::own_face_fluxes`] (with this pipe's own boundary neighbors).
    pub fn reconstruct_cells(
        &self,
        gas: &GasProperties,
        dt: f64,
        left_neighbor: PrimitiveState,
        right_neighbor: PrimitiveState,
    ) -> Vec<ReconstructedCell> {
        let cell_count = self.cell_count();
        (0..cell_count)
            .map(|i| {
                let before = if i == 0 { left_neighbor } else { self.state[i - 1].to_primitive(gas) };
                let center = self.state[i].to_primitive(gas);
                let after = if i == cell_count - 1 { right_neighbor } else { self.state[i + 1].to_primitive(gas) };
                reconstruct_cell(before, center, after, self.mesh.cells[i].width, dt, gas)
            })
            .collect()
    }

    /// All face fluxes for this pipe in isolation (length `cell_count + 1`),
    /// using this pipe's own boundary conditions at both ends. Convenient
    /// for pipes with no junction-coupled end; a
    /// [`crate::network::PipeNetwork`] with junctions instead calls
    /// [`Self::reconstruct_cells`] directly with the correct cross-pipe
    /// neighbor and combines both sides' reconstructed edges into a
    /// properly symmetric junction flux.
    pub fn own_face_fluxes(&self, gas: &GasProperties, dt: f64) -> Vec<Flux> {
        let left_neighbor = self.own_left_neighbor(gas);
        let right_neighbor = self.own_right_neighbor(gas);
        Self::assemble_face_fluxes(&self.reconstruct_cells(gas, dt, left_neighbor, right_neighbor), left_neighbor, right_neighbor, gas)
    }

    /// Assembles a full list of face fluxes (length `cell_count + 1`) from
    /// a pipe's already-reconstructed cells and its two end neighbors —
    /// shared logic between [`Self::own_face_fluxes`] and
    /// [`crate::network::PipeNetwork`], which calls this with
    /// junction-aware neighbors before overriding the junction-coupled end
    /// face(s) with a two-sided flux.
    pub fn assemble_face_fluxes(
        reconstructed_cells: &[ReconstructedCell],
        left_neighbor: PrimitiveState,
        right_neighbor: PrimitiveState,
        gas: &GasProperties,
    ) -> Vec<Flux> {
        let mut fluxes = Vec::with_capacity(reconstructed_cells.len() + 1);
        fluxes.push(hllc_flux(left_neighbor, reconstructed_cells[0].left_face, gas));
        for pair in reconstructed_cells.windows(2) {
            fluxes.push(hllc_flux(pair[0].right_face, pair[1].left_face, gas));
        }
        fluxes.push(hllc_flux(reconstructed_cells[reconstructed_cells.len() - 1].right_face, right_neighbor, gas));
        fluxes
    }

    /// Advances every cell by one explicit timestep, given the flux at
    /// every face (length `cell_count + 1`, first and last being the two
    /// end faces) and each cell's reconstructed state (needed for the
    /// taper source term below).
    ///
    /// The predictor half-step is already baked into how the face fluxes
    /// themselves were computed (see `reconstruction.rs`), so the
    /// conservative update here is still a single full-timestep
    /// forward-Euler-style step — see `solver.rs`'s doc comment for why
    /// that's sufficient for full second-order accuracy with this specific
    /// scheme.
    ///
    /// For a constant-diameter, wall-less pipe this reduces exactly to the
    /// original `dU/dt = -(F_right - F_left)/dx` update (face areas equal
    /// the cell area everywhere, so area-weighting and volume-dividing
    /// cancel, and the geometric/wall source terms are both identically
    /// zero). For a tapered pipe, fluxes are weighted by their own face's
    /// area before differencing, a geometric pressure-force source term
    /// `p̄ * (A_right_face - A_left_face)` is added to momentum (`p̄` being
    /// the average of this cell's two MUSCL-reconstructed face pressures —
    /// the choice that makes a uniform, at-rest tapered pipe generate
    /// exactly zero spurious velocity, see `tests/tapered_pipe_stays_at_rest.rs`),
    /// and the whole update is divided by cell *volume* (`area * width`)
    /// rather than just `width`. Wall friction/heat-transfer source terms
    /// (`source_terms::wall_sources`) are added the same way when `self.wall`
    /// is set.
    pub fn apply_face_fluxes(
        &mut self,
        face_fluxes: &[Flux],
        reconstructed_cells: &[ReconstructedCell],
        dt: f64,
        gas: &GasProperties,
    ) {
        debug_assert_eq!(face_fluxes.len(), self.state.len() + 1);
        debug_assert_eq!(reconstructed_cells.len(), self.state.len());

        for i in 0..self.state.len() {
            let cell = self.mesh.cells[i];
            let left_face_area = self.mesh.face_areas[i];
            let right_face_area = self.mesh.face_areas[i + 1];

            let area_weighted_divergence = face_fluxes[i + 1] * right_face_area - face_fluxes[i] * left_face_area;

            let mean_face_pressure =
                0.5 * (reconstructed_cells[i].left_face.pressure + reconstructed_cells[i].right_face.pressure);
            let geometric_momentum_source = mean_face_pressure * (right_face_area - left_face_area);

            let (wall_momentum_source, wall_energy_source) = match &self.wall {
                Some(wall) => {
                    let diameter = 2.0 * (cell.area / std::f64::consts::PI).sqrt();
                    source_terms::wall_sources(self.state[i].to_primitive(gas), wall, diameter, gas)
                }
                None => (0.0, 0.0),
            };

            let cell_volume = cell.area * cell.width;
            let dt_over_volume = dt / cell_volume;

            // The flux-divergence and geometric source terms are both in
            // force/power units (flux * area, or pressure * area) - a
            // total rate-of-change of the cell's *content* (U * volume),
            // so both need dividing by volume to get dU/dt. Wall
            // friction/heat-transfer sources are already expressed *per
            // unit volume* (see `source_terms::wall_sources`'s doc
            // comment) - i.e. already a dU/dt contribution - so they're
            // scaled by `dt` alone, not `dt_over_volume` again.
            self.state[i].mass -= area_weighted_divergence.mass * dt_over_volume;
            self.state[i].momentum -= area_weighted_divergence.momentum * dt_over_volume;
            self.state[i].momentum += geometric_momentum_source * dt_over_volume;
            self.state[i].momentum += wall_momentum_source * dt;
            self.state[i].energy -= area_weighted_divergence.energy * dt_over_volume;
            self.state[i].energy += wall_energy_source * dt;
        }
    }

    /// The largest timestep this pipe alone could take under the CFL
    /// condition, i.e. `min` over its cells of `cfl * cell_width / (|u| + c)`.
    ///
    /// A network of multiple pipes must use the *minimum* of this value
    /// across every pipe, not each pipe's own value independently — see
    /// [`crate::solver::step`]'s doc comment for why a shared timestep is
    /// required once pipes are junction-coupled.
    pub fn cfl_time_step(&self, gas: &GasProperties, cfl: f64) -> f64 {
        self.state
            .iter()
            .zip(self.mesh.cells.iter())
            .map(|(state, cell)| {
                let primitive = state.to_primitive(gas);
                let wave_speed = primitive.velocity.abs() + primitive.sound_speed(gas);
                cfl * cell.width / wave_speed
            })
            .fold(f64::INFINITY, f64::min)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary::BoundaryCondition;
    use crate::mesh::Mesh;

    fn air_at_rest(pressure: f64) -> PrimitiveState {
        PrimitiveState::from_pressure_temperature(pressure, 293.15, 0.0, &GasProperties::AIR)
    }

    #[test]
    fn own_face_fluxes_has_one_more_entry_than_cells() {
        let gas = GasProperties::AIR;
        let mesh = Mesh::uniform(1.0, 0.05, 0.01);
        let pipe = Pipe::uniform_initial_state(
            mesh,
            air_at_rest(150_000.0),
            &gas,
            BoundaryCondition::ClosedEnd,
            BoundaryCondition::ClosedEnd,
        );
        let fluxes = pipe.own_face_fluxes(&gas, 1e-6);
        assert_eq!(fluxes.len(), pipe.cell_count() + 1);
    }

    #[test]
    fn closed_closed_pipe_at_rest_has_zero_net_flux_everywhere() {
        // A uniform pipe at rest, closed both ends: every face sees
        // matching states either side after reconstruction (uniform state
        // has zero slope everywhere), giving zero mass flux throughout —
        // the mechanistic reason `tests/quiescent.rs` stays at rest exactly.
        let gas = GasProperties::AIR;
        let mesh = Mesh::uniform(1.0, 0.05, 0.01);
        let pipe = Pipe::uniform_initial_state(
            mesh,
            air_at_rest(150_000.0),
            &gas,
            BoundaryCondition::ClosedEnd,
            BoundaryCondition::ClosedEnd,
        );
        let fluxes = pipe.own_face_fluxes(&gas, 1e-6);
        for flux in &fluxes {
            assert!(flux.mass.abs() < 1e-9, "expected zero mass flux, got {}", flux.mass);
        }
    }
}
