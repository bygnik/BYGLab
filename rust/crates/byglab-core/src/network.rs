//! Multiple pipes connected by junctions.
//!
//! A same-area, same-gas junction between two pipe ends is treated as an
//! ordinary *internal* face of one continuous mesh: a single HLLC flux is
//! computed between the two pipes' adjoining boundary cells (or, with the
//! MUSCL-Hancock upgrade, their properly reconstructed boundary *faces*)
//! and applied to both sides. This is mathematically identical to an
//! internal face inside a single pipe (a numerical flux only depends on
//! the two adjacent states, not on which array they live in —
//! conservation just requires applying the same flux with opposite sign to
//! both neighbors), so no iterative junction solver is needed when area
//! and gas properties match on both sides.
//!
//! Getting second-order accuracy right AT a junction means each pipe's
//! reconstruction at its joined end must use the *other* pipe's real
//! boundary cell as its neighbor (not a boundary-condition-derived ghost),
//! and the shared flux must be built from *both* pipes' reconstructed
//! edges — using only one side's reconstruction (treating the other pipe's
//! cell as an unreconstructed ghost) would silently drop back to
//! first-order accuracy right at the junction. [`PipeNetwork::all_face_fluxes`]
//! does this in two passes: reconstruct every pipe with the correct
//! neighbors first, then resolve junction fluxes from both sides' results.

use crate::gas::{Flux, GasProperties, PrimitiveState};
use crate::pipe::Pipe;
use crate::reconstruction::ReconstructedCell;
use crate::riemann::hllc_flux;
use serde::{Deserialize, Serialize};

/// Which end of a pipe a [`PipeEndRef`] refers to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PipeEnd {
    Left,
    Right,
}

/// Identifies one specific end of one specific pipe within a
/// [`PipeNetwork`]'s `pipes` list.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct PipeEndRef {
    pub pipe_index: usize,
    pub end: PipeEnd,
}

/// A direct, same-area, same-gas connection between two pipe ends.
///
/// Only connects a pipe's `Right` end to another pipe's `Left` end (in
/// either order) — physically, gluing pipe A's right end to pipe B's left
/// end forms one continuous coordinate line where positions increase
/// monotonically through both pipes and "positive velocity" means the same
/// thing on both sides, so no sign convention needs to change. A `Left`-to-
/// `Left` or `Right`-to-`Right` junction would need one pipe's velocity
/// sign flipped to make physical sense (its local +x axis points the
/// "wrong way" relative to the joined pipe) — not needed by any current
/// use case, and not implemented; [`PipeNetwork`] panics with a clear
/// message if one is constructed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Junction {
    pub a: PipeEndRef,
    pub b: PipeEndRef,
}

/// A collection of pipes plus the junctions connecting them.
///
/// Owns the network-wide view needed for correct multi-pipe stepping:
/// junction-aware, second-order-accurate face fluxes (see
/// [`Self::all_face_fluxes`]) and a timestep that respects every cell in
/// every pipe (see [`Self::cfl_time_step`]) — see [`crate::solver::step`]
/// for why a single shared timestep across the whole network is required
/// once pipes are junction-coupled.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipeNetwork {
    pub pipes: Vec<Pipe>,
    pub junctions: Vec<Junction>,
}

impl PipeNetwork {
    /// A network containing just one pipe and no junctions.
    pub fn single_pipe(pipe: Pipe) -> Self {
        PipeNetwork { pipes: vec![pipe], junctions: Vec::new() }
    }

    /// The largest timestep the whole network could take under the CFL
    /// condition — the minimum of every individual pipe's
    /// [`Pipe::cfl_time_step`]. Always use this (never a single pipe's own
    /// value) once more than one pipe is involved.
    pub fn cfl_time_step(&self, gas: &GasProperties, cfl: f64) -> f64 {
        self.pipes.iter().map(|pipe| pipe.cfl_time_step(gas, cfl)).fold(f64::INFINITY, f64::min)
    }

    /// Advances every pipe in the network by one explicit timestep of size
    /// `dt`, resolving junction-coupled end faces with a shared,
    /// second-order-accurate flux instead of each pipe's own boundary
    /// condition.
    pub fn advance(&mut self, dt: f64, gas: &GasProperties) {
        let all_fluxes = self.all_face_fluxes(dt, gas);
        for (pipe, fluxes) in self.pipes.iter_mut().zip(all_fluxes.iter()) {
            pipe.apply_face_fluxes(fluxes, dt);
        }
    }

    /// Every pipe's full list of face fluxes (length `cell_count + 1`).
    ///
    /// Two passes: first, determine the correct neighbor state just past
    /// each pipe end (a pipe's own boundary condition by default, or the
    /// joined pipe's real boundary cell for a junction-coupled end) and
    /// reconstruct every pipe using those neighbors; then assemble each
    /// pipe's own fluxes, and finally overwrite junction-coupled end faces
    /// with a flux built from *both* connected pipes' reconstructed edges
    /// (see this module's doc comment for why that second step matters).
    fn all_face_fluxes(&self, dt: f64, gas: &GasProperties) -> Vec<Vec<Flux>> {
        let mut left_neighbors: Vec<PrimitiveState> =
            self.pipes.iter().map(|pipe| pipe.own_left_neighbor(gas)).collect();
        let mut right_neighbors: Vec<PrimitiveState> =
            self.pipes.iter().map(|pipe| pipe.own_right_neighbor(gas)).collect();

        for junction in &self.junctions {
            let (a, b) = self.validate_junction_orientation(junction);
            *Self::neighbor_slot(&mut left_neighbors, &mut right_neighbors, a) = self.end_state(b, gas);
            *Self::neighbor_slot(&mut left_neighbors, &mut right_neighbors, b) = self.end_state(a, gas);
        }

        let reconstructed: Vec<Vec<ReconstructedCell>> = self
            .pipes
            .iter()
            .enumerate()
            .map(|(i, pipe)| pipe.reconstruct_cells(gas, dt, left_neighbors[i], right_neighbors[i]))
            .collect();

        let mut fluxes: Vec<Vec<Flux>> = (0..self.pipes.len())
            .map(|i| Pipe::assemble_face_fluxes(&reconstructed[i], left_neighbors[i], right_neighbors[i], gas))
            .collect();

        for junction in &self.junctions {
            let (a, b) = self.validate_junction_orientation(junction);
            // a is the Right end, b is the Left end (validated above) - the
            // shared flux uses a's pipe's reconstructed RIGHT face and b's
            // pipe's reconstructed LEFT face, exactly like an ordinary
            // internal face between two adjacent cells.
            let a_edge = Self::edge_face(&reconstructed[a.pipe_index], a.end);
            let b_edge = Self::edge_face(&reconstructed[b.pipe_index], b.end);
            let shared_flux = hllc_flux(a_edge, b_edge, gas);
            Self::set_end_flux(&mut fluxes, a, shared_flux);
            Self::set_end_flux(&mut fluxes, b, shared_flux);
        }

        fluxes
    }

    /// Checks that `junction` connects a `Right` end to a `Left` end (in
    /// either order) and returns `(right_end, left_end)` — panics with a
    /// clear message otherwise (see [`Junction`]'s doc comment for why
    /// only this orientation is supported).
    fn validate_junction_orientation(&self, junction: &Junction) -> (PipeEndRef, PipeEndRef) {
        match (junction.a.end, junction.b.end) {
            (PipeEnd::Right, PipeEnd::Left) => (junction.a, junction.b),
            (PipeEnd::Left, PipeEnd::Right) => (junction.b, junction.a),
            _ => panic!(
                "Junction {{a: {:?}, b: {:?}}} connects two `{:?}` ends - only a `Right`-to-`Left` \
                 connection (in either order) is supported, since that's the only case where both \
                 pipes' velocity sign conventions already agree; see `Junction`'s doc comment.",
                junction.a, junction.b, junction.a.end
            ),
        }
    }

    /// The primitive gas state in the boundary cell at one specific pipe end.
    fn end_state(&self, end_ref: PipeEndRef, gas: &GasProperties) -> PrimitiveState {
        let pipe = &self.pipes[end_ref.pipe_index];
        match end_ref.end {
            PipeEnd::Left => pipe.left_boundary_cell_state(gas),
            PipeEnd::Right => pipe.right_boundary_cell_state(gas),
        }
    }

    /// The reconstructed face value at one end of a pipe's cell list — the
    /// face that actually touches a junction (a pipe's `Right` end touches
    /// its last cell's `right_face`; its `Left` end touches its first
    /// cell's `left_face`).
    fn edge_face(reconstructed: &[ReconstructedCell], end: PipeEnd) -> PrimitiveState {
        match end {
            PipeEnd::Left => reconstructed[0].left_face,
            PipeEnd::Right => reconstructed[reconstructed.len() - 1].right_face,
        }
    }

    /// Mutable reference into whichever of `left_neighbors`/`right_neighbors`
    /// corresponds to `end_ref`.
    fn neighbor_slot<'a>(
        left_neighbors: &'a mut [PrimitiveState],
        right_neighbors: &'a mut [PrimitiveState],
        end_ref: PipeEndRef,
    ) -> &'a mut PrimitiveState {
        match end_ref.end {
            PipeEnd::Left => &mut left_neighbors[end_ref.pipe_index],
            PipeEnd::Right => &mut right_neighbors[end_ref.pipe_index],
        }
    }

    /// Overwrites the face-flux entry for one specific pipe end within the
    /// per-pipe flux lists built by [`Self::all_face_fluxes`].
    fn set_end_flux(fluxes: &mut [Vec<Flux>], end_ref: PipeEndRef, flux: Flux) {
        let pipe_fluxes = &mut fluxes[end_ref.pipe_index];
        match end_ref.end {
            PipeEnd::Left => pipe_fluxes[0] = flux,
            PipeEnd::Right => {
                let last = pipe_fluxes.len() - 1;
                pipe_fluxes[last] = flux;
            }
        }
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

    fn two_pipe_network(pressure_a: f64, pressure_b: f64) -> PipeNetwork {
        let gas = GasProperties::AIR;
        let pipe_a = Pipe::uniform_initial_state(
            Mesh::uniform(1.0, 0.05, 0.01),
            air_at_rest(pressure_a),
            &gas,
            BoundaryCondition::ClosedEnd,
            BoundaryCondition::ClosedEnd,
        );
        let pipe_b = Pipe::uniform_initial_state(
            Mesh::uniform(1.0, 0.05, 0.01),
            air_at_rest(pressure_b),
            &gas,
            BoundaryCondition::ClosedEnd,
            BoundaryCondition::ClosedEnd,
        );
        PipeNetwork {
            pipes: vec![pipe_a, pipe_b],
            junctions: vec![Junction {
                a: PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
                b: PipeEndRef { pipe_index: 1, end: PipeEnd::Left },
            }],
        }
    }

    #[test]
    fn matched_states_across_a_junction_produce_zero_net_flux() {
        let gas = GasProperties::AIR;
        let network = two_pipe_network(150_000.0, 150_000.0);
        let fluxes = network.all_face_fluxes(1e-6, &gas);
        let junction_flux = fluxes[0][fluxes[0].len() - 1];
        assert!(junction_flux.mass.abs() < 1e-9);
    }

    #[test]
    #[should_panic(expected = "only a `Right`-to-`Left` connection")]
    fn left_to_left_junction_panics_with_a_clear_message() {
        let gas = GasProperties::AIR;
        let mut network = two_pipe_network(150_000.0, 150_000.0);
        network.junctions[0] = Junction {
            a: PipeEndRef { pipe_index: 0, end: PipeEnd::Left },
            b: PipeEndRef { pipe_index: 1, end: PipeEnd::Left },
        };
        network.all_face_fluxes(1e-6, &gas);
    }
}
