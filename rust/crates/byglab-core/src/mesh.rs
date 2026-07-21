//! Geometry of a 1D pipe discretized into finite-volume cells.

use serde::{Deserialize, Serialize};

/// A single finite-volume cell: where it sits along the pipe axis, how wide
/// it is, and the pipe's cross-sectional area there.
///
/// Area is stored per-cell even though every current use case has a
/// constant-diameter pipe (area uniform along its whole length) — this
/// keeps the door open for tapered/stepped ducts (intake horns, diffusers)
/// later without having to change the `Cell` layout, only how it's filled
/// in and how the flux update uses it.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Cell {
    /// Position of the cell's center along the pipe axis, meters, measured
    /// from the pipe's left (x=0) end.
    pub center: f64,
    /// Cell length along the pipe axis, meters. This is the "dx" in the
    /// finite-volume update.
    pub width: f64,
    /// Cross-sectional area at this cell, square meters.
    pub area: f64,
}

/// A 1D discretization of a pipe into a sequence of cells running from one
/// end to the other. Currently always uniform (equal width, equal area)
/// since every pipe modeled so far has constant diameter and a single
/// mesh-size target, matching how OpenWAM's own pipes are meshed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Mesh {
    pub cells: Vec<Cell>,
}

impl Mesh {
    /// Builds a uniform mesh for a constant-diameter pipe.
    ///
    /// The number of cells is `length / target_cell_size`, rounded to the
    /// nearest whole number (minimum 1) — the same convention OpenWAM uses
    /// (e.g. a 1 m pipe with a 0.01 m target cell size produces a 100-cell
    /// mesh; the reference cases in `benchmarks/openwam/` report 101 cells
    /// for the same nominal target due to how OpenWAM counts mesh points
    /// rather than cells — this solver counts cells directly, so expect a
    /// harmless off-by-one relative to OpenWAM's reported cell counts).
    pub fn uniform(length: f64, diameter: f64, target_cell_size: f64) -> Self {
        let cell_count = (length / target_cell_size).round().max(1.0) as usize;
        let width = length / cell_count as f64;
        let area = circle_area(diameter);

        let cells = (0..cell_count)
            .map(|i| Cell { center: (i as f64 + 0.5) * width, width, area })
            .collect();

        Mesh { cells }
    }

    /// Number of cells in this mesh.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
}

/// Cross-sectional area of a circular pipe from its diameter.
fn circle_area(diameter: f64) -> f64 {
    std::f64::consts::PI * (diameter / 2.0).powi(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_mesh_has_expected_cell_count_and_total_length() {
        let mesh = Mesh::uniform(1.0, 0.05, 0.01);
        assert_eq!(mesh.cell_count(), 100);
        let total_length: f64 = mesh.cells.iter().map(|c| c.width).sum();
        assert!((total_length - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cell_centers_are_evenly_spaced_and_span_the_pipe() {
        let mesh = Mesh::uniform(1.0, 0.05, 0.01);
        assert!((mesh.cells.first().unwrap().center - 0.005).abs() < 1e-12);
        assert!((mesh.cells.last().unwrap().center - 0.995).abs() < 1e-12);
    }

    #[test]
    fn area_matches_circle_area_of_given_diameter() {
        let mesh = Mesh::uniform(1.0, 0.05, 0.01);
        let expected_area = std::f64::consts::PI * 0.025 * 0.025;
        for cell in &mesh.cells {
            assert!((cell.area - expected_area).abs() < 1e-12);
        }
    }
}
