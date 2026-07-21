//! Geometry of a 1D pipe discretized into finite-volume cells, including
//! optional taper (linearly varying diameter along the pipe's length).

use serde::{Deserialize, Serialize};

/// A single finite-volume cell: where it sits along the pipe axis, how wide
/// it is, and the pipe's cross-sectional area at the cell's center.
///
/// A tapered pipe's area varies continuously, so this is only the area at
/// the cell's own center — see [`Mesh::face_areas`] for the (generally
/// different) areas at the cell's two boundaries, which the quasi-1D flux
/// update needs separately.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Cell {
    /// Position of the cell's center along the pipe axis, meters, measured
    /// from the pipe's left (x=0) end.
    pub center: f64,
    /// Cell length along the pipe axis, meters. This is the "dx" in the
    /// finite-volume update.
    pub width: f64,
    /// Cross-sectional area at the cell's center, square meters.
    pub area: f64,
}

/// A 1D discretization of a pipe into a sequence of equal-width cells
/// running from one end to the other, with diameter linearly interpolated
/// along the pipe's length — the same modeling convention OpenWAM uses for
/// conical/tapered sections. A constant-diameter pipe is just the special
/// case where the two end diameters match, which makes every area-change
/// term in the solver identically zero rather than a special code path.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Mesh {
    pub cells: Vec<Cell>,
    /// Cross-sectional area at each cell boundary, square meters.
    /// `face_areas[i]` is the area at the left boundary of `cells[i]` (and
    /// the right boundary of `cells[i-1]`); length is always
    /// `cells.len() + 1`. Computed directly from the diameter at each
    /// face's own position — not by averaging neighboring cell-center
    /// areas, which would only coincidentally agree for a linear taper.
    pub face_areas: Vec<f64>,
}

impl Mesh {
    /// Builds a mesh for a pipe whose diameter varies linearly from
    /// `diameter_left` (at x=0) to `diameter_right` (at x=length).
    ///
    /// The number of cells is `length / target_cell_size`, rounded to the
    /// nearest whole number (minimum 1) — the same convention OpenWAM uses
    /// (e.g. a 1 m pipe with a 0.01 m target cell size produces a 100-cell
    /// mesh; the reference cases in `benchmarks/openwam/` report 101 cells
    /// for the same nominal target due to how OpenWAM counts mesh points
    /// rather than cells — this solver counts cells directly, so expect a
    /// harmless off-by-one relative to OpenWAM's reported cell counts).
    pub fn tapered(length: f64, diameter_left: f64, diameter_right: f64, target_cell_size: f64) -> Self {
        let cell_count = (length / target_cell_size).round().max(1.0) as usize;
        let width = length / cell_count as f64;

        let diameter_at = |x: f64| -> f64 {
            let fraction = x / length;
            diameter_left + (diameter_right - diameter_left) * fraction
        };

        let cells = (0..cell_count)
            .map(|i| {
                let center = (i as f64 + 0.5) * width;
                Cell { center, width, area: circle_area(diameter_at(center)) }
            })
            .collect();

        let face_areas = (0..=cell_count).map(|i| circle_area(diameter_at(i as f64 * width))).collect();

        Mesh { cells, face_areas }
    }

    /// Builds a uniform mesh for a constant-diameter pipe — the taper-free
    /// special case of [`Mesh::tapered`].
    pub fn uniform(length: f64, diameter: f64, target_cell_size: f64) -> Self {
        Self::tapered(length, diameter, diameter, target_cell_size)
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

    #[test]
    fn uniform_mesh_has_constant_face_areas_equal_to_cell_areas() {
        let mesh = Mesh::uniform(1.0, 0.05, 0.01);
        let expected_area = std::f64::consts::PI * 0.025 * 0.025;
        assert_eq!(mesh.face_areas.len(), mesh.cell_count() + 1);
        for &face_area in &mesh.face_areas {
            assert!((face_area - expected_area).abs() < 1e-12);
        }
    }

    #[test]
    fn tapered_mesh_face_areas_match_end_diameters_exactly() {
        let mesh = Mesh::tapered(1.0, 0.05, 0.10, 0.01);
        let expected_left_area = circle_area(0.05);
        let expected_right_area = circle_area(0.10);
        assert!((mesh.face_areas.first().unwrap() - expected_left_area).abs() < 1e-12);
        assert!((mesh.face_areas.last().unwrap() - expected_right_area).abs() < 1e-12);
    }

    #[test]
    fn tapered_mesh_face_areas_increase_monotonically_for_a_widening_pipe() {
        let mesh = Mesh::tapered(1.0, 0.05, 0.10, 0.01);
        for window in mesh.face_areas.windows(2) {
            assert!(window[1] > window[0], "face areas should strictly increase along a widening taper");
        }
    }
}
