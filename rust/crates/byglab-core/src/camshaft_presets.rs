//! Real BMW S54 camshaft grind data, ingested directly from official
//! Schrick camshaft data sheets ("Nockenwellen Datenblatt") — not
//! estimated or approximated. Each constant below cites its exact source
//! PDF (schrick.com/media/.../<part-number>.pdf, all revised 2013-09-04
//! through 2013-09-11).
//!
//! **Crank angle convention**: these data sheets reference "0°" to the
//! *gas-exchange* TDC (the TDC between the exhaust and intake strokes),
//! not the firing/combustion TDC that `combustion.rs`'s Wiebe timing uses
//! — confirmed by cross-checking the exhaust lobe centerlines (e.g.
//! `0415A1960-01`'s -104° "Auslass Spreizung" matches exactly
//! `(open + close) / 2 = (-252 + 44) / 2`) against standard EVO/EVC
//! timing (72° BBDC / 44° ATDC relative to the *gas-exchange* TDC, both
//! physically ordinary figures — they would not be if "0" meant firing
//! TDC instead). This matters if this data is ever combined with the
//! combustion model in the same simulation (a ±360° shift would be
//! needed to align the two conventions) — not needed for a standalone
//! breathing/valve-flow simulation, since [`crate::cylinder::Cylinder::volume`]
//! is exactly 360°-periodic and doesn't care which TDC "0" refers to.
//!
//! **What's ingested vs. what isn't**: duration, valve lift, and installed
//! open/close timing (both spread options where the sheet publishes two)
//! are used directly to build a [`CamProfile`]. The manufacturer's own
//! `@1mm` figures are *not* used to build the profile — they're kept here
//! purely as an honest cross-check of how well this crate's idealized
//! versine shape (`camshaft.rs`) matches the real (more sophisticated,
//! ramped) cam lobe; see `camshaft_presets` tests for the measured
//! discrepancy, which is real and non-trivial, not hidden.

use crate::camshaft::CamProfile;

/// One published install option for a [`SchrickGrind`]: the lobe
/// separation angle and the resulting opening/closing crank angles (all
/// in Schrick's own "0 = gas-exchange TDC" convention, degrees), plus the
/// manufacturer's own `@1mm` cross-check figures.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SchrickSpread {
    pub lobe_separation_degrees: f64,
    pub opening_angle_degrees: f64,
    pub closing_angle_degrees: f64,
    /// Manufacturer-published opening angle at 1mm lift — not used to
    /// build the [`CamProfile`], only for the versine-approximation
    /// cross-check.
    pub opens_at_1mm_degrees: f64,
    pub closes_at_1mm_degrees: f64,
}

/// One real Schrick camshaft grind's published specifications.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SchrickGrind {
    pub part_number: &'static str,
    pub application: &'static str,
    /// True (zero-lift, seat-to-seat) duration, degrees.
    pub duration_degrees: f64,
    pub valve_lift_mm: f64,
    pub spreads: &'static [SchrickSpread],
}

impl SchrickGrind {
    /// Builds a [`CamProfile`] for the given published install option
    /// (index into `self.spreads`) — `opening_angle_radians` here is
    /// still in Schrick's "0 = gas-exchange TDC" convention (see this
    /// module's doc comment).
    pub fn cam_profile(&self, spread_index: usize) -> CamProfile {
        let spread = &self.spreads[spread_index];
        CamProfile {
            max_lift: self.valve_lift_mm / 1000.0,
            opening_angle_radians: spread.opening_angle_degrees.to_radians(),
            duration_radians: self.duration_degrees.to_radians(),
        }
    }
}

/// Intake, "low street". Source: schrick.com/media/11/7e/1c/1653300739/0415E1800-00.pdf
pub const SCHRICK_0415E1800_00: SchrickGrind = SchrickGrind {
    part_number: "0415E1800-00",
    application: "low street (intake)",
    duration_degrees: 280.0,
    valve_lift_mm: 12.25,
    spreads: &[
        SchrickSpread { lobe_separation_degrees: 132.0, opening_angle_degrees: -8.0, closing_angle_degrees: 272.0, opens_at_1mm_degrees: 7.0, closes_at_1mm_degrees: 256.0 },
        SchrickSpread { lobe_separation_degrees: 72.0, opening_angle_degrees: -68.0, closing_angle_degrees: 212.0, opens_at_1mm_degrees: -53.0, closes_at_1mm_degrees: 196.0 },
    ],
};

/// Intake, "low street" (higher-duration variant sharing the same
/// timing/spread family as [`SCHRICK_0415E1880_01`] but at 12.25mm
/// instead of 13.75mm lift). Source: schrick.com/media/48/ed/b4/1653300774/0415E1880-00.pdf
pub const SCHRICK_0415E1880_00: SchrickGrind = SchrickGrind {
    part_number: "0415E1880-00",
    application: "low street (intake)",
    duration_degrees: 288.0,
    valve_lift_mm: 12.25,
    spreads: &[
        SchrickSpread { lobe_separation_degrees: 132.0, opening_angle_degrees: -12.0, closing_angle_degrees: 276.0, opens_at_1mm_degrees: 3.0, closes_at_1mm_degrees: 262.0 },
        SchrickSpread { lobe_separation_degrees: 72.0, opening_angle_degrees: -72.0, closing_angle_degrees: 216.0, opens_at_1mm_degrees: -57.0, closes_at_1mm_degrees: 202.0 },
    ],
};

/// Intake, "medium street". Source: schrick.com/media/de/21/44/1653300680/0415E1040-00.pdf
pub const SCHRICK_0415E1040_00: SchrickGrind = SchrickGrind {
    part_number: "0415E1040-00",
    application: "medium street (intake)",
    duration_degrees: 304.0,
    valve_lift_mm: 12.25,
    spreads: &[
        SchrickSpread { lobe_separation_degrees: 132.0, opening_angle_degrees: -20.0, closing_angle_degrees: 284.0, opens_at_1mm_degrees: -2.0, closes_at_1mm_degrees: 266.0 },
        SchrickSpread { lobe_separation_degrees: 72.0, opening_angle_degrees: -80.0, closing_angle_degrees: 224.0, opens_at_1mm_degrees: -62.0, closes_at_1mm_degrees: 206.0 },
    ],
};

/// Intake, "high street", single fixed (non-adjustable-spread) install.
/// Source: schrick.com/media/89/b3/b7/1653300709/0415E1040-01.pdf
pub const SCHRICK_0415E1040_01: SchrickGrind = SchrickGrind {
    part_number: "0415E1040-01",
    application: "high street (intake)",
    duration_degrees: 304.0,
    valve_lift_mm: 12.25,
    spreads: &[SchrickSpread { lobe_separation_degrees: 104.0, opening_angle_degrees: -48.0, closing_angle_degrees: 256.0, opens_at_1mm_degrees: -30.0, closes_at_1mm_degrees: 238.0 }],
};

/// Intake, "high street", highest-lift grind. Source: schrick.com/media/f8/b0/b8/1653300821/0415E1880-01.pdf
pub const SCHRICK_0415E1880_01: SchrickGrind = SchrickGrind {
    part_number: "0415E1880-01",
    application: "high street (intake)",
    duration_degrees: 288.0,
    valve_lift_mm: 13.75,
    spreads: &[
        SchrickSpread { lobe_separation_degrees: 132.0, opening_angle_degrees: -12.0, closing_angle_degrees: 276.0, opens_at_1mm_degrees: 1.0, closes_at_1mm_degrees: 263.0 },
        SchrickSpread { lobe_separation_degrees: 72.0, opening_angle_degrees: -72.0, closing_angle_degrees: 216.0, opens_at_1mm_degrees: -59.0, closes_at_1mm_degrees: 203.0 },
    ],
};

/// Exhaust, "low street". Source: schrick.com/media/b9/ab/74/1653300553/0415A1800-00.pdf
pub const SCHRICK_0415A1800_00: SchrickGrind = SchrickGrind {
    part_number: "0415A1800-00",
    application: "low street (exhaust)",
    duration_degrees: 280.0,
    valve_lift_mm: 12.25,
    spreads: &[
        SchrickSpread { lobe_separation_degrees: -130.0, opening_angle_degrees: -270.0, closing_angle_degrees: 10.0, opens_at_1mm_degrees: -255.0, closes_at_1mm_degrees: -6.0 },
        SchrickSpread { lobe_separation_degrees: -85.0, opening_angle_degrees: -225.0, closing_angle_degrees: 55.0, opens_at_1mm_degrees: -210.0, closes_at_1mm_degrees: 39.0 },
    ],
};

/// Exhaust, "low street" (shorter-duration variant). Source: schrick.com/media/92/1e/dc/1653298433/0415A1720-00.pdf
/// (the sheet's own internal table header reads "0415 E1 720-00" — a
/// data-entry typo on Schrick's part, contradicted by the same sheet's
/// "Auslass"/"Exhaust profile" labels and part-number footer/URL; treated
/// as the exhaust profile it clearly is throughout).
pub const SCHRICK_0415A1720_00: SchrickGrind = SchrickGrind {
    part_number: "0415A1720-00",
    application: "low street (exhaust)",
    duration_degrees: 272.0,
    valve_lift_mm: 12.25,
    spreads: &[
        SchrickSpread { lobe_separation_degrees: -128.0, opening_angle_degrees: -264.0, closing_angle_degrees: 8.0, opens_at_1mm_degrees: -249.0, closes_at_1mm_degrees: -8.0 },
        SchrickSpread { lobe_separation_degrees: -83.0, opening_angle_degrees: -219.0, closing_angle_degrees: 53.0, opens_at_1mm_degrees: -204.0, closes_at_1mm_degrees: 37.0 },
    ],
};

/// Exhaust, "medium street". Source: schrick.com/media/e9/c6/c8/1653300620/0415A1960-00.pdf
pub const SCHRICK_0415A1960_00: SchrickGrind = SchrickGrind {
    part_number: "0415A1960-00",
    application: "medium street (exhaust)",
    duration_degrees: 296.0,
    valve_lift_mm: 12.25,
    spreads: &[
        SchrickSpread { lobe_separation_degrees: -130.0, opening_angle_degrees: -278.0, closing_angle_degrees: 18.0, opens_at_1mm_degrees: -260.0, closes_at_1mm_degrees: 0.0 },
        SchrickSpread { lobe_separation_degrees: -85.0, opening_angle_degrees: -233.0, closing_angle_degrees: 63.0, opens_at_1mm_degrees: -215.0, closes_at_1mm_degrees: 45.0 },
    ],
};

/// Exhaust, "high street". Source: schrick.com/media/be/04/b1/1653300585/0415A1800-01.pdf
pub const SCHRICK_0415A1800_01: SchrickGrind = SchrickGrind {
    part_number: "0415A1800-01",
    application: "high street (exhaust)",
    duration_degrees: 280.0,
    valve_lift_mm: 13.75,
    spreads: &[
        SchrickSpread { lobe_separation_degrees: -130.0, opening_angle_degrees: -270.0, closing_angle_degrees: 10.0, opens_at_1mm_degrees: -257.0, closes_at_1mm_degrees: -4.0 },
        SchrickSpread { lobe_separation_degrees: -85.0, opening_angle_degrees: -225.0, closing_angle_degrees: 55.0, opens_at_1mm_degrees: -212.0, closes_at_1mm_degrees: 41.0 },
    ],
};

/// Exhaust, "high street", highest-lift grind, single fixed-spread
/// install. Source: schrick.com/media/d9/ea/08/1653300653/0415A1960-01.pdf
pub const SCHRICK_0415A1960_01: SchrickGrind = SchrickGrind {
    part_number: "0415A1960-01",
    application: "high street (exhaust)",
    duration_degrees: 296.0,
    valve_lift_mm: 13.75,
    spreads: &[SchrickSpread { lobe_separation_degrees: -104.0, opening_angle_degrees: -252.0, closing_angle_degrees: 44.0, opens_at_1mm_degrees: -234.0, closes_at_1mm_degrees: 26.0 }],
};

/// All ten ingested grinds, for iterating in tests/tooling.
pub const ALL_SCHRICK_S54_GRINDS: &[SchrickGrind] = &[
    SCHRICK_0415E1800_00,
    SCHRICK_0415E1880_00,
    SCHRICK_0415E1040_00,
    SCHRICK_0415E1040_01,
    SCHRICK_0415E1880_01,
    SCHRICK_0415A1800_00,
    SCHRICK_0415A1720_00,
    SCHRICK_0415A1960_00,
    SCHRICK_0415A1800_01,
    SCHRICK_0415A1960_01,
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::camshaft::lift_at;

    #[test]
    fn every_grinds_duration_matches_its_own_open_close_angles() {
        // Self-consistency check: duration = close - open, exactly, for
        // every published spread option on every grind - catches a
        // transcription error in this module directly (independent of
        // trusting the source PDFs' own internal consistency, which was
        // already confirmed by hand while ingesting them).
        for grind in ALL_SCHRICK_S54_GRINDS {
            for spread in grind.spreads {
                let implied_duration = spread.closing_angle_degrees - spread.opening_angle_degrees;
                assert!(
                    (implied_duration - grind.duration_degrees).abs() < 1e-9,
                    "{}: spread (open={}, close={}) implies duration {implied_duration}, but grind duration is {}",
                    grind.part_number,
                    spread.opening_angle_degrees,
                    spread.closing_angle_degrees,
                    grind.duration_degrees
                );
            }
        }
    }

    #[test]
    fn every_grinds_lobe_separation_matches_its_own_open_close_midpoint() {
        // Self-consistency check: lobe separation = (open+close)/2,
        // exactly, matching how these figures are actually defined.
        for grind in ALL_SCHRICK_S54_GRINDS {
            for spread in grind.spreads {
                let implied_lsa = 0.5 * (spread.opening_angle_degrees + spread.closing_angle_degrees);
                assert!(
                    (implied_lsa - spread.lobe_separation_degrees).abs() < 1e-9,
                    "{}: spread (open={}, close={}) implies LSA {implied_lsa}, but published LSA is {}",
                    grind.part_number,
                    spread.opening_angle_degrees,
                    spread.closing_angle_degrees,
                    spread.lobe_separation_degrees
                );
            }
        }
    }

    #[test]
    fn versine_approximation_vs_published_at_1mm_timing_across_all_grinds() {
        // The versine model is fit to match the TRUE (zero-lift) duration
        // and opening angle exactly, by construction - this checks how
        // well that idealized shape reproduces each grind's manufacturer-
        // published "@1mm" timing, an honest, real comparison against a
        // genuine (more sophisticated, ramped) cam lobe, not assumed to
        // match closely.
        let lift_threshold = 0.001; // 1mm

        let mut max_discrepancy_degrees = 0.0_f64;
        for grind in ALL_SCHRICK_S54_GRINDS {
            for (spread_index, spread) in grind.spreads.iter().enumerate() {
                let profile = grind.cam_profile(spread_index);

                // Bisect on the model's own (already independently
                // validated in camshaft.rs's tests) forward lift_at
                // function to find where it crosses 1mm, on the rising side.
                let mut low = profile.opening_angle_radians;
                let mut high = profile.opening_angle_radians + 0.5 * profile.duration_radians;
                for _ in 0..100 {
                    let mid = 0.5 * (low + high);
                    if lift_at(&profile, mid) < lift_threshold {
                        low = mid;
                    } else {
                        high = mid;
                    }
                }
                let model_opens_at_1mm_degrees = (0.5 * (low + high)).to_degrees();
                let discrepancy_degrees = model_opens_at_1mm_degrees - spread.opens_at_1mm_degrees;
                max_discrepancy_degrees = max_discrepancy_degrees.max(discrepancy_degrees.abs());

                println!(
                    "{} spread#{spread_index} (LSA={}): model opens@1mm={model_opens_at_1mm_degrees:.2} deg, published={} deg, discrepancy={:.2} deg",
                    grind.part_number, spread.lobe_separation_degrees, spread.opens_at_1mm_degrees, discrepancy_degrees
                );
            }
        }

        println!("max |discrepancy| across all ingested grinds/spreads: {max_discrepancy_degrees:.2} deg");

        // Document, don't hide, the real deviation of the idealized
        // versine shape from actual (ramp-profiled) cam lobes - a
        // double-digit-degree discrepancy is expected here and
        // informative, not a bug: it's exactly why direct-lift-table
        // support (not yet implemented) would be needed for high-
        // fidelity valve-event timing. The versine is a first
        // approximation for overall duration/lift, not lobe shape.
        assert!(max_discrepancy_degrees < 25.0, "discrepancy {max_discrepancy_degrees:.2} deg is larger than expected even for an approximate versine fit");
    }
}
