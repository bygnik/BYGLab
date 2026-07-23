//! An N-way branch junction: several pipe ends meeting at one point,
//! generalizing [`crate::network::Junction`] (exactly 2 pipes, matching
//! area) to any number of pipes with any areas — the architecture needed
//! for a single intake/exhaust runner splitting into several (or several
//! merging into one).
//!
//! # The physics (re-derived from real ground truth, not invented)
//!
//! A real branch junction has no throat and no restriction, so it is NOT
//! the same physics as [`crate::valve_port`]'s compressible-orifice
//! coupling (an earlier design draft of this module made exactly that
//! mistake, and a design review caught it before any code shipped — there
//! is no choked-flow regime here at all, since there's nothing to choke
//! against). The correct closure, independently re-derived from OpenWAM's
//! actual N-way branch boundary (`TCCRamificacion.cpp`, used in this
//! repo's real 6-cylinder 6-2-1 exhaust header case) and Gordon Blair's
//! "Design and Simulation of Four-Stroke Engines" (Sec. 2.13, "Benson's
//! constant pressure model"), is a method-of-characteristics closure:
//!
//! - All connected legs are assumed to share one instantaneous static
//!   pressure `P_j` at the junction (Benson's core assumption).
//! - **Every leg — supplying or receiving — is resolved by its own
//!   Riemann-invariant/isentropic-expansion relation, using its own
//!   entropy, unconditionally.** A `Right` end's known (interior-sourced)
//!   characteristic is `J = u + 2c/(gamma-1)`; a `Left` end's is
//!   `J = u - 2c/(gamma-1)` (sign convention matches [`crate::valve_port`]'s
//!   `end_sign`: `Right = +1`, `Left = -1`). Given a trial `P_j`, boundary
//!   sound speed follows `c(P_j) = c_own * (P_j/P_own)^((gamma-1)/(2*gamma))`
//!   and boundary velocity follows directly from `J`.
//! - **Two earlier, more elaborate designs were tried and independently
//!   disproved before landing on this one** (both are worth recording,
//!   since the failure modes are non-obvious): (1) swapping a *different*
//!   leg's sound speed into a receiving leg's own `J` equation — this is
//!   provably wrong: with every leg starting at the same reference
//!   temperature (every validation case here), `J` for a `Right`-end leg
//!   and a `Left`-end leg with the same state are exact negatives of each
//!   other, which forces the naively-computed receiving velocity to have
//!   the *same sign* as the supplying leg's, so mass can never balance
//!   except degenerately. (2) resolving receiving legs via isentropic
//!   acceleration from a mass-flux-weighted mixed *stagnation* state
//!   (mirroring Gordon Blair's Eq. 2.14.12 energy mixing) — this is also
//!   wrong, for a different, purely geometric reason: accelerating from one
//!   shared stagnation state to the same static pressure gives every
//!   receiving leg the *same* velocity magnitude regardless of its own
//!   area, so total mass conservation collapses to requiring the supplying
//!   leg's area to exactly equal the sum of the receiving legs' areas — a
//!   coincidence, not a general solution. Both were caught by direct
//!   numerical scans of `total_mass_flow(P_j)` over the bracket, which
//!   either never crossed zero in the interior or crossed it only at a
//!   geometry-dependent accident. The much simpler "own entropy,
//!   unconditionally" closure above passed that same scan cleanly
//!   (monotonic, single interior root) for every case tried, including
//!   unequal areas and more than 2 legs — matching how Gordon Blair's own
//!   *simple* constant-pressure model (Sec. 2.13) doesn't mix entropy at
//!   all for the pressure/velocity solve either; entropy mixing only
//!   matters for tracking what composition ends up resident inside a
//!   receiving pipe afterward (a longer-timescale bookkeeping concern, out
//!   of scope for a single-step junction resolve), not for this instant's
//!   force/mass balance.
//!
//! `P_j` itself is found by a 1-D root-find (see [`solve_junction_pressure`])
//! enforcing exact mass conservation across the whole junction, since that
//! is literally the equation being solved.
//!
//! # A real, honest limitation (measured and documented, not hidden)
//!
//! This closure gives **exact mass conservation** (the equation solved
//! for) but not exact energy conservation in general — using the same
//! entropy-referenced Riemann invariant for every leg doesn't automatically
//! conserve stagnation enthalpy across the junction at finite amplitude
//! (only to `O((u/c)^2)`, the small-disturbance/linear-acoustics regime
//! already validated by `tests/acoustic_resonance.rs`). This is not an
//! implementation gap to engineer away - it is the same reason Gordon
//! Blair's own more-accurate model (Sec. 2.14, the non-isentropic
//! Bingham-Blair branch model) needs a *second*, separately-solved energy
//! equation and an angle-dependent pressure-loss term, well beyond this
//! module's deliberately narrow scope. See `tests/branch_junction_reflection.rs`
//! for the measured energy residual at both small and realistic amplitude.
//!
//! # Scope
//!
//! Lossless only (no discharge coefficient) - matches OpenWAM's own
//! `TCCRamificacion` and Blair's "constant pressure" model exactly. A
//! future per-leg loss coefficient would need to enter through the
//! momentum/energy balance the way Blair's angle-dependent `CL(theta)`
//! does, not as a naive multiplier the way `valve.rs`'s `Cd` works (there
//! is no orifice here) - its own future design pass, not a one-line
//! addition to this module.
//!
//! This module never learns what a [`crate::cylinder::Cylinder`] or
//! [`crate::valve::ValveGeometry`] is (mirrors the architectural boundary
//! already established for `valve_port.rs`), and `network.rs` never learns
//! what a [`BranchJunction`] is - it only consumes the
//! [`crate::network::ExternalPortFlux`] list [`resolve_branch_junction`]
//! produces.

use crate::gas::{GasProperties, PrimitiveState};
use crate::network::{ExternalPortFlux, PipeEnd, PipeEndRef, PipeNetwork};

/// Several pipe ends meeting at one point, sharing one instantaneous
/// static pressure (see the module doc comment). `Junction` remains the
/// right choice for an exact, same-area, 2-pipe connection (a genuine
/// 2-state HLLC solve, strictly more accurate than this approximate
/// closure) - use `BranchJunction` for `ends.len() != 2` or mismatched-area
/// cases, which as a side effect also covers the previously-unsupported
/// "genuine sudden area change" 2-pipe junction.
#[derive(Debug, Clone, PartialEq)]
pub struct BranchJunction {
    pub ends: Vec<PipeEndRef>,
}

fn end_sign(end: PipeEnd) -> f64 {
    match end {
        PipeEnd::Right => 1.0,
        PipeEnd::Left => -1.0,
    }
}

/// One pipe end's own state and derived quantities, fixed for the whole
/// resolve (independent of the trial junction pressure).
///
/// `pub(crate)` so `lossy_branch_junction.rs` can reuse this exact,
/// already-tested representation and its associated functions below rather
/// than risk a second, silently-diverging transcription of the same
/// isentropic relation.
pub(crate) struct Leg {
    pub(crate) pipe_index: usize,
    pub(crate) end: PipeEnd,
    pub(crate) sign: f64,
    pub(crate) area: f64,
    pub(crate) pressure: f64,
    pub(crate) sound_speed: f64,
    /// The leg's own known Riemann invariant: `u + sign*2*c/(gamma-1)`.
    pub(crate) riemann_invariant: f64,
}

pub(crate) fn collect_legs(network: &PipeNetwork, junction: &BranchJunction, gas: &GasProperties) -> Vec<Leg> {
    junction
        .ends
        .iter()
        .map(|end_ref| {
            let pipe = &network.pipes[end_ref.pipe_index];
            let (state, area) = match end_ref.end {
                PipeEnd::Left => (pipe.left_boundary_cell_state(gas), pipe.mesh.face_areas[0]),
                PipeEnd::Right => {
                    let face_areas = &pipe.mesh.face_areas;
                    (pipe.right_boundary_cell_state(gas), face_areas[face_areas.len() - 1])
                }
            };
            let sign = end_sign(end_ref.end);
            let sound_speed = state.sound_speed(gas);
            let riemann_invariant = state.velocity + sign * 2.0 * sound_speed / (gas.gamma - 1.0);
            Leg {
                pipe_index: end_ref.pipe_index,
                end: end_ref.end,
                sign,
                area,
                pressure: state.pressure,
                sound_speed,
                riemann_invariant,
            }
        })
        .collect()
}

/// One leg's resolved boundary state at a given trial junction pressure.
#[derive(Clone, Copy)]
pub(crate) struct LegSolution {
    pub(crate) velocity: f64,
    pub(crate) density: f64,
    /// Signed mass flow rate, positive = this leg supplies the junction.
    pub(crate) mass_flow_out: f64,
}

/// Every leg (supplying or receiving) is resolved the same way: isentropic
/// expansion/compression of its own gas from its own current state to the
/// trial junction pressure, using its own entropy and its own known
/// Riemann invariant. See the module doc comment for why this - and not a
/// mixed-entropy or mixed-stagnation reference - is the correct closure.
pub(crate) fn leg_solution(leg: &Leg, gas: &GasProperties, trial_pressure: f64) -> LegSolution {
    let sound_speed = leg.sound_speed * (trial_pressure / leg.pressure).powf((gas.gamma - 1.0) / (2.0 * gas.gamma));
    let velocity = leg.riemann_invariant - leg.sign * 2.0 * sound_speed / (gas.gamma - 1.0);
    let density = gas.gamma * trial_pressure / (sound_speed * sound_speed);
    let mass_flow_out = leg.sign * density * velocity * leg.area;
    LegSolution { velocity, density, mass_flow_out }
}

/// Resolves every leg at the given trial junction pressure.
pub(crate) fn resolve_all_legs(legs: &[Leg], gas: &GasProperties, trial_pressure: f64) -> Vec<LegSolution> {
    legs.iter().map(|leg| leg_solution(leg, gas, trial_pressure)).collect()
}

pub(crate) fn total_mass_flow(legs: &[Leg], gas: &GasProperties, trial_pressure: f64) -> f64 {
    resolve_all_legs(legs, gas, trial_pressure).iter().map(|s| s.mass_flow_out).sum()
}

/// Solves for the shared junction pressure via bisection: the equation
/// being solved is exact mass conservation, `sum(signed mass flow) = 0`.
///
/// Starts from the bracket `[min(leg pressures), max(leg pressures)]`,
/// which is a *proven* bracket when every leg's own velocity is zero (at
/// the low end, every other leg has strictly higher pressure and so
/// supplies non-negatively; symmetrically at the high end) - the common
/// case for every validation case in this crate (near-quiescent legs with
/// small perturbations). For legs with real nonzero velocity, the true
/// root can in principle sit slightly outside that raw pressure range (a
/// leg already blowing outward can keep supplying against a slightly
/// higher junction pressure) - handled by widening the bracket
/// defensively rather than assuming it never happens, and panicking with
/// a clear message if widening still fails to bracket a root, rather than
/// silently returning a wrong pressure.
pub(crate) fn solve_junction_pressure(legs: &[Leg], gas: &GasProperties) -> f64 {
    let min_pressure = legs.iter().map(|l| l.pressure).fold(f64::INFINITY, f64::min);
    let max_pressure = legs.iter().map(|l| l.pressure).fold(f64::NEG_INFINITY, f64::max);

    if (max_pressure - min_pressure) < 1e-9 * max_pressure {
        return 0.5 * (min_pressure + max_pressure);
    }

    let mut lo = min_pressure;
    let mut hi = max_pressure;
    let mut total_lo = total_mass_flow(legs, gas, lo);
    let mut total_hi = total_mass_flow(legs, gas, hi);

    let mut widen_attempts = 0;
    while total_lo.signum() == total_hi.signum() && widen_attempts < 20 {
        let span = hi - lo;
        lo = (lo - 0.5 * span).max(1.0);
        hi += 0.5 * span;
        total_lo = total_mass_flow(legs, gas, lo);
        total_hi = total_mass_flow(legs, gas, hi);
        widen_attempts += 1;
    }
    assert!(
        total_lo.signum() != total_hi.signum(),
        "branch junction: could not bracket a root for the shared junction pressure after widening \
         (total mass flow {total_lo:e} at {lo:e} Pa, {total_hi:e} at {hi:e} Pa)"
    );

    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if (hi - lo) < 1e-12 * mid {
            return mid;
        }
        let total_mid = total_mass_flow(legs, gas, mid);
        if total_mid == 0.0 {
            return mid;
        }
        if total_mid.signum() == total_lo.signum() {
            lo = mid;
            total_lo = total_mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Solves the branch junction and returns one [`ExternalPortFlux`] per
/// leg, ready to hand to [`PipeNetwork::advance_with_external_fluxes`] -
/// mirrors `valve_port.rs`'s own "compute once at the pre-step state,
/// apply once" pattern (no sub-stepping - see that module's doc comment
/// for why re-evaluating mid-step would break the same conservation
/// invariant this module relies on).
pub fn resolve_branch_junction(
    network: &PipeNetwork,
    junction: &BranchJunction,
    gas: &GasProperties,
) -> Vec<ExternalPortFlux> {
    assert!(junction.ends.len() >= 2, "a branch junction needs at least 2 pipe ends, got {}", junction.ends.len());

    let legs = collect_legs(network, junction, gas);
    let junction_pressure = solve_junction_pressure(&legs, gas);
    let solutions = resolve_all_legs(&legs, gas, junction_pressure);

    let typical_mass_flow_scale =
        solutions.iter().map(|s| s.mass_flow_out.abs()).fold(0.0, f64::max).max(1e-9);
    debug_assert!(
        solutions.iter().map(|s| s.mass_flow_out).sum::<f64>().abs() < 1e-6 * typical_mass_flow_scale,
        "branch junction mass conservation invariant violated at the converged root"
    );

    legs.iter()
        .zip(solutions)
        .map(|(leg, solution)| {
            let face_state =
                PrimitiveState { density: solution.density, velocity: solution.velocity, pressure: junction_pressure };
            let flux = face_state.to_conserved(gas).physical_flux(gas);
            ExternalPortFlux { end: PipeEndRef { pipe_index: leg.pipe_index, end: leg.end }, neighbor_state: face_state, flux }
        })
        .collect()
}

/// Advances a [`PipeNetwork`] containing one or more [`BranchJunction`]s by
/// one explicit, CFL-limited timestep, resolving every branch junction's
/// shared pressure once at the pre-step state (mirrors
/// [`crate::solver::step`], generalized to also carry branch-junction
/// fluxes through [`PipeNetwork::advance_with_external_fluxes`]). Returns
/// the `dt` actually taken.
pub fn step(network: &mut PipeNetwork, gas: &GasProperties, cfl: f64, junctions: &[BranchJunction]) -> f64 {
    let dt = network.cfl_time_step(gas, cfl);
    let external: Vec<ExternalPortFlux> =
        junctions.iter().flat_map(|junction| resolve_branch_junction(network, junction, gas)).collect();
    network.advance_with_external_fluxes(dt, gas, &external);
    dt
}

/// Repeatedly calls [`step`] until at least `t_end` seconds have elapsed.
/// Mirrors [`crate::solver::run_to_time`]'s pairing with
/// [`crate::solver::step`].
pub fn run_to_time(network: &mut PipeNetwork, gas: &GasProperties, cfl: f64, t_end: f64, junctions: &[BranchJunction]) -> f64 {
    let mut elapsed = 0.0;
    while elapsed < t_end {
        elapsed += step(network, gas, cfl, junctions);
    }
    elapsed
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary::BoundaryCondition;
    use crate::mesh::Mesh;
    use crate::pipe::Pipe;

    fn air_at_rest(pressure: f64) -> PrimitiveState {
        PrimitiveState::from_pressure_temperature(pressure, 293.15, 0.0, &GasProperties::AIR)
    }

    fn uniform_pipe(pressure: f64, diameter: f64) -> Pipe {
        Pipe::uniform_initial_state(
            Mesh::uniform(1.0, diameter, 0.01),
            air_at_rest(pressure),
            &GasProperties::AIR,
            BoundaryCondition::ClosedEnd,
            BoundaryCondition::ClosedEnd,
        )
    }

    #[test]
    fn matched_states_across_a_branch_produce_zero_net_flux_at_every_leg() {
        let gas = GasProperties::AIR;
        let network = PipeNetwork {
            pipes: vec![uniform_pipe(150_000.0, 0.05), uniform_pipe(150_000.0, 0.05), uniform_pipe(150_000.0, 0.05)],
            junctions: vec![],
        };
        let junction = BranchJunction {
            ends: vec![
                PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
                PipeEndRef { pipe_index: 1, end: PipeEnd::Left },
                PipeEndRef { pipe_index: 2, end: PipeEnd::Left },
            ],
        };
        let fluxes = resolve_branch_junction(&network, &junction, &gas);
        for flux in &fluxes {
            assert!(flux.flux.mass.abs() < 1e-9, "expected zero mass flux, got {}", flux.flux.mass);
        }
    }

    #[test]
    fn a_disturbed_leg_among_matched_legs_produces_outflow_from_it_and_inflow_to_the_others() {
        let gas = GasProperties::AIR;
        let network = PipeNetwork {
            pipes: vec![uniform_pipe(200_000.0, 0.05), uniform_pipe(150_000.0, 0.05), uniform_pipe(150_000.0, 0.05)],
            junctions: vec![],
        };
        let junction = BranchJunction {
            ends: vec![
                PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
                PipeEndRef { pipe_index: 1, end: PipeEnd::Left },
                PipeEndRef { pipe_index: 2, end: PipeEnd::Left },
            ],
        };
        let fluxes = resolve_branch_junction(&network, &junction, &gas);

        // `flux.mass` is always "quantity moving in the +x direction
        // through this face" (the pipe's own raw coordinate convention,
        // not an "outward positive" one) - at a Left end, +x points
        // FURTHER INTO the pipe, so a receiving (inflow) leg there shows
        // POSITIVE flux.mass, matching a Right end's own positive-outflow
        // convention (this mirrors the sign convention `valve_port.rs`
        // already established: positive flux at a Left face = inflow).
        // Leg 0 (Right end, higher pressure) supplies the junction
        // (outflow, positive); legs 1/2 (Left ends, lower pressure,
        // identical) both receive (inflow, also positive at a Left end),
        // split evenly by symmetry.
        assert!(fluxes[0].flux.mass > 0.0, "expected leg 0 to supply the junction");
        assert!(fluxes[1].flux.mass > 0.0, "expected leg 1 to receive from the junction");
        assert!(fluxes[2].flux.mass > 0.0, "expected leg 2 to receive from the junction");
        let relative_difference = (fluxes[1].flux.mass - fluxes[2].flux.mass).abs() / fluxes[1].flux.mass.abs();
        assert!(relative_difference < 1e-9, "identical legs should receive identically by symmetry");
    }

    #[test]
    fn total_mass_flow_is_exactly_zero_at_the_converged_root() {
        let gas = GasProperties::AIR;
        let network = PipeNetwork {
            pipes: vec![
                uniform_pipe(300_000.0, 0.04),
                uniform_pipe(100_000.0, 0.05),
                uniform_pipe(120_000.0, 0.06),
                uniform_pipe(90_000.0, 0.03),
            ],
            junctions: vec![],
        };
        let junction = BranchJunction {
            ends: vec![
                PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
                PipeEndRef { pipe_index: 1, end: PipeEnd::Left },
                PipeEndRef { pipe_index: 2, end: PipeEnd::Left },
                PipeEndRef { pipe_index: 3, end: PipeEnd::Right },
            ],
        };
        let fluxes = resolve_branch_junction(&network, &junction, &gas);

        let total: f64 = fluxes
            .iter()
            .map(|f| {
                let sign = end_sign(f.end.end);
                let area = network.pipes[f.end.pipe_index].mesh.face_areas[0].max(
                    network.pipes[f.end.pipe_index].mesh.face_areas
                        [network.pipes[f.end.pipe_index].mesh.face_areas.len() - 1],
                );
                sign * f.flux.mass * area
            })
            .sum();
        let scale: f64 = fluxes
            .iter()
            .map(|f| {
                let area = network.pipes[f.end.pipe_index].mesh.face_areas[0].max(
                    network.pipes[f.end.pipe_index].mesh.face_areas
                        [network.pipes[f.end.pipe_index].mesh.face_areas.len() - 1],
                );
                (f.flux.mass * area).abs()
            })
            .fold(0.0, f64::max);
        let relative_error = total.abs() / scale;
        println!("4-leg branch junction mass conservation: relative error {relative_error:e}");
        assert!(relative_error < 1e-9, "mass not conserved across the junction: relative error {relative_error:e}");
    }

    #[test]
    fn four_way_junction_with_uniform_area_matches_uniform_state_with_zero_flow() {
        // All 4 legs identical - the shared pressure should equal every
        // leg's own pressure exactly, with zero flow everywhere.
        let gas = GasProperties::AIR;
        let network = PipeNetwork {
            pipes: vec![uniform_pipe(180_000.0, 0.04); 4],
            junctions: vec![],
        };
        let junction = BranchJunction {
            ends: vec![
                PipeEndRef { pipe_index: 0, end: PipeEnd::Right },
                PipeEndRef { pipe_index: 1, end: PipeEnd::Left },
                PipeEndRef { pipe_index: 2, end: PipeEnd::Left },
                PipeEndRef { pipe_index: 3, end: PipeEnd::Left },
            ],
        };
        let fluxes = resolve_branch_junction(&network, &junction, &gas);
        for flux in &fluxes {
            assert!(flux.flux.mass.abs() < 1e-9);
            assert!((flux.flux.momentum - 180_000.0).abs() < 1e-3);
        }
    }
}
