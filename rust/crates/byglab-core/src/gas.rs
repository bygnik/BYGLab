//! Ideal-gas thermodynamics for a calorically perfect gas (constant ratio of
//! specific heats). This is the state representation the whole 1D solver is
//! built on: [`PrimitiveState`] for physically-intuitive quantities
//! (density, velocity, pressure), [`ConservedState`] for the quantities the
//! finite-volume update actually integrates (mass, momentum, energy
//! densities), and [`Flux`] for the rate at which those conserved quantities
//! pass through a point.
//!
//! `Flux` is a distinct type from `ConservedState` even though both are
//! three-component vectors in the same units-per-component order — a flux is
//! a *rate through a surface*, not a *state at a point*, and keeping them
//! separate lets the compiler catch a state accidentally used where a flux
//! was meant (or vice versa).

use serde::{Deserialize, Serialize};

/// Properties of the working fluid: ratio of specific heats and the
/// specific gas constant. Constant (not a function of temperature or
/// composition) — this models a calorically perfect gas, matching the
/// `nmGammaConstante` mode used throughout the OpenWAM reference cases this
/// solver is validated against.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GasProperties {
    /// Ratio of specific heats, cp/cv (dimensionless). 1.4 for diatomic
    /// gases like air at moderate temperatures.
    pub gamma: f64,
    /// Specific gas constant R = cp - cv, in J/(kg*K).
    pub gas_constant: f64,
}

impl GasProperties {
    /// Dry air at moderate temperature: gamma = 1.4, R = 287 J/(kg*K).
    /// Matches the gas model used in every OpenWAM reference case.
    pub const AIR: GasProperties = GasProperties { gamma: 1.4, gas_constant: 287.0 };
}

/// Flow state expressed as the physically-intuitive quantities: density,
/// velocity, and pressure. Convenient for specifying initial/boundary
/// conditions and for reading results, but NOT what the finite-volume
/// update integrates directly — see [`ConservedState`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PrimitiveState {
    /// kg/m^3
    pub density: f64,
    /// m/s, positive in the direction of increasing pipe-axis position
    pub velocity: f64,
    /// Pa (absolute, not gauge)
    pub pressure: f64,
}

impl PrimitiveState {
    /// Builds a state at rest (or with a given velocity) from pressure and
    /// temperature — the natural way to specify initial/boundary conditions
    /// for engine gas dynamics, where pressure and temperature are the
    /// measured/known quantities rather than density.
    ///
    /// Density follows from the ideal gas law: `density = pressure / (R * temperature_kelvin)`.
    pub fn from_pressure_temperature(
        pressure: f64,
        temperature_kelvin: f64,
        velocity: f64,
        gas: &GasProperties,
    ) -> Self {
        let density = pressure / (gas.gas_constant * temperature_kelvin);
        PrimitiveState { density, velocity, pressure }
    }

    /// Local speed of sound, `sqrt(gamma * pressure / density)`.
    ///
    /// Always well-defined and positive for a physical gas state (density
    /// and pressure are both positive), so this never risks a
    /// divide-by-zero the way computing it from temperature via density
    /// might if density were zero.
    pub fn sound_speed(&self, gas: &GasProperties) -> f64 {
        (gas.gamma * self.pressure / self.density).sqrt()
    }

    /// Temperature in Kelvin, back-derived from the ideal gas law.
    pub fn temperature_kelvin(&self, gas: &GasProperties) -> f64 {
        self.pressure / (gas.gas_constant * self.density)
    }

    /// Converts to the conserved (mass/momentum/energy density) form used
    /// internally by the finite-volume update.
    pub fn to_conserved(&self, gas: &GasProperties) -> ConservedState {
        let kinetic_energy_density = 0.5 * self.density * self.velocity * self.velocity;
        let internal_energy_density = self.pressure / (gas.gamma - 1.0);
        ConservedState {
            mass: self.density,
            momentum: self.density * self.velocity,
            energy: internal_energy_density + kinetic_energy_density,
        }
    }
}

/// Flow state expressed as the quantities a finite-volume scheme actually
/// conserves and integrates: mass, momentum, and total energy, all per unit
/// volume. This is what's stored per cell and updated each timestep.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ConservedState {
    /// Mass density, rho, kg/m^3. Same quantity as [`PrimitiveState::density`].
    pub mass: f64,
    /// Momentum density, rho*u, kg/(m^2*s).
    pub momentum: f64,
    /// Total energy density (internal + kinetic), J/m^3.
    pub energy: f64,
}

impl ConservedState {
    /// Converts to the primitive (density/velocity/pressure) form, the
    /// natural form for reading results or evaluating boundary conditions.
    pub fn to_primitive(&self, gas: &GasProperties) -> PrimitiveState {
        let velocity = self.momentum / self.mass;
        let kinetic_energy_density = 0.5 * self.mass * velocity * velocity;
        let internal_energy_density = self.energy - kinetic_energy_density;
        let pressure = internal_energy_density * (gas.gamma - 1.0);
        PrimitiveState { density: self.mass, velocity, pressure }
    }

    /// The physical flux `F(U)` of this state: the rate at which mass,
    /// momentum, and energy would pass through a stationary point if the
    /// flow here extended uniformly past it. This is what a Riemann solver
    /// approximates *across* a discontinuity between two different states;
    /// evaluated on a single state it's the exact analytic flux, used to
    /// check that a numerical flux function is consistent
    /// (`numerical_flux(U, U) == U.physical_flux()`).
    pub fn physical_flux(&self, gas: &GasProperties) -> Flux {
        let primitive = self.to_primitive(gas);
        Flux {
            mass: self.momentum,
            momentum: self.momentum * primitive.velocity + primitive.pressure,
            energy: primitive.velocity * (self.energy + primitive.pressure),
        }
    }
}

/// The rate at which mass, momentum, and energy pass through a point in
/// space — as opposed to [`ConservedState`], which is the amount present
/// *at* a point. Produced by [`crate::riemann::hllc_flux`] at cell
/// interfaces and by [`ConservedState::physical_flux`] for a single state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Flux {
    /// Mass flux, kg/(m^2*s).
    pub mass: f64,
    /// Momentum flux, N/m^2 (i.e. Pa, since momentum flux has units of
    /// pressure — this is `rho*u^2 + p`).
    pub momentum: f64,
    /// Energy flux, W/m^2.
    pub energy: f64,
}

impl std::ops::Sub for Flux {
    type Output = Flux;
    fn sub(self, rhs: Flux) -> Flux {
        Flux {
            mass: self.mass - rhs.mass,
            momentum: self.momentum - rhs.momentum,
            energy: self.energy - rhs.energy,
        }
    }
}

impl std::ops::Mul<f64> for Flux {
    type Output = Flux;
    fn mul(self, scalar: f64) -> Flux {
        Flux { mass: self.mass * scalar, momentum: self.momentum * scalar, energy: self.energy * scalar }
    }
}

impl ConservedState {
    /// Updates this state by one explicit finite-volume step, given the
    /// fluxes at the cell's left and right faces, the cell width, and the
    /// timestep. This is the core update formula of the whole solver:
    /// `dU/dt = -(F_right - F_left) / dx`, applied explicitly (forward Euler).
    pub fn advance(&mut self, left_face_flux: Flux, right_face_flux: Flux, cell_width: f64, dt: f64) {
        let divergence = right_face_flux - left_face_flux;
        let scale = dt / cell_width;
        self.mass -= divergence.mass * scale;
        self.momentum -= divergence.momentum * scale;
        self.energy -= divergence.energy * scale;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primitive_to_conserved_round_trips() {
        let gas = GasProperties::AIR;
        let original = PrimitiveState { density: 1.2, velocity: 15.0, pressure: 101_325.0 };
        let recovered = original.to_conserved(&gas).to_primitive(&gas);
        assert!((original.density - recovered.density).abs() < 1e-9);
        assert!((original.velocity - recovered.velocity).abs() < 1e-9);
        assert!((original.pressure - recovered.pressure).abs() < 1e-6);
    }

    #[test]
    fn from_pressure_temperature_matches_ideal_gas_law() {
        let gas = GasProperties::AIR;
        let state = PrimitiveState::from_pressure_temperature(101_325.0, 293.15, 0.0, &gas);
        // rho = p / (R T)
        let expected_density = 101_325.0 / (287.0 * 293.15);
        assert!((state.density - expected_density).abs() < 1e-9);
    }

    #[test]
    fn sound_speed_of_air_at_room_temperature_is_about_343() {
        let gas = GasProperties::AIR;
        let state = PrimitiveState::from_pressure_temperature(101_325.0, 293.15, 0.0, &gas);
        let c = state.sound_speed(&gas);
        assert!((c - 343.2).abs() < 0.5, "expected ~343.2 m/s, got {c}");
    }

    #[test]
    fn physical_flux_of_state_at_rest_is_pressure_only() {
        let gas = GasProperties::AIR;
        let state = PrimitiveState { density: 1.2, velocity: 0.0, pressure: 150_000.0 };
        let flux = state.to_conserved(&gas).physical_flux(&gas);
        assert_eq!(flux.mass, 0.0);
        assert!((flux.momentum - 150_000.0).abs() < 1e-6);
        assert_eq!(flux.energy, 0.0);
    }
}
