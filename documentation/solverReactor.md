
# solverReactor
## Overview

The module implements a one-dimensional steady-state reacting-flow solver based on a finite-volume discretization of species and temperature transport in a cylindrical reactor. It includes thermodynamic species models, a single-reaction kinetic model, mixture property evaluation, zone-based mesh generation, inlet and outlet boundary handling, a segregated steady-state solver, plotting utilities, and a factory function that builds a full reactor case from a context dictionary. [file:63]

The implementation is aimed at plug-flow-like simulations with axial transport only, constant mass flux, ideal-gas or constant-density mixture handling, and cellwise reaction and heating activation through mesh zones. The core solver advances species and temperature iteratively until a steady state is reached. [file:63]

## Module structure

The module is organized around these main components: [file:63]

- `Specie`: thermodynamic data and temperature-dependent heat capacity/enthalpy evaluation. [file:63]
- `Reaction`: kinetics, equilibrium, mass sources, heat release, and analytical derivatives. [file:63]
- `Mixture`: mixture molar mass, density, and heat capacity evaluation. [file:63]
- `domainSetup`: geometric container for the reactor diameter. [file:63]
- `Inlet` and `Outlet`: inlet boundary conditions and outlet result extraction. [file:63]
- `Zone` and `Mesh`: zone-wise reactor layout and 1D discretization. [file:63]
- `scalarField`: cell-centered scalar storage. [file:63]
- `solver`: the segregated steady-state reacting-flow solver. [file:63]
- `ReactorPlotter`: PNG plot generation for temperature and species. [file:63]
- `build_reactor_from_context()`: high-level case builder from a YAML-like context dictionary. [file:63]

## Constants and assumptions

The module defines the universal gas constant, a reference temperature, a reference pressure, and `PI` at module scope. These constants are used in thermodynamic and kinetic calculations across the whole solver. [file:63]

The most important modeling assumptions are: [file:63]

- One-dimensional axial transport only. [file:63]
- Steady-state solution. [file:63]
- Single chemical reaction currently supported. [file:63]
- Constant total mass flux computed from the inlet. [file:63]
- Ideal-gas density at reference pressure by default, with optional constant density. [file:63]
- Segregated solution of species and temperature rather than a fully coupled block solve. [file:63]

## Main classes

### `Specie`

`Specie` stores the thermodynamic description of one chemical species. It includes the species name, molar mass, enthalpy of formation, entropy of formation, and either a constant heat-capacity model or a polynomial heat-capacity model. [file:63]

Important features: [file:63]

- Automatic integer ID assignment through a class counter. [file:63]
- Construction from a dictionary with `Specie.from_dict()`. [file:63]
- `heatCapacity(T)` for constant or polynomial heat capacity evaluation. [file:63]
- `enthalpy(T)` for formation plus sensible enthalpy. [file:63]

The polynomial heat-capacity mode builds a NumPy polynomial and its antiderivative during initialization. This allows fast repeated evaluation of both `cp(T)` and sensible enthalpy over full cell arrays. [file:63]

### `Reaction`

`Reaction` represents the single reaction active in the solver. It stores stoichiometric coefficients, forward and reverse exponents, reversibility, Arrhenius parameters, and the participating species list. [file:63]

Key capabilities include: [file:63]

- Forward Arrhenius rate constant evaluation. [file:63]
- Equilibrium constant evaluation from reaction enthalpy and entropy. [file:63]
- Backward rate constant evaluation through `K_p` to `K_c` conversion. [file:63]
- Forward and backward mass-action rate evaluation. [file:63]
- Species mass source evaluation. [file:63]
- Reaction heat source evaluation. [file:63]
- Analytical derivatives with respect to temperature and concentration. [file:63]

#### Reaction source model

The reaction model uses concentrations in mol/m³, computed from local density, mass fraction, and molar mass. The net rate is `rateForward - rateBackward`, and mass sources are obtained by multiplying by stoichiometric coefficients and species molar masses. [file:63]

The heat source is computed from the net rate and the temperature-dependent reaction enthalpy. The module also supports an equilibrium mask that sets near-equilibrium cells to zero in selected source calculations to reduce stiffness. [file:63]

### `Mixture`

`Mixture` computes thermophysical mixture properties from the active species list and cellwise species mass fractions. It stores the selected density model and the species molar masses. [file:63]

Main methods: [file:63]

- `equivalentMolarMass(speciesFractions)` computes the mixture molar mass. [file:63]
- `idealGasDensity(T, speciesFractions)` computes density from the ideal gas law. [file:63]
- `mixtureHeatCapacity(T, speciesFractions)` computes mass-weighted heat capacity. [file:63]

By default, density is based on the ideal-gas law at `P_REF`, but the inlet logic can also use a constant density if `densityModel == 'const'`. [file:63]

### `domainSetup`

`domainSetup` is a lightweight container for the reactor geometry. In the current implementation, it mainly stores the reactor diameter, which is then used to compute cross-sectional area, inlet mass flow rate, and cell volumes. [file:63]

### `Inlet`

`Inlet` defines the inlet boundary condition in terms of axial velocity, temperature, and species mass fractions. It also contains `from_dict()` to reorder species according to the species list used in the solver and to normalize the inlet mass fractions so they sum to one. [file:63]

The method `inletValues(mixture, domain)` computes: [file:63]

- Inlet density. [file:63]
- Mass flow rate. [file:63]
- Temperature boundary value. [file:63]
- Ordered inlet species fraction vector. [file:63]

### `Outlet`

`Outlet` is a result container that stores the solution extracted from a selected cell, typically the last one. It stores temperature, mass fractions, density, velocity, mass flow rate, and molar concentrations. [file:63]

Useful methods include: [file:63]

- `fromSolver(slv, position=-1)` to build an outlet object from the current solver state. [file:63]
- `asDict(species=None)` to serialize output data. [file:63]
- `massFraction(name, species)` and `concentration(name, species)` for species-specific access. [file:63]

### `Zone`

`Zone` represents one axial segment of the reactor. It stores the zone length, a zone type string, and flags indicating whether heat and mass source terms are active in that segment. [file:63]

Each zone can therefore act as: [file:63]

- A reaction-active region. [file:63]
- An externally heated region. [file:63]
- A passive region. [file:63]

`Zone.from_dict()` maps YAML-style configuration to the internal object and copies `heatSource`, `massSource`, and `heatSourceValue` into the zone object. [file:63]

### `Mesh`

`Mesh` constructs a one-dimensional finite-volume mesh from an ordered list of zones and the domain diameter. Each zone is divided into uniformly sized cells based on the target mesh spacing. [file:63]

The method `meshCreate()` allocates and fills cell-level arrays such as: [file:63]

- `cell_centers` [file:63]
- `cell_sizes` [file:63]
- `cell_volumes` [file:63]
- `cell_zone_id` [file:63]
- `cell_zone_type` [file:63]
- `cell_heat_flag` [file:63]
- `cell_mass_flag` [file:63]
- `cell_heat_value` [file:63]

This provides the bridge between zone-level configuration and cellwise PDE assembly. `Mesh.from_dict()` reads the target sizing from the mesh configuration block and builds the final discretized mesh. [file:63]

### `scalarField`

`scalarField` is a thin wrapper around a cell-centered NumPy array. It is used for species, temperature, and velocity fields. [file:63]

The field is allocated by calling `fieldInitialize(mesh)`, which creates a zero array aligned with `mesh.cell_centers`. The solver then stacks species fields into a 2D matrix and keeps temperature and velocity as dedicated `scalarField` objects. [file:63]

### `solver`

`solver` is the main numerical engine of the module. It solves the species and energy equations on the 1D mesh for the given mixture, reaction, and inlet conditions. [file:63]

At construction, it: [file:63]

- Stores mesh, mixture, reaction, and inlet objects. [file:63]
- Computes the inlet mass flux from the inlet boundary condition. [file:63]
- Stacks species fields into a `(n_species, n_cells)` array. [file:63]
- Allocates temperature and velocity fields. [file:63]
- Initializes density, source arrays, residual arrays, and reaction-rate storage. [file:63]

## Solver workflow

### Initialization

`initializeCase()` assigns the inlet composition to all cells, sets the entire temperature field to the inlet temperature, and updates density. This gives a physically consistent plug-flow initial condition before starting the outer iterations. [file:63]

### Concentration evaluation

`concentrationArray()` converts mass fractions to molar concentrations using the current density field and species molar masses. The returned array has shape `(n_species, n_cells)`. [file:63]

### Density update

`update_density()` recomputes the cellwise density from current temperature and composition. It is called after composition or temperature changes so that all property calculations remain consistent. [file:63]

### Source-term evaluation

`sourcesEvaluation()` is one of the most important methods in the script. It evaluates reaction rates, computes species and heat source terms, applies zone masks, and under-relaxes the chemistry source terms between nonlinear iterations. [file:63]

Its workflow is: [file:63]

1. Reset source arrays. [file:63]
2. Compute concentrations, density, and temperature references. [file:63]
3. Read reaction-active and heat-active masks from the mesh. [file:63]
4. Add prescribed external heat from `cell_heat_value`. [file:63]
5. Compute forward and backward reaction rates. [file:63]
6. Build a near-equilibrium mask. [file:63]
7. Compute species mass sources and their derivatives. [file:63]
8. Compute reaction heat sources and `dQ/dT`. [file:63]
9. Apply source under-relaxation and store the final arrays. [file:63]

### Species equation assembly

`matrixSpecieEquationAssembly(specieIndex)` assembles a sparse lower-bidiagonal convection-reaction system for one species using `scipy.sparse.diags`. The source linearization follows the common split `S = S_U + S_P * Y`, with `S_P = min(dS/dY, 0)` to preserve stability. [file:63]

The boundary condition is imposed at the first cell through the inlet mass fraction. The method returns the sparse matrix `A` and right-hand side vector `b`. [file:63]

### Temperature equation assembly

`matrixTemperatureEquationAssembly()` builds the energy equation in the same finite-volume style, with convective transport based on `massFlux * cp_mix` and implicit source treatment from `dQ/dT`. [file:63]

Like the species equation, it returns a sparse matrix and a right-hand side. The inlet temperature is applied in the first control volume. [file:63]

### Residual evaluation

The module also provides explicit residual evaluators: [file:63]

- `specieScalarEquation()` computes the species residual matrix. [file:63]
- `heatEquation()` computes the energy residual vector. [file:63]

These are useful for diagnostics and convergence monitoring, even though the main outer loop currently monitors field change rather than residual norm. [file:63]

### Steady-state loop

`steadyState()` runs the segregated outer iteration. The procedure is: [file:63]

1. Update density and velocity. [file:63]
2. Evaluate source terms. [file:63]
3. Loop over all species and solve each species equation with `spsolve`. [file:63]
4. Clip species to `[0, 1]` and renormalize them so the cellwise sum is one. [file:63]
5. Update density and velocity again. [file:63]
6. Re-evaluate source terms. [file:63]
7. Solve the temperature equation with `spsolve`. [file:63]
8. Clip temperature between `temperatureClipLow` and `temperatureClipHigh`. [file:63]
9. Under-relax the temperature update. [file:63]
10. Measure `dY` and `dT`, print diagnostics, and stop if the convergence criterion is reached. [file:63]

On completion, the method stores `self.outlet = Outlet.fromSolver(self)` and returns the outlet object. [file:63]

## Plotting utilities

`ReactorPlotter` generates PNG output from a converged or partially converged solver state. It operates directly on the solver object and uses the mesh cell centers as the axial coordinate. [file:63]

Available methods: [file:63]

- `get_axis()` returns the axial coordinate array. [file:63]
- `save_temperature(path)` writes a temperature profile PNG. [file:63]
- `save_species(path)` writes a species mass-fraction profile PNG. [file:63]
- `save_all(path)` writes a side-by-side temperature and species summary PNG. [file:63]

This makes the solver suitable for automated case directories where each optimization sample can save its own diagnostic plots. [file:63]

## Factory function

### `build_reactor_from_context(case_ctx)`

This function builds a complete reactor simulation from a nested configuration dictionary. It resets `Specie.counter` and `Zone.counter`, creates species, creates the single reaction, builds the mixture, domain, zones, mesh, and inlet, initializes one `scalarField` per species, and then constructs the `solver` object. [file:63]

The function returns: [file:63]

- `slv`: the configured solver instance. [file:63]
- `species`: the ordered species list used across the case. [file:63]

This is the main integration point used by higher-level driver scripts such as optimization workflows. [file:63]

## Data layout

The module consistently uses NumPy arrays with parameters organized per species and per cell. [file:63]

| Quantity | Shape | Description |
| --- | --- | --- |
| `specieFields` | `(n_species, n_cells)` | Species mass fractions in all cells. [file:63] |
| `temperatureField.cellField` | `(n_cells,)` | Temperature field. [file:63] |
| `velocityField.cellField` | `(n_cells,)` | Velocity field. [file:63] |
| `density` | `(n_cells,)` | Mixture density. [file:63] |
| `massSources` | `(n_species, n_cells)` | Species source terms. [file:63] |
| `massSourcesDerivative` | `(n_species, n_cells)` | Species source Jacobian terms. [file:63] |
| `heatSources` | `(n_cells,)` | Total volumetric heat source. [file:63] |
| `heatSourcesDerivative` | `(n_cells,)` | Heat source derivative with respect to temperature. [file:63] |
| `reactionRates` | `(2, n_cells)` | Forward and backward rates. [file:63] |

## Minimal usage example

```python
from scripts.solverReactor import build_reactor_from_context, ReactorPlotter

case_ctx = {
    "mesh": {...},
    "inlet": {...},
    "chemistry": {...},
}

slv, species = build_reactor_from_context(case_ctx)
slv.initializeCase()
outlet = slv.steadyState(
    maxiter=500,
    relaxationFactorSpecie=0.4,
    relaxationFactorTemperature=0.4,
    convergenceCriteria=1e-6,
    temperatureClipLow=200.0,
    temperatureClipHigh=2000.0,
)

plotter = ReactorPlotter(slv)
plotter.save_temperature("temperature.png")
plotter.save_species("species.png")
plotter.save_all("profiles.png")

print(outlet.asDict(species=species))
```

## Practical notes

Several implementation details are especially relevant when using or extending this solver: [file:63]

- The current design supports exactly one reaction. `build_reactor_from_context()` raises an error if more than one reaction is present. [file:63]
- Reaction activity and external heating are controlled cellwise through the mesh zone flags. [file:63]
- Species are renormalized after each species sweep, which helps maintain physical mass fractions even under aggressive source terms. [file:63]
- The solver uses sparse matrices and `spsolve`, which is appropriate for the bidiagonal systems assembled here. [file:63]
- Temperature clipping is built into `steadyState()` and acts as a hard safeguard against runaway updates. [file:63]
- The outlet object stores both mass fractions and molar concentrations, making it convenient for optimization targets and postprocessing. [file:63]
