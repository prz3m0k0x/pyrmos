# main

## Overview

`main.py` is the orchestration layer of the project. It connects the YAML configuration files, the safe expression engine, the reactor solver, and the PSO optimizer into one executable optimization workflow.  

The script does not implement the reactor physics or the PSO algorithm itself. Instead, it prepares context dictionaries, applies optimization variables to those contexts, resolves dependent expressions, runs one reactor case per particle evaluation, extracts objective values, and finally launches the optimizer.  

## Role in the project

The script imports `UserExpression` from `scripts.usrExpr`, `PSOConfig` and `PSOOptimizer` from `scripts.PSOOPtimizer`, and `build_reactor_from_context`, `Outlet`, and `ReactorPlotter` from `scripts.solverReactor`. This makes `main.py` the integration point where configuration, chemistry simulation, and optimization are tied together.  

In practical terms, the script is responsible for: 

- Loading YAML-based reactor and PSO settings.
- Building a runtime context dictionary.
- Loading named user expressions and outlet expressions.
- Reapplying optimization variables to a fresh unresolved context for every particle.
- Resolving dependent expressions after particle application. 
- Running a reactor simulation for each optimization sample. 
- Evaluating derived output expressions and converting them to PSO objective values.  
- Creating case folders and writing YAML/PNG diagnostics.  

## Imported dependencies

The script uses these external and internal modules:  

| Import | Purpose |
| --- | --- |
| `copy` | Deep-copy of case contexts before expression resolution and particle application.   |
| `pathlib.Path` | File and directory handling.   |
| `yaml` | YAML loading and writing.   |
| `scripts.usrExpr.UserExpression` | Safe evaluation of symbolic expressions from configuration files.   |
| `scripts.PSOOPtimizer.PSOConfig` | PSO configuration container.   |
| `scripts.PSOOPtimizer.PSOOptimizer` | Optimization engine.   |
| `scripts.solverReactor.build_reactor_from_context` | Builds a solver instance from the resolved case context.   |
| `scripts.solverReactor.Outlet` | Fallback outlet extraction if the solver does not return one directly.   |
| `scripts.solverReactor.ReactorPlotter` | Writes reactor plots into the case directory.   |

## High-level workflow

The top-level workflow implemented in `main()` is:  

1. Build an unresolved base context from the `config` directory.  
2. Resolve a temporary copy of that context to safely construct the PSO configuration.  
3. Create an objective function that, for each particle, deep-copies the raw base context, applies particle values, resolves expressions, runs the reactor case, evaluates outputs, and returns objective values.  
4. Create a `PSOOptimizer` instance from random initialization.  
5. Run the optimizer and print the best particle at the end.  

This separation between unresolved template context and per-particle resolved context is a key design feature of the current script. It prevents dependent expressions from being permanently collapsed before optimization variables are updated.  

## Utility functions

### `load_yaml(path)`

This helper opens a YAML file, parses it using `yaml.safe_load`, and returns an empty dictionary if the file is empty. It is used as the base loader for all configuration blocks.  

### `_save_yaml(path, data)`

This helper writes a Python object to YAML. It ensures the parent directory exists and then uses `yaml.safe_dump(..., sort_keys=False)` to preserve a human-friendly key order in saved case files.  

This function is used to save generated case data such as `caseSetup.yaml`, `outlet.yaml`, and `reactorDebug.yaml`.  

### `load_expression_registry(path, expr_cls=UserExpression)`

This function loads a YAML mapping of expression names to expression strings and converts each string into a `UserExpression` object. It validates that the loaded YAML content is a dictionary and that every expression value is a string.  

The return value is a dictionary of the form:  

```python
{
    "expr1": UserExpression("..."),
    "expr2": UserExpression("..."),
}
```

It is used for both `userExpressions.yaml` and `outletConfig.yaml`.  

### `make_serializable_context(ctx)`

This helper deep-copies a context dictionary and removes the `expressions` and `outletExpressions` registries. These objects contain `UserExpression` instances, so removing them produces a cleaner YAML-serializable snapshot of the case.  

The function is used before saving `caseSetup.yaml` and before building `result_ctx` for output extraction.  

### `make_case_dir(base_dir, study_name)`

This function creates numbered case folders such as `cases/<study_name>1`, `cases/<study_name>2`, and so on. It increments the suffix until it finds a folder name that does not already exist.  

This gives each particle evaluation its own unique output directory.  

## Context construction

### `build_context(config_dir, expr_cls=UserExpression)`

This function loads the full project configuration from a directory. It reads:  

- `meshConfig.yaml`  
- `inletConfig.yaml`  
- `solverNumerics.yaml`  
- `speciesConfig.yaml`  
- `psoAlgorithm.yaml`  
- `userExpressions.yaml`  
- optionally `outletConfig.yaml`  

It then reorganizes these into a normalized runtime context dictionary with top-level keys:  

- `mesh`  
- `inlet`  
- `solver`  
- `chemistry`  
- `pso`  
- `expressions`  
- `outletExpressions`  

The inlet block is reshaped so that `diameter`, `velocity`, `temperature`, and `specie` are directly accessible in the runtime format expected by the solver builder and expression logic.  

## Expression resolution

### `resolve_value(value, root_ctx, expr_registry)`

This is the recursive worker used to replace symbolic expression names embedded inside the context with their numeric values. If `value` is a string that matches a key in the expression registry, it is evaluated immediately using the root context.  

If `value` is a dictionary or list, the function recursively resolves its children. It explicitly skips the `expressions` and `outletExpressions` keys to avoid trying to overwrite the expression registries themselves.  

### `resolve_expressions_in_context(ctx, max_passes=10)`

This function repeatedly calls `resolve_value()` until no more replacements occur or `max_passes` is exceeded. It exists because one expression may depend on another expression that becomes resolvable only after an earlier substitution pass.  

If the loop never stabilizes within the pass limit, the function raises `RuntimeError`. On success, it returns the mutated context with symbolic references replaced by numbers.  

### `evaluate_named_expressions(expr_registry, context)`

This function is different from `resolve_expressions_in_context()`. Instead of modifying the main context tree, it evaluates a registry of named expressions into a separate result dictionary.  

It repeatedly attempts to evaluate pending expressions, augmenting the available evaluation context with already-resolved results. This allows derived outputs to depend on each other as long as the dependency chain eventually resolves.  

If any expressions remain unresolved after the iteration budget, the function raises `RuntimeError`.  

## Dotted-path utilities

### `set_by_dotted_path(data, path, value)`

This helper assigns a value into a nested dictionary using a dotted key path such as `inlet.specie.so2`. It walks through the dictionary hierarchy and raises `KeyError` if an intermediate or final key does not exist.  

### `get_by_dotted_path(data, path)`

This helper reads a value from a nested dictionary using a dotted key path. It validates that each traversal step remains inside a dictionary and raises `KeyError` on invalid paths.  

These two helpers are central to parameter injection and output extraction.  

## Optimization-variable application

### `apply_particle_to_context(ctx, particle, parameter_defs)`

This function applies one PSO particle vector to the case context. It loops over the particle coordinates and the corresponding parameter definitions and writes each coordinate into the nested context using the parameter's `key` field.  

For example, if a parameter definition contains `key: inlet.specie.so2`, the corresponding particle value is written directly into `ctx["inlet"]["specie"]["so2"]`.  

## Objective extraction

### `extract_objectives_for_pso(result_ctx, output_defs)`

This function converts the solved reactor result into the response vector returned to the optimizer. It supports several output modes:  

- Values already present in `result_ctx["derived"]`.  
- Special-case outputs `expr3` and `expr4`.  
- General dotted-path extraction from the result context.  

The special cases are:  

- `expr3`: computes SO2 conversion as `(inlet_so2 - outlet_so2) / inlet_so2`.  
- `expr4`: computes temperature drop as `inlet_temperature - outlet_temperature`.  

Because the optimizer is written in minimization form, any output whose goal is `
maximize` is negated before being returned. Outputs marked as `minimize` are passed through unchanged.  

## Running one reactor case

### `run_case(case_ctx, case_dir)`

This function performs one full reactor simulation for one resolved case context. It builds the solver, initializes the case, reads numerical settings from the `solver` block, and calls `steadyState()` with the configured iteration limits, under-relaxation factors, residual criterion, and temperature clipping bounds.  

If the solver returns `None`, the function falls back to `Outlet.fromSolver(slv)`. After solving, it converts the outlet object to a labeled dictionary using `asDict(species=species)`.  

The function also writes a debug YAML file called `reactorDebug.yaml` containing mesh metadata and the outlet data. It then saves three plots into the case directory using `ReactorPlotter`:  

- `temperature.png`  
- `species.png`  
- `profiles.png`  

The return value is a compact outlet dictionary with temperature, species fractions, density, velocity, mass flow rate, and concentrations.  

## PSO configuration builder

### `make_pso_config_from_context(ctx)`

This function translates the `pso` block of the runtime context into a `PSOConfig` object. It extracts:  

- Global PSO settings from `ctx["pso"]["pso"]`.  
- Variable definitions from `ctx["pso"]["parameters"]`.  
- Response definitions from `ctx["pso"]["outputs"]`.  
- Optional linear constraints from `ctx["pso"].get("constraints", {})`.  

It builds lower and upper parameter bounds from each parameter definition and converts linear constraints into `constr_matrix`, `constr_lb`, and `constr_ub` arrays expected by `PSOConfig`.  

The returned configuration object contains swarm size and schedule settings such as `h_factor`, `max_iter`, `t_neighbors`, `w_init`, `w_finish`, `c1_init`, `c1_finish`, `c2_init`, `c2_finish`, and `v_max_factor`, along with bounds and constraints.  

## `main()` function

### Purpose

`main()` is the executable entry point of the script. It constructs the unresolved template context, derives the PSO configuration, defines the objective function closure, runs the optimizer, and prints the best result.  

### Detailed flow

#### 1. Build unresolved base context

The script starts with:  

```python
base_ctx_raw = build_context("config", UserExpression)
```

This preserves symbolic expressions in the template context rather than resolving them immediately.  

#### 2. Build a temporary resolved context for PSO settings

Next, `main()` creates a temporary resolved copy:  

```python
pso_ctx = resolve_expressions_in_context(copy.deepcopy(base_ctx_raw))
pso_cfg = make_pso_config_from_context(pso_ctx)
```

This allows PSO-related expressions to be evaluated without mutating the raw base template that will later be reused for each particle.  

#### 3. Define `objective_function(particle, iteration)`

The nested `objective_function` is the core callback passed to the optimizer. For every particle evaluation, it performs the following sequence:  

1. Deep-copy `base_ctx_raw` so that the new case starts from unresolved symbolic definitions.  
2. Apply the particle values into the copied context with `apply_particle_to_context(...)`.  
3. Print diagnostic information about the particle and selected inlet variables.  
4. Resolve all dependent expressions in the copied case context.  
5. Create a new numbered case directory under `cases/`.  
6. Save the resolved case setup to `caseSetup.yaml`.  
7. Run the reactor simulation through `run_case(...)`.  
8. Build `result_ctx` by combining the case input context and the reactor outlet data.  
9. Evaluate outlet expressions into `derived`.  
10. Save `outlet.yaml` containing both the raw outlet and the derived values.  
11. Convert the result to a PSO objective vector using `extract_objectives_for_pso(...)`.  

This design is especially important for dependent inlet expressions such as oxygen or nitrogen fractions defined in terms of another optimized variable. Because the context is copied first and resolved only after particle application, the dependent values are recalculated for every particle rather than frozen from startup.  

#### 4. Create and run the optimizer

After defining the objective function, `main()` creates the optimizer using:  

```python
optimizer = PSOOptimizer.from_random(pso_cfg, objective_function=objective_function)
swarm = optimizer.run()
```

After the optimization finishes, it reads and prints `swarm.global_best_position`.  

## Case output files

For each particle evaluation, the script writes a dedicated case directory under `cases/<study_name>N`. Inside that directory, the current script writes at least these artifacts:  

| File | Contents |
| --- | --- |
| `caseSetup.yaml` | Resolved input case used for the reactor run.   |
| `outlet.yaml` | Reactor outlet data plus evaluated derived outputs.   |
| `reactorDebug.yaml` | Mesh metadata and labeled outlet information.   |
| `temperature.png` | Axial temperature profile.   |
| `species.png` | Axial species mass-fraction profiles.   |
| `profiles.png` | Combined summary plot.   |

This output structure makes the optimization trace inspectable and reproducible after the run.  

## Data structures used by the script

The script relies on a nested context dictionary that acts as the shared data model between configuration loading, expression resolution, solver setup, and output processing.  

A simplified structure is:  

```python
ctx = {
    "mesh": {...},
    "inlet": {
        "diameter": ...,
        "velocity": ...,
        "temperature": ...,
        "specie": {...},
    },
    "solver": {...},
    "chemistry": {...},
    "pso": {...},
    "expressions": {...},
    "outletExpressions": {...},
}
```

The same structure is reused throughout the run, with particle values injected into selected locations and outlet data added later in a separate `result_ctx`.  

## Important design decision

One of the most important features of this version of `main.py` is that it keeps `base_ctx_raw` unresolved and only resolves expressions after particle values are applied inside the objective function. This prevents dependent symbolic inputs from being evaluated too early and becoming fixed before optimization begins.  

That design directly supports parameterizations where one optimized variable controls other inlet quantities through expressions.  

## Minimal execution example

```python
if __name__ == "__main__":
    main()
```

When executed, the script expects a `config/` directory containing the YAML files required by `build_context()`. It then launches the full reactor optimization workflow automatically.  

## Function summary

| Function | Purpose |
| --- | --- |
| `load_yaml` | Read one YAML file.   |
| `_save_yaml` | Write one YAML file and create parent folders if needed.   |
| `load_expression_registry` | Load named expressions as `UserExpression` objects.   |
| `make_serializable_context` | Remove non-serializable expression registries from a context snapshot.   |
| `make_case_dir` | Create a unique numbered case folder.   |
| `build_context` | Build the main runtime context from YAML configuration files.   |
| `resolve_value` | Recursive helper for symbolic-value replacement.   |
| `resolve_expressions_in_context` | Resolve embedded expressions inside the context.   |
| `evaluate_named_expressions` | Evaluate named derived expressions into a separate result dictionary.   |
| `set_by_dotted_path` | Write nested dictionary values via dotted keys.   |
| `get_by_dotted_path` | Read nested dictionary values via dotted keys.   |
| `apply_particle_to_context` | Inject one particle vector into the case context.   |
| `extract_objectives_for_pso` | Convert case results into PSO objective values.   |
| `run_case` | Execute one resolved reactor case and save diagnostics.   |
| `make_pso_config_from_context` | Convert context settings into `PSOConfig`.   |
| `main` | Launch the full optimization workflow.   |
