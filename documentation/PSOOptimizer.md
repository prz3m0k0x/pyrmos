# PSOOptimizer

GitHub-flavored Markdown documentation for the `PSOOPtimizer.py` module. It explains the purpose of the module, the role of each class, the expected data shapes, and the normal usage pattern. [file:60]

## Overview

The module implements Particle Swarm Optimization with support for both single-objective and multi-objective problems. Its main components are `PSOConfig`, `Swarm`, `HistoryLogger`, and `PSOOptimizer`. [file:60]

The optimizer evaluates particles through a user-provided objective function, updates velocities and positions, tracks personal bests, and maintains neighborhood-based best solutions. In the multi-objective case, it uses Das-Dennis weight vectors together with weighted Tchebycheff scalarization. [file:60]

## Main classes

### `PSOConfig`

`PSOConfig` is a dataclass that stores algorithm settings, parameter bounds, and optional linear constraints. It contains fields such as `h_factor`, `max_iter`, `n_params`, `n_responses`, `t_neighbors`, the PSO coefficients (`w`, `c1`, `c2`), and constraint definitions. [file:60]

During `__post_init__`, bounds and constraint arrays are converted to NumPy arrays and validated. Missing bounds are replaced with `-inf` and `inf`, and constraint dimensions are checked for consistency with `n_params`. [file:60]

#### Population size

The `pop_size()` method returns the number of particles:

- For `n_responses == 1`, it returns `h_factor`. [file:60]
- For `n_responses > 1`, it returns `comb(h_factor + n_responses - 1, n_responses - 1)`. [file:60]

This ties the swarm size to the number of Das-Dennis weight vectors used in the decomposition. [file:60]

#### Das-Dennis weights

The static method `das_dennis_weights(m, H)` generates normalized Das-Dennis weight vectors by recursively enumerating integer partitions and dividing by `H`. These weights are later used for scalarizing multi-objective responses. [file:60]

#### Velocity schedule

The `algorithm_velocity_parameters` property returns a linearly interpolated array of shape `(max_iter, 3)` between the initial and final values of `w`, `c1`, and `c2`. This allows the optimizer coefficients to evolve over the run. [file:60]

### `Swarm`

`Swarm` stores the dynamic PSO state. It is initialized from a `PSOConfig`, an array of initial particles, and the corresponding objective values. [file:60]

It keeps track of:

- Current particle positions. [file:60]
- Current velocities. [file:60]
- Das-Dennis decomposition weights. [file:60]
- Neighborhood topology. [file:60]
- Reference and nadir points. [file:60]
- Personal best positions and gains. [file:60]
- Neighborhood best positions and gains. [file:60]

#### Initialization

Particle positions are stored in an array of shape `(n_params, pop_size)`. The maximum velocity is defined as `v_max_factor * (x_ub - x_lb)` for each parameter, and initial velocities are sampled uniformly in the interval `[-v_max, v_max]`. [file:60]

#### Neighborhood topology

The swarm builds Das-Dennis weights, computes pairwise Euclidean distances between them, and assigns each weight vector a neighborhood of the `t_neighbors` nearest weights. This creates a decomposition-based local-best topology. [file:60]

#### Scalarization

The swarm uses:

- `z_ref`: componentwise minimum of gains. [file:60]
- `z_nad`: componentwise maximum of gains. [file:60]

The `_tcheby_scalars()` method normalizes deviations from `z_ref`, multiplies them by the decomposition weights, and takes the maximum over objectives. It then adds a penalty term for infeasible particles. [file:60]

#### Constraint penalty

The `_penalty()` method adds squared penalties for:

- Lower-bound violations. [file:60]
- Upper-bound violations. [file:60]
- Linear constraint violations of the form `constr_lb <= A x <= constr_ub`. [file:60]

These penalties are added directly to the scalarized objective value. [file:60]

#### Particle update

The `step(w, c1, c2)` method performs one PSO iteration:

1. Draw random matrices `r1` and `r2`. [file:60]
2. Update velocity using inertia, cognitive, and social terms. [file:60]
3. Clip velocity to `[-v_max, v_max]`. [file:60]
4. Update particle positions. [file:60]
5. Clip particle positions to `[x_lb, x_ub]`. [file:60]

The social term uses neighborhood-best positions stored in `gbest_positions`, not a single fully global best. [file:60]

#### Best updates

The `update_bests(new_gains)` method:

1. Updates `z_ref` and `z_nad`. [file:60]
2. Recomputes scalar values for the current particles. [file:60]
3. Replaces personal bests where the new scalar value is better. [file:60]
4. Rebuilds neighborhood-best data. [file:60]

#### Best-particle accessor

The `global_best_position` property returns the personal-best position associated with the minimum scalarized score in `pbest_scalars`. This gives a single representative best solution at the end of the run. [file:60]

### `HistoryLogger`

`HistoryLogger` records optimization history for later analysis. It stores:

- `best_gain_history`: best scalar value per iteration. [file:60]
- `particle_history`: full particle coordinates with shape `(n_params, pop_size, max_iter + 1)`. [file:60]

The `log()` method records the current particle state and best scalar value, and prints a progress line. The `save()` method writes:

- `particle_history.npy` [file:60]
- `gain_history.csv` [file:60]

### `PSOOptimizer`

`PSOOptimizer` is the top-level driver class. It connects the configuration, the objective function, the swarm state, and the logger. [file:60]

Its constructor takes:

- A `PSOConfig` instance. [file:60]
- An objective function with signature `Callable[[np.ndarray, int], np.ndarray]`. [file:60]
- An array of initial particles. [file:60]

During initialization, it evaluates the initial particles at iteration `0`, creates the swarm, creates the history logger, and logs the initial state. [file:60]

## Objective function interface

The objective function is called as:

```python
objective_function(particle, iteration)
```

Where:

- `particle` is a 1D NumPy array containing one particle. [file:60]
- `iteration` is the current iteration number. [file:60]

Return value rules:

- For a single-objective problem, return one scalar-like value or a one-element array. [file:60]
- For a multi-objective problem, return one value per response in a consistent order. [file:60]

The optimizer converts each returned value to `float` and stacks all particle responses column-wise. [file:60]

## Public workflow

### Random initialization

The `from_random()` classmethod generates uniformly random initial particles between `x_lb` and `x_ub`, then constructs a `PSOOptimizer` instance. [file:60]

### Evaluation

The private method `_evaluate()` loops over particles, calls the objective function particle by particle, converts each result to a NumPy array, and returns a gain matrix built with `np.column_stack`. [file:60]

### Run loop

The `run()` method:

1. Reads the iteration-dependent schedule from `algorithm_velocity_parameters`. [file:60]
2. Updates particle states with `swarm.step()`. [file:60]
3. Re-evaluates the swarm. [file:60]
4. Updates personal and neighborhood bests. [file:60]
5. Logs the current iteration. [file:60]
6. Returns the final `Swarm` object. [file:60]

## Data layout

The module uses a consistent array layout where parameters are the first axis and particles are the second axis. [file:60]

| Quantity | Shape | Description |
| --- | --- | --- |
| `particles` | `(n_params, pop_size)` | Current particle positions. [file:60] |
| `velocity` | `(n_params, pop_size)` | Current particle velocities. [file:60] |
| `initial_gains` | `(n_responses, pop_size)` | Initial objective values. [file:60] |
| `new_gains` | `(n_responses, pop_size)` | Objective values after an update. [file:60] |
| `weights` | `(pop_size, n_responses)` | Das-Dennis decomposition weights. [file:60] |
| `particle_history` | `(n_params, pop_size, max_iter + 1)` | Full saved particle history. [file:60] |

## Minimal example

```python
import numpy as np
from scripts.PSOOPtimizer import PSOConfig, PSOOptimizer

config = PSOConfig(
    h_factor=20,
    max_iter=30,
    n_params=2,
    n_responses=1,
    x_lb=[0.0, 300.0],
    x_ub=[1.0, 900.0],
)

def objective_function(particle: np.ndarray, iteration: int) -> np.ndarray:
    x1, x2 = particle
    value = (x1 - 0.3) ** 2 + ((x2 - 600.0) / 300.0) ** 2
    return np.array([value], dtype=float)

optimizer = PSOOptimizer.from_random(config, objective_function)
swarm = optimizer.run()

best_particle = swarm.global_best_position
print("Best particle:", best_particle)
```

## Output handling

The optimizer does not save history automatically. To export recorded data, call the logger explicitly: [file:60]

```python
from pathlib import Path

optimizer.logger.save(Path("results"))
```

This writes particle trajectories to `particle_history.npy` and best-gain history to `gain_history.csv`. [file:60]

## Practical notes

This implementation behaves like a decomposition-based local-best PSO with penalty handling for bounds and optional linear constraints. It is well suited for cases where the same code should support both scalar and vector objective functions. [file:60]

The `run()` method returns the final `Swarm` object rather than a scalar optimum. Downstream code should therefore read `swarm.global_best_position`, `swarm.pbest_positions`, or `swarm.pbest_gains` depending on what summary is needed. [file:60]