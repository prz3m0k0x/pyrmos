import yaml
from pathlib import Path
import copy
from scripts import PSOOPtimizer
from scripts import solverReactor
from scripts.usrExpr import UserExpression

def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f) or {}


def save_yaml(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)

def load_expression_registry(path, UserExpression):
    raw = load_yaml(path)
    return {name: UserExpression(expr) for name, expr in raw.items()}

def make_serializable_context(ctx):
    out = copy.deepcopy(ctx)
    out.pop("expressions", None)
    return out

def make_case_dir(base_dir, study_name):
    base_dir = Path(base_dir)
    base_dir.mkdir(parents=True, exist_ok=True)

    i = 1
    while True:
        case_dir = base_dir / f"{study_name}{i}"
        if not case_dir.exists():
            case_dir.mkdir(parents=True, exist_ok=False)
            return case_dir
        i += 1

def build_outlet_context(outlet_temperature, outlet_species):
    return {
        "temperature": outlet_temperature,
        "specie": outlet_species,
    }
        
def load_expression_registry(path, UserExpression):
    raw = load_yaml(path)
    return {name: UserExpression(expr) for name, expr in raw.items()}

def build_context(config_dir, UserExpression):
    config_dir = Path(config_dir)

    mesh_cfg = load_yaml(config_dir / "meshConfig.yaml")
    inlet_cfg = load_yaml(config_dir / "inletConfig.yaml")
    solver_cfg = load_yaml(config_dir / "solverNumerics.yaml")
    species_cfg = load_yaml(config_dir / "speciesConfig.yaml")
    pso_cfg = load_yaml(config_dir / "psoAlgorithm.yaml")
    expr_registry = load_expression_registry(config_dir / "userExpressions.yaml", UserExpression)

    ctx = {
        "mesh": mesh_cfg,
        "inlet": {
            "diameter": inlet_cfg["diameter"],
            "velocity": inlet_cfg["inletVelocity"],
            "temperature": inlet_cfg["inletTemperature"],
            "specie": inlet_cfg["speciesYInlet"],
        },
        "solver": solver_cfg,
        "chemistry": species_cfg,
        "pso": pso_cfg,
        "expressions": expr_registry,
    }

    ctx["zones"] = ctx["mesh"]["zones"]
    ctx["species"] = ctx["chemistry"]["species"]
    ctx["reactions"] = ctx["chemistry"]["reactions"]

    return ctx

def resolve_value(value, root_ctx, expr_registry):
    if isinstance(value, str) and value in expr_registry:
        return expr_registry[value].value(root_ctx), True

    if isinstance(value, dict):
        changed_any = False
        for key, subval in list(value.items()):
            if key == "expressions":
                continue
            newval, changed = resolve_value(subval, root_ctx, expr_registry)
            value[key] = newval
            changed_any = changed_any or changed
        return value, changed_any

    if isinstance(value, list):
        changed_any = False
        for i, item in enumerate(value):
            newval, changed = resolve_value(item, root_ctx, expr_registry)
            value[i] = newval
            changed_any = changed_any or changed
        return value, changed_any

    return value, False

def resolve_expressions_in_context(ctx, max_passes=10):
    expr_registry = ctx["expressions"]

    for _ in range(max_passes):
        _, changed = resolve_value(ctx, ctx, expr_registry)
        if not changed:
            break
    else:
        raise RuntimeError("Expression resolution exceeded max_passes")

    return ctx

def find_unresolved_expr_names(value, expr_registry, found=None):
    if found is None:
        found = []

    if isinstance(value, str) and value in expr_registry:
        found.append(value)
    elif isinstance(value, dict):
        for key, subval in value.items():
            if key == "expressions":
                continue
            find_unresolved_expr_names(subval, expr_registry, found)
    elif isinstance(value, list):
        for item in value:
            find_unresolved_expr_names(item, expr_registry, found)

    return found
def set_by_dotted_path(data, path, value):
    parts = path.split(".")
    node = data

    for part in parts[:-1]:
        if part not in node:
            raise KeyError(f"Unknown path while setting value: {path}")
        node = node[part]

    last = parts[-1]
    if last not in node:
        raise KeyError(f"Unknown final key while setting value: {path}")

    node[last] = value

def get_by_dotted_path(data, path):
    parts = path.split(".")
    node = data

    for part in parts:
        if not isinstance(node, dict):
            raise KeyError(f"Cannot descend into non-dict at: {part} for path {path}")
        if part not in node:
            raise KeyError(f"Unknown path: {path}")
        node = node[part]

    return node

def apply_particle_to_context(ctx, particle, parameter_defs):
    for x, param in zip(particle, parameter_defs):
        set_by_dotted_path(ctx, param["key"], float(x))

def evaluate_named_expressions(expr_registry, context):
    results = {}
    pending = dict(expr_registry)

    for _ in range(20):
        changed = False

        for name in list(pending.keys()):
            expr = pending[name]
            extended_context = dict(context)
            extended_context.update(results)

            try:
                results[name] = expr.value(extended_context)
                del pending[name]
                changed = True
            except (KeyError, ValueError, TypeError):
                pass

        if not changed:
            break

    if pending:
        raise RuntimeError(
            f"Could not resolve outlet expressions: {list(pending.keys())}"
        )

    return results

def run_single_case(base_ctx, particle, UserExpression, outlet_config_path, runs_root="cases"):
    case_ctx = copy.deepcopy(base_ctx)

    apply_particle_to_context(case_ctx, particle, case_ctx["pso"]["parameters"])
    case_ctx = resolve_expressions_in_context(case_ctx)

    case_dir = make_case_dir(runs_root, case_ctx["pso"]["study"]["name"])

    resolved_case = make_serializable_context(case_ctx)
    save_yaml(case_dir / "caseSetup.yaml", resolved_case)

    outlet_ctx = run_reactor_case(case_ctx, case_dir)

    full_context = copy.deepcopy(resolved_case)
    full_context["outlet"] = outlet_ctx

    outlet_expr_registry = load_expression_registry(outlet_config_path, UserExpression)
    outlet_expr_values = evaluate_named_expressions(outlet_expr_registry, full_context)

    outlet_data = {
        "outlet": outlet_ctx,
        "derived": outlet_expr_values,
    }
    save_yaml(case_dir / "outlet.yaml", outlet_data)

    result_ctx = copy.deepcopy(full_context)
    result_ctx["derived"] = outlet_expr_values

    return case_dir, result_ctx

def extract_objectives(result_ctx, output_defs):
    values = []

    for outdef in output_defs:
        key = outdef["key"]

        if key in result_ctx.get("derived", {}):
            val = result_ctx["derived"][key]
        else:
            val = get_by_dotted_path(result_ctx, key)

        values.append(float(val))

    return values
def evaluate_particle(base_ctx, particle, UserExpression, outlet_config_path, runs_root="cases"):
    case_dir, result_ctx = run_single_case(
        base_ctx=base_ctx,
        particle=particle,
        UserExpression=UserExpression,
        outlet_config_path=outlet_config_path,
        runs_root=runs_root,
    )

    output_defs = base_ctx["pso"]["outputs"]
    objectives = extract_objectives(result_ctx, output_defs)

    return objectives

def postprocess_case(
    ctx,
    UserExpression,
    outlet_temperature,
    outlet_species,
    outlet_config_path,
    runs_root="runs",
):
    study_name = ctx["pso"]["study"]["name"]
    case_dir = make_case_dir(runs_root, study_name)

    resolved_case = make_serializable_context(ctx)
    save_yaml(case_dir / "caseSetup.yaml", resolved_case)

    outlet_ctx = build_outlet_context(outlet_temperature, outlet_species)

    full_context = copy.deepcopy(resolved_case)
    full_context["outlet"] = outlet_ctx

    outlet_expr_registry = load_expression_registry(outlet_config_path, UserExpression)
    outlet_expr_values = evaluate_named_expressions(outlet_expr_registry, full_context)

    outlet_data = {
        "outlet": outlet_ctx,
        "derived": outlet_expr_values,
    }

    save_yaml(case_dir / "outlet.yaml", outlet_data)

    return case_dir, outlet_data

def main():
    base_ctx = build_context("config", UserExpression)
    base_ctx = resolve_expressions_in_context(base_ctx)

    pso_cfg = base_ctx["pso"]
    parameter_defs = pso_cfg["parameters"]

    bounds = [tuple(p["bounds"]) for p in parameter_defs]

    def objective_function(particle):
        return evaluate_particle(
            base_ctx=base_ctx,
            particle=particle,
            UserExpression=UserExpression,
            outlet_config_path="config/outletConfig.yaml",
            runs_root="cases",
        )

    optimizer = PSOOptimizer(
        bounds=bounds,
        config=pso_cfg["pso"],
        objective_function=objective_function,
    )

    best = optimizer.optimize()
    print(best)


ctx = build_context("config", UserExpression)
ctx = resolve_expressions_in_context(ctx)

outlet_temperature = 720.0
outlet_species = {
    "so2": 0.03,
    "so3": 0.07,
    "o2": 0.18,
    "n2": 0.72,
}

case_dir, outlet_data = postprocess_case(
    ctx=ctx,
    UserExpression=UserExpression,
    outlet_temperature=outlet_temperature,
    outlet_species=outlet_species,
    outlet_config_path="config/outletConfig.yaml",
    runs_root="cases",
)
