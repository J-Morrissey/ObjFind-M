# objfind_generalized.py
# Generalized ObjFind/CHO-style Pyomo formulation for arbitrary COBRA models

import warnings
import pandas as pd
import numpy as np
import cobra
from cobra.util.array import create_stoichiometric_matrix
from pyomo.environ import (
    ConcreteModel, Set, Var, NonNegativeReals, Binary, Objective, minimize, Constraint, SolverFactory
)
from objfindm_mapping import objfindm_mapping
from convert_to_irreversible import convert_to_irreversible

# =========================
# Config (edit as needed)
# =========================
CONFIG = {
    # --- Inputs ---
    # Works with SBML too; run_objfind will detect .xml and use cobra.io.read_sbml_model
    "model_path": "iCHO2441.xml",
    "exchange_file": "C13 Exchange Rates Quadmin.xlsx",   # rows: exchange rxn IDs, cols: experiments
    "intracellular_file": "Intracellular Reactions All.xlsx",
    "mapping_file": "iCHO2441 Mapping ObjFind-M.xlsx",
    "experiment": "ATTnegative",   # or "SV40positive", "SV40negative"

    # --- Objective hooks (OPTIONAL) ---
    # Add (v - measured)^2 terms for specific reactions, by ID, if you want
    "exometabolite_targets": {
        # "BIOMASS_CHO": <number>,
        # "PROD_IgG": <number>,
    },

    # --- Model overrides (OPTIONAL) ---
    "force_bounds": {
        "FECR": (0, 0),      # you clamp FECR in your script
        # If you want to mirror the other (currently commented) redox transport clamps, add them:
        # "NADPHtru": (0, 0),
        # "NADPHtxu": (0, 0),
        # "NADHtpu": (0, 0),
        # "NADHtru": (0, 0),
        # "NADPtru": (0, 0),
        # "NADPtxu": (0, 0),
        # "NADtm":   (0, 0),
        # "NADtn":   (0, 0),
        # "NADtru":  (0, 0),
        # "NADtx":   (0, 0),
    },

    # Set exchange reactions exactly to measured fluxes from the sheet (handles *_reverse if present)
    "fix_exchange_to_measurement": True,

    # Parsimonious cap (you used 100)
    "flux_sum_max": 100.0,

    # Solver (you can still set NEOS_EMAIL env var; this just picks the local solver name)
    "solver": "ipopt",

    # Output
    "output_excel": "ObjFind_iCHO2441_Results.xlsx",
}

# ====== Utilities ======
def rxn_index_map(model):
    """Reaction ID -> index (column in S)."""
    return {rxn.id: j for j, rxn in enumerate(model.reactions)}

def safe_get_idx(rxn_to_idx, rxn_id):
    """Return index or None if absent."""
    return rxn_to_idx.get(rxn_id, None)

def apply_forced_bounds(model, force_bounds: dict):
    present_ids = {r.id for r in model.reactions}
    for rid, (lb, ub) in force_bounds.items():
        if rid not in present_ids:
            warnings.warn(f"Forced bound skipped (rxn not found): {rid}")
            continue
        rxn = model.reactions.get_by_id(rid)
        rxn.bounds = (float(lb), float(ub))

def fix_exchanges_from_sheet(model, exchange_df, experiment):
    """Sets exchange reaction(s) to the measured flux. Supports split reversible naming (<id> and <id>_reverse)."""
    present_ids = {r.id for r in model.reactions}
    for rid in exchange_df.index:
        if experiment not in exchange_df.columns:
            raise ValueError(f"Experiment '{experiment}' not found in exchange data columns: {list(exchange_df.columns)}")
        val = float(exchange_df.loc[rid, experiment])
        # forward (export) is non-negative in irreversible models
        if rid in present_ids:
            model.reactions.get_by_id(rid).bounds = (max(0.0, val), max(0.0, val))
        # backward (import) often stored as <id>_reverse
        rev = f"{rid}_reverse"
        if rev in present_ids:
            model.reactions.get_by_id(rev).bounds = (max(0.0, -val), max(0.0, -val))


# ====== Main runner ======
def run_objfind(config: dict):
    # --- Load inputs ---
    GEM = cobra.io.load_json_model(config["model_path"]) if config["model_path"].lower().endswith(".json") \
        else cobra.io.read_sbml_model(config["model_path"])

    exchange_df = pd.read_excel(config["exchange_file"], index_col=0)
    intracellular_df = pd.read_excel(config["intracellular_file"], index_col=0)
    mapping_table = pd.read_excel(config["mapping_file"], index_col=0)
    experiment = config["experiment"]
    flux_sum_max = float(config["flux_sum_max"])

    # --- Normalize to irreversible model (network-agnostic) ---
    convert_to_irreversible(GEM)

    # --- Build stoichiometric matrix & index map ---
    S = create_stoichiometric_matrix(GEM, array_type='dense', dtype=float)
    m, n = S.shape
    rid_to_idx = rxn_index_map(GEM)

    # --- Apply exchange & forced bounds ---
    if config.get("fix_exchange_to_measurement", True):
        fix_exchanges_from_sheet(GEM, exchange_df, experiment)
    if config.get("force_bounds"):
        apply_forced_bounds(GEM, config["force_bounds"])

    # --- Build mapping artifacts (network-agnostic; plug in your mapper) ---
    mapping_art = objfindm_mapping(
        GEM=GEM,
        mapping_table=mapping_table,
        intracellular_data=intracellular_df,
        experiment=experiment
    )

    # Unpack with safe defaults
    lb = np.array(mapping_art.get("lowerbounds", [0.0]*n), dtype=float)
    ub = np.array(mapping_art.get("upperbounds", [1e6]*n), dtype=float)
    # clip to ≥ 0 because we model v ≥ 0 for irreversibles
    lb = np.maximum(lb, 0.0)
    ub = np.maximum(ub, 0.0)

    # Objective pieces
    ser_irrev_idx = mapping_art.get("series_irreversible_indices", [])
    ser_irrev_c13 = mapping_art.get("series_irreversible_c13", [])
    ser_irrev_coef = mapping_art.get("series_irreversible_obj_coeff", [])

    ser_rev_f_idx = mapping_art.get("series_reversible_forward_indices", [])
    ser_rev_b_idx = mapping_art.get("series_reversible_backward_indices", [])
    ser_rev_c13  = mapping_art.get("series_reversible_c13", [])
    ser_rev_coef = mapping_art.get("series_reversible_obj_coeff", [])

    parallel_reactions = mapping_art.get("parallel_reactions", [])  # list of dicts

    # --- Exometabolite fit targets (network agnostic) ---
    # E.g., {"Biomass": 0.05, "Product": 0.01}
    exo_targets = config.get("exometabolite_targets", {})
    exo_terms = []
    for r_id, measured in exo_targets.items():
        j = safe_get_idx(rid_to_idx, r_id)
        if j is None:
            warnings.warn(f"Exometabolite target skipped; reaction not in model: {r_id}")
            continue
        exo_terms.append((j, float(measured)))

    # --- Pyomo model ---
    model = ConcreteModel()
    model.i = Set(initialize=range(m))
    model.j = Set(initialize=range(n))

    # Decision vars
    model.v = Var(model.j, within=NonNegativeReals)   # irreversible fluxes ≥ 0
    model.c = Var(model.j)                            # Lagrange stationarity var
    model.u1 = Var(model.i)                           # duals for mass balance
    model.u2 = Var(model.j, within=NonNegativeReals)  # duals for upper bound
    model.u3 = Var(model.j, within=NonNegativeReals)  # duals for lower bound
    model.u4 = Var(within=NonNegativeReals)           # parsimonious scalar

    # (Binary & aux vars kept for parity with your template; not all are used below)
    model.delta  = Var(model.j, within=Binary)
    model.c_abs  = Var(model.j, within=NonNegativeReals)
    model.z      = Var(model.j, within=Binary)

    model.series_irrev_sum = Var()
    model.series_rev_sum = Var()
    model.parallel_sum = Var()
    model.exometabolomics_sum = Var()

    # Quadratic pieces
    if ser_irrev_idx:
        model.series_irreversible = Constraint(
            expr=model.series_irrev_sum == sum(
                coef * (model.v[j] - c13)**2
                for coef, j, c13 in zip(ser_irrev_coef, ser_irrev_idx, ser_irrev_c13)
            )
        )
    else:
        model.series_irreversible = Constraint(expr=model.series_irrev_sum == 0)

    if ser_rev_f_idx:
        model.series_reversible = Constraint(
            expr=model.series_rev_sum == sum(
                coef * (model.v[jf] - model.v[jb] - c13)**2
                for coef, jf, jb, c13 in zip(ser_rev_coef, ser_rev_f_idx, ser_rev_b_idx, ser_rev_c13)
            )
        )
    else:
        model.series_reversible = Constraint(expr=model.series_rev_sum == 0)

    if parallel_reactions:
        model.parallel = Constraint(
            expr=model.parallel_sum == sum(
                pr['MILP Coefficent'] * (
                    sum(model.v[idx] / dirn for idx, dirn in zip(pr['Model Indices'], pr['Model Directionality']))
                    - pr['MEDIAN']
                )**2
                for pr in parallel_reactions
            )
        )
    else:
        model.parallel = Constraint(expr=model.parallel_sum == 0)

    # Exometabolite fit terms (by reaction ID)
    if exo_terms:
        model.exometabolomics = Constraint(
            expr=model.exometabolomics_sum == sum(
                (model.v[j] - meas)**2 for j, meas in exo_terms
            )
        )
    else:
        model.exometabolomics = Constraint(expr=model.exometabolomics_sum == 0)

    # Objective
    model.objective = Objective(
        expr = model.series_irrev_sum + model.series_rev_sum + model.parallel_sum + model.exometabolomics_sum,
        sense = minimize
    )

    # Mass balance S v = 0
    def mass_balance_rule(model, i):
        return sum(S[i, j] * model.v[j] for j in model.j) == 0.0
    model.mass_balance = Constraint(model.i, rule=mass_balance_rule)

    # Bounds via dualized form
    def flux_lb(model, j): return model.v[j] >= lb[j]
    def flux_ub(model, j): return model.v[j] <= ub[j]
    model.flux_lb = Constraint(model.j, rule=flux_lb)
    model.flux_ub = Constraint(model.j, rule=flux_ub)

    # c in [0,1] (if needed for your dual relations)
    model.c_min = Constraint(model.j, rule=lambda m, j: m.c[j] >= 0)
    model.c_max = Constraint(model.j, rule=lambda m, j: m.c[j] <= 1)

    # Stationary Lagrangian: c[j] = u1^T*Scol + u2[j] - u3[j] + u4
    def stationary_lagrangian(model, j):
        return model.c[j] == sum(model.u1[i] * S[i, j] for i in model.i) + model.u2[j] - model.u3[j] + model.u4
    model.stationary_lagrangian = Constraint(model.j, rule=stationary_lagrangian)

    # Parsimonious cap: sum v <= flux_sum_max
    model.sum_flux = Constraint(expr=sum(model.v[j] for j in model.j) <= flux_sum_max)

    # Dual objective identity (if you’re enforcing strong duality)
    model.dual_obj = Constraint(
        expr=sum(model.c[j] * model.v[j] for j in model.j)
            == flux_sum_max * model.u4 + sum(model.u2[j] * ub[j] - model.u3[j] * lb[j] for j in model.j)
    )

    # Normalization of c (as in your script)
    model.sum_abs_c = Constraint(expr=sum(model.c[j] for j in model.j) == 1)

    # --- Solve ---
    solver = config.get("solver", "ipopt")
    opt = SolverFactory(solver)
    results = opt.solve(model)
    # Optional: check results.solver.termination_condition

    # --- Collect results ---
    rxn_list = [r.id for r in GEM.reactions]
    fluxes = [float(model.v[j]()) for j in range(n)]
    coeffs = [float(model.c[j]()) for j in range(n)]

    df_flux = pd.DataFrame({"Reaction": rxn_list, f"{experiment} Flux": fluxes})
    df_coef = pd.DataFrame({"Reaction": rxn_list, f"{experiment} Coefficient": coeffs})


    # Save
    with pd.ExcelWriter(config["output_excel"], engine="xlsxwriter") as writer:
        df_coef.to_excel(writer, sheet_name="Coefficients", index=False)
        df_flux.to_excel(writer, sheet_name="Fluxes", index=False)

    print(f"Solved. Objective = {float(model.objective())}")
    # Show top coefficients
    top_n = 10
    print("\nTop coefficients:")
    print(
        df_coef.sort_values(by=f"{experiment} Coefficient", ascending=False)
            .head(top_n)
            .to_string(index=False)
    )
    print(f"Results written to: {config['output_excel']}")
    return df_coef, df_flux, results


# ====== Run if executed directly ======
if __name__ == "__main__":
    run_objfind(CONFIG)
