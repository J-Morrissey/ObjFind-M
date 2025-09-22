from typing import Dict, List, Any
import pandas as pd

def objfindm_mapping(
    GEM,
    mapping_table: pd.DataFrame,
    intracellular_data: pd.DataFrame,
    experiment: str
) -> Dict[str, Any]:
    """
    Build all mapping artifacts needed by ObjFindCHO from a mapping table
    and intracellular (13C) data for a given experiment.

    Returns a dict with:
      - reaction_ids: List[str]
      - next_flux_constrained_reactions: List[str]
      - mapping_table_out: pd.DataFrame (mapping_table with LB/UB/MEDIAN & objective coeffs added)
      - lowerbounds_dict, upperbounds_dict, medians_dict: Dict[str, float]
      - parallel_reactions: List[dict] (each dict has Model Names/Indices, Directionality, LB/UB/MEDIAN, etc.)
      - series_irreversible_* and series_reversible_* lists for Pyomo objective construction:
            series_irreversible_names, series_irreversible_indices, series_irreversible_c13, series_irreversible_obj_coeff
            series_reversible_forward_names, series_reversible_forward_indices,
            series_reversible_backward_names, series_reversible_backward_indices,
            series_reversible_c13, series_reversible_obj_coeff
      - lowerbounds, upperbounds: Lists[float] in GEM reaction order
    """
    # Basic lists from model
    reaction_ids = [rxn.id for rxn in GEM.reactions]

    # Constrained reactions from mapping (exclude 'none')
    next_flux_constrained_reactions = [x for x in list(mapping_table['Model Name']) if x != 'none']

    # Identify all 13C reaction names present in the mapping (skip NaNs)
    all_reaction_names = [x for x in list(mapping_table['Enzyme']) if str(x) != 'nan']

    # Columns present in intracellular_data: *_LB, *_UB, *_MEDIAN.
    # Build a set of base reaction names that are bounded in the data.
    bounded_c13reactions = list(intracellular_data.columns)
    for suffix in ['_LB', '_MEDIAN', '_UB']:
        bounded_c13reactions = [i.replace(suffix, '') for i in bounded_c13reactions]
    bounded_c13reactions = list(dict.fromkeys(bounded_c13reactions))  # unique

    # Indices to split mapping_table into per-13C-reaction blocks
    reaction_indexes = []
    for name in all_reaction_names:
        reaction_indexes.append(list(mapping_table['Enzyme']).index(name))
    reaction_indexes.append(len(list(mapping_table['Enzyme'])))

    # Output containers mirroring your original code
    lowerbounds: List[Any] = []
    upperbounds: List[Any] = []
    medians: List[Any] = []
    milp_objective_coefficients: List[Any] = []
    parallel_reactions: List[dict] = []

    # Loop through each 13C reaction and compute per-row mapping info
    for i, name in enumerate(all_reaction_names):
        # Segment of mapping rows corresponding to this 13C reaction
        seg = mapping_table.iloc[reaction_indexes[i]:reaction_indexes[i + 1]][:]

        # Determine how many series ('s') and independent parallel groups exist
        pathways = list(seg['Pathway'])
        number_series = pathways.count('s')
        without_s = [x for x in pathways if x != 's']
        if number_series == len(pathways):
            number_independent_parallel = 0
        else:
            number_independent_parallel = max(1, len([k for k in range(len(without_s)) if without_s[k] != without_s[k - 1]]))

        # Coefficient distributes weight across series and parallel groupings
        milp_coeff = 1.0 / float(number_series + number_independent_parallel) if (number_series + number_independent_parallel) > 0 else 0.0

        # Track parallel block starts so we can group them
        parallel_indexes = [x for x in range(len(pathways)) if pathways[x] != pathways[x - 1]]
        parallel_indexes.insert(0, 'placeholder')

        parallel_names: List[str] = []
        parallel_directionality: List[float] = []
        parallel_reaction: dict = {}

        # Row-wise processing
        for j in range(len(seg['Model Name'])):
            model_name = seg.iloc[j]['Model Name']
            pathway = seg.iloc[j]['Pathway']
            direction = seg.iloc[j]['Direction']

            # Default to NaN rows if missing in data or not bounded in intracellular dataset
            if (name not in bounded_c13reactions) or (str(pathway) == 'nan') or (str(intracellular_data.get(f'{name}_LB', {}).get(experiment, 'nan')) == 'nan'):
                lowerbounds.append('nan')
                upperbounds.append('nan')
                medians.append('nan')
                milp_objective_coefficients.append('nan')
                continue

            # Pull LB/UB/MEDIAN for the 13C reaction and adjust for mapping direction
            c13_lb = intracellular_data[f'{name}_LB'][experiment] * (1.0 / float(direction))
            c13_ub = intracellular_data[f'{name}_UB'][experiment] * (1.0 / float(direction))
            c13_median = intracellular_data[f'{name}_MEDIAN'][experiment] * (1.0 / float(direction))

            if pathway == 's':
                # Series: direct mapping with directional correction and min/max ordering
                lowerbounds.append(min(c13_lb, c13_ub))
                upperbounds.append(max(c13_lb, c13_ub))
                medians.append(c13_median)
                milp_objective_coefficients.append(milp_coeff)

            # Parallel mapping: accumulate a group and later add a quadratic constraint on their sum
            if isinstance(pathway, str) and ('p' in pathway):
                # Start a new parallel group if this row marks a new block
                if j in parallel_indexes:
                    if parallel_reaction:  # push the previous group if not empty
                        parallel_reactions.append(parallel_reaction)
                    parallel_names = []
                    parallel_directionality = []
                    parallel_reaction = {}

                parallel_names.append(model_name)
                parallel_directionality.append(direction)
                parallel_reaction['Parallel Order'] = pathway
                parallel_reaction['C13 Name'] = name
                parallel_reaction['Model Names'] = parallel_names
                parallel_reaction['Model Directionality'] = parallel_directionality
                parallel_reaction['LB'] = intracellular_data[f'{name}_LB'][experiment]
                parallel_reaction['UB'] = intracellular_data[f'{name}_UB'][experiment]
                parallel_reaction['MEDIAN'] = intracellular_data[f'{name}_MEDIAN'][experiment]
                parallel_reaction['MILP Coefficent'] = milp_coeff

                # For each member of a parallel set, allow both directions when building bounds
                lowerbounds.append(min(c13_lb, c13_ub, -c13_lb, -c13_ub))
                upperbounds.append(max(c13_lb, c13_ub, -c13_lb, -c13_ub))
                medians.append(c13_median)
                milp_objective_coefficients.append(0)

            # If we hit the end of a parallel group, finalize the group entry
            if (('p' in str(pathway)) and ((parallel_indexes[-1] == 'placeholder') or (parallel_indexes[-1] == j))):
                if parallel_reaction and (parallel_reaction not in parallel_reactions):
                    parallel_reactions.append(parallel_reaction)

    # Build dicts keyed by model reaction name for LB/UB/MEDIAN
    lowerbounds_dict = dict(zip(next_flux_constrained_reactions, lowerbounds))
    upperbounds_dict = dict(zip(next_flux_constrained_reactions, upperbounds))
    medians_dict = dict(zip(next_flux_constrained_reactions, medians))

    # Attach columns back onto the mapping table for inspection/export
    mapping_table_out = mapping_table.copy()
    mapping_table_out['LB'] = lowerbounds
    mapping_table_out['UB'] = upperbounds
    mapping_table_out['MEDIAN'] = medians
    mapping_table_out['Objective Coefficient'] = milp_objective_coefficients

    # Model-level lower/upper bounds list (in GEM reaction order)
    model_lowerbounds = [rxn.lower_bound for rxn in GEM.reactions]
    model_upperbounds = [rxn.upper_bound for rxn in GEM.reactions]

    # Build series vs reversible lists for the Pyomo objective terms
    desired_reactions = list(mapping_table['Model Name'])
    series_irreversible_names = []
    series_reversible_forward_names = []
    series_reversible_backward_names = []
    series_irreversible_indices = []
    series_reversible_forward_indices = []
    series_reversible_backward_indices = []
    series_irreversible_c13 = []
    series_reversible_c13 = []
    series_irreversible_obj_coeff = []
    series_reversible_obj_coeff = []

    for i in range(len(desired_reactions)):
        rname = desired_reactions[i]
        c13 = medians[i]
        coeff = milp_objective_coefficients[i]

        if coeff == 'nan' or str(c13) == 'nan':
            continue

        if rname + '_reverse' in reaction_ids:
            series_reversible_forward_names.append(rname)
            series_reversible_forward_indices.append(reaction_ids.index(rname))
            series_reversible_backward_names.append(rname + '_reverse')
            series_reversible_backward_indices.append(reaction_ids.index(rname + '_reverse'))
            series_reversible_c13.append(c13)
            series_reversible_obj_coeff.append(coeff)
        else:
            series_irreversible_names.append(rname)
            series_irreversible_indices.append(reaction_ids.index(rname))
            series_irreversible_c13.append(c13)
            series_irreversible_obj_coeff.append(coeff)

    # For each parallel group, add reverse partners and indices
    for k in range(len(parallel_reactions)):
        model_names = list(parallel_reactions[k]["Model Names"])
        model_directionality = list(parallel_reactions[k]["Model Directionality"])
        # Add reverse if present and track opposite direction
        # iterate over a snapshot of current names to avoid infinite loop
        orig_len = len(model_names)
        for j in range(orig_len):
            name = model_names[j]
            if name + '_reverse' in reaction_ids:
                model_names.append(name + '_reverse')
                model_directionality.append(-model_directionality[j])
        parallel_reactions[k]["Model Names"] = model_names
        parallel_reactions[k]["Model Directionality"] = model_directionality
        parallel_reactions[k]["Model Indices"] = [reaction_ids.index(nm) for nm in model_names]

    return {
        "reaction_ids": reaction_ids,
        "next_flux_constrained_reactions": next_flux_constrained_reactions,
        "mapping_table_out": mapping_table_out,
        "lowerbounds_dict": lowerbounds_dict,
        "upperbounds_dict": upperbounds_dict,
        "medians_dict": medians_dict,
        "parallel_reactions": parallel_reactions,
        "series_irreversible_names": series_irreversible_names,
        "series_irreversible_indices": series_irreversible_indices,
        "series_irreversible_c13": series_irreversible_c13,
        "series_irreversible_obj_coeff": series_irreversible_obj_coeff,
        "series_reversible_forward_names": series_reversible_forward_names,
        "series_reversible_forward_indices": series_reversible_forward_indices,
        "series_reversible_backward_names": series_reversible_backward_names,
        "series_reversible_backward_indices": series_reversible_backward_indices,
        "series_reversible_c13": series_reversible_c13,
        "series_reversible_obj_coeff": series_reversible_obj_coeff,
        "lowerbounds": model_lowerbounds,
        "upperbounds": model_upperbounds,
    }