# TODO: 1) a parser for iqtree 2) a func to load all params in cogent3
#       3) there should be functions to reconstruct relevant parameters.
from cogent3.core.tree import PhyloNode

from phylo_limits.check_boundary import ParamRules
from phylo_limits.classify_matrix import ModelPsubs


def load_bounds(model_result) -> ...:
    """get all boundary settings in a model fit"""
    rules = model_result.lf.get_param_rules()
    bounds = dict()
    for rule in rules:
        if rule["par_name"] == "length":
            bounds["branch_upper"] = rule["upper"]
            bounds["branch_lower"] = rule["lower"]
            break

    for rule in rules:
        if (rule["par_name"] != "length") & (rule["par_name"] != "mprobs"):
            bounds["rate_upper"] = rule["upper"]
            bounds["rate_lower"] = rule["lower"]
            break
    return ...


def load_topo(model_result) -> PhyloNode:
    """get c3 tree"""
    return model_result.lf.tree


def load_psubs(model_result) -> ModelPsubs:
    """get psubs"""
    return ModelPsubs(source=model_result.source, psubs=model_result.lf.get_all_psubs())


def load_param_values(model_result) -> ParamRules:
    """get non-topology param values"""
    return ParamRules(
        source=model_result.source, params=model_result.lf.get_param_rules()
    )


def parsing_iqtree():
    """_summary_

    Returns:
        dict: _description_
    """
    pass
