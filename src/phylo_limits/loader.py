# TODO: 1) a parser for iqtree 2) a func to load all params in cogent3
#       3) there should be functions to reconstruct relevant parameters.


def load_bounds(lf) -> dict:
    """get all boundary settings in a model fit"""
    rules = lf.get_param_rules()
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
    return bounds


def load_topo(model_result) -> str:
    """get newick format tree"""
    return model_result.lf.tree.get_newick(with_node_names=True)


def load_param_values(model_result) -> list[dict]:
    """get non-topology param values"""
    return model_result.lf.get_param_rules()


def parsing_iqtree() -> dict:
    """_summary_

    Returns:
        dict: _description_
    """
    pass
