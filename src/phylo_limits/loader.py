# TODO: 1) a parser for iqtree 2) a func to load all params in cogent3
#       3) there should be functions to reconstruct relevant parameters.
import json

from cogent3.app.composable import define_app
from cogent3.app.result import model_result
from cogent3.app.typing import SerialisableType
from cogent3.core.tree import PhyloNode
from cogent3.evolve.parameter_controller import AlignmentLikelihoodFunction

from phylo_limits.check_boundary import diagonse
from phylo_limits.check_ident import has_valid_path
from phylo_limits.loader import load_bounds
from phylo_limits.matrix_class import classify_psubs
from phylo_limits.record import PhyloLimitRec


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


def parsing_iqtree():
    """_summary_

    Returns:
        dict: _description_
    """
    pass


@define_app
def generate_record(model_res: model_result, strictly=False) -> str:
    """record psubs classes, identifiability, boundary values and non-DLC projection.
    Args:
        "strictly" controls the sensitivity for Identity matrix (I); if false, treat I as DLC.
    Return:
        string in json format
    """
    tree = model_res.lf.tree.get_newick(with_node_names=True)
    psubs_class = classify_psubs(lf=model_res.lf, strictly=strictly)
    boundary_values = diagonse(
        statistics=model_res.lf.get_statistics(), bounds=load_bounds(model_res.lf)
    )
    bad_nodes = has_valid_path(lf=model_res.lf, strictly=strictly)
    identifiable = not bad_nodes  # if there are any bad nodes show up, unidentifiable

    rec = PhyloLimitRec(
        source=model_res.source,
        model_name=model_res.name,
        tree=tree,
        psubs_class=psubs_class,
        boundary_values=boundary_values,
        bad_nodes=bad_nodes,
        identifiable=identifiable,
    )

    return json.dumps(rec.to_rich_dict())
