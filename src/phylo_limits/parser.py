from cogent3.core.tree import PhyloNode

from phylo_limits.check_boundary import ParamRules
from phylo_limits.classify_matrix import ModelPsubs


def load_tree(model_result) -> PhyloNode:
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
