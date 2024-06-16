from cogent3.core.tree import PhyloNode
from cogent3.util.deserialise import deserialise_object

from phylo_limits.check_boundary import ParamRules
from phylo_limits.classify_matrix import ModelPsubs
from phylo_limits.parser import load_param_values, load_psubs, load_tree


def test_load_param_values():
    model_res = deserialise_object("data/eval_identifiability/unid_model_result.json")
    result = load_param_values(model_res)
    assert isinstance(result, ParamRules) == True


def test_load_psubs():
    model_res = deserialise_object("data/eval_identifiability/unid_model_result.json")
    result = load_psubs(model_res)
    assert isinstance(result, ModelPsubs) == True


def test_load_tree():
    model_res = deserialise_object("data/eval_identifiability/unid_model_result.json")
    result = load_tree(model_res)
    assert isinstance(result, PhyloNode) == True
