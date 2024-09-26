from cogent3.util.deserialise import deserialise_object

from phylo_limits.check_boundary import ParamRules
from phylo_limits.classify_matrix import ModelPsubs
from phylo_limits.record_phylo_limits import (
    PhyloLimitRec,
    generate_phylo_limit_record,
    load_param_values,
    load_psubs,
)


def test_load_param_values():
    model_res = deserialise_object("data/eval_identifiability/unid_model_result.json")
    result = load_param_values(model_res)
    assert isinstance(result, ParamRules) == True


def test_load_psubs():
    model_res = deserialise_object("data/eval_identifiability/unid_model_result.json")
    result = load_psubs(model_res)
    assert isinstance(result, ModelPsubs) == True


def test_generate_record():
    model_res = deserialise_object(
        "data/eval_identifiability/unid_model_result.json"
    )  # two I
    rec_app = generate_phylo_limit_record()  # default `strict` == F
    record = rec_app(model_res)
    assert isinstance(record, PhyloLimitRec) == True
    assert record.strict == False


def test_generate_record_strict_control():
    model_res = deserialise_object(
        "data/eval_identifiability/unid_model_result.json"
    )  # two I
    rec_app = generate_phylo_limit_record(strict=True)  # set `strict`
    record = rec_app(model_res)
    assert isinstance(record, PhyloLimitRec) == True
    assert record.strict == True
    assert record.identifiability == False
