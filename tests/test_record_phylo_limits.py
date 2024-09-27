import pathlib

from cogent3.util.deserialise import deserialise_object

from phylo_limits.check_boundary import BoundsViolation, ParamRules
from phylo_limits.classify_matrix import ModelMatrixCategories, ModelPsubs
from phylo_limits.eval_identifiability import IdentCheckRes
from phylo_limits.record_phylo_limits import (
    PhyloLimitRec,
    check_fit_boundary,
    check_identifiability,
    classify_model_psubs,
    generate_phylo_limit_record,
    load_param_values,
    load_psubs,
)


DATADIR = pathlib.Path(__file__).parent / "data"


def test_load_param_values():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    result = load_param_values(model_res)
    assert isinstance(result, ParamRules) == True


def test_load_psubs():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    result = load_psubs(model_res)
    assert isinstance(result, ModelPsubs) == True


def test_generate_record():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )  # two I
    rec_app = generate_phylo_limit_record()  # default `strict` == F
    record = rec_app(model_res)
    assert isinstance(record, PhyloLimitRec) == True
    assert record.strict == False


def test_generate_record_strict_control():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )  # two I
    rec_app = generate_phylo_limit_record(strict=True)  # set `strict`
    record = rec_app(model_res)
    assert isinstance(record, PhyloLimitRec) == True
    assert record.strict == True
    assert record.identifiability == False


def test_to_rich_dict_phylolimitrec():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    rec_app = generate_phylo_limit_record(strict=True)
    record = rec_app(model_res)
    result = record.to_rich_dict()
    assert isinstance(result, dict) == True
    assert all(
        k in result
        for k in [
            "model_name",
            "boundary_values",
            "ISCL_mcats",
            "source",
            "identifiability",
            "strict",
            "message",
            "version",
        ]
    )


def test_check_identifiability_detailed():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    check_app = check_identifiability(strict=True, details=True)
    record = check_app(model_res)
    assert isinstance(record, IdentCheckRes) == True


def test_check_identifiability():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    check_app = check_identifiability(strict=True, details=False)
    record = check_app(model_res)
    assert record == False


def test_classify_model_psubs():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    check_app = classify_model_psubs()
    res = check_app(model_res)
    assert isinstance(res, ModelMatrixCategories) == True


def test_check_fit_boundary():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    check_app = check_fit_boundary()
    res = check_app(model_res)
    assert isinstance(res, BoundsViolation) == True
