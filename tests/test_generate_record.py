from cogent3.util.deserialise import deserialise_object

from phylo_limits.record_limit import PhyloLimitRec, generate_record


def test_generate_record():
    model_res = deserialise_object(
        "data/eval_identifiability/unid_model_result.json"
    )  # two I
    rec_app = generate_record()  # default `strict` == F
    record = rec_app(model_res)
    assert isinstance(record, PhyloLimitRec) == True
    assert record.strict == False


def test_generate_record_strict_control():
    model_res = deserialise_object(
        "data/eval_identifiability/unid_model_result.json"
    )  # two I
    rec_app = generate_record(strict=True)  # set `strict`
    record = rec_app(model_res)
    assert isinstance(record, PhyloLimitRec) == True
    assert record.strict == True
    assert record.identifiability == False
