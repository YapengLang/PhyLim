import numpy
import pytest
from ylib.mles_filter import mles_discm, mles_within_bounds
from cogent3.util.deserialise import deserialise_object
from cogent3.app.result import model_result
from cogent3 import open_data_store

# todo: use built-in fixture tmpdir to store temp db
_eps = numpy.finfo(float).eps

rate_params = (
    numpy.full((3, 11), 10),  # fine params
    numpy.full((3, 11), 1e-6),  # exact lower bound
    numpy.full((3, 11), 1e-6 + _eps),  # close to lower bound
    numpy.full((3, 11), 200),  # exact upper bound
    numpy.full((3, 11), 200 - _eps),  # close to upper bound
)

branch_params = (
    numpy.full((3, 1), 10),  # fine params
    numpy.full((3, 1), 1e-6),  # exact lower bound
    numpy.full((3, 1), 1e-6 + _eps),  # close to lower bound
    numpy.full((3, 1), 50),  # exact upper bound
    numpy.full((3, 1), 50 - _eps),  # close to upper bound
)

params_to_try = (
    ((rate_params[0], branch_params[0]), True),  # nice
    ((rate_params[1], branch_params[0]), [1, 1, 0, 1]),  # rate is lower bound
    ((rate_params[2], branch_params[0]), [1, 1, 0, 1]),  # rate close to lower bound
    ((rate_params[0], branch_params[1]), [0, 1, 1, 1]),  # branch length is lower bound
    ((rate_params[1], branch_params[2]), [0, 1, 0, 1]),  # all close to lower bound
    ((rate_params[3], branch_params[4]), [1, 0, 1, 0]),  # all close to upper bound
    ((rate_params[0], branch_params[3]), [1, 0, 1, 1]),  # branch close to upper
    ((rate_params[3], branch_params[0]), [1, 1, 1, 0]),  # rate close to upper bound
    (
        (rate_params[4], branch_params[2]),
        [0, 1, 1, 0],
    ),  # branch close to lower, rate close to upper
)


model_result_to_try = (
    (
        deserialise_object("data/model_result_cases/case1.json"),
        model_result,
    ),
    (  # case 1 , well within
        deserialise_object("data/model_result_cases/case2.json"),
        "FAIL",
    ),
    (  # case 2 , branth length touch the lower bound
        deserialise_object("data/model_result_cases/case3.json"),
        "FAIL",
    ),
    (  # case 3 , rate param too close to the lower
        deserialise_object("data/model_result_cases/case4.json"),
        "FAIL",
    ),
    (  # case 4 , rate param touch the upper
        deserialise_object("data/model_result_cases/case5.json"),
        "FAIL",
    ),  # case 5 , branth length touch the upper
)

ds1 = open_data_store(
    "data/model_result_cases/tstcase_model_result_GTR_boundvio.sqlitedb",
    mode="r",
)

ds2 = open_data_store(
    "data/model_result_cases/tstcase_model_result_GN.sqlitedb",
    mode="r",
)

model_dsmemb_to_try = (
    # case 1 , close the rate lower
    (
        [memb for memb in ds1.completed if memb.unique_id == "159324_878_558178"][0],
        "FAIL",
    ),
    # case 2 , well within
    (
        [memb for memb in ds2.completed if memb.unique_id == "118291-119277"][0],
        model_result,
    ),
    # case 3 , close to branch lower
    (
        [memb for memb in ds2.completed if memb.unique_id == "103447-103710"][0],
        "FAIL",
    ),
)


@pytest.mark.parametrize("params", params_to_try)
def test_mles_discm(params):
    result = mles_discm(params[0][0], params[0][1])
    expected = params[1]

    assert numpy.array_equal(result, expected)


@pytest.mark.parametrize("results", model_dsmemb_to_try)
def test_mles_within_bounds(results):
    """test if the model results over the bound could be picked out"""
    # load out the test suite in json format for the model results above
    discm = mles_within_bounds()

    case = discm(results[0])
    if not isinstance(case, model_result):
        assert case.type == results[1]
        print(case.message)
    else:
        assert isinstance(case, model_result)
