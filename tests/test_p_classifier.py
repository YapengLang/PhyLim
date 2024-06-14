import numpy
import pytest

from phylo_limits.classify_matrix import (
    CHAINSAW,
    DLC,
    IDENTITY,
    LIMIT,
    SYMPATHETIC,
    ModelMatrixCategories,
    ModelPsubs,
    classify_psub,
    classify_psubs,
    is_chainsaw,
    is_dlc,
    is_identity,
    is_limit,
)


@pytest.mark.parametrize(
    "matrix_input, expected",
    [
        (numpy.eye(4), True),
        (numpy.ones((4, 4)), False),
    ],
)
def test_is_identity(matrix_input, expected):
    assert is_identity(matrix_input) == expected


@pytest.mark.repeat(10)
def test_is_dlc(make_dlc, repeat):
    for _ in range(repeat):
        assert is_dlc(make_dlc()) == True


@pytest.mark.repeat(10)
def test_is_dlc_close_elem(make_dlc, repeat):
    for _ in range(repeat):
        m = make_dlc()
        m[1, 0] = m[0, 0]
        assert is_dlc(m) == False


def test_is_dlc_all_rows_same(make_limit):
    assert is_dlc(make_limit()) == False


def test_is_dlc_by_chainsaw(make_chainsaw):
    assert is_dlc(make_chainsaw()) == False


@pytest.mark.repeat(10)
def test_is_chainsaw(make_chainsaw, repeat):
    for _ in range(repeat):
        assert is_chainsaw(make_chainsaw()) == True


@pytest.mark.repeat(10)
def test_is_chainsaw_two_rows_same(make_dlc, repeat):
    for _ in range(repeat):
        m = make_dlc()
        m[0] = m[1]
        assert is_chainsaw(m) == False


def test_is_chainsaw_all_rows_same(make_limit):
    assert is_chainsaw(make_limit()) == False


def test_is_chainsaw_by_dlc(make_dlc):
    assert is_chainsaw(make_dlc()) == False


def test_is_chainsaw_by_sym_empir():
    with open("data/matrices/sympathetic1.npy", "rb") as f:
        m = numpy.load(f)
    assert is_chainsaw(m) == False


def test_is_limit(make_limit):
    assert is_limit(make_limit()) == True


def test_is_limit_by_sym_empir():
    with open("data/matrices/sympathetic1.npy", "rb") as f:
        m = numpy.load(f)
    assert is_limit(m) == False


@pytest.mark.parametrize(
    "mtx_input, expected",
    [("make_limit", LIMIT), ("make_dlc", DLC), ("make_chainsaw", CHAINSAW)],
)
def test_classify_psub_by_fixtures(mtx_input, expected, request):
    mtx_fx = request.getfixturevalue(mtx_input)
    result = classify_psub(mtx_fx())
    assert result == expected


def test_classify_psub_by_sym_empir():
    with open("data/matrices/sympathetic1.npy", "rb") as f:
        m = numpy.load(f)
    assert classify_psub(m) == SYMPATHETIC


def test_classify_psub_by_I():
    assert classify_psub(numpy.eye(4)) == IDENTITY


def test_classify_psubs(make_dlc):
    from cogent3.util import dict_array

    psub = {(str("bar"),): dict_array.DictArray(make_dlc())}
    mpsubs = ModelPsubs(source="foo", psubs=psub)
    cls_app = classify_psubs()
    assert isinstance(cls_app(mpsubs), ModelMatrixCategories)
