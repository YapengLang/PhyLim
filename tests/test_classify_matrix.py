import pathlib

import numpy
import pytest

from cogent3.util.table import Table

from phylim.classify_matrix import (
    CHAINSAW,
    DLC,
    IDENTITY,
    LIMIT,
    SYMPATHETIC,
    ModelMatrixCategories,
    ModelPsubs,
    classify_matrix,
    classify_psub,
    is_chainsaw,
    is_dlc,
    is_identity,
    is_limit,
)


DATADIR = pathlib.Path(__file__).parent / "data"


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
        assert is_dlc(make_dlc())


@pytest.mark.repeat(10)
def test_is_dlc_close_elem(make_dlc, repeat):
    for _ in range(repeat):
        m = make_dlc()
        m[1, 0] = m[0, 0]
        assert not is_dlc(m)


def test_is_dlc_all_rows_same(make_limit):
    assert not is_dlc(make_limit())


def test_is_dlc_by_chainsaw(make_chainsaw):
    assert not is_dlc(make_chainsaw())


@pytest.mark.repeat(10)
def test_is_chainsaw(make_chainsaw, repeat):
    assert all(is_chainsaw(make_chainsaw()) for _ in range(repeat))


@pytest.mark.repeat(10)
def test_is_chainsaw_two_rows_same(make_dlc, repeat):
    for _ in range(repeat):
        m = make_dlc()
        m[0] = m[1]
        assert not is_chainsaw(m)


def test_is_chainsaw_all_rows_same(make_limit):
    assert not is_chainsaw(make_limit())


def test_is_chainsaw_by_dlc(make_dlc):
    assert not is_chainsaw(make_dlc())


def test_is_chainsaw_by_sym_empir():
    with open(f"{DATADIR}/matrices/sympathetic1.npy", "rb") as f:
        m = numpy.load(f)
    assert not is_chainsaw(m)


def test_is_limit(make_limit):
    assert is_limit(make_limit())


def test_is_limit_by_sym_empir():
    with open(f"{DATADIR}/matrices/sympathetic1.npy", "rb") as f:
        m = numpy.load(f)
    assert not is_limit(m)


@pytest.mark.parametrize(
    "mtx_input, expected",
    [("make_limit", LIMIT), ("make_dlc", DLC), ("make_chainsaw", CHAINSAW)],
)
def test_classify_psub_by_fixtures(mtx_input, expected, request):
    mtx_fx = request.getfixturevalue(mtx_input)
    result = classify_psub(mtx_fx())
    assert result == expected


def test_classify_psub_by_sym_empir():
    with open(f"{DATADIR}/matrices/sympathetic1.npy", "rb") as f:
        m = numpy.load(f)
    assert classify_psub(m) == SYMPATHETIC


def test_classify_psub_by_I():
    assert classify_psub(numpy.eye(4)) == IDENTITY


def test_to_rich_dict_modelmatrixcategories():
    labelled = ModelMatrixCategories(source="foo", mcats={("bar",): CHAINSAW})
    result = labelled.to_rich_dict()
    assert all(k in result for k in ["source", "mcats", "version"])


def test_to_table_modelmatrixcategories():
    labelled = ModelMatrixCategories(source="foo", mcats={("bar",): CHAINSAW})
    result = labelled.to_table()
    assert isinstance(result, Table)


def test_repr_html_modelmatrixcategories():
    labelled = ModelMatrixCategories(source="foo", mcats={("bar",): CHAINSAW})
    result = labelled._repr_html_()
    assert isinstance(result, str)


def test_classify_psubs(make_dlc):
    from cogent3.util import dict_array

    psub = {(str("bar"),): dict_array.DictArray(make_dlc())}
    mpsubs = ModelPsubs(source="foo", psubs=psub)
    assert isinstance(classify_matrix(mpsubs), ModelMatrixCategories)
