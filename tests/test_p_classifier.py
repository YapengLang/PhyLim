import numpy
import pytest

from phylo_limits.classify_matrix import (
    classify_psub,
    classify_psubs,
    is_chainsaw,
    is_dlc,
    is_identity,
    is_limit,
)


def test_is_identity():
    assert is_identity(numpy.eye(4)) == True


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


# def test_check_chainsaw_limit(repeat):
#     for i in range(repeat):
#         a = numpy.random.random(4)
#         a /= a.sum()
#         m = numpy.array([a - 1e-8] * 4)
#         assert is_chainsaw(p_matrix=m, p_limit=numpy.array([a] * 4)) == False


def test_check_chainsaw_close():  # TODO: store extreme cases into data folder.
    m = numpy.array(
        [
            [
                0.096737740606487349,
                0.412019393154277491,
                0.087722662266498261,
                0.403520203972736691,
            ],
            [
                0.084159672770039218,
                0.414528031407991659,
                0.086784264411416831,
                0.414528031410552167,
            ],
            [
                0.087722662266498261,
                0.403520203972736913,
                0.096737740606487418,
                0.412019393154277491,
            ],
            [
                0.086784264411416789,
                0.414528031410552333,
                0.084159672770039246,
                0.414528031407991437,
            ],
        ]
    )
    assert is_chainsaw(m) == False


# @pytest.mark.parametrize(
#     "matrix, expected",
#     [
#         (numpy.eye(4, dtype=float), True),
#         (numpy.ones((4, 4)), False),
#     ],
# )
# def test_check_I(matrix, expected):
#     assert is_I(matrix) == expected


# def test_check_limit():
#     a = numpy.random.random(4)
#     a /= a.sum()
#     m = numpy.array([a - 1e-8] * 4)
#     assert is_limit(p_matrix=m, p_limit=numpy.array([a] * 4)) == True
