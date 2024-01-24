import pytest
import numpy
from ylib.check_dlc import check_dlc, check_chainsaw


@pytest.mark.repeat(100)
def test_check_dlc(make_p):
    dlc = numpy.random.rand() > 0.5
    p = make_p(dlc)
    assert check_dlc(p) == dlc


@pytest.mark.repeat(100)
def test_check_chainsaw(make_p):
    m = make_p()
    row_indices = numpy.arange(4)
    numpy.random.shuffle(row_indices)
    dlc = not (row_indices == numpy.arange(4)).all()
    p = m[row_indices, :]
    assert check_chainsaw(p) == dlc


@pytest.mark.repeat(100)
def test_check_chainsaw_limit():
    a = numpy.random.random(4)
    a /= a.sum()
    m = numpy.array([a - 1e-8] * 4)
    assert check_chainsaw(p_matrix=m, p_limit=numpy.array([a] * 4)) == False


def test_check_chainsaw_close():
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
    assert check_chainsaw(m) == False
