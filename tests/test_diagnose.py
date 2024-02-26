import numpy
import pytest

from phylo_limits.diagnose import mles_within_bounds


_eps = numpy.finfo(float).eps

rate_params = (
    numpy.full((3, 11), 10),  # fine params
    numpy.full((3, 11), 1e-6),  # exact lower bound
    numpy.full((3, 11), 1e-6 + _eps),  # close to lower bound
    numpy.full((3, 11), 200),  # exact upper bound
    numpy.full((3, 11), 200 - _eps),  # close to upper bound
)

#TODO: decompose this test onto rate params only
@pytest.mark.parametrize("params", rate_params)
def test_mles_within_bounds(params):
    pass
