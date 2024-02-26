import numpy
import pytest

from cogent3.maths.optimisers import MaximumEvaluationsReached

from phylo_limits.delta_col import (
    EstimateMinColDelta,
    min_col_diff,
    min_diff_from_diag,
)


_eps = numpy.finfo(float).eps


def test_min_col_diff(make_p,repeat):
    for i in range(repeat):
        mask = numpy.ones((4, 4), bool)
        mask[numpy.diag_indices(4)] = False
        diag_indices = numpy.diag_indices(4)
        offdiag_indices = mask

        # randomly generate p matrix which has 50% probability to be dlc
        dlc = numpy.random.rand() > 0.5
        p = make_p(dlc)
        expected = p[0, 0] - p[0, 1]
        assert numpy.allclose(min_col_diff(p, diag_indices, offdiag_indices), expected)


def test_min_diff_from_diag(make_p,repeat):
    for i in range(repeat):
        mask = numpy.ones((4, 4), bool)
        mask[numpy.diag_indices(4)] = False
        diag_indices = numpy.diag_indices(4)
        offdiag_indices = mask

        dlc = numpy.random.rand() > 0.5
        p = make_p(dlc)
        expected = p[0, 0] - p[0, 1]
        assert numpy.allclose(
            min_diff_from_diag(p, diag_indices, offdiag_indices), [expected] * 4
        )


def test_opt_jc():
    jc = numpy.full((4, 4), 0.25)
    indices = numpy.diag_indices(4)
    jc[indices] = -0.75  # construct jukes-cantor rate m

    opt = EstimateMinColDelta(jc)
    opt.optimise()

    assert abs(0 - opt.fit) < _eps
