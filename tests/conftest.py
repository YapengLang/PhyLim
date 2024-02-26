import numpy
import pytest


@pytest.fixture
def repeat(request):
    # Get the repeat count from the test marker or default to 1
    repeat_count = request.node.get_closest_marker('repeat')
    if repeat_count is None:
        repeat_count = 1
    else:
        repeat_count = repeat_count.args[0]
    return repeat_count

@pytest.fixture()
def make_p():
    def _make_p(dlc=True):
        indices = numpy.diag_indices(4)
        rng = numpy.random.default_rng()
        v = rng.random(size=(2,))
        m = numpy.zeros((4, 4), dtype=float)
        diag, off_diag = v.max(), v.min() / 4
        if not dlc:
            diag, off_diag = off_diag, diag
        m[:] = off_diag
        m[indices] = diag
        m /= m.sum(axis=1)
        return m

    return _make_p


@pytest.fixture()
def make_q():
    rng = numpy.random.default_rng()
    q = rng.random(size=(4, 4))
    indices = numpy.diag_indices(4)
    q[indices] = 0
    q[indices] = -q.sum(axis=1)
    return q
