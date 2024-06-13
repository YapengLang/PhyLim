import numpy
import pytest


@pytest.fixture
def repeat(request):
    # Get the repeat count from the test marker or default to 1
    # refer to: https://docs.pytest.org/en/8.2.x/how-to/fixtures.html
    # #using-markers-to-pass-data-to-fixtures
    repeat_count = request.node.get_closest_marker("repeat")
    if repeat_count is None:
        repeat_count = 1
    else:
        repeat_count = repeat_count.args[0]
    return repeat_count


@pytest.fixture()
def make_dlc():
    def _make_dlc():
        indices = numpy.diag_indices(4)
        rng = numpy.random.default_rng()
        v = rng.random(size=(2,))  # two floats in [0, 1)
        m = numpy.zeros((4, 4), dtype=float)
        diag, off_diag = v.max(), v.min() / 4  # scale down the minimum
        m[:] = off_diag
        m[indices] = diag
        m /= m.sum(axis=1)
        return m

    return _make_dlc


@pytest.fixture()
def make_chainsaw(make_dlc):
    def _make_chainsaw():
        m = make_dlc()
        shuffled = numpy.arange(m.shape[0])
        original = shuffled.copy()
        while numpy.array_equal(shuffled, original):
            numpy.random.shuffle(shuffled)
        return m[shuffled, :]

    return _make_chainsaw


@pytest.fixture()
def make_limit():
    """create limit matrix with rows exactly the same"""

    def _make_limit():
        rng = numpy.random.default_rng()
        v = rng.random(4)
        v /= v.sum(axis=0)
        return numpy.tile(v, (4, 1))

    return _make_limit
