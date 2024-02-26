from typing import Any

import numpy
import scipy

from cogent3 import maths
from numpy.linalg import eig, inv  # , sum
from numpy.ma import dot as ma_dot
from numpy.ma.core import array, diag


expm = scipy.linalg.expm

_eps = numpy.finfo(float).eps

# create default features; a 4-state matrix
mask = numpy.ones((4, 4), bool)
mask[numpy.diag_indices(4)] = False
DIAG = numpy.diag_indices(4)
OFFDIAG = mask


def min_diff_from_diag(
    m: numpy.ndarray, diag_indices: numpy.ndarray, off_diag_indices: numpy.ndarray
) -> numpy.ndarray:
    """
    Args:
        m (numpy.ndarray): a matrix
        diag_indices (numpy.ndarray): diagonal indices
        off_diag_indices (numpy.ndarray): off-diagonal indices in boolean

    Des:
        compute difference for each column between diagonal and
        the largest off-diagonal element
    """
    return m[diag_indices] - m.max(axis=0, where=off_diag_indices, initial=m.min())


def min_col_diff(
    m: numpy.ndarray,
    diag_indices: numpy.ndarray = DIAG,
    off_diag_indices: numpy.ndarray = OFFDIAG,
) -> numpy.ndarray:
    return min_diff_from_diag(m, diag_indices, off_diag_indices).min()


class CalcMinColDeltaFit:
    """Represents the property of matrix. attributes: diag_indices, offdiag_indices"""

    def __init__(self, Q):
        self.Q = Q
        num_states = Q.shape[0]
        mask = numpy.ones(Q.shape, bool)
        mask[numpy.diag_indices(num_states)] = False
        self.diag_indices = numpy.diag_indices(num_states)
        self.offdiag_indices = mask

    def __call__(self, new_tau):
        # compute the DLC metric
        l_0 = -2 * _eps

        P = expm(self.Q * new_tau)
        col_delta = min_col_diff(P, self.diag_indices, self.offdiag_indices)
        return (
            l_0 - col_delta
        ) ** 2  # minimise the point to get where slightly passes the 0


class EstimateMinColDelta:
    def __init__(self, Q):
        self.calc = CalcMinColDeltaFit(Q)
        self.stat = 0.0  # tau

    def optimise(self, max_evaluations=1000):
        # force exit if the evaluation step over the limit
        try:
            tau = maths.optimisers.minimise(
                self.calc,
                xinit=(1,),  # the initial value
                bounds=((0,), (1000,)),  # [lower,upper] bounds for the parameter
                local=True,  # just local optimisation, not Simulated Annealing
                show_progress=False,
                limit_action="raise",
                max_evaluations=max_evaluations,
                return_calculator=True,
                tolerance=1e-8,
            )

        except ArithmeticError as err:
            print(err)

        self.stat = tau

        self.fit = min_col_diff(
            expm(self.calc.Q * tau), self.calc.diag_indices, self.calc.offdiag_indices
        )  # modified the self.fit as fitted delta_col


# 24/7 added check of limit distribution achieved todo: modify docs
def valid_symp(limit, to_check):
    return bool(numpy.allclose(limit, to_check))


def get_stat_pi_via_eigen(P, check_precision=True):
    """This code was provided by Gavin Huttley Obtain stationary distribution
    via EigenValue decomposition."""
    P = array(P).T
    eva, eve = eig(P)

    if check_precision:
        i = inv(eve)
        r = ma_dot(ma_dot(eve, diag(eva)), i).real
        if not numpy.allclose(P, r):
            raise ArithmeticError

    evect = eve[:, eva.round(10) == 1]
    stat_pi = evect / sum(evect.real)
    return stat_pi.flatten().real


def valid_bysect(delta):
    return delta <= -(2 * _eps)


def invalid_dlc(delta: float, limit, to_check) -> bool:
    # if this condition is ture, the dlc is violated.
    # the doubled machine precison makes sure is truly a dlc violation
    return bool(valid_bysect(delta) and not valid_symp(limit, to_check))
    # check the deviation #todo: check the avail of 0.5*_eps

