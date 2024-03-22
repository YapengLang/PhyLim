import warnings

from collections import Counter

import numpy

from cogent3 import maths
from cogent3.app.result import model_result
from cogent3.evolve.parameter_controller import AlignmentLikelihoodFunction
from cogent3.util.dict_array import DictArray
from numpy import allclose, eye, partition
from numpy.linalg import eig, inv  # , sum
from numpy.ma import dot as ma_dot
from numpy.ma.core import array, diag


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

#TODO: rename as is_..
def is_identity(p_matrix) -> bool:
    return allclose(p_matrix, eye(p_matrix.shape[0]))


def is_limit(p_matrix, p_limit) -> bool:
    """check if a given matrix is a limit matrix by definition"""
    return allclose(p_matrix, p_limit)


# modified, 1) if limit_m provided, the dlc should NOT equal to it, 2)the diag NOT allclose() to the largest offdiag #todo: add test
def is_dlc(p_matrix, p_limit=None) -> bool:
    if p_limit is not None and is_limit(p_matrix, p_limit):
        return False
    for i in range(len(p_matrix)):  # get column
        column = [p_matrix[j][i] for j in range(len(p_matrix))]
        diagonal = column[i]  # get diagonal value
        freq = Counter(column)
        if allclose(partition(column, -2)[-1], partition(column, -2)[-2]):
            return False
        if diagonal != max(column) or freq[max(column)] != 1:
            return False
    return True


# modified, 1) if limit_m provided, the chainsaw should NOT equal to it, 2) the diag NOT allclose() to the largest offdiag.
def is_chainsaw(p_matrix, p_limit=None) -> bool:
    """given a p is non-DLC"""  # could simplify a little bit more
    if is_dlc(p_matrix):
        return False
    if p_limit is not None and is_limit(p_matrix, p_limit):
        return False
    largest_row_index = {}
    for i in range(len(p_matrix)):  # get column
        column = [p_matrix[j][i] for j in range(len(p_matrix))]
        freq = Counter(column)
        if freq[max(column)] == 1 and not allclose(
            partition(column, -2)[-1], partition(column, -2)[-2]
        ):
            largest_row_index[i] = column.index(max(column))
        else:
            return False
    return len(set(largest_row_index.values())) == 4




def classify_psubs(
    lf:AlignmentLikelihoodFunction,
    strictly=True,
    label_L=True,
) -> dict:
    """get a dict for all psubs in a fitted likelihood function. add class of each matrix to the dict. 
    always warning users if Limit occurs.                                                                                           
    
    Args: `strictly` controls whether take I as DLC, if False, make it as DLC and warn it
               `label_L` if True, label a limit matrix as Limit instead of Sympathetic

    Returns:
        dict: with class of matrices
    """
    motif_probs=lf.get_motif_probs_by_node()
    psubs_dict = lf.get_all_psubs()
    model_name = lf.name
    
    new_dict = {}
    for key, value in psubs_dict.items():
        pi = (
            get_stat_pi_via_eigen(value)
            if model_name in {"ssGN", "GN"}
            else motif_probs[key[0]]
        )

        if is_identity(value):
            if strictly == True:
                new_dict[key] = {"value": value, "class": "Identity"}
            else:
                warnings.warn(
                    "Fit contains Identity matrix, which makes topology un-identifiable!"
                )
                new_dict[key] = {"value": value, "class": "DLC"}
            continue

        elif is_dlc(p_matrix=value, p_limit=numpy.array([pi, pi, pi, pi])):
            new_dict[key] = {"value": value, "class": "DLC"}
            continue

        elif is_chainsaw(
            p_matrix=value, p_limit=numpy.array([pi, pi, pi, pi])
        ):

            new_dict[key] = {"value": value, "class": "Chainsaw"}
            continue

        else:
            if is_limit(
                p_matrix=value, p_limit=numpy.array([pi, pi, pi, pi])
            ):
                warnings.warn("Fit contains Limit matrix!")
                if label_L == True:
                    new_dict[key] = {"value": value, "class": "Limit"}
                else:
                    new_dict[key] = {"value": value, "class": "Sympathetic"}
            else:
                new_dict[key] = {"value": value, "class": "Sympathetic"}
    return new_dict

