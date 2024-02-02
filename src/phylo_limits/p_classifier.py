import numpy
from cogent3.util.dict_array import DictArray
from cogent3.app.result import model_result
from phylo_limits import delta_col

import warnings

from collections import Counter
from numpy import allclose, partition, eye


def check_I(p_matrix) -> bool:
    return allclose(p_matrix, eye(p_matrix.shape[0]))


def check_limit(p_matrix, p_limit) -> bool:
    """check if a given matrix is a limit matrix by definition"""
    return bool(allclose(p_matrix, p_limit))


# modified, 1) if limit_m provided, the dlc should NOT equal to it, 2)the diag NOT allclose() to the largest offdiag #todo: add test
def check_dlc(p_matrix, p_limit=None) -> bool:
    if p_limit is not None and check_limit(p_matrix, p_limit):
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
def check_chainsaw(p_matrix, p_limit=None) -> bool:
    """given a p is non-DLC"""  # could simplify a little bit more
    if check_dlc(p_matrix):
        return False
    if p_limit is not None and check_limit(p_matrix, p_limit):
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




def check_all_psubs(
    model_res:model_result,
    strictly=True,
    label_L=False,
) -> dict:
    """get a dict for all psubs in a fitted likelihood function. add class of each matrix to the dict. always warning users if Limit occurs
    Parameter: `strictly` controls whether take I as DLC, if False, make it as DLC and warn it
               `label_L` if True, label a limit matrix as Limit instead of Sympathetic"""
    # motif_probs=lf.get_motif_probs_by_node()
    psubs_dict = model_res.lf.get_all_psubs()
    model_name = model_res.name
    motif_probs = model_res.lf.get_motif_probs_by_node()

    new_dict = {}
    for key, value in psubs_dict.items():
        pi = (
            delta_col.get_stat_pi_via_eigen(value)
            if model_name in {"ssGN", "GN"}
            else motif_probs[key[0]]
        )

        if check_I(value):
            if strictly == True:
                new_dict[key] = {"value": value, "class": "Identity"}
            else:
                warnings.warn(
                    "Fit contains Identity matrix, which makes topology un-identifiable!"
                )
                new_dict[key] = {"value": value, "class": "DLC"}
            continue

        elif check_dlc(p_matrix=value, p_limit=numpy.array([pi, pi, pi, pi])):
            new_dict[key] = {"value": value, "class": "DLC"}
            continue

        elif check_chainsaw(
            p_matrix=value, p_limit=numpy.array([pi, pi, pi, pi])
        ):

            new_dict[key] = {"value": value, "class": "Chainsaw"}
            continue

        else:
            if check_limit(
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
