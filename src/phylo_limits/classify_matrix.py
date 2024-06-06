import dataclasses

from enum import Enum

import numpy

from cogent3.app.composable import define_app
from cogent3.util.dict_array import DictArray
from numpy import allclose, eye, ndarray, tile


class MatrixCategory(Enum):
    identity = "identity"
    sympathetic = "sympathetic"
    limit = "limit"
    dlc = "DLC"
    chainsaw = "chainsaw"


IDENTITY = MatrixCategory.identity
SYMPATHETIC = MatrixCategory.sympathetic
LIMIT = MatrixCategory.limit
DLC = MatrixCategory.dlc
CHAINSAW = MatrixCategory.chainsaw


def is_identity(p_matrix: ndarray) -> bool:
    return allclose(p_matrix, eye(p_matrix.shape[0]))


def is_limit(p_matrix: ndarray) -> bool:
    """check if a given matrix is a Limit matrix, which all rows are same"""
    p_limit = tile(p_matrix[0], (4, 1))  # TODO; shape!
    return allclose(p_matrix, p_limit)


def is_dlc(p_matrix: ndarray) -> bool:
    """
    Judge whether the given matrix is DLC. IMPORTANT: whether it is in limit distribution does not matter,
    but the equality between the diagnoal and off-diag elements matters.
    """
    diags = numpy.diag(p_matrix)
    off_diags = p_matrix.T[~eye(4, dtype=bool)].reshape(  # TODO; shape!
        4, 3
    )  # take all off-diags as each column in p_matrix a vector
    for i in range(len(diags)):
        off_max = off_diags[i].max()
        diag = diags[i]
        if off_max >= diag or allclose(off_max, diag):
            return False
    return True


def is_chainsaw(p_matrix: ndarray) -> bool:
    """
    Judge whether the given matrix is a chainsaw. IMPORTANT: whether it is in limit distribution does
    not matter, but the equality between the diagnoal and off-diag elements matters.
    """
    max_indices = p_matrix.argmax(axis=0)
    unique = set(max_indices)
    if len(unique) != p_matrix.shape[0]:
        return False
    if (
        max_indices == numpy.arange(p_matrix.shape[0])
    ).all():  # make sure it's not DLC preliminarily
        return False
    return is_dlc(p_matrix[max_indices, :])


def classify_psub(p_matrix: ndarray) -> MatrixCategory:
    """Take a p_matrix and label it"""
    if is_identity(p_matrix):
        return IDENTITY
    elif is_limit(p_matrix):
        return LIMIT
    elif is_dlc(p_matrix):
        return DLC
    elif is_chainsaw(p_matrix):
        return CHAINSAW
    else:
        return SYMPATHETIC


@dataclasses.dataclass(slots=True)
class ModelPsubs:
    source: str
    psubs: dict[tuple[str], DictArray]

    def items(self):
        return self.psubs.items()


@dataclasses.dataclass(slots=True)
class ModelMatrixCategories:
    source: str
    mcats: dict[tuple[str], MatrixCategory]


@define_app
class classify_psubs:
    """labels all psubs in a given ModelPsubs object which has source info"""

    def main(self, psubs: ModelPsubs) -> ModelMatrixCategories:
        labelled_psubs_dict = {}
        for key, value in psubs.items():
            p_matrix = value.to_array()
            labelled_psubs_dict[key] = classify_psub(p_matrix)

        return ModelMatrixCategories(source=psubs.source, mcats=labelled_psubs_dict)