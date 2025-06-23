import dataclasses

from enum import Enum

import numpy

from cogent3.util.dict_array import DictArray
from cogent3.core.table import Table
from numpy import allclose, eye, ndarray, tile

from phylim._version import __version__


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
    p_limit = tile(p_matrix[0], (p_matrix.shape[0], 1))
    return allclose(p_matrix, p_limit)


def is_dlc(p_matrix: ndarray) -> bool:
    """
    Judge whether the given matrix is DLC. IMPORTANT: whether it is in limit distribution does not matter,
    but the equality between the diagnoal and off-diag elements matters.
    """
    diags = numpy.diag(p_matrix)
    off_diags = p_matrix.T[~eye(p_matrix.shape[0], dtype=bool)].reshape(
        p_matrix.shape[0], p_matrix.shape[0] - 1
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
    psubs: dict[tuple[str, ...], DictArray]

    def items(self):
        return self.psubs.items()


@dataclasses.dataclass(slots=True)
class ModelMatrixCategories:
    source: str
    mcats: dict[tuple[str, ...], MatrixCategory]

    def items(self):
        return self.mcats.items()

    def to_rich_dict(self) -> dict:
        return {
            "source": self.source,
            "mcats": {k: v.value for k, v in self.items()},
            "version": __version__,
        }

    def to_table(self) -> Table:
        headers = [
            "edge name",
            "matrix category",
        ]
        rows = []
        rows.extend([edge[0], mcat.value] for edge, mcat in self.items())
        return Table(
            header=headers, data=rows, title="Substitution Matrices Categories"
        )

    def _repr_html_(self) -> str:
        table = self.to_table()
        table.set_repr_policy(show_shape=False)
        return table._repr_html_()


def classify_matrix(psubs: ModelPsubs) -> ModelMatrixCategories:
    """labels all psubs in a given ModelPsubs object which has source info"""
    labelled_psubs_dict = {}
    for key, value in psubs.items():
        p_matrix = value.to_array()
        labelled_psubs_dict[key] = classify_psub(p_matrix)

    return ModelMatrixCategories(source=psubs.source, mcats=labelled_psubs_dict)
