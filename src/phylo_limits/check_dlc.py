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
