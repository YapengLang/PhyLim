from cogent3 import make_table, maths
from cogent3.util.table import Table
from cogent3.app.result import model_result
from cogent3.app.typing import SerialisableType
from cogent3.app.composable import define_app
import numpy

from phylo_limits.check_dlc import check_dlc, check_chainsaw
from phylo_limits.delta_col import (
    min_col_diff,
    CalcMinColDeltaFit,
    get_stat_pi_via_eigen,
    get_node_parent,
)

# input a list NotCompleted, extract bound violatio infor, and write as a tsv.
def get_bound_viol(notcomps: list) -> Table:
    rows = []
    for notcomp in notcomps:
        name = notcomp.source
        msg = notcomp.message
        msg = msg[1::2]
        row = [name] + list(msg)
        rows.append(row)

    header = ["name"] + ["within branch lower"] + ["within branch upper"]
    if len(rows[0]) == 5:
        header = header + ["within rate lower"] + ["within rate upper"]

    return make_table(header, data=rows)


# input model_result, validate matrix, and output list for only non-dlc.
@define_app
def get_dlc_viol(result: model_result) -> list:
    # unique_id = result.source
    edges = [k[0] for k in result.lf.get_all_psubs().keys()]
    rows = []
    for name in edges:
        P = result.lf.get_psub_for_edge(name=name)
        Q = result.lf.get_rate_matrix_for_edge(name=name)  # here the q is calibrated.
        t = result.lf.get_param_value(par_name="length", edge=name)
        p0 = result.lf.get_motif_probs_by_node()[
            get_node_parent(result=result, node=name)
        ]
        pi = (
            get_stat_pi_via_eigen(P) if result.name in ["ssGN", "GN"] else p0
        )  # get limit distribution
        ens = maths.matrix_exponential_integration.expected_number_subs(p0=p0, Q=Q, t=t)

        metric = CalcMinColDeltaFit(
            Q
        )  # get diag, off-diag indices. this line can be replaced.
        delta_col = min_col_diff(
            P.to_array(), metric.diag_indices, metric.offdiag_indices
        )
        if check_dlc(p_matrix=P, p_limit=numpy.array([pi, pi, pi, pi])):
            row = [name] + [t] + [ens] + [delta_col] + [0] + [0]
        elif check_chainsaw(p_matrix=P, p_limit=numpy.array([pi, pi, pi, pi])):
            row = [name] + [t] + [ens] + [delta_col] + [0] + [1]
        else:
            row = [name] + [t] + [ens] + [delta_col] + [1] + [0]
        rows.append(row)
    return rows


# input Cogent3.table, output algn name, edge, etc..
@define_app
def get_transit(table: Table) -> list:
    rows = []
    for i in range(table.shape[0]):
        if table[i, 4] == 1:
            row = table[i, :-1].tolist()[0]
            rows.append(row)
    return rows
