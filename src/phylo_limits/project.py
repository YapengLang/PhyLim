#TODO: this file will be filled with functions for projecting the fit: for given likelihood function, report the ens where DLC is broken.
#       (it should project only matrix who is DLC)
from cogent3.parse.table import load_delimited
from cogent3 import get_app, maths
from tqdm import tqdm

import numpy
import scipy

import click

from phylo_limits.check_dlc import check_chainsaw
from phylo_limits.delta_col import get_stat_pi_via_eigen


def num_saw_srch(t, q, p0, tau, model_name):
    start = float(tau) * t
    end = 1000  # change the search upper as 1000, step as 0.1
    step = 0.1

    interval_l = 0.0
    intervals = []
    for i in numpy.arange(start, end, step):
        p = scipy.linalg.expm(q.to_array() * i)
        pi = get_stat_pi_via_eigen(p) if model_name in ["ssGN", "GN"] else p0
        if check_chainsaw(p_matrix=p, p_limit=numpy.array([pi, pi, pi, pi])):
            if interval_l == 0.0:
                interval_l = i
            if i + step >= end:
                interval = [
                    maths.matrix_exponential_integration.expected_number_subs(
                        p0=p0, Q=q, t=interval_l
                    ),
                    maths.matrix_exponential_integration.expected_number_subs(
                        p0=p0, Q=q, t=end
                    ),
                ]

                intervals.append(interval)
        elif interval_l != 0.0:
            interval = [
                maths.matrix_exponential_integration.expected_number_subs(
                    p0=p0, Q=q, t=interval_l
                ),
                maths.matrix_exponential_integration.expected_number_subs(
                    p0=p0, Q=q, t=i - step
                ),
            ]
            interval_l = 0.0
            intervals.append(interval)
    return intervals


def find_saw(infile: click.Path, dict_of_members):
    loader = get_app("load_db")
    header, rows, title, legend = load_delimited(
        infile, header=False, sep="\t"
    )  # a tsv of transit

    tbl_rows = []

    for row in tqdm(rows[1:]):  # todo:modify
        memb = dict_of_members[row[0]]
        model = loader(memb)

        t = model.lf.get_param_value(par_name="length", edge=row[1])
        q = model.lf.get_rate_matrix_for_edge(name=row[1], calibrated=True)
        p0 = model.lf.get_motif_probs()
        tau = row[3]
        model_name = model.name
        if intervals := num_saw_srch(t, q, p0, tau, model_name):
            row_of_chainsaw = row + [intervals]
            tbl_rows.append(row_of_chainsaw)
    return tbl_rows



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
