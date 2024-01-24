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
