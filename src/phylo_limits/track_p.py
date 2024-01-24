from cogent3 import open_data_store, get_app, make_table, load_delimited
from phylo_limits.delta_col import tau_TX
from scipy.linalg import expm
import json

LOADER = get_app("load_db")


def load_chainsaws(M, where_ens):
    chainsaws = load_delimited(
        f"~/repos/YapengAnalysis/data/model_stat_microbial_18Jul_fixed/report_microbial_25Jul/chainsaws/{M}.tsv",
        header=False,
        sep="\t",
    )[1][1:]
    for row in chainsaws:
        intervals = json.loads(row[5])
        for interval in intervals:
            if (
                abs(where_ens - interval[0]) < 1e-5
            ):  # choose a proper ens which met the expected
                return row[0], row[1], interval[0]


def load_ndlcs(M, where_ens):
    ndlcs = load_delimited(
        f"~/repos/YapengAnalysis/data/model_stat_microbial_18Jul_fixed/report_microbial_25Jul/transit/{M}.tsv",
        header=False,
        sep="\t",
    )[1][1:]
    for row in ndlcs:
        if (
            abs(float(row[2]) - where_ens) < 1e-5
        ):  # choose a proper ens which met the expected
            return row[0], row[1], float(row[2])


def load_model_lf(ds, algn, edge):
    for i in ds.completed:
        if i.unique_id == algn:
            model = LOADER(i)
            q = model.lf.get_rate_matrix_for_edge(name=edge, calibrated=False)
            p0 = model.lf.get_motif_probs()
            return q, p0, model


def teleport(
    M, want_ndlc=False, want_chainsaw=False, where_ens=None, algn=None, edge=None
):
    """if algn and edge provided, then return p for given ens. Or choose a ndlc around ens provided"""
    if want_ndlc and want_chainsaw:
        raise ValueError("plz give right expectation")
    ds = open_data_store(
        f"~/repos/YapengAnalysis/data/model_results_microbial_2Jun/model_results_{M}_microBen2015.sqlitedb",
        mode="r",
    )
    if want_chainsaw:
        algn, edge, where_ens = load_chainsaws(M, where_ens)
    if want_ndlc:
        algn, edge, where_ens = load_ndlcs(M, where_ens)

    q, p0, model = load_model_lf(ds=ds, algn=algn, edge=edge)
    tau = tau_TX(p0=p0, q=q, ens=where_ens)
    p = expm(q.to_array() * tau)
    return algn, edge, p0, p, where_ens, tau, model
