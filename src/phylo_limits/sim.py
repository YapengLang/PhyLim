from cogent3 import (
    get_app,
    available_models,
    load_aligned_seqs,
    make_table,
    get_model,
    make_tree,
)
from cogent3.core.alignment import ArrayAlignment
from cogent3.evolve.parameter_controller import AlignmentLikelihoodFunction
from numpy import linalg
import numpy
from numpy import matmul as mm
from numpy import transpose as tt

ii = linalg.inv
BH = get_model("BH")


def make_psub(M):
    pdict = {}
    row = {}
    nucs = ["T", "C", "A", "G"]
    for i in [0, 1, 2, 3]:
        for j in [0, 1, 2, 3]:
            row[nucs[j]] = M[i, j]
        pdict[nucs[i]] = row
        row = {}
    return pdict


def apply_rules(psubs: dict, tips=None, motif=None):
    if tips is None:
        tips = ("a", "b", "c")
    if motif is None:
        motif = {"T": 0.25, "C": 0.25, "A": 0.25, "G": 0.25}
    rules_dict = []
    chrac = "("
    for i in range(len(tips)):
        chrac += f"{tips[i]})" if i + 1 == len(tips) else f"{tips[i]}, "
    tree = make_tree(chrac)
    param_rules = BH.make_likelihood_function(tree).get_param_rules()
    x = 0
    for rule in param_rules:
        if rule["par_name"] == "psubs":
            rule["edge"] = tips[x]
            rule["init"] = make_psub(psubs[tips[x]])
            x += 1
        if rule["par_name"] == "mprobs":
            rule["init"] = motif
        rules_dict.append(rule)

    lf = BH.make_likelihood_function(tree=tree)
    lf.apply_param_rules(rules_dict)
    return lf


def calc_G(Pa, Pb, Pc, Pim=numpy.array([0.25, 0.25, 0.25, 0.25]), gamma=0):
    """func to calc G in Chang, which uniquely determine the joint distribution on tri-terminal if DLC"""
    assert (
        abs(linalg.det(Pa) - 0) > 1e-5
        and abs(linalg.det(Pb) - 0) > 1e-5
        and abs(linalg.det(Pc) - 0) > 1e-5
    )
    Pia = Pim @ Pa
    Pib = Pim @ Pb
    Pic = Pim @ Pc

    Big_pia = numpy.diag(Pia)
    Big_pim = numpy.diag(Pim)
    Pam = ii(Big_pia) @ tt(Pa) @ Big_pim  # (1) in Chang
    Pab = Pam @ Pb
    Pab_gamma = Pam @ numpy.diag(Pc[:, gamma]) @ Pb
    return ii(Pab) @ Pab_gamma


def nor(matrix):
    sum_r = matrix.sum(axis=1)
    rows = []
    for i in [0, 1, 2, 3]:
        row = matrix[i] / sum_r[i]
        rows.append(row)
    return numpy.array(rows)


def simulator(model_name, tree, rules) -> AlignmentLikelihoodFunction:
    """apply a set of model, to generate sim data"""
    model = get_model(model_name)

    lf = model.make_likelihood_function(tree=tree)
    lf.apply_param_rules(rules)
    return lf
