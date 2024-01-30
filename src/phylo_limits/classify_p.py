import numpy
from cogent3.util.dict_array import DictArray
from cogent3.app.result import model_result
from phylo_limits import check_dlc, delta_col

import warnings

#TODO: this should get a cogent3 likelihood function as input.. the arguments are unnecessary. 
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

        if check_dlc.check_I(value):
            if strictly == True:
                new_dict[key] = {"value": value, "class": "Identity"}
            else:
                warnings.warn(
                    "Fit contains Identity matrix, which makes topology un-identifiable!"
                )
                new_dict[key] = {"value": value, "class": "DLC"}
            continue

        elif check_dlc.check_dlc(p_matrix=value, p_limit=numpy.array([pi, pi, pi, pi])):
            new_dict[key] = {"value": value, "class": "DLC"}
            continue

        elif check_dlc.check_chainsaw(
            p_matrix=value, p_limit=numpy.array([pi, pi, pi, pi])
        ):

            new_dict[key] = {"value": value, "class": "Chainsaw"}
            continue

        else:
            if check_dlc.check_limit(
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
