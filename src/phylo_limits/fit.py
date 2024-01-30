import numpy

from cogent3 import get_app
from cogent3.app.result import model_result
from cogent3.app.typing import AlignedSeqsType, SerialisableType


# parameters of model fitting
RATE_UPPER = 200
BRANCH_UPPER = 50
LOWER = 1e-6


def model_init(
    model, het:bool, tree_func = None, opt_args=None
):  
    """input: model you want to fit, and het setting, and tree (if needed)
        output: an cogent3 app
    """

    opt_args = opt_args or {}
    opt_args = {"tolerance": 1e-8, **opt_args}

    param_rules = [
        {
            "par_name": "length",
            "upper": BRANCH_UPPER,
        }
    ]

    if het == True:
        time_het = "max" if model not in {"JC69", "F81"} else None
    else:
        time_het = None

    if tree_func: # TODO: consider : 1) apply an interface to provide tree? 2) report invalid paralinear distance? 
        unique_trees = False
    else:
        unique_trees = True

    return get_app(
        "model",
        model,
        tree_func=tree_func,
        opt_args=opt_args,
        param_rules=param_rules,
        upper=RATE_UPPER,
        lower=LOWER,
        optimise_motif_probs=True,
        show_progress=False,
        time_het=time_het,
        unique_trees=unique_trees,
    )



def fit(algn:AlignedSeqsType, model:str, het=False) -> model_result:
    """input: a cogent3 alignment
        output: a cogent3 model result object 

        arguments: model: run evo on this model, het: time-het setting
        des: given an algn, run specified model on it, then return model result.
    """
    if algn.num_seqs < 4:
        tree_func = None
    else:
        dist_cal = get_app("fast_slow_dist", fast_calc="paralinear", moltype="dna")
        est_tree = get_app("quick_tree", drop_invalid=False)
        tree_func = dist_cal + est_tree
    
    modeller = model_init(model = model, het = het, tree_func = tree_func)
    return modeller(algn)



#TODO: # and also feature of fitting all alignments as a data store for command line using. 