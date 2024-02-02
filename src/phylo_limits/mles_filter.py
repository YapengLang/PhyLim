import numpy
import pickle
import json
from typing import Union
from cogent3 import get_app
from cogent3.app.composable import define_app, NotCompleted
from cogent3.app.result import model_result
from cogent3.app.typing import SerialisableType, AlignedSeqsType
from cogent3.app.evo import model
from cogent3.evolve.fast_distance import DistanceMatrix
from cogent3.app.data_store import DataMember

_eps = numpy.finfo(float).eps

RATE_UPPER = 200
BRANCH_UPPER = 50
LOWER = 1e-6



# a function to check if params are within the bound
def mles_discm(rate_params: numpy.array, branch_params: numpy.array):
    """check if any parameters are close to the bound.
       return True if all well within the bound, otherwise a list consists of 1 and 0 stand for Ture or False
       (in the order of [branch_min, branch_max, rate_min, rate_max]).

    Args:
        rate_params (numpy.array): array of exchangeability terms, 0 size if JC69 or F81
        branch_params (numpy.array): array of branch lengths

    Returns:
        boolean or list
    """
    discm = numpy.array(
        [
            abs(branch_params.min() - LOWER) > _eps,
            abs(branch_params.max() - BRANCH_UPPER) > _eps,
        ]
    )
    if rate_params.size > 0:
        discm = numpy.concatenate(
            (
                discm,
                [
                    abs(rate_params.min() - LOWER) > _eps,
                    abs(rate_params.max() - RATE_UPPER) > _eps,
                ],
            )
        )
    return True if all(discm) else discm * 1


# todo, test
class SortNotCompleted:
    def __init__(self, rate_array, branch_array, result):

        self.rate = rate_array
        self.branch = branch_array

        self.result = result

    def __call__(self):
        discm = mles_discm(self.rate, self.branch)
        if isinstance(discm, numpy.ndarray):
            return NotCompleted(
                "FAIL",
                "SortNotCompleted",
                str(discm),
                source=str(self.result.source),
            )


BOUNDS = {
    "rate": {"upper": RATE_UPPER, "lower": LOWER},
    "length": {"upper": BRANCH_UPPER, "lower": LOWER},
    "motif": {"upper": 1, "lower": 0},
}


@define_app
def mles_filter(memb: DataMember) -> SerialisableType:
    """post-checking fitted model, return NotCompleted or json"""
    loader = get_app("load_db")
    result = loader(memb)
    if isinstance(result, NotCompleted):
        return result

    glob_dict = {}
    tables = result.lf.get_statistics()
    for table in tables:
        if table.title == "global params":
            glob_dict = {
                k: {
                    "upper": abs(v - BOUNDS["rate"]["upper"]),
                    "lower": abs(v - BOUNDS["rate"]["lower"]),
                }
                for k, v in table.to_dict()[0].items()
            }

        elif table.title == "edge params":
            cols = set(table.columns) - {"edge", "parent"}
            edge_dict = {}
            for k, v in table.to_dict().items():
                edge_dict[v["edge"]] = {}
                for param in cols:
                    if param == "length":
                        edge_dict[v["edge"]][param] = {
                            "upper": abs(v[param] - BOUNDS["length"]["upper"]),
                            "lower": abs(v[param] - BOUNDS["length"]["lower"]),
                        }
                    else:
                        edge_dict[v["edge"]][param] = {
                            "upper": abs(v[param] - BOUNDS["rate"]["upper"]),
                            "lower": abs(v[param] - BOUNDS["rate"]["lower"]),
                        }

        elif table.title == "motif params":
            motif_dict = {
                k: {
                    "upper": abs(v - BOUNDS["motif"]["upper"]),
                    "lower": abs(v - BOUNDS["motif"]["lower"]),
                }
                for k, v in table.to_dict()[0].items()
            }

    res_dict = (
        {
            memb.unique_id: {
                "global params": glob_dict,
                "edge params": edge_dict,
                "motif params": motif_dict,
            }
        }
        if glob_dict
        else {memb.unique_id: {"edge params": edge_dict, "motif params": motif_dict}}
    )
    return json.dumps(res_dict)


def check_proxs(p_dict):
    for v in p_dict.values():
        if set(v.keys()) == {"upper", "lower"}:
            if any(
                [v["upper"] <= 1e-10, v["lower"] <= 1e-10]
            ):  # if any proximity to bound less or equal to 1e-10
                return False
        elif not check_proxs(v):
            return False
    return True


@define_app
def mles_within_bounds(memb: DataMember) -> Union[model_result, SerialisableType]:
    """validate fitted model, return NotCompleted(if any abs(params-bound) <= 10^{-10}) or model_result"""
    loader = get_app("load_db")
    result = loader(memb)
    if isinstance(result, NotCompleted):
        return result
    import json

    # filt = mles_filter()
    # params = json.loads(filt(memb))

    # if not check_proxs(params):
    #     return NotCompleted(
    #         "FAIL",
    #         "check_proxs",
    #         "the fit is too close to param bound",
    #         source=result.source,
    #     )

    # modifies to any RATE params ..
    exclude_cols = {"edge", "parent", "length"}
    tables = result.lf.get_statistics()
    for table in tables:
        if table.title in (
            "edge params",
            "global params",
        ):  # depend on time-het setting
            rate_arr = table[
                :, [c for c in table.columns if c not in exclude_cols]
            ].array
            break

    if rate_arr.size > 0:
        discm = numpy.array(
            [
                abs(rate_arr.min() - LOWER) <= 1e-10,
                abs(rate_arr.max() - RATE_UPPER) <= 1e-10,
            ]
        )

        if any(discm):
            return NotCompleted(
                "FAIL",
                "check_proxs",
                "the fit is too close to param bound",
                source=result.source,
            )

    return result
