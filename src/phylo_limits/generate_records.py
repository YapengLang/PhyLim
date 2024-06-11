import dataclasses
import json

from cogent3.app.composable import define_app
from cogent3.app.result import model_result

from phylo_limits import parser
from phylo_limits.check_boundary import get_bounds_violation
from phylo_limits.classify_matrix import classify_psubs
from phylo_limits.has_valid_path import IdentifiabilityCheck


clspsub_app = classify_psubs()
bound_app = get_bounds_violation()
Ident_app = IdentifiabilityCheck()


@dataclasses.dataclass(slots=True)
class PhyloLimitRec:
    """the record of phylogenetic limits"""

    source: str
    model_name: str
    boundary_values: dict
    identifiable: bool
    bad_nodes: set

    def to_rich_dict(self) -> dict:
        return {
            "source": self.source,
            "model_name": self.model_name,
            "boundary_values": self.boundary_values,
            "identifiable": self.identifiable,
            "bad_nodes": self.bad_nodes,
        }


@define_app
def generate_record(model_res: model_result, strictly=False) -> str:
    """record psubs classes, identifiability, boundary values and non-DLC projection.
    Args:
        "strictly" controls the sensitivity for Identity matrix (I); if false, treat I as DLC.
    Return:
        string in json format
    """
    tree = parser.load_tree(model_res)
    psubs = parser.load_psubs(model_res)
    params = parser.load_param_values(model_res)

    res = Ident_app(psubs=clspsub_app(psubs), tree=tree)

    rec = PhyloLimitRec(
        source=model_res.source,
        model_name=str(model_res.name),  # TODO
        boundary_values=bound_app(params).vio,
        identifiable=not bool(res),
        bad_nodes=res,
    )

    return json.dumps(rec.to_rich_dict())
