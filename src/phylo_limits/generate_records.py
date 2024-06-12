import dataclasses
import json

from typing import Union

from cogent3.app.composable import define_app
from cogent3.app.result import model_result

from phylo_limits import parser
from phylo_limits.check_boundary import get_bounds_violation
from phylo_limits.classify_matrix import MatrixCategory, classify_psubs
from phylo_limits.eval_identifiability import (
    IdentifiabilityCheck,
    IdentMsgTypes,
)


clspsub_app = classify_psubs()
bound_app = get_bounds_violation()


@dataclasses.dataclass(slots=True)
class PhyloLimitRec:
    """the record of phylogenetic limits"""

    source: str
    model_name: Union[str, None]
    boundary_values: Union[dict[str, float], None]
    mcats: dict[tuple[str], MatrixCategory]
    identifiability: bool
    strict: bool
    message: Union[set, None]
    message_type: Union[IdentMsgTypes, None]

    def to_rich_dict(self) -> dict:
        result = {
            "source": self.source,
            "model_name": "" if self.model_name is None else self.model_name,
            "boundary_values": (
                dict() if self.boundary_values is None else self.boundary_values
            ),
            "mcats": {k: v.value for k, v in self.mcats.items()},
            "identifiability": self.identifiability,
            "strict": self.strict,
            "message": {} if self.message is None else self.message,
            "message_type": self.message_type,
        }

        if info := result["message_type"]:
            result["message_type"] = info.value
        return result


@define_app
def generate_record(model_res: model_result, strict=False) -> str:
    """record psubs classes, identifiability, boundary values etc.
    Args:
        "strictly" controls the sensitivity for Identity matrix (I); if false, treat I as DLC.
    Return:
        string in json format
    """
    tree = parser.load_tree(model_res)
    psubs = parser.load_psubs(model_res)
    params = parser.load_param_values(model_res)

    Ident_app = IdentifiabilityCheck(strict=strict)
    psubs_labelled = clspsub_app(psubs)
    result = Ident_app(psubs=psubs_labelled, tree=tree)

    rec = PhyloLimitRec(
        source=model_res.source,
        model_name=model_res.name,
        boundary_values=bound_app(params).vio,
        mcats=psubs_labelled.mcats,
        identifiability=result.identifiability,
        strict=result.strict,
        message=result.message,
        message_type=result.message_type,
    )

    return json.dumps(rec.to_rich_dict())
