import dataclasses

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


@dataclasses.dataclass(slots=True)
class PhyloLimitRec:
    """the record of phylogenetic limits"""

    source: str
    model_name: Union[str, None]
    boundary_values: Union[list[dict], None]
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
                [] if self.boundary_values is None else self.boundary_values
            ),
            "mcats": {k[0]: v.value for k, v in self.mcats.items()},
            "identifiability": self.identifiability,
            "strict": self.strict,
            "message": [] if self.message is None else list(self.message),
            "message_type": self.message_type,
        }

        if info := result["message_type"]:
            result["message_type"] = info.value
        return result


@define_app
class generate_record:
    """record psubs classes, identifiability, boundary values etc.
    Args:
        "strict" controls the sensitivity for Identity matrix (I); if false, treat I as DLC.
    Return:
        string in json format
    """

    def __init__(self, strict=False) -> None:
        self.strict = strict

    def main(self, model_res: model_result) -> PhyloLimitRec:
        tree = parser.load_tree(model_res)
        psubs = parser.load_psubs(model_res)
        params = parser.load_param_values(model_res)
        clspsub_app = classify_psubs()
        bound_app = get_bounds_violation()

        Ident_app = IdentifiabilityCheck(strict=self.strict)
        psubs_labelled = clspsub_app(psubs)
        result = Ident_app(psubs_labelled, tree)

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
        return rec
