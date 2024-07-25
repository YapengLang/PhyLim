import dataclasses

from typing import Union

from cogent3.app.composable import define_app
from cogent3.app.result import model_result

from phylo_limits import parser
from phylo_limits.check_boundary import get_bounds_violation
from phylo_limits.classify_matrix import DLC, MatrixCategory, classify_psubs
from phylo_limits.eval_identifiability import (
    IdentCheckRes,
    IdentifiabilityCheck,
)


@dataclasses.dataclass(slots=True)
class PhyloLimitRec(IdentCheckRes):
    """the record of phylogenetic limits"""

    model_name: Union[str, None]
    boundary_values: Union[list[dict], None]
    ISCL_mcats: Union[dict[tuple[str], MatrixCategory], None]

    def to_rich_dict(self) -> dict:
        result = IdentCheckRes.to_rich_dict(self)
        result["model_name"] = self.model_name or ""
        result["boundary_values"] = self.boundary_values or []
        result["ISCL_mcats"] = {}
        if self.ISCL_mcats:
            result["ISCL_mcats"] = {k[0]: v.value for k, v in self.ISCL_mcats.items()}

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

        boundary_values = bound_app(params).vio
        psubs_labelled = clspsub_app(psubs)
        result = Ident_app(psubs_labelled, tree)

        return PhyloLimitRec(
            source=result.source,
            model_name=model_res.name,
            identifiability=result.identifiability,
            strict=result.strict,
            message=result.message,
            boundary_values=boundary_values,
            ISCL_mcats={k: v for k, v in psubs_labelled.items() if v is not DLC},
        )
