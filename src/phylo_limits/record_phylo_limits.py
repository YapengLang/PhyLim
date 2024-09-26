import dataclasses

from typing import Union

from cogent3.app.composable import define_app
from cogent3.app.result import model_result

from phylo_limits.__init__ import __version__
from phylo_limits.check_boundary import ParamRules, get_bounds_violation
from phylo_limits.classify_matrix import (
    DLC,
    MatrixCategory,
    ModelPsubs,
    classify_psubs,
)
from phylo_limits.eval_identifiability import (
    IdentCheckRes,
    IdentifiabilityCheck,
)


def load_psubs(model_result: model_result) -> ModelPsubs:
    """get psubs"""
    return ModelPsubs(source=model_result.source, psubs=model_result.lf.get_all_psubs())  # type: ignore


def load_param_values(model_result: model_result) -> ParamRules:
    """get non-topology param values"""
    return ParamRules(
        source=model_result.source, params=model_result.lf.get_param_rules()  # type: ignore
    )


@dataclasses.dataclass(slots=True)
class PhyloLimitRec(IdentCheckRes):
    """the record of phylogenetic limits"""

    model_name: Union[str, None]
    boundary_values: Union[list[dict], None]
    ISCL_mcats: Union[dict[tuple[str, ...], MatrixCategory], None]

    def to_rich_dict(self) -> dict:
        result = IdentCheckRes.to_rich_dict(self)
        result["model_name"] = self.model_name or ""
        result["boundary_values"] = self.boundary_values or []
        result["ISCL_mcats"] = {}
        if self.ISCL_mcats:
            result["ISCL_mcats"] = {k[0]: v.value for k, v in self.ISCL_mcats.items()}
        result["version"] = __version__
        return result


@define_app
class generate_phylo_limit_record:
    """record psubs classes, identifiability, boundary values etc.
    Args:
        "strict" controls the sensitivity for Identity matrix (I); if false, treat I as DLC.
    Return:
        PhyloLimitRec object
    """

    def __init__(self, strict: bool = False) -> None:
        self.strict = strict

    def main(self, model_result: model_result) -> PhyloLimitRec:
        tree = model_result.lf.tree  # type: ignore
        psubs = load_psubs(model_result)
        params = load_param_values(model_result)

        clspsub_app = classify_psubs()
        bound_app = get_bounds_violation()
        Ident_app = IdentifiabilityCheck(strict=self.strict)

        boundary_values = bound_app(params).vio
        psubs_labelled = clspsub_app(psubs)
        result = Ident_app(psubs_labelled, tree)

        return PhyloLimitRec(
            source=result.source,
            model_name=model_result.name,
            identifiability=result.identifiability,
            strict=result.strict,
            message=result.message,
            boundary_values=boundary_values,
            ISCL_mcats={k: v for k, v in psubs_labelled.items() if v is not DLC},
        )


@define_app
class eval_identifiability:
    """check the identifiability of a model.
    Args:
        strict: controls the sensitivity for Identity matrix (I); if false, treat I as DLC.
    Return:
        bool value of identifiability
    """

    def __init__(self, strict: bool = False) -> None:
        self.strict = strict

    def main(self, model_result: model_result) -> bool:
        tree = model_result.lf.tree  # type: ignore
        psubs = load_psubs(model_result)

        clspsub_app = classify_psubs()
        Ident_app = IdentifiabilityCheck(strict=self.strict)

        psubs_labelled = clspsub_app(psubs)
        result = Ident_app(psubs_labelled, tree)
        return result.identifiability
