import dataclasses

from typing import Union

from cogent3.app.composable import define_app
from cogent3.app.result import model_result
from cogent3.util.table import Table

from phylo_limits.__init__ import __version__
from phylo_limits.check_boundary import (
    BoundsViolation,
    ParamRules,
    check_boundary,
)
from phylo_limits.classify_matrix import (
    DLC,
    MatrixCategory,
    ModelMatrixCategories,
    ModelPsubs,
    classify_matrix,
)
from phylo_limits.eval_identifiability import (
    IdentCheckRes,
    eval_identifiability,
)


def load_psubs(model_result: model_result) -> ModelPsubs:
    """get psubs"""
    return ModelPsubs(source=model_result.source, psubs=model_result.lf.get_all_psubs())  # type: ignore


def load_param_values(model_result: model_result) -> ParamRules:
    """get non-topology param values"""
    return ParamRules(
        source=model_result.source, params=model_result.lf.get_param_rules()  # type: ignore
    )


@define_app
class check_fit_boundary:
    """check if there are any rate params proximity to the bounds as 1e-10.
    This value is important as two clusters of fits divided by the value.
    """

    def main(self, model_result: model_result) -> BoundsViolation:
        params = load_param_values(model_result)
        return check_boundary(params)


@define_app
class classify_model_psubs:
    """labels all psubs in a given ModelPsubs object which has source info"""

    def main(self, model_result: model_result) -> ModelMatrixCategories:
        psubs = load_psubs(model_result)
        return classify_matrix(psubs)


# a rich dataclass to store bound violations, ISCL matrices, etc., besides identifiability
@dataclasses.dataclass(slots=True)
class PhyloLimitRec:
    """the record of phylogenetic limits"""

    check: IdentCheckRes
    model_name: Union[str, None]
    boundary_values: Union[list[dict], None]
    nondlc_and_identity: Union[dict[tuple[str, ...], MatrixCategory], None]

    def to_rich_dict(self) -> dict:
        result = self.check.to_rich_dict()
        result["model_name"] = self.model_name or ""
        result["boundary_values"] = self.boundary_values or []
        result["nondlc_and_identity"] = {}
        if self.nondlc_and_identity:
            result["nondlc_and_identity"] = {
                k[0]: v.value for k, v in self.nondlc_and_identity.items()
            }
        result["version"] = __version__
        return result

    @property
    def is_identifiable(self) -> bool:
        return self.check.is_identifiable

    @property
    def has_BV(self) -> bool:
        return bool(self.boundary_values)

    @property
    def violation_type(self) -> str | None:
        return None if self.is_identifiable else self.check.violation_type.name

    def to_table(self) -> Table:
        headers = [
            "source",
            "model name",
            "identifiable",
            "has boundary values",
            "version",
        ]
        rows = [
            [
                self.check.source,
                self.model_name,
                self.is_identifiable,
                self.has_BV,
                __version__,
            ]
        ]

        return Table(header=headers, data=rows, title="Phylo Limits Record")

    def _repr_html_(self) -> str:
        table = self.to_table()
        table.set_repr_policy(show_shape=False)
        return table._repr_html_()


@define_app
class phylo_limits:
    """record psubs classes, identifiability, boundary values etc of a model_result.
    Args:
        "strict" controls the sensitivity for Identity matrix (I); if false,
        treat I as DLC.
    Return:
        PhyloLimitRec object
    """

    def __init__(self, strict: bool = False) -> None:
        self.strict = strict

    def main(self, model_result: model_result) -> PhyloLimitRec:
        tree = model_result.lf.tree  # type: ignore

        check_bound_app = check_fit_boundary()
        classify_psubs_app = classify_model_psubs()

        boundary_values = check_bound_app(model_result).vio
        psubs_labelled = classify_psubs_app(model_result)
        result = eval_identifiability(psubs_labelled, tree, self.strict)

        return PhyloLimitRec(
            check=result,
            model_name=model_result.name,
            boundary_values=boundary_values,
            nondlc_and_identity={
                k: v for k, v in psubs_labelled.items() if v is not DLC
            },
        )
