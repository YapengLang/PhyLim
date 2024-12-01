"""A library for checking the limits of phylogenetic tree estimation.
"""

from cogent3.app.composable import define_app
from cogent3.app.result import model_result

from phylim.check_boundary import BoundsViolation, ParamRules, check_boundary
from phylim.classify_matrix import (
    ModelMatrixCategories,
    ModelPsubs,
    classify_matrix,
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
