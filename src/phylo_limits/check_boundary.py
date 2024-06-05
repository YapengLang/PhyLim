import dataclasses

from cogent3.app.composable import define_app


@dataclasses.dataclass(slots=True)
class ParamRules:
    source: str
    params: list[dict]


@dataclasses.dataclass(slots=True)
class BoundsViolation:
    source: str
    vio: dict[str, float]

    def items(self):
        return self.vio.items()


@define_app
class get_bounds_violation:
    """check if there are any rate params proximity to the bounds as 1e-10.
    This value is important as two clusters of fits splitted by the value."""

    exclude_params = "length", "mprobs"

    def __init__(self) -> None:
        pass

    def main(self, params: ParamRules) -> BoundsViolation:
        vio = dict()
        list_of_params = params.params
        for param in list_of_params:
            if param["par_name"] not in self.exclude_params:
                if (abs(param["init"] - param["lower"]) <= 1e-10) or (
                    abs(param["init"] - param["upper"]) <= 1e-10
                ):
                    vio[param["par_name"]] = param["init"]
        return BoundsViolation(source=params.source, vio=vio)
