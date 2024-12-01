import dataclasses

from phylim._version import __version__


@dataclasses.dataclass(slots=True)
class ParamRules:
    source: str
    params: list[dict]


@dataclasses.dataclass(slots=True)
class BoundsViolation:
    source: str
    vio: list[dict]

    def to_rich_dict(self) -> dict:
        return {"source": self.source, "vio": self.vio, "version": __version__}


EXCLUDE_PARS = "length", "mprobs"


def check_boundary(params: ParamRules) -> BoundsViolation:
    """check if there are any rate params proximity to the bounds as 1e-10.
    This value is important as two clusters of fits divided by the value."""

    vio = []
    list_of_params = params.params
    for param in list_of_params:
        if param["par_name"] not in EXCLUDE_PARS:
            if (abs(param["init"] - param["lower"]) <= 1e-10) or (
                abs(param["init"] - param["upper"]) <= 1e-10
            ):
                vio.append(param)
    return BoundsViolation(source=params.source, vio=vio)
