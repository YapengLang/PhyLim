import pytest

from phylo_limits.check_boundary import (
    BoundsViolation,
    ParamRules,
    check_boundary,
)


@pytest.mark.parametrize(
    "param_input, expected",
    [
        (
            [
                {
                    "par_name": "mprobs",
                    "init": {
                        "T": 0.1816317033849437,
                        "C": 0.25091547179313317,
                        "A": 0.24240340225709273,
                        "G": 0.3250494225648303,
                    },
                    "lower": None,
                    "upper": None,
                },
                {
                    "par_name": "C/T",
                    "edge": "159324",
                    "init": 199.9999968157517,
                    "lower": 1e-06,
                    "upper": 200,
                },
                {
                    "par_name": "A/G",
                    "edge": "878",
                    "init": 1.0000000991577895e-06,
                    "lower": 1e-06,
                    "upper": 200,
                },
                {
                    "par_name": "length",
                    "edge": "878",
                    "init": 1e-06,
                    "lower": 1e-06,
                    "upper": 50,
                },
            ],
            [
                {
                    "par_name": "A/G",
                    "edge": "878",
                    "init": 1.0000000991577895e-06,
                    "lower": 1e-06,
                    "upper": 200,
                }
            ],
        )
    ],
)
def test_get_bounds_violation(param_input, expected):
    test_input = ParamRules(source="foo", params=param_input)
    result = check_boundary(test_input)
    assert isinstance(result, BoundsViolation)
    assert result.vio == expected
