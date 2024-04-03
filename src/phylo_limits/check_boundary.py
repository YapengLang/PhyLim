import numpy

from cogent3.evolve.parameter_controller import AlignmentLikelihoodFunction


# TODO: read from lf function


def mles_within_bounds(params, bounds) -> dict:
    """check if there are any rate params proximity to the bounds as 1e-10"""
    exclude_cols = {"edge", "parent", "length"}
    tables = params

    for table in tables:
        if table.title in (
            "edge params",
            "global params",
        ):  # depend on time-het setting
            rate_arr = table[
                :, [c for c in table.columns if c not in exclude_cols]
            ].array
            break

    if rate_arr.size > 0:
        discm = numpy.array(
            [
                abs(rate_arr.min() - bounds["rate_lower"]) <= 1e-10,
                abs(rate_arr.max() - bounds["rate_upper"]) <= 1e-10,
            ]
        )

        if any(discm):
            return {"rate min": rate_arr.min(), "rate max": rate_arr.max()}

    return {}


def diagonse(statistics, bounds) -> dict:
    """it will read a cogent3 likelihood function, then judge whether there are boundary values. the parameters in a lf
    are rates, branch lengths, motif probs, but we only care about the rate boundary values
    (as implemented in mles_filter.mles_within_bounds).

    Args:
        lf(AlignmentLikelihoodFunction)

    Returns:
        the maximum and minimum values of rate params if any of them close to the bounds as 1e-10
    """

    return mles_within_bounds(params=statistics, bounds=bounds)
