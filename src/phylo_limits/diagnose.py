from cogent3.app.result import model_result
from cogent3.app.typing import AlignedSeqsType, SerialisableType

from phylo_limits.mles_filter import mles_within_bounds


def diagonse(model_res:model_result) -> ...:
    """it will read a cogent3 model result, and judge whether there are boundary values. the parameters in a lf
    are rates, branch lengths, motif probs, but we only care about the rate boundary values 
    (as implemented in mles_filter.mles_within_bounds)  

    Args:
        model_res (model_result): _description_

    Returns:
        str: _description_
    """
    return mles_within_bounds(model_res)