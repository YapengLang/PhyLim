from pathlib import Path

import click

from scitrack import CachingLogger

import glob
import numpy
import scipy
from tqdm import tqdm
from cogent3 import (
    get_app,
    open_data_store,
    make_table,
    maths,
    make_aligned_seqs,
    load_aligned_seqs,
)
from cogent3.parse.table import load_delimited
from cogent3.app.typing import SerialisableType, AlignedSeqsType
from cogent3.app.composable import define_app, NotCompleted
from cogent3.app.data_store import DataMember

from phylo_limits.delta_col import get_sub_num
from phylo_limits.mles_filter import (
    make_model_app,
    mles_within_bounds,
    fit_core_E,
    fit_core_W,
    mles_filter,
)  # import filter
from phylo_limits.fit_stats import get_bound_viol, get_dlc_viol, get_transit
from phylo_limits.check_dlc import check_chainsaw
from phylo_limits.check_ident import check_ident, check_all_psubs
from phylo_limits.find_saw import find_saw
from phylo_limits.parse_bad_phylip import load_bad_phylip

from contextlib import redirect_stderr
import io

__author__ = "Yapeng Lang"
__copyright__ = "Copyright 2016-2023, Yapeng Lang"
__credits__ = ["Yapeng Lang"]
__license__ = "BSD"
__version__ = "2023.3.31"  # A DATE BASED VERSION
__maintainer__ = "Yapeng Lang"
__email__ = "Yapeng.Lang@anu.edu.au"
__status__ = "alpha"


LOGGER = CachingLogger()


@click.group()
@click.version_option(__version__)  # add version option
def main():
    """cl app for model fitting and result filtering"""
    pass


_verbose = click.option(
    "-v",
    "--verbose",
    count=True,
    help="is an integer indicating number of cl occurrences",
)


_outpath = click.option(
    "-o", "--outpath", type=Path, help="the input string will be cast to Path instance"
)


if __name__ == "__main__":
    main()
