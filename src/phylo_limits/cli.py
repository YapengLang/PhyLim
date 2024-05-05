from pathlib import Path

import click

from cogent3 import get_app, open_data_store
from scitrack import CachingLogger

from phylo_limits.record import generate_record


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


@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--inpath",
    required=True,
    type=click.Path(exists=True),
    help="the directory where your alignments stored",
)
@_outpath
@click.option(
    "-a",
    "--name",
    type=str,
    required=True,
    help="the model family you want to fit",
)
@click.option(
    "-h",
    "--het",
    type=str,
    required=True,
    help="time-homo or heter",
)
def ident_check(inpath, outpath, name, het):
    # open an input directory
    dstore = open_data_store(
        inpath,
        suffix="nexus",
        mode="r",
    )

    # open an output directory
    out_dstore = open_data_store(outpath, mode="w")
    out_dstore.unlock(force=True)

    # construct the process
    loader = get_app(
        "load_aligned", format="nexus", moltype="dna"
    )  # TODO: load model_results
    # model = fit(model=name, het=het)
    recorder = generate_record()
    writer = get_app("write_db", data_store=out_dstore)

    process = loader + recorder + writer
    process.apply_to(dstore[:5], show_progress=True, parallel=False)
    print(out_dstore.describe)


if __name__ == "__main__":
    main()
