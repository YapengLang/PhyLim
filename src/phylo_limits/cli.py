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
from cogent3.app.result import model_result

from phylo_limits.project import project
from phylo_limits.check_ident import check_ident
from phylo_limits.diagnose import diagonse
from phylo_limits.fit import fit

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


@define_app
class PhyloLimitStat:

    def __init__(self,unique_id:str,model_name:str,identifiability:bool,boundary_values:dict,projection:dict):
        self.unique_id = unique_id
        self.model_name = model_name
        self.identifiability = identifiability
        self.boundary_values = boundary_values
        self.projection = projection


    def to_rich_dict(self) -> dict:
        return ...  

@define_app
def app_for_all_infor(model_res:model_result) -> dict:
    
    all_infor = PhyloLimitStat()
    all_infor.unique_id = ...
    all_infor.model_name = model_res.name
    all_infor.identifiability = check_ident(model_res.lf)
    all_infor.boundary_values = diagonse(model_res)
    all_infor.projection = project(model_res)
    return all_infor.to_rich_dict()

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
    required=False,
    help="the model family you want to fit",
)
def ident(inpath, outpath, name):
    # open an input directory
    dstore = open_data_store(
            inpath,
            suffix="nexus",
            mode="r",
        )  
    
    # open an output directory   
    out_dstore = open_data_store(
        outpath, mode="w"
    ) 
    out_dstore.unlock(force=True)

    # construct the process
    loader = get_app("load_aligned", format="nexus", moltype="dna") 
    model = fit(model=name)
    process = loader + model
    #TODO: all the apps will communicate with json string 




if __name__ == "__main__":
    main()
