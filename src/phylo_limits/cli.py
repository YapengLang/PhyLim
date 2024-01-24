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
from phylo_limits.entro_meas import D_entropy_t_stat, only_varsites, de_gap
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


@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    default="/Users/yapenglang/repos/YapengAnalysis/data/doi_10.5061_dryad.g7g0n__v1/microbial",
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
@click.option(
    "-a",
    "--name",
    type=str,
    required=False,
    help="the model family you want to fit",
)
@click.option(
    "-A", "--allmodel", required=False, is_flag=True, help="auto fit all models"
)
@click.option(
    "-V",
    "--varonly",
    required=False,
    is_flag=True,
    help="only fit to infor sites as per Duchene 2022",
)
@click.option(
    "-D",
    "--degap",
    required=False,
    is_flag=True,
    help="only fit to degapped",
)
@click.option(
    "-P",
    "--purefit",
    required=False,
    is_flag=True,
    help="filter out bound violation fits",
)
def fit_micro(infile, outpath, overwrite, name, allmodel, varonly, degap, purefit):
    names = {"GTR", "ssGN", "GN"} if allmodel else {name}
    for n in names:
        name = n
        outfile = (
            f"{outpath}/model_result_{name}_varonly.sqlitedb"
            if varonly
            else f"{outpath}/model_result_{name}.sqlitedb"
        )
        # open data store for in- and out-put
        if overwrite:
            outfile = outpath
            outfile.unlink()
        out_dstore = open_data_store(
            outfile, mode="w"
        )  # specify a database to store the results
        out_dstore.unlock(force=True)

        # open ds where appiled to
        dstore = open_data_store(
            infile,
            suffix="nexus",
            mode="r",
        )  # open an input directory

        # construct the process
        loader = get_app("load_aligned", format="nexus", moltype="dna")
        model = make_model_app(name)
        param_filter = mles_within_bounds()
        writer = get_app("write_db", data_store=out_dstore)

        if varonly:
            raise KeyError("not yet been developed!")
            # get_var = only_varsites()
            # process = loader + get_var + model + param_filter + writer
        elif degap:
            no_gap = de_gap()
            process = loader + no_gap + model + param_filter + writer
        else:
            if purefit:
                process = loader + model + param_filter + writer
            else:
                process = loader + model + writer

        process.apply_to(
            dstore[:], show_progress=True, parallel=False
        )  # fine to change the selected number
        print(out_dstore.describe)


@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
@click.option(
    "-N",
    "--nprocs",
    type=int,
    default=1,
    help="number of processes",
)
@click.option(
    "-L",
    "--limit",
    type=int,
    default=0,
    help="number of data records to process",
)
@click.option(
    "-A", "--allmodel", required=False, is_flag=True, help="auto fit all models"
)
def calc_subs(infile, outpath, overwrite, nprocs, limit, allmodel):
    if allmodel:
        names = {"ssGN", "GN", "GTR"}
    # for n in names:
    # open data store for in- and output
    # if len(glob.glob(f"{infile}/model*result{n}*.sqlitedb")) == 1:
    #     dstore_path = glob.glob(f"{infile}/model*result*{n}*.sqlitedb")[0]
    # else:
    #     for path in glob.glob(f"{infile}/model*result*{n}*.sqlitedb"):
    #         if "ss" not in path:
    #             dstore_path = path
    dstore_path = infile
    dstore = open_data_store(dstore_path, mode="r", limit=limit or None)
    if overwrite:
        outpath.unlink()
    out_dstore = open_data_store(outpath, mode="w")
    out_dstore.unlock(force=True)

    # construct a calculating process
    # loader = get_app("load_db")
    filt = mles_within_bounds()
    calculator = get_sub_num()
    writer = get_app("write_db", data_store=out_dstore)
    # process = loader + calculator + writer
    process = filt + calculator + writer
    parallel = nprocs > 1
    process.apply_to(
        dstore,
        show_progress=True,
        parallel=parallel,
        par_kw=dict(max_workers=nprocs, chunksize=len(dstore) // nprocs),
    )
    print(out_dstore.describe)


@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
def get_notcomp(infile, outpath, overwrite):
    dstore = open_data_store(
        infile,
        mode="r",
    )  # open an input directory
    if overwrite:
        outpath.unlink()

    # construct the process
    loader = get_app("load_db")
    NotComps = dstore.not_completed
    notcomps = [loader(i) for i in NotComps]
    table = get_bound_viol(notcomps)
    table.write(outpath, format="tsv")


@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "-A", "--allmodel", required=False, is_flag=True, help="auto fit all models"
)
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
def get_ens(infile, outpath, overwrite, allmodel):
    # if allmodel:
    #     names = {"JC69", "F81", "HKY85", "TN93", "ssGN", "K80", "GN", "GTR"}
    # for n in names:
    # if len(glob.glob(f"{infile}/model*result{n}*.sqlitedb")) == 1:
    #     dstore_path = glob.glob(f"{infile}/model*result{n}*.sqlitedb")[0]
    # elif not glob.glob(f"{infile}/model*result{n}*.sqlitedb"):
    #     continue
    # else:
    #     for path in glob.glob(f"{infile}/model*result{n}*.sqlitedb"):
    #         if "ss" not in path:
    #             dstore_path = path
    dstore_path = infile
    # open data store
    dstore = open_data_store(
        dstore_path,
        mode="r",
    )
    if overwrite:
        outpath.unlink()

    # construct the process
    rows = []
    loader = get_app("load_db")
    # filt = mles_within_bounds()
    finder = get_dlc_viol()
    process = loader + finder
    # process = filt + finder

    for memb in tqdm(dstore.completed):
        algn_list = process(memb)
        if not isinstance(
            algn_list, NotCompleted
        ):  # some cases raise ValueError in limit matrix checking
            rows.extend(iter([[memb.unique_id] + row for row in algn_list]))

    header = (
        ["algn name"]
        + ["edge"]
        + ["length"]
        + ["ens"]
        + ["delta_col"]
        + ["sympathetic"]
        + ["chainsaw"]
    )

    table = make_table(header, data=rows)
    table.write(outpath, format="tsv")


@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
def find_transit(infile, outpath, overwrite):
    # open data store for stats tables
    dstore = open_data_store(
        infile,
        mode="r",
    )
    if overwrite:
        outpath.unlink()

    # construct the process of get true transition point
    rows = []
    loader = get_app("load_db")
    finder = get_transit()
    process = loader + finder

    for memb in dstore.completed:
        algn_list = process(memb)
        algn_list = [[memb.unique_id] + n for n in algn_list]
        rows.extend(iter(algn_list))

    header = (
        ["id"] + ["edge name"] + ["substitution number"] + ["tau_"] + ["delta column"]
    )

    table = make_table(header, data=rows)
    table.write(outpath, format="tsv")


# prototype of hunting chainsaw beyond the dlc violation point
@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
@click.option("-k", "--key", required=True, help="key is model name")
def dirty_saw(outpath, infile, overwrite, key):

    # open data store
    dstore = open_data_store(
        key,  # f"~/repos/YapengAnalysis/data/model_results_microbial_2Jun/model_results_{key}_microBen2015.sqlitedb",  # TODO: modify
        mode="r",
    )
    if overwrite:
        outpath.unlink()

    members = dstore.completed

    dict_of_members = {memb.unique_id: memb for memb in members}
    tbl_rows = find_saw(infile, dict_of_members)

    headers = (
        ["id"]
        + ["edge name"]
        + ["substitution number"]
        + ["tau_"]
        + ["delta column"]
        + ["interval"]
    )

    table = make_table(headers, data=tbl_rows)
    table.write(outpath, format="tsv")


@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    default="/Users/yapenglang/repos/YapengAnalysis/data/doi_10.5061_dryad.g7g0n__v1/microbial",
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
def calc_entro_D(outpath, infile, overwrite):
    """line to calc Duchene's entropy measure for an alignment"""
    ds = open_data_store(
        infile,
        suffix="nexus",
        mode="r",
    )  # open an input directory
    if overwrite:
        outpath.unlink()

    # construct the process
    loader = get_app("load_aligned", format="nexus", moltype="dna")
    calc_all = D_entropy_t_stat(only_varsites=False)
    calc_infor = D_entropy_t_stat(only_varsites=True)

    # process_all = loader + calc_all
    # process_infor = loader + calc_infor

    rows = []
    for memb in tqdm(ds.completed[:]):
        algn = loader(memb)
        result_all = calc_all(algn)
        result_infor = calc_infor(algn)
        row = [memb.unique_id] + [
            result_all[0],
            result_infor[0],
            result_all[1],
            result_infor[1],
        ]
        rows.append(row)

    header = (
        ["id"]
        + ["tstat for all sites"]
        + ["for only infor sites"]
        + ["pvalue for all"]
        + ["pvalue for infor"]
    )

    table = make_table(header, data=rows)
    table.write(outpath, format="tsv")


# TODO: delete
@define_app
def get_aligned_Espeland(
    dsmemb: DataMember, only_algn=False
) -> SerialisableType | AlignedSeqsType:
    """tackle with fasta format in Espeland et al, to get rid of space"""
    from cogent3.parse.fasta import MinimalFastaParser

    file = f"{str(dsmemb.data_store.source)}/{dsmemb.unique_id}"
    f = open(file)
    seqs = list(MinimalFastaParser(f))
    seqs = [
        (row[0], row[1].replace(" ", "").upper()) for row in seqs
    ]  # index long name to keep short
    algn = make_aligned_seqs(seqs, moltype="dna", source=dsmemb.unique_id[:-6])
    if only_algn:
        return algn
    try:
        matrix = algn.distance_matrix(calc="paralinear", drop_invalid=False)
    except ArithmeticError:
        print(f"{file} cannot be calced paralinear")

    return matrix.to_dict()


# TODO: modify the syntax to be general
@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--inpath",
    type=str,
    default="/Users/yapenglang/repos/YapengAnalysis/data/data_duchene2022/empiricalDataLoci/Espeland_2018/loci",
    help="data repo",
)
@click.option(
    "-o",
    "--outpath",
    type=str,
    default="/Users/yapenglang/repos/YapengAnalysis/data/data_duchene2022/empiricalDataLoci/Espeland_2018/paralinear_tree_8Sep/",
    help="data outpath",
)
def calc_para_Espeland(inpath, outpath):
    import pickle

    calc = get_aligned_Espeland()
    pathlist = Path(inpath).rglob("*.fasta")
    for path in pathlist:
        path_in_str = str(path)
        matrix = calc(path_in_str)
        with open(
            outpath + path_in_str.split("/")[-1][:-6] + ".pickle", "wb"
        ) as handle:
            pickle.dump(matrix, handle, protocol=pickle.HIGHEST_PROTOCOL)


@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    default="/Users/yapenglang/repos/YapengAnalysis/data/data_duchene2022/empiricalDataLoci/Espeland_2018/loci/",
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
@click.option(
    "-a",
    "--name",
    type=str,
    required=False,
    help="the model family you want to fit",
)
@click.option(
    "-A", "--allmodel", required=False, is_flag=True, help="auto fit all models"
)
@click.option(
    "-f",
    "--format",
    required=True,
    help="infile format indicated as suffix ",
    default="fasta",
)
def fit_timehomo(infile, outpath, overwrite, name, allmodel, format):
    names = {"GTR", "ssGN", "GN"} if allmodel else {name}
    formats = {
        "fa": "fasta",
        "fasta": "fasta",
        "nex": "nexus",
        "nexus": "nexus",
        "phylip": "phylip",
        "phy": "phylip",
    }
    for n in names:
        outfile = f"{outpath}/model_result_{n}.sqlitedb"
        # open data store for in- and out-put
        if overwrite:
            outfile = outpath
            outfile.unlink()
        out_dstore = open_data_store(
            outfile, mode="w"
        )  # specify a database to store the results
        out_dstore.unlock(force=True)

        # open ds where appiled to
        in_dstore = open_data_store(
            infile,
            suffix=format,
            mode="r",
        )  # open an input directory

        # construct the process
        dist_cal = get_app("fast_slow_dist", fast_calc="paralinear", moltype="dna")
        est_tree = get_app("quick_tree", drop_invalid=False)
        calc_tree = dist_cal + est_tree
        if formats[format] == "phylip":
            loader = load_bad_phylip()
        else:
            loader = get_app("load_aligned", moltype="dna", format=formats[format])
        modeller = get_app(
            "model",
            n,
            tree_func=calc_tree,
            upper=200,
            lower=1e-6,
            param_rules=[
                {
                    "par_name": "length",
                    "upper": 50,
                }
            ],
        )
        writer = get_app("write_db", data_store=out_dstore)
        app = loader + modeller + writer
        out_dstore = app.apply_to(
            in_dstore[:], show_progress=True
        )  # fine to change num
        print(out_dstore.describe)


@define_app
def pick_notcomp(dsmemb: DataMember) -> AlignedSeqsType:
    loader = get_app("load_db")
    uid = f"{dsmemb.unique_id}.nexus"
    if isinstance(loader(dsmemb), NotCompleted):
        return load_aligned_seqs(
            f"/Users/yapenglang/repos/YapengAnalysis/data/doi_10.5061_dryad.g7g0n__v1/microbial/{uid}",
            format="nexus",
            moltype="dna",
        )
    else:
        return NotCompleted(
            "FAIL",
            "pick_notcomp",
            "Done previously",
            source=uid,
        )


@main.command(no_args_is_help=True)
@_outpath
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
def fit_notcomp_cases(outpath, overwrite):  # infile is the existing model results

    names = {"GTR", "ssGN", "GN"}
    for name in names:
        outfile = f"{outpath}/model_result_{name}_boundvio.sqlitedb"
        inpath = f"/Users/yapenglang/repos/YapengAnalysis/data/model_results_microbial_2Jun/model_results_{name}_microBen2015.sqlitedb"
        # open data store for in- and out-put
        if overwrite:
            outfile = outpath
            outfile.unlink()
        out_dstore = open_data_store(
            outfile, mode="w"
        )  # specify a database to store the results
        out_dstore.unlock(force=True)

        dstore = open_data_store(
            inpath,
            mode="r",
        )

        # construct the process
        pick = pick_notcomp()
        model = make_model_app(name)
        writer = get_app("write_db", data_store=out_dstore)
        process = pick + model + writer

        process.apply_to(
            dstore[:], show_progress=True, parallel=False
        )  # fine to change the selected number
        print(out_dstore.describe)


@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
@click.option(
    "-a",
    "--name",
    type=str,
    required=False,
    help="",
)
def get_bounds(infile, outpath, overwrite, name):
    # names = {"GTR", "ssGN", "GN"}

    outfile = f"{outpath}/{name}_boundvio.sqlitedb"
    dstore_path = infile
    # open data store for in- and out-put
    if overwrite:
        outfile = outpath
        outfile.unlink()
    out_dstore = open_data_store(
        outfile, mode="w"
    )  # specify a database to store the results
    out_dstore.unlock(force=True)

    dstore = open_data_store(
        dstore_path,
        mode="r",
    )
    # process
    filt = mles_filter()
    writer = get_app("write_db", data_store=out_dstore)
    process = filt + writer

    process.apply_to(
        dstore[:], show_progress=True, parallel=False
    )  # fine to change the selected number
    print(out_dstore.describe)


@define_app
def tbl_breaking(memb: DataMember) -> SerialisableType:
    """get a memb, run ident. check on it. return a json with each row for a breaking node"""
    loader = get_app("load_db")
    model = loader(memb)
    uid = memb.unique_id
    nodes = check_ident(model.lf, strictly=False, which=True)
    rows = []
    for node in nodes:
        row = [uid, node]
        rows.append(row)
    header = ["id"] + ["bad node"]

    table = make_table(header, data=rows)
    return table.to_json()


@define_app
def tbl_ISCL(memb: DataMember) -> SerialisableType:
    """get a memb, run check all psubs. return a json with each row for a interesting node"""
    loader = get_app("load_db")
    model = loader(memb)
    uid = memb.unique_id
    lf = model.lf
    labelled_psubs = check_all_psubs(
        lf.get_all_psubs(),
        lf.name,
        lf.get_motif_probs_by_node(),
        strictly=True,
        label_L=True,
    )
    rows = []
    for k, v in labelled_psubs.items():
        if v["class"] != "DLC":
            row = [uid, k[0], v["class"]]
            rows.append(row)
    header = ["id"] + ["node"] + ["class"]

    table = make_table(header, data=rows)
    return table.to_json()


@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
@click.option(
    "-a",
    "--name",
    type=str,
    required=False,
    help="the model family you want to fit",
)
@click.option(
    "-A", "--allmodel", required=False, is_flag=True, help="auto fit all models"
)
@_outpath
def run_ident_check(infile, overwrite, allmodel, name, outpath):
    """a rich ident. check for generating two ds for one data set;
    one ds contains id breaking nodes, one ds contains all types of matrix other than dlc"""
    import os

    # os.mkdir(f"{infile}/identifiability")
    # os.mkdir(f"{infile}/identifiability/breaking")
    # os.mkdir(f"{infile}/identifiability/ISCL_psubs")
    names = {"GTR", "ssGN", "GN"} if allmodel else {name}

    for n in names:
        outfile_breaking = f"{outpath}/identifiability/breaking/{n}.sqlitedb"
        outfile_ISCL_psubs = f"{outpath}/identifiability/ISCL_psubs/{n}.sqlitedb"

        # open data store for in- and out-put
        if overwrite:
            outfile_breaking.unlink()
            outfile_ISCL_psubs.unlink()
        out_dstore_breaking = open_data_store(
            outfile_breaking,
            mode="w",
        )
        out_dstore_ISCL_psubs = open_data_store(
            outfile_ISCL_psubs,
            mode="w",
        )
        # open ds where appiled to
        # if len(glob.glob(f"{infile}/model*result{n}*.sqlitedb")) == 1:
        #     dstore_path = glob.glob(f"{infile}/model*result{n}*.sqlitedb")[0]
        # else:
        #     for path in glob.glob(f"{infile}/model*result{n}*.sqlitedb"):
        #         if "ss" not in path:
        #             dstore_path = path
        dstore_path = f"{infile}/model_result_{n}.sqlitedb" if allmodel else infile
        in_dstore = open_data_store(
            dstore_path,
            mode="r",
        )  # open an input directory

        # construct the process
        breaking = tbl_breaking()
        writer = get_app("write_db", data_store=out_dstore_breaking)
        process = breaking + writer
        process.apply_to(in_dstore, show_progress=True)
        print(out_dstore_breaking.describe)

        iscl = tbl_ISCL()
        writer = get_app("write_db", data_store=out_dstore_ISCL_psubs)
        process = iscl + writer
        process.apply_to(in_dstore, show_progress=True)
        print(out_dstore_ISCL_psubs.describe)


if __name__ == "__main__":
    main()
