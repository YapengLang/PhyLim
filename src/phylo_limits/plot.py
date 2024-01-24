# this module tends to extract infor from analysis
import json, numpy
import plotly.express as px
from cogent3 import load_delimited, open_data_store, get_app
from phylo_limits.entro_meas import de_gap
import ast
from tqdm import tqdm

COL_SEQ = px.colors.qualitative.D3

REPO = "/Users/yapenglang/repos/YapengAnalysis/data/doi_10.5061_dryad.g7g0n__v1/microbial"  # MICROBIAL_REPO

DS = open_data_store(
    REPO,
    suffix="nexus",
    mode="r",
)

ALGN_NAMES = [memb.unique_id for memb in DS.members]


class colours:
    def __init__(self):
        self.GTR = COL_SEQ[9]
        self.GN = "#8316B5"
        self.ssGN = "#F5802A"
        self.theme = COL_SEQ[0]
        self.JC69 = COL_SEQ[0]
        self.TN93 = COL_SEQ[2]
        self.K80 = COL_SEQ[3]
        self.F81 = COL_SEQ[7]
        self.HKY85 = COL_SEQ[5]
        self.red = "#D62728"


class notations:
    def __init__(self):
        self.coldel = r"$\LARGE\Delta_{\text{column}}$"
        self.mincoldel = r"$\LARGE\min(\Delta_{\text{column}})$"
        self.enstau_ = r"$\LARGE\mathrm{ENS}(\tau_{-})$"
        self.ens=r"$\LARGE\mathrm{ENS}$"
        self.tobs = r"$\LARGE\mathrm{t}_{obs}$"
        self.tvobs = r"$\LARGE\mathrm{t}_{var\text{_}obs}$"


def save_fig(fig, name, s=5, h=400, w=900):
    """save plotly fig into thesis repo, given name and params"""
    import plotly.io as pio

    pio.write_image(
        fig,
        f"/Users/yapenglang/Documents/master_project/thesis/figs/{name}.png",
        scale=s,
        height=h,
        width=w,
    )


def show_latex():
    from IPython.display import display, HTML
    import plotly

    ## Tomas Mazak's workaround
    plotly.offline.init_notebook_mode()
    display(
        HTML(
            '<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_SVG"></script>'
        )
    )


def set_font(fig):
    """set global font in a px fig object (not worked for embedded latex), auto show latex in jupyter"""
    show_latex()
    fig.update_layout(
        font=dict(
            size=25,
        ),
        font_family="Serif",
        template="seaborn"
    )


def odessy(all_data: list):
    """order sets by an intersection, in order to plot data
    assumption: input are {unique_id (algn) : datum}"""
    name_set = set()
    for data in all_data:
        name_set = name_set & data.keys() if name_set else data.keys()

    return [{name: data[name] for name in name_set} for data in all_data]


def get_deltas(infile):
    """load a tsv file including ens, tau, algn_name, delta all branches etc... It will return the all delta values in an algn"""
    rows = load_delimited(infile, header=False, sep="\t",)[
        1
    ][1:]

    dict_for_deltas = {}
    for row in rows:
        names = row[0]
        if names in dict_for_deltas:
            dict_for_deltas[names].append(row[-3])
        else:
            dict_for_deltas[names] = [row[-3]]
    return dict_for_deltas


def get_min_deltas(infile, way_to_tackle="microbial"):
    """load a tsv file including all delta_c for all branches. It will return the min in an algn"""
    dict_for_deltas = get_deltas(infile)
    minimal_delta = {}
    if way_to_tackle == "microbial":
        for names in tqdm(dict_for_deltas.keys()):
            uni_id = ast.literal_eval(names.replace(" ", ","))
            for ids in ALGN_NAMES:
                if all(name in ids for name in uni_id):
                    minimal_delta[ids] = numpy.min(
                        [float(i) for i in dict_for_deltas[names]]
                    )
    else:
        minimal_delta = {
            k: min(float(i) for i in v) for k, v in dict_for_deltas.items()
        }
    return minimal_delta


def get_ens(infile, based_on="min_delta", way_to_tackle="microbial"):
    """load a tsv file including ens, tau, algn_name, delta all branches etc... It will return all ens"""
    rows = load_delimited(infile, header=False, sep="\t",)[
        1
    ][1:]

    dict_for_ens = {}
    for row in rows:
        names = row[0]
        if names in dict_for_ens:
            dict_for_ens[names][0].append(float(row[3]))
            dict_for_ens[names][1].append(float(row[-3]))
        else:
            dict_for_ens[names] = [[float(row[3])], [float(row[-3])]]
    if based_on == "min_delta":
        result = {k: v[0][v[1].index(min(v[1]))] for k, v in dict_for_ens.items()}
    else:
        result = {k: v[0] for k, v in dict_for_ens.items()}

    output = {}
    if way_to_tackle == "microbial":
        for names in tqdm(result.keys()):
            uni_id = ast.literal_eval(names.replace(" ", ","))
            for ids in ALGN_NAMES:
                if all(name in ids for name in uni_id):
                    output[ids] = result[names]
    else:
        output = result
    return output


def get_seqlen(infile=None, based_on="degapped"):
    """get length of alignments for microB, by default"""
    dg = de_gap()
    loader = get_app("load_aligned", format="nexus", moltype="dna")
    return (
        {memb.unique_id: dg(loader(memb)).seq_len for memb in tqdm(DS.members)}
        if based_on == "degapped"
        else {memb.unique_id: loader(memb).seq_len for memb in tqdm(DS.members)}
    )


def get_mu(model_name: str, based_on="min_delta"):
    """get total mutation rate branch-wise for microB, by default"""
    loader = get_app("load_db")
    ds = open_data_store(
        f"/Users/yapenglang/repos/YapengAnalysis/data/entropy_on_microb_21Aug/model_results_microbial_degap_28Aug/model_result_{model_name}.sqlitedb",
        mode="r",
    )

    # infile means, if based on min_delta, first get min_delta, then get relevant mu for branch
    model_name = model_name.lower()
    rows = load_delimited(
        f"/Users/yapenglang/repos/YapengAnalysis/data/entropy_on_microb_21Aug/model_stat_microbial_degap_28Aug/ndlc_and_delta/{model_name}.tsv",
        header=False,
        sep="\t",
    )[1][1:]

    dict_for_ = {}
    for row in rows:
        names = row[0]
        if names in dict_for_:
            dict_for_[names][0].append(row[1])
            dict_for_[names][1].append(float(row[-3]))
        else:
            dict_for_[names] = [[row[1]], [float(row[-3])]]

    if based_on == "min_delta":
        result = {k: v[0][v[1].index(min(v[1]))] for k, v in dict_for_.items()}
    else:
        raise ValueError("need more func in this block!")

    output = {}
    for names in tqdm(result.keys()):
        uni_id = ast.literal_eval(names.replace(" ", ","))
        for memb in ds.completed:
            if all(name in memb.unique_id for name in uni_id):
                branch_name = result[names]
                model = loader(memb)
                try:
                    q = model.lf.get_all_rate_matrices()[(branch_name,)].to_array()
                except KeyError:
                    print(f"no rate m in fit{uni_id}")
                mu = -q[numpy.diag_indices(4)].sum()
                output[memb.unique_id] = mu
    return output


def to_bounds(ds_infile, ens_infile, bound = "upper", param = "rate"):
    """give a ds stored bounds prox, with specified a param, return y = ens(tau_), x = proximity to bounds"""
    ds = open_data_store(ds_infile, mode="r")
    loader = get_app("load_db")

    proximity_dict={}
    for memb in ds.completed:
        rec = loader(memb)
        uid = list(rec.keys())[0]
        edges = list(rec[uid]['edge params'].keys())

        proximity_dict[uid] = {}
        if param in ["rate", "length"]:
            for e in edges:
                edge_params = rec[uid]['edge params'][e]
                all_params = [edge_params.pop("length")[bound]]
                if param == "rate":
                    all_params = [row[bound] for row in list(edge_params.values())]

                proximity_dict[uid][e] = min(all_params)
        else:
            proximity_dict[uid] = min(
                row[bound] for row in list(rec[uid]["motif params"].values())
            )
    
    transits  = load_delimited(ens_infile, header=False, sep="\t",)[
        1
    ][1:]
    x=[]
    y=[]
    for r in transits:
        y.append(float(r[2]))
        if param in ["rate", "length"]:
            x.append(proximity_dict[r[0]][r[1]])
        else:
            x.append(proximity_dict[r[0]])
    return x, y
    
    