import pickle

from cogent3 import load_tree, make_tree
from cogent3.util.deserialise import deserialise_object
from numpy import array, eye

from phylo_limits import check_ident


a_identity = eye(4, dtype=float)
a_chainsaw = array(
    [
        [0.9001, 0.1, 0.1, 0.1],
        [0.9000001, 0.1, 0.9, 0.1],
        [0.1, 0.9, 0.1, 0.1],
        [0.1, 0.1, 0.1, 0.9],
    ]
)
a_sympathetic = array(
    [
        [0.9000000001, 0.1, 0.1, 0.1],
        [0.9000001, 0.1, 0.9, 0.1],
        [0.1, 0.9, 0.1, 0.1],
        [0.1, 0.1, 0.1, 0.9],
    ]
)
a_sympathetic2 = array(
    [
        [0.9001, 0.1, 0.1, 0.1],
        [0.1, 0.9, 0.1, 0.1],
        [0.1, 0.1, 0.9, 0.899999],
        [0.1, 0.1, 0.1, 0.9],
    ]
)
a_limit = array(
    [
        [
            0.138523008245210882,
            0.361476991754789756,
            0.138523008245208579,
            0.361476991754788590,
        ],
        [
            0.138523008245210077,
            0.361476991754790200,
            0.138523008245209300,
            0.361476991754789201,
        ],
        [
            0.138523008245208579,
            0.361476991754788590,
            0.138523008245210827,
            0.361476991754789589,
        ],
        [
            0.138523008245209328,
            0.361476991754789201,
            0.138523008245209994,
            0.361476991754790034,
        ],
    ]
)  # dlc-like
a_limit2 = array(
    [
        [
            0.086643076765260257,
            0.413356923234742824,
            0.086643076765260188,
            0.413356923234742490,
        ],
        [
            0.086643076765228907,
            0.413356923234774132,
            0.086643076765228991,
            0.413356923234774076,
        ],
        [
            0.086643076765260160,
            0.413356923234742546,
            0.086643076765260257,
            0.413356923234742601,
        ],
        [
            0.086643076765228977,
            0.413356923234774243,
            0.086643076765228935,
            0.413356923234773965,
        ],
    ]
)  # chainsaw-like



def test_cherry_picker_1():
    tree = load_tree("data/ident_check/test_tree.newick")
    tree.reassign_names({"edge.1": "inter.2", "edge.0": "inter.3"})

    psubs_dict = {
        tuple(sorted(("a", "inter.1"))): {"class": "Sympathetic"},
        tuple(sorted(("b", "inter.1"))): {"class": "Sympathetic"},
        tuple(sorted(("c", "inter.2"))): {"class": "DLC"},
        tuple(sorted(("d", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("e", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("inter.2", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.3", "inter.2"))): {"class": "Sympathetic"},
    }

    check_ident.LINKED_NODES = []

    check_ident.cherry_picker(tree, psubs_dict)

    assert check_ident.LINKED_NODES == ["inter.2", "inter.1"]


def test_cherry_picker_2():
    tree = load_tree("data/ident_check/test_tree.newick")
    tree.reassign_names({"edge.1": "inter.2", "edge.0": "inter.3"})

    psubs_dict = {
        tuple(sorted(("a", "inter.1"))): {"class": "Sympathetic"},
        tuple(sorted(("b", "inter.1"))): {"class": "Sympathetic"},
        tuple(sorted(("c", "inter.2"))): {"class": "Sympathetic"},
        tuple(sorted(("d", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("e", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("inter.2", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.3", "inter.2"))): {"class": "DLC"},
    }
    
    check_ident.LINKED_NODES = []

    check_ident.cherry_picker(tree,psubs_dict)

    assert check_ident.LINKED_NODES == ["inter.3", "inter.2", "inter.1"]


def test_cherry_picker_3():
    tree = load_tree("data/ident_check/test_tree.newick")
    tree.reassign_names({"edge.1": "inter.2", "edge.0": "inter.3"})

    psubs_dict = {
        tuple(sorted(("a", "inter.1"))): {"class": "Sympathetic"},
        tuple(sorted(("b", "inter.1"))): {"class": "Sympathetic"},
        tuple(sorted(("c", "inter.2"))): {"class": "Sympathetic"},
        tuple(sorted(("d", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("e", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("inter.2", "inter.1"))): {"class": "Sympathetic"},
        tuple(sorted(("inter.3", "inter.2"))): {"class": "DLC"},
    }
    
    check_ident.LINKED_NODES = []

    check_ident.cherry_picker(tree,psubs_dict)

    assert not check_ident.LINKED_NODES


def test_cherry_picker_4():
    tree = load_tree("data/ident_check/test_tree2.newick")
    tree.reassign_names(
        {
            "edge.2": "inter.2",
            "edge.3": "inter.3",
            "edge.5": "inter.4",
            "edge.0": "inter.5",
            "edge.1": "inter.6",
            "edge.4": "inter.7",
        }
    )

    psubs_dict = {
        tuple(sorted(("a", "inter.5"))): {"class": "Sympathetic"},
        tuple(sorted(("b", "inter.5"))): {"class": "Sympathetic"},
        tuple(sorted(("c", "inter.6"))): {"class": "DLC"},
        tuple(sorted(("d", "inter.6"))): {"class": "Sympathetic"},
        tuple(sorted(("e", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("f", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("g", "inter.7"))): {"class": "Sympathetic"},
        tuple(sorted(("h", "inter.7"))): {"class": "DLC"},
        tuple(sorted(("i", "inter.4"))): {"class": "DLC"},
        tuple(sorted(("inter.2", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.3", "inter.1"))): {"class": "Sympathetic"},
        tuple(sorted(("inter.4", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.5", "inter.2"))): {"class": "DLC"},
        tuple(sorted(("inter.6", "inter.2"))): {"class": "DLC"},
        tuple(sorted(("inter.7", "inter.4"))): {"class": "DLC"},
    }
    
    check_ident.LINKED_NODES = []

    check_ident.cherry_picker(tree,psubs_dict)

    assert check_ident.LINKED_NODES == ["inter.6", "inter.2", "inter.1"]


def test_cherry_picker_5():
    tree = load_tree("data/ident_check/test_tree2.newick")
    tree.reassign_names(
        {
            "edge.2": "inter.2",
            "edge.3": "inter.3",
            "edge.5": "inter.4",
            "edge.0": "inter.5",
            "edge.1": "inter.6",
            "edge.4": "inter.7",
        }
    )

    psubs_dict = {
        tuple(sorted(("a", "inter.5"))): {"class": "Sympathetic"},
        tuple(sorted(("b", "inter.5"))): {"class": "Sympathetic"},
        tuple(sorted(("c", "inter.6"))): {"class": "Sympathetic"},
        tuple(sorted(("d", "inter.6"))): {"class": "Sympathetic"},
        tuple(sorted(("e", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("f", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("g", "inter.7"))): {"class": "Sympathetic"},
        tuple(sorted(("h", "inter.7"))): {"class": "DLC"},
        tuple(sorted(("i", "inter.4"))): {"class": "Sympathetic"},
        tuple(sorted(("inter.2", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.3", "inter.1"))): {"class": "Sympathetic"},
        tuple(sorted(("inter.4", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.5", "inter.2"))): {"class": "DLC"},
        tuple(sorted(("inter.6", "inter.2"))): {"class": "DLC"},
        tuple(sorted(("inter.7", "inter.4"))): {"class": "DLC"},
    }
    
    check_ident.LINKED_NODES = []

    check_ident.cherry_picker(tree,psubs_dict)

    assert check_ident.LINKED_NODES == ["inter.7", "inter.4", "inter.1"]


def test_cherry_picker_6():
    tree = load_tree("data/ident_check/test_tree2.newick")
    tree.reassign_names(
        {
            "edge.2": "inter.2",
            "edge.3": "inter.3",
            "edge.5": "inter.4",
            "edge.0": "inter.5",
            "edge.1": "inter.6",
            "edge.4": "inter.7",
        }
    )

    psubs_dict = {
        tuple(sorted(("a", "inter.5"))): {"class": "Sympathetic"},
        tuple(sorted(("b", "inter.5"))): {"class": "Sympathetic"},
        tuple(sorted(("c", "inter.6"))): {"class": "Sympathetic"},
        tuple(sorted(("d", "inter.6"))): {"class": "Sympathetic"},
        tuple(sorted(("e", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("f", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("g", "inter.7"))): {"class": "Sympathetic"},
        tuple(sorted(("h", "inter.7"))): {"class": "DLC"},
        tuple(sorted(("i", "inter.4"))): {"class": "Sympathetic"},
        tuple(sorted(("inter.2", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.3", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.4", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.5", "inter.2"))): {"class": "DLC"},
        tuple(sorted(("inter.6", "inter.2"))): {"class": "DLC"},
        tuple(sorted(("inter.7", "inter.4"))): {"class": "Sympathetic"},
    }
    
    check_ident.LINKED_NODES = []

    check_ident.cherry_picker(tree,psubs_dict)

    assert check_ident.LINKED_NODES == ["inter.3", "inter.1"]


def test_cherry_picker_7():
    tree = load_tree("data/ident_check/test_tree2.newick")
    tree.reassign_names(
        {
            "edge.2": "inter.2",
            "edge.3": "inter.3",
            "edge.5": "inter.4",
            "edge.0": "inter.5",
            "edge.1": "inter.6",
            "edge.4": "inter.7",
        }
    )

    psubs_dict = {
        tuple(sorted(("a", "inter.5"))): {"class": "Sympathetic"},
        tuple(sorted(("b", "inter.5"))): {"class": "Sympathetic"},
        tuple(sorted(("c", "inter.6"))): {"class": "Sympathetic"},
        tuple(sorted(("d", "inter.6"))): {"class": "Sympathetic"},
        tuple(sorted(("e", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("f", "inter.3"))): {"class": "DLC"},
        tuple(sorted(("g", "inter.7"))): {"class": "Sympathetic"},
        tuple(sorted(("h", "inter.7"))): {"class": "DLC"},
        tuple(sorted(("i", "inter.4"))): {"class": "Sympathetic"},
        tuple(sorted(("inter.2", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.3", "inter.1"))): {"class": "Sympathetic"},
        tuple(sorted(("inter.4", "inter.1"))): {"class": "DLC"},
        tuple(sorted(("inter.5", "inter.2"))): {"class": "DLC"},
        tuple(sorted(("inter.6", "inter.2"))): {"class": "DLC"},
        tuple(sorted(("inter.7", "inter.4"))): {"class": "Sympathetic"},
    }
    
    check_ident.LINKED_NODES = []

    check_ident.cherry_picker(tree,psubs_dict)

    assert not check_ident.LINKED_NODES


# the below tests would check in case bifurc. nodes've been dropped
def test_reroot1():
    tree = check_ident.reroot("((A,(X,Y)edge2)edge0,(C,D)edge1)root;", "edge2")
    assert len(tree.get_node_names()) - len(tree.get_tip_names()) == 4


def test_reroot2():
    tree = check_ident.reroot("((A,B)edge.0,(C,D)edge.1)root;", "edge.1")
    assert len(tree.get_node_names()) - len(tree.get_tip_names()) == 3


# TODO: modularised these test below, this test should nail down if check_ident could precisely detect the unid node
def test_check_ident_core():
    with open("data/ident_check/unid_psubs.pkl", "rb") as f:
        psubs_dict = pickle.load(f)
    lf = deserialise_object("data/ident_check/unid_lf_case1.json")
    tree_str, DICT_PSUBS, renaming_projection = check_ident.rename(
        tree=lf.tree, psubs_dict=psubs_dict
    )
    tree = make_tree(tree_str)
    n = set(tree.get_node_names()) - set(
        tree.get_tip_names()
    )  # all internal nodes N, a set
    # loop for re-rooting
    bad_nodes = check_ident.check_ident_rerooting(n, tree_str, DICT_PSUBS,which=True)
    print(bad_nodes)
    print(tree_str)
    print(
        DICT_PSUBS[
            tuple(sorted((list(bad_nodes)[0], "ESAG_RS01755-gc0.47%")))
        ]["class"]
    )


def test_check_ident_core2():
    with open("data/ident_check/unid_psubs2.pkl", "rb") as f:
        psubs_dict = pickle.load(f)
    tree_str, DICT_PSUBS, renaming_projection = check_ident.rename(
        tree=make_tree("((A,B)edge.0,(C,D)edge.1);"), psubs_dict=psubs_dict
    )
    tree = make_tree(tree_str)
    n = set(tree.get_node_names()) - set(
        tree.get_tip_names()
    )  # all internal nodes N, a set
    # loop for re-rooting
    bad_nodes = check_ident.check_ident_rerooting(n, tree_str, DICT_PSUBS,which=True)
    print(bad_nodes)
    print(tree_str)
    print(DICT_PSUBS)


def test_check_ident_core3():
    # test the projection dict
    with open("data/ident_check/unid_psubs.pkl", "rb") as f:
        psubs_dict = pickle.load(f)
    lf = deserialise_object("data/ident_check/unid_lf_case1.json")
    tree_str, DICT_PSUBS, renaming_projection = check_ident.rename(
        tree=lf.tree, psubs_dict=psubs_dict
    )
    renaming_projection = {v: k for k, v in renaming_projection.items()}
    tree = make_tree(tree_str)
    n = set(tree.get_node_names()) - set(
        tree.get_tip_names()
    )  # all internal nodes N, a set
    # loop for re-rooting
    bad_nodes = check_ident.check_ident_rerooting(n, tree_str, DICT_PSUBS,which=True)
    for i in bad_nodes:
        print(renaming_projection[i])


def test_check_ident():
    lf = deserialise_object("data/ident_check/unid_lf_case1.json")
    assert not check_ident.has_valid_path(lf, which=False)


def test_check_ident2():
    lf = deserialise_object("data/ident_check/id_lf_case2.json")
    assert check_ident.has_valid_path(lf, which=False)


def test_check_ident3():
    # two I
    model = deserialise_object("data/ident_check/unid_model_case2.json")
    result = check_ident.has_valid_path(model.lf, which=True)
    assert result == {"A", "D"}


def test_check_ident4():
    lf = deserialise_object("data/ident_check/id_lf_case2.json")
    assert check_ident.has_valid_path(lf, which=True) == set()


def test_check_ident5():
    # all DLC
    lf = deserialise_object("data/ident_check/id_all_dlc_case.json")
    assert check_ident.has_valid_path(lf, which=True) == set()


def test_check_ident6():
    # all DLC
    lf = deserialise_object("data/ident_check/id_all_dlc_case.json")
    assert check_ident.has_valid_path(lf, which=False) == True
