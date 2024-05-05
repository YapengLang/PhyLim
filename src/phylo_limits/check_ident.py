import copy
import warnings

from io import StringIO

import numpy

from Bio import Phylo
from cogent3 import make_tree
from cogent3.app.composable import define_app
from cogent3.core.tree import PhyloNode
from cogent3.evolve.parameter_controller import AlignmentLikelihoodFunction
from cogent3.util.dict_array import DictArray

from phylo_limits.matrix_class import classify_psubs


def rename_inter_nodes(tree: PhyloNode, psubs_dict: dict) -> str | dict:
    """give name to every N in a PhyloNode object and psubs dict"""
    n_nodes = set(tree.get_node_names()) - set(tree.get_tip_names())
    result = {}
    k = 2
    for n in n_nodes:
        if n != "root":
            psubs_dict[(f"inter.{k}",)] = psubs_dict.pop((n,))
            result[n] = f"inter.{k}"
            k += 1
    tree.reassign_names(result)
    result["root"] = "inter.1"
    return (
        tree.get_newick(with_node_names=True).replace(";", "inter.1;"),
        psubs_dict,
        result,
    )


def rename_psubs_onpath(tree: PhyloNode, psubs_dict: dict) -> dict:
    """give name to psubs dict. Note: name in psubs dict is uniquely labeled by child and parent.
    name in tree and psubs dict should match"""
    name_rules = {
        (node.name,): tuple(sorted((node.parent.name, node.name)))
        for node in tree.iter_nontips()
    }
    for node in tree.iter_tips():
        name_rules[(node.name,)] = tuple(sorted((node.parent.name, node.name)))
    for key, value in name_rules.items():
        psubs_dict[value] = psubs_dict.pop(key)
    return psubs_dict


def rename(tree: PhyloNode, psubs_dict: dict) -> str | dict:
    """FIXED name in tree and psubs dict, with path indicated by two adjcent nodes. Note: renaming is randomly"""
    tree_copy = copy.deepcopy(tree)
    tree_str, renamed_psubs_dict, renaming_projection = rename_inter_nodes(
        tree_copy, psubs_dict
    )
    renamed_psubs_dict = rename_psubs_onpath(make_tree(tree_str), renamed_psubs_dict)
    return tree_str, renamed_psubs_dict, renaming_projection


def reroot(newick_tree: str, root: str) -> PhyloNode:
    """REROOT a newick tree str with the help of Biopython, keeping the name of each node. TODO: write this function by myself!"""
    handle = StringIO(newick_tree)
    tree = Phylo.read(handle, "newick")
    nontip_num = len(tree.get_nonterminals())
    old_root = tree.root.name
    tree.root_with_outgroup(root)
    tree_str = tree.format("newick").replace(":0.00000", "").replace("\n", "")

    # the below lines tackles with the dropped root in biopy if bifurcating; it should only run once in one lf
    if len(tree.get_nonterminals()) != nontip_num:
        child_on_bifurc = [
            c.get_newick(with_node_names=True)[:-1]
            for c in make_tree(newick_tree).children
        ]
        for achild in child_on_bifurc:
            if achild in tree_str:
                tree_str = tree_str.replace(achild, f"({achild}){old_root}")
                break

    return make_tree(tree_str)


# dict for all psubs with their class ie. an renamed output from check_all_psubs
# DICT_PSUBS = {}
# a dict stores all child-node pair, indicated all children below a node
LINKED_NODES = []


# note: cherry_picker has an alternative, which would track all DLC path from a root lead to tips. benchmark needed to make choice
def cherry_picker(tree: PhyloNode, DICT_PSUBS: dict) -> bool:
    for child in tree.children:
        if DICT_PSUBS[tuple(sorted((child.name, child.parent.name)))]["class"] == "DLC":
            result = True  # default before any recursion
            if child.children:  # check if the next subtree has internal
                result = cherry_picker(
                    child, DICT_PSUBS
                )  # store the result after each recursion
            if result:
                LINKED_NODES.append(tree.name)
                return True  # this break is necessary, it will stop other branch checking after we find a single way to a tip
    return False
    # this break also necessary, as in one stack, if loop over all branch but no dlc found,
    # we return to make sure the result will pass outside the stack to the upper level.


def check_bad_psubs(psubs_dict: dict) -> set:
    bad_nodes = [
        key[0]
        for key, value in psubs_dict.items()
        if value["class"] in ["Identity", "Chainsaw"]
    ]
    return set(bad_nodes)


def check_ident_rerooting(n, tree_str, DICT_PSUBS, which=True) -> bool or set:  # type: ignore
    """`which` if True, return every nodes broken or empty a set, otherwise, return False when a broken node occurs"""
    global LINKED_NODES
    fine_nodes = []
    if which == True:
        bad_nodes = []
        for node in n:
            if node in fine_nodes:  # this node meet cond. so check next one
                continue
            LINKED_NODES = []
            tree = reroot(newick_tree=tree_str, root=node)
            cherry_picker(tree, DICT_PSUBS)
            if not LINKED_NODES:
                bad_nodes += [node]
            else:
                fine_nodes += LINKED_NODES
        return set(bad_nodes)
    else:
        for node in n:
            if set(fine_nodes) == n:  # all nodes have meet condition
                return True
            if node in fine_nodes:  # this node meet cond. so check next one
                continue
            LINKED_NODES = []
            tree = reroot(newick_tree=tree_str, root=node)
            cherry_picker(tree, DICT_PSUBS)
            if not LINKED_NODES:
                return False
            else:
                fine_nodes += LINKED_NODES
        return True


def has_valid_path(
    lf: AlignmentLikelihoodFunction, strictly=True, which=True
) -> bool or set:  # type: ignore
    """for a given tree, whether rooted or not (not rooted at tips), find if: for every internal nodes denotes N,
    there is a path from that node to a tip meets all dlc Matrix on that path. otherwise unidentifiable.
    assumption: if any matrices are identity (if strictly) or chainsaw, unidentifiable. the tree labels are fixed. branches are labeled
    uniquely by adjcent nodes.

    Parameters: `strictly` controls the sensitivity for Identity matrix; if True, func will return False (or set{..}) whenever I occurs, otherwise, func
    will take I as DLC, and report it as warning

            `which` argu if True, func will report set. set{...} indicates which nodes are breaked as per original tree; if
    chainsaws and/or identities occur, return them discard any sympathetic issues. An empty set means ident.
    """

    # first, check if there are any I or Chainsaw
    psubs_dict = classify_psubs(
        lf,
        strictly=strictly,
        label_L=False,
    )
    bad_nodes = set()
    if all(
        value["class"] == "DLC" for value in psubs_dict.values()
    ):  # need this as most natural cases stop here
        return bad_nodes if which == True else True
    if len(bad_nodes := check_bad_psubs(psubs_dict)) > 0:  # stop if any I/C occurs
        return bad_nodes if which == True else False

    # then, do the general check
    # fixed names in tree, psubs_dict
    tree_str, DICT_PSUBS, renaming_projection = rename(
        tree=lf.tree, psubs_dict=psubs_dict
    )
    tree = make_tree(tree_str)
    n = set(tree.get_node_names()) - set(
        tree.get_tip_names()
    )  # all internal nodes N, a set
    bad_nodes = check_ident_rerooting(n, tree_str, DICT_PSUBS, which=which)
    if which == False or len(bad_nodes) < 1:  # if dont want nodes or get an empty set
        return bad_nodes
    renaming_projection = {v: k for k, v in renaming_projection.items()}
    return {renaming_projection[node] for node in bad_nodes}


@define_app
class IdentifiabilityCheck:
    """check the identifiability of a model.
    We make use of the feature of c3 that every tree is root-oriented,
    the question is: can every tip lead to the root with all DLC?"""

    def has_valid_path(self, psubs, tree):
        tree = make_tree(tree)
        tip_to_root_map = tree.child_parent_map()
