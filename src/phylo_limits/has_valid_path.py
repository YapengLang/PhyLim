from itertools import chain

from cogent3.app.composable import define_app
from cogent3.core.tree import PhyloNode

from phylo_limits.classify_matrix import SYMPATHETIC, ModelMatrixCategories


def trav_tip_to_root(tree: PhyloNode) -> list[list[str]]:
    """traverse from tip to root, return lists for tip to root path"""
    traversed = []

    for tip in tree.tips():
        ancest_list = [tip.name]
        for i in tip.ancestors():
            ancest_list.append(i.name)
            if any(i.name in path for path in traversed):
                break

        traversed.append(ancest_list)

    return traversed


def break_path(path: list[str], msyms: set) -> list[set]:
    """break the path by msyms(the set of sympathetics)"""
    split_paths = []
    linked = set()
    for i in range(len(path)):
        if path[i] in msyms:
            if linked:
                linked = linked | {path[i]}
                split_paths.append(linked)
            linked = set()
        else:
            linked = linked | {path[i]}
    if len(linked) > 1:
        split_paths.append(linked)

    return split_paths


def find_intersection(list_split_paths: list[set]) -> list[set]:
    """this function's from: https://stackoverflow.com/a/6800499"""
    reachable = list(map(set, list_split_paths))
    for i, v in enumerate(reachable):
        for j, k in enumerate(reachable[i + 1 :], i + 1):
            if v & k:
                reachable[i] = v.union(reachable.pop(j))
                return find_intersection(reachable)
    return reachable


def find_bad_nodes(reachable: list[set], tips: set, nodes: set) -> set:
    reachable = [x for x in reachable if x & tips]
    good_nodes = set(chain.from_iterable(reachable))
    return nodes - good_nodes


@define_app
class IdentifiabilityCheck:
    """check the identifiability of a model.
    We make use of the feature of c3 that every tree is root-oriented,
    the question is: can every tip lead to the root with all DLC?"""

    def main(self, psubs: ModelMatrixCategories, tree: PhyloNode) -> set:
        mcats = psubs.mcats
        msyms = {k for k, v in mcats.items() if v is SYMPATHETIC}
        tips = set(tree.get_tip_names())
        nodes = set(tree.get_node_names()) - tips
        traversed = trav_tip_to_root(tree)

        breaked_paths = []
        for i in traversed:
            breaked_paths = breaked_paths + break_path(i, msyms)

        bad_nodes = find_bad_nodes(find_intersection(breaked_paths), tips, nodes)
        return bad_nodes
