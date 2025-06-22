import dataclasses

from enum import Enum
from itertools import chain
from typing import Union

from cogent3.core.tree import PhyloNode

from phylim._version import __version__
from phylim.classify_matrix import (
    CHAINSAW,
    IDENTITY,
    LIMIT,
    SYMPATHETIC,
    MatrixCategory,
    ModelMatrixCategories,
)


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
    for item in path:
        if item in msyms:
            if linked:
                linked = linked | {item}
                split_paths.append(linked)
            linked = set()
        else:
            linked = linked | {item}
    if len(linked) > 1:
        split_paths.append(linked)

    return split_paths


def find_intersection(list_split_paths: list[set]) -> list[set]:
    """Take union of sets based on shared element(s).
    this function's referred to: https://stackoverflow.com/a/6800499."""
    reachable = list(map(set, list_split_paths))
    for i, v in enumerate(reachable):
        for j, k in enumerate(reachable[i + 1 :], i + 1):
            if v & k:
                reachable[i] = v | reachable.pop(j)
                return find_intersection(reachable)
    return reachable


def find_bad_nodes(reachable: list[set], tips: set, nodes: set) -> set:
    """retain nodes in each `reachable` set if there is a node in `tips`"""
    reachable = [x for x in reachable if x & tips]
    good_nodes = set(chain.from_iterable(reachable))
    return nodes - good_nodes


def eval_mcats(mcats: dict[tuple[str, ...], MatrixCategory], strict: bool) -> set:
    """return any chainsaws or identity matrices (depend on `strict`)"""
    bad_categories = {IDENTITY, CHAINSAW} if strict else {CHAINSAW}
    return {k[0] for k, v in mcats.items() if v in bad_categories}


def eval_paths(mcats: dict[tuple[str, ...], MatrixCategory], tree: PhyloNode) -> set:
    """if num of S = 1 or 0, return an empty set; if num of S >= 2, run the path validation algm,
    then return a set for bad nodes."""
    msyms = {k[0] for k, v in mcats.items() if v in {SYMPATHETIC, LIMIT}}
    if len(msyms) < 2:
        return set()
    tips = set(tree.get_tip_names())
    nodes = set(tree.get_node_names()) - tips
    traversed = trav_tip_to_root(tree)

    breaked_paths = []
    for i in traversed:
        breaked_paths = breaked_paths + break_path(i, msyms)

    return find_bad_nodes(find_intersection(breaked_paths), tips, nodes)


class ViolationType(Enum):
    bad_matrices = "bad_matrices"
    bad_nodes = "bad_nodes"
    none = "none"  # the model is identifiable


BADMTX = ViolationType.bad_matrices
BADNODES = ViolationType.bad_nodes
IDENTIFIABLE = ViolationType.none


@dataclasses.dataclass(slots=True)
class IdentCheckRes:
    source: str
    strict: bool
    names: Union[set[str], None]
    violation_type: ViolationType

    def to_rich_dict(self) -> dict:
        result = {
            "source": self.source,
            "strict": self.strict,
            "names": None,
            "violation_type": self.violation_type.name,
            "version": __version__,
        }
        if self.names:
            result["names"] = list(self.names)

        return result

    @property
    def is_identifiable(self) -> bool:
        return self.violation_type is IDENTIFIABLE


def eval_identifiability(
    psubs: ModelMatrixCategories, tree: PhyloNode, strict: bool
) -> IdentCheckRes:
    """check the identifiability of a model fit, provided tree and matrices categories.
    Args:
        strict: controls the sensitivity for Identity matrix (I); if false, treat I as DLC.
    """
    if bad_mtx_names := eval_mcats(psubs.mcats, strict=strict):
        return IdentCheckRes(
            source=psubs.source,
            strict=strict,
            names=bad_mtx_names,
            violation_type=BADMTX,
        )

    bad_node_names = eval_paths(psubs.mcats, tree)
    if bad_node_names:
        return IdentCheckRes(
            source=psubs.source,
            strict=strict,
            names=bad_node_names,
            violation_type=BADNODES,
        )

    return IdentCheckRes(
        source=psubs.source,
        strict=strict,
        names=None,
        violation_type=IDENTIFIABLE,
    )
