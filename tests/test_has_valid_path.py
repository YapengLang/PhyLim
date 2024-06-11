import pytest

from cogent3 import load_tree, make_tree

from phylo_limits.eval_identifiability import (
    break_path,
    find_intersection,
    trav_tip_to_root,
)


def test_trav_tip_to_root():
    tree = make_tree("((A,B)edge.0,(C,D)edge.1);")
    assert trav_tip_to_root(tree) == [
        ["A", "edge.0", "root"],
        ["B", "edge.0"],
        ["C", "edge.1", "root"],
        ["D", "edge.1"],
    ]


def test_trav_tip_to_root_adtree():
    tree = load_tree("data/ident_check/test_tree.newick")
    print(trav_tip_to_root(tree))


def test_break_path():
    path1, msyms1 = ["A", "edge.0", "root"], {"edge.0"}
    print(break_path(path1, msyms1))


def test_break_path2():
    path2, msyms2 = ["A", "edge.0", "root"], {"edge.0", "A"}
    print(break_path(path2, msyms2))


def test_break_path3():
    path2, msyms2 = ["A", "edge.0", "root"], {"A", "edge.0"}
    print(break_path(path2, msyms2))


def test_break_path4():
    path2, msyms2 = ["A", "edge.0", "root"], {"B", "edge.1"}
    print(break_path(path2, msyms2))


def test_break_path5():
    path2, msyms2 = ["A", "edge.0", "edge.1", "edge.2", "root"], {"B", "edge.1", "A"}
    print(break_path(path2, msyms2))

    # [{"edge.0", "edge.1"},{"edge.2", "root"}]


def test_break_path6():
    path2, msyms2 = ["A", "edge.0", "edge.1", "edge.2", "root"], {
        "B",
        "edge.1",
        "A",
        "edge.0",
    }
    print(break_path(path2, msyms2))

    # [{"edge.2", "root"}]


def test_break_path7():
    path2, msyms2 = ["A", "edge.0", "edge.1", "edge.2", "root"], {"edge.1", "root"}
    print(break_path(path2, msyms2))

    # [{"A", "edge.0", "edge.1"}, {"edge.2", "root"}]


def test_break_path8():
    path2, msyms2 = ["A", "edge.0", "edge.1", "edge.2", "root"], {
        "A",
        "edge.0",
    }
    print(break_path(path2, msyms2))

    # [{ "edge.1", "edge.2", "root"}]


def test_find_intersection():
    g = [
        {"1", "1.1", "0"},
        {"2", "1.1"},
        {"3", "1.3"},
        {"4", "1.3"},
        {"1.2", "0"},
        {"1.5", "1.4"},
        {"1.4", "1.6", "1.7"},
    ]
    print(find_intersection(g))


#     [{'0', '1', '1.1', '1.2', '2'},
#  {'1.3', '3', '4'},
#  {'1.4', '1.5', '1.6', '1.7'}]


def test_find_intersection_1():
    g = [
        {"1", "1.1", "0"},
        {"2", "1.1"},
        {"3", "1.3"},
        {"4", "1.3"},
        {"1.2", "0"},
        {"1.5", "1.4"},
        {"1.4", "1.6", "7"},
        {"1.2", "1.4"},
    ]
    print(find_intersection(g))


# [{'2', '1.4', '1.5', '7', '1.1', '1', '1.2', '0', '1.6'}, {'3', '4', '1.3'}]
