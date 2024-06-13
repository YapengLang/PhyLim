import pytest

from cogent3 import make_tree

from phylo_limits.eval_identifiability import (
    break_path,
    find_bad_nodes,
    find_intersection,
    trav_tip_to_root,
)


@pytest.mark.parametrize(
    "tree_input,expected",
    [
        (
            "((A,B)edge.0,(C,D)edge.1);",
            [
                ["A", "edge.0", "root"],
                ["B", "edge.0"],
                ["C", "edge.1", "root"],
                ["D", "edge.1"],
            ],
        ),
        (
            "(A,B,(C,(D,E)edge.0)edge.1);",
            [
                ["A", "root"],
                ["B", "root"],
                ["C", "edge.1", "root"],
                ["D", "edge.0", "edge.1"],
                ["E", "edge.0"],
            ],
        ),
    ],
)
def test_trav_tip_to_root(tree_input, expected):
    tree = make_tree(tree_input)
    assert trav_tip_to_root(tree) == expected


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ((["A", "edge.0", "root"], {"edge.0"}), [{"A", "edge.0"}]),
        ((["A", "edge.0", "root"], {"edge.0", "A"}), []),
        ((["A", "edge.0", "root"], {"A", "edge.0"}), []),
        ((["A", "edge.0", "root"], {"B", "edge.1"}), [{"A", "edge.0", "root"}]),
        (
            (["A", "edge.0", "edge.1", "edge.2", "root"], {"B", "edge.1", "A"}),
            [{"edge.0", "edge.1"}, {"edge.2", "root"}],
        ),
        (
            (
                ["A", "edge.0", "edge.1", "edge.2", "root"],
                {
                    "B",
                    "edge.1",
                    "A",
                    "edge.0",
                },
            ),
            [{"edge.2", "root"}],
        ),
        (
            (["A", "edge.0", "edge.1", "edge.2", "root"], {"edge.1", "root"}),
            [{"A", "edge.0", "edge.1"}, {"edge.2", "root"}],
        ),
        (
            (
                ["A", "edge.0", "edge.1", "edge.2", "root"],
                {
                    "A",
                    "edge.0",
                },
            ),
            [{"edge.1", "edge.2", "root"}],
        ),
    ],
)
def test_break_path(test_input, expected):
    assert break_path(*test_input) == expected


@pytest.mark.parametrize(
    "sets_input,expected",
    [
        (
            [
                {"1", "1.1", "0"},
                {"2", "1.1"},
                {"3", "1.3"},
                {"4", "1.3"},
                {"1.2", "0"},
                {"1.5", "1.4"},
                {"1.4", "1.6", "1.7"},
            ],
            [
                {"0", "1", "1.1", "1.2", "2"},
                {"1.3", "3", "4"},
                {"1.4", "1.5", "1.6", "1.7"},
            ],
        ),
        (
            [
                {"1", "1.1", "0"},
                {"2", "1.1"},
                {"3", "1.3"},
                {"4", "1.3"},
                {"1.2", "0"},
                {"1.5", "1.4"},
                {"1.4", "1.6", "7"},
                {"1.2", "1.4"},
            ],
            [
                {"2", "1.4", "1.5", "7", "1.1", "1", "1.2", "0", "1.6"},
                {"3", "4", "1.3"},
            ],
        ),
        ([{"1", "2"}, {"3", "4"}], [{"1", "2"}, {"3", "4"}]),
    ],
)
def test_find_intersection(sets_input, expected):
    assert (find_intersection(sets_input)) == expected


@pytest.mark.parametrize(
    "test_input,expected",
    [
        (([{"1", "2", "3"}, {"4", "5"}], {"1", "4"}, {"1", "2", "3", "4", "5"}), set()),
        (([{"1", "2", "3"}, {"4", "5"}], {"1"}, {"1", "2", "3", "4", "5"}), {"4", "5"}),
        (
            ([{"1", "2", "3"}, {"4", "5"}], {"1", "4"}, {"1", "2", "3", "4", "5", "6"}),
            {"6"},
        ),
        (
            ([{"1", "2", "3"}, {"4", "5"}], {"0"}, {"1", "2", "3", "4", "5"}),
            {"1", "2", "3", "4", "5"},
        ),
    ],
)
def test_find_bad_nodes(test_input, expected):
    assert find_bad_nodes(*test_input) == expected


def test_eval_mcats():
    """we assume all matrices are correctly labelled"""
    pass
