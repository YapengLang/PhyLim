import pytest

from cogent3 import make_tree

from phylo_limits.classify_matrix import (
    CHAINSAW,
    DLC,
    IDENTITY,
    LIMIT,
    SYMPATHETIC,
)
from phylo_limits.eval_identifiability import (
    BADMTX,
    BADNODES,
    IdentCheckRes,
    IdentifiabilityCheck,
    ModelMatrixCategories,
    break_path,
    eval_mcats,
    eval_paths,
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


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ({("1",): DLC, ("2",): DLC, ("3",): DLC}, set()),
        ({("1",): DLC, ("2",): DLC, ("3",): SYMPATHETIC}, set()),
        ({("1",): DLC, ("2",): DLC, ("3",): CHAINSAW}, {"3"}),
        ({("1",): DLC, ("2",): DLC, ("3",): IDENTITY}, set()),
        ({("1",): DLC, ("2",): IDENTITY, ("3",): CHAINSAW}, {"3"}),
        ({("1",): DLC, ("2",): DLC, ("3",): LIMIT}, set()),
    ],
)
def test_eval_mcats_not_strict(test_input, expected):
    """we assume all matrices are correctly labelled"""
    assert eval_mcats(test_input, strict=False) == expected


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ({("1",): DLC, ("2",): DLC, ("3",): DLC}, set()),
        ({("1",): DLC, ("2",): DLC, ("3",): SYMPATHETIC}, set()),
        ({("1",): DLC, ("2",): DLC, ("3",): CHAINSAW}, {"3"}),
        ({("1",): DLC, ("2",): DLC, ("3",): IDENTITY}, {"3"}),
        ({("1",): DLC, ("2",): IDENTITY, ("3",): CHAINSAW}, {"2", "3"}),
        ({("1",): DLC, ("2",): DLC, ("3",): LIMIT}, set()),
    ],
)
def test_eval_mcats_strict(test_input, expected):
    """we assume all matrices are correctly labelled"""
    assert eval_mcats(test_input, strict=True) == expected


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ({("1",): DLC, ("2",): DLC, ("3",): DLC}, set()),
        ({("1",): DLC, ("2",): DLC, ("3",): SYMPATHETIC}, set()),
        ({("1",): DLC, ("2",): DLC, ("3",): IDENTITY}, set()),
        ({("1",): DLC, ("2",): SYMPATHETIC, ("3",): SYMPATHETIC}, set()),
        ({("1",): SYMPATHETIC, ("2",): SYMPATHETIC, ("3",): SYMPATHETIC}, {"root"}),
        ({("1",): SYMPATHETIC, ("2",): SYMPATHETIC, ("3",): LIMIT}, {"root"}),
    ],
)
def test_eval_paths_starshape(test_input, expected):
    """we assume `eval_mcats` is good and no bad matrices in the input"""
    from cogent3 import make_tree

    tree = make_tree("(1,2,3);")
    assert eval_paths(test_input, tree) == expected


@pytest.mark.parametrize(
    "test_input,expected",
    [
        (
            {
                ("1",): DLC,
                ("2",): DLC,
                ("3",): DLC,
                ("4",): DLC,
                ("edge.0",): DLC,
                ("edge.1",): SYMPATHETIC,
            },
            set(),
        ),
        (
            {
                ("1",): DLC,
                ("2",): DLC,
                ("3",): DLC,
                ("4",): DLC,
                ("edge.0",): SYMPATHETIC,
                ("edge.1",): SYMPATHETIC,
            },
            {"root"},
        ),
        (
            {
                ("1",): DLC,
                ("2",): DLC,
                ("3",): SYMPATHETIC,
                ("4",): DLC,
                ("edge.0",): DLC,
                ("edge.1",): SYMPATHETIC,
            },
            set(),
        ),
        (
            {
                ("1",): DLC,
                ("2",): DLC,
                ("3",): SYMPATHETIC,
                ("4",): SYMPATHETIC,
                ("edge.0",): DLC,
                ("edge.1",): DLC,
            },
            set(),
        ),
        (
            {
                ("1",): DLC,
                ("2",): DLC,
                ("3",): SYMPATHETIC,
                ("4",): SYMPATHETIC,
                ("edge.0",): DLC,
                ("edge.1",): SYMPATHETIC,
            },
            {"edge.1"},
        ),
    ],
)
def test_eval_paths_smyt_quartet(test_input, expected):
    """we assume `eval_mcats` is good and no bad matrices in the input"""
    from cogent3 import make_tree

    tree = make_tree("((1,2)edge.0,(3,4)edge.1);")
    assert eval_paths(test_input, tree) == expected


@pytest.mark.parametrize(
    "test_input,expected",
    [
        (
            {
                ("1",): DLC,
                ("2",): DLC,
                ("3",): DLC,
                ("4",): DLC,
                ("5",): SYMPATHETIC,
                ("6",): SYMPATHETIC,
                ("7",): DLC,
                ("8",): DLC,
                ("edge.0",): DLC,
                ("edge.1",): SYMPATHETIC,
                ("edge.2",): DLC,
                ("edge.3",): DLC,
                ("edge.4",): SYMPATHETIC,
                ("edge.5",): DLC,
            },
            set(),
        ),
        (
            {
                ("1",): DLC,
                ("2",): DLC,
                ("3",): DLC,
                ("4",): DLC,
                ("5",): SYMPATHETIC,
                ("6",): SYMPATHETIC,
                ("7",): DLC,
                ("8",): DLC,
                ("edge.0",): DLC,
                ("edge.1",): SYMPATHETIC,
                ("edge.2",): DLC,
                ("edge.3",): DLC,
                ("edge.4",): SYMPATHETIC,
                ("edge.5",): SYMPATHETIC,
            },
            {"edge.3", "edge.4"},
        ),
    ],
)
def test_eval_paths_more_taxa(test_input, expected):
    from cogent3 import make_tree

    tree = make_tree(
        "((1,2)edge.0,((3,4)edge.1,((5,6)edge.3,(7,8)edge.5)edge.4)edge.2);"
    )
    assert eval_paths(test_input, tree) == expected


def test_identifiabilitycheck_badmtx():
    psubs = ModelMatrixCategories(
        source="foo", mcats={("1",): DLC, ("2",): DLC, ("3",): CHAINSAW}
    )
    tree = make_tree("(1,2,3);")
    ident_app = IdentifiabilityCheck(strict=False)
    result = ident_app(psubs, tree)
    assert isinstance(result, IdentCheckRes)
    assert result.message_type == BADMTX


def test_identifiabilitycheck_badnodes():
    psubs = ModelMatrixCategories(
        source="foo",
        mcats={("1",): SYMPATHETIC, ("2",): SYMPATHETIC, ("3",): SYMPATHETIC},
    )
    tree = make_tree("(1,2,3);")
    ident_app = IdentifiabilityCheck(strict=False)
    result = ident_app(psubs, tree)
    assert isinstance(result, IdentCheckRes)
    assert result.message_type == BADNODES
