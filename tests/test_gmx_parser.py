import pytest
from md_utils.simulation.gmx_parser import MDP, ParseError, parse_ndx

NDX_CONTENT = """\
[ System ]
1 2 3 4 5
[ Protein ]
1 2 3
[ Solvent ]
4 5
"""

MDP_CONTENT = """\
; minimal test MDP
integrator = md
nsteps = 50000
dt = 0.002
"""


# --- parse_ndx ---

def test_parse_ndx_groups(tmp_path):
    ndx = tmp_path / "test.ndx"
    ndx.write_text(NDX_CONTENT)
    groups, _ = parse_ndx(str(ndx))
    assert list(groups.keys()) == ["System", "Protein", "Solvent"]
    assert groups["System"] == [1, 2, 3, 4, 5]
    assert groups["Protein"] == [1, 2, 3]
    assert groups["Solvent"] == [4, 5]


def test_parse_ndx_group_str_contains_names(tmp_path):
    ndx = tmp_path / "test.ndx"
    ndx.write_text(NDX_CONTENT)
    _, group_str = parse_ndx(str(ndx))
    for name in ("System", "Protein", "Solvent"):
        assert name in group_str


def test_parse_ndx_empty_file(tmp_path):
    ndx = tmp_path / "empty.ndx"
    ndx.write_text("")
    groups, group_str = parse_ndx(str(ndx))
    assert groups == {}
    assert group_str == ""


# --- MDP._convert_to_numeric ---

@pytest.fixture
def empty_mdp():
    return MDP()


def test_convert_int(empty_mdp):
    result = empty_mdp._convert_to_numeric("50000")
    assert result == 50000
    assert type(result) is int


def test_convert_float(empty_mdp):
    assert empty_mdp._convert_to_numeric("0.002") == pytest.approx(0.002)


def test_convert_string_passthrough(empty_mdp):
    assert empty_mdp._convert_to_numeric("md") == "md"


def test_convert_list_of_ints(empty_mdp):
    assert empty_mdp._convert_to_numeric("1 2 3") == [1, 2, 3]


def test_convert_empty_string_returns_empty_string(empty_mdp):
    # Regression: previously returned [] which broke skipempty=True
    assert empty_mdp._convert_to_numeric("") == ""
    assert empty_mdp._convert_to_numeric("  ") == ""


def test_convert_non_string_passthrough(empty_mdp):
    assert empty_mdp._convert_to_numeric(42) == 42
    assert empty_mdp._convert_to_numeric([1, 2]) == [1, 2]


# --- MDP.read ---

def test_mdp_read_params(tmp_path):
    f = tmp_path / "test.mdp"
    f.write_text(MDP_CONTENT)
    mdp = MDP(str(f))
    assert mdp["integrator"] == "md"
    assert mdp["nsteps"] == 50000
    assert mdp["dt"] == pytest.approx(0.002)


def test_mdp_read_comment_preserved(tmp_path):
    f = tmp_path / "test.mdp"
    f.write_text(MDP_CONTENT)
    mdp = MDP(str(f))
    comment_values = [v for k, v in mdp.items() if k.startswith("C")]
    assert any("minimal test MDP" in c for c in comment_values)


def test_mdp_read_unknown_line_raises(tmp_path):
    f = tmp_path / "bad.mdp"
    f.write_text("this line has no equals sign\n")
    with pytest.raises(ParseError):
        MDP(str(f))


# --- MDP.write ---

def test_mdp_roundtrip(tmp_path):
    src = tmp_path / "in.mdp"
    out = tmp_path / "out.mdp"
    src.write_text(MDP_CONTENT)
    MDP(str(src)).write(str(out))
    mdp2 = MDP(str(out))
    assert mdp2["integrator"] == "md"
    assert mdp2["nsteps"] == 50000
    assert mdp2["dt"] == pytest.approx(0.002)


def test_mdp_write_skipempty(tmp_path):
    src = tmp_path / "in.mdp"
    out = tmp_path / "out.mdp"
    src.write_text("nsteps = 1000\nrefcoord-scaling = \n")
    MDP(str(src)).write(str(out), skipempty=True)
    content = out.read_text()
    assert "nsteps" in content
    assert "refcoord-scaling" not in content


def test_mdp_write_list_value(tmp_path):
    src = tmp_path / "in.mdp"
    out = tmp_path / "out.mdp"
    src.write_text("energygrps = Protein SOL\n")
    MDP(str(src)).write(str(out))
    assert "Protein SOL" in out.read_text()


def test_mdp_write_default_overwrites_input(tmp_path):
    src = tmp_path / "in.mdp"
    src.write_text("nsteps = 1000\n")
    mdp = MDP(str(src))
    mdp["nsteps"] = 2000
    mdp.write()
    assert MDP(str(src))["nsteps"] == 2000
