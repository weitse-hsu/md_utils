"""
Tests for CLI argument parsing.

All CLIs import general_utils (and some need MDAnalysis / pymol) at module
level, so tests are skipped automatically when those packages are absent.
They run in full locally once dependencies are installed.
"""
import argparse
import pytest

# ---------------------------------------------------------------------------
# Conditional imports — skip entire groups when deps are missing
# ---------------------------------------------------------------------------
map_residues = pytest.importorskip(
    "md_utils.cli.structure.map_residues",
    reason="general_utils not installed",
)
_parse_component = map_residues._parse_component
_parse_landmark  = map_residues._parse_landmark

# These only need general_utils (no GROMACS / PyMOL / MDAnalysis)
create_index_grps  = pytest.importorskip(
    "md_utils.cli.simulation.create_index_grps",  reason="general_utils not installed"
)
process_gmx_traj   = pytest.importorskip(
    "md_utils.cli.simulation.process_gmx_traj",   reason="general_utils not installed"
)
prep_simulation    = pytest.importorskip(
    "md_utils.cli.simulation.prep_simulation",     reason="general_utils not installed"
)

# These additionally need MDAnalysis / pymol — imported lazily inside the tests
# so that the three modules above are still tested even if these are missing.


# ---------------------------------------------------------------------------
# _parse_component
# ---------------------------------------------------------------------------

def test_parse_component_single_id():
    name, ids = _parse_component("ATP8B:PROA")
    assert name == "ATP8B"
    assert ids == ["PROA"]


def test_parse_component_multiple_ids():
    name, ids = _parse_component("ATP8B:PROA,PROC")
    assert name == "ATP8B"
    assert ids == ["PROA", "PROC"]


def test_parse_component_strips_whitespace():
    name, ids = _parse_component(" ATP8B : PROA , PROC ")
    assert name == "ATP8B"
    assert ids == ["PROA", "PROC"]


def test_parse_component_no_colon_raises():
    with pytest.raises(argparse.ArgumentTypeError, match="formatted as"):
        _parse_component("ATP8B-PROA")


def test_parse_component_empty_name_raises():
    with pytest.raises(argparse.ArgumentTypeError, match="name cannot be empty"):
        _parse_component(":PROA")


def test_parse_component_empty_ids_raises():
    with pytest.raises(argparse.ArgumentTypeError, match="No IDs provided"):
        _parse_component("ATP8B:")


# ---------------------------------------------------------------------------
# _parse_landmark
# ---------------------------------------------------------------------------

def test_parse_landmark_single_residue():
    name, tokens = _parse_landmark("ATP8B:100")
    assert name == "ATP8B"
    assert tokens == [100]


def test_parse_landmark_multiple_residues():
    name, tokens = _parse_landmark("ATP8B:100,200,300")
    assert tokens == [100, 200, 300]


def test_parse_landmark_range():
    name, tokens = _parse_landmark("ATP8B:118-132")
    assert tokens == [(118, 132)]


def test_parse_landmark_mixed():
    name, tokens = _parse_landmark("ATP8B:50,118-132,200")
    assert tokens == [50, (118, 132), 200]


def test_parse_landmark_range_start_gt_end_raises():
    with pytest.raises(argparse.ArgumentTypeError, match="start > end"):
        _parse_landmark("ATP8B:132-118")


def test_parse_landmark_invalid_residue_raises():
    with pytest.raises(argparse.ArgumentTypeError, match="Invalid residue"):
        _parse_landmark("ATP8B:abc")


def test_parse_landmark_no_colon_raises():
    with pytest.raises(argparse.ArgumentTypeError, match="formatted as"):
        _parse_landmark("ATP8B-100")


# ---------------------------------------------------------------------------
# map_residues initialize()
# ---------------------------------------------------------------------------

def test_map_residues_initialize_required_args():
    args = map_residues.initialize([
        "-p", "ref.pdb",
        "-g", "sim.gro",
        "-c", "ATP8B:PROA",
    ])
    assert args.pdb == "ref.pdb"
    assert args.gro == "sim.gro"
    assert args.components == [("ATP8B", ["PROA"])]
    assert args.output == "residue_mappings.tsv"   # default
    assert args.pdb_id_field == "segid"            # default


def test_map_residues_initialize_multiple_components():
    args = map_residues.initialize([
        "-p", "ref.pdb", "-g", "sim.gro",
        "-c", "ATP8B:PROA,PROC",
        "-c", "CDC50A:PROB",
    ])
    assert len(args.components) == 2
    assert args.components[1] == ("CDC50A", ["PROB"])


def test_map_residues_initialize_missing_required_raises():
    with pytest.raises(SystemExit):
        map_residues.initialize(["-p", "ref.pdb"])   # missing --gro and --component


def test_map_residues_initialize_landmark_parsed():
    args = map_residues.initialize([
        "-p", "ref.pdb", "-g", "sim.gro",
        "-c", "ATP8B:PROA",
        "--landmark", "ATP8B:118-132,200",
    ])
    assert args.landmarks == [("ATP8B", [(118, 132), 200])]


# ---------------------------------------------------------------------------
# create_index_grps initialize()
# ---------------------------------------------------------------------------

def test_create_index_grps_required_args():
    args = create_index_grps.initialize(["-f", "system.gro"])
    assert args.gro == "system.gro"
    assert args.output == "index.ndx"        # default
    assert args.selections == "selections.txt"  # default


def test_create_index_grps_missing_required_raises():
    with pytest.raises(SystemExit):
        create_index_grps.initialize([])


# ---------------------------------------------------------------------------
# process_gmx_traj initialize()
# ---------------------------------------------------------------------------

def test_process_gmx_traj_required_args():
    args = process_gmx_traj.initialize(["-i", "traj.xtc"])
    assert args.input == "traj.xtc"
    assert args.time_step == 200   # default


def test_process_gmx_traj_missing_required_raises():
    with pytest.raises(SystemExit):
        process_gmx_traj.initialize([])


# ---------------------------------------------------------------------------
# prep_simulation initialize()
# ---------------------------------------------------------------------------

def test_prep_simulation_required_args():
    args = prep_simulation.initialize(["-i", "inputs/"])
    assert args.input_dir == "inputs/"
    assert args.n_replicates == 3    # default


def test_prep_simulation_missing_required_raises():
    with pytest.raises(SystemExit):
        prep_simulation.initialize([])


# ---------------------------------------------------------------------------
# identify_pocket / identify_ss_residues — only imported if deps present
# ---------------------------------------------------------------------------

def test_identify_pocket_initialize():
    pytest.importorskip("MDAnalysis", reason="MDAnalysis not installed")
    identify_pocket = pytest.importorskip("md_utils.cli.structure.identify_pocket")
    args = identify_pocket.initialize(["-i", "complex.pdb"])
    assert args.input == "complex.pdb"
    assert args.resname == "LIG"    # default
    assert args.cutoff == 6.0       # default


def test_identify_ss_residues_initialize():
    pytest.importorskip("pymol", reason="pymol not installed", exc_type=ImportError)
    identify_ss = pytest.importorskip("md_utils.cli.structure.identify_ss_residues", exc_type=ImportError)
    args = identify_ss.initialize(["-i", "protein.pdb", "-s", "H"])
    assert args.input == "protein.pdb"
    assert args.ss_type == "H"
    assert args.min_length == 5    # default
