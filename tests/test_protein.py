import pytest
from md_utils.structure.protein import (
    AlignmentResult,
    GROResidue,
    PDBResidue,
    _parse_int_field,
    build_mapping_rows,
    check_monotonicity,
    check_nonoverlapping_mappings,
    convert_res_code,
    fitting_needleman_wunsch,
    read_gro_protein_residues,
    read_pdb_components,
    summarize_alignment,
)


# ---------------------------------------------------------------------------
# convert_res_code
# ---------------------------------------------------------------------------

def test_convert_3to1():
    assert convert_res_code("ALA") == "A"
    assert convert_res_code("GLU") == "E"
    assert convert_res_code("trp") == "W"   # case-insensitive


def test_convert_1to3():
    assert convert_res_code("A") == "ALA"
    assert convert_res_code("g") == "GLY"   # case-insensitive


def test_convert_unknown_returns_x():
    assert convert_res_code("XYZ") == "X"
    assert convert_res_code("Z") == "X"


def test_convert_invalid_length_raises():
    with pytest.raises(ValueError):
        convert_res_code("AL")    # 2-letter
    with pytest.raises(ValueError):
        convert_res_code("NALA")  # 4-letter N-terminal variant not accepted


# ---------------------------------------------------------------------------
# _parse_int_field
# ---------------------------------------------------------------------------

def test_parse_int_field():
    assert _parse_int_field("  42  ") == 42
    assert _parse_int_field("") is None
    assert _parse_int_field("   ") is None
    assert _parse_int_field("abc") is None


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _gro_atom(resid, resname, atomname="CA", atomnum=1):
    return f"{resid:5d}{resname:5s}{atomname:>5s}{atomnum:5d}   0.000   0.000   0.000\n"


def _pdb_atom(resname, chain, resid, icode=" ", segid="    "):
    """Builds a 76-char PDB ATOM line with key fields at their correct columns."""
    return (
        "ATOM      1  N   "              # [0:17]
        f"{resname:3s} {chain:1s}{resid:4d}{icode:1s}"  # [17:27]
        "      1.000   2.000   3.000  1.00  0.00      "  # [27:72]
        f"{segid:4s}\n"                 # [72:76]
    )


def _pdb_residue(resid, resname, aa, ordinal, component="A", pdb_id="A"):
    return PDBResidue(
        component=component, pdb_id=pdb_id,
        pdb_resid=resid, icode="", resname=resname, aa=aa,
        pdb_component_ordinal=ordinal,
    )


def _gro_residue(resid, resname, aa, protein_ordinal, file_ordinal=None):
    return GROResidue(
        gro_resid=resid, resname=resname, aa=aa,
        gro_protein_ordinal=protein_ordinal,
        gro_file_ordinal=file_ordinal if file_ordinal is not None else protein_ordinal,
    )


# ---------------------------------------------------------------------------
# read_gro_protein_residues
# ---------------------------------------------------------------------------

def test_read_gro_protein_residues(tmp_path):
    gro = tmp_path / "test.gro"
    gro.write_text(
        "Test system\n"
        "    4\n"
        + _gro_atom(1, "ALA", "N", 1)
        + _gro_atom(1, "ALA", "CA", 2)   # same residue, second atom → deduplicated
        + _gro_atom(2, "GLY", "N", 3)
        + _gro_atom(3, "SOL", "OW", 4)   # water → skipped
        + "   5.000   5.000   5.000\n"
    )
    residues = read_gro_protein_residues(str(gro))
    assert len(residues) == 2
    assert residues[0].resname == "ALA" and residues[0].aa == "A"
    assert residues[0].gro_protein_ordinal == 1
    assert residues[1].resname == "GLY" and residues[1].gro_protein_ordinal == 2


def test_read_gro_invalid_file_raises(tmp_path):
    gro = tmp_path / "bad.gro"
    gro.write_text("only one line\n")
    with pytest.raises(ValueError, match="valid GRO file"):
        read_gro_protein_residues(str(gro))


# ---------------------------------------------------------------------------
# read_pdb_components
# ---------------------------------------------------------------------------

def test_read_pdb_components_by_chain(tmp_path):
    pdb = tmp_path / "test.pdb"
    pdb.write_text(
        _pdb_atom("ALA", "A", 1)
        + _pdb_atom("ALA", "A", 1)    # duplicate atom for same residue → deduplicated
        + _pdb_atom("GLY", "A", 2)
        + _pdb_atom("HOH", "B", 100)  # water → skipped
    )
    components = read_pdb_components(str(pdb), {"chainA": ["A"]}, id_field="chain")
    assert len(components["chainA"]) == 2
    assert components["chainA"][0].resname == "ALA" and components["chainA"][0].pdb_resid == 1
    assert components["chainA"][1].resname == "GLY"


def test_read_pdb_components_by_segid(tmp_path):
    pdb = tmp_path / "test.pdb"
    pdb.write_text(
        _pdb_atom("ALA", "A", 1, segid="PROA")
        + _pdb_atom("GLU", "A", 2, segid="PROA")
    )
    components = read_pdb_components(str(pdb), {"A": ["PROA"]}, id_field="segid")
    assert len(components["A"]) == 2
    assert components["A"][1].aa == "E"


def test_read_pdb_components_invalid_id_field_raises(tmp_path):
    pdb = tmp_path / "test.pdb"
    pdb.write_text("")
    with pytest.raises(ValueError, match="id_field"):
        read_pdb_components(str(pdb), {}, id_field="nonsense")


# ---------------------------------------------------------------------------
# fitting_needleman_wunsch
# ---------------------------------------------------------------------------

def test_nw_exact_match():
    aln = fitting_needleman_wunsch("ACDEF", "ACDEF")
    assert aln.score == 5 * 2   # all 5 matches at match_score=2
    matched = [(q, t) for q, t in aln.pairs if q is not None and t is not None]
    assert len(matched) == 5


def test_nw_finds_subsequence():
    # "AK" should align to positions 1 and 2 inside "MAKE"
    aln = fitting_needleman_wunsch("AK", "MAKE")
    assert aln.gro_start_index == 1
    assert aln.gro_end_index == 2
    matched = [(q, t) for q, t in aln.pairs if q is not None and t is not None]
    assert len(matched) == 2


def test_nw_empty_inputs_raise():
    with pytest.raises(ValueError, match="Empty query"):
        fitting_needleman_wunsch("", "ACDEF")
    with pytest.raises(ValueError, match="Empty target"):
        fitting_needleman_wunsch("ACD", "")


# ---------------------------------------------------------------------------
# build_mapping_rows
# ---------------------------------------------------------------------------

def test_build_mapping_rows():
    pdb_res = [_pdb_residue(1, "ALA", "A", 1), _pdb_residue(2, "GLY", "G", 2)]
    gro_res = [_gro_residue(10, "ALA", "A", 1), _gro_residue(11, "GLY", "G", 2)]
    aln = AlignmentResult(score=4, pairs=[(0, 0), (1, 1), (None, None)],
                          gro_start_index=0, gro_end_index=1)
    rows = build_mapping_rows("A", pdb_res, gro_res, aln)
    assert len(rows) == 2   # (None, None) pair is skipped
    assert rows[0]["pdb_resid"] == 1 and rows[0]["gro_resid"] == 10
    assert rows[0]["aa_match"] is True
    assert rows[1]["pdb_resid"] == 2 and rows[1]["gro_resid"] == 11


# ---------------------------------------------------------------------------
# check_nonoverlapping_mappings
# ---------------------------------------------------------------------------

def _row(gro_ordinal, component="A", pdb_resid=1, gro_resid=10):
    return {"gro_protein_ordinal": gro_ordinal, "component": component,
            "pdb_segid": "A", "pdb_resid": pdb_resid,
            "gro_resid": gro_resid, "gro_resname": "ALA"}


def test_check_nonoverlapping_mappings_clean():
    check_nonoverlapping_mappings([_row(1), _row(2, pdb_resid=2, gro_resid=11)])


def test_check_nonoverlapping_mappings_overlap():
    with pytest.raises(RuntimeError, match="Non-unique GRO mapping"):
        check_nonoverlapping_mappings([_row(1, "A"), _row(1, "B")])


# ---------------------------------------------------------------------------
# check_monotonicity
# ---------------------------------------------------------------------------

def _mono_row(component, ordinal):
    return {"component": component, "gro_protein_ordinal": ordinal}


def test_check_monotonicity_valid(capsys):
    rows = [_mono_row("A", i) for i in (1, 2, 3)]
    check_monotonicity(rows, ["A"])   # should not raise


def test_check_monotonicity_invalid():
    rows = [_mono_row("A", i) for i in (1, 3, 2)]
    with pytest.raises(RuntimeError, match="monotonic"):
        check_monotonicity(rows, ["A"])


# ---------------------------------------------------------------------------
# summarize_alignment
# ---------------------------------------------------------------------------

def test_summarize_alignment_raises_on_low_mapped_fraction(capsys):
    pdb_res = [_pdb_residue(i, "ALA", "A", i) for i in range(1, 11)]  # 10 residues
    gro_res = [_gro_residue(i, "ALA", "A", i) for i in range(1, 11)]
    # 8/10 mapped → fraction=0.8 < default threshold 0.95
    pairs = [(i, i) for i in range(8)] + [(8, None), (9, None)]
    aln = AlignmentResult(score=16, pairs=pairs, gro_start_index=0, gro_end_index=7)
    with pytest.raises(RuntimeError, match="mapped"):
        summarize_alignment("A", pdb_res, gro_res, aln)


def test_summarize_alignment_raises_on_low_identity(capsys):
    pdb_res = [_pdb_residue(i, "ALA", "A", i) for i in range(1, 6)]   # 5 ALA
    gro_res = [_gro_residue(i, "GLY", "G", i) for i in range(1, 6)]   # 5 GLY (all mismatch)
    pairs = [(i, i) for i in range(5)]
    aln = AlignmentResult(score=-5, pairs=pairs, gro_start_index=0, gro_end_index=4)
    with pytest.raises(RuntimeError, match="identity"):
        summarize_alignment("A", pdb_res, gro_res, aln)
