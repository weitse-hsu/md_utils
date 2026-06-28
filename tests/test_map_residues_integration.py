"""
Golden-file integration test for map_residues.

Runs the CLI against protein-only copies of the ATP8B/CDC50A input files
(stored in tests/data/) and compares the output TSV byte-for-byte against a
committed golden file.  Because the input files are part of the repository,
the test runs in CI as well as locally.

To regenerate the golden file after an intentional behaviour change, run:

    map_residues \\
        -p tests/data/test_input.pdb \\
        -g tests/data/test_input.gro \\
        -c ATP8B:PROA,PROC,PROD \\
        -c CDC50A:PROB \\
        -o tests/data/map_residues_golden.tsv \\
        -l /dev/null
"""
import subprocess
from pathlib import Path


_DATA   = Path(__file__).parent / "data"
_PDB    = _DATA / "test_input.pdb"
_GRO    = _DATA / "test_input.gro"
_GOLDEN = _DATA / "map_residues_golden.tsv"

_COMPONENTS = ["-c", "ATP8B:PROA,PROC,PROD", "-c", "CDC50A:PROB"]


def _run(tmp_path):
    out_tsv = tmp_path / "output.tsv"
    result = subprocess.run(
        [
            "map_residues",
            "-p", str(_PDB),
            "-g", str(_GRO),
            *_COMPONENTS,
            "-o", str(out_tsv),
            "-l", str(tmp_path / "run.log"),
        ],
        capture_output=True,
        text=True,
    )
    return result, out_tsv


def test_map_residues_output_matches_golden(tmp_path):
    """Full end-to-end run produces a TSV identical to the committed golden file."""
    result, out_tsv = _run(tmp_path)
    assert result.returncode == 0, f"map_residues failed:\n{result.stderr}"
    assert out_tsv.read_text() == _GOLDEN.read_text(), (
        "Output TSV differs from the golden file. "
        "If this change is intentional, regenerate it — see the docstring "
        "at the top of this module."
    )


def test_map_residues_row_count(tmp_path):
    """1498 mapped residue pairs: 1171 ATP8B + 327 CDC50A."""
    _, out_tsv = _run(tmp_path)
    lines = out_tsv.read_text().splitlines()
    assert lines[0].startswith("component\t"), "Missing TSV header"
    assert len(lines) == 1499  # 1 header + 1498 data rows
