import subprocess
import pytest
from unittest.mock import MagicMock, patch
from md_utils.simulation.gmx_utils import run_gmx_cmd


def _success(stdout="output text"):
    mock = MagicMock()
    mock.returncode = 0
    mock.stdout = stdout
    return mock


def test_run_gmx_cmd_returns_code_and_stdout():
    with patch("subprocess.run", return_value=_success("done")):
        rc, out = run_gmx_cmd(["echo", "hello"], print_output=False)
    assert rc == 0
    assert out == "done"


def test_run_gmx_cmd_failure_raises_runtime_error():
    err = subprocess.CalledProcessError(1, ["gmx", "grompp"], output="fatal error")
    with patch("subprocess.run", side_effect=err):
        with pytest.raises(RuntimeError, match="gmx grompp failed with return code 1"):
            run_gmx_cmd(["gmx", "grompp", "-f", "md.mdp"], print_output=False)


def test_run_gmx_cmd_print_output(capsys):
    with patch("subprocess.run", return_value=_success("hello output")):
        run_gmx_cmd(["echo"], print_output=True)
    assert "hello output" in capsys.readouterr().out


def test_run_gmx_cmd_silent(capsys):
    with patch("subprocess.run", return_value=_success("hidden")):
        run_gmx_cmd(["echo"], print_output=False)
    assert capsys.readouterr().out == ""


def test_run_gmx_cmd_passes_prompt_input():
    with patch("subprocess.run", return_value=_success()) as mock_run:
        run_gmx_cmd(["gmx", "make_ndx"], prompt_input="0\nq\n", print_output=False)
    assert mock_run.call_args[1]["input"] == "0\nq\n"
