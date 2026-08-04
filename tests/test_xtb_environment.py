import os
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
LAB_ROOT = REPO_ROOT / "Laboratory"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(LAB_ROOT) not in sys.path:
    sys.path.insert(0, str(LAB_ROOT))

from OperatingTools import config_handler


def test_set_xtb_environment_uses_otool_xtb_next_to_orca(tmp_path, monkeypatch):
    orca_dir = tmp_path / "orca"
    orca_dir.mkdir()
    orca_exe = orca_dir / "orca"
    xtb_exe = orca_dir / "otool_xtb"
    orca_exe.touch()
    xtb_exe.touch()

    config = {"pathInfo": {"orcaExe": str(orca_exe)}}

    result = config_handler.set_xtb_environment(config)

    assert result == str(xtb_exe)
    assert os.environ["XTBEXE"] == str(xtb_exe)


def test_set_xtb_environment_fails_if_otool_xtb_is_missing(tmp_path, monkeypatch):
    orca_exe = tmp_path / "orca"
    orca_exe.touch()
    config = {"pathInfo": {"orcaExe": str(orca_exe)}}

    monkeypatch.delenv("XTBEXE", raising=False)
    with pytest.raises(FileNotFoundError, match="otool_xtb"):
        config_handler.set_xtb_environment(config)
