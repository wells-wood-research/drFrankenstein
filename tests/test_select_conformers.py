import sys
from pathlib import Path
from unittest.mock import patch

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
LAB_ROOT = REPO_ROOT / "Laboratory"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(LAB_ROOT))

from OperatingTools import select_conformers


def _config(selection):
    return {
        "conformerGenerationInfo": {"conformerSelection": selection},
        "miscInfo": {"conformerSelectionMethods": "ENERGY", "seed": 1818},
        "runtimeInfo": {
            "madeByConformers": {
                "conformerXyzs": ["c1.xyz", "c2.xyz"],
                "conformerEnergies": {"c1": 0.0, "c2": 1.0},
            }
        },
    }


def test_selection_uses_conformer_generation_info():
    config = _config("DIVERSE")
    with patch.object(select_conformers, "_select_conformers_by_pca_kmeans", return_value=["c2.xyz"]) as diverse:
        result = select_conformers.select_conformer_xyzs(config, 1)

    assert result == ["c2.xyz"]
    diverse.assert_called_once_with(config, 1)


def test_legacy_selection_method_is_rejected():
    with pytest.raises(ValueError, match="Unknown conformer selection method 'BOLTZMANN'"):
        select_conformers.select_conformer_xyzs(_config("BOLTZMANN"), 1)
