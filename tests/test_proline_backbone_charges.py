"""Regression tests for proline-type backbone charge enforcement.

Proline (and 3-/4-hydroxyproline) has a ring nitrogen with no backbone amide H.
`enforce_default_backbone_charges` used to assume a standard backbone (it indexed
`backboneAliases["H"]` and applied AMBER19_DEFAULT_BACKBONE_CHARGES["N"]), so a
correct proline config (which omits the "H" alias) raised KeyError/ValueError.
These tests pin the proline path and guard the standard path against regression.
"""

import os
import sys
import types
import tempfile
import unittest

import numpy as np  # noqa: F401  (imported by module under test)
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LAB_ROOT = os.path.join(REPO_ROOT, "Laboratory")
for _p in (REPO_ROOT, LAB_ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _stub_heavy_imports():
    """Charged_Assistant imports mdtraj / pdbUtils / project helpers at module load.
    enforce_default_backbone_charges needs none of them, so stub them out."""
    sys.modules.setdefault("mdtraj", types.ModuleType("mdtraj"))
    if "pdbUtils" not in sys.modules:
        pkg = types.ModuleType("pdbUtils")
        pkg.pdbUtils = types.ModuleType("pdbUtils.pdbUtils")
        sys.modules["pdbUtils"] = pkg
        sys.modules["pdbUtils.pdbUtils"] = pkg.pdbUtils
    ot = sys.modules.setdefault("OperatingTools", types.ModuleType("OperatingTools"))
    for sub in ("file_parsers", "symmetry_tool"):
        full = f"OperatingTools.{sub}"
        if full not in sys.modules:
            m = types.ModuleType(full)
            sys.modules[full] = m
            setattr(ot, sub, m)


_stub_heavy_imports()

from Laboratory.Experiments.Protocol_3_Charging import Charged_Assistant as CA


def _write_charges_csv(path, rows):
    # default index -> read back via index_col="Unnamed: 0", matching the module
    pd.DataFrame(rows, columns=["ATOM_NAME", "Charge"]).to_csv(path)
    return path


class TestProlineBackboneCharges(unittest.TestCase):
    def _config(self, chargeCsv, backboneAliases):
        return {
            "chargeFittingInfo": {"enforceDefaultBackboneCharges": True},
            "moleculeInfo": {"charge": 0, "backboneAliases": backboneAliases},
            "runtimeInfo": {"madeByCharges": {"chargesCsv": chargeCsv}},
        }

    def test_proline_without_backbone_H_uses_proline_charges(self):
        rows = [
            ["N", -0.30], ["CA", 0.10], ["HA", 0.05], ["C", 0.60], ["O", -0.55],
            ["CB", -0.10], ["HB2", 0.05], ["HB3", 0.05],
            ["CG", -0.10], ["HG2", 0.05], ["HG3", 0.05],
            ["CD", 0.00], ["HD2", 0.05], ["HD3", 0.05],
        ]
        with tempfile.TemporaryDirectory() as d:
            csv = _write_charges_csv(os.path.join(d, "charges.csv"), rows)
            # proline backbone aliases deliberately omit "H"
            aliases = {"N": ["N"], "CA": ["CA"], "HA": ["HA"], "C": ["C"], "O": ["O"]}

            # must not raise (previously KeyError/ValueError on the missing H)
            CA.enforce_default_backbone_charges(self._config(csv, aliases))

            out = pd.read_csv(csv, index_col="Unnamed: 0")
            byName = dict(zip(out["ATOM_NAME"], out["Charge"]))
            self.assertAlmostEqual(byName["N"], -0.2548, places=4)
            self.assertAlmostEqual(byName["CA"], -0.0266, places=4)
            self.assertAlmostEqual(byName["HA"], 0.0641, places=4)
            self.assertAlmostEqual(byName["C"], 0.5896, places=4)
            self.assertAlmostEqual(byName["O"], -0.5748, places=4)
            self.assertNotIn("H", byName)  # no spurious backbone amide H
            self.assertAlmostEqual(out["Charge"].sum(), 0.0, places=3)

    def test_standard_residue_with_H_unchanged(self):
        rows = [
            ["N", -0.30], ["H", 0.20], ["CA", 0.10], ["HA", 0.05], ["C", 0.60], ["O", -0.55],
            ["CB", -0.10], ["HB1", 0.05], ["HB2", 0.05], ["HB3", 0.05],
        ]
        with tempfile.TemporaryDirectory() as d:
            csv = _write_charges_csv(os.path.join(d, "charges.csv"), rows)
            aliases = {"N": ["N"], "H": ["H"], "CA": ["CA"], "HA": ["HA"], "C": ["C"], "O": ["O"]}

            CA.enforce_default_backbone_charges(self._config(csv, aliases))

            out = pd.read_csv(csv, index_col="Unnamed: 0")
            byName = dict(zip(out["ATOM_NAME"], out["Charge"]))
            self.assertAlmostEqual(byName["N"], -0.4157, places=4)
            self.assertAlmostEqual(byName["H"], 0.2719, places=4)
            self.assertAlmostEqual(byName["CA"], 0.0337, places=4)
            self.assertAlmostEqual(byName["C"], 0.5973, places=4)
            self.assertAlmostEqual(out["Charge"].sum(), 0.0, places=3)


if __name__ == "__main__":
    unittest.main()
