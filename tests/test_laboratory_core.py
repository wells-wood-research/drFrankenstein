import copy
import os
import sys
import types
import unittest

import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LAB_ROOT = os.path.join(REPO_ROOT, "Laboratory")
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
if LAB_ROOT not in sys.path:
    sys.path.insert(0, LAB_ROOT)

# Provide a lightweight local fallback for the external pdbUtils dependency.
if "pdbUtils" not in sys.modules:
    pdbutils_pkg = types.ModuleType("pdbUtils")
    pdbutils_mod = types.ModuleType("pdbUtils.pdbUtils")

    def _pdb2df(pdb_path):
        rows = []
        with open(pdb_path, "r") as handle:
            for line in handle:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                rows.append(
                    {
                        "ATOM_ID": int(line[6:11].strip()),
                        "ATOM_NAME": line[12:16].strip(),
                        "X": float(line[30:38].strip()),
                        "Y": float(line[38:46].strip()),
                        "Z": float(line[46:54].strip()),
                        "ELEMENT": line[76:78].strip() or line[12:16].strip()[0],
                    }
                )
        return pd.DataFrame(rows)

    def _df2pdb(df, out_path):
        with open(out_path, "w") as handle:
            for idx, row in df.iterrows():
                atom_id = int(row.get("ATOM_ID", idx + 1))
                atom_name = str(row["ATOM_NAME"])[:4].ljust(4)
                x = float(row["X"])
                y = float(row["Y"])
                z = float(row["Z"])
                element = str(row.get("ELEMENT", atom_name.strip()[0]))[:2].rjust(2)
                handle.write(
                    f"HETATM{atom_id:5d} {atom_name:4s} MOL A   1    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element}\n"
                )
            handle.write("TER\n")

    pdbutils_mod.pdb2df = _pdb2df
    pdbutils_mod.df2pdb = _df2pdb
    pdbutils_mod.pdbUtils = pdbutils_mod
    pdbutils_pkg.pdbUtils = pdbutils_mod
    sys.modules["pdbUtils"] = pdbutils_pkg
    sys.modules["pdbUtils.pdbUtils"] = pdbutils_mod

from Laboratory.OperatingTools import file_parsers
from Laboratory.OperatingTools import electron_checker
from Laboratory.OperatingTools import make_internal_coords
from Laboratory.OperatingTools import pdb_checker
from Laboratory.OperatingTools import set_config_defaults


PCY_ROOT = os.path.join(REPO_ROOT, "__PCY__")
PCY_PDB = os.path.join(PCY_ROOT, "01_termini_capping", "PCY_capped.pdb")
PCY_XYZ = os.path.join(PCY_ROOT, "03_GOAT_conformers", "PCY_conformer_10.xyz")
PCY_MOL2 = os.path.join(PCY_ROOT, "02_parameter_assembly", "PCY.mol2")


class TestFileParsers(unittest.TestCase):
    def test_xyz2df_parses_fixture(self):
        df = file_parsers.xyz2df(PCY_XYZ)
        self.assertEqual(len(df), 48)
        self.assertListEqual(list(df.columns), ["index", "element", "x", "y", "z"])
        self.assertEqual(df.iloc[0]["element"], "N")

    def test_parse_mol2_parses_fixture(self):
        atom_df, bond_df = file_parsers.parse_mol2(PCY_MOL2)
        self.assertEqual(len(atom_df), 48)
        self.assertGreater(len(bond_df), 0)
        self.assertListEqual(
            list(atom_df.columns),
            ["ATOM_ID", "ATOM_NAME", "X", "Y", "Z", "ATOM_TYPE", "RES_ID", "RES_NAME", "CHARGE"],
        )


class TestInternalCoords(unittest.TestCase):
    def setUp(self):
        from pdbUtils import pdbUtils

        self.pdb_df = pdbUtils.pdb2df(PCY_PDB)

    def test_pdb_to_graph_has_expected_nodes_and_edges(self):
        graph = make_internal_coords.pdb_to_graph(PCY_PDB)
        self.assertEqual(graph.number_of_nodes(), 48)
        self.assertGreater(graph.number_of_edges(), 40)
        self.assertIn("N", graph.nodes)
        self.assertIn("CA", graph.nodes)

    def test_find_dihedrals_returns_backbone_pattern(self):
        graph = make_internal_coords.pdb_to_graph(PCY_PDB)
        dihedrals = make_internal_coords.find_dihedrals(graph)
        self.assertGreater(len(dihedrals), 10)
        self.assertTrue(any(d[0] == "N" and d[1] == "CA" and d[2] == "C" for d in dihedrals))

    def test_geometry_calculations_return_reasonable_values(self):
        bond = make_internal_coords.calculate_bond_length(self.pdb_df, ["C", "O"])
        self.assertGreater(bond, 1.0)
        self.assertLess(bond, 1.4)

        angle = make_internal_coords.calculate_angle(self.pdb_df, ["N", "CA", "C"])
        self.assertGreater(angle, 90.0)
        self.assertLess(angle, 140.0)

        dihedral = make_internal_coords.calculate_dihedral(self.pdb_df, ["N", "CA", "C", "O"])
        self.assertGreaterEqual(dihedral, -180.0)
        self.assertLessEqual(dihedral, 180.0)

    def test_geometry_input_length_validation(self):
        with self.assertRaises(ValueError):
            make_internal_coords.calculate_bond_length(self.pdb_df, ["C"])
        with self.assertRaises(ValueError):
            make_internal_coords.calculate_angle(self.pdb_df, ["N", "CA"])
        with self.assertRaises(ValueError):
            make_internal_coords.calculate_dihedral(self.pdb_df, ["N", "CA", "C"])


class TestConfigDefaults(unittest.TestCase):
    def test_apply_defaults_sets_expected_optional_values(self):
        config = {
            "moleculeInfo": {"charge": 0, "multiplicity": 1, "moleculeName": "PCY"},
            "pathInfo": {"multiWfnDir": "/tmp", "orcaExe": "/etc/hosts"},
            "torsionScanInfo": {"scanMethod": "XTB2"},
            "chargeFittingInfo": {
                "chargeFittingProtocol": "RESP2",
                "optMethod": "R2SCAN-3c",
                "singlePointMethod": "wB97X-3c",
            },
            "parameterFittingInfo": {"forceField": "AMBER"},
            "miscInfo": {},
        }

        result = set_config_defaults.apply_defaults_and_validate(copy.deepcopy(config))
        self.assertIn("pathInfo", result)
        self.assertIn("inputDir", result["pathInfo"])
        self.assertTrue(result["pathInfo"]["outputDir"].endswith("PCY_FrankenParams"))
        self.assertIn("runScansOn", result["torsionScanInfo"])
        self.assertIn("phiPsi", result["torsionScanInfo"]["runScansOn"])


class TestPdbChecker(unittest.TestCase):
    def test_duplicate_atom_detection_raises(self):
        df = pd.DataFrame(
            [
                {"ATOM_ID": 1, "ATOM_NAME": "CA"},
                {"ATOM_ID": 2, "ATOM_NAME": "CB"},
                {"ATOM_ID": 3, "ATOM_NAME": "CA"},
            ]
        )
        with self.assertRaises(ValueError):
            pdb_checker.check_for_duplicate_atoms(df)

    def test_unique_atoms_pass(self):
        df = pd.DataFrame(
            [
                {"ATOM_ID": 1, "ATOM_NAME": "CA"},
                {"ATOM_ID": 2, "ATOM_NAME": "CB"},
                {"ATOM_ID": 3, "ATOM_NAME": "CG"},
            ]
        )
        pdb_checker.check_for_duplicate_atoms(df)


class TestElectronChecker(unittest.TestCase):
    def test_calculate_total_electrons_uses_charge_adjustment(self):
        df = pd.DataFrame(
            [
                {"ATOM_ID": 1, "ATOM_NAME": "N", "ELEMENT": "N"},
                {"ATOM_ID": 2, "ATOM_NAME": "H1", "ELEMENT": "H"},
                {"ATOM_ID": 3, "ATOM_NAME": "H2", "ELEMENT": "H"},
                {"ATOM_ID": 4, "ATOM_NAME": "H3", "ELEMENT": "H"},
            ]
        )
        self.assertEqual(electron_checker.calculate_total_electrons(df, charge=1), 9)

    def test_validate_charge_multiplicity_raises_for_both_odd(self):
        df = pd.DataFrame(
            [
                {"ATOM_ID": 1, "ATOM_NAME": "N", "ELEMENT": "N"},
                {"ATOM_ID": 2, "ATOM_NAME": "H1", "ELEMENT": "H"},
                {"ATOM_ID": 3, "ATOM_NAME": "H2", "ELEMENT": "H"},
                {"ATOM_ID": 4, "ATOM_NAME": "H3", "ELEMENT": "H"},
            ]
        )
        with self.assertRaisesRegex(ValueError, "both odd"):
            electron_checker.validate_charge_multiplicity_from_pdb(df, charge=1, multiplicity=1)

    def test_validate_charge_multiplicity_raises_for_parity_mismatch(self):
        df = pd.DataFrame(
            [
                {"ATOM_ID": 1, "ATOM_NAME": "O", "ELEMENT": "O"},
                {"ATOM_ID": 2, "ATOM_NAME": "H1", "ELEMENT": "H"},
                {"ATOM_ID": 3, "ATOM_NAME": "H2", "ELEMENT": "H"},
            ]
        )
        with self.assertRaisesRegex(ValueError, "parity mismatch"):
            electron_checker.validate_charge_multiplicity_from_pdb(df, charge=0, multiplicity=2)


if __name__ == "__main__":
    unittest.main()
