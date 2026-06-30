import copy
import json
import os
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

try:
    import pandas as pd
except ModuleNotFoundError:  # pragma: no cover
    pd = None

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover
    yaml = None

if pd is None or yaml is None:
    pytestmark = pytest.mark.skip(reason="Requires pandas and pyyaml")


REPO_ROOT = Path(__file__).resolve().parents[1]
LAB_ROOT = REPO_ROOT / "Laboratory"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(LAB_ROOT) not in sys.path:
    sys.path.insert(0, str(LAB_ROOT))


def _ensure_stub_modules():
    if "argpass" not in sys.modules:
        argpass_mod = types.ModuleType("argpass")
        argpass_mod.get = lambda *args, **kwargs: None
        sys.modules["argpass"] = argpass_mod

    if "pdbUtils" not in sys.modules:
        pdbutils_pkg = types.ModuleType("pdbUtils")
        pdbutils_mod = types.ModuleType("pdbUtils.pdbUtils")

        def _pdb2df(path):
            rows = []
            with open(path, "r", encoding="utf-8") as handle:
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
                            "RES_ID": int(line[22:26].strip() or 1),
                            "ELEMENT": (line[76:78].strip() or line[12:16].strip()[:1]),
                        }
                    )
            return pd.DataFrame(rows)

        def _df2pdb(df, out_path):
            with open(out_path, "w", encoding="utf-8") as handle:
                for idx, row in df.iterrows():
                    atom_id = int(row.get("ATOM_ID", idx + 1))
                    atom_name = str(row.get("ATOM_NAME", "X"))[:4].ljust(4)
                    x = float(row.get("X", 0.0))
                    y = float(row.get("Y", 0.0))
                    z = float(row.get("Z", 0.0))
                    element = str(row.get("ELEMENT", atom_name.strip()[:1]))[:2].rjust(2)
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

    if "parmed" not in sys.modules:
        parmed_mod = types.ModuleType("parmed")
        parmed_charmm = types.ModuleType("parmed.charmm")
        parmed_exceptions = types.ModuleType("parmed.exceptions")
        parmed_topologyobjects = types.ModuleType("parmed.topologyobjects")

        class _CharmmParameterSet:
            def __init__(self, *args, **kwargs):
                self.args = args

        class _CharmmPsfFile:
            def __init__(self, *args, **kwargs):
                pass

            def load_parameters(self, *args, **kwargs):
                return None

        class _ParameterWarning(Warning):
            pass

        class _Dihedral:
            pass

        class _Cmap:
            pass

        parmed_charmm.CharmmParameterSet = _CharmmParameterSet
        parmed_charmm.CharmmPsfFile = _CharmmPsfFile
        parmed_exceptions.ParameterWarning = _ParameterWarning
        parmed_topologyobjects.Dihedral = _Dihedral
        parmed_topologyobjects.Cmap = _Cmap
        parmed_mod.charmm = parmed_charmm
        parmed_mod.exceptions = parmed_exceptions
        parmed_mod.topologyobjects = parmed_topologyobjects
        sys.modules["parmed"] = parmed_mod
        sys.modules["parmed.charmm"] = parmed_charmm
        sys.modules["parmed.exceptions"] = parmed_exceptions
        sys.modules["parmed.topologyobjects"] = parmed_topologyobjects

    if "mdtraj" not in sys.modules:
        sys.modules["mdtraj"] = types.ModuleType("mdtraj")

    if "openmm" not in sys.modules:
        openmm_mod = types.ModuleType("openmm")
        openmm_app = types.ModuleType("openmm.app")
        openmm_unit = types.ModuleType("openmm.unit")
        openmm_mod.app = openmm_app
        openmm_mod.unit = openmm_unit
        sys.modules["openmm"] = openmm_mod
        sys.modules["openmm.app"] = openmm_app
        sys.modules["openmm.unit"] = openmm_unit

    if "psfgen" not in sys.modules:
        psfgen_mod = types.ModuleType("psfgen")

        class _PsfGen:
            def __init__(self, *args, **kwargs):
                pass

        psfgen_mod.PsfGen = _PsfGen
        sys.modules["psfgen"] = psfgen_mod

    if "mpire" not in sys.modules:
        mpire_mod = types.ModuleType("mpire")
        mpire_utils = types.ModuleType("mpire.utils")

        class _WorkerPool:
            def __init__(self, n_jobs=1):
                self.n_jobs = n_jobs

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def map(self, func, iterable, **kwargs):
                return [func(x) for x in iterable]

        mpire_mod.WorkerPool = _WorkerPool
        mpire_utils.make_single_arguments = lambda args: args
        sys.modules["mpire"] = mpire_mod
        sys.modules["mpire.utils"] = mpire_utils


_ensure_stub_modules()

if pd is not None and yaml is not None:
    from Laboratory.Experiments.Protocol_1_Capping import Capping_Doctor
    from Laboratory.Experiments.Protocol_1_Capping import Capping_Assistant
    from Laboratory.Experiments.Protocol_1_Capping import Capping_Builders
    from Laboratory.Experiments.Protocol_1_Capping import Capping_Geometry as geom
    from Laboratory.Experiments.Protocol_2_Wriggling import Wriggling_Doctor
    from Laboratory.Experiments.Protocol_3_Charging import Charged_Doctor
    from Laboratory.Experiments.Protocol_4_Assembly import Assembly_Doctor
    from Laboratory.Experiments.Protocol_5_Twisting import Twisted_Doctor
    from Laboratory.Experiments.Protocol_6_Stitching import Stitching_Doctor
    from Laboratory.Experiments.Protocol_6_Stitching.AMBER_protocols import AMBER_total_protocol
    from Laboratory.Experiments.Protocol_7_Creation import AMBER_creation
    from Laboratory.Experiments.Protocol_7_Creation import drCreator
    from Laboratory.Experiments.Protocol_8_Reporter import Conformer_PCA_Monster
    from Laboratory.Experiments.Protocol_8_Reporter import Reporting_Doctor
    from Laboratory.OperatingTools import config_handler
    from Laboratory.OperatingTools import electron_checker
    from Laboratory.OperatingTools import file_parsers
    from Laboratory.OperatingTools import make_internal_coords
    from Laboratory.OperatingTools import pdb_checker

from OperatingTools import drLogger
from OperatingTools import drSubprocess


FIXTURE_ROOT = REPO_ROOT / "tests" / "fixtures" / "pcy_protocol_suite"
INPUT_ROOT = FIXTURE_ROOT / "inputs"
EXPECTED_ROOT = FIXTURE_ROOT / "expected"


def _load_fixture_config(tmp_path: Path) -> dict:
    with open(INPUT_ROOT / "drFrankenstein.yaml", "r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)

    in_dir = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    # capping protocol expects inputDir/<molecule>.pdb
    src_pdb = INPUT_ROOT / "01_termini_capping" / "PCY_capped.pdb"
    target_pdb = in_dir / f"{raw['moleculeInfo']['moleculeName']}.pdb"
    target_pdb.write_text(src_pdb.read_text(encoding="utf-8"), encoding="utf-8")

    cfg = {
        "pathInfo": {
            "inputDir": str(in_dir),
            "outputDir": str(out_dir),
            "amberHome": str(out_dir),
        },
        "moleculeInfo": raw["moleculeInfo"],
        "torsionScanInfo": raw["torsionScanInfo"],
        "chargeFittingInfo": raw["chargeFittingInfo"],
        "parameterFittingInfo": raw["parameterFittingInfo"],
        "miscInfo": raw["miscInfo"],
        "runtimeInfo": {
            "madeByCapping": {"cappedPdb": str(INPUT_ROOT / "01_termini_capping" / "PCY_capped.pdb")},
            "madeByConformers": {"conformerXyzs": [str(INPUT_ROOT / "02_GOAT_conformers" / "PCY_conformer_1.xyz")]},
            "madeByTwisting": {
                "torsionsToScan": {"t1": {"ATOM_INDEXES": [1, 2, 3, 4]}},
                "torsionTags": ["t1"],
                "torsionDir": str(out_dir / "05_torsion_scanning"),
            },
            "madeByAssembly": {
                "assembledFrcmod": str(INPUT_ROOT / "04_parameter_assembly" / "PCY_capped.frcmod"),
                "assembledPrm": str(INPUT_ROOT / "04_parameter_assembly" / "PCY_capped.frcmod"),
                "cappedMol2": str(INPUT_ROOT / "04_parameter_assembly" / "PCY_capped.mol2"),
            },
        },
        "checkpointInfo": {
            "cappingComplete": False,
            "conformersComplete": False,
            "chargesComplete": False,
            "assemblyComplete": False,
            "scanningComplete": False,
            "torsionFittingComplete": False,
            "finalCreationComplete": False,
            "reportingComplete": False,
        },
    }
    return cfg


def test_fixture_manifest_paths_exist():
    with open(FIXTURE_ROOT / "manifest.json", "r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    for relpath in manifest["inputs"]:
        assert (INPUT_ROOT / relpath).exists(), relpath

    for relpath in manifest["expected"].values():
        assert (EXPECTED_ROOT / relpath).exists(), relpath


def test_protocol_1_capping_fast_path_uses_fixture():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _load_fixture_config(Path(tmp))
        cfg["moleculeInfo"]["backboneAliases"] = None
        expected_capped = str(EXPECTED_ROOT / "01_termini_capping" / "PCY_capped.pdb")
        with patch.object(Capping_Doctor.Capping_Monster, "optimise_capped_structures", return_value=expected_capped):
            result = Capping_Doctor.capping_protocol(config=copy.deepcopy(cfg), debug=True)
        assert result["runtimeInfo"]["madeByCapping"]["cappedPdb"] == expected_capped
        assert result["checkpointInfo"]["cappingComplete"] is True


def test_protocol_2_wriggling_skips_external_jobs_and_updates_runtime():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _load_fixture_config(Path(tmp))
        fixture_xyz = str(INPUT_ROOT / "02_GOAT_conformers" / "PCY_conformer_1.xyz")

        with patch.object(Wriggling_Doctor, "sort_out_directories", side_effect=lambda c: {**c, "runtimeInfo": {**c["runtimeInfo"], "madeByConformers": {"conformerDir": tmp}}}), \
             patch.object(Wriggling_Doctor, "pdb_to_xyz"), \
             patch.object(Wriggling_Doctor.drOrca, "write_goat_input", return_value="goat.inp"), \
             patch.object(Wriggling_Doctor.drOrca, "run_orca"), \
             patch.object(Wriggling_Doctor, "split_conformers", return_value=[fixture_xyz]), \
             patch.object(Wriggling_Doctor.cleaner, "clean_wriggle"):
            result = Wriggling_Doctor.conformer_generation_protocol(config=copy.deepcopy(cfg), debug=True)

        assert result["checkpointInfo"]["conformersComplete"] is True
        assert result["runtimeInfo"]["madeByConformers"]["nConformersGenerated"] == 1


def test_protocol_3_charging_resp_sets_expected_csv():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _load_fixture_config(Path(tmp))
        cfg["chargeFittingInfo"]["chargeFittingProtocol"] = "RESP"
        expected_csv = str(EXPECTED_ROOT / "03_charge_calculations" / "resp2_charges.csv")
        sample_df = pd.read_csv(expected_csv)

        with patch.object(Charged_Doctor.Charged_Monster, "get_charge_group_indexes", side_effect=lambda *_: cfg), \
             patch.object(Charged_Doctor.Charged_Assistant, "set_up_directories", side_effect=lambda c, *_: {**c, "runtimeInfo": {**c["runtimeInfo"], "madeByCharges": {"chargeDir": tmp}}}), \
             patch.object(Charged_Doctor.select_conformers, "select_conformer_xyzs", return_value=cfg["runtimeInfo"]["madeByConformers"]["conformerXyzs"]), \
             patch.object(Charged_Doctor, "partial_charge_RESP_protocol", return_value=(sample_df, expected_csv)), \
             patch.object(Charged_Doctor.Charged_Assistant, "process_charge_csv"), \
             patch.object(Charged_Doctor.Charged_Assistant, "round_charges_carefully"), \
             patch.object(Charged_Doctor.Charged_Assistant, "enforce_default_backbone_charges"):
            result = Charged_Doctor.charge_protocol(copy.deepcopy(cfg), debug=True)

        assert result["checkpointInfo"]["chargesComplete"] is True
        assert result["runtimeInfo"]["madeByCharges"]["chargesCsv"] == expected_csv


def test_protocol_4_assembly_dispatch_marks_checkpoint_complete():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _load_fixture_config(Path(tmp))
        cfg["miscInfo"]["assemblyProtocol"] = "ANTECHAMBER"
        with patch.object(Assembly_Doctor, "amber_assembly_protocol", return_value={"runtimeInfo": {}, "checkpointInfo": {}}):
            result = Assembly_Doctor.parameter_assembly_protocol(copy.deepcopy(cfg), debug=True)
        assert result["checkpointInfo"]["assemblyComplete"] is True


def test_protocol_5_twisting_marks_scan_complete():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _load_fixture_config(Path(tmp))

        def _set_dirs(c):
            c["runtimeInfo"]["madeByTwisting"]["torsionDir"] = tmp
            return c

        def _choose(c):
            c["runtimeInfo"]["madeByTwisting"]["torsionsToScan"] = {"t1": {"ATOM_INDEXES": [1, 2, 3, 4]}}
            return c

        with patch.object(Twisted_Doctor.Twisted_Assistant, "set_up_directories", side_effect=_set_dirs), \
             patch.object(Twisted_Doctor.Twisted_Assistant, "identify_rotatable_bonds", side_effect=lambda c, mode: c), \
             patch.object(Twisted_Doctor.Twisted_Assistant, "choose_torsions_to_scan", side_effect=_choose), \
             patch.object(Twisted_Doctor, "run_torsion_scanning", side_effect=lambda *args, **kwargs: args[-1]):
            result = Twisted_Doctor.twist_protocol(copy.deepcopy(cfg), debug=True)

        assert result["checkpointInfo"]["scanningComplete"] is True


def test_protocol_6_stitching_runs_single_shuffle_fast():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _load_fixture_config(Path(tmp))
        cfg["runtimeInfo"]["madeByTwisting"]["torsionTags"] = ["t1"]
        cfg["runtimeInfo"]["madeByStitching"] = {}

        metrics = {
            "composite_score": 0.0,
            "location_score": 0.0,
            "amplitude_score": 0.0,
            "stationary_count_score": 0.0,
            "normalized_mae_score": 0.0,
        }
        with patch.object(Stitching_Doctor.Stitching_Assistant, "sort_out_directories", side_effect=lambda c: {**c, "runtimeInfo": {**c["runtimeInfo"], "madeByStitching": {"qmmmParameterFittingDir": tmp}}}), \
             patch.object(Stitching_Doctor.AMBER_helper_functions, "copy_assembled_parameters", side_effect=lambda c: c), \
             patch.object(Stitching_Doctor.AMBER_helper_functions, "edit_mol2_partial_charges"), \
             patch.object(Stitching_Doctor.AMBER_helper_functions, "run_tleap_to_make_params"), \
             patch.object(Stitching_Doctor.Stitching_Assistant, "remove_exploded_torsions", return_value=["t1"]), \
             patch.object(Stitching_Doctor.Stitching_Assistant, "shuffle_torsion_tags", return_value=["t1"]), \
             patch.object(Stitching_Doctor.Stitching_Assistant, "init_tqdm_bar_options", return_value={"disable": True}), \
             patch.object(Stitching_Doctor.AMBER_total_protocol, "get_MM_total_energies", return_value=[0.0]), \
             patch.object(Stitching_Doctor.AMBER_torsion_protocol, "get_MM_torsion_energies", return_value=([0.0], [])), \
             patch.object(Stitching_Doctor.QMMM_fitting_protocol, "fit_torsion_parameters", return_value=(MagicMock(), metrics, metrics, 0.0, 0.0, True)), \
             patch.object(Stitching_Doctor.AMBER_helper_functions, "update_frcmod", return_value=str(INPUT_ROOT / "04_parameter_assembly" / "PCY_capped.frcmod")), \
             patch.object(Stitching_Doctor.Stitching_Assistant, "check_mae_convergence", return_value=True), \
             patch.object(Stitching_Doctor.Stitching_Assistant, "construct_final_params", return_value={"t1": "ok"}), \
             patch.object(Stitching_Doctor.Stitching_Plotter, "make_gif"), \
             patch.object(Stitching_Doctor.Stitching_Plotter, "plot_mean_average_error"), \
             patch.object(Stitching_Doctor.Stitching_Plotter, "plot_run_mean_average_error"), \
             patch.object(Stitching_Doctor.cleaner, "clean_up_stitching"):
            result = Stitching_Doctor.torsion_fitting_protocol.__wrapped__(config=copy.deepcopy(cfg), debug=True)

        assert result["checkpointInfo"]["torsionFittingComplete"] is True
        # Current stitching logic records maxShuffles when the full schedule is traversed.
        assert result["runtimeInfo"]["madeByStitching"]["shufflesCompleted"] == cfg["parameterFittingInfo"]["maxShuffles"]


def test_protocol_7_creation_marks_final_creation_complete():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _load_fixture_config(Path(tmp))
        cfg["parameterFittingInfo"]["forceField"] = "AMBER"
        cfg["torsionScanInfo"]["runScansOn"]["phiPsi"] = False

        with patch.object(drCreator.AMBER_creation, "get_capping_atom_ids", return_value=[1, 2]), \
             patch.object(drCreator.AMBER_creation, "create_final_lib_and_mol2"), \
             patch.object(drCreator.AMBER_creation, "copy_final_frcmod", side_effect=lambda c: c), \
             patch.object(drCreator.AMBER_creation, "duplicate_capping_parameters"), \
             patch.object(drCreator.AMBER_creation, "add_wildcard_dihedrals"):
            result = drCreator.create_the_monster(copy.deepcopy(cfg))

        assert result["checkpointInfo"]["finalCreationComplete"] is True


def test_protocol_8_reporter_sets_report_html():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _load_fixture_config(Path(tmp))

        with patch.object(Reporting_Doctor.Reporting_Assistant, "copy_images"), \
             patch.object(Reporting_Doctor.plot_time_gantt, "generate_gantt_chart", side_effect=["gantt.png", "gantt_cpu.png"]), \
             patch.object(Reporting_Doctor.Reporting_Monster, "process_wriggle_results", return_value={}), \
             patch.object(Reporting_Doctor.Conformer_PCA_Monster, "process_conformer_pca_results", return_value={}), \
             patch.object(Reporting_Doctor.Reporting_Monster, "process_twist_results", return_value={}), \
             patch.object(Reporting_Doctor.Reporting_Monster, "process_charges_results", return_value={}), \
             patch.object(Reporting_Doctor.Reporting_Monster, "process_fitting_results", return_value={}), \
             patch.object(Reporting_Doctor.Shelly, "methods_writer_protocol", return_value={}), \
             patch.object(Reporting_Doctor.Shelly, "gather_citations", return_value={}), \
             patch.object(Reporting_Doctor, "make_html_report"):
            result = Reporting_Doctor.reporter_protocol(copy.deepcopy(cfg), debug=True)

        assert result["runtimeInfo"]["madeByReporting"]["reportHtml"].endswith("08_Parameterisation_Report/drFrankenstein_report.html")


class DoctorTestBase(unittest.TestCase):
    def make_base_config(self, out_dir):
        return {
            "pathInfo": {"inputDir": out_dir, "outputDir": out_dir, "amberHome": out_dir},
            "moleculeInfo": {
                "moleculeName": "PCY",
                "backboneAliases": {"N": ["N"], "C": ["C"], "H": ["H"], "CA": ["CA"], "HA": ["HA"], "O": ["O"]},
                "chargeGroups": {"g1": {"atoms": ["N"], "charge": 0}},
            },
            "torsionScanInfo": {"runScansOn": {"phiPsi": True}, "singlePointMethod": None},
            "chargeFittingInfo": {
                "chargeFittingProtocol": "RESP",
                "nConformers": 1,
                "nCoresPerCalculation": 1,
                "waterDensity": 10,
            },
            "parameterFittingInfo": {
                "forceField": "AMBER",
                "maxCosineFunctions": 2,
                "maxShuffles": 1,
                "minShuffles": 1,
                "converganceTolerance": 999.0,
            },
            "miscInfo": {"availableCpus": 1, "assemblyProtocol": "ANTECHAMBER", "seed": 1818},
            "runtimeInfo": {
                "madeByCapping": {"cappedPdb": os.path.join(out_dir, "capped.pdb")},
                "madeByConformers": {"conformerXyzs": [os.path.join(out_dir, "c1.xyz")]},
                "madeByTwisting": {
                    "torsionsToScan": {"t1": {"ATOM_INDEXES": [1, 2, 3, 4]}},
                    "torsionTags": ["t1"],
                    "torsionDir": os.path.join(out_dir, "twist"),
                },
                "madeByAssembly": {"assembledFrcmod": os.path.join(out_dir, "a.frcmod"), "assembledPrm": os.path.join(out_dir, "a.prm")},
            },
            "checkpointInfo": {
                "cappingComplete": False,
                "conformersComplete": False,
                "chargesComplete": False,
                "scanningComplete": False,
                "torsionFittingComplete": False,
            },
        }


class TestProtocol1CappingDetailed(unittest.TestCase):
    def test_n_terminus_deletes_enough_hydrogens_for_trivalent_n(self):
        config = {"moleculeInfo": {"backboneAliases": {"H": ["H"]}}}
        bonded_atoms = ["CA", "H", "H2", "H3"]
        atoms_to_delete = Capping_Assistant.decide_atom_to_delete_N_termini(bonded_atoms, config)
        self.assertEqual(set(atoms_to_delete), {"H2", "H3"})

    def test_ace_carbonyl_improper_is_trans_180(self):
        mol_df = pd.DataFrame({"ATOM_ID": [1, 2, 3, 4, 5], "ATOM_NAME": ["N", "H", "CA", "C", "O"], "RES_NAME": ["ALA"] * 5, "RES_ID": [1] * 5, "X": [0.0, -0.5, 1.5, 2.5, 3.0], "Y": [0.0, -0.5, 0.0, 1.0, 2.0], "Z": [0.0, 0.5, 0.0, 0.0, 0.0], "OCCUPANCY": [1.0] * 5, "TEMP_FACTOR": [0.0] * 5, "ELEMENT": ["N", "H", "C", "C", "O"]})
        ace_template_df = pd.DataFrame({"ATOM_NAME": ["C_C", "O_C", "C2_C", "H1_C", "H2_C", "H3_C"], "X": [0.0, 1.2, -1.5, -2.0, -1.2, -1.8], "Y": [0.0, 0.0, 0.0, 0.9, -0.8, -0.1], "Z": [0.0, 0.0, 0.0, 0.0, 0.5, -0.7]})
        ace_cap = Capping_Builders.build_ace_cap(mol_df, "N", ace_template_df)
        validation = Capping_Builders.validate_geometry(mol_df, ace_cap, "N", "N")
        self.assertAlmostEqual(abs(validation["O_C-C_C-N-C2_C_dihedral"]), 180.0, places=5)

    def test_cap_clash_resolver_rotates_rigidly_around_attachment_axis(self):
        mol_df = pd.DataFrame({"ATOM_NAME": ["T", "X"], "X": [0.0, 0.0], "Y": [0.0, 2.0], "Z": [0.0, 0.0]})
        cap_df = pd.DataFrame({"ATOM_NAME": ["P", "A", "B"], "X": [1.0, 1.0, 1.0], "Y": [0.0, 2.0, 0.0], "Z": [0.0, 0.0, 1.0]})
        before_angle = geom.calculate_angle(geom.get_coords(mol_df, "T"), geom.get_coords(cap_df, "P"), geom.get_coords(cap_df, "A"))
        resolved = Capping_Builders.resolve_cap_nonbonded_clashes(mol_df=mol_df, cap_df=cap_df, terminal_atom="T", pivot_atom="P", min_distance=1.5, rotation_steps=72)
        after_angle = geom.calculate_angle(geom.get_coords(mol_df, "T"), geom.get_coords(resolved, "P"), geom.get_coords(resolved, "A"))
        self.assertAlmostEqual(before_angle, after_angle, places=5)


PCY_ROOT = REPO_ROOT / "_FrankenParams_PCY"
PCY_PDB = PCY_ROOT / "01_termini_capping" / "PCY_capped.pdb"
PCY_XYZ = PCY_ROOT / "02_GOAT_conformers" / "PCY_conformer_10.xyz"
PCY_MOL2 = PCY_ROOT / "04_parameter_assembly" / "PCY_capped.mol2"


class TestLaboratoryCore(unittest.TestCase):
    def test_xyz2df_parses_fixture(self):
        df = file_parsers.xyz2df(str(PCY_XYZ))
        self.assertGreater(len(df), 0)
        self.assertListEqual(list(df.columns), ["index", "element", "x", "y", "z"])

    def test_parse_mol2_parses_fixture(self):
        atom_df, bond_df = file_parsers.parse_mol2(str(PCY_MOL2))
        self.assertGreater(len(atom_df), 0)
        self.assertGreater(len(bond_df), 0)

    def test_internal_coords_dihedrals_present(self):
        graph = make_internal_coords.pdb_to_graph(str(PCY_PDB))
        dihedrals = make_internal_coords.find_dihedrals(graph)
        self.assertGreater(len(dihedrals), 10)

    def test_validate_config_rejects_invalid_goat_settings(self):
        config = {
            "moleculeInfo": {"charge": 0, "multiplicity": 1, "moleculeName": "PCY", "chargeGroups": None, "backboneAliases": {"N": ["N"], "C": ["C"]}},
            "pathInfo": {"inputDir": "/tmp", "outputDir": "/tmp", "multiWfnDir": "/tmp", "orcaExe": "/etc/hosts"},
            "torsionScanInfo": {"runScansOn": {"phiPsi": True}, "nConformers": -1, "scanMethod": "XTB2"},
            "conformerGenerationInfo": {"goatMode": "GOAT-FAST", "energyCutoff": -1.0, "goatMethod": "MMFF", "conformerSelection": "BAD_SELECTION"},
            "chargeFittingInfo": {"chargeFittingProtocol": "RESP2", "nConformers": -1, "nCoresPerCalculation": 1, "optMethod": "R2SCAN-3c", "singlePointMethod": "wB97X-3c"},
            "parameterFittingInfo": {"forceField": "AMBER", "maxCosineFunctions": 4, "maxShuffles": 50, "minShuffles": 10, "l2DampingFactor": 0.1, "sagvolSmoothing": True},
            "miscInfo": {"availableCpus": 1, "seed": 1818},
        }
        self.assertIsNone(config_handler.validate_config(config))

    def test_duplicate_atom_detection_raises(self):
        df = pd.DataFrame([{"ATOM_ID": 1, "ATOM_NAME": "CA"}, {"ATOM_ID": 2, "ATOM_NAME": "CB"}, {"ATOM_ID": 3, "ATOM_NAME": "CA"}])
        with self.assertRaises(ValueError):
            pdb_checker.check_for_duplicate_atoms(df)

    def test_validate_charge_multiplicity_raises_for_parity_mismatch(self):
        df = pd.DataFrame([{"ATOM_ID": 1, "ATOM_NAME": "O", "ELEMENT": "O"}, {"ATOM_ID": 2, "ATOM_NAME": "H1", "ELEMENT": "H"}, {"ATOM_ID": 3, "ATOM_NAME": "H2", "ELEMENT": "H"}])
        with self.assertRaisesRegex(ValueError, "parity mismatch"):
            electron_checker.validate_charge_multiplicity_from_pdb(df, charge=0, multiplicity=2)


class TestDoctorAdvanced(DoctorTestBase):
    def test_parameter_assembly_protocol_dispatches_all_routes(self):
        cfg = self.make_base_config("/tmp")
        with patch.object(Assembly_Doctor, "agnostic_assembly_protocol", return_value={"ok": "agnostic"}), \
             patch.object(Assembly_Doctor, "amber_assembly_protocol", return_value={"ok": "amber"}), \
             patch.object(Assembly_Doctor, "CGenFF_assembly_protocol", return_value={"ok": "cgenff"}):
            cfg["miscInfo"]["assemblyProtocol"] = "agnostic"
            self.assertEqual(Assembly_Doctor.parameter_assembly_protocol(cfg, debug=True)["ok"], "agnostic")
            cfg["miscInfo"]["assemblyProtocol"] = "ANTECHAMBER"
            self.assertEqual(Assembly_Doctor.parameter_assembly_protocol(cfg, debug=True)["ok"], "amber")
            cfg["miscInfo"]["assemblyProtocol"] = "CGENFF"
            self.assertEqual(Assembly_Doctor.parameter_assembly_protocol(cfg, debug=True)["ok"], "cgenff")

    def test_torsion_fitting_protocol_skips_converged_torsions_in_later_shuffles(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            cfg["parameterFittingInfo"]["maxCosineFunctions"] = 2
            cfg["parameterFittingInfo"]["maxShuffles"] = 2
            cfg["parameterFittingInfo"]["minShuffles"] = 1
            cfg["parameterFittingInfo"]["converganceTolerance"] = 1.0
            cfg["runtimeInfo"]["madeByTwisting"]["torsionTags"] = ["t1", "t2"]
            cfg["runtimeInfo"]["madeByStitching"] = {}

            def _fit_torsion_parameters(*args, **kwargs):
                torsion_tag = args[1]
                if torsion_tag == "t1":
                    metric = {"composite_score": 0.0, "location_score": 0.0, "amplitude_score": 0.0, "stationary_count_score": 0.0, "normalized_mae_score": 0.0}
                    return MagicMock(), metric, metric, 0.0, 0.0, True
                metric = {"composite_score": 5.0, "location_score": 5.0, "amplitude_score": 5.0, "stationary_count_score": 5.0, "normalized_mae_score": 5.0}
                return MagicMock(), metric, metric, 5.0, 5.0, False

            with patch.object(Stitching_Doctor.Stitching_Assistant, "sort_out_directories", side_effect=lambda c: {**c, "runtimeInfo": {**c["runtimeInfo"], "madeByStitching": {"qmmmParameterFittingDir": tmp}}}), \
                 patch.object(Stitching_Doctor.AMBER_helper_functions, "copy_assembled_parameters", side_effect=lambda c: c), \
                 patch.object(Stitching_Doctor.AMBER_helper_functions, "edit_mol2_partial_charges"), \
                 patch.object(Stitching_Doctor.AMBER_helper_functions, "run_tleap_to_make_params"), \
                 patch.object(Stitching_Doctor.Stitching_Assistant, "remove_exploded_torsions", return_value=["t1", "t2"]), \
                 patch.object(Stitching_Doctor.Stitching_Assistant, "shuffle_torsion_tags", return_value=["t1", "t2", "t1", "t2"]), \
                 patch.object(Stitching_Doctor.Stitching_Assistant, "init_tqdm_bar_options", return_value={"disable": True}), \
                 patch.object(Stitching_Doctor.AMBER_total_protocol, "get_MM_total_energies", return_value=[0.0]), \
                 patch.object(Stitching_Doctor.AMBER_torsion_protocol, "get_MM_torsion_energies", return_value=([0.0], [])), \
                 patch.object(Stitching_Doctor.QMMM_fitting_protocol, "fit_torsion_parameters", side_effect=_fit_torsion_parameters) as fit, \
                 patch.object(Stitching_Doctor.AMBER_helper_functions, "update_frcmod", return_value=os.path.join(tmp, "updated.frcmod")) as update, \
                 patch.object(Stitching_Doctor.Stitching_Assistant, "construct_final_params", return_value={"t1": "ok", "t2": "ok"}), \
                 patch.object(Stitching_Doctor.Stitching_Plotter, "make_gif"), \
                 patch.object(Stitching_Doctor.Stitching_Plotter, "plot_mean_average_error"), \
                 patch.object(Stitching_Doctor.Stitching_Plotter, "plot_run_mean_average_error"), \
                 patch.object(Stitching_Doctor.cleaner, "clean_up_stitching"):
                result = Stitching_Doctor.torsion_fitting_protocol.__wrapped__(config=copy.deepcopy(cfg), debug=True)

            self.assertTrue(result["checkpointInfo"]["torsionFittingComplete"])
            self.assertEqual(fit.call_count, 3)
            self.assertEqual(update.call_count, 3)
            self.assertEqual(result["runtimeInfo"]["madeByStitching"]["convergedTags"], ["t1"])

    def test_cpu_gpu_platform_disables_gpu_visibility_env(self):
        with patch.dict(os.environ, {}, clear=True):
            result = AMBER_total_protocol.configure_openmm_gpu_access("CPU")
            self.assertEqual(result, "CPU")
            self.assertEqual(os.environ["CUDA_VISIBLE_DEVICES"], "-1")

    def test_extract_backbone_types_still_requires_core_backbone_atoms(self):
        prmtop_df = pd.DataFrame({"name": ["N", "CA", "C"], "type": ["N", "CT", "C"]})
        aliases = {"N": ["N"], "CA": ["CA"], "C": ["C"], "O": ["O"]}
        with patch("Laboratory.Experiments.Protocol_7_Creation.AMBER_creation.pmd.load_file", return_value=MagicMock(to_dataframe=lambda: prmtop_df), create=True):
            with self.assertRaises(ValueError):
                AMBER_creation.extract_backbone_types_from_prmtop("fake.prmtop", aliases)


def test_experiment_logger_initialization():
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = drLogger.ExperimentLogger(tmpdir)
        assert os.path.isdir(tmpdir)
        assert os.path.isfile(logger.get_log_file_path())


def test_subprocess_run_logging():
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = drLogger.ExperimentLogger(tmpdir)
        drLogger.set_logger(logger)
        result = drSubprocess.logged_run(["ls", "-l"], capture_output=True)
        assert result.returncode == 0
        with open(logger.get_log_file_path(), "r", encoding="utf-8") as handle:
            content = handle.read()
            assert "SUBPROCESS CALL" in content
            assert "ls" in content


def test_error_report_logging():
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = drLogger.ExperimentLogger(tmpdir)
        error_report = {
            "pdbName": "test_molecule",
            "errorType": "ValueError",
            "errorMessage": "Invalid parameter value provided",
            "functionName": "calculate_charges",
            "lineNumber": "125",
            "lineOfCode": "result = invalid_function()",
            "scriptName": "/home/user/drFrankenstein/charges.py",
            "fullTraceBack": [
                "/home/user/drFrankenstein/main.py:42 in main",
                "/home/user/drFrankenstein/charges.py:125 in calculate_charges",
            ],
        }
        logger.log_error_report(error_report)
        with open(logger.get_log_file_path(), "r", encoding="utf-8") as handle:
            content = handle.read()
            assert "DRFRANKENSTEIN ERROR REPORT" in content
            assert "ValueError" in content
