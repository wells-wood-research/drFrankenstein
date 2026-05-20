import copy
import os
import sys
import tempfile
import types
import unittest
from unittest.mock import MagicMock, patch
import pandas as pd
import numpy as np


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LAB_ROOT = os.path.join(REPO_ROOT, "Laboratory")
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
if LAB_ROOT not in sys.path:
    sys.path.insert(0, LAB_ROOT)


# ---- test bootstrap stubs for optional/external deps ----
def _ensure_stub_modules():
    if "argpass" not in sys.modules:
        argpass_mod = types.ModuleType("argpass")

        def _get(*args, **kwargs):
            return None

        argpass_mod.get = _get
        sys.modules["argpass"] = argpass_mod

    if "pdbUtils" not in sys.modules:
        import pandas as _pd

        pdbutils_pkg = types.ModuleType("pdbUtils")
        pdbutils_mod = types.ModuleType("pdbUtils.pdbUtils")

        def _pdb2df(_):
            return _pd.DataFrame(
                [
                    {"ATOM_ID": 1, "ATOM_NAME": "N", "X": 0.0, "Y": 0.0, "Z": 0.0, "ELEMENT": "N", "RES_ID": 1},
                    {"ATOM_ID": 2, "ATOM_NAME": "CA", "X": 1.0, "Y": 0.0, "Z": 0.0, "ELEMENT": "C", "RES_ID": 1},
                    {"ATOM_ID": 3, "ATOM_NAME": "C", "X": 2.0, "Y": 0.0, "Z": 0.0, "ELEMENT": "C", "RES_ID": 1},
                ]
            )

        def _df2pdb(df, out_path):
            with open(out_path, "w") as handle:
                for i, row in df.iterrows():
                    handle.write(
                        f"HETATM{int(row.get('ATOM_ID', i+1)):5d} {str(row['ATOM_NAME'])[:4]:<4} MOL A   1"
                        f"{float(row.get('X', 0.0)):12.3f}{float(row.get('Y', 0.0)):8.3f}{float(row.get('Z', 0.0)):8.3f}"
                        "  1.00  0.00           C\\n"
                    )
                handle.write("TER\\n")

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

            @classmethod
            def from_structure(cls, _):
                return cls()

            def write(self, *args, **kwargs):
                return None

        class _CharmmPsfFile:
            def __init__(self, *args, **kwargs):
                pass

            def load_parameters(self, *args, **kwargs):
                return None

            def save(self, *args, **kwargs):
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

            def read_topology(self, *args, **kwargs):
                return None

            def add_segment(self, *args, **kwargs):
                return None

            def read_coords(self, *args, **kwargs):
                return None

            def write_psf(self, *args, **kwargs):
                return None

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

        def _make_single_arguments(args):
            return args

        mpire_mod.WorkerPool = _WorkerPool
        mpire_utils.make_single_arguments = _make_single_arguments
        sys.modules["mpire"] = mpire_mod
        sys.modules["mpire.utils"] = mpire_utils


_ensure_stub_modules()

from Laboratory.Experiments.Protocol_1_Capping import Capping_Doctor
from Laboratory.Experiments.Protocol_2_Wriggling import Wriggling_Doctor
from Laboratory.Experiments.Protocol_3_Charging import Charged_Doctor
from Laboratory.Experiments.Protocol_4_Assembly import Assembly_Doctor
from Laboratory.Experiments.Protocol_5_Twisting import Twisted_Doctor
from Laboratory.Experiments.Protocol_6_Stitching import Stitching_Doctor
from Laboratory.Experiments.Protocol_6_Stitching.AMBER_protocols import AMBER_total_protocol
from Laboratory.Experiments.Protocol_7_Creation import AMBER_creation
from Laboratory.Experiments.Protocol_8_Reporter import Reporting_Doctor


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


class TestCappingDoctor(DoctorTestBase):
    def test_capping_protocol_ligand_short_circuit(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            cfg["moleculeInfo"]["backboneAliases"] = None
            with patch.object(Capping_Doctor.Capping_Monster, "optimise_capped_structures", return_value="opt.pdb") as opt:
                result = Capping_Doctor.capping_protocol(config=copy.deepcopy(cfg))
            opt.assert_called_once()
            self.assertEqual(result["runtimeInfo"]["madeByCapping"]["cappedPdb"], "opt.pdb")
            self.assertTrue(result["checkpointInfo"]["cappingComplete"])


class TestWrigglingDoctor(DoctorTestBase):
    def test_conformer_generation_protocol_updates_runtime(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            out_xyz = os.path.join(tmp, "PCY_conformer_1.xyz")
            with open(out_xyz, "w") as f:
                f.write("2\n-1.0\nH 0 0 0\nH 0 0 1\n")

            with patch.object(Wriggling_Doctor, "sort_out_directories", side_effect=lambda c: {**c, "runtimeInfo": {**c["runtimeInfo"], "madeByConformers": {"conformerDir": tmp}}}), \
                 patch.object(Wriggling_Doctor, "pdb_to_xyz"), \
                 patch.object(Wriggling_Doctor.drOrca, "write_goat_input", return_value="goat.in"), \
                 patch.object(Wriggling_Doctor.drOrca, "run_orca"), \
                 patch.object(Wriggling_Doctor, "split_conformers", return_value=[out_xyz]), \
                 patch.object(Wriggling_Doctor.cleaner, "clean_wriggle"):
                result = Wriggling_Doctor.conformer_generation_protocol(config=copy.deepcopy(cfg))

            self.assertTrue(result["checkpointInfo"]["conformersComplete"])
            self.assertEqual(result["runtimeInfo"]["madeByConformers"]["nConformersGenerated"], 1)


class TestChargedDoctor(DoctorTestBase):
    def test_charge_protocol_resp_dispatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            cfg["chargeFittingInfo"]["chargeFittingProtocol"] = "RESP"
            with patch.object(Charged_Doctor.Charged_Monster, "get_charge_group_indexes", side_effect=lambda *_: cfg), \
                 patch.object(Charged_Doctor.Charged_Assistant, "set_up_directories", side_effect=lambda c, *_: {**c, "runtimeInfo": {**c["runtimeInfo"], "madeByCharges": {"chargeDir": tmp}}}), \
                 patch.object(Charged_Doctor.Twisted_Assistant, "get_conformer_xyzs", return_value=["c1.xyz"]), \
                 patch.object(Charged_Doctor, "partial_charge_RESP_protocol", return_value=(MagicMock(), os.path.join(tmp, "charges.csv"))), \
                 patch.object(Charged_Doctor.Charged_Assistant, "process_charge_csv"), \
                 patch.object(Charged_Doctor.Charged_Assistant, "round_charges_carefully"):
                result = Charged_Doctor.charge_protocol(copy.deepcopy(cfg), debug=True)

            self.assertTrue(result["checkpointInfo"]["chargesComplete"])
            self.assertTrue(result["runtimeInfo"]["madeByCharges"]["chargesCsv"].endswith("charges.csv"))


class TestAssemblyDoctor(DoctorTestBase):
    def test_parameter_assembly_protocol_dispatches(self):
        cfg = self.make_base_config("/tmp")
        with patch.object(Assembly_Doctor, "agnostic_assembly_protocol", return_value={"ok": "agnostic"}), \
             patch.object(Assembly_Doctor, "amber_assembly_protocol", return_value={"ok": "amber"}), \
             patch.object(Assembly_Doctor, "CGenFF_assembly_protocol", return_value={"ok": "cgenff"}):
            cfg["miscInfo"]["assemblyProtocol"] = "agnostic"
            self.assertEqual(Assembly_Doctor.parameter_assembly_protocol(cfg)["ok"], "agnostic")
            cfg["miscInfo"]["assemblyProtocol"] = "ANTECHAMBER"
            self.assertEqual(Assembly_Doctor.parameter_assembly_protocol(cfg)["ok"], "amber")
            cfg["miscInfo"]["assemblyProtocol"] = "CGENFF"
            self.assertEqual(Assembly_Doctor.parameter_assembly_protocol(cfg)["ok"], "cgenff")


class TestTwistedDoctor(DoctorTestBase):
    def test_twist_protocol_sets_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            cfg["runtimeInfo"]["madeByTwisting"] = {"torsionDir": tmp}

            def _set_dirs(c):
                c["runtimeInfo"]["madeByTwisting"]["torsionDir"] = tmp
                return c

            def _choose_torsions(c):
                c["runtimeInfo"]["madeByTwisting"]["torsionsToScan"] = {"t1": {"ATOM_INDEXES": [1, 2, 3, 4]}}
                return c

            with patch.object(Twisted_Doctor.Twisted_Assistant, "set_up_directories", side_effect=_set_dirs), \
                 patch.object(Twisted_Doctor.Twisted_Assistant, "identify_rotatable_bonds", side_effect=lambda c, mode: c), \
                 patch.object(Twisted_Doctor.Twisted_Assistant, "choose_torsions_to_scan", side_effect=_choose_torsions), \
                 patch.object(Twisted_Doctor, "run_torsion_scanning", side_effect=lambda *args, **kwargs: args[-1]):
                result = Twisted_Doctor.twist_protocol(copy.deepcopy(cfg))
            self.assertTrue(result["checkpointInfo"]["scanningComplete"])


class TestStitchingDoctor(DoctorTestBase):
    def test_profile_fit_score_prefers_matching_curve_shape(self):
        qm = np.array([0.0, 2.0, 0.0, 1.5, 0.0, 1.0, 0.0])
        same = qm.copy()
        shifted = np.array([0.0, 0.0, 2.0, 0.0, 1.5, 0.0, 1.0])
        flat = np.zeros_like(qm)

        same_score = Stitching_Doctor.Stitching_Assistant.calculate_profile_fit_score(qm, same, sampleSpacingDegrees=10)
        shifted_score = Stitching_Doctor.Stitching_Assistant.calculate_profile_fit_score(qm, shifted, sampleSpacingDegrees=10)
        flat_metrics = Stitching_Doctor.Stitching_Assistant.calculate_profile_fit_metrics(qm, flat, sampleSpacingDegrees=10)

        self.assertAlmostEqual(same_score, 0.0)
        self.assertGreater(shifted_score, same_score)
        self.assertEqual(flat_metrics["amplitude_score"], 1.0)
        self.assertEqual(flat_metrics["stationary_count_score"], 1.0)
        self.assertGreaterEqual(flat_metrics["composite_score"], 1.0)

    def test_torsion_fitting_protocol_converges_single_shuffle(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            cfg["runtimeInfo"]["madeByTwisting"]["torsionTags"] = ["t1"]
            cfg["runtimeInfo"]["madeByStitching"] = {}

            with patch.object(Stitching_Doctor.Stitching_Assistant, "sort_out_directories", side_effect=lambda c: {**c, "runtimeInfo": {**c["runtimeInfo"], "madeByStitching": {"qmmmParameterFittingDir": tmp}}}), \
                 patch.object(Stitching_Doctor.AMBER_helper_functions, "copy_assembled_parameters", side_effect=lambda c: c), \
                 patch.object(Stitching_Doctor.AMBER_helper_functions, "edit_mol2_partial_charges"), \
                 patch.object(Stitching_Doctor.AMBER_helper_functions, "run_tleap_to_make_params"), \
                 patch.object(Stitching_Doctor.Stitching_Assistant, "remove_exploded_torsions", return_value=["t1"]), \
                 patch.object(Stitching_Doctor.Stitching_Assistant, "shuffle_torsion_tags", return_value=["t1"]), \
                 patch.object(Stitching_Doctor.Stitching_Assistant, "init_tqdm_bar_options", return_value={"disable": True}), \
                 patch.object(Stitching_Doctor.AMBER_total_protocol, "get_MM_total_energies", return_value=[0.0]), \
                 patch.object(Stitching_Doctor.AMBER_torsion_protocol, "get_MM_torsion_energies", return_value=([0.0], [])), \
                  patch.object(Stitching_Doctor.QMMM_fitting_protocol, "fit_torsion_parameters", return_value=(MagicMock(), {"composite_score": 0.0, "location_score": 0.0, "amplitude_score": 0.0, "stationary_count_score": 0.0, "normalized_mae_score": 0.0}, {"composite_score": 0.0, "location_score": 0.0, "amplitude_score": 0.0, "stationary_count_score": 0.0, "normalized_mae_score": 0.0}, 0.0, 0.0, True)), \
                 patch.object(Stitching_Doctor.AMBER_helper_functions, "update_frcmod", return_value=os.path.join(tmp, "updated.frcmod")), \
                 patch.object(Stitching_Doctor.Stitching_Assistant, "check_mae_convergence", return_value=True), \
                 patch.object(Stitching_Doctor.Stitching_Assistant, "construct_final_params", return_value={"t1": "ok"}), \
                 patch.object(Stitching_Doctor.Stitching_Plotter, "make_gif"), \
                 patch.object(Stitching_Doctor.Stitching_Plotter, "plot_mean_average_error"), \
                 patch.object(Stitching_Doctor.Stitching_Plotter, "plot_run_mean_average_error"), \
                 patch.object(Stitching_Doctor.cleaner, "clean_up_stitching"):
                result = Stitching_Doctor.torsion_fitting_protocol.__wrapped__(config=copy.deepcopy(cfg), debug=True)

            self.assertTrue(result["checkpointInfo"]["torsionFittingComplete"])
            self.assertEqual(result["runtimeInfo"]["madeByStitching"]["shufflesCompleted"], 1)

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

    def test_clean_up_stitching_preserves_final_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            qmmm_dir = os.path.join(tmp, "qmmm")
            molecule_dir = os.path.join(tmp, "molecule")
            torsion_dir = os.path.join(qmmm_dir, "t1")
            os.makedirs(torsion_dir, exist_ok=True)
            os.makedirs(molecule_dir, exist_ok=True)
            cfg["miscInfo"]["cleanUpLevel"] = 3
            cfg["runtimeInfo"]["madeByStitching"] = {
                "qmmmParameterFittingDir": qmmm_dir,
                "moleculeParameterDir": molecule_dir,
                "shufflesCompleted": 5,
                "moleculeFrcmod": os.path.join(molecule_dir, "PCY_2.frcmod"),
            }

            final_png = os.path.join(torsion_dir, "fitting_shuffle_4_nCosines_2.png")
            extra_png = os.path.join(torsion_dir, "fitting_shuffle_1_nCosines_2.png")
            final_frcmod = os.path.join(molecule_dir, "PCY_2.frcmod")
            extra_frcmod = os.path.join(molecule_dir, "PCY_1.frcmod")
            for path in [final_png, extra_png, final_frcmod, extra_frcmod]:
                with open(path, "w") as handle:
                    handle.write("x")

            Stitching_Doctor.cleaner.clean_up_stitching(cfg)

            self.assertTrue(os.path.exists(final_png))
            self.assertTrue(os.path.exists(final_frcmod))
            self.assertFalse(os.path.exists(extra_png))
            self.assertFalse(os.path.exists(extra_frcmod))


class TestAmberTotalProtocol(DoctorTestBase):
    def test_minimisation_uses_current_torsion_indexes(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            scan_dir = os.path.join(tmp, "scan")
            os.makedirs(scan_dir, exist_ok=True)
            for idx in range(2):
                with open(os.path.join(scan_dir, f"orca_scan.{idx:03d}.xyz"), "w") as handle:
                    handle.write("1\ncomment\nH 0.0 0.0 0.0\n")

            torsion_indexes = [4, 5, 6, 7]
            cfg["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]["t1"]["ATOM_INDEXES"] = torsion_indexes

            with patch.object(AMBER_total_protocol.file_parsers, "convert_traj_xyz_to_pdb", return_value=["a.pdb", "b.pdb"]) as convert, \
                 patch.object(AMBER_total_protocol, "run_mm_constrained_em", return_value=[("0", 1.5), ("1", 2.5)]) as minimise, \
                 patch.object(AMBER_total_protocol.Stitching_Assistant, "get_scan_angles_from_orca_inp", return_value=[0, 10]), \
                 patch.object(AMBER_total_protocol.Stitching_Assistant, "rescale_angles_0_360", side_effect=lambda x: x):
                result = AMBER_total_protocol.get_singlepoint_energies_for_torsion_scan(
                    0, scan_dir, tmp, cfg["runtimeInfo"]["madeByCapping"]["cappedPdb"], "fake.prmtop", torsion_indexes, debug=True
                )

            convert.assert_called_once()
            minimise.assert_called_once_with(["a.pdb", "b.pdb"], "fake.prmtop", torsionIndexes=torsion_indexes)
            self.assertEqual(result["Energy"].tolist(), [1.5, 2.5])
            self.assertEqual(result["Angle"].tolist(), [0, 10])


class TestReportingDoctor(DoctorTestBase):
    def test_reporter_protocol_updates_runtime(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            with patch.object(Reporting_Doctor.Reporting_Assistant, "copy_images"), \
                 patch.object(Reporting_Doctor.plot_time_gantt, "generate_gantt_chart", return_value="gantt.png"), \
                 patch.object(Reporting_Doctor.Reporting_Monster, "process_wriggle_results", return_value={}), \
                 patch.object(Reporting_Doctor.Reporting_Monster, "process_twist_results", return_value={}), \
                 patch.object(Reporting_Doctor.Reporting_Monster, "process_charges_results", return_value={}), \
                 patch.object(Reporting_Doctor.Reporting_Monster, "process_fitting_results", return_value={}), \
                 patch.object(Reporting_Doctor.Shelly, "methods_writer_protocol", return_value={}), \
                 patch.object(Reporting_Doctor.Shelly, "gather_citations", return_value={}), \
                 patch.object(Reporting_Doctor, "make_html_report"):
                result = Reporting_Doctor.reporter_protocol(copy.deepcopy(cfg))

            self.assertIn("madeByReporting", result["runtimeInfo"])
            self.assertTrue(result["runtimeInfo"]["madeByReporting"]["reportHtml"].endswith("drFrankenstein_report.html"))

    def test_process_fitting_results_uses_new_shuffle_filename(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            reporter_dir = os.path.join(tmp, "report")
            images_dir = os.path.join(reporter_dir, "Images")
            qmmm_dir = os.path.join(tmp, "qmmm")
            torsion_dir = os.path.join(qmmm_dir, "t1")
            os.makedirs(torsion_dir, exist_ok=True)
            os.makedirs(images_dir, exist_ok=True)
            cfg["runtimeInfo"]["madeByReporting"] = {"reporterDir": reporter_dir, "imagesDir": images_dir}
            cfg["runtimeInfo"]["madeByStitching"] = {
                "qmmmParameterFittingDir": qmmm_dir,
                "shufflesCompleted": 25,
                "maxTorsions": 4,
                "finalParameters": {"t1": "params"},
            }
            cfg["parameterFittingInfo"]["sagvolSmoothing"] = None
            cfg["parameterFittingInfo"]["l2DampingFactor"] = 0.1

            final_png = os.path.join(torsion_dir, "fitting_shuffle_4_nCosines_4.png")
            mae_png = os.path.join(torsion_dir, "mean_average_error.png")
            gif = os.path.join(torsion_dir, "torsion_fitting.gif")
            all_mae = os.path.join(qmmm_dir, "run_mean_average_error.png")
            for path, content in [(final_png, "latest"), (mae_png, "mae"), (gif, "gif"), (all_mae, "all")]:
                with open(path, "w") as handle:
                    handle.write(content)

            result = Reporting_Doctor.Reporting_Monster.process_fitting_results(cfg)

            self.assertTrue(result["fittingImages"]["t1"]["fittingPng"].endswith("t1_final.png"))
            self.assertTrue(result["fittingImages"]["t1"]["fittingGif"].endswith("t1_fitting.gif"))
            self.assertTrue(os.path.exists(os.path.join(images_dir, "fitting_images", "t1_final.png")))
            with open(os.path.join(images_dir, "fitting_images", "t1_final.png")) as handle:
                self.assertEqual(handle.read(), "latest")

    def test_process_fitting_results_keeps_gif_and_png_separate(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self.make_base_config(tmp)
            reporter_dir = os.path.join(tmp, "report")
            images_dir = os.path.join(reporter_dir, "Images")
            qmmm_dir = os.path.join(tmp, "qmmm")
            torsion_dir = os.path.join(qmmm_dir, "t1")
            os.makedirs(torsion_dir, exist_ok=True)
            os.makedirs(images_dir, exist_ok=True)
            cfg["runtimeInfo"]["madeByReporting"] = {"reporterDir": reporter_dir, "imagesDir": images_dir}
            cfg["runtimeInfo"]["madeByStitching"] = {
                "qmmmParameterFittingDir": qmmm_dir,
                "shufflesCompleted": 25,
                "maxTorsions": 4,
                "finalParameters": {"t1": "params"},
            }
            cfg["parameterFittingInfo"]["sagvolSmoothing"] = None
            cfg["parameterFittingInfo"]["l2DampingFactor"] = 0.1

            gif = os.path.join(torsion_dir, "torsion_fitting.gif")
            mae_png = os.path.join(torsion_dir, "mean_average_error.png")
            all_mae = os.path.join(qmmm_dir, "run_mean_average_error.png")
            for path, content in [(gif, "gif"), (mae_png, "mae"), (all_mae, "all")]:
                with open(path, "w") as handle:
                    handle.write(content)

            result = Reporting_Doctor.Reporting_Monster.process_fitting_results(cfg)

            self.assertIsNone(result["fittingImages"]["t1"]["fittingPng"])
            self.assertTrue(result["fittingImages"]["t1"]["fittingGif"].endswith("t1_fitting.gif"))
            self.assertTrue(os.path.exists(os.path.join(images_dir, "fitting_images", "t1_fitting.gif")))


class TestAmberCreation(unittest.TestCase):
    def test_extract_backbone_types_allows_missing_optional_h_and_ha(self):
        prmtop_df = pd.DataFrame({"name": ["N", "CA", "C", "O"], "type": ["N", "CT", "C", "O"]})
        aliases = {"N": ["N"], "CA": ["CA"], "C": ["C"], "O": ["O"], "H": ["H"], "HA": ["HA"]}

        with patch("Laboratory.Experiments.Protocol_7_Creation.AMBER_creation.pmd.load_file", return_value=MagicMock(to_dataframe=lambda: prmtop_df), create=True):
            atom_types = AMBER_creation.extract_backbone_types_from_prmtop("fake.prmtop", aliases)

        self.assertEqual(atom_types["N"], "N")
        self.assertEqual(atom_types["CA"], "CT")
        self.assertNotIn("H", atom_types)
        self.assertNotIn("HA", atom_types)

    def test_extract_backbone_types_still_requires_core_backbone_atoms(self):
        prmtop_df = pd.DataFrame({"name": ["N", "CA", "C"], "type": ["N", "CT", "C"]})
        aliases = {"N": ["N"], "CA": ["CA"], "C": ["C"], "O": ["O"]}

        with patch("Laboratory.Experiments.Protocol_7_Creation.AMBER_creation.pmd.load_file", return_value=MagicMock(to_dataframe=lambda: prmtop_df), create=True):
            with self.assertRaises(ValueError):
                AMBER_creation.extract_backbone_types_from_prmtop("fake.prmtop", aliases)


if __name__ == "__main__":
    unittest.main()
