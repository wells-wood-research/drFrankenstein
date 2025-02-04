import os
from pathlib import Path
from typing import Optional

import rdkit
import psiresp

from parameteriser.prep import build_molecule_from_smiles
#from parameteriser.run_resp import run_psiresp


def in_notebook():
    # Simploe functtion to know if the code is executed from  nmotebook
    try:
        from IPython import get_ipython

        if "IPKernelApp" not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True


class FF_Genenerator:
    def __init__(self, input_string, config):
        # config=Optional(str | Path)):
        self.input_string: str = input_string

        mol, _, backbone_list, capping_list, sidechain_list = build_molecule_from_smiles(input_string)
        self.mol: rdkit.Molecule = mol
        self.backbone_list: list = backbone_list
        self.capping_list: list = capping_list
        self.sidechain_list: list = sidechain_list
        self.psiresp_job: psiresp.job.Job

        self.charges: list

        #self.conformers: 
        # config: str | Path

    def run_qm_resp(self):
        mol, _, backbone_list, capping_list = build_molecule_from_smiles(
            self.input_string
        )

        constraint_capping = run_psiresp(mol, backbone_list, capping_list)
        without_constraint = run_psiresp(mol, backbone_list, [])
        print(constraint_capping)
        atom_names = [
            f"{x}{i}"
            for i, x in enumerate([atom.GetSymbol() for atom in mol.GetAtoms()], 1)
        ]

    def run_resp(self, run_qm):
        constraint_capping = run_psiresp(
            self.mol, self.backbone_list, self.capping_list
        )
        print(f"running the QM now { run_qm}")
        if run_qm == True:
            os.chdir("resp_qm_calculations/single_point/")
            os.system("bash run_single_point.sh")

            constraint_capping = run_psiresp(
                self.mol, self.backbone_list, self.capping_list
            )
        self.charges = constraint_capping


    def get_resp_job_eva(self): 
        psirespmol = psiresp.Molecule.from_rdkit(self.mol)
        constraints = psiresp.ChargeConstraintOptions(symmetric_atoms_are_equivalent=True)
        constraints.add_charge_sum_constraint_for_molecule(psirespmol, charge=0, indices=self.capping_list)

        psirespmol.optimize_geometry=True

        geometry_options = psiresp.QMGeometryOptimizationOptions(
        method="hf",
        basis="6-31g*")

        esp_options = psiresp.QMEnergyOptions(
        method="hf",
        basis="6-31g*")

        job_multi = psiresp.Job(
            molecules=[psirespmol],
            charge_constraints=constraints,
            qm_optimization_options=geometry_options,
            qm_esp_options=esp_options,
            working_directory="resp_qm_calculations")
        
        self.psiresp_job = job_multi

    def run_resp_job_eva(self, resp_job):
        resp_job.generate_conformers()
        resp_job.generate_orientations()
        resp_job.run()

    def run_resp_optimisation_eva(self, resp_job):
        os.chdir("resp_qm_calculations/optimization/")
        os.system("bash run_optimization.sh")
        os.chdir("../../")
        resp_job.run()
    
    def run_qm_get_charges_eva(self, resp_job):
        os.chdir("resp_qm_calculations/single_point/")
        os.system("bash run_single_point.sh")
        os.chdir("../../")
        resp_job.run()
        self.charges = resp_job.molecules[0].stage_2_restrained_charges
        

    def produce_library(self):
        pass

    def prepare_dih_fit(self):
        pass

    def run_qm_dih_fit(self):
        pass

    def run_fit(self):
        pass
