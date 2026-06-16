
import os
from os import path as p
import pandas as pd
from pdbUtils.pdbUtils import pdb2df
from OperatingTools import electron_checker
from typing import Any
## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass

def check_pdb(config: dict) -> None:
    """Validate the input PDB and reject duplicate atom names."""
    inputDir = config["pathInfo"]["inputDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    molPdb = p.join(inputDir, f"{moleculeName}.pdb")
    pdbDf = pdb2df(molPdb)

    check_for_duplicate_atoms(pdbDf)


def validate_charge_multiplicity(config: dict, pdb_path: FilePath) -> None:
    """Verify that the PDB electron count matches the configured charge and multiplicity."""
    charge = config["moleculeInfo"]["charge"]
    multiplicity = config["moleculeInfo"]["multiplicity"]
    pdbDf = pdb2df(pdb_path)
    electron_checker.validate_charge_multiplicity_from_pdb(
        pdb_df=pdbDf,
        charge=charge,
        multiplicity=multiplicity
    )


def get_capped_pdb_path(config: dict) -> FilePath:
    """Return the capped PDB path recorded in runtimeInfo."""
    runtime_info = config.get("runtimeInfo", {})
    made_by_capping = runtime_info.get("madeByCapping", {})
    capped_pdb = made_by_capping.get("cappedPdb")
    if not capped_pdb:
        raise ValueError(
            "Capped PDB not found in runtimeInfo.madeByCapping.cappedPdb. "
            "Ensure the capping protocol has completed before running the electron check."
        )
    return capped_pdb


def check_for_duplicate_atoms(pdbDf: pd.DataFrame) -> None:
    """Raise an error if the PDB contains duplicate atom names."""
    atomNamesChecked = []
    duplicateAtomNames = []
    for _, row in pdbDf.iterrows():
        if row["ATOM_NAME"] in atomNamesChecked:
            duplicateAtomNames.append(row["ATOM_NAME"])
        else:
            atomNamesChecked.append(row["ATOM_NAME"])

    if len(duplicateAtomNames) > 0:
        duplicateReport = {}
        for atomName in duplicateAtomNames:
            duplicateIds = pdbDf[pdbDf["ATOM_NAME"] == atomName]["ATOM_ID"].tolist()
            duplicateReport[atomName] = duplicateIds
        raise ValueError(f"Duplicate atoms found in PDB file: {duplicateReport}")
    
