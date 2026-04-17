
import os
from os import path as p
import pandas as pd
from pdbUtils.pdbUtils import pdb2df
from OperatingTools import electron_checker
## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass

def check_pdb(config: dict) -> None:
    inputDir = config["pathInfo"]["inputDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    charge = config["moleculeInfo"]["charge"]
    multiplicity = config["moleculeInfo"]["multiplicity"]
    molPdb = p.join(inputDir, f"{moleculeName}.pdb")
    pdbDf = pdb2df(molPdb)

    check_for_duplicate_atoms(pdbDf)
    electron_checker.validate_charge_multiplicity_from_pdb(
        pdb_df=pdbDf,
        charge=charge,
        multiplicity=multiplicity
    )


def check_for_duplicate_atoms(pdbDf: pd.DataFrame):
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
    
