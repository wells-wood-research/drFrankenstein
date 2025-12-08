
import os
from os import path as p
import pandas as pd
from pdbUtils.pdbUtils import pdb2df
## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass

def check_pdb(config: dict) -> None:
    inputDir = config["pathInfo"]["inputDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    molPdb = p.join(inputDir, f"{moleculeName}.pdb")
    pdbDf = pdb2df(molPdb)

    check_for_duplicate_atoms(pdbDf)


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
    
