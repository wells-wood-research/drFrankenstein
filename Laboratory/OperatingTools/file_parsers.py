## BASIC IMPORTS ##
from os import path as p
import os
import pandas as pd


## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass


def xyz2df(xyzFile: FilePath) -> pd.DataFrame:
    """
    converts an xyz file to a pd.DataFrame

    Args:
        xyzFile (FilePath): input XYZ file

    Returns:
        xyzDf (pd.DataFrame): dataframe containing index, element and coords
    """

    xyzData = []
    atomIndex = 0
    with open(xyzFile, "r") as f:
        lines = f.readlines()[2:]
        for line in lines:
            atomIndex +=1
            lineData = line.split()
            xyzData.append({"index": atomIndex,
                            "element": lineData[0],
                              "x": float(lineData[1]), 
                              "y": float(lineData[2]),
                                "z": float(lineData[3])})
    return pd.DataFrame(xyzData)



def parse_mol2(mol2File):
    atomData = []
    bondData = []
    readingAtoms = False
    readingBonds = False

    with open(mol2File, "r") as mol2:
        for line in mol2:
            if line.strip() == "":
                continue
            if line.startswith("@<TRIPOS>ATOM"):
                readingAtoms=True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                readingBonds=True
                readingAtoms=False
                continue
            if line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                break
            if readingAtoms:
                atomData.append(line.split())
            elif readingBonds:
                bondData.append(line.split())


        
    atomDataColumns = ["ATOM_ID", "ATOM_NAME", "X", "Y", "Z", "ATOM_TYPE", "RES_ID", "RES_NAME", "CHARGE"]
    atomDf = pd.DataFrame(atomData, columns=atomDataColumns)
    bondDataColumns = ["BOND_ID", "ATOM_A_ID", "ATOM_B_ID", "BOND_ORDER"]
    bondDf = pd.DataFrame(bondData, columns=bondDataColumns)

    return atomDf, bondDf