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