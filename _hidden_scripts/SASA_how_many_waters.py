import os
from os import path as p
import pandas as pd
import mdtraj as md
from pdbUtils import pdbUtils   
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass
###########################################################
def inputs():
    conformerDir = "/home/esp/scriptDevelopment/drFrankenstein/07_TTC_outputs/03_GOAT_conformers"
    conformerXyzs = [p.join(conformerDir, file) for file in os.listdir(conformerDir) if file.endswith(".xyz") and "conformer" in file]
    cappedPdb = "/home/esp/scriptDevelopment/drFrankenstein/07_TTC_outputs/01_termini_capping/TCC_capped.pdb"
    outDir = "/home/esp/scriptDevelopment/drFrankenstein/_hidden_scripts/SASA_how_many_waters"

    return conformerXyzs, cappedPdb, outDir

###########################################################
def update_pdb_coords(inPdb: FilePath, xyzFile: FilePath, outPdb: FilePath) -> None:
    """
    updates PDB file with XYZ coords
    
    Args:
        inPdb (FilePath): input PDB file
        xyzFile (FilePath): input XYZ file
        outPdb (FilePath): output PDB file

    Returns:
        None (outPdb already defined!)
    """

    inDf = pdbUtils.pdb2df(inPdb)
    xyzDf = xyz2df(xyzFile)

    inDf["X"] = xyzDf["x"]
    inDf["Y"] = xyzDf["y"]
    inDf["Z"] = xyzDf["z"]

    pdbUtils.df2pdb(inDf, outPdb)

###########################################################
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

###########################################################
def convert_traj_xyz_to_pdb(trajXyzs: list[FilePath],
                             cappedPdb: FilePath,
                               fittingRoundDir: DirectoryPath) -> list[FilePath]:
    """
    Effectively converts a list of XYZ files to PDB
    In practice, updates a static PDB file with coords from a list of XYZ files

    Args:
        trajXyzs (list[FilePath]): a list of XYZ files
        cappedPdb (FilePath): a PDB file containing correct atom data
        fittingRoundDir (DirectoryPath): the output dir for this function
    
    """
    trajPdbs = []
    for i, trajXyz in enumerate(trajXyzs):
        trajPdb = p.join(fittingRoundDir, f"conformer_{i}.pdb")
        update_pdb_coords(cappedPdb, trajXyz, trajPdb)
        trajPdbs.append(trajPdb)

    return trajPdbs
###########################################################
def main():
    conformerXyzs, cappedPdb, outDir = inputs()
    os.makedirs(outDir, exist_ok=True)

    conformerPdbs = convert_traj_xyz_to_pdb(conformerXyzs, cappedPdb, outDir)

    totalSasas = []
    for pdbFile in conformerPdbs:
        traj = md.load(pdbFile)
        sasaPerAtom = md.shrake_rupley(traj)
        totalSasa = sasaPerAtom.sum()
        totalSasas.append(totalSasa)
    
    meanTotalSasa = sum(totalSasas)/len(totalSasas)

    nWatersPerNmSquared = 10

    nWaters = round(meanTotalSasa * nWatersPerNmSquared)

    print(nWaters)



###########################################################



if __name__ == "__main__":
    main()