## BASIC LIBRARIES ##
import os
from os import path as p
from subprocess import call, PIPE, run, STDOUT
import pandas as pd
import re
import numpy as np
from pdbUtils import pdbUtils
from shutil import move, copy
import random


## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass


from OperatingTools import file_parsers
# # 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
# def convert_traj_xyz_to_pdb(trajXyzs: list[FilePath],
#                              cappedPdb: FilePath,
#                                fittingRoundDir: DirectoryPath) -> list[FilePath]:
#     """
#     Effectively converts a list of XYZ files to PDB
#     In practice, updates a static PDB file with coords from a list of XYZ files

#     Args:
#         trajXyzs (list[FilePath]): a list of XYZ files
#         cappedPdb (FilePath): a PDB file containing correct atom data
#         fittingRoundDir (DirectoryPath): the output dir for this function
    
#     """
#     trajPdbs = []
#     for trajXyz in trajXyzs:
#         trajIndex = trajXyz.split(".")[1]
#         trajPdb = p.join(fittingRoundDir, f"orca.{trajIndex}.pdb")
#         update_pdb_coords(cappedPdb, trajXyz, trajPdb)
#         trajPdbs.append(trajPdb)

#     return trajPdbs




# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def shuffle_torsion_tags(torsion_tags: list[str]) -> list[str]:
    """
    Re-orders the torsion tags in a random order
    
    Args:
        torsion_tags (list[str]): the torsion tags associated with our torsions
    Returns:
        shuffled_torsion_tags: torsion_tags re-ordered
    """
    shuffledTorsionTags = torsion_tags # Variable name is already camelCase, matches return doc.
    random.shuffle(shuffledTorsionTags)

    return shuffledTorsionTags

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pdb2mol2(in_pdb: FilePath,
              out_mol2: FilePath,
                working_dir: DirectoryPath) -> None:
    """
    Uses antechamber to convert pdb to mol2
    Writes some unwanted temporary files
    TODO: clean these up

    Args:
        in_pdb (FilePath): input file in PDB format
        out_mol2 (FilePath): output file in MOL2 format
        working_dir (DirectoryPath): working directory (vital for cleanup)

    Returns:
        None (out_mol2 has already been defined!)
    
    """
    os.chdir(working_dir)

    ## set RES_ID to 1 for all atoms to keep antechamber happy
    pdbDf = pdbUtils.pdb2df(in_pdb)
    pdbDf["RES_ID"] = 1
    tmpPdb = p.join(working_dir, "tmp.pdb")
    pdbUtils.df2pdb(pdbDf, tmpPdb)
    ## get index, set path for antechamber to write outputs
    index: str = p.basename(in_pdb).split("_")[1].split(".")[0]
    antechamberOut: FilePath = p.join(working_dir, f"antechamber_{index}.out")
    ## run antechamber to create MOL2 file from PDB
    antechamberCommand: list = [
        "antechamber", "-i", tmpPdb, "-fi", "pdb", "-o", out_mol2,
        "-fo", "mol2", "-at", "gaff2", "-rn", "MOL", "-s", "2",
    ]
    with open(antechamberOut, 'w') as outfile:
        run(antechamberCommand, stdout=outfile, stderr=STDOUT)
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def edit_mol2_atom_types(in_mol2: FilePath,
                         out_mol2: FilePath) -> None: # Renamed function
    
    """
    Edits the atom types in a MOL2 file to make compatable with gaff2
    This may need to be added to in the case of new weird atom types
    from antechamber

    Args:
        in_mol2 (FilePath): input MOL2 file
        out_mol2 (FilePath): output MOL2 file

    Returns:
        None (out_mol2 is already defined!)
    """


    ## open in_mol2 for reading and out_mol2 for writing
    with open(in_mol2, 'r') as inMol2File, open(out_mol2, "w") as outMol2File: # Renamed file handlers
        ## loop through in_mol2 until we get to the atom section
        mol2Lines = inMol2File.readlines()
        isAtomLine = False
        for line in mol2Lines:
            if line.strip() == "@<TRIPOS>ATOM":
                isAtomLine = True
            elif line.strip() == "@<TRIPOS>BOND":
                isAtomLine = False
            elif isAtomLine:
                ## get atom name
                atomData = line.split()
                atomName = atomData[1]
                ## set alpha-carbon to CX rather than c3
                if atomName == "CA":
                    line = re.sub(r'(\s)c3(\s)', r'\1CX\2', line)
                ## other substitutions
                line = re.sub(r'(\s)os(\s)', r'\1o \2', line)
                line = re.sub(r'(\s)n7(\s)', r'\1ns\2', line)
                line = re.sub(r'(\s)h1(\s)', r'\1hc\2', line)
            outMol2File.write(line)
        ## loop through inMol2 until we get to the atom section
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def create_atom_type_map(in_mol2: FilePath) -> dict:
    """
    Creates a dict that maps atom names to atom types
    using a MOL2 file

    Args:
        in_mol2 (FilePath): input MOL2 file

    Returns:
        atom_type_map (dict): dict mapping atom names to atom types
    
    """
    ## init empty dict
    atomTypeMap = {}
    ## open in_mol2 for reading
    with open(in_mol2, 'r') as inMol2File: # Renamed file handler
        ## loop through in_mol2 until we get to the atoms section
        mol2Lines = inMol2File.readlines()
        isAtomLine = False
        for line in mol2Lines:
            if line.strip() == "@<TRIPOS>ATOM":
                isAtomLine = True
            elif line.strip() == "@<TRIPOS>BOND":
                isAtomLine = False
            elif isAtomLine:
                ## get atom name and type, add to dict
                atomData = line.split()
                atomTypeMap[atomData[1]] = atomData[5]
    return atomTypeMap

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def sort_out_directories(config: dict) -> dict:
    """
    Makes directories for the parameter fitting stage
    Adds these paths to the "pathInfo" section of config

    Args: 
        config (dict): contains all information needed for run
    
    Returns:
        config (dict): updated config with pathInfo of new dirs
    
    """
    outputDir: DirectoryPath = config["pathInfo"]["outputDir"]

    parameterFittingTopDir: DirectoryPath = p.join(outputDir, "06_parameter_fitting")
    os.makedirs(parameterFittingTopDir, exist_ok=True)
    config["runtimeInfo"]["madeByStitching"]["parameterFittingTopDir"] = parameterFittingTopDir

    moleculeParameterDir: DirectoryPath = p.join(parameterFittingTopDir, "molecule_parameters")
    os.makedirs(moleculeParameterDir, exist_ok=True)
    config["runtimeInfo"]["madeByStitching"]["moleculeParameterDir"] = moleculeParameterDir

    mmTotalCalculationDir: DirectoryPath = p.join(parameterFittingTopDir, "mm_total_energies")
    os.makedirs(mmTotalCalculationDir, exist_ok=True)
    config["runtimeInfo"]["madeByStitching"]["mmTotalCalculationDir"] = mmTotalCalculationDir

    qmmmFittingDir = p.join(parameterFittingTopDir, "qm-mm_parameter_fitting")
    config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"] = qmmmFittingDir
    os.makedirs(qmmmFittingDir,exist_ok=True)
    
    return config
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_completed_torsion_scan_dirs(config: dict, torsion_tag:str) -> List[DirectoryPath]:
    """
    Looks through ORCA scan directories 
    Works out whether scan completed
    Returns a list of completed scan dirs

    Args:
        config (dict): contains all information needed for run
        torsionTag (str): identifier for torsion (format "W-X-Y-Z")

    Returns:
        completedTorsionScanDirs (list): list of completed scan directories
    
    """
    ## init empty list to store scan dirs in
    completedTorsionScanDirs: list = []

    ## find dir where ORCA torsion scans were performed
    torsionTopDir: DirectoryPath = config["runtimeInfo"]["madeByTwisting"]["torsionDir"]

    torsionDir: DirectoryPath = p.join(torsionTopDir, f"torsion_{torsion_tag}")

    ## loop through batches of scans
    for conformerName in os.listdir(torsionDir): # conformerName is already camelCase
        if not "conformer" in conformerName:
            continue
        ## find conformer dir
        conformerDir: DirectoryPath = p.join(torsionDir, conformerName)
        for calculationName in os.listdir(conformerDir):
            ## get calculation dir (forwards and backwards)
            if  calculationName.endswith("forwards") or calculationName.endswith("backwards"):
                calculationDir: DirectoryPath = p.join(conformerDir, calculationName)
                orcaFinishedNormally = p.join(calculationDir, "ORCA_FINISHED_NORMALLY")
                if p.isfile(orcaFinishedNormally):
                    completedTorsionScanDirs.append(calculationDir)
    return completedTorsionScanDirs

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_scan_angles_from_orca_inp(scan_dir: DirectoryPath) -> np.ndarray:
    """
    Read through an ORCA input file and find start and end angles of a dihedral scan
    Creates a range between these

    Args:
        scanDir (DirectoryPath): location of orca_.inp

    Returns:
        angleRange (list): list of angles between start and end angles
    """
    ## find orca input
    orcaInp = p.join(scan_dir, "orca_scan.inp")
    if not p.isfile(orcaInp):
        raise FileNotFoundError(f"orca_scan.inp not found in {scan_dir}")
    ## open orca input for reading
    with open(orcaInp, "r") as f: # f is already short and conventional
        for line in f:
            ## find dihedral scan, get start, end angle
            if line.startswith("D"):
                lineData = line.split()
                startAngle = float(lineData[6][:-1])
                endAngle = float(lineData[7][:-1])
                break
    ## create range
    if startAngle < endAngle:
        stepSize = 10
    else:
        stepSize = -10
    angleRange = np.arange(startAngle, endAngle + stepSize, stepSize)
    
    return angleRange
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rescale_angles_0_360(angle: np.ndarray) -> np.ndarray:
    """
    puts angles on a -180 to +180 scale
    """
    angle = angle % 360  # First, reduce the angle to the 0-360 range
    # if angle > 180:
    #     angle -= 360  # Shift angles greater than 180 to the negative side
    return angle
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def merge_energy_dfs(energy_dfs: List[Optional[pd.DataFrame]]) -> pd.DataFrame:
    """
    merges energy DataFrames on the Angle column

    Args:
        energyDfs (list): list of DataFrames with Angle column
    
    Returns:
        mergedDf (pd.DataFrame): merged DataFrame
    """
    ##TODO: is this efficient?
    energyDfsValid = [df for df in energy_dfs if df is not None] # Renamed to avoid confusion

    if not energyDfsValid:
        # Or handle as an error: raise ValueError("No valid DataFrames to merge.")
        return pd.DataFrame(columns=['Angle']) 


    mergedDf = energyDfsValid[0][['Angle', 'Energy']].rename(
        columns={'Energy': 'Energy_0'}
    )

    for i, df in enumerate(energyDfsValid[1:], start=1):
        mergedDf = mergedDf.merge(
            df[['Angle', 'Energy']],
            on='Angle',
            how='outer',
            suffixes=('', f'_df{i}')
        ).rename(columns={'Energy': f'Energy_{i}'})
    mergedDf = mergedDf.groupby('Angle', as_index=False).first()

    return mergedDf
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def update_pdb_coords(in_pdb: FilePath, xyz_file: FilePath, out_pdb: FilePath) -> None:
    """
    updates PDB file with XYZ coords
    
    Args:
        in_pdb (FilePath): input PDB file
        xyz_file (FilePath): input XYZ file
        out_pdb (FilePath): output PDB file

    Returns:
        None (out_pdb already defined!)
    """

    inDf = pdbUtils.pdb2df(in_pdb)
    xyzDf = file_parsers.xyz2df(xyz_file)

    inDf["X"] = xyzDf["x"]
    inDf["Y"] = xyzDf["y"]
    inDf["Z"] = xyzDf["z"]

    pdbUtils.df2pdb(inDf, out_pdb)
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
