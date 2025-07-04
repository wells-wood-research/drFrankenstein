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
def remove_exploded_torsions(config: dict) -> None:
    """
    TODO: write this
    """

    torsionTags = []
    torsionsToScan = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]
    for torsionTag, torsionData in torsionsToScan.items():
        if torsionData["globalMinimaAngle"] is None:
            continue
        torsionTags.append(torsionTag)
    return torsionTags



# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def shuffle_torsion_tags(torsionTags: list[str]) -> list[str]:
    """
    Re-orders the torsion tags in a random order
    
    Args:
        torsionTags (list[str]): the torsion tags associated with our torsions
    Returns:
        shuffledTorsionTags: torsionTags re-ordered
    """
    shuffledTorsionTags = torsionTags
    random.shuffle(shuffledTorsionTags)

    return shuffledTorsionTags

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pdb2mol2(inPdb: FilePath,
              outMol2: FilePath,
                workingDir: DirectoryPath) -> None:
    """
    Uses antechamber to convert pdb to mol2
    Writes some unwanted temporary files
    TODO: clean these up

    Args:
        inPdb (FilePath): input file in PDB format
        outPdb (FilePath): output file in MOL2 format
        workingDir (DirectoryPath): working directory (vital for cleanup)

    Returns:
        None (outMol2 has already been defined!)
    
    """
    os.chdir(workingDir)

    ## set RES_ID to 1 for all atoms to keep antechamber happy
    pdbDf = pdbUtils.pdb2df(inPdb)
    pdbDf["RES_ID"] = 1
    tmpPdb = p.join(workingDir, "tmp.pdb")
    pdbUtils.df2pdb(pdbDf, tmpPdb)
    ## get index, set path for antechamber to write outputs
    index: str = p.basename(inPdb).split("_")[1].split(".")[0]
    antechamberOut: FilePath = p.join(workingDir, f"antechamber_{index}.out")
    ## run antechamber to create MOL2 file from PDB
    antechamberCommand: list = [
        "antechamber", "-i", tmpPdb, "-fi", "pdb", "-o", outMol2,
        "-fo", "mol2", "-at", "gaff2", "-rn", "MOL", "-s", "2",
    ]
    with open(antechamberOut, 'w') as outfile:
        run(antechamberCommand, stdout=outfile, stderr=STDOUT)
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def edit_mo2_atom_types(inMol2: FilePath,
                         outMol2: FilePath) -> None:
    
    """
    Edits the atom types in a MOL2 file to make compatable with gaff2
    This may need to be added to in the case of new weird atom types
    from antechamber

    Args:
        inMol2 (FilePath): input MOL2 file
        outMol2 (FilePath): output MOL2 file

    Returns:
        None (outMol2 is already defined!)
    """


    ## open inMol2 for reading and outMol2 for writing
    with open(inMol2, 'r') as inMol2, open(outMol2, "w") as outMol2:
        ## loop through inMol2 until we get to the atom section
        mol2Lines = inMol2.readlines()
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
            outMol2.write(line)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def create_atom_type_map(inMol2: FilePath) -> dict:
    """
    Creates a dict that maps atom names to atom types
    using a MOL2 file

    Args:
        inMol2 (FilePath): input MOL2 file

    Returns:
        atomTypeMape (dict): dict mapping atom names to atom types
    
    """
    ## init empty dict
    atomTypeMap = {}
    ## open inMol2 for reading
    with open(inMol2, 'r') as inMol2:
        ## loop through inMol2 until we get to the atoms section
        mol2Lines = inMol2.readlines()
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
def sort_out_directories(config) -> dict:
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
def get_completed_torsion_scan_dirs(config: dict, torsionTag:str):
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

    torsionDir: DirectoryPath = p.join(torsionTopDir, f"torsion_{torsionTag}")

    ## loop through batches of scans
    for conformerName in os.listdir(torsionDir):
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
def get_scan_angles_from_orca_inp(scanDir: DirectoryPath):
    """
    Read through an ORCA input file and find start and end angles of a dihedral scan
    Creates a range between these

    Args:
        scanDir (DirectoryPath): location of orca_.inp

    Returns:
        angleRange (list): list of angles between start and end angles
    """
    ## find orca input
    orcaInp = p.join(scanDir, "orca_scan.inp")
    if not p.isfile(orcaInp):
        raise FileNotFoundError(f"orca_scan.inp not found in {scanDir}")
    ## open orca input for reading
    with open(orcaInp, "r") as f:
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
def rescale_angles_0_360(angle: np.array) -> np.array:
    """
    puts angles on a -180 to +180 scale
    """
    angle = angle % 360  # First, reduce the angle to the 0-360 range
    # if angle > 180:
    #     angle -= 360  # Shift angles greater than 180 to the negative side
    return angle
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def merge_energy_dfs(energyDfs: list) -> pd.DataFrame:
    """
    merges energy DataFrames on the Angle column

    Args:
        energyDfs (list): list of DataFrames with Angle column
    
    Returns:
        mergedDf (pd.DataFrame): merged DataFrame
    """

    ##TODO: is this efficient?
    energyDfs = [df for df in energyDfs if not df is None]



    mergedDf = energyDfs[0][['Angle', 'Energy']].rename(
        columns={'Energy': 'Energy_0'}
    )

    for i, df in enumerate(energyDfs[1:], start=1):
        mergedDf = mergedDf.merge(
            df[['Angle', 'Energy']],
            on='Angle',
            how='outer',
            suffixes=('', f'_df{i}')
        ).rename(columns={'Energy': f'Energy_{i}'})
    mergedDf = mergedDf.groupby('Angle', as_index=False).first()

    return mergedDf
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
    xyzDf = file_parsers.xyz2df(xyzFile)

    inDf["X"] = xyzDf["x"]
    inDf["Y"] = xyzDf["y"]
    inDf["Z"] = xyzDf["z"]

    pdbUtils.df2pdb(inDf, outPdb)
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
