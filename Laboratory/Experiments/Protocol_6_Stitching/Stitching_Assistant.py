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
from functools import reduce
from scipy.signal import find_peaks


## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass


from OperatingTools import file_parsers

def construct_final_params(config: dict, paramFile: FilePath, param_extraction_function: object, torsionTags: list[str]) -> dict:
    """Extract final torsion parameters for all torsion tags."""
    finalParameters = {}
    for torsionTag in torsionTags:
        torsionParameters = param_extraction_function(paramFile, torsionTag, config)
        finalParameters[torsionTag] = torsionParameters


    return finalParameters


def score_flatline_status(
    scoresDf: pd.DataFrame,
    torsionTag: str = "All_Torsions",
    windowSize: int = 3,
    diffTolerance: float = 0.1,
) -> tuple[bool, bool]:
    """Return whether torsion and total scores have flatlined."""
    if "torsion_tag" in scoresDf.columns:
        scoresDf = scoresDf[scoresDf["torsion_tag"] == torsionTag]
    if scoresDf.empty:
        return False, False

    torsionCol = "fit_score_torsion" if "fit_score_torsion" in scoresDf.columns else "mae_torsion"
    totalCol = "fit_score_total" if "fit_score_total" in scoresDf.columns else "mae_total"

    diffTorsion = scoresDf[torsionCol].diff().abs().dropna()
    diffTotal = scoresDf[totalCol].diff().abs().dropna()

    torsionFlatLined = len(diffTorsion) >= windowSize and np.all(diffTorsion.tail(windowSize) < diffTolerance)
    totalFlatLined = len(diffTotal) >= windowSize and np.all(diffTotal.tail(windowSize) < diffTolerance)

    return torsionFlatLined, totalFlatLined


def check_scores_for_flatline(torsionTags: list[str], convergedTags: list[str], fittingScoresCsv: FilePath, windowSize: int = 3, diffTolerance: float = 0.1) -> list[str]:
    """Identify torsions whose scores have flatlined."""
    
    unconvergedTorsionTags = [tag for tag in torsionTags if tag not in convergedTags]

    scoresDf = pd.read_csv(fittingScoresCsv)

    flatlinedTorsions = []
    for torsionTag in unconvergedTorsionTags:
        torsionFlatLined, totalFlatLined = score_flatline_status(
            scoresDf,
            torsionTag=torsionTag,
            windowSize=windowSize,
            diffTolerance=diffTolerance,
        )
        if torsionFlatLined and totalFlatLined:
            flatlinedTorsions.append(torsionTag)
    
    return flatlinedTorsions


        



def rms_of_mae_dict(meanAverageErrors: dict[str, list[float]]) -> float:
    """Calculate the RMS of the MAE history."""
    rms = np.sqrt(sum(meanAverageErrors[torsion][-1] ** 2 for torsion in meanAverageErrors) / len(meanAverageErrors.keys()))
    return rms
def init_tqdm_bar_options() -> dict:
    """Return the shared tqdm configuration for stitching progress bars."""
    tqdmBarOptions = {
        "desc": f"\033[32mRunning Parameter Fitting\033[0m",
        "ascii": "-🗲→",    
        "colour": "yellow",
        "unit":  "scan",
        "dynamic_ncols": False, 
        "ncols": 102,
        "leave": True,
        "position": 0
    }
    return tqdmBarOptions

def shuffle_torsion_tags(torsionTags: list[str], maxShuffles: int, seed: int) -> list[str]:
    """Return repeated shuffled torsion tag orders."""
    shuffledTorsionTags = []
    baseTorsionTags = list(torsionTags)
    for i in range(maxShuffles):
        shuffleIndex = i + 1
        # Use a per-shuffle derived seed so each shuffle order is generated independently.
        derivedSeed = seed * shuffleIndex if seed != 0 else shuffleIndex
        workingTorsionTags = list(baseTorsionTags)
        random.Random(derivedSeed).shuffle(workingTorsionTags)
        shuffledTorsionTags.extend(workingTorsionTags)

    return shuffledTorsionTags

def check_mae_convergence(latestRmsMaeTorsion: float, latestRmsMaeTotal: float, converganceTolerance: float | None) -> bool:
    """Check whether MAE scores are below the convergence tolerance."""
    if converganceTolerance is None:
        return False
    else:
        return latestRmsMaeTorsion < converganceTolerance and latestRmsMaeTotal < converganceTolerance


def check_torsion_convergence(latestTorsionScore: float, latestTotalScore: float, converganceTolerance: float | None) -> bool:
    """Check whether torsion scores are below the convergence tolerance."""
    if converganceTolerance is None:
        return False
    return latestTorsionScore < converganceTolerance and latestTotalScore < converganceTolerance


def _stationary_point_angles(signal: np.ndarray, sampleSpacingDegrees: int = 10) -> np.ndarray:
    """Return the angular positions of stationary points in a signal."""
    signal = np.asarray(signal, dtype=float)
    maxima, _ = find_peaks(signal)
    minima, _ = find_peaks(-signal)
    stationaryIndexes = np.unique(np.concatenate([maxima, minima]))
    return stationaryIndexes.astype(float) * sampleSpacingDegrees


def _circular_angle_distance(angleA: float, angleB: float) -> float:
    """Return the shortest distance between two angles."""
    delta = abs(angleA - angleB) % 360.0
    return min(delta, 360.0 - delta)


def calculate_profile_fit_score(
    qmEnergy: np.ndarray,
    mmEnergy: np.ndarray,
    sampleSpacingDegrees: int = 10,
) -> float:
    """Return the composite score for a QM/MM profile comparison."""
    return calculate_profile_fit_metrics(qmEnergy, mmEnergy, sampleSpacingDegrees)["composite_score"]


def calculate_profile_fit_metrics(
    qmEnergy: np.ndarray,
    mmEnergy: np.ndarray,
    sampleSpacingDegrees: int = 10,
) -> dict:
    """Return the detailed metric breakdown for a QM/MM profile comparison."""
    qmEnergy = np.asarray(qmEnergy, dtype=float)
    mmEnergy = np.asarray(mmEnergy, dtype=float)
    if qmEnergy.shape != mmEnergy.shape:
        raise ValueError("QM and MM energy arrays must have the same shape.")

    qmEnergy = qmEnergy - qmEnergy.min()
    mmEnergy = mmEnergy - mmEnergy.min()
    qmAmplitude = float(np.ptp(qmEnergy))
    mmAmplitude = float(np.ptp(mmEnergy))
    energyScale = max(qmAmplitude, np.finfo(float).eps)
    qmIsFlat = qmAmplitude <= np.finfo(float).eps
    mmIsFlat = mmAmplitude <= np.finfo(float).eps

    qmStationaryAngles = _stationary_point_angles(qmEnergy, sampleSpacingDegrees)
    mmStationaryAngles = _stationary_point_angles(mmEnergy, sampleSpacingDegrees)

    if qmIsFlat and mmIsFlat:
        locationScore = 0.0
        amplitudeScore = 0.0
        stationaryCountScore = 0.0
        normalizedMaeScore = 0.0
    elif qmIsFlat != mmIsFlat:
        locationScore = 1.0
        amplitudeScore = 1.0
        stationaryCountScore = 1.0
        normalizedMaeScore = 1.0
    else:
        amplitudeScore = abs(qmAmplitude - mmAmplitude) / energyScale
        normalizedMaeScore = float(np.mean(np.abs(qmEnergy - mmEnergy)) / energyScale)

        if len(qmStationaryAngles) == 0 or len(mmStationaryAngles) == 0:
            locationScore = 1.0
            stationaryCountScore = 1.0
        else:
            matchedPoints = min(len(qmStationaryAngles), len(mmStationaryAngles))
            locationScore = float(
                np.mean(
                    [
                        _circular_angle_distance(qmAngle, mmAngle) / 180.0
                        for qmAngle, mmAngle in zip(qmStationaryAngles[:matchedPoints], mmStationaryAngles[:matchedPoints])
                    ]
                )
            )
            stationaryCountScore = abs(len(qmStationaryAngles) - len(mmStationaryAngles)) / max(
                len(qmStationaryAngles), len(mmStationaryAngles), 1
            )

    return {
        "composite_score": float(np.mean([locationScore, amplitudeScore, stationaryCountScore, normalizedMaeScore])),
        "location_score": locationScore,
        "amplitude_score": amplitudeScore,
        "stationary_count_score": stationaryCountScore,
        "normalized_mae_score": normalizedMaeScore,
        "qm_stationary_count": int(len(qmStationaryAngles)),
        "mm_stationary_count": int(len(mmStationaryAngles)),
        "qm_amplitude": qmAmplitude,
        "mm_amplitude": mmAmplitude,
    }

   
        



def remove_exploded_torsions(config: dict) -> list[str]:
    """Return torsion tags that remain after filtering exploded scans."""

    torsionTags = []
    torsionsToScan = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]
    for torsionTag, torsionData in torsionsToScan.items():
        if torsionData["globalMinimaAngle"] is None:
            continue
        torsionTags.append(torsionTag)
    return torsionTags



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
        
    ## clean up temporary files
    os.remove(tmpPdb)
    filesToRemove = [f for f in os.listdir(workingDir) if f.startswith("ANTECHAMBER")]
    for f in filesToRemove:
        os.remove(p.join(workingDir, f))
    return None
        
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
def get_completed_torsion_scan_dirs(config: dict, torsionTag: str) -> list[str]:
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
def get_scan_angles_from_orca_inp(scanDir: DirectoryPath) -> list[float]:
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
def merge_energy_dfs(energyDfs: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Merges a list of energy DataFrames efficiently by concatenating first.

    Args:
        energyDfs (list): A list of DataFrames, each with "Angle" and "Energy" columns.

    Returns:
        pd.DataFrame: A merged DataFrame with a single "Angle" column and
                      multiple "Energy_i" columns.
    """
    ## Remove empty or None DataFrames
    energyDfs = [df for df in energyDfs if df is not None and not df.empty]
    if not energyDfs:
        return pd.DataFrame()

    # Assign an ID to each DataFrame to track its origin, then concatenate.
    # The 'keys' argument in concat adds a new level to the index.
    combined_df = pd.concat(
        [df[["Angle", "Energy"]] for df in energyDfs],
        keys=range(len(energyDfs)),
        names=['df_index', 'original_index']
    ).reset_index(level='original_index', drop=True)

    # Perform a single groupby to get the first energy for each angle per original df
    grouped = combined_df.groupby(['df_index', 'Angle'])['Energy'].first()

    # Unstack the result to pivot the df_index into columns
    merged_df = grouped.unstack('df_index')

    # Rename columns to the desired "Energy_i" format
    merged_df.columns = [f"Energy_{i}" for i in merged_df.columns]

    # Reset the index to turn the 'Angle' index back into a column
    return merged_df.reset_index()

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



def _remove_converged_torsion_tags(fittingScoresCsv: FilePath, converganceTolerance: float | None, torsionTags: list[str], topN: int = 5) -> list[str]:
    """
    Reads through the MAE csv file to find torsions that have converged
    Removes these from the list of torsions to fit

    Args:
        fittingScoresCsv (FilePath): location of MAE csv file
        converganceTolerance (float): tolerance for convergence
        torsionTags (list): list of torsion tags to fit
    Returns:
        list: updated list of torsion tags without converged ones
    """
    maeDf = pd.read_csv(fittingScoresCsv, index_col=False, header=0)
    maxShuffle = maeDf["shuffle"].max()
    lastShuffleDf = maeDf[(maeDf["shuffle"] == maxShuffle) & (maeDf["torsion_tag"] != "All_Torsions")]
    ## sort by MAE total
    totalCol = "fit_score_total" if "fit_score_total" in lastShuffleDf.columns else "mae_total"
    lastShuffleDf = lastShuffleDf.sort_values(totalCol)
    ## get top N torsions by MAE total
    nonConvergedTorsions = lastShuffleDf[lastShuffleDf[totalCol] > converganceTolerance]["torsion_tag"].tolist()

    nonConvergedTorsions = list(set(nonConvergedTorsions))

    return nonConvergedTorsions
