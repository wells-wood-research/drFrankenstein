## BASIC IMPORTS ##
import os
from os import path as p
from pdbUtils import pdbUtils
from subprocess import call
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import glob 
import re
from shutil import rmtree

## MULTIPROCESSING AND LOADING BAR LIBRARIES ##
from mpire import WorkerPool
from mpire.utils import make_single_arguments

## RDKIT IMPORTS ##
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolToPDBFile
RDLogger.DisableLog('rdApp.warning')


## drFRANKENSTEIN LIBRARIES ##
import drPlotter
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

from typing import List, Tuple, Dict

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def twist_protocol(config):

    config = set_up_directories(config)

    ## identify rotatable bonds to be scanned
    rotatableBonds: Dict[List[Tuple[int,int,int,int]]] = identify_rotatable_bonds(config["moleculeInfo"]["cappedPdb"])

    ## loop over torsions that need to be scanned
    config["pathInfo"]["torsionDirs"] = []
    for rotatableBond in rotatableBonds:
        config = run_torsion_scanning(rotatableBond, config)


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_torsion_scanning(rotatableBond, config, debug=False) -> dict:
    print(f"Running torsion scans for {rotatableBond['atoms']}")
    ## make a top level dir for this torsion
    torsionId: str = "-".join(map(str, rotatableBond["atoms"])) 
    torsionDir = p.join(config["pathInfo"]["torsionTopDir"], f"torsion_{torsionId}" )
    os.makedirs(torsionDir, exist_ok=True)
    ## add to config
    config["pathInfo"]["torsionDirs"].append(torsionDir)
    ## init some variables
    meanAverageError = np.inf
    batchIndex = 1
    meanAverageErrors = []
    scanRawDfs = {}
    ## unpack torsionScanInfo
    conversionTolerance = config["torsionScanInfo"]["convergenceTolerance"]
    minScanBatches = config["torsionScanInfo"]["minScanBatches"]
    maxScanBatches = config["torsionScanInfo"]["minScanBatches"]
    ## INIT A WHILE LOOP - BREAK WHEN CONVERGENCE IS REACHED
    while (meanAverageError > conversionTolerance or batchIndex < minScanBatches) and (batchIndex + 1 != maxScanBatches):

        torsionBatchDir = p.join(torsionDir, f"batch_{batchIndex}")
        os.makedirs(torsionBatchDir, exist_ok=True)
        ## generate a batch of conformers
        conformerPdbs = gen_conformers(inputPdb=config["moleculeInfo"]["cappedPdb"],
                                        batchDir=torsionBatchDir,
                                        iteration=batchIndex,
                                        nConformers=config["hardwareInfo"]["nCores"])
        
        
        ## CREATE A DIR FOR THIS TORSION ANGLE
        torsionScanDir = p.join(torsionBatchDir, f"torsion_scans")
        os.makedirs(torsionScanDir, exist_ok=True)
        
        ## INIT EMPTY LIST TO STORE SCAN DATA
        scanDfs = []
        if debug:
            ## run in serial
            scanDfs = scan_in_serial(scanDfs, torsionScanDir, conformerPdbs, rotatableBond["indices"], batchIndex, config)
        else:
            # run torsion scans in paralell
            scanDfs = scan_in_parallel(scanDfs, torsionScanDir, conformerPdbs, rotatableBond["indices"], batchIndex, config)

        if not len(scanDfs) == 0:
            ## Merge scan data, calculate averages, rolling averages and mean average errors
            meanAverageError, mergedDf, scanAverageDf, rollingAverageDf = process_scan_data(scanDfs, torsionDir, batchIndex)
            meanAverageErrors.append(meanAverageError)
            scanRawDfs[batchIndex] = mergedDf
        ## STEP BATCH INDEX
        batchIndex += 1
    else:
        ## once convergence criteria has been met
        ## create a dataframe containing final scan energies
        finalScanEnergiesDf = pd.DataFrame()
        finalScanEnergiesDf["Angle"] = scanAverageDf["Angle"]
        finalScanEnergiesDf[torsionId] = rollingAverageDf.iloc[:, -1]
        ## write to file
        dataDir = p.join(torsionDir, "scan_data")
        finalScanEnergiesDf.to_csv(p.join(dataDir, "final_scan_energies.csv"), index=False)
        ## plot scan data
        drPlotter.plot_torsion_scans(torsionDir, scanRawDfs, scanAverageDf, rollingAverageDf, meanAverageErrors)

    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def set_up_directories(config: dict) -> dict:
    outputDir = config["pathInfo"]["outputDir"]
    torsionTopDir = p.join(outputDir, "02_torsion_scanning")
    os.makedirs(torsionTopDir, exist_ok=True)
    config["pathInfo"]["torsionTopDir"] = torsionTopDir

    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def scan_in_serial(scanDfs, torsionScanDir, conformerPdbs, torsionIndexes, batchIndex, config) -> List[pd.DataFrame]:
    argsList = [(conformerPdb, torsionScanDir,  torsionIndexes, config) for  conformerPdb in conformerPdbs]

    for args in argsList:
        forwardsDf, backwardsDf = do_the_twist_worker(args)
        if forwardsDf is not None and backwardsDf is not None:
            scanDfs.append(forwardsDf)
            scanDfs.append(backwardsDf)

    return scanDfs

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def scan_in_parallel(scanDfs, torsionScanDir, conformerPdbs, torsionIndexes, batchIndex, config) -> List[pd.DataFrame]:
    tqdmBarOptions = {
        "desc": f"\033[32mScanning Batch {batchIndex}\033[0m",
        "ascii": "-ϟ",  
        "colour": "yellow",
        "unit":  "scan",
        "dynamic_ncols": True
    }
    argsList = [(conformerPdb, torsionScanDir,  torsionIndexes, config) for  conformerPdb in conformerPdbs]

    with WorkerPool(n_jobs = config["hardwareInfo"]["nCores"]) as pool:
        results = pool.map(do_the_twist_worker,
                            make_single_arguments(argsList),
                              progress_bar=True,
                              iterable_len = len(argsList),
                              progress_bar_options=tqdmBarOptions)
    
    for forwardsDf, backwardsDf in results:
        if forwardsDf is not None and backwardsDf is not None:
            scanDfs.append(forwardsDf)
            scanDfs.append(backwardsDf)

    return scanDfs

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def process_scan_data(scanDfs: List[pd.DataFrame],
                       torsionTopDir: DirectoryPath,
                         batchIndex: int):
    
    ## merge scan dataframes
    scanDfs = [process_energy_outputs(scanDf) for scanDf in scanDfs]
    mergedDf = merge_scan_dfs(scanDfs)

    
    ## remove cols with large jumps in energy
    mergedDf = detect_jumps_in_data(mergedDf)
    ## make a dir to store output csv files
    dataDir = p.join(torsionTopDir, "scan_data")
    os.makedirs(dataDir, exist_ok=True)
    ## write to csv
    batchCsv = p.join(dataDir, f"batch_{batchIndex}.csv")
    mergedDf.to_csv(batchCsv, index=False)
    averagesCsv = p.join(dataDir, "scan_averages.csv")

    if batchIndex == 1:
        ## init the location of the averages csv
        scanAverageDf = pd.DataFrame()
        scanAverageDf["Angle"] = mergedDf["Angle"]
    else:   
        ## read data from averages csv
        scanAverageDf = pd.read_csv(averagesCsv, index_col=None)
    ## calculate averages
    scanAverageDf[f"Batch {batchIndex}"] = mergedDf.drop(columns="Angle").mean(axis=1)
    scanAverageDf.to_csv(averagesCsv, index=False)

    if batchIndex >= 2:
        rollingAveragesCsv = p.join(dataDir, "rolling_averages.csv")
        if p.exists(rollingAveragesCsv):
            rollingAverageDf = pd.read_csv(rollingAveragesCsv, index_col=None)
        else:
            rollingAverageDf = pd.DataFrame()
            rollingAverageDf["Angle"] = mergedDf["Angle"]
        rollingAverageDf[f"Batch {batchIndex}"] = scanAverageDf.drop(columns="Angle").mean(axis=1)
        rollingAverageDf.to_csv(rollingAveragesCsv, index=False)

    if batchIndex >= 3:
        meanAverageError = (rollingAverageDf[f"Batch {batchIndex}"] - rollingAverageDf[f"Batch {batchIndex - 1}"]).abs().mean()
        return meanAverageError, mergedDf, scanAverageDf, rollingAverageDf
    
    else:
        return np.inf, mergedDf, None, None

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def gen_conformers(inputPdb: FilePath,
                    batchDir: DirectoryPath,
                      iteration: int,
                        nConformers: int):
    
    ## make a new dir to store this batch of conformers
    
    conformerDir = p.join(batchDir, "conformers")
    os.makedirs(conformerDir, exist_ok=True)
    ## extract information for naming conformers
    molName = p.basename(inputPdb).split(".")[0]

    firstConformerIndex = (iteration) * nConformers - nConformers + 1

    mol = Chem.MolFromPDBFile(inputPdb, removeHs=False)
    confs = AllChem.EmbedMultipleConfs(mol, numConfs=nConformers, randomSeed=iteration)
    conformerPdbs = []
    for i in range(nConformers):    
        mol.SetProp("_Name", f"Conformer_{i+firstConformerIndex}")
        conformerPdb: FilePath = p.join(conformerDir, f"{molName}_conformer_{i+firstConformerIndex}.pdb")
        MolToPDBFile(mol, conformerPdb, confId=i)
        conformerPdbs.append(conformerPdb)
    return conformerPdbs

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_twist_worker(args):
    ## unpack args tuple
    conformerPdb, torsionScanDir,  torsionIndexes, config= args
    try:
        forwardsDf, backwardsDf = do_the_twist(conformerPdb, torsionScanDir, torsionIndexes, config)
        return forwardsDf, backwardsDf
    except Exception as e:
        # ## delete scan dir
        # conformerId = p.basename(conformerPdb).split(".")[0]
        # conformerScanDir = p.join(torsionScanDir, f"scans_conformer{conformerId}")
        # # rmtree(conformerScanDir)
        return None, None

        
        
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_twist(conformerPdb: FilePath,
                 torsionScanDir: DirectoryPath,
                 torsionIndexes: List[int],
                 config: dict) -> Tuple[pd.DataFrame]:

    conformerId = p.basename(conformerPdb).split(".")[0]

    conformerScanDir = p.join(torsionScanDir, f"scans_conformer{conformerId}")
    optXyz = run_optimisation_step(conformerPdb, torsionIndexes, conformerScanDir, conformerId, config)
    ## get angle of torsion in this conformer
    initalTorsionAngle = measure_current_torsion_angle(optXyz, torsionIndexes)

    forwardsScanDf, forwardsXyz = run_forwards_scan_step(optXyz, initalTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

    backwardsScanDf = run_backwards_scan_step(forwardsXyz, initalTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

    return forwardsScanDf, backwardsScanDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_optimisation_step(conformerPdb, torsionIndexes, conformerScanDir, conformerId, config):
            ## XTB-GFN2 OPTIMISATION ##
    optDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_optimise")
    os.makedirs(optDir, exist_ok=True)
    ## make an ORCA input file for optimisation
    optOrcaInput: FilePath = generate_orca_input(inputGeom=conformerPdb,
                                                  torsionIndexes=torsionIndexes,
                                                   outDir = optDir,
                                                    charge = config["moleculeInfo"]["charge"],
                                                      multiplicity = config["moleculeInfo"]["multiplicity"],
                                                        scanAngles = None)
    ## Run orca optimisation
    optOrcaOuput: FilePath = p.join(optDir, "orca.out")
    if not p.isfile(optOrcaOuput):
        run_orca(optOrcaInput, optOrcaOuput)
    optXyz = p.join(optDir, "orca.xyz")
    return optXyz

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_forwards_scan_step(optXyz, initalTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config):
    ## XTB-GFN2 FORWARDS SCAN ##
    ## do a forwards scan
    forwardsDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_forwards")
    os.makedirs(forwardsDir, exist_ok=True)
    forwardsScanAngle = initalTorsionAngle + 360
    forwardsScanText = f"{str(initalTorsionAngle)}, {str(forwardsScanAngle)}, {str(config['torsionScanInfo']['nScanSteps'])}"


    forwardsOrcaInput: FilePath = generate_orca_input(inputGeom=optXyz,
                                                      torsionIndexes=torsionIndexes,
                                                      outDir = forwardsDir,
                                                        charge = config["moleculeInfo"]["charge"],
                                                        multiplicity = config["moleculeInfo"]["multiplicity"],
                                                        scanAngles = forwardsScanText)
    forwardsOrcaOutput: FilePath = p.join(forwardsDir, "orca.out")
    if not p.isfile(forwardsOrcaOutput):
        run_orca(forwardsOrcaInput, forwardsOrcaOutput)
    forwardsXyz = find_final_xyz(forwardsDir)
    forwardsScanDf = read_data(forwardsDir)
    forwardsScanDf.columns = ["Angle", f"{conformerId}_f"]
    forwardsScanDf["Angle"] = forwardsScanDf["Angle"].apply(rescale_torsion_angles)
    forwardsScanDf = forwardsScanDf.sort_values(by="Angle", ascending=True)
    forwardsScanDf = average_duplicate_angles(forwardsScanDf)

    return forwardsScanDf, forwardsXyz
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_backwards_scan_step(optXyz, initalTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config):
    ## do a backwards scan
    backwardsDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_backwards")
    os.makedirs(backwardsDir, exist_ok=True)
    forwardsScanAngle = initalTorsionAngle + 360 

    backwardsScanText = f"{str(forwardsScanAngle)}, {str(initalTorsionAngle)}, {str(config['torsionScanInfo']['nScanSteps'])}"

    backwardsOrcaInput: FilePath = generate_orca_input(inputGeom = optXyz,
                                                       torsionIndexes = torsionIndexes,
                                                       outDir =backwardsDir,                                      
                                                        charge = config["moleculeInfo"]["charge"],
                                                        multiplicity = config["moleculeInfo"]["multiplicity"],
                                                        scanAngles = backwardsScanText)
    backwardsOrcaOutput: FilePath = p.join(backwardsDir, "orca.out")
    if not p.isfile(backwardsOrcaOutput):
        run_orca(backwardsOrcaInput, backwardsOrcaOutput)
    backwardsScanDf = read_data(backwardsDir)

    backwardsScanDf.columns = ["Angle", f"{conformerId}_b"]
    backwardsScanDf["Angle"] = backwardsScanDf["Angle"].apply(rescale_torsion_angles)
    backwardsScanDf = backwardsScanDf.sort_values(by="Angle", ascending=True)
    backwardsScanDf = average_duplicate_angles(backwardsScanDf)
    return backwardsScanDf  

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def average_duplicate_angles(df):
    return  df.groupby('Angle', as_index=False).mean()

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def merge_scan_dfs(scanDfs):
    mergedDf = pd.DataFrame()
    mergedDf["Angle"] = scanDfs[0]["Angle"]
    for scanDf in scanDfs:
        if scanDf.columns[1] not in mergedDf.columns:
            mergedDf[scanDf.columns[1]] = scanDf[scanDf.columns[1]]
    return mergedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def process_energy_outputs(scanDf):
    minEnergy = scanDf.iloc[:,1].min()
    scanDf.iloc[:,1] = scanDf.iloc[:,1] - minEnergy
    scanDf.iloc[:,1] = scanDf.iloc[:,1].apply(hartree_to_kcal_per_mol)
    return scanDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def hartree_to_kcal_per_mol(energy):
    return energy * 627.5095

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rescale_torsion_angles(angle):
    angle = angle % 360  # First, reduce the angle to the 0-360 range
    if angle > 180:
        angle -= 360  # Shift angles greater than 180 to the negative side
    return angle

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def find_final_xyz(conformertorsionScanDir):
    allXyzFiles = glob.glob(p.join(conformertorsionScanDir, "*.xyz"))
    scanXYZFiles = sorted([f for f in allXyzFiles if re.match(r'orca\.\d+\.xyz$', os.path.basename(f))])

    finalXyzFile = scanXYZFiles[-1]
    return finalXyzFile


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def measure_current_torsion_angle(conformerXyz, torsionIndexes):
    atomCoords = xyz_to_np_array(conformerXyz)
    torsionAngle = calculate_torsion_angle(atomCoords, torsionIndexes)
    roundedToTenAngle = round(torsionAngle, -1)
    return roundedToTenAngle

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_torsion_angle(coords, torsionIndexes):
    # Vectors between points
    b1 = coords[torsionIndexes[1]] - coords[torsionIndexes[0]]
    b2 = coords[torsionIndexes[2]] - coords[torsionIndexes[1]]
    b3 = coords[torsionIndexes[3]] - coords[torsionIndexes[2]]

    # Normal vectors to the planes
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    # Normalize the normal vectors
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)

    # Calculate the angle
    angleInradians = np.arctan2(np.dot(np.cross(n1, n2), b2 / np.linalg.norm(b2)), np.dot(n1, n2))

    # Convert from radians to degrees
    angleIdDegrees = np.degrees(angleInradians)

    return angleIdDegrees
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_torsion(df, torsionTopDir):
    plt.figure(figsize=(12, 8))
    
    # Plot each column as a scatter plot except 'Angle' and 'Average'
    for col in df.columns:
        if col not in ['Angle', 'Average']:
            try:
                plt.scatter(df['Angle'], df[col], label=col, cmap='plasma', alpha=0.7)
            except:
                pass
    # Plot 'Average' as a black line
    plt.plot(df['Angle'], df['Average'], color='black', linewidth=2, label='Average')
    
    # Add labels and legend
    plt.xlabel('Angle')
    plt.ylabel('Values')
    plt.title('Data Plot')
    plt.grid(True)

    # Save the plot
    plt.savefig(p.join(torsionTopDir, "torsion.png"))
    plt.close()
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def identify_rotatable_bonds(inputPdb) -> List[Tuple[int,int,int,int]]:
    # Load the molecule from a PDB file
    mol = Chem.MolFromPDBFile(inputPdb, removeHs=False)
    # Identify torsion angles for rotatable bonds
    torsionAngles = []
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        if bond.GetBondType() == Chem.BondType.SINGLE:
            atom2 = bond.GetBeginAtom()
            atom3 = bond.GetEndAtom()

            ## get atom names for begin and end atoms
            atom2Name = atom2.GetPDBResidueInfo().GetName().strip()
            atom3Name = atom3.GetPDBResidueInfo().GetName().strip()

            ## dont scan amide bonds
            if atom2Name in ["NN", "C"] and atom3Name in ["NN", "C"]:
                continue
            if atom2Name in ["N", "CC1"] and atom3Name in ["N", "CC1"]:
                continue

            
            if not (atom2.IsInRing() or atom3.IsInRing()):
                # Find neighboring atoms for torsion angle
                neighborsBegin = [a for a in atom2.GetNeighbors() if a.GetIdx() != atom3.GetIdx()]
                neighborsEnd = [a for a in atom3.GetNeighbors() if a.GetIdx() != atom2.GetIdx()]
                
                if neighborsBegin and neighborsEnd:
                    atom1 = neighborsBegin[0]
                    atom4 = neighborsEnd[0]
                    
                    atom1Name = atom1.GetPDBResidueInfo().GetName().strip()
                    atom4Name = atom4.GetPDBResidueInfo().GetName().strip()

                    ## dont scan bonds with non-polar hydrogens at either as atoms 1 or 4
                    if atom1Name.startswith("H"):
                        if atom2Name.startswith("C"):
                            continue
                    if atom4Name.startswith("H"):
                        if atom3Name.startswith("C"):
                            continue
                    ## add torsion data to list
                    torsionAngles.append({
                        'atoms': (atom1Name, atom2Name, atom3Name, atom4Name),
                        'indices': (atom1.GetIdx(), atom2.GetIdx(), atom3.GetIdx(), atom4.GetIdx())
                    })

    return torsionAngles

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def read_data(conformertorsionScanDir):
    scanDat: FilePath = p.join(conformertorsionScanDir, "orca.relaxscanact.dat")
    if not os.path.exists(scanDat):
        raise FileNotFoundError(f"Scan data not found in {conformertorsionScanDir}")
    scanDf: pd.DataFrame = pd.read_csv(scanDat, sep='\s+', header=None)
    return scanDf
    
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_orca(orcaInput, orcaOutput):
    orcaCommand = ["orca", orcaInput]
    with open(orcaOutput, 'w') as output_file:
        try:
            call(orcaCommand, stdout=output_file, stderr=output_file)
        except Exception as e:
            pass
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def generate_orca_input(inputGeom, torsionIndexes, outDir, charge, multiplicity, scanAngles):

    ## create orca input file
    orcaInputFile = p.join(outDir, "orca.inp")
    with open(orcaInputFile, "w") as f:
        ## METHOD
        f.write("! XTB2 Opt\n")
        ## SCAN (OPTIONAL)
        if not scanAngles is None:
            f.write("%geom Scan\n")
            torsionText: str = f"D {' '.join(map(str, torsionIndexes))} = {scanAngles}\n"
            f.write(torsionText)
            f.write("end\n")
            f.write("end\n")

        ## GEOMETRY
        ## if we have a pdb file, read coords and place in orca input file
        if p.splitext(inputGeom)[1] == ".pdb":
            conformerDf = pdbUtils.pdb2df(inputGeom)
            f.write(f"*xyz {charge} {multiplicity}\n\n")
            for index, row in conformerDf.iterrows():
                geomLine = f"{row['ELEMENT']} {row['X']} {row['Y']} {row['Z']}\n"
                f.write(geomLine)
            f.write("*")
        ## if we have an xyz file, put path in input line instead
        elif p.splitext(inputGeom)[1] == ".xyz":
            f.write(f"*xyzfile {charge} {multiplicity} {inputGeom}\n\n")

    return orcaInputFile
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def xyz_to_np_array(filePath):
    with open(filePath, 'r') as file:
        lines = file.readlines()

    # Skip the first two lines (atom count and comment)
    dataLines = lines[2:]

    # Parse the coordinates
    coords = []
    for line in dataLines:
        parts = line.split()
        if len(parts) >= 4:
            x, y, z = map(float, parts[1:4])
            coords.append([x, y, z])

    return np.array(coords)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def detect_jumps_in_data(df):
    # Calculate differences between consecutive values
    diffDf = df.drop(columns='Angle').diff().abs()
    
    # Identify columns with any difference greater than the threshold
    jumpyCols = diffDf.columns[((diffDf > 3).any())]
    
    # Drop these columns from the original DataFrame
    cleanDf = df.drop(columns=jumpyCols)
    
    return cleanDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_means(df):
    return df.drop(columns='Angle').mean(axis=1)

if __name__ == "__main__":
    twist_protocol()