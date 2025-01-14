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

## PLOTTING LIBRARIES ##
import matplotlib.pyplot as plt
import matplotlib.cm as cm

## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

from typing import List, Tuple, Dict
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def dummy_inputs() -> dict:
    config: dict = {
        "pathInfo": {
            "inputDir": "/home/esp/scriptDevelopment/drFrankenstein/dummy_input",
            "outputDir": "/home/esp/scriptDevelopment/drFrankenstein/outputs"
        },
        "moleculeInfo": {
            "charge": 0,
            "multiplicity": 1, 
            "moleculePdb": "/home/esp/scriptDevelopment/drFrankenstein/dummy_input/1-hexene.pdb"
        }, 
        "scanInfo": {
            "nScanSteps": 37,
            "convergenceTolerance": 0.1,
            "minBatches": 5,
            "maxBatches": 30
        },
        "hardwareInfo": {
            "nCores": 16
        }
    }
    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def twist_protocol(config):
    ## UNPACK CONFIG DICT ##
    ## PATH INFO
    inputDir: FilePath = config["pathInfo"]["inputDir"]
    outputDir: FilePath = config["pathInfo"]["outputDir"]
    ## MOLECULE INFO
    moleculeInfo: dict = config["moleculeInfo"]
    charge: int = moleculeInfo["charge"]
    multiplicity: int = moleculeInfo["multiplicity"]
    moleculePdb: FilePath = moleculeInfo["moleculePdb"]
    ## SCAN INFO
    scanInfo: dict = config["scanInfo"]
    nScanSteps: int = scanInfo["nScanSteps"]
    convergenceTolerance: float = scanInfo["convergenceTolerance"]
    minBatches: int = scanInfo["minBatches"]
    maxBatches: int = scanInfo["maxBatches"]
    ## HARDWARE INFO
    hardwareInfo: dict = config["hardwareInfo"]
    nCores: int = hardwareInfo["nCores"]

    ## make top level output dir
    os.makedirs(outputDir,exist_ok=True)

    ## identify rotatable bonds to be scanned
    rotatableBonds: Dict[List[Tuple[int,int,int,int]]] = identify_rotatable_bonds(moleculePdb, outputDir)

    ## loop over torsions that need to be scanned
    for rotatableBond in rotatableBonds:
        print(f"Running torsion scans for {rotatableBond['atoms']}")
        ## make a top level dir for this torsion
        torsionId: str = "-".join(map(str, rotatableBond["atoms"])) 
        torsionTopDir = p.join(outputDir, f"torsion_{torsionId}" )
        os.makedirs(torsionTopDir, exist_ok=True)

        ## init some variables
        meanAverageError = np.inf
        batchIndex = 1
        meanAverageErrors = []
        scanRawDfs = {}
        
        ## get torsion index
        torsionIndexes = rotatableBond["indices"]

        ## INIT A WHILE LOOP - BREAK WHEN CONVERGENCE IS REACHED
        while (meanAverageError > convergenceTolerance or batchIndex < minBatches) or (batchIndex + 1 == maxBatches):
            batchDir = p.join(torsionTopDir, f"batch_{batchIndex}")
            os.makedirs(batchDir, exist_ok=True)
            ## generate a batch of conformers
            conformerPdbs = gen_conformers(inputPdb=moleculePdb,
                                            batchDir=batchDir,
                                            iteration=batchIndex,
                                            nConformers=nCores)
            ## CREATE A DIR FOR THIS TORSION ANGLE
            torsionDir = p.join(batchDir, f"torsions")
            os.makedirs(torsionDir, exist_ok=True)
            ## INIT EMPTY LIST TO STORE SCAN DATA
            scanDfs = []
            ## CONSTRUCT ARGS LIST FOR MULTIPROCESSING
            argsList = [
                (conformerPdb, torsionDir, torsionIndexes, charge, multiplicity, nScanSteps)
                for  conformerPdb in conformerPdbs
            ]
            ## run torsion scans in paralell
            scanDfs = scan_in_paralell(scanDfs, argsList, nCores, batchIndex)
   

            meanAverageError, mergedDf, scanAverageDf, rollingAverageDf = process_scan_data(scanDfs, torsionTopDir, batchIndex)
            meanAverageErrors.append(meanAverageError)
            scanRawDfs[batchIndex] = mergedDf

            ## STEP BATCH INDEX
            batchIndex += 1
        else:
            finalScanEnergiesDf = pd.DataFrame()
            finalScanEnergiesDf["Angle"] = scanAverageDf["Angle"]
            finalScanEnergiesDf[torsionId] = rollingAverageDf.iloc[:, -1]

            dataDir = p.join(torsionTopDir, "scan_data")
            finalScanEnergiesDf.to_csv(p.join(dataDir, "final_scan_energies.csv"), index=False)
            



        plot_torsion_scans(torsionTopDir, scanRawDfs, scanAverageDf, rollingAverageDf, meanAverageErrors)
        exit()

def scan_in_paralell(scanDfs, argsList, nCores, batchIndex) -> List[pd.DataFrame]:


    tqdmBarOptions = {
        "desc": f"Scanning Batch {batchIndex}",
        "ascii": "-🗲",  # Use the character for the bar
        "colour": "yellow",
        "unit":  "scan",
        "dynamic_ncols": True
    }

    with WorkerPool(n_jobs=nCores) as pool:
        results = pool.map(process_conformer,
                            make_single_arguments(argsList)
                            , progress_bar=True,
                              iterable_len = len(argsList),
                              progress_bar_options=tqdmBarOptions)
    
    for forwardsDf, backwardsDf in results:
        if forwardsDf is not None and backwardsDf is not None:
            scanDfs.append(forwardsDf)
            scanDfs.append(backwardsDf)

    return scanDfs


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_torsion_scans(torsionTopDir, scanRawDfs, scanAverageDf, rollingAverageDf, meanAverageErrors):
    plotDir = p.join(torsionTopDir, "plots")
    os.makedirs(plotDir, exist_ok=True)

    plot_raw_data(scanRawDfs, plotDir)
    plot_average_data(scanAverageDf, plotDir)
    plot_rolling_average_data(rollingAverageDf, plotDir)
    plot_mean_average_data(meanAverageErrors, plotDir)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_raw_data(scanRawDfs, plotDir):
    colors = cm.YlGn(np.linspace(0, 1, len(scanRawDfs)))
    
    plt.figure(figsize=(12, 8))
    
    for (batchIndex, df), color in zip(scanRawDfs.items(), colors):
        for column in df.columns:
            if column != 'Angle':
                plt.scatter(df['Angle'], df[column], label=f'{batchIndex} - {column}', color=color)
    

    plt.xlabel('Angle')
    plt.ylabel('Relative Energy (kcal/mol)')
    plt.title('Raw Data')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(p.join(plotDir, "raw_data.png"))
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_average_data(scanAverageDf, plotDir):
    colors = cm.YlGn(np.linspace(0, 1, len(scanAverageDf.columns)-1))

    plt.figure(figsize=(12, 8))
    
    for column, color in zip(scanAverageDf.columns, colors):
        if column != 'Angle':
            plt.plot(scanAverageDf['Angle'], scanAverageDf[column], label=column, color=color)
    
    plt.xlabel('Angle')
    plt.ylabel('Relative Energy (kcal/mol)')
    plt.title('Average Data')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(p.join(plotDir, "average_data.png"))

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_rolling_average_data(rollingAverageDf, plotDir):
    colors = cm.YlGn(np.linspace(0, 1, len(rollingAverageDf.columns)-1))

    plt.figure(figsize=(12, 8))
    
    for column, color in zip(rollingAverageDf.columns, colors):
        if column != 'Angle':
            plt.plot(rollingAverageDf['Angle'], rollingAverageDf[column], label=column, color=color)
    
    plt.xlabel('Angle')
    plt.ylabel('Relative Energy (kcal/mol)')
    plt.title('Rolling Average Data')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(p.join(plotDir, "rolling_average_data.png"))
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_mean_average_data(meanAverageErrors, plotDir):
    plt.figure(figsize=(12, 8))
    
    plt.plot(meanAverageErrors)
    plt.ylabel('Relative Energy (kcal/mol)')
    plt.title('Mean Average Error')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(p.join(plotDir, "mean_average_error.png"))
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
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def process_conformer(args):
    conformerPdb, torsionDir, torsionIndex, charge, multiplicity, nScanSteps = args
    try:
        forwardsDf, backwardsDf = do_the_twist(
            torsionDir,
            conformerPdb,
            torsionIndex,
            charge,
            multiplicity,
            nScanSteps
        )
        return forwardsDf, backwardsDf
    except Exception as e:
        ## delete scan dir
        conformerId = p.basename(conformerPdb).split(".")[0]
        conformerScanDir = p.join(torsionDir, f"scans_conformer{conformerId}")
        # rmtree(conformerScanDir)
        return None, None

        
        
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_twist(torsionScanDir: DirectoryPath,
                 conformerPdb: FilePath,
                 torsionIndex: List[int],
                 charge: int,
                 multiplicity: int,
                 nScanSteps: int) -> Tuple[pd.DataFrame]:

    conformerId = p.basename(conformerPdb).split(".")[0]

    conformerScanDir = p.join(torsionScanDir, f"scans_conformer{conformerId}")
    ## do an initial optimisation with XTB
    optDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_optimise")
    os.makedirs(optDir, exist_ok=True)
    # extract xyz and element from input pdb
    optOrcaInput: FilePath = generate_orca_input(conformerPdb, torsionIndex, optDir, charge, multiplicity, scanAngles = None)
    optOrcaOuput: FilePath = p.join(optDir, "orca.out")
    if not p.isfile(optOrcaOuput):
        run_orca(optOrcaInput, optOrcaOuput)
    optXyz = p.join(optDir, "orca.xyz")

    initalTorsionAngle = measure_current_torsion_angle(optXyz, torsionIndex)

    ## do a forwards scan
    forwardsDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_forwards")
    os.makedirs(forwardsDir, exist_ok=True)
    forwardsScanAngle = initalTorsionAngle + 360
    forwardsScanText = f"{str(initalTorsionAngle)}, {str(forwardsScanAngle)}, {str(nScanSteps)}"

    forwardsOrcaInput: FilePath = generate_orca_input(optXyz, torsionIndex, forwardsDir, charge, multiplicity, scanAngles = forwardsScanText)
    forwardsOrcaOutput: FilePath = p.join(forwardsDir, "orca.out")
    if not p.isfile(forwardsOrcaOutput):
        run_orca(forwardsOrcaInput, forwardsOrcaOutput)
    forwardsXyz = find_final_xyz(forwardsDir)
    forwardsScanDf = read_data(forwardsDir)
    forwardsScanDf.columns = ["Angle", f"{conformerId}_f"]
    forwardsScanDf["Angle"] = forwardsScanDf["Angle"].apply(rescale_torsion_angles)
    forwardsScanDf = forwardsScanDf.sort_values(by="Angle", ascending=True)
    forwardsScanDf = average_duplicate_angles(forwardsScanDf)

    ## do a backwards scan
    backwardsDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_backwards")
    os.makedirs(backwardsDir, exist_ok=True)
    scanAngles_backwards = f"{str(forwardsScanAngle)}, {str(initalTorsionAngle)}, {str(nScanSteps)}"

    backwardsOrcaInput: FilePath = generate_orca_input(forwardsXyz, torsionIndex, backwardsDir, charge, multiplicity, scanAngles = scanAngles_backwards)
    backwardsOrcaOutput: FilePath = p.join(backwardsDir, "orca.out")
    if not p.isfile(backwardsOrcaOutput):
        run_orca(backwardsOrcaInput, backwardsOrcaOutput)
    backwardsScanDf = read_data(backwardsDir)

    backwardsScanDf.columns = ["Angle", f"{conformerId}_b"]
    backwardsScanDf["Angle"] = backwardsScanDf["Angle"].apply(rescale_torsion_angles)
    backwardsScanDf = backwardsScanDf.sort_values(by="Angle", ascending=True)
    backwardsScanDf = average_duplicate_angles(backwardsScanDf)
    return forwardsScanDf, backwardsScanDf  

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
def find_final_xyz(conformerTorsionDir):
    allXyzFiles = glob.glob(p.join(conformerTorsionDir, "*.xyz"))
    scanXYZFiles = sorted([f for f in allXyzFiles if re.match(r'orca\.\d+\.xyz$', os.path.basename(f))])

    finalXyzFile = scanXYZFiles[-1]
    return finalXyzFile


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def measure_current_torsion_angle(conformerXyz, torsionIndex):
    atomCoords = xyz_to_np_array(conformerXyz)
    torsionAngle = calculate_torsion_angle(atomCoords)
    roundedToTenAngle = round(torsionAngle, -1)
    return roundedToTenAngle

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_torsion_angle(coords):
    # Vectors between points
    b1 = coords[1] - coords[0]
    b2 = coords[2] - coords[1]
    b3 = coords[3] - coords[2]

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
def identify_rotatable_bonds(inputPdb, outputDir) -> List[Tuple[int,int,int,int]]:
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
            if atom2Name in ["N", "CC1"] and atom3Name in ["NN", "CC1"]:
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
                    
                    torsionAngles.append({
                        'atoms': (atom1Name, atom2Name, atom3Name, atom4Name),
                        'indices': (atom1.GetIdx(), atom2.GetIdx(), atom3.GetIdx(), atom4.GetIdx())
                    })

    

    return torsionAngles

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def read_data(conformerTorsionDir):
    scanDat: FilePath = p.join(conformerTorsionDir, "orca.relaxscanact.dat")
    if not os.path.exists(scanDat):
        raise FileNotFoundError(f"Scan data not found in {conformerTorsionDir}")
    scanDf: pd.DataFrame = pd.read_csv(scanDat, delim_whitespace=True, header=None)
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
def generate_orca_input(conformerGeom, torsionIndex, conformerTorsionDir, charge, multiplicity, scanAngles):

    ## create orca input file
    orcaInputFile = p.join(conformerTorsionDir, "orca.inp")
    with open(orcaInputFile, "w") as f:
        ## METHOD
        f.write("! XTB2 Opt\n")
        ## SCAN (OPTIONAL)
        if not scanAngles is None:
            f.write("%geom Scan\n")
            torsionText: str = f"D {' '.join(map(str, torsionIndex))} = {scanAngles}\n"
            f.write(torsionText)
            f.write("end\n")
            f.write("end\n")

        ## GEOMETRY
        ## if we have a pdb file, read coords and place in orca input file
        if p.splitext(conformerGeom)[1] == ".pdb":
            conformerDf = pdbUtils.pdb2df(conformerGeom)
            f.write(f"*xyz {charge} {multiplicity}\n\n")
            for index, row in conformerDf.iterrows():
                geomLine = f"{row['ELEMENT']} {row['X']} {row['Y']} {row['Z']}\n"
                f.write(geomLine)
            f.write("*")
        ## if we have an xyz file, put path in input line instead
        elif p.splitext(conformerGeom)[1] == ".xyz":
            f.write(f"*xyzfile {charge} {multiplicity} {conformerGeom}\n\n")

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

if __name__ == "__main__":
    main()