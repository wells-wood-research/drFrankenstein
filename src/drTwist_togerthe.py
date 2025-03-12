## BASIC IMPORTS ##
import os
from os import path as p
from pdbUtils import pdbUtils
from subprocess import call, PIPE
import pandas as pd
import yaml
import matplotlib
matplotlib.use('Agg')  # Set the Agg backend before importing pyplot
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import glob 
import re
from shutil import rmtree
from scipy.signal import argrelextrema

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
import drOrca
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

from typing import List, Tuple, Dict

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def twist_protocol(config):

#     config = set_up_directories(config)

#     ## identify rotatable bonds to be scanned
#     rotatableBonds: Dict[List[Tuple[int,int,int,int]]] = identify_rotatable_bonds(config["moleculeInfo"]["cappedPdb"])

#     ## loop over torsions that need to be scanned
#     config["pathInfo"]["torsionDirs"] = []
#     config["torsionScanInfo"]["torsionTags"] = []
#     config["torsionScanInfo"]["finalScanEnergies"] = {}
#     for rotatableBond in rotatableBonds:
#         config = run_torsion_scanning(rotatableBond, config)

#     config["checkpointInfo"]["scanningComplete"] = True
#     return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def run_torsion_scanning(rotatableBond, config, debug=False) -> dict:
#     ## make a top level dir for this torsion
#     torsionTag: str = "-".join(map(str, rotatableBond["atoms"])) 
#     torsionDir = p.join(config["pathInfo"]["torsionTopDir"], f"torsion_{torsionTag}" )
#     os.makedirs(torsionDir, exist_ok=True)
#     ## add to config
#     config["pathInfo"]["torsionDirs"].append(torsionDir)
#     config["torsionScanInfo"]["torsionTags"].append(torsionTag)

 
#     ## get conformer XYZ files
#     conformerXyzs = config["pathInfo"]["conformerXyzs"]
#     ## INIT EMPTY LIST TO STORE SCAN DATA
#     scanDfs = []
#     if debug:
#         ## run in serial
#         scanDfs = scan_in_serial(scanDfs, torsionDir, conformerXyzs, rotatableBond["indices"], config)
#     else:
#         # run torsion scans in paralell
#         scanDfs = scan_in_parallel(scanDfs, torsionDir, conformerXyzs, rotatableBond["indices"], torsionTag, config)

#     ## Merge scan data, calculate averages, rolling averages and mean average errors
#     finalScanEnergiesCsv  = process_scan_data(scanDfs, torsionDir, torsionTag)

#     config["torsionScanInfo"]["finalScanEnergies"][torsionTag] = finalScanEnergiesCsv   
#     return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def set_up_directories(config: dict) -> dict:
#     outputDir = config["pathInfo"]["outputDir"]
#     torsionTopDir = p.join(outputDir, "03_torsion_scanning")
#     os.makedirs(torsionTopDir, exist_ok=True)
#     config["pathInfo"]["torsionTopDir"] = torsionTopDir

#     return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def scan_in_serial(scanDfs, torsionScanDir, conformerXyzs, torsionIndexes, config) -> List[pd.DataFrame]:
#     argsList = [(conformerXyz, torsionScanDir,  torsionIndexes, config) for  conformerXyz in conformerXyzs]

#     for args in argsList:
#         forwardsDf, backwardsDf = do_the_twist_worker(args)
#         if forwardsDf is not None and backwardsDf is not None:
#             scanDfs.append(forwardsDf)
#             scanDfs.append(backwardsDf)

#     return scanDfs

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def scan_in_parallel(scanDfs, torsionScanDir, conformerXyzs, torsionIndexes, torsionTag, config) -> List[pd.DataFrame]:
#     tqdmBarOptions = {
#         "desc": f"\033[32mScanning Torsion {torsionTag}\033[0m",
#         "ascii": "-ϟ",  
#         "colour": "yellow",
#         "unit":  "scan",
#         "dynamic_ncols": True
#     }
#     argsList = [(conformerXyz, torsionScanDir,  torsionIndexes, config) for  conformerXyz in conformerXyzs]

#     ## save on cores 
#     nCores = min(len(argsList), config["hardwareInfo"]["nCores"])

#     with WorkerPool(n_jobs = nCores) as pool:
#         results = pool.map(do_the_twist_worker,
#                             make_single_arguments(argsList),
#                               progress_bar=True,
#                               iterable_len = len(argsList),
#                               progress_bar_options=tqdmBarOptions)
    
#     for forwardsDf, backwardsDf in results:
#         if forwardsDf is not None and backwardsDf is not None:
#             scanDfs.append(forwardsDf)
#             scanDfs.append(backwardsDf)

#     return scanDfs

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def process_scan_data(scanDfs: List[pd.DataFrame],
#                        torsionTopDir: DirectoryPath,
#                        torsionTag: str):

#     ## make a dir to store output csv files
#     dataDir = p.join(torsionTopDir, "scan_data")
#     os.makedirs(dataDir, exist_ok=True)

#     ## merge scan dataframes
#     # scanDfs = [process_energy_outputs(scanDf) for scanDf in scanDfs]
#     mergedDf = merge_scan_dfs(scanDfs)
 
#     ## remove cols with large jumps in energy
#     mergedDf = detect_jumps_in_data(mergedDf)

#     ## write to csv
#     mergedScanCsv = p.join(dataDir, f"scan_energies.csv")
#     mergedDf.to_csv(mergedScanCsv, index=False)

#     scanAverageDf = pd.DataFrame()
#     scanAverageDf["Angle"] = mergedDf["Angle"]

#     ## calculate averages
#     finalScanEnergiesCsv = p.join(dataDir, "final_scan_energies.csv")
#     scanAverageDf[torsionTag] = mergedDf.drop(columns="Angle").mean(axis=1)
#     scanAverageDf.to_csv(finalScanEnergiesCsv, index=False)

#     return finalScanEnergiesCsv


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def do_the_twist_worker(args):
#     ## unpack args tuple
#     conformerXyz, torsionScanDir,  torsionIndexes, config= args
#     try:
#         forwardsDf, backwardsDf = do_the_twist(conformerXyz, torsionScanDir, torsionIndexes, config)
#         return forwardsDf, backwardsDf
#     except Exception as e:
#         raise(e)
#         # ## delete scan dir
#         # conformerId = p.basename(conformerPdb).split(".")[0]
#         # conformerScanDir = p.join(torsionScanDir, f"scans_conformer{conformerId}")
#         # # rmtree(conformerScanDir)
#         return None, None

        
        
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def do_the_twist(conformerXyz: FilePath,
#                  torsionScanDir: DirectoryPath,
#                  torsionIndexes: List[int],
#                  config: dict) -> Tuple[pd.DataFrame]:

#     conformerId = p.basename(conformerXyz).split(".")[0]

#     conformerScanDir = p.join(torsionScanDir, f"scans_{conformerId}")
#     optXyz = run_optimisation_step(conformerXyz, conformerScanDir, conformerId, config)
#     ## get angle of torsion in this conformer
#     initialTorsionAngle = measure_current_torsion_angle(optXyz, torsionIndexes)

#     forwardsScanDf, forwardsXyz, forwardsDir = run_forwards_scan_step(optXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

#     backwardsScanDf, backwardsDir = run_backwards_scan_step(forwardsXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

#     ## if no QM method has been specified, return data as-is
#     if config["torsionScanInfo"]["singlePointMethod"] is None:
#         forwardsScanDf = process_energy_outputs(forwardsScanDf)
#         backwardsScanDf = process_energy_outputs(backwardsScanDf)
#         return forwardsScanDf, backwardsScanDf
    
#     ## otherwise apply the single-point protocol
#     forwardsSinglePointDf = run_singlepoints_on_scans(scanDir=forwardsDir,
#                               scanDf = forwardsScanDf, 
#                               outDir=conformerScanDir,
#                               conformerId=conformerId,
#                               config = config,
#                               tag = "forwards")
    

#     backwardsSinglePointDf = run_singlepoints_on_scans(scanDir=backwardsDir,
#                               scanDf = backwardsScanDf,
#                               outDir=conformerScanDir,
#                               conformerId=conformerId,
#                               config = config,
#                               tag = "backwards")
    
#     forwardsScanDf = process_energy_outputs(forwardsScanDf)
#     backwardsScanDf = process_energy_outputs(backwardsScanDf)

#     forwardsSinglePointDf = process_energy_outputs(forwardsSinglePointDf)
#     backwardsSinglePointDf = process_energy_outputs(backwardsSinglePointDf)



#     drPlotter.plot_scan_singlepoint_comparison(forwardsScanDf, forwardsSinglePointDf, conformerScanDir, tag="forwards")
#     drPlotter.plot_scan_singlepoint_comparison(backwardsScanDf, backwardsSinglePointDf, conformerScanDir, tag="backwards")

#     return forwardsSinglePointDf, backwardsSinglePointDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

# def rescale_and_sort_energy_angles(energyDf):
#     energyDf["Energy"] = energyDf["Energy"] - energyDf["Energy"].min()
#     energyDf["Angle"] = energyDf["Angle"].apply(rescale_torsion_angles)
#     energyDf = energyDf.sort_values(by="Angle", ascending=True)
#     return energyDf


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def run_singlepoints_on_scans(scanDir, scanDf, outDir, conformerId,  config, tag):

#     scanXyzs = find_scan_xyz_files(scanDir, expectedNumberOfFiles=config["torsionScanInfo"]["nScanSteps"])

#     singlePointsOn = config["torsionScanInfo"]["scanSinglePointsOn"]
#     if singlePointsOn == "all":
#         stationaryPointScanIndexes = scanDf["scan_index"].to_list()

#     elif singlePointsOn == "minMaxOnly":
#         stationaryPointsIndexes = find_local_extrema(scanDf["Energy"])
#         stationaryPointScanIndexes = scanDf.loc[stationaryPointsIndexes, "scan_index"].to_list()
#         scanXyzs = [scanXyz for scanXyz in scanXyzs if scanXyz.split(".")[1] in stationaryPointScanIndexes]

#     elif singlePointsOn == "minMaxMiddle":
#         stationaryPointsIndexes = find_local_extrema(scanDf["Energy"])
#         stationaryAndMidPointIndexes = add_mid_points(stationaryPointsIndexes)
#         stationaryPointScanIndexes = scanDf.loc[stationaryAndMidPointIndexes, "scan_index"].to_list()
#         scanXyzs = [scanXyz for scanXyz in scanXyzs if scanXyz.split(".")[1] in stationaryPointScanIndexes]

#     singlePointEnergies = {}
#     for scanXyz in scanXyzs:
#         scanId = scanXyz.split(".")[1]
#         scanDir = p.join(outDir, f"SP_{conformerId}_{tag}_{scanId}")
#         os.makedirs(scanDir, exist_ok=True)
#         spOrcaInput: FilePath = drOrca.make_orca_input_for_singlepoint(inputXyz=scanXyz,
#                                                                        outDir= scanDir,
#                                                                        moleculeInfo=config["moleculeInfo"],
#                                                                        qmMethod=config["torsionScanInfo"]["singlePointMethod"],
#                                                                        solvationMethod=config["torsionScanInfo"]["singlePointSolvationMethod"])
#         spOrcaOutput : FilePath = p.join(scanDir, "orca_sp.out")
#         if not p.isfile(spOrcaOutput):
#             run_orca(spOrcaInput, spOrcaOutput)
#         singlePointEnergy = read_singlepoint_energy(spOrcaOutput)
#         singlePointEnergies[scanId] = singlePointEnergy
#     singlePointEnergyDf = pd.DataFrame(list(singlePointEnergies.items()), columns = ["scan_index", "Energy"])

#     singlePointEnergyDf = singlePointEnergyDf.merge(
#         scanDf[["scan_index", "Angle"]], on="scan_index", how="left")

#     return singlePointEnergyDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def add_mid_points(indexes: np.array):
#     newIndexes = []
#     for i, indexA in enumerate(indexes[:-1]):
#         newIndexes.append(indexA)
#         indexB = indexes[i+1]
#         midPointIndex = (indexA + indexB) // 2
#         newIndexes.append(midPointIndex)
#     newIndexes.append(indexes[-1])
    
#     return newIndexes

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

# def find_local_extrema(energies):
#     # Convert the energies to a numpy array
#     energiesArray = energies.to_numpy()

#     # Extend the array to handle periodicity
#     extendedEnergies = np.concatenate(
#         (energiesArray[-1:], energiesArray, energiesArray[:1])
#     )

#     # Find local minima
#     localMinimaIndexes = argrelextrema(extendedEnergies, np.less)[0] - 1

#     # Find local maxima
#     localMaximaIndexes = argrelextrema(extendedEnergies, np.greater)[0] - 1

#     # Filter out-of-bound indices
#     localMinimaIndexes = localMinimaIndexes[
#         (localMinimaIndexes >= 0) & (localMinimaIndexes < len(energiesArray))
#     ]
#     localMaximaIndexes = localMaximaIndexes[
#         (localMaximaIndexes >= 0) & (localMaximaIndexes < len(energiesArray))
#     ]

#     combinedExtremaIndexes = sorted(
#         np.concatenate((localMinimaIndexes, localMaximaIndexes))
#     )
#     return [int(index) for index in combinedExtremaIndexes]

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def read_singlepoint_energy(spOrcaOutput):
#     with open(spOrcaOutput, "r") as f:
#         for line in f:
#             if line.startswith("FINAL SINGLE POINT ENERGY"):
#                 singlePointEnergy = float(line.split()[-1])
#                 return singlePointEnergy

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def run_forwards_scan_step(optXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config):
#     ## XTB-GFN2 FORWARDS SCAN ##
#     ## do a forwards scan
#     forwardsDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_forwards")
#     os.makedirs(forwardsDir, exist_ok=True)
#     forwardsScanAngle = initialTorsionAngle + 360
#     forwardsScanText = f"{str(initialTorsionAngle)}, {str(forwardsScanAngle)}, {str(config['torsionScanInfo']['nScanSteps'])}"


#     forwardsOrcaInput: FilePath = drOrca.make_orca_input_for_scan(inputXyz=optXyz,
#                                                       outDir = forwardsDir,
#                                                         moleculeInfo = config["moleculeInfo"],
#                                                         qmMethod=config["torsionScanInfo"]["scanMethod"],
#                                                         solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
#                                                         torsionIndexes=torsionIndexes,
#                                                         scanAngles = forwardsScanText)
#     forwardsOrcaOutput: FilePath = p.join(forwardsDir, "orca_scan.out")
#     if not p.isfile(forwardsOrcaOutput):
#         os.chdir(forwardsDir)
#         run_orca(forwardsOrcaInput, forwardsOrcaOutput)
#     forwardsXyz = find_final_xyz(forwardsDir)
#     forwardsScanDf = read_data(forwardsDir)
#     forwardsScanDf.columns = ["Angle", "Energy"]
#     ## add orca numbers 
#     forwardsScanDf["scan_index"] = [f"{i:03}" for i in range(1, 38)]

#     forwardsScanDf = take_min_duplicate_angles(forwardsScanDf)

#     ## find scan xyzFiles,
#     return forwardsScanDf, forwardsXyz, forwardsDir
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def run_backwards_scan_step(forwardsScanXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config):
#     ## do a backwards scan
#     backwardsDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_backwards")
#     os.makedirs(backwardsDir, exist_ok=True)
#     forwardsScanAngle = initialTorsionAngle + 360 

#     backwardsScanText = f"{str(forwardsScanAngle)}, {str(initialTorsionAngle)}, {str(config['torsionScanInfo']['nScanSteps'])}"

#     backwardsOrcaInput: FilePath = drOrca.make_orca_input_for_scan(inputXyz=forwardsScanXyz,
#                                                         outDir = backwardsDir,
#                                                         moleculeInfo = config["moleculeInfo"],
#                                                         qmMethod=config["torsionScanInfo"]["scanMethod"],
#                                                         solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
#                                                         torsionIndexes=torsionIndexes,
#                                                         scanAngles = backwardsScanText)


#     backwardsOrcaOutput: FilePath = p.join(backwardsDir, "orca_scan.out")
#     if not p.isfile(backwardsOrcaOutput):
#         run_orca(backwardsOrcaInput, backwardsOrcaOutput)
#     backwardsScanDf = read_data(backwardsDir)
#     backwardsScanDf.columns = ["Angle", "Energy"]

#     backwardsScanDf["scan_index"] = [f"{i:03}" for i in range(1, 38)]


#     backwardsScanDf = take_min_duplicate_angles(backwardsScanDf)


#     return backwardsScanDf, backwardsDir  

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

# def find_scan_xyz_files(scanDir: DirectoryPath, expectedNumberOfFiles: int):
#     scanXyzs = sorted([xyzFile for xyzFile in glob.glob(p.join(scanDir, 'orca_scan.[0-9][0-9][0-9].xyz'))])
#     if not len(scanXyzs) == expectedNumberOfFiles:
#         raise(ValueError(f"Number of scan xyz files ({len(scanXyzs)}) does not match expected ({expectedNumberOfFiles})"))

#     return scanXyzs

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

# def take_min_duplicate_angles(energyDf):
#     energyDf['Angle'] = energyDf['Angle'].round(0)
#     # Find the index of the minimum energy for each unique Angle
#     minEnergyIndexes = energyDf.groupby('Angle')["Energy"].idxmin()

#     minEnergyIndexes = minEnergyIndexes.dropna()

#     # Select the rows with the minimum energy for each Angle
#     return energyDf.loc[minEnergyIndexes].reset_index(drop=True)

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def merge_scan_dfs(scanDfs):
#     mergedDf = pd.DataFrame()
#     mergedDf["Angle"] = np.arange(0,360,10)
#     for colIndex, scanDf in enumerate(scanDfs):
#         mergedDf = mergedDf.merge(scanDf[["Angle", "Energy"]], on = "Angle", how= "left")
#         mergedDf.rename(columns={"Energy":f"Energy_{colIndex + 1}"}, inplace=True)
#     return mergedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def process_energy_outputs(energyDf):
#     energyDf = rescale_and_sort_energy_angles(energyDf)
#     energyDf = take_min_duplicate_angles(energyDf)

#     energyDf["Energy"] = energyDf["Energy"].apply(hartree_to_kcal_per_mol)
#     return energyDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def hartree_to_kcal_per_mol(energy):
    return energy * 627.5095

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rescale_torsion_angles(angle):
    angle = angle % 360  # reduce the angle to the 0-360 range
    return angle

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def find_final_xyz(conformertorsionScanDir):
#     allXyzFiles = glob.glob(p.join(conformertorsionScanDir, "*.xyz"))
#     scanXYZFiles = sorted([f for f in allXyzFiles if re.match(r'orca_scan\.\d+\.xyz$', os.path.basename(f))])

#     finalXyzFile = scanXYZFiles[-1]
#     return finalXyzFile

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def measure_current_torsion_angle(conformerXyz, torsionIndexes):
#     atomCoords = xyz_to_np_array(conformerXyz)
#     torsionAngle = calculate_torsion_angle(atomCoords, torsionIndexes)
#     roundedToTenAngle = round(torsionAngle, -1)
#     return roundedToTenAngle

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def calculate_torsion_angle(coords, torsionIndexes):
#     # Vectors between points
#     b1 = coords[torsionIndexes[1]] - coords[torsionIndexes[0]]
#     b2 = coords[torsionIndexes[2]] - coords[torsionIndexes[1]]
#     b3 = coords[torsionIndexes[3]] - coords[torsionIndexes[2]]

#     # Normal vectors to the planes
#     n1 = np.cross(b1, b2)
#     n2 = np.cross(b2, b3)

#     # Normalize the normal vectors
#     n1 /= np.linalg.norm(n1)
#     n2 /= np.linalg.norm(n2)

#     # Calculate the angle
#     angleInRadians = np.arctan2(np.dot(np.cross(n1, n2), b2 / np.linalg.norm(b2)), np.dot(n1, n2))

#     # Convert from radians to degrees
#     angleIdDegrees = np.degrees(angleInRadians)

#     return angleIdDegrees
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
# def identify_rotatable_bonds(inputPdb) -> List[Tuple[int,int,int,int]]:
#     # Load the molecule from a PDB file
#     mol = Chem.MolFromPDBFile(inputPdb, removeHs=False)
#     # Identify torsion angles for rotatable bonds
#     torsionAngles = []
#     for bond in mol.GetBonds():
#         if bond.IsInRing():
#             continue
#         if bond.GetBondType() == Chem.BondType.SINGLE:
#             atom2 = bond.GetBeginAtom()
#             atom3 = bond.GetEndAtom()

#             ## get atom names for begin and end atoms
#             atom2Name = atom2.GetPDBResidueInfo().GetName().strip()
#             atom3Name = atom3.GetPDBResidueInfo().GetName().strip()

#             ## dont scan amide bonds
#             if atom2Name in ["NN", "C"] and atom3Name in ["NN", "C"]:
#                 continue
#             if atom2Name in ["N", "CC1"] and atom3Name in ["N", "CC1"]:
#                 continue

            
#             if not (atom2.IsInRing() or atom3.IsInRing()):
#                 # Find neighboring atoms for torsion angle
#                 neighborsBegin = [a for a in atom2.GetNeighbors() if a.GetIdx() != atom3.GetIdx()]
#                 neighborsEnd = [a for a in atom3.GetNeighbors() if a.GetIdx() != atom2.GetIdx()]
                
#                 for atom1 in neighborsBegin:
#                     for atom4 in neighborsEnd:

#                         atom1Name = atom1.GetPDBResidueInfo().GetName().strip()
#                         atom4Name = atom4.GetPDBResidueInfo().GetName().strip()

#                         ## dont scan bonds with non-polar hydrogens at either as atoms 1 or 4
#                         if atom1Name.startswith("H"):
#                             if atom2Name.startswith("C"):
#                                 continue
#                         if atom4Name.startswith("H"):
#                             if atom3Name.startswith("C"):
#                                 continue
#                         ## add torsion data to list
#                         torsionAngles.append({
#                             'atoms': (atom1Name, atom2Name, atom3Name, atom4Name),
#                             'indices': (atom1.GetIdx(), atom2.GetIdx(), atom3.GetIdx(), atom4.GetIdx())
#                         })


#     return torsionAngles

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def read_data(conformertorsionScanDir):
#     scanDat: FilePath = p.join(conformertorsionScanDir, "orca_scan.relaxscanact.dat")
#     if not os.path.exists(scanDat):
#         raise FileNotFoundError(f"Scan data not found in {conformertorsionScanDir}")
#     scanDf: pd.DataFrame = pd.read_csv(scanDat, sep='\s+', header=None)
#     return scanDf
    
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def run_orca(orcaInput, orcaOutput):
#     orcaCommand = ["orca", orcaInput]
#     with open(orcaOutput, 'w') as output_file:
#         try:
#             call(orcaCommand, stdout=output_file, stderr=output_file)
#         except Exception as e:
#             raise(e)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def xyz_to_np_array(filePath):
#     with open(filePath, 'r') as file:
#         lines = file.readlines()

    # # Skip the first two lines (atom count and comment)
    # dataLines = lines[2:]

    # # Parse the coordinates
    # coords = []
    # for line in dataLines:
    #     parts = line.split()
    #     if len(parts) >= 4:
    #         x, y, z = map(float, parts[1:4])
    #         coords.append([x, y, z])

    # return np.array(coords)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def detect_jumps_in_data(df):
#     # Calculate differences between consecutive values
#     diffDf = df.drop(columns='Angle').diff().abs()
    
#     # Identify columns with any difference greater than the threshold
#     jumpyCols = diffDf.columns[((diffDf > 10).any())]
    
#     # Drop these columns from the original DataFrame
#     cleanDf = df.drop(columns=jumpyCols)
    
#     return cleanDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# def calculate_means(df):
#     return df.drop(columns='Angle').mean(axis=1)

# if __name__ == "__main__":
#     configYaml = "/home/esp/scriptDevelopment/drFrankenstein/NMH_outputs/drFrankenstein.yaml"
#     with open(configYaml, "r") as yamlFile:
#         config = yaml.safe_load(yamlFile)
#     twist_protocol(config)