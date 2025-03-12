import os
from os import path as p
import sys
import pandas as pd
## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import Dict, List, Tuple
## ADD SRC TO PATH ##
currentFilePath: FilePath = os.path.abspath(__file__)
currentDir: DirectoryPath = os.path.dirname(currentFilePath)
srcDir: DirectoryPath = os.path.dirname(currentDir)
sys.path.append(srcDir)

## drFRANKENSTEIN MODULES ##
import drOrca
from drTwist import Assistant

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_optimisation_step(conformerXyz, conformerScanDir, conformerId, config):
    optDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_opt")
    os.makedirs(optDir, exist_ok=True)
    ## make an ORCA input file for optimisation
    optOrcaInput: FilePath = drOrca.make_orca_input_for_opt(inputXyz=conformerXyz,
                                                   outDir = optDir,moleculeInfo=config["moleculeInfo"],
                                                   qmMethod=config["torsionScanInfo"]["scanMethod"],
                                                    solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"])
                                                    
    optOrcaOutput: FilePath = p.join(optDir, "orca_opt.out")
    if not p.isfile(optOrcaOutput):
        drOrca.run_orca(optOrcaInput, optOrcaOutput)
    optXyz = p.join(optDir, "orca_opt.xyz")
    return optXyz

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_forwards_scan_step(optXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config):
    ## XTB-GFN2 FORWARDS SCAN ##
    ## do a forwards scan
    forwardsDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_forwards")
    os.makedirs(forwardsDir, exist_ok=True)
    forwardsScanAngle = initialTorsionAngle + 360
    forwardsScanText = f"{str(initialTorsionAngle)}, {str(forwardsScanAngle)}, {str(config['torsionScanInfo']['nScanSteps'])}"


    forwardsOrcaInput: FilePath = drOrca.make_orca_input_for_scan(inputXyz=optXyz,
                                                      outDir = forwardsDir,
                                                        moleculeInfo = config["moleculeInfo"],
                                                        qmMethod=config["torsionScanInfo"]["scanMethod"],
                                                        solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
                                                        torsionIndexes=torsionIndexes,
                                                        scanAngles = forwardsScanText)
    forwardsOrcaOutput: FilePath = p.join(forwardsDir, "orca_scan.out")
    if not p.isfile(forwardsOrcaOutput):
        os.chdir(forwardsDir)
        drOrca.run_orca(forwardsOrcaInput, forwardsOrcaOutput)
    forwardsXyz = Assistant.find_final_xyz(forwardsDir)
    forwardsScanDf = Assistant.read_scan_energy_data(forwardsDir)
    forwardsScanDf.columns = ["Angle", "Energy"]
    ## add orca numbers 
    forwardsScanDf["scan_index"] = [f"{i:03}" for i in range(1, 38)]

    forwardsScanDf = Assistant.take_min_duplicate_angles(forwardsScanDf)

    ## find scan xyzFiles,
    return forwardsScanDf, forwardsXyz, forwardsDir
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def run_backwards_scan_step(forwardsScanXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config):
    ## do a backwards scan
    backwardsDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_backwards")
    os.makedirs(backwardsDir, exist_ok=True)
    forwardsScanAngle = initialTorsionAngle + 360 

    backwardsScanText = f"{str(forwardsScanAngle)}, {str(initialTorsionAngle)}, {str(config['torsionScanInfo']['nScanSteps'])}"

    backwardsOrcaInput: FilePath = drOrca.make_orca_input_for_scan(inputXyz=forwardsScanXyz,
                                                        outDir = backwardsDir,
                                                        moleculeInfo = config["moleculeInfo"],
                                                        qmMethod=config["torsionScanInfo"]["scanMethod"],
                                                        solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
                                                        torsionIndexes=torsionIndexes,
                                                        scanAngles = backwardsScanText)


    backwardsOrcaOutput: FilePath = p.join(backwardsDir, "orca_scan.out")
    if not p.isfile(backwardsOrcaOutput):
        drOrca.run_orca(backwardsOrcaInput, backwardsOrcaOutput)
    backwardsScanDf = Assistant.read_scan_energy_data(backwardsDir)
    backwardsScanDf.columns = ["Angle", "Energy"]

    backwardsScanDf["scan_index"] = [f"{i:03}" for i in range(1, 38)]


    backwardsScanDf = Assistant.take_min_duplicate_angles(backwardsScanDf)


    return backwardsScanDf, backwardsDir  
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_singlepoints_on_scans(scanDir, scanDf, outDir, conformerId,  config, tag):

    scanXyzs = Assistant.find_scan_xyz_files(scanDir, expectedNumberOfFiles=config["torsionScanInfo"]["nScanSteps"])

    singlePointsOn = config["torsionScanInfo"]["scanSinglePointsOn"]
    if singlePointsOn == "all":
        stationaryPointScanIndexes = scanDf["scan_index"].to_list()

    elif singlePointsOn == "minMaxOnly":
        stationaryPointsIndexes = Assistant.find_local_extrema(scanDf["Energy"])
        stationaryPointScanIndexes = scanDf.loc[stationaryPointsIndexes, "scan_index"].to_list()
        scanXyzs = [scanXyz for scanXyz in scanXyzs if scanXyz.split(".")[1] in stationaryPointScanIndexes]

    elif singlePointsOn == "minMaxMiddle":
        stationaryPointsIndexes = Assistant.find_local_extrema(scanDf["Energy"])
        stationaryAndMidPointIndexes = Assistant.add_mid_points(stationaryPointsIndexes)
        stationaryPointScanIndexes = scanDf.loc[stationaryAndMidPointIndexes, "scan_index"].to_list()
        scanXyzs = [scanXyz for scanXyz in scanXyzs if scanXyz.split(".")[1] in stationaryPointScanIndexes]

    singlePointEnergies = {}
    for scanXyz in scanXyzs:
        scanId = scanXyz.split(".")[1]
        scanDir = p.join(outDir, f"SP_{conformerId}_{tag}_{scanId}")
        os.makedirs(scanDir, exist_ok=True)
        spOrcaInput: FilePath = drOrca.make_orca_input_for_singlepoint(inputXyz=scanXyz,
                                                                       outDir= scanDir,
                                                                       moleculeInfo=config["moleculeInfo"],
                                                                       qmMethod=config["torsionScanInfo"]["singlePointMethod"],
                                                                       solvationMethod=config["torsionScanInfo"]["singlePointSolvationMethod"])
        spOrcaOutput : FilePath = p.join(scanDir, "orca_sp.out")
        if not p.isfile(spOrcaOutput):
            drOrca.run_orca(spOrcaInput, spOrcaOutput)
        singlePointEnergy = Assistant.read_singlepoint_energy(spOrcaOutput)
        singlePointEnergies[scanId] = singlePointEnergy
    singlePointEnergyDf = pd.DataFrame(list(singlePointEnergies.items()), columns = ["scan_index", "Energy"])

    singlePointEnergyDf = singlePointEnergyDf.merge(
        scanDf[["scan_index", "Angle"]], on="scan_index", how="left")

    return singlePointEnergyDf