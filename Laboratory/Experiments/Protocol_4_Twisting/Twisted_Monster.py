import os
from os import path as p
import pandas as pd

## drFRANKENSTEIN LIBRARIES ##
from OperatingTools import cleaner

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import Dict, List, Tuple


## drFRANKENSTEIN MODULES ##
from . import Twisted_Assistant
from OperatingTools import drOrca, Timer

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_optimisation_step(conformerXyz, conformerScanDir, conformerId, config):
    optDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_opt")
    os.makedirs(optDir, exist_ok=True)
    ## make an ORCA input file for optimisation
    optOrcaInput: FilePath = drOrca.make_orca_input_for_opt(inputXyz=conformerXyz,
                                                   outDir = optDir,moleculeInfo=config["moleculeInfo"],
                                                   qmMethod=config["torsionScanInfo"]["scanMethod"],
                                                    solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"])
    flags = [p.join(optOrcaInput, flag) for flag in ["ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"]]
                        
    if not any ([p.isfile(flag) for flag in flags]):
        optOrcaOutput: FilePath = p.join(optDir, "orca_opt.out")
        drOrca.run_orca(optOrcaInput, optOrcaOutput, config)
        Twisted_Assistant.create_orca_terminated_flag(optDir, optOrcaOutput)

    optXyz = p.join(optDir, "orca_opt.xyz")

    cleaner.clean_up_opt_dir(optDir, config)
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

    
    flags = [p.join(forwardsDir, flag) for flag in ["ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"]]

    if not any ([p.isfile(flag) for flag in flags]):
        os.chdir(forwardsDir)
        forwardsOrcaOutput: FilePath = p.join(forwardsDir, "orca_scan.out")
        drOrca.run_orca(forwardsOrcaInput, forwardsOrcaOutput, config)
        Twisted_Assistant.create_orca_terminated_flag(forwardsDir, forwardsOrcaOutput)
    
    forwardsXyz = Twisted_Assistant.find_final_xyz(forwardsDir)
    forwardsScanDf = Twisted_Assistant.read_scan_energy_data(forwardsDir)
    forwardsScanDf.columns = ["Angle", "Energy"]
    ## add orca numbers 
    forwardsScanDf["scan_index"] = [f"{i:03}" for i in range(1, 38)]

    forwardsScanDf = Twisted_Assistant.take_min_duplicate_angles(forwardsScanDf)
    cleaner.clean_up_scan_dir(forwardsDir, config)

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


    flags = [p.join(backwardsDir, flag) for flag in ["ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"]]

    if not any ([p.isfile(flag) for flag in flags]):
        os.chdir(backwardsDir)
        backwardsOrcaOutput: FilePath = p.join(backwardsDir, "orca_scan.out")
        drOrca.run_orca(backwardsOrcaInput, backwardsOrcaOutput, config)
        Twisted_Assistant.create_orca_terminated_flag(backwardsDir, backwardsOrcaOutput)

    backwardsScanDf = Twisted_Assistant.read_scan_energy_data(backwardsDir)
    backwardsScanDf.columns = ["Angle", "Energy"]

    backwardsScanDf["scan_index"] = [f"{i:03}" for i in range(1, 38)]

    backwardsScanDf = Twisted_Assistant.take_min_duplicate_angles(backwardsScanDf)

    cleaner.clean_up_scan_dir(backwardsDir, config)

    return backwardsScanDf, backwardsDir  
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_singlepoints_on_scans(scanDir, scanDf, conformerId,  config):

    scanXyzs = Twisted_Assistant.find_scan_xyz_files(scanDir, expectedNumberOfFiles=config["torsionScanInfo"]["nScanSteps"])
    singlePointsOn = config["torsionScanInfo"]["scanSinglePointsOn"]
    if singlePointsOn == "all":
        stationaryPointScanIndexes = scanDf["scan_index"].to_list()

    elif singlePointsOn == "minMaxOnly":
        stationaryPointsIndexes = Twisted_Assistant.find_local_extrema(scanDf["Energy"])
        stationaryPointScanIndexes = scanDf.loc[stationaryPointsIndexes, "scan_index"].to_list()
        scanXyzs = [scanXyz for scanXyz in scanXyzs if scanXyz.split(".")[1] in stationaryPointScanIndexes]

    elif singlePointsOn == "minMaxMiddle":
        stationaryPointsIndexes = Twisted_Assistant.find_local_extrema(scanDf["Energy"])
        stationaryAndMidPointIndexes = Twisted_Assistant.add_mid_points(stationaryPointsIndexes)
        stationaryPointScanIndexes = scanDf.loc[stationaryAndMidPointIndexes, "scan_index"].to_list()
        scanXyzs = [scanXyz for scanXyz in scanXyzs if scanXyz.split(".")[1] in stationaryPointScanIndexes]


    conformerScanDir = p.dirname(scanDir)
    tag = scanDir.split("_")[-1]
    spTopDir = p.join(conformerScanDir, f"SP_{tag}_{conformerId}")
    os.makedirs(spTopDir, exist_ok=True)


    singlePointEnergies = {}
    for scanXyz in scanXyzs:
        scanId = scanXyz.split(".")[1]
        spDir = p.join(spTopDir, f"SP_{conformerId}_{tag}_{scanId}")
        os.makedirs(spDir, exist_ok=True)
        spOrcaInput: FilePath = drOrca.make_orca_input_for_singlepoint(inputXyz=scanXyz,
                                                                       outDir= spDir,
                                                                       moleculeInfo=config["moleculeInfo"],
                                                                       qmMethod=config["torsionScanInfo"]["singlePointMethod"],
                                                                       solvationMethod=config["torsionScanInfo"]["singlePointSolvationMethod"])
        spOrcaOutput : FilePath = p.join(spDir, "orca_sp.out")
        if not p.isfile(spOrcaOutput):
            drOrca.run_orca(spOrcaInput, spOrcaOutput, config)
        singlePointEnergy = Twisted_Assistant.read_singlepoint_energy(spOrcaOutput)
        singlePointEnergies[scanId] = singlePointEnergy
        cleaner.clean_up_singlepoint_dir(spDir, config)
    singlePointEnergyDf = pd.DataFrame(list(singlePointEnergies.items()), columns = ["scan_index", "Energy"])

    singlePointEnergyDf = singlePointEnergyDf.merge(
        scanDf[["scan_index", "Angle"]], on="scan_index", how="left")
    singlePointEnergyDf = Twisted_Assistant.process_energy_outputs(singlePointEnergyDf)
    singlePointEnergyDf.to_csv(p.join(spTopDir, f"SP_{conformerId}_{tag}.csv"), index=False)

    
    return singlePointEnergyDf