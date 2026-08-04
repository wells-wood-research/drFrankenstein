import os
import fnmatch
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


def _orca_hit_optimisation_iteration_limit(orcaOutput: FilePath) -> bool:
    """Return whether ORCA stopped after exhausting optimisation iterations."""
    if not p.isfile(orcaOutput):
        return False

    failureMessage = "*** convergence criteria cannot be satisfied within * iterations ***"
    with open(orcaOutput, "r") as f:
        return any(fnmatch.fnmatch(line.strip().lower(), failureMessage) for line in f)


def _retry_crashed_scan_with_copt(inputXyz, scanDir, moleculeInfo, qmMethod,
                                  solvationMethod, torsionIndexes, scanAngles,
                                  nCpus, config) -> bool:
    """Retry a scan with COPT when its only detected failure is iteration exhaustion."""
    crashedFlag = p.join(scanDir, "ORCA_CRASHED")
    orcaOutput = p.join(scanDir, "orca_scan.out")
    if not (p.isfile(crashedFlag) and _orca_hit_optimisation_iteration_limit(orcaOutput)):
        return False

    # Preserve scan geometries and ancillary ORCA files; replace only this run's
    # input, output, and status so the COPT calculation is a clean retry.
    for filePath in (p.join(scanDir, "orca_scan.inp"), orcaOutput, crashedFlag):
        if p.isfile(filePath):
            os.remove(filePath)

    coptInput = drOrca.make_orca_input_for_scan(
        inputXyz=inputXyz, outDir=scanDir, moleculeInfo=moleculeInfo,
        qmMethod=qmMethod, solvationMethod=solvationMethod,
        torsionIndexes=torsionIndexes, scanAngles=scanAngles, nCpus=nCpus,
        fallback=True,
    )
    drOrca.run_orca(coptInput, orcaOutput, config)
    Twisted_Assistant.create_orca_terminated_flag(scanDir, orcaOutput)
    return True

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_optimisation_step(conformerXyz, conformerScanDir, conformerId, config):
    optDir: DirectoryPath = p.join(conformerScanDir, f"{conformerId}_opt")
    os.makedirs(optDir, exist_ok=True)
    nCoresPerCalculation = config["torsionScanInfo"].get("nCoresPerCalculation", 1)
    ## make an ORCA input file for optimisation
    optOrcaInput: FilePath = drOrca.make_orca_input_for_opt(inputXyz=conformerXyz,
                                                   outDir = optDir,moleculeInfo=config["moleculeInfo"],
                                                   qmMethod=config["torsionScanInfo"]["scanMethod"],
                                                    solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
                                                    nCpus=nCoresPerCalculation)
    finishedFlag = p.join(optDir, "ORCA_FINISHED_NORMALLY")
    crashedFlag = p.join(optDir, "ORCA_CRASHED")
    flags = [finishedFlag, crashedFlag]
                        
    if not any ([p.isfile(flag) for flag in flags]):
        optOrcaOutput: FilePath = p.join(optDir, "orca_opt.out")
        drOrca.run_orca(optOrcaInput, optOrcaOutput, config)
        Twisted_Assistant.create_orca_terminated_flag(optDir, optOrcaOutput)
    if p.isfile(crashedFlag):
        raise FileNotFoundError(f"ORCA crashed during optimisation for {conformerId}: {optDir}")

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
    forwardsScanText = f"{str(initialTorsionAngle)}, {str(forwardsScanAngle)}, 37"


    forwardsOrcaInput: FilePath = drOrca.make_orca_input_for_scan(inputXyz=optXyz,
                                                      outDir = forwardsDir,
                                                        moleculeInfo = config["moleculeInfo"],
                                                        qmMethod=config["torsionScanInfo"]["scanMethod"],
                                                        solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
                                                        torsionIndexes=torsionIndexes,
                                                        scanAngles = forwardsScanText,
                                                        nCpus=config["torsionScanInfo"].get("nCoresPerCalculation", 1))

    
    finishedFlag = p.join(forwardsDir, "ORCA_FINISHED_NORMALLY")
    crashedFlag = p.join(forwardsDir, "ORCA_CRASHED")
    flags = [finishedFlag, crashedFlag]

    if not any ([p.isfile(flag) for flag in flags]):
        os.chdir(forwardsDir)
        forwardsOrcaOutput: FilePath = p.join(forwardsDir, "orca_scan.out")
        drOrca.run_orca(forwardsOrcaInput, forwardsOrcaOutput, config)
        Twisted_Assistant.create_orca_terminated_flag(forwardsDir, forwardsOrcaOutput)
    _retry_crashed_scan_with_copt(
        inputXyz=optXyz, scanDir=forwardsDir, moleculeInfo=config["moleculeInfo"],
        qmMethod=config["torsionScanInfo"]["scanMethod"],
        solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
        torsionIndexes=torsionIndexes, scanAngles=forwardsScanText,
        nCpus=config["torsionScanInfo"].get("nCoresPerCalculation", 1), config=config,
    )
    if p.isfile(crashedFlag):
        raise FileNotFoundError(f"ORCA crashed during forwards scan for {conformerId}: {forwardsDir}")
    
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

    backwardsScanText = f"{str(forwardsScanAngle)}, {str(initialTorsionAngle)}, 37"

    backwardsOrcaInput: FilePath = drOrca.make_orca_input_for_scan(inputXyz=forwardsScanXyz,
                                                        outDir = backwardsDir,
                                                        moleculeInfo = config["moleculeInfo"],
                                                        qmMethod=config["torsionScanInfo"]["scanMethod"],
                                                        solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
                                                        torsionIndexes=torsionIndexes,
                                                        scanAngles = backwardsScanText,
                                                        nCpus=config["torsionScanInfo"].get("nCoresPerCalculation", 1))


    finishedFlag = p.join(backwardsDir, "ORCA_FINISHED_NORMALLY")
    crashedFlag = p.join(backwardsDir, "ORCA_CRASHED")
    flags = [finishedFlag, crashedFlag]

    if not any ([p.isfile(flag) for flag in flags]):
        os.chdir(backwardsDir)
        backwardsOrcaOutput: FilePath = p.join(backwardsDir, "orca_scan.out")
        drOrca.run_orca(backwardsOrcaInput, backwardsOrcaOutput, config)
        Twisted_Assistant.create_orca_terminated_flag(backwardsDir, backwardsOrcaOutput)
    _retry_crashed_scan_with_copt(
        inputXyz=forwardsScanXyz, scanDir=backwardsDir, moleculeInfo=config["moleculeInfo"],
        qmMethod=config["torsionScanInfo"]["scanMethod"],
        solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
        torsionIndexes=torsionIndexes, scanAngles=backwardsScanText,
        nCpus=config["torsionScanInfo"].get("nCoresPerCalculation", 1), config=config,
    )
    if p.isfile(crashedFlag):
        raise FileNotFoundError(f"ORCA crashed during backwards scan for {conformerId}: {backwardsDir}")

    backwardsScanDf = Twisted_Assistant.read_scan_energy_data(backwardsDir)
    backwardsScanDf.columns = ["Angle", "Energy"]

    backwardsScanDf["scan_index"] = [f"{i:03}" for i in range(1, 38)]

    backwardsScanDf = Twisted_Assistant.take_min_duplicate_angles(backwardsScanDf)

    cleaner.clean_up_scan_dir(backwardsDir, config)

    return backwardsScanDf, backwardsDir  
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_singlepoints_on_scans(scanDir, scanDf, conformerId,  config):

    scanXyzs = Twisted_Assistant.find_scan_xyz_files(scanDir, 37)

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
                                                                       solvationMethod=config["torsionScanInfo"]["singlePointSolvationMethod"],
                                                                       nCpus=config["torsionScanInfo"].get("nCoresPerCalculation", 1))
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
