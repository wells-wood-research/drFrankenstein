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
def run_optimization_step(conformer_xyz: FilePath, 
                           conformer_scan_dir: DirectoryPath, 
                           conformer_id: str, 
                           config: dict) -> FilePath:
    """
    Runs a geometry optimization for a given conformer.

    Args:
        conformer_xyz: Path to the conformer XYZ file.
        conformer_scan_dir: Directory for this conformer's scan steps.
        conformer_id: Identifier for the conformer.
        config: Configuration dictionary.

    Returns:
        Path to the optimized XYZ file.
    """
    optDir: DirectoryPath = p.join(conformer_scan_dir, f"{conformer_id}_opt")
    os.makedirs(optDir, exist_ok=True)
    ## make an ORCA input file for optimisation
    optOrcaInput: FilePath = drOrca.make_orca_input_for_opt(inputXyz=conformer_xyz,
                                                   outDir = optDir,moleculeInfo=config["moleculeInfo"],
                                                   qmMethod=config["torsionScanInfo"]["scanMethod"],
                                                    solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"])
    flags = [p.join(optDir, flag) for flag in ["ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"]] # Corrected path for flags
                        
    if not any ([p.isfile(flag) for flag in flags]):
        optOrcaOutput: FilePath = p.join(optDir, "orca_opt.out")
        drOrca.run_orca(optOrcaInput, optOrcaOutput, config)
        Twisted_Assistant.create_orca_terminated_flag(optDir, optOrcaOutput) # type: ignore

    optXyz = p.join(optDir, "orca_opt.xyz")

    cleaner.clean_up_opt_dir(optDir, config)
    return optXyz

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_forwards_scan_step(opt_xyz: FilePath, 
                            initial_torsion_angle: float, 
                            torsion_indexes: List[int], 
                            conformer_scan_dir: DirectoryPath, 
                            conformer_id: str, 
                            config: dict) -> Tuple[pd.DataFrame, FilePath, DirectoryPath]:
    """
    Runs the forward part of a torsion scan for a given conformer.

    Args:
        opt_xyz: Path to the optimized conformer XYZ file.
        initial_torsion_angle: The starting torsion angle.
        torsion_indexes: List of atom indexes defining the torsion.
        conformer_scan_dir: Directory for this conformer's scan steps.
        conformer_id: Identifier for the conformer.
        config: Configuration dictionary.

    Returns:
        A tuple containing the DataFrame of scan results, path to the final XYZ, and the forward scan directory.
    """
    ## XTB-GFN2 FORWARDS SCAN ##
    ## do a forwards scan
    forwardsDir: DirectoryPath = p.join(conformer_scan_dir, f"{conformer_id}_forwards")
    os.makedirs(forwardsDir, exist_ok=True)
    forwardsScanAngle = initial_torsion_angle + 360
    forwardsScanText = f"{str(initial_torsion_angle)}, {str(forwardsScanAngle)}, {str(config['torsionScanInfo']['nScanSteps'])}"


    forwardsOrcaInput: FilePath = drOrca.make_orca_input_for_scan(inputXyz=opt_xyz,
                                                      outDir = forwardsDir,
                                                        moleculeInfo = config["moleculeInfo"],
                                                        qmMethod=config["torsionScanInfo"]["scanMethod"],
                                                        solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
                                                        torsionIndexes=torsion_indexes,
                                                        scanAngles = forwardsScanText)

    
    flags = [p.join(forwardsDir, flag) for flag in ["ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"]]

    if not any ([p.isfile(flag) for flag in flags]):
        os.chdir(forwardsDir)
        forwardsOrcaOutput: FilePath = p.join(forwardsDir, "orca_scan.out")
        drOrca.run_orca(forwardsOrcaInput, forwardsOrcaOutput, config)
        Twisted_Assistant.create_orca_terminated_flag(forwardsDir, forwardsOrcaOutput) # type: ignore
    
    forwardsXyz = Twisted_Assistant.find_final_xyz(forwardsDir) # type: ignore
    forwardsScanDf = Twisted_Assistant.read_scan_energy_data(forwardsDir) # type: ignore
    forwardsScanDf.columns = ["Angle", "Energy"]
    ## add orca numbers 
    forwardsScanDf["scan_index"] = [f"{i:03}" for i in range(1, config['torsionScanInfo']['nScanSteps'] + 1)] # Corrected range

    forwardsScanDf = Twisted_Assistant.take_min_duplicate_angles(forwardsScanDf) # type: ignore
    cleaner.clean_up_scan_dir(forwardsDir, config)

    ## find scan xyzFiles,
    return forwardsScanDf, forwardsXyz, forwardsDir
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def run_backwards_scan_step(forwards_scan_xyz: FilePath, 
                             initial_torsion_angle: float, 
                             torsion_indexes: List[int], 
                             conformer_scan_dir: DirectoryPath, 
                             conformer_id: str, 
                             config: dict) -> Tuple[pd.DataFrame, DirectoryPath]:
    """
    Runs the backward part of a torsion scan for a given conformer.

    Args:
        forwards_scan_xyz: Path to the final XYZ from the forward scan.
        initial_torsion_angle: The starting torsion angle for the original forward scan.
        torsion_indexes: List of atom indexes defining the torsion.
        conformer_scan_dir: Directory for this conformer's scan steps.
        conformer_id: Identifier for the conformer.
        config: Configuration dictionary.

    Returns:
        A tuple containing the DataFrame of scan results and the backward scan directory.
    """
    ## do a backwards scan
    backwardsDir: DirectoryPath = p.join(conformer_scan_dir, f"{conformer_id}_backwards")
    os.makedirs(backwardsDir, exist_ok=True)
    forwardsScanAngle = initial_torsion_angle + 360 # This is the end point of the forward scan, start for backward

    backwardsScanText = f"{str(forwardsScanAngle)}, {str(initial_torsion_angle)}, {str(config['torsionScanInfo']['nScanSteps'])}"

    backwardsOrcaInput: FilePath = drOrca.make_orca_input_for_scan(inputXyz=forwards_scan_xyz,
                                                        outDir = backwardsDir,
                                                        moleculeInfo = config["moleculeInfo"],
                                                        qmMethod=config["torsionScanInfo"]["scanMethod"],
                                                        solvationMethod=config["torsionScanInfo"]["scanSolvationMethod"],
                                                        torsionIndexes=torsion_indexes,
                                                        scanAngles = backwardsScanText)


    flags = [p.join(backwardsDir, flag) for flag in ["ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"]]

    if not any ([p.isfile(flag) for flag in flags]):
        os.chdir(backwardsDir)
        backwardsOrcaOutput: FilePath = p.join(backwardsDir, "orca_scan.out")
        drOrca.run_orca(backwardsOrcaInput, backwardsOrcaOutput, config)
        Twisted_Assistant.create_orca_terminated_flag(backwardsDir, backwardsOrcaOutput) # type: ignore

    backwardsScanDf = Twisted_Assistant.read_scan_energy_data(backwardsDir) # type: ignore
    backwardsScanDf.columns = ["Angle", "Energy"]

    backwardsScanDf["scan_index"] = [f"{i:03}" for i in range(1, config['torsionScanInfo']['nScanSteps'] + 1)] # Corrected range

    backwardsScanDf = Twisted_Assistant.take_min_duplicate_angles(backwardsScanDf) # type: ignore

    cleaner.clean_up_scan_dir(backwardsDir, config)

    return backwardsScanDf, backwardsDir  
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_singlepoints_on_scans(scan_dir: DirectoryPath, 
                               scan_df: pd.DataFrame, 
                               conformer_id: str,  
                               config: dict) -> pd.DataFrame:
    """
    Runs single point energy calculations on selected geometries from a torsion scan.

    Args:
        scan_dir: Directory containing the scan XYZ files.
        scan_df: DataFrame with the energies and angles from the scan.
        conformer_id: Identifier for the conformer.
        config: Configuration dictionary.

    Returns:
        DataFrame with single point energies.
    """
    scanXyzs = Twisted_Assistant.find_scan_xyz_files(scan_dir, expectedNumberOfFiles=config["torsionScanInfo"]["nScanSteps"]) # type: ignore
    singlePointsOn = config["torsionScanInfo"]["scanSinglePointsOn"]
    if singlePointsOn == "all":
        stationaryPointScanIndexes = scan_df["scan_index"].to_list()

    elif singlePointsOn == "minMaxOnly":
        stationaryPointIndexes = Twisted_Assistant.find_local_extrema(scan_df["Energy"]) # type: ignore
        stationaryPointScanIndexes = scan_df.loc[stationaryPointIndexes, "scan_index"].to_list()
        scanXyzs = [scanXyz for scanXyz in scanXyzs if scanXyz.split(".")[1] in stationaryPointScanIndexes]

    elif singlePointsOn == "minMaxMiddle":
        stationaryPointIndexes = Twisted_Assistant.find_local_extrema(scan_df["Energy"]) # type: ignore
        stationaryAndMidPointIndexes = Twisted_Assistant.add_mid_points(stationaryPointIndexes) # type: ignore
        stationaryPointScanIndexes = scan_df.loc[stationaryAndMidPointIndexes, "scan_index"].to_list()
        scanXyzs = [scanXyz for scanXyz in scanXyzs if scanXyz.split(".")[1] in stationaryPointScanIndexes]


    conformerScanDir = p.dirname(scan_dir)
    tag = scan_dir.split("_")[-1]
    spTopDir = p.join(conformerScanDir, f"SP_{tag}_{conformer_id}")
    os.makedirs(spTopDir, exist_ok=True)


    singlePointEnergies = {}
    for scanXyz in scanXyzs:
        scanId = scanXyz.split(".")[1]
        spDir = p.join(spTopDir, f"SP_{conformer_id}_{tag}_{scanId}")
        os.makedirs(spDir, exist_ok=True)
        spOrcaInput: FilePath = drOrca.make_orca_input_for_singlepoint(inputXyz=scanXyz,
                                                                       outDir= spDir,
                                                                       moleculeInfo=config["moleculeInfo"],
                                                                       qmMethod=config["torsionScanInfo"]["singlePointMethod"],
                                                                       solvationMethod=config["torsionScanInfo"]["singlePointSolvationMethod"])
        spOrcaOutput : FilePath = p.join(spDir, "orca_sp.out")
        if not p.isfile(spOrcaOutput):
            drOrca.run_orca(spOrcaInput, spOrcaOutput, config)
        singlePointEnergy = Twisted_Assistant.read_singlepoint_energy(spOrcaOutput) # type: ignore
        singlePointEnergies[scanId] = singlePointEnergy
        cleaner.clean_up_singlepoint_dir(spDir, config)
    singlePointEnergyDf = pd.DataFrame(list(singlePointEnergies.items()), columns = ["scan_index", "Energy"])

    singlePointEnergyDf = singlePointEnergyDf.merge(
        scan_df[["scan_index", "Angle"]], on="scan_index", how="left")
    singlePointEnergyDf = Twisted_Assistant.process_energy_outputs(singlePointEnergyDf) # type: ignore
    singlePointEnergyDf.to_csv(p.join(spTopDir, f"SP_{conformer_id}_{tag}.csv"), index=False)
    return singlePointEnergyDf