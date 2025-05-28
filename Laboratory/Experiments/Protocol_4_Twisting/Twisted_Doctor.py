
import os
from os import path as p
import pandas as pd

## MULTIPROCESSING AND LOADING BAR LIBRARIES ##
from mpire import WorkerPool
from mpire.utils import make_single_arguments

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import Dict, List, Tuple


## drFRANKENSTEIN MODULES ##
from . import Twisted_Assistant
from . import Twisted_Monster
from . import Twisted_Plotter
from OperatingTools import drSplash
from OperatingTools import Timer
from typing import Optional # Added for Optional return types

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def twist_protocol(config: dict) -> dict:
    """
    Main protocol to perform torsion scanning for identified rotatable bonds.

    Args:
        config: Dictionary containing all run information and settings.

    Returns:
        The updated configuration dictionary with results from torsion scanning.
    """
    ## create an entry in runtimeInfo for twist
    config["runtimeInfo"]["madeByTwisting"] = {}
    config["runtimeInfo"]["madeByTwisting"]["torsionDirs"] = []
    config["runtimeInfo"]["madeByTwisting"]["torsionTags"] = []
    config["runtimeInfo"]["madeByTwisting"]["finalScanEnergies"] = {}
    config["runtimeInfo"]["madeByTwisting"]["scanEnergyData"] = {}

    config = Twisted_Assistant.set_up_directories(config)

    ## identify rotatable bonds to be scanned
    config = Twisted_Assistant.identify_rotatable_bonds(config, mode = config["parameterFittingInfo"]["forceField"])

    ## exclude backbone torsions if requested
    if config["torsionScanInfo"]["preserveBackboneTorsions"]:
        config = Twisted_Assistant.exclude_backbone_torsions(config)

    rotatableDihedrals = config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"]

    nRotatableBonds = len(rotatableDihedrals)
    config["runtimeInfo"]["madeByTwisting"]["nRotatableBonds"] = nRotatableBonds
    for torsionIndex, (torsionTag, torsionData) in enumerate(rotatableDihedrals.items()):
        config = run_torsion_scanning(torsionTag, torsionData, torsionIndex, nRotatableBonds, config)
    config["checkpointInfo"]["scanningComplete"] = True

    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_torsion_scanning(torsionTag: str,
                          torsionData: Dict[str, Dict[str, Tuple[str,...]]], # Corrected type hint
                            torsionIndex: int,
                              nRotatableBonds: int,
                                config: dict,
                                  debug: bool = False) -> dict:

    ## make a top level dir for this torsion
    torsionDir = p.join(config["runtimeInfo"]["madeByTwisting"]["torsionDir"], f"torsion_{torsionTag}" )
    os.makedirs(torsionDir, exist_ok=True)
    ## add to config
    config["runtimeInfo"]["madeByTwisting"]["torsionDirs"].append(torsionDir)
    config["runtimeInfo"]["madeByTwisting"]["torsionTags"].append(torsionTag)

    conformerXyzs = Twisted_Assistant.get_conformer_xyzs(config)
    drSplash.show_torsion_being_scanned(torsionTag, torsionIndex, nRotatableBonds)
    if debug:
        ## run in serial
        scanDfs, scanDirs = scan_in_serial(torsionDir, conformerXyzs, torsionData["ATOM_INDEXES"], config = config)
        if config["torsionScanInfo"]["singlePointMethod"] is None:
            singlePointDfs = None
        else:
            singlePointDfs = single_points_in_serial(scanDirs, scanDfs, torsionDir, torsionTag,  config = config)
    else:
        # run torsion scans in parallel
        scanDfs, scanDirs = scan_in_parallel(torsionDir, conformerXyzs, torsionData["ATOM_INDEXES"], torsionTag, config = config)
        ## run single point scans in parallel
        if config["torsionScanInfo"]["scanSinglePointsOn"] is None or config["torsionScanInfo"]["singlePointMethod"] == None:
            singlePointDfs = None
        else:
            singlePointDfs = single_points_in_parallel(scanDirs, scanDfs, torsionTag, config = config)
    ## Merge scan data, calculate averages, rolling averages and mean average errors
    scanEnergiesCsv, scanAveragesDf  = Twisted_Assistant.process_scan_data(scanDfs, torsionDir, torsionTag)
    if  config["torsionScanInfo"]["singlePointMethod"] is None:
        config["runtimeInfo"]["madeByTwisting"]["finalScanEnergies"][torsionTag] = scanEnergiesCsv
        singlePointAveragesDf = None
        config = Twisted_Assistant.gather_scan_data(scanAveragesDf, torsionTag, config)
    else:
        singlePointEnergiesCsv, singlePointAveragesDf = Twisted_Assistant.process_scan_data(singlePointDfs, torsionDir, torsionTag)
        config["runtimeInfo"]["madeByTwisting"]["finalScanEnergies"][torsionTag] = singlePointEnergiesCsv   
        config = Twisted_Assistant.gather_scan_data(singlePointAveragesDf, torsionTag, config)

    ## Plotting
    config = Twisted_Plotter.twist_plotting_protocol(scanDfs, scanAveragesDf, singlePointDfs, singlePointAveragesDf, torsionDir, torsionTag, config)
    

    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Scan Single Points", "TORSION_SCANS")
def single_points_in_serial(scanDirs: List[DirectoryPath], 
                             scanDfs: List[pd.DataFrame], 
                             torsionDir: DirectoryPath, 
                             torsionTag: str, 
                             config: dict) -> List[pd.DataFrame]:
    argsList = [(scanDir, scanDf, torsionDir, torsionTag, config) for scanDir, scanDf in zip(scanDirs, scanDfs)]
    singlePointDfs = []
    for args in argsList:
        singlePointDf = do_the_single_point_worker(args)
        if singlePointDf is not None:
            singlePointDfs.append(singlePointDf)
    return singlePointDfs

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Torsion Scanning", "TORSION_SCANS")
def scan_in_serial(torsionScanDir: DirectoryPath, 
                    conformerXyzs: List[FilePath], 
                    torsionIndexes: List[int], 
                    config: dict) -> Tuple[List[pd.DataFrame], List[DirectoryPath]]:
    argsList = [(conformerXyz, torsionScanDir,  torsionIndexes, config) for  conformerXyz in conformerXyzs]
    scanDfs = []
    scanDirs = []
    for args in argsList:
        scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir= do_the_twist_worker(args)
        if scanForwardsDf is not None and scanBackwardsDf is not None:
            scanDfs.append(scanForwardsDf)
            scanDfs.append(scanBackwardsDf)
            scanDirs.append(forwardsDir)
            scanDirs.append(backwardsDir)

    return scanDfs, scanDirs

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Scan Single Points", "TORSION_SCANS")
def single_points_in_parallel(scanDirs: List[DirectoryPath], 
                               scanDfs: List[pd.DataFrame], 
                               torsionTag: str, 
                               config: dict) -> List[pd.DataFrame]:
    argsList = [(scanDir, scanDf, torsionTag, config) for scanDir, scanDf in zip(scanDirs, scanDfs)]

    purpleText = "\033[35m"
    resetTextColor = "\033[0m"

    ## set up progress bar
    tqdmBarOptions = {
        "desc": f"{purpleText}Running Single-Points {resetTextColor}",
        "ascii": "-🗲→",    
        "colour": "yellow",
        "unit":  "scan",
        "ncols": 124,
        "dynamic_ncols": False
    }

    with WorkerPool(n_jobs = config["miscInfo"]["availableCpus"]) as pool:
        results = pool.map(do_the_single_point_worker,
                            make_single_arguments(argsList),
                              progress_bar=True,
                              iterable_len = len(argsList),
                              progress_bar_options=tqdmBarOptions)

    singlePointDfs = []
    for result in results:
        if result is not None:
            singlePointDfs.append(result)
    return singlePointDfs
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Torsion Scanning", "TORSION_SCANS")
def scan_in_parallel(torsionScanDir: DirectoryPath, 
                      conformerXyzs: List[FilePath], 
                      torsionIndexes: List[int], 
                      torsionTag: str, 
                      config: dict) -> Tuple[List[pd.DataFrame], List[DirectoryPath]]:
    greenText = "\033[32m"
    resetTextColor = "\033[0m"

    ## set up progress bar
    tqdmBarOptions = {
        "desc": f"{greenText}Scanning Torsion Angle{resetTextColor}",
        "ascii": "-🗲→",    
        "colour": "yellow",
        "unit":  "scan",
        "ncols": 124,
        "dynamic_ncols": False    }

    argsList = [(conformerXyz, torsionScanDir,  torsionIndexes, config) for  conformerXyz in conformerXyzs]

    ## save on cores 
    nCores = min(len(argsList), config["miscInfo"]["availableCpus"])

    with WorkerPool(n_jobs = nCores) as pool:
        results = pool.map(do_the_twist_worker,
                            make_single_arguments(argsList),
                              progress_bar=True,
                              iterable_len = len(argsList),
                              progress_bar_options=tqdmBarOptions)
    ## combine all dataframes into a big list
    scanDfs = []
    scanDirs = []
    for scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir in results:
        if scanForwardsDf is not None and scanBackwardsDf is not None:
            scanDfs.append(scanForwardsDf)
            scanDfs.append(scanBackwardsDf)
            scanDirs.append(forwardsDir)
            scanDirs.append(backwardsDir)

    return scanDfs, scanDirs

  

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_single_point_worker(args: tuple) -> Optional[pd.DataFrame]:
    """
    Worker function to perform single point calculations for a given scan.

    Args:
        args: A tuple containing (scanDir, scanDf, torsionTag, config).

    Returns:
        A pandas DataFrame with single point energies, or None if an error occurs.
    """
    scanDir, scanDf, torsionTag, config = args
    try:
        singlePointDf = Twisted_Monster.run_singlepoints_on_scans(scanDir=scanDir,
                                  scanDf = scanDf, 
                                  conformerId=torsionTag,
                                  config = config)
        return singlePointDf
    except Exception as e:
        # TODO: Log the exception e
        return None
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_twist_worker(args: tuple) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], Optional[DirectoryPath], Optional[DirectoryPath]]:
    """
    Worker function to perform torsion scanning for a given conformer.

    Args:
        args: A tuple containing (conformerXyz, torsionScanDir, torsionIndexes, config).

    Returns:
        A tuple containing (scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir),
        or (None, None, None, None) if an error occurs.
    """
    ## unpack args tuple
    conformerXyz, torsionScanDir,  torsionIndexes, config= args
    try:
        scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir  = do_the_twist(conformerXyz, torsionScanDir, torsionIndexes, config=config)
        return scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir
    except FileNotFoundError as e:
        # TODO: Log the exception e (this is fine TODO: make a custom error to look for - this is when a scan crashes)
        return None, None, None, None
    except Exception as e:
        ## for all other exceptions, raise them, this will be caught by the debugger in drFrankenstein.py
        raise(e)
    
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_twist(conformerXyz: FilePath,
                 torsionScanDir: DirectoryPath,
                 torsionIndexes: List[int],
                 config: dict) -> Tuple[pd.DataFrame, pd.DataFrame, DirectoryPath, DirectoryPath]:
    """
    Performs the forward and backward torsion scan for a given conformer.

    Args:
        conformerXyz: Path to the conformer XYZ file.
        torsionScanDir: Directory for the specific torsion scan.
        torsionIndexes: List of atom indexes defining the torsion.
        config: Configuration dictionary.

    Returns:
        A tuple containing (scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir).
    """
    conformerId = p.basename(conformerXyz).split(".")[0]

    conformerScanDir = p.join(torsionScanDir, f"scans_{conformerId}")
    optXyz = Twisted_Monster.run_optimization_step(conformerXyz, conformerScanDir, conformerId, config) # Updated call site
    ## get angle of torsion in this conformer
    initialTorsionAngle = Twisted_Assistant.measure_current_torsion_angle(optXyz, torsionIndexes) # type: ignore

    scanForwardsDf, forwardsXyz, forwardsDir = Twisted_Monster.run_forwards_scan_step(optXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

    scanBackwardsDf, backwardsDir = Twisted_Monster.run_backwards_scan_step(forwardsXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

    scanForwardsDf = Twisted_Assistant.process_energy_outputs(scanForwardsDf) # type: ignore
    scanBackwardsDf = Twisted_Assistant.process_energy_outputs(scanBackwardsDf) # type: ignore

    return  scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
