
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

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def twist_protocol(config):
    ## create an entry in runtimeInfo for twist
    """
    Main protocol for torsion scanning
    1. Create an entry in runtimeInfo for twist
    2. Set up directories for torsion scans
    3. Identify rotatable bonds to be scanned
    4. Choose torsions to be scanned
    5. Run torsion scanning for each rotatable bond
    6. Run Single points on scans
    6. Save torsion scan data

    Args:
        config (dict): the drFrankenstein config containing all run information

    Returns:
        config (dict): updated config
    """
    config["runtimeInfo"]["madeByTwisting"] = {}
    config["runtimeInfo"]["madeByTwisting"]["torsionDirs"] = []
    config["runtimeInfo"]["madeByTwisting"]["torsionTags"] = []
    config["runtimeInfo"]["madeByTwisting"]["finalScanEnergies"] = {}
    config["runtimeInfo"]["madeByTwisting"]["scanEnergyData"] = {}

    config = Twisted_Assistant.set_up_directories(config)

    ## identify rotatable bonds to be scanned
    config = Twisted_Assistant.identify_rotatable_bonds(config, mode = config["parameterFittingInfo"]["forceField"])

    config = Twisted_Assistant.choose_torsions_to_scan(config)

    torsionsToScan = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]

    nRotatableBonds = len(torsionsToScan)
    config["runtimeInfo"]["madeByTwisting"]["nRotatableBonds"] = nRotatableBonds
    for torsionIndex, (torsionTag, torsionData) in enumerate(torsionsToScan.items()):
        config = run_torsion_scanning(torsionTag, torsionData, torsionIndex, nRotatableBonds, config)
    config["checkpointInfo"]["scanningComplete"] = True

    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_torsion_scanning(torsionTag: str,
                          torsionData: dict[str:tuple[str], str:tuple[str], str:tuple[str]],
                            torsionIndex: int,
                              nRotatableBonds: int,
                                config: dict,
                                  debug: bool = False) -> dict:

    """
    Runs torsion scanning for a given torsion.

    This function is the main entry point for torsion scanning.
    It creates a top level directory for the torsion, runs torsion scanning in either serial or parallel,
    processes the scan data, and plots the results.

    Args:
        torsionTag (str): torsion tag
        torsionData (dict[str:tuple[str], str:tuple[str], str:tuple[str]]): torsion data
        torsionIndex (int): torsion index
        nRotatableBonds (int): number of rotatable bonds
        config (dict): the drFrankenstein config containing all run information
        debug (bool, optional): whether to run in debug mode. Defaults to False.

    Returns:
        dict: updated config
    """
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
@Timer.time_function("Torsion Scanning (serial)", "TORSION_SCANS")
def scan_in_serial(torsionScanDir, conformerXyzs, torsionIndexes, config) -> Tuple[pd.DataFrame, DirectoryPath]:
    """
    Runs torsion scanning for a given torsion in serial.

    This function runs torsion scanning in serial by calling do_the_twist_worker on each conformer.
    It appends the output of each scan to a list of dataframes and directories.

    Args:
        torsionScanDir (DirectoryPath): directory to store torsion scan data
        conformerXyzs (List[FilePath]): list of conformer xyz files
        torsionIndexes (List[int]): list of torsion indexes
        config (dict): the drFrankenstein config containing all run information

    Returns:
        Tuple[pd.DataFrame, DirectoryPath]: tuple of dataframes and directories
    """
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
@Timer.time_function("Torsion Scanning", "TORSION_SCANS")
def scan_in_parallel(torsionScanDir, conformerXyzs, torsionIndexes, torsionTag, config) -> Tuple[pd.DataFrame, DirectoryPath]:
    """
    Executes torsion scanning in parallel for a given torsion.

    This function utilizes multiprocessing to perform torsion scanning in parallel,
    distributing the workload across available CPU cores. The progress of the scanning
    is displayed using a progress bar.

    Args:
        torsionScanDir (DirectoryPath): Directory to store torsion scan data.
        conformerXyzs (List[FilePath]): List of conformer xyz files.
        torsionIndexes (List[int]): List of torsion indexes.
        torsionTag (str): Torsion tag.
        config (dict): The drFrankenstein config containing all run information.

    Returns:
        Tuple[List[pd.DataFrame], List[DirectoryPath]]: A tuple containing lists of
        dataframes with scan results and corresponding directories.
    """

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
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Scan Single Points (serial)", "TORSION_SCANS")
def single_points_in_serial(scanDirs, scanDfs, torsionDir, torsionTag, config):
    """
    Runs single-point energy calculations for a given torsion in serial.

    This function runs single-point energy calculations in serial by calling do_the_single_point_worker on each conformer.
    It appends the output of each scan to a list of dataframes.

    Args:
        scanDirs (List[DirectoryPath]): list of directories to store torsion scan data
        scanDfs (List[pd.DataFrame]): list of dataframes with scan results
        torsionDir (DirectoryPath): directory to store torsion scan data
        torsionTag (str): torsion tag
        config (dict): the drFrankenstein config containing all run information

    Returns:
        List[pd.DataFrame]: list of dataframes with single-point energy results
    """
    argsList = [(scanDir, scanDf, torsionTag, config) for scanDir, scanDf in zip(scanDirs, scanDfs)]
    singlePointDfs = []
    for args in argsList:
        singlePointDf = do_the_single_point_worker(args)
        if singlePointDf is not None:
            singlePointDfs.append(singlePointDf)
    return singlePointDfs
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Scan Single Points", "TORSION_SCANS")
def single_points_in_parallel(scanDirs, scanDfs, torsionTag, config):
    """
    Executes single-point energy calculations for torsions in parallel.

    This function utilizes multiprocessing to perform single-point energy
    calculations in parallel, distributing the workload across available
    CPU cores. The progress of the calculations is displayed using a
    progress bar.

    Args:
        scanDirs (List[DirectoryPath]): List of directories containing
        torsion scan data.
        scanDfs (List[pd.DataFrame]): List of dataframes with scan results.
        torsionTag (str): Identifier for the torsion being processed.
        config (dict): Configuration dictionary containing run information.

    Returns:
        List[pd.DataFrame]: A list of dataframes containing single-point 
        energy results for each torsion.
    """

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

  

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_single_point_worker(args):
    """
    Executes single-point energy calculations for a given torsion.

    This function runs single-point energy calculations using the provided scan directory,
    dataframe, torsion tag, and configuration. It returns the resulting dataframe
    containing the energy results. If an exception occurs during the execution,
    the function returns None.

    Args:
        args (tuple): A tuple containing the following:
            - scanDir (DirectoryPath): Directory to store scan data.
            - scanDf (pd.DataFrame): DataFrame with scan results.
            - torsionTag (str): Identifier for the torsion being processed.
            - config (dict): Configuration dictionary containing run information.

    Returns:
        pd.DataFrame or None: DataFrame containing single-point energy results,
        or None if an error occurs.
    """

    scanDir, scanDf, torsionTag, config = args
    try:
        singlePointDf = Twisted_Monster.run_singlepoints_on_scans(scanDir=scanDir,
                                  scanDf = scanDf, 
                                  conformerId=torsionTag,
                                  config = config)
        return singlePointDf
    except Exception as e:
        ...
        return None
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_twist_worker(args):
    """
    Executes the torsion scan for a given conformer.

    This function runs the torsion scan by calling do_the_twist with the provided
    arguments. It handles exceptions by returning None values if a FileNotFoundError
    occurs, indicating a scan crash, or raising other exceptions to be caught by the
    debugger in drFrankenstein.py.

    Args:
        args (tuple): A tuple containing the following:
            - conformerXyz (FilePath): Path to the conformer xyz file.
            - torsionScanDir (DirectoryPath): Directory to store torsion scan data.
            - torsionIndexes (List[int]): List of torsion indexes.
            - config (dict): Configuration dictionary containing run information.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, DirectoryPath, DirectoryPath]:
        - scanForwardsDf (pd.DataFrame): DataFrame with forward scan results.
        - scanBackwardsDf (pd.DataFrame): DataFrame with backward scan results.
        - forwardsDir (DirectoryPath): Directory for forward scan data.
        - backwardsDir (DirectoryPath): Directory for backward scan data.
        Returns None values if a scan crashes due to FileNotFoundError.
    """
    ## unpack args tuple
    conformerXyz, torsionScanDir,  torsionIndexes, config= args
    try:
        scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir  = do_the_twist(conformerXyz, torsionScanDir, torsionIndexes, config=config)
        return scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir
    except FileNotFoundError as e:
        ## this is fine TODO: make a custom error to look for - this is when a scan crashes
        return None, None, None, None
    except Exception as e:
        ## for all other exceptions, raise them, this will be caught by the debugger in drFrankenstein.py
        raise(e)
    
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_twist(conformerXyz: FilePath,
                 torsionScanDir: DirectoryPath,
                 torsionIndexes: List[int],
                 config: dict) -> Tuple[pd.DataFrame]:

    """
    Runs torsion scanning for a given conformer.

    Args:
        conformerXyz (FilePath): Path to the conformer xyz file.
        torsionScanDir (DirectoryPath): Directory to store torsion scan data.
        torsionIndexes (List[int]): List of torsion indexes.
        config (dict): Configuration dictionary containing run information.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, DirectoryPath, DirectoryPath]:
        - scanForwardsDf (pd.DataFrame): DataFrame with forward scan results.
        - scanBackwardsDf (pd.DataFrame): DataFrame with backward scan results.
        - forwardsDir (DirectoryPath): Directory for forward scan data.
        - backwardsDir (DirectoryPath): Directory for backward scan data.
        Returns None values if a scan crashes due to FileNotFoundError.
    """
    conformerId = p.basename(conformerXyz).split(".")[0]

    conformerScanDir = p.join(torsionScanDir, f"scans_{conformerId}")
    optXyz = Twisted_Monster.run_optimisation_step(conformerXyz, conformerScanDir, conformerId, config)
    ## get angle of torsion in this conformer
    initialTorsionAngle = Twisted_Assistant.measure_current_torsion_angle(optXyz, torsionIndexes)

    scanForwardsDf, forwardsXyz, forwardsDir = Twisted_Monster.run_forwards_scan_step(optXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

    if Twisted_Assistant.have_scans_exploded(forwardsDir, config):
        return None, None, None, None

    scanBackwardsDf, backwardsDir = Twisted_Monster.run_backwards_scan_step(forwardsXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

    if Twisted_Assistant.have_scans_exploded(backwardsDir, config):
        return None, None, None, None

    scanForwardsDf = Twisted_Assistant.process_energy_outputs(scanForwardsDf)
    scanBackwardsDf = Twisted_Assistant.process_energy_outputs(scanBackwardsDf)

    return  scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
