"""
This script contains:
 - twist_protocol (the main function for torsion scanning)
 - run_torsion_scanning (scanning for each rotatable bond)
 - scan_in_serial (non-multiprocessing, use for debug)
 - scan_in_parallel (multiprocessing with nice loading bars [IMPORTANT])
 - do_the_twist_worker(small helper for do_the_twist, has as try/except and unpacks args)

 and helper function
"""
import os
from os import path as p
import sys
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


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def twist_protocol(config):
    config = Twisted_Assistant.set_up_directories(config)

    ## identify rotatable bonds to be scanned
    rotatableBonds: Dict[List[Tuple[int,int,int,int]]] = Twisted_Assistant.identify_rotatable_bonds(config["moleculeInfo"]["cappedPdb"])

    ## loop over torsions that need to be scanned
    config["pathInfo"]["torsionDirs"] = []
    config["torsionScanInfo"]["torsionTags"] = []
    config["torsionScanInfo"]["finalScanEnergies"] = {}

    nRotatableBonds = len(rotatableBonds)
    for torsionIndex, rotatableBond in enumerate(rotatableBonds):
        config = run_torsion_scanning(rotatableBond, torsionIndex, nRotatableBonds, config)
    config["checkpointInfo"]["scanningComplete"] = True
    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_torsion_scanning(rotatableBond, torsionIndex, nRotatableBonds, config, debug=False) -> dict:

    ## make a top level dir for this torsion
    torsionTag: str = "-".join(map(str, rotatableBond["atoms"])) 
    torsionDir = p.join(config["pathInfo"]["torsionTopDir"], f"torsion_{torsionTag}" )
    os.makedirs(torsionDir, exist_ok=True)
    ## add to config
    config["pathInfo"]["torsionDirs"].append(torsionDir)
    config["torsionScanInfo"]["torsionTags"].append(torsionTag)

    conformerXyzs = Twisted_Assistant.get_conformer_xyzs(config)
    drSplash.show_torsion_being_scanned(torsionTag, torsionIndex, nRotatableBonds)
    if debug:
        ## run in serial
        scanDfs, scanDirs = scan_in_serial(torsionDir, conformerXyzs, rotatableBond["indices"], config)
        if config["torsionScanInfo"]["scanSinglePointsOn"] is None:
            singlePointDfs = None
        else:
            singlePointDfs = single_points_in_serial(scanDirs, scanDfs, torsionDir, config, torsionTag)
    else:
        # run torsion scans in parallel
        scanDfs, scanDirs = scan_in_parallel(torsionDir, conformerXyzs, rotatableBond["indices"], torsionTag, config)
        if config["torsionScanInfo"]["scanSinglePointsOn"] is None:
            singlePointDfs = None
        else:
            singlePointDfs = single_points_in_parallel(scanDirs, scanDfs, config, torsionTag)
    ## Merge scan data, calculate averages, rolling averages and mean average errors
    scanEnergiesCsv, scanAveragesDf  = Twisted_Assistant.process_scan_data(scanDfs, torsionDir, torsionTag)
    if  config["torsionScanInfo"]["singlePointMethod"] is None:
        config["torsionScanInfo"]["finalScanEnergies"][torsionTag] = scanEnergiesCsv
        singlePointAveragesDf = None
    else:
        singlePointEnergiesCsv, singlePointAveragesDf = Twisted_Assistant.process_scan_data(singlePointDfs, torsionDir, torsionTag)
        config["torsionScanInfo"]["finalScanEnergies"][torsionTag] = singlePointEnergiesCsv   


    ## Plotting
    Twisted_Plotter.twist_plotting_protocol(scanDfs, scanAveragesDf, singlePointDfs, singlePointAveragesDf, torsionDir, torsionTag, config)
    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def single_points_in_serial(scanDirs, scanDfs, torsionDir, config, torsionTag):
    argsList = [(scanDir, scanDf, torsionDir, torsionTag, config) for scanDir, scanDf in zip(scanDirs, scanDfs)]
    singlePointDfs = []
    for args in argsList:
        singlePointDf = do_the_single_point_worker(args)
        if singlePointDf is not None:
            singlePointDfs.append(singlePointDf)
    return singlePointDfs

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def scan_in_serial(torsionScanDir, conformerXyzs, torsionIndexes, config) -> Tuple[pd.DataFrame, DirectoryPath]:
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

def single_points_in_parallel(scanDirs, scanDfs, config, torsionTag):
    argsList = [(scanDir, scanDf, torsionTag, config) for scanDir, scanDf in zip(scanDirs, scanDfs)]
    tqdmBarOptions = {
        "desc": f"\033[32mRunning Single-Points \033[0m",
        "ascii": "-▒",  
        "colour": "magenta",
        "unit":  "scan",
        "dynamic_ncols": True
    }
    with WorkerPool(n_jobs = config["hardwareInfo"]["nCores"]) as pool:
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
def scan_in_parallel(torsionScanDir, conformerXyzs, torsionIndexes, torsionTag, config) -> Tuple[pd.DataFrame, DirectoryPath]:
    tqdmBarOptions = {
        "desc": f"\033[32mScanning Torsion Angle\033[0m",
        "ascii": "-▒",  
        "colour": "green",
        "unit":  "scan",
        "dynamic_ncols": True
    }


    argsList = [(conformerXyz, torsionScanDir,  torsionIndexes, config) for  conformerXyz in conformerXyzs]

    ## save on cores 
    nCores = min(len(argsList), config["hardwareInfo"]["nCores"])

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
def do_the_single_point_worker(args):
    scanDir, scanDf, torsionTag, config = args
    try:
        singlePointDf = Twisted_Monster.run_singlepoints_on_scans(scanDir=scanDir,
                                  scanDf = scanDf, 
                                  conformerId=torsionTag,
                                  config = config)
        return singlePointDf
    except Exception as e:
        print(e)
        return None
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_twist_worker(args):
    ## unpack args tuple
    conformerXyz, torsionScanDir,  torsionIndexes, config= args
    try:
        scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir  = do_the_twist(conformerXyz, torsionScanDir, torsionIndexes, config)
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

    conformerId = p.basename(conformerXyz).split(".")[0]

    conformerScanDir = p.join(torsionScanDir, f"scans_{conformerId}")
    optXyz = Twisted_Monster.run_optimisation_step(conformerXyz, conformerScanDir, conformerId, config)
    ## get angle of torsion in this conformer
    initialTorsionAngle = Twisted_Assistant.measure_current_torsion_angle(optXyz, torsionIndexes)

    scanForwardsDf, forwardsXyz, forwardsDir = Twisted_Monster.run_forwards_scan_step(optXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

    scanBackwardsDf, backwardsDir = Twisted_Monster.run_backwards_scan_step(forwardsXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

    scanForwardsDf = Twisted_Assistant.process_energy_outputs(scanForwardsDf)
    scanBackwardsDf = Twisted_Assistant.process_energy_outputs(scanBackwardsDf)

    return  scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
