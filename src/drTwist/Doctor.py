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

## ADD SRC TO PATH ##
currentFilePath: FilePath = os.path.abspath(__file__)
currentDir: DirectoryPath = os.path.dirname(currentFilePath)
srcDir: DirectoryPath = os.path.dirname(currentDir)
sys.path.append(srcDir)

from drTwist import Assistant, Monster, Plotter


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def twist_protocol(config):

    config = Assistant.set_up_directories(config)

    ## identify rotatable bonds to be scanned
    rotatableBonds: Dict[List[Tuple[int,int,int,int]]] = Assistant.identify_rotatable_bonds(config["moleculeInfo"]["cappedPdb"])

    ## loop over torsions that need to be scanned
    config["pathInfo"]["torsionDirs"] = []
    config["torsionScanInfo"]["torsionTags"] = []
    config["torsionScanInfo"]["finalScanEnergies"] = {}
    for rotatableBond in rotatableBonds:
        config = run_torsion_scanning(rotatableBond, config)

    config["checkpointInfo"]["scanningComplete"] = True
    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_torsion_scanning(rotatableBond, config, debug=False) -> dict:
    ## make a top level dir for this torsion
    torsionTag: str = "-".join(map(str, rotatableBond["atoms"])) 
    torsionDir = p.join(config["pathInfo"]["torsionTopDir"], f"torsion_{torsionTag}" )
    os.makedirs(torsionDir, exist_ok=True)
    ## add to config
    config["pathInfo"]["torsionDirs"].append(torsionDir)
    config["torsionScanInfo"]["torsionTags"].append(torsionTag)

 
    ## get conformer XYZ files
    conformerXyzs = config["pathInfo"]["conformerXyzs"]

    if debug:
        ## run in serial
        scanDfs, singlePointDfs = scan_in_serial(torsionDir, conformerXyzs, rotatableBond["indices"], config)
    else:
        # run torsion scans in paralell
        scanDfs, singlePointDfs = scan_in_parallel(torsionDir, conformerXyzs, rotatableBond["indices"], torsionTag, config)

    ## Merge scan data, calculate averages, rolling averages and mean average errors
    scanEnergiesCsv, scanAveragesDf  = Assistant.process_scan_data(scanDfs, torsionDir, torsionTag)
    if  config["torsionScanInfo"]["singlePointMethod"] is None:
        config["torsionScanInfo"]["finalScanEnergies"][torsionTag] = scanEnergiesCsv
        singlePointAveragesDf = None
    else:
        singlePointEnergiesCsv, singlePointAveragesDf = Assistant.process_scan_data(singlePointDfs, torsionDir, torsionTag)
        config["torsionScanInfo"]["finalScanEnergies"][torsionTag] = singlePointEnergiesCsv   


    ## Plotting
    Plotter.twist_plotting_protocol(scanDfs, scanAveragesDf, singlePointDfs, singlePointAveragesDf, torsionDir, torsionTag, config)
    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def scan_in_serial(torsionScanDir, conformerXyzs, torsionIndexes, config) -> List[pd.DataFrame]:
    argsList = [(conformerXyz, torsionScanDir,  torsionIndexes, config) for  conformerXyz in conformerXyzs]
    scanDfs = []
    singlePointDfs = []
    for args in argsList:
        scanForwardsDf, scanBackwardsDf, singlePointForwardsDf, singlePointBackwardsDf = do_the_twist_worker(args)
        if scanForwardsDf is not None and scanBackwardsDf is not None:
            scanDfs.append(scanForwardsDf)
            scanDfs.append(scanBackwardsDf)
        if singlePointForwardsDf is not None and singlePointBackwardsDf is not None:
            singlePointDfs.append(singlePointForwardsDf)
            singlePointDfs.append(singlePointBackwardsDf)

    return scanDfs, singlePointDfs
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def scan_in_parallel(torsionScanDir, conformerXyzs, torsionIndexes, torsionTag, config) -> List[pd.DataFrame]:
    tqdmBarOptions = {
        "desc": f"\033[32mScanning Torsion {torsionTag}\033[0m",
        "ascii": "-ϟ",  
        "colour": "yellow",
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
    scanDfs = []
    singlePointDfs = []
    for scanForwardsDf, scanBackwardsDf, singlePointForwardsDf, singlePointBackwardsDf  in results:
        if scanForwardsDf is not None and scanBackwardsDf is not None:
            scanDfs.append(scanForwardsDf)
            scanDfs.append(scanBackwardsDf)
        if singlePointForwardsDf is not None and singlePointBackwardsDf is not None:
            singlePointDfs.append(singlePointForwardsDf)
            singlePointDfs.append(singlePointBackwardsDf)
    return scanDfs, singlePointDfs
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_twist_worker(args):
    ## unpack args tuple
    conformerXyz, torsionScanDir,  torsionIndexes, config= args
    try:
        scanForwardsDf, scanBackwardsDf, singlePointForwardsDf, singlePointBackwardsDf  = do_the_twist(conformerXyz, torsionScanDir, torsionIndexes, config)
        return scanForwardsDf, scanBackwardsDf, singlePointForwardsDf, singlePointBackwardsDf 
    except FileNotFoundError as e:
        ## this is fine
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
    optXyz = Monster.run_optimisation_step(conformerXyz, conformerScanDir, conformerId, config)
    ## get angle of torsion in this conformer
    initialTorsionAngle = Assistant.measure_current_torsion_angle(optXyz, torsionIndexes)

    scanForwardsDf, forwardsXyz, forwardsDir = Monster.run_forwards_scan_step(optXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

    scanBackwardsDf, backwardsDir = Monster.run_backwards_scan_step(forwardsXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)

    ## if no QM method has been specified, return data as-is
    if config["torsionScanInfo"]["singlePointMethod"] is None:
        scanForwardsDf = Assistant.process_energy_outputs(scanForwardsDf)
        scanBackwardsDf = Assistant.process_energy_outputs(scanBackwardsDf)
        return scanForwardsDf, scanBackwardsDf, None, None
    
    ## otherwise apply the single-point protocol
    singlePointForwardsDf = Monster.run_singlepoints_on_scans(scanDir=forwardsDir,
                              scanDf = scanForwardsDf, 
                              outDir=conformerScanDir,
                              conformerId=conformerId,
                              config = config,
                              tag = "forwards")
    

    singlePointBackwardsDf = Monster.run_singlepoints_on_scans(scanDir=backwardsDir,
                              scanDf = scanBackwardsDf,
                              outDir=conformerScanDir,
                              conformerId=conformerId,
                              config = config,
                              tag = "backwards")
    
    scanForwardsDf = Assistant.process_energy_outputs(scanForwardsDf)
    scanBackwardsDf = Assistant.process_energy_outputs(scanBackwardsDf)

    singlePointForwardsDf = Assistant.process_energy_outputs(singlePointForwardsDf)
    singlePointBackwardsDf = Assistant.process_energy_outputs(singlePointBackwardsDf)

    return scanForwardsDf, scanBackwardsDf, singlePointForwardsDf, singlePointBackwardsDf
