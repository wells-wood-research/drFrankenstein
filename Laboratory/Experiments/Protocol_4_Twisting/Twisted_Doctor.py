import os
from os import path as p
import pandas as pd
from mpire import WorkerPool
from mpire.utils import make_single_arguments

class FilePath:
    pass

class DirectoryPath:
    pass
from typing import Dict, List, Tuple
from . import Twisted_Assistant
from . import Twisted_Monster
from . import Twisted_Plotter
from OperatingTools import drSplash
from OperatingTools import Timer

def twist_protocol(config):
    config['runtimeInfo']['madeByTwisting'] = {}
    config['runtimeInfo']['madeByTwisting']['torsionDirs'] = []
    config['runtimeInfo']['madeByTwisting']['torsionTags'] = []
    config['runtimeInfo']['madeByTwisting']['finalScanEnergies'] = {}
    config['runtimeInfo']['madeByTwisting']['scanEnergyData'] = {}
    config = Twisted_Assistant.set_up_directories(config)
    config = Twisted_Assistant.identify_rotatable_bonds(config, mode=config['parameterFittingInfo']['forceField'])
    if config['torsionScanInfo']['preserveBackboneTorsions']:
        config = Twisted_Assistant.exclude_backbone_torsions(config)
    rotatableDihedrals = config['runtimeInfo']['madeByTwisting']['rotatableDihedrals']
    nRotatableBonds = len(rotatableDihedrals)
    config['runtimeInfo']['madeByTwisting']['nRotatableBonds'] = nRotatableBonds
    for (torsionIndex, (torsionTag, torsionData)) in enumerate(rotatableDihedrals.items()):
        config = run_torsion_scanning(torsionTag, torsionData, torsionIndex, nRotatableBonds, config)
    config['checkpointInfo']['scanningComplete'] = True
    return config

def run_torsion_scanning(torsionTag: str, torsionData: dict[str:tuple[str], str:tuple[str], str:tuple[str]], torsionIndex: int, nRotatableBonds: int, config: dict, debug: bool=False) -> dict:
    torsionDir = p.join(config['runtimeInfo']['madeByTwisting']['torsionDir'], f'torsion_{torsionTag}')
    os.makedirs(torsionDir, exist_ok=True)
    config['runtimeInfo']['madeByTwisting']['torsionDirs'].append(torsionDir)
    config['runtimeInfo']['madeByTwisting']['torsionTags'].append(torsionTag)
    conformerXyzs = Twisted_Assistant.get_conformer_xyzs(config)
    drSplash.show_torsion_being_scanned(torsionTag, torsionIndex, nRotatableBonds)
    if debug:
        (scanDfs, scanDirs) = scan_in_serial(torsionDir, conformerXyzs, torsionData['ATOM_INDEXES'], config=config)
        if config['torsionScanInfo']['singlePointMethod'] is None:
            singlePointDfs = None
        else:
            singlePointDfs = single_points_in_serial(scanDirs, scanDfs, torsionDir, torsionTag, config=config)
    else:
        (scanDfs, scanDirs) = scan_in_parallel(torsionDir, conformerXyzs, torsionData['ATOM_INDEXES'], torsionTag, config=config)
        if config['torsionScanInfo']['scanSinglePointsOn'] is None or config['torsionScanInfo']['singlePointMethod'] == None:
            singlePointDfs = None
        else:
            singlePointDfs = single_points_in_parallel(scanDirs, scanDfs, torsionTag, config=config)
    (scanEnergiesCsv, scanAveragesDf) = Twisted_Assistant.process_scan_data(scanDfs, torsionDir, torsionTag)
    if config['torsionScanInfo']['singlePointMethod'] is None:
        config['runtimeInfo']['madeByTwisting']['finalScanEnergies'][torsionTag] = scanEnergiesCsv
        singlePointAveragesDf = None
        config = Twisted_Assistant.gather_scan_data(scanAveragesDf, torsionTag, config)
    else:
        (singlePointEnergiesCsv, singlePointAveragesDf) = Twisted_Assistant.process_scan_data(singlePointDfs, torsionDir, torsionTag)
        config['runtimeInfo']['madeByTwisting']['finalScanEnergies'][torsionTag] = singlePointEnergiesCsv
        config = Twisted_Assistant.gather_scan_data(singlePointAveragesDf, torsionTag, config)
    config = Twisted_Plotter.twist_plotting_protocol(scanDfs, scanAveragesDf, singlePointDfs, singlePointAveragesDf, torsionDir, torsionTag, config)
    return config

@Timer.time_function('Torsion Scanning (serial)', 'TORSION_SCANS')
def scan_in_serial(torsionScanDir, conformerXyzs, torsionIndexes, config) -> Tuple[pd.DataFrame, DirectoryPath]:
    argsList = [(conformerXyz, torsionScanDir, torsionIndexes, config) for conformerXyz in conformerXyzs]
    scanDfs = []
    scanDirs = []
    for args in argsList:
        (scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir) = do_the_twist_worker(args)
        if scanForwardsDf is not None and scanBackwardsDf is not None:
            scanDfs.append(scanForwardsDf)
            scanDfs.append(scanBackwardsDf)
            scanDirs.append(forwardsDir)
            scanDirs.append(backwardsDir)
    return (scanDfs, scanDirs)

@Timer.time_function('Torsion Scanning', 'TORSION_SCANS')
def scan_in_parallel(torsionScanDir, conformerXyzs, torsionIndexes, torsionTag, config) -> Tuple[pd.DataFrame, DirectoryPath]:
    greenText = '\x1b[32m'
    resetTextColor = '\x1b[0m'
    tqdmBarOptions = {'desc': f'{greenText}Scanning Torsion Angle{resetTextColor}', 'ascii': '-🗲→', 'colour': 'yellow', 'unit': 'scan', 'ncols': 124, 'dynamic_ncols': False}
    argsList = [(conformerXyz, torsionScanDir, torsionIndexes, config) for conformerXyz in conformerXyzs]
    nCores = min(len(argsList), config['miscInfo']['availableCpus'])
    with WorkerPool(n_jobs=nCores) as pool:
        results = pool.map(do_the_twist_worker, make_single_arguments(argsList), progress_bar=True, iterable_len=len(argsList), progress_bar_options=tqdmBarOptions)
    scanDfs = []
    scanDirs = []
    for (scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir) in results:
        if scanForwardsDf is not None and scanBackwardsDf is not None:
            scanDfs.append(scanForwardsDf)
            scanDfs.append(scanBackwardsDf)
            scanDirs.append(forwardsDir)
            scanDirs.append(backwardsDir)
    return (scanDfs, scanDirs)

@Timer.time_function('Scan Single Points (serial)', 'TORSION_SCANS')
def single_points_in_serial(scanDirs, scanDfs, torsionDir, torsionTag, config):
    argsList = [(scanDir, scanDf, torsionDir, torsionTag, config) for (scanDir, scanDf) in zip(scanDirs, scanDfs)]
    singlePointDfs = []
    for args in argsList:
        singlePointDf = do_the_single_point_worker(args)
        if singlePointDf is not None:
            singlePointDfs.append(singlePointDf)
    return singlePointDfs

@Timer.time_function('Scan Single Points', 'TORSION_SCANS')
def single_points_in_parallel(scanDirs, scanDfs, torsionTag, config):
    argsList = [(scanDir, scanDf, torsionTag, config) for (scanDir, scanDf) in zip(scanDirs, scanDfs)]
    purpleText = '\x1b[35m'
    resetTextColor = '\x1b[0m'
    tqdmBarOptions = {'desc': f'{purpleText}Running Single-Points {resetTextColor}', 'ascii': '-🗲→', 'colour': 'yellow', 'unit': 'scan', 'ncols': 124, 'dynamic_ncols': False}
    with WorkerPool(n_jobs=config['miscInfo']['availableCpus']) as pool:
        results = pool.map(do_the_single_point_worker, make_single_arguments(argsList), progress_bar=True, iterable_len=len(argsList), progress_bar_options=tqdmBarOptions)
    singlePointDfs = []
    for result in results:
        if result is not None:
            singlePointDfs.append(result)
    return singlePointDfs

def do_the_single_point_worker(args):
    (scanDir, scanDf, torsionTag, config) = args
    try:
        singlePointDf = Twisted_Monster.run_singlepoints_on_scans(scanDir=scanDir, scanDf=scanDf, conformerId=torsionTag, config=config)
        return singlePointDf
    except Exception as e:
        ...
        return None

def do_the_twist_worker(args):
    (conformerXyz, torsionScanDir, torsionIndexes, config) = args
    try:
        (scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir) = do_the_twist(conformerXyz, torsionScanDir, torsionIndexes, config=config)
        return (scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir)
    except FileNotFoundError as e:
        return (None, None, None, None)
    except Exception as e:
        raise e

def do_the_twist(conformerXyz: FilePath, torsionScanDir: DirectoryPath, torsionIndexes: List[int], config: dict) -> Tuple[pd.DataFrame]:
    conformerId = p.basename(conformerXyz).split('.')[0]
    conformerScanDir = p.join(torsionScanDir, f'scans_{conformerId}')
    optXyz = Twisted_Monster.run_optimisation_step(conformerXyz, conformerScanDir, conformerId, config)
    initialTorsionAngle = Twisted_Assistant.measure_current_torsion_angle(optXyz, torsionIndexes)
    (scanForwardsDf, forwardsXyz, forwardsDir) = Twisted_Monster.run_forwards_scan_step(optXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)
    (scanBackwardsDf, backwardsDir) = Twisted_Monster.run_backwards_scan_step(forwardsXyz, initialTorsionAngle, torsionIndexes, conformerScanDir, conformerId, config)
    scanForwardsDf = Twisted_Assistant.process_energy_outputs(scanForwardsDf)
    scanBackwardsDf = Twisted_Assistant.process_energy_outputs(scanBackwardsDf)
    return (scanForwardsDf, scanBackwardsDf, forwardsDir, backwardsDir)