## BASIC LIBRARIES ##
import os
from os import path as p
import pandas as pd
import numpy as np
import re
from shutil import rmtree
from scipy.signal import savgol_filter
from tqdm import tqdm

## MULTIPROCESSION LIBRARIES ##
from concurrent.futures import ProcessPoolExecutor
 ## MULTIPROCESSING AND LOADING BAR LIBRARIES ##
from mpire import WorkerPool
from mpire.utils import make_single_arguments

import multiprocessing

## OPENMM LIBRARIES
import openmm.app as app
import openmm as openmm
import  openmm.unit  as unit

## drFRANKENSTEIN LIBRARIES ##
from .. import Stitching_Assistant
from . import AMBER_helper_functions
from OperatingTools import file_parsers

from OperatingTools import drSplash
## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass


def get_MM_total_energies(config, torsionTag, debug=True):
    drSplash.show_getting_mm_total(torsionTag)

    mmTotalDir: DirectoryPath = config["runtimeInfo"]["madeByStitching"]["mmTotalCalculationDir"]
    cappedPdb: FilePath = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    completedTorsionScanDirs: list = Stitching_Assistant.get_completed_torsion_scan_dirs(config, torsionTag)
    
    moleculeFrcmod: FilePath = config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"]
    moleculePrmtop: FilePath = config["runtimeInfo"]["madeByStitching"]["moleculePrmtop"]

    ## read calculated partial charges
    chargesCsv: FilePath = config["runtimeInfo"]["madeByCharges"]["chargesCsv"]
    chargesDf: pd.DataFrame = pd.read_csv(chargesCsv, index_col="Unnamed: 0")
    ## round charge to 4 decimal places
    chargesDf["Charge"] = chargesDf["Charge"].round(4)

    torsionTotalDir = p.join(mmTotalDir, f"{torsionTag}")
    os.makedirs(torsionTotalDir, exist_ok=True)

    torsionFittingDir = p.join(torsionTotalDir, "fitting_data")
    os.makedirs(torsionFittingDir, exist_ok=True)

    if debug:
        singlePointEnergyDfs = run_serial(completedTorsionScanDirs, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, moleculePrmtop)
    else:
        singlePointEnergyDfs = run_parallel(completedTorsionScanDirs, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, moleculePrmtop, config, torsionTag)

    ## sort out data
    mergedEnergyDf = Stitching_Assistant.merge_energy_dfs(singlePointEnergyDfs)
    mergedEnergyDf["Mean_Energy"] = mergedEnergyDf.drop(columns="Angle").mean(axis=1)

    mergedEnergycsv = p.join(torsionFittingDir, "raw_energy_data.csv")
    mergedEnergyDf.to_csv(mergedEnergycsv, index=False)
        
    finalScanEnergiesDf = pd.DataFrame()
    finalScanEnergiesDf["Angle"] = mergedEnergyDf["Angle"]
    finalScanEnergiesDf[torsionTag] = mergedEnergyDf["Mean_Energy"]
    
    finalScanEnergiesDf[torsionTag] = finalScanEnergiesDf[torsionTag] - finalScanEnergiesDf[torsionTag].min()

    finalScanEnergiesDf["smoothedEnergy"] = savgol_filter(finalScanEnergiesDf[torsionTag], window_length=5, polyorder=2)
    finalScanEnergiesDf.to_csv(p.join(torsionFittingDir, "final_scan_energies.csv"), index=False)

    return  finalScanEnergiesDf["smoothedEnergy"].to_numpy()

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_serial(scanDirs, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, moleculePrmtop):
    argsList = [(scanIndex, scanDir, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, moleculePrmtop) for  scanIndex, scanDir in enumerate(scanDirs)]
    singlePointEnergyDfs = []
    for args in argsList:
        singlePointEnergyDf = get_singlepoint_energies_for_torsion_scan(*args)
        singlePointEnergyDfs.append(singlePointEnergyDf)
    return singlePointEnergyDfs


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_parallel(scanDirs, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, moleculePrmtop, config, torsionTag):
    argsList = [(scanIndex, scanDir, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, moleculePrmtop) for scanIndex, scanDir in enumerate(scanDirs)]

    # Use multiprocessing Pool without tqdm progress bar
    with multiprocessing.Pool() as pool:
        singlePointEnergyDfs = list(pool.imap(single_point_worker, argsList))
    
    return singlePointEnergyDfs

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def single_point_worker(args):
    scanIndex, scanDir, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, moleculePrmtop = args
    try:
        singlePointEnergyDf = get_singlepoint_energies_for_torsion_scan(scanIndex, scanDir, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, moleculePrmtop)
        return singlePointEnergyDf

    except Exception as e:
        raise(e)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_singlepoint_energies_for_torsion_scan(scanIndex, scanDir, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, moleculePrmtop, debug=False):
    fittingRoundDir = p.join(torsionTotalDir, f"fitting_round_{scanIndex+1}")
    os.makedirs(fittingRoundDir, exist_ok=True)


    trajXyzs = sorted([p.join(scanDir, file) 
                        for file in os.listdir(scanDir) 
                        if re.match(r'^orca_scan\.\d\d\d\.xyz$', file)])


    trajPdbs = file_parsers.convert_traj_xyz_to_pdb(trajXyzs, cappedPdb, fittingRoundDir)
    singlePointEnergies = run_mm_singlepoints(trajPdbs, moleculePrmtop)

    singlePointEnergyDf = pd.DataFrame(singlePointEnergies, columns=["TrajIndex", "Energy"])

    scanAngles = Stitching_Assistant.get_scan_angles_from_orca_inp(scanDir)
    if len(scanAngles) != len(singlePointEnergyDf):
        raise(ValueError(f"Number of scan angles ({len(scanAngles)}) does not match number of single point energies ({len(singlePointEnergyDf)})"))
    singlePointEnergyDf["Angle"] = scanAngles

    singlePointEnergyDf["Angle"] = singlePointEnergyDf["Angle"].apply(Stitching_Assistant.rescale_angles_0_360)

    if not debug:
        rmtree(fittingRoundDir)
    return singlePointEnergyDf

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_mm_singlepoints(trajPdbs: list, moleculePrmtop: FilePath) -> float:
    """
    Runs a singlepoint energy calculation at the MM level 
    using OpenMM

    Args:
        prmtop (FilePath): topology file for AMBER
        impcrd (FilePath): coordinate file for AMBER

    Returns:
        singlePointEnergy (float): energy of prmtop // xyz coords
    """
    # Load Amber files and create system
    prmtop: app.Topology = app.AmberPrmtopFile(moleculePrmtop)

    # Create the system.
    system: openmm.System = prmtop.createSystem(nonbondedMethod=app.NoCutoff,
                                                nonbondedCutoff=1 * unit.nanometer,
                                                constraints=None)

    integrator = openmm.LangevinIntegrator(300, 1/unit.picosecond,  0.0005*unit.picoseconds)
    platform = openmm.Platform.getPlatformByName('CPU')

    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    ## set coordinates of simulation 
    singlePointEnergies = []
    for trajPdb in trajPdbs:
        pdbFile = app.PDBFile(trajPdb)
        simulation.context.setPositions(pdbFile.positions)
        state: openmm.State = simulation.context.getState(getPositions=True, getEnergy=True)
        singlePointEnergy = state.getPotentialEnergy() / unit.kilocalories_per_mole

        trajIndex = trajPdb.split(".")[0].split("_")[1]
        singlePointEnergies.append((trajIndex,singlePointEnergy))

    return singlePointEnergies
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
