## BASIC LIBRARIES ##
import os
from os import path as p
import pandas as pd
import numpy as np
import re
from shutil import rmtree
from scipy.signal import savgol_filter

 ## MULTIPROCESSING AND LOADING BAR LIBRARIES ##
import multiprocessing 
from mpire.utils import make_single_arguments

## OPENMM LIBRARIES
import openmm.app as app
import openmm as openmm
import  openmm.unit  as unit

## drFRANKENSTEIN LIBRARIES ##
from .. import Stitching_Assistant
from OperatingTools import drSplash
from OperatingTools import file_parsers

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_MM_total_energies(config:dict,
                           torsionTag: str, 
                           moleculePrm: FilePath,
                             debug = True) -> np.ndarray:
    """
    Main protocol for calculating CHARMM energies for a torsion scan

    Args:
        config (dict): contains all info for run
        torsionTag (str): identifier for torsion of interest
        debug (bool): controls multiprocessing 
    Returns:
        smoothedEnergies (mp.ndarray): CHARMM energy
    
    """
    drSplash.show_getting_mm_total(torsionTag)

    mmTotalDir: DirectoryPath = config["runtimeInfo"]["madeByStitching"]["mmTotalCalculationDir"]
    completedTorsionScanDirs: list = Stitching_Assistant.get_completed_torsion_scan_dirs(config, torsionTag)

    torsionTotalDir = p.join(mmTotalDir, f"{torsionTag}")
    os.makedirs(torsionTotalDir, exist_ok=True)

    torsionFittingDir = p.join(torsionTotalDir, "fitting_data")
    os.makedirs(torsionFittingDir, exist_ok=True)

    if debug:
        singlePointEnergyDfs = run_serial(completedTorsionScanDirs, torsionTotalDir, moleculePrm, config)
    else:
        singlePointEnergyDfs = run_parallel(completedTorsionScanDirs, torsionTotalDir, moleculePrm, config)
    
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

    smoothedEnergies = finalScanEnergiesDf["smoothedEnergy"].to_numpy()

    return  smoothedEnergies

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_parallel(scanDirs, torsionTotalDir, moleculePrm, config):
    argsList = [(scanIndex, scanDir, torsionTotalDir, moleculePrm, config) for scanIndex, scanDir in enumerate(scanDirs)]

    with multiprocessing.Pool(processes=config["miscInfo"]["availableCpus"]) as pool:
        singlePointEnergyDfs = pool.starmap(single_point_worker, make_single_arguments(argsList))
    
    return singlePointEnergyDfs


def single_point_worker(args):
    scanIndex, scanDir, torsionTotalDir, moleculePrm, config = args
    try:
        singlePointEnergyDf = get_singlepoint_energies_for_torsion_scan(scanIndex, scanDir, torsionTotalDir, moleculePrm, config)
        return singlePointEnergyDf

    except Exception as e:
        raise(e)    

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_serial(scanDirs, torsionTotalDir, moleculePrm, config):
    argsList = [(scanIndex, scanDir, torsionTotalDir, moleculePrm, config) for  scanIndex, scanDir in enumerate(scanDirs)]
    singlePointEnergyDfs = []
    for args in argsList:
        singlePointEnergyDf = get_singlepoint_energies_for_torsion_scan(*args)
        singlePointEnergyDfs.append(singlePointEnergyDf)
    return singlePointEnergyDfs

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_singlepoint_energies_for_torsion_scan(scanIndex, scanDir, torsionTotalDir, moleculePrm, config, debug = False):
    ## unpack config
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]


    fittingRoundDir = p.join(torsionTotalDir, f"fitting_round_{scanIndex+1}")
    os.makedirs(fittingRoundDir, exist_ok=True)


    trajXyzs = sorted([p.join(scanDir, file) 
                        for file in os.listdir(scanDir) 
                        if re.match(r'^orca_scan\.\d\d\d\.xyz$', file)])


    trajPdbs = file_parsers.convert_traj_xyz_to_pdb(trajXyzs, cappedPdb, fittingRoundDir)
    singlePointEnergies = run_mm_singlepoints(trajPdbs, moleculePrm, config)


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

def run_mm_singlepoints(trajPdbs: list, moleculePrm: FilePath, config: dict) -> list[float]:
    """
    Runs a singlepoint energy calculation at the MM level 
    using OpenMM

    Args:
        trajPdbs (list[FilePath]): PDB files that lie along a torsion scan trajectory
        config (dicr)
    Returns:
        singlePointEnergies (list[float]): MM energy of torsion scan
    """

    ## unpack config
    moleculePsf = config["runtimeInfo"]["madeByStitching"]["moleculePsf"]
    moleculeRtf = config["runtimeInfo"]["madeByStitching"]["moleculeRtf"]

    # cgenffRtf = config["runtimeInfo"]["madeByStitching"]["cgenffRtf"]
    # cgenffPrm = config["runtimeInfo"]["madeByStitching"]["cgenffPrm"]


    ## load CHARMM parameters
    psf: app.Topology = app.CharmmPsfFile(moleculePsf)
    charmmParams = app.CharmmParameterSet(moleculeRtf, moleculePrm)

    # Create the system.
    system: openmm.System = psf.createSystem(charmmParams,
                                                   nonbondedMethod=app.NoCutoff,
                                                nonbondedCutoff=1 * unit.nanometer,
                                                constraints=None)
    ## set up intergrator
    integrator = openmm.LangevinIntegrator(300, 1/unit.picosecond,  0.0005*unit.picoseconds)
    ## use CPU platform - no point in GPU for this tiny system!
    platform = openmm.Platform.getPlatformByName('CPU')
    ## make simulation object
    simulation = app.Simulation(psf.topology, system, integrator, platform)
    ## set coordinates of simulation along torsion scan
    singlePointEnergies = []
    for trajPdb in trajPdbs:
        pdbFile = app.PDBFile(trajPdb)
        simulation.context.setPositions(pdbFile.positions)
        state: openmm.State = simulation.context.getState(getPositions=True, getEnergy=True)
        singlePointEnergy = state.getPotentialEnergy() / unit.kilocalories_per_mole

        trajIndex = trajPdb.split(".")[0].split("_")[1]
        singlePointEnergies.append((trajIndex,singlePointEnergy))

    return singlePointEnergies
