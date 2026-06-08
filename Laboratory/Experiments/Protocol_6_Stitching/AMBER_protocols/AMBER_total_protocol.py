## BASIC LIBRARIES ##
import os
from os import path as p
import pandas as pd
import numpy as np
import re
from shutil import rmtree
from scipy.signal import savgol_filter
from tqdm import tqdm
from typing import List, Tuple

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

GPU_DISABLE_ENV_VARS = (
    "CUDA_VISIBLE_DEVICES",
    "HIP_VISIBLE_DEVICES",
    "ROCR_VISIBLE_DEVICES",
    "GPU_DEVICE_ORDINAL",
)


def configure_openmm_gpu_access(gpuPlatform):
    if isinstance(gpuPlatform, str) and gpuPlatform.upper() == "CPU":
        for envVar in GPU_DISABLE_ENV_VARS:
            os.environ[envVar] = "-1"
        os.environ["OPENMM_DEFAULT_PLATFORM"] = "CPU"
        return "CPU"
    if isinstance(gpuPlatform, str):
        platformNames = {"CUDA": "CUDA", "HIP": "HIP", "OPENCL": "OpenCL"}
        return platformNames.get(gpuPlatform.upper(), gpuPlatform)
    return gpuPlatform


def get_MM_total_energies(config:dict, torsionTag, moleculeFrcmod, debug=False):
    drSplash.show_getting_mm_total(torsionTag)

    mmTotalDir: DirectoryPath = config["runtimeInfo"]["madeByStitching"]["mmTotalCalculationDir"]
    cappedPdb: FilePath = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    completedTorsionScanDirs: list = Stitching_Assistant.get_completed_torsion_scan_dirs(config, torsionTag)
    
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
    torsionIndexes = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"][torsionTag]["ATOM_INDEXES"]

    # Thread gpuPlatform from config (set earlier in set_config_defaults)
    gpuPlatform = None
    if isinstance(config, dict):
        gpuPlatform = config.get("miscInfo", {}).get("gpuPlatform", None)
    gpuPlatform = configure_openmm_gpu_access(gpuPlatform)


    if debug:
        singlePointEnergyDfs = run_serial(completedTorsionScanDirs, torsionTotalDir, cappedPdb, moleculePrmtop, torsionIndexes, gpuPlatform)
    else:
        singlePointEnergyDfs = run_parallel(completedTorsionScanDirs, torsionTotalDir, cappedPdb, moleculePrmtop, torsionIndexes, gpuPlatform)

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
def run_serial(scanDirs, torsionTotalDir, cappedPdb, moleculePrmtop, torsionIndexes, gpuPlatform=None):
    argsList = [(scanIndex, scanDir, torsionTotalDir, cappedPdb, moleculePrmtop, gpuPlatform) for  scanIndex, scanDir in enumerate(scanDirs)]
    singlePointEnergyDfs = []
    for args in argsList:
        singlePointEnergyDf = get_singlepoint_energies_for_torsion_scan(*args, torsionIndexes)
        singlePointEnergyDfs.append(singlePointEnergyDf)
    return singlePointEnergyDfs


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_parallel(scanDirs, torsionTotalDir, cappedPdb, moleculePrmtop, torsionIndexes, gpuPlatform=None):
    argsList = [(scanIndex, scanDir, torsionTotalDir, cappedPdb,  moleculePrmtop, gpuPlatform) for scanIndex, scanDir in enumerate(scanDirs)]

    # Use multiprocessing Pool without tqdm progress bar
    with multiprocessing.Pool() as pool:
        singlePointEnergyDfs = list(pool.imap(single_point_worker, [(args, torsionIndexes) for args in argsList]))
    
    return singlePointEnergyDfs

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def single_point_worker(args):
    (scanIndex, scanDir, torsionTotalDir, cappedPdb, moleculePrmtop, gpuPlatform), torsionIndexes = args
    try:
        singlePointEnergyDf = get_singlepoint_energies_for_torsion_scan(scanIndex, scanDir, torsionTotalDir, cappedPdb, moleculePrmtop, torsionIndexes, gpuPlatform)
        return singlePointEnergyDf

    except Exception as e:
        raise(e)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_singlepoint_energies_for_torsion_scan(scanIndex, scanDir, torsionTotalDir, cappedPdb, moleculePrmtop, torsionIndexes, gpuPlatform=None, debug=False):
    fittingRoundDir = p.join(torsionTotalDir, f"fitting_round_{scanIndex+1}")
    os.makedirs(fittingRoundDir, exist_ok=True)


    trajXyzs = sorted([p.join(scanDir, file) 
                        for file in os.listdir(scanDir) 
                        if re.match(r'^orca_scan\.\d\d\d\.xyz$', file)])


    trajPdbs = file_parsers.convert_traj_xyz_to_pdb(trajXyzs, cappedPdb, fittingRoundDir)
    minimisedEnergies = run_mm_constrained_em(trajPdbs, moleculePrmtop, torsionIndexes=torsionIndexes, gpuPlatform=gpuPlatform)

    singlePointEnergyDf = pd.DataFrame(minimisedEnergies, columns=["TrajIndex", "Energy"])

    scanAngles = np.array(Stitching_Assistant.get_scan_angles_from_orca_inp(scanDir))
    if len(scanAngles) != len(singlePointEnergyDf):
        raise(ValueError(f"Number of scan angles ({len(scanAngles)}) does not match number of single point energies ({len(singlePointEnergyDf)})"))

    # Ensure angles are numeric and in 0-360 range without per-row apply (avoids pandas concat issues)
    scanAngles = np.mod(scanAngles, 360)
    singlePointEnergyDf["Angle"] = scanAngles

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
    nParticles = simulation.context.getSystem().getNumParticles()
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
def run_mm_constrained_em(trajPdbs: List[FilePath], moleculePrmtop: FilePath, torsionIndexes: List[int], gpuPlatform: str | None = None) -> List[Tuple[str, float]]:
    """
    Runs a constrained energy minimization for each provided structure
    and returns the final potential energy.

    The minimization is constrained by freezing the positions of atoms
    whose indices are provided in torsionIndexes.

    Args:
        trajPdbs (List[FilePath]): A list of PDB file paths for the coordinate sets.
        moleculePrmtop (FilePath): The topology file path (AMBER prmtop).
        torsionIndexes (List[int]): A list of 0-based atom indices to constrain
                                     during the minimization.
        gpuPlatform (str|None): Preferred platform name ('CUDA','HIP','OpenCL','CPU') or None to use defaults.

    Returns:
        List[Tuple[str, float]]: A list of tuples, where each tuple contains the
                                 trajectory index (as a string) and its corresponding
                                 minimized potential energy in kcal/mol.
    """
    # 1. Load Amber files to create the base system
    prmtop = app.AmberPrmtopFile(str(moleculePrmtop))
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff,
                                 constraints=None) # Use None for constraints here to not conflict with our custom ones

    # 2. Apply constraints: Set the mass of specified atoms to zero.
    # OpenMM's minimizer and integrators will not move particles with zero mass.
    for atom_index in torsionIndexes:
        system.setParticleMass(atom_index, 0.0 * unit.dalton)

    # 3. Set up the simulation object
    # An integrator is required, but it's not used by the local energy minimizer.
    integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)

    gpuPlatform = configure_openmm_gpu_access(gpuPlatform)
    platform = openmm.Platform.getPlatformByName(gpuPlatform)


    simulation = app.Simulation(prmtop.topology, system, integrator, platform)

    minimizedEnergies = []
    # 4. Loop through each trajectory PDB file
    for trajPdb in trajPdbs:
        pdbFile = app.PDBFile(str(trajPdb))

        # Set the starting positions for this frame
        simulation.context.setPositions(pdbFile.positions)
        
        # 5. Run the energy minimization
        simulation.minimizeEnergy()

        # 6. Get the potential energy after minimization
        final_state = simulation.context.getState(getEnergy=True)
        minimizedEnergy = final_state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        
        # Extract a unique index from the filename (e.g., from "snap_123.pdb")
        try:
            trajIndex = str(trajPdb).split(".")[0].split("_")[-1]
        except IndexError:
            trajIndex = "unknown" # Fallback if filename format is different
            
        minimizedEnergies.append((trajIndex, minimizedEnergy))

    return minimizedEnergies
