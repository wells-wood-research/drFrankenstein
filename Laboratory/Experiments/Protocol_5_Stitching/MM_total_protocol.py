## BASIC LIBRARIES ##
import os
from os import path as p
import pandas as pd
import numpy as np
import re
from shutil import rmtree
from scipy.signal import savgol_filter


## MULTIPROCESSION LIBRARIES ##
from concurrent.futures import ProcessPoolExecutor
 ## MULTIPROCESSING AND LOADING BAR LIBRARIES ##
from mpire import WorkerPool
from mpire.utils import make_single_arguments

## drFRANKENSTEIN LIBRARIES ##
from . import Stitching_Assistant

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass


def get_MM_total_energies(config, torsionTag, debug=False):
    mmTotalDir: DirectoryPath = config["pathInfo"]["mmTotalCalculationDir"]
    cappedPdb: FilePath = config["moleculeInfo"]["cappedPdb"]
    completedTorsionScanDirs: list = Stitching_Assistant.get_completed_torsion_scan_dirs(config, torsionTag)
    
    moleculeFrcmod: FilePath = config["pathInfo"]["moleculeFrcmod"]

    ## read calculated partial charges
    chargesCsv: FilePath = config["chargeFittingInfo"]["chargesCsv"]
    chargesDf: pd.DataFrame = pd.read_csv(chargesCsv, index_col="Unnamed: 0")
    ## round charge to 4 decimal places
    chargesDf["Charge"] = chargesDf["Charge"].round(4)

    torsionTotalDir = p.join(mmTotalDir, f"{torsionTag}")
    os.makedirs(torsionTotalDir, exist_ok=True)

    torsionFittingDir = p.join(torsionTotalDir, "fitting_data")
    os.makedirs(torsionFittingDir, exist_ok=True)

    if debug:
        singlePointEnergyDfs = run_serial(completedTorsionScanDirs, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod)
    else:
        singlePointEnergyDfs = run_parallel(completedTorsionScanDirs, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, config, torsionTag)


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
def run_serial(scanDirs, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod):
    argsList = [(scanIndex, scanDir, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod) for  scanIndex, scanDir in enumerate(scanDirs)]
    singlePointEnergyDfs = []
    for args in argsList:
        singlePointEnergyDf = get_singlepoint_energies_for_torsion_scan(*args)
        singlePointEnergyDfs.append(singlePointEnergyDf)
    return singlePointEnergyDfs


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_parallel(scanDirs, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, config, torsionTag):

    tqdmBarOptions = {
        "desc": f"\033[32mFitting {torsionTag}\033[0m",
        "ascii": "-ϟ",  
        "colour": "yellow",
        "unit":  "scan",
        "dynamic_ncols": True
    }
    argsList = [(scanIndex, scanDir, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod) for  scanIndex, scanDir in enumerate(scanDirs)]

    with WorkerPool(n_jobs = config["hardwareInfo"]["nCores"]) as pool:
        singlePointEnergyDfs = pool.map(single_point_worker,
                            make_single_arguments(argsList),
                              progress_bar=True,
                              iterable_len = len(argsList),
                              progress_bar_options=tqdmBarOptions)
    return singlePointEnergyDfs

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def single_point_worker(args):
    scanIndex, scanDir, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod = args
    try:
        singlePointEnergyDf = get_singlepoint_energies_for_torsion_scan(scanIndex, scanDir, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod)
        return singlePointEnergyDf

    except Exception as e:
        raise(e)
        return None
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_singlepoint_energies_for_torsion_scan(scanIndex, scanDir, torsionTotalDir, cappedPdb, chargesDf, moleculeFrcmod, debug=False):
    fittingRoundDir = p.join(torsionTotalDir, f"fitting_round_{scanIndex+1}")
    os.makedirs(fittingRoundDir, exist_ok=True)


    trajXyzs = sorted([p.join(scanDir, file) 
                        for file in os.listdir(scanDir) 
                        if re.match(r'^orca_scan\.\d\d\d\.xyz$', file)])

    prmtop, inpcrd = make_prmtop_inpcrd(trajXyzs[0], fittingRoundDir, cappedPdb, chargesDf, moleculeFrcmod, debug)

    trajPdbs = convert_traj_xyz_to_pdb(trajXyzs, cappedPdb, fittingRoundDir)
    singlePointEnergies = Stitching_Assistant.run_mm_singlepoints(trajPdbs, prmtop, inpcrd)


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
def convert_traj_xyz_to_pdb(trajXyzs, cappedPdb, fittingRoundDir):
    trajPdbs = []
    for trajXyz in trajXyzs:
        trajIndex = trajXyz.split(".")[1]
        trajPdb = p.join(fittingRoundDir, f"orca.{trajIndex}.pdb")
        Stitching_Assistant.update_pdb_coords(cappedPdb, trajXyz, trajPdb)
        trajPdbs.append(trajPdb)

    return trajPdbs
##########################################################
def make_prmtop_inpcrd(trajXyz, fittingRoundDir, cappedPdb, chargesDf, moleculeFrcmod, debug=False):

    trajIndex = trajXyz.split(".")[1]

    os.makedirs(fittingRoundDir, exist_ok=True)
    os.chdir(fittingRoundDir)

    trajPdb = p.join(fittingRoundDir, f"orca_{trajIndex}.pdb")
    Stitching_Assistant.update_pdb_coords(cappedPdb, trajXyz, trajPdb)

    trajMol2 = p.join(fittingRoundDir, f"orca_{trajIndex}.mol2")
    Stitching_Assistant.pdb2mol2(trajPdb, trajMol2, fittingRoundDir)

    chargedMol2 = p.join(fittingRoundDir, f"orca_{trajIndex}_charged.mol2")
    Stitching_Assistant.edit_mol2_partial_charges(trajMol2, chargesDf, chargedMol2)

    renamedMol2 = p.join(fittingRoundDir, f"orca_{trajIndex}_charged_renamed.mol2")
    Stitching_Assistant.edit_mo2_atom_types(chargedMol2, renamedMol2)
    prmtop, inpcrd = Stitching_Assistant.make_prmtop_and_inpcrd(renamedMol2, moleculeFrcmod, fittingRoundDir, trajIndex)

    return prmtop, inpcrd

##########################################################

