## BASIC IMPORTS ##
from subprocess import call, DEVNULL 
from os import path as p
import os
from shutil import  copy
import pandas as pd
import numpy as np
import mdtraj as md
from pdbUtils import pdbUtils
## drFrankenstein LIBRARIES ##
from OperatingTools import file_parsers

## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass


def process_charge_csv(config: dict) -> dict:
    """
    Converts auto-generated charge csv into a format that can be used in the report protocol

    Args:
        config (dict): config containing all run information

    Returns:
        config (dict): updated config dict
    """
    ## unpack config ##
    chargeCsv = config["runtimeInfo"]["madeByCharges"]["chargesCsv"]
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]

    cappedDf = pdbUtils.pdb2df(cappedPdb)
    chargeDf = pd.read_csv(chargeCsv, index_col="Unnamed: 0")
    chargeDf["ATOM_NAME"] = cappedDf["ATOM_NAME"]

    chargeDf.to_csv(chargeCsv)

    return config




def find_tip3p_water_params() -> FilePath:

    """
    Looks through drFrankenstein's src directory to find TIP3P water params in orcaFF format

    Args:
        None
    Returns:
        nmePdb (FilePath): path to NME.pdb
        acePdb (FilePath): path to ACE.pdb
    """
    ## get location of this file
    chargeSrcDir = p.dirname(p.abspath(__file__))
    tip3pParams = p.join(chargeSrcDir, "TIP3P.ORCAFF.prms")
    return tip3pParams

def get_qm_atoms_for_solvated_system(xyzFile: FilePath, config:dict) -> dict:
    """
    Simple function that gets a contiguous list of atom indexes (starting with 0)
    for orca to use as a qm region. 
    """
    xyzDf = file_parsers.xyz2df(xyzFile)
    nAtoms = len(xyzDf.index)
    qmAtoms = "{" + "0:" + str(nAtoms-1) + "}"

    config["runtimeInfo"]["madeByCharges"]["qmAtoms"] = qmAtoms
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def how_many_waters_for_solvator(conformerXyzs: list[FilePath],
                                  conformerPdb: FilePath,
                                    outDir: DirectoryPath,
                                    nWatersPerNmSquared: int = 10) -> int:
    """
    Simple function that decides how many water molecules to add in a SOLVATOR calculation
    Currently we just do 2 x nHeavyAtoms
    TODO: This could be more detailed?
    
    Args: 
        xyzFile (FilePath): an XYZ file of our molecule
        multiplier (int): nWaters = nHeavyAtoms * multiplier

    Returns:
        nWaters (int):  how many water molecules to add in a SOLVATOR calculation
    """
    ## convert XYZs to PDBs so we can load them into mdtraj
    conformerPdbs = file_parsers.convert_traj_xyz_to_pdb(conformerXyzs, conformerPdb, outDir)

    totalSasas = []
    for pdbFile in conformerPdbs:
        traj = md.load(pdbFile)
        sasaPerAtom = md.shrake_rupley(traj)
        totalSasa = sasaPerAtom.sum()
        totalSasas.append(totalSasa)
    meanTotalSasa = sum(totalSasas)/len(totalSasas)

    nWaters = round(meanTotalSasa * nWatersPerNmSquared)

    ## clean up
    for pdbFile in conformerPdbs:
        os.remove(pdbFile)
    
    return nWaters
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def generate_charge_constraints_file(config: dict, outDir: DirectoryPath) -> FilePath:
    """
    Uses charge groups from config to generate charge constraints file to be used by MultiWFN

    Args:
        config (dict): config containing all run information
        outDir (DirectoryPath): path to output directory

    Returns:
        chargeConstraintsTxt (FilePath): path to charge constraints file
    
    """
    ## unpack config
    chargeGroups: dict = config["runtimeInfo"]["madeByCharges"]["chargeGroups"]


    chargeConstraintsTxt = p.join(outDir, "charge_group_constraints.txt")
    with open(chargeConstraintsTxt, "w") as f:
        for _, chargeGroupData in chargeGroups.items():
            indexes = chargeGroupData["indexes"]
            if len(indexes) == 0:
                continue
            charge = chargeGroupData["charge"]
            formattedIndexes = ",".join(map(str, indexes))
            f.write(f"{formattedIndexes} {str(charge)}\n")
          

    return chargeConstraintsTxt
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def generate_conformer_list_file(orcaCalculationsDir: DirectoryPath,
                                  chargeFittingDir: DirectoryPath) -> FilePath:
    """
    Reads through single point calculations for each conformer 
    and gets the single-point energies
    converts these into Boltzmann probabilities 
    and creates a conformer_list.txt file for MultiWFN

    Args:
        orcaCalculationsDir (DirectoryPath): path to dir containing single point calculations
        chargeFittingDir (DirectoryPath): path to dir for charge fitting calculations

    Returns:
        conformerListTxt (FilePath): path to conformer list file
    """

    singlePointData = {}
    for conformerName in os.listdir(orcaCalculationsDir):
        singlePointData[conformerName] = {}
        calculationDir = p.join(orcaCalculationsDir, conformerName)
        if not p.isdir:
            continue
        ## find ORCA.out file
        singlePointOutFile = [p.join(calculationDir, file) for file
                               in os.listdir(calculationDir)
                               if p.splitext(file)[1] == ".out"][0]
        singlePointEnergy = find_final_single_point_energy(singlePointOutFile)
       
        singlePointData[conformerName]["Energy"] = singlePointEnergy
        singlePointData[conformerName]["Path"] = p.join(chargeFittingDir, f"{conformerName}.molden.input")

    singlePointDf =pd.DataFrame.from_dict(singlePointData, orient='index')

    singlePointDf = calculate_boltzmann_probabilities(singlePointDf)

    conformerListTxt = p.join(chargeFittingDir, "conformer_list.txt")
    with open(conformerListTxt, "w") as f:
        for _, row in singlePointDf.iterrows():
            f.write(f"{row['Path']} {str(row['Probability'])}\n")

    return conformerListTxt

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_boltzmann_probabilities(df, temperature=298.15):
    boltzmannConstant = 0.0083144621  # kJ/(mol*K)
    energies = df['Energy'].values
    minEnergy = np.min(energies)
    deltaEnergies = energies - minEnergy
    boltzmannFactors = np.exp(-deltaEnergies / 
                              (boltzmannConstant * temperature))
    probabilities = boltzmannFactors / np.sum(boltzmannFactors)
    df['Probability'] = probabilities
    return df    

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def find_final_single_point_energy(outFilePath):
    with open(outFilePath, 'r') as file:
        for line in file:
            if "FINAL SINGLE POINT ENERGY" in line:
                return float(line.split()[-1])


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def apply_resp2_weighted_average(solvatedDf: pd.DataFrame,
                                  gasPhaseDf: pd.DataFrame,
                                    proportions: list[float] = [0.6, 0.4]) -> pd.DataFrame:
    """
    Takes charges that have been calculated with and without solvation
    and applies a weighted average to create RESP2 charges

    Args:
        solvatedDf (pd.DataFrame): DataFrame of charges calculated with solvation
        gasPhaseDf (pd.DataFrame): DataFrame of charges calculated in the gas phase
        proportions (list, optional): List of proportions to apply to solvatedDf and gasPhaseDf. Defaults to [0.6, 0.4]

    Returns:
       resp2Df (pd.DataFrame): DataFrame of RESP2 charges
    """

    resp2Df = pd.DataFrame()
    resp2Df[["atomIndex","atomElement"]] = solvatedDf[["atomIndex","atomElement"]]
    resp2Df["Charge"] = proportions[0] * solvatedDf["Charge"] + proportions[1] * gasPhaseDf["Charge"]
    return resp2Df
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def set_up_directories(config: dict, protocol: str) -> dict:
    """
    Creates a directory structure for charge calculations
    The structure depends on whether RESP or RESP2 is being run

    Args:
        config (dict): Dictionary containing all information needed
        protocol (str): either RESP or RESP2

    Returns:
        config (dict): Updated config dictionary
    """
    ## get pathInfo dict
    pathInfo: dict = config["pathInfo"]

    outputDir: DirectoryPath = pathInfo["outputDir"]
    ## charge dir - acts as a topDir for all charge-related processes
    chargeDir: DirectoryPath = p.join(outputDir, "05_charge_calculations")
    os.makedirs(chargeDir, exist_ok=True)

    #### FOR RESP FITTING, just run once ##
    if protocol == "RESP":
        orcaCalculationsDir: DirectoryPath = p.join(chargeDir, "01_orca_calculations")
        os.makedirs(orcaCalculationsDir, exist_ok=True)
        chargeFittingDir: DirectoryPath = p.join(chargeDir, "02_charge_fitting")
        os.makedirs(chargeFittingDir, exist_ok=True)
        config["runtimeInfo"]["madeByCharges"]["chargeDir"] = chargeFittingDir
        return config
    
    elif protocol == "RESP2":

        ## directory to run ORCA QM calculations
        solvatedDir: DirectoryPath = p.join(chargeDir, "RESP2_solvated")
        gasPhaseDir: DirectoryPath = p.join(chargeDir, "RESP2_gas_phase")

        ## update config
        config["runtimeInfo"]["madeByCharges"].update({
            "chargeDir": chargeDir,
            "solvatedDir": solvatedDir,
            "gasPhaseDir": gasPhaseDir,
        })
        return config
    elif protocol == "SOLVATOR":
        ## create output directories for orca single-point calculations and MultiWFN charge fitting
        solvatorDir = p.join(chargeDir, "01_solvator_calculations")
        os.makedirs(solvatorDir, exist_ok=True)
        qmmmOptDir = p.join(chargeDir, "02_QMMM_optimisations")
        os.makedirs(qmmmOptDir, exist_ok=True)
        qmmmSinglepointDir = p.join(chargeDir, "03_QMMM_singlepoints")
        os.makedirs(qmmmSinglepointDir, exist_ok=True)
        fittingDir = p.join(chargeDir, "04_charge_fitting")
        os.makedirs(fittingDir, exist_ok=True)
        ## update config
        config["runtimeInfo"]["madeByCharges"].update({
            "chargeDir": chargeDir,
            "solvatorDir": solvatorDir,
            "qmmmOptDir": qmmmOptDir,
            "qmmmSinglepointDir": qmmmSinglepointDir,
            "fittingDir": fittingDir
        })
        return config