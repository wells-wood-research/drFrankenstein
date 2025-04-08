## BASIC IMPORTS ##
from os import path as p
import os
import pandas as pd

## drFrankenstein LIBRARIES ##
from . import Charged_Monster
from . import Charged_Assistant

## MULTIPROCESSING AND LOADING BAR LIBRARIES ##c
from mpire import WorkerPool
from mpire.utils import make_single_arguments

## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

## drFRANKENSTEIN LIBRARIES ##
from OperatingTools import Timer
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function()
def charge_protocol(config: dict, debug: bool = False) -> dict:
    """
    Main protocol for charge fitting
    Runs either RESP or RESP2 style charge fitting methods
    Uses a combination of ORCA single-points and MultiWFN to calculate partial charges

    Args:
        config (dict): Config file for the charge fitting protocol
        debug (bool): Whether to run in debug mode (toggles multiprocessing)
    
    Returns:
        config (dict): updated config
    """
    ## create an entry in runtimeInfo for conformers
    config["runtimeInfo"]["madeByCharges"] = {}


    ## unpack config ##
    protocol = config["chargeFittingInfo"]["chargeFittingProtocol"]

    ## get indexes of charge groups (doesn't change so out of if block)
    config = Charged_Monster.get_charge_group_indexes(config["runtimeInfo"]["madeByCapping"]["cappedPdb"], config)

    ## set up directories for either RESP or RESP2 calculations
    config = Charged_Assistant.set_up_directories(config, protocol)

    ## For RESP protocol, just run charge fitting once
    if protocol == "RESP":
        _, chargesCsv = calculate_partial_charges(outDir = config["pathInfo"]["chargeDir"],
                                useSolvation = True,
                                 config = config,
                                  debug = debug)
        config["runtimeInfo"]["madeByCharges"]["chargesCsv"] = chargesCsv
        
    ## for RESP2 protocol, run charge fitting twice - once with solvation and once without
    elif protocol == "RESP2":
        print("RESP2: with solvation")
        solvatedDf, _ = calculate_partial_charges(outDir = config["runtimeInfo"]["madeByCharges"]["solvatedDir"],
                                useSolvation = True,
                                 config = config,
                                  debug = debug)
        print("RESP2: gas phase")
        gasPhaseDf, _ = calculate_partial_charges(outDir = config["runtimeInfo"]["madeByCharges"]["gasPhaseDir"],
                                useSolvation = False,
                                 config = config,
                                  debug = debug)
        

        resp2Df = Charged_Assistant.apply_resp2_weighted_average(solvatedDf, gasPhaseDf)
        resp2Csv = p.join(config["runtimeInfo"]["madeByCharges"]["chargeDir"], "resp2_charges.csv")
        resp2Df.to_csv(resp2Csv)
        config["runtimeInfo"]["madeByCharges"]["chargesCsv"] = resp2Csv

    ## update config with checkpoint flag
    config["checkpointInfo"]["chargesComplete"] = True
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def  calculate_partial_charges(outDir: DirectoryPath,
                                useSolvation: bool,
                                  config: dict,
                                    debug: bool = False) -> tuple[pd.DataFrame, FilePath]:
    """
    Sub-protocol that first runs ORCA calculations and then uses MultiWFN to calculate partial charges

    Args:
        outDir (DirectoryPath): Directory for output files
        useSolvation (bool): Whether to use solvation
        config (dict): dictionary containing all run information
        debug (bool): Whether to run in debug mode (toggles multiprocessing)

    Returns:
        chargesDf (pd.DataFrame): DataFrame containing partial charges (used in RESP2 protocol)
        chargesCsv (FilePath): Path to csv file containing partial charges (used in RESP protocol)

    """
    ## create output directories for orca single-point calculations and MultiWFN charge fitting
    orcaDir = p.join(outDir, "orca_calculations")
    os.makedirs(orcaDir, exist_ok=True)

    fittingDir = p.join(outDir, "charge_fitting")
    os.makedirs(fittingDir, exist_ok=True)
    ## run orca single-point calculations on conformers
    run_qm_calculations_for_charge_fitting(orcaDir=orcaDir,
                                            fittingDir=fittingDir,
                                               config=config,
                                                useSolvation=useSolvation,
                                                debug=debug)
    ## create input files for MultiWFN
    conformerListTxt = Charged_Assistant.generate_conformer_list_file(orcaDir, fittingDir)
    chargeConstraintsTxt = Charged_Assistant.generate_charge_constraints_file(config, fittingDir)
    ## run MultiWFN RESP charge fitting
    rawMultiWfnOutputs = Charged_Monster.run_charge_fitting(config, conformerListTxt, chargeConstraintsTxt, fittingDir)
    ## parse MultiWFN output, write to CSV file
    chargesDf = Charged_Monster.parse_multiwfn_output(rawMultiWfnOutputs)
    chargesCsv = p.join(fittingDir, "charges.csv")
    chargesDf.to_csv(chargesCsv)

    return chargesDf, chargesCsv
 
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_qm_calculations_for_charge_fitting(orcaDir: DirectoryPath,
                                           fittingDir: DirectoryPath,
                                                 config: dict,
                                                   useSolvation: bool = True,
                                                     debug: bool = False) -> None:
    """
    Runs ORCA single-point calculations on a set of conformers
    Handles multiprocessing and progress bar

    Args:
        orcaDir (DirectoryPath): Directory for output for orca calculations
        fittingDir (DirectoryPath): Directory for output for charge fitting
        config (dict): dictionary containing all run information
    Returns:
        None [we will find output files later]
    
    """
    
    ## unpack config ##
    moleculeInfo = config["moleculeInfo"]
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    nConformers = config["chargeFittingInfo"]["nConformers"]
    chargeFittingInfo = config["chargeFittingInfo"]

    ## decide how many conformers to sample
    if nConformers == -1 or nConformers > len(conformerXyzs):
        sampledConformerXyzs = conformerXyzs
        config["chargeFittingInfo"]["nConformers"] = len(conformerXyzs)
    else:
        sampledConformerXyzs = conformerXyzs[:nConformers]

    ## create an argument list for running charge calculations
    argsList = [(conformerXyz, orcaDir, fittingDir, chargeFittingInfo,  moleculeInfo, useSolvation, config) for conformerXyz in sampledConformerXyzs]
    
    ## set up progress bar
    tqdmBarOptions = {
        "desc": f"Running Charge Calculations",
        "ascii": "-ϟ",  
        "colour": "yellow",
        "unit":  "scan",
        "dynamic_ncols": True,
        "leave": True
    }
    ## run in serial
    if debug:
        for arg in argsList:
            Charged_Monster.run_orca_singlepoint_for_charge_calculations(arg)
    ## run in parallel
    else:
        nCpus = min(len(argsList), config["hardwareInfo"]["nCores"])
        with WorkerPool(n_jobs = nCpus) as pool:
            pool.map(Charged_Monster.run_orca_singlepoint_for_charge_calculations,
                    make_single_arguments(argsList),
                    progress_bar=True,
                    iterable_len = len(argsList),
                    progress_bar_options=tqdmBarOptions)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    raise NotImplementedError




