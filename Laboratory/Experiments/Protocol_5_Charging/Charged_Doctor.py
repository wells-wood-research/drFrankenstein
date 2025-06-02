## BASIC IMPORTS ##
from os import path as p
import os
import pandas as pd
from shutil import copy


## MULTIPROCESSING AND LOADING BAR LIBRARIES ##c
from mpire import WorkerPool
from mpire.utils import make_single_arguments

## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

## drFRANKENSTEIN LIBRARIES ##
from OperatingTools import Timer, cleaner
from . import Charged_Monster
from . import Charged_Assistant
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
        _, chargesCsv = partial_charge_RESP_protocol(outDir = config["runtimeInfo"]["madeByCharges"]["chargeDir"],
                                useSolvation = True,
                                 config = config,
                                  debug = debug)
        config["runtimeInfo"]["madeByCharges"]["chargesCsv"] = chargesCsv
        
    ## for RESP2 protocol, run charge fitting twice - once with solvation and once without
    elif protocol == "RESP2":
        print("RESP2: with solvation")
        solvatedDf, _ = partial_charge_RESP_protocol(outDir = config["runtimeInfo"]["madeByCharges"]["solvatedDir"],
                                useSolvation = True,
                                 config = config,
                                  debug = debug)
        print("RESP2: gas phase")
        gasPhaseDf, _ = partial_charge_RESP_protocol(outDir = config["runtimeInfo"]["madeByCharges"]["gasPhaseDir"],
                                useSolvation = False,
                                 config = config,
                                  debug = debug)
        resp2Df = Charged_Assistant.apply_resp2_weighted_average(solvatedDf, gasPhaseDf)
        resp2Csv = p.join(config["runtimeInfo"]["madeByCharges"]["chargeDir"], "resp2_charges.csv")
        resp2Df.to_csv(resp2Csv)
        config["runtimeInfo"]["madeByCharges"]["chargesCsv"] = resp2Csv
    ## for SOLVATOR protocol (CHARMM compatible)
    elif protocol == "SOLVATOR":
        chargesCsv, config = partial_charges_SOLVATOR_protocol(outDir = config["runtimeInfo"]["madeByCharges"]["chargeDir"],
                                                          config = config,
                                                          debug = debug)
        config["runtimeInfo"]["madeByCharges"]["chargesCsv"] = chargesCsv



    config = Charged_Assistant.process_charge_csv(config)

    ## clean up
    cleaner.clean_up_charges(config)

    ## update config with checkpoint flag
    config["checkpointInfo"]["chargesComplete"] = True

    print("CHARGES DONE!")
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def partial_charges_SOLVATOR_protocol(outDir: DirectoryPath,
                                      config: dict,
                                      debug: bool) -> tuple[pd.DataFrame, FilePath]:
    """
    Sub-protocol that:
        1. Runs ORCA SOLVATOR calculations to add explicit waters around the molecule
        2. Runs single-point calculations with ORCA to creates wavefunction files 
        3. Runs MultiWFN to assign partial charges to the molecule

    Args:
        outDir (DirectoryPath): Directory for output files
        config (dict): dictionary containing all run information
        debug (bool): Whether to run in debug mode (toggles multiprocessing)

    Returns:
        chargesCsv (FilePath): Path to csv file containing partial charges (used in RESP protocol)
        config (dict): updated config
    """
    ## unpack config
    qmmmSinglepointDir  = config["runtimeInfo"]["madeByCharges"]["qmmmSinglepointDir"]
    fittingDir  = config["runtimeInfo"]["madeByCharges"]["fittingDir"]

    ## run orca single-point calculations on conformers
    config = add_solvation_shell_with_SOLVATOR(config=config,
                                                debug=debug)

    ## solvate with tip3 waters, optimize and run single-point
    config = tip3p_qmmm_protocol(config=config, 
                                debug=debug)
    
    ## create input files for MultiWFN
    conformerListTxt = Charged_Assistant.generate_conformer_list_file(qmmmSinglepointDir, fittingDir)
    chargeConstraintsTxt = Charged_Assistant.generate_charge_constraints_file(config, fittingDir)
    symmetryConstraintsTxt = Charged_Assistant.generate_symmetry_constraints_file(config, fittingDir)

    ## run MultiWFN RESP charge fitting
    rawMultiWfnOutputs = Charged_Monster.run_charge_fitting(config = config,
                                                             conformerListTxt = conformerListTxt,
                                                               chargeConstraintsTxt = chargeConstraintsTxt,
                                                               symmetryConstraintsTxt = symmetryConstraintsTxt,
                                                               fittingDir = fittingDir)
    ## parse MultiWFN output, write to CSV file
    chargesDf = Charged_Monster.parse_multiwfn_output(rawMultiWfnOutputs)
    chargesCsv = p.join(fittingDir, "charges.csv")
    chargesDf.to_csv(chargesCsv)

    return chargesCsv, config



# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def tip3p_qmmm_protocol(config, debug):
    """
    Runs QM/MM calculations for charge calculations using TIP3P water molecules:
    1. Optimisation
    2. Single-Point Energy calculation
    3. converts output to molden files for MultiWFN
    """
    ## unpack config
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    chargeDir = config["runtimeInfo"]["madeByCharges"]["chargeDir"]

    ## make a dir for parameters
    qmmmParameterDir = p.join(chargeDir, "00_qm-mm_parameters")
    os.makedirs(qmmmParameterDir, exist_ok=True)
    config["runtimeInfo"]["madeByCharges"]["qmmmParameterDir"] = qmmmParameterDir
    ## copy over an unsolvated conformer file to make parameters with
    unsolvatedXyz = p.join(qmmmParameterDir, "unsolvated.xyz")
    copy(conformerXyzs[0], unsolvatedXyz)

    ## get qm atoms for QM/MM calculations
    config = Charged_Assistant.get_qm_atoms_for_solvated_system(unsolvatedXyz, config)
    ## create ORCAFF parameters with TIP3P waters
    config = Charged_Monster.create_orca_ff_parameters(unsolvatedXyz, config)
    ## optimise solvated geometries, update config
    qmmmOptXyzs = qmmm_opt_protocol_for_charges(config=config, debug=debug)
    config["runtimeInfo"]["madeByCharges"]["solvatedOptimisedXyzs"] = qmmmOptXyzs
    ## get singlepoint energies
    qmmm_singlepoint_protocol_for_charges(qmmmOptXyzs=qmmmOptXyzs, config=config, debug=debug)

    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("QM/MM Optimisations", "CHARGE_CALCULATIONS")
def qmmm_opt_protocol_for_charges(config, debug):
    ## unpack config
    solvatedXyzs = config["runtimeInfo"]["madeByCharges"]["solvatedXyzs"]
    qmmmOptDir = config["runtimeInfo"]["madeByCharges"]["qmmmOptDir"]

    ### QM-MM Optimisations
    qmmmOptArgsList = [(solvatedXyz, qmmmOptDir, config)
                     for solvatedXyz in solvatedXyzs]
    if debug: ## run serial if debug is on
        qmmmOptXyzs = []
        for qmmmOptArgs in qmmmOptArgsList:
            qmmmOptXyz = Charged_Monster.run_qmmm_opt(qmmmOptArgs)
            qmmmOptXyzs.append(qmmmOptXyz)
    else: ## run in parallel otherwise
        redText = "\033[31m"
        resetTextColor = "\033[0m"

        ## set up progress bar
        tqdmBarOptions = {
            "desc": f"{redText}Running QMMM Optimisations with TIP3P waters{resetTextColor}",
            "ascii": "-🗲→",    
            "colour": "yellow",
            "unit":  "scan",
            "dynamic_ncols": True,
        }
    ## run in parallel
        ## work out how many cpus to use
        nCoresPerCalculation = config["chargeFittingInfo"]["nCoresPerCalculation"]
        availableCpus = config["miscInfo"]["availableCpus"]

        nCpus = min((availableCpus // nCoresPerCalculation), len(qmmmOptArgsList))
        
        with WorkerPool(n_jobs = nCpus) as pool:
            qmmmOptXyzs = pool.map(Charged_Monster.run_qmmm_opt,
                    make_single_arguments(qmmmOptArgsList),
                    progress_bar=True,
                    iterable_len = len(qmmmOptArgsList),
                    progress_bar_options=tqdmBarOptions)
            
    return qmmmOptXyzs
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("QM/MM Single-Points", "CHARGE_CALCULATIONS")
def qmmm_singlepoint_protocol_for_charges(qmmmOptXyzs, config, debug):

    ## unpack config
    qmmmSinglepointDir = config["runtimeInfo"]["madeByCharges"]["qmmmSinglepointDir"]
    fittingDir = config["runtimeInfo"]["madeByCharges"]["fittingDir"]
    ### QM-MM single-points
    qmmmSinglepointArgsList = [(qmmmOptXyz, qmmmSinglepointDir, fittingDir, config)
                     for qmmmOptXyz in qmmmOptXyzs]
    if debug: ## run serial if debug is on
        for qmmmSinglepointArgs in qmmmSinglepointArgsList:
            Charged_Monster.run_qmmm_singlepoint(qmmmSinglepointArgs)
    else: ## run in parallel otherwise
        redText = "\033[31m"
        resetTextColor = "\033[0m"

        ## set up progress bar
        tqdmBarOptions = {
            "desc": f"{redText}Running QMMM Single-Points with TIP3P waters{resetTextColor}",
            "ascii": "-🗲→",    
            "colour": "yellow",
            "unit":  "scan",
            "dynamic_ncols": True,
        }
    ## run in parallel
        ## work out how many cpus to use
        nCoresPerCalculation = config["chargeFittingInfo"]["nCoresPerCalculation"]
        availableCpus = config["miscInfo"]["availableCpus"]
        ## set OMP_NUM_THREADS
        os.environ["OMP_NUM_THREADS"] = str(nCoresPerCalculation)

        nCpus = min((availableCpus // nCoresPerCalculation), len(qmmmSinglepointArgsList))        
        with WorkerPool(n_jobs = nCpus) as pool:        
            optXyzs = pool.map(Charged_Monster.run_qmmm_singlepoint,
                    make_single_arguments(qmmmSinglepointArgsList),
                    progress_bar=True,
                    iterable_len = len(qmmmSinglepointArgsList),
                    progress_bar_options=tqdmBarOptions)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Explicit Solvation", "CHARGE_CALCULATIONS")
def add_solvation_shell_with_SOLVATOR(config: dict,
                                     debug: bool = False) -> dict:
    """
    Runs ORCA single-point calculations on a set of conformers
    Handles multiprocessing and progress bar

    Args:
        orcaDir (DirectoryPath): Directory for output for orca calculations
        fittingDir (DirectoryPath): Directory for output for charge fitting
        config (dict): dictionary containing all run information
    Returns:
        config (dict): updated config
    """
    ## unpack config ##
    solvatorDir = config["runtimeInfo"]["madeByCharges"]["solvatorDir"]

    moleculeInfo = config["moleculeInfo"]
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    nConformers = config["chargeFittingInfo"]["nConformers"]
    chargeFittingInfo = config["chargeFittingInfo"]
    waterDensity = config["chargeFittingInfo"]["waterDensity"]

    nWaters = Charged_Assistant.how_many_waters_for_solvator(conformerXyzs, cappedPdb, solvatorDir, nWatersPerNmSquared=waterDensity)

    ## decide how many conformers to sample
    if nConformers == -1 or nConformers > len(conformerXyzs):
        sampledConformerXyzs = conformerXyzs
        config["chargeFittingInfo"]["nConformers"] = len(conformerXyzs)
    else:
        sampledConformerXyzs = conformerXyzs[:nConformers]

    ## create an argument list for running charge calculations
    solvatorArgList = [(conformerXyz, solvatorDir, chargeFittingInfo,  moleculeInfo, nWaters, config) for 
                conformerXyz in sampledConformerXyzs]
           
    ## run in solvation --> opt in serial
    solvatedXyzs = []
    if debug:
        for arg in solvatorArgList:
            solvatedXyz = Charged_Monster.run_orca_solvator_for_charge_calculations(arg)
            solvatedXyzs.append(solvatedXyz)
    else:
        redText = "\033[31m"
        resetTextColor = "\033[0m"

        ## set up progress bar
        tqdmBarOptions = {
            "desc": f"{redText}Running SOLVATOR for Charge Calculations{resetTextColor}",
            "ascii": "-🗲→",    
            "colour": "yellow",
            "unit":  "scan",
            "dynamic_ncols": True,
        }
    ## run in parallel
        ## work out how many cpus to use
        availableCpus = config["miscInfo"]["availableCpus"]
        nCoresPerCalculation = chargeFittingInfo["nCoresPerCalculation"]
        nCpus = min(availableCpus, len(solvatorArgList))

        ## set OMP_NUM_THREADS
        with WorkerPool(n_jobs = nCpus) as pool:        
            solvatedXyzs = pool.map(Charged_Monster.run_orca_solvator_for_charge_calculations,
                    make_single_arguments(solvatorArgList),
                    progress_bar=True,
                    iterable_len = len(solvatorArgList),
                    progress_bar_options=tqdmBarOptions)
            
    config["runtimeInfo"]["madeByCharges"]["solvatedXyzs"] = solvatedXyzs
    config["runtimeInfo"]["madeByCharges"]["nWaters"] = nWaters
    return config


    
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def  partial_charge_RESP_protocol(outDir: DirectoryPath,
                                useSolvation: bool,
                                  config: dict,
                                    debug: bool = False) -> tuple[pd.DataFrame, FilePath]:
    """
    Sub-protocol that first runs ORCA calculations and then uses MultiWFN to calculate partial charges
    Depending on the useSolvation arg, runs with or without implicit solvent
    Compatible with AMBER-style parameterisation methods

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
    orcaDir = p.join(outDir, "01_orca_calculations")
    os.makedirs(orcaDir, exist_ok=True)

    fittingDir = p.join(outDir, "02_charge_fitting")
    os.makedirs(fittingDir, exist_ok=True)

    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]

    ## run orca single-point calculations on conformers
    run_qm_calculations_for_RESP(conformerXyzs = conformerXyzs,
                                 orcaDir=orcaDir,
                                            fittingDir=fittingDir,
                                               config=config,
                                                useSolvation=useSolvation,
                                                debug=debug)
    ## create input files for MultiWFN
    conformerListTxt = Charged_Assistant.generate_conformer_list_file(orcaDir, fittingDir)
    chargeConstraintsTxt = Charged_Assistant.generate_charge_constraints_file(config, fittingDir)
    symmetryConstraintsTxt = Charged_Assistant.generate_symmetry_constraints_file(config, fittingDir)
    ## run MultiWFN RESP charge fitting
    rawMultiWfnOutputs = Charged_Monster.run_charge_fitting(config=config,
                                                            conformerListTxt= conformerListTxt,
                                                             chargeConstraintsTxt = chargeConstraintsTxt,
                                                             symmetryConstraintsTxt=symmetryConstraintsTxt,
                                                              fittingDir = fittingDir)
    ## parse MultiWFN output, write to CSV file
    chargesDf = Charged_Monster.parse_multiwfn_output(rawMultiWfnOutputs)
    chargesCsv = p.join(fittingDir, "charges.csv")
    chargesDf.to_csv(chargesCsv)

    return chargesDf, chargesCsv
 
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("QM Calculations", "CHARGE_CALCULATIONS")
def run_qm_calculations_for_RESP(conformerXyzs: list[FilePath],
                                    orcaDir: DirectoryPath,
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
    # conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
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
    

    redText = "\033[31m"
    resetTextColor = "\033[0m"

    ## set up progress bar
    tqdmBarOptions = {
        "desc": f"{redText}Running QM for Charge Calculations{resetTextColor}",
        "ascii": "-🗲→",    
        "colour": "yellow",
        "unit":  "scan",
        "dynamic_ncols": True,
    }

    ## run in serial
    if debug:
        for arg in argsList:
            Charged_Monster.run_orca_singlepoint_for_charge_calculations(arg)
    ## run in parallel
    else:
        nCpus = min(len(argsList), config["miscInfo"]["availableCpus"])
        with WorkerPool(n_jobs = nCpus) as pool:        

            pool.map(Charged_Monster.run_orca_singlepoint_for_charge_calculations,
                    make_single_arguments(argsList),
                    progress_bar=True,
                    iterable_len = len(argsList),
                    progress_bar_options=tqdmBarOptions)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    raise NotImplementedError




