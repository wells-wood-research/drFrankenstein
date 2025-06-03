from os import path as p
import os
import pandas as pd
from shutil import copy
from mpire import WorkerPool
from mpire.utils import make_single_arguments

class FilePath:
    pass

class DirectoryPath:
    pass
from OperatingTools import Timer, cleaner
from . import Charged_Monster
from . import Charged_Assistant

def charge_protocol(config: dict, debug: bool=False) -> dict:
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
    config['runtimeInfo']['madeByCharges'] = {}
    protocol = config['chargeFittingInfo']['chargeFittingProtocol']
    config = Charged_Monster.get_charge_group_indexes(config['runtimeInfo']['madeByCapping']['cappedPdb'], config)
    config = Charged_Assistant.set_up_directories(config, protocol)
    if protocol == 'RESP':
        (_, chargesCsv) = partial_charge_RESP_protocol(outDir=config['runtimeInfo']['madeByCharges']['chargeDir'], useSolvation=True, config=config, debug=debug)
        config['runtimeInfo']['madeByCharges']['chargesCsv'] = chargesCsv
    elif protocol == 'RESP2':
        print('RESP2: with solvation')
        (solvatedDf, _) = partial_charge_RESP_protocol(outDir=config['runtimeInfo']['madeByCharges']['solvatedDir'], useSolvation=True, config=config, debug=debug)
        print('RESP2: gas phase')
        (gasPhaseDf, _) = partial_charge_RESP_protocol(outDir=config['runtimeInfo']['madeByCharges']['gasPhaseDir'], useSolvation=False, config=config, debug=debug)
        resp2Df = Charged_Assistant.apply_resp2_weighted_average(solvatedDf, gasPhaseDf)
        resp2Csv = p.join(config['runtimeInfo']['madeByCharges']['chargeDir'], 'resp2_charges.csv')
        resp2Df.to_csv(resp2Csv)
        config['runtimeInfo']['madeByCharges']['chargesCsv'] = resp2Csv
    elif protocol == 'SOLVATOR':
        (chargesCsv, config) = partial_charges_SOLVATOR_protocol(outDir=config['runtimeInfo']['madeByCharges']['chargeDir'], config=config, debug=debug)
        config['runtimeInfo']['madeByCharges']['chargesCsv'] = chargesCsv
    config = Charged_Assistant.process_charge_csv(config)
    cleaner.clean_up_charges(config)
    config['checkpointInfo']['chargesComplete'] = True
    print('CHARGES DONE!')
    return config

def partial_charges_solvator_protocol(outDir: DirectoryPath, config: dict, debug: bool) -> tuple[pd.DataFrame, FilePath]:
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
    qmmmSinglepointDir = config['runtimeInfo']['madeByCharges']['qmmmSinglepointDir']
    fittingDir = config['runtimeInfo']['madeByCharges']['fittingDir']
    config = add_solvation_shell_with_SOLVATOR(config=config, debug=debug)
    config = tip3p_qmmm_protocol(config=config, debug=debug)
    conformerListTxt = Charged_Assistant.generate_conformer_list_file(qmmmSinglepointDir, fittingDir)
    chargeConstraintsTxt = Charged_Assistant.generate_charge_constraints_file(config, fittingDir)
    symmetryConstraintsTxt = Charged_Assistant.generate_symmetry_constraints_file(config, fittingDir)
    rawMultiWfnOutputs = Charged_Monster.run_charge_fitting(config=config, conformerListTxt=conformerListTxt, chargeConstraintsTxt=chargeConstraintsTxt, symmetryConstraintsTxt=symmetryConstraintsTxt, fittingDir=fittingDir)
    chargesDf = Charged_Monster.parse_multiwfn_output(rawMultiWfnOutputs)
    chargesCsv = p.join(fittingDir, 'charges.csv')
    chargesDf.to_csv(chargesCsv)
    return (chargesCsv, config)

def tip3p_qmmm_protocol(config, debug):
    """
    Runs QM/MM calculations for charge calculations using TIP3P water molecules:
    1. Optimisation
    2. Single-Point Energy calculation
    3. converts output to molden files for MultiWFN
    """
    conformerXyzs = config['runtimeInfo']['madeByConformers']['conformerXyzs']
    chargeDir = config['runtimeInfo']['madeByCharges']['chargeDir']
    qmmmParameterDir = p.join(chargeDir, '00_qm-mm_parameters')
    os.makedirs(qmmmParameterDir, exist_ok=True)
    config['runtimeInfo']['madeByCharges']['qmmmParameterDir'] = qmmmParameterDir
    unsolvatedXyz = p.join(qmmmParameterDir, 'unsolvated.xyz')
    copy(conformerXyzs[0], unsolvatedXyz)
    config = Charged_Assistant.get_qm_atoms_for_solvated_system(unsolvatedXyz, config)
    config = Charged_Monster.create_orca_ff_parameters(unsolvatedXyz, config)
    qmmmOptXyzs = qmmm_opt_protocol_for_charges(config=config, debug=debug)
    config['runtimeInfo']['madeByCharges']['solvatedOptimisedXyzs'] = qmmmOptXyzs
    qmmm_singlepoint_protocol_for_charges(qmmmOptXyzs=qmmmOptXyzs, config=config, debug=debug)
    return config

@Timer.time_function('QM/MM Optimisations', 'CHARGE_CALCULATIONS')
def qmmm_opt_protocol_for_charges(config, debug):
    solvatedXyzs = config['runtimeInfo']['madeByCharges']['solvatedXyzs']
    qmmmOptDir = config['runtimeInfo']['madeByCharges']['qmmmOptDir']
    qmmmOptArgsList = [(solvatedXyz, qmmmOptDir, config) for solvatedXyz in solvatedXyzs]
    if debug:
        qmmmOptXyzs = []
        for qmmmOptArgs in qmmmOptArgsList:
            qmmmOptXyz = Charged_Monster.run_qmmm_opt(qmmmOptArgs)
            qmmmOptXyzs.append(qmmmOptXyz)
    else:
        redText = '\x1b[31m'
        resetTextColor = '\x1b[0m'
        tqdmBarOptions = {'desc': f'{redText}Running QMMM Optimisations with TIP3P waters{resetTextColor}', 'ascii': '-🗲→', 'colour': 'yellow', 'unit': 'scan', 'dynamic_ncols': True}
        nCoresPerCalculation = config['chargeFittingInfo']['nCoresPerCalculation']
        availableCpus = config['miscInfo']['availableCpus']
        nCpus = min(availableCpus // nCoresPerCalculation, len(qmmmOptArgsList))
        with WorkerPool(n_jobs=nCpus) as pool:
            qmmmOptXyzs = pool.map(Charged_Monster.run_qmmm_opt, make_single_arguments(qmmmOptArgsList), progress_bar=True, iterable_len=len(qmmmOptArgsList), progress_bar_options=tqdmBarOptions)
    return qmmmOptXyzs

@Timer.time_function('QM/MM Single-Points', 'CHARGE_CALCULATIONS')
def qmmm_singlepoint_protocol_for_charges(qmmmOptXyzs, config, debug):
    qmmmSinglepointDir = config['runtimeInfo']['madeByCharges']['qmmmSinglepointDir']
    fittingDir = config['runtimeInfo']['madeByCharges']['fittingDir']
    qmmmSinglepointArgsList = [(qmmmOptXyz, qmmmSinglepointDir, fittingDir, config) for qmmmOptXyz in qmmmOptXyzs]
    if debug:
        for qmmmSinglepointArgs in qmmmSinglepointArgsList:
            Charged_Monster.run_qmmm_singlepoint(qmmmSinglepointArgs)
    else:
        redText = '\x1b[31m'
        resetTextColor = '\x1b[0m'
        tqdmBarOptions = {'desc': f'{redText}Running QMMM Single-Points with TIP3P waters{resetTextColor}', 'ascii': '-🗲→', 'colour': 'yellow', 'unit': 'scan', 'dynamic_ncols': True}
        nCoresPerCalculation = config['chargeFittingInfo']['nCoresPerCalculation']
        availableCpus = config['miscInfo']['availableCpus']
        os.environ['OMP_NUM_THREADS'] = str(nCoresPerCalculation)
        nCpus = min(availableCpus // nCoresPerCalculation, len(qmmmSinglepointArgsList))
        with WorkerPool(n_jobs=nCpus) as pool:
            optXyzs = pool.map(Charged_Monster.run_qmmm_singlepoint, make_single_arguments(qmmmSinglepointArgsList), progress_bar=True, iterable_len=len(qmmmSinglepointArgsList), progress_bar_options=tqdmBarOptions)

@Timer.time_function('Explicit Solvation', 'CHARGE_CALCULATIONS')
def add_solvation_shell_with_solvator(config: dict, debug: bool=False) -> dict:
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
    solvatorDir = config['runtimeInfo']['madeByCharges']['solvatorDir']
    moleculeInfo = config['moleculeInfo']
    conformerXyzs = config['runtimeInfo']['madeByConformers']['conformerXyzs']
    cappedPdb = config['runtimeInfo']['madeByCapping']['cappedPdb']
    nConformers = config['chargeFittingInfo']['nConformers']
    chargeFittingInfo = config['chargeFittingInfo']
    waterDensity = config['chargeFittingInfo']['waterDensity']
    nWaters = Charged_Assistant.how_many_waters_for_solvator(conformerXyzs, cappedPdb, solvatorDir, nWatersPerNmSquared=waterDensity)
    if nConformers == -1 or nConformers > len(conformerXyzs):
        sampledConformerXyzs = conformerXyzs
        config['chargeFittingInfo']['nConformers'] = len(conformerXyzs)
    else:
        sampledConformerXyzs = conformerXyzs[:nConformers]
    solvatorArgList = [(conformerXyz, solvatorDir, chargeFittingInfo, moleculeInfo, nWaters, config) for conformerXyz in sampledConformerXyzs]
    solvatedXyzs = []
    if debug:
        for arg in solvatorArgList:
            solvatedXyz = Charged_Monster.run_orca_solvator_for_charge_calculations(arg)
            solvatedXyzs.append(solvatedXyz)
    else:
        redText = '\x1b[31m'
        resetTextColor = '\x1b[0m'
        tqdmBarOptions = {'desc': f'{redText}Running SOLVATOR for Charge Calculations{resetTextColor}', 'ascii': '-🗲→', 'colour': 'yellow', 'unit': 'scan', 'dynamic_ncols': True}
        availableCpus = config['miscInfo']['availableCpus']
        nCoresPerCalculation = chargeFittingInfo['nCoresPerCalculation']
        nCpus = min(availableCpus, len(solvatorArgList))
        with WorkerPool(n_jobs=nCpus) as pool:
            solvatedXyzs = pool.map(Charged_Monster.run_orca_solvator_for_charge_calculations, make_single_arguments(solvatorArgList), progress_bar=True, iterable_len=len(solvatorArgList), progress_bar_options=tqdmBarOptions)
    config['runtimeInfo']['madeByCharges']['solvatedXyzs'] = solvatedXyzs
    config['runtimeInfo']['madeByCharges']['nWaters'] = nWaters
    return config

def partial_charge_resp_protocol(outDir: DirectoryPath, useSolvation: bool, config: dict, debug: bool=False) -> tuple[pd.DataFrame, FilePath]:
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
    orcaDir = p.join(outDir, '01_orca_calculations')
    os.makedirs(orcaDir, exist_ok=True)
    fittingDir = p.join(outDir, '02_charge_fitting')
    os.makedirs(fittingDir, exist_ok=True)
    conformerXyzs = config['runtimeInfo']['madeByConformers']['conformerXyzs']
    run_qm_calculations_for_RESP(conformerXyzs=conformerXyzs, orcaDir=orcaDir, fittingDir=fittingDir, config=config, useSolvation=useSolvation, debug=debug)
    conformerListTxt = Charged_Assistant.generate_conformer_list_file(orcaDir, fittingDir)
    chargeConstraintsTxt = Charged_Assistant.generate_charge_constraints_file(config, fittingDir)
    symmetryConstraintsTxt = Charged_Assistant.generate_symmetry_constraints_file(config, fittingDir)
    rawMultiWfnOutputs = Charged_Monster.run_charge_fitting(config=config, conformerListTxt=conformerListTxt, chargeConstraintsTxt=chargeConstraintsTxt, symmetryConstraintsTxt=symmetryConstraintsTxt, fittingDir=fittingDir)
    chargesDf = Charged_Monster.parse_multiwfn_output(rawMultiWfnOutputs)
    chargesCsv = p.join(fittingDir, 'charges.csv')
    chargesDf.to_csv(chargesCsv)
    return (chargesDf, chargesCsv)

@Timer.time_function('QM Calculations', 'CHARGE_CALCULATIONS')
def run_qm_calculations_for_resp(conformerXyzs: list[FilePath], orcaDir: DirectoryPath, fittingDir: DirectoryPath, config: dict, useSolvation: bool=True, debug: bool=False) -> None:
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
    moleculeInfo = config['moleculeInfo']
    nConformers = config['chargeFittingInfo']['nConformers']
    chargeFittingInfo = config['chargeFittingInfo']
    if nConformers == -1 or nConformers > len(conformerXyzs):
        sampledConformerXyzs = conformerXyzs
        config['chargeFittingInfo']['nConformers'] = len(conformerXyzs)
    else:
        sampledConformerXyzs = conformerXyzs[:nConformers]
    argsList = [(conformerXyz, orcaDir, fittingDir, chargeFittingInfo, moleculeInfo, useSolvation, config) for conformerXyz in sampledConformerXyzs]
    redText = '\x1b[31m'
    resetTextColor = '\x1b[0m'
    tqdmBarOptions = {'desc': f'{redText}Running QM for Charge Calculations{resetTextColor}', 'ascii': '-🗲→', 'colour': 'yellow', 'unit': 'scan', 'dynamic_ncols': True}
    if debug:
        for arg in argsList:
            Charged_Monster.run_orca_singlepoint_for_charge_calculations(arg)
    else:
        nCpus = min(len(argsList), config['miscInfo']['availableCpus'])
        with WorkerPool(n_jobs=nCpus) as pool:
            pool.map(Charged_Monster.run_orca_singlepoint_for_charge_calculations, make_single_arguments(argsList), progress_bar=True, iterable_len=len(argsList), progress_bar_options=tqdmBarOptions)
if __name__ == '__main__':
    raise NotImplementedError