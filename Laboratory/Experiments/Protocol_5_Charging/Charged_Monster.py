from os import path as p
import os
import glob
from shutil import copy, move
from subprocess import call, DEVNULL
from pdbUtils import pdbUtils
import pandas as pd
import pexpect
from tqdm import tqdm
from copy import deepcopy
import mdtraj as md
from OperatingTools import drOrca, cleaner, Timer
from Experiments.Protocol_1_Capping.Capping_Assistant import find_bonded_atoms
from . import Charged_Assistant

class FilePath:
    pass

class DirectoryPath:
    pass

def run_qmmm_opt(qmmmOptArgs):
    (solvatedXyz, singlepointDir, config) = qmmmOptArgs
    qmAtoms = config['runtimeInfo']['madeByCharges']['qmAtoms']
    solvatedParams = config['runtimeInfo']['madeByCharges']['solvatedParams']
    nCoresPerCalculation = config['chargeFittingInfo']['nCoresPerCalculation']
    conformerName = p.basename(solvatedXyz).split('.')[0]
    conformerQmmmOptDir = p.join(singlepointDir, f'{conformerName}')
    os.makedirs(conformerQmmmOptDir, exist_ok=True)
    qmmmOptOrcaInput = drOrca.make_orca_input_qmmm_opt(inputXyz=solvatedXyz, outDir=conformerQmmmOptDir, moleculeInfo=config['moleculeInfo'], qmMethod=config['chargeFittingInfo']['optMethod'], qmAtoms=qmAtoms, parameterFile=solvatedParams, nCpus=nCoresPerCalculation)
    qmmmOptOrcaOutput = p.join(conformerQmmmOptDir, 'QMMM_orca_opt.out')
    solvatedOptXyz = p.join(conformerQmmmOptDir, f'{conformerName}.xyz')
    if not p.isfile(qmmmOptOrcaOutput):
        drOrca.run_orca(qmmmOptOrcaInput, qmmmOptOrcaOutput, config)
        optXyz = p.join(conformerQmmmOptDir, 'QMMM_orca_opt.xyz')
        copy(optXyz, solvatedOptXyz)
        cleaner.clean_up_opt_dir(conformerQmmmOptDir, config)
    return solvatedOptXyz

def run_qmmm_singlepoint(qmmmSinglepointArgs):
    (solvatedOptXyz, singlepointDir, fittingDir, config) = qmmmSinglepointArgs
    qmAtoms = config['runtimeInfo']['madeByCharges']['qmAtoms']
    solvatedParams = config['runtimeInfo']['madeByCharges']['solvatedParams']
    conformerName = p.basename(solvatedOptXyz).split('.')[0]
    conformerQmmmSinglepointDir = p.join(singlepointDir, f'{conformerName}')
    os.makedirs(conformerQmmmSinglepointDir, exist_ok=True)
    qmmmSinglepointOrcaInput = drOrca.make_orca_input_qmmm_singlepoint(inputXyz=solvatedOptXyz, outDir=conformerQmmmSinglepointDir, moleculeInfo=config['moleculeInfo'], qmMethod=config['chargeFittingInfo']['singlePointMethod'], qmAtoms=qmAtoms, parameterFile=solvatedParams)
    qmmmSinglepointOrcaOutput = p.join(conformerQmmmSinglepointDir, 'QMMM_orca_sp.out')
    if not p.isfile(qmmmSinglepointOrcaOutput):
        drOrca.run_orca(qmmmSinglepointOrcaInput, qmmmSinglepointOrcaOutput, config)
    call(['orca_2mkl', p.join(conformerQmmmSinglepointDir, 'QMMM_orca_sp'), '-molden'], stdout=DEVNULL)
    singlePointMolden = p.join(conformerQmmmSinglepointDir, 'QMMM_orca_sp.molden.input')
    destMolden = p.join(fittingDir, f'{conformerName}.molden.input')
    copy(singlePointMolden, destMolden)

def create_orca_ff_parameters(unsolvatedXyz: FilePath, config: dict) -> dict:
    """
    Creates an ORCAFF.prms file for a solvated XYZ file
    1. Use orca_mm -makeff to make a forcefield file for just the unsolvated molecule
    2. Use orca_mm -repeatff to repeat a pre-made TIP3P water parameter * nWaters
    3. Ise orca_mm -mergeff to combine the results of steps 1 and 2   

    Args:
        unsolvatedXyz (FilePath): Path to a conformerXYZ moved to the qmmmParameters dir
        config (dict): Contains all run information
    Returns:
        config (dict): updated config
    """
    qmmmParameterDir = config['runtimeInfo']['madeByCharges']['qmmmParameterDir']
    nWaters = config['runtimeInfo']['madeByCharges']['nWaters']
    moleculeName = config['moleculeInfo']['moleculeName']
    makeffCommand = ['orca_mm', '-makeff', unsolvatedXyz]
    call(makeffCommand)
    unsolvatedParams = p.join(qmmmParameterDir, 'unsolvated.ORCAFF.prms')
    tip3pParamsSource = Charged_Assistant.find_tip3p_water_params()
    tip3pParams = p.join(qmmmParameterDir, 'TIP3P.ORCAFF.prms')
    copy(tip3pParamsSource, tip3pParams)
    repeatffCommand = ['orca_mm', '-repeatff', tip3pParams, str(nWaters)]
    call(repeatffCommand)
    watersParams = p.join(qmmmParameterDir, f'TIP3P_repeat{nWaters}.ORCAFF.prms')
    mergeffCommand = ['orca_mm', '-mergeff', unsolvatedParams, watersParams]
    call(mergeffCommand, stdout=DEVNULL)
    outParams = p.join(qmmmParameterDir, f'unsolvated_merged.ORCAFF.prms')
    solvatedParams = p.join(qmmmParameterDir, f'{moleculeName}_solvated.ORCAFF.prms')
    move(outParams, solvatedParams)
    config['runtimeInfo']['madeByCharges']['solvatedParams'] = solvatedParams
    return config

def run_orca_solvator_for_charge_calculations(args: tuple) -> FilePath:
    """
    Runs ORCA's SOLVATOR protocol to place explicit water molecules around the conformer of interest

    This is much faster than using the non-STOCHASTIC cluster mode in the solvator

    Args:
        args (tuple): to be unpacked, containing:
            TODO
    Returns:
        solvatedXyz (FilePath): XYZ file containing the solvated conformer
    
    """
    (conformerXyz, solvatorDir, chargeFittingInfo, moleculeInfo, nWaters, config) = args
    conformerName = p.basename(conformerXyz).split('.')[0]
    conformerSolvatorDir = p.join(solvatorDir, conformerName)
    os.makedirs(conformerSolvatorDir, exist_ok=True)
    solvatorOrcaInput = drOrca.make_orca_input_for_solvator(inputXyz=conformerXyz, outDir=conformerSolvatorDir, moleculeInfo=moleculeInfo, qmMethod=chargeFittingInfo['optMethod'], solvationMethod=chargeFittingInfo['optSolvationMethod'], nWaters=nWaters)
    solvatorOrcaOutput = p.join(conformerSolvatorDir, 'SOLVATOR_orca.out')
    solvatedXyz = p.join(conformerSolvatorDir, f'{conformerName}.xyz')
    if not p.isfile(solvatorOrcaOutput):
        drOrca.run_orca(solvatorOrcaInput, solvatorOrcaOutput, config)
        outXyz = p.join(conformerSolvatorDir, 'SOLVATOR_orca.solvator.xyz')
        os.rename(outXyz, solvatedXyz)
        cleaner.clean_up_solvator_dir(conformerSolvatorDir, config)
    return solvatedXyz

def run_orca_singlepoint_for_charge_calculations(args) -> None:
    (conformerXyz, orcaDir, fittingDir, chargeFittingInfo, moleculeInfo, useSolvation, config) = args
    conformerName = p.basename(conformerXyz).split('.')[0]
    conformerQmDir = p.join(orcaDir, conformerName)
    os.makedirs(conformerQmDir, exist_ok=True)
    if useSolvation:
        optSolvation = chargeFittingInfo['optSolvationMethod']
        singlePointSolvation = chargeFittingInfo['singlePointSolvationMethod']
    else:
        optSolvation = None
        singlePointSolvation = None
    orcaOptInput = drOrca.make_orca_input_for_opt(inputXyz=conformerXyz, outDir=conformerQmDir, moleculeInfo=moleculeInfo, qmMethod=chargeFittingInfo['optMethod'], solvationMethod=optSolvation)
    orcaOptOutput = p.join(conformerQmDir, 'orca_opt.out')
    if not p.isfile(orcaOptOutput):
        drOrca.run_orca(orcaOptInput, orcaOptOutput, config)
    optXyz = p.join(conformerQmDir, 'orca_opt.xyz')
    orcaSinglePointInput = drOrca.make_orca_input_for_singlepoint(inputXyz=optXyz, outDir=conformerQmDir, moleculeInfo=moleculeInfo, qmMethod=chargeFittingInfo['singlePointMethod'], solvationMethod=singlePointSolvation)
    orcaSinglePointOutput = p.join(conformerQmDir, 'orca_sp.out')
    if not p.isfile(orcaSinglePointOutput):
        drOrca.run_orca(orcaSinglePointInput, orcaSinglePointOutput, config)
    singlePointName = p.splitext(p.basename(orcaSinglePointInput))[0]
    call(['orca_2mkl', p.join(conformerQmDir, singlePointName), '-molden'], stdout=DEVNULL)
    singlePointMolden = p.join(conformerQmDir, 'orca_sp.molden.input')
    destMolden = p.join(fittingDir, f'{conformerName}.molden.input')
    cleaner.clean_up_singlepoint_dir(conformerQmDir, config, keepGbw=True)
    copy(singlePointMolden, destMolden)

def parse_multiwfn_output(rawMultiWfnOutputs: FilePath) -> pd.DataFrame:
    """
    Parses the output from Multiwfn to extract atomic charges.

    Args:
        rawMultiwfnOutputs (FilePath): Path to the file containing the raw Multiwfn output.

    Returns:
        chargesDf (pd.DataFrame) : A DataFrame containing columns 'atomIndex', 'atomElement', and 'Charge', 
                      representing the index, element, and charge of each atom, respectively.

    This function reads the Multiwfn output file, identifies the section after the 
    "Successfully converged!" line, and extracts atomic charges until the line containing 
    "Sum of charges:". The extracted data is structured into a DataFrame for further analysis.
    """
    with open(rawMultiWfnOutputs, 'r') as file:
        lines = file.readlines()
    startIndex = next((i for (i, line) in enumerate(lines) if 'Successfully converged!' in line)) + 3
    chargeDataLines = []
    for line in lines[startIndex:]:
        if 'Sum of charges:' in line.strip():
            break
        chargeDataLines.append(line.strip())
    data = []
    for line in chargeDataLines:
        parts = line.split(')')
        charge = float(parts[1].strip())
        parts = parts[0].split('(')
        atomIndex = int(parts[0])
        atomElement = parts[1].strip()
        data.append({'atomIndex': atomIndex, 'atomElement': atomElement, 'Charge': charge})
    chargesDf = pd.DataFrame(data)
    return chargesDf

def get_charge_group_indexes(pdbFile: FilePath, config: dict) -> dict:
    """
    Uses user-defined charge groups to construct a complete dictionary of atom indexes for each charge group
    Users do not need to specify hydrogen atoms, this function will find them automatically
    All non-specified heavy atoms and associated protons are assigned to a "left-overs" charge group
    Data is stored in the config dict and returned

    Args:
        pdbFile (FilePath): path to pdb file
        config (dict): drFrankenstein config

    Returns:
        config (dict): updated config
    """
    pdbDf = pdbUtils.pdb2df(pdbFile)
    chargeGroupsInput = config['moleculeInfo']['chargeGroups']
    chargeGroups = deepcopy(chargeGroupsInput)
    overallCharge = config['moleculeInfo']['charge']
    uncappedDf = pdbDf[~pdbDf['RES_NAME'].isin(['NME', 'ACE'])]
    userDefinedAtoms = []
    userDefinedCharge = 0
    for (chargeGroupName, chargeGroupData) in chargeGroupsInput.items():
        heavyChargeGroupAtoms = chargeGroupData['atoms']
        protonNames = []
        for heavyAtom in heavyChargeGroupAtoms:
            boundProtons = find_bonded_atoms(pdbDf, heavyAtom)
            protonNames.extend(boundProtons)
        protonNames = [atomName for atomName in protonNames if atomName.startswith('H')]
        chargeGroupAtoms = heavyChargeGroupAtoms + protonNames
        chargeGroupIndexes = uncappedDf[uncappedDf['ATOM_NAME'].isin(chargeGroupAtoms)]['ATOM_ID'].to_list()
        chargeGroups[chargeGroupName]['indexes'] = chargeGroupIndexes
        userDefinedCharge += chargeGroupData['charge']
        userDefinedAtoms.extend(chargeGroupAtoms)
    leftOverDf = uncappedDf[~uncappedDf['ATOM_NAME'].isin(userDefinedAtoms)]
    leftOverAtoms = leftOverDf['ATOM_NAME'].to_list()
    leftOverIndexes = leftOverDf['ATOM_ID'].to_list()
    leftOverCharge = overallCharge - userDefinedCharge
    if len(leftOverAtoms) > 0:
        chargeGroups['left-over-atoms'] = {'atoms': leftOverAtoms, 'charge': leftOverCharge, 'indexes': leftOverIndexes}
    for capResName in ['NME', 'ACE']:
        capDf = pdbDf[pdbDf['RES_NAME'] == capResName]
        for (capIndex, capResDf) in capDf.groupby('RES_ID'):
            chargeGroupAtoms = capResDf['ATOM_NAME'].to_list()
            chargeGroupIndexes = capResDf['ATOM_ID'].to_list()
            chargeGroupName = f'{capResName}_cap_{capIndex}'
            chargeGroups[chargeGroupName] = {'atoms': chargeGroupAtoms, 'charge': 0, 'indexes': chargeGroupIndexes}
    config['runtimeInfo']['madeByCharges']['chargeGroups'] = chargeGroups
    return config

@Timer.time_function('Charge Fitting', 'CHARGE_CALCULATIONS')
def run_charge_fitting(config: dict, conformerListTxt: FilePath, chargeConstraintsTxt: FilePath, symmetryConstraintsTxt: FilePath, fittingDir: DirectoryPath) -> dict:
    """
    Runs MultiWFN to calculate charges for a set of conformers
    Usually, MultiWFN_noGUI uses a silly command-line input style. 
    This function uses pexpext to automate this for RESP charge fitting

    Args:
        config (dict): drFrankenstein config
        conformerListTxt (FilePath): path to conformer list file
        chargeConstraintsTxt (FilePath): path to charge constraints file
        fittingDir (DirectoryPath): path to fitting directory
    Returns:
        config (dict): updated config
    """
    multiWfnDir = config['pathInfo']['multiWfnDir']
    nConformers = config['chargeFittingInfo']['nConformers']
    os.chdir(multiWfnDir)
    multiWfnExe = p.join(multiWfnDir, 'Multiwfn_noGUI')
    if not p.isfile(multiWfnExe):
        raise FileNotFoundError('Multiwfn_noGUI not found')
    moldenFile = glob.glob(p.join(fittingDir, '*.molden.input'))[0]
    outputFile = p.join(fittingDir, 'MultiWfn_raw_outputs.txt')
    with open(outputFile, 'w') as logFile:
        logFile.write('Starting charge fitting process...\n')
        logFile.flush()
        chargeFittingCommand = f'{multiWfnExe} {moldenFile}'
        child = pexpect.spawn(chargeFittingCommand, encoding='utf-8', logfile=logFile)
        child.expect('.*300.*\\n?\\r?')
        child.sendline('7')
        child.expect('.*20.*\\n?\\r?')
        child.sendline('18')
        child.expect('.*11.*\\n?\\r?')
        child.sendline('-1')
        child.expect('.*Input.*\\n?\\r?')
        child.sendline(conformerListTxt)
        child.expect('.*11.*\\n?\\r?')
        child.sendline('6')
        child.expect('.*1.*\\n?\\r?')
        child.sendline('1')
        child.expect('.*Input.*\\n?\\r?')
        child.sendline(chargeConstraintsTxt)
        child.expect('.*11.*\\n?\\r?')
        child.sendline('5')
        child.expect('.*1.*\\n?\\r?')
        child.sendline('1')
        child.expect('.*Input.*\\n?\\r?')
        child.sendline(symmetryConstraintsTxt)
        child.expect('.*11.*\\n?\\r?')
        child.sendline('2')
        orangeText = '\x1b[38;5;172m'
        resetTextColor = '\x1b[0m'
        tqdmBarOptions = {'desc': f'{orangeText}Performing Charge Fitting{resetTextColor}', 'ascii': '-ϟ→', 'colour': 'yellow', 'unit': 'calculations', 'dynamic_ncols': True, 'leave': True}
        totalTickProgress = nConformers * 100.0
        progressBar = tqdm(total=totalTickProgress, **tqdmBarOptions)
        try:
            while True:
                index = child.expect(['Progress: \\[.*?\\]\\s+(\\d+\\.\\d+) %', ' 11 Choose ESP type, current: Nuclear \\+ Electronic[\\n\\r]?', '.*\\(y/n\\)[\\n\\r]?', pexpect.TIMEOUT], timeout=5)
                if index == 0:
                    progressBar.update(1)
                elif index in [1, 2]:
                    break
                elif index == 3:
                    break
        except pexpect.ExceptionPexpect as e:
            raise e
        child.close()
    return outputFile