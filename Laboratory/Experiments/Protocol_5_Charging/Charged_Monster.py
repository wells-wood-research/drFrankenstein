
## BASIC IMPORTS ##
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
## drFrankenstein LIBRARIES ##
from OperatingTools import drOrca, cleaner, Timer
from Experiments.Protocol_1_Capping.Capping_Assistant import find_bonded_atoms
from . import Charged_Assistant
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_qmmm_opt(qmmmOptArgs):
    ## unpack args
    solvatedXyz, singlepointDir, config = qmmmOptArgs

    ## unpack config
    qmAtoms = config["runtimeInfo"]["madeByCharges"]["qmAtoms"]
    solvatedParams = config["runtimeInfo"]["madeByCharges"]["solvatedParams"]
    nCoresPerCalculation = config["chargeFittingInfo"]["nCoresPerCalculation"]
    ## get conformer name, make a dir for qmmm opt calculations
    conformerName = p.basename(solvatedXyz).split(".")[0]
    
    conformerQmmmOptDir = p.join(singlepointDir, f"{conformerName}")
    os.makedirs(conformerQmmmOptDir, exist_ok=True)

    qmmmOptOrcaInput = drOrca.make_orca_input_qmmm_opt(inputXyz = solvatedXyz,
                                                    outDir = conformerQmmmOptDir,
                                                    moleculeInfo = config["moleculeInfo"],
                                                    qmMethod = config["chargeFittingInfo"]["optMethod"],
                                                    qmAtoms=qmAtoms,
                                                    parameterFile=solvatedParams,
                                                    nCpus=nCoresPerCalculation)
    qmmmOptOrcaOutput = p.join(conformerQmmmOptDir, "QMMM_orca_opt.out")
    solvatedOptXyz = p.join(conformerQmmmOptDir, f"{conformerName}.xyz")

    if not p.isfile(qmmmOptOrcaOutput):
        drOrca.run_orca(qmmmOptOrcaInput, qmmmOptOrcaOutput, config)
        optXyz = p.join(conformerQmmmOptDir, "QMMM_orca_opt.xyz")
        copy(optXyz, solvatedOptXyz)
        cleaner.clean_up_opt_dir(conformerQmmmOptDir, config)

    return solvatedOptXyz

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_qmmm_singlepoint(qmmmSinglepointArgs):

    ## unpack args
    solvatedOptXyz, singlepointDir, fittingDir, config = qmmmSinglepointArgs 
    ## unpack config
    qmAtoms = config["runtimeInfo"]["madeByCharges"]["qmAtoms"]
    solvatedParams = config["runtimeInfo"]["madeByCharges"]["solvatedParams"]

    ## get conformer name, make a dir for qmmm opt calculations
    conformerName = p.basename(solvatedOptXyz).split(".")[0]
    
    conformerQmmmSinglepointDir = p.join(singlepointDir, f"{conformerName}")
    os.makedirs(conformerQmmmSinglepointDir, exist_ok=True)

    qmmmSinglepointOrcaInput = drOrca.make_orca_input_qmmm_singlepoint(inputXyz=solvatedOptXyz,
                                                                       outDir= conformerQmmmSinglepointDir,
                                                                       moleculeInfo= config["moleculeInfo"],
                                                                       qmMethod= config["chargeFittingInfo"]["singlePointMethod"],
                                                                       qmAtoms=qmAtoms,
                                                                       parameterFile=solvatedParams)
    qmmmSinglepointOrcaOutput = p.join(conformerQmmmSinglepointDir, "QMMM_orca_sp.out")
    if not p.isfile(qmmmSinglepointOrcaOutput):
        drOrca.run_orca(qmmmSinglepointOrcaInput, qmmmSinglepointOrcaOutput, config)
    
    ## create molden file for charge fitting TODO: move to its own func
    call(["orca_2mkl", p.join(conformerQmmmSinglepointDir,"QMMM_orca_sp"), "-molden"], stdout=DEVNULL)
    singlePointMolden = p.join(conformerQmmmSinglepointDir, "QMMM_orca_sp.molden.input")

    destMolden = p.join(fittingDir, f"{conformerName}.molden.input")

    copy(singlePointMolden, destMolden)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def create_orca_ff_parameters(unsolvatedXyz: FilePath, config:dict) -> dict:
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

    ## unpack config
    qmmmParameterDir = config["runtimeInfo"]["madeByCharges"]["qmmmParameterDir"]
    nWaters  = config["runtimeInfo"]["madeByCharges"]["nWaters"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    ## Create a dummy parameter set for the molecule without solvating waters
    makeffCommand = ["orca_mm", "-makeff", unsolvatedXyz]
    call(makeffCommand)
    unsolvatedParams = p.join(qmmmParameterDir, "unsolvated.ORCAFF.prms")


    ## Create a parameter file for nWaters TIP3P water molecules
    tip3pParamsSource = Charged_Assistant.find_tip3p_water_params()
    tip3pParams = p.join(qmmmParameterDir, "TIP3P.ORCAFF.prms")
    copy(tip3pParamsSource, tip3pParams)
    repeatffCommand = ["orca_mm", "-repeatff", tip3pParams, str(nWaters)]
    call(repeatffCommand)
    watersParams = p.join(qmmmParameterDir, f"TIP3P_repeat{nWaters}.ORCAFF.prms")

    ## Combine unsolvated parameters with water parameters
    mergeffCommand = ["orca_mm", "-mergeff", unsolvatedParams, watersParams]
    call(mergeffCommand, stdout=DEVNULL)
    outParams = p.join(qmmmParameterDir, f"unsolvated_merged.ORCAFF.prms")
    solvatedParams = p.join(qmmmParameterDir, f"{moleculeName}_solvated.ORCAFF.prms")

    move(outParams, solvatedParams)

    config["runtimeInfo"]["madeByCharges"]["solvatedParams"] = solvatedParams

    return config


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
    ## unpack config
    conformerXyz, solvatorDir, chargeFittingInfo,  moleculeInfo, nWaters, config = args
    ## get conformer name, make a dir for solvator + opt calculations
    conformerName = p.basename(conformerXyz).split(".")[0]
    
    conformerSolvatorDir = p.join(solvatorDir, conformerName)
    os.makedirs(conformerSolvatorDir, exist_ok=True)

    solvatorOrcaInput = drOrca.make_orca_input_for_solvator(inputXyz = conformerXyz,
                                                    outDir = conformerSolvatorDir,
                                                    moleculeInfo = moleculeInfo,
                                                    qmMethod = chargeFittingInfo["optMethod"],
                                                    solvationMethod = chargeFittingInfo["optSolvationMethod"],
                                                    nWaters= nWaters)
    solvatorOrcaOutput = p.join(conformerSolvatorDir, "SOLVATOR_orca.out")
    solvatedXyz = p.join(conformerSolvatorDir, f"{conformerName}.xyz")

    if not p.isfile(solvatorOrcaOutput):
        drOrca.run_orca(solvatorOrcaInput, solvatorOrcaOutput, config)
        outXyz = p.join(conformerSolvatorDir, "SOLVATOR_orca.solvator.xyz")
        ## rename the xyz file  
        os.rename(outXyz, solvatedXyz)

        cleaner.clean_up_solvator_dir(conformerSolvatorDir, config)

    return solvatedXyz

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_orca_singlepoint_for_charge_calculations(args) -> None: 
    conformerXyz, orcaDir, fittingDir, chargeFittingInfo,  moleculeInfo, useSolvation, config = args
    conformerName = p.basename(conformerXyz).split(".")[0]
    
    conformerQmDir = p.join(orcaDir, conformerName)
    os.makedirs(conformerQmDir, exist_ok=True)

    if useSolvation:
        optSolvation = chargeFittingInfo["optSolvationMethod"]
        singlePointSolvation = chargeFittingInfo["singlePointSolvationMethod"]
    else: 
        optSolvation = None
        singlePointSolvation = None

    orcaOptInput = drOrca.make_orca_input_for_opt(inputXyz = conformerXyz,
                                                    outDir = conformerQmDir,
                                                    moleculeInfo = moleculeInfo,
                                                    qmMethod = chargeFittingInfo["optMethod"],
                                                    solvationMethod = optSolvation)
    
    orcaOptOutput = p.join(conformerQmDir, "orca_opt.out")
    if not p.isfile(orcaOptOutput):
        drOrca.run_orca(orcaOptInput, orcaOptOutput, config)

    optXyz = p.join(conformerQmDir, "orca_opt.xyz")
    orcaSinglePointInput = drOrca.make_orca_input_for_singlepoint(inputXyz = optXyz,
                                                                outDir = conformerQmDir,
                                                                moleculeInfo = moleculeInfo,
                                                                qmMethod = chargeFittingInfo["singlePointMethod"],
                                                                solvationMethod = singlePointSolvation)
    
    orcaSinglePointOutput = p.join(conformerQmDir, "orca_sp.out")
    if not p.isfile(orcaSinglePointOutput):
        drOrca.run_orca(orcaSinglePointInput, orcaSinglePointOutput, config)

    singlePointName = p.splitext(p.basename(orcaSinglePointInput))[0]
    call(["orca_2mkl", p.join(conformerQmDir,singlePointName), "-molden"], stdout=DEVNULL)
    singlePointMolden = p.join(conformerQmDir, "orca_sp.molden.input")

    destMolden = p.join(fittingDir, f"{conformerName}.molden.input")

    cleaner.clean_up_singlepoint_dir(conformerQmDir, config, keepGbw=True)

    copy(singlePointMolden, destMolden)
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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

    # Find the index of the line after "Successfully converged!"
    start_index = next(i for i, line in enumerate(lines) 
                       if "Successfully converged!" in line) + 3

    # Extract the charge data lines
    chargeDataLines = []
    for line in lines[start_index:]:
        if "Sum of charges:" in line.strip():
            break
        chargeDataLines.append(line.strip())

    # Parse the charge data into a DataFrame
    data = []
    for line in chargeDataLines:
        parts = line.split(")")
        charge = float(parts[1].strip())
        parts = parts[0].split("(")
        atomIndex = int(parts[0])
        atomElement = parts[1].strip()

        data.append({'atomIndex': atomIndex, 'atomElement': atomElement, 'Charge': charge})

    chargesDf = pd.DataFrame(data)
    return chargesDf
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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

    chargeGroupsInput = config["moleculeInfo"]["chargeGroups"]
    chargeGroups = deepcopy(chargeGroupsInput)  # Use deepcopy to avoid modifying the original
    overallCharge = config["moleculeInfo"]["charge"]

    ## remove capping groups
    uncappedDf = pdbDf[~pdbDf["RES_NAME"].isin(["NME", "ACE"])]

    ## deal with user-defined charge groups
    userDefinedAtoms = []
    userDefinedCharge = 0
    for chargeGroupName, chargeGroupData in chargeGroupsInput.items():
        ## get heavy atom names
        heavyChargeGroupAtoms = chargeGroupData["atoms"]
        ## fill in hydrogen atom names
        protonNames = []
        for heavyAtom in heavyChargeGroupAtoms:
            boundProtons = find_bonded_atoms(pdbDf, heavyAtom)
            protonNames.extend(boundProtons)
        protonNames = [atomName for atomName in protonNames if atomName.startswith("H")] ##TODO: make this nicer
        chargeGroupAtoms = heavyChargeGroupAtoms + protonNames
        chargeGroupIndexes = uncappedDf[uncappedDf["ATOM_NAME"].isin(chargeGroupAtoms)]["ATOM_ID"].to_list()
        chargeGroups[chargeGroupName]["indexes"] = chargeGroupIndexes
        userDefinedCharge += chargeGroupData["charge"]
        userDefinedAtoms.extend(chargeGroupAtoms)

    ## deal with left-over atoms, these will all go in one group
    leftOverDf = uncappedDf[~uncappedDf["ATOM_NAME"].isin(userDefinedAtoms)]
    leftOverAtoms = leftOverDf["ATOM_NAME"].to_list()
    leftOverIndexes = leftOverDf["ATOM_ID"].to_list()
    leftOverCharge = overallCharge - userDefinedCharge 
    
    if len(leftOverAtoms) > 0:
        chargeGroups["left-over-atoms"] = {
            "atoms" : leftOverAtoms,
            "charge" : leftOverCharge,
            "indexes" : leftOverIndexes
        }
    ## deal with Terminal Caps 
    for capResName in ["NME", "ACE"]:
        capDf = pdbDf[pdbDf["RES_NAME"]==capResName]
        for capIndex, capResDf in capDf.groupby("RES_ID"):
            chargeGroupAtoms = capResDf["ATOM_NAME"].to_list()
            chargeGroupIndexes = capResDf["ATOM_ID"].to_list()
            chargeGroupName = f"{capResName}_cap_{capIndex}"
            chargeGroups[chargeGroupName] = {
                "atoms": chargeGroupAtoms,
                "charge": 0,
                "indexes": chargeGroupIndexes
            }

    config["runtimeInfo"]["madeByCharges"]["chargeGroups"] = chargeGroups

    return config
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Charge Fitting", "CHARGE_CALCULATIONS")
def run_charge_fitting(config: dict,
                        conformerListTxt: FilePath,
                          chargeConstraintsTxt: FilePath,
                          symmetryConstraintsTxt: FilePath,
                            fittingDir: DirectoryPath) -> dict:
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

    ## Unpack pathInfo
    multiWfnDir = config["pathInfo"]["multiWfnDir"]
    nConformers = config["chargeFittingInfo"]["nConformers"]

    ## Change dir to MultiWFN build dir so it can find settings.ini
    os.chdir(multiWfnDir)

    ## Get MultiWFN binary
    multiWfnExe = p.join(multiWfnDir, "Multiwfn_noGUI")
    if not p.isfile(multiWfnExe):
        raise FileNotFoundError("Multiwfn_noGUI not found")

    ## Get a molden file for the input command
    moldenFile = glob.glob(p.join(fittingDir, "*.molden.input"))[0]

    ## Define output file
    outputFile = p.join(fittingDir, "MultiWfn_raw_outputs.txt")

    ## Open a file object for continuous logging
    with open(outputFile, 'w') as logFile:
        # Clear file and write a starting message
        logFile.write("Starting charge fitting process...\n")
        logFile.flush()  # Ensure it’s written immediately

        ## Spawn the child process with logfile enabled
        chargeFittingCommand = f"{multiWfnExe} {moldenFile}"
        child = pexpect.spawn(chargeFittingCommand, encoding='utf-8', logfile=logFile)

        ####### GET TO RESP MAIN MENU #######
        child.expect(r'.*300.*\n?\r?')
        child.sendline("7")

        child.expect(r'.*20.*\n?\r?')
        child.sendline("18")

        ###### LOAD CONFORMERS #######
        child.expect(r'.*11.*\n?\r?')
        child.sendline("-1")

        child.expect(r'.*Input.*\n?\r?')
        child.sendline(conformerListTxt)

        ####### LOAD CHARGE CONSTRAINTS #######
        child.expect(r'.*11.*\n?\r?')
        child.sendline("6")

        child.expect(r'.*1.*\n?\r?')
        child.sendline("1")

        child.expect(r'.*Input.*\n?\r?')
        child.sendline(chargeConstraintsTxt)

        ####### LOAD SYMMETRY CONSTRAINTS #######
        child.expect(r'.*11.*\n?\r?')
        child.sendline("5")

        child.expect(r'.*1.*\n?\r?')
        child.sendline("1")

        child.expect(r'.*Input.*\n?\r?')
        child.sendline(symmetryConstraintsTxt)

        ####### RUN CALCULATIONS #######
        child.expect(r'.*11.*\n?\r?')
        child.sendline("2")

        # Monitor the output for progress updates
        orangeText = "\033[38;5;172m"
        resetTextColor = "\033[0m"
        # Initialize tqdm progress bar
        tqdmBarOptions = {
        "desc": f"{orangeText}Performing Charge Fitting{resetTextColor}",
        "ascii": "-ϟ→",  
        "colour": "yellow",
        "unit":  "calculations",
        "dynamic_ncols": True,
        "leave": True
            }


        totalTickProgress = nConformers * 100.00
        progress_bar = tqdm(total=totalTickProgress, **tqdmBarOptions)
        try:
            while True:
                # Expect either a progress line or the end of the process (EOF)
                index = child.expect([r"Progress: \[.*?\]\s+(\d+\.\d+) %",
                                      r" 11 Choose ESP type, current: Nuclear \+ Electronic[\n\r]?",
                                      r".*\(y/n\)[\n\r]?", pexpect.TIMEOUT], 
                                      timeout=5)
                if index == 0:  # Progress line matched
                    progress_bar.update(1)  # Update the tqdm bar by 1 tick
                
                elif index in [1, 2]:  # (process finished)
                    break
                
                elif index == 3:  # Timeout (no output for 5 seconds)
                    break

        except pexpect.ExceptionPexpect as e:
            raise e

        # Close the child process
        child.close()

    return outputFile
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
