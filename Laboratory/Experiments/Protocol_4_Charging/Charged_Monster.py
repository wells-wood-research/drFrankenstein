
## BASIC IMPORTS ##
from os import path as p
import os
import glob
from shutil import copy
from subprocess import call, DEVNULL
from pdbUtils import pdbUtils
import pandas as pd
import pexpect
from tqdm import tqdm
from copy import deepcopy
## drFrankenstein LIBRARIES ##
from OperatingTools import drOrca
from Experiments.Protocol_1_Capping.Capping_Assistant import find_bonded_atoms

## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

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
def run_charge_fitting(config: dict,
                        conformerListTxt: FilePath,
                          chargeConstraintsTxt: FilePath,
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

        ####### RUN CALCULATIONS #######
        child.expect(r'.*11.*\n?\r?')
        child.sendline("2")

        # Monitor the output for progress updates

        # Initialize tqdm progress bar
        tqdmBarOptions = {
        "desc": f"Performing Charge Fitting",
        "ascii": "-ϟ",  
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
                    continue

        except pexpect.ExceptionPexpect as e:
            raise e

        # Close the child process
        child.close()

    return outputFile
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
